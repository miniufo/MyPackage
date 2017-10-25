/**
 * @(#)GridDataFetcher.java	1.0 2013.02.17
 *
 * Copyright 2013 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import static miniufo.basic.InterpolationModel.linearInterpolation;
import static miniufo.basic.InterpolationModel.bilinearInterpolation;


/**
 * fetching data at any point from gridded data source
 * using a bi-linear or bi-cubic interpolation fashion
 *
 * @version 1.0, 2013.02.17
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GridDataFetcher{
	//
	private int tlev=-1;	// started   time-level for current buffer (start from 1)
	private int zlev=-1;	// started height-level for current buffer (start from 1)
	
	private float undef=Float.NaN;
	
	private float[] xdef =null;
	private float[] ydef =null;
	private float[] zdef =null;
	private float[] dxdef=null;
	private float[] dydef=null;
	private float[] dzdef=null;
	private float[] dtdef=null;
	private MDate[] tdef =null;
	
	private Region2D region=null;
	
	private DataDescriptor dd=null;
	
	private DataRead dr=null;
	
	
	/**
	 * constructor
	 */
	public GridDataFetcher(DataDescriptor dd){
		this.dd=dd;
		
		xdef=dd.getXDef().getSamples();	dxdef=dd.getDXDef();
		ydef=dd.getYDef().getSamples();	dydef=dd.getDYDef();
		zdef=dd.getZDef().getSamples();	dzdef=dd.getDZDef();
		tdef=dd.getTDef().getSamples();	dtdef=dd.getDTDef();
		
		region=dd.toRegion2D();
		
		dr=DataIOFactory.getDataRead(dd);
	}
	
	
	/*** getor and setor ***/
	public int getBufferTLevel(){ return tlev;}
	
	public int getBufferZLevle(){ return zlev;}
	
	
	/**
	 * prepare the horizontal slice buffer before fetching it
	 * 
	 * @param	vname	name of Variable
	 * @param	tstep	t level (start from 1)
	 * @param	zstep	z level (start from 1)
	 * 
	 * @return	xybuf	2-D buffer data
	 */
	public Variable prepareXYBuffer(String vname,int tstep,int zstep){
		return prepareXYBuffer(vname,tstep,zstep,0);
	}
	
	/**
	 * prepare the horizontal slice buffer before fetching it
	 * 
	 * @param	vname	name of Variable
	 * @param	tstep	t level (start from 1)
	 * @param	zstep	z level (start from 1)
	 * @param	expand	times of horizontal expanding to reduce undefined data
	 * 
	 * @return	xybuf	2-D buffer data
	 */
	public Variable prepareXYBuffer(String vname,int tstep,int zstep,int expand){
		if(tstep<1||tstep>dd.getTCount())
		throw new IllegalArgumentException("t-step ("+tstep+") should be in [1 "+dd.getTCount()+"]");
		
		if(zstep<1||zstep>dd.getZCount())
		throw new IllegalArgumentException("z-step ("+zstep+") should be in [1 "+dd.getZCount()+"]");
		
		Variable xybuf=new Variable(vname,true,new Range(1,1,dd.getYCount(),dd.getXCount()));
		
		Range r=xybuf.getRange();
		
		r.setTRange(tstep);
		r.setZRange(zstep);
		
		dr.readData(xybuf);
		
		tlev=tstep;
		zlev=zstep;
		
		expandXY(xybuf,expand);
		
		return xybuf;
	}
	
	
	/**
	 * prepare the XYZ buffer before fetching it
	 * 
	 * @param	vname	name of Variable
	 * @param	tstep	t level (start from 1)
	 * 
	 * @return	xyzbuf	3-D buffer data
	 */
	public Variable prepareXYZBuffer(String vname,int tstep){
		return prepareXYZBuffer(vname,tstep,1,dd.getZCount(),0);
	}
	
	/**
	 * prepare the XYZ buffer before fetching it
	 * 
	 * @param	vname	name of Variable
	 * @param	tstep	t level (start from 1)
	 * @param	zstr	z start level (start from 1)
	 * @param	zlen	length of z levels
	 * @param	expand	times of horizontal expanding to reduce undefined data
	 * 
	 * @return	xyzbuf	3-D buffer data
	 */
	public Variable prepareXYZBuffer(String vname,int tstep,int zstr,int zlen,int expand){
		if(tstep<1||tstep>dd.getTCount())
		throw new IllegalArgumentException("t-step should be in [1 "+dd.getTCount()+"]");
		if(zstr<1||zstr>dd.getZCount())
		throw new IllegalArgumentException("zstr should be in [1 "+dd.getZCount()+"]");
		if(zlen<1||zlen>dd.getZCount())
		throw new IllegalArgumentException("zlen should be in [1 "+dd.getZCount()+"]");
		
		Variable xyzbuf=new Variable(vname,true,new Range(1,zlen,dd.getYCount(),dd.getXCount()));
		
		xyzbuf.getRange().setTRange(tstep);
		
		int[] zrange=xyzbuf.getRange().getZRange();
		zrange[0]=zstr;
		zrange[1]=zstr+zlen-1;
		zrange[2]=zlen;
		
		dr.readData(xyzbuf);
		
		tlev=tstep;
		
		expandXY(xyzbuf,expand);
		
		return xyzbuf;
	}
	
	
	/**
	 * prepare the XY-T buffer before fetching it
	 * 
	 * @param	vname	name of Variable
	 * @param	zstep	level (start from 1)
	 * @param	expand	times of horizontal expanding to reduce undefined data
	 * 
	 * @return	xytbuf	3-D buffer data
	 */
	public Variable prepareXYTBuffer(String vname,int zstep){
		return prepareXYTBuffer(vname,zstep,1,dd.getTCount(),0);
	}
	
	/**
	 * prepare the XY-T buffer before fetching it
	 * 
	 * @param	vname	name of Variable
	 * @param	zstep	level (start from 1)
	 * @param	tstr	t start level (start from 1)
	 * @param	tlen	length of t levels
	 * @param	expand	times of horizontal expanding to reduce undefined data
	 * 
	 * @return	xytbuf	3-D buffer data
	 */
	public Variable prepareXYTBuffer(String vname,int zstep,int tstr,int tlen,int expand){
		if(zstep<1||zstep>dd.getZCount())
		throw new IllegalArgumentException("z-step should be in [1 "+dd.getZCount()+"]");
		if(tstr<1||tstr>dd.getTCount())
		throw new IllegalArgumentException("tstr should be in [1 "+dd.getTCount()+"]");
		if(tlen<1||tlen>dd.getTCount())
		throw new IllegalArgumentException("tlen should be in [1 "+dd.getTCount()+"]");
		
		Variable xytbuf=new Variable(vname,false,new Range(tlen,1,dd.getYCount(),dd.getXCount()));
		
		xytbuf.getRange().setZRange(zstep);
		
		int[] trange=xytbuf.getRange().getTRange();
		trange[0]=tstr;
		trange[1]=tstr+tlen-1;
		trange[2]=tlen;
		
		dr.readData(xytbuf);
		
		zlev=zstep;
		tlev=tstr;
		
		expandXYT(xytbuf,expand);
		
		return xytbuf;
	}
	
	
	/**
	 * fetch XY-slice data from prepared buffer
	 * 
	 * @param	xpos	x-position (unit of degree or m) of the point to be fetch
	 * @param	ypos	y-position (unit of degree or m) of the point to be fetch
	 * @param	xybuf	XY-slice buffer obtained from calling prepareXYSliceBuffer
	 */
	public float fetchXYBuffer(float xpos,float ypos,Variable xybuf){
		undef=xybuf.getUndef();
		
		float[][] xyb=xybuf.getData()[0][0];
		
		if(region.inRange(xpos,ypos)){
			int xtag=dd.getXLENum(xpos);
			int ytag=dd.getYLENum(ypos);
			int xend=xtag==dd.getXCount()-1?xtag:xtag+1;
			int yend=ytag==dd.getYCount()-1?ytag:ytag+1;
			
			float lt=xyb[yend][xtag];	if(lt==undef) return undef;
			float lb=xyb[ytag][xtag];	if(lb==undef) return undef;
			float rb=xyb[ytag][xend];	if(rb==undef) return undef;
			float rt=xyb[yend][xend];	if(rt==undef) return undef;
			
			float dx1=(xpos-xdef[xtag])/dxdef[xend-1];
			float dy1=(ypos-ydef[ytag])/dydef[yend-1];
			
			return bilinearInterpolation(lb,rb,lt,rt,dx1,dy1,undef);
			
		}else return undef;
	}
	
	public float fetchXYBufferPeriodicX(float xpos,float ypos,Variable xybuf){
		undef=xybuf.getUndef();
		
		float[][] xyb=xybuf.getData()[0][0];
		
		if(region.inYRange(ypos)){
			int xtag=dd.getXLENumPeriodicX(xpos);
			int ytag=dd.getYLENum(ypos);
			int xend=xtag==dd.getXCount()-1?0:xtag+1;
			int yend=ytag==dd.getYCount()-1?ytag:ytag+1;
			
			float lt=xyb[yend][xtag];	if(lt==undef) return undef;
			float lb=xyb[ytag][xtag];	if(lb==undef) return undef;
			float rb=xyb[ytag][xend];	if(rb==undef) return undef;
			float rt=xyb[yend][xend];	if(rt==undef) return undef;
			
			float dx1=(xpos-xdef[xtag])/dxdef[xend-1<0?0:xend-1];
			float dy1=(ypos-ydef[ytag])/dydef[yend-1];
			
			return bilinearInterpolation(lb,rb,lt,rt,dx1,dy1,undef);
			
		}else return undef;
	}
	
	
	/**
	 * fetch XYZ-slice data from prepared buffer
	 * 
	 * @param	xpos	x-position (unit of degree or m) of the point to be fetch
	 * @param	ypos	y-position (unit of degree or m) of the point to be fetch
	 * @param	zpos	z-position (unit of Pa or m) of the point to be fetch
	 * @param	xyzbuf	XYZ-slice buffer obtained from calling prepareXYZBuffer
	 */
	public float fetchXYZBuffer(float xpos,float ypos,float zpos,Variable xyzbuf){
		undef=xyzbuf.getUndef();
		
		float[][][] xyzb=xyzbuf.getData()[0];
		
		if(region.inRange(xpos,ypos)){
			int xtag=dd.getXLENum(xpos);
			int ytag=dd.getYLENum(ypos);
			int ztag=dd.getZLENum(zpos);	//int ztag=dd.getZLENum(lev+(lev==1000?0:0.1f));
			if(ztag==-1) throw new IllegalArgumentException(
				"zpos "+zpos+" outside ["+dd.getZDef().getMin()+","+dd.getZDef().getMax()+"]"
			);
			
			int xend=xtag==dd.getXCount()-1?xtag:xtag+1;
			int yend=ytag==dd.getYCount()-1?ytag:ytag+1;
			int zend=ztag==dd.getZCount()-1?ztag:ztag+1;
			
			float ltl=xyzb[ztag][yend][xtag];	if(ltl==undef) return undef;
			float lbl=xyzb[ztag][ytag][xtag];	if(lbl==undef) return undef;
			float rbl=xyzb[ztag][ytag][xend];	if(rbl==undef) return undef;
			float rtl=xyzb[ztag][yend][xend];	if(rtl==undef) return undef;
			
			float ltu=xyzb[zend][yend][xtag];	if(ltu==undef) return undef;
			float lbu=xyzb[zend][ytag][xtag];	if(lbu==undef) return undef;
			float rbu=xyzb[zend][ytag][xend];	if(rbu==undef) return undef;
			float rtu=xyzb[zend][yend][xend];	if(rtu==undef) return undef;
			
			float dx1=(xpos-xdef[xtag])/dxdef[xend-1];
			float dy1=(ypos-ydef[ytag])/dydef[yend-1];
			float dz1=(zpos-zdef[ztag])/dzdef[zend-1];
			
			float lower=bilinearInterpolation(lbl,rbl,ltl,rtl,dx1,dy1,undef);
			float upper=bilinearInterpolation(lbu,rbu,ltu,rtu,dx1,dy1,undef);
			
			return linearInterpolation(lower,upper,dz1,undef);
			
		}else return undef;
	}
	
	public float fetchXYZBufferPeriodicX(float xpos,float ypos,float zpos,Variable xyzbuf){
		undef=xyzbuf.getUndef();
		
		float[][][] xyzb=xyzbuf.getData()[0];
		
		if(region.inYRange(ypos)){
			int xtag=dd.getXLENumPeriodicX(xpos);
			int ytag=dd.getYLENum(ypos);
			int ztag=dd.getZLENum(zpos);
			if(ztag==-1) throw new IllegalArgumentException(
				"zpos "+zpos+" outside ["+dd.getZDef().getMin()+","+dd.getZDef().getMax()+"]"
			);
			
			int xend=xtag==dd.getXCount()-1?0:xtag+1;
			int yend=ytag==dd.getYCount()-1?ytag:ytag+1;
			int zend=ztag==dd.getZCount()-1?ztag:ztag+1;
			
			float ltl=xyzb[ztag][yend][xtag];	if(ltl==undef) return undef;
			float lbl=xyzb[ztag][ytag][xtag];	if(lbl==undef) return undef;
			float rbl=xyzb[ztag][ytag][xend];	if(rbl==undef) return undef;
			float rtl=xyzb[ztag][yend][xend];	if(rtl==undef) return undef;
			
			float ltu=xyzb[zend][yend][xtag];	if(ltu==undef) return undef;
			float lbu=xyzb[zend][ytag][xtag];	if(lbu==undef) return undef;
			float rbu=xyzb[zend][ytag][xend];	if(rbu==undef) return undef;
			float rtu=xyzb[zend][yend][xend];	if(rtu==undef) return undef;
			
			float dx1=(xpos-xdef[xtag])/dxdef[xend-1<0?0:xend-1];
			float dy1=(ypos-ydef[ytag])/dydef[yend-1];
			float dz1=(zpos-zdef[ztag])/dzdef[zend-1];
			
			float lower=bilinearInterpolation(lbl,rbl,ltl,rtl,dx1,dy1,undef);
			float upper=bilinearInterpolation(lbu,rbu,ltu,rtu,dx1,dy1,undef);
			
			return linearInterpolation(lower,upper,dz1);
			
		}else return undef;
	}
	
	
	/**
	 * fetch XYT-slice data from prepared buffer
	 * 
	 * @param	xpos	x-position (unit of degree or m) of the point to be fetch
	 * @param	ypos	y-position (unit of degree or m) of the point to be fetch
	 * @param	tim		time in long format
	 * @param	xytbuf	XYT-slice buffer obtained from calling prepareXYSliceBuffer
	 */
	public float fetchXYTBuffer(float xpos,float ypos,long tim,Variable xytbuf){
		undef=xytbuf.getUndef();
		
		float[][][] xytb=xytbuf.getData()[0];
		
		if(region.inRange(xpos,ypos)){
			int xtag=dd.getXLENum(xpos);
			int ytag=dd.getYLENum(ypos);
			int ttag=dd.getTLENum(tim)-tlev+1;
			int xend=xtag==dd.getXCount()-1?xtag:xtag+1;
			int yend=ytag==dd.getYCount()-1?ytag:ytag+1;
			int tend=ttag==xytb[0][0].length-1?ttag:ttag+1;
			
			if(ttag==-1) throw new IllegalArgumentException(
				"time "+tim+" outside ["+dd.getTDef().getFirst().getLongTime()+","+dd.getTDef().getLast().getLongTime()+"]"
			);
			
			float ltp=xytb[yend][xtag][ttag];	if(ltp==undef) return undef;
			float lbp=xytb[ytag][xtag][ttag];	if(lbp==undef) return undef;
			float rbp=xytb[ytag][xend][ttag];	if(rbp==undef) return undef;
			float rtp=xytb[yend][xend][ttag];	if(rtp==undef) return undef;
			
			float ltn=xytb[yend][xtag][tend];	if(ltn==undef) return undef;
			float lbn=xytb[ytag][xtag][tend];	if(lbn==undef) return undef;
			float rbn=xytb[ytag][xend][tend];	if(rbn==undef) return undef;
			float rtn=xytb[yend][xend][tend];	if(rtn==undef) return undef;
			
			float dx1=(xpos-xdef[xtag])/dxdef[xend-1];
			float dy1=(ypos-ydef[ytag])/dydef[yend-1];
			float dt1=new MDate(tim).getDT(tdef[ttag+tlev-1])/dtdef[0];
			
			float prev=bilinearInterpolation(lbp,rbp,ltp,rtp,dx1,dy1,undef);
			float next=bilinearInterpolation(lbn,rbn,ltn,rtn,dx1,dy1,undef);
			
			return linearInterpolation(prev,next,dt1);
			
		}else return undef;
	}
	
	public float fetchXYTBufferPeriodicX(float xpos,float ypos,long tim,Variable xytbuf){
		undef=xytbuf.getUndef();
		
		float[][][] xytb=xytbuf.getData()[0];
		
		if(region.inYRange(ypos)){
			int xtag=dd.getXLENumPeriodicX(xpos);
			int ytag=dd.getYLENum(ypos);
			int ttag=dd.getTLENum(tim)-tlev+1;
			int xend=xtag==dd.getXCount()-1?0:xtag+1;
			int yend=ytag==dd.getYCount()-1?ytag:ytag+1;
			int tend=ttag==xytb[0][0].length-1?ttag:ttag+1;
			
			if(ttag==-1) throw new IllegalArgumentException(
				"time "+tim+" outside ["+dd.getTDef().getFirst().getLongTime()+","+dd.getTDef().getLast().getLongTime()+"]"
			);
			
			float ltp=xytb[yend][xtag][ttag];	if(ltp==undef) return undef;
			float lbp=xytb[ytag][xtag][ttag];	if(lbp==undef) return undef;
			float rbp=xytb[ytag][xend][ttag];	if(rbp==undef) return undef;
			float rtp=xytb[yend][xend][ttag];	if(rtp==undef) return undef;
			
			float ltn=xytb[yend][xtag][tend];	if(ltn==undef) return undef;
			float lbn=xytb[ytag][xtag][tend];	if(lbn==undef) return undef;
			float rbn=xytb[ytag][xend][tend];	if(rbn==undef) return undef;
			float rtn=xytb[yend][xend][tend];	if(rtn==undef) return undef;
			
			float dx1=(xpos-xdef[xtag])/dxdef[xend-1<0?0:xend-1];
			float dy1=(ypos-ydef[ytag])/dydef[yend-1];
			float dt1=new MDate(tim).getDT(tdef[ttag+tlev-1])/dtdef[0];
			
			float prev=bilinearInterpolation(lbp,rbp,ltp,rtp,dx1,dy1,undef);
			float next=bilinearInterpolation(lbn,rbn,ltn,rtn,dx1,dy1,undef);
			
			return linearInterpolation(prev,next,dt1,undef);
			
		}else return undef;
	}
	
	
	/**
	 * close the DataRead
	 */
	public void closeFile(){ dr.closeFile();}
	
	
	/*** helper methods ***/
	private static void expandXY(Variable v,int iter){
		float undef=v.getUndef();
		
		for(int m=0;m<iter;m++)
		for(int l=0,L=v.getTCount();l<L;l++)
		for(int k=0,K=v.getZCount();k<K;k++){
			float[][] vdata=v.getData()[l][k];
			float[][] bdata=new float[v.getYCount()][v.getXCount()];
			
			for(int j=0,J=v.getYCount();j<J;j++) System.arraycopy(vdata[j],0,bdata[j],0,v.getXCount());
			
			int jtag=0,itag=0;
			
			for(int j=0,J=v.getYCount();j<J;j++)
			for(int i=0,I=v.getXCount();i<I;i++) if(bdata[j][i]==undef){
				int count=0; double sum=0;
				
				jtag=j+1; itag=i+1; if(i==I-1) itag=0;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				jtag=j  ; itag=i+1; if(i==I-1) itag=0;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				jtag=j-1; itag=i+1; if(i==I-1) itag=0;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				jtag=j+1; itag=i  ;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				jtag=j-1; itag=i  ;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				jtag=j+1; itag=i-1; if(i==0) itag=I-1;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				jtag=j  ; itag=i-1; if(i==0) itag=I-1;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				jtag=j-1; itag=i-1; if(i==0) itag=I-1;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag]!=undef){ sum+=bdata[jtag][itag]; count++;}
				
				if(count!=0) vdata[j][i]=(float)(sum/count);
			}
		}
	}
	
	private static void expandXYT(Variable v,int iter){
		float undef=v.getUndef();
		
		for(int m=0;m<iter;m++)
		for(int k=0,K=v.getZCount();k<K;k++){
			float[][][] vdata=v.getData()[k];
			float[][][] bdata=new float[v.getYCount()][v.getXCount()][v.getTCount()];
			
			for(int j=0,J=v.getYCount();j<J;j++)
			for(int i=0,I=v.getXCount();i<I;i++)
			System.arraycopy(vdata[j][i],0,bdata[j][i],0,v.getTCount());
			
			int jtag=0,itag=0;
			
			for(int j=0,J=v.getYCount();j<J;j++)
			for(int i=0,I=v.getXCount();i<I;i++)
			for(int l=0,L=v.getTCount();l<L;l++) if(bdata[j][i][l]==undef){
				int count=0; double sum=0;
				
				jtag=j+1; itag=i+1; if(i==I-1) itag=0;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				jtag=j  ; itag=i+1; if(i==I-1) itag=0;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				jtag=j-1; itag=i+1; if(i==I-1) itag=0;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				jtag=j+1; itag=i  ;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				jtag=j-1; itag=i  ;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				jtag=j+1; itag=i-1; if(i==0) itag=I-1;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				jtag=j  ; itag=i-1; if(i==0) itag=I-1;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				jtag=j-1; itag=i-1; if(i==0) itag=I-1;
				if(jtag<J&&jtag>=0&&bdata[jtag][itag][l]!=undef){ sum+=bdata[jtag][itag][l]; count++;}
				
				if(count!=0) vdata[j][i][l]=(float)(sum/count);
			}
		}
	}
	
	
	/** test
	public static void main(String argv[]){
		miniufo.diagnosis.DiagnosisFactory df=
		miniufo.diagnosis.DiagnosisFactory.parseFile("d:/Data/HYCOM/HYCOM_SCS_Clim_1993_2003.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		GridDataFetcher gdf=new GridDataFetcher(dd);
		
		Variable t1=df.getVariables(new Range("t(1,1)",dd),"t")[0]; t1.replaceUndefData(0);
		Variable t2=new Variable("t" ,true,new Range(1,dd.getZCount(),dd.getYCount(),dd.getXCount()));
		
		float[][][] data=gdf.prepareXYZBuffer(t2.getName(),1,1,dd.getZCount(),105);
		
		Variable t3=new Variable("t3",true,new Range(1,dd.getZCount(),dd.getYCount(),dd.getXCount()));
		
		t3.getData()[0]=data; t3.replaceUndefData(0);
		
		miniufo.io.DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/expandXY.dat");
		dw.writeData(dd,t1,t3); dw.closeFile();
	}*/
}
