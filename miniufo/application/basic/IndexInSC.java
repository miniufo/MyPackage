/**
 * @(#)IndexInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import miniufo.lagrangian.Typhoon;
import miniufo.util.TicToc;

import java.util.concurrent.TimeUnit;

import miniufo.application.EquationInSphericalCoordinate;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.basic.ArrayUtil;


/**
 * All kinds of climate index
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class IndexInSC extends EquationInSphericalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public IndexInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/**
     * calculate blocking index according to Pelly and Hoskins (2003)
     *
     * @param	hgt		geopotential height field
     * @param	xspan	grids span in x-direction, even number only
     * @param	yspan	grids span in y-direction, even number only
     *
     * @return	blocking index field
     */
	public Variable cBlockingIndexByPelly(Variable hgt,int xspan,int yspan){
		if(xspan<2||xspan%2!=0) throw new IllegalArgumentException("xspan should be even and larger than 1");
		if(yspan<2||yspan%2!=0) throw new IllegalArgumentException("yspan should be even and larger than 1");
		
		xspan/=2;	yspan/=2;	ystart=hgt.getRange().getYRange()[0];
		
		t=hgt.getTCount();	z=hgt.getZCount();	y=hgt.getYCount();	x=hgt.getXCount();
		
		float undef=hgt.getUndef();
		
		Variable bi=new Variable("bi",hgt);	bi.setValue(undef);
		bi.setCommentAndUnit("blocking index by Pelly and Hoskins (2003)");
		
		float[][][][] hgtdata=hgt.getData();	float[][][][] bidata=bi.getData();
		
		if(hgt.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=yspan;j<y-yspan;j++)
			for(int i=xspan;i<x-xspan;i++){
				bidata[l][k][j][i]=0;
				
				for(int jj=j-yspan;jj<=j;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[l][k][j][i]-=hgtdata[l][k][jj][ii];
				
				for(int jj=j;jj<=j+yspan;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[l][k][j][i]+=hgtdata[l][k][jj][ii];
				
				bidata[l][k][j][i]/=(yspan+1)*(xspan*2+1);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=yspan;j<y-yspan;j++)
			for(int i=xspan;i<x-xspan;i++){
				bidata[k][j][i][l]=0;
				
				for(int jj=j-yspan;jj<=j;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[k][j][i][l]-=hgtdata[k][jj][ii][l];
				
				for(int jj=j;jj<=j+yspan;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[k][j][i][l]+=hgtdata[k][jj][ii][l];
				
				bidata[k][j][i][l]/=(yspan+1)*(xspan*2+1);
			}
		}
		
		return bi;
	}
	
	
	/**
     * Calculate 2-dimensional (horizontal) index by assigning a properly averaged
     * diagnostics in cylindrical coordinates to each lat/lon grid.
     * 
     * Reference: (Qian et al. 2016, WAF)
     *
     * @param	dd			grid data descriptor
     * @param	rng			range of calculation, should be smaller than that of descriptor
     * @param	tr			a given Typhoon
     * @param	yinc		increment in radial direction (degree)
     * @param	ycount		ring count in cylindrical coordinates
     * @param	xcount		grid count in each ring
     * @param	idxNames	index names: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	2-dimensional horizontal indices
     */
	public static Variable[] c2DHorizontalIndex
	(DataDescriptor dd,String rng,Typhoon tr,float yinc,int ycount,int xcount,int str,int end,String... idxNames){
		Range r1=new Range(rng,dd);
		Range r2=new Range(tr.getTRange(),dd);
		
		MDate[] tdef=dd.getTDef().getSamples();
		
		if(r1.getTRange()[0]!=r2.getTRange()[0]) throw new IllegalArgumentException(
			"typhoon start time ("+tdef[r2.getTRange()[0]-1]+") is not the same as range-specified time ("+tdef[r1.getTRange()[0]-1]+")"
		);
		
		if(r1.getTRange()[1]!=r2.getTRange()[1]) throw new IllegalArgumentException(
			"typhoon end time ("+tdef[r2.getTRange()[1]-1]+") is not the same as range-specified time ("+tdef[r1.getTRange()[1]-1]+")"
		);
		
		return c2DHorizontalIndex(dd,rng,tr.getZonalVelocity(),tr.getMeridionalVelocity(),yinc,ycount,xcount,str,end,idxNames);
	}
	
	/**
     * Calculate 2-dimensional (horizontal) index by assigning a properly averaged
     * diagnostics in cylindrical coordinates to each lat/lon grid.
     * 
     * Reference: (Qian et al. 2016, WAF)
     *
     * @param	dd			grid data descriptor
     * @param	rng			range of calculation, should be smaller than that of descriptor
     * @param	yinc		increment in radial direction (degree)
     * @param	ycount		ring count in cylindrical coordinates
     * @param	xcount		grid count in each ring
     * @param	idxNames	index names: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	2-dimensional horizontal indices
     */
	public static Variable[] c2DHorizontalIndex
	(DataDescriptor dd,String rng,float yinc,int ycount,int xcount,int str,int end,String... idxNames){
		return c2DHorizontalIndex(dd,rng,0,0,yinc,ycount,xcount,str,end,idxNames);
	}
	
	/**
     * Calculate 2-dimensional (horizontal) index by assigning a properly averaged
     * diagnostics in cylindrical coordinates to each lat/lon grid.
     * 
     * Reference: (Qian et al. 2016, WAF)
     *
     * @param	dd			grid data descriptor
     * @param	rng			range of calculation, should be smaller than that of descriptor
     * @param	cu			zonal velocity of cylindrical coordinates
     * @param	cv			meridional velocity of cylindrical coordinates
     * @param	yinc		increment in radial direction (degree)
     * @param	ycount		ring count in cylindrical coordinates
     * @param	xcount		grid count in each ring
     * @param	idxNames	index names: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	2-dimensional horizontal indices
     */
	public static Variable[] c2DHorizontalIndex
	(DataDescriptor dd,String rng,float cu,float cv,float yinc,int ycount,int xcount,int str,int end,String... idxNames){
		Range r=new Range(rng,dd);
		
		int tlen=r.getTRange()[2];
		
		float[] cus=ArrayUtil.newArray(tlen,cu);
		float[] cvs=ArrayUtil.newArray(tlen,cv);
		
		return c2DHorizontalIndex(dd,rng,cus,cvs,yinc,ycount,xcount,str,end,idxNames);
	}
	
	
	/*** helper methods ***/
	
	/**
     * get the range without the lon and lat
     *
     * @return	the range String
     */
	private static String getTZBufferRange(String rng){
		String[] ss=rng.split(";");
		StringBuilder sb=new StringBuilder();
		
		for(String s:ss)
		if(s.startsWith("t")||s.startsWith("z")||s.startsWith("lev")){
			sb.append(s);	sb.append(";");
		}
		
		if(sb.length()!=0) sb.delete(sb.length()-1,sb.length());
		
		return sb.toString();
	}
	
	/**
     * Calculate 2-dimensional (horizontal) index by assigning a properly averaged
     * diagnostics in cylindrical coordinates to each lat/lon grid.
     * 
     * Reference: (Qian et al. 2016, WAF)
     *
     * @param	dd			grid data descriptor
     * @param	rng			range of calculation, should be smaller than that of descriptor
     * @param	cus			zonal velocities of cylindrical coordinates
     * @param	cvs			meridional velocities of cylindrical coordinates
     * @param	yinc		increment in radial direction (degree)
     * @param	ycount		ring count in cylindrical coordinates
     * @param	xcount		grid count in each ring
     * @param	str			start radial index (from 0)
     * @param	end			end   radial index (from 0)
     * @param	idxNames	index names: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	2-dimensional horizontal indices
     */
	private static Variable[] c2DHorizontalIndex
	(DataDescriptor dd,String rng,float[] cus,float[] cvs,float yinc,int ycount,int xcount,int str,int end,String... idxNames){
		if(idxNames.length==0) throw new IllegalArgumentException("index name required");
		
		boolean hasT=false;
		
		for(String s:idxNames) if(s.equals("htHFC")){ hasT=true; break;}
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		
		Range r   =new Range(rng,dd);
		Range rtmp=new Range(getTZBufferRange(rng),dd);
		
		int tlen=r.getTRange()[2];
		int outputInterval=Math.round(r.getYRange()[2]*r.getXRange()[2]/20f);
		
		if(cus.length!=tlen) throw new IllegalArgumentException("length of cus do not equal time count");
		if(cvs.length!=tlen) throw new IllegalArgumentException("length of cvs do not equal time count");
		
		Variable u=new Variable("u",rtmp);
		Variable v=new Variable("v",rtmp);
		Variable T=hasT?new Variable("T",rtmp):null;
		
		DataRead dr=DataIOFactory.getDataRead(dd);
		dr.setPrinting(false);
		if(hasT) dr.readData(u,v,T);
		else     dr.readData(u,v);
		dr.closeFile();
		
		int tc=idxNames.length;
		Variable[] idx=new Variable[tc];
		for(int m=0;m<tc;m++){
			idx[m]=new Variable(idxNames[m],r);
			idx[m].setUndef(dd.getUndef("u"));
			idx[m].setCommentAndUnit("2D horizontal index ("+idxNames[m]+")");
		}
		
		TicToc.tic("start computing index");
		
		for(int ystart=r.getYRange()[0]-1,jj=ystart,J=r.getYRange()[1],ptr=0;jj<J;jj++)
		for(int xstart=r.getXRange()[0]-1,ii=xstart,I=r.getXRange()[1];ii<I;ii++){ ptr++;
			float[] olons=ArrayUtil.newArray(tlen,dd.getXDef().getSamples()[ii]);
			float[] olats=ArrayUtil.newArray(tlen,dd.getYDef().getSamples()[jj]);
			
			CylindricalSpatialModel csm=new CylindricalSpatialModel(olons,olats,dd.getZDef().getSamples(),yinc,ycount,xcount);
			csm.setUWhole(cus); csm.setVWhole(cvs);
			
			Variable[] CIDX=cIndexAtPoint(csm,ssm,u,v,T,str,end,idxNames);
			
			for(int m=0;m<tc;m++){ // copy cylindrical result to the lat/lon grids
				float[][][][] edata=CIDX[m].getData();
				float[][][][] idata= idx[m].getData();
				
				if(idx[0].isTFirst()){
					for(int l=0,L=idx[0].getTCount();l<L;l++)
					for(int k=0,K=idx[0].getZCount();k<K;k++)
					idata[l][k][jj-ystart][ii-xstart]=edata[l][k][0][0];
					
				}else{
					for(int l=0,L=u.getTCount();l<L;l++)
					for(int k=0,K=u.getZCount();k<K;k++)
					idata[k][jj-ystart][ii-xstart][l]=edata[k][0][0][l];
				}
			}
			
			if(ptr%(outputInterval)==0) System.out.print(".");
		}
		
		TicToc.toc(TimeUnit.MINUTES);
		
		return idx;
	}
	
	/**
     * Calculate cylindrical-coordinate diagnostics for a single lat/lon grid.
     *
     * @param	csm			cylindrical spatial model
     * @param	ssm			spherical   spatial model
     * @param	u			gridded zonal wind
     * @param	v			gridded meridional wind
     * @param	T			gridded temperature
     * @param	idxNames	index names: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	radial-band averaged diagnostics
     */
	private static Variable[] cIndexAtPoint
	(CylindricalSpatialModel csm,SphericalSpatialModel ssm,Variable u,Variable v,Variable T,int str,int end,String... idxNames){
		int ic=idxNames.length;
		
		CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
		DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
		ThermoDynamicMethodsInCC tdm=new ThermoDynamicMethodsInCC(csm);
		
		Variable[] re=ct.reprojectToCylindrical(ct.transToCylindrical(u),ct.transToCylindrical(v));
		Variable Vr=re[1].copy();	// full radial velocity, not storm-relative
		Variable Tc=T==null?null:ct.transToCylindrical(T);
		
		dm.cStormRelativeAziRadVelocity(csm.getUWhole(),csm.getVWhole(),re[0],re[1]);
		Variable utm=re[0].anomalizeX();
		Variable uta=re[0]; re[1].anomalizeX();
		Variable vra=re[1];
		
		Variable[] CIDX=new Variable[ic];
		
		// for 0.3-deg interval
		//int str= 9,end=18;	// 300-600 km average
		//int str=12,end=21;	// 400-700 km average
		//int str=15,end=24;	// 500-800 km average
		
		for(int m=0;m<ic;m++){
			if(idxNames[m].equalsIgnoreCase("REFC")){
				CIDX[m]=dm.cREFC(uta,vra).averageAlong(Dimension.Y,str,end);
				
			}else if(idxNames[m].equalsIgnoreCase("PEFC")){
				Vr.anomalizeX(); // not storm-relative radial velocity
				CIDX[m]=dm.cPEFC(Vr).averageAlong(Dimension.Y,str,end);
				
			}else if(idxNames[m].equalsIgnoreCase("AEFC")){
				Variable gaz=dm.cAbsoluteAngularMomentumByJohnson(re[0]);
				gaz.anomalizeX();
				CIDX[m]=dm.cAEFC(gaz,vra).averageAlong(Dimension.Y,str,end);
				
			}else if(idxNames[m].equalsIgnoreCase("ISB")){
				CIDX[m]=dm.cMeanInertialStabilityByUT(utm).averageAlong(Dimension.Y,str,end);
				
			}else if(idxNames[m].equalsIgnoreCase("ETA")){
				CIDX[m]=dm.cMeanAbsoluteVorticity(utm).averageAlong(Dimension.Y,str,end);
				
			}else if(idxNames[m].equalsIgnoreCase("ULFI")){
				Vr.anomalizeX(); // not storm-relative radial velocity
				CIDX[m]=dm.cREFC(uta,vra).plusEq(dm.cPEFC(Vr)).divideEq(dm.cMeanAbsoluteVorticity(utm)).divideEq(86400f).averageAlong(Dimension.Y,str,end);
				
			}else if(idxNames[m].equalsIgnoreCase("htHFC")){
				Variable th=tdm.cPotentialTemperature(Tc);
				th.anomalizeX();
				CIDX[m]=dm.cEddyHeatHFC(th,vra).averageAlong(Dimension.Y,str,end);
				
			}else throw new IllegalArgumentException("invalid type for horizontal index: "+idxNames[m]);
		}
		
		return CIDX;
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			miniufo.diagnosis.DiagnosisFactory df=new miniufo.diagnosis.DiagnosisFactory(
				"E:/ERAInterim/200/2004.200.uv.nc"
			);
			DataDescriptor dd=df.getDataDescriptor();
			
			miniufo.io.DataWrite dw=miniufo.io.DataIOFactory.getDataWrite(
				dd,"E:/ERAInterim/200/PRI2004.dat"
			);
			
			int tt=1;
			while(tt<=1464){
				String r="t("+tt+","+(tt+23)+");lon(80,180);lat(6,60)";
				System.out.println(r);
				
				Variable[] vel=df.getVariables(new Range(r,dd),"u","v");
				
				Variable pri=IndexInSC.cPotentialReintensifyIndex(dd,r,0.2f,30,36);
				
				vel[0].setUndef(vel[1].getUndef());
				pri.setUndef(vel[1].getUndef());
				
				dw.writeData(pri);
				
				tt+=24;
			}
			
			dw.writeCtl(dd);	dw.closeFile();
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
