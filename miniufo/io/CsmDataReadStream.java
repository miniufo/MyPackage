/**
 * @(#)CsmDataReadStream.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import static miniufo.basic.InterpolationModel.bicubicLagrangeInterpolation;
import static miniufo.basic.InterpolationModel.bicubicPolynomialInterpolation;
import static miniufo.basic.InterpolationModel.bilinearInterpolation;


/**
 * used to read the commonest csm data file, including interpolation
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CsmDataReadStream implements DataRead,Print{
	//
	private boolean print=true;
	
	private Variable buf=null;
	
	private CsmDescriptor cd=null;
	
	private CtlDataReadStream cdrs=null;
	
	
	/**
     * constructor
     *
     * @param	cd	csm descriptor
     */
	public CsmDataReadStream(CsmDescriptor cd){
		this.cd=cd;
		
		if(!cd.hasCtl() ) throw new IllegalArgumentException("no corresponding ctl file");
		if(!cd.hasData()) throw new IllegalArgumentException("no corresponding data file");
		
		CtlDescriptor ctl=cd.getCtlDescriptor();
		
		cdrs=new CtlDataReadStream(ctl);
		
		if(!ctl.zLinear()) System.out.println("zdef is not linear in ctl");
			//throw new IllegalArgumentException("zdef is not linear in ctl");
		
		if(cd.getDZDef()[0]!=ctl.getDZDef()[0])
			throw new IllegalArgumentException("delta zdef not same in csm and ctl");
		
		if(cd.getDTDef()[0]!=ctl.getDTDef()[0])
			throw new IllegalArgumentException("delta tdef not same in csm and ctl");
		
		if(cd.getZCount()>ctl.getZCount())
			throw new IllegalArgumentException("z count beyond limit");
		
		if(cd.getTCount()>ctl.getTCount())
			throw new IllegalArgumentException("T count beyond limit");
		
		boolean b=false;
		for(int k=0;k<ctl.getZCount();k++)
		if(cd.getZDef().getSamples()[0]==ctl.getZDef().getSamples()[k]){ b=true; break;}
		if(!b) throw new IllegalArgumentException("zdef in csm is not exactly defined in ctl");
		
		b=false;
		for(int l=0;l<ctl.getTCount();l++)
		if(cd.getTimes()[0]==ctl.getTimes()[l]){ b=true; break;}
		if(!b) throw new IllegalArgumentException("tdef in csm is not exactly defined in ctl");
	}
	
	
	/**
	 * to read data from the specified file
	 *
     * @param	v	variable need to fill with data
     */ 
	public void readData(Variable... v){ readData(Type.CUBIC_P,v);}
	
	public void readData(Type type,Variable... v){
		if(print) System.out.print("\nStart reading ");
		
		initialBuffer(v[0].getName(),v[0].isTFirst());
		
		if(print) System.out.print(buf.getName()+" ");	readOne(v[0],type);
		
		for(int m=1;m<v.length;m++){
			if(!v[m].isAreaLike(v[m-1])) throw new IllegalArgumentException("dimension not same");
			
			buf.setName(v[m].getName());
			
			if(print) System.out.print(buf.getName()+" ");
			
			readOne(v[m],type);
		}
		
		if(print) System.out.println("data...\nFinish reading data.");
	}
	
	
	/**
	 * whether to print out
	 *
     * @param	print	print or disable print
     */ 
	public void setPrinting(boolean print){ this.print=print;}
	
	
	private void readOne(Variable v,Type type){
		int t=v.getTCount(),z=v.getZCount(),y=v.getYCount();
		
		CtlDescriptor ctl=cd.getCtlDescriptor();
		
		v.setUndef(ctl.getUndef(v.getName()));
		v.setComment(ctl.getVarComment(v.getName()));
		v.setUnit(ctl.getVarUnit(v.getName()));
		
		int[] trange=buf.getRange().getTRange();	int[] yrange=buf.getRange().getYRange();
		int[] zrange=buf.getRange().getZRange();	int[] xrange=buf.getRange().getXRange();
		
		for(int l=0;l<t;l++){
			/*** set trange ***/
			int vtstart=v.getRange().getTRange()[0];
			trange[0]=trange[1]=ctl.getTNum(cd.getTDef().getSamples()[vtstart-1+l])+1;
			
			/*** set x-y range ***/
			xrange[0]=ctl.getXLENum(ArrayUtil.getMin(cd.getLon()[vtstart-1+l][y-1]))-1;
			yrange[0]=ctl.getYLENum(ArrayUtil.getMin(cd.getLat()[vtstart-1+l][y-1]))-1;
			
			xrange[1]=xrange[2]+xrange[0]-1;
			yrange[1]=yrange[2]+yrange[0]-1;
			
			if(xrange[0]<0||xrange[1]>ctl.getXCount())
				throw new IllegalArgumentException("csm is beyond ctl xrange, olon = "+cd.getOLon()[l]);
			
			if(yrange[0]<0||yrange[1]>ctl.getYCount())
				throw new IllegalArgumentException("csm is beyond ctl yrange, olat = "+cd.getOLat()[l]);
			
			for(int k=0;k<z;k++){
				/*** set zrange ***/
				zrange[0]=zrange[1]=1+ctl.getZNum(
					cd.getZDef().getSamples()[v.getRange().getZRange()[0]-1+k]
				);
				
				/*** read buf ***/
				cdrs.readData(buf);
				
				/*** interpolation ***/
				switch(type){
				case CUBIC_P: interpolateCubicP(v,cd.getLon()[vtstart-1+l],cd.getLat()[vtstart-1+l],k,l); break;
				case CUBIC_L: interpolateCubicL(v,cd.getLon()[vtstart-1+l],cd.getLat()[vtstart-1+l],k,l); break;
				case LINEAR : interpolateLinear(v,cd.getLon()[vtstart-1+l],cd.getLat()[vtstart-1+l],k,l); break;
				default: throw new IllegalArgumentException("unsupported interpolation type");
				}
			}
		}
	}
	
	private void interpolateCubicP(Variable v,float[][] lon,float[][] lat,int k,int l){
		int x=v.getXCount(),xoffset=buf.getRange().getXRange()[0];
		int y=v.getYCount(),yoffset=buf.getRange().getYRange()[0];
		
		CtlDescriptor ctl=cd.getCtlDescriptor();
		
		float undef=ctl.getUndef(v.getName());
		
		float[] dlon=ctl.getDXDef();	float[] ydef=ctl.getYDef().getSamples();
		float[] dlat=ctl.getDYDef();	float[] xdef=ctl.getXDef().getSamples();
		
		float[][] bdata=buf.getData()[0][0];
		
		if(v.isTFirst()){
			float[][] vdata=v.getData()[l][k];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				int xtag=ctl.getXLENum(lon[j][i]);
				int ytag=ctl.getYLENum(lat[j][i]);
				
				float disx=lon[j][i]-xdef[xtag];	xtag-=xoffset;
				float disy=lat[j][i]-ydef[ytag];	ytag-=yoffset;
				
				vdata[j][i]=bicubicPolynomialInterpolation(
					bdata[ytag  ][xtag],bdata[ytag  ][xtag+1],bdata[ytag  ][xtag+2],bdata[ytag  ][xtag+3],
					bdata[ytag+1][xtag],bdata[ytag+1][xtag+1],bdata[ytag+1][xtag+2],bdata[ytag+1][xtag+3],
					bdata[ytag+2][xtag],bdata[ytag+2][xtag+1],bdata[ytag+2][xtag+2],bdata[ytag+2][xtag+3],
					bdata[ytag+3][xtag],bdata[ytag+3][xtag+1],bdata[ytag+3][xtag+2],bdata[ytag+3][xtag+3],
					disx/dlon[xtag],disy/dlat[ytag],undef
				);
			}
			
		}else{
			float[][][] vdata=v.getData()[k];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				int xtag=ctl.getXLENum(lon[j][i]);
				int ytag=ctl.getYLENum(lat[j][i]);
				
				float disx=lon[j][i]-xdef[xtag];	xtag-=xoffset;
				float disy=lat[j][i]-ydef[ytag];	ytag-=yoffset;
				
				vdata[j][i][l]=bicubicPolynomialInterpolation(
					bdata[ytag  ][xtag],bdata[ytag  ][xtag+1],bdata[ytag  ][xtag+2],bdata[ytag  ][xtag+3],
					bdata[ytag+1][xtag],bdata[ytag+1][xtag+1],bdata[ytag+1][xtag+2],bdata[ytag+1][xtag+3],
					bdata[ytag+2][xtag],bdata[ytag+2][xtag+1],bdata[ytag+2][xtag+2],bdata[ytag+2][xtag+3],
					bdata[ytag+3][xtag],bdata[ytag+3][xtag+1],bdata[ytag+3][xtag+2],bdata[ytag+3][xtag+3],
					disx/dlon[xtag],disy/dlat[ytag],undef
				);
			}
		}
	}
	
	private void interpolateCubicL(Variable v,float[][] lon,float[][] lat,int k,int l){
		int x=v.getXCount(),xoffset=buf.getRange().getXRange()[0];
		int y=v.getYCount(),yoffset=buf.getRange().getYRange()[0];
		
		CtlDescriptor ctl=cd.getCtlDescriptor();
		
		float undef=ctl.getUndef(v.getName());
		
		float[] dlon=ctl.getDXDef();	float[] ydef=ctl.getYDef().getSamples();
		float[] dlat=ctl.getDYDef();	float[] xdef=ctl.getXDef().getSamples();
		
		float[][] bdata=buf.getData()[0][0];
		
		if(v.isTFirst()){
			float[][] vdata=v.getData()[l][k];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				int xtag=ctl.getXLENum(lon[j][i]);
				int ytag=ctl.getYLENum(lat[j][i]);
				
				float disx=lon[j][i]-xdef[xtag];	xtag-=xoffset;
				float disy=lat[j][i]-ydef[ytag];	ytag-=yoffset;
				
				vdata[j][i]=bicubicLagrangeInterpolation(
					bdata[ytag  ][xtag],bdata[ytag  ][xtag+1],bdata[ytag  ][xtag+2],bdata[ytag  ][xtag+3],
					bdata[ytag+1][xtag],bdata[ytag+1][xtag+1],bdata[ytag+1][xtag+2],bdata[ytag+1][xtag+3],
					bdata[ytag+2][xtag],bdata[ytag+2][xtag+1],bdata[ytag+2][xtag+2],bdata[ytag+2][xtag+3],
					bdata[ytag+3][xtag],bdata[ytag+3][xtag+1],bdata[ytag+3][xtag+2],bdata[ytag+3][xtag+3],
					1f+disx/dlon[xtag],1f+disy/dlat[ytag],undef
				);
			}
			
		}else{
			float[][][] vdata=v.getData()[k];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				int xtag=ctl.getXLENum(lon[j][i]);
				int ytag=ctl.getYLENum(lat[j][i]);
				
				float disx=lon[j][i]-xdef[xtag];	xtag-=xoffset;
				float disy=lat[j][i]-ydef[ytag];	ytag-=yoffset;
				
				vdata[j][i][l]=bicubicLagrangeInterpolation(
					bdata[ytag  ][xtag],bdata[ytag  ][xtag+1],bdata[ytag  ][xtag+2],bdata[ytag  ][xtag+3],
					bdata[ytag+1][xtag],bdata[ytag+1][xtag+1],bdata[ytag+1][xtag+2],bdata[ytag+1][xtag+3],
					bdata[ytag+2][xtag],bdata[ytag+2][xtag+1],bdata[ytag+2][xtag+2],bdata[ytag+2][xtag+3],
					bdata[ytag+3][xtag],bdata[ytag+3][xtag+1],bdata[ytag+3][xtag+2],bdata[ytag+3][xtag+3],
					1f+disx/dlon[xtag],1f+disy/dlat[ytag],undef
				);
			}
		}
	}
	
	private void interpolateLinear(Variable v,float[][] lon,float[][] lat,int k,int l){
		int x=v.getXCount(),xoffset=buf.getRange().getXRange()[0];
		int y=v.getYCount(),yoffset=buf.getRange().getYRange()[0];
		
		CtlDescriptor ctl=cd.getCtlDescriptor();
		
		float undef=ctl.getUndef(v.getName());
		
		float[] dlon=ctl.getDXDef();	float[] ydef=ctl.getYDef().getSamples();
		float[] dlat=ctl.getDYDef();	float[] xdef=ctl.getXDef().getSamples();
		
		float[][] bdata=buf.getData()[0][0];
		
		if(v.isTFirst()){
			float[][] vdata=v.getData()[l][k];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				int xtag=ctl.getXLENum(lon[j][i]);
				int ytag=ctl.getYLENum(lat[j][i]);
				
				float disx=lon[j][i]-xdef[xtag];	xtag-=xoffset;
				float disy=lat[j][i]-ydef[ytag];	ytag-=yoffset;
				
				vdata[j][i]=bilinearInterpolation(
					bdata[ytag+1][xtag+1],bdata[ytag+1][xtag+2],
					bdata[ytag+2][xtag+1],bdata[ytag+2][xtag+2],
					disx/dlon[xtag],disy/dlat[ytag],undef
				);
			}
			
		}else{
			float[][][] vdata=v.getData()[k];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				int xtag=ctl.getXLENum(lon[j][i]);
				int ytag=ctl.getYLENum(lat[j][i]);
				
				float disx=lon[j][i]-xdef[xtag];	xtag-=xoffset;
				float disy=lat[j][i]-ydef[ytag];	ytag-=yoffset;
				
				vdata[j][i][l]=bilinearInterpolation(
					bdata[ytag+1][xtag+1],bdata[ytag+1][xtag+2],
					bdata[ytag+2][xtag+1],bdata[ytag+2][xtag+2],
					disx/dlon[xtag],disy/dlat[ytag],undef
				);
			}
		}
	}
	
	private void initialBuffer(String vname,boolean tFirst){
		CtlDescriptor ctl=cd.getCtlDescriptor();
		
		// find the northernmost time
		float maxLat=ArrayUtil.getMax(cd.getLat());
		
		int maxLatT=0;
		for(int l=0,L=cd.getTCount();l<L;l++)
		if(maxLat==cd.getLat()[l][cd.getYCount()-1][0]){
			maxLatT=l;	break;
		}
		
		float[][] lons=cd.getLon()[maxLatT];
		float[][] lats=cd.getLat()[maxLatT];
		
		int ymax=ctl.getYNum  (ArrayUtil.getMax(lats));
		int ymin=ctl.getYLENum(ArrayUtil.getMin(lats));
		int xmax=ctl.getXNum  (ArrayUtil.getMax(lons));
		int xmin=ctl.getXLENum(ArrayUtil.getMin(lons));
		
		if(ymax==-1) throw new IllegalArgumentException(ArrayUtil.getMax(lats)+" out of Y-range ["+ctl.getYDef().getFirst()+","+ctl.getYDef().getLast()+"]");
		if(ymin==-1) throw new IllegalArgumentException(ArrayUtil.getMin(lats)+" out of Y-range ["+ctl.getYDef().getFirst()+","+ctl.getYDef().getLast()+"]");
		if(xmax==-1) throw new IllegalArgumentException(ArrayUtil.getMax(lons)+" out of X-range ["+ctl.getXDef().getFirst()+","+ctl.getXDef().getLast()+"]");
		if(xmin==-1) throw new IllegalArgumentException(ArrayUtil.getMin(lons)+" out of X-range ["+ctl.getXDef().getFirst()+","+ctl.getXDef().getLast()+"]");
		
		Range brange=new Range(1,1,ymax-ymin+6,xmax-xmin+6);
		
		buf=new Variable(vname,true,brange);
	}
	
	
	/**
	 * close file method
     */
	public void closeFile(){
		if(cdrs!=null) cdrs.closeFile();
		
		cdrs=null;	buf=null;	cd=null;
	}
}
