/**
 * @(#)CoordinateTransformation.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.advanced;

import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.CylindricalSpatialModel;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static miniufo.basic.InterpolationModel.bicubicPolynomialInterpolation;


/**
 * transformation between cylindrical and spherical coordinates
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CoordinateTransformation{
	//
	private   SphericalSpatialModel ssm=null;
	private CylindricalSpatialModel csm=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		spacial model in spheral coordinate
     * @param	csm		spacial model in cylindrical coordinate
     */
	public CoordinateTransformation(SphericalSpatialModel ssm,CylindricalSpatialModel csm){
		this.ssm=ssm;	this.csm=csm;
	}
	
	
	/**
     * reproject a given variable (usually a vector)
     * from lat/lon coordinates to cylindrical coordinates
     *
     * @param	vx		a given Variable in x direction
     * @param	vy		a given Variable in y direction
     *
     * @return	result (usually a vector) of the transformation, [0] is tangential and [1] is radial
     */
	public Variable[] reprojectToCylindrical(Variable vx,Variable vy){
		if(!vx.isLike(vy)) throw new IllegalArgumentException();
		
		int t=vx.getTCount(),	z=vx.getZCount(),	y=csm.getYCount(),	x=csm.getXCount();
		
		float undef=vx.getUndef();
		float[][][] eta=csm.getEta();
		Variable[]   nv=new Variable[2];
		
		nv[0]=new Variable("ut",vx.isTFirst(),new Range(t,z,y,x));	nv[0].setUndef(undef);
		nv[1]=new Variable("vr",vx.isTFirst(),new Range(t,z,y,x));	nv[1].setUndef(undef);
		
		nv[0].setCommentAndUnit("tangential velocity (m s^-1)");
		nv[1].setCommentAndUnit("radial velocity (m s^-1)");
		
		Range r0=nv[0].getRange();
		Range r1=nv[1].getRange();
		Range ur=   vx.getRange();
		
		r0.setTRange(ur);	r0.setZRange(ur);	r0.setYRange(ur);	r0.setXRange(ur);
		r1.setTRange(ur);	r1.setZRange(ur);	r1.setYRange(ur);	r1.setXRange(ur);
		
		float[][][][]  v1data=vx.getData();
		float[][][][]  v2data=vy.getData();
		float[][][][] nv1data=nv[0].getData();
		float[][][][] nv2data=nv[1].getData();
		
		if(vx.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(v1data[l][k][j][i]!=undef&&v2data[l][k][j][i]!=undef){
					nv1data[l][k][j][i]=(float)(
						-v1data[l][k][j][i]*cos(eta[l][j][i])-v2data[l][k][j][i]*sin(eta[l][j][i])
					);
					
					nv2data[l][k][j][i]=(float)(
						-v1data[l][k][j][i]*sin(eta[l][j][i])+v2data[l][k][j][i]*cos(eta[l][j][i])
					);
					
				}else{
					nv1data[l][k][j][i]=undef;
					nv2data[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(v1data[k][j][i][l]!=undef&&v2data[k][j][i][l]!=undef){
					nv1data[k][j][i][l]=(float)(
						-v1data[k][j][i][l]*cos(eta[l][j][i])-v2data[k][j][i][l]*sin(eta[l][j][i])
					);
					
					nv2data[k][j][i][l]=(float)(
						-v1data[k][j][i][l]*sin(eta[l][j][i])+v2data[k][j][i][l]*cos(eta[l][j][i])
					);
					
				}else{
					nv1data[k][j][i][l]=undef;
					nv2data[k][j][i][l]=undef;
				}
			}
		}
		
		return nv;
	}
	
	
	/**
     * reproject a given variable (usually a vector)
     * from cylindrical coordinates to lat/lon coordinates
     *
     * @param	ut		a given Variable in tangential direction
     * @param	vr		a given Variable in radial direction
     *
     * @return	result (usually a vector) of the transformation, [0] is lon and [1] is lat
     */
	public Variable[] reprojectToLatLon(Variable ut,Variable vr){
		if(!ut.isLike(vr)) throw new IllegalArgumentException();
		
		int t=ut.getTCount(),	z=ut.getZCount(),	y=csm.getYCount(),	x=csm.getXCount();
		
		float undef=ut.getUndef();
		float[][][] eta=csm.getEta();
		Variable[]   nv=new Variable[2];
		
		nv[0]=new Variable("u",ut.isTFirst(),new Range(t,z,y,x));	nv[0].setUndef(undef);
		nv[1]=new Variable("v",ut.isTFirst(),new Range(t,z,y,x));	nv[1].setUndef(undef);
		
		nv[0].setCommentAndUnit("westerly");
		nv[1].setCommentAndUnit("southerly");
		
		Range r0=nv[0].getRange();
		Range r1=nv[1].getRange();
		Range ur=   ut.getRange();
		
		r0.setTRange(ur);	r0.setZRange(ur);	r0.setYRange(ur);	r0.setXRange(ur);
		r1.setTRange(ur);	r1.setZRange(ur);	r1.setYRange(ur);	r1.setXRange(ur);
		
		float[][][][]  v1data=ut.getData();
		float[][][][]  v2data=vr.getData();
		float[][][][] nv1data=nv[0].getData();
		float[][][][] nv2data=nv[1].getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(v1data[l][k][j][i]!=undef&&v2data[l][k][j][i]!=undef){
					nv1data[l][k][j][i]=(float)(
						-v1data[l][k][j][i]*cos(eta[l][j][i])-v2data[l][k][j][i]*sin(eta[l][j][i])
					);
					
					nv2data[l][k][j][i]=(float)(
						-v1data[l][k][j][i]*sin(eta[l][j][i])+v2data[l][k][j][i]*cos(eta[l][j][i])
					);
					
				}else{
					nv1data[l][k][j][i]=undef;
					nv2data[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(v1data[k][j][i][l]!=undef&&v2data[k][j][i][l]!=undef){
					nv1data[k][j][i][l]=(float)(
						-v1data[k][j][i][l]*cos(eta[l][j][i])-v2data[k][j][i][l]*sin(eta[l][j][i])
					);
					
					nv2data[k][j][i][l]=(float)(
						-v1data[k][j][i][l]*sin(eta[l][j][i])+v2data[k][j][i][l]*cos(eta[l][j][i])
					);
					
				}else{
					nv1data[k][j][i][l]=undef;
					nv2data[k][j][i][l]=undef;
				}
			}
		}
		
		return nv;
	}
	
	
	/**
     * transfered a given variable from lat/lon coordinates to cylindrical coordinates
     *
     * @param	v	a given Variable
     *
     * @return	result in cylindrical coordinates
     */
	public Variable transToCylindrical(Variable v){
		if(!ssm.isGlobal())
			throw new IllegalArgumentException("require a global model");
		
		if(!ssm.isAreaLike(v))
			throw new IllegalArgumentException("require a variable covering the global");
		
		if(isBeyondRange()) throw new IllegalArgumentException("csm is beyond the range of ssm");
		
		int tt=v.getTCount(),yy=csm.getYCount();
		int zz=v.getZCount(),xx=csm.getXCount();
		
		Variable re=new Variable(v.getName(),new Range(tt,zz,yy,xx));
		re.setComment(v.getComment());
		re.setUnit(v.getUnit());
		re.setUndef(v.getUndef());
		
		float[][][] lons=csm.getLon();	// radians
		float[][][] lats=csm.getLat();	// radians
		
		float[] xdef=ssm.getXDef().getSamples();
		float[] ydef=ssm.getYDef().getSamples();
		
		float undef=v.getUndef();
		float dlon=ssm.getXDef().getIncrements()[0];	// radians
		float dlat=ssm.getYDef().getIncrements()[0];	// radians
		
		if(v.isTFirst()){
			for(int l=0;l<tt;l++)
			for(int k=0;k<zz;k++){
				float[][] vdata= v.getData()[l][k];
				float[][] rdata=re.getData()[l][k];
				
				for(int j=0;j<yy;j++)
				for(int i=0;i<xx;i++){
					int xtags=ArrayUtil.getLEIdxIncre(xdef,lons[l][j][i]);
					int ytags=ArrayUtil.getLEIdxIncre(ydef,lats[l][j][i]);
					
					float disx=lons[l][j][i]-xdef[xtags];
					float disy=lats[l][j][i]-ydef[ytags];
					
					rdata[j][i]=bicubicPolynomialInterpolation(
						vdata[ytags-1][xtags-1],vdata[ytags-1][xtags],vdata[ytags-1][xtags+1],vdata[ytags-1][xtags+2],
						vdata[ytags  ][xtags-1],vdata[ytags  ][xtags],vdata[ytags  ][xtags+1],vdata[ytags  ][xtags+2],
						vdata[ytags+1][xtags-1],vdata[ytags+1][xtags],vdata[ytags+1][xtags+1],vdata[ytags+1][xtags+2],
						vdata[ytags+2][xtags-1],vdata[ytags+2][xtags],vdata[ytags+2][xtags+1],vdata[ytags+2][xtags+2],
						disx/dlon,disy/dlat,undef
					);
				}
			}
			
		}else{
			for(int k=0;k<zz;k++){
				float[][][] vdata= v.getData()[k];
				float[][][] rdata=re.getData()[k];
				
				for(int l=0;l<tt;l++)
				for(int j=0;j<yy;j++)
				for(int i=0;i<xx;i++){
					int xtags=ArrayUtil.getLEIdxIncre(xdef,lons[l][j][i]);
					int ytags=ArrayUtil.getLEIdxIncre(ydef,lats[l][j][i]);
					
					float disx=lons[l][j][i]-xdef[xtags];
					float disy=lats[l][j][i]-ydef[ytags];
					
					rdata[j][i][l]=bicubicPolynomialInterpolation(
						vdata[ytags-1][xtags-1][l],vdata[ytags-1][xtags][l],vdata[ytags-1][xtags+1][l],vdata[ytags-1][xtags+2][l],
						vdata[ytags  ][xtags-1][l],vdata[ytags  ][xtags][l],vdata[ytags  ][xtags+1][l],vdata[ytags  ][xtags+2][l],
						vdata[ytags+1][xtags-1][l],vdata[ytags+1][xtags][l],vdata[ytags+1][xtags+1][l],vdata[ytags+1][xtags+2][l],
						vdata[ytags+2][xtags-1][l],vdata[ytags+2][xtags][l],vdata[ytags+2][xtags+1][l],vdata[ytags+2][xtags+2][l],
						disx/dlon,disy/dlat,undef
					);
				}
			}
		}
		
		return re;
	}
	
	/**
     * transfered a given variable from lat/lon coordinates to cylindrical coordinates
     * the given variable should be invariant (tcount==1)
     *
     * @param	v	a given Variable
     *
     * @return	result in cylindrical coordinates
     */
	public Variable transToCylindricalInvariantly(Variable v){
		if(!ssm.isGlobal())
			throw new IllegalArgumentException("require a global model");
		
		if(!ssm.isAreaLike(v))
			throw new IllegalArgumentException("require a variable covering the global");
		
		if(isBeyondRange()) throw new IllegalArgumentException("csm is beyond the range of ssm");
		
		int tt=csm.getTCount(),yy=csm.getYCount();
		int zz=  v.getZCount(),xx=csm.getXCount();
		
		Variable re=new Variable(v.getName(),v.isTFirst(),new Range(tt,zz,yy,xx));
		re.setComment(v.getComment());
		re.setUnit(v.getUnit());
		re.setUndef(v.getUndef());
		
		float[][][] lons=csm.getLon();	// radians
		float[][][] lats=csm.getLat();	// radians
		
		float[] xdef=ssm.getXDef().getSamples();
		float[] ydef=ssm.getYDef().getSamples();
		
		float undef=v.getUndef();
		float dlon=ssm.getXDef().getIncrements()[0];	// radians
		float dlat=ssm.getYDef().getIncrements()[0];	// radians
		
		if(v.isTFirst()){
			for(int l=0;l<tt;l++)
			for(int k=0;k<zz;k++){
				float[][] vdata= v.getData()[0][k];
				float[][] rdata=re.getData()[l][k];
				
				for(int j=0;j<yy;j++)
				for(int i=0;i<xx;i++){
					int xtags=ArrayUtil.getLEIdxIncre(xdef,lons[l][j][i]);
					int ytags=ArrayUtil.getLEIdxIncre(ydef,lats[l][j][i]);
					
					float disx=lons[l][j][i]-xdef[xtags];
					float disy=lats[l][j][i]-ydef[ytags];
					
					rdata[j][i]=bicubicPolynomialInterpolation(
						vdata[ytags-1][xtags-1],vdata[ytags-1][xtags],vdata[ytags-1][xtags+1],vdata[ytags-1][xtags+2],
						vdata[ytags  ][xtags-1],vdata[ytags  ][xtags],vdata[ytags  ][xtags+1],vdata[ytags  ][xtags+2],
						vdata[ytags+1][xtags-1],vdata[ytags+1][xtags],vdata[ytags+1][xtags+1],vdata[ytags+1][xtags+2],
						vdata[ytags+2][xtags-1],vdata[ytags+2][xtags],vdata[ytags+2][xtags+1],vdata[ytags+2][xtags+2],
						disx/dlon,disy/dlat,undef
					);
				}
			}
			
		}else{
			for(int k=0;k<zz;k++){
				float[][][] vdata= v.getData()[k];
				float[][][] rdata=re.getData()[k];
				
				for(int l=0;l<tt;l++)
				for(int j=0;j<yy;j++)
				for(int i=0;i<xx;i++){
					int xtags=ArrayUtil.getLEIdxIncre(xdef,lons[l][j][i]);
					int ytags=ArrayUtil.getLEIdxIncre(ydef,lats[l][j][i]);
					
					float disx=lons[l][j][i]-xdef[xtags];
					float disy=lats[l][j][i]-ydef[ytags];
					
					rdata[j][i][l]=bicubicPolynomialInterpolation(
						vdata[ytags-1][xtags-1][0],vdata[ytags-1][xtags][0],vdata[ytags-1][xtags+1][0],vdata[ytags-1][xtags+2][0],
						vdata[ytags  ][xtags-1][0],vdata[ytags  ][xtags][0],vdata[ytags  ][xtags+1][0],vdata[ytags  ][xtags+2][0],
						vdata[ytags+1][xtags-1][0],vdata[ytags+1][xtags][0],vdata[ytags+1][xtags+1][0],vdata[ytags+1][xtags+2][0],
						vdata[ytags+2][xtags-1][0],vdata[ytags+2][xtags][0],vdata[ytags+2][xtags+1][0],vdata[ytags+2][xtags+2][0],
						disx/dlon,disy/dlat,undef
					);
				}
			}
		}
		
		return re;
	}
	
	
	/**
     * project a vector (u,v) onto along mean vector (um,vm) direction
     *
     * @param	u	zonal component of a vector
     * @param	v	meridional component of a vector
     * @param	um	zonal component of the reference vector
     * @param	vm	meridional component of the reference vector
     *
     * @return	re	projection in the mean vector coordinate
     */
	public static float[] projectToNaturalCoords(float u,float v,float um,float vm){
		float[] re=new float[2];
		
		double arg=Math.atan2(v,u)-Math.atan2(vm,um);
		double mod=Math.hypot(u,v);
		
		re[0]=(float)(mod*Math.cos(arg));
		re[1]=(float)(mod*Math.sin(arg));
		
		return re;
	}
	
	
	/*** helper methods ***/
	
	/**
     * whether the polar coordinate area is beyond the X-Y coordinate area
     *
     * @return	true or false
     */
	private boolean isBeyondRange(){
		int r_tag=csm.getYCount()-1;
		float min_lon=ssm.getXDef().getSamples()[0];
		float max_lon=ssm.getXDef().getSamples()[ssm.getXCount()-1];
		float min_lat=ssm.getYDef().getSamples()[0];
		float max_lat=ssm.getYDef().getSamples()[ssm.getYCount()-1];
		
		for(int l=0;l<csm.getTCount();l++)
		for(int i=0;i<csm.getXCount();i++){
			if(csm.getLon()[l][r_tag][i]<min_lon){
				System.out.println("  warning, invalid minLon: "+
					Math.toDegrees(csm.getLon()[l][r_tag][i])+"<"+
					Math.toDegrees(min_lon)
				);
				return true;
			}
			if(csm.getLon()[l][r_tag][i]>max_lon){
				System.out.println("  warning, invalid maxLon: "+
					Math.toDegrees(csm.getLon()[l][r_tag][i])+">"+
					Math.toDegrees(max_lon)
				);
				return true;
			}
			if(csm.getLat()[l][r_tag][i]<min_lat){
				System.out.println("  warning, invalid minLat: "+
					Math.toDegrees(csm.getLat()[l][r_tag][i])+"<"+
					Math.toDegrees(min_lat)
				);
				return true;
			}
			if(csm.getLat()[l][r_tag][i]>max_lat){
				System.out.println("  warning, invalid maxLat: "+
					Math.toDegrees(csm.getLat()[l][r_tag][i])+">"+
					Math.toDegrees(max_lat)
				);
				return true;
			}
		}
		
		return false;
	}
	
	
	/** test
	public static void main(String[] args){
		Complex co=new Complex(3,7);
		Complex cm=new Complex(5,4);
		
		Complex c=co.divide(cm.divide(cm.getMod()));
		
		float[] re=projectToNaturalCoords(co.getReal(),co.getImag(),cm.getReal(),cm.getImag());
		
		Complex c2=new Complex(re[0],re[1]);
		
		System.out.println(c);
		System.out.println(c2);
	}*/
}
