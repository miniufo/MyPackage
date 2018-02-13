/**
 * @(#)SphericalSpatialModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.SpatialCoordinate;
import miniufo.descriptor.TemporalCoordinate;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.tan;


/**
 * Spatial model in spherical-polar coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class SphericalSpatialModel extends SpatialModel{
	//
	private float[] f1  =null;	// 2*omega*sin(fai)
	private float[] f2  =null;	// 2*omega*cos(fai)
	private float[] Beta=null;	// Rossby arg, f2/r
	private float[] rcos=null;	// R * cos(lat)
	private float[] lsin=null;	// sin(lat)
	private float[] lcos=null;	// cos(lat)
	private float[] ltan=null;	// tan(lat)
	
	
	/**
     * constructor
     *
     * @param	pr	a given parser
     */
	public SphericalSpatialModel(DataDescriptor dd){ _SphericalSpatialModel(dd);}
	
	public SphericalSpatialModel(DataDescriptor dd,LevelType levType){
		this.levType=levType;
		_SphericalSpatialModel(dd);
	}
	
	
	/*** getor and setor ***/
	public float[] getF1(){ return f1;}
	
	public float[] getF2(){ return f2;}
	
	public float[] getRCos() { return rcos;}
	
	public float[] getLSin() { return lsin;}
	
	public float[] getLCos() { return lcos;}
	
	public float[] getLTan() { return ltan;}
	
	public float[] getBeta() { return Beta;}
	
	
	/**
     * Whether the given model is global.
     */
	public boolean isGlobal(){
		if(!periodicX) return false;
		
		double dlon=Math.toDegrees(xdef.getIncrements()[0]);
		double dlat=Math.toDegrees(ydef.getIncrements()[0]);
		
		if(Math.abs(360.0/xdef.length()-dlon)/dlon>1e-2f) return false;
		if(Math.abs(Math.toDegrees(ydef.getFirst())+90)>=dlat) return false;
		if(Math.abs(Math.toDegrees(ydef.getLast ())-90)>=dlat) return false;
		
		return true;
	}
	
	
	/**
     * used to print out
     */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		int tcount=tdef.length(),zcount=zdef.length();
		int ycount=ydef.length(),xcount=xdef.length();
		
		sb.append("\nSpherical Spatial Model:\n");
		
		sb.append("tdef:\n  from: ");
		sb.append(tdef.getSamples()[0]);		sb.append("\n  to  : ");
		sb.append(tdef.getSamples()[tcount-1]);	sb.append("\ttcount: ");		sb.append(tcount);
		
		sb.append("\nzdef:\n  from: ");
		sb.append(zdef.getSamples()[0]/100f);		sb.append("\n  to  : ");
		sb.append(zdef.getSamples()[zcount-1]/100f);sb.append("\tzcount: ");	sb.append(zcount);
		
		sb.append("\nydef:\n  from: ");
		sb.append((float)(ydef.getSamples()[0]*180/Math.PI));			sb.append("\n  to  : ");
		sb.append((float)(ydef.getSamples()[ycount-1]*180/Math.PI));	sb.append("\tycount: ");sb.append(ycount);
		
		sb.append("\nxdef:\n  from: ");
		sb.append((float)(xdef.getSamples()[0]*180/Math.PI));			sb.append("\n  to  : ");
		sb.append((float)(xdef.getSamples()[xcount-1]*180/Math.PI));	sb.append("\txcount: ");sb.append(xcount);
		
		return sb.toString();
	}
	
	
	/*** helper methods ***/
	private void _SphericalSpatialModel(DataDescriptor dd){
		if(dd instanceof CsmDescriptor||dd instanceof CtsDescriptor)
		throw new IllegalArgumentException("DataDescriptor cannot be CsmDescriptor and CtsDescriptor");
		
		tdef=(TemporalCoordinate)(dd.getTDef().clone());
		zdef=( SpatialCoordinate)(dd.getZDef().clone()); int zcount=zdef.length();
		ydef=( SpatialCoordinate)(dd.getYDef().clone()); int ycount=ydef.length();
		xdef=( SpatialCoordinate)(dd.getXDef().clone()); int xcount=xdef.length();
		
		f1  =new float[ycount];		rcos=new float[ycount];
		f2  =new float[ycount];		dxs =new float[ycount];
		Beta=new float[ycount];		lcos=new float[ycount];
		lsin=new float[ycount];		ltan=new float[ycount];
		
		/*** z-direction ***/
		if(levType==LevelType.PRESSURE)
		for(int k=0,K=zcount;k<K;k++) zdef.getSamples()[k]*=100;	// change to Pa
		
		/*** y-direction ***/
		for(int j=0;j<ycount;j++){
			float tmp=ydef.getSamples()[j]*=(float)(Math.PI/180);	// to radian
			lsin[j]=(float)sin(tmp);
			lcos[j]=(float)cos(tmp);
			ltan[j]=(float)tan(tmp);
			rcos[j]=EARTH_RADIUS*lcos[j];
			  f1[j]=2*EARTH_ROTATE_SPEED*lsin[j];
			  f2[j]=2*EARTH_ROTATE_SPEED*lcos[j];
			Beta[j]=f2[j]/EARTH_RADIUS;
		}
		
		if(levType==LevelType.PRESSURE) dz=dd.getDZDef()[0]*100;	// change to Pa
		else dz=dd.getDZDef()[0];
		
		dt=dd.getDTDef()[0];
		
		if(!ydef.isLinear()&&ycount!=1)
			dy=EARTH_RADIUS*ydef.getRange()/(ycount-1);
		else
			dy=EARTH_RADIUS*(float)(dd.getDYDef()[0]*Math.PI/180);
		
		if(!xdef.isLinear()&&xcount!=1)
			dx=EARTH_RADIUS*xdef.getRange()/(xcount-1);
		else
			dx=EARTH_RADIUS*(float)(dd.getDXDef()[0]*Math.PI/180);
		
		/*** x-direction ***/
		for(int i=0;i<xcount;i++) xdef.getSamples()[i]*=Math.PI/180;
		if(!xdef.isLinear()&&xcount!=1){
			for(int j=0;j<ycount;j++){
				dxs[j]=EARTH_RADIUS*lcos[j]*xdef.getRange()/(xcount-1);
				dxs[j]=dxs[j]<0?0:dxs[j];
			}
			
		}else for(int j=0;j<ycount;j++){
			dxs[j]=EARTH_RADIUS*lcos[j]*(float)(dd.getDXDef()[0]*Math.PI/180);
			dxs[j]=dxs[j]<0?0:dxs[j];
		}
		
		this.periodicX=dd.isPeriodicX();
		this.periodicY=dd.isPeriodicY();
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			SpheralSpacialModel ssm=new SpheralSpacialModel(
				new miniufo.descriptor.CtlDescriptor("d:/data/typhoon/catarina/catarina.ctl")
			);
			
			System.out.println(ssm);
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
