/**
 * @(#)SpatialCoordinate.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import java.util.Arrays;

import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel;
import miniufo.basic.InterpolationModel.Type;


/**
 * discrete sampled spatial coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class SpatialCoordinate extends Coordinate{
	//
	private float min=0;
	private float max=0;
	
	private float[] samples=null;
	
	
	/**
     * constructor
     *
     * param	samples		an array of samples
     */
	public SpatialCoordinate(String name,float[] samples){
		super(name);
		
		if(samples==null) throw new NullPointerException("samples are null");
		
		int length=samples.length;
		
		if(length<1) throw new IllegalArgumentException("sample length should be positive");
		else if(length>1&&samples[length-1]<samples[0]) isIncre=false;
		
		this.name=name;	this.samples=samples;
		
		float[] ex=ArrayUtil.getExtrema(this.samples);
		
		min=ex[0];	max=ex[1];
		
		if(length>2){
			float aveInc=(max-min)/(length-1);
			
			for(int i=1;i<length;i++)
			if(Math.abs(samples[i]-samples[i-1])/aveInc-1>5e-4f){ islinear=false; break;}
		}
	}
	
	
	/*** getor and setor ***/
	public int length(){ return samples.length;}
	
	public float getMin(){ return min;}
	
	public float getMax(){ return max;}
	
	public float getFirst(){ return samples[0];}
	
	public float getLast(){ return samples[samples.length-1];}
	
	public float getRange(){ return getLast()-getFirst();}
	
	public boolean isLike(Coordinate sc){
		if(!(sc instanceof SpatialCoordinate)) return false;
		if(islinear!=sc.islinear) return false;
		if(!name.equals(sc.name)) return false;
		if(!Arrays.equals(samples,((SpatialCoordinate)sc).samples)) return false;
		
		return true;
	}
	
	public float[] getSamples(){ return samples;}
	
	public float[] getIncrements(){
		int len=samples.length-1;
		
		if(len==0) return new float[]{1};
		
		float[] ds=new float[len];
		
		for(int i=0;i<len;i++) ds[i]=samples[i+1]-samples[i];
		
		return ds;
	}
	
	
	/**
     * clone method
     */
	public SpatialCoordinate clone(){
		SpatialCoordinate scd=(SpatialCoordinate)super.clone();
		
		scd.samples=samples.clone();
		
		return scd;
	}
	
	
	/**
     * Re-sample (interpolate) the coordinate using equally-spaced n-points
     * 
     * @param	n		equally-spaced n-points
     * @param	type	type of interpolation
     */
	public SpatialCoordinate resample(int n,Type type){
		if(n<=1) throw new IllegalArgumentException("n should be at least 2");
		
		float[] spl=InterpolationModel.interp1D(samples,n,type);
		
		if(type==Type.PERIODIC_LINEAR){
			int idx=0;
			
			float last=getLast();
			float incr=spl[1]-spl[0];
			
			for(int i=0,I=spl.length;i<I;i++) if(spl[i]+incr>last){ idx=i; break;}
			
			for(int i=idx+1,ii=1,I=spl.length;i<I;i++) spl[i]=spl[i-1]+ii*incr;
			
		}else if(type==Type.PERIODIC_CUBIC_L||type==Type.PERIODIC_CUBIC_P||type==Type.PERIODIC_SPLINE){
			throw new IllegalArgumentException("unsupported type for resampling");
		}
		
		return new SpatialCoordinate(name,spl);
	}
	
	
	/** test
	public static void main(String[] args){
		float[] orig=new float[]{0,1,2,3,4,5,6,7,8,9};
		
		SpatialCoordinate sc=new SpatialCoordinate("xdim",orig);
		
		SpatialCoordinate interp=sc.resample(19,Type.PERIODIC_LINEAR);
		
		System.out.println(Arrays.toString(interp.getSamples()));
	}*/
}
