/**
 * @(#)temporalCoordinate.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import java.util.Arrays;

import miniufo.diagnosis.MDate;


/**
 * discrete sampled temporal coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TemporalCoordinate extends Coordinate{
	//
	private MDate[] samples=null;
	
	
	/**
     * constructor
     *
     * param	samples		an array of samples
     */
	public TemporalCoordinate(String name,MDate[] samples){
		super(name);
		
		if(samples==null) throw new NullPointerException("samples are null");
		
		int length=samples.length;
		
		if(length<1) throw new IllegalArgumentException("sample length should be positive");
		
		this.samples=samples;
		
		for(int i=1,I=length-1;i<I;i++)
		if(samples[i].getDT(samples[i-1])!=samples[i].getDT(samples[i+1])){ islinear=false; break;}
	}
	
	public TemporalCoordinate(String name,MDate[] samples,boolean islinear){
		super(name);
		
		if(samples==null) throw new NullPointerException("samples are null");
		
		int length=samples.length;
		
		if(length<1) throw new IllegalArgumentException("sample length should be positive");
		
		this.samples=samples;	this.islinear=islinear;
	}
	
	
	/*** getor and setor ***/
	public int length(){ return samples.length;}
	
	public boolean isLike(Coordinate sc){
		if(!(sc instanceof TemporalCoordinate)) return false;
		if(islinear!=sc.islinear) return false;
		if(!name.equals(sc.name)) return false;
		if(!Arrays.equals(samples,((TemporalCoordinate)sc).samples)) return false;
		
		return true;
	}
	
	public MDate getFirst(){ return samples[0];}
	
	public MDate getLast(){ return samples[samples.length-1];}
	
	public  long[] getLongSamples(){
		long[] rs=new long[samples.length];
		
		for(int l=0,L=samples.length;l<L;l++) rs[l]=samples[l].getLongTime();
		
		return rs;
	}
	
	public int[] getIncrements(){
		int len=samples.length-1;
		
		if(len==0) throw new IllegalArgumentException("sample length is 1");
		
		int[] ds=new int[len];
		
		for(int i=0;i<len;i++) ds[i]=samples[i+1].getDT(samples[i]);
		
		return ds;
	}
	
	public MDate[] getSamples(){ return samples ;}
	
	
	/**
     * clone method
     */
	public TemporalCoordinate clone(){
		TemporalCoordinate cd=(TemporalCoordinate)super.clone();
		
		cd.samples=samples.clone();
		
		return cd;
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			
			
		}catch(Exception e){ e.printStackTrace();}
	}*/
}
