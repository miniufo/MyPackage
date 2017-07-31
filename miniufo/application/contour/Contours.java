/**
 * @(#)Contours.java	1.0 2017.06.09
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import miniufo.basic.ArrayUtil;
import static miniufo.application.statisticsModel.StatisticsApplication.cArithmeticMean;


/**
 * This class describe a series of contours used to setup a
 * tracer-contour-based coordinate space.
 *
 * @version 1.0, 2017.06.09
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Contours{
	//
	private boolean increSToN=true;	// increasing from south to north
	private boolean isLinear =true;	// equally-spaced contours
	
	private int      number=0;		// number of contours
	
	private double southVal=0;		// contour value at south boundary
	private double northVal=0;		// contour value at north boundary
	private double meanIntv=0;		// mean interval
	
	private double[] values=null;	// values of contours
	private double[] areas =null;	// area of each contour
	private double[] Yeqvs =null;	// equivalent Ys
	
	
	/**
	 * Compute the contours given the 2D array data and number of contours.
	 * Whether contours are increasing from south to north is determined by
	 * the data themselves.
	 * 
	 * @param	data	an 2D array data
	 * @param	undef	undefined value
	 * @param	numOfC	number of contours
	 */
	public Contours(float[][] data,float undef,int numOfC){
		float[] ex=ArrayUtil.getExtrema(data,undef);
		
		increSToN=cArithmeticMean(data[0],undef)<cArithmeticMean(data[data.length-1],undef);
		
		if(increSToN){ southVal=ex[0]; northVal=ex[1];}
		else         { southVal=ex[1]; northVal=ex[0];}
		
		number=numOfC;
		meanIntv=(northVal-southVal)/(number-1);
		
		check();
		
		values=new double[number];
		for(int c=0,I=number;c<I;c++) values[c]=southVal+meanIntv*c;
		
		isLinear=true;
	}
	
	/**
	 * Compute the contours given the south-most and north-most contours
	 * 
	 * @param	csouth	south-most contour
	 * @param	cnorth	north-most contour
	 * @param	numOfC	number of contours
	 */
	public Contours(float csouth,float cnorth,float inc){
		increSToN=csouth<cnorth;
		southVal =csouth;
		northVal =cnorth;
		number   =(int)Math.round((northVal-southVal)/inc)+1;
		meanIntv =inc;
		
		if(number-1!=Math.round((northVal-southVal)/inc))
		throw new IllegalArgumentException("overflow");
		
		check();
		
		values=new double[number];
		for(int c=0,I=number;c<I;c++) values[c]=southVal+inc*c;
		
		isLinear=true;
	}
	
	/**
	 * Compute the contours by specifying all the values
	 * from south to north
	 * 
	 * @param	cs	contours from south to north
	 */
	public Contours(float[] cs){
		increSToN=cs[0]<cs[cs.length-1];
		southVal =cs[0];
		northVal =cs[cs.length-1];
		number   =cs.length;
		meanIntv =(northVal-southVal)/(number-1.0);
		
		check();
		
		values=new double[number];
		for(int i=0,I=values.length;i<I;i++) values[i]=cs[i];
		
		for(int i=0,I=number-1;i<I;i++)
		if(Math.abs(values[i+1]-values[i]-meanIntv)/meanIntv>1e-3){ isLinear=false; break;}
	}
	
	/**
	 * Compute the contours by specifying all the values
	 * from south to north
	 * 
	 * @param	cs	contours from south to north
	 */
	public Contours(double[] cs){
		increSToN=cs[0]<cs[cs.length-1];
		southVal =cs[0];
		northVal =cs[cs.length-1];
		number   =cs.length;
		meanIntv =(northVal-southVal)/(number-1.0);
		
		check();
		
		values=cs.clone();
		
		for(int i=0,I=number-1;i<I;i++)
		if(Math.abs(values[i+1]-values[i]-meanIntv)/meanIntv>1e-3){ isLinear=false; break;}
	}
	
	
	/*** getor and setor ***/
	public boolean increaseSToN(){ return increSToN;}
	
	public boolean isLinear(){ return isLinear;}
	
	public int getContourNumber(){ return number;}
	
	public double getMaxValue(){ return increSToN?northVal:southVal;}
	
	public double getMinValue(){ return increSToN?southVal:northVal;}
	
	public double[] getEquivalentYs(){ return Yeqvs;}
	
	public double[] getValues(){ return values;}
	
	public double[] getAreas(){ return areas;}
	
	public void setAreas(double[] areas){ this.areas=areas;}
	
	public void setYEs(double[] Ye){ this.Yeqvs=Ye;}
	
	
	/**
     * used to print out
     */
	public String toString(){
		if(isLinear)
			return "contours are linear "+"["+southVal+" : " +meanIntv+" : "+northVal+"]";
		else
			return "contours are levels "+"["+southVal+" : ~"+meanIntv+" : "+northVal+"]";
	}
	
	
	/*** helper methods ***/
	private void check(){
		if(number<2)
		throw new IllegalArgumentException("the number of contour is at least 2");
		if(meanIntv==0)
		throw new IllegalArgumentException("mean increment is 0");
	}
	
	
	/*** test
	public static void main(String[] args){
		
	}*/
}
