/**
 * @(#)WindowFunction.java	1.0 2013.03.12
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.statistics.StatisticsUtil;


/**
 * window function
 *
 * @version 1.0, 2013.03.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class WindowFunction{
	
	/**
	 * prevent from instantiate
     */
	private WindowFunction(){}
	
	
	/**
	 * generates the N-point Hann/Hamming window using symmetric window sampling.
	 * Note that the first and last zero-weighted window samples are not included in hanning.
	 * equals to matlab code: w=hann(L)/hanning(L)
     *
     * @param	L	an positive integer
     *
     * @return	w	N-point symmetric Hanning window 
	 */
	public static float[] hann(int L){
		if(L<1) throw new IllegalArgumentException("L should be positive");
		
		float[] w=new float[L];
		
		for(int l=0;l<L;l++)
		w[l]=(float)(0.5*(1.0-Math.cos(2.0*Math.PI*l/(L-1))));
		
		return w;
	}
	
	public static float[] hanning(int L){
		if(L<1) throw new IllegalArgumentException("L should be positive");
		
		float[] w=new float[L];
		
		for(int l=0;l<L;l++)
		w[l]=(float)(0.5*(1.0-Math.cos(2.0*Math.PI*(l+1)/(L+1))));
		
		return w;
	}
	
	public static float[] hamming(int L){
		if(L<1) throw new IllegalArgumentException("L should be positive");
		
		float[] w=new float[L];
		
		for(int l=0;l<L;l++)
		w[l]=(float)(0.54-0.46*Math.cos(2.0*Math.PI*(l)/(L-1)));
		
		return w;
	}
	
	public static float[] rectangular(int L){
		if(L<1) throw new IllegalArgumentException("L should be positive");
		
		float[] w=new float[L];
		
		for(int l=0;l<L;l++) w[l]=1f;
		
		return w;
	}
	
	public static float[] tukey(int L,float R){
		if(L<1) throw new IllegalArgumentException("L should be positive");
		if(R>1||R<0) throw new IllegalArgumentException("R should be in range [0, 1]");
		
		float sec1=R*(L-1)/2f;
		float sec2=(L-1)*(1-R/2f);
		
		float[] w=new float[L];
		
		for(int l=0;l<sec1;l++)
		w[l]=(float)((1.0+Math.cos(Math.PI*(2.0*l/(R*(L-1))-1)))/2.0);
		
		for(int l=Math.round(sec1),LL=Math.round(sec2);l<=LL;l++)
		w[l]=1f;
		
		for(int l=Math.round(sec2)+1;l<L;l++)
		w[l]=(float)((1.0+Math.cos(Math.PI*(2.0*l/(R*(L-1))-2.0/R+1)))/2.0);
		
		return w;
	}
	
	public static float[] normalizedHann(int L){ return normalizing(hann(L));}
	
	public static float[] normalizedHanning(int L){ return normalizing(hanning(L));}
	
	public static float[] normalizedHamming(int L){ return normalizing(hamming(L));}
	
	
	/*** helper method ***/
	private static float[] normalizing(float[] win){
		float sum=StatisticsUtil.cArithmeticMean(win);
		
		for(int i=0,I=win.length;i<I;i++) win[i]/=sum;
		
		return win;
	}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(Arrays.toString(Tukey(16,0.8f)));
	}*/
}
