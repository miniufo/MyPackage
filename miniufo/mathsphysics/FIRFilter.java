/**
 * @(#)FIRFilter.java	1.0 2013.03.12
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.statistics.StatisticsUtil;


/**
 * finite impulse filter
 *
 * @version 1.0, 2013.03.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class FIRFilter{
	
	/**
	 * prevent from instantiate
     */
	private FIRFilter(){}
	
	
	/**
	 * digital filter function
	 * 
	 * @param	fc		cutoff frequency
	 * @param	win		window function
	 * @param	data	a time series to be filtered
	 * 
	 * @return	re		filtered result
	 */
	public static float[] windowFiltering(float fc,float[] win,float[] data){ return filter(fir1(fc,win),data);}
	
	
	/**
	 * Y = filter(B,X) filters the data in vector X with the
	 * filter described by vectors A and B to create the filtered
	 * data Y.  The filter is a "Direct Form II Transposed"
	 * implementation of the standard difference equation:
	 * y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
	 * same as matlab function y = filter(B,1,X);
	 * 
	 * @param	B		coefficient of filter
	 * @param	data	series to be filtered
	 * 
	 * @return	re		filtered result
	 */
	public static float[] filter(float[] B,float[] data){
		float[] re=new float[data.length];
		
		for(int i=0,I=data.length;i<I;i++){
			float sum=0;
			for(int j=0,J=B.length;j<J;j++)
			if(i-j>=0) sum+=B[j]*data[i-j];
			re[i]=sum;
		}
		
		return re;
	}
	
	/**
	 * B = fir1(Wn,wind(N+1)) designs an N'th order lowpass FIR digital
	 * filter and returns the filter coefficients in length N+1 vector B.
	 * similar to matlab command B = fir1(N,Wn,wind(N+1));
	 * 
	 * @param	fc	cutoff frequency between 0 < Wn < 1.0,
	 * 				with 1.0 corresponding to half the sample rate.
	 * @return	B	filter coefficients that are real and has linear phase.
	 * 				The normalized gain of the filter at Wn is -6 dB.
	 */
	public static float[] fir1(float fc,float[] win){
		int N=win.length;
		int odd=N%2;
		int nhlf=(int)Math.floor((N+1)/2f);
		
		if(fc>=1||fc<=0) throw new IllegalArgumentException("fc should be 0 <= fc <= 1");
		
		float[] B =new float[N];
		float[] xn=new float[nhlf];
		
		for(int i=0,ii=0;i<nhlf;i++,ii++) xn[i]=ii+0.5f*(1-odd);
		
		float[] tmp=sinc(fc,xn);
		
		for(int i=0,I=(int)Math.ceil(B.length/2f);i<I;i++){
			B[i]=tmp[I-i-1];
			B[B.length-1-i]=B[i];
		}
		
		for(int i=0;i<N;i++) B[i]*=win[i];
		
		float sum=Math.abs(StatisticsUtil.sum(B));
		
		for(int i=0;i<N;i++) B[i]/=sum;
		
		return B;
	}
	
	
	/*** helper method ***/
	/**
	 * Sin(c*pi*x)/(pi*x) function returns a matrix whose elements are
	 * the sinc of the elements of X, i.e.
	 *    y = sin(c1*pi*x)/(pi*x)    if x ~= 0
	 *      = 1                      if x == 0
	 * 
	 * @param	c	an coefficient
	 * @param	x	an array
	 * 
	 * @return	re	y = Sin(c1*pi*x)/(pi*x)
	 */
	private static float[] sinc(float c,float[] x){
		int len=x.length;
		
		float[] re=new float[len];
		
		for(int i=0;i<len;i++){
			double de=Math.PI*x[i];
			
			if(x[i]!=0) re[i]=(float)(Math.sin(c*de)/de);
			else re[i]=c;
		}
		
		return re;
	}
	
	
	/** test
	public static void main(String[] args){
		float[] B=new float[]{1,2,3};
		float[] data=new float[6];
		for(int i=0;i<data.length;i++) data[i]=(float)Math.sin(2*Math.PI*(1+0.2*i));
		float[] y=filter(B,data);
		
		System.out.println("B   : "+Arrays.toString(B));
		System.out.println("data: "+Arrays.toString(data));
		System.out.println("y   : "+Arrays.toString(y));
	}*/
}
