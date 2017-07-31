/**
 * @(#)MathsPhysicsUtil.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import java.util.Arrays;


/**
 * common methods concerning math-physics
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class MathsPhysicsUtil{
	
	/**
	 * prevent from instantiate
     */
	private MathsPhysicsUtil(){}
	
	
	/**
	 * Returns true if x is a power of 2
	 */
	public static boolean isPowerOf2(int x){
		// BITS_PER_WORD==32
		for(int i=1,y=2;i<32;i++,y<<=1)
		if(x==y) return true;
		
		return false;
	}
	
	
	/**
	 * Number of factor 2, n of 2s multiply like 2*2*...*2==pwrOf2
	 *
	 * @param	pwrOf2	a given number which contains n of 2s as factor
	 *
	 * @return	i		number of 2s
	 */
	public static int numberOf2(int pwrOf2){
		if(!isPowerOf2(pwrOf2))
		throw new IllegalArgumentException("the arguments is not power-of-2 number");
		
		for(int i=0;;i++)
		if((pwrOf2&(1<<i))>0) return i;
	}
	
	
	/**
	 * factoring a given number val.
	 *
	 * @param	val		a given number to factoring
	 *
	 * @return	a vector containing the prime factors of N
	 */
	public static int[] factor(int val){
		if(val<1) throw new IllegalArgumentException("val ("+val+") should be larger than 0");
		
		if(val==1) return new int[]{1};
		
		int limit=(int)Math.ceil(Math.sqrt(val));
		int index=0;
		
		int[] numArray=new int[limit];
		
		for(int i=2;i<=limit;i++)
		while(val%i==0){
			numArray[index++]=i;
			val/=i;
		}
		
		if(val!=1) numArray[index++]=val;
		
		return Arrays.copyOf(numArray,index);
	}
	
	
	/**
     * factorial function, = n!
     *
     * @param	n	an given integer
     *
     * @return	n!
     */
	public static int factorialInt(int n){
		if(n< 0) throw new IllegalArgumentException("n should be larger than 0");
		if(n>12) throw new IllegalArgumentException("result over flow");
		
		if(n==0||n==1) return 1;
		
		int re=1;
		
		for(int i=2;i<=n;i++) re*=i;
		
		return re;
	}
	
    public static int twoFactorialInt(int n){
		if(n< 0) throw new IllegalArgumentException("n should be larger than 0");
		if(n>19) throw new IllegalArgumentException("result over flow");
		
        int re=1;
        
        for(int i=n;i>=2;i-=2) re*=i;
        
        return re;
    }
    
    public static double factorialDouble(int n){
		if(n<0) throw new IllegalArgumentException("n should be larger than 0");
		
		if(n==0||n==1) return 1.0;
		
		double re=1;
		
		for(int i=2;i<=n;i++) re*=i;
		
		if(re<=0) throw new IllegalArgumentException("result over flow");
		
		return re;
	}
	
    public static double twoFactorialDouble(int n){
		if(n<0) throw new IllegalArgumentException("n should be larger than 0");
		
        double re=1;
        
        for(int i=n;i>=2;i-=2) re*=i;
		
		if(re<=0) throw new IllegalArgumentException("result over flow");
        
        return re;
    }
	
	
	/**
	 * gamma(x) = integral from 0 to infinite of t^(x-1) exp(-t) dt
     * The gamma function interpolates the factorial function.
     * For integer n, gamma(n+1) = n! (n factorial)
     *
     * @param	n	an given integer
     *
     * @return	gamma value
     */
	public static double gamma(int x){ return factorialDouble(x-1);}
	
	public static double gamma(double x){
		return Math.exp(org.apache.commons.math3.special.Gamma.logGamma(x));
	}
	
	
	/**
	 * convolution of data using the window function h (mask)
	 * equals to matlab code: y=conv(data,h,'full')
     *
     * @param	h		window function (normalized)
     * @param	data	a series of data
     *
     * @return	re		returns the full convolution
	 */
	public static float[] convolve(float[] data,float[] h){
		int hLen=h.length;
		int dLen=data.length+hLen-1;
		
		float[] re=new float[dLen];	// result
		float[] dp=new float[dLen];	// data padded with 0
		float[] hp=new float[dLen];	// window padded with 0
		
		System.arraycopy(data,0,dp,0,data.length);
		System.arraycopy(h   ,0,hp,0,hLen);
		
		for(int i=0;i<dLen;i++)
		for(int j=0,J=data.length;j<J;j++) if(i-j>=0) re[i]+=dp[j]*hp[i-j];
		
		return re;
	}
	
	/**
	 * convolution of data using the window function h (mask)
	 * equals to matlab code: y=conv(data,h,'same')
     *
     * @param	h		window function
     * @param	data	a series of data
     *
     * @return	re		the central part of the convolution that is the same size as A
	 */
	public static float[] convolveSame(float[] data,float[] h){
		float[] re=new float[data.length];
		
		System.arraycopy(convolve(data,h),h.length/2,re,0,data.length);
		
		return re;
	}
	
	
	/**
	 * compute area of a triangle using Heron formula given the lengths of three sides
     * area = sqrt( s * (s-s1) * (s-s2) * (s-s3) ) where s = (s1+s2+s3)/2
     * 
     * @param	s1	length of the first side
     * @param	s2	length of the second side
     * @param	s3	length of the third side
     *
     * @return	re	area
     */
	public static float cTriangleAreaUsingHeronFormula(float s1,float s2,float s3){
		float s=(s1+s2+s3)/2f;
		
		return (float)Math.sqrt(s*(s-s1)*(s-s2)*(s-s3));
	}
	
	public static double cTriangleAreaUsingHeronFormula(double s1,double s2,double s3){
		double s=(s1+s2+s3)/2.0;
		
		return Math.sqrt(s*(s-s1)*(s-s2)*(s-s3));
	}
	
	
	/**
	 * compute the normalized area (0-1) on a sphere bounded by two latitudes and longitudes
	 * following Matlab command areaquad.
	 * 
	 * @param	lon1	longitude of the southwest corner (Radian)
	 * @param	lat1	latitudes of the southwest corner (Radian)
	 * @param	lon2	longitude of the northeast corner (Radian)
	 * @param	lat2	latitudes of the northeast corner (Radian)
	 * 
	 * @return	area bounded by two latitudes and two longitudes
	 */
	public static float cAreaQuadByRadian(float lon1,float lat1,float lon2,float lat2){
		if(lon1==lon2||lat1==lat2) return 0;
		return (float)(Math.abs(lon1-lon2)*Math.abs(Math.sin(lat1)-Math.sin(lat2))/(4.0*Math.PI));
	}
	
	public static double cAreaQuadByRadian(double lon1,double lat1,double lon2,double lat2){
		if(lon1==lon2||lat1==lat2) return 0;
		return Math.abs(lon1-lon2)*Math.abs(Math.sin(lat1)-Math.sin(lat2))/(4.0*Math.PI);
	}
	
	
	/**
	 * compute Haversine function : hav(x) = (1-cos(x))/2 = sin^2(x/2)
	 */
	public static double Haversine(double x){ return (1.0-Math.cos(x))/2.0;}
	
	public static float Haversine(float x){ return (float)Haversine((double)x);}
	
	
	/*** test **
	public static void main(String[] args){
		for(int i=1;i<500;i++)
		System.out.println(Arrays.toString(factor(i)));
	}*/
}
