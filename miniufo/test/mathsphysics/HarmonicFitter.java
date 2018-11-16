/**
 * @(#)HarmonicFitter.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.mathsphysics;

import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import miniufo.basic.ArrayUtil;
import miniufo.mathsphysics.Fitter;


/**
 * Harmonic fitting class: y(t)=ym+[a*sin(omega*t)+b*cos(omega*t)]+n(t)
 * Reference: My notebook
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class HarmonicFitter extends Fitter{
	//
	private float a=0;		// coefficient of sin function
	private float b=0;		// coefficient of cos function
	private float amp=0;	// amplitude of harmonic, amp=sqrt(a*a+b*b)
	private float sph=0;	// initial phase for sin of harmonic
	private float omega=0;	// pulsation, omega=2*PI/T
	
	private double sumYSin=0;	// sum of   y*sin over all t
	private double sumYCos=0;	// sum of   y*cos over all t
	private double sumSin2=0;	// sum of sin*sin over all t
	private double sumCos2=0;	// sum of cos*cos over all t
	private double sumSC  =0;	// sum of sin*cos over all t
	
	private double[] sinV=null;	// buffer for sin(omega*l)
	private double[] cosV=null;	// buffer for cos(omega*l)
	
	
	/**
     * constructor
     *
     * @param	T	period of the harmonic to be fit (unit of x)
     * @param	x	array of points where function are evaluated
     */
	public HarmonicFitter(float T,float[] x){
		super(x);
		
		omega=(float)(2.0*PI/T);
		
		this.x=x;
		
		initTRelatedBuffer(x);
	}
	
	public HarmonicFitter(float T,int len){ this(T,ArrayUtil.newMonotonicArray(len,1));}
	
	
	/**
     * fitting, calculate the coefficients of the harmonic by providing y=f(t)
     * Reference: My notebook
     *
     * @param	y	a series of observations at each x
     */
	public void fit(float[] y){
		if(y.length!=len) throw new IllegalArgumentException("lengths not equal");
		
		initYRelatedBuffer(y);
		
		b=(float)((sumYCos*sumSin2-sumYSin*sumSC)/(sumSin2*sumCos2-sumSC*sumSC));
		a=(float)((sumYSin-b*sumSC)/sumSin2);
		
		amp=(float)Math.hypot(a,b);
		sph=(float)Math.atan2(b,a);
	}
	
	
	/*** getor and setor ***/
	public float getSinCoeff(){ return a;}
	
	public float getCosCoeff(){ return b;}
	
	public float getAmplitude(){ return amp;}
	
	public float getSinInitPhase(){ return sph;}
	
	public float getCosInitPhase(){ return sph-(float)(PI/2.0);}
	
	
	/**
     * calculate values, stored in array r
     *
     * @param	x	a given array of locations to be evaluated at
     */
	public float[] cValues(float[] x){
		float[] r=new float[x.length];
		
		for(int l=0,L=x.length;l<L;l++) r[l]=(float)(a*sin(omega*x[l])+b*cos(omega*x[l]));
		
		return r;
	}
	
	public float[] cValues(){
		float[] r=new float[len];
		
		for(int l=0;l<len;l++) r[l]=(float)(a*sinV[l]+b*cosV[l]);
		
		return r;
	}
	
	
	/*** helper methods ***/
	private void initTRelatedBuffer(float[] t){
		int len=t.length;
		
		sinV=new double[len];
		cosV=new double[len];
		
		for(int l=0;l<len;l++){
			sinV[l]=sin(omega*t[l]);
			cosV[l]=cos(omega*t[l]);
		}
		
		// clear all buffer
		sumSin2=0;sumCos2=0;sumSC=0;
		
		for(int l=0;l<len;l++){
			sumSin2+=sinV[l]*sinV[l];
			sumCos2+=cosV[l]*cosV[l];
			sumSC  +=sinV[l]*cosV[l];
		}
	}
	
	private void initYRelatedBuffer(float[] y){
		int len=y.length;
		// clear all buffer
		sumYSin=0;sumYCos=0;
		
		for(int l=0;l<len;l++){
			sumYSin+=y[l]*sinV[l];
			sumYCos+=y[l]*cosV[l];
		}
	}
	
	
	/** test
	public static void main(String[] arg){
		int len=400;
		
		float ym=30;
		float a1=-11;
		float b1=-13;
		float T1=63;
		float a2=4;
		float b2=5;
		float T2=21;
		
		java.util.Random r=new java.util.Random();
		
		float[] y=new float[len];
		float[] x=new float[len];
		float[] N=new float[len];
		
		for(int l=0;l<len;l++){
			N[l]=(float)r.nextGaussian();
			x[l]=l;
			y[l]=//ym+
				a1*(float)sin(2*PI/T1*l)+b1*(float)cos(2*PI/T1*l)+
				a2*(float)sin(2*PI/T2*l)+b2*(float)cos(2*PI/T2*l);//+
				//N[l];
		}
		
		HarmonicFitter hf1=new HarmonicFitter(T1,len); hf1.fit(y);
		HarmonicFitter hf2=new HarmonicFitter(T2,len); hf2.fit(y);
		
		System.out.println("a1:"+a1+"("+hf1.getSinCoeff()+")\tb1:"+b1+"("+hf1.getCosCoeff()+")\t");
		System.out.println("a2:"+a2+"("+hf2.getSinCoeff()+")\tb2:"+b2+"("+hf2.getCosCoeff()+")\t");
		
		for(int l=0;l<len;l++)
		System.out.println(
			y[l]+"\t"+
			(hf1.getAmplitude()*sin(2*PI/T1*l+hf1.getSinInitPhase()))+"\t"+
			(hf1.getAmplitude()*cos(2*PI/T1*l+hf1.getCosInitPhase()))+"\t"+
			(hf1.getSinCoeff()*sin(2*PI/T1*l)+hf1.getCosCoeff()*cos(2*PI/T1*l))+"\t"+
			(hf2.getAmplitude()*sin(2*PI/T2*l+hf2.getSinInitPhase()))+"\t"+
			(hf2.getAmplitude()*cos(2*PI/T2*l+hf2.getCosInitPhase()))+"\t"+
			(hf2.getSinCoeff()*sin(2*PI/T2*l)+hf2.getCosCoeff()*cos(2*PI/T2*l))
		);
	}*/
	/*
	public static void main(String[] arg){
		int len=1201;
		
		double Tk1=23.934469*3600; double wk1=2*Math.PI/Tk1; double Ak1=1.2; double Pk1=56;
		double To1=25.81934 *3600; double wo1=2*Math.PI/To1; double Ao1=2.3; double Po1=126;
		double Tp1=24.06589 *3600; double wp1=2*Math.PI/Tp1; double Ap1=3.1; double Pp1=176;
		double Tq1=26.86836 *3600; double wq1=2*Math.PI/Tq1; double Aq1=1.9; double Pq1=206;
		
		double Tm2=12.420601*3600; double wm2=2*Math.PI/Tm2; double Am2=3.1; double Pm2=346;
		double Ts2=12.0     *3600; double ws2=2*Math.PI/Ts2; double As2=2.8; double Ps2=286;
		double Tn2=12.658348*3600; double wn2=2*Math.PI/Tn2; double An2=4.1; double Pn2=136;
		double Tk2=11.967234*3600; double wk2=2*Math.PI/Tk2; double Ak2=1.8; double Pk2=71;
		
		java.util.Random r=new java.util.Random();
		
		float[] y=new float[len];
		float[] x=new float[len];
		float[] N=new float[len];
		
		for(int l=0;l<len;l++){
			N[l]=(float)r.nextGaussian();
			x[l]=l*3600;
			y[l]=5.3f+(float)(
				Ak1*sin(wk1*x[l]+Pk1*Math.PI/180)//+Ao1*sin(wo1*x[l]+Po1*Math.PI/180)+
				//Am2*sin(wm2*x[l]+Pm2*Math.PI/180)+As2*sin(ws2*x[l]+Ps2*Math.PI/180)
			);
		}
		
		HarmonicFitter hf1=new HarmonicFitter((float)Tk1,x); hf1.fit(y);
		
		System.out.println(hf1.getAmplitude()+"\t"+hf1.getCosInitPhase());
		
		
		for(int l=0;l<len;l++)
		System.out.println(
			y[l]+"\t"+
			(hf1.getAmplitude()*sin(2*PI/Tk1*l+hf1.getSinInitPhase()))+"\t"+
			(hf1.getAmplitude()*cos(2*PI/Tk1*l+hf1.getCosInitPhase()))+"\t"+
			(hf1.getSinCoeff()*sin(2*PI/Tk1*l)+hf1.getCosCoeff()*cos(2*PI/Tk1*l))
		);
	}*/
}
