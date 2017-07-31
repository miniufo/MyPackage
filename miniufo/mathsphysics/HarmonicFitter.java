/**
 * @(#)HarmonicFitter.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import miniufo.basic.ArrayUtil;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;


/**
 * Harmonic fitting class: y(t)=ym+a1*sin(o1*t)+b1*cos(o1*t)+a2*sin(o2*t)+b2*cos(o2*t)+...
 * Reference: My notebook
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class HarmonicFitter extends Fitter{
	//
	private float ym=0;			// fitted mean
	
	private  float[] Ts=null;	// periods of the given harmonics
	private  float[] as=null;	// coefficients of sin
	private  float[] bs=null;	// coefficients of cos
	private  float[] As=null;	// amplitudes
	private  float[] Ps=null;	// phases of sin
	private double[] Os=null;	// frequencies defined as Os = 2*PI/Ts
	
	private DMatrixRMaj A=null;
	private DMatrixRMaj X=null;
	
	
	/**
     * constructor
     *
     * @param	x	array of points where function are evaluated
     * @param	Ts	periods of the harmonics to be fit simultaneously (unit of x)
     */
	public HarmonicFitter(float[] x,float... Ts){
		super(x);
		
		this.Ts=Ts;
		
		as=new float[Ts.length];
		bs=new float[Ts.length];
		As=new float[Ts.length];
		Ps=new float[Ts.length];
		
		Os=new double[Ts.length];
		
		A=new DMatrixRMaj(len,Ts.length*2+1);
		X=new DMatrixRMaj(Ts.length*2+1,1);
		
		for(int i=0,I=Ts.length;i<I;i++) Os[i]=2.0*PI/Ts[i];
		
		for(int i=0,I=A.numRows;i<I;i++){
			A.set(i,0,1);
			
			for(int k=0,p=1,K=Ts.length;k<K;k++){
				A.set(i,p++,sin(Os[k]*x[i]));
				A.set(i,p++,cos(Os[k]*x[i]));
			}
		}
	}
	
	public HarmonicFitter(int len,float... Ts){ this(ArrayUtil.newMonotonousArray(len,1),Ts);}
	
	
	/**
     * fitting, calculate the coefficients of the harmonic by providing y=f(t)
     * Reference: My notebook
     *
     * @param	y	a series of observations at each x
     */
	public void fit(float[] y){
		if(y.length!=len) throw new IllegalArgumentException("lengths not equal");
		
		DMatrixRMaj Y=new DMatrixRMaj(len,1);
		
		for(int i=0;i<len;i++) Y.set(i,0,y[i]);
		
		boolean solved=CommonOps_DDRM.solve(A,Y,X);
		
		if(!solved) throw new IllegalArgumentException("cannot solve this");
		
		ym=(float)X.get(0,0);
		
		for(int i=0,p=1,I=Ts.length;i<I;i++){
			double a=X.get(p++,0);
			double b=X.get(p++,0);
			
			as[i]=(float)a;
			bs[i]=(float)b;
			As[i]=(float)Math.hypot(a,b);
			Ps[i]=(float)Math.atan2(b,a);
		}
	}
	
	
	/*** getor and setor ***/
	public float getMean(){ return ym;}
	
	public float[] getSinCoeffs(){ return as;}
	
	public float[] getCosCoeffs(){ return bs;}
	
	public float[] getAmplitudes(){ return As;}
	
	public float[] getSinInitPhases(){ return Ps;}
	
	public float[] getCosInitPhases(){
		float[] Pcos=Ps.clone();
		
		for(int i=0,I=Ps.length;i<I;i++) Pcos[i]=(float)Math.atan2(as[i],bs[i]);
		
		return Pcos;
	}
	
	
	/**
     * calculate values, stored in array r
     *
     * @param	x	a given array of locations to be evaluated at
     */
	public float[] cValues(float[] x){
		float[] r=new float[x.length];
		
		for(int l=0,L=x.length;l<L;l++){
			double tmp=0;
			
			for(int k=0,K=Ts.length;k<K;k++) tmp+=(as[k]*sin(Os[k]*x[l])+bs[k]*cos(Os[k]*x[l]));
			
			r[l]=(float)tmp;
		}
		
		return r;
	}
	
	public float[] cValues(){ return cValues(x);}
	
	
	/*** helper methods ***/
	
	
	/** test
	public static void main(String[] arg){
		int len=1201;
		
		int freq=1;
		
		double Tk1=23.934469*freq; double wk1=2*Math.PI/Tk1; double Ak1=1.2; double Pk1=56;
		double To1=25.81934 *freq; double wo1=2*Math.PI/To1; double Ao1=2.3; double Po1=126;
		double Tp1=24.06589 *freq; double wp1=2*Math.PI/Tp1; double Ap1=3.1; double Pp1=176;
		double Tq1=26.86836 *freq; double wq1=2*Math.PI/Tq1; double Aq1=1.9; double Pq1=206;
		
		double Tm2=12.420601*freq; double wm2=2*Math.PI/Tm2; double Am2=3.1; double Pm2=346;
		double Ts2=12.0     *freq; double ws2=2*Math.PI/Ts2; double As2=2.8; double Ps2=286;
		double Tn2=12.658348*freq; double wn2=2*Math.PI/Tn2; double An2=4.1; double Pn2=136;
		double Tk2=11.967234*freq; double wk2=2*Math.PI/Tk2; double Ak2=1.8; double Pk2=71;
		
		java.util.Random r=new java.util.Random();
		
		float[] y=new float[len];
		float[] x=new float[len];
		float[] N=new float[len];
		
		for(int l=0;l<len;l++){
			N[l]=(float)r.nextGaussian();
			x[l]=l*freq;
			y[l]=5.3f+(float)(
				Ak1*sin(wk1*x[l]+Pk1*Math.PI/180)+Ao1*sin(wo1*x[l]+Po1*Math.PI/180)+
				Am2*sin(wm2*x[l]+Pm2*Math.PI/180)+As2*sin(ws2*x[l]+Ps2*Math.PI/180)+N[l]
			);
		}
		
		HarmonicFitter2 hf1=new HarmonicFitter2(len,(float)Tk1,(float)To1,(float)Tm2,(float)Ts2); hf1.fit(y);
		
		float[] As=hf1.getAmplitudes();
		float[] Ps=hf1.getSinInitPhases();
		
		System.out.println(hf1.getMean());
		
		for(int i=0;i<Ps.length;i++){
			Ps[i]*=180/Math.PI;
			System.out.println(As[i]+"\t"+Ps[i]);
		}
		
		float[] yc=hf1.cValues();
		System.out.println();
		for(int l=0;l<len;l++) System.out.println(y[l]+"\t"+(yc[l]+hf1.getMean()));
	}*/
}
