/**
 * @(#)PolynomialFitter.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import org.ejml.simple.SimpleMatrix;
import miniufo.basic.ArrayUtil;


/**
 * Polynomial fitting class: 1+x^1+x^2+...+x^n
 * Reference: The Basis of Numerical Computation (by Jianhua Sheng)
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class PolynomialFitter extends Fitter{
	//
	private int npoly;
	
	private double[] coeff=null;	// the coefficients of a polynomial
	
	private double[][] ary=null;
	private double[][] tmp=null;
	
	
	/**
     * constructor
     *
     * @param	l	order of polynomial fitting
     * @param	x	array of points where function are evaluated
     */
	public PolynomialFitter(int l,float[] x){
		super(x);
		
		if(l<1)
		throw new IllegalArgumentException("order should be larger than or equal 1");
		
		if(l>3)
		System.out.println("Warning: the result may not be good since order is higher than 3");
		
		npoly=l+1;
		
		coeff=new double[npoly];
		ary  =new double[npoly][npoly];
		tmp  =new double[npoly][1];
	}
	
	public PolynomialFitter(int l,int len){ this(l,ArrayUtil.newMonotonousArray(len,1));}
	
	
	/*** getor and setor ***/
	public float getCoefficient(int order){ return (float)coeff[order];}
	
	
	/**
     * fitting, calculate the coefficients of the polynomial
     * by providing y=f(x)
     * Reference: The Basis of Numerical Computation (by Jianhua Sheng)
     *            Equation 2.77 at page 76.
     *
     * @param	y	a series of observations
     * @param	x	where observations are located at
     */
	public void fit(float[] y){
		clear(ary);	clear(tmp);
		
		if(len!=y.length) throw new IllegalArgumentException("lengths not equal");
		
		for(int j=0;j<npoly;j++){
			for(int i=j;i<npoly;i++)
			for(int k=0;k<len;k++) ary[j][i]+=Math.pow(x[k],i+j);
			
			for(int i=0;i<j;i++) ary[j][i]=ary[i][j];
		}
		
		for(int i=0;i<npoly;i++)
		for(int k=0;k<len;k++) tmp[i][0]+=Math.pow(x[k],i)*y[k];
		
		SimpleMatrix ml=new SimpleMatrix(ary);
		SimpleMatrix mr=new SimpleMatrix(tmp);
		SimpleMatrix re=ml.solve(mr);
		
		for(int i=0;i<npoly;i++) coeff[i]=re.get(i,0);
	}
	
	
	/**
     * calculate values, stored in array r
     *
     * @param	x	a given array of locations to be evaluated at
     */
	public float[] cValues(float[] x){
		if(coeff==null)	throw new NullPointerException("Please use fit(int n) first");
		
		float[] r=new float[x.length];
		
		for(int i=0,I=x.length;i<I;i++){
			r[i]=0;
			
			for(int j=0;j<npoly;j++) r[i]+=coeff[j]*Math.pow(x[i],j);
		}
		
		return r;
	}
	
	public float[] cValues(){ return cValues(x);}
	
	
	/*** helper methods ***/
	private void clear(double[][] data){
		for(int j=0,J=data.length;j<J;j++)
		for(int i=0,I=data[0].length;i<I;i++) data[j][i]=0;
	}
	
	
	/** test
	public static void main(String[] arg){
		float[] y={1, 1.284f, 1.6487f, 2.117f, 2.7183f};
		float[] x={1, 2f, 3f, 4f, 5f};
		
		PolynomialFitter pf=new PolynomialFitter(2,x);
		
		pf.fit(y);
		
		float[] re=pf.cValues(x);
		
		for(int i=0;i<re.length;i++) System.out.println(re[i]);
	}*/
}
