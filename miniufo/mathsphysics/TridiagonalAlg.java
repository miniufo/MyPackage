/**
 * @(#)TridiagonalAlg.java	1.0 2015.01.13
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;


/**
 * Using tridiagonal algorithm to solve the tridiagonal matrix of the form Ax=d
 * where A is tridiagonal matrix.
 * 
 * Reference:
 * http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 *
 * @version 1.0, 2015.01.13
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TridiagonalAlg{
	//
	private int n;				// length
	
	private double[] buf0=null;	// buffer with length n
	private double[] buf1=null;	// buffer with length n
	private double[] buf2=null;	// buffer with length n
	private double[] buf3=null;	// buffer with length n
	private double[] buf4=null;	// buffer with length n
	private double[] buf5=null;	// buffer with length n-1
	
	
	/**
     * constructor
     *
     * @param	length	for a length*length tridiagonal matrix
     */
	public TridiagonalAlg(int length){
		n=length;
		
		buf0=new double[n];	buf3=new double[n];
		buf1=new double[n];	buf4=new double[n];
		buf2=new double[n];	buf5=new double[n-1];
	}
	
	
	/*** getor and setor ***/
	public int getLength(){ return n;}
	
	
	/**
	 * trace method, solving equation of tri-diagonal matrix.
	 *
	 * @param	a	lower coefficients of the matrix
	 * @param	b	diagonal coefficients of the matrix
	 * @param	c	upper coefficients of the matrix
	 * @param	d	vector on the right of the equation
	 * @param	U	solution of the unknowns (output)
	 */
	public void trace(double[] a,double[] b,double[] c,double[] d,double[] U){
		int alen=a.length;
		
		if(a.length!=n-1) throw new IllegalArgumentException("length of given a ("+a.length+") should be "+(n-1));
		if(b.length!=n  ) throw new IllegalArgumentException("length of given b ("+b.length+") should be "+(n  ));
		if(c.length!=n-1) throw new IllegalArgumentException("length of given c ("+c.length+") should be "+(n-1));
		if(d.length!=n  ) throw new IllegalArgumentException("length of given d ("+d.length+") should be "+(n  ));
		
		buf5[0]=c[0]/b[0];
		buf0[0]=b[0];
		
		for(int i=1;i<alen;i++){
			buf0[i]=b[i]-a[i-1]*buf5[i-1];
			buf5[i]=c[i]/buf0[i];
		}
		
		buf0[alen]=b[alen]-a[n-2]*buf5[n-2];
		
		U[0]=d[0]/buf0[0];
		
		for(int i=1;i<n;i++) U[i]=(d[i]-a[i-1]*U[i-1])/buf0[i];
		
		// calculate U[i]
		for(int i=n-2;i>=0;i--) U[i]-=buf5[i]*U[i+1];
	}
	
	public double[] trace(double[] a,double[] b,double[] c,double[] d){
		double[] U=new double[n];
		
		trace(a,b,c,d,U);
		
		return U;
	}
	
	
	/**
	 * trace method, solving equation of cyclic tri-diagonal matrix
	 *
	 * @param	a	lower  coefficient of the matrix
	 * @param	b	diagonal coefficient of the matrix
	 * @param	c	upper coefficient of the matrix
	 * @param	d	vector on the right of the equation
	 * @param	a0	cyclic coefficient
	 * @param	cn	cyclic coefficient
	 * @param	U	solution of the unknowns (output)
	 */
	public void traceCyclic(double[] a,double[] b,double[] c,double[] d,double a0,double cn,double[] U){
		// calculate buf1
		buf4[n-1]=cn;
		trace(a,b,c,buf4,buf1);
		
		// calculate buf2
		buf4[n-1]=0;	buf4[0]=a0;
		trace(a,b,c,buf4,buf2);
		
		// calculate buf3
		buf4[n-1]=0;	buf4[0]=a0;
		trace(a,b,c,d,buf3);
		
		U[n-1]=((1+buf1[0])/buf1[n-1]*buf3[n-1]-buf3[0])/((1+buf1[0])*(1+buf2[n-1])/buf1[n-1]-buf2[0]);
		U[ 0 ]=(buf3[0]-buf2[0]*U[n-1])/(1+buf1[0]);
		
		for(int i=1,I=n-1;i<I;i++) U[i]=buf3[i]-buf1[i]*U[0]-buf2[i]*U[n-1];
	}
	
	public double[] traceCyclic(double[] a,double[] b,double[] c,double[] d,double a0,double cn){
		double[] U=new double[n];
		
		traceCyclic(a,b,c,d,a0,cn,U);
		
		return U;
	}
	
	
	/** test
	public static void main(String[] args){
		double[] a={2,2,0};
		double[] b={3,3,3,3};
		double[] c={0,1,1};
		double[] d={5,9,9,8};
		
		double[] U=new TridiagonalAlg(4).trace(a,b,c,d);
		
		for(int i=0;i<4;i++) System.out.println(U[i]);
	}*/
}
