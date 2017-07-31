/**
 * @(#)Spline.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.matrix;

import miniufo.basic.ArrayUtil;
import org.ejml.simple.SimpleMatrix;


/**
 * spline class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Spline{
	//
	private boolean incre=true;
	
	private float[] sx=null;	// given x vector
	private float[] sy=null;	// given y value vector refer to the x vector
	private float[] h =null;	// step length
	
	private float[] a =null;	// down-band coefficient in the matrix
	private float[] c =null;	// up-band coefficient int the matrix
	private float[] d =null;	// vector on the right of the equation
	
	private float[] m =null;	// derivation vector
	
	private float[] o3=null;	// coefficient vector of 3-order term
	private float[] o2=null;	// coefficient vector of 2-order term
	private float[] o1=null;	// coefficient vector of 1-order term
	
	private int slength=0;
	
	
	/**
     * constructor
     *
     * @param	sx	given x vector
     * @param	sy	given y value vector refer to the x vector
     */
	public Spline(float[] sx,float[] sy){
		slength=sx.length;
		
		if(slength!=sy.length) throw new IllegalArgumentException("lengths mismatch");
		if(slength<2) throw new IllegalArgumentException("array's length should be larger than 1");
		
		this.sx=sx;	this.sy=sy;
		
		incre=sx[slength-1]>sx[0];
		
		if(!incre){
			this.sx=sx.clone(); ArrayUtil.reverse(this.sx);
			this.sy=sy.clone(); ArrayUtil.reverse(this.sy);
		}
		
		// if no periodic trace, a c d lengths should be slength-3, slength-3, slength-2
		if(slength>2){
			a=new float[slength-2];
			c=new float[slength-2];
			d=new float[slength-1];
		}
		
		o1=new float[slength-1];
		o2=new float[slength-1];	m=new float[slength];
		o3=new float[slength-1];	h=new float[slength-1];
		
		/*** calculate heap ***/
		for(int i=0,I=slength-1;i<I;i++) h[i]=this.sx[i+1]-this.sx[i];
	}
	
	
	public void updateData(float[] sy){
		if(sy.length!=slength)
		throw new IllegalArgumentException("invalid length of array ("+sy.length+"), should be "+slength);
		
		if(incre) System.arraycopy(sy,0,this.sy,0,slength);
		else for(int i=0;i<slength;i++) this.sy[i]=sy[slength-1-i];
	}
	
	
	/**
     * Cubic spline interpolation with first-order derivatives
     * (complete or clamped) boundary condition.
     *
     * @param	c1	left boundary value
     * @param	c2	right boundary value
     */
	public void cubicSplineWith1stBC(float c1,float c2){
		if(slength==2){
			o1[0]=c1;
			
			double s2=h[0]*h[0],s3=s2*h[0];
			
			double[][] coef=new double[][]{
				{s3  ,s2    },
				{3*s2,2*h[0]}
			};
			
			double[][] y=new double[][]{
				{sy[1]-sy[0]-c1*h[0]},
				{c2-c1}
			};
			
			SimpleMatrix A=new SimpleMatrix(coef);
			SimpleMatrix b=new SimpleMatrix(y);
			
			b=A.solve(b);
			
			o3[0]=(float)b.get(0);
			o2[0]=(float)b.get(1);
			
			return;
		}
		
		/*** calculate a,c,d ***/
		for(int i=0,I=slength-3;i<I;i++) a[i]=h[i+2]/(h[i+1]+h[i+2]);
		
		for(int i=0,I=slength-3;i<I;i++) c[i]=h[i]/(h[i+1]+h[i]);
		
		for(int i=0,I=slength-2;i<I;i++) d[i]=3/(h[i]+h[i+1])*((sy[i+2]-sy[i+1])*h[i]/h[i+1]+(sy[i+1]-sy[i])*h[i+1]/h[i]);
		d[0]-=h[1]*c1/(h[0]+h[1]);	d[slength-3]-=h[slength-3]*c2/(h[slength-2]+h[slength-3]);
		
		/*** using trace method to get x vector ***/
		trace(a,c,d);	m[0]=c1;	m[slength-1]=c2;
		
		/*** calculate o3,o2,o1 ***/
		for(int i=0,I=slength-1;i<I;i++){
			o3[i]=(2*(sy[i]-sy[i+1])/h[i]+m[i]+m[i+1])/h[i]/h[i];
			o2[i]=(3*(sy[i+1]-sy[i])/h[i]-2*m[i]-m[i+1])/h[i];
			o1[i]=m[i];
		}
	}
	
	/**
     * cubic spline interpolation with second-order derivatives boundary condition
     *
     * @param	c1	left boundary value
     * @param	c2	right boundary value
     */
	public void cubicSplineWith2ndBC(float c1,float c2){
		if(slength==2){
			o2[0]=c1/2f;
			o3[0]=(c2-2*o2[0])/(6*h[0]);
			o1[0]=(sy[1]-sy[0])/h[0]-o3[0]*h[0]*h[0]-o2[0]*h[0];
			
			return;
		}
		
		/*** calculate a,c,d ***/
		for(int i=0,I=slength-3;i<I;i++) a[i]=h[i+1]/(h[i+1]+h[i+2]);
		
		for(int i=0,I=slength-3;i<I;i++) c[i]=h[i+1]/(h[i+1]+h[i]);
		
		for(int i=0,I=slength-2;i<I;i++) d[i]=6/(h[i]+h[i+1])*((sy[i+2]-sy[i+1])/h[i+1]-(sy[i+1]-sy[i])/h[i]);
		d[0]-=h[0]*c1/(h[0]+h[1]);	d[slength-3]-=h[slength-2]*c2/(h[slength-2]+h[slength-3]);
		
		/*** using trace method to get x vector ***/
		trace(a,c,d);	m[0]=c1;	m[slength-1]=c2;
		
		/*** calculate o3,o2,o1 ***/
		for(int i=0,I=slength-1;i<I;i++){
			o3[i]=(m[i+1]-m[i])/h[i]/6;
			o2[i]=m[i]/2;
			o1[i]=(sy[i+1]-sy[i])/h[i]-(m[i+1]+2*m[i])*h[i]/6;
		}
	}
	
	/**
     * cubic spline interpolation with natural boundary condition
     */
	public void cubicSplineWithNaturalBC(){ cubicSplineWith2ndBC(0,0);}
	
	/**
     * cubic spline interpolation with periodic boundary condition
     * 
     * @param	step	step length between x[0] and x[x.length-1]
     */
	public void cubicSplineWithPeriodicBC(float step){
		float c1=(sy[1]-sy[slength-1])/(step+h[0]);	// 1st derivative at start point
		float c2=(sy[0]-sy[slength-2])/(h[slength-2]+step);	// 1st derivative at end point
		
		cubicSplineWith1stBC(c1,c2);
	}
	
	
	/**
     * calculate values according to the given x
     *
     * @param	x	given x vector
     *
     * @return	y	result values vector
     */
	public float[] cValues(float[] x){
		float[] y=new float[x.length];
		
		cValues(x,y);
		
		return y;
	}
	
	/**
     * calculate values according to the given x
     *
     * @param	x	given x vector
     * @param	y	result values vector
     */
	public void cValues(float[] x,float[] y){
		int length=x.length;
		
		if(length!=y.length)
		throw new IllegalArgumentException("array lengths mismatch");
		
		for(int i=0;i<length;i++){
			int idx=0;
			
			idx=ArrayUtil.getLEIdxIncre(sx,x[i]);
			
			if(idx==slength-1){ y[i]=sy[slength-1]; continue;} // index out-of-boundary case
			
			float dx=x[i]-sx[idx];
			
			y[i]=o3[idx]*dx*dx*dx+o2[idx]*dx*dx+o1[idx]*dx+sy[idx];
		}
	}
	
	
	/**
     * use trace method to solve the tridiagonal matrix equation
     * this method has side-effect, values of array d would be changed in order to save memory
     *
     * @param	a	down-band coefficient in the matrix
     * @param	c	up-band coefficient int the matrix
     * @param	d	vector on the right of the equation
     */
	private void trace(float[] a,float[] c,float[] d){
		/*** calculate l ***/
		m[1]=2;	// l[0]=2
		
		for(int i=2;i<slength-1;i++) m[i]=2-a[i-2]*c[i-2]/m[i-1];	// x[i]=2-a[i-1]*c[i-1]/l[i-1];
		
		/*** calculate y ***/
		d[0]/=m[1];	//y[0]=d[0]/l[0];
		
		for(int i=1;i<slength-2;i++) d[i]=(d[i]-a[i-1]*d[i-1])/m[i+1];	// y[i]=(d[i]-a[i-1]*y[i-1])/l[i];
		
		/*** calculate x ***/
		m[slength-2]=d[slength-3];	// x[slength-3]=y[slength-3];
		
		for(int i=slength-3;i>0;i--) m[i]=d[i-1]-c[i-1]*m[i+1]/m[i];	// x[i]=y[i]-c[i]*x[i+1]/l[i];
	}
	
	
	/** test
	public static void main(String[] arg){
		float[] sx=new float[]{1,3,4,6,9,20};
		//float[] sx=new float[]{20,9,6,4,3,1};
		float[] sy=new float[]{3.5f,4.9f,26.8f,12f,21.3f,12.3f};
		
		Spline s1=new Spline(sx,sy);
		
		s1.cubicSplineWith1stBC(0.6f,2.2f);
		
		float[] re=s1.cValues(new float[]{1.6f,2.4f,5f,7f,9f,15f});
		System.out.println(Arrays.toString(re));
	}*/
}
