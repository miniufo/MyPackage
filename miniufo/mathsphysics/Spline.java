/**
 * @(#)Spline.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import miniufo.basic.ArrayUtil;


/**
 * spline class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Spline{
	//
	private int slength=0;
	
	private boolean incre=true;
	
	private double[] sx=null;	// given x vector
	private double[] sy=null;	// given y value vector refer to the x vector
	private double[] h =null;	// step length
	
	private double[] o3=null;	// coefficient vector of 3-order term
	private double[] o2=null;	// coefficient vector of 2-order term
	private double[] o1=null;	// coefficient vector of 1-order term
	
	private TridiagonalAlg trace=null;
	
	
	/**
     * constructor
     *
     * @param	sx	given x vector
     * @param	sy	given y value vector refer to the x vector
     */
	public Spline(float[] x,float[] y){
		slength=x.length;
		
		if(slength!=y.length) throw new IllegalArgumentException("lengths mismatch");
		if(slength<2) throw new IllegalArgumentException("array's length should be larger than 1");
		
		incre=x[slength-1]>x[0];
		
		sx=new double[slength  ];
		sy=new double[slength  ];
		h =new double[slength-1];
		o1=new double[slength-1];
		o2=new double[slength-1];
		o3=new double[slength-1];
		
		if(incre)
			for(int i=0;i<slength;i++){ sx[i]=x[i]; sy[i]=y[i];}
		else
			for(int i=0;i<slength;i++){ sx[slength-1-i]=x[i]; sy[slength-1-i]=y[i];}
		
		/*** calculate heap ***/
		for(int i=0,I=slength-1;i<I;i++) h[i]=sx[i+1]-sx[i];
		
		trace=new TridiagonalAlg(slength);
	}
	
	public Spline(double[] x,double[] y){
		slength=x.length;
		
		if(slength!=y.length) throw new IllegalArgumentException("lengths mismatch");
		if(slength<2) throw new IllegalArgumentException("array's length should be larger than 1");
		
		incre=x[slength-1]>x[0];
		
		sx=new double[slength  ];
		sy=new double[slength  ];
		h =new double[slength-1];
		o1=new double[slength-1];
		o2=new double[slength-1];
		o3=new double[slength-1];
		
		if(incre)
			for(int i=0;i<slength;i++){ sx[i]=x[i]; sy[i]=y[i];}
		else
			for(int i=0;i<slength;i++){ sx[slength-1-i]=x[i]; sy[slength-1-i]=y[i];}
		
		/*** calculate heap ***/
		for(int i=0,I=slength-1;i<I;i++) h[i]=sx[i+1]-sx[i];
		
		trace=new TridiagonalAlg(slength);
	}
	
	
	/**
     * changing data y without modifying grid info x
     *
     * @param	y	y = f(x), function values of x
     */
	public void updateData(float[] y){
		if(y.length!=slength)
		throw new IllegalArgumentException("invalid length of array ("+y.length+"), should be "+slength);
		
		if(incre) for(int i=0;i<slength;i++) sy[i]=y[i];
		else for(int i=0;i<slength;i++) sy[slength-1-i]=y[i];
	}
	
	public void updateData(double[] y){
		if(y.length!=slength)
		throw new IllegalArgumentException("invalid length of array ("+y.length+"), should be "+slength);
		
		if(incre) for(int i=0;i<slength;i++) sy[i]=y[i];
		else for(int i=0;i<slength;i++) sy[slength-1-i]=y[i];
	}
	
	
	/**
     * Cubic spline interpolation with first-order derivatives
     * (complete or clamped) boundary condition.
     *
     * @param	c1	left boundary value
     * @param	c2	right boundary value
     */
	public void cubicSplineWith1stBCMatrix(double c1,double c2){
		DMatrixRMaj A=new DMatrixRMaj(slength,slength);
		DMatrixRMaj d=new DMatrixRMaj(slength,1);
		DMatrixRMaj m=new DMatrixRMaj(slength,1);
		
		for(int i=1,I=slength-1;i<I;i++){
			A.set(i,i-1,h[i-1]/(h[i-1]+h[i]));	// miu
			A.set(i,i  ,2);						// diagnal
			A.set(i,i+1,h[i  ]/(h[i-1]+h[i]));	// lambda
			d.set(i,0,6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]));
		}
		
		A.set(0,0,2.0);
		A.set(0,1,1.0);
		d.set(0,0,6.0/h[0]*((sy[1]-sy[0])/h[0]-c1));
		
		A.set(slength-1,slength-2,1.0);
		A.set(slength-1,slength-1,2.0);
		d.set(slength-1,0,6.0/h[slength-2]*(c2-(sy[slength-1]-sy[slength-2])/h[slength-2]));
		
		CommonOps_DDRM.solve(A,d,m);
		
		computeCoefficients(m.data);
	}
	
	public void cubicSplineWith1stBC(double c1,double c2){
		double[] a=new double[slength-1];
		double[] b=new double[slength  ];
		double[] c=new double[slength-1];
		double[] d=new double[slength  ];
		
		for(int i=1,I=slength-1;i<I;i++){
			a[i-1]=h[i-1]/(h[i-1]+h[i]);	// lower-band
			b[i  ]=2;						// diagnal
			c[i  ]=h[i  ]/(h[i-1]+h[i]);	// upper-band
			d[i  ]=6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]);
		}
		
		b[0]=2.0;
		c[0]=1.0;
		d[0]=6.0/h[0]*((sy[1]-sy[0])/h[0]-c1);
		
		a[slength-2]=1.0;
		b[slength-1]=2.0;
		d[slength-1]=6.0/h[slength-2]*(c2-(sy[slength-1]-sy[slength-2])/h[slength-2]);
		
		if(trace.getLength()!=slength) trace=new TridiagonalAlg(slength);
		
		computeCoefficients(trace.trace(a,b,c,d));
	}
	
	/**
     * cubic spline interpolation with second-order derivatives boundary condition
     *
     * @param	c1	left boundary value
     * @param	c2	right boundary value
     */
	public void cubicSplineWith2ndBCMatrix(double c1,double c2){
		DMatrixRMaj A=new DMatrixRMaj(slength,slength);
		DMatrixRMaj d=new DMatrixRMaj(slength,1);
		DMatrixRMaj m=new DMatrixRMaj(slength,1);
		
		for(int i=1,I=slength-1;i<I;i++){
			A.set(i,i-1,h[i-1]/(h[i-1]+h[i]));	// miu
			A.set(i,i  ,2);						// diagnal
			A.set(i,i+1,h[i  ]/(h[i-1]+h[i]));	// lambda
			d.set(i,0,6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]));
		}
		
		A.set(0,0,1.0);
		d.set(0,0,c1);
		
		A.set(slength-1,slength-1,1.0);
		d.set(slength-1,0,c2);
		
		CommonOps_DDRM.solve(A,d,m);
		
		computeCoefficients(m.data);
	}
	
	public void cubicSplineWith2ndBC(double c1,double c2){
		double[] a=new double[slength-1];
		double[] b=new double[slength  ];
		double[] c=new double[slength-1];
		double[] d=new double[slength  ];
		
		for(int i=1,I=slength-1;i<I;i++){
			a[i-1]=h[i-1]/(h[i-1]+h[i]);	// lower-band
			b[i  ]=2;						// diagnal
			c[i  ]=h[i  ]/(h[i-1]+h[i]);	// upper-band
			d[i]=6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]);
		}
		
		b[0]=1.0;
		d[0]=c1;
		
		b[slength-1]=1.0;
		d[slength-1]=c2;
		
		if(trace.getLength()!=slength) trace=new TridiagonalAlg(slength);
		
		computeCoefficients(trace.trace(a,b,c,d));
	}
	
	/**
     * cubic spline interpolation with natural boundary condition
     */
	public void cubicSplineWithNaturalBC(){ cubicSplineWith2ndBC(0,0);}
	
	/**
     * Cubic spline interpolation with periodic boundary condition.
     * 
     * Periodic boundary condition means that the 1st and 2nd derivatives at the
     * end-points are equal while the values at the end-points could be different
     * 
     * @param	step	step length between x[0] and x[x.length-1]
     */
	public void cubicSplineWithPeriodicBCMatrix(){
		DMatrixRMaj A=new DMatrixRMaj(slength,slength);
		DMatrixRMaj d=new DMatrixRMaj(slength,1);
		DMatrixRMaj m=new DMatrixRMaj(slength,1);
		
		for(int i=1,I=slength-1;i<I;i++){
			A.set(i,i-1,h[i-1]/(h[i-1]+h[i]));	// miu
			A.set(i,i  ,2);						// diagnal
			A.set(i,i+1,h[i  ]/(h[i-1]+h[i]));	// lambda
			d.set(i,0,6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]));
		}
		
		A.set(0,0,1.0);
		A.set(0,slength-1,-1.0);
		
		A.set(slength-1,1,h[0]/(h[0]+h[slength-2]));
		A.set(slength-1,slength-2,h[slength-2]/(h[slength-2]+h[0]));
		A.set(slength-1,slength-1,2.0);
		d.set(slength-1,0,6.0/(h[slength-2]+h[0])*
			((sy[1]-sy[0])/h[0]-(sy[slength-1]-sy[slength-2])/h[slength-2]));
		
		CommonOps_DDRM.solve(A,d,m);
		
		computeCoefficients(m.data);
	}
	
	public void cubicSplineWithPeriodicBC(){
		if(slength<5){
			cubicSplineWithPeriodicBCMatrix();
			return;
		}
		
		double[] a=new double[slength-2];
		double[] b=new double[slength-1];
		double[] c=new double[slength-2];
		double[] d=new double[slength-1];
		
		for(int i=0,I=slength-2;i<I;i++){
			a[i]=i==I-1?h[I]/(h[I]+h[0]):h[i+1]/(h[i+1]+h[i+2]);	// lower-band
			b[i]=2;													// diagnal
			c[i]=h[i+1]/(h[i]+h[i+1]);								// upper-band
			d[i]=6.0/(h[i]+h[i+1])*((sy[i+2]-sy[i+1])/h[i+1]-(sy[i+1]-sy[i])/h[i]);
		}
		
		double a0=h[0]/(h[0]+h[1]);
		b[slength-2]=2.0;
		d[slength-2]=6.0/(h[slength-2]+h[0])*((sy[1]-sy[0])/h[0]-(sy[slength-1]-sy[slength-2])/h[slength-2]);
		double cn=h[0]/(h[0]+h[slength-2]);
		
		if(trace.getLength()!=slength-1) trace=new TridiagonalAlg(slength-1);
		
		double[] m=new double[slength];
		
		System.arraycopy(trace.traceCyclic(a,b,c,d,a0,cn),0,m,1,slength-1);
		
		m[0]=m[slength-1];
		
		computeCoefficients(m);
	}
	
	/**
     * cubic spline interpolation with periodic boundary condition
     * 
     * @param	step	step length between x[0] and x[x.length-1]
     */
	public void cubicSplineWithNotAKnotBCMatrix(){
		DMatrixRMaj A=new DMatrixRMaj(slength,slength);
		DMatrixRMaj d=new DMatrixRMaj(slength,1);
		DMatrixRMaj m=new DMatrixRMaj(slength,1);
		
		for(int i=1,I=slength-1;i<I;i++){
			A.set(i,i-1,h[i-1]/(h[i-1]+h[i]));	// miu
			A.set(i,i  ,2);						// diagnal
			A.set(i,i+1,h[i  ]/(h[i-1]+h[i]));	// lambda
			d.set(i,0,6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]));
		}
		
		A.set(0,0,-h[1]/(h[0]+h[1]));
		A.set(0,1,1.0);
		A.set(0,2,-h[0]/(h[0]+h[1]));
		
		A.set(slength-1,slength-3,-h[slength-2]/(h[slength-3]+h[slength-2]));
		A.set(slength-1,slength-2,1.0);
		A.set(slength-1,slength-1,-h[slength-3]/(h[slength-3]+h[slength-2]));
		
		CommonOps_DDRM.solve(A,d,m);
		
		computeCoefficients(m.data);
	}
	
	public void cubicSplineWithNotAKnotBC(){
		if(slength<5){
			cubicSplineWithNotAKnotBCMatrix();
			return;
		}
		
		if(h[0]!=h[1]){
			double[] a=new double[slength-1];
			double[] b=new double[slength  ];
			double[] c=new double[slength-1];
			double[] d=new double[slength  ];
			
			for(int i=1,I=slength-1;i<I;i++){
				a[i-1]=h[i-1]/(h[i-1]+h[i]);	// lower-band
				b[i  ]=2;						// diagnal
				c[i  ]=h[i  ]/(h[i-1]+h[i]);	// upper-band
				d[i  ]=6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]);
			}
			
			double c0_c1=-h[0]/h[1];
			
			b[0]=c0_c1*a[0]+h[1]/(h[0]+h[1]);
			c[0]=2.0*c0_c1-1.0;
			d[0]=c0_c1*d[1];
			
			double an_an_1=-h[slength-2]/h[slength-3];
			
			a[slength-2]=2.0*an_an_1-1.0;
			b[slength-1]=c[slength-2]*an_an_1+h[slength-3]/(h[slength-3]+h[slength-2]);
			d[slength-1]=an_an_1*d[slength-2];
			
			if(trace.getLength()!=slength) trace=new TridiagonalAlg(slength);
			
			computeCoefficients(trace.trace(a,b,c,d));
			
		}else{
			double[] a=new double[slength-2];
			double[] b=new double[slength-1];
			double[] c=new double[slength-2];
			double[] d=new double[slength-1];
			
			for(int i=1,I=slength-1;i<I;i++){
				if(i<slength-2) a[i-1]=h[i]/(h[i]+h[i+1]);	// lower-band
				b[i-1]=2;					// diagnal
				c[i-1]=h[i]/(h[i-1]+h[i]);	// upper-band
				d[i-1]=6.0/(h[i-1]+h[i])*((sy[i+1]-sy[i])/h[i]-(sy[i]-sy[i-1])/h[i-1]);
			}
			
			double m1=h[0]*d[0]/(2*h[0]+h[1]);
			
			b[0]=h[0]/(h[0]+h[1]);
			d[0]-=2.0*m1;
			d[1]-=h[1]/(h[1]+h[2])*m1;
			a[0]=0;
			
			double an_an_1=-h[slength-2]/h[slength-3];
			
			a[slength-3]=2.0*an_an_1-1.0;
			b[slength-2]=c[slength-3]*an_an_1+h[slength-3]/(h[slength-3]+h[slength-2]);
			c[slength-3]=h[slength-2]/(h[slength-3]+h[slength-2]);
			d[slength-3]=6.0/(h[slength-3]+h[slength-2])*((sy[slength-1]-sy[slength-2])/h[slength-2]-(sy[slength-2]-sy[slength-3])/h[slength-3]);
			d[slength-2]=an_an_1*d[slength-3];
			
			if(trace.getLength()!=slength-1) trace=new TridiagonalAlg(slength-1);
			
			double[] m=new double[slength];
			
			System.arraycopy(trace.trace(a,b,c,d),0,m,1,slength-1);
			
			m[0]=m[1]; m[1]=m1;
			
			computeCoefficients(m);
		}
	}
	
	
	/**
     * calculate values according to the given x
     *
     * @param	x		given x vector
     * @param	x		given x vector
     * @param	undef	undefined value
     *
     * @return	y		result values vector
     */
	public double[] cValues(double[] x,double undef){
		double[] y=new double[x.length];
		
		cValues(x,y,undef);
		
		return y;
	}
	
	public double[] cValues(double[] x){ return cValues(x,Double.NaN);}
	
	public float[] cValues(float[] x,float undef){
		float[] y=new float[x.length];
		
		cValues(x,y,undef);
		
		return y;
	}
	
	public float[] cValues(float[] x){ return cValues(x,Float.NaN);}
	
	
	/**
     * calculate values according to the given x.
     * if x outside the sx, then fill with undef.
     *
     * @param	x		given x vector
     * @param	y		result values vector
     * @param	undef	undefined value
     */
	public void cValues(double[] x,double[] y,double undef){
		int length=x.length;
		
		if(length!=y.length)
		throw new IllegalArgumentException("array lengths mismatch");
		
		boolean checkRange=!Double.isNaN(undef);
		
		for(int i=0;i<length;i++){
			int idx=0;
			
			if(checkRange)
			if(incre){
				if(x[i]<sx[0]||x[i]>sx[slength-1]){ y[i]=undef; continue;}
			}else{
				if(x[i]>sx[0]||x[i]<sx[slength-1]){ y[i]=undef; continue;}
			}
			
			idx=ArrayUtil.getLEIdxIncre(sx,x[i]);
			
			if(idx==slength-1){ y[i]=(float)sy[slength-1]; continue;} // last index case
			
			double dx=x[i]-sx[idx];
			
			y[i]=o3[idx]*dx*dx*dx+o2[idx]*dx*dx+o1[idx]*dx+sy[idx];
		}
	}
	
	public void cValues(float[] x,float[] y,float undef){
		// exactly the same as the above cValues(double[],double[], double)
		int length=x.length;
		
		if(length!=y.length)
		throw new IllegalArgumentException("array lengths mismatch");
		
		boolean checkRange=!Double.isNaN(undef);
		
		for(int i=0;i<length;i++){
			int idx=0;
			
			if(checkRange)
			if(incre){
				if(x[i]<sx[0]||x[i]>sx[slength-1]){ y[i]=undef; continue;}
			}else{
				if(x[i]>sx[0]||x[i]<sx[slength-1]){ y[i]=undef; continue;}
			}
			
			idx=ArrayUtil.getLEIdxIncre(sx,x[i]);
			
			if(idx==slength-1){ y[i]=(float)sy[slength-1]; continue;} // last index case
			
			double dx=x[i]-sx[idx];
			
			y[i]=(float)(o3[idx]*dx*dx*dx+o2[idx]*dx*dx+o1[idx]*dx+sy[idx]);
		}
	}
	
	public void cValues(double[] x,double[] y){ cValues(x,y,Double.NaN);}
	
	public void cValues(float[] x,float[] y){ cValues(x,y,Float.NaN);}
	
	
	/*** helper methods ***/
	private void computeCoefficients(double[] m){
		/*** calculate o3,o2,o1 ***/
		for(int i=0,I=slength-1;i<I;i++){
			o3[i]=(m[i+1]-m[i])/(6.0*h[i]);
			o2[i]=m[i]/2.0;
			o1[i]=(sy[i+1]-sy[i])/h[i]-h[i]/6.0*(m[i+1]+2.0*m[i]);
		}
	}
	
	
	/** test
	public static void main(String[] arg){
		float[] sx=ArrayUtil.newMonotonousArray(2000,1);
		//float[] sx=new float[]{20,9,6,4,3,1};
		float[] sy=ArrayUtil.newGaussianRandomSeries(2000,0,2);
		
		TicToc.tic("spline 1");
		Spline s1=new Spline(sx,sy);
		s1.cubicSplineWithNotAKnotBCMatrix();
		float[] re1=s1.cValues(new float[]{1.6f,2.4f,5f,7f,9f,15f});
		TicToc.toc(TimeUnit.MILLISECONDS);
		
		TicToc.tic("spline 2");
		Spline s2=new Spline(sx,sy);
		s2.cubicSplineWithNotAKnotBC();
		float[] re2=s2.cValues(new float[]{1.6f,2.4f,5f,7f,9f,15f});
		TicToc.toc(TimeUnit.MILLISECONDS);
		
		System.out.println(Arrays.toString(re1));
		System.out.println(Arrays.toString(re2));
	}*/
}
