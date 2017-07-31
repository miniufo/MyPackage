/**
 * @(#)TridiagonalSovler.java	1.0 2013.06.21
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;


/**
 * for solving tri-diagonal matrix
 *
 * @version 1.0, 2013.06.21
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TridiagonalSovler{
	
	// Prevent from construction
	private TridiagonalSovler(){}
	
	
	/**
	 * Solve the given tridiagonal system of linear equations. This method
	 * solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> = <I>b</I> where
	 * <I>A</I> is tridiagonal (<I>N</I> &gt;= 2). The form of <I>A</I> for the
	 * 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   0  ]
	 * A = [ f0  d1  e1  0  ]
	 *     [ 0   f1  d2  e2 ]
	 *     [ 0   0   f2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 2.
	 * @param  e  (input) Vector of super-diagonal elements. Length must be
	 *            <I>N</I>-1.
	 * @param  f  (input) Vector of sub-diagonal elements. Length must be
	 *            <I>N</I>-1.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solve(float[] d,float[] e,float[] f,float[] b,float[] x){
		// Verify preconditions.
		int N = d.length;
		
		if(N<2)
		throw new IllegalArgumentException("Tridiagonal.solve(): d.length = "+d.length+" illegal");
			
		if(e.length!=N-1)
		throw new IllegalArgumentException("Tridiagonal.solve(): e.length = "+e.length+" illegal");
		
		if(f.length!=N-1)
		throw new IllegalArgumentException("Tridiagonal.solve(): f.length = "+f.length+" illegal");
		
		if(b.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solve(): b.length = "+b.length+" illegal");
		
		if(x.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solve(): x.length = "+x.length+" illegal");
		
		// Working storage.
		float[] alpha = new float [N];
		float[] z = new float [N];
		
		// Elimination of sub-diagonal. alpha = new diagonal, z = new right hand
		// side.
		z[0]    =b[0];
		alpha[0]=d[0];
		
		if(alpha[0]==0.0) throw new IllegalArgumentException("Tridiagonal.solve(): Zero on diagonal");
		
		for(int i=1;i<N;i++){
			float t=f[i-1]/alpha[i-1];
			alpha[i]=d[i]-t*e[i-1];
			z[i]=b[i]-t*z[i-1];
			
			if (alpha[i]==0f) throw new IllegalArgumentException("Tridiagonal.solve(): Zero on diagonal");
		}
		
		// Back substitution.
		int Nminus1=N-1;
		x[Nminus1] =z[Nminus1]/alpha[Nminus1];
		
		for(int i=N-2;i>=0;i--) x[i]=(z[i]-e[i]*x[i+1])/alpha[i];
	}
	
	
	/**
	 * Solve the given symmetric tridiagonal system of linear equations. This
	 * method solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> =
	 * <I>b</I> where <I>A</I> is symmetric tridiagonal (<I>N</I> &gt;= 2). The
	 * form of <I>A</I> for the 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   0  ]
	 * A = [ e0  d1  e1  0  ]
	 *     [ 0   e1  d2  e2 ]
	 *     [ 0   0   e2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 2.
	 * @param  e  (input) Vector of off-diagonal elements. Length must be
	 *            <I>N</I>-1.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solveSymmetric(float[] d,float[] e,float[] b,float[] x){ solve(d,e,e,b,x);}
	
	
	/**
	 * Solve the given cyclic tridiagonal system of linear equations. This
	 * method solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> =
	 * <I>b</I> where <I>A</I> is cyclic tridiagonal (<I>N</I> &gt;= 3). The
	 * form of <I>A</I> for the 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   f3 ]
	 * A = [ f0  d1  e1  0  ]
	 *     [ 0   f1  d2  e2 ]
	 *     [ e3  0   f2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 3.
	 * @param  e  (input) Vector of super-diagonal elements. Length must be
	 *            <I>N</I>.
	 * @param  f  (input) Vector of sub-diagonal elements. Length must be
	 *            <I>N</I>.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solveCyclic(float[] d,float[] e,float[] f,float[] b,float[] x){
		// Verify preconditions.
		int N=d.length;
		if(N<3)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): d.length = "+d.length+" illegal");
		
		if(e.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): e.length = "+e.length+" illegal");
		
		if(f.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): f.length = "+f.length+" illegal");
		
		if(b.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): b.length = "+b.length+" illegal");
		
		if(x.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): x.length = "+x.length+" illegal");
		
		// Working storage.
		float[] alpha=new float[N];
		float[] zb   =new float[N];
		float[] zu   =new float[N];
		float[] w    =new float[N];

		// Elimination of sub-diagonal. alpha = new diagonal, zb = new right
		// hand side. A*q = zu.
		if(d[0]==0f||d[1]==0f)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): Zero on diagonal");
		
		zb[0]=b[0];
		
		float beta=-d[0];
		float q=1f-(e[0]*f[0])/(d[0]*d[1]);
		float abs_q_over_beta=Math.abs(q/beta);
		
		if(abs_q_over_beta<=0.5){
			
		}else if(abs_q_over_beta<1.0){
			beta *= 0.5;
		}else if(abs_q_over_beta<2.0){
			beta *= 2.0;
		}
		
		zu[0]=beta;
		alpha[0]=d[0]-beta;
		if(alpha[0]==0.0)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): Zero on diagonal");
		
		int Nminus1=N-1;
		for(int i=1;i<Nminus1;i++){
			float t=f[i-1]/alpha[i-1];
			
			alpha[i]=d[i]-t*e[i-1];
			zb[i]=b[i]-t*zb[i-1];
			zu[i]=-t*zu[i-1];
			
			if(alpha[i] == 0.0)
			throw new IllegalArgumentException("Tridiagonal.solveCyclic(): Zero on diagonal");
			
		}
		
		int Nminus2=N-2;
		float t=f[Nminus2]/alpha[Nminus2];
		alpha[Nminus1]=d[Nminus1]-e[Nminus1]*f[Nminus1]/beta-t*e[Nminus2];
		
		zb[Nminus1]=b[Nminus1]-t*zb[Nminus2];
		zu[Nminus1]=e[Nminus1]-t*zu[Nminus2];
		
		if(alpha[Nminus1]==0.0)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): Zero on diagonal");
		
		// Back substitution.
		w[Nminus1]=zu[Nminus1]/alpha[Nminus1];
		x[Nminus1]=zb[Nminus1]/alpha[Nminus1];
		
		for(int i=Nminus2;i>=0;i--){
			w[i]=(zu[i]-e[i]*w[i+1])/alpha[i];
			x[i]=(zb[i]-e[i]*x[i+1])/alpha[i];
		}
		
		// Sherman-Morrison to fix up from corner elements.
		float vw=w[0]+f[Nminus1]/beta*w[Nminus1]+1f;
		float vx=x[0]+f[Nminus1]/beta*x[Nminus1];
		
		if(vw==0.0)
		throw new IllegalArgumentException("Tridiagonal.solveCyclic(): Zero on diagonal");
		
		float vx_over_vw=vx/vw;
		
		for(int i=0;i<N;i++) x[i]-=vx_over_vw*w[i];
	}
	
	
	/**
	 * Solve the given symmetric cyclic tridiagonal system of linear equations.
	 * This method solves the general <I>N</I>-by-<I>N</I> system <I>Ax</I> =
	 * <I>b</I> where <I>A</I> is symmetric cyclic tridiagonal (<I>N</I> &gt;=
	 * 3). The form of <I>A</I> for the 4-by-4 case is:
	 * <PRE>
	 *     [ d0  e0  0   e3 ]
	 * A = [ e0  d1  e1  0  ]
	 *     [ 0   e1  d2  e2 ]
	 *     [ e3  0   e2  d3 ]
	 * </PRE>
	 *
	 * @param  d  (input) Vector of diagonal elements. Length <I>N</I> must be
	 *            &gt;= 3.
	 * @param  e  (input) Vector of off-diagonal elements. Length must be
	 *            <I>N</I>.
	 * @param  b  (input) Vector of right hand side elements. Length must be
	 *            <I>N</I>.
	 * @param  x  (output) Solution vector. Length must be <I>N</I>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any argument is the wrong length.
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the linear system cannot be solved.
	 */
	public static void solveSymmetricCyclic(float[] d,float[] e,float[] b,float[] x){
		// Verify preconditions.
		int N=d.length;
		if(N<3)
		throw new IllegalArgumentException("Tridiagonal.solveSymmetricCyclic(): d.length = "+d.length+" illegal");
		
		if(e.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solveSymmetricCyclic(): e.length = "+e.length+" illegal");
		
		if(b.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solveSymmetricCyclic(): b.length = "+b.length+" illegal");
		
		if(x.length!=N)
		throw new IllegalArgumentException("Tridiagonal.solveSymmetricCyclic(): x.length = "+x.length+" illegal");
		
		// Working storage.
		float sum=0;
		float[] alpha=new float [N];
		float[] gamma=new float [N];
		float[] delta=new float [N];
		float[] c    =new float [N];
		float[] z    =new float [N];
		
		// Factor.
		int Nminus1=N-1;
		int Nminus2=N-2;
		int Nminus3=N-3;
		
		if(d[0]==0.0)
		throw new IllegalArgumentException("Tridiagonal.solveSymmetricCyclic(): Zero on diagonal");
		
		alpha[0]=d[0];
		gamma[0]=e[0]/alpha[0];
		delta[0]=e[Nminus1]/alpha[0];
		
		sum+=alpha[0]*delta[0]*delta[0];
		
		for(int i=1;i<Nminus2;i++){
			alpha[i]=d[i]-e[i-1]*gamma[i-1];
			
			if(alpha[i]==0.0)
			throw new IllegalArgumentException("Tridiagonal.solveSymmetricCyclic(): Zero on diagonal");
			
			gamma[i]=e[i]/alpha[i];
			delta[i]=-delta[i-1]*e[i-1]/alpha[i];
			sum+=alpha[i]*delta[i]*delta[i];
		}
		
		alpha[Nminus2]= d[Nminus2]-e[Nminus3]*gamma[Nminus3];
		gamma[Nminus2]=(e[Nminus2]-e[Nminus3]*delta[Nminus3])/alpha[Nminus2];
		alpha[Nminus1]= d[Nminus1]-sum-alpha[Nminus2]*gamma[Nminus2]*gamma[Nminus2];
		
		// Update.
		z[0]=b[0];
		for(int i=1;i<Nminus1;i++) z[i]=b[i]-z[i-1]*gamma[i-1];
		sum=0.0f;
		for(int i=0;i<Nminus2;i++) sum+=delta[i]*z[i];
		z[Nminus1]=b[Nminus1]-sum-gamma[Nminus2]*z[Nminus2];
		for(int i=0;i<N;i++) c[i]=z[i]/alpha[i];
		
		// Back substitution.
		x[Nminus1]=c[Nminus1];
		x[Nminus2]=c[Nminus2]-gamma[Nminus2]*x[Nminus1];
		for(int i=Nminus3;i>=0;i--)x[i]=c[i]-gamma[i]*x[i+1]-delta[i]*x[Nminus1];
	}
	
	
	/** test
	public static void main(String[] args){
		float[] a={1,1,1,1};
		float[] b={3,3,3,3};
		float[] c={1,1,1,1};
		float[] d={4,4,4,4};
		
		float[] U=new float[4];
		
		solveCyclic(b,c,a,d,U);
		
		for(int i=0;i<4;i++) System.out.println(U[i]);
	}*/
}
