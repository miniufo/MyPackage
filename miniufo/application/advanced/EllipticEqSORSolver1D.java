/**
 * @(#)EllipticEqSORSolver1D.java	1.0 2018.04.08
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.advanced;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import miniufo.application.GeoFluidApplication;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.Print;
import miniufo.mathsphysics.TridiagonalAlg;


/**
 * A class for solving elliptic equation of the form:
 * 
 *   d     dS
 *   --( A --) = F
 *   dx    dx
 * 
 * Reference: my note on inverting 1D elliptic equation using SOR iteration
 *
 * @version 1.0, 2018.04.08
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EllipticEqSORSolver1D extends GeoFluidApplication implements Print{
	//
	private boolean print =true;	// boolean flag for printing
	private boolean setDim=false;
	private boolean setA  =false;
	
	private Variable APrime=null;	// elliptic coefficient A'
	
	private Dimension dim=null;
	
	public enum Dimension{X,Y}		// each spatial dimension is possible
	
	
	/**
	 * constructor
	 */
	public EllipticEqSORSolver1D(SpatialModel sm){ super(sm);}
	
	
	/*** getor and setor ***/
	public void setDimCombination(Dimension dim){ this.dim=dim; setDim=true;}
	
	public void setA(Variable APrime){
		this.APrime=APrime; setDim=true;
	}
	
	public void setPrinting(boolean print){ this.print=print;}
	
	
	/**
     * Inverting elliptic equation of the form:
     * d   dS    d   dS    d   dS    d   dS
     * --(A--) + --(B--) + --(B--) + --(C--) = F
     * dx  dx    dy  dx    dx  dy    dy  dy
     * for S (res) given F (frc) and first (x) and second (y) dimensions using SOR method.
     * 
     * If AC-BB > 0, a unique solution can be obtained.  Otherwise, no solution can be found.
     * 
     * @param	S		result on the left of the equation
     * @param	F		a given forcing function on the r.h.s. of the equation
     * @param	modify	whether modify coefficients to prevent from overflowing
     *
     * @exception	if r,v are not dimensionally the same
     */
	public void solve(Variable S,Variable F){
		if(F!=null) checkDimensions(S,F);
		assignSubDomainParams(S);
		
		if(print) System.out.println("Start solving elliptic equation using SOR with BCs ("+BCx+" "+BCy+")...");
		
		if(!setDim) throw new IllegalArgumentException("dimension not set");
		if(!setA  ) throw new IllegalArgumentException("elliptic coefficient A not set");
		
		List<Future<float[]>> results=new ArrayList<>();
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
		
	    if(S.isTFirst()){
	    	switch(dim){
			case X:{
				Params param=new Params(x,dx,S.getUndef(),BCx);
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					float[] Sbuf=S.getData()[l][k][j];
					float[] Fbuf=F==null?new float[x]:F.getData()[l][k][j];
					float[] Abuf=APrime.getData()[l][k][j];
					
					String info=tdef[l+tstart-1]+"\t"+zdef[k+zstart-1]/100+"hPa loops ";
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						return invertingOneSliceByTriDiagSolver(info,Abuf,Sbuf,Fbuf,param,overflow);
					}));
				}
				break;
				
			}case Y:{
				Params param=new Params(y,dy,S.getUndef(),BCy);
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					float[] Sbuf=new float[y]; float[][][] Sdata=S.getData()[l];
					float[] Fbuf=new float[y]; float[][][] Fdata=F==null?null:F.getData()[l];
					float[] Abuf=new float[y]; float[][][] Adata=APrime.getData()[l];
					
					for(int j=0;j<y;j++){
						Sbuf[j]=Sdata[k][j][i];
						Abuf[j]=Adata[k][j][i];
						if(F!=null) Fbuf[j]=Fdata[k][j][i];
					}
					
					String info=tdef[l+tstart-1]+" loops ";
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						return invertingOneSliceByTriDiagSolver(info,Abuf,Sbuf,Fbuf,param,overflow);
					}));
				}
				break;
			}	
			default: throw new IllegalArgumentException("unsupported dimension: "+dim);
			}
			
	    }else{
	    	switch(dim){
			case X:{
				Params param=new Params(x,dx,S.getUndef(),BCx);
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++){
					float[] Sbuf=new float[x]; float[][][] Sdata=S.getData()[k];
					float[] Fbuf=new float[x]; float[][][] Fdata=F==null?null:F.getData()[k];
					float[] Abuf=new float[x]; float[][][] Adata=APrime.getData()[k];
					
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						Sbuf[i]=Sdata[j][i][l];
						Abuf[i]=Adata[j][i][l];
						if(F!=null) Fbuf[i]=Fdata[j][i][l];
					}
					
					String info=tdef[l+tstart-1]+"\t"+zdef[k+zstart-1]/100+"hPa loops ";
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						return invertingOneSliceByTriDiagSolver(info,Abuf,Sbuf,Fbuf,param,overflow);
					}));
				}
				break;
				
			}case Y:{
				Params param=new Params(y,dy,S.getUndef(),BCy);
				for(int l=0;l<t;l++)
				for(int i=0;i<x;i++){
					float[] Sbuf=new float[y]; float[][][][] Sdata=S.getData();
					float[] Fbuf=new float[y]; float[][][][] Fdata=F==null?null:F.getData();
					float[] Abuf=new float[y]; float[][][][] Adata=APrime.getData();
					
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++){
						Sbuf[j]=Sdata[k][j][i][l];
						Abuf[j]=Adata[k][j][i][l];
						if(F!=null) Fbuf[j]=Fdata[k][j][i][l];
					}
					
					String info=tdef[l+tstart-1]+" loops ";
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						return invertingOneSliceByTriDiagSolver(info,Abuf,Sbuf,Fbuf,param,overflow);
					}));
				}
				break;
				
			}default: throw new IllegalArgumentException("unsupported dimention combinaiton: "+dim);
			}
		}
	    
		float[][][][] Sdata=S.getData();
    	switch(dim){
		case X:
			try{
				for(int l=0,ptr=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					float[] re=results.get(ptr++).get();
					
					if(!S.isTFirst())
					for(int i=0;i<x;i++) Sdata[k][j][i][l]=re[i];
				}
			}
			catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
			break;
			
		case Y:
			try{
				for(int l=0,ptr=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					float[] re=results.get(ptr++).get();
					
					if(S.isTFirst())
						for(int j=0;j<y;j++) Sdata[l][k][j][i]=re[j];
					else
						for(int j=0;j<y;j++) Sdata[k][j][i][l]=re[j];
				}
			}
			catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
			break;
			
		default: throw new IllegalArgumentException("unsupported dimension combinaiton: "+dim);
		}
	    
		if(print) System.out.println("Finished.");
	}
	
	
	/*** helper methods and classes ***/
	
	/**
	 * Inverting one slice (2D that constitute first and second dimensions)
	 * data using SOR iteration.
	 * 
	 * @param	info	used for printing
	 * @param	undef	undefined data value
	 * @param	A		coefficient for the first dimension
	 * @param	B		coefficient for the cross derivatives
	 * @param	C		coefficient for the second dimension
	 * @param	S		results of the SOR inversion
	 * @param	F		forcing function
	 */
	private float[] invertingOneSliceByTriDiagSolver(String info,float[] A,float[] S,float[] F,Params param,boolean[] overflow){
		int dimC=param.dimC;
		
		float delSqr=param.delSqr;
		
		BoundaryCondition dimBC=param.dimBC;
		
		if(dimBC!=BoundaryCondition.Fixed) throw new IllegalArgumentException("Only fixed boundary condition allowed");
		
		double[] a=new double[dimC-1];
		double[] b=new double[dimC  ];
		double[] c=new double[dimC-1];
		double[] d=new double[dimC  ];
		
		for(int i=1;i<dimC-1;i++){
			b[i]=-A[i-1]-A[i];
			a[i]=A[i-1];
			c[i]=A[i];
			d[i]=F[i]*delSqr;
		}
		
		b[0]=d[0]; b[dimC-1]=d[dimC-1];	// fixed boundary conditions
		
		TridiagonalAlg ta=new TridiagonalAlg(dimC);
		
		double[] U=ta.trace(a,b,c,d);
		
		for(int i=0;i<U.length;i++) S[i]=(float)U[i];
		
		double norm=cAbsMean(S);
		
		if(Double.isNaN(norm)||norm>1e9) overflow[0]=true;
		
		if(print) System.out.println(info+(overflow[0]?"   overflows!":""));
		
		return S;
	}
	
	private double cAbsMean(float[] data){
		double sum=0;
		
		for(int i=0,I=data.length;i<I;i++) sum+=Math.abs(data[i]);
		
		sum/=data.length;
		
		return sum;
	}
	
	
	/**
	 * parameters used for SOR iteration
	 */
	private static final class Params{
		//
		int dimC=0;		// grids in the first dimension
		
		float delSqr=0;	// squared increment
		
		BoundaryCondition dimBC=null;	// BC for the dimension
		
		/**
		 * constructor
		 * 
		 * @param dimC		grids in the first dimension
		 * @param delDim	increment in the first dimension
		 * @param undef		undefined data value
		 */
		public Params
		(int dimC,float delDim,float undef,BoundaryCondition dimBC){
			this.dimC =dimC;
			this.dimBC=dimBC;
			delSqr=delDim*delDim;
		}
	}
}
