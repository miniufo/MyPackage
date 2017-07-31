/**
 * @(#)EllipticEquationSolver.java	1.0 2015.03.12
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
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;


/**
 * a class for solving elliptic equation
 * Reference: my note on inverting elliptic equation using SOR iteration
 *
 * @version 1.0, 2015.03.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EllipticEqSORSolver extends GeoFluidApplication implements Print{
	//
	private int mxLoopCount=5000;	// max loop of the iteration
	
	private boolean print  =true;	// boolean flag for printing
	private boolean setComb=false;
	private boolean setABC =false;
	
	private double tolerance=1e-6f;	// tolerance for the iteration
	
	private Variable APrime=null;	// elliptic coefficient A'
	private Variable BPrime=null;	// elliptic coefficient B'
	private Variable CPrime=null;	// elliptic coefficient C'
	
	private DimCombination dimComb=null;
	
	public enum DimCombination{XY,YZ}	// horizontal or meridional planes
	
	
	/**
	 * constructor
	 */
	public EllipticEqSORSolver(SpatialModel sm){ super(sm);}
	
	
	/*** getor and setor ***/
	public void setDimCombination(DimCombination dimComb){ this.dimComb=dimComb; setComb=true;}
	
	public void setABC(Variable APrime,Variable BPrime,Variable CPrime){
		this.APrime=APrime;
		this.BPrime=BPrime;
		this.CPrime=CPrime; setABC=true;
	}
	
	public void setPrinting(boolean print){ this.print=print;}
	
	public void setMaxLoopCount(int count){
		if(count<1)
			throw new IllegalArgumentException("count should be greater than 1");
		
		mxLoopCount=count;
	}
	
	public void setTolerance(double tol){
		if(tol<=0)
			throw new IllegalArgumentException("tolerance should be positive");
		
		tolerance=tol;
	}
	
	
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
	public void solve(Variable S,Variable F,boolean modify){
		if(F!=null) checkDimensions(S,F);
		assignSubDomainParams(S);
		
		if(print) System.out.println("Start solving elliptic equation using SOR with BCs ("+BCx+" "+BCy+")...");
		
		if(!setComb) throw new IllegalArgumentException("dimension combination not set");
		if(!setABC ) throw new IllegalArgumentException("elliptic coefficients not set");
		
		List<Future<float[][]>> results=new ArrayList<>(y*(dimComb==DimCombination.XY?x:z));
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<float[][]> cs=new ExecutorCompletionService<>(es);
		
	    if(S.isTFirst()){
	    	switch(dimComb){
			case XY:{
				Params param=new Params(x,y,dx,dy,S.getUndef(),BCx,BCy);
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++){
					float[][] Sbuf=S.getData()[l][k];
					float[][] Fbuf=F==null?new float[y][x]:F.getData()[l][k];
					float[][] Abuf=APrime.getData()[l][k];
					float[][] Bbuf=BPrime.getData()[l][k];
					float[][] Cbuf=CPrime.getData()[l][k];
					
					String info=tdef[l+tstart-1]+"\t"+zdef[k+zstart-1]/100+"hPa loops ";
					final int ll=l;
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						float[][] re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
						if(modify&&overflow[0]){
							modifyCoeffAC(Abuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
							re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
							if(modify&&overflow[0]){
								modifyCoeffD(Abuf,Bbuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
								re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
							}
						}
						return re;
					}));
				}
				break;
				
			}case YZ:{
				Params param=new Params(y,z,dy,dz,S.getUndef(),BCy,BCz);
				for(int l=0;l<t;l++)
				for(int i=0;i<x;i++){
					float[][] Sbuf=new float[z][y]; float[][][] Sdata=S.getData()[l];
					float[][] Fbuf=new float[z][y]; float[][][] Fdata=F==null?null:F.getData()[l];
					float[][] Abuf=new float[z][y]; float[][][] Adata=APrime.getData()[l];
					float[][] Bbuf=new float[z][y]; float[][][] Bdata=BPrime.getData()[l];
					float[][] Cbuf=new float[z][y]; float[][][] Cdata=CPrime.getData()[l];
					
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++){
						Sbuf[k][j]=Sdata[k][j][i];
						Abuf[k][j]=Adata[k][j][i];
						Bbuf[k][j]=Bdata[k][j][i];
						Cbuf[k][j]=Cdata[k][j][i];
						if(F!=null) Fbuf[k][j]=Fdata[k][j][i];
					}
					
					String info=tdef[l+tstart-1]+" loops ";
					final int ll=l;
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						float[][] re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
						if(modify&&overflow[0]){
							modifyCoeffAC(Abuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
							re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
							if(modify&&overflow[0]){
								modifyCoeffD(Abuf,Bbuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
								re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
							}
						}
						return re;
					}));
				}
				break;
			}	
			default: throw new IllegalArgumentException("unsupported dimention combinaiton: "+dimComb);
			}
			
	    }else{
	    	switch(dimComb){
			case XY:{
				Params param=new Params(x,y,dx,dy,S.getUndef(),BCx,BCy);
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++){
					float[][] Sbuf=new float[y][x]; float[][][] Sdata=S.getData()[k];
					float[][] Fbuf=new float[y][x]; float[][][] Fdata=F==null?null:F.getData()[k];
					float[][] Abuf=new float[y][x]; float[][][] Adata=APrime.getData()[k];
					float[][] Bbuf=new float[y][x]; float[][][] Bdata=BPrime.getData()[k];
					float[][] Cbuf=new float[y][x]; float[][][] Cdata=CPrime.getData()[k];
					
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						Sbuf[j][i]=Sdata[j][i][l];
						Abuf[j][i]=Adata[j][i][l];
						Bbuf[j][i]=Bdata[j][i][l];
						Cbuf[j][i]=Cdata[j][i][l];
						if(F!=null) Fbuf[j][i]=Fdata[j][i][l];
					}
					
					String info=tdef[l+tstart-1]+"\t"+zdef[k+zstart-1]/100+"hPa loops ";
					final int ll=l;
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						float[][] re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
						if(modify&&overflow[0]){
							modifyCoeffAC(Abuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
							re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
							if(modify&&overflow[0]){
								modifyCoeffD(Abuf,Bbuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
								re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
							}
						}
						return re;
					}));
				}
				break;
				
			}case YZ:{
				Params param=new Params(y,z,dy,dz,S.getUndef(),BCy,BCz);
				for(int l=0;l<t;l++)
				for(int i=0;i<x;i++){
					float[][] Sbuf=new float[z][y]; float[][][][] Sdata=S.getData();
					float[][] Fbuf=new float[z][y]; float[][][][] Fdata=F==null?null:F.getData();
					float[][] Abuf=new float[z][y]; float[][][][] Adata=APrime.getData();
					float[][] Bbuf=new float[z][y]; float[][][][] Bdata=BPrime.getData();
					float[][] Cbuf=new float[z][y]; float[][][][] Cdata=CPrime.getData();
					
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++){
						Sbuf[k][j]=Sdata[k][j][i][l];
						Abuf[k][j]=Adata[k][j][i][l];
						Bbuf[k][j]=Bdata[k][j][i][l];
						Cbuf[k][j]=Cdata[k][j][i][l];
						if(F!=null) Fbuf[k][j]=Fdata[k][j][i][l];
					}
					
					String info=tdef[l+tstart-1]+" loops ";
					final int ll=l;
					results.add(cs.submit(()->{
						boolean[] overflow=new boolean[]{false};
						float[][] re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
						if(modify&&overflow[0]){
							modifyCoeffAC(Abuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
							re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
							if(modify&&overflow[0]){
								modifyCoeffD(Abuf,Bbuf,Cbuf,ll); overflow[0]=false; resetInner(Sbuf);
								re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Sbuf,Fbuf,param,overflow);
								if(overflow[0]) System.out.println(info+" still overflow after modification");
							}
						}
						return re;
					}));
				}
				break;
				
			}default: throw new IllegalArgumentException("unsupported dimention combinaiton: "+dimComb);
			}
		}
	    
		float[][][][] Sdata=S.getData();
    	switch(dimComb){
		case XY:
			try{
				for(int l=0,ptr=0;l<t;l++)
				for(int k=0;k<z;k++){
					float[][] re=results.get(ptr++).get();
					
					if(!S.isTFirst())
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++) Sdata[k][j][i][l]=re[j][i];
				}
			}
			catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
			break;
			
		case YZ:
			try{
				for(int l=0,ptr=0;l<t;l++)
				for(int i=0;i<x;i++){
					float[][] re=results.get(ptr++).get();
					
					if(S.isTFirst())
						for(int k=0;k<z;k++)
						for(int j=0;j<y;j++) Sdata[l][k][j][i]=re[k][j];
					else
						for(int k=0;k<z;k++)
						for(int j=0;j<y;j++) Sdata[k][j][i][l]=re[k][j];
				}
			}
			catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
			break;
			
		default: throw new IllegalArgumentException("unsupported dimention combinaiton: "+dimComb);
		}
	    
		if(print) System.out.println("Finished.");
	}
	
	public void solve(Variable S,Variable F){ solve(S,F,false);}
	
	
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
	private float[][] invertingOneSliceBySOR(String info,float[][] A,float[][] B,float[][] C,float[][] S,float[][] F,Params param,boolean[] overflow){
		int loop=0;	// current loop
		int dim1C=param.dim1C;
		int dim2C=param.dim2C;
		
		float ratioQtr=param.ratioQtr;
		float undef   =param.undef;
		float optArg  =param.optArg;
		float delD2Sqr=param.delD2Sqr;
		float ratioSqr=param.ratioSqr;
		
		double convSpd=0;	// convergent speed = delta norm(S) / norm(S)
		
		BoundaryCondition dim1BC=param.dim1BC;
		BoundaryCondition dim2BC=param.dim2BC;
		
		double normPrev=Double.MAX_VALUE;
		
		float[][] R=new float[dim2C][dim1C];
		
		do{
			// process the expand boundary condition i.e., expanding the inner grid to boundaries
			if(dim2BC==BoundaryCondition.Expanded){
				if(dim1BC==BoundaryCondition.Periodic){
					for(int i=0;i<x;i++){ S[0][i]=S[1][i]; S[y-1][i]=S[y-2][i];} // south and north boundaries
					
				}else{
					for(int i=1;i<x-1;i++){ S[0][i]=S[1][i]; S[y-1][i]=S[y-2][i];} // south and north boundaries
					for(int j=1;j<y-1;j++){ S[j][0]=S[j][1]; S[j][x-1]=S[j][x-2];} // east and west boundaries
					
					S[0  ][0  ]=S[1  ][1  ]; // southwest corner
					S[0  ][x-1]=S[1  ][x-2]; // southeast corner
					S[y-1][0  ]=S[y-2][1  ]; // northwest corner
					S[y-1][x-1]=S[y-2][x-2]; // northeast corner
				}
			}
			
			for(int j=1,J=dim2C-1;j<J;j++){
				if(dim1BC==BoundaryCondition.Periodic){
					// for the west boundary SOR i.e., i==0
					if(F[j][0]!=undef&&A[j][1]!=undef&&A[j][0]!=undef&&C[j+1][0]!=undef&&C[j][0]!=undef&&
					B[j][1]!=undef&&B[j][x-1]!=undef&&B[j+1][0]!=undef&&B[j-1][0]!=undef){
						R[j][0]=(
							(
								A[j][1  ]*(S[j  ][1  ]-S[j  ][0  ])-
								A[j][0  ]*(S[j  ][0  ]-S[j  ][x-1])
							)*ratioSqr+(
								B[j][  1]*(S[j+1][1  ]-S[j-1][1  ])-
								B[j][x-1]*(S[j+1][x-1]-S[j-1][x-1])
							)*ratioQtr+(
								B[j+1][0]*(S[j+1][1  ]-S[j+1][x-1])-
								B[j-1][0]*(S[j-1][1  ]-S[j-1][x-1])
							)*ratioQtr+(
								C[j+1][0]*(S[j+1][0  ]-S[j  ][0  ])-
								C[j  ][0]*(S[j  ][0  ]-S[j-1][0  ])
							)
						)-F[j][0]*delD2Sqr;
						
						R[j][0]*=optArg/((A[j][1]+A[j][0])*ratioSqr+(C[j+1][0]+C[j][0]));
						
						S[j][0]+=R[j][0];
					}
				}
				
				// process inner region SOR
				for(int i=1,I=dim1C-1;i<I;i++){
					if(F[j][i]==undef||A[j][i+1]==undef||A[j][i]==undef||C[j+1][i]==undef||C[j][i]==undef||
					B[j][i+1]==undef||B[j][i-1]==undef||B[j+1][i]==undef||B[j-1][i]==undef) continue;
					
					R[j][i]=(
						(
							A[j][i+1]*(S[j  ][i+1]-S[j  ][i  ])-
							A[j][i  ]*(S[j  ][i  ]-S[j  ][i-1])
						)*ratioSqr+(
							B[j][i+1]*(S[j+1][i+1]-S[j-1][i+1])-
							B[j][i-1]*(S[j+1][i-1]-S[j-1][i-1])
						)*ratioQtr+(
							B[j+1][i]*(S[j+1][i+1]-S[j+1][i-1])-
							B[j-1][i]*(S[j-1][i+1]-S[j-1][i-1])
						)*ratioQtr+(
							C[j+1][i]*(S[j+1][i  ]-S[j  ][i  ])-
							C[j  ][i]*(S[j  ][i  ]-S[j-1][i  ])
						)
					)-F[j][i]*delD2Sqr;
					
					R[j][i]*=optArg/((A[j][i+1]+A[j][i])*ratioSqr+(C[j+1][i]+C[j][i]));
					
					S[j][i]+=R[j][i];
				}
				
				if(dim1BC==BoundaryCondition.Periodic){
					// for the east boundary SOR i.e., i==x-1
					if(F[j][x-1]!=undef&&A[j][0]!=undef&&A[j][x-1]!=undef&&C[j+1][x-1]!=undef&&C[j][x-1]!=undef&&
					B[j][0]!=undef&&B[j][x-2]!=undef&&B[j+1][x-1]!=undef&&B[j-1][x-1]!=undef){
						R[j][x-1]=(
							(
								A[j  ][0  ]*(S[j  ][0  ]-S[j  ][x-1])-
								A[j  ][x-1]*(S[j  ][x-1]-S[j  ][x-2])
							)*ratioSqr+(
								B[j  ][0  ]*(S[j+1][0  ]-S[j-1][0  ])-
								B[j  ][x-2]*(S[j+1][x-2]-S[j-1][x-2])
							)*ratioQtr+(
								B[j+1][x-1]*(S[j+1][0  ]-S[j+1][x-2])-
								B[j-1][x-1]*(S[j-1][0  ]-S[j-1][x-2])
							)*ratioQtr+(
								C[j+1][x-1]*(S[j+1][x-1]-S[j  ][x-1])-
								C[j  ][x-1]*(S[j  ][x-1]-S[j-1][x-1])
							)
						)-F[j][x-1]*delD2Sqr;
						
						R[j][x-1]*=optArg/((A[j][0]+A[j][x-1])*ratioSqr+(C[j+1][x-1]+C[j][x-1]));
						
						S[j][x-1]+=R[j][x-1];
					}
				}
			}
			
			double norm=cAbsMean(S);
			
			if(Double.isNaN(norm)||norm>1e9){
				overflow[0]=true;
				break;
			}
			
			convSpd=Math.abs(norm-normPrev)/normPrev;
			
			if(convSpd<tolerance||loop>mxLoopCount) break;
			
			normPrev=norm; loop++;
			
		}while(true);
		
		if(print) System.out.println(info+String.format("%4d",loop)+" and tolerance is "+convSpd+(overflow[0]?"   overflows!":""));
		
		return S;
	}
	
	private double cAbsMean(float[][] data){
		double sum=0;
		
		for(int j=0,J=data.length;j<J;j++)
		for(int i=0,I=data[0].length;i<I;i++) sum+=Math.abs(data[j][i]);
		
		sum/=data.length*data[0].length;
		
		return sum;
	}
	
	private void resetInner(float[][] S){
		for(int j=1,J=S.length-1;j<J;j++)
		for(int i=1,I=S[0].length-1;i<I;i++) S[j][i]=0;
	}
	
	
	/**
	 * used for modifying the coefficients A, B, and C
	 */
	private void averageNeighboring(float[][] data){
		int z=data.length,y=data[0].length;
		
		float threshold=0;
		
		float[][] tmp=new float[z][y];
		
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++) tmp[k][j]=data[k][j];
		
		for(int k=1,K=z-1;k<K;k++)
		for(int j=1;j<y;j++) if(data[k][j]!=undef&&data[k][j]<=threshold){
			int idxJ=0,idxK=0,count=0;
			
			double sum=0;
			
			idxJ=j-1; idxK=k-1;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			idxJ=j-1; idxK=k;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			idxJ=j-1; idxK=k+1;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			idxJ=j; idxK=k-1;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			idxJ=j; idxK=k+1;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			idxJ=j+1; idxK=k-1;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			idxJ=j+1; idxK=k;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			idxJ=j+1; idxK=k+1;
			if(idxJ>=0&&idxJ<y&&idxK>=0&&idxK<z&&data[idxK][idxJ]!=undef&&data[idxK][idxJ]>threshold){
				sum+=data[idxK][idxJ]; count++;
			}
			
			if(count!=0) tmp[k][j]=(float)(sum/count);
		}
		
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++) data[k][j]=tmp[k][j];
	}
	
	private boolean hasNegative(float[][] data){
		int z=data.length,y=data[0].length;
		
		float threshold=0;
		
		for(int k=1,K=z-1;k<K;k++)
		for(int j=1;j<y;j++) if(data[k][j]!=undef&&data[k][j]<=threshold) return true;
		
		return false;
	}
	
	
	/**
	 * Modify the stability coefficients so that all values are larger than a threshold.
	 * 
	 * @param	v			a given variable, usually stability coefficients A, B, and C
	 * @param	threshold	values of the variable should be larger than this threshold	
	 */
	private void modifyCoeffAC(float[][] Adata,float[][] Cdata,int l){
		boolean modified=false;
		while(hasNegative(Adata)){ averageNeighboring(Adata); modified=true;}
		if(print&&modified) System.out.println("                          A has been modified");
		
		modified=false;
		while(hasNegative(Cdata)){ averageNeighboring(Cdata); modified=true;}
		if(print&&modified) System.out.println("                          C has been modified");
		
		float[][][][] A=APrime.getData();
		float[][][][] C=CPrime.getData();
		
		if(APrime.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				A[l][k][j][i]=Adata[k][j];
				C[l][k][j][i]=Cdata[k][j];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				A[k][j][i][l]=Adata[k][j];
				C[k][j][i][l]=Cdata[k][j];
			}
		}
	}
	
	private void modifyCoeffD(float[][] Adata,float[][] Bdata,float[][] Cdata,int l){
		int z=Adata.length,y=Adata[0].length;
		
		boolean modified=false;
		
		boolean[][] Acorr=new boolean[z][y];
		boolean[][] Ccorr=new boolean[z][y];
		
		float[][] tmpA=new float[z][y];
		float[][] tmpC=new float[z][y];
		
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++){
			tmpA[k][j]=Adata[k][j];
			tmpC[k][j]=Cdata[k][j];
		}
		
		for(int k=1,K=z-1;k<K;k++)
		for(int j=1,J=y-1;j<J;j++) if(Adata[k][j+1]!=undef&&Bdata[k][j]!=undef&&Cdata[k+1][j]!=undef){
			float D=Adata[k][j+1]*Cdata[k+1][j]-Bdata[k][j]*Bdata[k][j];
			
			if(D<=0){
				float scale=(float)Math.sqrt((Bdata[k][j]/Adata[k][j+1])*1.1f*(Bdata[k][j]/Cdata[k+1][j]));
				
				if(Float.isNaN(scale)){
					System.out.println(Adata[k][j+1]+"\t"+Bdata[k][j]+"\t"+Cdata[k+1][j]+"\t"+k+"\t"+j);
					throw new IllegalArgumentException("invalid scale");
				}
				
				Adata[k][j+1]=tmpA[k][j+1]*scale; Acorr[k][j+1]=true;
				Cdata[k+1][j]=tmpC[k+1][j]*scale; Ccorr[k+1][j]=true;
				
				modified=true;
			}
		}
		
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++){
			tmpA[k][j]=Adata[k][j];
			tmpC[k][j]=Cdata[k][j];
		}
		
		for(int k=1,K=z-1;k<K;k++)
		for(int j=1,J=y-1;j<J;j++){
			int count=0; double sum=0;
			if(Acorr[k-1][j-1]||Acorr[k-1][j]||Acorr[k-1][j+1]||Acorr[k][j-1]||Acorr[k][j+1]||Acorr[k+1][j-1]||Acorr[k+1][j]||Acorr[k+1][j+1]){
				if(tmpA[k-1][j-1]!=undef){ sum+=tmpA[k-1][j-1]; count++;}
				if(tmpA[k-1][j  ]!=undef){ sum+=tmpA[k-1][j  ]; count++;}
				if(tmpA[k-1][j+1]!=undef){ sum+=tmpA[k-1][j+1]; count++;}
				if(tmpA[k  ][j-1]!=undef){ sum+=tmpA[k  ][j-1]; count++;}
				if(tmpA[k  ][j+1]!=undef){ sum+=tmpA[k  ][j+1]; count++;}
				if(tmpA[k+1][j-1]!=undef){ sum+=tmpA[k+1][j-1]; count++;}
				if(tmpA[k+1][j  ]!=undef){ sum+=tmpA[k+1][j  ]; count++;}
				if(tmpA[k+1][j+1]!=undef){ sum+=tmpA[k+1][j+1]; count++;}
			}
			if(count!=0) Adata[k][j]=(float)(sum/count);
			
			count=0; sum=0;
			if(Ccorr[k-1][j-1]||Ccorr[k-1][j]||Ccorr[k-1][j+1]||Ccorr[k][j-1]||Ccorr[k][j+1]||Ccorr[k+1][j-1]||Ccorr[k+1][j]||Ccorr[k+1][j+1]){
				if(tmpC[k-1][j-1]!=undef){ sum+=tmpC[k-1][j-1]; count++;}
				if(tmpC[k-1][j  ]!=undef){ sum+=tmpC[k-1][j  ]; count++;}
				if(tmpC[k-1][j+1]!=undef){ sum+=tmpC[k-1][j+1]; count++;}
				if(tmpC[k  ][j-1]!=undef){ sum+=tmpC[k  ][j-1]; count++;}
				if(tmpC[k  ][j+1]!=undef){ sum+=tmpC[k  ][j+1]; count++;}
				if(tmpC[k+1][j-1]!=undef){ sum+=tmpC[k+1][j-1]; count++;}
				if(tmpC[k+1][j  ]!=undef){ sum+=tmpC[k+1][j  ]; count++;}
				if(tmpC[k+1][j+1]!=undef){ sum+=tmpC[k+1][j+1]; count++;}
			}
			if(count!=0) Cdata[k][j]=(float)(sum/count);
		}
		
		float[][][][] A=APrime.getData();
		float[][][][] C=CPrime.getData();
		
		if(APrime.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				A[l][k][j][i]=Adata[k][j];
				C[l][k][j][i]=Cdata[k][j];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				A[k][j][i][l]=Adata[k][j];
				C[k][j][i][l]=Cdata[k][j];
			}
		}
		
		if(print&&modified) System.out.println("                          D has been modified");
	}
	
	
	/**
	 * parameters used for SOR iteration
	 */
	private static final class Params{
		//
		int dim1C=0;	// grids in the first dimension
		int dim2C=0;	// grids in the second dimension
		
		float ratioQtr=0;	// ratio divided by 4
		float ratioSqr=0;	// square of the ratio
		float delD2Sqr=0;	// square of the increment in the second dimension
		float optArg=0;		// optimal argument for the SOR
		float undef=0;		// undefined value
		
		BoundaryCondition dim1BC=null;	// BC for the first dimension
		BoundaryCondition dim2BC=null;	// BC for the second dimension
		
		/**
		 * constructor
		 * 
		 * @param dim1C		grids in the first dimension
		 * @param dim2C		grids in the second dimension
		 * @param delDim1	increment in the first dimension
		 * @param delDim2	increment in the second dimension
		 * @param undef		undefined data value
		 */
		public Params
		(int dim1C,int dim2C,float delDim1,float delDim2,float undef,BoundaryCondition dim1BC,BoundaryCondition dim2BC){
			this.dim1C =dim1C;
			this.dim2C =dim2C;
			this.undef =undef;
			this.dim1BC=dim1BC;
			this.dim2BC=dim2BC;
			
			ratioQtr=delDim2/delDim1;
			ratioSqr=ratioQtr*ratioQtr;
			ratioQtr/=4f;
			delD2Sqr=delDim2*delDim2;
			
			float epsilon=(float)(pow(sin(PI/(2*dim1C+2)),2)+pow(sin(PI/(2*dim2C+2)),2));
			optArg=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		}
	}
}
