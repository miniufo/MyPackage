/**
 * @(#)GeneralEllipticEqSORSolver.java	1.0 2017.11.03
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
 * A class for solving elliptic equation in the general form:
 * 
 *   A Sxx + B Sxy + C Syy + D Sx + E Sy + F S + G = 0
 *   
 * Reference: my note on inverting elliptic equation using SOR iteration
 *
 * @version 1.0, 2017.11.03
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GeneralEllipticEqSORSolver extends GeoFluidApplication implements Print{
	//
	private int mxLoopCount=5000;	// max loop of the iteration
	
	private boolean print  =true;	// boolean flag for printing
	private boolean setComb=false;
	private boolean setCoef=false;
	
	private double tolerance=1e-6f;	// tolerance for the iteration
	
	private Variable A=null;	// elliptic coefficient A
	private Variable B=null;	// elliptic coefficient B
	private Variable C=null;	// elliptic coefficient C
	private Variable D=null;	// elliptic coefficient D
	private Variable E=null;	// elliptic coefficient E
	private Variable F=null;	// elliptic coefficient F
	private Variable G=null;	// elliptic coefficient G
	
	private DimCombination dimComb=null;
	
	public enum DimCombination{XY,YZ}	// horizontal or meridional planes
	
	
	/**
	 * constructor
	 */
	public GeneralEllipticEqSORSolver(SpatialModel sm){ super(sm);}
	
	
	/*** getor and setor ***/
	public void setDimCombination(DimCombination dimComb){ this.dimComb=dimComb; setComb=true;}
	
	public void setCoefficients(Variable A,Variable B,Variable C,Variable D,Variable E,Variable F,Variable G){
		this.A=A; this.D=D; this.G=G;
		this.B=B; this.E=E;
		this.C=C; this.F=F;
		setCoef=true;
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
     * 
     * A Sxx + B Sxy + C Syy + D Sx + E Sy + F S + G = 0
     * 
     * for S (res) given the first (x) and second (y) dimensions using SOR method.
     * 
     * If 4AC-BB > 0, a unique solution can be obtained.  Otherwise, no solution can be found.
     * 
     * @param	S		result on the left of the equation
     * @param	modify	whether modify coefficients to prevent from overflowing
     *
     * @exception	if r,v are not dimensionally the same
     */
	public void solve(Variable S,boolean modify){
		assignSubDomainParams(S);
		
		if(print) System.out.println("Start solving elliptic equation using SOR with BCs ("+BCx+" "+BCy+")...");
		
		if(!setComb) throw new IllegalArgumentException("dimension combination not set");
		if(!setCoef ) throw new IllegalArgumentException("elliptic coefficients not set");
		
		List<Future<float[][]>> results=new ArrayList<>(y*(dimComb==DimCombination.XY?x:z));
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<float[][]> cs=new ExecutorCompletionService<>(es);
		
	    if(S.isTFirst()){
	    	switch(dimComb){
			case XY:{
				throw new IllegalArgumentException("unsupported now");
				
			}case YZ:{
				Params param=new Params(y,z,dy,dz,S.getUndef(),BCy,BCz);
				for(int l=0;l<t;l++)
				for(int i=0;i<x;i++){
					float[][] zero=new float[z][y];
					
					float[][] Sbuf=new float[z][y]; float[][][] Sdata=S==null?null:S.getData()[l];
					float[][] Abuf=new float[z][y]; float[][][] Adata=A==null?null:A.getData()[l];
					float[][] Bbuf=new float[z][y]; float[][][] Bdata=B==null?null:B.getData()[l];
					float[][] Cbuf=new float[z][y]; float[][][] Cdata=C==null?null:C.getData()[l];
					float[][] Dbuf=new float[z][y]; float[][][] Ddata=D==null?null:D.getData()[l];
					float[][] Ebuf=new float[z][y]; float[][][] Edata=E==null?null:E.getData()[l];
					float[][] Fbuf=new float[z][y]; float[][][] Fdata=F==null?null:F.getData()[l];
					float[][] Gbuf=new float[z][y]; float[][][] Gdata=G==null?null:G.getData()[l];
					
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
						float[][] re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Dbuf,Ebuf,Fbuf,Gbuf,Sbuf,param,overflow);
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
				throw new IllegalArgumentException("unsupported now");
				
			}case YZ:{
				Params param=new Params(y,z,dy,dz,S.getUndef(),BCy,BCz);
				for(int l=0;l<t;l++)
				for(int i=0;i<x;i++){
					float[][] Sbuf=new float[z][y]; float[][][][] Sdata=S.getData();
					float[][] Abuf=new float[z][y]; float[][][][] Adata=A==null?null:A.getData();
					float[][] Bbuf=new float[z][y]; float[][][][] Bdata=B==null?null:B.getData();
					float[][] Cbuf=new float[z][y]; float[][][][] Cdata=C==null?null:C.getData();
					float[][] Dbuf=new float[z][y]; float[][][][] Ddata=D==null?null:D.getData();
					float[][] Ebuf=new float[z][y]; float[][][][] Edata=E==null?null:E.getData();
					float[][] Fbuf=new float[z][y]; float[][][][] Fdata=F==null?null:F.getData();
					float[][] Gbuf=new float[z][y]; float[][][][] Gdata=G==null?null:G.getData();
					
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
						float[][] re=invertingOneSliceBySOR(info,Abuf,Bbuf,Cbuf,Dbuf,Ebuf,Fbuf,Gbuf,Sbuf,param,overflow);
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
	
	public void solve(Variable S){ solve(S,false);}
	
	
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
	private float[][] invertingOneSliceBySOR(String info,float[][] A,float[][] B,float[][] C,float[][] D,float[][] E,float[][] F,float[][] G,float[][] S,Params param,boolean[] overflow){
		int loop=0;	// current loop
		int dim1C=param.dim1C;
		int dim2C=param.dim2C;
		
		float ratio =param.ratio;
		float undef =param.undef;
		float optArg=param.optArg;
		float delx1 =param.delx1;
		float delx2 =param.delx2;
		
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
					if(F[j][0]==undef&&A[j][1]==undef&&A[j][0]==undef&&C[j+1][0]==undef&&C[j][0]==undef&&
					B[j][1]==undef&&B[j][x-1]==undef&&B[j+1][0]==undef&&B[j-1][0]==undef) continue;
					
					R[j][0]=A[j][0]*(S[j][1]+S[j][x-1]-2f*S[j][0])*ratio+
							B[j][0]*(S[j+1][1]+S[j-1][x-1]-S[j-1][1]-S[j+1][x-1])/4f+
							C[j][0]*(S[j+1][0]+S[j-1][0]-2f*S[j][0])+
							D[j][0]*(S[j][1]-S[j][x-1])*delx2/2f+
							E[j][0]*(S[j][0]-S[j][0])*delx1/2f+
							(F[j][0]*S[j][0]+G[j][0])*delx1*delx2;
					
					R[j][0]*=optArg/(A[j][0]*ratio+C[j][0]/ratio)/2f;
					
					S[j][0]+=R[j][0];
				}
				
				// process inner region SOR
				for(int i=1,I=dim1C-1;i<I;i++){
					if(A[j][i]==undef||B[j][i]==undef||C[j][i]==undef||D[j][i]==undef||
					E[j][i]==undef||F[j][i]==undef||G[j][i]==undef) continue;
					
					R[j][i]=A[j][i]*(S[j][i+1]+S[j][i-1]-2f*S[j][i])*ratio+
							B[j][i]*(S[j+1][i+1]+S[j-1][i-1]-S[j-1][i+1]-S[j+1][i-1])/4f+
							C[j][i]*(S[j+1][i]+S[j-1][i]-2f*S[j][i])+
							D[j][i]*(S[j][i+1]-S[j][i-1])*delx2/2f+
							E[j][i]*(S[j][i]-S[j][i])*delx1/2f+
							(F[j][i]*S[j][i]+G[j][i])*delx1*delx2;
					
					R[j][i]*=optArg/(A[j][i]*ratio+C[j][i]/ratio)/2f;
					
					S[j][i]+=R[j][i];
				}
				
				if(dim1BC==BoundaryCondition.Periodic){
					// for the east boundary SOR i.e., i==x-1
					if(A[j][x-1]==undef||B[j][x-1]==undef||C[j][x-1]==undef||D[j][x-1]==undef||
					E[j][x-1]==undef||F[j][x-1]==undef||G[j][x-1]==undef) continue;
					
					R[j][x-1]=A[j][x-1]*(S[j][x-1+1]+S[j][x-1-1]-2f*S[j][x-1])*ratio+
							B[j][x-1]*(S[j+1][x-1+1]+S[j-1][x-1-1]-S[j-1][x-1+1]-S[j+1][x-1-1])/4f+
							C[j][x-1]*(S[j+1][x-1]+S[j-1][x-1]-2f*S[j][x-1])+
							D[j][x-1]*(S[j][x-1+1]-S[j][x-1-1])*delx2/2f+
							E[j][x-1]*(S[j][x-1]-S[j][x-1])*delx1/2f+
							(F[j][x-1]*S[j][x-1]+G[j][x-1])*delx1*delx2;
					
					R[j][x-1]*=optArg/(A[j][x-1]*ratio+C[j][x-1]/ratio)/2f;
					
					S[j][x-1]+=R[j][x-1];
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
	 * parameters used for SOR iteration
	 */
	private static final class Params{
		//
		int dim1C=0;	// grids in the first dimension
		int dim2C=0;	// grids in the second dimension
		
		float ratio =0;	// ratio divided by 4
		float delx1 =0;	// increment in the first  dimension
		float delx2 =0;	// increment in the second dimension
		float optArg=0;	// optimal argument for the SOR
		float undef =0;	// undefined value
		
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
		(int dim1C,int dim2C,float delx1,float delx2,float undef,BoundaryCondition dim1BC,BoundaryCondition dim2BC){
			this.dim1C =dim1C;
			this.dim2C =dim2C;
			this.delx1 =delx1;
			this.delx2 =delx2;
			this.undef =undef;
			this.dim1BC=dim1BC;
			this.dim2BC=dim2BC;
			
			ratio=delx2/delx1;
			
			float epsilon=(float)(pow(sin(PI/(2*dim1C+2)),2)+pow(sin(PI/(2*dim2C+2)),2));
			optArg=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		}
	}
}
