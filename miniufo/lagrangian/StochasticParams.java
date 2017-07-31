/**
 * @(#)StochasticParams.java	1.0 2014.08.27
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;


/**
 * all possible parameters for stochastic model
 *
 * @version 1.0, 2014.08.27
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class StochasticParams{
	//
	private static final int dim=2;	// only supporting 2D
	
	private int order=-1;
	
	private float[][] Tvel=null;	// velocity timescale tensor Tvij, for 2D case: [Tvxx Tvxy; Tvyx Tvyy]
	private float[][] Tacc=null;	// velocity timescale tensor Taij, for 2D case: [Taxx Taxy; Tayx Tayy]
	private float[][] Diff=null;	// diffusivity tensor Kij, for 2D case: [Kxx Kxy; Kyx Kyy]
	private float[][] VarV=null;	// velocity variance tensor Vij, for 2D case: [Vxx Vxy; Vyx Vyy]
	private float[][] VarA=null;	// acceleration variance tensor Aij, for 2D case: [Axx Axy; Ayx Ayy]
	private float[][] DGrd=null;	// diffusivity gradient tensor dKjk/dxi, for 2D case: [dKxx/dx dKxy/dy; dKyx/dx dKyy/dy]
	
	
	/**
	 * constructor
	 * 
	 * @param	dt		delta T for tracking (s)
	 * @param	diff	diffusivity tensor (m^2/s)
	 * @param	dGrd	diffusivity gradient tensor (m/s)
	 */
	public StochasticParams(float dt,float[][] diff,float[][] dGrd){
		if(dGrd.length!=dim) throw new IllegalArgumentException("only support 2D applications");
		
		checkDim(diff,dGrd);
		
		order=0;
		
		Diff=diff;
		Tvel=new float[][]{{dt/2f,dt/2f},{dt/2f,dt/2f}};
		Tacc=new float[2][2];
		VarV=new float[2][2];
		VarA=new float[2][2];
		DGrd=dGrd;
		
		for(int i=0;i<2;i++)
		for(int j=0;j<2;j++) VarV[i][j]=Diff[i][j]/Tvel[i][j];
	}
	
	public StochasticParams(float dt,float[][] diff){ this(dt,diff,new float[2][2]);}
	
	/**
	 * constructor
	 * 
	 * @param	tvel	velocity timescale tensor (day)
	 * @param	diff	diffusivity tensor (m^2/s)
	 * @param	dGrd	diffusivity gradient tensor (m/s)
	 */
	public StochasticParams(float[][] tvel,float[][] diff,float[][] dGrd){
		if(dGrd.length!=dim) throw new IllegalArgumentException("only support 2D applications");
		
		checkDim(tvel,diff,dGrd);
		
		order=1;
		
		Diff=diff;
		Tvel=new float[][]{tvel[0].clone(),tvel[1].clone()};
		Tacc=new float[2][2];
		VarV=new float[2][2];
		VarA=new float[2][2];
		DGrd=dGrd;
		
		for(int i=0;i<2;i++)
		for(int j=0;j<2;j++){
			Tvel[i][j]*=24f*60f*60f;	// changing unit to s
			VarV[i][j]=Diff[i][j]/Tvel[i][j];
		}
	}
	
	public StochasticParams(float[][] tvel,float[][] diff){ this(tvel,diff,new float[2][2]);}
	
	/**
	 * constructor
	 * 
	 * @param	tvel	velocity timescale tensor (s)
	 * @param	tacc	velocity timescale tensor (s)
	 * @param	diff	diffusivity tensor (m^2/s)
	 * @param	dGrd	diffusivity gradient tensor (m/s)
	 */
	public StochasticParams(float[][] tvel,float[][] tacc,float[][] diff,float[][] dGrd){
		if(dGrd.length!=dim) throw new IllegalArgumentException("only support 2D applications");
		
		checkDim(tvel,tacc,diff,dGrd);
		
		order=2;
		
		Diff=diff;
		Tvel=new float[][]{tvel[0].clone(),tvel[1].clone()};
		Tacc=new float[][]{tacc[0].clone(),tacc[1].clone()};
		VarV=new float[2][2];
		VarA=new float[2][2];
		DGrd=dGrd;
		
		for(int i=0;i<2;i++)
		for(int j=0;j<2;j++){
			Tvel[i][j]*=24f*60f*60f;	// changing unit to s
			Tacc[i][j]*=24f*60f*60f;	// changing unit to s
			VarV[i][j]=Diff[i][j]/Tvel[i][j];
			VarA[i][j]=VarV[i][j]/Tvel[i][j]/Tacc[i][j];
		}
	}
	
	
	/*** getor and setor ***/
	public int getOrder(){ return order;}
	
	/**
	 * get parameters
	 * 
	 * @param	i	zonal component, 1 or 2
	 * @param	j	meridional component, 1 or 2
	 */
	public float getTL(int i,int j){ checkIndex(i,j); return Tvel[i-1][j-1]+Tacc[i-1][j-1];}
	
	public float getTVel(int i,int j){ checkIndex(i,j); return Tvel[i-1][j-1];}
	
	public float getTAcc(int i,int j){ checkIndex(i,j); return Tacc[i-1][j-1];}
	
	public float getDiff(int i,int j){ checkIndex(i,j); return Diff[i-1][j-1];}
	
	public float getVarV(int i,int j){ checkIndex(i,j); return VarV[i-1][j-1];}
	
	public float getVarA(int i,int j){ checkIndex(i,j); return VarA[i-1][j-1];}
	
	public float getDGrd(int i,int j){ checkIndex(i,j); return DGrd[i-1][j-1];}
	
	
	/**
	 * used to print out
	 */
	public String toString(){
		return String.format(
			"              TL (day)    Tv (day)    Ta (day)   VarV (cm/s)^2   VarA (10^-12 cm/s^2)^2   Diff (m^2/s)\n"+
			"zonal:     %10.4f  %10.4f  %10.4f %13.4f           %14.4f    %13.4f\n"+
			"meridional:%10.4f  %10.4f  %10.4f %13.4f           %14.4f    %13.4f",
			getTL(1,1)/86400f,getTVel(1,1)/86400f,getTAcc(1,1)/86400f,getVarV(1,1)*1e4f,getVarA(1,1)*1e16f,getDiff(1,1),
			getTL(2,2)/86400f,getTVel(2,2)/86400f,getTAcc(2,2)/86400f,getVarV(2,2)*1e4f,getVarA(2,2)*1e16f,getDiff(2,2)
		);
	}
	
	
	/*** helper methods ***/
	private void checkIndex(int... indices){
		for(int idx:indices)
		if(idx<1||idx>2)
		throw new IllegalArgumentException("idx ("+idx+") should be 1 for x or 2 for y");
	}
	
	private void checkDim(float[][]... tensors){
		for(float[][] tensor:tensors)
		if(tensor.length!=dim||tensor[0].length!=dim)
		throw new IllegalArgumentException("invalid dimensions of the given tensor, should be "+dim);
	}
	
	
	/** test
	public static void main(String[] args){
		DenseMatrix64F diff=new DenseMatrix64F(new double[][]{{1000,500},{400,1000}});
		DenseMatrix64F Tvel=new DenseMatrix64F(new double[][]{{10,9},{8,10}});
		DenseMatrix64F Tacc=new DenseMatrix64F(new double[][]{{5,3},{2,5}});
		
		diff.print();
		Tvel.print();
		Tacc.print();
		
		System.out.println(new StochasticParams(86400,diff)+"\n");
		System.out.println(new StochasticParams(Tvel,diff)+"\n");
		System.out.println(new StochasticParams(Tvel,Tacc,diff)+"\n");
	}*/
}
