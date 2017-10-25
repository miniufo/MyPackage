/**
 * @(#)GlobalLaplaceEquation.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;

import miniufo.basic.InterpolationModel;
import miniufo.basic.InterpolationModel.Type;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;


/**
 * global laplace equation class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GlobalLaplaceEquationInSC extends LaplaceEquationInSC{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public GlobalLaplaceEquationInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isPeriodicX())
		throw new IllegalArgumentException("Not a zonal periodic model");
	}
	
	
	/**
     * to solve the Laplace Equation in the form of:
     * L(res)=frc
     * using SOR method
     *
     * @param	extBC	extend boundary condition or not
     * @param	res		result on the left of the equation
     * @param	frc		a given force on the right of the equation
     */
	public void solve(boolean extBC,int loop,Variable res,Variable frc){
		System.out.println("\nStart solving global Laplace equation...");
		
		if(frc!=null&&!res.isLike(frc))
			throw new IllegalArgumentException("dimensions not same");
		
		if(loop>0) mxLoopCount=loop;
		
		tstart=res.getRange().getTRange()[0];
		zstart=res.getRange().getZRange()[0];
		
		t=res.getTCount();	z=res.getZCount();
		y=res.getYCount();	x=res.getXCount();
		
        ExecutorService es=ConcurrentUtil.defaultExecutor();
        CompletionService<ExecutorService> cs=new ExecutorCompletionService<>(es);
		
	    if(res.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] rdata=res.getData()[l][k];
				float[][] fdata=(frc==null?null:frc.getData()[l][k]);
				
				String s=sm.getTDef().getSamples()[l+tstart-1]+" "+zdef[k+zstart-1]/100+"hPa loops ";
				
				// solve
				cs.submit(new SolverTFirst(extBC,s,rdata,fdata),null);
			}
			
		}else{
			for(int k=0;k<z;k++){
				float[][][] rdata=res.getData()[k];
				float[][][] fdata=(frc==null?null:frc.getData()[k]);
				
				for(int l=0;l<t;l++){
					String s=sm.getTDef().getSamples()[l+tstart-1]+" "+zdef[k+zstart-1]/100+"hPa loops ";
					
					// solve
					cs.submit(new SolverNonTFirst(l,extBC,s,rdata,fdata),null);
				}
			}
		}
	    
	    try{for(int i=0,I=z*t;i<I;i++) cs.take();}
	    catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
		
		System.out.println("Finish solving.");
	}
	
	public void solve(boolean extBC,Variable res,Variable frc){
		System.out.println("\nStart solving global Laplace equation...");
		
		if(frc!=null&&!res.isLike(frc))
			throw new IllegalArgumentException("dimensions not same");
		
		tstart=res.getRange().getTRange()[0];
		zstart=res.getRange().getZRange()[0];
		
		t=res.getTCount();	z=res.getZCount();
		y=res.getYCount();	x=res.getXCount();
		
        ExecutorService es=ConcurrentUtil.defaultExecutor();
        CompletionService<ExecutorService> cs=new ExecutorCompletionService<>(es);
		
	    if(res.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] rdata=res.getData()[l][k];
				float[][] fdata=(frc==null?null:frc.getData()[l][k]);
				
				String s=sm.getTDef().getSamples()[l+tstart-1]+" "+zdef[k+zstart-1]/100+"hPa loops ";
				
				// solve
				cs.submit(new SolverTFirst(extBC,s,rdata,fdata),null);
			}
			
		}else{
			for(int k=0;k<z;k++){
				float[][][] rdata=res.getData()[k];
				float[][][] fdata=(frc==null?null:frc.getData()[k]);
				
				for(int l=0;l<t;l++){
					String s=sm.getTDef().getSamples()[l+tstart-1]+" "+zdef[k+zstart-1]/100+"hPa loops ";
					
					// solve
					cs.submit(new SolverNonTFirst(l,extBC,s,rdata,fdata),null);
				}
			}
		}
	    
	    try{for(int i=0,I=z*t;i<I;i++) cs.take();}
	    catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
		
		System.out.println("Finish solving.");
	}
	
	public void solve(Variable res,Variable frc){ solve(false,res,frc);}
	
	
	// Solver for TFirst type
	private final class SolverTFirst implements Runnable{
		//
		private boolean extBC=false;	// default: not extends boundary condition
		
		private float[][] rdata=null;
		private float[][] fdata=null;
		
		private String info=null;
		
		private final int coarseLoop=100;
		private final int[] factors ={5,4,3,2};
		
		
		/**
		 * constructor
		 */
		public SolverTFirst(boolean extBC,String info,float[][] rdata,float[][] fdata){
			this.rdata=rdata;	this.extBC=extBC;
			this.fdata=fdata;	this.info =info;
		}
		
		
		public void run(){
			iterationCoarseGrid();
			iterationFineGrid();
		}
		
		
		private void iterationCoarseGrid(){
			int fctx=generateCoarseGridX();
			int fcty=generateCoarseGridY();
			int loop=coarseLoop/Math.min(fctx,fcty);
			
			int xx=0;				// counter
			int cx=x/fctx;			// x-count in coarse grid
			int cy=(y-1)/fcty+1;	// y-count in coarse grid
			
			float dy2=dy*dy*fcty*fcty;
			
			float[][] ccoe=new float[3][y-2];
			for(int j=fcty,J=y-1;j<J;j+=fcty){
				ccoe[0][j-fcty]=lcos[j]/(dxs[j]*dxs[j]*fctx*fctx);
				ccoe[1][j-fcty]=(float)cos((ydef[j]+ydef[j+fcty])/2);
				ccoe[2][j-fcty]=(float)cos((ydef[j]+ydef[j-fcty])/2);
			}
			
			do{
				for(int j=fcty,J=y-fcty;j<J;j+=fcty){
					/*** East boundary ***/
					float correct=(
						(rdata[j][fctx]-rdata[j][  0   ])-
						(rdata[j][0   ]-rdata[j][x-fctx])
						
					)*ccoe[0][j-fcty]+(
						ccoe[1][j-fcty]*(rdata[j+fcty][0]-rdata[j     ][0])-
						ccoe[2][j-fcty]*(rdata[j     ][0]-rdata[j-fcty][0])
						
					)/dy2-(fdata==null?0:fdata[j][0]);
					
					correct*=1/(ccoe[0][j-fcty]*2+(ccoe[1][j-fcty]+ccoe[2][j-fcty])/dy2);
					
					rdata[j][0]+=correct;
					
					/*** internal area ***/
					for(int i=fctx,I=x-fctx;i<I;i+=fctx){
						correct=(
							(rdata[j][i+fctx]-rdata[j][i     ])-
							(rdata[j][i     ]-rdata[j][i-fctx])
							
						)*ccoe[0][j-fcty]+(
							ccoe[1][j-fcty]*(rdata[j+fctx][i]-rdata[j     ][i])-
							ccoe[2][j-fcty]*(rdata[j     ][i]-rdata[j-fctx][i])
							
						)/dy2-(fdata==null?0:fdata[j][i]);
						
						correct*=1/(ccoe[0][j-fcty]*2+(ccoe[1][j-fcty]+ccoe[2][j-fcty])/dy2);
						
						rdata[j][i]+=correct;
					}
					
					/*** West boundary ***/
					correct=(
						(rdata[j][0     ]-rdata[j][x-fctx  ])-
						(rdata[j][x-fctx]-rdata[j][x-fctx*2])
						
					)*ccoe[0][j-fcty]+(
						ccoe[1][j-fcty]*(rdata[j+fcty][x-fctx]-rdata[j     ][x-fctx])-
						ccoe[2][j-fcty]*(rdata[j     ][x-fctx]-rdata[j-fcty][x-fctx])
						
					)/dy2-(fdata==null?0:fdata[j][x-fctx]);
					
					correct*=1/(ccoe[0][j-fcty]*2+(ccoe[1][j-fcty]+ccoe[2][j-fcty])/dy2);
					
					rdata[j][x-fctx]+=correct;
				}
				
				if(xx++>=loop) break;
				
			}while(true);
			
			// interpolate
			float[][] buf=new float[cy][cx];
			
			for(int j=0;j<cy;j++)
			for(int i=0;i<cx;i++) buf[j][i]=rdata[j*fcty][i*fctx];
			
			InterpolationModel.interp2D(buf,rdata,Type.PERIODIC_CUBIC_P,Type.CUBIC_P);
		}
		
		private void iterationFineGrid(){
			int xx=0;	float RP=0,dy2=dy*dy;
			
			float[][] coeff=new float[3][y-2];
			for(int j=1;j<y-1;j++){
				coeff[0][j-1]=lcos[j]/(dxs[j]*dxs[j]);
				coeff[1][j-1]=(float)cos((ydef[j]+ydef[j+1])/2);
				coeff[2][j-1]=(float)cos((ydef[j]+ydef[j-1])/2);
			}
			
			float epsilon=(float)(pow(sin(PI/(2*x+2)),2)+pow(sin(PI/(2*y+2)),2));
			float optimal=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
			
			do{
				for(int j=1,J=y-1;j<J;j++){
					/*** East boundary ***/
					float corrt=(
						(rdata[j][1]-rdata[j][0  ])-
						(rdata[j][0]-rdata[j][x-1])
						
					)*coeff[0][j-1]+(
						coeff[1][j-1]*(rdata[j+1][0]-rdata[j  ][0])-
						coeff[2][j-1]*(rdata[j  ][0]-rdata[j-1][0])
						
					)/dy2-(fdata==null?0:fdata[j][0]);
					
					corrt*=optimal/(coeff[0][j-1]*2+(coeff[1][j-1]+coeff[2][j-1])/dy2);
					
					rdata[j][0]+=corrt;
					
					corrt=corrt<0?-corrt:corrt;
					if(corrt>RP) RP=corrt;
					
					/*** internal area ***/
					for(int i=1,I=x-1;i<I;i++){
						corrt=(
							(rdata[j][i+1]-rdata[j][i  ])-
							(rdata[j][i  ]-rdata[j][i-1])
							
						)*coeff[0][j-1]+(
							coeff[1][j-1]*(rdata[j+1][i]-rdata[j  ][i])-
							coeff[2][j-1]*(rdata[j  ][i]-rdata[j-1][i])
							
						)/dy2-(fdata==null?0:fdata[j][i]);
						
						corrt*=optimal/(coeff[0][j-1]*2+(coeff[1][j-1]+coeff[2][j-1])/dy2);
						
						rdata[j][i]+=corrt;
						
						corrt=corrt<0?-corrt:corrt;
						if(corrt>RP) RP=corrt;
					}
					
					/*** West boundary ***/
					corrt=(
						(rdata[j][0  ]-rdata[j][x-1])-
						(rdata[j][x-1]-rdata[j][x-2])
						
					)*coeff[0][j-1]+(
						coeff[1][j-1]*(rdata[j+1][x-1]-rdata[j  ][x-1])-
						coeff[2][j-1]*(rdata[j  ][x-1]-rdata[j-1][x-1])
						
					)/dy2-(fdata==null?0:fdata[j][x-1]);
					
					corrt*=optimal/(coeff[0][j-1]*2+(coeff[1][j-1]+coeff[2][j-1])/dy2);
					
					rdata[j][x-1]+=corrt;
					
					corrt=corrt<0?-corrt:corrt;
					if(corrt>RP) RP=corrt;
				}
				
				if(RP<=tolerance||xx++>=mxLoopCount) break;	RP=0;
				
				if(extBC)
				for(int i=0;i<x;i++) rdata[y-1][i]=rdata[y-2][i];
				
			}while(true);
			
			if(xx<mxLoopCount)
				System.out.println(info+(--xx)+"  absolute err is "+RP);
			else
				System.out.println(info+(--xx)+" Max loop limits, absolute err is "+RP);
		}
		
		private int generateCoarseGridX(){
			for(int i:factors) if(x%i==0) return i;
			
			System.out.println(
				"Warning: in Multi-grid method\n"+
				"no suitable factor in x direction"
			);
			
			return 1;
		}
		
		private int generateCoarseGridY(){
			for(int i:factors) if((y-1)%i==0) return i;
			
			System.out.println(
				"Warning: in Multi-grid method\n"+
				"no suitable factor in y direction"
			);
			
			return 1;
		}
	}
	
	// solver for non-TFirst
	private final class SolverNonTFirst implements Runnable{
		//
		private int l;
		
		private boolean extBC=false;
		
		private float[][][] rdata=null;
		private float[][][] fdata=null;
		
		private String info=null;
		
		
		/**
		 * constructor
		 */
		public SolverNonTFirst(int l,boolean extBC,String info,float[][][] rdata,float[][][] fdata){
			this.rdata=rdata;	this.extBC=extBC;
			this.fdata=fdata;	this.info =info;	this.l=l;
		}
		
		
		public void run(){ iterationFineGrid();}
		
		
		private void iterationFineGrid(){
			int xx=0;	float RP=0,dy2=dy*dy;
			
			float[][] coeff=new float[3][y-2];
			for(int j=1;j<y-1;j++){
				coeff[0][j-1]=lcos[j]/(dxs[j]*dxs[j]);
				coeff[1][j-1]=(float)cos((ydef[j]+ydef[j+1])/2);
				coeff[2][j-1]=(float)cos((ydef[j]+ydef[j-1])/2);
			}
			
			float epsilon=(float)(pow(sin(PI/(2*x+2)),2)+pow(sin(PI/(2*y+2)),2));
			float optimal=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
			
			do{
				for(int j=1,J=y-1;j<J;j++){
					/*** East boundary ***/
					float R=(
						(rdata[j][1][l]-rdata[j][0  ][l])-
						(rdata[j][0][l]-rdata[j][x-1][l])
						
					)*coeff[0][j-1]+(
						coeff[1][j-1]*(rdata[j+1][0][l]-rdata[j  ][0][l])-
						coeff[2][j-1]*(rdata[j  ][0][l]-rdata[j-1][0][l])
						
					)/dy2-(fdata==null?0:fdata[j][0][l]);
					
					R*=optimal/(coeff[0][j-1]*2+(coeff[1][j-1]+coeff[2][j-1])/dy2);
					
					rdata[j][0][l]+=R;
					
					R=R<0?-R:R;
					if(R>RP) RP=R;
					
					/*** internal area ***/
					for(int i=1,I=x-1;i<I;i++){
						R=(
							(rdata[j][i+1][l]-rdata[j][i  ][l])-
							(rdata[j][i  ][l]-rdata[j][i-1][l])
							
						)*coeff[0][j-1]+(
							coeff[1][j-1]*(rdata[j+1][i][l]-rdata[j  ][i][l])-
							coeff[2][j-1]*(rdata[j  ][i][l]-rdata[j-1][i][l])
							
						)/dy2-(fdata==null?0:fdata[j][i][l]);
						
						R*=optimal/(coeff[0][j-1]*2+(coeff[1][j-1]+coeff[2][j-1])/dy2);
						
						rdata[j][i][l]+=R;
						
						R=R<0?-R:R;
						if(R>RP) RP=R;
					}
					
					/*** West boundary ***/
					R=(
						(rdata[j][0  ][l]-rdata[j][x-1][l])-
						(rdata[j][x-1][l]-rdata[j][x-2][l])
						
					)*coeff[0][j-1]+(
						coeff[1][j-1]*(rdata[j+1][x-1][l]-rdata[j  ][x-1][l])-
						coeff[2][j-1]*(rdata[j  ][x-1][l]-rdata[j-1][x-1][l])
						
					)/dy2-(fdata==null?0:fdata[j][x-1][l]);
					
					R*=optimal/(coeff[0][j-1]*2+(coeff[1][j-1]+coeff[2][j-1])/dy2);
					
					rdata[j][x-1][l]+=R;
					
					R=R<0?-R:R;
					if(R>RP) RP=R;
				}
				
				if(RP<=tolerance||xx++>=mxLoopCount) break;	RP=0;
				
				if(extBC)
				for(int i=0;i<x;i++) rdata[y-1][i]=rdata[y-2][i];
				
			}while(true);
			
			if(xx<mxLoopCount)
				System.out.println(info+(--xx)+"  absolute err is "+RP);
			else
				System.out.println(info+(--xx)+" Max loop limits, absolute err is "+RP);
		}
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			miniufo.descriptor.CtlDescriptor ctl=new 
			miniufo.descriptor.CtlDescriptor("D:/Data/DiagnosisVortex/Haima/Haima.ctl");
			
			miniufo.diagnosis.Range r=new miniufo.diagnosis.Range("t(1,3);z(35,37)",ctl);
			
			Variable u=new Variable("u",false,r);
			Variable v=new Variable("v",false,r);
			
			miniufo.io.DataRead dr=null;
			dr=miniufo.io.DataIOFactory.getDataRead(ctl);
			dr.readData(u);	dr.closeFile();
			dr=miniufo.io.DataIOFactory.getDataRead(ctl);
			dr.readData(v);	dr.closeFile();
			
			SphericalSpacialModel ssm=new SphericalSpacialModel(ctl);
			GlobalWindFieldInSC gwf=new GlobalWindFieldInSC(ssm);
			
			gwf.setThreadCount(2);
			
			Variable pf=gwf.cVelocityPotentialBySOR(u,v);
			Variable sf=gwf.cStreamFunctionBySOR(u,v);
			
			miniufo.io.DataWrite dw=miniufo.io.DataIOFactory.getDataWrite(ctl,"d:/pf.dat");
			dw.writeData(ctl,u,v,pf,sf);
			
		}catch(Exception e){ e.printStackTrace(); System.exit(0);}
	}*/
}
