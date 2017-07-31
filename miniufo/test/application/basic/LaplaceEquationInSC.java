/**
 * @(#)LaplaceEquation.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;

import miniufo.concurrent.ConcurrentUtil;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.application.EquationInSphericalCoordinate;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;


/**
 * Laplace equation class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class LaplaceEquationInSC extends EquationInSphericalCoordinate{
	//
	protected int threadCount=1;
	protected int mxLoopCount=5000;
	protected float tolerance=0.0001f;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public LaplaceEquationInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/**
     * to solve the Laplace Equation in the form:
     * L(res)=frc
     * using SOR method
     *
     * @param	res	result on the left of the equation
     * @param	frc	a given force on the right of the equation
     *
     * @exception	if r,v are not dimensionally the same
     */
	public void solve(boolean extBC,int loop,Variable res,Variable frc){
		System.out.println("\nStart solving laplace equation...");
		
		if(frc!=null&&!res.isLike(frc))
			throw new IllegalArgumentException("dimensions not same");
		
		if(loop>0) mxLoopCount=loop;
		
		tstart=res.getRange().getTRange()[0];
		zstart=res.getRange().getZRange()[0];
		ystart=res.getRange().getYRange()[0];
		
		t=res.getTCount();	z=res.getZCount();
		y=res.getYCount();	x=res.getXCount();
		
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<Void> cs=new ExecutorCompletionService<>(es);
		
	    if(res.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] rdata=res.getData()[l][k];
				float[][] fdata=(frc==null?null:frc.getData()[l][k]);
				
				String s=sm.getTDef().getSamples()[l+tstart-1]+" "+zdef[k+zstart-1]/100+"hPa loops ";
				
				// solve
				cs.submit(()->solvingTFirst(extBC,s,rdata,fdata),null);
			}
			
	    }else{
	    	for(int k=0;k<z;k++){
	    		float[][][] rdata=res.getData()[k];
				float[][][] fdata=(frc==null?null:frc.getData()[k]);
				
				for(int l=0;l<t;l++){
					String s=sm.getTDef().getSamples()[l+tstart-1]+" "+zdef[k+zstart-1]/100+"hPa loops ";
					final int ll=l;
					
					// solve
					cs.submit(()->solvingNonTFirst(ll,extBC,s,rdata,fdata),null);
				}
			}
		}
	    
	    try{for(int i=0,I=z*t;i<I;i++) cs.take();}
	    catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
	    
	    System.out.println("Finish solving.");
	}
	
	public void solve(Variable res,Variable frc){ solve(false,-1,res,frc);}
	
	public void solve(int loop,Variable res,Variable frc){ solve(false,loop,res,frc);}
	
	
	/*** getor and setor ***/
	public void setMaxLoopCount(int count){
		if(count<1)
			throw new IllegalArgumentException("count should be greater than 1");
		
		mxLoopCount=count;
	}
	
	public void setTolerance(float tol){
		if(tol<=0)
			throw new IllegalArgumentException("tolerance should be positive");
		
		tolerance=tol;
	}
	
	public void setThreadCount(int count){
		if(count<1)
			throw new IllegalArgumentException("Thread count should be positive");
		
		threadCount=count;
	}
	
	
	/*** helper methods ***/
	private void solvingTFirst(boolean extBC,String info,float[][] rdata,float[][] fdata){
		int xx=0;	float RP=0,ratio=dy/dx,dy2=dy*dy,ratio2=ratio*ratio;
		
		float[][] coeff=new float[3][y-2];
		for(int j=1;j<y-1;j++){
			coeff[1][j-1]=(float)cos((ydef[ystart-1+j]+ydef[ystart+j  ])/2);
			coeff[2][j-1]=(float)cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2);
		}
		
		float epsilon=(float)(pow(sin(PI/(2*x+2)),2)+pow(sin(PI/(2*y+2)),2));
		float optimal=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		
		do{
			if(extBC){ // extend the boundary
				for(int i=1;i<x-1;i++){
					rdata[0  ][i]=rdata[1  ][i];
					rdata[y-1][i]=rdata[y-2][i];
				}
				
				for(int j=1;j<y-1;j++){
					rdata[j][0  ]=rdata[j][1  ];
					rdata[j][x-1]=rdata[j][x-2];
				}
				
				rdata[0  ][0  ]=rdata[1  ][1  ];
				rdata[0  ][x-1]=rdata[1  ][x-2];
				rdata[y-1][0  ]=rdata[y-2][1  ];
				rdata[y-1][x-1]=rdata[y-2][x-2];
			}
			
			for(int j=1,J=y-1;j<J;j++)
			for(int i=1,I=x-1;i<I;i++){
				float correct=(
					1f/lcos[ystart-1+j]*(rdata[j][i+1]-rdata[j][i  ])-
					1f/lcos[ystart-1+j]*(rdata[j][i  ]-rdata[j][i-1])
					
				)*ratio2+(
					coeff[1][j-1]*(rdata[j+1][i]-rdata[j  ][i])-
					coeff[2][j-1]*(rdata[j  ][i]-rdata[j-1][i])
					
				)-(fdata==null?0:fdata[j][i]*dy2);
				
				correct*=optimal/(2f/lcos[ystart-1+j]*ratio2+(coeff[1][j-1]+coeff[2][j-1]));
				
				rdata[j][i]+=correct;
				
				correct=correct<0?-correct:correct;
				if(correct>RP) RP=correct;
			}
			
			if(RP<=tolerance||xx++>=mxLoopCount) break;	RP=0;
			
		}while(true);
		
		if(xx<mxLoopCount)
			System.out.println(info+(--xx)+"  absolute err is "+RP);
		else
			System.out.println(info+(--xx)+" Max loop limits, absolute err is "+RP);
	}
	
	private void solvingNonTFirst(int l,boolean extBC,String info,float[][][] rdata,float[][][] fdata){
		int xx=0;	float RP=0,dy2=dy*dy,ratio=dy/dx,ratio2=ratio*ratio;
		
		float[][] coeff=new float[3][y-2];
		for(int j=1;j<y-1;j++){
			coeff[1][j-1]=(float)cos((ydef[ystart-1+j]+ydef[ystart+j  ])/2);
			coeff[2][j-1]=(float)cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2);
		}
		
		float epsilon=(float)(pow(sin(PI/(2*x+2)),2)+pow(sin(PI/(2*y+2)),2));
		float optimal=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		
		do{
			if(extBC){
				for(int i=1;i<x-1;i++){
					rdata[0  ][i][l]=rdata[1  ][i][l];
					rdata[y-1][i][l]=rdata[y-2][i][l];
				}
				
				for(int j=1;j<y-1;j++){
					rdata[j][0  ][l]=rdata[j][1  ][l];
					rdata[j][x-1][l]=rdata[j][x-2][l];
				}
				
				rdata[0  ][0  ][l]=rdata[1  ][1  ][l];
				rdata[0  ][x-1][l]=rdata[1  ][x-2][l];
				rdata[y-1][0  ][l]=rdata[y-2][1  ][l];
				rdata[y-1][x-1][l]=rdata[y-2][x-2][l];
			}
			
			for(int j=1,J=y-1;j<J;j++)
			for(int i=1,I=x-1;i<I;i++){
				float correct=(
					1f/lcos[ystart-1+j]*(rdata[j][i+1][l]-rdata[j][i  ][l])-
					1f/lcos[ystart-1+j]*(rdata[j][i  ][l]-rdata[j][i-1][l])
					
				)*ratio2+(
					coeff[1][j-1]*(rdata[j+1][i][l]-rdata[j  ][i][l])-
					coeff[2][j-1]*(rdata[j  ][i][l]-rdata[j-1][i][l])
					
				)-(fdata==null?0:fdata[j][i][l]*dy2);
				
				correct*=optimal/(2f/lcos[ystart-1+j]*ratio2+(coeff[1][j-1]+coeff[2][j-1]));
				
				rdata[j][i][l]+=correct;
				
				correct=correct<0?-correct:correct;
				if(correct>RP) RP=correct;
			}
			
			if(RP<=tolerance||xx++>=mxLoopCount) break;	RP=0;
			
		}while(true);
		
		if(xx<mxLoopCount)
			System.out.println(info+(--xx)+"  absolute err is "+RP);
		else
			System.out.println(info+(--xx)+" Max loop limits, absolute err is "+RP);
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
	    	
	    }catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
	}*/
}
