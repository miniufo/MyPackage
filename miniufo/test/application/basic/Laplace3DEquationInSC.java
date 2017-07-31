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
public class Laplace3DEquationInSC extends EquationInSphericalCoordinate{
	//
	protected int mxLoopCount=1000;
	protected float tolerance=0.1f;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public Laplace3DEquationInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
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
	public void solve(boolean extBC,Variable res,Variable frc){
		System.out.println("\nStart solving laplace equation...");
		
		if(frc!=null&&!res.isLike(frc))
			throw new IllegalArgumentException("dimensions not same");
		
		tstart=res.getRange().getTRange()[0];
		zstart=res.getRange().getZRange()[0];
		ystart=res.getRange().getYRange()[0];
		
		t=res.getTCount();	z=res.getZCount();
		y=res.getYCount();	x=res.getXCount();
		
        ExecutorService es=ConcurrentUtil.defaultExecutor();
        CompletionService<Void> cs=new ExecutorCompletionService<>(es);
		
	    if(!res.isTFirst()) throw new IllegalArgumentException("T-first required");
	    
		for(int l=0;l<t;l++){
			float[][][] rdata=res.getData()[l];
			float[][][] fdata=(frc==null?null:frc.getData()[l]);
			
			String s=sm.getTDef().getSamples()[l+tstart-1]+" ";
			
			// solve
			cs.submit(()->solveTFirst(s,rdata,fdata),null);
		}
	    
	    try{for(int i=0,I=t;i<I;i++) cs.take();}
	    catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
		
		System.out.println("Finish solving.");
	}
	
	public void solve(Variable res,Variable frc){ solve(false,res,frc);}
	
	
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
	
	
	/*** helper methods ***/
	private void solveTFirst(String info,float[][][] rdata,float[][][] fdata){
		int xx=0;	float RP=0,dy2=dy*dy,dp2=dz*dz;
		
		float[] coeffx=new float[y-2];
		float[] coefy1=new float[y-2];
		float[] coefy2=new float[y-2];
		float[] coeffp=new float[y-2];
		
		for(int j=1;j<y-1;j++){
			coeffx[j-1]=lcos[ystart-1+j]/(dxs[ystart-1+j]*dxs[ystart-1+j]);
			coefy1[j-1]=(float)cos((ydef[ystart-1+j]+ydef[ystart+j  ])/2)/dy2;
			coefy2[j-1]=(float)cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2)/dy2;
			coeffp[j-1]=lcos[ystart-1+j]/dp2;
		}
		
		float epsilon=(float)(pow(sin(PI/(2*x+2)),2)+pow(sin(PI/(2*y+2)),2)+pow(sin(PI/(2*z+2)),2));
		float optimal=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		
		do{
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,J=y-1;j<J;j++)
			for(int i=1,I=x-1;i<I;i++){
				float correct=(
					(rdata[k][j][i+1]-rdata[k][j][i  ])-
					(rdata[k][j][i  ]-rdata[k][j][i-1])
					
				)*coeffx[j-1]+(
					coefy1[j-1]*(rdata[k][j+1][i]-rdata[k][j  ][i])-
					coefy2[j-1]*(rdata[k][j  ][i]-rdata[k][j-1][i])
					
				)+(
					(rdata[k+1][j][i]-rdata[k  ][j][i])-
					(rdata[k  ][j][i]-rdata[k-1][j][i])
					
				)*coeffp[j-1]-(fdata==null?0:fdata[k][j][i]*lcos[ystart-1+j]);
				
				correct*=optimal/(coeffx[j-1]*2+(coefy1[j-1]+coefy2[j-1])+coeffp[j-1]*2);
				
				rdata[k][j][i]+=correct;
				
				correct=correct<0?-correct:correct;
				if(correct>RP) RP=correct;
			}
			
			if(RP<=tolerance||xx++>=mxLoopCount) break;	RP=0;
			if(xx%100==0) System.out.println(xx);
			
		}while(true);
		
		if(xx<mxLoopCount)
			System.out.println(info+(--xx)+"  absolute err is "+RP);
		else
			System.out.println(info+(--xx)+" Max loop limits, absolute err is "+RP);
	}
	
	
	/** test
	public static void main(String[] arg){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		
		DynamicMethodsInSC dm=new DynamicMethodsInSC(ssm);
		Laplace3DEquationInSC le=new Laplace3DEquationInSC(ssm);
		
		Range r=new Range("t(1,1)",dd);
		
		Variable[] vs=df.getVariables(r,true,"u","v","w");
		Variable u=vs[0];
		Variable v=vs[1];
		Variable w=vs[2];
		
		Variable frc=new Variable("frc",u);
		Variable res=new Variable("res",u);
		
		//	frc.add3DDisturbance(180,91,19,2,10);
		for(int l=0,L=frc.getTCount();l<L;l++)
		for(int k=0,K=frc.getTCount();k<K;k++)
		for(int j=0,J=frc.getTCount();j<J;j++)
		for(int i=0,I=frc.getTCount();i<I;i++)
		frc.getData()[l][k][j][i]=2*(float)
		Math.exp(-(Math.sqrt((180-i)*(180-i)+(91-j)*(91-j)+(19-k)*(19-k))/10));
		
		le.solve(res,frc);
		
		Variable[] comp=dm.c3DGradient(res);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/test.dat");
		dw.writeData(dd,u,v,w,frc,res,comp[0],comp[1],comp[2]);	dw.closeFile();
	}*/
}
