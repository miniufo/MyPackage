/**
 * @(#)LSM0th.java	1.0 2014.08.29
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.function.Function;

import miniufo.descriptor.DataDescriptor;
import static java.lang.Math.sqrt;


/**
 * Zeroth-order stochastic model implemented using
 * 4th-order Runge-Kutta method
 *
 * @version 1.0, 2014.08.29
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LSM0th extends StochasticModel{
	
	/**
	 * constructor
	 * 
	 * @param	dtRatio		finite time interval for integration (s)
	 * @param	constMean	whether to use time-invariant mean flow
	 * @param	dd			DataDescriptor describe the domain
	 * @param	BCx			zonal boundary condition
	 * @param	BCy			meridional boundary condition
	 * @param	func		a function that map the current Record to StochasticParams
	 */
	public LSM0th(int dtRatio,boolean constMean,DataDescriptor dd,BCType BCx,BCType BCy,Function<Record,StochasticParams> func){
		super(dtRatio,constMean,dd,BCx,BCy,func);
	}
	
	public LSM0th(int dtRatio,DataDescriptor dd,BCType BCx,BCType BCy,Function<Record,StochasticParams> func){
		super(dtRatio,false,dd,BCx,BCy,func);
	}
	
	
	/*** getor and setor ***/
	public int getOrder(){ return 0;}
	
	
	/**
	 * integrate forward over dt using 4th-order RK method
	 * 
	 * @param	init	initial record
	 * @param	dt		delta time
	 */
	protected Record forwardRK4(Record init,float dt){
		//// compute background mean flow using 4th-order Runge-Kutta method ////
		Record r2=forwardByMean(init,init,dt/2f);	if(r2==null) return null;
		Record r3=forwardByMean(init,r2  ,dt/2f);	if(r3==null) return null;
		Record r4=forwardByMean(init,r3  ,dt   );	if(r4==null) return null;
		
		float[] velm=fetchVelocity(init.getTime(),init.getXPos(),init.getYPos());
		
		float velXk1=velm[0];
		float velYk1=velm[1];
		float velXk2=  r2.getDataValue(0);
		float velYk2=  r2.getDataValue(1);
		float velXk3=  r3.getDataValue(0);
		float velYk3=  r3.getDataValue(1);
		float velXk4=  r4.getDataValue(0);
		float velYk4=  r4.getDataValue(1);
		
		float velXRK4=(velXk1+2*velXk2+2*velXk3+velXk4)/6;
		float velYRK4=(velYk1+2*velYk2+2*velYk3+velYk4)/6;
		
		//// compute estimated position using mean background and previous random velocity ////
		float resX0=init.getDataValue(0)-velXk1;
		float resY0=init.getDataValue(1)-velYk1;
		
		r4.setData(0,velXRK4+resX0);
		r4.setData(1,velYRK4+resY0);
		
		Record esti=forwardByAll(init,r4,resX0,resY0,dt,false);
		
		//// compute final position using mean background and averaged random velocity ////
		float[] res1=getRandom(esti,null);
		
		r4.setData(0,velXRK4+res1[0]);
		r4.setData(1,velYRK4+res1[1]);
		
		r4=forwardByAll(init,r4,res1[0],res1[1],dt,true);
		
		return r4;
	}
	
	protected float[] getRandom(Record init,Record esti){
		StochasticParams sp0=func.apply(init); validateOrder(sp0);
		
		double Kxx=sp0.getDiff(1,1),Kxy=sp0.getDiff(1,2);
		double Kyx=sp0.getDiff(2,1),Kyy=sp0.getDiff(2,2);
		
		double Rx=rnd.nextGaussian(),Ry=rnd.nextGaussian();
		
		float DKxx=sp0.getDGrd(1,1),DKxy=sp0.getDGrd(1,2);
		float DKyx=sp0.getDGrd(2,1),DKyy=sp0.getDGrd(2,2);
		
		// compute inhomogeneous correction
		float velXCr=DKxx+DKxy;
		float velYCr=DKyx+DKyy;
		
		// compute random velocity, = sqrt(2*dt*Kxx)*N(0,1)/dt
		float resX1=(float)(sqrt(2.0/dt*Kxx)*Rx+sqrt(2.0/dt*Kxy)*Ry);
		float resY1=(float)(sqrt(2.0/dt*Kyx)*Rx+sqrt(2.0/dt*Kyy)*Ry);
		
		return new float[]{resX1+velXCr,resY1+velYCr,0,0};
	}
	
	protected float[] spinupRandom(Record init,int iter){ return getRandom(init,null);}
}
