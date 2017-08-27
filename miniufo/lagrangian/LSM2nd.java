/**
 * @(#)LSM2nd.java	1.0 2014.08.27
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.function.Function;
import miniufo.descriptor.DataDescriptor;
import static java.lang.Math.sqrt;


/**
 * Second-order stochastic model implemented according to the
 * reference: Rupolo 2007, JPO
 *
 * @version 1.0, 2014.08.27
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LSM2nd extends StochasticModel{
	
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
	public LSM2nd(int dtRatio,boolean constMean,DataDescriptor dd,BCType BCx,BCType BCy,Function<Record,StochasticParams> func){
		super(dtRatio,constMean,dd,BCx,BCy,func);
	}
	
	public LSM2nd(int dtRatio,DataDescriptor dd,BCType BCx,BCType BCy,Function<Record,StochasticParams> func){
		super(dtRatio,false,dd,BCx,BCy,func);
	}
	
	
	/*** getor and setor ***/
	public int getOrder(){ return 2;}
	
	
	/**
	 * integrate forward over dt using 4th-order RK method
	 * 
	 * @param	init	initial record
	 * @param	dt		delta time
	 */
	protected Record forwardRK4(Record init,float dt){
		Record r2=forwardByMean(init,init,dt/2f);	if(r2==null) return null;
		Record r3=forwardByMean(init,r2  ,dt/2f);	if(r3==null) return null;
		Record r4=forwardByMean(init,r3  ,dt   );	if(r4==null) return null;
		
		float[] velm=fetchVelocity(init.getTime(),init.getXPos(),init.getYPos());
		
		float velXk1=velm[0];
		float velYk1=velm[1];
		float velXk2=r2.getDataValue(0);
		float velYk2=r2.getDataValue(1);
		float velXk3=r3.getDataValue(0);
		float velYk3=r3.getDataValue(1);
		float velXk4=r4.getDataValue(0);
		float velYk4=r4.getDataValue(1);
		
		float velXRK4=(velXk1+2*velXk2+2*velXk3+velXk4)/6;
		float velYRK4=(velYk1+2*velYk2+2*velYk3+velYk4)/6;
		
		//// compute estimated position using mean background and previous random velocity ////
		float resX0=init.getDataValue(0)-velXk1;
		float resY0=init.getDataValue(1)-velYk1;
		
		r4.setData(0,velXRK4+resX0);
		r4.setData(1,velYRK4+resY0);
		
		Record esti=forwardByAll(init,r4,resX0,resY0,dt,false);	if(esti==null) return null;
		
		//// compute final position using mean background and averaged random velocity ////
		float[] res1=getRandom(init,esti);
		
		r4.setData(0,velXRK4+res1[0]); r4.setData(2,res1[2]);
		r4.setData(1,velYRK4+res1[1]); r4.setData(3,res1[3]);
		
		r4=forwardByAll(init,r4,res1[0],res1[1],dt,true);
		
		return r4;
	}
	
	protected float[] getRandom(Record init,Record esti){
		StochasticParams sp0=func.apply(init);	validateOrder(sp0);
		StochasticParams sp1=func.apply(esti);	validateOrder(sp1);
		
		float Kxx0 =sp0.getDiff(1,1),Kyy0 =sp0.getDiff(2,2);
		float DKxx0=sp0.getDGrd(1,1),DKyy0=sp0.getDGrd(2,2);
		float Taxx0=sp0.getTAcc(1,1),Tayy0=sp0.getTAcc(2,2);
		float Tvxx0=sp0.getTVel(1,1),Tvyy0=sp0.getTVel(2,2);
		float vXX0 =sp0.getVarV(1,1),vYY0 =sp0.getVarV(2,2);
		float vXX1 =sp1.getVarV(1,1),vYY1 =sp1.getVarV(2,2);
		
		float[] velm=fetchVelocity(init.getTime(),init.getXPos(),init.getYPos());
		
		float resX0=init.getDataValue(0)-velm[0];
		float resY0=init.getDataValue(1)-velm[1];
		float accX0=init.getDataValue(2);
		float accY0=init.getDataValue(3);
		
		//// compute inhomogeneous correction ////
		float velXCr=(
			//(DKxx/Txx0+DKxy/Txy0)*dt+
			DKxx0/Tvxx0*dt+
			//resX0/varXX0*(varXX1-varXX0)+resY0/varXY0*(varXX1-varXX0)+
			//resX0/varYX0*(varXY1-varXY0)+resY0/varYY0*(varXY1-varXY0)
			resX0/vXX0*(vXX1-vXX0)
		)/2f;
		
		float velYCr=(
			//(DKyx/Tyx0+DKyy/Tyy0)*dt+
			DKyy0/Tvyy0*dt+
			//resX0/varXX0*(varYX1-varYX0)+resY0/varXY0*(varYX1-varYX0)+
			//resX0/varYX0*(varYY1-varYY0)+resY0/varYY0*(varYY1-varYY0)
			resY0/vYY0*(vYY1-vYY0)
		)/2f;
		
		float accXCr=-resX0/Tvxx0/Taxx0*dt;
		float accYCr=-resY0/Tvyy0/Tayy0*dt;
		
		velXCr=velYCr=accXCr=accYCr=0;
		
		float daccX=(float)(sqrt(2.0*dt*Kxx0)/Taxx0/Tvxx0*rnd.nextGaussian());
		float daccY=(float)(sqrt(2.0*dt*Kyy0)/Tayy0/Tvyy0*rnd.nextGaussian());
		
		float accX1=accX0*(1-dt/Taxx0)+accXCr+daccX;
		float accY1=accY0*(1-dt/Tayy0)+accYCr+daccY;
		float resX1=resX0+((accX0+accX1)/2f-resX0/Tvxx0)*dt+velXCr;
		float resY1=resY0+((accY0+accY1)/2f-resY0/Tvyy0)*dt+velYCr;
		
		return new float[]{resX1,resY1,accX1,accY1};
	}
	
	protected float[] spinupRandom(Record init,int iter){
		StochasticParams sp=func.apply(init); validateOrder(sp);
		
		float Kxx =sp.getDiff(1,1),Kyy =sp.getDiff(2,2);
		float DKxx=sp.getDGrd(1,1),DKyy=sp.getDGrd(2,2);
		float Taxx=sp.getTAcc(1,1),Tayy=sp.getTAcc(2,2);
		float Tvxx=sp.getTVel(1,1),Tvyy=sp.getTVel(2,2);
		
		double coX=sqrt(2.0*dt*Kxx)/Tvxx/Taxx;
		double coY=sqrt(2.0*dt*Kyy)/Tvyy/Tayy;
		double vCX=DKxx/Tvxx/2.0*dt;
		double vCY=DKyy/Tvyy/2.0*dt;
		
		float accX0=0;
		float accY0=0;
		float resX0=0;
		float resY0=0;
		
		float accX1=(float)(accX0-accX0/Taxx*dt+coX*rnd.nextGaussian());
		float accY1=(float)(accY0-accY0/Tayy*dt+coY*rnd.nextGaussian());
		float resX1=(float)(resX0+((accX0+accX1)/2f-resX0/Tvxx)*dt+vCX);
		float resY1=(float)(resY0+((accY0+accY1)/2f-resY0/Tvyy)*dt+vCY);
		
		for(int i=1;i<iter;i++){
			accX0=accX1;	accY0=accY1;
			resX0=resX1;	resY0=resY1;
			
			accX1=(float)(accX0-accX0/Taxx*dt+coX*rnd.nextGaussian());
			accY1=(float)(accY0-accY0/Tayy*dt+coY*rnd.nextGaussian());
			resX1=(float)(resX0+((accX0+accX1)/2f-resX0/Tvxx)*dt+vCX);
			resY1=(float)(resY0+((accY0+accY1)/2f-resY0/Tvyy)*dt+vCY);
		}
		
		return new float[]{resX1,resY1,accX1,accY1};
	}
}
