/**
 * @(#)LSM1st.java	1.0 2014.08.27
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.function.Function;
import miniufo.descriptor.DataDescriptor;
import static java.lang.Math.sqrt;


/**
 * First-order stochastic model implemented according to the
 * reference: Oh et al. 2000, JGR
 *
 * @version 1.0, 2014.08.27
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LSM1st extends StochasticModel{
	
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
	public LSM1st(int dtRatio,boolean constMean,DataDescriptor dd,BCType BCx,BCType BCy,Function<Record,StochasticParams> func){
		super(dtRatio,constMean,dd,BCx,BCy,func);
	}
	
	public LSM1st(int dtRatio,DataDescriptor dd,BCType BCx,BCType BCy,Function<Record,StochasticParams> func){
		super(dtRatio,false,dd,BCx,BCy,func);
	}
	
	
	/*** getor and setor ***/
	public int getOrder(){ return 1;}
	
	
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
		float velXk2=r2.getData(VelX);
		float velYk2=r2.getData(VelY);
		float velXk3=r3.getData(VelX);
		float velYk3=r3.getData(VelY);
		float velXk4=r4.getData(VelX);
		float velYk4=r4.getData(VelY);
		
		float velXRK4=(velXk1+2*velXk2+2*velXk3+velXk4)/6;
		float velYRK4=(velYk1+2*velYk2+2*velYk3+velYk4)/6;
		
		//// compute estimated position using mean background and previous random velocity ////
		float resX0=init.getData(VelX)-velXk1;
		float resY0=init.getData(VelY)-velYk1;
		
		r4.setData(VelX,velXRK4+resX0);
		r4.setData(VelY,velYRK4+resY0);
		
		Record esti=forwardByAll(init,r4,resX0,resY0,dt,false);	if(esti==null) return null;
		
		//// compute final position using mean background and averaged random velocity ////
		float[] res1=getRandom(init,esti);
		
		r4.setData(VelX,velXRK4+res1[0]);
		r4.setData(VelY,velYRK4+res1[1]);
		
		r4=forwardByAll(init,r4,res1[0],res1[1],dt,true);
		
		return r4;
	}
	
	protected float[] getRandom(Record init,Record esti){
		StochasticParams sp0=func.apply(init);	validateOrder(sp0);
		StochasticParams sp1=func.apply(esti);	validateOrder(sp1);
		
		float Kxx0 =sp0.getDiff(1,1),Kyy0 =sp0.getDiff(2,2);
		float DKxx0=sp0.getDGrd(1,1),DKyy0=sp0.getDGrd(2,2);
		float Txx0 =sp0.getTVel(1,1),Tyy0 =sp0.getTVel(2,2);
		float vXX0 =sp0.getVarV(1,1),vYY0 =sp0.getVarV(2,2);
		float vXX1 =sp1.getVarV(1,1),vYY1 =sp1.getVarV(2,2);
		
		float[] velm=fetchVelocity(init.getTime(),init.getXPos(),init.getYPos());
		
		float resX0=init.getData(VelX)-velm[0];
		float resY0=init.getData(VelY)-velm[1];
		
		//// compute inhomogeneous correction ////
		float velXCr=(
			//(DKxx/Txx0+DKxy/Txy0)*dt+
			DKxx0/Txx0*dt+
			//resX0/varXX0*(varXX1-varXX0)+resY0/varXY0*(varXX1-varXX0)+
			//resX0/varYX0*(varXY1-varXY0)+resY0/varYY0*(varXY1-varXY0)
			resX0/vXX0*(vXX1-vXX0)
		)/2f;
		
		float velYCr=(
			//(DKyx/Tyx0+DKyy/Tyy0)*dt+
			DKyy0/Tyy0*dt+
			//resX0/varXX0*(varYX1-varYX0)+resY0/varXY0*(varYX1-varYX0)+
			//resX0/varYX0*(varYY1-varYY0)+resY0/varYY0*(varYY1-varYY0)
			resY0/vYY0*(vYY1-vYY0)
		)/2f;
		
		float dvelX=(float)(sqrt(2.0*dt*Kxx0)/Txx0*rnd.nextGaussian());
		float dvelY=(float)(sqrt(2.0*dt*Kyy0)/Tyy0*rnd.nextGaussian());
		
		float resX1=resX0*(1f-dt/Txx0)+velXCr+dvelX;
		float resY1=resY0*(1f-dt/Tyy0)+velYCr+dvelY;
		
		return new float[]{resX1,resY1,0,0};
	}
	
	protected float[] spinupRandom(Record init,int iter){
		StochasticParams sp=func.apply(init); validateOrder(sp);
		
		float Kxx =sp.getDiff(1,1),Kyy =sp.getDiff(2,2);
		float DKxx=sp.getDGrd(1,1),DKyy=sp.getDGrd(2,2);
		float Txx =sp.getTVel(1,1),Tyy =sp.getTVel(2,2);
		
		double coX=sqrt(2.0*dt*Kxx)/Txx;
		double coY=sqrt(2.0*dt*Kyy)/Tyy;
		double vCX=DKxx/Txx/2.0*dt;
		double vCY=DKyy/Tyy/2.0*dt;
		
		float resX0=0f;
		float resY0=0f;
		float resX1=(float)(resX0*(1.0-dt/Txx)+vCX+coX*rnd.nextGaussian());
		float resY1=(float)(resY0*(1.0-dt/Tyy)+vCY+coY*rnd.nextGaussian());
		
		for(int i=1;i<iter;i++){
			resX0=resX1;	resY0=resY1;
			
			resX1=(float)(resX0*(1.0-dt/Txx)+vCX+coX*rnd.nextGaussian());
			resY1=(float)(resY0*(1.0-dt/Tyy)+vCY+coY*rnd.nextGaussian());
		}
		
		return new float[]{resX1,resY1,0,0};
	}
}
