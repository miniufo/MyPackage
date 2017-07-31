/**
 * @(#)GlobalLinearHeatInducedModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import miniufo.application.EquationInSphericalCoordinate;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.mathsphysics.DiscreteFourier;


/**
 * Global Linear Heat-Induced Model
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GlobalLinearHeatInducedModel extends EquationInSphericalCoordinate{
	//
	private  int  x;				// x-grid count
	private  int  y;				// y-grid count
	private  int  b;
	private  int  dt=60;
	private  int  maxloop=3600/dt*24*12;// max loop of SOR
	private float epsilon=8e-6f;	// coefficient of dissipation
	
	private float[][] buff=null;	// buffer for smooth
	
	private Variable U=null;		// the basic state of u-wind
	private Variable V=null;		// the basic state of v-wind
	private Variable H=null;		// the basic state of geopotential
	private Variable Q=null;		// the heating source
	
	private Variable h=null;		// the disturbance of u-wind
	private Variable u=null;		// the disturbance of v-wind
	private Variable v=null;		// the disturbance of geopotential
	
	// fourier filter
	DiscreteFourier df=null;
	
	float[][] cosw=null;
	float[][] sinw=null;
	float[]   bufr=null;
	float[]    Ks =null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public GlobalLinearHeatInducedModel(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isZonalPeriodic())
		throw new IllegalArgumentException("Not a zonal periodic model");
		
		x=ssm.getXCount();	y=ssm.getYCount();
		
		df=new DiscreteFourier(x);
		cosw=df.getCosOmega();
		sinw=df.getSinOmega();
		
		buff=new float[y][x];
		
		bufr=new float[x];	b=6;
		
		Ks=new float[b];
		for(int i=0;i<b;i++) Ks[i]=i;
	}
	
	
	/*** getor and setor **/
	public void setBackgroundU(Variable U){ this.U=U;}
	
	public void setBackgroundV(Variable V){ this.V=V;}
	
	public void setBackgroundH(Variable H){ this.H=H;}
	
	public void setDisturbanceU(Variable u){ this.u=u;}
	
	public void setDisturbanceV(Variable v){ this.v=v;}
	
	public void setDisturbanceH(Variable h){ this.h=h;}
	
	public void setHeatingSource(Variable Q){ this.Q=Q;}
	
	public void setDissipationCoefficient(float e){ this.epsilon=e;}
	
	public Variable getBackgroundU(){ return U;}
	
	public Variable getBackgroundV(){ return V;}
	
	public Variable getBackgroundH(){ return H;}
	
	public Variable getBackgroundQ(){ return Q;}
	
	public Variable getDisturbanceU(){ return u;}
	
	public Variable getDisturbanceV(){ return v;}
	
	public Variable getDisturbanceH(){ return h;}
	
	
	private void smoothPolar(float[][] a){
		int length=a[0].length;
		
		for(int j=1;j<b;j++){
			float[] buf=a[j].clone();
			float ave=miniufo.statistics.StatisticsUtil.cArithmeticMean(buf);
			for(int i=0;i<length;i++) buf[i]-=ave;
			
			float[][] A=df.cFourierAB(buf);
			
			for(int i=0;i<length;i++){
				bufr[i]=0;	// initialized to 0
				
				for(int k=0;k<=Ks[j];k++)
				if(k!=0){
					bufr[i]+=A[0][k-1]*cosw[k-1][i]+A[1][k-1]*sinw[k-1][i];
					
				}else bufr[i]+=ave;
			}
			
			System.arraycopy(bufr,0,a[j],0,length);
		}
		
		for(int j=y-2;j>y-b-1;j--){
			float[] buf=a[j].clone();
			float ave=miniufo.statistics.StatisticsUtil.cArithmeticMean(buf);
			for(int i=0;i<length;i++) buf[i]-=ave;
			
			float[][] A=df.cFourierAB(buf);
			
			for(int i=0;i<length;i++){
				bufr[i]=0;	// initialized to 0
				
				for(int k=0;k<=Ks[y-1-j];k++)
				if(k!=0){
					bufr[i]+=A[0][k-1]*cosw[k-1][i]+A[1][k-1]*sinw[k-1][i];
					
				}else bufr[i]+=ave;
			}
			
			System.arraycopy(bufr,0,a[j],0,length);
		}
	}
	
	private void smooth(float[][] a){
		int length=a[0].length;
		
		for(int j=1;j<y-1;j++){
			for(int i=1;i<x-1;i++) buff[j][i]=a[j][i]/2+(a[j-1][i]+a[j+1][i]+a[j][i-1]+a[j][i+1])/8;
			
			buff[j][ 0 ]=a[j][ 0 ]/2+(a[j+1][0]+a[j-1][0]+a[j][x-1]+a[j][1])/8;
			buff[j][x-1]=a[j][x-1]/2+(a[j+1][x-1]+a[j-1][x-1]+a[j][x-2]+a[j][0])/8;
			
			System.arraycopy(buff[j],0,a[j],0,length);
		}
	}
	
	public void iteration(int interval,String output){
		if(!U.isTFirst()) throw new IllegalArgumentException("U should be T-first type");
		if(!V.isTFirst()) throw new IllegalArgumentException("V should be T-first type");
		if(!H.isTFirst()) throw new IllegalArgumentException("H should be T-first type");
		if(!Q.isTFirst()) throw new IllegalArgumentException("Q should be T-first type");
		if(!u.isTFirst()) throw new IllegalArgumentException("u should be T-first type");
		if(!v.isTFirst()) throw new IllegalArgumentException("v should be T-first type");
		if(!h.isTFirst()) throw new IllegalArgumentException("h should be T-first type");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(output);
		cdws.writeData(h);	cdws.writeData(u);	cdws.writeData(v);
		
		float[][] Qdata=Q.getData()[0][0];
		float[][] Udata=U.getData()[0][0];
		float[][] Vdata=V.getData()[0][0];
		float[][] Hdata=H.getData()[0][0];
		
		float[][][] udata=new float[2][y][x];
		float[][][] vdata=new float[2][y][x];
		float[][][] hdata=new float[2][y][x];
		
		int tag1=1,tag2=0;
		for(int l=0;l<=maxloop;l++){
			// predict h
			for(int j=1,J=y-1;j<J;j++){
				for(int i=1,I=x-1;i<I;i++){
					hdata[tag1][j][i]=hdata[tag2][j][i]+dt*(
						-Hdata[j][i]*(
							(udata[tag2][j][i+1]-udata[tag2][j][i-1])/(dxs[j]+dxs[j])+
							(vdata[tag2][j+1][i]-vdata[tag2][j-1][i])/(dy+dy)-
							vdata[tag2][j][i]*ltan[j]/EARTH_RADIUS
						)
						/**/
						-hdata[tag2][j][i]*(
							(Udata[j][i+1]-Udata[j][i-1])/(dxs[j]+dxs[j])+
							(Vdata[j+1][i]-Vdata[j-1][i])/(dy+dy)-
							Vdata[j][i]*ltan[j]/EARTH_RADIUS
						)
						-(
							udata[tag2][j][i]*(Hdata[j][i+1]-Hdata[j][i-1])+
							Udata[j][i]*(hdata[tag2][j][i+1]-hdata[tag2][j][i-1])
						)/(dxs[j]+dxs[j])
						-(
							vdata[tag2][j][i]*(Hdata[j+1][i]-Hdata[j-1][i])+
							Vdata[j][i]*(hdata[tag2][j+1][i]-hdata[tag2][j-1][i])
						)/(dy+dy)
						-epsilon*hdata[tag2][j][i]
						-Qdata[j][i]
					);
				}
				
				// left
				hdata[tag1][j][0]=hdata[tag2][j][0]+dt*(
					-Hdata[j][0]*(
						(udata[tag2][j][1]-udata[tag2][j][x-1])/(dxs[j]+dxs[j])+
						(vdata[tag2][j+1][0]-vdata[tag2][j-1][0])/(dy+dy)-
						vdata[tag2][j][0]*ltan[j]/EARTH_RADIUS
					)
					/**/
					-hdata[tag2][j][0]*(
						(Udata[j][1]-Udata[j][x-1])/(dxs[j]+dxs[j])+
						(Vdata[j+1][0]-Vdata[j-1][0])/(dy+dy)-
						Vdata[j][0]*ltan[j]/EARTH_RADIUS
					)
					-(
						udata[tag2][j][0]*(Hdata[j][1]-Hdata[j][x-1])+
						Udata[j][0]*(hdata[tag2][j][1]-hdata[tag2][j][x-1])
					)/(dxs[j]+dxs[j])
					-(
						vdata[tag2][j][0]*(Hdata[j+1][0]-Hdata[j-1][0])+
						Vdata[j][0]*(hdata[tag2][j+1][0]-hdata[tag2][j-1][0])
					)/(dy+dy)
					-epsilon*hdata[tag2][j][0]
					-Qdata[j][0]
				);
				
				// right
				hdata[tag1][j][x-1]=hdata[tag2][j][x-1]+dt*(
					-Hdata[j][x-1]*(
						(udata[tag2][j][0]-udata[tag2][j][x-2])/(dxs[j]+dxs[j])+
						(vdata[tag2][j+1][x-1]-vdata[tag2][j-1][x-1])/(dy+dy)-
						vdata[tag2][j][x-1]*ltan[j]/EARTH_RADIUS
					)
					/**/
					-hdata[tag2][j][x-1]*(
						(Udata[j][0]-Udata[j][x-2])/(dxs[j]+dxs[j])+
						(Vdata[j+1][x-1]-Vdata[j-1][x-1])/(dy+dy)-
						Vdata[j][x-1]*ltan[j]/EARTH_RADIUS
					)
					-(
						udata[tag2][j][x-1]*(Hdata[j][0]-Hdata[j][x-2])/+
						Udata[j][x-1]*(hdata[tag2][j][0]-hdata[tag2][j][x-2])
					)/(dxs[j]+dxs[j])
					-(
						vdata[tag2][j][x-1]*(Hdata[j+1][x-1]-Hdata[j-1][x-1])+
						Vdata[j][x-1]*(hdata[tag2][j+1][x-1]-hdata[tag2][j-1][x-1])
					)/(dy+dy)
					-epsilon*hdata[tag2][j][x-1]
					-Qdata[j][x-1]
				);
			}
			
			// u-wind
			for(int j=1,J=y-1;j<J;j++){
				for(int i=1,I=x-1;i<I;i++){
					udata[tag1][j][i]=udata[tag2][j][i]+dt*(
						+f1[j]*vdata[tag2][j][i]
						-(hdata[tag1][j][i+1]-hdata[tag1][j][i-1])/(dxs[j]+dxs[j])
						/**/
						-(Udata[j][i]*(udata[tag2][j][i+1]-udata[tag2][j][i-1]))/(dxs[j]+dxs[j])
						-(Vdata[j][i]*(udata[tag2][j+1][i]-udata[tag2][j-1][i]))/(dy+dy)
						+(
							Udata[j][i]*vdata[tag2][j][i]+Vdata[j][i]*udata[tag2][j][i]
						)*ltan[j]/EARTH_RADIUS
						-epsilon*udata[tag2][j][i]
					);
				}
				
				// left
				udata[tag1][j][0]=udata[tag2][j][0]+dt*(
					+f1[j]*vdata[tag2][j][0]
					-(hdata[tag1][j][1]-hdata[tag1][j][x-1])/(dxs[j]+dxs[j])
					/**/
					-(Udata[j][0]*(udata[tag2][j][1]-udata[tag2][j][x-1]))/(dxs[j]+dxs[j])
					-(Vdata[j][0]*(udata[tag2][j+1][0]-udata[tag2][j-1][0]))/(dy+dy)
					+(
						Udata[j][0]*vdata[tag2][j][0]+Vdata[j][0]*udata[tag2][j][0]
					)*ltan[j]/EARTH_RADIUS
					-epsilon*udata[tag2][j][0]
				);
				
				// right
				udata[tag1][j][x-1]=udata[tag2][j][x-1]+dt*(
					+f1[j]*vdata[tag2][j][x-1]
					-(hdata[tag1][j][0]-hdata[tag1][j][x-2])/(dxs[j]+dxs[j])
					/**/
					-(Udata[j][x-1]*(udata[tag2][j][0]-udata[tag2][j][x-2]))/(dxs[j]+dxs[j])
					-(Vdata[j][x-1]*(udata[tag2][j+1][x-1]-udata[tag2][j-1][x-1]))/(dy+dy)
					+(
						Udata[j][x-1]*vdata[tag2][j][x-1]+Vdata[j][x-1]*udata[tag2][j][x-1]
					)*ltan[j]/EARTH_RADIUS
					-epsilon*udata[tag2][j][x-1]
				);
			}
			
			// v-wind
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					vdata[tag1][j][i]=vdata[tag2][j][i]+dt*(
						-f1[j]*udata[tag2][j][i]
						-(hdata[tag1][j+1][i]-hdata[tag1][j-1][i])/(dy+dy)
						/**/
						-(Udata[j][i]*(vdata[tag2][j][i+1]-vdata[tag2][j][i-1]))/(dxs[j]+dxs[j])
						-(Vdata[j][i]*(vdata[tag2][j+1][i]-vdata[tag2][j-1][i]))/(dy+dy)
						-Udata[j][i]*ltan[j]/EARTH_RADIUS*udata[tag2][j][i]*2
						-epsilon*vdata[tag2][j][i]
					);
				}
				
				// left
				vdata[tag1][j][0]=vdata[tag2][j][0]+dt*(
					-f1[j]*udata[tag2][j][0]
					-(hdata[tag1][j+1][0]-hdata[tag1][j-1][0])/(dy+dy)
					/**/
					-(Udata[j][0]*(vdata[tag2][j][1]-vdata[tag2][j][x-1]))/(dxs[j]+dxs[j])
					-(Vdata[j][0]*(vdata[tag2][j+1][0]-vdata[tag2][j-1][0]))/(dy+dy)
					-Udata[j][0]*ltan[j]/EARTH_RADIUS*udata[tag2][j][0]*2
					-epsilon*vdata[tag2][j][0]
				);
				
				// right
				vdata[tag1][j][x-1]=vdata[tag2][j][x-1]+dt*(
					-f1[j]*udata[tag2][j][x-1]
					-(hdata[tag1][j+1][x-1]-hdata[tag1][j-1][x-1])/(dy+dy)
					/**/
					-(Udata[j][x-1]*(vdata[tag2][j][0]-vdata[tag2][j][x-2]))/(dxs[j]+dxs[j])
					-(Vdata[j][x-1]*(vdata[tag2][j+1][x-1]-vdata[tag2][j-1][x-1]))/(dy+dy)
					-Udata[j][x-1]*ltan[j]/EARTH_RADIUS*udata[tag2][j][x-1]*2
					-epsilon*vdata[tag2][j][x-1]
				);
			}
			
			if(l%5 ==0){ smoothPolar(hdata[tag1]); smoothPolar(udata[tag1]); smoothPolar(vdata[tag1]);}
			if(l%50==0){ smooth(hdata[tag1]); smooth(udata[tag1]); smooth(vdata[tag1]);}
			
			// output
			if(l%interval==0){
				h.getData()[0][0]=hdata[tag1];	cdws.writeData(h);
				u.getData()[0][0]=udata[tag1];	cdws.writeData(u);
				v.getData()[0][0]=vdata[tag1];	cdws.writeData(v);
				System.out.print(".");
			}
			
			tag1=(++tag1)%2;	tag2=(++tag2)%2;
		};
		
		System.out.println();
		
		cdws.closeFile();
	}
	
	
	/** test*/
	public static void main(String[] args){
		try{
			DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/GillModel/MeanFlow/Global.ctl");
			DataDescriptor dd=df.getDataDescriptor();
			SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
			GlobalLinearHeatInducedModel ghi=new GlobalLinearHeatInducedModel(ssm);
			
			Range r=new Range("",dd);
			
			Variable[] mean=df.getVariables(r,true,"hm","um","vm","olr");
			Variable H=mean[0];		H.multiplyEq(9.8f);
			Variable U=mean[1];
			Variable V=mean[2];
			Variable Q=mean[3];		Q.multiplyEq(-0.002f);
			
			for(int j=0;j<Q.getYCount();j++){
				if(j<20) for(int i=0;i<Q.getXCount();i++) Q.getData()[0][0][j][i]=0;
				if(j>52) for(int i=0;i<Q.getXCount();i++) Q.getData()[0][0][j][i]=0;
			}
			
			Variable h=new Variable("h"  ,true,r);
			Variable u=new Variable("u"  ,true,r);
			Variable v=new Variable("v"  ,true,r);
			
			ghi.setHeatingSource(Q);
			ghi.setBackgroundH(H);
			ghi.setBackgroundU(U);
			ghi.setBackgroundV(V);
			ghi.setDisturbanceH(h);
			ghi.setDisturbanceU(u);
			ghi.setDisturbanceV(v);
			
			ghi.iteration(3600/ghi.dt*6,"d:/Data/GillModel/MeanFlow/test.dat");
			
		}catch(Exception ex){ ex.printStackTrace();}
	}
}
