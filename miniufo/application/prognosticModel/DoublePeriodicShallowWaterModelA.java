/**
 * @(#)ShallowWaterModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import static java.lang.Math.sqrt;


/**
 * shallow water wave model
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DoublePeriodicShallowWaterModelA extends ShallowWaterModel2D{
	//
	private float deltaR=0;
	
	private float epsilonX=0;
	private float epsilonY=epsilonX;
	
	private float[][] Q=null;	// 2*y*x, diabatic heating
	private float[][] H=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public DoublePeriodicShallowWaterModelA(int dt,int y,int x){
		super(dt,y,x);
		
		Q=new float[y+2][x+2];
		H=new float[y+2][x+2];
		
		deltaR=50000*2/deltaT;
		
		float dtmax=deltaR/350/1.414f;
		
		System.out.println("Bounded shallow water model integration parameters:");
		System.out.println("  delta dt : "+deltaT+" second(s), < "+dtmax+" s");
		
		vd.getModel().xmin=1;
		vd.getModel().xmax=x;
		vd.getModel().ymin=1;
		vd.getModel().ymax=y;
	}
	
	
	public void integration(int step,int interval,String output){
		System.out.println("\nstart integration...");
		
		int tag1=1,tag2=0;
		
		copyToDispBuf(h[0]);
		vd.setData(dispBuf);
		
		for(int l=1;l<=step;l++){
			float[][] hres=h[tag1];	float[][] hsrc=h[tag2];
			float[][] ures=u[tag1];	float[][] usrc=u[tag2];
			float[][] vres=v[tag1];	float[][] vsrc=v[tag2];
			
			// predict h
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				hres[j][i]=hsrc[j][i]-Q[j][i]*deltaT-
				(usrc[j][i+1]*hsrc[j][i+1]-usrc[j][i-1]*hsrc[j][i-1])/deltaR-
				(vsrc[j+1][i]*hsrc[j+1][i]-vsrc[j-1][i]*hsrc[j-1][i])/deltaR+
				(hsrc[j][i+1]+hsrc[j][i-1]-2*hsrc[j][i])/deltaR/deltaR/deltaT*epsilonX+
				(hsrc[j+1][i]+hsrc[j-1][i]-2*hsrc[j][i])/deltaR/deltaR/deltaT*epsilonY;
			}
			
			updateMassBoundaries(hres,tag1);
			
			for(int j=0,J=y+2;j<J;j++)
			for(int i=0,I=x+2;i<I;i++) H[j][i]=(float)sqrt(hres[j][i]);
			
			// u-wind
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				ures[j][i]=usrc[j][i]*H[j][i]-(
					H[j][i]*(hres[j][i+1]-hres[j][i-1])
				)/deltaR-(
					(usrc[j][i+1]*usrc[j][i+1]*H[j][i+1]-usrc[j][i-1]*usrc[j][i-1]*H[j][i-1])+
					usrc[j][i]*(usrc[j][i+1]*H[j][i+1]-usrc[j][i-1]*H[j][i-1])
				)/deltaR/2-(
					(vsrc[j+1][i]*usrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*usrc[j-1][i]*H[j-1][i])+
					vsrc[j][i]*(usrc[j+1][i]*H[j+1][i]-usrc[j-1][i]*H[j-1][i])
				)/deltaR/2+epsilonX*H[j][i]*(
					usrc[j][i+1]+usrc[j][i-1]-2*usrc[j][i]
				)/deltaR/deltaR/deltaT+epsilonY*H[j][i]*(
					usrc[j+1][i]+usrc[j-1][i]-2*usrc[j][i]
				)/deltaR/deltaR/deltaT;
				
				ures[j][i]/=H[j][i];
				
				vres[j][i]=vsrc[j][i]*H[j][i]-(
					H[j][i]*(hres[j+1][i]-hres[j-1][i])
				)/deltaR-(
					(usrc[j][i+1]*vsrc[j][i+1]*H[j][i+1]-usrc[j][i-1]*vsrc[j][i-1]*H[j][i-1])+
					usrc[j][i]*(vsrc[j][i+1]*H[j][i+1]-vsrc[j][i-1]*H[j][i-1])
				)/deltaR/2-(
					(vsrc[j+1][i]*vsrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*vsrc[j-1][i]*H[j-1][i])+
					vsrc[j][i]*(vsrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*H[j-1][i])
				)/deltaR/2+epsilonX*H[j][i]*(
					vsrc[j][i+1]+vsrc[j][i-1]-2*vsrc[j][i]
				)/deltaR/deltaR/deltaT+epsilonY*H[j][i]*(
					vsrc[j+1][i]+vsrc[j-1][i]-2*vsrc[j][i]
				)/deltaR/deltaR/deltaT;
				
				vres[j][i]/=H[j][i];
			}
			
			updateWindBoundaries(ures,vres,tag1);
			
			// output and smooth if necessary
			if(l%interval==0){
				smooth5(h[tag1]); smooth5(u[tag1]); smooth5(v[tag1]);
				
				cTotalMassAndEnergy(tag1);
				
				System.out.println(
					"  output for "+l+" step, total mass and energy is: "+totalMass+"\t"+totalEnergy
				);
				
				copyToDispBuf(h[tag1]);
				vd.setData(dispBuf);
				vd.updateModel(50);
			}
			
			int tmp=tag1;	tag1=tag2;	tag2=tmp;
		}
		
		System.out.println("finish integration.");
	}
	
	
	public void addQDisturbance(int x0,int y0,float rad,float amp){
		addDisturbance(x0,y0,rad,amp,Q);
	}
	
	
	protected void cTotalMassAndEnergy(int tag){
		totalMass=0;	totalEnergy=0;	
		
		for(int j=1;j<=y;j++){
			float tmp0=0,tmp1=0;
			
			for(int i=1;i<=x;i++){
				tmp0+=h[tag][j][i];
				tmp1+=(h[tag][j][i]+u[tag][j][i]*u[tag][j][i]+v[tag][j][i]*v[tag][j][i])*h[tag][j][i]/2;
			}
			
			totalMass+=tmp0;	totalEnergy+=tmp1;
		}
	}
	
	protected void updateMassBoundaries(float[][] hres,int tag){
		periodicExtends(hres);
	}
	
	protected void updateWindBoundaries(float[][] ures,float[][] vres,int tag){
		periodicExtends(ures);
		periodicExtends(vres);
	}
	
	
	private void periodicExtends(float[][] a){
		for(int j=0,J=y+2;j<J;j++){
			a[j][ 0 ]=a[j][x];
			a[j][x+1]=a[j][1];
		}
		
		for(int i=0,I=x+2;i<I;i++){
			a[ 0 ][i]=a[y][i];
			a[y+1][i]=a[1][i];
		}
	}
	
	
	/** test*/
	public static void main(String[] args){
		try{
			DoublePeriodicShallowWaterModelA gswm=new DoublePeriodicShallowWaterModelA(60,192,192);
			
			gswm.initialMirrorSurface();
			gswm.addHDisturbance(80,60,30,3f);
			gswm.addHDisturbance(80,120,35,-3.5f);
			
			//gswm.integration(12*60,12*10,"d:/Gill.dat");
			
			gswm.integration(60*24*3,20,"d:/Gill.dat");
			
		}catch(Exception ex){ ex.printStackTrace();}
	}
}
