/**
 * @(#)ShallowWaterModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import miniufo.visualize.VisualizeData;
import static java.lang.Math.sqrt;


/**
 * bounded shallow water wave model
 * rigid boundaries in north and south, periodic boundaries in east and west
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class BoundedShallowWaterModelC{
	//
	private int y;
	private int x;
	
	private float deltaT=0;
	private float deltaR=0;
	
	private float[] M_E=null;	// [0] means mass and [1] means energy
	
	private float[][] H   =null;
	
	private float[][][] h=null;	// 2*y*x, fai=gz
	private float[][][] u=null;	// 2*y*x, westerly
	private float[][][] v=null;	// 2*y*x, southerly
	
	private VisualizeData vd=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public BoundedShallowWaterModelC(int dt,int y,int x){
		this.y=y;	this.x=x;	this.deltaT=dt;
		
		M_E=new float[2];
		
		h=new float[2][y+2][x+2];	H=new float[y+2][x+2];
		u=new float[2][y+2][x+2];	v=new float[2][y+2][x+2];
		
		deltaR=50000/deltaT;
		
		float dtmax=deltaR/350/1.414f;
		
		System.out.println("Bounded shallow water model integration parameters:");
		System.out.println("  delta dt : "+deltaT+" second(s), < "+dtmax+" s");
		
		vd =new VisualizeData(h[0]);
		vd.getModel().xmin=0;
		vd.getModel().xmax=x-1;
		vd.getModel().ymin=0;
		vd.getModel().ymax=y-1;
	}
	
	
	public void initialMirrorSurface(){
		for(int j=0,J=y+2;j<J;j++)
		for(int i=0,I=x+2;i<I;i++) h[0][j][i]=8000;
	}
	
	
	public void integration(int step,int interval,String output){
		System.out.println("\nstart integration...");
		
		int tag1=1,tag2=0;
		
		vd.setData(h[0]);
		
		for(int l=1;l<=step;l++){
			float[][] hres=h[tag1];	float[][] hsrc=h[tag2];
			float[][] ures=u[tag1];	float[][] usrc=u[tag2];
			float[][] vres=v[tag1];	float[][] vsrc=v[tag2];
			
			// predict h
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				hres[j][i]=hsrc[j][i]-
				(usrc[j][i]*(hsrc[j][i+1]+hsrc[j][i])-usrc[j][i-1]*(hsrc[j][i-1]+hsrc[j][i]))/deltaR/2-
				(vsrc[j][i]*(hsrc[j+1][i]+hsrc[j][i])-vsrc[j-1][i]*(hsrc[j-1][i]+hsrc[j][i]))/deltaR/2;
			}
			
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++) H[j][i]=(float)sqrt(hres[j][i]);
			
			// u-wind
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				ures[j][i]=usrc[j][i]*H[j][i]-(
					(H[j][i]+H[j][i+1])/2*(hres[j][i+1]-hres[j][i-1])
				)/deltaR-(
					(
						2*(usrc[j][i]+usrc[j][i+1])/2*(usrc[j][i]*(H[j][i]+H[j][i+1])/2+usrc[j][i+1]*(H[j][i+1]+H[j][i+2])/2)/2-
						2*(usrc[j][i]+usrc[j][i-1])/2*(usrc[j][i]*(H[j][i]+H[j][i+1])/2+usrc[j][i-1]*(H[j][i-1]+H[j][i  ])/2)/2
					)-
					usrc[j][i]*(H[j][i]+H[j][i+1])/2*(usrc[j][i+1]-usrc[j][i-1])/2
				)/deltaR/2-(
					(vsrc[j+1][i]*usrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*usrc[j-1][i]*H[j-1][i])+
					vsrc[j][i]*(usrc[j+1][i]*H[j+1][i]-usrc[j-1][i]*H[j-1][i])
				)/deltaR/2;
				
				ures[j][i]/=(H[j][i]+H[j][i+1])/2;
			}
			
			// v-wind
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				vres[j][i]=vsrc[j][i]*H[j][i]-(
						H[j][i]*(hres[j+1][i]-hres[j-1][i])
				)/deltaR-(
					(usrc[j][i+1]*vsrc[j][i+1]*H[j][i+1]-usrc[j][i-1]*vsrc[j][i-1]*H[j][i-1])+
					usrc[j][i]*(vsrc[j][i+1]*H[j][i+1]-vsrc[j][i-1]*H[j][i-1])
				)/deltaR/2-(
					(vsrc[j+1][i]*vsrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*vsrc[j-1][i]*H[j-1][i])+
					vsrc[j][i]*(vsrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*H[j-1][i])
				)/deltaR/2;
				
				vres[j][i]/=H[j][i];
			}
			
			// output and smooth if necessary
			if(l%interval==0){
				cTotalMassAndEnergy(tag1);
				
				System.out.println(
					"  output for "+l+" step, total mass and energy is: "+M_E[0]+"\t"+M_E[1]
				);
				
				vd.setData(h[tag1]);
				vd.updateModel(100);
			}
			
			int tmp=tag1;	tag1=tag2;	tag2=tmp;
		}
		
		System.out.println("finish integration.");
	}
	
	
	public void addHDisturbance(int x0,int y0,float rad,float amp){
		addDisturbance(x0,y0,rad,amp,h[0]);
	}
	
	public void addUDisturbance(int x0,int y0,float rad,float amp){
		addDisturbance(x0,y0,rad,amp,u[0]);
	}
	
	public void addVDisturbance(int x0,int y0,float rad,float amp){
		addDisturbance(x0,y0,rad,amp,v[0]);
	}
	
	private void addDisturbance(int x0,int y0,float rad,float amp,float[][] v){
		for(int j=0,J=y+2;j<J;j++)
		for(int i=0,I=x+2;i<I;i++)
		v[j][i]+=amp*(float)Math.exp(-Math.hypot(i-x0,j-y0)/rad);
	}
	
	private void cTotalMassAndEnergy(int tag){
		M_E[0]=0;	M_E[1]=0;	float tmp0=0,tmp1=0;
		
		for(int j=0;j<y;j++){
			for(int i=0;i<x;i++){
				tmp0+=h[tag][j][i];
				tmp1+=(h[tag][j][i]+u[tag][j][i]*u[tag][j][i]+v[tag][j][i]*v[tag][j][i])*h[tag][j][i]/2;
			}
			
			M_E[0]+=tmp0;	M_E[1]+=tmp1;
			
			tmp0=0;			tmp1=0;
		}
	}
	
	
	/*** getor and setor ***/
	public float[][] getInitialStateH(){ return h[0];}
	
	public float[][] getInitialStateU(){ return u[0];}
	
	public float[][] getInitialStateV(){ return v[0];}
	
	
	/** test*/
	public static void main(String[] args){
		try{
			BoundedShallowWaterModelC gswm=new BoundedShallowWaterModelC(60,100,100);
			
			gswm.initialMirrorSurface();
			//gswm.initialRossbyHaurwitz();
			
			//gswm.addQDisturbance(80,0,25,15,0.01f);
			gswm.addHDisturbance(50,50,10,3f);
			
			//gswm.integration(12*60,12*10,"d:/Gill.dat");
			
			gswm.integration(60*24*10,20,"d:/Gill.dat");
			
		}catch(Exception ex){ ex.printStackTrace();}
	}
}
