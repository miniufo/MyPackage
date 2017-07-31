/**
 * @(#)ShallowWaterModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import java.io.FileWriter;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


/**
 * shallow water wave model
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ShallowWaterModel1D{
	//
	private int x;
	
	private float deltaT=0;
	private float deltaX=0;
	private float ratio =0;
	
	private float[] Q  =null;	// 2*y*x, diabatic heating
	private float[] H  =null;
	private float[] buf=null;
	
	private double[] M_E=null;	// [0] means mass and [1] means energy
	
	private float[][] h=null;	// 2*x, fai=gz
	private float[][] u=null;	// 2*x, zonal wind
	
	Variable out=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public ShallowWaterModel1D(int dt,int x){
		this.x=x;	this.deltaT=dt;
		
		M_E=new double[2];
		buf=new float[x];
		
		Q=new float[x];	h=new float[2][x];
		H=new float[x];	u=new float[2][x];
		
		float dlon=(float)toRadians(360f/x);
		
		deltaX=EARTH_RADIUS*dlon;
		ratio =deltaT/(2*deltaX);
		
		float dtmax=deltaX/350;
		
		System.out.println("Global shallow water model integration parameters:");
		System.out.println("  delta time: "+deltaT+" second(s), < "+dtmax+" s");
		System.out.println("  delta lon : "+(float)toDegrees(dlon)+" degree");
		
		// set output variable state
		out=new Variable("out",true,new Range(1,1,1,1));
	}
	
	
	public void initialMirrorSurface(){ for(int i=0;i<x;i++) h[0][i]=8000;}
	
	
	public void integration(int step,int interval,String output){
		System.out.println("\nstart integration...");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(output);
		
		int tag1=1,tag2=0;
		
		float[][][][] odata=new float[1][1][1][];
		
		odata[0][0][0]=h[0];	out.setData(odata,true);
		
		// write initial condition
		odata[0][0][0]=h[0];	cdws.writeData(out);
		odata[0][0][0]=u[0];	cdws.writeData(out);
		
		for(int l=1;l<=step;l++){
			float[] hres=h[tag1];	float[] hsrc=h[tag2];
			float[] ures=u[tag1];	float[] usrc=u[tag2];
			
			/*************** predict h ***************/
			for(int i=1;i<x-1;i++)
			hres[i]=hsrc[i]-Q[i]*deltaT-(usrc[i+1]*hsrc[i+1]-usrc[i-1]*hsrc[i-1])*ratio;
			
			hres[0]=hsrc[0]-Q[0]*deltaT-(usrc[1]*hsrc[1]-usrc[x-1]*hsrc[x-1])*ratio;
			
			hres[x-1]=hsrc[x-1]-Q[x-1]*deltaT-(usrc[0]*hsrc[0]-usrc[x-2]*hsrc[x-2])*ratio;
			
			for(int i=0;i<x;i++) H[i]=(float)sqrt(hres[i]);
			
			/************ predict u-wind ************/
			for(int i=1;i<x-1;i++){
				ures[i]=usrc[i]*H[i]-(H[i]*(hres[i+1]-hres[i-1]))*ratio-(
					(usrc[i+1]*usrc[i+1]*H[i+1]-usrc[i-1]*usrc[i-1]*H[i-1])+
					usrc[i]*(usrc[i+1]*H[i+1]-usrc[i-1]*H[i-1])
				)*ratio/2;
				
				ures[i]/=H[i];
			}
			
			// left
			ures[0]=usrc[0]*H[0]-(H[0]*(hres[1]-hres[x-1]))*ratio-(
				(usrc[1]*usrc[1]*H[1]-usrc[x-1]*usrc[x-1]*H[x-1])+
				usrc[0]*(usrc[1]*H[1]-usrc[x-1]*H[x-1])
			)*ratio/2;
			
			ures[0]/=H[0];
			
			// right
			ures[x-1]=usrc[x-1]*H[x-1]-(H[x-1]*(hres[0]-hres[x-2]))*ratio-(
				(usrc[0]*usrc[0]*H[0]-usrc[x-2]*usrc[x-2]*H[x-2])+
				usrc[x-1]*(usrc[0]*H[0]-usrc[x-2]*H[x-2])
			)*ratio/2;
			
			ures[x-1]/=H[x-1];
			
			if(l%(interval/2)==0){ smooth(h[tag1]);	smooth(u[tag1]);}
			
			// output
			if(l%interval==0){
				cTotalMassAndEnergy(tag1);
				
				System.out.println(
					"  output for "+l+" step, total mass and energy is: "+M_E[0]+"\t"+M_E[1]
				);
				
				odata[0][0][0]=h[tag1];	cdws.writeData(out);
				odata[0][0][0]=u[tag1];	cdws.writeData(out);
			}
			
			// swap tag1 and tag2
			int tmp=tag1;	tag1=tag2;	tag2=tmp;
		}
		
		cdws.closeFile();
		
		try{
			// write ctl file
			StringBuffer sb=new StringBuffer();
			
			float dt=deltaT*interval/60;
			
			sb.append("dset ");	sb.append(output);	sb.append("\n");
			sb.append("undef -9999.0\n");
			sb.append("title shallow water model output\n");
			sb.append("xdef ");	sb.append(x);	sb.append(" linear   0 ");	sb.append(360f/x);	sb.append("\n");
			sb.append("ydef 1 linear 0 1\n");
			sb.append("zdef 1 levels 200\n");
			sb.append("tdef ");	sb.append(step/interval+1);	sb.append(" linear 00:00z1Jan2008 ");
			if(dt*interval<60){ sb.append((int)(dt*interval));	sb.append("mn\n");}
			else{ sb.append((int)(dt*interval/60));	sb.append("hr\n");}
			sb.append("vars 2\n");
			sb.append("h 0 99 geopotential\n");
			sb.append("u 0 99 u-wind\n");
			sb.append("endvars\n");
			
			FileWriter fw=new FileWriter(output.split("\\.")[0]+".ctl");
			
			fw.write(sb.toString());	fw.close();
			
		}catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
		
		System.out.println("finish integration.");
	}
	
	
	public void addQDisturbance(float lon0,float hlon,float amp){
		addDisturbance(lon0,hlon,amp,Q);
	}
	
	public void addHDisturbance(float lon0,float hlon,float amp){
		addDisturbance(lon0,hlon,amp,h[0]);
	}
	
	public void addUDisturbance(float lon0,float hlon,float amp){
		addDisturbance(lon0,hlon,amp,u[0]);
	}
	
	
	private void addDisturbance(float lon0,float hlon,float amp,float[] v){
		float dlon=360f/x;
		
		int istr=(int)((lon0-hlon)/dlon);
		int iend=(int)((lon0+hlon)/dlon);
		
		for(int i=istr;i<=iend;i++)
		v[i]+=amp*(float)(cos(PI*(i*dlon-lon0)/hlon)+1);
	}
	
	private void cTotalMassAndEnergy(int tag){
		M_E[0]=0;	M_E[1]=0;	double tmp0=0,tmp1=0;
		
		float[] tmpH=h[tag];
		float[] tmpU=u[tag];
		
		for(int i=0;i<x;i++){
			tmp0+=tmpH[i];
			tmp1+=(tmpH[i]+tmpU[i]*tmpU[i])*tmpH[i]/2;
		}
		
		M_E[0]+=tmp0;	M_E[1]+=tmp1;
	}
	
	private void smooth(float[] a){
		for(int i=1;i<x-1;i++) buf[i]=a[i]/2+(a[i-1]+a[i+1])/4;
		
		buf[ 0 ]=a[ 0 ]/2+(a[x-1]+a[1])/4;
		buf[x-1]=a[x-1]/2+(a[x-2]+a[0])/4;
	}
	
	
	/*** getor and setor ***/
	public float[] getInitialStateH(){ return h[0];}
	
	public float[] getInitialStateU(){ return u[0];}
	
	public float[] getInitialStateQ(){ return Q;}
	
	
	/** test*/
	public static void main(String[] args){
		try{
			ShallowWaterModel1D gswm=new ShallowWaterModel1D(20,144);
			
			gswm.initialMirrorSurface();
			
			//gswm.addQDisturbance(80,15,20,15,0.05f);
			gswm.addHDisturbance(140,50,1f);
			
			gswm.integration(3*60*24*10,60,"d:/Gill.dat");
			
		}catch(Exception ex){ ex.printStackTrace();}
	}
}
