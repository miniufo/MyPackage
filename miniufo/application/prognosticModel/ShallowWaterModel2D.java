/**
 * @(#)ShallowWaterModel2D.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import java.io.FileWriter;

import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.visualize.VisualizeData;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.diagnosis.SpatialModel.EARTH_ROTATE_SPEED;


/**
 * shallow water wave model
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class ShallowWaterModel2D{
	//
	protected int y;
	protected int x;
	
	protected boolean writeResult=false;	// true for IO, false for display
	
	protected float deltaT=0;
	protected float deltaY=0;
	protected float deltaX=0;
	
	protected float totalMass=0;
	protected float totalEnergy=0;
	
	protected float[][] smthBuf=null;
	protected float[][] dispBuf=null;
	
	protected float[][][] h=null;	// 2*y*x, fai=gz
	protected float[][][] u=null;	// 2*y*x, westerly
	protected float[][][] v=null;	// 2*y*x, southerly
	
	protected Variable out=null;
	protected VisualizeData vd=null;
	
	
	/**
     * constructor
     *
     * @param	dt	delta t in seconds
     * @param	y	y-grid count
     * @param	x	x-grid count
     */
	public ShallowWaterModel2D(int dt,int y,int x,boolean writeResult){
		this.y=y;	this.x=x;	this.deltaT=dt;
		
		smthBuf=new float[y][x];
		dispBuf=new float[y][x];
		
		h=new float[2][y+2][x+2];
		u=new float[2][y+2][x+2];
		v=new float[2][y+2][x+2];
		
		this.writeResult=writeResult;
		
		// set output variable state
		if(writeResult){
			out=new Variable("out",true,new Range(1,1,1,1));
		}else{
			vd =new VisualizeData(h[0]);
			vd.getModel().xmin=0;
			vd.getModel().xmax=360-360f/x;
			vd.getModel().ymin=-90;
			vd.getModel().ymax= 90;
		}
	}
	
	public ShallowWaterModel2D(int dt,int y,int x){
		this(dt,y,x,false);
	}
	
	
	public void initialRossbyHaurwitz(){
		float a=EARTH_RADIUS,o=7.848e-6f,K=7.848e-6f,fai0=9.8f*8000;
		float dlon=360f/x,dlat=180f/(y-1);
		
		for(int j=1;j<=y;j++){
			double fai=Math.toRadians(-90+(j-1)*dlat);
			
			for(int i=1;i<=x;i++){
				double lamda=Math.toRadians(i*dlon);
				
				double A=o*(2*EARTH_ROTATE_SPEED+o)*pow(cos(fai),2)/2+K*K*pow(cos(fai),8)*(
					5*pow(cos(fai),2)+26-32*pow(cos(fai),-2)
				)/4;
				
				double B=(EARTH_ROTATE_SPEED+o)*K/15*pow(cos(fai),4)*(26-25*pow(cos(fai),2));
				
				double C=K*K*pow(cos(fai),8)*(5*pow(cos(fai),2)-6)/4;
				
				h[0][j][i]=fai0+(float)(a*a*A+a*a*cos(4*lamda)*B+a*a*cos(8*lamda)*C);
				
				u[0][j][i]=(float)(a*o*cos(fai)+a*K*pow(cos(fai),3)*(4*pow(sin(fai),2)-pow(cos(fai),2))*cos(4*lamda));
				
				v[0][j][i]=(float)(-a*K*4*pow(cos(fai),3)*sin(fai)*sin(4*lamda));
			}
		}
	}
	
	public void initialMirrorSurface(){
		for(int j=0,J=y+2;j<J;j++)
		for(int i=0,I=x+2;i<I;i++) h[0][j][i]=8000;
	}
	
	
	/**
     * model integration
     *
     * @param	step		total steps to be integrated
     * @param	interval	interval for output or display
     * @param	output		output path
     */
	public abstract void integration(int step,int interval,String output);
	
	
	public void addHDisturbance(int x0,int y0,float rad,float amp){
		addDisturbance(x0,y0,rad,amp,h[0]);
	}
	
	public void addUDisturbance(int x0,int y0,float rad,float amp){
		addDisturbance(x0,y0,rad,amp,u[0]);
	}
	
	public void addVDisturbance(int x0,int y0,float rad,float amp){
		addDisturbance(x0,y0,rad,amp,v[0]);
	}
	
	public void addDisturbance(int x0,int y0,float rad,float amp,float[][] v){
		for(int j=1;j<=y;j++)
		for(int i=1;i<=x;i++)
		v[j][i]+=amp*(float)Math.exp(-Math.hypot(i-x0,j-y0)/rad);
	}
	
	
	protected abstract void cTotalMassAndEnergy(int tag);
	
	protected abstract void updateMassBoundaries(float[][] h,int tag);
	
	protected abstract void updateWindBoundaries(float[][] u,float[][] v,int tag);
	
	
	protected void smooth5(float[][] a){
		for(int j=1;j<=y;j++){
			for(int i=1;i<=x;i++){
				smthBuf[j-1][i-1]=a[j][i]/2+(a[j-1][i]+a[j+1][i]+a[j][i-1]+a[j][i+1])/8;
				//if(j==1&&i==1)
				//System.out.println(a[j][i]+"\t"+a[j-1][i]+"\t"+a[j+1][i]+"\t"+a[j][i-1]+"\t"+a[j][i+1]);
			}
			
			System.arraycopy(smthBuf[j-1],0,a[j],1,x);
		}
	}
	
	protected void copyToDispBuf(float[][] a){
		for(int j=1;j<=y;j++) System.arraycopy(a[j],1,dispBuf[j-1],0,x);
	}
	
	protected void writeCtl(int step,int interval,String output){
		// write ctl file
		StringBuffer sb=new StringBuffer();
		
		float dt=deltaT*interval/60;
		
		sb.append("dset ");	sb.append(output);	sb.append("\n");
		sb.append("undef -9999.0\n");
		sb.append("title shallow water model output\n");
		sb.append("xdef "); sb.append(x); sb.append(" linear   0 "); sb.append(360f/x); sb.append("\n");
		sb.append("ydef "); sb.append(y); sb.append(" linear -90 "); sb.append(180f/(y-1)); sb.append("\n");
		sb.append("zdef 1 levels 200\n");
		sb.append("tdef ");	sb.append(step/interval+1);	sb.append(" linear 00:00z1Jan2008 ");
		if(dt*interval<60){ sb.append((int)(dt*interval));	sb.append("mn\n");}
		else{ sb.append((int)(dt*interval/60));	sb.append("hr\n");}
		sb.append("vars 3\n");
		sb.append("h 0 99 geopotential\n");
		sb.append("u 0 99 u-wind\n");
		sb.append("v 0 99 v-wind\n");
		sb.append("endvars\n");
		
		try{
			FileWriter fw=new FileWriter(output.split("\\.")[0]+".ctl");
			
			fw.write(sb.toString());	fw.close();
			
		}catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
	}
}
