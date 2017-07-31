/**
 * @(#)ShallowWaterModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import static java.lang.Math.sqrt;
import miniufo.diagnosis.SpatialModel;
import miniufo.io.CtlDataWriteStream;


/**
 * shallow water wave model
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class BoundedShallowWaterModelA extends ShallowWaterModel2D{
	//
	private float f=2f*SpatialModel.EARTH_ROTATE_SPEED;
	private float deltaR=0;
	private float epsilonX=1e6f;
	private float epsilonY=epsilonX;
	
	private float[][] H=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public BoundedShallowWaterModelA(int dt,int y,int x,boolean writeResult){
		super(dt,y,x,writeResult);
		
		H=new float[y+2][x+2];
		
		deltaR=50000*2/deltaT;
		
		float dtmax=deltaR/350/1.414f;
		
		System.out.println("Bounded shallow water model integration parameters:");
		System.out.println("  delta dt : "+deltaT+" second(s), < "+dtmax+" s");
		
		if(!writeResult){
			vd.getModel().xmin=0;
			vd.getModel().xmax=x-1;
			vd.getModel().ymin=0;
			vd.getModel().ymax=y-1;
		}
	}
	
	public BoundedShallowWaterModelA(int dt,int y,int x){ this(dt,y,x,false);}
	
	
	public void integration(int step,int interval,String output){
		System.out.println("\nstart integration...");
		
		int tag1=1,tag2=0;
		
		float[][][][] odata=null;
		CtlDataWriteStream cdws=null;
		
		if(writeResult){
			odata=new float[1][1][y][x];
			
			out.setData(odata,true);
			
			cdws=new CtlDataWriteStream(output);
			
			copy(h[0],odata[0][0]);	cdws.writeData(out);
			copy(u[0],odata[0][0]);	cdws.writeData(out);
			copy(v[0],odata[0][0]);	cdws.writeData(out);
			
		}else{
			copyToDispBuf(h[0]);
			vd.setData(dispBuf);
		}
		
		for(int l=1;l<=step;l++){
			float[][] hres=h[tag1];	float[][] hsrc=h[tag2];
			float[][] ures=u[tag1];	float[][] usrc=u[tag2];
			float[][] vres=v[tag1];	float[][] vsrc=v[tag2];
			
			// predict h
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				hres[j][i]=hsrc[j][i]-
				(usrc[j][i+1]*hsrc[j][i+1]-usrc[j][i-1]*hsrc[j][i-1])/deltaR-
				(vsrc[j+1][i]*hsrc[j+1][i]-vsrc[j-1][i]*hsrc[j-1][i])/deltaR+
				(hsrc[j][i+1]+hsrc[j][i-1]-2*hsrc[j][i])/deltaR/deltaR/deltaT*epsilonX+
				(hsrc[j+1][i]+hsrc[j-1][i]-2*hsrc[j][i])/deltaR/deltaR/deltaT*epsilonY;
			}
			
			updateMassBoundaries(hres,tag2);
			
			for(int j=0,J=y+2;j<J;j++)
			for(int i=0,I=x+2;i<I;i++) H[j][i]=(float)sqrt(hres[j][i]);
			
			// wind
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				ures[j][i]=usrc[j][i]*H[j][i]-(
					H[j][i]*(hres[j][i+1]-hres[j][i-1])
				)/deltaR-f*vsrc[j][i]*H[j][i]-(
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
				)/deltaR+f*usrc[j][i]*H[j][i]-(
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
			
			updateWindBoundaries1(ures,vres,tag2);
			
			smooth5(h[tag1]); smooth5(u[tag1]); smooth5(v[tag1]);
			// output and smooth if necessary
			if(l%interval==0){
				
				cTotalMassAndEnergy(tag1);
				
				System.out.println(
					"  output for "+l+" step, total mass and energy is: "+totalMass+"\t"+totalEnergy
				);
				
				if(writeResult){
					copy(hres,odata[0][0]);	cdws.writeData(out);
					copy(ures,odata[0][0]);	cdws.writeData(out);
					copy(vres,odata[0][0]);	cdws.writeData(out);
					
				}else{
					copyToDispBuf(hres);
					vd.setData(dispBuf);
					vd.updateModel(5);
				}
			}
			
			int tmp=tag1;	tag1=tag2;	tag2=tmp;
		}
		
		if(writeResult){
			cdws.closeFile();
			writeCtl(step,interval,output);
		}
		
		System.out.println("finish integration.");
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
		float[][] hsrc=h[tag];
		float[][] usrc=u[tag];
		float[][] vsrc=v[tag];
		
		// deal with the south boundary
		for(int i=1;i<=x;i++){
			hres[0][i]=hsrc[0][i]-
			(vsrc[1][i]*hsrc[1][i]-vsrc[0][i]*hsrc[0][i])/deltaR/2;
		}
		
		// deal with the north boundary
		for(int i=1;i<=x;i++){
			hres[y+1][i]=hsrc[y+1][i]-
			(vsrc[y+1][i]*hsrc[y+1][i]-vsrc[y][i]*hsrc[y][i])/deltaR/2;
		}
		
		// deal with the west boundary
		for(int j=1;j<=y;j++){
			hres[j][0]=hsrc[j][0]-
			(usrc[j][1]*hsrc[j][1]-usrc[j][0]*hsrc[j][0])/deltaR/2;
		}
		
		// deal with the east boundary
		for(int j=1;j<=y;j++){
			hres[j][x+1]=hsrc[j][x+1]-
			(usrc[j][x+1]*hsrc[j][x+1]-usrc[j][x]*hsrc[j][x])/deltaR/2;
		}
		
		//periodicExtends(hres);
	}
	
	protected void updateWindBoundaries(float[][] ures,float[][] vres,int tag){
		// deal with the north boundaries
		for(int i=1;i<=x;i++){
			vres[y+1][i]=10f*(float)Math.sin(2.0*Math.PI*i/x);
		}
		
		//periodicExtends(ures);
		//periodicExtends(vres);
	}
	
	protected void updateWindBoundaries1(float[][] ures,float[][] vres,int l){
		// deal with the north boundaries
		for(int i=1;i<=x;i++){
			vres[y+1][i]=(1+l/10000f)*(float)Math.sin(2.0*Math.PI*i/x);
		}
	}
	
	
	private void copy(float[][] src,float[][] des){
		for(int j=1;j<=y;j++) System.arraycopy(src[j],1,des[j-1],0,x);
	}
	
	
	/** test*/
	public static void main(String[] args){
		try{
			BoundedShallowWaterModelA gswm=new BoundedShallowWaterModelA(60,140,140);
			
			gswm.initialMirrorSurface();
			//gswm.initialRossbyHaurwitz();
			
			//gswm.addQDisturbance(80,0,25,15,0.01f);
			//gswm.addHDisturbance(90,60,30,3f);
			//gswm.addHDisturbance(90,120,30,-3f);
			
			//gswm.integration(12*60,12*10,"d:/Gill.dat");
			
			gswm.integration(60*24*15,60*3,"d:/Gill.dat");
			
		}catch(Exception ex){ ex.printStackTrace();}
	}
}
