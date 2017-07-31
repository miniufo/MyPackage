/**
 * @(#)ShallowWaterModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import miniufo.io.CtlDataWriteStream;
import miniufo.mathsphysics.FastFourier;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.tan;
import static java.lang.Math.sqrt;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.diagnosis.SpatialModel.EARTH_ROTATE_SPEED;


/**
 * shallow water wave model
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GlobalShallowWaterModelA extends ShallowWaterModel2D{
	//
	private int b=0;
	
	private int[] Ks=null;	// fourier filter
	
	private float[] f  =null;
	private float[] cos=null;
	private float[] tan=null;
	
	private float[][] Q=null;	// 2*y*x, diabatic heating
	private float[][] H=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public GlobalShallowWaterModelA(int dt,int y,int x,boolean writeResult){
		super(dt,y,x,writeResult);
		
		Q=new float[y+2][x+2];
		H=new float[y+2][x+2];
		
		float lat =0;
		float dlon=360f/x;
		float dlat=180f/(y-1);
		
		dlon=(float)toRadians(dlon);
		dlat=(float)toRadians(dlat);
		
		deltaY=EARTH_RADIUS*dlat;
		deltaX=EARTH_RADIUS*dlon;
		
		f  =new float[y];
		cos=new float[y];
		tan=new float[y];
		
		for(int j=0;j<y;j++){
			lat=(float)(-PI/2+dlat*j);	if(lat<=Math.toRadians(-70)) b++;
			
			  cos[j]=(float)cos(lat);
			  tan[j]=(float)tan(lat);
			    f[j]=2*EARTH_ROTATE_SPEED*(float)sin(lat);
		}
		
		float dtmax=deltaX*cos[1]/350/1.414f;
		
		System.out.println("Global 2D shallow water model integration parameters:");
		System.out.println("  delta dt : "+deltaT+" second(s), CFL dt: "+dtmax+" s");
		System.out.println("  delta lat: "+(float)toDegrees(dlat)+" degree");
		System.out.println("  delta lon: "+(float)toDegrees(dlon)+" degree");
		
		Ks=new int[b];
		for(int i=0;i<b;i++) Ks[i]=i;
	}
	
	public GlobalShallowWaterModelA(int dt,int y,int x){ this(dt,y,x,false);}
	
	
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
			
			float[][] deltaH=new float[y+2][x+2];
			float[][] deltaU=new float[y+2][x+2];
			float[][] deltaV=new float[y+2][x+2];
			
			/******************** process h equation ********************/
			cContinuityH(deltaH,tag2);
			
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++)
			hres[j][i]=hsrc[j][i]-deltaH[j][i]*deltaT;
			
			updateMassBoundaries(hres,tag2);
			
			for(int j=1;j<=y;j++)
			for(int i=0,I=x+2;i<I;i++)
			H[j][i]=(float)sqrt(hres[j][i]);
			
			/***************** process u and v equations *****************/
			cPressureGradientX(deltaU,tag2);	cPressureGradientY(deltaV,tag2);
			
			cAdvectionX(deltaU,tag2);			cAdvectionY(deltaV,tag2);
			cCoriolisX(deltaU,tag2);			cCoriolisY(deltaV,tag2);
			
			for(int j=1;j<=y;j++)
			for(int i=1;i<=x;i++){
				ures[j][i]=usrc[j][i]-deltaU[j][i]*deltaT/H[j][i];
				vres[j][i]=vsrc[j][i]-deltaV[j][i]*deltaT/H[j][i];
			}
			
			updateWindBoundaries(ures,vres,tag2);
			
			
			/*************** output and smooth if necessary ****************/
			if(l%(interval/4)==0){
				smoothPolar(hres); smoothPolar(ures); smoothPolar(vres);
			}
			
			if(l%interval==0){
				cTotalMassAndEnergy(tag1);
				
				System.out.println(
					"  output for "+l+" steps, total mass and energy are: "+
					totalMass+"\t"+totalEnergy
				);
				
				smooth5(h[tag1]); smooth5(u[tag1]); smooth5(v[tag1]);
				
				if(writeResult){
					copy(hres,odata[0][0]);	cdws.writeData(out);
					copy(ures,odata[0][0]);	cdws.writeData(out);
					copy(vres,odata[0][0]);	cdws.writeData(out);
					
				}else{
					copyToDispBuf(hres);
					vd.setData(dispBuf);
					vd.updateModel(10);
				}
			}
			
			// swap tags
			int tmp=tag1;	tag1=tag2;	tag2=tmp;
		}
		
		if(writeResult){
			cdws.closeFile();
			writeCtl(step,interval,output);
		}
		
		System.out.println("finish integration.");
	}
	
	public void initialRossbyHaurwitz(){
		super.initialRossbyHaurwitz();
		
		periodicExtends(h[0]);
		periodicExtends(u[0]);
		periodicExtends(v[0]);
	}
	
	
	protected void cTotalMassAndEnergy(int tag){
		totalMass=0;	totalEnergy=0;
		
		for(int j=1;j<=y;j++){
			float tmp0=0,tmp1=0;
			
			for(int i=1;i<=x;i++){
				tmp0+=h[tag][j][i];
				tmp1+=(h[tag][j][i]+u[tag][j][i]*u[tag][j][i]+v[tag][j][i]*v[tag][j][i])*h[tag][j][i]/2;
			}
			
			tmp0*=cos[j-1];	tmp1*=cos[j-1];
			
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
	
	
	private void cContinuityH(float[][] re,int tag){
		float[][] usrc=u[tag];	float[][] hsrc=h[tag];
		float[][] vsrc=v[tag];
		
		for(int j=2;j<y;j++)
		for(int i=1;i<=x;i++){
			re[j][i]=	// continuity term
			(usrc[j][i+1]*hsrc[j][i+1]-usrc[j][i-1]*hsrc[j][i-1])/deltaX/2/cos[j-1]+
			(vsrc[j+1][i]*cos[j]*hsrc[j+1][i]-vsrc[j-1][i]*cos[j-2]*hsrc[j-1][i])/deltaY/2/cos[j-1]-
			
			// diabatic term
			Q[j][i];
		}
		
		// deal with the polars
		float tmpS=0,tmpN=0;
		for(int i=1;i<=x;i++){
			tmpS+=hsrc[ 2 ][i]*vsrc[ 2 ][i];
			tmpN+=hsrc[y-1][i]*vsrc[y-1][i];
		}
		
		tmpS*=4f/deltaY/2/x;	tmpN*=4f/deltaY/2/x;
		
		for(int i=1;i<=x;i++){
			re[1][i]= tmpS;
			re[y][i]=-tmpN;
		}
	}
	
	private void cAdvectionY(float[][] re,int tag){
		float[][] usrc=u[tag];
		float[][] vsrc=v[tag];
		
		for(int j=2;j<y;j++)
		for(int i=1;i<=x;i++)
		re[j][i]+=
		(	// zonal advection term
			(usrc[j][i+1]*vsrc[j][i+1]*H[j][i+1]-usrc[j][i-1]*vsrc[j][i-1]*H[j][i-1])+
			usrc[j][i]*(vsrc[j][i+1]*H[j][i+1]-vsrc[j][i-1]*H[j][i-1])
		)/deltaX/cos[j-1]/2+
		
		(	// meridional advection term
			(vsrc[j+1][i]*cos[j]*vsrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*cos[j-2]*vsrc[j-1][i]*H[j-1][i])/cos[j-1]+
			vsrc[j][i]*(vsrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*H[j-1][i])
		)/deltaY/2;
	}
	
	private void cCoriolisY(float[][] re,int tag){
		float[][] usrc=u[tag];
		
		for(int j=2;j<y;j++)
		for(int i=1;i<=x;i++)
		re[j][i]+=(f[j-1]+usrc[j][i]*tan[j-1]/EARTH_RADIUS)*usrc[j][i]*H[j][i];
	}
	
	private void cPressureGradientY(float[][] re,int tag){
		float[][] hsrc=h[tag];
		
		for(int j=2;j<y;j++)
		for(int i=1;i<=x;i++)
		re[j][i]+=(H[j][i]*(hsrc[j+1][i]-hsrc[j-1][i]))/deltaY/2;
	}
	
	private void cAdvectionX(float[][] re,int tag){
		float[][] usrc=u[tag];
		float[][] vsrc=v[tag];
		
		for(int j=2;j<y;j++)
		for(int i=1;i<=x;i++)
		re[j][i]+=
		(	// zonal advection term
			(usrc[j][i+1]*usrc[j][i+1]*H[j][i+1]-usrc[j][i-1]*usrc[j][i-1]*H[j][i-1])+
			usrc[j][i]*(usrc[j][i+1]*H[j][i+1]-usrc[j][i-1]*H[j][i-1])
		)/deltaX/cos[j-1]/2+
		
		(	// meridional advection term
			(vsrc[j+1][i]*cos[j]*usrc[j+1][i]*H[j+1][i]-vsrc[j-1][i]*cos[j-2]*usrc[j-1][i]*H[j-1][i])/cos[j-1]+
			vsrc[j][i]*(usrc[j+1][i]*H[j+1][i]-usrc[j-1][i]*H[j-1][i])
		)/deltaY/2;
	}
	
	private void cCoriolisX(float[][] re,int tag){
		float[][] usrc=u[tag];
		float[][] vsrc=v[tag];
		
		for(int j=2;j<y;j++)
		for(int i=1;i<=x;i++)
		re[j][i]-=(f[j-1]+usrc[j][i]*tan[j-1]/EARTH_RADIUS)*vsrc[j][i]*H[j][i];
	}
	
	private void cPressureGradientX(float[][] re,int tag){
		float[][] hsrc=h[tag];
		
		for(int j=2;j<y;j++)
		for(int i=1;i<=x;i++)
		re[j][i]+=(H[j][i]*(hsrc[j][i+1]-hsrc[j][i-1]))/deltaX/2/cos[j-1];
	}
	
	
	private void copy(float[][] src,float[][] des){
		for(int j=1;j<=y;j++) System.arraycopy(src[j],1,des[j-1],0,x);
	}
	
	private void smoothPolar(float[][] a){
		FastFourier ff=new FastFourier(x);
		
		float[] buf=new float[x];
		
		for(int j=2;j<b;j++){
			System.arraycopy(a[j],1,buf,0,x);
			ff.fftMixedRadix(buf);
			
			float[] re1=ff.getResultRealPartCopy();
			float[] im1=ff.getResultImagePartCopy();
			
			for(int k=Ks[j-1]+1,K=x-Ks[j-1];k<K;k++) re1[k]=im1[k]=0;
			
			ff.ifftMixedRadix(re1,im1);
			
			System.arraycopy(ff.getResultRealPart(),0,a[j],1,x);
		}
		
		for(int j=y-1,J=y-b+1;j>J;j--){
			System.arraycopy(a[j],1,buf,0,x);
			ff.fftMixedRadix(buf);
			
			float[] re1=ff.getResultRealPartCopy();
			float[] im1=ff.getResultImagePartCopy();
			
			for(int k=Ks[y-j]+1,K=x-Ks[y-j];k<K;k++) re1[k]=im1[k]=0;
			
			ff.ifftMixedRadix(re1,im1);
			
			System.arraycopy(ff.getResultRealPart(),0,a[j],1,x);
		}
		
		periodicExtends(a);
	}
	
	private void periodicExtends(float[][] a){
		for(int j=0,J=y+2;j<J;j++){
			a[j][ 0 ]=a[j][x];
			a[j][x+1]=a[j][1];
		}
	}
	
	
	/** test*/
	public static void main(String[] args){
		try{
			GlobalShallowWaterModelA gswm=new GlobalShallowWaterModelA(60,73,144);
			
			gswm.initialMirrorSurface();
			//gswm.initialRossbyHaurwitz();
			
			//gswm.addQDisturbance(80,0,25,15,0.01f);
			gswm.addHDisturbance(70,40,6,20f);
			
			//gswm.integration(12*60,12*10,"d:/Gill.dat");
			
			/*
			GlobalWindFieldInSC gwf=new GlobalWindFieldInSC(SphericalSpacialModel.Model2P5);
			
			Variable sf=new Variable("sf",new miniufo.diagnosis.Range(1,1,73,144));
			sf.setValue(100000000);
			sf.add3DDisturbance(72,37,0,-5000000,8);
			
			Variable[] wnd=gwf.cRotationalWind(sf);
			
			GlobalHGTInSC ghgt=new GlobalHGTInSC(SphericalSpacialModel.Model2P5);
			Variable h=ghgt.cBalancedGeopotentialBySH(wnd[0],wnd[1]);
			h.plusEq(6000);
			
			gswm.h[0]=h.getData()[0][0];
			gswm.u[0]=wnd[0].getData()[0][0];
			gswm.v[0]=wnd[1].getData()[0][0];*/
			
			gswm.integration(60*24*20,30,"d:/Gill.dat");
			
		}catch(Exception ex){ ex.printStackTrace();}
	}
}
