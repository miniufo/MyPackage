/**
 * @(#)ShallowWaterModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import java.io.File;
import java.io.FileWriter;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataReadStream;
import miniufo.io.CtlDataWriteStream;
import static java.lang.Math.PI;
import static java.lang.Math.pow;
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
public final class GlobalShallowWaterModel{
	//
	private int y;
	private int x;
	private int b;
	
	private float epsilon;
	private float deltaT=0;
	private float deltaY=0;
	
	private float[] M_E=null;	// 0 means mass and 1 means energy
	
	private float[] ff =null;
	private float[] cos=null;
	private float[] tan=null;
	private float[] deltaX=null;
	
	private float[][] buff=null;	// buffer for smooth
	//private float[][] Terr=null;	// 2*y*x, terrain
	private float[][] Q   =null;	// 2*y*x, diabatic heating
	
	private float[][] A1U=null;
	private float[][] A2V=null;
	private float[][] A3H=null;
	private float[][] A3h=null;
	
	private float[][] B1U=null;
	private float[][] B2V=null;
	//private float[][] B3H=null;
	private float[][] B3h=null;
	
	private float[][][] U=null;
	private float[][][] V=null;
	private float[][][] H=null;
	
	private float[][][] h=null;	// 2*y*x, fai=gz
	private float[][][] u=null;	// 2*y*x, westerly
	private float[][][] v=null;	// 2*y*x, southerly
	
	Variable out=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public GlobalShallowWaterModel(int dt,int y,int x){
		this.y=y;	this.x=x;	this.deltaT=dt;
		
		M_E=new float[2];
		
		A1U=new float[y][x];	B1U=new float[y][x];
		A2V=new float[y][x];	B2V=new float[y][x];
		A3H=new float[y][x];	//B3H=new float[y][x];
		A3h=new float[y][x];	B3h=new float[y][x];
		
		 Q  =new float[y][x];
		//Terr=new float[y][x];
		buff=new float[y][x];
		
		h=new float[2][y][x];	H=new float[2][y][x];
		u=new float[2][y][x];	U=new float[2][y][x];
		v=new float[2][y][x];	V=new float[2][y][x];
		
		float lat =0;
		float dlon=360f/x;
		float dlat=180f/(y-1);
		
		dlon=(float)toRadians(dlon);
		dlat=(float)toRadians(dlat);
		
		deltaY=EARTH_RADIUS*dlat*2;
		
		ff =new float[y];
		cos=new float[y];
		tan=new float[y];
		deltaX=new float[y];
		
		for(int j=0;j<y;j++){
			lat=(float)(-PI/2+dlat*j);	if(lat<=Math.toRadians(-80)) b++;
			
			  cos[j]=(float)cos(lat);
			  tan[j]=(float)tan(lat);
			   ff[j]=2*EARTH_ROTATE_SPEED*(float)sin(lat);
			deltaX[j]=EARTH_RADIUS*dlon*cos[j]*2;
		}
		
		float dtmin=deltaX[1]/350/1.414f;
		
		System.out.println("Global shallow water model integration parameters:");
		System.out.println("  delta dt : "+deltaT+" second(s), < "+dtmin+" s");
		System.out.println("  delta lat: "+(float)toDegrees(dlat)+" degree");
		System.out.println("  delta lon: "+(float)toDegrees(dlon)+" degree");
		
		// set output variable state
		out=new Variable("out",true,new Range(1,1,1,1));
	}
	
	
	public void initialRossbyHaurwitz(){
		float a=EARTH_RADIUS,o=7.848e-6f,K=7.848e-6f,fai0=9.8f*8000;
		float dlon=360f/x,dlat=180f/(y-1);
		
		double A,B,C;
		double fai,lamda;
		
		for(int j=0;j<y;j++){
			fai=toRadians(-90+j*dlat);
			
			for(int i=0;i<x;i++){
				lamda=toRadians(i*dlon);
				
				A=o*(2*EARTH_ROTATE_SPEED+o)*pow(cos(fai),2)/2+K*K*pow(cos(fai),8)*(
					5*pow(cos(fai),2)+26-32*pow(cos(fai),-2)
				)/4;
				
				B=(EARTH_ROTATE_SPEED+o)*K/15*pow(cos(fai),4)*(26-25*pow(cos(fai),2));
				
				C=K*K*pow(cos(fai),8)*(5*pow(cos(fai),2)-6)/4;
				
				h[0][j][i]=fai0+(float)(a*a*A+a*a*cos(4*lamda)*B+a*a*cos(8*lamda)*C);
				
				u[0][j][i]=(float)(a*o*cos(fai)+a*K*pow(cos(fai),3)*(4*pow(sin(fai),2)-pow(cos(fai),2))*cos(4*lamda));
				
				v[0][j][i]=(float)(-a*K*4*pow(cos(fai),3)*sin(fai)*sin(4*lamda));
			}
		}
	}
	
	public void initialMirrorSurface(){
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) h[0][j][i]=8000;
	}
	
	
	public void integration(int step,int interval,String output){
		System.out.println("\nstart integration...");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(output);
		
		int tag1=1,tag2=0;
		
		float[][][][] odata=new float[1][1][][];
		
		odata[0][0]=h[0];	out.setData(odata,true);
		
		odata[0][0]=h[0];	cdws.writeData(out);
		odata[0][0]=u[0];	cdws.writeData(out);
		odata[0][0]=v[0];	cdws.writeData(out);
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			H[tag2][j][i]=(float)sqrt(h[tag2][j][i]);
			U[tag2][j][i]=u[tag2][j][i]*H[tag2][j][i];
			V[tag2][j][i]=v[tag2][j][i]*H[tag2][j][i];
		}
		
		for(int l=1;l<=step;l++){
			// predict
			//A3(A3h,h[tag2],tag2);	B3(B3h,A3h,h[tag2],tag2);
			A3(A3H,H[tag2],tag2);	//B3(B3H,A3H,H[tag2],tag2);
			A1(A1U,U[tag2],tag2);	//B1(B1U,A1U,U[tag2],tag2);
			A2(A2V,V[tag2],tag2);	//B2(B2V,A2V,V[tag2],tag2);
			
			cEpsilon(tag2);	epsilon*=deltaT;
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				//h[tag1][j][i]=h[tag2][j][i]-(A3H[j][i]+epsilon*B3H[j][i])*deltaT;
				h[tag1][j][i]=h[tag2][j][i]-Q[j][i]-A3H[j][i]*deltaT;
				H[tag1][j][i]=(float)sqrt(h[tag1][j][i]);
			}
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				//U[tag1][j][i]=U[tag2][j][i]-(A1U[j][i]+epsilon*B1U[j][i])*deltaT;
				U[tag1][j][i]=U[tag2][j][i]-A1U[j][i]*deltaT;
				u[tag1][j][i]=U[tag1][j][i]/H[tag1][j][i];
			}
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				//V[tag1][j][i]=V[tag2][j][i]-(A2V[j][i]+epsilon*B2V[j][i])*deltaT;
				V[tag1][j][i]=V[tag2][j][i]-A2V[j][i]*deltaT;
				v[tag1][j][i]=V[tag1][j][i]/H[tag1][j][i];
			}
			
			if(l%20==0){
				smooth(h[tag1]);	smooth(u[tag1]);	smooth(v[tag1]);
				smooth(H[tag1]);	smooth(U[tag1]);	smooth(V[tag1]);
			}
			
			// output and smooth if necessary
			if(l%interval==0){
				
				cTotalMassAndEnergy(tag1);
				
				System.out.println(
					"  output for "+l+" step, total mass and energy is: "+M_E[0]+"\t"+M_E[1]
				);
				
				odata[0][0]=h[tag1];	cdws.writeData(out);
				odata[0][0]=u[tag1];	cdws.writeData(out);
				odata[0][0]=v[tag1];	cdws.writeData(out);
			}
			
			tag1=(++tag1)%2;	tag2=(++tag2)%2;
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
			sb.append("ydef ");	sb.append(y);	sb.append(" linear -90 ");	sb.append(180f/(y-1));	sb.append("\n");
			sb.append("zdef 1 levels 200\n");
			sb.append("tdef ");	sb.append(step/interval+1);	sb.append(" linear 00:00z1Jan2008 ");
			if(dt<60){ sb.append((int)dt);	sb.append("mn\n");}
			else{ sb.append((int)(dt/60));	sb.append("hr\n");}
			sb.append("vars 3\n");
			sb.append("h 0 99 geopotential\n");
			sb.append("u 0 99 u-wind\n");
			sb.append("v 0 99 v-wind\n");
			sb.append("endvars\n");
			
			FileWriter fw=new FileWriter(output.split("\\.")[0]+".ctl");
			
			fw.write(sb.toString());	fw.close();
			
		}catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
		
		System.out.println("finish integration.");
	}
	
	
	public void addQDisturbance(float lon0,float lat0,float hlon,float hlat,float amp){
		addDisturbance(lon0,lat0,hlon,hlat,amp,Q);
	}
	
	public void addHDisturbance(float lon0,float lat0,float hlon,float hlat,float amp){
		addDisturbance(lon0,lat0,hlon,hlat,amp,h[0]);
	}
	
	public void addUDisturbance(float lon0,float lat0,float hlon,float hlat,float amp){
		addDisturbance(lon0,lat0,hlon,hlat,amp,u[0]);
	}
	
	public void addVDisturbance(float lon0,float lat0,float hlon,float hlat,float amp){
		addDisturbance(lon0,lat0,hlon,hlat,amp,v[0]);
	}
	
	
	private void addDisturbance(float lon0,float lat0,float hlon,float hlat,float amp,float[][] v){
		float dlon=360f/x,dlat=180f/(y-1);
		
		int jstr=(int)((lat0+90-hlat)/dlat),istr=(int)((lon0-hlon)/dlon);
		int jend=(int)((lat0+90+hlat)/dlat),iend=(int)((lon0+hlon)/dlon);
		
		for(int j=jstr;j<=jend;j++){ double lat=-90+j*dlat;
		for(int i=istr;i<=iend;i++){ double lon=i*dlon;
			v[j][i]+=amp*(float)(cos(PI*(lat-lat0)/hlat/2)*cos(PI*(lon-lon0)/hlon/2));
		}}
	}
	
	private void smoothPolar(float[][] a){
		for(int j=1;j<b;j++){
			for(int i=1;i<x-1;i++) buff[j][i]=a[j][i]/2+(a[j][i-1]+a[j][i+1])/4;
			
			buff[j][ 0 ]=a[j][ 0 ]/2+(a[j][x-1]+a[j][1])/4;
			buff[j][x-1]=a[j][x-1]/2+(a[j][x-2]+a[j][0])/4;
			
			for(int i=0;i<x;i++) a[j][i]=buff[j][i];
			
			for(int i=2;i<x-2;i++) buff[j][i]=a[j][i]/2+(a[j][i-2]+a[j][i+2]+a[j][i-1]+a[j][i+1])/8;
			
			buff[j][ 0 ]=a[j][ 0 ]/2+(a[j][x-2]+a[j][2]+a[j][x-1]+a[j][1])/8;
			buff[j][x-1]=a[j][x-1]/2+(a[j][x-3]+a[j][1]+a[j][x-2]+a[j][0])/8;
			buff[j][ 1 ]=a[j][ 1 ]/2+(a[j][x-1]+a[j][3]+a[j][0]+a[j][2])/8;
			buff[j][x-2]=a[j][x-2]/2+(a[j][x-4]+a[j][0]+a[j][x-3]+a[j][x-1])/8;
			
			for(int i=0;i<x;i++) a[j][i]=buff[j][i];
		}
		
		for(int j=y-2;j>=y-b-1;j--){
			for(int i=1;i<x-1;i++) buff[j][i]=a[j][i]/2+(a[j][i-1]+a[j][i+1])/4;
			
			buff[j][ 0 ]=a[j][ 0 ]/2+(a[j][x-1]+a[j][1])/4;
			buff[j][x-1]=a[j][x-1]/2+(a[j][x-2]+a[j][0])/4;
			
			for(int i=0;i<x;i++) a[j][i]=buff[j][i];
			
			for(int i=2;i<x-2;i++) buff[j][i]=a[j][i]/2+(a[j][i-2]+a[j][i+2]+a[j][i-1]+a[j][i+1])/8;
			
			buff[j][ 0 ]=a[j][ 0 ]/2+(a[j][x-2]+a[j][2]+a[j][x-1]+a[j][1])/8;
			buff[j][x-1]=a[j][x-1]/2+(a[j][x-3]+a[j][1]+a[j][x-2]+a[j][0])/8;
			buff[j][ 1 ]=a[j][ 1 ]/2+(a[j][x-1]+a[j][3]+a[j][0]+a[j][2])/8;
			buff[j][x-2]=a[j][x-2]/2+(a[j][x-4]+a[j][0]+a[j][x-3]+a[j][x-1])/8;
			
			for(int i=0;i<x;i++) a[j][i]=buff[j][i];
		}
	}
	
	private void smooth(float[][] a){
		for(int j=1;j<y-1;j++){
			for(int i=1;i<x-1;i++) buff[j][i]=a[j][i]/2+(a[j-1][i]+a[j+1][i]+a[j][i-1]+a[j][i+1])/8;
			
			buff[j][ 0 ]=a[j][ 0 ]/2+(a[j+1][0]+a[j-1][0]+a[j][x-1]+a[j][1])/8;
			buff[j][x-1]=a[j][x-1]/2+(a[j+1][x-1]+a[j-1][x-1]+a[j][x-2]+a[j][0])/8;
			
			for(int i=0;i<x;i++) a[j][i]=buff[j][i];
		}
	}
	
	private void cTotalMassAndEnergy(int tag){
		M_E[0]=0;	M_E[1]=0;	float tmp0=0,tmp1=0;
		
		for(int j=0;j<y;j++){
			for(int i=0;i<x;i++){
				tmp0+=h[tag][j][i];
				tmp1+=(h[tag][j][i]+u[tag][j][i]*u[tag][j][i]+v[tag][j][i]*v[tag][j][i])*h[tag][j][i]/2;
			}
			
			tmp0*=cos[j];	tmp1*=cos[j];
			
			M_E[0]+=tmp0;	M_E[1]+=tmp1;
			
			tmp0=0;			tmp1=0;
		}
	}
	
	private float cNorm2(float[][] a,float[][] b){
		float re=0,tmp=0;
		
		for(int j=0;j<y;j++){
			for(int i=0;i<x;i++) tmp+=a[j][i]*b[j][i];
			
			tmp*=cos[j];	re+=tmp;	tmp=0;
		}
		
		return re;
	}
	
	private void A3(float[][] re,float[][] a,int tag){
		for(int j=1;j<y-1;j++){
			for(int i=1;i<x-1;i++){
				re[j][i]=
				(U[tag][j][i+1]*a[j][i+1]-U[tag][j][i-1]*a[j][i-1])/deltaX[j]+
				(V[tag][j+1][i]*cos[j+1]*a[j+1][i]-V[tag][j-1][i]*cos[j-1]*a[j-1][i])/deltaY/cos[j];
			}
			
			// left
			re[j][0]=
			(U[tag][j][1]*a[j][1]-U[tag][j][x-1]*a[j][x-1])/deltaX[j]+
			(V[tag][j+1][0]*cos[j+1]*a[j+1][0]-V[tag][j-1][0]*cos[j-1]*a[j-1][0])/deltaY/cos[j];
			
			// right
			re[j][x-1]=
			(U[tag][j][0]*a[j][0]-U[tag][j][x-2]*a[j][x-2])/deltaX[j]+
			(V[tag][j+1][x-1]*cos[j+1]*a[j+1][x-1]-V[tag][j-1][x-1]*cos[j-1]*a[j-1][x-1])/deltaY/cos[j];
		}
		
		// deal with the polars
		float tmpS=0,tmpN=0;
		for(int i=0;i<x;i++){
			tmpS+=a[ 1 ][i]*V[tag][ 1 ][i];
			tmpN+=a[y-2][i]*V[tag][y-2][i];
		}
		
		tmpS*=4f/deltaY/x;	tmpN*=4f/deltaY/x;
		
		for(int i=0;i<x;i++){
			re[ 0 ][i]=tmpS;
			re[y-1][i]=-tmpN;
		}
		
		smoothPolar(re);
	}
	
	private void A2(float[][] re,float[][] a,int tag){
		for(int j=1;j<y-1;j++){
			for(int i=1;i<x-1;i++){
				re[j][i]=(
					H[tag][j][i]*(h[tag][j+1][i]-h[tag][j-1][i])
				)/deltaY+(
					(u[tag][j][i+1]*a[j][i+1]-u[tag][j][i-1]*a[j][i-1])+
					u[tag][j][i]*(a[j][i+1]-a[j][i-1])
				)/deltaX[j]/2+(
					(v[tag][j+1][i]*cos[j+1]*a[j+1][i]-v[tag][j-1][i]*cos[j-1]*a[j-1][i])/cos[j]+
					v[tag][j][i]*(a[j+1][i]-a[j-1][i])
				)/deltaY/2+(ff[j]+u[tag][j][i]*tan[j]/EARTH_RADIUS)*U[tag][j][i];
			}
			
			// left
			re[j][0]=(
				H[tag][j][0]*(h[tag][j+1][0]-h[tag][j-1][0])
			)/deltaY+(
				(u[tag][j][1]*a[j][1]-u[tag][j][x-1]*a[j][x-1])+
				u[tag][j][0]*(a[j][1]-a[j][x-1])
			)/deltaX[j]/2+(
				(v[tag][j+1][0]*cos[j+1]*a[j+1][0]-v[tag][j-1][0]*cos[j-1]*a[j-1][0])/cos[j]+
				v[tag][j][0]*(a[j+1][0]-a[j-1][0])
			)/deltaY/2+(ff[j]+u[tag][j][0]*tan[j]/EARTH_RADIUS)*U[tag][j][0];
			
			// right
			re[j][x-1]=(
				H[tag][j][x-1]*(h[tag][j+1][x-1]-h[tag][j-1][x-1])
			)/deltaY+(
				(u[tag][j][0]*a[j][0]-u[tag][j][x-2]*a[j][x-2])+
				u[tag][j][x-1]*(a[j][0]-a[j][x-2])
			)/deltaX[j]/2+(
				(v[tag][j+1][x-1]*cos[j+1]*a[j+1][x-1]-v[tag][j-1][x-1]*cos[j-1]*a[j-1][x-1])/cos[j]+
				v[tag][j][x-1]*(a[j+1][x-1]-a[j-1][x-1])
			)/deltaY/2+(ff[j]+u[tag][j][x-1]*tan[j]/EARTH_RADIUS)*U[tag][j][x-1];
		}
		
		smoothPolar(re);
	}
	
	private void A1(float[][] re,float[][] a,int tag){
		for(int j=1;j<y-1;j++){
			for(int i=1;i<x-1;i++){
				re[j][i]=(
					H[tag][j][i]*(h[tag][j][i+1]-h[tag][j][i-1])
				)/deltaX[j]+(
					(u[tag][j][i+1]*a[j][i+1]-u[tag][j][i-1]*a[j][i-1])+
					u[tag][j][i]*(a[j][i+1]-a[j][i-1])
				)/deltaX[j]/2+(
					(v[tag][j+1][i]*cos[j+1]*a[j+1][i]-v[tag][j-1][i]*cos[j-1]*a[j-1][i])/cos[j]+
					v[tag][j][i]*(a[j+1][i]-a[j-1][i])
				)/deltaY/2-(ff[j]+u[tag][j][i]*tan[j]/EARTH_RADIUS)*V[tag][j][i];
			}
			
			// left
			re[j][0]=(
				H[tag][j][0]*(h[tag][j][1]-h[tag][j][x-1])
			)/deltaX[j]+(
				(u[tag][j][1]*a[j][1]-u[tag][j][x-1]*a[j][x-1])+
				u[tag][j][0]*(a[j][1]-a[j][x-1])
			)/deltaX[j]/2+(
				(v[tag][j+1][0]*cos[j+1]*a[j+1][0]-v[tag][j-1][0]*cos[j-1]*a[j-1][0])/cos[j]+
				v[tag][j][0]*(a[j+1][0]-a[j-1][0])
			)/deltaY/2-(ff[j]+u[tag][j][0]*tan[j]/EARTH_RADIUS)*V[tag][j][0];
			
			// right
			re[j][x-1]=(
				H[tag][j][x-1]*(h[tag][j][0]-h[tag][j][x-2])
			)/deltaX[j]+(
				(u[tag][j][0]*a[j][0]-u[tag][j][x-2]*a[j][x-2])+
				u[tag][j][x-1]*(a[j][0]-a[j][x-2])
			)/deltaX[j]/2+(
				(v[tag][j+1][x-1]*cos[j+1]*a[j+1][x-1]-v[tag][j-1][x-1]*cos[j-1]*a[j-1][x-1])/cos[j]+
				v[tag][j][x-1]*(a[j+1][x-1]-a[j-1][x-1])
			)/deltaY/2-(ff[j]+u[tag][j][x-1]*tan[j]/EARTH_RADIUS)*V[tag][j][x-1];
		}
		
		smoothPolar(re);
	}
	
	void B3(float[][] BF,float[][] AF,float[][] F,int tag){
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) buff[j][i]=F[j][i]-AF[j][i]*deltaT;
		
		A3(BF,buff,tag);
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) BF[j][i]=(BF[j][i]-AF[j][i])/deltaT;
	}
	
	void B2(float[][] BF,float[][] AF,float[][] F,int tag){
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) buff[j][i]=F[j][i]-AF[j][i]*deltaT;
		
		A2(BF,buff,tag);
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) BF[j][i]=(BF[j][i]-AF[j][i])/deltaT;
	}
	
	void B1(float[][] BF,float[][] AF,float[][] F,int tag){
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) buff[j][i]=F[j][i]-AF[j][i]*deltaT;
		
		A1(BF,buff,tag);
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) BF[j][i]=(BF[j][i]-AF[j][i])/deltaT;
	}
	
	private void cEpsilon(int tag){
		float K1,K2,K3,K,a,b,c,d;
		
		K=cNorm2(B1U,U[tag])+cNorm2(B2V,V[tag])+cNorm2(B3h,h[tag]);
		
		if(K==0) K=0.001f;
		
		a=cNorm2(A1U,A1U);	b=cNorm2(A2V,A2V);	c=cNorm2(A3H,A3H);	d=cNorm2(A3h,A3h);
		
		K1=(a+b+c)/K;
		K2=(cNorm2(B1U,A1U)+cNorm2(B2V,A2V)+cNorm2(B3h,A3H))/K;
		K3=(float)(sqrt(a*cNorm2(B1U,B1U))+sqrt(b*cNorm2(B2V,B2V))+sqrt(d*cNorm2(B3h,B3h)))/K;
		
		epsilon=K1/(float)(1-deltaT*K2+sqrt(Math.abs(pow(1-deltaT*K2,2)-pow(deltaT*K3,2))));
		
		System.out.println(epsilon);
	}
	
	
	/*** getor and setor ***/
	public float[][] getInitialStateH(){ return h[0];}
	
	public float[][] getInitialStateU(){ return u[0];}
	
	public float[][] getInitialStateV(){ return v[0];}
	
	public float[][] getInitialStateQ(){ return Q;}
	
	
	/** test*/
	public static void main(String[] args){
		try{
			int x=360,y=181;
			GlobalShallowWaterModel gswm=new GlobalShallowWaterModel(1,y,x);
			
			float[][] h=gswm.getInitialStateH();
			float[][] u=gswm.getInitialStateU();
			
			//gswm.initialRossbyHaurwitz();
			
			/*
			gswm.initialMirrorSurface();
			
			int yy=30,xx=70,r=30;
			for(int j=-yy;j<yy;j++)
			for(int i=xx-r;i<xx+r;i++)
			h[j+(y-1)/2][i]+=(float)(3*(cos(toRadians(j/(float)yy*180.0))+1)*(cos(toRadians((i-xx)/(float)r*180f))+1));
			*/
			/*
			float a=6371220/4,o=7.848e-6f,K=7.848e-6f;
			for(int j=0;j<181;j++){
				double fai=toRadians(-90+j);
				
				for(int i=0;i<360;i++)
				u[j][i]=(float)(a*o*cos(fai)+a*K*pow(cos(fai),3)*(4*pow(sin(fai),2)-pow(cos(fai),2)));
			}
			
			for(int j=40;j<141;j++){
				double fai=toRadians(-90+j);
				
				for(int i=0;i<360;i++)
				u[j][i]-=(float)(3*(cos(fai/50*180)+1));
			}*/
			
			CtlDescriptor ctl=new CtlDescriptor(new File("d:/data/TyphoonWesterly/SWM/init.ctl"));
			
			Range r=new Range("t(2,2);z(2,2)",ctl);
			
			Variable uwnd=new Variable("u",true,r);
			Variable hgt =new Variable("h",true,r);
			
			CtlDataReadStream cdrs=new CtlDataReadStream(ctl);
			cdrs.readData(uwnd,hgt);	cdrs.closeFile();
			
			float[][][][] udata=uwnd.getData();
			float[][][][] hdata=hgt.getData();
			
			float d=hdata[0][0][0][0]-12000;
			
			for(int j=0;j<uwnd.getYCount();j++)
			for(int i=0;i<uwnd.getXCount();i++){
				u[j][i]=udata[0][0][j][i];
				h[j][i]=hdata[0][0][j][i];
				
				h[j][i]-=d;
			}
			
			/* westerly
			for(int j=-20;j<20;j++)
			for(int i=40;i<180;i++){
				u[j+90][i]+=(float)(3.5*(cos(toRadians((j-0)/20.0*180.0))+1)*(cos(toRadians((i-110f)/70f*180f))+1));
			}*/
			
			/* easterly */
			for(int j=-20;j<20;j++)
			for(int i=110;i<190;i++){
				u[j+90][i]+=(float)(2*(cos(toRadians((j-0)/20.0*180.0))+1)*(cos(toRadians((i-150f)/40f*180f))+1));
			}
			
			gswm.integration(60*60*25,60*60,"d:/easterly_cc.dat");
			
		}catch(Exception ex){ ex.printStackTrace();}
	}
}
