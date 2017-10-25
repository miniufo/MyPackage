/**
 * @(#)GlobalGillModelInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import miniufo.application.EquationInSphericalCoordinate;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;


/**
 * response of geopotential field to a specified heating source
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GlobalGillModelInSC extends EquationInSphericalCoordinate{
	//
	private  int  maxloop=5000;		// max loop of SOR
	private float coeffic=8e-6f;	// coefficient of dissipation
	private float gheight=2000f;	// mean depth of the fluid (equals Cg^2)
	
	private float[] coef=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public GlobalGillModelInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isPeriodicX())
		throw new IllegalArgumentException("Not a zonal periodic model");
		
		x=ssm.getXCount();	y=ssm.getYCount();	coef=new float[y];
		
		for(int j=0;j<y;j++) coef[j]=coeffic*lsin[j];
	}
	
	
	/*** getor and setor ***/
	public int getMaxLoop(){ return maxloop;}
	
	public float getDissipationCoefficient(){ return coeffic;}
	
	public float getBasicGeopotential(){ return gheight;}
	
	public void setDissipationCoefficient(float coeff){
		if(coeff<0) throw new IllegalArgumentException("coeff should larger than 0");
		
		coeffic=coeff;
		
		for(int j=0;j<y;j++) coef[j]=coeffic*lcos[j];
	}
	
	public void setBasicGeopotential(float height){
		if(height<=0) throw new IllegalArgumentException("height should larger than 0");
		
		gheight=height;
	}
	
	
	/**
     * generate a heating source in a specified position (lon0,lat0)
     * with the radius in longitude (hlon) and latitude (hlat) and amplitude (amp)
     *
     * @param	Q		heating source variable
     * @param	lon0	center longitude of heating source 
     * @param	lat0	center latitude of heating source
     * @param	hlon	radius of heating source in longitude
     * @param	hlat	radius of heating source in latitude
     * @param	amp		amplitude (strength) of the heating source
     */
	public void addHeatingSource(Variable Q,float lon0,float lat0,float hlon,float hlat,float amp){
		if(Q.getXCount()!=x||Q.getYCount()!=y)
			throw new IllegalArgumentException("not a global variable");
		
		if(!Q.isTFirst()) throw new IllegalArgumentException("Q should be a T-First variable");
		
		if(Q.getTCount()!=1||Q.getZCount()!=1)
			throw new IllegalArgumentException("area-type variable require");
		
		float[][] Qdata=Q.getData()[0][0];
		
		float dlon=360f/x,dlat=180f/(y-1);
		
		int jstr=(int)((lat0+90-hlat)/dlat),istr=(int)((lon0-hlon)/dlon);
		int jend=(int)((lat0+90+hlat)/dlat),iend=(int)((lon0+hlon)/dlon);
		
		for(int j=jstr;j<=jend;j++){ double lat=-90+j*dlat;
		for(int i=istr;i<=iend;i++){ double lon=i*dlon;
			Qdata[j][i]+=amp*(float)(cos(PI*(lat-lat0)/hlat/2)*cos(PI*(lon-lon0)/hlon/2));
		}}
	}
	
	public void addHeatingSource2(Variable Q,float lon0,float lat0,float arg,float amp){
		if(Q.getXCount()!=x||Q.getYCount()!=y)
			throw new IllegalArgumentException("not a global variable");
		
		if(!Q.isTFirst()) throw new IllegalArgumentException("Q should be a T-First variable");
		
		if(Q.getTCount()!=1||Q.getZCount()!=1)
			throw new IllegalArgumentException("area-type variable require");
		
		if(arg<=1)
			throw new IllegalArgumentException("arg should be larger than 1");
		
		float[][] Qdata=Q.getData()[0][0];
		
		float dlon=360f/x,dlat=180f/(y-1);
		
		for(int j=0;j<y;j++){ double lat=-90+j*dlat,dislat=lat-lat0;dislat*=dislat;
		for(int i=0;i<x;i++){ double lon=i*dlon,dislon=lon-lon0;dislon*=dislon;
			Qdata[j][i]+=amp*(float)exp(-(dislat+dislon)/arg);
		}}
	}
	
	
	/**
     * the value of Q beyond the key region is set to zero
     *
     * @param	Q		heating source variable
     * @param	lon1	lower longitude of the key region
     * @param	lon2	higher longitude of the key region
     * @param	lat1	lower latitude of the key region
     * @param	lat2	higher latitude of the key region
     */
	public void setKeyRegion(Variable Q,float lon1,float lon2,float lat1,float lat2){
		if(lon1>=lon2||lat1>=lat2) throw new IllegalArgumentException("invalid region");
		
		if(Q.getXCount()!=x||Q.getYCount()!=y)
			throw new IllegalArgumentException("not a global variable");
		
		if(!Q.isTFirst()) throw new IllegalArgumentException("Q should be a T-First variable");
		
		if(Q.getTCount()!=1||Q.getZCount()!=1)
			throw new IllegalArgumentException("area-type variable require");
		
		float dlon=360f/x,dlat=180f/(y-1);
		
		float[][] Qdata=Q.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			float lon=i*dlon,lat=-90+j*dlat;
			if(lon<lon1||lon>lon2||lat<lat1||lat>lat2) Qdata[j][i]=0;
		}
	}
	
	
	/**
     * compute the wind fields according to geopotential field
     *
     * @param	h	geopotential field
     * 
     * @return	re	wind fields, re[0] is u-wind and re[1] is v-wind
     */
	public Variable[] cWindFields(Variable h){
		if(h.getXCount()!=x||h.getYCount()!=y)
			throw new IllegalArgumentException("not a global variable");
		
		if(!h.isTFirst()) throw new IllegalArgumentException("Q should be a T-First variable");
		
		if(h.getTCount()!=1||h.getZCount()!=1)
			throw new IllegalArgumentException("area-type variable require");
		
		Variable[] re=new Variable[2];
		re[0]=new Variable("u",h);
		re[1]=new Variable("v",h);
		
		float[][] udata=re[0].getData()[0][0];
		float[][] vdata=re[1].getData()[0][0];
		float[][] hdata=h.getData()[0][0];
		
		for(int j=1;j<y-1;j++){
			float de=coeffic*coeffic+f1[j]*f1[j];
			
			for(int i=1;i<x-1;i++){
				udata[j][i]=(
					- f1[j] *(hdata[j+1][i]-hdata[j-1][i])/(dy   +dy   )
					-coeffic*(hdata[j][i+1]-hdata[j][i-1])/(dxs[j]+dxs[j])
				)/de;
				
				vdata[j][i]=(
					  f1[j] *(hdata[j][i+1]-hdata[j][i-1])/(dxs[j]+dxs[j])
					-coeffic*(hdata[j+1][i]-hdata[j-1][i])/(dy   +dy   )
				)/de;
			}
			
			// east and west boundary
			udata[j][0]=(
				- f1[j] *(hdata[j+1][0]-hdata[j-1][0])/(dy   +dy   )
				-coeffic*(hdata[j  ][1]-hdata[j][x-1])/(dxs[j]+dxs[j])
			)/de;
			
			vdata[j][0]=(
				  f1[j] *(hdata[j  ][1]-hdata[j][x-1])/(dxs[j]+dxs[j])
				-coeffic*(hdata[j+1][0]-hdata[j-1][0])/(dy   +dy   )
			)/de;
			 
			udata[j][x-1]=(
				- f1[j] *(hdata[j+1][x-1]-hdata[j-1][x-1])/(dy   +dy   )
				-coeffic*(hdata[j  ][0  ]-hdata[j  ][x-2])/(dxs[j]+dxs[j])
			)/de;
			
			vdata[j][x-1]=(
				  f1[j] *(hdata[j  ][0  ]-hdata[j  ][x-2])/(dxs[j]+dxs[j])
				-coeffic*(hdata[j+1][x-1]-hdata[j-1][x-1])/(dy   +dy   )
			)/de;
		}
		
		return re;
	}
	
	
	/**
     * successive over relaxation iteration
     *
     * @param	err		err for stopping iteration
     * @param	h		geopotential field
     * @param	Q		heating source
     */
	public void SOR(float err,Variable h,Variable Q){
		if(!h.isLike(Q)) throw new IllegalArgumentException("dimensions not same");
		
		System.out.println("\nStart SORing...");
		
		y=Q.getYCount();	x=Q.getXCount();
		
		if(err<0) throw new IllegalArgumentException();
		
		tstart=h.getRange().getTRange()[0];
		
		float[][] hdata=h.getData()[0][0];
		float[][] Qdata=Q.getData()[0][0];
		
		if(err<=0) throw new IllegalArgumentException("Error limit should larger than 0");
		
		int xx=0;
		
		/*** calculate optimizing argument ***/
		float epsilon=(float)(pow(sin(PI/(2*x+2)),2)+pow(sin(PI/(2*y+2)),2));
		float opt_arg=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		opt_arg=1.3f;
		
		float coeff2=coeffic*coeffic,dy2=dy*dy;
		
		float[][] R=new float[y-2][x];
		
		do{
			float maxerr=0;
			
			for(int j=11;j<y-11;j++){
				float ff2=f1[j]*f1[j],dx2=dxs[j]*dxs[j];
				float A=coeffic*gheight/(coeff2+ff2);
				float B=- ( coeff2-ff2 )*Beta[j]*gheight/(coeff2+ff2)/(coeff2+ff2);
				float C=-2*coeffic*f1[j]*Beta[j]*gheight/(coeff2+ff2)/(coeff2+ff2);
				
				/*** East boundary ***/
				R[j-1][0]=A*(
					(
						(hdata[j  ][1]-hdata[j][0])-(hdata[j][0]-hdata[j][x-1])
						
					)/dx2+(
						(hdata[j+1][0]-hdata[j][0])-(hdata[j][0]-hdata[j-1][0])
						
					)/dy2
				)+(
					B*(hdata[j  ][1]-hdata[j][x-1])/(dxs[j]+dxs[j])+
					C*(hdata[j+1][0]-hdata[j-1][0])/(dy   +dy   )
				
				)-coeffic*hdata[j][0]-Qdata[j][0];
				
				R[j-1][0]*=opt_arg/(2*A/dx2+2*A/dy2-coeffic);
				
				float tmp=abs(R[j-1][0]);
				if(tmp>maxerr) maxerr=tmp;
				
				hdata[j][0]+=R[j-1][0];
				
				/*** internal area ***/
				for(int i=1;i<x-1;i++){
					R[j-1][i]=A*(
						(
							(hdata[j][i+1]-hdata[j][i])-(hdata[j][i]-hdata[j][i-1])
							
						)/dx2+(
							(hdata[j+1][i]-hdata[j][i])-(hdata[j][i]-hdata[j-1][i])
							
						)/dy2
					)+(
						B*(hdata[j][i+1]-hdata[j][i-1])/(dxs[j]+dxs[j])+
						C*(hdata[j+1][i]-hdata[j-1][i])/(dy   +dy   )
						
					)-coeffic*hdata[j][i]-Qdata[j][i];
					
					R[j-1][i]*=opt_arg/(2*A/dx2+2*A/dy2-coeffic);
					
					tmp=abs(R[j-1][i]);
					if(tmp>maxerr) maxerr=tmp;
					
					hdata[j][i]+=R[j-1][i];
				}
				
				/*** West boundary ***/
				R[j-1][x-1]=A*(
					(
						(hdata[j  ][0  ]-hdata[j][x-1])-(hdata[j][x-1]-hdata[j  ][x-2])
						
					)/dx2+(
						(hdata[j+1][x-1]-hdata[j][x-1])-(hdata[j][x-1]-hdata[j-1][x-1])
						
					)/dy2
				)+(
					B*(hdata[j  ][0  ]-hdata[j  ][x-2])/(dxs[j]+dxs[j])+
					C*(hdata[j+1][x-1]-hdata[j-1][x-1])/(dy   +dy   )
					
				)-coeffic*hdata[j][x-1]-Qdata[j][x-1];
				
				R[j-1][x-1]*=opt_arg/(2*A/dx2+2*A/dy2-coeffic);
				
				tmp=abs(R[j-1][x-1]);
				if(tmp>maxerr) maxerr=tmp;
				
				hdata[j][x-1]+=R[j-1][x-1];
			}
			
			if(maxerr<=err||xx++>=maxloop) break;
			
		}while(true);
		
		if(xx<maxloop) System.out.println("loops "+xx);
		else System.out.println("loops "+xx);
		
		System.out.println("Finish SORing.");
	}
	
	public void SORObs(float err,Variable h,Variable hm,Variable Q){
		if(!h.isLike(Q)) throw new IllegalArgumentException("dimensions not same");
		
		System.out.println("\nStart SORing...");
		
		y=Q.getYCount();	x=Q.getXCount();
		
		if(err<0) throw new IllegalArgumentException();
		
		tstart=h.getRange().getTRange()[0];
		
		float[][] hdata= h.getData()[0][0];
		float[][] mdata=hm.getData()[0][0];
		float[][] Qdata= Q.getData()[0][0];
		
		if(err<=0) throw new IllegalArgumentException("Error limit should larger than 0");
		
		int xx=0;
		
		/*** calculate optimizing argument ***/
		float epsilon=(float)(pow(sin(PI/(2*x+2)),2)+pow(sin(PI/(2*y+2)),2));
		float opt_arg=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		opt_arg=1.3f;
		
		float coeff2=coeffic*coeffic,dy2=dy*dy;
		
		do{
			float maxerr=0;
			
			for(int j=11;j<y-11;j++){
				float ff2=f1[j]*f1[j],dx2=dxs[j]*dxs[j];
				float A=coeffic*mdata[j][0]/(coeff2+ff2);
				float B=- ( coeff2-ff2 )*Beta[j]*mdata[j][0]/(coeff2+ff2)/(coeff2+ff2);
				float C=-2*coeffic*f1[j]*Beta[j]*mdata[j][0]/(coeff2+ff2)/(coeff2+ff2);
				
				/*** East boundary ***/
				float R=A*(
					(
						(hdata[j  ][1]-hdata[j][0])-(hdata[j][0]-hdata[j][x-1])
						
					)/dx2+(
						(hdata[j+1][0]-hdata[j][0])-(hdata[j][0]-hdata[j-1][0])
						
					)/dy2
				)+(
					B*(hdata[j  ][1]-hdata[j][x-1])/(dxs[j]+dxs[j])+
					C*(hdata[j+1][0]-hdata[j-1][0])/(dy+dy)
				
				)-coeffic*hdata[j][0]-Qdata[j][0];
				
				R*=opt_arg/(2*A/dx2+2*A/dy2-coeffic);
				
				float tmp=abs(R);
				if(tmp>maxerr) maxerr=tmp;
				
				hdata[j][0]+=R;
				
				/*** internal area ***/
				for(int i=1;i<x-1;i++){
					A=coeffic*mdata[j][i]/(coeff2+ff2);
					B=- ( coeff2-ff2 )*Beta[j]*mdata[j][i]/(coeff2+ff2)/(coeff2+ff2);
					C=-2*coeffic*f1[j]*Beta[j]*mdata[j][i]/(coeff2+ff2)/(coeff2+ff2);
					
					R=A*(
						(
							(hdata[j][i+1]-hdata[j][i])-(hdata[j][i]-hdata[j][i-1])
							
						)/dx2+(
							(hdata[j+1][i]-hdata[j][i])-(hdata[j][i]-hdata[j-1][i])
							
						)/dy2
					)+(
						B*(hdata[j][i+1]-hdata[j][i-1])/(dxs[j]+dxs[j])+
						C*(hdata[j+1][i]-hdata[j-1][i])/(dy   +dy   )
						
					)-coeffic*hdata[j][i]-Qdata[j][i];
					
					R*=opt_arg/(2*A/dx2+2*A/dy2-coeffic);
					
					tmp=abs(R);
					if(tmp>maxerr) maxerr=tmp;
					
					hdata[j][i]+=R;
				}
				
				/*** West boundary ***/
				A=coeffic*mdata[j][x-1]/(coeff2+ff2);
				B=- ( coeff2-ff2 )*Beta[j]*mdata[j][x-1]/(coeff2+ff2)/(coeff2+ff2);
				C=-2*coeffic*f1[j]*Beta[j]*mdata[j][x-1]/(coeff2+ff2)/(coeff2+ff2);
				
				R=A*(
					(
						(hdata[j  ][0  ]-hdata[j][x-1])-(hdata[j][x-1]-hdata[j  ][x-2])
						
					)/dx2+(
						(hdata[j+1][x-1]-hdata[j][x-1])-(hdata[j][x-1]-hdata[j-1][x-1])
						
					)/dy2
				)+(
					B*(hdata[j  ][0  ]-hdata[j  ][x-2])/(dxs[j]+dxs[j])+
					C*(hdata[j+1][x-1]-hdata[j-1][x-1])/(dy   +dy   )
					
				)-coeffic*hdata[j][x-1]-Qdata[j][x-1];
				
				R*=opt_arg/(2*A/dx2+2*A/dy2-coeffic);
				
				tmp=abs(R);
				if(tmp>maxerr) maxerr=tmp;
				
				hdata[j][x-1]+=R;
			}
			
			if(maxerr<=err||xx++>=maxloop) break;
			
		}while(true);
		
		if(xx<maxloop) System.out.println("loops "+xx);
		else System.out.println("loops "+xx);
		
		System.out.println("Finish SORing.");
	}
	
	public void SORObsLimitArea(float err,Variable h,Variable hm,Variable Q,int x1,int x2,int y1,int y2){
		if(!h.isLike(Q)) throw new IllegalArgumentException("dimensions not same");
		
		System.out.println("\nStart SORing...");
		
		y=Q.getYCount();	x=Q.getXCount();
		
		if(err<0) throw new IllegalArgumentException();
		
		tstart=h.getRange().getTRange()[0];
		
		float[][] hdata= h.getData()[0][0];
		float[][] mdata=hm.getData()[0][0];
		float[][] Qdata= Q.getData()[0][0];
		
		if(err<=0) throw new IllegalArgumentException("Error limit should larger than 0");
		
		int xx=0;
		
		/*** calculate optimizing argument ***/
		float epsilon=(float)(pow(sin(PI/(2*(x2-x1+1)+2)),2)+pow(sin(PI/(2*(y2-y1+1)+2)),2));
		float opt_arg=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		//System.out.println("opt_arg: "+opt_arg);
		opt_arg=1.3f;
		
		float coeff2=coeffic*coeffic,dy2=dy*dy;
		
		do{
			float maxerr=0;
			
			for(int j=y1;j<=y2;j++){
				float ff2=f1[j]*f1[j],dx2=dxs[j]*dxs[j];
				float A=coeffic*mdata[j][0]/(coeff2+ff2);
				float B=- ( coeff2-ff2 )*Beta[j]*mdata[j][0]/(coeff2+ff2)/(coeff2+ff2);
				float C=-2*coeffic*f1[j]*Beta[j]*mdata[j][0]/(coeff2+ff2)/(coeff2+ff2);
				
				/*** internal area ***/
				for(int i=x1;i<=x2;i++){
					A=coeffic*mdata[j][i]/(coeff2+ff2);
					B=- ( coeff2-ff2 )*Beta[j]*mdata[j][i]/(coeff2+ff2)/(coeff2+ff2);
					C=-2*coeffic*f1[j]*Beta[j]*mdata[j][i]/(coeff2+ff2)/(coeff2+ff2);
					
					float R=A*(
						(
							(hdata[j][i+1]-hdata[j][i])-(hdata[j][i]-hdata[j][i-1])
							
						)/dx2+(
							(hdata[j+1][i]-hdata[j][i])-(hdata[j][i]-hdata[j-1][i])
							
						)/dy2
					)+(
						B*(hdata[j][i+1]-hdata[j][i-1])/(dxs[j]+dxs[j])+
						C*(hdata[j+1][i]-hdata[j-1][i])/(dy   +dy   )
						
					)-coeffic*hdata[j][i]-Qdata[j][i];
					
					R*=opt_arg/(2*A/dx2+2*A/dy2-coeffic);
					
					float tmp=abs(R);
					if(tmp>maxerr) maxerr=tmp;
					
					hdata[j][i]+=R;
				}
			}
			
			if(maxerr<=err||xx++>=maxloop) break;
			
		}while(true);
		
		if(xx<maxloop) System.out.println("loops "+xx);
		else System.out.println("loops "+xx);
		
		System.out.println("Finish SORing.");
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			DataDescriptor dd=
			DescriptorFactory.getDescriptor("D:/Data/DiagnosisVortex/Catarina/catarina.ctl");
			
			SphericalSpacialModel ssm=new SphericalSpacialModel(dd);
			
			Range r=new Range("t(1,1);z(1,1)",dd);
			
			Variable Q=new Variable("olr",true,r);
			Variable h1=new Variable("h1",true,r);
			Variable h2=new Variable("h2",true,r);
			Variable h3=new Variable("h3",true,r);
			Variable h4=new Variable("h4",true,r);
			Variable h5=new Variable("h5",true,r);
			Variable h6=new Variable("h6",true,r);
			
			GlobalGillModelInSC ggm=new GlobalGillModelInSC(ssm);
			
			ggm.addHeatingSource2(Q,160,0,240,0.1f);
			
			ggm.setDissipationCoefficient(7e-5f); ggm.setBasicGeopotential(2000f); h1.setValue(2000f);
			ggm.SOREpsilon(0.1f,h1,Q); h1.setComment("coeff=7e-5, height=2000");
			
			ggm.setDissipationCoefficient(1e-5f); ggm.setBasicGeopotential(2000f); h2.setValue(2000f);
			ggm.SOREpsilon(0.1f,h2,Q); h2.setComment("coeff=1e-5, height=2000");
			
			ggm.setDissipationCoefficient(7e-6f); ggm.setBasicGeopotential(2000f); h3.setValue(2000f);
			ggm.SOREpsilon(0.1f,h3,Q); h3.setComment("coeff=7e-6, height=2000");
			
			ggm.setDissipationCoefficient(7e-5f); ggm.setBasicGeopotential(10000f); h4.setValue(10000f);
			ggm.SOREpsilon(0.1f,h4,Q); h4.setComment("coeff=7e-5, height=10000");
			
			ggm.setDissipationCoefficient(1e-5f); ggm.setBasicGeopotential(10000f); h5.setValue(10000f);
			ggm.SOREpsilon(0.1f,h5,Q); h5.setComment("coeff=1e-5, height=10000");
			
			ggm.setDissipationCoefficient(7e-6f); ggm.setBasicGeopotential(10000f); h6.setValue(10000f);
			ggm.SOREpsilon(0.1f,h6,Q); h6.setComment("coeff=7e-6, height=10000");
			
			Variable[] re1=ggm.cWindFields(h1); re1[0].setName("u1"); re1[1].setName("v1");
			Variable[] re2=ggm.cWindFields(h2); re2[0].setName("u2"); re2[1].setName("v2");
			Variable[] re3=ggm.cWindFields(h3); re3[0].setName("u3"); re3[1].setName("v3");
			Variable[] re4=ggm.cWindFields(h4); re4[0].setName("u4"); re4[1].setName("v4");
			Variable[] re5=ggm.cWindFields(h5); re5[0].setName("u5"); re5[1].setName("v5");
			Variable[] re6=ggm.cWindFields(h6); re6[0].setName("u6"); re6[1].setName("v6");
			
			miniufo.io.DataWrite dw=miniufo.io.DataIOFactory.getDataWrite(dd,"d:/Data/GillModel/GillTest/gill2.dat");
			dw.writeData(dd,Q,h1,h2,h3,h4,h5,h6,re1[0],re1[1],re2[0],re2[1],re3[0],re3[1],re4[0],re4[1],re5[0],re5[1],re6[0],re6[1]);
			dw.closeFile();
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
