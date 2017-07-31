/**
 * @(#)ThermoMethodsInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.application.EquationInSphericalCoordinate;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static miniufo.geophysics.atmos.ThermoDynamics.*;


/**
 * calculation of thermal methods in spheral coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class ThermoMethodsInSC extends EquationInSphericalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public ThermoMethodsInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/**
     * calculate specific humidity using vapor pressure
     *
     * @param	e	vapor pressure (Pa)
     *
     * @return	q	specific humidity (Kg/Kg)
     */
	public Variable cSpecificHumidity(Variable e){
		Variable q=new Variable("spfh",e);	q.setCommentAndUnit("specific humidity or mixing ratio (Kg Kg^-1)");
		
		zstart=e.getRange().getZRange()[0];
		
		t=e.getTCount();	z=e.getZCount();	y=e.getYCount();	x=e.getXCount();
		
		float[][][][]  edata= e.getData();
		float[][][][] redata=q.getData();
		
		float undef=e.getUndef();	q.setUndef(undef);
		
		if(e.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(edata[l][k][j][i]!=undef) redata[l][k][j][i]=
					Rd/Rv*edata[l][k][j][i]/(zdef[zstart-1+k]-(1-Rd/Rv)*edata[l][k][j][i]);
					
				else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(edata[k][j][i][l]!=undef) redata[k][j][i][l]=
					Rd/Rv*edata[k][j][i][l]/(zdef[zstart-1+k]-(1-Rd/Rv)*edata[k][j][i][l]);
					
				else redata[k][j][i][l]=undef;
		}
		
		return q;
	}
	
	/**
     * calculate specific humidity using relative humidity and temperature
     *
     * @param	T	temperature (K)
     * @param	rh	relative humidity (%)
     *
     * @return	q	specific humidity (Kg/Kg)
     */
	public Variable cSpecificHumidity(Variable T,Variable rh){
		if(!T.isLike(rh)) throw new IllegalArgumentException("dimensions not same");
		
		Variable q=new Variable("spfh",T);	q.setCommentAndUnit("specific humidity / mixing ratio (Kg Kg^-1)");
		
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float[][][][]  Tdata= T.getData();
		float[][][][] rhdata=rh.getData();
		float[][][][] redata=q.getData();
		
		float undef=T.getUndef();	q.setUndef(undef);
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&rhdata[l][k][j][i]!=undef){
					float Es=(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])));
					
					redata[l][k][j][i]=Rd/Rv*Es*rhdata[l][k][j][i]/(zdef[zstart-1+k]-(1-Rd/Rv)*Es)/100;
					
				}else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&rhdata[k][j][i][l]!=undef){
					float Es=(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])));
					
					redata[k][j][i][l]=Rd/Rv*Es*rhdata[k][j][i][l]/(zdef[zstart-1+k]-(1-Rd/Rv)*Es)/100;
					
				}else redata[k][j][i][l]=undef;
		}
		
		return q;
	}
	
	/**
     * calculate saturated specific humidity using temperature
     *
     * @param	T	temperature (K)
     *
     * @return	qs	saturated specific humidity (Kg/Kg)
     */
	public Variable cSaturatedSpecificHumidity(Variable T){
		Variable q=new Variable("sspfh",T);	q.setCommentAndUnit("saturated specific humidity / mixing ratio (Kg Kg^-1)");
		
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float[][][][]  Tdata= T.getData();
		float[][][][] redata=q.getData();
		
		float undef=T.getUndef();	q.setUndef(undef);
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef){
					float Es=(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])));
					
					redata[l][k][j][i]=Rd/Rv*Es/(zdef[zstart-1+k]-(1-Rd/Rv)*Es);
					
				}else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef){
					float Es=(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])));
					
					redata[k][j][i][l]=Rd/Rv*Es/(zdef[zstart-1+k]-(1-Rd/Rv)*Es);
					
				}else redata[k][j][i][l]=undef;
		}
		
		return q;
	}
	
	
	/**
     * calculate relative humidity using specific humidity and temperature
     *
     * @param	T	temperature (K)
     * @param	q	specific humidity (Kg/Kg)
     *
     * @return	re	relative humidity (%)
     */
	public Variable cRelativeHumidity(Variable T,Variable q){
		if(!T.isLike(q)) throw new IllegalArgumentException("dimensions not same");
		
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float undef=T.getUndef();
		
		Variable re=new Variable("rh",T);
		re.setUndef(undef);	re.setCommentAndUnit("relative humidity (%)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][]  qdata= q.getData();
		float[][][][] redata=re.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&qdata[l][k][j][i]!=undef)
					redata[l][k][j][i]=Rv/Rd/
						(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])))*
						qdata[l][k][j][i]*zdef[zstart-1+k]*100;
					
				else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&qdata[k][j][i][l]!=undef)
					redata[k][j][i][l]=Rv/Rd/
						(float)(100*exp(53.67957-6743.769/Tdata[k][j][i][l]-4.8451*log(Tdata[k][j][i][l])))*
						qdata[k][j][i][l]*zdef[zstart-1+k]*100;
					
				else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
    
	
	/**
     * calculate saturated vapor pressure, Emanuel recommanded in 'Atmospheric Convection' in 1994
     *
     * @param	T	air temperature (K)
     *
     * @return	Es	saturated vapor pressure (Pa)
     */
	public Variable cSaturatedVaporPressure(Variable T){
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		Variable Es=new Variable("Es",T);	Es.setCommentAndUnit("saturate vapor pressure (Pa)");
		
		float[][][][]  Tdata=null;	 Tdata= T.getData();
		float[][][][] Esdata=null;	Esdata=Es.getData();
		
		float undef=T.getUndef();	Es.setUndef(undef);
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef)
				Esdata[l][k][j][i]=
				(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])));
					
				else Esdata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef)
				Esdata[k][j][i][l]=
				(float)(100*exp(53.67957-6743.769/Tdata[k][j][i][l]-4.8451*log(Tdata[k][j][i][l])));
					
				else Esdata[k][j][i][l]=undef;
		}
		
		return Es;
	}
	
	/**
     * calculate saturated vapor pressure, Tetens-Magnus formulas
     *
     * @param	T		air temperature (K)
     * @param	type	0 means water surface, 1 means ice surface
     *
     * @return	Es		saturated vapor pressure (Pa)
     */
	public Variable cSaturatedVaporPressureOnSurface(Variable T,int type){
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		Variable Es=new Variable("Es",T);	Es.setCommentAndUnit("Surface saturated vapor pressure (Pa)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][] Esdata=Es.getData();
		
		float undef=T.getUndef();	Es.setUndef(undef);
		
		if(T.isTFirst()){
			if(type==0){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k][j][i]!=undef)
					Esdata[l][k][j][i]=E0*(float)exp(7.5f*(Tdata[l][k][j][i]-273.15f)/(Tdata[l][k][j][i]-35.85f));
						
					else Esdata[l][k][j][i]=undef;
				}
				
			}else if(type==1){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k][j][i]!=undef)
					Esdata[l][k][j][i]=E0*(float)exp(9.5f*(Tdata[l][k][j][i]-273.15f)/(Tdata[l][k][j][i]-7.65f));
						
					else Esdata[l][k][j][i]=undef;
				}
				
			}else throw new IllegalArgumentException(
				"unknown type for calculatiion of the surface specific humidity"
			);
			
		}else{
			if(type==0){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k][j][i][l]!=undef)
					Esdata[k][j][i][l]=E0*(float)exp(7.5f*(Tdata[k][j][i][l]-273.15f)/(Tdata[k][j][i][l]-35.85f));
						
					else Esdata[k][j][i][l]=undef;
				}
				
			}else if(type==1){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k][j][i][l]!=undef)
					Esdata[k][j][i][l]=E0*(float)exp(9.5f*(Tdata[l][k][j][i]-273.15f)/(Tdata[l][k][j][i]-7.65f));
						
					else Esdata[k][j][i][l]=undef;
				}
				
			}else throw new IllegalArgumentException(
				"unknown type for calculation of the surface specific humidity"
			);
		}
		
		return Es;
	}
	
	
	/**
     * calculate specific volume according to the equation of state
     *
     * @param	T	temperature of air (K)
     *
     * @return	a	specific volume (m^3/kg)
     */
    public Variable cSpecificVolume(Variable T){
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float undef=T.getUndef();
		Variable a=new Variable("a",T);	a.setCommentAndUnit("specific volume (m^3 Kg^-1)");
		
		float[][][][] Tdata=T.getData();
		float[][][][] adata=a.getData();
		
		if(T.isTFirst()){
			for(int k=0;k<z;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k][j][i]!=undef) adata[l][k][j][i]=tmp*Tdata[l][k][j][i];
					else adata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int k=0;k<z;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k][j][i][l]!=undef) adata[k][j][i][l]=tmp*Tdata[k][j][i][l];
					else adata[k][j][i][l]=undef;
				}
			}
		}
		
		return a;
     }
	
	/**
     * calculate static stability argument
     *
     * @param	T		temperature of air (K)
     *
     * @return	sigma	static stability argument (m^4 s^2 kg^-2)
     */
    public Variable cStaticStabilityArgByT(Variable T){
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float undef=T.getUndef();
		Variable sigma=new Variable("sigma",T);	sigma.setCommentAndUnit("static stability argument (m^4 s^2 kg^-2)");
		
		float[][][][] 	  Tdata=	T.getData();
		float[][][][] sigmadata=sigma.getData();
		
		if(T.isTFirst()){
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef)
						sigmadata[l][k][j][i]=(
							tmp*Tdata[l][k][j][i]/Cp-(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz)
						)*tmp;
						
					else sigmadata[l][k][j][i]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef)
					sigmadata[l][0][j][i]=(
						tmp*Tdata[l][0][j][i]/Cp-(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz
					)*tmp;
					
				else sigmadata[l][0][j][i]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef)
					sigmadata[l][z-1][j][i]=(
						tmp*Tdata[l][z-1][j][i]/Cp-(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz
					)*tmp;
					
				else sigmadata[l][z-1][j][i]=undef;
			}
			
		}else{
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef)
						sigmadata[k][j][i][l]=(
							tmp*Tdata[k][j][i][l]/Cp-(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz)
						)*tmp;
						
					else sigmadata[k][j][i][l]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef)
					sigmadata[0][j][i][l]=(
						tmp*Tdata[0][j][i][l]/Cp-(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz
					)*tmp;
					
				else sigmadata[0][j][i][l]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef)
					sigmadata[z-1][j][i][l]=(
						tmp*Tdata[z-1][j][i][l]/Cp-(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz
					)*tmp;
					
				else sigmadata[z-1][j][i][l]=undef;
			}
		}
		
		return sigma;
     }
    
    /**
     * calculate static stability argument
     *
     * @param	th		potential temperature (K)
     *
     * @return	sigma	static stability argument(m^4 s^2 kg^-2)
     */
    public Variable cStaticStabilityArgByPT(Variable th){
		zstart=th.getRange().getZRange()[0];
		
		t=th.getTCount();	z=th.getZCount();	y=th.getYCount();	x=th.getXCount();
		
		float undef=th.getUndef();
		Variable sigma=new Variable("sigma",th);	sigma.setCommentAndUnit("static stability argument (m^4 s^2 kg^-2)");
		
		float[][][][] 	  Tdata=   th.getData();
		float[][][][] sigmadata=sigma.getData();
		
		if(th.isTFirst()){
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef)
						sigmadata[l][k][j][i]=
							-tmp*(float)pow(zdef[zstart-1+k]/100000f,Rd/Cp)*
							(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz);
						
					else sigmadata[l][k][j][i]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef)
					sigmadata[l][0][j][i]=
						-tmp*(float)pow(zdef[zstart-1]/100000f,Rd/Cp)*
						(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
					
				else sigmadata[l][0][j][i]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef)
					sigmadata[l][z-1][j][i]=
						-tmp*(float)pow(zdef[zstart-2+z]/100000f,Rd/Cp)*
						(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
					
				else sigmadata[l][z-1][j][i]=undef;
			}
			
		}else{
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef)
						sigmadata[k][j][i][l]=
							-tmp*(float)pow(zdef[zstart-1+k]/100000f,Rd/Cp)*
							(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz);
						
					else sigmadata[k][j][i][l]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef)
					sigmadata[0][j][i][l]=
						-tmp*(float)pow(zdef[zstart-1]/100000f,Rd/Cp)*
						(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
					
				else sigmadata[0][j][i][l]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef)
					sigmadata[z-1][j][i][l]=
						-tmp*(float)pow(zdef[zstart-2+z]/100000f,Rd/Cp)*
						(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
					
				else sigmadata[z-1][j][i][l]=undef;
			}
		}
		
		return sigma;
     }
	
    
	/**
     * calculate dew point temperature
     *
     * @param	e	vapor pressure (Pa)
     *
     * @return	Td	dew point temperature (K)
     */
    public Variable cDewPointTemperature(Variable e){
		zstart=e.getRange().getZRange()[0];
		
		t=e.getTCount();	z=e.getZCount();	y=e.getYCount();	x=e.getXCount();
		
		float undef=e.getUndef();
		Variable Td=new Variable("Td",e);	Td.setCommentAndUnit("dew point temperature (K)");
		
		float[][][][]  edata= e.getData();
		float[][][][] Tddata=Td.getData();
		
		if(e.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(edata[l][k][j][i]!=undef){
					float tmp=(float)log(edata[l][k][j][i]/E0);
					
					Tddata[l][k][j][i]=241.88f*tmp/(17.558f-tmp)+273.15f;
					
				}else Tddata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(edata[k][j][i][l]!=undef){
					float tmp=(float)log(edata[k][j][i][l]/E0);
					
					Tddata[k][j][i][l]=241.88f*tmp/(17.558f-tmp)+273.15f;
					
				}else Tddata[k][j][i][l]=undef;
			}
		}
		
		return Td;
	}
	
	/**
     * calculate dew point temperature
     *
     * @param	T	temperature of air (K)
     * @param	rh	relative humidity (%)
     *
     * @return	Td	dew point temperature (K)
     */
    public Variable cDewPointTemperature(Variable T,Variable rh){
		if(!T.isLike(rh)) throw new IllegalArgumentException("dimensions not same");
		
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float undef=T.getUndef();
		Variable Td=new Variable("Td",T);	Td.setCommentAndUnit("dew point temperature (K)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][] rhdata=rh.getData();
		float[][][][] Tddata=Td.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][k][j][i]!=undef&&rhdata[l][k][j][i]!=undef){
					float tmp=(float)log(exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i]))*rhdata[l][k][j][i]/E0);
					
					Tddata[l][k][j][i]=241.88f*tmp/(17.558f-tmp)+273.15f;
					
				}else Tddata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[k][j][i][l]!=undef&&rhdata[k][j][i][l]!=undef){
					float tmp=(float)log(exp(53.67957-6743.769/Tdata[k][j][i][l]-4.8451*log(Tdata[k][j][i][l]))*rhdata[k][j][i][l]/E0);
					
					Tddata[k][j][i][l]=241.88f*tmp/(17.558f-tmp)+273.15f;
					
				}else Tddata[k][j][i][l]=undef;
			}
		}
		
		return Td;
	}
	
	/**
     * calculate temperature at lifting condensation level
     * according to Bolton (1980)
     *
     * @param	T	air temperature (K)
     * @param	rh	relative humidity (%)
     *
     * @return	Tc	temperature at lifting condensation level (K)
     */
	public Variable cLCLTemperature(Variable T,Variable rh){
		if(!T.isLike(rh)) throw new IllegalArgumentException("dimensions not same");
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		Variable Tc=new Variable("Tc",T);	Tc.setCommentAndUnit("temperature at lifting condensation level (K)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][] rhdata=rh.getData();
		float[][][][] Tcdata=Tc.getData();
		
		float undef=T.getUndef();	Tc.setUndef(undef);
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&rhdata[l][k][j][i]!=undef){
					//E=(float)(exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i]))*
					//rhdata[l][k][j][i]/100);
					
					//Tcdata[l][k][j][i]=55f+2840f/(float)(3.5*log(Tdata[l][k][j][i])-log(E)-4.805);
					
					Tcdata[l][k][j][i]=55+2840/(float)(
						8.3451*log(Tdata[l][k][j][i])-53.67957+
						6743.769/Tdata[l][k][j][i]-
						log(rhdata[l][k][j][i]/100)-4.805
					);
					
				}else Tcdata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&rhdata[k][j][i][l]!=undef)
					Tcdata[k][j][i][l]=55+2840/(float)(
						8.3451*log(Tdata[k][j][i][l])-53.67957+
						6743.769/Tdata[k][j][i][l]-
						log(rhdata[k][j][i][l]/100)-4.805
					);
					
				else Tcdata[k][j][i][l]=undef;
		}
		
		return Tc;
	}
	
	/**
     * calculate equivalent temperature, Te=T+Lq/Cp
     *
     * @param	T	air temperature (K)
     * @param	q	specific humidity (Kg/Kg)
     *
     * @return	Te	equivalent temperature (K)
     */
	public Variable cEquivalentTemperature(Variable T,Variable q){
		if(!T.isLike(q)) throw new IllegalArgumentException("dimensions not same");
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		Variable re=new Variable("Te",T);	re.setCommentAndUnit("equivalent temperature (K)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][]  qdata= q.getData();
		float[][][][] redata=re.getData();
		
		float undef=T.getUndef();	re.setUndef(undef);
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&qdata[l][k][j][i]!=undef)
					redata[l][k][j][i]=Tdata[l][k][j][i]+(L0-2327f*(Tdata[l][k][j][i]-273.15f))*qdata[l][k][j][i]/Cp;
					
				else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&qdata[k][j][i][l]!=undef)
					redata[k][j][i][l]=Tdata[k][j][i][l]+(L0-2327f*(Tdata[k][j][i][l]-273.15f))*qdata[k][j][i][l]/Cp;
					
				else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	/**
     * calculate virtual temperature, Tv=(1+0.61q)T
     *
     * @param	T	air temperature (K)
     * @param	q	specific humidity (Kg/Kg)
     *
     * @return	Tv	virtual temperature (K)
     */
	public Variable cVirtualTemperature(Variable T,Variable q){
		if(!T.isLike(q)) throw new IllegalArgumentException("dimensions not same");
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		Variable re=new Variable("Tv",T);	re.setCommentAndUnit("virtual temperature (K)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][]  qdata= q.getData();
		float[][][][] redata=re.getData();
		
		float undef=T.getUndef();	re.setUndef(undef);
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&qdata[l][k][j][i]!=undef)
					redata[l][k][j][i]=Tdata[l][k][j][i]*(1+0.607717f*qdata[l][k][j][i]);
					
				else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&qdata[k][j][i][l]!=undef)
					redata[k][j][i][l]=Tdata[k][j][i][l]*(1+0.607717f*qdata[k][j][i][l]);
					
				else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	/**
     * calculate potential temperature
     *
     * @param	T		temperature of air (K)
     *
     * @return	theta	potential temperature (K)
     */
    public Variable cPotentialTemperature(Variable T){
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float undef=T.getUndef();
		Variable theta=new Variable("th",T);	theta.setCommentAndUnit("potential temperature (K)");
		
		float[][][][] 	  Tdata=	T.getData();
		float[][][][] thetadata=theta.getData();
		
		if(T.isTFirst()){
			for(int k=0;k<z;k++){
				float exner=(float)pow((100000/zdef[zstart-1+k]),(Rd/Cp));
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k][j][i]!=undef) thetadata[l][k][j][i]=exner*Tdata[l][k][j][i];
					else thetadata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int k=0;k<z;k++){
				float exner=(float)pow((100000/zdef[zstart-1+k]),(Rd/Cp));
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k][j][i][l]!=undef) thetadata[k][j][i][l]=exner*Tdata[k][j][i][l];
					else thetadata[k][j][i][l]=undef;
				}
			}
		}
		
		return theta;
	}
	
	/**
     * calculate equivalent potential temperature
     * according to Bolton (1980)
     *
     * @param	T		air temperature (K)
     * @param	q		specific humidity or mixing ratio (Kg/Kg)
     * @param	Tc		temperature at lifting condensation level (K)
     *
     * @return	thetaE	equivalent potential temperature (K)
     */
	public Variable cEquivalentPotentialTemperature(Variable T,Variable q,Variable Tc){
		if(!T.isLike(q)||!T.isLike(Tc)) throw new IllegalArgumentException("dimensions not same");
		
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		Variable re=new Variable("The",T);	re.setCommentAndUnit("equivalent potential temperature (K)");
		
		float[][][][]  qdata= q.getData();
		float[][][][]  Tdata= T.getData();
		float[][][][] Tcdata=Tc.getData();
		float[][][][] redata=re.getData();
		
		float undef=T.getUndef();	re.setUndef(undef);
		
		if(T.isTFirst()){
			for(int k=0;k<z;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Tdata[l][k][j][i]!=undef&&qdata[l][k][j][i]!=undef&&Tcdata[l][k][j][i]!=undef)
				redata[l][k][j][i]=(float)(Tdata[l][k][j][i]*
				pow(100000f/zdef[zstart-1+k],Rd/Cp*(1-0.28f*qdata[l][k][j][i]))*exp(
					(3376f/Tcdata[l][k][j][i]-2.54f)*qdata[l][k][j][i]*(1+0.81f*qdata[l][k][j][i])
				));
				
			else redata[l][k][j][i]=undef;
			
		}else{
			for(int k=0;k<z;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Tdata[k][j][i][l]!=undef&&qdata[k][j][i][l]!=undef&&Tcdata[k][j][i][l]!=undef)
				redata[k][j][i][l]=(float)(Tdata[k][j][i][l]*
				pow(100000/zdef[zstart-1+k],Rd/Cp*(1-0.28f*qdata[k][j][i][l]))*exp(
					(3376f/Tcdata[k][j][i][l]-2.54f)*qdata[k][j][i][l]*(1+0.81f*qdata[k][j][i][l])
				));
				
			else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	/**
     * calculate temperature in cumulus cloud
     * according to Ding Y., 1989
     *
     * @param	T		air temperature (K)
     * @param	bottom	start index of cloud bottom
     * @param	top		end index of cloud top
     *
     * @return	Ts		cumulus cloud temperature (K)
     */
	public Variable cCumulusCloudTemperature(Variable T,int bottom,int top){
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();
		y=T.getYCount();	x=T.getXCount();
		
		zstart=T.getRange().getZRange()[0];
		
		if(z<2)
			throw new IllegalArgumentException("need a variable contains the whole air column");
		if(bottom<0)
			throw new IllegalArgumentException("invalid cloud bottom index");
		if(top>z)
			throw new IllegalArgumentException("invalid cloud top index");
		if(bottom>=top)
			throw new IllegalArgumentException("bottom index should be larger than top");
		
		Variable re=T.copy();
		re.setName("Ts");
		re.setCommentAndUnit("cumulus cloud temperature (K)");
		
		float[][][][] tdata= T.getData();
		float[][][][] rdata=re.getData();
		
		float undef=T.getUndef();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int k=bottom+1;k<=top;k++){
				float LEs=(L0-2327f*(tdata[l][k][j][i]-273.15f))*
				(float)(100*exp(53.67957-6743.769/tdata[l][k][j][i]-4.8451*log(tdata[l][k][j][i]))); // Es
				
				if(tdata[l][k][j][i]!=undef)
				rdata[l][k][j][i]=rdata[l][k-1][j][i]-0.2876f*rdata[l][k-1][j][i]/zdef[zstart-2+k]*
				(1+9.045f*LEs/zdef[zstart-2+k]/rdata[l][k-1][j][i])*(zdef[zstart-2+k]-zdef[zstart-1+k])/
				(1+17950f*LEs/zdef[zstart-2+k]/rdata[l][k-1][j][i]/rdata[l][k-1][j][i])*
				(1-rdata[l][k-1][j][i]/1300f);
				
				else rdata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int k=bottom+1;k<=top;k++){
				float LEs=(L0-2327f*(tdata[k][j][i][l]-273.15f))*
				(float)(100*exp(53.67957-6743.769/tdata[l][k][j][i]-4.8451*log(tdata[l][k][j][i]))); // Es
				
				if(tdata[k][j][i][l]!=undef)
				rdata[k][j][i][l]=rdata[k-1][j][i][l]-0.2876f*rdata[k-1][j][i][l]/zdef[zstart-2+k]*
				(1+9.045f*LEs/zdef[zstart-2+k]/rdata[k-1][j][i][l])*(zdef[zstart-2+k]-zdef[zstart-1+k])/
				(1+17950f*LEs/zdef[zstart-2+k]/rdata[k-1][j][i][l]/rdata[k-1][j][i][l])*
				(1-rdata[k-1][j][i][l]/1300f);
				
				else rdata[k][j][i][l]=undef;
			}
		}
		
		return re;
	}
	
	
	/**
     * calculate Exner function, pi=Cp*(p/p0)^(Rd/Cp)=Cp*T/theta
     * 
     * @param	T	any given variable for newing the result (K)
     *
     * @return	Exner function
     */
	public Variable cExnerFunction(Variable T){
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		Variable E=new Variable("exner",T);	E.setCommentAndUnit("Exner function, pi=Cp*(p/p0)^(Rd/Cp)=T/theta");
		
		float[][][][] Edata=E.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Edata[l][k][j][i]=Cp*(float)pow(zdef[k]/100000,Rd/Cp);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Edata[l][k][j][i]=Cp*(float)pow(zdef[k]/100000,Rd/Cp);
		}
		
		return E;
	}
	
	
	/**
     * calculate the maximum intensity of a storm according to
     * the underlying SST (DeMaria et al. 1993, JAS)
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity (m/s)
     */
	public Variable cMaxPotentialIntensity(Variable sst){
		zstart=sst.getRange().getZRange()[0];
		
		t=sst.getTCount();	z=sst.getZCount();	y=sst.getYCount();	x=sst.getXCount();
		
		float undef=sst.getUndef();
		Variable mpi=new Variable("mpi",sst);	mpi.setCommentAndUnit("maximum potential intensity (m s^-1)");
		
		float[][][][] sdata=sst.getData();
		float[][][][] mdata=mpi.getData();
		
		if(sst.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(sdata[l][k][j][i]!=undef){
					mdata[l][k][j][i]=(float)(0.51444*74.0*exp(0.2*(sdata[l][k][j][i]-273.15-25.0)));
					
				}else mdata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(sdata[k][j][i][l]!=undef){
					mdata[l][k][j][i]=(float)(0.51444*74.0*exp(0.2*(sdata[k][j][i][l]-273.15-25.0)));
					
				}else mdata[k][j][i][l]=undef;
			}
		}
		
		return mpi;
	}
	
	
	/** test
	public static void main(String[] args){
		miniufo.descriptor.CtlDescriptor ctl=new miniufo.descriptor.CtlDescriptor(
			"D:/Data/DiagnosisVortex/Durian/crt2.ctl"
		);
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(ctl);
		
		ThermoMethodsInSC tm=new ThermoMethodsInSC(ssm);
		ThermoDynamicMethodsInSC tdm=new ThermoDynamicMethodsInSC(ssm);
		
		miniufo.diagnosis.Range range=new miniufo.diagnosis.Range("t(1,1)",ctl);
		
		Variable rh=new Variable("rh",true,range);
		Variable u =new Variable("u" ,true,range);
		Variable v =new Variable("v" ,true,range);
		Variable w =new Variable("w" ,true,range);
		Variable T =new Variable("t" ,true,range);
		Variable q =new Variable("q" ,true,range);
		Variable td=new Variable("td",true,range);
		Variable th=new Variable("th",true,range);
		Variable pv=new Variable("pv",true,range);
		Variable the=new Variable("the",true,range);
		Variable vor=new Variable("vor",true,range);
		Variable hgt=new Variable("h",true,range);
		
		miniufo.io.CtlDataReadStream cdrs=new miniufo.io.CtlDataReadStream(ctl);
    	cdrs.readData(rh,u,v,w,T,q,td,th,pv,the,vor,hgt);	cdrs.closeFile();
    	
    	Variable a   =tm.cSpecificVolume(T);
    	Variable omega=((Variable)w.copy()).multiply(-GRAVITY_ACCERLERATION).divide(a);omega.setName("omega");
    	Variable  Tc =tm.cLCLTemperature(T,rh);
    	Variable spfh=tm.cSpecificHumidity(T,rh);
    	Variable  Td=tm.cDewPointTemperature(T,rh);	Td.setName("Td2");
    	Variable  pt=tm.cPotentialTemperature(T);		pt.setName("pt");
    	Variable ept=tm.cEquivalentPotentialTemperature(T,q,Tc);	ept.setName("ept");
    	Variable Te =tm.cEquivalentTemperature(T,q);
    	Variable Tv =tm.cVirtualTemperature(T,q);
    	Variable ept2=tm.cPotentialTemperature(Te);	ept2.setName("ept2");
    	Variable epv=tdm.cEquivalentPotentialVorticity(u,v,omega,th);	epv.setName("epv");
    	
    	miniufo.io.CtlDataWriteStream cdws=new miniufo.io.CtlDataWriteStream(
    		"D:/Data/DiagnosisVortex/Durian/test.dat"
    	);
    	cdws.writeData(ctl,rh,u,v,w,omega,T,Te,q,td,th,pv,the,Tc,Tv,spfh,Td,pt,ept,ept2,epv,vor,hgt);
    	cdws.closeFile();
	}*/
}
