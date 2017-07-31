/**
 * @(#)ThermoDynamicMethodsInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.geophysics.Empirical;
import miniufo.geophysics.ocean.SeaWater2;
import miniufo.application.EquationInSphericalCoordinate;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static miniufo.geophysics.atmos.ThermoDynamics.Cp;
import static miniufo.geophysics.atmos.ThermoDynamics.E0;
import static miniufo.geophysics.atmos.ThermoDynamics.L0;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;
import static miniufo.geophysics.atmos.ThermoDynamics.Rv;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.diagnosis.SpatialModel.GRAVITY_ACCERLERATION;


/**
 * Algorithms for thermo-dynamical methods in spherical coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class ThermoDynamicMethodsInSC extends EquationInSphericalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public ThermoDynamicMethodsInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
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
     * Calculate the maximum potential intensity of a storm according to the underlying SST.
     * This method is based on Merrill (1998, JAS) and DeMaria et al. (1993, JAS), applicable
     * to the Atlantic hurricanes.
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity in unit of surface wind speed (m s^-1)
     */
	public Variable cMPIAtl1(Variable sst){
		assignSubDomainParams(sst);
		
		Variable mpi=new Variable("mpi",sst);
		mpi.setCommentAndUnit("maximum potential intensity for Atlantic hurricanes (m s^-1)");
		mpi.setValue(undef);
		
		float[][][][] sdata=sst.getData();
		float[][][][] mdata=mpi.getData();
		
		if(sst.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(sdata[l][k][j][i]!=undef) mdata[l][k][j][i]=Empirical.cMPIAtl1(sdata[l][k][j][i]);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(sdata[k][j][i][l]!=undef) mdata[k][j][i][l]=Empirical.cMPIAtl1(sdata[k][j][i][l]);
		}
		
		return mpi;
	}
	
	/**
     * Calculate the maximum potential intensity of a storm according to the underlying SST.
     * This method is based on DeMaria et al. (1994, JC), applicable to the Atlantic hurricanes.
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity in unit of surface wind speed (m s^-1)
     */
	public Variable cMPIAtl2(Variable sst){
		assignSubDomainParams(sst);
		
		Variable mpi=new Variable("mpi",sst);
		mpi.setCommentAndUnit("maximum potential intensity over Atlantic Ocean basin (m s^-1)");
		mpi.setValue(undef);
		
		float[][][][] sdata=sst.getData();
		float[][][][] mdata=mpi.getData();
		
		if(sst.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(sdata[l][k][j][i]!=undef) mdata[l][k][j][i]=Empirical.cMPIAtl2(sdata[l][k][j][i]);
					
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(sdata[k][j][i][l]!=undef) mdata[k][j][i][l]=Empirical.cMPIAtl2(sdata[k][j][i][l]);
		}
		
		return mpi;
	}
	
	/**
     * Calculate the maximum potential intensity of a storm according to the underlying SST.
     * This method is based on Zeng (2007, MWR), applicable to TCs over western North Pacific.
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity in unit of surface wind speed (m s^-1)
     */
	public Variable cMPIWNP(Variable sst){
		assignSubDomainParams(sst);
		
		Variable mpi=new Variable("mpi",sst);
		mpi.setCommentAndUnit("maximum potential intensity for wester North Pacific typhoons (m s^-1)");
		mpi.setValue(undef);
		
		float[][][][] sdata=sst.getData();
		float[][][][] mdata=mpi.getData();
		
		if(sst.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(sdata[l][k][j][i]!=undef) mdata[l][k][j][i]=Empirical.cMPIWNP(sdata[l][k][j][i]);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(sdata[k][j][i][l]!=undef) mdata[k][j][i][l]=Empirical.cMPIWNP(sdata[k][j][i][l]);
		}
		
		return mpi;
	}
    
	
	/**
     * caculate diabatic heating rate, Cp*dT/dt+a*w=Q, units: W/kg or (m^2/s^3)
     *
     * @param	T	air temperature
     * @param	u	u-wind
     * @param	v	v-wind
     * @param	w	w-wind
     *
     * @return	Q	diabatic heating
     */
	public Variable cDiabaticHeatingRate(Variable T,Variable u,Variable v,Variable w){
		if(!T.isLike(u)||!T.isLike(v)||!T.isLike(w)) throw new IllegalArgumentException("dimensions not same");
		
		Variable Q=new Variable("Q",T);		float undef=T.getUndef();
		Q.setCommentAndUnit("diabatic heating rate: Q=Cp*dT/dt+a*omega (W kg^-1)");
		
		ystart=T.getRange().getYRange()[0];
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		float[][][][] wdata=w.getData();
		float[][][][] Tdata=T.getData();
		float[][][][] Qdata=Q.getData();
		
		if(T.isTFirst()){
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l+1][k][j][i]!=undef&&Tdata[l-1][k][j][i]!=undef)
					Qdata[l][k][j][i]=(Tdata[l+1][k][j][i]-Tdata[l-1][k][j][i])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[1  ][k][j][i]!=undef&&Tdata[0  ][k][j][i]!=undef)
					Qdata[0  ][k][j][i]=(Tdata[1  ][k][j][i]-Tdata[0  ][k][j][i])/dt;
				if(Tdata[t-1][k][j][i]!=undef&&Tdata[t-2][k][j][i]!=undef)
					Qdata[t-1][k][j][i]=(Tdata[t-1][k][j][i]-Tdata[t-2][k][j][i])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[l][k][j][i]!=undef&&Tdata[l][k][j][i+1]!=undef&&Tdata[l][k][j][i-1]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=udata[l][k][j][i]*(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
					
					else Qdata[l][k][j][i]=udata[l][k][j][i]*(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
				}
				
				if(vdata[l][k][j][i]!=undef&&Tdata[l][k][j+1][i]!=undef&&Tdata[l][k][j-1][i]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=vdata[l][k][j][i]*(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(dy*2);
					
					else Qdata[l][k][j][i]=vdata[l][k][j][i]*(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(dy*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[l][k][j][0]!=undef&&Tdata[l][k][j][1]!=undef&&Tdata[l][k][j][0]!=undef){
					if(Qdata[l][k][j][0]!=undef)
						Qdata[l][k][j][0]+=udata[l][k][j][0]*(Tdata[l][k][j][1]-Tdata[l][k][j][0])/dxs[ystart-1+j];
					
					else Qdata[l][k][j][0]=udata[l][k][j][0]*(Tdata[l][k][j][1]-Tdata[l][k][j][0])/dxs[ystart-1+j];
				}
				
				if(vdata[l][k][j][x-1]!=undef&&Tdata[l][k][j][x-1]!=undef&&Tdata[l][k][j][x-2]!=undef){
					if(Qdata[l][k][j][x-1]!=undef)
						Qdata[l][k][j][x-1]+=vdata[l][k][j][x-1]*(Tdata[l][k][j][x-1]-Tdata[l][k][j][x-2])/dxs[ystart-1+j];
					
					else Qdata[l][k][j][x-1]=vdata[l][k][j][x-1]*(Tdata[l][k][j][x-1]-Tdata[l][k][j][x-2])/dxs[ystart-1+j];
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[l][k][0][i]!=undef&&Tdata[l][k][0][i]!=undef&&Tdata[l][k][1][i]!=undef){
					if(Qdata[l][k][0][i]!=undef)
						Qdata[l][k][0][i]+=udata[l][k][0][i]*(Tdata[l][k][1][i]-Tdata[l][k][0][i])/dy;
					
					else Qdata[l][k][0][i]=udata[l][k][0][i]*(Tdata[l][k][1][i]-Tdata[l][k][0][i])/dy;
				}
				
				if(vdata[l][k][y-1][i]!=undef&&Tdata[l][k][y-1][i]!=undef&&Tdata[l][k][y-2][i]!=undef){
					if(Qdata[l][k][y-1][i]!=undef)
						Qdata[l][k][y-1][i]+=vdata[l][k][y-1][i]*(Tdata[l][k][y-1][i]-Tdata[l][k][y-2][i])/dy;
					
					else Qdata[l][k][y-1][i]=vdata[l][k][y-1][i]*(Tdata[l][k][y-1][i]-Tdata[l][k][y-2][i])/dy;
				}
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=wdata[l][k][j][i]*(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz);
					
					else Qdata[l][k][j][i]=wdata[l][k][j][i]*(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef&&Tdata[l][0][j][i]!=undef){
					if(Qdata[l][0][j][i]!=undef)
						Qdata[l][0][j][i]+=wdata[l][0][j][i]*(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
					
					else Qdata[l][0][j][i]=wdata[l][0][j][i]*(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
				}
				
				if(wdata[l][z-1][j][i]!=undef&&Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef){
					if(Qdata[l][z-1][j][i]!=undef)
						Qdata[l][z-1][j][i]+=wdata[l][z-1][j][i]*(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
					
					else Qdata[l][z-1][j][i]=wdata[l][z-1][j][i]*(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
				}
			}
			
			/*** -a*w ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&Tdata[l][k][j][i]!=undef){
					if(Qdata[l][k][j][i]!=undef){
						Qdata[l][k][j][i]*=Cp;
						Qdata[l][k][j][i]-=wdata[l][k][j][i]*Tdata[l][k][j][i]*Rd/zdef[zstart-1+k];
						
					}else Qdata[l][k][j][i]=-wdata[l][k][j][i]*Tdata[l][k][j][i]*Rd/zdef[zstart-1+k];
				}
			
		}else{
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l+1]!=undef&&Tdata[k][j][i][l-1]!=undef) Qdata[k][j][i][l]=(Tdata[k][j][i][l+1]-Tdata[k][j][i][l-1])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[k][j][i][1  ]!=undef&&Tdata[k][j][i][0  ]!=undef) Qdata[k][j][i][0  ]=(Tdata[k][j][i][1  ]-Tdata[k][j][i][0  ])/dt;
				if(Tdata[k][j][i][t-1]!=undef&&Tdata[k][j][i][t-2]!=undef) Qdata[k][j][i][t-1]=(Tdata[k][j][i][t-1]-Tdata[k][j][i][t-2])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[k][j][i][l]!=undef&&Tdata[k][j][i+1][l]!=undef&&Tdata[k][j][i-1][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=udata[k][j][i][l]*(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
					
					else Qdata[k][j][i][l]=udata[k][j][i][l]*(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
				}
				
				if(vdata[k][j][i][l]!=undef&&Tdata[k][j+1][i][l]!=undef&&Tdata[k][j-1][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=vdata[k][j][i][l]*(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(dy*2);
					
					else Qdata[k][j][i][l]=vdata[k][j][i][l]*(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(dy*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[k][j][0][l]!=undef&&Tdata[k][j][1][l]!=undef&&Tdata[k][j][0][l]!=undef){
					if(Qdata[k][j][0][l]!=undef)
						Qdata[k][j][0][l]+=udata[k][j][0][l]*(Tdata[k][j][1][l]-Tdata[k][j][0][l])/dxs[ystart-1+j];
					
					else Qdata[k][j][0][l]=udata[k][j][0][l]*(Tdata[k][j][1][l]-Tdata[k][j][0][l])/dxs[ystart-1+j];
				}
				
				if(vdata[k][j][x-1][l]!=undef&&Tdata[k][j][x-1][l]!=undef&&Tdata[k][j][x-2][l]!=undef){
					if(Qdata[k][j][x-1][l]!=undef)
						Qdata[k][j][x-1][l]+=vdata[k][j][x-1][l]*(Tdata[k][j][x-1][l]-Tdata[k][j][x-2][l])/dxs[ystart-1+j];
					
					else Qdata[k][j][x-1][l]=vdata[k][j][x-1][l]*(Tdata[k][j][x-1][l]-Tdata[k][j][x-2][l])/dxs[ystart-1+j];
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[k][0][i][l]!=undef&&Tdata[k][0][i][l]!=undef&&Tdata[k][1][i][l]!=undef){
					if(Qdata[k][0][i][l]!=undef)
						Qdata[k][0][i][l]+=udata[k][0][i][l]*(Tdata[k][1][i][l]-Tdata[k][0][i][l])/dy;
					
					else Qdata[k][0][i][l]=udata[k][0][i][l]*(Tdata[k][1][i][l]-Tdata[k][0][i][l])/dy;
				}
				
				if(vdata[k][y-1][i][l]!=undef&&Tdata[k][y-1][i][l]!=undef&&Tdata[k][y-2][i][l]!=undef){
					if(Qdata[k][y-1][i][l]!=undef)
						Qdata[k][y-1][i][l]+=vdata[k][y-1][i][l]*(Tdata[k][y-1][i][l]-Tdata[k][y-2][i][l])/dy;
					
					else Qdata[k][y-1][i][l]=vdata[k][y-1][i][l]*(Tdata[k][y-1][i][l]-Tdata[k][y-2][i][l])/dy;
				}
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=wdata[k][j][i][l]*(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz);
					
					else Qdata[k][j][i][l]=wdata[k][j][i][l]*(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef&&Tdata[0][j][i][l]!=undef){
					if(Qdata[0][j][i][l]!=undef)
						Qdata[0][j][i][l]+=wdata[0][j][i][l]*(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
					
					else Qdata[0][j][i][l]=wdata[0][j][i][l]*(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
				}
				
				if(wdata[z-1][j][i][l]!=undef&&Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef){
					if(Qdata[z-1][j][i][l]!=undef)
						Qdata[z-1][j][i][l]+=wdata[z-1][j][i][l]*(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
					
					else Qdata[z-1][j][i][l]=wdata[z-1][j][i][l]*(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
				}
			}
			
			/*** -a*w ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef){
						Qdata[k][j][i][l]*=Cp;
						Qdata[k][j][i][l]-=wdata[k][j][i][l]*Tdata[k][j][i][l]*Rd/zdef[zstart-1+k];
						
					}else Qdata[k][j][i][l]=-wdata[k][j][i][l]*Tdata[k][j][i][l]*Rd/zdef[zstart-1+k];
				}
		}
		
		return Q;
	}
	
	public Variable[] cDiabaticHeatingRateTerms(Variable T,Variable u,Variable v,Variable w){
		if(!T.isLike(u)||!T.isLike(v)||!T.isLike(w)) throw new IllegalArgumentException("dimensions not same");
		
		Variable[] Q=new Variable[4];		float undef=T.getUndef();
		Q[0]=new Variable("qloc",T); Q[0].setUndef(undef); Q[0].setCommentAndUnit("local partial term (W kg^-1)");
		Q[1]=new Variable("qadv",T); Q[1].setUndef(undef); Q[1].setCommentAndUnit("advection term (W kg^-1)");
		Q[2]=new Variable("qcon",T); Q[2].setUndef(undef); Q[2].setCommentAndUnit("convection term (W kg^-1)");
		Q[3]=new Variable("qaw" ,T); Q[3].setUndef(undef); Q[3].setCommentAndUnit("aw/Cp term (W kg^-1)");
		
		ystart=T.getRange().getYRange()[0];
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float[][][][] udata=u.getData();	float[][][][] Q0data=Q[0].getData();
		float[][][][] vdata=v.getData();	float[][][][] Q1data=Q[1].getData();
		float[][][][] wdata=w.getData();	float[][][][] Q2data=Q[2].getData();
		float[][][][] Tdata=T.getData();	float[][][][] Q3data=Q[3].getData();
		
		if(T.isTFirst()){
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l+1][k][j][i]!=undef&&Tdata[l-1][k][j][i]!=undef)
				Q0data[l][k][j][i]=(Tdata[l+1][k][j][i]-Tdata[l-1][k][j][i])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[1][k][j][i]!=undef&&Tdata[0][k][j][i]!=undef)
				Q0data[0][k][j][i]=(Tdata[1][k][j][i]-Tdata[0][k][j][i])/dt;
				
				if(Tdata[t-1][k][j][i]!=undef&&Tdata[t-2][k][j][i]!=undef)
				Q0data[t-1][k][j][i]=(Tdata[t-1][k][j][i]-Tdata[t-2][k][j][i])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[l][k][j][i]!=undef&&Tdata[l][k][j][i+1]!=undef&&Tdata[l][k][j][i-1]!=undef)
				Q1data[l][k][j][i]=udata[l][k][j][i]*(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
				
				if(vdata[l][k][j][i]!=undef&&Tdata[l][k][j+1][i]!=undef&&Tdata[l][k][j-1][i]!=undef)
				Q1data[l][k][j][i]+=vdata[l][k][j][i]*(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(dy*2);
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[l][k][j][0]!=undef&&Tdata[l][k][j][1]!=undef&&Tdata[l][k][j][0]!=undef)
				Q1data[l][k][j][0]=udata[l][k][j][0]*(Tdata[l][k][j][1]-Tdata[l][k][j][0])/dxs[ystart-1+j];
				
				if(vdata[l][k][j][x-1]!=undef&&Tdata[l][k][j][x-1]!=undef&&Tdata[l][k][j][x-2]!=undef)
				Q1data[l][k][j][x-1]=vdata[l][k][j][x-1]*(Tdata[l][k][j][x-1]-Tdata[l][k][j][x-2])/dxs[ystart-1+j];
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[l][k][0][i]!=undef&&Tdata[l][k][0][i]!=undef&&Tdata[l][k][1][i]!=undef)
				Q1data[l][k][0][i]=udata[l][k][0][i]*(Tdata[l][k][1][i]-Tdata[l][k][0][i])/dy;
				
				if(vdata[l][k][y-1][i]!=undef&&Tdata[l][k][y-1][i]!=undef&&Tdata[l][k][y-2][i]!=undef)
				Q1data[l][k][y-1][i]=vdata[l][k][y-1][i]*(Tdata[l][k][y-1][i]-Tdata[l][k][y-2][i])/dy;
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef)
				Q2data[l][k][j][i]=wdata[l][k][j][i]*(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz);
				
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef&&Tdata[l][0][j][i]!=undef)
				Q2data[l][0][j][i]=wdata[l][0][j][i]*(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
				
				if(wdata[l][z-1][j][i]!=undef&&Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef)
				Q2data[l][z-1][j][i]=wdata[l][z-1][j][i]*(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
			}
			
			/*** -a*w/Cp ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&Tdata[l][k][j][i]!=undef)
				Q3data[l][k][j][i]=-wdata[l][k][j][i]*Tdata[l][k][j][i]*Rd/zdef[zstart-1+k]/Cp;
			
		}else{
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l+1]!=undef&&Tdata[k][j][i][l-1]!=undef)
				Q0data[k][j][i][l]=(Tdata[k][j][i][l+1]-Tdata[k][j][i][l-1])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[k][j][i][1]!=undef&&Tdata[k][j][i][0]!=undef)
				Q0data[k][j][i][0]=(Tdata[k][j][i][1]-Tdata[k][j][i][0])/dt;
				
				if(Tdata[k][j][i][t-1]!=undef&&Tdata[k][j][i][t-2]!=undef)
				Q0data[k][j][i][t-1]=(Tdata[k][j][i][t-1]-Tdata[k][j][i][t-2])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[k][j][i][l]!=undef&&Tdata[k][j][i+1][l]!=undef&&Tdata[k][j][i-1][l]!=undef)
				Q1data[k][j][i][l]=udata[k][j][i][l]*(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
				
				if(vdata[k][j][i][l]!=undef&&Tdata[k][j+1][i][l]!=undef&&Tdata[k][j-1][i][l]!=undef)
				Q1data[k][j][i][l]+=vdata[k][j][i][l]*(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(dy*2);
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[k][j][0][l]!=undef&&Tdata[k][j][1][l]!=undef&&Tdata[k][j][0][l]!=undef)
				Q1data[k][j][0][l]=udata[k][j][0][l]*(Tdata[k][j][1][l]-Tdata[k][j][0][l])/dxs[ystart-1+j];
				
				if(vdata[k][j][x-1][l]!=undef&&Tdata[k][j][x-1][l]!=undef&&Tdata[k][j][x-2][l]!=undef)
				Q1data[k][j][x-1][l]=vdata[k][j][x-1][l]*(Tdata[k][j][x-1][l]-Tdata[k][j][x-2][l])/dxs[ystart-1+j];
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[k][0][i][l]!=undef&&Tdata[k][0][i][l]!=undef&&Tdata[k][1][i][l]!=undef)
				Q1data[k][0][i][l]=udata[k][0][i][l]*(Tdata[k][1][i][l]-Tdata[k][0][i][l])/dy;
				
				if(vdata[k][y-1][i][l]!=undef&&Tdata[k][y-1][i][l]!=undef&&Tdata[k][y-2][i][l]!=undef)
				Q1data[k][y-1][i][l]=vdata[k][y-1][i][l]*(Tdata[k][y-1][i][l]-Tdata[k][y-2][i][l])/dy;
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef)
				Q2data[k][j][i][l]=wdata[k][j][i][l]*(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz);
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef&&Tdata[0][j][i][l]!=undef)
				Q2data[0][j][i][l]=wdata[0][j][i][l]*(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
				
				if(wdata[z-1][j][i][l]!=undef&&Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef)
				Q2data[z-1][j][i][l]=wdata[z-1][j][i][l]*(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
			}
			
			/*** -a*w ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef)
				Q3data[k][j][i][l]=-wdata[k][j][i][l]*Tdata[k][j][i][l]*Rd/zdef[zstart-1+k]/Cp;
		}
		
		return Q;
	}
	
	
	/**
     * calculate latent heating due to condensation or evaporation,
     * L*dq/dt=Q, units: W/kg or m^2/s^3
     *
     * @param	q	specific humidity (kg/kg)
     * @param	u	uwind (m/s^-1)
     * @param	v	vwind (m/s^-1)
     * @param	w	omega (Pa/s^-1)
     * @param	T	Temperature (K)
     *
     * @return	re	latent heating
     */
	public Variable cLargeScaleLatentHeating(Variable q,Variable u,Variable v,Variable w,Variable T){
		if(!q.isLike(u)||!q.isLike(v)||!q.isLike(w)||!q.isLike(T))
			throw new IllegalArgumentException("dimensions not same");
		
		Variable re=new Variable("Hs",q);	re.setCommentAndUnit("large-scale latent heating (W kg^-1)");
		
		ystart=q.getRange().getYRange()[0];
		zstart=q.getRange().getZRange()[0];
		
		t=q.getTCount();	z=q.getZCount();	y=q.getYCount();	x=q.getXCount();
		
		float[][][][]  qdata= q.getData();
		float[][][][]  udata= u.getData();
		float[][][][]  vdata= v.getData();
		float[][][][]  wdata= w.getData();
		float[][][][]  Tdata= T.getData();
		float[][][][] redata=re.getData();
		
		float undef=q.getUndef();	re.setUndef(undef);
		
		if(q.isTFirst()){
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(qdata[l+1][k][j][i]!=undef&&qdata[l-1][k][j][i]!=undef)
					redata[l][k][j][i]=(qdata[l+1][k][j][i]-qdata[l-1][k][j][i])/(dt*2);
				else redata[l][k][j][i]=undef;
			
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[l][k][j][i]!=undef&&qdata[l][k][j][i+1]!=undef&&qdata[l][k][j][i-1]!=undef){
					if(redata[l][k][j][i]!=undef)
						redata[l][k][j][i]+=udata[l][k][j][i]*(qdata[l][k][j][i+1]-qdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
					
					else redata[l][k][j][i]=udata[l][k][j][i]*(qdata[l][k][j][i+1]-qdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
				}
				
				if(vdata[l][k][j][i]!=undef&&qdata[l][k][j+1][i]!=undef&&qdata[l][k][j-1][i]!=undef){
					if(redata[l][k][j][i]!=undef)
						redata[l][k][j][i]+=vdata[l][k][j][i]*(qdata[l][k][j+1][i]-qdata[l][k][j-1][i])/(dy*2);
					
					else redata[l][k][j][i]=vdata[l][k][j][i]*(qdata[l][k][j+1][i]-qdata[l][k][j-1][i])/(dy*2);
				}
			}
			
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&qdata[l][k+1][j][i]!=undef&&qdata[l][k-1][j][i]!=undef){
					if(redata[l][k][j][i]!=undef)
						redata[l][k][j][i]+=
							wdata[l][k][j][i]*(qdata[l][k+1][j][i]-qdata[l][k-1][j][i])/(dz+dz);
					
					else redata[l][k][j][i]=
						wdata[l][k][j][i]*(qdata[l][k+1][j][i]-qdata[l][k-1][j][i])/(dz+dz);
				}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)	// L=2500794-2370*t
				if(redata[l][k][j][i]!=undef&&Tdata[l][k][j][i]!=undef)
					redata[l][k][j][i]*=-L0+2327f*(Tdata[l][k][j][i]-273.15f);
			
		}else{
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(qdata[k][j][i][l+1]!=undef&&qdata[k][j][i][l-1]!=undef)
					redata[k][j][i][l]=(qdata[k][j][i][l+1]-qdata[k][j][i][l-1])/(dt*2);
				else redata[k][j][i][l]=undef;
			
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[k][j][i][l]!=undef&&qdata[k][j][i+1][l]!=undef&&qdata[k][j][i-1][l]!=undef){
					if(redata[k][j][i][l]!=undef)
						redata[k][j][i][l]+=udata[k][j][i][l]*(qdata[k][j][i+1][l]-qdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
					
					else redata[k][j][i][l]=udata[k][j][i][l]*(qdata[k][j][i+1][l]-qdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
				}
				
				if(vdata[k][j][i][l]!=undef&&qdata[k][j+1][i][l]!=undef&&qdata[k][j-1][i][l]!=undef){
					if(redata[k][j][i][l]!=undef)
						redata[k][j][i][l]+=vdata[k][j][i][l]*(qdata[k][j+1][i][l]-qdata[k][j-1][i][l])/(dy*2);
					
					else redata[k][j][i][l]=vdata[k][j][i][l]*(qdata[k][j+1][i][l]-qdata[k][j-1][i][l])/(dy*2);
				}
			}
			
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&qdata[k+1][j][i][l]!=undef&&qdata[k-1][j][i][l]!=undef){
					if(redata[k][j][i][l]!=undef)
						redata[k][j][i][l]+=
							wdata[k][j][i][l]*(qdata[k+1][j][i][l]-qdata[k-1][j][i][l])/(dz+dz);
					
					else redata[k][j][i][l]=
						wdata[k][j][i][l]*(qdata[k+1][j][i][l]-qdata[k-1][j][i][l])/(dz+dz);
				}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)	// L=2500794-2370*t
				if(redata[k][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef) redata[k][j][i][l]*=-L0+2327f*(Tdata[k][j][i][l]-273.15f);
		}
		
		return re;
	}
	
	/**
     * calculate latent heating, units: W/kg or m^2/s^3
     * reference: Yue and Shou, 2008, AAS
     *
     * @param	qs	saturated specific humidity
     * @param	w	omega
     * @param	T	Temperature
     *
     * @return	re	latent heating
     */
	public Variable cLargeScaleLatentHeating(Variable qs,Variable w,Variable T){
		if(!qs.isLike(w)||!qs.isLike(T)) throw new IllegalArgumentException("dimensions not same");
		
		Variable re=new Variable("Hs",qs);	re.setCommentAndUnit("large-scale latent heating (W kg^-1)");
		
		zstart=qs.getRange().getZRange()[0];
		
		t=qs.getTCount();	z=qs.getZCount();	y=qs.getYCount();	x=qs.getXCount();
		
		float[][][][]  qdata=qs.getData();
		float[][][][]  wdata= w.getData();
		float[][][][]  Tdata= T.getData();
		float[][][][] redata=re.getData();
		
		float undef=qs.getUndef();	re.setUndef(undef);
		
		if(qs.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(qdata[l][k][j][i]!=undef&&Tdata[l][k][j][i]!=undef&&wdata[l][k][j][i]!=undef){
					float tmp1=(float)(17.1543*(273.16f-36));
					float tmp2=Cp*(Tdata[l][k][j][i]-36)*(Tdata[l][k][j][i]-36);
					
					redata[l][k][j][i]=
					-(L0-2327f*(Tdata[l][k][j][i]-273.15f))*wdata[l][k][j][i]*qdata[l][k][j][i]*(
						tmp1*Rd*Tdata[l][k][j][i]*(1+0.61f*qdata[l][k][j][i])-tmp2
					)/(tmp1*(L0-2327f*(Tdata[l][k][j][i]-273.15f))*qdata[l][k][j][i]+tmp2)/zdef[zstart-1+k];
					
				}else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(qdata[k][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef&&wdata[k][j][i][l]!=undef){
					float tmp1=(float)(17.1543*(273.16f-36));
					float tmp2=Cp*(Tdata[k][j][i][l]-36)*(Tdata[k][j][i][l]-36);
					
					redata[k][j][i][l]=
					-(L0-2327f*(Tdata[k][j][i][l]-273.15f))*wdata[k][j][i][l]*qdata[k][j][i][l]*(
						tmp1*Rd*Tdata[k][j][i][l]*(1+0.61f*qdata[k][j][i][l])-tmp2
					)/(tmp1*(L0-2327f*(Tdata[k][j][i][l]-273.15f))*qdata[k][j][i][l]+tmp2)/zdef[zstart-1+k];
					
				}else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	/**
     * calculate latent heating due to cumulus convection,
     * reference: Ding Y. 1989, units: W/kg or m^2/s^3
     *
     * @param	u		u-wind
     * @param	v		v-wind
     * @param	w		vertical velocity in isobaric coordinate
     * @param	q		specific humidity
     * @param	T		Temperature
     * @param	bot		z-index of cloud bottom
     * @param	top		z-index of cloud top
     *
     * @return	re	cumulus latent heating
     */
	public Variable cCumulusLatentHeating
	(Variable u,Variable v,Variable w,Variable q,Variable T,int bot,int top){
		zstart=q.getRange().getZRange()[0];
		ystart=q.getRange().getYRange()[0];
		
		t=q.getTCount();	z=q.getZCount();
		y=q.getYCount();	x=q.getXCount();
		
		if(!q.isLike(T)||!q.isLike(u)||!q.isLike(v)||!q.isLike(w))
			throw new IllegalArgumentException("dimensions not same");
		if(z<2)
			throw new IllegalArgumentException("need a variable contains the whole air column");
		if(bot<0)
			throw new IllegalArgumentException("invalid cloud bottom index");
		if(top>z)
			throw new IllegalArgumentException("invalid cloud top index");
		if(bot>=top)
			throw new IllegalArgumentException("bottom index should be larger than top");
		
		Variable re=new Variable("Hc",q);
		re.setCommentAndUnit("cumulus-scale latent heating (W kg^-1)");
		
		float[] ts=new float[z];
		float[][][][]  udata= u.getData();
		float[][][][]  vdata= v.getData();
		float[][][][]  wdata= w.getData();
		float[][][][]  qdata= q.getData();
		float[][][][]  Tdata= T.getData();
		float[][][][] redata=re.getData();
		
		float undef=q.getUndef();	re.setUndef(undef);
		
		if(q.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				for(int k=0;k<z;k++) ts[k]=Tdata[l][k][j][i];	// copy
				
				for(int k=bot+1;k<=top;k++){	// compute Ts in cloud
					float LEs=(L0-2327f*(ts[k-1]-273.15f))*
					E0*(float)(Math.pow(273.0/ts[k-1],5.31)*Math.exp(25.22*(1-273.0/ts[k-1]))); // Es
					
					ts[k]=ts[k-1]-0.2876f*ts[k-1]/zdef[zstart-2+k]*
					(1+9.045f*LEs/zdef[zstart-2+k]/ts[k-1])*(zdef[zstart-2+k]-zdef[zstart-1+k])/
					(1+17950f*LEs/zdef[zstart-2+k]/ts[k-1]/ts[k-1])*
					(1-ts[k-1]/1300f);
				}
				
				float Tem=0;
				for(int k=bot+1;k<=top;k++)
				Tem+=(ts[k]-Tdata[l][k][j][i])*(zdef[zstart-2+k]-zdef[zstart-1+k])/2;
				
				float Vap=0;
				for(int k=bot+1;k<=top;k++)
				Vap+=(
					(qdata[l][k][j][i+1]*udata[l][k][j][i+1]-qdata[l][k][j][i-1]*udata[l][k][j][i-1])/
					(dxs[ystart-1+j]*2)+
					(qdata[l][k][j+1][i]*udata[l][k][j+1][i]-qdata[l][k][j-1][i]*udata[l][k][j-1][i])/
					(dy+dy)-
					(qdata[l][k][j  ][i]*vdata[l][k][j  ][i]*ltan[ystart-1+j]/EARTH_RADIUS)
					
				)*(zdef[zstart-2+k]-zdef[zstart-1+k])/2;
				
				Vap+=wdata[l][top][j][i]*qdata[l][top][j][i]-wdata[l][bot][j][i]*qdata[l][bot][j][i];
				
				for(int k=0;k<z;k++)
				redata[l][k][j][i]=(L0-2327f*(Tdata[l][k][j][i]-273.15f))*
				(float)Math.pow(100000f/zdef[k],Rd/Cp)*Vap/Tem*(ts[k]-Tdata[l][k][j][i]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				for(int k=0;k<z;k++) ts[k]=Tdata[k][j][i][l];	// copy
				
				for(int k=bot+1;k<=top;k++){	// compute Ts in cloud
					float LEs=(L0-2327f*(ts[k-1]-273.15f))*
					E0*(float)(Math.pow(273.0/ts[k-1],5.31)*Math.exp(25.22*(1-273.0/ts[k-1]))); // Es

					ts[k]=ts[k-1]-0.2876f*ts[k-1]/zdef[zstart-2+k]*
					(1+9.045f*LEs/zdef[zstart-2+k]/ts[k-1])*(zdef[zstart-2+k]-zdef[zstart-1+k])/
					(1+17950f*LEs/zdef[zstart-2+k]/ts[k-1]/ts[k-1])*
					(1-ts[k-1]/1300f);
				}
				
				float Tem=0;
				for(int k=bot+1;k<=top;k++)
				Tem+=(ts[k]-Tdata[k][j][i][l])*(zdef[zstart-2+k]-zdef[zstart-1+k])/2;
				
				float Vap=0;
				for(int k=bot+1;k<=top;k++)
				Vap+=(
					(qdata[k][j][i+1][l]*udata[k][j][i+1][l]-qdata[k][j][i-1][l]*udata[k][j][i-1][l])/
					(dxs[ystart-1+j]*2)+
					(qdata[k][j+1][i][l]*udata[k][j+1][i][l]-qdata[k][j-1][i][l]*udata[k][j-1][i][l])/
					(dy+dy)-
					(qdata[k][j  ][i][l]*vdata[k][j  ][i][l]*ltan[ystart-1+j]/EARTH_RADIUS)
					
				)*(zdef[zstart-2+k]-zdef[zstart-1+k])/2;
				
				Vap+=wdata[top][j][i][l]*qdata[top][j][i][l]-wdata[bot][j][i][l]*qdata[bot][j][i][l];
				
				for(int k=0;k<z;k++)
				redata[k][j][i][l]=(L0-2327f*(Tdata[k][j][i][l]-273.15f))*
				(float)Math.pow(100000f/zdef[k],Rd/Cp)*Vap/Tem*(ts[k]-Tdata[k][j][i][l]);
			}
		}
		
		return re;
	}
	
	
	/**
     * calculate (equivalent) potential vorticity
     *
     * @param	T	temperature (K)
     * @param	q	specific humidity (Kg/Kg)
     * @param	h	geopotential height (gpm)
     *
     * @return	mse		moist static energy (J/Kg)
     */
	public Variable cMoistStaticEnergy(Variable T,Variable q,Variable h){
		if(!T.isLike(q)||!T.isLike(h)) throw new IllegalArgumentException("dimensions not same");
		
		Variable mse=new Variable("mse",T);	mse.setCommentAndUnit("moist static energy (J kg^-1)");
		
		zstart=T.getRange().getZRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float[][][][] Tdata=T.getData();
		float[][][][] qdata=T.getData();
		float[][][][] hdata=h.getData();
		float[][][][] rdata=mse.getData();
		
		float undef=T.getUndef();	mse.setUndef(undef);
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(qdata[l][k][j][i]!=undef&&Tdata[l][k][j][i]!=undef&&hdata[l][k][j][i]!=undef){
				rdata[l][k][j][i]=Cp*Tdata[l][k][j][i]+GRAVITY_ACCERLERATION*hdata[l][k][j][i]+
					(L0-2327f*(Tdata[l][k][j][i]-273.15f))*qdata[l][k][j][i];
				
			}else rdata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(qdata[k][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef&&hdata[k][j][i][l]!=undef){
				rdata[k][j][i][l]=Cp*Tdata[k][j][i][l]+GRAVITY_ACCERLERATION*hdata[k][j][i][l]+
					(L0-2327f*(Tdata[k][j][i][l]-273.15f))*qdata[k][j][i][l];
				
			}else rdata[k][j][i][l]=undef;
		}
		
		return mse;
	}
	
	
	/**
     * calculate thermal wind using temperature
     *
     * @param	T	Temperature
     *
     * @return	re	thermal wind, [0] is u-component and [1] is v-component
     */
	public Variable[] cThermalWind(Variable T){
		zstart=T.getRange().getZRange()[0];
		ystart=T.getRange().getYRange()[0];
		
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		float undef=T.getUndef();
		Variable[] thm=new Variable[2];
		thm[0]=new Variable("thermu",T);	thm[0].setCommentAndUnit("zonal component of thermal wind");
		thm[1]=new Variable("thermv",T);	thm[1].setCommentAndUnit("meridinal component of thermal wind");
		
		float[][][][] tdata=     T.getData();
		float[][][][] udata=thm[0].getData();
		float[][][][] vdata=thm[1].getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++)
				for(int i=1;i<x-1;i++){
					if(tdata[l][k][j][i+1]!=undef&&tdata[l][k][j][i-1]!=undef)
						vdata[l][k][j][i]=Rd/zdef[zstart-1+k]/f1[ystart-1+j]*
						(tdata[l][k][j][i+1]-tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
						
					else vdata[l][k][j][i]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=0;i<x;i++){
					if(tdata[l][k][j+1][i]!=undef&&tdata[l][k][j-1][i]!=undef)
						udata[l][k][j][i]=-Rd/zdef[zstart-1+k]/f1[ystart-1+j]*
						(tdata[l][k][j+1][i]-tdata[l][k][j-1][i])/(dy*2);
						
					else udata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++)
				for(int i=1;i<x-1;i++){
					if(tdata[k][j][i+1][l]!=undef&&tdata[k][j][i-1][l]!=undef)
						vdata[k][j][i][l]=Rd/zdef[zstart-1+k]/f1[ystart-1+j]*
						(tdata[k][j][i+1][l]-tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
						
					else vdata[k][j][i][l]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=0;i<x;i++){
					if(tdata[k][j+1][i][l]!=undef&&tdata[k][j-1][i][l]!=undef)
						udata[k][j][i][l]=-Rd/zdef[zstart-1+k]/f1[ystart-1+j]*
						(tdata[k][j+1][i][l]-tdata[k][j-1][i][l])/(dy*2);
						
					else udata[k][j][i][l]=undef;
				}
			}
		}
		
		return thm;
	}
	
	/**
     * calculate water vapor flux
     *
     * @param	q	specific humidity
     * @param	u	uwind
     * @param	v	vwind
     *
     * @return	re	water vapor flux, [0] is x-direction, [1] is y-direction
     *
     * @exception if q u v are not dimensionally the same
     */
	public Variable[] cWaterVaporFlux(Variable q,Variable u,Variable v){
		if(!q.isLike(u)||!q.isLike(v)) throw new IllegalArgumentException("dimensions not same");
		
		Variable[] re=new Variable[2];	re[0]=new Variable("wvfx",q);	re[1]=new Variable("wvfy",q);
		
		re[0].setCommentAndUnit("water vapor flux in x-direction (m s^-1)");
		re[1].setCommentAndUnit("water vapor flux in y-direction (m s^-1)");
		
		t=q.getTCount();	z=q.getZCount();	y=q.getYCount();	x=q.getXCount();
		
		float[][][][]  qdata=null;	 qdata=    q.getData();
		float[][][][]  udata=null;	 udata=    u.getData();
		float[][][][]  vdata=null;	 vdata=    v.getData();
		float[][][][] r0data=null;	r0data=re[0].getData();
		float[][][][] r1data=null;	r1data=re[1].getData();
		
		float undef=q.getUndef();	re[0].setUndef(undef);	re[1].setUndef(undef);
		
		if(q.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(qdata[l][k][j][i]!=undef&&udata[l][k][j][i]!=undef){
					if(udata[l][k][j][i]!=undef) r0data[l][k][j][i]=qdata[l][k][j][i]*udata[l][k][j][i]/GRAVITY_ACCERLERATION;
					else r0data[l][k][j][i]=undef;
					
					if(vdata[l][k][j][i]!=undef) r1data[l][k][j][i]=qdata[l][k][j][i]*vdata[l][k][j][i]/GRAVITY_ACCERLERATION;
					else r1data[l][k][j][i]=undef;
					
				}else{
					r0data[l][k][j][i]=undef;
					r1data[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(qdata[k][j][i][l]!=undef&&udata[k][j][i][l]!=undef){
					if(udata[k][j][i][l]!=undef) r0data[k][j][i][l]=qdata[k][j][i][l]*udata[k][j][i][l]/GRAVITY_ACCERLERATION;
					else r0data[k][j][i][l]=undef;
					
					if(vdata[k][j][i][l]!=undef) r1data[k][j][i][l]=qdata[k][j][i][l]*vdata[k][j][i][l]/GRAVITY_ACCERLERATION;
					else r1data[k][j][i][l]=undef;
					
				}else{
					r0data[k][j][i][l]=undef;
					r1data[k][j][i][l]=undef;
				}
			}
		}
		
		qdata=null;	udata=null;	vdata=null;	r0data=null;	r1data=null;
		
		return re;
	}
	
	
	/**
     * Convert in-situ temperature to potential temperature referenced at 0 dbar.
     *
     * @param	T		in-situ temperature (degC)
     *
     * @return	PT		potential temperature (degC)
     */
	public Variable convertToTheta(Variable S,Variable T){
		assignSubDomainParams(T);
		
		Variable PT=new Variable("theta",T);
		PT.setCommentAndUnit("potential temperature referenced at 0 dbar (degC)");
		PT.setValue(undef);
		
		float[][][][] tdata= T.getData();
		float[][][][] sdata= S.getData();
		float[][][][] pdata=PT.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(tdata[l][k][j][i]!=undef)
				pdata[l][k][j][i]=(float)SeaWater2.theta(sdata[l][k][j][i],tdata[l][k][j][i],zdef[k],0);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(tdata[k][j][i][l]!=undef)
				pdata[k][j][i][l]=(float)SeaWater2.theta(sdata[k][j][i][l],tdata[k][j][i][l],zdef[k],0);
		}
		
		return PT;
	}
	
	/**
     * Convert potential temperature (referenced at 0 dbar) to in-situ temperature.
     *
     * @param	PT		potential temperature (degC)
     *
     * @return	T		in-situ temperature (degC)
     */
	public Variable convertToInSituT(Variable S,Variable PT){
		assignSubDomainParams(PT);
		
		Variable T=new Variable("T",PT);
		T.setCommentAndUnit("in-situ temperature (degC)");
		T.setValue(undef);
		
		float[][][][] tdata= T.getData();
		float[][][][] sdata= S.getData();
		float[][][][] pdata=PT.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(pdata[l][k][j][i]!=undef)
			tdata[l][k][j][i]=(float)SeaWater2.theta(sdata[l][k][j][i],pdata[l][k][j][i],0,zdef[k]);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(pdata[k][j][i][l]!=undef)
			tdata[k][j][i][l]=(float)SeaWater2.theta(sdata[k][j][i][l],pdata[k][j][i][l],0,zdef[k]);
		}
		
		return T;
	}
	
	
	/**
     * Compute potential density using seawater, assuming the vertical level is z (m).
     *
     * @param	S		salinity (PSU)
     * @param	PT		potential temperature (degC)
     *
     * @return	sgm		potential density sigmaT (kg/m^3)
     */
	public Variable cPotentialDensity(Variable S,Variable PT){
		assignSubDomainParams(PT);
		
		Variable sgm=new Variable("sigma",PT);
		sgm.setCommentAndUnit("potential density (kg/m^3)");
		sgm.setValue(undef);
		
		float[][][][] mdata=sgm.getData();
		float[][][][] sdata=  S.getData();
		float[][][][] tdata= PT.getData();
		
		if(sgm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(tdata[l][k][j][i]!=undef&&sdata[l][k][j][i]!=undef)
			mdata[l][k][j][i]=(float)SeaWater2.sigmat(sdata[l][k][j][i],tdata[l][k][j][i]);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(tdata[k][j][i][l]!=undef&&sdata[k][j][i][l]!=undef)
			mdata[k][j][i][l]=(float)SeaWater2.sigmat(sdata[k][j][i][l],tdata[k][j][i][l]);
		}
		
		return sgm;
	}
	
	
	/**
	 * 
	public Variable cDivergence(Variable Q,Variable T){
		if(!Q.isLike(T)) throw new IllegalArgumentException("dimensions not same");
		
		Variable div=new Variable("div",Q);
		
		zstart=T.getRange().getZRange()[0];
		
		t=Q.getTCount();	z=Q.getZCount();	y=Q.getYCount();	x=Q.getXCount();
		
		float[][][][] Qdata=Q.getData();
		float[][][][] Tdata=T.getData();
		float[][][][] ddata=div.getData();
		
		float undef=Q.getUndef(),Cv=Cp-Rd;
		
		if(Q.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Qdata[l][k-1][j][i]!=undef&&Qdata[l][k+1][j][i]!=undef&&Qdata[l][k][j][i]!=undef
				&& Tdata[l][k-1][j][i]!=undef&&Tdata[l][k+1][j][i]!=undef&&Tdata[l][k][j][i]!=undef){
					ddata[l][k][i][j]=-Qdata[l][k][j][i]/Tdata[l][k][j][i]/Cv-
						zdef[k]*(Qdata[l][k+1][j][i]-Qdata[l][k-1][j][i])/(dz[zstart-1+k]+dz[zstart-2+k])/Cv/Tdata[l][k][j][i]+
						zdef[k]*Qdata[l][k][j][i]/Tdata[l][k][j][i]/Tdata[l][k][j][i]/Cv*
						(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz[zstart-1+k]+dz[zstart-2+k]);
					
				}else ddata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Qdata[k-1][j][i][l]!=undef&&Qdata[k+1][j][i][l]!=undef&&Qdata[k][j][i][l]!=undef
				&& Tdata[k-1][j][i][l]!=undef&&Tdata[k+1][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef){
					ddata[k][i][j][l]=-Qdata[k][j][i][l]/Tdata[k][j][i][l]/Cv-
						zdef[k]*(Qdata[k+1][j][i][l]-Qdata[k-1][j][i][l])/(dz[zstart-1+k]+dz[zstart-2+k])/Cv/Tdata[k][j][i][l]+
						zdef[k]*Qdata[k][j][i][l]/Tdata[k][j][i][l]/Tdata[k][j][i][l]/Cv*
						(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz[zstart-1+k]+dz[zstart-2+k]);
					
				}else ddata[k][j][i][l]=undef;
			}
		}
		
		return div;
	}*/
	
	
	/** test
	public static void main(String[] args){
		miniufo.descriptor.CtlDescriptor ctl=new miniufo.descriptor.CtlDescriptor(
			"D:/Data/DiagnosisVortex/Durian/crt2.ctl"
		);
		
		SphericalSpacialModel ssm=new SphericalSpacialModel(ctl);
		
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
		
		miniufo.io.CtlDataReadStream cdrs=new miniufo.io.CtlDataReadStream(ctl);
	   	cdrs.readData(rh,u,v,w,T,q,td,th,pv,the,vor);	cdrs.closeFile();
	   	
	   	Variable a   =tm.cSpecificVolume(T);
	   	Variable omega=((Variable)w.clone()).multiply(-g).divide(a);omega.setName("omega");
	   	Variable  Tc =tm.cLCLTemperature(T,rh);
	   	Variable spfh=tm.cSpecificHumidity(T,rh);
	   	Variable  Td=tm.cDewPointTemperature(T,rh);	Td.setName("Td2");
	   	Variable  pt=tm.cPotentialTemperature(T);		pt.setName("pt");
	   	Variable ept=tm.cEquivalentPotentialTemperature(T,q,Tc);	ept.setName("ept");
	   	Variable epv1=tdm.cEquivalentPotentialVorticity(u,v,omega,th);	epv1.setName("epv1");
	   	Variable epv2=tdm.cEquivalentPotentialVorticity(u,v,omega,the);	epv2.setName("epv2");
	   	
	   	miniufo.io.CtlDataWriteStream cdws=new miniufo.io.CtlDataWriteStream(
	   		"D:/Data/DiagnosisVortex/Durian/test.dat"
	   	);
	   	cdws.writeData(ctl,rh,u,v,w,omega,T,q,td,th,pv,the,Tc,spfh,Td,pt,ept,epv1,epv2,vor);
	   	cdws.closeFile();
	}*/
}
