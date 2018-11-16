/**
 * @(#)VelocityFieldInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.application.EquationInSphericalCoordinate;
import miniufo.application.advanced.EllipticEqSORSolver2D.DimCombination;
import static java.lang.Math.abs;
import static miniufo.diagnosis.SpatialModel.REarth;
import static miniufo.diagnosis.SpatialModel.gEarth;


/**
 * Analysis methods of velocity fields in spherical coordinates.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class VelocityFieldInSC extends EquationInSphericalCoordinate{
	//
	protected int threadCount=1;
	
	
	/**
     * Constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public VelocityFieldInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/*** getor and setor ***/
	public void setThreadCount(int count){
		if(count<1)
			throw new IllegalArgumentException("Thread count should be positive");
		
		threadCount=count;
	}
	
	
	/**
     * Calculate geostrophic velocity, taking the singularity near the equator into account.
     * Reference: Lagerloef et al. 1999, JGR
     *
     * @param	h	geopotential height or sea surface height (m)
     *
     * @return	gv	geostrophic velocity (m s^-1), [0] is zonal component and [1] meridional component
     */
	public Variable[] cGeostrophicVelocity(Variable h){
		assignSubDomainParams(h);
		
		final double faiRf=Math.toRadians(2.2f);	// latitudinal length scale
		final double Campl=0.7f;					// constant amplitude scaling
		
		Variable[] gw=new Variable[2];
		gw[0]=new Variable("Ug",h); gw[0].setValue(undef);
		gw[1]=new Variable("Vg",h); gw[1].setValue(undef);
		gw[0].setCommentAndUnit("geostrophic velocity in x-direction (m s^-1)");
		gw[1].setCommentAndUnit("geostrophic velocity in y-direction (m s^-1)");
		
		float[][][][]  Ugdata=gw[0].getData();
		float[][][][]  Vgdata=gw[1].getData();
		float[][][][] hgtdata=    h.getData();
		
		if(h.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++){
				if(Math.abs(ydef[ystart-1+j])<=Math.toRadians(5)){	// within equatorial band (5S - 5N)
					double tmp=Math.exp(-Math.pow(ydef[ystart-1+j]/faiRf,2.0));
					
					// weights
					float Wb=(float)(Campl*tmp);
					float Wf=(float)(1-tmp);
					
					for(int i=1,I=x-1;i<I;i++){
						if(hgtdata[l][k][j+1][i]!=undef&&hgtdata[l][k][j-1][i]!=undef&&hgtdata[l][k][j][i]!=undef){
							float uf=-(hgtdata[l][k][j+1][i]-hgtdata[l][k][j-1][i])
									 *gEarth/(dy*2*f1[ystart-1+j]);
							float ub=-(hgtdata[l][k][j+1][i]+hgtdata[l][k][j-1][i]-2*hgtdata[l][k][j][i])
									 *gEarth/(dy*dy*Beta[ystart-1+j]);
							
							uf=Float.isInfinite(uf)?0:uf;
							Ugdata[l][k][j][i]=uf*Wf+ub*Wb;
						}
						
						if(hgtdata[l][k][j+1][i+1]!=undef&&hgtdata[l][k][j+1][i-1]!=undef&&hgtdata[l][k][j-1][i+1]!=undef&&hgtdata[l][k][j-1][i-1]!=undef){
							float vf=(hgtdata[l][k][j][i+1]-hgtdata[l][k][j][i-1])
									*gEarth/(dxs[ystart-1+j]*2*f1[ystart-1+j]);
							float vb=((hgtdata[l][k][j+1][i+1]-hgtdata[l][k][j+1][i-1])-(hgtdata[l][k][j-1][i+1]-hgtdata[l][k][j-1][i-1]))
									*gEarth/(dxs[ystart-1+j]*2*dy*2*Beta[ystart-1+j]);
							
							vf=Float.isInfinite(vf)?0:vf;
							Vgdata[l][k][j][i]=vf*Wf+vb*Wb;
						}
					}
					
				}else for(int i=1,I=x-1;i<I;i++){ // outside equatorial band
					if(hgtdata[l][k][j+1][i]!=undef&&hgtdata[l][k][j-1][i]!=undef)
					Ugdata[l][k][j][i]=-(hgtdata[l][k][j+1][i]-hgtdata[l][k][j-1][i])*gEarth/(dy*2*f1[ystart-1+j]);
					
					if(hgtdata[l][k][j][i+1]!=undef&&hgtdata[l][k][j][i-1]!=undef)
					Vgdata[l][k][j][i]=(hgtdata[l][k][j][i+1]-hgtdata[l][k][j][i-1])*gEarth/(dxs[ystart-1+j]*2*f1[ystart-1+j]);
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++){
				if(Math.abs(ydef[ystart-1+j])<=Math.toRadians(5)){	// within equatorial band (5S - 5N)
					double tmp=Math.exp(-Math.pow(ydef[ystart-1+j]/faiRf,2.0));
					
					// weights
					float Wb=(float)(Campl*tmp);
					float Wf=(float)(1-tmp);
					
					for(int i=1,I=x-1;i<I;i++){
						if(hgtdata[k][j+1][i][l]!=undef&&hgtdata[k][j-1][i][l]!=undef&&hgtdata[k][j][i][l]!=undef){
							float uf=-(hgtdata[k][j+1][i][l]-hgtdata[k][j-1][i][l])
									 *gEarth/(dy*2*f1[ystart-1+j]);
							float ub=-(hgtdata[k][j+1][i][l]+hgtdata[k][j-1][i][l]-2*hgtdata[k][j][i][l])
									 *gEarth/(dy*dy*Beta[ystart-1+j]);
							
							uf=Float.isInfinite(uf)?0:uf;
							Ugdata[k][j][i][l]=uf*Wf+ub*Wb;
						}
						
						if(hgtdata[k][j+1][i+1][l]!=undef&&hgtdata[k][j+1][i-1][l]!=undef&&hgtdata[k][j-1][i+1][l]!=undef&&hgtdata[k][j-1][i-1][l]!=undef){
							float vf=(hgtdata[k][j][i+1][l]-hgtdata[k][j][i-1][l])
									*gEarth/(dxs[ystart-1+j]*2*f1[ystart-1+j]);
							float vb=((hgtdata[k][j+1][i+1][l]-hgtdata[k][j+1][i-1][l])-(hgtdata[k][j-1][i+1][l]-hgtdata[k][j-1][i-1][l]))
									*gEarth/(dxs[ystart-1+j]*2*dy*2*Beta[ystart-1+j]);
							
							vf=Float.isInfinite(vf)?0:vf;
							Vgdata[k][j][i][l]=vf*Wf+vb*Wb;
						}
					}
					
				}else for(int i=1,I=x-1;i<I;i++){ // outside equatorial band
					if(hgtdata[k][j+1][i][l]!=undef&&hgtdata[k][j-1][i][l]!=undef)
					Ugdata[k][j][i][l]=-(hgtdata[k][j+1][i][l]-hgtdata[k][j-1][i][l])*gEarth/(dy*2*f1[ystart-1+j]);
					
					if(hgtdata[k][j][i+1][l]!=undef&&hgtdata[k][j][i-1][l]!=undef)
					Vgdata[k][j][i][l]=(hgtdata[k][j][i+1][l]-hgtdata[k][j][i-1][l])*gEarth/(dxs[ystart-1+j]*2*f1[ystart-1+j]);
				}
			}
		}
		
		return gw;
	}
	
	/**
     * Calculate rotational velocity using streamfunction.
     * 
     * @param	sf		stream function (m^2 s^-1)
     *
     * @return	rot		rotational velocity, [0] is in x direction while [1] is in y direction
     */
	public Variable[] cRotationalVelocity(Variable sf){
		assignSubDomainParams(sf);
		
		Variable[] rot=new Variable[2];
		rot[0]=new Variable("Usf",sf);	rot[0].setValue(undef);
		rot[1]=new Variable("Vsf",sf);	rot[1].setValue(undef);
		rot[0].setCommentAndUnit("rotational velocity in x-direction (m s^-1)");
		rot[1].setCommentAndUnit("rotational velocity in y-direction (m s^-1)");
		
		float[][][][] Usdata=rot[0].getData();
		float[][][][] Vsdata=rot[1].getData();
		float[][][][] sfdata=     sf.getData();
		
		if(sf.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			if(sfdata[l][k][j+1][i]!=undef&&sfdata[l][k][j-1][i]!=undef)
			Usdata[l][k][j][i]=-(sfdata[l][k][j+1][i]-sfdata[l][k][j-1][i])/(dy*2);
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=1;i<x-1;i++)
			if(sfdata[l][k][j][i+1]!=undef&&sfdata[l][k][j][i-1]!=undef)
			Vsdata[l][k][j][i]=(sfdata[l][k][j][i+1]-sfdata[l][k][j][i-1])/(dxs[ystart-1+j]+dxs[ystart-1+j]);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			if(sfdata[k][j+1][i][l]!=undef&&sfdata[k][j-1][i][l]!=undef)
			Usdata[k][j][i][l]=-(sfdata[k][j+1][i][l]-sfdata[k][j-1][i][l])/(dy*2);
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=1;i<x-1;i++)
			if(sfdata[k][j][i+1][l]!=undef&&sfdata[k][j][i-1][l]!=undef)
			Vsdata[k][j][i][l]=(sfdata[k][j][i+1][l]-sfdata[k][j][i-1][l])/(dxs[ystart-1+j]+dxs[ystart-1+j]);
		}
		
		return rot;
	}
	
	
	/**
     * Calculate Ekman current (Niiler et al. 2003, GRL; Ralph and Niller 1999, JPO).
     *
     * @param	u10		u-wind at 10 meter (m s^-1)
     * @param	v10		v-wind at 10 meter (m s^-1)
     *
     * @return	ec		Ekman velocity at 15m depth (m s^-1)
     */
	public Variable[] cEkmanVelocityByWind(Variable u10,Variable v10){
		checkDimensions(u10,v10);
		assignSubDomainParams(u10);
		
		Variable[] ec=new Variable[2];
		ec[0]=new Variable("Uek",u10);	ec[0].setValue(undef);
		ec[1]=new Variable("Vek",u10);	ec[1].setValue(undef);
		ec[0].setCommentAndUnit("zonal component of Ekman velocity at 15m depth (m s^-1)");
		ec[1].setCommentAndUnit("meridional component of Ekman velocity at 15m depth (m s^-1)");
		
		final float A=7e-7f;
		final double theta=Math.toRadians(54);
		
		float[][][][]  udata=u10.getData();
		float[][][][]  vdata=v10.getData();
		float[][][][] eudata=ec[0].getData();
		float[][][][] evdata=ec[1].getData();
		
		if(u10.isTFirst()){
			for(int j=0;j<y;j++){
				double the=ydef[ystart-1+j]<=0?theta:-theta;
				float coef=A/(float)Math.sqrt(Math.abs(f1[ystart-1+j]));
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++)
				if(udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef){
					eudata[l][k][j][i]=coef*
					(float)(Math.cos(the)*udata[l][k][j][i]-Math.sin(the)*vdata[l][k][j][i]);
					evdata[l][k][j][i]=coef*
					(float)(Math.cos(the)*vdata[l][k][j][i]+Math.sin(the)*udata[l][k][j][i]);
				}
			}
			
		}else{
			for(int j=0;j<y;j++){
				double the=ydef[ystart-1+j]<=0?theta:-theta;
				float coef=A/(float)Math.sqrt(Math.abs(f1[ystart-1+j]));
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++)
				if(udata[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef){
					eudata[k][j][i][l]=coef*
					(float)(Math.cos(the)*udata[k][j][i][l]-Math.sin(the)*vdata[k][j][i][l]);
					evdata[k][j][i][l]=coef*
					(float)(Math.cos(the)*vdata[k][j][i][l]+Math.sin(the)*udata[k][j][i][l]);
				}
			}
		}
		
		return ec;
	}
	
	/** unfinished... need validation
     * Calculate Ekman current, taking into account the equatorial region (25S - 25N).
     * Reference:
     *   van Meurs and Niiler 1997, JPO;
     *   Lagerloef et al. 1999, JGR;
     *   Sallee et al. 2008, JMR;
     *   Sudre and Morrow 2008, OD
     *
     * @param	taux	zonal component of wind stress (m s^-1)
     * @param	tauy	meridional component of wind stress (m s^-1)
     *
     * @return	ec		Ekman velocity at 15m depth (m s^-1)
     */
	public Variable[] cEkmanVelocityByStress1(Variable taux,Variable tauy){
 		checkDimensions(taux,tauy);
 		assignSubDomainParams(taux);
 		
 		Variable[] ec=new Variable[2];
 		ec[0]=new Variable("Uek",taux);	ec[0].setValue(undef);
 		ec[1]=new Variable("Vek",taux);	ec[1].setValue(undef);
 		ec[0].setCommentAndUnit("zonal component of Ekman velocity at 15m depth (m s^-1)");
 		ec[1].setCommentAndUnit("meridional component of Ekman velocity at 15m depth (m s^-1)");
 		
 		final float r=2.15e-4f;
 		final float hmd=32.5f;
 		final float row=1027f;
 		
 		float[][][][]  udata=taux.getData();
 		float[][][][]  vdata=tauy.getData();
 		float[][][][] eudata=ec[0].getData();
 		float[][][][] evdata=ec[1].getData();
 		
 		if(taux.isTFirst()){
 			for(int j=0;j<y;j++)
 			if(Math.abs(ydef[ystart-1+j])<=Math.toRadians(43)){ // within equatorial band (25S - 25N)
 				double theta=-Math.atan(f1[ystart-1+j]*hmd/r);
	 			
	 			float B=1f/(float)Math.sqrt((r*r+hmd*hmd*f1[ystart-1+j]*f1[ystart-1+j]))/row;
	 			float sint=(float)Math.sin(theta);
	 			float cost=(float)Math.cos(theta);
	 			
	 			for(int l=0;l<t;l++)
	 			for(int k=0;k<z;k++)
	 			for(int i=0;i<x;i++)
	 			if(udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef){
	 				eudata[l][k][j][i]=B*(cost*udata[l][k][j][i]-sint*vdata[l][k][j][i]);
	 				evdata[l][k][j][i]=B*(cost*vdata[l][k][j][i]+sint*udata[l][k][j][i]);
	 			}
	 			
 			}else{ // outside equatorial band
 				double theta=Math.toRadians(ydef[ystart-1+j]<=0?55:-55);
	 			
	 			float B=0.3f;
	 			float sint=(float)Math.sin(theta);
	 			float cost=(float)Math.cos(theta);
	 			
	 			for(int l=0;l<t;l++)
	 			for(int k=0;k<z;k++)
	 			for(int i=0;i<x;i++)
	 			if(udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef){
	 				eudata[l][k][j][i]=B*(cost*udata[l][k][j][i]-sint*vdata[l][k][j][i]);
	 				evdata[l][k][j][i]=B*(cost*vdata[l][k][j][i]+sint*udata[l][k][j][i]);
	 			}
 			}
 			
 		}else{
 			for(int j=0;j<y;j++)
 			if(Math.abs(ydef[ystart-1+j])<=Math.toRadians(43)){ // within equatorial band (25S - 25N)
 				double theta=-Math.atan(f1[ystart-1+j]*hmd/r);
	 			
	 			float B=1f/(float)Math.sqrt((r*r+hmd*hmd*f1[ystart-1+j]*f1[ystart-1+j]))/row;
	 			float sint=(float)Math.sin(theta);
	 			float cost=(float)Math.cos(theta);
	 			
	 			for(int l=0;l<t;l++)
	 			for(int k=0;k<z;k++)
	 			for(int i=0;i<x;i++)
	 			if(udata[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef){
	 				eudata[k][j][i][l]=B*(cost*udata[k][j][i][l]-sint*vdata[k][j][i][l]);
	 				evdata[k][j][i][l]=B*(cost*vdata[k][j][i][l]+sint*udata[k][j][i][l]);
	 			}
	 			
 			}else{ // outside equatorial band
 				double theta=Math.toRadians(ydef[ystart-1+j]<=0?55:-55);
	 			
	 			float B=0.3f;
	 			float sint=(float)Math.sin(theta);
	 			float cost=(float)Math.cos(theta);
	 			
	 			for(int l=0;l<t;l++)
	 			for(int k=0;k<z;k++)
	 			for(int i=0;i<x;i++)
	 			if(udata[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef){
	 				eudata[k][j][i][l]=B*(cost*udata[k][j][i][l]-sint*vdata[k][j][i][l]);
	 				evdata[k][j][i][l]=B*(cost*vdata[k][j][i][l]+sint*udata[k][j][i][l]);
	 			}
 			}
 		}
 		
 		return ec;
 	}
	
	/**
     * Calculate Ekman current (Rypina et al. 2012, JPO; Ralph and Niiler 1999, JPO).
     *
     * @param	taux	zonal component of wind stress (N m^-2)
     * @param	tauy	meridional component of wind stress (N m^-2)
     *
     * @return	ec		Ekman velocity at 15m depth (m s^-1)
     */
	public Variable[] cEkmanVelocityByStress2(Variable taux,Variable tauy){
 		checkDimensions(taux,tauy);
 		assignSubDomainParams(taux);
 		
 		Variable[] ec=new Variable[2];
 		ec[0]=new Variable("Uek",taux);	ec[0].setValue(undef);
 		ec[1]=new Variable("Vek",taux);	ec[1].setValue(undef);
 		ec[0].setCommentAndUnit("zonal component of Ekman velocity at 15m depth (m s^-1)");
 		ec[1].setCommentAndUnit("meridional component of Ekman velocity at 15m depth (m s^-1)");
 		
 		final double beta=0.065f;
 		final double row=1027f;
 		final double theta=Math.toRadians(55);
 		
 		float[][][][]  udata=taux.getData();
 		float[][][][]  vdata=tauy.getData();
 		float[][][][] eudata=ec[0].getData();
 		float[][][][] evdata=ec[1].getData();
 		
 		if(taux.isTFirst()){
 			for(int j=0;j<y;j++){
 				double the=ydef[ystart-1+j]<=0?theta:-theta;
 				
 				float coef=(float)(beta/Math.sqrt(Math.abs(f1[ystart-1+j])*row));
 				float sint=(float)Math.sin(the);
 				float cost=(float)Math.cos(the);
 				
 				for(int l=0;l<t;l++)
 				for(int k=0;k<z;k++)
 				for(int i=0;i<x;i++)
 				if(udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef){
 					eudata[l][k][j][i]=coef/(float)Math.sqrt(Math.hypot(udata[l][k][j][i],vdata[l][k][j][i]))*
 					(cost*udata[l][k][j][i]-sint*vdata[l][k][j][i]);
 					evdata[l][k][j][i]=coef/(float)Math.sqrt(Math.hypot(udata[l][k][j][i],vdata[l][k][j][i]))*
 					(cost*vdata[l][k][j][i]+sint*udata[l][k][j][i]);
 				}
 			}
 			
 		}else{
 			for(int j=0;j<y;j++){
 				double the=ydef[ystart-1+j]<=0?theta:-theta;
 				
 				float coef=(float)(beta/Math.sqrt(Math.abs(f1[ystart-1+j])*row));
 				float sint=(float)Math.sin(the);
 				float cost=(float)Math.cos(the);
 				
 				for(int l=0;l<t;l++)
 				for(int k=0;k<z;k++)
 				for(int i=0;i<x;i++)
 				if(udata[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef){
 					eudata[k][j][i][l]=coef/(float)Math.sqrt(Math.hypot(udata[k][j][i][l],vdata[k][j][i][l]))*
 					(cost*udata[k][j][i][l]-sint*vdata[k][j][i][l]);
 					evdata[k][j][i][l]=coef/(float)Math.sqrt(Math.hypot(udata[k][j][i][l],vdata[k][j][i][l]))*
 					(cost*vdata[k][j][i][l]+sint*udata[k][j][i][l]);
 				}
 			}
 		}
 		
 		return ec;
 	}
	
	
	/**
     * Endlich's method for calculating rotational velocity.
     *
     * @param	uu	zonal velocity (m s^-1)
     * @param	vv	meridional velocity (m s^-1)
     *
     * @return	rw	rotational velocity, [0] is x component and [1] is y component
     */
	public Variable[] cRotationalVelocityByEndlich(Variable uu,Variable vv){
		System.out.println("\nstart calculating rotational velocity...");
		
 		checkDimensions(uu,vv);
 		assignSubDomainParams(uu);
		
		Variable[] rw=new Variable[2];
		rw[0]=uu.copy(); rw[0].setName("us");
		rw[1]=vv.copy(); rw[1].setName("vs");
		
		rw[0].setCommentAndUnit("x-component of rotational velocity (m s^-1)");
		rw[1].setCommentAndUnit("y-component of rotational velocity (m s^-1)");
		
		float[][]    vordata=new float[y-2][x-2];
		float[][][][]  udata=rw[0].getData();
		float[][][][]  vdata=rw[1].getData();
		float[][][][] uudata=uu.getData();
		float[][][][] vvdata=vv.getData();
		
		if(uu.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0,cu=0,cv=0;
				float max=1,tmpmax=0,aveu0=0,avev0=0,aveu1=0,avev1=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef){ aveu0+=udata[l][k][j][i]; cu++;}
					if(vdata[l][k][j][i]!=undef){ avev0+=vdata[l][k][j][i]; cv++;}
				}
				
				if(cu!=0) aveu0/=cu; cu=0;		if(cv!=0) avev0/=cv; cv=0;
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(vvdata[l][k][j][i+1]!=undef&&vvdata[l][k][j][i-1]!=undef
					 &&uudata[l][k][j+1][i]!=undef&&uudata[l][k][j-1][i]!=undef)
					vordata[j-1][i-1]=
						(vvdata[l][k][j][i+1]-vvdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)-
						(uudata[l][k][j+1][i]-uudata[l][k][j-1][i])/(dy*2)+
						(uudata[l][k][j][i]/REarth)*ltan[ystart-1+j];
					
					else vordata[j-1][i-1]=undef;
				}
				
				do{
					for(int j=1;j<y-1;j++)
					for(int i=1;i<x-1;i++){
						if(vordata[j-1][i-1]!=undef
						 &&udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef){
							// calculate divergence
							float tmpdiv=
							(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
							(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
							(vdata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the divergence
							// multiplied by lcos to ensure the convergence of iteration
							float tmpu=-dxs[ystart-1+j]*tmpdiv/2*lcos[j];
							float tmpv=-dy*tmpdiv/2*lcos[j];
							
							// correct the u-velocity and v-velocity
							udata[l][k][j][i+1]+=tmpu;
							udata[l][k][j][i-1]-=tmpu;
							vdata[l][k][j+1][i]+=tmpv;
							vdata[l][k][j-1][i]-=tmpv;
							
							// calculate vorticity
							float tmpvor=
							(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)-
							(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)+
							(udata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the new vorticity
							// multiplied by lcos to ensure the convergence of iteration
							tmpu=dy*(tmpvor-vordata[j-1][i-1])/2*lcos[j];
							tmpv=-dxs[ystart-1+j]*(tmpvor-vordata[j-1][i-1])/2*lcos[j];
							
							// correct the u-velocity and v-velocity
							udata[l][k][j+1][i]+=tmpu;
							udata[l][k][j-1][i]-=tmpu;
							vdata[l][k][j][i+1]+=tmpv;
							vdata[l][k][j][i-1]-=tmpv;
							
							// calculate the new divergence
							tmpdiv=
							(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
							(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
							(vdata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							float tmp=abs(tmpdiv);
							if(tmp>tmpmax) tmpmax=tmp;
						}
					}
					
					if(max<1e-12f) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<9999);
				
				System.out.println(String.format("Max error of divergence is %.4e Adjust for %4d loops",max,count));
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef){ aveu1+=udata[l][k][j][i]; cu++;}
					if(vdata[l][k][j][i]!=undef){ avev1+=vdata[l][k][j][i]; cv++;}
				}
				
				if(cu!=0) aveu1/=cu;	if(cv!=0) avev1/=cv;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef) udata[l][k][j][i]+=aveu0-aveu1;
					if(vdata[l][k][j][i]!=undef) vdata[l][k][j][i]+=avev0-avev1;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(vordata[j-1][i-1]==undef){
					udata[l][k][j][i]=undef;
					vdata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0,cu=0,cv=0;
				float max=1,tmpmax=0,aveu0=0,avev0=0,aveu1=0,avev1=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[k][j][i][l]!=undef){ aveu0+=udata[k][j][i][l]; cu++;}
					if(vdata[k][j][i][l]!=undef){ avev0+=vdata[k][j][i][l]; cv++;}
				}
				
				if(cu!=0) aveu0/=cu; cu=0;		if(cv!=0) avev0/=cv; cv=0;
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(vvdata[k][j][i+1][l]!=undef&&vvdata[k][j][i-1][l]!=undef
					 &&uudata[k][j+1][i][l]!=undef&&uudata[k][j-1][i][l]!=undef)
					vordata[j-1][i-1]=
						(vvdata[k][j][i+1][l]-vvdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)-
						(uudata[k][j+1][i][l]-uudata[k][j-1][i][l])/(dy*2)+
						(uudata[k][j][i][l]/REarth)*ltan[ystart-1+j];
					
					else vordata[j-1][i-1]=undef;
				}
				
				do{
					for(int j=1;j<y-1;j++)
					for(int i=1;i<x-1;i++){
						if(vordata[j-1][i-1]!=undef
						 &&udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef){
							// calculate divergence
							float tmpdiv=
							(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
							(vdata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the divergence
							// multiplied by lcos to ensure the convergence of iteration
							float tmpu=-dxs[ystart-1+j]*tmpdiv/2*lcos[j];
							float tmpv=-dy*tmpdiv/2*lcos[j];
							
							// correct the u-velocity and v-velocity
							udata[k][j][i+1][l]+=tmpu;
							udata[k][j][i-1][l]-=tmpu;
							vdata[k][j+1][i][l]+=tmpv;
							vdata[k][j-1][i][l]-=tmpv;
							
							// calculate vorticity
							float tmpvor=
							(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)-
							(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2)+
							(udata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the new vorticity
							// multiplied by lcos to ensure the convergence of iteration
							tmpu=dy*(tmpvor-vordata[j-1][i-1])/2;
							tmpv=-dxs[ystart-1+j]*(tmpvor-vordata[j-1][i-1])/2;
							
							// correct the u-velocity and v-velocity
							udata[k][j+1][i][l]+=tmpu;
							udata[k][j-1][i][l]-=tmpu;
							vdata[k][j][i+1][l]+=tmpv;
							vdata[k][j][i-1][l]-=tmpv;
							
							// calculate the new divergence
							tmpdiv=
							(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
							(vdata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							float tmp=abs(tmpdiv);
							if(tmp>tmpmax) tmpmax=tmp;
						}
					}
					
					if(max<1e-12f) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<9999);
				
				System.out.println(String.format("Max error of divergence is %.4e Adjust for %4d loops",max,count));
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[k][j][i][l]!=undef){ aveu1+=udata[k][j][i][l]; cu++;}
					if(vdata[k][j][i][l]!=undef){ avev1+=vdata[k][j][i][l]; cv++;}
				}
				
				if(cu!=0) aveu1/=cu;	if(cv!=0) avev1/=cv;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[k][j][i][l]!=undef) udata[k][j][i][l]+=aveu0-aveu1;
					if(vdata[k][j][i][l]!=undef) vdata[k][j][i][l]+=avev0-avev1;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(vordata[j-1][i-1]==undef){
					udata[k][j][i][l]=undef;
					vdata[k][j][i][l]=undef;
				}
			}
		}
		
		System.out.println("finish calculating rotational velocity");
		
		return rw;
	}
	
	/**
     * Calculate rotational velocity using the vorticity, initial velocity field is given by u, v.
     *
     * @param	vor		vorticity of the observed velocity field (s^-1)
     * @param	 u		initial zonal velocity (m s^-1)
     * @param	 v		initial meridional velocity (m s^-1)
     *
     * @return	rw		rotational velocity, [0] is x component and [1] is y component
     */
	public Variable[] cRotationalVelocityByEndlich(Variable vor,Variable u,Variable v){
		System.out.println("\nstart calculating rotational velocity...");
		
 		checkDimensions(vor,u,v);
 		assignSubDomainParams(vor);
		
		Variable[] rw=new Variable[2];
		rw[0]=new Variable("us",vor);	rw[0].plus(u);
		rw[1]=new Variable("vs",vor);	rw[1].plus(v);
		
		rw[0].setCommentAndUnit("rotational velocity in x-direction (m s^-1)");
		rw[1].setCommentAndUnit("rotational velocity in y-direction (m s^-1)");
		
		float[][][][] vordata=vor.getData();
		float[][][][]  udata=rw[0].getData();
		float[][][][]  vdata=rw[1].getData();
		
		if(vor.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0,cu=0,cv=0;
				float max=1,tmpmax=0,aveu0=0,avev0=0,aveu1=0,avev1=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef){ aveu0+=udata[l][k][j][i]; cu++;}
					if(vdata[l][k][j][i]!=undef){ avev0+=vdata[l][k][j][i]; cv++;}
				}
				
				if(cu!=0) aveu0/=cu; cu=0;		if(cv!=0) avev0/=cv; cv=0;
				
				do{
					for(int j=1;j<y-1;j++)
					for(int i=1;i<x-1;i++){
						if(vordata[l][k][j][i]!=undef
						 &&udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef){
							// calculate divergence
							float tmpdiv=
							(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
							(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
							(vdata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the divergence
							float tmpu=-dxs[ystart-1+j]*tmpdiv/2;
							float tmpv=-dy*tmpdiv/2;
							
							// correct the u-velocity and v-velocity
							udata[l][k][j][i+1]+=tmpu;
							udata[l][k][j][i-1]-=tmpu;
							vdata[l][k][j+1][i]+=tmpv;
							vdata[l][k][j-1][i]-=tmpv;
							
							// calculate vorticity
							float tmpvor=
							(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)-
							(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)+
							(udata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the new vorticity
							tmpu=dy*(tmpvor-vordata[l][k][j][i])/2;
							tmpv=-dxs[ystart-1+j]*(tmpvor-vordata[l][k][j][i])/2;
							
							// correct the u-velocity and v-velocity
							udata[l][k][j+1][i]+=tmpu;
							udata[l][k][j-1][i]-=tmpu;
							vdata[l][k][j][i+1]+=tmpv;
							vdata[l][k][j][i-1]-=tmpv;
							
							// calculate the new divergence
							tmpdiv=
							(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
							(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
							(vdata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							float tmp=abs(tmpdiv);
							if(tmp>tmpmax) tmpmax=tmp;
						}
					}
					
					if(max<1e-12f) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<1999);
				
				System.out.println(String.format("Max error of vorticity is %.4e Adjust for %4d loops",max,count));
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef){ aveu1+=udata[l][k][j][i]; cu++;}
					if(vdata[l][k][j][i]!=undef){ avev1+=vdata[l][k][j][i]; cv++;}
				}
				
				if(cu!=0) aveu1/=cu;	if(cv!=0) avev1/=cv;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef) udata[l][k][j][i]+=aveu0-aveu1;
					if(vdata[l][k][j][i]!=undef) vdata[l][k][j][i]+=avev0-avev1;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(vordata[l][k][j][i]==undef){
					udata[l][k][j][i]=undef;
					vdata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0,cu=0,cv=0;
				float max=1,tmpmax=0,aveu0=0,avev0=0,aveu1=0,avev1=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[k][j][i][l]!=undef){ aveu0+=udata[k][j][i][l]; cu++;}
					if(vdata[k][j][i][l]!=undef){ avev0+=vdata[k][j][i][l]; cv++;}
				}
				
				if(cu!=0) aveu0/=cu; cu=0;		if(cv!=0) avev0/=cv; cv=0;
				
				do{
					for(int j=1;j<y-1;j++)
					for(int i=1;i<x-1;i++){
						if(vordata[l][k][j][i]!=undef
						 &&udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef){
							// calculate divergence
							float tmpdiv=
							(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
							(vdata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the divergence
							float tmpu=-dxs[ystart-1+j]*tmpdiv/2;
							float tmpv=-dy*tmpdiv/2;
							
							// correct the u-velocity and v-velocity
							udata[k][j][i+1][l]+=tmpu;
							udata[k][j][i-1][l]-=tmpu;
							vdata[k][j+1][i][l]+=tmpv;
							vdata[k][j-1][i][l]-=tmpv;
							
							// calculate vorticity
							float tmpvor=
							(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)-
							(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2)+
							(udata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the new vorticity
							tmpu=dy*(tmpvor-vordata[l][k][j][i])/2;
							tmpv=-dxs[ystart-1+j]*(tmpvor-vordata[l][k][j][i])/2;
							
							// correct the u-velocity and v-velocity
							udata[k][j+1][i][l]+=tmpu;
							udata[k][j-1][i][l]-=tmpu;
							vdata[k][j][i+1][l]+=tmpv;
							vdata[k][j][i-1][l]-=tmpv;
							
							// calculate the new divergence
							tmpdiv=
							(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
							(vdata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							float tmp=abs(tmpdiv);
							if(tmp>tmpmax) tmpmax=tmp;
						}
					}
					
					if(max<1e-12f) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<1999);
				
				System.out.println(String.format("Max error of vorticity is %.4e Adjust for %4d loops",max,count));
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[k][j][i][l]!=undef){ aveu1+=udata[k][j][i][l]; cu++;}
					if(vdata[k][j][i][l]!=undef){ avev1+=vdata[k][j][i][l]; cv++;}
				}
				
				if(cu!=0) aveu1/=cu;	if(cv!=0) avev1/=cv;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[k][j][i][l]!=undef) udata[k][j][i][l]+=aveu0-aveu1;
					if(vdata[k][j][i][l]!=undef) vdata[k][j][i][l]+=avev0-avev1;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(vordata[k][j][i][l]==undef){
					udata[k][j][i][l]=undef;
					vdata[k][j][i][l]=undef;
				}
			}
		}
		
		System.out.println("finish calculating rotational velocity");
		
		return rw;
	}
	
	/**
     * Endlich's method for calculating divergent velocity.
     * 
     * @param	u	original u-velocity (m s^-1)
     * @param	v	original v-velocity (m s^-1)
     *
     * @return	dw	divergent velocity, [0] is in x direction while [1] is in y direction
     */
	public Variable[] cDivergentVelocityByEndlich(Variable uu,Variable vv){
		System.out.println("\nstart calculating divergent velocity...");
		
 		checkDimensions(uu,vv);
 		assignSubDomainParams(uu);
 		
		Variable[] dw=new Variable[2];
		dw[0]=new Variable("ud",uu);	dw[0].setCommentAndUnit("divergent velocity in x-direction (m s^-1)");
		dw[1]=new Variable("vd",uu);	dw[1].setCommentAndUnit("divergent velocity in x-direction (m s^-1)");
		
		float[][]    divdata=new float[y-1][x-1];
		float[][][][]  udata=dw[0].getData();
		float[][][][]  vdata=dw[1].getData();
		float[][][][] uudata=uu.getData();
		float[][][][] vvdata=vv.getData();
		
		if(uu.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(uudata[l][k][j][i+1]!=undef&&uudata[l][k][j][i-1]!=undef&&vvdata[l][k][j+1][i]!=undef&&vvdata[l][k][j-1][i]!=undef)
						divdata[j-1][i-1]=
							(uudata[l][k][j][i+1]-uudata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
							(vvdata[l][k][j+1][i]-vvdata[l][k][j-1][i])/(dy*2)-
							(vvdata[l][k][j][i]/REarth)*ltan[ystart-1+j];
						
					else divdata[j-1][i-1]=undef;
				}
				
				int count=0,aveu=0,avev=0,cu=0,cv=0;
				float max=1,tmpmax=0;
				
				do{
					for(int j=1;j<y-1;j++)
					for(int i=1;i<x-1;i++){
						if(divdata[j-1][i-1]!=undef){
							// calculate vorticity
							float tmpvor=
							(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)-
							(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)+
							(udata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the new vorticity
							// multiplied by lcos to ensure the convergence of iteration
							float tmpu=dy*tmpvor/2*lcos[j];
							float tmpv=-dxs[ystart-1+j]*tmpvor/2*lcos[j];
							
							// correct the u-velocity and v-velocity
							udata[l][k][j+1][i]+=tmpu;
							udata[l][k][j-1][i]-=tmpu;
							vdata[l][k][j][i+1]+=tmpv;
							vdata[l][k][j][i-1]-=tmpv;
							
							// calculate divergence
							float tmpdiv=
							(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
							(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
							(vdata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the divergence
							// multiplied by lcos to ensure the convergence of iteration
							tmpu=-dxs[ystart-1+j]*(tmpdiv-divdata[j-1][i-1])/2*lcos[j];
							tmpv=-dy*(tmpdiv-divdata[j-1][i-1])/2*lcos[j];
							
							// correct the u-velocity and v-velocity
							udata[l][k][j][i+1]+=tmpu;
							udata[l][k][j][i-1]-=tmpu;
							vdata[l][k][j+1][i]+=tmpv;
							vdata[l][k][j-1][i]-=tmpv;
							
							// calculate the new divergence
							tmpvor=
							(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)-
							(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)+
							(udata[l][k][j][i]/REarth)*ltan[ystart-1+j];
							
							float tmp=abs(tmpvor);
							if(tmp>tmpmax) tmpmax=tmp;
							
						}
					}
					
					if(max<1e-13f) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<1999);
				
				System.out.println(String.format("Max error of vorticity is %.4e Adjust for %4d loops",max,count));
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef){ aveu+=udata[l][k][j][i]; cu++;}
					if(udata[l][k][j][i]!=undef){ avev+=vdata[l][k][j][i]; cv++;}
				}
				
				if(cu!=0) aveu/=cu;		if(cv!=0) avev/=cv;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef) udata[l][k][j][i]-=aveu;
					if(udata[l][k][j][i]!=undef) vdata[l][k][j][i]-=avev;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(divdata[j-1][i-1]==undef){
					udata[l][k][j][i]=undef;
					vdata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(uudata[k][j][i+1][l]!=undef&&uudata[k][j][i-1][l]!=undef&&vvdata[k][j+1][i][l]!=undef&&vvdata[k][j-1][i][l]!=undef)
						divdata[j-1][i-1]=
							(uudata[k][j][i+1][l]-uudata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
							(vvdata[k][j+1][i][l]-vvdata[k][j-1][i][l])/(dy*2)-
							(vvdata[k][j][i][l]/REarth)*ltan[ystart-1+j];
						
					else divdata[j-1][i-1]=undef;
				}
				
				int count=0,aveu=0,avev=0,cu=0,cv=0;
				float max=1,tmpmax=0;
				
				do{
					for(int j=1;j<y-1;j++)
					for(int i=1;i<x-1;i++){
						if(divdata[j-1][i-1]!=undef){
							// calculate vorticity
							float tmpvor=
							(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)-
							(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2)+
							(udata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the new vorticity
							// multiplied by lcos to ensure the convergence of iteration
							float tmpu=dy*tmpvor/2*lcos[j];
							float tmpv=-dxs[ystart-1+j]*tmpvor/2*lcos[j];
							
							// correct the u-velocity and v-velocity
							udata[k][j+1][i][l]+=tmpu;
							udata[k][j-1][i][l]-=tmpu;
							vdata[k][j][i+1][l]+=tmpv;
							vdata[k][j][i-1][l]-=tmpv;
							
							// calculate divergence
							float tmpdiv=
							(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
							(vdata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							// calculate the correct value according to the divergence
							// multiplied by lcos to ensure the convergence of iteration
							tmpu=-dxs[ystart-1+j]*(tmpdiv-divdata[j-1][i-1])/2*lcos[j];
							tmpv=-dy*(tmpdiv-divdata[j-1][i-1])/2*lcos[j];
							
							// correct the u-velocity and v-velocity
							udata[k][j][i+1][l]+=tmpu;
							udata[k][j][i-1][l]-=tmpu;
							vdata[k][j+1][i][l]+=tmpv;
							vdata[k][j-1][i][l]-=tmpv;
							
							
							// calculate the new divergence
							tmpdiv=
							(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
							(vdata[k][j][i][l]/REarth)*ltan[ystart-1+j];
							
							float tmp=abs(tmpdiv);
							if(tmp>tmpmax) tmpmax=tmp;
							
						}
					}
					
					if(max<1e-13f) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<1999);
				
				System.out.println(String.format("Max error of vorticity is %.4e Adjust for %4d loops",max,count));
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[k][j][i][l]!=undef){ aveu+=udata[k][j][i][l]; cu++;}
					if(udata[k][j][i][l]!=undef){ avev+=vdata[k][j][i][l]; cv++;}
				}
				
				if(cu!=0) aveu/=cu;		if(cv!=0) avev/=cv;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(udata[l][k][j][i]!=undef) udata[l][k][j][i]-=aveu;
					if(udata[l][k][j][i]!=undef) vdata[l][k][j][i]-=avev;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(divdata[j-1][i-1]==undef){
					udata[k][j][i][l]=undef;
					vdata[k][j][i][l]=undef;
				}
			}
		}
		
		System.out.println("finish calculating divergent velocity");
		
		return dw;
	}
	
	
	/**
     * Calculate the streamfunction using Endlich method.
     *
     * @param	u	non-divergence zonal velocity (m s^-1)
     * @param	v	non-divergence meridinoal velocity (m s^-1)
     *
     * @return	sf	streamfunction (m^2 s^-1)
     */
	public Variable cStreamFunctionByEndlich(Variable uu,Variable vv){
		System.out.println("\nstart calculating stream function...");
		
 		checkDimensions(uu,vv);
 		assignSubDomainParams(uu);
 		
		Variable sf=new Variable("sf",uu);
		sf.setCommentAndUnit("stream function (m^2 s^-1)");
		
		float[][][][]  udata=uu.getData();
		float[][][][]  vdata=vv.getData();
		float[][][][] sfdata=sf.getData();
		
		if(uu.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;
				float tmpmax=0,max=Float.MAX_VALUE,error=0.1f;
				
				do{
					for(int j=0;j<y;j++)
					for(int i=1;i<x-1;i++)
					if(vdata[l][k][j][i]!=uu.getUndef()){
						float dsf=dxs[ystart-1+j]*vdata[l][k][j][i]-(sfdata[l][k][j][i+1]-sfdata[l][k][j][i-1])/2;
						
						sfdata[l][k][j][i+1]+=dsf;
						sfdata[l][k][j][i-1]-=dsf;
						
						float tmp=abs(dsf);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					for(int j=1;j<y-1;j++)
					for(int i=0;i<x;i++)
					if(udata[l][k][j][i]!=uu.getUndef()){
						float dsf=-dy*udata[l][k][j][i]-(sfdata[l][k][j+1][i]-sfdata[l][k][j-1][i])/2;
						
						sfdata[l][k][j+1][i]+=dsf;
						sfdata[l][k][j-1][i]-=dsf;
						
						float tmp=abs(dsf);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					if(max<error||max==tmpmax) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<7999);
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
				if(udata[l][k][j][i]==undef) sfdata[l][k][j][i]=undef;
				
				System.out.println(String.format("Max adjustment is %13.4f, adjust for %4d loops",max,count));
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;
				float tmpmax=0,max=Float.MAX_VALUE,error=0.1f;
				
				do{
					for(int j=0;j<y;j++)
					for(int i=1;i<x-1;i++){
						float dsf=dxs[ystart-1+j]*vdata[k][j][i][l]-(sfdata[k][j][i+1][l]-sfdata[k][j][i-1][l])/2;
						
						sfdata[l][k][j][i+1]+=dsf;
						sfdata[l][k][j][i-1]-=dsf;
						
						float tmp=abs(dsf);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					for(int j=1;j<y-1;j++)
					for(int i=0;i<x;i++){
						float dsf=-dy*udata[k][j][i][l]-(sfdata[k][j+1][i][l]-sfdata[k][j-1][i][l])/2;
						
						sfdata[k][j+1][i][l]+=dsf;
						sfdata[k][j-1][i][l]-=dsf;
						
						float tmp=abs(dsf);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					if(max<error||max==tmpmax) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<7999);
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
				if(udata[k][j][i][l]==undef) sfdata[k][j][i][l]=undef;
				
				System.out.println(String.format("Max adjustment is %13.4f, adjust for %4d loops",max,count));
			}
		}
		
		System.out.println("finish calculating stream function");
		
		return sf;
	}
	
	/**
     * Calculate the velocity potential using Endlich method.
     *
     * @param	u	irrotational u-velocity (m s^-1)
     * @param	v	irrotational v-velocity (m s^-1)
     *
     * @return	vp	velocity potential (m^2 s^-1)
     */
	public Variable cVelocityPotentialByEndlich(Variable uu,Variable vv){
		System.out.println("\nstart calculating potential function....");
		
 		checkDimensions(uu,vv);
 		assignSubDomainParams(uu);
 		
		Variable[] dw=cDivergentVelocityByEndlich(uu,vv);
		Variable   vp=new Variable("vp",uu);
		vp.setCommentAndUnit("velocity potential (m^2 s^-1)");
		
		float[][][][]  udata=dw[0].getData();
		float[][][][]  vdata=dw[1].getData();
		float[][][][] vpdata=   vp.getData();
		
		if(uu.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;
				float tmpmax=0,max=Float.MAX_VALUE,error=0.1f;
				
				do{
					for(int j=0;j<y;j++)
					for(int i=1;i<x-1;i++)
					if(udata[l][k][j][i]!=uu.getUndef()){
						float dvp=dxs[ystart-1+j]*udata[l][k][j][i]-(vpdata[l][k][j][i+1]-vpdata[l][k][j][i-1])/2;
						
						vpdata[l][k][j][i+1]+=dvp;
						vpdata[l][k][j][i-1]-=dvp;
						
						float tmp=abs(dvp);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					for(int j=1;j<y-1;j++)
					for(int i=0;i<x;i++)
					if(vdata[l][k][j][i]!=uu.getUndef()){
						float dvp=dy*vdata[l][k][j][i]-(vpdata[l][k][j+1][i]-vpdata[l][k][j-1][i])/2;
						
						vpdata[l][k][j+1][i]+=dvp;
						vpdata[l][k][j-1][i]-=dvp;
						
						float tmp=abs(dvp);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					if(max<error||max==tmpmax) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<7999);
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
				if(udata[l][k][j][i]==undef) vpdata[l][k][j][i]=undef;
				
				System.out.println(String.format("Max adjustment is %11.4f, adjust for %4d loops",max,count));
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;
				float tmpmax=0,max=Float.MAX_VALUE,error=0.1f;
				
				do{
					for(int j=0;j<y;j++)
					for(int i=1;i<x-1;i++)
					if(udata[k][j][i][l]!=uu.getUndef()){
						float dvp=dxs[ystart-1+j]*udata[k][j][i][l]-(vpdata[k][j][i+1][l]-vpdata[k][j][i-1][l])/2;
						
						vpdata[l][k][j][i+1]+=dvp;
						vpdata[l][k][j][i-1]-=dvp;
						
						float tmp=abs(dvp);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					for(int j=1;j<y-1;j++)
					for(int i=0;i<x;i++)
					if(vdata[k][j][i][l]!=uu.getUndef()){
						float dvp=dy*vdata[k][j][i][l]-(vpdata[k][j+1][i][l]-vpdata[k][j-1][i][l])/2;
						
						vpdata[k][j+1][i][l]+=dvp;
						vpdata[k][j-1][i][l]-=dvp;
						
						float tmp=abs(dvp);	if(tmp>tmpmax) tmpmax=tmp;
					}
					
					if(max<error||max==tmpmax) break;	max=tmpmax;		tmpmax=0;
					
				}while(count++<7999);
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
				if(udata[k][j][i][l]==undef) vpdata[k][j][i][l]=undef;
				
				System.out.println(String.format("Max adjustment is %11.4f, adjust for %4d loops",max,count));
			}
		}
		
		System.out.println("finish calculating potential function");
		
		return vp;
	}
	
	/**
     * Calculate the stream function using SOR.
     *
     * @param	u	original u-velocity (m s^-1)
     * @param	v	original v-velocity (m s^-1)
     *
     * @return	sf	streamfunction
     */
	public Variable cStreamFunctionBySOR(Variable u,Variable v){
		SphericalSpatialModel ssm=(SphericalSpatialModel)(sm);
		DynamicMethodsInSC     dm=new DynamicMethodsInSC(ssm);
 		PoissonEquationInSC    pe=new PoissonEquationInSC(ssm);
 		
 		pe.setDimComBination(DimCombination.XY);
 		
 		Variable vor=dm.c2DVorticity(u,v);
 		Variable re =pe.invertingBySOR(vor);
 		
 		re.setName("sf");
 		re.setCommentAndUnit("streamfunction in X-Y plane (m^2 s^-1)");
		
		return re;
	}
	
	/**
     * Calculate the velocity potential using SOR.
     *
     * @param	u	original u-velocity (m s^-1)
     * @param	v	original v-velocity (m s^-1)
     *
     * @return	vp	velocity potential (m^2 s^-1)
     */
	public Variable cVelocityPotentialBySOR(Variable u,Variable v){
		SphericalSpatialModel ssm=(SphericalSpatialModel)(sm);
		DynamicMethodsInSC     dm=new DynamicMethodsInSC(ssm);
 		PoissonEquationInSC    pe=new PoissonEquationInSC(ssm);
 		
 		pe.setDimComBination(DimCombination.XY);
		
		Variable div=dm.c2DDivergence(u,v);
		Variable re =pe.invertingBySOR(div);
 		
 		re.setName("vp");
 		re.setCommentAndUnit("velocity potential in X-Y plane (m^2 s^-1)");
		
		return re;
	}
	
	
	/**
     * Calculate the stream function using SOR.
     *
     * @param	v	v-velocity (m s^-1)
     * @param	w	w-velocity (m s^-1)
     *
     * @return	sf	streamfunction (m^2 s^-1)
     */
	public Variable cYPStreamFunctionBySOR(Variable v,Variable w){
		SphericalSpatialModel ssm=(SphericalSpatialModel)(sm);
 		DynamicMethodsInSC     dm=new DynamicMethodsInSC (ssm);
 		PoissonEquationInSC    pe=new PoissonEquationInSC(ssm);
 		
 		pe.setDimComBination(DimCombination.YZ);
 		
 		Variable vor=dm.cYZVorticity(v,w);
 		Variable re =pe.invertingBySOR(vor);
 		
 		re.setName("sf");
 		re.setCommentAndUnit("streamfunction in Y-Z plane (m Pa s^-1)");
		
		return re;
	}
	
	
	/** test
	public static void main(String arg[]){
		DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/Aviso/msla.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		VelocityFieldInSC wf=new VelocityFieldInSC(ssm);
		
		Range r=new Range("",dd);
		
		Variable sla=df.getVariables(r,"sla")[0];
		
		Variable[] uv=wf.cGeostrophicVelocity(sla);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Data/Aviso/geo.dat");
		dw.writeData(dd,uv[0],uv[1],sla);	dw.closeFile();
	}*/
}
