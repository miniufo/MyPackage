/**
 * @(#)Dynamics.java	1.0 2015.01.07
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.geophysics.atmos;

import static miniufo.diagnosis.SpatialModel.gEarth;


/**
 * This class contains algorithems for the most basic diagnostic quantities
 * in geophysical fluids such as the atmosphere and the ocean, adopting the most
 * commomly used Cartesian coordinates.
 * 
 * These algorithems ignore particular partial difference schemes by assuming
 * the partial differential derivatives are all known.
 *
 * @version 1.0, 2015.01.07
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ThermoDynamics{
	//
	public static final float RHO  =1.225f;		// density of air at sea level at a temperature of 15 degree
	public static final float Cp   =1004.88f;	// thermal capacity of ideal gas (J kg^-1 K^-1)
	public static final float E0   =610.78f;	// saturated vapor pressure above pure water surface (Pa)
	public static final float Rd   =287.04f;	// constant of dry air (J kg^-1 K^-1)
	public static final float Rv   =461.52f;	// constant of vapor   (J kg^-1 K^-1)
	public static final float rd   =gEarth/Cp;	// vertical temperature decrease rate
	public static final float Cv   =Cp-Rd;		// (J kg^-1 K^-1)
	public static final float kappa=Rd/Cp;		// kappa
	public static final float L0   =2500794f;	// (J kg^-1)
	public static final float Pref =100000;		// 1000 hPa of pressure reference
	
	private static final float[] lvap={
		2.3823f,2.3945f,2.4062f,2.4183f,2.4300f,
		2.4418f,2.4535f,2.4656f,2.4774f,2.4891f,
		2.5008f,2.5247f,2.5494f,2.5749f,2.6030f,2.6348f
	};
	
	
	/**
	 * prevent from instantiate
	 */
	private ThermoDynamics(){}
	
	
	/**
     * calculate L in -Ldq/dt, as Prof. Yuan's, taking 'T<273.15f' into account
     * cL2 is a linear method which do not consider the 'T<273.15f' case
     *
     * @param	T	temperature (K)
     * 
     * @return	L	latent heating coefficient (J kg^-1)
     */
	public static float cL(float T){
		int itx,iarg;	float fctr;
		
		T-=273.15f;
		
		if(T>0){
			itx=(int)(T/5f);
			iarg=11-itx;
			fctr=(5f*(itx+1)-T)/5f;
			
		}else{
			itx=(int)(T/(-10f))+1;
			iarg=11+itx;
			
			if(iarg<=16) fctr=(10f*(1-itx)-T)/10f;
			else{
				iarg=16;
				fctr=(-50f-T)/10f+1f;
			}
		}
		
		return (lvap[iarg-2]+fctr*(lvap[iarg-1]-lvap[iarg-2]))*1.e6f;
	}
	
	/**
     * a linear method that does not consider the 'T<273.15f' case
     *
     * @param	T	temperature (K)
     * 
     * @return	L	latent heating coefficient (J kg^-1)
     */
	public static float cL2(float T){ return L0-2327f*(T-273.15f);}
	
	
	/**
     * Calculate 2D (horizontal) divergence.
     *
     * @param	Ux	x-gradient of x-velocity (s^-1)
     * @param	Vy	y-gradient of y-velocity (s^-1)
     *
     * @return	2D (horizontal) divergence (s^-1)
     */
	public static float cSpecificHumidity(){ return 1;}
	
	public static float cSaturatedSpecificHumidity(){ return 1;}
	
	public static float cRelativeHumidity(){ return 1;}
	
	public static float cSaturatedVaporPressure(){ return 1;}
	
	public static float cSaturatedVaporPressureOnSurface(){ return 1;}
	
	public static float cSpecificVolume(){ return 1;}
	
	public static float cStaticStabilityArgByT(){ return 1;}
	
	public static float cDewPointTemperature(){ return 1;}
	
	public static float cLCLTemperature(){ return 1;}
	
	public static float cEquivalentTemperature(){ return 1;}
	
	public static float cVirtualTemperature(float T,float q){
		return T*(1f+0.607717f*q);
	}
	
	/**
	 * Calculate potential temperature reference to 1000 hPa
	 * 
	 * @param	T	temperature (K)
	 * @param	p	pressure (Pa)
	 */
	public static float cPotentialTemperature(float T,float p){
		return (float)(T*Math.pow(Pref/p,kappa));
	}
	
	/**
	 * Calculate temperature from potential temperature referenced to 1000 hPa
	 * 
	 * @param	theta	potential temperature (K)
	 * @param	p		pressure (Pa)
	 */
	public static float cTemperature(float theta,float p){
		return (float)(theta*Math.pow(p/Pref,kappa));
	}
	
	public static float cEquivalentPotentialTemperature(){ return 1;}
	
	public static float cExnerFunction(float p){
		return (float)(Cp*Math.pow(p/Pref,Rd/Cp));
	}
}
