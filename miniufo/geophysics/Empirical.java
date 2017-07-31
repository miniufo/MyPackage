/**
 * @(#)Empirical.java	1.0 2016.10.12
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.geophysics;

import static java.lang.Math.exp;


/**
 * This class contains some empirical formula used in geophysical fluid studies.
 *
 * @version 1.0, 2016.10.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Empirical{
	//
	public enum Scheme{
		LP,	// Reference: Large and Pond (1981, JPO)
		DO,	// Reference: Donelan et al. (2004, GRL)
		LY	// Reference: Large and Yeager (2009, CD)
	};
	
	
	/**
	 * prevent from instantiate
	 */
	private Empirical(){}
	
	
	/**
     * Calculate bulk drag coefficient Cd under different wind speed.
     *
     * @param	wind	wind speed (m s^-1)
     * @param	scheme	which scheme to use
     *
     * @return	drag coefficient Cd
     */
	public static float cBulkDragCoefficient(float wind,Scheme scheme){
		switch(scheme){
		case LP: return cCdByLP(wind);
		case DO: return cCdByDO(wind);
		case LY: return cCdByLY(wind);
		default: throw new IllegalArgumentException("unsupported scheme:"+scheme);
		}
	}
	
	
	/**
	 * Reference: Large and Pond (1981, JPO)
	 * 
	 * @param wind	non-negative wind speed (m s^-1)
	 */
	public static float cCdByLP(float wind){
		if(wind<0) throw new IllegalArgumentException("wind speed should be non-negative");
		
		if(wind<11) return 1.2e-3f;
		else return 4.9e-4f+6.5e-5f*wind;
	}
	
	/**
	 * Reference: Donelan et al. (2004, GRL)
	 * 
	 * @param wind	non-negative wind speed (m s^-1)
	 */
	public static float cCdByDO(float wind){
		if(wind<0) throw new IllegalArgumentException("wind speed should be non-negative");
		
		if(wind<30) return 3.55e-4f+6.5e-5f*wind;
		else return 3.55e-4f+6.5e-5f*30f;
	}
	
	/**
	 * Reference: Large and Yeager (2009, CD)
	 * 
	 * @param wind	non-negative wind speed (m s^-1)
	 */
	public static float cCdByLY(float wind){
		if(wind<0) throw new IllegalArgumentException("wind speed should be non-negative");
		if(wind==0) return 0;
		
		double a1=2.7e-3;	double a2=1.42e-4;
		double a3=7.64e-5;	double a4=-3.14807e-13;
		
		if(wind<33) return (float)(a1/wind+a2+a3*wind+a4*Math.pow(wind,6.0));
		else return 2.34e-3f;
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
	public static float cMPIAtl1(float sst){ return (float)(0.51444*74.0*exp(0.2*(sst-273.15-25.0)));}
	
	/**
     * Calculate the maximum potential intensity of a storm according to the underlying SST.
     * This method is based on DeMaria et al. (1994, JC), applicable to the Atlantic hurricanes.
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity in unit of surface wind speed (m s^-1)
     */
	public static float cMPIAtl2(float sst){ return (float)(28.2+55.8*exp(0.1813*(sst-273.15-30)));}
	
	/**
     * Calculate the maximum potential intensity of a storm according to the underlying SST.
     * This method is based on Zeng (2007, MWR), applicable to TCs over western North Pacific.
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity in unit of surface wind speed (m s^-1)
     */
	public static float cMPIWNP(float sst){ return (float)(15.69+(98.03*exp(0.1806*(sst-30.0))));}
	
	
	/*** test**
	public static void main(String[] args){
		int len=50;
		
		float[] wind =new float[len];
		float[] tauLP=new float[len];
		float[] tauLY=new float[len];
		float[] tauLI=new float[len];
		float[] tauDO=new float[len];
		
		for(int i=0;i<len;i++){
			wind[i]=i;
			
			tauLP[i]=cCdByLP(wind[i]);
			tauLY[i]=cCdByLY(wind[i]);
			tauLI[i]=cCdByMIT(wind[i]);
			tauDO[i]=cCdByDO(wind[i]);
			
			System.out.println(wind[i]+"\t"+tauLP[i]+"\t"+tauLY[i]+"\t"+tauLI[i]+"\t"+tauDO[i]);
		}
	}*/
}
