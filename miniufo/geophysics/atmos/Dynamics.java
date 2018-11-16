/**
 * @(#)Dynamics.java	1.0 2018.09.27
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.geophysics.atmos;

import static miniufo.diagnosis.SpatialModel.gEarth;


/**
 * This class contains algorithms for the most basic diagnostic quantities
 * in geophysical fluids such as the atmosphere and the ocean.
 * 
 * These algorithms ignore particular partial difference schemes by assuming
 * the partial differential derivatives are all known.
 *
 * @version 1.0, 2018.09.27
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Dynamics{
	//
	
	/**
	 * prevent from instantiate
	 */
	private Dynamics(){}
	
	
	/**
	 * Calculate Montgomery streamfunction.
	 * 
	 * @param	z	geopotential height (m)
	 * @param	T	temperature (K)
	 */
	public static float cMontgomeryStreamfunction(float z,float T){ return (float)(ThermoDynamics.Cp*T+gEarth*z);}
	
	/**
	 * Calculate geopotential height.
	 * 
	 * @param	m	Montgomery streamfunction (m^2 s^-2)
	 * @param	T	temperature (K)
	 */
	public static float cGeopotentialHeight(float m,float T){ return (float)((m-ThermoDynamics.Cp*T)/gEarth);}
}
