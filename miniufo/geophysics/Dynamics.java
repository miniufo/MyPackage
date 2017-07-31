/**
 * @(#)Dynamics.java	1.0 2015.01.07
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.geophysics;


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
public final class Dynamics{
	
	/**
	 * prevent from instantiate
	 */
	private Dynamics(){}
	
	
	/**
     * Calculate 2D (horizontal) divergence.
     *
     * @param	Ux	x-gradient of x-velocity (s^-1)
     * @param	Vy	y-gradient of y-velocity (s^-1)
     *
     * @return	2D (horizontal) divergence (s^-1)
     */
	public float cDivergence2D(float Ux,float Vy){ return Ux+Vy;}
	
	/**
     * Calculate 3D divergence.
     *
     * @param	Ux	x-gradient of x-velocity (s^-1)
     * @param	Vy	y-gradient of y-velocity (s^-1)
     * @param	Wz	z-gradient of z-velocity (s^-1)
     *
     * @return	3D divergence (s^-1)
     */
	public float cDivergence3D(float Ux,float Vy,float Wz){ return Ux+Vy+Wz;}
	
	/**
     * Calculate 2D (horizontal) vorticity.
     *
     * @param	Vx	x-gradient of y-velocity (s^-1)
     * @param	Uy	y-gradient of x-velocity (s^-1)
     *
     * @return	2D (horizontal) vorticity (s^-1)
     */
	public float cVorticity2D(float Vx,float Uy){ return Vx-Uy;}
	
	/**
     * Calculate 3D vorticity.
     *
     * @param	Wy	y-gradient of z-velocity (s^-1)
     * @param	Vz	z-gradient of y-velocity (s^-1)
     * @param	Uz	z-gradient of x-velocity (s^-1)
     * @param	Wx	x-gradient of z-velocity (s^-1)
     * @param	Vx	x-gradient of y-velocity (s^-1)
     * @param	Uy	y-gradient of x-velocity (s^-1)
     *
     * @return	vor	3D vorticity vector (s^-1),
     * 				[0], [1] and [2] for x-, y- and z-components respectively
     */
	public float[] cVorticity3D(float Wy,float Vz,float Uz,float Wx,float Vx,float Uy){
		return new float[]{Wy-Vz, Uz-Wx, Vx-Uy};
	}
	
	/**
     * Calculate 2D (horizontal) tension strain (stretching deformation).
     *
     * @param	Ux	x-gradient of x-velocity (s^-1)
     * @param	Vy	y-gradient of y-velocity (s^-1)
     *
     * @return	2D (horizontal) tension strain (s^-1)
     */
	public float cTensionStrain2D(float Ux,float Vy){ return Ux-Vy;}
	
	/**
     * Calculate 2D (horizontal) shearing strain (shearing deformation).
     *
     * @param	Vx	x-gradient of y-velocity (s^-1)
     * @param	Uy	y-gradient of x-velocity (s^-1)
     *
     * @return	2D (horizontal) tension strain (s^-1)
     */
	public float cShearingStrain2D(float Vx,float Uy){ return Vx+Uy;}
	
	/**
     * Calculate 2D (horizontal) total strain (deformation rate).
     *
     * @param	ts	tension strain (stretching deformation) (s^-1)
     * @param	ss	shearing strain (shearing deformation) (s^-1)
     *
     * @return	2D (horizontal) total strain (s^-1)
     */
	public float cTotalStrain2D(float ts,float ss){ return (float)Math.hypot(ts,ss);}
	
	/**
     * Calculate 2D (horizontal) total strain angle defined between
     * the x axis and the axis of dilatation.
     * 
     * Reference: Spensberger and Spengler 2014, JAS
     *
     * @param	ts	tension strain (stretching deformation) (s^-1)
     * @param	ss	shearing strain (shearing deformation) (s^-1)
     *
     * @return	2D (horizontal) total strain angle (radian)
     */
	public float cTotalStrainAngle(float ts,float ss){ return (float)(Math.atan2(ss,ts)/2.0);}
}
