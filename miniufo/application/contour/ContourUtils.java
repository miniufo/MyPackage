/**
 * @(#)ContourUtils.java	1.0 2017.06.09
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;


/**
 * This class provides utility functions for the package.
 *
 * @version 1.0, 2017.06.09
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ContourUtils{
	
	/**
	 * compute the equivalent circle (minimum Le) given
	 * an area A = PI * R^2 (usually in spherical coordinates).
	 * 
	 * C = 2*PI*R = 2*PI*sqrt(A/PI) = 2*sqrt(A*PI)
	 * 
	 * @param	area	a given area, probably enclosed by a tracer contour (m^2)
	 */
	public static double cEquivalentCircleLength(float area){ return 2.0*Math.sqrt(area*Math.PI);}
	
	/**
	 * compute the equivalent circle radius R given
	 * an area A = PI * R^2 (usually in spherical coordinates).
	 * 
	 * @param	area	a given area, probably enclosed by a tracer contour (m^2)
	 */
	public static double cEquivalentCircleRadius(float area){ return Math.sqrt(area/Math.PI);}
	
	/**
	 * compute the equivalent Y (minimum Le) given
	 * an area A = X*Y (usually in Cartesian coordinates).
	 * 
	 * @param	area	a given area, probably enclosed by a tracer contour (m^2)
	 */
	public static double cEquivalentY(float area,float X){ return area/X;}
	
	
	/*** test
	public static void main(String[] args){
		
	}*/
}
