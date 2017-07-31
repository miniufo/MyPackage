/**
 * @(#)Function2D	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;


/**
 * define a 2-D function
 */
public interface Function2D{
	
	/**
	 * define a 2-Dimensional function like:
	 * Z = F (X , Y)
	 */
	public float map(float x,float y);
}
