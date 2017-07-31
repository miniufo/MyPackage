/**
 * @(#)Mapping.java	1.0 2014.03.29
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.basic;


/**
 * general mapping interface
 *
 * @version 1.0, 2014.03.29
 * @author  MiniUFO
 * @since   MDK1.0
 */
@FunctionalInterface
public interface Mapping{
	
	/**
     * Applies this function to given arguments
     *
     * @param	arg		an argument
     */
	float map(float... args);
}
