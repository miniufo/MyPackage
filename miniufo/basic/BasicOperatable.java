/**
 * @(#)Operatable.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.basic;


/**
 * an interface for all kinds of operations
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public interface BasicOperatable<T>{
	//
	public T plus(float f);
	public T plus(T o);
	
	public T plusEq(float f);
	public T plusEq(T o);
	
	public T minus(float f);
	public T minus(T o);
	
	public T minusEq(float f);
	public T minusEq(T o);
	
	public T multiply(float f);
	public T multiply(T o);
	
	public T multiplyEq(float f);
	public T multiplyEq(T o);
	
	public T divide(float f);
	public T divide(T o);
	
	public T divideEq(float f);
	public T divideEq(T o);
}
