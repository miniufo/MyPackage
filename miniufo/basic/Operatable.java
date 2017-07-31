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
public interface Operatable<T> extends BasicOperatable<T>{
	//
	public T pow(float f);
	public T pow(T o);
	
	public T powEq(float f);
	public T powEq(T o);
	
	public T abs();
	public T absEq();
	
	public T square();
	public T squareEq();
	
	public T sqrt();
	public T sqrtEq();
	
	public T exp();
	public T expEq();
	
	public T reciprocal();
	public T reciprocalEq();
}
