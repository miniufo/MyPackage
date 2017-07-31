/**
 * @(#)EllipticEquationInterface.java	1.0 2015.03.13
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application;

import miniufo.diagnosis.Variable;


/**
 * interface that a elliptic equation should implement
 * Reference: My Notes on inverting elliptical equation using SOR
 *
 * @version 1.0, 2015.03.13
 * @author  MiniUFO
 * @since   MDK1.0
 */
public interface EllipticEquationInterface{
	
	/**
	 * compute A defined at model grids
	 * 
	 * @param	v	a given variable probably used in the computation
	 * 
	 * @return	A	coefficient of the elliptic equation defined at model grids
	 */
	public Variable cA(Variable v);
	
	/**
	 * compute B defined at model grids
	 * 
	 * @param	v	a given variable probably used in the computation
	 * 
	 * @return	B	coefficient of the elliptic equation defined at model grids
	 */
	public Variable cB(Variable v);
	
	/**
	 * compute C defined at model grids
	 * 
	 * @param	v	a given variable probably used in the computation
	 * 
	 * @return	C	coefficient of the elliptic equation defined at model grids
	 */
	public Variable cC(Variable v);
	
	/**
	 * compute A' defined at semi-grids in the first dimension
	 * 
	 * @param	v	a given variable probably used in the computation
	 * 
	 * @return	A'	coefficient of the elliptic equation suitable for SOR
	 */
	public Variable cAPrime(Variable v);
	
	/**
	 * compute B' defined at model grids
	 * 
	 * @param	v	a given variable probably used in the computation
	 * 
	 * @return	B'	coefficient of the elliptic equation suitable for SOR
	 */
	public Variable cBPrime(Variable v);
	
	/**
	 * compute C' defined at semi-grids in the second dimension
	 * 
	 * @param	v	a given variable probably used in the computation
	 * 
	 * @return	C'	coefficient of the elliptic equation suitable for SOR
	 */
	public Variable cCPrime(Variable v);
}
