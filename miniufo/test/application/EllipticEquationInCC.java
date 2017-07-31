/**
 * @(#)EllipticEquationInCC.java	1.0 2015.03.13
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application;

import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInCylindricalCoordinate;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.Variable;


/**
 * basic class for elliptical equation in polar coordinate
 *
 * @version 1.0, 2015.03.13
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class EllipticEquationInCC
extends EquationInCylindricalCoordinate implements EllipticEquationInterface{
	//
	protected Variable A=null;
	protected Variable B=null;
	protected Variable C=null;
	
	/**
     * constructor
     *
     * @param	psm		initialized by a given spacial model in polar coordinate
     */
	public EllipticEquationInCC(CylindricalSpatialModel csm){ super(csm);}
}
