/**
 * @(#)EllipticEquationInSC.java	1.0 2015.03.13
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application;

import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInSphericalCoordinate;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.test.application.basic.LaplaceEquationInSC;


/**
 * basic class for elliptic equation in spherical coordinate
 *
 * @version 1.0, 2015.03.13
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class EllipticEquationInSC extends EquationInSphericalCoordinate
implements EllipticEquationInterface{
	//
	protected Variable A=null;
	protected Variable B=null;
	protected Variable C=null;
	
	protected LaplaceEquationInSC le=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a given spacial model in spheral coordinate
     */
	public EllipticEquationInSC(SphericalSpatialModel ssm){ super(ssm);}
}
