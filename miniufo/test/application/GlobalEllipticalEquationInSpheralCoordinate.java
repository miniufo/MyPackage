/**
 * @(#)EllipticalEquationInSpheralCoordinate.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application;

import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInSphericalCoordinate;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.test.application.basic.GlobalLaplaceEquationInSC;


/**
 * basic class for elliptical equation in spheral coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class GlobalEllipticalEquationInSpheralCoordinate extends EquationInSphericalCoordinate
implements EllipticEquationInterface{
	//
	protected Variable A=null;
	protected Variable B=null;
	protected Variable C=null;
	
	protected GlobalLaplaceEquationInSC le=null;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a given spacial model in spheral coordinate
     */
	public GlobalEllipticalEquationInSpheralCoordinate(SphericalSpatialModel ssm){ super(ssm);}
}
