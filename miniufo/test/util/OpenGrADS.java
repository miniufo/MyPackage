/**
 * @(#)OpenGrADS.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;


/**
 * use to implement some GrADS functions
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class OpenGrADS{
	
	/**
     * constructor prevent initialization
     */
	private OpenGrADS(){}
	
	
	/**
     * call GrADS to run gs in batch mode
     */
	public static void runGS(String gspath){
		Utility.runShell("opengrads -pbc '"+gspath+"'");
	}
	
	
	/** test
	public static void main(String[] args){
		runGS("e:/paint.gs");
	}*/
}
