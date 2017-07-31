/**
 * @(#)DataRead.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import miniufo.diagnosis.Variable;


/**
 * interface for I/O
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public interface DataRead extends Print{
	
	/**
	 * to read data from the specified file
	 *
     * @param	v	variable need to fill with data
     */ 
	public void readData(Variable... v);
	
	/**
	 * close file method
     */
	public void closeFile();
}
