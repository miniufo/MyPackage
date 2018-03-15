/**
 * @(#)AttachedMeta.java	1.0 2018.03.14
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;


/**
 * Used to describe the Meta-info of the attached Lagrangian data
 *
 * @version 1.0, 2018.03.14
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AttachedMeta{
	//
	public final  int  index;	// attached index
	public final String name;	// descriptive name
	
	
	/**
	 * constructor
	 */
	public AttachedMeta(String name,int index){
		this.name = name;
		this.index=index;
	}
	
	
	/**
     * used to print out
     */
	public String toString(){ return name+" ("+index+")";}
}
