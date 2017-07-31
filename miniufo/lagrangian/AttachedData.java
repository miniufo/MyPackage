/**
 * @(#)AttachedData.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;


/**
 * used to describe the data structure attached to one Record
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AttachedData{
	//
	private float value=0;		// value of the data
	
	private String name=null;	// descriptive name
	
	
	/**
	 * constructor
	 */
	public AttachedData(float value,String name){
		this.value=value;
		this.name =name;
	}
	
	
	/*** getor and setor ***/
	public float getValue(){ return value;}
	
	public String getName(){ return name;}
	
	public void setName(String name){ this.name=name;}
	
	public void setValue(float value){ this.value=value;}
	
	
	/**
     * used to print out
     */
	public String toString(){
		return name+": "+value;
	}
}
