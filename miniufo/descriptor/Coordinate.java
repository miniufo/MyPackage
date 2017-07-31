/**
 * @(#)Coordinate.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;


/**
 * discrete sampled coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class Coordinate implements Cloneable{
	//
	protected boolean islinear=true;
	protected boolean isIncre =true;
	
	protected  String   name  =null; 
	
	
	/**
     * constructor
     *
     * param	samples		an array of samples
     */
	public Coordinate(String name){ this.name=name;}
	
	
	/*** getor and setor ***/
	public boolean isLinear(){ return islinear;}
	
	public boolean isIncreasing(){ return isIncre;}
	
	public String getName(){ return name;}
	
	public abstract int length();
	
	public abstract boolean isLike(Coordinate cd);
	
	
	/**
     * clone method
     */
	public Object clone(){
		try{
			Coordinate cd=(Coordinate)super.clone();
			
			cd.islinear=islinear;
			
			cd.name=name;
			
			return cd;
			
	    }catch(CloneNotSupportedException ex){
		    // this shouldn't happen, since we are Cloneable
		    throw new InternalError();
	    }
	}
}
