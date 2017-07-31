/**
 * @(#)Fitter.java	1.0 2013.03.07
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.basic.ArrayUtil;


/**
 * fitter used to fit a series of data
 *
 * @version 1.0, 2013.03.07
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class Fitter{
	//
	protected int len=0;		// data length
	
	protected float[] x=null;	// array of points where function are evaluated
	
	
	/**
     * constructor
     *
     * @param	x	array of points where function are evaluated
     */
	public Fitter(float[] x){
		len=x.length;
		this.x=x;
		
		checkLength(len);
	}
	
	public Fitter(int len){
		this.len=len;
		
		checkLength(len);
		
		x=ArrayUtil.newMonotonousArray(len,1);
	}
	
	
	/**
     * fitting
     *
     * @param	y	a series of observations
     */
	public abstract void fit(float[] y);
	
	
	/**
     * calculate values at given points x
     *
     * @param	x	a given array of locations to be evaluated at
     */
	public abstract float[] cValues(float[] x);
	
	public abstract float[] cValues();
	
	
	/*** helper methods ***/
	private void checkLength(int len){
		if(len<2) throw new IllegalArgumentException("no enough data");
	}
}
