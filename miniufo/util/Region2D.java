/**
 * @(#)Region2D.java	1.0 2014.08.12
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;


/**
 * Define a 2D region.
 *
 * @version 1.0, 2014.08.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Region2D{
	//
	private float xmin=0;
	private float xmax=0;
	private float ymin=0;
	private float ymax=0;
	
	private String name=null;
	
	
	/**
	 * Constructor.
	 * 
	 * @param	xmin	west  bounded
	 * @param	ymin	south bounded
	 * @param	xmax	east  bounded
	 * @param	ymax	north bounded
	 */
	public Region2D(float xmin,float ymin,float xmax,float ymax,String name){
		if(xmax<xmin)
		throw new IllegalArgumentException("xmax ("+xmax+") should be larger than xmin ("+xmin+")");
		if(ymax<ymin)
		throw new IllegalArgumentException("ymax ("+ymax+") should be larger than ymin ("+ymin+")");
		
		this.xmin=xmin;
		this.xmax=xmax;
		this.ymin=ymin;
		this.ymax=ymax;
		this.name  =name;
	}
	
	public Region2D(float xmin,float ymin,float xmax,float ymax){ this(xmin,ymin,xmax,ymax,"noName");}
	
	
	/*** getor and setor ***/
	public float  getXMin(){ return xmin;}
	
	public float  getXMax(){ return xmax;}
	
	public float  getYMin(){ return ymin;}
	
	public float  getYMax(){ return ymax;}
	
	public String getName(){ return name;}
	
	
	/**
	 * whether a given point is in the region
	 * 
	 * @param	x	x-coordinate of a give location
	 * @param	y	y-coordinate of a give location
	 */
	public boolean inRange(float x,float y){ return inXRange(x)&&inYRange(y);}
	
	public boolean inXRange(float x){ return x>=xmin&&x<=xmax;}
	
	public boolean inYRange(float y){ return y>=ymin&&y<=ymax;}
	
	
	/**
	 * used to print out
	 */
	public String toString(){
		return String.format("%s: X[%8.3f ~ %8.3f]; Y[%8.3f ~ %8.3f]",name,xmin,xmax,ymin,ymax);
	}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(new Region(20,-40,120,30,"IO"));
	}*/
}
