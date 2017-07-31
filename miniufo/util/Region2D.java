/**
 * @(#)Region2D.java	1.0 2014.08.12
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;


/**
 * define a 2D region
 *
 * @version 1.0, 2014.08.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Region2D{
	//
	private float lonmin=0;
	private float lonmax=0;
	private float latmin=0;
	private float latmax=0;
	
	private String name =null;
	
	
	/**
	 * constructor
	 * 
	 * @param	lonmin	west-bounded of longitude
	 * @param	lonmax	east-bounded of longitude
	 * @param	latmin	south-bounded of latitude
	 * @param	latmax	north-bounded of latitude
	 */
	public Region2D(float lonmin,float latmin,float lonmax,float latmax,String name){
		if(lonmax<lonmin)
		throw new IllegalArgumentException("lonmax ("+lonmax+") should be larger than lonmin ("+lonmin+")");
		if(latmax<latmin)
		throw new IllegalArgumentException("latmax ("+latmax+") should be larger than latmin ("+latmin+")");
		
		this.lonmin=lonmin;
		this.lonmax=lonmax;
		this.latmin=latmin;
		this.latmax=latmax;
		this.name  =name;
	}
	
	public Region2D(float lonmin,float latmin,float lonmax,float latmax){
		this(lonmin,latmin,lonmax,latmax,"noName");
	}
	
	
	/*** getor and setor ***/
	public float getLonMin(){ return lonmin;}
	
	public float getLonMax(){ return lonmax;}
	
	public float getLatMin(){ return latmin;}
	
	public float getLatMax(){ return latmax;}
	
	public String getName(){ return name;}
	
	
	/**
	 * whether a given point is in the region
	 * 
	 * @param	lon		longitude of a give location
	 * @param	lat		latitude of a give location
	 */
	public boolean inRange(float lon,float lat){ return inXRange(lon)&&inYRange(lat);}
	
	public boolean inXRange(float lon){ return lon>=lonmin&&lon<=lonmax;}
	
	public boolean inYRange(float lat){ return lat>=latmin&&lat<=latmax;}
	
	
	/**
	 * used to print out
	 */
	public String toString(){
		return String.format("%s: Lon[%5.1f ~ %5.1f]; Lat[%5.1f ~ %5.1f]",name,lonmin,lonmax,latmin,latmax);
	}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(new Region(20,-40,120,30,"IO"));
	}*/
}
