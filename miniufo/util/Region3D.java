/**
 * @(#)Region3D.java	1.0 2017.05.12
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;


/**
 * define a 3D region
 *
 * @version 1.0, 2017.05.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Region3D{
	//
	private float lonmin=0;
	private float lonmax=0;
	private float latmin=0;
	private float latmax=0;
	private float levmin=0;
	private float levmax=0;
	
	private String name =null;
	
	
	/**
	 * constructor
	 * 
	 * @param	lonmin	west-bound of longitude
	 * @param	lonmax	east-bound of longitude
	 * @param	latmin	south-bound of latitude
	 * @param	latmax	north-bound of latitude
	 * @param	levmin	lower-bound of level
	 * @param	levmax	upper-bound of level
	 */
	public Region3D(float lonmin,float latmin,float levmin,float lonmax,float latmax,float levmax,String name){
		if(lonmax<lonmin)
		throw new IllegalArgumentException("lonmax ("+lonmax+") should be larger than lonmin ("+lonmin+")");
		if(latmax<latmin)
		throw new IllegalArgumentException("latmax ("+latmax+") should be larger than latmin ("+latmin+")");
		if(levmax<levmin)
		throw new IllegalArgumentException("levmax ("+levmax+") should be larger than levmin ("+levmin+")");
		
		this.lonmin=lonmin;
		this.lonmax=lonmax;
		this.latmin=latmin;
		this.latmax=latmax;
		this.levmin=levmin;
		this.levmax=levmax;
		this.name  =name;
	}
	
	public Region3D(float lonmin,float latmin,float levmin,float lonmax,float latmax,float levmax){
		this(lonmin,latmin,levmin,lonmax,latmax,levmax,"noName");
	}
	
	
	/*** getor and setor ***/
	public float getLonMin(){ return lonmin;}
	
	public float getLonMax(){ return lonmax;}
	
	public float getLatMin(){ return latmin;}
	
	public float getLatMax(){ return latmax;}
	
	public float getLevMin(){ return levmin;}
	
	public float getLevMax(){ return levmax;}
	
	public String getName(){ return name;}
	
	
	/**
	 * whether a given point is in the region
	 * 
	 * @param	lon		longitude of a give location
	 * @param	lat		latitude of a give location
	 */
	public boolean inRange(float lon,float lat,float lev){ return inXRange(lon)&&inYRange(lat)&&inZRange(lev);}
	
	public boolean inXRange(float lon){ return lon>=lonmin&&lon<=lonmax;}
	
	public boolean inYRange(float lat){ return lat>=latmin&&lat<=latmax;}
	
	public boolean inZRange(float lev){ return lev>=levmin&&lev<=levmax;}
	
	
	/**
	 * used to print out
	 */
	public String toString(){
		return String.format("%s: Lon[%5.1f ~ %5.1f]; Lat[%5.1f ~ %5.1f]; Lev[%5.1f ~ %5.1f]",name,lonmin,lonmax,latmin,latmax,levmin,levmax);
	}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(new Region(20,-40,120,30,"IO"));
	}*/
}
