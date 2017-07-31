/**
 * @(#)Record.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.io.Serializable;


/**
 * used to describe a record of a Lagrangian model at a specific time
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Record implements Serializable{
	//
	private static final long serialVersionUID = 1239922442535904638L;
	
	private int cycleNum=-1;	// designed for Argo
	
	private long  time=0;		// time
	
	private float lon =0;		// longitude
	private float lat =0;		// latitude
	
	private float[] data=null;	// values for attached data (e.g., uvel,vvel,temp,salt)
	
	public static final float undef=-9999.0f;	// default undefined value
	
	
	/**
	 * constructor
	 */
	public Record(long time){
		this.time=time;
		
		data=new float[2];
	}
	
	public Record(long time,float lon,float lat){
		this.time=time;
		this.lon =lon;
		this.lat =lat;
		
		data=new float[2];
	}
	
	public Record(long time,float lon,float lat,float uvel,float vvel){
		this.time=time;
		this.lon =lon;
		this.lat =lat;
		
		data=new float[2];
		data[0]=uvel;
		data[1]=vvel;
	}
	
	public Record(long time,float lon,float lat,int dataLen){
		this(time,lon,lat);
		
		data=new float[dataLen];
		
		for(int i=0;i<dataLen;i++) data[i]=undef;
	}
	
	public Record(Record r){
		this.time=r.time;
		this.lon =r.lon;
		this.lat =r.lat;
		
		if(r.data!=null) data=r.data.clone();
	}
	
	
	/*** getor and setor ***/
	public int getCycleNum(){ return cycleNum;}
	
	public int getDataLength(){ return data.length;}
	
	public long getTime(){ return time;}
	
	public float getLon(){ return lon;}
	
	public float getLat(){ return lat;}
	
	public float getDataValue(int idx){ return data[idx];}
	
	public float[] getDataValues(){ return data;}
	
	public void setCycleNum(int cycleNum){ this.cycleNum=cycleNum;}
	
	public void setLon(float lon){ this.lon=lon;}
	
	public void setLat(float lat){ this.lat=lat;}
	
	public void setTime(long time){ this.time=time;}
	
	public void setData(int idx,float value){ data[idx]=value;}
	
	public void addData(int idx,float value){ data[idx]+=value;}
	
	
	/**
     * used to print out
     */
	public String toString(){
		return String.format("%9.3f%10.3f%10.3f%10.3f%19s",lon,lat,data[0],data[1],time);
	}
}
