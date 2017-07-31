/**
 * @(#)MetaData.java	1.0 2013.02.27
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;


/**
 * used to describe the MetaData of GDP drifter
 *
 * @version 1.0, 2013.02.27
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class MetaData{
	//
	private int id=-1;		// AOML buoy identification number
	private int wmo=-1;		// World Meteorological Center buoy identification number
	private int expno=-1;	// Experiment number
	private int tdeath=0;	// Type of death:
							// 0: buoy still alive,
							// 1: buoy ran aground,
							// 2: picked up by vessel,
							// 3: stop transmitting,
							// 4: sporadic transmissions,
							// 5: bad batteries,
							// 6: inactive status
	
	private long dtime=0L;	// Deployment time
	private long etime=0L;	// End time
	private long undrg=0L;	// Data Drogue lost time
	
	private float dlon=0;	// Deployment longitude
	private float dlat=0;	// Deployment latitude
	private float elon=0;	// End longitude
	private float elat=0;	// End latitude
	
	private String type=null;	// Buoy type
	
	public static final long undefinedTime=19790101000000L;	// undefined time
	
	
	/**
	 * constructor
	 */
	public MetaData(String oneline){
		int off=0;//oneline.length()==119?0:1;
		
		id   =Integer.parseInt(oneline.substring(0     ,8 +off).trim());
		wmo  =Integer.parseInt(oneline.substring(8 +off,16+off).trim());
		expno=Integer.parseInt(oneline.substring(17+off,22+off).trim());
		type =oneline.substring(23+off,31+off).trim();
		
		dtime=Long.parseLong(
			oneline.substring(32+off,36+off)+oneline.substring(37+off,39+off)+oneline.substring(40+off,42+off)+
			oneline.substring(43+off,45+off)+oneline.substring(46+off,48+off)+"00"
		);
		dlat =Float.parseFloat(oneline.substring(49+off,56+off).trim());
		dlon =Float.parseFloat(oneline.substring(58+off,65+off).trim());
		
		etime=Long.parseLong(
			oneline.substring(66+off,70+off)+oneline.substring(71+off,73+off)+oneline.substring(74+off,76+off)+
			oneline.substring(77+off,79+off)+oneline.substring(80+off,82+off)+"00"
		);
		elat =Float.parseFloat(oneline.substring(83+off,90+off).trim());
		elon =Float.parseFloat(oneline.substring(92+off,99+off).trim());
		
		undrg=Long.parseLong(
			oneline.substring(100+off,104+off)+oneline.substring(105+off,107+off)+oneline.substring(108+off,110+off)+
			oneline.substring(111+off,113+off)+oneline.substring(114+off,116+off)+"00"
		);
		tdeath=Integer.parseInt(oneline.substring(117+off,119+off).trim());
	}
	
	
	/*** getor and setor ***/
	public int getID(){ return id;}
	
	public int getWMO(){ return wmo;}
	
	public int getExperimentNumber(){ return expno;}
	
	public int getTypeOfDeath(){ return tdeath;}
	
	public long getDeployTime(){ return dtime;}
	
	public long getEndTime(){ return etime;}
	
	public long getUndroguedTime(){ return undrg;}
	
	public float[] getDeployLocation(){ return new float[]{dlon,dlat};}
	
	public float[] getEndLocation(){ return new float[]{elon,elat};}
	
	public String getType(){ return type;}
	
	
	/**
     * used to print out
     */
	public String toString(){
		return String.format(
			"%8d%5s, deploy at [%6.2f, %6.2f] on %s, end at [%6.2f, %6.2f] on %s, undrogued on %s",
			id,type,dlon,dlat,dtime,elon,elat,etime,undrg
		);
	}
	
	
	/*** test**
	public static void main(String[] args){
		String s="7704821 32512   129    SVP   1987/07/20 20:09  -15.02   276.94"+
				 " 1988/12/19 14:19  -27.91   257.85 1988/12/19 14:19  3";
		MetaData md=new MetaData(s);
		System.out.println(md);
	}*/
}
