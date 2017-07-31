/**
 * @(#)DrogueOffData.java	1.0 2013.02.28
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;


/**
 * used to describe the MetaData of GDP drifter
 *
 * @version 1.0, 2013.02.28
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DrogueOffData{
	//
	private int id=-1;		// AOML buoy identification number
	
	private long otime=0L;	// original 
	private long atime=0L;	// automatic
	private long mtime=0L;	// manual
	
	
	/**
	 * constructor
	 */
	public DrogueOffData(String oneline){
		String[] tokens=oneline.trim().split("[\\s\\t]+");
		
		id   =Integer.parseInt(tokens[0].trim());
		otime=Long.parseLong(
			tokens[1].substring(0,4)+tokens[1].substring(5,7)+tokens[1].substring(8,10)+"000000"
		);
		atime=Long.parseLong(
			tokens[2].substring(0,4)+tokens[2].substring(5,7)+tokens[2].substring(8,10)+"000000"
		);
		mtime=Long.parseLong(
			tokens[3].substring(0,4)+tokens[3].substring(5,7)+tokens[3].substring(8,10)+"000000"
		);
	}
	
	
	/*** getor and setor ***/
	public int getID(){ return id;}
	
	public long getOriginalTime(){ return otime;}
	
	public long getAutomaticTime(){ return atime;}
	
	public long getManualTime(){ return mtime;}
	
	
	/**
     * used to print out
     */
	public String toString(){
		return String.format("%8d original:%s, automatic:%s, manual:%s",id,otime,atime,mtime);
	}
	
	
	/*** test**
	public static void main(String[] args){
		String s=" 46016	2010-05-24	2004-05-19	  2004-05-21	";
		DrogueOffData md=new DrogueOffData(s);
		System.out.println(md);
	}*/
}
