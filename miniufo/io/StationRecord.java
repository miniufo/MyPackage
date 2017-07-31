/**
 * @(#)StationRecord.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

/**
 * station record
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class StationRecord implements Cloneable{
	// header
	public int nlev;			// level count include surface variable and level dependent variable
	public int flag;			// flag = 1, flag = 0, 0 means no surface variable follows the header
	
	public float tim;			// -0.5 ~ 0.5, offset to the report time
	public float lat;			// latitude  (degree)
	public float lon;			// longitude (degree)
	
	// data
	public float[] sdata=null;	// surface data
	
	public float[][] ldata=null;	// levels data, ldata.length is z count, ldata[0].length is vcount+1
	
	public String   sid=null;
	
	
	/**
     * constructor
     */
	public StationRecord(){}
	
	/**
     * constructor, construct according to a given StationRecord
     *
     * @param	sr	a given StationRecord
     */
	public StationRecord(StationRecord sr){
		nlev=sr.nlev;
		flag=sr.flag;
		tim =sr.tim;
		lat =sr.lat;
		lon =sr.lon;
		sid =sr.sid;
	}
	
	
	/**
     * used to print out
     */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		sb.append("sid :\t");	sb.append(sid );	sb.append("\n");
		sb.append("tim :\t");	sb.append(tim );	sb.append("\n");
		sb.append("lon :\t");	sb.append(lon );	sb.append("\n");
		sb.append("lat :\t");	sb.append(lat );	sb.append("\n");
		sb.append("nlev:\t");	sb.append(nlev);	sb.append("\n");
		sb.append("flag:\t");	sb.append(flag);	sb.append("\n");
		
		return sb.toString();
	}
	
	/**
     * clone method
     */
	public Object clone(){
		try{
		    StationRecord sr=(StationRecord)super.clone();
			
			sr.nlev=nlev;
			sr.flag=flag;
			sr.tim =tim;
			sr.lat =lat;
			sr.lon =lon;
			sr.sid =sid;
			sr.sdata=sdata.clone();
			sr.ldata=new float[ldata.length][];
			
			for(int i=0;i<ldata.length;i++)
			sr.ldata[i]=ldata[i].clone();
			
			return sr;
			
	    }catch(CloneNotSupportedException ex){
		    // this shouldn't happen, since we are Cloneable
		    throw new InternalError();
	    }
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			float[] a={1,2,3};
			float[] b=a.clone();
			
			b[2]=6;
			
			for(int i=0;i<a.length;i++)
				System.out.println("a["+i+"]: "+a[i]+"\tb["+i+"]: "+b[i]);
			
			System.out.println(a+"\t"+b);
			
		}catch(Exception ex){
			ex.printStackTrace();
			System.exit(0);
		}
	}*/
}
