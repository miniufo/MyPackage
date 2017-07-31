/**
 * @(#)GrADSMapRecord.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

/**
 * Map record for GrADS
 *
 * @version 1.0, 02/01/2007
 * @author  MiniLynn
 * @since   MDK1.0
 */
public final class GrADSMapRecord implements Cloneable{
	// header
	public int recsty=1;		// record style : recsty = 1, means data record
	public int linsty=0;		//  line  style : linsty = 0, means coastline;
								// 				  linsty = 1, means national boundary;
								// 				  linsty = 2, means river;
								// 				  linsty = 3, means province boundary;
	
	public int num;				// total grid number : <= 255
	
	// data
	public float[] lon=null;	// latitude  (degree)
	public float[] lat=null;	// longitude (degree)
	
	
	/**
     * constructor
     */
	public GrADSMapRecord(int num){
		this.num=num;
		lon=new float[num];
		lat=new float[num];
	}
	
	
	/**
     * constructor, construct according to a given GradsMapRecord
     *
     * @param	gmr	a given GradsMapRecord
     */
	public GrADSMapRecord(GrADSMapRecord gmr){
		recsty=gmr.recsty;
		linsty=gmr.linsty;
		num =gmr.num;
		lon =new float[num];	System.arraycopy(gmr.lon,0,lon,0,num);
		lat =new float[num];	System.arraycopy(gmr.lat,0,lat,0,num);
	}
	
	
	/**
     * clone method
     */
	public Object clone(){
		try{
		    GrADSMapRecord gmr=(GrADSMapRecord)super.clone();
			
		    gmr.recsty=recsty;
		    gmr.linsty=linsty;
		    gmr.num =num;
		    gmr.lon =new float[num];	System.arraycopy(lon,0,gmr.lon,0,num);
		    gmr.lat =new float[num];	System.arraycopy(lat,0,gmr.lat,0,num);
			
			return gmr;
			
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
