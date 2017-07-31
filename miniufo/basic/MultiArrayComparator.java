/**
 * @(#)MultiArrayComparator.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.basic;

import java.util.Comparator;


/**
 * a class used to sorted the multiArray
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class MultiArrayComparator implements Comparator<Object>{
	//
	private int keyColumn=0;
	private int sortOrder=1;
	
	public enum sortOrders{ASCEND,DESCEND};
	
	
	/**
	 * constructor
	 */
	public MultiArrayComparator(){}
	
	public MultiArrayComparator(int keyColumn){ this.keyColumn=keyColumn;}
	
	public MultiArrayComparator(int keyColumn,sortOrders so){
		this.keyColumn=keyColumn;
		sortOrder=so.ordinal();
	}
	
	
	public int compare(Object a,Object b){
		if(a instanceof String[])
			return sortOrder*( (String[])a)[keyColumn].compareTo(((String[])b)[keyColumn]);
		else if(a instanceof int[])
			return sortOrder*(((   int[])a)[keyColumn]     -     ((   int[])b)[keyColumn]);
		else if(a instanceof long[])
			return sortOrder*(int)(((long[])a)[keyColumn]  -     ((  long[])b)[keyColumn]);
		else if(a instanceof float[]){
			float tmp=((float[])a)[keyColumn]-((float[])b)[keyColumn];
			
			if(tmp>=0) return sortOrder;
			else return -sortOrder;
			
		}else if(a instanceof double[]){
			double tmp=((double[])a)[keyColumn]-((double[])b)[keyColumn];
			
			if(tmp>=0) return sortOrder;
			else return -sortOrder;
			
		}else return 0;
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			String s="    26    23    23    23   -33  -120 -1000 -1000-32768-32768    22";
			
			long st=System.currentTimeMillis();
			String[] a=split(s,6);
			System.out.println(System.currentTimeMillis()-st);
			for(String t:a) System.out.println(t);
			
		}catch(Exception e){ e.printStackTrace();}
	}*/
}
