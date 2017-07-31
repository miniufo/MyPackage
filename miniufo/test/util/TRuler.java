/**
 * @(#)TRuler.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;

import java.util.ArrayList;
import miniufo.descriptor.DataDescriptor.DataType;
import miniufo.diagnosis.MDate;


/**
 * used to get array tags which meet some condition
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TRuler{
	//
	private int length;
	
	private long[] times=null;
	
	
	/**
     * constructor
     */
	public TRuler(MDate stime,int length,DataType type){
		this.length=length;
		
		times=new long[length];
		
		switch(type){
			case ANNUAL:
				for(int l=0;l<length;l++)
				times[l]=stime.addYears(l).getLongTime();
				break;
			case MONTHLY:
				for(int l=0;l<length;l++)
				times[l]=stime.addMonths(l).getLongTime();
				break;
			case DAILY:
				for(int l=0;l<length;l++)
				times[l]=stime.addDays(l).getLongTime();
				break;
			case DAILY4:
				for(int l=0;l<length;l++)
				times[l]=stime.addHours(6*l).getLongTime();
				break;
			default:
				throw new IllegalArgumentException("unsupported type: "+type);
		}
	}
	
	
	/**
     * get tags whose fields between start time and end time
     * i.e., get tags from 1998.1.1 to 1998.1.31
     *
     * @param	ss		start time
     * @param	ee		end time
     *
     * @return	tags	tags for arrays
     */
	public int[] getTags(String ss,String ee){
		String[] tmp1=ss.split("\\.");
		String[] tmp2=ee.split("\\.");
		
		MDate st=null;
		MDate et=null;
		
		if(tmp1.length==4&&tmp2.length==4){
			st=new MDate(
				Integer.parseInt(tmp1[0]),
				Integer.parseInt(tmp1[1]),
				Integer.parseInt(tmp1[2]),
				Integer.parseInt(tmp1[3]),0
			);
			
			et=new MDate(
				Integer.parseInt(tmp2[0]),
				Integer.parseInt(tmp2[1]),
				Integer.parseInt(tmp2[2]),
				Integer.parseInt(tmp2[3]),0
			);
			
		}else if(tmp1.length==3&&tmp2.length==3){
			st=new MDate(
				Integer.parseInt(tmp1[0]),
				Integer.parseInt(tmp1[1]),
				Integer.parseInt(tmp1[2])
			);
			
			et=new MDate(
				Integer.parseInt(tmp2[0]),
				Integer.parseInt(tmp2[1]),
				Integer.parseInt(tmp2[2])
			);
			
		}else if(tmp1.length==2&&tmp2.length==2){
			st=new MDate(Integer.parseInt(tmp1[0]),Integer.parseInt(tmp1[1]),1);
			
			et=new MDate(Integer.parseInt(tmp2[0]),Integer.parseInt(tmp2[1]),1);
			
		}else if(tmp1.length==1&&tmp2.length==1){
			st=new MDate(Integer.parseInt(tmp1[0]),1,1);
			
			et=new MDate(Integer.parseInt(tmp2[0]),1,1);
			
		}else throw new IllegalArgumentException(
			"Invalid stime String. Like this according to the data: 1998.1.1 or 1998.1 or 1998"
		);
		
		if(!st.before(et))
		throw new IllegalArgumentException("start time should be before the end time");
		
		long strlong=st.getLongTime();
		long endlong=et.getLongTime();
		
		ArrayList<Integer> al=new ArrayList<Integer>();
		
		for(int l=0;l<length;l++)
		if(times[l]>=strlong){
			if(times[l]<=endlong) al.add(l);
			else break;
		}
		
		int[] tags=new int[al.size()];
		
		for(int l=0,L=al.size();l<L;l++) tags[l]=al.get(l);
		
		return tags;
	}
	
	/**
     * get tags whose fields meet the condition
     *
     * @param	conditions	forms (no whitespace): year(1998,1998);month(2,3);day(1,1);hour(0,0)
     *
     * @return	tags	tags for arrays
     */
	public int[] getTags(String conditions){
		String[] conds=conditions.split(";");
		
		int[] yrRange=new int[2]; yrRange[0]=Integer.MIN_VALUE; yrRange[1]=Integer.MAX_VALUE;
		int[] moRange=new int[2]; moRange[0]=Integer.MIN_VALUE; moRange[1]=Integer.MAX_VALUE;
		int[] dyRange=new int[2]; dyRange[0]=Integer.MIN_VALUE; dyRange[1]=Integer.MAX_VALUE;
		int[] hrRange=new int[2]; hrRange[0]=Integer.MIN_VALUE; hrRange[1]=Integer.MAX_VALUE;
		
		for(String condition:conds){
			String[] tmp=condition.split("\\(");
			
			if(tmp.length!=2) throw new IllegalArgumentException(
				"invalid condition form\n" +
				"like this: year(1998,2000);month(2,3);day(5,10);hour(0,0)"
			);
			
			String[] vv=tmp[1].split(",");
			
			if(vv.length!=2) throw new IllegalArgumentException(
				"invalid condition form\n" +
				"like this: year(1998,2000);month(2,3);day(5,10);hour(0,0)"
			);
			
			if(tmp[0].equalsIgnoreCase("year")){
				yrRange[0]=Integer.parseInt(vv[0].substring(0,vv[0].length()));
				yrRange[1]=Integer.parseInt(vv[1].substring(0,vv[1].length()-1));
				
			}else if(tmp[0].equalsIgnoreCase("month")){
				moRange[0]=Integer.parseInt(vv[0].substring(0,vv[0].length()));
				moRange[1]=Integer.parseInt(vv[1].substring(0,vv[1].length()-1));
				
			}else if(tmp[0].equalsIgnoreCase("day")){
				dyRange[0]=Integer.parseInt(vv[0].substring(0,vv[0].length()));
				dyRange[1]=Integer.parseInt(vv[1].substring(0,vv[1].length()-1));
				
			}else if(tmp[0].equalsIgnoreCase("hour")){
				hrRange[0]=Integer.parseInt(vv[0].substring(0,vv[0].length()));
				hrRange[1]=Integer.parseInt(vv[1].substring(0,vv[1].length()-1));
			}
		}
		
		ArrayList<Integer> al=new ArrayList<Integer>();
		
		for(int l=0;l<length;l++){
			String tim=String.valueOf(times[l]);
			
			int yr=Integer.parseInt(tim.substring(0, 4));
			int mo=Integer.parseInt(tim.substring(4, 6));
			int dy=Integer.parseInt(tim.substring(6, 8));
			int hr=Integer.parseInt(tim.substring(8,10));
			
			if(
				yr>=yrRange[0]&&yr<=yrRange[1]&&
				mo>=moRange[0]&&mo<=moRange[1]&&
				dy>=dyRange[0]&&dy<=dyRange[1]&&
				hr>=hrRange[0]&&hr<=hrRange[1]
				
			) al.add(l);
		}
		
		int[] tags=new int[al.size()];
		
		for(int l=0,L=al.size();l<L;l++) tags[l]=al.get(l);
		
		return tags;
	}
	
	
	/** test
	public static void main(String[] args){
		TRuler tr=new TRuler(new MDate(1987,1,1),36524,DataType.DAILY4);
		
		String conditions="month(2,2);day(29,29);hour(0,0)";
		
		System.out.println(java.util.Arrays.toString(tr.getTags(conditions)));
	}*/
}
