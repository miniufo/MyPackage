/**
 * @(#)MDate.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.diagnosis;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.TimeZone;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.GregorianCalendar;


/**
 * used to describe the time in kinds of data type
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class MDate implements Cloneable{
	//
	public enum Months{JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC};
	
	private long time;	// YYYYMMDDHHMMSS
	
	private GregorianCalendar gc=new GregorianCalendar(utc);
	
	private static final TimeZone utc=TimeZone.getTimeZone("Etc/UTC");
	
	private static final SimpleDateFormat fmtS=new SimpleDateFormat("yyyy.MM.dd HH:mm:ss");
	private static final SimpleDateFormat fmtL=new SimpleDateFormat("yyyyMMddHHmmss");
	
	private static final Pattern p=Pattern.compile(
		"((\\d{2})?(:(\\d{2}))?z)?(\\d{1,2})?(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)(\\d{1,4})",
	Pattern.CASE_INSENSITIVE);
	
	static{
		synchronized(fmtS){fmtS.setTimeZone(utc);}
		synchronized(fmtL){fmtL.setTimeZone(utc);}
	}
	
	
	/**
     * constructor
     *
     * @param	grads_date_time		a string specified the date. example: hh:mmZddmmmyyyy (For GrADS ctl data type only)
     *
     * @return
     *
     * @exception
     */
	public MDate(String grads_date_time){
		Matcher m=p.matcher(grads_date_time.trim());
		
		String group=null;
		
		if(!m.find()) throw new IllegalArgumentException("grads date format incorrect!");
		
		if((group=m.group(4))==null) gc.set(Calendar.MINUTE,0);
		else gc.set(Calendar.MINUTE,Integer.parseInt(group));
		
		if((group=m.group(2))==null) gc.set(Calendar.HOUR_OF_DAY,0);
		else gc.set(Calendar.HOUR_OF_DAY,Integer.parseInt(group));
		
		if((group=m.group(5))==null) gc.set(Calendar.DATE,1);
		else gc.set(Calendar.DATE,Integer.parseInt(group));
		
		gc.set(Calendar.SECOND,0);
		gc.set(Calendar.MONTH,Months.valueOf(m.group(6).toUpperCase()).ordinal());
		gc.set(Calendar.YEAR,Integer.parseInt(m.group(7)));
		
		toLongTime();
	}
	
	public MDate(long time){
		String t=String.valueOf(time);
		
		if(t.length()==14){
			gc.set(Calendar.YEAR ,Integer.parseInt(t.substring( 0,4 )));
			gc.set(Calendar.MONTH,Integer.parseInt(t.substring( 4,6 ))-1);
			gc.set(Calendar.DATE ,Integer.parseInt(t.substring( 6,8 )));
			gc.set(Calendar.HOUR_OF_DAY,Integer.parseInt(t.substring( 8,10)));
			gc.set(Calendar.MINUTE,Integer.parseInt(t.substring(10,12)));
			gc.set(Calendar.SECOND,Integer.parseInt(t.substring(12,14)));
			
		}else if(t.length()==12){
			gc.set(Calendar.YEAR ,Integer.parseInt(t.substring( 0,4 )));
			gc.set(Calendar.MONTH,Integer.parseInt(t.substring( 4,6 ))-1);
			gc.set(Calendar.DATE ,Integer.parseInt(t.substring( 6,8 )));
			gc.set(Calendar.HOUR_OF_DAY,Integer.parseInt(t.substring( 8,10)));
			gc.set(Calendar.MINUTE,Integer.parseInt(t.substring(10,12)));
			gc.set(Calendar.SECOND,0);
			
		}else if(t.length()==10){
			gc.set(Calendar.YEAR ,Integer.parseInt(t.substring( 0,4 )));
			gc.set(Calendar.MONTH,Integer.parseInt(t.substring( 4,6 ))-1);
			gc.set(Calendar.DATE ,Integer.parseInt(t.substring( 6,8 )));
			gc.set(Calendar.HOUR_OF_DAY,Integer.parseInt(t.substring( 8,10)));
			gc.set(Calendar.MINUTE,0);
			gc.set(Calendar.SECOND,0);
			
		}else if(t.length()==8){
			gc.set(Calendar.YEAR ,Integer.parseInt(t.substring( 0,4 )));
			gc.set(Calendar.MONTH,Integer.parseInt(t.substring( 4,6 ))-1);
			gc.set(Calendar.DATE ,Integer.parseInt(t.substring( 6,8 )));
			gc.set(Calendar.HOUR_OF_DAY,0);
			gc.set(Calendar.MINUTE,0);
			gc.set(Calendar.SECOND,0);
			
		}else if(t.length()==6){
			gc.set(Calendar.YEAR ,Integer.parseInt(t.substring( 0,4 )));
			gc.set(Calendar.MONTH,Integer.parseInt(t.substring( 4,6 ))-1);
			gc.set(Calendar.DATE ,1);
			gc.set(Calendar.HOUR_OF_DAY,0);
			gc.set(Calendar.MINUTE,0);
			gc.set(Calendar.SECOND,0);
			
		}else if(t.length()==4){
			gc.set(Calendar.YEAR ,Integer.parseInt(t.substring( 0,4 )));
			gc.set(Calendar.MONTH,0);
			gc.set(Calendar.DATE ,1);
			gc.set(Calendar.HOUR_OF_DAY,0);
			gc.set(Calendar.MINUTE,0);
			gc.set(Calendar.SECOND,0);
			
		}else throw new IllegalArgumentException("length of time is invalid");
		
		toLongTime();
	}
	
	public MDate(int year,int month,int date,int hour,int minute,int second){
		gc.set(Calendar.YEAR,year);
		gc.set(Calendar.MONTH,month-1);
		gc.set(Calendar.DATE,date);
		gc.set(Calendar.HOUR_OF_DAY,hour);
		gc.set(Calendar.MINUTE,minute);
		gc.set(Calendar.SECOND,second);
		
		toLongTime();
	}
	
	public MDate(int year,int month,int date,int hour,int minute){
		this(year,month,date,hour,minute,0);
		
		toLongTime();
	}
	
	public MDate(int year,int month,int date){
		this(year,month,date,0,0,0);
		
		toLongTime();
	}
	
	public MDate(GregorianCalendar cldr){ cldr.setTimeZone(utc); this.gc=cldr; toLongTime();}
	
	public MDate(){ toLongTime();}
	
	
	/*** getor and setor ***/
	public int getDT(MDate md){
		long result=gc.getTimeInMillis()-md.gc.getTimeInMillis();
		
		result=result>=0?result/1000:-result/1000; // always positive in seconds
		
		int test=(int)result;
		
		if(test!=result) throw new IllegalArgumentException("overflow for too large dt");
		
		return test;
	}
	
	public int getYear()  { return gc.get(Calendar.YEAR);   }
	
	public int getMonth() { return gc.get(Calendar.MONTH)+1;}
	
	public int getDate()  { return gc.get(Calendar.DATE);   }
	
	public int getHour()  { return gc.get(Calendar.HOUR_OF_DAY);}
	
	public int getMinute(){ return gc.get(Calendar.MINUTE); }
	
	public int getSecond(){ return gc.get(Calendar.SECOND); }
	
	public int getDayOfYear(){ return gc.get(Calendar.DAY_OF_YEAR);}
	
	public long getLongTime(){ return time; }
	
	public GregorianCalendar getGregorianCalendar(){ return gc; }
	
	public void setYear(int year)    { gc.set(Calendar.YEAR,year);			toLongTime();}
	
	public void setMonth(int month)  { gc.set(Calendar.MONTH,month-1);		toLongTime();}
	
	public void setDate(int date)    { gc.set(Calendar.DATE,date);			toLongTime();}
	
	public void setHour(int hour)    { gc.set(Calendar.HOUR_OF_DAY,hour);	toLongTime();}
	
	public void setMinute(int minute){ gc.set(Calendar.MINUTE,minute);		toLongTime();}
	
	public void setSecond(int second){ gc.set(Calendar.SECOND,second);		toLongTime();}
	
	public void setDayOfYear(int doy){ gc.set(Calendar.DAY_OF_YEAR,doy);	toLongTime();}
	
	
	/**
     * add method
     *
     * @param	increment	a string specified the increment.
	 *						example: 3HR or 10MN or 1YR (For GrADS ctl data type only)
     *
     * @return	new MDate which is the add result
     */
	public MDate add(String increment){
		MDate newmd=(MDate)this.clone();
		
		int num=0;
		
		String lower_case_incre=increment.toLowerCase();
		
		if(lower_case_incre.indexOf("mn")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("mn",""));
			newmd.gc.add(Calendar.MINUTE,num);
			
		}else if(lower_case_incre.indexOf("hr")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("hr",""));
			newmd.gc.add(Calendar.HOUR,num);
			
		}else if(lower_case_incre.indexOf("dy")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("dy",""));
			newmd.gc.add(Calendar.DATE,num);
			
		}else if(lower_case_incre.indexOf("mo")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("mo",""));
			newmd.gc.add(Calendar.MONTH,num);
			
		}else if(lower_case_incre.indexOf("yr")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("yr",""));
			newmd.gc.add(Calendar.YEAR,num);
			
		}else throw new IllegalArgumentException("unsupported adding String format");
		
		newmd.toLongTime();
		
		return newmd;
	}
	
	/**
     * add method
     *
     * @param	field	which field to be added to
     * @param	amount	how much to add to
     *
     * @return	new MDate which is the add result
     */
	public MDate add(int field,int amount){
		MDate newmd=(MDate)this.clone();
		
		newmd.gc.add(field,amount);	newmd.toLongTime();
		
		return newmd;
	}
	
	public MDate addEq(int field,int amount){
		gc.add(field,amount);
		toLongTime();
		return this;
	}
	
	
	/**
     * whether the MDate is before the given MDate
     *
     * @param	md	given MDate
     *
     * @return	true or not
     */
	public boolean before(MDate md){ return time<md.time;}
	
	/**
     * whether the MDate is after the given MDate
     *
     * @param	md	given MDate
     *
     * @return	true or not
     */
	public boolean after(MDate md) { return time>md.time;}
	
	/**
     * whether the MDate is the same with the given MDate
     *
     * @param	md	given MDate
     *
     * @return	true or not
     */
	public boolean equals(Object md){
		if(time==((MDate)md).time) return true;
		else return false;
	}
	
	
	/**
     * Whether the given year is leap or not
     * A different implementation from GregorianCalendar's
     * Results are the same when year is larger than 1582
     *
     * @param	year	a give year
     *
     * @return	true or not
     */
	public static boolean isLeapYear(int year){
		if(year%4==0&&year%100!=0||year%400==0) return true;
		return false;
	}
	
	
	/**
     * convert to a string in GrADS form
     *
     * @return	a string in GrADS form
     */
	public String toGradsDate(){
		return  String.format("%02d",gc.get(Calendar.HOUR_OF_DAY))+":"+
				String.format("%02d",gc.get(Calendar.SECOND))+"Z"+
				String.format("%02d",gc.get(Calendar.DATE))+
				Months.values()[gc.get(Calendar.MONTH)]+gc.get(Calendar.YEAR);
	}
	
	
	/**
     * used to print out
     */
	public String toString(){
		synchronized(fmtL){
			return fmtS.format(gc.getTime());
		}
	}
	
	/**
     * clone method
     */
	public Object clone(){
		try{
			MDate md=null;	md=(MDate)super.clone();
			md.gc=(GregorianCalendar)gc.clone();
			md.time=time;
			
			return md;
			
	    }catch(CloneNotSupportedException ex){
	    	// this shouldn't happen, since we are Cloneable
		    throw new InternalError();
	    }
	}
	
	
	/*** helper methods ***/
	private void toLongTime(){
		synchronized(fmtL){
			time=Long.parseLong(fmtL.format(gc.getTime()));
		}
	}
	
	
	/** test
	public static void main(String[] arg){
		GregorianCalendar gc=new GregorianCalendar();
		
		for(int year=0;year<=10000;year++){
			if(gc.isLeapYear(year)!=isLeapYear(year))
			System.out.println(year+" is not same. gc: "+gc.isLeapYear(year)+"; MDate: "+isLeapYear(year));
		}
	}*/
}
