/**
 * @(#)MDate.java	1.0 2014.04.01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.time.temporal.ChronoUnit;
import java.util.Locale;


/**
 * used to describe the time in GrADS
 *
 * @version 1.0, 2014.04.01
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class MDate{
	//
	private long time;	// yyyyMMddHHmmss
	
	private LocalDateTime datetime=null;
	
	private static final DateTimeFormatter grads1=DateTimeFormatter.ofPattern("HH:mm'z'ddMMMyyyy"  ).withLocale(Locale.ENGLISH);
	private static final DateTimeFormatter longF =DateTimeFormatter.ofPattern("yyyyMMddHHmmss"     ).withLocale(Locale.ENGLISH);
	private static final DateTimeFormatter strF  =DateTimeFormatter.ofPattern("yyyy.MM.dd HH:mm:ss").withLocale(Locale.ENGLISH);
	
	
	/**
     * constructor
     *
     * @param	grads_date_time		a string specified the date (For GrADS ctl data type only)
     */
	public MDate(String grads_date_time){
		String input=grads_date_time.toLowerCase();
		
		switch(grads_date_time.length()){
			case 15:
			case 14: break;
			case 12:
			case 11: input=input.substring(0,2)+":00"+input.substring(2); break;
			case  9:
			case  8: input="00:00z"+input; break;
			case  7: input="00:00z01"+input; break;
			default: throw new IllegalArgumentException("grads date format incorrect: "+grads_date_time);
		}
		
		if(input.length()==15){
			String sub=input.substring(8,9); // first character of month
			input=input.replace(sub,sub.toUpperCase());
			
		}else{
			String sub=input.substring(6,8); // first character of month and character of day
			input=input.replace(sub,"0"+sub.toUpperCase());
		}
		
		datetime=LocalDateTime.parse(input,grads1);
		
		toLongTime();
	}
	
	public MDate(long time){
		String t=Long.toString(time);
		
		if(t.length()==14){
			datetime=LocalDateTime.parse(t,longF);
			this.time=time;
			
		}else throw new IllegalArgumentException("length of time is invalid");
	}
	
	public MDate(int year,int month,int date,int hour,int minute,int second){
		datetime=LocalDateTime.of(year,month,date,hour,minute,second);
		
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
	
	public MDate(){
		datetime=LocalDateTime.now();
		toLongTime();
	}
	
	private MDate(LocalDateTime dt){
		datetime=dt;
		toLongTime();
	}
	
	
	/*** getor and setor ***/
	public int getDT(MDate md){ return getDT(md,ChronoUnit.SECONDS);}
	
	public int getDT(MDate md,ChronoUnit unit){
		long result=datetime.until(md.datetime,unit);
		
		int test=(int)result;
		
		if(test!=result) throw new IllegalArgumentException("overflow for too large dt");
		
		return (test<0)?-test:test;
	}
	
	public int getYear()  { return datetime.getYear();}
	
	public int getMonth() { return datetime.getMonthValue();}
	
	public int getDate()  { return datetime.getDayOfMonth();}
	
	public int getHour()  { return datetime.getHour();}
	
	public int getMinute(){ return datetime.getMinute();}
	
	public int getSecond(){ return datetime.getSecond();}
	
	public int getDayOfYear(){ return datetime.getDayOfYear();}
	
	public long getLongTime(){ return time; }
	
	
	/**
     * add method
     *
     * @param	increment	a string specified the increment.
	 *						example: 3hr or 10mn or 1yr (For GrADS ctl increment type only)
     *
     * @return	new MDate which is the add result
     */
	public MDate add(String increment){
		LocalDateTime dt=null;
		
		int num=0;
		
		String lower_case_incre=increment.toLowerCase();
		
		if(lower_case_incre.indexOf("mn")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("mn",""));
			dt=datetime.plusMinutes(num);
			
		}else if(lower_case_incre.indexOf("hr")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("hr",""));
			dt=datetime.plusHours(num);
			
		}else if(lower_case_incre.indexOf("dy")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("dy",""));
			dt=datetime.plusDays(num);
			
		}else if(lower_case_incre.indexOf("mo")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("mo",""));
			dt=datetime.plusMonths(num);
			
		}else if(lower_case_incre.indexOf("yr")!=-1){
			num=Integer.parseInt(lower_case_incre.replace("yr",""));
			dt=datetime.plusYears(num);
			
		}else throw new IllegalArgumentException("unsupported adding String format");
		
		return new MDate(dt);
	}
	
	public MDate addYears(int amount){ return new MDate(datetime.plusYears(amount));}
	
	public MDate addMonths(int amount){ return new MDate(datetime.plusMonths(amount));}
	
	public MDate addDays(int amount){ return new MDate(datetime.plusDays(amount));}
	
	public MDate addHours(int amount){ return new MDate(datetime.plusHours(amount));}
	
	public MDate addMinutes(int amount){ return new MDate(datetime.plusMinutes(amount));}
	
	public MDate addSeconds(int amount){ return new MDate(datetime.plusSeconds(amount));}
	
	
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
     * convert to a string in GrADS form
     *
     * @return	a string in GrADS form
     */
	public String toGradsDate(){ return datetime.format(grads1);}
	
	public String toFormattedDate(DateTimeFormatter formatter){ return datetime.format(formatter);}
	
	/**
     * used to print out
     */
	public String toString(){ return datetime.format(strF);}
	
	
	/**
	 * return a copy
	 */
	public MDate copy(){ return new MDate(datetime);}
	
	
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
	
	
	/*** helper methods ***/
	private void toLongTime(){ time=Long.parseLong(datetime.format(longF));}
	
	
	/** test
	public static void main(String[] arg){
		int times=200000;
		String datetime="00Z04JUL2008";
		
		miniufo.test.diagnosis.MDateJoda m2=new miniufo.test.diagnosis.MDateJoda(datetime);
		MDate m1=new MDate(datetime);
		miniufo.test.diagnosis.MDate m3=new miniufo.test.diagnosis.MDate(datetime);
		
		miniufo.util.TicToc.tic("Testing MDate");
		for(int i=0;i<times;i++){
			MDate tmp=new MDate(m1.add("1dy").toGradsDate());
			m1=tmp;
		}
		miniufo.util.TicToc.toc();
		System.out.println(m1.toGradsDate());
		
		miniufo.util.TicToc.tic("Testing MDateJoda");
		for(int i=0;i<times;i++){
			miniufo.test.diagnosis.MDateJoda tmp=new miniufo.test.diagnosis.MDateJoda(m2.add("1dy").toGradsDate());
			m2=tmp;
		}
		miniufo.util.TicToc.toc();
		System.out.println(m2.toGradsDate());
		
		miniufo.util.TicToc.tic("Testing old MDate");
		for(int i=0;i<times;i++){
			miniufo.test.diagnosis.MDate tmp=new miniufo.test.diagnosis.MDate(m3.add("1dy").toGradsDate());
			m3=tmp;
		}
		miniufo.util.TicToc.toc();
		System.out.println(m3.toGradsDate());
	}*/
}
