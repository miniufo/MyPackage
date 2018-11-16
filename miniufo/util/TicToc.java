/**
 * @(#)TicToc.java	1.0 2013.08.20
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

import java.util.concurrent.TimeUnit;


/**
 * a class to mimic tic-toc commands in Matlab
 *
 * @version 1.0, 2013.08.20
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TicToc{
	//
	private static long savedTime=0;
	
	
	public static void tic(String message){
		System.out.print(message);
		System.out.print("...");
		System.out.flush();
		savedTime=System.nanoTime();
	}
	
	public static void tic(String format,Object... args){
		System.out.printf(format,args);
		System.out.print("...");
		System.out.flush();
		savedTime=System.nanoTime();
	}
	
	public static double toc(TimeUnit unit){
		long elapsedTime=System.nanoTime()-savedTime;
		double re=0;
		
		switch(unit){
		case MILLISECONDS: re=elapsedTime/1e6;		break;
		case SECONDS: re=elapsedTime/1e9;			break;
		case MINUTES: re=elapsedTime/1e9/60.0;		break;
		case HOURS  : re=elapsedTime/1e9/3600.0;	break;
		case DAYS   : re=elapsedTime/1e9/86400.0;	break;
		default: throw new IllegalArgumentException("unsupported TimeUnit: "+unit);
		}
		
		System.out.printf(" (%.3f "+unit+")\n",re);
		
		return re;
	}
	
	public static double toc(String s,TimeUnit unit){
		long elapsedTime=System.nanoTime()-savedTime;
		double re=0;
		
		switch(unit){
		case MILLISECONDS: re=elapsedTime/1e6;		break;
		case SECONDS: re=elapsedTime/1e9;			break;
		case MINUTES: re=elapsedTime/1e9/60.0;		break;
		case HOURS  : re=elapsedTime/1e9/3600.0;	break;
		case DAYS   : re=elapsedTime/1e9/86400.0;	break;
		default: throw new IllegalArgumentException("unsupported TimeUnit: "+unit);
		}
		
		System.out.printf(s+" (%.3f "+unit+")\n",re);
		
		return re;
	}
	
	
	/** test
	public static void main(String[] args){
		tic("counting to 10000000");
		int counter=0;
		for(int i=0;i<10000000;i++) counter++;
		toc();
	}*/
}
