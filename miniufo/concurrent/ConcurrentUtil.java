/**
 * @(#)ConcurrentUtil.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.concurrent;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;


/**
 * global concurrent environment
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ConcurrentUtil{
	//
	private static int threadCount=1;
	
	private static ExecutorService executor=null;
	
	
	/*** prevent from initialization ***/
	private ConcurrentUtil(){}
	
	
	/**
	 * initialize the default executor
	 */
	public static void initDefaultExecutor(int n){
		if(executor==null)
			executor=Executors.newFixedThreadPool(n);
		else
			((ThreadPoolExecutor)executor).setCorePoolSize(n);
		
		threadCount=n;
		
		System.out.println("default executor is initialized with "+n+" thread count");
	}
	
	
	/**
	 * get default executor globally given an integer
	 */
	public static ExecutorService defaultExecutor(){
		if(executor==null) initDefaultExecutor(1);
		return executor;
	}
	
	/**
	 * get thread count
	 */
	public static int threadCount(){ return threadCount;}
	
	
	/**
	 * shutdown the default executor
	 */
	public static void shutdown(){
		if(executor!=null) executor.shutdown();
		
		executor=null;
	}
}
