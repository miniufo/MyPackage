/**
 * @(#)Bootstrap.java	1.0 2013.11.27
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.statistics;

import java.lang.reflect.Array;
import java.util.Random;


/**
 * used for basic statistic method
 *
 * @version 1.0, 2013.11.27
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Bootstrap{
	//
	private static final Random rnd=new Random();
	
	
	/**
	 * constructor
	 */
	private Bootstrap(){}
	
	
	/**
	 * resample without allocating memory
	 */
	@SuppressWarnings("unchecked")
	public static <T> T[] resample(Class<T> ty,T[] samples){
		int size=samples.length;
		
		T[] resamples=(T[])Array.newInstance(ty,size);
		
		int[] idx=resampleIndex(size);
		
		for(int i=0;i<size;i++) resamples[i]=samples[idx[i]];
		
		return resamples;
	}
	
	
	/**
	 * resample the indices of original samples
	 */
	public static int[] resampleIndex(int size){
		int[] idx=new int[size];
		
		for(int i=0;i<size;i++) idx[i]=rnd.nextInt(size);
		
		return idx;
	}
	
	
	/** test
	public static void main(String[] args){
		
	}*/
}
