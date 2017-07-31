/**
 * @(#)Binner.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.statistics;

import miniufo.basic.ArrayUtil;


/**
 * used for binning the data
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Binner{
	
	/**
	 * prevent from construction
	 */
	private Binner(){}
	
	
	/**
     * Returns the distribution of data in bins with centers specified by bins.
     * This is similar to hist of Matlab
     *
     * @param	data	data series needs to be binned
     * @param	bins	bins that specified the centers of each bins
     *
     * @return	count	data count within each bins
     */
	public static int[] hist(float[] data,float undef,int length){
		return hist(data,undef,getBins(data,undef,length));
	}
	
	public static int[] hist(float[] data,float undef,float str,float inter,int length){
		return hist(data,undef,getBins(str,inter,length));
	}
	
	public static int[] hist(float[] data,float undef,float[] bins){
		int len=bins.length;
		
		int[] count=new int[len];
		
		for(int l=0,L=data.length;l<L;l++) if(data[l]!=undef){
			boolean isSet=false;
			
			for(int m=1;m<len;m++)
			if(data[l]<bins[m]-(bins[m]-bins[m-1])/2f){ count[m-1]++; isSet=true; break;}
			
			if(!isSet) count[len-1]++;
		}
		
		return count;
	}
	
	public static float[] getBins(float str,float inter,int length){
		float[] bins=new float[length];
		
		for(int i=0;i<length;i++) bins[i]=str+i*inter;
		
		return bins;
	}
	
	public static float[] getBins(float[] data,float undef,int length){
		float[] extremes=ArrayUtil.getExtrema(data,undef);
		
		float interval=(extremes[1]-extremes[0])/length;
		
		float[] bins=new float[length];
		
		bins[0]=extremes[0]+interval/2;
		
		for(int l=1;l<length;l++) bins[l]=bins[l-1]+interval;
		
		return bins;
	}
	
	
	/** test
	public static void main(String[] args){
		float[] u=TextReader.readColumn("d:/Data/GDP/SCS/EddySpeedPDF/dataU.txt",229527,false,1)[0];
		float[] v=TextReader.readColumn("d:/Data/GDP/SCS/EddySpeedPDF/dataV.txt",229527,false,1)[0];
		
		float std=StatisticsUtil.cStandardDeviation(u);
		
		int[] bins=hist(u,-6f*std,std,13);
		
		System.out.println("\n\n");
		
		for(int i:bins) System.out.println(i);
	}*/
}
