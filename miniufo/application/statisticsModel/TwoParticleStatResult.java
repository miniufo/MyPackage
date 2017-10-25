/**
 * @(#)TwoParticleStatResult.java	1.0 2017.08.28
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.io.FileWriter;
import java.io.IOException;


/**
 * Result of two-particle statistical results.
 *
 * @version 1.0, 2017.08.28
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TwoParticleStatResult{
	//
	int   dt=0;			// delta-T in seconds
	int tNum=0;			// maximum time lag
	
	int numOfPairs=0;	// total number of pairs at initial time
	
	  int[] num=null;	// number of pairs as a function of time
	
	float[] Dxx=null;	// x-dispersion as a function of time
	float[] Dyy=null;	// y-dispersion as a function of time
	float[] Dis=null;	// total dispersion as a function of time
	float[] Kxx=null;	// x-diffusivity as a function of time
	float[] Kyy=null;	// y-diffusivity as a function of time
	
	
	/**
	 * Constructor.
	 * 
	 * @param	tRad	maximum time lag
	 */
	public TwoParticleStatResult(int tNum){
		if(tNum<1) throw new IllegalArgumentException("tRad should be positive");
		
		this.tNum=tNum;
		
		num=new int  [tNum];
		Dxx=new float[tNum];
		Dyy=new float[tNum];
		Dis=new float[tNum];
		Kxx=new float[tNum];
		Kyy=new float[tNum];
	}
	
	
	/**
	 * write out to a specific file
	 */
	public void toFile(String fname){
		try(FileWriter fw=new FileWriter(fname)){
			fw.write(toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	/**
	 * used to print out
	 */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		sb.append(String.format(" initial pairs: %8d\n",num[0]));
		
		sb.append(" time(day)  samples   Dxx(km)   Dyy(km)   Dis(km)   Kxx(10^7cm^2/s)  Kyy(10^7cm^2/s)\n");
		
		for(int l=0,L=num.length;l<L;l++)
		sb.append(String.format(
			"   %5.1f    %6d   %7.3f   %7.3f   %7.3f      %7.3f          %7.3f\n",
			l*dt/86400f,num[l],Dxx[l]/1e6,Dyy[l]/1e6,Dis[l]/1e6,Kxx[l],Kyy[l]
		));
		
		return sb.toString();
	}
	
	
	/** test
	public static void main(String arg[]){
		System.out.println("[235,23]".replaceAll("\\[|\\]",""));
	}*/
}
