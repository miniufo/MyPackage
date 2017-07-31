/**
 * @(#)SignificanceTest.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.statistics;

import org.apache.commons.math3.distribution.TDistribution;


/**
 * used for significance test
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class SignificanceTest{
	
	/**
     * test the significance of correlation coefficient
     * using T-Distribution
     *
     * @param	corr	correlation coefficient
     * @param	P		confidence level, e.g., 0.95
     * @param	V		degree of freedom (samples)
     *
     * @return	bool	whether the result is significant
     */
	public static boolean testCorrelationCoefficient(float corr,float P,float V){
		if(corr>=getCriticalCorrelation(P,V)) return true;
		
		return false;
	}
	
	public static float getCriticalCorrelation(float P,float V){
		TDistribution td=new TDistribution(V);
		
		float x=0;
		
		x=(float)Math.abs(td.inverseCumulativeProbability((1-P)/2));
		
		return (float)Math.sqrt(x*x/(V-2+x*x));
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			System.out.println(getCriticalCorrelation(0.99f,21));
	    	
	    }catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
	}*/
}
