/**
 * @(#)ClimateStabilityTest.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.statistics;

import miniufo.statistics.StatisticsUtil;
import static java.lang.Math.sqrt;


/**
 * used to test climate stability
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ClimateStabilityTest extends StatisticsUtil{
	
	/**
     * get result of U-Test
     *
     * @param	arr_data	data in an array
     * @param	ave_exp		average as expected
     * @param	var_exp		standard deviation as expected
     *
     * @return	result
     */
	public static float cUTest(float[] data,float ave_exp,float var_exp){
		return (cArithmeticMean(data)-ave_exp)/var_exp*(float)sqrt(data.length);
	}
	
	public static double cUTest(double[] data,double ave_exp,double var_exp){
		return (cArithmeticMean(data)-ave_exp)/var_exp*sqrt(data.length);
	}
	
	
	/**
     * get result of U-Test
     *
     * @param	data1	data 1 in an array
     * @param	data2	data 2 in an array
     * @param	var1	variance 1 as expected
     * @param	var2	variance 2 as expected
     *
     * @return	result
     */
	public static float cUTest(float[] data1,float[] data2,float var1,float var2){
		return (cArithmeticMean(data1)-cArithmeticMean(data2))/(float)sqrt(var1/data1.length+var2/data2.length);
	}
	
	public static double cUTest(double[] data1,double[] data2,double var1,double var2){
		return (cArithmeticMean(data1)-cArithmeticMean(data2))/sqrt(var1/data1.length+var2/data2.length);
	}
	
	
	/**
     * get result of T-Test
     *
     * @param	data	data in an array
     * @param	ave_exp		average as expected
     *
     * @return	result
     */
	public static float cTTest(float[] data,float ave_exp){
		return (cArithmeticMean(data)-ave_exp)/cStandardDeviation(data)*(float)sqrt(data.length);
	}
	
	public static double cTTest(double[] data,double ave_exp){
		return (cArithmeticMean(data)-ave_exp)/cStandardDeviation(data)*sqrt(data.length);
	}
	
	
	/**
     * get result of T-Test
     *
     * @param	data1	data 1 in an array
     * @param	data2	data 2 in an array
     *
     * @return	result
     */
	public static float cTTest(float[] data1,float[] data2){
		return (cArithmeticMean(data1)-cArithmeticMean(data2))/
		(float)sqrt(cVariance(data1)/data1.length+cVariance(data2)/data2.length);
	}
	
	public static double cTTest(double[] data1,double[] data2){
		return (cArithmeticMean(data1)-cArithmeticMean(data2))/
		sqrt(cVariance(data1)/data1.length+cVariance(data2)/data2.length);
	}
	
	
	/**
     * get result of X2-Test
     *
     * @param	data	data in an array
     * @param	var_exp		variance as expected
     *
     * @return	result
     */
	public static float cXTest(float[] data,float var_exp){
		return (data.length-1)*cVariance(data)/var_exp;
	}
	
	public static double cXTest(double[] data,double var_exp){
		return (data.length-1)*cVariance(data)/var_exp;
	}
	
	
	/**
     * get result of X2-Test
     *
     * @param	data	data in an array
     * @param	ave_exp		average as expected
     * @param	var_exp		standard deviation as expected
     *
     * @return	result
     */
	public static float cXTest(float[] data,float ave_exp,float var_exp){
		float x=0,tmp=0;
		
		for(int i=0;i<data.length;i++){
			tmp=(data[i]-ave_exp)/var_exp;
			x+=tmp*tmp;
		}
		
		System.out.println(x);
		
		return (data.length-1)*cVariance(data)/var_exp;
	}
	
	public static double cXTest(double[] data,double ave_exp,double var_exp){
		double x=0,tmp=0;
		
		for(int i=0;i<data.length;i++){
			tmp=(data[i]-ave_exp)/var_exp;
			x+=tmp*tmp;
		}
		
		System.out.println(x);
		
		return (data.length-1)*cVariance(data)/var_exp;
	}
	
	
	/**
     * get result of F-Test
     *
     * @param	data1	data 1 in an array
     * @param	data2	data 2 in an array
     *
     * @return	result
     */
	public static float cFTest(float[] data1,float[] data2){
		return (data1.length*cVariance(data1)/(data1.length-1))/(data2.length*cVariance(data2)/(data2.length-1));
	}
	
	public static double cFTest(double[] data1,double[] data2){
		return (data1.length*cVariance(data1)/(data1.length-1))/(data2.length*cVariance(data2)/(data2.length-1));
	}
}
