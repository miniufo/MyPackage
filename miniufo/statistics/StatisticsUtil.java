/**
 * @(#)StatisticsUtiljava	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.statistics;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;


/**
 * used for basic statistic method
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class StatisticsUtil{
	
	/**
     * calculate arithmetic mean
     *
     * @param	data	data in an array
     *
     * @return	arithmetic mean
     */
	public static float sum(float[] data){
		float sum=0;
		
		for(float ei:data) sum+=ei;
		
		return sum;
	}
	
	public static double sum(double[] data){
		double sum=0;
		
		for(double ei:data) sum+=ei;
		
		return sum;
	}
	
	
	/**
     * calculate arithmetic mean
     *
     * @param	data	data in an array
     *
     * @return	arithmetic mean
     */
	public static float cArithmeticMean(float[] data){
		float ave=0;
		
		for(float ei:data) ave+=ei;
		
		return ave/data.length;
	}
	
	public static double cArithmeticMean(double[] data){
		double ave=0;
		
		for(double ei:data) ave+=ei;
		
		return ave/data.length;
	}
	
	
	/**
     * calculate variance weighted by n-1
     *
     * @param	data	data in an array
     *
     * @return	variance (n-1 degree of freedom)
     */
	public static float cVariance(float[] data){ return cVariance(data,true);}
	
	public static float cVariance(float[] data,boolean weiLengthMinusOne){
		float average=cArithmeticMean(data),variance=0;
		
		for(float ei:data){
			float tmp=ei-average;
			variance+=tmp*tmp;
		}
		
		if(data.length<2) throw new IllegalArgumentException("data length too short");
		
		if(weiLengthMinusOne)
			return variance/=data.length-1;
		else
			return variance/=data.length;
	}
	
	public static double cVariance(double[] data){ return cVariance(data,true);}
	
	public static double cVariance(double[] data,boolean weiLengthMinusOne){
		double average=cArithmeticMean(data),variance=0;
		
		for(double ei:data){
			double tmp=ei-average;
			variance+=tmp*tmp;
		}
		
		if(data.length<2) throw new IllegalArgumentException("data length too short");
		
		if(weiLengthMinusOne)
			return variance/=data.length-1;
		else
			return variance/=data.length;
	}
	
	
	/**
     * calculate standard deviation of a series
     *
     * @param	data	data in an array
     *
     * @return	standard deviation (n-1 degree of freedom)
     */
	public static float cStandardDeviation(float[] data){ return (float)sqrt(cVariance(data));}
	
	public static float cStandardDeviation(float[] data,boolean weiLengthMinusOne){
		return (float)sqrt(cVariance(data,weiLengthMinusOne));
	}
	
	public static double cStandardDeviation(double[] data){ return sqrt(cVariance(data));}
	
	public static double cStandardDeviation(double[] data,boolean weiLengthMinusOne){
		return sqrt(cVariance(data,weiLengthMinusOne));
	}
	
	
	/**
     * calculate standard error of the mean of a series, Se = std / sqrt(n)
     *
     * @param	data	data in an array
     *
     * @return	standard error
     */
	public static float cStandardError(float[] data){
		return cStandardDeviation(data)/(float)Math.sqrt(data.length);
	}
	
	public static double cStandardError(double[] data){
		return cStandardDeviation(data)/Math.sqrt(data.length);
	}
	
	
	/**
     * calculate skewness coefficient
     *
     * @param	data	data in an array
     *
     * @return	skewness coefficient
     */
	public static float cSkewnessCoefficient(float[] data){
		float sum=0;
		float ave=cArithmeticMean(data);
		float var=cStandardDeviation(data);
		
		for(int i=0,I=data.length;i<I;i++) sum+=pow((data[i]-ave)/var,3);
		
		return sum*(float)sqrt(1.0/(6*data.length));
	}
	
	public static double cSkewnessCoefficient(double[] data){
		double sum=0;
		double ave=cArithmeticMean(data);
		double var=cStandardDeviation(data);
		
		for(int i=0,I=data.length;i<I;i++) sum+=pow((data[i]-ave)/var,3);
		
		return sum*sqrt(1.0/(6*data.length));
	}
	
	
	/**
     * calculate kurtosis coefficient
     *
     * @param	data	data in an array
     *
     * @return	kurtosis coefficient
     */
	public static float cKurtosisCoefficient(float[] data){
		int len=data.length;
		
		float sum=0;
		float ave=cArithmeticMean(data);
		float var=cStandardDeviation(data);
		
		for(int i=0;i<len;i++) sum+=pow((data[i]-ave)/var,4);
		
		return (sum/len-3)*(float)sqrt(len/24f);
	}
	
	public static double cKurtosisCoefficient(double[] data){
		int len=data.length;
		
		double sum=0;
		double ave=cArithmeticMean(data);
		double var=cStandardDeviation(data);
		
		for(int i=0;i<len;i++) sum+=pow((data[i]-ave)/var,4);
		
		return (sum/len-3)*sqrt(len/24.0);
	}
	
	
	/**
     * Calculate root-mean-square error of two series
     *
     * @param	data1	data 1 in an array
     * @param	data2	data 2 in an array
     */
	public static float cRMSE(float[] data1,float[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		float sum=0;
		
		for(int i=0;i<len;i++){
			float diff=data1[i]-data2[i];
			sum+=diff*diff;
		}
		
		sum/=len;
		
		return (float)sqrt(sum);
	}
	
	public static double cRMSE(double[] data1,double[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		double sum=0;
		
		for(int i=0;i<len;i++){
			double diff=data1[i]-data2[i];
			sum+=diff*diff;
		}
		
		sum/=len;
		
		return sqrt(sum);
	}
	
	
	/**
     * Calculate 1-norm, summing absolute differences of between two arrays.
     *
     * @param	data1	data 1 in an array
     * @param	data2	data 2 in an array
     */
	public static float c1Norm(float[] data1,float[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		float sum=0;
		
		for(int i=0;i<len;i++) sum+=Math.abs(data1[i]-data2[i]);
		
		return (float)sum;
	}
	
	public static double c1Norm(double[] data1,double[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		double sum=0;
		
		for(int i=0;i<len;i++) sum+=Math.abs(data1[i]-data2[i]);
		
		return sum;
	}
	
	
	/**
     * calculate anomaly
     *
     * @param	data	data in an array
     *
     * @return	anomaly
     */
	public static float[] cAnomaly(float[] data){
		float average=cArithmeticMean(data);
		
		float[] variance=new float[data.length];
		
		for(int i=0,I=data.length;i<I;i++) variance[i]=data[i]-average;
		
		return variance;
	}
	
	public static double[] cAnomaly(double[] data){
		double average=cArithmeticMean(data);
		
		double[] variance=new double[data.length];
		
		for(int i=0,I=data.length;i<I;i++) variance[i]=data[i]-average;
		
		return variance;
	}
	
	
	/**
     * calculate covariance
     *
     * @param	data1	data in an array
     * @param	data2	data in another array
     *
     * @return	covariance
     */
	public static float cCovariance(float[] data1,float[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		float mean_a=cArithmeticMean(data1);
		float mean_b=cArithmeticMean(data2);
		float covariance=0;
		
		for(int i=0;i<len;i++) covariance+=(data1[i]-mean_a)*(data2[i]-mean_b);
		
		return covariance/=(len-1);
	}
	
	public static double cCovariance(double[] data1,double[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		double mean_a=cArithmeticMean(data1);
		double mean_b=cArithmeticMean(data2);
		double covariance=0;
		
		for(int i=0;i<len;i++) covariance+=(data1[i]-mean_a)*(data2[i]-mean_b);
		
		return covariance/=(len-1);
	}
	
	
	/**
     * calculate Euclidean distance between two vectors in N-dim space
     * i.e., ||A-B||=sqrt((a1-b1)^2+(a2-b2)^2+...+(an-bn)^2)
     *
     * @param	data1	data in an array
     * @param	data2	data in another array
     *
     * @return	dis		Euclidean distance in N-dim space
     */
	public static float cEuclideanDistance(float[] data1,float[] data2){
		if(data1.length!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		float dis=0;
		
		for(int i=0,I=data1.length;i<I;i++){
			float diff=data1[i]-data2[i];
			dis+=diff*diff;
		}
		
		return (float)sqrt(dis);
	}
	
	public static double cEuclideanDistance(double[] data1,double[] data2){
		if(data1.length!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		double dis=0;
		
		for(int i=0,I=data1.length;i<I;i++){
			double diff=data1[i]-data2[i];
			dis+=diff*diff;
		}
		
		return sqrt(dis);
	}
	
	
	/**
     * calculate correlation coefficient
     *
     * @param	data1	data in an array
     * @param	data2	data in another array
     *
     * @return	correlation coefficient
     */
	public static float cCorrelationCoefficient(float[] data1,float[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		float mean1=cArithmeticMean(data1);
		float mean2=cArithmeticMean(data2);
		float var1=0,var2=0;
		
		for(int i=0;i<len;i++){
			float tmp=data1[i]-mean1;
			var1+=tmp*tmp;
		}
		
		for(int i=0;i<len;i++){
			float tmp=data2[i]-mean2;
			var2+=tmp*tmp;
		}
		
		return cCorrelationCoefficient(data1,data2,mean1,mean2,var1,var2);
	}
	
	public static double cCorrelationCoefficient(double[] data1,double[] data2){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		double mean1=cArithmeticMean(data1);
		double mean2=cArithmeticMean(data2);
		double var1=0,var2=0;
		
		for(int i=0;i<len;i++){
			double tmp=data1[i]-mean1;
			var1+=tmp*tmp;
		}
		
		for(int i=0;i<len;i++){
			double tmp=data2[i]-mean2;
			var2+=tmp*tmp;
		}
		
		return cCorrelationCoefficient(data1,data2,mean1,mean2,var1,var2);
	}
	
	
	/**
     * calculate partial correlation coefficient r12_3
     *
     * @param	data1	time series
     * @param	data2	time series
     * @param	data3	time series
     *
     * @return	partial correlation coefficient
     */
	public static float cPartialCorrelationCoefficient(float[] data1,float[] data2,float[] data3){
		if(data1.length!=data2.length||data2.length!=data3.length)
			throw new IllegalArgumentException("array length not equal");
		
		float r12=cCorrelationCoefficient(data1,data2);
		float r13=cCorrelationCoefficient(data1,data3);
		float r23=cCorrelationCoefficient(data2,data3);
		
		float r12_3=(r12-r13*r23)/(float)sqrt((1-r13*r13)*(1-r23*r23));
		
		return r12_3;
	}
	
	public static double cPartialCorrelationCoefficient(double[] data1,double[] data2,double[] data3){
		if(data1.length!=data2.length||data2.length!=data3.length)
			throw new IllegalArgumentException("array length not equal");
		
		double r12=cCorrelationCoefficient(data1,data2);
		double r13=cCorrelationCoefficient(data1,data3);
		double r23=cCorrelationCoefficient(data2,data3);
		
		double r12_3=(r12-r13*r23)/sqrt((1-r13*r13)*(1-r23*r23));
		
		return r12_3;
	}
	
	
	/**
     * calculate unbiased correlation coefficient
     *
     * @param	data1	data in an array
     * @param	data2	data in another array
     *
     * @return	unbiased correlation coefficient
     */
	public static float cUnbiasedCorrelationCoefficient(float[] data1,float[] data2){
		float r=cCorrelationCoefficient(data1,data2);
		
		return r*(1+(1-r*r)/(2*(data1.length-4)));
	}
	
	public static double cUnbiasedCorrelationCoefficient(double[] data1,double[] data2){
		double r=cCorrelationCoefficient(data1,data2);
		
		return r*(1+(1-r*r)/(2*(data1.length-4)));
	}
	
	
	/**
     * calculate lead-lag correlation coefficient,
	 * when data2 is delayed by the given value 'delay' to data1
     *
     * @param	data1	data in an array
     * @param	data2	data in another array
     * @param	delay	time count that data2 is delayed to data1 (always positive)
     *
     * @return	correlation coefficient
     */
	public static float cLeadLagCorrelationCoefficient(float[] data1,float[] data2,int delay){
		int len=data1.length;
		
		if(len!=data2.length)
			throw new IllegalArgumentException("array length not equal");
		if(delay<0||delay>len-2)
			throw new IllegalArgumentException("invalid delay, should be in [0, "+(len-2)+"]");
		
		float var1=0,mean1=cArithmeticMean(data1);
		float var2=0,mean2=cArithmeticMean(data2);
		
		for(int i=0;i<len;i++){
			float tmp1=data1[i]-mean1;	var1+=tmp1*tmp1;
			float tmp2=data2[i]-mean2;	var2+=tmp2*tmp2;
		}
		
		var1*=(float)(len-delay)/len;
		var2*=(float)(len-delay)/len;
		
		float[] tmp1=new float[len];
		float[] tmp2=new float[len];
		
		for(int i=0;i<len;i++){
			tmp1[i]=mean1;
			tmp2[i]=mean2;
		}
		
		System.arraycopy(data1,delay,tmp1,0,len-delay);
		System.arraycopy(data2,    0,tmp2,0,len-delay);
		
		return cCorrelationCoefficient(tmp1,tmp2,mean1,mean2,var1,var2);
	}
	
	public static double cLeadLagCorrelationCoefficient(double[] data1,double[] data2,int delay){
		int len=data1.length;
		
		if(len!=data2.length)
			throw new IllegalArgumentException("array length not equal");
		if(delay<0||delay>len-2)
			throw new IllegalArgumentException("invalid delay, should be in [0, "+(len-2)+"]");
		
		double var1=0,mean1=cArithmeticMean(data1);
		double var2=0,mean2=cArithmeticMean(data2);
		
		for(int i=0;i<len;i++){
			double tmp1=data1[i]-mean1;	var1+=tmp1*tmp1;
			double tmp2=data2[i]-mean2;	var2+=tmp2*tmp2;
		}
		
		var1*=(double)(len-delay)/len;
		var2*=(double)(len-delay)/len;
		
		double[] tmp1=new double[len];
		double[] tmp2=new double[len];
		
		for(int i=0;i<len;i++){
			tmp1[i]=mean1;
			tmp2[i]=mean2;
		}
		
		System.arraycopy(data1,delay,tmp1,0,len-delay);
		System.arraycopy(data2,    0,tmp2,0,len-delay);
		
		return cCorrelationCoefficient(tmp1,tmp2,mean1,mean2,var1,var2);
	}
	
	
	/**
     * standardize
     *
     * @param	tdata	a given temporal series
     */
	public static void standardize(float[] tdata){
		int len=tdata.length;
		
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		float ave=cArithmeticMean(tdata);
		float var=0;
		
		for(int l=0;l<len;l++){
			tdata[l]-=ave;
			var+=(tdata[l]*tdata[l]);
		}
		
		var=(float)Math.sqrt(var/(len-1));
		
		for(int l=0;l<len;l++) tdata[l]/=var;
	}
	
	public static void standardize(double[] tdata){
		int len=tdata.length;
		
		if(len<2) throw new IllegalArgumentException("data length too short");
		
		double ave=cArithmeticMean(tdata);
		double var=0;
		
		for(int l=0;l<len;l++){
			tdata[l]-=ave;
			var+=(tdata[l]*tdata[l]);
		}
		
		var=Math.sqrt(var/(len-1));
		
		for(int l=0;l<len;l++) tdata[l]/=var;
	}
	
	
	/**
     * anomalize
     *
     * @param	tdata	a given temporal series
     */
	public static float anomalize(float[] tdata){
		float ave=cArithmeticMean(tdata);
		
		for(int l=0,L=tdata.length;l<L;l++) tdata[l]-=ave;
		
		return ave;
	}
	
	public static double anomalize(double[] tdata){
		double ave=cArithmeticMean(tdata);
		
		for(int l=0,L=tdata.length;l<L;l++) tdata[l]-=ave;
		
		return ave;
	}
	
	
	/**
     * T-Test of significance of mean value (Ref. Fengying Wei)
     *
     * @param	data1	a given series
     * @param	data2	a given series
     * 
     * @return	T variable for test
     */
	public static float TTest(float[] data1,float[] data2){
		int n1=data1.length;
		int n2=data2.length;
		
		float m1=cArithmeticMean(data1);
		float m2=cArithmeticMean(data2);
		
		float s2=(cVariance(data1)*(n1-1f)+cVariance(data2)*(n2-1f))/(n1+n2-2f);
		
		if(s2==0) return 0;
		else return (m1-m2)/(float)Math.sqrt(s2*(1f/n1+1f/n2));
	}
	
	
	/*** helper methods ***/
	/**
	 * compute correlation coefficient (r) by specifying the means and accumulated variances
	 * 
     * @param	data1	a given series
     * @param	data2	a given series
     * @param	mean1	mean of data1
     * @param	mean2	mean of data2
     * @param	var1	accumulated variance of data1
     * @param	var2	accumulated variance of data2
	 */
	private static float cCorrelationCoefficient
	(float[] data1,float[] data2,float mean1,float mean2,float var1,float var2){
		float covariance=0;
		
		for(int i=0,I=data1.length;i<I;i++) covariance+=(data1[i]-mean1)*(data2[i]-mean2);
		
		return covariance/(float)sqrt(var1*var2);
	}
	
	private static double cCorrelationCoefficient
	(double[] data1,double[] data2,double mean1,double mean2,double var1,double var2){
		double covariance=0;
		
		for(int i=0,I=data1.length;i<I;i++) covariance+=(data1[i]-mean1)*(data2[i]-mean2);
		
		return covariance/sqrt(var1*var2);
	}
	
	
	/** test
	public static void main(String[] args){
		final int len=160;
		final int delay=150;
		
		float[] a=new float[len];
		float[] b=new float[len];
		
		for(int i=0;i<len;i++){
			a[i]=(float)Math.random();
			b[i]=(float)Math.random();
		}
		
		float co1=cLeadLagCorrelationCoefficient(a,b,delay);
		float co2=cLeadLagCorrelationCoefficient2(a,b,delay);
		
		System.out.println(co1);
		System.out.println(co2);
		System.out.println(co1==co2);
	}*/
}
