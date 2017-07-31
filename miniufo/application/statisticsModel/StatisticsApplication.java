/**
 * @(#)StatisticsApplication.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import static java.lang.Math.sqrt;


/**
 * Basic class for statistics application
 * that takes into account undefined value.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class StatisticsApplication{
	
	/**
     * calculate arithmetic mean
     *
     * @param	data	data in an array
     * @param	undef	undefined value
     *
     * @return	arithmetic mean
     */
	public static float cArithmeticMean(float[] data,float undef){ return cArithmeticMean0(data,undef)[0];}
	
	public static double cArithmeticMean(double[] data,double undef){ return cArithmeticMean0(data,undef)[0];}
	
	
	/**
     * calculate variance
     *
     * @param	data	data in an array
     * @param	undef	undefined value
     *
     * @return	variance
     */
	public static float cVariance(float[] data,boolean weiLengthMinusOne,float undef){ return cVariance0(data,weiLengthMinusOne,undef)[0];}
	
	public static float cVariance(float[] data,float undef){ return cVariance(data,true,undef);}
	
	public static double cVariance(double[] data,boolean weiLengthMinusOne,double undef){ return cVariance0(data,weiLengthMinusOne,undef)[0];}
	
	public static double cVariance(double[] data,double undef){ return cVariance(data,true,undef);}
	
	
	/**
     * calculate range = max - min
     *
     * @param	data	data in an array
     * @param	undef	undefined value
     *
     * @return	variance
     */
	public static float cRange(float[] data,float undef){
		if(data.length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float ei:data)
		if(ei!=undef){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		if(re[0]!=Float.MAX_VALUE) return re[1]-re[0];
		else return undef;
	}
	
	public static double cRange(double[] data,double undef){
		if(data.length==0) throw new IllegalArgumentException("0-length of array");
		
		double[] re=new double[2];
		
		re[0]=Double.MAX_VALUE;
		re[1]=Double.MIN_VALUE;
		
		for(double ei:data)
		if(ei!=undef){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		if(re[0]!=Float.MAX_VALUE) return re[1]-re[0];
		else return undef;
	}
	
	
	/**
     * calculate standard deviation
     *
     * @param	data	data in an array
     * @param	undef	undefined value
     *
     * @return	standard deviation
     */
	public static float cStandardDeviation(float[] data,float undef){
		float var=cVariance(data,undef);
		
		if(var==undef) return undef;
		else return (float)sqrt(var);
	}
	
	public static double cStandardDeviation(double[] data,double undef){
		double var=cVariance(data,undef);
		
		if(var==undef) return undef;
		else return sqrt(var);
	}
	
	
	/**
     * Calculate standard error.
     *
     * @param	data	data in an array
     * @param	undef	undefined value
     *
     * @return	standard error
     */
	public static float cStandardError(float[] data,float undef){
		float[] re=cVariance0(data,true,undef);
		
		float var=re[0],count=re[1];
		
		if(var==undef) return undef;
		else return (float)sqrt(var/count);
	}
	
	public static double cStandardError(double[] data,double undef){
		double[] re=cVariance0(data,true,undef);
		
		double var=re[0],count=re[1];
		
		if(var==undef) return undef;
		else return sqrt(var/count);
	}
	
	
	/**
     * calculate anomaly
     *
     * @param	data	data in an array
     * @param	undef	undefined value
     *
     * @return	anomaly
     */
	public static float[] cAnomaly(float[] data,float undef){
		float[] variance=new float[data.length];
		
		cAnomaly(data,variance,undef);
		
		return variance;
	}
	
	public static double[] cAnomaly(double[] data,double undef){
		double[] variance=new double[data.length];
		
		cAnomaly(data,variance,undef);
		
		return variance;
	}

	public static void cAnomaly(float[] data,float[] anom,float undef){
		if(data.length!=anom.length)
			throw new IllegalArgumentException("array lengths not equal");
		
		float average=cArithmeticMean(data,undef);
		
		for(int i=0,I=data.length;i<I;i++)
		if(data[i]!=undef) anom[i]=(float)(data[i]-average);
	}
	
	public static void cAnomaly(double[] data,double[] anom,double undef){
		if(data.length!=anom.length)
			throw new IllegalArgumentException("array lengths not equal");
		
		double average=cArithmeticMean(data,undef);
		
		for(int i=0,I=data.length;i<I;i++)
		if(data[i]!=undef) anom[i]=(float)(data[i]-average);
	}
	
	
	/**
     * calculate covariance
     *
     * @param	data1	data in an array
     * @param	data2	data in another array
     * @param	undef	undefined value
     *
     * @return	covariance
     */
	public static float cCovariance(float[] data1,float[] data2,float undef){
		if(data1.length!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		float mean_a=cArithmeticMean(data1,undef);
		float mean_b=cArithmeticMean(data2,undef);
		float covariance=0; int count=0;
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data1[i]!=undef&&data2[i]!=undef){ covariance+=(data1[i]-mean_a)*(data2[i]-mean_b); count++;}
		
		return covariance/=count;
	}
	
	public static double cCovariance(double[] data1,double[] data2,double undef){
		if(data1.length!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		double mean_a=cArithmeticMean(data1,undef);
		double mean_b=cArithmeticMean(data2,undef);
		double covariance=0; int count=0;
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data1[i]!=undef&&data2[i]!=undef){ covariance+=(data1[i]-mean_a)*(data2[i]-mean_b); count++;}
		
		return covariance/=count;
	}
	
	
	/**
     * calculate correlation coefficient
     *
     * @param	data1	data in an array
     * @param	data2	data in another array
     * @param	undef	undefined value
     *
     * @return	correlation coefficient
     */
	public static float cCorrelationCoefficient(float[] data1,float[] data2,float undef){
		if(data1.length!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		float average_a=cArithmeticMean(data1,undef);
		float average_b=cArithmeticMean(data2,undef);
		float variance_a=0,variance_b=0,covariance=0;
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data1[i]!=undef){
			float tmp=data1[i]-average_a;
			variance_a+=tmp*tmp;
		}
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data2[i]!=undef){
			float tmp=data2[i]-average_b;
			variance_b+=tmp*tmp;
		}
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data1[i]!=undef&&data2[i]!=undef) covariance+=(data1[i]-average_a)*(data2[i]-average_b);
		
		if(variance_a!=undef&&variance_b!=undef)
			return covariance/(float)sqrt(variance_a*variance_b);
		else
			return undef;
	}
	
	public static double cCorrelationCoefficient(double[] data1,double[] data2,double undef){
		if(data1.length!=data2.length) throw new IllegalArgumentException("array length not equal");
		
		double average_a=cArithmeticMean(data1,undef);
		double average_b=cArithmeticMean(data2,undef);
		double variance_a=0,variance_b=0,covariance=0;
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data1[i]!=undef){
			double tmp=data1[i]-average_a;
			variance_a+=tmp*tmp;
		}
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data2[i]!=undef){
			double tmp=data2[i]-average_b;
			variance_b+=tmp*tmp;
		}
		
		for(int i=0,I=data1.length;i<I;i++)
		if(data1[i]!=undef&&data2[i]!=undef) covariance+=(data1[i]-average_a)*(data2[i]-average_b);
		
		if(variance_a!=undef&&variance_b!=undef)
			return covariance/sqrt(variance_a*variance_b);
		else
			return undef;
	}
	
	
	/**
     * standardize
     *
     * @param	tdata	a given temporal series
     * @param	undef	undefined value
     */
	public static void standardize(float[] tdata,float undef){
		int length=tdata.length,count=0;
		
		float ave=0,var=0,sum=0;
		
		for(int l=0;l<length;l++) if(tdata[l]!=undef){ sum+=tdata[l]; count++;}
		
		if(count>=2){
			ave=sum/count;
			
			for(int l=0;l<length;l++)
			if(tdata[l]!=undef){
				tdata[l]-=ave;
				var+=(tdata[l]*tdata[l]);
			}
			
			var=(float)Math.sqrt(var/(count-1));
			
			for(int l=0;l<length;l++) if(tdata[l]!=undef) tdata[l]/=var;
			
		}else{
			for(int l=0;l<length;l++) tdata[l]=undef;
		}
	}
	
	public static void standardize(double[] tdata,double undef){
		int length=tdata.length,count=0;
		
		double ave=0,var=0,sum=0;
		
		for(int l=0;l<length;l++) if(tdata[l]!=undef){ sum+=tdata[l]; count++;}
		
		if(count>=2){
			ave=sum/count;
			
			for(int l=0;l<length;l++)
			if(tdata[l]!=undef){
				tdata[l]-=ave;
				var+=(tdata[l]*tdata[l]);
			}
			
			var=Math.sqrt(var/(count-1));
			
			for(int l=0;l<length;l++) if(tdata[l]!=undef) tdata[l]/=var;
			
		}else{
			for(int l=0;l<length;l++) tdata[l]=undef;
		}
	}
	
	
	/**
     * anomalization
     *
     * @param	tdata	a given temporal series
     * @param	undef	undefined value
     * 
     * @return	ave		average of the data
     */
	public static float anomalize(float[] tdata,float undef){
		float ave=cArithmeticMean(tdata,undef);
		
		if(ave!=undef){
			for(int l=0,L=tdata.length;l<L;l++) if(tdata[l]!=undef) tdata[l]-=ave;
		}else{
			for(int l=0,L=tdata.length;l<L;l++) tdata[l]=undef;
		}
		
		return ave;
	}
	
	public static double anomalize(double[] tdata,double undef){
		double ave=cArithmeticMean(tdata,undef);
		
		if(ave!=undef){
			for(int l=0,L=tdata.length;l<L;l++) if(tdata[l]!=undef) tdata[l]-=ave;
		}else{
			for(int l=0,L=tdata.length;l<L;l++) tdata[l]=undef;
		}
		
		return ave;
	}
	
	
	/*** helper methods ***/
	private static float[] cArithmeticMean0(float[] data,float undef){
		float tmp_ave=0;
		int count=0;
		
		for(float f:data)
		if(f!=undef){ tmp_ave+=f; count++;}
		
		if(count!=0) return new float[]{tmp_ave/count,count};
		else return new float[]{undef,0};
	}
	
	private static double[] cArithmeticMean0(double[] data,double undef){
		double tmp_ave=0;
		int count=0;
		
		for(double f:data)
		if(f!=undef){ tmp_ave+=f; count++;}
		
		if(count!=0) return new double[]{tmp_ave/count,count};
		else return new double[]{undef,0};
	}
	
	private static float[] cVariance0(float[] data,boolean weiLengthMinusOne,float undef){
		float[] re=cArithmeticMean0(data,undef);
		float variance=0,average=re[0],count=re[1];
		
		for(float f:data)
		if(f!=undef){
			float tmp=f-average;
			variance+=tmp*tmp;
		}
		
		if(count<=1) return new float[]{undef,count};
		
		if(weiLengthMinusOne)
			return new float[]{variance/=count-1,count};
		else
			return new float[]{variance/=count,count};
	}
	
	private static double[] cVariance0(double[] data,boolean weiLengthMinusOne,double undef){
		double[] re=cArithmeticMean0(data,undef);
		double variance=0,average=re[0],count=re[1];
		
		for(double f:data)
		if(f!=undef){
			double tmp=f-average;
			variance+=tmp*tmp;
		}
		
		if(count<=1) return new double[]{undef,count};
		
		if(weiLengthMinusOne)
			return new double[]{variance/=count-1,count};
		else
			return new double[]{variance/=count,count};
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			
			
		}catch(Exception ex){ ex.printStackTrace();}
	}*/
}
