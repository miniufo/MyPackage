/**
 * @(#)GaussMarkovEstimator.java	1.0 2014.08.08
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.diagnosis.MDate;


/**
 * An estimator for the harmonic cycles as well as spatial variations
 * References:
 *   Lumpkin 2003, GRL
 * 	 Lumpkin and Johnson 2013, JGR
 *
 * @version 1.0, 2014.08.08
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class GaussMarkovEstimator{
	//
	protected boolean hasST=false;	// has spatial terms
	
	protected int len=0;
	
	protected double mean=0;		// mean of the data
	protected double dmn =0;		// delta from orginal bin mean, new mean = dmn + mean
	protected double var =0;		// variance of the data
	
	protected float[] disx=null;	// zonal distance from bin center (degree)
	protected float[] disy=null;	// meridional distance from bin center (degree)
	
	protected double[] times=null;	// times in double format (unit of year)
	protected double[] fs   =null;	// frequencies of cycle, including cycle-0 i.e., mean
	
	public enum AutoCorrType{T1st,T2nd,TCosExp}
	
	
	/**
	 * set arguments of frequencies and TL
	 * 
	 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
	 * @param	TL		Lagrangian integral timescale (unit: year)
	 */
	public void setFrequenciesAndTimescales(AutoCorrType act,float[] freqs,float... TSs){
		if(freqs.length<1)
		throw new IllegalArgumentException("number of cycles should be at least zero");
		
		if(TSs==null||TSs.length<1)
		throw new IllegalArgumentException("at least one timescale TL");
		
		if(TSs[0]<=0)
		throw new IllegalArgumentException("TL should be positive");
		
		for(int i=0;i<freqs.length;i++)
		if(freqs[i]<=0) throw new IllegalArgumentException("frequencies should be positive");
		
		fs=new double[freqs.length+1];
		
		fs[0]=0;	// time-invariant mean
		
		for(int i=1,I=freqs.length+1;i<I;i++) fs[i]=freqs[i-1];
		
		assignRxx();
		assignRnn(act,TSs);
		assignA();
	}
	
	/**
	 * Gauss-Markov estimation of cycles
	 */
	public abstract void estimateCycles(boolean computeError);
	
	
	/**
	 * reconstruct the mean and cycles by A * cycles
	 * 
	 * @return	re	matrix form of cycles: [xc1 xc2 ... xs1 xs2 ... x0]'
	 */
	public abstract float[] reconstructMeanCycles(long[] ltims,float[] dx,float[] dy);
	
	public abstract float[] reconstructMeanCycles(long[] ltims);
	
	/**
	 * compute the variance contributions of each component
	 * 
	 * @return	re	for fs = {0, 1, 2},
	 * 				[0] the mean squares of GM mean,
	 * 				[1] the mean squares of annual cycle,
	 * 				[2] the mean squares of semiannual cycle,
	 * 				[3] the mean squares of spatial variations,
	 * 				[4] the mean squares of residuals,
	 */
	public abstract float[] varianceContribution();
	
	
	/**
	 * reject results of decomposition according to error
	 */
	public abstract boolean rejectResidual(float threshold);
	
	public abstract boolean rejectAmplitude(float threshold);
	
	
	/*** getor and setor ***/
	
	/**
	 * get residual data after removing the mean and fs signals
	 * 
	 * @param	numOfCyc	number of cycles to be removed, started from 0, 
	 * 						2 corresponding to 0cycle (mean) + 1cycle + 2cycle (per year)
	 */
	public abstract float[] getResidualData(int numOfCyc);
	
	/**
	 * compute amplitudes [mean, fs1, fs2, ...]
	 */
	public float[] getCycleAmplitudes(){
		checkResult();
		
		float[] co=getCoefficients(true);
		float[] re=new float[fs.length];
		
		re[0]=co[0];
		
		for(int i=1,I=re.length;i<I;i++) re[i]=(float)Math.hypot(co[i],co[i+fs.length-1]);
		
		return re;
	}
	
	/**
	 * compute sin initial phase [fs1, fs2, ...]
	 */
	public float[] getCycleSinPhase(){
		checkResult();
		
		float[] co=getCoefficients(false);
		float[] re=new float[fs.length-1];
		
		for(int i=0,I=re.length;i<I;i++) re[i]=(float)Math.atan2(co[i+1],co[i+fs.length]);
		
		return re;
	}
	
	/**
	 * compute cos initial phase [fs1, fs2, ...]
	 */
	public float[] getCycleCosPhase(){
		checkResult();
		
		float[] co=getCoefficients(false);
		float[] re=new float[fs.length-1];
		
		for(int i=0,I=re.length;i<I;i++) re[i]=(float)(Math.atan2(co[i+1],co[i+fs.length])-Math.PI/2);
		
		return re;
	}
	
	/**
	 * compute estimated error for the unknowns (cycles)
	 * <=1 is OK while >1 is rejected
	 */
	public abstract float[] getCycleErrors();
	
	/**
	 * get coefficients of unknowns in float array
	 * 
	 * @return	re	matrix form of cycles: [x0 xc1 xc2 ... xs1 xs2 ... dx dx^2 dy dy^2 dxdy]'
	 */
	public abstract float[] getCoefficients(boolean addDataMean);
	
	/**
	 * set data without changing other parameters
	 */
	public abstract void setData(float[] data);
	
	
	/*** helper methods ***/
	protected abstract void assignRxx();
	
	protected abstract void assignRnn(AutoCorrType act,float... timescales);
	
	protected abstract void assignA();
	
	protected abstract void checkResult();
	
	
	protected static double cTd(double TL){
		if(TL<=0)
		throw new IllegalArgumentException("integral timescale should be positive");
		
		return Math.sqrt(Math.PI/2.0)/Math.exp(-0.5)*TL;
	}
	
	protected static double cTe(double Td){
		return 2.0*Math.sqrt(2.0)/Math.PI*Td;
	}
	
	protected static double cAutoCorr(AutoCorrType act,float dT,float... timescales){
		double ac=0;
		
		switch(act){
			case T1st: ac=cAutoCorr1stOrder(dT,timescales[0]); break;
			case T2nd: ac=cAutoCorr2ndOrder(dT,timescales[0],timescales[1]); break;
			case TCosExp: ac=cAutoCorrCosExp(dT,timescales[0]); break;
			default: throw new IllegalArgumentException("unsupported AutoCorrType: "+act);
		}
		
		return ac;
	}
	
	protected static double cAutoCorr1stOrder(double tlag,double TL){
		return Math.exp(-Math.abs(tlag)/TL);
	}
	
	protected static double cAutoCorr2ndOrder(double tlag,double Tv,double Ta){
		double ratio=Ta/Tv;
		return (Math.exp(-Math.abs(tlag)/Tv)-ratio*Math.exp(-Math.abs(tlag)/Ta))/(1-ratio);
	}
	
	protected static double cAutoCorrCosExp(double tlag,double TL){
		double Td=cTd(TL);
		double Te=cTe(Td);
		double tmp=tlag/Te;
		
		return Math.cos(Math.PI*tlag/(2.0*Td))*Math.exp(-tmp*tmp);
	}
	
	protected static double toDoubleTime(long time){
		MDate md=new MDate(time);
		
		double yr=md.getYear();
		double dy=md.getDayOfYear();
		double hr=md.getHour();
		double mn=md.getMinute();
		double se=md.getSecond();
		
		boolean leap=MDate.isLeapYear(md.getYear());
		
		return yr+(dy+(hr+(mn+se/60.0)/60.0)/24.0)/(leap?366.0:365.0);
	}
	
	protected static double cMeanSquares(double[] data){
		double re=0;
		
		for(double d:data) re+=d*d;
		
		return re/=data.length-1;
	}
	
	protected static double[] toDoubleTimes(long[] lt){
		double[] rt=new double[lt.length];
		
		for(int l=0,L=lt.length;l<L;l++) rt[l]=toDoubleTime(lt[l]);
		
		return rt;
	}
	
	
	/** test
	public static void main(String[] args){
		int len=1000;
		
		float[] y =new float[len];
		long[] t =new long[len];
		float[] dx=new float[len];
		float[] dy=new float[len];
		
		MDate str=new MDate();
		
		for(int i=0;i<len;i++){
			t [i]=str.addDays(i).getLongTime();
			dx[i]=(float)(Math.random()-0.5);
			dy[i]=(float)(Math.random()-0.5);
			y [i]=(float)(1.0+
				2.0*sin(2.0*PI*t[i]/10000000000.0+0.1)+1.5*sin(2.0*PI*t[i]/10000000000.0+1.3)+
				0.4*dx[i]*dx[i]+0.6*dx[i]+0.2*dy[i]*dy[i]+0.3*dy[i]+0.1*dx[i]*dy[i]
			);
		}
		
		GaussMarkovEstimator2Blas gme=new GaussMarkovEstimator2Blas(y,t,dx,dy);
		gme.setFrequenciesAndTL(new float[]{1,2},0.0137f);
		gme.estimateCycles(false);
		
		float[] re=gme.varianceContribution();
		
		for(int l=0;l<re.length;l++) System.out.println(re[l]);
	}*/
}
