/**
 * @(#)Fourier.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.basic.ArrayUtil;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.FDistribution;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.atan2;
import static miniufo.statistics.StatisticsUtil.cArithmeticMean;
import static miniufo.statistics.StatisticsUtil.cLeadLagCorrelationCoefficient;


/**
 * Fourier class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class PowerSpectrum{
	
	/**
	 * prevent from instantiate
     */
	private PowerSpectrum(){}
	
	
	/**
     * power spectral density estimate using fft algorithm
     * 
     * @param	data	data to estimate
     * @param	win		window function
     * @param	Fs		frequency of sampling
     * 
     * @return	re		re[0] is power spectral density, re[1] is periods in unit of 1/Fs
     */
	public static float[][] fftPSDEstimate(float[] data,float[] win,float Fs){
		int len=data.length;
		
		// windowing the data
		float[] windata=new float[len];
		for(int i=0;i<len;i++) windata[i]=data[i]*win[i];
		
		// get scaled spectrum
		float[] spectra=scaledSpectrum(FastFourier.fft(windata),win,Fs,true);
		
		// get periods
		float[] periods=new float[spectra.length];
		for(int i=0,I=periods.length;i<I;i++) periods[i]=len/(i*Fs);
		
		return new float[][]{spectra,periods};
	}
	
	public static float[][] fftPSDEstimate(float[] data){
		return fftPSDEstimate(data,WindowFunction.rectangular(data.length),1);
	}
	
	/**
     * discrete power spectral estimate
     * 
     * @param	data	data to estimate
     * 
     * @return	re		re[0] is power spectrum, re[1] is period
     */
	public static float[][] discretePSDEstimate(float[] data){
		int length=data.length;
		
		DiscreteFourier df=new DiscreteFourier(length);
		
		float[][] co=df.cFourierAB(data);
		
		float[][] re=new float[2][co[0].length];
		
		for(int k=0,K=co[0].length;k<K;k++){
			re[0][k]=(co[0][k]*co[0][k]+co[1][k]*co[1][k])/length;
			re[1][k]=length*1f/k;
		}
		
		scaleOneSide(re[0],length);
		
		return re;
	}
	
	/**
     * continuous power spectral estimate using lead/lag correlation (not checked!!!)
     *
     * @param	data	data to estimate
     * @param	m		max lag-time, if smaller, result will be smoother
     * 
     * @return	re		re[0] is power spectrum, re[1] is period, re[2] is red noise spectrum
     */
	public static float[][] continuousPSDEstimate(float[] data,int m){
		int length=data.length;
		
		if(length<=3) throw new IllegalArgumentException("data are too short for estimate");
		
		float[]   r =new float[m+1];
		float[][] re=new float[3][m+1];
		
		if(m>(float)length/3||m<(float)length/10)
			throw new IllegalArgumentException("m should be larger than n/10 and smaller than n/3");
		
		for(int i=0;i<=m;i++)
		r[i]=cLeadLagCorrelationCoefficient(data,data,i);
		
		float w0=(float)(PI/m);
		for(int k=0;k<=m;k++){
			re[0][k]=r[0];
			
			for(int i=1;i<m;i++) re[0][k]+=r[i]*(float)(1+cos(w0*i))*cos(w0*k*i);
			
			re[0][k]/=m;
		}
		
		re[0][0]*=0.5f;	re[0][m]*=0.5f;
		
		/*** calculate periods ***/
		re[1][0]=Float.MAX_VALUE;
		
		for(int i=1;i<=m;i++) re[1][i]=(float)(m<<1)/i;
		
		/*
		 * red or white noise check
		 * chi-square distribution with degrees of freedom v=(2*n-m/2)/m is used
		 * 95% confidence level is default
		 */
		float v =(2f*length-m/2f)/m;
		float r2=r[1]*r[1];
		float as=cArithmeticMean(re[0]);
		//float co=getChiSquare(0.95f,v)/v;
		float co=(float)new ChiSquaredDistribution(v).inverseCumulativeProbability(0.95)/v;
		
		if(r[1]>0){	// red noise check
			for(int i=0;i<=m;i++)
			re[2][i]=as*(1-r2)/(1+r2-2*r[1]*(float)cos(PI*i/m))*co;
			
		}else{		// white noise check
			for(int i=0;i<=m;i++)
			re[2][i]=as*co;
		}
		
		return re;
	}
	
	
	/**
     * discrete cross spectral denstiy estimate using fft algorithm
     *
     * @param	data1	series 1
     * @param	data2	series 2
     * @param	win		window function
     * @param	Fs		frequency of sampling
     * 
     * @return	re		re[0] is power spectrum, re[1] is period in unit of 1/Fs
     */
	public static float[][] fftCPSDEstimate(float[] data1,float[] data2,float[] win,float Fs){
		int len=data1.length;
		
		if(len!=data2.length) throw new IllegalArgumentException("lengths not equal");
		
		// windowing the data
		float[] windata1=new float[len];
		float[] windata2=new float[len];
		float squaredSum=0;
		for(int i=0;i<len;i++){
			windata1[i]=data1[i]*win[i];
			windata2[i]=data2[i]*win[i];
			squaredSum+=win[i]*win[i];
		}
		
		Complex[] r1=FastFourier.fft(windata1);
		Complex[] r2=FastFourier.fft(windata2);
		
		float[] p1=scaledSpectrum(r1,win,Fs,true);
		float[] p2=scaledSpectrum(r2,win,Fs,true);
		
		float[] C =new float[p1.length];	// coincident spectrum
		float[] Q =new float[p1.length];	// quadrature spectrum
		float[] R =new float[p1.length];	// coherence spectrum
		float[] L =new float[p1.length];	// phase lag spectrum
		float[] cy=new float[p1.length];	// cyclonic spectrum
		float[] ac=new float[p1.length];	// anticyclonic spectrum
		float[] Ts=new float[p1.length];	// periods
		
		for(int i=0,I=p1.length;i<I;i++){
			C[i]=( r1[i].getReal()*r2[i].getReal()+r1[i].getImag()*r2[i].getImag())/squaredSum/Fs;
			Q[i]=(-r1[i].getReal()*r2[i].getImag()+r1[i].getImag()*r2[i].getReal())/squaredSum/Fs;
		}
		
		for(int i=1,I=len%2==0?p1.length-1:p1.length;i<I;i++){ C[i]*=2; Q[i]*=2;}
		
		for(int i=0,I=p1.length;i<I;i++){
			R [i]=(C[i]*C[i]+Q[i]*Q[i])/(p1[i]*p2[i]);
			L [i]=(float)atan2(Q[i],C[i]);
			cy[i]=p1[i]+p2[i]-Q[i];
			ac[i]=p1[i]+p2[i]+Q[i];
			Ts[i]=(float)len/i;
		}
		
		return new float[][]{C,Q,R,L,cy,ac,Ts};
	}
	
	public static float[][] fftCPSDEstimate(float[] data1,float[] data2){
		return fftCPSDEstimate(data1,data2,WindowFunction.rectangular(data1.length),1);
	}
	
	/**
     * discrete cross spectral density estimate
     *
     * @param	data1	series 1
     * @param	data2	series 2
     * 
     * @return	re		re[0] is power spectrum, re[1] is period, re[2] is red noise spectrum
     */
	public static float[][] discreteCPSDEstimate(float[] data1,float[] data2){
		int len=data1.length;
		
		DiscreteFourier df=new DiscreteFourier(len);
		
		float[][] co1=df.cFourierAB(data1);
		float[][] co2=df.cFourierAB(data2);
		
		int length=co1[0].length;
		
		float[] puu=new float[length];	// auto-spectrum of u
		float[] pvv=new float[length];	// auto-spectrum of v
		float[] Cuv=new float[length];	// coincident spectrum
		float[] Quv=new float[length];	// quadrature spectrum
		float[] Ruv=new float[length];	// coherence spectrum
		float[] Luv=new float[length];	// phase lag spectrum
		float[] cy =new float[length];	// cyclonic coherence spectrum
		float[] ac =new float[length];	// anticyclonic coherence spectrum
		float[] Ts =new float[length];	// periods
		
		for(int i=0;i<length;i++){
			puu[i]=(co1[0][i]*co1[0][i]+co1[1][i]*co1[1][i])/len;
			pvv[i]=(co2[0][i]*co2[0][i]+co2[1][i]*co2[1][i])/len;
			Cuv[i]=(co1[0][i]*co2[0][i]+co1[1][i]*co2[1][i])/len;
			Quv[i]=(co1[0][i]*co2[1][i]-co2[0][i]*co1[1][i])/len;
		}
		
		scaleOneSide(puu,len);	scaleOneSide(pvv,len);
		scaleOneSide(Cuv,len);	scaleOneSide(Quv,len);
		
		for(int i=0;i<length;i++){
			Ruv[i]=(Cuv[i]*Cuv[i]+Quv[i]*Quv[i])/(puu[i]*pvv[i]);
			Luv[i]=(float)atan2(Quv[i],Cuv[i]);
			 cy[i]=puu[i]+pvv[i]-Quv[i];
			 ac[i]=puu[i]+pvv[i]+Quv[i];
			 Ts[i]=(float)len/i;
		}
		
		float[][] re=new float[7][];
		
		re[0]=Cuv;	re[1]=Quv;
		re[2]=Ruv;	re[3]=Luv;
		re[4]=cy;	re[5]=ac;
		re[6]=Ts;
		
		return re;
	}
	
	/**
     * continuous cross spectral estimate (not checked!!)
     *
     * @param	data1	series 1
     * @param	data2	series 2
     * @param	m		max lag-time, if smaller, result will be smoother
     * 
     * @return	re		re[0] is power spectrum, re[1] is period, re[2] is red noise spectrum
     */
	public static float[][] continuousCPSDEstimate(float[] data1,float[] data2,int m){
		int length=data1.length;
		
		if(length<=3) throw new IllegalArgumentException("data are too short for estimate");
		
		float[]   r1 =new float[m+1];
		float[]   r2 =new float[m+1];
		float[]   r12=new float[m+1];
		float[]   r21=new float[m+1];

		float[] px =new float[m+1];	// data1 spectrum
		float[] py =new float[m+1];	// data2 spectrum
		float[] pxy=new float[m+1];	// co-spectrum
		float[] qxy=new float[m+1];	// quad-spectrum
		float[] rxy=new float[m+1];	// coherence spectrum
		float[] cxy=new float[m+1];	// phase difference spectrum
		float[] lxy=new float[m+1];	// lag length spectrum
		float[] tim=new float[m+1];	// periodic
		float[] t95=new float[m+1];	// F-test
		float[] s95=new float[m+1];	// Goodman-test
		
		if(m>(float)length/3||m<(float)length/10)
			throw new IllegalArgumentException("m should be larger than n/10 and smaller than n/3");
		
		for(int i=0;i<=m;i++){
			 r1[i]=cLeadLagCorrelationCoefficient(data1,data1,i);
			 r2[i]=cLeadLagCorrelationCoefficient(data2,data2,i);
			r12[i]=cLeadLagCorrelationCoefficient(data1,data2,i);
			r21[i]=cLeadLagCorrelationCoefficient(data2,data1,i);
		}
		
		/*** calculate co-spectrum and quad-spectrum ***/
		float w0=(float)(PI/m);
		for(int k=0;k<=m;k++){
			px[k]=r1[0];	py[k]=r2[0];
			
			for(int i=1;i<m;i++){
				float tmp=(float)(1+cos(w0*i)*cos(w0*k*i));
				
				 px[k]+=r1[i]*tmp;
				 py[k]+=r2[i]*tmp;
				pxy[k]+=(r12[i]+r21[i])*tmp;
				qxy[k]+=(r12[i]-r21[i])*tmp;
			}
			
			 px[k]/=m;
			 py[k]/=m;
			pxy[k]=(r12[0]+pxy[k]/2)/m;
			qxy[k]=(qxy[k]/2)/m;
		}
		
		 px[0]/=2;	 px[m]/=2;
		 py[0]/=2;	 py[m]/=2;
		pxy[0]/=2;	pxy[m]/=2;
		qxy[0]/=2;	qxy[m]/=2;
		
		/*** calculate periods ***/
		tim[0]=Float.MAX_VALUE;
		for(int i=1;i<=m;i++) tim[i]=(float)(m<<1)/i;
		
		/*** calculate coherence spectrum and phase difference spectrum ***/
		for(int i=0;i<=m;i++){
			rxy[i]=(pxy[i]*pxy[i]+qxy[i]*qxy[i])/(px[i]*py[i]);
			cxy[i]=(float)atan2(qxy[i],pxy[i]);
			lxy[i]=cxy[i]*tim[i]/(float)(2*PI);
		}
		
		/*
		 * Significant test for continuous cross spectrum (95% confidence level)
		 * There are two test methods
		 * rxy951: F-test (=sqrt(F/(F+v-1)), where degree of freedom v=(2*n-m/2)/m)
		 *  	   F distribution here is with degrees of freedom numerator 2 and
		 *  	   degrees of freedom denominator v2=2*(v-1)
		 * rxy952: Goodman-test (=sqrt(1-a**(1/(v-1)), where a=0.05 is significant level)
		 *  	   The degrees of freedom v=(2*n-m/2)/m).
		 */
		float v  =(2f*length-m/2f)/m;
		float v2 =2*(v-1);
		//float f95=getF(0.95f,2,v2);
		float f95=(float)new FDistribution(2,v2).inverseCumulativeProbability(0.95);
		
		for(int i=0;i<=m;i++){
			t95[i]=f95/(f95+v-1);
			s95[i]=1-(float)pow(0.05,1f/(v-1));
		}
		
		float[][] re=new float[8][];
		
		re[0]=pxy;	re[1]=qxy;	re[2]=rxy;	re[3]=cxy;
		re[4]=lxy;	re[5]=tim;	re[6]=t95;	re[7]=s95;
		
		return re;
	}
	
	
	/**
     * compute variance-preserving spectrum
     *
     * @param	fftMod		mode of direct result from calling fft function
     * @param	oneSide		one-side or two-side spectrum
     * 
     * @return	re			variance-preserving spectrum
     */
	public static float[] variancePreservingSpectrum(Complex[] fft,float[] win,float Fs,boolean oneSided){
		int len=fft.length;
		
		if(len!=win.length) throw new IllegalArgumentException("lengths not equal");
		
		float winpow=0;								// squared sum of win
		float[] squaredMod=new float[fft.length];	// squared mode of fft
		
		for(int l=0;l<len;l++){
			winpow+=win[l]*win[l];
			squaredMod[l]=fft[l].getSquaredMod()*l/len;	// multiplied by frequency l/len to keep variance preserving
		}
		
		scaleSpectrum(squaredMod,winpow,Fs);
		
		return getSide(squaredMod,oneSided);
	}
	
	/**
     * compute scaled spectrum
     *
     * @param	fftMod		mode of direct result from calling fft function
     * @param	win			window function applied to original data, used to scale the spectrum
     * @param	Fs			frequency of sampling
     * @param	oneSide		one-side or two-side spectrum
     * 
     * @return	re			spectrum scaled by window power, equivalent to pwelch('psd') in matlab
     * 						unit is squared unit of original data per 1/Fs
     */
	public static float[] scaledSpectrum(Complex[] fft,float[] win,float Fs,boolean oneSided){
		int len=fft.length;
		
		if(len!=win.length) throw new IllegalArgumentException("lengths not equal");
		
		float winpow=0;								// squared sum of win
		float[] squaredMod=new float[fft.length];	// squared mode of fft
		
		for(int l=0;l<len;l++){
			winpow+=win[l]*win[l];
			squaredMod[l]=fft[l].getSquaredMod();
		}
		
		scaleSpectrum(squaredMod,winpow,Fs);
		
		return getSide(squaredMod,oneSided);
	}
	
	
	/**
     * generate a realization of series from a given two-sided power spectral density (PSD)
     *
     * @param	psd			a given power spectral density
     * @param	Fs			frequency of sampling
     * 
     * @return	re			one realization of time series
     */
	public static float[] genSeriesFromTwoSidedPSD(float[] psd,float Fs){
		int N=psd.length,half=0;
		boolean even=N%2==0;
		
		if(even) half=N/2-1;
		else half=(N-1)/2;
		
		float[] pha=ArrayUtil.newUniformRandomSeries(half,0,(float)(2.0*Math.PI));
		float[] amp=new float[N];
		
		for(int i=0;i<N;i++) amp[i]=(float)Math.sqrt(psd[i]*N*Fs);
		
		Complex[] whiteHalf =new Complex[half];
		Complex[] whiteNeg  =new Complex[half];
		Complex[] unit      =new Complex[]{new Complex(1,0)};
		Complex[] whiteNoise=null;
		
		for(int i=0;i<half;i++){
			whiteHalf[i       ]=Complex.polar(1,pha[i]);
			whiteNeg [half-i-1]=whiteHalf[i].conjugate(); // whiteNeg = flip(conj(whiteHalf))
		}
		
		if(even) whiteNoise=ArrayUtil.concatAll(Complex.class,unit,whiteHalf,unit,whiteNeg);
		else     whiteNoise=ArrayUtil.concatAll(Complex.class,unit,whiteHalf     ,whiteNeg);
		
		for(int i=0;i<N;i++) whiteNoise[i].multiplyEq(amp[i]);
		
		Complex[] data=FastFourier.ifft(whiteNoise);
		
		float[] re=new float[N];
		
		for(int i=0;i<N;i++) re[i]=data[i].getReal();
		
		return re;
	}
	
	/**
     * generate a realization of series from a given power spectral density (PSD)
     *
     * @param	psd			a given power spectral density
     * @param	Fs			frequency of sampling
     * @param	oneSide		the given PSD is one-sided or two-sided
     * 
     * @return	re			one realization of time series
     */
	public static float[] genSeriesFromPSD(float[] psd,float Fs,boolean oneSide){
		if(!oneSide) return genSeriesFromTwoSidedPSD(psd,Fs);
		else{
			// convert one-sided PSD to two-sided PSD
			float[] tmp=psd.clone();
			
			for(int i=1,I=psd.length;i<I;i++) tmp[i]/=2f;
			
			float[] tmprev=tmp.clone(); ArrayUtil.reverse(tmprev);
			
			float[] psdts=new float[tmp.length*2-1];
			
			System.arraycopy(tmp   ,0,psdts,         0,   tmp.length);
			System.arraycopy(tmprev,0,psdts,tmp.length,tmprev.length-1);
			
			return genSeriesFromTwoSidedPSD(psdts,Fs);
		}
	}
	
	
	/*** help methods ***/
	private static float[] getSide(float[] spectra,boolean oneSided){
		int len=spectra.length;
		int half=oneSided?(len%2==0?len/2+1:(len+1)/2):len;
		
		float[] re=new float[half];
		
		System.arraycopy(spectra,0,re,0,half);
		
		if(oneSided) scaleOneSide(re,len);
		
		return re;
	}
	
	private static void scaleOneSide(float[] oneside,int len){
		int half=len%2==0?len/2+1:(len+1)/2;
		
		if(half!=oneside.length) throw new IllegalArgumentException("invalid length");
		
		if(len%2==0) for(int l=1,L=half-1;l<L;l++) oneside[l]*=2f;
		else for(int l=1;l<half;l++) oneside[l]*=2f;
	}
	
	private static void scaleSpectrum(float[] squaredMod,float squaredSumOfWin,float Fs){
		for(int i=0,I=squaredMod.length;i<I;i++) squaredMod[i]/=squaredSumOfWin*Fs;
	}
	
	
	/** test
	public static void main(String[] args){
		int N=561;
		float Fs=2f;
		float dt=1f/Fs;
		
		boolean oneSide=true;
		int repeat=5000;
		
		if(!oneSide){
			float[] PSD =repeatMean(()->genTwoSidedPSD(dt,100,Fs,N),repeat);
			float[] PSD2=repeatMean(()->getTwoSidedPSD(genSeriesFromPSD(PSD,Fs,oneSide),Fs),repeat);
			
			for(int i=0;i<PSD.length;i++) System.out.println(PSD[i]+"\t"+PSD2[i]);
			
		}else{
			float[] PSD =repeatMean(()->genOneSidedPSD(dt,100,Fs,N),repeat);
			float[] PSD2=repeatMean(()->getOneSidedPSD(genSeriesFromPSD(PSD,Fs,oneSide),Fs),repeat);
			
			for(int i=0;i<PSD2.length;i++) System.out.println(PSD[i]+"\t"+PSD2[i]);
		}
	}
	
	static float[] repeatMean(ArrayGenerator ag,int repeat){
		float[] PSD=ag.generate();
		
		for(int i=1;i<repeat;i++){
			float[] tmp=ag.generate();
			
			for(int l=0,L=PSD.length;l<L;l++) PSD[l]+=tmp[l];
		}
		
		for(int l=0,L=PSD.length;l<L;l++) PSD[l]/=repeat;
		
		return PSD;
	}
	
	static float[] getTwoSidedPSD(float[] data,float Fs){
		Complex[] fft=FastFourier.fft(data);
		
		return scaledSpectrum(fft,WindowFunction.rectangular(data.length),Fs,false);
	}
	
	static float[] getOneSidedPSD(float[] data,float Fs){
		Complex[] fft=FastFourier.fft(data);
		
		return scaledSpectrum(fft,WindowFunction.rectangular(data.length),Fs,true);
	}
	
	static float[] genTwoSidedPSD(float dt,float T,float Fs,int N){
		Complex[] fft=FastFourier.fft(ArrayUtil.newAR1RandomSeries(N,dt,100,10));
		
		return scaledSpectrum(fft,WindowFunction.rectangular(N),Fs,false);
	}
	
	static float[] genOneSidedPSD(float dt,float T,float Fs,int N){
		Complex[] fft=FastFourier.fft(ArrayUtil.newAR1RandomSeries(N,dt,100,10));
		
		return scaledSpectrum(fft,WindowFunction.rectangular(N),Fs,true);
	}*/
}
