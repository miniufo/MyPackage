/**
 * @(#)FilterModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.statistics;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import miniufo.mathsphysics.Complex;
import miniufo.mathsphysics.FastFourier;
import miniufo.mathsphysics.HarmonicFitter;
import miniufo.mathsphysics.MathsPhysicsUtil;
import miniufo.mathsphysics.PolynomialFitter;
import miniufo.mathsphysics.WindowFunction;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.exp;
import static java.lang.Math.sqrt;


/**
 * Basic filter method for a series of data
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class FilterModel{
	
	// for Killworth transformation only
	private static final DMatrixRMaj invA=new DMatrixRMaj(12,12);
	
	static{
		double[] daysOfMonths={31,28,31,30,31,30,31,31,30,31,30,31};
		double[] e=new double[12];
		double[] g=new double[12];
		double[] f=new double[12];
		
		e[0]=daysOfMonths[0]/(daysOfMonths[11]+daysOfMonths[0])/4f;
		g[0]=daysOfMonths[0]/(daysOfMonths[1 ]+daysOfMonths[0])/4f;
		f[0]=1f-e[0]-g[0];
		
		for(int i=1;i<11;i++){
			e[i]=daysOfMonths[i]/(daysOfMonths[i-1]+daysOfMonths[i])/4f;
			g[i]=daysOfMonths[i]/(daysOfMonths[i+1]+daysOfMonths[i])/4f;
			f[i]=1f-e[i]-g[i];
		}
		
		e[11]=daysOfMonths[11]/(daysOfMonths[10]+daysOfMonths[11])/4f;
		g[11]=daysOfMonths[11]/(daysOfMonths[0 ]+daysOfMonths[11])/4f;
		f[11]=1f-e[11]-g[11];
		
		DMatrixRMaj A=new DMatrixRMaj(12,12);
		
		for(int i=0;i<12;i++) A.set(i,i  ,f[i]);
		for(int i=0;i<11;i++) A.set(i,i+1,g[i]);
		for(int i=1;i<12;i++) A.set(i,i-1,e[i]);
		
		A.set(0 ,11,e[0 ]);
		A.set(11,0 ,g[11]);
		
		if(!CommonOps_DDRM.invert(A,invA)) throw new IllegalArgumentException("A cannot be inverted");
	}
	
	
	/**
     * get the linear trend (slope) of a given data series
     *
     * @param	data	an array of data
     *
     * @return	linear trend using 1st-order polynomial fit
     */
	public static float getLinearTrend(float[] data){
		PolynomialFitter pf=new PolynomialFitter(1,data.length);
		pf.fit(data);
		return pf.getCoefficient(1);
	}
	
	
	/**
     * Get the running mean of the given data
     *
     * @param	data	an array of data
     * @param	points	points of data used to mean
     *
     * @return	result	running mean of the data
     */
	public static float[] runningMean(float[] data,int points){
		float[] result=new float[data.length];
		
		runningMean(data,result,points);
		
		return result;
	}
	
	/**
     * Get the running mean of the given data
     *
     * @param	data	an array of data
     * @param	points	points of data used to mean
     * @param	undef	undefined value
     *
     * @return	result	running mean of the data
     */
	public static float[] runningMean(float[] data,int points,float undef){
		float[] re=new float[data.length];
		
		runningMean(data,re,points,undef);
		
		return re;
	}
	
	/**
     * Get the running mean of the given data
     *
     * @param	data	an array of data
     * @param	result	result of running of data
     * @param	points	points of data used to mean
     *
     * @return	result	running mean of the data
     */
	public static void runningMean(float[] data,float[] result,int points){
		int length=data.length;
		
		if(points%2!=1||points>length)
		throw new IllegalArgumentException("points should be a odd number and smaller than data length");
		
		if(length!=result.length) throw new IllegalArgumentException("array lengths not equal");
		
		int start=points/2,end=length-1-start;
		
		for(int i=0;i<start;i++) result[i]=data[i];
		
		for(int i=start;i<=end;i++){
			float sum=0;
			
			for(int j=i-start,J=i-start+points;j<J;j++) sum+=data[j];
			
			result[i]=sum/points;
		}
		
		for(int i=end+1;i<length;i++) result[i]=data[i];
	}
	
	/**
     * Get the running mean of the given data
     *
     * @param	data	an array of data
     * @param	result	result of running of data
     * @param	points	points of data used to mean
     * @param	undef	undefined value
     *
     * @return	result	running mean of the data
     */
	public static void runningMean(float[] data,float[] result,int points,float undef){
		int length=data.length;
		
		if(points%2!=1||points>length)
		throw new IllegalArgumentException("points should be a odd number and smaller than data length");
		
		if(length!=result.length) throw new IllegalArgumentException("array lengths not equal");
		
		int rad=points/2;
		
		for(int i=0;i<length;i++){
			float sum=0; int count=0;
			
			for(int j=i-rad,J=i+rad;j<=J;j++)
			if(data[i]==undef){
				break;
			}else if(j>=0&&j<length&&data[j]!=undef){
				sum+=data[j]; count++;
			}
			
			if(count==0) result[i]=undef;
			else result[i]=sum/count;
		}
	}
	
	
	/**
     * filter using specific window
     *
     * @param	data	an array of data
     * @param	win		window function (e.g., hann, hanning, hamming...)
     *
     * @return	result	running mean of the data
     */
	public static float[] windowFilter(float[] data,float[] win){
		float sum=StatisticsUtil.sum(win);
		
		float[] normWin=win.clone();
		
		for(int i=0,I=win.length;i<I;i++) normWin[i]=win[i]/sum;	// normalize window
		
		return MathsPhysicsUtil.convolveSame(data,normWin);
	}
	
	/**
     * hanning filter
     *
     * @param	data	an array of data
     * @param	L		number of points of the window i.e. full width at half maximum (FWHM)
     *
     * @return	result	running mean of the data
     */
	public static float[] hanningFilter(float[] data,int L){
		return windowFilter(data,WindowFunction.hanning(L));
	}
	
	public static void hanningFilter(float[] data,float[] re,int L){
		float[] normWin=WindowFunction.hanning(L);
		
		float sum=StatisticsUtil.sum(normWin);
		
		for(int i=0,I=normWin.length;i<I;i++) normWin[i]/=sum;	// normalize window
		
		System.arraycopy(MathsPhysicsUtil.convolve(data,normWin),L/2,re,0,data.length);
	}
	
	
	/**
     * hamming filter
     *
     * @param	data	an array of data
     * @param	L		number of points of the window i.e. full width at half maximum (FWHM)
     *
     * @return	result	running mean of the data
     */
	public static float[] hammingFilter(float[] data,int L){
		return windowFilter(data,WindowFunction.hamming(L));
	}
	
	public static void hammingFilter(float[] data,float[] re,int L){
		float[] normWin=WindowFunction.hamming(L);
		
		float sum=StatisticsUtil.sum(normWin);
		
		for(int i=0,I=normWin.length;i<I;i++) normWin[i]/=sum;	// normalize window
		
		System.arraycopy(MathsPhysicsUtil.convolve(data,normWin),L/2,re,0,data.length);
	}
	
	
	/**
     * band filter of one order Butterworth
     *
     * @param	tdata	a given temporal series
     * @param	t1		lower bound of period
     * @param	t2		upper bound of period
     *
     * @return	the result of the filtration
     */
	public static float[] ButterworthFilter(float[] tdata,float t1,float t2){
		float[] result=new float[tdata.length];
		
		ButterworthFilter(tdata,result,t1,t2);
		
		return result;
	}
	
	public static void ButterworthFilter(float[] tdata,float[] result,float t1,float t2){
		if(tdata.length!=result.length) throw new IllegalArgumentException("array lengths not equal");
		if(t1<=0||t1>t2) throw new IllegalArgumentException("illegal band arguments");
		
		int t=tdata.length;
		float a,b1,b2,w1,w2,dQ,Q2,tmp1,tmp2,tmp3;
		
		w1=(float)(2*PI/t1);	w2=(float)(2*PI/t2);
		
		dQ=2*(float)abs(sin(w1)/(1+cos(w1))-sin(w2)/(1+cos(w2)));
		Q2=(float)(4*sin(w1)*sin(w2)/((1+cos(w1))*(1+cos(w2))));
		
		a =2*dQ/(4+2*dQ+Q2);
		b1=2*(Q2-4)/(4+2*dQ+Q2);
		b2=(4-2*dQ+Q2)/(4+2*dQ+Q2);
		
		result[0]=result[1]=0;
		for(int l=2;l<t;l++) result[l]=a*(tdata[l]-tdata[l-2])-b1*result[l-1]-b2*result[l-2];
		
		tmp1=result[t-3];	result[t-3]=a*(result[t-3]-result[t-1])-b1*result[t-2]-b2*result[t-1];
		tmp2=result[t-4];	result[t-4]=a*(result[t-4]-result[t-2])-b1*result[t-3]-b2*result[t-2];
		for(int l=t-5;l>=0;l--){
			tmp3=result[l];
			result[l]=a*(result[l]-tmp1)-b1*result[l+1]-b2*result[l+2];
			tmp1=tmp2;	tmp2=tmp3;
		}
	}
	
	
	/**
     * Fast-Fourier band filter
     *
     * @param	data	temporal series
     * @param	Ks		wavenumbers of harmonics to be kept, start from 0 and end at data.length/2
     *
     * @return	temporal series after being filtered
     */
	public static float[] FFTFilter(float[] data,int... Ks){
		float[] result=new float[data.length];
		
		FFTFilter(data,result,Ks);
		
		return result;
	}
	
	public static void FFTFilter(float[] data,float[] r,int... Ks){
		int len=data.length;
		
		if(len!=r.length) throw new IllegalArgumentException("array lengths not equal");
		
		boolean[] keep=new boolean[len];
		
		for(int K:Ks){
			if(K<0||K>len/2)
			throw new IllegalArgumentException("wavenumber "+K+" should be in [0, "+len/2+"]");
			
			keep[K]=true; keep[len-K-1]=true;
		}
		
		final Complex zero=new Complex(0,0);
		
		Complex[] c=FastFourier.fft(data);
		for(int i=0;i<len;i++) if(!keep[i]) c[i]=zero;
		
		Complex[] d=FastFourier.ifft(c);
		for(int i=0;i<len;i++) r[i]=d[i].getReal();
	}
	
	
	/**
     * Fourier filter using HarmonicFitter
     *
     * @param	data	temporal series
     * @param	Ts		periods of the harmonics need to be removed, in unit of data point
     *
     * @return	temporal series after being filtered
     */
	public static float[] FourierFilter(float[] data,float... Ts){
		float[] re=new float[data.length];
		
		FourierFilter(data,re,Ts);
		
		return re;
	}
	
	public static void FourierFilter(float[] data,float[] r,float... Ts){
		int len=data.length;
		
		if(len!=r.length) throw new IllegalArgumentException("array lengths not equal");
		
		float[] re=new float[len];
		
		System.arraycopy(data,0,r ,0,len);
		System.arraycopy(data,0,re,0,len);
		
		HarmonicFitter hf=new HarmonicFitter(len,Ts);
		
		hf.fit(re);
		
		re=hf.cValues();
		
		for(int l=0;l<len;l++) r[l]-=re[l];
		
		System.arraycopy(r,0,re,0,len);
	}
	
	
	/**
     * M-term Guassian-Type low-pass filter
     *
     * @param	data	a given temporal series
     * @param	m		the term number used to running mean, it must be an odd number.
     *
     * @return	the result of the filtration
     */
	public static float[] GuassFilter(float[] data,int m){
		float[] result=new float[data.length];
		
		GuassFilter(data,result,m);
		
		return result;
	}
	
	public static void GuassFilter(float[] data,float[] result,int m){
		int   n  =data.length;
		int   nl =(m-1)/2;
		float c=2.15f;	// a tunable parameter, generally, c>2.0
		float cgm=nl/c;	// variance of Guassian distribution
		
		if(n!=result.length) throw new IllegalArgumentException("array lengths not equal");
		
		float[] xw=new float[n+m];
		float[] ck=new float[m];
		
		for(int i=0     ;i<=nl    ;i++) xw[i]=data[0   ];
		for(int i=nl+1  ;i<=n+nl-2;i++) xw[i]=data[i-nl];
		for(int i=n+nl-1;i<=n+m-1 ;i++) xw[i]=data[n-1 ];
		
		ck[nl]=(float)(1/(cgm*sqrt(2*PI)));
		for(int i=1;i<=nl;i++){
			ck[nl+i]=ck[nl]*(float)exp(-i*i/(2*cgm*cgm));
			ck[nl-i]=ck[nl+i];
		}
		
		for(int i=0;i<n;i++){
			result[i]=0;
			
			for(int j=0;j<=2*nl;j++) result[i]+=ck[j]*xw[i+j];
		}
	}
	
	
	/**
     * filter the cycle of a given length
     *
     * @param	tdata	a given temporal series
     * @param	cycle	length of the cycle
     *
     * @return	re		the mean cycle
     */
	public static float[] cycleFilter(float[] tdata,int cycle){
		float[] means=new float[cycle];
		
		cycleFilter(tdata,means,cycle);
		
		return means;
	}
	
	public static void cycleFilter(float[] tdata,float[] means,int cycle){
		int len=tdata.length;
		
		if(len%cycle!=0)
			throw new IllegalArgumentException("length of data cannot be divided by cycle");
		if(means.length!=cycle)
			throw new IllegalArgumentException("length of means should be equal to cycle");
		
		for(int ll=0,de=len/cycle;ll<cycle;ll++){
			means[ll]=0;
			
			for(int l=ll;l<len;l+=cycle) means[ll]+=tdata[l];
			
			means[ll]/=de;
			
			for(int l=ll;l<len;l+=cycle) tdata[l]-=means[ll];
		}
	}
	
	
	/**
     * remove the linear trend of a series
     *
     * @param	data	an array of data
     * @param	result	array that used to store the result
     *
     * @return	result	result after removing linear trend
     */
	public static void removeLinearTrend(float[] data){
		PolynomialFitter pf=new PolynomialFitter(1,data.length);
		
		pf.fit(data);
		
		float[] trend=pf.cValues();	// get trend
		
		for(int i=0,I=data.length;i<I;i++) data[i]-=trend[i];
	}
	
	public static void removeLinearTrend(float[] data,float[] result){
		if(data.length!=result.length) throw new IllegalArgumentException("data lengths not equal");
		
		PolynomialFitter pf=new PolynomialFitter(1,data.length);
		
		pf.fit(data);
		
		float[] trend=pf.cValues();
		
		for(int i=0,I=data.length;i<I;i++) result[i]=data[i]-trend[i];
	}
	
	
	/**
     * Hilbert transform, PI/2 filter, refer to Page 138 of statistics book
     *
     * @param	tdata	a given temporal series
     * @param	L		filter length
     *
     * @return	re		result of the filter
     */
	public static float[] HilbertTransform(float[] tdata,int L){
		float[] re=new float[tdata.length];
		float[] hl=new float[(L<<1)+1];
		
		HilbertTransform(tdata,re,hl,L);
		
		return re;
	}
	
	public static void HilbertTransform(float[] tdata,float[] re,float[] hl,int L){
		if(L<7||L>25) throw new IllegalArgumentException("invalid filter length");
		if(tdata.length!=re.length) throw new IllegalArgumentException("array lengths not equal");
		
		int length=(L<<1)+1;
		if(hl.length!=length) throw new IllegalArgumentException("not valid hl.length");
		
		for(int i=0  ;i<L     ;i++){ double tmp=sin(PI*(i-L)/2.0);	hl[i]=(float)(2.0*tmp*tmp/PI/(i-L));}
		for(int i=L+1;i<length;i++){ double tmp=sin(PI*(i-L)/2.0);	hl[i]=(float)(2.0*tmp*tmp/PI/(i-L));}
		
		for(int l=L,LL=re.length-L;l<LL;l++)
		for(int i=0;i<length;i++) re[l]+=tdata[l-(i-L)]*hl[i];
		
		float mean=StatisticsUtil.cArithmeticMean(tdata);
		for(int i=0,I=tdata.length;i<I;i++) re[i]+=mean;
	}
	
	
	/**
     * Killworth transformation of monthly climatology for linear interpolation
     * Reference: Killworth 1996, JPO
     *
     * @param	path	path of file after interpolation
     * @param	type	type of interpolation, e.g. cubic, linear
     * @param	T		T count after interpolation
     */
	public static float[] KillworthTransform(float[] monthly){
		float[] re=new float[12];
		
		KillworthTransform(monthly,re);
		
		return re;
	}
	
	public static void KillworthTransform(float[] monthly,float[] re){
		if(monthly.length!=12)
		throw new IllegalArgumentException("length of monthly should be 12 (actually "+monthly.length+")");
		
		if(re.length!=12)
		throw new IllegalArgumentException("length of re should be 12 (actually "+re.length+")");
		
		DMatrixRMaj Y=new DMatrixRMaj(12,1);
		DMatrixRMaj X=new DMatrixRMaj(12,1,true,
			monthly[0],monthly[1],monthly[2],monthly[3],monthly[4 ],monthly[5 ],
			monthly[6],monthly[7],monthly[8],monthly[9],monthly[10],monthly[11]
		);
		
		CommonOps_DDRM.mult(invA,X,Y);
		
		for(int i=0;i<12;i++) re[i]=(float)Y.get(i);
	}
	
	
	/** test
	public static void main(String[] args){
		float[] a=new float[]{3,5,2,6,8,3,5,8,6,4,6,5};
		
		float[] re=runningMean(a,3,8);
		
		for(int i=0;i<12;i++){
			System.out.println(a[i]+"\t"+re[i]);
		}
	}*/
}
