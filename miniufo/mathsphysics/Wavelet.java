/**
 * @(#)Wavelet.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static miniufo.mathsphysics.MathsPhysicsUtil.gamma;
import static miniufo.mathsphysics.MathsPhysicsUtil.factorialDouble;


/**
 * wavelet class, base on the Morlet wavelet
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Wavelet{
	//
	private static final float  dt =1;		// data interval
	private static final float  s0 =2*dt;	// minimum scale, usually is 2*dt
	private static final float  dj =0.125f;	// scale resolution
	private static final float pad =0;		// padding number
	private static final float lag1=0.72f;	// lag-1 autocorrelation for red noise
	
	private int J;			// scale number
	private int N;			// number of points in the time series (after padding with 'pad')
	private int length;		// sample length, length <= N
	
	private float coi1;		// base for cone of influence
	private float para;		// the mother wavelet parameter
	private float ffctr;	// the ratio of Fourier period to scale
	
	private float dofmin;	// degrees of freedom with no smoothing
	//private float Cdelta;	// reconstruction factor
	private float gamFac;	// decorrelation factor for time  averaging
	//private float dj0;		// decorrelation factor for scale averaging
	
	private float[]  Sj =null;	// scales
	private float[] gws =null;	// global wavelet spectrum
	private float[] coi =null;	// cone of influence
	private float[] prd =null;	// period
	private float[] fftr=null;	// theoretical red-noise spectrum, Equation (16)
	
	private String mthr=null;	// String for mother wavelet: morlet, paul, dog
	
	private Complex[][]   wave  =null;	// wavelet coefficient
	private Complex[][] daughter=null;	// wavelet function
	
	
	/**
     * constructor
     *
     * @param	data	a time series of data for wavelet analysis
     */
	public Wavelet(float[] data,String mother){
		length=data.length;
		if(length<8) throw new IllegalArgumentException("not enough samples for analysis");
		
		/*** determine the length after padding and scale number ***/
		int ibase=(int)(log(length)/log(2)+0.4999);
		N=(int)pow(2,ibase+1);
		J=(int)(log(length*dt/s0)/log(2)/dj);	// scale number
		System.out.println("\nwavelet analysis...\npower-of-two is: "+ibase);
		
		/*** construct the wave number array, Equation (5) ***/
		float omega0=(float)(2*PI/(N*dt));
		float[] omegaK=new float[N];	omegaK[0]=0;
		
		for(int k=1;k<N/2+1;k++) omegaK[k]=k*omega0;
		for(int k=N/2+1;k<N;k++) omegaK[k]=-omegaK[N-k];
		
		/*** construct an array 'yfft' with padding for FFT, Equation (3) ***/
		Complex[] yfft=new Complex[N];
		for(int i=0;i<length;i++) yfft[i]=new Complex(data[i],0);
		for(int i=length;i<N;i++) yfft[i]=new Complex(pad    ,0);
		
		yfft=FastFourier.fft(yfft);
		
		/*** construct scale array ***/
		Sj=new float[J+1];
		for(int j=0;j<=J;j++) Sj[j]=s0*(float)pow(2,j*dj);
		
		/*** construct the daughter function ***/
		waveBases(mother,omegaK,Sj);
		
		/*** construct period array and wave array ***/
		prd =new float[J+1];
		fftr=new float[J+1];
		wave=new Complex[J+1][];
		
		/*** main wavelet loop ***/
		float variance=miniufo.statistics.StatisticsUtil.cVariance(data);
		for(int j=0;j<=J;j++){
			prd[j]=Sj[j]*ffctr;
			
			// theoretical red-noise spectrum as a function of period, Equation (16)
			fftr[j]=(float)(variance*(1-lag1*lag1)/(1-2*lag1*cos(2*PI*dt/prd[j])+lag1*lag1));
			
			for(int k=0;k<N;k++) daughter[j][k].multiplyEq(yfft[k]);
			
			wave[j]=FastFourier.ifft(daughter[j]);
		}
		
		// calculate global wavelet spectrum
		gws=new float[J+1];
		for(int j=0;j<=J;j++){
			for(int k=0;k<length;k++){ float tmp=wave[j][k].getMod(); gws[j]+=tmp*tmp;}
			gws[j]/=length;
		}
		
		
		// calculate cone of influence
		coi=new float[length];
		for(int i=0;i<length/2;i++){
			coi[i]=coi1*dt*i;
			coi[length-i-1]=coi[i];
		}
		if(length%2==1) coi[length/2]=coi1*dt*(length/2);
		
		/* reconstruction
		float Cdlt=0.776f;	// reconstruction factor
		float gama=2.32f;	// decorrelation factor for time averaging
		float dj0 =0.6f;	// factor for scale averaging
		float psi0=(float)pow(PI,-0.25f);
		rec=new float[N];
		for(int i=0;i<N;i++){
			for(int j=0;j<=J;j++) rec[i]+=wave[j][i].getReal()/Sj[j];
			
			rec[i]*=dj*(float)sqrt(dt)/Cdlt/psi0;
		}*/
	}
	
	
	// getor and setor
	public int getPadLength(){ return N;}
	
	public float[] getScales(){ return Sj;}
	
	public float[] getPeriods(){ return prd;}
	
	public float[] getConeOfInfluence(){ return coi;}
	
	public float[] getGlobalWaveletSpectrum(){ return gws ;}
	
	public Complex[][] getWaveCoefficient(){ return wave;}
	
	
	/**
	 * Do a regular chi-square test, i.e., Equation (18)
	 * 
	 * @param	level	confidence level, i.e., 0.95
	 * 
	 * @return	significance
	 */
	public float[] getChiSquareSignificance(float level){
		float[] sig=new float[J+1];
		
		for(int j=0;j<=J;j++)
		//sig[j]=fftr[j]*getChiSquare(level,dofmin)/dofmin;	// Equation (18)
		sig[j]=fftr[j]*(float)new ChiSquaredDistribution(dofmin).inverseCumulativeProbability(level)/dofmin;
		
		return sig;
	}
	
	/**
	 * Do a 'time-average' test, i.e., Equation (23)
	 * In this case, DOF should be set to NA
	 * The number of local wavelet spectra that were averaged together
	 * For the global wavelet spectrum, this would be NA=N,
	 * where N is the number of points in your time series
	 * 
	 * @param	level	confidence level, i.e., 0.95
	 * 
	 * @return	significance
	 */
	public float[] getGlobalSignificance(float level){
		float[] sig=new float[J+1];
		
		for(int j=0;j<=J;j++){
			float dof=N-Sj[j];
			
			if(dof<1) dof=1;
			
			// Equation (23)
			dof=dofmin*(float)sqrt(1+pow(dof*dt/gamFac/Sj[j],2));
			
			//sig[j]=fftr[j]*getChiSquare(level,dof)/dof;
			sig[j]=fftr[j]*(float)new ChiSquaredDistribution(dof).inverseCumulativeProbability(level)/dof;
		}
		
		return sig;
	}
	
	/**
	 * Do a 'scale-average' test, i.e., Equations (25) ~ (28)
	 * In this case, DOF should be set to s1 ~ s2,
	 * which gives the scale range that was averaged together
	 * 
	 * @param	level	confidence level, i.e., 0.95
	 * @param	s1		smaller scale
	 * @param	s2		larger scale
	 * 
	 * @return	significance
	public float getScaleAveragedSignificance(float level,float s1,float s2){
		if(s1<1||s1>s2) throw new IllegalArgumentException("s1 and s2 must satisfy: 1 < s1 < s2");
		
		int count=0;
		ArrayList<Integer> tag=new ArrayList<Integer>();
		
		for(int j=0;j<=J;j++)
		if(s1<=Sj[j]&&Sj[j]>=s2){ tag.add(j); count++;}
		
		if(count==0) throw new IllegalArgumentException("No valid scales between "+s1+" and "+s2);
		
		// Equation (25)
		float savg=0;
		for(int i=0;i<count;i++) savg+=1f/Sj[tag.get(i)];
		savg=1f/savg;
		
		// power-of-two midpoint
		float smid=(float)Math.exp((Math.log(s1)+Math.log(s2))/2);
		
		// Equation (28)
		float dof=(dofmin*count*savg/smid)*(float)sqrt(1+Math.pow(count*dj/dj0,2));
		
		float sig=-9999,totalFFT=0;
		for(int i=0;i<count;i++){
			totalFFT+=fftr[tag.get(i)]/Sj[tag.get(i)];	totalFFT*=savg;
			
			sig=(dj*dt/Cdelta/savg)*totalFFT*getChiSquare(level,dof)/dof;
		}
		
		return sig;
	}
	 */
	
	
	/**
     * get wavelet basis
     *
     * @param	mother	'Morlet', 'Paul' or 'Dog'
     * @param	omegaK	the Fourier frequencies at which to calculate the wavelet
     * @param	Sj		the wavelet scales
     */
	private void waveBases(String mother,float[] omegaK,float[] Sj){
		mthr=mother.toLowerCase();
		
		daughter=new Complex[J+1][N];
		
		if(mthr.equals("morlet")){
			para=6;
			ffctr=(float)(4*PI/(para+sqrt(2+para*para)));
			coi1=(float)(ffctr/sqrt(2));
			
			dofmin=2;	gamFac=2.32f;	//Cdelta=0.776f;	dj0=0.6f;	//psi0=(float)(pow(PI,0.25f));
			
			for(int j=0;j<=J;j++){
				float norm=(float)(sqrt(2*PI*Sj[j]/dt)*pow(PI,-0.25f));
				
				for(int k=0;k<N/2+1;k++){
					float expnt=-(float)(pow(Sj[j]*omegaK[k]-para,2)/2);
					daughter[j][k]=new Complex(norm*(float)exp(expnt),0);
				}
				
				for(int k=N/2+1;k<N;k++) daughter[j][k]=new Complex(0,0);
			}
			
		}else if(mthr.equals("paul")){
			para=4;
			ffctr=(float)(4*PI/(2*para+1));
			coi1=(float)(ffctr*sqrt(2));
			
			dofmin=2;	gamFac=1.17f;	//Cdelta=1.132f;	dj0=1.5f;	//psi0=1.079f;
			
			for(int j=0;j<=J;j++){
				float norm=(float)
				(sqrt(2*PI*Sj[j]/dt)*pow(2,para)/sqrt(para*factorialDouble((int)(2*para-1))));
				
				for(int k=0;k<N/2+1;k++){
					float expnt=-Sj[j]*omegaK[k];
					daughter[j][k]=new Complex((float)(norm*pow(-expnt,para)*exp(expnt)),0);
				}
				
				for(int k=N/2+1;k<N;k++) daughter[j][k]=new Complex(0,0);
			}
			
		}else if(mthr.equals("dog")){
			para=2;
			ffctr=(float)(2*PI*sqrt(2.0/(2*para+1)));
			coi1=(float)(ffctr/sqrt(2));
			
			// para=6 for [Cdelta=1.966f;	gamFac=1.37f;	dj0=0.97f;	psi0=0.884f;]
			dofmin=1;	gamFac=1.43f;	//Cdelta=3.541f;	dj0=1.4f;	//psi0=0.867f;	// para=2
			
			for(int j=0;j<=J;j++){
				Complex norm=new Complex(
					(float)(sqrt(2*PI*Sj[j]/dt)*sqrt(1/gamma(para+0.5f))),0
				);
				norm.multiplyEq(new Complex(0,1).pow(para));
				
				for(int k=0;k<N;k++){
					float sk=Sj[j]*omegaK[k];
					daughter[j][k]=norm.multiply((float)(pow(sk,para)*exp(-0.5f*sk*sk)));
				}
			}
			
		}else throw new IllegalArgumentException("mother should be 'Morlet', 'Paul' or 'Dog'");
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			System.out.println(new java.text.DecimalFormat(".000").format(123456.23456789f));
			System.out.println(String.format("%0.000s",123456.123456789f));
			
		}catch(Exception ex){ ex.printStackTrace();}
	}*/
}
