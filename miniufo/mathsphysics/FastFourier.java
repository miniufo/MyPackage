/**
 * @(#)Fourier.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;


/**
 * Fast Fourier class.
 * new for mixed-radix FFT
 * static for power-of-2 FFT
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class FastFourier{
	//
	private int N=0;				// length of N point FFT.
	private int NumofFactors=0;		// Number of factors of N.
	
	private int[] factors=null;		// Factors of N processed in the current stage.
	private int[] sofar  =null;		// Finished factors before the current stage.
	private int[] remain =null;		// Finished factors after the current stage.
	
	private float[] temRe   =null;	// Intermediate result of real part
	private float[] temIm   =null;	// Intermediate result of image part
	private float[] outputRe=null;	// real part of output
	private float[] outputIm=null;	// image part of output
	
	
	// Maximum numbers of factors allowed.
	private static final int MaxFactorsNumber=37;
	// Maximum factor of N.
	private static final int maxFactor=37;
	// 2*PI in float type
	private static final float TWO_PI=(float)(Math.PI*2.0);
	// cos2to3PI = cos(2*pi/3), using for 3 point FFT. cos(2*PI/3) is not -1.5
	private static final float cos2to3PI=-1.5000f;
	// sin2to3PI = sin(2*pi/3), using for 3 point FFT.
	private static final float sin2to3PI=8.6602540378444E-01f;
	// TwotoFivePI   = 2*pi/5.
	// c51, c52, c53, c54, c55 are used in fft5().
	// c51 =(cos(TwotoFivePI)+cos(2*TwotoFivePI))/2-1.
	private static final float c51=-1.25f;
	// c52 =(cos(TwotoFivePI)-cos(2*TwotoFivePI))/2.
	private static final float c52=5.5901699437495E-01f;
	// c53 = -sin(TwotoFivePI).
	private static final float c53=-9.5105651629515E-01f;
	// c54 =-(sin(TwotoFivePI)+sin(2*TwotoFivePI)).
	private static final float c54=-1.5388417685876E+00f;
	// c55 =(sin(TwotoFivePI)-sin(2*TwotoFivePI)).
	private static final float c55=3.6327126400268E-01f;
	// OnetoSqrt2 = 1/sqrt(2), used in fft8().
	private static final float OnetoSqrt2=7.0710678118655E-01f;
	
	
	/**
	 * constructor
     */
	public FastFourier(int N){
		this.N=N;
		outputRe=new float[N];
		outputIm=new float[N];
		
		factorize();
		
		// Allocate memory for intermediate result of FFT.
		temRe=new float[maxFactor];
		temIm=new float[maxFactor];
	}
	
	
	public void printFactors(){
		System.out.print("FFT Factors: [");
		for(int i=0,I=factors.length-1;i<I;i++) System.out.print(factors[i]+", ");
		System.out.print(factors[factors.length-1]);
		System.out.println("]");
	}
	
	
	/**
	 * mixed-radix (inverse) fast Fourier transform
	 */
	public void fftMixedRadix(float[] re,float[] im){
		if(re.length!=N||im.length!=N)
		throw new IllegalArgumentException("data lengths not equal");
		
		//clearBuffer(outputRe);	clearBuffer(outputIm);
		
		permute(re,im);
		
		for(int factorIndex=0;factorIndex<NumofFactors;factorIndex++)
		twiddle(factorIndex);
		
		//clearBuffer(temRe);		clearBuffer(temIm);
	}
	
	public void fftMixedRadix(float[] re){
		fftMixedRadix(re,new float[re.length]);
	}
	
	public void fftMixedRadix(Complex[] data){
		float[] re=new float[N];
		float[] im=new float[N];
		
		for(int i=0;i<N;i++){
			re[i]=data[i].getReal();
			im[i]=data[i].getImag();
		}
		
		fftMixedRadix(re,im);
	}
	
	public void ifftMixedRadix(float[] re,float[] im){
		float[] imNeg=new float[N];
		
		for(int i=0;i<N;i++) imNeg[i]=-im[i];
		
		fftMixedRadix(re,imNeg);
		
		for(int i=0;i<N;i++){
			outputRe[i]/= N;
			outputIm[i]/=-N;
		}
	}
	
	public void ifftMixedRadix(float[] re){
		fftMixedRadix(re);
		
		for(int i=0;i<N;i++){
			outputRe[i]/= N;
			outputIm[i]/=-N;
		}
	}
	
	public void ifftMixedRadix(Complex[] data){
		float[] re=new float[N];
		float[] im=new float[N];
		
		for(int i=0;i<N;i++){
			re[i]= data[i].getReal();
			im[i]=-data[i].getImag();
		}
		
		fftMixedRadix(re,im);
		
		for(int i=0;i<N;i++){
			outputRe[i]/= N;
			outputIm[i]/=-N;
		}
	}
	
	
	/**
	 * power-of-2 fast Fourier transform
	 */
	public void fftPowerOf2(float[] re,float[] im){
		if(re.length!=N||im.length!=N)
		throw new IllegalArgumentException("data lengths not equal");
		
		trf(re,im,-TWO_PI);
	}
	
	public void fftPowerOf2(float[] re){
		trf(re,new float[re.length],-TWO_PI);
	}
	
	public void fftPowerOf2(Complex[] data){
		float[] re=new float[N];
		float[] im=new float[N];
		
		for(int i=0;i<N;i++){
			re[i]=data[i].getReal();
			im[i]=data[i].getImag();
		}
		
		fftPowerOf2(re,im);
	}
	
	public void ifftPowerOf2(float[] re,float[] im){
		if(re.length!=N||im.length!=N)
		throw new IllegalArgumentException("data lengths not equal");
		
		trf(re,im,TWO_PI);
	}
	
	public void ifftPowerOf2(float[] re){
		trf(re,new float[re.length],TWO_PI);
	}
	
	public void ifftPowerOf2(Complex[] data){
		float[] re=new float[N];
		float[] im=new float[N];
		
		for(int i=0;i<N;i++){
			re[i]=data[i].getReal();
			im[i]=data[i].getImag();
		}
		
		ifftPowerOf2(re,im);
	}
	
	
	/*** getor and setor ***/
	public float[] getResultRealPart(){  return outputRe;}
	
	public float[] getResultImagePart(){ return outputIm;}
	
	public float[] getResultRealPartCopy(){  return outputRe.clone();}
	
	public float[] getResultImagePartCopy(){ return outputIm.clone();}
	
	public Complex[] getResult(){
		Complex[] re=new Complex[N];
		
		for(int i=0;i<N;i++)
		re[i]=new Complex(outputRe[i],outputIm[i]);
		
		return re;
	}
	
	
	/*** private method for mixed-radix fft ***/
	private void twiddle(int factorIndex){
		// Get factor data.
		int sofarRadix=sofar[factorIndex];
		int radix=factors[factorIndex];
		int remainRadix=remain[factorIndex];
		
		float tem;   // Temporary variable to do data exchange.
		
		float W=2*(float)Math.PI/(sofarRadix*radix);
		float cosW= (float)Math.cos(W);
		float sinW=-(float)Math.sin(W);
		
		float[] twiddleRe=new float[radix];
		float[] twiddleIm=new float[radix];
		
		float twRe=1.0f,twIm=0f;
		
		//Initialize twiddle addBk.address variables.
		int dataOffset=0,groupOffset=0,address=0;
		
		for(int dataNo=0;dataNo<sofarRadix;dataNo++){
			if(sofarRadix>1){
				twiddleRe[0]=1.0f;
				twiddleIm[0]=0.0f;
				twiddleRe[1]=twRe;
				twiddleIm[1]=twIm;
				
				for(int i=2;i<radix;i++){
					twiddleRe[i]=twRe*twiddleRe[i-1]-twIm*twiddleIm[i-1];
					twiddleIm[i]=twIm*twiddleRe[i-1]+twRe*twiddleIm[i-1];
				}
				
				tem =cosW*twRe-sinW*twIm;
				twIm=sinW*twRe+cosW*twIm;
				twRe=tem;
			}
			
			for(int groupNo=0;groupNo<remainRadix;groupNo++){
				if((sofarRadix>1)&&(dataNo>0)){
					temRe[0]=outputRe[address];
					temIm[0]=outputIm[address];
					
					int blockIndex=1;
					
					do{
						address=address+sofarRadix;
						temRe[blockIndex]=twiddleRe[blockIndex]*outputRe[address]-
						twiddleIm[blockIndex]*outputIm[address];
						temIm[blockIndex]=twiddleRe[blockIndex]*outputIm[address]+
						twiddleIm[blockIndex]*outputRe[address];
						
						blockIndex++;
					}while(blockIndex<radix);
					
				}else
					for(int i=0;i<radix;i++){
						temRe[i]=outputRe[address];
						temIm[i]=outputIm[address];
						address+=sofarRadix;
					}
				
				switch(radix){
				case 2:
					tem=temRe[0]+temRe[1];
					temRe[1]=temRe[0]-temRe[1];
					temRe[0]=tem;
					tem=temIm[0]+temIm[1];
					temIm[1]=temIm[0]-temIm[1];
					temIm[0]=tem;
					break;
				case 3:
					float t1Re=temRe[1]+temRe[2];
					float t1Im=temIm[1]+temIm[2];
					temRe[0]=temRe[0]+t1Re;
					temIm[0]=temIm[0]+t1Im;
					
					float m1Re=cos2to3PI*t1Re;
					float m1Im=cos2to3PI*t1Im;
					float m2Re=sin2to3PI*(temIm[1]-temIm[2]);
					float m2Im=sin2to3PI*(temRe[2]-temRe[1]);
					float s1Re=temRe[0]+m1Re;
					float s1Im=temIm[0]+m1Im;
					
					temRe[1]=s1Re+m2Re;
					temIm[1]=s1Im+m2Im;
					temRe[2]=s1Re-m2Re;
					temIm[2]=s1Im-m2Im;
					break;
				case 4:
					fft4(temRe, temIm);
					break;
				case 5:
					fft5(temRe, temIm);
					break;
				case 8:
					fft8();
					break;
				case 10:
					fft10();
					break;
				default:
					fftPrime(radix);
					break;
				}
				
				address=groupOffset;
				
				for(int i=0;i<radix;i++){
					outputRe[address]=temRe[i];
					outputIm[address]=temIm[i];
					address+=sofarRadix;
				}
				
				groupOffset+=sofarRadix*radix;
				address=groupOffset;
			}
			
			groupOffset=++dataOffset;
			address=groupOffset;
		}
	}
	
	private void fftPrime(int radix){
		// Initial WRe, WIm.
		float W=2*(float)Math.PI/radix;
		float cosW = (float)Math.cos(W);
		float sinW =-(float)Math.sin(W);
		float[] WRe=new float[radix];
		float[] WIm=new float[radix];
		
		WRe[0]=1;
		WIm[0]=0;
		WRe[1]=cosW;
		WIm[1]=sinW;
		
		for(int i=2;i<radix;i++){
			WRe[i]=cosW*WRe[i-1]-sinW*WIm[i-1];
			WIm[i]=sinW*WRe[i-1]+cosW*WIm[i-1];
		}
		
		// FFT of prime length data, using DFT, can be improved in the future.
		float rere, reim, imre, imim;
		int j, k;
		int max = (radix + 1) / 2;
		
		float[] tem1Re=new float[max];
		float[] tem1Im=new float[max];
		float[] tem2Re=new float[max];
		float[] tem2Im=new float[max];
		
		for(j=1;j<max;j++){
			tem1Re[j]=temRe[j]+temRe[radix - j];
			tem1Im[j]=temIm[j]-temIm[radix - j];
			tem2Re[j]=temRe[j]-temRe[radix - j];
			tem2Im[j]=temIm[j]+temIm[radix - j];
		}
		
		for(j=1;j<max;j++){
			temRe[j]=temRe[0];			temIm[j]=temIm[0];
			temRe[radix-j]=temRe[0];	temIm[radix-j]=temIm[0];
			
			k=j;
			
			for(int i=1;i<max;i++){
				rere=WRe[k]*tem1Re[i];
				imim=WIm[k]*tem1Im[i];
				reim=WRe[k]*tem2Im[i];
				imre=WIm[k]*tem2Re[i];
				
				temRe[radix-j]+=rere+imim;
				temIm[radix-j]+=reim-imre;
				temRe[j]+=rere-imim;
				temIm[j]+=reim+imre;
				
				k=k+j;
				
				if(k>=radix) k=k-radix;
			}
		}
		
		for(j=1;j<max;j++){
			temRe[0]=temRe[0]+tem1Re[j];
			temIm[0]=temIm[0]+tem2Im[j];
		}
	}
	
	private void permute(float[] re,float[] im){
		int[] count=new int[MaxFactorsNumber];
		int j;
		int k=0;
		
		for(int i=0;i<N-1;i++){
			outputRe[i]=re[k];
			outputIm[i]=im[k];
			
			j=0;
			k=k+remain[j];
			count[0]=count[0]+1;
			
			while(count[j]>=factors[j]){
				count[j]=0;
				k=k-(j==0?N:remain[j-1])+remain[j+1];
				j++;
				count[j]=count[j]+1;
			}
		}
		
		outputRe[N-1]=re[N-1];
		outputIm[N-1]=im[N-1];
	}
	
	private void factorize(){
		final int[] radices={2,3,4,5,8,10};
		
		int[] temFactors=new int[MaxFactorsNumber];
		
		// 1 - point FFT, no need to factorize N.
		if(N==1){
			temFactors[0]=1;
			NumofFactors=1;
		}
		
		// N - point FFT, N is needed to be factorized.
		int n=N;
		int index=0;    // index of temFactors.
		int i=radices.length - 1;
		
		while((n>1)&&(i>=0))
		if((n%radices[i])==0){
			n/=radices[i];
			temFactors[index++] = radices[i];
			
		}else i--;
		
		// Substitute 2x8 with 4x4.
		// index>0, in the case only one prime factor, such as N=263.
		if((index>0)&&(temFactors[index-1]==2))
		for(i=index-2;i>=0;i--)
		if(temFactors[i]==8){
			temFactors[index-1]=temFactors[i]=4;
			// break out of for loop, because only one '2' will exist in
			// temFactors, so only one substitutation is needed.
			break;
		}
		
		if(n>1){
			double K=Math.sqrt(n)+1;
			for(int k=2;k<K;k++)
			while((n%k)==0){
				n/=k;
				temFactors[index++] = k;
			}
			
			if(n>1) temFactors[index++]=n;
		}
		
		NumofFactors=index;
		//if(temFactors[NumofFactors-1] > 10)
		//   maxFactor = n;
		//else
		//   maxFactor = 10;
		
		// Inverse temFactors and store factors into factors[].
		factors=new int[NumofFactors];
		for(i=0;i<NumofFactors;i++)
		factors[i]=temFactors[NumofFactors-i-1];
		
		// Calculate sofar[], remain[].
		// sofar[]  : finished factors before the current stage.
		// factors[]: factors of N processed in the current stage.
		// remain[] : finished factors after the current stage.
		sofar =new int[NumofFactors];
		remain=new int[NumofFactors];
		
		remain[0]=N/factors[0];
		sofar[0]=1;
		for(i=1;i<NumofFactors;i++){
			sofar[i]=sofar[i-1]*factors[i-1];
			remain[i]=remain[i-1]/factors[i];
		}
	}
	
	private void fft4(float dataRe[],float dataIm[]){
		float t1Re,t1Im, t2Re,t2Im;
		float m2Re,m2Im, m3Re,m3Im;
		
		t1Re=dataRe[0]+dataRe[2];	t1Im=dataIm[0]+dataIm[2];
		t2Re=dataRe[1]+dataRe[3];	t2Im=dataIm[1]+dataIm[3];
		
		m2Re=dataRe[0]-dataRe[2];	m2Im=dataIm[0]-dataIm[2];
		m3Re=dataIm[1]-dataIm[3];	m3Im=dataRe[3]-dataRe[1];
		
		dataRe[0]=t1Re+t2Re;		dataIm[0]=t1Im+t2Im;
		dataRe[2]=t1Re-t2Re;		dataIm[2]=t1Im-t2Im;
		dataRe[1]=m2Re+m3Re;		dataIm[1]=m2Im+m3Im;
		dataRe[3]=m2Re-m3Re;		dataIm[3]=m2Im-m3Im;
	}
	
	private void fft5(float dataRe[],float dataIm[]){
		float t1Re,t1Im,t2Re,t2Im,t3Re,t3Im,t4Re,t4Im,t5Re,t5Im;
		float m1Re,m1Im,m2Re,m2Im,m3Re,m3Im,m4Re,m4Im,m5Re,m5Im;
		float s1Re,s1Im,s2Re,s2Im,s3Re,s3Im,s4Re,s4Im,s5Re,s5Im;
		
		t1Re=dataRe[1]+dataRe[4];	t1Im=dataIm[1]+dataIm[4];
		t2Re=dataRe[2]+dataRe[3];	t2Im=dataIm[2]+dataIm[3];
		t3Re=dataRe[1]-dataRe[4];	t3Im=dataIm[1]-dataIm[4];
		t4Re=dataRe[3]-dataRe[2];	t4Im=dataIm[3]-dataIm[2];
		
		t5Re=t1Re+t2Re;		t5Im=t1Im+t2Im;
		
		dataRe[0]=dataRe[0]+t5Re;
		dataIm[0]=dataIm[0]+t5Im;
		
		m1Re = c51 * t5Re;				m1Im = c51 * t5Im;
		m2Re = c52 * (t1Re - t2Re);		m2Im = c52 * (t1Im - t2Im);
		m3Re = -c53 * (t3Im + t4Im);	m3Im = c53 * (t3Re + t4Re);
		m4Re = -c54 * t4Im;				m4Im = c54 * t4Re;
		m5Re = -c55 * t3Im;				m5Im = c55 * t3Re;
		
		s3Re = m3Re - m4Re;			s3Im = m3Im - m4Im;
		s5Re = m3Re + m5Re;			s5Im = m3Im + m5Im;
		s1Re = dataRe[0] + m1Re;	s1Im = dataIm[0] + m1Im;
		s2Re = s1Re + m2Re;			s2Im = s1Im + m2Im;
		s4Re = s1Re - m2Re;			s4Im = s1Im - m2Im;
		
		dataRe[1] = s2Re + s3Re;	dataIm[1] = s2Im + s3Im;
		dataRe[2] = s4Re + s5Re;	dataIm[2] = s4Im + s5Im;
		dataRe[3] = s4Re - s5Re;	dataIm[3] = s4Im - s5Im;
		dataRe[4] = s2Re - s3Re;	dataIm[4] = s2Im - s3Im;
	}
	
	private void fft8(){
		float[] data1Re=new float[4];	float[] data1Im=new float[4];
		float[] data2Re=new float[4];	float[] data2Im=new float[4];
		
		// To improve the speed, use direct assaignment instead for loop here.
		data1Re[0]=temRe[0];	data2Re[0]=temRe[1];
		data1Re[1]=temRe[2];	data2Re[1]=temRe[3];
		data1Re[2]=temRe[4];	data2Re[2]=temRe[5];
		data1Re[3]=temRe[6];	data2Re[3]=temRe[7];
		
		data1Im[0]=temIm[0];	data2Im[0]=temIm[1];
		data1Im[1]=temIm[2];	data2Im[1]=temIm[3];
		data1Im[2]=temIm[4];	data2Im[2]=temIm[5];
		data1Im[3]=temIm[6];	data2Im[3]=temIm[7];
		
		fft4(data1Re, data1Im);	fft4(data2Re, data2Im);
		
		float tem =OnetoSqrt2*(data2Re[1]+data2Im[1]);
		data2Im[1]=OnetoSqrt2*(data2Im[1]-data2Re[1]);
		data2Re[1]=tem;
		tem=data2Im[2];
		data2Im[2]=-data2Re[2];
		data2Re[2]=tem;
		tem=OnetoSqrt2*(data2Im[3]-data2Re[3]);
		data2Im[3]=-OnetoSqrt2*(data2Re[3]+data2Im[3]);
		data2Re[3]=tem;
		
		temRe[0]=data1Re[0]+data2Re[0];	temIm[0]=data1Im[0]+data2Im[0];
		temRe[4]=data1Re[0]-data2Re[0];	temIm[4]=data1Im[0]-data2Im[0];
		temRe[1]=data1Re[1]+data2Re[1];	temIm[1]=data1Im[1]+data2Im[1];
		temRe[5]=data1Re[1]-data2Re[1];	temIm[5]=data1Im[1]-data2Im[1];
		temRe[2]=data1Re[2]+data2Re[2];	temIm[2]=data1Im[2]+data2Im[2];
		temRe[6]=data1Re[2]-data2Re[2];	temIm[6]=data1Im[2]-data2Im[2];
		temRe[3]=data1Re[3]+data2Re[3];	temIm[3]=data1Im[3]+data2Im[3];
		temRe[7]=data1Re[3]-data2Re[3];	temIm[7]=data1Im[3]-data2Im[3];
	}
	
	private void fft10(){
		float[] data1Re=new float[5];	float[] data1Im=new float[5];
		float[] data2Re=new float[5];	float[] data2Im=new float[5];
		
		// To improve the speed, use direct assaignment instead for loop here.
		data1Re[0]=temRe[0];	data1Im[0]=temIm[0];
		data2Re[0]=temRe[5];	data2Im[0]=temIm[5];
		data1Re[1]=temRe[2];	data1Im[1]=temIm[2];
		data2Re[1]=temRe[7];	data2Im[1]=temIm[7];
		data1Re[2]=temRe[4];	data1Im[2]=temIm[4];
		data2Re[2]=temRe[9];	data2Im[2]=temIm[9];
		data1Re[3]=temRe[6];	data1Im[3]=temIm[6];
		data2Re[3]=temRe[1];	data2Im[3]=temIm[1];
		data1Re[4]=temRe[8];	data1Im[4]=temIm[8];
		data2Re[4]=temRe[3];	data2Im[4]=temIm[3];
	    
		fft5(data1Re,data1Im);	fft5(data2Re,data2Im);

		temRe[0] = data1Re[0] + data2Re[0];	temIm[0] = data1Im[0] + data2Im[0];
		temRe[5] = data1Re[0] - data2Re[0];	temIm[5] = data1Im[0] - data2Im[0];
		temRe[6] = data1Re[1] + data2Re[1];	temIm[6] = data1Im[1] + data2Im[1];
		temRe[1] = data1Re[1] - data2Re[1];	temIm[1] = data1Im[1] - data2Im[1];
		temRe[2] = data1Re[2] + data2Re[2];	temIm[2] = data1Im[2] + data2Im[2];
		temRe[7] = data1Re[2] - data2Re[2];	temIm[7] = data1Im[2] - data2Im[2];
		temRe[8] = data1Re[3] + data2Re[3];	temIm[8] = data1Im[3] + data2Im[3];
		temRe[3] = data1Re[3] - data2Re[3];	temIm[3] = data1Im[3] - data2Im[3];
		temRe[4] = data1Re[4] + data2Re[4];	temIm[4] = data1Im[4] + data2Im[4];
		temRe[9] = data1Re[4] - data2Re[4];	temIm[9] = data1Im[4] - data2Im[4];
	}
	
	
	/**
	 * private method for power-of-2 fft
	 * results are usually a series of complex data.
	 * methods with 'component' suffix means that the results are returned in
	 * complex components, [0] is real part series and [1] is image part series
	 * of the complex result. 
	 *
	 * @param	data		temporal series
	 * @param	inverse		transform(TWO_PI) or inversetransform(-TWO_PI)
	 */
	private void trf(float[] re,float[] im,float inverse){
		int numBits=MathsPhysicsUtil.numberOf2(N);
		
		// Simultaneous data copy and bit-reversal ordering into output
		for(int i=0;i<N;i++){
			int j=reverseBits(i,numBits);
			outputRe[j]=re[i];
			outputIm[j]=im[i];
		}
		
		// FFT
		fftcc(outputRe,outputIm,inverse);
		
		if(inverse>0)
		for(int i=0;i<N;i++){ outputRe[i]/=N; outputIm[i]/=N;}
	}
	
	/**
	 * Common FFT code
	 *
	 * @param	arrayRe		real  parts of an array
	 * @param	arrayIm		image parts of an array
	 * @param	twoPi		TWO_PI for transform, -TWO_PI for inverse transform.
	 */
	private void fftcc(float[] arrayRe,float[] arrayIm,float twoPi){
		for(int blockSize=2,blockEnd=1;blockSize<=N;blockSize<<=1){
			float deltaAngle=twoPi/blockSize;
			float alpha=(float)Math.sin(0.5*deltaAngle);
			alpha*=2*alpha;
			float beta=(float)Math.sin(deltaAngle);
			
			for(int i=0;i<N;i+=blockSize){
				float angleRe=1;
				float angleIm=0;
				
				for(int j=i,n=0;n<blockEnd;j++,n++){
					int k=j+blockEnd;
					
					// tmp = angle*array[k]
					float tmpRe=angleRe*arrayRe[k]-angleIm*arrayIm[k];
					float tmpIm=angleRe*arrayIm[k]+angleIm*arrayRe[k];
					
					arrayRe[k]=arrayRe[j]-tmpRe;
					arrayIm[k]=arrayIm[j]-tmpIm;
					
					arrayRe[j]+=tmpRe;
					arrayIm[j]+=tmpIm;
					
					// angle = angle - (a-bi)*angle
					tmpRe=alpha*angleRe+beta*angleIm;
					tmpIm=alpha*angleIm-beta*angleRe;
					
					angleRe-=tmpRe;
					angleIm-=tmpIm;
				}
			}
			
			blockEnd=blockSize;
		}
	}
	
	/**
	 * Reverse bits, change 100 to 00(numBits==2), 001(numBits==3), 0010(numBits==4)
	 *
	 * @param	index		number to be reversed
	 * @param	numBits		length of the bits
	 *
	 * @return	rev			result of the reverse
	 */
	private int reverseBits(int index,int numBits){
		int i,rev;
		
		for(i=rev=0;i<numBits;i++){ rev=(rev<<1)|(index&1); index>>=1;}
		
		return rev;
	}
	
	
	/**
	 * FFT accessed in static way
	 */
	public static Complex[] fft(float[] re,float[] im){
		FastFourier ff=new FastFourier(re.length);
		
		if(MathsPhysicsUtil.isPowerOf2(re.length))
			ff.fftPowerOf2(re,im);
		else
			ff.fftMixedRadix(re,im);
		
		return ff.getResult();
	}
	
	public static Complex[] fft(float[] re){
		FastFourier ff=new FastFourier(re.length);
		
		if(MathsPhysicsUtil.isPowerOf2(re.length))
			ff.fftPowerOf2(re);
		else
			ff.fftMixedRadix(re);
		
		return ff.getResult();
	}
	
	public static Complex[] fft(Complex[] data){
		FastFourier ff=new FastFourier(data.length);
		
		if(MathsPhysicsUtil.isPowerOf2(data.length))
			ff.fftPowerOf2(data);
		else
			ff.fftMixedRadix(data);
		
		return ff.getResult();
	}
	
	public static Complex[] ifft(float[] re,float[] im){
		FastFourier ff=new FastFourier(re.length);
		
		if(MathsPhysicsUtil.isPowerOf2(re.length))
			ff.ifftPowerOf2(re,im);
		else
			ff.ifftMixedRadix(re,im);
		
		return ff.getResult();
	}
	
	public static Complex[] ifft(float[] re){
		FastFourier ff=new FastFourier(re.length);
		
		if(MathsPhysicsUtil.isPowerOf2(re.length))
			ff.ifftPowerOf2(re);
		else
			ff.ifftMixedRadix(re);
		
		return ff.getResult();
	}
	
	public static Complex[] ifft(Complex[] data){
		FastFourier ff=new FastFourier(data.length);
		
		if(MathsPhysicsUtil.isPowerOf2(data.length))
			ff.ifftPowerOf2(data);
		else
			ff.ifftMixedRadix(data);
		
		return ff.getResult();
	}
	
	
	/** test
	public static void main(String[] args){
		int len=105;
		float[] data=new float[len];
		
		for(int i=0;i<len;i++){
			data[i]=(float)(Math.sin(2*Math.PI*i/8)+Math.sin(2*Math.PI*i/20)+Math.sin(2*Math.PI*i/30))/3+
			(float)Math.random();
		}
		
		FastFourier ff=new FastFourier(len);
		
		ff.fftMixedRadix(data);
		float[] re1=ff.getResultRealPartCopy();
		float[] im1=ff.getResultImagePartCopy();
		
		ff.ifftMixedRadix(re1,im1);
		float[] data2=ff.getResultRealPart();
		
		Complex[] data3=FastFourier.ifft(FastFourier.fft(data));
		
		for(int i=0;i<len;i++)
		System.out.println(data[i]+"\t"+data2[i]+"\t"+data3[i].getReal());
		
		//for(int i=1;i<len/2;i++)
		//System.out.println(data[i]+"\t"+re2[i]+"\t"+re3[i]+"\t"+re1[i]+"\t"+im1[i]);
		//System.out.println(re1[i]+"\t"+re1[len-i]+"\t"+im1[i]+"\t"+-im1[len-i]);
	}*/
}
