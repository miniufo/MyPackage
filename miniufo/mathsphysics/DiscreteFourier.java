/**
 * @(#)Fourier.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import static java.lang.Math.PI;
import static java.lang.Math.atan;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;


/**
 * Fourier class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DiscreteFourier{
	//
	private int N   =0;			// length of the data
	private int half=0;			// half length of the data
	
	private float[][] sinw=null;	// sin(omega)
	private float[][] cosw=null;	// cos(omega)
	
	
	/**
	 * constructor
	 */
	public DiscreteFourier(int n){
		N=n;	half=n%2==0?n/2+1:(n+1)/2;
		
		if(N<=3) throw new IllegalArgumentException("data are too short for estimate");
		
		sinw=new float[half][N];
		cosw=new float[half][N];
		
		for(int k=0;k<half;k++){
			float o=(float)(2.0*PI*k/N);
			
			for(int i=0;i<N;i++){
				sinw[k][i]=(float)sin(o*i);
				cosw[k][i]=(float)cos(o*i);
			}
		}
	}
	
	
	/*** getor and setor ***/
	public float[][] getSinOmega(){ return sinw;}
	
	public float[][] getCosOmega(){ return cosw;}
	
	
	/**
     * calculate fourier coefficient
     *
     * @param	data	data to calculate
     * 
     * @return	re		re[0] is the Fourier coefficient of   sine component (a)
						re[1] is the Fourier coefficient of cosine component (b)
						re[2] is the amplitude spectrum = sqrt(a**2+b**2)/2 (c)
 						re[3] is the phase angle spectrum =atan(b/a) (cta)
						re[4] is the discrete power spectrum = (a**2+b**2)/2 (s)
     */
	public float[][] cFourierCoefficient(float[] data){
		if(data.length!=N) throw new IllegalArgumentException("data length is not "+N);
		
		float[][] re=new float[5][half];
		
		for(int k=0;k<half;k++){
			for(int i=0;i<N;i++){
				re[0][k]+=data[i]*cosw[k][i];
				re[1][k]+=data[i]*sinw[k][i];
			}
			
			re[0][k]*=2.0f/N;
			re[1][k]*=2.0f/N;

			re[2][k]=re[0][k]*re[0][k]+re[1][k]*re[1][k];
			re[4][k]=re[2][k];
			
			re[2][k]=(float)sqrt(re[2][k])/2;
			re[4][k]/=2;
			
			re[3][k]=(float)atan(re[1][k]/re[0][k]);
		}
		
		return re;
	}
	
	public double[][] cFourierCoefficient(double[] data){
		if(data.length!=N) throw new IllegalArgumentException("data length is not "+N);
		
		double[][] re=new double[5][half];
		
		for(int k=0;k<half;k++){
			for(int i=0;i<N;i++){
				re[0][k]+=data[i]*cosw[k][i];
				re[1][k]+=data[i]*sinw[k][i];
			}
			
			re[0][k]*=2.0/N;
			re[1][k]*=2.0/N;

			re[2][k]=re[0][k]*re[0][k]+re[1][k]*re[1][k];
			re[4][k]=re[2][k];
			
			re[2][k]=sqrt(re[2][k])/2.0;
			re[4][k]/=2.0;
			
			re[3][k]=atan(re[1][k]/re[0][k]);
		}
		
		return re;
	}
	
	
	/**
     * calculate fourier coefficient ak and bk
     *
     * @param	data	data to calculate
     * 
     * @return	re		re[0] is the Fourier coefficient of   sine component (a)
						re[1] is the Fourier coefficient of cosine component (b)
     */
	public float[][] cFourierAB(float[] data){
		if(data.length!=N) throw new IllegalArgumentException("data length is not "+N);
		
		float[][] re=new float[5][half];
		
		for(int k=0;k<half;k++)
		for(int i=0;i<N;i++){
			re[0][k]+=data[i]*cosw[k][i];
			re[1][k]+=data[i]*sinw[k][i];
		}
		
		return re;
	}
	
	public double[][] cFourierAB(double[] data){
		if(data.length!=N) throw new IllegalArgumentException("data length is not "+N);
		
		double[][] re=new double[5][half];
		
		for(int k=0;k<half;k++)
		for(int i=0;i<N;i++){
			re[0][k]+=data[i]*cosw[k][i];
			re[1][k]+=data[i]*sinw[k][i];
		}
		
		return re;
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			float PI=(float)Math.PI;
			float[] data=new float[512];
			
			int[] K=new int[]{1,2,3,4,5,6};
			
			for(int i=0;i<512;i++){
				for(int j=0;j<K.length;j++){
					data[i]+=(float)Math.sin(2*PI*1*K[j]*i/360);
				}
				
				data[i]+=100;	// ave
			}
			
			float[] re0=miniufo.statistics.FilterModel.FourierBandFilter(data,new int[]{0});
			float[] re1=miniufo.statistics.FilterModel.FourierBandFilter(data,new int[]{1});
			float[] re2=miniufo.statistics.FilterModel.FourierBandFilter(data,new int[]{2});
			float[] re3=miniufo.statistics.FilterModel.FourierBandFilter(data,new int[]{3});
			
			for(int i=0;i<512;i++) System.out.println(data[i]+"\t"+re0[i]+"\t"+re1[i]+"\t"+re2[i]+"\t"+re3[i]);
	    	
	    }catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
	}*/
}
