/**
 * @(#)EmpiricalModeDecomposition.java	1.0 2013.06.06
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import java.util.Arrays;
import java.util.Random;
import miniufo.basic.ArrayUtil;
import miniufo.statistics.StatisticsUtil;


/**
 * Empirical model decomposition
 *
 * @version 1.0, 2013.06.06
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EEMD{
	//
	private int len  =0;	// data length
	private int nmode=0;	// mode count
	private int itMax=10;	// maximum iteration
	private int enMax=100;	// ensemble number
	
	private float Nstd=0;	// std ratio
	
	private float[] dd=null;	// x-coordinates
	
	
	/**
	 * constructor
	 * 
	 * @param	len		length of data
	 * @param	nmode	number of modes to be decomposed
	 * @param	Nstd	ratio of the std. of the added noise and that of data
	 */
	public EEMD(int len,int nmode,float Nstd){
		if(len<3)   throw new IllegalArgumentException("length should be at least 3");
		if(nmode<1) throw new IllegalArgumentException("nmode should be at least 1");
		if(Nstd<0)  throw new IllegalArgumentException("Nstd should be non negative");
		
		this.len=len;
		this.nmode=nmode;
		this.Nstd=Nstd;
		
		dd=ArrayUtil.newMonotonousArray(len,1);
	}
	
	
	/*** getor and setor ***/
	public int getModes(){ return nmode;}
	
	public int getMaxIteration(){ return itMax;}
	
	public int getMaxEnsemble(){ return enMax;}
	
	public void setMaxIteration(int itMax){
		if(itMax<1) throw new IllegalArgumentException("itMax should be at least 1");
		this.itMax=itMax;
	}
	
	public void setMaxEnsemble(int enMax){
		if(enMax<1) throw new IllegalArgumentException("enMax should be at least 1");
		this.enMax=enMax;
	}
	
	
	/**
	 * decomposition of a given data
	 */
	public float[][] decomp(float[] data){
		float[][] allm=new float[nmode+2][len];
		
		decomp(data,allm);
		
		return allm;
	}
	
	public void decomp(float[] data,float[][] re){
		if(re.length!=nmode+2||re[0].length!=len)
		throw new IllegalArgumentException("invalid result sizes");
		
		if(data.length!=len)
		throw new IllegalArgumentException("invalid length of data");
		
		for(int j=0,J=nmode+2;j<J;j++)
		for(int i=0;i<len;i++) re[j][i]=0;
		
		// normalizing data
		float Ystd=StatisticsUtil.cStandardDeviation(data);
		for(int i=0;i<len;i++) data[i]/=Ystd;
		
		Random r=new Random();
		
		float[] R1=new float[len];
		float[] R2=new float[len];
		
		float[][] mode=new float[nmode+2][len];
		
		for(int m=0;m<enMax;m++){
			// noising data
			for(int i=0;i<len;i++){
				float rnd=(float)r.nextGaussian()*Nstd;
				R1[i]=data[i]+rnd;
				R2[i]=data[i]-rnd;
			}
			
			decompKernel(R1,dd,nmode,mode);
			
			for(int j=0,J=nmode+2;j<J;j++)
			for(int i=0;i<len;i++) re[j][i]+=mode[j][i];
			
			decompKernel(R2,dd,nmode,mode);
			
			for(int j=0,J=nmode+2;j<J;j++)
			for(int i=0;i<len;i++) re[j][i]+=mode[j][i];
		}
		
		for(int j=0,J=nmode+2;j<J;j++)
		for(int i=0;i<len;i++) re[j][i]*=Ystd/(enMax*2);
	}
	
	
	/*** helper methods ***/
	private static float[][] extremaMax(float[] data){
		int len =data.length;
		int cmax=1;
		
		// [0] are values and [1] are indices
		float[][] max=new float[2][len];
		
		max[0][0]=data[0]; max[1][0]=0;
		
		for(int i=1,I=len-1;i<I;i++)
		if(data[i-1]<=data[i]&&data[i]>=data[i+1]){
			max[0][cmax]=data[i];
			max[1][cmax]=i;
			cmax++;
		}
		
		max[0][cmax]=data[len-1]; max[1][cmax]=len-1; cmax++;
		
		if(cmax>=4){
			float slope=(max[0][1]-max[0][2])/(max[1][1]-max[1][2]);
			float tmp=slope*(max[1][0]-max[1][1])+max[0][1];
			if(tmp>max[0][0]) max[0][0]=tmp;
			
			slope=(max[0][cmax-2]-max[0][cmax-3])/(max[1][cmax-2]-max[1][cmax-3]);
			tmp=slope*(max[1][cmax-1]-max[1][cmax-2])+max[0][cmax-2];
			if(tmp>max[0][cmax-1]) max[0][cmax-1]=tmp;
		}
		
		float[] v=Arrays.copyOf(max[0],cmax);
		float[] i=Arrays.copyOf(max[1],cmax);
		
		return new float[][]{v,i};
	}
	
	private static float[][] extremaMin(float[] data){
		int len =data.length;
		int cmin=1;
		
		// [0] are values and [1] are indices
		float[][] min=new float[2][len];
		
		min[0][0]=data[0]; min[1][0]=0;
		
		for(int i=1,I=len-1;i<I;i++)
		if(data[i-1]>=data[i]&&data[i]<=data[i+1]){
			min[0][cmin]=data[i];
			min[1][cmin]=i;
			cmin++;
		}
		
		min[0][cmin]=data[len-1]; min[1][cmin]=len-1; cmin++;
		
		if(cmin>=4){
			float slope=(min[0][1]-min[0][2])/(min[1][1]-min[1][2]);
			float tmp=slope*(min[1][0]-min[1][1])+min[0][1];
			if(tmp<min[0][0]) min[0][0]=tmp;
			
			slope=(min[0][cmin-2]-min[0][cmin-3])/(min[1][cmin-2]-min[1][cmin-3]);
			tmp=slope*(min[1][cmin-1]-min[1][cmin-2])+min[0][cmin-2];
			if(tmp<min[0][cmin-1]) min[0][cmin-1]=tmp;
		}
		
		float[] v=Arrays.copyOf(min[0],cmin);
		float[] i=Arrays.copyOf(min[1],cmin);
		
		return new float[][]{v,i};
	}
	
	private void decompKernel(float[] data,float[] dd,int nmode,float[][] mode){
		int len=data.length;
		
		float[] resu=data.clone();
		
		System.arraycopy(data,0,mode[0],0,len);
		
		for(int i=0;i<nmode;i++){
			float[] tmp=resu.clone();
			
			for(int iter=0;iter<itMax;iter++){
				float[][] max=extremaMax(tmp);
				float[][] min=extremaMin(tmp);
				
				Spline s=null;
				
				s=new Spline(max[1],max[0]);
				s.cubicSplineWith1stBC(0,0);
				float[] upper=s.cValues(dd);
				
				s=new Spline(min[1],min[0]);
				s.cubicSplineWith1stBC(0,0);
				float[] lower=s.cValues(dd);
				
				for(int ii=0;ii<len;ii++) tmp[ii]-=(upper[ii]+lower[ii])/2f;
			}
			
			for(int ii=0;ii<len;ii++) resu[ii]-=tmp[ii];
			
			System.arraycopy(tmp,0,mode[i+1],0,len);
		}
		
		System.arraycopy(resu,0,mode[nmode+1],0,len);
	}
	
	
	/** test
	public static void main(String[] args){
		int nmode=4;
		int count=201;
		float[] data=TextReader.readColumn("d:/Matlab/EEMD/y.txt",count,false,1)[0];
		
		float[][] mode=new EEMD(count,nmode,0.02f).decomp(data);
		
		for(int i=0;i<count;i++){
			System.out.println(
				mode[0][i]+"\t"+
				mode[1][i]+"\t"+
				mode[2][i]+"\t"+
				mode[3][i]+"\t"+
				mode[4][i]+"\t"+
				mode[5][i]
			);
		}
	}*/
}
