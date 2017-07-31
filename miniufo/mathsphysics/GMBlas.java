/**
 * @(#)GMBlas.java	1.0 2014.08.08
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.basic.ArrayUtil;
import miniufo.statistics.StatisticsUtil;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;


/**
 * Gauss-Markov estimator implemented using BLAS
 *
 * @version 1.0, 2014.08.08
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GMBlas extends GaussMarkovEstimator{
	//
	private DoubleMatrix Y     =null;	// observations
	private DoubleMatrix A     =null;	// coefficients for the unknowns
	private DoubleMatrix AT    =null;	// transpose of matrix A
	private DoubleMatrix Rxx   =null;	// covariance matrix of the unknowns
	private DoubleMatrix Rnn   =null;	// variance structure of the eddy noise
	private DoubleMatrix Pxx   =null;	// error estimation
	private DoubleMatrix cycles=null;	// corresponding to the unknown matrix (Lumpkin and Johnson 2013, GRL)
										// fs=[1,2] corresponds to cycles=[x0 xc1 xc2 ... xs1 xs2 ...]
	
	
	/**
	 * constructor
	 * 
	 * @param	data	data of observations without any undefined values
	 */
	public GMBlas(float[] data,long[] ltims){
		len=data.length;
		
		if(len<2) throw new IllegalArgumentException("length should be larger than 1");
		if(len!=ltims.length) throw new IllegalArgumentException("lengths not equal");
		
		times=toDoubleTimes(ltims);
		hasST=false;
		
		Y=new DoubleMatrix(len,1);
		setData(data);
	}
	
	public GMBlas(float[] data,long[] ltims,float[] disx,float[] disy){
		len=data.length;
		
		if(len<2)
			throw new IllegalArgumentException("length should be larger than 1");
		
		if(len!=ltims.length||len!=disx.length||len!=disy.length)
			throw new IllegalArgumentException("lengths not equal");
		
		this.disx =disx;
		this.disy =disy;
		this.times=toDoubleTimes(ltims);
		this.hasST=true;
		
		Y=new DoubleMatrix(len,1);
		setData(data);
	}
	
	
	/**
	 * Gauss-Markov estimation of cycles
	 */
	public void estimateCycles(boolean computeError){
		if(fs==null) throw new IllegalArgumentException(
			"call setFrequenciesAndTL(double[] fs,double TL) first"
		);
		
		// RAT=Rxx*AT;
		// BLK=A*RAT+Rnn;
		// tmp=BLK^-1*Y;
		// cycles=RAT*tmp;
		DoubleMatrix RAT=Rxx.mmul(AT);
		DoubleMatrix BLK=A.mmul(RAT).addi(Rnn);
		DoubleMatrix tmp=Solve.solve(BLK,Y);
		cycles=RAT.mmul(tmp);
		
		if(computeError){
			// tmp=(A*RAT+Rnn)^-1*A
			// Pxx=Rxx-RAT*tmp*Rxx;
			tmp=Solve.solve(BLK,A);
			Pxx=RAT.mmul(tmp).mmul(Rxx).rsubi(Rxx);
		}
		
		dmn=cycles.get(0,0);
	}
	
	
	/**
	 * reconstruct the mean and cycles by A * cycles
	 * 
	 * @return	re	matrix form of cycles: [x0 xc1 xc2 ... xs1 xs2 ...]'
	 */
	public float[] reconstructMeanCycles(long[] ltims,float[] dx,float[] dy){
		checkResult();
		
		DoubleMatrix Aimg=assignAForm(toDoubleTimes(ltims),dx,dy);
		
		double[] tmp=Aimg.mmul(cycles).data;
		
		float[] re=new float[len];
		
		for(int i=0;i<len;i++) re[i]=(float)(tmp[i]+mean);
		
		return re;
	}
	
	public float[] reconstructMeanCycles(long[] ltims){
		if(hasST)
			return reconstructMeanCycles(ltims,new float[ltims.length],new float[ltims.length]);
		else
			return reconstructMeanCycles(ltims,null,null);
	}
	
	
	/**
	 * compute the variance contributions of each component
	 * 
	 * @return	re	for fs = {0, 1, 2},
	 * 				[0] the mean squares of GM mean,
	 * 				[1] the mean squares of annual cycle,
	 * 				[2] the mean squares of semiannual cycle,
	 * 				[3] the mean squares of spatial variations if any,
	 * 				[4] the mean squares of residuals,
	 */
	public float[] varianceContribution(){
		checkResult();
		
		int N=hasST? (2*fs.length+4) : (2*fs.length-1);	// No. of coefficients in cycles
		int C=hasST? fs.length : fs.length-1;			// No. of cycles excluding mean plus 1 for ST
		
		DoubleMatrix[] cycs=new DoubleMatrix[C];
		
		for(int i=0,I=cycs.length;i<I;i++) cycs[i]=new DoubleMatrix(N,1);
		
		// for cycles
		for(int i=1,I=fs.length-1;i<=I;i++){
			cycs[i-1].put(i  ,cycles.get(i  ));
			cycs[i-1].put(i+I,cycles.get(i+I));
		}
		
		// for spatial variations
		if(hasST){
			cycs[cycs.length-1].put(N-5,cycles.get(N-5));
			cycs[cycs.length-1].put(N-4,cycles.get(N-4));
			cycs[cycs.length-1].put(N-3,cycles.get(N-3));
			cycs[cycs.length-1].put(N-2,cycles.get(N-2));
			cycs[cycs.length-1].put(N-1,cycles.get(N-1));
		}
		
		double[][] tmp=new double[C][];
		
		for(int i=0,I=fs.length-1;i<I;i++) tmp[i]=A.mmul(cycs[i]).data;	// cycles
		if(hasST) tmp[tmp.length-1]=A.mmul(cycs[cycs.length-1]).data;	// spatial term
		
		float[] re=new float[C+2];
		
		re[0]=(float)((dmn+mean)*(dmn+mean));
		for(int i=0;i<C;i++) re[i+1]=(float)cMeanSquares(tmp[i]);
		re[re.length-1]=(float)cMeanSquares(A.mmul(cycles).rsubi(Y).data);
		
		return re;
	}
	
	
	/**
	 * reject residual of decomposition according to a threshold (Lumpkin and Johnson use 0.4)
	 */
	public boolean rejectResidual(float threshold){
		checkResult();
		
		// nrat=Y-A*cycles
		DoubleMatrix nrat=A.mmul(cycles).rsubi(Y);
		
		int count=0;
		for(int i=0,I=len;i<I;i++){
			nrat.put(i,0,Math.abs(nrat.get(i,0))/Math.sqrt(Rnn.get(i,i)));
			if(nrat.get(i,0)>1) count++;
		}
		
		if((float)count/len>threshold) return true;
		
		return false;
	}
	
	/**
	 * reject amplitude of decomposition according to a threshold (Lumpkin and Johnson use 1)
	 */
	public boolean rejectAmplitude(float threshold){
		checkResult();
		
		double[] xratio=new double[cycles.rows];
		
		for(int i=0,I=xratio.length;i<I;i++)
		xratio[i]=Math.abs(cycles.get(i,0))/Math.sqrt(Rxx.get(i,i));
		
		if(ArrayUtil.getMax(xratio)>threshold) return true;
		
		return false;
	}
	
	
	/*** getor and setor ***/
	
	/**
	 * get residual data after removing the mean and fs signals
	 * 
	 * @param	numOfCyc	number of cycles to be removed, started from 0, 
	 * 						2 corresponding to 0cycle (mean) + 1cycle + 2cycle (per year)
	 */
	public float[] getResidualData(int numOfCyc){
		checkResult();
		
		if(numOfCyc+1==fs.length){
			// res=Y-A*cycles
			DoubleMatrix res=A.mmul(cycles).rsubi(Y);
			
			float[] re=new float[res.rows];
			
			for(int i=0,I=re.length;i<I;i++) re[i]=(float)res.get(i);
			
			return re;
			
		}else{
			if(numOfCyc<0)
				throw new IllegalArgumentException("number of cycles to be removed should be at least 0 (mean)");
			if(numOfCyc>=fs.length)
				throw new IllegalArgumentException("only "+(fs.length-1)+" cycles can be removed");
			
			DoubleMatrix cycCopy=cycles.dup();
			
			int half=fs.length-1;
			
			for(int i=1,I=fs.length;i<I;i++)
			if(i>numOfCyc){
				cycCopy.put(i,0,0);
				cycCopy.put(i+half,0,0);
			}
			
			// res=Y-A*cycCopy
			DoubleMatrix res=A.mmul(cycCopy).rsubi(Y);
			
			float[] re=new float[res.rows];
			
			for(int i=0,I=re.length;i<I;i++) re[i]=(float)res.get(i);
			
			return re;
		}
	}
	
	/**
	 * compute estimated error for the unknowns (cycles)
	 * <=1 is OK while >1 is rejected
	 */
	public float[] getCycleErrors(){
		checkResult();
		
		float[] re=new float[Pxx.rows];
		
		for(int i=0,I=re.length;i<I;i++) re[i]=(float)Math.sqrt(Pxx.get(i,i));
		
		return re;
	}
	
	/**
	 * get coefficients of unknowns in float array
	 * 
	 * @return	re	matrix form of cycles: [x0 xc1 xc2 ... xs1 xs2 ... dx dx^2 dy dy^2 dxdy]'
	 */
	public float[] getCoefficients(boolean addDataMean){
		float[] re=new float[cycles.length];
		
		for(int i=0,I=cycles.length;i<I;i++) re[i]=(float)cycles.get(i);
		
		if(addDataMean) re[0]+=mean;
		
		return re;
	}
	
	/**
	 * set data without changing other parameters
	 */
	public void setData(float[] data){
		if(data.length!=len) throw new IllegalArgumentException("lengths not equal");
		
		mean=StatisticsUtil.cArithmeticMean(data);
		var =StatisticsUtil.cVariance(data);
		
		for(int j=0;j<len;j++) Y.put(j,0,data[j]-mean);
	}
	
	
	/*** helper methods ***/
	protected void assignRxx(){
		int N=hasST?2*fs.length+4:2*fs.length-1;
		
		Rxx=new DoubleMatrix(N,N);
		
		double rng=ArrayUtil.getRange(Y.data)/2;
		
		for(int i=0;i<N;i++) Rxx.put(i,i,rng*rng);
	}
	
	protected void assignRnn(AutoCorrType act,float... timescales){
		int N=times.length;
		
		double whitenoise_to_eddy=0.1;
		
		Rnn=new DoubleMatrix(N,N);
		
		for(int j=0;j<N;j++){
			for(int i=j;i<N;i++){
				float dT=(float)(times[i]-times[j]);
				Rnn.put(j,i,var*cAutoCorr(act,dT,timescales)*(1.0-whitenoise_to_eddy));
			}
			
			for(int i=0;i<j;i++) Rnn.put(j,i,Rnn.get(i,j));
			
			Rnn.put(j,j,Rnn.get(j,j)+whitenoise_to_eddy*var);
		}
	}
	
	protected void assignA(){
		if(hasST) A=assignAForm(times,disx,disy);
		else A=assignAForm(times,null,null);
		
		AT=A.transpose();
	}
	
	protected void checkResult(){
		if(cycles==null) throw new IllegalArgumentException("call estimateCycles(boolean) first");
	}
	
	protected DoubleMatrix assignAForm(double[] ntimes,float[] ndx,float[] ndy){
		int tlen=ntimes.length;
		
		if(hasST&&(tlen!=ndx.length||tlen!=ndy.length)) throw new IllegalArgumentException("lengths not equal");
		
		int N=hasST? (2*fs.length+4) : (2*fs.length-1);
		
		DoubleMatrix Aform=new DoubleMatrix(tlen,N);
		
		for(int j=0;j<tlen;j++){
			double tmp=2.0*Math.PI*ntimes[j];
			
			Aform.put(j,0,1);
			
			for(int ptr=1,I=fs.length-1;ptr<=I;ptr++){
				double angle=tmp*fs[ptr];
				
				Aform.put(j,ptr  ,Math.cos(angle));
				Aform.put(j,ptr+I,Math.sin(angle));
			}
			
			if(hasST){
				Aform.put(j,N-5,ndx[j]);
				Aform.put(j,N-4,ndx[j]*ndx[j]);
				Aform.put(j,N-3,ndy[j]);
				Aform.put(j,N-2,ndy[j]*ndy[j]);
				Aform.put(j,N-1,ndx[j]*ndy[j]);
			}
		}
		
		return Aform;
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
