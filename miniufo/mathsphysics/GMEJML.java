/**
 * @(#)GMEJML.java	1.0 2014.08.25
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.basic.ArrayUtil;
import miniufo.statistics.StatisticsUtil;
import org.ejml.data.DMatrixRMaj;
import static org.ejml.dense.row.CommonOps_DDRM.addEquals;
import static org.ejml.dense.row.CommonOps_DDRM.mult;
import static org.ejml.dense.row.CommonOps_DDRM.subtract;
import static org.ejml.dense.row.CommonOps_DDRM.solve;
import static org.ejml.dense.row.CommonOps_DDRM.transpose;


/**
 * Gauss-Markov estimator implemented using EJML
 *
 * @version 1.0, 2014.08.25
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GMEJML extends GaussMarkovEstimator{
	//
	private DMatrixRMaj Y     =null;	// observations
	private DMatrixRMaj A     =null;	// coefficients for the unknowns
	private DMatrixRMaj AT    =null;	// transpose of matrix A
	private DMatrixRMaj Rxx   =null;	// covariance matrix of the unknowns
	private DMatrixRMaj Rnn   =null;	// variance structure of the eddy noise
	private DMatrixRMaj Pxx   =null;	// error estimation
	private DMatrixRMaj cycles=null;	// corresponding to the unknown matrix (Lumpkin and Johnson 2013, GRL)
										// fs=[1,2] corresponds to cycles=[x0 xc1 xc2 ... xs1 xs2 ...]
	
	
	/**
	 * constructor
	 * 
	 * @param	data	data of observations without any undefined values
	 * @param	time	times corresponding to the data (unit: year)
	 */
	public GMEJML(float[] data,long[] ltims){
		len=data.length;
		
		if(len<2) throw new IllegalArgumentException("length should be larger than 1");
		if(len!=ltims.length) throw new IllegalArgumentException("lengths not equal");
		
		times=toDoubleTimes(ltims);
		hasST=false;
		
		Y=new DMatrixRMaj(len,1);
		setData(data);
	}
	
	public GMEJML(float[] data,long[] ltims,float[] disx,float[] disy){
		len=data.length;
		
		if(len<2)
			throw new IllegalArgumentException("length should be larger than 1");
		
		if(len!=ltims.length||len!=disx.length||len!=disy.length)
			throw new IllegalArgumentException("lengths not equal");
		
		this.disx =disx;
		this.disy =disy;
		this.times=toDoubleTimes(ltims);
		this.hasST=true;
		
		Y=new DMatrixRMaj(len,1);
		setData(data);
	}
	
	
	/**
	 * Gauss-Markov estimation of cycles
	 */
	public void estimateCycles(boolean computeError){
		if(fs==null) throw new IllegalArgumentException("call setFrequenciesAndTL(double[] fs,double TL) first");
		
		int N=Rxx.getNumCols();
		
		DMatrixRMaj tmp=new DMatrixRMaj(len,1);
		DMatrixRMaj RAT=new DMatrixRMaj(N,len);
		DMatrixRMaj BLK=new DMatrixRMaj(len,len);
		
		// RAT=Rxx*AT;
		// BLK=A*RAT+Rnn;
		// tmp=BLK^-1*Y;
		// cycles=RAT*tmp;
		mult(Rxx,AT,RAT);
		mult(A,RAT,BLK);
		addEquals(BLK,Rnn);
		solve(BLK,Y,tmp);
		mult(RAT,tmp,cycles);
		
		if(computeError){
			DMatrixRMaj tmp2=new DMatrixRMaj(len,N);
			DMatrixRMaj tmp3=new DMatrixRMaj(N,N);
			DMatrixRMaj tmp4=new DMatrixRMaj(N,N);
			// tmp2=(A*RAT+Rnn)^-1*A
			// Pxx=Rxx-RAT*tmp2*Rxx;
			solve(BLK,A,tmp2);
			mult(RAT,tmp2,tmp3);
			mult(tmp3,Rxx,tmp4);
			subtract(Rxx,tmp4,Pxx);
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
		
		DMatrixRMaj Aimg=assignAForm(toDoubleTimes(ltims),dx,dy);
		DMatrixRMaj res=new DMatrixRMaj(len,1);
		
		mult(Aimg,cycles,res);
		
		double[] tmp=res.data;
		
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
		
		DMatrixRMaj[] cycs=new DMatrixRMaj[C];
		
		for(int i=0,I=cycs.length;i<I;i++) cycs[i]=new DMatrixRMaj(N,1);
		
		// for cycles
		for(int i=1,I=fs.length-1;i<=I;i++){
			cycs[i-1].set(i  ,cycles.get(i  ));
			cycs[i-1].set(i+I,cycles.get(i+I));
		}
		
		// for spatial variations
		if(hasST){
			cycs[cycs.length-1].set(N-5,cycles.get(N-5));
			cycs[cycs.length-1].set(N-4,cycles.get(N-4));
			cycs[cycs.length-1].set(N-3,cycles.get(N-3));
			cycs[cycs.length-1].set(N-2,cycles.get(N-2));
			cycs[cycs.length-1].set(N-1,cycles.get(N-1));
		}
		
		double[][] tmp=new double[C][];
		
		DMatrixRMaj res=new DMatrixRMaj(len,1);
		
		for(int i=0,I=fs.length-1;i<I;i++){
			mult(A,cycs[i],res);
			tmp[i]=res.data;	// cycles
		}
		if(hasST){
			mult(A,cycs[cycs.length-1],res);
			tmp[tmp.length-1]=res.data;	// spatial term
		}
		
		float[] re=new float[C+2];
		
		re[0]=(float)((dmn+mean)*(dmn+mean));
		for(int i=0;i<C;i++) re[i+1]=(float)cMeanSquares(tmp[i]);
		
		mult(A,cycles,res);
		subtract(Y,res,res);
		
		re[re.length-1]=(float)cMeanSquares(res.data);
		
		return re;
	}
	
	
	/**
	 * reject residual of decomposition according to a threshold (Lumpkin and Johnson use 0.4)
	 */
	public boolean rejectResidual(float threshold){
		checkResult();
		
		// nrat=Y-A*cycles
		DMatrixRMaj nrat=new DMatrixRMaj(len,1);
		mult(A,cycles,nrat);
		subtract(Y,nrat,nrat);
		
		int count=0;
		for(int i=0,I=len;i<I;i++){
			nrat.set(i,0,Math.abs(nrat.get(i,0))/Math.sqrt(Rnn.get(i,i)));
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
		
		double[] xratio=new double[cycles.numRows];
		
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
			DMatrixRMaj res=new DMatrixRMaj(len,1);
			mult(A,cycles,res);
			subtract(Y,res,res);
			
			float[] re=new float[res.numRows];
			
			for(int i=0,I=re.length;i<I;i++) re[i]=(float)res.get(i);
			
			return re;
			
		}else{
			if(numOfCyc<0)
				throw new IllegalArgumentException("number of cycles to be removed should be at least 0 (mean)");
			if(numOfCyc>=fs.length)
				throw new IllegalArgumentException("only "+(fs.length-1)+" cycles can be removed");
			
			DMatrixRMaj cycCopy=cycles.copy();
			
			int half=fs.length-1;
			
			for(int i=1,I=fs.length;i<I;i++)
			if(i>numOfCyc){
				cycCopy.set(i,0,0);
				cycCopy.set(i+half,0,0);
			}
			
			// res=Y-A*cycCopy
			DMatrixRMaj res=new DMatrixRMaj(len,1);
			mult(A,cycCopy,res);
			subtract(Y,res,res);
			
			float[] re=new float[res.numRows];
			
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
		
		float[] re=new float[Pxx.numRows];
		
		for(int i=0,I=re.length;i<I;i++) re[i]=(float)Math.sqrt(Pxx.get(i,i));
		
		return re;
	}
	
	/**
	 * get coefficients of unknowns in float array
	 * 
	 * @return	re	matrix form of cycles: [x0 xc1 xc2 ... xs1 xs2 ... dx dx^2 dy dy^2 dxdy]'
	 */
	public float[] getCoefficients(boolean addDataMean){
		float[] re=new float[cycles.numRows];
		
		for(int i=0,I=cycles.numRows;i<I;i++) re[i]=(float)cycles.get(i);
		
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
		
		for(int j=0;j<len;j++) Y.set(j,0,data[j]-mean);
	}
	
	
	/*** helper methods ***/
	protected void assignRxx(){
		int N=hasST?2*fs.length+4:2*fs.length-1;
		
		Rxx=new DMatrixRMaj(N,N);
		
		double rng=ArrayUtil.getRange(Y.data)/2;
		
		for(int i=0;i<N;i++) Rxx.set(i,i,rng*rng);
	}
	
	protected void assignRnn(AutoCorrType act,float... timescales){
		int N=times.length;
		
		double whitenoise_to_eddy=0.1;
		
		Rnn=new DMatrixRMaj(N,N);
		
		for(int j=0;j<N;j++){
			for(int i=j;i<N;i++){
				float dT=(float)(times[i]-times[j]);
				Rnn.set(j,i,var*cAutoCorr(act,dT,timescales)*(1.0-whitenoise_to_eddy));
			}
			
			for(int i=0;i<j;i++) Rnn.set(j,i,Rnn.get(i,j));
			
			Rnn.set(j,j,Rnn.get(j,j)+whitenoise_to_eddy*var);
		}
	}
	
	protected void assignA(){
		if(hasST) A=assignAForm(times,disx,disy);
		else A=assignAForm(times,null,null);
		
		AT=transpose(A,null);
	}
	
	protected void checkResult(){
		if(cycles==null) throw new IllegalArgumentException("call estimateCycles(boolean) first");
	}
	
	protected DMatrixRMaj assignAForm(double[] ntimes,float[] ndx,float[] ndy){
		int tlen=ntimes.length;
		
		if(hasST&&(tlen!=ndx.length||tlen!=ndy.length)) throw new IllegalArgumentException("lengths not equal");
		
		int N=hasST? (2*fs.length+4) : (2*fs.length-1);
		
		DMatrixRMaj Aform=new DMatrixRMaj(tlen,N);
		
		for(int j=0;j<tlen;j++){
			double tmp=2.0*Math.PI*ntimes[j];
			
			Aform.set(j,0,1);
			
			for(int ptr=1,I=fs.length-1;ptr<=I;ptr++){
				double angle=tmp*fs[ptr];
				
				Aform.set(j,ptr  ,Math.cos(angle));
				Aform.set(j,ptr+I,Math.sin(angle));
			}
			
			if(hasST){
				Aform.set(j,N-5,ndx[j]);
				Aform.set(j,N-4,ndx[j]*ndx[j]);
				Aform.set(j,N-3,ndy[j]);
				Aform.set(j,N-2,ndy[j]*ndy[j]);
				Aform.set(j,N-1,ndx[j]*ndy[j]);
			}
		}
		
		return Aform;
	}
	
	
	/** test
	public static void main(String[] args){
		float[][] data=TextReader.readColumn("d:/data.txt",808,false,1,2,3,4,5);
		float[] yy=data[0];
		float[] dx=data[1];
		float[] dy=data[2];
		float[] tt=data[3];
		float[] soi=data[4];
		
		long str=System.nanoTime();
		GaussMarkovEstimator2 gme=new GaussMarkovEstimator2(yy,tt,dx,dy);
		gme.setFrequenciesAndTL(2,0.0137f);
		gme.estimateCycles();
		System.out.println("using: "+(System.nanoTime()-str)/1000000+" times");
		
		System.out.println(gme.mean+"\t"+gme.dmn);
	}*/
}
