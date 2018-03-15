/**
 * @(#)EulerianStatistics.java	1.0 2013.04.15
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.basic.ArrayUtil;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.AttachedMeta;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.mathsphysics.Complex;
import miniufo.mathsphysics.GMBlas;
import miniufo.mathsphysics.GaussMarkovEstimator;
import miniufo.mathsphysics.GaussMarkovEstimator.AutoCorrType;
import miniufo.statistics.CurrentEllipse;
import miniufo.statistics.StatisticsUtil;
import miniufo.util.TicToc;


/**
 * Eulerian statistics using GM2 with spatial terms
 *
 * @version 1.0, 2013.04.15
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EulerianStatistics extends SingleParticleStatistics{
	//
	private int attachedLen=0;
	
	private BinStatistics bs=null;
	
	
	/**
	 * constructor
	 * 
	 * @param	ls			a list of particle
	 * @param	dd			grid descriptor
	 * @param	dlon		delta lon for periodic boundary condition
	 * @param	hasST		has spatial terms or not
	 */
	public EulerianStatistics(List<? extends Particle> ls,DataDescriptor dd,boolean hasST){
		super(ls,dd);
		
		this.attachedLen=ls.get(0).getRecord(0).getDataLength();
		this.bs=new BinStatistics(dd,ls,hasST);
	}
	
	
	/**
     * calculate the eddy kinetic energy
     *
     * @return	re	EKE
     */
	public Variable cEKE(){
		int y=dd.getYCount(),x=dd.getXCount();
		
		Variable re=new Variable("eke",true,new Range(1,1,y,x));
		re.setCommentAndUnit("eddy kinetic energy");
		re.setUndef(Record.undef);
		
		float[][] rdata=re.getData()[0][0];
		float[][] eke  =BinStatistics.cEKE(bs.ad[0],bs.ad[1]);
		
		for(int j=0;j<y;j++) System.arraycopy(eke[j],0,rdata[j],0,x);
		
		return re;
	}
	
	public Variable cCount(){ return new BinningStatistics(dd).binningCount(ls);}
	
	
	/**
     * calculate array bias
     * Reference: Poulain 2001, JMS
     * 
     * @param	diff	diffusivity tensor (m^2/s)
     * @param	dd		DataDescriptor template
     *
     * @return	ab	array bias
     */
	public Variable[] cConcentrationAndArrayBias(float[][] diff,DataDescriptor dd){
		if(diff.length!=2||diff[0].length!=2)
		throw new IllegalArgumentException("diffusivity should be a 2-rank tensor");
		
		DynamicMethodsInSC dm=new DynamicMethodsInSC(new SphericalSpatialModel(dd));
		
		Variable conc=dm.cConcentration(cCount());
		Variable[] ab=dm.c2DGradient(conc);
		
		ab[0].divideEq(conc);
		ab[1].divideEq(conc);
		
		Variable ubias=ab[0].multiply(-diff[0][0]).plusEq(ab[1].multiply(-diff[0][1]));
		Variable vbias=ab[0].multiply(-diff[1][0]).plusEq(ab[1].multiply(-diff[1][1]));
		
		ubias.setName("ubias");	ubias.setCommentAndUnit("zonal component of array bias");
		vbias.setName("vbias");	vbias.setCommentAndUnit("meridional component of array bias");
		
		return new Variable[]{conc,ubias,vbias};
	}
	
	public Variable[] cConcentrationAndArrayBias(float[][] diff){ return cConcentrationAndArrayBias(diff,dd);}
	
	/**
     * calculate the variance ellipses
     *
     * @return	re	[0] is major axis, [1] is minor axis and [2] is theta
     */
	public Variable[] cVarianceEllipse(){
		int y=dd.getYCount(),x=dd.getXCount();
		
		Variable[] re=new Variable[3];
		
		re[0]=new Variable("major",true,new Range(1,1,y,x));
		re[1]=new Variable("minor",true,new Range(1,1,y,x));
		re[2]=new Variable("theta",true,new Range(1,1,y,x));
		
		re[0].setCommentAndUnit("variance along the major axis");
		re[1].setCommentAndUnit("variance along the minor axis");
		re[2].setCommentAndUnit("the direction of the axis of principal variability");
		
		re[0].setUndef(Record.undef);
		re[1].setUndef(Record.undef);
		re[2].setUndef(Record.undef);
		
		float[][] r0data=re[0].getData()[0][0];
		float[][] r1data=re[1].getData()[0][0];
		float[][] r2data=re[2].getData()[0][0];
		
		float[][][] ve=BinStatistics.cVarianceEllipse(bs.ad[0],bs.ad[1]);
		
		for(int j=0;j<y;j++){
			System.arraycopy(ve[0][j],0,r0data[j],0,x);
			System.arraycopy(ve[1][j],0,r1data[j],0,x);
			System.arraycopy(ve[2][j],0,r2data[j],0,x);
		}
		
		return re;
	}
	
	/**
     * calculate the means within a bin
     *
     * @return	re	[0] is attd0, [1] is attd1...
     */
	public Variable[] cMeansOfBins(){
		int y=dd.getYCount(),x=dd.getXCount();
		
		Variable[] re=new Variable[attachedLen];
		
		for(int m=0;m<attachedLen;m++){
			re[m]=new Variable("attd"+m,true,new Range(1,1,y,x));
			re[m].setCommentAndUnit("gridded pseudo-Eulerian mean of AttachedData "+m);
			re[m].setUndef(Record.undef);
		}
		
		float[][][] bm=bs.cMeansOfBins();
		
		for(int m=0;m<attachedLen;m++)
		for(int j=0;j<y;j++) System.arraycopy(bm[m][j],0,re[m].getData()[0][0][j],0,x);
		
		return re;
	}
	
	/**
     * reconstruct the mean and cycles by GM using resolution of a given grid
     * 
	 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
     * @param	TL		a given Lagrangian timescale (unit: year)
     * @param	grid	a grid template
     *
     * @return	re		[0] is attd0, [1] is attd1...
     */
	public Variable[] reconstructMeanCyclesByGM(float[] freqs,float TL,DataDescriptor grid){
		int ny=grid.getYCount(),nx=grid.getXCount(),nt=grid.getTCount();
		
		Variable[] re=new Variable[attachedLen];
		
		for(int m=0;m<attachedLen;m++){
			re[m]=new Variable("attd"+m,false,new Range(nt,1,ny,nx));
			re[m].setCommentAndUnit("reconstructed mean of AttachedData "+m);
			re[m].setUndef(Record.undef);
		}
		
		float[][][][] coeffs=bs.reconstructMeanCycles(freqs,TL);
		
		for(int m=0;m<attachedLen;m++){
			float[][][] data=re[m].getData()[0];
			
			for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++) System.arraycopy(coeffs[m][j][i],0,data[m][j][i],0,nt);
		}
		
		return re;
	}
	
	/**
     * calculate the standard deviations within a bin
     *
     * @return	re	[0] is attd0, [1] is attd1...
     */
	public Variable[] cSTDsOfBins(){
		int y=dd.getYCount(),x=dd.getXCount();
		
		Variable[] re=new Variable[attachedLen];
		
		for(int m=0;m<attachedLen;m++){
			re[m]=new Variable("attd"+m+"std",true,new Range(1,1,y,x));
			re[m].setCommentAndUnit("gridded STD of AttachedData "+m);
			re[m].setUndef(Record.undef);
		}
		
		float[][][] bstd=bs.cSTDsOfBins();
		
		for(int m=0;m<attachedLen;m++)
		for(int j=0;j<y;j++) System.arraycopy(bstd[m][j],0,re[m].getData()[0][0][j],0,x);
		
		return re;
	}
	
	/**
     * calculate seasonal sampling bias (Lumpkin 2003, GRL)
     *
     * @return	re	[0] is the amplitude and [1] is the phase of sampling bias
     */
	public Variable[] cSeasonalSamplingBias(){
		int y=dd.getYCount(),x=dd.getXCount();
		
		Variable[] re=new Variable[2];
		re[0]=new Variable("ampli",true,new Range(1,1,y,x));
		re[1]=new Variable("phase",true,new Range(1,1,y,x));
		
		re[0].setCommentAndUnit("amplitude");
		re[1].setCommentAndUnit("phase");
		
		re[0].setValue(Record.undef);	re[0].setUndef(Record.undef);
		re[1].setValue(Record.undef);	re[1].setUndef(Record.undef);
		
		float[][] adata=re[0].getData()[0][0];
		float[][] pdata=re[1].getData()[0][0];
		float[][][] bias=bs.cSeasonalSamplingBias();
		
		for(int j=0;j<y;j++){
			System.arraycopy(bias[0][j],0,adata[j],0,x);
			System.arraycopy(bias[1][j],0,pdata[j],0,x);
		}
		
		return re;
	}
	
	/**
     * compute variance contributions of each components
     * 
	 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
     * @param	TL		a given Lagrangian timescale (unit: year)
     *
     * @return	re		[0] is mean squares of , [1] is v and [2] is t
     */
	public Variable[][] cVarianceContribution(float[] freqs,float TL){
		int y=dd.getYCount(),x=dd.getXCount(),numOfCyc=freqs.length;
		
		Variable[][] re=new Variable[attachedLen][numOfCyc*2+5];
		
		for(int m=0;m<attachedLen;m++){
			re[m][0]=new Variable("attd"+m+"ms",true,new Range(1,1,y,x));
			re[m][0].setCommentAndUnit("mean squares of AttachedData "+m);
			re[m][0].setUndef(Record.undef);
			re[m][0].setValue(Record.undef);
			
			re[m][1]=new Variable("attd"+m+"gvar",true,new Range(1,1,y,x));
			re[m][1].setCommentAndUnit("GM variance of AttachedData "+m);
			re[m][1].setUndef(Record.undef);
			re[m][1].setValue(Record.undef);
			
			re[m][2]=new Variable("attd"+m+"var",true,new Range(1,1,y,x));
			re[m][2].setCommentAndUnit("variance of AttachedData "+m);
			re[m][2].setUndef(Record.undef);
			re[m][2].setValue(Record.undef);
			
			re[m][3]=new Variable("attd"+m+"bm",true,new Range(1,1,y,x));
			re[m][3].setCommentAndUnit("variance of bin mean for AttachedData "+m);
			re[m][3].setUndef(Record.undef);
			re[m][3].setValue(Record.undef);
			
			re[m][4]=new Variable("attd"+m+"gm",true,new Range(1,1,y,x));
			re[m][4].setCommentAndUnit("variance of GM mean for AttachedData "+m);
			re[m][4].setUndef(Record.undef);
			re[m][4].setValue(Record.undef);
			
			re[m][numOfCyc*2+3]=new Variable("attd"+m+"st",true,new Range(1,1,y,x));
			re[m][numOfCyc*2+3].setCommentAndUnit("variance of spatial terms for AttachedData "+m);
			re[m][numOfCyc*2+3].setUndef(Record.undef);
			re[m][numOfCyc*2+3].setValue(Record.undef);
			
			re[m][numOfCyc*2+4]=new Variable("attd"+m+"res",true,new Range(1,1,y,x));
			re[m][numOfCyc*2+4].setCommentAndUnit("variance of residuals for AttachedData "+m);
			re[m][numOfCyc*2+4].setUndef(Record.undef);
			re[m][numOfCyc*2+4].setValue(Record.undef);
		}
		
		for(int m=0;m<attachedLen;m++)
		for(int i=1;i<=numOfCyc;i++){
			re[m][i+4]=new Variable("attd"+m+"c"+i,true,new Range(1,1,y,x));
			re[m][i+4].setCommentAndUnit("variance of freq. "+i+" for AttachedData "+m);
			re[m][i+4].setUndef(Record.undef);
			re[m][i+4].setValue(Record.undef);
		}
		
		float[][][][] vc=bs.cVarianceContribution(freqs,TL);
		
		for(int m=0;m<attachedLen;m++)
		for(int l=0,L=re[0].length;l<L;l++){
			float[][] rdata=re[m][l].getData()[0][0];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) rdata[j][i]=vc[m][j][i][l];
		}
		
		return re;
	}
	
	/**
     * calculate the seasonal means
     *
     * @return	re	means of different data of different seasons, [seasons][AttachedData]
     */
	public Variable[][] cSeasonalMeans(int[][] seasons,AttachedMeta... meta){
		return new BinningStatistics(dd).binningSeasonalData(ls,seasons,meta);
	}
	
	/**
     * calculate the amplitudes of cycles and the estimated errors for the mean
     * 
	 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
     * @param	TL		a given Lagrangian timescale (unit: year)
     *
     * @return	re	amplitudes in Matrix form:
     * 			|attd0m  attd0_fs0_amp  attd0_fs1_amp  ...  attd0_fs0_pha  attd0_fs1_pha  ...  attd0error|
     * 			|attd1m  attd1_fs0_amp  attd1_fs1_amp  ...  attd0_fs0_pha  attd0_fs1_pha  ...  attd1error|
     * 			|attd2m  attd2_fs0_amp  attd2_fs1_amp  ...  attd0_fs0_pha  attd0_fs1_pha  ...  attd2error|
     * 			|attd3m  attd3_fs0_amp  attd3_fs1_amp  ...  attd0_fs0_pha  attd0_fs1_pha  ...  attd3error| ...
     */
	public Variable[][] cCycleAmplitudesAndPhases(float[] freqs,float TL){
		int y=dd.getYCount(),x=dd.getXCount(),numOfCyc=freqs.length;
		
		Variable[][] re=new Variable[attachedLen][numOfCyc*2+4];
		
		for(int m=0;m<attachedLen;m++){
			re[m][0]=new Variable("attd"+m+"m",true,new Range(1,1,y,x));
			re[m][0].setCommentAndUnit("estimated mean of AttachedData "+m);
			re[m][0].setUndef(Record.undef);
			re[m][0].setValue(Record.undef);
			
			re[m][numOfCyc*2+1]=new Variable("attd"+m+"err",true,new Range(1,1,y,x));
			re[m][numOfCyc*2+1].setCommentAndUnit("estimated error of mean of AttachedData "+m);
			re[m][numOfCyc*2+1].setUndef(Record.undef);
			re[m][numOfCyc*2+1].setValue(Record.undef);
			
			re[m][numOfCyc*2+2]=new Variable("attd"+m+"Amask",true,new Range(1,1,y,x));
			re[m][numOfCyc*2+2].setCommentAndUnit("amplitude mask of GM solution for AttachedData "+m);
			re[m][numOfCyc*2+2].setUndef(Record.undef);
			re[m][numOfCyc*2+2].setValue(Record.undef);
			
			re[m][numOfCyc*2+3]=new Variable("attd"+m+"Rmask",true,new Range(1,1,y,x));
			re[m][numOfCyc*2+3].setCommentAndUnit("residual mask of GM solution for AttachedData "+m);
			re[m][numOfCyc*2+3].setUndef(Record.undef);
			re[m][numOfCyc*2+3].setValue(Record.undef);
		}
		
		for(int m=0;m<attachedLen;m++){
			for(int i=1;i<=numOfCyc;i++){
				re[m][i]=new Variable("attd"+m+"amp"+i,true,new Range(1,1,y,x));
				re[m][i].setCommentAndUnit("estimated cycle amplitude of AttachedData "+m+" for freq. "+i);
				re[m][i].setUndef(Record.undef);
				re[m][i].setValue(Record.undef);
			}
			
			for(int i=numOfCyc+1;i<=2*numOfCyc;i++){
				re[m][i]=new Variable("attd"+m+"amp"+i,true,new Range(1,1,y,x));
				re[m][i].setCommentAndUnit("estimated cycle sin phase of AttachedData "+m+" for freq. "+i);
				re[m][i].setUndef(Record.undef);
				re[m][i].setValue(Record.undef);
			}
		}
		
		float[][][][] amp=bs.cCycleAmplitudesAndPhases(freqs,TL);
		
		for(int m=0;m<attachedLen;m++){
			for(int l=0;l<=numOfCyc;l++){
				float[][] rdata=re[m][l].getData()[0][0];
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) rdata[j][i]=amp[m][j][i][l];
			}
			
			for(int l=numOfCyc+1;l<=2*numOfCyc;l++){
				float[][] rdata=re[m][l].getData()[0][0];
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) rdata[j][i]=amp[m][j][i][l];
			}
			
			float[][] rdata=re[m][numOfCyc*2+1].getData()[0][0];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			rdata[j][i]=amp[m][j][i][numOfCyc*2+1];
			
			rdata=re[m][numOfCyc*2+2].getData()[0][0];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			rdata[j][i]=amp[m][j][i][numOfCyc*2+2];
			
			rdata=re[m][numOfCyc*2+3].getData()[0][0];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			rdata[j][i]=amp[m][j][i][numOfCyc*2+3];
		}
		
		return re;
	}
	
	/**
     * calculate the elliptic parameters of cycles
     * 
	 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
     * @param	TL		a given Lagrangian timescale (unit: year)
     *
	 * @return	re		re[0][3] are the elliptic parameters for annual cycle,
	 * 					re[1][3] that for semiannual cycle, ...
     */
	public Variable[][] cCycleEllipses(float[] freqs,float TL){
		int y=dd.getYCount(),x=dd.getXCount(),numOfCyc=freqs.length;
		
		Variable[][] re=new Variable[numOfCyc][4];
		
		for(int m=0;m<numOfCyc;m++){
			re[m][0]=new Variable("major"+m,true,new Range(1,1,y,x));
			re[m][0].setCommentAndUnit("major axe for "+m+" cycle");
			re[m][0].setUndef(Record.undef);
			re[m][0].setValue(Record.undef);
			
			re[m][1]=new Variable("minor"+m,true,new Range(1,1,y,x));
			re[m][1].setCommentAndUnit("minor axe for "+m+" cycle");
			re[m][1].setUndef(Record.undef);
			re[m][1].setValue(Record.undef);
			
			re[m][2]=new Variable("theta"+m,true,new Range(1,1,y,x));
			re[m][2].setCommentAndUnit("inclination for "+m+" cycle");
			re[m][2].setUndef(Record.undef);
			re[m][2].setValue(Record.undef);
			
			re[m][3]=new Variable("phase"+m,true,new Range(1,1,y,x));
			re[m][3].setCommentAndUnit("phase for "+m+" cycle");
			re[m][3].setUndef(Record.undef);
			re[m][3].setValue(Record.undef);
		}
		
		float[][][][] params=bs.cCycleEllipses(freqs,TL);
		
		for(int m=0;m<numOfCyc;m++){
			float[][] rdata=re[m][0].getData()[0][0];
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) rdata[j][i]=params[m][j][i][0];
			
			rdata=re[m][1].getData()[0][0];
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) rdata[j][i]=params[m][j][i][1];
			
			rdata=re[m][2].getData()[0][0];
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) rdata[j][i]=params[m][j][i][2];
			
			rdata=re[m][3].getData()[0][0];
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) rdata[j][i]=params[m][j][i][3];
		}
		
		return re;
	}
	
	
	public void removeMeansOfBins(){ bs.removeMeansOfBins();}
	
	public void removeSeasonalBinMean(int[][] seasons){
		AttachedMeta[] meta=new AttachedMeta[attachedLen];
		
		for(int i=0;i<attachedLen;i++) meta[i]=new AttachedMeta("",i);
		
		Variable[][] sm=cSeasonalMeans(seasons,meta);
		
		float[][][][] data=new float[seasons.length][attachedLen][][];
		
		for(int i=0,I=seasons.length;i<I;i++)
		for(int j=0;j<attachedLen;j++) data[i][j]=sm[i][j].getData()[0][0];
		
		for(Particle p:ls)
		for(int l=0,L=p.getTCount();l<L;l++){
			Record rec=p.getRecord(l);
			
			int mon=new MDate(rec.getTime()).getMonth();
			
			int tagX=dd.getXNum(rec.getXPos());
			int tagY=dd.getYNum(rec.getYPos());
			
			for(int s=0,S=seasons.length;s<S;s++)
			for(int m=0,M=seasons[s].length;m<M;m++)
			if(mon==seasons[s][m]){
				for(int i=0;i<attachedLen;i++)
				rec.setData(meta[i],rec.getDataValue(meta[i])-data[s][i][tagY][tagX]);
			}
		}
	}
	
	public void removeCyclesByGM(float[] freqs,float TL,int noCycRem){ bs.removeCyclesByGM(freqs,TL,noCycRem);}
	
	public void removeCyclesByGM(float[] freqs,float TL){ bs.removeCyclesByGM(freqs,TL,freqs.length);}
	
	
	public void normalizeByBins(){ bs.normalizeByBin();}
	
	public void normalizeByDividingSTD(){bs.normalizeByDividingSTD();}
	
	public void normalizeByGM(float[] freqs,float TL){ bs.normalizeByGM(freqs,TL);}
	
	
	/**
	 * Project the current in along- and cross-stream components.
	 * 
	 * @param	umean	time-mean u flow
	 * @param	vmean	time-mean v flow
	 * @param	us		along-stream component
	 * @param	vn		cross-stream component
	 */
	public void projectCurrentAlongAndCrossStream(Variable umean,Variable vmean,AttachedMeta us,AttachedMeta vn){
		float undef=umean.getUndef();
		
		float[][] udata=umean.getData()[0][0];
		float[][] vdata=vmean.getData()[0][0];
		
		for(Particle p:ls)
		for(int l=0,L=p.getTCount();l<L;l++){
			Record r=p.getRecord(l);
			
			int itag=dd.getXNum(r.getXPos());
			int jtag=dd.getYNum(r.getYPos());
			
			float u=r.getDataValue(Particle.UVEL);
			float v=r.getDataValue(Particle.VVEL);
			
			float um=udata[jtag][itag];
			float vm=vdata[jtag][itag];
			
			if(um!=undef&&vm!=undef){
				float[] ac=CoordinateTransformation.projectToNaturalCoords(u,v,um,vm);
				
				r.setData(us,ac[0]);
				r.setData(vn,ac[1]);
				
			}else{
				r.setData(us,Record.undef);
				r.setData(vn,Record.undef);
			}
		}
	}
	
	
	public void maskoutByBinObservation(int threshold){
		// for list data
		float[][] cdata=bs.cCount();
		
		for(Particle p:ls)
		for(int l=0,L=p.getTCount();l<L;l++){
			Record rec=p.getRecord(l);
			
			int tagX=dd.getXNum(rec.getXPos());
			int tagY=dd.getYNum(rec.getYPos());
			
			if(cdata[tagY][tagX]<threshold){
				rec.setData(Particle.UVEL,Record.undef);
				rec.setData(Particle.VVEL,Record.undef);
			}
		}
		
		// for buffer data
		for(int m=0;m<attachedLen;m++)
		for(int j=0,J=dd.getYCount();j<J;j++)
		for(int i=0,I=dd.getXCount();i<I;i++)
		if(cdata[j][i]<threshold)
		for(int l=0,L=Math.round(cdata[j][i]);l<L;l++) bs.ad[m][j][i][l]=Record.undef;
	}
	
	
	/**
     * write ellipse data for matlab
     *
     * @param	path		path to write file
     * @param	threshold	write out while count larger than this value
     * @param	vs			variables to output, one in a column
     */
	public void writeDataForMatlab(String path,int threshold,Variable... vs){
		try(BufferedWriter br=new BufferedWriter(new FileWriter(path))){
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[] xdef=dd.getXDef().getSamples();
			float[] ydef=dd.getYDef().getSamples();
			
			float[][] cdata=bs.cCount();
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(cdata[j][i]>threshold){
				StringBuilder sb=new StringBuilder();
				
				sb.append(xdef[i]+"  "+ydef[j]+"  ");
				
				for(int m=0,M=vs.length-1;m<M;m++) sb.append(vs[m].getData()[0][0][j][i]+"  ");
				sb.append(vs[vs.length-1].getData()[0][0][j][i]+"\n");
				
				br.write(sb.toString());
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public void writeDataForMatlab(String path,float[][] mask,Variable... vs){
		try(BufferedWriter br=new BufferedWriter(new FileWriter(path))){
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[] xdef=dd.getXDef().getSamples();
			float[] ydef=dd.getYDef().getSamples();
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(mask[j][i]>0){
				StringBuilder sb=new StringBuilder();
				
				sb.append(xdef[i]+"  "+ydef[j]+"  ");
				
				for(int m=0,M=vs.length-1;m<M;m++) sb.append(vs[m].getData()[0][0][j][i]+"  ");
				sb.append(vs[vs.length-1].getData()[0][0][j][i]+"\n");
				
				br.write(sb.toString());
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/*** helper class ***/
	private static final class BinStatistics{
		//
		boolean hasST=false;
		
		private int attachedLen=0;
		
		private int[][] ptrs=null;
		
		// dimension is [y][x][t]
		private long[][][] tims=null;
		
		private float[][][] disx=null;
		private float[][][] disy=null;
		
		// dimension is [dataCount][y][x][t]
		private float[][][][] ad=null;	// attached data
		
		private Variable count=null;
		
		private DataDescriptor dd=null;
		
		private AttachedMeta[] meta=null;
		
		private List<? extends Particle> ls=null;
		
		
		/**
		 * constructor
		 */
		public BinStatistics(DataDescriptor dd,List<? extends Particle> ls,boolean hasST){
			this.hasST=hasST;
			this.dd   =dd;
			this.ls   =ls;
			
			count=new BinningStatistics(dd).binningCount(ls);
			
			int y=count.getYCount(),x=count.getXCount();
			
			float[][] data=count.getData()[0][0];
			
			attachedLen=ls.get(0).getRecord(0).getDataLength();
			meta=ls.get(0).getAttachedMeta();
			
			ptrs=new  int[y][x];
			tims=new long[y][x][];
			ad  =new float[attachedLen][y][x][];
			
			if(hasST){
				disx=new float[y][x][];
				disy=new float[y][x][];
			}
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				int len=Math.round(data[j][i]);
				
				tims[j][i]=new long[len];
				
				if(hasST){
					disx[j][i]=new float[len];
					disy[j][i]=new float[len];
				}
				
				for(int m=0;m<attachedLen;m++) ad[m][j][i]=new float[len];
			}
			
			griddingIntoBuffer();
		}
		
		
		/*** getor and setor ***/
		public float[][] cCount(){ return count.getData()[0][0];}
		
		
		/**
		 * @return	re	re[0] is mean of u, re[1] that of v,
		 * 				re[2] that of AttachedData 0, re[3] that of AttachedData 1, ...
		 */
		public float[][][] cMeansOfBins(){
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[][][] rdata=new float[attachedLen][y][x];
			
			for(int m=0;m<attachedLen;m++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			rdata[m][j][i]=StatisticsApplication.cArithmeticMean(ad[m][j][i],Record.undef);
			
			return rdata;
		}
		
		/**
		 * @return	re	re[0] is std. of u, re[1] that of v,
		 * 				re[2] that of AttachedData 0, re[3] that of AttachedData 1, ...
		 */
		public float[][][] cSTDsOfBins(){
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[][][] rdata=new float[attachedLen][y][x];
			
			for(int m=0;m<attachedLen;m++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			rdata[m][j][i]=StatisticsApplication.cStandardDeviation(ad[m][j][i],Record.undef);
			
			return rdata;
		}
		
		/**
		 * @return	re	re[0] is amplitude, re[1] is phase (unit: degree)
		 */
		public float[][][] cSeasonalSamplingBias(){
			int y=dd.getYCount(),x=dd.getXCount();
			
			final Complex undefCPLX=new Complex(Record.undef,0);
			
			Complex[][] bias=new Complex[y][x];
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) bias[j][i]=new Complex(Record.undef,0);
			
			for(Particle p:ls)
			for(int i=0,I=p.getTCount();i<I;i++){
				Record r=p.getRecord(i);
				
				int itag=dd.getXNumPeriodicX(r.getXPos());
				int jtag=dd.getYNum(r.getYPos());
				
				MDate md=new MDate(r.getTime());
				
				boolean leap=MDate.isLeapYear(md.getYear());
				
				Complex c=Complex.polar(1f,md.getDayOfYear()*(float)(2*Math.PI)/(leap?366f:365f));
				
				if(undefCPLX.equals(bias[jtag][itag])) bias[jtag][itag].set(c);
				else bias[jtag][itag].plusEq(c);
			}
			
			float[][] adata=new float[y][x];
			float[][] pdata=new float[y][x];
			float[][] cdata=count.getData()[0][0];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(cdata[j][i]>1) bias[j][i].divideEq(cdata[j][i]);
				
				if(undefCPLX.equals(bias[j][i])){
					adata[j][i]=Record.undef;
					pdata[j][i]=Record.undef;
					
				}else{
					adata[j][i]=bias[j][i].getMod();
					pdata[j][i]=bias[j][i].getArg()*(float)(180.0/Math.PI);
					
					if(pdata[j][i]<0) pdata[j][i]+=360;
				}
			}
			
			return new float[][][]{adata,pdata};
		}
		
		/**
		 * reconstruct mean and cycles by GM
		 * 
		 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
		 * @param	TL		Lagrangian integral timescale (year)
		 * 
		 * @return	re	re[0][y][x][t] is for zonal velocity,
		 * 				re[1][y][x][t] for meridional velocity,
		 * 				re[2][y][x][t] for temperature...
		 */
		public float[][][][] reconstructMeanCycles(float[] freqs,float TL){
			TicToc.tic("  Gauss-Markov decomposition 2 for coefficients of spatial terms");
			
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[][][][] rdata=new float[attachedLen][y][x][];
			
			List<Future<float[]>> ls=new ArrayList<>(x*attachedLen);
			ExecutorService es=ConcurrentUtil.defaultExecutor();
			CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
			
			final float per=y/20f;
			float curr=0;
			for(int j=0;j<y;j++){
				for(int i=0;i<x;i++)
				for(int m=0;m<attachedLen;m++)
				ls.add(cs.submit(new SolverCoefficients(
					ad[m][j][i],tims[j][i],
					hasST?disx[j][i]:null,
					hasST?disy[j][i]:null,
					freqs,TL)));
				
				try{
					for(int i=0,ptr=0;i<x;i++)
					for(int m=0;m<attachedLen;m++) rdata[m][j][i]=ls.get(ptr++).get();
				}
				catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
				
				if(j-curr>=per){
					curr=j;
					System.out.print(".");
				}
				
				ls.clear();
				
				System.gc();
			}
			
			TicToc.toc(TimeUnit.MINUTES);
			
			return rdata;
		}
		
		/**
		 * compute cycle amplitudes using Gauss-Markov decomposition
		 * 
		 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
		 * @param	TL		Lagrangian integral timescale
		 * 
		 * @return	re		re[0][y][x][fs.length+1] is the amplitude for AttachedData 0,
		 * 					re[1][y][x][fs.length+1] that for AttachedData 1, ...
		 */
		public float[][][][] cCycleAmplitudesAndPhases(float[] freqs,float TL){
			TicToc.tic("  Gauss-Markov decomposition 2 for cycle amplitudes");
			
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[][][][] rdata=new float[attachedLen][y][x][];
			
			List<Future<float[]>> ls=new ArrayList<>(x*attachedLen);
			ExecutorService es=ConcurrentUtil.defaultExecutor();
			CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
			
			final float per=y/20f;
			float curr=0;
			for(int j=0;j<y;j++){
				for(int i=0;i<x;i++)
				for(int m=0;m<attachedLen;m++)
				ls.add(cs.submit(new SolverAmplitudeAndPhase(
					ad[m][j][i],tims[j][i],
					hasST?disx[j][i]:null,
					hasST?disy[j][i]:null,
					freqs,TL
				)));
				
				try{
					for(int i=0,ptr=0;i<x;i++)
					for(int m=0;m<attachedLen;m++) rdata[m][j][i]=ls.get(ptr++).get();
				}
				catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
				
				if(j-curr>=per){
					curr=j;
					System.out.print(".");
				}
				
				ls.clear();
				
				System.gc();
			}
			
			TicToc.toc(TimeUnit.MINUTES);
			
			return rdata;
		}
		
		/**
		 * compute cycle amplitudes using Gauss-Markov decomposition
		 * 
		 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
		 * @param	TL		Lagrangian integral timescale
		 * 
		 * @return	re	re[0][y][x][3] are the elliptic parameters for annual cycle,
		 * 				re[1][y][x][3] that for semiannual cycle, ...
		 */
		public float[][][][] cCycleEllipses(float[] freqs,float TL){
			TicToc.tic("  Gauss-Markov decomposition 2 for coefficients of spatial terms");
			
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[][][][] rdata=new float[attachedLen][y][x][];
			float[][][][] re=new float[freqs.length][y][x][];
			
			List<Future<float[]>> ls=new ArrayList<>(x*attachedLen);
			ExecutorService es=ConcurrentUtil.defaultExecutor();
			CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
			
			final float per=y/20f;
			float curr=0;
			for(int j=0;j<y;j++){
				for(int i=0;i<x;i++)
				for(int m=0;m<attachedLen;m++)
				ls.add(cs.submit(new SolverCoefficients(
					ad[m][j][i],tims[j][i],
					hasST?disx[j][i]:null,
					hasST?disy[j][i]:null,
					freqs,TL)));
				
				try{
					for(int i=0,ptr=0;i<x;i++)
					for(int m=0;m<attachedLen;m++) rdata[m][j][i]=ls.get(ptr++).get();
				}
				catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
				
				if(j-curr>=per){
					curr=j;
					System.out.print(".");
				}
				
				ls.clear();
				
				System.gc();
			}
			
			TicToc.toc(TimeUnit.MINUTES);
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(rdata[0][j][i][0]==Record.undef){
					for(int f=0,F=freqs.length;f<F;f++)
					re[f][j][i]=new float[]{Record.undef,Record.undef,Record.undef,Record.undef};
					
				}else{
					for(int f=0,F=freqs.length;f<F;f++){
						float ucos=rdata[0][j][i][f+1];
						float usin=rdata[0][j][i][f+freqs.length];
						float vcos=rdata[1][j][i][f+1];
						float vsin=rdata[1][j][i][f+freqs.length];
						
						re[f][j][i]=CurrentEllipse.cCurrentEllipseByCoeffs(usin,ucos,vsin,vcos);
					}
				}
			}
			
			return re;
		}
		
		/**
		 * compute variance contributions of each components
		 * 
		 * @param	freqs	frequencies in unit of year^-1, 1 for annual and 2 for semiannual cycles
		 * @param	TL		Lagrangian integral timescale
		 * 
		 * @return	re		re[0][y][x][fs.length+5] is the amplitude for AttachedData 0,
		 * 					re[1][y][x][fs.length+5] that for AttachedData 1, ...
		 */
		public float[][][][] cVarianceContribution(float[] freqs,float TL){
			TicToc.tic("  Gauss-Markov decomposition 2 for variance contribution");
			
			int y=dd.getYCount(),x=dd.getXCount();
			
			float[][][][] rdata=new float[attachedLen][y][x][];
			
			List<Future<float[]>> ls=new ArrayList<>(x*attachedLen);
			ExecutorService es=ConcurrentUtil.defaultExecutor();
			CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
			
			final float per=y/20f;
			float curr=0;
			for(int j=0;j<y;j++){
				for(int i=0;i<x;i++)
				for(int m=0;m<attachedLen;m++)
				ls.add(cs.submit(new SolverVariance(
					ad[m][j][i],tims[j][i],
					hasST?disx[j][i]:null,
					hasST?disy[j][i]:null,
					freqs,TL
				)));
				
				try{
					for(int i=0,ptr=0;i<x;i++)
					for(int m=0;m<attachedLen;m++) rdata[m][j][i]=ls.get(ptr++).get();
				}
				catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
				
				if(j-curr>=per){
					curr=j;
					System.out.print(".");
				}
				
				ls.clear();
				
				System.gc();
			}
			
			TicToc.toc(TimeUnit.MINUTES);
			
			return rdata;
		}
		
		
		/**
		 * @return	re	re[0] is major, re[1] is minor and re[2] is theta
		 */
		public static float[][][] cVarianceEllipse(float[][][] ueddy,float[][][] veddy){
			int y=ueddy.length,x=ueddy[0].length;
			
			float[][] major=new float[y][x];
			float[][] minor=new float[y][x];
			float[][] theta=new float[y][x];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] params=CurrentEllipse.cVarianceEllipse(ueddy[j][i],veddy[j][i]);
				
				major[j][i]=params[0];
				minor[j][i]=params[1];
				theta[j][i]=params[2];
			}
			
			return new float[][][]{major,minor,theta};
		}
		
		/**
		 * @return	re	eddy kinetic energy
		 */
		public static float[][] cEKE(float[][][] ueddy,float[][][] veddy){
			int y=ueddy.length,x=ueddy[0].length;
			
			float[][] eke=new float[y][x];
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float sum=0;	int count=0;
				
				for(int l=0,L=ueddy[j][i].length;l<L;l++)
				if(ueddy[j][i][l]!=Record.undef&&veddy[j][i][l]!=Record.undef){
					sum+=(ueddy[j][i][l]*ueddy[j][i][l]+veddy[j][i][l]*veddy[j][i][l])/2f;
					count++;
				}
				
				if(count!=0) eke[j][i]=sum/count;
				else eke[j][i]=Record.undef;
			}
			
			return eke;
		}
		
		
		/**
		 * get anomalies by removing means of bins
		 */
		public void removeMeansOfBins(){
			for(int m=0;m<attachedLen;m++)
			for(int j=0,J=dd.getYCount();j<J;j++)
			for(int i=0,I=dd.getXCount();i<I;i++)
			StatisticsApplication.anomalize(ad[m][j][i],Record.undef);
			
			bufToList();
		}
		
		/**
		 * get anomalies by removing means using GM decomposition
		 */
		public void removeCyclesByGM(float[] freqs,float TL,int noCycRem){
			if(noCycRem>freqs.length) throw new IllegalArgumentException(
				"no. of Cycle ("+noCycRem+") should be smaller than length of freqs ("+freqs.length+")"
			);
				
			TicToc.tic("  Gauss-Markov decomposition 2 for residuals");
			
			int y=dd.getYCount(),x=dd.getXCount();
			
			ExecutorService es=ConcurrentUtil.defaultExecutor();
			CompletionService<Void> cs=new ExecutorCompletionService<>(es);
			
			final float per=y/20f;
			float curr=0;
			for(int j=0;j<y;j++){
				for(int i=0;i<x;i++)
				for(int m=0;m<attachedLen;m++)
				cs.submit(new SolverResidual(
					ad[m][j][i],tims[j][i],
					hasST?disx[j][i]:null,
					hasST?disy[j][i]:null,
					freqs,TL,noCycRem
				),null);
				
				try{
					for(int i=0;i<x;i++)
					for(int m=0;m<attachedLen;m++) cs.take();
				}
				catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
				
				if(j-curr>=per){
					curr=j;
					System.out.print(".");
				}
				
				System.gc();
			}
			
			TicToc.toc(TimeUnit.MINUTES);
			
			bufToList();
		}
		
		/**
		 * normalize using local mean and std. of bins
		 */
		public void normalizeByBin(){
			for(int m=0;m<attachedLen;m++)
			for(int j=0,J=dd.getYCount();j<J;j++)
			for(int i=0,I=dd.getXCount();i<I;i++)
			StatisticsApplication.standardize(ad[m][j][i],Record.undef);
			
			bufToList();
		}
		
		/**
		 * normalize by dividing std. of bins
		 */
		public void normalizeByDividingSTD(){
			for(int m=0;m<attachedLen;m++)
			for(int j=0,J=dd.getYCount();j<J;j++)
			for(int i=0,I=dd.getXCount();i<I;i++) dividingSTD(ad[m][j][i],Record.undef);
			
			bufToList();
		}
		
		/**
		 * normalize using anomalies that derived by GM decomposition
		 */
		public void normalizeByGM(float[] freqs,float TL){
			removeCyclesByGM(freqs,TL,freqs.length);
			
			for(int m=0;m<attachedLen;m++)
			for(int j=0,J=dd.getYCount();j<J;j++)
			for(int i=0,I=dd.getXCount();i<I;i++)
			StatisticsApplication.standardize(ad[m][j][i],Record.undef);
			
			bufToList();
		}
		
		
		/*** helper methods ***/
		private DataPair getValidData(long[] tims,float[] data,float[] disx,float[] disy){
			int NC=data.length;
			int UC=0;
			
			for(int l=0;l<NC;l++) if(data[l]==Record.undef) UC++;
			
			if(UC==0) return new DataPair(tims,data,disx,disy);
			
			 long[] nt=new  long[NC-UC];
			float[] nd=new float[NC-UC];
			float[] nx=hasST?new float[NC-UC]:null;
			float[] ny=hasST?new float[NC-UC]:null;
			
			for(int l=0,ptr=0;l<NC;l++) if(data[l]!=Record.undef){
				nd[ptr]=data[l];
				nt[ptr]=tims[l];
				
				if(hasST){
					nx[ptr]=disx[l];
					ny[ptr]=disy[l];
				}
				
				ptr++;
			}
			
			return new DataPair(nt,nd,nx,ny);
		}
		
		private void dividingSTD(float[] tdata,float undef){
			float[] ori=tdata.clone();
			
			int length=tdata.length,count=0;
			
			double var=0;
			
			for(int l=0;l<length;l++) if(tdata[l]!=undef) count++;
			
			if(count>=2){
				for(int l=0;l<length;l++)
				if(tdata[l]!=undef) var+=tdata[l]*tdata[l];
				
				var=Math.sqrt(var/(count-1));
				
				for(int l=0;l<length;l++) if(tdata[l]!=undef) tdata[l]/=var;
				
				if(tdata[0]==0&&tdata[1]==0)
				for(int l=0;l<length;l++) System.out.println(Arrays.toString(ori));
				
			}else{
				for(int l=0;l<length;l++) tdata[l]=undef;
			}
		}
		
		private void bufToList(){
			int y=dd.getYCount(),x=dd.getXCount();
			
			int[][] tmpPtrs=new int[y][x];
			
			for(Particle p:ls)
			for(int l=0,L=p.getTCount();l<L;l++){
				Record r=p.getRecord(l);
				
				int itag=dd.getXNumPeriodicX(r.getXPos());
				int jtag=dd.getYNum(r.getYPos());
				
				for(int m=0;m<attachedLen;m++)
				r.setData(meta[m],ad[m][jtag][itag][tmpPtrs[jtag][itag]]);
				
				tmpPtrs[jtag][itag]++;
			}
		}
		
		private void griddingIntoBuffer(){
			float[] xdef=dd.getXDef().getSamples();
			float[] ydef=dd.getYDef().getSamples();
			
			for(Particle p:ls)
			for(int l=0,L=p.getTCount();l<L;l++){
				Record r=p.getRecord(l);
				
				int itag=dd.getXNumPeriodicX(r.getXPos());
				int jtag=dd.getYNum(r.getYPos());
				int ltag=ptrs[jtag][itag];
				
				for(int m=0;m<attachedLen;m++)
				ad[m][jtag][itag][ltag]=r.getDataValue(meta[m]);
				
				tims[jtag][itag][ltag]=r.getTime();
				
				if(hasST){
					disx[jtag][itag][ltag]=r.getXPos()-xdef[itag];
					disy[jtag][itag][ltag]=r.getYPos()-ydef[jtag];
				}
				
				ptrs[jtag][itag]++;
			}
		}
		
		
		/*** helper classes for multi-thread execution since GM decomposition is slow ***/
		private final class SolverAmplitudeAndPhase implements Callable<float[]>{
			//
			float TL=0;
			
			 long[] tims=null;
			float[] data=null;
			float[] disx=null;
			float[] disy=null;
			float[] freq=null;
			
			SolverAmplitudeAndPhase(float[] data,long[] tims,float[] disx,float[] disy,float[] freqs,float TL){
				this.data=data;
				this.tims=tims;
				this.disx=disx;
				this.disy=disy;
				this.freq=freqs;
				this.TL=TL;
			}
			
			public float[] call(){
				DataPair dp=getValidData(tims,data,disx,disy);
				
				 long[] nt=dp.tims;
				float[] nd=dp.data;
				float[] nx=dp.disx;
				float[] ny=dp.disy;
				
				int rlen=freq.length*2+4;
				
				if(nd.length<=1||ArrayUtil.getRange(nd)==0f||StatisticsUtil.cVariance(nd)==0f){
					// case where length of array is too short or
					// case where the matrix is singular
					float[] re=new float[rlen];
					
					for(int l=0,L=re.length;l<L;l++) re[l]=Record.undef;
					
					return re;
					
				}else{
					GaussMarkovEstimator gme=hasST?new GMBlas(nd,nt,nx,ny):new GMBlas(nd,nt);
					gme.setFrequenciesAndTimescales(AutoCorrType.TCosExp,freq,TL);
					gme.estimateCycles(true);
					
					float[] amp=gme.getCycleAmplitudes();
					float[] pha=gme.getCycleSinPhase();
					float[] err=gme.getCycleErrors();
					
					float[] re=new float[rlen];
					
					System.arraycopy(amp,0,re,0,amp.length);
					System.arraycopy(pha,0,re,amp.length,pha.length);
					re[re.length-3]=err[0];
					re[re.length-2]=gme.rejectAmplitude(1)?-1:1;
					re[re.length-1]=gme.rejectResidual(0.4f)?-1:1;
					
					return re;
				}
			}
		}
		
		private final class SolverVariance implements Callable<float[]>{
			//
			float TL=0;
			
			 long[] tims=null;
			float[] data=null;
			float[] disx=null;
			float[] disy=null;
			float[] freq=null;
			
			SolverVariance(float[] data,long[] tims,float[] disx,float[] disy,float[] freqs,float TL){
				this.data=data;
				this.tims=tims;
				this.disx=disx;
				this.disy=disy;
				this.freq=freqs;
				this.TL=TL;
			}
			
			public float[] call(){
				DataPair dp=getValidData(tims,data,disx,disy);
				
				 long[] nt=dp.tims;
				float[] nd=dp.data;
				float[] nx=dp.disx;
				float[] ny=dp.disy;
				
				int rlen=hasST?freq.length*2+1:freq.length*2;
				
				if(nd.length<=1||ArrayUtil.getRange(nd)==0f||StatisticsUtil.cVariance(nd)==0f){
					// case where length of array is too short or
					// case where the matrix is singular
					float[] re=new float[rlen];
					
					for(int l=0,L=re.length;l<L;l++) re[l]=Record.undef;
					
					return re;
					
				}else{
					GaussMarkovEstimator gme=hasST?new GMBlas(nd,nt,nx,ny):new GMBlas(nd,nt);
					gme.setFrequenciesAndTimescales(AutoCorrType.TCosExp,freq,TL);
					gme.estimateCycles(false);
					
					return gme.varianceContribution();
				}
			}
		}
		
		private final class SolverCoefficients implements Callable<float[]>{
			//
			float TL=0;
			
			 long[] tims=null;
			float[] data=null;
			float[] disx=null;
			float[] disy=null;
			float[] freq=null;
			
			SolverCoefficients(float[] data,long[] tims,float[] disx,float[] disy,float[] freqs,float TL){
				this.data=data;
				this.tims=tims;
				this.disx=disx;
				this.disy=disy;
				this.freq=freqs;
				this.TL=TL;
			}
			
			public float[] call(){
				DataPair dp=getValidData(tims,data,disx,disy);
				
				 long[] nt=dp.tims;
				float[] nd=dp.data;
				float[] nx=dp.disx;
				float[] ny=dp.disy;
				
				int rlen=hasST?freq.length*2+6:freq.length*2+1;
				
				if(nd.length<=1||ArrayUtil.getRange(nd)==0f||StatisticsUtil.cVariance(nd)==0f){
					// case where length of array is too short or
					// case where the matrix is singular
					float[] re=new float[rlen];
					
					for(int l=0,L=re.length;l<L;l++) re[l]=Record.undef;
					
					return re;
					
				}else{
					GaussMarkovEstimator gme=hasST?new GMBlas(nd,nt,nx,ny):new GMBlas(nd,nt);
					gme.setFrequenciesAndTimescales(AutoCorrType.TCosExp,freq,TL);
					gme.estimateCycles(false);
					
					return gme.getCoefficients(true);
				}
			}
		}
		
		private final class SolverResidual implements Runnable{
			//
			int noCycRem=0;
			float TL=0;
			
			 long[] tims=null;
			float[] data=null;
			float[] disx=null;
			float[] disy=null;
			float[] freq=null;
			
			SolverResidual(float[] data,long[] tims,float[] disx,float[] disy,float[] freqs,float TL,int noCycRem){
				this.data=data;
				this.tims=tims;
				this.disx=disx;
				this.disy=disy;
				this.freq=freqs;
				this.noCycRem=noCycRem;
				this.TL=TL;
			}
			
			public void run(){
				DataPair dp=getValidData(tims,data,disx,disy);
				
				 long[] nt=dp.tims;
				float[] nd=dp.data;
				float[] nx=dp.disx;
				float[] ny=dp.disy;
				
				if(nd.length<=2||ArrayUtil.getRange(nd)==0){
					// case where length of array is too short or
					// case where the matrix is singular
					for(int l=0,L=data.length;l<L;l++) data[l]=Record.undef;
					
				}else{
					GaussMarkovEstimator gme=hasST?new GMBlas(nd,nt,nx,ny):new GMBlas(nd,nt);
					gme.setFrequenciesAndTimescales(AutoCorrType.TCosExp,freq,TL);
					gme.estimateCycles(false);
					
					copyResidualBack(gme.getResidualData(noCycRem),data);
				}
			}
			
			void copyResidualBack(float[] residual,float[] dest){
				for(int i=0,ptr=0,I=dest.length;i<I;i++)
				if(dest[i]!=Record.undef){
					dest[i]=residual[ptr];
					ptr++;
				}
			}
		}
		
		private final static class DataPair{
			 long[] tims=null;
			float[] data=null;
			float[] disx=null;
			float[] disy=null;
			
			DataPair(long[] tims,float[] data,float[] disx,float[] disy){
				this.data=data;
				this.tims=tims;
				this.disx=disx;
				this.disy=disy;
			}
		}
	}
	
	
	/** test
	public static void main(String arg[]){
		float[] data=new float[]{3,4,1,2,4,6,8,5,3,2,5,8,5,0,4,1};
		 long[] tims=new  long[]{1,2,3,5,6,4,2,8,0,7,5,9,2,6,8,4};
		averagingSimultaneousData(data,tims);
	}*/
}
