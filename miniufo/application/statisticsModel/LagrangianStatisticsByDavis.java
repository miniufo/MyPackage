/**
 * @(#)LagrangianStatisticsByDavis.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.Predicate;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.AttachedMeta;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.StochasticModel;
import miniufo.util.Region2D;
import miniufo.util.TicToc;


/**
 * Lagrangian statistics
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LagrangianStatisticsByDavis extends SingleParticleStatistics{
	
	/**
	 * constructor, the lengths of Particles in the list should be the same
	 */
	public LagrangianStatisticsByDavis(List<? extends Particle> ls,DataDescriptor dd){ super(ls,dd);}
	
	
	/**
     * Compute the Lagrangian diffusivity using Davis' (1991) theory: <v'd'>.
     *
     * @param	cond	condition for a record to be the origin of pseudo-track
     * @param	tRad	maximum lead or lag count
     */
	public SingleParticleStatResult cStatisticsByDavisTheory(Predicate<Record> cond,int tRad){
		BinStatistics bd=new BinStatistics(tRad);
		
		bd.computeMean(cond);
		bd.computeAutoCovariance(cond);
		bd.computeDispersion(cond);
		bd.computeDiffByVD(cond);
		
		return bd.lsr;
	}
	
	/**
     * Compute the Lagrangian diffusivity using Taylor's (1922) theory: <v'v'>.
     *
     * @param	cond	condition for a record to be the origin of pseudo-track
     * @param	tRad	maximum lead or lag count
     */
	public SingleParticleStatResult cStatisticsByTaylorTheory(Predicate<Record> cond,int tRad){
		BinStatistics bd=new BinStatistics(tRad);
		
		bd.computeMean(cond);
		bd.computeAutoCovariance(cond);
		bd.computeDispersion(cond);
		bd.computeDiffByVV(cond);
		
		return bd.lsr;
	}
	
	/**
     * Compute the Lagrangian diffusivity using growth rate of dispersion: 0.5*d<d'd'>/dt.
     *
     * @param	cond	condition for a record to be the origin of pseudo-track
     * @param	tRad	maximum lead or lag count
     */
	public SingleParticleStatResult cStatisticsByDispersionTheory(Predicate<Record> cond,int tRad){
		BinStatistics bd=new BinStatistics(tRad);
		
		bd.computeMean(cond);
		bd.computeAutoCovariance(cond);
		bd.computeDispersion(cond);
		bd.computeDiffByDispDT(cond);
		
		return bd.lsr;
	}
	
	/**
     * Compute the Lagrangian diffusivity using three methods: <v'd'>, integral of <v'v'>, and 0.5*d<d'd'>/dt.
     *
     * @param	cond	condition for a record to be the origin of pseudo-track
     * @param	tRad	maximum lead or lag count
     */
	public SingleParticleStatResult[] cStatistics(Predicate<Record> cond,int tRad){
		BinStatistics bd1=new BinStatistics(tRad);
		
		bd1.computeMean(cond);
		bd1.computeAutoCovariance(cond);
		bd1.computeDispersion(cond);
		
		BinStatistics bd2=bd1.copy();
		BinStatistics bd3=bd1.copy();
		
		bd1.computeDiffByVD(cond);
		bd2.computeDiffByVV(cond);
		bd3.computeDiffByDispDT(cond);
		
		return new SingleParticleStatResult[]{bd1.lsr,bd2.lsr,bd3.lsr};
	}
	
	
	/**
     * Compute Lagrangian statistics maps using time mean within the given time lags [str, end].
     *
     * @param	tRad		maximum time lags
     * @param	bRad		radius for bins (degree)
     * @param	str			start index
     * @param	end			end index
     * @param	minTracks	minimum tracks within the bin, else undef value is assigned to result
     * 
     * @return	stats		[0,1] for [Kxx,Kyy], [2,3] for [Tx,Ty] and [4,5] for [Lx,Ly]
     */
	public Variable[] cMeanStatisticsMapByDavisTheory(int tRad,float bRad,int str,int end,int minTracks){
		return cStatisticsMapByDavisTheory(tRad,bRad,str,end,minTracks,true);
	}
	
	public Variable[] cMeanStatisticsMapByTaylorTheory(int tRad,float bRad,int str,int end,int minTracks){
		return cStatisticsMapByTaylorTheory(tRad,bRad,str,end,minTracks,true);
	}
	
	public Variable[] cMeanStatisticsMapByDispersionTheory(int tRad,float bRad,int str,int end,int minTracks){
		return cStatisticsMapByDispersionTheory(tRad,bRad,str,end,minTracks,true);
	}
	
	public Variable[][] cMeanStatisticsMap(int tRad,float bRad,int str,int end,int minTracks,AttachedMeta ypos){
		return cStatisticsMap(tRad,bRad,str,end,minTracks,true,ypos);
	}
	
	
	/**
     * Compute Lagrangian statistics maps using time mean within the given time lags [str, end].
     *
     * @param	tRad		maximum time lags
     * @param	bRad		radius for bins (degree)
     * @param	str			start index
     * @param	end			end index
     * @param	minTracks	minimum tracks within the bin, else undef value is assigned to result
     * 
     * @return	stats		[0,1] for [Kxx,Kyy], [2,3] for [Tx,Ty] and [4,5] for [Lx,Ly]
     */
	public Variable[] cMaxStatisticsMapByDavisTheory(int tRad,float bRad,int str,int end,int minTracks){
		return cStatisticsMapByDavisTheory(tRad,bRad,str,end,minTracks,false);
	}
	
	public Variable[] cMaxStatisticsMapByTaylorTheory(int tRad,float bRad,int str,int end,int minTracks){
		return cStatisticsMapByTaylorTheory(tRad,bRad,str,end,minTracks,false);
	}
	
	public Variable[] cMaxStatisticsMapByDispersionTheory(int tRad,float bRad,int str,int end,int minTracks){
		return cStatisticsMapByDispersionTheory(tRad,bRad,str,end,minTracks,false);
	}
	
	public Variable[][] cMaxStatisticsMap(int tRad,float bRad,int str,int end,int minTracks,AttachedMeta ypos){
		return cStatisticsMap(tRad,bRad,str,end,minTracks,false,ypos);
	}
	
	
	/*** helper method and class ***/
	
	/**
     * Compute Lagrangian statistics maps using mean/maximum value within the given time lags [str, end].
     *
     * @param	tRad		maximum time lags
     * @param	bRad		radius for bins (degree)
     * @param	str			start index
     * @param	end			end index
     * @param	minTracks	minimum tracks within the bin, else undef value is assigned to result
     * @param	ave			mean method or maximum method to be used
     * 
     * @return	stats		[0,1] for [Kxx,Kyy], [2,3] for [Tx,Ty] and [4,5] for [Lx,Ly]
     */
	private Variable[] cStatisticsMapByDavisTheory(int tRad,float bRad,int str,int end,int minTracks,boolean ave){
		float undef=dd.getUndef(null);
		
		Variable[] stats=new Variable[9];
		
		Range r=new Range("",dd);
		
		stats[0]=new Variable("Kxx",r);	stats[0].setUndef(undef);
		stats[1]=new Variable("Kyy",r);	stats[1].setUndef(undef);
		stats[2]=new Variable("Tx" ,r);	stats[2].setUndef(undef);
		stats[3]=new Variable("Ty" ,r);	stats[3].setUndef(undef);
		stats[4]=new Variable("Lx" ,r);	stats[4].setUndef(undef);
		stats[5]=new Variable("Ly" ,r);	stats[5].setUndef(undef);
		stats[6]=new Variable("K11",r);	stats[6].setUndef(undef);
		stats[7]=new Variable("K22",r);	stats[7].setUndef(undef);
		stats[8]=new Variable("ang",r);	stats[8].setUndef(undef);
		
		stats[0].setCommentAndUnit("zonal component of diffusivity (10^7 cm^2/s)");
		stats[1].setCommentAndUnit("meridional component of diffusivity (10^7 cm^2/s)");
		stats[2].setCommentAndUnit("zonal component of Lagrangian integral timescale (day)");
		stats[3].setCommentAndUnit("meridional component of Lagrangian integral timescale (day)");
		stats[4].setCommentAndUnit("zonal component of Lagrangian length scale (km)");
		stats[5].setCommentAndUnit("meridional component of Lagrangian length scale (km)");
		stats[6].setCommentAndUnit("major component of diffusivity (10^7 cm^2/s)");
		stats[7].setCommentAndUnit("minor component of diffusivity (10^7 cm^2/s)");
		stats[8].setCommentAndUnit("angle of major component (radian)");
		
		float[][] kxdata=stats[0].getData()[0][0];
		float[][] kydata=stats[1].getData()[0][0];
		float[][] txdata=stats[2].getData()[0][0];
		float[][] tydata=stats[3].getData()[0][0];
		float[][] lxdata=stats[4].getData()[0][0];
		float[][] lydata=stats[5].getData()[0][0];
		float[][] k1data=stats[6].getData()[0][0];
		float[][] k2data=stats[7].getData()[0][0];
		float[][] agdata=stats[8].getData()[0][0];
		
		float[] lons=dd.getXDef().getSamples();
		float[] lats=dd.getYDef().getSamples();
		
		List<Future<float[]>> ls=new ArrayList<>(dd.getXCount()-1);
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
		
		TicToc.tic("start computing statistics grid by grid");
		
		float per=dd.getYCount()/20f;
		float curr=0;
		
		for(int j=0,J=dd.getYCount()-1;j<J;j++){
			for(int i=0,I=dd.getXCount()-1;i<I;i++){
				final int itag=i;
				final int jtag=j;
				
				Predicate<Record> cond=rec->
				new Region2D(lons[itag]-bRad,lats[jtag]-bRad,lons[itag]+bRad,lats[jtag]+bRad).inRange(rec.getXPos(),rec.getYPos());
				
				if(ave)
					ls.add(cs.submit(()->cStatisticsByDavisTheory(cond,tRad).getMean(str,end,minTracks)));
				else
					ls.add(cs.submit(()->cStatisticsByDavisTheory(cond,tRad).getMax(str,end,minTracks)));
			}
			
			try{
				for(int i=0,ptr=0,I=dd.getXCount()-1;i<I;i++){
					float[] mean=ls.get(ptr++).get();
					
					if(mean!=null){
						kxdata[j][i]=mean[0];	kydata[j][i]=mean[1];
						txdata[j][i]=mean[2];	tydata[j][i]=mean[3];
						lxdata[j][i]=mean[4];	lydata[j][i]=mean[5];
						k1data[j][i]=mean[6];	k2data[j][i]=mean[7];
						agdata[j][i]=mean[8];
						
					}else{
						kxdata[j][i]=undef;	kydata[j][i]=undef;
						txdata[j][i]=undef;	tydata[j][i]=undef;
						lxdata[j][i]=undef;	lydata[j][i]=undef;
						k1data[j][i]=undef;	k2data[j][i]=undef;
						agdata[j][i]=undef;
					}
				}
			}
			catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
			catch(ExecutionException   e){ e.printStackTrace(); System.exit(0);}
			
			if(j-curr>=per){
				curr=j;
				System.out.print(".");
			}
			
			ls.clear();
		}
		
		TicToc.toc(TimeUnit.MINUTES);
		
		return stats;
	}
	
	private Variable[] cStatisticsMapByTaylorTheory(int tRad,float bRad,int str,int end,int minTracks,boolean ave){
		float undef=dd.getUndef(null);
		
		Variable[] stats=new Variable[9];
		
		Range r=new Range("",dd);
		
		stats[0]=new Variable("Kxx",r);	stats[0].setUndef(undef);
		stats[1]=new Variable("Kyy",r);	stats[1].setUndef(undef);
		stats[2]=new Variable("Tx" ,r);	stats[2].setUndef(undef);
		stats[3]=new Variable("Ty" ,r);	stats[3].setUndef(undef);
		stats[4]=new Variable("Lx" ,r);	stats[4].setUndef(undef);
		stats[5]=new Variable("Ly" ,r);	stats[5].setUndef(undef);
		stats[6]=new Variable("K11",r);	stats[6].setUndef(undef);
		stats[7]=new Variable("K22",r);	stats[7].setUndef(undef);
		stats[8]=new Variable("ang",r);	stats[8].setUndef(undef);
		
		stats[0].setCommentAndUnit("zonal component of diffusivity (10^7 cm^2/s)");
		stats[1].setCommentAndUnit("meridional component of diffusivity (10^7 cm^2/s)");
		stats[2].setCommentAndUnit("zonal component of Lagrangian integral timescale (day)");
		stats[3].setCommentAndUnit("meridional component of Lagrangian integral timescale (day)");
		stats[4].setCommentAndUnit("zonal component of Lagrangian length scale (km)");
		stats[5].setCommentAndUnit("meridional component of Lagrangian length scale (km)");
		stats[6].setCommentAndUnit("major component of diffusivity (10^7 cm^2/s)");
		stats[7].setCommentAndUnit("minor component of diffusivity (10^7 cm^2/s)");
		stats[8].setCommentAndUnit("angle of major component (radian)");
		
		float[][] kxdata=stats[0].getData()[0][0];
		float[][] kydata=stats[1].getData()[0][0];
		float[][] txdata=stats[2].getData()[0][0];
		float[][] tydata=stats[3].getData()[0][0];
		float[][] lxdata=stats[4].getData()[0][0];
		float[][] lydata=stats[5].getData()[0][0];
		float[][] k1data=stats[6].getData()[0][0];
		float[][] k2data=stats[7].getData()[0][0];
		float[][] agdata=stats[8].getData()[0][0];
		
		float[] lons=dd.getXDef().getSamples();
		float[] lats=dd.getYDef().getSamples();
		
		List<Future<float[]>> ls=new ArrayList<>(dd.getXCount()-1);
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
		
		TicToc.tic("start computing statistics grid by grid");
		
		float per=dd.getYCount()/20f;
		float curr=0;
		
		for(int j=0,J=dd.getYCount()-1;j<J;j++){
			for(int i=0,I=dd.getXCount()-1;i<I;i++){
				final int itag=i;
				final int jtag=j;
				
				Predicate<Record> cond=rec->
				new Region2D(lons[itag]-bRad,lats[jtag]-bRad,lons[itag]+bRad,lats[jtag]+bRad).inRange(rec.getXPos(),rec.getYPos());
				
				if(ave)
					ls.add(cs.submit(()->cStatisticsByTaylorTheory(cond,tRad).getMean(str,end,minTracks)));
				else
					ls.add(cs.submit(()->cStatisticsByTaylorTheory(cond,tRad).getMax(str,end,minTracks)));
			}
			
			try{
				for(int i=0,ptr=0,I=dd.getXCount()-1;i<I;i++){
					float[] mean=ls.get(ptr++).get();
					
					if(mean!=null){
						kxdata[j][i]=mean[0];	kydata[j][i]=mean[1];
						txdata[j][i]=mean[2];	tydata[j][i]=mean[3];
						lxdata[j][i]=mean[4];	lydata[j][i]=mean[5];
						k1data[j][i]=mean[6];	k2data[j][i]=mean[7];
						agdata[j][i]=mean[8];
						
					}else{
						kxdata[j][i]=undef;	kydata[j][i]=undef;
						txdata[j][i]=undef;	tydata[j][i]=undef;
						lxdata[j][i]=undef;	lydata[j][i]=undef;
						k1data[j][i]=undef;	k2data[j][i]=undef;
						agdata[j][i]=undef;
					}
				}
			}
			catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
			catch(ExecutionException   e){ e.printStackTrace(); System.exit(0);}
			
			if(j-curr>=per){
				curr=j;
				System.out.print(".");
			}
			
			ls.clear();
		}
		
		TicToc.toc(TimeUnit.MINUTES);
		
		return stats;
	}
	
	private Variable[] cStatisticsMapByDispersionTheory(int tRad,float bRad,int str,int end,int minTracks,boolean ave){
		float undef=dd.getUndef(null);
		
		Variable[] stats=new Variable[9];
		
		Range r=new Range("",dd);
		
		stats[0]=new Variable("Kxx",r);	stats[0].setUndef(undef);
		stats[1]=new Variable("Kyy",r);	stats[1].setUndef(undef);
		stats[2]=new Variable("Tx" ,r);	stats[2].setUndef(undef);
		stats[3]=new Variable("Ty" ,r);	stats[3].setUndef(undef);
		stats[4]=new Variable("Lx" ,r);	stats[4].setUndef(undef);
		stats[5]=new Variable("Ly" ,r);	stats[5].setUndef(undef);
		stats[6]=new Variable("K11",r);	stats[6].setUndef(undef);
		stats[7]=new Variable("K22",r);	stats[7].setUndef(undef);
		stats[8]=new Variable("ang",r);	stats[8].setUndef(undef);
		
		stats[0].setCommentAndUnit("zonal component of diffusivity (10^7 cm^2/s)");
		stats[1].setCommentAndUnit("meridional component of diffusivity (10^7 cm^2/s)");
		stats[2].setCommentAndUnit("zonal component of Lagrangian integral timescale (day)");
		stats[3].setCommentAndUnit("meridional component of Lagrangian integral timescale (day)");
		stats[4].setCommentAndUnit("zonal component of Lagrangian length scale (km)");
		stats[5].setCommentAndUnit("meridional component of Lagrangian length scale (km)");
		stats[6].setCommentAndUnit("major component of diffusivity (10^7 cm^2/s)");
		stats[7].setCommentAndUnit("minor component of diffusivity (10^7 cm^2/s)");
		stats[8].setCommentAndUnit("angle of major component (radian)");
		
		float[][] kxdata=stats[0].getData()[0][0];
		float[][] kydata=stats[1].getData()[0][0];
		float[][] txdata=stats[2].getData()[0][0];
		float[][] tydata=stats[3].getData()[0][0];
		float[][] lxdata=stats[4].getData()[0][0];
		float[][] lydata=stats[5].getData()[0][0];
		float[][] k1data=stats[6].getData()[0][0];
		float[][] k2data=stats[7].getData()[0][0];
		float[][] agdata=stats[8].getData()[0][0];
		
		float[] xdef=dd.getXDef().getSamples();
		float[] ydef=dd.getYDef().getSamples();
		
		List<Future<float[]>> ls=new ArrayList<>(dd.getXCount()-1);
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<float[]> cs=new ExecutorCompletionService<>(es);
		
		TicToc.tic("start computing statistics grid by grid");
		
		float per=dd.getYCount()/20f;
		float curr=0;
		
		for(int j=0,J=dd.getYCount()-1;j<J;j++){
			for(int i=0,I=dd.getXCount()-1;i<I;i++){
				final int itag=i;
				final int jtag=j;
				
				Predicate<Record> cond=rec->
				new Region2D(xdef[itag]-bRad,ydef[jtag]-bRad,xdef[itag]+bRad,ydef[jtag]+bRad).inRange(rec.getXPos(),rec.getYPos());
				
				if(ave)
					ls.add(cs.submit(()->cStatisticsByDispersionTheory(cond,tRad).getMean(str,end,minTracks)));
				else
					ls.add(cs.submit(()->cStatisticsByDispersionTheory(cond,tRad).getMax(str,end,minTracks)));
			}
			
			try{
				for(int i=0,ptr=0,I=dd.getXCount()-1;i<I;i++){
					float[] mean=ls.get(ptr++).get();
					
					if(mean!=null){
						kxdata[j][i]=mean[0];	kydata[j][i]=mean[1];
						txdata[j][i]=mean[2];	tydata[j][i]=mean[3];
						lxdata[j][i]=mean[4];	lydata[j][i]=mean[5];
						k1data[j][i]=mean[6];	k2data[j][i]=mean[7];
						agdata[j][i]=mean[8];
						
					}else{
						kxdata[j][i]=undef;	kydata[j][i]=undef;
						txdata[j][i]=undef;	tydata[j][i]=undef;
						lxdata[j][i]=undef;	lydata[j][i]=undef;
						k1data[j][i]=undef;	k2data[j][i]=undef;
						agdata[j][i]=undef;
					}
				}
			}
			catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
			catch(ExecutionException   e){ e.printStackTrace(); System.exit(0);}
			
			if(j-curr>=per){
				curr=j;
				System.out.print(".");
			}
			
			ls.clear();
		}
		
		TicToc.toc(TimeUnit.MINUTES);
		
		return stats;
	}
	
	private Variable[][] cStatisticsMap(int tRad,float bRad,int str,int end,int minTracks,boolean ave,AttachedMeta ypos){
		float undef=dd.getUndef(null);
		
		Variable[][] stats=new Variable[3][9];
		
		Range r=new Range("",dd);
		
		float[][][] kxdata=new float[3][][];
		float[][][] kydata=new float[3][][];
		float[][][] txdata=new float[3][][];
		float[][][] tydata=new float[3][][];
		float[][][] lxdata=new float[3][][];
		float[][][] lydata=new float[3][][];
		float[][][] k1data=new float[3][][];
		float[][][] k2data=new float[3][][];
		float[][][] agdata=new float[3][][];
		
		for(int i=0;i<3;i++){
			stats[i][0]=new Variable("Kxx"+(i+1),r);	stats[i][0].setUndef(undef);
			stats[i][1]=new Variable("Kyy"+(i+1),r);	stats[i][1].setUndef(undef);
			stats[i][2]=new Variable("Tx" +(i+1),r);	stats[i][2].setUndef(undef);
			stats[i][3]=new Variable("Ty" +(i+1),r);	stats[i][3].setUndef(undef);
			stats[i][4]=new Variable("Lx" +(i+1),r);	stats[i][4].setUndef(undef);
			stats[i][5]=new Variable("Ly" +(i+1),r);	stats[i][5].setUndef(undef);
			stats[i][6]=new Variable("K11"+(i+1),r);	stats[i][6].setUndef(undef);
			stats[i][7]=new Variable("K22"+(i+1),r);	stats[i][7].setUndef(undef);
			stats[i][8]=new Variable("ang"+(i+1),r);	stats[i][8].setUndef(undef);
			
			stats[i][0].setCommentAndUnit("zonal component of diffusivity (10^7 cm^2/s)");
			stats[i][1].setCommentAndUnit("meridional component of diffusivity (10^7 cm^2/s)");
			stats[i][2].setCommentAndUnit("zonal component of Lagrangian integral timescale (day)");
			stats[i][3].setCommentAndUnit("meridional component of Lagrangian integral timescale (day)");
			stats[i][4].setCommentAndUnit("zonal component of Lagrangian length scale (km)");
			stats[i][5].setCommentAndUnit("meridional component of Lagrangian length scale (km)");
			stats[i][6].setCommentAndUnit("major component of diffusivity (10^7 cm^2/s)");
			stats[i][7].setCommentAndUnit("minor component of diffusivity (10^7 cm^2/s)");
			stats[i][8].setCommentAndUnit("angle of major component (radian)");
			
			kxdata[i]=stats[i][0].getData()[0][0];
			kydata[i]=stats[i][1].getData()[0][0];
			txdata[i]=stats[i][2].getData()[0][0];
			tydata[i]=stats[i][3].getData()[0][0];
			lxdata[i]=stats[i][4].getData()[0][0];
			lydata[i]=stats[i][5].getData()[0][0];
			k1data[i]=stats[i][6].getData()[0][0];
			k2data[i]=stats[i][7].getData()[0][0];
			agdata[i]=stats[i][8].getData()[0][0];
		}
		
		float[] xdef=dd.getXDef().getSamples();
		float[] ydef=dd.getYDef().getSamples();
		
		List<Future<float[][]>> ls=new ArrayList<>(dd.getXCount()-1);
		ExecutorService es=ConcurrentUtil.defaultExecutor();
		CompletionService<float[][]> cs=new ExecutorCompletionService<>(es);
		
		TicToc.tic("start computing statistics grid by grid");
		
		float per=dd.getYCount()/20f;
		float curr=0;
		
		for(int j=0,J=dd.getYCount()-1;j<J;j++){
			for(int i=0,I=dd.getXCount()-1;i<I;i++){
				final int itag=i;
				final int jtag=j;
				
				Predicate<Record> cond=rec->
				new Region2D(xdef[itag]-bRad,ydef[jtag]-bRad,xdef[itag]+bRad,ydef[jtag]+bRad).inRange(rec.getXPos(),rec.getData(ypos));
				
				if(ave)
					ls.add(cs.submit(()->{
						SingleParticleStatResult[] re=cStatistics(cond,tRad);
						return new float[][]{
							re[0].getMean(str,end,minTracks),
							re[1].getMean(str,end,minTracks),
							re[2].getMean(str,end,minTracks)
						};
					}));
				else
					ls.add(cs.submit(()->{
						SingleParticleStatResult[] re=cStatistics(cond,tRad);
						return new float[][]{
							re[0].getMax(str,end,minTracks),
							re[1].getMax(str,end,minTracks),
							re[2].getMax(str,end,minTracks)
						};
					}));
			}
			
			try{
				for(int m=0;m<3;m++)
				for(int i=0,ptr=0,I=dd.getXCount()-1;i<I;i++){
					float[][] mean=ls.get(ptr++).get();
					
					if(mean[m]!=null){
						kxdata[m][j][i]=mean[m][0];	kydata[m][j][i]=mean[m][1];
						txdata[m][j][i]=mean[m][2];	tydata[m][j][i]=mean[m][3];
						lxdata[m][j][i]=mean[m][4];	lydata[m][j][i]=mean[m][5];
						k1data[m][j][i]=mean[m][6];	k2data[m][j][i]=mean[m][7];
						agdata[m][j][i]=mean[m][8];
						
					}else{
						kxdata[m][j][i]=undef;	kydata[m][j][i]=undef;
						txdata[m][j][i]=undef;	tydata[m][j][i]=undef;
						lxdata[m][j][i]=undef;	lydata[m][j][i]=undef;
						k1data[m][j][i]=undef;	k2data[m][j][i]=undef;
						agdata[m][j][i]=undef;
					}
				}
			}
			catch(InterruptedException|ExecutionException e){ e.printStackTrace(); System.exit(0);}
			
			if(j-curr>=per){
				curr=j;
				System.out.print(".");
			}
			
			ls.clear();
		}
		
		TicToc.toc(TimeUnit.MINUTES);
		
		return stats;
	}
	
	
	private final class BinStatistics{
		//
		private int pseudoTracks=0;	// number of pseudo-tracks, i.e. observations at lag 0
		private int noOfMaxLag  =0;	// observations at maximum lag
		private int noOfMinLag  =0;	// observations at minimum lag
		
		private SingleParticleStatResult lsr=null;
		
		
		/**
		 * Constructor.
		 */
		BinStatistics(int tRad){
			lsr=new SingleParticleStatResult(tRad);
			
			Particle p1=null;
			
			for(Particle p:ls) if(p.getTCount()>1) p1=p;
			
			if(p1==null) throw new IllegalArgumentException("no valid particle data");
			
			lsr.dt=new MDate(p1.getTime(1)).getDT(new MDate(p1.getTime(0)));
			
			//writeTracks("/lustre/home/qianyk/Data/GDP/Track/diff9regions/"+lon1+lat1);
		}
		
		
		/**
		 * Copy the current BinStatistics to a new one.
		 */
		BinStatistics copy(){
			BinStatistics bs=new BinStatistics();
			
			bs.pseudoTracks=pseudoTracks;
			bs.noOfMaxLag=noOfMaxLag;
			bs.noOfMinLag=noOfMinLag;
			bs.lsr=lsr.copy();
			
			return bs;
		}
		
		
		// compute Lagrangian mean as a function of tau
		void computeMean(Predicate<Record> cond){
			pseudoTracks=0;
			noOfMaxLag  =0;
			noOfMinLag  =0;
			
			int tRad=lsr.tRad;
			
			float umn=0,vmn=0;
			
			Averager av=new Averager(4,tRad);
			
			for(Particle p:ls){
				float[] uspd=p.getUVel();
				float[] vspd=p.getVVel();
				
				for(int l=0,L=p.getTCount();l<L;l++){
					float oX=p.getXPosition(l);
					float oY=p.getYPosition(l);
					
					// meet condition to be an origin of pseudo-track
					if(!cond.test(p.getRecord(l))) continue;
					if(uspd[l]==Record.undef) continue;
					
					umn+=uspd[l]; vmn+=vspd[l]; pseudoTracks++;
					
					for(int ll=0;ll<L;ll++){
						int tau=ll-l;
						
						if(tau>tRad||tau<-tRad) continue;
						if(uspd[ll]==Record.undef) continue;
						if(tau== tRad) noOfMaxLag++;
						if(tau==-tRad) noOfMinLag++;
						
						if(llpos){
							float dlon=(float)Math.toRadians(p.getXPosition(ll)-oX);
							float dlat=(float)Math.toRadians(p.getYPosition(ll)-oY);
							
							float disX=SpatialModel.REarth*dlon*(float)Math.cos(oY*Math.PI/180.0);
							float disY=SpatialModel.REarth*dlat;
							
							av.addSample(tau,uspd[ll],vspd[ll],disX,disY);
							
						}else{
							av.addSample(tau,uspd[ll],vspd[ll],p.getXPosition(ll)-oX,p.getYPosition(ll)-oY);
						}
					}
				}
			}
			
			lsr.umn=umn/pseudoTracks;
			lsr.vmn=vmn/pseudoTracks;
			
			lsr.pseudoTracks=pseudoTracks;
			lsr.noOfMaxLag  =noOfMaxLag;
			lsr.noOfMinLag  =noOfMinLag;
			
			av.average();
			
			for(int l=0,L=2*tRad+1;l<L;l++){
				lsr.num[l]=av.count[l];
				lsr.um [l]=(float)av.means[0][l];
				lsr.vm [l]=(float)av.means[1][l];
				lsr.DXm[l]=(float)av.means[2][l];
				lsr.DYm[l]=(float)av.means[3][l];
			}
		}
		
		// compute autocovariance of velocity Pxx, Pyy, Pxy and Pyx
		void computeAutoCovariance(Predicate<Record> cond){
			int tRad=lsr.tRad;
			
			boolean hasAcc=ls.get(0).getRecord(0).getDataLength()==4;
			
			// for Pxx Pxy Pyx Pyy, Qxx Qxy Qyx Qyy, ua va
			Averager av=new Averager(hasAcc?10:6,tRad);
			
			for(Particle p:ls){
				float[] uspd=p.getUVel();
				float[] vspd=p.getVVel();
				float[] accx=hasAcc?p.getAttachedData(StochasticModel.AccX):null;
				float[] accy=hasAcc?p.getAttachedData(StochasticModel.AccY):null;
				
				for(int l=0,L=p.getTCount();l<L;l++){
					// meet condition to be an origin of pseudo-track
					if(!cond.test(p.getRecord(l))) continue;
					if(uspd[l]==Record.undef) continue;
					
					float ua0=uspd[l]-lsr.um[tRad];	// u'(tau=0)
					float va0=vspd[l]-lsr.vm[tRad];	// v'(tau=0)
					
					for(int ll=0;ll<L;ll++){
						int tau=ll-l;
						
						if(tau>tRad||tau<-tRad) continue;
						if(uspd[ll]==Record.undef) continue;
						
						float ua=uspd[ll]-lsr.um[tRad+tau];	// u'(tau)
						float va=vspd[ll]-lsr.vm[tRad+tau];	// v'(tau)
						
						float Pxx=ua0*ua;	// u'(0)*u'(tau)
						float Pyy=va0*va;	// v'(0)*v'(tau)
						float Pxy=ua0*va;	// u'(0)*v'(tau)
						float Pyx=va0*ua;	// v'(0)*u'(tau)
						
						if(hasAcc){
							float Qxx=accx[l]*accx[ll];	// ax'(0)*ax'(tau)
							float Qyy=accy[l]*accy[ll];	// ay'(0)*ay'(tau)
							float Qxy=accx[l]*accy[ll];	// ax'(0)*ay'(tau)
							float Qyx=accy[l]*accx[ll];	// ay'(0)*ax'(tau)
							
							av.addSample(tau,Pxx,Pyy,Pxy,Pyx,ua,va,Qxx,Qyy,Qxy,Qyx);
							
						}else av.addSample(tau,Pxx,Pyy,Pxy,Pyx,ua,va);
					}
				}
			}
			
			av.average();
			
			for(int l=0,L=2*tRad+1;l<L;l++){
				lsr.num[l]=av.count[l];
				
				lsr.Pxx[l]=(float)av.means[0][l];
				lsr.Pyy[l]=(float)av.means[1][l];
				lsr.Pxy[l]=(float)av.means[2][l];
				lsr.Pyx[l]=(float)av.means[3][l];
				lsr.ua [l]=(float)av.means[4][l];
				lsr.va [l]=(float)av.means[5][l];
				
				if(hasAcc){
					lsr.Qxx[l]=(float)av.means[6][l];
					lsr.Qyy[l]=(float)av.means[7][l];
					lsr.Qxy[l]=(float)av.means[8][l];
					lsr.Qyx[l]=(float)av.means[9][l];
				}
			}
		}
		
		// compute dispersion by <d'd'> as a function of tau
		void computeDispersion(Predicate<Record> cond){
			int tRad=lsr.tRad;
			
			Averager av=new Averager(4,tRad);
			
			for(Particle p:ls){
				float[] uspd=p.getUVel();
				
				for(int l=0,L=p.getTCount();l<L;l++){
					float oX=p.getXPosition(l);
					float oY=p.getYPosition(l);
					
					// meet condition to be an origin of pseudo-track
					if(!cond.test(p.getRecord(l))) continue;
					if(uspd[l]==Record.undef) continue;
					
					for(int ll=0;ll<L;ll++){
						int tau=ll-l;
						
						if(tau>tRad||tau<-tRad) continue;
						if(uspd[ll]==Record.undef) continue;
						
						if(llpos){
							float dlon=(float)Math.toRadians(p.getXPosition(ll)-oX);
							float dlat=(float)Math.toRadians(p.getYPosition(ll)-oY);
							
							// dx'(tau)
							float disX=SpatialModel.REarth*dlon*(float)Math.cos(oY*Math.PI/180.0)-lsr.DXm[tRad+tau];
							// dy'(tau)
							float disY=SpatialModel.REarth*dlat-lsr.DYm[tRad+tau];
							
							float Dxx=disX*disX;	// dx'(tau)*dx'(tau)
							float Dyy=disY*disY;	// dy'(tau)*dy'(tau)
							float Dxy=disX*disY;	// dx'(tau)*dy'(tau)
							float Dyx=disY*disX;	// dy'(tau)*dx'(tau)
							
							av.addSample(tau,Dxx,Dyy,Dxy,Dyx);
							
						}else{
							// dx'(tau)
							float disX=p.getXPosition(ll)-oX-lsr.DXm[tRad+tau];
							// dy'(tau)
							float disY=p.getYPosition(ll)-oY-lsr.DYm[tRad+tau];
							
							float Dxx=disX*disX;	// dx'(tau)*dx'(tau)
							float Dyy=disY*disY;	// dy'(tau)*dy'(tau)
							float Dxy=disX*disY;	// dx'(tau)*dy'(tau)
							float Dyx=disY*disX;	// dy'(tau)*dx'(tau)
							
							av.addSample(tau,Dxx,Dyy,Dxy,Dyx);
						}
					}
				}
			}
			
			av.average();
			
			for(int l=0,L=2*tRad+1;l<L;l++){
				lsr.num[l]=av.count[l];
				
				lsr.Dxx[l]=(float)av.means[0][l];
				lsr.Dyy[l]=(float)av.means[1][l];
				lsr.Dxy[l]=(float)av.means[2][l];
				lsr.Dyx[l]=(float)av.means[3][l];
			}
		}
		
		// compute diffusivity by -<v'(0)d'(-tau)>
		void computeDiffByVD(Predicate<Record> cond){
			int tRad=lsr.tRad,count=0;
			
			float ucfd=0,vcfd=0;
			
			Averager av=new Averager(4,tRad);
			Averager av2=new Averager(2,tRad);
			
			for(Particle p:ls){
				float[] uspd=p.getUVel();
				float[] vspd=p.getVVel();
				
				float usqr=0,vsqr=0;
				
				for(int l=0,L=p.getTCount();l<L;l++){
					float oX=p.getXPosition(l);
					float oY=p.getYPosition(l);
					
					// meet condition to be an origin of pseudo-track
					if(!cond.test(p.getRecord(l))) continue;
					if(uspd[l]==Record.undef) continue;
					
					float ua0=uspd[l]-lsr.um[tRad];	// u'(tau=0)
					float va0=vspd[l]-lsr.vm[tRad];	// v'(tau=0)
					
					usqr+=ua0*ua0;
					vsqr+=va0*va0;	count++;
					
					for(int ll=0;ll<L;ll++){
						int tau=ll-l;
						
						if(tau>tRad||tau<-tRad) continue;
						if(uspd[ll]==Record.undef) continue;
						
						if(llpos){
							float dlon=(float)Math.toRadians(p.getXPosition(ll)-oX);
							float dlat=(float)Math.toRadians(p.getYPosition(ll)-oY);
							
							// dx'(tau)
							float disX=SpatialModel.REarth*dlon*(float)Math.cos(oY*Math.PI/180.0)-lsr.DXm[tRad+tau];
							// dy'(tau)
							float disY=SpatialModel.REarth*dlat-lsr.DYm[tRad+tau];
							
							float Kxx=-ua0*disX;	// -u'(0)*dx'(tau)
							float Kyy=-va0*disY;	// -v'(0)*dy'(tau)
							float Kxy=-ua0*disY;	// -u'(0)*dy'(tau)
							float Kyx=-va0*disX;	// -v'(0)*dx'(tau)
							
							av.addSample(-tau,Kxx,Kyy,Kxy,Kyx);
							av2.addSample(tau,disX,disY);
							
						}else{
							// dx'(tau)
							float disX=p.getXPosition(ll)-oX-lsr.DXm[tRad+tau];
							// dy'(tau)
							float disY=p.getYPosition(ll)-oY-lsr.DYm[tRad+tau];
							
							float Kxx=-ua0*disX;	// -u'(0)*dx'(tau)
							float Kyy=-va0*disY;	// -v'(0)*dy'(tau)
							float Kxy=-ua0*disY;	// -u'(0)*dy'(tau)
							float Kyx=-va0*disX;	// -v'(0)*dx'(tau)
							
							av.addSample(-tau,Kxx,Kyy,Kxy,Kyx);
							av2.addSample(tau,disX,disY);
						}
					}
				}
				
				ucfd+=usqr;
				vcfd+=vsqr;
			}
			
			av.average();
			av2.average();
			
			for(int l=0,L=2*tRad+1;l<L;l++){
				lsr.num[l]=av2.count[l];
				
				lsr.Kxx[l]=(float)av.means[0][l];
				lsr.Kyy[l]=(float)av.means[1][l];
				lsr.Kxy[l]=(float)av.means[2][l];
				lsr.Kyx[l]=(float)av.means[3][l];
				
				lsr.DXa[l]=(float)av2.means[0][l];
				lsr.DYa[l]=(float)av2.means[1][l];
			}
			
			principalAxeDecomp();
			
			float[] scales=cScales();
			
			lsr.Tu=scales[0];	lsr.Lu=scales[2];	lsr.Ku=scales[4];
			lsr.Tv=scales[1];	lsr.Lv=scales[3];	lsr.Kv=scales[5];
			
			ucfd/=(count-1);	ucfd=(float)Math.sqrt(ucfd);
			vcfd/=(count-1);	vcfd=(float)Math.sqrt(vcfd);
			
			lsr.ucfd*=2f/Math.sqrt(count/(2f*lsr.Tu*4f));
			lsr.vcfd*=2f/Math.sqrt(count/(2f*lsr.Tv*4f));
		}
		
		// compute diffusivity by integrating <v'(0)v'(tau)> from -tau to 0
		void computeDiffByVV(Predicate<Record> cond){
			int tRad=lsr.tRad,len=2*tRad+1;
			
			float ucfd=0,vcfd=0;
			
			for(int l=0,half=len/2;l<half;l++){
				lsr.Kxx[l]=0;
				lsr.Kyy[l]=0;
				lsr.Kxy[l]=0;
				lsr.Kyx[l]=0;
				
				for(int ll=half,LL=len-l;ll<LL;ll++){
					lsr.Kxx[l]+=lsr.Pxx[ll];
					lsr.Kyy[l]+=lsr.Pyy[ll];
					lsr.Kxy[l]+=lsr.Pxy[ll];
					lsr.Kyx[l]+=lsr.Pyx[ll];
				}
				
				lsr.Kxx[l]-=(lsr.Pxx[half]+lsr.Pxx[len-1-l])/2;
				lsr.Kyy[l]-=(lsr.Pyy[half]+lsr.Pyy[len-1-l])/2;
				lsr.Kxy[l]-=(lsr.Pxy[half]+lsr.Pxy[len-1-l])/2;
				lsr.Kyx[l]-=(lsr.Pyx[half]+lsr.Pyx[len-1-l])/2;
				
				lsr.Kxx[l]*=-lsr.dt;
				lsr.Kyy[l]*=-lsr.dt;
				lsr.Kxy[l]*=-lsr.dt;
				lsr.Kyx[l]*=-lsr.dt;
			}
			
			for(int half=len/2,l=half+1,L=len;l<L;l++){
				lsr.Kxx[l]=0;
				lsr.Kyy[l]=0;
				lsr.Kxy[l]=0;
				lsr.Kyx[l]=0;
				
				for(int ll=half;ll>=len-1-l;ll--){
					lsr.Kxx[l]+=lsr.Pxx[ll];
					lsr.Kyy[l]+=lsr.Pyy[ll];
					lsr.Kxy[l]+=lsr.Pxy[ll];
					lsr.Kyx[l]+=lsr.Pyx[ll];
				}
				
				lsr.Kxx[l]-=(lsr.Pxx[half]+lsr.Pxx[len-1-l])/2;
				lsr.Kyy[l]-=(lsr.Pyy[half]+lsr.Pyy[len-1-l])/2;
				lsr.Kxy[l]-=(lsr.Pxy[half]+lsr.Pxy[len-1-l])/2;
				lsr.Kyx[l]-=(lsr.Pyx[half]+lsr.Pyx[len-1-l])/2;
				
				lsr.Kxx[l]*=lsr.dt;
				lsr.Kyy[l]*=lsr.dt;
				lsr.Kxy[l]*=lsr.dt;
				lsr.Kyx[l]*=lsr.dt;
			}
			
			principalAxeDecomp();
			
			float[] scales=cScales();
			
			lsr.Tu=scales[0];	lsr.Lu=scales[2];	lsr.Ku=scales[4];
			lsr.Tv=scales[1];	lsr.Lv=scales[3];	lsr.Kv=scales[5];
			
			ucfd/=(lsr.pseudoTracks-1);	ucfd=(float)Math.sqrt(ucfd);
			vcfd/=(lsr.pseudoTracks-1);	vcfd=(float)Math.sqrt(vcfd);
			
			lsr.ucfd*=2f/Math.sqrt(lsr.pseudoTracks/(2f*lsr.Tu*4f));
			lsr.vcfd*=2f/Math.sqrt(lsr.pseudoTracks/(2f*lsr.Tv*4f));
		}
		
		/*// compute diffusivity by differencing dispersion
		void computeDiffByDispDT(Predicate<Record> cond){
			int tRad=lsr.tRad,len=2*tRad+1;
			
			float ucfd=0,vcfd=0;
			
			computeDispersion(cond);
			
			for(int l=0;l<tRad;l++){ // negative time lags
				lsr.Kxx[l]=(lsr.Dxx[l+1]-lsr.Dxx[l])/lsr.dt/2f;
				lsr.Kyy[l]=(lsr.Dyy[l+1]-lsr.Dyy[l])/lsr.dt/2f;
				lsr.Kxy[l]=(lsr.Dxy[l+1]-lsr.Dxy[l])/lsr.dt/2f;
				lsr.Kyx[l]=(lsr.Dyx[l+1]-lsr.Dyx[l])/lsr.dt/2f;
			}
			
			for(int l=tRad+1;l<len;l++){ // positive time lags
				lsr.Kxx[l]=(lsr.Dxx[l]-lsr.Dxx[l-1])/lsr.dt/2f;
				lsr.Kyy[l]=(lsr.Dyy[l]-lsr.Dyy[l-1])/lsr.dt/2f;
				lsr.Kxy[l]=(lsr.Dxy[l]-lsr.Dxy[l-1])/lsr.dt/2f;
				lsr.Kyx[l]=(lsr.Dyx[l]-lsr.Dyx[l-1])/lsr.dt/2f;
			}
			
			principalAxeDecomp();
			
			float[] scales=cScales();
			
			lsr.Tu=scales[0];	lsr.Lu=scales[2];	lsr.Ku=scales[4];
			lsr.Tv=scales[1];	lsr.Lv=scales[3];	lsr.Kv=scales[5];
			
			ucfd/=(lsr.pseudoTracks-1);	ucfd=(float)Math.sqrt(ucfd);
			vcfd/=(lsr.pseudoTracks-1);	vcfd=(float)Math.sqrt(vcfd);
		}
		
		*/
		void computeDiffByDispDT(Predicate<Record> cond){
			int tRad=lsr.tRad,len=2*tRad+1;
			
			float ucfd=0,vcfd=0;
			
			lsr.Kxx[0]=(lsr.Dxx[1]-lsr.Dxx[0])/lsr.dt/2f;
			lsr.Kyy[0]=(lsr.Dyy[1]-lsr.Dyy[0])/lsr.dt/2f;
			lsr.Kxy[0]=(lsr.Dxy[1]-lsr.Dxy[0])/lsr.dt/2f;
			lsr.Kyx[0]=(lsr.Dyx[1]-lsr.Dyx[0])/lsr.dt/2f;
			
			for(int l=1,L=len-1;l<L;l++){
				lsr.Kxx[l]=(lsr.Dxx[l+1]-lsr.Dxx[l-1])/(2f*lsr.dt)/2f;
				lsr.Kyy[l]=(lsr.Dyy[l+1]-lsr.Dyy[l-1])/(2f*lsr.dt)/2f;
				lsr.Kxy[l]=(lsr.Dxy[l+1]-lsr.Dxy[l-1])/(2f*lsr.dt)/2f;
				lsr.Kyx[l]=(lsr.Dyx[l+1]-lsr.Dyx[l-1])/(2f*lsr.dt)/2f;
			}
			
			int last=len-1;
			
			lsr.Kxx[last]=(lsr.Dxx[last]-lsr.Dxx[last-1])/lsr.dt/2f;
			lsr.Kyy[last]=(lsr.Dyy[last]-lsr.Dyy[last-1])/lsr.dt/2f;
			lsr.Kxy[last]=(lsr.Dxy[last]-lsr.Dxy[last-1])/lsr.dt/2f;
			lsr.Kyx[last]=(lsr.Dyx[last]-lsr.Dyx[last-1])/lsr.dt/2f;
			
			principalAxeDecomp();
			
			float[] scales=cScales();
			
			lsr.Tu=scales[0];	lsr.Lu=scales[2];	lsr.Ku=scales[4];
			lsr.Tv=scales[1];	lsr.Lv=scales[3];	lsr.Kv=scales[5];
			
			ucfd/=(lsr.pseudoTracks-1);	ucfd=(float)Math.sqrt(ucfd);
			vcfd/=(lsr.pseudoTracks-1);	vcfd=(float)Math.sqrt(vcfd);
		}
		
		/*
		void writeTracks(String path){
			for(Particle p:ls){
				int inidx=-1,outidx=-1;
				
				for(int l=0,L=p.getTCount();l<L;l++){
					float lon=p.getLongitude(l);
					float lat=p.getLatitude(l);
					
					if(inRange(lon,lat)){ inidx=l; break;}
				}
				
				for(int l=p.getTCount()-1;l>=0;l--){
					float lon=p.getLongitude(l);
					float lat=p.getLatitude(l);
					
					if(inRange(lon,lat)){ outidx=l; break;}
				}
				
				if(inidx!=-1&&outidx!=-1){
					int str=inidx-lsr.tRad<0?0:inidx-lsr.tRad;
					int end=outidx+lsr.tRad>p.getTCount()-1?p.getTCount()-1:outidx+lsr.tRad;
					
					Particle pp=p.subRecord(str,end-str);
					pp.toTrajectoryFile(path+p.getID());
				}
			}
		}*/
		
		float[] cScales(){
			int tRad=lsr.tRad;
			
			float[] Tx=lsr.Kxx.clone();
			float[] Ty=lsr.Kyy.clone();
			float[] Lx=lsr.Kxx.clone();
			float[] Ly=lsr.Kyy.clone();
			
			for(int l=0,L=Tx.length;l<L;l++){
				Tx[l]/=lsr.Pxx[tRad];
				Ty[l]/=lsr.Pyy[tRad];
				Lx[l]/=Math.sqrt(lsr.Pxx[tRad]);
				Ly[l]/=Math.sqrt(lsr.Pyy[tRad]);
			}
			
			return new float[]{
				getMaxInPositiveRange(Tx)/86400f,	// day
				getMaxInPositiveRange(Ty)/86400f,
				getMaxInPositiveRange(Lx)/1000f,	// km
				getMaxInPositiveRange(Ly)/1000f,
				getMaxInPositiveRange(lsr.Kxx)/1e3f,	// 1e7 cm^2/s
				getMaxInPositiveRange(lsr.Kyy)/1e3f		// 1e7 cm^2/s
			};
		}
		
		
		/*** helper methods and class ***/
		private float getMaxInPositiveRange(float[] data){
			float re=Float.MIN_VALUE;
			
			for(int I=data.length,i=I/2;i<I;i++)
			if(data[i]>re) re=data[i];
			
			return re;
		}
		
		private void principalAxeDecomp(){
			int len=2*lsr.tRad+1;
			
			for(int l=0,L=len;l<L;l++){
				double KxxS=lsr.Kxx[l];
				double KyyS=lsr.Kyy[l];
				double KxyS=(lsr.Kxy[l]+lsr.Kyx[l])/2.0;
				double tmp=(KxxS+KyyS+Math.sqrt((KxxS-KyyS)*(KxxS-KyyS)+4.0*KxyS*KxyS))/2.0;
				
				lsr.K11[l]=(float)tmp;
				lsr.K22[l]=(float)(KxxS+KyyS-tmp);
				lsr.ang[l]=(float)(Math.atan2(2*KxyS,KxxS-KyyS)/2.0);
			}
			
			// change K11 and K22 to ensure K22 is smaller in magnitude
			float m11=0,m22=0,tmp=0;
			
			// for negative time lag
			for(int l=0,L=len/2;l<L;l++){
				m11+=lsr.K11[l];
				m22+=lsr.K22[l];
			}
			
			if(Math.abs(m11)<Math.abs(m22))
			for(int l=0,L=len/2;l<L;l++){
				tmp=lsr.K22[l];
				lsr.K22[l]=lsr.K11[l];
				lsr.K11[l]=tmp;
			}
			
			m11=0;m22=0;
			// for positive time lag
			for(int l=len/2;l<len;l++){
				m11+=lsr.K11[l];
				m22+=lsr.K22[l];
			}
			
			if(Math.abs(m11)<Math.abs(m22))
			for(int l=len/2;l<len;l++){
				tmp=lsr.K22[l];
				lsr.K22[l]=lsr.K11[l];
				lsr.K11[l]=tmp;
			}
		}
		
		private BinStatistics(){}
	}
	
	private static final class Averager{
		private int vcount =0;
		private int tRad   =0;
		
		private int[] count=null;
		
		private double[][] means=null;
		
		
		Averager(int vcount,int tRad){
			this.tRad  =tRad;
			this.vcount=vcount;
			
			count=new int[2*tRad+1];
			means=new double[vcount][2*tRad+1];
		}
		
		void addSample(int tau,float... data){
			if(data.length!=vcount) throw new IllegalArgumentException("lengths not equal");
			
			int tag=tau+tRad;
			
			for(int i=0;i<vcount;i++) means[i][tag]+=data[i];
			
			count[tag]++;
		}
		
		void average(){
			for(int i=0;i<vcount;i++)
			for(int l=0,L=2*tRad+1;l<L;l++) means[i][l]/=count[l];
		}
	}
	
	
	/** test
	public static void main(String arg[]){
		
	}*/
}
