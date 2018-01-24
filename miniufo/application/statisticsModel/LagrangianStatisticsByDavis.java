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
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
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
     * compute the Lagrangian diffusivity using <v'd'>
     *
     * @param	cond	condition for a record to be the origin of pseudo-track
     * @param	tRad	maximum lead or lag count
     */
	public SingleParticleStatResult cStatisticsByDavisTheory1(Predicate<Record> cond,int tRad){
		BinStatistics bd=new BinStatistics(tRad);
		
		bd.computeMean(cond);
		bd.computeAutoCovariance(cond);
		bd.computeDispersion(cond);
		bd.computeDiffByVD(cond);
		
		return bd.lsr;
	}
	
	/**
     * compute the Lagrangian diffusivity using <v'v'>
     *
     * @param	cond	condition for a record to be the origin of pseudo-track
     * @param	tRad	maximum lead or lag count
     */
	public SingleParticleStatResult cStatisticsByDavisTheory2(Predicate<Record> cond,int tRad){
		BinStatistics bd=new BinStatistics(tRad);
		
		bd.computeMean(cond);
		bd.computeDispersion(cond);
		bd.computeDiffByVV(cond);
		
		return bd.lsr;
	}
	
	/**
     * compute the Lagrangian diffusivity using 0.5*d<d'd'>/dt
     *
     * @param	cond	condition for a record to be the origin of pseudo-track
     * @param	tRad	maximum lead or lag count
     */
	public SingleParticleStatResult cStatisticsByDavisTheory3(Predicate<Record> cond,int tRad){
		BinStatistics bd=new BinStatistics(tRad);
		
		bd.computeMean(cond);
		bd.computeAutoCovariance(cond);
		bd.computeDiffByDispDT(cond);
		
		return bd.lsr;
	}
	
	
	/**
     * compute Lagrangian statistics map
     *
     * @param	tRad		maximum time lags
     * @param	bRad		radius for bins (degree)
     * @param	str			start index
     * @param	end			end index
     * @param	minTracks	minimum tracks within the bin, else undef value is assigned to result
     * 
     * @return	stats	[0,1] for [Kxx,Kyy], [2,3] for [Tx,Ty] and [4,5] for [Lx,Ly]
     */
	public Variable[] cStatisticsMapByDavisTheory1(int tRad,float bRad,int str,int end,int minTracks){
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
				
				ls.add(cs.submit(()->cStatisticsByDavisTheory1(cond,tRad).getMean(str,end,minTracks)));
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
	
	public Variable[] cStatisticsMapByDavisTheory2(int tRad,float bRad,int str,int end,int minTracks){
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
				
				ls.add(cs.submit(()->cStatisticsByDavisTheory2(cond,tRad).getMean(str,end,minTracks)));
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
	
	public Variable[] cStatisticsMapByDavisTheory3(int tRad,float bRad,int str,int end,int minTracks){
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
				
				ls.add(cs.submit(()->cStatisticsByDavisTheory3(cond,tRad).getMean(str,end,minTracks)));
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
	
	
	/*** helper method and class ***/
	private final class BinStatistics{
		//
		private int pseudoTracks=0;	// number of pseudo-tracks, i.e. observations at lag 0
		private int noOfMaxLag  =0;	// observations at maximum lag
		private int noOfMinLag  =0;	// observations at minimum lag
		
		private SingleParticleStatResult lsr=null;
		
		
		/**
		 * constructor
		 */
		BinStatistics(int tRad){
			lsr=new SingleParticleStatResult(tRad);
			
			Particle p1=null;
			
			for(Particle p:ls) if(p.getTCount()>1) p1=p;
			
			if(p1==null) throw new IllegalArgumentException("no valid particle data");
			
			lsr.dt=new MDate(p1.getTime(1)).getDT(new MDate(p1.getTime(0)));
			
			//writeTracks("/lustre/home/qianyk/Data/GDP/Track/diff9regions/"+lon1+lat1);
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
							
							float disX=SpatialModel.EARTH_RADIUS*dlon*(float)Math.cos(oY*Math.PI/180.0);
							float disY=SpatialModel.EARTH_RADIUS*dlat;
							
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
				float[] accx=hasAcc?p.getAttachedData(2):null;
				float[] accy=hasAcc?p.getAttachedData(3):null;
				
				for(int l=0,L=p.getTCount();l<L;l++){
					// meet condition to be an origin of pseudo-track
					if(!cond.test(p.getRecord(l))) continue;
					if(uspd[l]==Record.undef) continue;
					
					float ua0=uspd[l];//-um[tRad];	// u'(tau=0)
					float va0=vspd[l];//-vm[tRad];	// v'(tau=0)
					
					for(int ll=0;ll<L;ll++){
						int tau=ll-l;
						
						if(tau>tRad||tau<-tRad) continue;
						if(uspd[ll]==Record.undef) continue;
						
						float ua=uspd[ll];	// u'(tau)
						float va=vspd[ll];	// v'(tau)
						
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
							float disX=SpatialModel.EARTH_RADIUS*dlon*(float)Math.cos(oY*Math.PI/180.0)-lsr.DXm[tRad+tau];
							// dy'(tau)
							float disY=SpatialModel.EARTH_RADIUS*dlat-lsr.DYm[tRad+tau];
							
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
					
					float ua0=uspd[l];//-um[tRad];	// u'(tau=0)
					float va0=vspd[l];//-vm[tRad];	// v'(tau=0)
					
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
							float disX=SpatialModel.EARTH_RADIUS*dlon*(float)Math.cos(oY*Math.PI/180.0)-lsr.DXm[tRad+tau];
							// dy'(tau)
							float disY=SpatialModel.EARTH_RADIUS*dlat-lsr.DYm[tRad+tau];
							
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
			
			computeAutoCovariance(cond);
			
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
		
		// compute diffusivity by differencing dispersion
		void computeDiffByDispDT(Predicate<Record> cond){
			int tRad=lsr.tRad,len=2*tRad+1;
			
			float ucfd=0,vcfd=0;
			
			computeDispersion(cond);
			
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
				double KxyS=(lsr.Kxy[l]+lsr.Kyx[l])/2f;
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
