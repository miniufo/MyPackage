/**
 * @(#)TwoParticleStatistics.java	1.0 2017.08.31
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.Predicate;
import miniufo.descriptor.DataDescriptor;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;


/**
 * Two-particle statistics.
 *
 * @version 1.0, 2017.08.31
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TwoParticleStatistics extends ParticleStatistics{
	
	/**
	 * prevent from initialization
	 */
	public TwoParticleStatistics(List<? extends Particle> ls,DataDescriptor dd){ super(ls,dd);}
	
	
	public List<ParticlePair> selectPairs(Predicate<ParticlePair> cond){
		int pcount=ls.size();
		int id=0;
		
		List<ParticlePair> res=new ArrayList<>();
		
		for(int p=0;p<pcount;p++)
		for(int pp=p+1;pp<pcount;pp++){
			ParticlePair pair=new ParticlePair(String.valueOf(id++),ls.get(p),ls.get(pp));
			
			if(cond.test(pair)) res.add(pair);
		}
		
		return res;
	}
	
	
	/**
	 * Compute relative dispersion.
	 * 
	 * @param	pairs	a list of particle pairs
	 */
	public TwoParticleStatResult cRelativeDispersion(List<ParticlePair> pairs){
		int tlen=pairs.get(0).getParticle1().getTCount();
		
		TwoParticleStatResult mpsr=new TwoParticleStatResult(tlen);
		
		mpsr.dt=pairs.get(0).getDeltaT();
		mpsr.numOfPairs=pairs.size();
		
		int  [] num=mpsr.num;
		float[] Dxx=mpsr.Dxx;
		float[] Dyy=mpsr.Dyy;
		float[] Dis=mpsr.Dis;
		float[] Kxx=mpsr.Kxx;
		float[] Kyy=mpsr.Kyy;
		
		for(ParticlePair pp:pairs)
		for(int l=0,L=pp.getTCount();l<L;l++){
			float dx=pp.cXDistance(l);
			float dy=pp.cYDistance(l);
			float ds=pp.cDistance(l);
			
			if(dx==Record.undef||dy==Record.undef||ds==Record.undef) continue;
			
			num[l]++;
			
			Dxx[l]+=dx;
			Dyy[l]+=dy;
			Dis[l]+=ds;
		}
		
		Kxx[0]=(Dxx[1]-Dxx[0])/mpsr.dt/2f;
		Kyy[0]=(Dyy[1]-Dyy[0])/mpsr.dt/2f;
		
		for(int l=1;l<tlen-1;l++){
			Kxx[l]=(Dxx[l+1]-Dxx[l-1])/mpsr.dt/4f;
			Kyy[l]=(Dyy[l+1]-Dyy[l-1])/mpsr.dt/4f;
		}
		
		Kxx[tlen-1]=(Dxx[tlen-1]-Dxx[tlen-2])/mpsr.dt/2f;
		Kyy[tlen-1]=(Dyy[tlen-1]-Dyy[tlen-2])/mpsr.dt/2f;
		
		return mpsr;
	}
	
	
	/**
	 * Re-sample the pairs using bootstrapping.
	 * 
	 * @param	sampleNum	how many pairs are needed to re-sample
	 * @param	pp			a list of ParticlePair
	 */
	public static List<ParticlePair> bootstrapping(int sampleNum,List<ParticlePair> pp){
		List<ParticlePair> spl=new ArrayList<>(sampleNum);
		
		Random rnd=new Random();
		
		for(int i=0;i<sampleNum;i++){
			// uniformly distributed between [0,sampleNum)
			int idx=rnd.nextInt(pp.size());
			
			spl.add(pp.get(idx));
		}
		
		return spl;
	}
	
	
	/*** getor and setor ***/
	
	
	/*** helper methods ***/
	
	
	/** test
	public static void main(String arg[]){
		
	}*/
}
