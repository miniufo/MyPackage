/**
 * @(#)ParticlePair.java	1.0 2017.08.28
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import miniufo.lagrangian.Particle;
import miniufo.test.diagnosis.MDate;


/**
 * Used to describe particle pair for computing relative dispersion.
 *
 * @version 1.0, 2017.08.28
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ParticlePair{
	//
	private int deltaT =0;	// in unit of second
	private int tcount =0;	// time count
	
	private String   id=null;
	
	private Particle p1=null;
	private Particle p2=null;
	
	
	/**
     * Constructor.
     *
     * @param	p1	the  first particle
     * @param	p2	the second particle
     */
	public ParticlePair(String id,Particle p1,Particle p2){
		tcount=p1.getTCount();
		
		if(p2.getTCount()!=tcount) throw new IllegalArgumentException("two particle should be of the same t-length");
		
		this.id    =id;
		this.p1    =p1;
		this.p2    =p2;
		this.deltaT=getDeltaTimeInSecond(p1);
		
		int deltaT2=getDeltaTimeInSecond(p2);
		
		if(deltaT!=deltaT2) throw new IllegalArgumentException(
			"delta-t for the two particles ("+deltaT+"!="+deltaT2+") are not the same\n"+p1
		);
	}
	
	
	/**
	 * Compute x-distance.
	 * 
	 * @param	ttag	t-index at which distance is computed
	 */
	public float cXDistance(int ttag){
		return p1.getRecord(ttag).cXDistanceTo(p2.getRecord(ttag),p1.isLatLonPosition());
	}
	
	/**
	 * Compute y-distance.
	 * 
	 * @param	ttag	t-index at which distance is computed
	 */
	public float cYDistance(int ttag){
		return p1.getRecord(ttag).cYDistanceTo(p2.getRecord(ttag),p1.isLatLonPosition());
	}
	
	/**
	 * Compute total distance.
	 * 
	 * @param	ttag	t-index at which distance is computed
	 */
	public float cDistance(int ttag){
		return p1.getRecord(ttag).cDistanceTo(p2.getRecord(ttag),p1.isLatLonPosition());
	}
	
	
	/**
	 * Compute x-distances for all time steps.
	 */
	public float[] cXDistances(){
		int t=p1.getTCount();
		
		float[] seps=new float[t];
		
		for(int l=0;l<t;l++) seps[l]=cXDistance(l);
		
		return seps;
	}
	
	/**
	 * Compute y-distances for all time steps.
	 */
	public float[] cYDistances(){
		int t=p1.getTCount();
		
		float[] seps=new float[t];
		
		for(int l=0;l<t;l++) seps[l]=cYDistance(l);
		
		return seps;
	}
	
	/**
	 * Compute total-distances for all time steps.
	 */
	public float[] cDistances(){
		int t=p1.getTCount();
		
		float[] seps=new float[t];
		
		for(int l=0;l<t;l++) seps[l]=cDistance(l);
		
		return seps;
	}
	
	
	/*** getor and setor ***/
	public int getDeltaT(){ return deltaT;}
	
	public int getTCount(){ return tcount;} 
	
	public String getID(){ return id;}
	
	public Particle getParticle1(){ return p1;}
	
	public Particle getParticle2(){ return p2;}
	
	
	/*** helper methods ***/
	private int getDeltaTimeInSecond(Particle p){
		if(p.getTCount()<2) throw new IllegalArgumentException("cannot get delta time as there are only one record");
		return new MDate(p.getTime(0)).getDT(new MDate(p.getTime(1)));
	}
	
	
	/**
	 * used to print out
	 */
	public String toString(){
		return String.format("particle pair %8s consists of P1 %8s and P2 %8s, initial separation: %9.5f km",id,p1.getID(),p2.getID(),cDistance(0)/1000f);
	}
}
