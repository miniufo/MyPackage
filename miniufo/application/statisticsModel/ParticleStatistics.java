/**
 * @(#)ParticleStatistics.java	1.0 2017.08.31
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.util.List;
import miniufo.descriptor.DataDescriptor;
import miniufo.lagrangian.Particle;


/**
 * Particle statistics.
 *
 * @version 1.0, 2017.08.31
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class ParticleStatistics{
	//
	protected boolean llpos=true;
	
	protected List<? extends Particle> ls=null;
	
	protected DataDescriptor dd=null;
	
	
	/**
	 * prevent from initialization
	 */
	public ParticleStatistics(List<? extends Particle> ls,DataDescriptor dd){
		this.ls=ls;
		this.dd=dd;
		this.llpos=ls.get(0).isLatLonPosition();
	}
	
	
	/*** getor and setor ***/
	public DataDescriptor getDataDescriptor(){ return dd;}
	
	
	/** test
	public static void main(String arg[]){
		
	}*/
}
