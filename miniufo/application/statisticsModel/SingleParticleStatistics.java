/**
 * @(#)SingleParticleStatistics.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import miniufo.descriptor.DataDescriptor;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;


/**
 * Single particle statistics.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class SingleParticleStatistics extends ParticleStatistics{
	
	/**
	 * prevent from initialization
	 */
	public SingleParticleStatistics(List<? extends Particle> ls,DataDescriptor dd){ super(ls,dd);}
	
	
	/**
	 * write distribution data
	 */
	public void writeMedianDistribution(String path){
		int len=ls.size();
		
		float[][] re=new float[2][len];
		
		for(int i=0;i<len;i++){
			Particle p=ls.get(i);
			
			Record r=p.getRecord(p.getMedianIndex());
			
			re[0][i]=r.getXPos();
			re[1][i]=r.getYPos();
		}
		
		writeDistribution(re[0],re[1],path);
	}
	
	public void writeInitialDistribution(String path){
		int len=ls.size();
		
		float[][] re=new float[2][len];
		
		for(int i=0;i<len;i++){
			float[] pos=ls.get(i).getInitialLocation();
			
			re[0][i]=pos[0];
			re[1][i]=pos[1];
		}
		
		writeDistribution(re[0],re[1],path);
	}
	
	
	/*** helper methods ***/
	private static void writeDistribution(float[] xpos,float[] ypos,String path){
		StringBuilder sb=new StringBuilder();
		
		for(int l=0,L=xpos.length;l<L;l++) sb.append(xpos[l]+" "+ypos[l]+"\n");
		
		try(FileWriter fw=new FileWriter(path)){
			fw.write(sb.toString());
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String arg[]){
		
	}*/
}
