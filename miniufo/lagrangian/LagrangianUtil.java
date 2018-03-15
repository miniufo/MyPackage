/**
 * @(#)LagrangianUtil.java	1.0 2014.08.07
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Stream;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.SpatialModel;


/**
 * utility for Lagrangian data
 *
 * @version 1.0, 2014.08.07
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LagrangianUtil{
	
	/**
	 * count all the records in the list
	 * 
	 * @param	ls	a give list
	 */
	public static int cTotalCount(List<? extends Particle> ls){
		int count=0;
		
		for(Particle p:ls) count+=p.getTCount();
		
		return count;
	}
	
	/**
	 * count the records that have velocity data
	 * 
	 * @param	ls	a give list
	 */
	public static int cDefinedCount(List<? extends Particle> ls){
		int count=0;
		
		for(Particle p:ls) for(int l=0,L=p.getTCount();l<L;l++)
		if(p.getRecord(l).getDataValue(Particle.UVEL)!=Record.undef) count++;
		
		return count;
	}
	
	/**
	 * compute all the record count in unit of drifter-year
	 * 
	 * @param	ls	a give list
	 */
	public static float cTotalDrifterYear(List<? extends Particle> ls){
		int dt=-1;
		
		for(Particle p:ls)
		if(p.getTCount()>1){
			dt=new MDate(p.getTime(1)).getDT(new MDate(p.getTime(0)));
			break;
		}
		
		if(dt==-1) throw new IllegalArgumentException("cannot determine delta-T");
		
		float ratio=dt/(3600f*24f*365f);
		
		return cTotalCount(ls)*ratio;
	}
	
	/**
	 * compute mean record count per drifter in unit of drifter-year
	 * 
	 * @param	ls	a give list
	 */
	public static float cMeanDrifterYear(List<? extends Particle> ls){
		return cTotalDrifterYear(ls)/ls.size();
	}
	
	/**
	 * compute all the record length in unit of km
	 * 
	 * @param	ls	a give list
	 */
	public static float cTotalTrajLength(List<? extends Particle> ls){
		float sum=0;
		
		for(Particle p:ls){
			float dis=0;
			
			float[] lons=p.getXPositions();
			float[] lats=p.getYPositions();
			
			for(int l=0,L=p.getTCount()-1;l<L;l++)
			dis+=SpatialModel.cSphericalDistanceByDegree(
				lons[l],lats[l],lons[l+1],lats[l+1]
			);
			
			sum+=dis;
		}
		
		return sum/1000;
	}
	
	/**
	 * compute mean the record length in unit of km
	 * 
	 * @param	ls	a give list
	 */
	public static float cMeanTrajLength(List<? extends Particle> ls){
		return cTotalTrajLength(ls)/ls.size();
	}
	
	
	/**
	 * change to a record stream
	 * 
	 * @param	ls	a list of particles
	 */
	public static Stream<Record> asRecordStream(List<? extends Particle> ls){
		Function<Particle,Stream<Record>> f=p->p.records.stream();
		
		return ls.stream().flatMap(f);
	}
	
	
	/**
	 * Output trajectories if the start point is in a given region
	 * 
	 * @param	ls		a list of particles
	 * @param	path	folder for output
	 * @param	wrtAtt	write attached data or not
	 * @param	cond	write to file if a particle meets this condition
	 */
	public static void writeTrajecories(List<? extends Particle> ls,String path,boolean wrtAtt,Predicate<Particle> cond){
		for(Particle p:ls) if(cond.test(p)) p.toTrajectoryFile(path,wrtAtt);
	}
	
	public static void writeTrajecories(List<? extends Particle> ls,String path){
		writeTrajecories(ls,path,false,p->true);
	}
	
	
	/**
	 * remove particles if a given condition is met
	 * 
	 * @param	ls		a give list
	 * @param	cond	a give condition based on particle
	 */
	public static void removeByCond(List<? extends Particle> ls,Predicate<Particle> cond){
		int total=ls.size();
		
		List<Particle> remove=new ArrayList<Particle>();
		
		for(Particle p:ls) if(cond.test(p)) remove.add(p);
		
		ls.removeAll(remove);
		
		int remain=ls.size();
		
		System.out.println("remain/total count: "+remain+"/"+total+" ("+(remain*100f/total)+"%)");
	}
	
	public static void removeLessThanT(List<? extends Particle> ls,int tcount){
		removeByCond(ls,p->p.getTCount()<tcount);
	}
	
	public static void removeNotEqualT(List<? extends Particle> ls,int tcount){
		removeByCond(ls,p->p.getTCount()!=tcount);
	}
	
	public static void removeLargeDeltaLon(List<? extends Particle> ls,float dlon){
		removeByCond(ls,p->Math.abs(p.getRecord(0).getXPos()-p.getRecord(p.getTCount()-1).getXPos())>dlon);
	}
	
	public static void removeLargeDeltaLat(List<? extends Particle> ls,float dlat){
		removeByCond(ls,p->Math.abs(p.getRecord(0).getYPos()-p.getRecord(p.getTCount()-1).getYPos())>dlat);
	}
	
	
	/**
	 * Write a list of Particle to binary file
	 * 
	 * @param	ls		a give list of Particle
	 * @param	fname	file name
	 */
	public static void writeAsBinaryFile(List<Particle> ls,String fname){
		try(ObjectOutputStream oos=new ObjectOutputStream(new FileOutputStream(fname))){
			oos.writeObject(ls);
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	/**
	 * read a list of Particle from binary file
	 * 
	 * @param	fname	file name
	 */
	@SuppressWarnings("unchecked")
	public static List<Particle> readFromBinaryFile(String fname){
		try(ObjectInputStream oos=new ObjectInputStream(new FileInputStream(fname))){
			return (List<Particle>)oos.readObject();
			
		}catch(Exception e){ e.printStackTrace(); System.exit(0);}
		
		throw new IllegalArgumentException("should not reach here");
	}
	
	
	/** test
	public static void main(String[] args){
		String path="d:/Data/dd.TXT/abc.cc.dd";
		
		System.out.println(getFileName(path));
		System.out.println(getFilePrefixName(path));
		System.out.println(getFileSuffixName(path));
		System.out.println(getCompleteFileNameWithoutExtension(path));
		System.out.println(getCompletePath(path));
	}*/
}
