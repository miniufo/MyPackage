/**
 * @(#)Particle.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.List;
import java.util.ArrayList;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.SpatialModel;
import static java.lang.Math.cos;
import static java.lang.Math.toRadians;


/**
 * used to describe a particle model in fluid
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class Particle implements Cloneable,Serializable{
	//
	private static final long serialVersionUID = -5698371935283702966L;

	protected boolean finished=false;
	
	protected String id=null;
	
	protected String[] attachedVars=null;	// names for attached data
	
	protected List<Record> records=null;
	
	
	/**
	 * constructor
	 * 
	 * @param	r	initial record
	 */
	public Particle(String id,int attLen){
		this(id,10,attLen);
	}
	
	public Particle(String id,int initSize,int attLen){
		this.id=id;
		
		attachedVars=new String[attLen];
		
		records=new ArrayList<>(initSize);
	}
	
	public Particle(String id,Record init,int attLen){
		this(id,10,attLen);
		
		records.add(init);
	}
	
	
	/*** getor and setor ***/
	public int getTCount(){ return records.size();}
	
	public int getMedianIndex(){ return (int)((records.size()+0.5f)/2);}
	
	public int getDataIndex(String name){
		for(int i=0,I=attachedVars.length;i<I;i++)
		if(name.equalsIgnoreCase(attachedVars[i])) return i;
		
		throw new IllegalArgumentException("cannot find attached data: "+name);
	}
	
	public boolean isFinished(){ return finished;}
	
	public boolean addRecord(Record r){ return records.add(r);}
	
	public boolean addAll(List<Record> ls){ return records.addAll(ls);}
	
	public boolean remove(Record r){ return records.remove(r);}
	
	public boolean removeAll(List<Record> al){ return records.removeAll(al);}
	
	public long getTime(int l){ return records.get(l).getTime();}
	
	public float getLongitude(int i){ return records.get(i).getLon();}
	
	public float getLatitude(int i){ return records.get(i).getLat();}
	
	public long[] getTimes(){
		int size=records.size();
		
		long[] times=new long[size];
		
		for(int l=0;l<size;l++) times[l]=records.get(l).getTime();
		
		return times;
	}
	
	public float[] getInitialLocation(){
		Record median=records.get(0);
		
		return new float[]{median.getLon(),median.getLat()};
	}
	
	public float[] getLongitudes(){
		int size=records.size();
		
		float[] lons=new float[size];
		
		for(int i=0;i<size;i++) lons[i]=records.get(i).getLon();
		
		return lons;
	}
	
	public float[] getLatitudes(){
		int size=records.size();
		
		float[] lats=new float[size];
		
		for(int i=0;i<size;i++) lats[i]=records.get(i).getLat();
		
		return lats;
	}
	
	public float[] getZonalVelocity(){ return getAttachedData(0);}
	
	public float[] getMeridionalVelocity(){ return getAttachedData(1);}
	
	public float[] getAttachedData(int idx){
		int size=records.size();
		
		float[] dr=new float[size];
		
		for(int i=0;i<size;i++) dr[i]=records.get(i).getDataValues()[idx];
		
		return dr;
	}
	
	public float[] getAttachedData(String name){ return getAttachedData(getDataIndex(name));}
	
	public String getID(){ return id;}
	
	public String[] getDataNames(){ return attachedVars;}
	
	public Record remove(int i){ return records.remove(i);}
	
	public Record getRecord(int idx){ return records.get(idx);}
	
	public void setAttachedDataNames(String... names){
		for(int i=0,I=names.length;i<I;i++) attachedVars[i]=names[i];
	}
	
	public Record setRecord(int idx,Record r){ return records.set(idx,r);}
	
	public void finish(){ finished=true;}
	
	
	/**
	 * sort the profiles according to the timestamp
	 */
	public void sort(){ records.sort((r1,r2)->Long.compare(r1.getTime(),r2.getTime()));}
	
	
	/**
	 * write the trajectory info to a file
	 * so that it could be used to draw figures
	 * 
	 * @param	path	path of the file, file name is this Particle's ID
	 */
	public void toTrajectoryFile(String path){
		StringBuilder buf=new StringBuilder();
		
		buf.append("* "+records.size()+" id: "+id+"\n");
		
		for(Record r:records){
			buf.append(r);
			buf.append("\n");
		}
		
		try(FileWriter fw=new FileWriter(new File(path+id+".txt"))){
			fw.write(buf.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/**
	 * compute velocity using centered difference of position
	 */
	public void cVelocityByPosition(){
		if(records.size()>1){
			float dt=new MDate(records.get(1).getTime()).getDT(new MDate(records.get(0).getTime()));
			
			// for the first record
			Record prev=records.get(0);
			Record next=records.get(1);
			Record pres=prev;
			
			float dlon=next.getLon()-prev.getLon();
			float dlat=next.getLat()-prev.getLat();
			float clat=(next.getLat()+prev.getLat())/2;
			
			float u=(float)(dlon/180.0*Math.PI*SpatialModel.EARTH_RADIUS/dt*cos(toRadians(clat)));
			float v=(float)(dlat/180.0*Math.PI*SpatialModel.EARTH_RADIUS/dt);
			
			pres.setData(0,u);
			pres.setData(1,v);
			
			// for all records between
			for(int l=1,L=records.size()-1;l<L;l++){
				prev=records.get(l-1);
				pres=records.get(l  );
				next=records.get(l+1);
				
				dlon=next.getLon()-prev.getLon();
				dlat=next.getLat()-prev.getLat();
				clat=(next.getLat()+prev.getLat())/2;
				
				u=(float)(dlon/180.0*Math.PI*SpatialModel.EARTH_RADIUS/(dt*2.0)*cos(toRadians(clat)));
				v=(float)(dlat/180.0*Math.PI*SpatialModel.EARTH_RADIUS/(dt*2.0));
				
				pres.setData(0,u);
				pres.setData(1,v);
			}
			
			// for the last record
			prev=records.get(records.size()-2);
			next=records.get(records.size()-1);
			pres=next;
			
			dlon=next.getLon()-prev.getLon();
			dlat=next.getLat()-prev.getLat();
			clat=(next.getLat()+prev.getLat())/2;
			
			u=(float)(dlon/180.0*Math.PI*SpatialModel.EARTH_RADIUS/dt*cos(toRadians(clat)));
			v=(float)(dlat/180.0*Math.PI*SpatialModel.EARTH_RADIUS/dt);
			
			pres.setData(0,u);
			pres.setData(1,v);
		}
	}
	
	
	/**
	 * split the float data into segments of a given size
	 */
	public Particle[] split(int size){
		if(size<1) throw new IllegalArgumentException("size should be positive");
		
		int len=records.size()/size;
		int lef=records.size()%size;
		
		if(lef!=0) len++;
		
		Particle[] re=new Particle[len];
		
		int ptr=0;
		
		for(int i=0,I=(lef!=0?len-1:len);i<I;i++){
			re[i]=new Particle(id,size);
			
			for(int j=0;j<size;j++)
			re[i].addRecord(records.get(ptr++));
		}
		
		if(lef!=0){
			re[len-1]=new Particle(id,lef);
			
			for(int j=0,J=lef;j<J;j++)
			re[len-1].addRecord(records.get(ptr++));
		}
		
		return re;
	}
	
	
	/**
	 * get a subList of records into a new Particle
	 * sharing the same memory
	 */
	public Particle subRecord(String newID,int str,int len){
		Particle p=new Particle(newID,len);
		
		for(int l=str,L=str+len;l<L;l++) p.records.add(records.get(l));
		
		//p.records=records.subList(str,str+len);
		
		return p;
	}
	
	public Particle subRecord(int str,int len){ return subRecord(id,str,len);}
	
	
	/**
     * delete records that not meet some requirement
     * 
     * @param	cond	an condition whether a record should be removed
     */
	public Particle removeRecords(Predicate<Record> cond){
		List<Record> remove=records.stream().filter(cond).collect(Collectors.toList());
		
		records.removeAll(remove);
		
		return this;
	}
	
	public Particle removeRecords(boolean[] del){
		int len=del.length;
		
		if(len!=records.size())
		throw new IllegalArgumentException("length of del array should be equal to count");
		
		List<Record> remove=new ArrayList<>();
		
		for(int l=0;l<len;l++) if(del[l]) remove.add(records.get(l));
		
		records.removeAll(remove);
		
		return this;
	}
	
	
	/**
     * clone method
     */
	public Object clone(){
		try{
		    Particle p=(Particle)super.clone();
		    
			p.finished=finished;
			p.id=id;
			p.records=new ArrayList<>(records.size());
			
			for(int i=0,I=records.size();i<I;i++)
		    p.addRecord(new Record(records.get(i)));
		    
			return p;
			
	    }catch(CloneNotSupportedException ex){
		    // this shouldn't happen, since we are Cloneable
		    throw new InternalError();
	    }
	}
	
	
	/**
     * used to print out
     */
	public String toString(){
		StringBuilder buf=new StringBuilder();
		
		buf.append(getClass().getSimpleName()+" id ("+id+") "+records.size()+" records:\n");
		buf.append("  lons(E-deg)  lats(N-deg) uspd(m/s) vspd(m/s) time(YYYYMMDDHHMMSS)\n");
		
		for(Record r:records){
			buf.append(r);
			buf.append("\n");
		}
		
		return buf.toString();
	}
	
	
	/** test
	public static void main(String[] arg){
		Particle p=new Particle(1,10);
		
		for(int l=0;l<10;l++)
		p.addRecord(new Record(l));
		
		System.out.println(p);
		System.out.println();
		System.out.println(p.subParticle(0,10));
	}*/
}
