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
import java.util.stream.Stream;

import miniufo.diagnosis.MDate;
import miniufo.diagnosis.SpatialModel;

import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.toRadians;


/**
 * Used to describe a Lagrangian particle in a fluid.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class Particle implements Cloneable,Serializable{
	//
	private static final long serialVersionUID = -5698371935283702966L;

	protected boolean finished=false;
	protected boolean llpos   =true;
	
	protected String id=null;
	
	protected AttachedMeta[] meta=null;	// names and indices for attached data
	
	protected List<Record> records=null;
	
	public static final AttachedMeta UVEL=new AttachedMeta("uvel",0);
	public static final AttachedMeta VVEL=new AttachedMeta("vvel",1);
	
	
	/**
	 * constructor
	 * 
	 * @param	r	initial record
	 */
	public Particle(String id,int attLen){
		this(id,10,attLen,true);
	}
	
	public Particle(String id,int initSize,int attLen,boolean llpos){
		this.id=id;
		
		meta=new AttachedMeta[attLen];
		
		records=new ArrayList<>(initSize);
		
		this.llpos=llpos;
	}
	
	public Particle(String id,Record init,int attLen){
		this(id,10,attLen,true);
		
		records.add(init);
	}
	
	
	/*** getor and setor ***/
	public int getTCount(){ return records.size();}
	
	public int getMedianIndex(){ return (int)((records.size()+0.5f)/2);}
	
	public int getDataIndex(String name){
		for(int i=0,I=meta.length;i<I;i++)
		if(name.equalsIgnoreCase(meta[i].name)) return meta[i].index;
		
		throw new IllegalArgumentException("cannot find attachedMeta: "+name);
	}
	
	public boolean isFinished(){ return finished;}
	
	public boolean isLatLonPosition(){ return llpos;}
	
	public boolean addRecord(Record r){ return records.add(r);}
	
	public boolean addAll(List<Record> ls){ return records.addAll(ls);}
	
	public boolean remove(Record r){ return records.remove(r);}
	
	public boolean removeAll(List<Record> al){ return records.removeAll(al);}
	
	public long getTime(int l){ return records.get(l).getTime();}
	
	public float getXPosition(int i){ return records.get(i).getXPos();}
	
	public float getYPosition(int i){ return records.get(i).getYPos();}
	
	public long[] getTimes(){
		int size=records.size();
		
		long[] times=new long[size];
		
		for(int l=0;l<size;l++) times[l]=records.get(l).getTime();
		
		return times;
	}
	
	public float[] getInitialLocation(){
		Record init=records.get(0);
		
		return new float[]{init.getXPos(),init.getYPos()};
	}
	
	public float[] getXPositions(){
		int size=records.size();
		
		float[] lons=new float[size];
		
		for(int i=0;i<size;i++) lons[i]=records.get(i).getXPos();
		
		return lons;
	}
	
	public float[] getYPositions(){
		int size=records.size();
		
		float[] lats=new float[size];
		
		for(int i=0;i<size;i++) lats[i]=records.get(i).getYPos();
		
		return lats;
	}
	
	public float[] getUVel(){ return getAttachedData(UVEL);}
	
	public float[] getVVel(){ return getAttachedData(VVEL);}
	
	public float[] getSpeeds(){
		int len=records.size();
		
		float[] xpos=getXPositions();
		float[] ypos=getYPositions();
		
		if(len>1){
			float dt=getDT(records.get(0).getTime(),records.get(1).getTime());
			
			float[] spds=new float[len];
			
			if(llpos){
				spds[0    ]=SpatialModel.cSphericalDistanceByDegree(xpos[0    ],ypos[0    ],xpos[1    ],ypos[1    ])/dt;
				for(int l=1,L=len-1;l<L;l++)
				spds[l    ]=SpatialModel.cSphericalDistanceByDegree(xpos[l-1  ],ypos[l-1  ],xpos[l+1  ],ypos[l+1  ])/(2f*dt);
				spds[len-1]=SpatialModel.cSphericalDistanceByDegree(xpos[len-2],ypos[len-2],xpos[len-1],ypos[len-1])/dt;
				
				return spds;
				
			}else{
				spds[0    ]=(float)(Math.hypot(xpos[0    ]-xpos[1    ],ypos[0    ]-ypos[1    ])/dt);
				for(int l=1,L=len-1;l<L;l++)
				spds[l    ]=(float)(Math.hypot(xpos[l-1  ]-xpos[l+1  ],ypos[l-1  ]-ypos[l+1  ])/(2f*dt));
				spds[len-1]=(float)(Math.hypot(xpos[len-2]-xpos[len-1],ypos[len-2]-ypos[len-1])/dt);
				
				return spds;
			}
			
		}else return new float[len];
	}
	
	public float[] getAttachedData(AttachedMeta meta){
		int size=records.size();
		
		float[] dr=new float[size];
		
		for(int i=0;i<size;i++) dr[i]=records.get(i).getDataValues()[meta.index];
		
		return dr;
	}
	
	public String getID(){ return id;}
	
	public Record remove(int i){ return records.remove(i);}
	
	public Record getRecord(int idx){ return records.get(idx);}
	
	public AttachedMeta[] getAttachedMeta(){ return meta;}
	
	public void setAttachedMeta(AttachedMeta... meta){
		for(int i=0,I=meta.length;i<I;i++) this.meta[i]=meta[i];
	}
	
	public Record setRecord(int idx,Record r){ return records.set(idx,r);}
	
	public void finish(){ finished=true;}
	
	public Stream<Record> stream(){ return records.stream();}
	
	
	/**
	 * sort the profiles according to the timestamp
	 */
	public void sort(){ records.sort((r1,r2)->Long.compare(r1.getTime(),r2.getTime()));}
	
	
	/**
	 * Write the trajectory info to a ASCII file.
	 * 
	 * @param	path	path of the file, file name is this Particle's ID
	 * @param	wrtAtt	write attached data or not
	 */
	public void toTrajectoryFile(String path,boolean wrtAtt){
		StringBuilder buf=new StringBuilder();
		
		buf.append("* "+records.size()+" id: "+id+"\n");
		
		llpos=isLatLonPosition();
		
		for(Record r:records){
			buf.append(r.toString(llpos,wrtAtt));
			buf.append("\n");
		}
		
		try(FileWriter fw=new FileWriter(new File(path+id+".txt"))){
			fw.write(buf.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public void toTrajectoryFile(String path){ toTrajectoryFile(path,false);}
	
	
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
			
			float dX= next.getXPos()-prev.getXPos();
			float dY= next.getYPos()-prev.getYPos();
			float mY=(next.getYPos()+prev.getYPos())/2f;
			
			llpos=isLatLonPosition();
			
			if(llpos){
				pres.setData(UVEL,(float)(dX/180.0*PI*EARTH_RADIUS/dt*cos(toRadians(mY))));
				pres.setData(VVEL,(float)(dY/180.0*PI*EARTH_RADIUS/dt));
				
			}else{
				pres.setData(UVEL,dX/dt);
				pres.setData(VVEL,dY/dt);
			}
			
			// for all records between
			for(int l=1,L=records.size()-1;l<L;l++){
				prev=records.get(l-1);
				pres=records.get(l  );
				next=records.get(l+1);
				
				dX= next.getXPos()-prev.getXPos();
				dY= next.getYPos()-prev.getYPos();
				mY=(next.getYPos()+prev.getYPos())/2f;
				
				if(llpos){
					pres.setData(UVEL,(float)(dX/180.0*PI*EARTH_RADIUS/(dt*2.0)*cos(toRadians(mY))));
					pres.setData(VVEL,(float)(dY/180.0*PI*EARTH_RADIUS/(dt*2.0)));
					
				}else{
					pres.setData(UVEL,dX/dt/2f);
					pres.setData(VVEL,dY/dt/2f);
				}
			}
			
			// for the last record
			prev=records.get(records.size()-2);
			next=records.get(records.size()-1);
			pres=next;
			
			dX= next.getXPos()-prev.getXPos();
			dY= next.getYPos()-prev.getYPos();
			mY=(next.getYPos()+prev.getYPos())/2f;
			
			if(llpos){
				pres.setData(UVEL,(float)(dX/180.0*PI*EARTH_RADIUS/dt*cos(toRadians(mY))));
				pres.setData(VVEL,(float)(dY/180.0*PI*EARTH_RADIUS/dt));
				
			}else{
				pres.setData(UVEL,dX/dt);
				pres.setData(VVEL,dY/dt);
			}
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
	 * Get a view of portion of this Particle records.
	 * 
	 * @param	str		start index (from 0)
	 * @param	len		length of records
	 */
	public void subRecords(int str,int len){ records=records.subList(str,str+len);}
	
	
	/**
     * remove records if they meet some requirement
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
			p.llpos=llpos;
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
		
		if(llpos)
			buf.append("  lons(E-deg)  lats(N-deg) time(YYYYMMDDHHMMSS) uspd(m/s) vspd(m/s)\n");
		else
			buf.append(" xpos(km)  ypos(km)   time(YYYYMMDDHHMMSS) uspd(m/s) vspd(m/s)\n");
		
		for(Record r:records){
			buf.append(r.toString(llpos,false));
			buf.append("\n");
		}
		
		return buf.toString();
	}
	
	
	/*** helper methods ***/
	private static int getDT(long t1,long t2){
		return new MDate(t1).getDT(new MDate(t2));
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
