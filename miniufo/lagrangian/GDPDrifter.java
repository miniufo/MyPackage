/**
 * @(#)GDPDrifter.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.ArrayList;
import java.util.List;
import miniufo.diagnosis.MDate;
import static miniufo.lagrangian.Record.undef;


/**
 * Global Drifter Program drifters
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GDPDrifter extends Particle{
	//
	private static final long serialVersionUID = -8586197613070975714L;
	
	private TYPE type=TYPE.SVP;
	
	public static enum TYPE{SVP,SVPB,SVPBS,SVPBW,SVPC,SVPI,SVPW}
	
	public static final AttachedMeta Temp  =new AttachedMeta("temp"   ,2);
	public static final AttachedMeta DrgOff=new AttachedMeta("drogOff",3);
	
	
	/**
	 * constructor
	 */
	public GDPDrifter(String id,int attLen){ super(id,attLen);}
	
	public GDPDrifter(String id,int size,int attLen){ super(id,size,attLen,true);}
	
	
	/*** getor and setor ***/
	public TYPE getType(){ return type;}
	
	
	/**
	 * remove undef velocity record at the start and the end of the records
	 */
	public void removeEndpointUndefVelocityRecords(){
		List<Record> removes=new ArrayList<Record>(2);
		
		Record r=null;
		float[] ad=null;
		
		r=records.get(0);
		ad=r.getDataValues();
		if(ad[0]==undef||ad[1]==undef) removes.add(r);
		
		r=records.get(records.size()-1);
		ad=r.getDataValues();
		if(ad[0]==undef||ad[1]==undef) removes.add(r);
		
		records.removeAll(removes);
	}
	
	
	/**
	 * whether the data has undef records
	 * 
	 * @param	idx		the index of AttachedData to be checked
	 */
	public boolean hasUndefRecords(AttachedMeta meta){
		for(Record r:records)
		if(r.getDataValue(meta)==undef) return true;
		
		return false;
	}
	
	/**
	 * whether the data has undef velocity records
	 */
	public boolean hasLargeVelocityRecords(float thre){
		for(Record r:records){
			float u=r.getDataValue(UVEL);
			float v=r.getDataValue(VVEL);
			
			if((u!=undef&&Math.abs(u)>thre)||(v!=undef&&Math.abs(v)>thre)) return true;
		}
		
		return false;
	}
	
	/**
	 * whether the data are continuous with 6-hour interval
	 */
	public boolean isContinuous(){
		for(int l=0,L=records.size()-1;l<L;l++){
			int dt=new MDate(records.get(l+1).getTime()).getDT(new MDate(records.get(l).getTime()));
			
			if(Math.abs(dt-21600)>1) return false;
		}
		
		return true;
	}
	
	
	/**
	 * split the float data into segments of a given size
	 */
	public GDPDrifter[] split(int size){
		if(size<1) throw new IllegalArgumentException("size should be positive");
		
		int len=records.size()/size;
		int lef=records.size()%size;
		
		if(lef!=0) len++;
		
		GDPDrifter[] re=new GDPDrifter[len];
		
		int ptr=0;
		
		for(int i=0,I=(lef!=0?len-1:len);i<I;i++){
			re[i]=new GDPDrifter(id+"_"+i,size,meta.length);
			
			re[i].setAttachedMeta(meta);
			
			for(int j=0;j<size;j++)
			re[i].addRecord(records.get(ptr++));
		}
		
		if(lef!=0){
			re[len-1]=new GDPDrifter(id+"_"+(len-1),lef,meta.length);
			
			re[len-1].setAttachedMeta(meta);
			
			for(int j=0,J=lef;j<J;j++)
			re[len-1].addRecord(records.get(ptr++));
		}
		
		return re;
	}
	
	public GDPDrifter[] splitByUndef(){
		boolean hasStr=false,hasEnd=false;
		
		int str=0,end=0;
		
		List<int[]> ls=new ArrayList<int[]>(10);
		
		Record r=records.get(0);
		if(r.getDataValue(UVEL)!=undef&&r.getDataValue(VVEL)!=undef){
			hasStr=true; str=0;
		}
		
		for(int l=1,L=getTCount();l<L;l++){
			r=records.get(l);
			boolean currDef=r.getDataValue(UVEL)!=undef&&r.getDataValue(VVEL)!=undef;
			r=records.get(l-1);
			boolean prevDef=r.getDataValue(UVEL)!=undef&&r.getDataValue(VVEL)!=undef;
			
			if(prevDef&&!currDef&&hasStr){
				hasEnd=true; end=l-1;
			}else if(!prevDef&&currDef){
				hasStr=true; str=l;
			}
			
			if(hasEnd&&hasStr){
				ls.add(new int[]{str,end});
				hasEnd=false;
				hasStr=false;
			}
			
			if(l==L-1&&hasStr&&currDef) ls.add(new int[]{str,L-1});
		}
		
		GDPDrifter[] re=new GDPDrifter[ls.size()];
		
		for(int i=0,I=ls.size();i<I;i++){
			int[] tags=ls.get(i);
			
			re[i]=subRecord(id+"_"+i,tags[0],tags[1]-tags[0]+1);
		}
		
		return re;
	}
	
	public GDPDrifter[] splitByDrogueOffDate(AttachedMeta meta){
		int dc=0,uc=0;
		
		for(int l=0,L=records.size();l<L;l++){
			float drogueState=records.get(l).getDataValue(meta);
			
			if(drogueState==1) dc++;
			else if(drogueState==-1) uc++;
			else throw new IllegalArgumentException("unknown drogue state, 1 for drogued and -1 for undrogued");
		}
		
		GDPDrifter[] re=new GDPDrifter[2];
		
		if(dc!=0) re[0]=subRecord(id+"d",0,dc);
		if(uc!=0) re[1]=subRecord(id+"u",dc,uc);
		
		return re;
	}
	
	
	/**
	 * get a subList of records into a new GDPDrifter
	 * sharing the same memory
	 */
	public GDPDrifter subRecord(String newID,int str,int len){
		GDPDrifter p=new GDPDrifter(newID,len,meta.length);
		
		for(int l=str,L=str+len;l<L;l++) p.records.add(records.get(l));
		
		//p.records=records.subList(str,str+len);
		
		p.setAttachedMeta(meta);
		
		return p;
	}
	
	public GDPDrifter subRecord(int str,int len){ return subRecord(id,str,len);}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(String.format("%d%03d",123456,1111));
	}*/
}
