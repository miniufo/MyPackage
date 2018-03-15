/**
 * @(#)AccessGDPDrifter.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.database;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import miniufo.diagnosis.MDate;
import miniufo.lagrangian.DrogueOffData;
import miniufo.lagrangian.GDPDrifter;
import miniufo.lagrangian.MetaData;
import miniufo.lagrangian.Record;
import miniufo.util.Region2D;


/**
 * parse GISST into binary
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AccessGDPDrifter{
	//
	private static final float GDPUndef=999.999f;
	
	
	/**
	 * prevent from construction
	 */
	private AccessGDPDrifter(){}
	
	
	/**
	 * collect the information from files into a list
	 * (if one snapshot of drifter is within a region, the whole track record is returned)
	 *
	 * @param	drftr		drifter records
	 * @param	filename	file name for data record
	 * @param	dataLen		length of the data attached to one single record
	 * @param	lon1		west longitude of the region (degree)
	 * @param	lon2		east longitude of the region (degree)
	 * @param	lat1		south latitude of the region (degree)
	 * @param	lat2		north latitude of the region (degree)
	 */
	public static void parseBasicGDPInfo(List<GDPDrifter> drftrs,String filename,int dataLen,Region2D region){
		if(dataLen<3) throw new IllegalArgumentException("dataLen should be at least 3 for [uvel,vvel,temp]");
		
		GDPDrifter prev=null;
		GDPDrifter curr=null;
		
		try(BufferedReader br=new BufferedReader(new FileReader(filename),8192*8)){
			String oneline=br.readLine();
			
			while(oneline!=null){
				String newID=oneline.substring(0,8).trim();
				
				if(prev!=null&&prev.getID().equalsIgnoreCase(newID)){
					prev.addRecord(constructRecord(oneline,dataLen));
					
				}else{
					curr=new GDPDrifter(newID,dataLen);
					curr.setAttachedMeta(GDPDrifter.UVEL,GDPDrifter.VVEL,GDPDrifter.Temp,GDPDrifter.DrgOff);
					curr.addRecord(constructRecord(oneline,dataLen));
					
					if(prev!=null&&inRange(prev,region)) drftrs.add(prev);
					
					prev=curr;
					curr=null;
				}
				
				oneline=br.readLine();
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public static void parseBasicGDPInfo(List<GDPDrifter> drftrs,String filename,Region2D region){
		parseBasicGDPInfo(drftrs,filename,4,region);
	}
	
	public static void parseBasicGDPInfo(List<GDPDrifter> drftrs,String filename,int dataLen){
		if(dataLen<3) throw new IllegalArgumentException("dataLen should be at least 3 for [uvel,vvel,temp]");
		
		try(BufferedReader br=new BufferedReader(new FileReader(filename),8192*8)){
			String oneline=br.readLine();
			
			while(oneline!=null){
				String newID=oneline.substring(0,8).trim();
				
				boolean inList=false;
				
				for(int i=drftrs.size()-1;i>=0;i--){
					GDPDrifter drftr=drftrs.get(i);
					
					if(drftr.getID().equalsIgnoreCase(newID)){
						inList=true;
						drftr.addRecord(constructRecord(oneline,dataLen));
						break;
					}
				}
				
				if(!inList){
					GDPDrifter dfr=new GDPDrifter(newID,dataLen);
					dfr.setAttachedMeta(GDPDrifter.UVEL,GDPDrifter.VVEL,GDPDrifter.Temp,GDPDrifter.DrgOff);
					dfr.addRecord(constructRecord(oneline,dataLen));
					drftrs.add(dfr);
				}
				
				oneline=br.readLine();
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public static void parseBasicGDPInfo(List<GDPDrifter> drftrs,String filename){
		parseBasicGDPInfo(drftrs,filename,4);
	}
	
	
	/**
	 * collect the MetaData from files into a list
	 *
	 * @param	ls		a list of MetaData returned as the result
	 * @param	filename	file name for data record
	 */
	public static void parseGDPMetaData(List<MetaData> ls,String filename){
		try(BufferedReader br=new BufferedReader(new FileReader(filename))){
			String oneline=br.readLine();
			
			while(oneline!=null){
				ls.add(new MetaData(oneline));
				
				oneline=br.readLine();
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	/**
	 * collect the Drogue-off data from files into a list
	 *
	 * @param	ls		a list of MetaData returned as the result
	 * @param	filename	file name for data record
	 */
	public static void parseDrogueOffData(List<DrogueOffData> ls,String filename){
		try(BufferedReader br=new BufferedReader(new FileReader(filename))){
			String oneline=br.readLine();
			oneline=br.readLine();
			oneline=br.readLine();
			oneline=br.readLine();
			
			while(oneline!=null){
				ls.add(new DrogueOffData(oneline));
				
				oneline=br.readLine();
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/**
	 * attach drogue data from meta-files, 1 for drogued and -1 for undrogued
	 *
	 * @param	drlst	list of drifters
	 * @param	mtlst	list of meta-data
	 */
	public static void attachDrogueData(List<GDPDrifter> drlst,List<MetaData> mtlst){
		for(GDPDrifter dr:drlst){
			int id=Integer.parseInt(dr.getID());
			
			boolean hasMD=false;
			
			for(MetaData md:mtlst)
			if(id==md.getID()){
				if(md.getUndroguedTime()==MetaData.undefinedTime){
					for(int l=0,L=dr.getTCount();l<L;l++) dr.getRecord(l).setData(GDPDrifter.DrgOff,1);
					
				}else for(int l=0,L=dr.getTCount();l<L;l++){
					long time=dr.getRecord(l).getTime();
					
					//if(time!=MetaData.undefinedTime){
						if(time<md.getUndroguedTime()) dr.getRecord(l).setData(GDPDrifter.DrgOff,1);
						else dr.getRecord(l).setData(GDPDrifter.DrgOff,-1);
						/*
						if(time>=(md.getDeployTime()-40000L)&&time<=(md.getEndTime()+50000L)){
							if(time<md.getUndroguedTime()) dr.getRecord(l).setData(3,1);
							else dr.getRecord(l).setData(3,-1);
							
						}else throw new IllegalArgumentException(
							"drifter ("+dr.getID()+") at "+time+" is not within [" +
							md.getDeployTime()+", "+md.getEndTime()+"]"
						);*/
						
					//}else throw new IllegalArgumentException("drifter ("+dr.getID()+") has undefined time");
				}
				
				hasMD=true;
				break;
			}
			
			if(!hasMD) throw new IllegalArgumentException("drifter ("+dr.getID()+") has no metadata");
		}
	}
	
	
	/**
	 * divided the drifter records into segments, and return those segments
	 * that are completely within a given region
	 * (all records are within the region)
	 *
	 * @param	drftr	drifter records
	 * @param	lon1	west longitude of the region (degree)
	 * @param	lon2	east longitude of the region (degree)
	 * @param	lat1	south latitude of the region (degree)
	 * @param	lat2	north latitude of the region (degree)
	 */
	public static GDPDrifter[] getRecordsWithinRegion(GDPDrifter drftr,Region2D region){
		int str=0,end=0;
		
		boolean hasStr=false,hasEnd=false;
		
		if(drftr.getTCount()==0) throw new IllegalArgumentException("no records of this drifter");
		
		float[] lons=drftr.getXPositions();
		float[] lats=drftr.getYPositions();
		
		List<int[]> ls=new ArrayList<int[]>(10);
		
		if(region.inRange(lons[0],lats[0])){
			hasStr=true; str=0;
		}
		
		for(int l=1,L=drftr.getTCount();l<L;l++){
			boolean currIn=region.inRange(lons[l  ],lats[l  ]);
			boolean prevIn=region.inRange(lons[l-1],lats[l-1]);
			
			if(prevIn&&!currIn&&hasStr){
				hasEnd=true; end=l-1;
			}else if(!prevIn&&currIn){
				hasStr=true; str=l;
			}
			
			if(hasEnd&&hasStr){
				ls.add(new int[]{str,end});
				hasEnd=false;
				hasStr=false;
			}
			
			if(l==L-1&&hasStr&&currIn) ls.add(new int[]{str,L-1});
		}
		
		GDPDrifter[] re=new GDPDrifter[ls.size()];
		
		String id=drftr.getID();
		for(int i=0,I=ls.size();i<I;i++){
			int[] tags=ls.get(i);
			
			if(ls.size()==1&&tags[0]==0&&tags[1]==(drftr.getTCount()-1))
				re[0]=drftr.subRecord(id,tags[0],tags[1]-tags[0]+1);
			else
				re[i]=drftr.subRecord(id+"_"+i,tags[0],tags[1]-tags[0]+1);
		}
		
		return re;
	}
	
	public static List<GDPDrifter> getRecordsWithinRegion(List<GDPDrifter> all,Region2D region){
		List<GDPDrifter> subls=new ArrayList<>();
		
		for(GDPDrifter drftr:all){
			GDPDrifter[] re=getRecordsWithinRegion(drftr,region);
			
			for(GDPDrifter p:re)
			if(p.getTCount()!=0) subls.add(p);
		}
		
		return subls;
	}
	
	/**
	 * get records within the specified temporal range
	 *
	 * @param	drftr	drifter records
	 * @param	tstr	start of the time
	 * @param	tend	end of the time
	 */
	public static List<GDPDrifter> getRecordsWithinRange(List<GDPDrifter> all,MDate str,MDate end){
		long tstr=str.getLongTime();
		long tend=end.getLongTime();
		
		if(tstr>tend) throw new IllegalArgumentException("T-start should be before T-end");
		
		List<GDPDrifter> res=new ArrayList<>();
		
		for(GDPDrifter drftr:all){
			int stag=-1,len=0;
			boolean lastIn=false;
			
			for(int l=0,L=drftr.getTCount();l<L;l++){
				long time=drftr.getRecord(l).getTime();
				
				boolean currIn=time>=tstr&&time<=tend;
				
				if(!lastIn&&currIn) stag=l;
				if(stag!=-1&&currIn) len++;
				
				lastIn=currIn;
			}
			
			if(stag!=-1) res.add(drftr.subRecord(stag,len));
		}
		
		return res;
	}
	
	
	/**
	 * helper methods
	 **/
	private static Record constructRecord(String oneline,int dataLen){
		int mo=Integer.parseInt(oneline.substring(10,13).trim());
		int dy=Integer.parseInt(oneline.substring(14,16).trim());
		int hr=Integer.parseInt(oneline.substring(17,19));
		int yr=Integer.parseInt(oneline.substring(21,25).trim());
		
		if(hr==0) hr=0;
		else if(hr==25) hr=6;
		else if(hr==50) hr=12;
		else if(hr==75) hr=18;
		else throw new IllegalArgumentException("invalid day ("+hr+"), should be [.000 .250 .500 .750]");
		
		Record r=new Record(
			new MDate(yr,mo,dy,hr,0).getLongTime(),
			Float.parseFloat(oneline.substring(38,45)),Float.parseFloat(oneline.substring(28,35)),
			dataLen
		);
		
		float u=Float.parseFloat(oneline.substring(57,65));
		float v=Float.parseFloat(oneline.substring(67,75));
		float t=Float.parseFloat(oneline.substring(47,56));
		
		r.setData(GDPDrifter.UVEL,u==GDPUndef?Record.undef:u/100f);
		r.setData(GDPDrifter.VVEL,v==GDPUndef?Record.undef:v/100f);
		r.setData(GDPDrifter.Temp,t==GDPUndef?Record.undef:t     );
		
		return r;
	}
	
	private static boolean inRange(GDPDrifter dft,Region2D region){
		for(int l=0,L=dft.getTCount();l<L;l++){
			Record r=dft.getRecord(l);
			
			float lon=r.getXPos();
			float lat=r.getYPos();
			
			if(region.inRange(lon,lat)) return true;
		}
		
		return false;
	}
	
	
	/** test
	public static void main(String[] args){
		String path="d:/Data/GDP/";
		
		List<MetaData> mdls=new ArrayList<>();
		//parseGDPMetaData(mdls,path+"dirfl_1_5000.dat");
		parseGDPMetaData(mdls,path+"dirfl_5001_sep12.dat");
		
		List<DrogueOffData> dfls=new ArrayList<>();
		parseDrogueOffData(dfls,path+"2013_02_14_drogueoff.dat");
		
	aa:	for(int l=0,L=mdls.size();l<L;l++){
			int id=mdls.get(l).getID();
			
			for(DrogueOffData dod:dfls){
				if(id==dod.getID()){
					System.out.println(
						mdls.get(l).getUndroguedTime()+"\t"+
						dod.getOriginalTime()+"\t"+
						dod.getAutomaticTime()+"\t"+
						dod.getManualTime()
					);
					
					continue aa;
				}
			}
			
			System.out.println(id+" has no drogue-off data");
		}
	}*/
}
