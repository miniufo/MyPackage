/**
 * @(#)AccessArgoNC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.database;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import ucar.ma2.Array;
import ucar.nc2.NetcdfFile;
import miniufo.diagnosis.MDate;
import miniufo.lagrangian.ArgoFloat;
import miniufo.lagrangian.Record;


/**
 * access Argo float data in NC format
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AccessArgoNC{
	//
	private static final float fillValue=999999.0f;
	private static final float ArgoUndef=99999.0f;
	
	private static final long  refLongTime=19500101000000L;
	private static final long undefLongTime=toLongTime(fillValue);
	
	private static final String  idname ="PLATFORM_NUMBER";
	private static final String cycleNum="CYCLE_NUMBER";
	private static final String lonname ="LONGITUDE";
	private static final String latname ="LATITUDE";
	private static final String timname ="JULD";
	
	private static final MDate  refTime =new MDate(refLongTime);
	
	
	/**
	 * prevent from construction
	 */
	private AccessArgoNC(){}
	
	
	/**
	 * Get basic info (time, lon, lat) from NC file into an ArrayList.
	 * The list may already contain some ArgoFloat.  If an ArgoFloat is
	 * already in the list, then the data will be appended to it.  Else
	 * the method would add a new ArgoFloat into the list.
	 * 
	 * @param	afs			an ArgoFloat list
	 * @param	filename	NC file name
	 */
	public static void parseBasicInfo(List<ArgoFloat> afs,String filename,int dataLen){
		NetcdfFile file=null;
		
		try{ file=NetcdfFile.open(filename);}
		catch(Exception e){ e.printStackTrace(); System.exit(0);}
		
		  int[] cycN=getCycleNumbers(file);
		float[] lons=getLons(file);	toEast(lons);
		float[] lats=getLats(file);
		long[] tims=getTimes(file);
		String[] ids=getIDs(file);
		
		if(ids.length==lons.length){	// for multiple ID floats data
			for(int i=0,I=ids.length;i<I;i++){
				boolean inList=false;
				
				for(int j=0;j<afs.size();j++){
					ArgoFloat af=afs.get(j);
					
					if(ids[i]==af.getID()){
						inList=true;
						if(lons[i]!=ArgoUndef&&lats[i]!=ArgoUndef)
						af.addRecord(new Record(tims[i],lons[i],lats[i],dataLen));
						break;
					}
				}
				
				if(!inList){
					ArgoFloat af=new ArgoFloat(ids[i],dataLen);
					af.setAttachedMeta(ArgoFloat.UVEL,ArgoFloat.VVEL);
					
					if(lons[i]!=ArgoUndef&&lats[i]!=ArgoUndef)
					af.addRecord(new Record(tims[i],lons[i],lats[i],dataLen));
					
					afs.add(af);
				}
			}
		}else{	// for single float multiple data, ids.length==1
			boolean inList=false;
			
			for(int j=0;j<afs.size();j++){
				ArgoFloat af=afs.get(j);
				
				if(ids[0]==af.getID()){
					inList=true;
					for(int i=0,I=lons.length;i<I;i++)
					if(tims[i]!=undefLongTime&&lons[i]!=ArgoUndef){
						Record r=new Record(tims[i],lons[i],lats[i],dataLen);
						r.setCycleNum(cycN[i]);
						af.addRecord(r);
					}
					break;
				}
			}
			
			if(!inList){
				ArgoFloat af=new ArgoFloat(ids[0],dataLen);
				af.setAttachedMeta(ArgoFloat.UVEL,ArgoFloat.VVEL);
				
				for(int i=0,I=lons.length;i<I;i++)
				if(tims[i]!=undefLongTime&&lons[i]!=ArgoUndef){
					Record r=new Record(tims[i],lons[i],lats[i],dataLen);
					r.setCycleNum(cycN[i]);
					af.addRecord(r);
				}
				
				afs.add(af);
			}
		}
		
		try{ file.close();}
		catch(Exception e){ e.printStackTrace(); System.exit(0);}
	}
	
	public static void parseBasicInfo(List<ArgoFloat> afs,String filename){
		parseBasicInfo(afs,filename,2);
	}
	
	/**
	 * whether the given list of ArgoFloat has two identical
	 * ArgoFloat (identical IDs) 
	 * 
	 * @param	afs			an ArgoFloat list
	 * @param	filename	NC file name
	 */
	public static boolean hasSameFloat(List<ArgoFloat> afs){
		for(int i=0,I=afs.size()-1;i<I;i++){
			String id=afs.get(i).getID();
			
			for(int j=i+1,J=afs.size();j<J;j++)
			if(id.equalsIgnoreCase(afs.get(j).getID())) return true;
		}
		
		return false;
	}
	
	public static void writeGS(List<ArgoFloat> afs,String path){
		StringBuffer sb=new StringBuffer();
		sb.append("'sdfopen d:/Data/uwnd.2010.4dl.500.nc'\n");
		sb.append("'enable print "+path+"trajectory.gmf'\n\n");
		sb.append("'setvpage 1.3 1.1 1 1'\n");
		sb.append("'setlopts 8 0.2 60 30'\n\n");
		sb.append("'set line 2 1 4'\n\n");
		for(ArgoFloat af:afs){
			sb.append("'tctrack uwnd "+path+af.getID()+".txt'\n");
		}
		sb.append("\n'draw title Argo trajectories'\n\n");
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		sb.append("'disable print'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		
		try{
			FileWriter fw=new FileWriter(new File(path+"trajectory.gs"));
			fw.write(sb.toString());	fw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	private static void toEast(float[] lons){
		for(int i=0,I=lons.length;i<I;i++)
		if(lons[i]<0) lons[i]=360f+lons[i];
	}
	
	private static int[] getCycleNumbers(NetcdfFile file){
		Array ary=readAll(cycleNum,file);
		
		return (int[])ary.get1DJavaArray(int.class);
	}
	
	
	/*** helper methods ***/
	private static long toLongTime(float daytime){
		// process date field
		int   intPart=(int)daytime;
		float decPart=daytime-(float)intPart;
		float tmpPart=0;
		
		MDate re=refTime.addDays(intPart);
		
		// process hour field
		tmpPart=decPart*24f;
		intPart=(int)tmpPart;
		decPart=tmpPart-(float)intPart;
		re.addHours(intPart);
		
		// process minute field
		tmpPart=decPart*60f;
		intPart=(int)tmpPart;
		decPart=tmpPart-(float)intPart;
		re.addMinutes(intPart);
		
		// process second field
		tmpPart=decPart*60f;
		intPart=(int)tmpPart;
		decPart=tmpPart-(float)intPart;
		re.addSeconds(intPart);
		
		return re.getLongTime();
	}
	
	private static float[] getLons(NetcdfFile file){
		Array ary=readAll(lonname,file);
		
		return (float[])ary.get1DJavaArray(float.class);
	}
	
	private static float[] getLats(NetcdfFile file){
		Array ary=readAll(latname,file);
		
		return (float[])ary.get1DJavaArray(float.class);
	}
	
	private static long[] getTimes(NetcdfFile file){
		Array ary=readAll(timname,file);
		
		float[] data=(float[])ary.get1DJavaArray(float.class);
		
		long[] time=new long[data.length];
		
		for(int i=0,I=data.length;i<I;i++)
		time[i]=toLongTime(data[i]);
		
		return time;
	}
	
	private static String[] getIDs(NetcdfFile file){
		Array ary=readAll(idname,file);
		
		if(ary.getRank()==2){
			char[][] data=(char[][])ary.copyToNDJavaArray();
			
			String[] ids=new String[data.length];
			
			for(int m=0,M=data.length;m<M;m++)
			ids[m]=new String(data[m]).trim();
			
			return ids;
			
		}else{
			char[] data=(char[])ary.copyToNDJavaArray();
			
			String[] ids=new String[1];
			
			ids[0]=new String(data).trim();
			
			return ids;
		}
	}
	
	private static Array readAll(String varname,NetcdfFile file){
		Array ary=null;
		try{ ary=file.findVariable(varname).read();}
		catch(IOException e){e.printStackTrace(); System.exit(0);}
		
		return ary;
	}
	
	
	/** test
	public static void main(String[] args){
		final String[] centers={"aoml","coriolis","csio","incois","jma","kma","kordi","meds"};
		//final String[] centers={"aoml"};
		
		StringBuffer sb=new StringBuffer();
		sb.append("'sdfopen d:/Data/uwnd.2010.4dl.500.nc'\n");
		sb.append("'enable print d:/Data/Argo/Traj/trajectory.gmf'\n\n");
		sb.append("'setvpage 1.3 1.1 1 1'\n");
		sb.append("'setlopts 8 0.2 60 30'\n\n");
		sb.append("'set line 2 1 1'\n\n");
		
		for(String center:centers){
			ArrayList<ArgoFloat> afs=new ArrayList<ArgoFloat>();
			
			File[] fs=new File("d:/Data/Argo/Traj/"+center).listFiles();
			
			System.out.print(String.format("%8s has %5d files:  ",center,fs.length));
			
			for(File f:fs) if(!f.isDirectory()) parseBasicInfo(afs,f.getAbsolutePath());
			
			int osize=afs.size();
			int nsize=0;
			
			for(ArgoFloat af:afs){
				if(af.getID()==2900586) continue;
				if(af.getID()==2900662) continue;
				if(af.getID()==3900465) continue;
				if(af.getID()==4901072) continue;
				if(af.getID()==4900541) continue;
				if(af.getID()==5900488) continue;
				if(af.getID()==5900943) continue;
				if(af.getID()==5901116) continue;
				if(af.getID()==5901286) continue;
				if(af.getID()==6900121) continue;
				if(af.getID()==6900358) continue;
				if(af.getID()==6900371) continue;
				if(af.getID()==6900384) continue;
				if(af.getID()==7900045) continue;
				if(af.getID()==7900046) continue;
				
				if(af.getRecordLength()<2) continue;
				
				af.sort();
				af.crossIDLToContinuousRecord();
				
				ArgoFloat af2=af.toDailyData(); if(af2.getRecordLength()<2) continue;
				af2.interpolateDailyPosition();
				af2.cDailyCurrentSpeed();
				af2.truncate(20070101000000L,20071231000000L);
				af2.crossIDLToDiscontinuousRecord();
				
				if(af2.getRecordLength()==365){
					af2.toTrajectoryFile("d:/Data/Argo/Traj/TXT/");
					sb.append("'tctrack uwnd d:/Data/Argo/Traj/TXT/"+af.getID()+".txt'\n");
					nsize++;
				}
			}
			
			System.gc();
			System.out.println(String.format("%5d/%5d (%6.2f%%)",nsize,osize,((float)nsize/osize*100f)));
		}
		
		sb.append("\n'draw title Argo trajectories'\n\n");
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		sb.append("'disable print'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		
		try{
			FileWriter fw=new FileWriter(new File("d:/Data/Argo/Traj/trajectory.gs"));
			fw.write(sb.toString());	fw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}*/
}
