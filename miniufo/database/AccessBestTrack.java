/**
 * @(#)AccessBestTrack.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.database;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.ArrayList;
import java.util.function.Predicate;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;


/**
 * access typhoon record files from the "best track" database
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AccessBestTrack{
	//
	public enum DataSets{CMA,GUAM,JMA,JTWC,NHC};
	
	public enum JTWCBasins{CP,EP,WP,IO,SH};
	
	
	/**
	 * prevent from construction
	 */
	private AccessBestTrack(){}
	
	
	/**
     * Write a list of typhoon to a file.
     *
     * @param	typhoons	a list of Typhoon
     */
	public static void recordsToFile(List<Typhoon> ls,String path){
		StringBuilder sb=new StringBuilder();
		
		for(Typhoon ty:ls) sb.append(ty.toString());
		
		try(FileWriter fw=new FileWriter(path)){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	/**
     * glue JTWC best-track source files (one TC per file) into a single file
     *
     * @param	basin	basin of JTWC data
     */
	public static void glueJTWCSourceFiles(JTWCBasins basin){
		String path="D:/Data/Typhoons/JTWC/";
		
		try(FileWriter fw=new FileWriter(path+"b"+basin.toString().toLowerCase()+".dat")){
			Files.list(Paths.get(path+"original/b"+basin.toString().toLowerCase()+"/")).forEach(p->{
				try{
					long count=Files.lines(p).count();
					
					StringBuilder sb=new StringBuilder();
					
					if(count!=0){
						sb.append(count+"\n");
						Files.lines(p).forEach(line->sb.append(line+"\n"));
						
					}else System.out.println("file "+p+" has 0 record");
					
					fw.write(sb.toString());
					
				}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			});
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/**
     * get Typhoons from "best track" dataset
     *
     * @param	path	path of the the dataset
     * @param	cond	condition for typhoon record typical form:
	 *					"name=chanchu;num=0001;id=0601;time=5jun1997-9jul1997;lon=80-270;lat=5-25"
     * @param	ds		dataset
     */
	public static List<Typhoon> getTyphoons(String path,String cond,DataSets ds){
		switch(ds){
		case CMA : return getTyphoonsFromCMA(path,cond);
		case JMA : return getTyphoonsFromJMA(path,cond);
		case JTWC: return getTyphoonsFromJTWC(path,cond);
		case NHC : return getTyphoonsFromNHC(path,cond);
		case GUAM: return getTyphoonsFromGUAM(path,cond);
		default  : throw new IllegalArgumentException("unsupported data set");
		}
	}
	
	public static List<Typhoon> getTyphoonsFromCMA(String path,String cond){
		System.out.println("Getting records from CMA ("+path+")");
		
		// to get all records from database
		List<Typhoon> all=new ArrayList<>();
		
		try(BufferedReader br=new BufferedReader(new FileReader(new File(path)))){
			int cc=0;
			while(true){
				String ln=br.readLine();
				if(ln==null) break;
				
				/********* for header *********/
				// Number of data lines
				int count=Integer.parseInt(ln.substring(12,15).trim());
				
				// tropical cyclone number ID
				//int  num =Integer.parseInt(ln.substring(16,20).trim());
				
				// Chinese number ID
				String id=ln.substring(21,25);
				
				// name of the storm, '(-)n' means sub-center and its number
				String name=ln.substring(30,50).trim();
				if(name.length()==0) name="nameless";
				
				if(name.indexOf("(-)")!=-1){
					for(int l=0;l<count;l++) br.readLine();	// skip sub-center records
					continue;
				}
				
				Typhoon typhoon=new Typhoon(id,name,count,5);
				
				/********* for records *********/
				for(int l=0;l<count;l++){
					ln=br.readLine();
					
					// date
					long time=Long.parseLong(ln.substring(0,10)+"0000");
					
					// intensity category according to "Chinese National Standard
					// for Grade of Tropical Cyclones", which was put in practice
					// since 15 June 2006
					// 0 weaker than TD or unknown
					// 1 tropical depression	(TC , 10.8-17.1 m/s)
					// 2 tropical storm			(TS , 17.2-24.4 m/s)
					// 3 severe tropical storm	(STS, 24.5-32.6 m/s)
					// 4 typhoon				(TY , 41.4-32.7 m/s)
					// 5 severe typhoon			(STY, 41.5-50.9 m/s)
					// 6 super typhoon			(superTC, >=51.0 m/s)
					// 9 extratropical transition cyclone (ET)
					TYPE type=TYPE.OTHERS;
					switch(Integer.parseInt(ln.substring(11,12))){
						case 1:
							type=TYPE.TD;		break;
						case 2:case 3:
							type=TYPE.TS;		break;
						case 4:case 5:case 6:
							type=TYPE.TY;		break;
						case 9:
							type=TYPE.EC;		break;
						case 0:
							type=TYPE.OTHERS;	break;
						default:
							throw new IllegalArgumentException("unknown type:"+ln.substring(11,12));
					}
					
					float lat=Float.parseFloat(ln.substring(13,16))/10f;	// latitude (0.1 degree)
					float lon=Float.parseFloat(ln.substring(17,21))/10f;	// longitude (0.1 degree)
					float prs=Float.parseFloat(ln.substring(22,26));		// minimum central pressure (hPa)
					float wnd=Float.parseFloat(ln.substring(31,34));		// maximum sustained wind speed (m/s)
					
					Record r=new Record(time,lon,lat,5);
					
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					r.setData(3,prs);
					r.setData(4,type.ordinal());
					
					typhoon.setAttachedDataNames("uvel","vvel","maxwspd","minslp","type");
					
					typhoon.addRecord(r);
				}
				
				all.add(typhoon);
				
				if(cc++%100==0) System.out.print(".");
			}
			
			System.out.println();
			
		}catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		// to select the records which meet the requires
		
		List<Typhoon> res=getTyphoons(all,cond);
		
		validate(res,r->true);
		
		return res;
	}
	
	public static List<Typhoon> getTyphoonsFromGUAM(String path,String cond){
		System.out.println("Getting records from GUAM ("+path+")");
		
		// typical form: "name=chanchu;num=0001;id=0601;time=5jun1997-9jul1997;lon=80-270;lat=5-25;"
		// area: ep  - eastern Pacific
		// area: cp  - central Pacific
		// area: wp  - western Pacific
		// area: io  - Indian ocean
		// area: sh  - southern hemisphere
		// area: all - all ocean basin
		
		List<Typhoon> all=new ArrayList<Typhoon>();
		
		try(BufferedReader br=new BufferedReader(new FileReader(new File(path)))){
			while(true){
				String ln=br.readLine();
				if(ln==null) break;
				
				int count=Integer.parseInt(ln);
				
				String id=null;
				
				List<Record> ls=new ArrayList<>(count);
				
				for(int i=0;i<count;i++){
					String oneline=br.readLine();
					
					String[] linesplit=oneline.split(",");
					
					//int num=Integer.parseInt(linesplit[1].trim());
					
					float lat=0;	String tmp=linesplit[6].trim();
					if(tmp.indexOf("N")!=-1) lat=Float.parseFloat(tmp.replace("N",""))/10f;
					else lat=Float.parseFloat(tmp.replace("S",""))/-10f;
					
					float lon=0;	tmp=linesplit[7].trim();
					if(tmp.indexOf("E")!=-1) lon=Float.parseFloat(tmp.replace("E",""))/10f;
					else lon=360f-Float.parseFloat(tmp.replace("W",""))/10f;
					
					// unit changes from knot to m/s
					float wnd=Float.parseFloat(linesplit[8].trim())*0.51444f;
					
					tmp=linesplit[2].trim();
					
					long time=Long.parseLong(tmp+"0000");
					
					if(i==0) id=tmp.substring(0,4)+linesplit[1].trim();
					
					Record r=new Record(time,lon,lat,3);
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					
					ls.add(r);
				}
				
				Typhoon ty=new Typhoon(id,"nameless",ls);
				ty.setAttachedDataNames("uvel","vvel","maxwspd");
				
				all.add(ty);
			}
			
		}catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		// to select the records which meet the requires
		List<Typhoon> res=getTyphoons(all,cond);
		
		validate(res,r->true);
		
		return res;
	}
	
	public static List<Typhoon> getTyphoonsFromJMA(String path,String cond){
		System.out.println("Getting records from JMA ("+path+")");
		
		// typical form: "name=chanchu;num=0001;id=0601;time=5jun1997-9jul1997;lon=80-270;lat=5-25;"
		
		List<Typhoon> all=new ArrayList<Typhoon>();
		
		try(BufferedReader br=new BufferedReader(new FileReader(new File(path)))){
			int cc=0;
			while(true){
				String ln=br.readLine();
				if(ln==null) break;
				
				/********* for header *********/
				// tropical cyclone number ID, Serial number ID of the storm
				// of intensity with maximum sustained wind speed of 28 kt (near gale) or greater
				//int num=Integer.parseInt(ln.substring(7,10).trim());
				
				// Number of data lines
				int count=Integer.parseInt(ln.substring(12,15).trim());
				
				// international number ID, replicate
				String id=ln.substring(21,25);
				
				// name of the storm
				String name=ln.substring(30,50).trim();
				if(name.length()==0) name="nameless";
				
				Typhoon typhoon=new Typhoon(id,name,count,5);
				
				for(int l=0;l<count;l++){
					ln=br.readLine();
					
					String tm=ln.substring(0,8)+"0000";
					
					long time=0;
					if(Integer.parseInt(tm.substring(0,2))>40) time=Long.parseLong("19"+tm);
					else time=Long.parseLong("20"+tm);
					
					// grade
					// 1 not used
					// 2 tropical depression	(TD)
					// 3 tropical storm			(TS)
					// 4 severe tropical storm	(STS)
					// 5 typhoon				(TY)
					// 6 extra-tropical cyclone	(L)
					// 7 just entering into the
					//	 responsible area of JMA
					// 8 not used
					// 9 tropical cyclone of TS intensity or higher
					TYPE type=TYPE.OTHERS;
					switch(Integer.parseInt(ln.substring(13,14))){
						case 3:case 4:
							type=TYPE.TS;
							break;
						case 5:case 9:
							type=TYPE.TY;
							break;
						case 6:
							type=TYPE.EC;
							break;
						case 2:
							type=TYPE.TD;
							break;
						case 1:case 7:case 8:
							type=TYPE.OTHERS;
							break;
						default:
							throw new IllegalArgumentException("unknown type:"+ln.substring(13,14));
					}
					
					float lat=Float.parseFloat(ln.substring(15,18))/10f;	// latitude (0.1 degree)
					float lon=Float.parseFloat(ln.substring(19,23))/10f;	// longitude (0.1 degree)
					float prs=Float.parseFloat(ln.substring(24,28).trim());	// minimum central pressure (hPa)
					float wnd=0;						// maximum sustained wind speed (knot changed to m/s)
					
					String wind=ln.substring(33,36).trim();
					if(wind.length()!=0) wnd=Float.parseFloat(wind)*0.51444f;
					
					Record r=new Record(time,lon,lat,5);
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					r.setData(3,prs);
					r.setData(4,type.ordinal());
					
					typhoon.setAttachedDataNames("uvel","vvel","maxwspd","minslp","type");
					
					typhoon.addRecord(r);
				}
				
				all.add(typhoon);
				
				if(cc++%100==0) System.out.print(".");
			}
			
			System.out.println();
			
		}catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		// to select the records which meet the requires
		List<Typhoon> res=getTyphoons(all,cond);
		
		validate(res,r->true);
		
		return res;
	}
	
	public static List<Typhoon> getTyphoonsFromJTWC(String path,String cond){
		System.out.println("Getting records from JTWC ("+path+")");
		
		// typical form: "name=chanchu;num=0001;id=0601;time=5jun1997-9jul1997;lon=80-270;lat=5-25;"
		// area: ep  - eastern Pacific
		// area: cp  - central Pacific
		// area: wp  - western Pacific
		// area: io  - Indian ocean
		// area: sh  - southern hemisphere
		// area: all - all ocean basin
		
		List<Typhoon> all=new ArrayList<Typhoon>();
		
		try(BufferedReader br=new BufferedReader(new FileReader(new File(path)))){
			int cc=0;
			while(true){
				String ln=br.readLine();
				if(ln==null) break;
				
				int count=Integer.parseInt(ln);
				
				String id=null;
				
				List<Record> ls=new ArrayList<>(count);
				
				for(int i=0;i<count;i++){
					String oneline=br.readLine();
					
					int num=Integer.parseInt(oneline.substring(3,6).trim());
					
					float lat=0;
					if(oneline.substring(38,39).equals("N"))
						lat=Float.parseFloat(oneline.substring(34,38).trim())/10f;
					else
						lat=Float.parseFloat(oneline.substring(34,38).trim())/-10f;
					
					float lon=0;
					if(oneline.substring(45,46).equals("E"))
						lon=Float.parseFloat(oneline.substring(40,45).trim())/10f;
					else
						lon=360f-Float.parseFloat(oneline.substring(40,45).trim())/10f;
					
					if(oneline.length()<51) System.out.println(oneline);
					
					// change unit from knot to m/s
					float wnd=Float.parseFloat(oneline.substring(47,51).trim())*0.51444f;
					float prs=0;
					
					if(oneline.length()>51) prs=Float.parseFloat(oneline.substring(52,57).trim());
					
					long time=Long.parseLong(oneline.substring(8,18)+"0000");
					
					if(i==0) id=oneline.substring(8 ,12)+String.format("%02d",num);
					
					Record r=new Record(time,lon,lat,4);
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					r.setData(3,prs);
					ls.add(r);
				}
				
				Typhoon ty=new Typhoon(id,"nameless",ls);
				ty.setAttachedDataNames("uvel","vvel","maxwspd","minslp");
				
				all.add(ty);
				
				if(cc++%100==0) System.out.print(".");
			}
			
			System.out.println();
			
		}catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		// to select the records which meet the requires
		List<Typhoon> res=getTyphoons(all,cond);
		
		validate(res,r->true);
		
		return res;
	}
	
	public static List<Typhoon> getTyphoonsFromNHC(String path,String cond){
		System.out.println("Getting records from NHC ("+path+")");
		
		// typical form: "name=chanchu;num=0001;id=0601;time=5jun1997-9jul1997;lon=80-270;lat=5-25;"
		// area: ep - eastern Pacific
		// area: at - Atlantic ocean basin
		
		List<Typhoon> all=new ArrayList<Typhoon>();
		
		try(BufferedReader br=new BufferedReader(new FileReader(new File(path)))){
			int cc=0;
			while(true){
				String head=br.readLine();
				if(head==null) break;
				
				int count=Integer.parseInt(head.substring(19,21).trim());
				int num  =Integer.parseInt(head.substring(22,24).trim());
				
				String year=head.substring(12,16).trim();
				String name=head.substring(35,46).trim();
				
				Typhoon typhoon=new Typhoon(year+num,name,count,4);
				
				for(int i=0;i<count;i++){
					String oneline=br.readLine();
					
					String month=oneline.substring(6,8 ).trim();
					String day  =oneline.substring(9,11).trim();
					
					String tmp=null;
					
					float lon=0,lat=0,wnd=0,prs=0;
					
					tmp=oneline.substring(12,15).trim();
					if(tmp.length()!=0) lat=Float.parseFloat(tmp)/10f; else lat=0;
					tmp=oneline.substring(15,19).trim();
					if(tmp.length()!=0) lon=360f-Float.parseFloat(tmp)/10f; else lon=0;
					tmp=oneline.substring(19,23).trim();
					if(tmp.length()!=0) wnd=Float.parseFloat(tmp)*0.51444f;  else wnd=0;
					tmp=oneline.substring(23,28).trim();
					if(tmp.length()!=0) prs=Float.parseFloat(tmp);  else prs=0;
					
					Record r=new Record(Long.parseLong(year+month+day+"000000"),lon,lat,4);
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					r.setData(3,prs);
					typhoon.addRecord(r);
					
					tmp=oneline.substring(29,32).trim();
					if(tmp.length()!=0) lat=Float.parseFloat(tmp)/10f; else lat=0;
					tmp=oneline.substring(32,36).trim();
					if(tmp.length()!=0) lon=360f-Float.parseFloat(tmp)/10f; else lon=0;
					tmp=oneline.substring(36,40).trim();
					if(tmp.length()!=0) wnd=Float.parseFloat(tmp)*0.51444f; else wnd=0;
					tmp=oneline.substring(40,45).trim();
					if(tmp.length()!=0) prs=Float.parseFloat(tmp); else prs=0;
					
					r=new Record(Long.parseLong(year+month+day+"060000"),lon,lat,4);
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					r.setData(3,prs);
					typhoon.addRecord(r);
					
					tmp=oneline.substring(46,49).trim();
					if(tmp.length()!=0) lat=Float.parseFloat(tmp)/10f; else lat=0;
					tmp=oneline.substring(49,53).trim();
					if(tmp.length()!=0) lon=360f-Float.parseFloat(tmp)/10f; else lon=0;
					tmp=oneline.substring(53,57).trim();
					if(tmp.length()!=0) wnd=Float.parseFloat(tmp)*0.51444f; else wnd=0;
					tmp=oneline.substring(57,62).trim();
					if(tmp.length()!=0) prs=Float.parseFloat(tmp); else prs=0;
					
					r=new Record(Long.parseLong(year+month+day+"120000"),lon,lat,4);
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					r.setData(3,prs);
					typhoon.addRecord(r);
					
					tmp=oneline.substring(63,66).trim();
					if(tmp.length()!=0) lat=Float.parseFloat(tmp)/10f; else lat=0;
					tmp=oneline.substring(66,70).trim();
					if(tmp.length()!=0) lon=360f-Float.parseFloat(tmp)/10f; else lon=0;
					tmp=oneline.substring(70,74).trim();
					if(tmp.length()!=0) wnd=Float.parseFloat(tmp)*0.51444f; else wnd=0;
					tmp=oneline.substring(74,79).trim();
					if(tmp.length()!=0) prs=Float.parseFloat(tmp); else prs=0;
					
					r=new Record(Long.parseLong(year+month+day+"180000"),lon,lat,4);
					r.setData(0,0  );
					r.setData(1,0  );
					r.setData(2,wnd);
					r.setData(3,prs);
					typhoon.addRecord(r);
				}
				
				typhoon.setAttachedDataNames("uvel","vvel","maxwspd","minslp");
				
				all.add(typhoon);
				
				br.readLine();
				
				if(cc++%100==0) System.out.print(".");
			}
			
			System.out.println();
			
		}catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		// to select the records which meet the requires
		List<Typhoon> res=getTyphoons(all,cond);
		
		validate(res,r->r.getDataValue(2)==0);
		
		return res;
	}
	
	public static List<Typhoon> getTyphoons(List<Typhoon> all,String cond){
		if(cond.equals("")) return all;
		
		List<Typhoon>   res =new ArrayList<Typhoon>();
		List<Typhoon> delrec=new ArrayList<Typhoon>();
		
		String[] cons=cond.split(";");
		
		for(int i=0;i<cons.length;i++){
			String[] conds=cons[i].split("=");
			
			if(conds.length!=2) throw new IllegalArgumentException(
				"invalid conditions form\n" +
				"typical form:\n" +
				"\"name=chanchu;id=0601;time=5jun1997-9jul1997;lon=80-270;lat=5-25;\""
			);
			
			if(conds[0].equals("name")){
				if(res.size()!=0){
					for(Typhoon one:res)
					if(!one.getName().equalsIgnoreCase(conds[1])) delrec.add(one);
					
					res.removeAll(delrec);
					
				}else{
					for(Typhoon one:all)
					if(one.getName().equalsIgnoreCase(conds[1])) res.add(one);
				}
				
			}else if(conds[0].equals("id")){
				if(res.size()!=0){
					for(Typhoon one:res)
					if(!one.getID().equals(conds[1])) delrec.add(one);
					
					res.removeAll(delrec);
					
				}else{
					for(Typhoon one:all)
					if(one.getID().equals(conds[1])) res.add(one);
				}
				
			}else if(conds[0].equals("time")){
				String[] range=conds[1].split("-");
				
				long str=new MDate(range[0]).getLongTime();
				long end=new MDate(range[1]).getLongTime();
				
				if(str>end) throw new IllegalArgumentException("invlid time range");
				
				if(res.size()!=0){
					for(Typhoon one:res){
						long tb=one.getTime(0);
						
						if(tb<str||tb>=end) delrec.add(one);
					}
					
					res.removeAll(delrec);
					
				}else{
					for(Typhoon one:all){
						long tb=one.getTime(0);
						
						if(tb>=str&&tb<end) res.add(one);
					}
				}
				
			}else throw new IllegalArgumentException(
				"invalid conditions form\n" +
				"typical form:\n" +
				"\"name=chanchu;num=0001;id=0601;time=5jun1997-9jul1997;lon=80-270;lat=5-25;\""
			);
		}
		
		return res;
	}
	
	
	/**
     * filter out the records whose wind speed < threshold
     *
     * @param	res			records in a list
     * @param	threshold	threshold for the wind speed (m s^-1)
	 *						records of wind < threshold will be removed
     */
	public static void maxWindFilter(List<Typhoon> res,float threshold){
		if(res.get(0).getDataNames().length<=2) throw new IllegalArgumentException("no maximum wind");
		
		for(int l=0;l<res.size();l++){
			Typhoon tr=res.get(l);
			
			// delete records do not reach the threshold
			if(tr.removeRecords(r->r.getDataValue(2)<threshold).getTCount()==0){
				res.remove(l);
				System.out.println("remove the "+l+"(th) record");
				l--;
			}
		}
	}
	
	/**
     * filter out the records whose minimum pressure > threshold
     *
     * @param	res			records in a list
     * @param	threshold	threshold for the wind speed (m s^-1)
	 *						records of wind < threshold will be removed
     */
	public static void minPressFilter(List<Typhoon> res,float threshold){
		if(res.get(0).getDataNames().length<=3) throw new IllegalArgumentException("no minimum pressure");
		
		for(int l=0;l<res.size();l++){
			Typhoon tr=res.get(l);
			
			if(tr.removeRecords(r->r.getDataValue(3)>threshold).getTCount()==0){
				res.remove(l);
				System.out.println("remove the "+l+"(th) record");
				l--;
			}
		}
	}
	
	
	/**
     * filter out the TCs that is birth (wind > threshold) outside a specific region
     *
     * @param	tr		typhoon record
     * @param	thre	threshold of the wind larger than which a TC is said to give birth
     * @param	lon1	west longitude for a region
     * @param	lat1	south latitude for a region
     * @param	lon2	east longitude for a region
     * @param	lat2	north latitude for a region
     */
	public static void birthInRegionFilter(List<Typhoon> res,float thre,float lon1,float lat1,float lon2,float lat2){
		for(int l=0;l<res.size();l++){
			Typhoon tr=res.get(l);
			
			if(!tr.birthInRegion(thre,lon1,lat1,lon2,lat2)){
				res.remove(l);
				l--;
			}
		}
	}
	
	
	/**
     * convert record to intensity variables
     *
     * @param	tr	a typhoon record
     * 
     * @return	vs	intensity variables
     */
	public static Variable[] toIntensityVariables(Typhoon tr){
		Variable pres=null;
		Variable wnds=null;
		
		float[] winds=tr.getWinds();
		float[] press=tr.getPressures();
		
		if(press!=null){
			pres=new Variable("prs"+tr.getName(),false,new Range(tr.getTCount(),1,1,1));
			pres.setCommentAndUnit("time series of "+tr.getName()+"'s minimum central pressure");
			System.arraycopy(press,0,pres.getData()[0][0][0],0,tr.getTCount());
			pres.setUndef(-9999);
		}
		
		if(winds!=null){
			wnds=new Variable("wnd"+tr.getName(),false,new Range(tr.getTCount(),1,1,1));
			wnds.setCommentAndUnit("time series of "+tr.getName()+"'s maximum sustained wind speed");
			System.arraycopy(winds,0,wnds.getData()[0][0][0],0,tr.getTCount());
			wnds.setUndef(-9999);
		}
		
		int count=0;
		if(pres!=null) count++;
		if(wnds!=null) count++;
		
		if(count==2){
			return new Variable[]{pres,wnds};
			
		}else if(count==1){
			if(pres!=null) return new Variable[]{pres};
			else return new Variable[]{wnds};
			
		}else throw new IllegalArgumentException("no intensity data");
	}
	
	
	/*** helper methods ***/
	private static void validate(List<Typhoon> res,Predicate<Record> cond){
		List<Typhoon> remove=new ArrayList<>();
		
		for(Typhoon ty:res){
			int count=ty.getTCount();
			
			float[] lon=ty.getXPositions();
			
			for(int i=1;i<count;i++) if(lon[i-1]-lon[i]>300) lon[i]+=360;
			
			boolean removeAll=false;
			
			if(removeIdenticalTimeRecords(ty)==count) removeAll=true;
			
			Predicate<Record> cond1=r->{
				int hour=(int)(r.getTime()%1000000L/10000L);
				return hour!=0&&hour!=6&&hour!=12&&hour!=18;
			};
			
			if(removeRecords(ty,cond1.and(cond))==count) removeAll=true;
			
			if(removeAll){
				remove.add(ty);
				System.out.println("remove the Typhoon ("+ty.getID()+") with no valid record");
			}
		}
		
		res.removeAll(remove);
		
		for(Typhoon tr:res) tr.cVelocityByPosition();
	}
	
	private static int removeRecords(Typhoon ty,Predicate<Record> cond){
		// stateless filtering
		int count =ty.getTCount();
		int remain=ty.removeRecords(cond).getTCount();
		return count-remain;
	}
	
	private static int removeIdenticalTimeRecords(Typhoon ty){
		// stateful filtering
		int count=ty.getTCount();
		long last=ty.getTime(0);
		
		boolean[] del=new boolean[count];
		
		for(int i=1;i<count;i++){
			long current=ty.getTime(i);
			
			if(current==last) del[i]=true;
			
			last=current;
		}
		
		int remain=ty.removeRecords(del).getTCount();
		
		return count-remain;
	}
	
	
	
	/** test
	public static void main(String[] args){
		AccessBestTrack.glueJTWCSourceFiles(JTWCBasins.IO);
		AccessBestTrack.glueJTWCSourceFiles(JTWCBasins.WP);
		AccessBestTrack.glueJTWCSourceFiles(JTWCBasins.SH);
	}*/
}
