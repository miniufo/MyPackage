/**
 * @(#)TyphoonRecord.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.ArrayList;
import java.util.List;

import miniufo.database.DataBaseUtil;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.SpatialModel;


/**
 * a class associated with the typhoon record
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Typhoon extends Particle{
	//
	private static final long serialVersionUID = 8912012454032659510L;
	
	private String name;
	
	public static enum TYPE{TD,TS,TY,EC,OTHERS};
	
	
	/*** Constructors ***/
	public Typhoon(String id,int attLen){ super(id,attLen);}
	
	public Typhoon(String id,int size,int attLen){ super(id,size,attLen);}
	
	public Typhoon(String id,String name,int size,int attLen){ super(id,size,attLen); this.name=name;}
	
	public Typhoon(String id,String name,List<Record> recs){
		super(id,recs.get(0).getDataLength());
		
		this.name=name;
		
		records.addAll(recs);
	}
	
	
	/*** getor and setor ***/
	public int getTag(long time){ // start from 0
		for(int l=0,L=records.size();l<L;l++)
		if(time==records.get(l).getTime()) return l;
		
		throw new IllegalArgumentException(
			"out of T-Range\n"+time+" is not in ["+
			records.get(0).getTime()+", "+records.get(records.size()-1).getTime()+"]"
		);
	}
	
	public float getACE(){
		float ace=0;	// accumulated cyclone energy [(m/s)^2]
		
		float[] wnds=getAttachedData(0);
		
		for(int l=0,L=records.size();l<L;l++) ace+=wnds[l]*wnds[l];
		
		return ace;
	}
	
	public float getPDI(){
		float pdi=0;	// power dissipation index [(m/s)^3]
		
		float[] wnds=getAttachedData(0);
		
		for(int l=0,L=records.size();l<L;l++) pdi+=wnds[l]*wnds[l]*wnds[l];
		
		return pdi;
	}
	
	public String getName(){ return name;}
	
	public String getTRange(){
		MDate md=new MDate(records.get(0).getTime());
		
		String yr1=String.valueOf(md.getYear());
		String mo1=String.format("%1$02d",md.getMonth());
		String dy1=String.format("%1$02d",md.getDate());
		String hr1=String.format("%1$02d",md.getHour());
		
		md=new MDate(records.get(records.size()-1).getTime());
		
		String yr2=String.valueOf(md.getYear());
		String mo2=String.format("%1$02d",md.getMonth());
		String dy2=String.format("%1$02d",md.getDate());
		String hr2=String.format("%1$02d",md.getHour());
		
		return "time("+yr1+"."+mo1+"."+dy1+"."+hr1+"," +
					   yr2+"."+mo2+"."+dy2+"."+hr2+")";
	}
	
	public long getBirthDateByWind(float windthreshold){
		for(int l=0,L=records.size();l<L;l++){
			Record r=records.get(l);
			
			if(r.getDataValue(2)>=windthreshold) return r.getTime();
		}
		
		throw new IllegalArgumentException("Typhoon never get stronger than: "+windthreshold+" m/s");
	}
	
	public  TYPE[] getTypes(){
		int len=records.size();
		
		TYPE[] types=new TYPE[len];
		
		for(int l=0;l<len;l++)
		types[l]=ordinalToType(Math.round(records.get(l).getDataValue(4)));
		
		return types;
	}
	
	public float[] getTranslatingSpeeds(){
		int len=records.size();
		
		float[] re=new float[len];
		
		float[] ux=getZonalVelocity();
		float[] vy=getMeridionalVelocity();
		
		for(int l=0;l<len;l++) re[l]=(float)Math.hypot(ux[l],vy[l]);
		
		return re;
	}
	
	public float[] getSpeeds(){
		int len=records.size();
		
		float[] lons=getLongitudes();
		float[] lats=getLatitudes();
		
		if(len>1){
			float dt=getDT(records.get(0).getTime(),records.get(1).getTime());
			
			float[] spds=new float[len];
			
			spds[0]=
			SpatialModel.cSphericalDistanceByDegree(lons[0],lats[0],lons[1],lats[1])/dt;
			
			for(int l=1,L=len-1;l<L;l++) spds[l]=
			SpatialModel.cSphericalDistanceByDegree(lons[l-1],lats[l-1],lons[l+1],lats[l+1])/(2*dt);
			
			spds[len-1]=
			SpatialModel.cSphericalDistanceByDegree(lons[len-2],lats[len-2],lons[len-1],lats[len-1])/dt;
			
			return spds;
			
		}else return new float[len];
	}
	
	public float[] getStormRelativeWinds(){
		int len=records.size();
		
		float[] spds=getSpeeds();
		float[] wnds=getAttachedData(0);
		
		for(int l=0;l<len;l++) wnds[l]-=spds[l];
		
		return wnds;
	}
	
	public float[] getWinds(){
		if(attachedVars.length<=2) return null;
		else return getAttachedData(2);
	}
	
	public float[] getPressures(){
		if(attachedVars.length<=3) return null;
		else return getAttachedData(3);
	}
	
	
	/**
	 * compute changes of the data (delta-data) using centered time differencing.
	 * The first and last are computed using one-sided time differencing.
	 */
	public static float[] getChangesByCentralDiff(float[] data){
		int N=data.length;
		
		if(N==1) return new float[1];
		
		float[] ch=new float[N];
		
		ch[0]=data[1]-data[0];
		
		for(int l=1,L=N-1;l<L;l++) ch[l]=data[l+1]-data[l-1];
		
		ch[N-1]=data[N-1]-data[N-2];
		
		return ch;
	}
	
	public static float[] getChangesByForwardDiff(float[] data,int interv){
		int N=data.length;
		
		float[] ch=new float[N];
		
		for(int l=0;l<N;l++) ch[l]=DataBaseUtil.undef;
		
		for(int l=0,L=N-interv;l<L;l++) ch[l]=data[l+interv]-data[l];
		
		return ch;
	}
	
	
	/**
     * whether the TC is birth (wind > threshold) within a specific region
     *
     * @param	tr		typhoon record
     * @param	thre	threshold of the wind larger than which a TC is said to give birth
     * @param	lon1	west longitude for a region
     * @param	lat1	south latitude for a region
     * @param	lon2	east longitude for a region
     * @param	lat2	north latitude for a region
     */
	public boolean birthInRegion(float thre,float lon1,float lat1,float lon2,float lat2){
		if(lon2<lon1) throw new IllegalArgumentException("invalid region: "+lon1+" > "+lon2);
		if(lat2<lat1) throw new IllegalArgumentException("invalid region: "+lat1+" > "+lat2);
		
		long time=getBirthDateByWind(thre);
		
		Record r=records.get(getTag(time));
		
		float lon=r.getLon();
		float lat=r.getLat();
		
		if(lon>=lon1&&lon<=lon2&&lat>=lat1&&lat<=lat2) return true;
		
		return false;
	}
	
	/**
     * insert records using linear interpolation
     * 
     * @param	insertNum	 number of records that inserted between two time slices
     */
	public Typhoon interpolateAlongT(int insertNum){
		if(insertNum<1)
		throw new IllegalArgumentException("insert number should be larger than 0");
		
		if(records.size()<2)
		throw new IllegalArgumentException("no enough records for insertion");
		
		MDate str=new MDate(records.get(0).getTime());
		
		int delim=insertNum+1;
		int dt =str.getDT(new MDate(records.get(1).getTime()));
		int ndt=dt/delim;
		int dlen=records.get(0).getDataLength();
		
		List<Record> ls=new ArrayList<>(records.size()+insertNum*(records.size()-1));
		
		for(int l=0,L=records.size()-1;l<L;l++){
			MDate tmp=str.addSeconds(ndt);
			
			Record rstr=records.get(l);
			Record rend=records.get(l+1);
			
			ls.add(rstr);
			
			float lonstr=rstr.getLon();	float lonend=rend.getLon();
			float latstr=rstr.getLat();	float latend=rend.getLat();
			
			float uvelstr=rstr.getDataValue(0);	float uvelend=rend.getDataValue(0);
			float vvelstr=rstr.getDataValue(1);	float vvelend=rend.getDataValue(1);
			
			float windstr=0;	float windend=0;
			float presstr=0;	float presend=0;
			
			if(dlen>2){ windstr=rstr.getDataValue(2); windend=rend.getDataValue(2);}
			if(dlen>3){ presstr=rstr.getDataValue(3); presend=rend.getDataValue(3);}
			
			for(int i=1;i<=insertNum;i++){
				long time=tmp.getLongTime();
				
				float lon =((delim-i)*lonstr+i*lonend)/delim;
				float lat =((delim-i)*latstr+i*latend)/delim;
				float uvel=((delim-i)*uvelstr+i*uvelend)/delim;
				float vvel=((delim-i)*vvelstr+i*vvelend)/delim;
				float wind=((delim-i)*windstr+i*windend)/delim;
				float pres=((delim-i)*presstr+i*presend)/delim;
				
				Record r=new Record(time,lon,lat,dlen);
				
				r.setData(0,uvel);
				r.setData(1,vvel);
				
				if(dlen>2) r.setData(2,wind);
				if(dlen>3) r.setData(3,pres);
				
				ls.add(r);
				
				tmp=tmp.addSeconds(ndt);
			}
			
			str=str.addSeconds(dt);
		}
		
		ls.add(records.get(records.size()-1));
		records=ls;
		
		return this;
	}
	
	/**
     * generate a CSM-type string
     * 
     * @param	xcount	grid count in azimuthal direction
     * @param	ycount	grid count in radial direction
     * @param	deltaY	grid interval in radial direction (degree)
     */
	public String toCSMString(int xcount,int ycount,int zcount,float deltaY,float deltaZ){
		return toCSMString("d:/ctl.ctl",xcount,ycount,zcount,deltaY,deltaZ,1000);
	}
	
	public String toCSMString(String ctlpath,int xcount,int ycount,int zcount,float deltaY,float deltaZ,float zstr){
		int len=records.size();
		
		long[] date=getTimes();
		float[] lons=getLongitudes();
		float[] lats=getLatitudes();
		float[] pres=getAttachedData(1);
		float[] wnds=getAttachedData(0);
		
		StringBuilder sb=new StringBuilder();
		
		sb.append("dpath "+ctlpath+"\n");
		sb.append("title "+name+" ("+id+") \n");
		sb.append("xdef "+xcount+"\n");
		sb.append("ydef "+ycount+" "+deltaY+"\n");
		sb.append("zdef "+zcount+" "+zstr+" "+deltaZ+"\n");
		sb.append(
			"tdef "+len+
			" "+new MDate(date[0]).toGradsDate()+
			" "+String.format("%.0f",(getDT(date[0],date[1])/3600f))+
			"hr\n"
		);
		sb.append("coords\n");
		for(int l=0;l<len;l++)
		sb.append(
			lons[l]+" "+lats[l]+" "+
			(pres==null?"":pres[l])+" "+
			(wnds==null?"":wnds[l])+"\n"
		);
		sb.append("endcoords\n");
		
		return sb.toString();
	}
	
	/**
     * used to print out
     */
	public String toString(){
		int len=records.get(0).getDataLength();
		
		float[] wnds=null;
		float[] pres=null;
		TYPE[]  type=null;
		
		if(len==5){
			wnds=getAttachedData(2);
			pres=getAttachedData(3);
			type=getTypes();
			
		}else if(len==4){
			wnds=getAttachedData(2);
			pres=getAttachedData(3);
			
		}else if(len==3){
			wnds=getAttachedData(2);
		}
		
		StringBuilder sb=new StringBuilder();
		
		sb.append(String.format("* %3d  %s  (%s)\n",getTCount(),name,id));
		
		for(int l=0,L=records.size();l<L;l++){
			Record r=records.get(l);
			
			sb.append(String.format("%5.1f  %4.1f  %14d",r.getLon(),r.getLat(),r.getTime()));
			
			if(len==5)      sb.append(String.format(" %7.3f  %4.0f  %s\n",wnds[l],pres[l],type[l]));
			else if(len==4) sb.append(String.format(" %7.3f  %4.0f\n",wnds[l],pres[l]));
			else if(len==3) sb.append(String.format(" %7.3f\n",wnds[l]));
		}
		
		return sb.toString();
	}
	
	
	/**
     * convert type data from flow to TYPE
     *
     * @param	type	type data in float format
     */
	public static TYPE getType(float type){ return ordinalToType(Math.round(type));}
	
	
	/*** helper methods ***/
	private static TYPE ordinalToType(int ordinal){
		TYPE[] types=TYPE.values();
		
		for(int i=0,I=types.length;i<I;i++)
		if(ordinal==types[i].ordinal()) return types[i];
		
		throw new IllegalArgumentException("no ordinal for number: "+ordinal);
	}
	
	private static int getDT(long t1,long t2){
		return new MDate(t1).getDT(new MDate(t2));
	}
	
	
	/** test
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons("D:/Data/Typhoons/JMA/JMA.txt","time=1Jan2011-31Dec2011",DataSets.JMA);
		
		int sum=0;
		for(Typhoon t:ls){ sum+=t.getTCount();}
		
		Stream<Record> s=Typhoon.toRecordStream(ls);
		
		System.out.println(sum+"\t"+s.count());
	}*/
}