/**
 * @(#)Typhoon.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.ArrayList;
import java.util.List;
import miniufo.diagnosis.MDate;


/**
 * A class of Typhoon representing best-track data
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
	
	public static final AttachedMeta Vmax=new AttachedMeta("Vmax",2);
	public static final AttachedMeta Pmin=new AttachedMeta("Pmin",3);
	public static final AttachedMeta Type=new AttachedMeta("Type",4);
	
	
	/*** Constructors ***/
	public Typhoon(String id,int attLen){ super(id,attLen);}
	
	public Typhoon(String id,int size,int attLen){ super(id,size,attLen,true);}
	
	public Typhoon(String id,String name,int size,int attLen){ super(id,size,attLen,true); this.name=name;}
	
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
		
		float[] wnds=getAttachedData(Vmax);
		
		for(int l=0,L=records.size();l<L;l++) ace+=wnds[l]*wnds[l];
		
		return ace;
	}
	
	public float getPDI(){
		float pdi=0;	// power dissipation index [(m/s)^3]
		
		float[] wnds=getAttachedData(Vmax);
		
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
			
			if(r.getData(Vmax)>=windthreshold) return r.getTime();
		}
		
		throw new IllegalArgumentException("Typhoon never get stronger than: "+windthreshold+" m/s");
	}
	
	public  TYPE[] getTypes(){
		int len=records.size();
		
		TYPE[] types=new TYPE[len];
		
		for(int l=0;l<len;l++)
		types[l]=ordinalToType(Math.round(records.get(l).getData(Type)));
		
		return types;
	}
	
	public float[] getTranslatingSpeeds(){
		int len=records.size();
		
		float[] re=new float[len];
		
		float[] ux=getUVel();
		float[] vy=getVVel();
		
		for(int l=0;l<len;l++) re[l]=(float)Math.hypot(ux[l],vy[l]);
		
		return re;
	}
	
	public float[] getStormRelativeWinds(){
		int len=records.size();
		
		float[] spds=getSpeeds();
		float[] wnds=getAttachedData(Vmax);
		
		for(int l=0;l<len;l++) wnds[l]-=spds[l];
		
		return wnds;
	}
	
	public float[] getWinds(){
		if(meta.length<=2) return null;
		else return getAttachedData(Vmax);
	}
	
	public float[] getPressures(){
		if(meta.length<=3) return null;
		else return getAttachedData(Pmin);
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
		
		float lon=r.getXPos();
		float lat=r.getYPos();
		
		if(lon>=lon1&&lon<=lon2&&lat>=lat1&&lat<=lat2) return true;
		
		return false;
	}
	
	/**
     * Insert records using linear interpolation.
     * 
     * @param	insertNum	 number of records that inserted between two time slices
     */
	public Typhoon interpolateAlongT(int insertNum){
		if(insertNum<1) return this;
		
		if(records.size()<2)
		throw new IllegalArgumentException("no enough records (2 at least) for insertion");
		
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
			
			float lonstr=rstr.getXPos();	float lonend=rend.getXPos();
			float latstr=rstr.getYPos();	float latend=rend.getYPos();
			
			float uvelstr=rstr.getData(UVEL);	float uvelend=rend.getData(UVEL);
			float vvelstr=rstr.getData(VVEL);	float vvelend=rend.getData(VVEL);
			
			float windstr=0;	float windend=0;
			float presstr=0;	float presend=0;
			
			if(dlen>2){ windstr=rstr.getData(Vmax); windend=rend.getData(Vmax);}
			if(dlen>3){ presstr=rstr.getData(Pmin); presend=rend.getData(Pmin);}
			
			for(int i=1;i<=insertNum;i++){
				long time=tmp.getLongTime();
				
				float lon =((delim-i)*lonstr+i*lonend)/delim;
				float lat =((delim-i)*latstr+i*latend)/delim;
				float uvel=((delim-i)*uvelstr+i*uvelend)/delim;
				float vvel=((delim-i)*vvelstr+i*vvelend)/delim;
				float wind=((delim-i)*windstr+i*windend)/delim;
				float pres=((delim-i)*presstr+i*presend)/delim;
				
				Record r=new Record(time,lon,lat,dlen);
				
				r.setData(UVEL,uvel);
				r.setData(VVEL,vvel);
				
				if(dlen>2) r.setData(Vmax,wind);
				if(dlen>3) r.setData(Pmin,pres);
				
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
     * Interpolate so that the time interval between records is deltaT.
     * 
     * @param	deltaT	 a given constant deltaT (second)
     */
	public Typhoon interpolateToDT(int deltaT){
		List<Record> ls=new ArrayList<>();
		
		ls.add(records.get(0));
		
		for(int i=1,I=records.size();i<I;i++){
			Record str=records.get(i-1);
			Record end=records.get(i  );
			
			ls.addAll(Record.interpolateToDT(str,end,deltaT));
			ls.add(end);
		}
		
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
		float[] lons=getXPositions();
		float[] lats=getYPositions();
		float[] pres=getAttachedData(Pmin);
		float[] wnds=getAttachedData(Vmax);
		
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
			wnds=getAttachedData(Vmax);
			pres=getAttachedData(Pmin);
			type=getTypes();
			
		}else if(len==4){
			wnds=getAttachedData(Vmax);
			pres=getAttachedData(Pmin);
			
		}else if(len==3){
			wnds=getAttachedData(Vmax);
		}
		
		StringBuilder sb=new StringBuilder();
		
		sb.append(String.format("* %3d  %s  (%s)\n",getTCount(),name,id));
		
		for(int l=0,L=records.size();l<L;l++){
			Record r=records.get(l);
			
			sb.append(String.format(
				"%6.2f  %5.2f %6.3f  %5.3f  %14d",
				r.getXPos(),r.getYPos(),r.getData(UVEL),r.getData(VVEL),r.getTime()
			));
			
			if(len==5)      sb.append(String.format(" %7.3f  %5.1f  %s\n",wnds[l],pres[l],type[l]));
			else if(len==4) sb.append(String.format(" %7.3f  %5.1f\n",wnds[l],pres[l]));
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
		
		return TYPE.OTHERS;
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
