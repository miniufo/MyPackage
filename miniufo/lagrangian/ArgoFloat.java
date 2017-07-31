/**
 * @(#)ArgoFloat.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.ArrayList;
import miniufo.diagnosis.MDate;
import miniufo.mathsphysics.Spline;
import static java.lang.Math.cos;
import static java.lang.Math.toRadians;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.lagrangian.Record.undef;


/**
 * a class representing an Argo float observation
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ArgoFloat extends Particle{
	//
	private static final long serialVersionUID = -5968083413710694473L;
	
	
	/**
	 * constructor
	 * 
	 * @param	id	identifier of this ArgoFloat
	 */
	public ArgoFloat(String id,int attLen){ super(id,attLen);}
	
	public ArgoFloat(String id,int size,int attLen){ super(id,size,attLen);}
	
	
	/**
	 * compute daily current speed using lon/lat coordinates
	 */
	public void cDailyCurrentSpeed(){
		int len=records.size();
		
		if(len<2) return;
		
		float[] uc=new float[len];
		float[] vc=new float[len];
		
		for(int i=1,I=len-1;i<I;i++){
			Record r1=records.get(i-1);
			Record r2=records.get(i+1);
			
			float lon1=r1.getLon(),lon2=r2.getLon();
			float lat1=r1.getLat(),lat2=r2.getLat();
			
			uc[i]=(float)toRadians(lon2-lon1)*(EARTH_RADIUS*(float)cos(toRadians((lat1+lat2)/2)))/(86400f*2);
			vc[i]=(float)toRadians(lat2-lat1)*(EARTH_RADIUS)/(86400f*2);
		}
		
		Record r1=records.get(0);
		Record r2=records.get(1);
		float lon1=r1.getLon(),lon2=r2.getLon();
		float lat1=r1.getLat(),lat2=r2.getLat();
		uc[0]=(float)toRadians(lon2-lon1)*(EARTH_RADIUS*(float)cos(toRadians((lat1+lat2)/2)))/86400f;
		vc[0]=(float)toRadians(lat2-lat1)*(EARTH_RADIUS)/86400f;
		
		r1=records.get(len-2);
		r2=records.get(len-1);
		lon1=r1.getLon();lon2=r2.getLon();
		lat1=r1.getLat();lat2=r2.getLat();
		uc[len-1]=(float)toRadians(lon2-lon1)*(EARTH_RADIUS*(float)cos(toRadians((lat1+lat2)/2)))/86400f;
		vc[len-1]=(float)toRadians(lat2-lat1)*(EARTH_RADIUS)/86400f;
		
		for(int i=0;i<len;i++){
			Record p=records.get(i);
			
			p.setData(0,uc[i]);
			p.setData(1,vc[i]);
		}
	}
	
	/**
	 * Interpolate the daily-mean data using spline
	 */
	public void interpolateDailyPosition(){
		if(records.size()<2){
			System.out.println(String.format("%10s has less than 2 records",id));
			return;
		}
		
		int[] tags=getLonLatValidTag();
		
		float[]   sx=tagToFloat(tags);
		float[][] sy=getLonLatValidData(tags);
		
		if(sx.length!=2){
			Spline slon=new Spline(sx,sy[0]);
			Spline slat=new Spline(sx,sy[1]);
			
			slon.cubicSplineWithNaturalBC();
			slat.cubicSplineWithNaturalBC();
			
			int size=records.size();
			float[] x=new float[size];
			float[] y=new float[size];
			
			for(int i=0;i<size;i++) x[i]=i;
			
			slon.cValues(x,y);
			for(int i=0;i<size;i++) records.get(i).setLon(y[i]);
			
			slat.cValues(x,y);
			for(int i=0;i<size;i++) records.get(i).setLat(y[i]);
		}
	}
	
	/**
	 * truncate the record that are out of t-range [str, end]
	 */
	public void truncate(long strtime,long endtime){
		ArrayList<Record> removes=new ArrayList<Record>();
		
		for(Record r:records){
			long t=r.getTime();
			
			if(t<strtime||t>endtime) removes.add(r);
		}
		
		records.removeAll(removes);
	}
	
	/**
	 * if float cross the international date line,
	 * its longitude will be added 360
	 */
	public void crossIDLToContinuousRecord(){
		boolean cross=false;
		
		if(records.size()<2) return;
		
		float lon1=records.get(0).getLon();
		
		if(lon1==undef) throw new IllegalArgumentException("first lon should not be undef");
		
		for(int i=1,I=records.size();i<I;i++){
			float lon2=records.get(i).getLon();
			
			if(lon2!=undef){
				if(Math.abs(lon2-lon1)>300){
					cross=true;
					break;
					
				}else{ lon1=lon2;}
			}
		}
		
		if(cross)
		for(Record r:records){
			float lon=r.getLon();
			if(lon!=undef&&lon<180) r.setLon(lon+360);
		}
	}
	
	/**
	 * prevent the longitude from being larger than 360
	 */
	public void crossIDLToDiscontinuousRecord(){
		for(Record r:records){
			float lon=r.getLon();
			if(lon!=undef&&lon>=360) r.setLon(lon-360);
		}
	}
	
	
	/**
     * whether the float's position is continous
     */
	public boolean isContinousPosition(){
		int[] tags=getLonLatValidTag();
		
		float[]  tag=tagToFloat(tags);
		float[][] ll=getLonLatValidData(tags);
		
		for(int i=0,I=tag.length-1;i<I;i++){
			float days=tag[i+1]-tag[i];
			float dlon=ll[0][i+1]-ll[0][i];
			float dlat=ll[1][i+1]-ll[1][i];
			
			if(Math.abs(dlon/days)>=1.4){
				System.out.println(String.format(
					"%10s moving too fast along longitude (lon[%4s]:%7.3f, dlon:%6.2f)",id,(int)tag[i+1],ll[0][i+1],dlon
				));
				return false;
			}
			
			if(Math.abs(dlat/days)>=1.4f){
				System.out.println(String.format(
					"%10s moving too fast along latitude  (lat[%4s]:%7.3f, dlat:%6.2f)",id,tag[i+1],ll[1][i+1],dlat
				));
				return false;
			}
			
			if(days>100){
				System.out.println(String.format(
					"%10s has no signal for more than 100 days (tag:%4s)",id,(int)tag[i+1]
				));
				return false;
			}
		}
		
		return true;
	}
	
	
	/**
	 * convert the float data to daily-mean data filled by undef
	 */
	public ArgoFloat toDailyData(){
		ArgoFloat re=new ArgoFloat(id,attachedVars.length);
		
		if(records.size()<1){
			System.out.println(String.format("%10s has no record",id));
			return null;
		}
		
		int[] tag=getDiffDayTag();
		
		for(int i=0,I=tag.length-1;i<I;i++){
			int  count=0;
			float olon=0,olat=0;
			
			for(int l=tag[i],L=tag[i+1];l<L;l++){
				Record r=records.get(l);
				
				float lon=r.getLon();
				float lat=r.getLat();
				
				if(lon!=undef&&lat!=undef){
					olon+=lon;
					olat+=lat;
					count++;
				}
			}
			
			if(count!=0){
				olon/=count;
				olat/=count;
				
			}else{
				olon=undef;
				olat=undef;
			}
			
			if(re.getTCount()>0){
				long curr=records.get(tag[i]).getTime();
				long last=re.getRecord(re.getTCount()-1).getTime();
				
				if(last>=curr) throw new IllegalArgumentException("sort profiles first");
				
				curr=curr-curr%1000000L;
				
				MDate tmp=new MDate(last);
				
				long str=System.currentTimeMillis();
				
				while(true){
					tmp=tmp.addDays(1);
					last=tmp.getLongTime();
					
					if(last!=curr){
						re.addRecord(new Record(last,undef,undef));
						
					}else{
						re.addRecord(new Record(curr,olon,olat));
						break;
					}
					
					if(System.currentTimeMillis()-str>10000L){
						System.out.println(String.format("%10s loops too long",id));
						return null;
					}
				}
				
			}else{
				long curr=records.get(tag[i]).getTime();
				
				curr=curr-curr%1000000L;
				
				re.addRecord(new Record(curr,olon,olat));
			}
		}
		
		return re;
	}
	
	
	/**
	 * split the float data into segments of a given size
	 */
	public ArgoFloat[] split(int size){
		if(size<1) throw new IllegalArgumentException("size should be positive");
		
		int len=records.size()/size;
		int lef=records.size()%size;
		
		if(lef!=0) len++;
		
		ArgoFloat[] re=new ArgoFloat[len];
		
		int ptr=0;
		
		for(int i=0,I=(lef!=0?len-1:len);i<I;i++){
			re[i]=new ArgoFloat(id,size);
			
			for(int j=0;j<size;j++)
			re[i].addRecord(records.get(ptr++));
		}
		
		if(lef!=0){
			re[len-1]=new ArgoFloat(id,lef);
			
			for(int j=0,J=lef;j<J;j++)
			re[len-1].addRecord(records.get(ptr++));
		}
		
		return re;
	}
	
	
	/*** helper methods ***/
	private float[] tagToFloat(int[] tags){
		float[] sx=new float[tags.length];
		
		for(int i=0,I=tags.length;i<I;i++)
		sx[i]=tags[i];
		
		return sx;
	}
	
	private int[] getLonLatValidTag(){
		int count=0;
		
		for(Record r:records){ if(r.getLon()!=undef) count++;}
		
		int[] tags=new int[count];
		
		for(int i=0,I=records.size(),tag=0;i<I;i++){
			Record r=records.get(i);
			
			if(r.getLon()!=undef) tags[tag++]=i;
		}
		
		return tags;
	}
	
	private float[][] getLonLatValidData(int[] tags){
		int count=tags.length;
		
		float[][] lonlat=new float[2][count];
		
		for(int i=0;i<count;i++){
			Record r=records.get(tags[i]);
			
			lonlat[0][i]=r.getLon();
			lonlat[1][i]=r.getLat();
		}
		
		return lonlat;
	}
	
	private int[] getDiffDayTag(){
		ArrayList<Integer> al=new ArrayList<Integer>();
		
		al.add(0);
		
		long last=records.get(0).getTime();
		last=last-last%1000000L;
		
		for(int i=1,I=records.size();i<I;i++){
			long curr=records.get(i).getTime();
			curr=curr-curr%1000000L;
			
			if(curr-last>=1000000L) al.add(i);
			
			last=curr;
		}
		
		al.add(records.size());
		
		int[] re=new int[al.size()];
		
		for(int i=0,I=al.size();i<I;i++) re[i]=al.get(i);
		
		return re;
	}
	
	
	/** test
	public static void main(String[] args){
		ArrayList<ArgoFloat> afs=new ArrayList<ArgoFloat>();
		AccessArgoNC.parseBasicInfo(afs,"D:/Data/Argo/Traj/aoml/19018_traj.nc");
		
		for(ArgoFloat af:afs){
			af.sort();
			af.crossIDLToContinuousRecord();
			
			ArgoFloat af2=af.toDailyData();
			if(af2!=null&&af2.isContinousPosition()){
				af2.interpolateDailyPosition();
				af2.cDailyCurrentSpeed();
				af2.crossIDLToDiscontinuousRecord();
				
				System.out.println(af2);
				ArgoFloat[] re=af2.split(10);
				
				for(ArgoFloat afre:re) System.out.println(afre+"\n");
			}
		}
	}*/
}
