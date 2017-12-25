/**
 * @(#)Record.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.LongStream;
import miniufo.basic.InterpolationModel;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.SpatialModel;


/**
 * Used to describe a record of a Lagrangian particle at a specific time.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Record implements Serializable{
	//
	private static final long serialVersionUID = 1239922442535904638L;
	
	private int cycleNum=-1;	// designed for Argo
	
	private long  time=0;		// time
	
	private float xpos=0;		// x-position
	private float ypos=0;		// y-position
	
	private float[] data=null;	// values for attached data (e.g., uvel,vvel,temp,salt)
	
	public static final float undef=-9999.0f;	// default undefined value
	
	
	/**
	 * constructor
	 */
	public Record(long time){
		this.time=time;
		
		data=new float[2];
	}
	
	public Record(long time,float xpos,float ypos){
		this.time=time;
		this.xpos=xpos;
		this.ypos=ypos;
		
		data=new float[2];
	}
	
	public Record(long time,float xpos,float ypos,float uvel,float vvel){
		this.time=time;
		this.xpos=xpos;
		this.ypos=ypos;
		
		data=new float[2];
		data[0]=uvel;
		data[1]=vvel;
	}
	
	public Record(long time,float xpos,float ypos,int dataLen){
		this(time,xpos,ypos);
		
		data=new float[dataLen];
		
		for(int i=0;i<dataLen;i++) data[i]=undef;
	}
	
	public Record(Record r){
		this.time=r.time;
		this.xpos=r.xpos;
		this.ypos=r.ypos;
		
		if(r.data!=null) data=r.data.clone();
	}
	
	
	/**
	 * Compute x-distance relative to the position of another record.
	 * 
	 * @param	r		another record
	 * @param	llpos	lat/lon position (degree) or Y/X position (m)
	 */
	public float cXDistanceTo(Record r,boolean llpos){
		if(  xpos==undef||  ypos==undef) return undef;
		if(r.xpos==undef||r.ypos==undef) return undef;
		
		if(llpos)	// in unit of degree (lat/lon grid)
			return (float)(SpatialModel.EARTH_RADIUS*(xpos-r.xpos)*Math.PI/180.0*Math.cos((ypos+r.ypos)/2.0));
		else		// in unit of meter (cartesian grid)
			return Math.abs(xpos-r.xpos);
	}
	
	/**
	 * Compute y-distance relative to the position of another record.
	 * 
	 * @param	r		another record
	 * @param	llpos	lat/lon position (degree) or Y/X position (m)
	 */
	public float cYDistanceTo(Record r,boolean llpos){
		if(  xpos==undef||  ypos==undef) return undef;
		if(r.xpos==undef||r.ypos==undef) return undef;
		
		if(llpos)	// in unit of degree (lat/lon grid)
			return (float)(SpatialModel.EARTH_RADIUS*(ypos-r.ypos)*Math.PI/180.0);
		else		// in unit of meter (cartesian grid)
			return Math.abs(ypos-r.ypos);
	}
	
	/**
	 * Compute total-distance relative to the position of another record.
	 * 
	 * @param	r		another record
	 * @param	llpos	lat/lon position (degree) or Y/X position (m)
	 */
	public float cDistanceTo(Record r,boolean llpos){
		if(  xpos==undef||  ypos==undef) return undef;
		if(r.xpos==undef||r.ypos==undef) return undef;
		
		if(llpos)	// in unit of degree (lat/lon grid)
			return SpatialModel.cSphericalDistanceByDegree(xpos,ypos,r.xpos,r.ypos);
		else		// in unit of meter (cartesian grid)
			return (float)Math.hypot(xpos-r.xpos,ypos-r.ypos);
	}
	
	
	/*** getor and setor ***/
	public int getCycleNum(){ return cycleNum;}
	
	public int getDataLength(){ return data.length;}
	
	public long getTime(){ return time;}
	
	public float getXPos(){ return xpos;}
	
	public float getYPos(){ return ypos;}
	
	public float getDataValue(int idx){ return data[idx];}
	
	public float[] getDataValues(){ return data;}
	
	public void setCycleNum(int cycleNum){ this.cycleNum=cycleNum;}
	
	public void setXPos(float xpos){ this.xpos=xpos;}
	
	public void setYPos(float ypos){ this.ypos=ypos;}
	
	public void setTime(long time){ this.time=time;}
	
	public void setData(int idx,float value){ data[idx]=value;}
	
	public void addData(int idx,float value){ data[idx]+=value;}
	
	
	/**
	 * Interpolate records between start and end.
	 * 
	 * @param	str		start record
	 * @param	end		end   record
	 * @param	times	times for interpolation
	 */
	public static List<Record> interpolateBetween(Record str,Record end,long... times){
		MDate ref=new MDate(str.time);
		
		float dt=ref.getDT(new MDate(end.time));
		
		return LongStream.of(times).mapToObj(time->{
			if(time<=str.time) throw new IllegalArgumentException("time is before str");
			if(time>=end.time) throw new IllegalArgumentException("time is after  end");
			
			float ratio=ref.getDT(new MDate(time))/dt;
			
			Record res=new Record(time,0,0,str.data.length);
			
			for(int i=0,I=str.data.length;i<I;i++)
			res.data[i]=InterpolationModel.linearInterpolation(str.data[i],end.data[i],ratio,undef);
			res.xpos    =InterpolationModel.linearInterpolation(str.xpos    ,end.xpos    ,ratio,undef);
			res.ypos    =InterpolationModel.linearInterpolation(str.ypos    ,end.ypos    ,ratio,undef);
			res.cycleNum=str.cycleNum;
			
			return res;
			
		}).collect(Collectors.toList());
	}
	
	/**
	 * Interpolate records to an interval of deltaT.
	 * 
	 * @param	str		start record
	 * @param	end		end   record
	 * @param	deltaT	interval of deltaT (seconds)
	 */
	public static List<Record> interpolateToDT(Record str,Record end,int deltaT){
		long dt=new MDate(str.time).getDT(new MDate(end.time));
		
		if(dt==deltaT  ) return new ArrayList<>();
		if(dt< deltaT  ) throw new IllegalArgumentException("deltaT should be larger than "+dt+" seconds");
		if(dt%deltaT!=0) throw new IllegalArgumentException("deltaT cannot divide "+dt);
		
		int len=(int)(dt/deltaT)-1;
		
		long[] times=new long[len];
		
		MDate ref=new MDate(str.time);
		
		for(int i=0;i<len;i++) times[i]=ref.addSeconds(deltaT*(i+1)).getLongTime();
		
		return interpolateBetween(str,end,times);
	}
	
	
	/**
     * used to print out
     */
	public String toString(){ return toStringSimple(true);}
	
	public String toString(boolean llpos,boolean wrtAtt){
		if(wrtAtt)
			return toStringFull(llpos);
		else
			return toStringSimple(llpos);
	}
	
	
	/*** helper methods ***/
	private String toStringSimple(boolean llpos){
		return String.format(
			"%9.3f%10.3f%19s%10.3f%10.3f",
			xpos/(llpos?1f:1000f),ypos/(llpos?1f:1000f),time,data[0],data[1]
		);
	}
	
	private String toStringFull(boolean llpos){
		StringBuilder format=new StringBuilder("%9.3f%10.3f%19s");
		
		Object[] os=new Object[data.length+3];
		
		os[0]=xpos/(llpos?1f:1000f);
		os[1]=ypos/(llpos?1f:1000f);
		os[2]=time;
		
		for(int i=0;i<data.length;i++){
			format.append("%10.3f");
			os[i+3]=data[i];
		}
		
		return String.format(format.toString(),os);
	}
	
	
	/*** test **
	public static void main(String[] args){
		Record r1=new Record(1,2,3,4,5);
		Record r2=new Record(r1);
		
		r1.time=3;
		
		System.out.println(r1.time+" "+r2.time);
	}*/
}
