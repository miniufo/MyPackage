/**
 * @(#)KMLDataWriter.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.IO;

import java.io.FileWriter;
import java.io.IOException;

import miniufo.io.FileWriteInterface;
import static java.lang.Math.sqrt;
import static java.lang.Math.atan2;
import static java.lang.Math.toDegrees;
import static miniufo.io.FileWriteInterface.Solution.*;
import static miniufo.diagnosis.SpatialModel.cLatLon;


/**
 * used to convert the data to kml format in Google earth
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class KMLDataWriter{
	//
	private boolean is_skip=true;	// whether to skip overwriting the file
	
	private float scale=0.04f;		// 1 m/s proportional to scale degree in Google earth
	private float arrow=scale*5;	// arrow size
	private float arrag=30f;		// arrow angle
	private float thick=2;			// thickness
	
	private float[] clev=null;
	
	private String name=null;
	private StringBuilder sb=null;
	
	private FileWriter fw=null;
	private FileWriteInterface fwi=null;
	
	// bbggrr
	private static String[] color={
		"c800a0","dc0082",
		"ff3c1e","ffa000",
		"c8c800","8cd200",
		"00dc00","32e6a0",
		"32dce6","2dafe6",
		"2882f0","3c3cfa","8200f0"
	};
	
	
	/**
     * constructor
     *
     * @param	path	path for writing
     */
	public KMLDataWriter(String path){
		fwi=new FileWriteInterface(path);
		this.name=fwi.getFile().getName();
		
		if(fwi.getFlag()!=SKIP){
			is_skip=false;
			
			try{
				switch(fwi.getFlag()){
				case RENAME:
					fw=new FileWriter(fwi.getParent()+fwi.getNewName());
					break;
					
				case OVERWRITE:
					fw=new FileWriter(path);
					break;
					
				case APPEND:
					fw=new FileWriter(path,true);
					break;

				default:
					break;
				}
				
		    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		}
		
		sb=new StringBuilder();	clev=new float[color.length+1];
	}
	
	
	/*** getor and setor ***/
	public void setScale(float f){ scale=f; arrow=scale*5;}
	
	public void setArrow(float f){ arrow=f;}
	
	public void setThick(float f){ thick=f;}
	
	public void setArrowAngle(float f){ arrag=f;}
	
	public void setRBRange(float min,float max){
		if(min>=max) throw new IllegalArgumentException("min should be smaller than max");
		
		float delta=(max-min)/color.length;
		
		for(int i=0,I=clev.length;i<I;i++) clev[i]=min+delta*i;
	}
	
	
	/*** add to buffer ***/
	public void addHeader(){
		sb.append("<?xml version=\"1.0\"?>\n");
		sb.append("<kml xmlns=\"http://earth.google.com/kml/2.0\">\n");
		sb.append("<Document>\n");
		sb.append("<name>"+name+"</name>\n");
		sb.append("<open>0</open>\n");
	}
	
	public void addName(String name){sb.append("<name>"+name+"</name>\n");}
	
	public void addOpen(boolean isopen){
		if(isopen)sb.append("<open>1</open>\n");
		else sb.append("<open>0</open>\n");
	}
	
	public void addFolder(StringBuilder content){
		sb.append("<folder>\n");
		sb.append(content);
		sb.append("</folder>\n");
	}
	
	public void addVector(float u,float v,float lon,float lat){
		float mag=(float)sqrt(u*u+v*v),lamda=(float)toDegrees(atan2(v,u))-90;
		
		String color=getColor(mag);	mag*=scale;
		
		if(lamda<0) lamda+=360;
		
		float[] L1=cLatLon(lon  ,lat  ,lamda,mag  );
		float[] L2=cLatLon(L1[0],L1[1],lamda-180-arrag,arrow);
		float[] L3=cLatLon(L1[0],L1[1],lamda-180+arrag,arrow);
		
		if(lon  >180) lon  -=360;
		if(L1[0]>180) L1[0]-=360;
		if(L2[0]>180) L2[0]-=360;
		if(L3[0]>180) L3[0]-=360;
		
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>ee"+color+"</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LineString>\n");
		sb.append("<coordinates>"+lon  +","+lat  +",0 "+L1[0]+","+L1[1]+",0</coordinates>\n");
		sb.append("</LineString>\n");
		sb.append("</Placemark>\n");
		
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>ee"+color+"</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LineString>\n");
		sb.append("<coordinates>"+L2[0]+","+L2[1]+",0 "+L1[0]+","+L1[1]+",0</coordinates>\n");
		sb.append("</LineString>\n");
		sb.append("</Placemark>\n");
		
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>ee"+color+"</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LineString>\n");
		sb.append("<coordinates>"+L3[0]+","+L3[1]+",0 "+L1[0]+","+L1[1]+",0</coordinates>\n");
		sb.append("</LineString>\n");
		sb.append("</Placemark>\n");
	}
	
	public void addLineString(float[] lon,float[] lat){
		int N=lon.length;
		
		if(lat.length!=N) throw new IllegalArgumentException("lengths not same");
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>eeffffff</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LineString>\n");
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++){
			sb.append(String.format("%.8f",lon[i])+",");
			sb.append(String.format("%.8f",lat[i])+",0\n");
		}
		sb.append("</coordinates>\n");
		sb.append("</LineString>\n");
		sb.append("</Placemark>\n");
	}
	
	public void addLineString(String[] lon,String[] lat){
		int N=lon.length;
		
		if(lat.length!=N) throw new IllegalArgumentException("lengths not same");
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>eeffffff</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LineString>\n");
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++){
			sb.append(lon[i]+",");
			sb.append(lat[i]+",0\n");
		}
		sb.append("</coordinates>\n");
		sb.append("</LineString>\n");
		sb.append("</Placemark>\n");
	}
	
	public void addLineString(String[] lonlat){
		int N=lonlat.length;
		
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>eeffffff</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LineString>\n");
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++)
		sb.append(lonlat[i]+",0\n");
		sb.append("</coordinates>\n");
		sb.append("</LineString>\n");
		sb.append("</Placemark>\n");
	}
	
	public void addLinearRing(float[] lon,float[] lat){
		int N=lon.length;
		
		if(lat.length!=N) throw new IllegalArgumentException("lengths not same");
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>ee000000</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LinearRing>\n");
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++){
			sb.append(String.format("%.8f",lon[i])+",");
			sb.append(String.format("%.8f",lat[i])+",0\n");
		}
		sb.append(String.format("%.8f",lon[0])+",");
		sb.append(String.format("%.8f",lat[0])+",0\n");
		sb.append("</coordinates>\n");
		sb.append("</LinearRing>\n");
		sb.append("</Placemark>\n");
	}
	
	public void addLinearRing(String[] lon,String[] lat){
		int N=lon.length;
		
		if(lat.length!=N) throw new IllegalArgumentException("lengths not same");
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>ee000000</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LinearRing>\n");
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++){
			sb.append(lon[i]+",");
			sb.append(lat[i]+",0\n");
		}
		sb.append(lon[0]+",");
		sb.append(lat[0]+",0\n");
		sb.append("</coordinates>\n");
		sb.append("</LinearRing>\n");
		sb.append("</Placemark>\n");
	}
	
	public void addLinearRing(String[] lonlat){
		int N=lonlat.length;
		
		sb.append("<Placemark>\n");
		sb.append("<Style>\n");
		sb.append("<LineStyle>\n");
		sb.append("<color>ee000000</color>\n");
		sb.append("<width>"+thick+"</width>\n");
		sb.append("</LineStyle>\n");
		sb.append("</Style>\n");
		sb.append("<LinearRing>\n");
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++)
		sb.append(lonlat[i]+",0\n");
		sb.append(lonlat[0]+",0\n");
		sb.append("</coordinates>\n");
		sb.append("</LinearRing>\n");
		sb.append("</Placemark>\n");
	}
	
	public void addTail(){
		sb.append("</Document>\n");
		sb.append("</kml>\n");
	}
	
	
	/*** write data ***/
	public void writeData(){
		try{
			if(!is_skip) fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public void closeFile(){
		is_skip=true;
		
		try{ if(fw!=null) fw.close();}
	    catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		fwi=null;	fw=null;
	}
	
	
	private String getColor(float mag){
		if(mag<clev[1]) return color[0];
		
		for(int i=2,I=clev.length-1;i<I;i++)
		if(mag<=clev[i]) return color[i-1];
		
		return color[color.length-1];
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			StringBuilder a=new StringBuilder("aaa");
			StringBuilder b=new StringBuilder("bbb");
			
			a.append(a);
			System.out.println(a);
			a.append("ccc");
			System.out.println(b);
			System.out.println(a);
			
		}catch(Exception e){ e.printStackTrace();}
	}*/
}
