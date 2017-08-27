/**
 * @(#)TCTracker.java	1.0 2014.12.31
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

import java.io.FileWriter;
import java.io.IOException;

import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import miniufo.lagrangian.Typhoon;


/**
 * TC tracker is used to find TC track and intensity in gridded data
 *
 * @version 1.0, 2014.12.31
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TCTracker{
	//
	private int tcount=0;
	
	private int[] itag=null;
	private int[] jtag=null;
	
	private float[] lons=null;
	private float[] lats=null;
	private float[] minp=null;
	
	private DataDescriptor dd=null;
	
	private Region2D region=null;
	
	
	/**
	 * constructor
	 */
	public TCTracker(DataDescriptor dd){
		this.dd=dd;
		this.region=dd.toRegion2D();
		
		tcount=dd.getTCount();
		
		itag=new int[tcount];
		jtag=new int[tcount];
		
		lons=new float[tcount];
		lats=new float[tcount];
		minp=new float[tcount];
	}
	
	
	/**
	 * Find track according to the minimum of sea-level pressure (SLP).
	 * Call this method first before calling other methods.
	 * 
	 * @param	vname	variable name, usually is SLP
	 * @param	slon	start longitude
	 * @param	slat	start latitude
	 * @param	sRad	search radius
	 */
	public void findTrackBySLP(String vname,float slon,float slat,int sRad){
		int x=dd.getXCount(),y=dd.getYCount();
		
		Variable var=new Variable(vname,false,new Range(tcount,1,y,x));
		
		DataRead dr=DataIOFactory.getDataRead(dd);
		dr.readData(var);dr.closeFile();
		
		float[][][][] vdata=var.getData();
		
		if(!region.inRange(slon,slat))
		throw new IllegalArgumentException("given location ("+slon+","+slat+") is not in the region "+region);
		
		int cX=dd.getXNum(slon);
		int cY=dd.getYNum(slat);
		
		for(int l=0;l<tcount;l++){
			minp[l]=Float.MAX_VALUE;
			
			for(int j=cY-sRad,J=cY+sRad;j<=J;j++)
			for(int i=cX-sRad,I=cX+sRad;i<=I;i++){
				if(i<0  ||j<0  ) continue;
				if(i>x-1||i>y-1) continue;
				
				if(vdata[0][j][i][l]<minp[l]){
					minp[l]=vdata[0][j][i][l];
					itag[l]=i;
					jtag[l]=j;
				}
			}
			
			cX=itag[l];	lons[l]=dd.getXDef().getSamples()[itag[l]];
			cY=jtag[l];	lats[l]=dd.getYDef().getSamples()[jtag[l]];
		}
	}
	
	public void findTrackBySLP(String vname,Typhoon ty,int sRad){
		int x=dd.getXCount(),y=dd.getYCount();
		
		if(tcount!=ty.getTCount()) throw new IllegalArgumentException("T-lengths not equal");
		
		Variable var=new Variable(vname,false,new Range(tcount,1,y,x));
		
		DataRead dr=DataIOFactory.getDataRead(dd);
		dr.readData(var);dr.closeFile();
		
		float[][][] vdata=var.getData()[0];
		
		for(int l=0;l<tcount;l++){
			minp[l]=Float.MAX_VALUE;
			
			float lon=ty.getXPosition(l);
			float lat=ty.getYPosition(l);
			
			if(!region.inRange(lon,lat))
			throw new IllegalArgumentException("given location ("+lon+","+lat+") is not in the region "+region);
			
			int cX=dd.getXNum(lon);
			int cY=dd.getYNum(lat);
			
			for(int j=cY-sRad,J=cY+sRad;j<=J;j++){ if(j<0||j>y-1) continue;
			for(int i=cX-sRad,I=cX+sRad;i<=I;i++){ if(i<0||i>x-1) continue;
				if(vdata[j][i][l]<minp[l]){
					minp[l]=vdata[j][i][l];
					itag[l]=i;
					jtag[l]=j;
				}
			}}
			
			if(minp[l]==Float.MAX_VALUE) throw new IllegalArgumentException("no valid comparison");
			
			lons[l]=dd.getXDef().getSamples()[itag[l]];
			lats[l]=dd.getYDef().getSamples()[jtag[l]];
		}
	}
	
	
	public Variable getDataAlongTrackAsVariable(String vname){
		Variable var=new Variable(vname,false,new Range(tcount,1,dd.getYCount(),dd.getXCount()));
		Variable res=new Variable(vname,false,new Range(tcount,1,1,1));
		
		DataRead dr=DataIOFactory.getDataRead(dd);
		dr.readData(var);dr.closeFile();
		
		float[] rdata=res.getData()[0][0][0];
		float[][][] vdata=var.getData()[0];
		
		for(int l=0;l<tcount;l++) rdata[l]=vdata[jtag[l]][itag[l]][l];
		
		return res;
	}
	
	
	/*** getor and setor ***/
	public int getTCount(){ return tcount;}
	
	public float[] getLongitudes(){ return lons;}
	
	public float[] getLatitudes(){ return lats;}
	
	public float[] getMinimumSLPs(){ return minp;}
	
	
	/**
	 * to track txt file which can be read by tctrack.gs
	 */
	public void toTrackFile(String path){
		StringBuilder sb=new StringBuilder();
		
		sb.append("*\t"+tcount+"\tnameless");
		for(int i=0;i<tcount;i++){
			sb.append(lons[i]+"\t");
			sb.append(lats[i]+"\t");
			sb.append(minp[i]+"\t\n");
		}
		
		try(FileWriter fw=new FileWriter(path)){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String argv[]){
		DataDescriptor dd=DiagnosisFactory.getDataDescriptor("d:/interim.ctl");
		
		TCTracer tt=new TCTracer(dd);
		tt.findTrack("slp",120,24,5);
		
		Variable slp=tt.getDataAlongTrack("slp");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream("d:/slp.dat");
		cdws.writeData(dd,slp);	cdws.closeFile();
		
		tt.toTrackFile("d:/track.txt");
	}*/
}
