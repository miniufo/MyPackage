/**
 * @(#)IntermediateFile.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.IO;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;
import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * Intermediate file format for WRF or MM5
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class IntermediateFile{
	//
	private int IFV  =5;	// intermediate-format version number for WPS (5), WRFSI (4), MM5 (3)
	private int NX   =0;	// Slab dimension in the X direction
	private int NY   =0;	// Slab dimension in the Y direction
	private int IPROJ=0;	// Flag denoting the projection
	
	private boolean IS_WIND_EARTH_REL=false;	// whether Lambert projected data has earth or model rotated winds
	
	private float DX   =0;	// Grid-spacing in x (km at TRUELAT1 (and TRUELAT2 as appropriate))
	private float DY   =0;	// Grid-spacing in y (km at TRUELAT1 (and TRUELAT2 as appropriate))
	private float XLVL =0;	// Pressure-level (Pa) of the data.
							// 200100 Pa indicates surface data; 201300 Pa indicates sea-level pressure
	private float NLATS=0;	// number of latitudes north of equator (for Gaussian grids)
	private float XFCST=0;	// The time, in format "YYYY-MM-DD_HH:mm:ss"
	private float XLONC=0;	// Center longitude of the projection
	private float TRUELAT1=0;	// Extra latitude (degrees north) used for defining Mercator,
								// Polar Stereographic, and Lambert conformal projections
	private float TRUELAT2=0;	// A second extra latitude (degrees north)
								// used for defining Lambert conformal projection
	private float STARTLAT=0;	// Starting latitude (degrees north)
	private float STARTLON=0;	// Starting longitude (degrees east)
	private float DELTALAT=0;	// Latitude increment (degrees) for lat/lon grid
	private float DELTALON=0;	// Longitude increment (degrees) for lat/lon grid
	private float EARTH_RADIUS=6367.47f;	// Radius of the earth (km)
	
	private long[] positions=null;	// start positions of slabs
	
	private byte[] HDATE     =new byte[24];	// The time, in format "YYYY-MM-DD_HH:mm:ss"
	private byte[] MAP_SOURCE=new byte[32];	// Source of data
	private byte[] STARTLOC  =new byte[8 ];	// Start location of data, "CENTER" or "SWCORNER". "SWCORNER" is typical
	private byte[] FIELD     =new byte[9 ]; // A field name.  Names with special meaning are described below
	private byte[] UNITS     =new byte[25];	// Units describing the field in the slab.
	private byte[] DESC      =new byte[46];	// Text description of the field in the slab.
	
	private String fname=null;		// full path and file name
	
	private StringBuilder info=null;// information of headers
	private StringBuilder dump=null;// information for dump
	
	private Format fmt=null;		// format of the file
	
	public enum Format{WPS, MM5, SI};
	
	
	/**
	 * constructor
	 *
     * @param	fname	complete file name
     */
 	public IntermediateFile(String fname){
		this.fname=fname;
		
		File f=new File(fname);
		
		if(!f.exists()) throw new IllegalArgumentException("File do not exist: "+fname);
		
		long flength=new File(fname).length();
		
		info=new StringBuilder();
		dump=new StringBuilder();
		
		List<Long> ls=new ArrayList<Long>();
		
		try{
			RandomAccessFile raf=new RandomAccessFile(fname,"r");
			FileChannel fc=raf.getChannel();
			
			int count=1;
			while(fc.position()<flength){
				readHeader(fc);
				
				ls.add(fc.position());
				
				fc.position(fc.position()+NX*NY*4+8);	// skip slab for next header
				
				if(count!=1) dump.append(varToString(count));
				else{
					dump.append(headerToString(false,count));
					dump.append(varToString(count));
				}
				
				info.append(headerToString(true,count));
				
				count++;
			}
			
			fc.close();	raf.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		positions=new long[ls.size()];
		for(int i=0,I=ls.size();i<I;i++) positions[i]=ls.get(i).longValue();
	}
	
 	
 	public void readToSlab(int count,float[][] slab){
 		if(slab.length!=NY||slab[0].length!=NX)
 		throw new IllegalArgumentException("invalid dimensions of slab, should be "+NY+" X "+NX);
 		
		ByteBuffer buf=readData(count);
		
		buf.getInt();	// skip 4
		
		for(int j=0;j<NY;j++)
		for(int i=0;i<NX;i++) slab[j][i]=buf.getFloat();
		
		ArrayUtil.reverse(slab);
	}
	
	public float[][] readToSlab(int count){
		ByteBuffer buf=readData(count);
		
		buf.getInt();	// skip 4
		
		float[][] slab=new float[NY][NX];
		
		for(int j=0;j<NY;j++)
		for(int i=0;i<NX;i++) slab[j][i]=buf.getFloat();
		
		ArrayUtil.reverse(slab);
		
		return slab;
	}
	
	public Variable readToVariable(int count){
		Variable v=new Variable("var",new Range(1,1,NY,NX));
		v.setUndef(-9999);
		v.setCommentAndUnit("var in intermediate file");
		
		readToSlab(count,v.getData()[0][0]);
		
		return v;
	}
	
	
	public void dump(){ System.out.println(dump.toString());}
	
	public String toString(){ return info.toString();}
	
	
	/*** helper methods ***/
	private StringBuilder headerToString(boolean appendVar,int num){
		StringBuilder sb=new StringBuilder();
		
		sb.append("Intermediate file: "+fname+" (format: ");
		
		if(IFV==5) sb.append("WPS)\n");
		else if(IFV==4) sb.append("SI)\n");
		else if(IFV==3) sb.append("MM5)\n");
		else throw new IllegalArgumentException("unsupported IFV: "+IFV);
		
		if(appendVar){
			String lvlinfo=null;
			
			if(XLVL==200100) lvlinfo="surface data, ";
			else if(XLVL==201300) lvlinfo="sea-level data, ";
			else lvlinfo=XLVL/100f+" hPa, ";
			
			sb.append("Field     : "+String.format("%3d ",num)+
				new String(FIELD).trim()+" ("+lvlinfo+
				new String(DESC).trim()+", unit:"+
				new String(UNITS).trim()+")\n"
			);
		}
		if(fmt==Format.WPS||fmt==Format.SI)
		sb.append("Source    : "+new String(MAP_SOURCE).trim()+"\n");
		sb.append("Time      : "+new String(HDATE).trim()+" (Fcst: "+XFCST+")\n");
		sb.append("Dimension : "+NX+" X "+NY+" ["+STARTLON+", "+STARTLAT+"] ");
		sb.append((fmt==Format.WPS||fmt==Format.SI)?new String(STARTLOC).trim():"");
		
		sb.append("\nProjection: "+IPROJ+" --");
		switch(IPROJ){
		case 0:	// cylindrical equidistant projection
			sb.append(" Cylindrical equidistant projection\n");
			sb.append("DLon/Dlat : "+DELTALON+"/"+DELTALAT+"\n");
			break;
			
		case 1:	// Mercator projection
			sb.append(" Mercator projection\n");
			sb.append("DX/DY     : "+DX+"/"+DY+"  TrueLat1: "+TRUELAT1+"\n");
			break;
			
		case 3:	// Lambert conformal projection
			sb.append(" Lambert conformal projection\n");
			sb.append("DX/DY     : "+DX+"/"+DY+"  TrueLat1/TrueLat2: "+TRUELAT1+"/"+TRUELAT2+"\n");
			break;
			
		case 4:	// Gaussian projection
			sb.append(" Gaussian projection\n");
			sb.append("DLon/Nlats: "+DELTALON+"/"+NLATS+"  Rearth: "+EARTH_RADIUS+"\n");
			break;
			
		case 5:	// Polar-stereographic projection
			sb.append(" Polar-stereographic projection\n");
			sb.append("DX/DY     : "+DX+"/"+DY+"  Xlonc: "+XLONC+"  TrueLat1: "+TRUELAT1+"\n");
			break;

		default: throw new IllegalArgumentException("Unsupported IPROJ: "+IPROJ);
		}
		
		if(fmt==Format.WPS)
		sb.append(IS_WIND_EARTH_REL?"wind is earth relative\n":"wind is model relative\n");
		
		sb.append("\n");
		
		return sb;
	}
	
	private StringBuilder varToString(int num){
		StringBuilder sb=new StringBuilder();
		
		String lvlinfo=null;
		
		if(XLVL==200100) lvlinfo="surface data, ";
		else if(XLVL==201300) lvlinfo="sea-level data, ";
		else lvlinfo=XLVL/100f+" hPa, ";
		
		sb.append("Field: "+String.format("%3d ",num)+
			new String(FIELD).trim()+" ("+lvlinfo+
			new String(DESC).trim()+", unit:"+
			new String(UNITS).trim()+")\n"
		);
		
		return sb;
	}
	
	private ByteBuffer readData(int count){
		ByteBuffer buf=ByteBuffer.allocate(NX*NY*4+8);
		buf.order(ByteOrder.BIG_ENDIAN);
		
		try{
			RandomAccessFile raf=new RandomAccessFile(fname,"r");
			FileChannel fc=raf.getChannel();
			
			fc.position(positions[count-1]);
			fc.read(buf);
			
			fc.close();
			raf.close();
			
		}catch(IOException e){ e.printStackTrace();System.exit(0);}
		
		buf.clear();
		
		return buf;
	}
	
	private void readHeader(FileChannel fc) throws IOException{
		////////// record 1 /////////////
		Record record1=new Record(4);
		fc.read(record1.buf);
		record1.buf.clear();
		getRecord1(record1.buf);
		
		////////// record 2 /////////////
		switch(IFV){
		case 5: fmt=Format.WPS; break;
		case 4: fmt=Format.SI ; break;
		case 3: fmt=Format.MM5; break;
		default: throw new IllegalArgumentException("unsupported IFV: "+IFV);
		}
		
		Record record2=null;
		
		if(fmt==Format.WPS||fmt==Format.SI){
			int len=4*5+HDATE.length+FIELD.length+UNITS.length+DESC.length+MAP_SOURCE.length;
			record2=new Record(len);
			
		}else{
			int len=4*5+HDATE.length+FIELD.length+UNITS.length+DESC.length;
			record2=new Record(len);
		}
		
		fc.read(record2.buf);
		record2.buf.clear();
		getRecord2(record2.buf);
		
		////////// record 3 /////////////
		Record record3=null;
		int len=0;
		switch(IPROJ){
		case 0:	// cylindrical equidistant projection
			len=4*4;
			if(fmt==Format.WPS) len+=STARTLOC.length+4;
			if(fmt==Format.SI ) len+=STARTLOC.length;
			break;
			
		case 1:	// Mercator projection
			len=4*5;
			if(fmt==Format.WPS) len+=STARTLOC.length+4;
			if(fmt==Format.SI ) len+=STARTLOC.length;
			break;
			
		case 3:	// Lambert conformal projection
			len=4*7;
			if(fmt==Format.WPS) len+=STARTLOC.length+4;
			if(fmt==Format.SI ) len+=STARTLOC.length;
			break;
			
		case 4:	// Gaussian projection
			len=4*4;
			if(fmt==Format.WPS) len+=STARTLOC.length+4;
			if(fmt==Format.SI ) len+=STARTLOC.length;
			break;
			
		case 5:	// Polar-stereographic projection
			len=4*6;
			if(fmt==Format.WPS) len+=STARTLOC.length+4;
			if(fmt==Format.SI ) len+=STARTLOC.length;
			break;
			
		default: throw new IllegalArgumentException("Unsupported IPROJ: "+IPROJ);
		}
		
		record3=new Record(len);
		fc.read(record3.buf);
		record3.buf.clear();
		getRecord3(record3.buf);
		
		////////// record 4 /////////////
		if(fmt==Format.WPS){
			Record record4=new Record(4);
			fc.read(record4.buf);
			record4.buf.clear();
			getRecord4(record4.buf);
		}
	}
	
	private void getRecord1(ByteBuffer buf){
		buf.getInt();	// skip 4
		IFV=buf.getInt();
	}
	
	private void getRecord2(ByteBuffer buf){
		buf.getInt();	// skip 4
		
		buf.get(HDATE);
		XFCST=buf.getFloat();
		if(fmt==Format.WPS||fmt==Format.SI) buf.get(MAP_SOURCE);
		buf.get(FIELD);
		buf.get(UNITS);
		buf.get(DESC );
		XLVL=buf.getFloat();
		NX=buf.getInt();
		NY=buf.getInt();
		IPROJ=buf.getInt();
	}
	
	private void getRecord3(ByteBuffer buf){
		buf.getInt();	// skip 4
		
		switch(IPROJ){
		case 0:	// cylindrical equidistant projection
			if(fmt==Format.WPS||fmt==Format.SI) buf.get(STARTLOC);
			STARTLAT=buf.getFloat();
			STARTLON=buf.getFloat();
			DELTALAT=buf.getFloat();
			DELTALON=buf.getFloat();
			if(fmt==Format.WPS) EARTH_RADIUS=buf.getFloat();
			break;
			
		case 1:	// Mercator projection
			if(fmt==Format.WPS||fmt==Format.SI) buf.get(STARTLOC);
			STARTLAT=buf.getFloat();
			STARTLON=buf.getFloat();
			DX=buf.getFloat();
			DY=buf.getFloat();
			TRUELAT1=buf.getFloat();
			if(fmt==Format.WPS) EARTH_RADIUS=buf.getFloat();
			break;
			
		case 3:	// Lambert conformal projection
			if(fmt==Format.WPS||fmt==Format.SI) buf.get(STARTLOC);
			STARTLAT=buf.getFloat();
			STARTLON=buf.getFloat();
			DX=buf.getFloat();
			DY=buf.getFloat();
			XLONC=buf.getFloat();
			TRUELAT1=buf.getFloat();
			TRUELAT2=buf.getFloat();
			if(fmt==Format.WPS) EARTH_RADIUS=buf.getFloat();
			break;
			
		case 4:	// Gaussian projection
			if(fmt==Format.WPS||fmt==Format.SI) buf.get(STARTLOC);
			STARTLAT=buf.getFloat();
			STARTLON=buf.getFloat();
			NLATS=buf.getFloat();
			DELTALON=buf.getFloat();
			if(fmt==Format.WPS) EARTH_RADIUS=buf.getFloat();
			break;
			
		case 5:	// Polar-stereographic projection
			if(fmt==Format.WPS||fmt==Format.SI) buf.get(STARTLOC);
			STARTLAT=buf.getFloat();
			STARTLON=buf.getFloat();
			DX=buf.getFloat();
			DY=buf.getFloat();
			XLONC=buf.getFloat();
			TRUELAT1=buf.getFloat();
			if(fmt==Format.WPS) EARTH_RADIUS=buf.getFloat();
			break;
			
		default: throw new IllegalArgumentException("Unsupported IPROJ: "+IPROJ);
		}
	}
	
	private void getRecord4(ByteBuffer buf){
		buf.getInt();	// skip 4
		IS_WIND_EARTH_REL=buf.getInt()==0?false:true;
	}
	
	private static final class Record{
		//
		int dlen=0;	// data  length
		
		ByteBuffer buf=null;
		
		public Record(int len){
			this.dlen=len;
			
			buf=ByteBuffer.allocate(dlen+8);
			buf.order(ByteOrder.BIG_ENDIAN);
		}
	}
	
	
	/** test
	public static void main(String[] args){
		//IntermediateFile imf=new IntermediateFile("d:/FILE_2010-10-19_00.GFS.SI");
		IntermediateFile imf=new IntermediateFile("d:/FILE_2004-09-11_00.interim.WPS");
		imf.dump();
		
		Variable v=imf.readToVariable(97);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream("d:/imf.dat");
		cdws.writeData(v);	cdws.closeFile();
	}*/
}
