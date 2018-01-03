/**
 * @(#)CsmDataWriteStream.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.FileWriter;
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import miniufo.io.FileWriteInterface;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Variable;
import static miniufo.io.FileWriteInterface.Solution.*;


/**
 * used to write the commonest csm data file (station type)
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CsmDataWriteStream extends DataWrite{
	//
	private boolean is_skip       =true;	// whether to skip overwriting the file
	private boolean print         =true;
	
	private String file_path      =null;
	private String file_name      =null;
	
	private FileChannel      fc   =null;	// file channel
	private FileOutputStream fos  =null;	// file access object
	private FileWriteInterface fwi=null;
	
	protected ArrayList<Var>  vars=null;
	
	
	/**
     * constructor
     *
     * @param	ctl_path		a string specified the ctl data file path
     * @param	byteswapped		whether need to wap byte
     */
	public CsmDataWriteStream(String file_path){
		fwi=new FileWriteInterface(file_path);
		this.file_path=fwi.getFile().getAbsolutePath();
		this.file_name=fwi.getFile().getName();
		
		if(fwi.getFlag()!=SKIP){
			is_skip=false;
			
			try{
				switch(fwi.getFlag()){
				case RENAME:
					fos=new FileOutputStream(fwi.getParent()+fwi.getNewName());
					break;
					
				case OVERWRITE:
					fos=new FileOutputStream(file_path);
					break;
					
				case APPEND:
					fos=new FileOutputStream(file_path,true);
					break;

				default:
					throw new IllegalArgumentException("unknown case");
				}
				
		    }catch(FileNotFoundException ex){ ex.printStackTrace(); System.exit(0);}
		}
		
		fc=fos.getChannel();
		vars=new ArrayList<Var>();
	}
	
	
	/**
     * write data in station data form
     *
     * @param	t_first		whether t dimension is the first dimension
     * @param	data		the storage of the data
     * @param	lon			longitudes of the stations
     * @param	lat			latitudes of the stations
     */
	public void writeData(DataDescriptor dd,Variable... v){
		if(!(dd instanceof CsmDescriptor))
			throw new IllegalArgumentException("CsmDescriptor only");
		
		CsmDescriptor cd=(CsmDescriptor)dd;
		
		if(!is_skip){
			if(print) System.out.println("\nStart writing station data...");
			
			write(cd.getLon(),cd.getLat(),cd.getZDef().getSamples(),v);
			
			/*** write ctl file ***/
			writeCtl(cd);
		    
		    if(print) System.out.println("Finish writing station data.");
		}
	}
	
	public void writeData(float[][][] lon,float[][][] lat,float[] zdef,Variable... v){
		if(!is_skip){
			if(print) System.out.println("\nStart writing station data...");
			write(lon,lat,zdef,v);
	    	if(print) System.out.println("Finish writing station data.");
		}
	}
	
	public void writeData(Variable... v){
		throw new IllegalArgumentException("use this method as: writeData(CsmDescriptor, Variable...)");
	}
	
	
	/**
     * Write data in station data form and all records are in the same time.
     * This function can be called consecutively to write records of different times.
     *
     * @param	srs		station records
     */
	public void writeSnapshot(StationRecord[] srs){
		if(srs.length<1) throw new IllegalArgumentException("no station record");
		
		for(StationRecord sr:srs) writeStationRecord(sr);
		writeStationRecordTerminator(srs[srs.length-1]);
	}
	
	
	/**
     * generate the ctl file of the given Variable according to the given CtlDescriptor
     *
     * @param	v	given Variable
     * @param	cd	given CtlDescriptor
     */
	public void writeCtl(DataDescriptor dd){ writeCtl(dd,null,null);}
	
	public void writeCtl(DataDescriptor dd,float[] zdef,String tinc){
		if(!(dd instanceof CsmDescriptor)) throw new IllegalArgumentException("CsmDescriptor only");
		
		CsmDescriptor cd=(CsmDescriptor)dd;
		
		String[] s=file_path.split("\\.");
		StringBuilder name=new StringBuilder();
		for(int i=0;i<s.length-1;i++) name.append(s[i]+".");
		String npath=name.append("ctl").toString();
		
		is_skip=true;	int vcount=vars.size();
		fwi=new FileWriteInterface(npath);
		Var last=vars.get(vcount-1);
		FileWriter fw=null;
		
		try{
			if(fwi.getFlag()!=SKIP){
				is_skip=false;
				
				switch(fwi.getFlag()){
				case RENAME   : fw=new FileWriter(fwi.getParent()+fwi.getNewName()); break;
				case OVERWRITE: fw=new FileWriter(npath); break;
				case APPEND   :
					System.out.println("Ctl file does not support appending");
					System.out.println("The program will terminate now");
					System.exit(0);
					break;
				default       : throw new IllegalArgumentException("unknown case");
				}
			}
	    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		if(!is_skip){
			StringBuilder sb=new StringBuilder();	int start=0;
			
			String mapname=file_name.substring(0,file_name.lastIndexOf("."))+".map";
			
			sb.append("dset ^");		sb.append(file_name);
			sb.append("\ntitle ");		sb.append(cd.getTitle());
			sb.append("\ndtype ");		sb.append("station");
			sb.append("\nstnmap ^");	sb.append(mapname);
			sb.append("\nundef ");		sb.append(cd.getUndef(null));
			
			start=last.tstart;
			sb.append("\ntdef ");
			sb.append(last.tcount);
			sb.append(" linear ");
			sb.append(cd.getTDef().getSamples()[start-1].toGradsDate());
			sb.append(" ");
			sb.append(tinc==null?cd.getTIncrement():tinc);
			
			sb.append("\nvars ");
			sb.append(vcount);
			sb.append("\n");
			
			for(Var v:vars){
				sb.append(String.format("%-9s",v.name));
				sb.append(" ");
				sb.append(v.zcount==1?0:v.zcount);
				sb.append(" 99 "+v.cmmt+"\n");
			}
			
			sb.append("endvars\n");
			
			try{ fw.write(sb.toString()); fw.close();}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		vars.clear();
	}
	
	
	/**
	 * whether to print out
	 *
     * @param	print	print or disable print
     */ 
	public void setPrinting(boolean print){ this.print=print;}
	
	
	/**
	 * close file method
     */
	public void closeFile(){
		is_skip=true;
		
		try{ if(fos!=null){ fc.close(); fos.close();}}
	    catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		fwi=null;	fos=null;
	}
	
	
	/*** helper methods ***/
	
	/**
     * write data in station data form
     *
     * @param	lon		longitudes of the stations
     * @param	lat		latitudes of the stations
     * @param	zdef	levels of the stations
     * @param	v		variables need to be written
     */
	private void write(float[][][] lon,float[][][] lat,float[] zdef,Variable... v){
		int vcount=v.length;
		
		int t=v[0].getTCount(),	z=v[vcount-1].getZCount(),	y=v[0].getYCount(),	x=v[0].getXCount();
		float undef=v[0].getUndef();
		
		if(lon.length!=lat.length)
			throw new IllegalArgumentException("t-dimensions not same");
		if(lon[0].length!=lat[0].length)
			throw new IllegalArgumentException("y-dimensions not same");
		if(lon[0][0].length!=lat[0][0].length)
			throw new IllegalArgumentException("x-dimensions not same");
		
		int std_id=500000,svc=0,lvc=0,tstart=v[0].getRange().getTRange()[0];
		
		if(v[0].getZCount()==1) svc++;
		else if(z!=1&&v[0].getZCount()==z) lvc++;
		
		for(int m=1;m<vcount;m++){
			if(v[m].getZCount()==1) svc++;
			else if(z!=1&&v[m].getZCount()==z) lvc++;
			
			if(v[m].getZCount()<v[m-1].getZCount())
			throw new IllegalArgumentException("invalid variable sequence");
			
			if(v[m].getTCount()!=t||v[m].getYCount()!=y||v[m].getXCount()!=x)
			throw new IllegalArgumentException("areas or times are not same");
			
			if(v[m].isTFirst()!=v[m-1].isTFirst())
			throw new IllegalArgumentException("tfirsts are not same");
			
			if(!Float.isNaN(undef)&&v[m].getUndef()!=undef)
			throw new IllegalArgumentException("undefs are not same");
			
			if(x!=lon[0][0].length||y!=lon[0].length||
				lon[0].length!=lat[0].length||lon[0][0].length!=lat[0][0].length
			)throw new IllegalArgumentException("dimensions not same");
		}
		
		if(svc+lvc!=vcount) throw new IllegalArgumentException("z counts are not valid");
		
		StationRecord sr=new StationRecord();
		sr.tim=0;	sr.flag=0;	sr.nlev=z;
		sr.sdata=new float[svc];
		if(lvc>0) sr.ldata=new float[z][lvc+1];
		if(svc>0){ sr.flag=1; if(lvc>0) sr.nlev++;}
		
		float[][][][][] data=new float[vcount][][][][];
		
		for(int m=0;m<vcount;m++){ data[m]=v[m].getData(); vars.add(new Var(v[m]));}
		
		if(v[0].isTFirst()){
			for(int l=0;l<t;l++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					sr.sid="  "+std_id;	std_id++;
					sr.lat=lat[tstart-1+l][j][i];
					sr.lon=lon[tstart-1+l][j][i];
					
					// surface data
					for(int m=0;m<svc;m++) sr.sdata[m]=data[m][l][0][j][i];
					
					// levels  data
					int zstart=v[0].getRange().getZRange()[0];
					if(lvc>0)
					for(int k=0;k<z;k++){
						sr.ldata[k][0]=zdef[zstart-1+k];
						for(int m=svc;m<vcount;m++) sr.ldata[k][m+1-svc]=data[m][l][k][j][i];
					}
					
					writeStationRecord(sr);
				}
				
				writeStationRecordTerminator(sr);	std_id=500000;
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					sr.sid="  "+std_id;	std_id++;
					sr.lat=lat[tstart-1+l][j][i];
					sr.lon=lon[tstart-1+l][j][i];
					
					// surface data
					for(int m=0;m<svc;m++) sr.sdata[m]=data[m][0][j][i][l];
					
					// levels  data
					int zstart=v[0].getRange().getZRange()[0];
					if(lvc>0)
					for(int k=0;k<z;k++){
						sr.ldata[k][0]=zdef[zstart-1+k];
						for(int m=svc;m<vcount;m++) sr.ldata[k][m+1-svc]=data[m][k][j][i][l];
					}
					
					writeStationRecord(sr);
				}
				
				writeStationRecordTerminator(sr);	std_id=500000;
			}
		}
	}
	
	/**
     * called by writeStationData to write a single station record
     *
     * @param	sr	a station record
     */
	private void writeStationRecord(StationRecord sr){
		int size=0;
		
		if(sr.ldata==null)
			size=sr.sid.length()+(5+sr.sdata.length)*4;
		else
			size=sr.sid.length()+(5+sr.sdata.length+sr.ldata.length*sr.ldata[0].length)*4;
		
		ByteBuffer buf=ByteBuffer.allocate(size);
		
		buf.order(ByteOrder.nativeOrder());
		buf.put(sr.sid.getBytes());	//sid
		buf.putFloat(sr.lat);		// lat
		buf.putFloat(sr.lon);		// lon
		buf.putFloat(sr.tim);		// tim
		buf.putInt(sr.nlev);		// nlev
		buf.putInt(sr.flag);		// nflag 1
		
		// surface data
		if(sr.flag==1){
		for(int i=0,I=sr.sdata.length;i<I;i++)
		buf.putFloat(sr.sdata[i]);}
		
		// levels data
		if(sr.nlev-sr.flag>=1)
		for(int i=0,I=sr.ldata.length;i<I;i++)
		for(int j=0,J=sr.ldata[0].length;j<J;j++)
		buf.putFloat(sr.ldata[i][j]);
		
		buf.flip();
		
		try{ fc.write(buf);}
		catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
	}
	
	/**
     * called by writeStationData to write a single station record terminator (do not write data)
     *
     * @param	sr	a station record
     */
	private void writeStationRecordTerminator(StationRecord sr){
		ByteBuffer buf=ByteBuffer.allocate(sr.sid.length()+5*4);
		
		buf.order(ByteOrder.nativeOrder());
		buf.put(sr.sid.getBytes());	//sid
		buf.putFloat(sr.lat);		//lat
		buf.putFloat(sr.lon);		//lon
		buf.putFloat(sr.tim);		//tim
		buf.putInt(0);				//nlev
		buf.putInt(1);				//nflag
		buf.flip();
		
		try{ fc.write(buf);}
		catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
	}
}
