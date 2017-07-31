/**
 * @(#)WriteGrADSMap.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

import java.io.FileOutputStream;
import java.io.IOException;
import miniufo.io.FileWriteInterface;
import static miniufo.io.FileWriteInterface.Solution.*;


/**
 * used to write the GrADS map
 *
 * @version 1.0, 02/01/2007
 * @author  MiniLynn
 * @since   MDK1.0
 */
public final class GrADSMapCreator{
	//
	private FileOutputStream   fos=null;	// file access object
	private FileWriteInterface fwi=null;
	
	
	/**
     * constructor
     *
     * @param	file_path		Map file path
     */
	public GrADSMapCreator(String file_path){
		fwi=new FileWriteInterface(file_path);
		
		if(fwi.getFlag()!=SKIP){
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
					throw new IllegalArgumentException("unsupported flag for "+fwi.getFlag());
				}
		    	
		    }catch(IOException ex1){ ex1.printStackTrace(); System.exit(0);}
		}
	}
	
	
	/**
     * write map records
     *
     * @param	GrADSMapRecord	GrADSMapRecord(s)
     */
	public void writeMap(GrADSMapRecord... gmr){
		for(int m=0;m<gmr.length;m++){
			int count=gmr[m].num;
			
			if(count>255)
				throw new IllegalArgumentException("The number of grids in one record exceeds 255 ("+count+")");
			
			byte[] buffer=new byte[3+6*count];
			
			//write one map record in bytes
			buffer[0]=(byte)gmr[m].recsty;
			buffer[1]=(byte)gmr[m].linsty;
			buffer[2]=(byte)gmr[m].num;
			
			for(int i=0;i<count;i++){
				buffer[i*6+0+3]=(byte)((int)(gmr[m].lon[i]*10000f)>>16);
				buffer[i*6+1+3]=(byte)((int)(gmr[m].lon[i]*10000f)>>8);
				buffer[i*6+2+3]=(byte)((int)(gmr[m].lon[i]*10000f));
				
				buffer[i*6+3+3]=(byte)((int)((gmr[m].lat[i]+90f)*10000f)>>16);
				buffer[i*6+4+3]=(byte)((int)((gmr[m].lat[i]+90f)*10000f)>>8);
				buffer[i*6+5+3]=(byte)((int)((gmr[m].lat[i]+90f)*10000f));
			}
			
			try{ fos.write(buffer);}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
	}
	
	
	/**
	 * close file method
     */
	public void closeFile(){
		try{ if(fos!=null) fos.close();}
	    catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String[] arg){		
		try{
			String[] s=("abc miniufo").split("miniufo");
			
			System.out.println(s[0]);
			System.out.println(s[1]);
			
	    }catch(Exception e){ e.printStackTrace();}
	}*/
}
