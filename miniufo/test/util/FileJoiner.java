/**
 * @(#)FileJoiner.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;

import java.io.File;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import miniufo.io.FileWriteInterface;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import static miniufo.io.FileWriteInterface.Solution.*;


/**
 * used to join kinds of data file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class FileJoiner{
	//
	protected boolean is_skip       =true;		// whether to skip writting file
	
	protected final int buffer_size =1048576*2;	// size of buffer (Byte) ==2MB
	
	protected byte[] byt_buffer     =null;		// buffer
	
	protected FileInputStream    fis=null;
	protected FileOutputStream   fos=null;
	protected FileWriteInterface fwi=null;
	
	
	/**
     * constructor
     *
     * @param	str_newFilePath		a string specified the file to write
     */
	public FileJoiner(String str_newFilePath){
		try{
			fwi=new FileWriteInterface(str_newFilePath);
			
			if(fwi.getFlag()!=SKIP){
				is_skip=false;
				
				switch(fwi.getFlag()){
				case RENAME:
					fos=new FileOutputStream(fwi.getParent()+fwi.getNewName());
					break;
					
				case OVERWRITE:
					fos=new FileOutputStream(str_newFilePath);
					break;
					
				case APPEND:
					fos=new FileOutputStream(str_newFilePath,true);
					break;

				default:
					throw new IllegalArgumentException("unsupported flag for "+fwi.getFlag());
				}
			}
	    	
	    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
	}
	
	
	/**
     * join data file just append one by one
     *
     * @param	arr_filePaths	strings specified the files to be joined
     */
	public void joinFile(String[] paths){
		if(paths.length<=1)
			throw new IllegalArgumentException("The files to be joined should be more than one");
		
		if(!is_skip){
			System.out.println("Start joining files...");
			
			try{
				File file =null;
		    	byt_buffer=new byte[buffer_size];
				
				for(int i=0;i<paths.length;i++){
					file=new File(paths[i]);
					fis =new FileInputStream(file);
					
					System.out.println("  With: "+file.getAbsolutePath());
					
					for(int j=0;j<(int)(file.length()/buffer_size);j++){
						fis.read(byt_buffer);
						fos.write(byt_buffer);
					}
					
					int tmp=(int)(file.length()%buffer_size);
					
					if(tmp!=0){
						byte[] ex_byt_buffer=null;	ex_byt_buffer=new byte[tmp];
						fis.read(ex_byt_buffer);
						fos.write(ex_byt_buffer);
						ex_byt_buffer=null;
					}
					
					fis.close();	fis=null;	file=null;
				}
				
				byt_buffer=null;
		    	
		    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
			
			System.out.println("Finish joining files.");
		}
	}
	
	/**
     * join different times file, the spacial dimension is the same
     *
     * @param	paths	strings specified the files to be joined
     */
	public void joinTFile(String[] paths){
		if(paths.length<=1)
			throw new IllegalArgumentException("The files to be joined should be more than one");
		
		if(!is_skip){
			System.out.println("Start joining T-files...");
			
			try{
				File file =null;
				byt_buffer=new byte[buffer_size];
				
				file=new File(paths[0]);
				long length=file.length();	file=null;
				
				for(int i=0;i<paths.length;i++){
					file=new File(paths[i]);
					fis =new FileInputStream(file);
					
					System.out.println("  With: "+file.getAbsolutePath());
					
					if(file.length()!=length) throw new IllegalArgumentException("file lengths not same");
					
					int tmp1=(int)(file.length()/buffer_size);
					int tmp2=(int)(file.length()%buffer_size);
					
					for(int j=0;j<tmp1;j++){
						fis.read(byt_buffer);
						fos.write(byt_buffer);
					}
					
					if(tmp2!=0){
						byte[] ex_byt_buffer=null;	ex_byt_buffer=new byte[tmp2];
						fis.read(ex_byt_buffer);
						fos.write(ex_byt_buffer);
						ex_byt_buffer=null;
					}
					
					fis.close();	fis=null;	file=null;
				}
			
				byt_buffer=null;
		    	
		    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
			
			System.out.println("Finish joining files.");
		}
	}
	
	/**
     * join different variable file, the spacial and temporal dimension is the same
     *
     * @param	paths	strings specified the files to be joined
     * @param	cpr				ctl parser used to specified the dimension
     */
	public void joinZFile(String[] paths,CtlDescriptor cd){
		if(paths.length<=1)
			throw new IllegalArgumentException("The files to be joined should be more than one");
		
		int tcount=cd.getTCount(),vcount=cd.getVCount();
		long onelevellength=cd.getOneLevelLength(),treclength=cd.getTRecLength(),filelength;
		
		if(!is_skip){
			System.out.println("Start joining files...");
			
			try{
				File file =null;	file=new File(paths[0]);	filelength=file.length();
				byt_buffer=new byte[(int)onelevellength];
				
				for(int i=1;i<paths.length;i++){
					file=new File(paths[i]);
					if(file.length()!=filelength) throw new IllegalArgumentException("file lengths not same");
					file=null;
				}
				
				for(int l=0;l<tcount;l++){
					System.out.println("  With: "+(l+1)+" time");
					
					for(int m=0;m<vcount;m++){
						for(int i=0;i<paths.length;i++){
							fis=new FileInputStream(new File(paths[i]));
							
							fis.skip(treclength*l);	fis.skip(onelevellength*m);
							
							fis.read(byt_buffer);
							fos.write(byt_buffer);
							
							fis.close();	fis=null;
						}
					}
				}
			
				byt_buffer=null;
		    	
		    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
			
			System.out.println("Finish joining files.");
		}
	}
	
	/**
     * join different variable file, the spacial and temporal dimension is the same
     *
     * @param	paths	strings specified the files to be joined
     * @param	var_length		length of the variable (bytes)
     * @param	t_count			time count
     */
	public void joinVFile(String[] paths,long var_length,int t_count){
		if(paths.length<=1)
			throw new IllegalArgumentException("The files to be joined should be more than one");
		
		if(!is_skip){
			System.out.println("Start joining files...");
			
			try{
				File file=null;
				byt_buffer=new byte[(int)var_length];	//buffer_size
				
				file=new File(paths[0]);
				long length=file.length();	file=null;
				
				for(int i=1;i<paths.length;i++){
					file=new File(paths[i]);
					if(file.length()!=length) throw new IllegalArgumentException("file lengths not same");
					file=null;
				}
				
				for(int l=0;l<t_count;l++){
					System.out.println("  With: "+(l+1)+" time");
					
					for(int j=0;j<paths.length;j++){
						file=new File(paths[j]);
						fis =new FileInputStream(file);
						fis.skip(var_length*l);
						fis.read(byt_buffer);
						fos.write(byt_buffer);
						fis.close();	fis=null;	file=null;
					}
				}
			
				byt_buffer=null;
		    	
		    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
			
			System.out.println("Finish joining files.");
		}
	}
	
	/**
     * join different variable file, the spacial and temporal dimension is the same
     *
     * @param	paths	strings specified the files to be joined
     * @param	cpr				ctl parser used to specified the dimension
     */
	public void joinVFile(String[] paths,CtlDescriptor cd){
		joinVFile(paths,cd.getTRecLength(),cd.getTCount());
	}
	
	/**
     * join different variable file, the area dimension and time dimension are the same
     * 
     * @param	cds		CtlDescriptors
     */
	public void joinFile(CtlDescriptor[] cds){
		System.out.println("Start joining files...");
		
		try{
			fos.close();	if(fwi.getFile().exists()) fwi.getFile().delete();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		int fcount=cds.length,tcount=cds[0].getTCount();
		
		for(int i=1;i<fcount;i++) if(!cds[i-1].isAreaLike(cds[i]))
		throw new IllegalArgumentException("areas of descriptors are not same");
		
		Range[] r=new Range[fcount];
		
		Variable[] v=new Variable[fcount];
		
		DataRead[] cdrs=new DataRead[fcount];
		
		for(int f=0;f<fcount;f++){
			   r[f]=new Range("t(1,1);z(1,1)",cds[f]);
			   v[f]=new Variable("",r[f]);
			cdrs[f]=DataIOFactory.getDataRead(cds[f]);
		}
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(fwi.getFile().getAbsolutePath());
		
		for(int l=0;l<tcount;l++){
			System.out.println("  With: "+(l+1)+" time");
			
			for(int f=0;f<fcount;f++){
				int vcount=cds[f].getVCount();	Range tmp=v[f].getRange();
				
				tmp.getTRange()[0]=tmp.getTRange()[1]=l+1;
				
				for(int m=0;m<vcount;m++){
					v[f].setName(cds[f].getVDef()[m].getName());
					
					int zcount=cds[f].getVDef()[m].getZCount();
					
					for(int k=0;k<zcount;k++){
						tmp.getZRange()[0]=tmp.getZRange()[1]=k+1;
						
						cdrs[f].readData(v[f]);	cdws.writeData(v[f]);
					}
				}
			}
		}
		
		cdws.closeFile();
		
		for(int f=0;f<fcount;f++) cdrs[f].closeFile();
		
		System.out.println("Finish joining files.");
	}
	
	
	/**
	 * close file method
     */
	public void closeFile(){
		try{
			is_skip=true;	byt_buffer=null;	fwi=null;
			if(fos!=null) fos.close();	fos=null;
			if(fis!=null) fis.close();	fis=null;
	    	
	    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String[] args){
		FileOutputStream fos1=new FileOutputStream("d:/test.txt");
		FileOutputStream fos2=new FileOutputStream("d:/test.txt");
		
		fos2.write(2);
		fos1.write(1);
		
		fos1.close();fos1.close();
		fos2.close();fos2.close();
	}*/
}
