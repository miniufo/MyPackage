/**
 * @(#)IOUtil.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.stream.Stream;


/**
 * utility for I/O
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class IOUtil{
	
	/**
	 * get line number of a text file
	 *
     * @param	completePath	complete path includes path and file name
	 */
	public static int getLineNumber(String completePath){
		try(LineNumberReader lnr=new LineNumberReader(new FileReader(completePath))){
			lnr.skip(Long.MAX_VALUE);
			return lnr.getLineNumber();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		throw new IllegalArgumentException("should not reach here");
	}
	
	
	/**
	 * Get file names.
	 *
     * @param	completePath	complete path includes path and file name
     */
	public static String getFileName(String completePath){ return new File(completePath).getName();}
	
	public static String getFilePrefixName(String completePath){
		String name=getFileName(completePath);
		
		int idx=name.lastIndexOf(".");
		
		if(idx==-1)
			return name;
		else
			return name.substring(0,name.lastIndexOf("."));
	}
	
	public static String getFileSuffixName(String completePath){
		String name=getFileName(completePath);
		
		int idx=name.lastIndexOf(".");
		
		if(idx==-1)
			return "";
		else
			return name.substring(name.lastIndexOf(".")+1);
	}
	
	public static String getCompleteFileNameWithoutExtension(String completePath){
		int idx=completePath.lastIndexOf(".");
		
		if(idx==-1)
			return completePath;
		else
			return completePath.substring(0,idx);
	}
	
	public static String getCompletePath(String completePath){
		int idx=completePath.lastIndexOf("/");
		
		if(idx==-1){
			idx=completePath.lastIndexOf("\\");
			
			if(idx==-1) return completePath;
			else return completePath.substring(0,idx)+"\\";
			
		}else return completePath.substring(0,idx)+"/";
	}
	
	
	/**
	 * Get file length in integer.
	 *
     * @param	fname	a given file name
     */
	public static int getFileLength(String fname){
		File f=new File(fname);
		
		if(!f.exists()) throw new IllegalArgumentException("file not found: "+fname);
		
		int fileLen=(int)f.length();
		
		if(fileLen!=f.length()) throw new IllegalArgumentException("file length overflow for integer");
		
		return fileLen;
	}
	
	
	public static void replaceContent(String fname,String oldContent,String newContent){
		String content=null;
		
		try(Stream<String> s=Files.lines(Paths.get(fname))){
			Optional<String> op=s.reduce((a,b)->a+"\n"+b);
			
			if(op.isPresent()) content=op.get()+"\n";
			else throw new IllegalArgumentException("no content");
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		content=content.replaceAll(oldContent,newContent);
		
		try(FileWriter fw=new FileWriter(fname)){ fw.write(content);}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String[] args){
		String path="d:/Data/dd.TXT/abc.cc.dd";
		
		System.out.println(getFileName(path));
		System.out.println(getFilePrefixName(path));
		System.out.println(getFileSuffixName(path));
		System.out.println(getCompleteFileNameWithoutExtension(path));
		System.out.println(getCompletePath(path));
	}*/
}
