/**
 * @(#)TextReader.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;


/**
 * read data from text file (ASCII format)
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TextReader{
	
	/**
	 * prevent from instantiate
	 */
	private TextReader(){}
	
	
	/**
	 * read data of the in the cols column and return them in float type
	 * 
	 * @param	fname		complete file name and path
	 * @param	skipTitle	whether skip the first line
	 * @param	cols		the colth column that need to return (start from 1)
	 */
	public static float[][] readColumnsF(String fname,boolean skipTitle,int... cols){
		int lines=IOUtil.getLineNumber(fname);
		
		if(skipTitle) lines--;
		
		float[][] data=new float[cols.length][lines];
		
		try(BufferedReader br=new BufferedReader(new InputStreamReader(new FileInputStream(fname),getCharSet(fname)))){
			String oneline=null;
			
			if(skipTitle) oneline=br.readLine();
			
			for(int j=0;j<lines;j++){
				oneline=br.readLine().trim();
				
				String[] ss=oneline.split("[\\s\\t]+");
				
				for(int i=0,I=cols.length;i<I;i++)
				data[i][j]=Float.parseFloat(ss[cols[i]-1]);
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		return data;
	}
	
	/**
	 * read data of the in the cols column and return them in String type
	 * 
	 * @param	fname		complete file name and path
	 * @param	skipTitle	whether skip the first line
	 * @param	cols		the colth column that need to return (start from 1)
	 */
	public static String[][] readColumnsS(String fname,boolean skipTitle,int... cols){
		int lines=IOUtil.getLineNumber(fname);
		
		if(skipTitle) lines--;
		
		String[][] data=new String[cols.length][lines];
		
		try(BufferedReader br=new BufferedReader(new InputStreamReader(new FileInputStream(fname),getCharSet(fname)))){
			String oneline=null;
			
			if(skipTitle) oneline=br.readLine();
			
			for(int j=0;j<lines;j++){
				oneline=br.readLine().trim();
				
				String[] ss=oneline.split("[\\s\\t]+");
				
				for(int i=0,I=cols.length;i<I;i++) data[i][j]=ss[cols[i]-1];
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		return data;
	}
	
	
	/*** helper methods ***/
	private static String getCharSet(String fname){
		int p=0;
		
		try(BufferedInputStream bis=new BufferedInputStream(new FileInputStream(fname))){
			p=(bis.read()<<8)+bis.read();
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		switch(p){
		case 0xefbb: return "UTF-8";
		case 0xfffe: return "Unicode";
		case 0xfeff: return "UTF-16BE";
		default    : return "GBK";
		}
	}
	
	
	/*** test **
	public static void main(String[] args){
		float[] data=readColumn("D:/Data/GDP/SCS/EddySpeedPDF/dataU.txt",1826,false,1)[0];
		
		for(float f:data)
		System.out.println(f);
	}*/
}
