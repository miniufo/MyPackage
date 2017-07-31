/**
 * @(#)AccessGISST.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.database;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import static miniufo.basic.ArrayUtil.splitByLength;


/**
 * parse GISST into binary
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AccessGISST{
	//
	private static int ygrid=180;
	private static int xgrid=360;
	
	private static float undef=-32768f;
	private static float cover=-1000f;
	
	
	/**
	 * prevent from construction
	 */
	private AccessGISST(){}
	
	
	//
	public static void parseGISST2Bin(String in,String out){
		try{
			Variable v=new Variable(null,new Range(1,1,ygrid,xgrid));
			
			float[][] vdata=v.getData()[0][0];
			
			Scanner sc=new Scanner(new File(in));
			
			CtlDataWriteStream cdws=new CtlDataWriteStream(out);
			
			while(sc.hasNextLine()){
				System.out.println(sc.nextLine());
				
				for(int j=0;j<ygrid;j++){
					String[] nums=splitByLength(sc.nextLine(),6);
					
					for(int i=0;i<180;i++){
						float f=Float.parseFloat(nums[i]);
						
						if(f==cover) vdata[ygrid-j-1][i+180]=undef;
						else if(f!=undef) vdata[ygrid-j-1][i+180]=f/100;
						else vdata[ygrid-j-1][i+180]=undef;
					}
					
					for(int i=180;i<xgrid;i++){
						float f=Float.parseFloat(nums[i]);
						
						if(f==cover) vdata[ygrid-j-1][i-180]=undef;
						else if(f!=undef) vdata[ygrid-j-1][i-180]=f/100;
						else vdata[ygrid-j-1][i-180]=undef;
					}
				}
				
				cdws.writeData(v);
			}
			
			sc.close();	cdws.closeFile();
		}
		catch(IOException e1){ e1.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1871-1880","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1881-1890","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1891-1900","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1901-1910","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1911-1920","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1921-1930","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1931-1940","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1941-1950","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1951-1960","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1961-1970","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1971-1980","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1981-1990","d:/Data/GISST.dat");
			parseGISST2Bin("d:/Data/GISST/GISST23B_SST_1991-2000","d:/Data/GISST.dat");
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
