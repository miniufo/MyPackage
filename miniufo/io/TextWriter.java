/**
 * @(#)TextWriter.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;


/**
 * read data from text file (ASCII format)
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TextWriter{
	
	/**
	 * write a slice of data in x-y plane
	 * 
	 * @param	fname		complete file name and path
	 * @param	data		data
	 * @param	format		format string
	 */
	public static void writeSlice(String fname,float[][] data,String format){
		int y=data.length;
		int x=data[0].length;
		
		try(BufferedWriter br=new BufferedWriter(new FileWriter(fname))){
			for(int j=0;j<y;j++){
				StringBuilder sb=new StringBuilder();
				
				for(int i=0;i<x;i++) sb.append(String.format(format,data[j][i]));
				
				sb.append("\n");
				
				br.write(sb.toString());
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/**
	 * write a slice of data in x-y plane
	 * 
	 * @param	fname		complete file name and path
	 * @param	data		data
	 * @param	format		format string
	 */
	public static void writeMultiColumnData(String fname,float[][] data){
		int y=data.length;
		int x=data[0].length;
		
		try(BufferedWriter br=new BufferedWriter(new FileWriter(fname))){
			for(int i=0;i<x;i++){
				StringBuilder sb=new StringBuilder();
				
				for(int j=0,J=y-1;j<J;j++) sb.append(data[j][i]+"  ");
				
				sb.append(data[y-1][i]+"\n");
				
				br.write(sb.toString());
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/*** test **
	public static void main(String[] args){
		miniufo.diagnosis.Variable u=miniufo.diagnosis.DiagnosisFactory.getVariables
		("d:/Data/NCEP/OriginalNC/uwnd.sig995.mon.mean.nc","lon(80,200);lat(0,60);t(1,768)","uwnd")[0];
		
		miniufo.diagnosis.Variable v=miniufo.diagnosis.DiagnosisFactory.getVariables
		("d:/Data/NCEP/OriginalNC/vwnd.sig995.mon.mean.nc","lon(80,200);lat(0,60);t(1,768)","vwnd")[0];
		
		miniufo.diagnosis.Variable um=u.anomalizeT();
		miniufo.diagnosis.Variable vm=v.anomalizeT();
		
		//FilterMethods.cycleFilter(u,12);
		//FilterMethods.cycleFilter(v,12);
		
		miniufo.diagnosis.Variable[] re=TurbulenceMethodsInSC.cVarianceEllipse(u,v);
		
		float[][][][] ra=re[0].getData();
		float[][][][] rb=re[1].getData();
		float[][][][] an=re[2].getData();
		float[][][][] uu=um.getData();
		float[][][][] vv=vm.getData();
		
		float[][] data=new float[7][u.getXCount()*u.getYCount()];
		
		for(int j=0,J=u.getYCount(),ptr=0;j<J;j++)
		for(int i=0,I=u.getXCount();i<I;i++){
			data[0][ptr]=ra[0][j][i][0];
			data[1][ptr]=rb[0][j][i][0];
			data[2][ptr]=an[0][j][i][0];
			data[3][ptr]=i*2.5f+80;
			data[4][ptr]=j*2.5f;
			data[5][ptr]=uu[0][j][i][0];
			data[6][ptr]=vv[0][j][i][0];
			
			ptr++;
		}
		
		writeMultiColumnData("d:/Data/NCEP/SurfWind/dataForMatlab.txt",data);
	}*/
}
