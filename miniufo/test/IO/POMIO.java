/**
 * @(#)POMIO.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.IO;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import miniufo.basic.ArrayUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * utility for I/O of Princeton Ocean Model
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class POMIO{
	//
	private static final String FORMAT="%13.5f";
	
	
	public static void writeUV
	(float[][] lons,float[][] lats,float[][] uvel,float[][] vvel,String fname){
		checkDim(lons,lats,uvel,vvel);
		
		BufferedWriter bw=null;
		
		try{
			bw=new BufferedWriter(new FileWriter(new File(fname)));
			
			bw.write(" LON\n");		write(lons,bw);
			bw.write(" LAT\n");		write(lats,bw);
			bw.write(" U-VEL\n");	write(uvel,bw);
			bw.write(" V-VEL\n");	write(vvel,bw);
			
			bw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public static void writeTS
	(float[][] lons,float[][] lats,float[] zdef,float[][][] temp,float[][][] salt,String fname){
		int z=zdef.length;
		
		if(z!=salt.length)
		throw new IllegalArgumentException("z-lengths of temp and salt are not the same");
		
		BufferedWriter bw=null;
		
		try{
			bw=new BufferedWriter(new FileWriter(new File(fname)));
			
			bw.write(" LEVELS\n");	write(zdef,bw);
			bw.write(" LON\n");		write(lons,bw);
			bw.write(" LAT\n");		write(lats,bw);
			
			for(int k=0;k<z;k++){
				checkDim(lons,lats,temp[k],salt[k]);
				
				bw.write(" T   "+(k+1)+"\n");	write(temp[k],bw);
				bw.write(" S   "+(k+1)+"\n");	write(salt[k],bw);
			}
			
			bw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public static void writeBath
	(float lon,float lat,float del,float[][] bath,String fname){
		BufferedWriter bw=null;
		
		try{
			bw=new BufferedWriter(new FileWriter(new File(fname)));
			
			bw.write(lon+"   "+lat+"   "+(1f/del)+"\n");
			write(bath,bw);
			
			bw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/**
	 * helper methods
	 * write out a 2D field of data
	 *
     * @param	data	a 2D field of data
     */ 
	private static void write(float[][] data,BufferedWriter bw){
		try{
			int y=data.length;
			int x=data[0].length;
			int bufLen=Integer.parseInt(FORMAT.substring(1,3))*x;
			
			for(int j=y-1;j>=0;j--){ // from north to south
				StringBuilder sb=new StringBuilder(bufLen);
				
				for(int i=0;i<x;i++)
				sb.append(String.format(FORMAT,data[j][i]));
				
				bw.write(sb.append("\n").toString());
			}
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	private static void write(float[] data,BufferedWriter bw){
		try{
			int x=data.length;
			int bufLen=Integer.parseInt(FORMAT.substring(1,3))*x+(x+1)/7;
			
			StringBuilder sb=new StringBuilder(bufLen);
			
			for(int i=0;i<x;i++){
				sb.append(String.format(FORMAT,data[i]));
				
				if((i+1)%15==0) sb.append("\n");
			}
			
			bw.write(sb.append("\n").toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	private static void checkDim(float[][] lons,float[][] lats,float[][] data1,float[][] data2){
		int y=lons.length;
		int x=lons[0].length;
		
		if(y!=lats.length||x!=lats[0].length)
		throw new IllegalArgumentException("lon/lat dimensions not same");
		
		if(y!=data1.length||x!=data1[0].length)
		throw new IllegalArgumentException("lon/data1 dimensions not same");
		
		if(y!=data2.length||x!=data2[0].length)
		throw new IllegalArgumentException("lon/data2 dimensions not same");
	}
	
	
	/** test*/
	public static void main(String[] args){
		String range="lon(90,145);lat(-10,30)";
		
		DiagnosisFactory df=DiagnosisFactory.parseFile("e:/SODA/Climatology/Climatology.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] ts =df.getVariables(new Range(range+";t(1,1)"       ,dd),"temp","salt");
		Variable[] tau=df.getVariables(new Range(range+";t(1,1);z(1,1)",dd),"taux","tauy");
		Variable bath =DiagnosisFactory.getVariables("e:/ETOPO/ETOPO5.ctl",range,"bath")[0];
		
		System.out.println(bath);
		
		ts[0].setDataMin(-100);	tau[0].setDataMin(-100);
		ts[1].setDataMin(-100);	tau[1].setDataMin(-100);
		
		int x=ts[0].getXCount(),y=ts[0].getYCount();
		
		float[] zdef=dd.getZDef().getSamples().clone();
		float[][] badata=bath.getData()[0][0];
		float[][] txdata=tau[0].getData()[0][0];
		float[][] tydata=tau[1].getData()[0][0];
		float[][][] tdata=ts[0].getData()[0];
		float[][][] sdata=ts[1].getData()[0];
		
		float[][] lons=new float[y][x];
		float[][] lats=new float[y][x];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			lons[j][i]= 90+i*0.5f;
			lats[j][i]=-10+j*0.5f;
		}
		
		ArrayUtil.reverse(zdef);
		ArrayUtil.reverse(tdata);
		ArrayUtil.reverse(sdata);
		
		writeTS(lons,lats,zdef,tdata,sdata,"d:/SODA_TS.dat");
		writeUV(lons,lats,txdata,tydata,"d:/SODA_UV.dat");
		writeBath(90,-10,0.083333333f,badata,"d:/SODA_BATH.dat");
	}
}
