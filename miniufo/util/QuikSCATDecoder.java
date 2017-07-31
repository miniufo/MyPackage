/**
 * @(#)QuikSCATDecoder.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.IOUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * used to decode the QuikSCAT wind data
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class QuikSCATDecoder{
	//
	private static final int y=353;	// y-grid count
	private static final int x=720;	// x-grid count
	
	
	/**
     * decode the data
     *
     * @param	inpath		input path
     * @param	outpath		output path
     * @param	writeCtl	whether to write the ctl file
     */
	public static void decode(String inpath,String outpath,boolean writeCtl){
		try{
			RandomAccessFile raf=new RandomAccessFile(inpath,"r");
			FileChannel fc=raf.getChannel();
			
			CtlDataWriteStream cdws=new CtlDataWriteStream(outpath);
			
			Variable u=new Variable("u",new Range(1,1,y,x));
			Variable v=new Variable("v",new Range(1,1,y,x));
			
			int reclen=(2*x*y+7)*4;
			int t=(int)raf.length()/reclen;
			
			ByteBuffer buf=ByteBuffer.allocate(reclen);
			buf.order(ByteOrder.LITTLE_ENDIAN);
			
			float[][] udata=u.getData()[0][0];
			float[][] vdata=v.getData()[0][0];
			
			long total=System.currentTimeMillis();
			
			for(int l=0;l<t;l++){
				fc.read(buf);	buf.position(16);	// skip 4*4 bytes
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) udata[j][i]=buf.getFloat();
				
				buf.position((6+x*y)<<2);	// skip 4*4 + x*y*4 + 2*4 bytes
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) vdata[j][i]=buf.getFloat();
				
				cdws.writeData(u);	cdws.writeData(v);
				
				buf.clear();
			}
			
			System.out.println("\ntotal: "+(System.currentTimeMillis()-total)/1000f);
			
			fc.close(); raf.close(); cdws.closeFile();
			
			// write ctl
			if(writeCtl) writeCTL(inpath,outpath);
		}
		catch(IOException e1){ e1.printStackTrace();System.exit(0);}
	}
	
	public static void decode(String inpath,String outpath){
		decode(inpath,outpath,true);
	}
	
	
	/**
     * write Ctl file
     */
	private static void writeCTL(String ipath,String opath){
		final String[] months={"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dev"};
		
		String ctlpath=IOUtil.getCompleteFileNameWithoutExtension(opath)+".ctl";
		
		FileWriter fw=null;
		
		try{ fw=new FileWriter(ctlpath);}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		final int[] days={
			31,28,31,30,31,30,31,31,30,31,30,31
		};
		
		int month=1,tcnt=1,year=2000;
		String[] dates=IOUtil.getFileName(ipath).split("\\.");
		
		if(dates.length>1){
			year =Integer.parseInt(dates[1].substring(0,4));
			month=Integer.parseInt(dates[1].substring(4,6));
		}
		
		tcnt=days[month-1]*4;
		
		try{
			fw.write(
				"dset ^"+IOUtil.getFileName(opath)+"\n"+
				"undef -9999.0\n"+
				"title QuikSCAT wind data\n"+
				"xdef 720 linear    0.5 0.5\n"+
				"ydef 353 linear  -88.0 0.5\n"+
				"zdef   1 linear 1000.0 1.0\n"+
				"tdef "+tcnt+" linear 06z01"+months[month-1]+year+" 6hr\n"+
				"vars 2\n"+
				"u 0 99 u-wind\n"+
				"v 0 99 v-wind\n"+
				"endvars\n"
			);
			
			fw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String[] args){
		for(int yy=2000;yy<=2009;yy++)
		QuikSCATDecoder.decode("//Lynn/QuikSCAT/uv."+yy+"05.bln","d:/lynn/uv"+yy+".dat",true);
	}*/
}
