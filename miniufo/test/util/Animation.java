/**
 * @(#)Animation.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.io.IOUtil;


/**
 * change a variable into an animation of GIF format
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Animation{
	//
	private int delay=10;		// default 10 ms
	
	private Variable vsrc=null;	// URL vector
	private Variable vdes=null;	// namelist
	
	
	/**
     * constructor
     */
	public Animation(Variable v,int N,Type type){
		vsrc=v;
		
		int c=0;
		
		if(v.getXCount()==1) c++;
		if(v.getYCount()==1) c++;
		if(v.getZCount()==1) c++;
		
		if(c==0&&c==3)
		throw new IllegalArgumentException("invalid dimensions of given variable");
		
		vdes=vsrc.interpolateT(N,type);
	}
	
	
	/**
     * generate GIF animation
     * 
     * @param	dd			data descriptor of the variable
     * @param	gifpath		path of GIF file
     */
	public void generateGIF(DataDescriptor dd,String gifpath){
		String name=IOUtil.getCompleteFileNameWithoutExtension(gifpath);
		
		String  gspath=name+".gs" ;
		String ctlpath=name+".ctl";
		String datpath=name+".dat";
		
		toDatFile(dd,datpath);
		generateGS(ctlpath,datpath,gspath,gifpath);
		
		OpenGrADS.runGS(gspath);
	}
	
	
	private void toDatFile(DataDescriptor dd,String path){
		DataWrite dw=DataIOFactory.getDataWrite(dd,path);
		dw.writeData(dd,vdes);	dw.closeFile();
	}
	
	private void generateGS(String ctlpath,String datpath,String gspath,String gifpath){
		String folder=new File(ctlpath).getParent()+File.separator+"_tmp_GIF_tmp_"+File.separator;
		
		StringBuilder sb=new StringBuilder();
		
		float min=ArrayUtil.getMin(vdes.getData());
		float max=ArrayUtil.getMax(vdes.getData());
		int cint=(int)((max-min)/235f);
		
		int rmin=(int)(min-min%cint);
		int rmax=(int)(max-max%cint);
		
		sb.append("'open "+ctlpath+"'\n");
		sb.append("'!mkdir "+folder+"'\n");
		sb.append("'set gxout shaded'\n");
		sb.append("'set mproj scaled'\n\n");
		sb.append("tt=1\n");
		sb.append("while(tt<="+vdes.getTCount()+")\n");
		sb.append("nn=100000+tt\n");
		sb.append("'set t 'tt\n");
		sb.append("'set grads off'\n");
		sb.append("'setNewRGBs'\n");			// run user-defined gs
		sb.append("'set rbrange "+rmin+" "+rmax+"'\n");
		sb.append("'set cint "+cint+"'\n");
		sb.append("'set ccols ");
		for(int ccols=254;ccols>=16;ccols--) sb.append(ccols+" ");
		sb.append("'\n");
		sb.append("'set clevs ");
		for(int levs=rmin,count=0;count<238;levs+=cint,count++)
		sb.append(levs+" ");
		sb.append("'\n");
		sb.append("'d "+vdes.getName()+"'\n");
		sb.append("'printim "+folder+"'nn'.png x640 y480 white'\n");
		sb.append("'c'\n");
		sb.append("tt=tt+1\n");
		sb.append("endwhile\n\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n\n");
		sb.append("'!convert -trim -delay "+delay+" "+folder+"*.png "+gifpath+"'\n");
		sb.append("'!rm -rf "+folder+"'\n");
		sb.append("'!rm  "+ctlpath+"'\n");
		sb.append("'!rm  "+datpath+"'\n");
		sb.append("'!rm  "+gspath+"'\n");
		sb.append("'quit'\n");
		
		try{
			FileWriter fw=new FileWriter(gspath);
			fw.write(sb.toString());	fw.close();
			
		}catch(IOException e){ e.printStackTrace();System.exit(0);}
	}
	
	
	/** test
	public static void main(String argv[]){
		miniufo.diagnosis.DiagnosisFactory df=new miniufo.diagnosis.DiagnosisFactory(
			"J:/DataI/NCEP1/Daily/Pressure/hgt/hgt.2008.nc"
		);
		DataDescriptor dd=df.getDataDescriptor();
		Variable hgt=df.getVariables(new Range("time(2008.1.1,2008.2.29);lev(600,600)",dd),"hgt")[0];
		
		Animation am=new Animation(hgt,237,Type.CUBIC);
		am.generateGIF(dd,"e:/hgt500.gif");
	}*/
}
