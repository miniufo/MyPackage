/**
 * @(#)CtsDataWriteStream.java	1.0 2017.07.04
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteOrder;
import java.text.DecimalFormat;
import miniufo.io.FileWriteInterface;
import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Variable;
import static miniufo.io.FileWriteInterface.Solution.*;


/**
 * Write the cts data file and corresponding cts file
 *
 * @version 1.0, 2017.07.04
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CtsDataWriteStream extends CtlDataWriteStream{
	
	/**
     * constructor
     *
     * @param	ctl_path		a string specified the ctl data file path
     * @param	byteswapped		whether need to wap byte
     */
	public CtsDataWriteStream(String file_path,ByteOrder order){ super(file_path,order);}
	
	/**
     * constructor
     *
     * @param	ctl_path		a string specified the ctl data file path
     */
	public CtsDataWriteStream(String file_path){ super(file_path,ByteOrder.nativeOrder());}
	
	
	/**
     * generate the ctl file of the given Variable according to the given CtlDescriptor
     *
     * @param	dd		data descriptor
     * @param	zdef	new zdef levels need to be written into new ctl file
     * @param	tinc	time incremence need to be written into new ctl file
     */
	public void writeCtl(DataDescriptor dd,float[] zdef,String tinc){
		if(al.size()==0) return;
		
		String[] s=file_path.split("\\.");
		StringBuilder name=new StringBuilder();
		for(int i=0;i<s.length-1;i++) name.append(s[i]+".");
		String ctlpath=name.append("cts").toString();	is_skip=true;
		
		FileWriter fw=null;
		fwi=new FileWriteInterface(ctlpath);
		
		try{
			if(fwi.getFlag()!=SKIP){
				is_skip=false;
				
				switch(fwi.getFlag()){
				case RENAME   : fw=new FileWriter(fwi.getParent()+fwi.getNewName()); break;
				case OVERWRITE: fw=new FileWriter(ctlpath); break;
				case APPEND   :
					System.out.println("Ctl file does not support appending");
					System.out.println("The program will terminate now");
					System.exit(0); break;
				default       : throw new IllegalArgumentException("unknown case");
				}
			}
	    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		if(!is_skip){
			StringBuilder sb=new StringBuilder();
			Variable last=al.get(al.size()-1);
			
			// pick up the variable has largest z-count
			for(Variable v:al) if(v.getZCount()>last.getZCount()) last=v;
			
			int start=0;
			
			sb.append("dset ^");	sb.append(file_name);
			if(order==ByteOrder.BIG_ENDIAN) sb.append("\noptions big_endian");
			sb.append("\nundef ");	sb.append(al.get(0).getUndef());
			sb.append("\ntitle ");	sb.append(dd.getTitle());
			
			start=last.getRange().getXRange()[0];
			if(dd.xLinear()){
				sb.append("\nxdef "+String.format("%4d",last.getXCount())+" linear ");
				sb.append(dd.getXDef().getSamples()[start-1]+" "+dd.getDXDef()[0]+"\n");
				
			}else{
				if(last.getXCount()!=1){
					sb.append("\nxdef "+String.format("%4d",last.getXCount())+" levels\n");
					for(int i=0;i<last.getXCount();i++)
					sb.append("\t"+dd.getXDef().getSamples()[start-1+i]+"\n");
					
				}else{
					sb.append("\nxdef    1 linear "+dd.getXDef().getSamples()[start-1]);
					sb.append(" "+dd.getDXDef()[0]+"\n");
				}
			}
			
			start=last.getRange().getYRange()[0];
			if(dd.yLinear()){
				sb.append("ydef "+String.format("%4d",last.getYCount())+" linear ");
				sb.append(dd.getYDef().getSamples()[start-1]+" "+dd.getDYDef()[0]+"\n");
				
			}else{
				if(last.getYCount()!=1){
					sb.append("ydef "+String.format("%4d",last.getYCount())+" levels\n");
					for(int i=0;i<last.getYCount();i++)
					sb.append("\t"+dd.getYDef().getSamples()[start-1+i]+"\n");
					
				}else{
					sb.append("ydef    1 linear "+dd.getYDef().getSamples()[start-1]);
					sb.append(" "+dd.getDYDef()[0]+"\n");
				}
			}
			
			if(zdef==null){
				start=last.getRange().getZRange()[0];
				if(dd.zLinear()){
					if(dd.getDZDef()==null||dd.getDZDef()[0]<0){
						sb.append("zdef "+String.format("%4d",last.getZCount())+" levels ");
						
						if(last.getZCount()<=dd.getZDef().length()){
							for(int i=0;i<last.getZCount();i++)
							sb.append(dd.getZDef().getSamples()[start-1+i]+" ");
							
						}else{
							float zstr=dd.getZDef().getSamples()[0];
							for(int i=0;i<last.getZCount();i++) sb.append(zstr-50*i+" ");
						}
						
						sb.append("\n");
						
					}else{
						sb.append("zdef "+String.format("%4d",last.getZCount())+" linear ");
						sb.append(dd.getZDef().getSamples()[start-1]+" "+dd.getDZDef()[0]+"\n");
					}
					
				}else{
					if(last.getZCount()!=1){
						if(dd.getZDef().length()!=1){
							sb.append("zdef "+String.format("%4d",last.getZCount())+" levels ");
							for(int i=0;i<last.getZCount();i++)
							sb.append(dd.getZDef().getSamples()[start-1+i]+" ");
							sb.append("\n");
							
						}else{
							sb.append("zdef "+String.format("%4d",last.getZCount())+" levels ");
							for(int i=0;i<last.getZCount();i++)
							sb.append(dd.getZDef().getSamples()[0]-10*i+" ");
							sb.append("\n");
						}
						
					}else sb.append("zdef    1 levels "+dd.getZDef().getSamples()[start-1]+"\n");
				}
				
			}else{
				sb.append("zdef "+zdef.length+" levels ");
				
				DecimalFormat df=new DecimalFormat(".000");
				for(int i=0;i<zdef.length;i++)
				sb.append(df.format(zdef[i])+" ");
				sb.append("\n");
			}
			
			start=last.getRange().getTRange()[0];
			
			sb.append("tdef "+String.format("%4d",last.getTCount())+" linear ");
			sb.append(dd.getTDef().getSamples()[start-1].toGradsDate());
			sb.append(" "+(tinc==null?dd.getTIncrement():tinc)+"\n");
			
			CtsDescriptor cd=(CtsDescriptor)dd;
			sb.append("f0   "+cd.getF0  ()+"\n");
			sb.append("Beta "+cd.getBeta()+"\n");
			
			sb.append(toStringBuilder(al,nm));
			
			try{ fw.write(sb.toString()); fw.close();}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		al.clear();
	}
	
	
	/** test
	public static void main(String[] arg){		
		try{
			ByteBuffer bb=ByteBuffer.allocate(10);
			System.out.println(bb.isDirect());
			
			bb=ByteBuffer.allocateDirect(10);
			System.out.println(bb.isDirect());
			
	    }catch(Exception e){ e.printStackTrace();}
	}*/
}
