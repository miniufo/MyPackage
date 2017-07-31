/**
 * @(#)CtlDataWriteStream.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.FileWriter;
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.text.DecimalFormat;
import java.util.ArrayList;
import miniufo.io.FileWriteInterface;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Variable;
import static miniufo.io.FileWriteInterface.Solution.*;


/**
 * used to write the commonest ctl data file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class CtlDataWriteStream implements DataWrite,Print{
	//
	protected boolean is_skip       =true;	// whether to skip overwriting the file
	protected boolean isUsed99      =false;
	protected boolean print         =true;
	
	protected String storageType    =null;
	protected String file_path      =null;
	protected String file_name      =null;
	
	protected ByteOrder order       =null;
	
	protected FileChannel         fc=null;	// file channel
	protected FileOutputStream   fos=null;	// file access object
	protected FileWriteInterface fwi=null;
	
	protected ArrayList<String>   nm=null;
	protected ArrayList<Variable> al=null;
	
	
	/**
     * constructor
     *
     * @param	ctl_path		a string specified the ctl data file path
     * @param	byteswapped		whether need to wap byte
     */
	public CtlDataWriteStream(String file_path,ByteOrder order){
		fwi=new FileWriteInterface(file_path);
		this.file_path=fwi.getFile().getAbsolutePath();
		this.file_name=fwi.getFile().getName();
		
		if(fwi.getFlag()!=SKIP){
			is_skip=false;
			
			try{
				switch(fwi.getFlag()){
				case RENAME:
					fos=new FileOutputStream(fwi.getParent()+"/"+fwi.getNewName());
					fc=fos.getChannel();
					break;
					
				case OVERWRITE:
					fos=new FileOutputStream(file_path);
					fc=fos.getChannel();
					break;
					
				case APPEND:
					fos=new FileOutputStream(file_path,true);
					fc=fos.getChannel();
					break;

				default:
					throw new IllegalArgumentException("unknown case");
				}
				
		    }catch(FileNotFoundException ex){ ex.printStackTrace(); System.exit(0);}
		}
		
		this.order=order;	fwi=null;
		
		nm=new ArrayList<String>();
		al=new ArrayList<Variable>();
	}
	
	/**
     * constructor
     *
     * @param	ctl_path		a string specified the ctl data file path
     */
	public CtlDataWriteStream(String file_path){ this(file_path,ByteOrder.nativeOrder());}
	
	
	/**
     * write data from array to local ctl data file
     *
     * @param	dd		data descriptor
     * @param	zdef	new zdef levels need to be written into new ctl file
     * @param	tinc	time incremence need to be written into new ctl file
     * @param	v		variables
     */
	public void writeData(Variable... v){
		if(!is_skip){
			if(v.length!=1&&print){
				System.out.print("\nStart writing ");
				
				for(int m=0;m<v.length;m++) System.out.print(v[m].getName()+" ");
				System.out.println("data...");
			}
			
			boolean tfirst=v[0].isTFirst();	float undef=v[0].getUndef();
			
			int tcount=v[0].getTCount();
			int ycount=v[0].getYCount();
			int xcount=v[0].getXCount();
			int vcount=v.length;
			
			for(int m=1;m<vcount;m++){
				if(v[m].getTCount()!=tcount)
					throw new IllegalArgumentException("time count not same for "+v[m].getName());
				if(v[m].getYCount()!=ycount||v[m].getXCount()!=xcount)
					throw new IllegalArgumentException("area count not same for "+v[m].getName());
				if(v[m].isTFirst()!=tfirst)
					throw new IllegalArgumentException("tfirst not same");
				if(!Float.isNaN(undef)&&v[m].getUndef()!=undef)
					throw new IllegalArgumentException("undef not same");
			}
			
			// size: one level grid-count * 4
			ByteBuffer buf=ByteBuffer.allocateDirect(ycount*xcount<<2);
			buf.order(order);
			
			try{
			    if(tfirst){
					for(int l=0;l<tcount;l++)
					for(int m=0;m<vcount;m++){
						int zcount=v[m].getZCount();
						float[][][][] vdata=v[m].getData();
						
						for(int k=0;k<zcount;k++){
							for(int j=0;j<ycount;j++)
							for(int i=0;i<xcount;i++)
							buf.putFloat(vdata[l][k][j][i]);
							
							buf.clear();	fc.write(buf);	buf.clear();
						}
					}
					
				}else{
					for(int l=0;l<tcount;l++)
					for(int m=0;m<vcount;m++){
						int zcount=v[m].getZCount();
						float[][][][] vdata=v[m].getData();
						
						for(int k=0;k<zcount;k++){
							for(int j=0;j<ycount;j++)
							for(int i=0;i<xcount;i++)
							buf.putFloat(vdata[k][j][i][l]);
							
							buf.clear();	fc.write(buf);	buf.clear();
						}
					}
				}
			    
		    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		    
			if(v.length==1){
				if(isUsed99) throw new UnsupportedOperationException(
					"writeData(Variable[]) and writeData(Variable) cannot be called in a same instance"
				);
				
				al.add(v[0]);	nm.add(v[0].getName()+"mini_ufo"+v[0].getCommentAndUnit());
				
				if(al.size()!=1){
					storageType="-1,20";
					
					if(!al.get(al.size()-1).isAreaLike(v[0]))
						throw new IllegalArgumentException("dimensions not same");
					
				}else storageType="99";
				
			}else{
				if(al.size()!=0) throw new UnsupportedOperationException(
					"writeData(Variable) and writeData(Variable[]) cannot be called in a same instance"
				);
				
				isUsed99=true;
				
				al.clear();	storageType="99";
				
				for(int i=0;i<v.length;i++){ al.add(v[i]);	nm.add(v[i].getName()+"mini_ufo"+v[i].getCommentAndUnit());}
				
				if(print) System.out.println("Finish writing data.");
			}
		}
	}
	
	public void writeData(DataDescriptor dd,Variable... v){
		writeData(v);	writeCtl(dd);
	}
	
	public void writeData(DataDescriptor dd,Variable[]... vas){
		for(Variable[] va:vas)
		for(Variable v:va) writeData(v);
		
		writeCtl(dd);
	}
	
	
	/**
     * generate the ctl file of the given Variable according to the given CtlDescriptor
     *
     * @param	dd		data descriptor
     */
	public void writeCtl(DataDescriptor dd){ writeCtl(dd,null,null);}
	
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
		String ctlpath=name.append("ctl").toString();	is_skip=true;
		
		FileWriter fw=null;
		fwi=new FileWriteInterface(ctlpath);
		
		try{
			if(fwi.getFlag()!=SKIP){
				is_skip=false;
				
				switch(fwi.getFlag()){
				case RENAME: fw=new FileWriter(fwi.getParent()+fwi.getNewName()); break;
				case OVERWRITE: fw=new FileWriter(ctlpath); break;
				case APPEND:
					System.out.println("Ctl file does not support appending");
					System.out.println("The program will terminate now");
					System.exit(0); break;
				default: throw new IllegalArgumentException("unknown case");
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
			
			sb.append("tdef ");
			sb.append(String.format("%4d",last.getTCount()));
			sb.append(" linear ");
			sb.append(dd.getTDef().getSamples()[start-1].toGradsDate());
			sb.append(" ");
			sb.append(tinc==null?dd.getTIncrement():tinc);
			
			sb.append("\n");
			sb.append(toStringBuilder(al,nm));
			
			try{ fw.write(sb.toString()); fw.close();}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		al.clear();
	}
	
	
	/**
	 * whether to print out
	 *
     * @param	print	print or disable print
     */ 
	public void setPrinting(boolean print){ this.print=print;}
	
	
	/**
	 * close file method
     */
	public void closeFile(){
		is_skip=true;
		
		al.clear();	al=null;
		
		try{ if(fos!=null){ fc.close(); fos.close();}}
	    catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		fwi=null;	fos=null;	file_path=null;
	}
	
	
	/*** change the info of an array of variables into a string ***/
	protected StringBuilder toStringBuilder(ArrayList<Variable> al,ArrayList<String> nm){
		StringBuilder sb=new StringBuilder();
		
		int length=al.size();
		
		sb.append("vars ");
		sb.append(String.format("%4d",length));
		sb.append("\n");
		for(int m=0;m<length;m++){
			Variable v=al.get(m);	String[] s=nm.get(m).split("mini_ufo");
			
			sb.append(String.format("%-11s",s[0]));
			sb.append(" ");
			sb.append(String.format("%3d",v.getZCount()==1?0:v.getZCount()));
			sb.append(" ");
			sb.append(String.format("%5s",storageType));
			sb.append(" ");
			if(s.length>1) sb.append((s[1].equals("")||s[1].equals("null"))?"comment":s[1]);
			else sb.append("comment");
			sb.append("\n");
		}
		
		sb.append("endvars\n");
		
		return sb;
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
