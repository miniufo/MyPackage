/**
 * @(#)CtlDataReadStream.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Range;


/**
 * used to read the binary ctl data file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CtlDataReadStream implements DataRead,Print{
	//
	private int t=0;
	private int z=0;
	private int y=0;
	private int x=0;
	
	private long skipZ=0;
	private long skipY=0;
	private long skipX=0;
	private long skipF=0;
	
	private boolean sequential  =false;
	private boolean print       =true;
	
	private ByteBuffer       buf=null;	// buffer fulfill with data in one time for process
	private CtlDescriptor     cd=null;
	private FileChannel      fc =null;
	private RandomAccessFile raf=null;
	
	
	/**
     * constructor
     *
     * @param	cd	ctl descriptor
     */
	public CtlDataReadStream(CtlDescriptor cd){
		this.cd=cd;	sequential=cd.isSequential();
		
		try{
			raf=new RandomAccessFile(cd.getDSet(),"r");
			
			if(sequential){
				if(raf.length()!=(cd.getTRecLength()+cd.getVCount()*2*4)*cd.getTCount())
				throw new IllegalArgumentException("length of data file is invalid");
				
			}else{
				if(raf.length()!=cd.getTRecLength()*cd.getTCount())
				throw new IllegalArgumentException(
					"length of data file is invalid:"+raf.length()+
					"(data, bytes), "+cd.getTRecLength()*cd.getTCount()+"(ctl, bytes)"
				);
			}
	    }
		catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		fc=raf.getChannel();
	}
	
	
	/**
	 * to read data from the specified file
	 *
     * @param	v	variable need to fill with data
     */ 
	public void readData(Variable... v){
		if(v.length!=1){
			if(print) System.out.print("\nStart reading ");
			
			for(int m=0;m<v.length;m++){
				if(print) System.out.print(v[m].getName()+" ");
				readOne(v[m]);
			}
			
			if(print) System.out.println("data...\nFinish reading data.");
			
		}else readOne(v[0]);
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
		try{ if(raf!=null){ fc.close();	raf.close();}}
		catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		
		t=z=y=x=0;	skipF=skipZ=skipY=skipX=0;
		
		buf=null;	raf=null;	cd=null;
	}
	
	
	/*** helper methods ***/
	private void readOne(Variable v){
		long one_level_length=cd.getOneLevelLength();
		
		boolean yrev=cd.isYRev();
		boolean zrev=cd.isZRev();
		
		Range range=v.getRange();
		
		int[] trange=range.getTRange();	int[] zrange=range.getZRange();
		int[] yrange=range.getYRange();	int[] xrange=range.getXRange();
		
		t=trange[2];	z=zrange[2];
		y=yrange[2];	x=xrange[2];
		
		int tcount=cd.getTCount();	int zcount=cd.getZCount();
		int ycount=cd.getYCount();	int xcount=cd.getXCount();
		
		if(t>tcount||z>zcount) throw new IllegalArgumentException("invalid range");
		if(y>ycount||x>xcount) throw new IllegalArgumentException("invalid range");
		
		/*** calculate skips ***/
		skipF=(xrange[0]-1)<<2;
		
		if(yrev) skipF+=xcount*(ycount-yrange[1])<<2;
		else skipF+=xcount*(yrange[0]-1)<<2;
		
		if(zrev) skipF+=one_level_length*(zcount-zrange[1]);
		else skipF+=one_level_length*(zrange[0]-1);
		
		skipF+=cd.getVarStartPosition(v.getName());
		
		if("99".equals(cd.getStorageType())){
			skipF+=cd.getTRecLength()*(trange[0]-1);
			skipZ=cd.getTRecLength()-one_level_length*zrange[2];
			
		}else{
			skipF+=one_level_length*cd.getVarZcount(v.getName())*(trange[0]-1);
			skipZ=one_level_length*(zcount-zrange[2]);
		}
		
		skipX=(xcount-xrange[2])<<2;
		skipY=(ycount-yrange[2])*xcount<<2;
		
		/*** buffer ***/
		buf=ByteBuffer.allocate(xrange[2]*4);
		buf.order(cd.getByteOrder());
		
		/*** start to read data ***/
		try{
			if(yrev&&!zrev){
				if(v.isTFirst()) readTFYRev(v.getData());
				else readYRev(v.getData());
				
			}else if(zrev&&!yrev){
				if(v.isTFirst()) readTFZRev(v.getData());
				else readZRev(v.getData());
				
			}else if(yrev&&zrev){
				if(v.isTFirst()) readTFYRevZRev(v.getData());
				else readYRevZRev(v.getData());
				
			}else{
				if(v.isTFirst()) readTF(v.getData());
				else read(v.getData());
			}
			
	    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
	    
	    v.setUndef(cd.getUndef(v.getName()));
	    v.setCommentAndUnit(cd.getVarCommentAndUnit(v.getName()));
	}
	
	private void read(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);
		
		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[0][0][i][0]=buf.getFloat();
		
		for(int j=1;j<y;j++){
			raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][j][i][0]=buf.getFloat();
		}
		
		for(int k=1;k<z;k++){	raf.seek(skipY+raf.getFilePointer());
		for(int j=0;j<y;j++){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[k][j][i][0]=buf.getFloat();
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){	raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=0;k<z;k++){	raf.seek(skipY+raf.getFilePointer());
		for(int j=0;j<y;j++){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[k][j][i][l]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
	
	private void readYRev(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);
		
		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[0][y-1][i][0]=buf.getFloat();
		
		for(int j=y-2;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			for(int i=0;i<x;i++)
			data[0][j][i][0]=buf.getFloat();
		}
		
		for(int k=1;k<z;k++){		raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[k][j][i][0]=buf.getFloat();
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){		raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=0;k<z;k++){		raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[k][j][i][l]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
	
	private void readZRev(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);
		
		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[z-1][0][i][0]=buf.getFloat();
		
		for(int j=1;j<y;j++){		raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[z-1][j][i][0]=buf.getFloat();
			
		}
		
		for(int k=z-2;k>=0;k--){	raf.seek(skipY);
		for(int j=0;j<y;j++){		raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[k][j][i][0]=buf.getFloat();
			
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){		raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=z-1;k>=0;k--){	raf.seek(skipY+raf.getFilePointer());
		for(int j=0;j<y;j++){		raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[k][j][i][l]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
	
	private void readYRevZRev(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);

		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[z-1][y-1][i][0]=buf.getFloat();
		
		for(int j=y-2;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[z-1][j][i][0]=buf.getFloat();
		}
		
		for(int k=z-2;k>=0;k--){	raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
		
			for(int i=0;i<x;i++)
			data[k][j][i][0]=buf.getFloat();
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){		raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=z-1;k>=0;k--){	raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[k][j][i][l]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
	
	private void readTF(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);
		
		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[0][0][0][i]=buf.getFloat();
		
		for(int j=1;j<y;j++){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][0][j][i]=buf.getFloat();
		}
		
		for(int k=1;k<z;k++){	raf.seek(skipY+raf.getFilePointer());
		for(int j=0;j<y;j++){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][k][j][i]=buf.getFloat();
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){	raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=0;k<z;k++){	raf.seek(skipY+raf.getFilePointer());
		for(int j=0;j<y;j++){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[l][k][j][i]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
	
	private void readTFYRev(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);
		
		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[0][0][y-1][i]=buf.getFloat();
		
		for(int j=y-2;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][0][j][i]=buf.getFloat();
		}
						
		for(int k=1;k<z;k++){		raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][k][j][i]=buf.getFloat();
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){		raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=0;k<z;k++){		raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[l][k][j][i]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
	
	private void readTFZRev(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);
		
		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[0][z-1][0][i]=buf.getFloat();
		
		for(int j=1;j<y;j++){		raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][z-1][j][i]=buf.getFloat();
		}
		
		for(int k=z-2;k>=0;k--){	raf.seek(skipY+raf.getFilePointer());
		for(int j=0;j<y;j++){		raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][k][j][i]=buf.getFloat();
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){		raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=z-1;k>=0;k--){	raf.seek(skipY+raf.getFilePointer());
		for(int j=0;j<y;j++){		raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[l][k][j][i]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
	
	private void readTFYRevZRev(float[][][][] data) throws IOException{
		raf.seek(skipF);
		
		if(sequential) raf.skipBytes(4);
		
		fc.read(buf); buf.clear();
		for(int i=0;i<x;i++)
		data[0][z-1][y-1][i]=buf.getFloat();
		
		for(int j=y-2;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][z-1][j][i]=buf.getFloat();
		}
		
		for(int k=z-2;k>=0;k--){	raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[0][k][j][i]=buf.getFloat();
		}}
		
		if(sequential) raf.skipBytes(4);
		
		for(int l=1;l<t;l++){		raf.seek(skipZ+raf.getFilePointer());	if(sequential) raf.skipBytes(4);
		for(int k=z-1;k>=0;k--){	raf.seek(skipY+raf.getFilePointer());
		for(int j=y-1;j>=0;j--){	raf.seek(skipX+raf.getFilePointer());
			buf.clear(); fc.read(buf); buf.clear();
			
			for(int i=0;i<x;i++)
			data[l][k][j][i]=buf.getFloat();
			
		}}if(sequential) raf.skipBytes(4);}
	}
}
