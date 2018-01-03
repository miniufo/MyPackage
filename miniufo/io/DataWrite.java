/**
 * @(#)DataWrite.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Variable;


/**
 * interface for I/O
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class DataWrite implements Print{
	
	/**
	 * to write data from the specified file
	 *
     * @param	v	variable need to fill with data
     */ 
	public abstract void writeData(Variable... v);
	
	public abstract void writeData(DataDescriptor dd,Variable... v);
	
	public abstract void writeCtl(DataDescriptor dd);
	
	public abstract void writeCtl(DataDescriptor dd,float[] zdef,String tinc);
	
	/**
	 * close file method
     */
	public abstract void closeFile();
	
	
	/*** helper class ***/
	protected static final class Var{
		//
		int tcount=0;		// t-count
		int zcount=0;		// z-count
		int ycount=0;		// y-count
		int xcount=0;		// x-count
		
		int tstart=0;
		int zstart=0;
		int ystart=0;
		int xstart=0;
		
		float undef=Float.NaN;	// undefined value
		
		String name=null;	// variable name
		String cmmt=null;	// comment
		
		//
		public Var(Variable v){
			this.tcount=v.getTCount();
			this.zcount=v.getZCount();
			this.ycount=v.getYCount();
			this.xcount=v.getXCount();
			
			this.tstart=v.getRange().getTRange()[0];
			this.zstart=v.getRange().getZRange()[0];
			this.ystart=v.getRange().getYRange()[0];
			this.xstart=v.getRange().getXRange()[0];
			
			this.undef =v.getUndef();
			this.name  =v.getName();
			this.cmmt  =v.getCommentAndUnit();
		}
		
		//
		public boolean isAreaLike(Var v){
			if(xcount!=v.xcount||ycount!=v.ycount) return false;
			return true;
		}
	}
	
}
