/**
 * @(#)Var.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import java.util.regex.Pattern;

/**
 * used to descripe the variable in the DataDescriptor file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class Var{
	//
	protected int tcount=1;
	protected int zcount=1;
	protected int ycount=1;
	protected int xcount=1;
	
	protected float undef=Float.NaN;
	
	protected String unit   =null;
	protected String vname  =null;
	protected String comment=null;
	
	protected static final Pattern unitPtn=Pattern.compile("\\([^\\(\\)]+?\\)");
	
	
	/*** getor and setor ***/
	public int getTCount(){ return tcount;}
	
	public int getZCount(){ return zcount;}
	
	public int getYCount(){ return ycount;}
	
	public int getXCount(){ return xcount;}
	
	public float getUndef(){ return undef;}
	
	public String getUnit(){ return unit ;}
	
	public String getName(){ return vname;}
	
	public String getComment(){ return comment;}
	
	public void setTCount(int t){ this.tcount=t;}
	
	public void setZCount(int z){ this.zcount=z;}
	
	public void setYCount(int y){ this.ycount=y;}
	
	public void setXCount(int x){ this.xcount=x;}
	
	public void setUndef(float a){ this.undef=a;}
	
	public void setUnit(String unit) { this.unit =unit ;}
	
	public void setName(String vname){ this.vname=vname;}
	
	public void setComment(String comment){ this.comment=comment;}
	
	
	/**
    * used to print out
    */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		sb.append("Varname:\t");	sb.append(vname);	sb.append("\n");
		sb.append("t-count:\t");	sb.append(tcount);	sb.append("\n");
		sb.append("z-count:\t");	sb.append(zcount);	sb.append("\n");
		sb.append("y-count:\t");	sb.append(ycount);	sb.append("\n");
		sb.append("x-count:\t");	sb.append(xcount);	sb.append("\n");
		sb.append("unit   :\t");	sb.append(unit);	sb.append("\n");
		sb.append("undef  :\t");	sb.append(undef);	sb.append("\n");
		sb.append("comment:\t");	sb.append(comment);	sb.append("\n");
		
		return sb.toString();
	}
}
