/**
 * @(#)DataComposite.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;

import miniufo.io.DataRead;
import miniufo.io.DataIOFactory;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * composition class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DataComposite{
	//
	private DataDescriptor dd=null;
	
	
	/**
     * constructor
     *
     * @param	src		source DataDescriptor
     */
	public DataComposite(DataDescriptor dd){ this.dd=dd;}
	
	
	/**
     * composite in T dimension
     *
     * @param	vname	name of a variable
     * @param	ranges	ranges for composite
     * 
     * @return	com		composite result, t-length is the same as each range
     */
	public Variable compositeT(String vname,String[] ranges){
		return compositeT(vname,initializeRanges(ranges));
	}
	
	public Variable compositeT(String vname,Range[] ranges){
		int rlength=ranges.length;
		
		int tcount=ranges[0].getTRange()[2];	int zcount=ranges[0].getZRange()[2];
		int ycount=ranges[0].getYRange()[2];	int xcount=ranges[0].getXRange()[2];
		
		for(int i=0;i<rlength;i++) if(
			ranges[i].getTRange()[2]!=tcount||ranges[i].getZRange()[2]!=zcount||
			ranges[i].getYRange()[2]!=ycount||ranges[i].getXRange()[2]!=xcount
		) throw new IllegalArgumentException("ranges are not dimensional the same");
		
		Variable tmp=new Variable(vname,new Range(tcount,zcount,ycount,xcount));
		Variable com=new Variable(vname,new Range(tcount,zcount,ycount,xcount));
		
		com.setCommentAndUnit("composite "+dd.getVarComment(vname));
		
		Range tmprg=tmp.getRange();	Range comrg=com.getRange();
		
		tmprg.setZRange(ranges[0]);
		tmprg.setYRange(ranges[0]);
		tmprg.setXRange(ranges[0]);
		
		comrg.setTRange(ranges[0]);	comrg.setZRange(ranges[0]);
		comrg.setYRange(ranges[0]);	comrg.setXRange(ranges[0]);
		
		DataRead cdrs=DataIOFactory.getDataRead(dd);
		
		for(int r=0;r<rlength;r++){
			tmprg.setTRange(ranges[r]);
			
			cdrs.readData(tmp);
			
			com.plusEq(tmp);
		}
		
		cdrs.closeFile();
		
		com.divideEq(rlength);	com.setUndef(dd.getUndef(""));
		
		return com;
	}
	
	
	/**
     * extract according to ranges before composite, t-size of ranges should be 1,
     * the T-length of result is the same as the length of ranges array
     *
     * @param	vname	name of a variable
     * @param	ranges	ranges for extraction
     * 
     * @return	res		result, t-length is the same as the count of ranges
     */
	public Variable extractT(String vname,String[] ranges){
		return extractT(vname,initializeRanges(ranges));
	}
	
	public Variable extractT(String vname,Range[] ranges){
		int rlength=ranges.length;
		
		int tcount=ranges[0].getTRange()[2];	int zcount=ranges[0].getZRange()[2];
		int ycount=ranges[0].getYRange()[2];	int xcount=ranges[0].getXRange()[2];
		
		if(tcount!=1) throw new IllegalArgumentException("tcount should be 1");
		
		for(int i=0;i<rlength;i++) if(
			ranges[i].getTRange()[2]!=tcount||ranges[i].getZRange()[2]!=zcount||
			ranges[i].getYRange()[2]!=ycount||ranges[i].getXRange()[2]!=xcount
		) throw new IllegalArgumentException("ranges are not dimensional the same");
		
		Variable tmp=new Variable(vname,true ,new Range(tcount ,zcount,ycount,xcount));
		Variable res=new Variable(vname,false,new Range(rlength,zcount,ycount,xcount));
		
		res.setUndef(dd.getUndef(""));
		res.setCommentAndUnit("extracted "+dd.getVarComment(vname));
		
		Range comrg=res.getRange();
		Range tmprg=tmp.getRange();
		
		DataRead cdrs=DataIOFactory.getDataRead(dd);
		
		for(int r=0;r<rlength;r++){
			tmprg.setZRange(ranges[r]);	tmprg.setYRange(ranges[r]);
			tmprg.setXRange(ranges[r]);	tmprg.setTRange(ranges[r]);
			
			cdrs.readData(tmp);
			
			float[][][]   tdata=tmp.getData()[0];
			float[][][][] cdata=res.getData();
			
			for(int k=0;k<zcount;k++)
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++) cdata[k][j][i][r]=tdata[k][j][i];
		}
		
		cdrs.closeFile();
		
		comrg.setZRange(ranges[0]);
		comrg.setYRange(ranges[0]);
		comrg.setXRange(ranges[0]);
		
		int[] tr=comrg.getTRange();
		tr[0]=ranges[0].getTRange()[0];
		tr[2]=ranges[0].getTRange()[2];
		tr[1]=tr[2]+tr[0]-1;
		
		return res;
	}
	
	
	/**
     * initialize ranges
     *
     * @param	src		source DataDescriptor
     * 
     * @return	r		ranges
     */
	private Range[] initializeRanges(String[] ranges){
		int N=ranges.length;
		
		Range[] r=new Range[N];
		
		for(int i=0;i<N;i++) r[i]=new Range(ranges[i],dd);
		
		return r;
	}
	
	
	/** test
	public static void main(String[] args){
		
	}*/
}
