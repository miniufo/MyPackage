/**
 * @(#)EOFResult.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * used to descripe the result of EOF analysis
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EOFResult{
	//
	private int modeCount=0;
	private int eigsCount=0;
	
	private float coeff=0;
	private float ratio=0;	// ratio of mcsum/svsum i.e., percentage of mcsum
	
	private Variable[] modes=null;
	private Variable[] times=null;
	private Variable[] contr=null;
	
	
	/**
     * constructor
     *
     * @param	modes	modes of EOF
     * @param	times	principle component of EOF
     * @param	eigen	eigenvalue of EOF
     */
	public EOFResult(int mc,int ec,Variable v){
		if(mc>ec)
		throw new IllegalArgumentException("eigen count should be larger than mode count");
		
		modeCount=mc;
		eigsCount=ec;
		
		coeff=(float)(Math.sqrt(2.0/v.getTCount()));
		
		modes=new Variable[modeCount];
		times=new Variable[modeCount];
		contr=new Variable[2];
		
		newVariables(v);
		setRanges(v);
	}
	
	
	/*** getor and setor ***/
	public int  getModeCount(){ return modeCount;}
	
	public int getEigenCount(){ return eigsCount;}
	
	public float getErrorCoefficient(){ return coeff;}
	
	public float getFirstFewRatio(){ return ratio;}
	
	public Variable[]  getModes(){ return modes;}
	
	public Variable[]  getTimes(){ return times;}
	
	public Variable[] getContributions(){ return contr;}
	
	
	void setRatio(float ratio){ this.ratio=ratio;}
	
	
	/*** helper methods ***/
	private void newVariables(Variable v){
		int t=v.getTCount();	int z=v.getZCount();
		int y=v.getYCount();	int x=v.getXCount();
		
		float undef=v.getUndef();	String name=v.getName();
		
		for(int m=0;m<modeCount;m++){
			modes[m]=new Variable(name+(m+1),true ,new Range(1,z,y,x));
			times[m]=new Variable(name+(m+1),false,new Range(t,z,1,1));
			
			modes[m].setCommentAndUnit((m+1)+" modes of "+name);
			times[m].setCommentAndUnit((m+1)+" time series of "+name);
			
			modes[m].setUndef(undef);	modes[m].setValue(undef);
			times[m].setUndef(undef);
		}
		
		contr[0]=new Variable("contrib",false,new Range(eigsCount,z,1,1));
		contr[1]=new Variable("conterr",false,new Range(eigsCount,z,1,1));
		
		contr[0].setUndef(undef);	contr[0].setCommentAndUnit("contributions (percentage)");
		contr[1].setUndef(undef);	contr[1].setCommentAndUnit("errors of contributions (percentage)");
	}
	
	private void setRanges(Variable v){
		// process the range
		Range range=v.getRange();
		
		for(int m=0,M=modes.length;m<M;m++){
			Range range0=modes[m].getRange();			Range range1=times[m].getRange();
			
			range0.setZRange(range.getZRange()[0]);		range1.setZRange(range.getZRange()[0]);
			range0.setTRange(range.getTRange()[0]);		range1.setTRange(range);
			range0.setXRange(range);					range1.setXRange(range.getXRange()[0]);
			range0.setYRange(range);					range1.setYRange(range.getYRange()[0]);
		}
		
		for(int i=0;i<2;i++){
			Range eirange=contr[i].getRange();
			
			eirange.setXRange(range.getXRange()[0]);
			eirange.setYRange(range.getYRange()[0]);
			eirange.setZRange(range.getZRange()[0]);
			
			eirange.getTRange()[0]=range.getTRange()[0];
			eirange.getTRange()[2]=eigsCount;
			eirange.getTRange()[1]=eigsCount-1+range.getTRange()[0];
		}
	}
	
	
	/**
	 * used to print out
	 */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		float[] cv=contr[0].getData()[0][0][0];
		float[] ce=contr[1].getData()[0][0][0];
		
		sb.append("\n Result of EOF:\n");
		for(int m=0;m<modeCount;m++){
			sb.append(" Contribution of ");
			sb.append(m+1);
			sb.append(" mode: ");
			sb.append(cv[m]);
			sb.append(" +/- ");
			sb.append(ce[m]);
			sb.append("\n");
		}
		
		return sb.toString();
	}
}
