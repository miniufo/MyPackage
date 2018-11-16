/**
 * @(#)Contours.java	1.0 2017.06.09
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel;


/**
 * This class describe a series of contours used to setup a
 * tracer-contour-based spatial coordinate.
 *
 * @version 1.0, 2017.06.09
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Contours{
	//
	private boolean increSToN=true;	// increasing from south to north
	private boolean isLinear =true;	// equally-spaced contours
	
	private int      number=0;		// number of contours
	
	private double southVal=0;		// contour value at south boundary
	private double northVal=0;		// contour value at north boundary
	private double meanIntv=0;		// mean interval
	
	private double[] values=null;	// values of contours
	
	private ContourTable ct=null;	// tabulated M(q), C(q), A(q) and Yeq(q)
	
	
	/**
	 * Compute the contours given the 2D array data and number of contours.
	 * Whether contours are increasing from south to north is determined by
	 * the data themselves.
	 * 
	 * @param	data		an 2D array data
	 * @param	undef		undefined value
	 * @param	numOfC		number of contours
	 * @param	increSToN	the contours are defined increasing from south to north
	 */
	public Contours(float[][] data,float undef,int numOfC,boolean increSToN){
		float[] ex=ArrayUtil.getExtrema(data,undef);
		
		this.increSToN=increSToN;
		
		if(increSToN){ southVal=ex[0]; northVal=ex[1];}
		else         { southVal=ex[1]; northVal=ex[0];}
		
		number=numOfC;
		meanIntv=(northVal-southVal)/(number-1);
		
		check();
		
		values=new double[number];
		for(int c=0,I=number;c<I;c++) values[c]=southVal+meanIntv*c;
		
		isLinear=true;
		
		ct=new ContourTable();
	}
	
	/**
	 * Compute the contours given the south-most and north-most contours.
	 * 
	 * @param	csouth	south-most contour
	 * @param	cnorth	north-most contour
	 * @param	inc		increment of contour
	 */
	public Contours(float csouth,float cnorth,float inc){
		increSToN=csouth<cnorth;
		southVal =csouth;
		northVal =cnorth;
		number   =(int)Math.round((northVal-southVal)/inc)+1;
		meanIntv =inc;
		
		if(number-1!=Math.round((northVal-southVal)/inc))
		throw new IllegalArgumentException("overflow");
		
		if( increSToN&&inc<0) throw new IllegalArgumentException("contour increment ("+inc+") should be positive");
		if(!increSToN&&inc>0) throw new IllegalArgumentException("contour increment ("+inc+") should be negative");
		
		check();
		
		values=new double[number];
		for(int c=0,I=number;c<I;c++) values[c]=southVal+inc*c;
		
		isLinear=true;
		
		ct=new ContourTable();
	}
	
	/**
	 * Compute the contours by specifying all the values from south to north
	 * 
	 * @param	cs	contours from south to north
	 */
	public Contours(float[] cs){
		increSToN=cs[0]<cs[cs.length-1];
		southVal =cs[0];
		northVal =cs[cs.length-1];
		number   =cs.length;
		meanIntv =(northVal-southVal)/(number-1.0);
		
		check();
		
		values=new double[number];
		for(int i=0,I=values.length;i<I;i++) values[i]=cs[i];
		
		for(int i=0,I=number-1;i<I;i++)
		if(Math.abs(values[i+1]-values[i]-meanIntv)/meanIntv>1e-3){ isLinear=false; break;}
		
		ct=new ContourTable();
	}
	
	/**
	 * Compute the contours by specifying all the values from south to north
	 * 
	 * @param	cs	contours from south to north
	 */
	public Contours(double[] cs){
		increSToN=cs[0]<cs[cs.length-1];
		southVal =cs[0];
		northVal =cs[cs.length-1];
		number   =cs.length;
		meanIntv =(northVal-southVal)/(number-1.0);
		
		check();
		
		values=cs.clone();
		
		for(int i=0,I=number-1;i<I;i++)
		if(Math.abs(values[i+1]-values[i]-meanIntv)/meanIntv>1e-3){ isLinear=false; break;}
		
		ct=new ContourTable();
	}
	
	
	/*** getor and setor ***/
	public boolean increaseSToN(){ return increSToN;}
	
	public boolean isLinear(){ return isLinear;}
	
	public int getContourNumber(){ return number;}
	
	public double getMaxValue(){ return increSToN?northVal:southVal;}
	
	public double getMinValue(){ return increSToN?southVal:northVal;}
	
	public double[] getValues(){ return values;}
	
	public double[] getMappedEquivalentYs(){ return ct.Yeqvs;}
	
	public double[] getMappedAreas(){ return ct.areas;}
	
	public double[] getMappedMass(){ return ct.mass;}
	
	public double[] getMappedCirculation(){ return ct.circu;}
	
	public void setContours(double[] vs){ this.values=vs;}
	
	public void setAreas(double[] areas){ ct.areas=areas;}
	
	public void setMass(double[] mass){ ct.mass=mass;}
	
	public void setCirculation(double[] circu){ ct.circu=circu;}
	
	public void setYEs(double[] Ye){ ct.Yeqvs=Ye;}
	
	
	/**
	 * Find a value of M corresponding to a value (qv) in q by linear interpolating the tabulated M(q).
	 * 
	 * @param	qv	a given contour value
	 */
	public double findArea(double qv){ return ct.findArrayValue(qv,values,ct.areas);}
	
	public double findMass(double qv){ return ct.findArrayValue(qv,values,ct.mass);}
	
	public double findCirculation(double qv){ return ct.findArrayValue(qv,values,ct.circu);}
	
	public double findYeq(double qv){ return ct.findArrayValue(qv,values,ct.Yeqvs);}
	
	public double findContourByMass(double M){ return ct.findArrayValue(M,ct.mass,values);}
	
	public double findContourByYeq(double Yeq){ return ct.findArrayValue(Yeq,ct.Yeqvs,values);}
	
	
	/**
	 * Find values of M corresponding to values of q (qvs) by linear interpolating the tabulated M(q).
	 * 
	 * @param	qvs	contour values
	 */
	public double[] findAreas(double[] qvs){ return ct.findArrayValues(qvs,values,ct.areas);}
	
	public double[] findMasses(double[] qvs){ return ct.findArrayValues(qvs,values,ct.mass);}
	
	public double[] findCirculations(double[] qvs){ return ct.findArrayValues(qvs,values,ct.circu);}
	
	public double[] findYeqs(double[] qvs){ return ct.findArrayValues(qvs,values,ct.Yeqvs);}
	
	public double[] findContoursByMass(double[] Ms){ return ct.findArrayValues(Ms,ct.mass,values);}
	
	public double[] findContoursByYeq(double[] Yeqs){ return ct.findArrayValues(Yeqs,ct.Yeqvs,values);}
	
	
	/**
     * used to print out
     */
	public String toString(){
		if(isLinear)
			return "contours are linear "+"["+southVal+" : " +meanIntv+" : "+northVal+"]";
		else
			return "contours are levels "+"["+southVal+" : ~"+meanIntv+" : "+northVal+"]";
	}
	
	
	/*** helper methods ***/
	private void check(){
		if(number<2)
		throw new IllegalArgumentException("the number of contour is at least 2");
		if(meanIntv==0)
		throw new IllegalArgumentException("mean increment is 0 for "+toString());
	}
	
	
	private final class ContourTable{
		//
		double[] areas=null;	// area of each contour
		double[] mass =null;	// mass within each contour
		double[] Yeqvs=null;	// equivalent Ys
		double[] circu=null;	// circulation
		
		
		/**
		 * Constructor.
		 */
		public ContourTable(){}
		
		
		/**
		 * Find a value of M corresponding to a value (qv) in q by linear interpolating the tabulated M(q).
		 * 
		 * @param	qv	a given value in array q
		 * @param	q	a given q array
		 * @param	M	a given M array
		 */
		public double findArrayValue(double qv,double[] q,double[] M){
			int len=q.length;
			
			if(len!=M.length) throw new IllegalArgumentException("lengths not equal");
			
			int idx=q[len-1]>q[0]?ArrayUtil.getLEIdxIncre(q,qv):ArrayUtil.getLEIdxDecre(q,qv);
			
			if(idx==-1) throw new IllegalArgumentException("value of "+qv+" is out of range ["+q[0]+", "+q[q.length-1]+"]");
			
			if(idx==len-1) return M[idx];
			
			return InterpolationModel.linearInterpolation(q[idx],q[idx+1],M[idx],M[idx+1],qv);
		}
		
		/**
		 * Find values of M corresponding to values of q (qvs) by linear interpolating the tabulated M(q).
		 * 
		 * @param	cntr	a given contour value
		 */
		public double[] findArrayValues(double[] qvs,double[] q,double[] M){
			int len=qvs.length;
			
			double[] re=new double[len];
			
			for(int i=0;i<len;i++) re[i]=findArrayValue(qvs[i],q,M);
			
			return re;
		}
		
	}
	
	
	/*** test
	public static void main(String[] args){
		
	}*/
}
