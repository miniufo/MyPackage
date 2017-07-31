/**
 * @(#)Variable.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import miniufo.basic.Operatable;
import miniufo.basic.InterpolationModel;
import miniufo.basic.InterpolationModel.Type;


/**
 * 4-dimensional variable for weather analysis
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Variable implements Operatable<Variable>{
	//
	private int t=0;
	private int z=0;
	private int y=0;
	private int x=0;
	
	private float undef		  =Float.NaN;
	
	private boolean has_undef =false;
	private boolean tfirst    =true;
	
	private String unit       =null;
	private String vname      =null;
	private String comment    =null;
	
	private Range range		  =null;
	
	private float[][][][] data=null;
	
	private static final Pattern unitPtn=Pattern.compile("\\(([^\\(\\).]*?)\\)$");
	
	public static enum Dimension{X,Y,Z,T}
	
	
	/**
     * constructor
     *
     * @param	vname	variable name
     * @param	tfirst	whether t dimension is the first dimension
     * @param	data	the storage
     */
	public Variable(String vname,boolean tfirst,float[][][][] data){
		if(tfirst){
			t=data.length;
			z=data[0].length;
			y=data[0][0].length;
			x=data[0][0][0].length;
			this.data=data;
			
		}else{
			z=data.length;
			y=data[0].length;
			x=data[0][0].length;
			t=data[0][0][0].length;
			this.data=data;
		}
		
		this.tfirst=tfirst;
		this.vname =vname;
		this.range =new Range(t,z,y,x);
	}
	
	/**
     * constructor
     *
     * @param	vname	variable name
     * @param	tfirst	whether t dimension is the first dimension
     * @param	r		data range
     */
	public Variable(String vname,boolean tfirst,Range range){
		this.vname =vname;
		this.range =range;
		this.tfirst=tfirst;
		
		t=range.getTRange()[2];	z=range.getZRange()[2];
		y=range.getYRange()[2];	x=range.getXRange()[2];
		
		if(tfirst) data=new float[t][z][y][x];
		else	   data=new float[z][y][x][t];
	}
	
	/**
     * constructor
     *
     * @param	vname	variable name
     * @param	r		data range
     */
	public Variable(String vname,Range range){
     	t=range.getTRange()[2];	z=range.getZRange()[2];
     	y=range.getYRange()[2];	x=range.getXRange()[2];
     	
     	this.range =range;
     	this.vname =vname;
     	this.tfirst=(t<=z);
     	
		if(tfirst) data=new float[t][z][y][x];
		else	   data=new float[z][y][x][t];
     }
     
	/**
     * constructor, construct according to the given variable
     *
     * @param	vname	variable name
     * @param	v		another variable
     */
	public Variable(String vname,Variable v){
		this.vname  =vname;
		this.tfirst =v.tfirst;
		this.range  =(Range)(v.range.clone());
		
		if(v.unit   !=null) this.unit   =v.unit;
		if(v.comment!=null) this.comment=v.comment;
		
		t=v.t;	z=v.z;	y=v.y;	x=v.x;
		
		undef=v.undef;	has_undef=v.has_undef;
		
		if(tfirst) data=new float[t][z][y][x];
		else	   data=new float[z][y][x][t];
	}
	
	
	/*** private constructor ***/
	private Variable(){}
	
	
	/**
     * joint variables into one without memory allocation
     * caution should be taken in the side effect
     * tFirst should be true
     *
     * @param	v	variables in an array
     */
	public static Variable joint(Variable... v){
		if(!v[0].isTFirst()) throw new IllegalArgumentException("variables should be t-first type");
		
		for(int i=1;i<v.length;i++)
		if(!v[0].isLike(v[i])) throw new IllegalArgumentException("dimension not same");
		
		Variable r=new Variable();
		
		r.t=v[0].t*v.length;
		r.z=v[0].z;
		r.y=v[0].y;
		r.x=v[0].x;

		r.tfirst	=true;
		r.undef		=v[0].undef;
		r.range		=(Range)(v[0].range.clone());
		r.has_undef	=v[0].has_undef;
		
		if(v[0].unit   !=null) r.unit   =v[0].unit   ;
		if(v[0].vname  !=null) r.vname  =v[0].vname  ;
		if(v[0].comment!=null) r.comment=v[0].comment;
		
		r.range.getTRange()[2]=r.t;	r.range.getTRange()[1]=r.range.getTRange()[0]-1+r.t;
		
		r.data=new float[r.t][][][];
		
		for(int l=0;l<v.length;l++)
		for(int ll=0;ll<v[l].t;ll++) r.data[v[0].t*l+ll]=v[l].data[ll];
		
		return r;
	}
	
	
	/*** getor and setor ***/
	public int getXCount(){ return x;}
	
	public int getYCount(){ return y;}
	
	public int getZCount(){ return z;}
	
	public int getTCount(){ return t;}
	
	public float getUndef(){ return undef;}
	
	public boolean isTFirst(){ return tfirst;}
	
	public Range  getRange() { return range; }
	
	public String getUnit(){ return unit; }
	
	public String getName(){ return vname;}
	
	public String getComment(){ return comment;}
	
	public String getCommentAndUnit(){
		if(unit!=null&&!unit.equals("")) return comment+" ("+unit+")";
		else return comment;
	}
	
	public float[][][][] getData(){ return data;}
	
	public void setRange(Range r){
		t=r.getTRange()[2];	z=r.getZRange()[2];
		y=r.getYRange()[2];	x=r.getXRange()[2];
		
		range.setTRange(r);	range.setZRange(r);
		range.setYRange(r);	range.setXRange(r);
	}
	
	public void setUndef(float undef){ this.undef=undef; has_undef=true;}
	
	public void setUnit(String unit) { this.unit =unit ;}
	
	public void setName(String vname){ this.vname=vname;}
	
	public void setComment(String comment){ this.comment=comment;}
	
	public void setCommentAndUnit(String commentAndUnit){
		Matcher m=unitPtn.matcher(commentAndUnit);
		
		if(m.find()){
			String unitInBracket=m.group();
			this.unit=unitInBracket.substring(1,unitInBracket.length()-1);
			this.comment=m.replaceAll("").trim();
			
		}else this.comment=commentAndUnit;
	}
	
	public void setData(float[][][][] data,boolean tfirst){
		this.data=data;
		
		if(tfirst){
			t=data.length;			z=data[0].length;
			y=data[0][0].length;	x=data[0][0][0].length;
			
		}else{
			t=data[0][0][0].length;	z=data.length;
			y=data[0].length;		x=data[0][0].length;
		}
		
		this.tfirst=tfirst;	int[] rr=null;
		
		rr=range.getTRange();	rr[2]=t;	rr[1]=rr[2]+rr[0]-1;
		rr=range.getZRange();	rr[2]=z;	rr[1]=rr[2]+rr[0]-1;
		rr=range.getYRange();	rr[2]=y;	rr[1]=rr[2]+rr[0]-1;
		rr=range.getXRange();	rr[2]=x;	rr[1]=rr[2]+rr[0]-1;
	}
	
	
	/**
     * plus a constant to the variable
     *
     * @param	v	a given constant
     */
	public Variable plus(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=data[l][k][j][i]+v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=data[k][j][i][l]+v;
		}
		
		return re;
	}
	
	public Variable plusEq(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(data[l][k][j][i]!=undef) data[l][k][j][i]+=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(data[k][j][i][l]!=undef) data[k][j][i][l]+=v;
		}
		
		return this;
	}
	
	/**
     * minus a constant from the variable
     *
     * @param	v	a given constant
     */
	public Variable minus(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=data[l][k][j][i]-v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=data[k][j][i][l]-v;
		}
		
		return re;
	}
	
	public Variable minusEq(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]-=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]-=v;
		}
		
		return this;
	}
	
	/**
     * multiply with a constant
     *
     * @param	v	a given constant
     */
	public Variable multiply(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=data[l][k][j][i]*v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=data[k][j][i][l]*v;
		}
		
		return re;
	}
	
	public Variable multiplyEq(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]*=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]*=v;
		}
		
		return this;
	}
	
	/**
     * divided by a constant
     *
     * @param	v	a given constant
     */
	public Variable divide(float v){
		if(v==0) throw new IllegalArgumentException("divided by zero");
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=data[l][k][j][i]/v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=data[k][j][i][l]/v;
		}
		
		return re;
	}
	
	public Variable divideEq(float v){
		if(v==0) throw new IllegalArgumentException("divided by zero");
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]/=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]/=v;
		}
		
		return this;
	}
	
	/**
     * power of a given float
     *
     * @param	v	a given constant
     */
	public Variable pow(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=(float)Math.pow(data[l][k][j][i],v);
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=(float)Math.pow(data[k][j][i][l],v);
		}
		
		return re;
	}
	
	public Variable powEq(float v){
		if(v==undef) System.out.println(" Warning: parameter equals undefined value");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]=(float)Math.pow(data[l][k][j][i],v);
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]=(float)Math.pow(data[k][j][i][l],v);
		}
		
		return this;
	}
	
	/**
     * exp
     */
	public Variable exp(){
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=(float)Math.exp(data[l][k][j][i]);
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=(float)Math.exp(data[k][j][i][l]);
		}
		
		return re;
	}
	
	public Variable expEq(){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]=(float)Math.exp(data[l][k][j][i]);
			else data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]=(float)Math.exp(data[k][j][i][l]);
			else data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	/**
     * abs
     */
	public Variable abs(){
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=Math.abs(data[l][k][j][i]);
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=Math.abs(data[k][j][i][l]);
		}
		
		return re;
	}
	
	public Variable absEq(){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]=Math.abs(data[l][k][j][i]);
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]=Math.abs(data[k][j][i][l]);
		}
		
		return this;
	}
	
	/**
     * square
     */
	public Variable square(){
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=data[l][k][j][i]*data[l][k][j][i];
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=data[k][j][i][l]*data[k][j][i][l];
		}
		
		return re;
	}
	
	public Variable squareEq(){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]*=data[l][k][j][i];
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]*=data[k][j][i][l];
		}
		
		return this;
	}
	
	/**
     * sqrt
     */
	public Variable sqrt(){
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=(float)(Math.sqrt(data[l][k][j][i]));
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=(float)(Math.sqrt(data[k][j][i][l]));
		}
		
		return re;
	}
	
	public Variable sqrtEq(){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]=(float)(Math.sqrt(data[l][k][j][i]));
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]=(float)(Math.sqrt(data[k][j][i][l]));
		}
		
		return this;
	}
	
	/**
     * reciprocal
     */
	public Variable reciprocal(){
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=1/data[l][k][j][i];
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=1/data[k][j][i][l];
		}
		
		return re;
	}
	
	public Variable reciprocalEq(){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]=1/data[l][k][j][i];
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]=1/data[k][j][i][l];
		}
		
		return this;
	}
	
	/**
     * logarithm (base e)
     */
	public Variable logarithm(){
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) rdata[l][k][j][i]=(float)(Math.log(data[l][k][j][i]));
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) rdata[k][j][i][l]=(float)(Math.log(data[k][j][i][l]));
		}
		
		return re;
	}
	
	public Variable logarithmEq(){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef) data[l][k][j][i]=(float)(Math.log(data[l][k][j][i]));
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef) data[k][j][i][l]=(float)(Math.log(data[k][j][i][l]));
		}
		
		return this;
	}
	
	
	/**
     * calculate hypotenuse
     */
	public Variable hypotenuse(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		Variable r=copy();
		
		float[][][][] rdata=r.getData();
		float[][][][] vdata=v.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef)
				rdata[l][k][j][i]=(float)(Math.hypot(data[l][k][j][i],vdata[l][k][j][i]));
			else
				rdata[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef)
				rdata[k][j][i][l]=(float)(Math.hypot(data[k][j][i][l],vdata[k][j][i][l]));
			else
				rdata[k][j][i][l]=undef;
		}
		
		return r;
	}
	
	public Variable hypotenuseEq(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		float[][][][] vdata=v.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef)
				data[l][k][j][i]=(float)(Math.hypot(data[l][k][j][i],vdata[l][k][j][i]));
			else
				data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef)
				data[k][j][i][l]=(float)(Math.hypot(data[k][j][i][l],vdata[k][j][i][l]));
			else
				data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	
	/**
     * plus a variable point to point
     *
     * @param	v	a given variable
     *
     * @exception	if this variable is not dimensionally the same as the given variable
     */
	public Variable plus(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[l][k][j][i]!=undef&&data[l][k][j][i]!=undef)
				rdata[l][k][j][i]=v.data[l][k][j][i]+data[l][k][j][i];
			else
				rdata[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[k][j][i][l]!=undef&&data[k][j][i][l]!=undef)
				rdata[k][j][i][l]=v.data[k][j][i][l]+data[k][j][i][l];
			else
				rdata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	public Variable plusEq(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[l][k][j][i]!=undef&&data[l][k][j][i]!=undef) data[l][k][j][i]+=v.data[l][k][j][i];
			else data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[k][j][i][l]!=undef&&data[k][j][i][l]!=undef) data[k][j][i][l]+=v.data[k][j][i][l];
			else data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	/**
     * minus a variable point to point
     *
     * @param	v	a given variable
     *
     * @exception	if this variable is not dimensionally the same as the given variable
     */
	public Variable minus(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[l][k][j][i]!=undef&&data[l][k][j][i]!=undef)
				rdata[l][k][j][i]=data[l][k][j][i]-v.data[l][k][j][i];
			else
				rdata[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[k][j][i][l]!=undef&&data[k][j][i][l]!=undef)
				rdata[k][j][i][l]=data[k][j][i][l]-v.data[k][j][i][l];
			else
				rdata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	public Variable minusEq(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[l][k][j][i]!=undef&&data[l][k][j][i]!=undef) data[l][k][j][i]-=v.data[l][k][j][i];
			else data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(v.data[k][j][i][l]!=undef&&data[k][j][i][l]!=undef) data[k][j][i][l]-=v.data[k][j][i][l];
			else data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	/**
     * multiply with a variable point to point
     *
     * @param	v	a given variable
     *
     * @exception	if this variable is not dimensionally the same as the given variable
     */
	public Variable multiply(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&v.data[l][k][j][i]!=undef)
				rdata[l][k][j][i]=data[l][k][j][i]*v.data[l][k][j][i];
			else
				rdata[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&v.data[k][j][i][l]!=undef)
				rdata[k][j][i][l]=data[k][j][i][l]*v.data[k][j][i][l];
			else
				rdata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	public Variable multiplyEq(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&v.data[l][k][j][i]!=undef) data[l][k][j][i]*=v.data[l][k][j][i];
			else data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&v.data[k][j][i][l]!=undef) data[k][j][i][l]*=v.data[k][j][i][l];
			else data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	/**
     * divided by a variable point to point
     *
     * @param	v	a given variable
     *
     * @exception	if this variable is not dimensionally the same as the given variable
     */
	public Variable divide(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&v.data[l][k][j][i]!=undef)
				rdata[l][k][j][i]=data[l][k][j][i]/v.data[l][k][j][i];
			else
				rdata[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&v.data[k][j][i][l]!=undef)
				rdata[k][j][i][l]=data[k][j][i][l]/v.data[k][j][i][l];
			else
				rdata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	public Variable divideEq(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&v.data[l][k][j][i]!=undef) data[l][k][j][i]/=v.data[l][k][j][i];
			else data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&v.data[k][j][i][l]!=undef) data[k][j][i][l]/=v.data[k][j][i][l];
			else data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	/**
     * power of a given Variable
     *
     * @param	v	a given Variable
     */
	public Variable pow(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		Variable re=copy();
		
		float[][][][] rdata=re.getData();
		float[][][][] vdata=v.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef)
				rdata[l][k][j][i]=(float)Math.pow(data[l][k][j][i],vdata[l][k][j][i]);
			else
				rdata[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef)
				rdata[k][j][i][l]=(float)Math.pow(data[k][j][i][l],vdata[k][j][i][l]);
			else
				rdata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	public Variable powEq(Variable v){
		if(!isLike(v)) throw new IllegalArgumentException("dimension not same");
		
		float[][][][] vdata=v.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef)
				data[l][k][j][i]=(float)Math.pow(data[l][k][j][i],vdata[l][k][j][i]);
			else
				data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef)
				data[k][j][i][l]=(float)Math.pow(data[k][j][i][l],vdata[k][j][i][l]);
			else
				data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	
	/**
     * set the data to a given value
     */
	public Variable setInner(float v){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++) data[l][k][j][i]=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++) data[k][j][i][l]=v;
		}
		
		return this;
	}
	
	public Variable setOuter(float v){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){ data[l][k][j][0]=v; data[l][k][j][y-1]=v;}
				for(int i=0;i<x;i++){ data[l][k][0][i]=v; data[l][k][x-1][i]=v;}
			}
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){ data[k][j][0][l]=v; data[k][j][y-1][l]=v;}
				for(int i=0;i<x;i++){ data[k][0][i][l]=v; data[k][x-1][i][l]=v;}
			}
		}
		
		return this;
	}
	
	public Variable setValue(float v){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) data[l][k][j][i]=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) data[k][j][i][l]=v;
		}
		
		return this;
	}
	
	/**
     * if data < v then set v to data
     */
	public Variable setDataMin(float v){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(data[l][k][j][i]<v) data[l][k][j][i]=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(data[k][j][i][l]<v) data[k][j][i][l]=v;
		}
		
		return this;
	}
	
	/**
     * if data > v then set v to data
     */
	public Variable setDataMax(float v){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(data[l][k][j][i]>v) data[l][k][j][i]=v;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(data[k][j][i][l]>v) data[k][j][i][l]=v;
		}
		
		return this;
	}
	
	/**
     * replace the undef to new undef
     */
	public Variable replaceUndefData(float newUndef){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[l][k][j][i]==undef||Float.isNaN(data[l][k][j][i])) data[l][k][j][i]=newUndef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(data[k][j][i][l]==undef||Float.isNaN(data[k][j][i][l])) data[k][j][i][l]=newUndef;
		}
		
		setUndef(newUndef);
		
		return this;
	}
	
	/**
     * change the NaN values to undefined values
     */
	public Variable changeNaNToUndef(){
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Float.isNaN(data[l][k][j][i])) data[l][k][j][i]=undef;
					
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Float.isNaN(data[k][j][i][l])) data[k][j][i][l]=undef;
		}
		
		return this;
	}
	
	
	/**
     * whether the dimension of the variable is the same with the given variable
     *
     * @param	v	a given variable
     *
     * @return	is the same or not
     */
	public boolean isLike(Variable v){
		if(x!=v.getXCount()||y!=v.getYCount()||z!=v.getZCount()||t!=v.getTCount()) return false;
		if(tfirst!=v.tfirst) return false;
		return true;
	}
	
	/**
     * whether the area of the variable is the same with the given variable
     *
     * @param	v	a given variable
     *
     * @return	is the same or not
     */
	public boolean isAreaLike(Variable v){
		if(x!=v.getXCount()||y!=v.getYCount()) return false;
		if(tfirst!=v.tfirst) return false;
		return true;
	}
	
	
	/**
     * average along one specific dimension
     * from start to end
     * 
     * @param	d		dimension (e.g., t, z, y, x)
     * @param	str		start tag
     * @param	end		end tag
     * @param	anom	whether anomalize the original data (true for side-effect)
     *
     * @return	re	new variable after averaging
     */
	public Variable averageAlong(Dimension d,int str,int end,boolean anom){
		int tlen=t;	int zlen=z;
		int ylen=y;	int xlen=x;
		
		switch(d){
			case T: tlen=1; break;
			case Z: zlen=1; break;
			case Y: ylen=1; break;
			case X: xlen=1; break;
			default: throw new IllegalArgumentException("unsupported dimension");
		}
		
		Range nr=new Range(tlen,zlen,ylen,xlen);
		
		Variable nv=new Variable(vname,tfirst,nr);	nv.setUndef(undef);
		if( unit  !=null) nv.setUnit(unit);
		if(comment!=null) nv.setComment(comment);
		
		float[][][][] nvdata=nv.getData();
		
		if(tfirst){
			switch(d){
			case T:
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					int count=0;	float sum=0;
					
					for(int l=str;l<=end;l++) if(data[l][k][j][i]!=undef){ sum+=data[l][k][j][i]; count++;}
					
					if(count!=0) nvdata[0][k][j][i]=sum/count;
					else nvdata[0][k][j][i]=undef;
					
					if(anom) for(int l=str;l<=end;l++)
					if(data[l][k][j][i]!=undef) data[l][k][j][i]-=nvdata[0][k][j][i];
				}
				
				nr.setTRange(range.getTRange()[0]);
				nr.setZRange(range);
				nr.setYRange(range);
				nr.setXRange(range);
				
				break;
				
			case Z:
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					int count=0;	float sum=0;
					
					for(int k=str;k<=end;k++) if(data[l][k][j][i]!=undef){ sum+=data[l][k][j][i]; count++;}
					
					if(count!=0) nvdata[l][0][j][i]=sum/count;
					else nvdata[l][0][j][i]=undef;
					
					if(anom) for(int k=str;k<=end;k++)
					if(data[l][k][j][i]!=undef) data[l][k][j][i]-=nvdata[l][0][j][i];
				}
				
				nr.setTRange(range);
				nr.setZRange(range.getZRange()[0]);
				nr.setYRange(range);
				nr.setXRange(range);
				
				break;
				
			case Y:
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					int count=0;	float sum=0;
					
					for(int j=str;j<=end;j++) if(data[l][k][j][i]!=undef){ sum+=data[l][k][j][i]; count++;}
					
					if(count!=0) nvdata[l][k][0][i]=sum/count;
					else nvdata[l][k][0][i]=undef;
					
					if(anom) for(int j=str;j<=end;j++)
					if(data[l][k][j][i]!=undef) data[l][k][j][i]-=nvdata[l][k][0][i];
				}
				
				nr.setTRange(range);
				nr.setZRange(range);
				nr.setYRange(range.getYRange()[0]);
				nr.setXRange(range);
				
				break;
				
			case X:
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;	float sum=0;
					
					for(int i=str;i<=end;i++) if(data[l][k][j][i]!=undef){ sum+=data[l][k][j][i]; count++;}
					
					if(count!=0) nvdata[l][k][j][0]=sum/count;
					else nvdata[l][k][j][0]=undef;
					
					if(anom) for(int i=str;i<=end;i++)
					if(data[l][k][j][i]!=undef) data[l][k][j][i]-=nvdata[l][k][j][0];
				}
				
				nr.setTRange(range);
				nr.setZRange(range);
				nr.setYRange(range);
				nr.setXRange(range.getXRange()[0]);
				
				break;
				
			default: throw new IllegalArgumentException("unsupported dimension");
			}
			
		}else{
			switch(d){
			case T:
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					int count=0;	float sum=0;
					
					for(int l=str;l<=end;l++) if(data[k][j][i][l]!=undef){ sum+=data[k][j][i][l]; count++;}
					
					if(count!=0) nvdata[k][j][i][0]=sum/count;
					else nvdata[k][j][i][0]=undef;
					
					if(anom) for(int l=str;l<=end;l++)
					if(data[k][j][i][l]!=undef) data[k][j][i][l]-=nvdata[k][j][i][0];
				}
				
				nr.setTRange(range.getTRange()[0]);
				nr.setZRange(range);
				nr.setYRange(range);
				nr.setXRange(range);
				
				break;
				
			case Z:
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					int count=0;	float sum=0;
					
					for(int k=str;k<=end;k++) if(data[k][j][i][l]!=undef){ sum+=data[k][j][i][l]; count++;}
					
					if(count!=0) nvdata[0][j][i][l]=sum/count;
					else nvdata[0][j][i][l]=undef;
					
					if(anom) for(int k=str;k<=end;k++)
					if(data[k][j][i][l]!=undef) data[k][j][i][l]-=nvdata[0][j][i][l];
				}
				
				nr.setTRange(range);
				nr.setZRange(range.getZRange()[0]);
				nr.setYRange(range);
				nr.setXRange(range);
				
				break;
				
			case Y:
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					int count=0;	float sum=0;
					
					for(int j=str;j<=end;j++) if(data[k][j][i][l]!=undef){ sum+=data[k][j][i][l]; count++;}
					
					if(count!=0) nvdata[k][0][i][l]=sum/count;
					else nvdata[k][0][i][l]=undef;
					
					if(anom) for(int j=str;j<=end;j++)
					if(data[k][j][i][l]!=undef) data[k][j][i][l]-=nvdata[k][0][i][l];
				}
				
				nr.setTRange(range);
				nr.setZRange(range);
				nr.setYRange(range.getYRange()[0]);
				nr.setXRange(range);
				
				break;
				
			case X:
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;	float sum=0;
					
					for(int i=str;i<=end;i++) if(data[k][j][i][l]!=undef){ sum+=data[k][j][i][l]; count++;}
					
					if(count!=0) nvdata[k][j][0][l]=sum/count;
					else nvdata[k][j][0][l]=undef;
					
					if(anom) for(int i=str;i<=end;i++)
					if(data[k][j][i][l]!=undef) data[k][j][i][l]-=nvdata[k][j][0][l];
				}
				
				nr.setTRange(range);
				nr.setZRange(range);
				nr.setYRange(range);
				nr.setXRange(range.getXRange()[0]);
				
				break;
				
			default: throw new IllegalArgumentException("unsupported dimension");
			}
		}
		
		return nv;
	}
	
	public Variable averageAlong(Dimension d,int str,int end){ return averageAlong(d,str,end,false);}
	
	public Variable averageAlong(Dimension d){
		switch(d){
			case T: return averageAlong(d,0,t-1,false);
			case Z: return averageAlong(d,0,z-1,false);
			case Y: return averageAlong(d,0,y-1,false);
			case X: return averageAlong(d,0,x-1,false);
			default: throw new IllegalArgumentException("unsupported dimension");
		}
	}
	
	
	/**
     * anomalize in t-direction and return mean value in t-direction
     *
     * @return	mean value in t-direction
     */
	public Variable anomalizeT(){ return averageAlong(Dimension.T,0,t-1,true);}
	
	/**
     * anomalize in z-direction and return mean value in z-direction
     *
     * @return	mean value in z-direction
     */
	public Variable anomalizeZ(){ return averageAlong(Dimension.Z,0,z-1,true);}
	
	/**
     * anomalize in y-direction and return mean value in y-direction
     *
     * @return	mean value in y-direction
     */
	public Variable anomalizeY(){ return averageAlong(Dimension.Y,0,y-1,true);}
	
	/**
     * anomalize in x-direction and return mean value in x-direction
     *
     * @return	mean value in x-direction
     */
	public Variable anomalizeX(){ return averageAlong(Dimension.X,0,x-1,true);}
	
	/**
     * anomalize in YX-plane and return mean value
     *
     * @return	mean value of area
     */
	public Variable anomalizeYX(){
		Range nr=new Range(t,z,1,1);
		
		Variable nv=new Variable(vname,tfirst,nr);	nv.setUndef(undef);
		if( unit  !=null) nv.setUnit(unit);
		if(comment!=null) nv.setCommentAndUnit(comment);
		
		float[][][][] nvdata=nv.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;	float sum=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) if(data[l][k][j][i]!=undef){ sum+=data[l][k][j][i]; count++;}
					
				if(count!=0) nvdata[l][k][0][0]=sum/count;
				else nvdata[l][k][0][0]=undef;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) if(data[l][k][j][i]!=undef) data[l][k][j][i]-=nvdata[l][k][0][0];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;	float sum=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) if(data[k][j][i][l]!=undef){ sum+=data[k][j][i][l]; count++;}
					
				if(count!=0) nvdata[k][0][0][l]=sum/count;
				else nvdata[k][0][0][l]=undef;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) if(data[k][j][i][l]!=undef) data[k][j][i][l]-=nvdata[k][0][0][l];
			}
		}
		
		nr.setTRange(range);	nr.setZRange(range);
		nr.getYRange()[0]=nr.getYRange()[1]=range.getYRange()[0];
		nr.getXRange()[0]=nr.getXRange()[1]=range.getXRange()[0];
		
		return nv;
	}
	
	/**
     * anomalize in ZY-plane and return mean value
     *
     * @return	mean value of area
     */
	public Variable anomalizeZY(){
		Range nr=new Range(t,1,1,x);
		
		Variable nv=new Variable(vname,tfirst,nr);	nv.setUndef(undef);
		if( unit  !=null) nv.setUnit(unit);
		if(comment!=null) nv.setCommentAndUnit(comment);
		
		float[][][][] nvdata=nv.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++){
				int count=0;	float sum=0;
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) if(data[l][k][j][i]!=undef){ sum+=data[l][k][j][i]; count++;}
				
				if(count!=0) nvdata[l][0][0][i]=sum/count;
				else nvdata[l][0][0][i]=undef;
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) if(data[l][k][j][i]!=undef) data[l][k][j][i]-=nvdata[l][0][0][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++){
				int count=0;	float sum=0;
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) if(data[k][j][i][l]!=undef){ sum+=data[k][j][i][l]; count++;}
				
				if(count!=0) nvdata[0][0][i][l]=sum/count;
				else nvdata[0][0][i][l]=undef;
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) if(data[k][j][i][l]!=undef) data[k][j][i][l]-=nvdata[0][0][i][l];
			}
		}
		
		nr.setTRange(range);	nr.setXRange(range);
		nr.getYRange()[0]=nr.getYRange()[1]=range.getYRange()[0];
		nr.getZRange()[0]=nr.getZRange()[1]=range.getZRange()[0];
		
		return nv;
	}
	
	/**
     * anomalize in ZX-plane and return mean value
     *
     * @return	mean value of area
     */
	public Variable anomalizeZX(){
		Range nr=new Range(t,1,y,1);
		
		Variable nv=new Variable(vname,tfirst,nr);	nv.setUndef(undef);
		if( unit  !=null) nv.setUnit(unit);
		if(comment!=null) nv.setCommentAndUnit(comment);
		
		float[][][][] nvdata=nv.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				int count=0;	float sum=0;
				
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++) if(data[l][k][j][i]!=undef){ sum+=data[l][k][j][i]; count++;}
				
				if(count!=0) nvdata[l][0][j][0]=sum/count;
				else nvdata[l][0][j][0]=undef;
				
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++) if(data[l][k][j][i]!=undef) data[l][k][j][i]-=nvdata[l][0][j][0];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				int count=0;	float sum=0;
				
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++) if(data[k][j][i][l]!=undef){ sum+=data[k][j][i][l]; count++;}
				
				if(count!=0) nvdata[0][j][0][l]=sum/count;
				else nvdata[0][j][0][l]=undef;
				
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++) if(data[k][j][i][l]!=undef) data[k][j][i][l]-=nvdata[0][j][0][l];
			}
		}
		
		nr.setTRange(range);	nr.setYRange(range);
		nr.getXRange()[0]=nr.getXRange()[1]=range.getXRange()[0];
		nr.getZRange()[0]=nr.getZRange()[1]=range.getZRange()[0];
		
		return nv;
	}
	
	
	/**
     * randomly sorted the T-data
     *
     * @param	v	a given Variable
     */
	public Variable randomT(){
		float[][] rand=new float[t][2];
		
		for(int l=0;l<t;l++){ rand[l][0]=(float)Math.random(); rand[l][1]=l;}
		
		java.util.Arrays.sort(rand,new miniufo.basic.MultiArrayComparator(0));
		
		if(tfirst){
			float[][][][] tmp=new float[t][][][];
			
			for(int l=0;l<t;l++)  tmp[l]=data[l];
			for(int l=0;l<t;l++) data[l]=tmp[(int)rand[l][1]];
			
		}else{
			float[] tmp=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				System.arraycopy(data[k][j][i],0,tmp,0,t);
				
				for(int l=0;l<t;l++) data[k][j][i][l]=tmp[(int)rand[l][1]];
			}
		}
		
		return this;
	}
	
	
	/**
     * interpolation along each dimension
     * 
     * @param	n		length in each dimension after interpolation
     * @param	type	type of interpolation
     */
	public Variable interpolateT(int n,Type type){
		if(n<2) throw new IllegalArgumentException("interpolate length should be at least 2");
		
		if(n==t) return copy();
		
		Variable res=new Variable(vname,tfirst,new Range(n,z,y,x));
		
		res.setComment(comment);
		res.setUnit(unit);
		res.setUndef(undef);
		
		float[][][][] rdata=res.getData();
		
		if(tfirst){
			float[] src=new float[t];
			float[] des=new float[n];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) src[l]=data[l][k][j][i];
				InterpolationModel.interp1D(src,des,type,undef);
				for(int l=0;l<n;l++) rdata[l][k][j][i]=des[l];
			}
			
		}else{
			float[] des=new float[n];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] src=data[k][j][i];
				InterpolationModel.interp1D(src,des,type,undef);
				System.arraycopy(des,0,rdata[k][j][i],0,n);
			}
		}
		
		Range nr=res.getRange();
		
		nr.setTRange(range); nr.getTRange()[2]=n;
		nr.setZRange(range);
		nr.setYRange(range);
		nr.setXRange(range);
		
		return res;
	}
	
	public Variable interpolateZ(int n,Type type){
		if(n==z) return copy();
		
		Variable res=new Variable(vname,tfirst,new Range(t,n,y,x));
		
		res.setComment(comment);
		res.setUnit(unit);
		res.setUndef(undef);
		
		float[][][][] rdata=res.getData();
		
		if(tfirst){
			float[] src=new float[z];
			float[] des=new float[n];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<z;k++) src[k]=data[l][k][j][i];
				InterpolationModel.interp1D(src,des,type);
				for(int k=0;k<n;k++) rdata[l][k][j][i]=des[k];
			}
			
		}else{
			float[] src=new float[z];
			float[] des=new float[n];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<z;k++) src[k]=rdata[k][j][i][l];
				InterpolationModel.interp1D(src,des,type);
				for(int k=0;k<n;k++) rdata[k][j][i][l]=des[k];
			}
		}
		
		Range nr=res.getRange();
		
		nr.setTRange(range);
		nr.setZRange(range); nr.getZRange()[2]=n;
		nr.setYRange(range);
		nr.setXRange(range);
		
		return res;
	}
	
	public Variable interpolateY(int n,Type type){
		if(n==y) return copy();
		
		Variable res=new Variable(vname,tfirst,new Range(t,z,n,x));
		
		res.setComment(comment);
		res.setUnit(unit);
		res.setUndef(undef);
		
		float[][][][] rdata=res.getData();
		
		if(tfirst){
			float[] src=new float[y];
			float[] des=new float[n];
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				for(int j=0;j<y;j++) src[j]=data[l][k][j][i];
				InterpolationModel.interp1D(src,des,type,undef);
				for(int j=0;j<n;j++) rdata[l][k][j][i]=des[j];
			}
			
		}else{
			float[] src=new float[z];
			float[] des=new float[n];
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				for(int j=0;j<y;j++) src[j]=rdata[k][j][i][l];
				InterpolationModel.interp1D(src,des,type,undef);
				for(int j=0;j<n;j++) rdata[k][j][i][l]=des[j];
			}
		}
		
		Range nr=res.getRange();
		
		nr.setTRange(range);
		nr.setZRange(range);
		nr.setYRange(range); nr.getYRange()[2]=n;
		nr.setXRange(range);
		
		return res;
	}
	
	public Variable interpolateX(int n,Type type){
		if(n==x) return copy();
		
		Variable res=new Variable(vname,tfirst,new Range(t,z,y,n));
		
		res.setComment(comment);
		res.setUnit(unit);
		res.setUndef(undef);
		
		float[][][][] rdata=res.getData();
		
		if(tfirst){
			float[] des=new float[n];
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				float[] src=data[l][k][j];
				InterpolationModel.interp1D(src,des,type,undef);
				System.arraycopy(des,0,rdata[l][k][j],0,n);
			}
			
		}else{
			float[] src=new float[x];
			float[] des=new float[n];
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				for(int i=0;i<x;i++) src[i]=rdata[k][j][i][l];
				InterpolationModel.interp1D(src,des,type,undef);
				for(int i=0;i<n;i++) rdata[k][j][i][l]=des[i];
			}
		}
		
		Range nr=res.getRange();
		
		nr.setTRange(range);
		nr.setZRange(range);
		nr.setYRange(range);
		nr.setXRange(range); nr.getXRange()[2]=n;
		
		return res;
	}
	
	/**
     * 2D interpolation
     * 
     * @param	yn		length in y-dimension after interpolation
     * @param	xn		length in x-dimension after interpolation
     * @param	ytype	type of interpolation in y-direction
     * @param	xtype	type of interpolation in x-direction
     */
	public Variable interpolateXY(int yn,int xn,Type ytype,Type xtype){
		if(xn==x&&yn==y) return copy();
		
		Variable res=new Variable(vname,tfirst,new Range(t,z,yn,xn));
		
		res.setComment(comment);
		res.setUnit(unit);
		res.setUndef(undef);
		
		float[][][][] rdata=res.getData();
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				float[][] src=data[l][k];
				InterpolationModel.interp2D(src,res.getData()[l][k],xtype,ytype,undef);
			}
			
		}else{
			float[][] src=new float[y ][x ];
			float[][] des=new float[yn][xn];
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) src[j][i]=rdata[k][j][i][l];
				InterpolationModel.interp2D(src,des,xtype,ytype,undef);
				for(int j=0;j<yn;j++)
				for(int i=0;i<xn;i++) rdata[k][j][i][l]=des[j][i];
			}
		}
		
		Range nr=res.getRange();
		
		nr.setTRange(range);
		nr.setZRange(range);
		nr.setYRange(range); nr.getYRange()[2]=yn;
		nr.setXRange(range); nr.getXRange()[2]=xn;
		
		return res;
	}
	
	
	/**
     * clone method
     */
	public Variable copy(){
	    Variable v=new Variable();
		
		v.x=x;	v.y=y;	v.z=z;	v.t=t;
		
		v.undef    =undef;
		v.tfirst   =tfirst;
		v.vname    =vname;
		v.range    =(Range)(range.clone());
		v.has_undef=has_undef;
		
		if(unit   !=null) v.unit   =unit;
		if(comment!=null) v.comment=comment;
		
		if(tfirst) v.data=new float[t][z][y][x];
		else	   v.data=new float[z][y][x][t];
		
		if(tfirst){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) System.arraycopy(data[l][k][j],0,v.data[l][k][j],0,x);
					
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) System.arraycopy(data[k][j][i],0,v.data[k][j][i],0,t);
		}
		
		return v;
	}
	
	/**
     * used to print out
     */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		sb.append("Varname:\t");	sb.append(vname);	sb.append("\n");
		sb.append("t-count:\t");	sb.append(t);		sb.append("\n");
		sb.append("z-count:\t");	sb.append(z);		sb.append("\n");
		sb.append("y-count:\t");	sb.append(y);		sb.append("\n");
		sb.append("x-count:\t");	sb.append(x);		sb.append("\n");
		sb.append("t-First:\t");	sb.append(tfirst);	sb.append("\n");
		sb.append("unit   :\t");	sb.append(unit);	sb.append("\n");
		sb.append("undef  :\t");	sb.append(undef);	sb.append("\n");
		sb.append("comment:\t");	sb.append(comment);	sb.append("\n");
		
		return sb.toString();
	}
	
	
	/** test
	public static void main(String[] args){
		Matcher m=unitPtn.matcher("[(sv)'M'] (kg m K^-1 s^-2)");
		
		if(m.find()){
			String unitInBracket=m.group();
			System.out.println("1 "+unitInBracket);
			
		};
		
		if(m.find()){
			String unitInBracket=m.group();
			System.out.println("2 "+unitInBracket);
		}
	}*/
}
