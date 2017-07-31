/**
 * @(#)StatisticsBasicAnalysisMethods.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import miniufo.mathsphysics.HarmonicFitter;
import miniufo.mathsphysics.PolynomialFitter;
import miniufo.statistics.StatisticsUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * basic methods of statistics
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class StatisticsBasicAnalysisMethods extends StatisticsApplication{
	
	/**
	 * prevent from instantiate
	 */
	private StatisticsBasicAnalysisMethods(){}
	
	
	/**
     * calculate areal average value
     *
     * @param	v	a given Variable
     *
     * @return	areal mean
     */
	public static Variable arealAverage(Variable v){
		float undef=v.getUndef();
		
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		Range nr=new Range(t,z,1,1);	Range r=v.getRange();
		
		Variable nv=new Variable(v.getName(),nr);
		nv.setComment("area-average of "+v.getComment());
		nv.setUnit(v.getUnit());
		nv.setUndef(undef);
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;	float sum=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
					if(vdata[l][k][j][i]!=undef){ sum+=vdata[l][k][j][i]; count++;}
					
				if(count!=0) nvdata[l][k][0][0]=sum/count;
				else nvdata[l][k][0][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int count=0;	float sum=0;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
					if(vdata[k][j][i][l]!=undef){ sum+=vdata[k][j][i][l]; count++;}
					
				if(count!=0) nvdata[k][0][0][l]=sum/count;
				else nvdata[k][0][0][l]=undef;
			}
		}
		
		nr.setTRange(r);	nr.setZRange(r);
		
		nr.getYRange()[0]=r.getYRange()[0];	nr.getYRange()[1]=r.getYRange()[0];
		nr.getXRange()[0]=r.getXRange()[0];	nr.getXRange()[1]=r.getXRange()[0];
		
		return nv;
	}
	
	/**
     * calculate spatial average value
     *
     * @param	v	a given Variable
     *
     * @return	spatial mean
     */
	public static Variable spatialAverage(Variable v){
		float undef=v.getUndef();
		
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		Range nr=new Range(t,1,1,1);
		Range  r=v.getRange();
		
		Variable nv=new Variable(v.getName(),nr);
		nv.setComment("space-average of "+v.getComment());
		nv.setUnit(v.getUnit());
		nv.setUndef(undef);
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++){
				int count=0;	float sum=0;
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
					if(vdata[l][k][j][i]!=undef){ sum+=vdata[l][k][j][i]; count++;}
					
				if(count!=0) nvdata[l][0][0][0]=sum/count;
				else nvdata[l][0][0][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++){
				int count=0;	float sum=0;
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++)
					if(vdata[k][j][i][l]!=undef){ sum+=vdata[k][j][i][l]; count++;}
					
				if(count!=0) nvdata[0][0][0][l]=sum/count;
				else nvdata[0][0][0][l]=undef;
			}
		}
		
		nr.getZRange()[0]=r.getZRange()[0];	nr.getZRange()[1]=r.getZRange()[0];
		nr.getYRange()[0]=r.getYRange()[0];	nr.getYRange()[1]=r.getYRange()[0];
		nr.getXRange()[0]=r.getXRange()[0];	nr.getXRange()[1]=r.getXRange()[0];
		
		nr.setTRange(r);
		
		return nv;
	}
	
	
	/**
     * calculate monthly mean using daily data (t=365 or 366)
     *
     * @param	v	a given Variable
     *
     * @return	re	monthly mean
     */
	public static Variable monthlyMean(Variable v){
		float undef=v.getUndef();
		
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		if(t!=365&&t!=366) throw new IllegalArgumentException("t length should be 365 or 366");
		
		Range nr=new Range(12,z,y,x);
		Range  r=v.getRange();
		
		Variable nv=new Variable(v.getName(),nr);
		nv.setComment("monthly mean of "+v.getComment());
		nv.setUnit(v.getUnit());
		nv.setUndef(undef);
		
		final int[][] lp={
			{  0, 30},{ 31, 58},{ 59, 89},{ 90,119},{120,150},{151,180},
			{181,211},{212,242},{243,272},{273,303},{304,333},{334,364}
		};
		final int[][] pp={
			{  0, 30},{ 31, 59},{ 60, 90},{ 91,120},{121,151},{152,181},
			{182,212},{213,243},{244,273},{274,304},{305,334},{335,365}
		};
		
		int[][] tag=null;
		if(v.getTCount()==366) tag=lp;
		else tag=pp;
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for (int ll=0;ll<12;ll++){
					int count=0;	float sum=0;
					
					for(int l=tag[ll][0];l<=tag[ll][1];l++)
						if(vdata[l][k][j][i]!=undef){ sum+=vdata[l][k][j][i]; count++;}
					
					if(count!=0) nvdata[ll][k][j][i]=sum/count;
					else nvdata[ll][k][j][i]=undef;
				}
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for (int ll=0;ll<12;ll++){
					int count=0;	float sum=0;

					for(int l=tag[ll][0];l<=tag[ll][1];l++)
						if(vdata[k][j][i][l]!=undef){ sum+=vdata[k][j][i][l]; count++;}
					
					if(count!=0) nvdata[k][j][i][ll]=sum/count;
					else nvdata[k][j][i][ll]=undef;
				}
			}
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=nr.getTRange()[0]+11;
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	
	/**
     * calculate one point correlation map
     *
     * @param	v		a given Variable
     * @param	xx		basic point given as xth point
     * @param	yy		basic point given as yth point
     *
     * @return	spacial atlas
     */
	public static Variable onePointCorrelationMap(Variable v,int xx,int yy){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float[] tmp1=new float[t];
		float[] tmp2=new float[t];
		
		Range  r=v.getRange();
		Range nr=new Range(1,z,y,x);
		
		Variable re=new Variable(v.getName(),true,nr);
		re.setComment("single-point correlation map of "+v.getComment());
		re.setUnit(v.getUnit());
		
		float[][][][] vdata= v.getData();
		float[][][][] rdata=re.getData();
		
		float undef=v.getUndef();
		
		/*** find the basic point ***/
		if(xx<0||xx>x) throw new IllegalArgumentException("Illegal x num");	xx--;
		if(yy<0||yy>y) throw new IllegalArgumentException("Illegal y num");	yy--;
		
		if(v.isTFirst()){
			for(int k=0;k<z;k++){
				for(int l=0;l<t;l++) tmp2[l]=vdata[l][k][yy][xx];
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					for(int l=0;l<t;l++) tmp1[l]=vdata[l][k][j][i];
					
					rdata[0][k][j][i]=cCorrelationCoefficient(tmp1,tmp2,undef);
				}
			}
		}else{
			for(int k=0;k<z;k++){
				tmp2=vdata[k][yy][xx];
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					tmp1=vdata[k][j][i];
					
					rdata[0][k][j][i]=cCorrelationCoefficient(tmp1,tmp2,undef);
				}
			}
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0];
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return re;
	}
	
	
	/**
     * calculate variance of every single point
     *
     * @param	v	a given Variable
     *
     * @return	re	variance
     */
	public static Variable cTVariance(Variable v){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float undef=v.getUndef();
		float[] tmp=new float[t];
		
		Range  r=v.getRange();
		Range nr=new Range(1,z,y,x);
		
		Variable re=new Variable(v.getName(),v.isTFirst(),nr);
		re.setComment("T-Variace of "+v.getComment());
		re.setUnit(v.getUnit());
		re.setUndef(undef);
		
		float[][][][] vdata= v.getData();
		float[][][][] rdata=re.getData();

		if(v.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) tmp[l]=vdata[l][k][j][i];
				
				rdata[0][k][j][i]=cVariance(tmp,undef);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) rdata[k][j][i][0]=cVariance(vdata[k][j][i],undef);
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0];
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return re;
	}
	
	public static Variable cTStandardDeviation(Variable v){
		Variable var=cTVariance(v);
		
		var.setCommentAndUnit("T-std. of "+v.getName());
		
		return var.sqrtEq();
	}
	
	
	/**
     * calculate the range in T dimension
     *
     * @param	v	a given Variable
     *
     * @return	re	range in T dimension
     */
	public static Variable cTRange(Variable v){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float undef=v.getUndef();
		float[] tmp=new float[t];
		
		Range  r=v.getRange();
		Range nr=new Range(1,z,y,x);
		
		Variable re=new Variable(v.getName()+"rng",v.isTFirst(),nr);
		re.setCommentAndUnit("T-Range of "+v.getName());
		re.setUndef(undef);
		
		float[][][][] vdata= v.getData();
		float[][][][] rdata=re.getData();

		if(v.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) tmp[l]=vdata[l][k][j][i];
				
				rdata[0][k][j][i]=cRange(tmp,undef);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) rdata[k][j][i][0]=cRange(vdata[k][j][i],undef);
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0];
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return re;
	}
	
	
	/**
     * polynomial fit along dimension T
     *
     * @param	v	a given Variable
     *
     * @return	re	fitting result
     */
	public static Variable polynomialFitAlongT(Variable v,int n){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float[] tmp=new float[t];
		
		Variable re=new Variable(v.getName()+"fit",v);
		re.setComment(n+" order polynomial fit of "+v.getComment());
		re.setUnit(v.getUnit());
		
		PolynomialFitter pf=new PolynomialFitter(n,t);
		
		float[][][][] vdata= v.getData();
		float[][][][] rdata=re.getData();
		
		if(v.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) tmp[l]=vdata[l][k][j][i];
				
				pf.fit(tmp);	float[] yy=pf.cValues();
				
				for(int l=0;l<t;l++) rdata[l][k][j][i]=yy[l];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				pf.fit(vdata[k][j][i]);
				System.arraycopy(pf.cValues(),0,rdata[k][j][i],0,t);
			}
		}
		
		return re;
	}
	
	
	/**
     * calculate the spatial similarity coefficient
     *
     * @param	v		a given Variable
     * @param	delay	delay time
     */
	public static Variable cSpatialSimilarity(Variable v,int delay){
		if(delay<0) throw new IllegalArgumentException("delay shoud not be negative");
		
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float undef=v.getUndef();
		
		Range nr=new Range(t-delay,z,1,1);
		Range  r=v.getRange();
		
		Variable nv=new Variable("spsm",nr);
		nv.setComment("spatial similarity of "+v.getComment());
		nv.setUnit(v.getUnit());
		nv.setUndef(undef);
		
		float[] buff1=new float[x*y];
		float[] buff2=new float[x*y];
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int k=0;k<z;k++)
			for(int l=0;l<t-delay;l++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){ buff1[j*x+i]=vdata[l][k][j][i]; buff2[j*x+i]=vdata[l+delay][k][j][i];}
				
				nvdata[l][k][0][0]=(float)StatisticsUtil.cCorrelationCoefficient(buff1,buff2);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int l=0;l<t-delay;l++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){ buff1[j*x+i]=vdata[k][j][i][l]; buff2[j*x+i]=vdata[k][j][i][l+delay];}
				
				nvdata[k][0][0][l]=(float)StatisticsApplication.cCorrelationCoefficient(buff1,buff2,undef);
			}
		}
		
		nr.getTRange()[0]=r.getTRange()[0];	nr.getTRange()[1]=r.getTRange()[0]+r.getTRange()[2]-1;
		nr.setYRange(r.getYRange()[0]);
		nr.setXRange(r.getXRange()[0]);
		nr.setZRange(r);
		
		return nv;
	}
	
	/**
     * calculate the spatial similarity coefficient
     *
     * @param	v1		a given Variable
     * @param	v2		a given Variable
     */
	public static Variable cSpatialSimilarity(Variable v1,Variable v2){
		if(!v1.isLike(v2)) throw new IllegalArgumentException("dimensions not same");
		
		int t=v1.getTCount(),	z=v1.getZCount(),	y=v1.getYCount(),	x=v1.getXCount();
		
		float undef=v1.getUndef();
		
		Range nr=new Range(t,z,1,1);
		Range  r=v1.getRange();
		
		Variable nv=new Variable("spsm",v1.isTFirst(),nr);
		nv.setComment("spatial similarity of "+v1.getComment()+" and "+v2.getComment());
		nv.setUnit("1");
		nv.setUndef(undef);
		
		float[] buff1=new float[x*y];
		float[] buff2=new float[x*y];
		
		float[][][][] v1data=v1.getData();
		float[][][][] v2data=v2.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v1.isTFirst()){
			for(int k=0;k<z;k++)
			for(int l=0;l<t;l++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){ buff1[j*x+i]=v1data[l][k][j][i]; buff2[j*x+i]=v2data[l][k][j][i];}
				
				nvdata[l][k][0][0]=(float)StatisticsApplication.cCorrelationCoefficient(buff1,buff2,undef);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int l=0;l<t;l++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){ buff1[j*x+i]=v1data[k][j][i][l]; buff2[j*x+i]=v2data[k][j][i][l];}
				
				nvdata[k][0][0][l]=(float)StatisticsApplication.cCorrelationCoefficient(buff1,buff2,undef);
			}
		}
		
		nr.getYRange()[0]=r.getYRange()[0];	nr.getYRange()[1]=r.getYRange()[0];
		nr.getXRange()[0]=r.getXRange()[0];	nr.getXRange()[1]=r.getXRange()[0];
		
		nr.setYRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	/**
     * calculate the vertical spatial similarity coefficient
     *
     * @param	v1		a given Variable
     * @param	v2		a given Variable
     */
	public static Variable cVerticalSpatialSimilarity(Variable v1,Variable v2){
		if(!v1.isLike(v2)) throw new IllegalArgumentException("dimensions not same");
		
		int t=v1.getTCount(),	z=v1.getZCount(),	y=v1.getYCount(),	x=v1.getXCount();
		
		float undef=v1.getUndef();
		
		Range nr=new Range(t,1,1,x);
		Range r=v1.getRange();
		
		Variable nv=new Variable("vspsm",v1.isTFirst(),nr);	nv.setUndef(undef);
		nv.setComment("vertical spatial similarity of "+v1.getComment()+" and "+v2.getComment());
		nv.setUnit("1");
		
		float[] buff1=new float[z*y];
		float[] buff2=new float[z*y];
		
		float[][][][] v1data=v1.getData();
		float[][][][] v2data=v2.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v1.isTFirst()){
			for(int l=0;l<t;l++)
				for(int i=0;i<x;i++){
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++){ buff1[k*y+j]=v1data[l][k][j][i]; buff2[k*y+j]=v2data[l][k][j][i];}
				
				nvdata[l][0][0][i]=(float)StatisticsApplication.cCorrelationCoefficient(buff1,buff2,undef);
			}
			
		}else{
			for(int l=0;l<t;l++)
				for(int i=0;i<x;i++){
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++){ buff1[k*y+j]=v1data[k][j][i][l]; buff2[k*y+j]=v2data[k][j][i][l];}
				
				nvdata[0][0][i][l]=(float)StatisticsApplication.cCorrelationCoefficient(buff1,buff2,undef);
			}
		}
		
		nr.getZRange()[0]=r.getZRange()[0];	nr.getZRange()[1]=r.getZRange()[0];
		nr.getYRange()[0]=r.getYRange()[0];	nr.getYRange()[1]=r.getYRange()[0];
		
		nr.setTRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	
	/**
     * calculate the correlation of temporal series and areal data
     *
     * @param	tv		a given Variable of time
     * @param	av		a given Variable of area
     *
     * @return	correlation of an area
     */
	public static Variable taco(Variable tv,Variable av){
		if(tv.getXCount()!=1||tv.getYCount()!=1||tv.getZCount()!=1)
			throw new IllegalArgumentException("temporal serial only");
		
		if(tv.getTCount()!=av.getTCount()||tv.isTFirst()!=av.isTFirst())
			throw new IllegalArgumentException("dimensions not same");
		
		float undef=tv.getUndef();
		
		int t=tv.getTCount(),	z=av.getZCount(),	y=av.getYCount(),	x=av.getXCount();
		
		Range nr=new Range(1,z,y,x);
		Range  r=av.getRange();
		
		Variable nv=new Variable("t"+av.getName(),tv.isTFirst(),nr);
		nv.setUndef(undef);	nv.setCommentAndUnit("spacial correlation map");
		
		float[][][][] tvdata=tv.getData();
		float[][][][] avdata=av.getData();
		float[][][][] nvdata=nv.getData();
		
		if(tv.isTFirst()){
			float[] tmp1=new float[t];
			float[] tmp2=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++){
					tmp1[l]=tvdata[l][0][0][0];
					tmp2[l]=avdata[l][k][j][i];
				}
				
				nvdata[0][k][j][i]=cCorrelationCoefficient(tmp1,tmp2,undef);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				nvdata[k][j][i][0]=cCorrelationCoefficient(tvdata[0][0][0],avdata[k][j][i],undef);
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0];
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	/**
     * calculate the lead or lag correlation of temporal series and areal data
     *
     * @param	tv		a given Variable of time
     * @param	av		a given Variable of area
     *
     * @return	correlation of an area
     */
	public static Variable cLeadLagTaco(Variable tv,Variable av,int delay){
		if(tv.getXCount()!=1||tv.getYCount()!=1||tv.getZCount()!=1)
			throw new IllegalArgumentException("temporal serial only");
		
		if(tv.getTCount()!=av.getTCount()||tv.isTFirst()!=av.isTFirst())
			throw new IllegalArgumentException("dimensions not same");
		
		float undef=tv.getUndef();
		
		int t=tv.getTCount(),	z=av.getZCount(),	y=av.getYCount(),	x=av.getXCount();
		
		Range nr=new Range(1,z,y,x),r=av.getRange();
		
		Variable nv=new Variable("t"+av.getName(),tv.isTFirst(),nr);
		nv.setUndef(undef);	nv.setCommentAndUnit("lead-lag spacial correlation map");
		
		float[][][][] tvdata=tv.getData();
		float[][][][] avdata=av.getData();
		float[][][][] nvdata=nv.getData();
		
		if(tv.isTFirst()){
			float[] tmp1=new float[t];
			float[] tmp2=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++){
					tmp1[l]=tvdata[l][0][0][0];
					tmp2[l]=avdata[l][k][j][i];
				}
				
				nvdata[0][k][j][i]=StatisticsUtil.cLeadLagCorrelationCoefficient(tmp1,tmp2,delay);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			nvdata[k][j][i][0]=StatisticsUtil.cLeadLagCorrelationCoefficient(tvdata[0][0][0],avdata[k][j][i],delay);
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0];
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	/**
     * calculate the self lead or lag correlation of a given variable
     *
     * @param	v		a given Variable
     * @param	delay	delays
     *
     * @return	correlation of an area
     */
	public static Variable cSelfLeadLagCorrelation(Variable v,int delay){
		float undef=v.getUndef();
		
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		Range nr=new Range(1,z,y,x),r=v.getRange();
		
		Variable nv=new Variable(v.getName()+"lag"+delay,v.isTFirst(),nr);
		nv.setUndef(undef);	nv.setCommentAndUnit("lead-lag spacial correlation map");
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			float[] tmp=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
		hh:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++)
				if(vdata[l][k][j][i]==undef){ nvdata[0][k][j][i]=undef; continue hh;}
				else tmp[l]=vdata[l][k][j][i];
				
				nvdata[0][k][j][i]=StatisticsUtil.cLeadLagCorrelationCoefficient(tmp,tmp,delay);
				
				if(Float.isNaN(nvdata[0][k][j][i])) nvdata[0][k][j][i]=undef;
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
		hh:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++)
				if(vdata[k][j][i][l]==undef){ nvdata[k][j][i][0]=undef; continue hh;}
				
				nvdata[k][j][i][0]=StatisticsUtil.cLeadLagCorrelationCoefficient(
					vdata[k][j][i],vdata[k][j][i],delay
				);
				
				if(Float.isNaN(nvdata[k][j][i][0])) nvdata[k][j][i][0]=undef;
			}
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0];
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	
	/**
     * calculate the skewness coefficient of a given Variable
     *
     * @param	pv	a given Variable
     *
     * @return	skewness coefficient
     */
	public static Variable cSkewnessCoefficient(Variable v){
		float undef=v.getUndef();
		
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		Range nr=new Range(1,z,y,x);
		Range  r=v.getRange();
		
		Variable nv=new Variable(v.getName(),nr);
		nv.setUndef(undef);	nv.setCommentAndUnit("skewness coefficient");
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] buf=new float[t];
				
				for(int l=0;l<t;l++) buf[l]=vdata[l][k][j][i];
				
				nvdata[0][k][j][i]=StatisticsUtil.cSkewnessCoefficient(buf);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				nvdata[k][j][i][0]=StatisticsUtil.cSkewnessCoefficient(vdata[k][j][i]);
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0];
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	
	/**
     * test the significance of difference of mean value (composite)
     *
     * @param	v1	composite variable 1
     * @param	v2	composite variable 2
     *
     * @return	T variable
     */
	public static Variable testDiffSignificance(Variable v1,Variable v2){
		if(!v1.isAreaLike(v2)) throw new IllegalArgumentException("not the same dimension");
		
		float undef=v1.getUndef();
		
		int t1=v1.getTCount();	int z=v1.getZCount();
		int t2=v2.getTCount();	int y=v1.getYCount();	int x=v1.getXCount();
		
		Range nr=new Range(1,z,y,x);
		Range r1=v1.getRange();
		
		Variable nv=new Variable("test"+v1.getName(),nr);
		nv.setUndef(undef);	nv.setCommentAndUnit("T-test Variable for "+v1.getName());
		
		float[][][][] v1data=v1.getData();
		float[][][][] v2data=v2.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v1.isTFirst()){
			float[] buf1=new float[t1];
			float[] buf2=new float[t2];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t1;l++) buf1[l]=v1data[l][k][j][i];
				for(int l=0;l<t2;l++) buf2[l]=v2data[l][k][j][i];
				
				nvdata[0][k][j][i]=StatisticsUtil.TTest(buf1,buf2);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			nvdata[k][j][i][0]=StatisticsUtil.TTest(v1data[k][j][i],v2data[k][j][i]);
		}
		
		nr.setTRange(r1.getTRange()[0]);
		nr.setZRange(r1);
		nr.setYRange(r1);
		nr.setXRange(r1);
		
		return nv;
	}
	
	
	/**
     * calculate the amplitudes given periods of harmonics along T-dimension
     *
     * @param	v	a given variable
     * 
     * @return	re	amplitudes
     */
	public static Variable[] cHarmonicAmplitudes(Variable v,float... Ts){
		int t=v.getTCount(),z=v.getZCount(),y=v.getYCount(),x=v.getXCount();
		
		float undef=v.getUndef();
		
		Variable[] re=new Variable[Ts.length];
		
		for(int m=0,M=Ts.length;m<M;m++){
			re[m]=new Variable(v.getName()+"amp"+m,v.isTFirst(),new Range(1,z,y,x));
			re[m].setCommentAndUnit("amplitude of harmonics of period "+Ts[m]);
			re[m].setUndef(undef);
		}
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			HarmonicFitter hf=new HarmonicFitter(t,Ts);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cc:		for(int i=0;i<x;i++){
				float[] buf=new float[t];
				
				for(int l=0;l<t;l++) buf[l]=vdata[l][k][j][i];
				
				if(buf[0]==undef){
					for(int m=0,M=Ts.length;m<M;m++) re[m].getData()[0][k][j][i]=undef;
					continue cc;
				}
				
				hf.fit(buf);
				for(int m=0,M=Ts.length;m<M;m++) re[m].getData()[0][k][j][i]=hf.getAmplitudes()[m];
			}
			
		}else{
			HarmonicFitter hf=new HarmonicFitter(t,Ts);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cc:		for(int i=0;i<x;i++){
				float[] buf=vdata[k][j][i];
				
				if(buf[0]==undef){
					for(int m=0,M=Ts.length;m<M;m++) re[m].getData()[k][j][i][0]=undef;
					continue cc;
				}
				
				hf.fit(buf);
				for(int m=0,M=Ts.length;m<M;m++) re[m].getData()[k][j][i][0]=hf.getAmplitudes()[m];
			}
		}
		
		for(int m=0,M=Ts.length;m<M;m++){
			re[m].getRange().setTRange(v.getRange().getTRange()[0]);
			re[m].getRange().setZRange(v.getRange());
			re[m].getRange().setYRange(v.getRange());
			re[m].getRange().setXRange(v.getRange());
		}
		
		return re;
	}
	
	/**
     * calculate the variance ellipses
     *
     * @param	ua	zonal component of eddy velocity
     * @param	va	meridional component of eddy velocity
     * 
     * @return	re	[0] is major axis, [1] is minor axis and [2] is theta
     */
	public static Variable[] cVarianceEllipse(Variable ua,Variable va){
		int t=ua.getTCount(),z=ua.getZCount(),y=ua.getYCount(),x=ua.getXCount();
		
		float undef=ua.getUndef();
		
		Variable[] re=new Variable[3];
		
		re[0]=new Variable("major",ua.isTFirst(),new Range(t,z,y,x));
		re[1]=new Variable("minor",ua.isTFirst(),new Range(t,z,y,x));
		re[2]=new Variable("theta",ua.isTFirst(),new Range(t,z,y,x));
		
		re[0].setCommentAndUnit("variance along the major axis");
		re[1].setCommentAndUnit("variance along the minor axis");
		re[2].setCommentAndUnit("the direction of the axis of principal variability");
		
		re[0].setUndef(undef);
		re[1].setUndef(undef);
		re[2].setUndef(undef);
		
		float[][][][] uadata=ua.getData();
		float[][][][] vadata=va.getData();
		float[][][][] r0data=re[0].getData();
		float[][][][] r1data=re[1].getData();
		float[][][][] r2data=re[2].getData();
		
		if(ua.isTFirst()){
			float[] ubuf=new float[t];
			float[] vbuf=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	brk:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++){
					ubuf[l]=uadata[l][k][j][i];
					vbuf[l]=vadata[l][k][j][i];
					
					if(ubuf[l]==undef||vbuf[l]==undef){
						r0data[0][k][j][i]=undef;
						r1data[0][k][j][i]=undef;
						r2data[0][k][j][i]=undef;
						continue brk;
					}
				}
				
				float uvar=squaredMean(ubuf);
				float vvar=squaredMean(vbuf);
				float cvar=productMean(ubuf,vbuf);
				
				r0data[0][k][j][i]=(uvar+vvar+(float)Math.sqrt((uvar-vvar)*(uvar-vvar)+4*cvar*cvar))/2;
				r1data[0][k][j][i]=(uvar+vvar)-r0data[0][k][j][i];
				r2data[0][k][j][i]=(float)Math.atan2(r0data[0][k][j][i]-uvar,cvar);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	brk:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++){
					if(uadata[k][j][i][l]==undef||vadata[k][j][i][l]==undef){
						r0data[k][j][i][0]=undef;
						r1data[k][j][i][0]=undef;
						r2data[k][j][i][0]=undef;
						continue brk;
					}
				}
				
				float uvar=squaredMean(uadata[k][j][i]);
				float vvar=squaredMean(vadata[k][j][i]);
				float cvar=productMean(uadata[k][j][i],vadata[k][j][i]);
				
				r0data[k][j][i][0]=(uvar+vvar+(float)Math.sqrt((uvar-vvar)*(uvar-vvar)+4*cvar*cvar))/2;
				r1data[k][j][i][0]=(uvar+vvar)-r0data[k][j][i][0];
				r2data[k][j][i][0]=(float)Math.atan2(r0data[k][j][i][0]-uvar,cvar);
			}
		}
		
		re[0].setRange(ua.getRange());	re[0].getRange().setTRange(ua.getRange().getTRange()[0]);
		re[1].setRange(ua.getRange());	re[1].getRange().setTRange(ua.getRange().getTRange()[0]);
		re[2].setRange(ua.getRange());	re[2].getRange().setTRange(ua.getRange().getTRange()[0]);
		
		return re;
	}
	
	
	/**
     * Calculate variance along X-dimension.
     *
     * @param	v	a given variable
     * 
     * @return	var	variance
     */
	public static Variable cVarianceAlongX(Variable v){
		float undef=v.getUndef();
		
		int t=v.getTCount(),z=v.getZCount(),y=v.getYCount(),x=v.getXCount();
		
		Range nr=new Range(t,z,y,1);	Range r=v.getRange();
		
		Variable nv=new Variable(v.getName(),v.isTFirst(),nr);
		nv.setComment("variance along X-dimension of "+v.getComment());
		nv.setUnit("("+v.getUnit()+")^2");
		nv.setUndef(undef);
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			nvdata[l][k][j][0]=cVariance(vdata[l][k][j],undef);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				float[] buf=new float[x];
				
				for(int i=0;i<x;i++) buf[i]=vdata[k][j][i][l];
				
				nvdata[k][j][0][l]=cVariance(buf,undef);
			}
		}
		
		nr.setTRange(r);
		nr.setZRange(r);
		nr.setYRange(r);
		nr.getXRange()[0]=r.getXRange()[0];
		nr.getXRange()[1]=r.getXRange()[0];
		
		return nv;
	}
	
	
	/*** helper methods ***/
	private static float squaredMean(float[] data){
		if(data.length==1||data.length==0) return 0;
		
		float average=StatisticsUtil.cArithmeticMean(data),variance=0;
		
		for(float ei:data){
			float tmp=ei-average;
			variance+=tmp*tmp;
		}
		
		return variance/data.length;
	}
	
	private static float productMean(float[] data1,float[] data2){
		if(data1.length==1||data1.length==0) return 0;
		
		float mean_a=StatisticsUtil.cArithmeticMean(data1);
		float mean_b=StatisticsUtil.cArithmeticMean(data2);
		float covariance=0;
		
		for(int i=0,I=data1.length;i<I;i++) covariance+=(data1[i]-mean_a)*(data2[i]-mean_b);
		
		return covariance/data1.length;
	}
	
	
	/** test
	public static void main(String arg[]){
		try{
			miniufo.descriptor.DataDescriptor dd=miniufo.descriptor.DescriptorFactory.getDescriptor("d:/anom.ctl");
			
			System.out.println(dd.getDZDef());
			
			miniufo.diagnosis.Range r=new miniufo.diagnosis.Range("",dd);
			
			Variable anom=new Variable("anom",r);
			
			miniufo.io.DataRead dr=miniufo.io.DataIOFactory.getDataRead(dd);
			dr.readData(anom);	dr.closeFile();
			
			Variable tend1=polynomialFit(anom,1);	tend1.setName("t1");
			Variable tend2=polynomialFit(anom,2);	tend2.setName("t2");
			Variable tend3=polynomialFit(anom,3);	tend3.setName("t3");
			Variable tend4=polynomialFit(anom,4);	tend4.setName("t4");
			
			miniufo.io.DataWrite dw=miniufo.io.DataIOFactory.getDataWrite(dd,"d:/re.dat");
			dw.writeData(dd,tend1,tend2,tend3,tend4,anom);
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
