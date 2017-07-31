/**
 * @(#)DataBaseUtil.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.database;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.LagrangianUtil;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;


/**
 * basic methods for database package
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DataBaseUtil{
	//
	public static final float undef=-9999.0f;
	
	
	/**
     * binning the types of records
     *
     * @param	dd		a DataDescriptor template used to describe the grids
     * @param	ls		records in a list
     * @param	dlon	delta lon for periodicX boundary condition
     * @param	wthre	wind threshold
     * @param	types	types that need to be binned
     * 
     * @return	re		type counts corresponding to types
     */
	public static Variable[] binningTCTypeCount(DataDescriptor dd,List<Typhoon> ls,float dlon,float wthre,TYPE... types){
		int tcount=dd.getTCount();
		
		if(types==null||types.length==0) throw new IllegalArgumentException("specify some types");
		
		Variable[] re=new Variable[types.length];
		
		for(int i=0,I=types.length;i<I;i++){
			re[i]=new Variable(types[i].name(),new Range(tcount,1,dd.getYCount(),dd.getXCount()));
			re[i].setCommentAndUnit("binned type count");
			re[i].setUndef(undef);
			
			TYPE t=types[i];
			
			Function<Record,Float> RcToF=r->1f;
			Predicate<Record> cond=r->{
				boolean isType=Typhoon.getType(r.getDataValue(4))==t;
				boolean largeW=r.getDataValue(2)>=wthre;
				
				return isType&&largeW;
			};
			
			binningDataMeanByCondition(dd,ls,RcToF,cond,re[i],dlon);
		}
		
		return re;
	}
	
	public static Variable[] binningTCTypeCount(DataDescriptor dd,List<Typhoon> ls,TYPE... types){
		return binningTCTypeCount(dd,ls,0,0,types);
	}
	
	// add by Lynn
	public static Variable binningTCACEByRange(DataDescriptor dd,List<Typhoon> ls,float dlon,float wndmin,float wndmax){
		int tcount=dd.getTCount();
		
		Variable re=new Variable("ace",new Range(tcount,1,dd.getYCount(),dd.getXCount()));
		re.setComment("accumulated cyclone energy");
		re.setUndef(undef);
		
		Function<Record,Float> RcToF=r->r.getDataValue(2)*r.getDataValue(2);
		Predicate<Record> cond=r->r.getDataValue(2)>=wndmin&&r.getDataValue(2)<wndmax;
		
		binningDataSumByCondition(dd,ls,RcToF,cond,re,dlon);
		
		return re;
	}

	// add by Lynn
	public static Variable binningTCCountByRange(DataDescriptor dd,List<Typhoon> ls,float dlon,float wndmin,float wndmax){
		int tcount=dd.getTCount();
		
		Variable re=new Variable("freq",new Range(tcount,1,dd.getYCount(),dd.getXCount()));
		re.setComment("accumulated cyclone energy");
		re.setUndef(undef);
		
		Function<Record,Float> RcToF=r->1f;
		Predicate<Record> cond=r->r.getDataValue(2)>=wndmin&&r.getDataValue(2)<wndmax;
		
		binningDataSumByCondition(dd,ls,RcToF,cond,re,dlon);
		
		return re;
	}
	
	public static Variable binningTCACE(DataDescriptor dd,List<Typhoon> ls,float dlon,float wthre){
		int tcount=dd.getTCount();
		
		Variable re=new Variable("ace",new Range(tcount,1,dd.getYCount(),dd.getXCount()));
		re.setCommentAndUnit("accumulated cyclone energy");
		re.setUndef(undef);
		
		Function<Record,Float> RcToF=r->r.getDataValue(2)*r.getDataValue(2);
		Predicate<Record> cond=r->r.getDataValue(2)>=wthre;
		
		binningDataMeanByCondition(dd,ls,RcToF,cond,re,dlon);
		
		return re;
	}
	
	public static Variable binningTCACE(DataDescriptor dd,List<Typhoon> ls){
		return binningTCACE(dd,ls,0,0);
	}
	
	public static Variable binningTCGenesisFrequency(DataDescriptor dd,List<Typhoon> ls,float dlon,float wthre){
		int tcount=dd.getTCount();
		
		Variable re=new Variable("gen",new Range(tcount,1,dd.getYCount(),dd.getXCount()));
		re.setCommentAndUnit("TC genesis frequency");
		re.setUndef(undef);
		
		Function<Record,Float> RcToF=r->1f;
		Predicate<Record> cond=r->r.getDataValue(2)>=wthre;
		
		binningDataMeanByCondition(dd,ls,RcToF,cond,re,dlon);
		
		return re;
	}
	
	public static Variable binningTCGenesisFrequency(DataDescriptor dd,List<Typhoon> ls){
		return binningTCGenesisFrequency(dd,ls,0,0);
	}
	
	
	/**
     * binning particle data
     *
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	dlon		delta lon for periodicX boundary condition
     * @param	maskUndef	mask the field with undef
     * @param	idx			indices of AttachedData
     * 
     * @return	re			length is equal to that of AttachedData
     */
	public static Variable[] binningData
	(DataDescriptor dd,List<? extends Particle> ls,float dlon,boolean maskUndef,int... idx){
		int tcount=dd.getTCount();
		
		if(idx==null||idx.length==0) throw new IllegalArgumentException("specify some indices of AttachedData");
		
		Variable[] re=new Variable[idx.length];
		
		String[] ads=ls.get(0).getDataNames();
		
		for(int i:idx){
			re[i]=new Variable("data"+i,new Range(tcount,1,dd.getYCount(),dd.getXCount()));
			re[i].setCommentAndUnit(ads[i]);
			re[i].setUndef(undef);
			
			Function<Record,Float> RcToF=r->r.getDataValue(i);
			Predicate<Record> cond=r->r.getDataValue(i)!=undef;
			
			binningDataMeanByCondition(dd,ls,RcToF,cond,re[i],dlon);
		}
		
		return re;
	}
	
	public static Variable[] binningData(DataDescriptor dd,List<? extends Particle> ls,int... idx){
		return binningData(dd,ls,0,true,idx);
	}
	
	/**
     * binning data by date
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	dlon		delta lon for periodicX boundary condition
     * @param	maskUndef	mask the field with undef
     * @param	idx			indices of AttachedData
     * @param	part		part of date in String form: ["year","month","date"]
     * @param	values		values of the date that would be binned
     * 
     * @return	re			data in those specific date
     */
	public static Variable[] binningDataByDate
	(DataDescriptor dd,List<? extends Particle> ls,float dlon,boolean maskUndef,int[] idx,String part,int... values){
		int tcount=dd.getTCount();
		
		if(tcount!=1) throw new IllegalArgumentException("T-count should be 1");
		if(idx==null||idx.length==0) throw new IllegalArgumentException("specify some indices of AttachedData");
		if(values==null||values.length==0) throw new IllegalArgumentException("specify some values of date");
		
		Variable[] re=new Variable[idx.length];
		
		String[] ads=ls.get(0).getDataNames();
		
		for(int vi:idx){
			re[vi]=new Variable("data"+vi,true,new Range(tcount,1,dd.getYCount(),dd.getXCount()));
			re[vi].setCommentAndUnit(ads[vi]);
			re[vi].setUndef(undef);
			
			Function<Record,Float> RcToF=r->r.getDataValue(vi);
			Predicate<Record> cond=r->{
				if(r.getDataValue(vi)==undef) return false;
				
				int value=getPartValue(r.getTime(),part);
				
				for(int v:values) if(v==value) return true;
				
				return false;
			};
			
			binningDataMeanByCondition(dd,ls,RcToF,cond,re[vi],dlon);
		}
		
		return re;
	}
	
	public static Variable[] binningDataByDate
	(DataDescriptor dd,List<? extends Particle> ls,boolean maskUndef,int[] idx,String part,int... values){
		return binningDataByDate(dd,ls,0,maskUndef,idx,part,values);
	}
	
	/**
     * binning data of different seasons
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	dlon		delta lon for periodicX boundary condition
     * @param	maskUndef	mask the field with undef
     * @param	season		e.g., int[][] season=new int[][]{
								{3, 4, 5},	{ 6,7,8},	spring, summer,
								{9,10,11},	{12,1,2}	autumn, winter
							};
     * @param	idx			indices of AttachedData
     * 
     * @return	re			means of different data of different seasons, [seasons][AttachedData]
     */
	public static Variable[][] binningSeasonalData
	(DataDescriptor dd,List<? extends Particle> ls,float dlon,boolean maskUndef,int[][] seasons,int... idx){
		int season=seasons.length;
		int sealen=seasons[0].length;
		
		for(int[] s:seasons)
		if(s.length!=sealen) throw new IllegalArgumentException("season lengths not equal");
		
		if(idx.length==0) throw new IllegalArgumentException("no idx specified");
		
		Variable[][] re=new Variable[season][idx.length];
		
		for(int i=0;i<season;i++){
			Variable[] tmp=binningDataByDate(dd,ls,dlon,maskUndef,idx,"month",seasons[i]);
			
			for(int j=0,J=tmp.length;j<J;j++){
				re[i][j]=tmp[j];
				re[i][j].setName(tmp[j].getName()+"_"+i);
			}
		}
		
		return re;
	}
	
	public static Variable[][] binningSeasonalData
	(DataDescriptor dd,List<? extends Particle> ls,boolean maskUndef,int[][] seasons,int... idx){
		return binningSeasonalData(dd,ls,0,maskUndef,seasons,idx);
	}
	
	
	/**
     * binning the data count
     *
     * @param	dd		a DataDescriptor template used to describe the grids
     * @param	ls		a list of data
     * @param	dlon	delta lon for periodicX boundary condition
     * 
     * @return	re		count of data within each bin
     */
	public static Variable binningCount(DataDescriptor dd,List<? extends Particle> ls,float dlon){
		Variable re=new Variable("count",new Range(dd.getTCount(),1,dd.getYCount(),dd.getXCount()));
		
		re.setCommentAndUnit("binned count");
		re.setUndef(undef);
		
		binningDataSumByCondition(dd,ls,r->1f,r->true,re,dlon);
		
		return re;
	}
	
	public static Variable binningCount(DataDescriptor dd,List<? extends Particle> ls){
		return binningCount(dd,ls,0);
	}
	
	public static Variable binningMedianCount(DataDescriptor dd,List<? extends Particle> ls,float dlon){
		int tcount=dd.getTCount();
		
		Variable re=new Variable("count",true,new Range(tcount,1,dd.getYCount(),dd.getXCount()));
		re.setCommentAndUnit("binned median-position count");
		re.setUndef(undef);
		
		Consumer<Particle> csmr=(tcount==1)?
			p->binningRecordTInvariantSum(dd,p.getRecord(p.getMedianIndex()),r->1f,r->true,null,re,dlon):
			p->binningRecordTVariantSum(dd,p.getRecord(p.getMedianIndex()),r->1f,r->true,null,re,dlon);
		
		ls.stream().forEach(csmr);
		
		return re;
	}
	
	public static Variable binningMedianCount(DataDescriptor dd,List<? extends Particle> ls){
		return binningMedianCount(dd,ls,0);
	}
	
	/**
     * binning count by date
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	dlon		delta lon for periodicX boundary condition
     * @param	maskUndef	mask the field with undef
     * @param	part		part of date in String form: ["year","month","date"]
     * @param	values		values of the date that would be binned
     * 
     * @return	re			count in those specific date
     */
	public static Variable binningCountByDate
	(DataDescriptor dd,List<? extends Particle> ls,float dlon,String part,int... values){
		int tcount=dd.getTCount();
		
		if(tcount!=1) throw new IllegalArgumentException("T-count should be 1");
		if(values==null||values.length==0) throw new IllegalArgumentException("specify some values of date");
		
		Predicate<Record> cond=r->{
			int value=getPartValue(r.getTime(),part);
			
			for(int v:values) if(v==value) return true;
			
			return false;
		};
		
		Variable re=new Variable("count",new Range(dd.getTCount(),1,dd.getYCount(),dd.getXCount()));
		
		re.setCommentAndUnit("binned count for "+part+" = "+Arrays.toString(values));
		re.setUndef(undef);
		
		binningDataSumByCondition(dd,ls,r->1f,cond,re,dlon);
		
		return re;
	}
	
	public static Variable binningCountByDate
	(DataDescriptor dd,List<? extends Particle> ls,String part,int... values){
		return binningCountByDate(dd,ls,0,part,values);
	}
	
	
	/**
     * binning count of different seasons
     * Spring [3, 4, 5]; Summer [6, 7, 8]; Autumn [9, 10, 11]; Winter [12, 1, 2]
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	dlon		delta lon for periodicX boundary condition
     * @param	season		e.g., int[][] season=new int[][]{
	 *							{3, 4, 5},	{ 6,7,8},	spring, summer,
	 *							{9,10,11},	{12,1,2}	autumn, winter
	 *						};
     * 
     * @return	re			count of different seasons
     */
	public static Variable[] binningSeasonalCount
	(DataDescriptor dd,List<? extends Particle> ls,float dlon,int[][] seasons){
		int sc=seasons.length;
		
		Variable[] cs=new Variable[sc];
		
		for(int i=0;i<sc;i++){
			cs[i]=binningCountByDate(dd,ls,dlon,"month",seasons[i]);
			cs[i].setName("count_mon"+(i+1));
		}
		
		return cs;
	}
	
	public static Variable[] binningSeasonalCount
	(DataDescriptor dd,List<? extends Particle> ls,int[][] seasons){
		return binningSeasonalCount(dd,ls,0,seasons);
	}
	
	
	/*** helper methods ***/
	
	/**
     * common interface for binning Record data
     * 
     * @param	dd		a DataDescriptor template used to describe the grids
     * @param	ls		a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     * 
     * @return	re		binning result return in a Variable
     */
	private static void binningDataMeanByCondition
	(DataDescriptor dd,List<? extends Particle> ls,Function<Record,Float> RcToF,Predicate<Record> cond,Variable var,float dlon){
		int tcount=dd.getTCount();
		
		if(tcount==1){
			int[][] count=new int[dd.getYCount()][dd.getXCount()];
			
			Consumer<Record> csmr=r->binningRecordTInvariantSum(dd,r,RcToF,cond,count,var,dlon);
			
			LagrangianUtil.asRecordStream(ls).forEach(csmr);
			
			meanTInvariant(count,true,var);
			
		}else{
			int[][][] count=new int[dd.getYCount()][dd.getXCount()][dd.getTCount()];
			
			Consumer<Record> csmr=r->binningRecordTVariantSum(dd,r,RcToF,cond,count,var,dlon);
			
			LagrangianUtil.asRecordStream(ls).forEach(csmr);
			
			meanTVariant(count,true,var);
		}
	}
	
	private static void binningDataSumByCondition
	(DataDescriptor dd,List<? extends Particle> ls,Function<Record,Float> RcToF,Predicate<Record> cond,Variable var,float dlon){
		int tcount=dd.getTCount();
		
		Consumer<Record> csmr=(tcount==1)?
			r->binningRecordTInvariantSum(dd,r,RcToF,cond,null,var,dlon):
			r->binningRecordTVariantSum(dd,r,RcToF,cond,null,var,dlon);
		
		LagrangianUtil.asRecordStream(ls).forEach(csmr);
	}
	
	/**
	 * binning a single Record accumulatively
	 * 
     * @param	dd		a DataDescriptor template used to describe the grids
     * @param	record	a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     * @param	count	count array for each grid
     * @param	var		variable need to store the result
     * @param	dlon	delta lon for periodicX boundary condition
	 */
	private static void binningRecordTInvariantSum
	(DataDescriptor dd,Record record,Function<Record,Float> RcToF,Predicate<Record> cond,int[][] count,Variable var,float dlon){
		float[][] rdata=var.getData()[0][0];
		
		float data=RcToF.apply(record);
		
		if(!isUndef(data)&&cond.test(record)){
			int itag=dd.getXNumPeriodicX(record.getLon(),dlon);
			int jtag=dd.getYNum(record.getLat());
			
			if(rdata[jtag][itag]!=undef) rdata[jtag][itag]+=data;
			else rdata[jtag][itag]=data;
			
			if(count!=null) count[jtag][itag]++;
		}
	}
	
	private static void binningRecordTVariantSum
	(DataDescriptor dd,Record record,Function<Record,Float> RcToF,Predicate<Record> cond,int[][][] count,Variable var,float dlon){
		float[][][] rdata=var.getData()[0];
		
		float data=RcToF.apply(record);
		
		if(!isUndef(data)&&cond.test(record)){
			int itag=dd.getXNumPeriodicX(record.getLon(),dlon);
			int jtag=dd.getYNum(record.getLat());
			int ltag=dd.getTNum(record.getTime());
			
			if(rdata[jtag][itag][ltag]!=undef) rdata[jtag][itag][ltag]+=data;
			else rdata[jtag][itag][ltag]=data;
			
			if(count!=null) count[jtag][itag][ltag]++;
		}
	}
	
	private static void meanTInvariant(int[][] count,boolean maskUndef,Variable var){
		float[][] rdata=var.getData()[0][0];
		
		for(int j=0,J=var.getYCount();j<J;j++)
		for(int i=0,I=var.getXCount();i<I;i++){
			float tmp=count[j][i];
			
			if(tmp!=0) rdata[j][i]/=tmp;
			else if(maskUndef) rdata[j][i]=undef;
		}
	}
	
	private static void meanTVariant(int[][][] count,boolean maskUndef,Variable var){
		float[][][] rdata=var.getData()[0];
		
		for(int j=0,J=var.getYCount();j<J;j++)
		for(int i=0,I=var.getXCount();i<I;i++)
		for(int l=0,L=var.getTCount();l<L;l++){
			float tmp=count[j][i][l];
			
			if(tmp!=0) rdata[j][i][l]/=tmp;
			else if(maskUndef) rdata[j][i][l]=undef;
		}
	}
	
	private static boolean isUndef(float data){
		if(data==undef) return true;
		
		return false;
	}
	
	private static int getPartValue(long time,String part){
		String partLower=part.toLowerCase();
		
		switch(partLower){
		case "year":
			return (int)(time/10000000000L);
		case "month":
			return (int)(time%10000000000L/100000000L);
		case "date":
			return (int)(time%100000000L/1000000L);
		default:
			throw new IllegalArgumentException(
				"unsupported String for \""+partLower+"\", should be [year, month, date]"
			);
		}
	}
	
	
	/** test
	public static void main(String[] args){
		long a=20111210063000L;
		boolean[] b=new boolean[1];
		
		System.out.println(getYear(a));
		System.out.println(getMonth(a));
		System.out.println(getDate(a));
		System.out.println(b[0]);
	}*/
}
