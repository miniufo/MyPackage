/**
 * @(#)BinningStatistics.java	1.0 2017.09.21
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.Stream;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.AttachedMeta;
import miniufo.lagrangian.LagrangianUtil;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;
import miniufo.statistics.StatisticsUtil;


/**
 * A general binning (grouping) class for Eulerian Statistics.
 *
 * @version 1.0, 2017.09.21
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class BinningStatistics{
	//
	private boolean maskUndef=false; // mask zero-count grids with undef
	
	private float gridSize=0;
	
	private DataDescriptor dd=null;
	
	public static final float undef=-9999.0f;
	
	
	/**
	 * Constructor
	 */
	public BinningStatistics(DataDescriptor dd){
		this.dd=dd;
		
		gridSize=(float)Math.toRadians(StatisticsUtil.cArithmeticMean(dd.getDYDef()))*SpatialModel.REarth;
	}
	
	
	/*** getor and setor ***/
	public boolean getMaskUndef(){ return maskUndef;}
	
	public DataDescriptor getDescriptor(){ return dd;}
	
	public void setMaskUndef(boolean maskUndef){ this.maskUndef=maskUndef;}
	
	
	
	/**
     * Binning the types of records.
     *
     * @param	ls		records in a list
     * @param	cond	condition that should be met (e.g., r->r.getDataValue(2)>17.2)
     * @param	types	types that need to be binned
     */
	public Variable[] binningTCTypeCount(List<Typhoon> ls,Predicate<Record> cond,TYPE... types){
		if(types==null||types.length==0) throw new IllegalArgumentException("specify some types");
		
		Variable[] re=new Variable[types.length];
		
		for(int i=0,I=types.length;i<I;i++){
			TYPE t=types[i];
			
			ToDoubleFunction<Record> RcToD=r->1.0;
			Predicate<Record> isType=r->Typhoon.getType(r.getData(Typhoon.Type))==t;
			
			re[i]=binningDataMeanByCondition(LagrangianUtil.asRecordStream(ls),RcToD,isType.and(cond));
			re[i].setName(t.name());
			re[i].setCommentAndUnit("binned type count");
			re[i].setUndef(undef);
		}
		
		return re;
	}
	
	/**
     * Binning the ACE of records.
     *
     * @param	ls		records in a list
     * @param	cond	condition that should be met (e.g., r->r.getDataValue(2)>17.2)
     * @param	types	types that need to be binned
     */
	public Variable binningTCACE(List<Typhoon> ls,Predicate<Record> cond){
		ToDoubleFunction<Record> RcToD=r->r.getData(Typhoon.Vmax)*r.getData(Typhoon.Vmax);
		
		Variable re=binningDataMeanByCondition(LagrangianUtil.asRecordStream(ls),RcToD,cond);
		re.setName("ace");
		re.setCommentAndUnit("accumulated cyclone energy");
		re.setUndef(undef);
		
		return re;
	}
	
	/**
     * Binning the genesis frequency of Typhoons.
     *
     * @param	ls		records in a list
     * @param	cond	condition that should be met (e.g., r->r.getDataValue(2)>17.2)
     * @param	types	types that need to be binned
     */
	public Variable binningTCGenesisFrequency(List<Typhoon> ls,float wthre){
		ToDoubleFunction<Record> RcToD=r->1.0;
		
		Function<Typhoon,Optional<Record>> genesisRec=ty->{
			Record re=null;
			
			for(int l=0,L=ty.getTCount();l<L;l++){
				Record r=ty.getRecord(l);
				if(r.getData(Typhoon.Vmax)>=wthre){ re=r; break;}
			}
			
			return Optional.ofNullable(re);
		};
		
		Variable re=binningDataMeanByCondition(ls.stream().map(genesisRec).filter(op->op.isPresent()).map(op->op.get()),RcToD,r->true);
		re.setName("gen");
		re.setCommentAndUnit("TC genesis frequency");
		re.setUndef(undef);
		
		return re;
	}
	
	
	
	/**
     * Binning particle attached data.
     *
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	maskUndef	mask the field with undef
     * @param	idx			indices of AttachedData
     */
	public Variable[] binningData(List<? extends Particle> ls,AttachedMeta... meta){
		if(meta==null||meta.length==0) throw new IllegalArgumentException("specify some AttachedMeta");
		
		int len=meta.length;
		
		Variable[] re=new Variable[len];
		
		for(int i=0;i<len;i++){
			final int ii=i;
			
			ToDoubleFunction<Record> RcToD=r->r.getData(meta[ii]);
			Predicate<Record> cond=r->r.getData(meta[ii])!=undef;
			
			re[i]=binningDataMeanByCondition(LagrangianUtil.asRecordStream(ls),RcToD,cond);
			re[i].setName("data"+i);
			re[i].setCommentAndUnit(meta[i].name);
			re[i].setUndef(undef);
		}
		
		return re;
	}
	
	/**
     * Binning data by date.
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	idx			indices of AttachedData
     * @param	part		part of date in String form: ["year","month","date"]
     * @param	values		values of the date that would be binned
     */
	public Variable[] binningDataByDate(List<? extends Particle> ls,AttachedMeta[] meta,String part,int... values){
		if(dd.getTCount()!=1) throw new IllegalArgumentException("T-count should be 1");
		if(meta==null||meta.length==0) throw new IllegalArgumentException("specify some indices of AttachedData");
		if(values==null||values.length==0) throw new IllegalArgumentException("specify some values of date");
		
		Variable[] re=new Variable[meta.length];
		
		for(int i=0,I=meta.length;i<I;i++){
			final int ii=i;
			
			ToDoubleFunction<Record> RcToD=r->r.getData(meta[ii]);
			Predicate<Record> cond=r->{
				if(r.getData(meta[ii])==undef) return false;
				
				int value=getPartValue(r.getTime(),part);
				
				for(int v:values) if(v==value) return true;
				
				return false;
			};
			
			re[i]=binningDataMeanByCondition(LagrangianUtil.asRecordStream(ls),RcToD,cond);
			re[i].setName("data");
			re[i].setCommentAndUnit(meta[i].name);
			re[i].setUndef(undef);
		}
		
		return re;
	}
	
	/**
     * Binning data of different seasons.
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	season		e.g., int[][] season=new int[][]{
								{3, 4, 5},	{ 6,7,8},	spring, summer,
								{9,10,11},	{12,1,2}	autumn, winter
							};
     * @param	idx			indices of AttachedData
     * 
     * @return	re			means of different data of different seasons, [seasons][AttachedData]
     */
	public Variable[][] binningSeasonalData(List<? extends Particle> ls,int[][] seasons,AttachedMeta... meta){
		int season=seasons.length;
		int sealen=seasons[0].length;
		
		for(int[] s:seasons)
		if(s.length!=sealen) throw new IllegalArgumentException("season lengths not equal");
		
		if(meta.length==0) throw new IllegalArgumentException("no idx specified");
		
		Variable[][] re=new Variable[season][meta.length];
		
		for(int i=0;i<season;i++){
			Variable[] tmp=binningDataByDate(ls,meta,"month",seasons[i]);
			
			for(int j=0,J=tmp.length;j<J;j++){
				re[i][j]=tmp[j];
				re[i][j].setName(tmp[j].getName()+"_"+i);
			}
		}
		
		return re;
	}
	
	/**
     * Binning count of different seasons.
     * Spring [3, 4, 5]; Summer [6, 7, 8]; Autumn [9, 10, 11]; Winter [12, 1, 2]
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	season		e.g., int[][] season=new int[][]{
	 *							{3, 4, 5},	{ 6,7,8},	spring, summer,
	 *							{9,10,11},	{12,1,2}	autumn, winter
	 *						};
     * 
     * @return	re			count of different seasons
     */
	public Variable[] binningSeasonalCount(List<? extends Particle> ls,int[][] seasons){
		int sc=seasons.length;
		
		Variable[] cs=new Variable[sc];
		
		for(int i=0;i<sc;i++){
			cs[i]=binningCountByDate(ls,"month",seasons[i]);
			cs[i].setName("count_mon"+(i+1));
		}
		
		return cs;
	}
	
	
	/**
     * Binning the data count.
     *
     * @param	dd		a DataDescriptor template used to describe the grids
     * @param	ls		a list of data
     */
	public Variable binningCount(List<? extends Particle> ls){
		Variable re=binningDataSumByCondition(LagrangianUtil.asRecordStream(ls),r->1.0,r->true);
		re.setName("count");
		re.setCommentAndUnit("binned count (1)");
		re.setUndef(undef);
		
		return re;
	}
	
	/**
     * Binning count by date.
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	part		part of date in String form: ["year","month","date"]
     * @param	values		values of the date that would be binned
     */
	public Variable binningCountByDate(List<? extends Particle> ls,String part,int... values){
		int tcount=dd.getTCount();
		
		if(tcount!=1) throw new IllegalArgumentException("T-count should be 1");
		if(values==null||values.length==0) throw new IllegalArgumentException("specify some values of date");
		
		Predicate<Record> cond=r->{
			int value=getPartValue(r.getTime(),part);
			
			for(int v:values) if(v==value) return true;
			
			return false;
		};
		
		Variable re=binningDataSumByCondition(LagrangianUtil.asRecordStream(ls),r->1.0,cond);
		
		re.setName("count");
		re.setCommentAndUnit("binned count for "+part+" = "+Arrays.toString(values));
		re.setUndef(undef);
		
		return re;
	}
	
	/**
     * Binning median count by date.
     * 
     * @param	dd			a DataDescriptor template used to describe the grids
     * @param	ls			a list of records
     * @param	part		part of date in String form: ["year","month","date"]
     * @param	values		values of the date that would be binned
     */
	public Variable binningMedianCount(List<? extends Particle> ls){
		Variable re=binningDataSumByCondition(ls.stream().map(p->p.getRecord(p.getMedianIndex())),r->1.0,r->true);
		re.setName("count");
		re.setCommentAndUnit("binned median-position count (1)");
		re.setUndef(undef);
		
		return re;
	}
	
	
	
	/**
     * Common interface for binning Record data and then averaging them.
     * 
     * @param	ls		a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     */
	public Variable binningDataMeanByCondition(Stream<Record> records,ToDoubleFunction<Record> RcToD,Predicate<Record> cond){
		Variable var=new Variable("ave",new Range(dd.getTCount(),1,dd.getYCount(),dd.getXCount()));
		var.setUndef(undef);
		
		if(dd.getTCount()==1){
			int[][] count=new int[dd.getYCount()][dd.getXCount()];
			
			Consumer<Record> csmr=r->countOneRecordTInvariantSum(r,RcToD,cond,count,var);
			
			records.forEach(csmr);
			
			meanTInvariant(count,var);
			
		}else{
			int[][][] count=new int[dd.getYCount()][dd.getXCount()][dd.getTCount()];
			
			Consumer<Record> csmr=r->countOneRecordTVariantSum(r,RcToD,cond,count,var);
			
			records.forEach(csmr);
			
			meanTVariant(count,var);
		}
		
		return var;
	}
	
	/**
     * Common interface for binning Record data and accumulating them.
     * 
     * @param	ls		a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     */
	public Variable binningDataSumByCondition(Stream<Record> records,ToDoubleFunction<Record> RcToD,Predicate<Record> cond){
		Variable var=new Variable("sum",new Range(dd.getTCount(),1,dd.getYCount(),dd.getXCount()));
		var.setUndef(undef);
		
		Consumer<Record> csmr=(dd.getTCount()==1)?
			r->countOneRecordTInvariantSum(r,RcToD,cond,null,var):
			r->countOneRecordTVariantSum  (r,RcToD,cond,null,var);
		
		records.forEach(csmr);
		
		return var;
	}
	
	/**
     * Common interface for binning maximum Record data.
     * 
     * @param	ls		a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     */
	public Variable binningDataMaxByCondition(Stream<Record> records,ToDoubleFunction<Record> RcToD,Predicate<Record> cond){
		Variable var=new Variable("max",new Range(dd.getTCount(),1,dd.getYCount(),dd.getXCount()));
		var.setUndef(undef);
		
		Consumer<Record> csmr=(dd.getTCount()==1)?
			r->countOneRecordTInvariantMax(r,RcToD,cond,var):
			r->countOneRecordTVariantMax  (r,RcToD,cond,var);
		
		records.forEach(csmr);
		
		return var;
	}
	
	
	/**
     * Common interface for spreading Record data using e-folding model.
     * 
     * @param	ls			a list of Particles containing Records
     * @param	efoldDis	e-folding distance (m)
     * @param	RcToF		a function that maps a Record to a Float as data
     * @param	cond		a condition that the Record is accepted
     */
	public Variable spreadingDataByCondition(Stream<Record> records,float efoldDis,ToDoubleFunction<Record> RcToD,Predicate<Record> cond){
		Variable var=new Variable("sum",new Range(dd.getTCount(),1,dd.getYCount(),dd.getXCount()));
		var.setUndef(undef);
		
		Consumer<Record> csmr=(dd.getTCount()==1)?
			r->spreadOneRecordTInvariantSum(r,efoldDis,RcToD,cond,var):
			r->spreadOneRecordTVariantSum  (r,efoldDis,RcToD,cond,var);
		
		records.forEach(csmr);
		
		return var;
	}
	
	/**
     * Common interface for spreading maximum Record data using e-folding model.
     * 
     * @param	ls			a list of Particles containing Records
     * @param	efoldDis	e-folding distance (m)
     * @param	RcToF		a function that maps a Record to a Float as data
     * @param	cond		a condition that the Record is accepted
     */
	public Variable spreadingDataMaxByCondition(Stream<Record> records,float efoldDis,ToDoubleFunction<Record> RcToD,Predicate<Record> cond){
		Variable var=new Variable("max",new Range(dd.getTCount(),1,dd.getYCount(),dd.getXCount()));
		var.setUndef(undef);
		
		Consumer<Record> csmr=(dd.getTCount()==1)?
			r->spreadOneRecordTInvariantMax(r,efoldDis,RcToD,cond,var):
			r->spreadOneRecordTVariantMax  (r,efoldDis,RcToD,cond,var);
		
		records.forEach(csmr);
		
		return var;
	}
	
	
	
	/*** helper methods ***/
	
	/**
	 * Counting a single Record accumulatively for t-invariant grid template.
	 * 
     * @param	record	a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     * @param	count	count array for each grid
     * @param	var		variable need to store the result
	 */
	private void countOneRecordTInvariantSum(Record record,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,int[][] count,Variable var){
		float[][] rdata=var.getData()[0][0];
		
		float data=(float)RcToD.applyAsDouble(record);
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			
			if(itag!=-1&&jtag!=-1){ // within the domain of grid template
				if(rdata[jtag][itag]!=undef) rdata[jtag][itag]+=data;
				else rdata[jtag][itag]=data;
				
				if(count!=null) count[jtag][itag]++;
			}
		}
	}
	
	/**
	 * Counting a single Record accumulatively for t-variant grid template.
	 * 
     * @param	record	a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     * @param	count	count array for each grid
     * @param	var		variable need to store the result
	 */
	private void countOneRecordTVariantSum(Record record,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,int[][][] count,Variable var){
		float[][][] rdata=var.getData()[0];
		
		float data=(float)RcToD.applyAsDouble(record);
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			int ltag=dd.getTNum(record.getTime());
			
			if(itag!=-1&&jtag!=-1&&ltag!=-1){ // within the domain of grid template
				if(rdata[jtag][itag][ltag]!=undef) rdata[jtag][itag][ltag]+=data;
				else rdata[jtag][itag][ltag]=data;
				
				if(count!=null) count[jtag][itag][ltag]++;
			}
		}
	}
	
	/**
	 * Spreading a single Record using e-folding model for t-invariant grid template.
	 * 
     * @param	record		a list of Particles containing Records
     * @param	efoldDis	e-folding distance (m)
     * @param	RcToF		a function that maps a Record to a Float as data
     * @param	cond		a condition that the Record is accepted
     * @param	var			variable need to store the result
	 */
	private void spreadOneRecordTInvariantSum(Record record,float efoldDis,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,Variable var){
		int x=dd.getXCount(),y=dd.getYCount(),spreadRad=Math.round(efoldDis*4f/gridSize)+1;
		
		double data=RcToD.applyAsDouble(record);
		
		float[] xdef=dd.getXDef().getSamples();
		float[] ydef=dd.getYDef().getSamples();
		
		float[][] rdata=var.getData()[0][0];
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			
			if(itag!=-1&&jtag!=-1) // within the domain of grid template
			for(int j=jtag-spreadRad,J=jtag+spreadRad;j<=J;j++) if(j>=0&&j<y){
			for(int i=itag-spreadRad,I=itag+spreadRad;i<=I;i++) if(i>=0&&i<x){
				float dis=SpatialModel.cSphericalDistanceByDegree(record.getXPos(),record.getYPos(),xdef[i],ydef[j]);
				float val=(float)(data*Math.exp(-Math.abs(dis/efoldDis)));
				
				if(  rdata[j][i]!=undef) rdata[j][i]+=val;
				else rdata[j][i] =val;
			}
			}
		}
	}
	
	/**
	 * Spreading a single Record using e-folding model for t-variant grid template.
	 * 
     * @param	record		a list of Particles containing Records
     * @param	efoldDis	e-folding distance (m)
     * @param	RcToF		a function that maps a Record to a Float as data
     * @param	cond		a condition that the Record is accepted
     * @param	var			variable need to store the result
	 */
	private void spreadOneRecordTVariantSum(Record record,float efoldDis,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,Variable var){
		int x=dd.getXCount(),y=dd.getYCount(),spreadRad=Math.round(efoldDis*4f/gridSize)+1;
		
		double data=RcToD.applyAsDouble(record);
		
		float[] xdef=dd.getXDef().getSamples();
		float[] ydef=dd.getYDef().getSamples();
		
		float[][][] rdata=var.getData()[0];
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			int ltag=dd.getTNum(record.getTime());
			
			if(itag!=-1&&jtag!=-1&&ltag!=-1) // within the domain of grid template
			for(int j=jtag-spreadRad,J=jtag+spreadRad;j<=J;j++) if(j>=0&&j<y){
			for(int i=itag-spreadRad,I=itag+spreadRad;i<=I;i++) if(i>=0&&i<x){
				float dis=SpatialModel.cSphericalDistanceByDegree(record.getXPos(),record.getYPos(),xdef[i],ydef[j]);
				float val=(float)(data*Math.exp(-Math.abs(dis/efoldDis)));
				
				if(  rdata[j][i][ltag]!=undef) rdata[j][i][ltag]+=val;
				else rdata[j][i][ltag] =val;
			}
			}
		}
	}
	
	
	/**
	 * Counting a single Record maximumly for t-invariant grid template.
	 * 
     * @param	record	a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     * @param	count	count array for each grid
     * @param	var		variable need to store the result
	 */
	private void countOneRecordTInvariantMax(Record record,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,Variable var){
		float[][] rdata=var.getData()[0][0];
		
		float data=(float)RcToD.applyAsDouble(record);
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			
			if(itag!=-1&&jtag!=-1){ // within the domain of grid template
				if(rdata[jtag][itag]!=undef) rdata[jtag][itag]=Math.max(rdata[jtag][itag],data);
				else rdata[jtag][itag]=data;
			}
		}
	}
	
	/**
	 * Counting a single Record maximumly for t-variant grid template.
	 * 
     * @param	record	a list of Particles containing Records
     * @param	RcToF	a function that maps a Record to a Float as data
     * @param	cond	a condition that the Record is accepted
     * @param	count	count array for each grid
     * @param	var		variable need to store the result
	 */
	private void countOneRecordTVariantMax(Record record,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,Variable var){
		float[][][] rdata=var.getData()[0];
		
		float data=(float)RcToD.applyAsDouble(record);
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			int ltag=dd.getTNum(record.getTime());
			
			if(itag!=-1&&jtag!=-1&&ltag!=-1){ // within the domain of grid template
				if(rdata[jtag][itag][ltag]!=undef) rdata[jtag][itag][ltag]=Math.max(rdata[jtag][itag][ltag],data);
				else rdata[jtag][itag][ltag]=data;
			}
		}
	}
	
	/**
	 * Spreading a single Record maximumly using e-folding model for t-invariant grid template.
	 * 
     * @param	record		a list of Particles containing Records
     * @param	efoldDis	e-folding distance (m)
     * @param	RcToF		a function that maps a Record to a Float as data
     * @param	cond		a condition that the Record is accepted
     * @param	var			variable need to store the result
	 */
	private void spreadOneRecordTInvariantMax(Record record,float efoldDis,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,Variable var){
		int x=dd.getXCount(),y=dd.getYCount(),spreadRad=Math.round(efoldDis*4f/gridSize)+1;
		
		double data=RcToD.applyAsDouble(record);
		
		float[] xdef=dd.getXDef().getSamples();
		float[] ydef=dd.getYDef().getSamples();
		
		float[][] rdata=var.getData()[0][0];
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			
			if(itag!=-1&&jtag!=-1) // within the domain of grid template
			for(int j=jtag-spreadRad,J=jtag+spreadRad;j<=J;j++) if(j>=0&&j<y){
			for(int i=itag-spreadRad,I=itag+spreadRad;i<=I;i++) if(i>=0&&i<x){
				float dis=SpatialModel.cSphericalDistanceByDegree(record.getXPos(),record.getYPos(),xdef[i],ydef[j]);
				float val=(float)(data*Math.exp(-Math.abs(dis/efoldDis)));
				
				if(  rdata[j][i]!=undef) rdata[j][i]=Math.max(rdata[j][i],val);
				else rdata[j][i] =val;
			}
			}
		}
	}
	
	/**
	 * Spreading a single Record maximumly using e-folding model for t-variant grid template.
	 * 
     * @param	record		a list of Particles containing Records
     * @param	efoldDis	e-folding distance (m)
     * @param	RcToF		a function that maps a Record to a Float as data
     * @param	cond		a condition that the Record is accepted
     * @param	var			variable need to store the result
	 */
	private void spreadOneRecordTVariantMax(Record record,float efoldDis,ToDoubleFunction<Record> RcToD,Predicate<Record> cond,Variable var){
		int x=dd.getXCount(),y=dd.getYCount(),spreadRad=Math.round(efoldDis*4f/gridSize)+1;
		
		double data=RcToD.applyAsDouble(record);
		
		float[] xdef=dd.getXDef().getSamples();
		float[] ydef=dd.getYDef().getSamples();
		
		float[][][] rdata=var.getData()[0];
		
		if(data!=undef&&cond.test(record)){
			int itag=dd.isPeriodicX()?dd.getXNumPeriodicX(record.getXPos()):dd.getXNum(record.getXPos());
			int jtag=dd.getYNum(record.getYPos());
			int ltag=dd.getTNum(record.getTime());
			
			if(itag!=-1&&jtag!=-1&&ltag!=-1) // within the domain of grid template
			for(int j=jtag-spreadRad,J=jtag+spreadRad;j<=J;j++) if(j>=0&&j<y){
			for(int i=itag-spreadRad,I=itag+spreadRad;i<=I;i++) if(i>=0&&i<x){
				float dis=SpatialModel.cSphericalDistanceByDegree(record.getXPos(),record.getYPos(),xdef[i],ydef[j]);
				float val=(float)(data*Math.exp(-Math.abs(dis/efoldDis)));
				
				if(  rdata[j][i][ltag]!=undef) rdata[j][i][ltag]=Math.max(rdata[j][i][ltag],val);
				else rdata[j][i][ltag] =val;
			}
			}
		}
	}
	
	
	/**
	 * Averaging the result.
	 * 
     * @param	count		count array for each grid
     * @param	maskUndef	mask those zero-count grid with undefined value
     * @param	var			variable need to store the result
	 */
	private void meanTInvariant(int[][] count,Variable var){
		float[][] rdata=var.getData()[0][0];
		
		for(int j=0,J=var.getYCount();j<J;j++)
		for(int i=0,I=var.getXCount();i<I;i++){
			float tmp=count[j][i];
			
			if(tmp!=0) rdata[j][i]/=tmp;
			else if(maskUndef) rdata[j][i]=undef;
		}
	}
	
	/**
	 * Averaging the result.
	 * 
     * @param	count		count array for each grid
     * @param	maskUndef	mask those zero-count grid with undefined value
     * @param	var			variable need to store the result
	 */
	private void meanTVariant(int[][][] count,Variable var){
		float[][][] rdata=var.getData()[0];
		
		for(int j=0,J=var.getYCount();j<J;j++)
		for(int i=0,I=var.getXCount();i<I;i++)
		for(int l=0,L=var.getTCount();l<L;l++){
			float tmp=count[j][i][l];
			
			if(tmp!=0) rdata[j][i][l]/=tmp;
			else if(maskUndef) rdata[j][i][l]=undef;
		}
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
