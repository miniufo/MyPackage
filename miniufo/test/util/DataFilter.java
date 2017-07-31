/**
 * @(#)DataMean.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;

import miniufo.application.statisticsModel.FilterMethods;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.DataDescriptor.DataType;
import miniufo.descriptor.Var;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import miniufo.io.DataWrite;


/**
 * mean process of the data
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DataFilter{
	//
	private DataDescriptor dd=null;
	
	
	/**
     * constructor
     *
     * @param	src		source DataDescriptor
     */
	public DataFilter(DataDescriptor src){
		dd=src;
		/*
		int z=0;
		
		for(Var v:dd.getVDef()) z+=v.getZCount();
		
		fileSize=dd.getXCount()*dd.getYCount()/1024f*dd.getTCount()/1024f*z/250f;
		
		System.out.println("\n size of filtering file: "+fileSize+" GB ("+fileSize/z+" GB per slice)");
		
		if(fileSize/z>1f) scratch=true;*/
	}
	
	
	/**
     * remove annual cycle
     *
     * @param	path	path of file after averaging
     */
	public void annualCycleFilter(String path){
		System.out.println("\nStart removing annual cycle...");
		
		DataRead  dr=DataIOFactory.getDataRead(dd);
		DataWrite dw=DataIOFactory.getDataWrite(dd,path);
		
		for(int m=0,M=dd.getVCount();m<M;m++){
			Var var=dd.getVDef()[m];
			
			int t=var.getTCount(),y=var.getYCount(),x=var.getXCount(),z=var.getZCount();
			
			Variable tmp=new Variable(var.getName(),false,new Range(t,1,y,x));
			tmp.setCommentAndUnit(var.getComment());
			tmp.setUndef(var.getUndef());
			
			for(int k=0;k<z;k++){
				int[] zrange=tmp.getRange().getZRange();
				
				zrange[0]=zrange[1]=k+1;
				
				dr.readData(tmp);
				
				DataType dt=getDataType(dd.getTIncrement());
				
				switch(dt){
				case MONTHLY:
					FilterMethods.cycleFilter(tmp,12);
					break;
				case SEASONAL:
					FilterMethods.cycleFilter(tmp,4);
					break;
				case DAILY:
					//FilterMethods.annualCycleFilterForDailyData(tmp,dd.getTDef().getSamples()[0].getYear());
					break;
				case DAILY4:
					//FilterMethods.annualCycleFilterForDaily4Data(tmp,dd.getTDef().getSamples()[0].getYear());
					break;

				default:
					throw new IllegalArgumentException("unsupported data type for annual cycle filter: "+dt);
				}
				
				dw.writeData(tmp);
			}
		}
		
		dw.closeFile();	dr.closeFile();
		
		System.out.println("Finish filtering.");
	}
	
	
	/*** helper methods and class ***/
	private static DataType getDataType(String incre){
		String unit=incre.substring(incre.length()-2,incre.length()).toLowerCase();
		
		if("yr".equals(unit)) return DataType.ANNUAL;
		else if("mo".equals(unit)) return DataType.MONTHLY;
		else if("dy".equals(unit)) return DataType.DAILY;
		else if("hr".equals(unit)) return DataType.DAILY4;
		else throw new IllegalArgumentException("unsupported data type for "+unit);
	}
	
	
	/** test
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/ERAInterim/Data.ctl");
		
		DataFilter dfilter=new DataFilter(df.getDataDescriptor());
		dfilter.annualCycleFilter("d:/Data/ERAInterim/filter.dat");
	}*/
}
