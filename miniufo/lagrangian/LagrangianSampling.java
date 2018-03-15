/**
 * @(#)LagrangianSampling.java	1.0 2018.03.14
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.IOUtil;
import miniufo.io.Print;
import miniufo.util.GridDataFetcher;
import miniufo.util.TicToc;


/**
 * Lagrangian (along-track) sampling of gridded data
 *
 * @version 1.0, 2018.03.14
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LagrangianSampling implements Print{
	//
	private int attachedLen=0;
	
	private boolean print=true;
	
	private List<? extends Particle> ps=null;
	
	
	/**
	 * Constructor
	 */
	public LagrangianSampling(List<? extends Particle> ps){
		this.ps=ps;
		
		attachedLen=ps.get(0).getRecord(0).getDataLength();
		
		LagrangianUtil.asRecordStream(ps).forEach(r->{
			if(r.getDataLength()!=attachedLen)
			throw new IllegalArgumentException("attached data lengths are not equal to "+attachedLen);
		});
	}
	
	
	/**
	 * Sample variables described in a ctl/cts file.
	 * 
	 * @param	ctl		a ctl file
	 * @param	meta	AttachedMeta that need to sample
	 */
	public void sampleVariables(String ctl,AttachedMeta... meta){
		int len=meta.length;
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(ctl);
		DataDescriptor dd=df.getDataDescriptor();
		
		int tLen=ps.get(0).getTCount();
		int tstr=dd.getTNum(ps.get(0).getTime(0     ))+1;
		int tend=dd.getTNum(ps.get(0).getTime(tLen-1))+1;
		
		GridDataFetcher gdf=new GridDataFetcher(dd);
		
		for(int i=0;i<len;i++){
			if(print) TicToc.tic("samping "+meta[i].name+" (attachedIndex: "+meta[i].index+") from t="+tstr+" to t="+tend);
			
			Variable buffer=gdf.prepareXYTBuffer(meta[i].name,1,tstr,tLen,0);
			
			buffer.replaceUndefData(Record.undef);
			
			final int ii=i;
			LagrangianUtil.asRecordStream(ps).forEach(r->{
				float v=gdf.fetchXYTBuffer(r.getXPos(),r.getYPos(),r.getTime(),buffer);
				
				r.setData(meta[ii],v);
			});
			
			for(Particle p:ps) p.setAttachedMeta(meta[i]);
			
			if(print) TicToc.toc(TimeUnit.SECONDS);
		}
		
		gdf.closeFile();
	}
	
	/**
	 * Sample variables at a instant time described in a ctl/cts file.
	 * 
	 * @param	ctl		a ctl file
	 * @param	tstep	time step (start from 1)
	 * @param	meta	AttachedMeta that need to sample
	 */
	public void sampleInstantVariables(String ctl,int tstep,AttachedMeta... meta){
		int len=meta.length;
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(ctl);
		DataDescriptor dd=df.getDataDescriptor();
		
		GridDataFetcher gdf=new GridDataFetcher(dd);
		
		for(int i=0;i<len;i++){
			if(print) TicToc.tic("samping "+meta[i].name+" (attachedIndex: "+meta[i].index+") at t="+tstep);
			
			Variable buffer=gdf.prepareXYBuffer(meta[i].name,tstep,1);
			
			buffer.replaceUndefData(Record.undef);
			
			final int ii=i;
			LagrangianUtil.asRecordStream(ps).forEach(r->{
				float v=gdf.fetchXYBuffer(r.getXPos(),r.getYPos(),buffer);
				
				r.setData(meta[ii],v);
			});
			
			for(Particle p:ps) p.setAttachedMeta(meta[i]);
			
			if(print) TicToc.toc(TimeUnit.SECONDS);
		}
		
		gdf.closeFile();
	}
	
	
	/**
	 * Write Lagrangian sampling data to a GrADS file (including ctl).
	 * 
	 * @param	out		file name for output
	 */
	public void toGrADSFile(String out){
		Variable[] vs=ps.stream().map(p->toTimeSeriesVariable(p)).toArray(Variable[]::new);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(out);
		cdws.writeData(vs); cdws.closeFile();
		
		StringBuilder sb=new StringBuilder();
		
		Particle p=ps.get(0); MDate tstr=new MDate(p.getTime(0));
		
		int dt=tstr.getDT(new MDate(p.getTime(1)));
		
		if(dt%86400!=0) throw new IllegalArgumentException("dt ("+dt+") cannot be divided by 86400");
		
		sb.append("dset ^"+IOUtil.getFileName(out)+"\n");
		sb.append("undef -9999\n");
		sb.append("title along-track sampled data\n");
		sb.append("xdef   1 linear 0 1\n");
		sb.append("ydef   1 linear 0 1\n");
		sb.append("zdef   "+vs[0].getZCount()+" linear 0 1\n");
		sb.append("tdef   "+vs[0].getTCount()+" linear "+tstr.toGradsDate()+" "+dt/86400+"dy\n");
		sb.append("vars "+vs.length+"\n");
		for(int i=0;i<vs.length;i++)
		sb.append(String.format("%-11s "+vs[i].getZCount()+" 99  %s\n","v"+(i+1),vs[i].getName()));
		sb.append("endvars\n");
		
		try(FileWriter fw=new FileWriter(IOUtil.getCompleteFileNameWithoutExtension(out)+".ctl")){ fw.write(sb.toString());}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/**
	 * whether to print out
	 *
     * @param	print	print or disable print
     */ 
	public void setPrinting(boolean print){ this.print=print;}
	
	
	/**
	 * Convert a particle data to a t-z Variable (x=y=1, a single point).
	 * 
	 * @param	p	a single particle
	 */
	private Variable toTimeSeriesVariable(Particle p){
		int t=p.getTCount(),z=p.getRecord(0).getDataLength();
		
		Variable v=new Variable("p"+p.getID(),false,new Range(t,z,1,1));
		
		float[][][][] vdata=v.getData();
		
		for(AttachedMeta meta:p.getAttachedMeta())
		for(int l=0;l<t;l++) vdata[meta.index][0][0][l]=p.getRecord(l).getDataValue(meta);
		
		return v;
	}
	
	
	/** test
	public static void main(String[] args){
		
	}*/
}
