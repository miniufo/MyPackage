/**
 * @(#)WaveNumberFrequencyAnalysis.java	1.0 2014.08.04
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import miniufo.diagnosis.Variable;
import miniufo.mathsphysics.Complex;
import miniufo.mathsphysics.FastFourier;


/**
 * wave number-frequency analysis
 *
 * @version 1.0, 2014.08.04
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class WaveNumberFrequencyAnalysis{
	//
	private int tlen=0;
	private int xlen=0;
	
	private Complex[][] res=null;
	
	private Variable v=null;
	
	
	/**
     * constructor
     */
	public WaveNumberFrequencyAnalysis(Variable v){
		checkDimension(v);
		
		tlen=v.getTCount();
		xlen=v.getXCount();
		
		res=new Complex[xlen][tlen];
		
		this.v=v;
	}
	
	
	public void windowingX(float[] window){
		if(window.length!=xlen)
		throw new IllegalArgumentException("window length ("+window.length+") should be "+xlen);
		
		float[][] vdata=v.getData()[0][0];
		
		for(int i=0;i<xlen;i++)
		for(int l=0;l<tlen;l++) vdata[i][l]*=window[i];
	}
	
	public void windowingT(float[] window){
		if(window.length!=tlen)
		throw new IllegalArgumentException("window length ("+window.length+") should be "+tlen);
		
		float[][] vdata=v.getData()[0][0];
		
		for(int i=0;i<xlen;i++)
		for(int l=0;l<tlen;l++) vdata[i][l]*=window[l];
	}
	
	
	public void transform(){
		Complex[]   tmp=new Complex[xlen];
		Complex[][] buf=new Complex[xlen][];
		
		float[][] vdata=v.getData()[0][0];
		
		for(int i=0;i<xlen;i++) buf[i]=FastFourier.fft(vdata[i]);
		
		for(int l=0;l<tlen;l++){
			for(int i=0;i<xlen;i++) tmp[i]=buf[i][l];
			
			Complex[] re=FastFourier.fft(tmp);
			
			for(int i=0;i<xlen;i++) res[i][l]=re[i];
		}
	}
	
	
	/*** getor and setor ***/
	public Variable getMode(){
		Variable mod=new Variable("mod"+v.getName(),v);
		
		float[][] mdata=mod.getData()[0][0];
		
		for(int i=0;i<xlen;i++)
		for(int l=0;l<tlen;l++){
			mdata[i][l]=res[i][l].getMod();
		}
		
		return mod;
	}
	
	public Variable getRealPart(){
		Variable re=new Variable("real"+v.getName(),v);
		
		float[][] rdata=re.getData()[0][0];
		
		for(int i=0;i<xlen;i++)
		for(int l=0;l<tlen;l++){
			rdata[i][l]=res[i][l].getReal();
		}
		
		return re;
	}
	
	public Variable getImagPart(){
		Variable im=new Variable("imag"+v.getName(),v);
		
		float[][] idata=im.getData()[0][0];
		
		for(int i=0;i<xlen;i++)
		for(int l=0;l<tlen;l++){
			idata[i][l]=res[i][l].getImag();
		}
		
		return im;
	}
	
	
	/*** helper methods ***/
	private void checkDimension(Variable v){
		if(v.getYCount()!=1)
		throw new IllegalArgumentException("y count should 1");
		
		if(v.getZCount()!=1)
		throw new IllegalArgumentException("z count should 1");
		
		if(v.getTCount()==1)
		throw new IllegalArgumentException("t count should be larger than 1");
		
		if(v.getXCount()==1)
		throw new IllegalArgumentException("x count should be larger than 1");
		
		if(v.isTFirst())
		throw new IllegalArgumentException("variable should not be T-first");
	}
	
	
	/** test
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/Validate/WaveNoFreq/hgt.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable hgt=df.getVariables(new Range("t(1,14600)",dd),"hgt")[0];
		
		float mean=hgt.averageAlong(Dimension.T).averageAlong(Dimension.X).getData()[0][0][0][0];
		
		System.out.println(mean);
		
		hgt.minusEq(mean);
		
		FilterMethods.FourierFilter(hgt,Dimension.T,1460);
		FilterMethods.FourierFilter(hgt,Dimension.T,730);
		FilterMethods.FourierFilter(hgt,Dimension.T,365);
		FilterMethods.removeLinearTrend(hgt);
		
		Variable h=new Variable("hgt",false,new Range(14000,1,1,hgt.getXCount()));
		
		for(int i=0;i<hgt.getXCount();i++)
		System.arraycopy(hgt.getData()[0][0][i],0,h.getData()[0][0][i],0,h.getTCount());
		
		WaveNumberFrequencyAnalysis wnfa=new WaveNumberFrequencyAnalysis(h);
		wnfa.windowingT(WindowFunction.tukey(h.getTCount(),0.05f));
		wnfa.transform();
		
		Variable mod =wnfa.getMode();
		Variable real=wnfa.getRealPart();
		Variable imag=wnfa.getImagPart();
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Data/Validate/WaveNoFreq/wave.dat");
		dw.writeData(dd,h,mod,real,imag);	dw.closeFile();
	}*/
}
