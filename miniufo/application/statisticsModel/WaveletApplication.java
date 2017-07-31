/**
 * @(#)WaveletApplication.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.mathsphysics.Complex;
import miniufo.mathsphysics.Wavelet;


/**
 * wavelet analysis application
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class WaveletApplication{
	//
	//private float variance;
	
	private float[] zdef=null;
	
	private Variable re =null;	// real part
	private Variable mo =null;	// mode
	private Variable gws=null;	// global wavelet spectrum
	private Variable coi=null;	// cone of influence
	
	private Wavelet  wl =null;
	
	
	/**
     * constructor
     */
	public WaveletApplication(Variable v){
		int t=v.getTCount(),z=v.getZCount(),y=v.getYCount(),x=v.getXCount();
		
		if(t<8) throw new IllegalArgumentException("T Dimension is too short");
		if(z!=y||z!=x||z!=1) throw new IllegalArgumentException("only accepted time series of data");
		
		/*** copy data into a buffer array 'data' ***/
		float[] data=null;	// buffer array
		float[][][][] vdata=v.getData();
		if(v.isTFirst()){
			data=new float[t];
			
			for(int l=0;l<t;l++)data[l]=vdata[l][0][0][0];
			
		}else data=vdata[0][0][0];
		
		for(int l=0;l<t;l++)
		if(data[l]==v.getUndef()) throw new IllegalArgumentException("time series contain undefined value");
		
		/*** get variance of the series ***/
		//miniufo.statistics.StatisticsModel.standardize(data);	// standardize
		//variance=miniufo.statistics.StatisticsModel.cVariance(data);
		
		/*** wavelet analysis ***/
		wl=new Wavelet(data,"morlet");
		Complex[][] wave=wl.getWaveCoefficient();
		
		/*** restore the result ***/
		int m=wave.length;
		Range range1=v.getRange();			Range range2=new Range(t,m,1,1);
		Range range3=new Range(t,1,1,1);	Range range4=new Range(1,m,1,1);
		
		re =new Variable("re" ,false,range2);	 re.setUndef(v.getUndef());
		mo =new Variable("mo" ,false,range2);	 mo.setUndef(v.getUndef());
		coi=new Variable("coi",false,range3);	coi.setUndef(v.getUndef());
		gws=new Variable("gws",false,range4);	gws.setUndef(v.getUndef());
		
		 re.setCommentAndUnit("real part of the wavelet");
		 mo.setCommentAndUnit("mode of the wavelet");
		coi.setCommentAndUnit("cone of influence");
		gws.setCommentAndUnit("global wavelet spectrum");
		
		float[][][][]  redata= re.getData();	float[][][][]  modata= mo.getData();
		float[][][][] coidata=coi.getData();	float[][][][] gwsdata=gws.getData();
		
		for(int l=0;l<t;l++)
		for(int k=0;k<m;k++){
			redata[k][0][0][l]=wave[m-k-1][l].getReal();// reverse and store
			modata[k][0][0][l]=wave[m-k-1][l].getMod();	// reverse and store
		}
		
		coidata[0][0][0]=wl.getConeOfInfluence();
		
		for(int k=0;k<m;k++) gwsdata[k][0][0][0]=wl.getGlobalWaveletSpectrum()[m-1-k];	// reverse and store
		
		/*** calculate zdef ***/
		zdef=wl.getPeriods().clone();
		miniufo.basic.ArrayUtil.reverse(zdef);	// reverse and store
		
		/*** set range ***/
		range2.setTRange(range1);	range2.setYRange(range1);	range2.setXRange(range1);
		range3.setTRange(range1);	range3.setYRange(range1);	range3.setXRange(range1);
		range4.setTRange(range1);	range4.setYRange(range1);	range4.setXRange(range1);
		
		range2.getZRange()[0]=range1.getZRange()[0];	range2.getZRange()[1]=m+1-range1.getZRange()[0];
		range3.getZRange()[0]=range1.getZRange()[0];	range3.getZRange()[1]=range1.getZRange()[0];
		range4.getZRange()[0]=range1.getZRange()[0];	range4.getZRange()[1]=m+1-range1.getZRange()[0];
	}
	
	
	/*** getor and setor ***/
	public float[]  getZDef(){ return zdef;}
	
	public Variable getReal(){ return re  ;}
	
	public Variable  getMod(){ return mo  ;}
	
	public Variable  getCOI(){ return coi ;}
	
	public Variable  getGWS(){ return gws ;}
	
	public Variable getCHISig(float level){
		// regular chi-square, > 1 means the power is significance
		Variable chisig=new Variable("chisig",mo);
		
		chisig.setCommentAndUnit("chi-square significance");
		
		float[] csig=wl.getChiSquareSignificance(level);
		float[][][][] sigdata=chisig.getData();
		
		for(int k=0,K=chisig.getZCount();k<K;k++)
		for(int l=0,L=chisig.getTCount();l<L;l++){
			float mod=mo.getData()[k][0][0][l];
			sigdata[k][0][0][l]=mod*mod/csig[K-1-k];	// reverse and store
		}
		
		return chisig;
	}
	
	public Variable getGWSSig(float level){
		Variable gwssig=new Variable("gwsig",gws);
		
		gwssig.setCommentAndUnit("global wavelet spectrum significance");
		
		float[] gsig=wl.getGlobalSignificance(level);
		float[][][][] sigdata=gwssig.getData();
		
		for(int k=0,K=gwssig.getZCount();k<K;k++)
			sigdata[k][0][0][0]=gsig[K-1-k];	// reverse and store
		
		return gwssig;
	}
	
	/**
	public Variable getSADSig(float level,float s1,float s2){
		if(s1<1||s1>s2) throw new IllegalArgumentException("s1 and s2 must satisfy: 1 < s1 < s2");
		
		int count=0;	float[] Sj=wl.getScales();
		ArrayList<Integer> tag=new ArrayList<Integer>();
		
		for(int j=0,J=zdef.length;j<J;j++)
		if(s1<=Sj[j]&&Sj[j]>=s2){ tag.add(j); count++;}
		
		if(count==0) throw new IllegalArgumentException("No valid scales between "+s1+" and "+s2);
		
		// scale averaged significance
		Variable sadsig=new Variable("sadsig",coi);
		
		sadsig.setComment("scale averaged spectrum significance");
		
		float ssig=wl.getScaleAveragedSignificance(level,s1,s2);
		float[][] scale_avg=new float[zdef.length][mo.getTCount()];
		float[][][][] sigdata=sadsig.getData();
		
		for(int k=0,K=mo.getZCount();k<K;k++)
		for(int l=0,L=mo.getTCount();l<L;l++){
			float mod=mo.getData()[k][0][0][l];
			scale_avg[k][l]=mod*mod/Sj[K-1-k];	// Equation (24)
		}
		
		for(int l=0,L=mo.getTCount();l<L;l++){
			for(int k=0;k<count;k++)
			sigdata[0][0][0][l]+=scale_avg[tag.get(k)][l];	// Equation (24)
			
			sigdata[0][0][0][l]*=variance;
		}
		
		return sadsig;
	}*/
	
	public void writeGS(String path){
		File f=new File(path);
		StringBuilder sb=new StringBuilder();
		
		sb.append("* open the real part and mode file\n");
		sb.append("'open "+path.replaceFirst("\\.gs",".ctl")+"'\n");
		sb.append("* open the global wavelet spectrum file\n");
		sb.append("'open "+f.getParent().replaceAll("\\\\","/")+"/gws.ctl'\n\n");
		
		sb.append("* draw real part\n");
		sb.append("'setvpage 2 1 2 1'\n");
		sb.append("'set t 1 "+mo.getTCount()+"'\n");
		sb.append("'set z "+(int)(1+mo.getZCount()/10f)+" "+mo.getZCount()+"'\n");
		sb.append("'set parea 1 8 0.5 4'\n");
		sb.append("'set xlopts 1 7 0.16'\n");
		sb.append("'set ylopts 1 7 0.16'\n");
		sb.append("'set cthick 10'\n");
		sb.append("'set zlog on'\n");
		sb.append("'set cthick 6'\n");
		sb.append("'set ylevs 4 8 16 32 64 128 256 512 1024 2048'\n");
		sb.append("'set gxout contour'\n");
		sb.append("'set clab off'\n");
		sb.append("'d re'\n");
		sb.append("* draw coi\n");
		sb.append("'set parea 1 8 0.5 4'\n");
		sb.append("'set z 1'\n");
		sb.append("'set grid off'\n");
		sb.append("'set yflip on'\n");
		sb.append("'set gxout line'\n");
		sb.append("'set cmark 0'\n");
		sb.append("'set ccolor 1'\n");
		sb.append("'set cthick 16'\n");
		sb.append("'set ylab off'\n");
		sb.append("'set cstyle 4'\n");
		sb.append("'d log(coi)'\n\n");
		
		sb.append("* draw mode\n");
		sb.append("'setvpage 2 1 1 1'\n");
		sb.append("'set t 1 "+mo.getTCount()+"'\n");
		sb.append("'set z "+(int)(1+mo.getZCount()/10f)+" "+mo.getZCount()+"'\n");
		sb.append("'set grid on'\n");
		sb.append("'set ylab on'\n");
		sb.append("'set parea 1 8 0.5 4'\n");
		sb.append("'set xlopts 1 7 0.16'\n");
		sb.append("'set ylopts 1 7 0.16'\n");
		sb.append("'set zlog on'\n");
		sb.append("'set ylevs 4 8 16 32 64 128 256 512 1024 2048'\n");
		sb.append("'set gxout shaded'\n");
		sb.append("'d mo'\n");
		sb.append("* draw sig > 1\n");
		sb.append("'set gxout contour'\n");
		sb.append("'set clevs 1'\n");
		sb.append("'set clab off'\n");
		sb.append("'set cthick 20'\n");
		sb.append("'d chisig'\n");
		sb.append("'set cmin 2'\n");
		sb.append("'set cthick 1'\n");
		sb.append("'set clab on'\n");
		sb.append("'set clskip 2'\n");
		sb.append("'set ccolor 1'\n");
		sb.append("'d mo'\n");
		sb.append("* draw coi\n");
		sb.append("'set parea 1 8 0.5 4'\n");
		sb.append("'set z 1'\n");
		sb.append("'set grid off'\n");
		sb.append("'set yflip on'\n");
		sb.append("'set gxout line'\n");
		sb.append("'set cmark 0'\n");
		sb.append("'set ccolor 1'\n");
		sb.append("'set cthick 16'\n");
		sb.append("'set ylab off'\n");
		sb.append("'set cstyle 4'\n");
		sb.append("'d log(coi)'\n\n");
		
		sb.append("* draw global wavelet significance\n");
		sb.append("'setvpage 2 1 1 1'\n");
		sb.append("'set parea 8.2 10.8 0.5 4'\n");
		sb.append("'set z "+(int)(1+mo.getZCount()/10f)+" "+mo.getZCount()+"'\n");
		sb.append("'set t 1'\n");
		sb.append("'set yflip off'\n");
		sb.append("'set ylab off'\n");
		sb.append("'set xlopts 1 7 0.16'\n");
		sb.append("'set cthick 12'\n");
		sb.append("'set zlog on'\n");
		sb.append("'set cmark 0'\n");
		sb.append("'set ccolor 1'\n");
		sb.append("'set cstyle 1'\n");
		sb.append("'set xlevs -3 -2 -1 0 1'\n");
		sb.append("'d log10(gws.2)'\n");
		sb.append("'set cmark 0'\n");
		sb.append("'set ccolor 2'\n");
		sb.append("'set cstyle 3'\n");
		sb.append("'d log10(gwsig.2)'\n\n");
		
		sb.append("'enable print "+path.replaceFirst("\\.gs",".gmf")+"'\n");
		sb.append("'print'\n");
		sb.append("'disable print'\n");
		sb.append("'close 2'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		
		try(FileWriter fw=new FileWriter(f,false)){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			miniufo.descriptor.NetCDFDescriptor ctl=
				new miniufo.descriptor.NetCDFDescriptor("D:/Data/GlobalTC/sst.mnmean.nc");
			
			Variable sst=new Variable("sst",new Range("lon(150,270);lat(-6,6);t(1,1848)",ctl));
			
			miniufo.io.NetCDFReadStream nc=new miniufo.io.NetCDFReadStream(ctl);
			nc.readData(sst);	nc.closeFile();
			
			//float[][][][] sstdata=sst.getData();
			//java.util.Scanner sn=new java.util.Scanner(new java.io.File(
			// "C:/MATLAB/R2008a/work/wavelet/sst_nino3.dat"
			//));
			//for(int i=0;i<504;i++) sstdata[0][0][0][i]=sn.nextFloat();
			
			Variable ssta=sst.anomalizeArea();
			miniufo.application.statisticsModel.FilterMethods.timeAveFilter(ssta,12);
			
			//for(float f:ssta.getData()[0][0][0]) System.out.println(f);
			
			WaveletApplication wa=new WaveletApplication(ssta);
			
			Variable re=wa.getReal();
			Variable mo=wa.getMod();
			Variable coi=wa.getCOI();
			Variable gws=wa.getGWS();
			Variable csig=wa.getCHISig(0.95f);
			Variable gsig=wa.getGWSSig(0.95f);
			
			miniufo.io.CtlDataWriteStream cdws=new miniufo.io.CtlDataWriteStream("d:/wlre.dat");
			cdws.writeData(coi,re,mo,csig);	cdws.writeCtl(ctl,wa.getZDef(),null);	cdws.closeFile();
			
			wa.writeGS("d:/wlre.gs");
			
			cdws=new miniufo.io.CtlDataWriteStream("d:/gws.dat");
			cdws.writeData(gws,gsig);	cdws.writeCtl(ctl,wa.getZDef(),null);	cdws.closeFile();
			
		}catch(Exception ex){
			ex.printStackTrace();
			System.exit(0);
		}
	}*/
}