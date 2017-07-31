/**
 * @(#)DataMean.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

import miniufo.io.DataRead;
import miniufo.io.DataWrite;
import miniufo.io.DataIOFactory;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.Var;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * mean process of the data
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DataMean{
	//
	private DataDescriptor ctl=null;
	
	
	/**
     * constructor
     *
     * @param	src		source DataDescriptor
     */
	public DataMean(DataDescriptor src){ ctl=src;}
	
	
	/**
     * pentad mean, daily data only
     *
     * @param	path	path of file after averaging
     */
	public void pentadMean(String path){
		System.out.println("\nStart data pentad averaging...");
		
		int t=ctl.getTCount(),z=ctl.getZCount(),y=ctl.getYCount(),x=ctl.getXCount(),v=ctl.getVCount();
		
		if(t<365)
			throw new IllegalArgumentException("daily data only, length should larger than 365");
		
		int syear=ctl.getTDef().getSamples()[0].getYear();
		int yearc=t/365;
		
		Var[] cv=ctl.getVDef();
		
		int[] dcounts=new int[yearc];
		int[] stags  =new int[yearc];
		
		float[] ltmp=new float[366];
		float[] ptmp=new float[365];
		
		if(MDate.isLeapYear(syear)) dcounts[0]=366;
		else dcounts[0]=365;
		
		/*** calculate the tags ***/
		for(int i=syear+1;i<=syear+yearc-1;i++)
		if(MDate.isLeapYear(i)){
			dcounts[i-syear]=366;
			stags[i-syear]=stags[i-1-syear]+365;
			
		}else{
			dcounts[i-syear]=365;
			if(MDate.isLeapYear(i-1)) stags[i-syear]=stags[i-1-syear]+366;
			else stags[i-syear]=stags[i-1-syear]+365;
		}
		
		final int[][] pp={
			{  0,  4},{  5,  9},{ 10, 14},{ 15, 19},{ 20, 24},{ 25, 30},
			{ 31, 35},{ 36, 40},{ 41, 45},{ 46, 50},{ 51, 55},{ 56, 58},
			{ 59, 63},{ 64, 68},{ 69, 73},{ 74, 78},{ 79, 83},{ 84, 89},
			{ 90, 94},{ 95, 99},{100,104},{105,109},{110,114},{115,119},
			{120,124},{125,129},{130,134},{135,139},{140,144},{145,150},
			{151,155},{156,160},{161,165},{166,170},{171,175},{176,180},
			{181,185},{186,190},{191,195},{196,200},{201,206},{207,211},
			{212,216},{217,221},{222,226},{227,231},{232,236},{237,242},
			{243,247},{248,252},{253,257},{258,262},{263,267},{268,272},
			{273,277},{278,282},{283,287},{288,292},{293,297},{298,303},
			{304,308},{309,313},{314,318},{319,323},{324,328},{329,333},
			{334,338},{339,343},{344,348},{349,353},{354,358},{359,364}
		};
		
		final int[][] lp={
			{  0,  4},{  5,  9},{ 10, 14},{ 15, 19},{ 20, 24},{ 25, 30},
			{ 31, 35},{ 36, 40},{ 41, 45},{ 46, 50},{ 51, 55},{ 56, 59},
			{ 60, 64},{ 65, 69},{ 70, 74},{ 75, 79},{ 80, 84},{ 85, 90},
			{ 91, 95},{ 96,100},{101,105},{106,110},{111,115},{116,120},
			{121,125},{126,130},{131,135},{136,140},{141,145},{146,151},
			{152,156},{157,161},{162,166},{167,171},{172,176},{177,181},
			{182,186},{187,191},{192,196},{197,201},{202,206},{207,212},
			{213,217},{218,222},{223,227},{228,232},{233,237},{238,243},
			{244,248},{249,253},{254,258},{259,263},{264,268},{269,273},
			{274,278},{279,283},{284,288},{289,293},{294,298},{299,304},
			{305,309},{310,314},{315,319},{320,324},{325,329},{330,334},
			{335,339},{340,344},{345,349},{350,354},{355,359},{360,365}
		};
		
		Range range=new Range(null,ctl);
		
		Variable buf=new Variable(cv[0].getName(),false,range);
		Variable res=new Variable("res",false,new Range(72*yearc,z,y,x));
		
		float[][][][] bdata=buf.getData();
		float[][][][] rdata=res.getData();
		
		System.out.println(" averaging...");
		
		DataRead  cdrs=DataIOFactory.getDataRead(ctl);
		DataWrite cdws=DataIOFactory.getDataWrite(ctl,path);
		
		for(int m=0;m<v;m++){
			System.out.println("    processing "+cv[m].getName());
			
			buf.setName(cv[m].getName());	buf.setUndef(cv[m].getUndef());
			res.setName(cv[m].getName());	res.setUndef(cv[m].getUndef());

			cdrs.readData(buf);
			res.setComment(buf.getComment());
			res.setUnit(buf.getUnit());
			
			for(int ll=0;ll<yearc;ll++){
				if(dcounts[ll]==365){
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						System.arraycopy(bdata[k][j][i],stags[ll],ptmp,0,365);
						
						for(int lm=0;lm<72;lm++){
							for(int lll=pp[lm][0];lll<=pp[lm][1];lll++) rdata[k][j][i][ll*72+lm]+=ptmp[lll];
							
							rdata[k][j][i][ll*72+lm]/=pp[lm][1]-pp[lm][0]+1;
						}
					}
					
				}else{
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						System.arraycopy(bdata[k][j][i],stags[ll],ltmp,0,366);
						
						for(int lm=0;lm<72;lm++){
							for(int lll=lp[lm][0];lll<=lp[lm][1];lll++) rdata[k][j][i][ll*72+lm]+=ltmp[lll];
							
							rdata[k][j][i][ll*72+lm]/=lp[lm][1]-lp[lm][0]+1;
						}
					}
				}
			}
			
			cdws.writeData(res);
		}

		cdws.writeCtl(ctl,null,"5dy");	cdrs.closeFile();	cdws.closeFile();
		
		System.out.println("Finish pentad averaging.");
	}
	
	
	/**
     * decadal mean, daily data only
     *
     * @param	path	path of file after averaging
     */
	public void decadalMean(String path){
		System.out.println("\nStart data decadal averaging...");
		
		int t=ctl.getTCount(),z=ctl.getZCount(),y=ctl.getYCount(),x=ctl.getXCount(),v=ctl.getVCount();
		
		if(t<365)
			throw new IllegalArgumentException("daily data only, length should larger than 365");
		
		int syear=ctl.getTDef().getSamples()[0].getYear();
		int yearc=t/365;
		
		Var[] cv=ctl.getVDef();
		
		int[] dcounts=new int[yearc];
		int[] stags  =new int[yearc];
		
		float[] ltmp=new float[366];
		float[] ptmp=new float[365];
		
		if(MDate.isLeapYear(syear)) dcounts[0]=366;
		else dcounts[0]=365;
		
		/*** calculate the tags ***/
		for(int i=syear+1;i<=syear+yearc-1;i++)
		if(MDate.isLeapYear(i)){
			dcounts[i-syear]=366;
			stags[i-syear]=stags[i-1-syear]+365;
			
		}else{
			dcounts[i-syear]=365;
			if(MDate.isLeapYear(i-1)) stags[i-syear]=stags[i-1-syear]+366;
			else stags[i-syear]=stags[i-1-syear]+365;
		}
		
		final int[][] pp={
			{  0,  9},{ 10, 19},{ 20, 30},
			{ 31, 40},{ 41, 50},{ 51, 58},
			{ 59, 68},{ 69, 78},{ 79, 89},
			{ 90, 99},{100,109},{110,119},
			{120,129},{130,139},{140,150},
			{151,160},{161,170},{171,180},
			{181,190},{191,200},{201,211},
			{212,221},{222,231},{232,242},
			{243,252},{253,262},{263,272},
			{273,282},{283,292},{293,303},
			{304,313},{314,323},{324,333},
			{334,343},{344,353},{354,364}
		};
		
		final int[][] lp={
			{  0,  9},{ 10, 19},{ 20, 30},
			{ 31, 40},{ 41, 50},{ 51, 59},
			{ 60, 69},{ 70, 79},{ 80, 90},
			{ 91,100},{101,110},{111,120},
			{121,130},{131,140},{141,151},
			{152,161},{162,171},{172,181},
			{182,191},{192,201},{202,212},
			{213,222},{223,232},{233,243},
			{244,253},{254,263},{264,273},
			{274,283},{284,293},{294,304},
			{305,314},{315,324},{325,334},
			{335,344},{345,354},{355,365}
		};
		
		Range range=new Range(null,ctl);
		
		Variable buf=new Variable(cv[0].getName(),false,range);
		Variable res=new Variable("res",false,new Range(36*yearc,z,y,x));
		
		float[][][][] bdata=buf.getData();
		float[][][][] rdata=res.getData();
		
		System.out.println(" averaging...");
		
		DataRead  cdrs=DataIOFactory.getDataRead(ctl);
		DataWrite cdws=DataIOFactory.getDataWrite(ctl,path);
		
		for(int m=0;m<v;m++){
			System.out.println("    processing "+cv[m].getName());
			
			buf.setName(cv[m].getName());	buf.setUndef(cv[m].getUndef());
			res.setName(cv[m].getName());	res.setUndef(cv[m].getUndef());

			cdrs.readData(buf);
			res.setComment(buf.getComment());
			res.setUnit(buf.getUnit());
			
			for(int ll=0;ll<yearc;ll++){
				if(dcounts[ll]==365){
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						System.arraycopy(bdata[k][j][i],stags[ll],ptmp,0,365);
						
						for(int lm=0;lm<36;lm++){
							for(int lll=pp[lm][0];lll<=pp[lm][1];lll++) rdata[k][j][i][ll*36+lm]+=ptmp[lll];
							
							rdata[k][j][i][ll*36+lm]/=pp[lm][1]-pp[lm][0]+1;
						}
					}
					
				}else{
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						System.arraycopy(bdata[k][j][i],stags[ll],ltmp,0,366);
						
						for(int lm=0;lm<36;lm++){
							for(int lll=lp[lm][0];lll<=lp[lm][1];lll++) rdata[k][j][i][ll*36+lm]+=ltmp[lll];
							
							rdata[k][j][i][ll*36+lm]/=lp[lm][1]-lp[lm][0]+1;
						}
					}
				}
			}
			
			cdws.writeData(res);
		}

		cdws.writeCtl(ctl,null,"10dy");	cdrs.closeFile();	cdws.closeFile();
		
		System.out.println("Finish decadal averaging.");
	}
	
	
	/**
     * monthly mean, daily data only
     *
     * @param	path	path of file after averaging
     */
	public void monthlyMean(String path){
		System.out.println("\nStart data monthly averaging...");
		
		int t=ctl.getTCount(),z=ctl.getZCount(),y=ctl.getYCount(),x=ctl.getXCount(),v=ctl.getVCount();
		
		if(t<365)
			throw new IllegalArgumentException("daily data only, length should larger than 365");
		
		int syear=ctl.getTDef().getSamples()[0].getYear();
		int yearc=t/365;
		
		Var[] cv=ctl.getVDef();
		
		int[] dcounts=new int[yearc];
		int[] stags  =new int[yearc];
		
		float[] ltmp=new float[366];
		float[] ptmp=new float[365];
		
		if(MDate.isLeapYear(syear)) dcounts[0]=366;
		else dcounts[0]=365;
		
		/*** calculate the tags ***/
		for(int i=syear+1;i<=syear+yearc-1;i++)
		if(MDate.isLeapYear(i)){
			dcounts[i-syear]=366;
			stags[i-syear]=stags[i-1-syear]+365;
			
		}else{
			dcounts[i-syear]=365;
			if(MDate.isLeapYear(i-1)) stags[i-syear]=stags[i-1-syear]+366;
			else stags[i-syear]=stags[i-1-syear]+365;
		}
		
		final int[][] pp={
			{  0, 30},{ 31, 58},{ 59, 89},{ 90,119},{120,150},{151,180},
			{181,211},{212,242},{243,272},{273,303},{304,333},{334,364}
		};
		
		final int[][] lp={
			{  0, 30},{ 31, 59},{ 60, 90},{ 91,120},{121,151},{152,181},
			{182,212},{213,243},{244,273},{274,304},{305,334},{335,365}
		};
		
		Range range=new Range(null,ctl);
		
		Variable buf=new Variable(cv[0].getName(),false,range);
		Variable res=new Variable("res",false,new Range(12*yearc,z,y,x));
		
		float[][][][] bdata=buf.getData();
		float[][][][] rdata=res.getData();
		
		System.out.println(" averaging...");
		
		DataRead  cdrs=DataIOFactory.getDataRead(ctl);
		DataWrite cdws=DataIOFactory.getDataWrite(ctl,path);
		
		for(int m=0;m<v;m++){
			System.out.println("    processing "+cv[m].getName());
			
			float undef=ctl.getUndef(cv[m].getName());
			
			buf.setName(cv[m].getName());	buf.setUndef(cv[m].getUndef());
			res.setName(cv[m].getName());	res.setUndef(cv[m].getUndef());
			
			cdrs.readData(buf);
			res.setComment(buf.getComment());
			res.setUnit(buf.getUnit());
			
			for(int ll=0;ll<yearc;ll++){
				if(dcounts[ll]==365){
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						System.arraycopy(bdata[k][j][i],stags[ll],ptmp,0,365);
						
						for(int lm=0;lm<12;lm++){
							int count=0;
							
							for(int lll=pp[lm][0];lll<=pp[lm][1];lll++)
							if(ptmp[lll]!=undef){ rdata[k][j][i][ll*12+lm]+=ptmp[lll]; count++;}
							
							rdata[k][j][i][ll*12+lm]/=count;
						}
					}
					
				}else{
					for(int k=0;k<z;k++)
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++){
						System.arraycopy(bdata[k][j][i],stags[ll],ltmp,0,366);
						
						for(int lm=0;lm<12;lm++){
							int count=0;
							
							for(int lll=lp[lm][0];lll<=lp[lm][1];lll++)
							if(ltmp[lll]!=undef){ rdata[k][j][i][ll*12+lm]+=ltmp[lll]; count++;}
							
							rdata[k][j][i][ll*12+lm]/=count;
						}
					}
				}
			}
			
			cdws.writeData(res);
		}
		
		cdws.writeCtl(ctl,null,"1mo");	cdrs.closeFile();	cdws.closeFile();
		
		System.out.println("Finish monthly averaging.");
	}
	
	
	/** test
	public static void main(String[] args){
		String var="h";
		
		DataMean dm=new DataMean(new CtlDescriptor("d:\\200\\"+var+".ctl"));
		
		dm.monthlyMean("d:\\200\\monthly\\"+var+"mly.dat");
	}*/
}
