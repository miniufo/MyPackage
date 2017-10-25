/**
 * @(#)CtsDescriptor.java	1.0 2017.06.01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import java.io.File;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Scanner;
import miniufo.diagnosis.MDate;
import static java.nio.ByteOrder.BIG_ENDIAN;
import static java.nio.ByteOrder.LITTLE_ENDIAN;


/**
 * Used to describe the control file (similar to ctl file for GrADS)
 * in Cartesian coordinates.
 *
 * @version 1.0, 2017.06.01
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CtsDescriptor extends CtlDescriptor{
	//
	private float  f0 =Float.NaN;
	private float beta=Float.NaN;
	
	
	/**
     * constructor
     *
     * param	content		content of the ctl file
     */
	public CtsDescriptor(String content){ super(); parse(content);}
	
	/**
     * constructor
     *
     * param	ctlFile		a ctl file
     */
	public CtsDescriptor(File ctsFile){
		super();
		
		descPath=ctsFile.getAbsolutePath().replace("\\","/");
		
		StringBuilder content=new StringBuilder();
		
		try{ Files.lines(ctsFile.toPath()).forEach(s->content.append(s+"\n"));}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		int idx=content.indexOf("^");
		if(idx!=-1) content.replace(idx,idx+1,ctsFile.getParentFile().getAbsolutePath().replace("\\","/")+"/");
		
		parse(content.toString());
	}
	
	
	/*** getor and setor ***/
	public float getF0(){ return f0;}
	
	public float getBeta(){ return beta;}
	
	
	/*** helper methods ***/
	private void parse(String content){
		float undef=Float.NaN;
		
		Scanner scnr_file=new Scanner(content);
		
		while(scnr_file.hasNextLine()){
			String oneline=scnr_file.nextLine();
			if(oneline.startsWith("*")) continue;	// skip comments
			
			if(!"".equals(oneline)){
				Scanner scnr_line=new Scanner(oneline);
				
				if(!scnr_line.hasNext()) break;
				
				String start_word=scnr_line.next().toLowerCase();
				
				     if(start_word.equals("title"  )) title=scnr_line.next();
				else if(start_word.equals("undef"  )) undef=Float.parseFloat(scnr_line.next());
				else if(start_word.equals("dset"   )) processDSet(scnr_line);
				else if(start_word.equals("options")) processOptions(scnr_line);
				else if(start_word.equals("xdef"   )) processXDef(scnr_line,scnr_file);
				else if(start_word.equals("ydef"   )) processYDef(scnr_line,scnr_file);
				else if(start_word.equals("zdef"   )) processZDef(scnr_line,scnr_file);
				else if(start_word.equals("tdef"   )) processTDef(scnr_line);
				else if(start_word.equals("f0"     )) f0  =Float.parseFloat(scnr_line.next());
				else if(start_word.equals("beta"   )) beta=Float.parseFloat(scnr_line.next());
				else if(start_word.equals("vars"   )) processVars(scnr_line,scnr_file,undef);
				
				scnr_line.close();
			}
		}
		
		scnr_file.close();
		
		/*** check information ***/
		if(    tdef==null) throw new IllegalArgumentException("missing tdef information");
		if(    zdef==null) throw new IllegalArgumentException("missing zdef information");
		if(    vdef==null) throw new IllegalArgumentException("missing vdef information");
		if(dsetPath==null) throw new IllegalArgumentException("missing data set path"   );
		
		if(Float.isNaN( f0 )) throw new IllegalArgumentException("missing f0"  );
		if(Float.isNaN(beta)) throw new IllegalArgumentException("missing beta");
		
		postProcess();
	}
	
	private void processDSet(Scanner sc){
		dsetPath=sc.next();
		hasData=Files.exists(Paths.get(dsetPath));
	}
	
	private void processOptions(Scanner sc){
		String tmp_options=sc.nextLine().toLowerCase();
		
		if(tmp_options.indexOf("yrev"            )!=-1) yrev=true;
		if(tmp_options.indexOf("zrev"            )!=-1) zrev=true;
		if(tmp_options.indexOf("template"        )!=-1) template=true;
		if(tmp_options.indexOf("sequential"      )!=-1) sequential=true;
		if(tmp_options.indexOf("big_endian"      )!=-1) order=ByteOrder.BIG_ENDIAN;
		if(tmp_options.indexOf("byteswapped"     )!=-1) order=(order==BIG_ENDIAN?LITTLE_ENDIAN:BIG_ENDIAN);
		if(tmp_options.indexOf("365_day_calendar")!=-1) cal_365_day=true;
	}
	
	private void processXDef(Scanner scL,Scanner scF){
		int xcount=Integer.parseInt(scL.next());
		
		if(xcount==0) throw new IllegalArgumentException("x count is 0");
		
		float[] df=new float[xcount];	String mapping=scL.next().toLowerCase();
		
		float incre=-1.0f;
		
		if(mapping.equals("linear")){
			df[0]=Float.parseFloat(scL.next());
			incre=Float.parseFloat(scL.next());
			
			if(incre<=0) throw new IllegalArgumentException("Illegal (negative) xdef increment");
			
			for(int i=1;i<xcount;i++) df[i]=df[0]+i*incre;
			
			if(xcount==1){ dxdef=new float[1]; dxdef[0]=incre;}
			
			xdef=new SpatialCoordinate("xdef",df);
			
		}else if(mapping.equals("levels")){
			if(scL.hasNext()) for(int i=0;i<xcount;i++) df[i]=Float.parseFloat(scL.next()    );
			else              for(int i=0;i<xcount;i++) df[i]=Float.parseFloat(scF.nextLine());
			
			xdef=new SpatialCoordinate("xdef",df);
			
		}else throw new IllegalArgumentException("unsupported xdef mapping");
	}
	
	private void processYDef(Scanner scL,Scanner scF){
		int ycount=Integer.parseInt(scL.next());
		
		if(ycount==0) throw new IllegalArgumentException("y count is 0");
		
		float[] df=new float[ycount];	String mapping=scL.next().toLowerCase();
		
		float incre=1.0f;
		
		if(mapping.equals("linear")){
			df[0]=Float.parseFloat(scL.next());
			incre=Float.parseFloat(scL.next());
			
			if(incre<=0) throw new IllegalArgumentException("Illegal (negative) ydef increment");
			
			for(int j=1;j<ycount;j++) df[j]=df[0]+j*incre;
			
			if(ycount==1){ dydef=new float[1]; dydef[0]=incre;}
			
			ydef=new SpatialCoordinate("ydef",df);
			
		}else if(mapping.equals("levels")){
			if(scL.hasNext()) for(int j=0;j<ycount;j++) df[j]=Float.parseFloat(scL.next()    );
			else              for(int j=0;j<ycount;j++) df[j]=Float.parseFloat(scF.nextLine());
			
			ydef=new SpatialCoordinate("ydef",df);
			
		}else throw new IllegalArgumentException("unsupported ydef mapping");
	}
	
	private void processZDef(Scanner scL,Scanner scF){
		int zcount=Integer.parseInt(scL.next());
		
		if(zcount==0) throw new IllegalArgumentException("z count is 0");
		
		float[] df=new float[zcount];	String mapping=scL.next().toLowerCase();
		
		float incre=50;
		
		if(mapping.equals("linear")){
			df[0]=Float.parseFloat(scL.next());
			incre=Float.parseFloat(scL.next());
			
			if(incre<=0) throw new IllegalArgumentException("Illegal (negative) zdef increment");
			
			for(int i=1;i<zcount;i++) df[i]=df[0]+i*incre;
			
			if(zcount==1){ dzdef=new float[1]; dzdef[0]=incre;}
			
			zdef=new SpatialCoordinate("zdef",df);
			
		}else if(mapping.equals("levels")){
			if(scL.hasNext()) for(int i=0;i<zcount;i++) df[i]=Float.parseFloat(scL.next()    );
			else              for(int i=0;i<zcount;i++) df[i]=Float.parseFloat(scF.nextLine());
			
			if(zcount==1){ dzdef=new float[1]; dzdef[0]=incre;}
			
			zdef=new SpatialCoordinate("zdef",df);
			
		}else throw new IllegalArgumentException("unsupported zdef mapping");
	}
	
	private void processTDef(Scanner sc){
		int tcount=Integer.parseInt(sc.next());
		
		if(tcount==0) throw new IllegalArgumentException("t count is 0");
		
		MDate[] df=new MDate[tcount];
		
		String mapping=sc.next().toLowerCase();
		
		if(mapping.equals("linear")){
			df[0]=new MDate(sc.next());
			tincrement=sc.next();
			
			for(int i=1;i<tcount;i++) df[i]=df[i-1].add(tincrement);
			
			tdef=new TemporalCoordinate("tdef",df,true);
			
		}else throw new IllegalArgumentException("unsupported tdef mapping");
		
		dtdef=new float[1];
		
		tincrement=tincrement.toLowerCase();
		
		if(tcount>1) dtdef[0]=df[1].getDT(df[0]);
		else switch(IncreType.valueOf(tincrement.substring(tincrement.length()-2))){
			case mn: dtdef[0]=Integer.parseInt(tincrement.replace("mn",""))*60;           break;
			case hr: dtdef[0]=Integer.parseInt(tincrement.replace("hr",""))*60*60;        break;
			case dy: dtdef[0]=Integer.parseInt(tincrement.replace("dy",""))*60*60*24;     break;
			case mo: dtdef[0]=Integer.parseInt(tincrement.replace("mo",""))*60*60*24*30;  break;
			case yr: dtdef[0]=Integer.parseInt(tincrement.replace("yr",""))*60*60*24*365; break;
			default: throw new IllegalArgumentException("unsupported time increment type");
		}
	}
	
	private void processVars(Scanner scL,Scanner scF,float undef){
		if(tdef==null||zdef==null||ydef==null||xdef==null)
			throw new NullPointerException("vdef should be after x,y,z,t definition in ctl file");
		
		if(Float.isNaN(undef))
			throw new IllegalArgumentException("invalid undef/no undef value before vdef");
		
		int tcount=tdef.length();
		int ycount=ydef.length();
		int xcount=xdef.length();
		
		one_level_length=xcount*ycount<<2;
		
		vcount=Integer.parseInt(scL.next());
		vdef=new CtlVar[vcount];
		
		/*** first var ***/
		String vline=scF.nextLine();
		while(vline.equals("")||vline==null) vline=scF.nextLine();	// skip space line
		
		vdef[0]=new CtlVar(vline);
		
		vdef[0].setTCount(tcount);	vdef[0].setYCount(ycount);
		vdef[0].setXCount(xcount);	vdef[0].setUndef(undef);
		
		CtlVar tmp=(CtlVar)vdef[0];
		
		tmp.setIndex(0);	tmp.setStartPosition(0);
		
		totalZcount+=vdef[0].getZCount();
		
		storageType=tmp.getStorageType();
		
		boolean t1=false,t2=false;
		
		/*** vars remain ***/
		for(int i=1;i<vcount;i++){
			vline=scF.nextLine();
			
			while(vline.equals("")||vline==null) vline=scF.nextLine();	// skip space line
			
			vdef[i]=new CtlVar(vline);	tmp=(CtlVar)vdef[i];
			vdef[i].setTCount(tcount);	vdef[i].setYCount(ycount);
			vdef[i].setXCount(xcount);	vdef[i].setUndef(undef);
			
			tmp.setIndex(i);
			
			if(tmp.getStorageType().equals("99")){
				tmp.setStartPosition(vdef[i-1].getZCount()*one_level_length+((CtlVar)(vdef[i-1])).getStartPosition());
				t1=true;
			}else{
				tmp.setStartPosition(vdef[i-1].getZCount()*one_level_length*tcount+((CtlVar)(vdef[i-1])).getStartPosition());
				t2=true;
			}
			
			if(t1==t2) throw new UnsupportedClassVersionError("storage type should be same");
			
			totalZcount+=vdef[i].getZCount();
		}
		
		t_rec_length=one_level_length*totalZcount;
	}
	
	
	/** test
	public static void main(String arg[]){
		CtsDescriptor cts=new CtsDescriptor(new File("D:/Data/MITgcm/barotropicDG/BetaCartRL/Stat.cts"));
		System.out.println(cts.getF0());
		System.out.println(cts.getBeta());
		System.out.println(cts.getVarCommentAndUnit("tr1"));
	}*/
}
