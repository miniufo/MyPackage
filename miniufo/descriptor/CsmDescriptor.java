/**
 * @(#)CsmDescriptor.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Scanner;
import miniufo.diagnosis.MDate;
import static java.lang.Math.toRadians;
import static miniufo.diagnosis.SpatialModel.cLatLons;


/**
 * A class used to describe the cylindrical spatial model.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CsmDescriptor extends DataDescriptor{
	//
	private float[]   olon =null;	// original x-coordinate (degree)
	private float[]   olat =null;	// original y-coordinate (degree)
	private float[]   prs  =null;	// minimum pressure
	private float[]   wnd  =null;	// maximum wind speed
	
	private float[][][] lon=null;	// units: degree, lon[t][y][x]
	private float[][][] lat=null;	// units: degree, lat[t][y][x]
	private float[][][] eta=null;	// units: degree, eta[t][y][x], angle between radial direction and north (counter-clockwise)
	
	private boolean hasCtl=false;
	
	private CtlDescriptor ctl=null;
	
	private enum IncreType{mn,hr,dy,mo,yr};
	
	
	/**
     * constructor
     *
     * param	content		content of the csm file
     */
	public CsmDescriptor(String content){ parse(content);}
	
	/**
     * constructor
     *
     * param	csmFile		a csm file
     */
	public CsmDescriptor(File csmFile){
		descPath=csmFile.getAbsolutePath().replace("\\","/");
		
		StringBuilder content=new StringBuilder();
		
		try{ Files.lines(csmFile.toPath()).forEach(s->content.append(s+"\n"));}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		int idx=content.indexOf("^");
		if(idx!=-1)
		content.replace(idx,idx+1,csmFile.getParentFile().getAbsolutePath().replace("\\","/")+"/");
		
		parse(content.toString());
	}
	
	
	/*** getor and setor ***/
	public float getUndef(String vname){
		if(hasCtl) return ctl.getUndef(null);
		else return Float.NaN;
	}
	
	public float[] getMinPres(){ return prs;}
	
	public float[] getMaxWind(){ return wnd;}
	
	public float[] getOLon(){ return olon;}
	
	public float[] getOLat(){ return olat;}
	
	public float[][][] getLon(){ return lon;}
	
	public float[][][] getLat(){ return lat;}
	
	public float[][][] getEta(){ return eta;}
	
	public boolean hasCtl(){ return hasCtl;}
	
	public CtlDescriptor getCtlDescriptor(){ return ctl;}
	
	public void setCtlDescriptor(CtlDescriptor ctl){
		hasCtl=true;
		
		this.ctl=ctl;
		
		hasData=this.ctl.hasData;
		
		/*** set descriptor status ***/
		dsetPath=this.ctl.getDSet();	yrev=this.ctl.yrev;
		vdef    =this.ctl.getVDef();	zrev=this.ctl.zrev;	vcount=this.ctl.getVCount();
	}
	
	
	/**
     * parse content in terms of a string
     */
	private void parse(String content){
		Scanner scnr_content=new Scanner(content);
		
		while(scnr_content.hasNextLine()){
			String oneline=scnr_content.nextLine();
			if(oneline.startsWith("*")) continue;	// skip comments
			
			if(!"".equals(oneline)){
				Scanner scnr_line=new Scanner(oneline);
				
				if(!scnr_line.hasNext()) break;
				
				String start_word=scnr_line.next().toLowerCase();
				
				if(start_word.equals("title")) title=scnr_line.next();
				else if(start_word.equals("dpath" )) processDPath (scnr_line);
				else if(start_word.equals("xdef"  )) processXDef  (scnr_line);
				else if(start_word.equals("ydef"  )) processYDef  (scnr_line);
				else if(start_word.equals("zdef"  )) processZDef  (scnr_line);
				else if(start_word.equals("tdef"  )) processTDef  (scnr_line);
				else if(start_word.equals("coords")) processCoords(scnr_content);
				
				scnr_line.close();
			}
		}
		
		scnr_content.close();
		
		/*** check information ***/
		if(    tdef==null) throw new IllegalArgumentException("missing tdef information");
		if(    ydef==null) throw new IllegalArgumentException("missing ydef information");
		if(    xdef==null) throw new IllegalArgumentException("missing xdef information");
		if(dsetPath==null) throw new IllegalArgumentException("missing descriptor path" );
		
		/*** calculate lon, lat and eta ***/
		int tcount=tdef.length();	eta=new float[tcount][][];
		int ycount=ydef.length();	lon=new float[tcount][][];
		int xcount=xdef.length();	lat=new float[tcount][][];
		
		for(int l=0;l<tcount;l++){
			float[][][] re=cLatLons(
				toRadians(olon[l]),toRadians(olat[l]),
				toRadians(ydef.getIncrements()[0]),ycount,xcount
			);
			
			lon[l]=re[0];	lat[l]=re[1];	eta[l]=re[2];
			
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++){
				lon[l][j][i]=(float)(Math.toDegrees(lon[l][j][i]));
				lat[l][j][i]=(float)(Math.toDegrees(lat[l][j][i]));
				eta[l][j][i]=(float)(Math.toDegrees(eta[l][j][i]));
			}
		}
		
		periodicX=true;
	}
	
	private void processDPath(Scanner sc){
		dsetPath=sc.next();
		
		hasCtl=Files.exists(Paths.get(dsetPath));
		
		if(hasCtl){
			if(ctl==null) ctl=new CtlDescriptor(new File(dsetPath));
			
			hasData=ctl.hasData;
			
			/*** set descriptor status ***/
			yrev=ctl.yrev;	vdef  =ctl.getVDef();
			zrev=ctl.zrev;	vcount=ctl.getVCount();
			
		}else if(ctl==null){ hasCtl=false; hasData=false;}
	}
	
	private void processXDef(Scanner sc){
		int xcount=Integer.parseInt(sc.next());
		
		if(xcount<=0) throw new IllegalArgumentException("x count should be positive");
		
		dxdef=new float[]{360.0f/xcount};
		
		float[] df=new float[xcount];
		for(int i=0;i<xcount;i++) df[i]=i*dxdef[0];
		
		xdef=new SpatialCoordinate("xdef",df);
	}
	
	private void processYDef(Scanner sc){
		int ycount=Integer.parseInt(sc.next());
		
		if(ycount<=0) throw new IllegalArgumentException("y count should be positive");
		
		dydef=new float[]{Float.parseFloat(sc.next())};
		
		if(dydef[0]<=0) throw new IllegalArgumentException("Illegal (negative) ydef increment");
		
		float[] df=new float[ycount];
		for(int j=0;j<ycount;j++) df[j]=j*dydef[0];
		
		ydef=new SpatialCoordinate("ydef",df);
	}
	
	private void processZDef(Scanner sc){
		int zcount=Integer.parseInt(sc.next());
		
		if(zcount<=0) throw new IllegalArgumentException("z count should be positive");
		
		float[] df=new float[zcount];
		
		String next=sc.next();
		
		if(next.equalsIgnoreCase("levels")){
			dzdef=new float[zcount-1];
			for(int k=0;k<zcount;k++) df[k]=Float.parseFloat(sc.next());
			for(int k=1;k<zcount;k++) dzdef[k-1]=df[k]-df[k-1];
			
		}else{
			df[0]=Float.parseFloat(next);
			dzdef=new float[]{Float.parseFloat(sc.next())};
			
			for(int i=1;i<zcount;i++) df[i]=df[i-1]+dzdef[0];
			//for(int i=1;i<zcount-1;i++) dzdef[i]=dzdef[0];
		}
		
		zdef=new SpatialCoordinate("zdef",df);
	}
	
	private void processTDef(Scanner sc){
		int tcount=Integer.parseInt(sc.next());
		
		if(tcount==0) throw new IllegalArgumentException("t count is 0");
		
		dtdef=new float[1];
		times=new long[tcount+1];
		MDate[] df=new MDate[tcount];
		
		df[0]=new MDate(sc.next());
		tincrement=sc.next().toLowerCase();
		
		for(int i=1;i<tcount;i++) df[i]=df[i-1].add(tincrement);
		
		if(tcount>1) dtdef[0]=df[1].getDT(df[0]);
		else{
			switch(IncreType.valueOf(tincrement.substring(tincrement.length()-2))){
			case mn: dtdef[0]=Integer.parseInt(tincrement.replace("mn",""))*60;           break;
			case hr: dtdef[0]=Integer.parseInt(tincrement.replace("hr",""))*60*60;        break;
			case dy: dtdef[0]=Integer.parseInt(tincrement.replace("dy",""))*60*60*24;     break;
			case mo: dtdef[0]=Integer.parseInt(tincrement.replace("mo",""))*60*60*24*30;  break;
			case yr: dtdef[0]=Integer.parseInt(tincrement.replace("yr",""))*60*60*24*365; break;
			default: throw new IllegalArgumentException("unsupported time increment type");
			}
		}
		
		for(int i=0;i<tcount;i++) times[i]=df[i].getLongTime();
		times[tcount]=df[tcount-1].add(tincrement).getLongTime();
		
		tdef=new TemporalCoordinate("tdef",df,true);
	}
	
	private void processCoords(Scanner sc){
		if(tdef==null) throw new NullPointerException("tdef do not define before coords");
		
		olon=new float[tdef.length()]; wnd=new float[tdef.length()];
		olat=new float[tdef.length()]; prs=new float[tdef.length()];
		
		boolean hasWnd=false,hasPrs=false;
		
		for(int l=0,L=tdef.length();l<L;l++){
			Scanner scnr_tmp=new Scanner(sc.nextLine());
			
			olon[l]=Float.parseFloat(scnr_tmp.next());
			olat[l]=Float.parseFloat(scnr_tmp.next());
			
			if(scnr_tmp.hasNext()){ prs[l]=Float.parseFloat(scnr_tmp.next()); hasPrs=true;}
			if(scnr_tmp.hasNext()){ wnd[l]=Float.parseFloat(scnr_tmp.next()); hasWnd=true;}
			
			scnr_tmp.close();
		}
		
		if(!"endcoords".equalsIgnoreCase(sc.nextLine()))
			throw new IllegalArgumentException("no endcoords or invlid count of coords");
		
		if(hasWnd) wnd=null;
		if(hasPrs) prs=null;
	}
	
	
	/** test
	public static void main(String arg[]){
		CsmDescriptor cd=new CsmDescriptor(new File("D:/Data/DiagnosisVortex/Haima/Haima.csm"));
		
		System.out.println(cd.descPath);
		System.out.println(cd.dsetPath);
		
		System.out.println(cd);
	}*/
}
