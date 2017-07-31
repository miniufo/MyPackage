/**
 * @(#)DataInterpolation.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import miniufo.io.DataWrite;
import miniufo.io.FileWriteInterface;
import miniufo.io.IOUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.Var;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.mathsphysics.Spline;
import static miniufo.io.FileWriteInterface.Solution.*;
import static miniufo.basic.InterpolationModel.Type;
import static miniufo.basic.InterpolationModel.interp1D;
import static miniufo.basic.InterpolationModel.interp2D;
import static miniufo.basic.InterpolationModel.bilinearInterpolation;
import static miniufo.basic.InterpolationModel.bicubicPolynomialInterpolation;
import static miniufo.basic.InterpolationModel.bicubicLagrangeInterpolation;


/**
 * interpolation of a data file into a new one
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DataInterpolation{
	//
	private DataDescriptor srcdata=null;
	
	
	/**
     * constructor
     *
     * @param	src		source DataDescriptor
     */
	public DataInterpolation(DataDescriptor src){ srcdata=src;}
	
	
	/**
     * horizontal interpolation
     *
     * @param	path	path of file after interpolation
     * @param	type	type of interpolation, e.g. cubic, linear
     * @param	y		y count
     * @param	x		x count
     */
	public void horizontalInterp(String path,Type zonalType,Type meridionalType,int y,int x){
		System.out.println("Start horizontal interpolating...");
		
		Range srange=new Range(1,1,srcdata.getYCount(),srcdata.getXCount());
		Range drange=new Range(1,1,y,x);
		
		Variable srcv=new Variable(null,srange);
		Variable desv=new Variable(null,drange);
		
		DataRead  cdrs=DataIOFactory.getDataRead(srcdata);
		DataWrite cdws=DataIOFactory.getDataWrite(srcdata,path);
		
		float[][] srcvdata=srcv.getData()[0][0];
		float[][] desvdata=desv.getData()[0][0];
		
		for(int l=1;l<=srcdata.getTCount();l++){
			System.out.println("  Processing "+srcdata.getTDef().getSamples()[l-1].toString()+" ...");
			
			int[] p=srcv.getRange().getTRange();	p[0]=l;p[1]=l;
			
			for(int m=0;m<srcdata.getVCount();m++){
				srcv.setName(srcdata.getVDef()[m].getName());
				desv.setName(srcv.getName());
				
				float undef=srcdata.getUndef(srcdata.getVDef()[m].getName());
				
				for(int k=1;k<=srcdata.getVDef()[m].getZCount();k++){
					p=srcv.getRange().getZRange();	p[0]=k;p[1]=k;
					
					cdrs.readData(srcv);
					interp2D(srcvdata,desvdata,zonalType,meridionalType,undef);
					cdws.writeData(desv);
				}
			}
		}
		
		cdrs.closeFile();	cdws.closeFile();
		
		// write ctl
		FileWriteInterface fwi=new FileWriteInterface(IOUtil.getCompleteFileNameWithoutExtension(path)+".ctl");
		
		if(fwi.getFlag()!=SKIP){
			FileWriter fw=null;
			Scanner sn=null;
			
			try{
				switch(fwi.getFlag()){
				case RENAME   : fw=new FileWriter(fwi.getParent()+fwi.getNewName()); break;
				case OVERWRITE: fw=new FileWriter(fwi.getFile());                    break;
				case APPEND   : fw=new FileWriter(fwi.getFile(),true);               break;
				default       : throw new IllegalArgumentException("unsupported flag for "+fwi.getFlag());
				}
				
				sn=new Scanner(new File(srcdata.getPath()));
			
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			
			StringBuilder sb=new StringBuilder();
			
			while(sn.hasNextLine()){
				String line=sn.nextLine();
				
				if(line.startsWith("options")) continue;
				else if(line.startsWith("xdef")){
					sb.append("xdef ");	sb.append(x);	sb.append(" linear ");
					
					Scanner tmp=new Scanner(line);
					
					tmp.next();	int xc=Integer.parseInt(tmp.next());
					tmp.next();	sb.append(tmp.next());
					
					sb.append(" ");
					
					if(zonalType==Type.PERIODIC_CUBIC_L||zonalType==Type.PERIODIC_CUBIC_P||zonalType==Type.PERIODIC_SPLINE)
						sb.append(Float.valueOf(tmp.next())*xc/x);
					else
						sb.append(Float.valueOf(tmp.next())*(xc-1)/(x-1));
					
					sb.append("\n");	tmp.close();
					
				}else if(line.startsWith("ydef")){
					sb.append("ydef ");	sb.append(y);	sb.append(" linear ");
					
					Scanner tmp=new Scanner(line);
					
					tmp.next();	int yc=Integer.parseInt(tmp.next());
					tmp.next();	sb.append(tmp.next());
					
					sb.append(" "); sb.append(Float.parseFloat(tmp.next())*(yc-1)/(y-1));
					
					sb.append("\n");	tmp.close();
					
				}else if(line.startsWith("dset")){
					sb.append("dset ^"+IOUtil.getFileName(path));
					sb.append("\n");
					
				}else{
					sb.append(line);
					sb.append("\n");
				}
			}
			
			sn.close();
			
			try{ fw.write(sb.toString());	fw.close();}
			catch(IOException e){  e.printStackTrace(); System.exit(0);}
		}
		
		System.out.println("Finish horizontal interpolating.");
	}
	
	public void horizontalInterp(String path,Type type,float[] coordsY,float[] coordsX){
		System.out.println("Start horizontal interpolating...");
		
		int sy=srcdata.getYCount(),sx=srcdata.getXCount();
		int y=coordsY.length,x=coordsX.length;
		
		Range srange=new Range(1,1,sy,sx);
		Range drange=new Range(1,1,y,x);
		
		Variable srcv=new Variable(null,srange);
		Variable desv=new Variable(null,drange);
		
		DataRead  cdrs=DataIOFactory.getDataRead(srcdata);
		DataWrite cdws=DataIOFactory.getDataWrite(srcdata,path);
		
		float[][] srcvdata=srcv.getData()[0][0];
		float[][] desvdata=desv.getData()[0][0];
		
		float[] srcxdef=srcdata.getXDef().getSamples();
		float[] srcydef=srcdata.getYDef().getSamples();
		
		switch(type){
		case LINEAR:
			for(int l=1;l<=srcdata.getTCount();l++){
				System.out.println("  Processing "+srcdata.getTDef().getSamples()[l-1].toString()+" ...");
				
				int[] p=srcv.getRange().getTRange();	p[0]=l;p[1]=l;
				
				for(int m=0;m<srcdata.getVCount();m++){
					srcv.setName(srcdata.getVDef()[m].getName());
					desv.setName(srcv.getName());
					
					for(int k=1;k<=srcdata.getVDef()[m].getZCount();k++){
						p=srcv.getRange().getZRange();	p[0]=k;p[1]=k;
						
						cdrs.readData(srcv); float undef=srcv.getUndef(); desv.setValue(undef);
						
						for(int j=0;j<y;j++){
							int tagj=srcdata.getYLENum(coordsY[j]);
							
							if(tagj>=0) for(int i=0;i<x;i++){
								int tagi=srcdata.getXLENum(coordsX[i]);
								
								if(tagi>=0) desvdata[j][i]=bilinearInterpolation(
									srcvdata[tagj  ][tagi],srcvdata[tagj  ][tagi+1],
									srcvdata[tagj+1][tagi],srcvdata[tagj+1][tagi+1],
									(coordsX[i]-srcxdef[tagi])/(srcxdef[tagi+1]-srcxdef[tagi]),
									(coordsY[j]-srcydef[tagj])/(srcydef[tagj+1]-srcydef[tagj]),undef
								);
							}
						}
						
						cdws.writeData(desv);
					}
				}
			}
			
			break;
			
		case CUBIC_P:
			for(int l=1;l<=srcdata.getTCount();l++){
				System.out.println("  Processing "+srcdata.getTDef().getSamples()[l-1].toString()+" ...");
				
				int[] p=srcv.getRange().getTRange();	p[0]=l;p[1]=l;
				
				for(int m=0;m<srcdata.getVCount();m++){
					srcv.setName(srcdata.getVDef()[m].getName());
					desv.setName(srcv.getName());
					
					for(int k=1;k<=srcdata.getVDef()[m].getZCount();k++){
						p=srcv.getRange().getZRange();	p[0]=k;p[1]=k;
						
						cdrs.readData(srcv); float undef=srcv.getUndef();
						
						for(int j=0;j<y;j++){
							int tagj=srcdata.getYLENum(coordsY[j]);
							
							if(tagj>=0) for(int i=0;i<x;i++){
								int tagi=srcdata.getXLENum(coordsX[i]);
								
								if(tagi>=0) desvdata[j][i]=bicubicPolynomialInterpolation(
									srcvdata[tagj-1][tagi-1],srcvdata[tagj-1][tagi],srcvdata[tagj-1][tagi+1],srcvdata[tagj-1][tagi+2],
									srcvdata[tagj  ][tagi-1],srcvdata[tagj  ][tagi],srcvdata[tagj  ][tagi+1],srcvdata[tagj  ][tagi+2],
									srcvdata[tagj+1][tagi-1],srcvdata[tagj+1][tagi],srcvdata[tagj+1][tagi+1],srcvdata[tagj+1][tagi+2],
									srcvdata[tagj+2][tagi-1],srcvdata[tagj+2][tagi],srcvdata[tagj+2][tagi+1],srcvdata[tagj+2][tagi+2],
									(coordsX[i]-srcxdef[tagi])/(srcxdef[tagi+1]-srcxdef[tagi]),
									(coordsY[j]-srcydef[tagj])/(srcydef[tagj+1]-srcydef[tagj]),undef
								);
							}
						}
						
						cdws.writeData(desv);
					}
				}
			}
			
			break;
			
		case CUBIC_L:
			for(int l=1;l<=srcdata.getTCount();l++){
				System.out.println("  Processing "+srcdata.getTDef().getSamples()[l-1].toString()+" ...");
				
				int[] p=srcv.getRange().getTRange();	p[0]=l;p[1]=l;
				
				for(int m=0;m<srcdata.getVCount();m++){
					srcv.setName(srcdata.getVDef()[m].getName());
					desv.setName(srcv.getName());
					
					for(int k=1;k<=srcdata.getVDef()[m].getZCount();k++){
						p=srcv.getRange().getZRange();	p[0]=k;p[1]=k;
						
						cdrs.readData(srcv); float undef=srcv.getUndef();
						
						for(int j=0;j<y;j++){
							int tagj=srcdata.getYLENum(coordsY[j]);
							
							if(tagj>=0) for(int i=0;i<x;i++){
								int tagi=srcdata.getXLENum(coordsX[i]);
								
								if(tagi>=0) desvdata[j][i]=bicubicLagrangeInterpolation(
									srcvdata[tagj-1][tagi-1],srcvdata[tagj-1][tagi],srcvdata[tagj-1][tagi+1],srcvdata[tagj-1][tagi+2],
									srcvdata[tagj  ][tagi-1],srcvdata[tagj  ][tagi],srcvdata[tagj  ][tagi+1],srcvdata[tagj  ][tagi+2],
									srcvdata[tagj+1][tagi-1],srcvdata[tagj+1][tagi],srcvdata[tagj+1][tagi+1],srcvdata[tagj+1][tagi+2],
									srcvdata[tagj+2][tagi-1],srcvdata[tagj+2][tagi],srcvdata[tagj+2][tagi+1],srcvdata[tagj+2][tagi+2],
									(coordsX[i]-srcxdef[tagi])/(srcxdef[tagi+1]-srcxdef[tagi]),
									(coordsY[j]-srcydef[tagj])/(srcydef[tagj+1]-srcydef[tagj]),undef
								);
							}
						}
						
						cdws.writeData(desv);
					}
				}
			}
			
			break;
			
		default: throw new UnsupportedOperationException("unsupported type for "+type);
		}
		
		cdrs.closeFile();	cdws.closeFile();
		
		System.out.println("Finish horizontal interpolating.");
	}
	
	
	/**
     * vertical interpolation (natural spline)
     *
     * @param	path	path of file after interpolation
     * @param	z		z count after interpolation
     * @param	vs		variables which need linear interpolation (e.g., relative humidity)
     */
	public void verticalInterp(String path,int z,String... vs){
		System.out.println("Start vertical interpolating...");
		
		int sx=srcdata.getXCount(),sy=srcdata.getYCount(),sz=srcdata.getZCount();
		int tcount=srcdata.getTCount();
		
		// destinated data for writting
		Variable desv=new Variable(null,true,new Range(1,z,sy,sx));
		
		DataWrite cdws=DataIOFactory.getDataWrite(srcdata,path);
		
		float[] zdefS=srcdata.getZDef().getSamples();
		float[] ydef =srcdata.getYDef().getSamples();
		float[] xdef =srcdata.getXDef().getSamples();
		float[] zdefD=new float[z];
		
		zdefD[0]=zdefS[0];
		float increS=(zdefS[sz-1]-zdefS[0])/(zdefS.length-1);
		float increD=(zdefS[sz-1]-zdefS[0])/(z-1);
		for(int k=1;k<z;k++) zdefD[k]=zdefD[0]+k*increD;
		
		Spline spln=new Spline(zdefS,new float[sz]);
		GridDataFetcher gdf=new GridDataFetcher(srcdata);
		
		for(int l=0;l<tcount;l++){
			System.out.println("  processing "+srcdata.getTDef().getSamples()[l].toString()+" ...");
			
			for(Var vv:srcdata.getVDef())
			if(vv.getZCount()!=1){
				float[][][] buffer=gdf.prepareXYZBuffer(vv.getName(),l+1);
				float[][][] dedata=desv.getData()[0];
				
				if(inNameList(vv.getName(),vs)){ // linear interpolation
					for(int j=0;j<sy;j++){ float lat=ydef[j];
					for(int i=0;i<sx;i++){ float lon=xdef[i];
						for(int k=0;k<z;k++) dedata[k][j][i]=gdf.fetchXYZBuffer(lon,lat,zdefD[k],buffer);
					}}
					
				}else{ // spline interpolation with 1st BC
					for(int j=0;j<sy;j++)
					for(int i=0;i<sx;i++){
						float[] bin =new float[sz];
						float[] bout=new float[z ];
						
						for(int k=0;k<sz;k++) bin[k]=buffer[k][j][i];
						
						spln.updateData(bin);
						if(increS>0)
							spln.cubicSplineWith1stBC((bin[1]-bin[0])/(zdefS[1]-zdefS[0]),(bin[sz-1]-bin[sz-2])/(zdefS[sz-1]-zdefS[sz-2]));
						else
							spln.cubicSplineWith1stBC((bin[sz-1]-bin[sz-2])/(zdefS[sz-1]-zdefS[sz-2]),(bin[1]-bin[0])/(zdefS[1]-zdefS[0]));
						spln.cValues(zdefD,bout);
						
						for(int k=0;k<z;k++) dedata[k][j][i]=bout[k];
					}
				}
				
				cdws.writeData(desv);
				
			}else{ // for single-level data
				float[][] buf=gdf.prepareXYBuffer(vv.getName(),l+1,1);
				
				cdws.writeData(new Variable("",true,new float[][][][]{{buf}}));
			}
		}
		
		cdws.closeFile();
		
		// write ctl
		FileWriteInterface fwi=new FileWriteInterface(IOUtil.getCompleteFileNameWithoutExtension(path)+".ctl");
		
		FileWriter fw=null;	Scanner sn=null;
		
		if(fwi.getFlag()!=SKIP){
			StringBuilder sb=new StringBuilder();
			
			try{
				switch(fwi.getFlag()){
				case RENAME:
					fw=new FileWriter(fwi.getParent()+fwi.getNewName());
					break;
					
				case OVERWRITE:
					fw=new FileWriter(fwi.getFile());
					break;
					
				case APPEND:
					fw=new FileWriter(fwi.getFile(),true);
					break;

				default:
					throw new IllegalArgumentException("unsupported for "+fwi.getFlag());
				}
				
				sn=new Scanner(new File(srcdata.getPath()));
				
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			
			while(sn.hasNextLine()){
				String line=sn.nextLine();
				
				if(line.startsWith("dset")){
					sb.append("dset ^");	sb.append(IOUtil.getFileName(path));
					sb.append("\n");
				
				}else if(line.startsWith("zdef")){
					sb.append("zdef ");	sb.append(z);
					
					if(zdefS[0]-zdefS[1]>0){
						sb.append(" levels ");
						for(int k=0;k<z;k++){
							sb.append(zdefS[0]+increD*k);
							sb.append(" ");
						}
						sb.append("\n");
						
					}else{
						sb.append(" linear ");
						Scanner tmp=new Scanner(line);
						tmp.next();	tmp.next();	tmp.next();
						sb.append(tmp.next());
						sb.append(" ");	sb.append(increD);
						sb.append("\n");	tmp.close();
					}
				
				}else if(line.startsWith("vars")){
					sb.append(line);	sb.append("\n");
					
					Scanner tmp=new Scanner(line);
					
					tmp.next();	int lc=Integer.parseInt(tmp.next());
					
					tmp.close();
					
					for(int i=0;i<lc;i++){
						line=sn.nextLine();
						
						tmp=new Scanner(line);
						
						sb.append(String.format("%-9s",tmp.next()));
						
						if(!tmp.next().equals("0")){
							sb.append(String.format("%3d ",z));
							sb.append(tmp.nextLine()+"\n");
							
						}else{ sb.append("  0 "+tmp.nextLine()+"\n");}
					}
				}
				else if(line.startsWith("options")) continue;
				else sb.append(line+"\n");
			}
			
			sn.close();
			
			try{ fw.write(sb.toString());	fw.close();}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		System.out.println("Finish vertical interpolating.");
	}
	
	
	/**
     * temporal interpolation
     *
     * @param	path	path of file after interpolation
     * @param	type	type of interpolation, e.g. cubic, linear
     * @param	T		T count after interpolation
     */
	public void temporalInterp(String path,Type type,int T){
		DataRead  dr=DataIOFactory.getDataRead(srcdata);
		DataWrite dw=DataIOFactory.getDataWrite(srcdata,path);
		
		Range rone=new Range("z(1,1)",srcdata);
		Range rmul=new Range("",srcdata);
		
		for(Var v:srcdata.getVDef()){
			Range r=null;
			
			if(v.getZCount()!=1) r=rmul;
			else r=rone;
			
			Variable tmp=new Variable(v.getName(),r);
			
			dr.readData(tmp);
			
			Variable res=tmp.interpolateT(T,type);
			
			dw.writeData(res);
		}
		
		dw.writeCtl(srcdata);
		
		dw.closeFile();	dr.closeFile();
	}
	
	
	/**
     * Global Gaussian grid to even lat/lon grid interpolation, meridional interpolation is natural spline
     *
     * @param	path		path of file after interpolation
     * @param	zonalType	type of zonal interpolation
     * @param	y			y count
     * @param	x			x count
     */
	public void GaussianToEvenGridInterp(String path,Type zonalType,int y,int x){
		System.out.println("Start gauss to grid interpolating...");
		
		int sy=srcdata.getYCount();
		
		Range srange=new Range(1,1,sy,srcdata.getXCount());
		Range drange=new Range(1,1,y,x);
		
		Variable srcv=new Variable(null,srange);
		Variable desv=new Variable(null,drange);
		
		float[]   tmp2=new float[y ];
		float[]    nzd=new float[y];
		float[]   ydef=new float[sy+2];
		float[][] buff=new float[sy][x];
		
		ydef[0]=-90;	ydef[sy+1]=90;
		
		for(int j=0;j<sy;j++) ydef[j+1]=srcdata.getYDef().getSamples()[j];
		for(int j=0;j<y ;j++) nzd[j]=-90+180f/(y-1)*j;
		
		DataRead  cdrs=DataIOFactory.getDataRead(srcdata);
		DataWrite cdws=DataIOFactory.getDataWrite(srcdata,path);
		
		float[][] srcvdata=srcv.getData()[0][0];
		float[][] desvdata=desv.getData()[0][0];
		
		for(int l=1;l<=srcdata.getTCount();l++){
			System.out.println("  Processing "+srcdata.getTDef().getSamples()[l-1].toString()+" ...");
			
			int[] p=srcv.getRange().getTRange();	p[0]=l;p[1]=l;
			
			for(int m=0;m<srcdata.getVCount();m++){
				srcv.setName(srcdata.getVDef()[m].getName());
				desv.setName(srcv.getName());
				
				float undef=srcdata.getUndef(srcdata.getVDef()[m].getName());
				
				for(int k=1;k<=srcdata.getVDef()[m].getZCount();k++){
					p=srcv.getRange().getZRange();	p[0]=k;p[1]=k;
					
					cdrs.readData(srcv);
					
					for(int j=0;j<sy;j++) interp1D(srcvdata[j],buff[j],zonalType,undef);
					
					for(int i=0;i<x;i++){
						int defCount=0;
						for(int j=0;j<sy;j++) if(buff[j][i]!=undef) defCount++;
						
						float[] tmp1=new float[defCount+2];
						float[] nydf=new float[defCount+2];
						
						int ptr=1;
						for(int j=0;j<sy;j++) if(buff[j][i]!=undef){
							nydf[ptr]=ydef[j+1];
							tmp1[ptr]=buff[j][i];
							ptr++;
						}
						nydf[0]=ydef[0];	nydf[defCount+1]=ydef[ydef.length-1];
						tmp1[0]=tmp1[1];	tmp1[tmp1.length-1]=tmp1[tmp1.length-2];
						
						Spline spln=new Spline(nydf,tmp1);
						
						spln.cubicSplineWithNaturalBC();
						
						spln.cValues(nzd,tmp2);
						
						for(int j=0;j<y;j++) desvdata[j][i]=tmp2[j];
					}
					
					cdws.writeData(desv);
				}
			}
		}
		
		cdrs.closeFile();	cdws.closeFile();
		
		// write ctl
		FileWriteInterface fwi=new FileWriteInterface(IOUtil.getCompleteFileNameWithoutExtension(path)+".ctl");
		
		FileWriter fw=null;
		
		if(fwi.getFlag()!=SKIP){
			StringBuilder sb=new StringBuilder();
			
			try{
				switch(fwi.getFlag()){
				case RENAME:
					fw=new FileWriter(fwi.getParent()+fwi.getNewName());
					break;
					
				case OVERWRITE:
					fw=new FileWriter(fwi.getFile());
					break;
					
				case APPEND:
					fw=new FileWriter(fwi.getFile(),true);
					break;
	
				default:
					throw new IllegalArgumentException("unsupported flag for "+fwi.getFlag());
				}
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			
			sb.append("dset ^"+IOUtil.getFileName(path)+"\n");
			sb.append("title "+srcdata.getTitle()+"\n");
			sb.append("undef "+srcdata.getUndef(srcdata.getVDef()[0].getName())+"\n");
			sb.append("xdef "+x+" linear "+srcdata.getXDef().getSamples()[0]+" "+360f/x+"\n");
			sb.append("ydef "+y+" linear -90 "+180f/(y-1)+"\n");
			sb.append("zdef "+srcdata.getZCount());
			
			if(srcdata.zLinear())
				sb.append(" linear "+srcdata.getZDef().getSamples()[0]+" "+srcdata.getDZDef()[0]+"\n");
			else{
				sb.append(" levels ");
				
				for(int k=0;k<srcdata.getZCount();k++) sb.append(srcdata.getZDef().getSamples()[k]+" ");
				
				sb.append("\n");
			}
			
			sb.append(
				"tdef "+srcdata.getTCount()+" linear "+
				srcdata.getTDef().getSamples()[0].toGradsDate()+" "+
				srcdata.getTIncrement()+"\n"
			);
			
			sb.append("vars "+srcdata.getVCount()+"\n");
			
			for(int v=0;v<srcdata.getVCount();v++){
				int zc=srcdata.getVDef()[v].getZCount();
				
				sb.append(
					srcdata.getVDef()[v].getName()+" "+
					(zc==1?0:zc)+" 99 "+
					srcdata.getVDef()[v].getComment()+"\n"
				);
			}
			
			sb.append("endvars\n");
			
			try{ fw.write(sb.toString());	fw.close();}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		System.out.println("Finish gauss to grid interpolating.");
	}
	
	
	/*** helper methods ***/
	private boolean inNameList(String name,String[] list){
		if(list==null) return false;
		
		for(String s:list)
		if(name.equalsIgnoreCase(s)) return true;
		
		return false;
	}
	
	
	/** test
	public static void main(String[] args){
		DataInterpolation di=new DataInterpolation(DiagnosisFactory.parseFile(
			"D:/Data/ERAInterim/Keff/PV/PV.ctl"
		).getDataDescriptor());
	}*/
}
