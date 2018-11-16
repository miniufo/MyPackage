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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import miniufo.io.DataWrite;
import miniufo.io.FileWriteInterface;
import miniufo.io.IOUtil;
import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.Var;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.geophysics.atmos.ThermoDynamics;
import miniufo.mathsphysics.Spline;
import static miniufo.io.FileWriteInterface.Solution.*;
import static miniufo.basic.InterpolationModel.Type;
import static miniufo.basic.InterpolationModel.interp1D;
import static miniufo.basic.InterpolationModel.interp2D;
import static miniufo.basic.InterpolationModel.bilinearInterpolation;
import static miniufo.basic.InterpolationModel.bicubicPolynomialInterpolation;
import static miniufo.basic.InterpolationModel.bicubicLagrangeInterpolation;
import static miniufo.diagnosis.SpatialModel.gEarth;
import static miniufo.geophysics.atmos.ThermoDynamics.kappa;
import static miniufo.geophysics.atmos.ThermoDynamics.Pref;


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
	
	private static final double P0k=Math.pow(Pref,kappa);
	
	
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
				Variable buffer=gdf.prepareXYZBuffer(vv.getName(),l+1);
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
						
						for(int k=0;k<sz;k++) bin[k]=buffer.getData()[0][k][j][i];
						
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
				Variable buf=gdf.prepareXYBuffer(vv.getName(),l+1,1);
				
				cdws.writeData(buf);
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
		
		dw.closeFile();	dr.closeFile();
		
		// write ctl
		FileWriteInterface fwi=new FileWriteInterface(IOUtil.getCompleteFileNameWithoutExtension(path)+".ctl");
		
		FileWriter fw=null;	Scanner sn=null;
		
		if(fwi.getFlag()!=SKIP){
			StringBuilder sb=new StringBuilder();
			
			try{
				switch(fwi.getFlag()){
				case RENAME   : fw=new FileWriter(fwi.getParent()+fwi.getNewName()); break;
				case OVERWRITE: fw=new FileWriter(fwi.getFile());                    break;
				case APPEND   : fw=new FileWriter(fwi.getFile(),true);               break;
				default       : throw new IllegalArgumentException("unsupported for "+fwi.getFlag());
				}
				
				sn=new Scanner(new File(srcdata.getPath()));
				
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			
			while(sn.hasNextLine()){
				String line=sn.nextLine();
				
				if(line.toLowerCase().startsWith("dset")){
					sb.append("dset ^");	sb.append(IOUtil.getFileName(path));
					sb.append("\n");
				
				}else if(line.toLowerCase().startsWith("tdef")){
					sb.append("tdef "+T+" linear ");
					
					int  dt=Math.round(srcdata.getDTDef()[0]);
					int ndt=dt*(srcdata.getTCount()-1)/(T-1)/3600;
					
					Scanner tmp=new Scanner(line);
					tmp.next();	tmp.next();	tmp.next();
					sb.append(tmp.next()+" "+ndt+"hr\n");
					tmp.close();
				
				}else if(line.toLowerCase().startsWith("vars")){
					sb.append(line+"\n");
					
					Scanner tmp=new Scanner(line);
					
					tmp.next();	int lc=Integer.parseInt(tmp.next());
					
					tmp.close();
					
					for(int i=0;i<lc;i++){
						line=sn.nextLine();
						
						tmp=new Scanner(line);
						
						sb.append(String.format("%-9s",tmp.next())+"  "+tmp.next()+"  -1,20  ");
						tmp.next();
						sb.append(tmp.nextLine()+"\n");
					}
				}
				else if(line.startsWith("options")) continue;
				else sb.append(line+"\n");
			}
			
			sn.close();
			
			try{ fw.write(sb.toString());	fw.close();}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
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
	
	
	/**
     * Data from isobaric levels are interpolated to isentropic levels, assuming a
     * linear dependence of temperature on log(p) (method c in Ziv and Alpert 1994, JAM).
     * Other variables are assumed to vary linearly with potential temperature and will
     * be linearly interpolated to the new isentropic levels.
     *
     * @param	path	path for interpolated file
     * @param	TPrs	temperature at pressure levels (K)
     * @param	ZPrs	geopotential at pressure levels (m^2 s^-2)
     * @param	Zsfc	surface geopotential (m^2 s^-2)
     * @param	ptDes	destinated theta levels for interpolation
     * @param	vPrs	variables at pressure levels that need to be interpolated
     */
	public void isobaricToIsentropicInterp(String path,String TPrs,String ZPrs,String Zsfc,float[] ptDes,String... vPrs){
		System.out.println("Start isobaric to isentropic interpolating...");
		
		int sx=srcdata.getXCount(),sy=srcdata.getYCount(),sz=srcdata.getZCount(),st=srcdata.getTCount();
		int z=ptDes.length,M=vPrs.length;
		
		// destinated data for writting, [vs..., sgm, p]
		Variable[] desv=new Variable[M+4];
		float[][][][] dedata=new float[M+4][][][];
		for(int m=0;m<M+4;m++){
			desv[m]=new Variable(null,true,new Range(1,z,sy,sx));
			dedata[m]=desv[m].getData()[0];
		}
		
		float[] zdef  =srcdata.getZDef().getSamples().clone();
		float[] lnPSrc=new float[sz];
		for(int k=0;k<sz;k++){
			zdef[k]*=100f; // change hPa to Pa
			lnPSrc[k]=(float)Math.log(zdef[k]);
		}
		
		GridDataFetcher gdf=new GridDataFetcher(srcdata);
		DataWrite cdws=DataIOFactory.getDataWrite(srcdata,path);
		
		for(int l=0;l<st;l++){
			System.out.println("  processing "+srcdata.getTDef().getSamples()[l].toString()+" ...");
			
			// preparing data
			Variable t  =gdf.prepareXYZBuffer(TPrs,l+1  );
			Variable Z  =gdf.prepareXYZBuffer(ZPrs,l+1  );
			Variable Zsf=gdf.prepareXYBuffer (Zsfc,l+1,1);
			Variable[] vsBuf=new Variable[M];
			
			float[][][]   tpdata=t.getData()[0];
			float[][][]   ZZdata=Z.getData()[0];
			float[][][][] bfData=new float[M][][][];
			
			for(int m=0;m<M;m++){
				 vsBuf[m]=gdf.prepareXYZBuffer(vPrs[m],l+1);
				bfData[m]=vsBuf[m].getData()[0];
			}
			
			for(int j=0;j<sy;j++)
			for(int i=0;i<sx;i++){	// loop for each grid
				// preparing 1D vertical arrays
				float[]   tpSrc=new float[sz];	// temperature
				float[]   zzSrc=new float[sz];	// geopotential
				float[]   sfSrc=new float[M ];
				float[][] vsSrc=new float[M ][sz];
				
				// store 1D vertical data
				for(int k=0;k<sz;k++){
					tpSrc[k]=tpdata[k][j][i];
					zzSrc[k]=ZZdata[k][j][i];
					for(int m=0;m<M;m++) vsSrc[m][k]=bfData[m][k][j][i];
				}
				
				for(int m=0;m<M;m++) sfSrc[m]=vsSrc[m][0];
				
				OnePointZMask pnt=new OnePointZMask(zdef[0],tpSrc[0],Zsf.getData()[0][0][j][i],t.getUndef(),zdef,lnPSrc,tpSrc,zzSrc);
				
				float[][] vDes=interpFromPress2Thetas(pnt,ptDes,sfSrc,vsSrc);
				
				for(int m=0;m<M+4;m++)
				for(int k=0;k<z;k++) dedata[m][k][j][i]=vDes[m][k];
			}
			
			for(int m=0;m<M+4;m++) cdws.writeData(desv[m]);
		}
		
		cdws.closeFile();
		
		// write ctl
		FileWriteInterface fwi=new FileWriteInterface(IOUtil.getCompleteFileNameWithoutExtension(path)+".ctl");
		FileWriter fw=null;	Scanner sn=null;
		
		if(fwi.getFlag()!=SKIP){
			StringBuilder sb=new StringBuilder();
			
			try{
				switch(fwi.getFlag()){
				case RENAME   : fw=new FileWriter(fwi.getParent()+fwi.getNewName()); break;
				case OVERWRITE: fw=new FileWriter(fwi.getFile()); break;
				case APPEND   : fw=new FileWriter(fwi.getFile(),true); break;
				default       : throw new IllegalArgumentException("unsupported for "+fwi.getFlag());
				}
				sn=new Scanner(new File(srcdata.getPath()));
				
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			
			while(sn.hasNextLine()){
				String line=sn.nextLine();
				
				if(line.startsWith("dset")){
					sb.append("dset ^"+IOUtil.getFileName(path)+"\n");
					
				}else if(line.startsWith("zdef")){
					sb.append("zdef "+z+" levels ");
					for(int k=0;k<z;k++) sb.append(ptDes[k]+" ");
					sb.append("\n");
				
				}else if(line.startsWith("vars")){
					sb.append("vars "+(vPrs.length+4)+"\n");
					
					Scanner tmp=new Scanner(line);
					tmp.next();	int lc=Integer.parseInt(tmp.next());
					tmp.close();
					
					List<String> vars=new ArrayList<>();
					
					for(int i=0;i<lc;i++) vars.add(sn.nextLine());
					
					for(int m=0;m<M;m++) for(String v:vars)
					if(v.startsWith(vPrs[m]+" ")){ sb.append(v.replace(sz+"",z+"")+"\n"); break;}
					
					sb.append(String.format("%-11s %3d %5s %2s\n","sgm",z,99,"isentropic density (kg m^-2 s^-1)"));
					sb.append(String.format("%-11s %3d %5s %2s\n","p"  ,z,99,"pressure (Pa)"));
					sb.append(String.format("%-11s %3d %5s %2s\n","z"  ,z,99,"geopotential (m^2 s^-2)"));
					sb.append(String.format("%-11s %3d %5s %2s\n","M"  ,z,99,"Montgomery streamfunction (m^2 s^-2)"));
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
     * Data from isobaric level are interpolated to isentropic levels, assuming a
     * linear dependence of temperature on lnp (method c in Ziv and Alpert 1994, JAM).
     * Other variables are assumed to vary linearly with potential temperature and will
     * be linearly interpolated to the new isentropic levels.
     * 
     * Surface variables sfcT or sfcP are used to mask those grid points under the ground,
     * using the extended definitions of Andrews (1983, JAS).
     *
     * @param	path	path for interpolated file
     * @param	TPrs	temperature at pressure levels (K)
     * @param	ZPrs	geopotential at pressure levels (m^2 s^-2)
     * @param	Tsfc	surface temperature (K)
     * @param	Psfc	surface pressure (Pa)
     * @param	Zsfc	surface geopotential (m^2 s^-2)
     * @param	ptDes	destinated theta levels for interpolation (K)
     * @param	vSfc	surface variables corresponding to vPrs for masking
     * @param	vPrs	variables at pressure levels that need to be interpolated
     */
	public void isobaricToIsentropicInterp(String path,String TPrs,String ZPrs,String Tsfc,String Psfc,String Zsfc,float[] ptDes,String[] vSfc,String[] vPrs){
		System.out.println("Start isobaric to isentropic interpolating with surface masking...");
		
		int sx=srcdata.getXCount(),sy=srcdata.getYCount(),sz=srcdata.getZCount(),st=srcdata.getTCount();
		int z=ptDes.length,M=vPrs.length;
		
		if(M!=vSfc.length) throw new IllegalArgumentException("vSfc ("+vSfc.length+") should have the same length as vPrs ("+vPrs.length+")");
		
		// destinated data for writting, [vPrs..., p, sgm]
		Variable[] desv=new Variable[M+4];
		float[][][][] dedata=new float[M+4][][][];
		for(int m=0;m<M+4;m++){
			desv[m]=new Variable(null,true,new Range(1,z,sy,sx));
			dedata[m]=desv[m].getData()[0];
		}
		
		float[] zdef  =srcdata.getZDef().getSamples().clone();
		float[] lnPSrc=new float[sz];
		for(int k=0;k<sz;k++){
			zdef[k]*=100f; // change hPa to Pa
			lnPSrc[k]=(float)Math.log(zdef[k]);
		}
		
		GridDataFetcher gdf=new GridDataFetcher(srcdata);
		DataWrite cdws=DataIOFactory.getDataWrite(srcdata,path);
		
		for(int l=0;l<st;l++){
			System.out.println("  processing "+srcdata.getTDef().getSamples()[l].toString()+" ...");
			
			// preparing data
			Variable t  =gdf.prepareXYZBuffer(TPrs,l+1  );
			Variable Z  =gdf.prepareXYZBuffer(ZPrs,l+1  );
			Variable Tsf=gdf.prepareXYBuffer (Tsfc,l+1,1);
			Variable Psf=gdf.prepareXYBuffer (Psfc,l+1,1);
			Variable Zsf=gdf.prepareXYBuffer (Zsfc,l+1,1);
			Variable[] sfBuf=new Variable[M];
			Variable[] vsBuf=new Variable[M];
			
			float[][][]   tpdata=t.getData()[0];
			float[][][]   ZZdata=Z.getData()[0];
			float[][][]   sfData=new float[M][][];
			float[][][][] vsData=new float[M][][][];
			
			for(int m=0;m<M;m++){
				sfBuf[m]=gdf.prepareXYBuffer (vPrs[m],l+1,1);
				vsBuf[m]=gdf.prepareXYZBuffer(vPrs[m],l+1  );
				
				sfData[m]=sfBuf[m].getData()[0][0];
				vsData[m]=vsBuf[m].getData()[0];
			}
			
			for(int j=0;j<sy;j++)
			for(int i=0;i<sx;i++){	// loop for each grid
				// preparing 1D vertical arrays
				float[]   tpSrc=new float[sz];	// temperature
				float[]   zzSrc=new float[sz];	// geopotential
				float[]   sfSrc=new float[M ];
				float[][] vsSrc=new float[M ][sz];
				
				for(int m=0;m<M;m++) sfSrc[m]=sfData[m][j][i];
				
				// store 1D vertical data
				for(int k=0;k<sz;k++){
					tpSrc[k]=tpdata[k][j][i];
					zzSrc[k]=ZZdata[k][j][i];
					for(int m=0;m<M;m++) vsSrc[m][k]=vsData[m][k][j][i];
				}
				
				OnePointZMask pnt=new OnePointZMask(
					Psf.getData()[0][0][j][i],
					Tsf.getData()[0][0][j][i],
					Zsf.getData()[0][0][j][i],
					t.getUndef(),zdef,lnPSrc,tpSrc,zzSrc
				);
				
				float[][] vDes=interpFromPress2Thetas(pnt,ptDes,sfSrc,vsSrc);
				
				for(int m=0;m<M+4;m++)
				for(int k=0;k<z;k++) dedata[m][k][j][i]=vDes[m][k];
			}
			
			for(int m=0;m<M+4;m++) cdws.writeData(desv[m]);
		}
		
		cdws.closeFile();
		
		// write ctl
		FileWriteInterface fwi=new FileWriteInterface(IOUtil.getCompleteFileNameWithoutExtension(path)+".ctl");
		FileWriter fw=null;	Scanner sn=null;
		
		if(fwi.getFlag()!=SKIP){
			StringBuilder sb=new StringBuilder();
			
			try{
				switch(fwi.getFlag()){
				case RENAME   : fw=new FileWriter(fwi.getParent()+fwi.getNewName()); break;
				case OVERWRITE: fw=new FileWriter(fwi.getFile());                    break;
				case APPEND   : fw=new FileWriter(fwi.getFile(),true);               break;
				default       : throw new IllegalArgumentException("unsupported for "+fwi.getFlag());
				}
				sn=new Scanner(new File(srcdata.getPath()));
				
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
			
			while(sn.hasNextLine()){
				String line=sn.nextLine();
				
				if(line.startsWith("dset")){
					sb.append("dset ^"+IOUtil.getFileName(path)+"\n");
					
				}else if(line.startsWith("zdef")){
					sb.append("zdef "+z+" levels ");
					for(int k=0;k<z;k++) sb.append(ptDes[k]+" ");
					sb.append("\n");
				
				}else if(line.startsWith("vars")){
					sb.append("vars "+(vPrs.length+4)+"\n");
					
					Scanner tmp=new Scanner(line);
					tmp.next();	int lc=Integer.parseInt(tmp.next());
					tmp.close();
					
					List<String> vars=new ArrayList<>();
					
					for(int i=0;i<lc;i++) vars.add(sn.nextLine());
					
					for(int m=0;m<M;m++) for(String v:vars)
					if(v.startsWith(vPrs[m]+" ")){ sb.append(v.replace(sz+"",z+"")+"\n"); break;}
					
					sb.append(String.format("%-11s %3d %5s %2s\n","sgm",z,99,"isentropic density (kg m^-2 s^-1)"));
					sb.append(String.format("%-11s %3d %5s %2s\n","p",z,99,"pressure (Pa)"));
					sb.append(String.format("%-11s %3d %5s %2s\n","z",z,99,"geopotential (m^2 s^-2)"));
					sb.append(String.format("%-11s %3d %5s %2s\n","M",z,99,"Montgomery streamfunction (m^2 s^-2)"));
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
	
	
	/*** helper methods ***/
	
	/**
	 * Whether a variable of a given name is in a given list.
	 * 
	 * @param	name	a given variable name
	 * @param	list	a given name list
	 */
	private static boolean inNameList(String name,String[] list){
		if(list==null) return false;
		
		for(String s:list)
		if(name.equalsIgnoreCase(s)) return true;
		
		return false;
	}
	
	
	/**
	 * Using the Newton-Rapson iteration to solve Eq. (9) of 
	 * Ziv and Alpert (1994, JAM) and find lnp for a given theta.
	 * 
	 * @param	iniGuess	initial guess of lnp (logPa)
	 * @param	theta		a given theta
	 * @param	a			coefficient in Eq. (8) of Ziv and Alpert (1994, JAM)
	 * @param	b			coefficient in Eq. (8) of Ziv and Alpert (1994, JAM)
	 * @param	undef		undefined value
	 */
	private static float iterateToFindLogP(float iniGuess,float theta,float a,float b,float undef){
		if(iniGuess==undef) return undef;
		
		double exnerInv=P0k*Math.exp(-kappa*iniGuess);
		double t=a*iniGuess+b;
		
		// Newton-Rapson interation
		double f=theta-t*exnerInv;
		double fp=exnerInv*(kappa*t-a);
		
		return (float)(iniGuess-(f/fp));
	}
	
	
	/**
	 * Interpolate from pressure levels to specified theta levels with surface masks.
	 * 
	 * @param	pnt		OnePoint data needed for column interpolation
	 * @param	ptDes	destined array of potential temperature (K)
	 * @param	sfcV	all surface variable corresponding to vSrcs
	 * @param	vPres	all variable sources to be interpolated
	 */
	private static float[][] interpFromPress2Thetas(OnePointZMask pnt,float[] ptDes,float[] sfcV,float[][] vPres){
		int zD =ptDes.length;
		int vC =vPres.length;
		int idx=pnt.idx;
		int sL =pnt.pres.length;
		
		float undef=pnt.undef;
		
		float[][] vDes=new float[vC+4][zD];
		float[] sgmDes=vDes[vC  ];
		float[] lnPDes=vDes[vC+1];
		float[] geoDes=vDes[vC+2];
		float[] monDes=vDes[vC+3];
		
		for(int k=0;k<zD;k++){
			if(ptDes[k]>pnt.theta[sL-1]){// mask with undefined value
				sgmDes[k]=undef;
				lnPDes[k]=undef;
				geoDes[k]=undef;
				monDes[k]=undef;
				for(int m=0;m<vC;m++) vDes[m][k]=undef;
				continue;
			}
			
			if(ptDes[k]==pnt.theta[sL-1]){// upper most point
				sgmDes[k]=pnt.sigma[sL-1];
				lnPDes[k]=pnt.pres[sL-1];
				geoDes[k]=pnt.geopt[sL-1];
				monDes[k]=ThermoDynamics.cTemperature(ptDes[k],pnt.pres[sL-1])*ThermoDynamics.Cp+pnt.geopt[sL-1];
				for(int m=0;m<vC;m++) vDes[m][k]=vPres[m][sL-1];
				continue;
			}
			
			int idxC=ArrayUtil.getLEIdxIncre(pnt.theta,ptDes[k]);
			
			if(idxC==-1||idxC<idx){// mask with surface value
				sgmDes[k]=0;
				lnPDes[k]=pnt.pSfc;
				geoDes[k]=pnt.zSfc;
				monDes[k]=ThermoDynamics.cTemperature(ptDes[k],pnt.pSfc)*ThermoDynamics.Cp+pnt.zSfc;
				for(int m=0;m<vC;m++) vDes[m][k]=sfcV[m];
				continue;
				
			}else if(idxC==idx){
				if(ptDes[k]<pnt.ptSfc){// mask with surface value
					sgmDes[k]=0;
					lnPDes[k]=pnt.pSfc;
					geoDes[k]=pnt.zSfc;
					monDes[k]=ThermoDynamics.cTemperature(ptDes[k],pnt.pSfc)*ThermoDynamics.Cp+pnt.zSfc;
					for(int m=0;m<vC;m++) vDes[m][k]=sfcV[m];
					continue;
					
				}else{// interpolate between surface and idxC+1
					float a=(pnt.temp[idxC+1]-pnt.tSfc)/(pnt.lnPrs[idxC+1]-pnt.lPSfc);
					float b= pnt.temp[idxC+1]-pnt.lnPrs[idxC+1]*a;
					
					lnPDes[k]=(pnt.lnPrs[idxC+1]+pnt.lPSfc)/2f;	// initial guess
					
					for(int i=0;i<50;i++){
						float tmp=lnPDes[k];
						lnPDes[k]=iterateToFindLogP(lnPDes[k],ptDes[k],a,b,undef);
						tmp-=lnPDes[k];
						if(tmp==undef||Math.abs(tmp)<1e-6) break;
					}
					
					sgmDes[k]=-(pnt.pres[idxC+1]-pnt.pSfc)/(pnt.theta[idxC+1]-pnt.ptSfc)/gEarth;
					lnPDes[k]=(float)Math.exp(lnPDes[k]);
					geoDes[k]=InterpolationModel.linearInterpolation(pnt.ptSfc,pnt.theta[idxC+1],pnt.zSfc,pnt.geopt[idxC+1],ptDes[k],undef);
					monDes[k]=ThermoDynamics.cTemperature(ptDes[k],lnPDes[k])*ThermoDynamics.Cp+geoDes[k];
					
					for(int m=0;m<vC;m++) vDes[m][k]=InterpolationModel.linearInterpolation(
						pnt.theta[idxC],pnt.theta[idxC+1],sfcV[m],vPres[m][idxC+1],ptDes[k],undef
					);
				}
				
			}else{// interpolate between idxC and idxC+1
				float a=(pnt.temp[idxC+1]-pnt.temp[idxC])/(pnt.lnPrs[idxC+1]-pnt.lnPrs[idxC]);
				float b= pnt.temp[idxC+1]-pnt.lnPrs[idxC+1]*a;
				
				lnPDes[k]=(pnt.lnPrs[idxC+1]+pnt.lnPrs[idxC])/2f;	// initial guess
				
				for(int i=0;i<50;i++){
					float tmp=lnPDes[k];
					lnPDes[k]=iterateToFindLogP(lnPDes[k],ptDes[k],a,b,undef);
					tmp-=lnPDes[k];
					if(tmp==undef||Math.abs(tmp)<1e-6) break;
				}
				
				sgmDes[k]=InterpolationModel.linearInterpolation(pnt.theta[idxC],pnt.theta[idxC+1],pnt.sigma[idxC],pnt.sigma[idxC+1],ptDes[k],undef);
				lnPDes[k]=(float)Math.exp(lnPDes[k]);
				geoDes[k]=InterpolationModel.linearInterpolation(pnt.theta[idxC],pnt.theta[idxC+1],pnt.geopt[idxC],pnt.geopt[idxC+1],ptDes[k],undef);
				monDes[k]=ThermoDynamics.cTemperature(ptDes[k],lnPDes[k])*ThermoDynamics.Cp+geoDes[k];
				
				for(int m=0;m<vC;m++) vDes[m][k]=InterpolationModel.linearInterpolation(
					pnt.theta[idxC],pnt.theta[idxC+1],vPres[m][idxC],vPres[m][idxC+1],ptDes[k],undef
				);
			}
		}
		
		final float threshold=800;
		for(int k=0;k<zD;k++)
		if(sgmDes[k]>threshold){
			if(sgmDes[k+1]<threshold) sgmDes[k]=sgmDes[k+1];
			if(sgmDes[k+1]>threshold&&sgmDes[k+2]<threshold) sgmDes[k]=sgmDes[k+1]=sgmDes[k+2];
		}
		
		return vDes;
	}
	
	
	/**
	 * Help to collect those data needed for interpolate a single column.
	 */
	private static final class OnePointZMask{
		//
		int       idx=-1;	// index that pointed to the last underground point
		
		float   undef=0;	// undefined value
		float    pSfc=0;	// surface pressure (Pa)
		float    tSfc=0;	// surface temperature (K)
		float    zSfc=0;	// surface geopotential (m^2 s^-2)
		float   ptSfc=0;	// surface potential temperature (K)
		float   lPSfc=0;	// surface natural log of pressure
		
		float[] pres =null;	// pressure at each pressure level (Pa)
		float[] temp =null;	// temperature at each pressure level (K)
		float[] geopt=null;	// geopotential at each pressure level (m^2 s^-2)
		float[] theta=null;	// potential temperature at each pressure level (K)
		float[] lnPrs=null;	// natural log of pressure at each pressure level
		float[] sigma=null;	// -g*d(p)/d(theta) at each pressure level, pre-calculated before interpolation
		
		
		/**
		 * Constructor.
		 */
		public OnePointZMask(float pSfc,float tSfc,float zSfc,float undef,float[] pres,float[] lnPrs,float[] temp,float[] geopt){
			int z=pres.length;
			
			this.undef=undef; this.pSfc=pSfc;
			this.tSfc =tSfc ; this.zSfc=zSfc;
			this.lnPrs=lnPrs; this.pres=pres;
			this.geopt=geopt; this.temp=temp;
			
			lPSfc=(float)Math.log(pSfc);
			ptSfc=ThermoDynamics.cPotentialTemperature(tSfc,pSfc);
			
			theta=new float[z];
			sigma=new float[z];
			
			for(int k=0;k<z;k++) theta[k]=ThermoDynamics.cPotentialTemperature(temp[k],pres[k]);
			
			Arrays.sort(theta);
			
			//mono=monotonic(theta);
			
			idx=ArrayUtil.getLEIdxIncre(geopt,zSfc);
			
			interpIsentropicDensity();
		}
		
		
		/**
		 * Check the potential temperature at pressure levels to be monotonic increasing with height.
		private boolean monotonic(float[] theta){
			for(int k=idx+1,K=theta.length-1;k<K;k++)
			if(theta[k+1]<theta[k]) return false;
			
			return true;
		}
		 */
		
		/**
		 * Calculate isentropic density at irregular pressure levels.
		 */
		private void interpIsentropicDensity(){
			int z=theta.length;
			
			float[] prsHalf=new float[z+1];	// half grid pressure
			float[] denHalf=new float[z+1];	// half grid density
			
			for(int k=1;k<z;k++){
				prsHalf[k]= (pres[k]+pres[k-1])/2f;
				denHalf[k]=-(pres[k]-pres[k-1])/(theta[k]-theta[k-1])/gEarth;
			}
			
			prsHalf[0]=pres[0  ]; denHalf[0]=denHalf[1];
			prsHalf[z]=pres[z-1]; denHalf[z]=denHalf[z-1];
			
			// interpolate back to original pressure levels and stored in isenDen array
			InterpolationModel.interp1D(prsHalf,denHalf,pres,sigma,Type.LINEAR,undef,false);
		}
	}
	
	
	/** test
	public static void main(String[] args){System.out.println(ThermoDynamics.kappa);
		DataInterpolation di=new DataInterpolation(DiagnosisFactory.parseFile(
			"D:/Data/ERAInterim/BKGState/OriginalData/DataPrs.ctl"
		).getDataDescriptor());
		
		di.isobaricToIsentropicInterp(
			"D:/Data/ERAInterim/BKGState/OriginalData/DataPTInterp.dat","t","z","zsfc",
			new float[]{265,275,285,300,315,330,350,370,395,430,475,530,600,700,850},new String[]{"u","v"});
		
		di.isobaricToIsentropicInterp(
			"D:/Data/ERAInterim/BKGState/OriginalData/DataPTInterpMask.dat","t","z","tsfc","psfc","zsfc",
			new float[]{265,275,285,300,315,330,350,370,395,430,475,530,600,700,850},
			new String[]{"u10","v10"},new String[]{"u","v"});
	}*/
}
