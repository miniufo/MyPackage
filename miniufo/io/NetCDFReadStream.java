/**
 * @(#)NetCDFReadStream.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.IOException;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.NetcdfFile;
import miniufo.descriptor.NetCDFDescriptor;
import miniufo.diagnosis.Range;


/**
 * used to read data from NetCDF file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class NetCDFReadStream implements DataRead,Print{
	//
	private int t=0;
	private int z=0;
	private int y=0;
	private int x=0;
	
	private boolean print=true;
	
	private NetcdfFile       nc=null;
	private NetCDFDescriptor nd=null;
	
	
	/**
     * constructor
     *
     * @param	nd	NetCDF descriptor
     */
	public NetCDFReadStream(NetCDFDescriptor nd){
		this.nd=nd;
		
		try{ nc=NetcdfFile.open(nd.getDSet());}
		catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
	}
	
	
	/**
	 * to read data from the specified file
	 *
     * @param	v	variable need to fill with data
     */ 
	public void readData(miniufo.diagnosis.Variable... vs){
		float offset=0;	float factor=1;	
		
		if(print) System.out.print("\nStart reading ");
		
		for(int m=0;m<vs.length;m++){
			if(print) System.out.print(vs[m].getName()+" ");
			
			Range range=vs[m].getRange();
			
			float[][][][] data=vs[m].getData();
			
			int[] trange=range.getTRange();	t=vs[m].getTCount();
			int[] zrange=range.getZRange();	z=vs[m].getZCount();
			int[] yrange=range.getYRange();	y=vs[m].getYCount();
			int[] xrange=range.getXRange();	x=vs[m].getXCount();
			
			ucar.nc2.Variable v=nc.findVariable(vs[m].getName());
			
			if(v==null)
			throw new NullPointerException("Cannot find "+vs[m].getName()+" in "+nd.getDSet());
			
			Attribute attr=null;
			String comment=null;
			if((attr=v.findAttribute("add_offset"))!=null)
				offset=attr.getNumericValue().floatValue();
			
			if((attr=v.findAttribute("scale_factor"))!=null)
				factor=attr.getNumericValue().floatValue();
			
			if((attr=v.findAttribute("missing_value"))!=null)
				vs[m].setUndef(attr.getNumericValue().floatValue());
			
			if((attr=v.findAttribute("long_name"))!=null)
				comment=attr.getStringValue();
			
			if((attr=v.findAttribute("units"))!=null)
				comment+=" ("+attr.getStringValue()+")";
			
			vs[m].setCommentAndUnit(comment);
			
			try{
		    	if(nd.getDimensionCount()==4){
					int[] ix=new int[4];
					
					final int[] shp={1,1,1,1};
				    
				    ix[0]=trange[0]-1;
				    ix[1]=zrange[0]-1;
				    ix[2]=nd.isYRev()?nd.getYCount()-yrange[0]:yrange[0]-1;
				    ix[3]=xrange[0]-1;
				    
				    if(vs[m].isTFirst()){
				    	if(nd.isYRev()){
						    for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[l][k][j][i]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[3]++;
							}
							ix[3]=xrange[0]-1;					ix[2]--;
							}
							ix[2]=nd.getYCount()-yrange[0];		ix[1]++;
							}
							ix[1]=zrange[0]-1;					ix[0]++;
							}
							
				    	}else{
				    		for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[l][k][j][i]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[3]++;
							}
							ix[3]=xrange[0]-1;					ix[2]++;
							}
							ix[2]=yrange[0]-1;					ix[1]++;
							}
							ix[1]=zrange[0]-1;					ix[0]++;
							}
				    	}
						
				    }else{
				    	if(nd.isYRev()){
							for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[k][j][i][l]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[3]++;
							}
							ix[3]=xrange[0]-1;					ix[2]--;
							}
							ix[2]=nd.getYCount()-yrange[0];		ix[1]++;
							}
							ix[1]=zrange[0]-1;					ix[0]++;
							}
							
				    	}else{
							for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[k][j][i][l]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[3]++;
							}
							ix[3]=xrange[0]-1;					ix[2]++;
							}
							ix[2]=yrange[0]-1;					ix[1]++;
							}
							ix[1]=zrange[0]-1;					ix[0]++;
							}
				    	}
				    }
					
				}else if(nd.getDimensionCount()==3){
					int[] ix=new int[3];
					
					final int[] shp={1,1,1};
					
					if(vs[m].isTFirst()){
						if(!nd.hasTDef()){
						    ix[0]=zrange[0]-1;	ix[1]=nd.getYCount()-yrange[0];	ix[2]=xrange[0]-1;
						    
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[0][k][j][i]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]++;
							}
							ix[2]=xrange[0]-1;					ix[1]--;
							}
							ix[1]=nd.getYCount()-yrange[0];		ix[0]++;
							}
							
						}else if(!nd.hasZDef()){
						    ix[0]=trange[0]-1;	ix[1]=nd.getYCount()-yrange[0];	ix[2]=xrange[0]-1;
						    
							for(int l=0;l<t;l++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[l][0][j][i]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]++;
							}
							ix[2]=xrange[0]-1;					ix[1]--;
							}
							ix[1]=nd.getYCount()-yrange[0];		ix[0]++;
							}
							
						}else if(!nd.hasYDef()){
						    ix[0]=trange[0]-1;	ix[1]=zrange[0]-1;	ix[2]=xrange[0]-1;
				    		
							for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int i=0;i<x;i++){
								data[l][k][0][i]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]++;
							}
							ix[2]=xrange[0]-1;	ix[1]++;
							}
							ix[1]=zrange[0]-1;	ix[0]++;
							}
							
						}else{
						    ix[0]=trange[0]-1;	ix[1]=zrange[0]-1;	ix[2]=nd.getYCount()-yrange[0];
				    		
							for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
								data[l][k][j][0]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]--;
							}
							ix[2]=nd.getYCount()-yrange[0];		ix[1]++;
							}
							ix[1]=zrange[0]-1;					ix[0]++;
							}
						}
						
					}else{
						if(!nd.hasTDef()){
						    ix[0]=zrange[0]-1;	ix[1]=nd.getYCount()-yrange[0];	ix[2]=xrange[0]-1;
				    		
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[k][j][i][0]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]++;
							}
							ix[2]=xrange[0]-1;					ix[1]--;
							}
							ix[1]=nd.getYCount()-yrange[0];	ix[0]++;
							}
							
						}else if(!nd.hasZDef()){
						    ix[0]=trange[0]-1;	ix[1]=nd.getYCount()-yrange[0];	ix[2]=xrange[0]-1;
				    		
							for(int l=0;l<t;l++){
							for(int j=0;j<y;j++){
							for(int i=0;i<x;i++){
								data[0][j][i][l]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]++;
							}
							ix[2]=xrange[0]-1;					ix[1]--;
							}
							ix[1]=nd.getYCount()-yrange[0];		ix[0]++;
							}
							
						}else if(!nd.hasYDef()){
						    ix[0]=trange[0]-1;	ix[1]=zrange[0]-1;	ix[2]=xrange[0]-1;
				    		
							for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int i=0;i<x;i++){
								data[k][0][i][l]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]++;
							}
							ix[2]=xrange[0]-1;	ix[1]++;
							}
							ix[1]=zrange[0]-1;	ix[0]++;
							}
							
						}else{
						    ix[0]=trange[0]-1;	ix[1]=zrange[0]-1;	ix[2]=nd.getYCount()-yrange[0];
				    		
							for(int l=0;l<t;l++){
							for(int k=0;k<z;k++){
							for(int j=0;j<y;j++){
								data[k][j][0][l]=v.read(ix,shp).getFloat(0)*factor+offset;	ix[2]--;
							}
							ix[2]=nd.getYCount()-yrange[0];		ix[1]++;
							}
							ix[1]=zrange[0]-1;	ix[0]++;
							}
						}
					}
					
				}else throw new IllegalArgumentException("Invalid dimension count");
				
		    }catch(IOException ex1){
		    	ex1.printStackTrace();	System.exit(0);
		    }catch(IllegalArgumentException ex2){
		    	ex2.printStackTrace();	System.exit(0);
		    }catch(InvalidRangeException ex3){
		    	ex3.printStackTrace();	System.exit(0);
		    }
		    
	    	vs[m].setUndef(nd.getUndef(vs[m].getName()));
		}
		
		if(print) System.out.println("data...\nFinish reading.");
	}
	
	
	/**
	 * whether to print out
	 *
     * @param	print	print or disable print
     */ 
	public void setPrinting(boolean print){ this.print=print;}
	
	
	/**
	 * close file method
     */
	public void closeFile(){
		if(nc!=null){
			try{ nc.close();}
			catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		}
	}
	
	
	/** test
	public static void main(String[] args){
		// Open the file.
		String filename = "//lynn/dataI/NCEP1/Daily4/Surface/shtfl/shtfl.sfc.gauss.2004.nc";
		NetcdfFile dataFile = null;
		
		try{
			dataFile = NetcdfFile.open(filename, null);
			
			// Get the latitude and longitude Variables.
			Variable latVar = dataFile.findVariable("lat");
			
			if (latVar == null) {
				System.out.println("Cant find Variable latitude");
				return;
			}
			
			Variable lonVar = dataFile.findVariable("lon");
			
			if (lonVar == null) {
				System.out.println("Cant find Variable longitude");
				return;
			}
			
			// Get the lat/lon data from the file.
			ArrayFloat.D1 latArray;
			ArrayFloat.D1 lonArray;
			
			System.out.println(latVar.read(new int[]{0},new int[]{1}).getFloat(0));
			
			latArray = (ArrayFloat.D1) latVar.read(new int[]{0},new int[]{3});
			lonArray = (ArrayFloat.D1) lonVar.read(new int[]{0},new int[]{3});
			
			System.out.println(latArray);
			System.out.println(lonArray);
			
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
		
		System.out.println("*** SUCCESS reading example file "+filename);
	}*/
}
