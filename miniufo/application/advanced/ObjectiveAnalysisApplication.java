/**
 * @(#)ObjectiveAnalysisApplication.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.advanced;

import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.ObjectiveAnalysis;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;


/**
 * objective analysis application
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ObjectiveAnalysisApplication{
	//
	private CtlDescriptor ctl=null;	// descriptor of the grids
	private CsmDescriptor csm=null;	// descriptor of the stns
	
	private static final float[] dradii={10,7,4,2,1};	// default analysis radii
	
	
	/**
     * constructor
     *
     * @param	dd		a data descriptor
     */
	public ObjectiveAnalysisApplication(CtlDescriptor ctl,CsmDescriptor csm){
		this.ctl=ctl;	this.csm=csm;
		
		if(!ctl.xLinear()||!ctl.yLinear())
			throw new IllegalArgumentException("not linear xdef and ydef");
	}
	
	
	/**
     * objective analysis by Cressman, like oacres in GrADS
     *
     * @param	grid		grid variable
     * @param	station		station variable
     * @param	radii		radii for analysis, in grid space unit
     * 						default: 10, 7, 4, 2, 1
     */
	public void oacres(Variable grid,Variable station,float... radii){
		System.out.println("\nStart objective analysis...");
		
		if(grid.isTFirst()!=station.isTFirst())
			throw new IllegalArgumentException("TFirst not same between grid and station variables");
		
		int tcount=grid.getTCount(),zcount=grid.getZCount();
		int ycount=grid.getYCount(),xcount=grid.getXCount();
		
		int sxcount=station.getXCount(),sycount=station.getYCount();
		int scount =sxcount*sycount;
		
		float[] r=radii;
		
		if(station.getTCount()!=tcount||station.getZCount()!=zcount)
			throw new IllegalArgumentException("t and z of grids are not the same as stations");
		
		if(ctl.getTCount()!=tcount||ctl.getZCount()!=zcount)
			throw new IllegalArgumentException("t and z of grids are not the same as the ctl");
		
		if(r==null) r=dradii.clone();
		else if(r.length>30||r.length==0){
			System.out.println("  Warning: invalid radii, use default (10,7,4,2,1)");
			r=dradii.clone();
		}
		
		float dgrid=(float)Math.toRadians(ctl.getDYDef()[0])*SpatialModel.EARTH_RADIUS;
		for(int i=0;i<r.length;i++) r[i]*=dgrid;
		
		grid.setUndef(station.getUndef());
		
		float[][][] slats=csm.getLat();	float[][][][] gdata=grid.getData();
		float[][][] slons=csm.getLon();	float[][][][] sdata=station.getData();
		
		float[] glats=ctl.getYDef().getSamples();
		float[] glons=ctl.getXDef().getSamples();
		
		float[][]   stn=new float[3][scount];
		float[][][] grd=new float[3][ycount][xcount];
		ObjectiveAnalysis oa=new ObjectiveAnalysis(grd,stn,station.getUndef());
		
		if(grid.isTFirst()){
			/*** prepare for the grid data ***/
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++){
				grd[0][j][i]=glons[i];
				grd[1][j][i]=glats[j];
			}
			
			for(int l=0;l<tcount;l++)
			for(int k=0;k<zcount;k++){
				System.out.println("  analyzing for t="+(l+1)+" z="+(k+1)+"...");
				
				/*** prepare for the station data ***/
				for(int j=0;j<sycount;j++)
				for(int i=0;i<sxcount;i++){
					stn[0][j*sxcount+i]=slons[l][j][i];
					stn[1][j*sxcount+i]=slats[l][j][i];
					stn[2][j*sxcount+i]=sdata[l][k][j][i];
				}
				
				oa.SCMByCressman(r);
				
				/*** storing grid data ***/
				for(int j=0;j<ycount;j++)
				for(int i=0;i<xcount;i++) gdata[l][k][j][i]=grd[2][j][i];
			}
			
		}else{
			/*** prepare for the grid data ***/
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++){
				grd[0][j][i]=glons[i];
				grd[1][j][i]=glats[j];
			}
			
			for(int l=0;l<tcount;l++)
			for(int k=0;k<zcount;k++){
				System.out.println("  analyzing for t="+(l+1)+" z="+(k+1)+"...");
				
				/*** prepare for the station data ***/
				for(int j=0;j<sycount;j++)
				for(int i=0;i<sxcount;i++){
					stn[0][j*sxcount+i]=slons[l][j][i];
					stn[1][j*sxcount+i]=slats[l][j][i];
					stn[2][j*sxcount+i]=sdata[k][j][i][l];
				}
				
				oa.SCMByCressman(r);
				
				/*** storing grid data ***/
				for(int j=0;j<ycount;j++)
				for(int i=0;i<xcount;i++) gdata[k][j][i][l]=grd[2][j][i];
			}
		}
		
		System.out.println("Finished.");
	}
	
	public void oacres(Variable grid,Variable station){ oacres(grid,station,null);}
	
	
	/** test
	public static void main(String[] args){
		CsmDescriptor dd=new CsmDescriptor("d:/Data/DiagnosisVortex/Tokage/tokage.csm");
		
		Range sr=new Range("z(1,1);t(1,1)",dd);
		
		Variable u=new Variable("u",sr);
		Variable v=new Variable("v",sr);
		
		DataRead dr=DataIOFactory.getDataRead(dd);
		dr.readData(u,v);	dr.closeFile();
		
		CtlDescriptor ctl=new CtlDescriptor("d:/cressman/grid.ctl");
		
		Range gr=new Range("",ctl);
		
		Variable gu=new Variable("u",gr);
		Variable gv=new Variable("v",gr);
		
		ObjectiveAnalysisApplication oaa=new ObjectiveAnalysisApplication(ctl,dd);
		
		oaa.oacres(gu,u);
		oaa.oacres(gv,v);
		
		DataWrite dw=DataIOFactory.getDataWrite(ctl,"d:/cressman/oa.dat");
		dw.writeData(ctl,gu,gv);	dw.closeFile();
	}*/
}
