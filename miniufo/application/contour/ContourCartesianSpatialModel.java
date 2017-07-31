/**
 * @(#)ContourCartesianSpatialModel.java	1.0 2017.07.02
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import miniufo.application.GeoFluidApplication.BoundaryCondition;
import miniufo.application.basic.DynamicMethodsInCTS;
import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.diagnosis.Variable;


/**
 * This class contains the contour-related algorithms in Cartesian coordinates.
 *
 * @version 1.0, 2017.07.02
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ContourCartesianSpatialModel extends ContourSpatialModel{
	
	/**
     * constructor
     *
     * @param	dd		a DataDescriptor (lat/lon descriptor)
     */
	public ContourCartesianSpatialModel(DataDescriptor dd){
		super(dd);
		
		if(!(dd instanceof CtsDescriptor))
		throw new IllegalArgumentException("DataDescriptor is not a CtsDescriptor");
		
		rngX=dd.getXDef().getRange();
		rngY=dd.getYDef().getRange();
		
		if(BCx==BoundaryCondition.Periodic) rngX=dd.getXDef().getIncrements()[0]*dd.getXCount();
		
		areaT=rngX*rngY;
	}
	
	
	/*** helper methods and classes ***/
	
	/**
	 * compute the area element dS for a possibly-refined grid
	 * 
	 * @param	idx		zonal index
	 * @param	idy		meridional index
	 */
	protected double computeDS(int idx,int idy){
		int X=xdefS.length;
		int Y=ydefS.length;
		
		float dX=dxdefS[0]; // assume grids are uniform
		float dY=dydefS[0]; // assume grids are uniform
		
		//if(idx<0||idx>=X) throw new IllegalArgumentException("idx should be in [0, "+X+")");
		//if(idy<0||idy>=Y) throw new IllegalArgumentException("idy should be in [0, "+Y+")");
		
		double dx=0,dy=0;
		
		switch(BCx){
		case Fixed:
			if(idx==0||idx==X-1) dx=dX/2.0;
			else dx=dX;
			break;
		case Periodic:
			dx=dX;
			break;
		default: throw new IllegalArgumentException("unsupported BCx: "+BCx);
		}
		
		switch(BCy){
		case Fixed:
			if(idy==0||idy==Y-1) dy=dY/2.0;
			else dy=dY;
			break;
		case Periodic:
			dy=dY;
			break;
		default: throw new IllegalArgumentException("unsupported BCx: "+BCx);
		}
		
		return dx*dy;
	}
	
	/**
	 * Compute the equivalent Y given an area enclosed by a tracer contour.
	 * 
	 * For a rectangle domain, using the following formula:
	 * A = Y0 - area / X0
	 * 
	 * @param	area	a given area (m^2)
	 * @param	undef	undefined value
	 */
	protected double cEquivalentY(double area,float undef){
		if(area==undef) return undef;
		else return rngY-area/rngX;
	}
	
	/**
     * Compute squared gradient in Cartesian coordinates,
     * i.e., |grad(tracer)|^2
     * 
     * @return	grd		squared gradient (a scalar)
     */
	protected Variable cSquaredTracerGradient(){
		CartesianSpatialModel csm=new CartesianSpatialModel((CtsDescriptor)dd);
		DynamicMethodsInCTS dm=new DynamicMethodsInCTS(csm);
		
		dm.setBCofX(BCx);
		
		Variable[] grd=dm.c2DGradient(tracer);
		
		grd[0].squareEq();
		grd[1].squareEq();
		
		Variable grdmag2=grd[0].plusEq(grd[1]);
		
		grdmag2.setName("grd2");
		grdmag2.setCommentAndUnit("squared gradient of "+tracer.getName()+" ("+tracer.getUnit()+"^2 m^-2)");
		
		return grdmag2;
	}
	
	
	/*** test **
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/ERAInterim/Keff/PV/GRDSqr.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] vs=df.getVariables(new Range("t(1,1)",dd),"pv","grd2pv");
		
		Variable pv =vs[0].multiplyEq(1e-6f);
		Variable grd=vs[1];
		
		ContourCoords cc=new ContourCoords(dd);
		
		Contours[][] ctss=cc.getContourInfo(pv,51);System.out.println(ctss[0][0]);
		
		Variable area=cc.integrateOverContour(null,pv,ctss);
		Variable grd2=cc.integrateOverContour(grd,pv,ctss);
		Variable grdA=cc.cGradientWRTArea(area,ctss);
		Variable Le2 =cc.cEquivalentLengthSquare(grd2,area,ctss);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"D:/Data/ERAInterim/Keff/PV/re.dat");
		dw.writeData(dd,area,grd2,grdA,Le2); dw.closeFile();
	}*/
}
