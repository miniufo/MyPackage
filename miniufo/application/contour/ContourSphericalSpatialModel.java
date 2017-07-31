/**
 * @(#)ContourSphericalSpatialModel.java	1.0 2017.07.02
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import miniufo.application.GeoFluidApplication.BoundaryCondition;
import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.mathsphysics.MathsPhysicsUtil;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


/**
 * This class contains the contour-related algorithms in spherical coordinates.
 *
 * @version 1.0, 2015.06.19
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ContourSphericalSpatialModel extends ContourSpatialModel{
	//
	private double sinN=0;	// sin(N) where N is the latitude of the north boundary
	private double sinS=0;	// sin(S) where S is the latitude of the south boundary
	
	
	/**
     * constructor
     *
     * @param	dd		a DataDescriptor (lat/lon descriptor)
     */
	public ContourSphericalSpatialModel(DataDescriptor dd){
		super(dd);
		
		if(dd instanceof CtsDescriptor)
		throw new IllegalArgumentException("DataDescriptor is a CtsDescriptor");
		
		if(dd.isPeriodicX()){ rngX=360.0; BCx=BoundaryCondition.Periodic;}
		else rngX=dd.getXDef().getRange();
		
		rngY=dd.getYDef().getRange();
		
		if(rngX>360) rngX=360;
		if(rngY>180) rngY=180;
		
		sinN =Math.sin(Math.toRadians(dd.getYDef().getLast ()));
		sinS =Math.sin(Math.toRadians(dd.getYDef().getFirst()));
		areaT=Math.toRadians(rngX)*(sinN-sinS)*EARTH_RADIUS*EARTH_RADIUS;
	}
	
	
	/*** helper methods and classes ***/
	
	/**
	 * compute the area element dS for a possibly-refined grid
	 * 
	 * @param	idx		zonal index
	 * @param	idy		meridional index
	 */
	protected double computeDS(int idx,int idy){
		float[] xdef=xdefS; float[] dlon=dxdefS;
		float[] ydef=ydefS; float[] dlat=dydefS;
		
		//if(idx<0||idx>=xdef.length) throw new IllegalArgumentException("idx should be in [0, "+xdef.length+")");
		//if(idy<0||idy>=ydef.length) throw new IllegalArgumentException("idy should be in [0, "+ydef.length+")");
		
		double lon1,lat1,lon2,lat2;
		
		switch(BCx){
		case Fixed:
			if(idx==0){
				lon1=xdef[idx];
				lon2=xdef[idx]+dlon[idx]/2.0;
			}else if(idx==xdef.length-1){
				lon1=xdef[idx]-dlon[idx-1]/2.0;
				lon2=xdef[idx];
			}else{
				lon1=xdef[idx]-dlon[idx-1]/2.0;
				lon2=xdef[idx]+dlon[idx]/2.0;
			}
			break;
			
		case Periodic:
			if(idx==0){
				lon1=xdef[idx]-dlon[idx]/2.0;
				lon2=xdef[idx]+dlon[idx]/2.0;
			}else if(idx==xdef.length-1){
				lon1=xdef[idx]-dlon[idx-1]/2.0;
				lon2=xdef[idx]+dlon[idx-1]/2.0;
			}else{
				lon1=xdef[idx]-dlon[idx-1]/2.0;
				lon2=xdef[idx]+dlon[idx]/2.0;
			}
			break;

		default: throw new IllegalArgumentException("unsupported BCx: "+BCx);
		}
		
		switch(BCy){
		case Fixed:
			if(idy==0){ // south pole
				lat1=ydef[idy];
				//lat2=ydef[idy]+dlat[idy]/2.0;
				lat2=ydef[idy]+dlat[idy];
			}else if(idy==ydef.length-1){ // north pole
				lat1=ydef[idy];
				//lat1=ydef[idy]-dlat[idy-1]/2.0;
				lat2=ydef[idy];
			}else{
				//lat1=ydef[idy]-dlat[idy-1]/2.0;
				//lat2=ydef[idy]+dlat[idy]/2.0;
				lat1=ydef[idy];
				lat2=ydef[idy]+dlat[idy];
			}
			break;
			
		case Periodic:
			if(idy==0){
				lat1=ydef[idy]-dlat[idy]/2.0;
				lat2=ydef[idy]+dlat[idy]/2.0;
			}else if(idy==ydef.length-1){
				lat1=ydef[idy]-dlat[idy-1]/2.0;
				lat2=ydef[idy]+dlat[idy-1]/2.0;
			}else{
				lat1=ydef[idy]-dlat[idy-1]/2.0;
				lat2=ydef[idy]+dlat[idy]/2.0;
			}
			break;
			
		default: throw new IllegalArgumentException("unsupported BCx: "+BCx);
		}
		
		return (4.0*Math.PI*EARTH_RADIUS*EARTH_RADIUS)*MathsPhysicsUtil.cAreaQuadByRadian(
			Math.toRadians(lon1),Math.toRadians(lat1),
			Math.toRadians(lon2),Math.toRadians(lat2)
		);
	}
	
	/**
	 * Compute the equivalent latitude (degree) given an area of a regional domain.
	 * 
	 * For global area, using the following formula:
	 * A = 2 * PI * R^2 * [1 - sin(lat)]
	 * 
	 * @param	area	a given area (m^2)
	 * @param	undef	undefined value
	 */
	protected double cEquivalentY(double area,float undef){
		if(area==undef) return undef;
		else return Math.toDegrees(Math.asin(sinN-area/(Math.toRadians(rngX)*EARTH_RADIUS*EARTH_RADIUS)));
	}
	
	/**
     * Compute squared gradient in latitude/longitude coordinates,
     * i.e., |grad(tracer)|^2
     * 
     * @return	grd		squared gradient (a scalar)
     */
	protected Variable cSquaredTracerGradient(){
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		DynamicMethodsInSC dm=new DynamicMethodsInSC(ssm);
		
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
