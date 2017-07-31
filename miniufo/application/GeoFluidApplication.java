/**
 * @(#)GeoFluidApplication.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application;

import miniufo.diagnosis.MDate;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;


/**
 * application of weather analysis
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class GeoFluidApplication{
	// four dimensions of the sub-domain built on a domain described by a ctl file
	protected int t=0;
	protected int z=0;
	protected int y=0;
	protected int x=0;
	
	// sub-domain offsets from the complete domain
	protected int tstart=1;
	protected int zstart=1;
	protected int ystart=1;
	protected int xstart=1;
	
	//
	protected float undef=Float.NaN;
	
	// delta elements in four dimensions
	protected float dt;				// unit: s
	protected float dz;				// unit: m or Pa
	protected float dy;				// unit: m
	protected float dx;				// unit: m
	
	// grid definitions
	protected MDate[] tdef=null;	// times (UTC)
	protected float[] zdef=null;	// levels(Pa)
	protected float[] ydef=null;	// latitude (radian)
	protected float[] xdef=null;	// longitude(radian)
	
	protected SpatialModel sm=null;
	
	protected BoundaryCondition BCt=BoundaryCondition.Fixed;
	protected BoundaryCondition BCz=BoundaryCondition.Fixed;
	protected BoundaryCondition BCy=BoundaryCondition.Fixed;
	protected BoundaryCondition BCx=BoundaryCondition.Fixed;
	
	public enum BoundaryCondition{
		Fixed,		// values are fixed and the 1st-order derivative is approximated by 1st-order forward/backward difference (default)
		Periodic,	// both values and the derivatives are periodic
		Expanded	// values at the BC are assigned with those at nearest inner grids
	};
	
	
	/**
     * constructor
     *
     * @param	sm	spacial model used for the calculation
     */
	public GeoFluidApplication(SpatialModel sm){
		this.sm=sm;
		
		dt=sm.getDT(); tdef=sm.getTDef().getSamples();
		dz=sm.getDZ(); zdef=sm.getZDef().getSamples();
		dy=sm.getDY(); ydef=sm.getYDef().getSamples();
		dx=sm.getDX(); xdef=sm.getXDef().getSamples();
	}
	
	
	/*** getor and setor ***/
	public SpatialModel getSpatialModel(){ return sm;}
	
	public BoundaryCondition getBCofT(){ return BCt;}
	
	public BoundaryCondition getBCofZ(){ return BCz;}
	
	public BoundaryCondition getBCofY(){ return BCy;}
	
	public BoundaryCondition getBCofX(){ return BCx;}
	
	public void setBCofT(BoundaryCondition BCt){ this.BCt=BCt;}
	
	public void setBCofZ(BoundaryCondition BCz){ this.BCz=BCz;}
	
	public void setBCofY(BoundaryCondition BCy){ this.BCy=BCy;}
	
	public void setBCofX(BoundaryCondition BCx){ this.BCx=BCx;}
	
	
	/*** helper methods ***/
	protected void assignSubDomainParams(Variable v){
		undef=v.getUndef();
		
		t=v.getTCount(); tstart=v.getRange().getTRange()[0];
		z=v.getZCount(); zstart=v.getRange().getZRange()[0];
		y=v.getYCount(); ystart=v.getRange().getYRange()[0];
		x=v.getXCount(); xstart=v.getRange().getXRange()[0];
	}
	
	protected void checkDimensions(Variable v1,Variable... vs){
		for(Variable v:vs) if(!v1.isLike(v))
		throw new IllegalArgumentException("dimensions not same for:\n"+v1+"\n"+v);
	}
	
	
	/**
     * Standard 2nd-order RK integrator for solving dData/dh = deriv,
     * where derivatives are estimated using finite difference,
     * started with a given initial value.
     *
     * @param	init		initial value to start with
     * @param	deriv		1st-order derivatives
     * @param	h			step
     * @param	increIdx	whether the index is increment or decrement
     */
	protected float[] RK2(float init,float[] deriv,float h,boolean increIdx){
		int len=deriv.length;
		
		if(h==0) throw new IllegalArgumentException("h should be non-zero");
		
		float[] re=deriv.clone();
		
		re[0]=re[len-1]=init;
		
		if(increIdx) for(int i=0,I=len-1;i<I;i++){
			re[i+1]=re[i]+h*(deriv[i]+deriv[i+1])/2f;
			
		}else for(int i=len-1;i>=1;i--){
			re[i-1]=re[i]+h*(deriv[i]+deriv[i-1])/2f;
		}
		
		return re;
	}
}
