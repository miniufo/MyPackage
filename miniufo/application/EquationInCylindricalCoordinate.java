/**
 * @(#)EquationInCylindricalCoordinate.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application;

import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.tan;
import static miniufo.diagnosis.SpatialModel.REarth;


/**
 * equation application in cylindrical coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class EquationInCylindricalCoordinate extends GeoFluidApplication{
	//
	protected float[] dxs =null;	// dx scaled by sin(beta)
	protected float[] rs  =null;	// radial distances = R * sin(beta)
	protected float[] bsin=null;	// sin(beta)
	protected float[] bcos=null;	// cos(beta)
	
	
	/**
     * constructor
     *
     * @param	psm		a given spacial model in polar coordinate
     */
	public EquationInCylindricalCoordinate(CylindricalSpatialModel csm){
		super(csm);
		
		dxs =csm.getDXs();
		bsin=csm.getBSin();
		bcos=csm.getBCos();
		
		rs=new float[csm.getYCount()];
		
		for(int j=0,J=csm.getYCount();j<J;j++) rs[j]=REarth*bsin[j];
		
		BCx=BoundaryCondition.Periodic;
	}
	
	
	/**
	 * the compute derivative (gradients) of a variable in a specific dimension
	 * 
	 * @param	var		a given variable
	 * @param	dim		specified dimension
	 */
	public Variable cDerivative(Variable var,Dimension dim){
		assignSubDomainParams(var);
		
		Variable der=new Variable("grd"+var.getName(),var);
		der.setCommentAndUnit("gradient of "+var.getName()+" along "+dim+" dimension");
		der.setValue(undef);
		
		float[][][][] vdata=var.getData();
		float[][][][] ddata=der.getData();
		
		switch(dim){
		case X:{
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					if(BCx==BoundaryCondition.Periodic){
						if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][x-1]!=undef)
						ddata[l][k][j][0  ]=(vdata[l][k][j][1]-vdata[l][k][j][x-1])/(dxs[ystart-1+j]*2);
						if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-2]!=undef)
						ddata[l][k][j][x-1]=(vdata[l][k][j][0]-vdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
					}
					
					for(int i=1;i<x-1;i++) if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef)
					ddata[l][k][j][i]=(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					if(BCx==BoundaryCondition.Periodic){
						if(vdata[k][j][1][l]!=undef&&vdata[k][j][x-1][l]!=undef)
						ddata[k][j][0  ][l]=(vdata[k][j][1][l]-vdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
						if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-2][l]!=undef)
						ddata[k][j][x-1][l]=(vdata[k][j][0][l]-vdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
					}
					
					for(int i=1;i<x-1;i++) if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef)
					ddata[k][j][i][l]=(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
				}
			}
			break;
		}
		case Y:{
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					if(BCy==BoundaryCondition.Fixed){
						if(vdata[l][k][  1][i]!=undef&&vdata[l][k][  0][i]!=undef)
						ddata[l][k][0  ][i]=(vdata[l][k][  1][i]-vdata[l][k][  0][i])/dy;
						if(vdata[l][k][y-1][i]!=undef&&vdata[l][k][y-2][i]!=undef)
						ddata[l][k][y-1][i]=(vdata[l][k][y-1][i]-vdata[l][k][y-2][i])/dy;
					}
					
					for(int j=1;j<y-1;j++) if(vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
					ddata[l][k][j][i]=(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					if(BCy==BoundaryCondition.Fixed){
						if(vdata[k][  1][i][l]!=undef&&vdata[k][  0][i][l]!=undef)
						ddata[k][0  ][i][l]=(vdata[k][  1][i][l]-vdata[k][  0][i][l])/dy;
						if(vdata[k][y-1][i][l]!=undef&&vdata[k][y-2][i][l]!=undef)
						ddata[k][y-1][i][l]=(vdata[k][y-1][i][l]-vdata[k][y-2][i][l])/dy;
					}
					
					for(int j=0;j<y;j++) if(vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
					ddata[k][j][i][l]=(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2);
				}
			}
			break;
		}
		case Z:{
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCz==BoundaryCondition.Fixed){
						if(vdata[l][  1][j][i]!=undef&&vdata[l][  0][j][i]!=undef)
						ddata[l][0  ][j][i]=(vdata[l][  1][j][i]-vdata[l][  0][j][i])/dz;
						if(vdata[l][z-1][j][i]!=undef&&vdata[l][z-2][j][i]!=undef)
						ddata[l][z-1][j][i]=(vdata[l][z-1][j][i]-vdata[l][z-2][j][i])/dz;
					}
					
					for(int k=1;k<z-1;k++) if(vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef)
					ddata[l][k][j][i]=(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCz==BoundaryCondition.Fixed){
						if(vdata[  1][j][i][l]!=undef&&vdata[  0][j][i][l]!=undef)
						ddata[0  ][j][i][l]=(vdata[  1][j][i][l]-vdata[  0][j][i][l])/dz;
						if(vdata[z-1][j][i][l]!=undef&&vdata[z-2][j][i][l]!=undef)
						ddata[z-1][j][i][l]=(vdata[z-1][j][i][l]-vdata[z-2][j][i][l])/dz;
					}
					
					for(int k=1;k<z-1;k++) if(vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef)
					ddata[k][j][i][l]=(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz*2);
				}
			}
			break;
		}
		case T:{
			if(var.isTFirst()){
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCt==BoundaryCondition.Fixed){
						if(vdata[  1][k][j][i]!=undef&&vdata[  0][k][j][i]!=undef)
						ddata[0  ][k][j][i]=(vdata[  1][k][j][i]-vdata[  0][k][j][i])/dt;
						if(vdata[t-1][k][j][i]!=undef&&vdata[t-2][k][j][i]!=undef)
						ddata[t-1][k][j][i]=(vdata[t-1][k][j][i]-vdata[t-2][k][j][i])/dt;
					}
					
					for(int l=1;l<t-1;l++) if(vdata[l+1][k][j][i]!=undef&&vdata[l-1][k][j][i]!=undef)
					ddata[l][k][j][i]=(vdata[l+1][k][j][i]-vdata[l-1][k][j][i])/(dt*2);
				}
				
			}else{
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCt==BoundaryCondition.Fixed){
						if(vdata[k][j][i][  1]!=undef&&vdata[k][j][i][  0]!=undef)
						ddata[k][j][i][0  ]=(vdata[k][j][i][  1]-vdata[k][j][i][  0])/dt;
						if(vdata[k][j][i][t-1]!=undef&&vdata[k][j][i][t-2]!=undef)
						ddata[k][j][i][t-1]=(vdata[k][j][i][t-1]-vdata[k][j][i][t-2])/dt;
					}
					
					for(int l=1;l<t-1;l++) if(vdata[k][j][i][l+1]!=undef&&vdata[k][j][i][l-1]!=undef)
					ddata[k][j][i][l]=(vdata[k][j][i][l+1]-vdata[k][j][i][l-1])/(dt*2);
				}
			}
			break;
		}
		default: throw new IllegalArgumentException("unsupported dimension: "+dim);
		}
		
		return der;
	}
	
	
	/**
     * weighting a variable with sin(beta)
     *
     * @param	v	a given variable
     */
	public Variable weightBSinEq(Variable v){
		assignSubDomainParams(v);
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(vdata[l][k][j][i]!=undef) vdata[l][k][j][i]*=bsin[j];
			
		}else{
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(vdata[k][j][i][l]!=undef) vdata[k][j][i][l]*=bsin[j];
		}
		
		return v;
	}
	
	public Variable weightBSin(Variable v){
		assignSubDomainParams(v);
		
		Variable nv=v.copy();
		
		float[][][][] ndata=nv.getData();
		
		if(v.isTFirst()){
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[l][k][j][i]!=undef) ndata[l][k][j][i]*=bsin[j];
			
		}else{
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[k][j][i][l]!=undef) ndata[k][j][i][l]*=bsin[j];
		}
		
		return nv;
	}
	
	
	/**
     * de-weighting a variable with sin(beta)
     *
     * @param	v	a given variable
     */
	public Variable deWeightBSinEq(Variable v){
		assignSubDomainParams(v);
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vdata[l][k][j][i]!=undef&&bsin[j]!=0) vdata[l][k][j][i]/=bsin[j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vdata[k][j][i][l]!=undef&&bsin[j]!=0) vdata[k][j][i][l]/=bsin[j];
		}
		
		return v;
	}
	
	public Variable deWeightBSin(Variable v){
		assignSubDomainParams(v);
		
		Variable nv=v.copy();
		
		float[][][][] ndata=nv.getData();
		
		if(v.isTFirst()){
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[l][k][j][i]!=undef){
				if(j!=0) ndata[l][k][j][i]/=bsin[j];
				else ndata[k][j][i][l]=0;
			}
			
		}else{
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[k][j][i][l]!=undef){
				if(j!=0) ndata[k][j][i][l]/=bsin[j];
				else ndata[k][j][i][l]=0;
			}
		}
		
		return nv;
	}
	
	
	/**
	 * Calculate storm-relative lat/lon velocity by subtracting homogeneous
	 * velocity given by cu and cv i.e., every grid in cylindrical coordinates
	 * translates at the same speed.
	 * 
	 * Notice that this method gives a azimuthally-mean storm-relative velocity
	 * slightly different from the azimuthally-mean original velocity.  Thus it is
	 * recommended to use cStormRaltiveAziRadVelocity() instead of using this one.
	 * 
	 * @param	cu	moving speed of the system in zonal direction (m s^-1)
	 * @param	cv	moving speed of the system in meridional direction (m s^-1)
	 * @param	u	zonal wind speed in cylindrical coordinates (m s^-1)
	 * @param	v	meridional wind speed in cylindrical coordinates (m s^-1)
	 */
	public void cStormRelativeLatLonVelocity(float[] cu,float[] cv,Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		if(t!=cu.length||t!=cv.length)
		throw new IllegalArgumentException("tlengths are not the same");
		
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		
		if(u.isTFirst())
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				udata[l][k][j][i]-=cu[l];
				vdata[l][k][j][i]-=cv[l];
			}
		else
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				udata[k][j][i][l]-=cu[l];
				vdata[k][j][i][l]-=cv[l];
			}
	}
	
	public void cStormRelativeLatLonVelocity(Variable u,Variable v){
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		cStormRelativeLatLonVelocity(csm.getUWhole(),csm.getVWhole(),u,v);
	}
	
	
	/**
     * Calculate storm-relative azi/rad velocity by subtracting homogeneous
	 * velocity given by cu and cv i.e., every grid in cylindrical coordinates
	 * translates at the same speed.
     *
     * @param	cu	moving speed of the system in zonal direction (m s^-1)
     * @param	cv	moving speed of the system in meridional direction (m s^-1)
     * @param	ut	tangential wind speed (m s^-1)
     * @param	vr	radial wind speed (m s^-1)
     */
    public void cStormRelativeAziRadVelocity(float[] cu,float[] cv,Variable ut,Variable vr){
    	checkDimensions(ut,vr);
    	assignSubDomainParams(ut);
		
		if(t!=cu.length||t!=cv.length)
		throw new IllegalArgumentException("tlengths are not the same");
		
		float[] xdef=((CylindricalSpatialModel)sm).getXDef().getSamples();
    	float[][][][] utdata=ut.getData();
    	float[][][][] vrdata=vr.getData();
    	
    	if(ut.isTFirst()) for(int l=0;l<t;l++){
	    	float cd=(float)(Math.atan2(cv[l],cu[l])-PI/2.0); // angle start from north, counter-clockwise
	    	float cs=(float)Math.hypot(cu[l],cv[l]);
	    	
	    	for(int i=0;i<x;i++){
	    		float angle=cd-xdef[i];
	    		float wtan=cs*(float)sin(angle);
	    		float wnor=cs*(float)cos(angle);
	    		
		   		for(int k=0;k<z;k++)
		   		for(int j=0;j<y;j++){
		   			utdata[l][k][j][i]-=wtan;
		   			vrdata[l][k][j][i]-=wnor;
		   		}
	    	}
	    }
    	else for(int l=0;l<t;l++){
	    	float cd=(float)(Math.atan2(cv[l],cu[l])-PI/2.0); // angle start from north, counter-clockwise
	    	float cs=(float)Math.hypot(cu[l],cv[l]);
	    	
	    	for(int i=0;i<x;i++){
	    		float angle=cd-xdef[i];
	    		float wtan=cs*(float)sin(angle);
	    		float wnor=cs*(float)cos(angle);
	    		
		   		for(int k=0;k<z;k++)
		   		for(int j=0;j<y;j++){
		   			utdata[k][j][i][l]-=wtan;
		   			vrdata[k][j][i][l]-=wnor;
		   		}
	    	}
	    }
    }
	
    public void cStormRelativeAziRadVelocity(Variable ut,Variable vr){
    	CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
    	
    	cStormRelativeAziRadVelocity(csm.getUWhole(),csm.getVWhole(),ut,vr);
    }
    
	
	/*** helper methods ***/
	protected void assignSubDomainParams(Variable v){
		super.assignSubDomainParams(v);
		
		if(xstart!=1||ystart!=1)
		throw new IllegalArgumentException("invalid sub-domain parameters. xstart and ystart should be 1");
	}
	
	
	/**
     * Compute translating velocity of each point in cylindrical coordinate
     * on a sphere given the velocity at the center of the coordinate.
     *
     * @param	u0		zonal component of velocity at cylindrical coordinate center (m s^-1)
     * @param	v0		meridional component of velocity at cylindrical coordinate center (m s^-1)
     * @param	lon0	longitude of the cylindrical coordinate center (radian)
     * @param	lat0	latitude of the cylindrical coordinate center (radian)
     * @param	lambda	azimuths in cylindrical coordinate (radian), start from north and counter-clockwise
     * @param	radius	radial angles along the great circle of the Earth (radian)
     * @param	lons	longitudes for cylindrical points (radian)
     * @param	lats	latitudes for cylindrical points (radian)
     *
     * @return	re		translating velocity, [0] is zonal and [1] is meridional
     */
	private static float[][][] cInhomoTranslatingVelocityByAnaly
	(float u0,float v0,float lon0,float lat0,float[] lambda,float[] radius,float[][] lons,float[][] lats){
		int y=radius.length,x=lambda.length;
		
		float[][][] re=new float[2][y][x];
		
		float[][] u=re[0];
		float[][] v=re[1];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			double alpha=lambda[i];
			double beta =radius[j];
			double lon  =lons[j][i];
			double lat  =lats[j][i];
			
			v[j][i]=(float)((cos(lat0)*cos(beta)-sin(lat0)*sin(beta)*cos(alpha))/cos(lat)*v0);
			u[j][i]=(float)(u0*cos(lat)/cos(lat0)-sin(beta)*sin(alpha)*tan(lat)/cos(lon0-lon)*v[j][i]);
		}
		
		return re;
	}
	
	/**
     * Compute translating velocity of each point in cylindrical coordinate analytically
     *
     * @return	re		translating velocity, [0] is zonal and [1] is meridional
     */
	public Variable[] cInhomoTranslatingVelocityByAnaly(){
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		int t=csm.getTCount(),y=csm.getYCount(),x=csm.getXCount();
		
		Variable[] re=new Variable[2];
		re[0]=new Variable("uwhole",false,new Range(t,1,y,x));
		re[1]=new Variable("vwhole",false,new Range(t,1,y,x));
		re[0].setComment("zonal translating velocity of the cylindrical coordinate (m s^-1)");
		re[1].setComment("meridional translating velocity of the cylindrical coordinate (m s^-1)");
		
		float[][][] udata=re[0].getData()[0];
		float[][][] vdata=re[1].getData()[0];
		
		for(int l=0;l<t;l++){
			float[][][] tmp=cInhomoTranslatingVelocityByAnaly(
				csm.getUWhole()[l]        ,csm.getVWhole()[l],
				csm.getOLon()[l]          ,csm.getOLat()[l],
				csm.getXDef().getSamples(),csm.getYDef().getSamples(),
				csm.getLon()[l]           ,csm.getLat()[l]
			);
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				udata[j][i][l]=tmp[0][j][i];
				vdata[j][i][l]=tmp[1][j][i];
			}
		}
		
		return re;
	}
	
	/**
     * Compute translating velocity of each point in cylindrical coordinate
     * using central finite difference of the lat/lon coordinates.
     * 
     * Only used to validate cInhomoTranslatingVelocityByAnaly().
     *
     * @return	re		translating velocity, [0] is zonal and [1] is meridional
     */
	public Variable[] cInhomoTranslatingVelocityByDiff(){
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		int t=csm.getTCount(),y=csm.getYCount(),x=csm.getXCount();
		
		Variable[] re=new Variable[2];
		re[0]=new Variable("uwdiff",false,new Range(t,1,y,x));
		re[1]=new Variable("vwdiff",false,new Range(t,1,y,x));
		re[0].setComment("zonal translating velocity of the cylindrical coordinate (m s^-1)");
		re[1].setComment("meridional translating velocity of the cylindrical coordinate (m s^-1)");
		
		float dt=csm.getDT();
		
		float[][][]   lon=csm.getLon();
		float[][][]   lat=csm.getLat();
		float[][][] udata=re[0].getData()[0];
		float[][][] vdata=re[1].getData()[0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			for(int l=1;l<t-1;l++){
				udata[j][i][l]=(float)((lon[l+1][j][i]-lon[l-1][j][i])*REarth*cos(lat[l][j][i])/(dt*2));
				vdata[j][i][l]=(float)((lat[l+1][j][i]-lat[l-1][j][i])*REarth/(dt*2));
			}
			
			udata[j][i][0]=(float)((lon[1][j][i]-lon[0][j][i])*REarth*cos(lat[0][j][i])/dt);
			vdata[j][i][0]=(float)((lat[1][j][i]-lat[0][j][i])*REarth/dt);
			
			udata[j][i][t-1]=(float)((lon[t-1][j][i]-lon[t-2][j][i])*REarth*cos(lat[t-2][j][i])/dt);
			vdata[j][i][t-1]=(float)((lat[t-1][j][i]-lat[t-2][j][i])*REarth/dt);
		}
		
		return re;
	}
	
	public Variable[] cHomoTranslatingVelocity(){
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		int t=csm.getTCount(),y=csm.getYCount(),x=csm.getXCount();
		
		Variable[] re=new Variable[2];
		re[0]=new Variable("uconst",false,new Range(t,1,y,x));
		re[1]=new Variable("vconst",false,new Range(t,1,y,x));
		re[0].setComment("zonal translating velocity of the cylindrical coordinate (m s^-1)");
		re[1].setComment("meridional translating velocity of the cylindrical coordinate (m s^-1)");
		
		float[] cu=csm.getUWhole();
		float[] cv=csm.getVWhole();
		
		float[][][] udata=re[0].getData()[0];
		float[][][] vdata=re[1].getData()[0];
		
		for(int l=0;l<t;l++)
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			udata[j][i][l]=cu[l];
			vdata[j][i][l]=cv[l];
		}
		
		return re;
	}
}
