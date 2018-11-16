/**
 * @(#)SpatialModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.asin;
import static java.lang.Math.acos;
import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;
import static java.lang.Math.toDegrees;
import miniufo.descriptor.SpatialCoordinate;
import miniufo.descriptor.TemporalCoordinate;
import miniufo.mathsphysics.MathsPhysicsUtil;


/**
 * Abstract class for spatial model.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class SpatialModel{
	//
	protected boolean periodicX=false;
	protected boolean periodicY=false;
	
	protected TemporalCoordinate tdef=null;
	protected  SpatialCoordinate zdef=null;
	protected  SpatialCoordinate ydef=null;
	protected  SpatialCoordinate xdef=null;
	
	protected float dt;	// unit: s, length is usually 1
	protected float dz;	// unit: depend on LevelType
	protected float dy;	// unit: m
	protected float dx;	// unit: m
	
	protected LevelType levType=LevelType.PRESSURE;
	
	protected float[] dxs=null;	// dx scaled by lcos[j]
	
	public static final float gEarth    =9.80665f;	// gravitational acceleration
	public static final float omegaEarth=7.292e-5f;	// rotating speed of earth
	public static final float REarth    =6371200f;	// radius of the earth (meter), consistent with GrADS
	
	public enum LevelType{GEOMETRIC,PRESSURE,ISENTROPIC,SIGMA,ISOPYCNAL};
	
	
	/*** getor and setor ***/
	public int getTCount(){ return tdef.length();}
	
	public int getZCount(){ return zdef.length();}
	
	public int getYCount(){ return ydef.length();}
	
	public int getXCount(){ return xdef.length();}
	
	public boolean isLinearModel(){
		return
		tdef.isLinear()&&zdef.isLinear()&&
		ydef.isLinear()&&xdef.isLinear();
	}
	
	public boolean isPeriodicX(){ return periodicX;}
	
	public boolean isPeriodicY(){ return periodicY;}
	
	public float getDT(){ return dt;}
	
	public float getDY(){ return dy;}
	
	public float getDZ(){ return dz;}
	
	public float getDX(){ return dx;}
	
	public float[] getDXs(){ return dxs;}
	
	public LevelType getLevelType(){ return levType;}
	
	public TemporalCoordinate getTDef(){ return tdef;}
	
	public  SpatialCoordinate getZDef(){ return zdef;}
	
	public  SpatialCoordinate getYDef(){ return ydef;}
	
	public  SpatialCoordinate getXDef(){ return xdef;}
	
	
	/**
     * whether the model is spacially the same as the given Variable
     *
     * @param	v	a given Variable
     */
	public boolean isDimLike(Variable v){
		if(tdef.length()!=v.getTCount()) return false;
		if(zdef.length()!=v.getZCount()) return false;
		if(ydef.length()!=v.getYCount()) return false;
		if(xdef.length()!=v.getXCount()) return false;
		
		return true;
	}
	
	/**
     * whether the model is planarily the same as the given Variable
     *
     * @param	v	a given Variable
     */
	public boolean isAreaLike(Variable v){
		//if(tdef.length()!=v.getTCount()) return false;
		if(ydef.length()!=v.getYCount()) return false;
		if(xdef.length()!=v.getXCount()) return false;
		
		return true;
	}
	
    
    /**
     * calculate spherical (great-circle) distance between two given points
     *
     * @param	lon1	longitude of the first point (degree)
     * @param	lat1	latitude  of the first point (degree)
     * @param	lon2	longitude of the second point (degree)
     * @param	lat2	latitude  of the second point (degree)
     *
     * @return	spherical distance (m)
     */
    public static float cSphericalDistanceByDegree(float lon1,float lat1,float lon2,float lat2){
    	return cSphericalDistanceByRadian(
    		(float)toRadians(lon1),(float)toRadians(lat1),(float)toRadians(lon2),(float)toRadians(lat2)
    	);
    }
    
    /**
     * calculate spherical (great-circle) distance between two given points
     *
     * @param	lon1	longitude of the first point (radian)
     * @param	lat1	latitude  of the first point (radian)
     * @param	lon2	longitude of the second point (radian)
     * @param	lat2	latitude  of the second point (radian)
     *
     * @return	spherical distance
     */
    public static float cSphericalDistanceByRadian(float lon1,float lat1,float lon2,float lat2){
    	if(lon1==lon2&&lat1==lat2)
    		return 0;
    	else{
    		// Algorithm from Matlab, using Haversine formula
    		// faster than R * acos( cos(lat1)*cos(lat2)*cos(lon1-lon2)+sin(lat1)*sin(lat2) )
    		double a=((1.0-cos(lat2-lat1))+cos(lat1)*cos(lat2)*(1.0-cos(lon2-lon1)))/2.0;
    		
    		if(a<0) a=0;
    		if(a>1) a=1;
    		
    		return (float)(REarth*2.0*Math.atan2(sqrt(a),sqrt(1-a)));
    	}
    }
    
	
	/**
	 * compute the Area on a sphere bounded by two latitudes and longitudes
	 * following Matlab command areaquad.
	 * 
	 * @param	lon1	longitude of the southwest corner (Radian)
	 * @param	lat1	latitudes of the southwest corner (Radian)
	 * @param	lon2	longitude of the northeast corner (Radian)
	 * @param	lat2	latitudes of the northeast corner (Radian)
	 */
	public static float cAreaQuadByDegree(float lon1,float lat1,float lon2,float lat2){
		return (float)(MathsPhysicsUtil.cAreaQuadByRadian(
			Math.toRadians(lon1),Math.toRadians(lat1),
			Math.toRadians(lon2),Math.toRadians(lat2)
		)*4.0*Math.PI*REarth*REarth);
	}
	
	public static float cAreaQuadByRadian(float lon1,float lat1,float lon2,float lat2){
		return (float)MathsPhysicsUtil.cAreaQuadByRadian((double)lon1,lat1,lon2,lat2);
	}
	
	
    /**
     * Calculate the lat lon of a point given the cylindrical
     * coordinates radius and lamda relevance to a fix point
     *
     * @param	olon	longitude of a fix point (degree)
     * @param	olat	latitude  of a fix point (degree)
     * @param	lamda	counter-clockwise away from the North (degree)
     * @param	r		radius away from the fix point (degree)
     *
     * @return	re		re[0] is longitude, re[1] is latitude and [2] is eta in degrees
     */
    public static float[] cLatLon(double olon,double olat,double lambda,double r){
    	float[] re=new float[3];
    	
    	olon  =toRadians(olon);
    	olat  =toRadians(olat);
    	lambda=toRadians(lambda);
    	r     =toRadians(r);
    	
		double tmp=asin(sin(olat)*cos(r)+cos(olat)*sin(r)*cos(lambda));
		
		re[1]=(float)toDegrees(tmp);
		
		double dlambda=asin(sin(r)*sin(lambda)/cos(tmp));
		
		tmp=olon-dlambda;
		
		if(tmp<0)     tmp+=PI*2;
		if(tmp>=PI*2) tmp-=PI*2;
		
		re[0]=(float)toDegrees(tmp);
		
		if(lambda<=PI) re[2]=(float)toDegrees(PI-acos(-cos(lambda)*cos(dlambda)+sin(lambda)*sin(dlambda)*sin(olat)));
		else           re[2]=(float)toDegrees(PI+acos(-cos(lambda)*cos(dlambda)+sin(lambda)*sin(dlambda)*sin(olat)));
    	
    	return re;
    }
    
	/**
     * calculate the lat/lon of cylindrical grids
     *
     * @param	lon0	central longitude of cylindrical coordinate (in radian)
     * @param	lat0	central latitude of cylindrical coordinate (in radian)
     * @param	yinc	delta angle in radial direction (in radian)
     * @param	ycount	number of rings in cylindrical coordinate
     * @param	xcount	number of grids in one ring
     *
     * @return	lats/lons, [0] is lons, [1] is lats and [2] is etas in radians
     */
	public static float[][][] cLatLons(double olon,double olat,double yinc,int ycount,int xcount){
		double xinc=2*PI/xcount;
		
		float[][][] re=new float[3][ycount][xcount];
		
		for(int j=0;j<ycount;j++){
			int half=xcount/2;
			
			/*** when alpha is less than or equals PI/2 ***/
			for(int i=0;i<=half;i++){
				re[1][j][i]=(float)asin(
					sin(olat)*cos(yinc*j)+
					cos(olat)*sin(yinc*j)*cos(xinc*i)
				);
				
				double dlambda=asin(sin(yinc*j)*sin(xinc*i)/cos(re[1][j][i]));
				
				re[0][j][i]=(float)(olon-dlambda);
				
				if(re[0][j][i]<0)     re[0][j][i]+=(float)PI*2;
				if(re[0][j][i]>=PI*2) re[0][j][i]-=(float)PI*2;
				
				re[2][j][i]=(float)(PI-
					acos(-cos(xinc*i)*cos(dlambda)+sin(xinc*i)*sin(dlambda)*sin(olat))
				);
			}
			
			/*** when alpha is larger than PI/2 ***/
			for(int i=half+1;i<xcount;i++){
				re[1][j][i]=(float)asin(
					sin(olat)*cos(yinc*j)+
					cos(olat)*sin(yinc*j)*cos(xinc*i)
				);
				
				double dlambda=asin(sin(yinc*j)*sin(xinc*i)/cos(re[1][j][i]));
				
				re[0][j][i]=(float)(olon-dlambda);
				
				if(re[0][j][i]<0)     re[0][j][i]+=(float)PI*2;
				if(re[0][j][i]>=PI*2) re[0][j][i]-=(float)PI*2;
				
				re[2][j][i]=(float)(PI+
					acos(-cos(xinc*i)*cos(dlambda)+sin(xinc*i)*sin(dlambda)*sin(olat))
				);
			}
		}
		
		return re;
	}
	
    
	/** test
	public static void main(String[] args){
		float[] pos=cLatLon(120,0,1,90);
		System.out.println(pos[0]+"\t"+pos[1]);
	}*/
}
