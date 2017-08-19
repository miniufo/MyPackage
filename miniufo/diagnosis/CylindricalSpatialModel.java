/**
 * @(#)CylindricalSpacialModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.SpatialCoordinate;
import miniufo.descriptor.TemporalCoordinate;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.toRadians;


/**
 * Spatial model in Cylindrical coordinates
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CylindricalSpatialModel extends SpatialModel{
	//
	private float[]   olon=null;	// original x-coordinate (radian)
	private float[]   olat=null;	// original y-coordinate (radian)
	
	private float[]   bsin=null;	// sin(beta)
	private float[]   bcos=null;	// cos(beta)
	private float[] uwhole=null;	// velocity in x-direction while associated with the model as a whole
	private float[] vwhole=null;	// velocity in y-direction while associated with the model as a whole
	
	private float[][][] eta=null;	// angle between radial direction and local meridian (radian)
	private float[][][] lon=null;	// (radian)
	private float[][][] lat=null;	// (radian)
	
	
	/**
     * constructor
     *
     * @param	cd		csm descriptor
     */
	public CylindricalSpatialModel(CsmDescriptor cd){ _CylindricalSpatialModel(cd);}
	
	public CylindricalSpatialModel(CsmDescriptor cd,LevelType levType){
		this.levType=levType;
		_CylindricalSpatialModel(cd);
	}
	
	/**
     * constructor
     *
     * @param	olon	central longitudes (degree)
     * @param	olat	central latitudes (degree)
     * @param	yinc	increment in radial direction (degree)
     * @param	ycount	ring count in cylindrical coordinates
     * @param	xcount	grid count in each ring
     */
	public CylindricalSpatialModel(float[] olons,float[] olats,float[] levs,float yinc,int ycount,int xcount){
		int tcount=olons.length;
		
		if(olats.length!=tcount)
			throw new IllegalArgumentException("lengths of lons and lats are not the same");
		
		olon=olons.clone();	for(int l=0;l<tcount;l++) olon[l]=(float)toRadians(olons[l]);
		olat=olats.clone();	for(int l=0;l<tcount;l++) olat[l]=(float)toRadians(olats[l]);
		
		yinc=(float)toRadians(yinc);
		
		boolean same=true;
		for(int l=1;l<tcount;l++)
		if(olon[0]!=olon[l]||olat[0]!=olat[l]){ same=false; break;}
		
		if(!same){
			lon=new float[tcount][][];
			lat=new float[tcount][][];
			eta=new float[tcount][][];
			
			for(int l=0;l<tcount;l++){
				float[][][] re=cLatLons(olon[l],olat[l],yinc,ycount,xcount);
				
				lon[l]=re[0];
				lat[l]=re[1];
				eta[l]=re[2];
			}
			
		}else{
			float[][][] re=cLatLons(olon[0],olat[0],yinc,ycount,xcount);
			
			lon=new float[tcount][ycount][];
			lat=new float[tcount][ycount][];
			eta=new float[tcount][ycount][];
			
			for(int l=0;l<tcount;l++)
			for(int j=0;j<ycount;j++){
				lon[l][j]=re[0][j].clone();
				lat[l][j]=re[1][j].clone();
				eta[l][j]=re[2][j].clone();
			}
		}
		
		float[] xdf=new float[xcount];
		float[] ydf=new float[ycount];
		
		double xinc=2*PI/xcount;
		for(int i=0;i<xcount;i++) xdf[i]=(float)(i*xinc);
		for(int j=0;j<ycount;j++) ydf[j]=(float)(j*yinc);
		
		MDate[] tdf=new MDate[tcount];
		for(int l=0;l<tcount;l++) tdf[l]=new MDate();
		
		float[] zlevs=levs.clone();
		
		if(levType==LevelType.PRESSURE) for(int k=0;k<zlevs.length;k++) zlevs[k]*=100;
		
		tdef=new TemporalCoordinate("tdef",tdf);
		zdef=new  SpatialCoordinate("zdef",zlevs);
		ydef=new  SpatialCoordinate("ydef",ydf);
		xdef=new  SpatialCoordinate("xdef",xdf);
		
		bsin=new float[ycount];	dxs=new float[ycount];
		bcos=new float[ycount];
		
		dt=0;	dy=EARTH_RADIUS*yinc;
		
		for(int j=0;j<ycount;j++)
		dxs[j]=(float)(EARTH_RADIUS*sin(ydef.getSamples()[j])*xinc);
		
		for(int j=0;j<ycount;j++){
			bsin[j]=(float)sin(ydef.getSamples()[j]);
			bcos[j]=(float)cos(ydef.getSamples()[j]);
		}
		
		/*** calculate velocity of the model which moved as a whole ***/
		uwhole=new float[tcount];	vwhole=new float[tcount];
		if(tcount==1){ uwhole[0]=0;	vwhole[0]=0;}
		else{
			for(int l=1;l<tcount-1;l++){
				uwhole[l]=(float)((olon[l+1]-olon[l-1])*EARTH_RADIUS*cos((olat[l+1]+olat[l-1])/2f)/(dt*2f));
				vwhole[l]=(float)((olat[l+1]-olat[l-1])*EARTH_RADIUS/(dt*2f));
			}
			
			uwhole[0]=(float)((olon[1]-olon[0])*EARTH_RADIUS*cos((olat[1]+olat[0])/2f)/dt);
			vwhole[0]=(float)((olat[1]-olat[0])*EARTH_RADIUS/dt);
			
			uwhole[tcount-1]=(float)((olon[tcount-1]-olon[tcount-2])*EARTH_RADIUS*cos((olat[tcount-1]+olat[tcount-2])/2f)/dt);
			vwhole[tcount-1]=(float)((olat[tcount-1]-olat[tcount-2])*EARTH_RADIUS/dt);
		}
	}
	
	
	/*** getor and setor ***/
	public float[] getBSin(){ return bsin;}
	
	public float[] getBCos(){ return bcos;}
	
	public float[] getOLon(){ return olon;}
	
	public float[] getOLat(){ return olat;}
	
	public float[] getUWhole(){ return uwhole;}
	
	public float[] getVWhole(){ return vwhole;}
	
	public float[][][] getEta(){ return eta  ;}
	
	public float[][][] getLon(){ return lon  ;}
	
	public float[][][] getLat(){ return lat  ;}
	
	public void setUWhole(float[] u){
		if(uwhole.length==u.length) uwhole=u;
		else throw new IllegalArgumentException("lengths are not the same");
	}
	
	public void setVWhole(float[] v){
		if(vwhole.length==v.length) vwhole=v;
		else throw new IllegalArgumentException("lengths are not the same");
	}
	
	
	/**
     * whether the given point is in the area
     *
     * @param	lon		longitude (degree) of the given point
     * @param	lat		latitude  (degree) of the given point
     * @param	t		time of the model
     *
     * @return	true or false
     */
	public boolean isInArea(float lon,float lat,int t){
		lon=(float)toRadians(lon);
		lat=(float)toRadians(lat);
		
		if(cSphericalDistanceByRadian(lon,lat,olon[t-1],olat[t-1])<=dy*(ydef.length()-1))
			return true;
		
		return false;
	}
	
	/**
     * whether the given point is in the area
     *
     * @param	lon		longitude (degree) of the given point
     * @param	lat		latitude  (degree) of the given point
     *
     * @return	true or false
     */
	public boolean isInArea(float lon,float lat){ return isInArea(lon,lat,0);}
	
	
	/**
     * used to print out
     */
	public String toString(){
		StringBuffer sb=new StringBuffer();
		
		int tcount=tdef.length(),zcount=zdef.length();
		int ycount=ydef.length(),xcount=xdef.length();
		
		sb.append("\nCylindricalSpacialModel Spacial Model:");
		
		sb.append("\ntdef ("+tcount+"):\n  from: ");
		sb.append(tdef.getSamples()[0]);		sb.append("\n  to  : ");
		sb.append(tdef.getSamples()[tcount-1]);
		
		sb.append("\nzdef ("+zcount+"):\n  from: ");
		sb.append(zdef.getSamples()[0]/100);	sb.append("\n  to  : ");
		sb.append(zdef.getSamples()[zcount-1]/100);
		
		sb.append("\nydef count: ");	sb.append(ycount);		sb.append("      dy:  ");
		sb.append((float)(ydef.getIncrements()[0]*180/PI));	sb.append("\tdegree");
		
		sb.append("\nxdef count: ");	sb.append(xcount);		sb.append("      dx:  ");
		sb.append((float)(xdef.getIncrements()[0]*180/PI));	sb.append("\tdegree");
		
		return sb.toString();
	}
	
	
	/*** helper methods***/
	private void _CylindricalSpatialModel(CsmDescriptor cd){
		tdef=(TemporalCoordinate)(cd.getTDef().clone());
		zdef= (SpatialCoordinate)(cd.getZDef().clone());
		ydef= (SpatialCoordinate)(cd.getYDef().clone());
		xdef= (SpatialCoordinate)(cd.getXDef().clone());
		
		int tcount=tdef.length(),zcount=zdef.length();
		int ycount=ydef.length(),xcount=xdef.length();
		
		olon=new float[tcount];
		olat=new float[tcount];
		
		lon=new float[tcount][ycount][xcount];
		lat=new float[tcount][ycount][xcount];
		eta=new float[tcount][ycount][xcount];
		
		  bsin=new float[ycount];	dxs =new float[ycount];
		  bcos=new float[ycount];
		uwhole=new float[tcount]; 	vwhole=new float[tcount];
		
		/*** process data ***/
		for(int l=0;l<tcount;l++){
			olon[l]=(float)toRadians(cd.getOLon()[l]);
			olat[l]=(float)toRadians(cd.getOLat()[l]);
			
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++){
				lon[l][j][i]=(float)toRadians(cd.getLon()[l][j][i]);
				lat[l][j][i]=(float)toRadians(cd.getLat()[l][j][i]);
				eta[l][j][i]=(float)toRadians(cd.getEta()[l][j][i]);
			}
		}
		
		if(levType==LevelType.PRESSURE)
		for(int k=0;k<zcount;k++) zdef.getSamples()[k]*=100;	// to Pa
		for(int j=0;j<ycount;j++) ydef.getSamples()[j]*=PI/180;	// to radian
		for(int i=0;i<xcount;i++) xdef.getSamples()[i]*=PI/180;	// to radian
		
		if(levType==LevelType.PRESSURE) dz=cd.getDZDef()[0]*100;// to Pa
		else dz=cd.getDZDef()[0];
		
		dt=cd.getDTDef()[0];
		dy=(float)(EARTH_RADIUS*PI*cd.getDYDef()[0]/180.0);
		
		for(int j=0;j<ycount;j++)
		dxs[j]=(float)(EARTH_RADIUS*sin(ydef.getSamples()[j])*PI*cd.getDXDef()[0]/180.0);
		
		/*** calculate velocity of the model which moved as a whole ***/
		if(tcount==1){ uwhole[0]=0;	vwhole[0]=0;}
		else{
			for(int l=1;l<tcount-1;l++){
				uwhole[l]=(float)((olon[l+1]-olon[l-1])*EARTH_RADIUS*cos((olat[l+1]+olat[l-1])/2f)/(dt*2f));
				vwhole[l]=(float)((olat[l+1]-olat[l-1])*EARTH_RADIUS/(dt*2f));
			}
			
			uwhole[0]=(float)((olon[1]-olon[0])*EARTH_RADIUS*cos((olat[1]+olat[0])/2f)/dt);
			vwhole[0]=(float)((olat[1]-olat[0])*EARTH_RADIUS/dt);
			
			uwhole[tcount-1]=(float)((olon[tcount-1]-olon[tcount-2])*EARTH_RADIUS*cos((olat[tcount-1]+olat[tcount-2])/2f)/dt);
			vwhole[tcount-1]=(float)((olat[tcount-1]-olat[tcount-2])*EARTH_RADIUS/dt);
		}
		
		for(int j=0;j<ycount;j++){
			bsin[j]=(float)sin(ydef.getSamples()[j]);
			bcos[j]=(float)cos(ydef.getSamples()[j]);
		}
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
