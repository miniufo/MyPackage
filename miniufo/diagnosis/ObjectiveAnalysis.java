/**
 * @(#)ObjectiveAnalysis.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import java.util.List;
import java.util.ArrayList;
import static java.lang.Math.abs;
import static miniufo.diagnosis.SpatialModel.cSphericalDistanceByDegree;


/**
 * objective analysis
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class ObjectiveAnalysis{
	//
	private int ycount  =0;
	private int xcount  =0;
	private int stncount=0;
	
	private float undef =0;
	
	private float[]		buf=null;	// length is stncount
	private float[][]   stn=null;	// a 3*stncount array, [0] is lons, [1] is lats and [2] is values
	private float[][][] grd=null;	// a 3*ycount*xcount array, [0] is lons, [1] is lats and [2] is values
	
	private List<Integer> stntag=new ArrayList<Integer>();	// station tags in the radius
	private List< Float > stndis=new ArrayList< Float >();	// distances of the stations
	private List< Float > stnwei=new ArrayList< Float >();	// weight coefficient of the stations
	
	private List< int[] > grdtag=new ArrayList< int[] >();	// grid tags in the radius
	private List< Float > grddis=new ArrayList< Float >();	// distances of the grids
	private List< Float > grdwei=new ArrayList< Float >();	// weight coefficient of the grids
	
	
	/**
     * constructor
     *
     * @param	stn		stations data, stn[0] is lons, stn[1] lats and stn[2] data
     */
	public ObjectiveAnalysis(float[][][] grids,float[][] stations,float undef){
		if(stations.length!=3)
		throw new IllegalArgumentException("stations data should be an array as stn[3][stncount]");
		
		if(grids.length!=3)
		throw new IllegalArgumentException("grids data should be an array as grd[3][ycount][xcount]");
		
		  ycount=grids[0].length;
		  xcount=grids[0][0].length;
		stncount=stations[0].length;
		
		stn=stations;	grd=grids;	this.undef=undef;
		
		buf=new float[stncount];
	}
	
	
	/*** getor and setor ***/
	public int getStationCount(){ return stncount;}
	
	public float[][] getStations(){ return stn;}
	
	public float[][][] getGrids(){ return grd;}
	
	
	/**
     * successive correction method
     *
     * @param	radii	radii of scanning (m)
     */
	public void SCMByCressman(float... radii){
		/*** first guess, initialized by mean of the stations ***/
		initialization();
		
		for(float rad:radii){
			float radInGridUnit=
			(float)Math.toDegrees(rad/SpatialModel.EARTH_RADIUS)/(grd[1][1][0]-grd[1][0][0]);
			System.out.println("start analyzing for rad = "+Math.round(radInGridUnit)+" grid space");
			
			/*** compute the observation increment ***/
			for(int i=0;i<stncount;i++){
				findGrids(stn[0][i],stn[1][i],rad);
				cGridWeightCoefficientsByCressman(rad);
				
				buf[i]=interpGrid();
			}
			
			float rad2=rad*rad;
			/*** correct the grid ***/
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++){
				findStations(grd[0][j][i],grd[1][j][i],rad);
				cStationWeightCoefficientsByCressman(rad);
				
				float dens=stntag.size()/rad2;
				
				if(dens>=1E-11f&&grd[2][j][i]!=undef) grd[2][j][i]+=interpStation();
				else if(rad==radii[0]) grd[2][j][i]=undef;
			}
		}
	}
	
	/**
     * successive correction method (unfinished)
     *
     * @param	B		argument B
     * @param	C		argument C
     */
	public void SCMByBarnes(float... radii){
		/*** first guess, initialized by mean of the stations ***/
		initialization();
		
		for(float rad:radii){
			float rad2=rad*rad;
			
			/*** compute the observation increment ***/
			for(int i=0;i<stncount;i++){
				findGrids(stn[0][i],stn[1][i],rad);
				cGridWeightCoefficientsByBarnes(rad);
				
				buf[i]=interpGrid();
			}
			
			/*** correct the grid ***/
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++){
				findStations(grd[0][j][i],grd[1][j][i],rad);
				cStationWeightCoefficientsByBarnes(rad);
				
				float dens=stntag.size()/rad2;
				
				if(dens<1E-11f) grd[2][j][i]=undef;
				else grd[2][j][i]+=interpStation();
			}
		}
	}
	
	
	/**
     * initial field of grids by mean of stations
     *
     * @param	grid	grid data, grid[0][j][i] is lons[j][i],
     * 					grid[1][j][i] lats[j][i] and grid[2][j][i] data[j][i]
     */
	private void initialization(){
		int count=0;
		
		float ave=0;
		for(int i=0;i<stncount;i++)
		if(stn[2][i]!=undef){ ave+=stn[2][i]; count++;}
		if(count!=0) ave/=count;
		
		for(int j=0;j<ycount;j++)
		for(int i=0;i<xcount;i++) grd[2][j][i]=ave;
	}
	
	
	/**
     * find stations in the radius around a grid
     *
     * @param	lon		longitude of the station (degree)
     * @param	lat		latitude  of the station (degree)
     * @param	rad		radius of the analysis (m)
     */
	private void findStations(float lon,float lat,float rad){
		stntag.clear();	stndis.clear();
		
		// scanning the stations in the radius
		float radDegree=(float)Math.toDegrees(rad/SpatialModel.EARTH_RADIUS);
		for(int i=0;i<stncount;i++)
		if(abs(lon-stn[0][i])<radDegree&&abs(lat-stn[1][i])<radDegree){
			float sdis=cSphericalDistanceByDegree(stn[0][i],stn[1][i],lon,lat);
			
			if(sdis<=rad){ stntag.add(i);	stndis.add(sdis);}
		}
	}
	
	/**
     * find grid in the radius around a station
     *
     * @param	lon		longitude of the grid (degree)
     * @param	lat		latitude  of the grid (degree)
     * @param	rad		radius of the analysis (m)
     */
	private void findGrids(float lon,float lat,float rad){
		grdtag.clear();	grddis.clear();
		
		// scanning the stations in the radius
		float radDegree=(float)Math.toDegrees(rad/SpatialModel.EARTH_RADIUS);
		for(int j=0;j<ycount;j++)
		for(int i=0;i<xcount;i++)
		if(abs(lon-grd[0][j][i])<radDegree&&abs(lat-grd[1][j][i])<radDegree){
			float sdis=cSphericalDistanceByDegree(grd[0][j][i],grd[1][j][i],lon,lat);
			if(sdis<=rad){ grdtag.add(new int[]{j,i});	grddis.add(sdis);}
		}
	}
	
	
	/**
     * compute the weight coefficients
     * 
     * @param	rad		radius of the analysis (m)
     */
	private void cStationWeightCoefficientsByCressman(float rad){
		stnwei.clear();	float rad2=rad*rad;
		
		for(int i=0,I=stntag.size();i<I;i++){
			float dis2=stndis.get(i);	dis2*=dis2;
			
			stnwei.add((rad2-dis2)/(rad2+dis2));
		}
	}
	
	private void cStationWeightCoefficientsByBarnes(float rad){
		stnwei.clear();	float rad2=rad*rad;
		
		for(int i=0,I=stntag.size();i<I;i++){
			float dis2=stndis.get(i);	dis2*=dis2;
			
			stnwei.add((float)Math.exp(-dis2/(2*rad2)));
		}
	}
	
	private void cGridWeightCoefficientsByCressman(float rad){
		grdwei.clear();	float rad2=rad*rad;
		
		for(int i=0,I=grdtag.size();i<I;i++){
			float dis2=grddis.get(i);	dis2*=dis2;
			
			grdwei.add((rad2-dis2)/(rad2+dis2));
		}
	}
	
	private void cGridWeightCoefficientsByBarnes(float rad){
		grdwei.clear();	float rad2=rad*rad;
		
		for(int i=0,I=grdtag.size();i<I;i++){
			float dis2=grddis.get(i);	dis2*=dis2;
			
			grdwei.add((float)Math.exp(-dis2/(2*rad2)));
		}
	}
	
	
	/**
     * interpolate one grid
     * 
     * @param	f	value of the station interpolated by grids
     */
	private float interpStation(){
		float sum=0,wsum=0;
		
		for(int i=0,I=stntag.size();i<I;i++){
			float w=stnwei.get(i);
			int tag=stntag.get(i);
			
			sum+=w*(stn[2][tag]-buf[tag]);
			wsum+=w;
		}
		
		return sum/wsum;
	}
	
	/**
     * interpolate one station
     */
	private float interpGrid(){
		float sum=0,wsum=0;
		
		for(int i=0,I=grdtag.size();i<I;i++){
			float w=grdwei.get(i);
			int[] t=grdtag.get(i);
			
			sum+=w*grd[2][t[0]][t[1]];
			wsum+=w;
		}
		
		return sum/wsum;
	}
	
	
	/** test 
	public static void main(String[] args){
		try{
			float[][] stn={
				{180,190,180,190},
				{ 30, 30, 20, 20},
				{  1,  4,  4,  9}
			};
			
			float[][][] grid={
				{{181}},
				{{ 31}},
				{{  1}}
			};
			
			ObjectiveAnalysis oa=new ObjectiveAnalysis(stn);
			
		}catch(Exception e){ e.printStackTrace();}
	}*/
}
