/**
 * @(#)SurfacePlotModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;


/**
 * template for surfacePlot
 */
public final class SurfacePlotModel{
	//
	public enum PlotMode{WIREFRAME,NORENDER,SPECTRUM,GRAYSCALE,DUALSHADE}
	
	public PlotMode plotMode=PlotMode.SPECTRUM;
	
	public int calcDivisionX=20;
	public int calcDivisionY=20;
	public int dispDivision =Math.min(calcDivisionY,calcDivisionX);
	
	public boolean isMesh        =true;
	public boolean isBoxed       =true;
	public boolean isScaleBox    =true;
	public boolean isDisplayZ    =true;
	public boolean isDisplayXY   =true;
	public boolean isDisplayGrids=true;
	
	public long phase=0;		// for phase animation
	
	public float xmin=-1;
	public float xmax= 1;
	public float ymin=-1;
	public float ymax= 1;
	public float zmin=-1;
	public float zmax= 1;
	
	public float[][] data=null;
	
	public String xlabel="X";
	public String ylabel="Y";
	public String zlabel="Z";
	
	public Function2D function=null;
	
	
	/**
	 * return the mapping of the 2D-function
	 */
	public float calculateZ(float x,float y){ return function.map(x,y);}
	
	public float getXIncrement(){ return (xmax-xmin)/(float)calcDivisionX;}
	
	public float getYIncrement(){ return (ymax-ymin)/(float)calcDivisionY;}
}
