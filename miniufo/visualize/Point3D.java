/**
 * @(#)Point3D.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import java.awt.Point;


/**
 * The class <code>SurfaceVertex</code> represents a vertex in 3D space.
 */
public final class Point3D{
	//
	public float x,y,z;	// coordinates
	
	private int    project_index;
	private Point  projection;
	
	private static int   master_project_index=0;	// over 4 billion times to reset
	private static float zmin,zmax;
	private static float zfactor;
	
	private static Projector projector;
	
	
	/**
	 * The constructor of <code>SurfaceVertex</code>.
	 * The x and y coordinated must be in normalized form, i.e: in the range -10 .. +10.
	 *
	 * @param ix the x coordinate
	 * @param iy the y coordinate
	 * @param iz the z coordinate
	 */
	Point3D(float ix,float iy,float iz){
		x=ix;	y=iy;	z=iz;
		project_index=master_project_index-1;
	}
	
	
	/**
	 * Determines whether this vertex is invalid, i.e has invalid coordinates value.
	 *
	 * @return <code>true</code> if this vertex is invalid
	 */
	public final boolean isInvalid(){ return Float.isNaN(z);}
	
	/**
	 * Gets the 2D projection of the vertex.
	 *
	 * @return the 2D projection
	 */
	public final Point projection(){
		if(project_index!=master_project_index)
		projection=projector.project(x,y,(z-zmin)*zfactor-10);
		project_index=master_project_index;
		
		return projection;
	}
	
	/**
	 * Transforms coordinate values to fit the scaling factor of the
	 * projector. This routine is only used for transforming center of projection
	 * in Surface Plotter.
	 */
	public final void transform(){
		x=x/projector.getXScaling();
		y=y/projector.getYScaling();
		z=(zmax-zmin)*(z/projector.getZScaling()+10)/20+zmin;
	}
	
	
	/**
	 * Invalidates all vertices. This will force the projector
	 * to recalculate vertex projection.
	 */
	public static void invalidate(){ master_project_index++;}
	
	/**
	 * Sets the projector to project this vertex.
	 *
	 * @param projector the projector
	 */
	public static void setProjector(Projector projector){ Point3D.projector=projector;}
	
	/**
	 * Sets the minimum and maximum value of z range.
	 * This values is used to compute a factor to normalized
	 * z values into the range -10 .. +10.
	 *
	 * @param zmin the minimum z
	 * @param zmax the maximum z 
	 */
	public static void setZRange(float zmin,float zmax){
		Point3D.zmin=zmin;
		Point3D.zmax=zmax;
		
		zfactor=20/(zmax-zmin);
	}
}
