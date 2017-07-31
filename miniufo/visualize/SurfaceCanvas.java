/**
 * @(#)SurfaceCanvas.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import java.awt.Font;
import java.awt.Color;
import java.awt.Event;
import java.awt.Image;
import java.awt.Label;
import java.awt.Point;
import java.awt.Canvas;
import java.awt.Graphics;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.PrintGraphics;
import static miniufo.visualize.SurfacePlotModel.PlotMode;


/**
 * The class <code>SurfaceCanvas</code> is responsible
 * for the generation of surface images and user mouse events handling. 
 *
 * @author  Yanto Suryono
 */
public final class SurfaceCanvas extends Canvas{
	//
	private static final long serialVersionUID = 1082939544452440193L;
	
	private int cvBufWidth, cvBufHeight;	// canvas size
	private int printwidth, printheight;	// print size
	private int calc_divisionsX;
	private int calc_divisionsY;
	
	private boolean image_drawn;			// image drawn flag
	private boolean data_available;			// data availability flag
	private boolean printing;				// printing flag
	private boolean isBoxed;
	private boolean isMesh;
	private boolean isScaleBox;
	private boolean isDisplayXY;
	private boolean isDisplayZ;
	private boolean isDisplayGrids;
	
	private float color;					// color of surface
	private float xmin,xmax;
	private float ymin,ymax;
	private float zmin,zmax;
	
	private Point3D cop;					// center of projection
	private Point3D[] vertex;				// vertices array
	private Projector projector;			// the projector
	private Image bufferImage;				// the backing buffer
	private Graphics bufferGraphics;		// the graphics context of backing buffer
	private SurfacePlotModel model;			// the model

	// constants
	private static final int TOP = 0;
	private static final int CENTER = 1;
	
	
	/**
	 * The constructor of <code>SurfaceCanvas</code>
	 *
	 * @param model The model provides data about the surface
	 */
	public SurfaceCanvas(){
		super();
		projector = new Projector();
		projector.setDistance(65);
		projector.set2DScaling(11);
		projector.setRotationAngle(220);
		projector.setElevationAngle(25);
		Point3D.setProjector(projector);
	}
	
	
	public void setModel(SurfacePlotModel model){
		this.model = model;
		isBoxed = model.isBoxed;
		isMesh = model.isMesh;
		isScaleBox = model.isScaleBox;
		isDisplayXY = model.isDisplayXY;
		isDisplayZ = model.isDisplayZ;
		isDisplayGrids = model.isDisplayGrids;
		calc_divisionsX = model.calcDivisionX;
		calc_divisionsY = model.calcDivisionY;
		bufferImage = null;
		bufferGraphics = null;
		image_drawn = false;
		data_available = false;
		printing = false;
		cvBufWidth = cvBufHeight = -1;
		if(bufferImage != null) bufferImage.flush();
		bufferImage = null;
		
		renderSurface(model.data);
	}

	/**
	 * Destroys the internal image. It will force <code>SurfaceCanvas</code>
	 * to regenerate all images when the <code>paint</code> method is called.
	 */
	public void destroyImage(){ image_drawn=false;}

	/**
	 * Sets the x and y ranges of calculated surface vertices.
	 * The ranges will not affect surface appearance. They affect axes
	 * scale appearance.
	 *
	 * @param xmin the minimum x
	 * @param xmax the maximum x
	 * @param ymin the minimum y
	 * @param ymax the maximum y
	 */
	public void setRanges(float xmin,float xmax,float ymin,float ymax){
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
	}
	
	/**
	 * Sets the data availability flag. If this flag is <code>false</code>,
	 * <code>SurfaceCanvas</code> will not generate any surface image, even
	 * if the data is available. But it is the programmer's responsiblity
	 * to set this flag to <code>false</code> when data is not available.
	 *
	 * @param avail the availability flag
	 */
	public void setDataAvailability(boolean avail){ data_available=avail;}
	
	/**
	 * Sets the new vertices array of surface.
	 *
	 * @param vertex the new vertices array
	 * @see   #getValuesArray
	 */
	public void setValuesArray(Point3D[] vertex){ this.vertex=vertex;}
	
	/**
	 * Gets the current vertices array.
	 *
	 * @return current vertices array
	 * @see    #setValuesArray
	 */
	public Point3D[] getValuesArray(){
		if (!data_available) return null;
		
		return vertex;
	}
	
	/**
	 * Gets the current x, y, and z ranges.
	 *
	 * @return array of x,y, and z ranges in order of
	 *         xmin, xmax, ymin, ymax, zmin, zmax
	 */
	public float[] getRanges(){
		float[] ranges=new float[6];
		
		ranges[0]=xmin;
		ranges[1]=xmax;
		ranges[2]=ymin;
		ranges[3]=ymax;
		ranges[4]=zmin;
		ranges[5]=zmax;
		
		return ranges;
	}
	
	
	private int click_x, click_y;    // previous mouse cursor position
	
	/**
	 * <code>mouseDown</code> event handler. Sets internal tracking variables
	 * for dragging operations.
	 *
	 * @param e the event
	 * @param x the x coordinate of cursor
	 * @param y the y coordinate of cursor
	 */
	@Override
	public boolean mouseDown(Event e,int x,int y){
		click_x = x;
		click_y = y;
		return true;
	}
	
	/**
	 * <code>mouseDrag<code> event handler. Tracks dragging operations.
	 * Checks the delay regeneration flag and does proper actions.
	 *
	 * @param e the event
	 * @param x the x coordinate of cursor
	 */
	@Override
	public boolean mouseDrag(Event e,int x,int y){
		if (e.controlDown()){
			projector.set2D_xTranslation(
				projector.get2D_xTranslation() + (x - click_x));
			projector.set2D_yTranslation(
				projector.get2D_yTranslation() + (y - click_y));
			
		}else{
			if (e.shiftDown()){
				float new_value = projector.get2DScaling() + (y - click_y) * 0.5f;
				
				if (new_value > 60.0f) new_value = 60.0f;
				
				if (new_value < 2.0f) new_value = 2.0f;
				
				projector.set2DScaling(new_value);
				
			}else{
				float new_value = projector.getRotationAngle() + (x - click_x);
				
				while (new_value > 360) new_value -= 360;
				
				while (new_value < 0) new_value += 360;
				
				projector.setRotationAngle(new_value);
				
				new_value = projector.getElevationAngle() + (y - click_y);
				
				if (new_value > 90)new_value = 90;
				else if (new_value < 0) new_value = 0;
				
				projector.setElevationAngle(new_value);
			}
		}
		
		image_drawn = false;
		repaint();
		click_x = x;
		click_y = y;
		return true;
	}
	
	private void renderSurface(float[][] data){
		final float xi = model.xmin;
		final float yi = model.ymin;
		final float xx = model.xmax;
		final float yx = model.ymax;

		final float stepx = (xx - xi) / calc_divisionsX;
		final float stepy = (yx - yi) / calc_divisionsY;

		final float xfactor = 20 / (xx - xi);
		final float yfactor = 20 / (yx - yi);

		setRanges(xi, xx, yi, yx);
		setDataAvailability(false);
		destroyImage();
		
		Point3D[] tmpVertices = new Point3D[(calc_divisionsX+1)*(calc_divisionsY+1)];

		float max = Float.NaN;
		float min = Float.NaN;

		float x = xi;
		for (int i=0,k=0;i<=calc_divisionsX;i++){
			float y = yi;
			
			for(int j=0; j<=calc_divisionsY;j++){
				float v = data[j][i];
				
				if(Float.isInfinite(v)) v = Float.NaN;
				
				if(!Float.isNaN(v)){
					if(Float.isNaN(max) || (v > max)) max = v;
					else if(Float.isNaN(min) || (v < min))min = v;
				}
				
				tmpVertices[k++] = new Point3D(
					(x - xi) * xfactor - 10,
					(y - yi) * yfactor - 10, v
				);
				
				y += stepy;
			}
			
			x += stepx;
		}

		setValuesArray(tmpVertices);
		setDataAvailability(true);
		repaint();
	}
	
	public Image getImage(){ return bufferImage;}
	
	/**
	 * Paints surface. Creates surface plot, contour plot, or density plot
	 * based on current vertices array, contour plot flag, and density plot
	 * flag. If no data is available, creates image of base plane and axes.
	 *
	 * @param g the graphics context to paint
	 * @see   #setContour
	 * @see   #setDensity
	 * @see   #setValuesArray
	 * @see   #setDataAvailability
	 */
	public void paint(Graphics g){
		if((getBounds().width<=0)||(getBounds().height<=0)) return;
		
		// Initialize Buffer Image
		if(bufferImage == null
			|| (getBounds().width != cvBufWidth)
			|| (getBounds().height!= cvBufHeight)){
			projector.setProjectionArea(new Rectangle(0, 0, getBounds().width, getBounds().height));
			
			image_drawn = false;
			
			if(bufferImage!=null) bufferImage.flush();
			bufferImage=createImage(getBounds().width,getBounds().height);
			
			if(bufferGraphics!=null) bufferGraphics.dispose();
			bufferGraphics=bufferImage.getGraphics();
			
			cvBufWidth =getBounds().width;
			cvBufHeight=getBounds().height;
		}
		
		if(g instanceof PrintGraphics){
			// modifies variables
			Graphics savedgc=bufferGraphics;
			bufferGraphics = g;
			
			Dimension pagedimension=((PrintGraphics)g).getPrintJob().getPageDimension();
			
			printwidth =pagedimension.width;
			printheight=cvBufHeight*printwidth/cvBufWidth;
			
			if(printheight > pagedimension.height){
				printheight= pagedimension.height;
				printwidth = cvBufWidth*printheight/cvBufHeight;
			}
			
			float savedscalingfactor=projector.get2DScaling();
			projector.setProjectionArea(new Rectangle(0,0,printwidth,printheight));
			projector.set2DScaling(savedscalingfactor*printwidth/cvBufWidth);
			
			bufferGraphics.clipRect(0,0,printwidth,printheight);
			
			// starts printing
			if(!data_available) drawBoxGridsTicksLabels(bufferGraphics);
			else{
				int fontsize=(int)(Math.round(projector.get2DScaling()*0.8));
				bufferGraphics.setFont(new Font("Arial",Font.BOLD,fontsize));
				
				Point3D.invalidate();
				
				// surface plot
				if(model.plotMode==PlotMode.WIREFRAME) plotWireframe();
				else plotSurface();
				
				if(isBoxed) drawBoundingBox();
			}
			
			bufferGraphics.drawRect(0,0,printwidth-1,printheight-1);
			
			// restores variables
			projector.set2DScaling(savedscalingfactor);
			projector.setProjectionArea(new Rectangle(0,0,getBounds().width,getBounds().height));
			
			bufferGraphics = savedgc;
			
		}else{
			if(image_drawn&&(bufferImage!=null)) g.drawImage(bufferImage,0,0,this);
			else{
				if(data_available){
					int fontsize=(int)(Math.round(projector.get2DScaling()*0.8));
					
					if(bufferGraphics!=null)
					bufferGraphics.setFont(new Font("Arial",Font.BOLD,fontsize));
					
					Point3D.invalidate();
					
					// surface plot
					if(model.plotMode==PlotMode.WIREFRAME) plotWireframe();
					else plotSurface();
					
					if(isBoxed) drawBoundingBox();
					
					image_drawn=true;
					
					g.drawImage(bufferImage,0,0,this);
					
				}else{
					g.setColor(Color.lightGray);
					g.fillRect(0,0,getBounds().width,getBounds().height);
					drawBoxGridsTicksLabels(g);
				}
			}
		}
	}
	
	/**
	 * Updates image. Just call the <code>paint</code> method to
	 * avoid flickers.
	 * do not erase, just paint.
	 *
	 * @param g the graphics context to update
	 * @see   #paint
	 */
	public void update(Graphics g){ paint(g);}
	
	
	/********* Private methods begin here *********/
	private int factor_x, factor_y;   // conversion factors
	private int t_x, t_y, t_z;        // determines ticks density
	
	private float color_factor;
	private final int[] poly_x = new int[9];
	private final int[] poly_y = new int[9];
	private Color line_color;
	
	/**
	 * Draws the bounding box of surface.
	 */
	private void drawBoundingBox(){
		Point startingpoint=projector.project(factor_x*10,factor_y*10,10);
		if(bufferGraphics!=null)
		bufferGraphics.setColor(Color.black);
		
		Point projection=null;
		projection=projector.project(-factor_x*10, factor_y*10, 10);
		if(bufferGraphics!=null)
		bufferGraphics.drawLine(startingpoint.x,startingpoint.y,projection.x,projection.y);
		
		projection=projector.project( factor_x*10,-factor_y*10, 10);
		if(bufferGraphics!=null)
		bufferGraphics.drawLine(startingpoint.x,startingpoint.y,projection.x,projection.y);
		
		projection=projector.project( factor_x*10, factor_y*10,-10);
		if(bufferGraphics!=null)
		bufferGraphics.drawLine(startingpoint.x,startingpoint.y,projection.x,projection.y);
	}
	
	/**
	 * Draws the base plane. The base plane is the x-y plane.
	 *
	 * @param g the graphics context to draw.
	 * @param x used to retrieve x coordinates of drawn plane from this method.
	 * @param y used to retrieve y coordinates of drawn plane from this method.
	 */
	private void drawBase(Graphics g,int[] x,int[] y){
		Point projection=null;
		
		projection=projector.project(-10,-10,-10);
		x[0]=projection.x;
		y[0]=projection.y;
		projection=projector.project(-10, 10,-10);
		x[1]=projection.x;
		y[1]=projection.y;
		projection=projector.project( 10, 10,-10);
		x[2]=projection.x;
		y[2]=projection.y;
		projection=projector.project( 10,-10,-10);
		x[3]=projection.x;
		y[3]=projection.y;
		
		x[4]=x[0];
		y[4]=y[0];
		
		if(model!=null){
			if(model.plotMode!=PlotMode.WIREFRAME){
				if(model.plotMode==PlotMode.NORENDER&&g!=null) g.setColor(Color.lightGray);
				else if(g!=null) g.setColor(new Color(192,220,192));
				
				if(g!=null) g.fillPolygon(x,y,4);
			}
			
			if(g!=null) g.setColor(Color.black);
			if(g!=null) g.drawPolygon(x,y,5);
		}
	}
	
	/**
	 * Draws non-surface parts, i.e: bounding box, axis grids, axis ticks,
	 * axis labels, base plane.
	 *
	 * @param g         the graphics context to draw
	 * @param draw_axes if <code>true</code>, only draws base plane and z axis
	 */
	private void drawBoxGridsTicksLabels(Graphics g){
		if(projector==null) return;
		
		boolean x_left = false, y_left = false;

		int[] x=new int[5];
		int[] y=new int[5];

		factor_x = factor_y = 1;
		Point projection = projector.project(0, 0, -10);
		x[0] = projection.x;
		projection = projector.project(10.5f, 0, -10);
		y_left = projection.x > x[0];
		int i = projection.y;
		projection = projector.project(-10.5f, 0, -10);
		
		if (projection.y > i){
			factor_x = -1;
			y_left = projection.x > x[0];
		}
		
		projection = projector.project(0, 10.5f, -10);
		x_left = projection.x > x[0];
		i = projection.y;
		projection = projector.project(0, -10.5f, -10);
		
		if (projection.y > i){
			factor_y = -1;
			x_left = projection.x > x[0];
		}
		
		setAxesScale();
		drawBase(g, x, y);

		if (isBoxed){
			projection = projector.project(-factor_x*10,-factor_y*10,-10);
			x[0] = projection.x;
			y[0] = projection.y;
			projection = projector.project(-factor_x*10,-factor_y*10, 10);
			x[1] = projection.x;
			y[1] = projection.y;
			projection = projector.project( factor_x*10,-factor_y*10, 10);
			x[2] = projection.x;
			y[2] = projection.y;
			projection = projector.project( factor_x*10,-factor_y*10,-10);
			x[3] = projection.x;
			y[3] = projection.y;
			x[4] = x[0];
			y[4] = y[0];

			if (model.plotMode != PlotMode.WIREFRAME){
				if (model.plotMode == PlotMode.NORENDER&&g!=null) g.setColor(Color.lightGray);
				else if(g!=null) g.setColor(new Color(192,220,192));
				
				if(g!=null) g.fillPolygon(x, y, 4);
			}
			
			if(g!=null) g.setColor(Color.black);
			if(g!=null) g.drawPolygon(x, y, 5);

			projection = projector.project(-factor_x * 10, factor_y * 10, 10);
			x[2] = projection.x;
			y[2] = projection.y;
			projection = projector.project(-factor_x * 10, factor_y * 10, -10);
			x[3] = projection.x;
			y[3] = projection.y;
			x[4] = x[0];
			y[4] = y[0];

			if(model.plotMode!=PlotMode.WIREFRAME){
				if (model.plotMode==PlotMode.NORENDER&&g!=null) g.setColor(Color.lightGray);
				else if(g!=null) g.setColor(new Color(192, 220, 192));
				
				if(g!=null) g.fillPolygon(x, y, 4);
			}
			
			if(g!=null) g.setColor(Color.black);
			if(g!=null) g.drawPolygon(x,y,5);
			
		}else{
			if(isDisplayZ){
				projection = projector.project(factor_x * 10, -factor_y * 10, -10);
				x[0] = projection.x;
				y[0] = projection.y;
				projection = projector.project(factor_x * 10, -factor_y * 10, 10);
				g.drawLine(x[0], y[0], projection.x, projection.y);

				projection = projector.project(-factor_x * 10, factor_y * 10, -10);
				x[0] = projection.x;
				y[0] = projection.y;
				projection = projector.project(-factor_x * 10, factor_y * 10, 10);
				g.drawLine(x[0], y[0], projection.x, projection.y);
			}
		}

		for (i = -9; i <= 9; i++){
			if (isDisplayXY || isDisplayGrids){
				if (!isDisplayGrids || (i % (t_y / 2) == 0) || isDisplayXY){
					if (isDisplayGrids && (i % t_y == 0))
						projection = projector.project(-factor_x * 10, i, -10);
					else{
						if(i%t_y!= 0) projection = projector.project(factor_x*9.8f,i,-10);
						else projection = projector.project(factor_x * 9.5f, i, -10);
					}
					
					Point tickpos = projector.project(factor_x * 10, i, -10);
					if(g!=null) g.drawLine(projection.x, projection.y, tickpos.x, tickpos.y);
					
					if ((i % t_y == 0) && isDisplayXY){
						tickpos = projector.project(factor_x * 10.5f, i, -10);
						if (y_left)
							drawNumber(g, tickpos.x, tickpos.y,
									 (float) ((double) (i + 10) / 20 * (ymax - ymin) + ymin),
									 Label.LEFT, TOP);
						else
							drawNumber(g, tickpos.x, tickpos.y,
									 (float) ((double) (i + 10) / 20 * (ymax - ymin) + ymin),
									 Label.RIGHT, TOP);
					}
				}
				
				if(!isDisplayGrids || (i % (t_x / 2) == 0) || isDisplayXY){
					if(isDisplayGrids&&(i%t_x==0)) projection=projector.project(i,-factor_y*10,-10);
					else{
						if(i%t_x!= 0) projection = projector.project(i,factor_y*9.8f,-10);
						else projection = projector.project(i, factor_y * 9.5f, -10);
					}
					
					Point tickpos = projector.project(i, factor_y * 10, -10);
					if(g!=null) g.drawLine(projection.x, projection.y, tickpos.x, tickpos.y);
					
					if ((i % t_x == 0) && isDisplayXY){
						tickpos = projector.project(i, factor_y * 10.5f, -10);
						
						if (x_left)
							drawNumber(g, tickpos.x, tickpos.y,
									 (float) ((double) (i + 10) / 20 * (xmax - xmin) + xmin),
									 Label.LEFT, TOP);
						else
							drawNumber(g, tickpos.x, tickpos.y,
									 (float) ((double) (i + 10) / 20 * (xmax - xmin) + xmin),
									 Label.RIGHT, TOP);
					}
				}
			}
			
			// z grids and ticks
			if (isDisplayZ || (isDisplayGrids && isBoxed)){
				if (!isDisplayGrids || (i % (t_z / 2) == 0) || isDisplayZ){
					Point tickpos=null;
					if (isBoxed && isDisplayGrids && (i % t_z == 0)){
						projection = projector.project(-factor_x * 10, -factor_y * 10, i);
						tickpos = projector.project(-factor_x * 10, factor_y * 10, i);
						
					}else{
						if(i%t_z==0) projection = projector.project(-factor_x*10,factor_y*9.5f,i);
						else projection = projector.project(-factor_x * 10, factor_y * 9.8f, i);
						
						tickpos = projector.project(-factor_x * 10, factor_y * 10, i);
					}
					
					if(g!=null) g.drawLine(projection.x, projection.y, tickpos.x, tickpos.y);
					
					if(isDisplayZ){
						tickpos = projector.project(-factor_x * 10, factor_y * 10.5f, i);
						
						if(i%t_z==0){
							if(x_left)
								drawNumber(g, tickpos.x, tickpos.y,
										 (float) ((double) (i + 10) / 20 * (zmax - zmin) + zmin),
										 Label.LEFT, CENTER);
							else
								drawNumber(g, tickpos.x, tickpos.y,
										 (float) ((double) (i + 10) / 20 * (zmax - zmin) + zmin),
										 Label.RIGHT, CENTER);
						}
					}
					
					if (isDisplayGrids && isBoxed && (i % t_z == 0)){
						projection = projector.project(-factor_x * 10, -factor_y * 10, i);
						tickpos = projector.project(factor_x * 10, -factor_y * 10, i);
						
					}else{
						if (i % t_z == 0) projection = projector.project(factor_x * 9.5f, -factor_y * 10, i);
						else projection = projector.project(factor_x * 9.8f, -factor_y * 10, i);
						
						tickpos = projector.project(factor_x * 10, -factor_y * 10, i);
					}
					
					if(g!=null) g.drawLine(projection.x, projection.y, tickpos.x, tickpos.y);
					
					if(isDisplayZ){
						tickpos = projector.project(factor_x * 10.5f, -factor_y * 10, i);
						
						if(i%t_z==0){
							if (y_left)
								drawNumber(g, tickpos.x, tickpos.y,
										 (float) ((double) (i + 10) / 20 * (zmax - zmin) + zmin),
										 Label.LEFT, CENTER);
							else
								drawNumber(g, tickpos.x, tickpos.y,
										 (float) ((double) (i + 10) / 20 * (zmax - zmin) + zmin),
										 Label.RIGHT, CENTER);
						}
					}
					
					if(isDisplayGrids && isBoxed){
						if(i%t_y==0){
							projection = projector.project(-factor_x * 10, i, -10);
							tickpos = projector.project(-factor_x * 10, i, 10);
							if(g!=null) g.drawLine(projection.x, projection.y, tickpos.x, tickpos.y);
						}
						
						if(i%t_x==0){
							projection = projector.project(i, -factor_y * 10, -10);
							tickpos = projector.project(i, -factor_y * 10, 10);
							if(g!=null) g.drawLine(projection.x, projection.y, tickpos.x, tickpos.y);
						}
					}
				}
			}
		}
		
		if (isDisplayXY){
			Point tickpos=null;
			tickpos = projector.project(0, factor_y * 14, -10);
			drawString(g, tickpos.x, tickpos.y, model.xlabel, Label.CENTER, TOP);
			tickpos = projector.project(factor_x * 14, 0, -10);
			drawString(g, tickpos.x, tickpos.y, model.ylabel, Label.CENTER, TOP);
		}
		
		if(isDisplayZ){
			Point tickpos=projector.project(-factor_x * 10, factor_y * 14, 0);
			drawString(g, tickpos.x, tickpos.y, model.zlabel, Label.CENTER, TOP);
		}
	}
	
	/**
	 * Draws string at the specified coordinates with the specified alignment.
	 *
	 * @param g       graphics context to draw
	 * @param x       the x coordinate
	 * @param y       the y coordinate
	 * @param s       the string to draw
	 * @param x_align the alignment in x direction
	 * @param y_align the alignment in y direction
	 */
	private void drawString(Graphics g, int x, int y,
								 String s, int x_align, int y_align)
	{
		switch (y_align)
		{
			case TOP:
				if(g!=null) y += g.getFontMetrics(g.getFont()).getAscent();
				break;
			case CENTER:
				if(g!=null) y += g.getFontMetrics(g.getFont()).getAscent() / 2;
				break;
		}
		switch (x_align)
		{
			case Label.LEFT:
				if(g!=null) g.drawString(s, x, y);
				break;
			case Label.RIGHT:
				if(g!=null) g.drawString(s, x - g.getFontMetrics(
					g.getFont()).stringWidth(s), y);
				break;
			case Label.CENTER:
				if(g!=null) g.drawString(s, x - g.getFontMetrics(
					g.getFont()).stringWidth(s) / 2, y);
				break;
		}
	}
	
	/**Plotting routines and methods begin here ***/
	

	/**
	 * Draws float at the specified coordinates with the specified alignment.
	 *
	 * @param g       graphics context to draw
	 * @param x       the x coordinate
	 * @param y       the y coordinate
	 * @param f       the float to draw
	 * @param x_align the alignment in x direction
	 * @param y_align the alignment in y direction
	 */
	private void drawNumber(Graphics g,int x,int y,float f,int x_align,int y_align){
		String s = Float.toString(f);
		drawString(g, x, y, s, x_align, y_align);
	}
	
	/**
	 * Sets the axes scaling factor. Computes the proper axis lengths
	 * based on the ratio of variable ranges. The axis lengths will
	 * also affect the size of bounding box.
	 */
	private void setAxesScale(){
		float divisor;
		int longest;

		if(!isScaleBox){
			projector.setScaling(1f);
			t_x = t_y = t_z = 4;
			return;
		}

		float scale_x = xmax - xmin;
		float scale_y = ymax - ymin;
		float scale_z = zmax - zmin;

		if (scale_x < scale_y){
			if (scale_y < scale_z){
				longest = 3;
				divisor = scale_z;
				
			}else{
				longest = 2;
				divisor = scale_y;
			}
			
		}else{
			if (scale_x < scale_z){
				longest = 3;
				divisor = scale_z;
				
			}else{
				longest = 1;
				divisor = scale_x;
			}
		}
		scale_x /= divisor;
		scale_y /= divisor;
		scale_z /= divisor;

		if ((scale_x < 0.2f) || (scale_y < 0.2f) && (scale_z < 0.2f)){
			switch (longest){
				case 1:
					if (scale_y < scale_z){
						scale_y /= scale_z;
						scale_z = 1.0f;
						
					}else{
						scale_z /= scale_y;
						scale_y = 1.0f;
					}
					break;
				case 2:
					if (scale_x < scale_z){
						scale_x /= scale_z;
						scale_z = 1.0f;
						
					}else{
						scale_z /= scale_x;
						scale_x = 1.0f;
					}
					break;
				case 3:
					if (scale_y < scale_x){
						scale_y /= scale_x;
						scale_x = 1.0f;
						
					}else{
						scale_x /= scale_y;
						scale_y = 1.0f;
					}
					break;
			}
		}
		
		if(scale_x < 0.2f) scale_x = 1.0f;
		projector.setXScaling(scale_x);
		if(scale_y < 0.2f) scale_y = 1.0f;
		projector.setYScaling(scale_y);
		if(scale_z < 0.2f) scale_z = 1.0f;
		projector.setZScaling(scale_z);

		if(scale_x < 0.5f) t_x = 8;
		else t_x = 4;
		if(scale_y < 0.5f) t_y = 8;
		else t_y = 4;
		if(scale_z < 0.5f) t_z = 8;
		else t_z = 4;
	}
	
	/**
	 * Determines whether a plane is plottable, i.e: does not have
	 * invalid vertex.
	 *
	 * @return <code>true</code> if the plane is plottable,
	 *         <code>false</code> otherwise
	 * @param values vertices array of the plane
	 */
	private static boolean isPointsValid(Point3D[] values){
		return (!values[0].isInvalid()
			&&  !values[1].isInvalid()
			&&  !values[2].isInvalid()
			&&  !values[3].isInvalid());
	}
	
	private void plotPlane(Point3D[] vertex,int verticescount){
		Point projection;
		
		int count = 0;
		float z = 0.0f;
		
		boolean low1 = (vertex[0].z < zmin);
		boolean valid1 = !low1 && (vertex[0].z <= zmax);
		for(int loop=0,index=1;loop<verticescount;loop++){
			boolean low2 = (vertex[index].z < zmin);
			boolean valid2 = !low2 && (vertex[index].z <= zmax);
			
			if((valid1||valid2)||(low1^low2)){
				float result;
				
				if (!valid1){
					if(low1) result = zmin;
					else result = zmax;
					
					float ratio = (result - vertex[index].z) / (vertex[loop].z - vertex[index].z);
					float new_x = ratio * (vertex[loop].x - vertex[index].x) + vertex[index].x;
					float new_y = ratio * (vertex[loop].y - vertex[index].y) + vertex[index].y;
					
					if (low1) projection = projector.project(new_x, new_y, -10);
					else projection = projector.project(new_x, new_y, 10);
					
					poly_x[count] = projection.x;
					poly_y[count] = projection.y;
					count++;
					z += result;
				}
				
				if (valid2){
					projection = vertex[index].projection();
					if(projection!=null){
						poly_x[count] = projection.x;
						poly_y[count] = projection.y;
						count++;
						z += vertex[index].z;
					}
				}else{
					if(low2) result = zmin;
					else result = zmax;
					
					float ratio = (result - vertex[loop].z) / (vertex[index].z - vertex[loop].z);
					float new_x = ratio * (vertex[index].x - vertex[loop].x) + vertex[loop].x;
					float new_y = ratio * (vertex[index].y - vertex[loop].y) + vertex[loop].y;
					
					if (low2)projection = projector.project(new_x, new_y, -10);
					else projection = projector.project(new_x, new_y, 10);
					
					poly_x[count] = projection.x;
					poly_y[count] = projection.y;
					count++;
					z += result;
				}
			}
			
			if(++index==verticescount) index=0;
			
			valid1 = valid2;
			low1 = low2;
		}
		
		line_color = Color.darkGray;
		if(count>0){
			switch (model.plotMode){
				case NORENDER:
					if(bufferGraphics!=null) bufferGraphics.setColor(Color.lightGray);
					break;
				case SPECTRUM:
					z = 0.8f - (z / count - zmin) * color_factor;
					if(bufferGraphics!=null)
					bufferGraphics.setColor(Color.getHSBColor(z, 1.0f, 1.0f));
					break;
				case GRAYSCALE:
					z = (z / count - zmin) * color_factor;
					if(bufferGraphics!=null) bufferGraphics.setColor(Color.getHSBColor(0, 0, z));
					if (z < 0.3f) line_color = new Color(0.6f, 0.6f, 0.6f);
					break;
				case DUALSHADE:
					z = (z / count - zmin) * color_factor + 0.4f;
					if(bufferGraphics!=null) bufferGraphics.setColor(Color.getHSBColor(color, 0.7f, z));
					break;
				default:
					break;
			}
			
			if(bufferGraphics!=null)
			bufferGraphics.fillPolygon(poly_x, poly_y, count);
			
			if(bufferGraphics!=null)
			bufferGraphics.setColor(line_color);
			
			if(isMesh||(model.plotMode==SurfacePlotModel.PlotMode.NORENDER)){
				poly_x[count] = poly_x[0];
				poly_y[count] = poly_y[0];
				count++;
				if(bufferGraphics!=null)
				bufferGraphics.drawPolygon(poly_x, poly_y, count);
			}
		}
	}
	
	private void plotArea(int start_lx,int start_ly,int end_lx,int end_ly,int sx,int sy){
		Point3D values1[] = new Point3D[4];
		
		start_lx *= calc_divisionsY + 1;
		sx *= calc_divisionsY + 1;
		end_lx *= calc_divisionsY + 1;
		
		int lx = start_lx;
		int ly = start_ly;
		
		while(ly != end_ly){
			values1[1] = vertex[lx + ly];
			values1[2] = vertex[lx + ly + sy];
			while(lx != end_lx){
				values1[0] = values1[1];
				values1[1] = vertex[lx + sx + ly];
				values1[3] = values1[2];
				values1[2] = vertex[lx + sx + ly + sy];
				
				if(model.plotMode==SurfacePlotModel.PlotMode.DUALSHADE) color=0.2f;
				if(isPointsValid(values1)) plotPlane(values1, 4);
				
				lx += sx;
			}
			
			ly += sy;
			lx = start_lx;
		}
	}
	
	private void plotSurface(){
		image_drawn = false;
		
		float zi=model.zmin;	zmin = zi;
		float zx=model.zmax;	zmax = zx;
		
		color_factor = 0.8f / (zmax - zmin);
		if(model.plotMode==SurfacePlotModel.PlotMode.DUALSHADE) color_factor*=0.6f/0.8f;

		if (!printing){
			if(bufferGraphics!=null) bufferGraphics.setColor(Color.lightGray);
			if(bufferGraphics!=null) bufferGraphics.fillRect(0, 0, getBounds().width, getBounds().height);
		}
		
		drawBoxGridsTicksLabels(bufferGraphics);
		
		Point3D.setZRange(zmin, zmax);
		
		// direction test
		float distance=projector.getDistance()*projector.getCosElevationAngle();

		// cop : center of projection
		cop=new Point3D(
			distance * projector.getSinRotationAngle(),
			distance * projector.getCosRotationAngle(),
			projector.getDistance()*projector.getSinElevationAngle()
		);
		
		cop.transform();
		
		int start_lx,sx,end_lx;
		int start_ly,sy,end_ly;

		int plot_density = model.dispDivision;
		int multiple_factor = 1;//calc_divisionsX / plot_density;
		
		if(cop.x>0){
			start_lx = 0;
			end_lx = calc_divisionsX;
			sx = multiple_factor;
		}else{
			start_lx = calc_divisionsX;
			end_lx = 0;
			sx = -multiple_factor;
		}
		
		if(cop.y>0){
			start_ly = 0;
			end_ly = calc_divisionsY;
			sy = multiple_factor;
		}else{
			start_ly = calc_divisionsY;
			end_ly = 0;
			sy = -multiple_factor;
		}

		if ((cop.x > 10) || (cop.x < -10)){
			if ((cop.y > 10) || (cop.y < -10)) plotArea(start_lx, start_ly, end_lx, end_ly, sx, sy);
			else{	// split in y direction
				int split_y = (int) ((cop.y + 10) * plot_density / 20) * multiple_factor;
				plotArea(start_lx, 0, end_lx, split_y, sx, multiple_factor);
				plotArea(start_lx, calc_divisionsY, end_lx, split_y, sx, -multiple_factor);
			}
		}else{
			if ((cop.y > 10) || (cop.y < -10)){   // split in x direction
				int split_x = (int) ((cop.x + 10) * plot_density / 20) * multiple_factor;
				plotArea(0, start_ly, split_x, end_ly, multiple_factor, sy);
				plotArea(calc_divisionsX, start_ly, split_x, end_ly, -multiple_factor, sy);
			}else{    // split in both x and y directions
				int split_x = (int) ((cop.x + 10) * plot_density / 20) * multiple_factor;
				int split_y = (int) ((cop.y + 10) * plot_density / 20) * multiple_factor;
				plotArea(0, 0, split_x, split_y, multiple_factor, multiple_factor);
				plotArea(0, calc_divisionsY, split_x, split_y, multiple_factor, -multiple_factor);
				plotArea(calc_divisionsX, 0, split_x, split_y, -multiple_factor, multiple_factor);
				plotArea(calc_divisionsX, calc_divisionsY, split_x, split_y,
						 -multiple_factor, -multiple_factor);
			}
		}
	}
	
	private void plotWireframe(){
		float lx = 0, ly = 0, lastz = 0;
		boolean error, lasterror, invalid;
		
		image_drawn = false;
		
		Point lastproj=new Point(0,0);
		Point projection = new Point(0, 0);
		
		int plot_density = model.dispDivision;
		int multiple_factor = calc_divisionsX / plot_density;
		
		zmin = model.zmin;
		zmax = model.zmax;
		
		if(!printing){
			if(bufferGraphics!=null) bufferGraphics.setColor(Color.lightGray);
			if(bufferGraphics!=null) bufferGraphics.fillRect(0, 0, getBounds().width, getBounds().height);
		}

		drawBoxGridsTicksLabels(bufferGraphics);
		if(bufferGraphics!=null) bufferGraphics.setColor(Color.black);

		Point3D.setZRange(zmin, zmax);
		
		int counter = 0;

		// plot - x direction
		for (int i=0,k=0;i <= calc_divisionsX;i++){
			lasterror = true;
			
			if (counter == 0){
				for (int j=0;j <= calc_divisionsY;j++){
					float z = vertex[k].z;
					invalid = Float.isNaN(z);
					if (!invalid){
						if (z < zmin){
							error = true;
							float ratio = (zmin - lastz) / (z - lastz);
							projection = projector.project(ratio * (vertex[k].x - lx) + lx,
														   ratio * (vertex[k].y - ly) + ly, -10);
						}else{
							if (z > zmax){
								error = true;
								float ratio = (zmax - lastz) / (z - lastz);
								projection = projector.project(ratio * (vertex[k].x - lx) + lx,
															   ratio * (vertex[k].y - ly) + ly, 10);
							}else{
								error = false;
								projection = vertex[k].projection();
							}
						}
						
						if(lasterror && (!error) && (j != 0)){
							if (lastz > zmax){
								float ratio = (zmax - z) / (lastz - z);
								lastproj = projector.project(
									ratio * (lx - vertex[k].x) + vertex[k].x,
									ratio * (ly - vertex[k].y) + vertex[k].y, 10);
							}else{
								if (lastz < zmin){
									float ratio = (zmin - z) / (lastz - z);
									lastproj = projector.project(
										ratio * (lx - vertex[k].x) + vertex[k].x,
										ratio * (ly - vertex[k].y) + vertex[k].y, -10);
								}
							}
							
						}else invalid = error && lasterror;
						
					}else error = true;
					
					if(!invalid&&(j!=0)&&bufferGraphics!=null)
					bufferGraphics.drawLine(lastproj.x,lastproj.y,projection.x,projection.y);
					
					lastproj = projection;
					lasterror = error;
					lx = vertex[k].x;
					ly = vertex[k].y;
					lastz = z;
					k++;
				}
				
			}else k += calc_divisionsX + 1;
			
			counter = (counter + 1) % multiple_factor;
		}
		
		// plot - y direction
		counter = 0;
		
		for(int j=0,k=0;j <= calc_divisionsY;j++){
			lasterror = true;
			
			if(counter==0)
			for(int i=0;i <= calc_divisionsX;i++){
				float z = vertex[k].z;
				invalid = Float.isNaN(z);
				
				if (!invalid){
					if (z < zmin){
						error = true;
						float ratio = (zmin - lastz) / (z - lastz);
						projection = projector.project(ratio * (vertex[k].x - lx) + lx,
													   ratio * (vertex[k].y - ly) + ly, -10);
					}else{
						if (z > zmax){
							error = true;
							float ratio = (zmax - lastz) / (z - lastz);
							projection = projector.project(ratio * (vertex[k].x - lx) + lx,
														   ratio * (vertex[k].y - ly) + ly, 10);
						}else{
							error = false;
							projection = vertex[k].projection();
						}
					}
					
					if (lasterror && (!error) && (i != 0)){
						if (lastz > zmax){
							float ratio = (zmax - z) / (lastz - z);
							lastproj = projector.project(
								ratio * (lx - vertex[k].x) + vertex[k].x,
								ratio * (ly - vertex[k].y) + vertex[k].y, 10);
						}else{
							if (lastz < zmin){
								float ratio = (zmin - z) / (lastz - z);
								lastproj = projector.project(
									ratio * (lx - vertex[k].x) + vertex[k].x,
									ratio * (ly - vertex[k].y) + vertex[k].y, -10);
							}
						}
						
					}else invalid = error && lasterror;
					
				}else error = true;
				
				if (!invalid && (i != 0) && bufferGraphics!=null)
				bufferGraphics.drawLine(lastproj.x, lastproj.y, projection.x, projection.y);
				
				lastproj = projection;
				lasterror = error;
				lx = vertex[k].x;
				ly = vertex[k].y;
				lastz = z;
				k += calc_divisionsX + 1;
			}
			
			k = j;
			counter = (counter + 1) % multiple_factor;
		}
	}
}
