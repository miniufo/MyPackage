/**
 * @(#)VisualizeDataInJOGL.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import miniufo.diagnosis.SpatialModel;
import java.awt.BorderLayout;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.media.opengl.GL2;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLCapabilities;
import javax.media.opengl.GLProfile;
import javax.media.opengl.awt.GLCanvas;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;


/**
 * visualize an array of 2-D discrete data
 */
public final class VisualizeInJOGLTest{
	//
	private JFrame   frame =null;
	private GLCanvas canvas=null;// RGB range from [0, 1]
	private Handler  handle=null;
	
	
	/**
	 * constructor
	 */
	public VisualizeInJOGLTest(){
		frame =new JFrame("JOGL Visualization");
		canvas=new GLCanvas(new GLCapabilities(GLProfile.get(GLProfile.GL2)));
		
		handle=new Handler(canvas);
		
		canvas.addMouseMotionListener(handle);
		canvas.addMouseListener(handle);
		canvas.addMouseWheelListener(handle);
		canvas.addGLEventListener(handle);
		canvas.addKeyListener(handle);
		
		frame.addWindowListener(new WindowAdapter(){
            public void windowClosing(WindowEvent e){
        		SwingUtilities.invokeLater(()->{
        			frame.dispose();
        			System.exit(0);
        		});
            }
		});
		frame.add(canvas,BorderLayout.CENTER);
		frame.pack();
		frame.setSize(800,600);
		frame.setVisible(true);
		
		canvas.requestFocus();
	}
	
	
	/*** helper methods and classes ***/
	static final class Handler extends GLHandler{
		//
		private static final int resolution=60;
		
		private static final double r=2.0;
		
		private static final Point O =new Point( 0, 0, 0);
		private static final Point X =new Point(2.5f,0,0);
		private static final Point Y =new Point(0,2.5f,0);
		private static final Point Z =new Point(0,0,2.5f);
		
		private static final Point[] arcXY=getXYArc(resolution);
		private static final Point[] arcXZ=getXZArc(resolution);
		private static final Point[] arcYZ=getYZArc(resolution);
		
		//
		private double lat0=35;
		private double lon0=45;
		private double beta=20;
		private double alfa=350;
		private double lat1=0;
		private double lon1=0;
		
		private Point r0=new Point(
			r*Math.cos(Math.toRadians(lat0))*Math.sin(Math.toRadians(lon0)),
			r*Math.sin(Math.toRadians(lat0)),
			r*Math.cos(Math.toRadians(lat0))*Math.cos(Math.toRadians(lon0))
		);
		
		private Point r1=getPoint1();
		private Point rL=getPointL(0.4);
		private Point rK=getPointK(0.4);
		private Point rM=getPointM(0.4);
		
		private Point[] arcLon0=getLon0Arc(resolution);
		private Point[] arcLat0=getLat0Arc(resolution);
		private Point[] arcRad =getRadialArc(resolution);
		private Point[] cylindC=getCylindCircle(resolution);
		private Point[] cylind2=getCylindCircle2(resolution);
		
		
		//
		public Handler(GLCanvas canvas){ super(canvas);}
		
		void drawLine(GL2 gl,Point p1,Point p2){
			gl.glVertex3f((float)p1.x,(float)p1.y,(float)p1.z);
			gl.glVertex3f((float)p2.x,(float)p2.y,(float)p2.z);
		}
		
		void drawPoint(GL2 gl,Point p){ gl.glVertex3f((float)p.x,(float)p.y,(float)p.z);}
		
		
		static Point[] getXYArc(int segs){
			Point[] re=new Point[segs+1];
			for(int i=0;i<=segs;i++){
				double a=Math.PI/2.0*i/segs;
				re[i]=new Point(r*Math.cos(a),r*Math.sin(a),0.0);
			}
			return re;
		}
		
		static Point[] getXZArc(int segs){
			Point[] re=new Point[segs+1];
			for(int i=0;i<=segs;i++){
				double a=Math.PI/2.0*i/segs;
				re[i]=new Point(r*Math.cos(a),0.0,r*Math.sin(a));
			}
			return re;
		}
		
		static Point[] getYZArc(int segs){
			Point[] re=new Point[segs+1];
			for(int i=0;i<=segs;i++){
				double a=Math.PI/2.0*i/segs;
				re[i]=new Point(0.0,r*Math.sin(a),r*Math.cos(a));
			}
			return re;
		}
		
		Point[] getLon0Arc(int segs){
			Point[] re=new Point[segs+1];
			for(int i=0;i<=segs;i++){
				double a=Math.PI/2.0*i/segs;
				re[i]=new Point(
					r*Math.sin(a)*Math.cos(Math.toRadians(lon0)),
					r*Math.cos(a),
					r*Math.sin(a)*Math.sin(Math.toRadians(lon0))
				);
			}
			return re;
		}
		
		Point[] getLat0Arc(int segs){
			Point[] re=new Point[segs+1];
			for(int i=0;i<=segs;i++){
				double a=Math.PI/2.0*i/segs;
				re[i]=new Point(
					r*Math.cos(a)*Math.cos(Math.toRadians(lat0)),
					r*Math.sin(Math.toRadians(lat0)),
					r*Math.sin(a)*Math.cos(Math.toRadians(lat0))
				);
			}
			return re;
		}
		
		Point[] getRadialArc(int segs){
			Point[] re=new Point[segs+1];
			for(int i=0;i<segs;i++){
				float[] re2=SpatialModel.cLatLon(lon0,lat0,360-alfa,i*beta/segs);
				double lon=(float)Math.toRadians(re2[0]);
				double lat=(float)Math.toRadians(re2[1]);
				re[i]=new Point(
					r*Math.cos(lat)*Math.cos(lon),
					r*Math.sin(lat),
					r*Math.cos(lat)*Math.sin(lon)
				);
			}
			re[segs]=re[0];
			return re;
		}
		
		Point[] getCylindCircle(int segs){
			Point[] re=new Point[segs+1];
			float[][][] ll=SpatialModel.cLatLons(Math.toRadians(lon0),Math.toRadians(lat0),Math.toRadians(beta),2,segs);
			
			for(int i=0;i<segs;i++) re[i]=new Point(
				r*Math.cos(ll[1][1][i])*Math.cos(ll[0][1][i]),
				r*Math.sin(ll[1][1][i]),
				r*Math.cos(ll[1][1][i])*Math.sin(ll[0][1][i])
			);
			
			re[segs]=re[0];
			return re;
		}
		
		Point[] getCylindCircle2(int segs){
			Point[] re=new Point[segs+1];
			float[][][] ll=SpatialModel.cLatLons(Math.toRadians(lon0),Math.toRadians(lat0),Math.toRadians(beta*0.95),2,segs);
			
			for(int i=0;i<segs;i++) re[i]=new Point(
				r*1.06*Math.cos(ll[1][1][i])*Math.cos(ll[0][1][i]),
				r*1.06*Math.sin(ll[1][1][i]),
				r*1.06*Math.cos(ll[1][1][i])*Math.sin(ll[0][1][i])
			);
			
			re[segs]=re[0];
			return re;
		}
		
		Point getPoint1(){
			float[] ll=SpatialModel.cLatLon(lon0,lat0,360-alfa,beta);
			lon1=ll[0];
			lat1=ll[1];
			double lontmp=Math.toRadians(lon1);
			double lattmp=Math.toRadians(lat1);
			return new Point(
				r*Math.cos(lattmp)*Math.cos(lontmp),
				r*Math.sin(lattmp),
				r*Math.cos(lattmp)*Math.sin(lontmp)
			);
		}
		
		Point getPointL(double mod){
			double l=r0.y*r1.z-r0.z*r1.y;
			double m=r0.z*r1.x-r0.x*r1.z;
			double n=r0.x*r1.y-r0.y*r1.x;
			
			double ratio=mod/Math.sqrt(l*l+m*m+n*n);
			
			return new Point(l*ratio+r1.x,m*ratio+r1.y,n*ratio+r1.z);
		}
		
		Point getPointK(double mod){
			double ratio=(mod+r)/r;
			return new Point(r1.x*ratio,r1.y*ratio,r1.z*ratio);
		}
		
		Point getPointM(double mod){
			double Kx=rK.x, Lx=rL.x-r1.x;
			double Ky=rK.y, Ly=rL.y-r1.y;
			double Kz=rK.z, Lz=rL.z-r1.z;
			
			double x=Ly*Kz-Lz*Ky;
			double y=Lz*Kx-Lx*Kz;
			double z=Lx*Ky-Ly*Kx;
			
			double ratio=mod/Math.sqrt(x*x+y*y+z*z);
			
			return new Point(x*ratio+r1.x,y*ratio+r1.y,z*ratio+r1.z);
		}
		
		double getDistance(Point p1,Point p2){
			double dx=p1.x-p2.x; dx*=dx;
			double dy=p1.y-p2.y; dy*=dy;
			double dz=p1.z-p2.z; dz*=dz;
			return Math.sqrt(dx+dy+dz);
		}
		
		void drawAll(GL2 gl,int segs){
			// coordinate points
			gl.glBegin(GL2.GL_POINTS);
			drawPoint(gl,O);
			drawPoint(gl,X);
			drawPoint(gl,Y);
			drawPoint(gl,Z);
			drawPoint(gl,r0);
			drawPoint(gl,r1);
			drawPoint(gl,rL);
			drawPoint(gl,rK);
			drawPoint(gl,rM);
			gl.glEnd();
			
			gl.glBegin(GL2.GL_LINES);
			// coordinate lines
			drawLine(gl,O,X);
			drawLine(gl,O,Y);
			drawLine(gl,O,Z);
			drawLine(gl,O,r0);
			drawLine(gl,O,r1);
			drawLine(gl,r1,rL);
			drawLine(gl,r1,rK);
			drawLine(gl,r1,rM);
			
			for(int i=0;i<segs;i++) drawLine(gl,arcXY[i],arcXY[i+1]);	// for X-Y arc
			for(int i=0;i<segs;i++) drawLine(gl,arcXZ[i],arcXZ[i+1]);	// for X-Z arc
			for(int i=0;i<segs;i++) drawLine(gl,arcYZ[i],arcYZ[i+1]);	// for Y-Z arc
			
			for(int i=0;i<segs;i++) drawLine(gl,arcLon0[i],arcLon0[i+1]);	// for longitude lambda0 arc
			for(int i=0;i<segs;i++) drawLine(gl,arcLat0[i],arcLat0[i+1]);	// for latitude fai0 arc
			for(int i=0;i<segs;i++) drawLine(gl,arcRad [i],arcRad [i+1]);	// for radial arc
			for(int i=0;i<segs;i++) drawLine(gl,cylindC[i],cylindC[i+1]);	// for cylindrical circle
			for(int i=0;i<segs;i++) drawLine(gl,cylind2[i],cylind2[i+1]);	// for cylindrical circle
			gl.glEnd();
		}
		
		
		/*** overwrite methods ***/
		public void display(GLAutoDrawable gld){
			super.display(gld);
			
			GL2 gl=gld.getGL().getGL2();
			
			gl.glLineWidth(1.6f);
			gl.glPointSize(5f);
			
			drawAll(gl,60);
		}
		
		public void keyPressed(KeyEvent e){
			int code=e.getKeyCode();
			
			if(     code==KeyEvent.VK_LEFT ) alfa++;
			else if(code==KeyEvent.VK_RIGHT) alfa--;
			
			r1=getPoint1();
			rL=getPointL(0.4);
			rK=getPointK(0.4);
			rM=getPointM(0.4);
			
			arcRad =getRadialArc(resolution);
			cylindC=getCylindCircle(resolution);
			cylind2=getCylindCircle2(resolution);
			
			canvas.repaint();
		}
	}
	
	
	private static final class Point{
		double x=0;
		double y=0;
		double z=0;
		
		public Point(double x,double y,double z){
			this.x=x;
			this.y=y;
			this.z=z;
		}
	}
	
	
	/** test*/
	public static void main(String args[]){
		new VisualizeInJOGLTest();
    }
}
