/**
 * @(#)GLHandler	1.0 2017.09.19
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import java.awt.Component;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import javax.media.opengl.GL2;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.awt.GLCanvas;
import javax.media.opengl.glu.GLU;
import com.jogamp.opengl.util.gl2.GLUT;


/**
 * A handler template for JOGL display
 *
 * @version 1.0, 2017.09.19
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class GLHandler implements GLEventListener,KeyListener,MouseMotionListener,MouseListener,MouseWheelListener{
	//
	protected int   pos_x =0;
	protected int   pos_y =0;
	
	protected float angleZ=0;
	protected float angleY=0;
	protected float angleX=0;
	
	protected float transX=0;
	protected float transY=0;
	protected float transZ=0;
	
	protected GLU      glu   =null;
	protected GLUT     glut  =null;
	protected GLCanvas canvas=null;
	
	
	/**
	 * Constructor
	 */
	public GLHandler(GLCanvas canvas){
		this.canvas=canvas;
		this.glu =new GLU();
		this.glut=new GLUT();
	}
	
	
	/*** overwrite methods ***/
	public void init(GLAutoDrawable gld){
		GL2 gl=gld.getGL().getGL2();
		
		gl.glEnable(GL2.GL_DEPTH_TEST);
		
		gl.glEnableClientState(GL2.GL_VERTEX_ARRAY);
		gl.glEnableClientState(GL2.GL_COLOR_ARRAY);
		
		gl.glMatrixMode(GL2.GL_PROJECTION);
		gl.glLoadIdentity();
		
		double w=((Component)gld).getWidth();
		double h=((Component)gld).getHeight();
		double aspect=w/h;
		
		glu.gluPerspective(20,aspect,0.1,20);
	}
	
	public void display(GLAutoDrawable gld){
		GL2 gl=gld.getGL().getGL2();
		
		gl.glClearColor(1f,1f,1f,1);
        gl.glClear(GL2.GL_COLOR_BUFFER_BIT | GL2.GL_DEPTH_BUFFER_BIT);
		gl.glColor3f(0,0,0);
		
		gl.glMatrixMode(GL2.GL_MODELVIEW);
		gl.glLoadIdentity();
		
		glu.gluLookAt(0,0,15, 0,0,0, 0,1,0);
		
		gl.glTranslatef(transX,transY,transZ);
		
		gl.glRotatef(angleX,1,0,0);
		gl.glRotatef(angleY,0,1,0);
		gl.glRotatef(angleZ,0,0,1);
	}
	
	public void reshape(GLAutoDrawable gld,int x,int y,int width,int height){ canvas.repaint();}
	
	public void dispose(GLAutoDrawable gld){}
	
	
	/**
	 * Keyboard events.
	 */
	public void keyPressed(KeyEvent e){}
	
	public void keyReleased(KeyEvent e){}
	
	public void keyTyped(KeyEvent e){
		float incre=e.getModifiers()==1?10f:2f;
		
		switch(e.getKeyChar()){
		case KeyEvent.VK_D+32:
		case KeyEvent.VK_D:
			angleY+=incre; if(angleY>360) angleY-=360;
			break;
		case KeyEvent.VK_A+32:
		case KeyEvent.VK_A:
			angleY-=incre; if(angleY<0) angleY+=360;
			break;
		case KeyEvent.VK_W+32:
		case KeyEvent.VK_W:
			angleX+=incre; if(angleX>90) angleX=90;
			break;
		case KeyEvent.VK_S+32:
		case KeyEvent.VK_S:
			angleX-=incre; if(angleX<0) angleX=0;
			break;
		default: break;
		}
		
		canvas.repaint();
	}
	
	
	/**
	 * Mouse events.
	 */
	public void mouseDragged(MouseEvent e){
		int x=e.getX(),dx=x-pos_x;
		int y=e.getY(),dy=y-pos_y;
		
		if(e.isControlDown()){
			transX+=dx/100f;
			transY-=dy/100f;
			
		}else{
			angleY+=dx;
			angleX+=dy;
			
			if(angleY>360) angleY-=360;
			else if(angleY<0) angleY+=360;
			
			if(angleX>90) angleX=90;
			if(angleX<0) angleX=0;
		}
		
		canvas.repaint();
		pos_x=x;
		pos_y=y;
	}
	
	public void mouseMoved(MouseEvent e){}
	
	public void mouseClicked(MouseEvent e){}
	
	public void mouseEntered(MouseEvent e){}
	
	public void mouseExited(MouseEvent e){}
	
	public void mousePressed(MouseEvent e){ pos_x=e.getX(); pos_y=e.getY();}
	
	public void mouseReleased(MouseEvent e){}
	
	public void mouseWheelMoved(MouseWheelEvent e){
		transZ-=e.getWheelRotation()*0.1f;
		
		if(transZ<-5) transZ=-5;
		if(transZ> 3) transZ= 3;
		
		canvas.repaint();
	}
}
