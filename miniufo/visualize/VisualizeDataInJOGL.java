/**
 * @(#)VisualizeDataInJOGL.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import com.jogamp.opengl.util.GLBuffers;
import com.jogamp.opengl.util.gl2.GLUT;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.concurrent.TimeUnit;
import javax.media.opengl.GL2;
import javax.media.opengl.GL3;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLCapabilities;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.GLProfile;
import javax.media.opengl.awt.GLCanvas;
import javax.media.opengl.glu.GLU;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;


/**
 * visualize an array of 2-D discrete data
 */
public final class VisualizeDataInJOGL{
	//
	private JFrame   frame =null;
	private GLCanvas canvas=null;// RGB range from [0, 1]
	private Primitives prmt=null;
	private Handler  handle=null;
	
	private static final float[][] GRADIENT_COLORS=new float[239][3];
	
	static{	// initialize the gradient colors
		int r=255,g=0,b=140;
		
		for(int i=  0;i< 10;i++){ b-=10; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
		for(int i= 10;i< 18;i++){ b-=5 ; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
		for(int i= 18;i< 69;i++){ g+=5 ; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
		for(int i= 69;i<120;i++){ r-=5 ; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
		for(int i=120;i<171;i++){ b+=5 ; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
		for(int i=171;i<222;i++){ g-=5 ; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
		for(int i=222;i<229;i++){ r+=5 ; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
		for(int i=229;i<239;i++){ r+=10; GRADIENT_COLORS[i][0]=r/255f; GRADIENT_COLORS[i][1]=g/255f; GRADIENT_COLORS[i][2]=b/255f;}
	}
	
	
	/**
	 * constructor
	 */
	public VisualizeDataInJOGL(float[][] data){
		frame =new JFrame("JOGL Visualization");
		canvas=new GLCanvas(new GLCapabilities(GLProfile.get(GLProfile.GL2)));
		prmt=new Primitives(data,2);
		
		handle=new Handler();
		
		handle.setBuffer(prmt);
		
		canvas.addMouseMotionListener(handle);
		canvas.addMouseListener(handle);
		canvas.addMouseWheelListener(handle);
		canvas.addGLEventListener(handle);
		canvas.addKeyListener(handle);
		
		frame.addWindowListener(new WindowAdapter(){
            public void windowClosing(WindowEvent e){
                runExit();
            }
		});
		frame.add(canvas,BorderLayout.CENTER);
		frame.pack();
		frame.setSize(800,600);
		frame.setVisible(true);
		
		canvas.requestFocus();
	}
	
	
	public void setData(float[][] data){
		prmt.setData(data);
		handle.setBuffer(prmt);
		canvas.repaint();
	}
	
	
	/*** helper methods and classes ***/
	private void runExit(){
		SwingUtilities.invokeLater(new Runnable(){
			public void run(){
				frame.dispose();
				System.exit(0);
			}
		});
	}
	
	
	private final class Handler
	implements GLEventListener,KeyListener,MouseMotionListener,MouseListener,MouseWheelListener{
		//
		private int pos_x,pos_y;
		
		private float angleZ=0;
		private float angleY=0;
		private float angleX=0;
		
		private float translateX=0;
		private float translateY=0;
		private float translateZ=0;
		
		private IntBuffer    indices=null;
		private FloatBuffer vertices=null;
		private FloatBuffer colors  =null;
		
		
		//
		public Handler(){}
		
		
		/*** getor and setor ***/
		public void setBuffer(Primitives prmt){
			this.vertices=prmt.getVertices();
			this.colors  =prmt.getColors();
			this.indices =prmt.getIndices();
		}
		
		//
		public void init(GLAutoDrawable gld){
			GL2 gl=gld.getGL().getGL2();
			GLU glu=new GLU();
			
			gl.glEnable(GL2.GL_DEPTH_TEST);
			
			gl.glEnableClientState(GL2.GL_VERTEX_ARRAY);
			gl.glEnableClientState(GL2.GL_COLOR_ARRAY);
			gl.glVertexPointer(3, GL2.GL_FLOAT, 0,vertices);
			gl.glColorPointer( 3, GL2.GL_FLOAT, 0, colors );
			
			gl.glMatrixMode(GL2.GL_PROJECTION);
			gl.glLoadIdentity();
			
			double w=((Component)gld).getWidth();
			double h=((Component)gld).getHeight();
			double aspect=w/h;
			
			glu.gluPerspective(60,aspect,0.1,20);
		}
		
		public void display(GLAutoDrawable gld){
			GL2 gl=gld.getGL().getGL2();
			GLU glu=new GLU();
			GLUT glut=new GLUT();
			
			gl.glClearColor(0.8f,0.8f,0.8f,1);
	        gl.glClear(GL2.GL_COLOR_BUFFER_BIT | GL2.GL_DEPTH_BUFFER_BIT);
			gl.glColor3f(0,0,0);
			
			gl.glMatrixMode(GL2.GL_MODELVIEW);
			gl.glLoadIdentity();
			
			glu.gluLookAt(0,0,4,  0,0,0,  0,1,0);
			
			gl.glTranslatef(translateX,translateY,translateZ);
			
			gl.glRotatef(angleX,1,0,0);
			gl.glRotatef(angleY,0,1,0);
			gl.glRotatef(angleZ,0,0,1);
			
	        gl.glDrawElements(GL2.GL_QUADS,indices.capacity(),GL3.GL_UNSIGNED_INT,indices);
			glut.glutWireCube(2);
		}
		
		public void reshape(GLAutoDrawable gld,int x,int y,int width,int height){canvas.repaint();}
		
		public void dispose(GLAutoDrawable gld){}
		
		//
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
		
		//
		public void mouseDragged(MouseEvent e){
			int x=e.getX(),dx=x-pos_x;
			int y=e.getY(),dy=y-pos_y;
			
			if(e.isControlDown()){
				translateX+=dx/100f;
				translateY-=dy/100f;
				
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
			translateZ-=e.getWheelRotation()*0.1f;
			
			if(translateZ<-5) translateZ=-5;
			if(translateZ> 3) translateZ= 3;
			
			canvas.repaint();
		}
	}
	
	private final class Primitives{
		//
		private int y;
		private int x;
		
		private float min=0;
		private float max=0;
		private float scale=0;
		
		private int[] indices=null;
		
		private float[][] data=null;
		
		private IntBuffer   idx=null;	// indices for primitives, generate once
		private FloatBuffer clr=null;	// define colors, change if data change
		private FloatBuffer vtc=null;	// define vertices, change if data change
		
		
		/**
		 * constructor
		 */
		public Primitives(float[][] data,float scale){
			y=data.length;
			x=data[0].length;
			
			if(x<2) throw new IllegalArgumentException("x-length should be larger than 1");
			if(y<2) throw new IllegalArgumentException("y-length should be larger than 1");
			
			this.scale=scale;
			this.data=new float[y][];
			
			for(int j=0;j<y;j++) this.data[j]=data[j].clone();
			
			scale(this.data,scale);
			
			indices=new int[(y-1)*(x-1)*4];
			
			clr=GLBuffers.newDirectFloatBuffer(x*y*3);
			vtc=GLBuffers.newDirectFloatBuffer(x*y*3);
			
			idx=genIndices();
			genColors();
			genVertices(-1,1,-1,1);
		}
		
		
		/*** getor and setor ***/
		public IntBuffer getIndices(){ return idx;}
		
		public FloatBuffer getColors(){ return clr;}
		
		public FloatBuffer getVertices(){ return vtc;}
		
		
		public void setData(float[][] data){
			for(int j=0;j<y;j++) System.arraycopy(data[j],0,this.data[j],0,x);
			
			scale(this.data,scale);
			
			genColors();
			updateVertices();
		}
		
		
		/*** helper methods ***/
		private IntBuffer genIndices(){
			for(int j=0,J=y-1,ptr=0;j<J;j++)
			for(int i=0,I=x-1;i<I;i++){
				indices[ptr++]= j   *x+i  ;
				indices[ptr++]= j   *x+i+1;
				indices[ptr++]=(j+1)*x+i+1;
				indices[ptr++]=(j+1)*x+i  ;
			}
			
			return IntBuffer.wrap(indices);
		}
		
		private void genVertices(float xstr,float xend,float ystr,float yend){
			float xinc=(xend-xstr)/(x-1);
			float yinc=(yend-ystr)/(y-1);
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				vtc.put(xstr+xinc*i);
				vtc.put(data[j][i]);
				vtc.put(ystr+yinc*j);
			}
			
			vtc.rewind();
		}
		
		private void updateVertices(){
			int pos=1;
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				vtc.put(pos,data[j][i]);
				pos+=3;
			}
			
			vtc.rewind();
		}
		
		private void genColors(){
			float step=(max-min)/(GRADIENT_COLORS.length-1);
			
			if(step>1e-11f){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					int tag=(int)((max-data[j][i])/step);
					
					float[] color=GRADIENT_COLORS[tag];
					
					clr.put(color[0]);
					clr.put(color[1]);
					clr.put(color[2]);
				}
				
			}else{
				float[] color=GRADIENT_COLORS[GRADIENT_COLORS.length/2];
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					clr.put(color[0]);
					clr.put(color[1]);
					clr.put(color[2]);
				}
			}
			
			clr.rewind();
		}
		
		private void scale(float[][] data,float scale){
			float ave=0;
			
			max=Float.MIN_VALUE;	min=Float.MAX_VALUE;
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float tmp=data[j][i];
				
				ave+=tmp;
				if(tmp>max) max=tmp;
				if(tmp<min) min=tmp;
			}
			
			ave/=(x*y);
			
			float de=(max-min)*scale;
			
			if(de<1e-11f) return;
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				data[j][i]-=ave;
				data[j][i]/=de;
			}
			
			min-=ave;	max-=ave;
			min/=de;	max/=de;
		}
	}
	
	
	/** test*/
	public static void main(String args[]){
		int len=101;
		float[][] data=new float[len][len];
		
		VisualizeDataInJOGL de=new VisualizeDataInJOGL(data);
		
		for(int l=0;l<1000;l++){
			for(int j=0;j<len;j++)
			for(int i=0;i<len;i++)
			data[j][i]=(float)(Math.sin(i*Math.cos(Math.toRadians(l))/3f)*Math.cos(j/6f))/2;
			
			de.setData(data);
			
			try{ TimeUnit.MILLISECONDS.sleep(40);}
			catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
		}
    }
}
