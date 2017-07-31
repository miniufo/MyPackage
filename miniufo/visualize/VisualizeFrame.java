/**
 * @(#)VisualizeFrame.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import java.awt.Image;
import java.awt.FlowLayout;
import java.awt.BorderLayout;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.SwingUtilities;


/**
 * GUI to show the discrete 2-Dim data
 */
public abstract class VisualizeFrame extends javax.swing.JFrame{
	//
	private static final long serialVersionUID=-3359784694978970547L;
	
	protected boolean   finishRender=true;
	protected SurfacePlotModel model=null;
	protected SurfaceCanvas   canvas=null;
	
	
	/**
	 * constructor
	 */
	public VisualizeFrame(){
		getContentPane().setLayout(new BorderLayout());
		
		JPanel southPanel=new JPanel(new FlowLayout(FlowLayout.CENTER,25,5));
		
		southPanel.add(new JLabel("Rotate: Mouse Click & Drag"));
		southPanel.add(new JLabel("Zoom: Shift Key + Mouse Click & Drag"));
		southPanel.add(new JLabel("Move: Control Key + Mouse Click & Drag"));
		
		getContentPane().add(southPanel,BorderLayout.SOUTH);
		
		canvas=new SurfaceCanvas();
		getContentPane().add(canvas,BorderLayout.CENTER);
		
		setSize(700,500);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    }
	
	
	/*** getor and setor ***/
	public boolean finishRendering(){ return finishRender;}
	
	public Image getImage(){ return canvas.getImage();}
	
	public SurfacePlotModel getModel(){	return model;}
	
	public void setModel(SurfacePlotModel model){ this.model=model;}
	
	
	/**
	 * update the model data and then sleep for a while (sleep)
	 */
	public void updateModel(long sleep){
		while(!finishRender){
			try{Thread.sleep(5);}
			catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
		}
		
		finishRender=false;
		//new Thread(new Updater(sleep)).start();
		SwingUtilities.invokeLater(new Updater(sleep));
	}
	
	
	/**
	 * using another thread to update the model
	 */
	private class Updater implements Runnable{
		//
		private long sleeptime=0;
		
		public Updater(long t){ this.sleeptime=t;}
		
		//
		public void run(){
			canvas.setModel(model);
			try{Thread.sleep(sleeptime+5);}
			catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
			finishRender=true;
		}
	}
}
