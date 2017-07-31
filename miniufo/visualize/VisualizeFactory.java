/**
 * @(#)VisualizeFactory.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.descriptor.DataDescriptor;


/**
 * factory method for visualization
 */
public final class VisualizeFactory{
	//
	private boolean loop=false;
	
	private DataDescriptor dd=null;
	
	
	/**
	 * constructor
	 */
	public VisualizeFactory(DataDescriptor dd){ this.dd=dd;}
	
	
	/**
	 * display a horizontal field in a variable specified by t and z
	 * t and z start from 0
	 * 
	 * @param	v	variable contains a horizontal field
	 * @param	t	specify the t dimension
	 * @param	z	specify the z dimension
	 */
	public void display(Variable v,int t,int z){
		int xcount=v.getXCount();
		int ycount=v.getYCount();
		
		float[][] data=new float[ycount][xcount];
		
		// prepare the data for display
		if(v.isTFirst()){
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++)
			data[ycount-j-1][xcount-i-1]=v.getData()[t][z][j][i];
			
		}else{
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++)
			data[ycount-j-1][xcount-i-1]=v.getData()[z][j][i][t];
		}
		
		VisualizeData vd=new VisualizeData(data);
		setVisualRange(vd,v);
		vd.updateModel(30);
	}
	
	public void animate(final Variable v,final int z,final long interval){
		loop=true;
		
		final int xcount=v.getXCount();
		final int ycount=v.getYCount();
		
		final float[][] data=new float[ycount][xcount];
		
		final VisualizeData vd=new VisualizeData(data);
		
		// prepare the data for display
		if(v.isTFirst()){
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++)
			data[j][i]=v.getData()[0][z][j][i];
			
		}else{
			for(int j=0;j<ycount;j++)
			for(int i=0;i<xcount;i++)
			data[j][i]=v.getData()[z][j][i][0];
		}
		
		setVisualRange(vd,v);
		
		new Thread(){
			public void run(){
				while(loop){
					for(int l=0,L=v.getTCount();l<L;l++){
						// prepare the data for display
						if(v.isTFirst()){
							for(int j=0;j<ycount;j++)
							for(int i=0;i<xcount;i++)
							data[ycount-j-1][xcount-i-1]=v.getData()[l][z][j][i];
							
						}else{
							for(int j=0;j<ycount;j++)
							for(int i=0;i<xcount;i++)
							data[ycount-j-1][xcount-i-1]=v.getData()[z][j][i][l];
						}
						
						vd.updateModel(interval);
					}
					
					vd.updateModel(interval+200);
				}
			}
		}.start();
	}
	
	public void stopAnimation(){ loop=false;}
	
	
	/**
	 * set visualize range
	 * 
	 * @param	vd	VisualizeData to be set
	 * @param	r	range of the variable
	 */
	private void setVisualRange(VisualizeData vd,Variable v){
		Range range=v.getRange();
		SurfacePlotModel pm=vd.model;
		
		pm.xmax=dd.getXDef().getSamples()[range.getXRange()[0]-1];
		pm.xmin=dd.getXDef().getSamples()[range.getXRange()[1]-1];
		
		pm.ymax=dd.getYDef().getSamples()[range.getYRange()[0]-1];
		pm.ymin=dd.getYDef().getSamples()[range.getYRange()[1]-1];
		
		float zmax=ArrayUtil.getMax(v.getData());
		float zmin=ArrayUtil.getMin(v.getData());
		float zint=(zmax-zmin)/8;
		
		pm.zmax=zmax+zint;
		pm.zmin=zmin-zint;
	}
	
	
	/** test*/
	public static void main(String args[]){
		miniufo.diagnosis.DiagnosisFactory df=
		miniufo.diagnosis.DiagnosisFactory.parseFile("D:/Data/DiagnosisVortex/Haima/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable REFC=df.getVariables(new Range("lev(200,200);t(1,1)",dd),"h")[0];
		
		VisualizeFactory vf=new VisualizeFactory(dd);
		//vf.animate(REFC,0,100);
		
		//try{ Thread.sleep(20000);}
		//catch(Exception e){e.printStackTrace(); System.exit(0);}
		
		//vf.stopAnimation();
		REFC.setValue(1000000);
		//REFC.add3DDisturbance(180,90,0,-5000,10);
		vf.display(REFC,0,0);
    }
}
