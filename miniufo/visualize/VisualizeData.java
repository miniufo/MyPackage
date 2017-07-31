/**
 * @(#)VisualizeData.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import miniufo.basic.ArrayUtil;


/**
 * visualize an array of 2-D discrete data
 */
public final class VisualizeData extends VisualizeFrame{
	//
	private static final long serialVersionUID = 4412909891321313417L;
	
	
	/**
	 * constructor
	 */
	public VisualizeData(float[][] data){
		super();
		
		model=new SurfacePlotModel();
		
		setData(data);
		
		setVisible(true);
	}
	
	
	public void setData(float[][] data){
		model.data=data;
		
		model.calcDivisionY=data.length-1;
		model.calcDivisionX=data[0].length-1;
		
		float zmax=ArrayUtil.getMax(data);
		float zmin=ArrayUtil.getMin(data);
		float zint=(zmax-zmin)/6;
		
		if(model.zmax<zmax){
			model.zmax=zmax+zint;
			model.zmin=zmin-zint;
		}
		
		if(model.zmin>zmin){
			model.zmax=zmax+zint;
			model.zmin=zmin-zint;
		}
	}
	
	
	/** test
	public static void main(String args[]){
		int len=31;
		float[][] data=new float[len][len];
		
		VisualizeData de=new VisualizeData(data);
		
		for(int l=0;l<1000;l++){
			for(int j=0;j<len;j++)
			for(int i=0;i<len;i++)
			data[j][i]=(float)(Math.sin(i*Math.cos(Math.toRadians(l))/3f)*Math.cos(j/6f))/2;
			
			de.updateModel(50);
		}
    }*/
}
