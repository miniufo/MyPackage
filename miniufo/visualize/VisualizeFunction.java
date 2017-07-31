/**
 * @(#)VisualizeFunction.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;


/**
 * visualize a 2-D function
 */
public final class VisualizeFunction extends VisualizeFrame{
	//
	private static final long serialVersionUID = 4415434635925455964L;
	
	private boolean loop=false;
	
	private float[][] bckp=null;
	
	
	/**
	 * constructor
	 */
	public VisualizeFunction(Function2D func){
		super();
		
		model=new SurfacePlotModel();
		
		model.function=func;
		
		calculateAndUpdate();
		
		cloneData();
		
		setVisible(true);
	}
    
	
	/**
	 * change display domain
	 * 
	 * @param	xstr	start point of domain along x-direction
	 * @param	xend	end point of domain along x-direction
	 * @param	xcount	grid count along x-direction
	 * @param	ystr	start point of domain along y-direction
	 * @param	yend	end point of domain along y-direction
	 * @param	ycount	grid count along y-direction
	 */
	public void changeDomain(float xstr,float xend,int xcount,float ystr,float yend,int ycount){
		model.xmin=xstr;
		model.ymin=ystr;
		model.xmax=xend;
		model.ymax=yend;
		model.calcDivisionX=xcount;
		model.calcDivisionY=ycount;
		
		calculateAndUpdate();
		cloneData();
	}
	
	
	/**
	 * for animation
	 */
	public void animateAmplitude(final long interval){
		loop=true;
		
		new Thread(){
			public void run(){
				long phase=model.phase;
				
				float[][] data=model.data;
				
				while(loop){
					phase+=2;
					model.phase=phase;
					float amp=(float)Math.cos(Math.toRadians(phase));
					
					for(int j=0,J=data.length;j<J;j++)
					for(int i=0,I=data[0].length;i<I;i++)
					data[j][i]=bckp[j][i]*amp;
					
					updateModel(interval);
				}
			}
		}.start();
	}
	
	public void stopAnimation(){ loop=false;}
	
	
	/**
	 * calculate data and update model display
	 */
	private void calculateAndUpdate(){
		model.data=new float[model.calcDivisionY+1][model.calcDivisionX+1];
		
		float xstr=model.xmin,xinc=model.getXIncrement();
		float ystr=model.ymin,yinc=model.getYIncrement();
		
		float zmax=-Float.MAX_VALUE;
		float zmin= Float.MAX_VALUE;
		
		float[][] tmp=model.data;
		
		for(int j=0,J=model.calcDivisionY;j<=J;j++)
		for(int i=0,I=model.calcDivisionX;i<=I;i++){
			model.data[j][i]=model.calculateZ(xstr+i*xinc,ystr+j*yinc);
			
			if(tmp[j][i]>zmax) zmax=tmp[j][i];
			if(tmp[j][i]<zmin) zmin=tmp[j][i];
		}
		
		float interval=zmax-zmin;
		
		model.zmax=zmax+interval/5;
		model.zmin=zmin-interval/5;
		
		updateModel(100);
	}
	
	private void cloneData(){
		bckp=new float[model.data.length][];
		
		for(int i=0,I=model.data.length;i<I;i++)
		bckp[i]=model.data[i].clone();
	}
	
	
	/** test
	public static void main(String args[]){
		VisualizeFunction vf=new VisualizeFunction(
			new Function2D(){
				public float map(float x,float y){
					return (float)Math.sqrt(Math.hypot(x,y))/2;
				}
			}
		);
		
		vf.model.zmin=-1;
		
		vf.animateAmplitude(20);
		
		try{ Thread.sleep(5000);}
		catch(Exception e){e.printStackTrace(); System.exit(0);}
		
		vf.stopAnimation();
	}*/
}
