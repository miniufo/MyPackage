/**
 * @(#)DiffusionModel1D.java	1.0 2018.07.31
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.prognosticModel;

import java.io.FileWriter;
import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.mathsphysics.TridiagonalAlg;


/**
 * 1D diffusion model
 *
 * @version 1.0, 2018.07.31
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DiffusionModel1D{
	//
	private int t;				// timesteps to be integrated, excluding initial state
	private int y;				// y-grid points
	
	private float deltaT=0;		// time step used to forward the model
	private float deltaY=0;		// y-grid spacing (m)
	private float ratio =0;		// deltaT/(2*deltaY*deltaY)
	private float trMass=0;		// total tracer mass
	private float widVar=0;		// variance of initial width of Gaussian tracer distribution (m^2)
	private float y0    =0;		// initial y-position of Gaussian tracer distribution (m)
	
	private float[][] k=null;	// t*y, diffusivity (m^2 s^-1)
	private float[][] q=null;	// 2*y, tracer (%)
	
	
	/**
     * Constructor.
     *
     * @param	t	how many timesteps to be integrated, excluding initial state
     * @param	dt	delta t in second for integration (second)
     * @param	y	how many y-grid points
     * @param	dy	delta y (meter)
     */
	public DiffusionModel1D(int t,int dt,int y,float dy){
		this.t=t;	this.deltaT=dt;
		this.y=y;	this.deltaY=dy;
		
		k=new float[t+1][y];
		q=new float[2][y];
		
		ratio=deltaT/(2f*deltaY*deltaY);
	}
	
	
	/**
     * Initialize tracer distribution using a prescribed tracer variable.
     * Only the first time record and y-direction are used.
     *
     * @param	tr	tracer distribution
     */
	public void initializeTracer(Variable tr){
		if(tr.getYCount()!=y) throw new IllegalArgumentException("y-count not the same");
		
		float[][][][] tdata=tr.getData();
		
		if(tr.isTFirst())
			for(int j=0;j<y;j++) q[0][j]=tdata[0][0][j][0];
		else
			for(int j=0;j<y;j++) q[0][j]=tdata[0][j][0][0];
	}
	
	/**
     * Initialize Gaussian-like tracer distribution.
     *
     * @param	ypos	center of the Gaussian distribution (m)
     * @param	widVar	variance of y-direction spreading (m^2)
     */
	public void initializeGaussianTracer(float ypos,float widVar){
		this.widVar=widVar;
		this.y0=ypos;
		
		for(int j=0;j<y;j++){
			float ydef=deltaY*j;
			double tmp1=ydef-ypos; tmp1*=tmp1;
			double tmp2=2.0*widVar;
			q[0][j]=(float)(100.0*Math.exp(-tmp1/tmp2)/Math.sqrt(Math.PI*tmp2));
		}
	}
	
	/**
     * Initialize diffusivity using a prescribed variable.
     *
     * @param	kappa	diffusivity variable (m^2 s^-1)
     */
	public void initializeDiffusivity(Variable kappa){
		int tk=kappa.getTCount();
		
		if(kappa.getYCount()!=y) throw new IllegalArgumentException("y-count not the same");
		if(tk!=t+1) throw new IllegalArgumentException("t-count of kappa ("+tk+") is not the same as t+1 ("+(t+1)+")");
		
		float[][][][] kdata=kappa.getData();
		
		if(kappa.isTFirst())
			for(int j=0;j<y;j++)
			for(int l=0;l<tk;l++) k[l][j]=kdata[l][0][j][0];
		else
			for(int j=0;j<y;j++)
			for(int l=0;l<tk;l++) k[l][j]=kdata[0][j][0][l];
	}
	
	/**
     * Initialize diffusivity using a prescribed constant.
     *
     * @param	kappa	diffusivity variable (m^2 s^-1)
     */
	public void initializeDiffusivity(float kappa){
		for(int l=0;l<=t;l++)
		for(int j=0;j< y;j++) k[l][j]=kappa;
	}
	
	
	/**
     * Integrate the 1D diffusion model forward.
     *
     * @param	interval	output frequency (s)
     * @param	output		output file name
     */
	public void integration(int interval,String output){
		if(interval%deltaT!=0) throw new IllegalArgumentException(
			"interval ("+interval+") cannot be divided by deltaT ("+deltaT+")"
		);
		
		float t0=checkStability(ArrayUtil.getMax(k));
		
		System.out.println("\nstart integration...");
		
		int tag1=1,tag2=0;
		
		float[][][][] odata=new float[1][1][1][];
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(output);
		
		TridiagonalAlg tri=new TridiagonalAlg(y);
		
		Variable out=new Variable("out",true,new Range(1,1,1,1));
		
		// write initial condition
		odata[0][0][0]=q[0];
		out.setData(odata,true);
		cdws.writeData(out);
		cdws.writeData(out);
		
		for(int l=1;l<=t;l++){
			float[] qres=q[tag1];
			float[] qsrc=q[tag2];
			float[] anal=new float[y];
			
			/*************** predict q ***************/
			double[] a=new double[y-1];
			double[] b=new double[y  ];
			double[] c=new double[y-1];
			double[] d=new double[y  ];
			
			a[y-2]=-1;
			for(int j=0,J=y-2;j<J;j++){
				a[j]=-ratio*(k[l][j]+k[l][j+1])/2f;
			}
			
			b[0]=1; b[y-1]=1;
			for(int j=1,J=y-1;j<J;j++){
				b[j]=1+ratio*(k[l][j-1]+k[l][j])/2f+ratio*(k[l][j]+k[l][j+1])/2f;
			}
			
			c[0]=-1;
			for(int j=1,J=y-1;j<J;j++){
				c[j]=-ratio*(k[l][j]+k[l][j+1])/2f;
			}
			
			d[0]=0; d[y-1]=0;
			for(int j=1,J=y-1;j<J;j++){
				d[j]=ratio*(k[l-1][j-1]+k[l-1][j])/2f*qsrc[j-1]+
					 ratio*(k[l-1][j]+k[l-1][j+1])/2f*qsrc[j+1]+
					(1-ratio*(k[l-1][j-1]+k[l-1][j])/2f-ratio*(k[l-1][j]+k[l-1][j+1])/2f)*qsrc[j];
			}
			
			double[] re=tri.trace(a,b,c,d);
			
			float tt=l*deltaT;
			for(int j=0;j<y;j++){
				float ydef=deltaY*j;
				double tmp1=ydef-y0; tmp1*=tmp1;
				
				qres[j]=(float)re[j];
				anal[j]=(float)(100.0*
					Math.sqrt(1.0/((2.0*Math.PI)*(widVar+2.0*k[l][j]*tt)))*Math.exp(-(tmp1)/(4.0*k[l][j]*(t0+tt)))
				);
			}
			
			// output q
			if(l%(interval/deltaT)==0){
				cTotalTracer(tag1);
				
				System.out.println(
					"  output "+l+" step, total tracer mass is: "+trMass
				);
				
				odata[0][0][0]=q[tag1];	cdws.writeData(out);
				odata[0][0][0]=anal;	cdws.writeData(out);
			}
			
			// swap tag1 and tag2
			int tmp=tag1;	tag1=tag2;	tag2=tmp;
		}
		
		cdws.closeFile();
		
		
		// write ctl file
		StringBuffer sb=new StringBuffer();
		
		float dt=interval/60;
		
		sb.append("dset ");	sb.append(output);	sb.append("\n");
		sb.append("undef -9999.0\n");
		sb.append("title 1D diffusion model output\n");
		sb.append("xdef   1 linear 0 0.05\n");
		sb.append("ydef "+y+" linear -10 0.05\n");
		sb.append("zdef   1 levels 1\n");
		sb.append("tdef ");	sb.append(Math.round(t*deltaT/interval)+1);	sb.append(" linear 00:00Z16dec2002 ");
		
		if(dt<60)       { sb.append((int)(dt)     ); sb.append("mn\n");}
		else if(dt<1440){ sb.append((int)(dt/60)  ); sb.append("hr\n");}
		else            { sb.append((int)(dt/1440)); sb.append("dy\n");}
		
		sb.append("vars 2\n");
		sb.append("q  0 99 tracer concentration\n");
		sb.append("qa 0 99 analytical solution\n");
		sb.append("endvars\n");
		
		try(FileWriter fw=new FileWriter(output.split("\\.")[0]+".ctl")){
			fw.write(sb.toString());
		}catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
		
		System.out.println("finish integration.");
	}
	
	
	/*** getor and setor ***/
	public float[] getInitialTracer(){ return q[0];}
	
	
	/*** helper methods ***/
	private void cTotalTracer(int tag){
		double tmp=0;
		
		float[] tmpTr=q[tag];
		
		for(int i=0;i<y;i++) tmp+=tmpTr[i];
		
		trMass=(float)(tmp*deltaY);
	}
	
	private float checkStability(float maxKappa){
		float criterion=ratio*2f*maxKappa;
		
		double t0=widVar/(2.0*maxKappa);
		
		System.out.println("Beta-plane diffusion model integration parameters:");
		System.out.println("  criterion: "+criterion+" (linear), should be < 0.5");
		System.out.println("  t0       : "+t0+" (constant diffusivity)");
		System.out.println("  delta Y  : "+deltaY/1000f+" km");
		
		return (float)t0;
	}
	
	
	/** test
	public static void main(String[] args){
		int steps=240000;
		int deltaT=3600;	// 1 hr
		int ycount=401 ;
		int deltaY=5000;	// 5 km
		
		DiffusionModel1D diffModel=new DiffusionModel1D(steps,deltaT,ycount,deltaY);
		
		diffModel.initializeGaussianTracer(1000e3f,2f*10000f*86400f);
		diffModel.initializeDiffusivity(10000);
		
		diffModel.integration(864000,"/public/home/qianyk/Data/DiffModel.dat");
	}*/
}
