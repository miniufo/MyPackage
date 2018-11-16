/**
 * @(#)ReferenceStat1D.java	1.0 2018.04.08
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.mathsphysics.TridiagonalAlg;
import miniufo.statistics.StatisticsUtil;


/**
 * Calculate 1D reference state for barotropic flow or single-layer isentropic flow.
 * Reference: Thuburn and Lagneau (1999, JAS)
 *
 * @version 1.0, 2018.04.08
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ReferenceState1D{
	//
	private Variable Mref =null;	//     mass    in the transformed contour-based space
	private Variable Cref =null;	// circulation in the transformed contour-based space
	private Variable Qref =null;	//  PV values  in the transformed contour-based space
	private Variable Uref =null;	// zonal wind  in the transformed contour-based space
	private Variable Href =null;	//  thickness  in the transformed contour-based space
	
	private Contours[][] cntrs=null;
	private Params      params=null;
	
	private ContourSpatialModel csm=null;	// defines contours
	
	
	/**
	 * Constructor
	 * 
	 * @param	csm		contour spatial model that defines contours
	 * @param	sm		spatial model that defines transformed Y-space
	 */
	public ReferenceState1D(ContourSpatialModel csm,int ygrid){
		this.csm=csm;
		
		params=new Params(ygrid);
	}
	
	
	/**
	 * Initialized by PV contours, mass per area, and absolute vorticity.
	 * 
	 * @param	PVXY		PV in original coordinate (m^2 K s^-1 kg^-1)
	 * @param	rhoXY		mass per area in original coordinate (kg m^-2)
	 * @param	zetaXY		absolute vorticity in original coordinate (s^-1)
	 * @param	numOfC		number of contours used to in contour coordinate
     * @param	resRatio	resolution ratio (>=1) for on-the-fly interpolation
	 * @param	adjCtr		adjust contours or not
	 */
	public void initializeByPV(Variable PVXY,Variable rhoXY,Variable zetaXY,int numOfC,int resRatio,boolean adjCtr){
		csm.initContourByTracer(PVXY,numOfC,resRatio,true,adjCtr);
		csm.getAreasBoundedByContour();
		csm.cContourEnclosedMass(rhoXY);
		csm.cContourEnclosedCirculation(zetaXY);
		
		cntrs=csm.getContours();
		
		if(cntrs==null)
		throw new IllegalArgumentException("ContourSpatialModel is not initialized by tracer");
	}
	
	
	/**
	 * Solve to get the reference state of mass and circulation.
	 */
	public void solve(){
		int t=csm.dd.getTCount(),z=csm.dd.getZCount(),y=params.yc,x=1;
		
		Mref=new Variable("Mref",new Range(t,z,y,x));
		Cref=new Variable("Cref",new Range(t,z,y,x));
		Qref=new Variable("Qref",new Range(t,z,y,x));
		Uref=new Variable("Uref",new Range(t,z,y,x));
		Href=new Variable("Href",new Range(t,z,y,x));
		
		float[][][][] mdata=Mref.getData();
		float[][][][] cdata=Cref.getData();
		float[][][][] qdata=Qref.getData();
		float[][][][] udata=Uref.getData();
		float[][][][] hdata=Href.getData();
		
		if(Mref.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				Contours cs=cntrs[k][l];
				
				double[] tmp=cs.getMappedMass();
				
				tmp[tmp.length-1]=0;// ensure minimum mass is zero
				
				double[] M=initializeM(y,ArrayUtil.getMax(cs.getMappedMass()));
				
				double acc0=StatisticsUtil.sum(M);
				
				// maximum loop is 200 and tolerance is 1e-13
				for(int i=0;i<200;i++){
					double[] prev=M.clone();
					params.iterateOnce(cs,M);
					double acc=StatisticsUtil.c1Norm(prev,M);
					
					if(acc/acc0<1e-13) break;
					if(i==199) System.out.println("maximum loop of 200 reached at t="+(l+1)+", z="+(k+1));
				}
				
				double[] Q=cs.findContoursByMass(M);
				double[] C=cs.findCirculations(Q);
				double[] U=params.cZonalWind(C);
				double[] h=params.cDensity(M);
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					mdata[l][k][j][i]=(float)M[j];
					cdata[l][k][j][i]=(float)C[j];
					qdata[l][k][j][i]=(float)Q[j];
					udata[l][k][j][i]=(float)U[j];
					hdata[l][k][j][i]=(float)h[j];
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				Contours cs=cntrs[k][l];
				
				double[] M=initializeM(y,ArrayUtil.getMax(cs.getMappedMass()));
				
				double acc0=StatisticsUtil.sum(M);
				
				// maximum loop is 200 and tolerance is 1e-13
				for(int i=0;i<200;i++){
					double[] prev=M.clone();
					params.iterateOnce(cs,M);
					double acc=StatisticsUtil.c1Norm(prev,M);
					
					if(acc/acc0<1e-13) break;
				}
				
				double[] Q=cs.findContoursByMass(M);
				double[] C=cs.findCirculations(Q);
				double[] U=params.cZonalWind(C);
				double[] h=params.cDensity(M);
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					mdata[k][j][i][l]=(float)M[j];
					cdata[k][j][i][l]=(float)C[j];
					qdata[k][j][i][l]=(float)Q[j];
					udata[k][j][i][l]=(float)U[j];
					hdata[k][j][i][l]=(float)h[j];
				}
			}
		}
	}
	
	
	/*** getor and setor ***/
	public Variable getMassRef(){ return Mref;}
	
	public Variable getCirculationRef(){ return Cref;}
	
	public Variable getPVRef(){ return Qref;}
	
	public Variable getURef(){ return Uref;}
	
	public Variable getHRef(){ return Href;}
	
	
	/*** helper methods and classes ***/
	private double[] initializeM(int ycount,double Mmax){
		double[] m=new double[ycount];
		
		double incY=180.0/(ycount-1.0);
		
		for(int j=0;j<ycount;j++){
			double tmp=Math.sin(Math.toRadians(-90.0+j*incY));
			m[j]=(1.0-tmp)*Mmax/2.0;
		}
		
		return m;
	}
	
	
	private final class Params{
		//
		private static final double twoPI =2.0*Math.PI;
		private static final double omega2=SpatialModel.omegaEarth*SpatialModel.omegaEarth;
		
		private int yc=0;			// grids in the Y-dimension
		
		private double delY  =0;	// Y-increment (m)
		private double delSqr=0;	// squared increment
		
		private double[] Ys  =null;	// values of Ys (degree)
		private double[] cos =null;	// cosine of latitude at half grid
		private double[] acos=null;	// a * cos(Lat)
		private double[] asg =null;	// a * sin(Lat) / g
		
		private TridiagonalAlg solver=null;
		
		
		/**
		 * Constructor
		 * 
		 * @param	ycount	ycount from South Pole to North Pole
		 * @param	cs		contours
		 */
		public Params(int ycount){
			this.yc=ycount;
			this.Ys    =new double[ycount];
			this.cos   =new double[ycount];
			this.acos  =new double[ycount];
			this.asg   =new double[ycount];
			
			solver=new TridiagonalAlg(ycount);
			
			double incY=180.0/(ycount-1.0);
			
			for(int j=0;j<ycount;j++){
				  Ys[j]=-90+j*incY;
				 cos[j]=Math.cos(Math.toRadians(Ys[j]+incY/2.0));
				acos[j]=Math.cos(Math.toRadians(Ys[j]))*SpatialModel.REarth;
				 asg[j]=SpatialModel.REarth*Math.sin(Math.toRadians(Ys[j]))/SpatialModel.gEarth;
			}
			
			delY=Math.toRadians(incY)*SpatialModel.REarth;
			delSqr=delY*delY;
		}
		
		
		public void iterateOnce(Contours cs,double[] M){
			double[] Q=cs.findContoursByMass(M);
			double[] C=cs.findCirculations(Q);
			
			double[] a=new double[yc-1];
			double[] b=new double[yc  ];
			double[] c=new double[yc-1];
			double[] d=new double[yc  ];
			
			for(int j=1,J=yc-1;j<J;j++){
				b[j  ]=-1.0/cos[j-1]-1.0/cos[j]-asg[j]/Math.PI/Math.pow(acos[j],3.0)*C[j]*Q[j]*delSqr;
				a[j-1]= 1.0/cos[j-1];
				c[j  ]= 1.0/cos[j];
				d[j  ]=asg[j]*((C[j]*C[j]-2.0*C[j]*Q[j]*M[j])/twoPI/Math.pow(acos[j],3.0)-twoPI*omega2*acos[j])*delSqr;
			}
			
			b[0]=1;    b[yc-1]=1;		// fixed boundary conditions
			d[0]=M[0]; d[yc-1]=M[yc-1];	// fixed boundary conditions
			
			solver.trace(a,b,c,d,M);
		}
		
		public double[] cDensity(double[] M){
			double[] re=new double[yc];
			
			re[0]=-(M[1]-M[0])/delY/twoPI/((acos[0]+acos[1])/2.0);
			for(int j=1;j<yc-1;j++) re[j]=-(M[j+1]-M[j-1])/(2.0*delY)/twoPI/acos[j];
			re[yc-1]=-(M[yc-1]-M[yc-2])/delY/twoPI/((acos[yc-1]+acos[yc-2])/2.0);
			
			return re;
		}
		
		public double[] cZonalWind(double[] C){
			double[] re=new double[yc];
			
			re[0]=0; re[yc-1]=0;
			for(int j=1;j<yc-1;j++) re[j]=C[j]/twoPI/acos[j]-SpatialModel.omegaEarth*acos[j];
			
			return re;
		}
	}
	
	
	/*** test
	public static void main(String[] args){
		
	}*/
}
