/**
 * @(#)ReferenceStat2D.java	1.0 2018.09.11
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import java.util.Arrays;
import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.geophysics.atmos.Dynamics;
import miniufo.geophysics.atmos.ThermoDynamics;
import miniufo.mathsphysics.Integrator1D;


/**
 * Calculate 2D reference state for multi-layer isentropic flow.
 *
 * @version 1.0, 2018.09.11
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ReferenceState2D{
	//
	private int YGrid=0;
	
	private Variable pm  =null;
	
	private Variable Qref=null;	//  PV values  in the transformed contour-based space
	private Variable Uref=null;	// zonal wind  in the transformed contour-based space
	private Variable Pref=null;	//  pressure   in the transformed contour-based space
	private Variable Dref=null;	//   density   in the transformed contour-based space
	private Variable Mref=null;	//     mass    in the transformed contour-based space
	private Variable Cref=null;	// circulation in the transformed contour-based space
	
	private Contours[][] cntrs=null;
	private Params      params=null;
	
	private ContourSpatialModel csm=null;	// defines contours
	
	
	/**
	 * Constructor
	 * 
	 * @param	csm		contour spatial model that defines contours
	 * @param	YGrid	Y-Grid count in the contour space
	 */
	public ReferenceState2D(ContourSpatialModel csm,int YGrid){
		this.csm=csm;
		this.YGrid=YGrid;
	}
	
	
	/**
	 * Initialized by PV contours, mass per area, and absolute vorticity.
	 * 
	 * @param	PVXY		PV in original coordinate (m^2 K s^-1 kg^-1)
	 * @param	rhoXY		mass per area in original coordinate (kg m^-2)
	 * @param	zetaXY		absolute vorticity in original coordinate (s^-1)
	 * @param	pm			zonal mean pressure for upper/lower boundary conditions (Pa)
	 * @param	numOfC		number of contours used to in contour coordinate
     * @param	resRatio	resolution ratio (>=1) for on-the-fly interpolation
	 * @param	adjCtr		adjust contours or not
	 */
	public void initializeByPV(Variable PVXY,Variable rhoXY,Variable zetaXY,Variable pm,int numOfC,int resRatio,boolean adjCtr){
		this.pm=pm;
		
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
		int t=pm.getTCount(),z=csm.dd.getZCount(),y=YGrid,x=1;
		
		Qref=new Variable("Qref",false,new Range(t,z,y,x));
		Uref=new Variable("Uref",false,new Range(t,z,y,x));
		Pref=new Variable("Pref",false,new Range(t,z,y,x));
		Dref=new Variable("Dref",false,new Range(t,z,y,x));
		Mref=new Variable("Mref",false,new Range(t,z,y,x));
		Cref=new Variable("Cref",false,new Range(t,z,y,x));
		
		float[][][][] qdata=Qref.getData();
		float[][][][] udata=Uref.getData();
		float[][][][] pdata=Pref.getData();
		float[][][][] ddata=Dref.getData();
		float[][][][] mdata=Mref.getData();
		float[][][][] cdata=Cref.getData();
		
		for(int l=0;l<t;l++){
			Contours[] cs=new Contours[z];
			
			double[] pTop=new double[y];
			double[] pSfc=new double[y];
			
			for(int k=0;k<z;k++) cs[k]=cntrs[k][l];
			for(int j=0;j<y;j++){
				pTop[j]=pm.getData()[l][z-1][j][0];
				pSfc[j]=pm.getData()[l][  0][j][0];
			}
			
			params=new Params(z,y,csm.dd.getZDef().getIncrements()[0],pTop,pSfc,cs);
			params.firstGuess();
			
			for(int i=0;i<8;i++) params.outerLoop();
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				qdata[k][j][i][l]=(float)params.Q [k][j];
				udata[k][j][i][l]=(float)params.U [k][j];
				pdata[k][j][i][l]=(float)params.P [k][j];
				ddata[k][j][i][l]=(float)params.sg[k][j];
				mdata[k][j][i][l]=(float)params.M2[k][j];
				cdata[k][j][i][l]=(float)params.C2[k][j];
			}
		}
	}
	
	
	/*** getor and setor ***/
	public Variable getMassRef(){ return Mref;}
	
	public Variable getCirculationRef(){ return Cref;}
	
	public Variable getPVRef(){ return Qref;}
	
	public Variable getURef(){ return Uref;}
	
	public Variable getPressureRef(){ return Pref;}
	
	public Variable getDensityRef(){ return Dref;}
	
	public Variable cMontgomeryStreamfunction(){
		int t=Pref.getTCount(),z=Pref.getZCount(),y=Pref.getYCount();
		
		Variable mont=new Variable("mont",Pref);
		
		float[]       theta=csm.dd.getZDef().getSamples();
		float[][][][] mdata=mont.getData();
		float[][][][] pdata=Pref.getData();
		
		for(int l=0;l<t;l++)
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++){
			float T=ThermoDynamics.cTemperature(theta[k],pdata[k][j][0][l]);
			float Z=Dynamics.cGeopotentialHeight(mdata[k][j][0][l],T);
			mdata[k][j][0][l]=Dynamics.cMontgomeryStreamfunction(Z,T);
		}
		
		return mont;
	}
	
	public Variable cGeopotential(){
		int t=Pref.getTCount(),z=Pref.getZCount(),y=Pref.getYCount();
		
		Variable geop=new Variable("geop",Pref);
		
		float[]       theta=csm.dd.getZDef().getSamples();
		float[][][][] gdata=geop.getData();
		float[][][][] pdata=Pref.getData();
		
		for(int l=0;l<t;l++)
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++){
			float T=ThermoDynamics.cTemperature(theta[k],pdata[k][j][0][l]);
			gdata[k][j][0][l]=Dynamics.cGeopotentialHeight(gdata[k][j][0][l],T);
		}
		
		return geop;
	}
	
	
	/*** helper methods and classes ***/
	
	
	private static final class Params{
		//
		private static final double twoPI  =2.0*Math.PI;
		private static final double twoPIa2=2.0*Math.PI*SpatialModel.REarth*SpatialModel.REarth;
		private static final double PIRg   =Math.PI*ThermoDynamics.Rd*SpatialModel.gEarth;
		private static final double p0Kappa=Math.pow(ThermoDynamics.Pref,ThermoDynamics.kappa);
		
		private int mxLoop   =600;	// max loop of the iteration
		private int yc       =0;	// grids in the Y-dimension
		private int zc       =0;	// grids in the z-dimension
		
		private double delY  =0;	// dY (1)
		private double delZ  =0;	// dTheta (K)
		private double delSqr=0;	// (dY/dTheta)^2
		private double optArg=0;	// optimal argument for SOR
		private double tol   =1e-12f;// tolerance for the iteration
		
		private double[] Ys  =null;	// values of Ys (degree)
		private double[] Zs  =null;	// values of thetas (K)
		private double[] sinQ=null;	// abs(sin(Lat))^0.25
		private double[] sin =null;	// sin(lat)
		private double[] cos =null;	// cos(Lat)
		private double[] fCor=null;	// 2 * Omega * sin(lat)
		
		private double[]   pTop=null;	// pressure at the top boundary (Pa)
		private double[]   pSfc=null;	// pressure at the surface boundary (Pa)
		
		private double[][] ud=null;	// update of Yeq
		private double[][] Ym=null;	// running mean of Yeq
		private double[][] Ye=null;	// equivalent Ys (degree)
		private double[][] Q =null;	// PV contours
		private double[][] M2=null;	// mass of 2D
		private double[][] C2=null;	// circulation of 2D
		private double[][] M3=null;	// mass of 3D
		private double[][] C3=null;	// circulation of 3D
		private double[][] P =null;	// pressure
		private double[][] U =null;	// zonal wind
		private double[][] sg=null;	// isentropic density
		
		private Contours[] cs=null;
		
		
		/**
		 * Constructor.
		 * 
		 * @param	zc		z-count
		 * @param	yc		y-count from South Pole to North Pole
		 * @param	delZ	vertical interval (K)
		 * @param	pTop	pressure at the top boundary
		 * @param	pSfc	pressure at the surface boundary
		 * @param	cs		contours at a time instant
		 */
		public Params(int zc,int yc,float delZ,double[] pTop,double[] pSfc,Contours[] cs){
			this.yc=yc;
			this.zc=zc; this.pTop=pTop;
			this.cs=cs; this.pSfc=pSfc;
			
			Ys  =new double[yc];		Zs  =new double[zc];
			cos =new double[yc];		sin =new double[yc];
			fCor=new double[yc];		sinQ=new double[yc];
			Ym  =new double[zc][yc];
			ud  =new double[zc][yc];	sg  =new double[zc][yc];
			Ye  =new double[zc][yc];	Q   =new double[zc][yc];
			M2  =new double[zc][yc];	C2  =new double[zc][yc];
			M3  =new double[zc][yc];	C3  =new double[zc][yc];
			P   =new double[zc][yc];	U   =new double[zc][yc];
			
			for(int j=0;j<yc;j++){
				  Ys[j]=-1.0+2.0*j/(yc-1.0);	// sin(lat)
				 sin[j]=Ys[j];
				fCor[j]=2.0*SpatialModel.omegaEarth*sin[j];
				 cos[j]=Math.sqrt(1.0-sin[j]*sin[j]);
				
				double lat=Math.toDegrees(Math.asin(Ys[j]));
				double factor=1;
				if(Math.abs(lat)<20) factor=0.1+0.9*Math.abs(lat)/20.0;
				sinQ[j]=Math.pow(Math.abs(sin[j]),0.7)*Math.pow(factor,2.0);
			}
			
			for(int k=0;k<zc;k++){
				Zs[k]=265.0+15.0*k;
				
				for(int j=0;j<yc;j++){ Ye[k][j]=Ys[j]; Ym[k][j]=Ys[j];}
				
				double[] tmp=cs[k].getMappedEquivalentYs();
				
				for(int j=0;j<tmp.length;j++) tmp[j]=Math.sin(Math.toRadians(tmp[j]));
			}
			
			delY=2.0/(yc-1.0);
			this.delZ=delZ;
			delSqr=delY/delZ; delSqr*=delSqr;
			
			double eps=Math.pow(Math.sin(Math.PI/(2.0*yc+2.0)),2.0)+Math.pow(Math.sin(Math.PI/(2.0*zc+2.0)),2.0);
			optArg=2.0/(1.0+Math.sqrt(eps*(2.0-eps)));
		}
		
		
		public void firstGuess(){
			mapToRegularYs();
			
			for(int k=0;k<zc;k++){
				System.arraycopy(M2[k],0,M3[k],0,yc);
				System.arraycopy(C2[k],0,C3[k],0,yc);
				
				cs[k].setContours(Q[k]);
				cs[k].setYEs(Ye[k]);
				cs[k].setMass(M2[k]);
				cs[k].setCirculation(C2[k]);
			}
			
			cDensity();
			cPressure();
			cZonalWind();
		}
		
		public void outerLoop(){
			mapToRegularYs();
			for(int k=0;k<zc;k++){
				System.arraycopy(M2[k],0,M3[k],0,yc);
				System.arraycopy(C2[k],0,C3[k],0,yc);
			}
			inversion2();
			cDensity();
			cPressure();
			cZonalWind();
			correctCT();
		}
		
		public void mapToRegularYs(){
			for(int k=0;k<zc;k++){
				double[] tmp;
				tmp=cs[k].findContoursByYeq(Ye[k]);	System.arraycopy(tmp,0, Q[k],0,yc);
				tmp=cs[k].findMasses(Q[k]);			System.arraycopy(tmp,0,M2[k],0,yc);
				tmp=cs[k].findCirculations(Q[k]);	System.arraycopy(tmp,0,C2[k],0,yc);
			}
		}
		
		public void inversion(){
			int loop=0;
			double convSpd=0;	// convergent speed = delta norm(S) / norm(S)
			double normPrev=Double.MAX_VALUE;
			
			double[][] dM=new double[zc][yc];
			
			//System.out.println("optArg: "+optArg);//1.8958787389964686
			optArg=1.9252;
			
			do{
				//for(int k=1,K=zc-1;k<K;k++)
				for(int k=zc-2;k>=1;k--)
				for(int j=1,J=yc-1;j<J;j++){
					double a=Math.pow((P[k][j-1]+P[k][j  ])/2.0,ThermoDynamics.kappa-1.0);
					double c=Math.pow((P[k][j  ]+P[k][j+1])/2.0,ThermoDynamics.kappa-1.0);
					double b=-(a+c);
					
					double coeff=delSqr*p0Kappa/PIRg*sin[j]/Math.pow(cos[j],4.0);
					
					double d=coeff*C2[k-1][j]*Q[k-1][j];
					double e=coeff*C2[k  ][j]*Q[k  ][j]*-2.0;
					double f=coeff*C2[k+1][j]*Q[k+1][j];
					
					double g=-a*M2[k][j-1]-b*M2[k][j]-c*M2[k][j+1]-coeff*(C2[k-1][j]*C2[k-1][j]-2.0*C2[k][j]*C2[k][j]+C2[k+1][j]*C2[k+1][j])/2.0;
					
					double res=a*dM[k][j-1]+b*dM[k][j]+c*dM[k][j+1]+d*dM[k-1][j]+e*dM[k][j]+f*dM[k+1][j]-g;
					
					dM[k][j]+=optArg*res/(-b-e);
				}
				
				double norm=cAbsMean(dM);
				
				if(Double.isNaN(norm)||norm>1e18){
					System.out.println("overflow: "+norm);
					break;
				}
				
				convSpd=Math.abs(norm-normPrev)/normPrev;
				
				if(convSpd<tol||loop>mxLoop) break;
				
				normPrev=norm; loop++;
				
				//System.out.println(loop+"\t"+norm+"\t"+convSpd);
				
			}while(true);
			
			System.out.println("loop: "+loop);
			
			for(int k=0;k<zc;k++){
				for(int j=1;j<yc-1;j++){
					M2[k][j]=M3[k][j]+dM[k][j];
					C2[k][j]=C3[k][j]+Q[k][j]*dM[k][j];
				}
				
				int idx=ArrayUtil.nonMonoDecreIdx(M2[k]);
				if(idx!=-1) throw new IllegalArgumentException("[k,j]: ["+k+", "+idx+"] "+Arrays.toString(M2[k]));
			}
		}
		
		public void inversion2(){
			int loop=0;
			double convSpd=0;	// convergent speed = delta norm(S) / norm(S)
			double normPrev=Double.MAX_VALUE;
			
			double[][] dM=new double[zc][yc];
			
			for(int k=0;k<zc;k++)
			for(int j=0;j<yc;j++) dM[k][j]=M2[k][j];
			
			//System.out.println("optArg: "+optArg);//1.8958787389964686
			optArg=1.9252;
			
			do{
				for(int k=1,K=zc-1;k<K;k++)
				for(int j=1,J=yc-1;j<J;j++){
					double a=Math.pow((P[k][j-1]+P[k][j  ])/2.0,ThermoDynamics.kappa-1.0);
					double c=Math.pow((P[k][j  ]+P[k][j+1])/2.0,ThermoDynamics.kappa-1.0);
					double b=-(a+c);
					
					double coeff=delSqr*p0Kappa/PIRg*sin[j]/Math.pow(cos[j],4.0);
					
					double d=coeff*C2[k-1][j]*Q[k-1][j];
					double e=coeff*C2[k  ][j]*Q[k  ][j]*-2.0;
					double f=coeff*C2[k+1][j]*Q[k+1][j];
					
					double g=(d*M2[k-1][j]+e*M2[k][j]+f*M2[k+1][j])-coeff*(C2[k-1][j]*C2[k-1][j]-2.0*C2[k][j]*C2[k][j]+C2[k+1][j]*C2[k+1][j])/2.0;
					
					double res=a*dM[k][j-1]+b*dM[k][j]+c*dM[k][j+1]+d*dM[k-1][j]+e*dM[k][j]+f*dM[k+1][j]-g;
					
					dM[k][j]+=optArg*res/(-b-e);
				}
				
				double norm=cAbsMean(dM);
				
				if(Double.isNaN(norm)||norm>1e18){
					System.out.println("overflow: "+norm);
					break;
				}
				
				convSpd=Math.abs(norm-normPrev)/normPrev;
				
				if(convSpd<tol||loop>mxLoop) break;
				
				normPrev=norm; loop++;
				
				//System.out.println(loop+"\t"+norm+"\t"+convSpd);
				
			}while(true);
			
			System.out.println("loop: "+loop);
			
			for(int k=0;k<zc;k++){
				for(int j=1;j<yc-1;j++){
					double delM=dM[k][j]-M2[k][j];
					M2[k][j]=dM[k][j];
					C2[k][j]=C3[k][j]+Q[k][j]*delM;
				}
				
				int idx=ArrayUtil.nonMonoDecreIdx(M2[k]);
				if(idx!=-1) throw new IllegalArgumentException("[k,j]: ["+k+", "+idx+"] "+Arrays.toString(M2[k]));
			}
		}
		
		public void cDensity(){
			for(int k=0;k<zc;k++)
			for(int j=0;j<yc-1;j++) sg[k][j]=-(M2[k][j+1]-M2[k][j])/delY/twoPIa2;
		}
		
		public void cPressure(){
			double[] integrand=new double[zc];
			double[] tmp=null;
			
			for(int k=0;k<zc;k++) integrand[k]=sg[k][0]*SpatialModel.gEarth;
			tmp=Integrator1D.integrateBackward(-delZ,integrand,pTop[0]);
			for(int k=0;k<zc;k++) P[k][0]=tmp[k];
			
			for(int j=1;j<yc-1;j++){ // density is defined at half grid of Y
				for(int k=0;k<zc;k++) integrand[k]=(sg[k][j]+sg[k][j-1])/2.0*SpatialModel.gEarth;
				tmp=Integrator1D.integrateBackward(-delZ,integrand,pTop[j]);
				for(int k=0;k<zc;k++) P[k][j]=tmp[k];
			}
			
			for(int k=0;k<zc;k++) integrand[k]=sg[k][yc-2]*SpatialModel.gEarth;
			tmp=Integrator1D.integrateBackward(-delZ,integrand,pTop[yc-1]);
			for(int k=0;k<zc;k++) P[k][yc-1]=tmp[k];
		}
		
		public void cZonalWind(){
			for(int k=0;k<zc;k++)
			for(int j=1;j<yc-1;j++) U[k][j]=C2[k][j]/twoPI/cos[j]/SpatialModel.REarth-SpatialModel.omegaEarth*SpatialModel.REarth*cos[j];
		}
		
		public void correctCT(){
			final double Nr=4.0;
			
			for(int k=0;k<zc;k++){
				for(int j=1,J=yc-1;j<J;j++){
					double dYdM= (Ye[k][j+1]-Ye[k][j-1])/(M3[k][j+1]-M3[k][j-1]);
					double dYdC= (Ye[k][j+1]-Ye[k][j-1])/(C3[k][j+1]-C3[k][j-1]);
					double tmpM=-(M2[k][j]-M3[k][j])*dYdM*sinQ[j];
					double tmpC=-(C2[k][j]-C3[k][j])*dYdC*sinQ[j];
					double old =Ye[k][j];
					
					Ye[k][j]=(Ye[k][j]+ud[k][j]/2.0+(tmpM+tmpC)/4.0+Ym[k][j]/Nr)*Nr/(Nr+1);
					ud[k][j]= Ye[k][j]-old;
					Ym[k][j]=(Ym[k][j]*(Nr-1)+Ye[k][j])/Nr;
				}
				
				cs[k].setYEs(Ye[k]);
				
				int tmp=ArrayUtil.nonMonoIncreIdx(Ye[k]);
				if(tmp!=-1) throw new IllegalArgumentException("non-monotonic Ye["+k+","+(int)tmp+"]: "+Arrays.toString(Ye[k]));
			}
		}
		
		public double cAbsMean(double[][] data){
			double sum=0;
			
			for(int j=0,J=data.length;j<J;j++)
			for(int i=0,I=data[0].length;i<I;i++) sum+=Math.abs(data[j][i]);
			
			sum/=data.length*data[0].length;
			
			return sum;
		}
	}
	
	
	/*** test
	public static void main(String[] args){
		
	}*/
}
