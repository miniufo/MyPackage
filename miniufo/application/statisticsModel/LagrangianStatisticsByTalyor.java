/**
 * @(#)LagrangianStatistics.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.util.List;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.Particle;
import miniufo.lagrangian.Record;
import miniufo.statistics.StatisticsUtil;


/**
 * Lagrangian statistics
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LagrangianStatisticsByTalyor extends SingleParticleStatistics{
	//
	private Variable[] vels=null;
	
	private static final BisectionSolver    BS=new BisectionSolver();
	private static final PrescribedFunction PF=new PrescribedFunction();
	
	
	/**
	 * constructor, the lengths of Particles in the list should be the same
	 */
	public LagrangianStatisticsByTalyor(List<? extends Particle> ls,DataDescriptor dd){
		super(ls,dd);
		
		vels=toVelocityVariables(ls);
	}
	
	
	/**
     * compute the Lagrangian statistics
     *
     * @param	maxlag	maximum lag time in the computation of lead-lag correlation
     * @param	dt		delta T of data interval (unit: s)
     */
	public Variable[] cStatisticsByTaylorTheory(int maxLag,int dt){
		Variable Ru=cAutoCorrelationFunction(vels[0],maxLag);
		Variable Rv=cAutoCorrelationFunction(vels[1],maxLag);
		
		Variable varu=StatisticsBasicAnalysisMethods.cTVariance(vels[0]); varu.setName("varu");
		Variable varv=StatisticsBasicAnalysisMethods.cTVariance(vels[1]); varv.setName("varv");
		
		Variable stdu=varu.sqrt();	stdu.setName("stdu");
		Variable stdv=varv.sqrt();	stdv.setName("stdv");
		
		Variable Tu=cTimescaleByIntegrateTo1stZeroCross(Ru,dt);	Tu.setName("Tu");
		Variable Tv=cTimescaleByIntegrateTo1stZeroCross(Rv,dt);	Tv.setName("Tv");
		
		Variable Lu=stdu.multiply(Tu);	Lu.setName("Lu");
		Variable Lv=stdv.multiply(Tv);	Lv.setName("Lv");
		
		Variable KHu=varu.multiply(Tu); KHu.setName("KHu");	KHu.setCommentAndUnit("eddy diffusive coefficient (m^2 s^-1)");
		Variable KHv=varv.multiply(Tv); KHv.setName("KHv");	KHv.setCommentAndUnit("eddy diffusive coefficient (m^2 s^-1)");
		
		stdu.multiplyEq(100);	// change unit to cm/s
		stdv.multiplyEq(100);
		
		Tu.divideEq(86400);		// change unit to day
		Tv.divideEq(86400);
		
		Lu.divideEq(1000);		// change unit to km
		Lv.divideEq(1000);
		
		KHu.divideEq(1000);		// change unit to 10^3 m^2/s
		KHv.divideEq(1000);
		
		return new Variable[]{varu,varv,stdu,stdv,Tu,Tv,Lu,Lv,KHu,KHv};
	}
	
	
	/**
     * remove the Lagrangian mean and trend
     */
	public void removeLagrangianMean(){
		vels[0].anomalizeT();
		vels[1].anomalizeT();
	}
	
	public void removeLagrangianTrend(){
		FilterMethods.removeLinearTrend(vels[0]);
		FilterMethods.removeLinearTrend(vels[1]);
	}
	
	
	/**
     * gridding the statistics variable to 2D variable by median position
     *
     * @param	maxlag	maximum lag time in the computation of lead-lag correlation
     * @param	dt		delta T of data interval (unit: s)
     */
	public Variable[] binningMeanByMedianPosition(Variable... vs){
		int len=ls.size();
		
		float[][] pos=new float[2][len];
		
		for(int i=0;i<len;i++){
			Particle p=ls.get(i);
			
			Record r=p.getRecord(p.getMedianIndex());
			
			pos[0][i]=r.getLon();
			pos[1][i]=r.getLat();
		}
		
		Variable[] re=null;//DataBaseUtil.binningMean(dd,pos[0],pos[1],false,vs);
		
		return re;
	}
	
	
	/**
     * calculate the autocorrelation function
     *
     * @param	v		a given variable
     * @param	maxlag	maximum lag time in the computation of lead-lag correlation
     */
	public static Variable cAutoCorrelationFunction(Variable vel,int maxLag){
		int t=vel.getTCount(),z=vel.getZCount();
		int y=vel.getYCount(),x=vel.getXCount();
		
		float undef=vel.getUndef();
		
		Range nr=new Range(maxLag+1,z,y,x);
		Range  r=vel.getRange();
		
		Variable R=new Variable("R_"+vel.getName(),vel.isTFirst(),nr);
		R.setCommentAndUnit("autocorrelation function of "+vel.getName());
		R.setUndef(undef);
		
		float[][][][] vdata=vel.getData();
		float[][][][] rdata=R.getData();
		
		if(vel.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cont:	for(int i=0;i<x;i++){
				float[] buf1=new float[t];
				
				for(int l=0;l<t;l++){
					buf1[l]=vdata[l][k][j][i];
					
					if(buf1[l]==undef){
						for(int lag=0;lag<=maxLag;lag++) rdata[lag][k][j][i]=undef;
						continue cont;
					}
				}
				
				float[] buf2=buf1.clone();
				
				rdata[0][k][j][i]=1;
				for(int lag=1;lag<=maxLag;lag++)
				rdata[lag][k][j][i]=StatisticsUtil.
				cLeadLagCorrelationCoefficient(buf1,buf2,lag);
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cont:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) if(vdata[k][j][i][l]==undef){
					for(int lag=0;lag<=maxLag;lag++) rdata[k][j][i][lag]=undef;
					continue cont;
				}
				
				rdata[k][j][i][0]=1;
				for(int lag=1;lag<=maxLag;lag++)
				rdata[k][j][i][lag]=StatisticsUtil.
				cLeadLagCorrelationCoefficient(vdata[k][j][i],vdata[k][j][i],lag);
			}
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0]+nr.getTRange()[2]-1;
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return R;
	}
	
	public static Variable cVelocityCovariance(Variable u1,Variable u2,int maxLag){
		int t=u1.getTCount(),z=u1.getZCount();
		int y=u1.getYCount(),x=u1.getXCount();
		
		float undef=u1.getUndef();
		
		Range nr=new Range(maxLag+1,z,y,x);
		Range  r=u1.getRange();
		
		Variable R=new Variable("P"+u1.getName()+u2.getName(),u1.isTFirst(),nr);
		R.setCommentAndUnit("velocity covariance of "+u1.getName()+" and "+u2.getName());
		R.setUndef(undef);
		
		float[][][][] u1data=u1.getData();
		float[][][][] u2data=u1.getData();
		float[][][][]  rdata=R.getData();
		
		int[] count=new int[maxLag];
		
		if(u1.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cont:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++)
				if(u1data[l][k][j][i]==undef||u2data[l][k][j][i]==undef){
					for(int lag=0;lag<=maxLag;lag++) rdata[lag][k][j][i]=undef;
					continue cont;
				}
				
				for(int lag=0;lag<=maxLag;lag++){
					rdata[lag][k][j][i]+=u1data[0][k][j][i]*u2data[lag][k][j][i];
					count[lag]++;
				}
			}
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int lag=0;lag<=maxLag;lag++) rdata[lag][k][j][i]/=count[lag];
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cont:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++)
				if(u1data[k][j][i][l]==undef||u2data[k][j][i][l]==undef){
					for(int lag=0;lag<=maxLag;lag++) rdata[k][j][i][lag]=undef;
					continue cont;
				}
				
				for(int lag=0;lag<=maxLag;lag++){
					rdata[k][j][i][lag]+=u1data[k][j][i][0]*u2data[k][j][i][lag];
					count[lag]++;
				}
			}
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int lag=0;lag<=maxLag;lag++) rdata[k][j][i][lag]/=count[lag];
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=r.getTRange()[0]+nr.getTRange()[2]-1;
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return R;
	}
	
	
	/**
     * calculate the characteristic time scale by integrating the autocorrelation function
     * using trapezoidal approximation along t-dimension
     *
     * @param	R		autocorrelation function
     * @param	lim		upper limit of integration along t-dimension (start from 0)
     */
	public static Variable cTimescaleByIntegrateTo1stZeroCross(Variable R,float dt){
		int t=R.getTCount(),z=R.getZCount();
		int y=R.getYCount(),x=R.getXCount();
		
		float undef=R.getUndef();
		
		Range nr=new Range(1,z,y,x);
		Range  r=R.getRange();
		
		Variable T=new Variable("ets_"+R.getName(),R.isTFirst(),nr);
		T.setCommentAndUnit("eddy timescale of "+R.getName()+" (s^-1)");
		T.setUndef(undef);
		
		float[][][][] vdata=R.getData();
		float[][][][] tdata=T.getData();
		
		if(R.isTFirst()){
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				float[] buf=new float[t];
				
				for(int l=0;l<t;l++) buf[l]=vdata[l][k][j][i];
				
				int lim=get1stZeroCrossing(buf);
				
				tdata[0][k][j][i]=trapezoidalIntegrate(buf,lim,dt,undef);
			}
			
		}else{
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				int lim=get1stZeroCrossing(vdata[k][j][i]);
				
				tdata[k][j][i][0]=trapezoidalIntegrate(vdata[k][j][i],lim,dt,undef);
			}
		}
		
		nr.setTRange(r.getTRange()[0]);
		nr.setZRange(r);
		nr.setYRange(r);
		nr.setXRange(r);
		
		return T;
	}
	
	public static Variable cTimescaleByFitting(Variable R,float dt){
		int t=R.getTCount(),z=R.getZCount();
		int y=R.getYCount(),x=R.getXCount();
		
		float undef=R.getUndef();
		
		Range nr=new Range(1,z,y,x);
		Range  r=R.getRange();
		
		Variable T=new Variable("ets_"+R.getName(),R.isTFirst(),nr);
		T.setCommentAndUnit("eddy timescale of "+R.getName()+" (s^-1)");
		T.setUndef(undef);
		
		float[][][][] tdata=T.getData();
		float[][][][] rdata=R.getData();
		
		if(R.isTFirst()){
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				float[] buf=new float[t];
				
				for(int l=0;l<t;l++) buf[l]=rdata[l][k][j][i];
				
				double TD=get1stZeroCrossing(buf)*dt;
				double tauE=cTauE(buf,TD,dt);
				double pow=Math.PI*tauE/4.0/TD;
				
				tdata[0][k][j][i]=(float)(Math.sqrt(Math.PI)/2.0*tauE*Math.exp(-pow*pow));
			}
			
		}else{
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				double TD =get1stZeroCrossing(rdata[k][j][i])*dt;
				double tauE=cTauE(rdata[k][j][i],TD,dt);
				double pow=Math.PI*tauE/4.0/TD;
				tdata[k][j][i][0]=(float)(Math.sqrt(Math.PI)/2.0*tauE*Math.exp(-pow*pow));
			}
		}
		
		nr.setTRange(r.getTRange()[0]);
		nr.setZRange(r);
		nr.setYRange(r);
		nr.setXRange(r);
		
		return T;
	}
	
	public static Variable cTimescalesByIntegration(Variable R,float dt){
		int t=R.getTCount(),z=R.getZCount();
		int y=R.getYCount(),x=R.getXCount();
		
		float undef=R.getUndef();
		
		Variable T=new Variable("TL",R);
		T.setCommentAndUnit("eddy characteristic timescale (s^-1)");
		
		float[][][][] vdata=R.getData();
		float[][][][] tdata=T.getData();
		
		if(R.isTFirst()){
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				float[] buf=new float[t];
				
				for(int l=0;l<t;l++) buf[l]=vdata[l][k][j][i];
				
				for(int l=0;l<t;l++)
				tdata[l][k][j][i]=trapezoidalIntegrate(buf,l,dt,undef);
			}
			
		}else{
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				for(int l=0;l<t;l++)
				tdata[k][j][i][l]=trapezoidalIntegrate(vdata[k][j][i],l,dt,undef);
			}
		}
		
		return T;
	}
	
	
	/**
     * change the list of data to Variable (y=1,z=1, x!=1,t!=1)
     */
	public static Variable[] toVelocityVariables(List<? extends Particle> particles){
		int x=particles.size();
		int t=particles.get(0).getTCount();
		
		for(Particle p:particles)
		if(p.getTCount()!=t)
		throw new IllegalArgumentException("invalid particle size: "+p.getTCount()+", should be "+t);
		
		Variable[] re=new Variable[2];
		
		re[0]=new Variable("uL",true,new Range(t,1,1,x));
		re[1]=new Variable("vL",true,new Range(t,1,1,x));
		
		re[0].setUndef(-9999);	re[0].setCommentAndUnit("lagrangian u data");
		re[1].setUndef(-9999);	re[1].setCommentAndUnit("lagrangian v data");
		
		float[][][][] udata=re[0].getData();
		float[][][][] vdata=re[1].getData();
		
		for(int i=0;i<x;i++){
			Particle p=particles.get(i);
			
			for(int l=0;l<t;l++){
				Record r=p.getRecord(l);
				
				udata[l][0][0][i]=r.getDataValue(0);
				vdata[l][0][0][i]=r.getDataValue(1);
			}
		}
		
		return re;
	}
	
	
	/*** help methods and class ***/
	private static float trapezoidalIntegrate(float[] data,int lim,float dt,float undef){
		if(lim==0) return 0;
		
		int len=data.length;
		
		for(int l=0;l<len;l++) if(data[l]==undef) return undef;
		
		float re=(data[0]+data[lim])/2;
		
		for(int l=1,L=lim;l<L;l++) re+=data[l];
		
		return re*dt;
	}
	
	private static int get1stZeroCrossing(float[] data){
		for(int l=1,L=data.length;l<L;l++) if(data[l]<0) return l-1;
		
		throw new IllegalArgumentException("no zero crossing");
	}
	
	private static double cTauE(float[] R,double TD,double dt){
		int len=R.length;
		
		double tauE=Double.NaN;
		
		PF.len=len;
		PF.TD =TD;
		PF.dt =dt;
		PF.R  =R;
		
		double str=TD*0.05;
		double end=TD*4;
		
		if(PF.value(str)*PF.value(end)>=0) throw new IllegalArgumentException("no solution");
		
		try{ tauE=BS.solve(300,PF,str,end);}
		catch(Exception e){ e.printStackTrace(); System.exit(0);}
		
		double[] RR=cR(len,TD,tauE,3600*6);
		for(double d:RR) System.out.println(d);
		
		System.out.println("TD, tauE and result:"+TD+"\t"+tauE+"\t"+PF.value(tauE));
		System.out.println("Eulerian dis:"+cErr(R,RR));
		
		return tauE;
	}
	
	private static double[] cR(int N,double TD,double tauE,double dt){
		double[] re=new double[N];
		
		for(int l=0;l<N;l++){
			double pow=l*dt/tauE;
			re[l]=Math.cos(Math.PI*l*dt/2.0/TD)*Math.exp(-pow*pow);
		}
		
		return re;
	}
	
	private static double cErr(float[] R1,double[] R2){
		int len=R1.length;
		double sum=0;
		
		for(int l=0;l<len;l++){
			double dis=R1[l]-R2[l];
			sum+=dis*dis;
		}
		
		return sum;
	}
	
	
	static final class PrescribedFunction implements UnivariateFunction{
		//
		 int   len=0;
		double TD =0;
		double dt =0;
		
		float[] R =null;
		
		//
		PrescribedFunction(){}
		
		public double value(double tauE){
			double re=0;
			
			for(int l=0;l<len;l++){
				double pow =l*dt/tauE;
				double tau2=pow*pow;
				double Rstr=Math.cos(Math.PI*l*dt/2.0/TD)*Math.exp(-tau2);
				re+=Rstr*tau2/tauE*(Rstr-R[l]);
			}
			
			return re;
		}
	}
	
	
	/** test
	public static void main(String arg[]){
		DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/ERAInterim/Data.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable u=df.getVariables(new Range("lon(150,150);lat(15,15);lev(850,850);t(1,1000)",dd),"u")[0];
		
		FilterMethods.removeLinearTrend(u);
		
		Variable R=TurbulenceMethodsInSC.cAutoCorrelationFunction(u,120);
		
		for(int l=0;l<R.getTCount();l++) System.out.println(R.getData()[0][0][0][l]);
		System.out.println("\n");
		Variable T1=TurbulenceMethodsInSC.cTimescaleByIntegrateTo1stZeroCross(R,3600*6);
		Variable T2=TurbulenceMethodsInSC.cTimescaleByFitting(R,3600*6);
		
		System.out.println(T1.getData()[0][0][0][0]);
		System.out.println(T2.getData()[0][0][0][0]);
	}*/
}
