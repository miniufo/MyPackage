/**
 * @(#)DynamicMethodsInCC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import java.util.function.Function;
import miniufo.application.EquationInCylindricalCoordinate;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.statistics.StatisticsUtil;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.diagnosis.SpatialModel.EARTH_ROTATE_SPEED;
import static miniufo.diagnosis.SpatialModel.GRAVITY_ACCERLERATION;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;


/**
 * basic analysis methods in cylindrical coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class DynamicMethodsInCC extends EquationInCylindricalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spatial model in cylindrical coordinates
     */
	public DynamicMethodsInCC(CylindricalSpatialModel csm){ super(csm);}
	
	
	/**
     * Calculate radial averaging (area-weighted)
     *
     * @param	v	a given Variable
     * @param	rs	the start of the ring count
     * @param	re	the end of the ring count (re >= rs >= 0)
     *
     * @return	radial mean
     */
	public Variable cRadialAverage(Variable v,int rs,int re){
		if(re<rs||rs<0) throw new IllegalArgumentException("invalid ring band (re >= rs >= 0)");
		
		assignSubDomainParams(v);
		
		Range nr=new Range(t,z,1,x);
		Range r =v.getRange();
		
		Variable nv=new Variable(v.getName(),v.isTFirst(),nr);
		nv.setUndef(undef);
		nv.setComment(v.getComment());
		nv.setUnit(v.getUnit());
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				float weight=0;	float sum=0;
				
				for(int j=rs;j<=re;j++)
				if(vdata[l][k][j][i]!=undef){ sum+=bsin[j]*vdata[l][k][j][i]; weight+=bsin[j];}
					
				if(weight!=0) nvdata[l][k][0][i]=sum/weight;
				else nvdata[l][k][0][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				float weight=0;	float sum=0;
				
				for(int j=rs;j<=re;j++)
				if(vdata[k][j][i][l]!=undef){ sum+=bsin[j]*vdata[k][j][i][l]; weight+=bsin[j];}
					
				if(weight!=0) nvdata[k][0][i][l]=sum/weight;
				else nvdata[k][0][i][l]=undef;
			}
		}
		
		nr.setTRange(r);
		nr.setZRange(r);
		nr.setXRange(r);
		nr.setYRange(r.getYRange()[0]);
		
		return nv;
	}
	
	
	/**
     * Calculate absolute angular momentum by Johnson and Downey (1975, Part II, MWR)
     *
     * @param	ut		tangential wind (m s^-1)
     *
     * @return	aam		absolute angular momentum (m^2 s^-1)
     */
	public Variable cAbsoluteAngularMomentumByJohnson(Variable ut){
		assignSubDomainParams(ut);
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable aam=new Variable("aam",ut);
		aam.setValue(undef);
		aam.setCommentAndUnit("absolute angular momentum (m^2 s^-1) by Johnson and Downey (1975, Part II, MWR)");
		
		float[]    fi0=csm.getOLat();
		float[]   beta=csm.getYDef().getSamples();
		float[][][] fi=csm.getLat();
		
		float[][][][] udata= ut.getData();
		float[][][][] adata=aam.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double Fi=
				(EARTH_ROTATE_SPEED*(sin(fi[l][j][i])+sin(fi0[l]))*(1.0-cos(beta[j]))*EARTH_RADIUS*EARTH_RADIUS);
				
				for(int k=0;k<z;k++)
				if(udata[l][k][j][i]!=undef) adata[l][k][j][i]=(float)(udata[l][k][j][i]*rs[j]+Fi);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double Fi=
				(EARTH_ROTATE_SPEED*(sin(fi[l][j][i])+sin(fi0[l]))*(1.0-cos(beta[j]))*EARTH_RADIUS*EARTH_RADIUS);
				
				for(int k=0;k<z;k++)
				if(udata[k][j][i][l]!=undef) adata[k][j][i][l]=(float)(udata[k][j][i][l]*rs[j]+Fi);
			}
		}
		
		return aam;
	}
	
	/**
     * Calculate absolute angular momentum by Johnson and Downey (1975, Part II, MWR)
     * with small-angle approximation (SAA).  See Yuan and Jian (2002, AAS)
     * This algorithm has been shown to be not consistent with the traditional one as well as without SAA.
     *
     * @param	ut		tangential wind (m s^-1)
     *
     * @return	aam		absolute angular momentum (m^2 s^-1)
     */
	public Variable cAbsoluteAngularMomentumByJohnsonWithSAA(Variable ut){
		assignSubDomainParams(ut);
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable aam=new Variable("aam",ut);
		aam.setValue(undef);
		aam.setCommentAndUnit("absolute angular momentum (m^2 s^-1) by Johnson and Downey (1975, Part II, MWR) with SAA");
		
		float[]    fi0=csm.getOLat();
		float[][][] fi=csm.getLat();
		
		float[][][][] udata= ut.getData();
		float[][][][] adata=aam.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double Fi=EARTH_ROTATE_SPEED*(sin(fi[l][j][i])+sin(fi0[l]))*rs[j]*rs[j]/2.0;
				
				for(int k=0;k<z;k++)
				if(udata[l][k][j][i]!=undef) adata[l][k][j][i]=(float)(udata[l][k][j][i]*rs[j]+Fi);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double Fi=EARTH_ROTATE_SPEED*(sin(fi[l][j][i])+sin(fi0[l]))*rs[j]*rs[j]/2.0;
				
				for(int k=0;k<z;k++)
				if(udata[k][j][i][l]!=undef) adata[k][j][i][l]=(float)(udata[k][j][i][l]*rs[j]+Fi);
			}
		}
		
		return aam;
	}
	
	/**
     * Calculate traditional absolute angular momentum (e.g., Holland 1983, QJRMS)
     *
     * @param	ut		tangential wind (m s^-1)
     *
     * @return	aam		absolute angular momentum (m^2 s^-1)
     */
	public Variable cAbsoluteAngularMomentum(Variable ut){
		assignSubDomainParams(ut);
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable aam=new Variable("aam",ut);
		aam.setCommentAndUnit("traditional absolute angular momentum (e.g., Holland 1983, QJRMS) (m^2 s^-1)");
		
		float[]         fi0=csm.getOLat();
		float[][][][] udata= ut.getData();
		float[][][][] adata=aam.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double Fi=EARTH_ROTATE_SPEED*sin(fi0[l])*rs[j]*rs[j];
				
				for(int k=0;k<z;k++)
				if(udata[l][k][j][i]!=undef) adata[l][k][j][i]=(float)(udata[l][k][j][i]*rs[j]+Fi);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double Fi=EARTH_ROTATE_SPEED*sin(fi0[l])*rs[j]*rs[j];
				
				for(int k=0;k<z;k++)
				if(udata[k][j][i][l]!=undef) adata[k][j][i][l]=(float)(udata[k][j][i][l]*rs[j]+Fi);
			}
		}
		
		return aam;
	}
	
	/**
     * Calculate absolute angular velocity = ut/r+f/2
     *
     * @param	ut		tangential wind (m s^-1)
     *
     * @return	absw	absolute angular velocity (m s^-1)
     */
	public Variable cAbsoluteAngularVelocity(Variable ut){
		assignSubDomainParams(ut);
		
		Variable absw=new Variable("absw",ut);
		absw.setValue(undef);
		absw.setCommentAndUnit("absolute angular velocity (s^-1)");
		
		float[][][]      fi=((CylindricalSpatialModel)sm).getLat();
		float[][][][] udata=  ut.getData();
		float[][][][] wdata=absw.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float f1=2f*EARTH_ROTATE_SPEED*(float)sin(fi[l][j][i]);
				
				for(int k=0;k<z;k++)
				if(udata[l][k][j][i]!=undef) wdata[l][k][j][i]=udata[l][k][j][i]/rs[j]+f1/2f;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float f1=2f*EARTH_ROTATE_SPEED*(float)sin(fi[l][j][i]);
				
				for(int k=0;k<z;k++)
				if(udata[k][j][i][l]!=undef) wdata[k][j][i][l]=udata[k][j][i][l]/rs[j]+f1/2f;
			}
		}
		
		return absw;
	}
	
	/**
     * Calculate vertical component of absolute voticity eta = f+(d(ut*r)/dr-d(vr)/da)/r
     *
     * @param	ut		tangential wind (m s^-1)
     * @param	vr		radial wind (m s^-1)
     *
     * @return	eta		vertical component of absolute vorticity (m s^-1)
     */
	public Variable cAbsoluteVorticity(Variable ut,Variable vr){
		checkDimensions(ut,vr);
		assignSubDomainParams(ut);
		
		Variable eta=new Variable("eta",ut);
		eta.setCommentAndUnit("vertical component of absolute voticity (s^-1)");
		
		float[][][]    lats=((CylindricalSpatialModel)sm).getLat();
		float[][][][] udata= ut.getData();
		float[][][][] vdata= vr.getData();
		float[][][][] edata=eta.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				float f1=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][0][0]);
				float um=0; int count=0;
				for(int i=0;i<x;i++) if(udata[l][k][1][i]!=undef){ um+=udata[l][k][1][i]; count++;}
				
				if(count!=0) for(int i=0;i<x;i++) edata[l][k][0][i]=f1+2*um/count/rs[1];
				else for(int i=0;i<x;i++) edata[l][k][0][i]=undef;
				
				for(int j=1;j<y-1;j++){
					/*** i==0 ***/
					if(udata[l][k][j][0]!=undef) edata[l][k][j][0]=
						2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][0])+
						(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2)+udata[l][k][j][0]/rs[j]-
						(vdata[l][k][j  ][1]-vdata[l][k][j][x-1])/(dxs[j]*2);
					else edata[l][k][j][0]=undef;
					
					for(int i=1;i<x-1;i++)
					if(udata[l][k][j][i]!=undef) edata[l][k][j][i]=
						2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i])+
						(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)+udata[l][k][j][i]/rs[j]-
						(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[j]*2);
					else edata[l][k][j][i]=undef;
					
					/*** i==x-1 ***/
					if(udata[l][k][j][x-1]!=undef) edata[l][k][j][x-1]=
						2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][x-1])+
						(udata[l][k][j+1][x-1]-udata[l][k][j-1][x-1])/(dy*2)+udata[l][k][j][x-1]/rs[j]-
						(vdata[l][k][j  ][0  ]-vdata[l][k][j  ][x-2])/(dxs[j]*2);
					else edata[l][k][j][x-1]=undef;
				}
				
				/*** j==y-1 ***/
				if(udata[l][k][y-1][0]!=undef) edata[l][k][y-1][0]=
					2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][0])+
					(udata[l][k][y-1][0]-udata[l][k][y-2][0  ])/dy+udata[l][k][y-1][0]/rs[y-1]-
					(vdata[l][k][y-1][1]-vdata[l][k][y-1][x-1])/(dxs[y-1]*2);
				else edata[l][k][y-1][0]=undef;
				
				for(int i=1;i<x-1;i++)
				if(udata[l][k][y-1][i]!=undef) edata[l][k][y-1][i]=
					2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][i])+
					(udata[l][k][y-1][i  ]-udata[l][k][y-2][i  ])/dy+udata[l][k][y-1][i]/rs[y-1]-
					(vdata[l][k][y-1][i+1]-vdata[l][k][y-1][i-1])/(dxs[y-1]*2);
				else edata[l][k][y-1][i]=undef;
				
				if(udata[l][k][y-1][x-1]!=undef) edata[l][k][y-1][x-1]=
					2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][x-1])+
					(udata[l][k][y-1][x-1]-udata[l][k][y-2][x-1])/dy+udata[l][k][y-1][x-1]/rs[y-1]-
					(vdata[l][k][y-1][0]-vdata[l][k][y-1][x-2])/(dxs[y-1]*2);
				else edata[l][k][y-1][x-1]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				float f1=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][0][0]);
				float um=0; int count=0;
				for(int i=0;i<x;i++) if(udata[k][1][i][l]!=undef){ um+=udata[k][1][i][l]; count++;}
				
				if(count!=0) for(int i=0;i<x;i++) edata[k][0][i][l]=f1+2*um/count/rs[1];
				else for(int i=0;i<x;i++) edata[k][0][i][l]=undef;
				
				for(int j=1;j<y-1;j++){
					/*** i==0 ***/
					if(udata[k][j][0][l]!=undef) edata[k][j][0][l]=
						2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][0])+
						(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2)+udata[k][j][0][l]/rs[j]-
						(vdata[k][j  ][1][l]-vdata[k][j][x-1][l])/(dxs[j]*2);
					else edata[k][j][0][l]=undef;
					
					for(int i=1;i<x-1;i++)
					if(udata[k][j][i][l]!=undef) edata[k][j][i][l]=
						2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i])+
						(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2)+udata[k][j][i][l]/rs[j]-
						(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[j]*2);
					else edata[k][j][i][l]=undef;
					
					/*** i==x-1 ***/
					if(udata[k][j][x-1][l]!=undef) edata[k][j][x-1][l]=
						2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][x-1])+
						(udata[k][j+1][x-1][l]-udata[k][j-1][x-1][l])/(dy*2)+udata[k][j][x-1][l]/rs[j]-
						(vdata[k][j  ][0  ][l]-vdata[k][j  ][x-2][l])/(dxs[j]*2);
					else edata[k][j][x-1][l]=undef;
				}
				
				/*** j==y-1 ***/
				if(udata[k][y-1][0][l]!=undef) edata[k][y-1][0][l]=
					2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][0])+
					(udata[k][y-1][0][l]-udata[k][y-2][0  ][l])/dy+udata[k][y-1][0][l]/rs[y-1]-
					(vdata[k][y-1][1][l]-vdata[k][y-1][x-1][l])/(dxs[y-1]*2);
				else edata[k][y-1][0][l]=undef;
				
				for(int i=1;i<x-1;i++)
				if(udata[k][y-1][i][l]!=undef) edata[k][y-1][i][l]=
					2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][i])+
					(udata[k][y-1][i  ][l]-udata[k][y-2][i  ][l])/dy+udata[k][y-1][i][l]/rs[y-1]-
					(vdata[k][y-1][i+1][l]-vdata[k][y-1][i-1][l])/(dxs[y-1]*2);
				else edata[k][y-1][i][l]=undef;
				
				if(udata[k][y-1][x-1][l]!=undef) edata[k][y-1][x-1][l]=
					2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][x-1])+
					(udata[k][y-1][x-1][l]-udata[k][y-2][x-1][l])/dy+udata[k][y-1][x-1][l]/rs[y-1]-
					(vdata[k][y-1][0][l]-vdata[k][y-1][x-2][l])/(dxs[y-1]*2);
				else edata[k][y-1][x-1][l]=undef;
			}
		}
		
		return eta;
	}
	
	/**
     * Calculate relative angular momentum
     *
     * @param	ut		tangential wind (m s^-1)
     *
     * @return	ram		relative angular momentum (m^2 s^-1)
     */
	public Variable cRelativeAngularMomentum(Variable ut){
		assignSubDomainParams(ut);
		
		Variable ram=new Variable("ram",ut);
		ram.setCommentAndUnit("relative angular momentum (m^2 s^-1)");
		
		float[][][][] udata= ut.getData();
		float[][][][] rdata=ram.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			if(udata[l][k][j][i]!=undef) rdata[l][k][j][i]=udata[l][k][j][i]*rs[j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int k=0;k<z;k++)
			if(udata[k][j][i][l]!=undef) rdata[k][j][i][l]=udata[k][j][i][l]*rs[j];
		}
		
		return ram;
	}
	
	/**
     * Calculate local inertial stability = (f+2ut/r)*eta
     *
     * @param	ut		tangential wind (m s^-1)
     * @param	vr		radial wind (m s^-1)
     *
     * @return	iner	local inertial stability (s^-2)
     */
	public Variable cInertialStability(Variable ut,Variable vr){
		Variable iner=cAbsoluteVorticity(ut,vr);
		
		iner.setName("iner");
		iner.setCommentAndUnit("local inertial stability (s^-2)");
		
		iner.multiplyEq(cAbsoluteAngularVelocity(ut).multiplyEq(2f));
		
		return iner;
	}
	
	/**
     * Calculate local inertial stability = (f+2ut/r)*eta normalized by f^2
     *
     * @param	ut		tangential wind (m s^-1)
     * @param	vr		radial wind (m s^-1)
     *
     * @return	eta		normalized local inertial stability
     */
	public Variable cInertialStabilityNorm(Variable ut,Variable vr){
		Variable iner=cInertialStability(ut,vr);
		
		iner.setName("inerN");
		iner.setCommentAndUnit("normalized local inertial stability (1)");
		
		float[][][]    lats=((CylindricalSpatialModel)sm).getLat();
		float[][][][] idata=iner.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double f2=2f*EARTH_ROTATE_SPEED*sin(lats[l][j][i]); f2*=f2;
				
				if(f2!=0) for(int k=0;k<z;k++) idata[l][k][j][i]/=f2;
				else for(int k=0;k<z;k++) idata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double f2=2f*EARTH_ROTATE_SPEED*sin(lats[l][j][i]); f2*=f2;
				
				if(f2!=0) for(int k=0;k<z;k++) idata[k][j][i][l]/=f2;
				else for(int k=0;k<z;k++) idata[k][j][i][l]=undef;
			}
		}
		
		return iner;
	}
	
	/**
	 * Calculate gradient wind from gradient wind balance
	 * relationship locally according to Yuan and Jian (2002, AAS).
	 * It is emphasized that the local gradient wind balance is not usually valid.
	 * 
	 * Reference: Yuan and Jian 2002, AAS
	 * 
	 * @param	fai		geopotential (m^2 s^-2)
	 * 
	 * @return	gw		gradient wind (m s^-1)
	 */
	public Variable cGradientWindByJohnson(Variable fai){
		assignSubDomainParams(fai);
		
		Variable gw=new Variable("gw",fai);
		gw.setValue(undef);
		gw.setCommentAndUnit("local gradient wind (m s^-1) computed following Yuan and Jian (2002, AAS)");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		float[]    olats=csm.getOLat();
		float[][][] lats=csm.getLat();
		
		float[][][][] gdata= gw.getData();
		float[][][][] fdata=fai.getData();
		
		if(fai.isTFirst()){
			for(int l=0;l<t;l++){
				float[] faim=new float[y];
				
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int k=0;k<z;k++)
					for(int i=0;i<x;i++)
					if(fdata[l][k][j][i]!=undef){ faim[j]+=fdata[l][k][j][i]; count++;}
					
					if(count!=0) faim[j]/=count;
					else faim[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=0;i<x;i++){
					double tmp=EARTH_ROTATE_SPEED*(sin(olats[l])+sin(lats[l][j][i]))*rs[j]/2.0;
					
					for(int k=0;k<z;k++){
						float re=(float)(sqrt(
							((fdata[l][k][j+1][i]-fdata[l][k][j-1][i])-(faim[j+1]-faim[j-1]))/(dy+dy)*rs[j]+tmp*tmp
						)-tmp);
						
						if(!Float.isNaN(re)) gdata[l][k][j][i]=re;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				float[] faim=new float[y];
				
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int k=0;k<z;k++)
					for(int i=0;i<x;i++)
					if(fdata[k][j][i][l]!=undef){ faim[j]+=fdata[k][j][i][l]; count++;}
					
					if(count!=0) faim[j]/=count;
					else faim[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				for(int i=0;i<x;i++){
					double tmp=EARTH_ROTATE_SPEED*(sin(olats[l])+sin(lats[l][j][i]))*rs[j]/2;
					
					for(int k=0;k<z;k++){
						float re=(float)(sqrt(
							((fdata[k][j+1][i][l]-fdata[k][j-1][i][l])-(faim[j+1]-faim[j-1]))/(dy+dy)*rs[j]+tmp*tmp
						)-tmp);
						
						if(!Float.isNaN(re)) gdata[k][j][i][l]=re;
					}
				}
			}
		}
		
		return gw;
    }
    
	/**
	 * Calculate traditional gradient wind from local gradient wind balance.
	 * It is emphasized that the local gradient wind balance is not always valid.
	 * 
	 * @param	fai		geopotential (m^2 s^-2)
	 * 
	 * @return	gw		gradient wind (m s^-1)
	 */
    public Variable cGradientWind(Variable fai){
    	assignSubDomainParams(fai);
    	
		Variable gw=new Variable("gw",fai);
		gw.setValue(undef);
		gw.setCommentAndUnit("local traditional gradient wind (m s^-1)");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		float[][][]    lats=csm.getLat();
		float[][][][] gdata= gw.getData();
		float[][][][] fdata=fai.getData();
		
		if(fai.isTFirst()){
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int j=1;j<y-1;j++){
				double tmp=EARTH_ROTATE_SPEED*sin(lats[l][j][i])*rs[j];
				
				for(int k=0;k<z;k++){
					float re=(float)(-tmp+sqrt(
						tmp*tmp+rs[j]*(fdata[l][k][j+1][i]-fdata[l][k][j-1][i])/(dy+dy)
					));
					
					if(!Float.isNaN(re)) gdata[l][k][j][i]=re;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int j=1;j<y-1;j++){
				double tmp=EARTH_ROTATE_SPEED*sin(lats[l][j][i])*rs[j];
				
				for(int k=0;k<z;k++){
					float re=(float)(-tmp+sqrt(
						tmp*tmp+rs[j]*(fdata[k][j+1][i][l]-fdata[k][j-1][i][l])/(dy+dy)
					));
					
					if(!Float.isNaN(re)) gdata[k][j][i][l]=re;
				}
			}
		}
		
		return gw;
    }
    
	/**
     * Calculate sum of centrifugal and Coriolis force = ut^2/r+f*ut
     *
     * @param	ut		tangential wind (m s^-1)
     *
     * @return	ccf		centrifugal and Coriolis force (m s^-2)
     */
	public Variable cCentrifugalCoriolisForce(Variable ut){
		assignSubDomainParams(ut);
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable ccf=new Variable("ccf",ut);
		ccf.setValue(undef);
		ccf.setCommentAndUnit("sum of centrifugal and Coriolis force (m s^-2)");
		
		float[][][]      fi=csm.getLat();
		float[][][][] udata= ut.getData();
		float[][][][] wdata=ccf.getData();
		
		if(ut.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double f1=2.0*EARTH_ROTATE_SPEED*sin(fi[l][j][i]);
				
				for(int k=0;k<z;k++) if(udata[l][k][j][i]!=undef&&rs[j]!=0)
				wdata[l][k][j][i]=(float)(udata[l][k][j][i]/rs[j]+f1)*udata[l][k][j][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				double f1=2.0*EARTH_ROTATE_SPEED*sin(fi[l][j][i]);
				
				for(int k=0;k<z;k++) if(udata[k][j][i][l]!=undef&&rs[j]!=0)
				wdata[k][j][i][l]=(float)(udata[k][j][i][l]/rs[j]+f1)*udata[k][j][i][l];
			}
		}
		
		return ccf;
	}
	
	/**
     * Calculate hydrostatic vertical velocity omega
     *
     * @param	w	vertical component of velocity (m s^-1)
     * @param	T	temperature (K)
     *
     * @return	o	static vertical velocity omega (Pa s^-1)
     */
    public Variable cHydrostaticOmega(Variable w,Variable T){
    	checkDimensions(w,T);
    	assignSubDomainParams(w);
    	
		Variable o=new Variable("omega",w);
		o.setUndef(undef);
		o.setCommentAndUnit("static vertical velocity omega (Pa s^-1)");
		
		float[][][][] wdata=w.getData();
		float[][][][] odata=o.getData();
		float[][][][] Tdata=T.getData();
		
		if(w.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[l][k][j][i]!=undef) odata[l][k][j][i]=
					-zdef[zstart-1+k]*GRAVITY_ACCERLERATION*wdata[l][k][j][i]/Rd/Tdata[l][k][j][i];
				else odata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[k][j][i][l]!=undef) odata[k][j][i][l]=
					-zdef[zstart-1+k]*GRAVITY_ACCERLERATION*wdata[k][j][i][l]/Rd/Tdata[k][j][i][l];
				else odata[k][j][i][l]=undef;
			}
		}
		
		return o;
     }
    
	/**
     * Calculate geopotential using hydrostatic balance relationship
     *
     * @param	T		temperature (K)
     * @param	h		geopotential (m^2 s^-2) for upper boundary values
     *
     * @return	geop	hydrostatic balanced geopotential (m^2 s^-2)
     */
    public Variable cHydrostaticGeopotential(Variable T,Variable h){
    	checkDimensions(h,T);
    	assignSubDomainParams(T);
    	
		Variable geop=h.copy();
		geop.setCommentAndUnit("hydrostatic balanced geopotential (m^2 s^-2)");
		
		float[][][][] gdata=geop.getData();
		float[][][][] Tdata=   T.getData();
		float[][][][] hdata=   h.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] deriv=new float[z];
				
				for(int k=0;k<z;k++) deriv[k]=-Rd*Tdata[l][k][j][i]/zdef[k];
				
				float[] re=RK2(hdata[l][z-1][j][i],deriv,-dz,false);
				
				for(int k=0;k<z;k++) gdata[l][k][j][i]=re[k];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] deriv=new float[z];
				
				for(int k=0;k<z;k++) deriv[k]=-Rd*Tdata[k][j][i][l]/zdef[k];
				
				float[] re=RK2(hdata[z-1][j][i][l],deriv,-dz,false);
				
				for(int k=0;k<z;k++) gdata[k][j][i][l]=re[k];
			}
		}
		
		return geop;
     }
    
	/**
     * Calculate temperature using hydrostatic balance relationship
     *
     * @param	h	geopotential (m^2 s^-2) for outer boundary values
     *
     * @return	T	temperature (K)
     */
    public Variable cHydrostaticTemperature(Variable h){
    	assignSubDomainParams(h);
    	
		Variable T=new Variable("T",h);
		T.setValue(undef);
		T.setCommentAndUnit("hydrostatic balanced temperature (K)");
		
		float[][][][] Tdata=T.getData();
		float[][][][] hdata=h.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int k=1,K=z-1;k<K;k++)
			Tdata[l][k][j][i]=-zdef[k]*(hdata[l][k+1][j][i]-hdata[l][k-1][j][i])/(dz+dz)/Rd;
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int k=1,K=z-1;k<K;k++)
			Tdata[k][j][i][l]=-zdef[k]*(hdata[k+1][j][i][l]-hdata[k-1][j][i][l])/(dz+dz)/Rd;
		}
		
		return T;
     }
    
	/**
     * Calculate the convergence of eddy angular momentum flux
     * following Molinari and Vollaro (1990), unit: m s^-1 day^-1
     *
     * @param	ua	storm-relative tangential wind anomaly (m s^-1)
     * @param	va	storm-relative radial wind anomaly (m s^-1)
     *
     * @return	EFC	convergence of eddy angular momentum flux (m s^-1 day^-1)
     */
	public Variable cREFC(Variable ua,Variable va){
		checkDimensions(ua,va);
		assignSubDomainParams(ua);
		
		Variable EFC=new Variable("refc",ua.isTFirst(),new Range(t,z,y,1));
		EFC.setUndef(undef);
		EFC.setCommentAndUnit("convergence of relative eddy angular momentum flux (m s^-1 day^-1)");
		
		float[][][][] efdata=EFC.getData();
		float[][][][] uadata= ua.getData();
		float[][][][] vadata= va.getData();
		
		if(EFC.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] buf=new float[y];
				
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(uadata[l][k][j][i]!=undef&&vadata[l][k][j][i]!=undef){
						buf[j]+=uadata[l][k][j][i]*vadata[l][k][j][i];	count++;
					}
					
					if(count!=0) buf[j]/=count;
					else buf[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef)
					efdata[l][k][j][0]=-60*60*24*(
						buf[j+1]*bsin[j+1]*bsin[j+1]-buf[j-1]*bsin[j-1]*bsin[j-1]
					)/dy/bsin[j]/bsin[j]/2;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] buf=new float[y];
				
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(uadata[k][j][i][l]!=undef&&vadata[k][j][i][l]!=undef){
						buf[j]+=uadata[k][j][i][l]*vadata[k][j][i][l];	count++;
					}
					
					if(count!=0) buf[j]/=count;
					else buf[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef)
					efdata[k][j][0][l]=-60*60*24*(
						buf[j+1]*bsin[j+1]*bsin[j+1]-buf[j-1]*bsin[j-1]*bsin[j-1]
					)/dy/bsin[j]/bsin[j]/2;
			}
		}
		
		EFC.getRange().setTRange(ua.getRange());
		EFC.getRange().setYRange(ua.getRange());
		EFC.getRange().setXRange(ua.getRange().getXRange()[0]);
		EFC.getRange().setZRange(ua.getRange());
		
		return EFC;
	}
	
	/**
     * Calculate localized the convergence of eddy angular momentum flux,
     * following Molinari and Vollaro (1990), unit: m s^-1 day^-1
     * except no tangential averaging
     *
     * @param	ua	storm-relative tangential wind anomaly (m s^-1)
     * @param	va	storm-relative radial wind anomaly (m s^-1)
     *
     * @return	EFC	localized convergence of eddy angular momentum flux (m s^-1 day^-1)
     */
	public Variable cLocalREFC(Variable ua,Variable va){
		checkDimensions(ua,va);
		assignSubDomainParams(ua);
		
		Variable EFC=new Variable("lcefc",ua.isTFirst(),new Range(t,z,y,x));
		EFC.setUndef(undef);
		EFC.setCommentAndUnit("localized convergence of relative eddy angular momentum flux (m s^-1 day^-1)");
		
		float[][][][] efdata=EFC.getData();
		float[][][][] uadata= ua.getData();
		float[][][][] vadata= va.getData();
		
		if(EFC.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				for(int j=1;j<y-1;j++)
				if(uadata[l][k][j+1][i]!=undef&&uadata[l][k][j-1][i]!=undef)
					efdata[l][k][j][i]=-60*60*24*(
						uadata[l][k][j+1][i]*vadata[l][k][j+1][i]*bsin[j+1]*bsin[j+1]-
						uadata[l][k][j-1][i]*vadata[l][k][j-1][i]*bsin[j-1]*bsin[j-1]
					)/dy/bsin[j]/bsin[j]/2;
				
				if(uadata[l][k][y-1][i]!=undef&&uadata[l][k][y-2][i]!=undef)
					efdata[l][k][y-1][i]=-60*60*24*(
						uadata[l][k][y-1][i]*vadata[l][k][y-1][i]*bsin[y-1]*bsin[y-1]-
						uadata[l][k][y-2][i]*vadata[l][k][y-2][i]*bsin[y-2]*bsin[y-2]
					)/dy/bsin[y-1]/bsin[y-1];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				for(int j=1;j<y-1;j++)
				if(uadata[k][j+1][i][l]!=undef&&uadata[k][j-1][i][l]!=undef)
					efdata[k][j][i][l]=-60*60*24*(
						uadata[k][j+1][i][l]*vadata[k][j+1][i][l]*bsin[j+1]*bsin[j+1]-
						uadata[k][j-1][i][l]*vadata[k][j-1][i][l]*bsin[j-1]*bsin[j-1]
					)/dy/bsin[j]/bsin[j]/2;
				
				if(uadata[k][y-1][i][l]!=undef&&uadata[k][y-2][i][l]!=undef)
					efdata[k][y-1][i][l]=-60*60*24*(
						uadata[k][y-1][i][l]*vadata[k][y-1][i][l]*bsin[y-1]*bsin[y-1]-
						uadata[k][y-2][i][l]*vadata[k][y-2][i][l]*bsin[y-2]*bsin[y-2]
					)/dy/bsin[y-1]/bsin[y-1];
			}
		}
		
		EFC.getRange().setTRange(ua.getRange());
		EFC.getRange().setYRange(ua.getRange());
		EFC.getRange().setXRange(ua.getRange().getXRange()[0]);
		EFC.getRange().setZRange(ua.getRange());
		
		return EFC;
	}
	
	/**
     * Calculate the planetary eddy angular momentum flux
     * following Molinari and Vollaro (1990), unit: m s^-1 day^-1
     *
     * @param	va		radial wind anomaly (m s^-1), (not storm-relative)
     *
     * @return	PEFC	planetary eddy angular momentum flux (m s^-1 day^-1)
     */
	public Variable cPEFC(Variable va){
		assignSubDomainParams(va);
		
		Variable PEFC=new Variable("pefc",va.isTFirst(),new Range(t,z,y,1));
		PEFC.setUndef(undef);
		PEFC.setCommentAndUnit("planetary eddy angular momentum flux (m s^-1 day^-1)");
		
		float[][][]     lats=((CylindricalSpatialModel)sm).getLat();
		float[][][][] efdata=PEFC.getData();
		float[][][][] vadata=  va.getData();
		
		if(PEFC.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					float[] f=new float[x];
					
					for(int i=0;i<x;i++)
					f[i]=2*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i]);
					
					int count=0;	float fm=StatisticsUtil.cArithmeticMean(f);
					for(int i=0;i<x;i++)
					if(vadata[l][k][j][i]!=undef){
						efdata[l][k][j][0]+=-(f[i]-fm)*vadata[l][k][j][i];
						count++;
					}
					
					if(count!=0) efdata[l][k][j][0]*=60f*60f*24f/count;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					float[] f=new float[x];
					
					for(int i=0;i<x;i++)
					f[i]=2*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i]);
					
					int count=0;	float fm=StatisticsUtil.cArithmeticMean(f);
					for(int i=0;i<x;i++)
					if(vadata[k][j][i][l]!=undef){
						efdata[k][j][0][l]+=-(f[i]-fm)*vadata[k][j][i][l];
						count++;
					}
					
					if(count!=0) efdata[k][j][0][l]*=60f*60f*24f/count;
				}
			}
		}
		
		PEFC.getRange().setTRange(va.getRange());
		PEFC.getRange().setYRange(va.getRange());
		PEFC.getRange().setXRange(va.getRange().getXRange()[0]);
		PEFC.getRange().setZRange(va.getRange());
		
		return PEFC;
	}
	
	/**
     * Calculate the convergence of absolute eddy angular momentum flux
     * following Yuan and Jian (2002, AAS), unit: m s^-1 day^-1
     *
     * @param	aa		absolute angular momentum anomaly (m^2 s^-1)
     * @param	va		storm-relative radial wind anomaly (m s^-1)
     *
     * @return	AEFC	absolute convergence of eddy angular momentum flux (m s^-1 day^-1)
     */
	public Variable cAEFC(Variable aa,Variable va){
		checkDimensions(aa,va);
		assignSubDomainParams(aa);
		
		Variable AEFC=new Variable("aefc",aa.isTFirst(),new Range(t,z,y,1));
		AEFC.setUndef(undef);
		AEFC.setCommentAndUnit("convergence of absolute eddy angular momentum flux (m s^-1 day^-1)");
		
		float[][][][] efdata=AEFC.getData();
		float[][][][] aadata= aa.getData();
		float[][][][] vadata= va.getData();
		
		if(AEFC.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] buf=new float[y];
				
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(aadata[l][k][j][i]!=undef&&vadata[l][k][j][i]!=undef){
						buf[j]+=aadata[l][k][j][i]*vadata[l][k][j][i];	count++;
					}
					
					if(count!=0) buf[j]/=count;
					else buf[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef)
					efdata[l][k][j][0]=-60*60*24*(buf[j+1]*bsin[j+1]-buf[j-1]*bsin[j-1])/(dy*2)/bsin[j]/rs[j];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] buf=new float[y];
				
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(aadata[k][j][i][l]!=undef&&vadata[k][j][i][l]!=undef){
						buf[j]+=aadata[k][j][i][l]*vadata[k][j][i][l];	count++;
					}
					
					if(count!=0) buf[j]/=count;
					else buf[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef)
					efdata[k][j][0][l]=-60*60*24*(buf[j+1]*bsin[j+1]-buf[j-1]*bsin[j-1])/(dy*2)/bsin[j]/rs[j];
			}
		}
		
		AEFC.getRange().setTRange(aa.getRange());
		AEFC.getRange().setYRange(aa.getRange());
		AEFC.getRange().setXRange(aa.getRange().getXRange()[0]);
		AEFC.getRange().setZRange(aa.getRange());
		
		return AEFC;
	}
	
	/**
     * Calculate the horizontal flux convergence of eddy heat, unit: K day^-1
     *
     * @param	ta		potential temperature anomaly (K)
     * @param	va		storm-relative radial wind anomaly (m s^-1)
     *
     * @return	htHFC	convergence of eddy heat flux (K day^-1)
     */
	public Variable cEddyHeatHFC(Variable ta,Variable va){
		checkDimensions(ta,va);
		assignSubDomainParams(ta);
		
		Variable htHFC=new Variable("htHFC",ta.isTFirst(),new Range(t,z,y,1));
		htHFC.setUndef(undef);
		htHFC.setCommentAndUnit("convergence of eddy heat flux (K day^-1)");
		
		float[][][][] efdata=htHFC.getData();
		float[][][][] tadata=   ta.getData();
		float[][][][] vadata=   va.getData();
		
		if(htHFC.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] buf=new float[y];
				
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(tadata[l][k][j][i]!=undef&&vadata[l][k][j][i]!=undef){
						buf[j]+=tadata[l][k][j][i]*vadata[l][k][j][i];	count++;
					}
					
					if(count!=0) buf[j]/=count;
					else buf[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef)
				efdata[l][k][j][0]=-60*60*24*(buf[j+1]*bsin[j+1]-buf[j-1]*bsin[j-1])/dy/bsin[j]/2;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] buf=new float[y];
				
				/*** Calculate average ***/
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(tadata[k][j][i][l]!=undef&&vadata[k][j][i][l]!=undef){
						buf[j]+=tadata[k][j][i][l]*vadata[k][j][i][l];	count++;
					}
					
					if(count!=0) buf[j]/=count;
					else buf[j]=undef;
				}
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef)
				efdata[k][j][0][l]=-60*60*24*(buf[j+1]*bsin[j+1]-buf[j-1]*bsin[j-1])/dy/bsin[j]/2;
			}
		}
		
		htHFC.getRange().setTRange(ta.getRange());
		htHFC.getRange().setYRange(ta.getRange());
		htHFC.getRange().setXRange(ta.getRange().getXRange()[0]);
		htHFC.getRange().setZRange(ta.getRange());
		
		return htHFC;
	}
	
	/**
     * Calculate radial-vertical (YZ) plane divergence.
     *
     * @param	vecY	azimuthal-mean radial component of a vector
     * @param	vecZ	azimuthal-mean vertical component of a vector
     *
     * @return	divYZ	divergence in YZ plane (vecY * m^-1 or vecZ * Pa^-1)
     */
	public Variable cYZDivergence(Variable vecY,Variable vecZ){
		assignSubDomainParams(vecY);
		checkDimensions(vecY,vecZ);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable divYZ=new Variable("divYZ",vecY);
		divYZ.setCommentAndUnit("divergence in YZ plane (vecY * m^-1 or vecZ * Pa^-1)");
		divYZ.setValue(undef);
		
		float[][][][] dvdata=divYZ.getData();
		float[][][][] vYdata= vecY.getData();
		float[][][][] vZdata= vecZ.getData();
		
		if(divYZ.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,Y=y-1;j<Y;j++)
			if(vYdata[l][k][j+1][0]!=undef&&vYdata[l][k][j-1][0]!=undef&&vZdata[l][k+1][j][0]!=undef&&vZdata[l][k-1][j][0]!=undef){
				dvdata[l][k][j][0]=
				(vYdata[l][k][j+1][0]*bsin[j+1]-vYdata[l][k][j-1][0]*bsin[j-1])/(2f*dy)/bsin[j]+
				(vZdata[l][k+1][j][0]-vZdata[l][k-1][j][0])/(2f*dz);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,Y=y-1;j<Y;j++)
			if(vYdata[k][j+1][0][l]!=undef&&vYdata[k][j-1][0][l]!=undef&&vZdata[k+1][j][0][l]!=undef&&vZdata[k-1][j][0][l]!=undef){
				dvdata[k][j][0][l]=
				(vYdata[k][j+1][0][l]*bsin[j+1]-vYdata[k][j-1][0][l]*bsin[j-1])/(2f*dy)/bsin[j]+
				(vZdata[k+1][j][0][l]-vZdata[k-1][j][0][l])/(2f*dz);
			}
		}
		
		return divYZ;
	}
	
	
	/**
     * Calculate vertical wind shear between two levels i.e., 200 - 850 hPa
     *
     * @param	u	zonal wind (m s^-1)
     * @param	v	meridional wind (m s^-1)
     *
     * @return	vws	[0] is u-shear, [1] is v-shear
     */
	public Variable[] cVerticalWindShear(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		if(z!=2) throw new IllegalArgumentException("z count should be 2 for difference");
		
		Variable[] vws=new Variable[2];
		vws[0]=new Variable("vwsu",u.isTFirst(),new Range(t,1,y,x)); vws[0].setUndef(undef);
		vws[1]=new Variable("vwsv",u.isTFirst(),new Range(t,1,y,x)); vws[1].setUndef(undef);
		vws[0].setCommentAndUnit("vertical wind shear in zonal direction (m s^-1)");
		vws[1].setCommentAndUnit("vertical wind shear in meridional direction (m s^-1)");
		
		float[][][][] s0data=vws[0].getData();
		float[][][][] s1data=vws[1].getData();
		float[][][][] udata =     u.getData();
		float[][][][] vdata =     v.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(udata[l][1][j][i]!=undef&&udata[l][0][j][i]!=undef)
					s0data[l][0][j][i]=udata[l][1][j][i]-udata[l][0][j][i];
				
				if(vdata[l][1][j][i]!=undef&&vdata[l][0][j][i]!=undef)
					s1data[l][0][j][i]=vdata[l][1][j][i]-vdata[l][0][j][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(udata[1][j][i][l]!=undef&&udata[0][j][i][l]!=undef)
					s0data[0][j][i][l]=udata[1][j][i][l]-udata[0][j][i][l];
				
				if(vdata[1][j][i][l]!=undef&&vdata[0][j][i][l]!=undef)
					s1data[0][j][i][l]=vdata[1][j][i][l]-vdata[0][j][i][l];
			}
		}
		
		vws[0].getRange().setTRange(u.getRange());
		vws[1].getRange().setTRange(u.getRange());
		vws[0].getRange().setYRange(u.getRange());
		vws[1].getRange().setYRange(u.getRange());
		vws[0].getRange().setXRange(u.getRange());
		vws[1].getRange().setXRange(u.getRange());
		vws[0].getRange().setZRange(u.getRange().getZRange()[0]);
		vws[1].getRange().setZRange(u.getRange().getZRange()[0]);
		
		return vws;
	}
	
	/**
     * Calculate the surface frictional stress rho*Cd*|V10|*V10 (z==1 only) using
     * bulk formula and a constant air density of 1.225 kg / m^3.
     *
     * @param	Cdfc	a function that maps the wind speed to Cd
     * @param	ut10	tangential wind at 10 meter (m s^-1)
     * @param	vr10	radial wind at 10 meter (m s^-1)
     *
     * @return	tau		frictional stress (kg m^-1 s^-2 or Pa), [0] is tangential and [1] is radial
     */
	public Variable[] cSurfaceFrictionalStress(Function<Float,Float> Cdfc,Variable ut10,Variable vr10){
		assignSubDomainParams(ut10);
		checkDimensions(ut10,vr10);
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		Variable[] tau=new Variable[2];
		tau[0]=new Variable("taux",ut10);
		tau[1]=new Variable("tauy",ut10);
		tau[0].setCommentAndUnit("surface frictional stress rho*Cd*|V10|*ut10 (kg m^-1 s^-2 or Pa)");
		tau[1].setCommentAndUnit("surface frictional stress rho*Cd*|V10|*vr10 (kg m^-1 s^-2 or Pa)");
		tau[0].setValue(undef);
		tau[1].setValue(undef);
		
		// At sea level and at 15°„C air has a density of approximately 1.225 kg/m^3
		// (0.001225 g/cm^3, 0.0023769 slug/ft^3, 0.0765 lbm/ft^3) according to ISA
		// (International Standard Atmosphere).  See wikipedia.
		final float rho=1.225f;
		float[][][][] utdata=ut10.getData();
		float[][][][] vrdata=vr10.getData();
		float[][][][] usdata=tau[0].getData();
		float[][][][] vsdata=tau[1].getData();
		
		if(ut10.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(utdata[l][0][j][i]!=undef&&vrdata[l][0][j][i]!=undef){
				float mag=(float)Math.hypot(utdata[l][0][j][i],vrdata[l][0][j][i]);
				
				usdata[l][0][j][i]=rho*Cdfc.apply(mag)*mag*utdata[l][0][j][i];
				vsdata[l][0][j][i]=rho*Cdfc.apply(mag)*mag*vrdata[l][0][j][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(utdata[0][j][i][l]!=undef&&vrdata[0][j][i][l]!=undef){
				float mag=(float)Math.hypot(utdata[0][j][i][l],vrdata[0][j][i][l]);
				
				usdata[0][j][i][l]=rho*Cdfc.apply(mag)*mag*utdata[0][j][i][l];
				vsdata[0][j][i][l]=rho*Cdfc.apply(mag)*mag*vrdata[0][j][i][l];
			}
		}
		
		return tau;
	}
	
	/**
     * Calculate the surface friction (z==1 only) defined as
     * -1/rho*dtau/dz or g*dtau/dp.
     *
     * @param	taux	tangential wind stress (kg m^-1 s^-2 or Pa)
     * @param	tauy	radial wind stress (kg m^-1 s^-2 or Pa)
     *
     * @return	fsfc	surface friction (m s^-2), [0] is tangential and [1] is radial
     */
	public Variable[] cSurfaceFriction(Variable taux,Variable tauy,Variable pblh){
		assignSubDomainParams(taux);
		checkDimensions(taux,tauy,pblh);
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		Variable[] fsfc=new Variable[2];
		fsfc[0]=new Variable("fsfcx",taux);
		fsfc[1]=new Variable("fsfcy",tauy);
		fsfc[0].setCommentAndUnit("surface tangential friction -1/rho*dtaux/dz (m s^-2)");
		fsfc[1].setCommentAndUnit("surface radial friction -1/rho*dtauy/dz (m s^-2)");
		fsfc[0].setValue(undef);
		fsfc[1].setValue(undef);
		
		// At sea level and at 15 ÔøΩÔøΩC air has a density of approximately 1.225 kg/m^3
		// (0.001225 g/cm^3, 0.0023769 slug/ft^3, 0.0765 lbm/ft^3) according to ISA
		// (International Standard Atmosphere).  See wikipedia.
		final float rho=1.225f;
		float[][][][] txdata=taux.getData();
		float[][][][] tydata=tauy.getData();
		float[][][][] pldata=pblh.getData();
		float[][][][] usdata=fsfc[0].getData();
		float[][][][] vsdata=fsfc[1].getData();
		
		if(taux.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(txdata[l][0][j][i]!=undef&&tydata[l][0][j][i]!=undef){
				usdata[l][0][j][i]=-1f/rho*txdata[l][0][j][i]/pldata[l][0][j][i];
				vsdata[l][0][j][i]=-1f/rho*tydata[l][0][j][i]/pldata[l][0][j][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(txdata[0][j][i][l]!=undef&&tydata[0][j][i][l]!=undef){
				usdata[0][j][i][l]=-1f/rho*txdata[0][j][i][l]/pldata[0][j][i][l];
				vsdata[0][j][i][l]=-1f/rho*tydata[0][j][i][l]/pldata[0][j][i][l];
			}
		}
		
		return fsfc;
	}
	
	/**
	 * Calculate balanced tangential wind, temperature, and geopotential fields using
	 * gradient-wind and hydrostatic balance relationships under spherical geometry.
     *
	 * @param	aamm	azimuthal-mean absolute angular momentum (m^2 s^-1)
     * @param	Tm		azimuthal-mean temperature (K)
     * @param	hm		azimuthal-mean geopotential (m^2 s^-2) for upper-outer boundary values
     *
     * @return	blc		balanced fields, [0] is temperature (K), [1] is geopotential (m^2 s^-2),
     * 					and [2] is gradient wind (m s^-1)
     */
	public Variable[] cGradientHydrostaticBalancedFieldsByT(Variable aamm,Variable Tm,Variable hm){
		assignSubDomainParams(aamm);
		checkDimensions(aamm,Tm,hm);
		
		Variable[] blc=new Variable[3];
		
		Variable hbc=cHydrostaticGeopotential(Tm,hm);	// only needs its outer boundary values
		blc[0]=cMeanBalancedTemperatureByCurvedTWB(aamm,Tm);
		blc[1]=cMeanGeopotentialByCurvedGWB(aamm,hbc);
		blc[2]=cMeanGradientWindByCurvedGWB(blc[1]);
		
		blc[0].setName("Tb"   ); blc[0].setCommentAndUnit("gradient- and hydrostatic-balanced temperature (K)");
		blc[1].setName("geopb"); blc[1].setCommentAndUnit("gradient- and hydrostatic-balanced geopotential (m^2 s^-2)");
		blc[2].setName("gdwb" ); blc[2].setCommentAndUnit("gradient- and hydrostatic-balanced gradient wind (m s^-1)");
		
		return blc;
	}
	
	/**
	 * Calculate balanced tangential wind, temperature, and geopotential fields using
	 * gradient-wind and hydrostatic balance relationships under spherical geometry.
     *
	 * @param	aamm	azimuthal-mean absolute angular momentum (m^2 s^-1)
     * @param	hm		azimuthal-mean geopotential (m^2 s^-2)
     *
     * @return	blc		balanced fields, [0] is temperature (K), [1] is geopotential (m^2 s^-2),
     * 					and [2] is gradient wind (m s^-1)
     */
	public Variable[] cGradientHydrostaticBalancedFieldsByH(Variable aamm,Variable hm){
		assignSubDomainParams(aamm);
		checkDimensions(aamm,hm);
		
		Variable[] blc=new Variable[3];
		
		Variable Tbc=cHydrostaticTemperature(hm);	// only needs its outer boundary values
		blc[0]=cMeanBalancedTemperatureByCurvedTWB(aamm,Tbc);
		blc[1]=cMeanGeopotentialByCurvedGWB(aamm,hm);
		blc[2]=cMeanGradientWindByCurvedGWB(blc[1]);
		
		blc[0].setName("Tb"   ); blc[0].setCommentAndUnit("gradient- and hydrostatic-balanced temperature (K)");
		blc[1].setName("geopb"); blc[1].setCommentAndUnit("gradient- and hydrostatic-balanced geopotential (m^2 s^-2)");
		blc[2].setName("gdwb" ); blc[2].setCommentAndUnit("gradient- and hydrostatic-balanced gradient wind (m s^-1)");
		
		return blc;
	}
	
	
	/**
     * Calculate azimuthal mean absolute angular momentum by Johnson and Downey (1975, Part II, MWR)
     *
     * @param	utm		azimuthal mean tangential wind (m s^-1)
     *
     * @return	am		azimuthal mean absolute angular momentum (m^2 s^-1)
     */
	public Variable cMeanAbsoluteAngularMomentumByJohnson(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable aamm=new Variable("aamm",utm);
		aamm.setValue(undef);
		aamm.setCommentAndUnit("azimuthal mean absolute angular momentum (m^2 s^-1) by Johnson and Downey (1975, MWR)");
		
		float[]   beta=csm.getYDef().getSamples();
		float[]    fi0=csm.getOLat();
		float[][][] fi=csm.getLat();
		
		float[][][][] udata= utm.getData();
		float[][][][] adata=aamm.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				double Fi=0;
				
				for(int i=0,I=sm.getXCount();i<I;i++)
				Fi+=EARTH_ROTATE_SPEED*(sin(fi[l][j][i])+sin(fi0[l]))*(1.0-cos(beta[j]))*EARTH_RADIUS*EARTH_RADIUS;
				Fi/=sm.getXCount();
				
				for(int k=0;k<z;k++)
				if(udata[l][k][j][0]!=undef) adata[l][k][j][0]=udata[l][k][j][0]*rs[j]+(float)Fi;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				double Fi=0;
				
				for(int i=0,I=sm.getXCount();i<I;i++)
				Fi+=EARTH_ROTATE_SPEED*(sin(fi[l][j][i])+sin(fi0[l]))*(1.0-cos(beta[j]))*EARTH_RADIUS*EARTH_RADIUS;
				Fi/=sm.getXCount();
				
				for(int k=0;k<z;k++)
				if(udata[k][j][0][l]!=undef) adata[k][j][0][l]=udata[k][j][0][l]*rs[j]+(float)Fi;
			}
		}
		
		return aamm;
	}
	
	/**
     * Calculate azimuthal-mean traditional absolute angular momentum
     *
     * @param	utm		azimuthal mean tangential wind (m s^-1)
     *
     * @return	aammm	azimuthal mean absolute angular momentum (m^2 s^-1)
     */
	public Variable cMeanAbsoluteAngularMomentum(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable aamm=new Variable("aamm",utm);
		aamm.setValue(undef);
		aamm.setCommentAndUnit("azimuthal mean traditional absolute angular momentum (m^2 s^-1)");
		
		float[]        olat=((CylindricalSpatialModel)sm).getOLat();
		float[][][][] udata= utm.getData();
		float[][][][] adata=aamm.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				double Fi=EARTH_ROTATE_SPEED*sin(olat[l])*rs[j]*rs[j];
				
				for(int k=0;k<z;k++) if(udata[l][k][j][0]!=undef)
				adata[l][k][j][0]=(float)(udata[l][k][j][0]*rs[j]+Fi);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				double Fi=EARTH_ROTATE_SPEED*sin(olat[l])*rs[j]*rs[j];
				
				for(int k=0;k<z;k++) if(udata[k][j][0][l]!=undef)
				adata[k][j][0][l]=(float)(udata[k][j][0][l]*rs[j]+Fi);
			}
		}
		
		return aamm;
	}
	
    /**
	 * Calculate azimuthal-mean absolute angular momentum using
	 * gradient wind balance relationship under spherical geometry.
	 * 
	 * @param	fm		azimuthal-mean geopotential (m^2 s^-2)
	 * 
	 * @return	aamm	azimuthal-mean absolute angular momentum (m^2 s^-1)
	 */
    public Variable cMeanAbsoluteAngularMomentumByCurvedGWB(Variable fm){
    	assignSubDomainParams(fm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable aamm=new Variable("aamm",fm);
		aamm.setValue(undef);
		aamm.setCommentAndUnit("azimuthal-mean absolute angular momentum (m^2 s^-1)");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		float[]       olats= csm.getOLat();
		float[][][][] adata=aamm.getData();
		float[][][][] fdata=  fm.getData();
		
		if(fm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				double tmp=EARTH_ROTATE_SPEED*sin(olats[l])*rs[j];
				
				for(int k=0;k<z;k++)
				adata[l][k][j][0]=(float)(sqrt(
					((fdata[l][k][j+1][0]-fdata[l][k][j-1][0])/(dy+dy)/bcos[j]*rs[j]+tmp*tmp)
				)*rs[j]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				double tmp=EARTH_ROTATE_SPEED*sin(olats[l])*rs[j];
				
				for(int k=0;k<z;k++)
				adata[k][j][0][l]=(float)(sqrt(
					((fdata[k][j+1][0][l]-fdata[k][j-1][0][l])/(dy+dy)/bcos[j]*rs[j]+tmp*tmp)
				)*rs[j]);
			}
		}
		
		return aamm;
    }
    
	/**
     * Calculate azimuthal mean absolute angular velocity = ut/r+f/2
     *
     * @param	utm		azimuthal mean tangential wind (m s^-1)
     *
     * @return	absw	absolute angular velocity (s^-1)
     */
	public Variable cMeanAbsoluteAngularVelocity(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable absw=new Variable("absw",utm);
		absw.setValue(undef);
		absw.setCommentAndUnit("absolute angular velocity (s^-1)");
		
		float[][][]      fi= csm.getLat();
		float[][][][] udata= utm.getData();
		float[][][][] wdata=absw.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				float f1=0;
				
				for(int i=0,I=sm.getXCount();i<I;i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(fi[l][j][i]);
				f1/=sm.getXCount();
				
				for(int k=0;k<z;k++)
				if(udata[l][k][j][0]!=undef&&rs[j]!=0) wdata[l][k][j][0]=udata[l][k][j][0]/rs[j]+f1/2f;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				float f1=0;
				
				for(int i=0,I=sm.getXCount();i<I;i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(fi[l][j][i]);
				f1/=sm.getXCount();
				
				for(int k=0;k<z;k++)
				if(udata[k][j][0][l]!=undef&&rs[j]!=0) wdata[k][j][0][l]=udata[k][j][0][l]/rs[j]+f1/2f;
			}
		}
		
		return absw;
	}
	
	/**
     * Calculate azimuthal mean vertical component of absolute vorticity = fm+d(utm*r)/dr/r
     *
     * @param	utm		azimuthal mean tangential wind (m s^-1)
     *
     * @return	etam	vertical component of absolute vorticity (s^-1)
     */
	public Variable cMeanAbsoluteVorticity(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable eta=new Variable("etam",utm);
		eta.setValue(undef);
		eta.setCommentAndUnit("azimuthal mean vertical component of absolute voticity (s^-1)");
		
		float[][][]    lats=((CylindricalSpatialModel)sm).getLat();
		float[][][][] udata=utm.getData();
		float[][][][] edata=eta.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				float f1=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][0][0]);
				if(udata[l][k][1][0]!=0) edata[l][k][0][0]=f1+2*udata[l][k][1][0]/rs[1];
				
				for(int j=1;j<y-1;j++){
					f1=0;
					for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i]);
					f1/=sm.getXCount();
					
					if(udata[l][k][j][0]!=undef) edata[l][k][j][0]=f1+
						(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2)+udata[l][k][j][0]/rs[j];
				}
				
				/*** j==y-1 ***/
				f1=0;
				for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][i]);
				f1/=sm.getXCount();
				
				if(udata[l][k][y-1][0]!=undef) edata[l][k][y-1][0]=f1+
					(udata[l][k][y-1][0]-udata[l][k][y-2][0])/dy+udata[l][k][y-1][0]/rs[y-1];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				float f1=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][0][0]);
				if(udata[k][1][0][l]!=0) edata[k][0][0][l]=f1+2*udata[k][1][0][l]/rs[1];
				
				for(int j=1;j<y-1;j++){
					f1=0;
					for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i]);
					f1/=sm.getXCount();
					
					if(udata[k][j][0][l]!=undef) edata[k][j][0][l]=f1+
						(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2)+udata[k][j][0][l]/rs[j];
				}
				
				/*** j==y-1 ***/
				f1=0;
				for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][i]);
				f1/=sm.getXCount();
				
				if(udata[k][y-1][0][l]!=undef) edata[k][y-1][0][l]=f1+
					(udata[k][y-1][0][l]-udata[k][y-2][0][l])/dy+udata[k][y-1][0][l]/rs[y-1];
			}
		}
		
		return eta;
	}
	
	/**
     * Calculate azimuthal mean vertical component of relative vorticity = d(utm*r)/dr/r
     *
     * @param	utm		azimuthal mean tangential wind (m s^-1)
     *
     * @return	zeta	vertical component of relative vorticity (s^-1)
     */
	public Variable cMeanRelativeVorticity(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable zeta=new Variable("zetam",utm);
		zeta.setValue(undef);
		zeta.setCommentAndUnit("azimuthal mean vertical component of relative voticity (s^-1)");
		
		float[][][][] udata= utm.getData();
		float[][][][] edata=zeta.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				if(udata[l][k][1][0]!=0) edata[l][k][0][0]=2*udata[l][k][1][0]/rs[1];
				
				for(int j=1;j<y-1;j++)
				if(udata[l][k][j][0]!=undef) edata[l][k][j][0]=
					(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2)+udata[l][k][j][0]/rs[j];
				
				/*** j==y-1 ***/
				if(udata[l][k][y-1][0]!=undef) edata[l][k][y-1][0]=
					(udata[l][k][y-1][0]-udata[l][k][y-2][0])/dy+udata[l][k][y-1][0]/rs[y-1];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				if(udata[k][1][0][l]!=0) edata[k][0][0][l]=2*udata[k][1][0][l]/rs[1];
				
				for(int j=1;j<y-1;j++)
				if(udata[k][j][0][l]!=undef) edata[k][j][0][l]=
					(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2)+udata[k][j][0][l]/rs[j];
				
				/*** j==y-1 ***/
				if(udata[k][y-1][0][l]!=undef) edata[k][y-1][0][l]=
					(udata[k][y-1][0][l]-udata[k][y-2][0][l])/dy+udata[k][y-1][0][l]/rs[y-1];
			}
		}
		
		return zeta;
	}
	
	/**
     * Calculate azimuthal mean vertical component of planetary vorticity
     *
     * @return	etam	vertical component of absolute vorticity (m s^-1)
     * 					only used to allocate the memory
     */
	public Variable cMeanPlanetaryVorticity(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable fm=new Variable("fm",utm);
		fm.setCommentAndUnit("azimuthal mean vertical component of planetary voticity (s^-1)");
		
		float[][][]    lats=((CylindricalSpatialModel)sm).getLat();
		float[][][][] edata=fm.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				edata[l][k][0][0]=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][0][0]);
				
				for(int j=1;j<y-1;j++){
					float f1=0;
					for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i]);
					edata[l][k][j][0]=f1/sm.getXCount();
				}
				
				/*** j==y-1 ***/
				float f1=0;
				for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][i]);
				edata[l][k][y-1][0]=f1/sm.getXCount();
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** j==0 ***/
				edata[k][0][0][l]=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][0][0]);
				
				for(int j=1;j<y-1;j++){
					float f1=0;
					for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][j][i]);
					edata[k][j][0][l]=f1/sm.getXCount();
				}
				
				/*** j==y-1 ***/
				float f1=0;
				for(int i=0;i<sm.getXCount();i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(lats[l][y-1][i]);
				edata[k][y-1][0][l]=f1/sm.getXCount();
			}
		}
		
		return fm;
	}
	
	/**
	 * Calculate azimuthal mean gradient wind from gradient wind balance
	 * relationship according to Johnson and Downey (1975, Part II, MWR).
	 * 
	 * Reference: Yuan and Jian 2002, AAS
	 * 
	 * @param	fm		azimuthal mean geopotential (m^2 s^-2)
	 * 
	 * @return	gwm		azimuthal mean gradient wind (m s^-1)
	 */
	public Variable cMeanGradientWindByJohnson(Variable fm){
		assignSubDomainParams(fm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable gwm=new Variable("gwm",fm);
		gwm.setValue(undef);
		gwm.setCommentAndUnit("azimuthal mean gradient wind (m s^-1) by Johnson and Downey (1975, Part II, MWR)");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		float[]       olats=csm.getOLat();
		float[][][][] gdata=gwm.getData();
		float[][][][] fdata= fm.getData();
		
		if(fm.isTFirst()){
			for(int l=0;l<t;l++){
				float[] faim=new float[y];
				
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int k=0;k<z;k++)
					if(fdata[l][k][j][0]!=undef){ faim[j]+=fdata[l][k][j][0]; count++;}
					
					if(count!=0) faim[j]/=count;
					else faim[j]=undef;
				}
				
				for(int j=1;j<y-1;j++){
					double tmp=2.0*EARTH_ROTATE_SPEED*sin(olats[l])*rs[j]/2.0;
					
					for(int k=0;k<z;k++){
						float re=(float)(sqrt(
							((fdata[l][k][j+1][0]-fdata[l][k][j-1][0])-(faim[j+1]-faim[j-1]))/(dy+dy)*rs[j]+tmp*tmp
						)-tmp);
						
						if(!Float.isNaN(re)) gdata[l][k][j][0]=re;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				float[] faim=new float[y];
				
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int k=0;k<z;k++)
					if(fdata[k][j][0][l]!=undef){
						faim[j]+=fdata[k][j][0][l];
						count++;
					}
					
					if(count!=0) faim[j]/=count;
					else faim[j]=undef;
				}
				
				for(int j=1;j<y-1;j++){
					double tmp=2.0*EARTH_ROTATE_SPEED*sin(olats[l])*rs[j]/2.0;
					
					for(int k=0;k<z;k++){
						float re=(float)(sqrt(
							((fdata[k][j+1][0][l]-fdata[k][j-1][0][l])-(faim[j+1]-faim[j-1]))/(dy+dy)*rs[j]+tmp*tmp
						)-tmp);
						
						if(!Float.isNaN(re)) gdata[k][j][0][l]=re;
					}
				}
			}
		}
		
		return gwm;
    }
	
	/**
	 * Calculate azimuthal-mean traditional gradient wind under spherical geometry
	 * 
	 * @param	fm		azimuthal-mean geopotential (m^2 s^-2)
	 * 
	 * @return	gw		azimuthal-mean gradient wind (m s^-1)
	 */
    public Variable cMeanGradientWindByCurvedGWB(Variable fm){
    	assignSubDomainParams(fm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable gwm=new Variable("gwm",fm);
		gwm.setValue(undef);
		gwm.setCommentAndUnit("azimuthal-mean traditional gradient wind (m s^-1)");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		float[]       olats=csm.getOLat();
		float[][][][] gdata=gwm.getData();
		float[][][][] fdata= fm.getData();
		
		if(fm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				double tmp=2.0*EARTH_ROTATE_SPEED*sin(olats[l])*rs[j]/2.0;
				
				for(int k=0;k<z;k++)
				gdata[l][k][j][0]=(float)(sqrt(
					(fdata[l][k][j+1][0]-fdata[l][k][j-1][0])/(dy+dy)/bcos[j]*rs[j]+tmp*tmp
				)-tmp);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				double tmp=2.0*EARTH_ROTATE_SPEED*sin(olats[l])*rs[j]/2.0;
				
				for(int k=0;k<z;k++)
				gdata[k][j][0][l]=(float)(sqrt(
					(fdata[k][j+1][0][l]-fdata[k][j-1][0][l])/(dy+dy)/bcos[j]*rs[j]+tmp*tmp
				)-tmp);
			}
		}
		
		return gwm;
    }
    
	/**
	 * Calculate azimuthal-mean tangential wind speed using
	 * traditional definition of absolute angular momentum.
	 * 
	 * @param	aamm	azimuthal-mean absolute angular momentum (m^2 s^-1)
	 * 
	 * @return	ut		azimuthal-mean tangential wind (m s^-1)
	 */
    public Variable cMeanTangentialWind(Variable aamm){
    	assignSubDomainParams(aamm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable utm=new Variable("utm",aamm);
		utm.setValue(undef);
		utm.setCommentAndUnit("azimuthal-mean tangential wind (m s^-1)");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		float[]       olats= csm.getOLat();
		float[][][][] udata= utm.getData();
		float[][][][] adata=aamm.getData();
		
		if(aamm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1;j<y;j++){
				double tmp=EARTH_ROTATE_SPEED*sin(olats[l])*rs[j];
				
				for(int k=0;k<z;k++) if(adata[l][k][j][0]!=undef)
				udata[l][k][j][0]=(float)(adata[l][k][j][0]/rs[j]-tmp);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1;j<y;j++){
				double tmp=EARTH_ROTATE_SPEED*sin(olats[l])*rs[j];
				
				for(int k=0;k<z;k++) if(adata[k][j][0][l]!=undef)
				udata[k][j][0][l]=(float)(adata[k][j][0][l]/rs[j]-tmp);
			}
		}
		
		return utm;
    }
	
	/**
     * Calculate the balanced temperature field using thermal-wind
     * balance relationship under spherical geometry.
     *
     * @param	aamm	azimuthal-mean absolute angular momentum (m^2 s^-1)
     * @param	Tm		azimuthal-mean temperature (K)
     *
     * @return	Tmb		balanced temperature (K)
     */
	public Variable cMeanBalancedTemperatureByCurvedTWB(Variable aamm,Variable Tm){
		assignSubDomainParams(aamm);
		checkDimensions(aamm,Tm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable Tmb=Tm.copy();
		Tmb.setName("Tb");
		Tmb.setCommentAndUnit("thermal-wind balanced temperature (K) under spherical geometry");
		
		float[][][][] Tbdata= Tmb.getData();
		float[][][][] Tmdata=  Tm.getData();
		float[][][][] amdata=aamm.getData();
		
		if(aamm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++){
				float[] deriv=new float[y];
				
				////// asignment //////
				for(int j=0;j<y;j++) deriv[j]=j==0?0:(float)
				(-zdef[k]*bcos[j]/Rd*2.0*amdata[l][k][j][0]/rs[j]/rs[j]/rs[j]*(amdata[l][k+1][j][0]-amdata[l][k-1][j][0])/(dz*2.0));
				
				float[] re=RK2(Tmdata[l][k][y-1][0],deriv,-dy,false);
				
				for(int j=0;j<y;j++) Tbdata[l][k][j][0]=re[j];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++){
				float[] deriv=new float[y];
				
				////// asignment //////
				for(int j=0;j<y;j++) deriv[j]=j==0?0:(float)
				(-zdef[k]*bcos[j]/Rd*2.0*amdata[k][j][0][l]/rs[j]/rs[j]/rs[j]*(amdata[k+1][j][0][l]-amdata[k-1][j][0][l])/(dz*2.0));
				
				float[] re=RK2(Tmdata[k][y-1][0][l],deriv,-dy,false);
				
				for(int j=0;j<y;j++) Tbdata[k][j][0][l]=re[j];
			}
		}
		
		return Tmb;
	}
	
	/**
	 * Calculate azimuthal-mean geopotential using gradient-wind balance
	 * relationship under spherical geometry.
	 * 
	 * @param	aamm	azimuthal-mean absolute angular momentum (m^2 s^-1)
	 * @param	hm		azimuthal-mean geopotential (m^2 s^-1) for outer boundary values
	 * 
	 * @return	geop	azimuthal-mean geopotential (m s^-1) satisfying gradient-wind balance
	 */
    public Variable cMeanGeopotentialByCurvedGWB(Variable aamm,Variable hm){
    	assignSubDomainParams(aamm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable geop=hm.copy();
		geop.setCommentAndUnit("azimuthal-mean geopotential (m^2 s^-2)");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		float[]       olats= csm.getOLat();
		float[][][][] gdata=geop.getData();
		float[][][][] adata=aamm.getData();
		
		if(aamm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] deriv=new float[y];
				
				for(int j=1;j<y;j++){
					double tmp=EARTH_ROTATE_SPEED*sin(olats[l]);
					
					deriv[j]=(float)(Math.pow(adata[l][k][j][0]/rs[j],2.0)/rs[j]-tmp*tmp*rs[j])*bcos[j];
				}
				
				float[] re=RK2(gdata[l][k][y-1][0],deriv,-dy,false);
				
				for(int j=0;j<y;j++) gdata[l][k][j][0]=re[j];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] deriv=new float[y];
				
				for(int j=1;j<y;j++){
					double tmp=EARTH_ROTATE_SPEED*sin(olats[l]);
					
					deriv[j]=(float)(Math.pow(adata[k][j][0][l]/rs[j],2.0)/rs[j]-tmp*tmp*rs[j])*bcos[j];
				}
				
				float[] re=RK2(gdata[k][y-1][0][l],deriv,-dy,false);
				
				for(int j=0;j<y;j++) gdata[k][j][0][l]=re[j];
			}
		}
		
		return geop;
    }
    
	/**
     * Calculate azimuthal mean of the sum of centrifugal and Coriolis force = ut^2/r+f*ut
     *
     * @param	utm		azimuthal mean tangential wind (m s^-1)
     *
     * @return	summ	azimuthal mean of the sum of centrifugal and Coriolis force (m s^-2)
     */
	public Variable cMeanCentrifugalCoriolisForce(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable ccfm=new Variable("ccfm",utm);
		ccfm.setValue(undef);
		ccfm.setCommentAndUnit("azimuthal mean of the sum of centrifugal and Coriolis force (m s^-2)");
		
		float[]       olats=csm.getOLat();
		float[][][][] udata=utm.getData();
		float[][][][] wdata=ccfm.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++){
				double f1=2f*EARTH_ROTATE_SPEED*sin(olats[l]);
				
				for(int j=0;j<y;j++)
				for(int k=0;k<z;k++) if(udata[l][k][j][0]!=undef&&rs[j]!=0)
				wdata[l][k][j][0]=(udata[l][k][j][0]/rs[j]+(float)f1)*udata[l][k][j][0];
			}
			
		}else{
			for(int l=0;l<t;l++){
				double f1=2f*EARTH_ROTATE_SPEED*sin(olats[l]);
				
				for(int j=0;j<y;j++)
				for(int k=0;k<z;k++) if(udata[k][j][0][l]!=undef&&rs[j]!=0)
				wdata[k][j][0][l]=(udata[k][j][0][l]/rs[j]+(float)f1)*udata[k][j][0][l];
			}
		}
		
		return ccfm;
	}
	
	/**
     * Calculate azimuthal mean inertial stability = (f+2ut/r)*(f+dut/dr+ut/r)
     * Reference: Molinari et al. (1990, JAS)
     *
     * @param	utm		azimuthal mean tangential wind (m s^-1)
     *
     * @return	ism		azimuthal mean inertial stability (s^-2)
     */
	public Variable cMeanInertialStabilityByUT(Variable utm){
		assignSubDomainParams(utm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		Variable ism=new Variable("ism",utm);
		ism.setValue(undef);
		ism.setCommentAndUnit("azimuthal mean of inertial stability (s^-2) by tangential wind");
		
		float[][][]      fi=csm.getLat();
		float[][][][] udata=utm.getData();
		float[][][][] idata=ism.getData();
		
		if(utm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				float f1=0;
				
				for(int i=0,I=sm.getXCount();i<I;i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(fi[l][j][i]);
				f1/=sm.getXCount();
				
				for(int k=0;k<z;k++) if(udata[l][k][j][0]!=undef) idata[l][k][j][0]=
				(f1+2f*udata[l][k][j][0]/rs[j])*(f1+(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(2*dy)+udata[l][k][j][0]/rs[j]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				float f1=0;
				
				for(int i=0,I=sm.getXCount();i<I;i++) f1+=2f*EARTH_ROTATE_SPEED*(float)sin(fi[l][j][i]);
				f1/=sm.getXCount();
				
				for(int k=0;k<z;k++) if(udata[k][j][0][l]!=undef) idata[k][j][0][l]=
				(f1+2f*udata[k][j][0][l]/rs[j])*(f1+(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(2*dy)+udata[k][j][0][l]/rs[j]);
			}
		}
		
		return ism;
	}
	
	/**
     * Calculate azimuthal mean inertial stability = dM^2/dr/r^3
     * Reference: Yuan and Jian (2002, AAS)
     *
     * @param	aamm	azimuthal mean absolute angular momentum (m^2 s^-1)
     *
     * @return	ism		azimuthal mean inertial stability (s^-2)
     */
	public Variable cMeanInertialStabilityByAM(Variable aamm){
		assignSubDomainParams(aamm);
    	
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable ism=new Variable("ism",aamm);
		ism.setValue(undef);
		ism.setCommentAndUnit("azimuthal mean of inertial stability (s^-2) by angular momentum");
		
		float[][][][] mdata=aamm.getData();
		float[][][][] idata=ism.getData();
		
		if(aamm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			if(mdata[l][k][j][0]!=undef) idata[l][k][j][0]=2*mdata[l][k][j][0]*
				(mdata[l][k][j+1][0]-mdata[l][k][j-1][0])/(dy*2)/rs[j]/rs[j]/rs[j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			if(mdata[k][j][0][l]!=undef) idata[k][j][0][l]=2*mdata[k][j][0][l]*
				(mdata[k][j+1][0][l]-mdata[k][j-1][0][l])/(dy*2)/rs[j]/rs[j]/rs[j];
		}
		
		return ism;
	}
	
	
	/**
     * Calculate the balanced temperature field.
     *
     * @param	Tm		azimuthal-mean temperature or potential temperature (K)
     */
	public void cDeviationFromOuterBoundary(Variable v){
		assignSubDomainParams(v);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) vdata[l][k][j][0]-=vdata[l][k][y-1][0];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) vdata[k][j][0][l]-=vdata[k][y-1][0][l];
		}
	}
	
	
	/** test
	public static void main(String arg[]){
		String name="Lupit";
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/"+name+"/"+name+".csm");
		CsmDescriptor csd=(CsmDescriptor)df.getDataDescriptor();
		CtlDescriptor ctl=csd.getCtlDescriptor();
		
		CylindricalSpacialModel csm=new CylindricalSpacialModel(csd);
		
		DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
		CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpacialModel(ctl),csm);
		
		Variable[] vs=df.getVariables(new Range("",csd),"u","v");
		Variable[] utvr=ct.reprojectToCylindrical(vs[0],vs[1]);
		dm.cRelativeVelocity(utvr[0],utvr[1]);
		
		Variable ut=utvr[0];
		Variable vr=utvr[1];
		Variable gaz=dm.cAbsoluteAngularMomentum(ut);
		
		ut.anomalizeX();
		vr.anomalizeX();
		
		Variable REFC=dm.cREFC(ut,vr);			REFC.setName("REFC");
		Variable EAMA=dm.cEAMAIndex(gaz,vr);	EAMA.setName("EAMA");
		Variable FFct=dm.cFFIndex(gaz,vr);		FFct.setName("FFCT");
		Variable FFbs=dm.cFFbsinIndex(gaz,vr);	FFbs.setName("FFbs");
		
		DataWrite dw=DataIOFactory.getDataWrite(ctl,"d:/FF.dat");
		dw.writeData(ctl,REFC,EAMA,FFct,FFbs);	dw.closeFile();
	}*/
}
