/**
 * @(#)EliassenModelInCC.java	1.0 2015.03.13
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.diagnosticModel;

import java.util.Arrays;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInCylindricalCoordinate;
import miniufo.basic.ArrayUtil;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static miniufo.geophysics.atmos.ThermoDynamics.Cp;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;
import static miniufo.geophysics.atmos.ThermoDynamics.kp;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.diagnosis.SpatialModel.EARTH_ROTATE_SPEED;


/**
 * Eliassen's balanced vortex model in cylindrical coordinate
 * Reference: Molinari and Vollaro (1990, JAS); Yuan and Jian (2002, AAS)
 *
 * @version 1.0, 2015.03.13
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EliassenModelInCC extends EquationInCylindricalCoordinate implements EllipticEquationInterface{
	//
	public static final float Cd=1.5e-3f;
	public static final float Kh=5.0e+3f;	// horizontal eddy viscous coefficient (m^2 s^-1)
	public static final float Kv=1.0e-3f;	// vertical eddy viscous coefficient in p-coordinate (Pa^2 s^-1)
	
	private short[][][] sfcIdx=null;		// vertical indices of the surface topography
	
	private float[] cs=null;		// central speed of the model (for compute the tilting)
	private float[] cd=null;		// central direction of the speed (for compute the tilting)
	
	private float[][] wtan=null;	// whole moving velocity projected to tangential direction (for compute the tilting)
	private float[][] wnor=null;	// whole moving velocity projected to radial direction (for compute the tilting)
	
	
	/**
     * constructor
     *
     * @param	csm		initialized by spatial model in Cylindrical coordinate
     */
	public EliassenModelInCC(CylindricalSpatialModel csm){
		super(csm);
		
		cWhole();
		
		sfcIdx=new short[csm.getTCount()][csm.getYCount()][csm.getXCount()];
	}
	
	/**
     * constructor
     *
     * @param	csm		initialized by spatial model in Cylindrical coordinate
     * @param	sfp		surface pressure
     */
	public EliassenModelInCC(CylindricalSpatialModel csm,Variable sfp){
		super(csm);
		
		cWhole();
		
		sfcIdx=new short[csm.getTCount()][csm.getYCount()][csm.getXCount()];
		
		markTopographyIdx(sfp);
	}
	
	
	/**
     * Assign surface friction (m s^-2) or heating rate (K s^-1) to pressure levels.
     *
     * @param	sfc		surface friction (m s^-2) or heating rate (K s^-1)
     * @param	sfh		surface geopotential height (m)
     * @param	pblh	planetary boundary layer height (m)
     * @param	hgt		geopotential height at each pressure levels (m)
     *
     * @return	sfcL	surface friction (m s^-2) or heating rate (K s^-1) at pressure levels
     */
	public Variable assignSurfaceToLevels(Variable sfc,Variable sfh,Variable pblh,Variable hgt){
		assignSubDomainParams(sfc);
		checkDimensions(sfc,sfh,pblh);
		
		if(!sfc.isAreaLike(hgt)) throw new IllegalArgumentException("dimensions not same");
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		z=hgt.getZCount();
		
		if(zdef[z-1]-zdef[0]>=0) throw new IllegalArgumentException("zdef is not decreasing upward (e.g., pressure)");
		
		Variable sfcL=new Variable(sfc.getName()+"L",hgt);
		sfcL.setCommentAndUnit("surface quantity at pressure levels");
		
		float[][][][] sfcdata = sfc.getData();
		float[][][][] hgtdata = hgt.getData();
		float[][][][] sfhdata = sfh.getData();
		float[][][][] sfLdata =sfcL.getData();
		float[][][][] pblhdata=pblh.getData();
		
		if(sfc.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] hlevs=new float[z];
				
				for(int k=0;k<z;k++) hlevs[k]=hgtdata[l][k][j][i];
				
				float pbl=pblhdata[l][0][j][i]==-9.99e8f?500:pblhdata[l][0][j][i];
				
				if(pbl!=pblhdata[l][0][j][i]) System.out.println(" warning: plbh data has undefined values");
				
				int bltIdx=ArrayUtil.getIdxIncre(hlevs,pbl+sfhdata[l][0][j][i]); // boundary layer top index
				
				if(bltIdx<0) bltIdx=0;
				
				if(bltIdx<sfcIdx[l][j][i]) throw new IllegalArgumentException("boundary layer height is below surface");
				
				for(int k2=0;k2<=bltIdx;k2++) sfLdata[l][k2][j][i]=sfcdata[l][0][j][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] hlevs=new float[z];
				
				for(int k=0;k<z;k++) hlevs[k]=hgtdata[k][j][i][l];
				
				float pbl=pblhdata[0][j][i][l]<-100?500:pblhdata[0][j][i][l];
				
				//if(pbl!=pblhdata[0][j][i][l]) System.out.println(" warning: plbh data has undefined values");
				
				int bltIdx=ArrayUtil.getIdxIncre(hlevs,pbl+sfhdata[0][j][i][l]); // boundary layer top index
				
				if(bltIdx<0) bltIdx=0;
				
				if(bltIdx<sfcIdx[l][j][i]){
					System.out.println(pbl+"\t"+pblhdata[0][j][i][l]+"\t"+sfhdata[0][j][i][l]);
					System.out.println(Arrays.toString(hlevs));
					System.out.println();
					throw new IllegalArgumentException("boundary layer height is below surface");
				}
				
				for(int k2=0;k2<=bltIdx;k2++) sfLdata[k2][j][i][l]=sfcdata[0][j][i][l];
			}
		}
		
		return sfcL;
	}
	
	/**
     * Calculate the horizontal viscous friction in free atmosphere
     * due to horizontal momentum dissipation: KH * (ï¿½ï¿½^2-1/r^2)ut.
     *
     * @param	utm		tangential-mean tangential wind (m s^-1)
     *
     * @return	frHVism	tangential-mean horizontal viscous friction (m s^-2)
     */
	public Variable cHorizontalViscousFriction(Variable utm){
		assignSubDomainParams(utm);
		
		Variable frHVism=new Variable("frHVism",utm);
		frHVism.setCommentAndUnit("horizontal viscous friction in free atmosphere (m s^-2)");
		frHVism.setValue(undef);
		
		float[][][][] udata=    utm.getData();
		float[][][][] fdata=frHVism.getData();
		
		if(frHVism.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(udata[l][k][j+1][0]!=undef&&udata[l][k][j-1][0]!=undef&&udata[l][k][j][0]!=undef)
			fdata[l][k][j][0]=
			(udata[l][k][j+1][0]+udata[l][k][j-1][0]-2f*udata[l][k][j][0])/dy/dy+
			(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(2f*dy)/rs[j]-
			udata[l][k][j][0]/rs[j]/rs[j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(udata[k][j+1][0][l]!=undef&&udata[k][j-1][0][l]!=undef&&udata[k][j][0][l]!=undef)
			fdata[k][j][0][l]=
			(udata[k][j+1][0][l]+udata[k][j-1][0][l]-2f*udata[k][j][0][l])/dy/dy+
			(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(2f*dy)/rs[j]-
			udata[k][j][0][l]/rs[j]/rs[j];
		}
		
		return frHVism;
	}
	
	/**
     * Calculate the vertical viscous friction in free atmosphere
     * due to vertical momentum dissipation: 1/g * KH * du/dp.
     *
     * @param	utm		tangential-mean tangential wind (m s^-1)
     *
     * @return	frVisVm	tangential friction (m s^-2)
     */
	public Variable cVerticalViscousFriction(Variable utm){
		assignSubDomainParams(utm);
		
		Variable frVisVm=new Variable("frVVism",utm);
		frVisVm.setCommentAndUnit("vertical viscous friction in free atmosphere (m s^-2)");
		frVisVm.setValue(undef);
		
		float[][][][] udata=  utm.getData();
		float[][][][] fdata=frVisVm.getData();
		
		if(frVisVm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int k=1,K=z-1;k<K;k++)
			if(udata[l][k+1][j][0]!=undef&&udata[l][k-1][j][0]!=undef&&udata[l][k][j][0]!=undef)
			fdata[l][k][j][0]=Kv*(udata[l][k+1][j][0]+udata[l][k-1][j][0]-2f*udata[l][k][j][0])/dz/dz;
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int k=1,K=z-1;k<K;k++)
			if(udata[k+1][j][0][l]!=undef&&udata[k-1][j][0][l]!=undef&&udata[k][j][0][l]!=undef)
			fdata[k][j][0][l]=Kv*(udata[k+1][j][0][l]+udata[k-1][j][0][l]-2f*udata[k][j][0][l])/dz/dz;
		}
		
		return frVisVm;
	}
	
	
	/**
     * Calculate force due to diabatic heating.
     *
     * @param	Qm	tangential-mean heating rate (W kg^-1 or m^2 s^-3)
     *
     * @return	f	force (m^2 kg^-1 s^-1)
     */
	public Variable cDiabaticHeatingForce(Variable Qm){
		assignSubDomainParams(Qm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("dhtFor",Qm);
		f.setCommentAndUnit("force due to diabatic heating (m^2 kg^-1 s^-1)");
		f.setValue(undef);
		
		float[][][][] fdata= f.getData();
		float[][][][] Qdata=Qm.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// core region (j==0), not used in inversion
				if(Qdata[l][k][0][0]!=undef&&Qdata[l][k][1][0]!=undef)
				fdata[l][k][0][0]=(Qdata[l][k][1][0]-Qdata[l][k][0][0])/dy*Rd/Cp/zdef[k];
				
				// inner region
				for(int j=1;j<y-1;j++)
				if(Qdata[l][k][j+1][0]!=undef&&Qdata[l][k][j-1][0]!=undef)
				fdata[l][k][j][0]=(Qdata[l][k][j+1][0]-Qdata[l][k][j-1][0])/(2f*dy)*Rd/Cp/zdef[k];
				
				// outer region (j==y-1), not used in inversion
				if(Qdata[l][k][y-1][0]!=undef&&Qdata[l][k][y-2][0]!=undef)
				fdata[l][k][y-1][0]=(Qdata[l][k][y-1][0]-Qdata[l][k][y-2][0])/dy*Rd/Cp/zdef[k];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// core region (j==0), not used in inversion
				if(Qdata[k][0][0][l]!=undef&&Qdata[k][1][0][l]!=undef)
				fdata[k][0][0][l]=(Qdata[k][1][0][l]-Qdata[k][0][0][l])/dy*Rd/Cp/zdef[k];
				
				// inner region
				for(int j=1;j<y-1;j++)
				if(Qdata[k][j+1][0][l]!=undef&&Qdata[k][j-1][0][l]!=undef)
				fdata[k][j][0][l]=(Qdata[k][j+1][0][l]-Qdata[k][j-1][0][l])/(2f*dy)*Rd/Cp/zdef[k];
				
				// outer region (j==y-1), not used in inversion
				if(Qdata[k][y-1][0][l]!=undef&&Qdata[k][y-2][0][l]!=undef)
				fdata[k][y-1][0][l]=(Qdata[k][y-1][0][l]-Qdata[k][y-2][0][l])/dy*Rd/Cp/zdef[k];
			}
		}
		
		return f;
	}
	
	/**
     * Calculate force due to eddy heat horizontal flux convergence (HFC)
     *
     * @param	tava	eddy heat horizontal flux = [Î¸'v'] (K m s^-1)
     *
     * @return	f		force (m^2 kg^-1 s^-1)
     */
	public Variable cEddyHeatHFCForce(Variable tava){
		assignSubDomainParams(tava);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("htHFCFor",tava);
		f.setCommentAndUnit("force due to eddy heat horizontal flux convergence (m^2 kg^-1 s^-1)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] tvdata=tava.getData();
		
		if(tava.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(tvdata[l][k][j+1][0]!=undef&&tvdata[l][k][j][0]!=undef&&tvdata[l][k][j-1][0]!=undef)
			fdata[l][k][j][0]=-(
				(tvdata[l][k][j+1][0]*bsin[j+1]-tvdata[l][k][j][0]*bsin[j])
				/((bsin[j+1]+bsin[j])/2f)-
				(tvdata[l][k][j][0]*bsin[j]-tvdata[l][k][j-1][0]*bsin[j-1])
				/((bsin[j]+bsin[j-1])/2f)
			)/(dy*dy)*Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(tvdata[k][j+1][0][l]!=undef&&tvdata[k][j][0][l]!=undef&&tvdata[k][j-1][0][l]!=undef)
			fdata[k][j][0][l]=-(
				(tvdata[k][j+1][0][l]*bsin[j+1]-tvdata[k][j][0][l]*bsin[j])
				/((bsin[j+1]+bsin[j])/2f)-
				(tvdata[k][j][0][l]*bsin[j]-tvdata[k][j-1][0][l]*bsin[j-1])
				/((bsin[j]+bsin[j-1])/2f)
			)/(dy*dy)*Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
		}
		
		return f;
	}
	
	/**
     * Calculate force due to eddy heat vertical flux convergence (VFC)
     *
     * @param	tawa	eddy heat vertical flux = [Î¸'w'] (K Pa s^-1)
     *
     * @return	f		force (m^2 kg^-1 s^-1)
     */
	public Variable cEddyHeatVFCForce(Variable tawa){
		assignSubDomainParams(tawa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("htVFCFor",tawa);
		f.setCommentAndUnit("force due to eddy heat vertical flux convergence (m^2 kg^-1 s^-1)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] twdata=tawa.getData();
		
		if(tawa.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(twdata[l][k+1][j+1][0]!=undef&&twdata[l][k-1][j+1][0]!=undef
			 &&twdata[l][k+1][j-1][0]!=undef&&twdata[l][k-1][j-1][0]!=undef)
			fdata[l][k][j][0]=-(
				(twdata[l][k+1][j+1][0]-twdata[l][k-1][j+1][0])/(2f*dz)-
				(twdata[l][k+1][j-1][0]-twdata[l][k-1][j-1][0])/(2f*dz)
			)/(2f*dy)*Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(twdata[k+1][j+1][0][l]!=undef&&twdata[k-1][j+1][0][l]!=undef
			 &&twdata[k+1][j-1][0][l]!=undef&&twdata[k-1][j-1][0][l]!=undef)
			fdata[k][j][0][l]=-(
				(twdata[k+1][j+1][0][l]-twdata[k-1][j+1][0][l])/(2f*dz)-
				(twdata[k+1][j-1][0][l]-twdata[k-1][j-1][0][l])/(2f*dz)
			)/(2f*dy)*Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
		}
		
		return f;
	}
	
	/**
     * Calculate force due to eddy absolute angular momentum (AAM)
     * horizontal flux convergence (HFC).
     *
     * @param	gm		tangential-mean absolute angular momentum (m^2 s^-1)
     * @param	gava	eddy AAM horizontal flux = [g'v'] (m^3 s^-2)
     *
     * @return	f		force (m^2 kg^-1 s^-1)
     */
	public Variable cEddyAAMHFCForce(Variable gm,Variable gava){
		assignSubDomainParams(gm);
		checkDimensions(gm,gava);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("aamHFCFor",gm);
		f.setCommentAndUnit("force due to eddy absolute angular momentum horizontal flux convergence (m^2 kg^-1 s^-1)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] gmdata=  gm.getData();
		float[][][][] gvdata=gava.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(gmdata[l][k+1][j+1][0]!=undef&&gmdata[l][k][j][0]!=undef&&gmdata[l][k-1][j][0]!=undef&&
			gvdata[l][k+1][j][0]!=undef&&gvdata[l][k][j][0]!=undef&&gvdata[l][k-1][j][0]!=undef)
			fdata[l][k][j][0]=-2f*(
				(gvdata[l][k+1][j+1][0]*bsin[j+1]-gvdata[l][k+1][j-1][0]*bsin[j-1])*
				gmdata[l][k+1][j][0]*bcos[j]/bsin[j]/(2f*dy)-
				(gvdata[l][k-1][j+1][0]*bsin[j+1]-gvdata[l][k-1][j-1][0]*bsin[j-1])*
				gmdata[l][k-1][j][0]*bcos[j]/bsin[j]/(2f*dy)
			)/(float)pow(rs[j],3.0)/(2f*dz);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(gmdata[k+1][j+1][0][l]!=undef&&gmdata[k][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&
			gvdata[k+1][j][0][l]!=undef&&gvdata[k][j][0][l]!=undef&&gvdata[k-1][j][0][l]!=undef)
			fdata[k][j][0][l]=-2f*(
				(gvdata[k+1][j+1][0][l]*bsin[j+1]-gvdata[k+1][j-1][0][l]*bsin[j-1])*
				gmdata[k+1][j][0][l]*bcos[j]/bsin[j]/(2f*dy)-
				(gvdata[k-1][j+1][0][l]*bsin[j+1]-gvdata[k-1][j-1][0][l]*bsin[j-1])*
				gmdata[k-1][j][0][l]*bcos[j]/bsin[j]/(2f*dy)
			)/(float)pow(rs[j],3.0)/(2f*dz);
		}
		
		return f;
	}
	
	/**
     * Calculate force due to eddy absolute angular momentum (AAM)
     * vertical flux convergence (VFC).
     *
     * @param	gm		tangential-mean absolute angular momentum (m^2 s^-1)
     * @param	gawa	eddy AAM vertical flux = [g'w'] (m^2 Pa s^-2)
     *
     * @return	f		force (m^2 kg^-1 s^-1)
     */
	public Variable cEddyAAMVFCForce(Variable gm,Variable gawa){
		assignSubDomainParams(gm);
		checkDimensions(gm,gawa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("aamVFCFor",gm);
		f.setCommentAndUnit("force due to eddy absolute angular momentum vertical flux convergence (m^2 kg^-1 s^-1)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] gmdata=  gm.getData();
		float[][][][] gwdata=gawa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1;j<y;j++)
			if(gmdata[l][k+1][j][0]!=undef&&gmdata[l][k][j][0]!=undef&&gmdata[l][k-1][j][0]!=undef&&
			gwdata[l][k+1][j][0]!=undef&&gwdata[l][k][j][0]!=undef&&gwdata[l][k-1][j][0]!=undef)
			fdata[l][k][j][0]=-2f*(
				(gmdata[l][k+1][j][0]+gmdata[l][k][j][0])/2f*(gwdata[l][k+1][j][0]-gwdata[l][k][j][0])-
				(gmdata[l][k-1][j][0]+gmdata[l][k][j][0])/2f*(gwdata[l][k][j][0]-gwdata[l][k-1][j][0])
			)*bcos[j]/(float)pow(rs[j],3)/(dz*dz);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1;j<y;j++)
			if(gmdata[k+1][j][0][l]!=undef&&gmdata[k][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&
			gwdata[k+1][j][0][l]!=undef&&gwdata[k][j][0][l]!=undef&&gwdata[k-1][j][0][l]!=undef)
			fdata[k][j][0][l]=-2f*(
				(gmdata[k+1][j][0][l]+gmdata[k][j][0][l])/2f*(gwdata[k+1][j][0][l]-gwdata[k][j][0][l])-
				(gmdata[k-1][j][0][l]+gmdata[k][j][0][l])/2f*(gwdata[k][j][0][l]-gwdata[k-1][j][0][l])
			)*bcos[j]/(float)pow(rs[j],3)/(dz*dz);
		}
		
		return f;
	}
	
	/**
     * Calculate frictional torque.
     *
     * @param	gm	tangential-mean absolute angular momentum (m^2 s^-1)
     * @param	fm	tangential-mean tangential friction (m s^-2)
     *
     * @return	f	force (m^2 kg^-1 s^-1)
     */
	public Variable cFrictionalTorqueForce(Variable gm,Variable fm){
		assignSubDomainParams(gm);
		
		if(!fm.isAreaLike(gm)) throw new IllegalArgumentException("areal not same for:\n"+fm+"\n"+gm);
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("frtqFor",gm);
		f.setCommentAndUnit("force due to frictional torque (m^2 kg^-1 s^-1)");
		
		float[][][][] ftdata= f.getData();
		float[][][][] frdata=fm.getData();
		float[][][][] gmdata=gm.getData();
		
		if(gm.isTFirst()){
			for(int j=1;j<y-1;j++){
				float tmp=(float)Math.pow(rs[j],2);
				
				for(int l=0;l<t;l++)
				for(int k=1;k<z-1;k++)
				if(frdata[l][k-1][j][0]!=undef&&frdata[l][k+1][j][0]!=undef&&gmdata[l][k][j][0]!=undef)
				ftdata[l][k][j][0]=
				2f*(frdata[l][k+1][j][0]*gmdata[l][k+1][j][0]-frdata[l][k-1][j][0]*gmdata[l][k-1][j][0])*bcos[j]/(2f*dz)/tmp;
			}
			
		}else{
			for(int j=1;j<y-1;j++){
				float tmp=(float)Math.pow(rs[j],2);
				
				for(int l=0;l<t;l++)
				for(int k=1;k<z-1;k++)
				if(frdata[k-1][j][0][l]!=undef&&frdata[k+1][j][0][l]!=undef&&gmdata[k][j][0][l]!=undef)
				ftdata[k][j][0][l]=
				2f*(frdata[k+1][j][0][l]*gmdata[k+1][j][0][l]-frdata[k-1][j][0][l]*gmdata[k-1][j][0][l])*bcos[j]/(2f*dz)/tmp;
			}
		}
		
		return f;
	}
	
	/**
     * Calculate tilting force.
     *
     * @param	ut	tangential wind component (not storm-relative) (m s^-1)
     * @param	vr	radial wind component (not storm-relative) (m s^-1)
     * @param	gm	tangential-mean absolute angular momentum (m^2 s^-1)
     *
     * @return	f	force (m^2 kg^-1 s^-1)
     */
	public Variable cTiltingForce(Variable gm,Variable ut,Variable vr){
		assignSubDomainParams(ut);
		checkDimensions(ut,vr);
		
		if(gm.getXCount()!=1)
			throw new IllegalArgumentException("x-direction of 1st param should be averaged");
		
		Variable f=new Variable("tiltFor",gm.isTFirst(),new Range(t,z,y,1));
		f.setCommentAndUnit("force due to tilting effect (m^2 kg^-1 s^-1)");
		f.setUndef(undef);
		f.setValue(undef);
		
		float til2,til3,til4,til5;
    	float[] olat=((CylindricalSpatialModel)sm).getOLat();
    	float[][] buf=new float[z][y];
		float[][][][]  fdata = f.getData();
		float[][][][] gmdata=gm.getData();
		float[][][][] utdata=ut.getData();
		float[][][][] vrdata=vr.getData();

		if(f.isTFirst()){
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(utdata[l][k][j][i]!=undef&&vrdata[l][k][j][i]!=undef){
						til3=(1f-(float)cos(ydef[j]))
							*(wtan[l][i]*vrdata[l][k][j][i]+wnor[l][i]*utdata[l][k][j][i]);
						
						til4=-vrdata[l][k][j][i]*bsin[j]*(float)(
							cos(ydef[j])*sin(olat[l])+bsin[j]*cos(olat[l])*cos(xdef[i])-sin(olat[l])
						);
						
						til5=(1f-(float)cos(ydef[j]))*(float)(
							(utdata[l][k][j][i]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[l][k][j][i]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i])
						);
						
						buf[k][j]+=til3+EARTH_RADIUS*EARTH_ROTATE_SPEED*(til4+til5);  count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					
					til2=-EARTH_ROTATE_SPEED*cs[l]*rs[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2f;
					buf[k][j]+=til2;
				}

				/*** partial p ***/
				for(int j=1;j<y;j++)
				for(int k=1;k<z-1;k++)
				if(gmdata[l][k+1][j][0]!=undef&&gmdata[l][k-1][j][0]!=undef&&buf[k+1][j]!=undef&&buf[k-1][j]!=undef)
					fdata[l][k][j][0]=
					2f*(gmdata[l][k+1][j][0]*buf[k+1][j]-gmdata[l][k-1][j][0]*buf[k-1][j])
					*bcos[j]/(float)pow(rs[j],3.0)/(dz*2f);
				
				/*** the first level ***/
				for(int j=1;j<y;j++)
				if(gmdata[l][1][j][0]!=undef&&gmdata[l][0][j][0]!=undef&&buf[1][j]!=undef&&buf[0][j]!=undef)
					fdata[l][0][j][0]=
					2f*(gmdata[l][1][j][0]*buf[1][j]-gmdata[l][0][j][0]*buf[0][j])
					*bcos[j]/(float)pow(rs[j],3.0)/dz;
				
				/*** the last level ***/
				for(int j=1;j<y;j++)
				if(gmdata[l][z-1][j][0]!=undef&&gmdata[l][z-2][j][0]!=undef&&buf[z-1][j]!=undef&&buf[z-2][j]!=undef)
					fdata[l][z-1][j][0]=
					2f*(gmdata[l][z-1][j][0]*buf[z-1][j]-gmdata[l][z-2][j][0]*buf[z-2][j])
					*bcos[j]/(float)pow(rs[j],3.0)/dz;
				
				/*** reset buf ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=0;
			}
			
		}else{
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(utdata[k][j][i][l]!=undef&&vrdata[k][j][i][l]!=undef){
						til3=(1f-(float)cos(ydef[j]))
							*(wtan[l][i]*vrdata[k][j][i][l]+wnor[l][i]*utdata[k][j][i][l]);
						
						til4=-vrdata[k][j][i][l]*bsin[j]
						    *(float)(cos(ydef[j])*sin(olat[l])
						    +bsin[j]*cos(olat[l])*cos(xdef[i])
						    -sin(olat[l]));
						
						til5=(1f-(float)cos(ydef[j]))*(float)(
							(utdata[k][j][i][l]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[k][j][i][l]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i]));
						
						buf[k][j]+=til3+EARTH_RADIUS*EARTH_ROTATE_SPEED*(til4+til5);  count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					
					til2=-EARTH_ROTATE_SPEED*cs[l]*rs[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2f;
					buf[k][j]+=til2;
				}

				/*** partial p ***/
				for(int j=1;j<y;j++)
				for(int k=1;k<z-1;k++)
				if(gmdata[k+1][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&buf[k+1][j]!=undef&&buf[k-1][j]!=undef)
					fdata[k][j][0][l]=
					2f*(gmdata[k+1][j][0][l]*buf[k+1][j]-gmdata[k-1][j][0][l]*buf[k-1][j])
					*bcos[j]/(float)pow(rs[j],3.0)/(dz*2f);
				
				/*** the first level ***/
				for(int j=1;j<y;j++)
				if(gmdata[1][j][0][l]!=undef&&gmdata[0][j][0][l]!=undef&&buf[1][j]!=undef&&buf[0][j]!=undef)
					fdata[0][j][0][l]=
					2f*(gmdata[1][j][0][l]*buf[1][j]-gmdata[0][j][0][l]*buf[0][j])
					*bcos[j]/(float)pow(rs[j],3.0)/dz;
				
				/*** the last level ***/
				for(int j=1;j<y;j++)
				if(gmdata[z-1][j][0][l]!=undef&&gmdata[z-2][j][0][l]!=undef&&buf[z-1][j]!=undef&&buf[z-2][j]!=undef)
					fdata[z-1][j][0][l]=
					2f*(gmdata[z-1][j][0][l]*buf[z-1][j]-gmdata[z-2][j][0][l]*buf[z-2][j])
					*bcos[j]/(float)pow(rs[j],3.0)/dz;
				
				/*** reset buf ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=0;
			}
		}
		
		return f;
	}
	
	
	/**
     * Calculate forcing factor due to heat functions,
     * i.e., FF = heat * R * ¦Ð / p.
     *
     * @param	heat	heat function e.g., eddy heat HFC or VFC (K s^-1)
     *
     * @return	ff		heat forcing factor (m^3 kg^-1 s^-1)
     */
	public Variable cHeatFF(Variable heat){
		assignSubDomainParams(heat);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable ff=new Variable(heat.getName()+"FF",heat);
		ff.setCommentAndUnit("forcing factor due to "+heat.getName()+" (m^3 kg^-1 s^-1)");
		ff.setValue(undef);
		
		float[][][][] fdata=  ff.getData();
		float[][][][] hdata=heat.getData();
		
		if(ff.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(hdata[l][k][j][0]!=undef) fdata[l][k][j][0]=
			hdata[l][k][j][0]*Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(hdata[k][j][0][l]!=undef) fdata[k][j][0][l]=
			hdata[k][j][0][l]*Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
		}
		
		return ff;
	}
	
	/**
     * Calculate diabatic heating rate d¦È/dt = Qm / (Cp * ¦Ð).
     *
     * @param	Qm	tangential-mean heating rate (W kg^-1 or m^2 s^-3)
     *
     * @return	f	diabatic heating rate (K s^-1)
     */
    public Variable cDiabaticHeatingRate(Variable Qm){
		assignSubDomainParams(Qm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("dhr",Qm);
		f.setCommentAndUnit("diabatic heating rate (K s^-1)");
		f.setValue(undef);
		
		float[][][][] fdata= f.getData();
		float[][][][] Qdata=Qm.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) if(Qdata[l][k][j][0]!=undef)
			fdata[l][k][j][0]=Qdata[l][k][j][0]/Cp/(float)pow(zdef[k]/100000.0,kp);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) if(Qdata[k][j][0][l]!=undef)
			fdata[k][j][0][l]=Qdata[k][j][0][l]/Cp/(float)pow(zdef[k]/100000.0,kp);
		}
		
		return f;
	}
    
	/**
     * Calculate eddy heat horizontal flux convergence (HFC).
     *
     * @param	tavam	eddy heat horizontal flux = [¦È'v'] (K m s^-1)
     *
     * @return	f		eddy heat HFC (K s^-1)
     */
	public Variable cEddyHeatHFC(Variable tavam){
		if(tavam.getXCount()!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable htHFC=cRadialDiv(tavam).multiplyEq(-1f);
		htHFC.setName("htHFC");
		htHFC.setCommentAndUnit("eddy heat horizontal flux convergence (K s^-1)");
		
		return htHFC;
	}
	
	/**
     * Calculate eddy heat vertical flux convergence (VFC).
     *
     * @param	tawam	eddy heat vertical flux = [¦È'¦Ø'] (K Pa s^-1)
     *
     * @return	f		eddy heat VHC (K s^-1)
     */
	public Variable cEddyHeatVFC(Variable tawam){
		if(tawam.getXCount()!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable htVFC=cVerticalDiv(tawam).multiplyEq(-1f);
		htVFC.setName("htVFC");
		htVFC.setCommentAndUnit("eddy heat vertical flux convergence (K s^-1)");
		
		return htVFC;
	}
	
	
	/**
     * Calculate forcing factor due to momentum functions,
     * i.e., FF = mome * 2 * gm * sin(lat) / rcos ^ 3.
     *
     * @param	gm		tangential-mean absolute angular momentum (m^2 s^-1)
     * @param	mome	momentum function e.g., eddy AAM HFC or VFC (m^2 s^-2)
     *
     * @return	ff		momentum forcing factor (m s^-3)
     */
	public Variable cMomentumFF(Variable gm,Variable mome){
		checkDimensions(gm,mome);
		assignSubDomainParams(mome);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable ff=new Variable(mome.getName()+"FF",mome);
		ff.setCommentAndUnit("forcing factor due to "+mome.getName()+" (m s^-3)");
		ff.setValue(undef);
		
		float[][][][] fdata=  ff.getData();
		float[][][][] gdata=  gm.getData();
		float[][][][] hdata=mome.getData();
		
		if(ff.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(hdata[l][k][j][0]!=undef) fdata[l][k][j][0]=
			hdata[l][k][j][0]*2f*gdata[l][k][j][0]*bcos[j]/(float)pow(rs[j],3.0);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(hdata[k][j][0][l]!=undef) fdata[k][j][0][l]=
			hdata[k][j][0][l]*2f*gdata[k][j][0][l]*bcos[j]/(float)pow(rs[j],3.0);
		}
		
		return ff;
	}
	
	/**
     * Calculate tangential frictional torque.
     *
     * @param	fm	tangential-mean tangential friction (m s^-2)
     *
     * @return	f	force (m^2 s^-2)
     */
	public Variable cFrictionalTorque(Variable fm){
		assignSubDomainParams(fm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable ft=new Variable("frtq",fm);
		ft.setCommentAndUnit("frictional torque (m^2 s^-2)");
		ft.setValue(undef);
		
		float[][][][] ftdata=ft.getData();
		float[][][][] frdata=fm.getData();
		
		if(fm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(frdata[l][k][j][0]!=undef) ftdata[l][k][j][0]=frdata[l][k][j][0]*rs[j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(frdata[k][j][0][l]!=undef) ftdata[k][j][0][l]=frdata[k][j][0][l]*rs[j];
		}
		
		return ft;
	}
	
	/**
     * Calculate eddy absolute angular momentum (AAM) horizontal flux convergence (HFC).
     *
     * @param	gava	eddy AAM horizontal flux = [g'v'] (m^3 s^-2)
     *
     * @return	f		eddy AAM HFC (m^2 s^-2)
     */
	public Variable cEddyAAMHFC(Variable gava){
		if(gava.getXCount()!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable aamHFC=cRadialDiv(gava).multiplyEq(-1f);
		aamHFC.setName("aamHFC");
		aamHFC.setCommentAndUnit("eddy absolute angular momentum horizontal flux convergence (m^2 s^-2)");
		
		return aamHFC;
	}
	
	/**
     * Calculate eddy absolute angular momentum (AAM) vertical flux convergence (VFC).
     *
     * @param	gawa	eddy AAM vertical flux = [g'w'] (m^2 Pa s^-2)
     *
     * @return	f		eddy AAM VFC (m^2 s^-2)
     */
	public Variable cEddyAAMVFC(Variable gawa){
		if(gawa.getXCount()!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable aamVFC=cVerticalDiv(gawa).multiplyEq(-1f);
		aamVFC.setName("aamVFC");
		aamVFC.setCommentAndUnit("eddy absolute angular momentum vertical flux convergence (m^2 s^-2)");
		
		return aamVFC;
	}
	
	/**
     * Calculate tilting effect.
     *
     * @param	gm	tangential-mean absolute angular momentum (m^2 s^-1)
     * @param	ut	tangential wind component (not storm-relative) (m s^-1)
     * @param	vr	radial wind component (not storm-relative) (m s^-1)
     *
     * @return	f	tilting effect (m^2 s^-2)
     */
	public Variable cTilting(Variable ut,Variable vr){
		assignSubDomainParams(ut);
		checkDimensions(ut,vr);
		
		Variable f=new Variable("tilt",ut.isTFirst(),new Range(t,z,y,1));
		f.setCommentAndUnit("tilting effect (m^2 s^-2)");
		f.setUndef(undef);
		f.setValue(undef);
		
    	float[] olat=((CylindricalSpatialModel)sm).getOLat();
		float[][][][]  fdata= f.getData();
		float[][][][] utdata=ut.getData();
		float[][][][] vrdata=vr.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++){
		    	float[][] buf=new float[z][y];
		    	
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(utdata[l][k][j][i]!=undef&&vrdata[l][k][j][i]!=undef){
						float til3=(1f-(float)cos(ydef[j]))*(wtan[l][i]*vrdata[l][k][j][i]+wnor[l][i]*utdata[l][k][j][i]);
						float til4=-vrdata[l][k][j][i]*bsin[j]*(float)(cos(ydef[j])*sin(olat[l])+bsin[j]*cos(olat[l])*cos(xdef[i])-sin(olat[l]));
						float til5=(1f-(float)cos(ydef[j]))*(float)(
							(utdata[l][k][j][i]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[l][k][j][i]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i])
						);
						buf[k][j]+=til3+EARTH_RADIUS*EARTH_ROTATE_SPEED*(til4+til5);  count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					
					float til2=-EARTH_ROTATE_SPEED*cs[l]*rs[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2;
					buf[k][j]+=til2;  count=0;
				}
				
				/*** assign ***/
				for(int j=1;j<y;j++)
				for(int k=0;k<z;k++) fdata[l][k][j][0]=buf[k][j];
			}
			
		}else{
			for(int l=0;l<t;l++){
		    	float[][] buf=new float[z][y];
		    	
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(utdata[k][j][i][l]!=undef&&vrdata[k][j][i][l]!=undef){
						float til3=(1f-(float)cos(ydef[j]))*(wtan[l][i]*vrdata[k][j][i][l]+wnor[l][i]*utdata[k][j][i][l]);
						float til4=-vrdata[k][j][i][l]*bsin[j]*(float)(cos(ydef[j])*sin(olat[l])+bsin[j]*cos(olat[l])*cos(xdef[i])-sin(olat[l]));
						float til5=(1f-(float)cos(ydef[j]))*(float)(
							(utdata[k][j][i][l]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[k][j][i][l]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i])
						);
						buf[k][j]+=til3+EARTH_RADIUS*EARTH_ROTATE_SPEED*(til4+til5);  count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					
					float til2=-EARTH_ROTATE_SPEED*cs[l]*rs[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2;
					buf[k][j]+=til2;  count=0;
				}
				
				/*** assign ***/
				for(int j=1;j<y;j++)
				for(int k=0;k<z;k++) fdata[k][j][0][l]=buf[k][j];
			}
		}
		
		return f;
	}
	
	
	/**
     * Calculate eddy-induced streamfunction.
     *
     * @param	tava	eddy heat horizontal flux = [¦È'v'] (K m s^-1)
     * @param	thm		azimuthal-mean potential temperature (K)
     *
     * @return	sfe		eddy-induced streamfunction (Pa m s^-1)
     */
	public Variable cEddyInducedStreamfunction(Variable tava,Variable thm){
		assignSubDomainParams(tava);
		checkDimensions(tava,thm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable sfe=new Variable("sfe",tava);
		sfe.setCommentAndUnit("eddy-induced streamfunction (Pa m s^-1)");
		sfe.setValue(undef);
		
		float[][][][] sfdata= sfe.getData();
		float[][][][] tvdata=tava.getData();
		float[][][][] tmdata= thm.getData();
		
		float threshold=1e-1f;
		
		if(tava.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0)
				if(tvdata[l][0][j][0]!=undef&&tmdata[l][0][j][0]!=undef&&tmdata[l][1][j][0]!=undef){
					float diffZ=tmdata[l][1][j][0]-tmdata[l][0][j][0];
					if(Math.abs(diffZ)>=threshold) sfdata[l][0][j][0]=tvdata[l][0][j][0]/diffZ*dz;
				}
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(tvdata[l][k][j][0]!=undef&&tmdata[l][k-1][j][0]!=undef&&tmdata[l][k+1][j][0]!=undef){
					float diffZ=tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0];
					if(Math.abs(diffZ)>=threshold) sfdata[l][k][j][0]=tvdata[l][k][j][0]/diffZ*(dz*2f);
				}
				
				// upper boundary (k==z-1)
				if(tvdata[l][z-1][j][0]!=undef&&tmdata[l][z-1][j][0]!=undef&&tmdata[l][z-2][j][0]!=undef){
					float diffZ=tmdata[l][z-1][j][0]-tmdata[l][z-2][j][0];
					if(Math.abs(diffZ)>=threshold) sfdata[l][z-1][j][0]=tvdata[l][z-1][j][0]/diffZ*dz;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0)
				if(tvdata[0][j][0][l]!=undef&&tmdata[0][j][0][l]!=undef&&tmdata[1][j][0][l]!=undef){
					float diffZ=tmdata[1][j][0][l]-tmdata[0][j][0][l];
					if(Math.abs(diffZ)>=threshold) sfdata[0][j][0][l]=tvdata[0][j][0][l]/diffZ*dz;
				}
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(tvdata[k][j][0][l]!=undef&&tmdata[k-1][j][0][l]!=undef&&tmdata[k+1][j][0][l]!=undef){
					float diffZ=tmdata[k+1][j][0][l]-tmdata[k-1][j][0][l];
					if(Math.abs(diffZ)>=threshold) sfdata[k][j][0][l]=tvdata[k][j][0][l]/diffZ*(dz*2f);
				}
				
				// upper boundary (k==z-1)
				if(tvdata[z-1][j][0][l]!=undef&&tmdata[z-1][j][0][l]!=undef&&tmdata[z-2][j][0][l]!=undef){
					float diffZ=tmdata[z-1][j][0][l]-tmdata[z-2][j][0][l];
					if(Math.abs(diffZ)>=threshold) sfdata[z-1][j][0][l]=tvdata[z-1][j][0][l]/diffZ*dz;
				}
			}
		}
		
		return sfe;
	}
	
	/**
     * Calculate coordinate-independent eddy-induced streamfunction.
     * Reference: Andrews and McIntyre 1978, JAS; Holton 1981, JGR
     *
     * @param	tava	eddy heat horizontal flux = [¦È'v'] (K m s^-1)
     * @param	tawa	eddy heat  vertical  flux = [¦È'¦Ø'] (K Pa s^-1)
     * @param	thm		azimuthal-mean potential temperature (K)
     *
     * @return	sfe		eddy-induced streamfunction (Pa m s^-1)
     */
	public Variable cCoordinateIndependentEddyInducedStreamfunction(Variable tava,Variable tawa,Variable thm){
		assignSubDomainParams(tava);
		checkDimensions(tava,tawa,thm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable sfe=new Variable("sfe",tava);
		sfe.setCommentAndUnit("coordinate-independent eddy-induced streamfunction (Pa m s^-1)");
		sfe.setValue(undef);
		
		float[][][][] sfdata= sfe.getData();
		float[][][][] tvdata=tava.getData();
		float[][][][] twdata=tawa.getData();
		float[][][][] tmdata= thm.getData();
		
		float threshold=0.05f;
		
		if(tava.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1,J=y-1;j<J;j++){
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(tvdata[l][k][j][0]!=undef&&tmdata[l][k-1][j][0]!=undef&&tmdata[l][k+1][j][0]!=undef
				&&tmdata[l][k][j-1][0]!=undef&&tmdata[l][k][j+1][0]!=undef){
					double diffZ=((double)tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0])/(dz*2.0);
					double diffY=((double)tmdata[l][k][j+1][0]-tmdata[l][k][j-1][0])/(dy*2.0);
					double gradT2=diffZ*diffZ+diffY*diffY;
					
					if(Math.abs(gradT2)>=threshold)
					sfdata[l][k][j][0]=(float)((tvdata[l][k][j][0]*diffZ-twdata[l][k][j][0]*diffY)/gradT2);
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1,J=y-1;j<J;j++){
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(tvdata[k][j][0][l]!=undef&&tmdata[k-1][j][0][l]!=undef&&tmdata[k+1][j][0][l]!=undef
				&&tmdata[k][j-1][0][l]!=undef&&tmdata[k][j+1][0][l]!=undef){
					double diffZ=((double)tmdata[k+1][j][0][l]-tmdata[k-1][j][0][l])/(dz*2.0);
					double diffY=((double)tmdata[k][j+1][0][l]-tmdata[k][j-1][0][l])/(dy*2.0);
					double gradT2=diffZ*diffZ+diffY*diffY;
					
					if(Math.abs(gradT2)>=threshold)
					sfdata[k][j][0][l]=(float)((tvdata[k][j][0][l]*diffZ-twdata[k][j][0][l]*diffY)/gradT2);
				}
			}
		}
		
		return sfe;
	}
	
	
    /**
     * Calculate v and w, [0] is in y direction while [1] is in z direction
     *
     * @param	sf	streamfunction in the radial-vertical plane (m Pa s^-1)
     *
     * @return	vw	[0] is radial velocity and [1] is vertical velocity
     */
	public Variable[] cVW(Variable sf){
		assignSubDomainParams(sf);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		Variable[] vw=new Variable[2];
		
		vw[0]=new Variable("vs",sf);	vw[0].setValue(undef);
		vw[1]=new Variable("ws",sf);	vw[1].setValue(undef);
		vw[0].setCommentAndUnit("simulated radial velocity (m s^-1)");
		vw[1].setCommentAndUnit("simulated vertical velocity (Pa s^-1)");
		
		float[][][][] sfdata=   sf.getData();
		float[][][][] vsdata=vw[0].getData();
		float[][][][] wsdata=vw[1].getData();
		
		if(sf.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=1,K=z-1;k<K;k++)
				for(int j=1;j<y;j++)
				if(sfdata[l][k+1][j][0]!=undef&&sfdata[l][k-1][j][0]!=undef)
				vsdata[l][k][j][0]=(sfdata[l][k+1][j][0]-sfdata[l][k-1][j][0])/bsin[j]/(dz*2);
				
				/*** top and bottom boundary **
				for(int j=1;j<y;j++){
					if(sfdata[l][z-1][j][0]!=undef&&sfdata[l][z-2][j][0]!=undef)
					vsdata[l][z-1][j][0]=2*(sfdata[l][z-1][j][0]-sfdata[l][z-2][j][0])/bsin[j]/dz[0]-vsdata[l][z-2][j][0];
					
					if(sfdata[l][0][j][0]!=undef&&sfdata[l][1][j][0]!=undef)
					vsdata[l][0][j][0]=2*(sfdata[l][1][j][0]-sfdata[l][0][j][0])/bsin[j]/dz[0]-vsdata[l][1][j][0];
				}*/
				
				for(int k=0;k<z;k++)
				for(int j=1,J=y-1;j<J;j++)
				if(sfdata[l][k][j+1][0]!=undef&&sfdata[l][k][j-1][0]!=undef)
				wsdata[l][k][j][0]=-(sfdata[l][k][j+1][0]-sfdata[l][k][j-1][0])/bsin[j]/(dy*2);
				
				/*** left and right boundary **
				for(int k=0;k<z;k++){
					if(sfdata[l][k][y-1][0]!=undef&&sfdata[l][k][y-2][0]!=undef)
					wsdata[l][k][y-1][0]=-2*(sfdata[l][k][y-1][0]-sfdata[l][k][y-2][0])/
					(float)sin((ydef[y-1]+ydef[y-2])/2)/dy-wsdata[l][k][y-2][0];
					
					if(sfdata[l][k][1][0]!=undef&&sfdata[l][k][0][0]!=undef)
					wsdata[l][k][0][0]=-2*(sfdata[l][k][1][0]-sfdata[l][k][0][0])/
					(float)sin((ydef[0]+ydef[1])/2)/dy-wsdata[l][k][1][0];
				}*/
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=1,K=z-1;k<K;k++)
				for(int j=1;j<y;j++)
				if(sfdata[k+1][j][0][l]!=undef&&sfdata[k-1][j][0][l]!=undef)
				vsdata[k][j][0][l]=(sfdata[k+1][j][0][l]-sfdata[k-1][j][0][l])/bsin[j]/(dz*2);
				
				/*** top and bottom boundary **
				for(int j=1;j<y;j++){
					if(sfdata[z-1][j][0][l]!=undef&&sfdata[z-2][j][0][l]!=undef)
					vsdata[z-1][j][0][l]=2*(sfdata[z-1][j][0][l]-sfdata[z-2][j][0][l])/bsin[j]/dz[0]-vsdata[z-2][j][0][l];
					
					if(sfdata[0][j][0][l]!=undef&&sfdata[1][j][0][l]!=undef)
					vsdata[0][j][0][l]=2*(sfdata[1][j][0][l]-sfdata[0][j][0][l])/bsin[j]/dz[0]-vsdata[1][j][0][l];
				}*/
				
				for(int k=0;k<z;k++)
				for(int j=1,J=y-1;j<J;j++)
				if(sfdata[k][j+1][0][l]!=undef&&sfdata[k][j-1][0][l]!=undef)
				wsdata[k][j][0][l]=-(sfdata[k][j+1][0][l]-sfdata[k][j-1][0][l])/bsin[j]/(dy*2);
				
				/*** left and right boundary **
				for(int k=0;k<z;k++){
					if(sfdata[k][y-1][0][l]!=undef&&sfdata[k][y-2][0][l]!=undef)
					wsdata[k][y-1][0][l]=-2*(sfdata[k][y-1][0][l]-sfdata[k][y-2][0][l])/
					(float)sin((ydef[y-1]+ydef[y-2])/2)/dy-wsdata[k][y-2][0][l];
					
					if(sfdata[k][1][0][l]!=undef&&sfdata[k][0][0][l]!=undef)
					wsdata[k][0][0][l]=-2*(sfdata[k][1][0][l]-sfdata[k][0][0][l])/
					(float)sin((ydef[0]+ydef[1])/2)/dy-wsdata[k][1][0][l];
				}*/
			}
		}
		
		return vw;
	}
	
	/**
     * Calculate eddy-induced velocity according to eddy-induced streamfunction.
     *
     * @param	sfe		eddy-induced streamfunction (Pa m s^-2)
     *
     * @return	vedy	eddy-induced velocity, [0] is meridional (m s^-1) and [1] is vertical (Pa s^-1)
     */
	public Variable[] cEddyInducedVelocity(Variable sfe){
		assignSubDomainParams(sfe);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable[] vedy=new Variable[2];
		vedy[0]=new Variable("vedy",sfe); vedy[0].setValue(undef);
		vedy[1]=new Variable("wedy",sfe); vedy[1].setValue(undef);
		vedy[0].setCommentAndUnit("meridional component of eddy-induced velocity (m s^-1)");
		vedy[1].setCommentAndUnit("vertical component of eddy-induced velocity (Pa s^-1)");
		
		float[][][][] vdata=vedy[0].getData();
		float[][][][] wdata=vedy[1].getData();
		float[][][][] sdata=    sfe.getData();
		
		if(sfe.isTFirst()){
			for(int l=0;l<t;l++){
				/*** meridional component ***/
				for(int j=0;j<y;j++){
					// lower boundary (k==0), not used in inversion
					if(sdata[l][0][j][0]!=undef&&sdata[l][1][j][0]!=undef)
					vdata[l][0][j][0]=(sdata[l][1][j][0]-sdata[l][0][j][0])/dz;
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(sdata[l][k-1][j][0]!=undef&&sdata[l][k+1][j][0]!=undef)
					vdata[l][k][j][0]=(sdata[l][k+1][j][0]-sdata[l][k-1][j][0])/(dz*2f);
					
					// upper boundary (k==z-1), not used in inversion
					if(sdata[l][z-1][j][0]!=undef&&sdata[l][z-2][j][0]!=undef)
					vdata[l][z-1][j][0]=(sdata[l][z-1][j][0]-sdata[l][z-2][j][0])/dz;
				}
				
				/*** vertical component ***/
				for(int k=0;k<z;k++){
					// south boundary (j==0), not used in inversion
					if(sdata[l][k][0][0]!=undef&&sdata[l][k][1][0]!=undef)
					wdata[l][k][0][0]=-(sdata[l][k][1][0]*bsin[1]-sdata[l][k][0][0]*bsin[0])
					/(dy*2f)/((bsin[0]+bsin[1])/2f);
					
					// inner region, centered difference
					for(int j=1,J=y-1;j<J;j++)
					if(sdata[l][k][j-1][0]!=undef&&sdata[l][k][j+1][0]!=undef)
					wdata[l][k][j][0]=-(sdata[l][k][j+1][0]*bsin[j+1]-sdata[l][k][j-1][0]*bsin[j-1])
					/(dy*2f)/bsin[j];
					
					// north boundary (j==y-1), not used in inversion
					if(sdata[l][k][y-1][0]!=undef&&sdata[l][k][y-2][0]!=undef)
					wdata[l][k][y-1][0]=-(sdata[l][k][y-1][0]*bsin[y-1]-sdata[l][k][y-2][0]*bsin[y-2])
					/(dy*2f)/((bsin[y-1]+bsin[y-2])/2f);
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				/*** meridional component ***/
				for(int j=0;j<y;j++){
					// lower boundary (k==0), not used in inversion
					if(sdata[0][j][0][l]!=undef&&sdata[1][j][0][l]!=undef)
					vdata[0][j][0][l]=(sdata[1][j][0][l]-sdata[0][j][0][l])/dz;
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(sdata[k-1][j][0][l]!=undef&&sdata[k+1][j][0][l]!=undef)
					vdata[k][j][0][l]=(sdata[k+1][j][0][l]-sdata[k-1][j][0][l])/(dz*2f);
					
					// upper boundary (k==z-1), not used in inversion
					if(sdata[z-1][j][0][l]!=undef&&sdata[z-2][j][0][l]!=undef)
					vdata[z-1][j][0][l]=(sdata[z-1][j][0][l]-sdata[z-2][j][0][l])/dz;
				}
				
				/*** vertical component ***/
				for(int k=0;k<z;k++){
					// south boundary (j==0), not used in inversion
					if(sdata[k][0][0][l]!=undef&&sdata[k][1][0][l]!=undef)
					wdata[k][0][0][l]=-(sdata[k][1][0][l]*bsin[1]-sdata[k][0][0][l]*bsin[0])
					/(dy*2f)/((bsin[0]+bsin[1])/2f);
					
					// inner region, centered difference
					for(int j=1,J=y-1;j<J;j++)
					if(sdata[k][j-1][0][l]!=undef&&sdata[k][j+1][0][l]!=undef)
					wdata[k][j][0][l]=-(sdata[k][j+1][0][l]*bsin[j+1]-sdata[k][j-1][0][l]*bsin[j-1])
					/(dy*2f)/bsin[j];
					
					// north boundary (j==y-1), not used in inversion
					if(sdata[k][y-1][0][l]!=undef&&sdata[k][y-2][0][l]!=undef)
					wdata[k][y-1][0][l]=-(sdata[k][y-1][0][l]*bsin[y-1]-sdata[k][y-2][0][l]*bsin[y-2])
					/(dy*2f)/((bsin[y-1]+bsin[y-2])/2f);
				}
			}
		}
		
		return vedy;
	}
	
	/**
     * Calculate Eliassen-Palm (EP) vector.
     *
     * @param	tava	eddy heat horizontal flux = [Î¸'v'] (K m s^-1)
     * @param	gava	eddy AAM  horizontal flux = [g'v'] (m^3 s^-2)
     * @param	gawa	eddy AAM vertical flux = [g'w'] (m^2 Pa s^-2)
     * @param	gm		azimuthal-mean absolute angular momentum (m^2 s^-1)
     * @param	thm		azimuthal-mean potential temperature (K)
     * @param	rZToY	ratio of Z to Y: Z/Y in the plot used to scale the vector components (e.g., 0.8)
     *
     * @return	EP		EP vector, [0] meridional (m^3 s^-2) and [1] vertical components (Pa m^2 s^-2);
     * 					[2] and [3] are the same as [0] and [1] except for visual scaling.
     */
	public Variable[] cEPVector(Variable tava,Variable gava,Variable gawa,Variable gm,Variable thm,float rZToY){
		assignSubDomainParams(tava);
		checkDimensions(tava,gava,gawa,gm,thm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable[] EP=new Variable[4];
		EP[0]=new Variable("EPy",tava); EP[0].setValue(undef);
		EP[1]=new Variable("EPz",tava); EP[1].setValue(undef);
		EP[0].setCommentAndUnit("meridional component of Eliassen-Palm vector (m^3 s^-2)");
		EP[1].setCommentAndUnit("vertical component of Eliassen-Palm vector (Pa m^2 s^-2)");
		
		float[][][][] mcdata=EP[0].getData();
		float[][][][] vcdata=EP[1].getData();
		float[][][][] tvdata= tava.getData();
		float[][][][] gvdata= gava.getData();
		float[][][][] gwdata= gawa.getData();
		float[][][][] gmdata=   gm.getData();
		float[][][][] tmdata=  thm.getData();
		
		if(tava.isTFirst()){
			for(int l=0;l<t;l++){
				/*** meridional component ***/
				for(int j=0;j<y;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[l][0][j][0]!=undef&&gwdata[l][0][j][0]!=undef&&tvdata[l][0][j][0]!=undef)
					mcdata[l][0][j][0]=-gvdata[l][0][j][0]+tvdata[l][0][j][0]*
					(gmdata[l][1][j][0]-gmdata[l][0][j][0])/(tmdata[l][1][j][0]-tmdata[l][0][j][0]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[l][k][j][0]!=undef&&gwdata[l][k][j][0]!=undef&&tvdata[l][k][j][0]!=undef)
					mcdata[l][k][j][0]=-gvdata[l][k][j][0]+tvdata[l][k][j][0]*
					(gmdata[l][k+1][j][0]-gmdata[l][k-1][j][0])/(tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[l][z-1][j][0]!=undef&&gwdata[l][z-1][j][0]!=undef&&tvdata[l][z-1][j][0]!=undef)
					mcdata[l][z-1][j][0]=-gvdata[l][z-1][j][0]+tvdata[l][z-1][j][0]*
					(gmdata[l][z-1][j][0]-gmdata[l][z-2][j][0])/(tmdata[l][z-1][j][0]-tmdata[l][z-2][j][0]);
				}
				
				/*** vertical component ***/
				// lower south boundary (k==0, j==0), not used in inversion
				if(gvdata[l][0][0][0]!=undef&&gwdata[l][0][0][0]!=undef&&tvdata[l][0][0][0]!=undef)
				vcdata[l][0][0][0]=-gwdata[l][0][0][0]-tvdata[l][0][0][0]*dz/dy*
				(gmdata[l][0][1][0]-gmdata[l][0][0][0])/(tmdata[l][1][0][0]-tmdata[l][0][0][0]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[l][k][0][0]!=undef&&gwdata[l][k][0][0]!=undef&&tvdata[l][k][0][0]!=undef)
				vcdata[l][k][0][0]=-gwdata[l][k][0][0]-tvdata[l][k][0][0]*2f*dz/dy*
				(gmdata[l][k][1][0]-gmdata[l][k][0][0])/(tmdata[l][k+1][0][0]-tmdata[l][k-1][0][0]);
				
				// upper south boundary (k==z-1, j==0), not used in inversion
				if(gvdata[l][z-1][0][0]!=undef&&gwdata[l][z-1][0][0]!=undef&&tvdata[l][z-1][0][0]!=undef)
				vcdata[l][z-1][0][0]=-gwdata[l][z-1][0][0]-tvdata[l][z-1][0][0]*dz/dy*
				(gmdata[l][z-1][1][0]-gmdata[l][z-1][0][0])/(tmdata[l][z-1][0][0]-tmdata[l][z-2][0][0]);
				
				for(int j=1,J=y-1;j<J;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[l][0][j][0]!=undef&&gwdata[l][0][j][0]!=undef&&tvdata[l][0][j][0]!=undef)
					vcdata[l][0][j][0]=-gwdata[l][0][j][0]-tvdata[l][0][j][0]*dz/dy/2f*
					(gmdata[l][0][j+1][0]-gmdata[l][0][j-1][0])/(tmdata[l][1][j][0]-tmdata[l][0][j][0]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[l][k][j][0]!=undef&&gwdata[l][k][j][0]!=undef&&tvdata[l][k][j][0]!=undef)
					vcdata[l][k][j][0]=-gwdata[l][k][j][0]-tvdata[l][k][j][0]*dz/dy*
					(gmdata[l][k][j+1][0]-gmdata[l][k][j-1][0])/(tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[l][z-1][j][0]!=undef&&gwdata[l][z-1][j][0]!=undef&&tvdata[l][z-1][j][0]!=undef)
					vcdata[l][z-1][j][0]=-gwdata[l][z-1][j][0]-tvdata[l][z-1][j][0]*dz/dy/2f*
					(gmdata[l][z-1][j+1][0]-gmdata[l][z-1][j-1][0])/(tmdata[l][z-1][j][0]-tmdata[l][z-2][j][0]);
				}
				
				// lower north boundary (k==0,j==y-1), not used in inversion
				if(gvdata[l][0][y-1][0]!=undef&&gwdata[l][0][y-1][0]!=undef&&tvdata[l][0][y-1][0]!=undef)
				vcdata[l][0][y-1][0]=-gwdata[l][0][y-1][0]-tvdata[l][0][y-1][0]*dz/dy*
				(gmdata[l][0][y-1][0]-gmdata[l][0][y-2][0])/(tmdata[l][1][y-1][0]-tmdata[l][0][y-1][0]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[l][k][y-1][0]!=undef&&gwdata[l][k][y-1][0]!=undef&&tvdata[l][k][y-1][0]!=undef)
				vcdata[l][k][y-1][0]=-gwdata[l][k][y-1][0]-tvdata[l][k][y-1][0]*2f*dz/dy*
				(gmdata[l][k][y-1][0]-gmdata[l][k][y-2][0])/(tmdata[l][k+1][y-1][0]-tmdata[l][k-1][y-1][0]);
				
				// upper north boundary (k==z-1,j==y-1), not used in inversion
				if(gvdata[l][z-1][y-1][0]!=undef&&gwdata[l][z-1][y-1][0]!=undef&&tvdata[l][z-1][y-1][0]!=undef)
				vcdata[l][z-1][y-1][0]=-gwdata[l][z-1][y-1][0]-tvdata[l][z-1][y-1][0]*dz/dy*
				(gmdata[l][z-1][y-1][0]-gmdata[l][z-1][y-2][0])/(tmdata[l][z-1][y-1][0]-tmdata[l][z-2][y-1][0]);
			}
			
		}else{
			for(int l=0;l<t;l++){
				/*** meridional component ***/
				for(int j=0;j<y;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[0][j][0][l]!=undef&&gwdata[0][j][0][l]!=undef&&tvdata[0][j][0][l]!=undef)
					mcdata[0][j][0][l]=-gvdata[0][j][0][l]+tvdata[0][j][0][l]*
					(gmdata[1][j][0][l]-gmdata[0][j][0][l])/(tmdata[1][j][0][l]-tmdata[0][j][0][l]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[k][j][0][l]!=undef&&gwdata[k][j][0][l]!=undef&&tvdata[k][j][0][l]!=undef)
					mcdata[k][j][0][l]=-gvdata[k][j][0][l]+tvdata[k][j][0][l]*
					(gmdata[k+1][j][0][l]-gmdata[k-1][j][0][l])/(tmdata[k+1][j][0][l]-tmdata[k-1][j][0][l]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[z-1][j][0][l]!=undef&&gwdata[z-1][j][0][l]!=undef&&tvdata[z-1][j][0][l]!=undef)
					mcdata[z-1][j][0][l]=-gvdata[z-1][j][0][l]+tvdata[z-1][j][0][l]*
					(gmdata[z-1][j][0][l]-gmdata[z-2][j][0][l])/(tmdata[z-1][j][0][l]-tmdata[z-2][j][0][l]);
				}
				
				/*** vertical component ***/
				// lower south boundary (k==0, j==0), not used in inversion
				if(gvdata[0][0][0][l]!=undef&&gwdata[0][0][0][l]!=undef&&tvdata[0][0][0][l]!=undef)
				vcdata[0][0][0][l]=-gwdata[0][0][0][l]-tvdata[0][0][0][l]*dz/dy*
				(gmdata[0][1][0][l]-gmdata[0][0][0][l])/(tmdata[1][0][0][l]-tmdata[0][0][0][l]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[k][0][0][l]!=undef&&gwdata[k][0][0][l]!=undef&&tvdata[k][0][0][l]!=undef)
				vcdata[k][0][0][l]=-gwdata[k][0][0][l]-tvdata[k][0][0][l]*2f*dz/dy*
				(gmdata[k][1][0][l]-gmdata[k][0][0][l])/(tmdata[k+1][0][0][l]-tmdata[k-1][0][0][l]);
				
				// upper south boundary (k==z-1, j==0), not used in inversion
				if(gvdata[z-1][0][0][l]!=undef&&gwdata[z-1][0][0][l]!=undef&&tvdata[z-1][0][0][l]!=undef)
				vcdata[z-1][0][0][l]=-gwdata[z-1][0][0][l]-tvdata[z-1][0][0][l]*dz/dy*
				(gmdata[z-1][1][0][l]-gmdata[z-1][0][0][l])/(tmdata[z-1][0][0][l]-tmdata[z-2][0][0][l]);
				
				for(int j=1,J=y-1;j<J;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[0][j][0][l]!=undef&&gwdata[0][j][0][l]!=undef&&tvdata[0][j][0][l]!=undef)
					vcdata[0][j][0][l]=-gwdata[0][j][0][l]-tvdata[0][j][0][l]*dz/dy/2f*
					(gmdata[0][j+1][0][l]-gmdata[0][j-1][0][l])/(tmdata[1][j][0][l]-tmdata[0][j][0][l]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[k][j][0][l]!=undef&&gwdata[k][j][0][l]!=undef&&tvdata[k][j][0][l]!=undef)
					vcdata[k][j][0][l]=-gwdata[k][j][0][l]-tvdata[k][j][0][l]*dz/dy*
					(gmdata[k][j+1][0][l]-gmdata[k][j-1][0][l])/(tmdata[k+1][j][0][l]-tmdata[k-1][j][0][l]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[z-1][j][0][l]!=undef&&gwdata[z-1][j][0][l]!=undef&&tvdata[z-1][j][0][l]!=undef)
					vcdata[z-1][j][0][l]=-gwdata[z-1][j][0][l]-tvdata[z-1][j][0][l]*dz/dy/2f*
					(gmdata[z-1][j+1][0][l]-gmdata[z-1][j-1][0][l])/(tmdata[z-1][j][0][l]-tmdata[z-2][j][0][l]);
				}
				
				// lower north boundary (k==0,j==y-1), not used in inversion
				if(gvdata[0][y-1][0][l]!=undef&&gwdata[0][y-1][0][l]!=undef&&tvdata[0][y-1][0][l]!=undef)
				vcdata[0][y-1][0][l]=-gwdata[0][y-1][0][l]-tvdata[0][y-1][0][l]*dz/dy*
				(gmdata[0][y-1][0][l]-gmdata[0][y-2][0][l])/(tmdata[1][y-1][0][l]-tmdata[0][y-1][0][l]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[k][y-1][0][l]!=undef&&gwdata[k][y-1][0][l]!=undef&&tvdata[k][y-1][0][l]!=undef)
				vcdata[k][y-1][0][l]=-gwdata[k][y-1][0][l]-tvdata[k][y-1][0][l]*2f*dz/dy*
				(gmdata[k][y-1][0][l]-gmdata[k][y-2][0][l])/(tmdata[k+1][y-1][0][l]-tmdata[k-1][y-1][0][l]);
				
				// upper north boundary (k==z-1,j==y-1), not used in inversion
				if(gvdata[z-1][y-1][0][l]!=undef&&gwdata[z-1][y-1][0][l]!=undef&&tvdata[z-1][y-1][0][l]!=undef)
				vcdata[z-1][y-1][0][l]=-gwdata[z-1][y-1][0][l]-tvdata[z-1][y-1][0][l]*dz/dy*
				(gmdata[z-1][y-1][0][l]-gmdata[z-1][y-2][0][l])/(tmdata[z-1][y-1][0][l]-tmdata[z-2][y-1][0][l]);
			}
		}
		
		EP[2]=EP[0].copy(); EP[2].setName("EPyS");
		EP[3]=EP[1].copy(); EP[3].setName("EPzS");
		EP[2].setCommentAndUnit("scaled EPy by multiplying actual length representing 1m (m^3 s^-2)");
		EP[3].setCommentAndUnit("scaled EPz by multiplying actual length representing 1Pa (Pa m^2 s^-2)");
		
		float[][][][] msdata=EP[2].getData();
		float[][][][] vsdata=EP[3].getData();
		
		float denomZ=Math.abs(zdef[z-1]-zdef[0])/rZToY;
		float denomY=Math.abs(ydef[y-1]-ydef[0])*EARTH_RADIUS;
		
		//System.out.println(dz);
		//System.out.println(dy);
		
		//System.out.println(denomZ+"\t"+denomY+"\t"+rZToY);
		
		if(tava.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(mcdata[l][k][j][0]!=undef) msdata[l][k][j][0]=mcdata[l][k][j][0]/denomY;
				if(vcdata[l][k][j][0]!=undef) vsdata[l][k][j][0]=vcdata[l][k][j][0]/denomZ;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(mcdata[k][j][0][l]!=undef) msdata[k][j][0][l]=mcdata[k][j][0][l]/denomY;
				if(vcdata[k][j][0][l]!=undef) vsdata[k][j][0][l]=vcdata[k][j][0][l]/denomZ;
			}
		}
		
		return EP;
	}
	
	
	/**
     * Calculate radial-vertical (YZ) plane divergence components and total divergence.
     *
     * @param	vecY	azimuthal-mean   radial component of a vector
     * @param	vecZ	azimuthal-mean vertical component of a vector
     *
     * @return	div		divergence in YZ plane (vecY * m^-1 or vecZ * Pa^-1)
     * 					[0] is y-component, [1] z-component, and [2] total (sum of [0] amd [1]) divergence
     */
	public Variable[] cYZDivergence(Variable vecY,Variable vecZ){
		assignSubDomainParams(vecY);
		checkDimensions(vecY,vecZ);
		
		Variable divY=cRadialDiv(vecY);
		Variable divZ=cVerticalDiv(vecZ);
		Variable div =divY.plus(divZ);
		
		divY.setName("divY"); divY.setComment("y-comp divergence in YZ plane (vecY * m^-1 or vecZ * Pa^-1)");
		divZ.setName("divZ"); divZ.setComment("z-comp divergence in YZ plane (vecY * m^-1 or vecZ * Pa^-1)");
		div.setName("divYZ");  div.setComment("total  divergence in YZ plane (vecY * m^-1 or vecZ * Pa^-1)");
		
		return new Variable[]{divY,divZ,div};
	}
	
	
	/**
     * Calculate pseudo density sigma = -dp/d¦È/g.
     *
     * @param	pres	pressure (Pa)
     *
     * @return	sigma	pseudo density (kg m^-2 K^-1)
     */
	public Variable cPseudoDensity(Variable pres){
		assignSubDomainParams(pres);
		
		Variable sigma=new Variable("sigma",pres);
		sigma.setCommentAndUnit("pseudo density (kg m^-2 K^-1)");
		sigma.setValue(undef);
		
		float[][][][] sdata=sigma.getData();
		float[][][][] pdata= pres.getData();
		
		if(pres.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				sdata[l][0][j][i]=-(pdata[l][1][j][i]-pdata[l][0][j][i])/dz/SpatialModel.GRAVITY_ACCERLERATION;
				for(int k=1,K=z-1;k<K;k++) sdata[l][k][j][i]=-(pdata[l][k+1][j][i]-pdata[l][k-1][j][i])/(2f*dz)/SpatialModel.GRAVITY_ACCERLERATION;
				sdata[l][z-1][j][i]=-(pdata[l][z-1][j][i]-pdata[l][z-2][j][i])/dz/SpatialModel.GRAVITY_ACCERLERATION;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				sdata[0][j][i][l]=-(pdata[1][j][i][l]-pdata[0][j][i][l])/dz/SpatialModel.GRAVITY_ACCERLERATION;
				for(int k=1,K=z-1;k<K;k++) sdata[k][j][i][l]=-(pdata[k+1][j][i][l]-pdata[k-1][j][i][l])/(2f*dz)/SpatialModel.GRAVITY_ACCERLERATION;
				sdata[z-1][j][i][l]=-(pdata[z-1][j][i][l]-pdata[z-2][j][i][l])/dz/SpatialModel.GRAVITY_ACCERLERATION;
			}
		}
		
		return sigma;
	}
	
	/**
     * Calculate the pseudo-density-weighted zonal mean and corresponding
     * anomaly (side-effect on the given var).
     *
     * @param	var		a given variable
     *
     * @return	den		pseudo density (kg m^-2 K^-1)
     */
	public Variable weightedAnomalizeX(Variable var,Variable den){
		assignSubDomainParams(var);
		checkDimensions(var,den);
		
		Variable vmm=new Variable(var.getName()+"mm",var.isTFirst(),new Range(t,z,y,1));
		vmm.setCommentAndUnit("pseudo-density-weighted zonal mean of "+var.getCommentAndUnit());
		vmm.setUndef(undef);
		vmm.setValue(undef);
		
		float[][][][] vrdata=var.getData();
		float[][][][] vmdata=vmm.getData();
		float[][][][] dsdata=den.getData();
		
		if(var.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				double sumVar=0,sumDen=0;
				
				for(int i=0;i<x;i++)
				if(vrdata[l][k][j][i]!=undef&&dsdata[l][k][j][i]!=undef){
					sumVar+=vrdata[l][k][j][i]*dsdata[l][k][j][i];
					sumDen+=dsdata[l][k][j][i];
				}
				
				if(sumDen!=0) vmdata[l][k][j][0]=(float)(sumVar/sumDen);
				
				for(int i=0;i<x;i++) vrdata[l][k][j][i]-=vmdata[l][k][j][0];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				double sumVar=0,sumDen=0;
				
				for(int i=0;i<x;i++)
				if(vrdata[k][j][i][l]!=undef&&dsdata[k][j][i][l]!=undef){
					sumVar+=vrdata[k][j][i][l]*dsdata[k][j][i][l];
					sumDen+=dsdata[k][j][i][l];
				}
				
				if(sumDen!=0) vmdata[k][j][0][l]=(float)(sumVar/sumDen);
				
				for(int i=0;i<x;i++) vrdata[k][j][i][l]-=vmdata[k][j][0][l];
			}
		}
		
		return vmm;
	}
	
	/**
     * Calculate pseudo-density-weighted eddy flux of A and B = [¦ÒA*B*] = [¦Ò]<A*B*>.
     *
     * @param	Ama		a given variable
     * @param	Bma		a given variable
     * @param	sigma	pseudo density (kg m^-2 K^-1)
     *
     * @return	eddy flux of A and B defined as [¦ÒA*B*] or [¦Ò]<A*B*>
     */
	public Variable cMassWeightedEddyFlux(Variable Ama,Variable Bma,Variable sigma){
		assignSubDomainParams(Ama);
		checkDimensions(Ama,Bma,sigma);
		
		Variable flx=sigma.multiply(Ama).multiplyEq(Bma).anomalizeX();
		
		flx.setName("flx");
		flx.setComment("pseudo-density-weighted eddy flux of A and B defined as [¦ÒA*B*] or [¦Ò]<A*B*>");
		
		return flx;
	}
	
	/**
     * Calculate azimuthal derivative of Montgomery streamfunction.
     *
     * @param	mont	Montgomery streamfunction (m^2 s^-2)
     *
     * @return	Mlambda	azimuthal derivative of Montgomery streamfunction (m^2 s^-2)
     */
	public Variable cMontLambda(Variable mont){
		assignSubDomainParams(mont);
		
		Variable Mld=new Variable("Mld",mont);
		Mld.setCommentAndUnit("azimuthal derivative of Montgomery streamfunction (m^2 s^-2)");
		Mld.setValue(undef);
		
		float delLambda=xdef[1]-xdef[0];
		float[][][][] sdata=  Mld.getData();
		float[][][][] mdata=mont.getData();
		
		if(mont.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				sdata[l][k][j][0]=(mdata[l][k][j][1]-mdata[l][k][j][x-1])/(delLambda*2);
				for(int i=1,I=x-1;i<I;i++) sdata[l][k][j][i]=(mdata[l][k][j][i+1]-mdata[l][k][j][i-1])/(delLambda*2);
				sdata[l][k][j][x-1]=(mdata[l][k][j][0]-mdata[l][k][j][x-2])/(delLambda*2);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				sdata[k][j][0][l]=(mdata[k][j][1][l]-mdata[k][j][x-1][l])/(delLambda*2);
				for(int i=1,I=x-1;i<I;i++) sdata[k][j][i][l]=(mdata[k][j][i+1][l]-mdata[k][j][i-1][l])/(delLambda*2);
				sdata[k][j][x-1][l]=(mdata[k][j][0][l]-mdata[k][j][x-2][l])/(delLambda*2);
			}
		}
		
		return Mld;
	}
	
	/**
     * Calculate small-amplitude Eliassen-Palm (EP) vector in isentropic coordinate.
     *
     * @param	sVaMa	eddy horizontal absolute angular momentum flux = [(sv)'M'] (kg m K^-1 s^-2)
     * @param	paMa	eddy horizontal heat flux = [p'M¦Ë'] (Pa m^2 s^-2)
     * @param	rZToY	ratio of Z to Y: Z/Y in the plot used to scale the vector components (e.g., 0.8)
     *
     * @return	EP		EP vector, [0] meridional (kg m K^-1 s^-2) and [1] vertical components (Pa m);
     * 					[2] and [3] are the same as [0] and [1] except for visual scaling according
     * 					to Edmon et al. (1980, JAS).
     */
	public Variable[] cIsentropicEPVectorByM(Variable sVaMa,Variable paMa,float rZToY){
		assignSubDomainParams(sVaMa);
		checkDimensions(sVaMa,paMa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable[] EP=new Variable[4];
		EP[0]=new Variable("EPyM",sVaMa); EP[0].setValue(undef);
		EP[1]=new Variable("EPzM",sVaMa); EP[1].setValue(undef);
		EP[0].setCommentAndUnit("meridional component of Eliassen-Palm vector EPy (kg m K^-1 s^-2)");
		EP[1].setCommentAndUnit("vertical component of Eliassen-Palm vector EPz (Pa m)");
		
		float[][][][] mcdata=EP[0].getData();
		float[][][][] vcdata=EP[1].getData();
		float[][][][] sadata=sVaMa.getData();
		float[][][][] pmdata= paMa.getData();
		
		if(sVaMa.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(sadata[l][k][j][0]!=undef) mcdata[l][k][j][0]=-sadata[l][k][j][0];
				if(pmdata[l][k][j][0]!=undef) vcdata[l][k][j][0]= pmdata[l][k][j][0]/SpatialModel.GRAVITY_ACCERLERATION;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(sadata[k][j][0][l]!=undef) mcdata[k][j][0][l]=-sadata[k][j][0][l];
				if(pmdata[k][j][0][l]!=undef) vcdata[k][j][0][l]= pmdata[k][j][0][l]/SpatialModel.GRAVITY_ACCERLERATION;
			}
		}
		
		EP[2]=EP[0].copy(); EP[2].setName("EPySM");
		EP[3]=EP[1].copy(); EP[3].setName("EPzSM");
		EP[2].setCommentAndUnit("scaled EPy by multiplying actual length representing 1m (kg m K^-1 s^-2)");
		EP[3].setCommentAndUnit("scaled EPz by multiplying actual length representing 1K (Pa m)");
		
		float[][][][] msdata=EP[2].getData();
		float[][][][] vsdata=EP[3].getData();
		
		float denomZ=Math.abs(zdef[z-1]-zdef[0])/rZToY;
		float denomY=Math.abs(ydef[y-1]-ydef[0])*EARTH_RADIUS;
		
		if(sVaMa.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(mcdata[l][k][j][0]!=undef) msdata[l][k][j][0]=mcdata[l][k][j][0]/denomY;
				if(vcdata[l][k][j][0]!=undef) vsdata[l][k][j][0]=vcdata[l][k][j][0]/denomZ;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(mcdata[k][j][0][l]!=undef) msdata[k][j][0][l]=mcdata[k][j][0][l]/denomY;
				if(vcdata[k][j][0][l]!=undef) vsdata[k][j][0][l]=vcdata[k][j][0][l]/denomZ;
			}
		}
		
		return EP;
	}
	
	/**
     * Calculate small-amplitude Eliassen-Palm (EP) vector in isentropic coordinate.
     *
     * @param	sVaUa	eddy horizontal absolute angular momentum flux = [(sv)'U'] (kg K^-1 s^-2)
     * @param	paMa	eddy horizontal heat flux = [p'M¦Ë'] (Pa m^2 s^-2)
     * @param	rZToY	ratio of Z to Y: Z/Y in the plot used to scale the vector components (e.g., 0.8)
     *
     * @return	EP		EP vector, [0] meridional (kg m K^-1 s^-2) and [1] vertical components (Pa m);
     * 					[2] and [3] are the same as [0] and [1] except for visual scaling according
     * 					to Edmon et al. (1980, JAS).
     */
	public Variable[] cIsentropicEPVectorByU(Variable sVaUa,Variable paMa,float rZToY){
		assignSubDomainParams(sVaUa);
		checkDimensions(sVaUa,paMa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable[] EP=new Variable[4];
		EP[0]=new Variable("EPyU",sVaUa); EP[0].setValue(undef);
		EP[1]=new Variable("EPzU",sVaUa); EP[1].setValue(undef);
		EP[0].setCommentAndUnit("meridional component of Eliassen-Palm vector EPy (kg m K^-1 s^-2)");
		EP[1].setCommentAndUnit("vertical component of Eliassen-Palm vector EPz (Pa m)");
		
		float[][][][] mcdata=EP[0].getData();
		float[][][][] vcdata=EP[1].getData();
		float[][][][] sadata=sVaUa.getData();
		float[][][][] pmdata= paMa.getData();
		
		if(sVaUa.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(sadata[l][k][j][0]!=undef) mcdata[l][k][j][0]=-sadata[l][k][j][0]*rs[j];
				if(pmdata[l][k][j][0]!=undef) vcdata[l][k][j][0]= pmdata[l][k][j][0]/SpatialModel.GRAVITY_ACCERLERATION;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(sadata[k][j][0][l]!=undef) mcdata[k][j][0][l]=-sadata[k][j][0][l]*rs[j];
				if(pmdata[k][j][0][l]!=undef) vcdata[k][j][0][l]= pmdata[k][j][0][l]/SpatialModel.GRAVITY_ACCERLERATION;
			}
		}
		
		EP[2]=EP[0].copy(); EP[2].setName("EPySU");
		EP[3]=EP[1].copy(); EP[3].setName("EPzSU");
		EP[2].setCommentAndUnit("scaled EPy by multiplying actual length representing 1m (kg m K^-1 s^-2)");
		EP[3].setCommentAndUnit("scaled EPz by multiplying actual length representing 1K (Pa m)");
		
		float[][][][] msdata=EP[2].getData();
		float[][][][] vsdata=EP[3].getData();
		
		float denomZ=Math.abs(zdef[z-1]-zdef[0])/rZToY;
		float denomY=Math.abs(ydef[y-1]-ydef[0])*EARTH_RADIUS;
		
		if(sVaUa.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(mcdata[l][k][j][0]!=undef) msdata[l][k][j][0]=mcdata[l][k][j][0]/denomY;
				if(vcdata[l][k][j][0]!=undef) vsdata[l][k][j][0]=vcdata[l][k][j][0]/denomZ;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(mcdata[k][j][0][l]!=undef) msdata[k][j][0][l]=mcdata[k][j][0][l]/denomY;
				if(vcdata[k][j][0][l]!=undef) vsdata[k][j][0][l]=vcdata[k][j][0][l]/denomZ;
			}
		}
		
		return EP;
	}
	
	/**
     * Calculate finite-amplitude Eliassen-Palm (EP) vector in isentropic coordinate.
     *
     * @param	sVmaMma	eddy absolute angular momentum flux = [sv*M*] (kg m K^-1 s^-2)
     * @param	paMa	eddy horizontal (adiabatic) heat flux = [p'M¦Ë'] (Pa m^2 s^-2)
     * @param	sWmaMma	eddy vertical (diabatic) absolute angular momentum flux = [sw*M*] (kg s^-2)
     * @param	rZToY	ratio of Z to Y: Z/Y in the plot used to scale the vector components (e.g., 0.8)
     *
     * @return	EP		EP vector, [0] meridional (kg m K^-1 s^-2) and [1] vertical components (Pa m);
     * 					[2] and [3] are the same as [0] and [1] except for visual scaling according
     * 					to Edmon et al. (1980, JAS).
     */
	public Variable[] cIsentropicFAEPVectorByM(Variable sVmaMma,Variable paMa,Variable sWmaMma,float rZToY){
		assignSubDomainParams(sVmaMma);
		checkDimensions(sVmaMma,paMa,sWmaMma);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable[] EP=new Variable[6];
		EP[0]=new Variable("EPmyM" ,sVmaMma); EP[0].setValue(undef);
		EP[1]=new Variable("EPmz1M",sVmaMma); EP[1].setValue(undef);
		EP[2]=new Variable("EPmz2M",sVmaMma); EP[2].setValue(undef);
		EP[0].setCommentAndUnit("meridional component of EP vector EPy -[sv*M*] (kg m K^-1 s^-2)");
		EP[1].setCommentAndUnit("adiabatic vertical component of EP vector EPz  [p'M¦Ë']/g (kg s^-2)");
		EP[2].setCommentAndUnit(" diabatic vertical component of EP vector EPz -[sw*M*]   (kg s^-2)");
		
		float[][][][] mcdata =EP[0].getData();
		float[][][][] vc1data=EP[1].getData();
		float[][][][] vc2data=EP[2].getData();
		float[][][][] svudata=sVmaMma.getData();
		float[][][][]  pmdata=   paMa.getData();
		float[][][][] swudata=sWmaMma.getData();
		
		if(sVmaMma.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(svudata[l][k][j][0]!=undef)  mcdata[l][k][j][0]=-svudata[l][k][j][0];
				if( pmdata[l][k][j][0]!=undef) vc1data[l][k][j][0]=  pmdata[l][k][j][0]/SpatialModel.GRAVITY_ACCERLERATION;
				if(swudata[l][k][j][0]!=undef) vc2data[l][k][j][0]=-swudata[l][k][j][0];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(svudata[k][j][0][l]!=undef)  mcdata[k][j][0][l]=-svudata[k][j][0][l];
				if( pmdata[k][j][0][l]!=undef) vc1data[k][j][0][l]=  pmdata[k][j][0][l]/SpatialModel.GRAVITY_ACCERLERATION;
				if(swudata[k][j][0][l]!=undef) vc2data[k][j][0][l]=-swudata[k][j][0][l];
			}
		}
		
		EP[3]=EP[0].copy(); EP[3].setName("EPmySM");
		EP[4]=EP[1].copy(); EP[4].setName("EPmz1SM");
		EP[5]=EP[2].copy(); EP[5].setName("EPmz2SM");
		EP[3].setCommentAndUnit("scaled EPy by multiplying actual length representing 1m (kg m K^-1 s^-2)");
		EP[4].setCommentAndUnit("scaled adiabatic EPz by multiplying actual length representing 1K (Pa m)");
		EP[5].setCommentAndUnit("scaled  diabatic EPz by multiplying actual length representing 1K (Pa m)");
		
		float[][][][] msdata =EP[3].getData();
		float[][][][] vs1data=EP[4].getData();
		float[][][][] vs2data=EP[5].getData();
		
		float denomZ=Math.abs(zdef[z-1]-zdef[0])/rZToY;
		float denomY=Math.abs(ydef[y-1]-ydef[0])*EARTH_RADIUS;
		
		if(sVmaMma.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if( mcdata[l][k][j][0]!=undef)  msdata[l][k][j][0]= mcdata[l][k][j][0]/denomY;
				if(vc1data[l][k][j][0]!=undef) vs1data[l][k][j][0]=vc1data[l][k][j][0]/denomZ;
				if(vc2data[l][k][j][0]!=undef) vs2data[l][k][j][0]=vc2data[l][k][j][0]/denomZ;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if( mcdata[k][j][0][l]!=undef)  msdata[k][j][0][l]= mcdata[k][j][0][l]/denomY;
				if(vc1data[k][j][0][l]!=undef) vs1data[k][j][0][l]=vc1data[k][j][0][l]/denomZ;
				if(vc2data[k][j][0][l]!=undef) vs2data[k][j][0][l]=vc2data[k][j][0][l]/denomZ;
			}
		}
		
		return EP;
	}
	
	/**
     * Calculate finite-amplitude Eliassen-Palm (EP) vector in isentropic coordinate.
     *
     * @param	sVmaUma	eddy absolute angular momentum flux = [sv*U*] (kg K^-1 s^-2)
     * @param	paMa	eddy horizontal (adiabatic) heat flux = [p'M¦Ë'] (Pa m^2 s^-2)
     * @param	sWmaUma	eddy vertical (diabatic) absolute angular momentum flux = [sw*U*] (kg s^-2)
     * @param	rZToY	ratio of Z to Y: Z/Y in the plot used to scale the vector components (e.g., 0.8)
     *
     * @return	EP		EP vector, [0] meridional (kg m K^-1 s^-2) and [1] vertical components (Pa m);
     * 					[2] and [3] are the same as [0] and [1] except for visual scaling according
     * 					to Edmon et al. (1980, JAS).
     */
	public Variable[] cIsentropicFAEPVectorByU(Variable sVmaUma,Variable paMa,Variable sWmaUma,float rZToY){
		assignSubDomainParams(sVmaUma);
		checkDimensions(sVmaUma,paMa,sWmaUma);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable[] EP=new Variable[6];
		EP[0]=new Variable("EPmyU" ,sVmaUma); EP[0].setValue(undef);
		EP[1]=new Variable("EPmz1U",sVmaUma); EP[1].setValue(undef);
		EP[2]=new Variable("EPmz2U",sVmaUma); EP[2].setValue(undef);
		EP[0].setCommentAndUnit("meridional component of EP vector EPy -[sv*U*]r (kg m K^-1 s^-2)");
		EP[1].setCommentAndUnit("adiabatic vertical component of EP vector EPz  [p'M¦Ë']/g (kg s^-2)");
		EP[2].setCommentAndUnit(" diabatic vertical component of EP vector EPz -[sw*U*]r  (kg s^-2)");
		
		float[][][][] mcdata =EP[0].getData();
		float[][][][] vc1data=EP[1].getData();
		float[][][][] vc2data=EP[2].getData();
		float[][][][] svudata=sVmaUma.getData();
		float[][][][]  pmdata=   paMa.getData();
		float[][][][] swudata=sWmaUma.getData();
		
		if(sVmaUma.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(svudata[l][k][j][0]!=undef)  mcdata[l][k][j][0]=-svudata[l][k][j][0]*rs[j];
				if( pmdata[l][k][j][0]!=undef) vc1data[l][k][j][0]=  pmdata[l][k][j][0]/SpatialModel.GRAVITY_ACCERLERATION;
				if(swudata[l][k][j][0]!=undef) vc2data[l][k][j][0]=-swudata[l][k][j][0]*rs[j];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(svudata[k][j][0][l]!=undef)  mcdata[k][j][0][l]=-svudata[k][j][0][l]*rs[j];
				if( pmdata[k][j][0][l]!=undef) vc1data[k][j][0][l]=  pmdata[k][j][0][l]/SpatialModel.GRAVITY_ACCERLERATION;
				if(swudata[k][j][0][l]!=undef) vc2data[k][j][0][l]=-swudata[k][j][0][l]*rs[j];
			}
		}
		
		EP[3]=EP[0].copy(); EP[3].setName("EPmySU");
		EP[4]=EP[1].copy(); EP[4].setName("EPmz1SU");
		EP[5]=EP[2].copy(); EP[5].setName("EPmz2SU");
		EP[3].setCommentAndUnit("scaled EPy by multiplying actual length representing 1m (kg m K^-1 s^-2)");
		EP[4].setCommentAndUnit("scaled adiabatic EPz by multiplying actual length representing 1K (Pa m)");
		EP[5].setCommentAndUnit("scaled  diabatic EPz by multiplying actual length representing 1K (Pa m)");
		
		float[][][][] msdata =EP[3].getData();
		float[][][][] vs1data=EP[4].getData();
		float[][][][] vs2data=EP[5].getData();
		
		float denomZ=Math.abs(zdef[z-1]-zdef[0])/rZToY;
		float denomY=Math.abs(ydef[y-1]-ydef[0])*EARTH_RADIUS;
		
		if(sVmaUma.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if( mcdata[l][k][j][0]!=undef)  msdata[l][k][j][0]= mcdata[l][k][j][0]/denomY;
				if(vc1data[l][k][j][0]!=undef) vs1data[l][k][j][0]=vc1data[l][k][j][0]/denomZ;
				if(vc2data[l][k][j][0]!=undef) vs2data[l][k][j][0]=vc2data[l][k][j][0]/denomZ;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if( mcdata[k][j][0][l]!=undef)  msdata[k][j][0][l]= mcdata[k][j][0][l]/denomY;
				if(vc1data[k][j][0][l]!=undef) vs1data[k][j][0][l]=vc1data[k][j][0][l]/denomZ;
				if(vc2data[k][j][0][l]!=undef) vs2data[k][j][0][l]=vc2data[k][j][0][l]/denomZ;
			}
		}
		
		return EP;
	}
	
	
	/**
     * implement the methods in EllipticalInterface
     */
	public Variable cAPrime(Variable sm){
		assignSubDomainParams(sm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable A=new Variable("Ap",sm);
		A.setCommentAndUnit("coefficient A' of elliptic equation (m^4 s^2 kg^-1)");
		A.setValue(undef);
		
		float[][][][] Adata= A.getData();
		float[][][][] Sdata=sm.getData();
		
		if(sm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			if(Sdata[l][k][j][0]!=undef)
			Adata[l][k][j][0]=(float)((Sdata[l][k][j-1][0]+Sdata[l][k][j][0])/2.0/sin((ydef[j-1]+ydef[j])/2.0));
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			if(Sdata[k][j][0][l]!=undef)
			Adata[k][j][0][l]=(float)((Sdata[k][j-1][0][l]+Sdata[k][j][0][l])/2.0/sin((ydef[j-1]+ydef[j])/2.0));
		}
		
		return A;
	}
	
	public Variable cBPrime(Variable thm){
		assignSubDomainParams(thm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable B=new Variable("Bp",thm);
		B.setCommentAndUnit("coefficient B' of elliptic equation (m^2 kg^-1)");
		B.setValue(undef);
		
		float[][][][] Bdata= B.getData();
		float[][][][] tdata=thm.getData();
		
		if(thm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				if(tdata[l][k][j+1][0]!=undef&&tdata[l][k][j-1][0]!=undef)
				Bdata[l][k][j][0]=(tdata[l][k][j+1][0]-tdata[l][k][j-1][0])/(2f*dy)/bsin[j]*
				Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
				
				if(tdata[l][k][y-1][0]!=undef&&tdata[l][k][y-2][0]!=undef)
				Bdata[l][k][y-1][0]=2f*(tdata[l][k][y-1][0]-tdata[l][k][y-2][0])/dy/(float)sin((ydef[y-2]+ydef[y-1])/2.0)*
				Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k]-Bdata[l][k][y-2][0];
				
				if(tdata[l][k][0][0]!=undef&&tdata[l][k][1][0]!=undef)
				Bdata[l][k][0][0]=2f*(tdata[l][k][1][0]-tdata[l][k][0][0])/dy/(float)sin(ydef[1]/2)*
				Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k]-Bdata[l][k][1][0];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				if(tdata[k][j+1][0][l]!=undef&&tdata[k][j-1][0][l]!=undef)
				Bdata[k][j][0][l]=(tdata[k][j+1][0][l]-tdata[k][j-1][0][l])/bsin[j]/(2f*dy)*
				Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k];
				
				if(tdata[k][y-1][0][l]!=undef&&tdata[k][y-2][0][l]!=undef)
				Bdata[k][y-1][0][l]=2f*(tdata[k][y-1][0][l]-tdata[k][y-2][0][l])/(float)sin((ydef[y-2]+ydef[y-1])/2.0)/dy*
				Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k]-Bdata[k][y-2][0][l];
				
				if(tdata[k][0][0][l]!=undef&&tdata[k][1][0][l]!=undef)
				Bdata[k][0][0][l]=2f*(tdata[k][1][0][l]-tdata[k][0][0][l])/dy/(float)sin(ydef[1]/2)*
				Rd*(float)pow(zdef[k]/100000.0,kp)/zdef[k]-Bdata[k][1][0][l];
			}
		}
		
		return B;
	}
	
	public Variable cCPrime(Variable gm){
		assignSubDomainParams(gm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable C=new Variable("Cp",gm);
		C.setCommentAndUnit("coefficient C' of elliptic equation (s^-2)");
		C.setValue(undef);
		
		float[]         buf=new float[y];
		float[][][][] Cdata= C.getData();
		float[][][][] gdata=gm.getData();
		
		if(gm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++){
				for(int j=1;j<y;j++)
				if(gdata[l][k-1][j][0]!=undef&&gdata[l][k][j][0]!=undef)
				buf[j]=(gdata[l][k-1][j][0]+gdata[l][k][j][0])/2f;
				
				for(int j=1,J=y-1;j<J;j++)
				if(buf[j]!=undef&&buf[j+1]!=undef&&buf[j-1]!=undef)
				Cdata[l][k][j][0]=bcos[j]*2f*buf[j]*(buf[j+1]-buf[j-1])/(2f*dy*rs[j]*rs[j]*rs[j]*bsin[j]);
				
				if(buf[y-1]!=undef&&buf[y-2]!=undef)
				Cdata[l][k][y-1][0]=bcos[y-1]*(buf[y-1]+buf[y-2])*(buf[y-1]-buf[y-2])/(dy*rs[y-1]*rs[y-1]*rs[y-1]*bsin[y-1]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++){
				for(int j=1;j<y;j++)
				if(gdata[k-1][j][0][l]!=undef&&gdata[k][j][0][l]!=undef)
				buf[j]=(gdata[k-1][j][0][l]+gdata[k][j][0][l])/2f;
				
				for(int j=1,J=y-1;j<J;j++)
				if(buf[j]!=undef&&buf[j+1]!=undef&&buf[j-1]!=undef)
				Cdata[k][j][0][l]=bcos[j]*2f*buf[j]*(buf[j+1]-buf[j-1])/(2*dy*rs[j]*rs[j]*rs[j]*bsin[j]);
				
				if(buf[y-1]!=undef&&buf[y-2]!=undef)
				Cdata[k][y-1][0][l]=bcos[y-1]*(buf[y-1]+buf[y-2])*(buf[y-1]-buf[y-2])/(dy*rs[y-1]*rs[y-1]*rs[y-1]*bsin[y-1]);
			}
		}
		
		return C;
	}
	
	public Variable cA(Variable sm){
		assignSubDomainParams(sm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable A=new Variable("Aa",sm);
		A.setCommentAndUnit("coefficient A of elliptical equation (m^4 s^2 kg^-1)");
		A.setValue(undef);
		
		float[][][][] Adata= A.getData();
		float[][][][] Sdata=sm.getData();
		
		if(sm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			if(Sdata[l][k][j][0]!=undef) Adata[l][k][j][0]=Sdata[l][k][j][0]/bsin[j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			if(Sdata[k][j][0][l]!=undef) Adata[k][j][0][l]=Sdata[k][j][0][l]/bsin[j];
		}
		
		return A;
	}
	
	public Variable cB(Variable thm){
		Variable B=cBPrime(thm);
		
		B.setName("Bb");
		B.setCommentAndUnit("coefficient B of elliptic equation (m^2 kg^-1)");
		
		return B;
	}
	
	public Variable cC(Variable gm){
		assignSubDomainParams(gm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable C=new Variable("Cc",gm);
		C.setCommentAndUnit("coefficient C of elliptical equation (s^-2)");
		C.setValue(undef);
		
		float[][][][] Cdata= C.getData();
		float[][][][] gdata=gm.getData();
		
		if(gm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				if(gdata[l][k][j+1][0]!=undef&&gdata[l][k][j][0]!=undef&&gdata[l][k][j-1][0]!=undef)
				Cdata[l][k][j][0]=bcos[j]*2f*gdata[l][k][j][0]*(gdata[l][k][j+1][0]-gdata[l][k][j-1][0])/(2*dy*rs[j]*rs[j]*rs[j]*bsin[j]);
				
				if(gdata[l][k][y-1][0]!=undef&&gdata[l][k][y-2][0]!=undef)
				Cdata[l][k][y-1][0]=bcos[y-1]*(gdata[l][k][y-1][0]+gdata[l][k][y-2][0])*
				(gdata[l][k][y-1][0]-gdata[l][k][y-2][0])/(dy*rs[y-1]*rs[y-1]*rs[y-1]*bsin[y-1]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				if(gdata[k][j+1][0][l]!=undef&&gdata[k][j][0][l]!=undef&&gdata[k][j-1][0][l]!=undef)
				Cdata[k][j][0][l]=bcos[j]*2f*gdata[k][j][0][l]*(gdata[k][j+1][0][l]-gdata[k][j-1][0][l])/(2*dy*rs[j]*rs[j]*rs[j]*bsin[j]);
				
				if(gdata[k][y-1][0][l]!=undef&&gdata[k][y-2][0][l]!=undef)
				Cdata[k][y-1][0][l]=bcos[y-1]*(gdata[k][y-1][0][l]+gdata[k][y-2][0][l])*
				(gdata[k][y-1][0][l]-gdata[k][y-2][0][l])/(dy*rs[y-1]*rs[y-1]*rs[y-1]*bsin[y-1]);
			}
		}
		
		return C;
	}
	
	
	/**
     * Initialize the boundaries of streamfunction using vm and wm.
     *
     * @param	sf	streamfunction (Pa m s^-1)
     * @param	vm	tangential-mean radial wind (m s^-1)
     * @param	wm	tangential-mean vertical velocity (Pa s^-1)
     */
	public void initialSFBoundary(Variable sf,Variable vm,Variable wm){
		checkDimensions(sf,vm,wm);
		assignSubDomainParams(sf);
		
		if(x!=1) throw new IllegalArgumentException("x-count should be 1 (tangential averaged)");
		
		float[][][][]  vdata=vm.getData();
		float[][][][]  wdata=wm.getData();
		float[][][][] sfdata=sf.getData();
		
		if(sf.isTFirst()){
			for(int l=0;l<t;l++){
				float err1=0,err2=0,max1=0,max2=0,vsf=0,wsf=0,ratio=0;
				/**/
				// caculate lateral boundries using vm
				for(int k=1;k<z;k++){
					sfdata[l][k][y-1][0]=sfdata[l][k-1][y-1][0]+
					dz*(vdata[l][k][y-1][0]+vdata[l][k-1][y-1][0])/2*bsin[y-1];
					
					if(Math.abs(sfdata[l][k][y-1][0])>max1) max1=sfdata[l][k][y-1][0];
				}
				
				vsf=sfdata[l][z-1][y-1][0];
				/**/
				// caculate lateral boundries using wm
				for(int j=1;j<y;j++){
					sfdata[l][z-1][j][0]=sfdata[l][z-1][j-1][0]-
					dy*(wdata[l][z-1][j-1][0]+wdata[l][z-1][j][0])/2*(float)sin((ydef[j]+ydef[j-1])/2);
					
					if(Math.abs(sfdata[l][z-1][j][0])>max2) max2=sfdata[l][z-1][j][0];
				}
				
				wsf=sfdata[l][z-1][y-1][0];
				
				// calculate err
				ratio=max1/(max1+max2);
				
				err1=(vsf-wsf)*ratio/(z-1);
				err2=(vsf-wsf)*(1-ratio)/(y-1);
				
				// correct lateral boundries
				for(int k=1;k<z-1;k++) sfdata[l][k][y-1][0]-=err1*k;
				for(int j=1;j<y  ;j++) sfdata[l][z-1][j][0]+=err2*j;
			}
			
		}else{
			for(int l=0;l<t;l++){
				float err1=0,err2=0,max1=0,max2=0,vsf=0,wsf=0,ratio=0;
				/**/
				// caculate lateral boundries using vm
				for(int k=1;k<z;k++){
					sfdata[k][y-1][0][l]=sfdata[k-1][y-1][0][l]+
					dz*(vdata[k][y-1][0][l]+vdata[k-1][y-1][0][l])/2*bsin[y-1];
					
					if(Math.abs(sfdata[k][y-1][0][l])>max1) max1=sfdata[k][y-1][0][l];
				}
				
				vsf=sfdata[z-1][y-1][0][l];
				/**/
				// caculate lateral boundries using wm
				for(int j=1;j<y;j++){
					sfdata[z-1][j][0][l]=sfdata[z-1][j-1][0][l]-
					dy*(wdata[z-1][j-1][0][l]+wdata[z-1][j][0][l])/2*(float)sin((ydef[j]+ydef[j-1])/2);
					
					if(Math.abs(sfdata[z-1][j][0][l])>max2) max2=sfdata[z-1][j][0][l];
				}
				
				wsf=sfdata[z-1][y-1][0][l];
				
				// calculate err
				ratio=max1/(max1+max2);
				
				err1=(vsf-wsf)*ratio/(z-1);
				err2=(vsf-wsf)*(1-ratio)/(y-1);
				
				// correct lateral boundries
				for(int k=1;k<z-1;k++) sfdata[k][y-1][0][l]-=err1*k;
				for(int j=1;j<y  ;j++) sfdata[z-1][j][0][l]+=err2*j;
			}
		}
	}
	
	
	/*** helper methods ***/
	
	/**
     * Calculate radial-component of divergence in the YZ-plane.
     *
     * @param	vecY	y-component of a vector in the YZ-plane
     */
	private Variable cRadialDiv(Variable vecY){
		assignSubDomainParams(vecY);
		
		Variable divY=new Variable("divY",vecY); divY.setValue(undef);
		divY.setCommentAndUnit("y-comp divergence in YZ plane (vecY * m^-1)");
		
		float[][][][] dvdata=divY.getData();
		float[][][][] vYdata=vecY.getData();
		
		if(vecY.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				// core boundary (j==0)
				if(vYdata[l][k][0][i]!=undef&&vYdata[l][k][1][i]!=undef)
				dvdata[l][k][0][i]=(vYdata[l][k][1][i]*bsin[1]-vYdata[l][k][0][i]*bsin[0])
				/dy/((bsin[1]+bsin[0])/2f);
				
				// inner region
				for(int j=1,Y=y-1;j<Y;j++)
				if(vYdata[l][k][j+1][i]!=undef&&vYdata[l][k][j-1][i]!=undef)
				dvdata[l][k][j][i]=(vYdata[l][k][j+1][i]*bsin[j+1]-vYdata[l][k][j-1][i]*bsin[j-1])/(2f*dy)/bsin[j];
				
				// outer boundary (j==y-1)
				if(vYdata[l][k][y-1][i]!=undef&&vYdata[l][k][y-2][i]!=undef)
				dvdata[l][k][y-1][i]=(vYdata[l][k][y-1][i]*bsin[y-1]-vYdata[l][k][y-2][i]*bsin[y-2])
				/dy/((bsin[y-1]+bsin[y-2])/2f);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				// core boundary (j==0)
				if(vYdata[k][0][i][l]!=undef&&vYdata[k][1][i][l]!=undef)
				dvdata[k][0][i][l]=(vYdata[k][1][i][l]*bsin[1]-vYdata[k][0][i][l]*bsin[0])
				/dy/((bsin[1]+bsin[0])/2f);
				
				// inner region
				for(int j=1,Y=y-1;j<Y;j++)
				if(vYdata[k][j+1][i][l]!=undef&&vYdata[k][j-1][i][l]!=undef)
				dvdata[k][j][i][l]=(vYdata[k][j+1][i][l]*bsin[j+1]-vYdata[k][j-1][i][l]*bsin[j-1])/(2f*dy)/bsin[j];
				
				// outer boundary (j==y-1)
				if(vYdata[k][y-1][i][l]!=undef&&vYdata[k][y-2][i][l]!=undef)
				dvdata[k][y-1][i][l]=(vYdata[k][y-1][i][l]*bsin[y-1]-vYdata[k][y-2][i][l]*bsin[y-2])
				/dy/((bsin[y-1]+bsin[y-2])/2f);
			}
		}
		
		return divY;
	}
	
	/**
     * Calculate vertical-component of divergence in the YZ-plane.
     *
     * @param	vecZ	z-component of a vector in the YZ-plane
     */
	private Variable cVerticalDiv(Variable vecZ){
		assignSubDomainParams(vecZ);
		
		Variable div=new Variable("divZ",vecZ); div.setValue(undef);
		div.setCommentAndUnit("z-comp divergence in YZ plane (vecZ * Pa^-1)");
		
		float[][][][] dvdata= div.getData();
		float[][][][] vZdata=vecZ.getData();
		
		if(vecZ.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(z>1){
				// lower boundary (k==0)
				if(vZdata[l][0][j][i]!=undef&&vZdata[l][1][j][i]!=undef)
				dvdata[l][0][j][i]=(vZdata[l][1][j][i]-vZdata[l][0][j][i])/dz;
				
				// inner region
				for(int k=1,K=z-1;k<K;k++)
				if(vZdata[l][k+1][j][i]!=undef&&vZdata[l][k-1][j][i]!=undef)
				dvdata[l][k][j][i]=(vZdata[l][k+1][j][i]-vZdata[l][k-1][j][i])/(2f*dz);
				
				// upper boundary (k==z-1)
				if(vZdata[l][z-1][j][i]!=undef&&vZdata[l][z-2][j][i]!=undef)
				dvdata[l][z-1][j][i]=(vZdata[l][z-1][j][i]-vZdata[l][z-2][j][i])/dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(z>1){
				// lower boundary (k==0)
				if(vZdata[0][j][i][l]!=undef&&vZdata[1][j][i][l]!=undef)
				dvdata[0][j][i][l]=(vZdata[1][j][i][l]-vZdata[0][j][i][l])/dz;
				
				// inner region
				for(int k=1,K=z-1;k<K;k++)
				if(vZdata[k+1][j][i][l]!=undef&&vZdata[k-1][j][i][l]!=undef)
				dvdata[k][j][i][l]=(vZdata[k+1][j][i][l]-vZdata[k-1][j][i][l])/(2f*dz);
				
				// upper boundary (k==z-1)
				if(vZdata[z-1][j][i][l]!=undef&&vZdata[z-2][j][i][l]!=undef)
				dvdata[z-1][j][i][l]=(vZdata[z-1][j][i][l]-vZdata[z-2][j][i][l])/dz;
			}
		}
		
		return div;
	}
	
	
	/**
     * Mark the vertical indices of surface topography.
     *
     * @param	sfp		surface pressure
     */
	private void markTopographyIdx(Variable sfp){
		assignSubDomainParams(sfp);
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		float[][][][] sfpdata=sfp.getData();
		
		float tmp=zdef[zdef.length-1]-zdef[0];
		
		if(tmp==0) throw new IllegalArgumentException("zdef is only one level or not monototic");
		
		boolean incre=tmp>0?true:false;
		
		if(sfp.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				sfcIdx[l][j][i]=incre?
				(short)ArrayUtil.getIdxIncre(zdef,sfpdata[l][0][j][i]):
				(short)ArrayUtil.getIdxDecre(zdef,sfpdata[l][0][j][i]);
				if(sfcIdx[l][j][i]<0) sfcIdx[l][j][i]=0;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				sfcIdx[l][j][i]=incre?
				(short)ArrayUtil.getIdxIncre(zdef,sfpdata[0][j][i][l]):
				(short)ArrayUtil.getIdxDecre(zdef,sfpdata[0][j][i][l]);
				if(sfcIdx[l][j][i]<0) sfcIdx[l][j][i]=0;
			}
		}
	}
	
	/**
     * Calculate the translating velocity of the coordinates.
     */
	private void cWhole(){
		CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
		
		int t=csm.getTCount(),x=csm.getXCount();
		
		float[] cu=csm.getUWhole();
		float[] cv=csm.getVWhole();
		
		cs=new float[t];	// central speed of the coordinates
		cd=new float[t];	// central direction of the speed
		
		wtan=new float[t][x];
		wnor=new float[t][x];
		
		for(int l=0;l<t;l++){
			cs[l]=(float)Math.hypot(cu[l],cv[l]);
			cd[l]=(float)(Math.atan2(cv[l],cu[l])-PI/2); // angle start from north, counter-clockwise
			
			for(int i=0;i<x;i++){
				float angle=cd[l]-xdef[i];
				wtan[l][i]=cs[l]*(float)sin(angle);
				wnor[l][i]=cs[l]*(float)cos(angle);
			}
		}
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			
		}catch(Exception ex){ ex.printStackTrace();}
	}*/
}
