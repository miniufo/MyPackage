/**
 * @(#)EliassenModelInSC.java	1.0 2015.11.11
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.diagnosticModel;

import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInSphericalCoordinate;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static miniufo.diagnosis.SpatialModel.REarth;
import static miniufo.diagnosis.SpatialModel.omegaEarth;
import static miniufo.geophysics.atmos.ThermoDynamics.Cp;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;
import static miniufo.geophysics.atmos.ThermoDynamics.kappa;


/**
 * Eliassen model in Spherical coordinate (Hadley model)
 *
 * @version 1.0, 2015.11.11
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EliassenModelInSC extends EquationInSphericalCoordinate
implements EllipticEquationInterface{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spatial model in spherical coordinate
     */
	public EliassenModelInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isPeriodicX()  ) throw new IllegalArgumentException("spatial model should be zonally periodic");
		if(!ssm.isLinearModel()) throw new IllegalArgumentException("EliassenModelInSC requires a linear model");
	}
	
	
	/**
     * Calculate absolute angular momentum about the rotating axis of the earth.
     *
     * @param	u	zonal wind (m s^-1)
     *
     * @return	ga	absolute angular momentum (m^2 s^-1)
     */
	public Variable cAbsoluteAngularMomentum(Variable u){
		assignSubDomainParams(u);
		
		Variable ga=new Variable("ga",u);
		ga.setCommentAndUnit("absolute angular momentum (m^2 s^-1)");
		ga.setValue(undef);
		
		float[][][][] Mdata=ga.getData();
		float[][][][] udata= u.getData();
		
		if(u.isTFirst()){
			for(int j=0;j<y;j++){
				float r=rcos[ystart-1+j];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++)
				Mdata[l][k][j][i]=(udata[l][k][j][i]+omegaEarth*r)*r;
			}
			
		}else{
			for(int j=0;j<y;j++){
				float r=rcos[ystart-1+j];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++)
				Mdata[k][j][i][l]=(udata[k][j][i][l]+omegaEarth*r)*r;
			}
		}
		
		return ga;
	}
	
	
	/**
     * Calculate force due to eddy heat horizontal flux convergence (HFC)
     *
     * @param	tava	eddy heat horizontal flux = [ï¿½ï¿½'v'] (K m s^-1)
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
			fdata[l][k][j][0]=(
				(tvdata[l][k][j+1][0]*lcos[ystart-1+j+1]-tvdata[l][k][j][0]*lcos[ystart-1+j])
				/((lcos[ystart-1+j+1]+lcos[ystart-1+j])/2f)-
				(tvdata[l][k][j][0]*lcos[ystart-1+j]-tvdata[l][k][j-1][0]*lcos[ystart-1+j-1])
				/((lcos[ystart-1+j]+lcos[ystart-1+j-1])/2f)
			)/(dy*dy)*Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(tvdata[k][j+1][0][l]!=undef&&tvdata[k][j][0][l]!=undef&&tvdata[k][j-1][0][l]!=undef)
			fdata[k][j][0][l]=(
				(tvdata[k][j+1][0][l]*lcos[ystart-1+j+1]-tvdata[k][j][0][l]*lcos[ystart-1+j])
				/((lcos[ystart-1+j+1]+lcos[ystart-1+j])/2f)-
				(tvdata[k][j][0][l]*lcos[ystart-1+j]-tvdata[k][j-1][0][l]*lcos[ystart-1+j-1])
				/((lcos[ystart-1+j]+lcos[ystart-1+j-1])/2f)
			)/(dy*dy)*Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
		}
		
		return f;
	}
	
	/**
     * Calculate force due to eddy heat vertical flux convergence (VFC)
     *
     * @param	tawa	eddy heat vertical flux = [ï¿½ï¿½'w'] (K Pa s^-1)
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
			fdata[l][k][j][0]=(
				(twdata[l][k+1][j+1][0]-twdata[l][k-1][j+1][0])/(2f*dz)-
				(twdata[l][k+1][j-1][0]-twdata[l][k-1][j-1][0])/(2f*dz)
			)/(2f*dy)*Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(twdata[k+1][j+1][0][l]!=undef&&twdata[k-1][j+1][0][l]!=undef
			 &&twdata[k+1][j-1][0][l]!=undef&&twdata[k-1][j-1][0][l]!=undef)
			fdata[k][j][0][l]=(
				(twdata[k+1][j+1][0][l]-twdata[k-1][j+1][0][l])/(2f*dz)-
				(twdata[k+1][j-1][0][l]-twdata[k-1][j-1][0][l])/(2f*dz)
			)/(2f*dy)*Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
		}
		
		return f;
	}
	
	/**
     * Calculate forcing factor due to heating functions,
     * i.e., FF = heat * R * ï¿½ï¿½ / p.
     *
     * @param	heat	heating function e.g., eddy heat HFC or VFC (K s^-1)
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
			hdata[l][k][j][0]*Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(hdata[k][j][0][l]!=undef) fdata[k][j][0][l]=
			hdata[k][j][0][l]*Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
		}
		
		return ff;
	}
	
	/**
     * Calculate diabatic heating rate d¦È/dt = Qm / (Cp * ¦Ð).
     *
     * @param	Qm	zonally averaged heating rate (W kg^-1 or m^2 s^-3)
     *
     * @return	f	diabatic heating rate (K s^-1)
     */
    public Variable cDiabaticHeating(Variable Qm){
		assignSubDomainParams(Qm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("dht",Qm);
		f.setCommentAndUnit("diabatic heating rate (K s^-1)");
		f.setValue(undef);
		
		float[][][][] fdata= f.getData();
		float[][][][] Qdata=Qm.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) if(Qdata[l][k][j][0]!=undef)
			fdata[l][k][j][0]=-Qdata[l][k][j][0]/Cp/(float)pow(zdef[zstart-1+k]/100000.0,kappa);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) if(Qdata[k][j][0][l]!=undef)
			fdata[k][j][0][l]=-Qdata[k][j][0][l]/Cp/(float)pow(zdef[zstart-1+k]/100000.0,kappa);
		}
		
		return f;
	}
    
	/**
     * Calculate eddy heat horizontal flux convergence (HFC).
     *
     * @param	tava	eddy heat horizontal flux = [ï¿½ï¿½'v'] (K m s^-1)
     *
     * @return	f		eddy heat HFC (K s^-1)
     */
	public Variable cEddyHeatHFC(Variable tava){
		assignSubDomainParams(tava);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("htHFC",tava);
		f.setCommentAndUnit("eddy heat horizontal flux convergence (K s^-1)");
		f.setValue(undef);
		
		float[][][][] fdata=   f.getData();
		float[][][][] vdata=tava.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used in inversion
				if(vdata[l][k][0][0]!=undef&&vdata[l][k][1][0]!=undef)
				fdata[l][k][0][0]=-(vdata[l][k][1][0]*lcos[ystart-1+1]-vdata[l][k][0][0]*lcos[ystart-1])
				/dy/((lcos[ystart-1+1]+lcos[ystart-1])/2f);
				
				// inner region
				for(int j=1;j<y-1;j++)
				if(vdata[l][k][j+1][0]!=undef&&vdata[l][k][j-1][0]!=undef) fdata[l][k][j][0]=-
				(vdata[l][k][j+1][0]*lcos[ystart-1+j+1]-vdata[l][k][j-1][0]*lcos[ystart-1+j-1])/(2f*dy)/lcos[ystart-1+j];
				
				// north boundary (j==y-1), not used in inversion
				if(vdata[l][k][y-1][0]!=undef&&vdata[l][k][y-2][0]!=undef)
				fdata[l][k][y-1][0]=-(vdata[l][k][y-1][0]*lcos[ystart-1+y-1]-vdata[l][k][y-2][0]*lcos[ystart-1+y-2])
				/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used in inversion
				if(vdata[k][0][0][l]!=undef&&vdata[k][1][0][l]!=undef)
				fdata[k][0][0][l]=-(vdata[k][1][0][l]*lcos[ystart-1+1]-vdata[k][0][0][l]*lcos[ystart-1])
				/dy/((lcos[ystart-1+1]+lcos[ystart-1])/2f);
				
				// inner region
				for(int j=1;j<y-1;j++)
				if(vdata[k][j+1][0][l]!=undef&&vdata[k][j-1][0][l]!=undef) fdata[k][j][0][l]=-
				(vdata[k][j+1][0][l]*lcos[ystart-1+j+1]-vdata[k][j-1][0][l]*lcos[ystart-1+j-1])/(2f*dy)/lcos[ystart-1+j];
				
				// north boundary (j==y-1), not used in inversion
				if(vdata[k][y-1][0][l]!=undef&&vdata[k][y-2][0][l]!=undef)
				fdata[k][y-1][0][l]=-(vdata[k][y-1][0][l]*lcos[ystart-1+y-1]-vdata[k][y-2][0][l]*lcos[ystart-1+y-2])
				/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f);
			}
		}
		
		return f;
	}
	
	/**
     * Calculate eddy heat vertical flux convergence (VFC).
     *
     * @param	tawa	eddy heat vertical flux = [ï¿½ï¿½'ï¿½ï¿½'] (K Pa s^-1)
     *
     * @return	f		eddy heat VHC (K s^-1)
     */
	public Variable cEddyHeatVFC(Variable tawa){
		assignSubDomainParams(tawa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("htVFC",tawa);
		f.setCommentAndUnit("eddy heat vertical flux convergence (K s^-1)");
		f.setValue(undef);
		
		float[][][][] fdata=   f.getData();
		float[][][][] wdata=tawa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(wdata[l][0][j][0]!=undef&&wdata[l][1][j][0]!=undef)
				fdata[l][0][j][0]=-(wdata[l][1][j][0]-wdata[l][0][j][0])/dz;
				
				// inner region
				for(int k=1,K=z-1;k<K;k++)
				if(wdata[l][k+1][j][0]!=undef&&wdata[l][k-1][j][0]!=undef)
				fdata[l][k][j][0]=-(wdata[l][k+1][j][0]-wdata[l][k-1][j][0])/(2f*dz);
				
				// upper boundary (k==z-1), not used in inversion
				if(wdata[l][z-1][j][0]!=undef&&wdata[l][z-2][j][0]!=undef)
				fdata[l][z-1][j][0]=-(wdata[l][z-1][j][0]-wdata[l][z-2][j][0])/dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(wdata[0][j][0][l]!=undef&&wdata[1][j][0][l]!=undef)
				fdata[0][j][0][l]=-(wdata[1][j][0][l]-wdata[0][j][0][l])/dz;
				
				// inner region
				for(int k=1,K=z-1;k<K;k++)
				if(wdata[k+1][j][0][l]!=undef&&wdata[k-1][j][0][l]!=undef)
				fdata[k][j][0][l]=-(wdata[k+1][j][0][l]-wdata[k-1][j][0][l])/(2f*dz);
				
				// upper boundary (k==z-1), not used in inversion
				if(wdata[z-1][j][0][l]!=undef&&wdata[z-2][j][0][l]!=undef)
				fdata[z-1][j][0][l]=-(wdata[z-1][j][0][l]-wdata[z-2][j][0][l])/dz;
			}
		}
		
		return f;
	}
	
	/**
     * Calculate eddy heat EP flux term.
     *
     * @param	tava	eddy heat horizontal flux = [ï¿½ï¿½'v'] (K m s^-1)
     * @param	tawa	eddy heat vertical flux = [ï¿½ï¿½'ï¿½ï¿½'] (K Pa s^-1)
     * @param	thm		zonal mean potential temperature (K)
     *
     * @return	f		eddy heat VHC (K s^-1)
     */
	public Variable cEddyHeatEP(Variable tava,Variable tawa,Variable thm){
		assignSubDomainParams(tawa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("htEPF",tawa);
		f.setCommentAndUnit("eddy heat EP flux (K s^-1)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] tvdata=tava.getData();
		float[][][][] twdata=tawa.getData();
		float[][][][] tmdata= thm.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1,J=y-1;j<J;j++){
				// lower boundary (k==0), not used in inversion
				
				// inner region
				for(int k=1,K=z-1;k<K;k++)
				if(tvdata[l][k][j][0]!=undef&&twdata[l][k][j][0]!=undef&&tmdata[l][k+1][j][0]!=undef&&
				tmdata[l][k-1][j][0]!=undef&&tmdata[l][k][j+1][0]!=undef&&tmdata[l][k][j-1][0]!=undef)
				fdata[l][k][j][0]=-twdata[l][k][j][0]-tvdata[l][k][j][0]*dz/dy*
				(tmdata[l][k][j+1][0]-tmdata[l][k][j-1][0])/(tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0]);
				
				// upper boundary (k==z-1), not used in inversion
			}
			
		}else{
		}
		
		return f;
	}
	
	
	/**
     * Calculate force due to eddy absolute angular momentum (AAM)
     * horizontal flux convergence (HFC).
     *
     * @param	gm		zonally averaged absolute angular momentum (m^2 s^-1)
     * @param	gava	eddy AAM horizontal flux = [g'v'] (m^3 s^-2)
     *
     * @return	f		force (m^2 kg^-1 s^-1)
     */
	public Variable cEddyAAMHFCForce(Variable gm,Variable gava){
		assignSubDomainParams(gm);
		checkDimensions(gm,gava);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("aamVFCFor",gm);
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
			fdata[l][k][j][0]=-2f*lsin[ystart-1+j]*(
				(gvdata[l][k+1][j+1][0]*lcos[ystart-1+j+1]-gvdata[l][k+1][j-1][0]*lcos[ystart-1+j-1])*
				gmdata[l][k+1][j][0]/lcos[ystart-1+j]/(2f*dy)-
				(gvdata[l][k-1][j+1][0]*lcos[ystart-1+j+1]-gvdata[l][k-1][j-1][0]*lcos[ystart-1+j-1])*
				gmdata[l][k-1][j][0]/lcos[ystart-1+j]/(2f*dy)
			)/(float)pow(rcos[ystart-1+j],3.0)/(2f*dz);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=1,J=y-1;j<J;j++)
			if(gmdata[k+1][j+1][0][l]!=undef&&gmdata[k][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&
			gvdata[k+1][j][0][l]!=undef&&gvdata[k][j][0][l]!=undef&&gvdata[k-1][j][0][l]!=undef)
			fdata[k][j][0][l]=-2f*lsin[ystart-1+j]*(
				(gvdata[k+1][j+1][0][l]*lcos[ystart-1+j+1]-gvdata[k+1][j-1][0][l]*lcos[ystart-1+j-1])*
				gmdata[k+1][j][0][l]/lcos[ystart-1+j]/(2f*dy)-
				(gvdata[k-1][j+1][0][l]*lcos[ystart-1+j+1]-gvdata[k-1][j-1][0][l]*lcos[ystart-1+j-1])*
				gmdata[k-1][j][0][l]/lcos[ystart-1+j]/(2f*dy)
			)/(float)pow(rcos[ystart-1+j],3.0)/(2f*dz);
		}
		
		return f;
	}
	
	/**
     * Calculate force due to eddy absolute angular momentum (AAM)
     * vertical flux convergence (HFC).
     *
     * @param	gm		zonally averaged absolute angular momentum (m^2 s^-1)
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
			for(int j=0;j<y;j++)
			if(gmdata[l][k+1][j][0]!=undef&&gmdata[l][k][j][0]!=undef&&gmdata[l][k-1][j][0]!=undef&&
			gwdata[l][k+1][j][0]!=undef&&gwdata[l][k][j][0]!=undef&&gwdata[l][k-1][j][0]!=undef)
			fdata[l][k][j][0]=-2f*lsin[ystart-1+j]*(
				(gmdata[l][k+1][j][0]+gmdata[l][k][j][0])/2f*(gwdata[l][k+1][j][0]-gwdata[l][k][j][0])-
				(gmdata[l][k-1][j][0]+gmdata[l][k][j][0])/2f*(gwdata[l][k][j][0]-gwdata[l][k-1][j][0])
			)/(float)pow(rcos[ystart-1+j],3.0)/(dz*dz);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1,K=z-1;k<K;k++)
			for(int j=0;j<y;j++)
			if(gmdata[k+1][j][0][l]!=undef&&gmdata[k][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&
			gwdata[k+1][j][0][l]!=undef&&gwdata[k][j][0][l]!=undef&&gwdata[k-1][j][0][l]!=undef)
			fdata[k][j][0][l]=-2f*lsin[ystart-1+j]*(
				(gmdata[k+1][j][0][l]+gmdata[k][j][0][l])/2f*(gwdata[k+1][j][0][l]-gwdata[k][j][0][l])-
				(gmdata[k-1][j][0][l]+gmdata[k][j][0][l])/2f*(gwdata[k][j][0][l]-gwdata[k-1][j][0][l])
			)/(float)pow(rcos[ystart-1+j],3.0)/(dz*dz);
		}
		
		return f;
	}
	
	/**
     * Calculate forcing factor due to momentum functions,
     * i.e., FF = mome * 2 * gm * sin(lat) / rcos ^ 3.
     *
     * @param	gm		zonal mean absolute angular momentum (m^2 s^-1)
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
			hdata[l][k][j][0]*2f*gdata[l][k][j][0]*lsin[ystart-1+j]/(float)pow(rcos[ystart-1+j],3.0);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(hdata[k][j][0][l]!=undef) fdata[k][j][0][l]=
			hdata[k][j][0][l]*2f*gdata[k][j][0][l]*lsin[ystart-1+j]/(float)pow(rcos[ystart-1+j],3.0);
		}
		
		return ff;
	}
	
	/**
     * Calculate zonal frictional torque.
     *
     * @param	fm	zonally averaged zonal friction (m s^-2)
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
			if(frdata[l][k][j][0]!=undef) ftdata[l][k][j][0]*=rcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(frdata[k][j][0][l]!=undef) ftdata[k][j][0][l]*=rcos[ystart-1+j];
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
		assignSubDomainParams(gava);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("aamHFC",gava);
		f.setCommentAndUnit("eddy absolute angular momentum horizontal flux convergence (m^2 s^-2)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] gvdata=gava.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used for inversion
				if(gvdata[l][k][0][0]!=undef) fdata[l][k][0][0]=
				-(gvdata[l][k][0][0]*lcos[ystart-1]-gvdata[l][k][1][0]*lcos[ystart-1+1])
				/dy/((lcos[ystart-1]+lcos[ystart-1+1])/2f);
				
				// inner region, centered difference
				for(int j=1;j<y-1;j++)
				if(gvdata[l][k][j+1][0]!=undef&&gvdata[l][k][j-1][0]!=undef) fdata[l][k][j][0]=
				-(gvdata[l][k][j+1][0]*lcos[ystart-1+j+1]-gvdata[l][k][j-1][0]*lcos[ystart-1+j-1])
				/(dy*2f)/(lcos[ystart-1+j]);
				
				// north boundary (j==y-1), not used for inversion
				if(gvdata[l][k][y-2][0]!=undef) fdata[l][k][y-1][0]=
				-(gvdata[l][k][y-1][0]*lcos[ystart-1+y-1]-gvdata[l][k][y-2][0]*lcos[ystart-1+y-2])
				/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used for inversion
				if(gvdata[k][0][0][l]!=undef) fdata[k][0][0][l]=
				-(gvdata[k][0][0][l]*lcos[ystart-1]-gvdata[k][1][0][l]*lcos[ystart-1+1])
				/dy/((lcos[ystart-1]+lcos[ystart-1+1])/2f);
				
				// inner region, centered difference
				for(int j=1;j<y-1;j++)
				if(gvdata[k][j+1][0][l]!=undef&&gvdata[k][j-1][0][l]!=undef) fdata[k][j][0][l]=
				-(gvdata[k][j+1][0][l]*lcos[ystart-1+j+1]-gvdata[k][j-1][0][l]*lcos[ystart-1+j-1])
				/(dy*2f)/(lcos[ystart-1+j]);
				
				// north boundary (j==y-1), not used for inversion
				if(gvdata[k][y-2][0][l]!=undef) fdata[k][y-1][0][l]=
				-(gvdata[k][y-1][0][l]*lcos[ystart-1+y-1]-gvdata[k][y-2][0][l]*lcos[ystart-1+y-2])
				/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f);
			}
		}
		
		return f;
	}
	
	/**
     * Calculate eddy absolute angular momentum (AAM) vertical flux convergence (VFC).
     *
     * @param	gawa	eddy AAM vertical flux = [g'w'] (m^2 Pa s^-2)
     *
     * @return	f		eddy AAM VFC (m^2 s^-2)
     */
	public Variable cEddyAAMVFC(Variable gawa){
		assignSubDomainParams(gawa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("aamVFC",gawa);
		f.setCommentAndUnit("eddy absolute angular momentum vertical flux convergence (m^2 s^-2)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] gadata=gawa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(gadata[l][1][j][0]!=undef&&gadata[l][0][j][0]!=undef)
				fdata[l][0][j][0]=-(gadata[l][1][j][0]-gadata[l][0][j][0])/dz;
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gadata[l][k+1][j][0]!=undef&&gadata[l][k-1][j][0]!=undef)
				fdata[l][k][j][0]=-(gadata[l][k+1][j][0]-gadata[l][k-1][j][0])/(dz*2f);
				
				// upper boundary (k==z-1), not used in inversion
				if(gadata[l][z-1][j][0]!=undef&&gadata[l][z-2][j][0]!=undef)
				fdata[l][z-1][j][0]=-(gadata[l][z-1][j][0]-gadata[l][z-2][j][0])/dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(gadata[1][j][0][l]!=undef&&gadata[0][j][0][l]!=undef)
				fdata[0][j][0][l]=-(gadata[1][j][0][l]-gadata[0][j][0][l])/dz;
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gadata[k+1][j][0][l]!=undef&&gadata[k-1][j][0][l]!=undef)
				fdata[k][j][0][l]=-(gadata[k+1][j][0][l]-gadata[k-1][j][0][l])/(dz*2f);
				
				// upper boundary (k==z-1), not used in inversion
				if(gadata[z-1][j][0][l]!=undef&&gadata[z-2][j][0][l]!=undef)
				fdata[z-1][j][0][l]=-(gadata[z-1][j][0][l]-gadata[z-2][j][0][l])/dz;
			}
		}
		
		return f;
	}
	
	/**
     * Calculate Eliassen-Palm (EP) horizontal flux convergence (HFC).
     *
     * @param	tava	eddy heat horizontal flux = [ï¿½ï¿½'v'] (K m s^-1)
     * @param	gava	eddy AAM  horizontal flux = [g'v'] (m^3 s^-2)
     * @param	gm		zonal mean absolute angular momentum (m^2 s^-1)
     * @param	thm		zonal mean potential temperature (K)
     *
     * @return	f		forcing factor (m^2 s^-2)
     */
	public Variable cEPHFC(Variable gava){
		assignSubDomainParams(gava);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("aamHFC",gava);
		f.setCommentAndUnit("eddy absolute angular momentum horizontal flux convergence (m^2 s^-2)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] gvdata=gava.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used for inversion
				if(gvdata[l][k][0][0]!=undef) fdata[l][k][0][0]=
				-(gvdata[l][k][0][0]*lcos[ystart-1]-gvdata[l][k][1][0]*lcos[ystart-1+1])
				/dy/((lcos[ystart-1]+lcos[ystart-1+1])/2f);
				
				// inner region, centered difference
				for(int j=1;j<y-1;j++)
				if(gvdata[l][k][j+1][0]!=undef&&gvdata[l][k][j-1][0]!=undef) fdata[l][k][j][0]=
				-(gvdata[l][k][j+1][0]*lcos[ystart-1+j+1]-gvdata[l][k][j-1][0]*lcos[ystart-1+j-1])
				/(dy*2f)/(lcos[ystart-1+j]);
				
				// north boundary (j==y-1), not used for inversion
				if(gvdata[l][k][y-2][0]!=undef) fdata[l][k][y-1][0]=
				-(gvdata[l][k][y-1][0]*lcos[ystart-1+y-1]-gvdata[l][k][y-2][0]*lcos[ystart-1+y-2])
				/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-1])/2f);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used for inversion
				if(gvdata[k][0][0][l]!=undef) fdata[k][0][0][l]=
				-(gvdata[k][0][0][l]*lcos[ystart-1]-gvdata[k][1][0][l]*lcos[ystart-1+1])
				/dy/((lcos[ystart-1]+lcos[ystart-1+1])/2f);
				
				// inner region, centered difference
				for(int j=1;j<y-1;j++)
				if(gvdata[k][j+1][0][l]!=undef&&gvdata[k][j-1][0][l]!=undef) fdata[k][j][0][l]=
				-(gvdata[k][j+1][0][l]*lcos[ystart-1+j+1]-gvdata[k][j-1][0][l]*lcos[ystart-1+j-1])
				/(dy*2f)/(lcos[ystart-1+j]);
				
				// north boundary (j==y-1), not used for inversion
				if(gvdata[k][y-2][0][l]!=undef) fdata[k][y-1][0][l]=
				-(gvdata[k][y-1][0][l]*lcos[ystart-1+y-1]-gvdata[k][y-2][0][l]*lcos[ystart-1+y-2])
				/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-1])/2f);
			}
		}
		
		return f;
	}
	
	/**
     * Calculate Eliassen-Palm (EP) vertical flux convergence (VFC).
     *
     * @param	tava	eddy heat horizontal flux = [ï¿½ï¿½'v'] (K m s^-1)
     * @param	gawa	eddy AAM vertical flux = [g'w'] (m^2 Pa s^-2)
     * @param	gm		zonal mean absolute angular momentum (m^2 s^-1)
     * @param	thm		zonal mean potential temperature (K)
     *
     * @return	f		force (m^2 s^-2)
     */
	public Variable cEPVFC(Variable gawa){
		assignSubDomainParams(gawa);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("aamVFC",gawa);
		f.setCommentAndUnit("eddy absolute angular momentum vertical flux convergence (m^2 s^-2)");
		f.setValue(undef);
		
		float[][][][]  fdata=   f.getData();
		float[][][][] gadata=gawa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(gadata[l][1][j][0]!=undef&&gadata[l][0][j][0]!=undef)
				fdata[l][0][j][0]=-(gadata[l][1][j][0]-gadata[l][0][j][0])/dz;
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gadata[l][k+1][j][0]!=undef&&gadata[l][k-1][j][0]!=undef)
				fdata[l][k][j][0]=-(gadata[l][k+1][j][0]-gadata[l][k-1][j][0])/(dz*2f);
				
				// upper boundary (k==z-1), not used in inversion
				if(gadata[l][z-1][j][0]!=undef&&gadata[l][z-2][j][0]!=undef)
				fdata[l][z-1][j][0]=-(gadata[l][z-1][j][0]-gadata[l][z-2][j][0])/dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(gadata[1][j][0][l]!=undef&&gadata[0][j][0][l]!=undef)
				fdata[0][j][0][l]=-(gadata[1][j][0][l]-gadata[0][j][0][l])/dz;
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gadata[k+1][j][0][l]!=undef&&gadata[k-1][j][0][l]!=undef)
				fdata[k][j][0][l]=-(gadata[k+1][j][0][l]-gadata[k-1][j][0][l])/(dz*2f);
				
				// upper boundary (k==z-1), not used in inversion
				if(gadata[z-1][j][0][l]!=undef&&gadata[z-2][j][0][l]!=undef)
				fdata[z-1][j][0][l]=-(gadata[z-1][j][0][l]-gadata[z-2][j][0][l])/dz;
			}
		}
		
		return f;
	}
	
	
	/**
     * Calculate eddy-induced streamfunction = [¦È'v']/[¦È]p.
     *
     * @param	tava	eddy heat horizontal flux = [¦È'v'] (K m s^-1)
     * @param	thm		zonal-mean potential temperature (K)
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
		
		if(tava.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(tvdata[l][0][j][0]!=undef&&tmdata[l][0][j][0]!=undef&&tmdata[l][1][j][0]!=undef)
				sfdata[l][0][j][0]=tvdata[l][0][j][0]/(tmdata[l][1][j][0]-tmdata[l][0][j][0])*dz;
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(tvdata[l][k][j][0]!=undef&&tmdata[l][k-1][j][0]!=undef&&tmdata[l][k+1][j][0]!=undef)
				sfdata[l][k][j][0]=tvdata[l][k][j][0]/(tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0])*(dz*2f);
				
				// upper boundary (k==z-1), not used in inversion
				if(tvdata[l][z-1][j][0]!=undef&&tmdata[l][z-1][j][0]!=undef&&tmdata[l][z-2][j][0]!=undef)
				sfdata[l][z-1][j][0]=tvdata[l][z-1][j][0]/(tmdata[l][z-1][j][0]-tmdata[l][z-2][j][0])*dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				// lower boundary (k==0), not used in inversion
				if(tvdata[0][j][0][l]!=undef&&tmdata[0][j][0][l]!=undef&&tmdata[1][j][0][l]!=undef)
				sfdata[0][j][0][l]=tvdata[0][j][0][l]/(tmdata[1][j][0][l]-tmdata[0][j][0][l])*dz;
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(tvdata[k][j][0][l]!=undef&&tmdata[k-1][j][0][l]!=undef&&tmdata[k+1][j][0][l]!=undef)
				sfdata[k][j][0][l]=tvdata[k][j][0][l]/(tmdata[k+1][j][0][l]-tmdata[k-1][j][0][l])*(dz*2f);
				
				// upper boundary (k==z-1), not used in inversion
				if(tvdata[z-1][j][0][l]!=undef&&tmdata[z-1][j][0][l]!=undef&&tmdata[z-2][j][0][l]!=undef)
				sfdata[z-1][j][0][l]=tvdata[z-1][j][0][l]/(tmdata[z-1][j][0][l]-tmdata[z-2][j][0][l])*dz;
			}
		}
		
		return sfe;
	}
	
	/**
     * Calculate eddy-induced velocity.
     *
     * @param	sfe		eddy-induced streamfunction (m^2 s^-2)
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
					wdata[l][k][0][0]=-(sdata[l][k][1][0]*lcos[ystart-1+1]-sdata[l][k][0][0]*lcos[ystart-1])
					/(dy*2f)/((lcos[ystart-1]+lcos[ystart-1+1])/2f);
					
					// inner region, centered difference
					for(int j=1,J=y-1;j<J;j++)
					if(sdata[l][k][j-1][0]!=undef&&sdata[l][k][j+1][0]!=undef)
					wdata[l][k][j][0]=-(sdata[l][k][j+1][0]*lcos[ystart-1+j+1]-sdata[l][k][j-1][0]*lcos[ystart-1+j-1])
					/(dy*2f)/lcos[ystart-1+j];
					
					// north boundary (j==y-1), not used in inversion
					if(sdata[l][k][y-1][0]!=undef&&sdata[l][k][y-2][0]!=undef)
					wdata[l][k][y-1][0]=-(sdata[l][k][y-1][0]*lcos[ystart-1+y-1]-sdata[l][k][y-2][0]*lcos[ystart-1+y-2])
					/(dy*2f)/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f);
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
					wdata[k][0][0][l]=-(sdata[k][1][0][l]*lcos[ystart-1+1]-sdata[k][0][0][l]*lcos[ystart-1])
					/(dy*2f)/((lcos[ystart-1]+lcos[ystart-1+1])/2f);
					
					// inner region, centered difference
					for(int j=1,J=y-1;j<J;j++)
					if(sdata[k][j-1][0][l]!=undef&&sdata[k][j+1][0][l]!=undef)
					wdata[k][j][0][l]=-(sdata[k][j+1][0][l]*lcos[ystart-1+j+1]-sdata[k][j-1][0][l]*lcos[ystart-1+j-1])
					/(dy*2f)/lcos[ystart-1+j];
					
					// north boundary (j==y-1), not used in inversion
					if(sdata[k][y-1][0][l]!=undef&&sdata[k][y-2][0][l]!=undef)
					wdata[k][y-1][0][l]=-(sdata[k][y-1][0][l]*lcos[ystart-1+y-1]-sdata[k][y-2][0][l]*lcos[ystart-1+y-2])
					/(dy*2f)/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f);
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
     * @param	gm		zonal mean absolute angular momentum (m^2 s^-1)
     * @param	thm		zonal mean potential temperature (K)
     * @param	rZToY	ratio of Z to Y: Z/Y in the plot used to scale the vector components (e.g., 0.8)
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
		float[][][][] tadata= tava.getData();
		float[][][][] gvdata= gava.getData();
		float[][][][] gwdata= gawa.getData();
		float[][][][] gmdata=   gm.getData();
		float[][][][] tmdata=  thm.getData();
		
		if(tava.isTFirst()){
			for(int l=0;l<t;l++){
				/*** meridional component ***/
				for(int j=0;j<y;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[l][0][j][0]!=undef&&gwdata[l][0][j][0]!=undef&&tadata[l][0][j][0]!=undef)
					mcdata[l][0][j][0]=-gvdata[l][0][j][0]+tadata[l][0][j][0]*
					(gmdata[l][1][j][0]-gmdata[l][0][j][0])/(tmdata[l][1][j][0]-tmdata[l][0][j][0]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[l][k][j][0]!=undef&&gwdata[l][k][j][0]!=undef&&tadata[l][k][j][0]!=undef)
					mcdata[l][k][j][0]=-gvdata[l][k][j][0]+tadata[l][k][j][0]*
					(gmdata[l][k+1][j][0]-gmdata[l][k-1][j][0])/(tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[l][z-1][j][0]!=undef&&gwdata[l][z-1][j][0]!=undef&&tadata[l][z-1][j][0]!=undef)
					mcdata[l][z-1][j][0]=-gvdata[l][z-1][j][0]+tadata[l][z-1][j][0]*
					(gmdata[l][z-1][j][0]-gmdata[l][z-2][j][0])/(tmdata[l][z-1][j][0]-tmdata[l][z-2][j][0]);
				}
				
				/*** vertical component ***/
				// lower south boundary (k==0, j==0), not used in inversion
				if(gvdata[l][0][0][0]!=undef&&gwdata[l][0][0][0]!=undef&&tadata[l][0][0][0]!=undef)
				vcdata[l][0][0][0]=-gwdata[l][0][0][0]-tadata[l][0][0][0]*dz/dy*
				(gmdata[l][0][1][0]-gmdata[l][0][0][0])/(tmdata[l][1][0][0]-tmdata[l][0][0][0]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[l][k][0][0]!=undef&&gwdata[l][k][0][0]!=undef&&tadata[l][k][0][0]!=undef)
				vcdata[l][k][0][0]=-gwdata[l][k][0][0]-tadata[l][k][0][0]*2f*dz/dy*
				(gmdata[l][k][1][0]-gmdata[l][k][0][0])/(tmdata[l][k+1][0][0]-tmdata[l][k-1][0][0]);
				
				// upper south boundary (k==z-1, j==0), not used in inversion
				if(gvdata[l][z-1][0][0]!=undef&&gwdata[l][z-1][0][0]!=undef&&tadata[l][z-1][0][0]!=undef)
				vcdata[l][z-1][0][0]=-gwdata[l][z-1][0][0]-tadata[l][z-1][0][0]*dz/dy*
				(gmdata[l][z-1][1][0]-gmdata[l][z-1][0][0])/(tmdata[l][z-1][0][0]-tmdata[l][z-2][0][0]);
				
				for(int j=1,J=y-1;j<J;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[l][0][j][0]!=undef&&gwdata[l][0][j][0]!=undef&&tadata[l][0][j][0]!=undef)
					vcdata[l][0][j][0]=-gwdata[l][0][j][0]-tadata[l][0][j][0]*dz/dy/2f*
					(gmdata[l][0][j+1][0]-gmdata[l][0][j-1][0])/(tmdata[l][1][j][0]-tmdata[l][0][j][0]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[l][k][j][0]!=undef&&gwdata[l][k][j][0]!=undef&&tadata[l][k][j][0]!=undef)
					vcdata[l][k][j][0]=-gwdata[l][k][j][0]-tadata[l][k][j][0]*dz/dy*
					(gmdata[l][k][j+1][0]-gmdata[l][k][j-1][0])/(tmdata[l][k+1][j][0]-tmdata[l][k-1][j][0]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[l][z-1][j][0]!=undef&&gwdata[l][z-1][j][0]!=undef&&tadata[l][z-1][j][0]!=undef)
					vcdata[l][z-1][j][0]=-gwdata[l][z-1][j][0]-tadata[l][z-1][j][0]*dz/dy/2f*
					(gmdata[l][z-1][j+1][0]-gmdata[l][z-1][j-1][0])/(tmdata[l][z-1][j][0]-tmdata[l][z-2][j][0]);
				}
				
				// lower north boundary (k==0,j==y-1), not used in inversion
				if(gvdata[l][0][y-1][0]!=undef&&gwdata[l][0][y-1][0]!=undef&&tadata[l][0][y-1][0]!=undef)
				vcdata[l][0][y-1][0]=-gwdata[l][0][y-1][0]-tadata[l][0][y-1][0]*dz/dy*
				(gmdata[l][0][y-1][0]-gmdata[l][0][y-2][0])/(tmdata[l][1][y-1][0]-tmdata[l][0][y-1][0]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[l][k][y-1][0]!=undef&&gwdata[l][k][y-1][0]!=undef&&tadata[l][k][y-1][0]!=undef)
				vcdata[l][k][y-1][0]=-gwdata[l][k][y-1][0]-tadata[l][k][y-1][0]*2f*dz/dy*
				(gmdata[l][k][y-1][0]-gmdata[l][k][y-2][0])/(tmdata[l][k+1][y-1][0]-tmdata[l][k-1][y-1][0]);
				
				// upper north boundary (k==z-1,j==y-1), not used in inversion
				if(gvdata[l][z-1][y-1][0]!=undef&&gwdata[l][z-1][y-1][0]!=undef&&tadata[l][z-1][y-1][0]!=undef)
				vcdata[l][z-1][y-1][0]=-gwdata[l][z-1][y-1][0]-tadata[l][z-1][y-1][0]*dz/dy*
				(gmdata[l][z-1][y-1][0]-gmdata[l][z-1][y-2][0])/(tmdata[l][z-1][y-1][0]-tmdata[l][z-2][y-1][0]);
			}
			
		}else{
			for(int l=0;l<t;l++){
				/*** meridional component ***/
				for(int j=0;j<y;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[0][j][0][l]!=undef&&gwdata[0][j][0][l]!=undef&&tadata[0][j][0][l]!=undef)
					mcdata[0][j][0][l]=-gvdata[0][j][0][l]+tadata[0][j][0][l]*
					(gmdata[1][j][0][l]-gmdata[0][j][0][l])/(tmdata[1][j][0][l]-tmdata[0][j][0][l]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[k][j][0][l]!=undef&&gwdata[k][j][0][l]!=undef&&tadata[k][j][0][l]!=undef)
					mcdata[k][j][0][l]=-gvdata[k][j][0][l]+tadata[k][j][0][l]*
					(gmdata[k+1][j][0][l]-gmdata[k-1][j][0][l])/(tmdata[k+1][j][0][l]-tmdata[k-1][j][0][l]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[z-1][j][0][l]!=undef&&gwdata[z-1][j][0][l]!=undef&&tadata[z-1][j][0][l]!=undef)
					mcdata[z-1][j][0][l]=-gvdata[z-1][j][0][l]+tadata[z-1][j][0][l]*
					(gmdata[z-1][j][0][l]-gmdata[z-2][j][0][l])/(tmdata[z-1][j][0][l]-tmdata[z-2][j][0][l]);
				}
				
				/*** vertical component ***/
				// lower south boundary (k==0, j==0), not used in inversion
				if(gvdata[0][0][0][l]!=undef&&gwdata[0][0][0][l]!=undef&&tadata[0][0][0][l]!=undef)
				vcdata[0][0][0][l]=-gwdata[0][0][0][l]-tadata[0][0][0][l]*dz/dy*
				(gmdata[0][1][0][l]-gmdata[0][0][0][l])/(tmdata[1][0][0][l]-tmdata[0][0][0][l]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[k][0][0][l]!=undef&&gwdata[k][0][0][l]!=undef&&tadata[k][0][0][l]!=undef)
				vcdata[k][0][0][l]=-gwdata[k][0][0][l]-tadata[k][0][0][l]*2f*dz/dy*
				(gmdata[k][1][0][l]-gmdata[k][0][0][l])/(tmdata[k+1][0][0][l]-tmdata[k-1][0][0][l]);
				
				// upper south boundary (k==z-1, j==0), not used in inversion
				if(gvdata[z-1][0][0][l]!=undef&&gwdata[z-1][0][0][l]!=undef&&tadata[z-1][0][0][l]!=undef)
				vcdata[z-1][0][0][l]=-gwdata[z-1][0][0][l]-tadata[z-1][0][0][l]*dz/dy*
				(gmdata[z-1][1][0][l]-gmdata[z-1][0][0][l])/(tmdata[z-1][0][0][l]-tmdata[z-2][0][0][l]);
				
				for(int j=1,J=y-1;j<J;j++){
					// lower boundary (k==0), not used in inversion
					if(gvdata[0][j][0][l]!=undef&&gwdata[0][j][0][l]!=undef&&tadata[0][j][0][l]!=undef)
					vcdata[0][j][0][l]=-gwdata[0][j][0][l]-tadata[0][j][0][l]*dz/dy/2f*
					(gmdata[0][j+1][0][l]-gmdata[0][j-1][0][l])/(tmdata[1][j][0][l]-tmdata[0][j][0][l]);
					
					// inner region, centered difference
					for(int k=1,K=z-1;k<K;k++)
					if(gvdata[k][j][0][l]!=undef&&gwdata[k][j][0][l]!=undef&&tadata[k][j][0][l]!=undef)
					vcdata[k][j][0][l]=-gwdata[k][j][0][l]-tadata[k][j][0][l]*dz/dy*
					(gmdata[k][j+1][0][l]-gmdata[k][j-1][0][l])/(tmdata[k+1][j][0][l]-tmdata[k-1][j][0][l]);
					
					// upper boundary (k==z-1), not used in inversion
					if(gvdata[z-1][j][0][l]!=undef&&gwdata[z-1][j][0][l]!=undef&&tadata[z-1][j][0][l]!=undef)
					vcdata[z-1][j][0][l]=-gwdata[z-1][j][0][l]-tadata[z-1][j][0][l]*dz/dy/2f*
					(gmdata[z-1][j+1][0][l]-gmdata[z-1][j-1][0][l])/(tmdata[z-1][j][0][l]-tmdata[z-2][j][0][l]);
				}
				
				// lower north boundary (k==0,j==y-1), not used in inversion
				if(gvdata[0][y-1][0][l]!=undef&&gwdata[0][y-1][0][l]!=undef&&tadata[0][y-1][0][l]!=undef)
				vcdata[0][y-1][0][l]=-gwdata[0][y-1][0][l]-tadata[0][y-1][0][l]*dz/dy*
				(gmdata[0][y-1][0][l]-gmdata[0][y-2][0][l])/(tmdata[1][y-1][0][l]-tmdata[0][y-1][0][l]);
				
				// inner region, centered difference
				for(int k=1,K=z-1;k<K;k++)
				if(gvdata[k][y-1][0][l]!=undef&&gwdata[k][y-1][0][l]!=undef&&tadata[k][y-1][0][l]!=undef)
				vcdata[k][y-1][0][l]=-gwdata[k][y-1][0][l]-tadata[k][y-1][0][l]*2f*dz/dy*
				(gmdata[k][y-1][0][l]-gmdata[k][y-2][0][l])/(tmdata[k+1][y-1][0][l]-tmdata[k-1][y-1][0][l]);
				
				// upper north boundary (k==z-1,j==y-1), not used in inversion
				if(gvdata[z-1][y-1][0][l]!=undef&&gwdata[z-1][y-1][0][l]!=undef&&tadata[z-1][y-1][0][l]!=undef)
				vcdata[z-1][y-1][0][l]=-gwdata[z-1][y-1][0][l]-tadata[z-1][y-1][0][l]*dz/dy*
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
		float denomY=Math.abs(ydef[y-1]-ydef[0])*REarth;
		
		System.out.println(dz);
		System.out.println(dy);
		
		System.out.println(denomZ+"\t"+denomY+"\t"+rZToY);
		
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
     * implement the methods in EllipticalInterface
     */
	public Variable cAPrime(Variable sigm){
		assignSubDomainParams(sigm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable A=new Variable("Ap",sigm);
		A.setCommentAndUnit("coefficient A' of elliptic equation");
		A.setValue(undef);
		
		float[][][][] Adata=   A.getData();
		float[][][][] Sdata=sigm.getData();
		
		if(sigm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(Sdata[l][k][j][0]!=undef)
			Adata[l][k][j][0]=(float)((Sdata[l][k][j-1][0]+Sdata[l][k][j][0])/2.0/cos((ydef[j-1]+ydef[j])/2.0));
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(Sdata[k][j][0][l]!=undef)
			Adata[k][j][0][l]=(float)((Sdata[k][j-1][0][l]+Sdata[k][j][0][l])/2.0/cos((ydef[j-1]+ydef[j])/2.0));
		}
		
		return A;
	}
    
	public Variable cBPrime(Variable thm){
		assignSubDomainParams(thm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable B=new Variable("Bp",thm);
		B.setCommentAndUnit("coefficient B' of elliptic equation");
		B.setValue(undef);
		
		float[][][][] Bdata=  B.getData();
		float[][][][] tdata=thm.getData();
		
		if(thm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used in inversion
				if(tdata[l][k][0][0]!=undef&&tdata[l][k][1][0]!=undef)
				Bdata[l][k][0][0]=(tdata[l][k][1][0]-tdata[l][k][0][0])/dy/((lcos[ystart-1]+lcos[ystart-1+1])/2f)*
				Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
				
				// inner region
				for(int j=1,J=y-1;j<J;j++)
				if(tdata[l][k][j+1][0]!=undef&&tdata[l][k][j-1][0]!=undef)
				Bdata[l][k][j][0]=(tdata[l][k][j+1][0]-tdata[l][k][j-1][0])/(2f*dy)/lcos[ystart-1+j]*
				Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
				
				// north boundary (j==y-1), not used in inversion
				if(tdata[l][k][y-1][0]!=undef&&tdata[l][k][y-2][0]!=undef)
				Bdata[l][k][y-1][0]=(tdata[l][k][y-1][0]-tdata[l][k][y-2][0])/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f)*
				Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used in inversion
				if(tdata[k][0][0][l]!=undef&&tdata[k][1][0][l]!=undef)
				Bdata[k][0][0][l]=(tdata[k][1][0][l]-tdata[k][0][0][l])/dy/((lcos[ystart-1]+lcos[ystart-1+1])/2f)*
				Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
				
				// inner region
				for(int j=1,J=y-1;j<J;j++)
				if(tdata[k][j+1][0][l]!=undef&&tdata[k][j-1][0][l]!=undef)
				Bdata[k][j][0][l]=(tdata[k][j+1][0][l]-tdata[k][j-1][0][l])/(2f*dy)/lcos[ystart-1+j]*
				Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
				
				// north boundary (j==y-1), not used in inversion
				if(tdata[k][y-1][0][l]!=undef&&tdata[k][y-2][0][l]!=undef)
				Bdata[k][y-1][0][l]=(tdata[k][y-1][0][l]-tdata[k][y-2][0][l])/dy/((lcos[ystart-1+y-1]+lcos[ystart-1+y-2])/2f)*
				Rd*(float)pow(zdef[zstart-1+k]/100000.0,kappa)/zdef[zstart-1+k];
			}
		}
		
		return B;
	}
	
	public Variable cCPrime(Variable Mm){
		assignSubDomainParams(Mm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable C=new Variable("Cp",Mm);
		C.setCommentAndUnit("coefficient C' of elliptic equation");
		C.setValue(undef);
		
		float[]         buf=new float[y];
		float[][][][] Cdata= C.getData();
		float[][][][] Mdata=Mm.getData();
		
		if(Mm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++){
				for(int j=0;j<y;j++)
				if(Mdata[l][k-1][j][0]!=undef&&Mdata[l][k][j][0]!=undef)
				buf[j]=(Mdata[l][k-1][j][0]+Mdata[l][k][j][0])/2f;
				
				// south boundary (j==0), not used in inversion
				if(buf[0]!=undef&&buf[1]!=undef)
				Cdata[l][k][0][0]=-2f*buf[0]/(float)pow(rcos[ystart-1],3.0)*
				(buf[1]-buf[0])/dy*ltan[ystart-1];
				
				// inner region
				for(int j=1,J=y-1;j<J;j++) if(buf[j]!=undef&&buf[j+1]!=undef&&buf[j-1]!=undef)
				Cdata[l][k][j][0]=-2f*buf[j]/(float)pow(rcos[ystart-1+j],3.0)*
				(buf[j+1]-buf[j-1])/(2f*dy)*ltan[ystart-1+j];
				
				// north boundary (j==y-1), not used in inversion
				if(buf[y-1]!=undef&&buf[y-2]!=undef)
				Cdata[l][k][y-1][0]=-2f*buf[y-1]/(float)pow(rcos[ystart-1+y-1],3.0)*
				(buf[y-1]-buf[y-2])/dy*ltan[ystart-1+y-1];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++){
				for(int j=0;j<y;j++)
				if(Mdata[k-1][j][0][l]!=undef&&Mdata[k][j][0][l]!=undef)
				buf[j]=(Mdata[k-1][j][0][l]+Mdata[k][j][0][l])/2f;
				
				// south boundary (j==0), not used in inversion
				if(buf[0]!=undef&&buf[1]!=undef)
				Cdata[k][0][0][l]=-2f*buf[0]/(float)pow(rcos[ystart-1],3.0)*
				(buf[1]-buf[0])/dy*ltan[ystart-1];
				
				// inner region
				for(int j=1,J=y-1;j<J;j++) if(buf[j]!=undef&&buf[j+1]!=undef&&buf[j-1]!=undef)
				Cdata[k][j][0][l]=-2f*buf[j]/(float)pow(rcos[ystart-1+j],3.0)*
				(buf[j+1]-buf[j-1])/(2f*dy)*ltan[ystart-1+j];
				
				// north boundary (j==y-1), not used in inversion
				if(buf[y-1]!=undef&&buf[y-2]!=undef)
				Cdata[k][y-1][0][l]=-2f*buf[y-1]/(float)pow(rcos[ystart-1+y-1],3.0)*
				(buf[y-1]-buf[y-2])/dy*ltan[ystart-1+y-1];
			}
		}
		
		return C;
	}
    
    public Variable cA(Variable sigm){
		assignSubDomainParams(sigm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable A=new Variable("A",sigm);
		A.setCommentAndUnit("coefficient A of elliptic equation");
		A.setValue(undef);
		
		float[][][][] Adata=   A.getData();
		float[][][][] Sdata=sigm.getData();
		
		if(sigm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(Sdata[l][k][j][0]!=undef&&lcos[ystart-1+j]!=0) Adata[l][k][j][0]=Sdata[l][k][j][0]/lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			if(Sdata[k][j][0][l]!=undef&&lcos[ystart-1+j]!=0) Adata[k][j][0][l]=Sdata[k][j][0][l]/lcos[ystart-1+j];
		}
		
		return A;
	}
	
	public Variable cB(Variable am){
		Variable B=cBPrime(am);
		
		B.setName("Bb");
		
		return B;
	}
	
	public Variable cC(Variable Mm){
		assignSubDomainParams(Mm);
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable C=new Variable("C",Mm);
		C.setCommentAndUnit("coefficient C of elliptic equation");
		C.setValue(undef);
		
		float[][][][] Cdata= C.getData();
		float[][][][] Mdata=Mm.getData();
		
		if(Mm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used in inversion
				if(Mdata[l][k][0][0]!=undef&&Mdata[l][k][1][0]!=undef)
				Cdata[l][k][0][0]=-2f*Mdata[l][k][0][0]/(float)pow(rcos[ystart-1],3.0)*ltan[ystart-1]*
				(Mdata[l][k][1][0]-Mdata[l][k][0][0])/dy;
				
				// inner region
				for(int j=1,J=y-1;j<J;j++)
				if(Mdata[l][k][j][0]!=undef&&Mdata[l][k][j+1][0]!=undef&&Mdata[l][k][j-1][0]!=undef)
				Cdata[l][k][j][0]=-2f*Mdata[l][k][j][0]/(float)pow(rcos[ystart-1+j],3.0)*ltan[ystart-1+j]*
				(Mdata[l][k][j+1][0]-Mdata[l][k][j-1][0])/(dy+dy);
				
				// north boundary (j==y-1), not used in inversion
				if(Mdata[l][k][y-1][0]!=undef&&Mdata[l][k][y-2][0]!=undef)
				Cdata[l][k][y-1][0]=-2f*Mdata[l][k][y-1][0]/(float)pow(rcos[ystart-1+y-1],3.0)*ltan[ystart-1+y-1]*
				(Mdata[l][k][y-1][0]-Mdata[l][k][y-2][0])/dy;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				// south boundary (j==0), not used in inversion
				if(Mdata[k][0][0][l]!=undef&&Mdata[k][1][0][l]!=undef)
				Cdata[k][0][0][l]=-2f*Mdata[k][0][0][l]/(float)pow(rcos[ystart-1],3.0)*ltan[ystart-1]*
				(Mdata[k][1][0][l]-Mdata[k][0][0][l])/dy;
				
				// inner region
				for(int j=1,J=y-1;j<J;j++)
				if(Mdata[k][j][0][l]!=undef&&Mdata[k][j+1][0][l]!=undef&&Mdata[k][j-1][0][l]!=undef)
				Cdata[k][j][0][l]=-2f*Mdata[k][j][0][l]/(float)pow(rcos[ystart-1+j],3.0)*ltan[ystart-1+j]*
				(Mdata[k][j+1][0][l]-Mdata[k][j-1][0][l])/(dy+dy);
				
				// north boundary (j==y-1), not used in inversion
				if(Mdata[k][y-1][0][l]!=undef&&Mdata[k][y-2][0][l]!=undef)
				Cdata[k][y-1][0][l]=-2f*Mdata[k][y-1][0][l]/(float)pow(rcos[ystart-1+y-1],3.0)*ltan[ystart-1+y-1]*
				(Mdata[k][y-1][0][l]-Mdata[k][y-2][0][l])/dy;
			}
		}
		
		return C;
	}
	
	
	/*** helper methods ***/
	
	
	/** test
	public static void main(String[] args){
		try{
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
