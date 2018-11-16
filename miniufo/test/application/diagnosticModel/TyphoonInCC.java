/**
 * @(#)TyphoonInCC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.diagnosticModel;

import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.io.Print;
import miniufo.test.application.EllipticEquationInCC;
import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static miniufo.geophysics.atmos.ThermoDynamics.cL;
import static miniufo.geophysics.atmos.ThermoDynamics.Cp;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;
import static miniufo.basic.ArrayUtil.getAbsMax;
import static miniufo.diagnosis.SpatialModel.cSphericalDistanceByRadian;
import static miniufo.diagnosis.SpatialModel.REarth;
import static miniufo.diagnosis.SpatialModel.omegaEarth;


/**
 * Typhoon model in cylindrical coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class TyphoonInCC extends EllipticEquationInCC implements Print{
	//
	public static final float Cd=1e-3f;
	public static final float Kh=5e+3f;	// m^2 s^-1
	public static final float Kv=1e+3f;	// vertical eddy coefficient in p-coordinate
	
	private boolean print=true;
	
	private short[][][] tags=null;
	
	private float[] cs=null;	// central speed of the model
	private float[] cd=null;	// central direction of the speed
	
	private float[][] wtan=null;	// whole moving velocity projected to tangential direction
	private float[][] wnor=null;	// whole moving velocity projected to radial direction
	
	
	/**
     * constructor
     *
     * @param	csm		initialized by spacial model in Cylindrical coordinate
     */
	public TyphoonInCC(CylindricalSpatialModel csm){
		super(csm);
		
		tags=new short[csm.getTCount()][csm.getYCount()][csm.getXCount()];
	}
	
	/**
     * constructor
     *
     * @param	csm		initialized by spacial model in Cylindrical coordinate
     * @param	sfp		surface pressure
     */
	public TyphoonInCC(CylindricalSpatialModel csm,Variable sfp){ super(csm); cTopography(sfp);}
	
	
	/**
     * calculate the tangential friction
     *
     * @param	ut10	tangential wind at 10 meter
     * @param	vr10	radial wind at 10 meter
     * @param	pblh	planet boundary layer height
     * @param	sfh		surface potential height
     * @param	hgt		potential height
     *
     * @return	fric	tangential friction
     */
	public Variable cSurfaceFriction(Variable ut10,Variable vr10,Variable pblh,Variable sfh,Variable hgt){
		if(!ut10.isLike(vr10)||!ut10.isLike(pblh)||!ut10.isLike(sfh)||!ut10.isAreaLike(hgt))
		throw new IllegalArgumentException("dimensions not same");
		
		t=ut10.getTCount();	z=ut10.getZCount();	y=ut10.getYCount();	x=ut10.getXCount();
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		z=hgt.getZCount();
		
		Variable fsfc=new Variable("fsfc",hgt);	fsfc.setCommentAndUnit("surface friction");
		
		float undef=hgt.getUndef();	fsfc.setUndef(undef);
		
		float[][][][] ut10data=ut10.getData();
		float[][][][] vr10data=vr10.getData();
		float[][][][] pblhdata=pblh.getData();
		float[][][][] hgtdata = hgt.getData();
		float[][][][] sfhdata = sfh.getData();
		float[][][][] fricdata=fsfc.getData();
		
		if(ut10.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) fricdata[l][k][j][i]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[l][k][j][i]-sfhdata[l][0][j][i];
					
					if(deltz>pblhdata[l][0][j][i]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[l][0][j][i])>=(pblhdata[l][0][j][i]-hgtdata[l][k-1][j][i]))
							kmax--;
						
						float tmp=-Cd*(float)Math.hypot(ut10data[l][0][j][i],vr10data[l][0][j][i])*
						ut10data[l][0][j][i]/pblhdata[l][0][j][i];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) fricdata[l][k2][j][i]=tmp;
						
						break;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) fricdata[k][j][i][l]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[k][j][i][l]-sfhdata[0][j][i][l];
					
					if(deltz>pblhdata[0][j][i][l]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[0][j][i][l])>=(pblhdata[0][j][i][l]-hgtdata[k-1][j][i][l]))
							kmax--;
						
						float tmp=-Cd*(float)Math.hypot(ut10data[0][j][i][l],vr10data[0][j][i][l])*
						ut10data[0][j][i][l]/pblhdata[0][j][i][l];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) fricdata[k2][j][i][l]=tmp;
						
						break;
					}
				}
			}
		}
		
		return fsfc;
	}
	
	/**
     * calculate the tangential friction in free atmosphere
     * due to vertical momentum dissipation
     *
     * @param	ut		tangential wind
     *
     * @return	fric	tangential friction
     */
	public Variable cFriction(Variable ut){
		t=ut.getTCount();	z=ut.getZCount();	y=ut.getYCount();	x=ut.getXCount();
		
		Variable fric=new Variable("fric",ut);	fric.setCommentAndUnit("friction in free atmosphere");
		
		float[] buf1=new float[z-1];
		float[][][][] udata=  ut.getData();
		float[][][][] fdata=fric.getData();
		
		if(fric.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<z-1;k++)
				buf1[k]=Kv*(zdef[k+1]+zdef[k])/2/zdef[0]*(udata[l][k+1][j][i]-udata[l][k][j][i])/dz;
				
				for(int k=1;k<z-1;k++)
				fdata[l][k][j][i]=(buf1[k]-buf1[k-1])/dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<z-1;k++)
				buf1[k]=Kv*(zdef[k+1]+zdef[k])/2/zdef[0]*(udata[k+1][j][i][l]-udata[k][j][i][l])/dz;
				
				for(int k=1;k<z-1;k++)
				fdata[k][j][i][l]=(buf1[k]-buf1[k-1])/dz;
			}
		}
		
		return fric;
	}
	
	/**
     * calculate the sensible/latent heat flux using surface heat flux data
     *
     * @param	flux	flux data
     * @param	T2m	 	temperature at 2 meter
     * @param	T		temperature
     * @param	pblh	planet boundary layer height
     * @param	sfh		surface potential height
     * @param	hgt		potential height
     *
     * @return	htfl	heat flux
     */
	public Variable cHeatFlux(Variable flux,Variable T2m,Variable pblh,Variable sfh,Variable hgt){
		if(!flux.isLike(T2m)||!flux.isLike(pblh)||!flux.isLike(T2m)||!flux.isLike(sfh)||!flux.isAreaLike(hgt))
		throw new IllegalArgumentException("dimensions not same");
		
		t=flux.getTCount();	z=flux.getZCount();	y=flux.getYCount();	x=flux.getXCount();
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		z=hgt.getZCount();
		
		Variable htfl=new Variable("htfl",hgt);	htfl.setCommentAndUnit("heat flux");
		
		float undef=hgt.getUndef();	htfl.setUndef(undef);
		
		float[][][][] fluxdata=flux.getData();
		float[][][][] T2mdata = T2m.getData();
		float[][][][] pblhdata=pblh.getData();
		float[][][][] hgtdata = hgt.getData();
		float[][][][] sfhdata = sfh.getData();
		float[][][][] htfldata=htfl.getData();
		
		if(flux.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) htfldata[l][k][j][i]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[l][k][j][i]-sfhdata[l][0][j][i];
					
					if(deltz>pblhdata[l][0][j][i]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[l][0][j][i])>=(pblhdata[l][0][j][i]-hgtdata[l][k-1][j][i]))
							kmax--;
						
						float tmp=-Rd*T2mdata[l][0][j][i]*fluxdata[l][0][j][i]/zdef[tags[l][j][i]]/pblhdata[l][0][j][i];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) htfldata[l][k2][j][i]=tmp;
						
						break;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) htfldata[k][j][i][l]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[k][j][i][l]-sfhdata[0][j][i][l];
					
					if(deltz>pblhdata[0][j][i][l]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[0][j][i][l])>=(pblhdata[0][j][i][l]-hgtdata[k-1][j][i][l]))
							kmax--;
						
						float tmp=-Rd*T2mdata[0][j][i][l]*fluxdata[0][j][i][l]/zdef[tags[l][j][i]]/pblhdata[0][j][i][l];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) htfldata[k2][j][i][l]=tmp;
						
						break;
					}
				}
			}
		}
		
		return htfl;
	}
	
	/**
     * calculate the sensible heat flux
     *
     * @param	ut10	tangential wind at 10 meter
     * @param	vr10	radial wind at 10 meter
     * @param	T2m	 	temperature at 2 meter
     * @param	T		temperature
     * @param	pblh	planet boundary layer height
     * @param	sfh		surface potential height
     * @param	hgt		potential height
     *
     * @return	shfl	sensible heat flux
     */
	public Variable cSensibleHeatFlux
	(Variable ut10,Variable vr10,Variable T2m,Variable T,Variable pblh,Variable sfh,Variable hgt){
		if(!ut10.isLike(vr10)||!ut10.isLike(pblh)||!ut10.isLike(T2m)||!ut10.isLike(sfh)||!ut10.isAreaLike(hgt))
		throw new IllegalArgumentException("dimensions not same");
		
		t=ut10.getTCount();	z=ut10.getZCount();	y=ut10.getYCount();	x=ut10.getXCount();
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		z=hgt.getZCount();
		
		Variable shfl=new Variable("shfl",hgt);	shfl.setCommentAndUnit("sensible heat flux");
		
		float undef=hgt.getUndef();	shfl.setUndef(undef);
		
		float[][][][] ut10data=ut10.getData();
		float[][][][] vr10data=vr10.getData();
		float[][][][] T2mdata = T2m.getData();
		float[][][][] Tdata   =   T.getData();
		float[][][][] pblhdata=pblh.getData();
		float[][][][] hgtdata = hgt.getData();
		float[][][][] sfhdata = sfh.getData();
		float[][][][] sflxdata=shfl.getData();
		
		if(ut10.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) sflxdata[l][k][j][i]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[l][k][j][i]-sfhdata[l][0][j][i];
					
					if(deltz>pblhdata[l][0][j][i]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[l][0][j][i])>=(pblhdata[l][0][j][i]-hgtdata[l][k-1][j][i]))
							kmax--;
						
						float tmp=-Cd*(float)Math.hypot(ut10data[l][0][j][i],vr10data[l][0][j][i])*
							Cp*(T2mdata[l][0][j][i]-Tdata[l][tags[l][j][i]][j][i])/pblhdata[l][0][j][i];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) sflxdata[l][k2][j][i]=tmp;
						
						break;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) sflxdata[k][j][i][l]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[k][j][i][l]-sfhdata[0][j][i][l];
					
					if(deltz>pblhdata[0][j][i][l]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[0][j][i][l])>=(pblhdata[0][j][i][l]-hgtdata[k-1][j][i][l]))
							kmax--;
						
						float tmp=-Cd*(float)Math.hypot(ut10data[0][j][i][l],vr10data[0][j][i][l])*
							Cp*(T2mdata[0][j][i][l]-Tdata[tags[l][j][i]][j][i][l])/pblhdata[0][j][i][l];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) sflxdata[k2][j][i][l]=tmp;
						
						break;
					}
				}
			}
		}
		
		return shfl;
	}
	
	/**
     * calculate the latent heat flux
     *
     * @param	ut10	tangential wind at 10 meter
     * @param	vr10	radial wind at 10 meter
     * @param	q2m	 	specific humidity at 2 meter
     * @param	q		specific humidity
     * @param	T		temperature
     * @param	pblh	planet boundary layer height
     * @param	sfh		surface potential height
     * @param	hgt		potential height
     *
     * @return	lhfl	latent heat flux
     */
	public Variable cLatentHeatFlux
	(Variable ut10,Variable vr10,Variable q2m,Variable q,Variable T,Variable pblh,Variable sfh,Variable hgt){
		if(!ut10.isLike(vr10)||!ut10.isLike(pblh)||!ut10.isLike(q2m)||!ut10.isLike(sfh)||!ut10.isAreaLike(hgt))
		throw new IllegalArgumentException("dimensions not same");
		
		t=ut10.getTCount();	z=ut10.getZCount();	y=ut10.getYCount();	x=ut10.getXCount();
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		z=hgt.getZCount();
		
		Variable lhfl=new Variable("lhfl",hgt);	lhfl.setCommentAndUnit("latent heat flux");
		
		float undef=hgt.getUndef();	lhfl.setUndef(undef);
		
		float[][][][] ut10data=ut10.getData();
		float[][][][] vr10data=vr10.getData();
		float[][][][] q2mdata = q2m.getData();
		float[][][][] qdata   =   q.getData();
		float[][][][] Tdata   =   T.getData();
		float[][][][] pblhdata=pblh.getData();
		float[][][][] hgtdata = hgt.getData();
		float[][][][] sfhdata = sfh.getData();
		float[][][][] lhfldata=lhfl.getData();
		
		if(ut10.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) lhfldata[l][k][j][i]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[l][k][j][i]-sfhdata[l][0][j][i];
					
					if(deltz>pblhdata[l][0][j][i]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[l][0][j][i])>=(pblhdata[l][0][j][i]-hgtdata[l][k-1][j][i]))
							kmax--;
						
						float tmp=-Cd*(float)Math.hypot(ut10data[l][0][j][i],vr10data[l][0][j][i])*
							cL(Tdata[l][tags[l][j][i]][j][i])*(q2mdata[l][0][j][i]-qdata[l][tags[l][j][i]][j][i])/pblhdata[l][0][j][i];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) lhfldata[l][k2][j][i]=tmp;
						
						break;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int k=0;k<tags[l][j][i];k++) lhfldata[k][j][i][l]=undef;
				
				for(int k=tags[l][j][i];k<z;k++){
					float deltz=hgtdata[k][j][i][l]-sfhdata[0][j][i][l];
					
					if(deltz>pblhdata[0][j][i][l]){
						int kmax=k;
						
						if(kmax>0&&(deltz-pblhdata[0][j][i][l])>=(pblhdata[0][j][i][l]-hgtdata[k-1][j][i][l]))
							kmax--;
						
						float tmp=-Cd*(float)Math.hypot(ut10data[0][j][i][l],vr10data[0][j][i][l])*
							cL(Tdata[tags[l][j][i]][j][i][l])*(q2mdata[0][j][i][l]-qdata[tags[l][j][i]][j][i][l])/pblhdata[0][j][i][l];
						
						for(int k2=tags[l][j][i];k2<=kmax;k2++) lhfldata[k2][j][i][l]=tmp;
						
						break;
					}
				}
			}
		}
		
		return lhfl;
	}
	
	
	/**
     * calculate force due to diabatic heating
     *
     * @param	Qm	tangentially averaged heating ratio
     *
     * @return	f	force
     */
	public Variable cDiabaticHeatingForce(Variable Qm){
		t=Qm.getTCount();	z=Qm.getZCount();	y=Qm.getYCount();	x=Qm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("dhf",Qm);	f.setCommentAndUnit("forcing of diabatic heating");
		
		float undef=Qm.getUndef();
		
		float[][][][] fdata= f.getData();
		float[][][][] Qdata=Qm.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++){
					if(Qdata[l][k][j+1][0]!=undef&&Qdata[l][k][j-1][0]!=undef)
						fdata[l][k][j][0]=(Qdata[l][k][j+1][0]-Qdata[l][k][j-1][0])/(2*dy)*Rd/Cp/zdef[k];
						
					else fdata[l][k][j][0]=undef;
				}
				
				/*** j==0 ***/
				if(Qdata[l][k][0][0]!=undef&&Qdata[l][k][1][0]!=undef)
					fdata[l][k][0][0]=(Qdata[l][k][1][0]-Qdata[l][k][0][0])/dy*Rd/Cp/zdef[k];
					
				else fdata[l][k][0][0]=undef;
				
				/*** j==y-1 ***/
				if(Qdata[l][k][y-1][0]!=undef&&Qdata[l][k][y-2][0]!=undef)
					fdata[l][k][y-1][0]=(Qdata[l][k][y-1][0]-Qdata[l][k][y-2][0])/dy*Rd/Cp/zdef[k];
					
				else fdata[l][k][y-1][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++){
					if(Qdata[k][j+1][0][l]!=undef&&Qdata[k][j-1][0][l]!=undef)
						fdata[k][j][0][l]=(Qdata[k][j+1][0][l]-Qdata[k][j-1][0][l])/(2*dy)*Rd/Cp/zdef[k];
						
					else fdata[k][j][0][l]=undef;
				}
				
				/*** j==0 ***/
				if(Qdata[k][0][0][l]!=undef&&Qdata[k][1][0][l]!=undef)
					fdata[k][0][0][l]=(Qdata[k][1][0][l]-Qdata[k][0][0][l])/dy*Rd/Cp/zdef[k];
					
				else fdata[k][0][0][l]=undef;
				
				/*** j==y-1 ***/
				if(Qdata[k][y-1][0][l]!=undef&&Qdata[k][y-2][0][l]!=undef)
					fdata[k][y-1][0][l]=(Qdata[k][y-1][0][l]-Qdata[k][y-2][0][l])/dy*Rd/Cp/zdef[k];
					
				else fdata[k][y-1][0][l]=undef;
			}
		}
		
		return f;
	}
	
	/**
     * calculate force due to eddy heat advection
     *
     * @param	va	radial wind anomaly
     * @param	aa	specific volume anomaly
     *
     * @return	f	force
     */
	public Variable cEddyHeatAdvectionForce(Variable va,Variable aa){
		if(!va.isLike(aa)) throw new IllegalArgumentException("dimensions not same");
		
		int count=0;		float undef=va.getUndef();
		t=va.getTCount();	z=va.getZCount();	y=va.getYCount();	x=va.getXCount();
		
		Variable f=new Variable("edyhaf",va.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing of eddy heat advection");
		
		float[]		  buf  =new float[y];
		float[][][][] fdata= f.getData();
		float[][][][] vdata=va.getData();
		float[][][][] adata=aa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++){
					for(int i=0;i<x;i++)
						if(adata[l][k][j+1][i]!=undef&&adata[l][k][j-1][i]!=undef&&vdata[l][k][j][i]!=undef){
							count++;
							fdata[l][k][j][0]+=(adata[l][k][j+1][i]-adata[l][k][j-1][i])*vdata[l][k][j][i];
						}
					
					if(count!=0) fdata[l][k][j][0]/=(count*dy*2);
					else fdata[l][k][j][0]=undef;
					
					count=0;
				}
				
				/*** first ring ***/
				for(int i=0;i<x;i++)
					if(adata[l][k][0][i]!=undef&&adata[l][k][1][i]!=undef&&vdata[l][k][0][i]!=undef){
						count++;
						fdata[l][k][0][0]+=(adata[l][k][1][i]-adata[l][k][0][i])*vdata[l][k][0][i];
					}
				
				if(count!=0) fdata[l][k][0][0]/=(count*dy);
				else fdata[l][k][0][0]=undef;
				
				count=0;
				
				/*** last ring ***/
				for(int i=0;i<x;i++)
					if(adata[l][k][y-1][i]!=undef&&adata[l][k][y-2][i]!=undef&&vdata[l][k][y-1][i]!=undef){
						count++;
						fdata[l][k][y-1][0]+=(adata[l][k][y-1][i]-adata[l][k][y-2][i])*vdata[l][k][y-1][i];
					}
				
				if(count!=0) fdata[l][k][y-1][0]/=(count*dy);
				else fdata[l][k][y-1][0]=undef;
				
				count=0;
				
				/*** partial Beta again ***/
				for(int j=0;j<y;j++) buf[j]=fdata[l][k][j][0];
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef) fdata[l][k][j  ][0]=-(buf[j+1]-buf[j-1])/(2*dy);
				else fdata[l][k][j  ][0]=undef;
				
				/*** deal with the first ring ***/
				if(buf[1  ]!=undef&&buf[0  ]!=undef) fdata[l][k][0  ][0]=-(buf[1]-buf[0])*2/dy;
				else fdata[l][k][0  ][0]=undef;
				
				/*** deal with the last ring ***/
				if(buf[y-1]!=undef&&buf[y-2]!=undef) fdata[l][k][y-1][0]=-(buf[y-1]-buf[y-2])*2/dy;
				else fdata[l][k][y-1][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++){
					for(int i=0;i<x;i++)
						if(adata[k][j+1][i][l]!=undef&&adata[k][j-1][i][l]!=undef&&vdata[k][j][i][l]!=undef){
							count++;
							fdata[k][j][0][l]+=(adata[k][j+1][i][l]-adata[k][j-1][i][l])*vdata[k][j][i][l];
						}
					
					if(count!=0) fdata[k][j][0][l]/=(2*dy*count);
					else fdata[k][j][0][l]=undef;
					
					count=0;
				}
				
				/*** deal with the first ring ***/
				for(int i=0;i<x;i++)
					if(adata[k][y-1][i][l]!=undef&&adata[k][y-2][i][l]!=undef&&vdata[k][y-1][i][l]!=undef){
						count++;
						fdata[k][y-1][0][l]+=(adata[k][y-1][i][l]-adata[k][y-2][i][l])*vdata[k][y-1][i][l];
					}
				
				if(count!=0) fdata[k][y-1][0][l]/=dy*count;
				else fdata[k][y-1][0][l]=undef;
				
				count=0;
				
				/*** deal with the last ring ***/
				for(int i=0;i<x;i++)
					if(adata[k][0][i][l]!=undef&&adata[k][1][i][l]!=undef&&vdata[k][0][i][l]!=undef){
						count++;
						fdata[k][0][0][l]+=(adata[k][1][i][l]-adata[k][0][i][l])*vdata[k][0][i][l];
					}
				
				if(count!=0) fdata[k][0][0][l]/=dy*count;
				else fdata[k][0][0][l]=undef;
				
				count=0;
				
				/*** partial Beta again ***/
				for(int j=0;j<y;j++) buf[j]=fdata[k][j][0][l];
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef) fdata[k][j][0][l]=-(buf[j+1]-buf[j-1])/(2*dy);
				else fdata[k][j][0][l]=undef;
				
				/*** deal with the first ring ***/
				if(buf[1]!=undef&&buf[0]!=undef) fdata[k][0][0][l]=-(buf[1]-buf[0])*2/dy;
				else fdata[k][0][0][l]=undef;
				
				/*** deal with the last ring ***/
				if(buf[y-1]!=undef&&buf[y-2]!=undef) fdata[k][y-1][0][l]=-(buf[y-1]-buf[y-2])*2/dy;
				else fdata[k][y-1][0][l]=undef;
			}
		}
		
		return f;
	}
	
	/**
     * calculate force due to adiabatic heating
     *
     * @param	segmaa	static stability anomally
     * @param	omegaa	vertical velocity anomally
     *
     * @return	f	force
     */
	public Variable cAdiabaticHeatingForce(Variable segmaa,Variable omegaa){
		if(!segmaa.isLike(omegaa)) throw new IllegalArgumentException("dimensions not same");
		
		float undef=segmaa.getUndef();
		t=segmaa.getTCount();	z=segmaa.getZCount();	y=segmaa.getYCount();	x=segmaa.getXCount();
		
		Variable f=new Variable("adhf",segmaa.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing of adiabatic heating");
		
		float[]		  		buf  =new float[y];
		float[][][][] 	  fdata=f.getData();
		float[][][][] segmadata=segmaa.getData();
		float[][][][] omegadata=omegaa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
						if(segmadata[l][k][j][i]!=undef&&omegadata[l][k][j][i]!=undef){
							count++;
							fdata[l][k][j][0]+=segmadata[l][k][j][i]*omegadata[l][k][j][i];
						}
					
					if(count!=0) fdata[l][k][j][0]/=count;
					else fdata[l][k][j][0]=undef;
				}
				
				/*** partial Beta ***/
				for(int j=0;j<y;j++) buf[j]=fdata[l][k][j][0];
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef) fdata[l][k][j][0]=(buf[j+1]-buf[j-1])/(2*dy);
				else fdata[l][k][j][0]=undef;
				
				if(buf[1]!=undef&&buf[0]!=undef) fdata[l][k][0][0]=(buf[1]-buf[0])/dy;
				else fdata[l][k][0  ][0]=undef;
				if(buf[y-1]!=undef&&buf[y-2]!=undef) fdata[l][k][y-1][0]=(buf[y-1]-buf[y-2])/dy;
				else fdata[l][k][y-1][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
						if(segmadata[k][j][i][l]!=undef&&omegadata[k][j][i][l]!=undef){
							count++;
							fdata[k][j][0][l]+=segmadata[k][j][i][l]*omegadata[k][j][i][l];
						}
					
					if(count!=0) fdata[k][j][0][l]/=count;
					else fdata[k][j][0][l]=undef;
				}
				
				/*** partial Beta ***/
				for(int j=0;j<y;j++) buf[j]=fdata[k][j][0][l];
				
				for(int j=1;j<y-1;j++)
				if(buf[j+1]!=undef&&buf[j-1]!=undef) fdata[k][j][0][l]=(buf[j+1]-buf[j-1])/(2*dy);
				else fdata[k][j][0][l]=undef;
				
				if(buf[1]!=undef&&buf[0]!=undef) fdata[k][0][0][l]=(buf[1]-buf[0])/dy;
				else fdata[k][0  ][0][l]=undef;
				if(buf[y-1]!=undef&&buf[y-2]!=undef) fdata[k][y-1][0][l]=(buf[y-1]-buf[y-2])/dy;
				else fdata[k][y-1][0][l]=undef;
			}
		}
		
		return f;
	}
	
	/**
     * calculate force due to eddy angular momentum horizontal transport
     *
     * @param	gazm	tangential averaged angular momentum
     * @param	gaza	angular momentum anomaly
     * @param	va		vwind anomaly
     *
     * @return	f		force
     *
     * @exception		if gazm,gaza,va are not dimensionally the same
     */
	public Variable cEddyAngMomenHorizontalTransportForce(Variable gazm,Variable gaza,Variable va){
		t=gaza.getTCount();	z=gaza.getZCount();	y=gaza.getYCount();	x=gaza.getXCount();
		
		if(gazm.getXCount()!=1)
			throw new IllegalArgumentException("x-direction of 1st param should be averaged");
		if(!gaza.isLike(va))
			throw new IllegalArgumentException("dimensions not same");
		if(gazm.getYCount()!=va.getYCount()||gazm.getZCount()!=va.getZCount()||gazm.getTCount()!=va.getTCount())
			throw new IllegalArgumentException("dimensions not same");
		
		float undef=gazm.getUndef();
		
		Variable f=new Variable("edyadvf",gazm.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing of eddy angular momentum advection");
		
		float[][]	  	buf =new float[z][y];
		float[][][][]  fdata= f.getData();
		float[][][][] gmdata=gazm.getData();
		float[][][][] gadata=gaza.getData();
		float[][][][] vadata=va.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(gadata[l][k][j][i]!=undef&&vadata[l][k][j][i]!=undef){
						buf[k][j]+=gadata[l][k][j][i]*vadata[l][k][j][i];	count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					else buf[k][j]=undef;
				}
				
				/*** calculate one level ***/
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++)
				if(buf[k][j]!=undef&&buf[k][j+1]!=undef&&buf[k][j-1]!=undef&&gmdata[l][k][j][0]!=undef)
					fdata[l][k][j][0]=gmdata[l][k][j][0]*(buf[k][j+1]*bsin[j+1]-buf[k][j-1]*bsin[j-1])/dy
					/(REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
					
				else fdata[l][k][j][0]=undef;
				
				for(int k=0;k<z;k++){
					fdata[l][k][0][0]=undef;
					
					if(buf[k][y-1]!=undef&&buf[k][y-2]!=undef&&gmdata[l][k][y-1][0]!=undef)
						fdata[l][k][y-1][0]=2*gmdata[l][k][y-1][0]*(buf[k][y-1]*bsin[y-1]-buf[k][y-2]*bsin[y-2])/dy
						/(REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else fdata[l][k][y-1][0]=undef;
				}
				
				/*** stored to buf ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=fdata[l][k][j][0];
				
				/*** partial p ***/
				for(int k=1;k<z-1;k++)
				for(int j=1;j<y;j++)
					if(buf[k+1][j]!=undef&&buf[k-1][j]!=undef) fdata[l][k][j][0]=-(buf[k+1][j]-buf[k-1][j])/(2*dz);
					else fdata[l][k][j][0]=undef;
				
				for(int j=1;j<y;j++){
					if(buf[0][j]!=undef&&buf[1][j]!=undef) fdata[l][0][j][0]=-(buf[1][j]-buf[0][j])/dz;
					else fdata[l][0][j][0]=undef;
					
					if(buf[z-1][j]!=undef&&buf[z-2][j]!=undef) fdata[l][z-1][j][0]=-(buf[z-1][j]-buf[z-2][j])/dz;
					else fdata[l][z-1][j][0]=undef;
				}
				
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
					if(gadata[k][j][i][l]!=undef&&vadata[k][j][i][l]!=undef){
						buf[k][j]+=gadata[k][j][i][l]*vadata[k][j][i][l];	count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					else buf[k][j]=undef;
				}
				
				/*** calculate one level ***/
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++)
				if(buf[k][j]!=undef&&buf[k][j+1]!=undef&&buf[k][j-1]!=undef&&gmdata[k][j][0][l]!=undef)
					fdata[k][j][0][l]=gmdata[k][j][0][l]*(buf[k][j+1]*bsin[j+1]-buf[k][j-1]*bsin[j-1])/dy
					/(REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
					
				else fdata[k][j][0][l]=undef;
				
				for(int k=0;k<z;k++){
					fdata[k][0][0][l]=undef;
					
					if(buf[k][y-1]!=undef&&buf[k][y-2]!=undef&&gmdata[k][y-1][0][l]!=undef)
						fdata[k][y-1][0][l]=2*gmdata[k][y-1][0][l]*(buf[k][y-1]*bsin[y-1]-buf[k][y-2]*bsin[y-2])/dy
						/(REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else fdata[k][y-1][0][l]=undef;
				}
				
				/*** stored to buf ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=fdata[k][j][0][l];
				
				/*** partial p ***/
				for(int k=1;k<z-1;k++)
				for(int j=1;j<y;j++)
					if(buf[k+1][j]!=undef&&buf[k-1][j]!=undef) fdata[k][j][0][l]=-(buf[k+1][j]-buf[k-1][j])/(2*dz);
					else fdata[k][j][0][l]=undef;
				
				for(int j=1;j<y;j++){
					if(buf[0][j]!=undef&&buf[1][j]!=undef) fdata[0][j][0][l]=-(buf[1][j]-buf[0][j])/dz;
					else fdata[0][j][0][l]=undef;
					
					if(buf[z-1][j]!=undef&&buf[z-2][j]!=undef) fdata[z-1][j][0][l]=-(buf[z-1][j]-buf[z-2][j])/dz;
					else fdata[z-1][j][0][l]=undef;
				}
				
				/*** reset buf ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=0;
			}
		}
		
		return f;
	}
	
	/**
     * calculate force due to eddy angular momentum vertical transport
     *
     * @param	gazm		tangential averaged angular momentum
     * @param	gaza		angular momentum anomaly
     * @param	omegaa	    vertical velocity anomaly
     *
     * @return	f			force
     *
     * @exception			if um,ua,omegaa are not dimensionally the same
     */
	public Variable cEddyAngMomenVerticalTransportForce(Variable gazm,Variable gaza,Variable omegaa){
		t=gaza.getTCount();	z=gaza.getZCount();	y=gaza.getYCount();	x=gaza.getXCount();
		
		if(gazm.getXCount()!=1)
			throw new IllegalArgumentException("x-direction of 1st param should be averaged");
		if(!gaza.isLike(omegaa))
			throw new IllegalArgumentException("dimensions not same");
		if(gazm.getYCount()!=gaza.getYCount()||gazm.getZCount()!=gaza.getZCount()||gazm.getTCount()!=gaza.getTCount())
			throw new IllegalArgumentException("dimensions not same");
		
		float undef=gazm.getUndef();
		
		Variable f=null;	f=new Variable("edyconf",gazm.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing of eddy angular momentum convection");
		
		float[][]	  	buf  =new float[z][y];
		float[][][][]  fdata = f.getData();
		float[][][][] gmdata=gazm.getData();
		float[][][][] gadata=gaza.getData();
		float[][][][] oadata=omegaa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(gadata[l][k][j][i]!=undef&&oadata[l][k][j][i]!=undef){
						buf[k][j]+=gadata[l][k][j][i]*oadata[l][k][j][i];	count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					else buf[k][j]=undef;
				}
				
				/*** partial p ***/
				for(int k=1;k<z-1;k++){
					fdata[l][k][0][0]=undef;
					
					for(int j=1;j<y;j++){
						if(gmdata[l][k+1][j][0]!=undef&&gmdata[l][k][j][0]!=undef&&gmdata[l][k-1][j][0]!=undef&&
						buf[k+1][j]!=undef&&buf[k][j]!=undef&&buf[k-1][j]!=undef)
							fdata[l][k][j][0]=-2*(
								(gmdata[l][k+1][j][0]+gmdata[l][k][j][0])*(buf[k+1][j]-buf[k][j])-
								(gmdata[l][k-1][j][0]+gmdata[l][k][j][0])*(buf[k][j]-buf[k-1][j])
							)/(float)pow(bsin[j]*REarth,3)/(dz*dz*2);
						
						else fdata[l][k][j][0]=undef;
					}
				}
				
				fdata[l][0][0][0]=undef;	fdata[l][z-1][0][0]=undef;
				
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
					if(gadata[k][j][i][l]!=undef&&oadata[k][j][i][l]!=undef){
						buf[k][j]+=gadata[k][j][i][l]*oadata[k][j][i][l];	count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					else buf[k][j]=undef;
				}
				
				/*** partial p ***/
				for(int k=1;k<z-1;k++){
					fdata[k][0][0][l]=undef;
					
					for(int j=1;j<y;j++){
						if(gmdata[k+1][j][0][l]!=undef&&gmdata[k][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&
						buf[k+1][j]!=undef&&buf[k][j]!=undef&&buf[k-1][j]!=undef)
							fdata[k][j][0][l]=-2*(
								(gmdata[k+1][j][0][l]+gmdata[k][j][0][l])*(buf[k+1][j]-buf[k][j])-
								(gmdata[k-1][j][0][l]+gmdata[k][j][0][l])*(buf[k][j]-buf[k-1][j])
							)/(float)pow(bsin[j]*REarth,3)/(dz*dz*2);
							
						else fdata[k][j][0][l]=undef;
					}
				}
				
				fdata[0][0][0][l]=undef;	fdata[z-1][0][0][l]=undef;
				
				/*** reset buf ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=0;
			}
		}
		
		return f;
	}
	
	/**
     * calculate frictional torque
     *
     * @param	gazm	tangential averaged angular momentum
     * @param	frm		tangential averaged friction
     *
     * @return	f	force
     */
	public Variable cFrictionalTorqueForce(Variable frm,Variable gazm){
		t=frm.getTCount();	z=frm.getZCount();	y=frm.getYCount();	x=frm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		if(!frm.isLike(gazm)) throw new IllegalArgumentException("dimensions not same");
		
		float undef=frm.getUndef();
		
		Variable ft=new Variable("frictf",frm);
		ft.setCommentAndUnit("forcing of frictional torque");
		
		float[][][][] ftdata=  ft.getData();
		float[][][][] frdata= frm.getData();
		float[][][][] gmdata=gazm.getData();
		
		if(frm.isTFirst()){
			for(int j=1;j<y-1;j++){
				float tmp=(float)Math.pow(REarth*bsin[j],2);
				
				for(int l=0;l<t;l++)
				for(int k=1;k<z-1;k++)
				if(frdata[l][k-1][j][0]!=undef&&frdata[l][k+1][j][0]!=undef&&gmdata[l][k][j][0]!=undef)
					ftdata[l][k][j][0]=
					2*(frdata[l][k+1][j][0]*gmdata[l][k+1][j][0]-frdata[l][k-1][j][0]*gmdata[l][k-1][j][0])/(2*dz)/tmp;
			}
			
		}else{
			for(int j=1;j<y-1;j++){
				float tmp=(float)Math.pow(REarth*bsin[j],2);
				
				for(int l=0;l<t;l++)
				for(int k=1;k<z-1;k++)
				if(frdata[k-1][j][0][l]!=undef&&frdata[k+1][j][0][l]!=undef&&gmdata[k][j][0][l]!=undef)
					ftdata[k][j][0][l]=
					2*(frdata[k+1][j][0][l]*gmdata[k+1][j][0][l]-frdata[k-1][j][0][l]*gmdata[k-1][j][0][l])/(2*dz)/tmp;
			}
		}
		
		return ft;
	}
	
	/**
     * calculate tilting force
     *
     * @param	ut		tangential wind component
     * @param	vr		radial wind component
     * @param	gazm	tangential averaged angular momentum
     *
     * @return	f	force
     */
	public Variable cTiltingForce(Variable gazm,Variable ut,Variable vr){
		t=ut.getTCount();	z=ut.getZCount();	y=ut.getYCount();	x=ut.getXCount();
		
		if(gazm.getXCount()!=1)
			throw new IllegalArgumentException("x-direction of 1st param should be averaged");
		if(!ut.isLike(vr))
			throw new IllegalArgumentException("dimensions not same");
		if(gazm.getYCount()!=ut.getYCount()||gazm.getZCount()!=ut.getZCount()||gazm.getTCount()!=ut.getTCount())
			throw new IllegalArgumentException("dimensions not same");
		
		int count=0;	float undef=gazm.getUndef();
		
		Variable f=new Variable("tiltf",gazm.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing of tilting effect");
		
		float til2,til3,til4,til5;
    	float[] olat=((CylindricalSpatialModel)sm).getOLat();
    	float[][] buf=new float[z][y];
		float[][][][]  fdata = f.getData();
		float[][][][] gmdata=gazm.getData();
		float[][][][] utdata=ut.getData();
		float[][][][] vrdata=vr.getData();

		if(f.isTFirst()){
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					for(int i=0;i<x;i++)
					if(utdata[l][k][j][i]!=undef&&vrdata[l][k][j][i]!=undef){
						til3=(1-(float)cos(ydef[j]))
							*(wtan[l][i]*vrdata[l][k][j][i]+wnor[l][i]*utdata[l][k][j][i]);
						
						til4=-vrdata[l][k][j][i]*bsin[j]*(float)(
							cos(ydef[j])*sin(olat[l])+bsin[j]*cos(olat[l])*cos(xdef[i])-sin(olat[l])
						);
						
						til5=(1-(float)cos(ydef[j]))*(float)(
							(utdata[l][k][j][i]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[l][k][j][i]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i])
						);
						
						buf[k][j]+=til3+REarth*omegaEarth*(til4+til5);
						
						count++;
						
					}else buf[k][j]=undef;
					
					if(count!=0) buf[k][j]/=count;
					
					til2=-REarth*omegaEarth*cs[l]*bsin[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2;
					buf[k][j]+=til2;
					count=0;
				}

				/*** partial p ***/
				for(int j=1;j<y;j++)
				for(int k=1;k<z-1;k++)
				if(gmdata[l][k+1][j][0]!=undef&&gmdata[l][k-1][j][0]!=undef&&buf[k+1][j]!=undef&&buf[k-1][j]!=undef)
					fdata[l][k][j][0]=
					2*(gmdata[l][k+1][j][0]*buf[k+1][j]-gmdata[l][k-1][j][0]*buf[k-1][j])
					/(float)pow(bsin[j]*REarth,3)/(dz*2);
				
				else fdata[k][j][0][l]=undef;
				
				/*** the first level ***/
				for(int j=1;j<y;j++)
				if(gmdata[l][1][j][0]!=undef&&gmdata[l][0][j][0]!=undef&&buf[1][j]!=undef&&buf[0][j]!=undef)
					fdata[l][0][j][0]=
					2*(gmdata[l][1][j][0]*buf[1][j]-gmdata[l][0][j][0]*buf[0][j])
					/(float)pow(bsin[j]*REarth,3)/dz;
				
				else fdata[l][0][j][0]=undef;
				
				/*** the last level ***/
				for(int j=1;j<y;j++)
				if(gmdata[l][z-1][j][0]!=undef&&gmdata[l][z-2][j][0]!=undef&&buf[z-1][j]!=undef&&buf[z-2][j]!=undef)
					fdata[l][z-1][j][0]=
					2*(gmdata[l][z-1][j][0]*buf[z-1][j]-gmdata[l][z-2][j][0]*buf[z-2][j])
					/(float)pow(bsin[j]*REarth,3)/dz;
				
				else fdata[l][z-1][j][0]=undef;
				
				/*** reset buf ***/
				for(int k=0;k<z;k++){
					for(int j=0;j<y;j++) buf[k][j]=0;
					
					fdata[l][k][0][0]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					for(int i=0;i<x;i++)
					if(utdata[k][j][i][l]!=undef&&vrdata[k][j][i][l]!=undef){
						til3=(1-(float)cos(ydef[j]))
							*(wtan[l][i]*vrdata[k][j][i][l]+wnor[l][i]*utdata[k][j][i][l]);
						
						til4=-vrdata[k][j][i][l]*bsin[j]
						    *(float)(cos(ydef[j])*sin(olat[l])
						    +bsin[j]*cos(olat[l])*cos(xdef[i])
						    -sin(olat[l]));
						
						til5=(1-(float)cos(ydef[j]))*(float)(
							(utdata[k][j][i][l]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[k][j][i][l]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i]));
						
						buf[k][j]+=til3+REarth*omegaEarth*(til4+til5);
						
						count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					
					til2=-REarth*omegaEarth*cs[l]*bsin[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2;
					buf[k][j]+=til2;
					count=0;
				}

				/*** partial p ***/
				for(int j=1;j<y;j++)
				for(int k=1;k<z-1;k++)
				if(gmdata[k+1][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&buf[k+1][j]!=undef&&buf[k-1][j]!=undef)
					fdata[k][j][0][l]=
					2*(gmdata[k+1][j][0][l]*buf[k+1][j]-gmdata[k-1][j][0][l]*buf[k-1][j])
					/(float)pow(bsin[j]*REarth,3)/(dz*2);
				
				else fdata[k][j][0][l]=undef;
				
				/*** the first level ***/
				for(int j=1;j<y;j++)
				if(gmdata[1][j][0][l]!=undef&&gmdata[0][j][0][l]!=undef&&buf[1][j]!=undef&&buf[0][j]!=undef)
					fdata[0][j][0][l]=
					2*(gmdata[1][j][0][l]*buf[1][j]-gmdata[0][j][0][l]*buf[0][j])
					/(float)pow(bsin[j]*REarth,3)/dz;
				
				else fdata[0][j][0][l]=undef;
				
				/*** the last level ***/
				for(int j=1;j<y;j++)
				if(gmdata[z-1][j][0][l]!=undef&&gmdata[z-2][j][0][l]!=undef&&buf[z-1][j]!=undef&&buf[z-2][j]!=undef)
					fdata[z-1][j][0][l]=
					2*(gmdata[z-1][j][0][l]*buf[z-1][j]-gmdata[z-2][j][0][l]*buf[z-2][j])
					/(float)pow(bsin[j]*REarth,3)/dz;
				
				else fdata[z-1][j][0][l]=undef;
				
				/*** reset buf ***/
				for(int k=0;k<z;k++){
					for(int j=0;j<y;j++) buf[k][j]=0;
					
					fdata[k][0][0][l]=undef;
				}
			}
		}
		
		return f;
	}
	
	
	/**
     * calculate force due to diabatic heating
     *
     * @param	Qm	tangentially averaged heating ratio
     *
     * @return	f	force
     */
	public Variable cDiabaticHeatingFF(Variable Qm){
		t=Qm.getTCount();	z=Qm.getZCount();	y=Qm.getYCount();	x=Qm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable f=new Variable("dhFF",Qm);	f.setCommentAndUnit("forcing factor of diabatic heating");
		
		float[][][][] fdata= f.getData();
		float[][][][] Qdata=Qm.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) fdata[l][k][j][0]=Qdata[l][k][j][0]*Rd/Cp/zdef[k];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++) fdata[k][j][0][l]=Qdata[k][j][0][l]*Rd/Cp/zdef[k];
		}
		
		return f;
	}
	
	/**
     * calculate force due to eddy heat advection
     *
     * @param	va	radial wind anomaly
     * @param	aa	specific volume anomaly
     *
     * @return	f	force
     */
	public Variable cEddyHeatAdvectionFF(Variable va,Variable aa){
		if(!va.isLike(aa)) throw new IllegalArgumentException("dimensions not same");
		
		int count=0;		float undef=va.getUndef();
		t=va.getTCount();	z=va.getZCount();	y=va.getYCount();	x=va.getXCount();
		
		Variable f=new Variable("edyhaFF",va.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing factor of eddy heat advection");
		
		float[][][][] fdata= f.getData();
		float[][][][] vdata=va.getData();
		float[][][][] adata=aa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++){
					for(int i=0;i<x;i++)
						if(adata[l][k][j+1][i]!=undef&&adata[l][k][j-1][i]!=undef&&vdata[l][k][j][i]!=undef){
							count++;
							fdata[l][k][j][0]+=(adata[l][k][j+1][i]-adata[l][k][j-1][i])*vdata[l][k][j][i];
						}
					
					if(count!=0) fdata[l][k][j][0]/=-(count*dy*2);
					else fdata[l][k][j][0]=undef;
					
					count=0;
				}
				
				/*** first ring ***/
				for(int i=0;i<x;i++)
					if(adata[l][k][0][i]!=undef&&adata[l][k][1][i]!=undef&&vdata[l][k][0][i]!=undef){
						count++;
						fdata[l][k][0][0]+=(adata[l][k][1][i]-adata[l][k][0][i])*vdata[l][k][0][i];
					}
				
				if(count!=0) fdata[l][k][0][0]/=-(count*dy);
				else fdata[l][k][0][0]=undef;
				
				count=0;
				
				/*** last ring ***/
				for(int i=0;i<x;i++)
					if(adata[l][k][y-1][i]!=undef&&adata[l][k][y-2][i]!=undef&&vdata[l][k][y-1][i]!=undef){
						count++;
						fdata[l][k][y-1][0]+=(adata[l][k][y-1][i]-adata[l][k][y-2][i])*vdata[l][k][y-1][i];
					}
				
				if(count!=0) fdata[l][k][y-1][0]/=-(count*dy);
				else fdata[l][k][y-1][0]=undef;
				
				count=0;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++){
					for(int i=0;i<x;i++)
						if(adata[k][j+1][i][l]!=undef&&adata[k][j-1][i][l]!=undef&&vdata[k][j][i][l]!=undef){
							count++;
							fdata[k][j][0][l]+=(adata[k][j+1][i][l]-adata[k][j-1][i][l])*vdata[k][j][i][l];
						}
					
					if(count!=0) fdata[k][j][0][l]/=-(2*dy*count);
					else fdata[k][j][0][l]=undef;
					
					count=0;
				}
				
				/*** deal with the first ring ***/
				for(int i=0;i<x;i++)
					if(adata[k][y-1][i][l]!=undef&&adata[k][y-2][i][l]!=undef&&vdata[k][y-1][i][l]!=undef){
						count++;
						fdata[k][y-1][0][l]+=(adata[k][y-1][i][l]-adata[k][y-2][i][l])*vdata[k][y-1][i][l];
					}
				
				if(count!=0) fdata[k][y-1][0][l]/=-dy*count;
				else fdata[k][y-1][0][l]=undef;
				
				count=0;
				
				/*** deal with the last ring ***/
				for(int i=0;i<x;i++)
					if(adata[k][0][i][l]!=undef&&adata[k][1][i][l]!=undef&&vdata[k][0][i][l]!=undef){
						count++;
						fdata[k][0][0][l]+=(adata[k][1][i][l]-adata[k][0][i][l])*vdata[k][0][i][l];
					}
				
				if(count!=0) fdata[k][0][0][l]/=-dy*count;
				else fdata[k][0][0][l]=undef;
				
				count=0;
			}
		}
		
		return f;
	}
	
	/**
     * calculate force due to adiabatic heating
     *
     * @param	segmaa	static stability anomally
     * @param	omegaa	vertical velocity anomally
     *
     * @return	f	force
     */
	public Variable cAdiabaticHeatingFF(Variable sigmaa,Variable omegaa){
		if(!sigmaa.isLike(omegaa)) throw new IllegalArgumentException("dimensions not same");
		
		float undef=sigmaa.getUndef();
		t=sigmaa.getTCount();	z=sigmaa.getZCount();	y=sigmaa.getYCount();	x=sigmaa.getXCount();
		
		Variable f=new Variable("adhFF",sigmaa.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing factor of adiabatic heating");
		
		float[][][][] 	  fdata=f.getData();
		float[][][][] sigmadata=sigmaa.getData();
		float[][][][] omegadata=omegaa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(sigmadata[l][k][j][i]!=undef&&omegadata[l][k][j][i]!=undef){
						count++;
						fdata[l][k][j][0]+=sigmadata[l][k][j][i]*omegadata[l][k][j][i];
					}
					
					if(count!=0) fdata[l][k][j][0]/=count;
					else fdata[l][k][j][0]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(sigmadata[k][j][i][l]!=undef&&omegadata[k][j][i][l]!=undef){
						count++;
						fdata[k][j][0][l]+=sigmadata[k][j][i][l]*omegadata[k][j][i][l];
					}
					
					if(count!=0) fdata[k][j][0][l]/=count;
					else fdata[k][j][0][l]=undef;
				}
			}
		}
		
		return f;
	}
	
	/**
     * calculate force due to eddy angular momentum horizontal transport
     *
     * @param	gazm	tangential averaged angular momentum
     * @param	gaza	angular momentum anomaly
     * @param	va		vwind anomaly
     *
     * @return	f		force
     */
	public Variable cEddyAngMomenHorizontalTransportFF(Variable gazm,Variable gaza,Variable va){
		t=gaza.getTCount();	z=gaza.getZCount();	y=gaza.getYCount();	x=gaza.getXCount();
		
		if(gazm.getXCount()!=1)
			throw new IllegalArgumentException("x-direction of 1st param should be averaged");
		if(!gaza.isLike(va))
			throw new IllegalArgumentException("dimensions not same");
		if(gazm.getYCount()!=va.getYCount()||gazm.getZCount()!=va.getZCount()||gazm.getTCount()!=va.getTCount())
			throw new IllegalArgumentException("dimensions not same");
		
		float undef=gazm.getUndef();
		
		Variable f=new Variable("edyadvFF",gazm.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing factor of eddy angular momentum horizontal transport");
		
		float[][]	  	gava =new float[z][y];
		float[][][][]  fdata = f.getData();
		float[][][][] gmdata=gazm.getData();
		float[][][][] gadata=gaza.getData();
		float[][][][] vadata=va.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(gadata[l][k][j][i]!=undef&&vadata[l][k][j][i]!=undef){
						gava[k][j]+=gadata[l][k][j][i]*vadata[l][k][j][i];	count++;
					}
					
					if(count!=0) gava[k][j]/=count;
					else gava[k][j]=undef;
				}
				
				/*** calculate one level ***/
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++)
				if(gava[k][j]!=undef&&gava[k][j+1]!=undef&&gava[k][j-1]!=undef&&gmdata[l][k][j][0]!=undef)
					fdata[l][k][j][0]=-gmdata[l][k][j][0]*(gava[k][j+1]*bsin[j+1]-gava[k][j-1]*bsin[j-1])/dy
					/(REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
				
				else fdata[l][k][j][0]=undef;
				
				for(int k=0;k<z;k++){
					fdata[l][k][0][0]=undef;
					
					if(gava[k][y-1]!=undef&&gava[k][y-2]!=undef&&gmdata[l][k][y-1][0]!=undef)
						fdata[l][k][y-1][0]=-2*gmdata[l][k][y-1][0]*(gava[k][y-1]*bsin[y-1]-gava[k][y-2]*bsin[y-2])/dy
						/(REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else fdata[l][k][y-1][0]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(gadata[k][j][i][l]!=undef&&vadata[k][j][i][l]!=undef){
						gava[k][j]+=gadata[k][j][i][l]*vadata[k][j][i][l];	count++;
					}
					
					if(count!=0) gava[k][j]/=count;
					else gava[k][j]=undef;
				}
				
				/*** calculate one level ***/
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++)
				if(gava[k][j]!=undef&&gava[k][j+1]!=undef&&gava[k][j-1]!=undef&&gmdata[k][j][0][l]!=undef)
					fdata[k][j][0][l]=-gmdata[k][j][0][l]*(gava[k][j+1]*bsin[j+1]-gava[k][j-1]*bsin[j-1])/dy
					/(REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
				
				else fdata[k][j][0][l]=undef;
				
				for(int k=0;k<z;k++){
					fdata[k][0][0][l]=undef;
					
					if(gava[k][y-1]!=undef&&gava[k][y-2]!=undef&&gmdata[k][y-1][0][l]!=undef)
						fdata[k][y-1][0][l]=-2*gmdata[k][y-1][0][l]*(gava[k][y-1]*bsin[y-1]-gava[k][y-2]*bsin[y-2])/dy
						/(REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else fdata[k][y-1][0][l]=undef;
				}
			}
		}
		
		return f;
	}
	
	/**
     * calculate force due to eddy angular momentum vertical transport
     *
     * @param	gazm		tangential averaged angular momentum
     * @param	gaza		angular momentum anomaly
     * @param	omegaa	    vertical velocity anomaly
     *
     * @return	f			force
     *
     * @exception			if um,ua,omegaa are not dimensionally the same
     */
	public Variable cEddyAngMomenVerticalTransportFF(Variable gazm,Variable gaza,Variable omegaa){
		t=gaza.getTCount();	z=gaza.getZCount();	y=gaza.getYCount();	x=gaza.getXCount();
		
		if(gazm.getXCount()!=1)
			throw new IllegalArgumentException("x-direction of 1st param should be averaged");
		if(!gaza.isLike(omegaa))
			throw new IllegalArgumentException("dimensions not same");
		if(gazm.getYCount()!=gaza.getYCount()||gazm.getZCount()!=gaza.getZCount()||gazm.getTCount()!=gaza.getTCount())
			throw new IllegalArgumentException("dimensions not same");
		
		float undef=gazm.getUndef();
		
		Variable f=new Variable("edyconFF",gazm.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing factor of eddy angular momentum vertical transport");
		
		float[][]	  	buf =new float[z][y];
		float[][][][]  fdata= f.getData();
		float[][][][] gmdata=gazm.getData();
		float[][][][] gadata=gaza.getData();
		float[][][][] oadata=omegaa.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(gadata[l][k][j][i]!=undef&&oadata[l][k][j][i]!=undef){
						buf[k][j]+=gadata[l][k][j][i]*oadata[l][k][j][i];	count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					else buf[k][j]=undef;
				}
				
				/*** partial p ***/
				for(int k=1;k<z-1;k++){
					fdata[l][k][0][0]=undef;
					
					for(int j=1;j<y;j++){
						if(gmdata[l][k+1][j][0]!=undef&&gmdata[l][k][j][0]!=undef&&gmdata[l][k-1][j][0]!=undef&&
						buf[k+1][j]!=undef&&buf[k][j]!=undef&&buf[k-1][j]!=undef)
							fdata[l][k][j][0]=-2*(
								gmdata[l][k][j][0]*(buf[k+1][j]-buf[k-1][j])
							)/(float)pow(bsin[j]*REarth,3)/(dz*2);
						
						else fdata[l][k][j][0]=undef;
					}
				}
				
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
					if(gadata[k][j][i][l]!=undef&&oadata[k][j][i][l]!=undef){
						buf[k][j]+=gadata[k][j][i][l]*oadata[k][j][i][l];	count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					else buf[k][j]=undef;
				}
				
				/*** partial p ***/
				for(int k=1;k<z-1;k++){
					fdata[k][0][0][l]=undef;
					
					for(int j=1;j<y;j++){
						if(gmdata[k+1][j][0][l]!=undef&&gmdata[k][j][0][l]!=undef&&gmdata[k-1][j][0][l]!=undef&&
						buf[k+1][j]!=undef&&buf[k][j]!=undef&&buf[k-1][j]!=undef)
							fdata[k][j][0][l]=-2*(
								gmdata[k][j][0][l]*(buf[k+1][j]-buf[k-1][j])
							)/(float)pow(bsin[j]*REarth,3)/(dz*2);
							
						else fdata[k][j][0][l]=undef;
					}
				}
				
				/*** reset buf ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++) buf[k][j]=0;
			}
		}
		
		return f;
	}
	
	/**
     * calculate frictional torque factor
     *
     * @param	gazm	tangential averaged angular momentum
     * @param	frm		tangential averaged friction
     *
     * @return	f	force
     */
	public Variable cFrictionalTorqueFF(Variable frm,Variable gazm){
		t=frm.getTCount();	z=frm.getZCount();	y=frm.getYCount();	x=frm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		if(!frm.isLike(gazm)) throw new IllegalArgumentException("dimensions not same");
		
		float undef=frm.getUndef();
		
		Variable ft=new Variable("frictFF",frm);
		ft.setCommentAndUnit("forcing factor of frictional torque");
		
		float[][][][] ftdata=  ft.getData();
		float[][][][] frdata= frm.getData();
		float[][][][] gmdata=gazm.getData();
		
		if(frm.isTFirst()){
			for(int j=1;j<y-1;j++){
				float tmp=(float)Math.pow(REarth*bsin[j],2);
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				if(frdata[l][k][j][0]!=undef&&gmdata[l][k][j][0]!=undef)
				ftdata[l][k][j][0]=2*gmdata[l][k][j][0]*frdata[l][k][j][0]/tmp;
			}
			
		}else{
			for(int j=1;j<y-1;j++){
				float tmp=(float)Math.pow(REarth*bsin[j],2);
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				if(frdata[k][j][0][l]!=undef&&gmdata[k][j][0][l]!=undef)
				ftdata[k][j][0][l]=2*gmdata[k][j][0][l]*frdata[k][j][0][l]/tmp;
			}
		}
		
		return ft;
	}
	
	/**
     * calculate tilting force factor
     *
     * @param	ut		tangential wind component
     * @param	vr		radial wind component
     * @param	gazm	tangential averaged angular momentum
     *
     * @return	f	force
     */
	public Variable cTiltingFF(Variable gazm,Variable ut,Variable vr){
		t=ut.getTCount();	z=ut.getZCount();	y=ut.getYCount();	x=ut.getXCount();
		
		if(gazm.getXCount()!=1)
			throw new IllegalArgumentException("x-direction of 1st param should be averaged");
		if(!ut.isLike(vr))
			throw new IllegalArgumentException("dimensions not same");
		if(gazm.getYCount()!=ut.getYCount()||gazm.getZCount()!=ut.getZCount()||gazm.getTCount()!=ut.getTCount())
			throw new IllegalArgumentException("dimensions not same");
		
		int count=0;	float undef=gazm.getUndef();
		
		Variable f=new Variable("tiltFF",gazm.isTFirst(),new Range(t,z,y,1));
		f.setUndef(undef);	f.setCommentAndUnit("forcing factor of tilting effect");
		
		float til2,til3,til4,til5;
    	float[] olat=((CylindricalSpatialModel)sm).getOLat();
    	float[][] buf=new float[z][y];
		float[][][][]  fdata = f.getData();
		float[][][][] gmdata=gazm.getData();
		float[][][][] utdata=ut.getData();
		float[][][][] vrdata=vr.getData();
		
		if(f.isTFirst()){
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					for(int i=0;i<x;i++)
					if(utdata[l][k][j][i]!=undef&&vrdata[l][k][j][i]!=undef){
						til3=(1-(float)cos(ydef[j]))
							*(wtan[l][i]*vrdata[l][k][j][i]+wnor[l][i]*utdata[l][k][j][i]);
						
						til4=-vrdata[l][k][j][i]*bsin[j]*(float)(
							cos(ydef[j])*sin(olat[l])+bsin[j]*cos(olat[l])*cos(xdef[i])-sin(olat[l])
						);
						
						til5=(1-(float)cos(ydef[j]))*(float)(
							(utdata[l][k][j][i]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[l][k][j][i]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i])
						);
						
						buf[k][j]+=til3+REarth*omegaEarth*(til4+til5);
						
						count++;
						
					}else buf[k][j]=undef;
					
					if(count!=0) buf[k][j]/=count;
					
					til2=-REarth*omegaEarth*cs[l]*bsin[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2;
					buf[k][j]+=til2;
					count=0;
				}
				
				/*** assign ***/
				for(int j=1;j<y;j++)
				for(int k=0;k<z;k++)
				if(gmdata[l][k][j][0]!=undef&&buf[k][j]!=undef)
					fdata[l][k][j][0]=2*gmdata[l][k][j][0]*buf[k][j]/(float)pow(bsin[j]*REarth,3);
				else
					fdata[l][k][j][0]=undef;
				
				/*** reset buf ***/
				for(int k=0;k<z;k++){
					for(int j=0;j<y;j++) buf[k][j]=0;
					
					fdata[l][k][0][0]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				/*** calculate average ***/
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					for(int i=0;i<x;i++)
					if(utdata[k][j][i][l]!=undef&&vrdata[k][j][i][l]!=undef){
						til3=(1-(float)cos(ydef[j]))
							*(wtan[l][i]*vrdata[k][j][i][l]+wnor[l][i]*utdata[k][j][i][l]);
						
						til4=-vrdata[k][j][i][l]*bsin[j]
						    *(float)(cos(ydef[j])*sin(olat[l])
						    +bsin[j]*cos(olat[l])*cos(xdef[i])
						    -sin(olat[l]));
						
						til5=(1-(float)cos(ydef[j]))*(float)(
							(utdata[k][j][i][l]-wtan[l][i])*cos(olat[l])*sin(xdef[i])
							+vrdata[k][j][i][l]*(cos(ydef[j])*cos(olat[l])*cos(xdef[i])-bsin[j]*sin(olat[l]))
							-wnor[l][i]*cos(olat[l])*cos(xdef[i]));
						
						buf[k][j]+=til3+REarth*omegaEarth*(til4+til5);
						
						count++;
					}
					
					if(count!=0) buf[k][j]/=count;
					
					til2=-REarth*omegaEarth*cs[l]*bsin[j]*bsin[j]*(float)(cos(olat[l])*cos(cd[l]))/2;
					buf[k][j]+=til2;
					count=0;
				}

				/*** assign ***/
				for(int j=1;j<y;j++)
				for(int k=0;k<z;k++)
				if(gmdata[k][j][0][l]!=undef&&buf[k][j]!=undef)
					fdata[k][j][0][l]=2*gmdata[k][j][0][l]*buf[k][j]/(float)pow(bsin[j]*REarth,3);
				else
					fdata[k][j][0][l]=undef;
				
				/*** reset buf ***/
				for(int k=0;k<z;k++){
					for(int j=0;j<y;j++) buf[k][j]=0;
					
					fdata[k][0][0][l]=undef;
				}
			}
		}
		
		return f;
	}
	
	
	/**
     * implement the methods in EllipticalInterface
     */
	public Variable cAPrime(Variable sigmam){
		t=sigmam.getTCount();	z=sigmam.getZCount();	y=sigmam.getYCount();	x=sigmam.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		A=new Variable("Aa",sigmam);	float undef=sigmam.getUndef();
		A.setCommentAndUnit("coefficient A of elliptical equation");
		
		float[][][][] Adata=     A.getData();
		float[][][][] Sdata=sigmam.getData();
		
		if(sigmam.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y;j++)
					if(Sdata[l][k][j][0]!=undef)
						Adata[l][k][j][0]=
						(float)((Sdata[l][k][j-1][0]+Sdata[l][k][j][0])/2/sin((ydef[j-1]+ydef[j])/2));
					
					else Adata[l][k][j][0]=undef;
				
				Adata[l][k][0][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y;j++)
					if(Sdata[k][j][0][l]!=undef)
						Adata[k][j][0][l]=
						(float)((Sdata[k][j-1][0][l]+Sdata[k][j][0][l])/2/sin((ydef[j-1]+ydef[j])/2));
						
					else Adata[k][j][0][l]=undef;
				
				Adata[k][0][0][l]=undef;
			}
		}
		
		return A;
	}
	
	public Variable cBPrime(Variable am){
		t=am.getTCount();	z=am.getZCount();	y=am.getYCount();	x=am.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		B=new Variable("Bb",am);	float undef=am.getUndef();
		B.setCommentAndUnit("coefficient B of elliptical equation");
		
		float[][][][] Bdata= B.getData();
		float[][][][] adata=am.getData();
		
		if(am.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++){
					for(int j=1;j<y-1;j++)
					if(adata[l][k][j+1][0]!=undef&&adata[l][k][j-1][0]!=undef)
						Bdata[l][k][j][0]=(adata[l][k][j+1][0]-adata[l][k][j-1][0])/bsin[j]/(2*dy);
					
					else Bdata[l][k][j][0]=undef;
					
					if(adata[l][k][y-1][0]!=undef&&adata[l][k][y-2][0]!=undef)
						Bdata[l][k][y-1][0]=2*(adata[l][k][y-1][0]-adata[l][k][y-2][0])/
							(float)sin((ydef[y-2]+ydef[y-1])/2)/dy-Bdata[l][k][y-2][0];
					
					else Bdata[l][k][y-1][0]=undef;
					
					if(adata[l][k][0][0]!=undef&&adata[l][k][1][0]!=undef)
						Bdata[l][k][0][0]=2*(adata[l][k][1][0]-adata[l][k][0][0])/dy/(float)sin(ydef[1]/2)-Bdata[l][k][1][0];
					
					else Bdata[l][k][0][0]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++){
					for(int j=1;j<y-1;j++)
					if(adata[k][j+1][0][l]!=undef&&adata[k][j-1][0][l]!=undef)
						Bdata[k][j][0][l]=(adata[k][j+1][0][l]-adata[k][j-1][0][l])/bsin[j]/(2*dy);
					
					else Bdata[k][j][0][l]=undef;
					
					if(adata[k][y-1][0][l]!=undef&&adata[k][y-2][0][l]!=undef)
						Bdata[k][y-1][0][l]=2*(adata[k][y-1][0][l]-adata[k][y-2][0][l])/
							(float)sin((ydef[y-2]+ydef[y-1])/2)/dy-Bdata[k][y-2][0][l];
					
					else Bdata[k][y-1][0][l]=undef;
					
					if(adata[k][0][0][l]!=undef&&adata[k][1][0][l]!=undef)
						Bdata[k][0][0][l]=2*(adata[k][1][0][l]-adata[k][0][0][l])/dy/(float)sin(ydef[1]/2)-Bdata[k][1][0][l];
					
					else Bdata[k][0][0][l]=undef;
				}
			}
		}
		
		return B;
	}
	
	public Variable cCPrime(Variable gazm){
		t=gazm.getTCount();	z=gazm.getZCount();	y=gazm.getYCount();	x=gazm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		C=new Variable("Cc",gazm);	float undef=gazm.getUndef();
		C.setCommentAndUnit("coefficient C of elliptical equation");
		
		float[][]    buffer=new float[z-1][y];
		float[][][][] Cdata = C.getData();
		float[][][][] gdata =gazm.getData();
		
		if(gazm.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=1;k<z;k++)
				for(int j=1;j<y;j++){
					if(gdata[l][k-1][j][0]!=undef&&gdata[l][k][j][0]!=undef)
						buffer[k-1][j]=(gdata[l][k-1][j][0]+gdata[l][k][j][0])*(gdata[l][k-1][j][0]+gdata[l][k][j][0])/4;
						
					else buffer[k-1][j]=undef;
				}
				
				for(int k=1;k<z;k++)
				for(int j=1;j<y-1;j++){
					if(buffer[k-1][j]!=undef&&buffer[k-1][j+1]!=undef&&buffer[k-1][j-1]!=undef)
						Cdata[l][k][j][0]=(buffer[k-1][j+1]-buffer[k-1][j-1])/
							(2*dy*REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
					
					else Cdata[l][k][j][0]=undef;
				}
				
				for(int k=1;k<z;k++){
					if(buffer[k-1][y-1]!=undef&&buffer[k-1][y-2]!=undef)
						Cdata[l][k][y-1][0]=(buffer[k-1][y-1]-buffer[k-1][y-2])/
						(dy*REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else Cdata[l][k][y-1][0]=undef;
					
					Cdata[l][k][0][0]=undef;
				}
				
				for(int j=0;j<y;j++) Cdata[l][0][j][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=1;k<z;k++)
					for(int j=1;j<y;j++){
						if(gdata[k-1][j][0][l]!=undef&&gdata[k][j][0][l]!=undef)
							buffer[k-1][j]=(gdata[k-1][j][0][l]+gdata[k][j][0][l])*(gdata[k-1][j][0][l]+gdata[k][j][0][l])/4;
							
						else buffer[k-1][j]=undef;
					}
					
					for(int k=1;k<z;k++)
					for(int j=1;j<y-1;j++){
						if(buffer[k-1][j]!=undef&&buffer[k-1][j+1]!=undef&&buffer[k-1][j-1]!=undef)
							Cdata[k][j][0][l]=(buffer[k-1][j+1]-buffer[k-1][j-1])/
								(2*dy*REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
						
						else Cdata[k][j][0][l]=undef;
					}
					
					for(int k=1;k<z;k++){
						if(buffer[k-1][y-1]!=undef&&buffer[k-1][y-2]!=undef)
							Cdata[k][y-1][0][l]=(buffer[k-1][y-1]-buffer[k-1][y-2])/
							(dy*REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
						
						else Cdata[k][y-1][0][l]=undef;
						
						Cdata[k][0][0][l]=undef;
					}
					
					for(int j=0;j<y;j++) Cdata[0][j][0][l]=undef;
			}
		}
		
		return C;
	}
	
	public Variable cA(Variable sigmam){
		t=sigmam.getTCount();	z=sigmam.getZCount();	y=sigmam.getYCount();	x=sigmam.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable A=new Variable("Aa",sigmam);	float undef=sigmam.getUndef();
		A.setCommentAndUnit("coefficient A of elliptical equation");
		
		float[][][][] Adata=     A.getData();
		float[][][][] Sdata=sigmam.getData();
		
		if(sigmam.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y;j++)
				if(Sdata[l][k][j][0]!=undef) Adata[l][k][j][0]=Sdata[l][k][j][0]/bsin[j];
				else Adata[l][k][j][0]=undef;
				
				Adata[l][k][0][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y;j++)
				if(Sdata[k][j][0][l]!=undef)
				Adata[k][j][0][l]=Sdata[k][j][0][l]/bsin[j];
				else Adata[k][j][0][l]=undef;
				
				Adata[k][0][0][l]=undef;
			}
		}
		
		return A;
	}
	
	public Variable cB(Variable am){
		if(B==null) cBPrime(am);
		
		return B.copy();
	}
	
	public Variable cC(Variable gazm){
		t=gazm.getTCount();	z=gazm.getZCount();	y=gazm.getYCount();	x=gazm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable C=new Variable("Cc",gazm);	float undef=gazm.getUndef();
		C.setCommentAndUnit("coefficient C of elliptical equation");
		
		float[][]    buffer=new float[z][y];
		float[][][][] Cdata = C.getData();
		float[][][][] gdata =gazm.getData();
		
		if(gazm.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++)
				for(int j=1;j<y;j++){
					if(gdata[l][k][j][0]!=undef)
						buffer[k][j]=gdata[l][k][j][0]*gdata[l][k][j][0];
						
					else buffer[k][j]=undef;
				}
				
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++){
					if(buffer[k][j+1]!=undef&&buffer[k][j-1]!=undef)
						Cdata[l][k][j][0]=(buffer[k][j+1]-buffer[k][j-1])/
							(2*dy*REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
					
					else Cdata[l][k][j][0]=undef;
				}
				
				for(int k=0;k<z;k++){
					if(buffer[k][y-1]!=undef&&buffer[k][y-2]!=undef)
						Cdata[l][k][y-1][0]=(buffer[k][y-1]-buffer[k][y-2])/
						(dy*REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else Cdata[l][k][y-1][0]=undef;
					
					Cdata[l][k][0][0]=undef;
				}
				
				for(int j=0;j<y;j++) Cdata[l][0][j][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++)
					for(int j=1;j<y;j++){
						if(gdata[k][j][0][l]!=undef)
							buffer[k][j]=gdata[k][j][0][l]*gdata[k][j][0][l];
							
						else buffer[k][j]=undef;
					}
					
					for(int k=0;k<z;k++)
					for(int j=1;j<y-1;j++){
						if(buffer[k][j+1]!=undef&&buffer[k][j-1]!=undef)
							Cdata[k][j][0][l]=(buffer[k][j+1]-buffer[k][j-1])/
								(2*dy*REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]*bsin[j]);
						
						else Cdata[k][j][0][l]=undef;
					}
					
					for(int k=0;k<z;k++){
						if(buffer[k][y-1]!=undef&&buffer[k][y-2]!=undef)
							Cdata[k][y-1][0][l]=(buffer[k][y-1]-buffer[k][y-2])/
							(dy*REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]*bsin[y-1]);
						
						else Cdata[k][y-1][0][l]=undef;
						
						Cdata[k][0][0][l]=undef;
					}
					
					for(int j=0;j<y;j++) Cdata[0][j][0][l]=undef;
			}
		}
		
		return C;
	}
	
	public Variable cBSinB(Variable am){
		t=am.getTCount();	z=am.getZCount();	y=am.getYCount();	x=am.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable B=new Variable("BsinB",am);	float undef=am.getUndef();
		B.setCommentAndUnit("coefficient B multiplied by sinB");
		
		float[][][][] Bdata= B.getData();
		float[][][][] adata=am.getData();
		
		if(am.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++){
					for(int j=1;j<y-1;j++)
					if(adata[l][k][j+1][0]!=undef&&adata[l][k][j-1][0]!=undef)
						Bdata[l][k][j][0]=(adata[l][k][j+1][0]-adata[l][k][j-1][0])/(2*dy);
					
					else Bdata[l][k][j][0]=undef;
					
					if(adata[l][k][y-1][0]!=undef&&adata[l][k][y-2][0]!=undef)
						Bdata[l][k][y-1][0]=2*(adata[l][k][y-1][0]-adata[l][k][y-2][0])/dy-Bdata[l][k][y-2][0];
					
					else Bdata[l][k][y-1][0]=undef;
					
					if(adata[l][k][0][0]!=undef&&adata[l][k][1][0]!=undef)
						Bdata[l][k][0][0]=2*(adata[l][k][1][0]-adata[l][k][0][0])/dy-Bdata[l][k][1][0];
					
					else Bdata[l][k][0][0]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++){
					for(int j=1;j<y-1;j++)
					if(adata[k][j+1][0][l]!=undef&&adata[k][j-1][0][l]!=undef)
						Bdata[k][j][0][l]=(adata[k][j+1][0][l]-adata[k][j-1][0][l])/(2*dy);
					
					else Bdata[k][j][0][l]=undef;
					
					if(adata[k][y-1][0][l]!=undef&&adata[k][y-2][0][l]!=undef)
						Bdata[k][y-1][0][l]=2*(adata[k][y-1][0][l]-adata[k][y-2][0][l])/dy-Bdata[k][y-2][0][l];
					
					else Bdata[k][y-1][0][l]=undef;
					
					if(adata[k][0][0][l]!=undef&&adata[k][1][0][l]!=undef)
						Bdata[k][0][0][l]=2*(adata[k][1][0][l]-adata[k][0][0][l])/dy-Bdata[k][1][0][l];
					
					else Bdata[k][0][0][l]=undef;
				}
			}
		}
		
		return B;
	}
	
	public Variable cCSinB(Variable gazm){
		t=gazm.getTCount();	z=gazm.getZCount();	y=gazm.getYCount();	x=gazm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		
		Variable C=new Variable("CsinB",gazm);	float undef=gazm.getUndef();
		C.setCommentAndUnit("coefficient C multiplied by sinB");
		
		float[][]    buffer=new float[z][y];
		float[][][][] Cdata = C.getData();
		float[][][][] gdata =gazm.getData();
		
		if(gazm.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++)
				for(int j=1;j<y;j++){
					if(gdata[l][k][j][0]!=undef)
						buffer[k][j]=gdata[l][k][j][0]*gdata[l][k][j][0];
						
					else buffer[k][j]=undef;
				}
				
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++){
					if(buffer[k][j+1]!=undef&&buffer[k][j-1]!=undef)
						Cdata[l][k][j][0]=(buffer[k][j+1]-buffer[k][j-1])/
							(2*dy*REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]);
					
					else Cdata[l][k][j][0]=undef;
				}
				
				for(int k=0;k<z;k++){
					if(buffer[k][y-1]!=undef&&buffer[k][y-2]!=undef)
						Cdata[l][k][y-1][0]=(buffer[k][y-1]-buffer[k][y-2])/
						(dy*REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else Cdata[l][k][y-1][0]=undef;
					
					Cdata[l][k][0][0]=undef;
				}
				
				for(int j=0;j<y;j++) Cdata[l][0][j][0]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++)
				for(int j=1;j<y;j++){
					if(gdata[k][j][0][l]!=undef)
						buffer[k][j]=gdata[k][j][0][l]*gdata[k][j][0][l];
						
					else buffer[k-1][j]=undef;
				}
				
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++){
					if(buffer[k][j+1]!=undef&&buffer[k][j-1]!=undef)
						Cdata[k][j][0][l]=(buffer[k][j+1]-buffer[k][j-1])/
							(2*dy*REarth*REarth*REarth*bsin[j]*bsin[j]*bsin[j]);
					
					else Cdata[k][j][0][l]=undef;
				}
				
				for(int k=0;k<z;k++){
					if(buffer[k][y-1]!=undef&&buffer[k][y-2]!=undef)
						Cdata[k][y-1][0][l]=(buffer[k][y-1]-buffer[k][y-2])/
						(dy*REarth*REarth*REarth*bsin[y-1]*bsin[y-1]*bsin[y-1]);
					
					else Cdata[k][y-1][0][l]=undef;
					
					Cdata[k][0][0][l]=undef;
				}
				
				for(int j=0;j<y;j++) Cdata[0][j][0][l]=undef;
			}
		}
		
		return C;
	}
	
	
	/**
     * initialize the boundaries of stream function using vm and wm
     *
     * @param	sf	stream function
     * @param	vm	tangential averaged radial wind
     * @param	wm	tangential averaged vertical velocity
     */
	public void initialSFBoundary(Variable sf,Variable vm,Variable wm){
		t=sf.getTCount();	z=sf.getZCount();	y=sf.getYCount();	x=sf.getXCount();
		
		if(!sf.isLike(vm)||!sf.isLike(wm)) throw new IllegalArgumentException("dimensions not same");
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
	
	
	/**
     * compute the left hand side of the equation
     *
     * @param	gm	tangentially averaged absolute angular momentum
     * @param	vm	tangentially averaged radial wind
     * @param	wm	tangentially averaged vertical velocity
     * @param	sm	tangentially averaged static stability
     * @param	am	tangentially averaged specific volume
     * 
     * @return	re	[0] is term 1 and [1] is term 2 before radial partial difference
     * 				[2] is term 3 and [3] is term 4 before vertical partial difference
     */
	public Variable[] cLeft(Variable gm,Variable vm,Variable wm,Variable sm,Variable am){
		t=gm.getTCount();	z=gm.getZCount();	y=gm.getYCount();	x=gm.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		if(!gm.isLike(vm)||!gm.isLike(wm)||!gm.isLike(sm)||!gm.isLike(am))
			throw new IllegalArgumentException("dimensions not same");
		
		float undef=gm.getUndef();
		
		Range range=new Range(t,z,y,1);
		
		Variable[] re=new Variable[4];
		re[0]=new Variable("lp1",gm.isTFirst(),range);
		re[1]=new Variable("lp2",gm.isTFirst(),range);
		re[2]=new Variable("lp3",gm.isTFirst(),range);
		re[3]=new Variable("lp4",gm.isTFirst(),range);
		
		re[0].setUndef(undef); re[0].setCommentAndUnit("term 1"); re[0].setValue(undef);
		re[1].setUndef(undef); re[1].setCommentAndUnit("term 2"); re[1].setValue(undef);
		re[2].setUndef(undef); re[2].setCommentAndUnit("term 3"); re[2].setValue(undef);
		re[3].setUndef(undef); re[3].setCommentAndUnit("term 4"); re[3].setValue(undef);
		
		float[][][][] lp1data=re[0].getData();
		float[][][][] lp2data=re[1].getData();
		float[][][][] lp3data=re[2].getData();
		float[][][][] lp4data=re[3].getData();
		float[][][][]  gmdata=gm.getData();
		float[][][][]  vmdata=vm.getData();
		float[][][][]  wmdata=wm.getData();
		float[][][][]  smdata=sm.getData();
		float[][][][]  amdata=am.getData();
		
		if(gm.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++){
				float gm2=gmdata[l][k][j+1][0]*gmdata[l][k][j+1][0];
				float gm1=gmdata[l][k][j-1][0]*gmdata[l][k][j-1][0];
				float fct=(float)Math.pow(REarth*bsin[j],3);
				
				lp1data[l][k][j][0]=-smdata[l][k][j][0]* wmdata[l][k][j][0];
				lp2data[l][k][j][0]= vmdata[l][k][j][0]*(amdata[l][k][j+1][0]-amdata[l][k][j-1][0])/(2*dy);
				lp3data[l][k][j][0]=-wmdata[l][k][j][0]*(amdata[l][k][j+1][0]-amdata[l][k][j-1][0])/(2*dy);
				lp4data[l][k][j][0]= vmdata[l][k][j][0]*(gm2-gm1)/(2*dy)/fct;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++){
				float gm2=gmdata[k][j+1][0][l]*gmdata[k][j+1][0][l];
				float gm1=gmdata[k][j-1][0][l]*gmdata[k][j-1][0][l];
				float fct=(float)Math.pow(REarth*bsin[j],3);
				
				lp1data[k][j][0][l]=-smdata[k][j][0][l]* wmdata[k][j][0][l];
				lp2data[k][j][0][l]= vmdata[k][j][0][l]*(amdata[k][j+1][0][l]-amdata[k][j-1][0][l])/(2*dy);
				lp3data[k][j][0][l]=-wmdata[k][j][0][l]*(amdata[k][j+1][0][l]-amdata[k][j-1][0][l])/(2*dy);
				lp4data[k][j][0][l]= vmdata[k][j][0][l]*(gm2-gm1)/(2*dy)/fct;
			}
		}
		
		return re;
	}
	
	
	/**
     * implement the methods in SORInterface
     */
	public void sor(float tol,int loop,Variable sf){
		if(print) System.out.println("\nStart SORing...");
		t=sf.getTCount();	z=sf.getZCount();	y=sf.getYCount();	x=sf.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be averaged");
		if(tol<0) throw new IllegalArgumentException();
		if(!A.isLike(B)||!A.isLike(C)) throw new IllegalArgumentException("Dimension not the same");
		
		tstart=sf.getRange().getTRange()[0];
		
		float[][][][]  Adata= A.getData();
		float[][][][]  Bdata= B.getData();
		float[][][][]  Cdata= C.getData();
		float[][][][] sfdata=sf.getData();
		
		MDate[] times=sm.getTDef().getSamples();
		
		int xx=0;	// loop record
		float max_error,relative_e,b,tmp_max=Float.MIN_VALUE;
		float factor2=dz/dy,factor1=factor2*factor2;	factor2/=4;
		
		/*** calculate optimizing argument ***/
		float epsilon=(float)(pow(sin(PI/(2*y+2)),2)+pow(sin(PI/(2*z+2)),2));
		float opt_arg=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		
		float[][] R=new float[z-2][y-2];
		
		final int mcount=loop>0?loop:5999;	// max loop count
		
		if(sf.isTFirst()){
			for(int l=0;l<t;l++){
				do{
					for(int k=1;k<z-1;k++)
					for(int j=1;j<y-1;j++){
						b=(sfdata[l][k][j][0]<0)?(-sfdata[l][k][j][0]):sfdata[l][k][j][0];
						if(b>tmp_max) tmp_max=b;
						
						R[k-1][j-1]=(
								Adata[l][k][j+1][0]*(sfdata[l][k  ][j+1][0]-sfdata[l][k  ][j  ][0])-
								Adata[l][k][j  ][0]*(sfdata[l][k  ][j  ][0]-sfdata[l][k  ][j-1][0])
							)*factor1+(
								Bdata[l][k][j+1][0]*(sfdata[l][k+1][j+1][0]-sfdata[l][k-1][j+1][0])-
								Bdata[l][k][j-1][0]*(sfdata[l][k+1][j-1][0]-sfdata[l][k-1][j-1][0])
							)*factor2+(
								Bdata[l][k+1][j][0]*(sfdata[l][k+1][j+1][0]-sfdata[l][k+1][j-1][0])-
								Bdata[l][k-1][j][0]*(sfdata[l][k-1][j+1][0]-sfdata[l][k-1][j-1][0])
							)*factor2+(
								Cdata[l][k+1][j][0]*(sfdata[l][k+1][j  ][0]-sfdata[l][k  ][j  ][0])-
								Cdata[l][k  ][j][0]*(sfdata[l][k  ][j  ][0]-sfdata[l][k-1][j  ][0])
							);
						
						R[k-1][j-1]*=opt_arg/((Adata[l][k][j+1][0]+Adata[l][k][j][0])*factor1+(Cdata[l][k+1][j][0]+Cdata[l][k][j][0]));
						
						sfdata[l][k][j][0]+=R[k-1][j-1];
					}
					
					max_error =getAbsMax(R);
					relative_e=max_error/tmp_max;	xx++;
					
				}while(xx<=mcount&&relative_e>=tol);
				
				if(print) System.out.println(times[l+tstart-1]+"\tloops "+xx+"  RE is "+relative_e);
				xx=0;	tmp_max=Float.MIN_VALUE;
			}
			
		}else{
			for(int l=0;l<t;l++){
				do{
					for(int k=1;k<z-1;k++)
					for(int j=1;j<y-1;j++){
						b=(sfdata[k][j][0][l]<0)?(-sfdata[k][j][0][l]):sfdata[k][j][0][l];
						if(b>tmp_max) tmp_max=b;
						
						R[k-1][j-1]=(
								Adata[k][j+1][0][l]*(sfdata[k  ][j+1][0][l]-sfdata[k  ][j  ][0][l])-
								Adata[k][j  ][0][l]*(sfdata[k  ][j  ][0][l]-sfdata[k  ][j-1][0][l])
							)*factor1+(
								Bdata[k][j+1][0][l]*(sfdata[k+1][j+1][0][l]-sfdata[k-1][j+1][0][l])-
								Bdata[k][j-1][0][l]*(sfdata[k+1][j-1][0][l]-sfdata[k-1][j-1][0][l])
							)*factor2+(
								Bdata[k+1][j][0][l]*(sfdata[k+1][j+1][0][l]-sfdata[k+1][j-1][0][l])-
								Bdata[k-1][j][0][l]*(sfdata[k-1][j+1][0][l]-sfdata[k-1][j-1][0][l])
							)*factor2+(
								Cdata[k+1][j][0][l]*(sfdata[k+1][j  ][0][l]-sfdata[k  ][j  ][0][l])-
								Cdata[k  ][j][0][l]*(sfdata[k  ][j  ][0][l]-sfdata[k-1][j  ][0][l])
							);
						
						R[k-1][j-1]*=opt_arg/((Adata[k][j+1][0][l]+Adata[k][j][0][l])*factor1+(Cdata[k+1][j][0][l]+Cdata[k][j][0][l]));
						
						sfdata[k][j][0][l]+=R[k-1][j-1];
					}
					
					max_error =getAbsMax(R);
					relative_e=max_error/tmp_max;	xx++;
					
				}while(xx<=mcount&&relative_e>=tol);
				
				if(print) System.out.println(times[l+tstart-1]+"\tloops "+xx+"  RE is "+relative_e);
				xx=0;	tmp_max=Float.MIN_VALUE;
			}
		}
		
		if(print) System.out.println("Finish SORing.");
	}
	
	public void sor(float tol,int loop,Variable sf,Variable F){
		if(print) System.out.println("\nStart SORing...");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		if(!A.isLike(B)||!A.isLike(C)||!A.isLike(F)) throw new IllegalArgumentException("dimensions not same");
		if(tol<0||x!=1) throw new IllegalArgumentException();
		
		tstart=sf.getRange().getTRange()[0];
		
		float undef=F.getUndef();
		float[][][][]  Fdata= F.getData();
		float[][][][]  Adata= A.getData();
		float[][][][]  Bdata= B.getData();
		float[][][][]  Cdata= C.getData();
		float[][][][] sfdata=sf.getData();
		
		MDate[] times=sm.getTDef().getSamples();
		
		int xx=0;	// loop record
		float max_error,relative_e=2,b,tmp_max=Float.MIN_VALUE,dp2=dz*dz;
		float factor2=dz/dy,factor1=factor2*factor2;	factor2/=4;
		
		/*** calculate optimizing argument ***/
		float epsilon=(float)(pow(sin(PI/(2*y+2)),2)+pow(sin(PI/(2*z+2)),2));
		float opt_arg=(float)(2/(1+sqrt((2-epsilon)*epsilon)));
		
		float[][] R=new float[z-2][y-2];
		
		final int mcount=loop>0?loop:5999;	// max loop count
		
		if(sf.isTFirst()){
			for(int l=0;l<t;l++){
				do{
					for(int k=1;k<z-1;k++)
					for(int j=1;j<y-1;j++){
						if(Fdata[l][k][j][0]==undef||Adata[l][k][j+1][0]==undef||Adata[l][k][j][0]==undef||
						Bdata[l][k][j+1][0]==undef||Bdata[l][k][j-1][0]==undef||Bdata[l][k+1][j][0]==undef||
						Bdata[l][k-1][j][0]==undef||Cdata[l][k+1][j][0]==undef||Cdata[l][k  ][j][0]==undef)
						continue;
						
						b=(sfdata[l][k][j][0]<0)?(-sfdata[l][k][j][0]):sfdata[l][k][j][0];
						if(b>tmp_max) tmp_max=b;
						
						R[k-1][j-1]=(
							(
								Adata[l][k][j+1][0]*(sfdata[l][k  ][j+1][0]-sfdata[l][k  ][j  ][0])-
								Adata[l][k][j  ][0]*(sfdata[l][k  ][j  ][0]-sfdata[l][k  ][j-1][0])
							)*factor1+(
								Bdata[l][k][j+1][0]*(sfdata[l][k+1][j+1][0]-sfdata[l][k-1][j+1][0])-
								Bdata[l][k][j-1][0]*(sfdata[l][k+1][j-1][0]-sfdata[l][k-1][j-1][0])
							)*factor2+(
								Bdata[l][k+1][j][0]*(sfdata[l][k+1][j+1][0]-sfdata[l][k+1][j-1][0])-
								Bdata[l][k-1][j][0]*(sfdata[l][k-1][j+1][0]-sfdata[l][k-1][j-1][0])
							)*factor2+(
								Cdata[l][k+1][j][0]*(sfdata[l][k+1][j  ][0]-sfdata[l][k  ][j  ][0])-
								Cdata[l][k  ][j][0]*(sfdata[l][k  ][j  ][0]-sfdata[l][k-1][j  ][0])
							)
						)-Fdata[l][k][j][0]*dp2;
						
						R[k-1][j-1]*=opt_arg/((Adata[l][k][j+1][0]+Adata[l][k][j][0])*factor1+(Cdata[l][k+1][j][0]+Cdata[l][k][j][0]));
						
						sfdata[l][k][j][0]+=R[k-1][j-1];
					}
					
					max_error=getAbsMax(R);
					if(tmp_max!=0) relative_e=max_error/tmp_max;	xx++;
					
				}while(xx<=mcount&&relative_e>=tol);
				
				if(print) System.out.println(times[l+tstart-1]+"\tloops "+xx+"  RE is "+relative_e);
				xx=0;	tmp_max=Float.MIN_VALUE;
			}
			
		}else{
			for(int l=0;l<t;l++){
				do{
					for(int k=1;k<z-1;k++)
					for(int j=1;j<y-1;j++){
						if(Fdata[k][j][0][l]==undef||Adata[k][j+1][0][l]==undef||Adata[k][j][0][l]==undef||
						Bdata[k][j+1][0][l]==undef||Bdata[k][j-1][0][l]==undef||Bdata[k+1][j][0][l]==undef||
						Bdata[k-1][j][0][l]==undef||Cdata[k+1][j][0][l]==undef||Cdata[k  ][j][0][l]==undef)
						continue;
						
						b=(sfdata[k][j][0][l]<0)?(-sfdata[k][j][0][l]):sfdata[k][j][0][l];
						if(b>tmp_max) tmp_max=b;
						
						R[k-1][j-1]=(
							(
								Adata[k][j+1][0][l]*(sfdata[k  ][j+1][0][l]-sfdata[k  ][j  ][0][l])-
								Adata[k][j  ][0][l]*(sfdata[k  ][j  ][0][l]-sfdata[k  ][j-1][0][l])
							)*factor1+(
								Bdata[k][j+1][0][l]*(sfdata[k+1][j+1][0][l]-sfdata[k-1][j+1][0][l])-
								Bdata[k][j-1][0][l]*(sfdata[k+1][j-1][0][l]-sfdata[k-1][j-1][0][l])
							)*factor2+(
								Bdata[k+1][j][0][l]*(sfdata[k+1][j+1][0][l]-sfdata[k+1][j-1][0][l])-
								Bdata[k-1][j][0][l]*(sfdata[k-1][j+1][0][l]-sfdata[k-1][j-1][0][l])
							)*factor2+(
								Cdata[k+1][j][0][l]*(sfdata[k+1][j  ][0][l]-sfdata[k  ][j  ][0][l])-
								Cdata[k  ][j][0][l]*(sfdata[k  ][j  ][0][l]-sfdata[k-1][j  ][0][l])
							)
						)-Fdata[k][j][0][l]*dp2;
						
						R[k-1][j-1]*=opt_arg/((Adata[k][j+1][0][l]+Adata[k][j][0][l])*factor1+(Cdata[k+1][j][0][l]+Cdata[k][j][0][l]));
						
						sfdata[k][j][0][l]+=R[k-1][j-1];
					}
					
					max_error =getAbsMax(R);
					if(tmp_max!=0) relative_e=max_error/tmp_max;	xx++;
					
				}while(xx<=mcount);
				
				if(print) System.out.println(times[l+tstart-1]+"\tloops "+xx+"  RE is "+relative_e);
				xx=0;	tmp_max=Float.MIN_VALUE;
			}
		}
		
		if(print) System.out.println("Finish SORing.");
	}
	
	
	/**
     * calculate relative velocity according to the moving speed of the coordinate
     *
     * @param		ut	tangential wind speed
     * @param		vr	radial wind speed
     */
    public void cRelativeVelocity(Variable ut,Variable vr){
    	if(!ut.isLike(vr)) throw new IllegalArgumentException("dimensions not same");
    	
    	CylindricalSpatialModel csm=(CylindricalSpatialModel)sm;
    	
		t=ut.getTCount();	z=ut.getZCount();	y=ut.getYCount();	x=ut.getXCount();
		tstart=ut.getRange().getTRange()[0];
    	
    	float[] olon=csm.getOLon();	float[] cu=csm.getUWhole();
    	float[] olat=csm.getOLat();	float[] cv=csm.getVWhole();
    	
    	float[][][][] utdata=ut.getData();
    	float[][][][] vrdata=vr.getData();
    	
    	cs=new float[t];	// central speed of the model
    	cd=new float[t];	// central direction of the speed
    	
    	wtan=new float[t][x];
    	wnor=new float[t][x];
    	
    	if(t==1) cs[0]=0;
    	else{
    		for(int l=1;l<t-1;l++) cs[l]=cSphericalDistanceByRadian(
    			olon[tstart-2+l],olat[tstart-2+l],olon[tstart+l],olat[tstart+l]
    		)/(csm.getDT()*2);
    		
    		cs[0  ]=cSphericalDistanceByRadian(
    			olon[tstart-1  ],olat[tstart-1  ],olon[tstart    ],olat[tstart    ]
    		)/csm.getDT();
    		
    		cs[t-1]=cSphericalDistanceByRadian(
    			olon[tstart-3+t],olat[tstart-3+t],olon[tstart-2+t],olat[tstart-2+t]
    		)/csm.getDT();
    	}
    	
    	if(ut.isTFirst()){
	    	for(int l=0;l<t;l++){
	    		cd[l]=(float)(Math.atan2(cv[tstart-1+l],cu[tstart-1+l])-PI/2); // angle start from north, counter-clockwise
	    		
	    		for(int i=0;i<x;i++){
	    			float angle=cd[l]-xdef[i];
	    			wtan[l][i]=cs[l]*(float)sin(angle);
	    			wnor[l][i]=cs[l]*(float)cos(angle);
	    			
		    		for(int k=0;k<z;k++)
		    		for(int j=0;j<y;j++){
		    			utdata[l][k][j][i]-=wtan[l][i];
		    			vrdata[l][k][j][i]-=wnor[l][i];
		    		}
	    		}
	    	}
    		
    	}else{
	    	for(int l=0;l<t;l++){
	    		cd[l]=(float)(Math.atan2(cv[tstart-1+l],cu[tstart-1+l])-PI/2); // angle start from north, counter-clockwise
	    		
	    		for(int i=0;i<x;i++){
	    			float angle=cd[l]-xdef[i];
	    			wtan[l][i]=cs[l]*(float)sin(angle);
	    			wnor[l][i]=cs[l]*(float)cos(angle);
	    			
		    		for(int k=0;k<z;k++)
		    		for(int j=0;j<y;j++){
		    			utdata[k][j][i][l]-=wtan[l][i];
		    			vrdata[k][j][i][l]-=wnor[l][i];
		    		}
	    		}
	    	}
    	}
    }
	
    
    /**
     * calculate v and w, a vector, [0] is in y direction while [1] is in z direction
     *
     * @param	sf	stream function
     *
     * @return	vw (vector)
     */
	public Variable[] cVW(Variable sf){
		t=sf.getTCount();	z=sf.getZCount();	y=sf.getYCount();	x=sf.getXCount();
		
		if(x!=1) throw new IllegalArgumentException("x-direction should be only one point");
		
		float undef=sf.getUndef();
		Variable[] vw=new Variable[2];
		
		vw[0]=new Variable("vs",sf); vw[1]=new Variable("ws",sf);
		
		vw[0].setValue(undef);	vw[0].setCommentAndUnit("simulated radial velocity (m s^-1)");
		vw[1].setValue(undef);	vw[1].setCommentAndUnit("simulated vertical velocity (Pa s^-1)");
		
		float[][][][] sfdata=   sf.getData();
		float[][][][] vsdata=vw[0].getData();
		float[][][][] wsdata=vw[1].getData();
		
		if(sf.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=1;k<z-1;k++)
				for(int j=1;j<y;j++){
					if(sfdata[l][k+1][j][0]!=undef&&sfdata[l][k-1][j][0]!=undef)
						vsdata[l][k][j][0]=(sfdata[l][k+1][j][0]-sfdata[l][k-1][j][0])/bsin[j]/(dz*2);
						
					else vsdata[l][k][j][0]=undef;
				}
				
				/*** top and bottom boundary **
				for(int j=1;j<y;j++){
					if(sfdata[l][z-1][j][0]!=undef&&sfdata[l][z-2][j][0]!=undef)
						vsdata[l][z-1][j][0]=2*(sfdata[l][z-1][j][0]-sfdata[l][z-2][j][0])/bsin[j]/dz[0]-vsdata[l][z-2][j][0];
						
					else vsdata[l][z-1][j][0]=undef;
					
					if(sfdata[l][0][j][0]!=undef&&sfdata[l][1][j][0]!=undef)
						vsdata[l][0][j][0]=2*(sfdata[l][1][j][0]-sfdata[l][0][j][0])/bsin[j]/dz[0]-vsdata[l][1][j][0];
						
					else vsdata[l][0][j][0]=undef;
				}*/
				
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++){
					if(sfdata[l][k][j+1][0]!=undef&&sfdata[l][k][j-1][0]!=undef)
						wsdata[l][k][j][0]=-(sfdata[l][k][j+1][0]-sfdata[l][k][j-1][0])/bsin[j]/(dy*2);
						
					else wsdata[l][k][j][0]=undef;
				}
				
				/*** left and right boundary **
				for(int k=0;k<z;k++){
					if(sfdata[l][k][y-1][0]!=undef&&sfdata[l][k][y-2][0]!=undef)
						wsdata[l][k][y-1][0]=-2*(sfdata[l][k][y-1][0]-sfdata[l][k][y-2][0])/
							(float)sin((ydef[y-1]+ydef[y-2])/2)/dy-wsdata[l][k][y-2][0];
						
					else wsdata[l][k][y-1][0]=undef;
					
					if(sfdata[l][k][1][0]!=undef&&sfdata[l][k][0][0]!=undef)
						wsdata[l][k][0][0]=-2*(sfdata[l][k][1][0]-sfdata[l][k][0][0])/
							(float)sin((ydef[0]+ydef[1])/2)/dy-wsdata[l][k][1][0];
						
					else wsdata[l][k][0][0]=undef;
				}*/
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=1;k<z-1;k++)
				for(int j=1;j<y;j++){
					if(sfdata[k+1][j][0][l]!=undef&&sfdata[k-1][j][0][l]!=undef)
						vsdata[k][j][0][l]=(sfdata[k+1][j][0][l]-sfdata[k-1][j][0][l])/bsin[j]/(dz*2);
						
					else vsdata[k][j][0][l]=undef;
				}
				
				/*** top and bottom boundary **
				for(int j=1;j<y;j++){
					if(sfdata[z-1][j][0][l]!=undef&&sfdata[z-2][j][0][l]!=undef)
						vsdata[z-1][j][0][l]=2*(sfdata[z-1][j][0][l]-sfdata[z-2][j][0][l])/bsin[j]/dz[0]-vsdata[z-2][j][0][l];
						
					else vsdata[z-1][j][0][l]=undef;
					
					if(sfdata[0][j][0][l]!=undef&&sfdata[1][j][0][l]!=undef)
						vsdata[0][j][0][l]=2*(sfdata[1][j][0][l]-sfdata[0][j][0][l])/bsin[j]/dz[0]-vsdata[1][j][0][l];
						
					else vsdata[0][j][0][l]=undef;
				}*/
				
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++){
					if(sfdata[k][j+1][0][l]!=undef&&sfdata[k][j-1][0][l]!=undef)
						wsdata[k][j][0][l]=-(sfdata[k][j+1][0][l]-sfdata[k][j-1][0][l])/bsin[j]/(dy*2);
						
					else wsdata[k][j][0][l]=undef;
				}
				
				/*** left and right boundary **
				for(int k=0;k<z;k++){
					if(sfdata[k][y-1][0][l]!=undef&&sfdata[k][y-2][0][l]!=undef)
						wsdata[k][y-1][0][l]=-2*(sfdata[k][y-1][0][l]-sfdata[k][y-2][0][l])/
							(float)sin((ydef[y-1]+ydef[y-2])/2)/dy-wsdata[k][y-2][0][l];
						
					else wsdata[k][y-1][0][l]=undef;
					
					if(sfdata[k][1][0][l]!=undef&&sfdata[k][0][0][l]!=undef)
						wsdata[k][0][0][l]=-2*(sfdata[k][1][0][l]-sfdata[k][0][0][l])/
							(float)sin((ydef[0]+ydef[1])/2)/dy-wsdata[k][1][0][l];
						
					else wsdata[k][0][0][l]=undef;
				}*/
			}
		}
		
		return vw;
	}
	
	
	/**
	 * whether to print out
	 *
     * @param	print	print or disable print
     */ 
	public void setPrinting(boolean print){ this.print=print;}
	
	
	/**
     * mark the level
     *
     * @param	sfp		surface pressure
     */
	private void cTopography(Variable sfp){
		t=sfp.getTCount();	z=sfp.getZCount();	y=sfp.getYCount();	x=sfp.getXCount();
		
		if(z!=1) throw new IllegalArgumentException("surface data only");
		
		tags=new short[t][y][x];
		
		float[][][][] sfpdata=sfp.getData();
		
		float tmp=zdef[zdef.length-1]-zdef[0];
		
		if(tmp==0) throw new IllegalArgumentException("zdef is only one level or not monototic");
		
		boolean incre=tmp>0?true:false;
		
		if(sfp.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				tags[l][j][i]=incre?
				(short)ArrayUtil.getIdxIncre(zdef,sfpdata[l][0][j][i]):
				(short)ArrayUtil.getIdxDecre(zdef,sfpdata[l][0][j][i]);
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				tags[l][j][i]=incre?
				(short)ArrayUtil.getIdxIncre(zdef,sfpdata[0][j][i][l]):
				(short)ArrayUtil.getIdxDecre(zdef,sfpdata[0][j][i][l]);
		}
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
