/**
 * @(#)PotentialVorticityInSC.java	1.0 2014.07.22
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.application.EquationInSphericalCoordinate;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;
import static miniufo.geophysics.atmos.ThermoDynamics.Cp;
import static miniufo.diagnosis.SpatialModel.GRAVITY_ACCERLERATION;


/**
 * potential vorticity in spheral coordinate
 *
 * @version 1.0, 2014.07.22
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class PotentialVorticityInSC extends EquationInSphericalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public PotentialVorticityInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/**
     * calculate quasi-geostrophic potential vorticity terms
     *
     * @param	fai		geopotential (m^2/s^2)
     * @param	T		temperature(K)
     *
     * @return	QGPV	quasi-geostrophic potential vorticity (1 PVU = 10^-6 K m^2 Kg^-1 s^-1)
     */
	public Variable[] cQGPVTerms(Variable fai,Variable T){
		if(!fai.isLike(T)) throw new IllegalArgumentException("dimensions not same");
		
		ystart=fai.getRange().getYRange()[0];
		
		t=fai.getTCount();	z=fai.getZCount();	y=fai.getYCount();	x=fai.getXCount();
		
		Variable[] pv=new Variable[3];
		pv[0]=new Variable("fterm",fai);	pv[0].setCommentAndUnit("planetary vorticity term (10^-6 K m^2 Kg^-1 s^-1)");
		pv[1]=new Variable("hterm",fai);	pv[1].setCommentAndUnit("relative vorticity term (10^-6 K m^2 Kg^-1 s^-1)");
		pv[2]=new Variable("vterm",fai);	pv[2].setCommentAndUnit("vertical term (10^-6 K m^2 Kg^-1 s^-1)");

		float[][][][]  fdata=  fai.getData();
		float[][][][]  Tdata=    T.getData();
		float[][][][] q1data=pv[0].getData();
		float[][][][] q2data=pv[1].getData();
		float[][][][] q3data=pv[2].getData();
		
		float undef=fai.getUndef();
		
		if(fai.isTFirst()){
			for(int l=0;l<t;l++){
				// term 1
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) q1data[l][k][j][i]=f1[ystart-1+j];
				
				// term 2
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(fdata[l][k][j][i+1]!=undef&&fdata[l][k][j][i-1]!=undef
				&&fdata[l][k][j+1][i]!=undef&&fdata[l][k][j-1][i]!=undef)
					q2data[l][k][j][i]=(
						(fdata[l][k][j][i+1]+fdata[l][k][j][i-1]-2*fdata[l][k][j][i])/(dxs[ystart-1+j]*dxs[ystart-1+j])+
						(fdata[l][k][j+1][i]+fdata[l][k][j-1][i]-2*fdata[l][k][j][i])/(dy*dy)
					)/f1[ystart-1+j];
				else
					q2data[l][k][j][i]=undef;
				
				// term 3
				for(int k=1;k<z-1;k++){
					float tmp=Rd/zdef[zstart-1+k];
					
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++)
					if(fdata[l][k+1][j][i]!=undef&&fdata[l][k-1][j][i]!=undef&&
					Tdata[l][k][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef&&Tdata[l][k+1][j][i]!=undef){
						float thetaU=(float)(Tdata[l][k+1][j][i]*pow(zdef[zstart-1+k+1]/100000,Rd/Cp));
						float thetaM=(float)(Tdata[l][k  ][j][i]*pow(zdef[zstart-1+k  ]/100000,Rd/Cp));
						float thetaD=(float)(Tdata[l][k-1][j][i]*pow(zdef[zstart-1+k-1]/100000,Rd/Cp));
						
						q3data[l][k][j][i]=f1[ystart-1+j]*(
							(fdata[l][k+1][j][i]-fdata[l][k][j][i])/
							(tmp*(Tdata[l][k+1][j][i]+Tdata[l][k][j][i])/2*(float)log(thetaU/thetaM))-
							(fdata[l][k][j][i]-fdata[l][k-1][j][i])/
							(tmp*(Tdata[l][k][j][i]+Tdata[l][k-1][j][i])/2*(float)log(thetaM/thetaD))
						)/((zdef[zstart-1+k]+zdef[zstart-2+k])/2);
						
					}else q3data[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++){
				// term 1
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) q1data[k][j][i][l]=f1[j];
				
				// term 2
				for(int k=0;k<z;k++)
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(fdata[k][j][i+1][l]!=undef&&fdata[k][j][i-1][l]!=undef
				&&fdata[k][j+1][i][l]!=undef&&fdata[k][j-1][i][l]!=undef)
					q2data[k][j][i][l]=(
						(fdata[k][j][i+1][l]+fdata[k][j][i-1][l]-2*fdata[k][j][i][l])/(dxs[ystart-1+j]*dxs[ystart-1+j])+
						(fdata[k][j+1][i][l]+fdata[k][j-1][i][l]-2*fdata[k][j][i][l])/(dy*dy)
					)/f1[ystart-1+j];
				else
					q2data[k][j][i][l]=undef;
				
				// term 3
				for(int k=1;k<z-1;k++){
					float tmp=Rd/zdef[zstart-1+k];
					
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++)
					if(fdata[k+1][j][i][l]!=undef&&fdata[k-1][j][i][l]!=undef&&
					Tdata[k][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef&&Tdata[k+1][j][i][l]!=undef){
						float thetaU=(float)(Tdata[k+1][j][i][l]*pow(zdef[zstart-1+k+1]/100000,Rd/Cp));
						float thetaM=(float)(Tdata[k  ][j][i][l]*pow(zdef[zstart-1+k  ]/100000,Rd/Cp));
						float thetaD=(float)(Tdata[k-1][j][i][l]*pow(zdef[zstart-1+k-1]/100000,Rd/Cp));
						
						q3data[k][j][i][l]=f1[ystart-1+j]*(
							(fdata[k+1][j][i][l]-fdata[k][j][i][l])/
							(tmp*(Tdata[k+1][j][i][l]+Tdata[k][j][i][l])/2*(float)log(thetaU/thetaM))-
							(fdata[k][j][i][l]-fdata[k-1][j][i][l])/
							(tmp*(Tdata[k][j][i][l]+Tdata[k-1][j][i][l])/2*(float)log(thetaM/thetaD))
						)/((zdef[zstart-1+k]+zdef[zstart-2+k])/2);
						
					}else q3data[k][j][i][l]=undef;
				}
			}
		}
		
		return pv;
	}
	
	
	/**
     * calculate (Ertel) potential vorticity
     *
     * @param	u		u-wind (m/s)
     * @param	v		v-wind (m/s)
     * @param	w		omega (Pa/s)
     * @param	Th		(equivalent) potential temperature (K)
     *
     * @return	epv		(equivalent) Ertel potential vorticity (1 PVU = 10^-6 K m^2 Kg^-1 s^-1)
     */
	public Variable cErtelPV(Variable u,Variable v,Variable w,Variable Th){
		if(!u.isLike(v)||!u.isLike(w)||!u.isLike(Th)) throw new IllegalArgumentException("dimensions not same");
		
		ystart=u.getRange().getYRange()[0];
		zstart=u.getRange().getZRange()[0];
		
		t=u.getTCount();	z=u.getZCount();	y=u.getYCount();	x=u.getXCount();
		
		float undef=u.getUndef();
		
		Variable epv=new Variable("epv",u);
		epv.setUndef(undef);	epv.setCommentAndUnit("Ertel potential vorticity (1 PVU = 10^-6 K m^2 Kg^-1 s^-1)");
		
		float[][][][] udata=  u.getData();
		float[][][][] vdata=  v.getData();
		float[][][][] wdata=  w.getData();
		float[][][][] Tdata= Th.getData();
		float[][][][] edata=epv.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				if(vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef
				&& Tdata[l][k][j][i+1]!=undef&&Tdata[l][k][j][i-1]!=undef
				&& wdata[l][k][j+1][i]!=undef&&wdata[l][k][j-1][i]!=undef){
					edata[l][k][j][i]=(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(2*dxs[ystart-1+j])*(
						(wdata[l][k][j+1][i]-wdata[l][k][j-1][i])/(2*dy)+
						(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz)
					);
				}
				
				if(udata[l][k+1][j][i]!=undef&&udata[l][k-1][j][i]!=undef
				&& Tdata[l][k][j+1][i]!=undef&&Tdata[l][k][j-1][i]!=undef
				&& wdata[l][k][j][i+1]!=undef&&wdata[l][k][j][i-1]!=undef&&edata[l][k][j][i]!=undef){
					edata[l][k][j][i]+=(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(2*dy)*(
						f2[ystart-1+j]+
						(udata[l][k+1][j][i]-udata[l][k-1][j][i])/(dz+dz)+
						(wdata[l][k][j][i+1]-wdata[l][k][j][i-1])/(dxs[ystart-2+j]+dxs[ystart-1+j])
					);
				}
				
				if(udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef
				&& vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef
				&& Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef&&edata[l][k][j][i]!=undef){
					edata[l][k][j][i]+=(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz)*(
						f1[ystart-1+j]+
						(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(2*dxs[ystart-1+j])-
						(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(2*dy)
					);
				}
				
				if(edata[l][k][j][i]!=undef) edata[l][k][j][i]*=-GRAVITY_ACCERLERATION*1e6f;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				if(vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef
				&& Tdata[k][j][i+1][l]!=undef&&Tdata[k][j][i-1][l]!=undef
				&& wdata[k][j+1][i][l]!=undef&&wdata[k][j-1][i][l]!=undef){
					edata[k][j][i][l]=(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(2*dxs[ystart-1+j])*(
						(wdata[k][j+1][i][l]-wdata[k][j-1][i][l])/(2*dy)+
						(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz)
					);
				}
				
				if(udata[k+1][j][i][l]!=undef&&udata[k-1][j][i][l]!=undef
				&& Tdata[k][j+1][i][l]!=undef&&Tdata[k][j-1][i][l]!=undef
				&& wdata[k][j][i+1][l]!=undef&&wdata[k][j][i-1][l]!=undef&&edata[k][j][i][l]!=undef){
					edata[k][j][i][l]+=(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(2*dy)*(
						f2[ystart-1+j]+
						(udata[k+1][j][i][l]-udata[k-1][j][i][l])/(dz+dz)+
						(wdata[k][j][i+1][l]-wdata[k][j][i-1][l])/(dxs[ystart-2+j]+dxs[ystart-1+j])
					);
				}
				
				if(udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef
				&& vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef
				&& Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef&&edata[k][j][i][l]!=undef){
					edata[k][j][i][l]+=(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz)*(
						f1[ystart-1+j]+
						(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(2*dxs[ystart-1+j])-
						(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(2*dy)
					);
				}
				
				if(edata[k][j][i][l]!=undef) edata[k][j][i][l]*=-GRAVITY_ACCERLERATION*1e6f;
			}
		}
		
		return epv;
	}
	
	
	/** test
	public static void main(String[] arg){
		DiagnosisFactory df=DiagnosisFactory.parseFile("/lustre/home/qianyk/Data/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		SphericalSpacialModel ssm=new SphericalSpacialModel(dd);
		PotentialVorticityInSC pv=new PotentialVorticityInSC(ssm);
		ThermoMethodsInSC tm=new ThermoMethodsInSC(ssm);
		
		Range r=new Range("lon(90,200);lat(10,60);t(3,3)",dd);
		
		Variable[] vs=df.getVariables(r,"h","T","u","v","w");
		
		Variable fai=vs[0].multiply(9.8f);
		Variable T  =vs[1];
		Variable u  =vs[2];
		Variable v  =vs[3];
		Variable w  =vs[4];
		Variable th =tm.cPotentialTemperature(T);
		Variable sig1=tm.cStaticStabilityArgByPT(th);	sig1.setName("sig1");
		Variable sig2=tm.cStaticStabilityArgByT(T);		sig2.setName("sig2");
		
		Variable[] terms=pv.cQGPVTerms(fai,T);
		Variable epv=pv.cErtelPV(u,v,w,th);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"/lustre/home/qianyk/Data/QGPV.dat");
		dw.writeData(dd,fai,T,u,v,w,th,epv,sig1,sig2,terms[0],terms[1],terms[2]);	dw.closeFile();
	}*/
}
