/**
 * @(#)ThermoMethodsInCC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.application.EquationInCylindricalCoordinate;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static miniufo.geophysics.atmos.ThermoDynamics.*;


/**
 * calculation of heat quantity in spherical coordinate (parametric)
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ThermoMethodsInCC extends EquationInCylindricalCoordinate{
	
	/**
     * constructor
     *
     * @param	csm		initialized by a spacial model in cylindrical coordinate
     */
	public ThermoMethodsInCC(CylindricalSpatialModel csm){ super(csm);}
	
	
	/**
     * calculate latent heating due to condensation or evaporation, -L*dq/dt=Q, units: W/kg
     *
     * @param	q	specific humidity (kg/kg)
     * @param	u	uwind (m/s)
     * @param	v	vwind (m/s)
     * @param	w	omega (Pa/s)
     * @param	T	Temperature (K)
     *
     * @return	re	latent heating (W/kg)
     */
	public Variable cLatentHeating(Variable q,Variable u,Variable v,Variable w,Variable T){
		checkDimensions(q,u,v,w,T);
		assignSubDomainParams(q);
		
		Variable re=new Variable("lh",q);
		re.setCommentAndUnit("latent heating: -L*dq/dt (W kg^-1)");
		
		float[][][][]  qdata= q.getData();
		float[][][][]  udata= u.getData();
		float[][][][]  vdata= v.getData();
		float[][][][]  wdata= w.getData();
		float[][][][]  Tdata= T.getData();
		float[][][][] redata=re.getData();
		
		if(q.isTFirst()){
			/*** partial t ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(qdata[l+1][k][j][i]!=undef&&qdata[l-1][k][j][i]!=undef) redata[l][k][j][i]=(qdata[l+1][k][j][i]-qdata[l-1][k][j][i])/(dt*2);
				else redata[l][k][j][i]=undef;
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(qdata[1][k][j][i]!=undef&&qdata[0][k][j][i]!=undef) redata[0][k][j][i]=(qdata[1][k][j][i]-qdata[0][k][j][i])/dt;
				else redata[0][k][j][i]=undef;
				
				if(qdata[t-1][k][j][i]!=undef&&qdata[t-2][k][j][i]!=undef) redata[t-1][k][j][i]=(qdata[t-1][k][j][i]-qdata[t-2][k][j][i])/dt;
				else redata[t-1][k][j][i]=undef;
			}
			
			/*** partial x ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			for(int i=1;i<x-1;i++)
				if(udata[l][k][j][i]!=undef&&qdata[l][k][j][i+1]!=undef&&qdata[l][k][j][i-1]!=undef){
					if(redata[l][k][j][i]!=undef)
						redata[l][k][j][i]+=udata[l][k][j][i]*(qdata[l][k][j][i+1]-qdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
					
					else redata[l][k][j][i]=udata[l][k][j][i]*(qdata[l][k][j][i+1]-qdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
				}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++){
				if(udata[l][k][j][0]!=undef&&qdata[l][k][j][1]!=undef&&qdata[l][k][j][x-1]!=undef){
					if(redata[l][k][j][0]!=undef)
						redata[l][k][j][0]+=udata[l][k][j][0]*(qdata[l][k][j][1]-qdata[l][k][j][x-1])/(dxs[ystart-1+j]*2);
					
					else redata[l][k][j][0]=udata[l][k][j][0]*(qdata[l][k][j][1]-qdata[l][k][j][x-1])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[l][k][j][x-1]!=undef&&qdata[l][k][j][0]!=undef&&qdata[l][k][j][x-2]!=undef){
					if(redata[l][k][j][x-1]!=undef)
						redata[l][k][j][x-1]+=udata[l][k][j][x-1]*(qdata[l][k][j][0]-qdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
					
					else redata[l][k][j][x-1]=udata[l][k][j][x-1]*(qdata[l][k][j][0]-qdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
				}
			}
			
			/*** partial y ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			for(int j=1;j<y-1;j++)
				if(vdata[l][k][j][i]!=undef&&qdata[l][k][j+1][i]!=undef&&qdata[l][k][j-1][i]!=undef){
					if(redata[l][k][j][i]!=undef)
						redata[l][k][j][i]+=vdata[l][k][j][i]*(qdata[l][k][j+1][i]-qdata[l][k][j-1][i])/(dy*2);
					
					else redata[l][k][j][i]=vdata[l][k][j][i]*(qdata[l][k][j+1][i]-qdata[l][k][j-1][i])/(dy*2);
				}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(vdata[l][k][0][i]!=undef&&qdata[l][k][1][i]!=undef&&qdata[l][k][0][i]!=undef){
					if(redata[l][k][0][i]!=undef)
						redata[l][k][0][i]+=vdata[l][k][0][i]*(qdata[l][k][1][i]-qdata[l][k][0][i])/dy;
					
					else redata[l][k][0][i]=vdata[l][k][0][i]*(qdata[l][k][1][i]-qdata[l][k][0][i])/dy;
				}
				
				if(vdata[l][k][y-1][i]!=undef&&qdata[l][k][y-2][i]!=undef&&qdata[l][k][y-1][i]!=undef){
					if(redata[l][k][y-1][i]!=undef)
						redata[l][k][y-1][i]+=vdata[l][k][y-1][i]*(qdata[l][k][y-1][i]-qdata[l][k][y-2][i])/dy;
					
					else redata[l][k][y-1][i]=vdata[l][k][y-1][i]*(qdata[l][k][y-1][i]-qdata[l][k][y-2][i])/dy;
				}
			}
			
			/*** partial p ***/
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			for(int k=1;k<z-1;k++)
				if(wdata[l][k][j][i]!=undef&&qdata[l][k+1][j][i]!=undef&&qdata[l][k-1][j][i]!=undef){
					if(redata[l][k][j][i]!=undef)
						redata[l][k][j][i]+=wdata[l][k][j][i]*(qdata[l][k+1][j][i]-qdata[l][k-1][j][i])/(dz*2);
					
					else redata[l][k][j][i]=wdata[l][k][j][i]*(qdata[l][k+1][j][i]-qdata[l][k-1][j][i])/(dz*2);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[l][0][j][i]!=undef&&qdata[l][1][j][i]!=undef&&qdata[l][0][j][i]!=undef){
					if(redata[l][0][j][i]!=undef)
						redata[l][0][j][i]+=wdata[l][0][j][i]*(qdata[l][1][j][i]-qdata[l][0][j][i])/dz;
					
					else redata[l][0][j][i]=wdata[l][0][j][i]*(qdata[l][1][j][i]-qdata[l][0][j][i])/dz;
				}
				
				if(wdata[l][z-1][j][i]!=undef&&qdata[l][z-1][j][i]!=undef&&qdata[l][z-2][j][i]!=undef){
					if(redata[l][z-1][j][i]!=undef)
						redata[l][z-1][j][i]+=wdata[l][z-1][j][i]*(qdata[l][z-1][j][i]-qdata[l][z-2][j][i])/dz;
					
					else redata[l][z-1][j][i]=wdata[l][z-1][j][i]*(qdata[l][z-1][j][i]-qdata[l][z-2][j][i])/dz;
				}
			}
			
			/*** calculate latent heat ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)	// L=2500794-2370*t
				if(redata[l][k][j][i]!=undef&&Tdata[l][k][j][i]!=undef) redata[l][k][j][i]*=-cL(Tdata[l][k][j][i]);
			
		}else{
			/*** partial t ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(qdata[k][j][i][l+1]!=undef&&qdata[k][j][i][l-1]!=undef) redata[k][j][i][l]=(qdata[k][j][i][l+1]-qdata[k][j][i][l-1])/(dt*2);
				else redata[k][j][i][l]=undef;
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(qdata[k][j][i][1]!=undef&&qdata[k][j][i][0]!=undef) redata[k][j][i][0]=(qdata[k][j][i][1]-qdata[k][j][i][0])/dt;
				else redata[0][k][j][i]=undef;
				
				if(qdata[k][j][i][t-1]!=undef&&qdata[k][j][i][t-2]!=undef) redata[k][j][i][t-1]=(qdata[k][j][i][t-1]-qdata[k][j][i][t-2])/dt;
				else redata[k][j][i][t-1]=undef;
			}
			
			/*** partial x ***/
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
				if(udata[k][j][i][l]!=undef&&qdata[k][j][i+1][l]!=undef&&qdata[k][j][i-1][l]!=undef){
					if(redata[k][j][i][l]!=undef)
						redata[k][j][i][l]+=udata[k][j][i][l]*(qdata[k][j][i+1][l]-qdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
					
					else redata[k][j][i][l]=udata[k][j][i][l]*(qdata[k][j][i+1][l]-qdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
				}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++){
				if(udata[k][j][0][l]!=undef&&qdata[k][j][1][l]!=undef&&qdata[k][j][x-1][l]!=undef){
					if(redata[k][j][0][l]!=undef)
						redata[k][j][0][l]+=udata[k][j][0][l]*(qdata[k][j][1][l]-qdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
					
					else redata[k][j][0][l]=udata[k][j][0][l]*(qdata[k][j][1][l]-qdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[k][j][x-1][l]!=undef&&qdata[k][j][0][l]!=undef&&qdata[k][j][x-2][l]!=undef){
					if(redata[k][j][x-1][l]!=undef)
						redata[k][j][x-1][l]+=udata[k][j][x-1][l]*(qdata[k][j][0][l]-qdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
					
					else redata[k][j][x-1][l]=udata[k][j][x-1][l]*(qdata[k][j][0][l]-qdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
				}
			}
			
			/*** partial y ***/
			for(int j=1;j<y-1;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
				if(vdata[k][j][i][l]!=undef&&qdata[k][j+1][i][l]!=undef&&qdata[k][j-1][i][l]!=undef){
					if(redata[k][j][i][l]!=undef)
						redata[k][j][i][l]+=vdata[k][j][i][l]*(qdata[k][j+1][i][l]-qdata[k][j-1][i][l])/(dy*2);
					
					else redata[k][j][i][l]=vdata[k][j][i][l]*(qdata[k][j+1][i][l]-qdata[k][j-1][i][l])/(dy*2);
				}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(vdata[k][0][i][l]!=undef&&qdata[k][1][i][l]!=undef&&qdata[k][0][i][l]!=undef){
					if(redata[k][0][i][l]!=undef)
						redata[k][0][i][l]+=vdata[k][0][i][l]*(qdata[k][1][i][l]-qdata[k][0][i][l])/dy;
					
					else redata[k][0][i][l]=vdata[k][0][i][l]*(qdata[k][1][i][l]-qdata[k][0][i][l])/dy;
				}
				
				if(vdata[k][y-1][i][l]!=undef&&qdata[k][y-2][i][l]!=undef&&qdata[k][y-1][i][l]!=undef){
					if(redata[k][y-1][i][l]!=undef)
						redata[k][y-1][i][l]+=vdata[k][y-1][i][l]*(qdata[k][y-1][i][l]-qdata[k][y-2][i][l])/dy;
					
					else redata[k][y-1][i][l]=vdata[k][y-1][i][l]*(qdata[k][y-1][i][l]-qdata[k][y-2][i][l])/dy;
				}
			}
			
			/*** partial p ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&qdata[k+1][j][i][l]!=undef&&qdata[k-1][j][i][l]!=undef){
					if(redata[k][j][i][l]!=undef)
						redata[k][j][i][l]+=wdata[k][j][i][l]*(qdata[k+1][j][i][l]-qdata[k-1][j][i][l])/(dz*2);
					
					else redata[k][j][i][l]=wdata[k][j][i][l]*(qdata[k+1][j][i][l]-qdata[k-1][j][i][l])/(dz*2);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[0][j][i][l]!=undef&&qdata[1][j][i][l]!=undef&&qdata[0][j][i][l]!=undef){
					if(redata[0][j][i][l]!=undef)
						redata[0][j][i][l]+=wdata[0][j][i][l]*(qdata[1][j][i][l]-qdata[0][j][i][l])/dz;
					
					else redata[0][j][i][l]=wdata[0][j][i][l]*(qdata[1][j][i][l]-qdata[0][j][i][l])/dz;
				}
				
				if(wdata[z-1][j][i][l]!=undef&&qdata[z-1][j][i][l]!=undef&&qdata[z-2][j][i][l]!=undef){
					if(redata[z-1][j][i][l]!=undef)
						redata[z-1][j][i][l]+=wdata[z-1][j][i][l]*(qdata[z-1][j][i][l]-qdata[z-2][j][i][l])/dz;
					
					else redata[z-1][j][i][l]=wdata[z-1][j][i][l]*(qdata[z-1][j][i][l]-qdata[z-2][j][i][l])/dz;
				}
			}
			
			/*** calculate latent heat ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)	// L=2500794-2370*t
				if(redata[k][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef) redata[k][j][i][l]*=-cL(Tdata[k][j][i][l]);
		}
		
		return re;
	}
	
	/**
     * calculate diabatic heating rate, Cp*dT/dt+a*w=Q, units: W/kg
     *
     * @param	T	air temperature (K)
     * @param	u	u-wind (m/s)
     * @param	v	v-wind (m/s)
     * @param	w	w-wind (Pa/s)
     *
     * @return	Q	diabatic heating (W/kg)
     */
	public Variable cDiabaticHeatingRateByT(Variable T,Variable u,Variable v,Variable w){
		checkDimensions(u,v,w,T);
		assignSubDomainParams(u);
		
		Variable Q=new Variable("Qt",T);
		Q.setCommentAndUnit("diabatic heating rate: Q=Cp*dT/dt+a*omega (W kg^-1)");
		
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		float[][][][] wdata=w.getData();
		float[][][][] Tdata=T.getData();
		float[][][][] Qdata=Q.getData();
		
		if(T.isTFirst()){
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l+1][k][j][i]!=undef&&Tdata[l-1][k][j][i]!=undef) Qdata[l][k][j][i]=(Tdata[l+1][k][j][i]-Tdata[l-1][k][j][i])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[1  ][k][j][i]!=undef&&Tdata[0  ][k][j][i]!=undef) Qdata[0  ][k][j][i]=(Tdata[1  ][k][j][i]-Tdata[0  ][k][j][i])/dt;
				if(Tdata[t-1][k][j][i]!=undef&&Tdata[t-2][k][j][i]!=undef) Qdata[t-1][k][j][i]=(Tdata[t-1][k][j][i]-Tdata[t-2][k][j][i])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[l][k][j][i]!=undef&&Tdata[l][k][j][i+1]!=undef&&Tdata[l][k][j][i-1]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=udata[l][k][j][i]*(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
					
					else Qdata[l][k][j][i]=udata[l][k][j][i]*(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[l][k][j][i]!=undef&&Tdata[l][k][j+1][i]!=undef&&Tdata[l][k][j-1][i]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=vdata[l][k][j][i]*(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(dy*2);
					
					else Qdata[l][k][j][i]=vdata[l][k][j][i]*(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(dy*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[l][k][j][0]!=undef&&Tdata[l][k][j][1]!=undef&&Tdata[l][k][j][x-1]!=undef){
					if(Qdata[l][k][j][0]!=undef)
						Qdata[l][k][j][0]+=udata[l][k][j][0]*(Tdata[l][k][j][1]-Tdata[l][k][j][x-1])/dxs[ystart-1+j];
					
					else Qdata[l][k][j][0]=udata[l][k][j][0]*(Tdata[l][k][j][1]-Tdata[l][k][j][x-1])/dxs[ystart-1+j];
				}
				
				if(udata[l][k][j][x-1]!=undef&&Tdata[l][k][j][0]!=undef&&Tdata[l][k][j][x-2]!=undef){
					if(Qdata[l][k][j][x-1]!=undef)
						Qdata[l][k][j][x-1]+=vdata[l][k][j][x-1]*(Tdata[l][k][j][0]-Tdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
					
					else Qdata[l][k][j][x-1]=vdata[l][k][j][x-1]*(Tdata[l][k][j][0]-Tdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[l][k][0][i]!=undef&&Tdata[l][k][0][i]!=undef&&Tdata[l][k][1][i]!=undef){
					if(Qdata[l][k][0][i]!=undef)
						Qdata[l][k][0][i]+=udata[l][k][0][i]*(Tdata[l][k][1][i]-Tdata[l][k][0][i])/dy;
					
					else Qdata[l][k][0][i]=udata[l][k][0][i]*(Tdata[l][k][1][i]-Tdata[l][k][0][i])/dy;
				}
				
				if(udata[l][k][y-1][i]!=undef&&Tdata[l][k][y-1][i]!=undef&&Tdata[l][k][y-2][i]!=undef){
					if(Qdata[l][k][y-1][i]!=undef)
						Qdata[l][k][y-1][i]+=vdata[l][k][y-1][i]*(Tdata[l][k][y-1][i]-Tdata[l][k][y-2][i])/dy;
					
					else Qdata[l][k][y-1][i]=vdata[l][k][y-1][i]*(Tdata[l][k][y-1][i]-Tdata[l][k][y-2][i])/dy;
				}
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=wdata[l][k][j][i]*(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz*2);
					
					else Qdata[l][k][j][i]=wdata[l][k][j][i]*(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz*2);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef&&Tdata[l][0][j][i]!=undef){
					if(Qdata[l][0][j][i]!=undef)
						Qdata[l][0][j][i]+=wdata[l][0][j][i]*(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
					
					else Qdata[l][0][j][i]=wdata[l][0][j][i]*(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
				}
				
				if(wdata[l][z-1][j][i]!=undef&&Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef){
					if(Qdata[l][z-1][j][i]!=undef)
						Qdata[l][z-1][j][i]+=wdata[l][z-1][j][i]*(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
					
					else Qdata[l][z-1][j][i]=wdata[l][z-1][j][i]*(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
				}
			}
			
			/*** a*w ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&Tdata[l][k][j][i]!=undef){
					if(Qdata[l][k][j][i]!=undef){
						Qdata[l][k][j][i]*=Cp;
						Qdata[l][k][j][i]-=wdata[l][k][j][i]*Tdata[l][k][j][i]*Rd/zdef[zstart-1+k];
						
					}else Qdata[l][k][j][i]=-wdata[l][k][j][i]*Tdata[l][k][j][i]*Rd/zdef[zstart-1+k];
				}
			
		}else{
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l+1]!=undef&&Tdata[k][j][i][l-1]!=undef) Qdata[k][j][i][l]=(Tdata[k][j][i][l+1]-Tdata[k][j][i][l-1])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[k][j][i][1  ]!=undef&&Tdata[k][j][i][0  ]!=undef) Qdata[k][j][i][0  ]=(Tdata[k][j][i][1  ]-Tdata[k][j][i][0  ])/dt;
				if(Tdata[k][j][i][t-1]!=undef&&Tdata[k][j][i][t-2]!=undef) Qdata[k][j][i][t-1]=(Tdata[k][j][i][t-1]-Tdata[k][j][i][t-2])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[k][j][i][l]!=undef&&Tdata[k][j][i+1][l]!=undef&&Tdata[k][j][i-1][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=udata[k][j][i][l]*(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
					
					else Qdata[k][j][i][l]=udata[k][j][i][l]*(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[k][j][i][l]!=undef&&Tdata[k][j+1][i][l]!=undef&&Tdata[k][j-1][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=vdata[k][j][i][l]*(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(dy*2);
					
					else Qdata[k][j][i][l]=vdata[k][j][i][l]*(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(dy*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[k][j][0][l]!=undef&&Tdata[k][j][1][l]!=undef&&Tdata[k][j][x-1][l]!=undef){
					if(Qdata[k][j][0][l]!=undef)
						Qdata[k][j][0][l]+=udata[k][j][0][l]*(Tdata[k][j][1][l]-Tdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
					
					else Qdata[k][j][0][l]=udata[k][j][0][l]*(Tdata[k][j][1][l]-Tdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[k][j][x-1][l]!=undef&&Tdata[k][j][0][l]!=undef&&Tdata[k][j][x-2][l]!=undef){
					if(Qdata[k][j][x-1][l]!=undef)
						Qdata[k][j][x-1][l]+=vdata[k][j][x-1][l]*(Tdata[k][j][0][l]-Tdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
					
					else Qdata[k][j][x-1][l]=vdata[k][j][x-1][l]*(Tdata[k][j][0][l]-Tdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[k][0][i][l]!=undef&&Tdata[k][0][i][l]!=undef&&Tdata[k][1][i][l]!=undef){
					if(Qdata[k][0][i][l]!=undef)
						Qdata[k][0][i][l]+=udata[k][0][i][l]*(Tdata[k][1][i][l]-Tdata[k][0][i][l])/dy;
					
					else Qdata[k][0][i][l]=udata[k][0][i][l]*(Tdata[k][1][i][l]-Tdata[k][0][i][l])/dy;
				}
				
				if(udata[k][y-1][i][l]!=undef&&Tdata[k][y-1][i][l]!=undef&&Tdata[k][y-2][i][l]!=undef){
					if(Qdata[k][y-1][i][l]!=undef)
						Qdata[k][y-1][i][l]+=vdata[k][y-1][i][l]*(Tdata[k][y-1][i][l]-Tdata[k][y-2][i][l])/dy;
					
					else Qdata[k][y-1][i][l]=vdata[k][y-1][i][l]*(Tdata[k][y-1][i][l]-Tdata[k][y-2][i][l])/dy;
				}
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=wdata[k][j][i][l]*(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz*2);
					
					else Qdata[k][j][i][l]=wdata[k][j][i][l]*(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz*2);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef&&Tdata[0][j][i][l]!=undef){
					if(Qdata[0][j][i][l]!=undef)
						Qdata[0][j][i][l]+=wdata[0][j][i][l]*(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
					
					else Qdata[0][j][i][l]=wdata[0][j][i][l]*(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
				}
				
				if(wdata[z-1][j][i][l]!=undef&&Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef){
					if(Qdata[z-1][j][i][l]!=undef)
						Qdata[z-1][j][i][l]+=wdata[z-1][j][i][l]*(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
					
					else Qdata[z-1][j][i][l]=wdata[z-1][j][i][l]*(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
				}
			}
			
			/*** a*w ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&Tdata[k][j][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef){
						Qdata[k][j][i][l]*=Cp;
						Qdata[k][j][i][l]-=wdata[k][j][i][l]*Tdata[k][j][i][l]*Rd/zdef[zstart-1+k];
						
					}else Qdata[k][j][i][l]=-wdata[k][j][i][l]*Tdata[k][j][i][l]*Rd/zdef[zstart-1+k];
				}
		}
		
		return Q;
	}
	
	/**
     * calculate diabatic heating rate, ¦Ð*Cp*d¦È/dt=Q, units: W/kg
     *
     * @param	Th	potential temperature (K)
     * @param	u	u-wind (m/s)
     * @param	v	v-wind (m/s)
     * @param	w	w-wind (Pa/s)
     *
     * @return	Q	diabatic heating (W/kg)
     */
	public Variable cDiabaticHeatingRateByPT(Variable Th,Variable u,Variable v,Variable w){
		checkDimensions(u,v,w,Th);
		assignSubDomainParams(u);
		
		Variable Q=new Variable("Qpt",Th);
		Q.setCommentAndUnit("diabatic heating rate: Q=¦Ð*Cp*d¦È/dt+a*omega (W kg^-1)");
		
		float[][][][] udata= u.getData();
		float[][][][] vdata= v.getData();
		float[][][][] wdata= w.getData();
		float[][][][] Tdata=Th.getData();
		float[][][][] Qdata= Q.getData();
		
		if(Th.isTFirst()){
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l+1][k][j][i]!=undef&&Tdata[l-1][k][j][i]!=undef) Qdata[l][k][j][i]=(Tdata[l+1][k][j][i]-Tdata[l-1][k][j][i])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[1  ][k][j][i]!=undef&&Tdata[0  ][k][j][i]!=undef) Qdata[0  ][k][j][i]=(Tdata[1  ][k][j][i]-Tdata[0  ][k][j][i])/dt;
				if(Tdata[t-1][k][j][i]!=undef&&Tdata[t-2][k][j][i]!=undef) Qdata[t-1][k][j][i]=(Tdata[t-1][k][j][i]-Tdata[t-2][k][j][i])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[l][k][j][i]!=undef&&Tdata[l][k][j][i+1]!=undef&&Tdata[l][k][j][i-1]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=udata[l][k][j][i]*(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
					
					else Qdata[l][k][j][i]=udata[l][k][j][i]*(Tdata[l][k][j][i+1]-Tdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[l][k][j][i]!=undef&&Tdata[l][k][j+1][i]!=undef&&Tdata[l][k][j-1][i]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=vdata[l][k][j][i]*(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(dy*2);
					
					else Qdata[l][k][j][i]=vdata[l][k][j][i]*(Tdata[l][k][j+1][i]-Tdata[l][k][j-1][i])/(dy*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[l][k][j][0]!=undef&&Tdata[l][k][j][1]!=undef&&Tdata[l][k][j][x-1]!=undef){
					if(Qdata[l][k][j][0]!=undef)
						Qdata[l][k][j][0]+=udata[l][k][j][0]*(Tdata[l][k][j][1]-Tdata[l][k][j][x-1])/dxs[ystart-1+j];
					
					else Qdata[l][k][j][0]=udata[l][k][j][0]*(Tdata[l][k][j][1]-Tdata[l][k][j][x-1])/dxs[ystart-1+j];
				}
				
				if(udata[l][k][j][x-1]!=undef&&Tdata[l][k][j][0]!=undef&&Tdata[l][k][j][x-2]!=undef){
					if(Qdata[l][k][j][x-1]!=undef)
						Qdata[l][k][j][x-1]+=vdata[l][k][j][x-1]*(Tdata[l][k][j][0]-Tdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
					
					else Qdata[l][k][j][x-1]=vdata[l][k][j][x-1]*(Tdata[l][k][j][0]-Tdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[l][k][0][i]!=undef&&Tdata[l][k][0][i]!=undef&&Tdata[l][k][1][i]!=undef){
					if(Qdata[l][k][0][i]!=undef)
						Qdata[l][k][0][i]+=udata[l][k][0][i]*(Tdata[l][k][1][i]-Tdata[l][k][0][i])/dy;
					
					else Qdata[l][k][0][i]=udata[l][k][0][i]*(Tdata[l][k][1][i]-Tdata[l][k][0][i])/dy;
				}
				
				if(udata[l][k][y-1][i]!=undef&&Tdata[l][k][y-1][i]!=undef&&Tdata[l][k][y-2][i]!=undef){
					if(Qdata[l][k][y-1][i]!=undef)
						Qdata[l][k][y-1][i]+=vdata[l][k][y-1][i]*(Tdata[l][k][y-1][i]-Tdata[l][k][y-2][i])/dy;
					
					else Qdata[l][k][y-1][i]=vdata[l][k][y-1][i]*(Tdata[l][k][y-1][i]-Tdata[l][k][y-2][i])/dy;
				}
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef){
					if(Qdata[l][k][j][i]!=undef)
						Qdata[l][k][j][i]+=wdata[l][k][j][i]*(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz*2);
					
					else Qdata[l][k][j][i]=wdata[l][k][j][i]*(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz*2);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef&&Tdata[l][0][j][i]!=undef){
					if(Qdata[l][0][j][i]!=undef)
						Qdata[l][0][j][i]+=wdata[l][0][j][i]*(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
					
					else Qdata[l][0][j][i]=wdata[l][0][j][i]*(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
				}
				
				if(wdata[l][z-1][j][i]!=undef&&Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef){
					if(Qdata[l][z-1][j][i]!=undef)
						Qdata[l][z-1][j][i]+=wdata[l][z-1][j][i]*(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
					
					else Qdata[l][z-1][j][i]=wdata[l][z-1][j][i]*(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
				}
			}
			
			/*** scaled by ¦ÐCp ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Qdata[l][k][j][i]!=undef) Qdata[l][k][j][i]*=Cp*Math.pow(zdef[k]/100000,kp);
					
			
		}else{
			/*** local partial term ***/
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l+1]!=undef&&Tdata[k][j][i][l-1]!=undef) Qdata[k][j][i][l]=(Tdata[k][j][i][l+1]-Tdata[k][j][i][l-1])/(dt*2);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[k][j][i][1  ]!=undef&&Tdata[k][j][i][0  ]!=undef) Qdata[k][j][i][0  ]=(Tdata[k][j][i][1  ]-Tdata[k][j][i][0  ])/dt;
				if(Tdata[k][j][i][t-1]!=undef&&Tdata[k][j][i][t-2]!=undef) Qdata[k][j][i][t-1]=(Tdata[k][j][i][t-1]-Tdata[k][j][i][t-2])/dt;
			}
			
			/*** advectiion term ***/
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				if(udata[k][j][i][l]!=undef&&Tdata[k][j][i+1][l]!=undef&&Tdata[k][j][i-1][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=udata[k][j][i][l]*(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
					
					else Qdata[k][j][i][l]=udata[k][j][i][l]*(Tdata[k][j][i+1][l]-Tdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[k][j][i][l]!=undef&&Tdata[k][j+1][i][l]!=undef&&Tdata[k][j-1][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=vdata[k][j][i][l]*(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(dy*2);
					
					else Qdata[k][j][i][l]=vdata[k][j][i][l]*(Tdata[k][j+1][i][l]-Tdata[k][j-1][i][l])/(dy*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(udata[k][j][0][l]!=undef&&Tdata[k][j][1][l]!=undef&&Tdata[k][j][x-1][l]!=undef){
					if(Qdata[k][j][0][l]!=undef)
						Qdata[k][j][0][l]+=udata[k][j][0][l]*(Tdata[k][j][1][l]-Tdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
					
					else Qdata[k][j][0][l]=udata[k][j][0][l]*(Tdata[k][j][1][l]-Tdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
				}
				
				if(udata[k][j][x-1][l]!=undef&&Tdata[k][j][0][l]!=undef&&Tdata[k][j][x-2][l]!=undef){
					if(Qdata[k][j][x-1][l]!=undef)
						Qdata[k][j][x-1][l]+=vdata[k][j][x-1][l]*(Tdata[k][j][0][l]-Tdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
					
					else Qdata[k][j][x-1][l]=vdata[k][j][x-1][l]*(Tdata[k][j][0][l]-Tdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
				}
			}
			
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++){
				if(udata[k][0][i][l]!=undef&&Tdata[k][0][i][l]!=undef&&Tdata[k][1][i][l]!=undef){
					if(Qdata[k][0][i][l]!=undef)
						Qdata[k][0][i][l]+=udata[k][0][i][l]*(Tdata[k][1][i][l]-Tdata[k][0][i][l])/dy;
					
					else Qdata[k][0][i][l]=udata[k][0][i][l]*(Tdata[k][1][i][l]-Tdata[k][0][i][l])/dy;
				}
				
				if(udata[k][y-1][i][l]!=undef&&Tdata[k][y-1][i][l]!=undef&&Tdata[k][y-2][i][l]!=undef){
					if(Qdata[k][y-1][i][l]!=undef)
						Qdata[k][y-1][i][l]+=vdata[k][y-1][i][l]*(Tdata[k][y-1][i][l]-Tdata[k][y-2][i][l])/dy;
					
					else Qdata[k][y-1][i][l]=vdata[k][y-1][i][l]*(Tdata[k][y-1][i][l]-Tdata[k][y-2][i][l])/dy;
				}
			}
			
			/*** convective term ***/
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef){
					if(Qdata[k][j][i][l]!=undef)
						Qdata[k][j][i][l]+=wdata[k][j][i][l]*(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz*2);
					
					else Qdata[k][j][i][l]=wdata[k][j][i][l]*(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz*2);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef&&Tdata[0][j][i][l]!=undef){
					if(Qdata[0][j][i][l]!=undef)
						Qdata[0][j][i][l]+=wdata[0][j][i][l]*(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
					
					else Qdata[0][j][i][l]=wdata[0][j][i][l]*(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
				}
				
				if(wdata[z-1][j][i][l]!=undef&&Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef){
					if(Qdata[z-1][j][i][l]!=undef)
						Qdata[z-1][j][i][l]+=wdata[z-1][j][i][l]*(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
					
					else Qdata[z-1][j][i][l]=wdata[z-1][j][i][l]*(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
				}
			}
			
			/*** scaled by ¦ÐCp ***/
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Qdata[k][j][i][l]!=undef) Qdata[k][j][i][l]*=Cp*Math.pow(zdef[k]/100000,kp);
		}
		
		return Q;
	}
	
	
	/**
     * calculate static stability argument
     *
     * @param	T		temperature of air (K)
     *
     * @return	sigma	static stability argument (m^4 s^2 kg^-1)
     */
	public Variable cStaticStabilityArgByT(Variable T){
		assignSubDomainParams(T);
		
		Variable sigma=new Variable("sigmaT",T);
		sigma.setCommentAndUnit("static stability argument (m^4 s^2 kg^-1)");
		
		float[][][][] 	  Tdata=	T.getData();
		float[][][][] sigmadata=sigma.getData();
		
		if(T.isTFirst()){
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef)
						sigmadata[l][k][j][i]=(
							tmp*Tdata[l][k][j][i]/Cp-(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz)
						)*tmp;
						
					else sigmadata[l][k][j][i]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef)
					sigmadata[l][0][j][i]=(
						tmp*Tdata[l][0][j][i]/Cp-(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz
					)*tmp;
					
				else sigmadata[l][0][j][i]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef)
					sigmadata[l][z-1][j][i]=(
						tmp*Tdata[l][z-1][j][i]/Cp-(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz
					)*tmp;
					
				else sigmadata[l][z-1][j][i]=undef;
			}
			
		}else{
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef)
						sigmadata[k][j][i][l]=(
							tmp*Tdata[k][j][i][l]/Cp-(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz)
						)*tmp;
						
					else sigmadata[k][j][i][l]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef)
					sigmadata[0][j][i][l]=(
						tmp*Tdata[0][j][i][l]/Cp-(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz
					)*tmp;
					
				else sigmadata[0][j][i][l]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef)
					sigmadata[z-1][j][i][l]=(
						tmp*Tdata[z-1][j][i][l]/Cp-(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz
					)*tmp;
					
				else sigmadata[z-1][j][i][l]=undef;
			}
		}
		
		return sigma;
	}
	
    /**
     * calculate static stability argument
     *
     * @param	th		potential temperature (K)
     *
     * @return	sigma	static stability argument (m^4 s^2 kg^-1)
     */
    public Variable cStaticStabilityArgByPT(Variable th){
		assignSubDomainParams(th);
		
		Variable sigma=new Variable("sigmaPT",th);
		sigma.setCommentAndUnit("static stability argument (m^4 s^2 kg^-1)");
		
		float[][][][] 	  Tdata=   th.getData();
		float[][][][] sigmadata=sigma.getData();
		
		if(th.isTFirst()){
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k+1][j][i]!=undef&&Tdata[l][k-1][j][i]!=undef)
						sigmadata[l][k][j][i]=
							-tmp*(float)pow(zdef[zstart-1+k]/100000f,kp)*
							(Tdata[l][k+1][j][i]-Tdata[l][k-1][j][i])/(dz+dz);
						
					else sigmadata[l][k][j][i]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][0][j][i]!=undef&&Tdata[l][1][j][i]!=undef)
					sigmadata[l][0][j][i]=
						-tmp*(float)pow(zdef[zstart-1]/100000f,kp)*
						(Tdata[l][1][j][i]-Tdata[l][0][j][i])/dz;
					
				else sigmadata[l][0][j][i]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[l][z-1][j][i]!=undef&&Tdata[l][z-2][j][i]!=undef)
					sigmadata[l][z-1][j][i]=
						-tmp*(float)pow(zdef[zstart-2+z]/100000f,kp)*
						(Tdata[l][z-1][j][i]-Tdata[l][z-2][j][i])/dz;
					
				else sigmadata[l][z-1][j][i]=undef;
			}
			
		}else{
			for(int k=1;k<z-1;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k+1][j][i][l]!=undef&&Tdata[k-1][j][i][l]!=undef)
						sigmadata[k][j][i][l]=
							-tmp*(float)pow(zdef[zstart-1+k]/100000f,kp)*
							(Tdata[k+1][j][i][l]-Tdata[k-1][j][i][l])/(dz+dz);
						
					else sigmadata[k][j][i][l]=undef;
				}
			}
			
			/*** k==0 ***/
			float tmp=Rd/zdef[zstart-1];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[0][j][i][l]!=undef&&Tdata[1][j][i][l]!=undef)
					sigmadata[0][j][i][l]=
						-tmp*(float)pow(zdef[zstart-1]/100000f,kp)*
						(Tdata[1][j][i][l]-Tdata[0][j][i][l])/dz;
					
				else sigmadata[0][j][i][l]=undef;
			}
			
			/*** k==z-1 ***/
			tmp=Rd/zdef[zstart-2+z];
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(Tdata[z-1][j][i][l]!=undef&&Tdata[z-2][j][i][l]!=undef)
					sigmadata[z-1][j][i][l]=
						-tmp*(float)pow(zdef[zstart-2+z]/100000f,kp)*
						(Tdata[z-1][j][i][l]-Tdata[z-2][j][i][l])/dz;
					
				else sigmadata[z-1][j][i][l]=undef;
			}
		}
		
		return sigma;
     }
	
    
	/**
     * calculate temperature at lifting condensation level
     *
     * @param	T	air temperature (K)
     * @param	rh	relative humidity (%)
     *
     * @return	Tc	temperature at lifting condensation level (K)
     */
    public Variable cLCLTemperature(Variable T,Variable rh){
		checkDimensions(T,rh);
		assignSubDomainParams(T);
		
		Variable Tc=new Variable("Tc",T);
		Tc.setCommentAndUnit("temperature at lifting condensation level (K)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][] rhdata=rh.getData();
		float[][][][] Tcdata=Tc.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&rhdata[l][k][j][i]!=undef){
					//tmp=(float)(exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i]))
					//*rhdata[l][k][j][i]/100);
					
					//Tcdata[l][k][j][i]=55f+2840f/(float)(3.5*log(Tdata[l][k][j][i])-log(tmp)-4.805);
					
					if(rhdata[l][k][j][i]<=0) rhdata[l][k][j][i]=0.1f;
					
					Tcdata[l][k][j][i]=55+2840/(float)(
					8.3451*log(Tdata[l][k][j][i])-53.67957+6743.769/Tdata[l][k][j][i]-log(rhdata[l][k][j][i]/100)-4.805
					);
					
				}else Tcdata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&rhdata[k][j][i][l]!=undef){
					if(rhdata[k][j][i][l]<=0) rhdata[k][j][i][l]=0.1f;
					
					Tcdata[k][j][i][l]=55+2840/(float)(
					8.3451*log(Tdata[k][j][i][l])-53.67957+6743.769/Tdata[k][j][i][l]-log(rhdata[k][j][i][l]/100)-4.805
					);
					
				}else Tcdata[k][j][i][l]=undef;
		}
		
		return Tc;
	}
    
	/**
     * calculate potential temperature
     *
     * @param	T		temperature of air (K)
     *
     * @return	theta	potential temperature (K)
     */
    public Variable cPotentialTemperature(Variable T){
		assignSubDomainParams(T);
		
		Variable theta=new Variable("th",T);
		theta.setCommentAndUnit("potential temperature (K)");
		
		float[][][][] 	  Tdata=	T.getData();
		float[][][][] thetadata=theta.getData();
		
		if(T.isTFirst()){
			for(int k=0;k<z;k++){
				float tmp=(float)pow((100000/zdef[zstart-1+k]),(Rd/Cp));
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k][j][i]!=undef) thetadata[l][k][j][i]=tmp*Tdata[l][k][j][i];
					else thetadata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int k=0;k<z;k++){
				float tmp=(float)pow((100000/zdef[zstart-1+k]),(Rd/Cp));
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k][j][i][l]!=undef) thetadata[k][j][i][l]=tmp*Tdata[k][j][i][l];
					else thetadata[k][j][i][l]=undef;
				}
			}
		}
		
		Tdata=null;	thetadata=null;
		
		return theta;
	}
    
	/**
     * calculate equivalent potential temperature
     * according to Bolton (1980)
     *
     * @param	T	air temperature (K)
     * @param	q	specific humidity or mixing ratio (Kg/Kg)
     * @param	Tc	temperature at lifting condensation level (K)
     *
     * @return	re	equivalent potential temperature (K)
     */
	public Variable cEquivalentPotentialTemperature(Variable T,Variable q,Variable Tc){
		checkDimensions(T,q,Tc);
		assignSubDomainParams(T);
		
		Variable re=new Variable("The",T);
		re.setCommentAndUnit("equivalent potential temperature (K)");
		
		float[][][][]  qdata= q.getData();
		float[][][][]  Tdata= T.getData();
		float[][][][] Tcdata=Tc.getData();
		float[][][][] redata=re.getData();
		
		if(T.isTFirst()){
			for(int k=0;k<z;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Tdata[l][k][j][i]!=undef&&qdata[l][k][j][i]!=undef&&Tcdata[l][k][j][i]!=undef)
				redata[l][k][j][i]=(float)(Tdata[l][k][j][i]*
				pow(100000f/zdef[zstart-1+k],Rd/Cp*(1-0.28f*qdata[l][k][j][i]))*exp(
					(3376f/Tcdata[l][k][j][i]-2.54f)*qdata[l][k][j][i]*(1+0.81f*qdata[l][k][j][i])
				));
				
			else redata[l][k][j][i]=undef;
			
		}else{
			for(int k=0;k<z;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(Tdata[k][j][i][l]!=undef&&qdata[k][j][i][l]!=undef&&Tcdata[k][j][i][l]!=undef)
				redata[k][j][i][l]=(float)(Tdata[k][j][i][l]*
				pow(100000/zdef[zstart-1+k],Rd/Cp*(1-0.28f*qdata[k][j][i][l]))*exp(
					(3376f/Tcdata[k][j][i][l]-2.54f)*qdata[k][j][i][l]*(1+0.81f*qdata[k][j][i][l])
				));
				
			else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
    
	
	/**
     * calculate specific humidity using relative humidity and temperature
     *
     * @param	T	temperature (K)
     * @param	rh	relative humidity (%)
     *
     * @return	re	specific humidity (Kg/Kg)
     */
	public Variable cSpecificHumidity(Variable T,Variable rh){
		checkDimensions(T,rh);
		assignSubDomainParams(T);
		
		Variable re=new Variable("spfh",T);
		re.setCommentAndUnit("specific humidity or mixing ratio (Kg Kg^-1)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][] rhdata=rh.getData();
		float[][][][] redata=re.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&rhdata[l][k][j][i]!=undef)
					redata[l][k][j][i]=Rd/Rv*
						(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])))*
						rhdata[l][k][j][i]/zdef[zstart-1+k]/100;
					
				else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&rhdata[k][j][i][l]!=undef)
					redata[k][j][i][l]=Rd/Rv*
						(float)(100*exp(53.67957-6743.769/Tdata[k][j][i][l]-4.8451*log(Tdata[k][j][i][l])))*
						rhdata[k][j][i][l]/zdef[zstart-1+k]/100;
					
				else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
	
	/**
     * calculate relative humidity using specific humidity and temperature
     *
     * @param	T	temperature (K)
     * @param	q	specific humidity (Kg/Kg)
     *
     * @return	re	relative humidity (%)
     */
	public Variable cRelativeHumidity(Variable T,Variable q){
		checkDimensions(T,q);
		assignSubDomainParams(T);
		
		Variable re=new Variable("rh",T);
		re.setCommentAndUnit("relative humidity (%)");
		
		float[][][][]  Tdata= T.getData();
		float[][][][]  qdata= q.getData();
		float[][][][] redata=re.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[l][k][j][i]!=undef&&qdata[l][k][j][i]!=undef)
					redata[l][k][j][i]=Rv/Rd/
						(float)(100*exp(53.67957-6743.769/Tdata[l][k][j][i]-4.8451*log(Tdata[l][k][j][i])))*
						qdata[l][k][j][i]*zdef[zstart-1+k]*100;
					
				else redata[l][k][j][i]=undef;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(Tdata[k][j][i][l]!=undef&&qdata[k][j][i][l]!=undef)
					redata[k][j][i][l]=Rv/Rd/
						(float)(100*exp(53.67957-6743.769/Tdata[k][j][i][l]-4.8451*log(Tdata[k][j][i][l])))*
						qdata[k][j][i][l]*zdef[zstart-1+k]*100;
					
				else redata[k][j][i][l]=undef;
		}
		
		return re;
	}
    
	/**
     * calculate specific volume according to the equation of state
     *
     * @param	T	temperature of air (K)
     *
     * @return	a	specific volume (m^3/kg)
     */
    public Variable cSpecificVolume(Variable T){
		assignSubDomainParams(T);
		
		Variable a=new Variable("a",T);
		a.setCommentAndUnit("specific volume (m^3 Kg^-1)");
		
		float[][][][] Tdata=T.getData();
		float[][][][] adata=a.getData();
		
		if(T.isTFirst()){
			for(int k=0;k<z;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[l][k][j][i]!=undef) adata[l][k][j][i]=tmp*Tdata[l][k][j][i];
					else adata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int k=0;k<z;k++){
				float tmp=Rd/zdef[zstart-1+k];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(Tdata[k][j][i][l]!=undef) adata[k][j][i][l]=tmp*Tdata[k][j][i][l];
					else adata[k][j][i][l]=undef;
				}
			}
		}
		
		return a;
     }
    
    
	/**
     * calculate the maximum intensity of a storm according to
     * the underlying SST (Merrill 1988, JAS; DeMaria et al. 1993, JAS)
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity (m/s)
     */
	public Variable cMaximumPotentialIntensity1(Variable sst){
		assignSubDomainParams(sst);
		
		Variable mpi=new Variable("mpi",sst);
		mpi.setCommentAndUnit("maximum potential intensity (m s^-1), (Merrill 1988, JAS; DeMaria et al. 1993, JAS)");
		
		float[][][][] sdata=sst.getData();
		float[][][][] mdata=mpi.getData();
		
		if(sst.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(sdata[l][k][j][i]!=undef){
					mdata[l][k][j][i]=(float)(0.51444*74.0*exp(0.2*(sdata[l][k][j][i]-273.15-25.0)));
					
				}else mdata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(sdata[k][j][i][l]!=undef){
					mdata[l][k][j][i]=(float)(0.51444*74.0*exp(0.2*(sdata[k][j][i][l]-273.15-25.0)));
					
				}else mdata[k][j][i][l]=undef;
			}
		}
		
		return mpi;
	}
	
	/**
     * calculate the maximum intensity of a storm according to
     * the underlying SST over the Atlantic Ocean basin
     * (DeMaria and Kaplan 1994, JC)
     *
     * @param	sst		underlying sea surface temperature (K)
     *
     * @return	mpi		maximum potential intensity (m/s)
     */
	public Variable cMaximumPotentialIntensity2(Variable sst){
		assignSubDomainParams(sst);
		
		Variable mpi=new Variable("mpi",sst);
		mpi.setCommentAndUnit(
			"maximum potential intensity (m s^-1) over Atlantic Ocean basin,"+
			" (DeMaria and Kaplan 1994, JC)"
		);
		
		float[][][][] sdata=sst.getData();
		float[][][][] mdata=mpi.getData();
		
		if(sst.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(sdata[l][k][j][i]!=undef){
					mdata[l][k][j][i]=(float)(28.2+55.8*exp(0.1813*(sdata[l][k][j][i]-273.15-30)));
					
				}else mdata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(sdata[k][j][i][l]!=undef){
					mdata[k][j][i][l]=(float)(28.2+55.8*exp(0.1813*(sdata[k][j][i][l]-273.15-30)));
					
				}else mdata[k][j][i][l]=undef;
			}
		}
		
		return mpi;
	}
    
    
	/** test
	public static void main(String[] args){
		for(int t=100;t>=-60;t-=16){
			System.out.println(t+"\t"+cL(t));
		}
	}*/
}
