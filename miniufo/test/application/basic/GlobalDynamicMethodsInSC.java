/**
 * @(#)GlobalDynamicMethodsInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


/**
 * basic analysis methods in global spheral coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GlobalDynamicMethodsInSC extends DynamicMethodsInSC{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public GlobalDynamicMethodsInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isZonalPeriodic())
		throw new IllegalArgumentException("Not a zonal periodic model");
	}
	
	
	/**
     * calculate divergence under the condition of west-east periodic boundary
     *
     * @param	u	u-wind
     * @param	v	v-wind
     *
     * @return	divergence
     *
     * @exception	if u,v are not dimensionlly the same
     */
	public Variable c2DDivergence(Variable u,Variable v){
		if(!u.isLike(v)) throw new IllegalArgumentException("dimensions not same");
		
		t=u.getTCount();	z=u.getZCount();	y=u.getYCount();	x=u.getXCount();
		
		float undef=u.getUndef();
		Variable div=new Variable("div",u);	div.setCommentAndUnit("divergence (s^-1)");
		
		float[][][][]   udata=u.getData();
		float[][][][]   vdata=v.getData();
		float[][][][] divdata=div.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					/********************* div *********************/
					if(udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
						divdata[l][k][j][i]=
							(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[j]*2)+
							(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
							vdata[l][k][j][i]*ltan[j]/EARTH_RADIUS;
						
					else divdata[l][k][j][i]=undef;
				}
				
				/*** left and right boundry of div ***/
				if(udata[l][k][j][1]!=undef&&udata[l][k][j][x-1]!=undef&&vdata[l][k][j+1][0]!=undef&&vdata[l][k][j-1][0]!=undef)
					divdata[l][k][j][0]=
						(udata[l][k][j  ][1]-udata[l][k][j  ][x-1])/(dxs[j]*2)+
						(vdata[l][k][j+1][0]-vdata[l][k][j-1][0])/(dy*2)-
						vdata[l][k][j][0]*ltan[j]/EARTH_RADIUS;
					
				else divdata[l][k][j][0]=undef;
				
				if(udata[l][k][j][0]!=undef&&udata[l][k][j][x-2]!=undef&&vdata[l][k][j+1][x-1]!=undef&&vdata[l][k][j-1][x-1]!=undef)
					divdata[l][k][j][x-1]=
						(udata[l][k][j  ][0]-udata[l][k][j  ][x-2])/(dxs[j]*2)+
						(vdata[l][k][j+1][x-1]-vdata[l][k][j-1][x-1])/(dy*2)-
						vdata[l][k][j][x-1]*ltan[j]/EARTH_RADIUS;
					
				else divdata[l][k][j][x-1]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					/********************* div *********************/
					if(udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
						divdata[k][j][i][l]=
							(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[j]*2)+
							(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
							(vdata[k][j][i][l]/EARTH_RADIUS)*ltan[j];
						
					else divdata[k][j][i][l]=undef;
				}
				
				/*** left and right boundry of div ***/
				if(udata[k][j][1][l]!=undef&&udata[k][j][x-1][l]!=undef&&vdata[k][j+1][0][l]!=undef&&vdata[k][j-1][0][l]!=undef)
					divdata[k][j][0][l]=
						(udata[k][j  ][1][l]-udata[k][j][x-1][l])/(dxs[j]*2)+
						(vdata[k][j+1][0][l]-vdata[k][j-1][0][l])/(dy*2)-
						vdata[k][j][0][l]*ltan[j]/EARTH_RADIUS;
					
				else divdata[k][j][0][l]=undef;
				
				if(udata[k][j][0][l]!=undef&&udata[k][j][x-2][l]!=undef&&vdata[k][j+1][x-1][l]!=undef&&vdata[k][j-1][x-1][l]!=undef)
					divdata[k][j][x-1][l]=
						(udata[k][j  ][0  ][l]-udata[k][j  ][x-2][l])/(dxs[j]*2)+
						(vdata[k][j+1][x-1][l]-vdata[k][j-1][x-1][l])/(dy*2)-
						vdata[k][j][x-1][l]*ltan[j]/EARTH_RADIUS;
					
				else divdata[k][j][x-1][l]=undef;
			}
		}
		
		return div;
	}
	
	/**
     * calculate vorticity under the condition of west-east periodic boundary
     *
     * @param	u	u-wind
     * @param	v	v-wind
     *
     * @return	vorticity
     *
     * @exception	if u,v are not dimensionlly the same
     */
	public Variable c2DVorticity(Variable u,Variable v){
		if(!u.isLike(v)) throw new IllegalArgumentException("dimensions not same");
		
		t=u.getTCount();	z=u.getZCount();	y=u.getYCount();	x=u.getXCount();
		
		float undef=u.getUndef();
		
		Variable vor=new Variable("vor",u);	vor.setCommentAndUnit("vorticity (s^-1)");
		
		float[][][][]   udata=u.getData();
		float[][][][]   vdata=v.getData();
		float[][][][] vordata=vor.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					/********************* vor *********************/
					if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef&&udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef)
						vordata[l][k][j][i]=
							(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[j]*2)-
							(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)+
							udata[l][k][j][i]*ltan[j]/EARTH_RADIUS;
						
					else vordata[l][k][j][i]=undef;
				}
				
				/*** left and right boundry of vor ***/
				if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][x-1]!=undef&&udata[l][k][j+1][0]!=undef&&udata[l][k][j-1][0]!=undef)
					vordata[l][k][j][0]=
						(vdata[l][k][j][1]-vdata[l][k][j][x-1])/(dxs[j]*2)-
						(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2)+
						udata[l][k][j][0]*ltan[j]/EARTH_RADIUS;
					
				else vordata[l][k][j][0]=undef;
				
				if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-2]!=undef&&udata[l][k][j+1][x-1]!=undef&&udata[l][k][j-1][x-1]!=undef)
					vordata[l][k][j][x-1]=
						(vdata[l][k][j  ][0  ]-vdata[l][k][j  ][x-2])/(dxs[j]*2)-
						(udata[l][k][j+1][x-1]-udata[l][k][j-1][x-1])/(dy*2)+
						udata[l][k][j][x-1]*ltan[j]/EARTH_RADIUS;
					
				else vordata[l][k][j][x-1]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					/********************* vor *********************/
					if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef&&udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef)
						vordata[k][j][i][l]=
							(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[j]*2)-
							(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2)+
							udata[k][j][i][l]*ltan[j]/EARTH_RADIUS;
						
					else vordata[k][j][i][l]=undef;
				}
				
				/*** left and right boundry of vor ***/
				if(vdata[k][j][1][l]!=undef&&vdata[k][j][x-1][l]!=undef&&udata[k][j+1][0][l]!=undef&&udata[k][j-1][0][l]!=undef)
					vordata[k][j][0][l]=
						(vdata[k][j  ][1][l]-vdata[k][j  ][x-1][l])/(dxs[j]*2)-
						(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2)+
						udata[k][j][0][l]*ltan[j]/EARTH_RADIUS;
					
				else vordata[k][j][0][l]=undef;
				
				if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-2][l]!=undef&&udata[k][j+1][x-1][l]!=undef&&udata[k][j-1][x-1][l]!=undef)
					vordata[k][j][x-1][l]=
						(vdata[k][j  ][0  ][l]-vdata[k][j  ][x-2][l])/(dxs[j]*2)-
						(udata[k][j+1][x-1][l]-udata[k][j-1][x-1][l])/(dy*2)+
						udata[k][j][x-1][l]*ltan[j]/EARTH_RADIUS;
					
				else vordata[k][j][x-1][l]=undef;
			}
		}
		
		return vor;
	}
	
	
	/**
     * calculate advectiion term
     *
     * @param	a	a given variable
     * @param	u	u-wind
     * @param	v	v-wind
     *
     * @return advectiion term
     *
     * @exception	if a,u,v are not dimensionlly the same
     */
	public Variable cAdvectionTerm(Variable a,Variable u,Variable v){
		if(!a.isLike(u)||!a.isLike(v)) throw new IllegalArgumentException("dimensions not same");
		
		t=v.getTCount();	z=v.getZCount();	y=v.getYCount();	x=v.getXCount();
		
		float undef=v.getUndef();
		Variable at=new Variable("at",v);	at.setCommentAndUnit("advection term");
		
		float[][][][]  adata= a.getData();
		float[][][][]  udata= u.getData();
		float[][][][]  vdata= v.getData();
		float[][][][] atdata=at.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** u*da/dx ***/
				for(int j=0;j<y;j++){
					for(int i=1;i<x-1;i++){
						if(udata[l][k][j][i]!=undef&&adata[l][k][j][i+1]!=undef&&adata[l][k][j][i-1]!=undef)
							atdata[l][k][j][i]=udata[l][k][j][i]*(adata[l][k][j][i+1]-adata[l][k][j][i-1])/(dxs[j]*2);
							
						else atdata[l][k][j][i]=undef;
					}
					
					if(udata[l][k][j][0]!=undef&&adata[l][k][j][1]!=undef&&adata[l][k][j][x-1]!=undef)
						atdata[l][k][j][0]=udata[l][k][j][0]*(adata[l][k][j][1]-adata[l][k][j][x-1])/(dxs[j]*2);
						
					else atdata[l][k][j][0]=undef;
					
					if(udata[l][k][j][x-1]!=undef&&adata[l][k][j][0]!=undef&&adata[l][k][j][x-2]!=undef)
						atdata[l][k][j][x-1]=udata[l][k][j][x-1]*(adata[l][k][j][0]-adata[l][k][j][x-2])/(dxs[j]*2);
						
					else atdata[l][k][j][x-1]=undef;
				}
				
				/*** v*da/dy ***/
				for(int i=0;i<x;i++){
					for(int j=1;j<y-1;j++){
						if(vdata[l][k][j][i]!=undef&&adata[l][k][j+1][i]!=undef&&adata[l][k][j-1][i]!=undef){
							if(atdata[l][k][j][i]!=undef) atdata[l][k][j][i]+=vdata[l][k][j][i]*(adata[l][k][j+1][i]-adata[l][k][j-1][i])/(dy*2);
							else atdata[l][k][j][i]=vdata[l][k][j][i]*(adata[l][k][j+1][i]-adata[l][k][j-1][i])/(dy*2);
						}
					}
					
					if(vdata[l][k][0][i]!=undef&&adata[l][k][0][i]!=undef&&adata[l][k][1][i]!=undef){
						if(atdata[l][k][0][i]!=undef) atdata[l][k][0][i]+=vdata[l][k][0][i]*(adata[l][k][1][i]-adata[l][k][0][i])/dy;
						else atdata[l][k][0][i]=vdata[l][k][0][i]*(adata[l][k][1][i]-adata[l][k][0][i])/dy;
					}
					
					if(vdata[l][k][y-1][i]!=undef&&adata[l][k][y-1][i]!=undef&&adata[l][k][y-2][i]!=undef){
						if(atdata[l][k][y-1][i]!=undef) atdata[l][k][y-1][i]+=vdata[l][k][y-1][i]*(adata[l][k][y-1][i]-adata[l][k][y-2][i])/dy;
						else atdata[l][k][y-1][i]=vdata[l][k][y-1][i]*(adata[l][k][y-1][i]-adata[l][k][y-2][i])/dy;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** u*da/dx ***/
				for(int j=0;j<y;j++){
					for(int i=1;i<x-1;i++){
						if(udata[k][j][i][l]!=undef&&adata[k][j][i+1][l]!=undef&&adata[k][j][i-1][l]!=undef)
							atdata[k][j][i][l]=udata[k][j][i][l]*(adata[k][j][i+1][l]-adata[k][j][i-1][l])/(dxs[j]*2);
							
						else atdata[k][j][i][l]=undef;
					}
					
					if(udata[k][j][0][l]!=undef&&adata[k][j][1][l]!=undef&&adata[k][j][x-1][l]!=undef)
						atdata[k][j][0][l]=udata[k][j][0][l]*(adata[k][j][1][l]-adata[k][j][x-1][l])/(dxs[j]*2);
						
					else atdata[k][j][0][l]=undef;
					
					if(udata[k][j][x-1][l]!=undef&&adata[k][j][0][l]!=undef&&adata[k][j][x-2][l]!=undef)
						atdata[k][j][x-1][l]=udata[k][j][x-1][l]*(adata[k][j][0][l]-adata[k][j][x-2][l])/(dxs[j]*2);
						
					else atdata[k][j][x-1][l]=undef;
				}
				
				/*** v*da/dy ***/
				for(int i=0;i<x;i++){
					for(int j=1;j<y-1;j++){
						if(vdata[k][j][i][l]!=undef&&adata[k][j+1][i][l]!=undef&&adata[k][j-1][i][l]!=undef){
							if(atdata[l][k][j][i]!=undef) atdata[k][j][i][l]+=vdata[k][j][i][l]*(adata[k][j+1][i][l]-adata[k][j-1][i][l])/(dy*2);
							else atdata[k][j][i][l]=vdata[k][j][i][l]*(adata[k][j+1][i][l]-adata[k][j-1][i][l])/(dy*2);
						}
					}
					
					if(vdata[k][0][i][l]!=undef&&adata[k][0][i][l]!=undef&&adata[k][1][i][l]!=undef){
						if(atdata[k][0][i][l]!=undef) atdata[k][0][i][l]+=vdata[k][0][i][l]*(adata[k][1][i][l]-adata[k][0][i][l])/dy;
						else atdata[k][0][i][l]=vdata[k][0][i][l]*(adata[k][1][i][l]-adata[k][0][i][l])/dy;
					}
					
					if(vdata[k][y-1][i][l]!=undef&&adata[k][y-1][i][l]!=undef&&adata[k][y-2][i][l]!=undef){
						if(atdata[k][y-1][i][l]!=undef) atdata[k][y-1][i][l]+=vdata[k][y-1][i][l]*(adata[k][y-1][i][l]-adata[k][y-2][i][l])/dy;
						else atdata[k][y-1][i][l]=vdata[k][y-1][i][l]*(adata[k][y-1][i][l]-adata[k][y-2][i][l])/dy;
					}
				}
			}
		}
		
		return at;
	}
	
	
	/**
     * calculate 2D wave activity flux (Plumb 1985, JAS)
     *
     * @param	ua		u-wind anomaly
     * @param	va		v-wind anomaly
     * @param	ha		geopotential height anomaly
     *
     * @return	re	re[0] is zonal component and re[1] is meridional component
     *
     * @exception	if ua,va,ha are not dimensionally the same
     */
	public Variable[] cWaveActivityFlux(Variable ua,Variable va,Variable ha){
		if(!ua.isLike(va)||!ua.isLike(ha)) throw new IllegalArgumentException("dimensions not same");
		
		zstart=ua.getRange().getZRange()[0];
		
		t=ua.getTCount();	z=ua.getZCount();	y=ua.getYCount();	x=ua.getXCount();
		
		float undef=ua.getUndef();
		
		Variable[] flx=new Variable[2];
		flx[0]=new Variable("Fsx",ua);
		flx[1]=new Variable("Fsy",ua);
		flx[0].setUndef(undef);	flx[0].setCommentAndUnit("zonal component of wave activity flux");
		flx[1].setUndef(undef);	flx[1].setCommentAndUnit("meridonal component of wave activity flux");
		
		float[][][][]  udata=ua.getData();
		float[][][][]  vdata=va.getData();
		float[][][][]  hdata=ha.getData();
		float[][][][] f0data=flx[0].getData();
		float[][][][] f1data=flx[1].getData();
		
		if(ua.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++){
				for(int i=1,I=x-1;i<I;i++){
					if(vdata[l][k][j][i-1]!=undef&&vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i]!=undef
					 &&hdata[l][k][j][i-1]!=undef&&hdata[l][k][j][i+1]!=undef&&f1[j]!=0)
						f0data[l][k][j][i]=zdef[zstart-1+k]*lcos[j]*(vdata[l][k][j][i]*vdata[l][k][j][i]-
						(vdata[l][k][j][i+1]*hdata[l][k][j][i+1]-vdata[l][k][j][i-1]*hdata[l][k][j][i-1])/(2*f1[j]*dxs[j]));
					else
						f0data[l][k][j][i]=undef;
					
					if(udata[l][k][j][i-1]!=undef&&udata[l][k][j][i+1]!=undef&&udata[l][k][j][i]!=undef
					 &&hdata[l][k][j][i-1]!=undef&&hdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i]!=undef&&f1[j]!=0)
						f1data[l][k][j][i]=zdef[zstart-1+k]*lcos[j]*(-udata[l][k][j][i]*vdata[l][k][j][i]+
						(udata[l][k][j][i+1]*hdata[l][k][j][i+1]-udata[l][k][j][i-1]*hdata[l][k][j][i-1])/(2*f1[j]*dxs[j]));
					else
						f1data[l][k][j][i]=undef;
				}
				
				if(vdata[l][k][j][x-1]!=undef&&vdata[l][k][j][1]!=undef&&vdata[l][k][j][0]!=undef
				 &&hdata[l][k][j][x-1]!=undef&&hdata[l][k][j][1]!=undef&&f1[j]!=0)
					f0data[l][k][j][0]=zdef[zstart-1+k]*lcos[j]*(vdata[l][k][j][0]*vdata[l][k][j][0]-
					(vdata[l][k][j][1]*hdata[l][k][j][1]-vdata[l][k][j][x-1]*hdata[l][k][j][x-1])/(2*f1[j]*dxs[j]));
				else
					f0data[l][k][j][0]=undef;
				
				if(udata[l][k][j][x-1]!=undef&&udata[l][k][j][1]!=undef&&udata[l][k][j][0]!=undef
				 &&hdata[l][k][j][x-1]!=undef&&hdata[l][k][j][1]!=undef&&vdata[l][k][j][0]!=undef&&f1[j]!=0)
					f1data[l][k][j][0]=zdef[zstart-1+k]*lcos[j]*(-udata[l][k][j][0]*vdata[l][k][j][0]+
					(udata[l][k][j][1]*hdata[l][k][j][1]-udata[l][k][j][x-1]*hdata[l][k][j][x-1])/(2*f1[j]*dxs[j]));
				else
					f1data[l][k][j][0]=undef;
				
				if(vdata[l][k][j][x-2]!=undef&&vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-1]!=undef
				 &&hdata[l][k][j][x-2]!=undef&&hdata[l][k][j][0]!=undef&&f1[j]!=0)
					f0data[l][k][j][x-1]=zdef[zstart-1+k]*lcos[j]*(vdata[l][k][j][x-1]*vdata[l][k][j][x-1]-
					(vdata[l][k][j][0]*hdata[l][k][j][0]-vdata[l][k][j][x-2]*hdata[l][k][j][x-2])/(2*f1[j]*dxs[j]));
				else
					f0data[l][k][j][x-1]=undef;

				if(udata[l][k][j][x-2]!=undef&&udata[l][k][j][0]!=undef&&udata[l][k][j][x-1]!=undef
				 &&hdata[l][k][j][x-2]!=undef&&hdata[l][k][j][0]!=undef&&vdata[l][k][j][x-1]!=undef&&f1[j]!=0)
					f1data[l][k][j][x-1]=zdef[zstart-1+k]*lcos[j]*(-udata[l][k][j][x-1]*vdata[l][k][j][x-1]+
					(udata[l][k][j][0]*hdata[l][k][j][0]-udata[l][k][j][x-2]*hdata[l][k][j][x-2])/(2*f1[j]*dxs[j]));
				else
					f1data[l][k][j][x-1]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1,J=y-1;j<J;j++){
				for(int i=1,I=x-1;i<I;i++){
					if(vdata[k][j][i-1][l]!=undef&&vdata[k][j][i+1][l]!=undef&&vdata[k][j][i][l]!=undef
					 &&hdata[k][j][i-1][l]!=undef&&hdata[k][j][i+1][l]!=undef&&f1[j]!=0)
						f0data[k][j][i][l]=zdef[zstart-1+k]*lcos[j]*(vdata[k][j][i][l]*vdata[k][j][i][l]-
						(vdata[k][j][i+1][l]*hdata[k][j][i+1][l]-vdata[k][j][i-1][l]*hdata[k][j][i-1][l])/(2*f1[j]*dxs[j]));
					else
						f0data[k][j][i][l]=undef;
					
					if(udata[k][j][i-1][l]!=undef&&udata[k][j][i+1][l]!=undef&&udata[k][j][i][l]!=undef
					 &&hdata[k][j][i-1][l]!=undef&&hdata[k][j][i+1][l]!=undef&&vdata[k][j][i][l]!=undef&&f1[j]!=0)
						f1data[k][j][i][l]=zdef[zstart-1+k]*lcos[j]*(-udata[k][j][i][l]*vdata[k][j][i][l]+
						(udata[k][j][i+1][l]*hdata[k][j][i+1][l]-udata[k][j][i-1][l]*hdata[k][j][i-1][l])/(2*f1[j]*dxs[j]));
					else
						f1data[k][j][i][l]=undef;
				}
				
				if(vdata[k][j][x-1][l]!=undef&&vdata[k][j][1][l]!=undef&&vdata[k][j][0][l]!=undef
				 &&hdata[k][j][x-1][l]!=undef&&hdata[k][j][1][l]!=undef&&f1[j]!=0)
					f0data[k][j][0][l]=zdef[zstart-1+k]*lcos[j]*(vdata[k][j][0][l]*vdata[k][j][0][l]-
					(vdata[k][j][1][l]*hdata[k][j][1][l]-vdata[k][j][x-1][l]*hdata[k][j][x-1][l])/(2*f1[j]*dxs[j]));
				else
					f0data[k][j][0][l]=undef;
				
				if(udata[k][j][x-1][l]!=undef&&udata[k][j][1][l]!=undef&&udata[k][j][0][l]!=undef
				 &&hdata[k][j][x-1][l]!=undef&&hdata[k][j][1][l]!=undef&&vdata[k][j][0][l]!=undef&&f1[j]!=0)
					f1data[k][j][0][l]=zdef[zstart-1+k]*lcos[j]*(-udata[k][j][0][l]*vdata[k][j][0][l]+
					(udata[k][j][1][l]*hdata[k][j][1][l]-udata[k][j][x-1][l]*hdata[k][j][x-1][l])/(2*f1[j]*dxs[j]));
				else
					f1data[k][j][0][l]=undef;
				
				if(vdata[k][j][x-2][l]!=undef&&vdata[k][j][0][l]!=undef&&vdata[k][j][x-1][l]!=undef
				 &&hdata[k][j][x-2][l]!=undef&&hdata[k][j][0][l]!=undef&&f1[j]!=0)
					f0data[k][j][x-1][l]=zdef[zstart-1+k]*lcos[j]*(vdata[k][j][x-1][l]*vdata[k][j][x-1][l]-
					(vdata[k][j][0][l]*hdata[k][j][0][l]-vdata[k][j][x-2][l]*hdata[k][j][x-2][l])/(2*f1[j]*dxs[j]));
				else
					f0data[k][j][x-1][l]=undef;

				if(udata[k][j][x-2][l]!=undef&&udata[k][j][0][l]!=undef&&udata[k][j][x-1][l]!=undef
				 &&hdata[k][j][x-2][l]!=undef&&hdata[k][j][0][l]!=undef&&vdata[k][j][x-1][l]!=undef&&f1[j]!=0)
					f1data[k][j][x-1][l]=zdef[zstart-1+k]*lcos[j]*(-udata[k][j][x-1][l]*vdata[k][j][x-1][l]+
					(udata[k][j][0][l]*hdata[k][j][0][l]-udata[k][j][x-2][l]*hdata[k][j][x-2][l])/(2*f1[j]*dxs[j]));
				else
					f1data[k][j][x-1][l]=undef;
			}
		}
		
		return flx;
	}
	
	
	/**
     * calculate gradient, a vector, [0] is in x direction while [1] is in y direction
     *
     * @param	v	a given Variable
     *
     * @return	gradient (vector)
     */
	public Variable[] cGradient(Variable v){
		t=v.getTCount();	z=v.getZCount();	y=v.getYCount();	x=v.getXCount();
		
		float undef=v.getUndef();
		Variable[] grd=new Variable[2];
		grd[0]=new Variable("grdx",v);	grd[0].setCommentAndUnit("gradient of "+v.getName()+" in x-direction");
		grd[1]=new Variable("grdy",v);	grd[1].setCommentAndUnit("gradient of "+v.getName()+" in y-direction");
		
		float[][][][]    vdata=     v.getData();
		float[][][][] grdxdata=grd[0].getData();
		float[][][][] grdydata=grd[1].getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){
					for(int i=1;i<x-1;i++){
						if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef)
							grdxdata[l][k][j][i]=(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[j]*2);
							
						else grdxdata[l][k][j][i]=undef;
					}
					
					if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][x-1]!=undef)
						grdxdata[l][k][j][0]=(vdata[l][k][j][1]-vdata[l][k][j][x-1])/(dxs[j]*2);
						
					else grdxdata[l][k][j][0]=undef;
					
					if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-2]!=undef)
						grdxdata[l][k][j][x-1]=(vdata[l][k][j][0]-vdata[l][k][j][x-2])/(dxs[j]*2);
						
					else grdxdata[l][k][j][x-1]=undef;
				}
				
				for(int i=0;i<x;i++){
					for(int j=1;j<y-1;j++){
						if(vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
							grdydata[l][k][j][i]=(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2);
							
						else grdydata[l][k][j][i]=undef;
					}
					
					if(vdata[l][k][1][i]!=undef&&vdata[l][k][0][i]!=undef)
						grdydata[l][k][0][i]=(vdata[l][k][1][i]-vdata[l][k][0][i])/dy;
						
					else grdydata[l][k][0][i]=undef;
					
					if(vdata[l][k][y-1][i]!=undef&&vdata[l][k][y-2][i]!=undef)
						grdydata[l][k][y-1][i]=(vdata[l][k][y-1][i]-vdata[l][k][y-2][i])/dy;
						
					else grdydata[l][k][y-1][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++){
					for(int i=1;i<x-1;i++){
						if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef)
							grdxdata[k][j][i][l]=(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[j]*2);
							
						else grdxdata[k][j][i][l]=undef;
					}
					
					if(vdata[k][j][1][l]!=undef&&vdata[k][j][x-1][l]!=undef)
						grdxdata[k][j][0][l]=(vdata[k][j][1][l]-vdata[k][j][x-1][l])/(dxs[j]*2);
						
					else grdxdata[k][j][0][l]=undef;
					
					if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-2][l]!=undef)
						grdxdata[k][j][x-1][l]=(vdata[k][j][0][l]-vdata[k][j][x-2][l])/(dxs[j]*2);
						
					else grdxdata[k][j][x-1][l]=undef;
				}
				
				for(int i=0;i<x;i++){
					for(int j=1;j<y-1;j++){
						if(vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
							grdydata[k][j][i][l]=(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2);
							
						else grdydata[k][j][i][l]=undef;
					}
					
					if(vdata[k][1][i][l]!=undef&&vdata[k][0][i][l]!=undef)
						grdydata[k][0][i][l]=(vdata[k][1][i][l]-vdata[k][0][i][l])/dy;
						
					else grdydata[k][0][i][l]=undef;
					
					if(vdata[k][y-1][i][l]!=undef&&vdata[k][y-2][i][l]!=undef)
						grdydata[k][y-1][i][l]=(vdata[k][y-1][i][l]-vdata[k][y-2][i][l])/dy;
						
					else grdydata[k][y-1][i][l]=undef;
				}
			}
		}
		
		return grd;
	}
	
	/**
     * calculate the laplace result of a given Variable
     *
     * @param	v	a given Variable
     *
     * @return	laplace of the given Variable
     */
	public Variable cLaplace(Variable v){
		t=v.getTCount();	z=v.getZCount();	y=v.getYCount();	x=v.getXCount();
		
		float undef=v.getUndef();
		Variable la=new Variable("la",v);	la.setCommentAndUnit("laplace of "+v.getName());
		
		float[][][][]  vdata=v.getData();
		float[][][][] ladata=la.getData();
		
	    if(la.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** inner area ***/
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i]!=undef&&
					vdata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
						ladata[l][k][j][i]=
						(vdata[l][k][j][i+1]+vdata[l][k][j][i-1]-2*vdata[l][k][j][i])/(dxs[j]*dxs[j])+
						(vdata[l][k][j+1][i]+vdata[l][k][j-1][i]-2*vdata[l][k][j][i])/(dy*dy)-
						(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(2*dy)*ltan[j]/EARTH_RADIUS;
				}
				
				/*** west and east boundary ***/
				for(int j=1;j<y-1;j++){
					if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][0]!=undef&&
					vdata[l][k][j][x-1]!=undef&&vdata[l][k][j+1][0]!=undef&&vdata[l][k][j-1][0]!=undef)
						ladata[l][k][j][0]=
						(vdata[l][k][j][1]+vdata[l][k][j][x-1]-2*vdata[l][k][j][0])/(dxs[j]*dxs[j])+
						(vdata[l][k][j+1][0]+vdata[l][k][j-1][0]-2*vdata[l][k][j][0])/(dy*dy)-
						(vdata[l][k][j+1][0]-vdata[l][k][j-1][0])/(2*dy)*ltan[j]/EARTH_RADIUS;
					
					if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-1]!=undef&&
					vdata[l][k][j][x-2]!=undef&&vdata[l][k][j+1][x-1]!=undef&&vdata[l][k][j-1][x-1]!=undef)
						ladata[l][k][j][x-1]=
						(vdata[l][k][j][0]+vdata[l][k][j][x-2]-2*vdata[l][k][j][x-1])/(dxs[j]*dxs[j])+
						(vdata[l][k][j+1][x-1]+vdata[l][k][j-1][x-1]-2*vdata[l][k][j][x-1])/(dy*dy)-
						(vdata[l][k][j+1][x-1]-vdata[l][k][j-1][x-1])/(2*dy)*ltan[j]/EARTH_RADIUS;
				}
				
				/*** south and north boundary ***/
				for(int i=0;i<x;i++){
					ladata[l][k][0][i]=undef;
					ladata[l][k][y-1][i]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** inner area ***/
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i][l]!=undef&&
					vdata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
						ladata[k][j][i][l]=
						(vdata[k][j][i+1][l]+vdata[k][j][i-1][l]-2*vdata[k][j][i][l])/(dxs[j]*dxs[j])+
						(vdata[k][j+1][i][l]+vdata[k][j-1][i][l]-2*vdata[k][j][i][l])/(dy*dy)-
						(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(2*dy)*ltan[j]/EARTH_RADIUS;
				}
				
				/*** west and east boundary ***/
				for(int j=1;j<y-1;j++){
					if(vdata[k][j][1][l]!=undef&&vdata[k][j][0][l]!=undef&&
					vdata[k][j][x-1][l]!=undef&&vdata[k][j+1][0][l]!=undef&&vdata[k][j-1][0][l]!=undef)
						ladata[k][j][0][l]=
						(vdata[k][j][1][l]+vdata[k][j][x-1][l]-2*vdata[k][j][0][l])/(dxs[j]*dxs[j])+
						(vdata[k][j+1][0][l]+vdata[k][j-1][0][l]-2*vdata[k][j][0][l])/(dy*dy)-
						(vdata[k][j+1][0][l]-vdata[k][j-1][0][l])/(2*dy)*ltan[j]/EARTH_RADIUS;
					
					if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-1][l]!=undef&&
					vdata[k][j][x-2][l]!=undef&&vdata[k][j+1][x-1][l]!=undef&&vdata[k][j-1][x-1][l]!=undef)
						ladata[k][j][x-1][l]=
						(vdata[k][j][0][l]+vdata[k][j][x-2][l]-2*vdata[k][j][x-1][l])/(dxs[j]*dxs[j])+
						(vdata[k][j+1][x-1][l]+vdata[k][j-1][x-1][l]-2*vdata[k][j][x-1][l])/(dy*dy)-
						(vdata[k][j+1][x-1][l]-vdata[k][j-1][x-1][l])/(2*dy)*ltan[j]/EARTH_RADIUS;
				}
				
				/*** south and north boundary ***/
				for(int i=0;i<x;i++){
					ladata[k][0][i][l]=undef;
					ladata[k][y-1][i][l]=undef;
				}
			}
		}
		
		return la;
	}
	
	
	/** test
	public static void main(String arg[]){
		miniufo.descriptor.CtlDescriptor cd=
		new miniufo.descriptor.CtlDescriptor("D:/data/typhoon/chanchu/chanchu2.ctl");
		miniufo.diagnosis.SphericalSpacialModel ssm=new miniufo.diagnosis.SphericalSpacialModel(cd);
		DynamicMethodsInSC am=new DynamicMethodsInSC(ssm);
		
		miniufo.diagnosis.Range range=new miniufo.diagnosis.Range("",cd);
		
		miniufo.io.CtlDataReadStream cdrs=new miniufo.io.CtlDataReadStream(cd);
		
		Variable q=new Variable("q",range);
		Variable u=new Variable("u",range);
		Variable v=new Variable("v",range);
		Variable w=new Variable("w",range);
		
		cdrs.readData(u,v,w,q);
		
		Variable lp=am.cLocalPartialTerm(q);
		Variable at=am.cAdvectionTerm(q,u,v);
		Variable ct=am.cConvectiveTerm(q,w);
		
		miniufo.io.CtlDataWriteStream cdws=
		new miniufo.io.CtlDataWriteStream("D:/data/typhoon/chanchu/ac2.dat");
		
		cdws.writeData(cd,u,v,w,q,lp,at,ct);
		
		cdrs.closeFile();	cdws.closeFile();
	}*/
}