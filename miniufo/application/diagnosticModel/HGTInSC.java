/**
 * @(#)HGTInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.diagnosticModel;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInSphericalCoordinate;
import static java.lang.Math.cos;
import static miniufo.diagnosis.SpatialModel.REarth;


/**
 * geopotential height diagnosis equation in spheral coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class HGTInSC extends EquationInSphericalCoordinate implements EllipticEquationInterface{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public HGTInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isLinearModel()) System.out.println("\nNot a equal-space model");
	}
	
	
	/**
     * calculate forces one by one
     */
	public Variable cTerm1(Variable div){
		assignSubDomainParams(div);
		
		Variable F=new Variable("F1",div);	F.setCommentAndUnit("tendency of divergence");
		
		float[][][][] Fdata=  F.getData();
		float[][][][] ddata=div.getData();
		
		if(F.isTFirst()){
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			Fdata[l][k][j][i]=-(ddata[l+1][k][j][i]-ddata[l-1][k][j][i])/(2*dt)*lcos[ystart-1+j];
			
		}else{
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			Fdata[k][j][i][l]=-(ddata[k][j][i][l+1]-ddata[k][j][i][l-1])/(2*dt)*lcos[ystart-1+j];
		}
		
		return F;
	}
	
	public Variable cTerm2(Variable u){
		assignSubDomainParams(u);
		
		Variable F=new Variable("F2",u);	F.setCommentAndUnit("zonal-jet term");
		
		float[][][][] Fdata=F.getData();
		float[][][][] udata=u.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[l][k][j][0]=-(
						(udata[l][k][j][1]+udata[l][k][j][0  ])/2f*
						(udata[l][k][j][1]-udata[l][k][j][0  ])/dxs[ystart-1+j]-
						(udata[l][k][j][0]+udata[l][k][j][x-1])/2f*
						(udata[l][k][j][0]-udata[l][k][j][x-1])/dxs[ystart-1+j]
					)/dx;
					
					/*** East boundary ***/
					Fdata[l][k][j][x-1]=-(
						(udata[l][k][j][0  ]+udata[l][k][j][x-1])/2f*
						(udata[l][k][j][0  ]-udata[l][k][j][x-1])/dxs[ystart-1+j]-
						(udata[l][k][j][x-1]+udata[l][k][j][x-2])/2f*
						(udata[l][k][j][x-1]-udata[l][k][j][x-2])/dxs[ystart-1+j]
					)/dx;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=-(
					(udata[l][k][j][i+1]+udata[l][k][j][i])/2f*
					(udata[l][k][j][i+1]-udata[l][k][j][i])/dxs[ystart-1+j]-
					(udata[l][k][j][i]+udata[l][k][j][i-1])/2f*
					(udata[l][k][j][i]-udata[l][k][j][i-1])/dxs[ystart-1+j]
				)/dx;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[k][j][0][l]=-(
						(udata[k][j][1][l]+udata[k][j][0  ][l])/2f*
						(udata[k][j][1][l]-udata[k][j][0  ][l])/dxs[ystart-1+j]-
						(udata[k][j][0][l]+udata[k][j][x-1][l])/2f*
						(udata[k][j][0][l]-udata[k][j][x-1][l])/dxs[ystart-1+j]
					)/dx;
					
					/*** East boundary ***/
					Fdata[k][j][x-1][l]=-(
						(udata[k][j][0  ][l]+udata[k][j][x-1][l])/2f*
						(udata[k][j][0  ][l]-udata[k][j][x-1][l])/dxs[ystart-1+j]-
						(udata[k][j][x-1][l]+udata[k][j][x-2][l])/2f*
						(udata[k][j][x-1][l]-udata[k][j][x-2][l])/dxs[ystart-1+j]
					)/dx;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=-(
					(udata[k][j][i+1][l]+udata[k][j][i][l])/2f*
					(udata[k][j][i+1][l]-udata[k][j][i][l])/dxs[ystart-1+j]-
					(udata[k][j][i][l]+udata[k][j][i-1][l])/2f*
					(udata[k][j][i][l]-udata[k][j][i-1][l])/dxs[ystart-1+j]
				)/dx;
			}
		}
		
		return F;
	}
	
	public Variable cTerm3(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable F=new Variable("F3",u);	F.setCommentAndUnit("saddle pattern flow");
		
		float[][][][] Fdata=F.getData();
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[l][k][j][0]=-(
						udata[l][k][j+1][0]*(vdata[l][k][j+1][1]-vdata[l][k][j+1][x-1])-
						udata[l][k][j-1][0]*(vdata[l][k][j-1][1]-vdata[l][k][j-1][x-1])
					)/(dy*dx*4);
					
					/*** East boundary ***/
					Fdata[l][k][j][x-1]=-(
						udata[l][k][j+1][x-1]*(vdata[l][k][j+1][0]-vdata[l][k][j+1][x-2])-
						udata[l][k][j-1][x-1]*(vdata[l][k][j-1][0]-vdata[l][k][j-1][x-2])
					)/(dy*dx*4);
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=-(
					udata[l][k][j+1][i]*(vdata[l][k][j+1][i+1]-vdata[l][k][j+1][i-1])-
					udata[l][k][j-1][i]*(vdata[l][k][j-1][i+1]-vdata[l][k][j-1][i-1])
				)/(dy*dx*4);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[k][j][0][l]=-(
						udata[k][j+1][0][l]*(vdata[k][j+1][1][l]-vdata[k][j+1][x-1][l])-
						udata[k][j-1][0][l]*(vdata[k][j-1][1][l]-vdata[k][j-1][x-1][l])
					)/(dy*dx*4);
					
					/*** East boundary ***/
					Fdata[k][j][x-1][l]=-(
						udata[k][j+1][x-1][l]*(vdata[k][j+1][0][l]-vdata[k][j+1][x-2][l])-
						udata[k][j-1][x-1][l]*(vdata[k][j-1][0][l]-vdata[k][j-1][x-2][l])
					)/(dy*dx*4);
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=-(
					udata[k][j+1][i][l]*(vdata[k][j+1][i+1][l]-vdata[k][j+1][i-1][l])-
					udata[k][j-1][i][l]*(vdata[k][j-1][i+1][l]-vdata[k][j-1][i-1][l])
				)/(dy*dx*4);
			}
		}
		
		return F;
	}
	
	public Variable cTerm4(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable F=new Variable("F4",u);	F.setCommentAndUnit("saddle pattern flow");
		
		float[][][][] Fdata=F.getData();
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[l][k][j][0]=-(
						vdata[l][k][j][1  ]*(udata[l][k][j+1][1  ]-udata[l][k][j-1][1  ])-
						vdata[l][k][j][x-1]*(udata[l][k][j+1][x-1]-udata[l][k][j-1][x-1])
					)/(dx*dy*4);
					
					/*** East boundary ***/
					Fdata[l][k][j][x-1]=-(
						vdata[l][k][j][0  ]*(udata[l][k][j+1][0  ]-udata[l][k][j-1][0  ])-
						vdata[l][k][j][x-2]*(udata[l][k][j+1][x-2]-udata[l][k][j-1][x-2])
					)/(dx*dy*4);
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=-(
					vdata[l][k][j][i+1]*(udata[l][k][j+1][i+1]-udata[l][k][j-1][i+1])-
					vdata[l][k][j][i-1]*(udata[l][k][j+1][i-1]-udata[l][k][j-1][i-1])
				)/(dx*dy*4);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[k][j][0][l]=-(
						vdata[k][j][1  ][l]*(udata[k][j+1][1  ][l]-udata[k][j-1][1  ][l])-
						vdata[k][j][x-1][l]*(udata[k][j+1][x-1][l]-udata[k][j-1][x-1][l])
					)/(dx*dy*4);
					
					/*** East boundary ***/
					Fdata[k][j][x-1][l]=-(
						vdata[k][j][0  ][l]*(udata[k][j+1][0  ][l]-udata[k][j-1][0  ][l])-
						vdata[k][j][x-2][l]*(udata[k][j+1][x-2][l]-udata[k][j-1][x-2][l])
					)/(dx*dy*4);
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=-(
					vdata[k][j][i+1][l]*(udata[k][j+1][i+1][l]-udata[k][j-1][i+1][l])-
					vdata[k][j][i-1][l]*(udata[k][j+1][i-1][l]-udata[k][j-1][i-1][l])
				)/(dx*dy*4);
			}
		}
		
		return F;
	}
	
	public Variable cTerm5(Variable v){
		assignSubDomainParams(v);
		
		Variable F=new Variable("F5",v);	F.setCommentAndUnit("meridional-jet term");
		
		float[][][][] Fdata=F.getData();
		float[][][][] vdata=v.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			Fdata[l][k][j][i]=-(
				(vdata[l][k][j+1][i]+vdata[l][k][j][i])/2f*
				(float)cos((ydef[ystart+j  ]+ydef[ystart-1+j])/2f)*
				(vdata[l][k][j+1][i]-vdata[l][k][j][i])/dy-
				(vdata[l][k][j][i]+vdata[l][k][j-1][i])/2f*
				(float)cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)*
				(vdata[l][k][j][i]-vdata[l][k][j-1][i])/dy
			)/dy;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			Fdata[k][j][i][l]=-(
				(vdata[k][j+1][i][l]+vdata[k][j][i][l])/2f*
				(float)cos((ydef[ystart+j]+ydef[ystart-1+j])/2f)*
				(vdata[k][j+1][i][l]-vdata[k][j][i][l])/dy-
				(vdata[k][j][i][l]+vdata[k][j-1][i][l])/2f*
				(float)cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)*
				(vdata[k][j][i][l]-vdata[k][j-1][i][l])/dy
			)/dy;
		}
		
		return F;
	}
	
	public Variable cTerm6(Variable w,Variable div){
		checkDimensions(w,div);
		assignSubDomainParams(w);
		
		Variable F=new Variable("F6",w);	F.setCommentAndUnit("convection of divergence");
		
		float[][][][] Fdata=  F.getData();
		float[][][][] ddata=div.getData();
		float[][][][] odata=  w.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			Fdata[l][k][j][i]=-odata[l][k][j][i]*lcos[ystart-1+j]*
			(ddata[l][k+1][j][i]-ddata[l][k-1][j][i])/(dz+dz);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			Fdata[k][j][i][l]=-odata[k][j][i][l]*lcos[ystart-1+j]*
			(ddata[k+1][j][i][l]-ddata[k-1][j][i][l])/(dz+dz);
		}
		
		return F;
	}
	
	public Variable cTerm7(Variable u,Variable w){
		checkDimensions(w,u);
		assignSubDomainParams(u);
			
		Variable F=new Variable("F7",u);	F.setCommentAndUnit("Walker circulation term");
		
		float[][][][] Fdata=F.getData();
		float[][][][] udata=u.getData();
		float[][][][] odata=w.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[l][k][j][0]=-
						(odata[l][k  ][j][1]-odata[l][k][j][x-1])/(dx+dx)*
						(udata[l][k+1][j][0]-udata[l][k-1][j][0])/(dz+dz);
					
					/*** East boundary ***/
					Fdata[l][k][j][x-1]=-
						(odata[l][k  ][j][0  ]-odata[l][k  ][j][x-2])/(dx+dx)*
						(udata[l][k+1][j][x-1]-udata[l][k-1][j][x-1])/(dz+dz);
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=-
					(odata[l][k][j][i+1]-odata[l][k][j][i-1])/(dx+dx)*
					(udata[l][k+1][j][i]-udata[l][k-1][j][i])/(dz+dz);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[k][j][0][l]=-
						(odata[k  ][j][1][l]-odata[k][j][x-1][l])/(dx+dx)*
						(udata[k+1][j][0][l]-udata[k-1][j][0][l])/(dz+dz);
					
					/*** East boundary ***/
					Fdata[k][j][x-1][l]=-
						(odata[k  ][j][0  ][l]-odata[k  ][j][x-2][l])/(dx+dx)*
						(udata[k+1][j][x-1][l]-udata[k-1][j][x-1][l])/(dz+dz);
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=-
					(odata[k][j][i+1][l]-odata[k][j][i-1][l])/(dx+dx)*
					(udata[k+1][j][i][l]-udata[k-1][j][i][l])/(dz+dz);
			}
		}
		
		return F;
	}
	
	public Variable cTerm8(Variable v,Variable w){
		checkDimensions(w,v);
		assignSubDomainParams(v);
		
		Variable F=new Variable("F8",v);	F.setCommentAndUnit("Hadley circulation term");
		
		float[][][][] Fdata=F.getData();
		float[][][][] vdata=v.getData();
		float[][][][] odata=w.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			Fdata[l][k][j][i]=-lcos[ystart-1+j]*
				(odata[l][k][j+1][i]-odata[l][k][j-1][i])/(dy+dy)*
				(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			Fdata[k][j][i][l]=-lcos[ystart-1+j]*
				(odata[k][j+1][i][l]-odata[k][j-1][i][l])/(dy+dy)*
				(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz);
		}
		
		return F;
	}
	
	public Variable cTerm9(Variable vor){
		assignSubDomainParams(vor);
		
		Variable F=new Variable("F9",vor);	F.setCommentAndUnit("geostropic term");
		
		float[][][][] Fdata=  F.getData();
		float[][][][] vdata=vor.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			Fdata[l][k][j][i]=f1[ystart-1+j]*vdata[l][k][j][i]*lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			Fdata[k][j][i][l]=f1[ystart-1+j]*vdata[k][j][i][l]*lcos[ystart-1+j];
		}
		
		return F;
	}
	
	public Variable cTerm10(Variable u){
		assignSubDomainParams(u);
		
		Variable F=new Variable("F10",u);	F.setCommentAndUnit("planetary vorticity term");
		
		float[][][][] Fdata=F.getData();
		float[][][][] udata=u.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Fdata[l][k][j][i]=
			-udata[l][k][j][i]*Beta[ystart-1+j]*lcos[ystart-1+j];
			//-udata[l][k][j][i]*(Beta[ystart-1+j]-f1[ystart-1+j]*ltan[ystart-1+j]/EARTH_RADIUS)*lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){ Fdata[k][j][i][l]=
			-udata[k][j][i][l]*Beta[ystart-1+j]*lcos[ystart-1+j];
			}
		}
		
		return F;
	}
	
	public Variable cTerm11(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable F=new Variable("F11",u);	F.setCommentAndUnit("spherical term");
		
		float[][][][] Fdata=F.getData();
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[l][k][j][0]=(
						udata[l][k][j][1  ]*vdata[l][k][j][1  ]-
						udata[l][k][j][x-1]*vdata[l][k][j][x-1]
					)*ltan[ystart-1+j]/(dx*2)/REarth;
					
					/*** East boundary ***/
					Fdata[l][k][j][x-1]=(
						udata[l][k][j][0  ]*vdata[l][k][j][0  ]-
						udata[l][k][j][x-2]*vdata[l][k][j][x-2]
					)*ltan[ystart-1+j]/(dx*2)/REarth;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=(
					udata[l][k][j][i+1]*vdata[l][k][j][i+1]-
					udata[l][k][j][i-1]*vdata[l][k][j][i-1]
				)*ltan[ystart-1+j]/(dx*2)/REarth;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** West boundary ***/
					Fdata[k][j][0][l]=(
						udata[k][j][1  ][l]*vdata[k][j][1  ][l]-
						udata[k][j][x-1][l]*vdata[k][j][x-1][l]
					)*ltan[ystart-1+j]/(dx*2)/REarth;
					
					/*** East boundary ***/
					Fdata[k][j][x-1][l]=(
						udata[k][j][0  ][l]*vdata[k][j][0  ][l]-
						udata[k][j][x-2][l]*vdata[k][j][x-2][l]
					)*ltan[ystart-1+j]/(dx*2)/REarth;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=(
					udata[k][j][i+1][l]*vdata[k][j][i+1][l]-
					udata[k][j][i-1][l]*vdata[k][j][i-1][l]
				)*ltan[ystart-1+j]/(dx*2)/REarth;
			}
		}
		
		return F;
	}
	
	public Variable cTerm12(Variable u){
		assignSubDomainParams(u);
		
		Variable F=new Variable("F12",u);	F.setCommentAndUnit("spherical term");
		
		float[][][][] Fdata=F.getData();
		float[][][][] udata=u.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			Fdata[l][k][j][i]=-(
				udata[l][k][j+1][i]*udata[l][k][j+1][i]*lsin[ystart  +j]-
				udata[l][k][j-1][i]*udata[l][k][j-1][i]*lsin[ystart-2+j]
			)/(dy*2)/REarth;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
			Fdata[k][j][i][l]=-(
				udata[k][j+1][i][l]*udata[k][j+1][i][l]*lsin[ystart  +j]-
				udata[k][j-1][i][l]*udata[k][j-1][i][l]*lsin[ystart-2+j]
			)/(dy*2)/REarth;
		}
		
		return F;
	}
	
	
	/**
     * calculate all forces
     */
	public Variable cAllTerms(Variable u,Variable v,Variable w){
		DynamicMethodsInSC bm=new DynamicMethodsInSC((SphericalSpatialModel)sm);
		
		bm.setBCofX(BCx);
		bm.setBCofY(BCy);
		
		Variable div=bm.c2DDivergence(u,v);
		Variable vor=bm.c2DVorticity(u,v);
		
		return cAllTerms(u,v,w,div,vor);
	}
	
	public Variable cAllTerms(Variable u,Variable v,Variable w,Variable div,Variable vor){
		checkDimensions(u,v,w,div,vor);
		
		Variable all=cTerm1(div)
			.plusEq(cTerm2(u)).plusEq(cTerm3(u,v)).plusEq(cTerm4(u,v)).plusEq(cTerm5(v))
			.plusEq(cTerm6(w,div)).plusEq(cTerm7(u,w)).plusEq(cTerm8(v,w))
			.plusEq(cTerm9(vor)).plusEq(cTerm10(u))
			.plusEq(cTerm11(u,v)).plusEq(cTerm12(u));
		
		return all;
	}
	
	
	/**
     * implement the methods in EllipticEquationInterface
     */
	public Variable cAPrime(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Ap",v); A.setCommentAndUnit("coefficient A' of elliptic equation");
		
		float[][][][] Adata=A.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[l][k][j][i]=1f/lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=1f/lcos[ystart-1+j];
		}
		
		return A;
	}
	
	public Variable cBPrime(Variable v){
		Variable B=new Variable("Bp",v);	B.setCommentAndUnit("coefficient B' of elliptic equation");
		return B;
	}
	
	public Variable cCPrime(Variable v){
		assignSubDomainParams(v);
		
		Variable C=new Variable("Cp",v); C.setCommentAndUnit("coefficient C' of elliptic equation");
		
		float[][][][] Cdata=C.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			for(int i=0;i<x;i++) Cdata[l][k][j][i]=(float)cos((ydef[ystart-2+j]+ydef[ystart-1+j])/2f);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			for(int i=0;i<x;i++) Cdata[k][j][i][l]=(float)cos((ydef[ystart-2+j]+ydef[ystart-1+j])/2f);
		}
		
		return C;
	}
	
	public Variable cA(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Aa",v); A.setCommentAndUnit("coefficient A of elliptical equation");
		
		float[][][][] Adata=A.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[l][k][j][i]=lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=lcos[ystart-1+j];
		}
		
		return A;
	}
	
	public Variable cB(Variable v){
		Variable B=new Variable("Bb",v);	B.setCommentAndUnit("coefficient B of elliptic equation");
		return B;
	}
	
	public Variable cC(Variable v){
		Variable C=cAPrime(v);
		C.setName("Cc"); C.setCommentAndUnit("coefficient C of elliptic equation");
		
		return C;
	}
	
	
	/** test
	public static void main(String[] args){
		ConcurrentUtil.initDefaultExecutor(1);
		try{
			miniufo.descriptor.CtlDescriptor cd=
			new miniufo.descriptor.CtlDescriptor(new java.io.File("D:/data/DiagnosisVortex/Haima/Haima.ctl"));
			miniufo.diagnosis.SphericalSpatialModel ssm=new miniufo.diagnosis.SphericalSpatialModel(cd);
			
			HGTInSC mHGT=new HGTInSC(ssm);mHGT.setBCofX(BoundaryCondition.Periodic);
			DynamicMethodsInSC gdm=new DynamicMethodsInSC(ssm); gdm.setBCofX(BoundaryCondition.Periodic);
			EllipticEqSORSolver ees=new EllipticEqSORSolver(ssm);
			GlobalLaplaceEquationInSC le=new GlobalLaplaceEquationInSC(ssm);
			SphericalHarmonicExpansion she=new SphericalHarmonicExpansion(ssm);
			
			miniufo.diagnosis.Range range=
			new miniufo.diagnosis.Range("t(2,2);z(36,36)",cd);
			//new miniufo.diagnosis.Range("lon(80,180);lat(10,60);t(1,1);z(36,36)",cd);
			
			Variable u=new Variable("u",false,range);
			Variable v=new Variable("v",false,range);
			//Variable w=new Variable("w",false,range);
			//Variable t=new Variable("t",false,range);
			Variable h=new Variable("h",false,range);
			
			miniufo.io.CtlDataReadStream  cdrs=new miniufo.io.CtlDataReadStream(cd);
			
			cdrs.readData(u,v,h);	h.multiplyEq(9.8f);
			
			Variable div=gdm.c2DDivergence(u,v);
			Variable vor=gdm.c2DVorticity(u,v);
			
			//Variable h1 =(Variable)(h.copy());	h1.setName("h1");	h1.setValue(0);
			Variable h2 =(Variable)(h.copy());	h2.setName("h2");	h2.setValue(0);
			Variable h3 =(Variable)(h.copy());	h3.setName("h3");	h3.setValue(0);
			Variable h4 =(Variable)(h.copy());	h4.setName("h4");	h4.setValue(0);
			Variable h5 =(Variable)(h.copy());	h5.setName("h5");	h5.setValue(0);
			//Variable h6 =(Variable)(h.copy());	h6.setName("h6");	h6.setValue(0);
			//Variable h7 =(Variable)(h.copy());	h7.setName("h7");	h7.setValue(0);
			//Variable h8 =(Variable)(h.copy());	h8.setName("h8");	h8.setValue(0);
			Variable h9 =(Variable)(h.copy());	h9.setName("h9");	h9.setValue(0);
			Variable h10=(Variable)(h.copy());	h10.setName("h10");	h10.setValue(0);
			Variable h11=(Variable)(h.copy());	h11.setName("h11");	h11.setValue(0);
			Variable h12=(Variable)(h.copy());	h12.setName("h12");	h12.setValue(0);
			Variable hb =(Variable)(h.copy());	hb.setName("hb");	hb.setInner(0);
			Variable hh =(Variable)(h.copy());	hh.setName("hh");	hh.setInner(0);
			
			//Variable hh1 =(Variable)(h.copy());	hh1.setName("hh1");		hh1.setValue(0);
			Variable hh2 =(Variable)(h.copy());	hh2.setName("hh2");		hh2.setValue(0);
			Variable hh3 =(Variable)(h.copy());	hh3.setName("hh3");		hh3.setValue(0);
			Variable hh4 =(Variable)(h.copy());	hh4.setName("hh4");		hh4.setValue(0);
			Variable hh5 =(Variable)(h.copy());	hh5.setName("hh5");		hh5.setValue(0);
			//Variable hh6 =(Variable)(h.copy());	hh6.setName("hh6");		hh6.setValue(0);
			//Variable hh7 =(Variable)(h.copy());	hh7.setName("hh7");		hh7.setValue(0);
			//Variable hh8 =(Variable)(h.copy());	hh8.setName("hh8");		hh8.setValue(0);
			Variable hh9 =(Variable)(h.copy());	hh9.setName("hh9");		hh9.setValue(0);
			Variable hh10=(Variable)(h.copy());	hh10.setName("hh10");	hh10.setValue(0);
			Variable hh11=(Variable)(h.copy());	hh11.setName("hh11");	hh11.setValue(0);
			Variable hh12=(Variable)(h.copy());	hh12.setName("hh12");	hh12.setValue(0);
			Variable hhb =(Variable)(h.copy());	hhb.setName("hhb");		hhb.setInner(0);
			Variable hhh =(Variable)(h.copy());	hhh.setName("hhh");		hhh.setInner(0);
			
			//Variable f1 =mHGT.cTerm1(div);
			Variable f2 =mHGT.cTerm2(u);
			Variable f3 =mHGT.cTerm3(u,v);
			Variable f4 =mHGT.cTerm4(u,v);
			Variable f5 =mHGT.cTerm5(v);
			//Variable f6 =mHGT.cTerm6(w,div);
			//Variable f7 =mHGT.cTerm7(u,w);
			//Variable f8 =mHGT.cTerm8(v,w);
			Variable f9 =mHGT.cTerm9(vor);
			Variable f10=mHGT.cTerm10(u);
			Variable f11=mHGT.cTerm11(u,v);
			Variable f12=mHGT.cTerm12(u);
			Variable ff =f2.plus(f3).plusEq(f4).plusEq(f5).plusEq(f9).plusEq(f10).plusEq(f11).plusEq(f12);
			//Variable ff =f1.plus(f2).plusEq(f3).plusEq(f4).plusEq(f5).plusEq(f6).plusEq(f7).plusEq(f8).plusEq(f9).plusEq(f10).plusEq(f11).plusEq(f12);
			
			//le.solve(false,5000,hh1,f1);	hh1.divideEq(9.8f);
			le.solve(false,5000,hh2,f2);	hh2.divideEq(9.8f);
			le.solve(false,5000,hh3,f3); 	hh3.divideEq(9.8f);
			le.solve(false,5000,hh4,f4); 	hh4.divideEq(9.8f);
			le.solve(false,5000,hh5,f5); 	hh5.divideEq(9.8f);
			//le.solve(false,5000,hh6,f6); 	hh6.divideEq(9.8f);
			//le.solve(false,5000,hh7,f7); 	hh7.divideEq(9.8f);
			//le.solve(false,5000,hh8,f8);	hh8.divideEq(9.8f);
			le.solve(false,5000,hh9,f9); 	hh9.divideEq(9.8f);
			le.solve(false,5000,hh10,f10); 	hh10.divideEq(9.8f);
			le.solve(false,5000,hh11,f11); 	hh11.divideEq(9.8f);
			le.solve(false,5000,hh12,f12);	hh12.divideEq(9.8f);
			le.solve(false,5000,hhh,ff);	hhh.divideEq(9.8f);
			le.solve(false,5000,hhb,null);	hhb.divideEq(9.8f);
			
			ees.setDimCombination(DimCombination.XY);
			ees.setBCofX(BoundaryCondition.Periodic);
			ees.setABC(mHGT.cAPrime(vor),mHGT.cBPrime(vor),mHGT.cCPrime(vor));
			
			//ees.solve(h1,f1);		h1.divideEq(9.8f);
			ees.solve(h2,f2);		h2.divideEq(9.8f);
			ees.solve(h3,f3);		h3.divideEq(9.8f);
			ees.solve(h4,f4);		h4.divideEq(9.8f);
			ees.solve(h5,f5);		h5.divideEq(9.8f);
			//ees.solve(h6,f6);		h6.divideEq(9.8f);
			//ees.solve(h7,f7);		h7.divideEq(9.8f);
			//ees.solve(h8,f8);		h8.divideEq(9.8f);
			ees.solve(h9,f9);		h9.divideEq(9.8f);
			ees.solve(h10,f10);		h10.divideEq(9.8f);
			ees.solve(h11,f11);		h11.divideEq(9.8f);
			ees.solve(h12,f12);		h12.divideEq(9.8f);
			ees.solve(hh,ff);		hh.divideEq(9.8f);
			ees.solve(hb,null);		hb.divideEq(9.8f);
			
			mHGT.deWeightLCos(f2);
			mHGT.deWeightLCos(f3);
			mHGT.deWeightLCos(f4);
			mHGT.deWeightLCos(f5);
			mHGT.deWeightLCos(f9);
			mHGT.deWeightLCos(f10);
			mHGT.deWeightLCos(f11);
			mHGT.deWeightLCos(f12);
			
			she.setM(180);
			Variable hs2 =she.solvePoissonEquation(f2 ); hs2.setName("hs2"); hs2.divideEq(9.8f);
			Variable hs3 =she.solvePoissonEquation(f3 ); hs3.setName("hs3"); hs3.divideEq(9.8f);
			Variable hs4 =she.solvePoissonEquation(f4 ); hs4.setName("hs4"); hs4.divideEq(9.8f);
			Variable hs5 =she.solvePoissonEquation(f5 ); hs5.setName("hs5"); hs5.divideEq(9.8f);
			Variable hs9 =she.solvePoissonEquation(f9 ); hs9.setName("hs9"); hs9.divideEq(9.8f);
			Variable hs10=she.solvePoissonEquation(f10); hs10.setName("hs10"); hs10.divideEq(9.8f);
			Variable hs11=she.solvePoissonEquation(f11); hs11.setName("hs11"); hs11.divideEq(9.8f);
			Variable hs12=she.solvePoissonEquation(f12); hs12.setName("hs12"); hs12.divideEq(9.8f);
			Variable hs  =she.solvePoissonEquation(ff ); hs.setName("hs"); hs.divideEq(9.8f);
			
			h.divideEq(9.8f);
			
			miniufo.io.CtlDataWriteStream cdws=new miniufo.io.CtlDataWriteStream("d:/hgt2.dat");
			cdws.writeData(cd,u,v,h,div,vor,
				hb,hh,h2,h3,h4,h5,h9,h10,h11,h12,
				hhb,hhh,hh2,hh3,hh4,hh5,hh9,hh10,hh11,hh12,
				hs,hs2,hs3,hs4,hs5,hs9,hs10,hs11,hs12,
				f2,f3,f4,f5,f9,f10,f11,f12);
			cdrs.closeFile();	cdws.closeFile();
			
	    }catch(Exception ex){ ex.printStackTrace();}
		ConcurrentUtil.shutdown();
	}*/
}
