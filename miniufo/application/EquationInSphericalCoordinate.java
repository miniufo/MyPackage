/**
 * @(#)EquationInSpheralCoordinate.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable.Dimension;


/**
 * equation application in spherical coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class EquationInSphericalCoordinate extends GeoFluidApplication{
	//
	protected float[] dxs =null;	// dx scaled by cos(lat)
	protected float[] f1  =null;	// 2*omega*sin(lat)
	protected float[] f2  =null;	// 2*omega*cos(lat)
	protected float[] rcos=null;	// R * cos(lat)
	protected float[] lsin=null;	// sin(lat)
	protected float[] lcos=null;	// cos(lat)
	protected float[] ltan=null;	// tan(lat)
	protected float[] Beta=null;	// Rossby arg
	
	
	/**
     * constructor
     *
     * @param	ssm		a given spatial model in spherical coordinate
     */
	public EquationInSphericalCoordinate(SphericalSpatialModel ssm){
		super(ssm);
		
		dxs =ssm.getDXs();
		f1  =ssm.getF1();
		f2  =ssm.getF2();
		rcos=ssm.getRCos();
		lsin=ssm.getLSin();
		lcos=ssm.getLCos();
		ltan=ssm.getLTan();
		Beta=ssm.getBeta();
	}
	
	
	/**
	 * compute derivative (gradients) of a variable in a specific dimension
	 * 
	 * @param	var		a given variable
	 * @param	dim		specified dimension
	 */
	public Variable cDerivative(Variable var,Dimension dim){
		assignSubDomainParams(var);
		
		Variable der=new Variable("grd"+var.getName(),var);
		der.setCommentAndUnit("gradient of "+var.getName()+" along "+dim+" dimension");
		der.setValue(undef);
		
		float[][][][] vdata=var.getData();
		float[][][][] ddata=der.getData();
		
		switch(dim){
		case X: if(x!=1){
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					if(dxs[ystart-1+j]!=0){
						if(BCx==BoundaryCondition.Periodic){
							if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][x-1]!=undef)
							ddata[l][k][j][0  ]=(vdata[l][k][j][1]-vdata[l][k][j][x-1])/(dxs[ystart-1+j]*2);
							if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-2]!=undef)
							ddata[l][k][j][x-1]=(vdata[l][k][j][0]-vdata[l][k][j][x-2])/(dxs[ystart-1+j]*2);
							
						}else if(BCx==BoundaryCondition.Fixed){
							if(vdata[l][k][j][1  ]!=undef&&vdata[l][k][j][0  ]!=undef)
							ddata[l][k][j][0  ]=(vdata[l][k][j][1  ]-vdata[l][k][j][0  ])/dxs[ystart-1+j];
							if(vdata[l][k][j][x-1]!=undef&&vdata[l][k][j][x-2]!=undef)
							ddata[l][k][j][x-1]=(vdata[l][k][j][x-1]-vdata[l][k][j][x-2])/dxs[ystart-1+j];
						}
						
						for(int i=1;i<x-1;i++) if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef)
						ddata[l][k][j][i]=(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2);
						
					}else{
						for(int i=0;i<x;i++) ddata[l][k][j][i]=0f;
					}
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					if(dxs[ystart-1+j]!=0){
						if(BCx==BoundaryCondition.Periodic){
							if(vdata[k][j][1][l]!=undef&&vdata[k][j][x-1][l]!=undef)
							ddata[k][j][0  ][l]=(vdata[k][j][1][l]-vdata[k][j][x-1][l])/(dxs[ystart-1+j]*2);
							if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-2][l]!=undef)
							ddata[k][j][x-1][l]=(vdata[k][j][0][l]-vdata[k][j][x-2][l])/(dxs[ystart-1+j]*2);
							
						}else if(BCx==BoundaryCondition.Fixed){
							if(vdata[k][j][1  ][l]!=undef&&vdata[k][j][0  ][l]!=undef)
							ddata[k][j][0  ][l]=(vdata[k][j][1  ][l]-vdata[k][j][0  ][l])/dxs[ystart-1+j];
							if(vdata[k][j][x-1][l]!=undef&&vdata[k][j][x-2][l]!=undef)
							ddata[k][j][x-1][l]=(vdata[k][j][x-1][l]-vdata[k][j][x-2][l])/dxs[ystart-1+j];
						}
						
						for(int i=1;i<x-1;i++) if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef)
						ddata[k][j][i][l]=(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[ystart-1+j]*2);
						
					}else{
						for(int i=0;i<x;i++) ddata[k][j][i][l]=0f;
					}
				}
			}
			break;
		}
		case Y: if(y!=1){
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					if(BCy==BoundaryCondition.Fixed){
						if(vdata[l][k][1  ][i]!=undef&&vdata[l][k][0  ][i]!=undef)
						ddata[l][k][0  ][i]=(vdata[l][k][1  ][i]-vdata[l][k][0  ][i])/dy;
						if(vdata[l][k][y-1][i]!=undef&&vdata[l][k][y-2][i]!=undef)
						ddata[l][k][y-1][i]=(vdata[l][k][y-1][i]-vdata[l][k][y-2][i])/dy;
					}
					
					for(int j=1;j<y-1;j++) if(vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
					ddata[l][k][j][i]=(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					if(BCy==BoundaryCondition.Fixed){
						if(vdata[k][1  ][i][l]!=undef&&vdata[k][0  ][i][l]!=undef)
						ddata[k][0  ][i][l]=(vdata[k][1  ][i][l]-vdata[k][0  ][i][l])/dy;
						if(vdata[k][y-1][i][l]!=undef&&vdata[k][y-2][i][l]!=undef)
						ddata[k][y-1][i][l]=(vdata[k][y-1][i][l]-vdata[k][y-2][i][l])/dy;
					}
					
					for(int j=1;j<y-1;j++) if(vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
					ddata[k][j][i][l]=(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2);
				}
			}
			break;
		}
		case Z: if(z!=1){
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCz==BoundaryCondition.Fixed){
						if(vdata[l][1  ][j][i]!=undef&&vdata[l][0  ][j][i]!=undef)
						ddata[l][0  ][j][i]=(vdata[l][1  ][j][i]-vdata[l][0  ][j][i])/dz;
						if(vdata[l][z-1][j][i]!=undef&&vdata[l][z-2][j][i]!=undef)
						ddata[l][z-1][j][i]=(vdata[l][z-1][j][i]-vdata[l][z-2][j][i])/dz;
					}
					
					for(int k=1;k<z-1;k++) if(vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef)
					ddata[l][k][j][i]=(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCz==BoundaryCondition.Fixed){
						if(vdata[1  ][j][i][l]!=undef&&vdata[0  ][j][i][l]!=undef)
						ddata[0  ][j][i][l]=(vdata[1  ][j][i][l]-vdata[0  ][j][i][l])/dz;
						if(vdata[z-1][j][i][l]!=undef&&vdata[z-2][j][i][l]!=undef)
						ddata[z-1][j][i][l]=(vdata[z-1][j][i][l]-vdata[z-2][j][i][l])/dz;
					}
					
					for(int k=1;k<z-1;k++) if(vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef)
					ddata[k][j][i][l]=(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz*2);
				}
			}
			break;
		}
		case T: if(t!=1){
			if(var.isTFirst()){
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCt==BoundaryCondition.Fixed){
						if(vdata[1  ][k][j][i]!=undef&&vdata[0  ][k][j][i]!=undef)
						ddata[0  ][k][j][i]=(vdata[1  ][k][j][i]-vdata[0  ][k][j][i])/dt;
						if(vdata[t-1][k][j][i]!=undef&&vdata[t-2][k][j][i]!=undef)
						ddata[t-1][k][j][i]=(vdata[t-1][k][j][i]-vdata[t-2][k][j][i])/dt;
					}
					
					for(int l=1;l<t-1;l++) if(vdata[l+1][k][j][i]!=undef&&vdata[l-1][k][j][i]!=undef)
					ddata[l][k][j][i]=(vdata[l+1][k][j][i]-vdata[l-1][k][j][i])/(dt*2);
				}
				
			}else{
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCt==BoundaryCondition.Fixed){
						if(vdata[k][j][i][1  ]!=undef&&vdata[k][j][i][0  ]!=undef)
						ddata[k][j][i][0  ]=(vdata[k][j][i][1  ]-vdata[k][j][i][0  ])/dt;
						if(vdata[k][j][i][t-1]!=undef&&vdata[k][j][i][t-2]!=undef)
						ddata[k][j][i][t-1]=(vdata[k][j][i][t-1]-vdata[k][j][i][t-2])/dt;
					}
					
					for(int l=1;l<t-1;l++) if(vdata[k][j][i][l+1]!=undef&&vdata[k][j][i][l-1]!=undef)
					ddata[k][j][i][l]=(vdata[k][j][i][l+1]-vdata[k][j][i][l-1])/(dt*2);
				}
			}
			break;
		}
		default: throw new IllegalArgumentException("unsupported dimension: "+dim);
		}
		
		der.changeNaNToUndef();
		
		return der;
	}
	
	/**
     * weighting a variable with cos(latitude)
     *
     * @param	v	a given variable
     */
	public Variable weightLCosEq(Variable v){
		assignSubDomainParams(v);
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(vdata[l][k][j][i]!=undef) vdata[l][k][j][i]*=lcos[ystart-1+j];
			
		}else{
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(vdata[k][j][i][l]!=undef) vdata[k][j][i][l]*=lcos[ystart-1+j];
		}
		
		return v;
	}
	
	public Variable weightLCos(Variable v){
		assignSubDomainParams(v);
		
		Variable nv=v.copy();
		
		float[][][][] ndata=nv.getData();
		
		if(v.isTFirst()){
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[l][k][j][i]!=undef) ndata[l][k][j][i]*=lcos[ystart-1+j];
			
		}else{
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[k][j][i][l]!=undef) ndata[k][j][i][l]*=lcos[ystart-1+j];
		}
		
		return nv;
	}
	
	/**
     * de-weighting a variable with cos(latitude)
     *
     * @param	v	a given variable
     */
	public Variable deWeightLCosEq(Variable v){
		assignSubDomainParams(v);
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vdata[l][k][j][i]!=undef&&lcos[ystart-1+j]!=0) vdata[l][k][j][i]/=lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vdata[k][j][i][l]!=undef&&lcos[ystart-1+j]!=0) vdata[k][j][i][l]/=lcos[ystart-1+j];
		}
		
		return v;
	}
	
	public Variable deWeightLCos(Variable v){
		assignSubDomainParams(v);
		
		Variable nv=v.copy();
		
		float[][][][] ndata=nv.getData();
		
		if(v.isTFirst()){
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[l][k][j][i]!=undef) ndata[l][k][j][i]/=lcos[ystart-1+j];
			
		}else{
			for(int j=0;j<y;j++)
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++)
			if(ndata[k][j][i][l]!=undef) ndata[k][j][i][l]/=lcos[ystart-1+j];
		}
		
		return nv;
	}
	
}
