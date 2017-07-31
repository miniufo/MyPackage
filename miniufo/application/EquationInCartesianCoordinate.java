/**
 * @(#)EquationInCartesianCoordinate.java	1.0 2017.06.21
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.diagnosis.Variable.Dimension;


/**
 * Equation application in Cartesian coordinates.
 *
 * @version 1.0, 2017.06.21
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class EquationInCartesianCoordinate extends GeoFluidApplication{
	//
	protected float beta=0;			// a constant = df/dy
	
	protected float[] fCor=null;	// = f0 + beta*y
	
	
	/**
     * constructor
     *
     * @param	csm		a given spatial model in Cartesian coordinates
     */
	public EquationInCartesianCoordinate(CartesianSpatialModel csm){
		super(csm);
		
		beta=csm.getBeta();
		fCor=csm.getFCor();
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
		case X:{
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					if(BCx==BoundaryCondition.Periodic){
						if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][x-1]!=undef)
						ddata[l][k][j][0  ]=(vdata[l][k][j][1]-vdata[l][k][j][x-1])/(dx*2f);
						if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-2]!=undef)
						ddata[l][k][j][x-1]=(vdata[l][k][j][0]-vdata[l][k][j][x-2])/(dx*2f);
						
					}else if(BCx==BoundaryCondition.Fixed){
						if(vdata[l][k][j][1  ]!=undef&&vdata[l][k][j][0  ]!=undef)
						ddata[l][k][j][0  ]=(vdata[l][k][j][1  ]-vdata[l][k][j][0  ])/dx;
						if(vdata[l][k][j][x-1]!=undef&&vdata[l][k][j][x-2]!=undef)
						ddata[l][k][j][x-1]=(vdata[l][k][j][x-1]-vdata[l][k][j][x-2])/dx;
					}
					
					for(int i=1;i<x-1;i++)
					if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef)
					ddata[l][k][j][i]=(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dx*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					if(BCx==BoundaryCondition.Periodic){
						if(vdata[k][j][1][l]!=undef&&vdata[k][j][x-1][l]!=undef)
						ddata[k][j][0  ][l]=(vdata[k][j][1][l]-vdata[k][j][x-1][l])/(dx*2f);
						if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-2][l]!=undef)
						ddata[k][j][x-1][l]=(vdata[k][j][0][l]-vdata[k][j][x-2][l])/(dx*2f);
						
					}else if(BCx==BoundaryCondition.Fixed){
						if(vdata[k][j][1  ][l]!=undef&&vdata[k][j][0  ][l]!=undef)
						ddata[k][j][0  ][l]=(vdata[k][j][1  ][l]-vdata[k][j][0  ][l])/dx;
						if(vdata[k][j][x-1][l]!=undef&&vdata[k][j][x-2][l]!=undef)
						ddata[k][j][x-1][l]=(vdata[k][j][x-1][l]-vdata[k][j][x-2][l])/dx;
					}
					
					for(int i=1;i<x-1;i++)
					if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef)
					ddata[k][j][i][l]=(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dx*2f);
				}
			}
			break;
		}
		case Y:{
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					if(BCy==BoundaryCondition.Periodic){
						if(vdata[l][k][1][i]!=undef&&vdata[l][k][y-1][i]!=undef)
						ddata[l][k][0  ][i]=(vdata[l][k][1][i]-vdata[l][k][y-1][i])/(dy*2);
						if(vdata[l][k][0][i]!=undef&&vdata[l][k][y-2][i]!=undef)
						ddata[l][k][y-1][i]=(vdata[l][k][0][i]-vdata[l][k][y-2][i])/(dy*2);
						
					}else if(BCy==BoundaryCondition.Fixed){
						if(vdata[l][k][1][i  ]!=undef&&vdata[l][k][0  ][i]!=undef)
						ddata[l][k][0  ][i]=(vdata[l][k][1  ][i]-vdata[l][k][0  ][i])/dy;
						if(vdata[l][k][1][y-1]!=undef&&vdata[l][k][y-2][i]!=undef)
						ddata[l][k][y-1][i]=(vdata[l][k][y-1][i]-vdata[l][k][y-2][i])/dy;
					}
					
					for(int j=1;j<y-1;j++)
					if(vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
					ddata[l][k][j][i]=(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					if(BCy==BoundaryCondition.Periodic){
						if(vdata[k][1][i][l]!=undef&&vdata[k][y-1][i][l]!=undef)
						ddata[k][0  ][i][l]=(vdata[k][1][i][l]-vdata[k][y-1][i][l])/(dy*2);
						if(vdata[k][0][i][l]!=undef&&vdata[k][y-2][i][l]!=undef)
						ddata[k][y-1][i][l]=(vdata[k][0][i][l]-vdata[k][y-2][i][l])/(dy*2);
						
					}else if(BCy==BoundaryCondition.Fixed){
						if(vdata[k][1  ][i][l]!=undef&&vdata[k][0  ][i][l]!=undef)
						ddata[k][0  ][i][l]=(vdata[k][1  ][i][l]-vdata[k][0  ][i][l])/dy;
						if(vdata[k][y-1][i][l]!=undef&&vdata[k][y-2][i][l]!=undef)
						ddata[k][y-1][i][l]=(vdata[k][y-1][i][l]-vdata[k][y-2][i][l])/dy;
					}
					
					for(int j=1;j<y-1;j++)
					if(vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
					ddata[k][j][i][l]=(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2);
				}
			}
			break;
		}
		case Z:{
			if(var.isTFirst()){
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCz==BoundaryCondition.Periodic){
						if(vdata[l][1][j][i]!=undef&&vdata[l][z-1][j][i]!=undef)
						ddata[l][0  ][j][i]=(vdata[l][1][j][i]-vdata[l][z-1][j][i])/(dz*2);
						if(vdata[l][0][j][i]!=undef&&vdata[l][z-2][j][i]!=undef)
						ddata[l][z-1][j][i]=(vdata[l][0][j][i]-vdata[l][z-2][j][i])/(dz*2);
						
					}else if(BCz==BoundaryCondition.Fixed){
						if(vdata[l][1  ][j][i]!=undef&&vdata[l][0  ][j][i]!=undef)
						ddata[l][0  ][j][i]=(vdata[l][1  ][j][i]-vdata[l][0  ][j][i])/dz;
						if(vdata[l][z-1][j][i]!=undef&&vdata[l][z-2][j][i]!=undef)
						ddata[l][z-1][j][i]=(vdata[l][z-1][j][i]-vdata[l][z-2][j][i])/dz;
					}
					
					for(int k=1;k<z-1;k++)
					if(vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef)
					ddata[l][k][j][i]=(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz*2);
				}
				
			}else{
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCz==BoundaryCondition.Periodic){
						if(vdata[1][j][i][l]!=undef&&vdata[z-1][j][i][l]!=undef)
						ddata[0  ][j][i][l]=(vdata[1][j][i][l]-vdata[z-1][j][i][l])/(dz*2);
						if(vdata[0][j][i][l]!=undef&&vdata[z-2][j][i][l]!=undef)
						ddata[z-1][j][i][l]=(vdata[0][j][i][l]-vdata[z-2][j][i][l])/(dz*2);
						
					}else if(BCz==BoundaryCondition.Fixed){
						if(vdata[1  ][j][i][l]!=undef&&vdata[0  ][j][i][l]!=undef)
						ddata[0  ][j][i][l]=(vdata[1  ][j][i][l]-vdata[0  ][j][i][l])/dz;
						if(vdata[z-1][j][i][l]!=undef&&vdata[z-2][j][i][l]!=undef)
						ddata[z-1][j][i][l]=(vdata[z-1][j][i][l]-vdata[z-2][j][i][l])/dz;
					}
					
					for(int k=1;k<z-1;k++)
					if(vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef)
					ddata[k][j][i][l]=(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz*2);
				}
			}
			break;
		}
		case T:{
			if(var.isTFirst()){
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCt==BoundaryCondition.Periodic){
						if(vdata[1][k][j][i]!=undef&&vdata[t-1][k][j][i]!=undef)
						ddata[0  ][k][j][i]=(vdata[1][k][j][i]-vdata[t-1][k][j][i])/(dt*2);
						if(vdata[0][k][j][i]!=undef&&vdata[t-2][k][j][i]!=undef)
						ddata[t-1][k][j][i]=(vdata[0][k][j][i]-vdata[t-2][k][j][i])/(dt*2);
						
					}else if(BCt==BoundaryCondition.Fixed){
						if(vdata[1  ][k][j][i]!=undef&&vdata[0  ][k][j][i]!=undef)
						ddata[0  ][k][j][i]=(vdata[1  ][k][j][i]-vdata[0  ][k][j][i])/dt;
						if(vdata[t-1][k][j][i]!=undef&&vdata[t-2][k][j][i]!=undef)
						ddata[t-1][k][j][i]=(vdata[t-1][k][j][i]-vdata[t-2][k][j][i])/dt;
					}
					
					for(int l=1;l<t-1;l++)
					if(vdata[l+1][k][j][i]!=undef&&vdata[l-1][k][j][i]!=undef)
					ddata[l][k][j][i]=(vdata[l+1][k][j][i]-vdata[l-1][k][j][i])/(dt*2);
				}
				
			}else{
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					if(BCt==BoundaryCondition.Periodic){
						if(vdata[k][j][i][1]!=undef&&vdata[k][j][i][t-1]!=undef)
						ddata[k][j][i][0  ]=(vdata[k][j][i][1]-vdata[k][j][i][t-1])/(dt*2);
						if(vdata[k][j][i][0]!=undef&&vdata[k][j][i][t-2]!=undef)
						ddata[k][j][i][t-1]=(vdata[k][j][i][0]-vdata[k][j][i][t-2])/(dt*2);
						
					}else if(BCt==BoundaryCondition.Fixed){
						if(vdata[k][j][i][1  ]!=undef&&vdata[k][j][i][0  ]!=undef)
						ddata[k][j][i][0  ]=(vdata[k][j][i][1  ]-vdata[k][j][i][0  ])/dt;
						if(vdata[k][j][i][t-1]!=undef&&vdata[k][j][i][t-2]!=undef)
						ddata[k][j][i][t-1]=(vdata[k][j][i][t-1]-vdata[k][j][i][t-2])/dt;
					}
					
					for(int l=1;l<t-1;l++)
					if(vdata[k][j][i][l+1]!=undef&&vdata[k][j][i][l-1]!=undef)
					ddata[k][j][i][l]=(vdata[k][j][i][l+1]-vdata[k][j][i][l-1])/(dt*2);
				}
			}
			break;
		}
		default: throw new IllegalArgumentException("unsupported dimension: "+dim);
		}
		
		return der;
	}
}
