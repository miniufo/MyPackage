/**
 * @(#)DynamicMethodsInCTS.java	1.0 2017.06.21
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.application.EquationInCartesianCoordinate;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.geophysics.Empirical;
import miniufo.geophysics.Empirical.Scheme;
import miniufo.geophysics.atmos.ThermoDynamics;
import static miniufo.diagnosis.SpatialModel.gEarth;


/**
 * Basic analysis methods in Cartesian coordinates
 *
 * @version 1.0, 2017.06.21
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class DynamicMethodsInCTS extends EquationInCartesianCoordinate{
	//
	public static final float HORCON_PARAM=0.2f;
	
	
	/**
     * Constructor
     *
     * @param	csm		initialized by a spatial model in Cartesian coordinates
     */
	public DynamicMethodsInCTS(CartesianSpatialModel csm){ super(csm);}
	
	
	/**
     * Calculate horizontal divergence.
     *
     * @param	u	u-wind (m s^-1)
     * @param	v	v-wind (m s^-1)
     *
     * @return	divergence (s^-1)
     */
	public Variable c2DDivergence(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable div=new Variable("div",u);
		div.setValue(undef);
		div.setCommentAndUnit("divergence (s^-1)");
		
		float[][][][]   udata=u.getData();
		float[][][][]   vdata=v.getData();
		float[][][][] divdata=div.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east boundary ***/
					if(udata[l][k][j][1]!=undef&&udata[l][k][j][x-1]!=undef&&vdata[l][k][j+1][0]!=undef&&vdata[l][k][j-1][0]!=undef)
						divdata[l][k][j][0]=
						(udata[l][k][j  ][1]-udata[l][k][j][x-1])/(dx*2f)+
						(vdata[l][k][j+1][0]-vdata[l][k][j-1][0])/(dy*2f);
					
					/*** west boundary ***/
					if(udata[l][k][j][0]!=undef&&udata[l][k][j][x-2]!=undef&&vdata[l][k][j+1][x-1]!=undef&&vdata[l][k][j-1][x-1]!=undef)
						divdata[l][k][j][x-1]=
						(udata[l][k][j  ][0  ]-udata[l][k][j  ][x-2])/(dx*2f)+
						(vdata[l][k][j+1][x-1]-vdata[l][k][j-1][x-1])/(dy*2f);
				}
				
				for(int i=1;i<x-1;i++)
				if(udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
					divdata[l][k][j][i]=
					(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dx*2f)+
					(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2f);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east boundary ***/
					if(udata[k][j][1][l]!=undef&&udata[k][j][x-1][l]!=undef&&vdata[k][j+1][0][l]!=undef&&vdata[k][j-1][0][l]!=undef)
						divdata[k][j][0][l]=
						(udata[k][j  ][1][l]-udata[k][j][x-1][l])/(dx*2f)+
						(vdata[k][j+1][0][l]-vdata[k][j-1][0][l])/(dy*2f);
					
					/*** west boundary ***/
					if(udata[k][j][0][l]!=undef&&udata[k][j][x-2][l]!=undef&&vdata[k][j+1][x-1][l]!=undef&&vdata[k][j-1][x-1][l]!=undef)
						divdata[k][j][x-1][l]=
						(udata[k][j  ][0  ][l]-udata[k][j  ][x-2][l])/(dx*2f)+
						(vdata[k][j+1][x-1][l]-vdata[k][j-1][x-1][l])/(dy*2f);
				}
				
				for(int i=1;i<x-1;i++)
				if(udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
					divdata[k][j][i][l]=
					(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dx*2f)+
					(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2f);
			}
		}
		
		return div;
	}
	
	public Variable c3DDivergence(Variable u,Variable v,Variable w){
		checkDimensions(u,v,w);
		assignSubDomainParams(u);
		
		Variable div=new Variable("div",u);
		div.setCommentAndUnit("divergence (s^-1)");
		div.setValue(undef);
		
		float[][][][]   udata=u.getData();
		float[][][][]   vdata=v.getData();
		float[][][][]   wdata=w.getData();
		float[][][][] divdata=div.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++) divdata[l][k][j][i]=
			(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dx*2f)+
			(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2f)+
			(wdata[l][k+1][j][i]-wdata[l][k-1][j][i])/(dz*2f);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++) divdata[k][j][i][l]=
			(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dx*2f)+
			(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2f)+
			(wdata[k+1][j][i][l]-wdata[k-1][j][i][l])/(dz*2f);
		}
		
		return div;
	}
	
	/**
     * Calculate vertical component of vorticity.
     *
     * @param	u	u-wind (m s^-1)
     * @param	v	v-wind (m s^-1)
     *
     * @return	vorticity (s^-1)
     */
	public Variable c2DVorticity(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable vor=new Variable("vor",u);
		vor.setCommentAndUnit("vorticity (s^-1)");
		vor.setValue(undef);
		
		float[][][][]   udata=u.getData();
		float[][][][]   vdata=v.getData();
		float[][][][] vordata=vor.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east boundary ***/
					if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][x-1]!=undef&&udata[l][k][j+1][0]!=undef&&udata[l][k][j-1][0]!=undef)
						vordata[l][k][j][0]=
						(vdata[l][k][j][1]-vdata[l][k][j][x-1])/(dx*2f)-
						(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2f);
					
					/*** west boundary ***/
					if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-2]!=undef&&udata[l][k][j+1][x-1]!=undef&&udata[l][k][j-1][x-1]!=undef)
						vordata[l][k][j][x-1]=
						(vdata[l][k][j][0]-vdata[l][k][j][x-2])/(dx*2f)-
						(udata[l][k][j+1][x-1]-udata[l][k][j-1][x-1])/(dy*2f);
				}
				
				for(int i=1;i<x-1;i++)
				if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef&&udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef)
					vordata[l][k][j][i]=
					(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dx*2f)-
					(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2f);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east boundary ***/
					if(vdata[k][j][1][l]!=undef&&vdata[k][j][x-1][l]!=undef&&udata[k][j+1][0][l]!=undef&&udata[k][j-1][0][l]!=undef)
						vordata[k][j][0][l]=
						(vdata[k][j  ][1][l]-vdata[k][j][x-1][l])/(dx*2f)-
						(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2f);
					
					/*** west boundary ***/
					if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-2][l]!=undef&&udata[k][j+1][x-1][l]!=undef&&udata[k][j-1][x-1][l]!=undef)
						vordata[k][j][x-1][l]=
						(vdata[k][j][0][l]-vdata[k][j][x-2][l])/(dx*2f)-
						(udata[k][j+1][x-1][l]-udata[k][j-1][x-1][l])/(dy*2f);
				}
				
				for(int i=1;i<x-1;i++)
				if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef&&udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef)
					vordata[k][j][i][l]=
					(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dx*2f)-
					(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2f);
			}
		}
		
		return vor;
	}
	
	/**
     * Calculate the divergence in Y-Z plane.
     *
     * @param	v	v-wind (m s^-1)
     * @param	w	w-wind (Pa s^-1)
     *
     * @return	div divergence (s^-1)
     */
	public Variable cYZDivergence(Variable v,Variable w){
		checkDimensions(v,w);
		assignSubDomainParams(v);
		
		Variable div=new Variable("div",v);
		div.setCommentAndUnit("divergence in Y-Z section (s^-1)");
		div.setValue(undef);
		
		float[][][][]   vdata=v.getData();
		float[][][][]   wdata=w.getData();
		float[][][][] divdata=div.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			if(vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef&&wdata[l][k+1][j][i]!=undef&&wdata[l][k-1][j][i]!=undef)
				divdata[l][k][j][i]=
				(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2f)+
				(wdata[l][k+1][j][i]-wdata[l][k-1][j][i])/(dz*2f);
			
		}else{
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			if(vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef&&wdata[k+1][j][i][l]!=undef&&wdata[k-1][j][i][l]!=undef)
				divdata[k][j][i][l]=
				(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2f)+
				(wdata[k+1][j][i][l]-wdata[k-1][j][i][l])/(dz*2f);
		}
		
		return div;
	}
	
	/**
     * Calculate the vorticity in Y-Z plane.
     *
     * @param	v	v-wind (m s^-1)
     * @param	w	w-wind (m s^-1)
     *
     * @return vor	vorticity (s^-1)
     */
	public Variable cYZVorticity(Variable v,Variable w){
		checkDimensions(v,w);
		assignSubDomainParams(v);
		
		Variable vor=new Variable("vor",v);
		vor.setCommentAndUnit("vorticity in Y-P section (s^-1)");
		vor.setValue(undef);
		
		float[][][][]   vdata=v.getData();
		float[][][][]   wdata=w.getData();
		float[][][][] vordata=vor.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			if(wdata[l][k][j+1][i]!=undef&&wdata[l][k][j-1][i]!=undef&&vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef){
				vordata[l][k][j][i]=
				(wdata[l][k][j+1][i]-wdata[l][k][j-1][i])/(dy*2f)-
				(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			if(wdata[k][j+1][i][l]!=undef&&wdata[k][j-1][i][l]!=undef&&vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef){
				vordata[k][j][i][l]=
				(wdata[k][j+1][i][l]-wdata[k][j-1][i][l])/(dy*2f)-
				(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz);
			}
		}
		
		return vor;
	}
	
	/**
     * Calculate the Laplacian in Y-Z plane of a given Variable
     *
     * @param	v	a given Variable
     *
     * @return	F	Laplacian in Y-Z plane of the given Variable
     */
	public Variable cYZLaplacian(Variable v){
		assignSubDomainParams(v);
		
		Variable F=new Variable("lp",v);
		F.setCommentAndUnit("Y-Z Laplacian of "+v.getName());
		F.setValue(undef);
		
		float[][][][] Fdata=F.getData();
		float[][][][] gdata=v.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++){
				for(int k=1;k<z-1;k++)
				for(int j=1;j<y-1;j++)
				Fdata[l][k][j][i]=(
					(gdata[l][k][j+1][i]-gdata[l][k][j][i])-
					(gdata[l][k][j][i]-gdata[l][k][j-1][i])
				)/dy/dy+(
					(gdata[l][k+1][j][i]-gdata[l][k][j][i])-
					(gdata[l][k][j][i]-gdata[l][k-1][j][i])
				)/dz/dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++){
				for(int k=1;k<z-1;k++)
				for(int j=1;j<y-1;j++)
				Fdata[k][j][i][l]=(
					(gdata[k][j+1][i][l]-gdata[k][j][i][l])-
					(gdata[k][j][i][l]-gdata[k][j-1][i][l])
				)/dy/dy+(
					(gdata[k+1][j][i][l]-gdata[k][j][i][l])-
					(gdata[k][j][i][l]-gdata[k-1][j][i][l])
				)/dz/dz;
			}
		}
		
		return F;
	}
	/**
     * Calculate horizontal tension strain (stretching deformation).
     *
     * @param	v	v-wind (m s^-1)
     * @param	w	w-wind (m s^-1)
     *
     * @return	str	tension strain (s^-1)
     */
	public Variable c2DTensionStrain(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable str=new Variable("tstrn",u);
		str.setCommentAndUnit("tension strain (s^-1)");
		str.setValue(undef);
		
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		float[][][][] sdata=str.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef){
					sdata[l][k][j][i]=
					(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dx*2f)-
					(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2f);
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(udata[l][k][j][1]!=undef&&udata[l][k][j][0]!=undef&&vdata[l][k][j+1][0]!=undef&&vdata[l][k][j-1][0]!=undef)
						sdata[l][k][j][0]=
						(udata[l][k][j  ][1]-udata[l][k][j  ][0])/dx-
						(vdata[l][k][j+1][0]-vdata[l][k][j-1][0])/(dy*2f);
					
					if(udata[l][k][j][x-1]!=undef&&udata[l][k][j][x-2]!=undef&&vdata[l][k][j+1][x-1]!=undef&&vdata[l][k][j-1][x-1]!=undef)
						sdata[l][k][j][x-1]=
						(udata[l][k][j  ][x-1]-udata[l][k][j  ][x-2])/dx-
						(vdata[l][k][j+1][x-1]-vdata[l][k][j-1][x-1])/(dy*2f);
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(udata[l][k][0][i+1]!=undef&&udata[l][k][0][i-1]!=undef&&vdata[l][k][1][i]!=undef&&vdata[l][k][0][i]!=undef)
						sdata[l][k][0][i]=
						(udata[l][k][0][i+1]-udata[l][k][0][i-1])/(dx*2f)-
						(vdata[l][k][1][i  ]-vdata[l][k][0][i  ])/dy;
					
					if(udata[l][k][y-1][i+1]!=undef&&udata[l][k][y-1][i-1]!=undef&&vdata[l][k][y-1][i]!=undef&&vdata[l][k][y-2][i]!=undef)
						sdata[l][k][y-1][i]=
						(udata[l][k][y-1][i+1]-udata[l][k][y-1][i-1])/(dx*2f)-
						(vdata[l][k][y-2][i  ]-vdata[l][k][y-1][i  ])/dy;
				}
				
				/*** corner points ***/
				if(udata[l][k][0][0]!=undef&&udata[l][k][0][1]!=undef&&vdata[l][k][0][0]!=undef&&vdata[l][k][1][0]!=undef)
					sdata[l][k][0][0]=
					(udata[l][k][0][1]-udata[l][k][0][0])/dx-
					(vdata[l][k][1][0]-vdata[l][k][0][0])/dy;
				
				if(udata[l][k][y-1][0]!=undef&&udata[l][k][y-1][1]!=undef&&vdata[l][k][y-1][0]!=undef&&vdata[l][k][y-2][0]!=undef)
					sdata[l][k][y-1][0]=
					(udata[l][k][y-1][1]-udata[l][k][y-1][0])/dx-
					(vdata[l][k][y-1][0]-vdata[l][k][y-2][0])/dy;
				
				if(udata[l][k][y-1][x-1]!=undef&&udata[l][k][y-1][x-2]!=undef&&vdata[l][k][y-1][x-1]!=undef&&vdata[l][k][y-2][x-1]!=undef)
					sdata[l][k][y-1][x-1]=
					(udata[l][k][y-1][x-1]-udata[l][k][y-1][x-2])/dx-
					(vdata[l][k][y-1][x-1]-vdata[l][k][y-2][x-1])/dy;
				
				if(udata[l][k][0][x-1]!=undef&&udata[l][k][0][x-2]!=undef&&vdata[l][k][0][x-1]!=undef&&vdata[l][k][1][x-1]!=undef)
					sdata[l][k][0][x-1]=
					(udata[l][k][0][x-1]-udata[l][k][0][x-2])/dx-
					(vdata[l][k][1][x-1]-vdata[l][k][0][x-1])/dy;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef){
					sdata[k][j][i][l]=
					(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dx*2f)-
					(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2f);
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(udata[k][j][1][l]!=undef&&udata[k][j][0][l]!=undef&&vdata[k][j+1][0][l]!=undef&&vdata[k][j-1][0][l]!=undef)
						sdata[k][j][0][l]=
						(udata[k][j  ][1][l]-udata[k][j  ][0][l])/dx-
						(vdata[k][j+1][0][l]-vdata[k][j-1][0][l])/(dy*2f);
					
					if(udata[k][j][x-1][l]!=undef&&udata[k][j][x-2][l]!=undef&&vdata[k][j+1][x-1][l]!=undef&&vdata[k][j-1][x-1][l]!=undef)
						sdata[k][j][x-1][l]=
						(udata[k][j  ][x-1][l]-udata[k][j  ][x-2][l])/dx-
						(vdata[k][j+1][x-1][l]-vdata[k][j-1][x-1][l])/(dy*2f);
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(udata[k][0][i+1][l]!=undef&&udata[k][0][i-1][l]!=undef&&vdata[k][1][i][l]!=undef&&vdata[k][0][i][l]!=undef)
						sdata[k][0][i][l]=
						(udata[k][0][i+1][l]-udata[k][0][i-1][l])/(dx*2f)-
						(vdata[k][1][i  ][l]-vdata[k][0][i  ][l])/dy;
					
					if(udata[k][y-1][i+1][l]!=undef&&udata[k][y-1][i-1][l]!=undef&&vdata[k][y-1][i][l]!=undef&&vdata[k][y-2][i][l]!=undef)
						sdata[k][y-1][i][l]=
						(udata[k][y-1][i+1][l]-udata[k][y-1][i-1][l])/(dx*2f)-
						(vdata[k][y-2][i  ][l]-vdata[k][y-1][i  ][l])/dy;
				}
				
				/*** corner points ***/
				if(udata[k][0][0][l]!=undef&&udata[k][0][1][l]!=undef&&vdata[k][0][0][l]!=undef&&vdata[k][1][0][l]!=undef)
					sdata[k][0][0][l]=
					(udata[k][0][1][l]-udata[k][0][0][l])/dx-
					(vdata[k][1][0][l]-vdata[k][0][0][l])/dy;
				
				if(udata[k][y-1][0][l]!=undef&&udata[k][y-1][1][l]!=undef&&vdata[k][y-1][0][l]!=undef&&vdata[k][y-2][0][l]!=undef)
					sdata[k][y-1][0][l]=
					(udata[k][y-1][1][l]-udata[k][y-1][0][l])/dx-
					(vdata[k][y-1][0][l]-vdata[k][y-2][0][l])/dy;
				
				if(udata[k][y-1][x-1][l]!=undef&&udata[k][y-1][x-2][l]!=undef&&vdata[k][y-1][x-1][l]!=undef&&vdata[k][y-2][x-1][l]!=undef)
					sdata[k][y-1][x-1][l]=
					(udata[k][y-1][x-1][l]-udata[k][y-1][x-2][l])/dx-
					(vdata[k][y-1][x-1][l]-vdata[k][y-2][x-1][l])/dy;
				
				if(udata[k][0][x-1][l]!=undef&&udata[k][0][x-2][l]!=undef&&vdata[k][0][x-1][l]!=undef&&vdata[k][1][x-1][l]!=undef)
					sdata[k][0][x-1][l]=
					(udata[k][0][x-1][l]-udata[k][0][x-2][l])/dx-
					(vdata[k][1][x-1][l]-vdata[k][0][x-1][l])/dy;
			}
		}
		
		return str;
	}
	
	/**
     * Calculate horizontal shearing strain (shearing deformation).
     *
     * @param	u	u-wind (m s^-1)
     * @param	v	v-wind (m s^-1)
     *
     * @return	shr	shearing strain (s^-1)
     */
	public Variable c2DShearingStrain(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable shr=new Variable("sstrn",u);
		shr.setCommentAndUnit("shearing strain (s^-1)");
		shr.setValue(undef);
		
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		float[][][][] sdata=shr.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef&&udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef){
					sdata[l][k][j][i]=
					(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dx*2f)+
					(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2f);
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][0]!=undef&&udata[l][k][j+1][0]!=undef&&udata[l][k][j-1][0]!=undef)
						sdata[l][k][j][0]=
						(vdata[l][k][j][1]-vdata[l][k][j][0])/dx+
						(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2f);
					
					if(vdata[l][k][j][x-1]!=undef&&vdata[l][k][j][x-2]!=undef&&udata[l][k][j+1][x-1]!=undef&&udata[l][k][j-1][x-1]!=undef)
						sdata[l][k][j][x-1]=
						(vdata[l][k][j  ][x-1]-vdata[l][k][j  ][x-2])/dx+
						(udata[l][k][j+1][x-1]-udata[l][k][j-1][x-1])/(dy*2f);
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(vdata[l][k][0][i+1]!=undef&&vdata[l][k][0][i-1]!=undef&&udata[l][k][1][i]!=undef&&udata[l][k][0][i]!=undef)
						sdata[l][k][0][i]=
						(vdata[l][k][0][i+1]-vdata[l][k][0][i-1])/(dx*2f)+
						(udata[l][k][1][i  ]-udata[l][k][0][i  ])/dy;
					
					if(vdata[l][k][y-1][i+1]!=undef&&vdata[l][k][y-1][i-1]!=undef&&udata[l][k][y-1][i]!=undef&&udata[l][k][y-2][i]!=undef)
						sdata[l][k][y-1][i]=
						(vdata[l][k][y-1][i+1]-vdata[l][k][y-1][i-1])/(dx*2f)+
						(udata[l][k][y-1][i]-udata[l][k][y-2][i])/dy;
				}
				
				/*** corner points ***/
				if(vdata[l][k][0][1]!=undef&&vdata[l][k][0][0]!=undef&&udata[l][k][1][0]!=undef&&udata[l][k][0][0]!=undef)
					sdata[l][k][0][0]=
					(vdata[l][k][0][1]-vdata[l][k][0][0])/dx+
					(udata[l][k][1][0]-udata[l][k][0][0])/dy;
				
				if(vdata[l][k][0][x-1]!=undef&&vdata[l][k][0][x-2]!=undef&&udata[l][k][1][x-1]!=undef&&udata[l][k][0][x-1]!=undef)
					sdata[l][k][0][x-1]=
					(vdata[l][k][0][x-1]-vdata[l][k][0][x-2])/dx+
					(udata[l][k][1][x-1]-udata[l][k][0][x-1])/dy;
				
				if(vdata[l][k][y-1][x-1]!=undef&&vdata[l][k][y-1][x-2]!=undef&&udata[l][k][y-1][x-1]!=undef&&udata[l][k][y-2][x-1]!=undef)
					sdata[l][k][y-1][x-1]=
					(vdata[l][k][y-1][x-1]-vdata[l][k][y-1][x-2])/dx+
					(udata[l][k][y-1][x-1]-udata[l][k][y-2][x-1])/dy;
				
				if(vdata[l][k][y-1][1]!=undef&&vdata[l][k][y-1][0]!=undef&&udata[l][k][y-1][0]!=undef&&udata[l][k][y-2][0]!=undef)
					sdata[l][k][y-1][0]=
					(vdata[l][k][y-1][1]-vdata[l][k][y-1][0])/dx+
					(udata[l][k][y-1][0]-udata[l][k][y-2][0])/dy;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef&&udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef){
					sdata[k][j][i][l]=
					(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dx*2f)+
					(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2f);
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(vdata[k][j][1][l]!=undef&&vdata[k][j][0][l]!=undef&&udata[k][j+1][0][l]!=undef&&udata[k][j-1][0][l]!=undef)
						sdata[k][j][0][l]=
						(vdata[k][j  ][1][l]-vdata[k][j  ][0][l])/dx+
						(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2f);
					
					if(vdata[k][j][x-1][l]!=undef&&vdata[k][j][x-2][l]!=undef&&udata[k][j+1][x-1][l]!=undef&&udata[k][j-1][x-1][l]!=undef)
						sdata[k][j][x-1][l]=
						(vdata[k][j  ][x-1][l]-vdata[k][j  ][x-2][l])/dx+
						(udata[k][j+1][x-1][l]-udata[k][j-1][x-1][l])/(dy*2f);
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(vdata[k][0][i+1][l]!=undef&&vdata[k][0][i-1][l]!=undef&&udata[k][1][i][l]!=undef&&udata[k][0][i][l]!=undef)
						sdata[k][0][i][l]=
						(vdata[k][0][i+1][l]-vdata[k][0][i-1][l])/(dx*2f)+
						(udata[k][1][i  ][l]-udata[k][0][i  ][l])/dy;
					
					if(vdata[k][y-1][i+1][l]!=undef&&vdata[k][y-1][i-1][l]!=undef&&udata[k][y-1][i][l]!=undef&&udata[k][y-2][i][l]!=undef)
						sdata[k][y-1][i][l]=
						(vdata[k][y-1][i+1][l]-vdata[k][y-1][i-1][l])/(dx*2f)+
						(udata[k][y-1][i][l]-udata[k][y-2][i][l])/dy;
				}
				
				/*** corner points ***/
				if(vdata[k][0][1][l]!=undef&&vdata[k][0][0][l]!=undef&&udata[k][1][0][l]!=undef&&udata[k][0][0][l]!=undef)
					sdata[k][0][0][l]=
					(vdata[k][0][1][l]-vdata[k][0][0][l])/dx+
					(udata[k][1][0][l]-udata[k][0][0][l])/dy;
				
				if(vdata[k][0][x-1][l]!=undef&&vdata[k][0][x-2][l]!=undef&&udata[k][1][x-1][l]!=undef&&udata[k][0][x-1][l]!=undef)
					sdata[k][0][x-1][l]=
					(vdata[k][0][x-1][l]-vdata[k][0][x-2][l])/dx+
					(udata[k][1][x-1][l]-udata[k][0][x-1][l])/dy;
				
				if(vdata[k][y-1][x-1][l]!=undef&&vdata[k][y-1][x-2][l]!=undef&&udata[k][y-1][x-1][l]!=undef&&udata[k][y-2][x-1][l]!=undef)
					sdata[k][y-1][x-1][l]=
					(vdata[k][y-1][x-1][l]-vdata[k][y-1][x-2][l])/dx+
					(udata[k][y-1][x-1][l]-udata[k][y-2][x-1][l])/dy;
				
				if(vdata[k][y-1][1][l]!=undef&&vdata[k][y-1][0][l]!=undef&&udata[k][y-1][0][l]!=undef&&udata[k][y-2][0][l]!=undef)
					sdata[k][y-1][0][l]=
					(vdata[k][y-1][1][l]-vdata[k][y-1][0][l])/dx+
					(udata[k][y-1][0][l]-udata[k][y-2][0][l])/dy;
			}
		}
		
		return shr;
	}
	
	/**
     * Calculate Smagorinsky viscosity.
     *
     * @param	u	u-wind (m s^-1)
     * @param	v	v-wind (m s^-1)
     *
     * @return	vis	viscosity (s^-1)
     */
	public Variable cSmagorinskyViscosity(Variable u,Variable v){
		checkDimensions(u,v);
		assignSubDomainParams(u);
		
		Variable vis=new Variable("vis",u);
		vis.setCommentAndUnit("Smagorinsky viscosity (m^2 s^-1)");
		vis.setValue(undef);
		
		float[][][][]   udata=u.getData();
		float[][][][]   vdata=v.getData();
		float[][][][] vsdata=vis.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
				if(udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef
				 &&udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef
				 &&vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef
				 &&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef
				){
					vsdata[l][k][j][i]=HORCON_PARAM*dx*dy/2*(float)Math.sqrt(
						Math.pow((udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dx+dx),2)+
						Math.pow((vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy+dy),2)+
						Math.pow((vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dx+dx)+
						(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy+dy),2)/2
					);
				}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
				if(udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef
				 &&udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef
				 &&vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef
				 &&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef
				){
					vsdata[k][j][i][l]=HORCON_PARAM*dx*dy/2*(float)Math.sqrt(
						Math.pow((udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dx+dx),2)+
						Math.pow((vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy+dy),2)+
						Math.pow((vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dx+dx)+
						(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy+dy),2)/2
					);
				}
		}
		
		return vis;
	}
	
	/**
     * Calculate vertical component of potential vorticity.
     *
     * @param	vor	vertical component of vorticity (s^-1)
     * @param	Th	potential temperature (K)
     *
     * @return	PV	potential vorticity (m^2 K s^-1 kg^-1)
     */
	public Variable cVerticalPotentialVorticity(Variable vor,Variable Th){
		checkDimensions(vor,Th);
		assignSubDomainParams(vor);
		
		Variable pv=new Variable("pv",vor);
		pv.setCommentAndUnit("potential vorticity (m^2 K s^-1 kg^-1)");
		pv.setValue(undef);
		
		float[][][][] vordata=vor.getData();
		float[][][][]  Tedata=Th.getData();
		float[][][][]  pvdata=pv.getData();
		
		if(vor.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vordata[l][k][j][i]!=undef&&Tedata[l][k+1][j][i]!=undef&&Tedata[l][k-1][j][i]!=undef){
				pvdata[l][k][j][i]=-gEarth*vordata[l][k][j][i]*
				(Tedata[l][k+1][j][i]-Tedata[l][k-1][j][i])/(dz+dz);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vordata[k][j][i][l]!=undef&&Tedata[k+1][j][i][l]!=undef&&Tedata[k-1][j][i][l]!=undef){
				pvdata[k][j][i][l]=-gEarth*vordata[k][j][i][l]*
				(Tedata[k+1][j][i][l]-Tedata[k-1][j][i][l])/(dz+dz);
			}
		}
		
		return pv;
	}
	
	
	/**
     * Calculate advection term.
     *
     * @param	a	a given variable
     * @param	u	u-wind (m s^-1)
     * @param	v	v-wind (m s^-1)
     *
     * @return adv	advection term
     */
	public Variable cAdvectionTerm(Variable a,Variable u,Variable v){
		checkDimensions(a,u,v);
		assignSubDomainParams(u);
		
		Variable adv=new Variable(a.getName()+"adv",v);
		adv.setCommentAndUnit("advection term");
		adv.setValue(undef);
		
		float[][][][]  adata= a.getData();
		float[][][][]  udata= u.getData();
		float[][][][]  vdata= v.getData();
		float[][][][] atdata=adv.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			if(udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef&&adata[l][k][j][i+1]!=undef
			 &&adata[l][k][j][i-1]!=undef&&adata[l][k][j+1][i]!=undef&&adata[l][k][j-1][i]!=undef){
				atdata[l][k][j][i]=
				udata[l][k][j][i]*(adata[l][k][j][i+1]-adata[l][k][j][i-1])/(dx*2f)+
				vdata[l][k][j][i]*(adata[l][k][j+1][i]-adata[l][k][j-1][i])/(dy*2f);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			if(udata[k][j][i][l]!=undef&&adata[k][j][i+1][l]!=undef&&adata[k][j][i-1][l]!=undef
			 &&vdata[k][j][i][l]!=undef&&adata[k][j+1][i][l]!=undef&&adata[k][j-1][i][l]!=undef){
				atdata[k][j][i][l]=
				udata[k][j][i][l]*(adata[k][j][i+1][l]-adata[k][j][i-1][l])/(dx*2f)+
				vdata[k][j][i][l]*(adata[k][j+1][i][l]-adata[k][j-1][i][l])/(dy*2f);
			}
		}
		
		return adv;
	}
	
	/**
     * Calculate convection term.
     *
     * @param	v	a given variable
     * @param	w	w-wind  (m s^-1 or Pa s^-1)
     *
     * @return	cov convection term
     */
	public Variable cConvectiveTerm(Variable v,Variable w){
		checkDimensions(v,w);
		assignSubDomainParams(v);
		
		Variable cov=new Variable(v.getName()+"con",v);	cov.setCommentAndUnit("convection term");
		
		float[][][][]  vdata=  v.getData();
		float[][][][]  wdata=  w.getData();
		float[][][][] ctdata=cov.getData();
		
		if(v.isTFirst()){
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[l][k][j][i]!=undef&&vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef){
					if(ctdata[l][k][j][i]!=undef)
						ctdata[l][k][j][i]=wdata[l][k][j][i]*(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz);
					
					else ctdata[l][k][j][i]=wdata[l][k][j][i]*(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[l][0][j][i]!=undef&&vdata[l][1][j][i]!=undef&&vdata[l][0][j][i]!=undef){
					if(ctdata[l][0][j][i]!=undef)
						ctdata[l][0][j][i]+=wdata[l][0][j][i]*(vdata[l][1][j][i]-vdata[l][0][j][i])/dz;
					
					else ctdata[l][0][j][i]=wdata[l][0][j][i]*(vdata[l][1][j][i]-vdata[l][0][j][i])/dz;
				}
				
				if(wdata[l][z-1][j][i]!=undef&&vdata[l][z-1][j][i]!=undef&&vdata[l][z-2][j][i]!=undef){
					if(ctdata[l][z-1][j][i]!=undef)
						ctdata[l][z-1][j][i]+=wdata[l][z-1][j][i]*(vdata[l][z-1][j][i]-vdata[l][z-2][j][i])/dz;
					
					else ctdata[l][z-1][j][i]=wdata[l][z-1][j][i]*(vdata[l][z-1][j][i]-vdata[l][z-2][j][i])/dz;
				}
			}
			
		}else{
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
				if(wdata[k][j][i][l]!=undef&&vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef){
					if(ctdata[k][j][i][l]!=undef)
						ctdata[k][j][i][l]=wdata[k][j][i][l]*(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz);
					
					else ctdata[k][j][i][l]=wdata[k][j][i][l]*(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz);
				}
			
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(wdata[0][j][i][l]!=undef&&vdata[1][j][i][l]!=undef&&vdata[0][j][i][l]!=undef){
					if(ctdata[0][j][i][l]!=undef)
						ctdata[0][j][i][l]+=wdata[0][j][i][l]*(vdata[1][j][i][l]-vdata[0][j][i][l])/dz;
					
					else ctdata[0][j][i][l]=wdata[0][j][i][l]*(vdata[1][j][i][l]-vdata[0][j][i][l])/dz;
				}
				
				if(wdata[z-1][j][i][l]!=undef&&vdata[z-1][j][i][l]!=undef&&vdata[z-2][j][i][l]!=undef){
					if(ctdata[z-1][j][i][l]!=undef)
						ctdata[z-1][j][i][l]+=wdata[z-1][j][i][l]*(vdata[z-1][j][i][l]-vdata[z-2][j][i][l])/dz;
					
					else ctdata[z-1][j][i][l]=wdata[z-1][j][i][l]*(vdata[z-1][j][i][l]-vdata[z-2][j][i][l])/dz;
				}
			}
		}
		
		return cov;
	}
	
	/**
     * Calculate the Laplacian of a given Variable.
     *
     * @param	v	a given Variable
     *
     * @return	F	Laplacian of the given Variable
     */
	public Variable cLaplacian(Variable v){
		assignSubDomainParams(v);
		
		Variable F=new Variable("lp",v);
		F.setCommentAndUnit("Laplacian of "+v.getName());
		F.setValue(undef);
		
		float[][][][] Fdata=F.getData();
		float[][][][] gdata=v.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east BC with i=0 ***/
					Fdata[l][k][j][0]=(
						(gdata[l][k][j][1]-gdata[l][k][j][0  ])-
						(gdata[l][k][j][0]-gdata[l][k][j][x-1])
					)/dx/dx+(
						(gdata[l][k][j+1][0]-gdata[l][k][j][0])-
						(gdata[l][k][j][0]-gdata[l][k][j-1][0])
					)/dy/dy;
					
					/*** east BC with i=x-1 ***/
					Fdata[l][k][j][x-1]=(
						(gdata[l][k][j][0  ]-gdata[l][k][j][x-1])-
						(gdata[l][k][j][x-1]-gdata[l][k][j][x-2])
					)/dx/dx+(
						(gdata[l][k][j+1][x-1]-gdata[l][k][j][x-1])-
						(gdata[l][k][j][x-1]-gdata[l][k][j-1][x-1])
					)/dy/dy;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=(
					(gdata[l][k][j][i+1]-gdata[l][k][j][i])-
					(gdata[l][k][j][i]-gdata[l][k][j][i-1])
				)/dx/dx+(
					(gdata[l][k][j+1][i]-gdata[l][k][j][i])-
					(gdata[l][k][j][i]-gdata[l][k][j-1][i])
				)/dy/dy;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east BC with i=0 ***/
					Fdata[k][j][0][l]=(
						(gdata[k][j][1][l]-gdata[k][j][0  ][l])-
						(gdata[k][j][0][l]-gdata[k][j][x-1][l])
					)/dx/dx+(
						(gdata[k][j+1][0][l]-gdata[k][j][0][l])-
						(gdata[k][j][0][l]-gdata[k][j-1][0][l])
					)/dy/dy;
					
					/*** east BC with i=x-1 ***/
					Fdata[k][j][x-1][l]=(
						(gdata[k][j][0  ][l]-gdata[k][j][x-1][l])-
						(gdata[k][j][x-1][l]-gdata[k][j][x-2][l])
					)/dx/dx+(
						(gdata[k][j+1][x-1][l]-gdata[k][j][x-1][l])-
						(gdata[k][j][x-1][l]-gdata[k][j-1][x-1][l])
					)/dy/dy;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=(
					(gdata[k][j][i+1][l]-gdata[k][j][i][l])-
					(gdata[k][j][i][l]-gdata[k][j][i-1][l])
				)/dx/dx+(
					(gdata[k][j+1][i][l]-gdata[k][j][i][l])-
					(gdata[k][j][i][l]-gdata[k][j-1][i][l])
				)/dy/dy;
			}
		}
		
		return F;
	}
	
	/**
     * Calculate the concentration given number of observations.
     *
     * @param	count	number of observations (1)
	 * @param	dlon	delta lon for periodic boundary condition
     *
     * @return	conc	concentration (m^-2)
     */
	public Variable cConcentration(Variable count,float dlon){
		assignSubDomainParams(count);
		
		Variable conc=new Variable("conc",count);
		conc.setCommentAndUnit("concentration");
		
		float[][][][] ctdata=count.getData();
		float[][][][] ccdata=conc.getData();
		
	    if(conc.isTFirst()){
			for(int j=0;j<y;j++){
				float deltaY=(j==0||j==y-1)?dy/2f:dy;
				
				for(int i=0;i<x;i++){
					float deltaX=(i==0||i==x-1)?dx/2f:dx;
					
					for(int l=0;l<t;l++)
					for(int k=0;k<z;k++)
					if(ccdata[l][k][j][i]!=undef) ccdata[l][k][j][i]=ctdata[l][k][j][i]/(deltaY*deltaX);
					else ccdata[l][k][j][i]=undef;
				}
			}
			
		}else{
			for(int j=0;j<y;j++){
				float deltaY=(j==0||j==y-1)?dy/2f:dy;
				
				for(int i=0;i<x;i++){
					float deltaX=(i==0||i==x-1)?dx/2f:dx;
					
					for(int l=0;l<t;l++)
					for(int k=0;k<z;k++)
					if(ccdata[k][j][i][l]!=undef) ccdata[k][j][i][l]=ctdata[k][j][i][l]/(deltaY*deltaX);
					else ccdata[k][j][i][l]=undef;
				}
			}
		}
		
		return conc;
	}
	
	/**
     * Calculate helicity.
     *
     * @param	u	u-wind (m s^-1)
     * @param	v	v-wind (m s^-1)
     * @param	w	w-wind (m s^-1)
     *
     * @return	h	helicity ()
     */
	public Variable cHelicity(Variable u,Variable v,Variable w){
		checkDimensions(u,v,w);
		assignSubDomainParams(u);
		
		Variable h=new Variable("hel",u);
		h.setCommentAndUnit("helicity ()");
		h.setValue(undef);
		
		float[][][][] hdata=h.getData();
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		float[][][][] wdata=w.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			if(
				udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef&&wdata[l][k][j][i]!=undef&&
				udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef&&udata[l][k+1][j][i]!=undef&&udata[l][k-1][j][i]!=undef&&
				vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef&&vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef&&
				wdata[l][k][j][i+1]!=undef&&wdata[l][k][j][i-1]!=undef&&wdata[l][k][j+1][i]!=undef&&wdata[l][k][j-1][i]!=undef
			){
				hdata[l][k][j][i]=
				udata[l][k][j][i]*(
					(wdata[l][k][j+1][i]-wdata[l][k][j-1][i])/(2f*dy)-
					(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz)
				)+
				vdata[l][k][j][i]*(
					(udata[l][k+1][j][i]-udata[l][k-1][j][i])/(dz+dz)-
					(wdata[l][k][j][i+1]-wdata[l][k][j][i-1])/(2f*dx)
				)+
				wdata[l][k][j][i]*(
					(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(2f*dx)-
					(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(2f*dy)
				);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++)
			if(
				udata[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef&&wdata[k][j][i][l]!=undef&&
				udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef&&udata[k+1][j][i][l]!=undef&&udata[k-1][j][i][l]!=undef&&
				vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef&&vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef&&
				wdata[k][j][i+1][l]!=undef&&wdata[k][j][i-1][l]!=undef&&wdata[k][j+1][i][l]!=undef&&wdata[k][j-1][i][l]!=undef
			){
				hdata[k][j][i][l]=
				udata[k][j][i][l]*(
					(wdata[k][j+1][i][l]-wdata[k][j-1][i][l])/(2f*dy)-
					(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz)
				)+
				vdata[k][j][i][l]*(
					(udata[k+1][j][i][l]-udata[k-1][j][i][l])/(dz+dz)-
					(wdata[k][j][i+1][l]-wdata[k][j][i-1][l])/(2f*dx)
				)+
				wdata[k][j][i][l]*(
					(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(2f*dx)-
					(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(2f*dy)
				);
			}
		}
		
		return h;
	}
	
	/**
     * Calculate vertical velocity according to the continuity equation.
     *
     * @param	div		horizontal divergence (s^-1)
     *
     * @return	w		vertical velocity in pressure coordinate (Pa s^-1)
     */
	public Variable cOmega(Variable div){
		assignSubDomainParams(div);
		
		Variable w=new Variable("omega",div);
		w.setCommentAndUnit("vertical velocity (Pa s^-1)");
		w.setValue(undef);
		
		float[][][][] wdata  =w.getData();
		float[][][][] divdata=div.getData();
		
		if(div.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(divdata[l][k-1][j][i]!=undef&&divdata[l][k][j][i]!=undef){
				if(wdata[l][k-1][j][i]!=undef)
					wdata[l][k][j][i]=wdata[l][k-1][j][i]+
					(divdata[l][k][j][i]+divdata[l][k-1][j][i])/2*dz;
				else
					wdata[l][k][j][i]=
					(divdata[l][k][j][i]+divdata[l][k-1][j][i])/2*dz;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(divdata[i][j][k-1][l]!=undef&&divdata[k][j][i][l]!=undef){
				if(wdata[k-1][j][i][l]!=undef)
					wdata[k][j][i][l]=wdata[k-1][j][i][l]+
					(divdata[k][j][i][l]+divdata[k-1][j][i][l])/2*dz;
				else
					wdata[k][j][i][l]=
					(divdata[k][j][i][l]+divdata[k-1][j][i][l])/2*dz;
			}
		}
		
		return w;
	}
	
	
	/**
     * Correct vertical velocity through error redistribution.
     *
     * @param	correct_type	different type of correction
     * @param	w				vertical velocity (Pa s^-1)
     */
	public void correctOmega(int correct_type,Variable w){
		assignSubDomainParams(w);
		
		SphericalSpatialModel ssm=(SphericalSpatialModel)sm;
		
		float[] max=new float[t];
		float[]  a =new float[t];
		float[][][] tmp_top_w=new float[t][y][x];
		float[][][][] wdata=w.getData();
		
		if(w.isTFirst()){
			for(int l=0;l<t;l++){
				a[l]=1;
				max[l]=Float.MIN_VALUE;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					float b=(wdata[l][z-1][j][i]<0)?(-wdata[l][z-1][j][i]):wdata[l][z-1][j][i];
					if(b>max[l]) max[l]=b;
				}
				
				while(max[l]>=0.1){ a[l]*=0.1f; max[l]*=a[l];}
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					tmp_top_w[l][j][i]=wdata[l][z-1][j][i];
					wdata[l][z-1][j][i]*=a[l];
				}
			}
			
			float[] sum={0f,0f};
			switch(correct_type){
				case 1:
					for(int k=0;k<z-1;k++) sum[1]+=dz;
					
					for(int l=0;l<t;l++){
						for(int k=0;k<z-1;k++){
							sum[0]+=dz*100;
							
							for(int j=0;j<y;j++)
							for(int i=0;i<x;i++)
							wdata[l][k][j][i]-=(1-a[l])*tmp_top_w[l][j][i]*sum[0]/
								(zdef[zstart-1]-zdef[zstart-2+ssm.getZCount()]);
						}
						
						sum[0]=0;
					}
				
					break;
					
				case 2:
					for(int k=0;k<z-1;k++) sum[1]+=k*dz;
					
					for(int l=0;l<t;l++){
						for(int k=0;k<z-1;k++){
							sum[0]+=k*dz;
							
							for(int j=0;j<y;j++)
							for(int i=0;i<x;i++)
							wdata[l][k][j][i]-=(1-a[l])*tmp_top_w[l][j][i]*sum[0]/sum[1];
						}
							
						sum[0]=0;
					}
					
					break;
				
				default:
					for(int k=0;k<z-1;k++) sum[1]+=(1-2*k)*dz;
					
					for(int l=0;l<t;l++){
						for(int k=0;k<z-1;k++){
							sum[0]+=(1-2*k)*dz;
							
							for(int j=0;j<y;j++)
							for(int i=0;i<x;i++)
							wdata[l][k][j][i]-=(1-a[l])*tmp_top_w[l][j][i]*sum[0]/sum[1];
						}
						
						sum[0]=0;
					}
			}
			
		}else{
			for(int l=0;l<t;l++){
				a[l]=1;
				max[l]=Float.MIN_VALUE;
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					float b=(wdata[z-1][j][i][l]<0)?(-wdata[z-1][j][i][l]):wdata[z-1][j][i][l];
					if(b>max[l]) max[l]=b;
				}
				
				while(max[l]>=0.1){ a[l]*=0.1f; max[l]*=a[l];}
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					tmp_top_w[j][i][l]=wdata[z-1][j][i][l];
					wdata[z-1][j][i][l]*=a[l];
				}
			}
			
			float[] sum={0f,0f};
			switch(correct_type){
				case 1:
					for(int k=0;k<z-1;k++) sum[1]+=dz;
					
					for(int l=0;l<t;l++){
						for(int k=0;k<z-1;k++){
							sum[0]+=dz;
							
							for(int j=0;j<y;j++)
							for(int i=0;i<x;i++)
							wdata[k][j][i][l]-=(1-a[l])*tmp_top_w[j][i][l]*sum[0]/
								(zdef[zstart-1]-zdef[zstart-2+ssm.getZCount()]);
						}
						
						sum[0]=0;
					}
				
					break;
					
				case 2:
					for(int k=0;k<z-1;k++) sum[1]+=(k+1)*dz;
					
					for(int l=0;l<t;l++){
						for(int k=0;k<z-1;k++){
							sum[0]+=(k+1)*dz;
							
							for(int j=0;j<y;j++)
							for(int i=0;i<x;i++)
							wdata[k][j][i][l]-=(1-a[l])*tmp_top_w[j][i][l]*sum[0]/sum[1];
						}
							
						sum[0]=0;
					}
					
					break;
				
				default:
					for(int k=0;k<z-1;k++) sum[1]+=(1-2*(k+1))*dz;
					
					for(int l=0;l<t;l++){
						for(int k=0;k<z-1;k++){
							sum[0]+=(1-2*(k+1))*dz;
							
							for(int j=0;j<y;j++)
							for(int i=0;i<x;i++)
							wdata[k][j][i][l]-=(1-a[l])*tmp_top_w[j][i][l]*sum[0]/sum[1];
						}
						
						sum[0]=0;
					}
			}
		}
	}
	
	
	/**
     * Calculate the magnitude of horizontal gradient.
     *
     * @param	v		a given Variable
     *
     * @return	grd		gradient magnitude
     */
	public Variable c2DGradientMagnitude(Variable v){
		Variable grdx=cDerivative(v,Dimension.X);
		Variable grdy=cDerivative(v,Dimension.Y);
		
		Variable mag=grdx.hypotenuseEq(grdy);
		
		mag.setName("grdmag");
		mag.setCommentAndUnit("magnitude of "+v.getName()+" gradient ("+v.getUnit()+" m^-1)");
		
		return mag;
	}
	
	/**
     * Calculate horizontal gradient, [0] is x- and [1] y-components.
     *
     * @param	v		a given Variable
     *
     * @return	grd		2D gradient (vector)
     */
	public Variable[] c2DGradient(Variable v){
		Variable grdx=cDerivative(v,Dimension.X);
		Variable grdy=cDerivative(v,Dimension.Y);
		
		grdx.setCommentAndUnit("gradient of "+v.getName()+" in x-direction");
		grdy.setCommentAndUnit("gradient of "+v.getName()+" in y-direction");
		
		return new Variable[]{grdx,grdy};
	}
	
	/**
     * Calculate unit vector of horizontal gradient, [0] is x- and [1] y-components.
     *
     * @param	v		a given Variable
     *
     * @return	grd		2D unit gradient vector
     */
	public Variable[] c2DUnitGradient(Variable v){
		assignSubDomainParams(v);
		
		Variable grdx=cDerivative(v,Dimension.X);
		Variable grdy=cDerivative(v,Dimension.Y);
		
		float[][][][] xdata=grdx.getData();
		float[][][][] ydata=grdy.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(xdata[l][k][j][i]!=undef&&ydata[l][k][j][i]!=undef){
				double mag=Math.hypot(xdata[l][k][j][i],ydata[l][k][j][i]);
				xdata[l][k][j][i]/=mag;
				ydata[l][k][j][i]/=mag;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(xdata[k][j][i][l]!=undef&&ydata[k][j][i][l]!=undef){
				double mag=Math.hypot(xdata[k][j][i][l],ydata[k][j][i][l]);
				xdata[k][j][i][l]/=mag;
				ydata[k][j][i][l]/=mag;
			}
		}
		
		grdx.setCommentAndUnit("gradient of "+v.getName()+" in x-direction");
		grdy.setCommentAndUnit("gradient of "+v.getName()+" in y-direction");
		
		return new Variable[]{grdx,grdy};
	}
	
	/**
     * Calculate 3D gradient, [0] is x-, [1] y- and [2] z-components.
     *
     * @param	v		a given Variable
     *
     * @return	grd		3D gradient (vector)
     */
	public Variable[] c3DGradient(Variable v){
		Variable grdx=cDerivative(v,Dimension.X);
		Variable grdy=cDerivative(v,Dimension.Y);
		Variable grdz=cDerivative(v,Dimension.Z);
		
		grdx.setCommentAndUnit("gradient of "+v.getName()+" in x-direction");
		grdy.setCommentAndUnit("gradient of "+v.getName()+" in y-direction");
		grdz.setCommentAndUnit("gradient of "+v.getName()+" in z-direction");
		
		return new Variable[]{grdx,grdy,grdz};
	}
	
	
	/**
     * Calculate wind stress using the 10-m wind.
     *
     * @param	uwnd	zonal wind (m s^-1)
     * @param	vwnd	meridional wind (m s^-1)
     * @param	scheme	which scheme is used [LP, DO, LY]
     *
     * @return	gradient (vector)
     */
	public Variable[] cWindStress(Variable uwnd,Variable vwnd,Scheme scheme){
		checkDimensions(uwnd,vwnd);
		assignSubDomainParams(uwnd);
		
		Variable taux=new Variable("taux",uwnd);
		Variable tauy=new Variable("tauy",vwnd);
		taux.setValue(undef);	taux.setCommentAndUnit("zonal wind stress (N)");
		tauy.setValue(undef);	tauy.setCommentAndUnit("meridianal wind stress (N)");
		
		float[][][][]  udata=uwnd.getData();
		float[][][][]  vdata=vwnd.getData();
		float[][][][]  xdata=taux.getData();
		float[][][][]  ydata=tauy.getData();
		
		if(uwnd.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float wndSpd=(float)Math.hypot(udata[l][k][j][i],vdata[l][k][j][i]);
				float Cd=Empirical.cBulkDragCoefficient(wndSpd,scheme);
				float tmp=ThermoDynamics.RHO*Cd*wndSpd;
				
				xdata[l][k][j][i]=tmp*udata[l][k][j][i];
				ydata[l][k][j][i]=tmp*vdata[l][k][j][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float wndSpd=(float)Math.hypot(udata[k][j][i][l],vdata[k][j][i][l]);
				float Cd=Empirical.cBulkDragCoefficient(wndSpd,scheme);
				float tmp=ThermoDynamics.RHO*Cd*wndSpd;
				
				xdata[k][j][i][l]=tmp*udata[k][j][i][l];
				ydata[k][j][i][l]=tmp*vdata[k][j][i][l];
			}
		}
		
		return new Variable[]{taux,tauy};
	}
	
	
	/** test
	public static void main(String arg[]){
		DiagnosisFactory dfu=DiagnosisFactory.parseFile("d:/Data/NCEP/OriginalNC/uwnd.sig995.mon.mean.nc");
		DataDescriptor dd=dfu.getDataDescriptor();
		
		Variable u=dfu.getVariables(new Range("t(1,60)",dd),"uwnd")[0];
		Variable v=DiagnosisFactory.getVariables("d:/Data/NCEP/OriginalNC/vwnd.sig995.mon.mean.nc","t(1,60)","vwnd")[0];
		
		SphericalSpacialModel ssm=new SphericalSpacialModel(dd);
		DynamicMethodsInSC dm=new DynamicMethodsInSC(ssm);
		
		Variable[] ek=dm.cEkmanCurrent(u,v);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Ekman.dat");
		dw.writeData(dd,u,v,ek[0],ek[1]);	dw.closeFile();
	}*/
}
