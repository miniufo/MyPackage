/**
 * @(#)DynamicMethodsInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.application.EquationInSphericalCoordinate;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.geophysics.Empirical;
import miniufo.geophysics.Empirical.Scheme;
import miniufo.geophysics.atmos.ThermoDynamics;
import static miniufo.diagnosis.SpatialModel.GRAVITY_ACCERLERATION;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


/**
 * Basic analysis methods in spherical coordinates.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class DynamicMethodsInSC extends EquationInSphericalCoordinate{
	//
	public static final float HORCON_PARAM=0.2f;
	
	
	/**
     * Constructor
     *
     * @param	ssm		initialized by a spatial model in spherical coordinates
     */
	public DynamicMethodsInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
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
		
		Variable div=new Variable("div",u);	div.setCommentAndUnit("divergence (s^-1)");
		
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
							(udata[l][k][j][1]-udata[l][k][j][x-1])/(dxs[ystart-1+j]*2)+
							(vdata[l][k][j+1][0]-vdata[l][k][j-1][0])/(dy*2)-
							vdata[l][k][j][0]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else divdata[l][k][j][0]=undef;
					
					/*** west boundary ***/
					if(udata[l][k][j][0]!=undef&&udata[l][k][j][x-2]!=undef&&vdata[l][k][j+1][x-1]!=undef&&vdata[l][k][j-1][x-1]!=undef)
						divdata[l][k][j][x-1]=
							(udata[l][k][j][0]-udata[l][k][j][x-2])/(dxs[ystart-1+j]*2)+
							(vdata[l][k][j+1][x-1]-vdata[l][k][j-1][x-1])/(dy*2)-
							vdata[l][k][j][x-1]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else divdata[l][k][j][x-1]=undef;
				}
				
				for(int i=1;i<x-1;i++)
				if(udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
					divdata[l][k][j][i]=
						(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
						(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
						vdata[l][k][j][i]*ltan[ystart-1+j]/EARTH_RADIUS;
					
				else divdata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east boundary ***/
					if(udata[k][j][1][l]!=undef&&udata[k][j][x-1][l]!=undef&&vdata[k][j+1][0][l]!=undef&&vdata[k][j-1][0][l]!=undef)
						divdata[k][j][0][l]=
							(udata[k][j][1][l]-udata[k][j][x-1][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][0][l]-vdata[k][j-1][0][l])/(dy*2)-
							vdata[k][j][0][l]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else divdata[l][k][j][0]=undef;
					
					/*** west boundary ***/
					if(udata[k][j][0][l]!=undef&&udata[k][j][x-2][l]!=undef&&vdata[k][j+1][x-1][l]!=undef&&vdata[k][j-1][x-1][l]!=undef)
						divdata[k][j][x-1][l]=
							(udata[k][j][0][l]-udata[k][j][x-2][l])/(dxs[ystart-1+j]*2)+
							(vdata[k][j+1][x-1][l]-vdata[k][j-1][x-1][l])/(dy*2)-
							vdata[k][j][x-1][l]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else divdata[l][k][j][x-1]=undef;
				}
				
				for(int i=1;i<x-1;i++)
				if(udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
					divdata[k][j][i][l]=
						(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
						(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
						(vdata[k][j][i][l]/EARTH_RADIUS)*ltan[ystart-1+j];
					
				else divdata[k][j][i][l]=undef;
			}
		}
		
		return div;
	}
	
	public Variable c3DDivergence(Variable u,Variable v,Variable w){
		checkDimensions(u,v,w);
		assignSubDomainParams(u);
		
		Variable div=new Variable("div",u);	div.setCommentAndUnit("divergence (s^-1)");
		
		float[][][][]   udata=u.getData();
		float[][][][]   vdata=v.getData();
		float[][][][]   wdata=w.getData();
		float[][][][] divdata=div.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				divdata[l][k][j][i]=
					(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
					(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)-
					vdata[l][k][j][i]*ltan[ystart-1+j]/EARTH_RADIUS+
					(wdata[l][k+1][j][i]-wdata[l][k-1][j][i])/(dz*2);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				divdata[k][j][i][l]=
					(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
					(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)-
					(vdata[k][j][i][l]/EARTH_RADIUS)*ltan[ystart-1+j]+
					(wdata[k+1][j][i][l]-wdata[k-1][j][i][l])/(dz*2);
			}
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
		
		Variable vor=new Variable("vor",u);	vor.setCommentAndUnit("vorticity (s^-1)");
		
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
							(vdata[l][k][j][1]-vdata[l][k][j][x-1])/(dxs[ystart-1+j]*2)-
							(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2)+
							udata[l][k][j][0]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else vordata[l][k][j][0]=undef;
					
					/*** west boundary ***/
					if(vdata[l][k][j][0]!=undef&&vdata[l][k][j][x-2]!=undef&&udata[l][k][j+1][x-1]!=undef&&udata[l][k][j-1][x-1]!=undef)
						vordata[l][k][j][x-1]=
							(vdata[l][k][j][0]-vdata[l][k][j][x-2])/(dxs[ystart-1+j]*2)-
							(udata[l][k][j+1][x-1]-udata[l][k][j-1][x-1])/(dy*2)+
							udata[l][k][j][x-1]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else vordata[l][k][j][x-1]=undef;
				}
				
				for(int i=1;i<x-1;i++)
				if(vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef&&udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef)
					vordata[l][k][j][i]=
						(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)-
						(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)+
						udata[l][k][j][i]*ltan[ystart-1+j]/EARTH_RADIUS;
					
				else vordata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				if(BCx==BoundaryCondition.Periodic){
					/*** east boundary ***/
					if(vdata[k][j][1][l]!=undef&&vdata[k][j][x-1][l]!=undef&&udata[k][j+1][0][l]!=undef&&udata[k][j-1][0][l]!=undef)
						vordata[k][j][0][l]=
							(vdata[k][j][1][l]-vdata[k][j][x-1][l])/(dxs[ystart-1+j]*2)-
							(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2)+
							udata[k][j][0][l]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else vordata[l][k][j][0]=undef;
					
					/*** west boundary ***/
					if(vdata[k][j][0][l]!=undef&&vdata[k][j][x-2][l]!=undef&&udata[k][j+1][x-1][l]!=undef&&udata[k][j-1][x-1][l]!=undef)
						vordata[k][j][x-1][l]=
							(vdata[k][j][0][l]-vdata[k][j][x-2][l])/(dxs[ystart-1+j]*2)-
							(udata[k][j+1][x-1][l]-udata[k][j-1][x-1][l])/(dy*2)+
							udata[k][j][x-1][l]*ltan[ystart-1+j]/EARTH_RADIUS;
						
					else vordata[l][k][j][x-1]=undef;
				}
				
				for(int i=1;i<x-1;i++)
				if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef&&udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef)
					vordata[k][j][i][l]=
						(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)-
						(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2)+
						udata[k][j][i][l]*ltan[ystart-1+j]/EARTH_RADIUS;
					
				else vordata[k][j][i][l]=undef;
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
		
		Variable div=new Variable("div",v);	div.setCommentAndUnit("divergence in Y-Z section");
		
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
				(vdata[l][k][j+1][i]*lcos[ystart+j]-vdata[l][k][j-1][i]*lcos[ystart-2+j])/(dy*2f)/lcos[ystart-1+j]+
				(wdata[l][k+1][j][i]-wdata[l][k-1][j][i])/(dz*2f);
			
		}else{
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			if(vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef&&wdata[k+1][j][i][l]!=undef&&wdata[k-1][j][i][l]!=undef)
				divdata[k][j][i][l]=
				(vdata[k][j+1][i][l]*lcos[ystart+j]-vdata[k][j-1][i][l]*lcos[ystart-2+j])/(dy*2f)/lcos[ystart-1+j]+
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
		
		Variable vor=new Variable("vor",v);	vor.setCommentAndUnit("vorticity in Y-P section");
		
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
				(wdata[l][k][j+1][i]-wdata[l][k][j-1][i])/(dy+dy)-
				(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz)+
				wdata[l][k][j][i]*ltan[ystart-1+j]/EARTH_RADIUS;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int i=0;i<x;i++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			if(wdata[k][j+1][i][l]!=undef&&wdata[k][j-1][i][l]!=undef&&vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef){
				vordata[k][j][i][l]=
				(wdata[k][j+1][i][l]-wdata[k][j-1][i][l])/(dy+dy)-
				(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz)+
				wdata[k][j][i][l]*ltan[ystart-1+j]/EARTH_RADIUS;
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
					(gdata[l][k][j+1][i]-gdata[l][k][j][i])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
					(gdata[l][k][j][i]-gdata[l][k][j-1][i])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
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
					(gdata[k][j+1][i][l]-gdata[k][j][i][l])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
					(gdata[k][j][i][l]-gdata[k][j-1][i][l])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
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
				for(int i=1;i<x-1;i++){
					if(udata[l][k][j][i+1]!=undef&&udata[l][k][j][i-1]!=undef&&vdata[l][k][j+1][i]!=undef&&vdata[l][k][j-1][i]!=undef)
						sdata[l][k][j][i]=
						(udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[ystart-1+j]*2)-
						(vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy*2)+
						vdata[l][k][j][i]*ltan[ystart-1+j]/EARTH_RADIUS;
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(udata[l][k][j][1]!=undef&&udata[l][k][j][0]!=undef&&vdata[l][k][j+1][0]!=undef&&vdata[l][k][j-1][0]!=undef)
						sdata[l][k][j][0]=
						(udata[l][k][j  ][1]-udata[l][k][j  ][0])/dxs[ystart-1+j]-
						(vdata[l][k][j+1][0]-vdata[l][k][j-1][0])/(dy*2)+
						vdata[l][k][j][0]*ltan[ystart-1+j]/EARTH_RADIUS;
					
					if(udata[l][k][j][x-1]!=undef&&udata[l][k][j][x-2]!=undef&&vdata[l][k][j+1][x-1]!=undef&&vdata[l][k][j-1][x-1]!=undef)
						sdata[l][k][j][x-1]=
						(udata[l][k][j  ][x-1]-udata[l][k][j  ][x-2])/dxs[ystart-1+j]-
						(vdata[l][k][j+1][x-1]-vdata[l][k][j-1][x-1])/(dy*2)+
						vdata[l][k][j][x-1]*ltan[ystart-1+j]/EARTH_RADIUS;
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(udata[l][k][0][i+1]!=undef&&udata[l][k][0][i-1]!=undef&&vdata[l][k][1][i]!=undef&&vdata[l][k][0][i]!=undef)
						sdata[l][k][0][i]=
						(udata[l][k][0][i+1]-udata[l][k][0][i-1])/(dxs[ystart-1]*2)-
						(vdata[l][k][1][i  ]-vdata[l][k][0][i  ])/dy+
						vdata[l][k][0][i]*ltan[ystart-1]/EARTH_RADIUS;
					
					if(udata[l][k][y-1][i+1]!=undef&&udata[l][k][y-1][i-1]!=undef&&vdata[l][k][y-1][i]!=undef&&vdata[l][k][y-2][i]!=undef)
						sdata[l][k][y-1][i]=
						(udata[l][k][y-1][i+1]-udata[l][k][y-1][i-1])/(dxs[ystart-2+y]*2)-
						(vdata[l][k][y-2][i  ]-vdata[l][k][y-1][i  ])/dy+
						vdata[l][k][y-1][i]*ltan[ystart-2+y]/EARTH_RADIUS;
				}
				
				/*** corner points ***/{
					if(udata[l][k][0][0]!=undef&&udata[l][k][0][1]!=undef&&vdata[l][k][0][0]!=undef&&vdata[l][k][1][0]!=undef)
						sdata[l][k][0][0]=
						(udata[l][k][0][1]-udata[l][k][0][0])/dxs[ystart-1]-
						(vdata[l][k][1][0]-vdata[l][k][0][0])/dy+
						vdata[l][k][0][0]*ltan[ystart-1]/EARTH_RADIUS;
					
					if(udata[l][k][y-1][0]!=undef&&udata[l][k][y-1][1]!=undef&&vdata[l][k][y-1][0]!=undef&&vdata[l][k][y-2][0]!=undef)
						sdata[l][k][y-1][0]=
						(udata[l][k][y-1][1]-udata[l][k][y-1][0])/dxs[ystart-2+y]-
						(vdata[l][k][y-1][0]-vdata[l][k][y-2][0])/dy+
						vdata[l][k][y-1][0]*ltan[ystart-2+y]/EARTH_RADIUS;
					
					if(udata[l][k][y-1][x-1]!=undef&&udata[l][k][y-1][x-2]!=undef&&vdata[l][k][y-1][x-1]!=undef&&vdata[l][k][y-2][x-1]!=undef)
						sdata[l][k][y-1][x-1]=
						(udata[l][k][y-1][x-1]-udata[l][k][y-1][x-2])/dxs[ystart-2+y]-
						(vdata[l][k][y-1][x-1]-vdata[l][k][y-2][x-1])/dy+
						vdata[l][k][y-1][x-1]*ltan[ystart-2+y]/EARTH_RADIUS;
					
					if(udata[l][k][0][x-1]!=undef&&udata[l][k][0][x-2]!=undef&&vdata[l][k][0][x-1]!=undef&&vdata[l][k][1][x-1]!=undef)
						sdata[l][k][0][x-1]=
						(udata[l][k][0][x-1]-udata[l][k][0][x-2])/dxs[ystart-1]-
						(vdata[l][k][1][x-1]-vdata[l][k][0][x-1])/dy+
						vdata[l][k][0][x-1]*ltan[ystart-1]/EARTH_RADIUS;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(udata[k][j][i+1][l]!=undef&&udata[k][j][i-1][l]!=undef&&vdata[k][j+1][i][l]!=undef&&vdata[k][j-1][i][l]!=undef)
						sdata[k][j][i][l]=
						(udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[ystart-1+j]*2)-
						(vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy*2)+
						(vdata[k][j][i][l]/EARTH_RADIUS)*ltan[ystart-1+j];
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(udata[k][j][1][l]!=undef&&udata[k][j][0][l]!=undef&&vdata[k][j+1][0][l]!=undef&&vdata[k][j-1][0][l]!=undef)
						sdata[k][j][0][l]=
						(udata[k][j  ][1][l]-udata[k][j  ][0][l])/dxs[ystart-1+j]-
						(vdata[k][j+1][0][l]-vdata[k][j-1][0][l])/(dy*2)+
						vdata[k][j][0][l]*ltan[ystart-1+j]/EARTH_RADIUS;
					
					if(udata[k][j][x-1][l]!=undef&&udata[k][j][x-2][l]!=undef&&vdata[k][j+1][x-1][l]!=undef&&vdata[k][j-1][x-1][l]!=undef)
						sdata[k][j][x-1][l]=
						(udata[k][j  ][x-1][l]-udata[k][j  ][x-2][l])/dxs[ystart-1+j]-
						(vdata[k][j+1][x-1][l]-vdata[k][j-1][x-1][l])/(dy*2)+
						vdata[k][j][x-1][l]*ltan[ystart-1+j]/EARTH_RADIUS;
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(udata[k][0][i+1][l]!=undef&&udata[k][0][i-1][l]!=undef&&vdata[k][1][i][l]!=undef&&vdata[k][0][i][l]!=undef)
						sdata[k][0][i][l]=
						(udata[k][0][i+1][l]-udata[k][0][i-1][l])/(dxs[ystart-1]*2)-
						(vdata[k][1][i  ][l]-vdata[k][0][i  ][l])/dy+
						vdata[k][0][i][l]*ltan[ystart-1]/EARTH_RADIUS;
					
					if(udata[k][y-1][i+1][l]!=undef&&udata[k][y-1][i-1][l]!=undef&&vdata[k][y-1][i][l]!=undef&&vdata[k][y-2][i][l]!=undef)
						sdata[k][y-1][i][l]=
						(udata[k][y-1][i+1][l]-udata[k][y-1][i-1][l])/(dxs[ystart-2+y]*2)-
						(vdata[k][y-2][i  ][l]-vdata[k][y-1][i  ][l])/dy+
						vdata[k][y-1][i][l]*ltan[ystart-2+y]/EARTH_RADIUS;
				}
				
				/*** corner points ***/
				if(udata[k][0][0][l]!=undef&&udata[k][0][1][l]!=undef&&vdata[k][0][0][l]!=undef&&vdata[k][1][0][l]!=undef)
					sdata[k][0][0][l]=
					(udata[k][0][1][l]-udata[k][0][0][l])/dxs[ystart-1]-
					(vdata[k][1][0][l]-vdata[k][0][0][l])/dy+
					vdata[k][0][0][l]*ltan[ystart-1]/EARTH_RADIUS;
				
				if(udata[k][y-1][0][l]!=undef&&udata[k][y-1][1][l]!=undef&&vdata[k][y-1][0][l]!=undef&&vdata[k][y-2][0][l]!=undef)
					sdata[k][y-1][0][l]=
					(udata[k][y-1][1][l]-udata[k][y-1][0][l])/dxs[ystart-2+y]-
					(vdata[k][y-1][0][l]-vdata[k][y-2][0][l])/dy+
					vdata[k][y-1][0][l]*ltan[ystart-2+y]/EARTH_RADIUS;
				
				if(udata[k][y-1][x-1][l]!=undef&&udata[k][y-1][x-2][l]!=undef&&vdata[k][y-1][x-1][l]!=undef&&vdata[k][y-2][x-1][l]!=undef)
					sdata[k][y-1][x-1][l]=
					(udata[k][y-1][x-1][l]-udata[k][y-1][x-2][l])/dxs[ystart-2+y]-
					(vdata[k][y-1][x-1][l]-vdata[k][y-2][x-1][l])/dy+
					vdata[k][y-1][x-1][l]*ltan[ystart-2+y]/EARTH_RADIUS;
				
				if(udata[k][0][x-1][l]!=undef&&udata[k][0][x-2][l]!=undef&&vdata[k][0][x-1][l]!=undef&&vdata[k][1][x-1][l]!=undef)
					sdata[k][0][x-1][l]=
					(udata[k][0][x-1][l]-udata[k][0][x-2][l])/dxs[ystart-1]-
					(vdata[k][1][x-1][l]-vdata[k][0][x-1][l])/dy+
					vdata[k][0][x-1][l]*ltan[ystart-1]/EARTH_RADIUS;
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
					(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
					(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(dy*2)-
					udata[l][k][j][i]*ltan[ystart-1+j]/EARTH_RADIUS;
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(vdata[l][k][j][1]!=undef&&vdata[l][k][j][0]!=undef&&udata[l][k][j+1][0]!=undef&&udata[l][k][j-1][0]!=undef)
						sdata[l][k][j][0]=
						(vdata[l][k][j][1]-vdata[l][k][j][0])/dxs[ystart-1+j]+
						(udata[l][k][j+1][0]-udata[l][k][j-1][0])/(dy*2)-
						udata[l][k][j][0]*ltan[ystart-1+j]/EARTH_RADIUS;
					
					if(vdata[l][k][j][x-1]!=undef&&vdata[l][k][j][x-2]!=undef&&udata[l][k][j+1][x-1]!=undef&&udata[l][k][j-1][x-1]!=undef)
						sdata[l][k][j][x-1]=
						(vdata[l][k][j  ][x-1]-vdata[l][k][j  ][x-2])/dxs[ystart-1+j]+
						(udata[l][k][j+1][x-1]-udata[l][k][j-1][x-1])/(dy*2)-
						udata[l][k][j][x-1]*ltan[ystart-1+j]/EARTH_RADIUS;
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(vdata[l][k][0][i+1]!=undef&&vdata[l][k][0][i-1]!=undef&&udata[l][k][1][i]!=undef&&udata[l][k][0][i]!=undef)
						sdata[l][k][0][i]=
						(vdata[l][k][0][i+1]-vdata[l][k][0][i-1])/(dxs[ystart-1]*2)+
						(udata[l][k][1][i  ]-udata[l][k][0][i  ])/dy-
						udata[l][k][0][i]*ltan[ystart-1]/EARTH_RADIUS;
					
					if(vdata[l][k][y-1][i+1]!=undef&&vdata[l][k][y-1][i-1]!=undef&&udata[l][k][y-1][i]!=undef&&udata[l][k][y-2][i]!=undef)
						sdata[l][k][y-1][i]=
						(vdata[l][k][y-1][i+1]-vdata[l][k][y-1][i-1])/(dxs[ystart-2+y]*2)+
						(udata[l][k][y-1][i]-udata[l][k][y-2][i])/dy-
						udata[l][k][y-1][i]*ltan[ystart-2+y]/EARTH_RADIUS;
				}
				
				/*** corner points ***/
				if(vdata[l][k][0][1]!=undef&&vdata[l][k][0][0]!=undef&&udata[l][k][1][0]!=undef&&udata[l][k][0][0]!=undef)
					sdata[l][k][0][0]=
					(vdata[l][k][0][1]-vdata[l][k][0][0])/dxs[ystart-1]+
					(udata[l][k][1][0]-udata[l][k][0][0])/dy-
					udata[l][k][0][0]*ltan[ystart-1]/EARTH_RADIUS;
				
				if(vdata[l][k][0][x-1]!=undef&&vdata[l][k][0][x-2]!=undef&&udata[l][k][1][x-1]!=undef&&udata[l][k][0][x-1]!=undef)
					sdata[l][k][0][x-1]=
					(vdata[l][k][0][x-1]-vdata[l][k][0][x-2])/dxs[ystart-1]+
					(udata[l][k][1][x-1]-udata[l][k][0][x-1])/dy-
					udata[l][k][0][x-1]*ltan[ystart-1]/EARTH_RADIUS;
				
				if(vdata[l][k][y-1][x-1]!=undef&&vdata[l][k][y-1][x-2]!=undef&&udata[l][k][y-1][x-1]!=undef&&udata[l][k][y-2][x-1]!=undef)
					sdata[l][k][y-1][x-1]=
					(vdata[l][k][y-1][x-1]-vdata[l][k][y-1][x-2])/dxs[ystart-2+y]+
					(udata[l][k][y-1][x-1]-udata[l][k][y-2][x-1])/dy-
					udata[l][k][y-1][x-1]*ltan[ystart-2+y]/EARTH_RADIUS;
				
				if(vdata[l][k][y-1][1]!=undef&&vdata[l][k][y-1][0]!=undef&&udata[l][k][y-1][0]!=undef&&udata[l][k][y-2][0]!=undef)
					sdata[l][k][y-1][0]=
					(vdata[l][k][y-1][1]-vdata[l][k][y-1][0])/dxs[ystart-2+y]+
					(udata[l][k][y-1][0]-udata[l][k][y-2][0])/dy-
					udata[l][k][y-1][0]*ltan[ystart-2+y]/EARTH_RADIUS;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef&&udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef)
						sdata[k][j][i][l]=
						(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
						(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(dy*2)-
						udata[k][j][i][l]*ltan[ystart-1+j]/EARTH_RADIUS;
				}
				
				/*** left and right boundry ***/
				for(int j=1;j<y-1;j++){
					if(vdata[k][j][1][l]!=undef&&vdata[k][j][0][l]!=undef&&udata[k][j+1][0][l]!=undef&&udata[k][j-1][0][l]!=undef)
						sdata[k][j][0][l]=
						(vdata[k][j  ][1][l]-vdata[k][j  ][0][l])/dxs[ystart-1+j]+
						(udata[k][j+1][0][l]-udata[k][j-1][0][l])/(dy*2)-
						udata[k][j][0][l]*ltan[ystart-1+j]/EARTH_RADIUS;
					
					if(vdata[k][j][x-1][l]!=undef&&vdata[k][j][x-2][l]!=undef&&udata[k][j+1][x-1][l]!=undef&&udata[k][j-1][x-1][l]!=undef)
						sdata[k][j][x-1][l]=
						(vdata[k][j  ][x-1][l]-vdata[k][j  ][x-2][l])/dxs[ystart-1+j]+
						(udata[k][j+1][x-1][l]-udata[k][j-1][x-1][l])/(dy*2)-
						udata[k][j][x-1][l]*ltan[ystart-1+j]/EARTH_RADIUS;
				}
				
				/*** top and bottom boundry ***/
				for(int i=1;i<x-1;i++){
					if(vdata[k][0][i+1][l]!=undef&&vdata[k][0][i-1][l]!=undef&&udata[k][1][i][l]!=undef&&udata[k][0][i][l]!=undef)
						sdata[k][0][i][l]=
						(vdata[k][0][i+1][l]-vdata[k][0][i-1][l])/(dxs[ystart-1]*2)+
						(udata[k][1][i  ][l]-udata[k][0][i  ][l])/dy-
						udata[k][0][i][l]*ltan[ystart-1]/EARTH_RADIUS;
					
					if(vdata[k][y-1][i+1][l]!=undef&&vdata[k][y-1][i-1][l]!=undef&&udata[k][y-1][i][l]!=undef&&udata[k][y-2][i][l]!=undef)
						sdata[k][y-1][i][l]=
						(vdata[k][y-1][i+1][l]-vdata[k][y-1][i-1][l])/(dxs[ystart-2+y]*2)+
						(udata[k][y-1][i][l]-udata[k][y-2][i][l])/dy-
						udata[k][y-1][i][l]*ltan[ystart-2+y]/EARTH_RADIUS;
				}
				
				/*** corner points ***/
				if(vdata[k][0][1][l]!=undef&&vdata[k][0][0][l]!=undef&&udata[k][1][0][l]!=undef&&udata[k][0][0][l]!=undef)
					sdata[k][0][0][l]=
					(vdata[k][0][1][l]-vdata[k][0][0][l])/dxs[ystart-1]+
					(udata[k][1][0][l]-udata[k][0][0][l])/dy-
					udata[k][0][0][l]*ltan[ystart-1]/EARTH_RADIUS;
				
				if(vdata[k][0][x-1][l]!=undef&&vdata[k][0][x-2][l]!=undef&&udata[k][1][x-1][l]!=undef&&udata[k][0][x-1][l]!=undef)
					sdata[k][0][x-1][l]=
					(vdata[k][0][x-1][l]-vdata[k][0][x-2][l])/dxs[ystart-1]+
					(udata[k][1][x-1][l]-udata[k][0][x-1][l])/dy-
					udata[k][0][x-1][l]*ltan[ystart-1]/EARTH_RADIUS;
				
				if(vdata[k][y-1][x-1][l]!=undef&&vdata[k][y-1][x-2][l]!=undef&&udata[k][y-1][x-1][l]!=undef&&udata[k][y-2][x-1][l]!=undef)
					sdata[k][y-1][x-1][l]=
					(vdata[k][y-1][x-1][l]-vdata[k][y-1][x-2][l])/dxs[ystart-2+y]+
					(udata[k][y-1][x-1][l]-udata[k][y-2][x-1][l])/dy-
					udata[k][y-1][x-1][l]*ltan[ystart-2+y]/EARTH_RADIUS;
				
				if(vdata[k][y-1][1][l]!=undef&&vdata[k][y-1][0][l]!=undef&&udata[k][y-1][0][l]!=undef&&udata[k][y-2][0][l]!=undef)
					sdata[k][y-1][0][l]=
					(vdata[k][y-1][1][l]-vdata[k][y-1][0][l])/dxs[ystart-2+y]+
					(udata[k][y-1][0][l]-udata[k][y-2][0][l])/dy-
					udata[k][y-1][0][l]*ltan[ystart-2+y]/EARTH_RADIUS;
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
				vsdata[l][k][j][i]=HORCON_PARAM*dxs[j]*dy/2*(float)Math.sqrt(
					Math.pow((udata[l][k][j][i+1]-udata[l][k][j][i-1])/(dxs[j]+dxs[j]),2)+
					Math.pow((vdata[l][k][j+1][i]-vdata[l][k][j-1][i])/(dy   +dy   ),2)+
					Math.pow((vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(dxs[j]+dxs[j])+
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
				vsdata[k][j][i][l]=HORCON_PARAM*dxs[j]*dy/2*(float)Math.sqrt(
					Math.pow((udata[k][j][i+1][l]-udata[k][j][i-1][l])/(dxs[j]+dxs[j]),2)+
					Math.pow((vdata[k][j+1][i][l]-vdata[k][j-1][i][l])/(dy   +dy   ),2)+
					Math.pow((vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(dxs[j]+dxs[j])+
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
		pv.setCommentAndUnit("potential vorticity");
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
				pvdata[l][k][j][i]=-GRAVITY_ACCERLERATION*vordata[l][k][j][i]*
				(Tedata[l][k+1][j][i]-Tedata[l][k-1][j][i])/(dz+dz);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vordata[k][j][i][l]!=undef&&Tedata[k+1][j][i][l]!=undef&&Tedata[k-1][j][i][l]!=undef){
				pvdata[k][j][i][l]=-GRAVITY_ACCERLERATION*vordata[k][j][i][l]*
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
			for(int k=0;k<z;k++){
				/*** internal area ***/
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef&&adata[l][k][j][i+1]!=undef
					 &&adata[l][k][j][i-1]!=undef&&adata[l][k][j+1][i]!=undef&&adata[l][k][j-1][i]!=undef)
						atdata[l][k][j][i]=
						udata[l][k][j][i]*(adata[l][k][j][i+1]-adata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
						vdata[l][k][j][i]*(adata[l][k][j+1][i]-adata[l][k][j-1][i])/(dy*2);
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** internal area ***/
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++){
					if(udata[k][j][i][l]!=undef&&adata[k][j][i+1][l]!=undef&&adata[k][j][i-1][l]!=undef
					 &&vdata[k][j][i][l]!=undef&&adata[k][j+1][i][l]!=undef&&adata[k][j-1][i][l]!=undef)
						atdata[k][j][i][l]=
						udata[k][j][i][l]*(adata[k][j][i+1][l]-adata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
						vdata[k][j][i][l]*(adata[k][j+1][i][l]-adata[k][j-1][i][l])/(dy*2);
				}
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
		
		Variable cov=new Variable(v.getName()+"con",v);
		cov.setCommentAndUnit("convection term");
		
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
		
		Variable F=new Variable("lp",v);	F.setCommentAndUnit("Laplacian of "+v.getName());
		
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
					)/dx/dxs[ystart-1+j]+(
						(gdata[l][k][j+1][0]-gdata[l][k][j][0])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
						(gdata[l][k][j][0]-gdata[l][k][j-1][0])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
					)/dy/dy;
					
					/*** east BC with i=x-1 ***/
					Fdata[l][k][j][x-1]=(
						(gdata[l][k][j][0  ]-gdata[l][k][j][x-1])-
						(gdata[l][k][j][x-1]-gdata[l][k][j][x-2])
					)/dx/dxs[ystart-1+j]+(
						(gdata[l][k][j+1][x-1]-gdata[l][k][j][x-1])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
						(gdata[l][k][j][x-1]-gdata[l][k][j-1][x-1])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
					)/dy/dy;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=(
					(gdata[l][k][j][i+1]-gdata[l][k][j][i])-
					(gdata[l][k][j][i]-gdata[l][k][j][i-1])
				)/dx/dxs[ystart-1+j]+(
					(gdata[l][k][j+1][i]-gdata[l][k][j][i])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
					(gdata[l][k][j][i]-gdata[l][k][j-1][i])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
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
					)/dx/dxs[ystart-1+j]+(
						(gdata[k][j+1][0][l]-gdata[k][j][0][l])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
						(gdata[k][j][0][l]-gdata[k][j-1][0][l])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
					)/dy/dy;
					
					/*** east BC with i=x-1 ***/
					Fdata[k][j][x-1][l]=(
						(gdata[k][j][0  ][l]-gdata[k][j][x-1][l])-
						(gdata[k][j][x-1][l]-gdata[k][j][x-2][l])
					)/dx/dxs[ystart-1+j]+(
						(gdata[k][j+1][x-1][l]-gdata[k][j][x-1][l])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
						(gdata[k][j][x-1][l]-gdata[k][j-1][x-1][l])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
					)/dy/dy;
				}
				
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=(
					(gdata[k][j][i+1][l]-gdata[k][j][i][l])-
					(gdata[k][j][i][l]-gdata[k][j][i-1][l])
				)/dx/dxs[ystart-1+j]+(
					(gdata[k][j+1][i][l]-gdata[k][j][i][l])*(float)Math.cos((ydef[ystart  +j]+ydef[ystart-1+j])/2f)-
					(gdata[k][j][i][l]-gdata[k][j-1][i][l])*(float)Math.cos((ydef[ystart-1+j]+ydef[ystart-2+j])/2f)
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
	public Variable cConcentration(Variable count){
		assignSubDomainParams(count);
		
		Variable conc=new Variable("conc",count);
		conc.setCommentAndUnit("concentration");
		
		float[][][][] ctdata=count.getData();
		float[][][][] ccdata=conc.getData();
		
	    if(conc.isTFirst()){
			for(int j=0;j<y;j++){
				float deltaY=(j==0||j==y-1)?dy/2f:dy;
				
				for(int i=0;i<x;i++){
					float deltaX=dxs[j];
					
					if(!sm.isPeriodicX()&&(i==0||i==x-1)) deltaX/=2f;
					
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
					float deltaX=dxs[j];
					
					if(!sm.isPeriodicX()&&(i==0||i==x-1)) deltaX/=2f;
					
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
		
		Variable h=new Variable("hel",u);	h.setCommentAndUnit("helicity ()");
		
		float[][][][] hdata=h.getData();
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		float[][][][] wdata=w.getData();
		
		if(u.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				if(
					udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef&&wdata[l][k][j][i]!=undef&&
					udata[l][k][j+1][i]!=undef&&udata[l][k][j-1][i]!=undef&&udata[l][k+1][j][i]!=undef&&udata[l][k-1][j][i]!=undef&&
					vdata[l][k][j][i+1]!=undef&&vdata[l][k][j][i-1]!=undef&&vdata[l][k+1][j][i]!=undef&&vdata[l][k-1][j][i]!=undef&&
					wdata[l][k][j][i+1]!=undef&&wdata[l][k][j][i-1]!=undef&&wdata[l][k][j+1][i]!=undef&&wdata[l][k][j-1][i]!=undef
				)
					hdata[l][k][j][i]=
						udata[l][k][j][i]*(
							(wdata[l][k][j+1][i]-wdata[l][k][j-1][i])/(2*dy)-
							(vdata[l][k+1][j][i]-vdata[l][k-1][j][i])/(dz+dz)
						)+
						vdata[l][k][j][i]*(
							(udata[l][k+1][j][i]-udata[l][k-1][j][i])/(dz+dz)-
							(wdata[l][k][j][i+1]-wdata[l][k][j][i-1])/(2*dxs[ystart-1+j])
						)+
						wdata[l][k][j][i]*(
							(vdata[l][k][j][i+1]-vdata[l][k][j][i-1])/(2*dxs[ystart-1+j])-
							(udata[l][k][j+1][i]-udata[l][k][j-1][i])/(2*dy)
						);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				if(
					udata[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef&&wdata[k][j][i][l]!=undef&&
					udata[k][j+1][i][l]!=undef&&udata[k][j-1][i][l]!=undef&&udata[k+1][j][i][l]!=undef&&udata[k-1][j][i][l]!=undef&&
					vdata[k][j][i+1][l]!=undef&&vdata[k][j][i-1][l]!=undef&&vdata[k+1][j][i][l]!=undef&&vdata[k-1][j][i][l]!=undef&&
					wdata[k][j][i+1][l]!=undef&&wdata[k][j][i-1][l]!=undef&&wdata[k][j+1][i][l]!=undef&&wdata[k][j-1][i][l]!=undef
				)
					hdata[k][j][i][l]=
						udata[k][j][i][l]*(
							(wdata[k][j+1][i][l]-wdata[k][j-1][i][l])/(2*dy)-
							(vdata[k+1][j][i][l]-vdata[k-1][j][i][l])/(dz+dz)
						)+
						vdata[k][j][i][l]*(
							(udata[k+1][j][i][l]-udata[k-1][j][i][l])/(dz+dz)-
							(wdata[k][j][i+1][l]-wdata[k][j][i-1][l])/(2*dxs[ystart-1+j])
						)+
						wdata[k][j][i][l]*(
							(vdata[k][j][i+1][l]-vdata[k][j][i-1][l])/(2*dxs[ystart-1+j])-
							(udata[k][j+1][i][l]-udata[k][j-1][i][l])/(2*dy)
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
		
		Variable w=new Variable("omega",div);	w.setCommentAndUnit("vertical velocity (Pa s^-1)");
		
		float[][][][] wdata  =w.getData();
		float[][][][] divdata=div.getData();
		
		if(div.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(divdata[l][k-1][j][i]!=undef&&divdata[l][k][j][i]!=undef)
					if(wdata[l][k-1][j][i]!=undef)
						wdata[l][k][j][i]=wdata[l][k-1][j][i]+
						(divdata[l][k][j][i]+divdata[l][k-1][j][i])/2*dz;
					else
						wdata[l][k][j][i]=
						(divdata[l][k][j][i]+divdata[l][k-1][j][i])/2*dz;
						
				else wdata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=1;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(divdata[i][j][k-1][l]!=undef&&divdata[k][j][i][l]!=undef)
					if(wdata[k-1][j][i][l]!=undef)
						wdata[k][j][i][l]=wdata[k-1][j][i][l]+
						(divdata[k][j][i][l]+divdata[k-1][j][i][l])/2*dz;
					else
						wdata[k][j][i][l]=
						(divdata[k][j][i][l]+divdata[k-1][j][i][l])/2*dz;
				
				else wdata[k][j][i][l]=undef;
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
     * Calculate Eliassen-Palm flux vector.
     *
     * @param	ua		u-wind anomaly (m s^-1)
     * @param	va		v-wind anomaly (m s^-1)
     * @param	tha		potential temperature anomaly (K)
     * @param	thm		potential temperature zonal mean (K)
     *
     * @return	re		re[0] is meridional component and re[1] is vertical component
     */
	public Variable[] cEPFlux(Variable ua,Variable va,Variable tha,Variable thm){
		if(thm.getXCount()!=1) throw new IllegalArgumentException("invalid x-count");
		
		checkDimensions(ua,va,tha);
		assignSubDomainParams(ua);
		
		Variable[] flx=new Variable[2];
		flx[0]=new Variable("EPfy",thm);
		flx[1]=new Variable("EPfp",thm);
		flx[0].setUndef(undef);	flx[0].setCommentAndUnit("EP-Flux meridional component");
		flx[1].setUndef(undef);	flx[1].setCommentAndUnit("EP-Flux vertical component");
		
		float[][][][]  udata= ua.getData();
		float[][][][]  vdata= va.getData();
		float[][][][]  adata=tha.getData();
		float[][][][]  mdata=thm.getData();
		float[][][][] f0data=flx[0].getData();
		float[][][][] f1data=flx[1].getData();
		
		float[][] tmp=new float[z][y];
		
		if(ua.isTFirst()){
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(udata[l][k][j][i]!=undef&&vdata[l][k][j][i]!=undef){
						f0data[l][k][j][0]+=udata[l][k][j][i]*vdata[l][k][j][i];
						count++;
					}
					
					if(count!=0) f0data[l][k][j][0]/=-count;
					
					count=0;
					
					for(int i=0;i<x;i++)
					if(vdata[l][k][j][i]!=undef&&adata[l][k][j][i]!=undef){
						tmp[k][j]+=vdata[l][k][j][i]*adata[l][k][j][i];
						count++;
					}
					
					if(count!=0) tmp[k][j]/=count/f1[ystart-1+j];
				}
				
				for(int j=0;j<y;j++){
					f1data[l][0][j][0]=
					tmp[0][j]*dz/(mdata[l][1][j][0]-mdata[l][0][j][0]);
					
					for(int k=1;k<z-1;k++){
						f1data[l][k][j][0]=
						tmp[k][j]*(dz+dz)/(mdata[l][k+1][j][0]-mdata[l][k-1][j][0]);
					}
					
					f1data[l][z-1][j][0]=
					tmp[z-1][j]*dz/(mdata[l][z-1][j][0]-mdata[l][z-2][j][0]);
				}
				
				for(int j=0;j<y;j++)
				for(int k=0;k<z;k++) tmp[k][j]=0;
			}
			
		}else{
			for(int l=0;l<t;l++){
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					int count=0;
					
					for(int i=0;i<x;i++)
					if(udata[k][j][i][l]!=undef&&vdata[k][j][i][l]!=undef){
						f0data[k][j][0][l]+=udata[k][j][i][l]*vdata[k][j][i][l];
						count++;
					}
					
					if(count!=0) f0data[k][j][0][l]/=-count;
					
					count=0;
					
					for(int i=0;i<x;i++)
					if(vdata[k][j][i][l]!=undef&&adata[k][j][i][l]!=undef){
						tmp[k][j]+=vdata[k][j][i][l]*adata[k][j][i][l];
						count++;
					}
					
					if(count!=0) tmp[k][j]/=count/f1[ystart-1+j];
				}
				
				for(int j=0;j<y;j++){
					f1data[0][j][0][l]=
					tmp[0][j]*dz/(mdata[1][j][0][l]-mdata[0][j][0][l]);
					
					for(int k=1;k<z-1;k++){
						f1data[k][j][0][l]=
						tmp[k][j]*(dz+dz)/(mdata[k+1][j][0][l]-mdata[k-1][j][0][l]);
					}
					
					f1data[z-1][j][0][l]=
					tmp[z-1][j]*dz/(mdata[z-1][j][0][l]-mdata[z-2][j][0][l]);
				}
				
				for(int j=0;j<y;j++)
				for(int k=0;k<z;k++) tmp[k][j]=0;
			}
		}
		
		return flx;
	}
	
	
	/**
     * Calculate 2D wave activity flux.
     * Reference: Plumb 1985 JAS
     *
     * @param	ua	u-wind anomaly (m s^-1)
     * @param	va	v-wind anomaly (m s^-1)
     * @param	ha	geopotential height anomaly (m)
     *
     * @return	re	re[0] is zonal component and re[1] is meridional component
     */
	public Variable[] cWaveActivityFlux(Variable ua,Variable va,Variable ha){
		checkDimensions(ua,va,ha);
		assignSubDomainParams(ua);
		
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
		
		Variable taux=new Variable("taux",uwnd); taux.setUndef(undef);
		Variable tauy=new Variable("tauy",vwnd); tauy.setUndef(undef);
		
		taux.setCommentAndUnit("zonal wind stress (N)");
		tauy.setCommentAndUnit("meridianal wind stress (N)");
		
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
