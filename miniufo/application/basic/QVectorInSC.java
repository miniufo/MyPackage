/**
 * @(#)QVectorInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.application.EquationInSphericalCoordinate;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;


/**
 * QVector in spheral coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class QVectorInSC extends EquationInSphericalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public QVectorInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/**
     * calculate quasi-geostrophic QVector, [0] is x direction and [1] is y direction
     *
     * @param	Ug		geostropic-wind in x direction
     * @param	Vg		geostropic-wind in y direction
     * @param	air		temperature(K)
     *
     * @return	quasi-geostrophic QVector
     *
     * @exception	if Ug,Vg,air are not dimensionlly the same
     */
	public Variable[] cQuasiGeostrophicQVector(Variable Ug,Variable Vg,Variable air){
		if(!Ug.isLike(Vg)||!Ug.isLike(air)) throw new IllegalArgumentException("dimensions not same");
		
		ystart=Ug.getRange().getYRange()[0];
		
		t=air.getTCount();	z=air.getZCount();	y=air.getYCount();	x=air.getXCount();
		
		Variable[] qv=new Variable[2];
		qv[0]=new Variable("qvx",Ug);	qv[0].setCommentAndUnit("quasi-geostrophic QVector in x-direction");
		qv[1]=new Variable("qvy",Ug);	qv[1].setCommentAndUnit("quasi-geostrophic QVector in y-direction");

		float[][][][] airdata=air.getData();
		float[][][][]  Ugdata= Ug.getData();	float[][][][]  Vgdata= Vg.getData();
		float[][][][] qvxdata=qv[0].getData();	float[][][][] qvydata=qv[1].getData();
		
		float undef=Ug.getUndef();
		
		if(Ug.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				if(Ugdata[l][k][j][i+1]!=undef&& Ugdata[l][k][j][i-1]!=undef&&
				  airdata[l][k][j][i+1]!=undef&&airdata[l][k][j][i-1]!=undef&&
				   Vgdata[l][k][j][i+1]!=undef&& Vgdata[l][k][j][i-1]!=undef&&
				  airdata[l][k][j+1][i]!=undef&&airdata[l][k][j-1][i]!=undef){
					
					qvxdata[l][k][j][i]=-Rd/zdef[k]*(
						( Ugdata[l][k][j][i+1]- Ugdata[l][k][j][i-1])/dxs[ystart-1+j]*
						(airdata[l][k][j][i+1]-airdata[l][k][j][i-1])/dxs[ystart-1+j]+
						( Vgdata[l][k][j][i+1]- Vgdata[l][k][j][i-1])/dxs[ystart-1+j]*
						(airdata[l][k][j+1][i]-airdata[l][k][j-1][i])/dy
					)/4;
					
				}else qvxdata[l][k][j][i]=undef;
				
				if(Ugdata[l][k][j+1][i]!=undef&& Ugdata[l][k][j-1][i]!=undef&&
				  airdata[l][k][j][i+1]!=undef&&airdata[l][k][j][i-1]!=undef&&
				   Vgdata[l][k][j+1][i]!=undef&& Vgdata[l][k][j-1][i]!=undef&&
				  airdata[l][k][j+1][i]!=undef&&airdata[l][k][j-1][i]!=undef){
					
					qvydata[l][k][j][i]=-Rd/zdef[k]*(
						( Ugdata[l][k][j+1][i]- Ugdata[l][k][j-1][i])/dy*
						(airdata[l][k][j][i+1]-airdata[l][k][j][i-1])/dxs[ystart-1+j]+
						( Vgdata[l][k][j+1][i]- Vgdata[l][k][j-1][i])/dy*
						(airdata[l][k][j+1][i]-airdata[l][k][j-1][i])/dy
					)/4;
					
				}else qvydata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				if(Ugdata[l][k][j][i+1]!=undef&& Ugdata[l][k][j][i-1]!=undef&&
				  airdata[l][k][j][i+1]!=undef&&airdata[l][k][j][i-1]!=undef&&
				   Vgdata[l][k][j][i+1]!=undef&& Vgdata[l][k][j][i-1]!=undef&&
				  airdata[l][k][j+1][i]!=undef&&airdata[l][k][j-1][i]!=undef){
					
					qvxdata[k][j][i][l]=-Rd/zdef[k]*(
						( Ugdata[k][j][i+1][l]- Ugdata[k][j][i-1][l])/dxs[ystart-1+j]*
						(airdata[k][j][i+1][l]-airdata[k][j][i-1][l])/dxs[ystart-1+j]+
						( Vgdata[k][j][i+1][l]- Vgdata[k][j][i-1][l])/dxs[ystart-1+j]*
						(airdata[k][j+1][i][l]-airdata[k][j-1][i][l])/dy
					)/4;
					
				}else qvxdata[k][j][i][l]=undef;
				
				if(Ugdata[k][j+1][i][l]!=undef&& Ugdata[k][j-1][i][l]!=undef&&
				  airdata[k][j][i+1][l]!=undef&&airdata[k][j][i-1][l]!=undef&&
				   Vgdata[k][j+1][i][l]!=undef&& Vgdata[k][j-1][i][l]!=undef&&
				  airdata[k][j+1][i][l]!=undef&&airdata[k][j-1][i][l]!=undef){
					
					qvydata[k][j][i][l]=-Rd/zdef[k]*(
						( Ugdata[k][j+1][i][l]- Ugdata[k][j-1][i][l])/dy*
						(airdata[k][j][i+1][l]-airdata[k][j][i-1][l])/dxs[ystart+j]+
						( Vgdata[k][j+1][i][l]- Vgdata[k][j-1][i][l])/dy*
						(airdata[k][j+1][i][l]-airdata[k][j-1][i][l])/dy
					)/4;
					
				}else qvydata[k][j][i][l]=undef;
			}
		}
		
		return qv;
	}
	
	
	/**
     * calculate quasi-geostrophic QVector, [0] is x direction and [1] is y direction
     *
     * @param	hgt		geopotential height
     * @param	air		temperature(K)
     *
     * @return	quasi-geostrophic QVector
     *
     * @exception	if Ug,Vg,air are not dimensionlly the same
     */
	public Variable[] cQuasiGeostrophicQVector(Variable hgt,Variable air){
		if(!hgt.isLike(air)) throw new IllegalArgumentException("dimensions not same");
		
		ystart=hgt.getRange().getYRange()[0];
		
		t=air.getTCount();	z=air.getZCount();	y=air.getYCount();	x=air.getXCount();
		
		Variable[] qv=new Variable[2];
		qv[0]=new Variable("qvx",hgt);	qv[0].setCommentAndUnit("quasi-geostrophic QVector in x-direction");
		qv[1]=new Variable("qvy",hgt);	qv[1].setCommentAndUnit("quasi-geostrophic QVector in y-direction");
		
		float[][][][] hgtdata=hgt.getData();	float[][][][] airdata=air.getData();
		float[][][][] qvxdata=qv[0].getData();	float[][][][] qvydata=qv[1].getData();
		
		float undef=hgt.getUndef();
		
		if(hgt.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=2;i<x-2;i++){
				if(hgtdata[l][k][j  ][i  ]!=undef&&hgtdata[l][k][j+1][i-1]!=undef&&
				   hgtdata[l][k][j-1][i+1]!=undef&&hgtdata[l][k][j+1][i+1]!=undef&&
				   hgtdata[l][k][j-1][i-1]!=undef&&hgtdata[l][k][j  ][i+2]!=undef&&
				   hgtdata[l][k][j  ][i-1]!=undef&&airdata[l][k][j  ][i+1]!=undef&&
				   airdata[l][k][j  ][i-1]!=undef&&airdata[l][k][j+1][i  ]!=undef&&
				   airdata[l][k][j-1][i  ]!=undef){
					
					qvxdata[l][k][j][i]=-Rd/zdef[k]*(
						(hgtdata[l][k][j+1][i-1]+hgtdata[l][k][j-1][i+1]-hgtdata[l][k][j-1][i-1]-hgtdata[l][k][j+1][i+1])*
						(airdata[l][k][j  ][i+1]-airdata[l][k][j  ][i-1])/dxs[ystart-1+j]+
						(hgtdata[l][k][j  ][i+2]+hgtdata[l][k][j  ][i-2]-hgtdata[l][k][j][i]*2)*
						(airdata[l][k][j+1][i  ]-airdata[l][k][j-1][i  ])/dy
					)/dxs[ystart-1+j]/dxs[ystart-1+j]/dy/f1[j]/8;
					
				}else qvxdata[l][k][j][i]=undef;
				
				if(hgtdata[l][k][j  ][i  ]!=undef&&hgtdata[l][k][j+2][i  ]!=undef&&
				   hgtdata[l][k][j-2][i  ]!=undef&&hgtdata[l][k][j+1][i+1]!=undef&&
				   hgtdata[l][k][j+1][i-1]!=undef&&hgtdata[l][k][j-1][i+1]!=undef&&
				   hgtdata[l][k][j-1][i-1]!=undef&&airdata[l][k][j  ][i+1]!=undef&&
				   airdata[l][k][j  ][i-1]!=undef&&airdata[l][k][j+1][i  ]!=undef&&
				   airdata[l][k][j-1][i  ]!=undef){
					
					qvydata[l][k][j][i]=-Rd/zdef[k]*(
						((hgtdata[l][k][j  ][i  ]-hgtdata[l][k][j+2][i  ])/f1[j+1]+(hgtdata[l][k][j  ][i  ]-hgtdata[l][k][j-2][i  ])/f1[j-1])*
						( airdata[l][k][j  ][i+1]-airdata[l][k][j  ][i-1])+
						((hgtdata[l][k][j+1][i+1]-hgtdata[l][k][j+1][i-1])/f1[j+1]+(hgtdata[l][k][j-1][i+1]-hgtdata[l][k][j-1][i-1])/f1[j-1])*
						( airdata[l][k][j+1][i  ]-airdata[l][k][j-1][i  ])
					)/dy/dy/dxs[ystart-1+j]/8;
					
				}else qvydata[l][k][j][i]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=2;i<x-2;i++){
				if(hgtdata[k][j  ][i  ][l]!=undef&&hgtdata[k][j+1][i-1][l]!=undef&&
				   hgtdata[k][j-1][i+1][l]!=undef&&hgtdata[k][j+1][i+1][l]!=undef&&
				   hgtdata[k][j-1][i-1][l]!=undef&&hgtdata[k][j  ][i+2][l]!=undef&&
				   hgtdata[k][j  ][i-1][l]!=undef&&airdata[k][j  ][i+1][l]!=undef&&
				   airdata[k][j  ][i-1][l]!=undef&&airdata[k][j+1][i  ][l]!=undef&&
				   airdata[k][j-1][i  ][l]!=undef){
					
					qvxdata[k][j][i][l]=-Rd/zdef[k]*(
						(hgtdata[k][j+1][i-1][l]+hgtdata[k][j-1][i+1][l]-hgtdata[k][j-1][i-1][l]-hgtdata[k][j+1][i+1][l])*
						(airdata[k][j  ][i+1][l]-airdata[k][j  ][i-1][l])/dxs[ystart-1+j]+
						(hgtdata[k][j  ][i+2][l]+hgtdata[k][j  ][i-2][l]-hgtdata[k][j][i][l]*2)*
						(airdata[k][j+1][i  ][l]-airdata[k][j-1][i  ][l])/dy
					)/dxs[ystart-1+j]/dxs[ystart-1+j]/dy/f1[j]/8;
					
				}else qvxdata[k][j][i][l]=undef;
				
				if(hgtdata[k][j  ][i  ][l]!=undef&&hgtdata[k][j+2][i  ][l]!=undef&&
				   hgtdata[k][j-2][i  ][l]!=undef&&hgtdata[k][j+1][i+1][l]!=undef&&
				   hgtdata[k][j+1][i-1][l]!=undef&&hgtdata[k][j-1][i+1][l]!=undef&&
				   hgtdata[k][j-1][i-1][l]!=undef&&airdata[k][j  ][i+1][l]!=undef&&
				   airdata[k][j  ][i-1][l]!=undef&&airdata[k][j+1][i  ][l]!=undef&&
				   airdata[k][j-1][i  ][l]!=undef){
					
					qvydata[k][j][i][l]=-Rd/zdef[k]*(
						((hgtdata[k][j  ][i  ][l]-hgtdata[k][j+2][i  ][l])/f1[j+1]+(hgtdata[k][j  ][i  ][l]-hgtdata[k][j-2][i  ][l])/f1[j-1])*
						( airdata[k][j  ][i+1][l]-airdata[k][j  ][i-1][l])+
						((hgtdata[k][j+1][i+1][l]-hgtdata[k][j+1][i-1][l])/f1[j+1]+(hgtdata[k][j-1][i+1][l]-hgtdata[k][j-1][i-1][l])/f1[j-1])*
						( airdata[k][j+1][i  ][l]-airdata[k][j-1][i  ][l])
					)/dy/dy/dxs[ystart-1+j]/8;
					
				}else qvydata[k][j][i][l]=undef;
			}
		}
		
		return qv;
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}