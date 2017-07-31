/**
 * @(#)GlobalVelocityField.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.basic;

import miniufo.application.basic.VelocityFieldInSC;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import static miniufo.diagnosis.SpatialModel.GRAVITY_ACCERLERATION;


/**
 * analysis methods of global velocity fields in spherical coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GlobalVelocityFieldInSC extends VelocityFieldInSC{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public GlobalVelocityFieldInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isZonalPeriodic())
		throw new IllegalArgumentException("Not a zonal periodic model");
	}
	
	
	/**
     * calculate geostrophic velocity,
     * [0] is in x direction while [1] is in y direction
     *
     * @param	hgt		geopotential height
     *
     * @return	geostrophic velocity (vector)
     */
	public Variable[] cGeostrophicVelocity(Variable hgt){
		t=hgt.getTCount();	z=hgt.getZCount();	y=hgt.getYCount();	x=hgt.getXCount();
		
		float undef=hgt.getUndef();
		Variable[] geow=new Variable[2];
		geow[0]=new Variable("Ug",hgt);	geow[0].setUndef(undef);
		geow[1]=new Variable("Vg",hgt);	geow[1].setUndef(undef);
		
		geow[0].setCommentAndUnit("geostrophic Velocity in x-direction (m s^-1)");
		geow[1].setCommentAndUnit("geostrophic Velocity in y-direction (m s^-1)");
		
		float[][][][]  Ugdata=geow[0].getData();
		float[][][][]  Vgdata=geow[1].getData();
		float[][][][] hgtdata=    hgt.getData();
		
		if(hgt.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** Ug ***/
				for(int i=0;i<x;i++){
					for(int j=1;j<y-1;j++){
						if(hgtdata[l][k][j+1][i]!=undef&&hgtdata[l][k][j-1][i]!=undef)
							Ugdata[l][k][j][i]=-(hgtdata[l][k][j+1][i]-hgtdata[l][k][j-1][i])*
							GRAVITY_ACCERLERATION/(dy*2*f1[j]);
							
						else Ugdata[l][k][j][i]=undef;
					}
					
					if(hgtdata[l][k][1][i]!=undef&&hgtdata[l][k][0][i]!=undef)
						Ugdata[l][k][0][i]=-(hgtdata[l][k][1][i]-hgtdata[l][k][0][i])*
						GRAVITY_ACCERLERATION/(dy*f1[0]);
						
					else Ugdata[l][k][0][i]=undef;
					
					if(hgtdata[l][k][y-1][i]!=undef&&hgtdata[l][k][y-2][i]!=undef)
						Ugdata[l][k][y-1][i]=-(hgtdata[l][k][y-1][i]-hgtdata[l][k][y-2][i])*
						GRAVITY_ACCERLERATION/(dy*f1[y-1]);
						
					else Ugdata[l][k][y-1][i]=undef;
				}
				
				/*** Vg ***/
				for(int j=0;j<y;j++){
					for(int i=1;i<x-1;i++){
						if(hgtdata[l][k][j][i+1]!=undef&&hgtdata[l][k][j][i-1]!=undef)
							Vgdata[l][k][j][i]=(hgtdata[l][k][j][i+1]-hgtdata[l][k][j][i-1])*
							GRAVITY_ACCERLERATION/(dxs[j]*2*f1[j]);
						
						else Vgdata[l][k][j][i]=undef;
					}
					
					if(hgtdata[l][k][j][1]!=undef&&hgtdata[l][k][j][x-1]!=undef)
						Vgdata[l][k][j][0]=(hgtdata[l][k][j][1]-hgtdata[l][k][j][x-1])*
						GRAVITY_ACCERLERATION/(dxs[j]*2*f1[j]);
					
					else Vgdata[l][k][j][0]=undef;
					
					if(hgtdata[l][k][j][0]!=undef&&hgtdata[l][k][j][x-2]!=undef)
						Vgdata[l][k][j][x-1]=(hgtdata[l][k][j][0]-hgtdata[l][k][j][x-2])*
						GRAVITY_ACCERLERATION/(dxs[j]*2*f1[j]);
					
					else Vgdata[l][k][j][x-1]=undef;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				/*** Ug ***/
				for(int i=0;i<x;i++){
					for(int j=1;j<y-1;j++){
						if(hgtdata[k][j+1][i][l]!=undef&&hgtdata[k][j-1][i][l]!=undef)
							Ugdata[k][j][i][l]=-(hgtdata[k][j+1][i][l]-hgtdata[k][j-1][i][l])*
							GRAVITY_ACCERLERATION/(dy*2*f1[j]);
							
						else Ugdata[k][j][i][l]=undef;
					}
					
					if(hgtdata[k][1][i][l]!=undef&&hgtdata[k][0][i][l]!=undef)
						Ugdata[k][0][i][l]=-(hgtdata[k][1][i][l]-hgtdata[k][0][i][l])*
						GRAVITY_ACCERLERATION/(dy*f1[0]);
						
					else Ugdata[k][0][i][l]=undef;
					
					if(hgtdata[k][y-1][i][l]!=undef&&hgtdata[k][y-2][i][l]!=undef)
						Ugdata[k][y-1][i][l]=-(hgtdata[k][y-1][i][l]-hgtdata[k][y-2][i][l])*
						GRAVITY_ACCERLERATION/(dy*f1[y-1]);
						
					else Ugdata[k][y-1][i][l]=undef;
				}
				
				/*** Vg ***/
				for(int j=0;j<y;j++){
					for(int i=1;i<x-1;i++){
						if(hgtdata[k][j][i+1][l]!=undef&&hgtdata[k][j][i-1][l]!=undef)
							Vgdata[k][j][i][l]=(hgtdata[k][j][i+1][l]-hgtdata[k][j][i-1][l])*
							GRAVITY_ACCERLERATION/(dxs[j]*2*f1[j]);
						
						else Vgdata[k][j][i][l]=undef;
					}
					
					if(hgtdata[k][j][1][l]!=undef&&hgtdata[k][j][x-1][l]!=undef)
						Vgdata[k][j][0][l]=(hgtdata[k][j][1][l]-hgtdata[k][j][x-1][l])*
						GRAVITY_ACCERLERATION/(dxs[j]*2*f1[j]);
					
					else Vgdata[k][j][0][l]=undef;
					
					if(hgtdata[k][j][0][l]!=undef&&hgtdata[k][j][x-2][l]!=undef)
						Vgdata[k][j][x-1][l]=(hgtdata[k][j][0][l]-hgtdata[k][j][x-2][l])*
						GRAVITY_ACCERLERATION/(dxs[j]*2*f1[j]);
					
					else Vgdata[k][j][x-1][l]=undef;
				}
			}
		}
		
		return geow;
	}
	
	/**
     * calculate Velocity vector by using stream function,
     * [0] is in x direction while [1] is in y direction
     *
     * @param	sf	stream function
     *
     * @return	velocity (vector)
     */
	public Variable[] cRotationalVelocity(Variable sf){
		t=sf.getTCount();	z=sf.getZCount();	y=sf.getYCount();	x=sf.getXCount();
		
		float undef=sf.getUndef();
		Variable[] geow=new Variable[2];
		geow[0]=new Variable("Usf",sf);	geow[0].setCommentAndUnit("rotational Velocity in y-direction (m s^-1)");
		geow[1]=new Variable("Vsf",sf);	geow[1].setCommentAndUnit("rotational Velocity in y-direction (m s^-1)");
		
		float[][][][] Usdata=geow[0].getData();
		float[][][][] Vsdata=geow[1].getData();
		float[][][][] sfdata=     sf.getData();
		
		if(sf.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				// east boundary
				if(sfdata[l][k][j+1][0]!=undef&&sfdata[l][k][j-1][0]!=undef)
					Usdata[l][k][j][0]=-(sfdata[l][k][j+1][0]-sfdata[l][k][j-1][0])/(dy*2);
					
				else Usdata[l][k][j][0]=undef;
				
				if(sfdata[l][k][j][1]!=undef&&sfdata[l][k][j][x-1]!=undef)
					Vsdata[l][k][j][0]=(sfdata[l][k][j][1]-sfdata[l][k][j][x-1])/(dxs[j]+dxs[j]);
				
				else Vsdata[l][k][j][0]=undef;
				
				
				for(int i=1;i<x-1;i++){
					if(sfdata[l][k][j+1][i]!=undef&&sfdata[l][k][j-1][i]!=undef)
						Usdata[l][k][j][i]=-(sfdata[l][k][j+1][i]-sfdata[l][k][j-1][i])/(dy*2);
						
					else Usdata[l][k][j][i]=undef;
					
					if(sfdata[l][k][j][i+1]!=undef&&sfdata[l][k][j][i-1]!=undef)
						Vsdata[l][k][j][i]=(sfdata[l][k][j][i+1]-sfdata[l][k][j][i-1])/(dxs[j]+dxs[j]);
					
					else Vsdata[l][k][j][i]=undef;
				}
				
				// west boundary
				if(sfdata[l][k][j+1][x-1]!=undef&&sfdata[l][k][j-1][x-1]!=undef)
					Usdata[l][k][j][x-1]=-(sfdata[l][k][j+1][x-1]-sfdata[l][k][j-1][x-1])/(dy*2);
					
				else Usdata[l][k][j][x-1]=undef;
				
				if(sfdata[l][k][j][0]!=undef&&sfdata[l][k][j][x-2]!=undef)
					Vsdata[l][k][j][x-1]=(sfdata[l][k][j][0]-sfdata[l][k][j][x-2])/(dxs[j]+dxs[j]);
				
				else Vsdata[l][k][j][x-1]=undef;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				// east boundary
				if(sfdata[k][j+1][0][l]!=undef&&sfdata[k][j-1][0][l]!=undef)
					Usdata[k][j][0][l]=-(sfdata[k][j+1][0][l]-sfdata[k][j-1][0][l])/(dy*2);
					
				else Usdata[k][j][0][l]=undef;
				
				if(sfdata[k][j][1][l]!=undef&&sfdata[k][j][x-1][l]!=undef)
					Vsdata[k][j][0][l]=(sfdata[k][j][1][l]-sfdata[k][j][x-1][l])/(dxs[j]+dxs[j]);
				
				else Vsdata[l][k][j][0]=undef;
				
				
				for(int i=1;i<x-1;i++){
					if(sfdata[k][j+1][i][l]!=undef&&sfdata[k][j-1][i][l]!=undef)
						Usdata[k][j][i][l]=-(sfdata[k][j+1][i][l]-sfdata[k][j-1][i][l])/(dy*2);
						
					else Usdata[k][j][i][l]=undef;
					
					if(sfdata[k][j][i+1][l]!=undef&&sfdata[k][j][i-1][l]!=undef)
						Vsdata[k][j][i][l]=(sfdata[k][j][i+1][l]-sfdata[k][j][i-1][l])/(dxs[j]+dxs[j]);
					
					else Vsdata[l][k][j][i]=undef;
				}
				
				// west boundary
				if(sfdata[k][j+1][x-1][l]!=undef&&sfdata[k][j-1][x-1][l]!=undef)
					Usdata[k][j][x-1][l]=-(sfdata[k][j+1][x-1][l]-sfdata[k][j-1][x-1][l])/(dy*2);
					
				else Usdata[k][j][x-1][l]=undef;
				
				if(sfdata[k][j][0][l]!=undef&&sfdata[k][j][x-2][l]!=undef)
					Vsdata[k][j][x-1][l]=(sfdata[k][j][0][l]-sfdata[k][j][x-2][l])/(dxs[j]+dxs[j]);
				
				else Vsdata[l][k][j][x-1]=undef;
			}
		}
		
		return geow;
	}
	
	
	/**
     * calculate the stream function using SOR
     *
     * @param	u	u-velocity
     * @param	v	v-velocity
     *
     * @return	stream function
     */
	public Variable cStreamFunctionBySOR(Variable u,Variable v){
		SphericalSpatialModel 	  ssm=(SphericalSpatialModel)(sm);
		GlobalDynamicMethodsInSC  gdm=new GlobalDynamicMethodsInSC(ssm);
		GlobalLaplaceEquationInSC gle=new GlobalLaplaceEquationInSC(ssm);
		
		Variable vor=gdm.c2DVorticity(u,v);
		Variable sf=new Variable("sf",u);
		sf.setUndef(u.getUndef());	sf.setCommentAndUnit("stream function (m^2 s^-1)");
		
		gle.setMaxLoopCount(3000);
		gle.setTolerance(500);
		gle.setThreadCount(threadCount);
		gle.solve(sf,vor);
	    
	    return sf;
	}
	
	/**
     * calculate the potential function using SOR
     *
     * @param	u	u-velocity
     * @param	v	v-velocity
     *
     * @return	velocity potential
     */
	public Variable cVelocityPotentialBySOR(Variable u,Variable v){
		SphericalSpatialModel 	  ssm=(SphericalSpatialModel)(sm);
		GlobalDynamicMethodsInSC  gdm=new GlobalDynamicMethodsInSC(ssm);
		GlobalLaplaceEquationInSC gle=new GlobalLaplaceEquationInSC(ssm);
		
		Variable div=gdm.c2DDivergence(u,v);
		Variable pf=new Variable("pf",u);
		pf.setUndef(u.getUndef());	pf.setCommentAndUnit("potential function (m^2 s^-1)");
		
		gle.setMaxLoopCount(3000);
		gle.setTolerance(100);
		gle.setThreadCount(threadCount);
		gle.solve(pf,div);
	    
	    return pf;
	}
	
	
	/** test
	public static void main(String arg[]){
		try{
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
