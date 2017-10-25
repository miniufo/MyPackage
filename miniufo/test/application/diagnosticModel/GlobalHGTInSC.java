/**
 * @(#)GlobalHGTInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.diagnosticModel;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.test.application.GlobalEllipticalEquationInSpheralCoordinate;
import miniufo.test.application.basic.GlobalDynamicMethodsInSC;
import miniufo.test.application.basic.GlobalLaplaceEquationInSC;
import miniufo.application.basic.SphericalHarmonicExpansion;
import static java.lang.Math.cos;
import static miniufo.geophysics.atmos.ThermoDynamics.Rd;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.diagnosis.SpatialModel.GRAVITY_ACCERLERATION;


/**
 * global geopotential height diagnosis equation in spheral coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GlobalHGTInSC extends GlobalEllipticalEquationInSpheralCoordinate{
	//
	protected int threadCount=Runtime.getRuntime().availableProcessors()-1;
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public GlobalHGTInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isPeriodicX())
		throw new IllegalArgumentException("Not a zonal periodic model");
		
		le=new GlobalLaplaceEquationInSC(ssm);
	}
	
	
	/*** getor and setor ***/
	public void setThreadCount(int count){
		if(count<1)
			throw new IllegalArgumentException("tolerance should be positive");
		
		threadCount=count;
	}
	
	
	/**
     * calculate forces one by one
     */
	public Variable cTerm1(Variable div){
		Variable F=new Variable("F1",div);	F.setCommentAndUnit("local partial of divergence");
		
		cTerm1(F,div);
		
		return F;
	}
	
	public Variable cTerm2(Variable uwnd){
		Variable F=new Variable("F2",uwnd);	F.setCommentAndUnit("x-jet term");
		
		cTerm2(F,uwnd);
		
		return F;
	}
	
	public Variable cTerm3(Variable uwnd,Variable vwnd){
		Variable F=new Variable("F3",uwnd);	F.setCommentAndUnit("saddle pattern flow");
		
		cTerm3(F,uwnd,vwnd);
		
		return F;
	}

	public Variable cTerm4(Variable uwnd,Variable vwnd){
		Variable F=new Variable("F4",uwnd);	F.setCommentAndUnit("saddle pattern flow");
		
		cTerm4(F,uwnd,vwnd);
		
		return F;
	}
	
	public Variable cTerm5(Variable vwnd){
		Variable F=new Variable("F5",vwnd);	F.setCommentAndUnit("y-jet term");
		
		cTerm5(F,vwnd);
		
		return F;
	}
	
	public Variable cTerm6(Variable omega,Variable div){
		Variable F=new Variable("F6",omega);	F.setCommentAndUnit("convection of divergence");
		
		cTerm6(F,omega,div);
		
		return F;
	}
	
	public Variable cTerm7(Variable uwnd,Variable omega){
		Variable F=new Variable("F7",uwnd);	F.setCommentAndUnit("Walker circulation term");
		
		cTerm7(F,uwnd,omega);
		
		return F;
	}
	
	public Variable cTerm8(Variable vwnd,Variable omega){
		Variable F=new Variable("F8",vwnd);	F.setCommentAndUnit("Hadley circulation term");
		
		cTerm8(F,vwnd,omega);
		
		return F;
	}
	
	public Variable cTerm77(Variable uwnd,Variable omega){
		Variable F=new Variable("F7",uwnd);	F.setCommentAndUnit("saddle pattern in X-P section");
		
		cTerm77(F,uwnd,omega);
		
		return F;
	}
	
	public Variable cTerm88(Variable vwnd,Variable omega){
		Variable F=new Variable("F8",vwnd);	F.setCommentAndUnit("saddle pattern in Y-P section");
		
		cTerm88(F,vwnd,omega);
		
		return F;
	}
	
	public Variable cTerm9(Variable vor){
		Variable F=new Variable("F9",vor);	F.setCommentAndUnit("geostropic term");
		
		cTerm9(F,vor);
		
		return F;
	}
	
	public Variable cTerm10(Variable uwnd){
		Variable F=new Variable("F10",uwnd);	F.setCommentAndUnit("planetary vorticity term");
		
		cTerm10(F,uwnd);
		
		return F;
	}
	
	public Variable cTerm11(Variable uwnd,Variable vwnd){
		Variable F=new Variable("F11",uwnd);	F.setCommentAndUnit("spherical term");
		
		cTerm11(F,uwnd,vwnd);
		
		return F;
	}
	
	public Variable cTerm12(Variable uwnd){
		Variable F=new Variable("F12",uwnd);	F.setCommentAndUnit("spherical term");
		
		cTerm12(F,uwnd);
		
		return F;
	}
	
	public Variable cTerm13(Variable omega,Variable div,Variable air){
		Variable F=new Variable("F13",omega);	F.setCommentAndUnit("thermal term");
		
		cTerm13(F,omega,div,air);
		
		return F;
	}
	
	public Variable cTerm14(Variable uwnd,Variable vwnd,Variable omega,Variable air){
		Variable F=new Variable("F14",uwnd);	F.setCommentAndUnit("thermal term");
		
		cTerm14(F,uwnd,vwnd,omega,air);
		
		return F;
	}
	
	public Variable cTerm15(Variable omega,Variable air){
		Variable F=new Variable("F15",omega);	F.setCommentAndUnit("thermal term");
		
		cTerm15(F,omega,air);
		
		return F;
	}
	
	/**
     * calculate all forces
     */
	public Variable cAllTerms(Variable uwnd,Variable vwnd,Variable omega,Variable div,Variable vor,Variable air){
		Variable F=new Variable("F",omega);	F.setCommentAndUnit("all forcing term");
		
		cAllTerms(F,uwnd,vwnd,omega,div,vor,air);
		
		return F;
	}
	
	
	/**
     * calculate forces one by one
     */
	public void cTerm1(Variable F,Variable div){
		if(!F.isLike(div)) throw new IllegalArgumentException("dimensions not same");
		if(dt<=0) throw new IllegalArgumentException("delta t should be larger than 0");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]   Fdata=  F.getData();
		float[][][][] divdata=div.getData();
		
		if(F.isTFirst()){
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[l][k][j][i]=-(divdata[l+1][k][j][i]-divdata[l-1][k][j][i])/(2*dt);	//2*dt
			
		}else{
			for(int l=1;l<t-1;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[k][j][i][l]=-(divdata[k][j][i][l+1]-divdata[k][j][i][l-1])/(2*dt);	//2*dt
		}
	}
	
	public void cTerm2(Variable F,Variable uwnd){
		if(!F.isLike(uwnd)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]    Fdata=   F.getData();
		float[][][][] uwnddata=uwnd.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[l][k][j][i]=-(
						(uwnddata[l][k][j][i+1]+uwnddata[l][k][j][i])*(uwnddata[l][k][j][i+1]-uwnddata[l][k][j][i])-
						(uwnddata[l][k][j][i]+uwnddata[l][k][j][i-1])*(uwnddata[l][k][j][i]-uwnddata[l][k][j][i-1])
					)/(dxs[j]*dxs[j]*2);
				}
				
				/*** East boundary ***/
				Fdata[l][k][j][0]=-(
					(uwnddata[l][k][j][1]+uwnddata[l][k][j][0])*(uwnddata[l][k][j][1]-uwnddata[l][k][j][0])-
					(uwnddata[l][k][j][0]+uwnddata[l][k][j][x-1])*(uwnddata[l][k][j][0]-uwnddata[l][k][j][x-1])
				)/(dxs[j]*dxs[j]*2);
				
				/*** West boundary ***/
				Fdata[l][k][j][x-1]=-(
					(uwnddata[l][k][j][0]+uwnddata[l][k][j][x-1])*(uwnddata[l][k][j][0]-uwnddata[l][k][j][x-1])-
					(uwnddata[l][k][j][x-1]+uwnddata[l][k][j][x-2])*(uwnddata[l][k][j][x-1]-uwnddata[l][k][j][x-2])
				)/(dxs[j]*dxs[j]*2);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[k][j][i][l]=-(
						(uwnddata[k][j][i+1][l]+uwnddata[k][j][i][l])*(uwnddata[k][j][i+1][l]-uwnddata[k][j][i][l])-
						(uwnddata[k][j][i][l]+uwnddata[k][j][i-1][l])*(uwnddata[k][j][i][l]-uwnddata[k][j][i-1][l])
					)/(dxs[j]*dxs[j]*2);
				}
				
				/*** East boundary ***/
				Fdata[k][j][0][l]=-(
					(uwnddata[k][j][1][l]+uwnddata[k][j][0][l])*(uwnddata[k][j][1][l]-uwnddata[k][j][0][l])-
					(uwnddata[k][j][0][l]+uwnddata[k][j][x-1][l])*(uwnddata[k][j][0][l]-uwnddata[k][j][x-1][l])
				)/(dxs[j]*dxs[j]*2);
				
				/*** West boundary ***/
				Fdata[k][j][x-1][l]=-(
					(uwnddata[k][j][0][l]+uwnddata[k][j][x-1][l])*(uwnddata[k][j][0][l]-uwnddata[k][j][x-1][l])-
					(uwnddata[k][j][x-1][l]+uwnddata[k][j][x-2][l])*(uwnddata[k][j][x-1][l]-uwnddata[k][j][x-2][l])
				)/(dxs[j]*dxs[j]*2);
			}
		}
	}
	
	public void cTerm3(Variable F,Variable uwnd,Variable vwnd){
		if(!F.isLike(uwnd)||!F.isLike(vwnd)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]    Fdata=   F.getData();
		float[][][][] uwnddata=uwnd.getData();
		float[][][][] vwnddata=vwnd.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++)
				Fdata[l][k][j][i]=-(
					uwnddata[l][k][j+1][i]*(vwnddata[l][k][j+1][i+1]-vwnddata[l][k][j+1][i-1])-
					uwnddata[l][k][j-1][i]*(vwnddata[l][k][j-1][i+1]-vwnddata[l][k][j-1][i-1])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				
				/*** East boundary ***/
				Fdata[l][k][j][0]=-(
					 uwnddata[l][k][j+1][0]*(vwnddata[l][k][j+1][1]-vwnddata[l][k][j+1][x-1])-
					 uwnddata[l][k][j-1][0]*(vwnddata[l][k][j-1][1]-vwnddata[l][k][j-1][x-1])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				
				/*** West boundary ***/
				Fdata[l][k][j][x-1]=-(
					 uwnddata[l][k][j+1][x-1]*(vwnddata[l][k][j+1][0]-vwnddata[l][k][j+1][x-2])-
					 uwnddata[l][k][j-1][x-1]*(vwnddata[l][k][j-1][0]-vwnddata[l][k][j-1][x-2])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++)
				Fdata[k][j][i][l]=-(
					uwnddata[k][j+1][i][l]*(vwnddata[k][j+1][i+1][l]-vwnddata[k][j+1][i-1][l])-
					uwnddata[k][j-1][i][l]*(vwnddata[k][j-1][i+1][l]-vwnddata[k][j-1][i-1][l])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				
				/*** East boundary ***/
				Fdata[k][j][0][l]=-(
					 uwnddata[k][j+1][0][l]*(vwnddata[k][j+1][1][l]-vwnddata[k][j+1][x-1][l])-
					 uwnddata[k][j-1][0][l]*(vwnddata[k][j-1][1][l]-vwnddata[k][j-1][x-1][l])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				
				/*** West boundary ***/
				Fdata[k][j][x-1][l]=-(
					 uwnddata[k][j+1][x-1][l]*(vwnddata[k][j+1][0][l]-vwnddata[k][j+1][x-2][l])-
					 uwnddata[k][j-1][x-1][l]*(vwnddata[k][j-1][0][l]-vwnddata[k][j-1][x-2][l])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
			}
		}
	}
	
	public void cTerm4(Variable F,Variable uwnd,Variable vwnd){
		if(!F.isLike(uwnd)||!F.isLike(vwnd)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]    Fdata=   F.getData();
		float[][][][] uwnddata=uwnd.getData();
		float[][][][] vwnddata=vwnd.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[l][k][j][i]=-(
						vwnddata[l][k][j][i+1]*(uwnddata[l][k][j+1][i+1]-uwnddata[l][k][j-1][i+1])-
						vwnddata[l][k][j][i-1]*(uwnddata[l][k][j+1][i-1]-uwnddata[l][k][j-1][i-1])
					)/(dxs[j]*dy*4);
				}
				
				/*** East boundary ***/
				Fdata[l][k][j][0]=-(
					vwnddata[l][k][j][  1]*(uwnddata[l][k][j+1][  1]-uwnddata[l][k][j-1][  1])-
					vwnddata[l][k][j][x-1]*(uwnddata[l][k][j+1][x-1]-uwnddata[l][k][j-1][x-1])
				)/(dxs[j]*dy*4);
				
				/*** West boundary ***/
				Fdata[l][k][j][x-1]=-(
					vwnddata[l][k][j][  0]*(uwnddata[l][k][j+1][  0]-uwnddata[l][k][j-1][  0])-
					vwnddata[l][k][j][x-2]*(uwnddata[l][k][j+1][x-2]-uwnddata[l][k][j-1][x-2])
				)/(dxs[j]*dy*4);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[k][j][i][l]=-(
						vwnddata[k][j][i+1][l]*(uwnddata[k][j+1][i+1][l]-uwnddata[k][j-1][i+1][l])-
						vwnddata[k][j][i-1][l]*(uwnddata[k][j+1][i-1][l]-uwnddata[k][j-1][i-1][l])
					)/(dxs[j]*dy*4);
				}
				
				/*** East boundary ***/
				Fdata[k][j][0][l]=-(
					vwnddata[k][j][  1][l]*(uwnddata[k][j+1][  1][l]-uwnddata[k][j-1][  1][l])-
					vwnddata[k][j][x-1][l]*(uwnddata[k][j+1][x-1][l]-uwnddata[k][j-1][x-1][l])
				)/(dxs[j]*dy*4);
				
				/*** West boundary ***/
				Fdata[k][j][x-1][l]=-(
					vwnddata[k][j][  0][l]*(uwnddata[k][j+1][  0][l]-uwnddata[k][j-1][  0][l])-
					vwnddata[k][j][x-2][l]*(uwnddata[k][j+1][x-2][l]-uwnddata[k][j-1][x-2][l])
				)/(dxs[j]*dy*4);
			}
		}
	}
	
	public void cTerm5(Variable F,Variable vwnd){
		if(!F.isLike(vwnd)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]    Fdata=   F.getData();
		float[][][][] vwnddata=vwnd.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++){
				Fdata[l][k][j][i]=-(
					(vwnddata[l][k][j+1][i]+vwnddata[l][k][j][i])*
					(float)cos((ydef[j+1]+ydef[j])/2)*
					(vwnddata[l][k][j+1][i]-vwnddata[l][k][j][i])
					-
					(vwnddata[l][k][j][i]+vwnddata[l][k][j-1][i])*
					(float)cos((ydef[j]+ydef[j-1])/2)*
					(vwnddata[l][k][j][i]-vwnddata[l][k][j-1][i])
				)/(dy*dy*2*lcos[j]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++){
				Fdata[k][j][i][l]=-(
					(vwnddata[k][j+1][i][l]+vwnddata[k][j][i][l])*(float)cos((ydef[j+1]+ydef[j])/2)*
					(vwnddata[k][j+1][i][l]-vwnddata[k][j][i][l])-
					(vwnddata[k][j][i][l]+vwnddata[k][j-1][i][l])*(float)cos((ydef[j]+ydef[j-1])/2)*
					(vwnddata[k][j][i][l]-vwnddata[k][j-1][i][l])
				)/(dy*dy*2*lcos[j]);
			}
		}
	}
	
	public void cTerm6(Variable F,Variable omega,Variable div){
		if(!F.isLike(omega)||!F.isLike(div)) throw new IllegalArgumentException("dimensions not same");
		
		zstart=F.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]   divdata=  div.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[l][k][j][i]=-omegadata[l][k][j][i]*
				(divdata[l][k+1][j][i]-divdata[l][k-1][j][i])/(dz+dz);
			
		}else{
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[k][j][i][l]=-omegadata[k][j][i][l]*
				(divdata[k+1][j][i][l]-divdata[k-1][j][i][l])/(dz+dz);
		}
	}
	
	public void cTerm77(Variable F,Variable uwnd,Variable omega){
		if(!F.isLike(omega)||!F.isLike(uwnd)) throw new IllegalArgumentException("dimensions not same");
		
		zstart=F.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]     udata= uwnd.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				for(int i=1;i<x-1;i++) Fdata[l][k][j][i]=-(
					(udata[l][k+1][j][i+1]-udata[l][k-1][j][i+1])*omegadata[l][k][j][i+1]-
					(udata[l][k+1][j][i-1]-udata[l][k-1][j][i-1])*omegadata[l][k][j][i-1]
				)/(4*dz*dxs[j]);
				
				Fdata[l][k][j][0]=-(
					(udata[l][k+1][j][  1]-udata[l][k-1][j][  1])*omegadata[l][k][j][  1]-
					(udata[l][k+1][j][x-1]-udata[l][k-1][j][x-1])*omegadata[l][k][j][x-1]
				)/(4*dz*dxs[j]);
				
				Fdata[l][k][j][x-1]=-(
					(udata[l][k+1][j][  0]-udata[l][k-1][j][  0])*omegadata[l][k][j][  0]-
					(udata[l][k+1][j][x-2]-udata[l][k-1][j][x-2])*omegadata[l][k][j][x-2]
				)/(4*dz*dxs[j]);
			}
			
		}else{
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++){
				for(int i=1;i<x-1;i++) Fdata[k][j][i][l]=-(
					(udata[k+1][j][i+1][l]-udata[k-1][j][i+1][l])*omegadata[k][j][i+1][l]-
					(udata[k+1][j][i-1][l]-udata[k-1][j][i-1][l])*omegadata[k][j][i-1][l]
				)/(4*dz*dxs[j]);
				
				Fdata[l][k][j][0]=-(
					(udata[k+1][j][  1][l]-udata[k-1][j][  1][l])*omegadata[k][j][  1][l]-
					(udata[k+1][j][x-1][l]-udata[k-1][j][x-1][l])*omegadata[k][j][x-1][l]
				)/(4*dz*dxs[j]);
				
				Fdata[l][k][j][x-1]=-(
					(udata[k+1][j][  0][l]-udata[k-1][j][  0][l])*omegadata[k][j][  0][l]-
					(udata[k+1][j][x-2][l]-udata[k-1][j][x-2][l])*omegadata[k][j][x-2][l]
				)/(4*dz*dxs[j]);
			}
		}
	}
	
	public void cTerm7(Variable F,Variable uwnd,Variable omega){
		if(!F.isLike(uwnd)||!F.isLike(omega)) throw new IllegalArgumentException("dimensions not same");

		zstart=F.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]  uwnddata= uwnd.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[l][k][j][i]=-
						(omegadata[l][k][j][i+1]-omegadata[l][k][j][i-1])/(dxs[j]*2)*
						( uwnddata[l][k+1][j][i]- uwnddata[l][k-1][j][i])/(dz+dz);
				}
				
				/*** East boundary ***/
				Fdata[l][k][j][0]=-
					(omegadata[l][k][j][1]-omegadata[l][k][j][x-1])/(dxs[j]*2)*
					(uwnddata[l][k+1][j][0]-uwnddata[l][k-1][j][0])/(dz+dz);
				
				/*** West boundary ***/
				Fdata[l][k][j][x-1]=-
					(omegadata[l][k][j][0]-omegadata[l][k][j][x-2])/(dxs[j]*2)*
					(uwnddata[l][k+1][j][x-1]-uwnddata[l][k-1][j][x-1])/(dz+dz);
			}
			
		}else{
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[k][j][i][l]=-
						(omegadata[k][j][i+1][l]-omegadata[k][j][i-1][l])/(dxs[j]*2)*
						( uwnddata[k+1][j][i][l]- uwnddata[k-1][j][i][l])/(dz+dz);
				}
				
				/*** East boundary ***/
				Fdata[k][j][0][l]=-
					(omegadata[k][j][1][l]-omegadata[k][j][x-1][l])/(dxs[j]*2)*
					(uwnddata[k+1][j][0][l]-uwnddata[k-1][j][0][l])/(dz+dz);
				
				/*** West boundary ***/
				Fdata[k][j][x-1][l]=-
					(omegadata[k][j][0][l]-omegadata[k][j][x-2][l])/(dxs[j]*2)*
					(uwnddata[k+1][j][x-1][l]-uwnddata[k-1][j][x-1][l])/(dz+dz);
			}
		}
	}
	
	public void cTerm88(Variable F,Variable vwnd,Variable omega){
		if(!F.isLike(vwnd)||!F.isLike(omega)) throw new IllegalArgumentException("dimensions not same");

		zstart=F.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]     vdata= vwnd.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				for(int i=0;i<x;i++) Fdata[l][k][j][i]=-(
					(vdata[l][k+1][j+1][i]-vdata[l][k-1][j+1][i])*omegadata[l][k][j+1][i]-
					(vdata[l][k+1][j-1][i]-vdata[l][k-1][j-1][i])*omegadata[l][k][j-1][i]
				)/(4*dz*dy);
			}
			
		}else{
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++){
				for(int i=0;i<x;i++) Fdata[k][j][i][l]=-(
					(vdata[k+1][j+1][i][l]-vdata[k-1][j+1][i][l])*omegadata[k][j+1][i][l]-
					(vdata[k+1][j-1][i][l]-vdata[k-1][j-1][i][l])*omegadata[k][j-1][i][l]
				)/(4*dz*dy);
			}
		}
	}
	
	public void cTerm8(Variable F,Variable vwnd,Variable omega){
		if(!F.isLike(vwnd)||!F.isLike(omega)) throw new IllegalArgumentException("dimensions not same");

		zstart=F.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]  vwnddata= vwnd.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++){
				Fdata[l][k][j][i]=-
					(omegadata[l][k][j+1][i]-omegadata[l][k][j-1][i])/(dy*2)*
					( vwnddata[l][k+1][j][i]- vwnddata[l][k-1][j][i])/(dz+dz);
			}
			
		}else{
			for(int k=1;k<z-1;k++)
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++){
				Fdata[k][j][i][l]=-
					(omegadata[k][j+1][i][l]-omegadata[k][j-1][i][l])/(dy*2)*
					( vwnddata[k+1][j][i][l]- vwnddata[k-1][j][i][l])/(dz+dz);
			}
		}
	}
	
	public void cTerm9(Variable F,Variable vor){
		if(!F.isLike(vor)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]   Fdata=  F.getData();
		float[][][][] vordata=vor.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[l][k][j][i]=f1[j]*vordata[l][k][j][i];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[k][j][i][l]=f1[j]*vordata[k][j][i][l];
		}
	}
	
	public void cTerm10(Variable F,Variable uwnd){
		if(!F.isLike(uwnd)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]    Fdata=   F.getData();
		float[][][][] uwnddata=uwnd.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[l][k][j][i]=-uwnddata[l][k][j][i]*f2[j]/EARTH_RADIUS;
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++) Fdata[k][j][i][l]=-uwnddata[k][j][i][l]*f2[j]/EARTH_RADIUS;
		}
	}
	
	public void cTerm11(Variable F,Variable uwnd,Variable vwnd){
		if(!F.isLike(uwnd)||!F.isLike(vwnd)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]    Fdata=   F.getData();
		float[][][][] uwnddata=uwnd.getData();
		float[][][][] vwnddata=vwnd.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[l][k][j][i]=(
						uwnddata[l][k][j][i+1]*vwnddata[l][k][j][i+1]-
						uwnddata[l][k][j][i-1]*vwnddata[l][k][j][i-1]
					)*ltan[j]/(dxs[j]*2)/EARTH_RADIUS;
				}
				
				/*** East boundary ***/
				Fdata[l][k][j][0]=(
					uwnddata[l][k][j][1]*vwnddata[l][k][j][1]-
					uwnddata[l][k][j][x-1]*vwnddata[l][k][j][x-1]
				)*ltan[j]/(dxs[j]*2)/EARTH_RADIUS;
				
				/*** West boundary ***/
				Fdata[l][k][j][x-1]=(
					uwnddata[l][k][j][0]*vwnddata[l][k][j][0]-
					uwnddata[l][k][j][x-2]*vwnddata[l][k][j][x-2]
				)*ltan[j]/(dxs[j]*2)/EARTH_RADIUS;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[k][j][i][l]=(
						uwnddata[k][j][i+1][l]*vwnddata[k][j][i+1][l]-
						uwnddata[k][j][i-1][l]*vwnddata[k][j][i-1][l]
					)*ltan[j]/(dxs[j]*2)/EARTH_RADIUS;
				}
				
				/*** East boundary ***/
				Fdata[k][j][0][l]=(
					uwnddata[k][j][1][l]*vwnddata[k][j][1][l]-
					uwnddata[k][j][x-1][l]*vwnddata[k][j][x-1][l]
				)*ltan[j]/(dxs[j]*2)/EARTH_RADIUS;
				
				/*** West boundary ***/
				Fdata[k][j][x-1][l]=(
					uwnddata[k][j][0][l]*vwnddata[k][j][0][l]-
					uwnddata[k][j][x-2][l]*vwnddata[k][j][x-2][l]
				)*ltan[j]/(dxs[j]*2)/EARTH_RADIUS;
			}
		}
	}
	
	public void cTerm12(Variable F,Variable uwnd){
		if(!F.isLike(uwnd)) throw new IllegalArgumentException("dimensions not same");
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]    Fdata=   F.getData();
		float[][][][] uwnddata=uwnd.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++){
				Fdata[l][k][j][i]=-(
					uwnddata[l][k][j+1][i]*uwnddata[l][k][j+1][i]*lsin[j+1]-
					uwnddata[l][k][j-1][i]*uwnddata[l][k][j-1][i]*lsin[j-1]
				)/(dy*2)/lcos[j]/EARTH_RADIUS;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++){
				Fdata[k][j][i][l]=-(
					uwnddata[k][j+1][i][l]*uwnddata[k][j+1][i][l]*lsin[j+1]-
					uwnddata[k][j-1][i][l]*uwnddata[k][j-1][i][l]*lsin[j-1]
				)/(dy*2)/lcos[j]/EARTH_RADIUS;
			}
		}
	}
	
	public void cTerm13(Variable F,Variable omega,Variable div,Variable air){
		if(!F.isLike(omega)||!F.isLike(div)||!F.isLike(air)) throw new IllegalArgumentException("dimensions not same");

		zstart=F.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]   divdata=  div.getData();
		float[][][][]   airdata=  air.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
				Fdata[l][k][j][i]=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
					airdata[l][k][j][i]*omegadata[l][k][j][i]*divdata[l][k][j][i];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++)
			for(int i=0;i<x;i++)
				Fdata[k][j][i][l]=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
					airdata[k][j][i][l]*omegadata[k][j][i][l]*divdata[k][j][i][l];
		}
	}
	
	public void cTerm14(Variable F,Variable uwnd,Variable vwnd,Variable omega,Variable air){
		if(!F.isLike(uwnd)||!F.isLike(vwnd)||!F.isLike(omega)||!F.isLike(air)) throw new IllegalArgumentException("dimensions not same");

		zstart=F.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]   airdata=  air.getData();
		float[][][][]  uwnddata= uwnd.getData();
		float[][][][]  vwnddata= vwnd.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[l][k][j][i]=(
						(omegadata[l][k][j][i+1]*airdata[l][k][j][i+1]-omegadata[l][k][j][i-1]*airdata[l][k][j][i-1])*
						uwnddata[l][k][j][i]/(dxs[j]*2)+
						(omegadata[l][k][j+1][i]*airdata[l][k][j+1][i]-omegadata[l][k][j-1][i]*airdata[l][k][j-1][i])*
						vwnddata[l][k][j][i]/(dy*2)
					)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				}
				
				/*** East boundary ***/
				Fdata[l][k][j][0]=(
					(omegadata[l][k][j][1]*airdata[l][k][j][1]-omegadata[l][k][j][x-1]*airdata[l][k][j][x-1])*
					uwnddata[l][k][j][0]/(dxs[j]*2)+
					(omegadata[l][k][j+1][0]*airdata[l][k][j+1][0]-omegadata[l][k][j-1][0]*airdata[l][k][j-1][0])*
					vwnddata[l][k][j][0]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				
				/*** West boundary ***/
				Fdata[l][k][j][x-1]=(
					(omegadata[l][k][j][0]*airdata[l][k][j][0]-omegadata[l][k][j][x-2]*airdata[l][k][j][x-2])*
					uwnddata[l][k][j][x-1]/(dxs[j]*2)+
					(omegadata[l][k][j+1][x-1]*airdata[l][k][j+1][x-1]-omegadata[l][k][j-1][x-1]*airdata[l][k][j-1][x-1])*
					vwnddata[l][k][j][x-1]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[k][j][i][l]=(
						(omegadata[k][j][i+1][l]*airdata[k][j][i+1][l]-omegadata[k][j][i-1][l]*airdata[k][j][i-1][l])*
						uwnddata[k][j][i][l]/(dxs[j]*2)+
						(omegadata[k][j+1][i][l]*airdata[k][j+1][i][l]-omegadata[k][j-1][i][l]*airdata[k][j-1][i][l])*
						vwnddata[k][j][i][l]/(dy*2)
					)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				}
				
				/*** East boundary ***/
				Fdata[k][j][0][l]=(
					(omegadata[k][j][1][l]*airdata[k][j][1][l]-omegadata[k][j][x-1][l]*airdata[k][j][x-1][l])*
					uwnddata[k][j][0][l]/(dxs[j]*2)+
					(omegadata[k][j+1][0][l]*airdata[k][j+1][0][l]-omegadata[k][j-1][0][l]*airdata[k][j-1][0][l])*
					vwnddata[k][j][0][l]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				
				/*** West boundary ***/
				Fdata[k][j][x-1][l]=(
					(omegadata[k][j][0][l]*airdata[k][j][0][l]-omegadata[k][j][x-2][l]*airdata[k][j][x-2][l])*
					uwnddata[k][j][x-1][l]/(dxs[j]*2)+
					(omegadata[k][j+1][x-1][l]*airdata[k][j+1][x-1][l]-omegadata[k][j-1][x-1][l]*airdata[k][j-1][x-1][l])*
					vwnddata[k][j][x-1][l]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
			}
		}
	}
	
	public void cTerm15(Variable F,Variable omega,Variable air){
		if(!F.isLike(omega)||!F.isLike(air)) throw new IllegalArgumentException("dimensions not same");
		
		zstart=omega.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]   airdata=  air.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[l][k][j][i]=(
						(omegadata[l][k][j][i+1]*airdata[l][k][j][i+1]-omegadata[l][k][j][i-1]*airdata[l][k][j][i-1])*
						f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
					);
				}
				
				/*** East boundary ***/
				Fdata[l][k][j][0]=(
					(omegadata[l][k][j][1]*airdata[l][k][j][1]-omegadata[l][k][j][x-1]*airdata[l][k][j][x-1])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
				
				/*** West boundary ***/
				Fdata[l][k][j][x-1]=(
					(omegadata[l][k][j][0]*airdata[l][k][j][0]-omegadata[l][k][j][x-2]*airdata[l][k][j][x-2])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					Fdata[k][j][i][l]=(
						(omegadata[k][j][i+1][l]*airdata[k][j][i+1][l]-omegadata[k][j][i-1][l]*airdata[k][j][i-1][l])*
						f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
					);
				}
				
				/*** East boundary ***/
				Fdata[k][j][0][l]=(
					(omegadata[k][j][1][l]*airdata[k][j][1][l]-omegadata[k][j][x-1][l]*airdata[k][j][x-1][l])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
				
				/*** West boundary ***/
				Fdata[k][j][x-1][l]=(
					(omegadata[k][j][0][l]*airdata[k][j][0][l]-omegadata[k][j][x-2][l]*airdata[k][j][x-2][l])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
			}
		}
	}
	
	/**
     * calculate all forces
     */
	public void cAllTerms(Variable F,Variable uwnd,Variable vwnd,
	Variable omega,Variable div,Variable vor,Variable air){
		if(!F.isLike(uwnd)||!F.isLike(vwnd)||!F.isLike(omega)||!F.isLike(div)
		 ||!F.isLike(vor) ||!F.isLike(air)) throw new IllegalArgumentException("dimensions not same");
		if(dt<=0) throw new IllegalArgumentException();
		
		tstart=uwnd.getRange().getTRange()[0];
		zstart=uwnd.getRange().getZRange()[0];
		
		t=F.getTCount();	z=F.getZCount();	y=F.getYCount();	x=F.getXCount();
		
		float[][][][]     Fdata=    F.getData();
		float[][][][]   divdata=  div.getData();
		float[][][][]   vordata=  vor.getData();
		float[][][][]   airdata=  air.getData();
		float[][][][]  uwnddata= uwnd.getData();
		float[][][][]  vwnddata= vwnd.getData();
		float[][][][] omegadata=omega.getData();
		
		if(F.isTFirst()){
			for(int k=1;k<z-1;k++)
			for(int l=1;l<t-1;l++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					/************************ inside area ************************/
					Fdata[l][k][j][i]=-(divdata[l+1][k][j][i]-divdata[l-1][k][j][i])/(2*dt);	//2*dt
					/**/
					Fdata[l][k][j][i]+=-(
						(uwnddata[l][k][j][i+1]+uwnddata[l][k][j][i]  )*
						(uwnddata[l][k][j][i+1]-uwnddata[l][k][j][i]  )/dxs[j]-
						(uwnddata[l][k][j][i]  +uwnddata[l][k][j][i-1])*
						(uwnddata[l][k][j][i]  -uwnddata[l][k][j][i-1])/dxs[j]
					)/(dxs[j]*2);
					/**/
					Fdata[l][k][j][i]+=-(
						vwnddata[l][k][j][i+1]*(uwnddata[l][k][j+1][i+1]-uwnddata[l][k][j-1][i+1])/dy-
						vwnddata[l][k][j][i-1]*(uwnddata[l][k][j+1][i-1]-uwnddata[l][k][j-1][i-1])/dy
					)/dxs[j];
					/**/
					Fdata[l][k][j][i]+=-(
						 uwnddata[l][k][j+1][i]*(vwnddata[l][k][j+1][i+1]-vwnddata[l][k][j+1][i-1])-
						 uwnddata[l][k][j-1][i]*(vwnddata[l][k][j-1][i+1]-vwnddata[l][k][j-1][i-1])
					)/(dy*dxs[(y+1)/2]*lcos[j]*4);
					/**/
					Fdata[l][k][j][i]+=-(
						(vwnddata[l][k][j+1][i]+vwnddata[l][k][j][i])*(float)cos((ydef[j+1]+ydef[j])/2)*
						(vwnddata[l][k][j+1][i]-vwnddata[l][k][j][i])/dy-
						(vwnddata[l][k][j][i]+vwnddata[l][k][j-1][i])*(float)cos((ydef[j-1]+ydef[j])/2)*
						(vwnddata[l][k][j][i]-vwnddata[l][k][j-1][i])/dy
					)/(dy*2*lcos[j]);
					
					Fdata[l][k][j][i]+=f1[j]*vordata[l][k][j][i];
					
					Fdata[l][k][j][i]+=-uwnddata[l][k][j][i]*f2[j]/EARTH_RADIUS;
					/**/
					Fdata[l][k][j][i]+=(
						uwnddata[l][k][j][i+1]*vwnddata[l][k][j][i+1]*ltan[j]-
						uwnddata[l][k][j][i-1]*vwnddata[l][k][j][i-1]*ltan[j]
					)/(dxs[j]*2)/EARTH_RADIUS;
					/**/
					Fdata[l][k][j][i]+=-(
						uwnddata[l][k][j+1][i]*uwnddata[l][k][j+1][i]*lsin[j]-
						uwnddata[l][k][j-1][i]*uwnddata[l][k][j-1][i]*lsin[j]
					)/(dy*2)/lcos[j]/EARTH_RADIUS;
					/**/
					Fdata[l][k][j][i]+=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
						airdata[l][k][j][i]*omegadata[l][k][j][i]*divdata[l][k][j][i];
					/**/
					Fdata[l][k][j][i]+=(
						(omegadata[l][k][j][i+1]*airdata[l][k][j][i+1]-omegadata[l][k][j][i-1]*airdata[l][k][j][i-1])*
						uwnddata[l][k][j][i]/(dxs[j]*2)+
						(omegadata[l][k][j+1][i]*airdata[l][k][j+1][i]-omegadata[l][k][j-1][i]*airdata[l][k][j-1][i])*
						vwnddata[l][k][j][i]/(dy*2)
					)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
					/**/
					Fdata[l][k][j][i]+=(
						(omegadata[l][k][j][i+1]*airdata[l][k][j][i+1]-omegadata[l][k][j][i-1]*airdata[l][k][j][i-1])*
						f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
					);
					/**/
					Fdata[l][k][j][i]+=-omegadata[l][k][j][i]*(divdata[l][k+1][j][i]-divdata[l][k-1][j][i])/(dz+dz);
					/**/
					Fdata[l][k][j][i]+=-
						(omegadata[l][k][j][i+1]-omegadata[l][k][j][i-1])/(dxs[j]*2)*
						( uwnddata[l][k+1][j][i]- uwnddata[l][k-1][j][i])/(dz+dz);
					/**/
					Fdata[l][k][j][i]+=-
						(omegadata[l][k][j+1][i]-omegadata[l][k][j-1][i])/(dy*2)*
						( vwnddata[l][k+1][j][i]- vwnddata[l][k-1][j][i])/(dz+dz);
				}
				
				/************************ East boundary ************************/
				Fdata[l][k][j][0]=-(divdata[l+1][k][j][0]-divdata[l-1][k][j][0])/(2*dt);	//2*dt
				/**/
				Fdata[l][k][j][0]+=-(
					(uwnddata[l][k][j][1]+uwnddata[l][k][j][0]  )*
					(uwnddata[l][k][j][1]-uwnddata[l][k][j][0]  )/dxs[j]-
					(uwnddata[l][k][j][0]+uwnddata[l][k][j][x-1])*
					(uwnddata[l][k][j][0]-uwnddata[l][k][j][x-1])/dxs[j]
				)/(dxs[j]*2);
				/**/
				Fdata[l][k][j][0]+=-(
					vwnddata[l][k][j][1]*(uwnddata[l][k][j+1][1]-uwnddata[l][k][j-1][1])/dy-
					vwnddata[l][k][j][x-1]*(uwnddata[l][k][j+1][x-1]-uwnddata[l][k][j-1][x-1])/dy
				)/dxs[j];
				/**/
				Fdata[l][k][j][0]+=-(
					 uwnddata[l][k][j+1][0]*(vwnddata[l][k][j+1][1]-vwnddata[l][k][j+1][x-1])-
					 uwnddata[l][k][j-1][0]*(vwnddata[l][k][j-1][1]-vwnddata[l][k][j-1][x-1])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				/**/
				Fdata[l][k][j][0]+=-(
					(vwnddata[l][k][j+1][0]+vwnddata[l][k][j][0])*(float)cos((ydef[j+1]+ydef[j])/2)*
					(vwnddata[l][k][j+1][0]-vwnddata[l][k][j][0])/dy-
					(vwnddata[l][k][j][0]+vwnddata[l][k][j-1][0])*(float)cos((ydef[j-1]+ydef[j])/2)*
					(vwnddata[l][k][j][0]-vwnddata[l][k][j-1][0])/dy
				)/(dy*2*lcos[j]);
				
				Fdata[l][k][j][0]+=f1[j]*vordata[l][k][j][0];
				
				Fdata[l][k][j][0]+=-uwnddata[l][k][j][0]*f2[j]/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][0]+=(
					uwnddata[l][k][j][1]*vwnddata[l][k][j][1]*ltan[j]-
					uwnddata[l][k][j][x-1]*vwnddata[l][k][j][x-1]*ltan[j]
				)/(dxs[j]*2)/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][0]+=-(
					uwnddata[l][k][j+1][0]*uwnddata[l][k][j+1][0]*lsin[j]-
					uwnddata[l][k][j-1][0]*uwnddata[l][k][j-1][0]*lsin[j]
				)/(dy*2)/lcos[j]/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][0]+=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
					airdata[l][k][j][0]*omegadata[l][k][j][0]*divdata[l][k][j][0];
				/**/
				Fdata[l][k][j][0]+=(
					(omegadata[l][k][j][1]*airdata[l][k][j][1]-omegadata[l][k][j][x-1]*airdata[l][k][j][x-1])*
					uwnddata[l][k][j][0]/(dxs[j]*2)+
					(omegadata[l][k][j+1][0]*airdata[l][k][j+1][0]-omegadata[l][k][j-1][0]*airdata[l][k][j-1][0])*
					vwnddata[l][k][j][0]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][0]+=(
					(omegadata[l][k][j][1]*airdata[l][k][j][1]-omegadata[l][k][j][x-1]*airdata[l][k][j][x-1])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
				/**/
				Fdata[l][k][j][0]+=-omegadata[l][k][j][0]*(divdata[l][k+1][j][0]-divdata[l][k-1][j][0])/(dz+dz);
				/**/
				Fdata[l][k][j][0]+=-
					(omegadata[l][k][j][1]-omegadata[l][k][j][x-1])/(dxs[j]*2)*
					( uwnddata[l][k+1][j][0]- uwnddata[l][k-1][j][0])/(dz+dz);
				/**/
				Fdata[l][k][j][0]+=-
					(omegadata[l][k][j+1][0]-omegadata[l][k][j-1][0])/(dy*2)*
					( vwnddata[l][k+1][j][0]- vwnddata[l][k-1][j][0])/(dz+dz);
				
				/************************ West boundary ************************/
				Fdata[l][k][j][x-1]=-(divdata[l+1][k][j][x-1]-divdata[l-1][k][j][x-1])/(2*dt);	//2*dt
				/**/
				Fdata[l][k][j][x-1]+=-(
					(uwnddata[l][k][j][0]+uwnddata[l][k][j][x-1]  )*
					(uwnddata[l][k][j][0]-uwnddata[l][k][j][x-1]  )/dxs[j]-
					(uwnddata[l][k][j][x-1]  +uwnddata[l][k][j][x-2])*
					(uwnddata[l][k][j][x-1]  -uwnddata[l][k][j][x-2])/dxs[j]
				)/(dxs[j]*2);
				/**/
				Fdata[l][k][j][x-1]+=-(
					vwnddata[l][k][j][0]*(uwnddata[l][k][j+1][0]-uwnddata[l][k][j-1][0])/dy-
					vwnddata[l][k][j][x-2]*(uwnddata[l][k][j+1][x-2]-uwnddata[l][k][j-1][x-2])/dy
				)/dxs[j];
				/**/
				Fdata[l][k][j][x-1]+=-(
					 uwnddata[l][k][j+1][x-1]*(vwnddata[l][k][j+1][0]-vwnddata[l][k][j+1][x-2])-
					 uwnddata[l][k][j-1][x-1]*(vwnddata[l][k][j-1][0]-vwnddata[l][k][j-1][x-2])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				/**/
				Fdata[l][k][j][x-1]+=-(
					(vwnddata[l][k][j+1][x-1]+vwnddata[l][k][j][x-1])*(float)cos((ydef[j+1]+ydef[j])/2)*
					(vwnddata[l][k][j+1][x-1]-vwnddata[l][k][j][x-1])/dy-
					(vwnddata[l][k][j][x-1]+vwnddata[l][k][j-1][x-1])*(float)cos((ydef[j-1]+ydef[j])/2)*
					(vwnddata[l][k][j][x-1]-vwnddata[l][k][j-1][x-1])/dy
				)/(dy*2*lcos[j]);
				
				Fdata[l][k][j][x-1]+=f1[j]*vordata[l][k][j][x-1];
				
				Fdata[l][k][j][x-1]+=-uwnddata[l][k][j][x-1]*f2[j]/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][x-1]+=(
					uwnddata[l][k][j][0]*vwnddata[l][k][j][0]*ltan[j]-
					uwnddata[l][k][j][x-2]*vwnddata[l][k][j][x-2]*ltan[j]
				)/(dxs[j]*2)/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][x-1]+=-(
					uwnddata[l][k][j+1][x-1]*uwnddata[l][k][j+1][x-1]*lsin[j]-
					uwnddata[l][k][j-1][x-1]*uwnddata[l][k][j-1][x-1]*lsin[j]
				)/(dy*2)/lcos[j]/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][x-1]+=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
					airdata[l][k][j][x-1]*omegadata[l][k][j][x-1]*divdata[l][k][j][x-1];
				/**/
				Fdata[l][k][j][x-1]+=(
					(omegadata[l][k][j][0]*airdata[l][k][j][0]-omegadata[l][k][j][x-2]*airdata[l][k][j][x-2])*
					uwnddata[l][k][j][x-1]/(dxs[j]*2)+
					(omegadata[l][k][j+1][x-1]*airdata[l][k][j+1][x-1]-omegadata[l][k][j-1][x-1]*airdata[l][k][j-1][x-1])*
					vwnddata[l][k][j][x-1]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				/**/
				Fdata[l][k][j][x-1]+=(
					(omegadata[l][k][j][0]*airdata[l][k][j][0]-omegadata[l][k][j][x-2]*airdata[l][k][j][x-2])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
				/**/
				Fdata[l][k][j][x-1]+=-omegadata[l][k][j][x-1]*(divdata[l][k+1][j][x-1]-divdata[l][k-1][j][x-1])/(dz+dz);
				/**/
				Fdata[l][k][j][x-1]+=-
					(omegadata[l][k][j][0]-omegadata[l][k][j][x-2])/(dxs[j]*2)*
					( uwnddata[l][k+1][j][x-1]- uwnddata[l][k-1][j][x-1])/(dz+dz);
				/**/
				Fdata[l][k][j][x-1]+=-
					(omegadata[l][k][j+1][x-1]-omegadata[l][k][j-1][x-1])/(dy*2)*
					( vwnddata[l][k+1][j][x-1]- vwnddata[l][k-1][j][x-1])/(dz+dz);
			}
			
		}else{
			for(int k=1;k<z-1;k++)
			for(int l=1;l<t-1;l++)
			for(int j=1;j<y-1;j++){
				for(int i=1;i<x-1;i++){
					//
					Fdata[k][j][i][l]=-(divdata[k][j][i][l+1]-divdata[k][j][i][l-1])/(2*dt);	//2*dt
					/**/
					Fdata[k][j][i][l]+=-(
						(uwnddata[k][j][i+1][l]+uwnddata[k][j][i][l])*(uwnddata[k][j][i+1][l]-uwnddata[k][j][i][l])/dxs[j]-
						(uwnddata[k][j][i][l]+uwnddata[k][j][i-1][l])*(uwnddata[k][j][i][l]-uwnddata[k][j][i-1][l])/dxs[j]
					)/(dxs[j]*2);
					/**/
					Fdata[k][j][i][l]+=-(
						vwnddata[k][j][i+1][l]*(uwnddata[k][j+1][i+1][l]-uwnddata[k][j-1][i+1][l])/dy-
						vwnddata[k][j][i-1][l]*(uwnddata[k][j+1][i-1][l]-uwnddata[k][j-1][i-1][l])/dy
					)/dxs[j];
					/**/
					Fdata[k][j][i][l]+=-(
						 uwnddata[k][j+1][i][l]*(vwnddata[k][j+1][i+1][l]-vwnddata[k][j+1][i-1][l])-
						 uwnddata[k][j-1][i][l]*(vwnddata[k][j-1][i+1][l]-vwnddata[k][j-1][i-1][l])
					)/(dy*dxs[(y+1)/2]*lcos[j]*4);
					/**/
					Fdata[k][j][i][l]+=-(
						(vwnddata[k][j+1][i][l]+vwnddata[k][j][i][l])*(float)cos((ydef[j+1]+ydef[j])/2)*
						(vwnddata[k][j+1][i][l]-vwnddata[k][j][i][l])/dy-
						(vwnddata[k][j][i][l]+vwnddata[k][j-1][i][l])*(float)cos((ydef[j-1]+ydef[j])/2)*
						(vwnddata[k][j][i][l]-vwnddata[k][j-1][i][l])/dy
					)/(dy*2*lcos[j]);
					
					Fdata[k][j][i][l]+=f1[j]*vordata[k][j][i][l];
					
					Fdata[k][j][i][l]+=-uwnddata[k][j][i][l]*f2[j]/EARTH_RADIUS;
					/**/
					Fdata[k][j][i][l]+=(
						uwnddata[k][j][i+1][l]*vwnddata[k][j][i+1][l]*ltan[j]-
						uwnddata[k][j][i-1][l]*vwnddata[k][j][i-1][l]*ltan[j]
					)/(dxs[j]*2)/EARTH_RADIUS;
					/**/
					Fdata[k][j][i][l]+=-(
						uwnddata[k][j+1][i][l]*uwnddata[k][j+1][i][l]*lsin[j]-
						uwnddata[k][j-1][i][l]*uwnddata[k][j-1][i][l]*lsin[j]
					)/(dy*2)/lcos[j]/EARTH_RADIUS;
					/**/
					Fdata[k][j][i][l]+=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
						airdata[k][j][i][l]*omegadata[k][j][i][l]*divdata[k][j][i][l];
					/**/
					Fdata[k][j][i][l]+=(
						(omegadata[k][j][i+1][l]*airdata[k][j][i+1][l]-omegadata[k][j][i-1][l]*airdata[k][j][i-1][l])*
						uwnddata[k][j][i][l]/(dxs[j]*2)+
						(omegadata[k][j+1][i][l]*airdata[k][j+1][i][l]-omegadata[k][j-1][i][l]*airdata[k][j-1][i][l])*
						vwnddata[k][j][i][l]/(dy*2)
					)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
					/**/
					Fdata[k][j][i][l]+=(
						(omegadata[k][j][i+1][l]*airdata[k][j][i+1][l]-omegadata[k][j][i-1][l]*airdata[k][j][i-1][l])*
						f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
					);
					/**/
					Fdata[k][j][i][l]+=-omegadata[k][j][i][l]*(divdata[k+1][j][i][l]-divdata[k-1][j][i][l])/(dz+dz);
					/**/
					Fdata[k][j][i][l]+=-
						(omegadata[k][j][i+1][l]-omegadata[k][j][i-1][l])/(dxs[j]*2)*
						( uwnddata[k+1][j][i][l]- uwnddata[k-1][j][i][l])/(dz+dz);
					/**/
					Fdata[k][j][i][l]+=-
						(omegadata[k][j+1][i][l]-omegadata[k][j-1][i][l])/(dy*2)*
						( vwnddata[k+1][j][i][l]- vwnddata[k-1][j][i][l])/(dz+dz);
				}
				
				/************************ East boundary ************************/
				Fdata[k][j][0][l]=-(divdata[k][j][0][l+1]-divdata[k][j][0][l-1])/(2*dt);	//2*dt
				/**/
				Fdata[k][j][0][l]+=-(
					(uwnddata[k][j][1][l]+uwnddata[k][j][0][l]  )*
					(uwnddata[k][j][1][l]-uwnddata[k][j][0][l]  )/dxs[j]-
					(uwnddata[k][j][0][l]+uwnddata[k][j][x-1][l])*
					(uwnddata[k][j][0][l]-uwnddata[k][j][x-1][l])/dxs[j]
				)/(dxs[j]*2);
				/**/
				Fdata[k][j][0][l]+=-(
					vwnddata[k][j][1][l]*(uwnddata[k][j+1][1][l]-uwnddata[k][j-1][1][l])/dy-
					vwnddata[k][j][x-1][l]*(uwnddata[k][j+1][x-1][l]-uwnddata[k][j-1][x-1][l])/dy
				)/dxs[j];
				/**/
				Fdata[k][j][0][l]+=-(
					 uwnddata[k][j+1][0][l]*(vwnddata[k][j+1][1][l]-vwnddata[k][j+1][x-1][l])-
					 uwnddata[k][j-1][0][l]*(vwnddata[k][j-1][1][l]-vwnddata[k][j-1][x-1][l])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				/**/
				Fdata[k][j][0][l]+=-(
					(vwnddata[k][j+1][0][l]+vwnddata[k][j][0][l])*(float)cos((ydef[j+1]+ydef[j])/2)*
					(vwnddata[k][j+1][0][l]-vwnddata[k][j][0][l])/dy-
					(vwnddata[k][j][0][l]+vwnddata[k][j-1][0][l])*(float)cos((ydef[j-1]+ydef[j])/2)*
					(vwnddata[k][j][0][l]-vwnddata[k][j-1][0][l])/dy
				)/(dy*2*lcos[j]);
				
				Fdata[k][j][0][l]+=f1[j]*vordata[k][j][0][l];
				
				Fdata[k][j][0][l]+=-uwnddata[k][j][0][l]*f2[j]/EARTH_RADIUS;
				/**/
				Fdata[k][j][0][l]+=(
					uwnddata[k][j][1][l]*vwnddata[k][j][1][l]*ltan[j]-
					uwnddata[k][j][x-1][l]*vwnddata[k][j][x-1][l]*ltan[j]
				)/(dxs[j]*2)/EARTH_RADIUS;
				/**/
				Fdata[k][j][0][l]+=-(
					uwnddata[k][j+1][0][l]*uwnddata[k][j+1][0][l]*lsin[j]-
					uwnddata[k][j-1][0][l]*uwnddata[k][j-1][0][l]*lsin[j]
				)/(dy*2)/lcos[j]/EARTH_RADIUS;
				/**/
				Fdata[k][j][0][l]+=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
					airdata[k][j][0][l]*omegadata[k][j][0][l]*divdata[k][j][0][l];
				/**/
				Fdata[k][j][0][l]+=(
					(omegadata[k][j][1][l]*airdata[k][j][1][l]-omegadata[k][j][x-1][l]*airdata[k][j][x-1][l])*
					uwnddata[k][j][0][l]/(dxs[j]*2)+
					(omegadata[k][j+1][0][l]*airdata[k][j+1][0][l]-omegadata[k][j-1][0][l]*airdata[k][j-1][0][l])*
					vwnddata[k][j][0][l]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				/**/
				Fdata[k][j][0][l]+=(
					(omegadata[k][j][1][l]*airdata[k][j][1][l]-omegadata[k][j][x-1][l]*airdata[k][j][x-1][l])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
				/**/
				Fdata[k][j][0][l]+=-omegadata[k][j][0][l]*(divdata[k+1][j][0][l]-divdata[k-1][j][0][l])/(dz+dz);
				/**/
				Fdata[k][j][0][l]+=-
					(omegadata[k][j][1][l]-omegadata[k][j][x-1][l])/(dxs[j]*2)*
					(uwnddata[k+1][j][0][l]-uwnddata[k-1][j][0][l])/(dz+dz);
				/**/
				Fdata[k][j][0][l]+=-
					(omegadata[k][j+1][0][l]-omegadata[k][j-1][0][l])/(dy*2)*
					( vwnddata[k+1][j][0][l]- vwnddata[k-1][j][0][l])/(dz+dz);
				
				/************************ West boundary ************************/
				Fdata[k][j][x-1][l]=-(divdata[k][j][x-1][l+1]-divdata[k][j][x-1][l-1])/(2*dt);	//2*dt
				/**/
				Fdata[k][j][x-1][l]+=-(
					(uwnddata[k][j][0][l]+uwnddata[k][j][x-1][l]  )*
					(uwnddata[k][j][0][l]-uwnddata[k][j][x-1][l]  )/dxs[j]-
					(uwnddata[k][j][x-1][l]+uwnddata[k][j][x-2][l])*
					(uwnddata[k][j][x-1][l]-uwnddata[k][j][x-2][l])/dxs[j]
				)/(dxs[j]*2);
				/**/
				Fdata[k][j][x-1][l]+=-(
					vwnddata[k][j][0][l]*(uwnddata[k][j+1][0][l]-uwnddata[k][j-1][0][l])/dy-
					vwnddata[k][j][x-2][l]*(uwnddata[k][j+1][x-2][l]-uwnddata[k][j-1][x-2][l])/dy
				)/dxs[j];
				/**/
				Fdata[k][j][x-1][l]+=-(
					 uwnddata[k][j+1][x-1][l]*(vwnddata[k][j+1][0][l]-vwnddata[k][j+1][x-2][l])-
					 uwnddata[k][j-1][x-1][l]*(vwnddata[k][j-1][0][l]-vwnddata[k][j-1][x-2][l])
				)/(dy*dxs[(y+1)/2]*lcos[j]*4);
				/**/
				Fdata[k][j][x-1][l]+=-(
					(vwnddata[k][j+1][x-1][l]+vwnddata[k][j][x-1][l])*(float)cos((ydef[j+1]+ydef[j])/2)*
					(vwnddata[k][j+1][x-1][l]-vwnddata[k][j][x-1][l])/dy-
					(vwnddata[k][j][x-1][l]+vwnddata[k][j-1][x-1][l])*(float)cos((ydef[j-1]+ydef[j])/2)*
					(vwnddata[k][j][x-1][l]-vwnddata[k][j-1][x-1][l])/dy
				)/(dy*2*lcos[j]);
				
				Fdata[k][j][x-1][l]+=f1[j]*vordata[k][j][x-1][l];
				
				Fdata[k][j][x-1][l]+=-uwnddata[k][j][x-1][l]*f2[j]/EARTH_RADIUS;
				/**/
				Fdata[k][j][x-1][l]+=(
					uwnddata[k][j][0][l]*vwnddata[k][j][0][l]*ltan[j]-
					uwnddata[k][j][x-2][l]*vwnddata[k][j][x-2][l]*ltan[j]
				)/(dxs[j]*2)/EARTH_RADIUS;
				/**/
				Fdata[k][j][x-1][l]+=-(
					uwnddata[k][j+1][x-1][l]*uwnddata[k][j+1][x-1][l]*lsin[j]-
					uwnddata[k][j-1][x-1][l]*uwnddata[k][j-1][x-1][l]*lsin[j]
				)/(dy*2)/lcos[j]/EARTH_RADIUS;
				/**/
				Fdata[k][j][x-1][l]+=Rd/GRAVITY_ACCERLERATION/EARTH_RADIUS/zdef[zstart-1+k]*
					airdata[k][j][x-1][l]*omegadata[k][j][x-1][l]*divdata[k][j][x-1][l];
				/**/
				Fdata[k][j][x-1][l]+=(
					(omegadata[k][j][0][l]*airdata[k][j][0][l]-omegadata[k][j][x-2][l]*airdata[k][j][x-2][l])*
					uwnddata[k][j][x-1][l]/(dxs[j]*2)+
					(omegadata[k][j+1][x-1][l]*airdata[k][j+1][x-1][l]-omegadata[k][j-1][x-1][l]*airdata[k][j-1][x-1][l])*
					vwnddata[k][j][x-1][l]/(dy*2)
				)*Rd/zdef[zstart-1+k]/GRAVITY_ACCERLERATION/EARTH_RADIUS;
				/**/
				Fdata[k][j][x-1][l]+=(
					(omegadata[k][j][0][l]*airdata[k][j][0][l]-omegadata[k][j][x-2][l]*airdata[k][j][x-2][l])*
					f2[j]*Rd/GRAVITY_ACCERLERATION/zdef[zstart-1+k]/(dxs[j]*2)
				);
				/**/
				Fdata[k][j][x-1][l]+=-omegadata[k][j][x-1][l]*(divdata[k+1][j][x-1][l]-divdata[k-1][j][x-1][l])/(dz+dz);
				/**/
				Fdata[k][j][x-1][l]+=-
					(omegadata[k][j][0][l]-omegadata[k][j][x-2][l])/(dxs[j]*2)*
					(uwnddata[k+1][j][x-1][l]-uwnddata[k-1][j][x-1][l])/(dz+dz);
				/**/
				Fdata[k][j][x-1][l]+=-
					(omegadata[k][j+1][x-1][l]-omegadata[k][j-1][x-1][l])/(dy*2)*
					( vwnddata[k+1][j][x-1][l]- vwnddata[k-1][j][x-1][l])/(dz+dz);
			}
		}
	}
	
	
	/**
     * calculate balanced geopotential
     * 
     * @param	u		u-wind
     * @param	v		v-wind
     * 
     * @return	gpt		geopotential
     */
	public Variable cBalancedGeopotentialBySOR(Variable u,Variable v){
		GlobalDynamicMethodsInSC gdm=
		new GlobalDynamicMethodsInSC((SphericalSpatialModel)sm);
		
		Variable vor=gdm.c2DVorticity(u,v);
		Variable gpt=new Variable("h",vor);
		gpt.setValue(5000);	gpt.setInner(0);
		
		Variable F=cTerm2(u);
		F.plusEq(cTerm3(u,v));
		F.plusEq(cTerm4(u,v));
		F.plusEq(cTerm5(v));
		F.plusEq(cTerm9(vor));
		F.plusEq(cTerm10(u));
		F.plusEq(cTerm11(u,v));
		F.plusEq(cTerm12(u));
		
		//sor(0.00001f,5000,true,gpt,F);
		
		return gpt;
	}
	
	public Variable cBalancedGeopotentialBySH(Variable u,Variable v){
		GlobalDynamicMethodsInSC gdm=
		new GlobalDynamicMethodsInSC((SphericalSpatialModel)sm);
		
		SphericalHarmonicExpansion she=
		new SphericalHarmonicExpansion((SphericalSpatialModel)sm);
		
		she.setM(sm.getYCount()-1);
		
		Variable vor=gdm.c2DVorticity(u,v);
		
		Variable F=cTerm2(u);
		F.plusEq(cTerm3(u,v));
		F.plusEq(cTerm4(u,v));
		F.plusEq(cTerm5(v));
		F.plusEq(cTerm9(vor));
		F.plusEq(cTerm10(u));
		F.plusEq(cTerm11(u,v));
		F.plusEq(cTerm12(u));
		
		return she.solvePoissonEquation(F);
	}
	
	
	/**
     * implement the methods in EllipticalInterface
     */
	public Variable cAPrime(Variable v){
		ystart=v.getRange().getYRange()[0];
		
		t=v.getTCount();	z=v.getZCount();	y=v.getYCount();	x=v.getXCount();
		
		A=new Variable("Ap",v);	A.setCommentAndUnit("coefficient A' of elliptic equation");
		
		float[][][][] Adata=A.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=1f/lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=1f/lcos[ystart-1+j];
		}
		
		return A;
	}
	
	public Variable cBPrime(Variable v){
		B=new Variable("Bp",v);	B.setCommentAndUnit("coefficient B' of elliptic equation");
		return B;
	}
	
	public Variable cCPrime(Variable v){
		ystart=v.getRange().getYRange()[0];
		
		t=v.getTCount();	z=v.getZCount();	y=v.getYCount();	x=v.getXCount();
		
		C=new Variable("Cp",v);	A.setCommentAndUnit("coefficient C' of elliptic equation");
		
		float[][][][] Adata=A.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y-1;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=(float)cos((ydef[ystart-1+j]+ydef[ystart-1+j+1])/2f);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y-1;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=(float)cos((ydef[ystart-1+j]+ydef[ystart-1+j+1])/2f);
		}
		
		return C;
	}
	
	public Variable cA(Variable T){
		t=T.getTCount();	z=T.getZCount();	y=T.getYCount();	x=T.getXCount();
		
		A=new Variable("Aa",T);
		
		float[][][][] Adata=A.getData();
		
		if(T.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++) Adata[l][k][j][i]=lcos[j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=1;j<y-1;j++)
			for(int k=0;k<z;k++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=lcos[j];
		}
		
		return A;
	}
	
	public Variable cB(Variable T){ return new Variable("Bb",T);}
	
	public Variable cC(Variable uwnd){
		C=cA(uwnd);	C.setName("Cc");
		
		return C;
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			miniufo.descriptor.CtlDescriptor ctl=
			new miniufo.descriptor.CtlDescriptor(new File("D:/Data/DiagnosisVortex/Haima/Haima.ctl"));
			SphericalSpacialModel ssm=new SphericalSpacialModel(ctl);
			GlobalHGTInSC gh=new GlobalHGTInSC(ssm);gh.setThreadCount(1);
			miniufo.application.basic.GlobalDynamicMethodsInSC gdm=
			new miniufo.application.basic.GlobalDynamicMethodsInSC(ssm);
			
			miniufo.diagnosis.Range range=new miniufo.diagnosis.Range("z(35,37);t(1,3)",ctl);
			
			Variable h=new Variable("h",true,range);
			Variable u=new Variable("u",true,range);
			Variable v=new Variable("v",true,range);
			
			Variable h2 =(Variable)h.clone();	h2.setName("h2");
			Variable h3 =(Variable)h.clone();	h3.setName("h3");
			Variable h4 =(Variable)h.clone();	h4.setName("h4");
			Variable h5 =(Variable)h.clone();	h5.setName("h5");
			Variable h9 =(Variable)h.clone();	h9.setName("h9");
			Variable h10=(Variable)h.clone();	h10.setName("h10");
			Variable h12=(Variable)h.clone();	h12.setName("h12");
			
			miniufo.io.DataRead dr=miniufo.io.DataIOFactory.getDataRead(ctl);
			dr.readData(h,u,v);	dr.closeFile();
			
			h.multiplyEq(9.8f);
			
			float undef=h.getUndef();
			
			Variable hb=(Variable)h.clone();	hb.setName("hb");	hb.setInner(0);	hb.setUndef(undef);
			Variable hh=(Variable)h.clone();	hh.setName("hh");	hh.setInner(0);	hh.setUndef(undef);
			
			Variable vor=gdm.cVorticity(u,v);
			
			Variable f2 =gh.cTerm2(u);		f2.setUndef(undef);		h2.setUndef(undef);
			Variable f3 =gh.cTerm3(u,v);	f3.setUndef(undef);		h3.setUndef(undef);
			Variable f4 =gh.cTerm4(u,v);	f4.setUndef(undef);		h4.setUndef(undef);
			Variable f5 =gh.cTerm5(v);		f5.setUndef(undef);		h5.setUndef(undef);
			Variable f9 =gh.cTerm9(vor);	f9.setUndef(undef);		h9.setUndef(undef);
			Variable f10=gh.cTerm10(u);		f10.setUndef(undef);	h10.setUndef(undef);
			Variable f12=gh.cTerm12(u);		f12.setUndef(undef);	h12.setUndef(undef);
			
			Variable ff=f2.plus(f3).plusEq(f4).plusEq(f5).plusEq(f9).plusEq(f10).plusEq(f12);
			ff.setName("ff");
			
			//gh.sor(0.000001f,h2 ,f2 );	h2.divideEq(9.8f);
			//gh.sor(0.000001f,h3 ,f3 );	h3.divideEq(9.8f);
			//gh.sor(0.000001f,h4 ,f4 );	h4.divideEq(9.8f);
			//gh.sor(0.000001f,h5 ,f5 );	h5.divideEq(9.8f);
			//gh.sor(0.000001f,h9 ,f9 );	h9.divideEq(9.8f);
			//gh.sor(0.000001f,h10,f10);	h10.divideEq(9.8f);
			//gh.sor(0.000001f,h12,f12);	h12.divideEq(9.8f);
			//gh.sor(0.000001f,hh ,ff);		hh.divideEq(9.8f);
			//gh.sor(0.000001f,hb ,null);	hb.divideEq(9.8f);
			
			h.divideEq(9.8f);
			
			miniufo.io.CtlDataWriteStream cdws=
			new miniufo.io.CtlDataWriteStream("d:/hgt.dat");
			cdws.writeData(ctl,u,v,h,hh,hb,h2,h3,h4,h5,h9,h10,h12,ff,f2,f3,f4,f5,f9,f10,f12);
			cdws.closeFile();
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
