/**
 * @(#)BalanceInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application;

import miniufo.application.EquationInSphericalCoordinate;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;


/**
 * balance of equation
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class BalanceInSC extends EquationInSphericalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public BalanceInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isLinearModel()) throw new IllegalArgumentException("Wave equation requires a linear spacial model");
	}
	
	
	public Variable cLeft(Variable K,Variable u,Variable v,Variable w){
		if(!K.isLike(u)||!K.isLike(v)||!K.isLike(w)) throw new IllegalArgumentException("dimensions not same");
		
		ystart=K.getRange().getYRange()[0];
		t=K.getTCount();	y=K.getYCount();	x=K.getXCount();	z=K.getZCount();
		
		Variable A=new Variable("l",K);
		
		float[][][][] Adata=A.getData();
		float[][][][] Kdata=K.getData();
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		float[][][][] wdata=w.getData();
		
		if(K.isTFirst()){
			for(int l=1;l<t-1;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				Adata[l][k][j][i]=(Kdata[l+1][k][j][i]-Kdata[l-1][k][j][i])/(dt*2)+
				udata[l][k][j][i]*(Kdata[l][k][j][i+1]-Kdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
				vdata[l][k][j][i]*(Kdata[l][k][j+1][i]-Kdata[l][k][j-1][i])/(dy*2)+
				wdata[l][k][j][i]*(Kdata[l][k+1][j][i]-Kdata[l][k-1][j][i])/(dz*2);
			}
			
		}else{
			for(int l=1;l<t-1;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				Adata[k][j][i][l]=(Kdata[k][j][i][l+1]-Kdata[k][j][i][l-1])/(dt*2)+
				udata[k][j][i][l]*(Kdata[k][j][i+1][l]-Kdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
				vdata[k][j][i][l]*(Kdata[k][j+1][i][l]-Kdata[k][j-1][i][l])/(dy*2)+
				wdata[k][j][i][l]*(Kdata[k+1][j][i][l]-Kdata[k-1][j][i][l])/(dz*2);
			}
		}
		
		return A;
	}
	
	public Variable cRight(Variable h,Variable u,Variable v){
		if(!h.isLike(u)||!h.isLike(v)) throw new IllegalArgumentException("dimensions not same");
		
		ystart=h.getRange().getYRange()[0];
		t=h.getTCount();	y=h.getYCount();	x=h.getXCount();	z=h.getZCount();
		
		Variable A=new Variable("r",h);
		
		float[][][][] Adata=A.getData();
		float[][][][] hdata=h.getData();
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		
		if(h.isTFirst()){
			for(int l=1;l<t-1;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				Adata[l][k][j][i]=
				udata[l][k][j][i]*(hdata[l][k][j][i+1]-hdata[l][k][j][i-1])/(dxs[ystart-1+j]*2)+
				vdata[l][k][j][i]*(hdata[l][k][j+1][i]-hdata[l][k][j-1][i])/(dy*2);
			}
			
		}else{
			for(int l=1;l<t-1;l++)
			for(int k=1;k<z-1;k++)
			for(int j=1;j<y-1;j++)
			for(int i=1;i<x-1;i++){
				Adata[k][j][i][l]=
				udata[k][j][i][l]*(hdata[k][j][i+1][l]-hdata[k][j][i-1][l])/(dxs[ystart-1+j]*2)+
				vdata[k][j][i][l]*(hdata[k][j+1][i][l]-hdata[k][j-1][i][l])/(dy*2);
			}
		}
		
		return A;
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			miniufo.data.CtlDescriptor ctl=new miniufo.data.CtlDescriptor("d:\\data\\data.ctl");
			
			SpheralSpacialModel ssm=new SpheralSpacialModel(ctl);
			BalanceInSC blc=new BalanceInSC(ssm);
			
			miniufo.util.Range range=new miniufo.util.Range("",ctl);
			
			Variable h=new Variable("h",range);
			Variable u=new Variable("u",range);
			Variable v=new Variable("v",range);
			Variable w=new Variable("w",range);
			
			Variable K=(Variable)(u.clone());	K.setName("K");
			
			miniufo.io.CtlDataReadStream cdrs=new miniufo.io.CtlDataReadStream(ctl);
			
			cdrs.readData(h,u,v,w);	cdrs.closeFile();
			
			
			K.square().plus(((Variable)(v.clone())).square()).divide(2);
			
			Variable l=blc.cLeft(K,u,v,w);
			Variable r=blc.cRight(h,u,v);
			
			
			miniufo.io.CtlDataWriteStream cdws=new miniufo.io.CtlDataWriteStream("d:\\data\\lr.dat");
			
			cdws.writeData(ctl,l,r);
			
			cdws.closeFile();
	    	
	    }catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
	}*/
}
