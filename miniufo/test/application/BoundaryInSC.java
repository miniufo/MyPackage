/**
 * @(#)BoundaryInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application;

import miniufo.application.EquationInSphericalCoordinate;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import static miniufo.geophysics.atmos.ThermoDynamics.Cp;
import static miniufo.geophysics.atmos.ThermoDynamics.L0;


/**
 * boundary scheme in spheral coordinate
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class BoundaryInSC extends EquationInSphericalCoordinate{
	//
	private static final float Ce=1.5e-3f;
	private static final float Ch=1.5e-3f;
	
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public BoundaryInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/**
     * calculate the distance from pressure level to the ground, < 0 means under the surface
     *
     * @param	hgt		potential height
     * @param	hsfc	potential height of surface
     *
     * @return	plh		pressure level height
     *
     * @exception	if hgt,hsfc are dimensionlly not the same
     */
	public Variable cPressureLevelHeight(Variable hgt,Variable hsfc){
		if(!hgt.isAreaLike(hsfc)) throw new IllegalArgumentException("dimensions not same");
		
		t=hgt.getTCount();	z=hgt.getZCount();	y=hgt.getYCount();	x=hgt.getXCount();
		
		float undef=hgt.getUndef();
		Variable plh=new Variable("plh",hgt);
		
		float[][][][]  plhdata= plh.getData();
		float[][][][]  hgtdata= hgt.getData();
		float[][][][] hsfcdata=hsfc.getData();
		
		if(hgt.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(hgtdata[l][k][j][i]!=undef) plhdata[l][k][j][i]=hgtdata[l][k][j][i]-hsfcdata[l][0][j][i];
				else hgtdata[l][k][j][i]=undef;
				
				if(plhdata[l][k][j][i]<0) plhdata[l][k][j][i]=undef;
				else if(plhdata[l][k][j][i]==0) plhdata[l][k][j][i]=2;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				if(hgtdata[k][j][i][l]!=undef) plhdata[k][j][i][l]=hgtdata[k][j][i][l]-hsfcdata[0][j][i][l];
				else hgtdata[k][j][i][l]=undef;
				
				if(plhdata[k][j][i][l]<0) plhdata[k][j][i][l]=undef;
				else if(plhdata[k][j][i][l]==0) plhdata[k][j][i][l]=2;
			}
		}
		
		return plh;
	}
	
	
	/**
     * calculate sensible heat flux
     *
     * @param	wsfc	wind speed of surface
     * @param	Tsfc	temperature of surface
     * @param	Tair	temperature of air
     * @param	hsfc	height of surface
     * @param	hgt		geopotential height
     * @param	pblh	planetary boundary layer height
     *
     * @return	re	sensible heat flux
     *
     * @exception if wsfc Tsfc hsfc are not dimensionally the same or Tair hgt are not dimensionally the same
     */
	public Variable cSensibleHeatFlux(Variable wsfc,Variable Tsfc,Variable T2m,Variable pblh){
		if(wsfc.getZCount()!=1) throw new IllegalArgumentException("surface variable only one level");
		if(!wsfc.isLike(Tsfc)||!wsfc.isLike(pblh)) throw new IllegalArgumentException("dimensions not same");
		
		Variable re=new Variable("hfs",wsfc);
		
		t=wsfc.getTCount();	y=wsfc.getYCount();	x=wsfc.getXCount();
		
		float[][][][]   redata=  re.getData();
		float[][][][]  T2mdata= T2m.getData();
		float[][][][] wsfcdata=wsfc.getData();
		float[][][][] Tsfcdata=Tsfc.getData();
		float[][][][] pblhdata=pblh.getData();
		
		float undef=wsfc.getUndef();	re.setUndef(undef);
		
		if(wsfc.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				redata[l][0][j][i]=
					Cp*Ch*wsfcdata[l][0][j][i]*(Tsfcdata[l][0][j][i]-T2mdata[l][0][j][i])/pblhdata[l][0][j][i];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				redata[0][j][i][l]=
					Cp*Ch*wsfcdata[0][j][i][l]*(Tsfcdata[0][j][i][l]-T2mdata[0][j][i][l])/pblhdata[0][j][i][l];
			}
		}
		
		return re;
	}
	
	/**
     * calculate latent heat flux
     *
     * @param	wsfc	wind speed of surface
     * @param	qsfc	specific humidity of surface
     * @param	qair	specific humidity of air
     * @param	Tair	temperature of air
     * @param	hsfc	height of surface
     * @param	hgt		geopotential height
     * @param	pblh	planetary boundary layer height
     *
     * @return	re	latent heat flux
     *
     * @exception if wsfc qsfc hsfc are not dimensionally the same or qair Tair hgt are not dimensionally the same
     */
	public Variable cLatentHeatFlux(Variable wsfc,Variable q2m,Variable qair,Variable Tair,Variable hsfc,Variable pblh,Variable hgt){
		if(wsfc.getZCount()!=1) throw new IllegalArgumentException("surface variable only one level");
		if(!wsfc.isLike(q2m)||!wsfc.isLike(hsfc)) throw new IllegalArgumentException("dimensions not same");
		if(!wsfc.isAreaLike(qair)||!qair.isLike(Tair)||!qair.isLike(hgt)) throw new IllegalArgumentException("dimensions not same");
		
		Variable re=new Variable("hfl",wsfc);
		
		t=Tair.getTCount();	z=Tair.getZCount();	y=Tair.getYCount();	x=Tair.getXCount();
		
		float[][][][]   redata=  re.getData();
		float[][][][]  hgtdata= hgt.getData();
		float[][][][]  q2mdata= q2m.getData();
		float[][][][] wsfcdata=wsfc.getData();
		float[][][][] qairdata=qair.getData();
		float[][][][] Tairdata=Tair.getData();
		float[][][][] hsfcdata=hsfc.getData();
		float[][][][] pblhdata=pblh.getData();
		
		float undef=wsfc.getUndef(),tmp,dh;	re.setUndef(undef);
		
		if(wsfc.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				tmp=pblhdata[l][0][j][i]+hsfcdata[l][0][j][i];
				
				for(int k=1;k<z;k++){
					dh=hgtdata[l][k][j][i]-hgtdata[l][k-1][j][i];
					
					if(hgtdata[l][k][j][i]>tmp){
						if(tmp-hgtdata[l][k-1][j][i]<dh/2)
							redata[l][0][j][i]=Ce*(L0-2327f*(Tairdata[l][k-1][j][i]-273.15f))*wsfcdata[l][0][j][i]*
								(q2mdata[l][0][j][i]-qairdata[l][k-1][j][i])/pblhdata[l][0][j][i];
						else
							redata[l][0][j][i]=Ce*(L0-2327f*(Tairdata[l][k][j][i]-273.15f))*wsfcdata[l][0][j][i]*
								(q2mdata[l][0][j][i]-qairdata[l][k][j][i])/pblhdata[l][0][j][i];
						
						break;
					}
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				tmp=pblhdata[0][j][i][l]+hsfcdata[0][j][i][l];
				
				for(int k=1;k<z;k++){
					dh=hgtdata[k][j][i][l]-hgtdata[k-1][j][i][l];
					
					if(hgtdata[k][j][i][l]>tmp){
						if(tmp-hgtdata[k-1][j][i][l]<dh/2)
							redata[0][j][i][l]=Ce*(L0-2327f*(Tairdata[k-1][j][i][l]-273.15f))*wsfcdata[0][j][i][l]*
								(q2mdata[0][j][i][l]-qairdata[k-1][j][i][l])/pblhdata[0][j][i][l];
						else
							redata[0][j][i][l]=Ce*(L0-2327f*(Tairdata[k][j][i][l]-273.15f))*wsfcdata[0][j][i][l]*
								(q2mdata[0][j][i][l]-qairdata[k][j][i][l])/pblhdata[0][j][i][l];
						
						break;
					}
				}
			}
		}
		
		return re;
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			miniufo.util.CtlDescriptor cd=new miniufo.util.CtlDescriptor("D:\\data\\Typhoon\\chanchu\\chanchu3.ctl");
			SpheralSpacialModel ssm=new SpheralSpacialModel(cd);
			BoundaryInSC tm=new BoundaryInSC(ssm);
			
			miniufo.util.Range range1=new miniufo.util.Range(null,cd);
			miniufo.util.Range range2=new miniufo.util.Range("z(1,1)",cd);
			
			Variable    q=new Variable("q"   ,true,range1);
			Variable    T=new Variable("t"   ,true,range1);
			Variable  hgt=new Variable("hgt" ,true,range1);
			//Variable  t2m=new Variable("t2m" ,true,range2);
			Variable  q2m=new Variable("spfh2m" ,true,range2);
			Variable u10m=new Variable("u10m",true,range2);
			Variable v10m=new Variable("v10m",true,range2);
			//Variable Tsfc=new Variable("tsfc",true,range2);
			Variable hsfc=new Variable("hsfc",true,range2);
			Variable pblh=new Variable("pblh",true,range2);
			
			miniufo.io.CtlDataReadStream cdrs=new miniufo.io.CtlDataReadStream(cd);
			
			miniufo.io.CtlDataWriteStream cdws=new miniufo.io.CtlDataWriteStream("D:\\data\\Typhoon\\chanchu\\flux.dat");
			
			cdrs.readData(hgt,q,T,q2m,u10m,v10m,hsfc,pblh);
			
			u10m.square().add(v10m.square()).sqrt();
			
			//Variable hfs=tm.cSensibleHeatFlux(u10m,Tsfc,t2m,pblh);			
			Variable hfs=tm.cLatentHeatFlux(u10m,q2m,q,T,hsfc,pblh,hgt);
			
			cdws.writeData(cd,u10m,q2m,q,T,hsfc,pblh,hgt,hfs);
			
			cdrs.closeFile();	cdws.closeFile();
	    	
	    }catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
	}*/
}
