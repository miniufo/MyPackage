/**
 * @(#)KeffInCTS.java	1.0 2017.07.02
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import miniufo.diagnosis.Variable;


/**
 * This class contains the contour-related algorithms in Cartesian coordinates.
 * All the variables returned are defined in contour coordinates.
 *
 * @version 1.0, 2017.07.02
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class KeffInCTS extends Keffective{
	
	/**
     * Constructor
     *
     * @param	cssm	contour-Cartesian spatial model
     */
	public KeffInCTS(ContourCartesianSpatialModel ccsm){ super(ccsm);}
	
	
	/**
     * Compute minimum length square.
     */
	public Variable cMinimumLengthSquare(){
		int C=cntrs[0][0].getContourNumber();
		int t=areas.getTCount();
		int z=areas.getZCount();
		
		Variable Lmin2=new Variable("Lmin2",areas);
		Lmin2.setCommentAndUnit("squared minimum length (m^2)");
		
		double rngX=csm.getXRange();
		
		if(areas.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] ldata=Lmin2.getData()[l][k][0];
				
				for(int c=0;c<C;c++) ldata[c]=(float)(rngX*rngX);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] ldata=Lmin2.getData()[k][0];
				
				for(int c=0;c<C;c++) ldata[c][l]=(float)(rngX*rngX);
			}
		}
		
		return Lmin2;
	}
	
	
	/*** helper methods ***/
	protected double ydefTransform(double ydef){ return ydef;}
	
	
	/*** test **
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/ERAInterim/Keff/PV/GRDSqr.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] vs=df.getVariables(new Range("t(1,1)",dd),"pv","grd2pv");
		
		Variable pv =vs[0].multiplyEq(1e-6f);
		Variable grd=vs[1];
		
		ContourCoords cc=new ContourCoords(dd);
		
		Contours[][] ctss=cc.getContourInfo(pv,51);System.out.println(ctss[0][0]);
		
		Variable area=cc.integrateOverContour(null,pv,ctss);
		Variable grd2=cc.integrateOverContour(grd,pv,ctss);
		Variable grdA=cc.cGradientWRTArea(area,ctss);
		Variable Le2 =cc.cEquivalentLengthSquare(grd2,area,ctss);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"D:/Data/ERAInterim/Keff/PV/re.dat");
		dw.writeData(dd,area,grd2,grdA,Le2); dw.closeFile();
	}*/
}
