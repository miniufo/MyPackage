/**
 * @(#)Keffective.java	1.0 2017.07.02
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import miniufo.basic.InterpolationModel;
import miniufo.diagnosis.Variable;


/**
 * This class abstract the contour-related algorithms in contour coordinates.
 * All the variables returned are defined in contour coordinates.
 *
 * @version 1.0, 2017.07.02
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class Keffective{
	//
	protected Variable tracer=null;
	protected Variable areas =null;
	
	protected Contours[][] cntrs=null;
	
	protected ContourSpatialModel csm=null;
	
	
	/**
     * Constructor
     *
     * @param	cssm	contour spatial model
     */
	public Keffective(ContourSpatialModel csm){
		this.csm   =csm;
		this.areas =csm.getAreasBoundedByContour();
		this.tracer=csm.getTracer();
		this.cntrs =csm.getContours();
		
		if(cntrs==null)
		throw new IllegalArgumentException("ContourSpatialModel is not initialized by tracer");
	}
	
	
	/**
     * Compute gradient of a given variable with respect to area.
     * 
     * @param	v		a given variable
     * 
     * @return	grdA	dv/dA
     */
	public Variable cGradientWRTArea(Variable v){
		int C=cntrs[0][0].getContourNumber();
		int t=v.getTCount();
		int z=v.getZCount();
		
		if(v.getYCount()!=1||v.getXCount()!=C)
		throw new IllegalArgumentException(v.getName()+" is not in contour coordinates");
		
		Variable grdA=new Variable("grdA"+v.getName(),v);
		grdA.setCommentAndUnit("gradient of "+v.getName()+" wrt area ("+v.getUnit()+" m^-2)");
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[] areas=cntrs[k][l].getMappedAreas();
				 float[] vdata= v.getData()[l][k][0];
				 float[] gdata=grdA.getData()[l][k][0];
				 
				// lower end
				gdata[0]=(float)((vdata[1]-vdata[0])/(areas[1]-areas[0]));
				
				// internal
				for(int c=1;c<C-1;c++){
					gdata[c]=(float)InterpolationModel.linearInterpolation(
						(vdata[c  ]-vdata[c-1])/(areas[c  ]-areas[c-1]),
						(vdata[c+1]-vdata[c  ])/(areas[c+1]-areas[c  ]),
						(areas[c  ]-areas[c-1])/(areas[c+1]-areas[c-1])
					);
				}
				
				// upper end
				gdata[C-1]=(float)((vdata[C-1]-vdata[C-2])/(areas[C-1]-areas[C-2]));
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[]  areas=cntrs[k][l].getMappedAreas();
				float[][] vdata= v.getData()[k][0];
				float[][] ldata=grdA.getData()[k][0];
				
				// lower end
				ldata[0][l]=(float)((vdata[1][l]-vdata[0][l])/(areas[1]-areas[0]));
				
				// internal
				for(int c=1;c<C-1;c++){
					ldata[c][l]=(float)InterpolationModel.linearInterpolation(
						(vdata[c  ][l]-vdata[c-1][l])/(areas[c  ]-areas[c-1]),
						(vdata[c+1][l]-vdata[c  ][l])/(areas[c+1]-areas[c  ]),
						(areas[c  ]   -areas[c-1]   )/(areas[c+1]-areas[c-1])
					);
				}
				
				// upper end
				ldata[C-1][l]=(float)((vdata[C-1][l]-vdata[C-2][l])/(areas[C-1]-areas[C-2]));
			}
		}
		
		return grdA;
	}
	
	/**
     * Compute gradient of the tracer with respect to area.
     * 
     * return	grdA	dq/dA
     */
	public Variable cGradientWRTArea(){
		int C=cntrs[0][0].getContourNumber();
		int t=areas.getTCount();
		int z=areas.getZCount();
		
		Variable grdA=new Variable("grdA"+tracer.getName(),areas);
		grdA.setCommentAndUnit("gradient of "+tracer.getName()+" wrt area ("+tracer.getUnit()+" m^-2)");
		
		if(areas.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[] areas=cntrs[k][l].getMappedAreas();
				double[] conts=cntrs[k][l].getValues();
				 float[] gdata=grdA.getData()[l][k][0];
				 
				// lower end
				gdata[0]=(float)((conts[1]-conts[0])/(areas[1]-areas[0]));
				
				// internal
				for(int c=1;c<C-1;c++){
					gdata[c]=(float)InterpolationModel.linearInterpolation(
						(conts[c  ]-conts[c-1])/(areas[c  ]-areas[c-1]),
						(conts[c+1]-conts[c  ])/(areas[c+1]-areas[c  ]),
						(areas[c  ]-areas[c-1])/(areas[c+1]-areas[c-1])
					);
				}
				
				// upper end
				gdata[C-1]=(float)((conts[C-1]-conts[C-2])/(areas[C-1]-areas[C-2]));
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[]  areas=cntrs[k][l].getMappedAreas();
				double[]  conts=cntrs[k][l].getValues();
				float[][] ldata=grdA.getData()[k][0];
				
				// lower end
				ldata[0][l]=(float)(conts[1]-conts[0]/(areas[1]-areas[0]));
				
				// internal
				for(int c=1;c<C-1;c++){
					ldata[c][l]=(float)InterpolationModel.linearInterpolation(
						(conts[c  ]-conts[c-1])/(areas[c  ]-areas[c-1]),
						(conts[c+1]-conts[c  ])/(areas[c+1]-areas[c  ]),
						(areas[c  ]-areas[c-1])/(areas[c+1]-areas[c-1])
					);
				}
				
				// upper end
				ldata[C-1][l]=(float)(conts[C-1]-conts[C-2]/(areas[C-1]-areas[C-2]));
			}
		}
		
		return grdA;
	}
	
	
	/**
     * Compute equivalent length square:
     * Le^2 = aveAlgQ / (dq/dA)^2
     *
     * @param	aveAlgQ		mean squared gradient along contours
     */
	public Variable cEquivalentLengthSquare(Variable aveAlgQ,Variable dqdA){
		Variable Le2=aveAlgQ.divide(dqdA).divideEq(dqdA);
		Le2.setName("Le2");
		Le2.setCommentAndUnit("squared equivalent length (m^2)");
		
		return Le2;
	}
	
	/**
     * Compute minimum length square.
     */
	public abstract Variable cMinimumLengthSquare();
	
	/**
     * Compute gradient with respect to equivalent Y: dq/dye.
     */
	public Variable cDqDye(){
		int C=cntrs[0][0].getContourNumber();
		int t=tracer.getTCount();
		int z=tracer.getZCount();
		
		Variable dqdye=new Variable("dqdye",areas);
		dqdye.setCommentAndUnit("derivative of q wrt equivalent Y: dq/dye ("+tracer.getUnit()+" m^-1)");
		
		if(areas.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[] ldata=dqdye.getData()[l][k][0];
				
				double[] vals=cntrs[k][l].getValues();
				double[] ydef=cntrs[k][l].getMappedEquivalentYs().clone();
				
				for(int j=0,J=ydef.length;j<J;j++) ydef[j]=ydefTransform(ydef[j]);
				
				// forward difference
				ldata[0]=(float)((vals[1]-vals[0])/(ydef[1]-ydef[0]));
				
				// centered difference
				for(int c=1;c<C-1;c++)
				ldata[c]=(float)InterpolationModel.linearInterpolation(
					(vals[c  ]-vals[c-1])/(ydef[c  ]-ydef[c-1]),
					(vals[c+1]-vals[c  ])/(ydef[c+1]-ydef[c  ]),
					(ydef[c  ]-ydef[c-1])/(ydef[c+1]-ydef[c-1])
				);
				
				// backward difference
				ldata[C-1]=(float)((vals[C-1]-vals[C-2])/(ydef[C-1]-ydef[C-2]));
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] ldata=dqdye.getData()[k][0];
				
				double[] vals=cntrs[k][l].getValues();
				double[] ydef=cntrs[k][l].getMappedEquivalentYs().clone();
				
				for(int j=0,J=ydef.length;j<J;j++) ydef[j]=ydefTransform(ydef[j]);
				
				// forward difference
				ldata[0][l]=(float)((vals[1]-vals[0])/(ydef[1]-ydef[0]));
				
				// centered difference
				for(int c=1;c<C-1;c++)
				ldata[c][l]=(float)InterpolationModel.linearInterpolation(
					(vals[c  ]-vals[c-1])/(ydef[c  ]-ydef[c-1]),
					(vals[c+1]-vals[c  ])/(ydef[c+1]-ydef[c  ]),
					(ydef[c  ]-ydef[c-1])/(ydef[c+1]-ydef[c-1])
				);
				
				// backward difference
				ldata[C-1][l]=(float)((vals[C-1]-vals[C-2])/(ydef[C-1]-ydef[C-2]));
			}
		}
		
		return dqdye;
	}
	
	/**
     * Compute normalized Keff = Le2/Lmin2 = <|grad(q)|^2> / (dq/dye)^2
     *
     * @param	aveAlgQ		mean squared gradient along contours
     * @param	dqdye		gradient with respect to equivalent Y
     */
	public Variable cNormalizedKeff(Variable aveAlgQ,Variable dqdye){
		Variable nkeff=aveAlgQ.copy();
		
		nkeff.divideEq(dqdye).divideEq(dqdye);
		
		nkeff.setName("nkeff");
		nkeff.setCommentAndUnit("normalized effective diffusivity (1)");
		
		return nkeff;
	}
	
	
	/*** helper methods ***/
	protected abstract double ydefTransform(double ydef);
	
	
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
