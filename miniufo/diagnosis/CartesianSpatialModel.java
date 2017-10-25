/**
 * @(#)CartesianSpatialModel.java	1.0 2017.06.01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.SpatialCoordinate;
import miniufo.descriptor.TemporalCoordinate;


/**
 * Spatial model in Cartesian coordinates
 *
 * @version 1.0, 2017.06.01
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CartesianSpatialModel extends SpatialModel{
	//
	private float beta;			// a constant = df/dy
	
	private float[] fCor=null;	// Coriolis parameter
	
	
	/**
     * constructor
     *
     * @param	pr	a given parser
     */
	public CartesianSpatialModel(DataDescriptor dd){ _CartesianSpatialModel(dd);}
	
	public CartesianSpatialModel(DataDescriptor dd,LevelType levType){
		this.levType=levType;
		_CartesianSpatialModel(dd);
	}
	
	
	/*** getor and setor ***/
	public float getBeta() { return beta;}
	
	public float[] getFCor(){ return fCor;}
	
	
	/**
     * used to print out
     */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		int tcount=tdef.length(),zcount=zdef.length();
		int ycount=ydef.length(),xcount=xdef.length();
		
		sb.append("\nCartesian Spatial Model:\n");
		
		sb.append("tdef:\n  from: ");
		sb.append(tdef.getSamples()[0]);		sb.append("\n  to  : ");
		sb.append(tdef.getSamples()[tcount-1]);	sb.append("\ttcount: ");	sb.append(tcount);
		
		sb.append("\nzdef:\n  from: ");
		sb.append(zdef.getSamples()[0]/100f);		sb.append("\n  to  : ");
		sb.append(zdef.getSamples()[zcount-1]/100f);sb.append("\tzcount: ");sb.append(zcount);
		
		sb.append("\nydef:\n  from: ");
		sb.append(ydef.getSamples()[0]);		sb.append("\n  to  : ");
		sb.append(ydef.getSamples()[ycount-1]);	sb.append("\tycount: ");	sb.append(ycount);
		
		sb.append("\nxdef:\n  from: ");
		sb.append(xdef.getSamples()[0]);		sb.append("\n  to  : ");
		sb.append(xdef.getSamples()[xcount-1]);	sb.append("\txcount: ");	sb.append(xcount);
		
		return sb.toString();
	}
	
	
	/*** helper methods ***/
	private void _CartesianSpatialModel(DataDescriptor dd){
		if(!(dd instanceof CtsDescriptor))
		throw new IllegalArgumentException("DataDescriptor should be CtsDescriptor");
		
		CtsDescriptor cd=(CtsDescriptor)dd;
		
		tdef=(TemporalCoordinate)(dd.getTDef().clone());
		zdef=( SpatialCoordinate)(dd.getZDef().clone()); int zcount=zdef.length();
		ydef=( SpatialCoordinate)(dd.getYDef().clone()); int ycount=ydef.length();
		xdef=( SpatialCoordinate)(dd.getXDef().clone());
		
		beta=cd.getBeta();
		dxs =new float[]{xdef.getIncrements()[0]};
		fCor=new float[ycount];
		
		/*** z-direction ***/
		if(levType==LevelType.PRESSURE)
		for(int k=0,K=zcount;k<K;k++) zdef.getSamples()[k]*=100;	// change to Pa
		
		/*** y-direction ***/
		for(int j=0;j<ycount;j++){
			float tmp=ydef.getSamples()[j];
			fCor[j]=cd.getF0()+beta*tmp;
		}
		
		if(levType==LevelType.PRESSURE) dz=dd.getDZDef()[0]*100;	// change to Pa
		else dz=dd.getDZDef()[0];
		
		dt=dd.getDTDef()[0];
		dy=ydef.getIncrements()[0];
		dx=xdef.getIncrements()[0];
		
		this.periodicX=dd.isPeriodicX();
		this.periodicY=dd.isPeriodicY();
	}
	
	
	/** test
	public static void main(String[] arg){
		CartesianSpatialModel csm=new CartesianSpatialModel(
			new miniufo.descriptor.CtsDescriptor(new File("d:/Data/MITgcm/barotropicDG/BetaCartRL/Stat.cts"))
		);
		
		System.out.println(csm);
	}*/
}
