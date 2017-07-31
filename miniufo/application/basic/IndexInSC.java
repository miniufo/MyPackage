/**
 * @(#)IndexInSC.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.lagrangian.Typhoon;
import miniufo.application.EquationInSphericalCoordinate;
import miniufo.application.advanced.CoordinateTransformation;


/**
 * All kinds of climate index
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class IndexInSC extends EquationInSphericalCoordinate{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by spacial model in spheral coordinate
     */
	public IndexInSC(SphericalSpatialModel ssm){ super(ssm);}
	
	
	/**
     * calculate blocking index according to Pelly and Hoskins (2003)
     *
     * @param	hgt		geopotential height field
     * @param	xspan	grids span in x-direction, even number only
     * @param	yspan	grids span in y-direction, even number only
     *
     * @return	blocking index field
     */
	public Variable cBlockingIndexByPelly(Variable hgt,int xspan,int yspan){
		if(xspan<2||xspan%2!=0) throw new IllegalArgumentException("xspan should be even and larger than 1");
		if(yspan<2||yspan%2!=0) throw new IllegalArgumentException("yspan should be even and larger than 1");
		
		xspan/=2;	yspan/=2;	ystart=hgt.getRange().getYRange()[0];
		
		t=hgt.getTCount();	z=hgt.getZCount();	y=hgt.getYCount();	x=hgt.getXCount();
		
		float undef=hgt.getUndef();
		
		Variable bi=new Variable("bi",hgt);	bi.setValue(undef);
		bi.setCommentAndUnit("blocking index by Pelly and Hoskins (2003)");
		
		float[][][][] hgtdata=hgt.getData();	float[][][][] bidata=bi.getData();
		
		if(hgt.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=yspan;j<y-yspan;j++)
			for(int i=xspan;i<x-xspan;i++){
				bidata[l][k][j][i]=0;
				
				for(int jj=j-yspan;jj<=j;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[l][k][j][i]-=hgtdata[l][k][jj][ii];
				
				for(int jj=j;jj<=j+yspan;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[l][k][j][i]+=hgtdata[l][k][jj][ii];
				
				bidata[l][k][j][i]/=(yspan+1)*(xspan*2+1);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=yspan;j<y-yspan;j++)
			for(int i=xspan;i<x-xspan;i++){
				bidata[k][j][i][l]=0;
				
				for(int jj=j-yspan;jj<=j;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[k][j][i][l]-=hgtdata[k][jj][ii][l];
				
				for(int jj=j;jj<=j+yspan;jj++)
				for(int ii=i-xspan;ii<=i+xspan;ii++) bidata[k][j][i][l]+=hgtdata[k][jj][ii][l];
				
				bidata[k][j][i][l]/=(yspan+1)*(xspan*2+1);
			}
		}
		
		return bi;
	}
	
	
	/**
     * Calculate 2-dimensional horizontal index in
     * storm-relative framework (Molinari and Vollaro 1990, JAS)
     *
     * @param	dd		global descriptor
     * @param	rng		range for index
     * @param	yinc	increment in radial direction (degree)
     * @param	ycount	ring count in cylindrical coordinates
     * @param	xcount	grid count in each ring
     * @param	types	types of index: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	2-dimensional horizontal REFC index
     */
	public static Variable[] c2DHorizontalIndex(DataDescriptor dd,String rng,Typhoon tr,float yinc,int ycount,int xcount,String... types){
		if(types.length==0) throw new IllegalArgumentException("type required");
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		
		if(!ssm.isGlobal())
			throw new IllegalArgumentException("require a global model");
		
		Range r=new Range(rng,dd);
		Range rtmp=new Range(getGlobalBufferRange(rng),dd);
		
		Variable u=new Variable("u",rtmp);
		Variable v=new Variable("v",rtmp);
		Variable T=new Variable("T",rtmp);
		
		miniufo.io.DataRead dr=miniufo.io.DataIOFactory.getDataRead(dd);
		dr.setPrinting(false);
		dr.readData(u,v,T);	dr.closeFile();
		
		int tc=types.length;
		Variable[] idx=new Variable[tc];
		for(int m=0;m<tc;m++){
			idx[m]=new Variable(types[m],r);
			idx[m].setUndef(dd.getUndef("u"));
			idx[m].setCommentAndUnit("2D horizontal index ("+types[m]+")");
		}
		
		float[] cu=new float[idx[0].getTCount()];
		float[] cv=new float[idx[0].getTCount()];
		
		int tstart=tr.getTag(dd.getTDef().getSamples()[r.getTRange()[0]-1].getLongTime());
		for(int l=0,L=idx[0].getTCount();l<L;l++){
			cu[l]=tr.getZonalVelocity()[tstart+l];
			cv[l]=tr.getMeridionalVelocity()[tstart+l];
		}
		
		for(int ystart=r.getYRange()[0]-1,jj=ystart,J=r.getYRange()[1];jj<J;jj++)
		for(int xstart=r.getXRange()[0]-1,ii=xstart,I=r.getXRange()[1];ii<I;ii++){
			float[] olons=new float[idx[0].getTCount()];
			float[] olats=new float[idx[0].getTCount()];
			
			for(int l=0,L=idx[0].getTCount();l<L;l++){
				olons[l]=(float)Math.toDegrees(ssm.getXDef().getSamples()[ii]);
				olats[l]=(float)Math.toDegrees(ssm.getYDef().getSamples()[jj]);
			}
			
			CylindricalSpatialModel csm=new CylindricalSpatialModel(olons,olats,dd.getZDef().getSamples(),yinc,ycount,xcount);
			CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			ThermoDynamicMethodsInCC tdm=new ThermoDynamicMethodsInCC(csm);
			
			Variable[] re=ct.reprojectToCylindrical(ct.transToCylindrical(u),ct.transToCylindrical(v));
			Variable Vr=re[1].copy();	// full radial velocity, not storm-relative
			Variable Tc=ct.transToCylindrical(T);
			
			dm.cStormRelativeAziRadVelocity(cu,cv,re[0],re[1]);
			Variable utm=re[0].anomalizeX();
			Variable uta=re[0]; re[1].anomalizeX();
			Variable vra=re[1];
			
			Variable[] CIDX=new Variable[tc];
			
			// for 0.3-deg interval
			int str=9,end=18;	// 300-600 km average
			//int str=12,end=21;	// 400-700 km average
			//int str=15,end=24;	// 500-800 km average
			
			for(int m=0;m<tc;m++){
				if(types[m].equalsIgnoreCase("REFC")){
					CIDX[m]=dm.cREFC(uta,vra).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("PEFC")){
					Vr.anomalizeX(); // not storm-relative radial velocity
					CIDX[m]=dm.cPEFC(Vr).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("AEFC")){
					Variable gaz=dm.cAbsoluteAngularMomentumByJohnson(re[0]);
					gaz.anomalizeX();
					CIDX[m]=dm.cAEFC(gaz,vra).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("ISB")){
					CIDX[m]=dm.cMeanInertialStabilityByUT(utm).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("ETA")){
					CIDX[m]=dm.cMeanAbsoluteVorticity(utm).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("ULFI")){
					Vr.anomalizeX(); // not storm-relative radial velocity
					CIDX[m]=dm.cREFC(uta,vra).plusEq(dm.cPEFC(Vr)).divideEq(dm.cMeanAbsoluteVorticity(utm)).divideEq(86400f).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("htHFC")){
					Variable th=tdm.cPotentialTemperature(Tc);
					th.anomalizeX();
					CIDX[m]=dm.cEddyHeatHFC(th,vra).averageAlong(Dimension.Y,str,end);
					
				}else throw new IllegalArgumentException("invalid type for horizontal index: "+types[m]);
				
				float[][][][] edata=CIDX[m].getData();
				float[][][][] idata= idx[m].getData();
				
				if(idx[m].isTFirst()){
					for(int l=0,L=idx[m].getTCount();l<L;l++)
					for(int k=0,K=idx[m].getZCount();k<K;k++)
					idata[l][k][jj-ystart][ii-xstart]=edata[l][k][0][0];
					
				}else{
					for(int l=0,L=idx[m].getTCount();l<L;l++)
					for(int k=0,K=idx[m].getZCount();k<K;k++)
					idata[k][jj-ystart][ii-xstart][l]=edata[k][0][0][l];
				}
			}
			
			if(ii==xstart) System.out.print(".");
		}
		
		System.out.println();
		
		return idx;
	}
	
	/**
     * calculate 2-dimensional horizontal index
     * (Molinari and Vollaro 1990, JAS)
     *
     * @param	dd		global descriptor
     * @param	rng		range for index
     * @param	yinc	increment in radial direction (degree)
     * @param	ycount	ring count in cylindrical coordinates
     * @param	xcount	grid count in each ring
     * @param	types	types of index: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	2-dimensional horizontal PEFC index
     */
	public static Variable[] c2DHorizontalIndex(DataDescriptor dd,String rng,float yinc,int ycount,int xcount,String... types){
		return c2DHorizontalIndex(dd,rng,0,0,yinc,ycount,xcount,types);
	}
	
	/**
     * calculate 2-dimensional horizontal index
     * (Molinari and Vollaro 1990, JAS)
     *
     * @param	dd		global descriptor
     * @param	rng		range for index
     * @param	cu		zonal velocity of cylindrical coordinates
     * @param	cv		meridional velocity of cylindrical coordinates
     * @param	yinc	increment in radial direction (degree)
     * @param	ycount	ring count in cylindrical coordinates
     * @param	xcount	grid count in each ring
     * @param	types	types of index: REFC, PEFC, AEFC, EAMA, FFCT, FFBS, ISB, ETA, ULFI
     *
     * @return	2-dimensional horizontal PEFC index
     */
	public static Variable[] c2DHorizontalIndex
	(DataDescriptor dd,String rng,float cu,float cv,float yinc,int ycount,int xcount,String... types){
		if(types.length==0) throw new IllegalArgumentException("type required");
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		
		if(!ssm.isGlobal())
			throw new IllegalArgumentException("require a global model");
		
		Range r   =new Range(rng,dd);
		Range rtmp=new Range(getGlobalBufferRange(rng),dd);
		
		Variable u=new Variable("u",rtmp);
		Variable v=new Variable("v",rtmp);
		Variable T=new Variable("T",rtmp);
		
		miniufo.io.DataRead dr=miniufo.io.DataIOFactory.getDataRead(dd);
		dr.setPrinting(false);
		dr.readData(u,v,T);	dr.closeFile();
		
		int tc=types.length;
		Variable[] idx=new Variable[tc];
		for(int m=0;m<tc;m++){
			idx[m]=new Variable(types[m],r);
			idx[m].setUndef(dd.getUndef("u"));
			idx[m].setCommentAndUnit("2D horizontal index ("+types[m]+")");
		}
		
		float[] cus=new float[idx[0].getTCount()];
		float[] cvs=new float[idx[0].getTCount()];
		
		for(int l=0,L=idx[0].getTCount();l<L;l++){ cus[l]=cu; cvs[l]=cv;}
		
		for(int ystart=r.getYRange()[0]-1,jj=ystart,J=r.getYRange()[1];jj<J;jj++)
		for(int xstart=r.getXRange()[0]-1,ii=xstart,I=r.getXRange()[1];ii<I;ii++){
			float[] olons=new float[idx[0].getTCount()];
			float[] olats=new float[idx[0].getTCount()];
			
			for(int l=0,L=idx[0].getTCount();l<L;l++){
				olons[l]=(float)Math.toDegrees(ssm.getXDef().getSamples()[ii]);
				olats[l]=(float)Math.toDegrees(ssm.getYDef().getSamples()[jj]);
			}
			
			CylindricalSpatialModel csm=new CylindricalSpatialModel(olons,olats,dd.getZDef().getSamples(),yinc,ycount,xcount);
			CoordinateTransformation ct=new CoordinateTransformation(ssm,csm);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			ThermoDynamicMethodsInCC tdm=new ThermoDynamicMethodsInCC(csm);
			
			Variable[] re=ct.reprojectToCylindrical(ct.transToCylindrical(u),ct.transToCylindrical(v));
			Variable Vr=re[1].copy(); // full radial velocity (not storm-relative)
			Variable Tc=ct.transToCylindrical(T);
			
			dm.cStormRelativeAziRadVelocity(cus,cvs,re[0],re[1]);
			Variable utm=re[0].anomalizeX();
			Variable uta=re[0]; re[1].anomalizeX();
			Variable vra=re[1];
			
			Variable[] CIDX=new Variable[tc];
			
			// for 0.3-deg interval
			int str=9,end=18;		// 300-600 km average
			//int str=12,end=21;	// 400-700 km average
			//int str=15,end=24;	// 500-800 km average
			
			for(int m=0;m<tc;m++){
				if(types[m].equalsIgnoreCase("REFC")){
					CIDX[m]=dm.cREFC(uta,vra).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("PEFC")){
					Vr.anomalizeX(); // not storm-relative radial velocity
					CIDX[m]=dm.cPEFC(Vr).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("AEFC")){
					Variable gaz=dm.cAbsoluteAngularMomentumByJohnson(re[0]);
					gaz.anomalizeX();
					CIDX[m]=dm.cAEFC(gaz,vra).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("ISB")){
					CIDX[m]=dm.cMeanInertialStabilityByUT(utm).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("ETA")){
					CIDX[m]=dm.cMeanAbsoluteVorticity(utm).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("ULFI")){
					Vr.anomalizeX(); // not storm-relative radial velocity
					CIDX[m]=dm.cREFC(uta,vra).plusEq(dm.cPEFC(Vr)).divideEq(dm.cMeanAbsoluteVorticity(utm)).divideEq(86400f).averageAlong(Dimension.Y,str,end);
					
				}else if(types[m].equalsIgnoreCase("htHFC")){
					Variable th=tdm.cPotentialTemperature(Tc);
					th.anomalizeX();
					CIDX[m]=dm.cEddyHeatHFC(th,vra).averageAlong(Dimension.Y,str,end);
					
				}else throw new IllegalArgumentException("invalid type for horizontal index");
				
				float[][][][] edata=CIDX[m].getData();
				float[][][][] idata= idx[m].getData();
				
				if(idx[0].isTFirst()){
					for(int l=0,L=idx[0].getTCount();l<L;l++)
					for(int k=0,K=idx[0].getZCount();k<K;k++)
					idata[l][k][jj-ystart][ii-xstart]=edata[l][k][0][0];
					
				}else{
					for(int l=0,L=u.getTCount();l<L;l++)
					for(int k=0,K=u.getZCount();k<K;k++)
					idata[k][jj-ystart][ii-xstart][l]=edata[k][0][0][l];
				}
			}
			
			if(ii==xstart) System.out.print(".");
		}
		
		System.out.println();
		
		return idx;
	}
	
	
	/**
     * get the range without the lon and lat
     *
     * @return	the range String
     */
	private static String getGlobalBufferRange(String rng){
		String[] ss=rng.split(";");
		StringBuilder sb=new StringBuilder();
		
		for(String s:ss)
		if(s.startsWith("t")||s.startsWith("z")||s.startsWith("lev")){
			sb.append(s);	sb.append(";");
		}
		
		if(sb.length()!=0) sb.delete(sb.length()-1,sb.length());
		
		return sb.toString();
	}
	
	
	/** test
	public static void main(String[] arg){
		try{
			miniufo.diagnosis.DiagnosisFactory df=new miniufo.diagnosis.DiagnosisFactory(
				"E:/ERAInterim/200/2004.200.uv.nc"
			);
			DataDescriptor dd=df.getDataDescriptor();
			
			miniufo.io.DataWrite dw=miniufo.io.DataIOFactory.getDataWrite(
				dd,"E:/ERAInterim/200/PRI2004.dat"
			);
			
			int tt=1;
			while(tt<=1464){
				String r="t("+tt+","+(tt+23)+");lon(80,180);lat(6,60)";
				System.out.println(r);
				
				Variable[] vel=df.getVariables(new Range(r,dd),"u","v");
				
				Variable pri=IndexInSC.cPotentialReintensifyIndex(dd,r,0.2f,30,36);
				
				vel[0].setUndef(vel[1].getUndef());
				pri.setUndef(vel[1].getUndef());
				
				dw.writeData(pri);
				
				tt+=24;
			}
			
			dw.writeCtl(dd);	dw.closeFile();
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
