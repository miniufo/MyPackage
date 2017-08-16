/**
 * @(#)ContourSpatialModel.java	1.0 2017.07.02
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.contour;

import java.util.Arrays;

import miniufo.application.GeoFluidApplication.BoundaryCondition;
import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.SpatialCoordinate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * This class contains the contour-related algorithms abstract for various coordinates.
 *
 * @version 1.0, 2017.07.02
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class ContourSpatialModel{
	//
	protected int resRatio=1;	// resolution ratio, 1 for original resolution, 2 for halved, ...
	protected double rngX =0;	// delta X =  east - west
	protected double rngY =0;	// delta Y = north - south
	protected double areaT=0;	// total area of the region in DataDescriptor (m^2)
	
	protected float[]  xdefS=null;	// xdef scaled by resRatio
	protected float[]  ydefS=null;	// ydef scaled by resRatio
	protected float[] dxdefS=null;	// delta-xdef scaled by resRatio
	protected float[] dydefS=null;	// delta-ydef scaled by resRatio
	
	protected Contours[][] cntrs=null;	// contours information
	
	protected Variable tracer   =null;	// a given tracer (X-Y coordinates)
	protected Variable grdxy2   =null;	// squares of tracer gradient (X-Y coordinates)
	protected Variable areas    =null;	// areas enclosed by each contour (contour coordinate)
	
	protected DataDescriptor dd =null;	// describing the grids of the tracer data
	
	protected BoundaryCondition BCy=BoundaryCondition.Fixed;	// BCy for y
	protected BoundaryCondition BCx=BoundaryCondition.Fixed;	// BCx for x
	
	
	/**
     * constructor
     *
     * @param	dd	a DataDescriptor (X-Y descriptor)
     */
	public ContourSpatialModel(DataDescriptor dd){
		this.dd=dd;
		
		if(dd.isPeriodicX()) BCx=BoundaryCondition.Periodic;
	}
	
	
	/**
     * Initialization of tracer contour coordinates.
     * Call this method first before calling other methods
     *
     * @param	trace		a given tracer variable
     * @param	numOfC		number of contours
     * @param	resRatio	resolution ratio (>=1) for on-the-fly interpolation
     * @param	adjCtr		adjust contours or not
     */
	public void initContourByTracer(Variable tracer,int numOfC,int resRatio,boolean adjCtr){
		setResolutionRatio(resRatio);
		
		this.tracer=tracer;
		this.cntrs =newContours(numOfC);
		this.grdxy2=cSquaredTracerGradient();
		
		// to compute the area within each contour and check
		integrateWithinContour1(adjCtr);
		
		this.areas=toAreaVariable();
	}
	
	public void initContourByTracer(Variable tracer,int numOfC){
		initContourByTracer(tracer,numOfC,1,false);
	}
	
	/**
     * Initialization of tracer contour coordinates.
     * Call this method first before calling other methods
     *
     * @param	trace		a given tracer variable
     * @param	csouth		south-most contour
     * @param	cnorth		north-most contour
     * @param	inc			increment of contour from south to north (<0 if decreasing from south to north)
     * @param	resRatio	resolution ratio (>=1) for on-the-fly interpolation
     * @param	adjCtr		adjust contours or not
     */
	public void initContourByTracer(Variable tracer,float csouth,float cnorth,float inc,int resRatio,boolean adjCtr){
		setResolutionRatio(resRatio);
		
		this.tracer=tracer;
		this.cntrs =newContours(csouth,cnorth,inc);
		this.grdxy2=cSquaredTracerGradient();
		
		// to compute the area within each contour and check
		integrateWithinContour1(adjCtr);
		
		this.areas=toAreaVariable();
	}
	
	public void initContourByTracer(Variable tracer,float csouth,float cnorth,float inc){
		initContourByTracer(tracer,csouth,cnorth,inc,1,false);
	}
	
	
	/**
     * Set resolution ratio so that the tracer field could be refined.
     *
     * @param	resRatio	resolution ratio (>=1)
     */
	public void setResolutionRatio(int resRatio){
		if(resRatio<1) throw new IllegalArgumentException("resRatio should be at least 1");
		
		this.resRatio=resRatio;
		
		int x=dd.getXCount(),nx=(x-1)*resRatio+1;
		int y=dd.getYCount(),ny=(y-1)*resRatio+1;
		
		Type xType=Type.LINEAR;
		Type yType=Type.LINEAR;
		
		if(BCx==BoundaryCondition.Periodic){ xType=Type.PERIODIC_LINEAR; nx=x*resRatio;}
		if(BCy==BoundaryCondition.Periodic){ yType=Type.PERIODIC_LINEAR; ny=y*resRatio;}
		
		SpatialCoordinate xres=dd.getXDef().resample(nx,xType);
		SpatialCoordinate yres=dd.getYDef().resample(ny,yType);
		
		 xdefS=xres.getSamples();
		 ydefS=yres.getSamples();
		dxdefS=xres.getIncrements();
		dydefS=yres.getIncrements();
	}
	
	/**
	 * Print out the equivalent Ys for ctl.
	 * 
	 * @param	ttag	index of time
	 * @param	ztag	index of z-level
	 */
	public void printEquivalentYs(int ttag,int ztag){
		Contours cts=cntrs[ztag][ttag];
		
		for(int c=0,C=cts.getContourNumber();c<C;c++)
		System.out.println(String.format("%12.6f %12.7f %12.5g",
			cts.getEquivalentYs()[c],cts.getAreas()[c]/areaT,cts.getValues()[c]
		));
	}
	
	
	/**
     * Interpolate the equivalent Y (Ye) to a regular
     * equal-spaced Y of a given resolution.
     * 
     * For grids outside the Ye range, undefined value is filled.
     *
     * @param	v		a variable needs to be interpolated
     * @param	ycount	number of Ys within the domain
     * @param	type	type of interpolation
     */
	public Variable interpolatedToYs(Variable v,int ycount,Type type){
		float[] ydef=dd.getYDef().getSamples();
		
		double interval=rngY/(ycount-1);
		double[] dy=new double[ycount]; // destination y grids
		
		for(int j=0;j<ycount;j++) dy[j]=ydef[0]+interval*j;
		
		if(dy[ycount-1]>ydef[ydef.length-1]) dy[ycount-1]=ydef[ydef.length-1];
		
		return interpolatedToYs(v,dy,type);
	}
	
	/**
     * Interpolate the equivalent Y (Ye) to a prescribed Ys.
     * 
     * For grids outside the Ye range, undefined value is filled.
     *
     * @param	v		a variable needs to be interpolated
     * @param	dy		destination Ys to be interpolated to
     * @param	type	type of interpolation
     */
	public Variable interpolatedToYs(Variable v,double[] dy,Type type){
		int ycount=dy.length;
		int t=v.getTCount(),z=v.getZCount();
		int C=v.getXCount(),y=v.getYCount();
		
		if(y    !=1) throw new IllegalArgumentException("not an area-coordinate variable");
		if(ycount<2) throw new IllegalArgumentException("not enough knots for interpolation");
		if(dy[0]>dy[ycount-1]) throw new IllegalArgumentException("dy should be increasing");
		
		float undef=v.getUndef();
		
		Variable nv=new Variable(v.getName(),v.isTFirst(),new Range(t,z,ycount,1));
		nv.setCommentAndUnit(v.getCommentAndUnit());
		nv.setUndef(undef);
		nv.setValue(undef);
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[] sy=cntrs[k][l].getEquivalentYs();	// source y grids
				double[] sf=new double[sy.length];			// source functional value from y grids sf = f(sy)
				
				for(int c=0;c<C;c++) sf[c]=v.getData()[l][k][0][c];
				
				double[][] re=getValid(sy,sf,undef);
				
				double[] df=InterpolationModel.interp1D(re[0],re[1],dy,type,undef);
				
				float[][] ndata=nv.getData()[l][k];
				
				for(int j=0;j<ycount;j++) ndata[j][0]=(float)df[j];
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[] sy=cntrs[k][l].getEquivalentYs();	// source y grids
				double[] sf   =new double[C];				// source functional value from y grids sf = f(sy)
				
				for(int c=0;c<C;c++) sf[c]=v.getData()[k][0][c][l];
				
				double[][] re=getValid(sy,sf,undef);
				
				double[] df=InterpolationModel.interp1D(re[0],re[1],dy,type,undef);
				
				float[][][] ndata=nv.getData()[k];
				
				for(int j=0;j<ycount;j++) ndata[j][0][l]=(float)df[j];
			}
		}
		
		nv.getRange().setTRange(v.getRange());
		nv.getRange().setZRange(v.getRange());
		
		return nv;
	}
	
	/**
     * Areal-integrator that performed over areas enclosed by specific contours:
     * re = integral(v dS) within each contour
     *
     * @param	v	a given integrand variable
     */
	public Variable integrateWithinContour(Variable v){
		if(!v.isLike(tracer))
		throw new IllegalArgumentException("dimensions not same for:\n"+v+"\n"+tracer);
		
		if(tracer.getRange().getYRange()[0]!=1||tracer.getRange().getXRange()[0]!=1)
		throw new IllegalArgumentException("using tracer over the entire domain");
		
		float undef=v.getUndef();
		
		int t=tracer.getTCount();
		int z=tracer.getZCount();
		int y=ydefS.length;
		int x=xdefS.length;
		
		Type xType=Type.LINEAR;
		Type yType=Type.LINEAR;
		
		if(BCx==BoundaryCondition.Periodic) xType=Type.PERIODIC_LINEAR;
		if(BCy==BoundaryCondition.Periodic) yType=Type.PERIODIC_LINEAR;
		
		Variable re=new Variable(v.getName()+"ct",v.isTFirst(),new Range(t,z,1,cntrs[0][0].getContourNumber()));
		re.setCommentAndUnit("integration of "+v.getName()+" within contours ("+v.getUnit()+" m^2)");
		re.setUndef(undef);
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int C=cntrs[k][l].getContourNumber();
				
				double[] cts =cntrs[k][l].getValues();
				double[] rbuf=new double[C];
				
				float[][] vbuff=     v.getData()[l][k];
				float[][] tbuff=tracer.getData()[l][k];
				float[][] vdata=InterpolationModel.interp2D(vbuff,x,y,xType,yType,undef);
				float[][] tdata=InterpolationModel.interp2D(tbuff,x,y,xType,yType,undef);
				float[] extreme=ArrayUtil.getExtrema(tdata);
				
				if(cntrs[k][l].increaseSToN()){
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++) if(tdata[j][i]!=undef&&vdata[j][i]!=undef){
						//double dS=computeDS(xstart-1+i,ystart-1+j);
						double dS=computeDS(i,j);
						
						if(tdata[j][i]>=cts[0]) rbuf[0]+=vdata[j][i]*dS;
						for(int c=1;c<C;c++) if(tdata[j][i]>cts[c]) rbuf[c]+=vdata[j][i]*dS;
					}
				}else{
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++) if(tdata[j][i]!=undef&&vdata[j][i]!=undef){
						//double dS=computeDS(xstart-1+i,ystart-1+j);
						double dS=computeDS(i,j);
						
						if(tdata[j][i]<=cts[0]) rbuf[0]+=vdata[j][i]*dS;
						for(int c=1;c<C;c++) if(tdata[j][i]<cts[c]) rbuf[c]+=vdata[j][i]*dS;
					}
				}
				
				checkDataAtEnds(rbuf,cts,undef,extreme);
				
				float[] rdata=re.getData()[l][k][0];
				for(int c=0;c<C;c++) rdata[c]=(float)(rbuf[c]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				int C=cntrs[k][l].getContourNumber();
				
				double[]  cts=cntrs[k][l].getValues();
				double[] rbuf=new double[C];
				
				float[][] vbuff=new float[tracer.getYCount()][tracer.getXCount()];
				float[][] tbuff=new float[tracer.getYCount()][tracer.getXCount()];
				
				for(int j=0,J=tracer.getYCount();j<J;j++)
				for(int i=0,I=tracer.getXCount();i<I;i++){
					vbuff[j][i]=     v.getData()[k][j][i][l];
					tbuff[j][i]=tracer.getData()[k][j][i][l];
				}
				
				float[][] vdata=InterpolationModel.interp2D(vbuff,x,y,xType,yType,undef);
				float[][] tdata=InterpolationModel.interp2D(tbuff,x,y,xType,yType,undef);
				float[] extreme=ArrayUtil.getExtrema(tdata);
				
				if(cntrs[k][l].increaseSToN()){
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++) if(tdata[j][i]!=undef&&vdata[j][i]!=undef){
						//double dS=computeDS(xstart-1+i,ystart-1+j);
						double dS=computeDS(i,j);
						
						if(tdata[j][i]>=cts[0]) rbuf[0]+=vdata[j][i]*dS;
						for(int c=1;c<C;c++) if(tdata[j][i]>cts[c]) rbuf[c]+=vdata[j][i]*dS;
					}
				}else{
					for(int j=0;j<y;j++)
					for(int i=0;i<x;i++) if(tdata[j][i]!=undef&&vdata[j][i]!=undef){
						//double dS=computeDS(xstart-1+i,ystart-1+j);
						double dS=computeDS(i,j);
						
						if(tdata[j][i]<=cts[0]) rbuf[0]+=vdata[j][i]*dS;
						for(int c=1;c<C;c++) if(tdata[j][i]<cts[c]) rbuf[c]+=vdata[j][i]*dS;
					}
				}
				
				checkDataAtEnds(rbuf,cts,undef,extreme);
				
				float[][] rdata=re.getData()[k][0];
				for(int c=0;c<C;c++) rdata[c][l]=(float)(rbuf[c]);
			}
		}
		
		return re;
	}
	
	
	/*** getor and setor ***/
	public double getXRange(){ return rngX;}
	
	public double getYRange(){ return rngY;}
	
	public double getTotalArea(){ return areaT;}
	
	public Variable getAreasBoundedByContour(){ return areas;}
	
	public Variable getTracer(){ return tracer;}
	
	public Variable getSquaredTracerGradient(){ return grdxy2;}
	
	public Contours[][] getContours(){ return cntrs;}
	
	
	/*** helper methods and classes ***/
	
	/**
	 * compute the area element dS for a grid
	 * 
	 * @param	idx		zonal index
	 * @param	idy		meridional index
	 */
	protected abstract double computeDS(int idx,int idy);
	
	
	/**
	 * Compute the equivalent Ys given areas enclosed by tracer contours.
	 * 
	 * @param	areas	a series of area (m^2)
	 * @param	undef	undefined value
	 */
	protected double[] cEquivalentYs(double[] areas,float undef){
		int L=areas.length;
		
		double[] Ys=new double[L];
		
		for(int l=0;l<L;l++) Ys[l]=cEquivalentY(areas[l],undef);
		
		return Ys;
	}
	
	/**
	 * Compute the equivalent Y given the area enclosed by a tracer contour.
	 * 
	 * @param	area	area (m^2)
	 * @param	undef	undefined value
	 */
	protected abstract double cEquivalentY(double area,float undef);
	
	/**
     * Areal-integrator that performed over areas enclosed by specific contours
     * with 1 as integrand: area = integral(1 dS) within each contour
     * 
     * @param	selfAdjust	whether to automatically adjust contour interval each time
     * 						if false, then adjust interval when non-monotonic area-contour occurs
     */
	protected void integrateWithinContour1(boolean selfAdjust){
		if(tracer.getRange().getYRange()[0]!=1||tracer.getRange().getXRange()[0]!=1)
		throw new IllegalArgumentException("using tracer over the entire domain");
		
		float undef=tracer.getUndef();
		
		Type xType=Type.LINEAR;
		Type yType=Type.LINEAR;
		
		if(BCx==BoundaryCondition.Periodic) xType=Type.PERIODIC_LINEAR;
		if(BCy==BoundaryCondition.Periodic) yType=Type.PERIODIC_LINEAR;
		
		if(tracer.isTFirst()){
			for(int l=0,L=tracer.getTCount();l<L;l++)
			for(int k=0,K=tracer.getZCount();k<K;k++){
				float[][] tbuff=tracer.getData()[l][k];
				float[][] tdata=InterpolationModel.interp2D(tbuff,xdefS.length,ydefS.length,xType,yType,undef);
				
				String info=mappingArea(cntrs[k][l],tdata);
				
				if(info!=null||selfAdjust){
					cntrs[k][l]=adjustContours(cntrs[k][l],undef);
					
					info=mappingArea(cntrs[k][l],tdata);
					
					if(info!=null) throw new IllegalArgumentException(
						"at ttag ("+l+") and ztag ("+k+") "+
						"still found non-monotonic contour-area relation after adjustment\n"+info
					);
				}
			}
			
		}else{
			for(int l=0,L=tracer.getTCount();l<L;l++)
			for(int k=0,K=tracer.getZCount();k<K;k++){
				float[][] tbuff=new float[tracer.getYCount()][tracer.getXCount()];
				
				for(int j=0,J=tracer.getYCount();j<J;j++)
				for(int i=0,I=tracer.getXCount();i<I;i++) tbuff[j][i]=tracer.getData()[k][j][i][l];
				
				float[][] tdata=InterpolationModel.interp2D(tbuff,xdefS.length,ydefS.length,xType,yType,undef);
				
				String info=mappingArea(cntrs[k][l],tdata);
				
				if(l==0&&k==0){
					System.out.println(cntrs[0][0]);
					System.out.println(Arrays.toString(cntrs[0][0].getValues()));
					System.out.println(Arrays.toString(cntrs[0][0].getAreas()));
				}
				
				if(info!=null||selfAdjust){
					System.out.println("adjusting");
					cntrs[k][l]=adjustContours(cntrs[k][l],undef);
					
					info=mappingArea(cntrs[k][l],tdata);
					
					if(info!=null) throw new IllegalArgumentException(
						"at ttag ("+l+") and ztag ("+k+") "+
						"still found non-monotonic contour-area relation after adjustment\n"+info
					);
				}
			}
		}
	}
	
	/**
     * Mapping each contour monotonically to a area.
     * 
     * @param	ct		contours
     * @param	tdata	2D x-y field after possible interpolation
     */
	private String mappingArea(Contours ct,float[][] tdata){
		int C=ct.getContourNumber();
		int x=xdefS.length;
		int y=ydefS.length;
		
		float undef=tracer.getUndef();
		
		double[] cts =ct.getValues();
		double[] rbuf=new double[C];
		
		float[] extreme=ArrayUtil.getExtrema(tdata,undef);
		
		if(ct.increaseSToN()){
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(tdata[j][i]!=undef){
				//double dS=computeDS(xstart-1+i,ystart-1+j);
				double dS=computeDS(i,j);
				
				if(tdata[j][i]>=cts[0]) rbuf[0]+=dS;
				for(int c=1;c<C;c++) if(tdata[j][i]>cts[c]) rbuf[c]+=dS;
			}
			
		}else{
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(tdata[j][i]!=undef){
				//double dS=computeDS(xstart-1+i,ystart-1+j);
				double dS=computeDS(i,j);
				
				if(tdata[j][i]<=cts[0]) rbuf[0]+=dS;
				for(int c=1;c<C;c++) if(tdata[j][i]<cts[c]) rbuf[c]+=dS;
			}
		}
		
		checkDataAtEnds(rbuf,cts,undef,extreme);
		String info=checkArea(rbuf,cts,undef,extreme);
		
		ct.setAreas(rbuf);
		ct.setYEs(cEquivalentYs(rbuf,undef));
		
		return info;
	}
	
	
	/**
	 * Get valid data for interpolation.
	 * 
	 * @param	sx		source x
	 * @param	sy		source y
	 * @param	undef	undefined value
	 */
	private double[][] getValid(double[] sx,double[] sy,float undef){
		if(sx.length!=sy.length)
		throw new IllegalArgumentException("lengths of sx and sy are not equal");
		
		int len=sx.length,elen=1;
		
		double last=undef;
		for(int i=0;i<len;i++) if(sx[i]!=undef){ last=sx[i]; break;}
		
		if(last==undef) throw new IllegalArgumentException("no valid data");
		
		for(int i=0;i<len;i++) if(sx[i]!=undef&&sx[i]!=last){ elen++; last=sx[i];}
		
		if(elen<1) throw new IllegalArgumentException("no valid data");
		
		double[][] re=new double[2][elen];
		
		last=undef;
		for(int i=0,ii=0;i<len;i++) if(sx[i]!=undef&&sx[i]!=last){
			re[0][ii]=sx[i];
			re[1][ii]=sy[i];
			last=sx[i];
			ii++;
		}
		
		return re;
	}
	
	/**
     * Ensure that the areas are monotonic.
     */
	private String checkArea(double[] areas,double[] values,float undef,float[] extreme){
		// get the first valid area (may close to the total area)
		int idx=0;
		for(int c=0,C=areas.length;c<C;c++) if(areas[c]!=undef){ idx=c; break;}
		
		if(areas[idx]>=areaT) areas[idx]=areaT;
		else{
			double dlat=cEquivalentY(areas[idx],undef)-dd.getYDef().getFirst();
			if(dlat/ArrayUtil.getMin(dd.getDYDef())<5e-4) areas[idx]=areaT;
		}
		
		for(int c=1,C=areas.length;c<C;c++) if(areas[c]!=undef&&areas[c]==areas[c-1])
		return 
			"areas of "+tracer.getName()+" ["+extreme[0]+", "+extreme[1]+"] "+
			"enclosed by contours ("+values[c-1]+","+values[c]+") are not monotonic ("+areas[c]+"=="+areas[c-1]+"),\n"+
			"please try reducing the number of contours ("+C+")";
		
		return null;
	}
	
	private void checkDataAtEnds(double[] data,double[] values,float undef,float[] extreme){
		double incre=(extreme[1]-extreme[0])/(values.length-1.0);
		
		for(int c=0,C=data.length;c<C;c++){
			if(values[c]<extreme[0]&&Math.abs(values[c]-extreme[0])/incre>2e-3) data[c]=undef;
			if(values[c]>extreme[1]&&Math.abs(values[c]-extreme[1])/incre>2e-3) data[c]=undef;
		}
	}
	
	
	/**
     * Get contours information from a given variable by specifying
     * number of contour, equally spaced between cmin(z,t) and cmax(z,t)
     */
	private Contours[][] newContours(int numOfC){
		float undef=tracer.getUndef();
		
		int t=tracer.getTCount(),z=tracer.getZCount();
		int y=tracer.getYCount(),x=tracer.getXCount();
		
		Contours[][] ctss=new Contours[z][t];
		
		if(tracer.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] vdata=tracer.getData()[l][k];
				
				ctss[k][l]=new Contours(vdata,undef,numOfC);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] bdata=new float[y][x];
				float[][][] vdata=tracer.getData()[k];
				
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) bdata[j][i]=vdata[j][i][l];
				
				ctss[k][l]=new Contours(bdata,undef,numOfC);
			}
		}
		
		return ctss;
	}
	
	/**
     * Get contours information from a given variable by specifying
     * the csouth and cnorth contours, as well as increment.
     */
	private Contours[][] newContours(float csouth,float cnorth,float inc){
		int t=tracer.getTCount();
		int z=tracer.getZCount();
		
		Contours[][] ctss=new Contours[z][t];
		
		for(int l=0;l<t;l++)
		for(int k=0;k<z;k++) ctss[k][l]=new Contours(csouth,cnorth,inc);
		
		return ctss;
	}
	
	/**
     * Get contours information from a given variable by specifying
     * number of contour, adjusting the interval so that the equivalent
     * Y is approximately equally-spaced.
     */
	private Contours adjustContours(Contours cs,float undef){
		double min=cs.getMinValue();
		double max=cs.getMaxValue();
		
		double[] values=cs.getValues();
		double[] areas =cs.getAreas();
		
		double[][] valid=getValid(values,areas,undef);
		
		int vLen=valid[0].length;
		
		boolean hasStr=false;
		boolean hasEnd=false;
		
		if(valid[1][0] ==areaT) hasStr=true;
		if(valid[1][vLen-1]==0) hasEnd=true;
		
		if(hasStr&&hasEnd){
			values=valid[0];
			areas =valid[1];
			
		}else if(hasStr&&!hasEnd){
			values=new double[vLen+1];
			areas =new double[vLen+1];
			
			System.arraycopy(valid[0],0,values,0,vLen);
			System.arraycopy(valid[1],0, areas,0,vLen);
			
			values[vLen]=max;
			 areas[vLen]=0;
			
		}else if(!hasStr&&hasEnd){
			values=new double[vLen+1];
			areas =new double[vLen+1];
			
			System.arraycopy(valid[0],0,values,1,vLen);
			System.arraycopy(valid[1],0, areas,1,vLen);
			
			values[0]=min;
			 areas[0]=areaT;
			
		}else{
			values=new double[vLen+2];
			areas =new double[vLen+2];
			
			System.arraycopy(valid[0],0,values,1,vLen);
			System.arraycopy(valid[1],0, areas,1,vLen);
			
			values[0]=min;   values[vLen]=max;
			 areas[0]=areaT;  areas[vLen]=0;
		}
		
		double[] newAs=new double[cs.getContourNumber()];
		
		double inc=areaT/(cs.getContourNumber()-1.0);
		
		for(int i=0,I=newAs.length-1;i<I;i++) newAs[i]=areaT-i*inc;
		newAs[newAs.length-1]=0; // ensure no roundoff error
		
		double[] newCs=InterpolationModel.interp1D(areas,values,newAs,Type.LINEAR);
		
		Contours ncs=new Contours(newCs);
		ncs.setAreas(newAs);
		ncs.setYEs(cEquivalentYs(newAs,undef));
		
		return ncs;
	}
	
	
	/**
     * Compute squared gradient in X-Y coordinates, i.e., |grad(tracer)|^2
     * 
     * @return	grd		squared gradient (a scalar)
     */
	protected abstract Variable cSquaredTracerGradient();
	
	/**
     * Collect the area data into a Variable in contour coordinate.
     */
	private Variable toAreaVariable(){
		float undef=tracer.getUndef();
		
		int t=tracer.getTCount();
		int z=tracer.getZCount();
		int C=cntrs[0][0].getContourNumber();
		
		Variable area=new Variable("area",tracer.isTFirst(),new Range(t,z,1,C));
		area.setCommentAndUnit("area within contours of "+tracer.getName()+" (m^2)");
		area.setUndef(undef);
		area.setValue(undef);
		
		if(tracer.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[] areas=cntrs[k][l].getAreas();
				 float[] adata=area.getData()[l][k][0];
				
				for(int c=0;c<C;c++) adata[c]=(float)(areas[c]);
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				double[] areas=cntrs[k][l].getAreas();
				float[][] adata=area.getData()[k][0];
				
				for(int c=0;c<C;c++) adata[c][l]=(float)(areas[c]);
			}
		}
		
		return area;
	}
	
	
	/*** test **
	public static void main(String[] args){
		float undef=6;
		double[] sx=new double[]{1,2,2,3,4,4,5,5,6};
		double[] sy=new double[]{3,4,4,5,6,6,7,7,8};
		
		double[][] re=getValid(sx,sy,undef);
		
		System.out.println(Arrays.toString(re[0]));
		System.out.println(Arrays.toString(re[1]));
	}*/
}