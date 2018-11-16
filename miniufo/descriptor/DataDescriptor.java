/**
 * @(#)DataDescriptor.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.MDate;
import miniufo.util.Region2D;
import miniufo.util.Region3D;


/**
 * Used to describe the 4D binary data
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class DataDescriptor{
	//
	protected boolean periodicX=false;
	protected boolean periodicY=false;
	
	protected TemporalCoordinate tdef=null;
	protected  SpatialCoordinate zdef=null;
	protected  SpatialCoordinate ydef=null;
	protected  SpatialCoordinate xdef=null;
	
	protected int vcount;
	
	protected boolean zrev   =false;
	protected boolean yrev   =false;
	protected boolean hasData=false;
	
	protected static boolean print  =true;
	
	protected float[] dtdef=null;	// unit: s, length is usually 1
	protected float[] dzdef=null;	// unit: hPa
	protected float[] dydef=null;	// unit: degree
	protected float[] dxdef=null;	// unit: degree
	
	protected  long[] times=null;	// MDate to long, length = tcount+1
	protected   Var[] vdef =null;
	
	protected String dsetPath  =null;	// path of dset
	protected String descPath  =null;	// path of descriptor
	protected String title     =null;	// title
	protected String tincrement=null;	// for writing ctl
	
	protected DataType dataType=DataType.OTHERS;
	
	protected enum IncreType{mn,hr,dy,mo,yr};
	
	public enum DataType{ANNUAL,SEASONAL,MONTHLY,DAILY,DAILY4,HOURLY,OTHERS}
	
	
	/*** getor and setor ***/
	public int getTNum(MDate md){ return getTNum(md.getLongTime());}
	
	public int getTNum(long tt){ return ArrayUtil.getIdxIncre(times,tt);}
	
	public int getZNum(float zpos){
		if(zdef.isIncre) return ArrayUtil.getIdxIncre(zdef.getSamples(),zpos);
		else return ArrayUtil.getIdxDecre(zdef.getSamples(),zpos);
	}
	
	public int getYNum(float ypos){ return ArrayUtil.getIdxIncre(ydef.getSamples(),ypos);}
	
	public int getXNum(float xpos){ return ArrayUtil.getIdxIncre(xdef.getSamples(),xpos);}
	
	public int getTLENum(MDate md){ return getTLENum(md.getLongTime());}
	
	public int getTLENum(long tt){ return ArrayUtil.getLEIdxIncre(times,tt);}
	
	public int getZLENum(float zpos){
		if(zdef.isIncre) return ArrayUtil.getLEIdxIncre(zdef.getSamples(),zpos);
		else return ArrayUtil.getLEIdxDecre(zdef.getSamples(),zpos);
	}
	
	public int getYLENum(float ypos){ return ArrayUtil.getLEIdxIncre(ydef.getSamples(),ypos);}
	
	public int getXLENum(float xpos){ return ArrayUtil.getLEIdxIncre(xdef.getSamples(),xpos);}
	
	public int getXNumPeriodicX(float xpos){
		float delx=dxdef[0];
		
		float lmin=xdef.getMin();
		float lmax=xdef.getMax();
		float lrep=lmax+delx;	// physically, lrep equals to lmin
		float wide=lrep-lmin;
		
		// set lon into [lmin, lrep)
		while(xpos>=lrep) xpos-=wide;
		while(xpos< lmin) xpos+=wide;
		
		
		if(xpos>lmax){
			// process the range (lmax, lrep)
			if(xpos<=lmax+delx/2f) return xdef.length()-1;
			else return 0;
			
		}else{
			// process the range [lmin, lmax]
			return getXNum(xpos);
		}
	}
	
	public int getXLENumPeriodicX(float xpos){
		float delx=dxdef[0];
		
		float lmin=xdef.getMin();
		float lmax=xdef.getMax();
		float lrep=lmax+delx;	// physically, lrep equals to lmin
		float wide=lrep-lmin;
		
		// set lon into [lmin, lrep)
		while(xpos>=lrep) xpos-=wide;
		while(xpos< lmin) xpos+=wide;
		
		if(xpos>lmax) return xdef.length()-1;	// process the range (lmax, lrep)
		else return getXLENum(xpos);				// process the range [lmin, lmax]
	}
	
	public int getTCount(){ return tdef.length();}
	
	public int getZCount(){ return zdef.length();}
	
	public int getYCount(){ return ydef.length();}
	
	public int getXCount(){ return xdef.length();}
	
	public int getVCount(){ return vcount ;}
	
	public int getVarZcount(String vname){
		for(int i=0,I=vdef.length;i<I;i++)
			if(vname.equals(vdef[i].getName())) return vdef[i].getZCount();
		
		throw new IllegalArgumentException("cannot find "+vname+" in "+descPath);
	}
	
	public boolean  isZRev(){ return zrev;}
	
	public boolean  isYRev(){ return yrev;}
	
	public boolean tLinear(){ return tdef.islinear;}
	
	public boolean zLinear(){ return zdef.islinear;}
	
	public boolean yLinear(){ return ydef.islinear;}
	
	public boolean xLinear(){ return xdef.islinear;}
	
	public boolean hasData(){ return hasData;}
	
	public boolean isPeriodicX(){ return periodicX;}
	
	public boolean isPeriodicY(){ return periodicY;}
	
	public  long[] getTimes(){return times;}
	
	public float[] getDTDef(){return dtdef;}
	
	public float[] getDZDef(){return dzdef;}
	
	public float[] getDYDef(){return dydef;}
	
	public float[] getDXDef(){return dxdef;}
	
	public   Var[] getVDef(){ return vdef ;}
	
	public String[] getVarNames(){
		int len=vdef.length;
		
		String[] names=new String[len];
		
		for(int i=0;i<len;i++) names[i]=vdef[i].vname;
		
		return names;
	}
	
	public String getTitle(){ return title;}
	
	public String  getDSet(){ return dsetPath;}
	
	public String  getPath(){ return descPath;}
	
	public String getTIncrement(){ return tincrement;}
	
	public String getVarUnit(String vname){
		for(int i=0;i<vdef.length;i++)
		if(vname.equalsIgnoreCase(vdef[i].getName())) return vdef[i].unit;
		
		throw new IllegalArgumentException("cannot find "+vname+" in "+descPath);
	}
	
	public String getVarComment(String vname){
		for(int i=0;i<vdef.length;i++)
		if(vname.equalsIgnoreCase(vdef[i].getName())) return vdef[i].comment;
		
		throw new IllegalArgumentException("cannot find "+vname+" in "+descPath);
	}
	
	public String getVarCommentAndUnit(String vname){
		for(int i=0;i<vdef.length;i++)
		if(vname.equalsIgnoreCase(vdef[i].getName())) return vdef[i].comment+" ("+vdef[i].unit+")";
		
		throw new IllegalArgumentException("cannot find "+vname+" in "+descPath);
	}
	
	public TemporalCoordinate getTDef(){ return tdef;}
	
	public  SpatialCoordinate getZDef(){ return zdef;}
	
	public  SpatialCoordinate getYDef(){ return ydef;}
	
	public  SpatialCoordinate getXDef(){ return xdef;}
	
	public DataType getDataType(){ return dataType;}
	
	public abstract float getUndef(String vname);
	
	public void setPeriodicX(boolean periodicX){ this.periodicX=periodicX;}
	
	public void setPeriodicY(boolean periodicY){ this.periodicY=periodicY;}
	
	
	/**
     * Post process the data.
     */
	protected void postProcess(){
		if(zdef.length()>1) dzdef=zdef.getIncrements();
		if(ydef.length()>1) dydef=ydef.getIncrements();
		if(xdef.length()>1) dxdef=xdef.getIncrements();
		
		int tcount=tdef.length();
		times=new long[tcount+1];	MDate[] tdf=tdef.getSamples();
		
		for(int i=0;i<tcount;i++) times[i]=tdf[i].getLongTime();
		times[tcount]=tdf[tcount-1].add(tincrement).getLongTime();
		
		     if("1yr".equalsIgnoreCase(tincrement)) dataType=DataType.ANNUAL  ;
		else if("3mo".equalsIgnoreCase(tincrement)) dataType=DataType.SEASONAL;
		else if("1mo".equalsIgnoreCase(tincrement)) dataType=DataType.MONTHLY ;
		else if("1dy".equalsIgnoreCase(tincrement)) dataType=DataType.DAILY   ;
		else if("6hr".equalsIgnoreCase(tincrement)) dataType=DataType.DAILY4  ;
		else if("1hr".equalsIgnoreCase(tincrement)) dataType=DataType.HOURLY  ;
	}
	
	/**
     * whether the x-dim is periodic given the total x-range
     */
	protected boolean isPeriodic(SpatialCoordinate sc,float range){
		// not physically but generally true
		if(!sc.islinear) return false;
		
		float delta=sc.getIncrements()[0];
		float start=sc.getLast()+delta-range;
		
		if(Math.abs((start-sc.getFirst())/delta)>1e-4) return false;
		
		return true;
	}
	
	/**
     * whether the dimension of the descriptor is the same with the given one
     */
	public boolean isLikes(DataDescriptor dd){
		if(tdef.length()!=dd.tdef.length()) return false;
		if(zdef.length()!=dd.zdef.length()) return false;
		if(ydef.length()!=dd.ydef.length()) return false;
		if(xdef.length()!=dd.xdef.length()) return false;
		
		return true;
	}
	
	/**
     * whether the area of the descriptor is the same with the given one
     */
	public boolean isAreaLike(DataDescriptor dd){
		if(tdef.length()!=dd.tdef.length()) return false;
		if(ydef.length()!=dd.ydef.length()) return false;
		if(xdef.length()!=dd.xdef.length()) return false;
		
		return true;
	}
	
	/**
     * Whether the given point is in the domain area.
     *
     * @param	lon		longitude (degree) of the given point
     * @param	lat		latitude  (degree) of the given point
     */
	public boolean isInArea(float xcoord,float ycoord){
		float x0=xdef.getFirst(),xN=xdef.getLast();
		float y0=ydef.getFirst(),yN=ydef.getLast();
		
		if(!periodicX&&(xcoord<x0||xcoord>xN)) return false;
		if(!periodicY&&(ycoord<y0||ycoord>yN)) return false;
		
		return true;
	}
	
	
	/**
     * used to print out
     */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		sb.append(  "  title:  ");	sb.append(title);
		sb.append("\n   path:  ");	sb.append(descPath);
		sb.append("\n   dset:  ");	sb.append(dsetPath);
		
		/*** print tdef ***/
		sb.append("\n\n   tdef:  ");
		
        sb.append("{\n\tStart time : ");
        sb.append(tdef.getSamples()[0]);
 		
 		sb.append("\n\t End  time : ");
 		sb.append(tdef.getSamples()[tdef.length()-1]);
 		
 		sb.append("\n\t time count: ");
 		sb.append(tdef.length());
 		
 		sb.append(" }");
		
		/*** print zdef ***/
		sb.append("\n\n   zdef:  ");
		
        sb.append("{ ");
        sb.append(zdef.getSamples()[0]);
 		
        for(int i=1,I=zdef.length();i<I;i++){
            sb.append(", ");
            sb.append(zdef.getSamples()[i]);
        }
 		
        sb.append(" }");
        
		/*** print ydef ***/
		sb.append("\n\n   ydef:  ");
		
        sb.append("{ ");
        sb.append(ydef.getSamples()[0]);
 		
        for(int i=1,I=ydef.length();i<I;i++){
            sb.append(", ");
            sb.append(ydef.getSamples()[i]);
        }
 		
        sb.append(" }");
        
		/*** print xdef ***/
		sb.append("\n\n   xdef:  ");
		
        sb.append("{ ");
        sb.append(xdef.getSamples()[0]);
 		
        for(int i=1,I=xdef.length();i<I;i++){
            sb.append(", ");
            sb.append(xdef.getSamples()[i]);
        }
 		
        sb.append(" }");
		
		sb.append("\n\n   vdef:  ");	sb.append(vcount);
		
		for(int i=0;i<vcount;i++){ sb.append("\n");	sb.append(vdef[i].toString());}
		
		return sb.toString();
	}
	
	/**
     * convert to a 2D region
     */
	public Region2D toRegion2D(){
		float[] lons=ArrayUtil.getExtrema(xdef.getSamples());
		float[] lats=ArrayUtil.getExtrema(ydef.getSamples());
		
		return new Region2D(lons[0],lats[0],lons[1],lats[1]);
	}
	
	/**
     * convert to a 3D region
     */
	public Region3D toRegion3D(){
		float[] lons=ArrayUtil.getExtrema(xdef.getSamples());
		float[] lats=ArrayUtil.getExtrema(ydef.getSamples());
		float[] levs=ArrayUtil.getExtrema(zdef.getSamples());
		
		return new Region3D(lons[0],lats[0],levs[0],lons[1],lats[1],levs[1]);
	}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(DiagnosisFactory.DF1.getDataDescriptor().getXNumPeriodicX(-360f));System.exit(0);
		
		System.out.println(DiagnosisFactory.DF10.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DF5.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DF2P5.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DF2.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DF1P5.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DF1.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DFHalf.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DFQuater.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DFThird.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DFT42.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DFT106.getDataDescriptor().isPeriodicX());
		System.out.println(DiagnosisFactory.DFT170.getDataDescriptor().isPeriodicX());
	}*/
}
