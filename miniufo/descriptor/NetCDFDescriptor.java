/**
 * @(#)NetCDFDescriptor.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import ucar.nc2.Variable;
import ucar.nc2.NetcdfFile;
import miniufo.diagnosis.MDate;
import static miniufo.basic.ArrayUtil.reverse;


/**
 * used to describe the data of NetCDF format
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class NetCDFDescriptor extends DataDescriptor{
	//
	private int dim_count    =0;
	
	private boolean has_tdef=false;
	private boolean has_zdef=false;
	private boolean has_ydef=false;
	private boolean has_xdef=false;
	
	private String dumpinfo =null;
	
	private NetcdfFile nc   =null;
	
	
	/**
     * constructor
     *
     * @param	cdfPath
     *
     * @exception	if I/O error or invalid NetCDF content encountered
     */
	public NetCDFDescriptor(String cdfPath) throws IOException{
		nc=NetcdfFile.open(cdfPath);	// open it readonly
		
		hasData=true;
		
		dumpinfo=nc.toString();
		
		dsetPath=nc.getLocation();	descPath=dsetPath;
		
		Variable vrbl=nc.findVariable("lon");	if(vrbl==null) vrbl=nc.findVariable("longitude");
		if(vrbl!=null){
			xdef=new SpatialCoordinate("longitude",(float[])vrbl.read().get1DJavaArray(float.class));
			
			dim_count++;	has_xdef=true;
			
			if(xdef.length()==1){ dxdef=new float[1];	dxdef[0]=1;}
			
		}else{ dxdef=new float[1];	dxdef[0]=1;}
		
		vrbl=nc.findVariable("lat");	if(vrbl==null) vrbl=nc.findVariable("latitude");
		if(vrbl!=null){
			ydef=new SpatialCoordinate("latitude",(float[])vrbl.read().get1DJavaArray(float.class));
			
			dim_count++;	has_ydef=true;	float[] df=ydef.getSamples();
			
			if(ydef.length()==1){ dydef=new float[1];	dydef[0]=1;}
			
			if(df[0]>df[ydef.length()-1]){ reverse(df);	yrev=true;}
			
		}else{ dydef=new float[1];	dydef[0]=1;}
		
		vrbl=nc.findVariable("level");	if(vrbl==null) vrbl=nc.findVariable("levelist");
		if(vrbl!=null){
			switch(vrbl.getDataType()){
			case FLOAT:
				zdef=new SpatialCoordinate("level",(float[])vrbl.read().get1DJavaArray(float.class));
				break;
			case INT:
				int[] zdata=(int[])vrbl.read().get1DJavaArray(int.class);
				float[] fdata=new float[zdata.length];
				
				for(int i=0,I=zdata.length;i<I;i++) fdata[i]=zdata[i];
				
				zdef=new SpatialCoordinate("level",fdata);
				break;

			default:
				throw new IllegalArgumentException("unsupported level type");
			}
			
			dim_count++;	has_zdef=true;
			
			if(zdef.length()==1){ dzdef=new float[1];	dzdef[0]=50;}
			
		}else{ dzdef=new float[1];	dzdef[0]=50;}
		
		vrbl=nc.findVariable("time");	if(vrbl==null) vrbl=nc.findVariable("date");
		if(vrbl!=null){
			int tcount=vrbl.getShape(0);	MDate[] df=new MDate[tcount];
			dim_count++;	has_tdef=true;
			
			String[] arrayT=vrbl.findAttribute("units").getStringValue().split(" ");
			String   units =arrayT[0].toLowerCase();
			String[] sdate =arrayT[2].split("-");
			String[] stime =arrayT.length==4?arrayT[3].split(":"):null;
			
			MDate startT=null;
			if(stime!=null)
				startT=new MDate(
					Integer.parseInt(sdate[0]),Integer.parseInt(sdate[1]),Integer.parseInt(sdate[2]),
					Integer.parseInt(stime[0]),Integer.parseInt(stime[1])
				);
			else
				startT=new MDate(
					Integer.parseInt(sdate[0]),Integer.parseInt(sdate[1]),Integer.parseInt(sdate[2])
				);
			
			int offset=0;
			if(vrbl.findAttribute("offset")!=null)
				offset=(int)(vrbl.findAttribute("offset").getNumericValue().floatValue());

			double[] dbl_arr=null;
			switch(vrbl.getDataType()){
			case DOUBLE:
				dbl_arr=(double[])vrbl.read().get1DJavaArray(double.class);
				break;
			case FLOAT:
				dbl_arr=new double[tcount];
				float[] arrf=(float[])vrbl.read().get1DJavaArray(float.class);
				for(int i=0;i<tcount;i++) dbl_arr[i]=arrf[i];
				break;
			case INT:
				dbl_arr=new double[tcount];
				int[] arri=(int[])vrbl.read().get1DJavaArray(int.class);
				for(int i=0;i<tcount;i++) dbl_arr[i]=arri[i];
				break;
			case LONG:
				dbl_arr=new double[tcount];
				long[] arrl=(long[])vrbl.read().get1DJavaArray(long.class);
				for(int i=0;i<tcount;i++) dbl_arr[i]=arrl[i];
				break;

			default:
				throw new IllegalArgumentException("the type of time is "+vrbl.getDataType());
			}
			
			dtdef=new float[1];
			
			if(vrbl.findAttribute("delta_t")!=null) dtdef[0]=getDeltaT(vrbl.findAttribute("delta_t").getStringValue());
			else dtdef[0]=getDeltaT(units);
			
			if(vrbl.findAttribute("delta_t")==null&&tcount!=0&&tcount!=1&&dtdef[0]!=dbl_arr[1]-dbl_arr[0]){
				dtdef[0]=(float)((dbl_arr[1]-dbl_arr[0])*dtdef[0]);
				System.out.println("warning for debug in nc");
				tincrement=(int)(dbl_arr[1]-dbl_arr[0])+tincrement.substring(tincrement.length()-2);
			}
			
			if(units.equals("hours")||units.equals("hour")){
				df[0]=startT.addHours((int)dbl_arr[0]+offset);
				for(int i=1;i<tcount;i++) df[i]=startT.addHours((int)dbl_arr[i]+offset);
				
			}else if(units.equals("days")||units.equals("day")){
				df[0]=startT.addDays((int)dbl_arr[0]+offset);
				for(int i=1;i<tcount;i++) df[i]=startT.addDays((int)dbl_arr[i]+offset);
				
			}else if(units.equals("months")||units.equals("month")){
				df[0]=startT.addMonths((int)dbl_arr[0]+offset);
				for(int i=1;i<tcount;i++) df[i]=startT.addMonths((int)dbl_arr[i]+offset);
				
			}else if(units.equals("minutes")||units.equals("minute")){
				df[0]=startT.addMinutes((int)dbl_arr[0]+offset);
				for(int i=1;i<tcount;i++) df[i]=startT.addMinutes((int)dbl_arr[i]+offset);
				
			}else if(units.equals("years")||units.equals("year")){
				df[0]=startT.addYears((int)dbl_arr[0]+offset);
				for(int i=1;i<tcount;i++) df[i]=startT.addYears((int)dbl_arr[i]+offset);
				
			}else throw new IllegalArgumentException("units is "+units);
			
			tdef=new TemporalCoordinate("time",df,true);
		}
		
		if(nc.findGlobalAttributeIgnoreCase("title")!=null)
		title=nc.findGlobalAttributeIgnoreCase("title").getStringValue();
		
		if(!has_tdef){
			System.out.println("\n  Warning: no tdef in "+descPath);
			dtdef=new float[1]; dtdef[0]=3600;
			tdef=new TemporalCoordinate("time",new MDate[]{new MDate(1970,1,1)},true);	has_tdef=true;
		}
		
		if(!has_zdef){
			System.out.println("\n  Warning: no zdef in "+descPath);
			zdef=new SpatialCoordinate("level",new float[]{1000});
		}
		
		if(!has_ydef){
			System.out.println("\n  Warning: no ydef in "+descPath);
			ydef=new SpatialCoordinate("latitude",new float[]{0});
		}
		
		if(!has_xdef){
			System.out.println("\n  Warning: no xdef in "+descPath);
			xdef=new SpatialCoordinate("longitude",new float[]{0});
		}
		
		vcount=nc.getVariables().size()-dim_count;	vdef=new NetCDFVar[vcount];
		
		int i=0;
		for(Variable vv:nc.getVariables())
		if(!vv.isCoordinateVariable()) vdef[i++]=new NetCDFVar(vv);
		
		postProcess();	nc.close();
	}
	
	
	/*** getor and setor ***/
	public int getDimensionCount(){ return dim_count;}
	
	public float getUndef(String vname){
		for(int m=0;m<vcount;m++) if(vdef[m].getName().equals(vname)) return vdef[m].getUndef();
		
		if(vname==null) return vdef[0].getUndef();
		else throw new IllegalArgumentException("Cannot find "+vname+" in "+descPath);
	}
	
	public boolean hasTDef(){ return has_tdef;}
	
	public boolean hasZDef(){ return has_zdef;}
	
	public boolean hasYDef(){ return has_ydef;}
	
	public boolean hasXDef(){ return has_xdef;}
	
	public NetcdfFile getNetCDFFile(){ return nc;}
	
	
	/**
     * dump all file information
     */
	public void dump(){ System.out.println(dumpinfo);}
	
	public static void dump(String path){
		try{
			NetcdfFile nf=NetcdfFile.open(path);
			System.out.println(nf);
			nf.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/**
     * get delta t from String to float
     *
     * @param	deltaT	delta time in String form
     *
     * @return	dt		delta time in float form
     */
	private float getDeltaT(String deltaT){
		if(deltaT.equals("hours"       )){ tincrement="1hr"; return 60*60;       }
		else if(deltaT.equals("days"   )){ tincrement="1dy"; return 60*60*24;    }
		else if(deltaT.equals("minutes")){ tincrement="1mn"; return 60;          }
		else if(deltaT.equals("years"  )){ tincrement="1yr"; return 60*60*24*365;}
		
		final Matcher m=
		Pattern.compile("(\\d{4})-(\\d{2})-(\\d{2}) (\\d{2}):(\\d{2}):(\\d{2})").matcher(deltaT);
		
		if(!m.find()) throw new IllegalArgumentException("no matched delta T in "+descPath);
		
		int yr,mo,dy,hr,mn,se;
		
		yr=Integer.parseInt(m.group(1)); if(yr!=0){ tincrement=yr+"yr"; return yr*60*60*24*365;}
		mo=Integer.parseInt(m.group(2)); if(mo!=0){ tincrement=mo+"mo"; return mo*60*60*24*31 ;} // probably no use
		dy=Integer.parseInt(m.group(3)); if(dy!=0){ tincrement=dy+"dy"; return dy*60*60*24    ;}
		hr=Integer.parseInt(m.group(4)); if(hr!=0){ tincrement=hr+"hr"; return hr*60*60       ;}
		mn=Integer.parseInt(m.group(5)); if(mn!=0){ tincrement=mn+"mn"; return mn*60          ;}
		se=Integer.parseInt(m.group(6)); if(se!=0){ tincrement=se+"se"; return se             ;}
		
		System.out.println("  Warning: Delta T is not defined, 1 hour is selected");
		
		tincrement="1hr";
		
		return 60*60;
	}
	
	
	/** test
	public static void main(String arg[]){
		try{
			NetCDFDescriptor nc=new NetCDFDescriptor("//lynn/DataI/NCEP1/Monthly/Surface/land.nc");
			nc.dump();
			System.out.println(nc);
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
