/**
 * @(#)Range.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import java.util.regex.Pattern;
import java.util.regex.Matcher;

import miniufo.descriptor.DataDescriptor;


/**
 * Describe the range in the spatial model, only recording the indices.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Range implements Cloneable{
	//
	private int[] t_range=null;	// start from 1, [0] is start index
	private int[] z_range=null;	//               [1] is end   index
	private int[] y_range=null;	//               [2] is length
	private int[] x_range=null;	//
	
	// like   x(1,144);y(1,73);z(1,20);t(1,10)
	// or     lon(0,360);lat(-90,90);lev(1000,100);time(1998.1.1,1998.12.31);
	// RegExp (t|z|y|x|time|lev|lat|lon)\((.+?),(.+?)\)
	private static Pattern pr=Pattern.compile(
		"(t|z|y|x|time|lev|lat|lon)\\((.+?),(.+?)\\)",
	Pattern.CASE_INSENSITIVE);
	
	// 1998.1.1.12.00  year.month.date.hour.minute
	// RegExp (\d{1,4})\.(\d{1,2})(\.(\d{1,2}))?(\.(\d{1,2}))?(\.(\d{1,2}))?
	private static Pattern pt=Pattern.compile(
			"(\\d{1,4})\\.(\\d{1,2})(\\.(\\d{1,2}))?(\\.(\\d{1,2}))?(\\.(\\d{1,2}))?",
	Pattern.CASE_INSENSITIVE);
	
	
	/**
     * constructor
     *
     * @param	sformations		formation in string form
     * @param	dd				DataDescriptor
     */
	public Range(String sformations,DataDescriptor dd){
		boolean t_isdef=false,z_isdef=false,y_isdef=false,x_isdef=false;
		
		Matcher m=pr.matcher(sformations.replace(" ",""));
		
		String[][] mtchs=new String[8][4];
		
		int mc=0;
		while(m.find()){ for(int i=0;i<4;i++) mtchs[mc][i]=m.group(i);	mc++;}
		
		for(int i=0;i<mc;i++){
			if(mtchs[i][1].equals("x")){
				if(x_isdef) throw new IllegalArgumentException("x has already been defined");
				
				x_range=new int[3];
				x_range[0]=Integer.parseInt(mtchs[i][2]);
				x_range[1]=Integer.parseInt(mtchs[i][3]);
				x_isdef=true;
				
			}else if(mtchs[i][1].equals("y")){
				if(y_isdef) throw new IllegalArgumentException("y has already been defined");
				
				y_range=new int[3];
				y_range[0]=Integer.parseInt(mtchs[i][2]);
				y_range[1]=Integer.parseInt(mtchs[i][3]);
				y_isdef=true;
				
			}else if(mtchs[i][1].equals("z")){
				if(z_isdef) throw new IllegalArgumentException("z has already been defined");
				
				z_range=new int[3];
				z_range[0]=Integer.parseInt(mtchs[i][2]);
				z_range[1]=Integer.parseInt(mtchs[i][3]);
				z_isdef=true;
				
			}else if(mtchs[i][1].equals("t")){
				if(t_isdef) throw new IllegalArgumentException("t has already been defined");
				
				t_range=new int[3];
				t_range[0]=Integer.parseInt(mtchs[i][2]);
				t_range[1]=Integer.parseInt(mtchs[i][3]);
				t_isdef=true;
				
			}else if(mtchs[i][1].equals("lon")){
				if(x_isdef) throw new IllegalArgumentException("x has already been defined");
				
				x_range=new int[3];
				x_range[0]=dd.getXNum(Float.parseFloat(mtchs[i][2]))+1;
				x_range[1]=dd.getXNum(Float.parseFloat(mtchs[i][3]))+1;
				x_isdef=true;
				
			}else if(mtchs[i][1].equals("lat")){
				if(y_isdef) throw new IllegalArgumentException("y has already been defined");
				
				y_range=new int[3];
				y_range[0]=dd.getYNum(Float.parseFloat(mtchs[i][2]))+1;
				y_range[1]=dd.getYNum(Float.parseFloat(mtchs[i][3]))+1;
				y_isdef=true;
				
			}else if(mtchs[i][1].equals("lev")){
				if(z_isdef) throw new IllegalArgumentException("z has already been defined");
				
				z_range=new int[3];
				z_range[0]=dd.getZNum(Float.parseFloat(mtchs[i][2]))+1;
				z_range[1]=dd.getZNum(Float.parseFloat(mtchs[i][3]))+1;
				z_isdef=true;
				
			}else{
				if(t_isdef) throw new IllegalArgumentException("t has already been defined");
				
				t_range=new int[3];
				
				Matcher mst=pt.matcher(mtchs[i][2]);	String[] stime=null;
				Matcher met=pt.matcher(mtchs[i][3]);	String[] etime=null;
				
				if(mst.find()&&met.find()){
					stime=new String[mst.groupCount()+1];
					etime=new String[mst.groupCount()+1];
					
					if(stime.length!=etime.length)
					throw new IllegalArgumentException("time format are not same");
					
					for(int ii=0;ii<=mst.groupCount();ii++){
						stime[ii]=mst.group(ii);
						etime[ii]=met.group(ii);
					}
				}
				
				t_range[0]=dd.getTNum(new MDate(
					Integer.parseInt(stime[1]),
					Integer.parseInt(stime[2]),
					stime[4]==null?1:Integer.parseInt(stime[4]),
					stime[6]==null?0:Integer.parseInt(stime[6]),
					stime[8]==null?0:Integer.parseInt(stime[8]),
				0))+1;
				
				t_range[1]=dd.getTNum(new MDate(
					Integer.parseInt(etime[1]),
					Integer.parseInt(etime[2]),
					etime[4]==null?1:Integer.parseInt(etime[4]),
					etime[6]==null?0:Integer.parseInt(etime[6]),
					etime[8]==null?0:Integer.parseInt(etime[8]),
				0))+1;
				
				t_isdef=true;
			}
		}
		
		if(!x_isdef){
			x_range=new int[3];
			x_range[0]=1;	x_range[1]=dd.getXDef().length();		x_isdef=true;
		}
		
		if(!y_isdef){
			y_range=new int[3];
			y_range[0]=1;	y_range[1]=dd.getYDef().length();		y_isdef=true;
		}
		
		if(!z_isdef){
			z_range=new int[3];
			z_range[0]=1;	z_range[1]=dd.getZDef().length();		z_isdef=true;
		}
		
		if(!t_isdef){
			t_range=new int[3];
			t_range[0]=1;	t_range[1]=dd.getTDef().length();		t_isdef=true;
		}
		
		t_range[2]=t_range[1]-t_range[0]+1;
		z_range[2]=z_range[1]-z_range[0]+1;
		y_range[2]=y_range[1]-y_range[0]+1;
		x_range[2]=x_range[1]-x_range[0]+1;
		
		if(x_range[2]<=0||y_range[2]<=0||z_range[2]<=0||t_range[2]<=0)
			throw new IllegalArgumentException("invalid range");
		
		if(x_range[0]<=0||y_range[0]<=0||z_range[0]<=0||t_range[0]<=0)
			throw new IllegalArgumentException("invalid range");
		
		if(t_range[1]>dd.getTDef().length()||z_range[1]>dd.getZDef().length())
			throw new IllegalArgumentException("range beyond limit");
		
		if(y_range[1]>dd.getYDef().length()||x_range[1]>dd.getXDef().length())
			throw new IllegalArgumentException("range beyond limit");
	}
	
	/**
     * constructor
     *
     * @param	t	t-count
     * @param	z	z-count
     * @param	y	y-count
     * @param	x	x-count
     */
	public Range(int t,int z,int y,int x){
		if(t<1) throw new IllegalArgumentException("t count should larger than 1");
		if(z<1) throw new IllegalArgumentException("z count should larger than 1");
		if(y<1) throw new IllegalArgumentException("y count should larger than 1");
		if(x<1) throw new IllegalArgumentException("x count should larger than 1");
		
		t_range=new int[3];	z_range=new int[3];
		y_range=new int[3];	x_range=new int[3];
		
		t_range[0]=z_range[0]=y_range[0]=x_range[0]=1;
		
		t_range[2]=t;	z_range[2]=z;
		y_range[2]=y;	x_range[2]=x;
		
		t_range[1]=t;	z_range[1]=z;
		y_range[1]=y;	x_range[1]=x;
	}
	
	
	/*** getor and setor ***/
	public int[] getTRange(){ return t_range;}
	
	public int[] getZRange(){ return z_range;}
	
	public int[] getYRange(){ return y_range;}
	
	public int[] getXRange(){ return x_range;}
	
	public void setTRange(int r){ t_range[0]=t_range[1]=r; t_range[2]=1;}
	
	public void setZRange(int r){ z_range[0]=z_range[1]=r; z_range[2]=1;}
	
	public void setYRange(int r){ y_range[0]=y_range[1]=r; y_range[2]=1;}
	
	public void setXRange(int r){ x_range[0]=x_range[1]=r; x_range[2]=1;}
	
	public void setTRange(Range r){ System.arraycopy(r.t_range,0,t_range,0,3);}
	
	public void setZRange(Range r){ System.arraycopy(r.z_range,0,z_range,0,3);}
	
	public void setYRange(Range r){ System.arraycopy(r.y_range,0,y_range,0,3);}
	
	public void setXRange(Range r){ System.arraycopy(r.x_range,0,x_range,0,3);}
	
	
	/**
     * used to print out
     */
	public String toString(){
		return  "Range:\n"+
				" T Range:\t"+t_range[0]+"  --  "+t_range[1]+"\t"+t_range[2]+"\n"+
				" Z Range:\t"+z_range[0]+"  --  "+z_range[1]+"\t"+z_range[2]+"\n"+
				" Y Range:\t"+y_range[0]+"  --  "+y_range[1]+"\t"+y_range[2]+"\n"+
				" X Range:\t"+x_range[0]+"  --  "+x_range[1]+"\t"+x_range[2];
	}
	
	
	/**
     * clone method
     */
    public Object clone(){
		try{
		    Range r=null;	r=(Range)super.clone();
			
			r.t_range=new int[3];	r.z_range=new int[3];
			r.y_range=new int[3];	r.x_range=new int[3];
			
			System.arraycopy(t_range,0,r.t_range,0,3);
			System.arraycopy(z_range,0,r.z_range,0,3);
			System.arraycopy(y_range,0,r.y_range,0,3);
			System.arraycopy(x_range,0,r.x_range,0,3);
			
			return r;
			
	    }catch(CloneNotSupportedException ex){
		    // this shouldn't happen, since we are Cloneable
		    throw new InternalError();
	    }
    }
	
	
	/** test
	public static void main(String[] args){
		try{
			String r="time(2006.5.14,2006.5.16);lon(120,250);lat(20,50);lev(900,500)";
			
			CtlDescriptor ctl=new CtlDescriptor("E:/typhoon/chanchu/chanchu2.ctl");
			
			long start=0;
			start=System.currentTimeMillis();
			CopyOfRange r1=new CopyOfRange(r,ctl);
			System.out.println("new: "+(System.currentTimeMillis()-start));
			
			start=System.currentTimeMillis();
			Range r2=new Range(r,ctl);
			System.out.println("old: "+(System.currentTimeMillis()-start));
			
			System.out.println(r1+"\n"+r2);
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
