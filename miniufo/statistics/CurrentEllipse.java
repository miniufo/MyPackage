/**
 * @(#)CurrentEllipse.java	1.0 2014.06.12
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.statistics;

import static java.lang.Math.abs;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import miniufo.mathsphysics.Complex;


/**
 * used for basic statistic method
 *
 * @version 1.0, 2014.06.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CurrentEllipse{
	//
	
	/**
	 * constructor
	 */
	private CurrentEllipse(){}
	
	
	/**
     * calculate parameters of variance ellipse
     * 
     * @param	ueddy	zonal component of eddy velocity time series
     * @param	veddy	meridional component of eddy velocity time series
     *
     * @return	re		major, minor axes and orientation (theta)
     */
	public static float[] cVarianceEllipse(float[] ueddy,float[] veddy){
		float uvar=squaredMean(ueddy);
		float vvar=squaredMean(veddy);
		float cvar=productMean(ueddy,veddy);
		
		float major=(uvar+vvar+(float)sqrt((uvar-vvar)*(uvar-vvar)+4*cvar*cvar))/2;
		float minor=(uvar+vvar)-major;
		float theta=(float)atan2(major-uvar,cvar);
		
		float[] re=new float[]{major,minor,theta};
		
		return re;
	}
	
	
	/**
     * calculate parameters of current ellipse given the form:
     * u = U * cos(t - fai)
     * v = V * cos(t - sai)
     * 
     * @param	uamp	zonal amplitude U
     * @param	upha	zonal phase fai
     * @param	vamp	meridional amplitude V
     * @param	vpha	meridional phase sai
     *
     * @return	re		major and minor axes, theta (inclination), and phase
     */
	public static float[] cCurrentEllipseByCoeffs(float usin,float ucos,float vsin,float vcos){
		float uamp=(float)Math.hypot(ucos,usin);
		float vamp=(float)Math.hypot(vcos,vsin);
		float upha=(float)(Math.PI/2.0-Math.atan2(ucos,usin));
		float vpha=(float)(Math.PI/2.0-Math.atan2(vcos,vsin));
		
		return cCurrentEllipseByAmps(uamp,upha,vamp,vpha);
	}
	
	
	/**
     * calculate parameters of current ellipse given the form:
     * u = U * cos(t - fai)
     * v = V * cos(t - sai)
     * 
     * @param	uamp	zonal amplitude U
     * @param	upha	zonal phase fai (radian)
     * @param	vamp	meridional amplitude V
     * @param	vpha	meridional phase sai (radian)
     *
     * @return	re		major and minor axes, theta (inclination), and phase
     */
	public static float[] cCurrentEllipseByAmps(float uamp,float upha,float vamp,float vpha){
		Complex u=Complex.polar(uamp,-upha);
		Complex v=Complex.polar(vamp,-vpha);
		
		Complex wp=new Complex(0f, 1f).multiplyEq(v).plusEq(u).divideEq(2f);
		Complex wm=new Complex(0f,-1f).multiplyEq(v).plusEq(u).conjugateEq().divideEq(2f);
		
		float Wp=wp.getMod();
		float Wm=wm.getMod();
		float THETAp=wp.getArg();
		float THETAm=wm.getArg();
		
		float major=Wp+Wm;
		float minor=abs(Wp-Wm);
		float theta=(THETAm+THETAp)/2f;
		float phase=(THETAm-THETAp)/2f;
		
		float[] re=new float[]{major,minor,theta,phase};
		
		return re;
	}
	
	public static float[] cCurrentEllipseByAmps2(float uamp,float upha,float vamp,float vpha){
		float uv=uamp*vamp,uu=uamp*uamp,vv=vamp*vamp;
		
		float cosdpha=(float)cos(upha-vpha);
		
		float theta=(float)atan2(2.0*uv*cosdpha,uu-vv)/2f;
		
		float sintheta=(float)sin(theta);
		float costheta=(float)cos(theta);
		
		float nu=uv*(float)sin(vpha-upha); nu*=nu;
		float co=uv*(float)(sin(2.0*theta)*cosdpha);
		
		float major=(float)sqrt(nu/(vv*costheta*costheta+uu*sintheta*sintheta-co));
		float minor=(float)sqrt(nu/(vv*sintheta*sintheta+uu*costheta*costheta+co));
		
		float[] re=new float[]{major,minor,theta};
		
		return re;
	}
	
	
	/*** helper methods ***/
	private static float squaredMean(float[] data){
		if(data.length==1||data.length==0) return 0;
		
		float average=StatisticsUtil.cArithmeticMean(data),variance=0;
		
		for(float ei:data){
			float tmp=ei-average;
			variance+=tmp*tmp;
		}
		
		return variance/data.length;
	}
	
	private static float productMean(float[] data1,float[] data2){
		if(data1.length==1||data1.length==0) return 0;
		
		float mean_a=StatisticsUtil.cArithmeticMean(data1);
		float mean_b=StatisticsUtil.cArithmeticMean(data2);
		float covariance=0;
		
		for(int i=0,I=data1.length;i<I;i++) covariance+=(data1[i]-mean_a)*(data2[i]-mean_b);
		
		return covariance/data1.length;
	}
	
	
	/** test
	public static void main(String[] args){
		float[] re2=null,re3=null;
		
		TicToc.tic("testing cCurrentEllipse2");
		for(int i=0;i<100000;i++) re2=cCurrentEllipse2(2,i,1,1);
		TicToc.toc();
		
		TicToc.tic("testing cCurrentEllipse3");
		for(int i=0;i<100000;i++) re3=cCurrentEllipse3(2,i,1,1);
		TicToc.toc();
	}*/
}
