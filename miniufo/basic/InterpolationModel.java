/**
 * @(#)InterpolationModel.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.basic;

import java.util.ArrayList;
import java.util.List;
import miniufo.mathsphysics.Spline;
import miniufo.statistics.FilterModel;


/**
 * interpolation methods
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class InterpolationModel{
	//
	public enum Type{LINEAR,CUBIC_L,CUBIC_P,SPLINE,PERIODIC_LINEAR,PERIODIC_CUBIC_L,PERIODIC_CUBIC_P,PERIODIC_SPLINE};
	
	
	/**
	 * prevent from instantiate
	 */
	private InterpolationModel(){}
	
	
	/**
     * Linear interpolation y=f(x) using two point values.
     * 
     * dis within [0,1] for interpolation and else for extrapolation
     *
     * @param	y0		value of the first point
     * @param	y1		value of the second point
     * @param	dis		ratio of the distance defined as (x-x0)/(x1-x0)
     */
	public static float linearInterpolation(float y0,float y1,float dis){ return y0+(y1-y0)*dis;}
	
	public static double linearInterpolation(double y0,double y1,double dis){ return y0+(y1-y0)*dis;}
	
	
	/**
     * Linear interpolation y=f(x) using two point values.
     * 
     * dis within [0,1] for interpolation and else for extrapolation
     *
     * @param	y0		value of the first point
     * @param	y1		value of the second point
     * @param	dis		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	undef	undefined value
     */
	public static float linearInterpolation(float y0,float y1,float dis,float undef){
		if(y0!=undef){
			if(y1!=undef) return y0+(y1-y0)*dis;
			else return y0;
		
		}else{
			if(y1!=undef) return y1;
			else return undef;
		}
	}
	
	public static double linearInterpolation(double y0,double y1,double dis,double undef){
		if(y0!=undef){
			if(y1!=undef) return y0+(y1-y0)*dis;
			else return y0;
		
		}else{
			if(y1!=undef) return y1;
			else return undef;
		}
	}
	
	
	/**
     * Quadratic interpolation y=f(x) using three equal-space values
     * 
     * No x-coordinates as it is assumed that x are evenly distributed
     * 
     * dis within [0,2] for interpolation and else for extrapolation
     * 
     * Equivalent to nLagrangeInterpolation(new float[]{0,1,2},new float[]{x0,x1,x2},dis)
     *
     * @param	y0		value of the first point
     * @param	y1		value of the second point
     * @param	y2		value of the third point
     * @param	dis		ratio of the distance defined as (x-x0)/(x1-x0)
     */
	public static float quadraticLagrangeInterpolation(float y0,float y1,float y2,float dis){
		float re=y0;
		re+=dis*(-1.5*y0+2.0*y1-0.5*y2);
		re+=dis*dis*(0.5*y0-y1+0.5*y2);
		return re;
	}
	
	public static double quadraticLagrangeInterpolation(double y0,double y1,double y2,double dis){
		double re=y0;
		re+=dis*(-1.5*y0+2.0*y1-0.5*y2);
		re+=dis*dis*(0.5*y0-y1+0.5*y2);
		return re;
	}
	
	/**
     * Quadratic interpolation y=f(x) using three equal-space values
     * 
     * No x-coordinates as it is assumed that x are evenly distributed
     * 
     * dis within [0,2] for interpolation and else for extrapolation
     * 
     * Equivalent to nLagrangeInterpolation(new float[]{0,1,2},new float[]{x0,x1,x2},dis)
     *
     * @param	y0		value of the first point
     * @param	y1		value of the second point
     * @param	y2		value of the third point
     * @param	dis		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	undef	undefined value
     */
	public static float quadraticLagrangeInterpolation(float y0,float y1,float y2,float dis,float undef){
		if(y0!=undef){
			if(y1!=undef){
				if(y2!=undef) return quadraticLagrangeInterpolation(y0,y1,y2,dis);
				else{
					if(dis<=1) return linearInterpolation(y0,y1,dis);
					else return undef;
				}
				
			}else{
				if(y2!=undef) return linearInterpolation(y0,y2,dis/2f);
				else return undef;
			}
			
		}else{
			if(y1!=undef){
				if(y2!=undef){
					if(dis>1) return linearInterpolation(y1,y2,dis-1f);
					else return undef;
					
				}else return undef;
				
			}else return undef;
		}
	}
	
	public static double quadraticLagrangeInterpolation(double y0,double y1,double y2,double dis,double undef){
		if(y0!=undef){
			if(y1!=undef){
				if(y2!=undef) return quadraticLagrangeInterpolation(y0,y1,y2,dis);
				else{
					if(dis<=1) return linearInterpolation(y0,y1,dis);
					else return undef;
				}
				
			}else{
				if(y2!=undef) return linearInterpolation(y0,y2,dis/2.0);
				else return undef;
			}
			
		}else{
			if(y1!=undef){
				if(y2!=undef){
					if(dis>1) return linearInterpolation(y1,y2,dis-1.0);
					else return undef;
					
				}else return undef;
				
			}else return undef;
		}
	}
	
	
	/**
     * Cubic Lagrange interpolation y=f(x) using four equal-space values.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
     * 
     * dis within [0,3] for interpolation and else for extrapolation
     * 
     * equivalent to nLagrangeInterpolation(new float[]{0,1,2,3},new float[]{y0,y1,y2,y3},dis)
     *
     * @param	y0		value of the first point
     * @param	y1		value of the second point
     * @param	y2		value of the third point
     * @param	y3		value of the forth point
     * @param	dis		ratio of the distance defined as (x-x0)/(x1-x0)
     *
     * @return	re		result of interpolation
     */
	public static float cubicLagrangeInterpolation(float y0,float y1,float y2,float y3,float dis){
		// float re=y0;
		// re+=dis*(11.0f*y0/(-6.0f)+y1*3.0f-3.0f*y2/2.0f+y3/3.0f);
		// re+=dis*dis*(y0-5.0f*y1/2.0f+2.0f*y2-y3/2.0f);
		// re+=dis*dis*dis*((y3-y0)/6.0f+(y1-y2)/2.0f);
		// return re;
		float re=dis*(-11f*y0+18f*y1-9f*y2+2f*y3);
		re+=dis*dis*(6f*y0-15f*y1+12f*y2-3f*y3);
		re+=dis*dis*dis*((y3-y0)+3f*(y1-y2));
		re/=6f;
		re+=y0;
		return re;
	}
	
	public static double cubicLagrangeInterpolation(double y0,double y1,double y2,double y3,double dis){
		double re=dis*(-11.0*y0+18.0*y1-9.0*y2+2.0*y3);
		re+=dis*dis*(6.0*y0-15.0*y1+12.0*y2-3.0*y3);
		re+=dis*dis*dis*((y3-y0)+3.0*(y1-y2));
		re/=6.0;
		re+=y0;
		return re;
	}
	
	/**
     * Cubic Lagrange interpolation y=f(x) using four equal-space values.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
     * 
     * dis within [0,3] for interpolation and else for extrapolation
     * 
     * equivalent to nLagrangeInterpolation(new float[]{0,1,2,3},new float[]{y0,y1,y2,y3},dis)
     * 
     * @param	y0		value of the first point
     * @param	y1		value of the second point
     * @param	y2		value of the third point
     * @param	y3		value of the forth point
     * @param	dis		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	undef	undefined value
     */
	public static float cubicLagrangeInterpolation(float y0,float y1,float y2,float y3,float dis,float undef){
		if(y0!=undef){
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef) return cubicLagrangeInterpolation(y0,y1,y2,y3,dis); // 0 1 2 3
					else return quadraticLagrangeInterpolation(y0,y1,y2,dis); // 0 1 2 #
					
				}else{
					if(y3!=undef){
						float[] x={0,1,3};	float[] y={y0,y1,y3};
						
						return nLagrangeInterpolation(x,y,dis); // 0 1 # 3
						
					}else{ // 0 1 # #
						if(dis<=1) return linearInterpolation(y0,y1,dis);
						else return undef;
					}
				}
				
			}else{
				if(y2!=undef){
					if(y3!=undef){	// 0 # 2 3
						float[] x={0,2,3};	float[] y={y0,y2,y3};
						
						return nLagrangeInterpolation(x,y,dis);
						
					}else{ // 0 # 2 #
						if(dis<=2) return linearInterpolation(y0,y2,dis/2f);
						else return undef;
					}
					
				}else{
					if(y3!=undef) return linearInterpolation(y0,y3,dis/3f); // 0 # # 3
					else return undef; // 0 # # #
				}
			}
			
		}else{
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef){	// # 1 2 3
						if(dis>=1) return quadraticLagrangeInterpolation(y1,y2,y3,dis-1);
						else return undef;
						
					}else{ // # 1 2 #
						if(dis<1) return undef;
						if(dis>2) return undef;
						return linearInterpolation(y1,y2,dis-1f);
					}
					
				}else{
					if(y3!=undef){ // # 1 # 3
						if(dis>=1) return linearInterpolation(y1,y3,(dis-1f)/2f);
						else return undef;
						
					}else return undef; // # 1 # #
				}
				
			}else return undef; // # # # #
		}
	}
	
	public static double cubicLagrangeInterpolation(double y0,double y1,double y2,double y3,double dis,double undef){
		if(y0!=undef){
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef) return cubicLagrangeInterpolation(y0,y1,y2,y3,dis); // 0 1 2 3
					else return quadraticLagrangeInterpolation(y0,y1,y2,dis); // 0 1 2 #
					
				}else{
					if(y3!=undef){
						double[] x={0,1,3};	double[] y={y0,y1,y3};
						
						return nLagrangeInterpolation(x,y,dis); // 0 1 # 3
						
					}else{ // 0 1 # #
						if(dis<=1) return linearInterpolation(y0,y1,dis);
						else return undef;
					}
				}
				
			}else{
				if(y2!=undef){
					if(y3!=undef){	// 0 # 2 3
						double[] x={0,2,3};	double[] y={y0,y2,y3};
						
						return nLagrangeInterpolation(x,y,dis);
						
					}else{ // 0 # 2 #
						if(dis<=2) return linearInterpolation(y0,y2,dis/2.0);
						else return undef;
					}
					
				}else{
					if(y3!=undef) return linearInterpolation(y0,y3,dis/3.0); // 0 # # 3
					else return undef; // 0 # # #
				}
			}
			
		}else{
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef){	// # 1 2 3
						if(dis>=1) return quadraticLagrangeInterpolation(y1,y2,y3,dis-1.0);
						else return undef;
						
					}else{ // # 1 2 #
						if(dis<1) return undef;
						if(dis>2) return undef;
						return linearInterpolation(y1,y2,dis-1.0);
					}
					
				}else{
					if(y3!=undef){ // # 1 # 3
						if(dis>=1) return linearInterpolation(y1,y3,(dis-1.0)/2.0);
						else return undef;
						
					}else return undef; // # 1 # #
				}
				
			}else return undef; // # # # #
		}
	}
	
	
	/**
     * Cubic polynomial interpolation y=f(x) using four equal-space values.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
	 * 
	 * Note that dis should be within [0,1] and f(0)=x1 and f(1)=x2, or else for extrapolation
	 * 
	 * x0 and x3 are used to approximate the first derivatives at x1 and x2 through 2nd central difference.
     *
     * @param	y0		value of the first point
     * @param	y1		value of the second point
     * @param	y2		value of the third point
     * @param	y3		value of the forth point
     * @param	dis		ratio of the distance defined as (x-x1)/(x2-x1)
     */
	public static float cubicPolynomialInterpolation(float y0,float y1,float y2,float y3,float dis){
		float re=y1;
		re+=dis*(y2-y0)*0.5f;
		re+=dis*dis*(y0-2.5f*y1+2f*y2-0.5f*y3);
		re+=dis*dis*dis*(0.5f*(y3-y0)+1.5f*(y1-y2));
		return re;
	}
	
	public static double cubicPolynomialInterpolation(double y0,double y1,double y2,double y3,double dis){
		double re=y1;
		re+=dis*(y2-y0)*0.5;
		re+=dis*dis*(y0-2.5*y1+2.0*y2-0.5*y3);
		re+=dis*dis*dis*(0.5*(y3-y0)+1.5*(y1-y2));
		return re;
	}
	
	/**
     * Cubic polynomial interpolation y=f(x) using four equal-space values.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
	 * 
	 * Note that dis should be within [0,1] and f(0)=x1 and f(1)=x2, or else for extrapolation
	 * 
	 * x0 and x3 are used to approximate the first derivatives at x1 and x2 through 2nd central difference.
     *
     * @param	x0		value of the first point
     * @param	x1		value of the second point
     * @param	x2		value of the third point
     * @param	x3		value of the forth point
     * @param	dis		ratio of the distance defined as (x-x1)/(x2-x1)
     * @param	undef	undefined value
     */
	public static float cubicPolynomialInterpolation(float y0,float y1,float y2,float y3,float dis,float undef){
		if(y0!=undef){
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef) return cubicPolynomialInterpolation(y0,y1,y2,y3,dis); // 0 1 2 3
					else return cubicPolynomialInterpolation(y0,y1,y2,y2,dis); // 0 1 2 #
					
				}else{
					if(y3!=undef){
						float[] x={0,1,3};	float[] y={y0,y1,y3};
						
						return nLagrangeInterpolation(x,y,dis+1); // 0 1 # 3
						
					}else return undef; // 0 1 # #
				}
				
			}else{
				if(y2!=undef){
					if(y3!=undef){	// 0 # 2 3
						float[] x={0,2,3};	float[] y={y0,y2,y3};
						
						return nLagrangeInterpolation(x,y,dis+1);
						
					}else{ // 0 # 2 #
						if(dis<=2) return linearInterpolation(y0,y2,(dis+1f)/2f);
						else return undef;
					}
					
				}else{
					if(y3!=undef) return linearInterpolation(y0,y3,(dis+1f)/3f); // 0 # # 3
					else return undef; // 0 # # #
				}
			}
			
		}else{
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef) return cubicPolynomialInterpolation(y1,y1,y2,y3,dis); // # 1 2 3
					else return linearInterpolation(y1,y2,dis); // # 1 2 #
					
				}else{
					if(y3!=undef) return linearInterpolation(y1,y3,dis/2f); // # 1 # 3
					else return undef; // # 1 # #
				}
				
			}else return undef; // # # # #
		}
	}
	
	public static double cubicPolynomialInterpolation(double y0,double y1,double y2,double y3,double dis,double undef){
		if(y0!=undef){
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef) return cubicPolynomialInterpolation(y0,y1,y2,y3,dis); // 0 1 2 3
					else return cubicPolynomialInterpolation(y0,y1,y2,y2,dis); // 0 1 2 #
					
				}else{
					if(y3!=undef){
						double[] x={0,1,3};	double[] y={y0,y1,y3};
						
						return nLagrangeInterpolation(x,y,dis+1.0); // 0 1 # 3
						
					}else return undef; // 0 1 # #
				}
				
			}else{
				if(y2!=undef){
					if(y3!=undef){	// 0 # 2 3
						double[] x={0,2,3};	double[] y={y0,y2,y3};
						
						return nLagrangeInterpolation(x,y,dis+1.0);
						
					}else{ // 0 # 2 #
						if(dis<=2) return linearInterpolation(y0,y2,(dis+1.0)/2.0);
						else return undef;
					}
					
				}else{
					if(y3!=undef) return linearInterpolation(y0,y3,(dis+1.0)/3.0); // 0 # # 3
					else return undef; // 0 # # #
				}
			}
			
		}else{
			if(y1!=undef){
				if(y2!=undef){
					if(y3!=undef) return cubicPolynomialInterpolation(y1,y1,y2,y3,dis); // # 1 2 3
					else return linearInterpolation(y1,y2,dis); // # 1 2 #
					
				}else{
					if(y3!=undef) return linearInterpolation(y1,y3,dis/2.0); // # 1 # 3
					else return undef; // # 1 # #
				}
				
			}else return undef; // # # # #
		}
	}
	
	
	/**
     * Bi-linear interpolation z = f(x,y) using neighboring four point.
     * Reference: Lee et al. 2005 OSJ
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     * 
     * dx1 and dy1 within [0,1] for interpolation and else for extrapolation
     *
     * @param	z11		lower-left point
     * @param	z12		lower-right point
     * @param	z21		upper-left point
     * @param	z22		upper-right point
     * @param	dx1		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	dy1		ratio of the distance defined as (y-y0)/(y1-y0)
     */
	public static float bilinearInterpolation(float z11,float z12,float z21,float z22,float dx1,float dy1){
		float dx2=1f-dx1,dy2=1f-dy1;
		
		return z11*dx2*dy2+z12*dx1*dy2+z21*dx2*dy1+z22*dx1*dy1;
	}
	
	public static double bilinearInterpolation(double z11,double z12,double z21,double z22,double dx1,double dy1){
		double dx2=1.0-dx1,dy2=1.0-dy1;
		
		return z11*dx2*dy2+z12*dx1*dy2+z21*dx2*dy1+z22*dx1*dy1;
	}
	
	/**
     * Bi-linear interpolation z = f(x,y) using neighboring four point.
     * Reference: Li et al. 2005 OSJ
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     * 
     * dx1 and dy1 within [0,1] for interpolation and else for extrapolation
     *
     * @param	z11		lower-left point
     * @param	z12		lower-right point
     * @param	z21		upper-left point
     * @param	z22		upper-right point
     * @param	dx1		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	dy1		ratio of the distance defined as (y-y0)/(y1-y0)
     * @param	undef	undefined value
     */
	public static float bilinearInterpolation(float z11,float z12,float z21,float z22,float dx1,float dy1,float undef){
		if(z11==undef) return undef;
		if(z12==undef) return undef;
		if(z21==undef) return undef;
		if(z22==undef) return undef;
		
		return bilinearInterpolation(z11,z12,z21,z22,dx1,dy1);
	}
	
	public static double bilinearInterpolation(double z11,double z12,double z21,double z22,double dx1,double dy1,double undef){
		if(z11==undef) return undef;
		if(z12==undef) return undef;
		if(z21==undef) return undef;
		if(z22==undef) return undef;
		
		return bilinearInterpolation(z11,z12,z21,z22,dx1,dy1);
	}
	
	
	/**
     * Bi-cubic Lagrange interpolation z = f(x,y) using neighboring 16 point.
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     *
     * @param	z00		lower-left point
     * @param	z01		lower-1 point
     * @param	z02		lower-2 point
     * @param	z03		lower-right point
     * @param	z10		1-left point
     * @param	z11		1-1 point
     * @param	z12		1-2 point
     * @param	z13		1-right point
     * @param	z20		2-left point
     * @param	z21		2-1 point
     * @param	z22		2-2 point
     * @param	z23		2-right point
     * @param	z30		upper-left point
     * @param	z31		upper-1 point
     * @param	z32		upper-2 point
     * @param	z33		upper-right point
     * @param	dx		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	dy		ratio of the distance defined as (y-y0)/(y1-y0)
     */
	public static float bicubicLagrangeInterpolation(
		float z00,float z01,float z02,float z03,float z10,float z11,float z12,float z13,
		float z20,float z21,float z22,float z23,float z30,float z31,float z32,float z33,
		float dx,float dy
	){
		float buf0=cubicLagrangeInterpolation(z00,z01,z02,z03,dx);
		float buf1=cubicLagrangeInterpolation(z10,z11,z12,z13,dx);
		float buf2=cubicLagrangeInterpolation(z20,z21,z22,z23,dx);
		float buf3=cubicLagrangeInterpolation(z30,z31,z32,z33,dx);
		
		return cubicLagrangeInterpolation(buf0,buf1,buf2,buf3,dy);
	}
	
	public static double bicubicLagrangeInterpolation(
		double z00,double z01,double z02,double z03,double z10,double z11,double z12,double z13,
		double z20,double z21,double z22,double z23,double z30,double z31,double z32,double z33,
		double dx,double dy
	){
		double buf0=cubicLagrangeInterpolation(z00,z01,z02,z03,dx);
		double buf1=cubicLagrangeInterpolation(z10,z11,z12,z13,dx);
		double buf2=cubicLagrangeInterpolation(z20,z21,z22,z23,dx);
		double buf3=cubicLagrangeInterpolation(z30,z31,z32,z33,dx);
		
		return cubicLagrangeInterpolation(buf0,buf1,buf2,buf3,dy);
	}
	
	/**
     * Bi-cubic Lagrange interpolation z = f(x,y) using neighboring 16 point.
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     *
     * @param	z00		lower-left point
     * @param	z01		lower-1 point
     * @param	z02		lower-2 point
     * @param	z03		lower-right point
     * @param	z10		1-left point
     * @param	z11		1-1 point
     * @param	z12		1-2 point
     * @param	z13		1-right point
     * @param	z20		2-left point
     * @param	z21		2-1 point
     * @param	z22		2-2 point
     * @param	z23		2-right point
     * @param	z30		upper-left point
     * @param	z31		upper-1 point
     * @param	z32		upper-2 point
     * @param	z33		upper-right point
     * @param	dx		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	dy		ratio of the distance defined as (y-y0)/(y1-y0)
     * @param	undef	undefined value
     */
	public static float bicubicLagrangeInterpolation(
		float z00,float z01,float z02,float z03, float z10,float z11,float z12,float z13,
		float z20,float z21,float z22,float z23, float z30,float z31,float z32,float z33,
		float dx,float dy,float undef
	){
		float buf0=cubicLagrangeInterpolation(z00,z01,z02,z03,dx,undef);
		float buf1=cubicLagrangeInterpolation(z10,z11,z12,z13,dx,undef);
		float buf2=cubicLagrangeInterpolation(z20,z21,z22,z23,dx,undef);
		float buf3=cubicLagrangeInterpolation(z30,z31,z32,z33,dx,undef);
		
		return cubicLagrangeInterpolation(buf0,buf1,buf2,buf3,dy,undef);
	}
	
	public static double bicubicLagrangeInterpolation(
		double z00,double z01,double z02,double z03,double z10,double z11,double z12,double z13,
		double z20,double z21,double z22,double z23,double z30,double z31,double z32,double z33,
		double dx,double dy,double undef
	){
		double buf0=cubicLagrangeInterpolation(z00,z01,z02,z03,dx,undef);
		double buf1=cubicLagrangeInterpolation(z10,z11,z12,z13,dx,undef);
		double buf2=cubicLagrangeInterpolation(z20,z21,z22,z23,dx,undef);
		double buf3=cubicLagrangeInterpolation(z30,z31,z32,z33,dx,undef);
		
		return cubicLagrangeInterpolation(buf0,buf1,buf2,buf3,dy,undef);
	}
	
	
	/**
     * Bi-cubic polynomial interpolation z = f(x,y) using neighboring 16 point.
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     *
     * @param	z00		lower-left point
     * @param	z01		lower-1 point
     * @param	z02		lower-2 point
     * @param	z03		lower-right point
     * @param	z10		1-left point
     * @param	z11		1-1 point
     * @param	z12		1-2 point
     * @param	z13		1-right point
     * @param	z20		2-left point
     * @param	z21		2-1 point
     * @param	z22		2-2 point
     * @param	z23		2-right point
     * @param	z30		upper-left point
     * @param	z31		upper-1 point
     * @param	z32		upper-2 point
     * @param	z33		upper-right point
     * @param	dx		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	dy		ratio of the distance defined as (y-y0)/(y1-y0)
     */
	public static float bicubicPolynomialInterpolation(
		float z00,float z01,float z02,float z03,float z10,float z11,float z12,float z13,
		float z20,float z21,float z22,float z23,float z30,float z31,float z32,float z33,
		float dx,float dy
	){
		float buf0=cubicPolynomialInterpolation(z00,z01,z02,z03,dx);
		float buf1=cubicPolynomialInterpolation(z10,z11,z12,z13,dx);
		float buf2=cubicPolynomialInterpolation(z20,z21,z22,z23,dx);
		float buf3=cubicPolynomialInterpolation(z30,z31,z32,z33,dx);
		
		return cubicPolynomialInterpolation(buf0,buf1,buf2,buf3,dy);
	}
	
	public static double bicubicPolynomialInterpolation(
		double z00,double z01,double z02,double z03,double z10,double z11,double z12,double z13,
		double z20,double z21,double z22,double z23,double z30,double z31,double z32,double z33,
		double dx,double dy
	){
		double buf0=cubicPolynomialInterpolation(z00,z01,z02,z03,dx);
		double buf1=cubicPolynomialInterpolation(z10,z11,z12,z13,dx);
		double buf2=cubicPolynomialInterpolation(z20,z21,z22,z23,dx);
		double buf3=cubicPolynomialInterpolation(z30,z31,z32,z33,dx);
		
		return cubicPolynomialInterpolation(buf0,buf1,buf2,buf3,dy);
	}
	
	/**
     * Bi-cubic polynomial interpolation z = f(x,y) using neighboring 16 point.
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     *
     * @param	z00		lower-left point
     * @param	z01		lower-1 point
     * @param	z02		lower-2 point
     * @param	z03		lower-right point
     * @param	z10		1-left point
     * @param	z11		1-1 point
     * @param	z12		1-2 point
     * @param	z13		1-right point
     * @param	z20		2-left point
     * @param	z21		2-1 point
     * @param	z22		2-2 point
     * @param	z23		2-right point
     * @param	z30		upper-left point
     * @param	z31		upper-1 point
     * @param	z32		upper-2 point
     * @param	z33		upper-right point
     * @param	dx		ratio of the distance defined as (x-x0)/(x1-x0)
     * @param	dy		ratio of the distance defined as (y-y0)/(y1-y0)
     * @param	undef	undefined value
     */
	public static float bicubicPolynomialInterpolation(
		float z00,float z01,float z02,float z03,float z10,float z11,float z12,float z13,
		float z20,float z21,float z22,float z23,float z30,float z31,float z32,float z33,
		float dx,float dy,float undef
	){
		float buf0=cubicPolynomialInterpolation(z00,z01,z02,z03,dx,undef);
		float buf1=cubicPolynomialInterpolation(z10,z11,z12,z13,dx,undef);
		float buf2=cubicPolynomialInterpolation(z20,z21,z22,z23,dx,undef);
		float buf3=cubicPolynomialInterpolation(z30,z31,z32,z33,dx,undef);
		
		return cubicPolynomialInterpolation(buf0,buf1,buf2,buf3,dy,undef);
	}
	
	public static double bicubicPolynomialInterpolation(
		double z00,double z01,double z02,double z03,double z10,double z11,double z12,double z13,
		double z20,double z21,double z22,double z23,double z30,double z31,double z32,double z33,
		double dx,double dy,double undef
	){
		double buf0=cubicPolynomialInterpolation(z00,z01,z02,z03,dx,undef);
		double buf1=cubicPolynomialInterpolation(z10,z11,z12,z13,dx,undef);
		double buf2=cubicPolynomialInterpolation(z20,z21,z22,z23,dx,undef);
		double buf3=cubicPolynomialInterpolation(z30,z31,z32,z33,dx,undef);
		
		return cubicPolynomialInterpolation(buf0,buf1,buf2,buf3,dy,undef);
	}
	
	
	/**
     * Linear interpolation y=f(x) given two data pairs y0=f(x0) and y1=f(x1) interpolate for yv=f(xv).
     *
     * @param	x0		1st x points
     * @param	x1		2nd x points
     * @param	y0		1st y points
     * @param	y1		2nd y points
     * @param	xv		value of x to be evaluated
     *
     * @return	re		result of interpolation at xv
     */
	public static float linearInterpolation(float x0,float x1,float y0,float y1,float xv){
		return linearInterpolation(y0,y1,(xv-x0)/(x1-x0));
	}
	
	public static double linearInterpolation(double x0,double x1,double y0,double y1,double xv){
		return linearInterpolation(y0,y1,(xv-x0)/(x1-x0));
	}
	
	
	/**
     * Linear interpolation y=f(x) given two data pairs y0=f(x0) and y1=f(x1) interpolate for yv=f(xv).
     *
     * @param	x0		1st x points
     * @param	x1		2nd x points
     * @param	y0		1st y points
     * @param	y1		2nd y points
     * @param	xv		value of x to be evaluated
     * @param	undef	undefined value
     *
     * @return	re		result of interpolation at xv
     */
	public static float linearInterpolation(float x0,float x1,float y0,float y1,float xv,float undef){
		if(y0==undef||y1==undef) return undef;
		return linearInterpolation(y0,y1,(xv-x0)/(x1-x0));
	}
	
	public static double linearInterpolation(double x0,double x1,double y0,double y1,double xv,double undef){
		if(y0==undef||y1==undef) return undef;
		return linearInterpolation(y0,y1,(xv-x0)/(x1-x0));
	}
	
	
	/**
     * Cubic Lagrangian interpolation y=f(x) given four data pairs
     * y0=f(x0), y1=f(x1), y2=f(x2) and y3=f(x3) interpolate for yv=f(xv).
     * 
     * This is suitable for the case that four data points are not evenly distributed.
     * 
     * xv within [x0,x3] for interpolation and else for extrapolation
     *
     * @param	x0		1st x points
     * @param	x1		2nd x points
     * @param	x2		3rd x points
     * @param	x3		4th x points
     * @param	y0		1st y points
     * @param	y1		2nd y points
     * @param	y2		3rd y points
     * @param	y3		4th y points
     * @param	xv		value of x to be evaluated
     *
     * @return	re		result of interpolation at xv
     */
	public static float cubicLagrangeInterpolation
	(float x0,float x1,float x2,float x3,float y0,float y1,float y2,float y3,float xv){
		return nLagrangeInterpolation(new float[]{x0,x1,x2,x3},new float[]{y0,y1,y2,y3},xv);
	}
	
	public static double cubicLagrangeInterpolation
	(double x0,double x1,double x2,double x3,double y0,double y1,double y2,double y3,double xv){
		return nLagrangeInterpolation(new double[]{x0,x1,x2,x3},new double[]{y0,y1,y2,y3},xv);
	}
	
	/**
     * Cubic polynomial interpolation y=f(x) given four data pairs
     * y0=f(x0), y1=f(x1), y2=f(x2) and y3=f(x3) interpolate for yv=f(xv).
     * 
     * This is suitable for the case that four data points are not evenly distributed.
     * 
     * xv within [x1,x2] for interpolation and else for extrapolation
     *
     * @param	x0		1st x points
     * @param	x1		2nd x points
     * @param	x2		3rd x points
     * @param	x3		4th x points
     * @param	y0		1st y points
     * @param	y1		2nd y points
     * @param	y2		3rd y points
     * @param	y3		4th y points
     * @param	xv		value of x to be evaluated
     *
     * @return	re		result of interpolation at xv
     */
	public static float cubicPolynomialInterpolation
	(float x0,float x1,float x2,float x3,float y0,float y1,float y2,float y3,float xv){
		float dx0=x1-x0,dx1=x2-x1,dx2=x3-x2,dx=xv-x1;
		float dy0=y1-y0,dy1=y2-y1,dy2=y3-y2;
		float ratio=dx/dx1;
		
		float z1=linearInterpolation(dy0/dx0,dy1/dx1,dx0/(dx0+dx1));
		float z2=linearInterpolation(dy1/dx1,dy2/dx2,dx1/(dx1+dx2));
		
		float d=y1;
		float c=z1;
		float b=3f*dy1-2f*z1*dx1-z2*dx1;
		float a=(z1+z2)*dx1-2f*dy1;
		
		float re=d;
		
		re+=c*dx;
		re+=b*ratio*ratio;
		re+=a*ratio*ratio*ratio;
		
		return re;
	}
	
	public static double cubicPolynomialInterpolation
	(double x0,double x1,double x2,double x3,double y0,double y1,double y2,double y3,double xv){
		double dx0=x1-x0,dx1=x2-x1,dx2=x3-x2,dx=xv-x1;
		double dy0=y1-y0,dy1=y2-y1,dy2=y3-y2;
		double ratio=dx/dx1;
		
		double z1=linearInterpolation(dy0/dx0,dy1/dx1,dx0/(dx0+dx1));
		double z2=linearInterpolation(dy1/dx1,dy2/dx2,dx1/(dx1+dx2));
		
		double d=y1;
		double c=z1;
		double b=3.0*dy1-2.0*z1*dx1-z2*dx1;
		double a=(z1+z2)*dx1-2.0*dy1;
		
		double re=d;
		
		re+=c*dx;
		re+=b*ratio*ratio;
		re+=a*ratio*ratio*ratio;
		
		return re;
	}
	
	/**
     * Lagrangian interpolation
     *
     * @param	x		an array of x points
     * @param	y		an array of y points
     * @param	xv		value within x array for interpolation, else for extrapolation
     *
     * @return	re		result of interpolation
     */
	public static float nLagrangeInterpolation(float[] x,float[] y,float xv){
		int n=x.length;	float re=0;
		
		if(n!=y.length) throw new IllegalArgumentException("Array lengths are not the same");
		
		for(int k=0;k<n;k++){
			float l=1;
			
			for(int i=0;i<n;i++)
			if(i!=k) l*=(xv-x[i])/(x[k]-x[i]);
			
			re+=l*y[k];
		}
		
		return re;
	}
	
	public static double nLagrangeInterpolation(double[] x,double[] y,double xv){
		int n=x.length;	double re=0;
		
		if(n!=y.length) throw new IllegalArgumentException("Array lengths are not the same");
		
		for(int k=0;k<n;k++){
			double l=1.0;
			
			for(int i=0;i<n;i++)
			if(i!=k) l*=(xv-x[i])/(x[k]-x[i]);
			
			re+=l*y[k];
		}
		
		return re;
	}
	
	
	/**
     * Interpolate an equally spaced data into a new one.  The length of the
     * new array can be shorter than that of the original one.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
     * 
     * Notice that if Type is not a periodic one, then the end-points of the new
     * array equal to those of the original one.
     * 
     * If Type is periodic, then the first points of both array are equal while the
     * last points are generally not.
     *
     * @param	src		array of data to be interpolated
     * @param	des		result after interpolation
     * @param	type	type of interpolation in enum Type
     */
	public static void interp1D(float[] src,float[] des,Type t){ interp1D(src,des,t,Float.NaN);}
	
	public static void interp1D(double[] src,double[] des,Type t){ interp1D(src,des,t,Double.NaN);}
	
	/**
     * Interpolate an equally spaced data into a new one.  The length of the
     * new array can be shorter than that of the original one.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
     * 
     * Notice that if Type is not a periodic one, then the end-points of the new
     * array equal to those of the original one.
     * 
     * If Type is periodic, then the first points of both array are equal while the
     * last points are generally not.
     *
     * @param	src		array of data to be interpolated
     * @param	des		result after interpolation
     * @param	type	type of interpolation in enum Type
     * @param	undef	undefined value
     */
	public static void interp1D(float[] src,float[] des,Type t,float undef){
		switch(t){
		case LINEAR : interp1DLinear(src,des,undef); break;
		case CUBIC_L: interp1DCubicL(src,des,undef); break;
		case CUBIC_P: interp1DCubicP(src,des,undef); break;
		case SPLINE : interp1DSpline(src,des,undef); break;
		case PERIODIC_LINEAR : interp1DPeriodicLinear(src,des,undef); break;
		case PERIODIC_CUBIC_L: interp1DPeriodicCubicL(src,des,undef); break;
		case PERIODIC_CUBIC_P: interp1DPeriodicCubicP(src,des,undef); break;
		case PERIODIC_SPLINE : interp1DPeriodicSpline(src,des,undef); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type for interp1D:\t"+t);
		}
	}
	
	public static void interp1D(double[] src,double[] des,Type t,double undef){
		switch(t){
		case LINEAR : interp1DLinear(src,des,undef); break;
		case CUBIC_L: interp1DCubicL(src,des,undef); break;
		case CUBIC_P: interp1DCubicP(src,des,undef); break;
		case SPLINE : interp1DSpline(src,des,undef); break;
		case PERIODIC_LINEAR : interp1DPeriodicLinear(src,des,undef); break;
		case PERIODIC_CUBIC_L: interp1DPeriodicCubicL(src,des,undef); break;
		case PERIODIC_CUBIC_P: interp1DPeriodicCubicP(src,des,undef); break;
		case PERIODIC_SPLINE : interp1DPeriodicSpline(src,des,undef); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type for interp1D:\t"+t);
		}
	}
	
	/**
     * Interpolate an equally spaced data into a new one.  The length of the
     * new array can be shorter than that of the original one.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
     * 
     * Notice that if Type is not a periodic one, then the end-points of the new
     * array equal to those of the original one.
     * 
     * If Type is periodic, then the first points of both array are equal while the
     * last points are generally not.
     *
     * @param	src		array of data to be interpolated
     * @param	T		length of the new array after interpolation
     * @param	type	type of interpolation in enum Type
     */
	public static float[] interp1D(float[] src,int T,Type t){ return interp1D(src,T,t,Float.NaN);}
	
	public static double[] interp1D(double[] src,int T,Type t){ return interp1D(src,T,t,Double.NaN);}
	
	/**
     * Interpolate an equally spaced data into a new one.  The length of the
     * new array can be shorter than that of the original one.
     * 
     * No x-coordinates as it is assumed that x are evenly distributed.
     * 
     * Notice that if Type is not a periodic one, then the end-points of the new
     * array equal to those of the original one.
     * 
     * If Type is periodic, then the first points of both array are equal while the
     * last points are generally not.
     *
     * @param	src		array of data to be interpolated
     * @param	T		length of the new array after interpolation
     * @param	type	type of interpolation in enum Type
     * @param	undef	undefined value
     */
	public static float[] interp1D(float[] src,int T,Type t,float undef){
		float[] des=new float[T];
		
		interp1D(src,des,t,undef);
		
		return des;
	}
	
	public static double[] interp1D(double[] src,int T,Type t,double undef){
		double[] des=new double[T];
		
		interp1D(src,des,t,undef);
		
		return des;
	}
	
	
	/**
     * Interpolate a 2D data into another one.
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     *
     * @param	src		array of data to be interpolated
     * @param	des		result after interpolation
     * @param	xType	type of interpolation in enum Type for src[j]
     * @param	yType	type of interpolation in enum Type for src
     */
	public static void interp2D(float[][] src,float[][] des,Type xType,Type yType){
		interp2D(src,des,xType,yType,Float.NaN);
	}
	
	public static void interp2D(double[][] src,double[][] des,Type xType,Type yType){
		interp2D(src,des,xType,yType,Double.NaN);
	}
	
	public static float[][] interp2D(float[][] src,int x,int y,Type xType,Type yType){
		float[][] des=new float[y][x];
		
		interp2D(src,des,xType,yType);
		
		return des;
	}
	
	public static double[][] interp2D(double[][] src,int x,int y,Type xType,Type yType){
		double[][] des=new double[y][x];
		
		interp2D(src,des,xType,yType);
		
		return des;
	}
	
	
	/**
     * interpolate a 2D data into another one.
     * 
     * No x and y coordinates as it is assumed that they are evenly distributed.
     *
     * @param	src		array of data to be interpolated
     * @param	des		result after interpolation
     * @param	xType	type of interpolation in enum Type for src[j]
     * @param	yType	type of interpolation in enum Type for src
     * @param	undef	undefined value
     */
	public static void interp2D(float[][] src,float[][] des,Type xType,Type yType,float undef){
		int sy=src.length,dy=des.length,dx=des[0].length;
		
		float[]   tmp1=new float[sy];
		float[]   tmp2=new float[dy];
		float[][] buff=new float[sy][dx];
		
		for(int j=0;j<sy;j++) interp1D(src[j],buff[j],xType,undef);
		
		for(int i=0;i<dx;i++){
			for(int j=0;j<sy;j++) tmp1[j]=buff[j][i];
			interp1D(tmp1,tmp2,yType,undef);
			for(int j=0;j<dy;j++) des[j][i]=tmp2[j];
		}
	}
	
	public static void interp2D(double[][] src,double[][] des,Type xType,Type yType,double undef){
		int sy=src.length,dy=des.length,dx=des[0].length;
		
		double[]   tmp1=new double[sy];
		double[]   tmp2=new double[dy];
		double[][] buff=new double[sy][dx];
		
		for(int j=0;j<sy;j++) interp1D(src[j],buff[j],xType,undef);
		
		for(int i=0;i<dx;i++){
			for(int j=0;j<sy;j++) tmp1[j]=buff[j][i];
			interp1D(tmp1,tmp2,yType,undef);
			for(int j=0;j<dy;j++) des[j][i]=tmp2[j];
		}
	}
	
	public static float[][] interp2D(float[][] src,int x,int y,Type xType,Type yType,float undef){
		float[][] des=new float[y][x];
		
		interp2D(src,des,xType,yType,undef);
		
		return des;
	}
	
	public static double[][] interp2D(double[][] src,int x,int y,Type xType,Type yType,double undef){
		double[][] des=new double[y][x];
		
		interp2D(src,des,xType,yType,undef);
		
		return des;
	}
	
	
	/**
     * Interpolate a data (probably not evenly distributed) specified as sy=f(sx) into a new one.
     * The length of the new array can be shorter than that of the original one.
     * 
     * Notice that Type can not be periodic-like as data may be unevenly distributed.
     *
     * @param	sx		array of data to be interpolated
     * @param	sy		array of data to be interpolated
     * @param	des		result after interpolation
     * @param	type	type of interpolation in enum Type
     */
	public static void interp1D(float[] sx,float[] sy,float[] dx,float[] dy,Type t){ interp1D(sx,sy,dx,dy,t,Float.NaN);}
	
	public static void interp1D(double[] sx,double[] sy,double[] dx,double[] dy,Type t){ interp1D(sx,sy,dx,dy,t,Double.NaN);}
	
	/**
     * Interpolate a data (probably not evenly distributed) specified as sy=f(sx) into a new one.
     * The length of the new array can be shorter than that of the original one.
     * 
     * Notice that Type can not be periodic-like as data may be unevenly distributed.
     *
     * @param	sx		array of data to be interpolated
     * @param	sy		array of data to be interpolated
     * @param	des		result after interpolation
     * @param	type	type of interpolation in enum Type
     * @param	undef	undefined value
     */
	public static void interp1D(float[] sx,float[] sy,float[] dx,float[] dy,Type t,float undef){
		switch(t){
		case LINEAR : interp1DLinearU(sx,sy,dx,dy,undef); break;
		case CUBIC_L: interp1DCubicLU(sx,sy,dx,dy,undef); break;
		case CUBIC_P: interp1DCubicPU(sx,sy,dx,dy,undef); break;
		case SPLINE : interp1DSplineU(sx,sy,dx,dy,undef); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type for interp1D:\t"+t);
		}
	}
	
	public static void interp1D(double[] sx,double[] sy,double[] dx,double[] dy,Type t,double undef){
		switch(t){
		case LINEAR : interp1DLinearU(sx,sy,dx,dy,undef); break;
		case CUBIC_L: interp1DCubicLU(sx,sy,dx,dy,undef); break;
		case CUBIC_P: interp1DCubicPU(sx,sy,dx,dy,undef); break;
		case SPLINE : interp1DSplineU(sx,sy,dx,dy,undef); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type for interp1D:\t"+t);
		}
	}
	
	
	/**
     * Interpolate a data (probably not evenly distributed) specified as sy=f(sx) into a new one.
     * The length of the new array can be shorter than that of the original one.
     * 
     * Notice that Type can not be periodic-like as data may be unevenly distributed.
     *
     * @param	sx		array of data to be interpolated
     * @param	sy		array of data to be interpolated
     * @param	dx		x-points needs to be evaluated
     * @param	t		type of interpolation in enum Type
     * 
     * @return	dy		results of the interpolation
     */
	public static float[] interp1D(float[] sx,float[] sy,float[] dx,Type t){
		float[] dy=new float[dx.length];
		interp1D(sx,sy,dx,dy,t,Float.NaN);
		return dy;
	}
	
	public static double[] interp1D(double[] sx,double[] sy,double[] dx,Type t){
		double[] dy=new double[dx.length];
		interp1D(sx,sy,dx,dy,t,Double.NaN);
		return dy;
	}
	
	/**
     * Interpolate a data (probably not evenly distributed) specified as sy=f(sx) into a new one.
     * The length of the new array can be shorter than that of the original one.
     * 
     * Notice that Type can not be periodic-like as data may be unevenly distributed.
     *
     * @param	sx		array of data to be interpolated
     * @param	sy		array of data to be interpolated
     * @param	dx		x-points needs to be evaluated
     * @param	t		type of interpolation in enum Type
     * @param	undef	undefined value
     * 
     * @param	dy		result after interpolation
     */
	public static float[] interp1D(float[] sx,float[] sy,float[] dx,Type t,float undef){
		float[] dy=new float[dx.length];
		
		switch(t){
		case LINEAR : interp1DLinearU(sx,sy,dx,dy,undef); break;
		case CUBIC_L: interp1DCubicLU(sx,sy,dx,dy,undef); break;
		case CUBIC_P: interp1DCubicPU(sx,sy,dx,dy,undef); break;
		case SPLINE : interp1DSplineU(sx,sy,dx,dy,undef); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type for interp1D:\t"+t);
		}
		
		return dy;
	}
	
	public static double[] interp1D(double[] sx,double[] sy,double[] dx,Type t,double undef){
		double[] dy=new double[dx.length];
		
		switch(t){
		case LINEAR : interp1DLinearU(sx,sy,dx,dy,undef); break;
		case CUBIC_L: interp1DCubicLU(sx,sy,dx,dy,undef); break;
		case CUBIC_P: interp1DCubicPU(sx,sy,dx,dy,undef); break;
		case SPLINE : interp1DSplineU(sx,sy,dx,dy,undef); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type for interp1D:\t"+t);
		}
		
		return dy;
	}
	
	
	/**
     * fill the undefined values in a given periodic array using linear interpolation
     *
     * @param	src		array of data to be interpolated
     * @param	undef	undefined value
     */
	public static void fillUndefInPeriodicData(float[] src,Type t,float undef){
		switch(t){
		case LINEAR: linearFillUndef(src,undef,false); break;
		case SPLINE: splineFillUndef(src,undef,false); break;
		case PERIODIC_LINEAR: linearFillUndef(src,undef,true); break;
		case PERIODIC_SPLINE: splineFillUndef(src,undef,true); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type:\t"+t);
		}
	}
	
	public static void fillUndefInPeriodicData(double[] src,Type t,double undef){
		switch(t){
		case LINEAR: linearFillUndef(src,undef,false); break;
		case SPLINE: splineFillUndef(src,undef,false); break;
		case PERIODIC_LINEAR: linearFillUndef(src,undef,true); break;
		case PERIODIC_SPLINE: splineFillUndef(src,undef,true); break;
		default: throw new UnsupportedOperationException("unsupported interpolation type:\t"+t);
		}
	}
	
	
	/**
     * cubic spline interpolation with first boundary condition
     *
     * @param	sx	source x points
     * @param	sy	source y value, sy=f(sx[i])
     * @param	dx	destined x point
     *
     * @return	r	result of the given x point, r[i]=f(dx[i])
     */
	public static float[] cubicSplineWith1stBC(float[] sx,float[] sy,float[] dx,float c1,float c2){
		Spline spln=new Spline(sx,sy);
		
		spln.cubicSplineWith1stBC(c1,c2);
		
		return spln.cValues(dx);
	}
	
	public static double[] cubicSplineWith1stBC(double[] sx,double[] sy,double[] dx,double c1,double c2){
		Spline spln=new Spline(sx,sy);
		
		spln.cubicSplineWith1stBC(c1,c2);
		
		return spln.cValues(dx);
	}
	
	/**
     * cubic spline interpolation with first boundary condition
     *
     * @param	sx	source x points
     * @param	sy	source y value, sy[i]=f(sx[i])
     * @param	dx	destined x point
     *
     * @return	r	result of the given x point, r[i]=f(dx[i])
     */
	public static float[] cubicSplineWith2ndBC(float[] sx,float[] sy,float[] dx,float c1,float c2){
		Spline spln=new Spline(sx,sy);
		
		spln.cubicSplineWith2ndBC(c1,c2);
		
		return spln.cValues(dx);
	}
	
	public static double[] cubicSplineWith2ndBC(double[] sx,double[] sy,double[] dx,double c1,double c2){
		Spline spln=new Spline(sx,sy);
		
		spln.cubicSplineWith2ndBC(c1,c2);
		
		return spln.cValues(dx);
	}
	
	
	/**
     * Linearly interpolate the monthly climatology to daily climatology using Killworth method
     * 
     * Reference: Killworth 1996, JPO
     *
     * @param	monthly		monthly climatology
     *
     * @return	daily		daily climatology
     */
	public static float[] monthlyClimToDailyKillworth(float[] monthly){
		float[] buf=new float[366];
		
		interp1D(FilterModel.KillworthTransform(monthly),buf,Type.PERIODIC_LINEAR);
		
		float[] daily=new float[365];
		
		// offset half month
		for(int i=15;i<365;i++) daily[i]=buf[i-15];
		for(int i=0 ;i<15 ;i++) daily[i]=buf[350+i];
		
		return daily;
	}
	
	/**
     * Linearly interpolate the monthly climatology to daily climatology
     *
     * @param	monthly		monthly climatology
     *
     * @return	daily		daily climatology
     */
	public static float[] monthlyClimToDailyLinearly(float[] monthly){
		float[] buf=new float[366];
		
		interp1D(monthly,buf,Type.PERIODIC_LINEAR);
		
		float[] daily=new float[365];
		
		// offset half month
		for(int i=15;i<365;i++) daily[i]=buf[i-15];
		for(int i=0 ;i<15 ;i++) daily[i]=buf[350+i];
		
		return daily;
	}
	
	/**
     * Stepwisely interpolate the monthly climatology to daily climatology
     *
     * @param	monthly		monthly climatology
     *
     * @return	daily		daily climatology
     */
	public static float[] monthlyClimToDailyStepwise(float[] monthly){
		float[] daily=new float[365];
		
		for(int i=0  ;i<31 ;i++) daily[i]=monthly[0];
		for(int i=31 ;i<59 ;i++) daily[i]=monthly[1];
		for(int i=59 ;i<90 ;i++) daily[i]=monthly[2];
		for(int i=90 ;i<120;i++) daily[i]=monthly[3];
		for(int i=120;i<151;i++) daily[i]=monthly[4];
		for(int i=151;i<181;i++) daily[i]=monthly[5];
		for(int i=181;i<212;i++) daily[i]=monthly[6];
		for(int i=212;i<243;i++) daily[i]=monthly[7];
		for(int i=243;i<273;i++) daily[i]=monthly[8];
		for(int i=273;i<304;i++) daily[i]=monthly[9];
		for(int i=304;i<334;i++) daily[i]=monthly[10];
		for(int i=334;i<365;i++) daily[i]=monthly[11];
		
		return daily;
	}
	
	
	/*** helper methods ***/
	private static void linearFillUndef(float[] src,float undef,boolean periodic){
		int x=src.length;
		
		List<int[]> lst=new ArrayList<int[]>();
		
		int str=0,end=0;
		boolean setStr=false,setEnd=false;
		
		// find piecewise undefined values
		for(int i=0,I=x-1;i<I;i++){
			if(setEnd==false&&src[i]!=undef&&src[i+1]==undef){ str=i  ; setStr=true;}
			if(setStr==true &&src[i]==undef&&src[i+1]!=undef){ end=i+1; setEnd=true;}
			
			if(setStr&&setEnd){
				lst.add(new int[]{str,end});
				
				setStr=false; setEnd=false;
			}
		}
		
		// piecewise linear interpolation
		for(int[] tags:lst){
			int len=tags[1]-tags[0]+1;
			
			for(int i=tags[0]+1;i<tags[1];i++)
			src[i]=linearInterpolation(src[tags[0]],src[tags[1]],(float)(i-tags[0])/(len-1));
		}
		
		if(periodic){
			str=-1;	end=x;
			
			// process the undefined values between end-points
			for(int i=0;i<x;i++){
				if(str==-1&&src[i    ]!=undef) str=i;
				if(end== x&&src[x-i-1]!=undef) end=x-i-1;
			}
			
			float[] buf=new float[str+x-end+1];
			
			for(int i=0,I=buf.length;i<I;i++){
				int tag=end+i;
				
				if(tag>=x) tag-=x;
				
				buf[i]=src[tag];
			}
			
			for(int i=1,I=buf.length-1;i<I;i++)
			buf[i]=linearInterpolation(buf[0],buf[I],(float)i/I);
			
			for(int i=0,I=buf.length;i<I;i++){
				int tag=end+i;
				
				if(tag>=x) tag-=x;
				
				src[tag]=buf[i];
			}
		}
	}
	
	private static void linearFillUndef(double[] src,double undef,boolean periodic){
		int x=src.length;
		
		List<int[]> lst=new ArrayList<int[]>();
		
		int str=0,end=0;
		boolean setStr=false,setEnd=false;
		
		// find piecewise undefined values
		for(int i=0,I=x-1;i<I;i++){
			if(setEnd==false&&src[i]!=undef&&src[i+1]==undef){ str=i  ; setStr=true;}
			if(setStr==true &&src[i]==undef&&src[i+1]!=undef){ end=i+1; setEnd=true;}
			
			if(setStr&&setEnd){
				lst.add(new int[]{str,end});
				
				setStr=false; setEnd=false;
			}
		}
		
		// piecewise linear interpolation
		for(int[] tags:lst){
			int len=tags[1]-tags[0]+1;
			
			for(int i=tags[0]+1;i<tags[1];i++)
			src[i]=linearInterpolation(src[tags[0]],src[tags[1]],(float)(i-tags[0])/(len-1));
		}
		
		if(periodic){
			str=-1;	end=x;
			
			// process the undefined values between end-points
			for(int i=0;i<x;i++){
				if(str==-1&&src[i    ]!=undef) str=i;
				if(end== x&&src[x-i-1]!=undef) end=x-i-1;
			}
			
			double[] buf=new double[str+x-end+1];
			
			for(int i=0,I=buf.length;i<I;i++){
				int tag=end+i;
				
				if(tag>=x) tag-=x;
				
				buf[i]=src[tag];
			}
			
			for(int i=1,I=buf.length-1;i<I;i++)
			buf[i]=linearInterpolation(buf[0],buf[I],(double)i/I);
			
			for(int i=0,I=buf.length;i<I;i++){
				int tag=end+i;
				
				if(tag>=x) tag-=x;
				
				src[tag]=buf[i];
			}
		}
	}
	
	
	private static void splineFillUndef(float[] src,float undef,boolean periodic){
		float[] buf=new float[src.length];
		
		int shift=0;
		
		if(periodic) shift=shiftToStartAndEndWithDef(src,buf,undef);
		
		int x=buf.length;
		int defCount=0;
		
		float[] sx=new float[x];
		
		for(int i=0;i<x;i++){
			sx[i]=i;
			
			if(buf[i]!=undef) defCount++;
		}
		
		float[] tmpx=new float[defCount];
		float[] tmpy=new float[defCount];
		
		for(int i=0,ptr=0;i<x;i++) if(buf[i]!=undef){
			tmpx[ptr]=i;
			tmpy[ptr]=buf[i];
			ptr++;
		}
		
		Spline spln=new Spline(tmpx,tmpy);
		
		if(periodic)
			spln.cubicSplineWithPeriodicBC();
		else
			spln.cubicSplineWithNaturalBC();
		
		spln.cValues(sx,buf);
		
		for(int i=0;i<x;i++){
			int tag=i-shift;
			
			if(tag<0) tag+=x;
			
			src[i]=buf[tag];
		}
	}
	
	private static void splineFillUndef(double[] src,double undef,boolean periodic){
		double[] buf=new double[src.length];
		
		int shift=0;
		
		if(periodic) shift=shiftToStartAndEndWithDef(src,buf,undef);
		
		int x=buf.length;
		int defCount=0;
		
		double[] sx=new double[x];
		
		for(int i=0;i<x;i++){
			sx[i]=i;
			
			if(buf[i]!=undef) defCount++;
		}
		
		double[] tmpx=new double[defCount];
		double[] tmpy=new double[defCount];
		
		for(int i=0,ptr=0;i<x;i++) if(buf[i]!=undef){
			tmpx[ptr]=i;
			tmpy[ptr]=buf[i];
			ptr++;
		}
		
		Spline spln=new Spline(tmpx,tmpy);
		
		if(periodic)
			spln.cubicSplineWithPeriodicBC();
		else
			spln.cubicSplineWithNaturalBC();
		
		spln.cValues(sx,buf);
		
		for(int i=0;i<x;i++){
			int tag=i-shift;
			
			if(tag<0) tag+=x;
			
			src[i]=buf[tag];
		}
	}
	
	
	private static int shiftToStartAndEndWithDef(float[] src,float[] re,float undef){
		int x=src.length;
		
		int str=0;
		
		if(src[0]!=undef&&src[x-1]!=undef){
			System.arraycopy(src,0,re,0,x);
			
			return str;
		}
		
		boolean match=false;
		for(int i=0,I=x-1;i<I;i++)
		if(src[i]!=undef&&src[i+1]!=undef){ str=i+1; match=true; break;}
		
		if(!match)
		throw new IllegalArgumentException("no at least two consecutive defined values in the array");
		
		for(int i=0;i<x;i++){
			int tag=str+i;
			
			if(tag>=x) tag-=x;
			
			re[i]=src[tag];
		}
		
		return str;
	}
	
	private static int shiftToStartAndEndWithDef(double[] src,double[] re,double undef){
		int x=src.length;
		
		int str=0;
		
		if(src[0]!=undef&&src[x-1]!=undef){
			System.arraycopy(src,0,re,0,x);
			
			return str;
		}
		
		boolean match=false;
		for(int i=0,I=x-1;i<I;i++)
		if(src[i]!=undef&&src[i+1]!=undef){ str=i+1; match=true; break;}
		
		if(!match)
		throw new IllegalArgumentException("no at least two consecutive defined values in the array");
		
		for(int i=0;i<x;i++){
			int tag=str+i;
			
			if(tag>=x) tag-=x;
			
			re[i]=src[tag];
		}
		
		return str;
	}
	
	
	private static void interp1DLinear(float[] src,float[] des,float undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<2) throw new IllegalArgumentException(
			"source array is too short (at least 2) for linear interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for linear interpolation"
		);
		
		float ratio=(float)(slen-1)/(dlen-1);
		
		des[0]=src[0];	des[dlen-1]=src[slen-1];
		
		for(int i=1,I=dlen-1;i<I;i++){
			float dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Float.isNaN(undef)?
			linearInterpolation(src[tag],src[tag+1],dis):
			linearInterpolation(src[tag],src[tag+1],dis,undef);
		}
	}
	
	private static void interp1DLinear(double[] src,double[] des,double undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<2) throw new IllegalArgumentException(
			"source array is too short (at least 2) for linear interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for linear interpolation"
		);
		
		double ratio=(double)(slen-1.0)/(dlen-1.0);
		
		des[0]=src[0];	des[dlen-1]=src[slen-1];
		
		for(int i=1,I=dlen-1;i<I;i++){
			double dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Double.isNaN(undef)?
			linearInterpolation(src[tag],src[tag+1],dis):
			linearInterpolation(src[tag],src[tag+1],dis,undef);
		}
	}
	
	
	private static void interp1DCubicL(float[] src,float[] des,float undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for cubic Lagrange interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for cubic Lagrange interpolation"
		);
		
		float ratio=(float)(slen-1)/(dlen-1);
		
		des[0]=src[0];	des[dlen-1]=src[slen-1];
		
		// count for points pairs located between
		// src[0] and src[1], or src[slen-2] and src[slen-1]
		int secount=(dlen-1)/(slen-1);
		
		for(int i=1;i<=secount;i++)
		if(Float.isNaN(undef)){
			des[i]=
			cubicLagrangeInterpolation(src[0],src[1],src[2],src[3],ratio*i);
			des[dlen-1-i]=
			cubicLagrangeInterpolation(src[slen-4],src[slen-3],src[slen-2],src[slen-1],3-ratio*i);
			
		}else{
			des[i]=
			cubicLagrangeInterpolation(src[0],src[1],src[2],src[3],ratio*i,undef);
			des[dlen-1-i]=
			cubicLagrangeInterpolation(src[slen-4],src[slen-3],src[slen-2],src[slen-1],3-ratio*i,undef);
		}
		
		for(int i=secount+1;i<dlen-secount-1;i++){
			float dis=i*ratio;  int tag=(int)dis;  dis-=tag-1;
			des[i]=Float.isNaN(undef)?
			cubicLagrangeInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis):
			cubicLagrangeInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis,undef);
		}
	}
	
	private static void interp1DCubicL(double[] src,double[] des,double undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for cubic Lagrange interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for cubic Lagrange interpolation"
		);
		
		double ratio=(double)(slen-1.0)/(dlen-1.0);
		
		des[0]=src[0];	des[dlen-1]=src[slen-1];
		
		// count for points pairs located between
		// src[0] and src[1], or src[slen-2] and src[slen-1]
		int secount=(dlen-1)/(slen-1);
		
		for(int i=1;i<=secount;i++)
		if(Double.isNaN(undef)){
			des[i]=
			cubicLagrangeInterpolation(src[0],src[1],src[2],src[3],ratio*i);
			des[dlen-1-i]=
			cubicLagrangeInterpolation(src[slen-4],src[slen-3],src[slen-2],src[slen-1],3-ratio*i);
			
		}else{
			des[i]=
			cubicLagrangeInterpolation(src[0],src[1],src[2],src[3],ratio*i,undef);
			des[dlen-1-i]=
			cubicLagrangeInterpolation(src[slen-4],src[slen-3],src[slen-2],src[slen-1],3-ratio*i,undef);
		}
		
		for(int i=secount+1;i<dlen-secount-1;i++){
			double dis=i*ratio;  int tag=(int)dis;  dis-=tag-1;
			des[i]=Double.isNaN(undef)?
			cubicLagrangeInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis):
			cubicLagrangeInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis,undef);
		}
	}
	
	
	private static void interp1DCubicP(float[] src,float[] des,float undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for cubic polynomial interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for cubic polynomial interpolation"
		);
		
		float ratio=(float)(slen-1)/(dlen-1);
		
		des[0]=src[0];	des[dlen-1]=src[slen-1];
		
		// count for points pairs located between
		// src[0] and src[1], or src[slen-2] and src[slen-1]
		int secount=(dlen-1)/(slen-1);
		
		for(int i=1;i<=secount;i++)
		if(Float.isNaN(undef)){
			des[i]=
			cubicPolynomialInterpolation(src[0],src[0],src[1],src[2],ratio*i);
			des[dlen-1-i]=
			cubicPolynomialInterpolation(src[slen-3],src[slen-2],src[slen-1],src[slen-1],1-ratio*i);
			
		}else{
			des[i]=
			cubicPolynomialInterpolation(src[0],src[0],src[1],src[2],ratio*i,undef);
			des[dlen-1-i]=
			cubicPolynomialInterpolation(src[slen-3],src[slen-2],src[slen-1],src[slen-1],1-ratio*i,undef);
		}
		
		for(int i=secount+1,I=dlen-secount-1;i<I;i++){
			float dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Float.isNaN(undef)?
			cubicPolynomialInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis):
			cubicPolynomialInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis,undef);
		}
	}
	
	private static void interp1DCubicP(double[] src,double[] des,double undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for cubic polynomial interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for cubic polynomial interpolation"
		);
		
		double ratio=(double)(slen-1.0)/(dlen-1.0);
		
		des[0]=src[0];	des[dlen-1]=src[slen-1];
		
		// count for points pairs located between
		// src[0] and src[1], or src[slen-2] and src[slen-1]
		int secount=(dlen-1)/(slen-1);
		
		for(int i=1;i<=secount;i++)
		if(Double.isNaN(undef)){
			des[i]=
			cubicPolynomialInterpolation(src[0],src[0],src[1],src[2],ratio*i);
			des[dlen-1-i]=
			cubicPolynomialInterpolation(src[slen-3],src[slen-2],src[slen-1],src[slen-1],1-ratio*i);
			
		}else{
			des[i]=
			cubicPolynomialInterpolation(src[0],src[0],src[1],src[2],ratio*i,undef);
			des[dlen-1-i]=
			cubicPolynomialInterpolation(src[slen-3],src[slen-2],src[slen-1],src[slen-1],1-ratio*i,undef);
		}
		
		for(int i=secount+1,I=dlen-secount-1;i<I;i++){
			double dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Double.isNaN(undef)?
			cubicPolynomialInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis):
			cubicPolynomialInterpolation(src[tag-1],src[tag],src[tag+1],src[tag+2],dis,undef);
		}
	}
	
	
	private static void interp1DSpline(float[] src,float[] des,float undef){
		if(!Float.isNaN(undef))
		throw new IllegalArgumentException("unsupported spline interpolation with undefined value");
		
		int slen=src.length,dlen=des.length;
		
		float ratio=(float)(slen-1)/(dlen-1);
		
		float[] sx=new float[slen];
		float[] dx=new float[dlen];
		
		for(int i=0;i<slen;i++) sx[i]=i;
		for(int i=0;i<dlen;i++) dx[i]=i*ratio;
		
		Spline s=new Spline(sx,src);
		s.cubicSplineWithNaturalBC();
		s.cValues(dx,des);
	}
	
	private static void interp1DSpline(double[] src,double[] des,double undef){
		if(!Double.isNaN(undef))
		throw new IllegalArgumentException("unsupported spline interpolation with undefined value");
		
		int slen=src.length,dlen=des.length;
		
		double ratio=(double)(slen-1.0)/(dlen-1.0);
		
		double[] sx=new double[slen];
		double[] dx=new double[dlen];
		
		for(int i=0;i<slen;i++) sx[i]=i;
		for(int i=0;i<dlen;i++) dx[i]=i*ratio;
		
		Spline s=new Spline(sx,src);
		s.cubicSplineWithNaturalBC();
		s.cValues(dx,des);
	}
	
	
	private static void interp1DPeriodicLinear(float[] src,float[] des,float undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<2) throw new IllegalArgumentException(
			"source array is too short (at least 2) for periodic linear interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for periodic linear interpolation"
		);
		
		float ratio=(float)slen/dlen;
		
		des[0]=src[0];
		
		for(int i=1;i<dlen;i++){
			float dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Float.isNaN(undef)?
			linearInterpolation(src[tag],tag+1>=slen?src[0]:src[tag+1],dis):
			linearInterpolation(src[tag],tag+1>=slen?src[0]:src[tag+1],dis,undef);
		}
	}
	
	private static void interp1DPeriodicLinear(double[] src,double[] des,double undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<2) throw new IllegalArgumentException(
			"source array is too short (at least 2) for periodic linear interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for periodic linear interpolation"
		);
		
		double ratio=(double)slen/dlen;
		
		des[0]=src[0];
		
		for(int i=1;i<dlen;i++){
			double dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Double.isNaN(undef)?
			linearInterpolation(src[tag],tag+1>=slen?src[0]:src[tag+1],dis):
			linearInterpolation(src[tag],tag+1>=slen?src[0]:src[tag+1],dis,undef);
		}
	}
	
	
	private static void interp1DPeriodicCubicL(float[] src,float[] des,float undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for periodic cubic Lagrange interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for periodic cubic Lagrange interpolation"
		);
		
		float ratio=(float)slen/dlen;
		
		des[0]=src[0];
		
		for(int i=1;i<dlen;i++){
			float dis=i*ratio;  int tag=(int)dis;  dis-=tag-1;
			
			des[i]=Float.isNaN(undef)?
			cubicLagrangeInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis
			):
			cubicLagrangeInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis,undef
			);
		}
	}
	
	private static void interp1DPeriodicCubicL(double[] src,double[] des,double undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for periodic cubic Lagrange interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for periodic cubic Lagrange interpolation"
		);
		
		float ratio=(float)slen/dlen;
		
		des[0]=src[0];
		
		for(int i=1;i<dlen;i++){
			double dis=i*ratio;  int tag=(int)dis;  dis-=tag-1;
			
			des[i]=Double.isNaN(undef)?
			cubicLagrangeInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis
			):
			cubicLagrangeInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis,undef
			);
		}
	}
	
	
	private static void interp1DPeriodicCubicP(float[] src,float[] des,float undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for periodic cubic polynomial interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for periodic cubic polynomial interpolation"
		);
		
		float ratio=(float)slen/dlen;
		
		des[0]=src[0];
		
		for(int i=1;i<dlen;i++){
			float dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Float.isNaN(undef)?
			cubicPolynomialInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis
			):
			cubicPolynomialInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis,undef
			);
		}
	}
	
	private static void interp1DPeriodicCubicP(double[] src,double[] des,double undef){
		int slen=src.length,dlen=des.length;
		
		if(slen<4) throw new IllegalArgumentException(
			"source array is too short (at least 4) for periodic cubic polynomial interpolation"
		);
		
		if(dlen<2) throw new IllegalArgumentException(
			"destinated array is too short (at least 2) for periodic cubic polynomial interpolation"
		);
		
		double ratio=(double)slen/dlen;
		
		des[0]=src[0];
		
		for(int i=1;i<dlen;i++){
			double dis=i*ratio;  int tag=(int)dis;  dis-=tag;
			des[i]=Double.isNaN(undef)?
			cubicPolynomialInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis
			):
			cubicPolynomialInterpolation(
				tag-1<0?src[slen-1]:src[tag-1],src[tag],
				tag+1>=slen?src[tag+1-slen]:src[tag+1],
				tag+2>=slen?src[tag+2-slen]:src[tag+2],dis,undef
			);
		}
	}
	
	
	private static void interp1DPeriodicSpline(float[] src,float[] des,float undef){
		if(!Float.isNaN(undef))
		throw new IllegalArgumentException("unsupported periodic spline interpolation with undefined value");
		
		float[] extSrc=extendPeriodicData(src);
		
		int slen=extSrc.length,dlen=des.length;
		
		float ratio=(float)(slen-1)/dlen;
		
		float[] sx=new float[slen];
		float[] dx=new float[dlen];
		
		for(int i=0;i<slen;i++) sx[i]=i;
		for(int i=0;i<dlen;i++) dx[i]=i*ratio;
		
		Spline s=new Spline(sx,extSrc);
		s.cubicSplineWithPeriodicBC();
		s.cValues(dx,des);
	}
	
	private static void interp1DPeriodicSpline(double[] src,double[] des,double undef){
		if(!Double.isNaN(undef))
		throw new IllegalArgumentException("unsupported periodic spline interpolation with undefined value");
		
		double[] extSrc=extendPeriodicData(src);
		
		int slen=extSrc.length,dlen=des.length;
		
		double ratio=(double)(slen-1.0)/dlen;
		
		double[] sx=new double[slen];
		double[] dx=new double[dlen];
		
		for(int i=0;i<slen;i++) sx[i]=i;
		for(int i=0;i<dlen;i++) dx[i]=i*ratio;
		
		Spline s=new Spline(sx,extSrc);
		s.cubicSplineWithPeriodicBC();
		s.cValues(dx,des);
	}
	
	
	private static void interp1DLinearU(float[] sx,float[] sy,float[] dx,float[] dy,float undef){
		int slen=sx.length;
		
		boolean incre=sx[slen-1]-sx[0]>0;
		
		for(int i=0,I=dx.length;i<I;i++){
			float x=dx[i];
			
			if(incre){ if(x>sx[slen-1]||x<sx[0]){ dy[i]=undef; continue;}}
			else{      if(x>sx[0]||x<sx[slen-1]){ dy[i]=undef; continue;}}
			
			int idx=0;
			if(incre) idx=ArrayUtil.getLEIdxIncre(sx,x);
			else      idx=ArrayUtil.getLEIdxDecre(sx,x);
			
			if(x==sx[idx]){ dy[i]=sy[idx]; continue;}
			
			dy[i]=linearInterpolation(sy[idx],sy[idx+1],(x-sx[idx])/(sx[idx+1]-sx[idx]));
		}
	}
	
	private static void interp1DLinearU(double[] sx,double[] sy,double[] dx,double[] dy,double undef){
		int slen=sx.length;
		
		boolean incre=sx[slen-1]-sx[0]>0;
		
		for(int i=0,I=dx.length;i<I;i++){
			double x=dx[i];
			
			if(incre){ if(x>sx[slen-1]||x<sx[0]){ dy[i]=undef; continue;}}
			else{      if(x>sx[0]||x<sx[slen-1]){ dy[i]=undef; continue;}}
			
			int idx=0;
			if(incre) idx=ArrayUtil.getLEIdxIncre(sx,x);
			else      idx=ArrayUtil.getLEIdxDecre(sx,x);
			
			if(x==sx[idx]){ dy[i]=sy[idx]; continue;}
			
			dy[i]=linearInterpolation(sy[idx],sy[idx+1],(x-sx[idx])/(sx[idx+1]-sx[idx]));
		}
	}
	
	
	private static void interp1DCubicLU(float[] sx,float[] sy,float[] dx,float[] dy,float undef){
		int slen=sx.length;
		
		boolean incre=sx[slen-1]-sx[0]>0;
		
		for(int i=0,I=dx.length;i<I;i++){
			float x=dx[i];
			
			if(incre){ if(x>sx[slen-1]||x<sx[0]){ dy[i]=undef; continue;}}
			else{      if(x>sx[0]||x<sx[slen-1]){ dy[i]=undef; continue;}}
			
			int idx=0;
			if(incre) idx=ArrayUtil.getLEIdxIncre(sx,x);
			else      idx=ArrayUtil.getLEIdxDecre(sx,x);
			
			if(x==sx[idx]){ dy[i]=sy[idx]; continue;}
			
			if(idx==0){ dy[i]=nLagrangeInterpolation(
				new float[]{sx[0],sx[1],sx[2],sx[3]},
				new float[]{sy[0],sy[1],sy[2],sy[3]},x
			); continue;}
			
			if(idx==slen-2){ dy[i]=nLagrangeInterpolation(
				new float[]{sx[slen-4],sx[slen-3],sx[idx],sx[slen-1]},
				new float[]{sy[slen-4],sy[slen-3],sy[idx],sy[slen-1]},x
			); continue;}
			
			dy[i]=nLagrangeInterpolation(
				new float[]{sx[idx-1],sx[idx],sx[idx+1],sx[idx+2]},
				new float[]{sy[idx-1],sy[idx],sy[idx+1],sy[idx+2]},x
			);
		}
	}
	
	private static void interp1DCubicLU(double[] sx,double[] sy,double[] dx,double[] dy,double undef){
		int slen=sx.length;
		
		boolean incre=sx[slen-1]-sx[0]>0;
		
		for(int i=0,I=dx.length;i<I;i++){
			double x=dx[i];
			
			if(incre){ if(x>sx[slen-1]||x<sx[0]){ dy[i]=undef; continue;}}
			else{      if(x>sx[0]||x<sx[slen-1]){ dy[i]=undef; continue;}}
			
			int idx=0;
			if(incre) idx=ArrayUtil.getLEIdxIncre(sx,x);
			else      idx=ArrayUtil.getLEIdxDecre(sx,x);
			
			if(x==sx[idx]){ dy[i]=sy[idx]; continue;}
			
			if(idx==0){ dy[i]=nLagrangeInterpolation(
				new double[]{sx[0],sx[1],sx[2],sx[3]},
				new double[]{sy[0],sy[1],sy[2],sy[3]},x
			); continue;}
			
			if(idx==slen-2){ dy[i]=nLagrangeInterpolation(
				new double[]{sx[slen-4],sx[slen-3],sx[idx],sx[slen-1]},
				new double[]{sy[slen-4],sy[slen-3],sy[idx],sy[slen-1]},x
			); continue;}
			
			dy[i]=nLagrangeInterpolation(
				new double[]{sx[idx-1],sx[idx],sx[idx+1],sx[idx+2]},
				new double[]{sy[idx-1],sy[idx],sy[idx+1],sy[idx+2]},x
			);
		}
	}
	
	
	private static void interp1DCubicPU(float[] sx,float[] sy,float[] dx,float[] dy,float undef){
		int slen=sx.length;
		
		boolean incre=sx[slen-1]-sx[0]>0;
		
		for(int i=0,I=dx.length;i<I;i++){
			float x=dx[i];
			
			if(incre){ if(x>sx[slen-1]||x<sx[0]){ dy[i]=undef; continue;}}
			else{      if(x>sx[0]||x<sx[slen-1]){ dy[i]=undef; continue;}}
			
			int idx=0;
			if(incre) idx=ArrayUtil.getLEIdxIncre(sx,x);
			else      idx=ArrayUtil.getLEIdxDecre(sx,x);
			
			if(x==sx[idx]){ dy[i]=sy[idx]; continue;}
			
			if(idx==0){ dy[i]=cubicPolynomialInterpolation(
				sx[0],sx[0],sx[1],sx[2],
				sy[0],sy[0],sy[1],sy[2],x
			); continue;}
			
			if(idx==slen-2){ dy[i]=cubicPolynomialInterpolation(
				sx[slen-3],sx[idx],sx[slen-1],sx[slen-1],
				sy[slen-3],sy[idx],sy[slen-1],sy[slen-1],x
			); continue;}
			
			dy[i]=cubicPolynomialInterpolation(
				sx[idx-1],sx[idx],sx[idx+1],sx[idx+2],
				sy[idx-1],sy[idx],sy[idx+1],sy[idx+2],x
			);
		}
	}
	
	private static void interp1DCubicPU(double[] sx,double[] sy,double[] dx,double[] dy,double undef){
		int slen=sx.length;
		
		boolean incre=sx[slen-1]-sx[0]>0;
		
		for(int i=0,I=dx.length;i<I;i++){
			double x=dx[i];
			
			if(incre){ if(x>sx[slen-1]||x<sx[0]){ dy[i]=undef; continue;}}
			else{      if(x>sx[0]||x<sx[slen-1]){ dy[i]=undef; continue;}}
			
			int idx=0;
			if(incre) idx=ArrayUtil.getLEIdxIncre(sx,x);
			else      idx=ArrayUtil.getLEIdxDecre(sx,x);
			
			if(x==sx[idx]){ dy[i]=sy[idx]; continue;}
			
			if(idx==0){ dy[i]=cubicPolynomialInterpolation(
				sx[0],sx[0],sx[1],sx[2],
				sy[0],sy[0],sy[1],sy[2],x
			); continue;}
			
			if(idx==slen-2){ dy[i]=cubicPolynomialInterpolation(
				sx[slen-3],sx[idx],sx[slen-1],sx[slen-1],
				sy[slen-3],sy[idx],sy[slen-1],sy[slen-1],x
			); continue;}
			
			dy[i]=cubicPolynomialInterpolation(
				sx[idx-1],sx[idx],sx[idx+1],sx[idx+2],
				sy[idx-1],sy[idx],sy[idx+1],sy[idx+2],x
			);
		}
	}
	
	
	private static void interp1DSplineU(float[] sx,float[] sy,float[] dx,float[] dy,float undef){
		if(!Float.isNaN(undef))
		throw new IllegalArgumentException("unsupported spline interpolation with undefined value");
		
		Spline spl=new Spline(sx,sy);
		spl.cubicSplineWithNaturalBC();
		spl.cValues(dx,dy,undef);
	}
	
	private static void interp1DSplineU(double[] sx,double[] sy,double[] dx,double[] dy,double undef){
		if(!Double.isNaN(undef))
		throw new IllegalArgumentException("unsupported spline interpolation with undefined value");
		
		Spline spl=new Spline(sx,sy);
		spl.cubicSplineWithNaturalBC();
		spl.cValues(dx,dy,undef);
	}
	
	
	private static float[] extendPeriodicData(float[] src){
		float[] re=new float[src.length+1];
		
		System.arraycopy(src,0,re,0,src.length);
		
		re[src.length]=src[0];
		
		return re;
	}
	
	private static double[] extendPeriodicData(double[] src){
		double[] re=new double[src.length+1];
		
		System.arraycopy(src,0,re,0,src.length);
		
		re[src.length]=src[0];
		
		return re;
	}
	
	
	/** test
	public static void main(String[] arg){
		int len=401;
		float[] sx=new float[]{1,3,3.6f,4,8};
		float[] sy=new float[]{2,1,0,3,6};
		
		float[] dx=new float[len];
		
		for(int i=0;i<len;i++) dx[i]=3+i/400f;
		
		float[] des1=new float[len];
		float[] des2=new float[len];
		float[] des3=new float[len];
		float[] des4=new float[len];
		
		interp1D(sx,sy,dx,des1,Type.LINEAR);
		interp1D(sx,sy,dx,des2,Type.CUBIC_L);
		interp1D(sx,sy,dx,des3,Type.CUBIC_P);
		interp1D(sx,sy,dx,des4,Type.SPLINE);
		
		for(int i=0;i<len;i++){
			System.out.println(des1[i]+"\t"+des2[i]+"\t"+des3[i]+"\t"+des4[i]);
		}
	}*/
}
