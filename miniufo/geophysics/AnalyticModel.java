/**
 * @(#)AnalyticModel.java	1.0 2015.12.19
 * 
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.geophysics;

import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.tanh;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;
import static miniufo.diagnosis.SpatialModel.EARTH_ROTATE_SPEED;


/**
 * This class contains some analytic model used in geophysical fluid studies.
 *
 * @version 1.0, 2015.12.19
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AnalyticModel{
	
	/**
	 * prevent from instantiate
	 */
	private AnalyticModel(){}
	
	
	/**
     * Calculate a Bickley jet with standing waves on it.
     * 
     * Reference: Jayne and Marotzke (2002, JPO)
     *
     * @param	Samp	amplitude of streamfunction (m^2 s^-1)
     * @param	Wjet	half-width of the jet (degree)
     * @param	Mamp	meander amplitude (degree)
     * @param	Mwvl	meander wavelength (degree)
     * @param	Mprd	meander period (day)
     * @param	x		x-coordinate (degree)
     * @param	y		y-coordinate (degree)
     * @param	t		time (day)
     *
     * @return	streamfunction of the jet (m^2 s^-1)
     */
	public double cBickleyJetWithStandingWave(double Samp,double Wjet,double Mamp,double Mwvl,double Mprd,double x,double y,double t){
		return Samp*tanh((y/Wjet)+Mamp/Wjet*sin(2.0*PI*x/Mwvl)*sin(2.0*PI*t/Mprd));
	}
	
	/**
     * Calculate a Bickley jet with moving waves on it.
     * 
     * Reference: Jayne and Marotzke (2002, JPO)
     *
     * @param	Samp	amplitude of streamfunction (m^2 s^-1)
     * @param	Wjet	half-width of the jet (degree)
     * @param	Ldec	zonal decaying length of the wave field (degree)
     * @param	Mamp	meander amplitude (degree)
     * @param	Mwvl	meander wavelength (degree)
     * @param	Mprd	meander period (day)
     * @param	x		x-coordinate (degree)
     * @param	y		y-coordinate (degree)
     * @param	t		time (day)
     *
     * @return	streamfunction of the jet (m^2 s^-1)
     */
	public double cBickleyJetWithMovingWave(double Samp,double Wjet,double Ldec,double Mamp,double Mwvl,double Mprd,double x,double y,double t){
		return Samp*tanh((y/Wjet)+Mamp/Wjet*exp(-pow(x/Ldec,2.0))*sin(2.0*PI*x/Mwvl-2.0*PI*t/Mprd));
	}
	
	/**
     * Calculate Rossby-Haurwitz wave.
     * 
     * Reference: 
     *
     * @param	lamda	longitude (radian)
     * @param	fai		latitude (radian)
     * @param	Mamp	meander amplitude (degree)
     * @param	Mwvl	meander wavelength (degree)
     * @param	Mprd	meander period (day)
     * @param	x		x-coordinate (degree)
     * @param	y		y-coordinate (degree)
     * @param	t		time (day)
     *
     * @return	RH wave, [0] is geopotential (m^2 s^-2), [1] u- and [2] v- velocity components (m s^-1)
     */
	public double[] cRossbyHaurwitzWave(double lamda,double fai){
		double a=EARTH_RADIUS,o=7.848e-6,K=7.848e-6,fai0=9.8*8000.0;
		
		double A=o*(2.0*EARTH_ROTATE_SPEED+o)*pow(cos(fai),2.0)/2.0+K*K*pow(cos(fai),8.0)*(
			5.0*pow(cos(fai),2.0)+26.0-32.0*pow(cos(fai),-2.0)
		)/4.0;
		
		double B=(EARTH_ROTATE_SPEED+o)*K/15*pow(cos(fai),4.0)*(26.0-25.0*pow(cos(fai),2.0));
		
		double C=K*K*pow(cos(fai),8.0)*(5.0*pow(cos(fai),2.0)-6.0)/4.0;
		
		double h=fai0+(a*a*A+a*a*cos(4.0*lamda)*B+a*a*cos(8.0*lamda)*C);
		double u=a*o*cos(fai)+a*K*pow(cos(fai),3.0)*(4*pow(sin(fai),2.0)-pow(cos(fai),2.0))*cos(4.0*lamda);
		double v=-a*K*4.0*pow(cos(fai),3.0)*sin(fai)*sin(4.0*lamda);
		
		return new double[]{h,u,v};
	}
}
