/**
 * @(#)SphericalHarmonic.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import static java.lang.Math.PI;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;


/**
 * spherical harmonic class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class SphericalHarmonic{
	//
	private int l;
	private int m;
	
	private AssociatedLegendre al=null;
	
	
	/**
     * constructor
     *
     * @param	l	degree
     * @param	m	order
     */
	public SphericalHarmonic(){ this(0,0); al=new AssociatedLegendre();}
	
	public SphericalHarmonic(int l,int m){
		this.l=l; this.m=m;
		
		al=new AssociatedLegendre(l,Math.abs(m));
	}
	
	
	/*** getor and setor ***/
	public int getL(){ return l;}
	
	public int getM(){ return m;}
	
	public void setL(int l){ this.l=l; al.setL(l);}
	
	public void setM(int m){ this.m=m; al.setM(m);}
	
	public void setLM(int l,int m){ this.l=l; this.m=m; al.setLM(l,Math.abs(m));}
	
	
	/**
     * evaluate the spherical harmonic
     *
     * @param	lamda	longitude (in radians)
     * @param	fai		latitude  (in radians)
     */
	public Complex eval(float lamda,float fai){
		int absm=Math.abs(m);
		float sign=(absm%2==1)?-1:1;
		float factor=(float)(
			sign*sqrt((2*l+1)/(4.0*PI)*ff(l,absm))*al.evalDouble(sin(fai))
		);
		
		Complex retval=new Complex(0,absm*lamda).expEq().multiplyEq(factor);
		
		if(m<0) return retval.conjugateEq();
		
		return retval.multiply(sign);
    }
	
	
	/**
     * evaluate factorialDouble(l-m)/factorialDouble(l+m)
     *
     * @param	l	degree
     * @param	m	order
     */
	private static double ff(int l,int m){
		double re=1.0;
		
		for(int i=l-m+1,I=l+m;i<=I;i++) re*=i;
		
		if(re<=0) throw new IllegalArgumentException("result over flow");
		
		return 1.0/re;
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			SphericalHarmonic sh=new SphericalHarmonic();
			
			int l=10;
			float lamda=(float)Math.toRadians(20);
			float fai  =(float)Math.toRadians(10);
			
			for(int m=-10;m<=0;m++){
				sh.setLM(l,m);
				
				System.out.print(sh.eval(lamda,fai)+"\t");
				
				sh.setLM(l,-m);
				
				System.out.println(sh.eval(lamda,fai));
			}
			
		}catch(Exception e){ e.printStackTrace();}
	}*/
}
