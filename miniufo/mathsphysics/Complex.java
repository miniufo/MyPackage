/**
 * @(#)Complex.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import miniufo.basic.Operatable;


/**
 * variable of weather analysis
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Complex implements Cloneable,Operatable<Complex>{
	//
	private float re;
	private float im;
	
	
	/**
	 * Constructor
	 *
	 * @param	x	the real value of a complex number.
	 * @param	y	the imaginary value of a complex number.
	 */
	public Complex(float x,float y){ re=x;	im=y;}

	/**
	 * Constructor
	 *
	 * @param	c	another Complex
	 */
	public Complex(Complex c){ re=c.re;	im=c.im;}
	
	public Complex(){ re=im=0;}
	
	/**
	 * Creates a complex number with the given modulus and argument.
	 *
	 * @param	mod 	the modulus  of a complex number.
	 * @param	arg 	the argument of a complex number.
	 */
	public static Complex polar(float mod,float arg){
		return new Complex(mod*(float)Math.cos(arg),mod*(float)Math.sin(arg));
	}
	
	
	/*** getor and setor ***/
	public float getReal(){ return re;}
	
	public float getImag(){ return im;}
	
	public float  getMod(){ return cMod(re,im);}
	
	public float  getArg(){ return cArg(re,im);}
	
	public float getSquaredMod(){ return re*re+im*im;}
	
	public void setReal(float re){ this.re=re;}
	
	public void setImag(float im){ this.im=im;}
	
	public void setRI(float re,float im){ this.re=re;	this.im=im;}
	
	public void setMA(float mod,float arg){
		re=mod*(float)Math.cos(arg);
		im=mod*(float)Math.sin(arg);
	}
	
	public void set(Complex c){ this.re=c.re;	this.im=c.im;}
	
	
	/**
	 * Returns true if the modulus of this complex number is within the zero tolerance.
	 */
	public boolean isZero(){ return cMod(re,im)==0;}
	
	/**
	 * Returns true if either the real or imaginary part is NaN.
	 */
	public boolean isNaN(){ return Float.isNaN(re)||Float.isNaN(im);}
	
	/**
	 * Returns true if either the real or imaginary part is infinite.
	 */
	public boolean isInfinite(){
		return re==Float.POSITIVE_INFINITY||re==Float.NEGATIVE_INFINITY
		     ||im==Float.POSITIVE_INFINITY||im==Float.NEGATIVE_INFINITY;
	}
	
	
	/**
	 * Returns the inverse of this complex number.
	 */
	public Complex reciprocal(){
		float denominator,real,imag;
		
		if(Math.abs(re)<Math.abs(im)){
			real=re/im;
			imag=-1;
			denominator=re*real+im;
			
		}else{
			real=1;
			imag=-im/re;
			denominator=re-im*imag;
		}
		
		return new Complex(real/denominator,imag/denominator);
	}
	
	public Complex reciprocalEq(){
		float denominator,real,imag;
		
		if(Math.abs(re)<Math.abs(im)){
			real=re/im;
			imag=-1;
			denominator=re*real+im;
			
		}else{
			real=1;
			imag=-im/re;
			denominator=re-im*imag;
		}
		
		re=real/denominator;	im=imag/denominator;
		
		return this;
	}
	
	/**
	 * Returns the complex conjugate of this complex number.
	 */
	public Complex conjugate(){ return new Complex(re,-im);}
	
	public Complex conjugateEq(){ im=-im;	return this;}
	
	/**
	 * Returns the addition of this complex number and another.
	 *
	 * @param z a complex number.
	 */
	public Complex plus(Complex z){ return new Complex(re+z.re,im+z.im);}
	
	public Complex plus(float real){ return new Complex(real+re,im);}
	
	public Complex plusEq(Complex z){
		re+=z.re;			im+=z.im;
		
		return this;
	}
	
	public Complex plusEq(float real){
		re+=real;
		
		return this;
	}
	
	/**
	 * Returns the subtraction of this complex number by another.
	 *
	 * @param z a complex number.
	 */
	public Complex minus(Complex z){ return new Complex(re-z.re,im-z.im);}
	
	public Complex minus(float real){ return new Complex(re-real,im);}
	
	public Complex minusEq(Complex z){
		re-=z.re; 			im-=z.im;
		
		return this;
	}
	
	public Complex minusEq(float real){
		re-=real;
		
		return this;
	}
	
	/**
	 * Returns the multiplication of this complex number and another.
	 *
	 * @param z a complex number.
	 */
	public Complex multiply(Complex z){ return new Complex(re*z.re-im*z.im,re*z.im+im*z.re);}
	
	public Complex multiply(float x){ return new Complex(re*x,im*x);}
	
	public Complex multiplyEq(Complex z){
		float real=re;
		
		re=real*z.re-im*z.im; im=real*z.im+im*z.re;
		
		return this;
	}
	
	public Complex multiplyEq(float x) {
		re*=x;	im*=x;
		
		return this;
	}
	
	/**
	 * Returns the division of this complex number by another.
	 *
	 * @param z a complex number.
	 */
	public Complex divide(Complex z){
		if(Math.abs(z.re)<Math.abs(z.im)){
			float a=z.re/z.im;
			float denominator=z.re*a+z.im;
			
			return new Complex((re*a+im)/denominator,(im*a-re)/denominator);
			
		}else{
			float a=z.im/z.re;
			float denominator=z.re+z.im*a;
			
			return new Complex((im*a+re)/denominator,(im-re*a)/denominator);
		}
	}
	
	public Complex divide(float x){ return new Complex(re/x,im/x);}
	
	public Complex divideEq(Complex z){
		float denominator,a;
		
		if(Math.abs(z.re)<Math.abs(z.im)){
			a=z.re/z.im;
			denominator=z.re*a+z.im;
			re=re*a+im;
			im=im*a-re;
			
		}else{
			a=z.im/z.re;
			denominator=z.re+z.im*a;
			re=re+im*a;
			im=im-re*a;
		}
		
		re/=denominator;	im/=denominator;
		
		return this;
	}
	
	public Complex divideEq(float x){
		re/=x;	im/=x;
		
		return this;
	}
	
	/**
	 * Returns this complex number raised to the power of another.
	 *
	 * @param z a complex number.
	 */
	public Complex pow(Complex z){
		float Mod=cMod(re,im);	float Arg=cArg(re,im);
		
		return polar((float)(Math.pow(Mod,z.re)/Math.exp(Arg*z.im)),Arg*z.re+(float)Math.log(Mod)*z.im);
	}
	
	public Complex pow(float x){ return polar((float)Math.pow(cMod(re,im),x),cArg(re,im)*x);}
	
	public Complex powEq(Complex z){
		float thisMod=cMod(re,im);
		
		float mod=(float)(Math.pow(thisMod,z.re)/Math.exp(cArg(re,im)*z.im));
		float arg=cArg(re,im)*z.re+(float)Math.log(thisMod)*z.im;
		
		re=mod*(float)Math.cos(arg);	im=mod*(float)Math.sin(arg);
		
		return this;
	}
	
	public Complex powEq(float x){
		float mod=(float)Math.pow(cMod(re,im),x);
		float arg=cArg(re,im)*x;
		
		re=mod*(float)Math.cos(arg);
		im=mod*(float)Math.sin(arg);
		
		return this;
	}
	
	/**
	 * Returns the square of this complex number.
	 */
	public Complex square(){ return new Complex(re*re-im*im,2*re*im);}
	
	public Complex squareEq(){
		float real=re;
		
		re=re*re-im*im; im*=2*real;
		
		return this;
	}
	
	/**
	 * Returns the square root of this complex number.
	 */
	public Complex sqrt(){
		float mod=(float)Math.sqrt(cMod(re,im));
		float arg=cArg(re,im)/2;
		
		return polar(mod,arg);
	}
	
	public Complex sqrtEq(){
		float mod=(float)Math.sqrt(cMod(re,im));
		float arg=cArg(re,im)/2;
		
		re=mod*(float)Math.cos(arg);	im=mod*(float)Math.sin(arg);
		
		return this;
	}
	
	/**
	 * Returns the exponent of this complex number.
	 */
	public Complex exp(){
		Complex result=new Complex(this);
		
		result.setRI(
			(float)(Math.exp(re)*Math.cos(im)),
			(float)(Math.exp(re)*Math.sin(im))
		);
		
		return result;
	}
	
	public Complex expEq(){
		float retmp=re;
		
		re=(float)(Math.exp(retmp)*Math.cos(im));
		im=(float)(Math.exp(retmp)*Math.sin(im));
		
		return this;
	}
	
	/**
	 * Returns the abstract value (mod) of this complex number.
	 */
	public Complex abs(){ return new Complex(cMod(re,im),0);}
	
	public Complex absEq(){
		setRI(cMod(re,im),0);
		
		return this;
	}
	
	
	/**
     * used to test whether this complex number equals another
     */
	public boolean equals(Object o){
		if(o instanceof Complex){
			Complex c=(Complex)o;
			
			float mod=cMod(re,im);
			float arg=cArg(re,im);
			
			if(mod==0) return mod==c.getMod();
			else return mod==c.getMod()&&arg==c.getArg();
			
		}else return false;
	}
	
	/**
     * used to print out
     */
	public String toString(){ return re+(im>=0?"+":"")+im+"i";}
	
	
	/** Make a deep copy of a Complex*/
	public Object clone(){
		try{
			Complex c=null;	c=(Complex)super.clone();
			
			c.im=im;	c.re=re;
			
			return c;
			
		}catch(CloneNotSupportedException ex){
			// this shouldn't happen, since we are Cloneable
			throw new InternalError();
		}
	}
	
	
	/**
	 * calculate modulus
	 *
	 * @param	rep		real  part
	 * @param	imp		image part
	 */
	private float cMod(float rep,float imp){
		return (float)Math.hypot(rep,imp);
	}
	
	/**
	 * calculate argument
	 *
	 * @param	rep		real  part
	 * @param	imp		image part
	 */
	private float cArg(float rep,float imp){ return (float)Math.atan2(imp,rep);}
	
	
	/** test
	public static void main(String[] args){
		try{
			Complex c=new Complex(2,9);
			float r=10;
			float i=-0.01f;
			System.out.println(c.cMod(r,i));
			System.out.println(c.cMod2(r,i));
	    	
	    }catch(Exception ex){ ex.printStackTrace(); System.exit(0);}
	}*/
}
