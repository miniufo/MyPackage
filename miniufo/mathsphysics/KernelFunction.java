/**
 * @(#)KernelFunction.java	1.0 2018.08.03
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;


/**
 * Kernel function.
 *
 * @version 1.0, 2018.08.03
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class KernelFunction{
	
	/**
	 * prevent from instantiate
     */
	private KernelFunction(){}
	
	
	/**
	 * Generates the kernel bump along a coordinate definition.
     *
     * @param	xdef	coordinate definition
     * @param	xpos	x-position in the xdef
     * @param	widVar	variance (width) of the Gaussian bump
	 */
	public static float[] gaussian(float[] xdef,float xpos,float widVar){
		int L=xdef.length;
		
		if(L<2) throw new IllegalArgumentException("length of xdef should be at least 2");
		
		float[] ker=new float[L];
		
		for(int l=0;l<L;l++){
			double tmp=xpos-xdef[l]; tmp*=tmp;
			ker[l]=(float)(Math.exp(-tmp/(2.0*widVar))/Math.sqrt(Math.PI*2.0*widVar));
		}
		
		return ker;
	}
	
	
	/** test
	public static void main(String[] args){
		float slope=2;
		float[] bump=gaussian(ArrayUtil.newMonotonousArray(100,slope),40.5f,25);
		
		System.out.println(StatisticsUtil.sum(bump)*slope);
		
		for(int i=0;i<bump.length;i++) System.out.println(bump[i]);
	}*/
}
