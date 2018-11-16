/**
 * @(#)Integrator1D.java	1.0 2018.09.14
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;


/**
 * One-dimensional integrator using Trapezoid method.
 *
 * @version 1.0, 2018.09.14
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Integrator1D{
	
	/**
	 * Prevent from constructing
	 */
	private Integrator1D(){}
	
	
	/**
	 * Integrate the discrete integrand from init to fnl, using the given step.
	 * 
	 * @param	step		dx
	 * @param	integrand	discrete integrand (derivatives)
	 * @param	init		initial value at first point re[0]
	 * @param	fnl			final value at last point re[end], used for correction
	 */
	public static float[] integrateForward(float step,float[] integrand,float init,float fnl){
		float[] re=integrateForward(step,integrand,init);
		
		correctForward(re,fnl);
		
		return re;
	}
	
	public static float[] integrateForward(float step,float[] integrand,float init){
		int len=integrand.length;
		
		float[] re=new float[len]; re[0]=init;
		
		for(int i=1;i<len;i++) re[i]=re[i-1]+(integrand[i-1]+integrand[i])/2f*step;
		
		return re;
	}
	
	public static double[] integrateForward(double step,double[] integrand,double init,double fnl){
		double[] re=integrateForward(step,integrand,init);
		
		correctForward(re,fnl);
		
		return re;
	}
	
	public static double[] integrateForward(double step,double[] integrand,double init){
		int len=integrand.length;
		
		double[] re=new double[len]; re[0]=init;
		
		for(int i=1;i<len;i++) re[i]=re[i-1]+(integrand[i-1]+integrand[i])/2.0*step;
		
		return re;
	}
	
	
	/**
	 * Integrate the discrete integrand from init to fnl, using the given step.
	 * 
	 * @param	step		dx
	 * @param	integrand	discrete integrand (derivatives)
	 * @param	init		initial value at first point re[0]
	 * @param	fnl			final value at last point re[end], used for correction
	 */
	public static float[] integrateBackward(float step,float[] integrand,float init,float fnl){
		float[] re=integrateBackward(step,integrand,init);
		
		correctBackward(re,fnl);
		
		return re;
	}
	
	public static float[] integrateBackward(float step,float[] integrand,float init){
		int len=integrand.length;
		
		float[] re=new float[len]; re[len-1]=init;
		
		for(int i=len-2;i>=0;i++) re[i]=re[i+1]-(integrand[i]+integrand[i+1])/2f*step;
		
		return re;
	}
	
	public static double[] integrateBackward(double step,double[] integrand,double init,double fnl){
		double[] re=integrateBackward(step,integrand,init);
		
		correctBackward(re,fnl);
		
		return re;
	}
	
	public static double[] integrateBackward(double step,double[] integrand,double init){
		int len=integrand.length;
		
		double[] re=new double[len]; re[len-1]=init;
		
		for(int i=len-2;i>=0;i--) re[i]=re[i+1]-(integrand[i]+integrand[i+1])/2.0*step;
		
		return re;
	}
	
	
	
	/*** helper methods ***/
	
	/**
	 * Correct the integration given the value at last point,
	 * so that after correction we have re[end]=endVal.
	 * 
	 * @param	re		result of integration
	 * @param	endVal	value at the last point
	 */
	private static void correctForward(float[] re,float endVal){
		int len=re.length;
		
		float corr=(re[len-1]-endVal)/(len-1);
		
		for(int i=1;i<len;i++) re[i]-=corr*(float)i;
	}
	
	private static void correctForward(double[] re,double endVal){
		int len=re.length;
		
		double corr=(re[len-1]-endVal)/(len-1);
		
		for(int i=1;i<len;i++) re[i]-=corr*(double)i;
	}
	
	/**
	 * Correct the integration given the value at first point,
	 * so that after correction we have re[0]=endVal.
	 * 
	 * @param	re		result of integration
	 * @param	endVal	value at the first point
	 */
	private static void correctBackward(float[] re,float endVal){
		int len=re.length;
		
		float corr=(re[0]-endVal)/(len-1);
		
		for(int i=len-2;i>=0;i--) re[i]-=corr*(float)(len-i-1);
	}
	
	private static void correctBackward(double[] re,double endVal){
		int len=re.length;
		
		double corr=(re[0]-endVal)/(len-1);
		
		for(int i=len-2;i>=0;i--) re[i]-=corr*(double)(len-i-1);
	}
	
	
	/** test
	public static void main(String[] args){
		int nmode=4;
		int count=201;
		float[] data=TextReader.readColumn("d:/Matlab/EEMD/y.txt",count,false,1)[0];
		
		float[][] mode=new EEMD(count,nmode,0.02f).decomp(data);
		
		for(int i=0;i<count;i++){
			System.out.println(
				mode[0][i]+"\t"+
				mode[1][i]+"\t"+
				mode[2][i]+"\t"+
				mode[3][i]+"\t"+
				mode[4][i]+"\t"+
				mode[5][i]
			);
		}
	}*/
}
