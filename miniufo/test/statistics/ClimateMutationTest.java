/**
 * @(#)ClimateMutationTest.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.statistics;

import miniufo.statistics.StatisticsUtil;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;


/**
 * used to test climate mutation
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ClimateMutationTest extends StatisticsUtil{
	
	/**
     * slipped T test
     *
     * @param	arr_data		data in an array
     * @param	slipped_length	slip length
     *
     * @return	t				result
     */
	public static float[] cSlippedTTest(float[] data,int slipLen){
		int length=data.length;
		if(slipLen<2||length<2*slipLen) throw new IllegalArgumentException();
		
		float ave1,ave2,var1,var2;
		float[] t  =null;	t  =new float[length-2*slipLen+1];
		float[] tmp=null;	tmp=new float[slipLen];
		
		for(int i=0;i<=length-2*slipLen;i++){
			for(int j=0;j<slipLen;j++) tmp[j]=data[i+j];
			
			ave1=(float)cArithmeticMean(tmp);
			var1=(float)cVariance(tmp);
			
			for(int j=0;j<slipLen;j++) tmp[j]=data[i+j+slipLen];
			
			ave2=(float)cArithmeticMean(tmp);
			var2=(float)cVariance(tmp);
			
			t[i]=(ave1-ave2)/(float)sqrt(slipLen*(var1+var2)/(2*slipLen-2))/(float)sqrt(2f/slipLen);
		}
		
		tmp=null;
		
		return t;
	}
	
	public static float[] cslipLenpedTTest(double[] data,int slipLen){
		int length=data.length;
		if(slipLen<2||length<2*slipLen) throw new IllegalArgumentException();
		
		float ave1,ave2,var1,var2;
		float[] t  =null;	t  =new float[length-2*slipLen+1];
		float[] tmp=null;	tmp=new float[slipLen];
		
		for(int i=0;i<=length-2*slipLen;i++){
			for(int j=0;j<slipLen;j++) tmp[j]=(float)data[i+j];
			
			ave1=(float)cArithmeticMean(tmp);
			var1=(float)cVariance(tmp);
			
			for(int j=0;j<slipLen;j++) tmp[j]=(float)data[i+j+slipLen];
			
			ave2=(float)cArithmeticMean(tmp);
			var2=(float)cVariance(tmp);
			
			t[i]=(ave1-ave2)/(float)sqrt(slipLen*(var1+var2)/(2*slipLen-2))*(float)sqrt(2/slipLen);
		}
		
		tmp=null;
		
		return t;
	}
	
	
	/**
     * Cramer test
     *
     * @param	data	data in an array
     * @param	sub_length	length
     *
     * @return	t			result
     */
	public static float[] cCramerTest(float[] data,int subLen){
		int length=data.length;
		if(subLen<2||length<=subLen) throw new IllegalArgumentException();
		
		float ave,ave1,var,temp;
		float[] t  =null;	t  =new float[length-subLen+1];
		float[] tmp=null;	tmp=new float[subLen];
		
		ave=(float)cArithmeticMean(data);
		var=(float)cVariance(data);
		
		for(int i=0;i<=length-subLen;i++){
			for(int j=0;j<subLen;j++) tmp[j]=data[i+j];
			
			ave1=(float)cArithmeticMean(tmp);
			temp=(ave1-ave)/var;
			
			t[i]=temp*(float)sqrt((subLen*(length-2))/(length-subLen*(1+temp*temp)));
		}
		
		tmp=null;
		
		return t;
	}
	
	public static float[] cCramerTest(double[] data,int subLen){
		int length=data.length;
		if(subLen<2||length<=subLen) throw new IllegalArgumentException();
		
		float ave,ave1,var,temp;
		float[] t  =null;	t  =new float[length-subLen+1];
		float[] tmp=null;	tmp=new float[subLen];
		
		ave=(float)cArithmeticMean(data);
		var=(float)cVariance(data);
		
		for(int i=0;i<=length-subLen;i++){
			for(int j=0;j<subLen;j++) tmp[j]=(float)data[i+j];
			
			ave1=(float)cArithmeticMean(tmp);
			temp=(ave1-ave)/var;
			
			t[i]=temp*(float)sqrt((subLen*(length-2))/(length-subLen*(1+temp*temp)));
		}
		
		tmp=null;
		
		return t;
	}
	
	
	/**
     * Yamamoto test
     *
     * @param	data		data in an array
     * @param	slipLen	slipLen length
     *
     * @return	t				result
     */
	public static float[] cYamamotoTest(float[] data,int slipLen){
		int length=data.length;
		if(slipLen<2||length<2*slipLen) throw new IllegalArgumentException();
		
		float ave1,ave2,var1,var2;
		float[] t  =null;	t  =new float[length-2*slipLen+1];
		float[] tmp=null;	tmp=new float[slipLen];
		
		for(int i=0;i<=length-2*slipLen;i++){
			for(int j=0;j<slipLen;j++) tmp[j]=data[i+j];
			
			ave1=(float)cArithmeticMean(tmp);
			var1=(float)cVariance(tmp);
			
			for(int j=0;j<slipLen;j++) tmp[j]=data[i+j+slipLen];
			
			ave2=(float)cArithmeticMean(tmp);
			var2=(float)cVariance(tmp);
			
			t[i]=abs(ave1-ave2)/(var1+var2);
		}
		
		tmp=null;
		
		return t;
	}
	
	public static float[] cYamamotoTest(double[] data,int slipLen){
		int length=data.length;
		if(slipLen<2||length<2*slipLen) throw new IllegalArgumentException();
		
		float ave1,ave2,var1,var2;
		float[] t  =null;	t  =new float[length-2*slipLen+1];
		float[] tmp=null;	tmp=new float[slipLen];
		
		for(int i=0;i<=length-2*slipLen;i++){
			for(int j=0;j<slipLen;j++) tmp[j]=(float)data[i+j];
			
			ave1=(float)cArithmeticMean(tmp);
			var1=(float)cVariance(tmp);
			
			for(int j=0;j<slipLen;j++) tmp[j]=(float)data[i+j+slipLen];
			
			ave2=(float)cArithmeticMean(tmp);
			var2=(float)cVariance(tmp);
			
			t[i]=abs(ave1-ave2)/(var1+var2);
		}
		
		tmp=null;
		
		return t;
	}
	
	
	/**
     * Mann-Kendall test
     *
     * @param	data	data in an array
     *
     * @return	t			result
     */
	public static float[] cMannKendallTest(float[] data){
		int length=data.length;
		if(length<3) throw new IllegalArgumentException();
		
		float r=0,s=0,Es=0,Vs=0;
		float[] t=null;	t=new float[length];
		
		for(int i=0;i<length;i++){
			for(int j=0;j<i;j++) if(data[i]>data[j]) r++;
			
			for(int j=0;j<=i;j++) s+=r;
			
			Es=(float)(i+1)*(i+2)/4;
			Vs=(float)i*(i+1)*(2*(i+1)+5)/72;
			
			t[i]=(s-Es)/(float)sqrt(Vs);
			
			r=0;	s=0;
		}
		
		t[0]=0;
		
		return t;
	}
	
	
	
	/** test
	public static void main(String[] args){
		try{
			Scanner sn=new Scanner(new File("D:\\CourseWare\\Statistics\\homework\\t654.txt"));
			float[] data=null; data=new float[50];
			
			for(int i=0;i<data.length;i++) data[i]=Float.parseFloat(sn.next());
			
			sn.close();
			
			float[] UF=cMannKendallTest(data);
			
			BufferedWriter bw=new BufferedWriter(new FileWriter("D:\\CourseWare\\Statistics\\homework\\t.txt"));
			//for(int i=0;i<length;i++) bw.write(-UF[length-1-i]+" ");
			for(int i=0;i<UF.length;i++) bw.write(UF[i]+" ");
			bw.close();
	    	
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}