/**
 * @(#)ArrayUtil.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.basic;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Random;


/**
 * a class associated with the array operation
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class ArrayUtil{
	
	/**
	 * prevent from instantiate
	 */
	private ArrayUtil(){}
	
	
	/**
     * return a concatenated array
     *
     * @param	first	first array
     * @param	rests	rest arrays or elements
     */
	@SuppressWarnings("unchecked")
	public static <T> T[] concatAll(Class<T> ty,T[]... rests){
		int len=0,off=0;
		
		for(T[] array:rests) len+=array.length;
		
		T[] re=(T[])Array.newInstance(ty,len);
		
		for(T[] array:rests){
			System.arraycopy(array,0,re,off,array.length);
			off+=array.length;
		}
		
		return re;
	}
	
	@SuppressWarnings("unchecked")
	public static <T> T[] concatAll(Class<T> ty,T[] first,T sec,T... rests){
		int len=rests.length+1+first.length;
		
		T[] re=(T[])Array.newInstance(ty,len);
		
		System.arraycopy(first,0,re,0,first.length);
		
		re[first.length]=sec;
		
		System.arraycopy(rests,0,re,first.length+1,rests.length);
		
		return re;
	}
	
	
	/**
     * return a monotonous array
     *
     * @param	len		array length
     * @param	slope	slope
     */
	public static float[] newMonotonousArray(int len,float slope){
		float[] re=new float[len];
		
		for(int l=0;l<len;l++) re[l]=slope*l;
		
		return re;
	}
	
	
	/**
     * return a harmonic series with specified parameters in the form:
     * y(t)=A*sin(2*pi/T*t+fai)
     *
     * @param	len		length of data
     * @param	A		amplitude
     * @param	T		period
     * @param	fai		initial phase
     */
	public static float[] newHarmonicSeriesWithAmp(int len,float A,float T,float fai){
		float[] re=new float[len];
		
		double omega=2.0*Math.PI/T;
		
		for(int l=0;l<len;l++) re[l]=A*(float)Math.sin(omega*l+fai);
		
		return re;
	}
	
	/**
     * return a harmonic series with specified parameters in the form:
     * y(t)=a*sin(2*pi/T*t)+b*cos(2*pi/T*t)
     *
     * @param	len		length of data
     * @param	a		coefficient of sin
     * @param	b		coefficient of cos
     * @param	T		period
     */
	public static float[] newHarmonicSeriesWithCoeff(int len,float a,float b,float T){
		float[] re=new float[len];
		
		double omega=2.0*Math.PI/T;
		
		for(int l=0;l<len;l++){
			double fai=omega*l;
			
			re[l]=(float)(a*Math.sin(fai)+b*Math.cos(fai));
		}
		
		return re;
	}
	
	
	/**
     * return a Gaussian-distributed random time series with mean and amplitude specified
     *
     * @param	len		length of data
     * @param	mean	mean of the data
     * @param	amp		amplitude of the data
     */
	public static float[] newGaussianRandomSeries(int len,float mean,float amp){
		Random r=new Random();
		
		float[] re=new float[len];
		
		for(int l=0;l<len;l++) re[l]=mean+amp*(float)r.nextGaussian();
		
		return re;
	}
	
	/**
     * return a Uniform-distributed random time series with mean and amplitude specified
     *
     * @param	len		length of data
     * @param	str		the start of the range of distribution e.g., 0
     * @param	range	range of data, for (R1, R2), range = R2 - R1
     */
	public static float[] newUniformRandomSeries(int len,float str,float range){
		Random r=new Random();
		
		float[] re=new float[len];
		
		for(int l=0;l<len;l++) re[l]=str+range*r.nextFloat();
		
		return re;
	}
	
	/**
     * return a 1st-order autoregressive time series AR(1)
     *
     * @param	len		length of data
     * @param	dt		sampling interval
     * @param	T		decaying timescale
     * @param	K		diffusivity = variance * T
     */
	public static float[] newAR1RandomSeries(int len,float dt,float T,float K){
		Random r=new Random();
		
		float[] re=new float[len];
		
		float spinup=0;
		
		for(int l=0;l<50;l++) spinup=(1f-dt/T)*spinup+(float)(Math.sqrt(2.0*dt*K)/T*r.nextGaussian());
		
		re[0]=spinup;
		
		for(int l=1;l<len;l++) re[l]=(1f-dt/T)*re[l-1]+(float)(Math.sqrt(2.0*dt*K)/T*r.nextGaussian());
		
		return re;
	}
	
	
	/**
	 * mirror-symmetric extension of a series of data to a given length
	 */
	public static float[] mirrorSymmetricExtention(float[] data,int length){
		float[] re=new float[length];
		
		mirrorSymmetricExtention(data,re);
		
		return re;
	}
	
	public static void mirrorSymmetricExtention(float[] data,float[] result){
		int dl=data.length,dl2=dl*2-2;
		int rl=result.length;
		
		if(rl<=dl)
			throw new IllegalArgumentException("result's length should be larger than data's");
		
		for(int i=0;i<dl;i++) result[i]=data[i];
		
		for(int i=dl;i<rl;i++) result[i]=data[dl2-i];
	}
	
	
	/**
	 * periodic extension of a series of data to a given length
	 */
	public static float[] periodicExtention(float[] data,int length){
		float[] re=new float[length];
		
		periodicExtention(data,re);
		
		return re;
	}
	
	public static void periodicExtention(float[] data,float[] result){
		int dl=data.length;
		int rl=result.length;
		
		if(rl<=dl)
			throw new IllegalArgumentException("result's length should be larger than data's");
		
		for(int i=0;i<rl;i++) result[i]=data[i%dl];
	}
	
	
	/**
	 * Build up a string from a given array in the following form:
	 * prefix data[0] delimiter data[1] delimiter ... delimiter data[end] suffix
	 * 
	 * @param	data	a given array
	 * @param	pre		prefix
	 * @param	del		delimiter
	 * @param	suf		suffix
	 */
	public static String allToString(float[] data,String pre,String del,String suf){
		StringBuilder sb=new StringBuilder();
		
		sb.append(pre);
		
		for(int i=0,I=data.length-1;i<I;i++) sb.append(data[i]+del);
		
		sb.append(data[data.length-1]+suf);
		
		return sb.toString();
	}
	
	public static String allToString(double[] data,String pre,String del,String suf){
		StringBuilder sb=new StringBuilder();
		
		sb.append(pre);
		
		for(int i=0,I=data.length-1;i<I;i++) sb.append(data[i]+del);
		
		sb.append(data[data.length-1]+suf);
		
		return sb.toString();
	}
	
	
	/**
	 * split a String into String array with fixed length
	 */
	public static String[] splitByLength(String s,int splitLen){
		if(splitLen<1) throw new IllegalArgumentException("split length "+splitLen+" should be positive");
		if(s==null||"".equals(s)) throw new IllegalArgumentException("empty string");
		
		int arrayNum=0;
		
		if(s.length()%splitLen==0) arrayNum=s.length()/splitLen;
		else arrayNum=s.length()/splitLen+1;
		
		String[] strArray=new String[arrayNum];
		
		for(int i=0,start=0,I=arrayNum-1;i<I;i++,start+=splitLen)
		strArray[i]=s.substring(start,start+splitLen);
		
		strArray[arrayNum-1]=s.substring(splitLen*(arrayNum-1));
		
		return strArray;
	}
	
	public static String[] splitByCount(String s,int splitCount){
		if(splitCount<1)
			throw new IllegalArgumentException("split count "+splitCount+" should be positive");
		if(s==null||"".equals(s))
			throw new IllegalArgumentException("empty string");
		if(s.length()%splitCount!=0)
			throw new IllegalArgumentException("string length cannot be divided by count");
		
		return splitByLength(s,s.length()/splitCount);
	}
	
	public static void splitByLength(String s,int splitLen,String[] strArray){
		if(splitLen<1) throw new IllegalArgumentException("split length "+splitLen+" should be positive");
		if(s==null||"".equals(s)) throw new IllegalArgumentException("empty string");
		
		int arrayNum=0;
		
		if(s.length()%splitLen==0) arrayNum=s.length()/splitLen;
		else arrayNum=s.length()/splitLen+1;
		
		if(arrayNum!=strArray.length) throw new IllegalArgumentException("lengths not equal");
		
		for(int i=0,start=0,I=arrayNum-1;i<I;i++,start+=splitLen)
		strArray[i]=s.substring(start,start+splitLen);
		
		strArray[arrayNum-1]=s.substring(splitLen*(arrayNum-1));
	}
	
	public static void splitByCount(String s,int splitCount,String[] strArray){
		if(splitCount<1)
			throw new IllegalArgumentException("split count "+splitCount+" should be positive");
		if(s==null||"".equals(s))
			throw new IllegalArgumentException("empty string");
		if(s.length()%splitCount!=0)
			throw new IllegalArgumentException("string length cannot be divided by count");
		
		splitByLength(s,s.length()/splitCount,strArray);
	}
	
	
	/**
	 * Get the index of the largest element which is smaller than a given value in an increasing array
	 * 
	 * Return -1 if value < data[0] or value > data[data.length-1]
	 * 
	 * @param	data	a given increasing array
	 * @param	value	a specified value
	 */
	public static int getLEIdxIncre(int[] data,int value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		if(idx<0) idx=-idx-2;
		
		return idx;
	}
	
	public static int getLEIdxIncre(short[] data,short value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		if(idx<0) idx=-idx-2;
		
		return idx;
	}
	
	public static int getLEIdxIncre(long[] data,long value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		if(idx<0) idx=-idx-2;
		
		return idx;
	}
	
	public static int getLEIdxIncre(float[] data,float value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		if(idx<0) idx=-idx-2;
		
		return idx;
	}
	
	public static int getLEIdxIncre(double[] data,double value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		if(idx<0) idx=-idx-2;
		
		return idx;
	}
	
	
	/**
	 * Get the index of the largest element which is smaller than a given value in an decreasing array
	 * 
	 * Return -1 if value > data[0] or value < data[data.length-1]
	 * 
	 * @param	data	a given decreasing array
	 * @param	value	a specified value
	 */
	public static int getLEIdxDecre(int[] data,int value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		if(idx<0) idx=-idx-1;
		
		return idx;
	}
	
	public static int getLEIdxDecre(short[] data,short value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		if(idx<0) idx=-idx-1;
		
		return idx;
	}
	
	public static int getLEIdxDecre(long[] data,long value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		if(idx<0) idx=-idx-1;
		
		return idx;
	}
	
	public static int getLEIdxDecre(float[] data,float value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		if(idx<0) idx=-idx-1;
		
		return idx;
	}
	
	public static int getLEIdxDecre(double[] data,double value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		if(idx<0) idx=-idx-1;
		
		return idx;
	}
	
	
	/**
	 * Get the index of the element which is closest to a given value in an increasing array
	 * 
	 * Return -1 if value < data[0] or value > data[data.length-1]
	 * 
	 * @param	data	a given increasing array
	 * @param	value	a specified value
	 */
	public static int getIdxIncre(int[] data,int value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		
		if(idx<0){
			idx=-idx-2;
			
			if(value-data[idx]<=data[idx+1]-value) return idx;
			else return idx+1;
		}
		
		return idx;
	}
	
	public static int getIdxIncre(short[] data,short value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		
		if(idx<0){
			idx=-idx-2;
			
			if(value-data[idx]<=data[idx+1]-value) return idx;
			else return idx+1;
		}
		
		return idx;
	}
	
	public static int getIdxIncre(long[] data,long value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		
		if(idx<0){
			idx=-idx-2;
			
			if(value-data[idx]<=data[idx+1]-value) return idx;
			else return idx+1;
		}
		
		return idx;
	}
	
	public static int getIdxIncre(float[] data,float value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		
		if(idx<0){
			idx=-idx-2;
			
			if(value-data[idx]<=data[idx+1]-value) return idx;
			else return idx+1;
		}
		
		return idx;
	}
	
	public static int getIdxIncre(double[] data,double value){
		if(value<data[0]||value>data[data.length-1]) return -1;
		
		int idx=Arrays.binarySearch(data,value);
		
		if(idx<0){
			idx=-idx-2;
			
			if(value-data[idx]<=data[idx+1]-value) return idx;
			else return idx+1;
		}
		
		return idx;
	}
	
	
	/**
	 * Get the index of the element which is closest to a given value in an decreasing array
	 * 
	 * Return -1 if value > data[0] or value < data[data.length-1]
	 * 
	 * @param	data	a given increasing array
	 * @param	value	a specified value
	 */
	public static int getIdxDecre(int[] data,int value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		
		if(idx<0){
			idx=-idx;
			
			if(value-data[idx-1]>=data[idx]-value) return idx-1;
			else return idx;
		}
		
		return idx;
	}
	
	public static int getIdxDecre(short[] data,short value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		
		if(idx<0){
			idx=-idx;
			
			if(value-data[idx-1]>=data[idx]-value) return idx-1;
			else return idx;
		}
		
		return idx;
	}
	
	public static int getIdxDecre(long[] data,long value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		
		if(idx<0){
			idx=-idx;
			
			if(value-data[idx-1]>=data[idx]-value) return idx-1;
			else return idx;
		}
		
		return idx;
	}
	
	public static int getIdxDecre(float[] data,float value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		
		if(idx<0){
			idx=-idx;
			
			if(value-data[idx-1]>=data[idx]-value) return idx-1;
			else return idx;
		}
		
		return idx;
	}
	
	public static int getIdxDecre(double[] data,double value){
		if(value>data[0]||value<data[data.length-1]) return -1;
		
		int idx=binarySearchDecre(data,value);
		
		if(idx<0){
			idx=-idx;
			
			if(value-data[idx-1]>=data[idx]-value) return idx-1;
			else return idx;
		}
		
		return idx;
	}
	
	
	/**
	 * deep comparison of two arrays
	 */
	public static boolean arrayEquals(int[] a1,int[] a2){
		if(a1==a2) return true;
		
		int length=a1.length;
		
		if(length!=a2.length) return false;
		
		for(int i=0;i<length;i++)
		if(a1[i]!=a2[i]) return false;
		
		return true;
	}
	
	public static boolean arrayEquals(long[] a1,long[] a2){
		if(a1==a2) return true;
		
		int length=a1.length;
		
		if(length!=a2.length) return false;
		
		for(int i=0;i<length;i++)
		if(a1[i]!=a2[i]) return false;
		
		return true;
	}
	
	public static boolean arrayEquals(float[] a1,float[] a2){
		if(a1==a2) return true;
		
		int length=a1.length;
		
		if(length!=a2.length) return false;
		
		for(int i=0;i<length;i++)
		if(a1[i]!=a2[i]) return false;
		
		return true;
	}
	
	public static boolean arrayEquals(double[] a1,double[] a2){
		if(a1==a2) return true;
		
		int length=a1.length;
		
		if(length!=a2.length) return false;
		
		for(int i=0;i<length;i++)
		if(a1[i]!=a2[i]) return false;
		
		return true;
	}
	
	
	/**
	 * reverse the sequence of the elements in an array
	 */
	public static void reverse(int[] array){
		for(int i=0,I=array.length,hf=I/2;i<hf;i++){
			int tmp=array[i];
			
			array[i]=array[I-1-i];
			array[I-1-i]=tmp;
		}
	}
	
	public static void reverse(long[] array){
		for(int i=0,I=array.length,hf=I/2;i<hf;i++){
			long tmp=array[i];
			
			array[i]=array[I-1-i];
			array[I-1-i]=tmp;
		}
	}
	
	public static void reverse(float[] array){
		for(int i=0,I=array.length,hf=I/2;i<hf;i++){
			float tmp=array[i];
			
			array[i]=array[I-1-i];
			array[I-1-i]=tmp;
		}
	}
	
	public static void reverse(double[] array){
		for(int i=0,I=array.length,hf=I/2;i<hf;i++){
			double tmp=array[i];
			
			array[i]=array[I-1-i];
			array[I-1-i]=tmp;
		}
	}
	
	public static void reverse(Object[] array){
		for(int i=0,I=array.length,hf=I/2;i<hf;i++){
			Object tmp=array[i];
			
			array[i]=array[I-1-i];
			array[I-1-i]=tmp;
		}
	}
	
	
	/**
     * get maximum of absolute values in a array
     *
     * @param	arr_data	an array
     *
     * @return	the maximum value
     *
     * @exception	if array length is zero
     */
	public static int getAbsMax(int[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=0;
		
		for(int ei:array){
			int b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static int getAbsMax(int[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=0;
		
		for(int[] ej:array)
		for(int ei:ej){
			int b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static int getAbsMax(int[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=0;
		
		for(int[][] ek:array)
		for(int[] ej:ek)
		for(int ei:ej){
			int b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static int getAbsMax(int[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=0;
		
		for(int[][][] el:array)
		for(int[][] ek:el)
		for(int[] ej:ek)
		for(int ei:ej){
			int b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	public static long getAbsMax(long[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=0;
		
		for(long ei:array){
			long b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static long getAbsMax(long[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=0;
		
		for(long[] ej:array)
		for(long ei:ej){
			long b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static long getAbsMax(long[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=0;
		
		for(long[][] ek:array)
		for(long[] ej:ek)
		for(long ei:ej){
			long b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static long getAbsMax(long[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=0;
		
		for(long[][][] el:array)
		for(long[][] ek:el)
		for(long[] ej:ek)
		for(long ei:ej){
			long b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	public static float getAbsMax(float[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=0;
		
		for(float ei:array){
			float b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static float getAbsMax(float[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=0;
		
		for(float[] ej:array)
		for(float ei:ej){
			float b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static float getAbsMax(float[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=0;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej){
			float b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static float getAbsMax(float[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=0;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej){
			float b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	public static double getAbsMax(double[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=0;
		
		for(double ei:array){
			double b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static double getAbsMax(double[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=0;
		
		for(double[] ej:array)
		for(double ei:ej){
			double b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static double getAbsMax(double[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=0;
		
		for(double[][] ek:array)
		for(double[] ej:ek)
		for(double ei:ej){
			double b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static double getAbsMax(double[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=0;
		
		for(double[][][] el:array)
		for(double[][] ek:el)
		for(double[] ej:ek)
		for(double ei:ej){
			double b=(ei<0)?(-ei):ei;
			if(b>tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	/**
     * get maximum of values in a array
     *
     * @param	array	an array
     *
     * @return	the maximum value
     *
     * @exception	if array length is zero
     */
	public static int getMax(int[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MIN_VALUE;
		
		for(int ei:array) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static int getMax(int[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MIN_VALUE;
		
		for(int[] ej:array)
		for(int ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static int getMax(int[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MIN_VALUE;
		
		for(int[][] ek:array)
		for(int[] ej:ek)
		for(int ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static int getMax(int[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MIN_VALUE;
		
		for(int[][][] el:array)
		for(int[][] ek:el)
		for(int[] ej:ek)
		for(int ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static long getMax(long[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MIN_VALUE;
		
		for(long ei:array) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static long getMax(long[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MIN_VALUE;
		
		for(long[] ej:array)
		for(long ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static long getMax(long[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MIN_VALUE;
		
		for(long[][] ek:array)
		for(long[] ej:ek)
		for(long ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static long getMax(long[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MIN_VALUE;
		
		for(long[][][] el:array)
		for(long[][] ek:el)
		for(long[] ej:ek)
		for(long ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static float getMax(float[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float ei:array) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMax(float[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMax(float[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMax(float[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static float getMax(float[] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float ei:array) if(ei!=undef&&ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMax(float[][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej) if(ei!=undef&&ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMax(float[][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef&&ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMax(float[][][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=-Float.MAX_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef&&ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static double getMax(double[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=-Double.MAX_VALUE;
		
		for(double ei:array) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static double getMax(double[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=-Double.MAX_VALUE;
		
		for(double[] ej:array)
		for(double ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static double getMax(double[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=-Double.MAX_VALUE;
		
		for(double[][] ek:array)
		for(double[] ej:ek)
		for(double ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	public static double getMax(double[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=-Double.MAX_VALUE;
		
		for(double[][][] el:array)
		for(double[][] ek:el)
		for(double[] ej:ek)
		for(double ei:ej) if(ei>tmp) tmp=ei;
		
		return tmp;
	}
	
	
	/**
     * get minimum of absolute value in a array
     *
     * @param	arr_data	an array
     *
     * @return	the minimum value
     *
     * @exception	if array length is zero
     */
	public static int getAbsMin(int[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int ei:array){
			int b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
				
		return tmp;
	}

	public static int getAbsMin(int[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int[] ej:array)
		for(int ei:ej){
			int b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static int getAbsMin(int[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int[][] ek:array)
		for(int[] ej:ek)
		for(int ei:ej){
			int b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static int getAbsMin(int[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int[][][] el:array)
		for(int[][] ek:el)
		for(int[] ej:ek)
		for(int ei:ej){
			int b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	public static long getAbsMin(long[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long ei:array){
			long b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static long getAbsMin(long[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long[] ej:array)
		for(long ei:ej){
			long b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static long getAbsMin(long[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long[][] ek:array)
		for(long[] ej:ek)
		for(long ei:ej){
			long b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static long getAbsMin(long[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long[][][] el:array)
		for(long[][] ek:el)
		for(long[] ej:ek)
		for(long ei:ej){
			long b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	public static float getAbsMin(float[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float ei:array){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static float getAbsMin(float[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static float getAbsMin(float[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static float getAbsMin(float[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	public static float getAbsMin(float[] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float ei:array) if(ei!=undef){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static float getAbsMin(float[][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej) if(ei!=undef){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static float getAbsMin(float[][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static float getAbsMin(float[][][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef){
			float b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	public static double getAbsMin(double[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double ei:array){
			double b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
				
		return tmp;
	}
	
	public static double getAbsMin(double[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double[] ej:array)
		for(double ei:ej){
			double b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
			
		return tmp;
	}
	
	public static double getAbsMin(double[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double[][] ek:array)
		for(double[] ej:ek)
		for(double ei:ej){
			double b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	public static double getAbsMin(double[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double[][][] el:array)
		for(double[][] ek:el)
		for(double[] ej:ek)
		for(double ei:ej){
			double b=(ei<0)?(-ei):ei;
			if(b<tmp) tmp=b;
		}
		
		return tmp;
	}
	
	
	/**
     * get maximum of values in a array
     *
     * @param	array	an array
     *
     * @return	the maximum value
     *
     * @exception	if array length is zero
     */
	public static int getMin(int[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int ei:array) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static int getMin(int[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int[] ej:array)
		for(int ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static int getMin(int[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int[][] ek:array)
		for(int[] ej:ek)
		for(int ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static int getMin(int[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int tmp=Integer.MAX_VALUE;
		
		for(int[][][] el:array)
		for(int[][] ek:el)
		for(int[] ej:ek)
		for(int ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static long getMin(long[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long ei:array) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static long getMin(long[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long[] ej:array)
		for(long ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static long getMin(long[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long[][] ek:array)
		for(long[] ej:ek)
		for(long ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static long getMin(long[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long tmp=Long.MAX_VALUE;
		
		for(long[][][] el:array)
		for(long[][] ek:el)
		for(long[] ej:ek)
		for(long ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static float getMin(float[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float ei:array) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMin(float[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMin(float[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMin(float[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static float getMin(float[] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float ei:array) if(ei!=undef&&ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMin(float[][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej) if(ei!=undef&&ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMin(float[][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef&&ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static float getMin(float[][][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float tmp=Float.MAX_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef&&ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	
	public static double getMin(double[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double ei:array) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static double getMin(double[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double[] ej:array)
		for(double ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static double getMin(double[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double[][] ek:array)
		for(double[] ej:ek)
		for(double ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	public static double getMin(double[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double tmp=Double.MAX_VALUE;
		
		for(double[][][] el:array)
		for(double[][] ek:el)
		for(double[] ej:ek)
		for(double ei:ej) if(ei<tmp) tmp=ei;
		
		return tmp;
	}
	
	
	/**
     * get the range of the array (max-min)
     *
     * @param	array	an array
     *
     * @return	range of the array
     */
	public static int getRange(int[] array){
		int[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static int getRange(int[][] array){
		int[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static int getRange(int[][][] array){
		int[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static int getRange(int[][][][] array){
		int[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	
	public static long getRange(long[] array){
		long[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static long getRange(long[][] array){
		long[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static long getRange(long[][][] array){
		long[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static long getRange(long[][][][] array){
		long[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	
	public static float getRange(float[] array){
		float[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static float getRange(float[][] array){
		float[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static float getRange(float[][][] array){
		float[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static float getRange(float[][][][] array){
		float[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	
	public static double getRange(double[] array){
		double[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static double getRange(double[][] array){
		double[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static double getRange(double[][][] array){
		double[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	public static double getRange(double[][][][] array){
		double[] ex=getExtrema(array);
		
		return ex[1]-ex[0];
	}
	
	
	/**
     * get minimum and maximum of values in a array
     *
     * @param	array	an array
     *
     * @return	the extreme values, [0] is minimum and [1] is maximum
     */
	public static int[] getExtrema(int[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		int[] re=new int[2];
		
		re[0]=Integer.MAX_VALUE;
		re[1]=Integer.MIN_VALUE;
		
		for(int ei:array){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static int[] getExtrema(int[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int[] re=new int[2];
		
		re[0]=Integer.MAX_VALUE;
		re[1]=Integer.MIN_VALUE;
		
		for(int[] ej:array)
		for(int ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static int[] getExtrema(int[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int[] re=new int[2];
		
		re[0]=Integer.MAX_VALUE;
		re[1]=Integer.MIN_VALUE;
		
		for(int[][] ek:array)
		for(int[] ej:ek)
		for(int ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static int[] getExtrema(int[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		int[] re=new int[2];
		
		re[0]=Integer.MAX_VALUE;
		re[1]=Integer.MIN_VALUE;
		
		for(int[][][] el:array)
		for(int[][] ek:el)
		for(int[] ej:ek)
		for(int ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	
	public static long[] getExtrema(long[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		long[] re=new long[2];
		
		re[0]=Long.MAX_VALUE;
		re[1]=Long.MIN_VALUE;
		
		for(long ei:array){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static long[] getExtrema(long[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long[] re=new long[2];
		
		re[0]=Long.MAX_VALUE;
		re[1]=Long.MIN_VALUE;
		
		for(long[] ej:array)
		for(long ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static long[] getExtrema(long[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long[] re=new long[2];
		
		re[0]=Long.MAX_VALUE;
		re[1]=Long.MIN_VALUE;
		
		for(long[][] ek:array)
		for(long[] ej:ek)
		for(long ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static long[] getExtrema(long[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		long[] re=new long[2];
		
		re[0]=Long.MAX_VALUE;
		re[1]=Long.MIN_VALUE;
		
		for(long[][][] el:array)
		for(long[][] ek:el)
		for(long[] ej:ek)
		for(long ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}	
	
	
	public static float[] getExtrema(float[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float ei:array){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static float[] getExtrema(float[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static float[] getExtrema(float[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static float[] getExtrema(float[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	
	public static float[] getExtrema(float[] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float ei:array) if(ei!=undef){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static float[] getExtrema(float[][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float[] ej:array)
		for(float ei:ej) if(ei!=undef){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static float[] getExtrema(float[][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float[][] ek:array)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static float[] getExtrema(float[][][][] array,float undef){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		float[] re=new float[2];
		
		re[0]=Float.MAX_VALUE;
		re[1]=Float.MIN_VALUE;
		
		for(float[][][] el:array)
		for(float[][] ek:el)
		for(float[] ej:ek)
		for(float ei:ej) if(ei!=undef){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	
	public static double[] getExtrema(double[] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		
		double[] re=new double[2];
		
		re[0]=Double.MAX_VALUE;
		re[1]=Double.MIN_VALUE;
		
		for(double ei:array){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static double[] getExtrema(double[][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double[] re=new double[2];
		
		re[0]=Double.MAX_VALUE;
		re[1]=Double.MIN_VALUE;
		
		for(double[] ej:array)
		for(double ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static double[] getExtrema(double[][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double[] re=new double[2];
		
		re[0]=Double.MAX_VALUE;
		re[1]=Double.MIN_VALUE;
		
		for(double[][] ek:array)
		for(double[] ej:ek)
		for(double ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	public static double[] getExtrema(double[][][][] array){
		if(array.length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0].length==0) throw new IllegalArgumentException("0-length of array");
		if(array[0][0][0].length==0) throw new IllegalArgumentException("0-length of array");
		
		double[] re=new double[2];
		
		re[0]=Double.MAX_VALUE;
		re[1]=Double.MIN_VALUE;
		
		for(double[][][] el:array)
		for(double[][] ek:el)
		for(double[] ej:ek)
		for(double ei:ej){
			if(ei<re[0]) re[0]=ei;
			if(ei>re[1]) re[1]=ei;
		}
		
		return re;
	}
	
	
	/**
     * compress data from float to short or from double to integer
     * using scale and offset
     *
     * @param	data	an array of data
     *
     * @return	re		in short type after compression
     */
	public static short[] compressToShort(float[] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		short[] re=new short[data.length];
		
		for(int i=0,I=data.length;i<I;i++)
		re[i]=compressToShort(data[i],param[0],param[1]);
		
		return re;
	}
	
	public static short[] compressToShort(double[] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		short[] re=new short[data.length];
		
		for(int i=0,I=data.length;i<I;i++)
		re[i]=compressToShort(data[i],param[0],param[1]);
		
		return re;
	}
	
	public static int[] compressToInteger(double[] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		int[] re=new int[data.length];
		
		for(int i=0,I=data.length;i<I;i++)
		re[i]=compressToInteger(data[i],param[0],param[1]);
		
		return re;
	}
	
	
	public static short[][] compressToShort(float[][] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		int J=data.length,I=data[0].length;
		
		short[][] re=new short[J][I];
		
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[j][i]=compressToShort(data[j][i],param[0],param[1]);
		
		return re;
	}
	
	public static short[][] compressToShort(double[][] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		int J=data.length,I=data[0].length;
		
		short[][] re=new short[J][I];
		
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[j][i]=compressToShort(data[j][i],param[0],param[1]);
		
		return re;
	}
	
	public static int[][] compressToInteger(double[][] data){
		double[] param=getScaleAndOffsetToInteger(getExtrema(data));
		
		int J=data.length,I=data[0].length;
		
		int[][] re=new int[J][I];
		
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[j][i]=compressToInteger(data[j][i],param[0],param[1]);
		
		return re;
	}

	
	public static short[][][] compressToShort(float[][][] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		int K=data.length,J=data[0].length,I=data[0][0].length;
		
		short[][][] re=new short[K][J][I];
		
		for(int k=0;k<K;k++)
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[k][j][i]=compressToShort(data[k][j][i],param[0],param[1]);
		
		return re;
	}
	
	public static short[][][] compressToShort(double[][][] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		int K=data.length,J=data[0].length,I=data[0][0].length;
		
		short[][][] re=new short[K][J][I];
		
		for(int k=0;k<K;k++)
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[k][j][i]=compressToShort(data[k][j][i],param[0],param[1]);
		
		return re;
	}
	
	public static int[][][] compressToInteger(double[][][] data){
		double[] param=getScaleAndOffsetToInteger(getExtrema(data));
		
		int K=data.length,J=data[0].length,I=data[0][0].length;
		
		int[][][] re=new int[K][J][I];
		
		for(int k=0;k<K;k++)
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[k][j][i]=compressToInteger(data[k][j][i],param[0],param[1]);
		
		return re;
	}
	

	public static short[][][][] compressToShort(float[][][][] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		int L=data.length,K=data[0].length,J=data[0][0].length,I=data[0][0][0].length;
		
		short[][][][] re=new short[L][K][J][I];
		
		for(int l=0;l<L;l++)
		for(int k=0;k<K;k++)
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[l][k][j][i]=compressToShort(data[l][k][j][i],param[0],param[1]);
		
		return re;
	}
	
	public static short[][][][] compressToShort(double[][][][] data){
		double[] param=getScaleAndOffsetToShort(getExtrema(data));
		
		int L=data.length,K=data[0].length,J=data[0][0].length,I=data[0][0][0].length;
		
		short[][][][] re=new short[L][K][J][I];
		
		for(int l=0;l<L;l++)
		for(int k=0;k<K;k++)
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[l][k][j][i]=compressToShort(data[l][k][j][i],param[0],param[1]);
		
		return re;
	}
	
	public static int[][][][] compressToInteger(double[][][][] data){
		double[] param=getScaleAndOffsetToInteger(getExtrema(data));
		
		int L=data.length,K=data[0].length,J=data[0][0].length,I=data[0][0][0].length;
		
		int[][][][] re=new int[L][K][J][I];
		
		for(int l=0;l<L;l++)
		for(int k=0;k<K;k++)
		for(int j=0;j<J;j++)
		for(int i=0;i<I;i++)
		re[l][k][j][i]=compressToInteger(data[l][k][j][i],param[0],param[1]);
		
		return re;
	}
	
	
	/*** helper methods ***/
	
	/**
     * compress data from double to short or to integer
     * using scale and offset
     * 
     * Reference:
     * http://www.unidata.ucar.edu/software/netcdf/docs/BestPractices.html#Packed%20Data%20Values
     *
     * @param	data	data value
     *
     * @return	re		result in short format
     */
	private static short compressToShort(double data,double scale,double offset){
		long tmp=Math.round((data-offset)/scale);
		short re=(short)tmp;
		
		if(re!=tmp) throw new IllegalArgumentException("overflow when compress data");
		
		return re;
	}
	
	private static int compressToInteger(double data,double scale,double offset){
		long tmp=Math.round((data-offset)/scale);
		int re=(int)tmp;
		
		if(re!=tmp) throw new IllegalArgumentException("overflow when compress data");
		
		return re;
	}
	
	
	/**
     * compute scale and offset for compressing data
     * from float to short or from double to integer (short)
     * 
     * Reference:
     * http://www.unidata.ucar.edu/software/netcdf/docs/BestPractices.html#Packed%20Data%20Values
     *
     * @param	ex	extreme values, [0] is min and [1] is max
     *
     * @return	re	[0] is scale_factor and [1] is add_offset
     */
	private static double[] getScaleAndOffsetToShort(float[] ex){
		final short bits=16;
		
		int de=2<<bits-1;
		
		double scale_factor=(double)(ex[1]-ex[0])/de;
		double add_offset  =ex[0]+de*scale_factor;
		
		return new double[]{scale_factor,add_offset};
	}
	
	private static double[] getScaleAndOffsetToShort(double[] ex){
		final short bits=16;
		
		int de=2<<bits-1;
		
		double scale_factor=(ex[1]-ex[0])/de;
		double add_offset  =ex[0]+de*scale_factor;
		
		return new double[]{scale_factor,add_offset};
	}
	
	private static double[] getScaleAndOffsetToInteger(double[] ex){
		final short bits=32;
		
		long de=2<<bits-1;
		
		double scale_factor=(ex[1]-ex[0])/de;
		double add_offset  =ex[0]+de*scale_factor;
		
		return new double[]{scale_factor,add_offset};
	}
	
	
	/**
     * similar to Arrays.binarySearch() but the array is in decreasing order
     */
	private static int binarySearchDecre(int[] a,int key){
		int low=0;
		int high=a.length-1;
		
		while(low<=high){
			int mid=(low+high)>>> 1;
			int midVal=a[mid];
			
			if(midVal>key) low=mid+1;		// Neither val is NaN, thisVal is smaller
			else if(midVal<key) high=mid-1;	// Neither val is NaN, thisVal is larger
			else return mid;	// key found
		}
		
		return -low;  // key not found.
	}
	
	private static int binarySearchDecre(short[] a,short key){
		int low=0;
		int high=a.length-1;
		
		while(low<=high){
			int mid=(low+high)>>> 1;
			int midVal=a[mid];
			
			if(midVal>key) low=mid+1;		// Neither val is NaN, thisVal is smaller
			else if(midVal<key) high=mid-1;	// Neither val is NaN, thisVal is larger
			else return mid;	// key found
		}
		
		return -low;  // key not found.
	}
	
	private static int binarySearchDecre(long[] a,long key){
		int low=0;
		int high=a.length-1;
		
		while(low<=high){
			int mid=(low+high)>>> 1;
			long midVal=a[mid];
			
			if(midVal>key) low=mid+1;		// Neither val is NaN, thisVal is smaller
			else if(midVal<key) high=mid-1;	// Neither val is NaN, thisVal is larger
			else return mid;	// key found
		}
		
		return -low;  // key not found.
	}
	
	private static int binarySearchDecre(float[] a,float key){
		int low=0;
		int high=a.length-1;
		
		while(low<=high){
			int mid=(low+high)>>> 1;
			float midVal=a[mid];
			
			if(midVal>key) low=mid+1;		// Neither val is NaN, thisVal is smaller
			else if(midVal<key) high=mid-1;	// Neither val is NaN, thisVal is larger
			else{
				int midBits=Float.floatToIntBits(midVal);
				int keyBits=Float.floatToIntBits(key);
				
				if(midBits==keyBits) return mid;		// Values are equal, Key found
				else if(midBits<keyBits) low = mid+1;	// (-0.0, 0.0) or (!NaN, NaN)
				else high=mid-1;						// (0.0, -0.0) or (NaN, !NaN)
			}
		}
		
		return -low;  // key not found.
	}
	
	private static int binarySearchDecre(double[] a,double key){
		int low=0;
		int high=a.length-1;
		
		while(low<=high){
			int mid=(low+high)>>> 1;
			double midVal=a[mid];
			
			if(midVal>key) low=mid+1;		// Neither val is NaN, thisVal is smaller
			else if(midVal<key) high=mid-1;	// Neither val is NaN, thisVal is larger
			else{
				long midBits=Double.doubleToLongBits(midVal);
				long keyBits=Double.doubleToLongBits(key);
				
				if(midBits==keyBits) return mid;		// Values are equal, Key found
				else if(midBits<keyBits) low = mid+1;	// (-0.0, 0.0) or (!NaN, NaN)
				else high=mid-1;						// (0.0, -0.0) or (NaN, !NaN)
			}
		}
		
		return -low;  // key not found.
	}
	
	
	/** test
	public static void main(String[] args){
		double[] array=new double[]{4,6,7,9,12,17,100};
		//double[] array=new double[]{100,17,12,9,7,6,4};
		
		System.out.println(getLEIdxIncre(array,100.1));
	}*/
}
