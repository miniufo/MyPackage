/**
 * @(#)Indexing.java	1.0 2018.03.23
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.basic;

import java.util.regex.Pattern;


/**
 * Indexing the Java array as Matlab or Python.
 *
 * @version 1.0, 2018.03.23
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Indexing{
	//
	public static final Pattern colon=Pattern.compile(":"  ,Pattern.CASE_INSENSITIVE);
	public static final Pattern comma=Pattern.compile(","  ,Pattern.CASE_INSENSITIVE);
	public static final Pattern space=Pattern.compile("\\s",Pattern.CASE_INSENSITIVE);
	
	
	/**
	 * Prevent from instantiation
	 */
	private Indexing(){}
	
	
	/**
	 * Indexing a 1D array similar to Matlab or Python.
	 *
	 * @param	data	1D array
	 * @param	index	eg, ":", "1:2:end", "end-1", "[1 3 5]", ...
	 */
	public static int[] get(int[] data,String index){
		int[] idx=getIndices(data.length,space.matcher(index).replaceAll(""));
		
		int[] re=new int[idx.length];
		
		for(int i=0,I=idx.length;i<I;i++) re[i]=data[idx[i]];
		
		return re;
	}
	
	public static long[] get(long[] data,String index){
		int[] idx=getIndices(data.length,space.matcher(index).replaceAll(""));
		
		long[] re=new long[idx.length];
		
		for(int i=0,I=idx.length;i<I;i++) re[i]=data[idx[i]];
		
		return re;
	}
	
	public static float[] get(float[] data,String index){
		int[] idx=getIndices(data.length,space.matcher(index).replaceAll(""));
		
		float[] re=new float[idx.length];
		
		for(int i=0,I=idx.length;i<I;i++) re[i]=data[idx[i]];
		
		return re;
	}
	
	public static double[] get(double[] data,String index){
		int[] idx=getIndices(data.length,space.matcher(index).replaceAll(""));
		
		double[] re=new double[idx.length];
		
		for(int i=0,I=idx.length;i<I;i++) re[i]=data[idx[i]];
		
		return re;
	}
	
	
	/**
	 * Indexing a 2D array similar to Matlab or Python.
	 *
	 * @param	data	1D array
	 * @param	index	eg, ":,1:2:end", "end-1,[1 3 5]", ...
	 */
	public static int[][] get(int[][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=2) throw new IllegalArgumentException("index ("+index+") should contain only one comma");
		
		int[] idx1=getIndices(data[0].length,tokens[1]);
		int[] idx2=getIndices(data.length,tokens[0]);
		
		int[][] re=new int[idx2.length][idx1.length];
		
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[j][i]=data[idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static long[][] get(long[][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=2) throw new IllegalArgumentException("index ("+index+") should contain only one comma");
		
		int[] idx1=getIndices(data[0].length,tokens[1]);
		int[] idx2=getIndices(data.length,tokens[0]);
		
		long[][] re=new long[idx2.length][idx1.length];
		
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[j][i]=data[idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static float[][] get(float[][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=2) throw new IllegalArgumentException("index ("+index+") should contain only one comma");
		
		int[] idx1=getIndices(data[0].length,tokens[1]);
		int[] idx2=getIndices(data.length,tokens[0]);
		
		float[][] re=new float[idx2.length][idx1.length];
		
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[j][i]=data[idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static double[][] get(double[][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=2) throw new IllegalArgumentException("index ("+index+") should contain only one comma");
		
		int[] idx1=getIndices(data[0].length,tokens[1]);
		int[] idx2=getIndices(data.length,tokens[0]);
		
		double[][] re=new double[idx2.length][idx1.length];
		
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[j][i]=data[idx2[j]][idx1[i]];
		
		return re;
	}
	
	
	/**
	 * Indexing a 3D array similar to Matlab or Python.
	 *
	 * @param	data	1D array
	 * @param	index	eg, ":,1:2:end,[1 3 5]", ...
	 */
	public static int[][][] get(int[][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=3) throw new IllegalArgumentException("index ("+index+") should contain only two comma");
		
		int[] idx1=getIndices(data[0][0].length,tokens[2]);
		int[] idx2=getIndices(data[0].length,tokens[1]);
		int[] idx3=getIndices(data.length,tokens[0]);
		
		int[][][] re=new int[idx3.length][idx2.length][idx1.length];
		
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[k][j][i]=data[idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static long[][][] get(long[][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=3) throw new IllegalArgumentException("index ("+index+") should contain only two comma");
		
		int[] idx1=getIndices(data[0][0].length,tokens[2]);
		int[] idx2=getIndices(data[0].length,tokens[1]);
		int[] idx3=getIndices(data.length,tokens[0]);
		
		long[][][] re=new long[idx3.length][idx2.length][idx1.length];
		
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[k][j][i]=data[idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static float[][][] get(float[][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=3) throw new IllegalArgumentException("index ("+index+") should contain only two comma");
		
		int[] idx1=getIndices(data[0][0].length,tokens[2]);
		int[] idx2=getIndices(data[0].length,tokens[1]);
		int[] idx3=getIndices(data.length,tokens[0]);
		
		float[][][] re=new float[idx3.length][idx2.length][idx1.length];
		
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[k][j][i]=data[idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static double[][][] get(double[][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=3) throw new IllegalArgumentException("index ("+index+") should contain only two comma");
		
		int[] idx1=getIndices(data[0][0].length,tokens[2]);
		int[] idx2=getIndices(data[0].length,tokens[1]);
		int[] idx3=getIndices(data.length,tokens[0]);
		
		double[][][] re=new double[idx3.length][idx2.length][idx1.length];
		
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[k][j][i]=data[idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	
	/**
	 * Indexing a 4D array similar to Matlab or Python.
	 *
	 * @param	data	1D array
	 * @param	index	eg, ":,1:2:end,[1 3 5],end-1", ...
	 */
	public static int[][][][] get(int[][][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=4) throw new IllegalArgumentException("index ("+index+") should contain only three comma");
		
		int[] idx1=getIndices(data[0][0][0].length,tokens[3]);
		int[] idx2=getIndices(data[0][0].length,tokens[2]);
		int[] idx3=getIndices(data[0].length,tokens[1]);
		int[] idx4=getIndices(data.length,tokens[0]);
		
		int[][][][] re=new int[idx4.length][idx3.length][idx2.length][idx1.length];
		
		for(int l=0,L=idx3.length;l<L;l++)
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[l][k][j][i]=data[idx4[l]][idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static long[][][][] get(long[][][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=4) throw new IllegalArgumentException("index ("+index+") should contain only three comma");
		
		int[] idx1=getIndices(data[0][0][0].length,tokens[3]);
		int[] idx2=getIndices(data[0][0].length,tokens[2]);
		int[] idx3=getIndices(data[0].length,tokens[1]);
		int[] idx4=getIndices(data.length,tokens[0]);
		
		long[][][][] re=new long[idx4.length][idx3.length][idx2.length][idx1.length];
		
		for(int l=0,L=idx3.length;l<L;l++)
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[l][k][j][i]=data[idx4[l]][idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static float[][][][] get(float[][][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=4) throw new IllegalArgumentException("index ("+index+") should contain only three comma");
		
		int[] idx1=getIndices(data[0][0][0].length,tokens[3]);
		int[] idx2=getIndices(data[0][0].length,tokens[2]);
		int[] idx3=getIndices(data[0].length,tokens[1]);
		int[] idx4=getIndices(data.length,tokens[0]);
		
		float[][][][] re=new float[idx4.length][idx3.length][idx2.length][idx1.length];
		
		for(int l=0,L=idx3.length;l<L;l++)
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[l][k][j][i]=data[idx4[l]][idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	public static double[][][][] get(double[][][][] data,String index){
		String[] tokens=comma.split(space.matcher(index).replaceAll(""));
		
		if(tokens.length!=4) throw new IllegalArgumentException("index ("+index+") should contain only three comma");
		
		int[] idx1=getIndices(data[0][0][0].length,tokens[3]);
		int[] idx2=getIndices(data[0][0].length,tokens[2]);
		int[] idx3=getIndices(data[0].length,tokens[1]);
		int[] idx4=getIndices(data.length,tokens[0]);
		
		double[][][][] re=new double[idx4.length][idx3.length][idx2.length][idx1.length];
		
		for(int l=0,L=idx3.length;l<L;l++)
		for(int k=0,K=idx3.length;k<K;k++)
		for(int j=0,J=idx2.length;j<J;j++)
		for(int i=0,I=idx1.length;i<I;i++) re[l][k][j][i]=data[idx4[l]][idx3[k]][idx2[j]][idx1[i]];
		
		return re;
	}
	
	
	/**
	 * Get indices along one dimension.
	 *
	 * @param	dlen	total length of the data (1D array)
	 * @param	index	eg, ":", "1:2:end", "[1 3 5]", "end-1", ...
	 */
	public static int[] getIndices(int dlen,String index){
		if(index.indexOf(":")==-1)
			return parseLevelIndex (dlen,index);
		else
			return parseLinearIndex(dlen,index);
	}
	
	
	/*** helper methods ***/
	
	/**
	 * Parse a single index along one dimension.
	 * 
	 * Return an array specifying start-index, stride, and end-index.
	 * 
	 * @param	dlen	length of a given array
	 * @param	oneIdx	one index such as "1:end", "2:3:end", ":"...
	 */
	private static int[] parseLinearIndex(int dlen,String oneIdx){
		if(":".equals(oneIdx)||"".equals(oneIdx)){
			int[] idx=new int[dlen];
			
			for(int i=0;i<dlen;i++) idx[i]=i;
			
			return idx;
		}
		
		int strIdx=0;
		int stride=1;
		int endIdx=dlen-1;
		
		String[] parts=colon.split(oneIdx,-1);
		
		int len=parts.length;
		
		switch(len){
		case 1:
			strIdx=Integer.parseInt(parts[0]);
			endIdx=strIdx;
			break;
			
		case 2:
			strIdx=Integer.parseInt(parts[0]);
			endIdx=parseEndIdx(dlen,parts[1]);
			
			break;
			
		case 3:
			strIdx=Integer.parseInt(parts[0]);
			stride=Integer.parseInt(parts[1]);
			endIdx=parseEndIdx(dlen,parts[2]);
			
			break;
			
		default: throw new IllegalArgumentException("invalid index: "+oneIdx);
		}
		
		int rlen=(endIdx-strIdx)/stride+1;
		
		int[] idx=new int[rlen];
		
		for(int i=strIdx,ii=0;i<=endIdx;i+=stride,ii++) idx[ii]=i;
		
		return idx;
	}
	
	/**
	 * Parse a single index along one dimension.
	 * 
	 * Return an array specifying every index (in square bracket if more
	 * than two) that are required.
	 * 
	 * @param	dlen	length of a given array
	 * @param	oneIdx	one index such as "1", "[2 3 end]"...
	 */
	private static int[] parseLevelIndex(int dlen,String oneIdx){
		int lBracket=oneIdx.indexOf("[");
		int rBracket=oneIdx.indexOf("]");
		
		int[] idx=null;
		
		if(lBracket!=-1){
			String[] tokens=space.split(oneIdx.substring(lBracket+1,rBracket));
			
			idx=new int[tokens.length];
			
			for(int i=0,I=tokens.length;i<I;i++){
				if(tokens[i].indexOf("end")!=-1)
					idx[i]=parseEndIdx(dlen,tokens[i]);
				else
					idx[i]=Integer.parseInt(tokens[i]);
			}
			
		}else{
			if(oneIdx.indexOf("end")!=-1)
				idx=new int[]{parseEndIdx(dlen,oneIdx)};
			else
				idx=new int[]{Integer.parseInt(oneIdx)};
		}
		
		return idx;
	}
	
	/**
	 * Parse the end index along one dimension.
	 * 
	 * @param	dlen	length of a given array
	 * @param	last	last index such as "end", "end-3"...
	 */
	private static int parseEndIdx(int dlen,String last){
		int endPos=last.indexOf("end");
		int endIdx;
		
		if(endPos!=-1){
			String remain=last.substring(endPos+3).trim();
			
			endIdx=dlen-1;
			
			if(!"".equals(remain)) endIdx+=Integer.parseInt(remain);
			
		}else endIdx=Integer.parseInt(last);
		
		return endIdx;
	}
	
	
	/** test 
	public static void main(String[] args){
		float[][] data=new float[][]{
			{ 1, 2, 3, 4, 5, 6, 7},
			{ 2, 3, 4, 5, 6, 7, 8},
			{ 3, 4, 5, 6, 7, 8, 9},
			{ 4, 5, 6, 7, 8, 9,10}
		};
		
		float[][] sub=Indexing.get(data,":,1:2:end-1");
		
		for(float[] a:sub) System.out.println(Arrays.toString(a));
	}*/
}
