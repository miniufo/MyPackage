/**
 * @(#)BWLabel.java	1.0 2016.12.27
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import java.util.Arrays;
import org.ejml.data.FMatrixRMaj;


/**
 * Labelling connected components in binary 2D array
 *
 * @version 1.0, 2013.03.12
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class BWLabel{
	//
	private int y=0;
	private int x=0;
	
	private int maxLabel=0;
	
	private float[][] data=null;
	private float[][] labs=null;
	
	
	/**
	 * constructor
	 * 
	 * @param	data	2D array containing 1 or 0 only
     */
	public BWLabel(float[][] data){
		this.data=data;
		this.y=data.length;
		this.x=data[0].length;
		
		labs=new float[y][x];
		
		checkBWMatrix(data);
	}
	
	
	/**
	 * connected-components labelling
	 * 
	 * @param	n	connectivity of 4 or 8
     */
	public void connComponentlabelling(int n){
		maxLabel=0; clearLabel();
		
		if(n==4){
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(labs[j][i]==0&&data[j][i]==1) dealBWLabel4(j,i,++maxLabel);
			
		}else if(n==8){
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(labs[j][i]==0&&data[j][i]==1) dealBWLabel8(j,i,++maxLabel);
			
		}else throw new IllegalArgumentException("n ("+n+") should be 4 or 8");
	}
	
	
	/*** getor and setor ***/
	public int getMaxLabel(){ return maxLabel;}
	
	public int[] getLabelsSortedByCount(){
		if(maxLabel==0) return new int[]{0,0};
		
		int[] labels=new int[maxLabel];
		
		int[][] tmp=new int[maxLabel][2];
		
		for(int l=1;l<=maxLabel;l++){
			int count=0;
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(labs[j][i]==l) count++;
			
			tmp[l-1][0]=l; tmp[l-1][1]=count;
		}
		
		Arrays.sort(tmp,(int[] a,int[] b)->{ return Integer.compare(b[1],a[1]);});
		
		for(int l=0;l<maxLabel;l++) labels[l]=tmp[l][0];
		
		return labels;
	}
	
	public float[][] getLabelData(){ return labs;}
	
	
	/*** helper method ***/
	private void checkBWMatrix(float[][] data){
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(data[j][i]!=0&&data[j][i]!=1)
		throw new IllegalArgumentException("data should contain 0 and 1 only");
	}
	
	private void clearLabel(){
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) labs[j][i]=0;
	}
	
	private void dealBWLabel4(int j,int i,int label){
		if(i<0||i==x) return;
		if(j<0||j==y) return;
		if(labs[j][i]!=0||data[j][i]!=1) return;
		
		labs[j][i]=label;
		
		dealBWLabel4(j+1,i,label);
		dealBWLabel4(j,i+1,label);
		dealBWLabel4(j-1,i,label);
		dealBWLabel4(j,i-1,label);
	}
	
	private void dealBWLabel8(int j,int i,int label){
		if(i<0||i==x) return;
		if(j<0||j==y) return;
		if(labs[j][i]!=0||data[j][i]!=1) return;
		
		labs[j][i]=label;
		
		dealBWLabel8(j+1,i  ,label);
		dealBWLabel8(j  ,i+1,label);
		dealBWLabel8(j-1,i  ,label);
		dealBWLabel8(j  ,i-1,label);
		dealBWLabel8(j-1,i-1,label);
		dealBWLabel8(j-1,i+1,label);
		dealBWLabel8(j+1,i-1,label);
		dealBWLabel8(j+1,i+1,label);
	}
	
	
	/**
	 * used to print out the 
	 */
	public void print(String format){
		FMatrixRMaj d=new FMatrixRMaj(data);
		FMatrixRMaj l=new FMatrixRMaj(labs);
		
		d.print(format);
		l.print(format);
	}
	
	
	/** test
	public static void main(String[] args){
		int y=8,x=2;
		
		float[][] data=new float[y][x];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) data[j][i]=Math.round(Math.random());
		
		BWLabel bw=new BWLabel(data);
		
		bw.connComponentlabelling(4); bw.print("%2.0f");
		System.out.println("there are "+bw.getMaxLabel()+" labels for 4\n\n");
		//System.out.println(Arrays.toString(bw.getLabelsInDescendingNumber()));
		
		bw.connComponentlabelling(8); bw.print("%2.0f"); System.out.println("there are "+bw.getMaxLabel()+" labels for 8");
		//System.out.println(Arrays.toString(bw.getLabelsInDescendingNumber()));
	}*/
}
