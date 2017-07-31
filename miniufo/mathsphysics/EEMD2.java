/**
 * @(#)EmpiricalModeDecomposition.java	1.0 2013.06.06
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;


/**
 * Empirical model decomposition
 *
 * @version 1.0, 2013.06.06
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EEMD2{
	//
	private int xcount=0;
	private int ycount=0;
	private int nmode =0;
	
	private float Nstd=0;
	
	private EEMD emd=null;
	
	
	/**
	 * constructor
	 * 
	 * @param	len		length of data
	 * @param	nmode	number of modes to be decomposed
	 * @param	Nstd	ratio of the std. of the added noise and that of data
	 */
	public EEMD2(int xcount,int ycount,int nmode,float Nstd){
		if(xcount<3) throw new IllegalArgumentException("x-count should be at least 3");
		if(ycount<3) throw new IllegalArgumentException("y-count should be at least 3");
		
		this.xcount=xcount;
		this.ycount=ycount;
		this.nmode =nmode;
		this.Nstd  =Nstd;
	}
	
	
	/*** getor and setor ***/
	public int getModes(){ return emd.getModes();}
	
	public int getMaxIteration(){ return emd.getMaxIteration();}
	
	public int getMaxEnsemble(){ return emd.getMaxEnsemble();}
	
	public void setMaxIteration(int itMax){ emd.setMaxIteration(itMax);}
	
	public void setMaxEnsemble(int enMax){ emd.setMaxEnsemble(enMax);}
	
	
	/**
	 * decomposition of a given data
	 */
	public float[][][] decomp(float[][] data){
		float[][][] allm=new float[nmode+1][ycount][xcount];
		
		decomp(data,allm);
		
		return allm;
	}
	
	public void decomp(float[][] data,float[][][] re){
		if(re.length!=nmode+1||re[0].length!=ycount||re[0][0].length!=xcount)
		throw new IllegalArgumentException("invalid result sizes");
		
		float[][] tmpx=new float[nmode+2][xcount];
		float[][] tmpy=new float[nmode+2][ycount];
		float[][][][] tmp=new float[nmode+1][nmode+1][ycount][xcount];
		
		// decompose in the first direction
		emd=new EEMD(xcount,nmode,Nstd);
		for(int j=0;j<ycount;j++){System.out.println("j: "+j);
			emd.decomp(data[j],tmpx);
			
			for(int m=0,M=nmode+1;m<M;m++) System.arraycopy(tmpx[m+1],0,re[m][j],0,xcount);
		}
		
		// decompose in the second direction
		emd=new EEMD(ycount,nmode,Nstd);
		float[] buf=new float[ycount];
		for(int i=0;i<xcount;i++){System.out.println("i: "+i);
		for(int m=0,M=nmode+1;m<M;m++){
			for(int j=0;j<ycount;j++) buf[j]=re[m][j][i];
			
			emd.decomp(buf,tmpy);
			
			for(int j=0;j<ycount;j++)
			for(int mm=0,MM=nmode+1;mm<MM;mm++) tmp[mm][m][j][i]=tmpy[mm+1][j];
		}}
		
		// combine modes
		for(int j=0;j<ycount;j++)
		for(int i=0;i<xcount;i++){
			for(int m=0,M=nmode+1;m<M;m++){
				re[m][j][i]=0;
				
				for(int mm=m,MM=nmode+1;mm<MM;mm++)
				re[m][j][i]+=tmp[m][mm][j][i]+tmp[mm][m][j][i];
				
				re[m][j][i]-=tmp[m][m][j][i];
			}
		}
	}
	
	
	/** test
	public static void main(String[] args){
		int xcount=240;
		int ycount=155;
		
		int[] cols=new int[ycount];
		for(int l=0;l<ycount;l++) cols[l]=l+1;
		
		float[][] data=TextReader.readColumnsF("d:/uwrf.txt",false,cols);
		
		EEMD2 emd=new EEMD2(xcount,ycount,5,0.02f);
		float[][][] re=emd.decomp(data);
		
		Variable[] modes=new Variable[re.length];
		
		for(int i=0,I=re.length;i<I;i++){
			modes[i]=new Variable("m"+i,new Range(1,1,ycount,xcount));
			modes[i].setCommentAndUnit("mode "+i);
			
			float[][] mdata=modes[i].getData()[0][0];
			
			for(int j=0;j<ycount;j++)
			System.arraycopy(re[i][j],0,mdata[j],0,xcount);
		}
		
		float[][] mdata=modes[re.length].getData()[0][0];
		
		for(int j=0;j<ycount;j++)
		System.arraycopy(data[j],0,mdata[j],0,xcount);
		
		CtlDataWriteStream dw=new CtlDataWriteStream("d:/EEMD.dat");
		dw.writeData(modes);
		dw.closeFile();
	}*/
}
