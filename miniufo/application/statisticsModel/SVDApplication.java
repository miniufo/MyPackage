/**
 * @(#)SVDApplication.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.SingularValueDecomposition_F64;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;


/**
 * SVD application
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class SVDApplication extends StatisticsApplication{
	
	/**
     * singular value decomposition
     *
     * @param	v1	a given Variable
     * @param	v2	another given Variable
     * @param	mc	mode count which need to output
     *
     * @return	re	[0][m] stands for the mth modes of v1, [1][m] stands for the mth temporal serial of v1
     *				[2][m] stands for the mth modes of v2, [3][m] stands for the mth temporal serial of v2
     */
	public static Variable[][] SVD(Variable v1,Variable v2,int mc){
		System.out.println("\nStart singular value decomposition (SVD) analysis...");
		
		if(v1.getTCount()!=v2.getTCount()) throw new IllegalArgumentException("dimensions not same");
		if(v1.getZCount()!=v2.getZCount()) throw new IllegalArgumentException("dimensions not same");
		
		int y1=v1.getYCount(), x1=v1.getXCount(), z=v1.getZCount();
		int y2=v2.getYCount(), x2=v2.getXCount(), t=v1.getTCount();
		int p,q,cc=0,g1=x1*y1,g2=x2*y2,scount=0,ucount=0;
		
		if(z > 1) throw new IllegalArgumentException("only one level allowed when using this method");
		if(t < 2) throw new IllegalArgumentException("t-dimension length is not enough");
		if(mc<=0) throw new IllegalArgumentException("mode count should be larger than 0");
		if(v1.isTFirst()||v2.isTFirst()) throw new IllegalArgumentException("T is first dimension of array");
		
		float undef=v1.getUndef(),svsum=0;
		
		float[][][][] data1=v1.getData();
		float[][][][] data2=v2.getData();
		
		// record the undef grid
		for(int j=0;j<y1;j++)
		for(int i=0;i<x1;i++) if(data1[0][j][i][0]==undef) ucount++;
		
		p=g1-ucount;	ucount=0;
		
		for(int j=0;j<y2;j++)
		for(int i=0;i<x2;i++) if(data2[0][j][i][0]==undef) ucount++;
		
		q=g2-ucount;
		
		if(mc>Math.min(p,q))
			throw new IllegalArgumentException("output modes are more than the matrix could produce");
		
		//
		Variable[][] re=new Variable[4][mc];
		
		for(int m=0;m<mc;m++){
			re[0][m]=new Variable(v1.getName()+(m+1),false,new Range(1,1,y1,x1));	re[0][m].setValue(undef);
			re[1][m]=new Variable(v1.getName()+(m+1),false,new Range(t,1,1,1));
			re[2][m]=new Variable(v2.getName()+(m+1),false,new Range(1,1,y2,x2));	re[2][m].setValue(undef);
			re[3][m]=new Variable(v2.getName()+(m+1),false,new Range(t,1,1,1));
			
			re[0][m].setCommentAndUnit("SVD "+(m+1)+" modes of "+v1.getName());
			re[1][m].setCommentAndUnit("SVD "+(m+1)+" time series of "+v1.getName());
			re[2][m].setCommentAndUnit("SVD "+(m+1)+" modes of "+v2.getName());
			re[3][m].setCommentAndUnit("SVD "+(m+1)+" time series of "+v2.getName());
			
			re[0][m].setUndef(undef);	re[1][m].setUndef(undef);
			re[2][m].setUndef(undef);	re[3][m].setUndef(undef);
		}
		
		float[][] mtrix1=new float[p][];
		float[][] mtrix2=new float[q][];
		
		// change 4D array to 2D array
		for(int j=0;j<y1;j++)
		for(int i=0;i<x1;i++)
		if(data1[0][j][i][0]!=undef) mtrix1[cc++]=data1[0][j][i];
		
		for(int j=0;j<t;j++)
		for(int i=0;i<p;i++)
		if(mtrix1[i][j]==undef) throw new IllegalArgumentException("Matrix contains undef values");
		
		cc=0;
		
		for(int j=0;j<y2;j++)
		for(int i=0;i<x2;i++)
		if(data2[0][j][i][0]!=undef) mtrix2[cc++]=data2[0][j][i];
		
		for(int j=0;j<t;j++)
		for(int i=0;i<q;i++)
		if(mtrix2[i][j]==undef) throw new IllegalArgumentException("Matrix contains undef values");
		
		// calculate the covariance matrix
		double[][] Sdata=new double[p][q];
		
		for(int i=0;i<p;i++)
		for(int j=0;j<q;j++){
			for(int k=0;k<t;k++) Sdata[i][j]+=mtrix1[i][k]*mtrix2[j][k];
			Sdata[i][j]/=t;
		}
		
		// start svd
		SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(p,q,true,true,true);
		svdF.decompose(new DMatrixRMaj(Sdata));
		
		DMatrixRMaj U=svdF.getU(null,false);
		DMatrixRMaj V=svdF.getV(null,false);
		
		double[] sv=svdF.getSingularValues();
		
		for(int i=0;i<sv.length;i++) if(sv[i]!=0){ svsum+=sv[i]*sv[i];	scount++;}
		if(scount<mc) throw new IllegalArgumentException("Not enough singular values for output");
		
		float[][][][][] r11data=new float[mc][][][][];
		float[][][][][] r12data=new float[mc][][][][];
		float[][][][][] r21data=new float[mc][][][][];
		float[][][][][] r22data=new float[mc][][][][];
		
		for(int m=0;m<mc;m++){
			r11data[m]=re[0][m].getData();	r12data[m]=re[1][m].getData();
			r21data[m]=re[2][m].getData();	r22data[m]=re[3][m].getData();
		}
		
		// A=U.transpose().multiply(V1);	B=V.transpose().multiply(V2);
		for(int i=0;i<mc;i++){
			for(int j=0;j<t;j++){
				for(int k=0;k<p;k++) r12data[i][0][0][0][j]+=U.get(k,i)*mtrix1[k][j];
				for(int k=0;k<q;k++) r22data[i][0][0][0][j]+=V.get(k,i)*mtrix2[k][j];
			}
			
			// coefficient for Heterogeneous correlation map
			/**/
			double c1;
			
			c1=sv[i]/miniufo.statistics.StatisticsUtil.cStandardDeviation(r22data[i][0][0][0]);
			for(int k=0,K=U.numRows;k<K;k++) U.set(k,i,U.get(k,i)*c1);
			
			c1=sv[i]/miniufo.statistics.StatisticsUtil.cStandardDeviation(r12data[i][0][0][0]);
			for(int k=0,K=V.numRows;k<K;k++) V.set(k,i,V.get(k,i)*c1);
		}
		
		// restore mode data
		cc=0;
		for(int j=0;j<y1;j++)
		for(int i=0;i<x1;i++){
			if(data1[0][j][i][0]!=undef){
				for(int k=0;k<mc;k++) r11data[k][0][j][i][0]=(float)U.get(cc,k);
				
				cc++;
			}
		}
		
		cc=0;
		for(int j=0;j<y2;j++)
		for(int i=0;i<x2;i++){
			if(data2[0][j][i][0]!=undef){
				for(int k=0;k<mc;k++) r21data[k][0][j][i][0]=(float)V.get(cc,k);
				
				cc++;
			}
		}
		
		// output contribution of each mode
		StringBuilder sb=new StringBuilder();
		
		sb.append("\nResult of SVD:\n");
		for(int i=0;i<mc;i++){
			sb.append("Contribution of ");
			sb.append(i+1);
			sb.append(" pair of modes: ");
			sb.append(sv[i]*sv[i]/svsum+"\t"+sv[i]);
			sb.append("\n");
		}
		
		System.out.println(sb.toString());
		
		Range range1=v1.getRange();
		Range range2=v2.getRange();
		
		for(int m=0;m<mc;m++){
			Range range11=re[0][m].getRange();			Range range12=re[1][m].getRange();
			Range range21=re[2][m].getRange();			Range range22=re[3][m].getRange();
			
			// process the range
			range11.setZRange(range1.getZRange()[0]);	range12.setZRange(range1.getZRange()[0]);
			range11.setTRange(range1.getTRange()[0]);	range12.setTRange(range1);
			range11.setXRange(range1);					range12.setXRange(range1.getXRange()[0]);
			range11.setYRange(range1);					range12.setYRange(range1.getYRange()[0]);
			
			range21.setZRange(range2.getZRange()[0]);	range22.setZRange(range2.getZRange()[0]);
			range21.setTRange(range2.getTRange()[0]);	range22.setTRange(range2);
			range21.setXRange(range2);					range22.setXRange(range2.getXRange()[0]);
			range21.setYRange(range2);					range22.setYRange(range2.getYRange()[0]);
		}
		
		System.out.println("Finish SVD.");
		
		return re;
	}
	
	
	/** test
	public static void main(String arg[]){
		try{
			
			
	    }catch(Exception ex){ ex.printStackTrace();}
	}*/
}
