/**
 * @(#)EOFApplication.java	1.0 07/02/01
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
import miniufo.mathsphysics.VerticalModeDecomposition;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.pow;
import static java.lang.Math.atan2;


/**
 * EOF application
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class EOFApplication extends StatisticsApplication{
	
	/**
     * constructor
     */
	private EOFApplication(){}
	
	
	/**
     * empirical orthogonal decomposition
     *
     * @param	v	a given Variable
     * @param	mc	mode count need to be restored
     * @param	ec	eigenvalue count need to be restored
     * 
     * @return	result	result of EOF
     */
	public static EOFResult EOF(Variable v,int mc,int ec){
		EOFResult result=new EOFResult(mc,ec,v);
		
		int t=v.getTCount(),y=v.getYCount(),x=v.getXCount();
		
		float undef=v.getUndef();	float[][][] data=v.getData()[0];
		
		int ge=cDefinedGridCount(v);
		
		checkData(v,mc,ge);
		
		float[][] mtrix=new float[ge][];	// 2D matrix
		DMatrixRMaj Vmtrx=null;			// modes
		
		// change 4D array to 2D array
		for(int j=0,cc=0;j<y;j++)
		for(int i=0;i<x;i++) if(data[j][i][0]!=undef) mtrix[cc++]=data[j][i];
		
		checkUndef(mtrix,undef);
		
		double[] sv=null;
		
		// if m > t, then change the time and space
		if(t>=ge){
			// calculate the Matrix of covariance
			double[][] Sdata=new double[ge][ge];
			
			cCovMatrix1(mtrix,Sdata);
			
			SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(ge,ge,true,true,true);
			if(!svdF.decompose(new DMatrixRMaj(Sdata))) throw new IllegalArgumentException("matrix cannot be decomposed");
			Vmtrx=svdF.getU(null,false);
			sv=svdF.getSingularValues();
			
			restoreContribution(cRatio(sv,ec),sv,result);
			
			for(int m=0;m<mc;m++)
			for(int i=0;i<ge;i++) Vmtrx.set(i,m,Vmtrx.get(i,m)*Math.sqrt(sv[m]/t));
			
		}else{
			// calculate the Matrix of covariance
			double[][] Sdata=new double[t][t];
			
			cCovMatrix2(mtrix,Sdata);
			
			SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(ge,ge,true,true,true);
			if(!svdF.decompose(new DMatrixRMaj(Sdata))) throw new IllegalArgumentException("matrix cannot be decomposed");
			sv=svdF.getSingularValues();
			
			// Vmtrx=A.multiply(svd.getU());
			DMatrixRMaj U=svdF.getU(null,false);
			
			restoreContribution(cRatio(sv,ec),sv,result);
			
			Vmtrx=new DMatrixRMaj(ge,mc);
			
			for(int g=0;g<ge;g++)
			for(int m=0;m<mc;m++){
				for(int l=0;l<t;l++) Vmtrx.set(g,m,Vmtrx.get(g,m)+mtrix[g][l]*U.get(l,m));
				Vmtrx.set(g,m,Vmtrx.get(g,m)/Math.sqrt(t));
			}
		}
		
		// restore temporal data
		for(int m=0;m<mc;m++){
			float[][][] r1data=result.getTimes()[m].getData()[0];
			
			// T=V.transpose().multiply(A);	restore temporal data to re[1]
			for(int l=0;l<t;l++){
				for(int g=0;g<ge;g++) r1data[0][0][l]+=Vmtrx.get(g,m)*mtrix[g][l];
				
				r1data[0][0][l]/=sv[m]/t;
			}
		}
		
		// restore mode data to re[0]
		for(int j=0,cc=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int m=0;m<mc;m++) result.getModes()[m].getData()[0][0][j][i]=(float)Vmtrx.get(cc,m);
			
			cc++;
		}
		
		System.out.println(result);
		
		return result;
	}
	
	/**
     * singular spectrum analysis
     *
     * @param	v	a given Variable
     * @param	mc	mode count need to be restored
     * @param	ec	eigenvalue count need to be restored
     * @param	mm	length of matrix in spatial-dimension
     * 
     * @return	result	result of SSA
     */
	public static EOFResult SSA(Variable v,int mc,int ec,int mm){
		EOFResult result=new EOFResult(mc,ec,v);
		
		int t=v.getTCount(),y=v.getYCount(),x=v.getXCount();
		int nn=t+1-mm;
		
		if( nn<=t/2f ) throw new IllegalArgumentException("invalid length of matrix in spatial-dimension");
		if( mm>= nn  ) throw new IllegalArgumentException("spatial grid is more than temporal grid");
		if(x!=1||y!=1) throw new IllegalArgumentException("only time series allowed");
		
		float undef=v.getUndef();	float[][][][] data=v.getData();
		
		float[][] mtrix=new float[mm][nn];	// 2D matrix
		double[]   sqsv =new double[mc];	// sqrt of the singular values
		double[]   sv   =null;				// singular values
		
		DMatrixRMaj Vmtrx=null;			// modes
		
		int cc=0;
		// change 4D array to 2D array
		for(int i=0;i<mm;i++) System.arraycopy(data[0][0][0],i,mtrix[i],0,nn);
		
		checkUndef(mtrix,undef);
		
		// calculate the Matrix of covariance
		double[][] Sdata=new double[mm][mm];
		
		cCovMatrix1(mtrix,Sdata);
		
		SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(mm,mm,true,true,true);
		if(!svdF.decompose(new DMatrixRMaj(Sdata))) throw new IllegalArgumentException("matrix cannot be decomposed");
		Vmtrx=svdF.getU(null,false);
		sv=svdF.getSingularValues();
		
		restoreContribution(cRatio(sv,ec),sv,result);
		
		for(int m=0;m<mc;m++){
			sqsv[m]=Math.sqrt(sv[m]/nn);
			
			for(int i=0;i<mm;i++) Vmtrx.set(i,m,Vmtrx.get(i,m)*sqsv[m]);
		}
		
		// restore temporal data
		for(int m=0;m<mc;m++){
			float[][][] r1data=result.getTimes()[m].getData()[0];
			
			// T=V.transpose().multiply(A);	restore temporal data to re[1]
			for(int l=0;l<nn;l++){
				for(int g=0;g<mm;g++) r1data[0][0][l]+=Vmtrx.get(g,m)*mtrix[g][l];
				
				r1data[0][0][l]/=sqsv[m]*sqsv[m];
			}
		}
		
		// restore mode data to re[0]
		cc=0;
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[0][j][i][0]!=undef){
			for(int m=0;m<mc;m++) result.getModes()[m].getData()[0][0][j][i]=(float)Vmtrx.get(cc,m);
			
			cc++;
		}
		
		System.out.println(result);
		
		return result;
	}
	
	/**
     * rotated empirical orthogonal decomposition
     *
     * @param	v	a given Variable
     * @param	mc	mode count need to be restored
     * @param	ec	eigenvalue count need to be restored
     * 
     * @return	result	result of REOF
     */
	public static EOFResult REOF(Variable v,int mc,int ec){
		EOFResult result=new EOFResult(mc,ec,v);
		
		int t=v.getTCount(),y=v.getYCount(),x=v.getXCount();
		
		float undef=v.getUndef();	float[][][] data=v.getData()[0];
		
		int ge=cDefinedGridCount(v);
		
		checkData(v,mc,ge);
		
		double[]     sv=null;				// singular value
		float[][] mtrix=new float[ge][];	// 2D matrix
		DMatrixRMaj Vmtrx=null;			// modes
		
		int cc=0;
		// change 4D array to 2D array
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(data[j][i][0]!=undef) mtrix[cc++]=data[j][i];
		
		checkUndef(mtrix,undef);
		
		// if m > t, then change the time and space
		if(t>=ge){
			// calculate the Matrix of covariance
			double[][] Sdata=new double[ge][ge];
			
			cCovMatrix1(mtrix,Sdata);
			
			SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(ge,ge,true,true,true);
			if(!svdF.decompose(new DMatrixRMaj(Sdata))) throw new IllegalArgumentException("matrix cannot be decomposed");
			Vmtrx=svdF.getU(null,false);
			sv=svdF.getSingularValues();
			
			// initial the factor loading matrix
			for(int m=0;m<mc;m++)
			for(int g=0;g<ge;g++) Vmtrx.set(g,m,Vmtrx.get(g,m)*Math.sqrt(sv[m]/t));
			
		}else{
			// calculate the Matrix of covariance
			double[][] Sdata=new double[t][t];
			
			cCovMatrix2(mtrix,Sdata);
			
			SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(ge,ge,true,true,true);
			if(!svdF.decompose(new DMatrixRMaj(Sdata))) throw new IllegalArgumentException("matrix cannot be decomposed");
			sv=svdF.getSingularValues();
			
			// V=A.multiply(svd.getU());
			DMatrixRMaj U=svdF.getU(null,false);
			
			Vmtrx=new DMatrixRMaj(ge,mc);
			
			for(int g=0;g<ge;g++)
			for(int m=0;m<mc;m++){
				for(int l=0;l<t;l++) Vmtrx.set(g,m,Vmtrx.get(g,m)+mtrix[g][l]*U.get(l,m));
				Vmtrx.set(g,m,Vmtrx.get(g,m)/Math.sqrt(t));
			}
		}
		
		System.out.println("\nStart rotating...");
		
		// calculate the common degree
		double[] h2=new double[ge];
		
		for(int g=0;g<ge;g++)
		for(int m=0;m<mc;m++) h2[g]+=Vmtrx.get(g,m)*Vmtrx.get(g,m);
		
		// calculate the sum of the variance contribution by per common degree
		double g1=0,g2=0,s2o=0;
		
		// start varimax rotating
		double ds2=0,s2=0,ccc=0,ddd=0,hhh=0,eee=0,theta;	int count=0;	final int mcount=30;
		double[] wl=new double[ge];	double[] wk=new double[ge];
		double[] ur=new double[ge];	double[] vr=new double[ge];
		
		do{
			// one cycle
			for(int ll=0;ll<mc;ll++)
			for(int kk=ll+1;kk<mc;kk++){
				// calculate theta
				for(int i=0;i<ge;i++){
					ur[i]=(Vmtrx.get(i,ll)*Vmtrx.get(i,ll)-Vmtrx.get(i,kk)*Vmtrx.get(i,kk))/h2[i];
					vr[i]=2*Vmtrx.get(i,ll)*Vmtrx.get(i,kk)/h2[i];
					
					ccc+=ur[i]*ur[i]-vr[i]*vr[i];
					ddd+=ur[i]*vr[i];
					hhh+=ur[i];
					eee+=vr[i];
				}
				
				theta=atan2(2*ddd-2*hhh*eee/ge,ccc-(hhh*hhh-eee*eee)/ge)/4;
				
				ccc=0;	ddd=0;	hhh=0;	eee=0;
				
				double sinth=sin(theta);
				double costh=cos(theta);
				
				// calculate new Vdata with theta
				for(int i=0;i<ge;i++){
					wl[i]=Vmtrx.get(i,ll);
					wk[i]=Vmtrx.get(i,kk);
					
					Vmtrx.set(i,ll, wl[i]*costh+wk[i]*sinth);
					Vmtrx.set(i,kk,-wl[i]*sinth+wk[i]*costh);
				}
			}
			
			// re-calculate the sum of the variance contribution by per common degree
			for(int m=0;m<mc;m++){
				g1=0;	g2=0;
				
				for(int i=0;i<ge;i++) g1+=pow(Vmtrx.get(i,m),4)/(h2[i]*h2[i]);
				for(int i=0;i<ge;i++) g2+=Vmtrx.get(i,m)*Vmtrx.get(i,m)/h2[i];
				
				g1/=ge;	g2*=g2/(ge*ge);
				
				s2+=g1-g2;
			}
			
			if(++count>=mcount) break;
			
			ds2=s2-s2o;	s2o=s2;
			
		}while(ds2>0.1f);
		
		System.out.print("rotated times: "+count);
		if(count==mcount) System.out.println(" (max loop count)");
		else System.out.println();
		
		// calculate the explained variance and output the contribution
		double[] st=new double[mc];
		
		for(int m=0;m<mc;m++)
		for(int g=0;g<ge;g++) st[m]+=Vmtrx.get(g,m)*Vmtrx.get(g,m);
		
		restoreContribution(cRatio(sv,ec),st,result);
		
		// storing temporal data
		for(int m=0;m<mc;m++){
			float[][][][] r1data=result.getTimes()[m].getData();
			
			// T=V.transpose().multiply(A);	restore temporal data to re[1]
			for(int l=0;l<t;l++){
				for(int g=0;g<ge;g++)
				r1data[0][0][0][l]+=Vmtrx.get(g,m)*mtrix[g][l];
				
				r1data[0][0][0][l]/=sv[m]/t;
			}
		}
		
		// restore mode data to re[0]
		cc=0;
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int m=0;m<mc;m++) result.getModes()[m].getData()[0][0][j][i]=(float)Vmtrx.get(cc,m);
			
			cc++;
		}
		
		System.out.println(result);
		
		return result;
	}
	
	/**
     * extended empirical orthogonal decomposition
     *
     * @param	v	a given Variable
     * @param	mc	mode count need to be restored
     * @param	ec	eigenvalue count need to be restored
     * @param	lag	lag time
     * 
     * @return	result	result of EEOF
     */
	public static EOFResult[] EEOF(Variable v,int mc,int ec,int lag){
		String vname=v.getName();
		EOFResult[] results=new EOFResult[3];
		
		for(int i=0;i<3;i++){
			v.setName(vname+i);
			results[i]=new EOFResult(mc,ec,v);
		}
		
		v.setName(vname);
		
		int t=v.getTCount(),y=v.getYCount(),x=v.getXCount();
		int cc=0,ttt=t-2*lag;
		int ge=cDefinedGridCount(v),ge3=ge*3;
		
		checkData(v,mc,ge);
		
		float undef=v.getUndef();
		
		if(mc>Math.min(ge,t))
		throw new IllegalArgumentException("output modes are more than the matrix could produce");
		
		float[][]  mtrix=new float[ge3][ttt];
		float[][][] data=v.getData()[0];
		
		// change 4D array to 2D array
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int l=0;l<ttt;l++) mtrix[cc][l]=data[j][i][l];
			cc++;
		}
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int l=lag;l<ttt+lag;l++) mtrix[cc][l-lag]=data[j][i][l];
			cc++;
		}
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int l=lag*2;l<ttt+2*lag;l++) mtrix[cc][l-lag*2]=data[j][i][l];
			cc++;
		}
		
		// check undef
		for(int i=0;i<ge3;i++)
		for(int j=0;j<ttt;j++)
		if(mtrix[i][j]==undef) throw new IllegalArgumentException("Matrix contains undef values");
		
		double[] sv=null;
		DMatrixRMaj Vmtrx=null;			// modes
		
		// if m > t, then change the time and space
		if(ttt>=ge3){
			// calculate the Matrix of covariance
			double[][] Sdata=new double[ge3][ge3];
			
			System.out.println("start calculating the covariance matrix");
			for(int i=0;i<ge3;i++)
			for(int j=i;j<ge3;j++)
			for(int k=0;k<ttt;k++) Sdata[i][j]+=mtrix[i][k]*mtrix[j][k];
			
			for(int i=1;i<ge3;i++)
			for(int j=0;j<i;j++) Sdata[i][j]=Sdata[j][i];
			System.out.println("finish calculating the covariance matrix");
			
			SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(ge3,ge3,true,true,true);
			if(!svdF.decompose(new DMatrixRMaj(Sdata))) throw new IllegalArgumentException("matrix cannot be decomposed");
			Vmtrx=svdF.getU(null,false);
			sv=svdF.getSingularValues();
			
			for(int i=0;i<3;i++)
			restoreContribution(cRatio(sv,ec),sv,results[i]);
			
			for(int k=0;k<mc;k++)
			if(sv[k]!=0)
			for(int i=0;i<ge;i++) Vmtrx.set(i,k,Vmtrx.get(i,k)*Math.sqrt(sv[k]/t));
			
		}else{
			// calculate the Matrix of covariance
			double[][] Sdata=new double[ttt][ttt];
			
			System.out.println("start calculating the covariance matrix");
			for(int i=0;i<ttt;i++)
			for(int j=i;j<ttt;j++)
			for(int k=0;k<ge3;k++) Sdata[i][j]+=mtrix[k][i]*mtrix[k][j];
			
			for(int i=1;i<ttt;i++)
			for(int j=0;j<i;j++) Sdata[i][j]=Sdata[j][i];
			System.out.println("finish calculating the covariance matrix");
			
			SingularValueDecomposition_F64<DMatrixRMaj> svdF=DecompositionFactory_DDRM.svd(ttt,ttt,true,true,true);
			if(!svdF.decompose(new DMatrixRMaj(Sdata))) throw new IllegalArgumentException("matrix cannot be decomposed");
			sv=svdF.getSingularValues();
			
			for(int i=0;i<3;i++)
			restoreContribution(cRatio(sv,ec),sv,results[i]);
			
			// V=A.multiply(svd.getU());
			DMatrixRMaj U=svdF.getU(null,false);
			Vmtrx=new DMatrixRMaj(ge,mc);
			
			for(int i=0;i<ge3;i++)
			for(int j=0;j<mc;j++)
			for(int k=0;k<ttt;k++){
				Vmtrx.set(i,j,Vmtrx.get(i,j)+mtrix[i][k]*U.get(k,j));
			}
			
			for(int k=0;k<mc;k++)
			if(sv[k]!=0)
			for(int i=0;i<ge;i++) Vmtrx.set(i,k,Vmtrx.get(i,k)/Math.sqrt(t));
		}
		
		// T=V.transpose().multiply(A);	restore temporal data to re[1]
		for(int i=0;i<mc;i++){
			for(int j=0;j<ttt;j++){
				for(int k=0;k<ge;k++)
				results[0].getTimes()[i].getData()[0][0][0][j]+=Vmtrx.get(k,i)*mtrix[k][j];
				
				results[0].getTimes()[i].getData()[0][0][0][j]/=sv[i]/t;
			}
			
			for(int j=0;j<ttt;j++){
				for(int k=ge+1;k<2*ge;k++)
				results[1].getTimes()[i].getData()[0][0][0][j]+=Vmtrx.get(k,i)*mtrix[k][j];
				
				results[1].getTimes()[i].getData()[0][0][0][j]/=sv[i]/t;
			}
			
			for(int j=0;j<ttt;j++){
				for(int k=2*ge;k<ge3;k++)
				results[2].getTimes()[i].getData()[0][0][0][j]+=Vmtrx.get(k,i)*mtrix[k][j];
				
				results[2].getTimes()[i].getData()[0][0][0][j]/=sv[i]/t;
			}
		}
		
		// restore mode data to re[0]
		cc=0;
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int k=0;k<mc;k++) results[0].getModes()[k].getData()[0][0][j][i]=(float)Vmtrx.get(cc,k);
			
			cc++;
		}

		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int k=0;k<mc;k++) results[1].getModes()[k].getData()[0][0][j][i]=(float)Vmtrx.get(cc,k);
			
			cc++;
		}

		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int k=0;k<mc;k++) results[2].getModes()[k].getData()[0][0][j][i]=(float)Vmtrx.get(cc,k);
			
			cc++;
		}
		
		System.out.println(results[0]);
		
		return results;
	}
	
	/**
     * complex empirical orthogonal decomposition
     *
     * @param	v	a given Variable
     * @param	mc	mode count need to be restored
     * @param	ec	eigenvalue count need to be restored
     * 
     * @return	result	result of CEOF
     
	public static EOFResult[] CEOF(Variable v,int mc,int ec){
		String vname=v.getName();
		EOFResult[] results=new EOFResult[2];
		
		v.setName(vname+"am");	results[0]=new EOFResult(mc,ec,v);
		v.setName(vname+"ph");	results[1]=new EOFResult(mc,ec,v);
		
		v.setName(vname);
		
		int t=v.getTCount(),y=v.getYCount(),x=v.getXCount();
		int cc=0;
		int ge=cNonUndefGridCount(v);
		
		checkData(v,mc,ge);
		
		float undef=v.getUndef();
		
		float[][][] data=v.getData()[0];
		
		if(mc>Math.min(ge,t))
		throw new IllegalArgumentException("output modes are more than the matrix could produce");
		
		float[] real=null;			Complex[][] Vdata=null;
		float[] imag=new float[t];	Complex[][] mtrix=new Complex[ge][t];
		
		// change 4D array to 2D array
		for(int i=0;i<x;i++)
		for(int j=0;j<y;j++)
		if(data[j][i][0]!=undef){
			real=data[j][i];
			
			cOrthogonalData(real,imag);
			
			for(int l=0;l<t;l++)
			mtrix[cc][l]=new Complex(data[j][i][l],imag[l]);
			
			cc++;
		}
		
		for(int j=0;j<t;j++)
		for(int i=0;i<ge;i++)
		if(mtrix[i][j].getReal()==undef) throw new IllegalArgumentException("Matrix contains undef values");
		
		// if m > t, then change the time and space
		if(t>=ge){
			// calculate the Matrix of covariance
			Complex[][] Sdata=new Complex[ge][ge];
			
			for(int i=0;i<ge;i++)
			for(int j=i;j<ge;j++){
				Sdata[i][j]=new Complex(0,0);
				
				for(int k=0;k<t;k++) Sdata[i][j].plusEq(mtrix[i][k].conjugate().multiplyEq(mtrix[j][k]));
			}
			
			for(int i=1;i<ge;i++)
			for(int j=0;j<i;j++) Sdata[i][j]=Sdata[j][i].conjugate();
			
			ComplexSingularValueDecomposition csvd=new ComplexSingularValueDecomposition(Sdata);
			
			Vdata=csvd.getU();
			
			float[] sv=csvd.getSingularValues();
			
			for(int i=0;i<2;i++)
			restoreContribution(cRatio(sv,ec),sv,results[i]);
			
		}else{
			// calculate the Matrix of covariance
			Complex[][] Sdata=new Complex[t][t];
			
			for(int i=0;i<t;i++)
			for(int j=i;j<t;j++){
				Sdata[i][j]=new Complex(0,0);
				
				for(int k=0;k<ge;k++) Sdata[i][j].plusEq(mtrix[k][i].conjugate().multiplyEq(mtrix[k][j]));
			}
			
			for(int i=1;i<t;i++)
			for(int j=0;j<i;j++) Sdata[i][j]=Sdata[j][i].conjugate();
			
			ComplexSingularValueDecomposition csvd=new ComplexSingularValueDecomposition(Sdata);
			
			float[] sv=csvd.getSingularValues();
			
			for(int i=0;i<2;i++)
			restoreContribution(cRatio(sv,ec),sv,results[i]);
			
			// V=A.multiply(csvd.getU());
			Complex[][] Udata=csvd.getU();
			
			Vdata=new Complex[ge][mc];
			
			for(int i=0;i<ge;i++)
			for(int j=0;j<mc;j++){
				Vdata[i][j]=new Complex(0,0);
				
				for(int k=0;k<t;k++) Vdata[i][j].plusEq(mtrix[i][k].multiply(Udata[k][j]));
			}
			
			for(int k=0;k<mc;k++)
			if(sv[k]!=0){
				float tmp=(float)(Math.sqrt(sv[k]));
				
				for(int i=0;i<ge;i++) Vdata[i][k].divideEq(tmp);
			}
		}
		
		// T=V.transpose().multiply(A);	restore temporal data to re[1]
		Complex tmp=new Complex(0,0);
		for(int i=0;i<mc;i++)
		for(int j=0;j<t;j++){
			tmp.setRI(0,0);
			for(int k=0;k<ge;k++)
			tmp.plusEq(Vdata[k][i].multiply(mtrix[k][j]));
			
			results[0].getTimes()[i].getData()[0][0][0][j]=tmp.conjugate().multiplyEq(tmp).sqrtEq().getReal();
			results[1].getTimes()[i].getData()[0][0][0][j]=tmp.getArg();
		}
		
		// restore mode data to re[0]
		cc=0;
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(data[j][i][0]!=undef){
			for(int k=0;k<mc;k++){
				results[0].getModes()[k].getData()[0][0][j][i]=Vdata[cc][k].conjugate().multiplyEq(Vdata[cc][k]).sqrtEq().getReal();
				results[1].getModes()[k].getData()[0][0][j][i]=Vdata[cc][k].getArg();
			}
			
			cc++;
		}
		
		System.out.println(results[0]);
		
		return results;
	}*/
	
	
	/**
     * reconstruct the original variable
     *
     * @param	mod		mode of EOF result
     * @param	tim		time series of EOF result
     *
     * @return	rec		reconstructed variable
     */
	public static Variable reconstruct(Variable mod,Variable tim){
		int t=mod.getTCount(),z=mod.getZCount(),y=mod.getYCount(),x=mod.getXCount();
		
		if(t!=1)
			throw new IllegalArgumentException("t-count of mod should only be 1");
		if(tim.getZCount()!=z)
			throw new IllegalArgumentException("z-counts of mod and tim are not the same");
		if(tim.getYCount()!=1&&tim.getXCount()!=1)
			throw new IllegalArgumentException("tim should contain only one grid");
		if(mod.isTFirst()||tim.isTFirst())
			throw new IllegalArgumentException("T-First error");
		
		t=tim.getTCount();
		
		Range recrange=new Range(t,z,y,x);
		
		Variable rec=new Variable(mod.getName(),recrange);
		rec.setCommentAndUnit("EOF reconstruction of "+mod.getName());
		
		float[][][][] recdata=rec.getData();
		float[][][][] moddata=mod.getData();
		float[][][][] timdata=tim.getData();
		
		for(int l=0;l<t;l++)
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) recdata[k][j][i][l]+=moddata[0][k][j][i]*timdata[k][0][0][l];
		
		Range modrange=mod.getRange();	Range timrange=tim.getRange();
		
		recrange.setYRange(modrange);	recrange.setTRange(timrange);
		recrange.setXRange(modrange);	recrange.setZRange(timrange);
		
		return rec;
	}
	
	/**
     * vertical decompose a given variable into normal modes
     *
     * @param	N2		Brunt-Vaisala frequency (could be static stability)
     * @param	v		variable to be decomposed
     * @param	mc		mode count need to be restored
     * @param	dz		delta z (could be pressure)
     *
     * @return	re		modes for decomposition of v
     */
	public static Variable[] verticalDecompose(Variable N2,Variable v,int mc,float dz){
		if(!N2.isLike(v)) throw new IllegalArgumentException("dimensions not same");
		
		int t=N2.getTCount(),z=N2.getZCount(),y=N2.getYCount(),x=N2.getXCount();
		
		VerticalModeDecomposition vmd=new VerticalModeDecomposition(z);
		
		Variable[] re=new Variable[mc];
		
		for(int i=0;i<mc;i++){
			re[i]=new Variable(v.getName()+"m"+i,v);
			re[i].setCommentAndUnit(i+" vertical mode of "+v.getName());
		}
		
		float[][][][] Ndata=N2.getData();
		float[][][][] vdata=v.getData();
		
		if(N2.isTFirst()){
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] Npf=new float[z-1];
				float[] vpf=new float[z];
				
				for(int k=0,K=z-1;k<K;k++){
					Npf[k]=Ndata[l][k][j][i];
					vpf[k]=vdata[l][k][j][i];
				}
				
				vpf[z-1]=vdata[l][z-1][j][i];
				
				vmd.updateNSquare(Npf,dz);
				
				float[][] modes=vmd.decompose(vpf,mc);
				
				for(int m=0;m<mc;m++){
					float[][][] rdata=re[m].getData()[l];
					
					for(int k=0;k<z;k++) rdata[k][j][i]=modes[m][k];
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] Npf=new float[z];
				float[] vpf=new float[z];
				
				for(int k=0,K=z-1;k<K;k++){
					Npf[k]=Ndata[k][j][i][l];
					vpf[k]=vdata[k][j][i][l];
				}
				
				vpf[z-1]=vdata[z-1][j][i][l];
				
				vmd.updateNSquare(Npf,dz);
				
				float[][] modes=vmd.decompose(vpf,mc);
				
				for(int m=0;m<mc;m++){
					float[][][][] rdata=re[m].getData();
					
					for(int k=0;k<z;k++) rdata[k][j][i][l]=modes[m][k];
				}
			}
		}
		
		return re;
	}
	
	
	/*** private method ***/
	private static int cDefinedGridCount(Variable v){
		int y=v.getYCount();	int ucount =0;
		int x=v.getXCount();	float undef=v.getUndef();
		
		float[][][] data=v.getData()[0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(data[j][i][0]==undef) ucount++;
		
		return x*y-ucount;
	}
	
	private static double cRatio(double[] sv,int ec){
		return cAccumulatedSV(sv,ec)/cAccumulatedSV(sv,sv.length);
	}
	
	private static double cAccumulatedSV(double[] sv,int length){
		float sum=0;
		
		for(int i=0;i<length;i++) if(sv[i]!=0) sum+=sv[i];
		
		return sum;
	}
	
	/*
	private static void cOrthogonalData(float[] real,float[] imag){
		Complex[] co=FastFourier.fft(real);
		
		int length=real.length;
		float o,oi,mean=StatisticsUtil.cArithmeticMean(real);
		
		for(int k=0;k<length/2;k++){
			o=(float)(2*PI*k/length);
			
			for(int i=0;i<length;i++){
				oi=o*i;
				
				//imag[i]+=co[k].getReal()*(float)cos(oi)+co[k].getImag()*(float)sin(oi);
				imag[i]+=co[k].getReal()*(float)sin(oi)-co[k].getImag()*(float)cos(oi);
			}
		}
		
		//for(int i=0;i<length;i++) imag[i]=imag[i]/length*2-mean;
		for(int i=0;i<length;i++) imag[i]=imag[i]/length*2+mean;
	}*/
	
	private static void checkData(Variable v,int mc,int ge){
		if(v.getZCount()!=1)
			throw new IllegalArgumentException("z-count should be 1");
		if(mc<=0)
			throw new IllegalArgumentException("mode count should be larger than 0");
		if(v.isTFirst())
			throw new IllegalArgumentException("T should not be the first dimension");
		if(v.getTCount()<3)
			throw new IllegalArgumentException("t-length is too small (<2)");
		if(mc>Math.min(ge,v.getTCount()))
			throw new IllegalArgumentException("output modes are more than the matrix could produce");
	}
	
	private static void checkUndef(float[][] data,float undef){
		for(int j=0,J=data.length;j<J;j++)
		for(int i=0,I=data[0].length;i<I;i++)
		if(data[j][i]==undef) throw new IllegalArgumentException("Matrix contains undef values");
	}
	
	private static void cCovMatrix1(float[][] mtrix,double[][] Sdata){
		int ge=mtrix.length;	int t=mtrix[0].length;
		
		// t>=ge
		for(int g=0;g<ge;g++)
		for(int h=g;h<ge;h++)
		for(int k=0;k<t;k++) Sdata[g][h]+=mtrix[g][k]*mtrix[h][k];
		
		for(int g=1;g<ge;g++)
		for(int h=0;h<g;h++) Sdata[g][h]=Sdata[h][g];
	}
	
	private static void cCovMatrix2(float[][] mtrix,double[][] Sdata){
		int ge=mtrix.length;	int t=mtrix[0].length;
		
		// t<=ge
		for(int l=0;l<t;l++)
		for(int m=l;m<t;m++)
		for(int g=0;g<ge;g++) Sdata[l][m]+=mtrix[g][l]*mtrix[g][m];
		
		for(int l=1;l<t;l++)
		for(int m=0;m<l;m++) Sdata[l][m]=Sdata[m][l];
	}
	
	private static void restoreContribution(double ratio,double[] sv,EOFResult result){
		int tcount=result.getContributions()[0].getTCount();
		
		double svsum=cAccumulatedSV(sv,tcount);
		float coeff=result.getErrorCoefficient();
		
		if(tcount>sv.length)
		throw new IllegalArgumentException("Not enough singular values for restore");
		
		result.setRatio((float)ratio);
		
		for(int m=0;m<2;m++){
			float[] sd=result.getContributions()[m].getData()[0][0][0];
			
			for(int i=0;i<tcount;i++)
			sd[i]=(float)(sv[i]/svsum*ratio);
			
			if(m==1)
			for(int i=0,I=sd.length;i<I;i++) sd[i]*=coeff;
		}
	}
	
	
	/** test
	public static void main(String arg[]){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/sst.mnmean.nc");
		DataDescriptor dd=df.getDataDescriptor();
		
		Range range=new Range("t(1297,1896);lon(120,260);lat(-30,30)",dd);
		
		Variable v=df.getVariables(range,"sst")[0];
		
		FilterMethods.cycleFilter(v,12);
		
		EOFResult re=EOF(v,20,20);
		
		miniufo.io.CtlDataWriteStream cdws=null;
		cdws=new miniufo.io.CtlDataWriteStream("D:/mod.dat");
		cdws.writeData(dd,re.getModes());	cdws.closeFile();
		
		cdws=new miniufo.io.CtlDataWriteStream("D:/tim.dat");
		cdws.writeData(dd,re.getTimes());	cdws.closeFile();
		
		cdws=new miniufo.io.CtlDataWriteStream("D:/contr.dat");
		cdws.writeData(dd,re.getContributions());	cdws.closeFile();
	}*/
}
