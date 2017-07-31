/**
 * @(#)VerticalModeDecomposition.java	1.0 2014.07.24
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import java.util.Arrays;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;


/**
 * Decomposite a vertical profile into vertical modes
 *
 * @version 1.0, 2014.07.24
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class VerticalModeDecomposition{
	//
	private boolean hasN2=false;
	
	private int len=0;
	
	private float[]   eigvalues =null;
	private float[][] eigvectors=null;
	
	private DMatrixRMaj A=null;
	
	private EigenDecomposition_F64<DMatrixRMaj> ed=null;
	
	
	/**
	 * constructor
	 * 
	 * @param	len		vertical data length
	 */
	public VerticalModeDecomposition(int len){
		if(len<3) throw new IllegalArgumentException("length should be at least 3");
		
		this.len=len;
		
		eigvalues =new float[len];
		eigvectors=new float[len][];
		
		A=new DMatrixRMaj(len,len);
		
		ed=DecompositionFactory_DDRM.eig(len,true);
	}
	
	
	/**
	 * Construct a matrix by setting Brunt-Vaisala frequency.
	 * Then find the eigenvalues and eigenvectors of the matrix,
	 * which is the difference approximation to the linear operator:
	 * L = -d(1/N^2*d( )/dz)/dz
	 * 
	 * Since central difference scheme is used, N2 is defined at half grids
	 * and thus its length should be one less than the matrix size
	 * 
	 * @param	N2	profile of Brunt-Vaisala frequency (could be static stability)
	 * @param	dz	delta z (could be pressure)
	 */
	public void updateNSquare(float[] N2,float dz){
		if(N2.length!=len-1) throw new IllegalArgumentException(
			"N2 is defined at half grids so its length ("+N2.length+")"+
			"should be one less than the matrix size ("+len+")"
		);
		
		A.set(0,0, 1.0/(dz*N2[0]*dz));	A.set(len-1,len-1, 1.0/((dz*N2[len-2]*dz)));
		A.set(0,1,-1.0/(dz*N2[0]*dz));	A.set(len-1,len-2,-1.0/((dz*N2[len-2]*dz)));
		
		for(int i=1,I=len-1;i<I;i++){
			A.set(i,i-1,-1.0/(dz*N2[i-1]*dz));
			A.set(i,i  , 1.0/(dz*N2[i-1]*dz)+1.0/(dz*N2[i]*dz));
			A.set(i,i+1,-1.0/(dz*N2[i]*dz));
		}
		
		decomposition();
		ascendingOrder();
		
		hasN2=true;
	}
	
	public float[][] decompose(float[] var,int modes){
		if(!hasN2) throw new IllegalArgumentException("call updateNSquare() first");
		if(var.length!=len) throw new IllegalArgumentException("lengths not equal");
		
		DMatrixRMaj data=new DMatrixRMaj(len,1);
		
		for(int i=0;i<len;i++) data.set(i,var[i]);
		
		float[] coeff=getCoefficients(data);
		
		float[][] re=new float[modes][len];
		
		for(int i=0;i<modes;i++)
		for(int j=0;j<len;j++) re[i][j]=coeff[i]*eigvectors[i][j];
		
		return re;
	}
	
	
	/*** helper methods ***/
	private void decomposition(){
		if(!ed.decompose(A))
			throw new IllegalArgumentException("cannot decompose the matrix");
		
		if(ed.getNumberOfEigenvalues()!=len)
			throw new IllegalArgumentException("number of eigenvalues do not equal matrix size");
		
		if(ed.inputModified())
			throw new IllegalArgumentException("input matrix has been modified");
	}
	
	private void ascendingOrder(){
		float[]   eigvlu=getAllEigenValues();
		float[][] eigvct=getAllEigenVectors();
		
		System.arraycopy(eigvlu,0,eigvalues,0,len);
		
		Arrays.sort(eigvalues);
		
		for(int i=0;i<len;i++){
			int idx=0;
			
			for(int j=0;j<len;j++) if(eigvalues[i]==eigvlu[j]){ idx=j; break;}
			
			
			eigvectors[i]=eigvct[idx];
			
			float tmp=eigvectors[i][0];
			
			for(int k=0;k<len;k++) eigvectors[i][k]/=tmp;
		}
	}
	
	private float[] getCoefficients(DMatrixRMaj data){
		DMatrixRMaj inv=new DMatrixRMaj(len,len);
		DMatrixRMaj coe=new DMatrixRMaj(len,1);
		
		for(int i=0;i<len;i++)
		for(int j=0;j<len;j++) inv.set(i,j,eigvectors[j][i]);
		
		if(!CommonOps_DDRM.invert(inv))
		throw new IllegalArgumentException("cannot invert the eigenvector matrix");
		
		CommonOps_DDRM.mult(inv,data,coe);
		
		return double2Float(coe.data);
	}
	
	private float[] getAllEigenValues(){
		float[] re=new float[len];
		
		for(int i=0;i<len;i++) re[i]=(float)ed.getEigenvalue(i).getReal();
		
		return re;
	}
	
	private float[][] getAllEigenVectors(){
		float[][] re=new float[len][];
		
		for(int i=0;i<len;i++) re[i]=double2Float(ed.getEigenVector(i).data);
		
		return re;
	}
	
	private float[] double2Float(double[] array){
		int len=array.length;
		
		float[] re=new float[len];
		
		for(int i=0;i<len;i++) re[i]=(float)array[i];
		
		return re;
	}
	
	
	/** test
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("/lustre/home/qianyk/Data/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		SphericalSpacialModel ssm=new SphericalSpacialModel(dd);
		ThermoMethodsInSC tm=new ThermoMethodsInSC(ssm);
		
		Variable T=df.getVariables(new Range("lon(120,120);lat(20,20);t(1,1)",dd),"t")[0];
		Variable S=tm.cStaticStabilityArgByT(T);
		
		float[] N2=new float[S.getZCount()-1];
		
		for(int i=0;i<N2.length;i++) N2[i]=S.getData()[0][i][0][0];
		
		VerticalModeDecomposition vmd=new VerticalModeDecomposition(S.getZCount());
		vmd.updateNSquare(N2,dd.getDZDef()[0]);
		
		float[] Tdata=new float[S.getZCount()];
		
		for(int k=0;k<S.getZCount();k++) Tdata[k]=T.getData()[0][k][0][0];
		
		float[][] vecs=vmd.decompose(Tdata,S.getZCount());
		
		Variable[] ms=new Variable[S.getZCount()];
		
		for(int i=0;i<ms.length;i++){
			ms[i]=new Variable("m"+(i+1),S);
			
			for(int k=0;k<S.getZCount();k++) ms[i].getData()[0][k][0][0]=vecs[i][k];
		}
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"/lustre/home/qianyk/Data/VMD.dat");
		dw.writeData(dd,ArrayUtil.concatAll(Variable.class,new Variable[]{T,S},ms));
		dw.closeFile();
	}*/
}
