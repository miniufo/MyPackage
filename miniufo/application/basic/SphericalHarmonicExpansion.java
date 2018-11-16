/**
 * @(#)SphericalHarmonicExpansion.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.application.GeoFluidApplication;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.mathsphysics.AssociatedLegendre;
import miniufo.mathsphysics.FastFourier;
import static miniufo.diagnosis.SpatialModel.REarth;


/**
 * spherical harmonic expansion
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class SphericalHarmonicExpansion extends GeoFluidApplication{
	//
	private boolean setM=false;
	
	private int M=0;
	
	private float dlat=0;
	
	private float[] lcos=null;
	private float[] lsin=null;
	
	private float[][][] plm=null;
	
	private FastFourier ff=null;
	
	
	/**
     * constructor
     *
     * @param	sm	spacial model used for the calculation
     */
	public SphericalHarmonicExpansion(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isGlobal()) throw new IllegalArgumentException("Not a global spatial model");
		
		y=ssm.getYCount();	lcos=ssm.getLCos();
		x=ssm.getXCount();	lsin=ssm.getLSin();	M=x/2;
		
		plm=new float[y][][];
		
		dlat=ssm.getYDef().getIncrements()[0];
		
		ff=new FastFourier(x);
	}
	
	
	/*** getor and setor ***/
	public int getM(){ return M;}
	
	public void setM(int M){
		if(M<2)
		throw new IllegalArgumentException("M should larger than 1");
		
		if(M>x/2)
		throw new IllegalArgumentException("M should not larger than x/2");
		
		this.M=M;
		
		for(int j=0;j<y;j++){
			plm[j]=AssociatedLegendre.legendreT(M,lsin[j]);
			
			for(int m=0;m<plm[j].length;m++)
			for(int n=0;n<plm[j][m].length;n++)
			plm[j][m][n]/=(float)Math.sqrt(2);	// normalization
		}
		
		setM=true;
	}
	
	
	/**
	 * transform a variable into spectral space
	 * 
	 * @param v	a give variable
	 * 
	 * @return	the spectrum coefficient,
	 * 			[0] is the real part and [1] is the image part
	 */
	public Variable[] cSpectrumCoefficient(Variable v){
		if(y!=v.getYCount()||x!=v.getXCount())
		throw new IllegalArgumentException("not a global variable or dimensions not the same");
		
		checkM();
		
		t=v.getTCount();	z=v.getZCount();
		
		Range nr=new Range(t,z,M+1,M+1);	Range r=v.getRange();
		
		Variable[] co=new Variable[2];	float undef=v.getUndef();
		co[0]=new Variable("shre"+v.getName(),v.isTFirst(),nr);	co[0].setUndef(undef);
		co[1]=new Variable("shim"+v.getName(),v.isTFirst(),nr);	co[1].setUndef(undef);
		
		co[0].setCommentAndUnit("real part of spherical harmonics expansion coefficiences");
		co[1].setCommentAndUnit("image part ofspherical harmonics expansion coefficiences");
		
		float[][][][] rdata=co[0].getData();
		float[][][][] idata=co[1].getData();
		float[][][][] vdata=    v.getData();
		long tt=System.nanoTime();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] fftRe=new float[y][];
				float[][] fftIm=new float[y][];
				
				// FFT along jth latitude, start from South Pole.
				// transformation from grid space to spectrum space
				for(int j=0;j<y;j++){
					ff.fftMixedRadix(vdata[l][k][j]);
					fftRe[j]=ff.getResultRealPartCopy();
					fftIm[j]=ff.getResultImagePartCopy();
				}
				
				for(int m=0;m<=M;m++){
				for(int n=m;n<=M;n++){
					float sumRe=0,sumIm=0;
					
					for(int j=0;j<y;j++){
						float tmp=plm[j][n][m]*dlat*lcos[j];
						sumRe+=tmp*fftRe[j][m];
						sumIm+=tmp*fftIm[j][m];
					}
					
					rdata[l][k][n][m]=sumRe;
					idata[l][k][n][m]=sumIm;
				}}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] fftRe=new float[y][];
				float[][] fftIm=new float[y][];
				
				// FFT along jth latitude, start from South Pole.
				// transformation from grid space to spectrum space
				for(int j=0;j<y;j++){
					float[]  fftBuf=new float[x];
					for(int i=0;i<x;i++) fftBuf[i]=vdata[k][j][i][l];
					ff.fftMixedRadix(fftBuf);
					fftRe[j]=ff.getResultRealPartCopy();
					fftIm[j]=ff.getResultImagePartCopy();
				}
				
				for(int m=0;m<=M;m++){
				for(int n=m;n<=M;n++){
					float sumRe=0,sumIm=0;
					
					for(int j=0;j<y;j++){
						float tmp=plm[j][n][m]*dlat*lcos[j];
						sumRe+=tmp*fftRe[j][m];
						sumIm+=tmp*fftIm[j][m];
					}
					
					rdata[k][n][m][l]=sumRe;
					idata[k][n][m][l]=sumIm;
				}}
			}
		}
		System.out.println("cSpectrumCoefficient: "+(System.nanoTime()-tt)/1000000000.0+" sec");
		
		nr.setTRange(r);	nr.setZRange(r);
		
		nr.getYRange()[0]=r.getYRange()[0];	nr.getYRange()[1]=r.getYRange()[0]+M-1;
		nr.getXRange()[0]=r.getXRange()[0];	nr.getXRange()[1]=r.getXRange()[0]+2*M;
		
		return co;
	}
	
	/**
	 * reconstruct a variable from its spectral coefficient
	 * 
	 * @param re	real part of the coefficient
	 * @param im	image part of the coefficient
	 * 
	 * @return	the reconstructed variable (triangle truncated by M)
	 */
	public Variable reconstruct(Variable re,Variable im){
		if(!re.isLike(im)) throw new IllegalArgumentException("dimensions not same");
		
		checkM();
		
		t=re.getTCount();	z=re.getZCount();
		
		Range nr=new Range(t,z,y,x);	Range r=re.getRange();
		
		Variable rc=new Variable(re.getName().replace("shre","rec"),re.isTFirst(),nr);
		rc.setUndef(re.getUndef());
		rc.setCommentAndUnit("reconstructed "+re.getName().replace("shre",""));
		
		float[][][][] rdata=re.getData();
		float[][][][] idata=im.getData();
		float[][][][] cdata=rc.getData();
		long tt=System.nanoTime();
		
		if(re.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] fftRe=new float[y][M+1];
				float[][] fftIm=new float[y][M+1];
				
				for(int m=0;m<=M;m++){
				for(int n=m;n<=M;n++){
					float rr=rdata[l][k][n][m];
					float ii=idata[l][k][n][m];
					
					for(int j=0;j<y;j++){
						float tmp=plm[j][n][m];
						
						fftRe[j][m]+=tmp*rr;
						fftIm[j][m]+=tmp*ii;
					}
				}}
				
				// inverse FFT along each latitude
				// transform from spetrum space to grid space
				// only real part are needed (stored)
				for(int j=0;j<y;j++){
					float[] tmpRe=new float[x];
					float[] tmpIm=new float[x];
					// padding fftRe (fftIm) to tmpRe (tmpIm)
					// to match the length of x
					tmpRe[0]=fftRe[j][0];
					for(int i=1;i<=M;i++){
						tmpRe[i]=fftRe[j][i];
						tmpIm[i]=fftIm[j][i];
						tmpRe[x-i]= fftRe[j][i];
						tmpIm[x-i]=-fftIm[j][i];
					}
					
					ff.ifftMixedRadix(tmpRe,tmpIm);
					
					System.arraycopy(ff.getResultRealPart(),0,cdata[l][k][j],0,x);
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] fftRe=new float[y][M+1];
				float[][] fftIm=new float[y][M+1];
				
				for(int m=0;m<=M;m++){
				for(int n=m;n<=M;n++){
					float rr=rdata[k][n][m][l];
					float ii=idata[k][n][m][l];
					
					for(int j=0;j<y;j++){
						float tmp=plm[j][n][m];
						
						fftRe[j][m]+=tmp*rr;
						fftIm[j][m]+=tmp*ii;
					}
				}}
				
				// inverse FFT along each latitude
				// transform from spetrum space to grid space
				// only real part are needed (stored)
				for(int j=0;j<y;j++){
					float[] tmpRe=new float[x];
					float[] tmpIm=new float[x];
					// padding fftRe (fftIm) to tmpRe (tmpIm)
					// to match the length of x
					tmpRe[0]=fftRe[j][0];
					for(int i=1;i<=M;i++){
						tmpRe[i]=fftRe[j][i];
						tmpIm[i]=fftIm[j][i];
						tmpRe[x-i]= fftRe[j][i];
						tmpIm[x-i]=-fftIm[j][i];
					}
					
					ff.ifftMixedRadix(tmpRe,tmpIm);
					float[] rec=ff.getResultRealPart();
					
					for(int i=0;i<x;i++)
					cdata[k][j][i][l]=rec[i];
				}
			}
		}
		System.out.println("reconstruct: "+(System.nanoTime()-tt)/1000000000.0+" sec");
		
		nr.setTRange(r);	nr.setZRange(r);
		
		nr.getYRange()[0]=r.getYRange()[0];	nr.getYRange()[1]=r.getYRange()[0]+y-1;
		nr.getXRange()[0]=r.getXRange()[0];	nr.getXRange()[1]=r.getXRange()[0]+x-1;
		
		return rc;
	}
	
	
	public Variable solvePoissonEquation(Variable v){
		Variable[] co=cSpectrumCoefficient(v);
		return reconstructL(co[0],co[1]);
	}
	
	
	public Variable reconstructL(Variable re,Variable im){
		if(!re.isLike(im))
		throw new IllegalArgumentException("dimensions not same");
		
		checkM();
		
		t=re.getTCount();	z=re.getZCount();
		
		Range nr=new Range(t,z,y,x);	Range r=re.getRange();
		
		Variable rc=new Variable(re.getName().replace("shre","rec"),re.isTFirst(),nr);
		rc.setUndef(re.getUndef());
		rc.setCommentAndUnit("reconstructed "+re.getName().replace("shre",""));
		
		float[][][][] rdata=re.getData();
		float[][][][] idata=im.getData();
		float[][][][] cdata=rc.getData();
		long tt=System.nanoTime();
		
		if(re.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] fftRe=new float[y][M+1];
				float[][] fftIm=new float[y][M+1];
				
				for(int m=0;m<=M;m++)
				for(int n=m;n<=M;n++)
				if(n!=0){	// skip n==0
					float co=-n*(n+1)/REarth/REarth;
					float rr=rdata[l][k][n][m];
					float ii=idata[l][k][n][m];
					
					for(int j=0;j<y;j++){
						float tmp=plm[j][n][m]/co;
						
						fftRe[j][m]+=tmp*rr;
						fftIm[j][m]+=tmp*ii;
					}
				}
				
				// inverse FFT along each latitude
				// transform from spetrum space to grid space
				// only real part are needed (stored)
				for(int j=0;j<y;j++){
					float[] tmpRe=new float[x];
					float[] tmpIm=new float[x];
					// padding fftRe (fftIm) to tmpRe (tmpIm)
					// with 0 to match the length of x
					tmpRe[0]=fftRe[j][0];
					for(int i=1;i<=M;i++){
						tmpRe[i]=fftRe[j][i];
						tmpIm[i]=fftIm[j][i];
						tmpRe[x-i]= fftRe[j][i];
						tmpIm[x-i]=-fftIm[j][i];
					}
					
					ff.ifftMixedRadix(tmpRe,tmpIm);
					
					System.arraycopy(ff.getResultRealPart(),0,cdata[l][k][j],0,x);
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				float[][] fftRe=new float[y][M+1];
				float[][] fftIm=new float[y][M+1];
				
				for(int m=0;m<=M;m++)
				for(int n=m;n<=M;n++)
				if(n!=0){	// skip n==0
					float co=-n*(n+1)/REarth/REarth;
					float rr=rdata[k][n][m][l];
					float ii=idata[k][n][m][l];
					
					for(int j=0;j<y;j++){
						float tmp=plm[j][n][m]/co;
						
						fftRe[j][m]+=tmp*rr;
						fftIm[j][m]+=tmp*ii;
					}
				}
				
				// inverse FFT along each latitude
				// transform from spetrum space to grid space
				// only real part are needed (stored)
				for(int j=0;j<y;j++){
					float[] tmpRe=new float[x];
					float[] tmpIm=new float[x];
					// padding fftRe (fftIm) to tmpRe (tmpIm)
					// with 0 to match the length of x
					tmpRe[0]=fftRe[j][0];
					for(int i=1;i<=M;i++){
						tmpRe[i]=fftRe[j][i];
						tmpIm[i]=fftIm[j][i];
						tmpRe[x-i]= fftRe[j][i];
						tmpIm[x-i]=-fftIm[j][i];
					}
					
					ff.ifftMixedRadix(tmpRe,tmpIm);
					float[] rec=ff.getResultRealPart();
					
					for(int i=0;i<x;i++)
					cdata[k][j][i][l]=rec[i];
				}
			}
		}
		System.out.println("reconstructL: "+(System.nanoTime()-tt)/1000000000.0+" sec");
		
		nr.setTRange(r);	nr.setZRange(r);
		
		nr.getYRange()[0]=r.getYRange()[0];	nr.getYRange()[1]=r.getYRange()[0]+y-1;
		nr.getXRange()[0]=r.getXRange()[0];	nr.getXRange()[1]=r.getXRange()[0]+x-1;
		
		return rc;
	}
	
	
	/**
	 * check whether the triangle truncation number M has been set
	 */
	private void checkM(){
		if(!setM) throw new IllegalArgumentException(
			"set triangle truncation number M first\n" +
			"use this.setM(int) before any computational calls"
		);
	}
	
	
	/** test
	public static void main(String[] args){
		miniufo.descriptor.DataDescriptor ctl=miniufo.diagnosis.DiagnosisFactory.getDataDescriptor(
			"d:/Data/DiagnosisVortex/Linfa/Linfa.ctl"
			//"d:/Spec/grid.ctl"
			//"d:/noPole.ctl"
		);
		
		SphericalSpacialModel ssm=new SphericalSpacialModel(ctl);
		
		Range r=new Range("t(1,1);lev(200,200)",ctl);
		
		Variable u=new Variable("t",false,r);
		
		SphericalHarmonicExpansion sh=new SphericalHarmonicExpansion(ssm);
		
		miniufo.io.DataRead dr=miniufo.io.DataIOFactory.getDataRead(ctl);
		dr.readData(u);	dr.closeFile();
		
		//miniufo.io.DataWrite dw2=null;
		//dw2=miniufo.io.DataIOFactory.getDataWrite(ctl,"d:/noPole.dat");
		//dw2.writeData(ctl,u);System.exit(0);
		
		sh.setM(180);
		
		//Variable[] re=sh.cSpectrumCoefficient(u);
		//Variable  urc=sh.reconstruct(re[0],re[1]);
		
		Variable[] re2=sh.cSpectrumCoefficient(u);
		Variable  urc2=sh.reconstruct(re2[0],re2[1]);urc2.setName("urec2");
		
		miniufo.io.DataWrite dw=null;
		dw=miniufo.io.DataIOFactory.getDataWrite(ctl,"d:/shF.dat");
		dw.writeData(ctl,u,urc2);
	}*/
}
