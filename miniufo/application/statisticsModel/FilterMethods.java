/**
 * @(#)FilterMethods.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.mathsphysics.PolynomialFitter;
import miniufo.statistics.FilterModel;


/**
 * Basic filter methods for Variable
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class FilterMethods extends StatisticsApplication{
	
	/**
	 * prevent from instantiate
	 */
	private FilterMethods(){}
	
	
	/**
     * removing the cycle of a given length and returning the mean cycle
     *
     * @param	v			a given Variable
     * @param	cycle		length of cycle that need to be removed
     * 
     * @return	nv			mean cycle
     */
	public static Variable cycleFilter(Variable v,int cycle){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float undef=v.getUndef();
		
		Range nr=new Range(cycle,z,y,x);	Range r=v.getRange();
		
		Variable nv=new Variable(v.getName(),v.isTFirst(),nr);
		nv.setComment(v.getComment());
		nv.setUnit(v.getUnit());
		nv.setUndef(undef);
		nv.setValue(undef);
		
		float[][][][]  vdata= v.getData();
		float[][][][] nvdata=nv.getData();
		
		if(v.isTFirst()){
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	bk:		for(int i=0;i<x;i++){
				float[] buf=new float[t];
				
				for(int l=0;l<t;l++){
					buf[l]=vdata[l][k][j][i];
					
					if(buf[l]==undef){
						for(int ll=0;ll<t    ;ll++)  vdata[l][k][j][i]=undef;
						for(int ll=0;ll<cycle;ll++) nvdata[l][k][j][i]=undef;
						
						continue bk;
					}
				}
				
				float[] means=FilterModel.cycleFilter(buf,cycle);
				
				for(int l=0;l<t    ;l++)  vdata[l][k][j][i]=buf[l];
				for(int l=0;l<cycle;l++) nvdata[l][k][j][i]=means[l];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	bk:		for(int i=0;i<x;i++){
				for(int l=0;l<t;l++)
				if(vdata[k][j][i][l]==undef){
					for(int ll=0;ll<t    ;ll++)  vdata[k][j][i][ll]=undef;
					for(int ll=0;ll<cycle;ll++) nvdata[k][j][i][ll]=undef;
					
					continue bk;
				}
				
				FilterModel.cycleFilter(vdata[k][j][i],nvdata[k][j][i],cycle);
			}
		}
		
		nr.getTRange()[0]=r.getTRange()[0];
		nr.getTRange()[1]=nr.getTRange()[0]-1+cycle;
		
		nr.setZRange(r);	nr.setYRange(r);	nr.setXRange(r);
		
		return nv;
	}
	
	
	/**
     * monthly mean, daily data only
     *
     * @param	v		a give variable in daily format
     * @param	syear	start year of the data
     */
	public static Variable monthlyMean(Variable v,int syear){
		int t=v.getTCount(),z=v.getZCount(),y=v.getYCount(),x=v.getXCount();
		
		if(t<365)
			throw new IllegalArgumentException("daily data only, length should larger than 365");
		
		int yearc=t/365;	float undef=v.getUndef();
		
		int[] dcounts=new int[yearc];
		int[] stags  =new int[yearc];
		
		float[] ltmp=new float[366];
		float[] ptmp=new float[365];
		
		if(MDate.isLeapYear(syear)) dcounts[0]=366;
		else dcounts[0]=365;
		
		/*** calculate the tags ***/
		for(int i=syear+1;i<=syear+yearc-1;i++)
		if(MDate.isLeapYear(i)){
			dcounts[i-syear]=366;
			stags[i-syear]=stags[i-1-syear]+365;
			
		}else{
			dcounts[i-syear]=365;
			if(MDate.isLeapYear(i-1)) stags[i-syear]=stags[i-1-syear]+366;
			else stags[i-syear]=stags[i-1-syear]+365;
		}
		
		if(stags[yearc-1]+dcounts[yearc-1]!=t) throw new IllegalArgumentException("invalid data length");
		
		final int[][] pp={
			{  0, 30},{ 31, 58},{ 59, 89},{ 90,119},{120,150},{151,180},
			{181,211},{212,242},{243,272},{273,303},{304,333},{334,364}
		};
		
		final int[][] lp={
			{  0, 30},{ 31, 59},{ 60, 90},{ 91,120},{121,151},{152,181},
			{182,212},{213,243},{244,273},{274,304},{305,334},{335,365}
		};
		
		Variable m=new Variable(v.getName(),false,new Range(12*yearc,z,y,x));
		
		float[][][][] vdata=v.getData();
		float[][][][] mdata=m.getData();
		
		for(int ll=0;ll<yearc;ll++){
			if(dcounts[ll]==365){
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					System.arraycopy(vdata[k][j][i],stags[ll],ptmp,0,365);
					
					for(int lm=0;lm<12;lm++){
						int count=0;
						
						for(int lll=pp[lm][0];lll<=pp[lm][1];lll++)
						if(ptmp[lll]!=undef){ mdata[k][j][i][ll*12+lm]+=ptmp[lll]; count++;}
						
						mdata[k][j][i][ll*12+lm]/=count;
					}
				}
				
			}else{
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					System.arraycopy(vdata[k][j][i],stags[ll],ltmp,0,366);
					
					for(int lm=0;lm<12;lm++){
						int count=0;
						
						for(int lll=lp[lm][0];lll<=lp[lm][1];lll++)
						if(ltmp[lll]!=undef){ mdata[k][j][i][ll*12+lm]+=ltmp[lll]; count++;}
						
						mdata[k][j][i][ll*12+lm]/=count;
					}
				}
			}
		}
		
		return m;
	}
	
	
	/**
     * standardization of the temporal series of a given Variable
     *
     * @param	v	a given Variable
     */
	public static void temporalStandardization(Variable v){
		float undef=v.getUndef();
		
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			float[] tmp=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) tmp[l]=vdata[l][k][j][i];
				
				standardize(tmp,undef);

				for(int l=0;l<t;l++) vdata[l][k][j][i]=tmp[l];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) standardize(vdata[k][j][i],undef);
		}
	}
	
	
	/**
     * removal of trend along T-dimension
     *
     * @param	v		a given Variable
     * @param	order	1 for linear trend, 2 for quadratic trend...
     */
	public static void removeTrend(Variable v,int order){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float undef=v.getUndef();
		float[][][][] vdata=v.getData();
		
		PolynomialFitter pf=new PolynomialFitter(order,t);
		
		if(v.isTFirst()){
			float[] tmp=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cont:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++){
					tmp[l]=vdata[l][k][j][i];
					
					if(tmp[l]==undef){
						for(int ll=0;ll<t;ll++) vdata[ll][k][j][i]=undef;
						continue cont;
					}
				}
				
				pf.fit(tmp);
				
				float[] fitting=pf.cValues();

				for(int l=0;l<t;l++) vdata[l][k][j][i]-=fitting[l];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
	cont:	for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) if(vdata[k][j][i][l]==undef){
					for(int ll=0;ll<t;ll++) vdata[k][j][i][ll]=undef;
					continue cont;
				}
				
				pf.fit(vdata[k][j][i]);
				
				float[] fitting=pf.cValues();

				for(int l=0;l<t;l++) vdata[k][j][i][l]-=fitting[l];
			}
		}
	}
	
	public static void removeLinearTrend(Variable v){ removeTrend(v,1);}
	
	public static void removeQuadraticTrend(Variable v){ removeTrend(v,2);}
	
	
	/**
     * Remove the signal of variable ix from variable vv by the maximum covariance method,
     * Reference: Soon-Il An (2003, JC).
     * Notice that after calling this method, vv will change (no new memory allocated).
     *
     * @param	vv	a given Variable contains the signal of i
     * @param	ix	just a time series, such as index of ENSO
     */
	public static void maxCovarianceFilter(Variable vv,Variable ix){
		if(ix.getXCount()!=1||ix.getYCount()!=1||ix.getZCount()!=1)
			throw new IllegalArgumentException("index should be a time series only");
		
		if(vv.getTCount()!=ix.getTCount())
			throw new IllegalArgumentException("t-lengths are not same");
		
		float undef=vv.getUndef();
		
		int t=vv.getTCount(),	z=vv.getZCount(),	y=vv.getYCount(),	x=vv.getXCount();
		
		float[][][][] vdata=vv.getData();
		float[][][][] idata=ix.getData();
		
		if(vv.isTFirst()){
			float[] idx=new float[t];
			float[] tmp=new float[t];
			
			for(int l=0;l<t;l++) idx[l]=idata[l][0][0][0];
			
			float var=cVariance(idx,undef);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vdata[0][k][j][i]!=undef){
				for(int l=0;l<t;l++) tmp[l]=vdata[l][k][j][i];
				
				float cov=cCovariance(tmp,idx,undef);
				
				for(int l=0;l<t;l++) vdata[l][k][j][i]-=idx[l]*cov/var;
			}
			
		}else{
			float var=cVariance(idata[0][0][0],undef);
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(vdata[k][j][i][0]!=undef){
				float cov=cCovariance(vdata[k][j][i],idata[0][0][0],undef);
				
				for(int l=0;l<t;l++) vdata[k][j][i][l]-=idata[0][0][0][l]*cov/var;
			}
		}
	}
	
	
	/**
     * Running mean along T-dimension
     *
     * @param	v	a given variable
     * @param	p	points of running length (should be odd)
     */
	public static void TRunningMean(Variable v,int p){
		int t=v.getTCount(),z=v.getZCount();
		int y=v.getYCount(),x=v.getXCount();
		
		float undef=v.getUndef();
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			float[] buf=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) buf[l]=vdata[l][k][j][i];
				
				float[] re=FilterModel.runningMean(buf,p,undef);
				
				for(int l=0;l<t;l++) vdata[l][k][j][i]=re[l];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				float[] re=FilterModel.runningMean(vdata[k][j][i],p,undef);
				
				System.arraycopy(re,0,vdata[k][j][i],0,t);
			}
		}
	}
	
	
	/**
     * Band filter of one order Butterworth
     *
     * @param	v	a given Variable
     * @param	t1	lower boundry of time
     * @param	t2	upper boundry of time
     *
     * @return	re	result of the filtering
     */
	public static void ButterworthFilter(Variable v,float t1,float t2){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float[][][][] vdata=v.getData();
		
		float[] buf1=new float[t];
		
		if(v.isTFirst()){
			float[] buf2=new float[t];
			
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				for(int l=0;l<t;l++) buf1[l]=vdata[l][j][k][i];
				
				FilterModel.ButterworthFilter(buf1,buf2,t1,t2);
				
				for(int l=0;l<t;l++) vdata[l][j][k][i]=buf2[l];
			}
			
		}else{
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++){
				FilterModel.ButterworthFilter(vdata[k][j][i],buf1,t1,t2);
				
				System.arraycopy(buf1,0,vdata[k][j][i],0,t);
			}
		}
	}
	
	
	/**
     * FFT filter along dimension D
     *
     * @param	v	a given Variable
     * @param	D	Dimension
     * @param	Ks	wavenumbers that are removed after filering, start from 1 and end at length/2
     */
	public static void FFTFilter(Variable v,Dimension D,int... Ks){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float[][][][] vdata=v.getData();
		
		switch(D){
		case T:
			if(v.isTFirst()){
				float[] buf=new float[t];
				float[] res=new float[t];
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					for(int l=0;l<t;l++) buf[l]=vdata[l][k][j][i];
					FilterModel.FFTFilter(buf,res,Ks);
					for(int l=0;l<t;l++) vdata[l][k][j][i]=res[l];
				}
				
			}else{
				float[] res=new float[t];
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					FilterModel.FFTFilter(vdata[k][j][i],res,Ks);
					System.arraycopy(res,0,vdata[k][j][i],0,t);
				}
			}
			break;
			
		case X:
			if(v.isTFirst()){
				float[] res=new float[x];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					FilterModel.FFTFilter(vdata[l][k][j],res,Ks);
					System.arraycopy(res,0,vdata[l][k][j],0,x);
				}
				
			}else{
				float[] buf=new float[x];
				float[] res=new float[x];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					for(int i=0;i<x;i++) buf[i]=vdata[k][j][i][l];
					FilterModel.FFTFilter(buf,res,Ks);
					for(int i=0;i<x;i++) vdata[k][j][i][l]=res[i];
				}
			}
			break;
			
		case Y:
			if(v.isTFirst()){
				float[] buf=new float[y];
				float[] res=new float[y];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					for(int j=0;j<y;j++) buf[j]=vdata[l][k][j][i];
					FilterModel.FFTFilter(buf,res,Ks);
					for(int j=0;j<y;j++) vdata[l][k][j][i]=res[j];
				}
				
			}else{
				float[] buf=new float[y];
				float[] res=new float[y];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					for(int j=0;j<y;j++) buf[j]=vdata[k][j][i][l];
					FilterModel.FFTFilter(buf,res,Ks);
					for(int j=0;j<y;j++) vdata[k][j][i][l]=res[j];
				}
			}
			break;
			
		case Z:
			if(v.isTFirst()){
				float[] buf=new float[z];
				float[] res=new float[z];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					for(int k=0;k<z;k++) buf[k]=vdata[l][k][j][i];
					FilterModel.FFTFilter(buf,res,Ks);
					for(int k=0;k<z;k++) vdata[l][k][j][i]=res[k];
				}
				
			}else{
				float[] buf=new float[z];
				float[] res=new float[z];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					for(int k=0;k<z;k++) buf[k]=vdata[k][j][i][l];
					FilterModel.FFTFilter(buf,res,Ks);
					for(int k=0;k<z;k++) vdata[k][j][i][l]=res[k];
				}
			}
			break;

		default: throw new IllegalArgumentException("not supported dimension: "+D);
		}
	}
	
	/**
     * Fourier filter along dimension D
     *
     * @param	v	a given Variable
     * @param	D	Dimension
     * @param	Ts	periods of harmonics
     */
	public static void FourierFilter(Variable v,Dimension D,float... Ts){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float[][][][] vdata=v.getData();
		
		switch(D){
		case T:
			if(v.isTFirst()){
				float[] buf=new float[t];
				float[] res=new float[t];
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					for(int l=0;l<t;l++) buf[l]=vdata[l][k][j][i];
					FilterModel.FourierFilter(buf,res,Ts);
					for(int l=0;l<t;l++) vdata[l][k][j][i]=res[l];
				}
				
			}else{
				float[] res=new float[t];
				
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					FilterModel.FourierFilter(vdata[k][j][i],res,Ts);
					System.arraycopy(res,0,vdata[k][j][i],0,t);
				}
			}
			break;
			
		case X:
			if(v.isTFirst()){
				float[] res=new float[x];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					FilterModel.FourierFilter(vdata[l][k][j],res,Ts);
					System.arraycopy(res,0,vdata[l][k][j],0,x);
				}
				
			}else{
				float[] buf=new float[x];
				float[] res=new float[x];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int j=0;j<y;j++){
					for(int i=0;i<x;i++) buf[i]=vdata[k][j][i][l];
					FilterModel.FourierFilter(buf,res,Ts);
					for(int i=0;i<x;i++) vdata[k][j][i][l]=res[i];
				}
			}
			break;
			
		case Y:
			if(v.isTFirst()){
				float[] buf=new float[y];
				float[] res=new float[y];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					for(int j=0;j<y;j++) buf[j]=vdata[l][k][j][i];
					FilterModel.FourierFilter(buf,res,Ts);
					for(int j=0;j<y;j++) vdata[l][k][j][i]=res[j];
				}
				
			}else{
				float[] buf=new float[y];
				float[] res=new float[y];
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++){
					for(int j=0;j<y;j++) buf[j]=vdata[k][j][i][l];
					FilterModel.FourierFilter(buf,res,Ts);
					for(int j=0;j<y;j++) vdata[k][j][i][l]=res[j];
				}
			}
			break;
			
		case Z:
			if(v.isTFirst()){
				float[] buf=new float[z];
				float[] res=new float[z];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					for(int k=0;k<z;k++) buf[k]=vdata[l][k][j][i];
					FilterModel.FourierFilter(buf,res,Ts);
					for(int k=0;k<z;k++) vdata[l][k][j][i]=res[k];
				}
				
			}else{
				float[] buf=new float[z];
				float[] res=new float[z];
				
				for(int l=0;l<t;l++)
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++){
					for(int k=0;k<z;k++) buf[k]=vdata[k][j][i][l];
					FilterModel.FourierFilter(buf,res,Ts);
					for(int k=0;k<z;k++) vdata[k][j][i][l]=res[k];
				}
			}
			break;

		default: throw new IllegalArgumentException("not supported dimension: "+D);
		}
	}
	
	
	/**
     * Shuman-Shapiro's 2D spatial filter
     *
     * @param	v	a given Variable
     * @param	s1	argument in x direction
     * @param	s2	argument in x direction
     * @param	s3	argument in y direction
     * @param	s4	argument in y direction
     */
	public static void ShumanShapiroFilter(Variable v,float s1,float s2,float s3,float s4){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float[][][][] vdata=v.getData();
		
		float[][] buf=new float[y][x];
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++) System.arraycopy(vdata[l][k][j],0,buf[j],0,x);
				
				for(int j=2;j<y-2;j++)
				for(int i=2;i<x-2;i++){
					vdata[l][k][j][i]=
						((1-s1)*(1-s2)+s1*s2/2)*((1-s3)*(1-s4)+s3*s4/2)*buf[j][i]+
						0.5f*(s1*(1-s2)+s2*(1-s1))*((1-s3)*(1-s4)+s3*s4/2)*(buf[j][i-1]+buf[j][i+1])+
						0.5f*((1-s1)*(1-s2)+s1*s2/2)*(s3*(1-s4)+s4*(1-s3))*(buf[j-1][i]+buf[j+1][i])+
						0.25f*(s1*(1-s2)+s2*(1-s1))*(s3*(1-s4)+s4*(1-s3))*(buf[j-1][i-1]+buf[j-1][i+1]+buf[j+1][i-1]+buf[j+1][i+1])+
						0.25f*s1*s2*((1-s3)*(1-s4)+s3*s4/2)*(buf[j][i-2]+buf[j][i+2])+
						0.25f*s3*s4*((1-s1)*(1-s2)+s1*s2/2)*(buf[j-2][i]+buf[j+2][i])+
						0.125f*s1*s2*(s3*(1-s4)+s4*(1-s3))*(buf[j-1][i-2]+buf[j-1][i+2]+buf[j+1][i-2]+buf[j+1][i+2])+
						0.125f*s3*s4*(s1*(1-s2)+s2*(1-s1))*(buf[j-2][i-1]+buf[j-2][i+1]+buf[j+2][i-1]+buf[j+2][i+1])+
						0.0625f*s1*s2*s3*s4*(buf[j-2][i-2]+buf[j-2][i+2]+buf[j+2][i-2]+buf[j+2][i+2]);
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) buf[j][i]=vdata[k][j][i][l];
				
				for(int j=2;j<y-2;j++)
				for(int i=2;i<x-2;i++){
					vdata[k][j][i][l]=
						((1-s1)*(1-s2)+s1*s2/2)*((1-s3)*(1-s4)+s3*s4/2)*buf[j][i]+
						0.5f*(s1*(1-s2)+s2*(1-s1))*((1-s3)*(1-s4)+s3*s4/2)*(buf[j][i-1]+buf[j][i+1])+
						0.5f*((1-s1)*(1-s2)+s1*s2/2)*(s3*(1-s4)+s4*(1-s3))*(buf[j-1][i]+buf[j+1][i])+
						0.25f*(s1*(1-s2)+s2*(1-s1))*(s3*(1-s4)+s4*(1-s3))*(buf[j-1][i-1]+buf[j-1][i+1]+buf[j+1][i-1]+buf[j+1][i+1])+
						0.25f*s1*s2*((1-s3)*(1-s4)+s3*s4/2)*(buf[j][i-2]+buf[j][i+2])+
						0.25f*s3*s4*((1-s1)*(1-s2)+s1*s2/2)*(buf[j-2][i]+buf[j+2][i])+
						0.125f*s1*s2*(s3*(1-s4)+s4*(1-s3))*(buf[j-1][i-2]+buf[j-1][i+2]+buf[j+1][i-2]+buf[j+1][i+2])+
						0.125f*s3*s4*(s1*(1-s2)+s2*(1-s1))*(buf[j-2][i-1]+buf[j-2][i+1]+buf[j+2][i-1]+buf[j+2][i+1])+
						0.0625f*s1*s2*s3*s4*(buf[j-2][i-2]+buf[j-2][i+2]+buf[j+2][i-2]+buf[j+2][i+2]);
				}
			}
		}
	}
	
	
	/**
     * 9-point smooth operator, similar to OpenGrads
     *
     * @param	v	a given Variable
     */
	public static void smooth9(Variable v){
		int t=v.getTCount(),	z=v.getZCount(),	y=v.getYCount(),	x=v.getXCount();
		
		float undef=v.getUndef();
		
		float[][][][] vdata=v.getData();
		
		float[][] buf=new float[y][x];
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++) System.arraycopy(vdata[l][k][j],0,buf[j],0,x);
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(buf[j][i]!=undef){
					float s=buf[j][i],w=1;
					
					if(buf[j+1][i  ]!=undef){ s+=buf[j+1][i  ]*0.5f; w+=0.5f;}
					if(buf[j  ][i+1]!=undef){ s+=buf[j  ][i+1]*0.5f; w+=0.5f;}
					if(buf[j-1][i  ]!=undef){ s+=buf[j-1][i  ]*0.5f; w+=0.5f;}
					if(buf[j  ][i-1]!=undef){ s+=buf[j  ][i-1]*0.5f; w+=0.5f;}
					if(buf[j+1][i-1]!=undef){ s+=buf[j+1][i-1]*0.3f; w+=0.3f;}
					if(buf[j-1][i-1]!=undef){ s+=buf[j-1][i-1]*0.3f; w+=0.3f;}
					if(buf[j+1][i+1]!=undef){ s+=buf[j+1][i+1]*0.3f; w+=0.3f;}
					if(buf[j-1][i+1]!=undef){ s+=buf[j-1][i+1]*0.3f; w+=0.3f;}
					
					vdata[l][k][j][i]=s/w;
				}
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++)
				for(int i=0;i<x;i++) buf[j][i]=vdata[k][j][i][l];
				
				for(int j=1;j<y-1;j++)
				for(int i=1;i<x-1;i++)
				if(buf[j][i]!=undef){
					float s=buf[j][i],w=1;
					
					if(buf[j+1][i  ]!=undef){ s+=buf[j+1][i  ]*0.5f; w+=0.5f;}
					if(buf[j  ][i+1]!=undef){ s+=buf[j  ][i+1]*0.5f; w+=0.5f;}
					if(buf[j-1][i  ]!=undef){ s+=buf[j-1][i  ]*0.5f; w+=0.5f;}
					if(buf[j  ][i-1]!=undef){ s+=buf[j  ][i-1]*0.5f; w+=0.5f;}
					if(buf[j+1][i-1]!=undef){ s+=buf[j+1][i-1]*0.3f; w+=0.3f;}
					if(buf[j-1][i-1]!=undef){ s+=buf[j-1][i-1]*0.3f; w+=0.3f;}
					if(buf[j+1][i+1]!=undef){ s+=buf[j+1][i+1]*0.3f; w+=0.3f;}
					if(buf[j-1][i+1]!=undef){ s+=buf[j-1][i+1]*0.3f; w+=0.3f;}
					
					vdata[k][j][i][l]=s/w;
				}
			}
		}
	}
	
	
	/** test
	public static void main(String arg[]){
		int len=324;
		
		float ym=10;
		float slp=0.01f;
		float ramp=0.2f;
		float a1=6, a2=3;
		float b1=8, b2=4;
		float T1=6, T2=24;
		
		float[] H1=miniufo.basic.ArrayUtil.newHarmonicSeriesWithCoeff(len,a1,b1,T1);
		float[] H2=miniufo.basic.ArrayUtil.newHarmonicSeriesWithCoeff(len,a2,b2,T2);
		float[] N =miniufo.basic.ArrayUtil.newGaussianRandomSeries(len,ym,ramp);
		
		float[] y=new float[len];
		
		for(int l=0;l<len;l++) y[l]=slp*l+H1[l]+H2[l]+N[l];
		
		float[] re=FilterModel.FourFilter(y,T1,T2);
		
		for(int l=0;l<len;l++)
		System.out.println(y[l]+"\t"+re[l]+"\t"+(N[l]+slp*l));
	}*/
}
