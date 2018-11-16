/**
 * @(#)SingleParticleStatResult.java	1.0 2013.07.21
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.statisticsModel;

import java.io.FileWriter;
import java.io.IOException;


/**
 * Result of single-particle statistics.
 *
 * @version 1.0, 2013.07.21
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class SingleParticleStatResult{
	//
	int    dt=0;	// delta-T in seconds
	int  tRad=0;	// maximum time lag
	
	int pseudoTracks=0;	// total data
	int noOfMaxLag  =0;	// observations at maximum positive lag
	int noOfMinLag  =0;	// observations at minimum negative lag
	
	float umn =0;	// bin mean
	float vmn =0;	// bin mean
	float ucfd=0;	// confidence interval = 2 * Usigma / sqrt(Nu*)
	float vcfd=0;	// confidence interval = 2 * Vsigma / sqrt(Nv*)
	float Ku  =0;	// zonal diffusivity
	float Kv  =0;	// meridional diffusivity
	float Tu  =0;	// zonal timescale
	float Tv  =0;	// meridional timescale
	float Lu  =0;	// zonal length scale
	float Lv  =0;	// meridional length scale
	
	  int[] num=null;	// mean as a function of tau
	
	float[]  um=null;	// Lagrangian mean
	float[]  vm=null;
	float[]  ua=null;	// departure
	float[]  va=null;
	
	float[] mlon=null;	// averaged longitude for all pseudo-tracks
	float[] mlat=null;	// averaged latitude for all pseudo-tracks
	
	float[] DXm=null;	// Lagrangian mean
	float[] DYm=null;
	float[] DXa=null;	// departure
	float[] DYa=null;
	
	float[] Pxx=null;	// autocovariance of velocity as a function of tau
	float[] Pyy=null;
	float[] Pxy=null;
	float[] Pyx=null;
	
	float[] Qxx=null;	// autocovariance of acceleration as a function of tau
	float[] Qyy=null;
	float[] Qxy=null;
	float[] Qyx=null;
	
	float[] Kxx=null;	// diffusivity as a function of tau
	float[] Kyy=null;
	float[] Kxy=null;
	float[] Kyx=null;
	
	float[] K11=null;	// principle-axe-decomposed diffusivity as a function of tau
	float[] K22=null;
	float[] ang=null;	// angle of principle diffusivity
	
	float[] Dxx=null;	// dispersion as a function of tau
	float[] Dyy=null;
	float[] Dxy=null;
	float[] Dyx=null;
	
	public static enum LagType{NegativeLag,PositiveLag,All}
	
	
	/**
	 * constructor
	 * 
	 * @param	lon1	start longitude for the bin
	 * @param	lat1	start latitude for the bin
	 * @param	lon2	end longitude for the bin
	 * @param	lat2	end latitude for the bin
	 * @param	tRad	maximum time lag
	 */
	public SingleParticleStatResult(int tRad){
		if(tRad<1) throw new IllegalArgumentException("tRad should be positive");
		
		this.tRad=tRad;
		
		num=new int  [2*tRad+1];
		um =new float[2*tRad+1];	ua =new float[2*tRad+1];
		vm =new float[2*tRad+1];	va =new float[2*tRad+1];
		
		DXm=new float[2*tRad+1];	DXa=new float[2*tRad+1];
		DYm=new float[2*tRad+1];	DYa=new float[2*tRad+1];
		
		Pxx=new float[2*tRad+1];	Pyy=new float[2*tRad+1];
		Pxy=new float[2*tRad+1];	Pyx=new float[2*tRad+1];
		
		Qxx=new float[2*tRad+1];	Qyy=new float[2*tRad+1];
		Qxy=new float[2*tRad+1];	Qyx=new float[2*tRad+1];
		
		Kxx=new float[2*tRad+1];	Kyy=new float[2*tRad+1];
		Kxy=new float[2*tRad+1];	Kyx=new float[2*tRad+1];
		K11=new float[2*tRad+1];	K22=new float[2*tRad+1];
		ang=new float[2*tRad+1];
		
		Dxx=new float[2*tRad+1];	Dyy=new float[2*tRad+1];
		Dxy=new float[2*tRad+1];	Dyx=new float[2*tRad+1];
		
		mlon=new float[2*tRad+1];	mlat=new float[2*tRad+1];
	}
	
	
	/**
	 * Copy the current result to a new one.
	 */
	public SingleParticleStatResult copy(){
		SingleParticleStatResult spsr=new SingleParticleStatResult(tRad);
		
		spsr.dt  =dt;
		spsr.umn =umn;
		spsr.vmn =vmn;
		spsr.ucfd=ucfd;
		spsr.vcfd=vcfd;
		spsr.Ku  =Ku;
		spsr.Kv  =Kv;
		spsr.Tu  =Tu;
		spsr.Tv  =Tv;
		spsr.Lu  =Lu;
		spsr.Lv  =Lv;
		spsr.pseudoTracks=pseudoTracks;
		spsr.noOfMaxLag  =noOfMaxLag;
		spsr.noOfMinLag  =noOfMinLag;
		
		if(num !=null) spsr.num =num.clone();
		if(um  !=null) spsr.um  = um.clone();
		if(vm  !=null) spsr.vm  = vm.clone();
		if(ua  !=null) spsr.ua  = ua.clone();
		if(va  !=null) spsr.va  = va.clone();
		if(mlon!=null) spsr.mlon=mlon.clone();
		if(mlat!=null) spsr.mlat=mlat.clone();
		
		if(DXm !=null) spsr.DXm =DXm.clone();
		if(DYm !=null) spsr.DYm =DYm.clone();
		if(DXa !=null) spsr.DXa =DXa.clone();
		if(DYa !=null) spsr.DYa =DYa.clone();
		
		if(Pxx !=null) spsr.Pxx =Pxx.clone();
		if(Pyy !=null) spsr.Pyy =Pyy.clone();
		if(Pxy !=null) spsr.Pxy =Pxy.clone();
		if(Pyx !=null) spsr.Pyx =Pyx.clone();
		
		if(Qxx !=null) spsr.Qxx =Qxx.clone();
		if(Qyy !=null) spsr.Qyy =Qyy.clone();
		if(Qxy !=null) spsr.Qxy =Qxy.clone();
		if(Qyx !=null) spsr.Qyx =Qyx.clone();
		
		if(Kxx !=null) spsr.Kxx =Kxx.clone();
		if(Kyy !=null) spsr.Kyy =Kyy.clone();
		if(Kxy !=null) spsr.Kxy =Kxy.clone();
		if(Kyx !=null) spsr.Kyx =Kyx.clone();
		
		if(K11 !=null) spsr.K11 =K11.clone();
		if(K22 !=null) spsr.K22 =K22.clone();
		if(ang !=null) spsr.ang =ang.clone();
		
		if(Dxx !=null) spsr.Dxx =Dxx.clone();
		if(Dyy !=null) spsr.Dyy =Dyy.clone();
		if(Dxy !=null) spsr.Dxy =Dxy.clone();
		if(Dyx !=null) spsr.Dyx =Dyx.clone();
		
		return spsr;
	}
	
	
	/**
	 * Average the result from str to end lag for both positive and negtive time lags.
	 * 
	 * @param	str			start index
	 * @param	end			end index
	 * @param	minTracks	minimum number of pseudotrack
	 * @param	type		mean of positive lags or negative lags or all
	 * 
	 * @return	mean		[0,1] for [Kxx,Kyy], [2,3] for [Tx,Ty] and [4,5] for [Lx,Ly]
	 */
	public float[] getMean(int str,int end,int minTracks){ return getMean(str,end,minTracks,LagType.All);}
	
	public float[] getMean(int str,int end,int minTracks,LagType type){
		switch(type){
		case All: return getMeanAll(str,end,minTracks);
		case NegativeLag: return getMeanNegative(str,end,minTracks);
		case PositiveLag: return getMeanPositive(str,end,minTracks);
		default: throw new IllegalArgumentException("unsupported LagType: "+type);
		}
	}
	
	/**
	 * Get the maximum result from str to end lag for both positive and negtive time lags.
	 * 
	 * @param	str			start index
	 * @param	end			end index
	 * @param	minTracks	minimum number of pseudotrack
	 * @param	type		maximum of positive lags or negative lags or all
	 * 
	 * @return	max			[0,1] for [Kxx,Kyy], [2,3] for [Tx,Ty] and [4,5] for [Lx,Ly]
	 */
	public float[] getMax(int str,int end,int minTracks){ return getMax(str,end,minTracks,LagType.All);}
	
	public float[] getMax(int str,int end,int minTracks,LagType type){
		switch(type){
		case All: return getMaxAll(str,end,minTracks);
		case NegativeLag: return getMaxNegative(str,end,minTracks);
		case PositiveLag: return getMaxPositive(str,end,minTracks);
		default: throw new IllegalArgumentException("unsupported LagType: "+type);
		}
	}
	
	
	/**
	 * write out to a specific file
	 */
	public void toFile(String fname){
		try(FileWriter fw=new FileWriter(fname)){
			fw.write(toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	/**
	 * used to print out
	 */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		
		sb.append(String.format(
			" total count:%8d (%7.2f drifter-years), "+
			"%8d for max lag and %8d for min lag\n",
			pseudoTracks,pseudoTracks/(4f*365.2f),noOfMaxLag,noOfMinLag
		));
		
		sb.append(String.format(
			" um: %6.3f +- %6.3f (std., cm/s), Tu: %5.2f (days), Lu: %6.2f (km), Ku: %7.3f (10^7 cm^2/s)\n"+
			" vm: %6.3f +- %6.3f (std., cm/s), Tv: %5.2f (days), Lv: %6.2f (km), Kv: %7.3f (10^7 cm^2/s)\n",
			umn*100,ucfd*100,Tu,Lu,Ku,
			vmn*100,vcfd*100,Tv,Lv,Kv
		));
		
		sb.append(
			"  timelag samples   "+
			"um     vm      ua     va      "+
			"DXm     DYm    DXa     DYa    "+
			"Pxx     Pyy     Pxy     Pyx     "+
			"Qxx     Qyy     Qxy     Qyx     "+
			"Kxx     Kyy     Kxy     Kyx     K11     K22    Angle   AngM     "+
			"Dxx     Dyy      Dxy     Dyx    "+
			"mlon    mlat\n"
		);
		
		for(int l=0,L=num.length;l<L;l++)
		sb.append(String.format(
			" %7.2f %8d "+
			"%6.3f %6.3f %6.3f %6.3f " +
			"%7.2f %7.2f %7.2f %7.2f "+
			"%7.3f %7.3f %7.3f %7.3f "+
			"%7.3f %7.3f %7.3f %7.3f "+
			"%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f "+
			"%7.2f %7.2f %7.2f %7.2f "+
			"%7.2f %7.2f\n",
			(l-tRad)/(86400f/dt),num[l],
			um[l]*100f,vm[l]*100f,ua[l]*100f,va[l]*100f,
			DXm[l]/1000f,DYm[l]/1000f,DXa[l]/1000f,DYa[l]/1000f,
			Pxx[l]*1e4f,Pyy[l]*1e4f,Pxy[l]*1e4f,Pyx[l]*1e4f,
			Qxx[l]*1e15f,Qyy[l]*1e15f,Qxy[l]*1e15f,Qyx[l]*1e15f,
			Kxx[l]/1e3f,Kyy[l]/1e3f,Kxy[l]/1e3f,Kyx[l]/1e3f,K11[l]/1e3f,K22[l]/1e3f,Math.toDegrees(ang[l]),(Kyx[l]-Kxy[l])/1e3f,
			Dxx[l]/1e6f,Dyy[l]/1e6f,Dxy[l]/1e6f,Dyx[l]/1e6f,
			mlon[l],mlat[l]
		));
		
		return sb.toString();
	}
	
	
	/*** helper methods ***/
	private void checkIndices(int str,int end){
		if(str<0||end>tRad||str>end) throw new IllegalArgumentException("0 <= str <= end <= tRad");
	}
	
	private float getAbsMax(float[] data,int str,int end){
		float re=0;
		
		for(int i=str;i<=end;i++) if(!Float.isNaN(data[i])){
			float abs=Math.abs(data[i]);
			
			if(abs>re) re=abs;
		}
		
		return re;
	}
	
	
	private float[] getMeanAll(int str,int end,int minTracks){
		checkIndices(str,end);
		
		if(pseudoTracks<minTracks||noOfMaxLag<minTracks*0.2f||noOfMinLag<minTracks*0.2f){
			System.out.println("no enough tracks pseudoTracks("+pseudoTracks+") noOfMaxLag("+noOfMaxLag+") noOfMinLag("+noOfMinLag+") minTracks("+minTracks+")");
			return null;
		}
		
		int cK=0;
		
		float KxxM=0,KyyM=0,K11M=0,K22M=0,AngM=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(Kxx[i])){ KxxM+=Kxx[i]; cK++;}
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(Kxx[i])){ KxxM-=Kxx[i]; cK++;}
		
		KxxM/=cK;	cK=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(Kyy[i])){ KyyM+=Kyy[i]; cK++;}
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(Kyy[i])){ KyyM-=Kyy[i]; cK++;}
		
		KyyM/=cK;	cK=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(K11[i])){ K11M+=K11[i]; cK++;}
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(K11[i])){ K11M-=K11[i]; cK++;}
		
		K11M/=cK;	cK=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(K22[i])){ K22M+=K22[i]; cK++;}
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(K22[i])){ K22M-=K22[i]; cK++;}
		
		K22M/=cK;	cK=0;
		
		float angP=0,angN=0;	int cP=0,cN=0;
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(ang[i])){ angP+=ang[i]; cP++;}
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(ang[i])){ angN+=ang[i]; cN++;}
		
		angP/=cP;
		angN/=cN;
		
		if(angN!=0){
			angN+=Math.PI/2.0;
			if(angN>Math.PI/2.0) angN-=Math.PI;
		}
		
		AngM=(angP+angN)/2f;
		
		float[] mean=new float[9];
		
		mean[0]=KxxM/1e3f;
		mean[1]=KyyM/1e3f;
		
		mean[2]=KxxM/Pxx[tRad]/86400f;
		mean[3]=KyyM/Pyy[tRad]/86400f;
		
		mean[4]=KxxM/(float)Math.sqrt(Pxx[tRad])/1000f;
		mean[5]=KyyM/(float)Math.sqrt(Pyy[tRad])/1000f;
		
		mean[6]=K11M/1e3f;
		mean[7]=K22M/1e3f;
		mean[8]=AngM;
		
		return mean;
	}
	
	private float[] getMeanPositive(int str,int end,int minTracks){
		checkIndices(str,end);
		
		if(pseudoTracks<minTracks||noOfMaxLag<minTracks*0.2f){
			System.out.println("no enough tracks pseudoTracks("+pseudoTracks+") noOfMaxLag("+noOfMaxLag+") minTracks("+minTracks+")");
			return null;
		}
		
		int cK=0;
		
		float KxxM=0,KyyM=0,K11M=0,K22M=0,AngM=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(Kxx[i])){ KxxM+=Kxx[i]; cK++;}
		
		KxxM/=cK;	cK=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(Kyy[i])){ KyyM+=Kyy[i]; cK++;}
		
		KyyM/=cK;	cK=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(K11[i])){ K11M+=K11[i]; cK++;}
		
		K11M/=cK;	cK=0;
		
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(K22[i])){ K22M+=K22[i]; cK++;}
		
		K22M/=cK;	cK=0;
		
		float angP=0;	int cP=0;
		for(int i=tRad+str,I=tRad+end;i<=I;i++) if(!Float.isNaN(ang[i])){ angP+=ang[i]; cP++;}
		
		angP/=cP;
		
		AngM=angP;
		
		float[] mean=new float[9];
		
		mean[0]=KxxM/1e3f;
		mean[1]=KyyM/1e3f;
		
		mean[2]=KxxM/Pxx[tRad]/86400f;
		mean[3]=KyyM/Pyy[tRad]/86400f;
		
		mean[4]=KxxM/(float)Math.sqrt(Pxx[tRad])/1000f;
		mean[5]=KyyM/(float)Math.sqrt(Pyy[tRad])/1000f;
		
		mean[6]=K11M/1e3f;
		mean[7]=K22M/1e3f;
		mean[8]=AngM;
		
		return mean;
	}
	
	private float[] getMeanNegative(int str,int end,int minTracks){
		checkIndices(str,end);
		
		if(pseudoTracks<minTracks||noOfMinLag<minTracks*0.2f){
			System.out.println("no enough tracks pseudoTracks("+pseudoTracks+") noOfMinLag("+noOfMinLag+") minTracks("+minTracks+")");
			return null;
		}
		
		int cK=0;
		
		float KxxM=0,KyyM=0,K11M=0,K22M=0,AngM=0;
		
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(Kxx[i])){ KxxM-=Kxx[i]; cK++;}
		
		KxxM/=cK;	cK=0;
		
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(Kyy[i])){ KyyM-=Kyy[i]; cK++;}
		
		KyyM/=cK;	cK=0;
		
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(K11[i])){ K11M-=K11[i]; cK++;}
		
		K11M/=cK;	cK=0;
		
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(K22[i])){ K22M-=K22[i]; cK++;}
		
		K22M/=cK;	cK=0;
		
		float angN=0;	int cN=0;
		for(int i=tRad-end,I=tRad-str;i<=I;i++) if(!Float.isNaN(ang[i])){ angN+=ang[i]; cN++;}
		
		angN/=cN;
		
		if(angN!=0){
			angN+=Math.PI/2.0;
			if(angN>Math.PI/2.0) angN-=Math.PI;
		}
		
		AngM=angN;
		
		float[] mean=new float[9];
		
		mean[0]=KxxM/1e3f;
		mean[1]=KyyM/1e3f;
		
		mean[2]=KxxM/Pxx[tRad]/86400f;
		mean[3]=KyyM/Pyy[tRad]/86400f;
		
		mean[4]=KxxM/(float)Math.sqrt(Pxx[tRad])/1000f;
		mean[5]=KyyM/(float)Math.sqrt(Pyy[tRad])/1000f;
		
		mean[6]=K11M/1e3f;
		mean[7]=K22M/1e3f;
		mean[8]=AngM;
		
		return mean;
	}
	
	
	private float[] getMaxAll(int str,int end,int minTracks){
		checkIndices(str,end);
		
		if(pseudoTracks<minTracks||noOfMaxLag<minTracks*0.2f||noOfMinLag<minTracks*0.2f){
			System.out.println("no enough tracks pseudoTracks("+pseudoTracks+") noOfMaxLag("+noOfMaxLag+") noOfMinLag("+noOfMinLag+") minTracks("+minTracks+")");
			return null;
		}
		
		int cK=0;
		
		float KxxM=0,KyyM=0,K11M=0,K22M=0,AngM=0,tmpMax=0;
		
		tmpMax=getAbsMax(Kxx,tRad+str,tRad+end); if(tmpMax!=0){ KxxM+=tmpMax; cK++;}
		tmpMax=getAbsMax(Kxx,tRad-end,tRad-str); if(tmpMax!=0){ KxxM+=tmpMax; cK++;}
		
		if(cK!=0) KxxM/=cK;	cK=0;
		
		tmpMax=getAbsMax(Kyy,tRad+str,tRad+end); if(tmpMax!=0){ KyyM+=tmpMax; cK++;}
		tmpMax=getAbsMax(Kyy,tRad-end,tRad-str); if(tmpMax!=0){ KyyM+=tmpMax; cK++;}
		
		if(cK!=0) KyyM/=cK;	cK=0;
		
		tmpMax=getAbsMax(K11,tRad+str,tRad+end); if(tmpMax!=0){ K11M+=tmpMax; cK++;}
		tmpMax=getAbsMax(K11,tRad-end,tRad-str); if(tmpMax!=0){ K11M+=tmpMax; cK++;}
		
		if(cK!=0) K11M/=cK;	cK=0;
		
		tmpMax=getAbsMax(K22,tRad+str,tRad+end); if(tmpMax!=0){ K22M+=tmpMax; cK++;}
		tmpMax=getAbsMax(K22,tRad-end,tRad-str); if(tmpMax!=0){ K22M+=tmpMax; cK++;}
		
		if(cK!=0) K22M/=cK;	cK=0;
		
		float angP=0,angN=0;	int cP=0,cN=0;
		tmpMax=getAbsMax(ang,tRad+str,tRad+end); if(tmpMax!=0){ angP+=tmpMax; cP++;}
		tmpMax=getAbsMax(ang,tRad-end,tRad-str); if(tmpMax!=0){ angN+=tmpMax; cN++;}
		
		if(cP!=0) angP/=cP;
		if(cN!=0) angN/=cN;
		
		if(angN!=0){
			angN+=Math.PI/2.0;
			if(angN>Math.PI/2.0) angN-=Math.PI;
		}
		
		AngM=(angP+angN)/2f;
		
		float[] extremes=new float[9];
		
		extremes[0]=KxxM/1e3f;
		extremes[1]=KyyM/1e3f;
		
		extremes[2]=KxxM/Pxx[tRad]/86400f;
		extremes[3]=KyyM/Pyy[tRad]/86400f;
		
		extremes[4]=KxxM/(float)Math.sqrt(Pxx[tRad])/1000f;
		extremes[5]=KyyM/(float)Math.sqrt(Pyy[tRad])/1000f;
		
		extremes[6]=K11M/1e3f;
		extremes[7]=K22M/1e3f;
		extremes[8]=AngM;
		
		return extremes;
	}
	
	private float[] getMaxPositive(int str,int end,int minTracks){
		checkIndices(str,end);
		
		if(pseudoTracks<minTracks||noOfMaxLag<minTracks*0.2f){
			System.out.println("no enough tracks pseudoTracks("+pseudoTracks+") noOfMaxLag("+noOfMaxLag+") minTracks("+minTracks+")");
			return null;
		}
		
		int cK=0;
		
		float KxxM=0,KyyM=0,K11M=0,K22M=0,AngM=0,tmpMax=0;
		
		tmpMax=getAbsMax(Kxx,tRad+str,tRad+end); if(tmpMax!=0){ KxxM+=tmpMax; cK++;}
		
		if(cK!=0) KxxM/=cK;	cK=0;
		
		tmpMax=getAbsMax(Kyy,tRad+str,tRad+end); if(tmpMax!=0){ KyyM+=tmpMax; cK++;}
		
		if(cK!=0) KyyM/=cK;	cK=0;
		
		tmpMax=getAbsMax(K11,tRad+str,tRad+end); if(tmpMax!=0){ K11M+=tmpMax; cK++;}
		
		if(cK!=0) K11M/=cK;	cK=0;
		
		tmpMax=getAbsMax(K22,tRad+str,tRad+end); if(tmpMax!=0){ K22M+=tmpMax; cK++;}
		
		if(cK!=0) K22M/=cK;	cK=0;
		
		float angP=0;	int cP=0;
		tmpMax=getAbsMax(ang,tRad+str,tRad+end); if(tmpMax!=0){ angP+=tmpMax; cP++;}
		
		if(cP!=0) angP/=cP;
		
		AngM=angP;
		
		float[] extremes=new float[9];
		
		extremes[0]=KxxM/1e3f;
		extremes[1]=KyyM/1e3f;
		
		extremes[2]=KxxM/Pxx[tRad]/86400f;
		extremes[3]=KyyM/Pyy[tRad]/86400f;
		
		extremes[4]=KxxM/(float)Math.sqrt(Pxx[tRad])/1000f;
		extremes[5]=KyyM/(float)Math.sqrt(Pyy[tRad])/1000f;
		
		extremes[6]=K11M/1e3f;
		extremes[7]=K22M/1e3f;
		extremes[8]=AngM;
		
		return extremes;
	}
	
	private float[] getMaxNegative(int str,int end,int minTracks){
		checkIndices(str,end);
		
		if(pseudoTracks<minTracks||noOfMinLag<minTracks*0.2f){
			System.out.println("no enough tracks pseudoTracks("+pseudoTracks+") noOfMinLag("+noOfMinLag+") minTracks("+minTracks+")");
			return null;
		}
		
		int cK=0;
		
		float KxxM=0,KyyM=0,K11M=0,K22M=0,AngM=0,tmpMax=0;
		
		tmpMax=getAbsMax(Kxx,tRad-end,tRad-str); if(tmpMax!=0){ KxxM+=tmpMax; cK++;}
		
		if(cK!=0) KxxM/=cK;	cK=0;
		
		tmpMax=getAbsMax(Kyy,tRad-end,tRad-str); if(tmpMax!=0){ KyyM+=tmpMax; cK++;}
		
		if(cK!=0) KyyM/=cK;	cK=0;
		
		tmpMax=getAbsMax(K11,tRad-end,tRad-str); if(tmpMax!=0){ K11M+=tmpMax; cK++;}
		
		if(cK!=0) K11M/=cK;	cK=0;
		
		tmpMax=getAbsMax(K22,tRad-end,tRad-str); if(tmpMax!=0){ K22M+=tmpMax; cK++;}
		
		if(cK!=0) K22M/=cK;	cK=0;
		
		float angN=0;	int cN=0;
		tmpMax=getAbsMax(ang,tRad-end,tRad-str); if(tmpMax!=0){ angN+=tmpMax; cN++;}
		
		if(cN!=0) angN/=cN;
		
		if(angN!=0){
			angN+=Math.PI/2.0;
			if(angN>Math.PI/2.0) angN-=Math.PI;
		}
		
		AngM=angN;
		
		float[] extremes=new float[9];
		
		extremes[0]=KxxM/1e3f;
		extremes[1]=KyyM/1e3f;
		
		extremes[2]=KxxM/Pxx[tRad]/86400f;
		extremes[3]=KyyM/Pyy[tRad]/86400f;
		
		extremes[4]=KxxM/(float)Math.sqrt(Pxx[tRad])/1000f;
		extremes[5]=KyyM/(float)Math.sqrt(Pyy[tRad])/1000f;
		
		extremes[6]=K11M/1e3f;
		extremes[7]=K22M/1e3f;
		extremes[8]=AngM;
		
		return extremes;
	}
	
	
	/** test
	public static void main(String arg[]){
		System.out.println("[235,23]".replaceAll("\\[|\\]",""));
	}*/
}
