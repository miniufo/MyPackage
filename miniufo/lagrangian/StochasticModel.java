/**
 * @(#)StochasticModel.java	1.0 2014.05.10
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.lagrangian;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import miniufo.concurrent.ConcurrentUtil;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Variable;
import miniufo.util.GridDataFetcher;
import miniufo.util.Region2D;
import miniufo.util.TicToc;
import static java.lang.Math.cos;
import static java.lang.Math.toDegrees;
import static java.lang.Math.toRadians;
import static miniufo.diagnosis.SpatialModel.EARTH_RADIUS;


/**
 * tracking a particle's trajectory using stochastic method
 *
 * @version 1.0, 2014.05.10
 * @author  MiniUFO
 * @since   MDK1.0
 */
public abstract class StochasticModel{
	//
	protected int dtRatio=0;
	
	protected boolean constMean=false;
	
	protected float dt=0;		// dt = deltaT/dtRatio, deltaT is defined in DataDescriptor
	protected float undef=0;	// undefined value
	
	protected Variable ubuf=null;
	protected Variable vbuf=null;
	
	protected BCType BCx=BCType.Landing;
	protected BCType BCy=BCType.Landing;
	
	protected Region2D region=null;
	protected DataDescriptor dd=null;
	protected GridDataFetcher gdf=null;
	
	protected Function<Record,StochasticParams> func=null;
	
	protected final Random rnd=new Random();
	
	public enum BCType{Periodic,Reflected,Landing}
	
	
	/**
	 * constructor
	 * 
	 * @param	dtRatio		delta time ratio, dt = deltaT / dtRatio, deltaT defined in DataDescriptor
	 * @param	constMean	whether to use time-invariant mean flow
	 * @param	dd			DataDescriptor describe the domain
	 * @param	BCx			zonal boundary condition
	 * @param	BCy			meridional boundary condition
	 * @param	func		a function that map the current Record to StochasticParams
	 */
	public StochasticModel
	(int dtRatio,boolean constMean,DataDescriptor dd,BCType BCx,BCType BCy,Function<Record,StochasticParams> func){
		if(dtRatio<1) throw new IllegalArgumentException("dtRatio should be at least 1");
		
		if(dd instanceof CsmDescriptor)
		throw new IllegalArgumentException("unsupport for CsmDescriptor");
		
		this.dt=dd.getDTDef()[0]/dtRatio;
		this.undef=dd.getUndef(null);
		this.dtRatio=dtRatio;
		this.constMean=constMean;
		this.dd=dd;
		this.BCx=BCx;
		this.BCy=BCy;
		this.func=func;
		this.region=dd.toRegion2D();
		
		gdf=new GridDataFetcher(dd);
	}
	
	
	public abstract int getOrder();
	
	
	/**
     * set the spatially temporally variated velocity fields
     *
     * @param	u		u name defined in DataDescriptor
     * @param	v		v name defined in DataDescriptor
     * @param	tstr	start time in DataDescriptor (start from 1)
     */
	public void setVelocityBuffer(String u,String v,int tstr){
		if(constMean){
			ubuf=gdf.prepareXYBuffer(u,tstr,1);
			vbuf=gdf.prepareXYBuffer(v,tstr,1);
			
		}else{
			ubuf=gdf.prepareXYTBuffer(u,1,tstr,2,0);
			vbuf=gdf.prepareXYTBuffer(v,1,tstr,2,0);
		}
	}
	
	
	/**
     * integrate all particles forwarding over deltaT defined in DataDescriptor
     *
     * @param	ps			a list of particles
	 * @param	appendRec	append the record (false for updating record)
     */
	public void integrateForward(List<Particle> ps,boolean appendRec){ for(Particle p:ps) forwardDeltaT(p,appendRec);}
	
	public void integrateForward(List<Particle> ps){ for(Particle p:ps) forwardDeltaT(p,false);}
	
	
	/**
     * deploy a particle initially located at (lon,lat) with a velocity specified
     * at that point (linear interpolation fashion), adding a random velocity and
     * acceleration fluctuations depending on a given stochastic model
     *
     * @param	lon			longitude (degree)
     * @param	lat			latitude (degree)
     * @param	initlen		initial length of particle
     * @param	params		stochastic parameters, including diffusivity and timescale
     */
	public Particle deployAt(String id,float lon,float lat){
		return deployAt(id,lon,lat,10);
	}
	
	public Particle deployAt(String id,float lon,float lat,int initLen){
		if(ubuf==null||vbuf==null)
		throw new IllegalArgumentException("call setVelocityBuffer() first");
		
		long time=dd.getTDef().getSamples()[gdf.getBufferTLevel()-1].getLongTime();
		float undef=dd.getUndef(null);
		float U=0,V=0;
		
		if(constMean){
			U=gdf.fetchXYBuffer(lon,lat,ubuf); if(U==undef) return null;	// zonal mean advective velocity
			V=gdf.fetchXYBuffer(lon,lat,vbuf); if(V==undef) return null;	// meridional mean advective velocity
			
		}else{
			U=gdf.fetchXYTBuffer(lon,lat,time,ubuf); if(U==undef) return null;	// zonal mean advective velocity
			V=gdf.fetchXYTBuffer(lon,lat,time,vbuf); if(V==undef) return null;	// meridional mean advective velocity
		}
		
		StochasticParams params=func.apply(new Record(time,lon,lat)); validateOrder(params);
		
		Record init=new Record(time,lon,lat,4);
		
		float[] rnds=spinupRandom(init,200);
		
		init.setData(0,U+rnds[0]);	// velX+resX
		init.setData(1,V+rnds[1]);	// velY+resY
		init.setData(2,rnds[2]);	// accX
		init.setData(3,rnds[3]);	// accY
		
		Particle p=new Particle(id,initLen,4,true);
		p.setAttachedDataNames("velX","velY","accX","accY");
		p.addRecord(init);
		
		return p;
	}
	
	
	/**
	 * deploy a patch of particles
	 * 
	 * @param	r			region to deploy
	 * @param	del			spacing of deployment in both directions (degree)
	 * @param	ensemble	number of particles deployed at the same point
	 * @param	initLen		initial length of the particle (pre-allocating memory)
	 * @param	func		function that maps the current Record to StochasticParams
	 */
	public List<Particle> deployPatch(Region2D r,float del,int ensemble,int initLen){
		int idx=0;
		
		List<Particle> ps=new ArrayList<>();
		
		int yc=Math.round((r.getYMax()-r.getYMin())/del);
		int xc=Math.round((r.getXMax()-r.getXMin())/del);
		
		if(((r.getYMax()-r.getYMin())%del)/del>0.99f) yc++;
		if(((r.getXMax()-r.getXMin())%del)/del>0.99f) xc++;
		
		for(int j=0;j<yc;j++)
		for(int i=0;i<xc;i++)
		for(int m=0;m<ensemble;m++){
			float lon=r.getXMin()+del*i;
			float lat=r.getYMin()+del*j;
			
			Particle p=deployAt(""+(1000000+idx++),lon,lat,initLen);
			
			if(p!=null) ps.add(p);
		}
		
		if(ps.size()>1000000)
		throw new IllegalArgumentException("more than 1,000,000 particles released (ID overflow)");
		
		return ps;
	}
	
	public List<Particle> deployPatch(Region2D r,float del,int ensemble){
		return deployPatch(r,del,ensemble,10);
	}
	
	
	/**
	 * simulate a list of particles
	 * 
	 * @param	ls			particles in a list
	 * @param 	intLen		length of integration (steps)
	 * @param	appendRec	append the record (false for updating record)
	 * @param	threads		number of threads
	 */
	public void simulateParticles(List<Particle> ls,String u,String v,int intLen,boolean appendRec,int threads){
		if(intLen<1) throw new IllegalArgumentException("steps should be positive");
		
		TicToc.tic("start tracking "+ls.size()+" particles");
		
		if(constMean){
			setVelocityBuffer(u,v,gdf.getBufferTLevel());
			
			int len=ls.size();
			int count=len/threads;
			int left =len%threads;
			
			int[] tags=new int[threads+1];
			
			for(int i=1;i<=left;i++) tags[i]=tags[i-1]+count+1;
			for(int i=left+1;i<threads;i++) tags[i]=tags[i-1]+count;
			tags[threads]=len;
			
			ExecutorService es=ConcurrentUtil.defaultExecutor();
			CompletionService<Void> cs=new ExecutorCompletionService<>(es);
			
			for(int i=0;i<threads;i++){
				final int ii=i;
				cs.submit(()->{for(int l=1;l<intLen;l++) integrateForward(ls.subList(tags[ii],tags[ii+1]),appendRec);},null);
			}
		    
		    try{for(int i=0;i<threads;i++) cs.take();}
		    catch(InterruptedException e){ e.printStackTrace(); System.exit(0);}
			
		}else for(int l=1;l<intLen;l++){
			if(threads!=1)
			throw new IllegalArgumentException("parallel computing not supported for non-constant mean");
			
			setVelocityBuffer(u,v,l);
			integrateForward(ls,appendRec);
		}
		
		TicToc.toc(TimeUnit.MINUTES);
	}
	
	public void simulateParticles(List<Particle> ls,String u,String v,int intLen){
		simulateParticles(ls,u,v,intLen,true,1);
	}
	
	
	/**
	 * shut down the tracking and close files
	 */
	public void shutdown(){ gdf.closeFile();}
	
	
	/*** helper methods ***/
	
	/**
	 * integrate forward over delta-T defined in DataDescriptor
	 * 
	 * @param	p			one particle
	 * @param	appendRec	append the record (false for updating record)
	 */
	protected void forwardDeltaT(Particle p,boolean appendRec){
		if(p.isFinished()) return;
		
		Record last=p.getRecord(p.getTCount()-1);
		
		for(int l=0;l<dtRatio;l++){
			Record now=forwardDT(last);
			
			// validate record
			if(now==null){ p.finish(); break;}
			
			last=now;
		}
		
		if(!p.isFinished()){
			if(appendRec) p.addRecord(last);
			else p.setRecord(p.getTCount()-1,last);
		}
	}
	
	/**
	 * integrate forward over dt and return a new Record
	 * null would be return if no valid Record
	 * 
	 * @param	record	initial record
	 * @param	l		time step (t-tag)
	 */
	protected Record forwardDT(Record init){ return forwardRK4(init,dt);}
	
	
	protected abstract Record forwardRK4(Record init,float dt);
	
	protected abstract float[] getRandom(Record init,Record esti);
	
	protected abstract float[] spinupRandom(Record init,int iter);
	
	
	/**
	 * Forward dt using position and time of Record init and the velocity of Record esti.
	 * The residual velocity and acceleration of Record esti would be check if encounting a BC.
	 * The returned Record contains new position and checked residual velocity and acceleration.
	 * 
	 * 
	 * @param	init	initial record
	 * @param	esti	estimated record
	 * @param	dt		forward time (s)
	 */
	protected Record forwardByMean(Record init,Record esti,float dt){
		if(init==null||esti==null) return null;
		
		float[] velm=fetchVelocity(esti.getTime(),esti.getXPos(),esti.getYPos());
		
		float lon0=init.getXPos();
		float lat0=init.getYPos();
		float dlon=(float)toDegrees(velm[0]*dt/(EARTH_RADIUS*cos(toRadians(lat0))));
		float dlat=(float)toDegrees(velm[1]*dt/EARTH_RADIUS);
		float lon1=lon0+dlon;
		float lat1=lat0+dlat;
		
		float[] res1=processBCs(lon1,lat1,0,0,0,0); if(res1==null) return null;
		
		Record re=new Record(init);
		
		long tim1=new MDate(init.getTime()).addSeconds(Math.round(dt)).getLongTime();
		
		float[] vel1=fetchVelocity(tim1,res1[0],res1[1]);
		
		if(vel1==null) return null;
		
		re.setTime(tim1);
		re.setXPos(res1[0]); re.setData(0,vel1[0]);
		re.setYPos(res1[1]); re.setData(1,vel1[1]);
		
		return re;
	}
	
	protected Record forwardByAll(Record init,Record esti,float resX,float resY,float dt,boolean cVel){
		if(init==null||esti==null) return null;
		
		float lon0=init.getXPos();
		float lat0=init.getYPos();
		
		float velXe=esti.getDataValue(0);
		float velYe=esti.getDataValue(1);
		float accXe=esti.getDataValue(2);
		float accYe=esti.getDataValue(3);
		
		float dlon=(float)toDegrees(velXe*dt/(EARTH_RADIUS*cos(toRadians(lat0))));
		float dlat=(float)toDegrees(velYe*dt/EARTH_RADIUS);
		float lon1=lon0+dlon;
		float lat1=lat0+dlat;
		
		float[] res1=processBCs(lon1,lat1,resX,resY,accXe,accYe); if(res1==null) return null;
		
		Record re=new Record(init);
		
		// new position
		re.setXPos(res1[0]);
		re.setYPos(res1[1]);
		
		if(cVel){ // new background mean and checked residual velocities
			long tim1=new MDate(init.getTime()).addSeconds(Math.round(dt)).getLongTime();
			
			float[] vel1=fetchVelocity(tim1,res1[0],res1[1]);
			
			if(vel1==null) return null;
			
			re.setTime(tim1);
			re.setData(0,vel1[0]+res1[2]); re.setData(2,res1[4]);
			re.setData(1,vel1[1]+res1[3]); re.setData(3,res1[5]);
		}
		
		return re;
	}
	
	protected float[] processBCs(float lon,float lat,float resX,float resY,float accX,float accY){
		if(!region.inXRange(lon)) switch(BCx){
			case Landing: return null;
			case Reflected:{
				float lonmax=region.getXMax();
				float lonmin=region.getXMin();
				
				if(lon>lonmax) lon=lonmax-(lon-lonmax);
				else lon=lonmin+(lonmin-lon);
				
				resX=-resX;
				accX=-accX;
				
				break;
			}
			case Periodic:{
				float lonmax=region.getXMax()+dd.getDXDef()[0];
				float lonmin=region.getXMin();
				
				if(lon>lonmax) lon-=lonmax-lonmin;
				else if(lon<lonmin) lon+=lonmax-lonmin;
				
				break;
			}
			default: throw new IllegalArgumentException("unsupported BoundaryType: "+BCx);
		}
		
		if(!region.inYRange(lat)) switch(BCy){
			case Landing: return null;
			case Reflected:{
				float latmax=region.getYMax();
				float latmin=region.getYMin();
				
				if(lat>latmax) lat=latmax-(lat-latmax);
				else lat=latmin+(latmin-lat);
				
				resY=-resY;
				accY=-accY;
				
				break;
			}
			case Periodic:
			default: throw new IllegalArgumentException("unsupported BoundaryType: "+BCy);
		}
		
		return new float[]{lon,lat,resX,resY,accX,accY};
	}
	
	protected float[] fetchVelocity(long tim,float lon,float lat){
		float velX1=0,velY1=0;
		
		if(constMean){
			if(BCx==BCType.Periodic){
				velX1=gdf.fetchXYBufferPeriodicX(lon,lat,ubuf); if(velX1==undef) return null;
				velY1=gdf.fetchXYBufferPeriodicX(lon,lat,vbuf); if(velY1==undef) return null;
				
			}else{
				velX1=gdf.fetchXYBuffer(lon,lat,ubuf); if(velX1==undef) return null;
				velY1=gdf.fetchXYBuffer(lon,lat,vbuf); if(velY1==undef) return null;
			}
		}else{
			if(BCx==BCType.Periodic){
				velX1=gdf.fetchXYTBufferPeriodicX(lon,lat,tim,ubuf); if(velX1==undef) return null;
				velY1=gdf.fetchXYTBufferPeriodicX(lon,lat,tim,vbuf); if(velY1==undef) return null;
				
			}else{
				velX1=gdf.fetchXYTBuffer(lon,lat,tim,ubuf); if(velX1==undef) return null;
				velY1=gdf.fetchXYTBuffer(lon,lat,tim,vbuf); if(velY1==undef) return null;
			}
		}
		
		return new float[]{velX1,velY1};
	}
	
	protected void validateOrder(StochasticParams sp){
		if(sp.getOrder()!=getOrder()) throw new IllegalArgumentException(
			"Order ("+sp.getOrder()+") of StochasticParams should be "+getOrder()
		);
	}
}
