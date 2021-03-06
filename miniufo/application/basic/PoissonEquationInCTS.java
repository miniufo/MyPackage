/**
 * @(#)PoissonEquationInCTS.java	1.0 2017.07.05
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.application.basic;

import miniufo.diagnosis.Variable;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInCartesianCoordinate;
import miniufo.application.advanced.EllipticEqSORSolver2D;
import miniufo.application.advanced.EllipticEqSORSolver2D.DimCombination;


/**
 * Poisson equation within horizontal plane in Cartesian coordinates.
 *
 * @version 1.0, 2017.07.05
 * @author  MiniUFO
 * @since   MDK1.0
 */
public class PoissonEquationInCTS extends EquationInCartesianCoordinate implements EllipticEquationInterface{
	//
	private DimCombination dimC=DimCombination.XY;
	
	
	/**
     * Constructor
     *
     * @param	ssm		initialized by spatial model in spherical coordinate
     */
	public PoissonEquationInCTS(CartesianSpatialModel csm){
		super(csm);
		
		if(!csm.isLinearModel()) System.out.println("\nNot a equal-space model");
	}
	
	public PoissonEquationInCTS(CartesianSpatialModel csm,DimCombination dimC){
		this(csm);
		
		this.dimC=dimC;
	}
	
	
	/*** getor and setor ***/
	public DimCombination getDimCombination(){ return dimC;}
	
	public void setDimComBination(DimCombination dimC){ this.dimC=dimC;}
	
	
	/**
     * Implement the methods in EllipticEquationInterface.
     */
	public Variable cAPrime(Variable v){
		if(dimC==DimCombination.XY) return cAPrimeXY(v);
		else return cAPrimeYZ(v);
	}
	
	public Variable cBPrime(Variable v){
		if(dimC==DimCombination.XY) return cBPrimeXY(v);
		else return cBPrimeYZ(v);
	}
	
	public Variable cCPrime(Variable v){
		if(dimC==DimCombination.XY) return cCPrimeXY(v);
		else return cCPrimeYZ(v);
	}
	
	public Variable cA(Variable v){
		if(dimC==DimCombination.XY) return cAXY(v);
		else return cAYZ(v);
	}
	
	public Variable cB(Variable v){
		if(dimC==DimCombination.XY) return cBXY(v);
		else return cBYZ(v);
	}
	
	public Variable cC(Variable v){
		if(dimC==DimCombination.XY) return cCXY(v);
		else return cCYZ(v);
	}
	
	
	/**
	 * Inverting the Poisson equation using SOR iteration.
	 * 
	 * @param	F	forcing on the r.h.s. of Poisson equation
	 */
	public Variable invertingBySOR(Variable F){
		EllipticEqSORSolver2D ees=new EllipticEqSORSolver2D(sm);
		
		Variable SS=new Variable("S",F);
		
		ees.setBCofX(BCx);
		ees.setBCofY(BCy);
		ees.setBCofZ(BCz);
		ees.setDimCombination(dimC);
		ees.setABC(cAPrime(F),cBPrime(F),cCPrime(F));
		ees.solve(SS,F);
		
		return SS;
	}
	
	
	/*** helper methods ***/
	private Variable cAPrimeXY(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Ap",v);
		A.setCommentAndUnit("coefficient A' of elliptic equation");
		A.setValue(1f);
		
		return A;
	}
	
	private Variable cAPrimeYZ(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Ap",v);
		A.setCommentAndUnit("coefficient A' of elliptic equation");
		A.setValue(1f);
		
		return A;
	}
	
	private Variable cBPrimeXY(Variable v){
		Variable B=new Variable("Bp",v);
		B.setCommentAndUnit("coefficient B' of elliptic equation");
		return B;
	}
	
	private Variable cBPrimeYZ(Variable v){
		Variable B=new Variable("Bp",v);
		B.setCommentAndUnit("coefficient B' of elliptic equation");
		return B;
	}
	
	private Variable cCPrimeXY(Variable v){
		assignSubDomainParams(v);
		
		Variable C=new Variable("Cp",v);
		C.setCommentAndUnit("coefficient C' of elliptic equation");
		C.setValue(1f);
		
		return C;
	}
	
	private Variable cCPrimeYZ(Variable v){
		assignSubDomainParams(v);
		
		Variable C=new Variable("Cp",v);
		C.setCommentAndUnit("coefficient C' of elliptic equation");
		C.setValue(1f);
		
		return C;
	}
	
	private Variable cAXY(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Aa",v);
		A.setCommentAndUnit("coefficient A of elliptical equation");
		A.setValue(1f);
		
		return A;
	}
	
	private Variable cAYZ(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Aa",v);
		A.setCommentAndUnit("coefficient A of elliptical equation");
		A.setValue(1f);
		
		return A;
	}
	
	private Variable cBXY(Variable v){
		Variable B=new Variable("Bb",v);
		B.setCommentAndUnit("coefficient B of elliptic equation");
		return B;
	}
	
	private Variable cBYZ(Variable v){
		Variable B=new Variable("Bb",v);
		B.setCommentAndUnit("coefficient B of elliptic equation");
		return B;
	}
	
	private Variable cCXY(Variable v){
		assignSubDomainParams(v);
		
		Variable C=new Variable("Cc",v);
		C.setCommentAndUnit("coefficient C of elliptic equation");
		C.setValue(1f);
		
		return C;
	}
	
	private Variable cCYZ(Variable v){
		assignSubDomainParams(v);
		
		Variable C=new Variable("Cc",v);
		C.setCommentAndUnit("coefficient C of elliptic equation");
		C.setValue(1f);
		
		return C;
	}
	
	
	/** test
	public static void main(String[] arg){
		ConcurrentUtil.initDefaultExecutor(1);
		
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] vs=df.getVariables(new Range("t(1,1);z(36,36)",dd),"u","v");
		
		Variable u=vs[0];
		Variable v=vs[1];
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		
		DynamicMethodsInSC dm=new DynamicMethodsInSC(ssm); dm.setBCofX(BoundaryCondition.Periodic);
		PoissonEquationInSC pe=new PoissonEquationInSC(ssm);
		GlobalVelocityFieldInSC gvf=new GlobalVelocityFieldInSC(ssm);
		
		Variable div=dm.c2DDivergence(u,v);
		Variable vor=dm.c2DVorticity(u,v);
		
		Variable sf1=pe.invertingBySOR(vor); sf1.setName("sf1");
		Variable sf2=pe.invertingBySH(vor); sf2.setName("sf2");
		
		Variable pf1=pe.invertingBySOR(div); pf1.setName("pf1");
		Variable pf2=pe.invertingBySH(div); pf2.setName("pf2");
		
		Variable[] velR1=gvf.cRotationalVelocity(sf1); velR1[0].setName("ur1"); velR1[1].setName("vr1");
		Variable[] velR2=gvf.cRotationalVelocity(sf2); velR2[0].setName("ur2"); velR2[1].setName("vr2");
		
		Variable[] velD1=dm.c2DGradient(pf1); velD1[0].setName("ud1"); velD1[1].setName("vd1");
		Variable[] velD2=dm.c2DGradient(pf2); velD2[0].setName("ud2"); velD2[1].setName("vd2");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/sf.dat");
		dw.writeData(dd,u,v,div,vor,sf1,sf2,pf1,pf2,velR1[0],velR1[1],velR2[0],velR2[1],velD1[0],velD1[1],velD2[0],velD2[1]);
		dw.closeFile();
		
		ConcurrentUtil.shutdown();
	}*/
}
