/**
 * @(#)InverseHGTInSC.java	1.0 2015.03.15
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.application.diagnosticModel;

import miniufo.concurrent.ConcurrentUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.test.application.basic.GlobalVelocityFieldInSC;
import miniufo.application.advanced.EllipticEqSORSolver;
import miniufo.application.advanced.EllipticEqSORSolver.DimCombination;
import miniufo.application.basic.DynamicMethodsInSC;
import miniufo.application.basic.SphericalHarmonicExpansion;
import miniufo.application.EllipticEquationInterface;
import miniufo.application.EquationInSphericalCoordinate;
import static java.lang.Math.cos;


/**
 * inverting wind field using geopotential in spheral coordinate
 *
 * @version 1.0, 2015.03.15
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class InverseHGTInSC extends EquationInSphericalCoordinate implements EllipticEquationInterface{
	
	/**
     * constructor
     *
     * @param	ssm		initialized by a spacial model in spheral coordinate
     */
	public InverseHGTInSC(SphericalSpatialModel ssm){
		super(ssm);
		
		if(!ssm.isLinearModel()) System.out.println("\nNot a equal-space model");
	}
	
	
	/**
     * inverting non-divergent flow using geopotential
     * 
     * @param	gpt		geopotential (m^2 s^-2)
     */
	public Variable invertingNondivergentFlow(Variable gpt){
		EllipticEqSORSolver ees=new EllipticEqSORSolver(sm);
		DynamicMethodsInSC dm=new DynamicMethodsInSC((SphericalSpatialModel)sm);
		
		ees.setBCofX(BCx); dm.setBCofX(BCx);
		ees.setBCofY(BCy); dm.setBCofY(BCy);
		
		Variable F=dm.cLaplacian(gpt);
		Variable S=new Variable("sf",gpt);
		S.setCommentAndUnit("geostrophic streamfunction (m^2/s)");
		
		ees.setDimCombination(DimCombination.XY);
		ees.setABC(cAPrime(F),cBPrime(F),cCPrime(F));
		ees.solve(S,F);
		
		return S;
	}
	
	
	/**
     * implement the methods in EllipticEquationInterface
     */
	public Variable cAPrime(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Ap",v); A.setCommentAndUnit("coefficient A' of elliptic equation");
		
		float[][][][] Adata=A.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[l][k][j][i]=f1[ystart-1+j]/lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=f1[ystart-1+j]/lcos[ystart-1+j];
		}
		
		return A;
	}
	
	public Variable cBPrime(Variable v){
		Variable B=new Variable("Bp",v);	B.setCommentAndUnit("coefficient B' of elliptic equation");
		return B;
	}
	
	public Variable cCPrime(Variable v){
		assignSubDomainParams(v);
		
		Variable C=new Variable("Cp",v); C.setCommentAndUnit("coefficient C' of elliptic equation");
		
		float[][][][] Cdata=C.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			for(int i=0;i<x;i++) Cdata[l][k][j][i]=
			(f1[ystart-2+j]+f1[ystart-1+j])/2f*(float)cos((ydef[ystart-2+j]+ydef[ystart-1+j])/2f);
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=1;j<y;j++)
			for(int i=0;i<x;i++) Cdata[k][j][i][l]=
			(f1[ystart-2+j]+f1[ystart-1+j])/2f*(float)cos((ydef[ystart-2+j]+ydef[ystart-1+j])/2f);
		}
		
		return C;
	}
	
	public Variable cA(Variable v){
		assignSubDomainParams(v);
		
		Variable A=new Variable("Aa",v); A.setCommentAndUnit("coefficient A of elliptical equation");
		
		float[][][][] Adata=A.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[l][k][j][i]=f1[ystart-1+j]*lcos[ystart-1+j];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) Adata[k][j][i][l]=f1[ystart-1+j]*lcos[ystart-1+j];
		}
		
		return A;
	}
	
	public Variable cB(Variable v){
		Variable B=new Variable("Bb",v);	B.setCommentAndUnit("coefficient B of elliptic equation");
		return B;
	}
	
	public Variable cC(Variable v){
		Variable C=cAPrime(v);
		C.setName("Cc"); C.setCommentAndUnit("coefficient C of elliptic equation");
		
		return C;
	}
	
	
	/** test*/
	public static void main(String[] args){
		ConcurrentUtil.initDefaultExecutor(1);
		
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] vs=df.getVariables(new Range("t(1,1);z(36,36)",dd),"u","v","h");
		
		Variable u=vs[0];
		Variable v=vs[1];
		Variable h=vs[2]; h.multiplyEq(9.8f);
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		
		DynamicMethodsInSC dm=new DynamicMethodsInSC(ssm);dm.setBCofX(BoundaryCondition.Periodic);
		InverseHGTInSC iHGT=new InverseHGTInSC(ssm);iHGT.setBCofX(BoundaryCondition.Periodic);iHGT.setBCofY(BoundaryCondition.Expanded);
		SphericalHarmonicExpansion she=new SphericalHarmonicExpansion(ssm); she.setM(dd.getYCount()-1);
		GlobalVelocityFieldInSC gvf=new GlobalVelocityFieldInSC(ssm);
		
		Variable vor=dm.c2DVorticity(u,v);
		
		Variable lpc=dm.cLaplacian(h);
		Variable frc=lpc.copy(); dm.weightLCos(lpc); frc.setName("frc");
		
		Variable sf1=she.solvePoissonEquation(vor); sf1.setName("sf1");
		Variable sf2=iHGT.invertingNondivergentFlow(h); sf2.setName("sf2");
		
		Variable[] velR1=gvf.cRotationalVelocity(sf1); velR1[0].setName("ur1"); velR1[1].setName("vr1");
		Variable[] velR2=gvf.cRotationalVelocity(sf2); velR2[0].setName("ur2"); velR2[1].setName("vr2");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/sf.dat");
		dw.writeData(dd,u,v,h,lpc,frc,vor,sf1,sf2,velR1[0],velR1[1],velR2[0],velR2[1]);
		dw.closeFile();
		
		ConcurrentUtil.shutdown();
	}
}
