/**
 * @(#)AssociatedLegendre.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.mathsphysics;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static miniufo.mathsphysics.MathsPhysicsUtil.twoFactorialDouble;


/**
 * associated Legendre class
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class AssociatedLegendre{
	//
	private int l;		// degree
	private int m;		// order
	
	
	/**
     * constructor
     *
     * @param	l	degree
     * @param	m	order
     */
	public AssociatedLegendre(int l,int m){
		if(m<0||m>l)
		throw new IllegalArgumentException("m and l should satisfy: 0 <= m < l");
		
		this.l=l; this.m=m;
	}
	
	public AssociatedLegendre(){ this(0,0);}
	
	
	/*** getor and setor ***/
	public int getL(){ return l;}
	
	public int getM(){ return m;}
	
	public void setL(int l){ this.l=l;}
	
	public void setM(int m){ this.m=m;}
	
	public void setLM(int l,int m){ this.l=l; this.m=m;}
	
	
	/*** recurrence formula ***/
	public float eval(float u){
		if(l==m&&m==0) return 1f;
		else{
			float last=0f,current=1f;
			
			if(m!=0) current=(float)(twoFactorialDouble(2*m-1)*pow((1-u*u),0.5*m));
			
			if(l!=m)
			for(int k=m;k<l;k++){
				float next=((2*k+1)*u*current-(k+m)*last)/(k+1-m);
				last=current;
				current=next;
			}
			
			return (m%2==1)?-current:current;
		}
	}
	
	public double evalDouble(double u){
		if(l==m&&m==0) return 1.0;
		else{
			double last=0.0,current=1.0;
			
			if(m!=0) current=twoFactorialDouble(2*m-1)*pow((1-u*u),0.5*m);
			
			if(l!=m)
			for(int k=m;k<l;k++){
				double next=((2*k+1)*u*current-(k+m)*last)/(k+1-m);
				last=current;
				current=next;
			}
			
			return (m%2==1)?-current:current;
		}
	}
	
	
	/**
     * normalized (by sqrt((2l+1)*(l-m)!/(l+m)!)) Legendre functions
     * (triangle truncation)
     *
     * @param	L	max degree
     * @param	u	Legendre function argument
     */
	public static float[][] legendreT(int L,float u){
		float[][] le=new float[L+1][];
		for(int i=0,I=L+1;i<I;i++) le[i]=new float[i+1];
		
		legendreT(L,u,le);
		
		return le;
	}
	
	public static void legendreT(int L,float u,float[][] le){
		L++;
		
		if(L<2) throw new IllegalArgumentException("L should be larger than 1");
		if(le.length!=L) throw new IllegalArgumentException("Invalid array dimension: L");
		for(int i=0;i<L;i++)
		if(le[i].length<=i) throw new IllegalArgumentException("Invalid array dimension: M");
		
		// first few
		le[0][0]=1;
		le[1][0]=(float)(sqrt(3)*u);
		le[1][1]=(float)(sqrt(1.5)*sqrt(1-u*u));
		le[2][1]=(float)sqrt(5.0)*u*le[1][1];
		
		// calculated P0N (3.30)
		for(int i=1,I=L-1;i<I;i++){
			double i2=i<<1;
			
			le[i+1][0]=(float)(
				sqrt((i2+1.0)*(i2+3.0))/(i+1.0)*u*le[i  ][0]-
				sqrt((i2+3.0)/(i2-1.0))/(i+1.0)*i*le[i-1][0]
			);
		}
		
		// calculated P1N (3.31)
		for(int i=2,I=L-1;i<I;i++){
			double i2=i<<1;
			
			le[i+1][1]=(float)(
				sqrt((i2+1.0)*(i2+3.0)/i/(i+2.0)) * u * le[i  ][1]-
				sqrt((i-1)*(i+1)*(i2+3)/i/(i+2)/(i2-1))*le[i-1][1]
			);
		}
		
		// calculated Pmn (3.29)
		for(int j=2;j<L;j++){
			double j2=j<<1;
			
			// Pm,m-1 == 0
			le[j][j]=(float)(
				sqrt((j2+1)*(j2-1)*(j2-3)/(j2-3)/j2/(j2-2))*le[j-2][j-2]-
				sqrt((j2+1)*(j2-1)/(j2-1)/j2/(j2-2))*u*le[j-1][j-2]
			);
			
			for(int i=j+1;i<L;i++){
				double i2=i<<1;
				
				le[i][j]=(float)(
					sqrt((i2+1)*(i+j-1)*(i+j-3)/(i2-3)/(i+j)/(i+j-2))*le[i-2][j-2]  -
					sqrt((i2+1)*(i+j-1)*(i-j+1)/(i2-1)/(i+j)/(i+j-2))*le[i-1][j-2]*u+
					sqrt((i2+1)*(i-j)/(i2-1)/(i+j))*le[i-1][j  ]*u
				);
			}
		}
	}
	
	
    /** test
    public static void main(String[] args){
    	try{
    		int N=5;
    		float[][] le=legendreT(N,0.2f);
    		
    		for(int i=N-1;i>=0;i--){
    			for(int j=0;j<=i;j++) System.out.print("P"+i+""+j+":"+le[i][j]+"\t");
    			System.out.println();
    		}
    		
    	}catch(Exception e){ e.printStackTrace();}
    }*/
}
