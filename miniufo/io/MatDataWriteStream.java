/**
 * @(#)CtlDataWriteStream.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import miniufo.diagnosis.Variable;


/**
 * used to write the commonest ctl data file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class MatDataWriteStream{
	
	/**
     * prevent from initialization
     */
	private MatDataWriteStream(){}
	
	
	/**
     * write data from array to local .mat file
     *
     * @param	vs		variables
     */
	public static void writeData(String path,Variable... vs){
		List<MLArray> ls=new ArrayList<>();
		
		for(Variable v:vs){
			int t=v.getTCount();
			int z=v.getZCount();
			int y=v.getYCount();
			int x=v.getXCount();
			
			MLDouble mld=null;
			
			float[][][][] vdata=v.getData();
			
			if(v.isTFirst()){
				mld=new MLDouble(v.getName(),new int[]{y,x,z,t});
				
				int idx=0;
				
				for(int l=0;l<t;l++)
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++)
				for(int j=0;j<y;j++) mld.set((double)vdata[l][k][j][i],idx++);
				
			}else{
				mld=new MLDouble(v.getName(),new int[]{t,y,x,z});
				
				int idx=0;
				
				for(int k=0;k<z;k++)
				for(int i=0;i<x;i++)
				for(int j=0;j<y;j++)
				for(int l=0;l<t;l++) mld.set((double)vdata[k][j][i][l],idx++);
			}
			
			ls.add(mld);
		}
		
		try{ new MatFileWriter(path,ls);}
		catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	/** test
	public static void main(String[] arg){
		writeData("d:/test.mat",DiagnosisFactory.getVariables("d:/Plumb/flx.ctl","",true,"uwnd","vwnd"));
	}*/
}
