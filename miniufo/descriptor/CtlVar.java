/**
 * @(#)CtlVar.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import java.util.Scanner;
import java.util.regex.Matcher;


/**
 * used to descripe the variable in the CtlDescriptor file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class CtlVar extends Var{
	//
	private int index;		// start with 0
	
	private long start_pos;
	
	private String storage=null;
	
	
	/**
	 * constructor
	 *
	 * @param	oneline		Variable in String form
	 */
	public CtlVar(String oneline){
		try(Scanner sn=new Scanner(oneline)){
			vname=sn.next();
			zcount=Integer.parseInt(sn.next());	zcount=(zcount==0?1:zcount);
			storage=sn.next();
			
			if(!"99".equals(storage)&&!"-1,20".equals(storage))
				throw new IllegalArgumentException("only '99' or '-1,20' supported");
			
			String commentAndUnit=sn.nextLine();
			Matcher m=unitPtn.matcher(commentAndUnit);
			
			if(m.find()){
				String unitInBracket=m.group();
				unit=unitInBracket.substring(1,unitInBracket.length()-1);
				comment=m.replaceAll("").trim();
				
			}else comment=commentAndUnit.trim();
		}
	}
	
	
	/*** getor and setor ***/
	public  int getIndex(){ return index;}
	
	public long getStartPosition(){ return start_pos;}
	
	public String getStorageType(){ return storage  ;}
	
	public void setIndex(int index){ this.index=index;}
	
	public void setStartPosition(long start_pos){ this.start_pos=start_pos;}
}
