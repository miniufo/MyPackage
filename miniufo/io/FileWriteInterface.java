/**
 * @(#)FileWriteInterface.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import java.io.File;
import java.io.IOException;


/**
 * used to interact with the user when writting a file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class FileWriteInterface{
	//
	private File   file   =null;
	private String newname=null;	// store the new name of the file if flag is 2
	private Solution  flag=Solution.OVERWRITE;
	
	public enum Solution{OVERWRITE,SKIP,RENAME,APPEND}
	
	
	/**
     * constructor
     *
     * @param	path	a string specified the file path
     */
	public FileWriteInterface(String path){
		file=new File(path);
		confirm();
	}
	
	
	/*** getor and setor ***/
	public Solution getFlag(){ return flag;}
	
	public File getFile(){ return file;}
	
	public String getNewName(){ return newname;}
	
	public String getParent(){ return file.getParent();}
	
	
	/**
     * interaction when the file is existed
     */
	private synchronized void confirm(){
		if(file.exists()){
			try{
			    byte[] flag=new byte[1];
				System.out.println("Warning:");
				System.out.println(file.getAbsolutePath()+" already exists!\nWould you like to continue?");
				System.in.skip(System.in.available());
				
			bb:	while(true){
					System.out.println("Input A(a) to append the file");
					System.out.println("Input O(o) to overwrite the existing file");
					System.out.println("Input R(r) to rename the file to be written(30chars)");
					System.out.println("Input S(s) to skip writing this file");
					System.out.println("Input T(t) to terminate the program without writing");
					System.in.skip(System.in.available());
					System.in.read(flag);
					
					if(System.in.available()==2||System.in.available()==1){
						switch(flag[0]){
						case 0x4F:
						case 0x6F:
							this.flag=Solution.OVERWRITE;	break bb;
							
						case 0x53:
						case 0x73:
							this.flag=Solution.SKIP     ;	break bb;
							
						case 0x41:
						case 0x61:
							this.flag=Solution.APPEND   ;	break bb;
							
						case 0x52:
						case 0x72:
							byte[] byt=new byte[64];
							
							System.in.skip(System.in.available());
							System.out.println("Please enter the new file name:(30 chars)");
							System.in.read(byt);
							newname=new String(byt);
							newname=newname.trim();
							
							this.flag=Solution.RENAME   ;	 break bb;
							
						case 0x54:
						case 0x74:
							System.out.println("The program terminates");
							System.exit(0);
						}
					}
				}
		    
		    }catch(IOException ex){ ex.printStackTrace(); System.exit(0);}
		}
	}
	
	
	/** test
	public static void main(String[] args){
		try{
			FileWriteInterface fwi=new FileWriteInterface("d:/a.dat");
			System.out.println(fwi.getParent()+fwi.getNewName());
			System.out.println(fwi.getParent()+"\\"+fwi.getNewName());
			System.out.println(fwi.getParent()+"/"+fwi.getNewName());
			
	    }catch(Exception e){ e.printStackTrace();}
	}*/
}
