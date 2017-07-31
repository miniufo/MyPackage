/**
 * @(#)Utility.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;


/**
 * utility implemented in java
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Utility{
	
	/**
	 * prevent from instantiate
     */
	private Utility(){}
	
	
	/**
	 * run shell/cmd command and print the standard output
	 * 
	 * @param	command		command used to run in shell/cmd
     */
	public static void runShell(String command){
		try{
			BufferedReader br=new BufferedReader(
				new InputStreamReader(
					Runtime.getRuntime().exec(command).getInputStream()
				)
			);
			String str=null;
			
			System.out.println();
			while((str=br.readLine())!=null) System.out.println(str);
			System.out.println();
			
			br.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	/**
	 * run shell/cmd command and return the standard output
	 * 
	 * @param	command		command used to run in shell/cmd
     */
	public static String getShellOutput(String command){
		StringBuilder sb=new StringBuilder();
		
		try{
			BufferedReader br=new BufferedReader(
				new InputStreamReader(
					Runtime.getRuntime().exec(command).getInputStream()
				)
			);
			
			String str=null;
			
			while((str=br.readLine())!=null) sb.append(str+"\n");
			
			br.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		return sb.toString();
	}
	
	
	/** test
	public static void main(String[] arg){
		System.out.println(getShellOutput("ls"));
	}*/
}
