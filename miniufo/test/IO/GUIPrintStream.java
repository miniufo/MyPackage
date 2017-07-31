/**
 * @(#)GUIPrintStream.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.IO;

import java.io.PrintStream;
import javax.swing.SwingUtilities;
import javax.swing.text.JTextComponent;


/**
 * redirect System.io to a GUI text component
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class GUIPrintStream extends PrintStream{
	//
	private int threshold=1024;
	
	private JTextComponent dest=null;
	private StringBuilder buffer=null;
	
	/**
     * constructor
     *
     * @param	component	a GUI component for redirection
     */
	public GUIPrintStream(JTextComponent component){
		super(System.out);
		
		dest=component;
		buffer=new StringBuilder();
	}
	
	
	/*** getor and setor ***/
	public int getSize(){ return threshold;}
	
	public void setSize(int size){ threshold=size;}
	
	
	/** 
	 * @Override
	 */
	public void write(byte[] buf,int off,int len){
		final String message=new String(buf,off,len);
		
		SwingUtilities.invokeLater(new Runnable(){
			public void run(){
				buffer.append(message);
				
				int length=buffer.length();
				int start=0;
				
				while(buffer.length()-start>threshold){
					start=buffer.indexOf("\n",start+1);
				}
				
				if(length>threshold) buffer.delete(0,start+1);
				
				dest.setText(buffer.toString());
				dest.setCaretPosition(dest.getCaretPosition());
			}
		});
	}
}
