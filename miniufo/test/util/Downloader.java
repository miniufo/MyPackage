/**
 * @(#)Downloader.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;

import java.io.BufferedInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.List;
import java.util.ArrayList;
import miniufo.io.FileWriteInterface;
import static miniufo.io.FileWriteInterface.Solution.*;


/**
 * mean process of the data
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class Downloader{
	//
	private static final int buffer_size=8096;	// buffer size
	
	private List<String> vdownload=null;		// URL vector
	private List<String> vfilelist=null;		// namelist
	
	
	/**
     * constructor
     */
	public Downloader(){
		vdownload=new ArrayList<String>();
		vfilelist=new ArrayList<String>();
	}
	
	
	/**
	 * clear items
	 */
	public void clearList(){
		vdownload.clear();
		vfilelist.clear();
	}
	
	/**
	 * add one item
	 * 
	 * @param	url		string
	 * @param	name	string
	 */
	public void addItem(String url,String name){
		vdownload.add(url);
		vfilelist.add(name);
	}
	
	/**
	 * download one by one according to the list
	 */
	public void downloadByList(){
		for(int i=0;i<vdownload.size();i++){
			String url =vdownload.get(i);
			String name=vfilelist.get(i);
			
			saveToFile(url, name);
		}
	}
	
	
	/**
	 * setting proxy
	 *
	 * @param	proxy		String
	 * @param	proxyport	String
	 */
	public void setProxyServer(String proxy,String proxyport){
		System.getProperties().put("proxyset", "true");
		System.getProperties().put("proxyhost",proxy);
		System.getProperties().put("proxyport",proxyport);
	}
	
	/**
	 * setting username and password
	 *
	 * @param uid String
	 * @param pwd String
	 */
	public void setAuthenticator(String uid,String pwd){
		//Authenticator.setdefault(new Myauthenticator(uid,pwd));
	}
	
	
	/**
	 * save url to local file
	 * 
	 * @param	desturl		String
	 * @param	filename	String
	 */
	public void saveToFile(String desturl,String filename){
		System.out.println("\nrequesting [" + desturl + "]...");
		
		try{
			URL url=new URL(desturl);
			HttpURLConnection httpurl=(HttpURLConnection)url.openConnection();
			
			httpurl.connect();
			
			BufferedInputStream bis=new BufferedInputStream(httpurl.getInputStream());
			
			boolean is_skip=true;	String name=filename;
			FileOutputStream fos=null;
			FileWriteInterface fwi=new FileWriteInterface(filename);
			
			if(fwi.getFlag()!=SKIP){
				is_skip=false;
				
				try{
					switch(fwi.getFlag()){
					case RENAME:
						name=fwi.getParent()+fwi.getNewName();
						fos=new FileOutputStream(name);
						break;
						
					case OVERWRITE:
						fos=new FileOutputStream(filename);
						break;
						
					case APPEND:
						fos=new FileOutputStream(filename,true);
						break;

					default:
						break;
					}
					
			    }catch(FileNotFoundException ex){ ex.printStackTrace(); System.exit(0);}
			}
			
			if(!is_skip){
				System.out.print("saving as [" + name + "]...\t");
				
				byte[] buf=new byte[buffer_size]; int size=0;
				while((size=bis.read(buf))!=-1) fos.write(buf,0,size);
				
				fos.close();	bis.close();	httpurl.disconnect();
			}
			
		}catch(IOException e){ System.out.print(e+"\ndo not ");}
		
		System.out.println("finished.");
	}
	
	
	/** test
	public static void main(String argv[]){
		miniufo.diagnosis.MDate md=new miniufo.diagnosis.MDate(2008,12,1,0,0,0);
		
		for(int i=220;i<=91*4;i++){
			String lt=String.valueOf(
				md.add(java.util.Calendar.HOUR_OF_DAY,i*6).getLongTime()
			);
			
			String yy=lt.substring(0,4);
			String tt=lt.substring(2,10);
			
			Downloader dl=new Downloader();
			//dl.setProxyServer("202.116.64.226","3128");
			dl.saveToFile(
				"http://agora.ex.nii.ac.jp/digital-typhoon/globe/color/"+yy+"/512x512/MTS1"+tt+".globe.1.jpg",
				"d:/animation/daily/MTS1"+tt+".globe.1.jpg"
			);
		}
	}*/
}
