/**
 * @(#)Shaded2D.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.visualize;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.image.BufferedImage;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.SwingConstants;
import miniufo.basic.ArrayUtil;


/**
 * Shaded view of an array of 2-D discrete data
 */
public final class Shaded2D extends JFrame{
	//
	private int mul=1;
	private int x  =1;
	private int y  =1;
	
	private int[][] colors=null;
	
	private float[][] data=null;
	
	private BufferedImage bi=null;
	
	private static final long serialVersionUID=4412909891321313417L;
	
	public static final int[] GRADIENTCOLORS=gradientColors();
	
	
	/**
	 * constructor
	 */
	public Shaded2D(int multiple,float[][] data){
		this.data=data;
		
		y=data.length;	x=data[0].length;	mul=multiple;
		
		dataToColors();
		
		getContentPane().setLayout(new BorderLayout());
		
		bi=new BufferedImage(x*mul,y*mul,BufferedImage.TYPE_4BYTE_ABGR);
		
		JLabel label=new JLabel(new ImageIcon(bi),SwingConstants.CENTER);
		getContentPane().add(label,BorderLayout.CENTER);
		
		setSize(x*mul+30,y*mul+50);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		
		setVisible(true);
	}
	
	
	/*** getor and setor ***/
	public float[][] getData(){ return data;}
	
	public BufferedImage getBufferedImage(){ return bi;}
	
	
	public void repaint(){
		dataToColors();
		
		for(int j=0,J=y*mul;j<J;j++){
			int jtag=j/mul;
			
			for(int i=0,I=x*mul;i<I;i++)
			bi.setRGB(i,j,colors[jtag][i/mul]);
			
		}
		
		super.repaint();
	}
	
	public void smoothColor5(){
		int[][] buf=colors.clone();
		
		for(int j=1,J=y-1;j<J;j++)
		for(int i=1,I=x-1;i<I;i++){
			int cc=colors[j][i];
			int rc=(cc>>16)&0xFF;
			int gc=(cc>>8 )&0xFF;
			int bc=(cc    )&0xFF;
			
			int cl=colors[j][i-1];
			int rl=(cl>>16)&0xFF;
			int gl=(cl>>8 )&0xFF;
			int bl=(cl    )&0xFF;
			
			int cr=colors[j][i+1];
			int rr=(cr>>16)&0xFF;
			int gr=(cr>>8 )&0xFF;
			int br=(cr    )&0xFF;
			
			int cu=colors[j+1][i];
			int ru=(cu>>16)&0xFF;
			int gu=(cu>>8 )&0xFF;
			int bu=(cu    )&0xFF;
			
			int cd=colors[j-1][i];
			int rd=(cd>>16)&0xFF;
			int gd=(cd>>8 )&0xFF;
			int bd=(cd    )&0xFF;
			
			int mr=(int)((rc*4+(rl+rr+ru+rd)+0.5f)/8);
			int mg=(int)((gc*4+(gl+gr+gu+gd)+0.5f)/8);
			int mb=(int)((bc*4+(bl+br+bu+bd)+0.5f)/8);
			
			buf[j][i]=getRBGColor(mr,mg,mb);
		}
		
		for(int j=0;j<y;j++) System.arraycopy(buf[j],0,colors[j],0,x);
	}
	
	public void smoothData5(){
		float[][] buf=data.clone();
		
		for(int j=1,J=y-1;j<J;j++)
		for(int i=1,I=x-1;i<I;i++)
		buf[j][i]=data[j][i]/2+(data[j+1][i]+data[j-1][i]+data[j][i+1]+data[j][i-1])/8;
		
		for(int j=0;j<y;j++) System.arraycopy(buf[j],0,data[j],0,x);
	}
	
	
	private static int getRBGColor(int r,int g,int b){
		return ((0xFF)<<24)|((r&0xFF)<<16)|((g&0xFF)<<8)|((b&0xFF));
	}
	
	
	private void dataToColors(){
		colors=new int[y][x];
		
		float max=ArrayUtil.getMax(data);
		float min=ArrayUtil.getMin(data);
		
		float step=(max-min)/(GRADIENTCOLORS.length-1);
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			int tag=(int)((max-data[j][i])/step);
			
			colors[y-j-1][i]=GRADIENTCOLORS[tag];
		}
	}
	
	private static int[] gradientColors(){
		int[] cs=new int[239];
		
		int r=255,g=0,b=140;
		
		for(int i=  0;i< 10;i++) cs[i]=new Color(r,g,b-=10).getRGB();
		for(int i= 10;i< 18;i++) cs[i]=new Color(r,g,b-=5).getRGB();
		for(int i= 18;i< 69;i++) cs[i]=new Color(r,g+=5,b).getRGB();
		for(int i= 69;i<120;i++) cs[i]=new Color(r-=5,g,b).getRGB();
		for(int i=120;i<171;i++) cs[i]=new Color(r,g,b+=5).getRGB();
		for(int i=171;i<222;i++) cs[i]=new Color(r,g-=5,b).getRGB();
		for(int i=222;i<229;i++) cs[i]=new Color(r+=5,g,b).getRGB();
		for(int i=229;i<239;i++) cs[i]=new Color(r+=10,g,b).getRGB();
		
		return cs;
	}
	
	
	/** test
	public static void main(String args[]){
		Variable u=DiagnosisFactory.getVariables("d:/Data/sst.mnmean.nc","t(1,1)","sst")[0];
		
		Shaded2D sd=new Shaded2D(5,u.getData()[0][0]);
		sd.repaint();
		
		for(int i=0;i<10;i++){
			try{ Thread.sleep(2000L);}
			catch(Exception e){ e.printStackTrace(); System.exit(0);}
			
			sd.smoothData5();
			sd.repaint();
		}
    }*/
}
