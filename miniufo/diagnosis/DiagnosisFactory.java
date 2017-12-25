/**
 * @(#)DiagnosisFactory.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.diagnosis;

import java.io.File;
import java.io.IOException;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import miniufo.io.DataRead;
import miniufo.io.DataIOFactory;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.NetCDFDescriptor;


/**
 * Wrapper of a DataDescriptor and associated variables.
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DiagnosisFactory{
	//
	private static boolean print=true;
	
	private DataDescriptor dd=null;
	
	public static final DiagnosisFactory DF10=
	DiagnosisFactory.parseContent(
		"dset ^Model10\n"+
		"title 10-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 36 linear   0 10\n"+
		"ydef 19 linear -90 10\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DF5=
	DiagnosisFactory.parseContent(
		"dset ^Model5\n"+
		"title 10-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 72 linear   0 5\n"+
		"ydef 37 linear -90 5\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DF2P5=
	DiagnosisFactory.parseContent(
		"dset ^Model2P5\n"+
		"title 2.5-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 144 linear   0 2.5\n"+
		"ydef  73 linear -90 2.5\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DF2=
	DiagnosisFactory.parseContent(
		"dset ^Model2\n"+
		"title 2-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 180 linear   0 2\n"+
		"ydef  90 linear -89 2\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DF1P5=
	DiagnosisFactory.parseContent(
		"dset ^Model1P5\n"+
		"title 1.5-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 240 linear   0 1.5\n"+
		"ydef 121 linear -90 1.5\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DF1=
	DiagnosisFactory.parseContent(
		"dset ^Model1\n"+
		"title 1-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 360 linear   0 1\n"+
		"ydef 181 linear -90 1\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DFHalf=
	DiagnosisFactory.parseContent(
		"dset ^ModelHalf\n"+
		"title 0.5-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 720 linear   0 0.5\n"+
		"ydef 361 linear -90 0.5\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DFThird=
			DiagnosisFactory.parseContent(
				"dset ^ModelHalf\n"+
				"title 1/3-deg resolution model\n"+
				"undef -9.99e8\n"+
				"xdef 1080 linear   0 0.3333333333\n"+
				"ydef  541 linear -90 0.3333333333\n"+
				"zdef   1 levels 0 1\n"+
				"tdef   1 linear 00z01Jan2000 1dy\n"+
				"vars 1\n"+
				"test 1 99 test variable\n"+
				"endvars\n"
			);
	
	public static final DiagnosisFactory DFQuater=
	DiagnosisFactory.parseContent(
		"dset ^ModelQuater\n"+
		"title 0.25-deg resolution model\n"+
		"undef -9.99e8\n"+
		"xdef 1440 linear   0 0.25\n"+
		"ydef  721 linear -90 0.25\n"+
		"zdef    1 levels 0 1\n"+
		"tdef    1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DFT42=
	DiagnosisFactory.parseContent(
		"dset ^ModelHalf\n"+
		"title T42 gridded model\n"+
		"undef -9.99e8\n"+
		"xdef 128 linear 0 2.8125\n"+
		"ydef  64 levels\n"+
		 "-87.86385\n"+
		 "-85.09651\n"+
		 "-82.31290\n"+
		 "-79.52561\n"+
		 "-76.73689\n"+
		 "-73.94751\n"+
		 "-71.15775\n"+
		 "-68.36776\n"+
		 "-65.57761\n"+
		 "-62.78735\n"+
		 "-59.99702\n"+
		 "-57.20663\n"+
		 "-54.41620\n"+
		 "-51.62574\n"+
		 "-48.83524\n"+
		 "-46.04473\n"+
		 "-43.25419\n"+
		 "-40.46365\n"+
		 "-37.67309\n"+
		 "-34.88252\n"+
		 "-32.09195\n"+
		 "-29.30136\n"+
		 "-26.51077\n"+
		 "-23.72017\n"+
		 "-20.92958\n"+
		 "-18.13897\n"+
		 "-15.34837\n"+
		 "-12.55776\n"+
		 "-9.767150\n"+
		 "-6.976533\n"+
		 "-4.185923\n"+
		 "-1.395312\n"+
		 " 1.395312\n"+
		 " 4.185923\n"+
		 " 6.976533\n"+
		 " 9.767150\n"+
		 " 12.55776\n"+
		 " 15.34837\n"+
		 " 18.13897\n"+
		 " 20.92958\n"+
		 " 23.72017\n"+
		 " 26.51077\n"+
		 " 29.30136\n"+
		 " 32.09195\n"+
		 " 34.88252\n"+
		 " 37.67309\n"+
		 " 40.46365\n"+
		 " 43.25419\n"+
		 " 46.04473\n"+
		 " 48.83524\n"+
		 " 51.62574\n"+
		 " 54.41620\n"+
		 " 57.20663\n"+
		 " 59.99702\n"+
		 " 62.78735\n"+
		 " 65.57761\n"+
		 " 68.36776\n"+
		 " 71.15775\n"+
		 " 73.94751\n"+
		 " 76.73689\n"+
		 " 79.52561\n"+
		 " 82.31290\n"+
		 " 85.09651\n"+
		 " 87.86385\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DFT106=
	DiagnosisFactory.parseContent(
		"dset ^ModelHalf\n"+
		"title T106 gridded model\n"+
		"undef -9.99e8\n"+
		"xdef 320 linear 0 1.125\n"+
		"ydef 160 levels\n"+
		 "-89.14157\n"+
		 "-88.02940\n"+
		 "-86.91079\n"+
		 "-85.79062\n"+
		 "-84.66994\n"+
		 "-83.54895\n"+
		 "-82.42783\n"+
		 "-81.30659\n"+
		 "-80.18530\n"+
		 "-79.06399\n"+
		 "-77.94262\n"+
		 "-76.82124\n"+
		 "-75.69984\n"+
		 "-74.57844\n"+
		 "-73.45701\n"+
		 "-72.33559\n"+
		 "-71.21414\n"+
		 "-70.09268\n"+
		 "-68.97125\n"+
		 "-67.84978\n"+
		 "-66.72832\n"+
		 "-65.60687\n"+
		 "-64.48540\n"+
		 "-63.36394\n"+
		 "-62.24246\n"+
		 "-61.12099\n"+
		 "-59.99952\n"+
		 "-58.87804\n"+
		 "-57.75657\n"+
		 "-56.63509\n"+
		 "-55.51361\n"+
		 "-54.39214\n"+
		 "-53.27066\n"+
		 "-52.14918\n"+
		 "-51.02769\n"+
		 "-49.90621\n"+
		 "-48.78473\n"+
		 "-47.66325\n"+
		 "-46.54177\n"+
		 "-45.42028\n"+
		 "-44.29879\n"+
		 "-43.17731\n"+
		 "-42.05583\n"+
		 "-40.93434\n"+
		 "-39.81285\n"+
		 "-38.69137\n"+
		 "-37.56989\n"+
		 "-36.44839\n"+
		 "-35.32691\n"+
		 "-34.20541\n"+
		 "-33.08393\n"+
		 "-31.96244\n"+
		 "-30.84095\n"+
		 "-29.71946\n"+
		 "-28.59798\n"+
		 "-27.47649\n"+
		 "-26.35500\n"+
		 "-25.23351\n"+
		 "-24.11202\n"+
		 "-22.99053\n"+
		 "-21.86904\n"+
		 "-20.74756\n"+
		 "-19.62608\n"+
		 "-18.50458\n"+
		 "-17.38309\n"+
		 "-16.26160\n"+
		 "-15.14011\n"+
		 "-14.01863\n"+
		 "-12.89714\n"+
		 "-11.77564\n"+
		 "-10.65415\n"+
		 "-9.532662\n"+
		 "-8.411180\n"+
		 "-7.289683\n"+
		 "-6.168194\n"+
		 "-5.046710\n"+
		 "-3.925221\n"+
		 "-2.803717\n"+
		 "-1.682235\n"+
		 " 0.5607517\n"+
		 " 0.5607517\n"+
		 " 1.682235\n"+
		 " 2.803717\n"+
		 " 3.925221\n"+
		 " 5.046710\n"+
		 " 6.168194\n"+
		 " 7.289683\n"+
		 " 8.411180\n"+
		 " 9.532662\n"+
		 " 10.65415\n"+
		 " 11.77564\n"+
		 " 12.89714\n"+
		 " 14.01863\n"+
		 " 15.14011\n"+
		 " 16.26160\n"+
		 " 17.38309\n"+
		 " 18.50458\n"+
		 " 19.62608\n"+
		 " 20.74756\n"+
		 " 21.86904\n"+
		 " 22.99053\n"+
		 " 24.11202\n"+
		 " 25.23351\n"+
		 " 26.35500\n"+
		 " 27.47649\n"+
		 " 28.59798\n"+
		 " 29.71946\n"+
		 " 30.84095\n"+
		 " 31.96244\n"+
		 " 33.08393\n"+
		 " 34.20541\n"+
		 " 35.32691\n"+
		 " 36.44839\n"+
		 " 37.56989\n"+
		 " 38.69137\n"+
		 " 39.81285\n"+
		 " 40.93434\n"+
		 " 42.05583\n"+
		 " 43.17731\n"+
		 " 44.29879\n"+
		 " 45.42028\n"+
		 " 46.54177\n"+
		 " 47.66325\n"+
		 " 48.78473\n"+
		 " 49.90621\n"+
		 " 51.02769\n"+
		 " 52.14918\n"+
		 " 53.27066\n"+
		 " 54.39214\n"+
		 " 55.51361\n"+
		 " 56.63509\n"+
		 " 57.75657\n"+
		 " 58.87804\n"+
		 " 59.99952\n"+
		 " 61.12099\n"+
		 " 62.24246\n"+
		 " 63.36394\n"+
		 " 64.48540\n"+
		 " 65.60687\n"+
		 " 66.72832\n"+
		 " 67.84978\n"+
		 " 68.97125\n"+
		 " 70.09268\n"+
		 " 71.21414\n"+
		 " 72.33559\n"+
		 " 73.45701\n"+
		 " 74.57844\n"+
		 " 75.69984\n"+
		 " 76.82124\n"+
		 " 77.94262\n"+
		 " 79.06399\n"+
		 " 80.18530\n"+
		 " 81.30659\n"+
		 " 82.42783\n"+
		 " 83.54895\n"+
		 " 84.66994\n"+
		 " 85.79062\n"+
		 " 86.91079\n"+
		 " 88.02940\n"+
		 " 89.14157\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	public static final DiagnosisFactory DFT170=
	DiagnosisFactory.parseContent(
		"dset ^ModelHalf\n"+
		"title T170 gridded model\n"+
		"undef -9.99e8\n"+
		"xdef 512 linear 0 0.703125\n"+
		"ydef 256 levels\n"+
		  "-89.46294\n"+
		  "-88.76694\n"+
		  "-88.06700\n"+
		  "-87.36604\n"+
		  "-86.66481\n"+
		  "-85.96337\n"+
		  "-85.26186\n"+
		  "-84.56026\n"+
		  "-83.85863\n"+
		  "-83.15699\n"+
		  "-82.45533\n"+
		  "-81.75364\n"+
		  "-81.05194\n"+
		  "-80.35023\n"+
		  "-79.64852\n"+
		  "-78.94681\n"+
		  "-78.24510\n"+
		  "-77.54337\n"+
		  "-76.84164\n"+
		  "-76.13990\n"+
		  "-75.43818\n"+
		  "-74.73644\n"+
		  "-74.03470\n"+
		  "-73.33298\n"+
		  "-72.63124\n"+
		  "-71.92949\n"+
		  "-71.22775\n"+
		  "-70.52600\n"+
		  "-69.82426\n"+
		  "-69.12252\n"+
		  "-68.42078\n"+
		  "-67.71903\n"+
		  "-67.01729\n"+
		  "-66.31554\n"+
		  "-65.61380\n"+
		  "-64.91206\n"+
		  "-64.21030\n"+
		  "-63.50855\n"+
		  "-62.80680\n"+
		  "-62.10506\n"+
		  "-61.40331\n"+
		  "-60.70155\n"+
		  "-59.99981\n"+
		  "-59.29806\n"+
		  "-58.59632\n"+
		  "-57.89456\n"+
		  "-57.19282\n"+
		  "-56.49107\n"+
		  "-55.78931\n"+
		  "-55.08756\n"+
		  "-54.38581\n"+
		  "-53.68406\n"+
		  "-52.98231\n"+
		  "-52.28057\n"+
		  "-51.57881\n"+
		  "-50.87706\n"+
		  "-50.17530\n"+
		  "-49.47356\n"+
		  "-48.77181\n"+
		  "-48.07006\n"+
		  "-47.36831\n"+
		  "-46.66655\n"+
		  "-45.96480\n"+
		  "-45.26305\n"+
		  "-44.56130\n"+
		  "-43.85955\n"+
		  "-43.15779\n"+
		  "-42.45604\n"+
		  "-41.75429\n"+
		  "-41.05254\n"+
		  "-40.35078\n"+
		  "-39.64904\n"+
		  "-38.94728\n"+
		  "-38.24553\n"+
		  "-37.54378\n"+
		  "-36.84203\n"+
		  "-36.14027\n"+
		  "-35.43852\n"+
		  "-34.73676\n"+
		  "-34.03502\n"+
		  "-33.33326\n"+
		  "-32.63151\n"+
		  "-31.92976\n"+
		  "-31.22801\n"+
		  "-30.52625\n"+
		  "-29.82450\n"+
		  "-29.12274\n"+
		  "-28.42100\n"+
		  "-27.71924\n"+
		  "-27.01749\n"+
		  "-26.31573\n"+
		  "-25.61399\n"+
		  "-24.91223\n"+
		  "-24.21048\n"+
		  "-23.50872\n"+
		  "-22.80698\n"+
		  "-22.10522\n"+
		  "-21.40347\n"+
		  "-20.70172\n"+
		  "-19.99995\n"+
		  "-19.29821\n"+
		  "-18.59645\n"+
		  "-17.89471\n"+
		  "-17.19295\n"+
		  "-16.49120\n"+
		  "-15.78944\n"+
		  "-15.08769\n"+
		  "-14.38594\n"+
		  "-13.68418\n"+
		  "-12.98243\n"+
		  "-12.28068\n"+
		  "-11.57893\n"+
		  "-10.87717\n"+
		  "-10.17542\n"+
		  "-9.473663\n"+
		  "-8.771915\n"+
		  "-8.070162\n"+
		  "-7.368414\n"+
		  "-6.666654\n"+
		  "-5.964906\n"+
		  "-5.263159\n"+
		  "-4.561391\n"+
		  "-3.859637\n"+
		  "-3.157890\n"+
		  "-2.456136\n"+
		  "-1.754382\n"+
		  "-1.052635\n"+
		  "-0.3508805\n"+
		  " 0.3508805\n"+
		  "1.052635\n"+
		  "1.754382\n"+
		  "2.456136\n"+
		  "3.157890\n"+
		  "3.859637\n"+
		  "4.561391\n"+
		  "5.263159\n"+
		  "5.964906\n"+
		  "6.666654\n"+
		  "7.368414\n"+
		  "8.070162\n"+
		  "8.771915\n"+
		  "9.473663\n"+
		  "10.17542\n"+
		  "10.87717\n"+
		  "11.57893\n"+
		  "12.28068\n"+
		  "12.98243\n"+
		  "13.68418\n"+
		  "14.38594\n"+
		  "15.08769\n"+
		  "15.78944\n"+
		  "16.49120\n"+
		  "17.19295\n"+
		  "17.89471\n"+
		  "18.59645\n"+
		  "19.29821\n"+
		  "19.99995\n"+
		  "20.70172\n"+
		  "21.40347\n"+
		  "22.10522\n"+
		  "22.80698\n"+
		  "23.50872\n"+
		  "24.21048\n"+
		  "24.91223\n"+
		  "25.61399\n"+
		  "26.31573\n"+
		  "27.01749\n"+
		  "27.71924\n"+
		  "28.42100\n"+
		  "29.12274\n"+
		  "29.82450\n"+
		  "30.52625\n"+
		  "31.22801\n"+
		  "31.92976\n"+
		  "32.63151\n"+
		  "33.33326\n"+
		  "34.03502\n"+
		  "34.73676\n"+
		  "35.43852\n"+
		  "36.14027\n"+
		  "36.84203\n"+
		  "37.54378\n"+
		  "38.24553\n"+
		  "38.94728\n"+
		  "39.64904\n"+
		  "40.35078\n"+
		  "41.05254\n"+
		  "41.75429\n"+
		  "42.45604\n"+
		  "43.15779\n"+
		  "43.85955\n"+
		  "44.56130\n"+
		  "45.26305\n"+
		  "45.96480\n"+
		  "46.66655\n"+
		  "47.36831\n"+
		  "48.07006\n"+
		  "48.77181\n"+
		  "49.47356\n"+
		  "50.17530\n"+
		  "50.87706\n"+
		  "51.57881\n"+
		  "52.28057\n"+
		  "52.98231\n"+
		  "53.68406\n"+
		  "54.38581\n"+
		  "55.08756\n"+
		  "55.78931\n"+
		  "56.49107\n"+
		  "57.19282\n"+
		  "57.89456\n"+
		  "58.59632\n"+
		  "59.29806\n"+
		  "59.99981\n"+
		  "60.70155\n"+
		  "61.40331\n"+
		  "62.10506\n"+
		  "62.80680\n"+
		  "63.50855\n"+
		  "64.21030\n"+
		  "64.91206\n"+
		  "65.61380\n"+
		  "66.31554\n"+
		  "67.01729\n"+
		  "67.71903\n"+
		  "68.42078\n"+
		  "69.12252\n"+
		  "69.82426\n"+
		  "70.52600\n"+
		  "71.22775\n"+
		  "71.92949\n"+
		  "72.63124\n"+
		  "73.33298\n"+
		  "74.03470\n"+
		  "74.73644\n"+
		  "75.43818\n"+
		  "76.13990\n"+
		  "76.84164\n"+
		  "77.54337\n"+
		  "78.24510\n"+
		  "78.94681\n"+
		  "79.64852\n"+
		  "80.35023\n"+
		  "81.05194\n"+
		  "81.75364\n"+
		  "82.45533\n"+
		  "83.15699\n"+
		  "83.85863\n"+
		  "84.56026\n"+
		  "85.26186\n"+
		  "85.96337\n"+
		  "86.66481\n"+
		  "87.36604\n"+
		  "88.06700\n"+
		  "88.76694\n"+
		  "89.46294\n"+
		"zdef   1 levels 0 1\n"+
		"tdef   1 linear 00z01Jan2000 1dy\n"+
		"vars 1\n"+
		"test 1 99 test variable\n"+
		"endvars\n"
	);
	
	
	/**
	 * constructor
	 */
	private DiagnosisFactory(){}
	
	
	public static DiagnosisFactory parseFile(String filePath){
		DiagnosisFactory df=new DiagnosisFactory();
		
		if(filePath.endsWith(".ctl"))
			df.dd=new CtlDescriptor(new File(filePath));
		
		else if(filePath.endsWith(".csm"))
			df.dd=new CsmDescriptor(new File(filePath));
		
		else if(filePath.endsWith(".cts"))
			df.dd=new CtsDescriptor(new File(filePath));
		
		else if(filePath.endsWith(".nc" )||filePath.endsWith(".cdf")){
			try{df.dd=new NetCDFDescriptor(filePath);}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		else throw new IllegalArgumentException("unsupported DataDescriptor type");
		
		return df;
	}
	
	public static DiagnosisFactory parseContent(String content){
		DiagnosisFactory df=new DiagnosisFactory();
		
		if(content.toLowerCase().indexOf("\nvars ")!=-1)
			df.dd=new CtlDescriptor(content);
		
		else if(content.toLowerCase().indexOf("\ncoords")!=-1)
			df.dd=new CsmDescriptor(content);
		
		else if(content.toLowerCase().indexOf("\nf0")!=-1)
			df.dd=new CtsDescriptor(content);
			
		else throw new IllegalArgumentException("unsupported DataDescriptor type");
		
		return df;
	}
	
	
	/**
	 * Get variables from tstr to tend steps and return them as a stream of variable[].
	 * This may be useful for long-time dataset analysis to avoid out-of-memory error.
	 * 
	 * @param	tstr	start time step (from 1), inclusive
	 * @param	tend	end   time step (from 1), inclusive
	 * @param	v1		the first variable name
	 * @param	others	other names of variables
	 */
	public Stream<Variable[]> getVariablesTimeByTime(int tstr,int tend,String v1,String... others){
		if(tstr<1   ) throw new IllegalArgumentException("tstr should be >= 1");
		if(tstr>tend) throw new IllegalArgumentException("tstr should not be larger than tend");
		
		String[] vnames=new String[others.length+1];
		
		vnames[0]=v1;
		System.arraycopy(others,0,vnames,1,others.length);
		
		return IntStream.range(tstr,tend+1).sequential().mapToObj(l->{
			return getVariables(new Range("t("+l+","+l+")",dd),vnames);
		});
	}
	
	public Stream<Variable[]> getVariablesTimeByTime(String v1,String... others){
		return getVariablesTimeByTime(1,dd.getTCount(),v1,others);
	}
	
	/**
	 * Get a variable from tstr to tend steps and return it as a stream of variable.
	 * This may be useful for long-time dataset analysis to avoid out-of-memory error.
	 * 
	 * @param	tstr	start time step (from 1), inclusive
	 * @param	tend	end   time step (from 1), inclusive
	 * @param	vname	the variable name
	 */
	public Stream<Variable> getVariableTimeByTime(int tstr,int tend,String vname){
		if(tstr<0   ) throw new IllegalArgumentException("tstr should be > 0");
		if(tstr>tend) throw new IllegalArgumentException("tstr should not be larger than tend");
		
		return IntStream.range(tstr,tend+1).sequential().mapToObj(l->{
			return getVariables(new Range("t("+l+","+l+")",dd),vname)[0];
		});
	}
	
	public Stream<Variable> getVariableTimeByTime(String vname){
		return getVariableTimeByTime(1,dd.getTCount(),vname);
	}
	
	
	/*** getor and setor ***/
	public DataDescriptor getDataDescriptor(){ return dd;}
	
	public Variable[] getVariables(Range r,String... names){
		int count=names.length;
		
		Variable[] v=new Variable[count];
		
		for(int i=0;i<count;i++) v[i]=new Variable(names[i],r);
		
		DataRead dr=DataIOFactory.getDataRead(dd);
		dr.setPrinting(print);
		dr.readData(v);	dr.closeFile();
		
		return v;
	}
	
	public Variable[] getVariables(Range r,boolean tfirst,String... names){
		int count=names.length;
		
		Variable[] v=new Variable[count];
		
		for(int i=0;i<count;i++) v[i]=new Variable(names[i],tfirst,r);
		
		DataRead dr=DataIOFactory.getDataRead(dd);
		dr.setPrinting(print);
		dr.readData(v);	dr.closeFile();
		
		return v;
	}
	
	public void setPrinting(boolean print){ DiagnosisFactory.print=print;}
	
	
	/**
	 * static factory method
	 */
	public static DataDescriptor getDataDescriptor(String path){
		return parseFile(path).dd;
	}
	
	public static Variable[] getVariables(String path,String range,String... names){
		DiagnosisFactory df=parseFile(path);
		
		return df.getVariables(new Range(range,df.getDataDescriptor()),names);
	}
	
	public static Variable[] getVariables(String path,String range,boolean tfirst,String... names){
		DiagnosisFactory df=parseFile(path);
		
		return df.getVariables(new Range(range,df.getDataDescriptor()),tfirst,names);
	}
	
	
	/** test
	public static void main(String[] args){
		System.out.println(DFT170.dd);
	}*/
}
