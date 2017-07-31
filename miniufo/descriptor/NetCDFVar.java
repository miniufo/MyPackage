/**
 * @(#)NetCDFVar.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.descriptor;

import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.Variable;


/**
 * used to descripe the variable in the NetCDFDescriptor file
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class NetCDFVar extends Var{
	
	/**
	 * variable in the NetCDFDescriptor file
	 *
	 * @param	v	NetCDF Variable
	 */
	public NetCDFVar(Variable v){
		for(Dimension dms:v.getDimensions()){
			if(dms.getFullName().equals("time" )){ tcount=dms.getLength();	continue;}
			if(dms.getFullName().equals("level")){ zcount=dms.getLength();	continue;}
			if(dms.getFullName().equals("lat"  )){ ycount=dms.getLength();	continue;}
			if(dms.getFullName().equals("lon"  )){ xcount=dms.getLength();	continue;}
		}
		
		Attribute attrb=null;	float factor=1,offset=0;
		attrb=v.findAttribute("scale_factor");	if(attrb!=null) factor=attrb.getNumericValue().floatValue();
		attrb=v.findAttribute("add_offset");		if(attrb!=null) offset=attrb.getNumericValue().floatValue();
		attrb=v.findAttribute("units");			if(attrb!=null) unit  =attrb.getStringValue();
		
		if(v.findAttribute("missing_value")!=null)
			undef=v.findAttribute("missing_value").getNumericValue().floatValue()*factor+offset;
		
		vname=v.getFullName();
		
		attrb=v.findAttribute("var_desc");	if(attrb!=null) comment=attrb.getStringValue().replace("\n"," ");
	}
}
