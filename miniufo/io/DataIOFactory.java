/**
 * @(#)DataIOFactory.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.io;

import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.descriptor.CtsDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.descriptor.NetCDFDescriptor;


/**
 * factory for I/O
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class DataIOFactory{
	
	/**
	 * prevent from instantiate
     */
	private DataIOFactory(){}
	
	
	/**
	 * static factory method
     */ 
	public static DataRead getDataRead(DataDescriptor dd){
		if(dd instanceof    CtlDescriptor) return new CtlDataReadStream((   CtlDescriptor)dd);
		if(dd instanceof    CtsDescriptor) return new CtlDataReadStream((   CtsDescriptor)dd);
		if(dd instanceof    CsmDescriptor) return new CsmDataReadStream((   CsmDescriptor)dd);
		if(dd instanceof NetCDFDescriptor) return new  NetCDFReadStream((NetCDFDescriptor)dd);
		
		throw new IllegalArgumentException("unsupported DataDescriptor type");
	}
    
	public static DataWrite getDataWrite(DataDescriptor dd,String path){
		if(dd instanceof CsmDescriptor) return new CsmDataWriteStream(path);
		else if(dd instanceof CtsDescriptor) return new CtsDataWriteStream(path);
		else return new CtlDataWriteStream(path);
	}
}
