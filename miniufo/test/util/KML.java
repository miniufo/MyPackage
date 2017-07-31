/**
 * @(#)KML.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.util;


/**
 * used to great the data to kml format in Google earth
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class KML{
	
	//
	public static StringBuilder newKml(StringBuilder content){
		StringBuilder sb=new StringBuilder();
		sb.append("<?xml version=\"1.0\"?>\n");
		sb.append("<kml xmlns=\"http://earth.google.com/kml/2.0\">\n");
		sb.append(content);
		sb.append("</kml>\n");
		return sb;
	}
	
	public static StringBuilder newDocument(StringBuilder content,String name,boolean isopen){
		StringBuilder sb=new StringBuilder();
		sb.append("<Document>\n");
		if(name!=null) sb.append(KML.newName(name));
		sb.append(KML.newOpen(isopen));
		sb.append(content);
		sb.append("</Document>\n");
		return sb;
	}
	
	public static StringBuilder newFolder(StringBuilder content,String name,boolean isopen){
		StringBuilder sb=new StringBuilder();
		sb.append("<Folder>\n");
		if(name!=null) sb.append(KML.newName(name));
		sb.append(KML.newOpen(isopen));
		sb.append(content);
		sb.append("</Folder>\n");
		return sb;
	}
	
	public static StringBuilder newPlacemark(StringBuilder content,String name){
		StringBuilder sb=new StringBuilder();
		sb.append("<Placemark>\n");
		if(name!=null) sb.append(KML.newName(name));
		sb.append(content);
		sb.append("</Placemark>\n");
		return sb;
	}
	
	public static StringBuilder newStyle(StringBuilder stylecontent,String id){
		StringBuilder sb=new StringBuilder();
		if(id!=null) sb.append("<style id=\""+id+"\">\n");
		else sb.append("<Style>\n");
		sb.append(stylecontent);
		sb.append("</Style>\n");
		return sb;
	}
	
	public static StringBuilder newLineString(String[] lon,String[] lat,String[] alt,StringBuilder altitudeMode){
		StringBuilder sb=new StringBuilder();
		int N=lon.length;
		
		if(lat.length!=N||alt.length!=N) throw new IllegalArgumentException("lengths not same");
		
		sb.append("<LineString>\n");
		if(altitudeMode!=null)	sb.append(altitudeMode);
		
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++){
			sb.append(lon[i]+",");
			sb.append(lat[i]+",");
			sb.append(alt[i]+",\n");
		}
		sb.append("</coordinates>\n");
		
		sb.append("</LineString>\n");
		
		return sb;
	}
		
	public static StringBuilder newLinearRing(String[] lon,String[] lat,String[] alt,StringBuilder altitudeMode){
		StringBuilder sb=new StringBuilder();
		int N=lon.length;
		
		if(lat.length!=N||alt.length!=N) throw new IllegalArgumentException("lengths not same");
		if(!lon[0].equals(lon[N-1])||!lat[0].equals(lat[N-1])) throw new IllegalArgumentException("illegal data");
		
		sb.append("<LinearRing>\n");
		if(altitudeMode!=null)	sb.append(altitudeMode);
		
		sb.append("<coordinates>\n");
		for(int i=0;i<N;i++){
			sb.append(lon[i]+",");
			sb.append(lat[i]+",");
			sb.append(alt[i]+",\n");
		}
		sb.append("</coordinates>\n");
		
		sb.append("</LinearRing>\n");
		
		return sb;
	}
	
    /** newPolygon
    *
    * @param	outerBoundary	necessarily, an element of newLinerRing
    * @param	innerBoundary	unnecessarily, elements of newLinerRing
    */
	public static StringBuilder newPolygon(StringBuilder outerBoundary,StringBuilder altitudeMode,StringBuilder... innerBoundary){
		StringBuilder sb=new StringBuilder();
		
		sb.append("<Polygon>\n");
		if(innerBoundary!=null)	sb.append(altitudeMode);
		
		sb.append("<outerBoundaryIs>\n");
		sb.append(outerBoundary);
		sb.append("</outerBoundaryIs>\n");
		
		if(innerBoundary!=null){
			for(int i=0;i<innerBoundary.length;i++) {
				sb.append("<innerBoundaryIs>\n");
				sb.append(innerBoundary[i]);
				sb.append("</innerBoundaryIs>\n");
			}
		}
		
		sb.append("</Polygon>\n");
		
		return sb;
	}
	
	
	public static StringBuilder newName(String name){
		StringBuilder sb=new StringBuilder();
		sb.append("<name>");
		sb.append(name);
		sb.append("</name>\n");
		return sb;
	}
	
	public static StringBuilder newOpen(boolean isopen){
		StringBuilder sb=new StringBuilder();
		if(isopen) sb.append("<open>1</open>\n");
		else sb.append("<open>0</open>\n");
		return sb;
	}
	
	public static StringBuilder newDescription(String description){
		StringBuilder sb=new StringBuilder();
		sb.append("<description>");
		sb.append(description);
		sb.append("</description>\n");
		return sb;
	}
		
	public static StringBuilder newStyleUrl(String id){
		StringBuilder sb=new StringBuilder();
		sb.append("<styleUrl>#"+id+"</styleUrl>\n");
		return sb;
	}
		
	public static StringBuilder newLineStyle(String width,String color){
		StringBuilder sb=new StringBuilder();
		sb.append("<LineStyle>\n");
		if(color!=null) sb.append("<color>"+color+"</color>\n");
		sb.append("<width>"+width+"</width>\n");
		sb.append("</LineStyle>\n");
		return sb;
	}
	
	public static StringBuilder newPolyStyle(boolean fill,boolean outline,String width,String color){
		StringBuilder sb=new StringBuilder();
		
		sb.append("<PolyStyle>\n");
		if(fill) sb.append("<fill>1</fill>\n");
		else sb.append("<fill>0</fill>\n");
		if(outline) sb.append("<outline>1</outline>\n");
		else sb.append("<outline>1</outline>\n");
		if(color!=null) sb.append("<color>"+color+"</color>\n");
		if(width!=null) sb.append("<width>"+width+"</width>\n");
		sb.append("</LineStyle>\n");
		return sb;
	}
	
    /** newAltitudeMode
    *
    * @param	mode	mode of altitude, can be set to clampToGround, relativeToGround, absolute  
    */
	public static StringBuilder newAltitudeMode(String mode){
		StringBuilder sb=new StringBuilder();
		sb.append("<extrude>1</extrude>\n");
		sb.append("<tessellate>1</tessellate>\n");
		sb.append("<altitudeMode>"+mode+"</altitudeMode>\n");
		return sb;
	}
	
    /** newRegion
    *
    * @param	LatLonAltBox	coordinates and altitudes,
    * 							the order must be north, south, east, west, minAltitude, maxAltitude
    * @param	Lod				Level of Detail, at lest contains one data of minLodPixels,
    * 							others can be maxLodPixels, minFadeExtent, maxFadeExtent
    */
	public static StringBuilder newRegion(String[] LatLonAltBox,String... Lod){
		StringBuilder sb=new StringBuilder();
		
		if(LatLonAltBox.length!=6) throw new IllegalArgumentException("Length of data must be 6");
		sb.append("<Region>\n");
		sb.append("<LatLonAltBox>\n");
		sb.append("<north>"+LatLonAltBox[0]+"</north>\n");
		sb.append("<south>"+LatLonAltBox[1]+"</south>\n");
		sb.append("<east>"+LatLonAltBox[2]+"</east>\n");
		sb.append("<west>"+LatLonAltBox[3]+"</west>\n");
		sb.append("<minAltitude>"+LatLonAltBox[4]+"</minAltitude>\n");
		sb.append("<maxAltitude>"+LatLonAltBox[5]+"</maxAltitude>\n");
		sb.append("</LatLonAltBox>\n");
		
		if(Lod.length>4) throw new IllegalArgumentException("Length of data must be less than 4");
		sb.append("<Lod>\n");
		sb.append("<minLodPixels>"+Lod[0]+"</minLodPixels>\n");
		if(Lod.length==2){
			sb.append("<maxLodPixels>"+Lod[1]+"</maxLodPixels>\n");
		}
		if(Lod.length==3){
			sb.append("<maxLodPixels>"+Lod[1]+"</maxLodPixels>\n");
			sb.append("<minFadeExtent>"+Lod[2]+"</minFadeExtent>\n");
		}
		if(Lod.length==4){
			sb.append("<maxLodPixels>"+Lod[1]+"</maxLodPixels>\n");
			sb.append("<minFadeExtent>"+Lod[2]+"</minFadeExtent>\n");
			sb.append("<maxFadeExtent>"+Lod[3]+"</maxFadeExtent>\n");
		}
		sb.append("</Lod>\n");
		sb.append("</Region>\n");
		
		return sb;		
	}
	
	/** test
	public static void main(String[] arg){
		long l=System.currentTimeMillis();
		WGrib wg=new WGrib("//LYNN/Data/Typhoon/200605/grib2006050100","d:/Data/Typhoon/CHANCHU/test.dat");
		System.out.println((System.currentTimeMillis()-l)/1000f);
		
		l=System.currentTimeMillis();
		wg.outputAll();
		System.out.println((System.currentTimeMillis()-l)/1000f);
		
		wg.closeFile();
	}*/
}
