/**
 * @(#)LittleRReport.java	1.0 07/02/01
 *
 * Copyright 2007 MiniUFO, All rights reserved.
 * MiniUFO Studio. Use is subject to license terms.
 */
package miniufo.test.IO;


/**
 * little-R record format
 *
 * @version 1.0, 02/01/2007
 * @author  MiniUFO
 * @since   MDK1.0
 */
public final class LittleRReport{
	//
	private HeaderRecord hr=null;
	private DataRecord   dr=null;
	private EndRecord    er=null;
	
	
	/**
     * constructor
     */
	public LittleRReport(HeaderRecord hr,DataRecord dr,EndRecord er){
		this.hr=hr;System.out.println(this.hr);
		this.dr=dr;System.out.println(this.dr);
		this.er=er;System.out.println(this.er);
	}
	
	
	static final class HeaderRecord{
		int sut;
		int julian;
		int seq_num;
		int num_dups;
		int num_vld_fld;
		int num_error;
		int num_warning;
		
		boolean bogus;
		boolean discard;
		boolean is_sound;
		
		float latitude;
		float longitude;
		float elevation;
		
		String id;
		String name;
		String source;
		String platform;
		String date_char;
		
		QCRecord slp;
		QCRecord sst;
		QCRecord psfc;
		QCRecord t_max;
		QCRecord t_min;
		QCRecord precip;
		QCRecord ceiling;
		QCRecord ref_pres;
		QCRecord ground_t;
		QCRecord p_tend03;
		QCRecord p_tend24;
		QCRecord cloud_cvr;
		QCRecord t_min_night;
	}
	
	static final class DataRecord{
		QCRecord pressure;
		QCRecord height;
		QCRecord temperature;
		QCRecord dew_point;
		QCRecord speed;
		QCRecord direction;
		QCRecord u;
		QCRecord v;
		QCRecord rh;
		QCRecord thickness;
	}
	
	static final class EndRecord{
		int num_vld_fld;
		int num_error;
		int num_warning;
	}
	
	static final class QCRecord{
		//
		float data;
		int qc;
		
		
		/**
		 * constructor
		 */
		QCRecord(float data,int qc){
			this.data=data;
			this.qc  =qc;
		}
	}
	
	static enum QC{
		PressureInterpolatedFromFirstGuessHeight(1),
		TemperatureAndDewPointBoth(4),
		WindSpeedAndDirectionBoth(5),
		WindSpeedNegative(6),
		WindDirectionLT0OrGT360(7),
		LevelVerticallyInterpolated(8),
		ValueVerticallyExtrapolatedFromSingleLevel(9),
		SignOfTemperatureReversed(10),
		SuperadiabaticLevelDetected(11),
		VerticalSpikeInWindSpeedOrDirection(12),
		ConvectiveAdjustmentAppliedToTemperatureField(13),
		NoNeighboringObservationsForBuddyCheck(14),
		FailsErrorMaximumTest(15),
		FailsBuddyTest(16),
		ObservationOutsideOfDomainDetectedByQC(17);
		
		final int number;
		
		QC(int num){ number=1<<num;}
		
		public String toString(){ return super.toString()+" ("+number+")";}
	}
	
	
	// namelist class
	static class Record1{	// analysis times
		int start_year;		// 4-digit year of the starting time to process
		int start_month;	// 2-digit month of the starting time to process
		int start_day;		// 2-digit day of the starting time to process
		int start_hour;		// 2-digit hour of the starting time to process
		int end_year;		// 4-digit year of the ending time to process
		int end_month;		// 2-digit month of the ending time to process
		int end_day;		// 2-digit day of the ending time to process
		int end_hour;		// 2-digit hour of the ending time to process
		int interval;		// time interval (s) between consecutive times to process
	}
	
	static class Record2{	// filenames (may include directory information)
		// filename of the first-guess fields, there is only a single name
		String fg_filename;
		// filename(s) of the observation files, one required for each time period to run through the objective analysis
		String obs_filename;
		// filename(s) of the observation files to be used for the surface analyses option (only when .F4D=.TRUE.)
		String sfc_obs_filename;
	}
	
	static class Record3{ // concern space allocated within the program for observationstly
		// anticipated maximum number of reports per time period
		int max_number_of_obs;
		// T/F flag allows the user to decide the severity of
		// not having enough space to store all of the available observations
		boolean fatal_if_excedd_max_obs;
	}
	
	static class Record4{	// QC options
		// check the difference between the first-guess and the observation
		boolean qc_test_error_max;
		// check the difference between a single observation and neighboring observations
		boolean qc_test_buddy;
		// check for vertical spikes in temperature, dew point, wind speed and wind direction
		boolean qc_test_vert_consistency;
		// remove any super-adiabatic lapse rate in a sounding by conservation of dry static energy
		boolean qc_test_convective_adj;
	}
	
	
	/** test*/
	public static void main(String[] args){
		for(QC c:QC.values()){
			System.out.println(c.toString());
		}
	}
}
