
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.0.0 - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 2.0.0 
 
 GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "times.h"
#include "geotop_common.h"
#include "inputKeywords.h"

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

//void set_time_step(PAR *par, TIMES *times){
  void set_time_step(Par *par, Times *times){
	
	double JDold;
	static double JDnew;
	static long line_interp, posix;

	if(times->time==0) posix=0;
	
	if(par->tsteps_from_file==1){
		
		if(times->time==0){
			line_interp=0;
		//	JDnew=par->init_date->co[geotop::common::Variables::i_sim]-1.0;
			JDnew=par->init_date[geotop::common::Variables::i_sim]-1.0;
		}
		
		JDold=JDnew;	
	//	JDnew=convert_tfromstart_JDfrom0(times->time, par->init_date->co[geotop::common::Variables::i_sim]);
		JDnew=convert_tfromstart_JDfrom0(times->time, par->init_date[geotop::common::Variables::i_sim]);
		
		if(floor(JDold)!=floor(JDnew)){
		//	time_no_interp(0, &line_interp, times->Dt_vector, times->Dt_matrix, times->numlinesDt_matrix, GTConst::max_cols_time_steps_file+1, 0,
		//				 par->init_date->co[geotop::common::Variables::i_sim]+times->time/86400.);
			time_no_interp(0, &line_interp, times->Dt_vector, times->Dt_matrix, times->numlinesDt_matrix, GTConst::max_cols_time_steps_file+1, 0,
								 par->init_date[geotop::common::Variables::i_sim]+times->time/86400.);
			posix=0;
		}
	}
	
	if( (long)times->Dt_vector[1+posix]== geotop::input::gDoubleNoValue ) posix=0;
	par->Dt=times->Dt_vector[1+posix];
	posix++;
		
//	if(par->plot_discharge_with_Dt_integration->co[geotop::common::Variables::i_sim]==1) par->Dtplot_discharge->co[geotop::common::Variables::i_sim]=par->Dt;
	if(par->plot_discharge_with_Dt_integration[geotop::common::Variables::i_sim]==1) par->Dtplot_discharge[geotop::common::Variables::i_sim]=par->Dt;
//	if(par->plot_point_with_Dt_integration->co[geotop::common::Variables::i_sim]==1) par->Dtplot_point->co[geotop::common::Variables::i_sim]=par->Dt;
	if(par->plot_point_with_Dt_integration[geotop::common::Variables::i_sim]==1) par->Dtplot_point[geotop::common::Variables::i_sim]=par->Dt;
//	if(par->plot_basin_with_Dt_integration->co[geotop::common::Variables::i_sim]==1) par->Dtplot_basin->co[geotop::common::Variables::i_sim]=par->Dt;
	if(par->plot_basin_with_Dt_integration[geotop::common::Variables::i_sim]==1) par->Dtplot_basin[geotop::common::Variables::i_sim]=par->Dt;

}
	
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

short is_leap(long y){

	short i=0;

	if(fmod(y,4.0)==0.0){
		i=1;
		if(fmod(y,100.0)==0.0){
			if(fmod(y,400)!=0.0) i=0;
		}
	}

	return(i);
}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*Calculates the fractional Julian day from 01/01/0000 0:00 (set at Julian Day = 1) given the year and Julian Day 
 of the year (from 1-Jan 0:00 of the year). This is the reference used in Matlab*/

double convert_JDandYear_JDfrom0(double JD, long year){
	
	long days;
	FILE *f;
	
	if (year < 0){
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: Not Possible to handle dates with year earlier than 0\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	

		days = 0;
	}else if (year == 0) {
		days = 1;//offset(01/01/0000 corrsepond to JDfrom0=1)
	}else {
		days = 1;//offset
		days += 366;//days of year 0
		days += (year-1)*365 + floor((year-1)/4.) - floor((year-1)/100.) + floor((year-1)/400.);
	}
	
	return(JD+(double)days);
}	
	
void convert_JDfrom0_JDandYear(double JDfrom0, double *JD, long *year){
	
	FILE *f;
	
	if (JDfrom0 < 1) {
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: Not Possible to handle dates with year earlier than 0\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
		
	}else if (JDfrom0 < 367) {
		*JD = JDfrom0 - 1.;
		*year = 0;
	}else {
		*year = floor( (JDfrom0-1.) / 365.25 );
		*JD = JDfrom0 - convert_JDandYear_JDfrom0(0.0, *year);
		while (*JD<0 || *JD>=365+is_leap(*year)) {
			if (*JD<0) {
				*year = *year - 1;
				*JD = *JD + 365 + is_leap(*year);
			}else {
				*JD = *JD - 365 - is_leap(*year);
				*year = *year + 1;
			}
		}
	}
	
}

double convert_JDfrom0_JD(double JDfrom0){
	
	double JD;
	long year;
	FILE *f;
	
	if (JDfrom0 < 1) {
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: Not Possible to handle dates with year earlier than 0\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
		JD = 0.0;
	}else if (JDfrom0 < 367) {
		JD = JDfrom0 - 1.;
	}else {
		year = floor( (JDfrom0-1.) / 365.25 );
		JD = JDfrom0 - convert_JDandYear_JDfrom0(0.0, year);
		while (JD<0 || JD>=365+is_leap(year)) {
			if (JD<0) {
				year -= 1;
				JD += (365 + is_leap(year));
			}else {
				JD -= (365 + is_leap(year));
				year += 1;
			}
		}
	}
	
	return JD;
	
}

/*==================================================================================================================*/

void convert_JDandYear_daymonthhourmin(double JD, long year, long *d, long *m, long *h, long *min){

	long JDint=floor(JD)+1;//integer Julian day
	double frac=JD-floor(JD);
	long days_in_months_in_leap[ ]={0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	long days_in_months_in_non_leap[ ]={0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
			
	*m = 1;//initialization
	
	if(is_leap(year) == 1){
		while (JDint > days_in_months_in_leap[*m]){
			JDint -= days_in_months_in_leap[*m];
			*m = (*m) + 1;
		}
	}else {
		while (JDint > days_in_months_in_non_leap[*m]){
			JDint -= days_in_months_in_non_leap[*m];
			*m = (*m) + 1;
		}		
	}
	
	*d = JDint;
	
	*min= ((int)floor(frac*((double)24.0*60.0) + 0.5)) % 60;
	*h = (int) floor(((((double)1440.0)*frac-(double)*min)/(double)60.0) + 0.5);
		
}

double convert_daymonthyearhourmin_JD(long d, long m, long y, long h, long min){
	
	double JD = 0.0;
	
	long i;
	long days_in_months_in_leap[ ]={0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	long days_in_months_in_non_leap[ ]={0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	
	if (is_leap(y) == 1) {
		for (i=1; i<=m-1; i++) {
			JD += (double)days_in_months_in_leap[i];
		}	
	}else {
		for (i=1; i<=m-1; i++) {
			JD += (double)days_in_months_in_non_leap[i];
		}			
	}

	JD += (double)d - 1.;
	JD += (double)h/24. + (double)min/60.;
	
	return JD;
	
}

/*==================================================================================================================*/

void convert_dateeur12_daymonthyearhourmin(double date, long *day, long *month, long *year, long *hour, long *min){
	
	FILE *f;
		
	*day = floor(date/1.E10);
	*month = floor(date/1.E8 - (*day)*1.E2);
	*year = floor(date/1.E4 - (*day)*1.E6 - (*month)*1.E4);
	*hour = floor(date/1.E2 - (*day)*1.E8 - (*month)*1.E6 - (*year)*1.E2 );
	*min = floor(date/1.E0 - (*day)*1.E10 - (*month)*1.E8 - (*year)*1.E4 - (*hour)*1.E2);
	
	if (*day<1 || *day>31 || *month<1 || *month>12 || *year<1700 || *year>2900 || *hour<0 || *hour>23 || *min<0 || *min>59) {
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: the date %12.0f cannot be read\n",date);
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}
}
	
double convert_daymonthyearhourmin_dateeur12(long day, long month, long year, long hour, long min){

	return(day*1.E10 + month*1.E8 + year*1.E4 + hour*1.E2 + min*1.E0);
	
}

/*==================================================================================================================*/

void convert_dateeur12_JDandYear(double date, double *JD, long *year){
	
	long i, day, month, hour, min;
	long days_in_months_in_leap[ ]={0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	long days_in_months_in_non_leap[ ]={0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	
	convert_dateeur12_daymonthyearhourmin(date, &day, &month, year, &hour, &min);

	*JD = 0.0;
	if (is_leap(*year) == 1) {
		for (i=1; i<month; i++) {
			*JD = (*JD) + (double)days_in_months_in_leap[i];
		}
	}else {
		for (i=1; i<month; i++) {
			*JD = (*JD) + (double)days_in_months_in_non_leap[i];
		}
	}
	
	*JD = (*JD) + (double)(day - 1);
	*JD = (*JD) + (double)hour/24.;
	*JD = (*JD) + (double)min/(24.*60.);
		
}
	
double convert_JDandYear_dateeur12(double JD, long year){
	
	long d, m, h, min;
	
	convert_JDandYear_daymonthhourmin(JD, year, &d, &m, &h, &min);
	return convert_daymonthyearhourmin_dateeur12(d, m, year, h, min);

}
		   
/*==================================================================================================================*/

double convert_dateeur12_JDfrom0(double date){
	
	double JD;
	long year;
	
	convert_dateeur12_JDandYear(date, &JD, &year);
	return convert_JDandYear_JDfrom0(JD, year);
}

double convert_JDfrom0_dateeur12(double JDfrom0){

	double JD;
	long year;

	convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
	return convert_JDandYear_dateeur12(JD, year);
}

/*==================================================================================================================*/

double convert_tfromstart_JDfrom0(double t, double JDfrom0_start){

	return JDfrom0_start + t/86400.;

}

double convert_JDfrom0_tfromstart(double JDfrom0, double JDfrom0_start){

	return (JDfrom0 - JDfrom0_start)*86400.;

}

/*==================================================================================================================*/
