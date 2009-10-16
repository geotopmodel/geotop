
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/






#include "keywords_file.h"
#include "struct.geotop.09375.h"
#include "times.h"
#include "constant.h"

extern STRINGBIN *files;
extern char *error_file_name;



/**********************************************************************************************************/
/**********************************************************************************************************/
/*update of the vector with times controls and counters:*/
/**********************************************************************************************************/
/**********************************************************************************************************/
void updates_times(TIMES *times, PAR *par){

	long i;
	FILE *f;
	double tmin, tmax, N;
	short occurring;

	/* for global output each dt pixel cancel the following if */
	if(times->i_pixel==times->n_pixel) {
		times->i_pixel=0;
		times->count_super_printed++;
	}
	times->i_pixel+=1;
	if(times->i_basin==times->n_basin) times->i_basin=0;
	times->i_basin+=1;
	date_time(times->time, par->year0, par->JD0, 0.0, &(times->JD), &(times->DD), &(times->MM), &(times->AAAA), &(times->hh), &(times->mm));
	if(par->JD_plots->co[1]>=0){
		occurring=0;
		for(i=1;i<=(long)(par->JD_plots->nh/2.);i++){
			get_time( &tmin, par->JD_plots->co[2*i-1], times->AAAA, par->JD0, par->year0 );
			get_time( &tmax, par->JD_plots->co[2*i]  , times->AAAA, par->JD0, par->year0 );
			if( floor((tmax-tmin)/(times->n_plot*par->Dt)) != (tmax-tmin)/(times->n_plot*par->Dt) ){
				N=floor(((tmax-tmin)/(times->n_plot*par->Dt)));
				tmax=tmin+N*times->n_plot*par->Dt;
			}

			if(times->time>=tmin && times->time-par->Dt<tmin){
				times->i_plot=0;
				times->nt_plot=1;
				times->d_plot=floor(par->JD_plots->co[2*i-1])+1;
			}

			if(times->time>=tmin && times->time<tmax){
				if(times->i_plot==times->n_plot){
					times->i_plot=0;
					times->nt_plot+=1;
				}
				times->i_plot+=1;
				occurring=1;
			}
		}
		if(occurring==0) times->i_plot=0;

	}else{
		times->n_plot=1;
	}

	date_time(times->time+0.5*par->Dt, par->year0, par->JD0, 0.0, &(times->JD), &(times->DD), &(times->MM), &(times->AAAA), &(times->hh), &(times->mm));

	if (times->time>=times->TH*3600.0){
		f=fopen(error_file_name,"a");
		fprintf(f,"\nTime table [s]: \nt_E=%10.2f   t_Wv=%10.2f   t_Wh=%10.2f   t_write=%10.2f\n", times->egy, times->vert_wb, times->horiz_wb, times->writeout);
		fclose(f);
	}

}



/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void date_time(double t, long y0, double JDstart, double delay, double *JD, long *d, long *m, long *y, long *h, long *min){
	/* Author: Stefano Endrizzi   Year:
	* function that calculates the current date (year, month, day, hour, minute)
	* Input:
	* 			t: current time (sec) from the start of the simulation
	* 			y0: year of the beginning of the simulation
	* 		JDstart: Julian day of the beginning of the simulation
	* 		 delay:
	* Output:
	* 			JD: julian day of the time of plot
	* 			d: day of the time of plot
	* 			m: month of the time of plot
	* 			y: year of the time of plot
	* 			h: hour of the time of plot
	* 		  min: minute of the time of plot
	* comment: Matteo Dall'Amico, May 2009 */
	short i;
	long JDint;

	i=is_leap(y0); // 1 is leap, 0 is not leap
	*JD=JDstart+delay+t/86400.0;
	*y=y0;

	while(*JD>=365+i){
		*y+=1;
		*JD-=(365+i);
		i=is_leap(*y);
	}

	JDint=floor(*JD)+1;

	if(JDint<=31){
		*d=JDint;
		*m=1;
	}else if(JDint<=31+28+i){
		*d=JDint-31;
		*m=2;
	}else if(JDint<=31+28+i+31){
		*d=JDint-31-28-i;
		*m=3;
	}else if(JDint<=31+28+i+31+30){
		*d=JDint-31-28-i-31;
		*m=4;
	}else if(JDint<=31+28+i+31+30+31){
		*d=JDint-31-28-i-31-30;
		*m=5;
	}else if(JDint<=31+28+i+31+30+31+30){
		*d=JDint-31-28-i-31-30-31;
		*m=6;
	}else if(JDint<=31+28+i+31+30+31+30+31){
		*d=JDint-31-28-i-31-30-31-30;
		*m=7;
	}else if(JDint<=31+28+i+31+30+31+30+31+31){
		*d=JDint-31-28-i-31-30-31-30-31;
		*m=8;
	}else if(JDint<=31+28+i+31+30+31+30+31+31+30){
		*d=JDint-31-28-i-31-30-31-30-31-31;
		*m=9;
	}else if(JDint<=31+28+i+31+30+31+30+31+31+30+31){
		*d=JDint-31-28-i-31-30-31-30-31-31-30;
		*m=10;
	}else if(JDint<=31+28+i+31+30+31+30+31+31+30+31+30){
		*d=JDint-31-28-i-31-30-31-30-31-31-30-31;
		*m=11;
	}else{
		*d=JDint-31-28-i-31-30-31-30-31-31-30-31-30;
		*m=12;
	}
	// before was as in the commented lines
	//*h=floor((*JD+1-JDint)*24);// was before
	//*min=floor((*JD+1-JDint)*1440-*h*60);// was before
	// correction by Thomas Egger, see email of 6/10/2009
	double frac = *JD - floor(*JD);
	*min= ((int)floor(frac*((double)24.0*60.0) + 0.5)) % 60;
	*h = (int) floor(((((double)1440.0)*frac-(double)*min)/(double)60.0) + 0.5);
}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/


void time_conversion(double JD01, long Y01, double t1, double JD02, long Y02, double *t2){
	/* Author:    Year:
	 * function that calculates the current time [second] based on a new origin given by the beginning of the
	 * first datum of the meteo file.
	 * Input:
	 * 			JD01: julian day of the beginning of the simulation
	 * 			Y01: year of the beginning of the simulation
	 * 			t1: current time counter of the simulation
	 * 			JD02: julian day of the first datum of the meteo station
	 * 			Y02: year of the first datum of the meteo station
	 * Output:
	 * 			t2:  current time [second] based on a new origin given by the beginning of the first datum of the meteo file
	 * t2<0: the simulation started prior of the meteo data and the current time counter hasn't reached the beginning of the meteo data yet.
	 * t2>0: the simulation started after of the meteo data or the current time counter has reached the beginning of the meteo data
	 * comment: Matteo Dall'Amico, April 2009 */
	double dJD;
	long y;

	y=Y02;
	dJD=JD01-JD02;// the difference between the decimal julian day of the current day and that of the meteo station

	do{
		if(y<Y01){// the year of the first datum of the meteo station < the current year of the simulation
			dJD+=(365.0+is_leap(y));
			y++;
		}else if(y>Y01){
			dJD-=(365.0+is_leap(y-1));
			y--;
		}
	}while(y!=Y01);

	*t2=t1+dJD*86400;		//s

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void get_time(double *t, double JD1, long y1, double JD0, long y0){

	double dJD;
	long y;

	dJD=JD1-JD0;

	for(y=y0;y<y1;y++){
		dJD+=(365.0+is_leap(y));
	}

	if(dJD<0) dJD+=(365.0+is_leap(y1));		//in case it is negative, it is considered the following year

	*t=dJD*86400.0;

	if(*t<0) t_error("Negative time");

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

long daysfrom0(long year){

	long days;

	days=year*365 + floor(year/4.0) - floor(year/100.0) + floor(year/400.0) + 2;

	return(days);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
