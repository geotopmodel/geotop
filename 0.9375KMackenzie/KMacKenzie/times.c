
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion MacLavagna

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 MacLavagna.
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



#include "struct.geotop.09375.h"
#include "times.h"

#include "constant.h"
#include "keywords_file.h"
#include "pedo.funct.h"
#include "write_dem.h"
#include "t_datamanipulation.h"
#include "geo_statistic.09375.h"
#include "t_random.h"
#include "networks.h"
#include "t_utilities.h"

extern STRINGBIN *files;
extern T_INIT *UV;



/**********************************************************************************************************/
/**********************************************************************************************************/
/*update of the vector with times controls and counters:*/
/**********************************************************************************************************/
/**********************************************************************************************************/
void updates_times(TIMES *times, PAR *par){

	long i, N;
	FILE *f;
	double tmin, tmax;
	short occurring;

	/* for global output each dt pixel cancel the following if */
	if(times->i_pixel==times->n_pixel) times->i_pixel=0;
	times->i_pixel+=1;
	if(times->i_basin==times->n_basin) times->i_basin=0;
	times->i_basin+=1;
	date_time(times->time, par->year0, par->JD0, 0.0, &(times->JD), &(times->DD), &(times->MM), &(times->AAAA), &(times->hh), &(times->mm));
	if(par->JD_plots->co[1]!=0){
		occurring=0;
		for(i=1;i<=par->JD_plots->nh;i++){
			tmin=get_time( (double)(par->JD_plots->co[i]-1), times->AAAA, par->JD0, par->year0 );
			tmax=get_time( (double)(par->JD_plots->co[i]  ), times->AAAA, par->JD0, par->year0 );
			if(fmod(tmax-tmin,times->n_plot*par->Dt)!=0.0){
				N=floor(tmax-tmin/(times->n_plot*par->Dt))+1;
				tmax=tmin+N*times->n_plot*par->Dt;
			}
			//printf("i:%ld JD:%ld t:%f tmin:%f tmax:%f\n",i,par->JD_plots->co[i],times->time,tmin,tmax);
			if(times->time>=tmin && times->time-par->Dt<tmin){
				times->i_plot=0;
				times->nt_plot=1;
				times->d_plot=par->JD_plots->co[i];
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

	date_time(times->time+par->Dt, par->year0, par->JD0, 0.0, &(times->JD), &(times->DD), &(times->MM), &(times->AAAA), &(times->hh), &(times->mm));

	if (times->time>=times->TH*3600.0){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"\nTime table [s]: \nt_E=%10.2f   t_Wv=%10.2f   t_Wh=%10.2f   t_write=%10.2f\n", times->egy, times->vert_wb, times->horiz_wb, times->writeout);
		fclose(f);
	}

}



/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void date_time(double t, long y0, double JDstart, double delay, double *JD, long *d, long *m, long *y, long *h, long *min){

	short i;
	long JDint;

	i=is_leap(y0);
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

	*h=floor((*JD+1-JDint)*24);
	*min=floor((*JD+1-JDint)*1440-*h*60);

}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/


void time_conversion(double JD01, long Y01, double t1, double JD02, long Y02, double *t2){

	double dJD;
	long y;

	y=Y02;
	dJD=JD01-JD02;

	do{
		if(y<Y01){
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

double get_time(double JD1, long y1, double JD0, long y0){

	double t, dJD;
	long y;

	dJD=JD1-JD0;

	for(y=y0;y<y1;y++){
		dJD+=(365.0+is_leap(y));
	}

	if(dJD<0) dJD+=(365.0+is_leap(y1));		//in case it is negative, it is considered the following year

	t=dJD*86400.0;

	if(t<0) t_error("Negative time");

	return(t);
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
