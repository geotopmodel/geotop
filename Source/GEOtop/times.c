
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */


#include "struct.geotop.h"
#include "times.h"
#include "constant.h"
#include "keywords_file.h"

extern STRINGBIN *files;



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
	if(times->i_pixel==times->n_pixel) times->i_pixel=0;
	times->i_pixel+=1;
	if(times->i_discharge==times->n_discharge) times->i_discharge=0;
	times->i_discharge+=1;
	if(times->i_basin==times->n_basin) times->i_basin=0;
	times->i_basin+=1;
	date_time(times->time, par->year0, par->JD0, 0.0, &(times->JD), &(times->day), &(times->month), &(times->year), &(times->hour), &(times->min));
	if(par->JD_plots->co[1]>=0){
		occurring=0;
		for(i=1;i<=(long)(par->JD_plots->nh/2.);i++){
			get_time( &tmin, par->JD_plots->co[2*i-1], times->year, par->JD0, par->year0 );
			get_time( &tmax, par->JD_plots->co[2*i]  , times->year, par->JD0, par->year0 );
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

	date_time(times->time+0.5*par->Dt, par->year0, par->JD0, 0.0, &(times->JD), &(times->day), &(times->month), &(times->year), &(times->hour), &(times->min));

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
	double frac;

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
	if(*min<0) *min=0;
	
	frac = *JD - floor(*JD);
	*min= ((int)floor(frac*((double)24.0*60.0) + 0.5)) % 60;
	*h = (int) floor(((((double)1440.0)*frac-(double)*min)/(double)60.0) + 0.5);
	

}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/


void time_conversion(double JD01, long Y01, double t1, double JD02, long Y02, double *t2){

	double dJD;
	long y;

	//printf("JD01:%f Y01:%ld t1:%f JD02:%f Y02:%ld\n",JD01,Y01,t1,JD02,Y02);

	y=Y02;
	dJD=JD01-JD02;

	//printf("y:%ld dJD:%f\n",y,dJD);

	do{
		if(y<Y01){
			dJD+=(365.0+is_leap(y));
			y++;
		}else if(y>Y01){
			dJD-=(365.0+is_leap(y-1));
			y--;
		}
	}while(y!=Y01);

	//printf("y:%ld Y01:%ld\n",y,Y01);

	*t2=t1+dJD*86400;		//s

	//printf("dJD:%f\n",dJD);

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
