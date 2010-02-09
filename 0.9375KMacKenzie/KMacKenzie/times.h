
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
    
    
#ifndef __TIMES_H__
#define __TIMES_H__    

#include "keywords_file.h"
#include "struct.geotop.09375.h"
#include "constant.h"
    
/****************************************************************************************************/
/* The dates, times and counters of the simulation are updated:                                      */
/****************************************************************************************************/
void updates_times(TIMES *times, PAR *par);

void date_time(double t, long y0, double JDstart, double delay, double *JD, long *d, long *m, long *y, long *h, long *min);

void time_conversion(double JD01, long Y01, double t1, double JD02, long Y02, double *t2);

void get_time(double *t, double JD1, long y1, double JD0, long y0);

short is_leap(long y);

long daysfrom0(long year);

double convert_date_to_julian(int year, int month, int day, int hour, int minute); 

#endif
