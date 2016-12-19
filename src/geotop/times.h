
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1 
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
    
    
/****************************************************************************************************/
/* The dates, times and counters of the simulation are updated:                                      */
/****************************************************************************************************/
#ifndef TIMES_H
#define TIMES_H
#include "struct.geotop.h"
#include "constants.h"
#include "meteo.h"
#include "meteodata.h"
#include <exception>

class InvalidDateException : std::exception
{
    public:
        const char* what()
        {
            return "Invalid date.";
        }
};

void set_time_step(Par *par, Times *times);

short is_leap(long year);

double convert_JDandYear_JDfrom0(double JD, long year);

void convert_JDfrom0_JDandYear(double JDfrom0, double *JD, long *year);

double convert_JDfrom0_JD(double JDfrom0);

void convert_JDandYear_daymonthhourmin(double JD, long year, long *d, long *m, long *h, long *min);

double convert_daymonthyearhourmin_JD(long d, long m, long y, long h, long min);

void convert_dateeur12_daymonthyearhourmin(double date, long *day, long *month, long *year, long *hour, long *min);

double convert_daymonthyearhourmin_dateeur12(long day, long month, long year, long hour, long min);

void convert_dateeur12_JDandYear(double date, double *JD, long *year);

double convert_JDandYear_dateeur12(double JD, long year);

double convert_dateeur12_JDfrom0(double date);

double convert_JDfrom0_dateeur12(double JDfrom0);

double convert_tfromstart_JDfrom0(double t, double JDfrom0_start);

double convert_JDfrom0_tfromstart(double JDfrom0, double JDfrom0_start);
#endif
