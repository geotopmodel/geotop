
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 2.0.0
 
 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
    
    
/****************************************************************************************************/
/* The dates, times and counters of the simulation are updated:                                      */
/****************************************************************************************************/

void set_time_step(PAR *par, TIMES *times);

short is_leap(long y);


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
