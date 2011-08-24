
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

void meteo_interp(short flag, short update, long *istart, double *out, double **data, long nlines, long ncols, long col_date, double tbeg, double tend);
void meteo_interp2(short flag, long *istart, double *out, double **data, long nlines, long ncols, long col_date, double tbeg, double tend);
short fixing_dates(long imeteo, double **data, double ST, double STstat, long nlines, long date12col, long JDfrom0col);
short fill_wind_speed(double **data, long nlines, long Wspeed, long Wdir, long Wx, long Wy, char *HeaderWx, char *HeaderWy);
short fill_Tdew(long imeteo, DOUBLEVECTOR *Z, double **data, long nlines, long RH, long Tair, long Tairdew, char *HeaderTdew, double RHmin);
void rewrite_meteo_files(double **meteo, long meteolines, char **header, char *name, short added_JD, short added_wind, short added_cloudiness, short added_Tdew);
double integrate_meas(short flag, double t, long i, double **data, long col, long col_date);
long find_line_data(short flag, double t, long ibeg, double **data, long col_date, long nlines, short *a);
double time_in_JDfrom0(short flag, long i, long col, double **data);
long find_station(long metvar, long nstat, double **var);
double **read_horizon(short a, long i, char *name, char **ColDescr, long *num_lines, FILE *flog);
