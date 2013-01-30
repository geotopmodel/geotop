
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225-9 'Moab' - 24 Aug 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225-9 'Moab'
 
 GEOtop 1.225-9 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225-9 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
#ifndef OUTPUT_H
#define OUTPUT_H

#include "struct.geotop.h"
#include "pedo.funct.h"
#include "../libraries/geomorphology/networks.h"
#include "../libraries/ascii/rw_maps.h"
#include "constants.h"
#include "../libraries/ascii/extensions.h"
#include "times.h"
#include "energy.balance.h"
#include "input.h"
#include "../libraries/ascii/tabs.h"
#include "vegetation.h"
#include "tables.h"
#include "snow.h"
#include "../libraries/ascii/init.h"
#include "water.balance.h"


#include <errno.h>
#include <time.h>

extern long number_novalue, number_absent;
extern char *string_novalue;

extern char *WORKING_DIRECTORY;

extern T_INIT *UV;
extern char **files, *logfile;
extern long Nl, Nr, Nc;

extern double t_meteo, t_energy, t_water, t_sub, t_sup, t_out, t_blowingsnow;

extern double **odpnt, **odp, *odbsn, *odb;
extern long *opnt, nopnt, *obsn, nobsn, *osnw, nosnw;
extern long *oglc, noglc, *osl, nosl;
extern char **hpnt, **hbsn, **hsnw, **hglc, **hsl;

extern char *keywords_num[num_par_number] , *keywords_char[num_par_char];

extern FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav, *ffice, *fficeav, *ffsnow, *ffglac;

extern long i_sim, i_run;

extern time_t start_time;
extern double elapsed_time, elapsed_time_start, cum_time, max_time;




void write_output(TIMES *times,WATER *wat,CHANNEL *cnet,PAR *par,TOPO *top,LAND *land,SOIL *sl,ENERGY *egy,SNOW *snow,GLACIER *glac,METEO *met);

void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac);

void write_soil_output(long i, long iname, double init_date, double JDfrom0, double JD, long day, long month, long year, long hour, long minute, DOUBLEVECTOR *n, SOIL *sl, PAR *par, double psimin, double cosslope);

void write_soil_file(long lmin, long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, double *var, DOUBLEVECTOR *n, double *dz, double cosslope);

void write_soil_header(FILE *f, DOUBLEVECTOR *n, double *dz);

void plot(char *name, long i_plot, DOUBLEVECTOR *V, short format, long **J);

double interpolate_soil(long lmin, double h, long max, double *Dz, double *Q);

double interpolate_soil2(long lmin, double h, long max, double *Dz, DOUBLEMATRIX *Q, long i);

void write_tensorseries_soil(long lmin, char *suf, char *filename, short type, short format, DOUBLEMATRIX *T, DOUBLEVECTOR *n, long **J, LONGMATRIX *RC, double *dz, DOUBLEMATRIX *slope, short vertical);

void fill_output_vectors(double Dt, double W, ENERGY *egy, SNOW *snow, GLACIER *glac, WATER *wat, METEO *met, PAR *par, TIMES *time, TOPO *top, SOIL *sl);

void print_run_average(SOIL *sl, TOPO *top, PAR *par);

void init_run(SOIL *sl, PAR *par);

void end_period_1D(SOIL *sl, TOPO *top, PAR *par);

void change_grid(long previous_sim, long next_sim, PAR *par, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet);
#endif
