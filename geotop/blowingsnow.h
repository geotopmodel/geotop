
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef BLOWINGSNOW_H
#define BLOWINGSNOW_H

#include "input.h"
#include "constants.h"
#include "struct.geotop.h"
#include "snow.h"
#include "PBSM.h"
#include "meteo.h"
#include "vegetation.h"
#include "energy.balance.h"
#include "meteodata.h"

extern long number_novalue, number_absent;

//extern T_INIT *UV;
extern TInit *UV;

//extern char *logfile;
extern std::string logfile;

extern long Nl, Nr, Nc;
//extern char *WORKING_DIRECTORY;
extern std::string WORKING_DIRECTORY;
extern double **odp;

extern long i_sim;

//void windtrans_snow(SNOW *snow, METEO *met, WATER *wat, LAND *land, TOPO *top, PAR *par, double t0);
void windtrans_snow(Snow *snow, Meteo *met, Water *wat, Land *land, Topo *top, Par *par, double t0);

//void set_inhomogeneous_fetch(SNOW *snow, METEO *met, LAND *land, PAR *par, TOPO *top, short *yes);
void set_inhomogeneous_fetch(Snow *snow, Meteo *met, Land *land, Par *par, Topo *top, short *yes);

//void set_windtrans_snow(double Dt, double t, SNOW *snow, METEO *met, LAND *land, PAR *par, FILE *f);
void set_windtrans_snow(double Dt, double t, Snow *snow, Meteo *met, Land *land, Par *par, FILE *f);

//void print_windtrans_snow(double Dt, SNOW *snow, PAR *par, TOPO *top, METEO *met, DOUBLEMATRIX *LC);
void print_windtrans_snow(double Dt, Snow *snow, Par *par, Topo *top, Meteo *met, GeoMatrix<double>& LC);

void extend_topography(DOUBLEMATRIX *M, double novalue);

void extend_topography_row(DOUBLEMATRIX *M, double novalue);

void extend_topography_column(DOUBLEMATRIX *M, double novalue);

void find_the_nearest(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

void find_the_nearest_row(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

void find_the_nearest_column(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

short no_novalue(long r, long c, DOUBLEMATRIX *M, double novalue, long *rr, long *cc);

void set_no_value(DOUBLEMATRIX *M, DOUBLEMATRIX *N, double undef);

#endif
