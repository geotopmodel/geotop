
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 1.225-15
 
 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

void windtrans_snow(SNOW *snow, METEO *met, WATER *wat, LAND *land, TOPO *top, PAR *par, double t0);

void set_inhomogeneous_fetch(SNOW *snow, METEO *met, LAND *land, PAR *par, TOPO *top, short *yes);

void set_windtrans_snow(double Dt, double t, SNOW *snow, METEO *met, LAND *land, PAR *par, FILE *f);

void print_windtrans_snow(double Dt, SNOW *snow, PAR *par, TOPO *top, METEO *met, DOUBLEMATRIX *LC);

void extend_topography(DOUBLEMATRIX *M, double novalue);

void extend_topography_row(DOUBLEMATRIX *M, double novalue);

void extend_topography_column(DOUBLEMATRIX *M, double novalue);

void find_the_nearest(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

void find_the_nearest_row(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

void find_the_nearest_column(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

short no_novalue(long r, long c, DOUBLEMATRIX *M, double novalue, long *rr, long *cc);

void set_no_value(DOUBLEMATRIX *M, DOUBLEMATRIX *N, double undef);

void wind_packing(SNOW *snow, PAR *par, long r, long c, double Dt);
