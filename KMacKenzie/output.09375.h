
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


#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "recovery.h"
#include "keywords_file.h"
#include "struct.geotop.09375.h"
#include "vegetation.h"
#include "pedo.funct.h"
#include "geo_statistic.09375.h"
#include "networks.h"
#include "rw_maps.h"
#include "constant.h"
#include "extensions.h"
#include "times.h"
#include "energy.balance.h"
#include "tabs.h"

#include "frost_table.h"


/* These variables have a global scope */
extern double wt0_basin; /*mean intercepted precipitation [mm] in the previous output-basin Dt*/
extern double Ssup;       /*supercial Storage of water in all the basin [mm]*/
extern double Ssub;      /*subsuperficial Storage of water in all the basin [mm]*/
extern double Rout;      /*sum of the output flows from the last output-basin for unit of area[mm]*/
extern double R_G;
extern double S_ch0;     /*wat in the channel at the start of a basin-output step-time*/
extern double S_ch1;     /*wat in the channel at the end of a basin-output step-time*/
extern double Qsub_ch, Qsup_ch ,Q_G; /*averaged output flows*/
extern double SWE_previous;
extern double GWE_previous;
extern double Smelt;   /*Snow melt [mm] during the time interval*/
extern double Ssubl;   /*Snow sublimation [mm] during the time interval*/
extern double Sevap;   /*Snow evaporation [mm] during the time interval*/
extern double Gmelt;   /*Glacier melt [mm] during the time interval*/
extern double Gsubl;   /*Glacier sublimation [mm] during the time interval*/
extern double Gevap;   /*Glacier evaporation [mm] during the time interval*/
extern double Smelt_previous;
extern double Ssubl_previous;
extern double Sevap_previous;
extern double Gmelt_previous;
extern double Gsubl_previous;
extern double Gevap_previous;
extern long isavings;
/* End variables that have a global scope */


/****************************************************************************************************/
/* write_output: stampa i dati degli output variabili nel tempo                                     */
/****************************************************************************************************/
void write_output(TIMES *times,WATER *wat,CHANNEL *cnet,PAR *par,TOPO *top,LAND *land,SOIL *sl,ENERGY *egy,SNOW *snow,GLACIER *glac,METEO *met);



/****************************************************************************************************/
/* All the structs and substructs of the simulation are deallocated:                                */
/****************************************************************************************************/
void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times);

void dealloc_meteostations(METEO_STATIONS *st);


/*==================================================================================================================*/
void write_date(FILE *f, long day, long month, long year, long hour, long min);

void plot(char *name, long JD, long y, long i, DOUBLEMATRIX *M, short format);

void write_init_condit(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac);

void write_soil_output(long n, long i, double t, double dt, long y0, double JD0, LONGMATRIX *rc, SOIL *sl, double psimin, double Esoil);

void initmatrix(double val, DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue);

void initlongmatrix(long val, LONGMATRIX *destination, DOUBLEMATRIX *origin, double novalue);

void inittensor(double val, DOUBLETENSOR *destination, DOUBLEMATRIX *origin, double novalue);

double interp_value(double E, double N, DOUBLEMATRIX *M, DOUBLEMATRIX *Z);
long row1(double N, long nrows, long i, T_INIT *UV);
long col1(double E, long ncols, long i, T_INIT *UV);

void write_output_superfast(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac, METEO *met);

void write_supertensor(long r, long c,short sy,double Dt_output, short n,double t, double dt, long y0, double JD0, char* file, double* Dz, TIMES *times, double *** out);

#endif
