
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion KMackenzie

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 KMackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOtop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/




/****************************************************************************************************/
/* write_output: stampa i dati degli output variabili nel tempo                                     */
/****************************************************************************************************/
void write_output(TIMES *times,WATER *wat,CHANNEL *cnet,PAR *par,TOPO *top,LAND *land,SOIL *sl,ENERGY *egy,SNOW *snow,GLACIER *glac,METEO *met);



/****************************************************************************************************/
/* All the structs and substructs of the simulation are deallocated:                                */
/****************************************************************************************************/
void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, LISTON *liston);

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
