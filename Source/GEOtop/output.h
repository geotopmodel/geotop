
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
    
    
    
/****************************************************************************************************/
/* write_output: stampa i dati degli output variabili nel tempo                                     */
/****************************************************************************************************/
void write_output(TIMES *times,WATER *wat,CHANNEL *cnet,PAR *par,TOPO *top,LAND *land,SOIL *sl,ENERGY *egy,SNOW *snow,GLACIER *glac,METEO *met);



/****************************************************************************************************/
/* All the structs and substructs of the simulation are deallocated:                                */
/****************************************************************************************************/
void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met);

void dealloc_meteostations(METEO_STATIONS *st);


/*==================================================================================================================*/
void plot(char *name, long JD, long y, long i, DOUBLEMATRIX *M, short format);

void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac);

void write_soil_output(long n, long i, double t, double dt, long y0, double JD0, LONGMATRIX *rc, SOIL *sl, double psimin);

void initmatrix(double val, DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue);

void initlongmatrix(long val, LONGMATRIX *destination, DOUBLEMATRIX *origin, double novalue);

void inittensor(double val, DOUBLETENSOR *destination, DOUBLEMATRIX *origin, double novalue);

double interp_value(double E, double N, DOUBLEMATRIX *M, DOUBLEMATRIX *Z);
long row1(double N, long nrows, long i, T_INIT *UV);
long col1(double E, long ncols, long i, T_INIT *UV);
