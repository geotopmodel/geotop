
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
#ifndef RW_MAPS_H
#define RW_MAPS_H
#include "../fluidturtle/turtle.h"
#include "import_ascii.h"
#include "write_ascii.h"
#include "extensions.h"
#include "../fluidturtle/t_datamanipulation.h"
#include "../fluidturtle/t_utilities.h"
#include "../fluidturtle/tensor3D.h"
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//BASE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

SHORTMATRIX *copyshort_doublematrix(DOUBLEMATRIX *M);

LONGMATRIX *copylong_doublematrix(DOUBLEMATRIX *M);

DOUBLEMATRIX *copydouble_shortmatrix(SHORTMATRIX *S);

DOUBLEMATRIX *copydouble_longmatrix(LONGMATRIX *L);

DOUBLEMATRIX *copydoublematrix_const(double c0, DOUBLEMATRIX *Mref, double NOVALUE);

DOUBLEMATRIX *multiplydoublematrix(double f, DOUBLEMATRIX *Mref, double NOVALUE);

void build_doubletensor(DOUBLETENSOR *T, DOUBLEMATRIX *M, long l);

DOUBLEMATRIX *extract_doublematrix(DOUBLETENSOR *T, long l);

DOUBLEMATRIX *extract_fromtensor(DOUBLETENSOR *T, long l);

DOUBLETENSOR *build_frommatrix(DOUBLEMATRIX *M, long l, long lmax);

void write_frommatrix(long l, DOUBLEMATRIX *M, DOUBLETENSOR *T);

void fmultiplydoublematrix(DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double f, double novalue);

void assignnovalue(DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//UTILITITY subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_suffix(char *suffix, long i, short start);

short existing_file(char *name);
short existing_file_wext(char *name, char *extension);
short existing_file_woext(char *name);

char *namefile_i(char *name, long i);
char *namefile_i_we(char *name, long i);
char *namefile_i_we2(char *name, long i);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//READ subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

DOUBLEMATRIX *read_map(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref, double no_value);

DOUBLEVECTOR *read_map_vector(short type, char *namefile, DOUBLEMATRIX *mask, T_INIT *grid, double no_value, LONGMATRIX *rc);

DOUBLEMATRIX *read_mapseries(long i, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref, double no_value);

DOUBLETENSOR *read_tensor(long nl, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref, double no_value);

DOUBLETENSOR *read_maptensor(long i, long lmax, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref, double no_value);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//WRITE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_map(char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV, long novalue);

void write_map_vector(char *filename, short type, short format, DOUBLEVECTOR *V, T_INIT *UV, long novalue, long **j, long nr, long nc);

void write_mapseries(long i, char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV, long novalue);

void write_tensorseries(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);

void write_tensorseries_vector(short a, long l, long i, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc);

void rename_tensorseries(short a, long l, long i, char *filename);

void rename_map(char *filename);

void write_tensorseries2(char *suf, long l, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);

void write_tensorseries2_vector(char *suf, long l, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc);

void write_tensorseries3(char *suffix, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);

void write_tensorseries3_vector(char *suffix, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc);
#endif
