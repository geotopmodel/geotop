
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
short existing_file_text(char *name);

char *namefile_i(char *name, long i);
char *namefile_i_we(char *name, long i);
char *namefile_i_we2(char *name, long i);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//READ subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

DOUBLEMATRIX *read_map(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref);

DOUBLEMATRIX *read_mapseries(long i, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref);

DOUBLETENSOR *read_tensor(long nl, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref);

DOUBLETENSOR *read_maptensor(long i, long lmax, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//WRITE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_map(char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV);

void write_mapseries(long i, char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV);

void write_tensorseries(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV);

void write_tensorseries_bis(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV);

void write_tensorseries2(long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV);