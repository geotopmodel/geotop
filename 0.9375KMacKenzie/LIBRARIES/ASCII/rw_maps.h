
/* STATEMENT:

ASCII-GEOtop LIBRARIES

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon

 LICENSE:

 This file is part of ASCII-GEOtop LIBRARIES
 ASCII-GEOtop LIBRARIES is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
    
    
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