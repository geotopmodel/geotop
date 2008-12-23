
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
    
    
    
#define end_vector -9.0E+5
#define end_vector_long 999999
/*----------------------------------------------------------------------------------------------------------*/
char **readline_textarray(FILE *f, long offset);

/*----------------------------------------------------------------------------------------------------------*/
void readline_array(FILE *f, double *a, long offset, long ncol, double ndef, short *endoffile);

/*----------------------------------------------------------------------------------------------------------*/
double decod(char *ch, long n, double ndef);

/*----------------------------------------------------------------------------------------------------------*/
long dim1(double *a);
long dim1l(long *a);
long dim2(double **a);
long dim_string(char *a);
long dim_vect_strings(char **a);

/*----------------------------------------------------------------------------------------------------------*/

double **alloc2(long n, long m);
double *alloc1(long n);
long **alloc_long2(long n);
long *alloc_long1(long n);


/*----------------------------------------------------------------------------------------------------------*/

short compare_strings(char *a, char *b);