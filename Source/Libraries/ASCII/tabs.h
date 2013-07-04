
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