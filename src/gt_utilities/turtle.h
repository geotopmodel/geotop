/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#ifndef TURTLE_H
#define TURTLE_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "../geotop/datastructs.h"
//#define element co //to inhertit old nomunclature for FluidTurle

#define isDynamic 1         /*This number will be used to mark a dynamic allocates
                              quantity */
#define MAXBUFFERSIZE 50000 /* Many functions uses a buffer to allocate data. A
                             safe habit is to limit its size. */

#define MAXLABELSIZE 1024    /* The same as above */

#define BUFFERINCREMENT 256 /*Buffers memory is increases when necessary of this
							 storage capability till MAXBUFFERSIZE*/

#define LABELINCREMENT 64   /*The same as above */

#define NR_END    1         /* Required by the Numerical Recipes' allocation routines
							*/

#define FREE_ARG  char*     /* The same as above */

#define OK 1

#define NOK -1

#define WOK 0

#define MAX_KEYWORD_LENGTH 10

#define PRINT 1
#define NOPRINT 0         /*With the above PRINT is used to set a printint option
                           in some functions */

#define NL 1              /* Numerical Recipes allocation routines allow to have
							 arbitrary  subscripts for vector and matrixes. The
							 fluid turtle library restrict this freedom by setting
						     their  lower value to NL  */
						     
						     
#define TEST 1

#define NOTEST 0

#define SQRT2  1.414213562373095

/**-------------------------------------------------------------------
VECTORS: essentially the same types as in Numerical Recipes, except
that char and long are not unsigned as there. We have vector for short, 
int,long , char, float, double.  
---------------------------------------------------------------------*/

/**-------------------------------------------------------------------
MATRIXES: The same types as for vectors except for char:
SHORTMATRIX, INTMATRIX, LONGMATRIX, FLOATMATRIX, DOUBLEMATRIX
---------------------------------------------------------------------*/



/**-------------------------------------------------------------------
BINS: A BIN is a collection of vectors of different lengths. 
We have BINS for short, long and double, a BIN of char vectors is 
here called STRINGBIN. Thus, we have:
SHORTBIN, INTBIN, LONGBIN, DOUBLEBIN, STRING
---------------------------------------------------------------------*/



/* For the Numerical Recipes Random Generator  ran1 */
//#define PI 3.141592654
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349
#define NR_END 1
/* For the Numerical Recipes Random Generator  ran2 */

#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/IM1)
#define IMM1 (IM1-1)
#define IAA1 40014
#define IAA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2  3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX  (1.0 - EPS)

#define THRESH 0
#define ITOL 3
#define TOL 0.00001
#define ITMAX 1000
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define EPS1 1.0e-14
#define TINY 1.0E-20


#endif
