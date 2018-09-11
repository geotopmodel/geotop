#include "turtle.h"


/*------------------------------------------------------------------------

Standard NR allocation routines for vectors:
for documentation refer to turtle.h

--------------------------------------------------------------------------*/

float *vector(long nl, long nh)

/* allocate a float vector  with subscript range v[nl ....nh] */
{

  float *v;

  v = (float *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(float)));

  if (!v) t_error("allocation failure in fvector()");

  return v - nl + NR_END;

}

/*------------------------------------------------------------------------

Standard NR allocation routines for matrixes:
for documentation refer to turtle.h

--------------------------------------------------------------------------*/

float **matrix(long nrl, long nrh, long ncl, long nch)

/* Allocate a float matrix  with subscript range v[nrl ....nrh][ncl ....nch] */
{

  long i, rows = nrh - nrl + 1, cols = nch - ncl + 1;

  float **m;

  /* Allocate array of pointers to arrays */

  m = (float **) malloc((size_t) ((rows + NR_END) * sizeof(float *)));

  if (!m) t_error("Allocation failure 1 in matrix()");

  m += NR_END;

  m -= nrl;

  /* Allocate rows and set pointers to them   */

  m[nrl] = (float *) malloc((size_t) ((rows * cols + NR_END) * sizeof(float)));

  if (!m) t_error("Allocation failure 2 in matrix()");

  m[nrl] += NR_END;

  m[nrl] -= ncl;

  /* Returns pointer to array of pointers to rows */

  for (i = nrl + 1; i <= nrh; i++) {

    m[i] = m[i - 1] + cols;

  }

  return m;

}

