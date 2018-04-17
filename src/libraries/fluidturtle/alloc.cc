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


/*-----------------------------------------------------------------------*/


long *lvector(long nl, long nh)


/* allocate an integer vector  with subscript range v[nl ....nh] */


{


  long *v;


  v = (long *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(long)));


  if (!v) t_error("allocation failure in fvector()");


  return v - nl + NR_END;


}


/*-----------------------------------------------------------------------*/


short *svector(long nl, long nh)


/* allocate an integer vector  with subscript range v[nl ....nh] */


{


  short *v;


  v = (short *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(short)));


  if (!v) t_error("allocation failure in fvector()");


  return v - nl + NR_END;


}


/*------------------------------------------------------------------------


Standard NR allocation routines for matrixes:


for documentation refer to turtle.h


--------------------------------------------------------------------------*/





float **matrix(long nrl, long nrh, long ncl, long nch)


/* Allocate a float matrix  with subscript range v[nrl ....nrh][ncl ....nrh] */


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


/*-----------------------------------------------------------------------*/


short **smatrix(long nrl, long nrh, long ncl, long nch)


/* Allocate an int matrix  with subscript range v[nrl ....nrh][ncl ....nrh] */


{


  long i, rows = nrh - nrl + 1, cols = nch - ncl + 1;


  short **m;





  /* Allocate array of pointers to arrays */





  m = (short **) malloc((size_t) ((rows + NR_END) * sizeof(short *)));


  if (!m) t_error("Allocation failure 1 in matrix()");


  m += NR_END;


  m -= nrl;





  /* Allocate rows and set pointers to them   */





  m[nrl] = (short *) malloc((size_t) ((rows * cols + NR_END) * sizeof(short)));


  if (!m) t_error("Allocation failure 2 in matrix()");


  m[nrl] += NR_END;


  m[nrl] -= ncl;


  /* Returns pointer to array of pointers to rows */


  for (i = nrl + 1; i <= nrh; i++) {


    m[i] = m[i - 1] + cols;


  }


  return m;


}


/*-----------------------------------------------------------------------*/


long **lmatrix(long nrl, long nrh, long ncl, long nch)


/* Allocate a long matrix  with subscript range v[nrl ....nrh][ncl ....nrh] */


{


  long i, rows = nrh - nrl + 1, cols = nch - ncl + 1;


  long **m;





  /* Allocate array of pointers to arrays */





  m = (long **) malloc((size_t) ((rows + NR_END) * sizeof(long *)));


  if (!m) t_error("Allocation failure 1 in matrix()");


  m += NR_END;


  m -= nrl;





  /* Allocate rows and set pointers to them   */





  m[nrl] = (long *) malloc((size_t) ((rows * cols + NR_END) * sizeof(long)));


  if (!m) t_error("Allocation failure 2 in matrix()");


  m[nrl] += NR_END;


  m[nrl] -= ncl;


  /* Returns pointer to array of pointers to rows */


  for (i = nrl + 1; i <= nrh; i++) {


    m[i] = m[i - 1] + cols;


  }


  return m;


}


/*-----------------------------------------------------------------------*/


double **dmatrix(long nrl, long nrh, long ncl, long nch)


/* Allocate a long matrix  with subscript range v[nrl ....nrh][ncl ....nrh] */


{


  long i, rows = nrh - nrl + 1, cols = nch - ncl + 1;


  double **m;





  /* Allocate array of pointers to arrays */





  m = (double **) malloc((size_t) ((rows + NR_END) * sizeof(double *)));


  if (!m) t_error("Allocation failure 1 in matrix()");


  m += NR_END;


  m -= nrl;





  /* Allocate rows and set pointers to them   */





  m[nrl] = (double *) malloc((size_t) ((rows * cols + NR_END) * sizeof(double)));


  if (!m) t_error("Allocation failure 2 in matrix()");


  m[nrl] += NR_END;


  m[nrl] -= ncl;


  /* Returns pointer to array of pointers to rows */


  for (i = nrl + 1; i <= nrh; i++) {


    m[i] = m[i - 1] + cols;


  }


  return m;


}








/*------------------------------------------------------------------------


Wrappers for vectors and matrixes


--------------------------------------------------------------------------*/





/*-----------------------------------------------------------------------*/





SHORTVECTOR *new_shortvector(long nh) {


  SHORTVECTOR *m;


  m = (SHORTVECTOR *) malloc(sizeof(SHORTVECTOR));


  if (!m) t_error("allocation failure in SHORTVECTOR()");


  m->isdynamic = isDynamic;


  m->nl = NL;


  m->nh = nh;


  m->co = svector(1, nh);


  return m;


}


/*-----------------------------------------------------------------------*/


FLOATVECTOR *new_floatvector(long nh) {


  FLOATVECTOR *m;


  m = (FLOATVECTOR *) malloc(sizeof(FLOATVECTOR));


  if (!m) t_error("allocation failure in FLOATVECTOR()");


  m->isdynamic = isDynamic;


  m->nl = NL;


  m->nh = nh;


  m->co = vector(m->nl, nh);


  return m;


}


/*-----------------------------------------------------------------------*/


LONGVECTOR *new_longvector(long nh) {


  LONGVECTOR *m;


  m = (LONGVECTOR *) malloc(sizeof(LONGVECTOR));


  if (!m) t_error("allocation failure in LONGVECTOR()");


  m->isdynamic = isDynamic;


  m->nl = NL;


  m->nh = nh;


  m->co = lvector(m->nl, nh);


  return m;


}


/*-----------------------------------------------------------------------*/


SHORTMATRIX *new_shortmatrix(long nrh, long nch) {


  SHORTMATRIX *m;


  m = (SHORTMATRIX *) malloc(sizeof(SHORTMATRIX));


  if (!m) t_error("allocation failure in SHORTMATRIX()");


  m->isdynamic = isDynamic;


  m->nrl = NL;


  m->nrh = nrh;


  m->ncl = NL;


  m->nch = nch;


  m->co = smatrix(1, nrh, 1, nch);


  return m;


}


/*-----------------------------------------------------------------------*/


LONGMATRIX *new_longmatrix(long nrh, long nch) {


  LONGMATRIX *m;


  m = (LONGMATRIX *) malloc(sizeof(LONGMATRIX));


  if (!m) t_error("allocation failure in LONGMATRIX()");


  m->isdynamic = isDynamic;


  m->nrl = NL;


  m->nrh = nrh;


  m->ncl = NL;


  m->nch = nch;


  m->co = lmatrix(1, nrh, 1, nch);


  return m;


}


/*-----------------------------------------------------------------------*/


DOUBLEMATRIX *new_doublematrix(long nrh, long nch) {


  DOUBLEMATRIX *m;


  m = (DOUBLEMATRIX *) malloc(sizeof(DOUBLEMATRIX));


  if (!m) t_error("allocation failure in new_doublematrix()");


  m->isdynamic = isDynamic;


  m->nrl = NL;


  m->nrh = nrh;


  m->ncl = NL;


  m->nch = nch;


  m->co = dmatrix(1, nrh, 1, nch);


  return m;


}


/*-----------------------------------------------------------------------*/


DOUBLEMATRIX *new_doublematrix0_(long nrh, long nch) {


  DOUBLEMATRIX *m;


  m = (DOUBLEMATRIX *) malloc(sizeof(DOUBLEMATRIX));


  if (!m) t_error("allocation failure in new_doublematrix()");


  m->isdynamic = isDynamic;


  m->nrl = 0;


  m->nrh = nrh;


  m->ncl = NL;


  m->nch = nch;


  m->co = dmatrix(m->nrl, m->nrh, m->ncl, m->nch);


  return m;


}

/*-----------------------------------------------------------------------*/


DOUBLEMATRIX *new_doublematrix_0(long nrh, long nch) {


  DOUBLEMATRIX *m;


  m = (DOUBLEMATRIX *) malloc(sizeof(DOUBLEMATRIX));


  if (!m) t_error("allocation failure in new_doublematrix()");


  m->isdynamic = isDynamic;


  m->nrl = NL;


  m->nrh = nrh;


  m->ncl = 0;


  m->nch = nch;


  m->co = dmatrix(m->nrl, m->nrh, m->ncl, m->nch);


  return m;


}


/*------------------------------------------------------------------------


Memory deallocation routines


--------------------------------------------------------------------------*/





void free_svector(short *v, long nl) {


  free((FREE_ARG) (v + nl - NR_END));


}


/*-----------------------------------------------------------------------*/


void free_vector(float *v, long nl) {


  free((FREE_ARG) (v + nl - NR_END));


}


/*-----------------------------------------------------------------------*/


void free_lvector(long *v, long nl) {


  free((FREE_ARG) (v + nl - NR_END));


}





/*-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/





void free_smatrix(short **m, long nrl, long ncl) {


  free((FREE_ARG) (m[nrl] + ncl - NR_END));


  free((FREE_ARG) (m + nrl - NR_END));


}


/*-----------------------------------------------------------------------*/


void free_lmatrix(long **m, long nrl, long ncl) {


  free((FREE_ARG) (m[nrl] + ncl - NR_END));


  free((FREE_ARG) (m + nrl - NR_END));


}


/*-----------------------------------------------------------------------*/


void free_dmatrix(double **m, long nrl, long ncl) {


  free((FREE_ARG) (m[nrl] + ncl - NR_END));


  free((FREE_ARG) (m + nrl - NR_END));


}


/*-----------------------------------------------------------------------*/





void free_shortvector(SHORTVECTOR *v) {


  if (v == NULL || v->co == NULL) {


    t_error("This shortvector was never allocated");


  } else if (v->isdynamic == 1) {


    free_svector(v->co, NL);


    v->isdynamic = v->nl = v->nh = -1;


    free(v);


    return;


  } else {


    printf("\nWarning::An attemp was made to free a non dynamic vector 1\n");


  }


}


/*-----------------------------------------------------------------------*/





void free_floatvector(FLOATVECTOR *v) {


  if (v == NULL || v->co == NULL) {


    t_error("This floatvector was never allocated");


  } else if (v->isdynamic == 1) {


    free_vector(v->co, NL);


    v->isdynamic = v->nl = v->nh = -1;


    free(v);


    return;


  } else {


    printf("\nWarning::An attemp was made to free a non dynamic vector 3\n");


  }


}


/*-----------------------------------------------------------------------*/





void free_longvector(LONGVECTOR *v) {


  if (v == NULL || v->co == NULL) {


    t_error("This longvector was never allocated");


  } else if (v->isdynamic == 1) {


    free_lvector(v->co, NL);


    v->isdynamic = v->nl = v->nh = -1;


    free(v);


    return;


  } else {


    printf("\nWarning::An attemp was made to free a non dynamic vector 4\n");


  }


}


/*-----------------------------------------------------------------------*/





void free_shortmatrix(SHORTMATRIX *m) {


  if (m == NULL || m->co == NULL) {


    t_error("This matrix was never allocated");


  } else if (m->isdynamic == 1) {


    free_smatrix(m->co, NL, NL);


    m->isdynamic = m->nrl = m->ncl = m->nrh = m->nch = -1;


    free(m);


    return;


  } else {


    printf("\nWarning::An attemp was made to free a non dynamic matrix\n");


  }


}


/*-----------------------------------------------------------------------*/





void free_longmatrix(LONGMATRIX *m) {


  if (m == NULL || m->co == NULL) {


    t_error("This matrix was never allocated");


  } else if (m->isdynamic == 1) {


    free_lmatrix(m->co, NL, NL);


    m->isdynamic = m->nrl = m->ncl = m->nrh = m->nch = -1;


    free(m);


    return;


  } else {


    printf("\nWarning::An attemp was made to free a non dynamic matrix\n");


  }


}


/*-----------------------------------------------------------------------*/





void free_doublematrix(DOUBLEMATRIX *m) {

  if (m == NULL || m->co == NULL) {
    t_error("This matrix was never allocated");
  } else if (m->isdynamic == 1) {
    free_dmatrix(m->co, m->nrl, m->ncl);
    m->isdynamic = m->nrl = m->ncl = m->nrh = m->nch = -1;
    free(m);
    return;
  } else {
    printf("\nWarning::An attemp was made to free a non dynamic matrix\n");
  }


}








