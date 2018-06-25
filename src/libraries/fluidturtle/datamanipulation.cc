#include "turtle.h"
#include "t_datamanipulation.h"
//#include "t_statistics.h"



/*-------------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
void initialize_shortmatrix(SHORTMATRIX *L, short sign) {

  long i, j;

  if (L != NULL) {
    if (L->isdynamic == 1) {
      for (i = 1; i <= L->nrh; i++) {
        for (j = 1; j <= L->nch; j++) {
          L->co[i][j] = sign;
        }
      }
    } else {
      t_error("This matrix was no properly allocated");
    }
  } else {
    t_error("A null matrix was addressed");
  }
}


/*---------------------------------------------------------------------------*/
void initialize_longmatrix(LONGMATRIX *L, long sign) {

  long i, j;

  if (L != NULL) {
    if (L->isdynamic == 1) {
      for (i = 1; i <= L->nrh; i++) {
        for (j = 1; j <= L->nch; j++) {
          L->co[i][j] = sign;
        }
      }
    } else {
      t_error("This matrix was no properly allocated");
    }
  } else {
    t_error("A null matrix was addressed");
  }
}


/*---------------------------------------------------------------------------*/
void initialize_doublematrix(DOUBLEMATRIX *L, double sign) {

  long i, j;

  if (L != NULL) {
    if (L->isdynamic == 1) {
      for (i = L->nrl; i <= L->nrh; i++) {
        for (j = L->ncl; j <= L->nch; j++) {
          L->co[i][j] = sign;
        }
      }
    } else {
      t_error("This matrix was no properly allocated");
    }
  } else {
    t_error("A null matrix was addressed");
  }
}


/*--------------------------------------------------------------------------*/

void copy_shortmatrix(SHORTMATRIX *origin, SHORTMATRIX *destination) {

  long i, j;

  if (origin == NULL || destination == NULL || origin->co == NULL
      || destination->co == NULL) {

    t_error("A matrix was not allocated");

  } else if (origin->isdynamic != 1 || destination->isdynamic != 1 || origin->nrh < 1
             || destination->nrh < 1 || origin->nch < 1 || destination->nch < 1) {

    t_error("A matrix was not allocated properly");

  } else if (origin->nrh != destination->nrh
             || origin->nch != destination->nch) {

    t_error("The matrixes do not have the same dimensions");

  }

  for (i = 1; i <= origin->nrh; i++) {
    for (j = 1; j <= origin->nch; j++) {

      destination->co[i][j] = origin->co[i][j];

    }
  }

}


/*--------------------------------------------------------------------------*/

void copy_doublematrix(DOUBLEMATRIX *origin, DOUBLEMATRIX *destination) {

  long i, j;

  if (origin == NULL || destination == NULL || origin->co == NULL
      || destination->co == NULL) {

    t_error("A matrix was not allocated");

  } else if (origin->isdynamic != 1 || destination->isdynamic != 1 || origin->nrh < 1
             || destination->nrh < 1 || origin->nch < 1 || destination->nch < 1) {

    t_error("A matrix was not allocated properly");

  } else if (origin->nrh != destination->nrh
             || origin->nch != destination->nch) {

    t_error("The matrixes do not have the same dimensions");

  }

  for (i = origin->nrl; i <= origin->nrh; i++) {
    for (j = origin->ncl; j <= origin->nch; j++) {

      destination->co[i][j] = origin->co[i][j];

    }
  }

}




