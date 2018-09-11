#include "turtle.h"
#include "t_datamanipulation.h"
//#include "t_statistics.h"


/*--------------------------------------------------------------------------*/

void copy_doublematrix(Matrix<double> *origin, Matrix<double> *destination) {

  long i, j;

//  if (origin == nullptr || destination == nullptr || origin->co == nullptr
//      || destination->co == nullptr) {
//
//    t_error("A matrix was not allocated");
//
//  } else if (origin->isdynamic != 1 || destination->isdynamic != 1 || origin->nrh < 1
//             || destination->nrh < 1 || origin->nch < 1 || destination->nch < 1) {
//
//    t_error("A matrix was not allocated properly");
//
//  } else if (origin->nrh != destination->nrh
//             || origin->nch != destination->nch) {
//
//    t_error("The matrixes do not have the same dimensions");
//
//  }

  for (i = origin->nrl; i <= origin->nrh; i++) {
    for (j = origin->ncl; j <= origin->nch; j++) {

      (*destination)(i,j) = (*origin)(i,j);

    }
  }
}




