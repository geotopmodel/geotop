
/* TO DO
 * \author Emanuele Cordano
 */

#ifdef USE_NETCDF

#ifndef NCGT_UTILITIES_H
#define NCGT_UTILITIES_H

#include "../libraries/fluidturtle/turtle.h"
//#include <netcdf.h>
#include "gt_utilities.h"

//char *copy_stringnames(const char *);

int rotate180_y_doublematrix(DOUBLEMATRIX *);
int rotate180_y_doublematrix(GeoMatrix<double>& M);
int rotate180_y_doubletensor(DOUBLETENSOR *);
int rotate180_y_doubletensor(GeoTensor<double>& M);
int rotate180_y_floatmatrix(FLOATMATRIX *);


int rotate180_y_longmatrix(LONGMATRIX *);
int rotate180_y_intmatrix(INTMATRIX *);
int rotate180_y_shortmatrix(SHORTMATRIX *);

int invert_order_doublevector(DOUBLEVECTOR *);
int invert_order_longvector(LONGVECTOR *);

#endif
#endif
