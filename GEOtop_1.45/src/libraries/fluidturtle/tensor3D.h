#ifndef TENSOR3D_H
#define TENSOR3D_H
#include "turtle.h"
void initialize_doubletensor(DOUBLETENSOR *L,double sign);
void copy_doubletensor(DOUBLETENSOR *origin,DOUBLETENSOR *destination);
DOUBLETENSOR *new_doubletensor_flexlayer(long ndl,long ndh,long nrh,long nch);
DOUBLETENSOR *new_doubletensor0(long ndh,long nrh,long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long ncl, long ndl);
void free_doubletensor( DOUBLETENSOR *m);
#endif
