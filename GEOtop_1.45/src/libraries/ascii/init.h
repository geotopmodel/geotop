#ifndef INIT_H
#define INIT_H
#include "../fluidturtle/turtle.h"
void initmatrix(double val, DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue);

void initlongmatrix(long val, LONGMATRIX *destination, DOUBLEMATRIX *origin, double novalue);

void inittensor(double val, DOUBLETENSOR *destination, DOUBLEMATRIX *origin, double novalue);
#endif
