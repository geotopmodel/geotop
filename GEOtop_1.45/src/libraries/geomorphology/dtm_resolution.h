#ifndef DTM_RESOLUTION_H
#define DTM_RESOLUTION_H
#include "../fluidturtle/turtle.h"

void reduce_resolution(long n, DOUBLEMATRIX *DTM1, DOUBLEMATRIX *DTM2, T_INIT *UV1, T_INIT *UV2);

void amplify_resolution(long n, long Rr, long Rc, DOUBLEMATRIX *DTM1, DOUBLEMATRIX *DTM2, T_INIT *UV1, T_INIT *UV2);
#endif
