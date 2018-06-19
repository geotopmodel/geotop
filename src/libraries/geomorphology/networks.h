#define MAXELEVATION 10000


/**

Name:  randomflow

Synopsis:  void randomflow(long ,long , SHORTMATRIX *);


Version:  0.9

Description: given the matrix of flowing directions and a position
in the matrix it returns a new flowing direction selected at random

Authors & date: Riccardo Rigon, November 1997


Inputs: 1) 2) the position; 3) The pointer to the flowing direction matrix

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997




*/


void initialize_longmatrix(LONGMATRIX *,long );


