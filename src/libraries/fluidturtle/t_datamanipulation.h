#define RATIO 1000

/**

Name: initialize_

Version: 1.0

Synopsis:
void initialize_longmatrix(LONGMATRIX *, long);


Description:  It initialize the matrix or vector with  the specified value

Inputs:  1) The pointer to the structure to initalize; 2) the value used for initialization

Examples:  Variogram.c


Authors & Date: Riccardo Rigon, October 1997.

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c


*/


void initialize_longvector(LONGVECTOR *, long);
void initialize_shortvector(SHORTVECTOR *,short );

void initialize_longmatrix(LONGMATRIX *, long);
void initialize_shortmatrix(SHORTMATRIX *,short );

void initialize_doublematrix(DOUBLEMATRIX *,double );


void copy_shortmatrix(SHORTMATRIX *,SHORTMATRIX *);

void copy_doublematrix(DOUBLEMATRIX *,DOUBLEMATRIX *);

