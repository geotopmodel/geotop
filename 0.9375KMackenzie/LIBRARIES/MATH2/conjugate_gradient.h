/*
 * coniugate_gradient.h
 *
 *  Created on: Nov 7, 2008
 *      Author: ecor
 */



long conjugate_gradient_search(long icnt, double epsilon,  DOUBLEVECTOR *x, DOUBLEVECTOR *b, int (* funz)(DOUBLEVECTOR *y,DOUBLEVECTOR *x));
int linear_comb_doublevector(DOUBLEVECTOR *result,DOUBLEVECTOR *a, DOUBLEVECTOR *b, double ca, double cb);
