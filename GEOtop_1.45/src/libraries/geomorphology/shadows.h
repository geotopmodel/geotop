#ifndef SHADOWS_H
#define SHADOWS_H
#include "../fluidturtle/turtle.h"

/*----------------------------------------------------------------------------------------------------------*/
/*Used to calculate shadow matrix*/
/*----------------------------------------------------------------------------------------------------------*/


void Orizzonte1(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
void Orizzonte2(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
void Orizzonte3(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
void Orizzonte4(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
void Orizzonte5(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
void Orizzonte6(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
void Orizzonte7(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
void Orizzonte8(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue);
				
#endif
