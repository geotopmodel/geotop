#ifndef TENSOR3D_H
#define TENSOR3D_H

#include "turtle.h"
#include "tensors3D.h"
#include "../../geotop/constants.h"
#include <vector>
#include "../../geotop/struct.geotop.h"
#include "../../geotop/datastructs.h"

void initialize_doubletensor(DOUBLETENSOR *L,double sign);
void copy_doubletensor(DOUBLETENSOR *origin,DOUBLETENSOR *destination);
DOUBLETENSOR *new_doubletensor_flexlayer(long ndl,long ndh,long nrh,long nch);
DOUBLETENSOR *new_doubletensor0(long ndh,long nrh,long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long ncl, long ndl);
void free_doubletensor( DOUBLETENSOR *m);


/*==============fucntions copied from utilities.h ===========*/


#define   FL     512

void stop_execution(void);
double Fmin(double a, double b);
long Fminlong(long a, long b);
double Fmax(double a, double b);
long Fmaxlong(long a, long b);

/*===============================functions copied from util_math.c========================================*/

//short tridiag2(short a, long r, long c, long nbeg, long nend, DOUBLEVECTOR *ld, DOUBLEVECTOR *d, DOUBLEVECTOR *ud, DOUBLEVECTOR *b, DOUBLEVECTOR *e);
short tridiag2(short a, long r, long c, long nbeg, long nend, const GeoVector<double>& ld, const GeoVector<double>& d, const GeoVector<double>& ud, const GeoVector<double>& b, GeoVector<double>& e)
;
//double norm_inf(DOUBLEVECTOR *V, long nbeg, long nend);
double norm_inf(const GeoVector<double>& V, long nbeg, long nend);
//double norm_2(DOUBLEVECTOR *V, long nbeg, long nend);
double norm_2(const GeoVector<double>& V, long nbeg, long nend);
void Cramer_rule(double A, double B, double C, double D, double E, double F, double *x, double *y);
double minimize_merit_function(double res0, double lambda1, double res1, double lambda2, double res2);
double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom);
double adaptiveSimpsons(double (*f)(double), double a, double b, double epsilon, int maxRecursionDepth);
double adaptiveSimpsonsAux2(double (*f)(double x, void *p), void *arg, double a, double b, double epsilon,
							double S, double fa, double fb, double fc, int bottom);
double adaptiveSimpsons2(double (*f)(double x, void *p), void *arg, double a, double b, double epsilon, int maxRecursionDepth);



#endif
