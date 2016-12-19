#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "turtle.h"
#include "../geotop/constants.h"
#include <vector>


/*==============functions copied from utilities.h ===========*/


#define   FL     512

void stop_execution(void);
double Fmin(double a, double b);
long Fminlong(long a, long b);
double Fmax(double a, double b);
long Fmaxlong(long a, long b);

short tridiag2(long nbeg, long nend, const GeoVector<double>& ld, const GeoVector<double>& d, const GeoVector<double>& ud, const GeoVector<double>& b, GeoVector<double>& e);

double norm_inf(const GeoVector<double>& V, long nbeg, long nend);
double norm_1(const GeoVector<double>& V, long nbeg, long nend);
double norm_2(const GeoVector<double>& V, long nbeg, long nend);
void Cramer_rule(double A, double B, double C, double D, double E, double F, double *x, double *y);
double minimize_merit_function(double res0, double lambda1, double res1, double lambda2, double res2);
double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom);
double adaptiveSimpsons(double (*f)(double), double a, double b, double epsilon, int maxRecursionDepth);
double adaptiveSimpsonsAux2(double (*f)(double x, void *p), void *arg, double a, double b, double epsilon,
							double S, double fa, double fb, double fc, int bottom);
double adaptiveSimpsons2(double (*f)(double x, void *p), void *arg, double a, double b, double epsilon, int maxRecursionDepth);

void product_matrix_using_lower_part_by_vector_plus_vector(double k, GeoVector<double>& out, const GeoVector<double>& y, const GeoVector<double>& x, const GeoVector<long>& Li, const GeoVector<long>& Lp, GeoVector<double>& Lx);

void product_using_only_strict_lower_diagonal_part(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<long>& Li, const GeoVector<long>& Lp, GeoVector<double>& Lx);

long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(double tol_rel, double tol_min, double tol_max, GeoVector<double>& x,
                                                          const GeoVector<double>& b, GeoVector<double>& y, const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);

short tridiag(short a, long r, long c, long nx, const GeoVector<double>& diag_inf, const GeoVector<double>& diag, const GeoVector<double>& diag_sup, const GeoVector<double>& b, GeoVector<double>& e);

void get_diag_strict_lower_matrix_plus_identity_by_vector(GeoVector<double>& diag, GeoVector<double>& udiag, const GeoVector<double>& y,const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);
void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<double>& y,const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);

double product(const GeoVector<double>& a, const GeoVector<double>& b);




#endif
