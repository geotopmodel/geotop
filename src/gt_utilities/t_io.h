#ifndef T_IO_H
#define T_IO_H
#include "turtle.h"
#include "math_utils.h"
#include <vector>


void multipass_topofilter(long ntimes, GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n);
void find_min_max(GeoMatrix<double>& M, long novalue, double *max, double *min);

double norm_1(const GeoVector<double>& V, long nbeg, long nend);

/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of part in which you want divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti       */

void sky_view_factor(GeoMatrix<double>& sky, long N, TInit *UV, GeoMatrix<double>& input, GeoMatrix<short>& convess, long novalue);
void find_aspect(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M);

void topofilter(GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n);
void find_slope(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef);
void find_max_slope(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M);


#endif
