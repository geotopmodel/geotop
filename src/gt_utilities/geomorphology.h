/* STATEMENT:

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)

 Copyright (c), 2016 - GEOtop Foundation

 This file is part of GEOtop 2.1

 GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

 GEOtop 2.1  is distributed as a free software in the hope to create and support
 a community of developers and users that constructively interact. If you just
 use the code, please give feedback to the authors and the community. Any way
 you use the model, may be the most trivial one, is significantly helpful for
 the future development of the GEOtop model. Any feedback will be highly
 appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#ifndef T_IO_H
#define T_IO_H
#include "turtle.h"
#include "math_utils.h"
#include <vector>

void multipass_topofilter(long ntimes,
                          GeoMatrix<double> &Zin,
                          GeoMatrix<double> &Zout,
                          long novalue,
                          long n);
void find_min_max(GeoMatrix<double> &M, long novalue, double *max,
                  double *min);

/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of part in which you want divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti       */

void sky_view_factor(GeoMatrix<double> &sky,
                     size_t N,
                     TInit *UV,
                     GeoMatrix<double> &input,
                     GeoMatrix<short> &convess,
                     long novalue);
void find_aspect(GeoMatrix<double> &topo,
                 GeoMatrix<double> &dzdx,
                 GeoMatrix<double> &dzdy,
                 long undef,
                 GeoMatrix<double> &M);

void topofilter(GeoMatrix<double> &Zin,
                GeoMatrix<double> &Zout,
                long novalue,
                long n);
void find_slope(double deltax,
                double deltay,
                GeoMatrix<double> &topo,
                GeoMatrix<double> &dzdx,
                GeoMatrix<double> &dzdy,
                long undef);
void find_max_slope(GeoMatrix<double> &topo,
                    GeoMatrix<double> &dzdx,
                    GeoMatrix<double> &dzdy,
                    long undef,
                    GeoMatrix<double> &M);

#endif
