
/* STATEMENT:
 overloaded function
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

#ifndef OUTPUT_H
#define OUTPUT_H

#include "struct.geotop.h"
#include "pedo.funct.h"
#include "../libraries/ascii/rw_maps.h"
#include "constants.h"
#include "times.h"
#include "energy.balance.h"
#include "input.h"
#include "../libraries/ascii/tabs.h"
#include "vegetation.h"
#include "tables.h"
#include "snow.h"

#include <time.h>

void write_output(Times *times,
                  Water *wat,
                  Channel *cnet,
                  Par *par,
                  Topo *top,
                  Land *land,
                  Soil *sl,
                  Energy *egy,
                  Snow *snow,
                  Glacier *glac,
                  Meteo *met);

void write_output_headers(long n,
                          Times *times,
                          Water *wat,
                          Par *par,
                          Topo *top,
                          Land *land,
                          Soil *sl,
                          Energy *egy,
                          Snow *snow,
                          Glacier *glac);

void write_soil_output(long i,
                       long iname,
                       double init_date,
                       double JDfrom0,
                       double JD,
                       long day,
                       long month,
                       long year,
                       long hour,
                       long minute,
                       const GeoVector<double> &n,
                       Soil *sl,
                       Par *par,
                       double psimin,
                       double cosslope);

void write_soil_file(long lmin,
                     long i,
                     FILE *f,
                     long d,
                     long m,
                     long y,
                     long h,
                     long mi,
                     double JDfrom0,
                     double JDfrom0init,
                     double *var,
                     const GeoVector<double> &n,
                     double *dz,
                     double cosslope);
void write_soil_file(long lmin,
                     long i,
                     FILE *f,
                     long d,
                     long m,
                     long y,
                     long h,
                     long mi,
                     double JDfrom0,
                     double JDfrom0init,
                     const GeoMatrix<double> &var,
                     long row,
                     const GeoVector<double> &n,
                     const GeoTensor<double> &dz,
                     double cosslope);

void write_soil_header(FILE *f,
                       const GeoVector<double> &n,
                       const GeoTensor<double> &dz);

void plot(std::string name,
          long i_plot,
          const GeoVector<double> &V,
          short format,
          long **J);

double interpolate_soil(long lmin,
                        double h,
                        long max,
                        const GeoTensor<double> &Dz,
                        const GeoMatrix<double> &Q,
                        const long &row);

double interpolate_soil2(long lmin,
                         double h,
                         long max,
                         const GeoTensor<double> &Dz,
                         GeoMatrix<double> &Q,
                         long i);

void write_tensorseries_soil(long lmin,
                             std::string suf,
                             std::string filename,
                             short type,
                             short format,
                             GeoMatrix<double> &T,
                             const GeoVector<double> &n,
                             long **J,
                             GeoMatrix<long> &RC,
                             GeoTensor<double> &dz,
                             GeoMatrix<double> &slope,
                             short vertical);
void fill_output_vectors(double Dt,
                         double W,
                         Energy *egy,
                         Snow *snow,
                         Glacier *glac,
                         Water *wat,
                         Meteo *met,
                         Par *par,
                         Times *time,
                         Topo *top,
                         Soil *sl);

void write_snow_file(short choice,
                     long IDpoint,
                     long r,
                     long c,
                     long actual_snow_layer_numb,
                     long max_snow_layer_numb,
                     FILE *f,
                     long d,
                     long m,
                     long y,
                     long h,
                     long mi,
                     double JDfrom0,
                     double JDfrom0init,
                     const GeoVector<double> &plot_depth,
                     const GeoTensor<double> &dz,
                     const GeoTensor<double> &var_to_print,
                     double cosslope);
void write_snow_output(long i,
                       long iname,
                       long r,
                       long c,
                       double init_date,
                       double JDfrom0,
                       long day,
                       long month,
                       long year,
                       long hour,
                       long minute,
                       const GeoVector<double> &plot_depth,
                       Statevar3D *Snow,
                       Par *par,
                       double cosslope);

void write_snow_header(FILE *f,
                       const GeoVector<double> &plot_depth,
                       size_t max_snow_layer);

#endif
