/*

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 - 31 December 2016

 Copyright (c), 2016 - GEOtop Foundation

 This file is part of GEOtop 2.1

 GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

 GEOtop 2.1 is distributed as a free software in the hope to create and support
 a community of developers and users that constructively interact. If you just
 use the code, please give feedback to GEOtop Foundation and the community. Any
 way you use the model, may be the most trivial one, is significantly helpful
 for the future development of the GEOtop model. Any feedback will be highly
 appreciated.

 */

/*
 * @brief Energy Data implementation
 */

#include "meteo_class.h"

// MB 3.1.17
/*Meteo::Meteo(double novalue, size_t Z0nrows, size_t Z0ncols)
{
   tau_cl_map = GeoMatrix<double>(Z0nrows, Z0ncols, novalue);
   tau_cl_av_map = GeoMatrix<double>(Z0nrows, Z0ncols, novalue);
   tau_cl_map_yes = GeoMatrix<short>(Z0nrows, Z0ncols, (short)novalue);
   tau_cl_av_map_yes = GeoMatrix<short>(Z0nrows, Z0ncols, (short)novalue);
}*/

Meteo::Meteo(size_t nrows, size_t ncols, double Vmin, double Pa)
{
  Tgrid = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, 5.);
  Pgrid = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, Pa);
  Vgrid = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, Vmin);
  Vdir = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, 0.);
  RHgrid = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, 0.7);
  ILWRgrid = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, 0);
  tau_cloud = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, 1.0);
  tau_cloud_av = GeoMatrix<double>((int)nrows + 1, (int)ncols + 1, 1.0);
  tau_cloud_yes = GeoMatrix<short>((int)nrows + 1, (int)ncols + 1, 0);
  tau_cloud_av_yes = GeoMatrix<short>((int)nrows + 1, (int)ncols + 1, 0);
}

// FIXME: Horrible hack needed to cope with legacy code structure
void Meteo::allocate_data(double novalue,
                          size_t nrows,
                          size_t ncols,
                          size_t /*Z0nrows*/,
                          size_t /*Z0ncols*/,
                          double Vmin,
                          double Pa)
{
  //   tau_cl_map.resize(Z0nrows, Z0ncols, novalue);
  //   tau_cl_av_map.resize(Z0nrows, Z0ncols, novalue);
  //   tau_cl_map_yes.resize(Z0nrows, Z0ncols, (short)novalue);
  //   tau_cl_av_map_yes.resize(Z0nrows, Z0ncols, (short)novalue);

  Tgrid.resize(nrows + 1, ncols + 1, 5.);
  Pgrid.resize(nrows + 1, ncols + 1, Pa);
  Vgrid.resize(nrows + 1, ncols + 1, Vmin);
  Vdir.resize(nrows + 1, ncols + 1, 0.);
  RHgrid.resize(nrows + 1, ncols + 1, 0.7);
  ILWRgrid.resize(nrows + 1, ncols + 1, novalue);
  tau_cloud.resize(nrows + 1, ncols + 1, 1.0);
  tau_cloud_av.resize(nrows + 1, ncols + 1, 1.0);
  tau_cloud_yes.resize(nrows + 1, ncols + 1, 0);
  tau_cloud_av_yes.resize(nrows + 1, ncols + 1, 0);
}
