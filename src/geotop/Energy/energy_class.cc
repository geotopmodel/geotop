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

#include "energy_class.h"

Energy::Energy(double /*novalue*/, size_t total_pixel)
{
  Rn_mean = GeoVector<double>(total_pixel + 1, 0.);
  LWin_mean = GeoVector<double>(total_pixel + 1, 0.);
  LW_mean = GeoVector<double>(total_pixel + 1, 0.);
  SW_mean = GeoVector<double>(total_pixel + 1, 0.);
  ET_mean = GeoVector<double>(total_pixel + 1, 0.);
  H_mean = GeoVector<double>(total_pixel + 1, 0.);
  SEB_mean = GeoVector<double>(total_pixel + 1, 0.);
  Ts_mean = GeoVector<double>(total_pixel + 1, 0.);
  Rswdown_mean = GeoVector<double>(total_pixel + 1, 0.);
  Rswbeam_mean = GeoVector<double>(total_pixel + 1, 0.);

  nDt_shadow = GeoVector<long>(total_pixel + 1, 0);
  nDt_sun = GeoVector<long>(total_pixel + 1, 0);
  Rn = GeoVector<double>(total_pixel + 1);
  LWin = GeoVector<double>(total_pixel + 1);
  LW = GeoVector<double>(total_pixel + 1);
  SW = GeoVector<double>(total_pixel + 1);
  LE = GeoVector<double>(total_pixel + 1);
  H = GeoVector<double>(total_pixel + 1);
  G = GeoVector<double>(total_pixel + 1);
  Ts = GeoVector<double>(total_pixel + 1);
  SWin = GeoVector<double>(total_pixel + 1);
  SWinb = GeoVector<double>(total_pixel + 1);
  shad = GeoVector<short>(total_pixel + 1, 0);
}

// FIXME: Horrible hack needed to cope with legacy code structure
void Energy::allocate_data(double /*novalue*/, size_t total_pixel)
{
  Rn_mean.resize(total_pixel + 1, 0.);
  LWin_mean.resize(total_pixel + 1, 0.);
  LW_mean.resize(total_pixel + 1, 0.);
  SW_mean.resize(total_pixel + 1, 0.);
  ET_mean.resize(total_pixel + 1, 0.);
  H_mean.resize(total_pixel + 1, 0.);
  SEB_mean.resize(total_pixel + 1, 0.);
  Ts_mean.resize(total_pixel + 1, 0.);
  Rswdown_mean.resize(total_pixel + 1, 0.);
  Rswbeam_mean.resize(total_pixel + 1, 0.);

  nDt_shadow.resize(total_pixel + 1, 0.);
  nDt_sun.resize(total_pixel + 1, 0.);
  Rn.resize(total_pixel + 1);
  LWin.resize(total_pixel + 1);
  LW.resize(total_pixel + 1);
  SW.resize(total_pixel + 1);
  LE.resize(total_pixel + 1);
  H.resize(total_pixel + 1);
  G.resize(total_pixel + 1);
  Ts.resize(total_pixel + 1);
  SWin.resize(total_pixel + 1);
  SWinb.resize(total_pixel + 1);
  shad.resize(total_pixel + 1, 0.);
}
