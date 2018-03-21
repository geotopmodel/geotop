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
 * @brief Water Data implementation
 */

#include "water_class.h"

Water::Water(double novalue, size_t nrows, size_t ncols, size_t total_pixel)
{
  // output variables
  PrecTot = GeoMatrix<double>(nrows + 1, ncols + 1, 0.);
  Pnet = GeoMatrix<double>(nrows + 1, ncols + 1, 0.);
  HN = GeoMatrix<double>(nrows + 1, ncols + 1, 0.);
  PrTOT_mean = GeoVector<double>(total_pixel + 1, 0.);
  PrSNW_mean = GeoVector<double>(total_pixel + 1, 0.);
  Pt = GeoVector<double>(total_pixel + 1, novalue);
  Ps = GeoVector<double>(total_pixel + 1, novalue);
  // computational variables
  h_sup = GeoVector<double>(total_pixel + 1, 0.);
  // error = GeoMatrix<double>(nlayers, total_pixel + 1, novalue);
  // Lx = GeoVector<double>(total_pixel + 1, novalue);
  // H0 = GeoVector<double>(total_pixel + 1, novalue);
  // H1 = GeoVector<double>(total_pixel + 1, novalue);
  // dH = GeoVector<double>(total_pixel + 1, novalue);
  // B = GeoVector<double>(total_pixel + 1, novalue);
  // f = GeoVector<double>(total_pixel + 1, novalue);
  // df = GeoVector<double>(total_pixel + 1, novalue);
  // Klat = GeoMatrix<double>(nlayers, total_pixel + 1, 0.);
  // Kbottom = GeoMatrix<double>(nlayers, total_pixel + 1, 0.);
}

// FIXME: Horrible hack needed to cope with legacy code structure
void Water::allocate_data(double novalue,
                          size_t nrows,
                          size_t ncols,
                          size_t total_pixel)
{
  // output variables
  PrecTot.resize(nrows + 1, ncols + 1, 0.);
  Pnet.resize(nrows + 1, ncols + 1, 0.);
  HN.resize(nrows + 1, ncols + 1, 0.);
  PrTOT_mean.resize(total_pixel + 1, 0.);
  PrSNW_mean.resize(total_pixel + 1, 0.);
  Pt.resize(total_pixel + 1, novalue);
  Ps.resize(total_pixel + 1, novalue);
  // computational variables
  h_sup.resize(total_pixel + 1, 0.);
  // error.resize(nlayers, total_pixel + 1, novalue);
  // Lx.resize(total_pixel + 1, novalue);
  // H0.resize(total_pixel + 1, novalue);
  // H1.resize(total_pixel + 1, novalue);
  // dH.resize(total_pixel + 1, novalue);
  // B.resize(total_pixel + 1, novalue);
  // f.resize(total_pixel + 1, novalue);
  // df.resize(total_pixel + 1, novalue);

  // Klat.resize(nlayers, total_pixel + 1, 0.);
  // Kbottom.resize(nlayers, total_pixel + 1, 0.);

  Voutlandsub = 0.;
  Voutlandsup = 0.;
  Voutbottom = 0.;
}
