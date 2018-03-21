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

/**
 * @brief Snow and Glacier State Variables implementation
 * @date September 2014
 */

#include "statevar.h"

Statevar3D::Statevar3D(double novalue,
                       size_t layers,
                       size_t rows,
                       size_t columns)
{
  type = GeoMatrix<short>(rows + 1, columns + 1, 2);
  lnum = GeoMatrix<long>(rows + 1, columns + 1, 0);

  // FIXME: GeoTensor's constructor assumes that the order of indices is:
  // Number of Rows, Number of Columns and Number of Layers BUT the rest of the
  // code  assumes that the order is layers, rows and then columns.
  Dzl = GeoTensor<double>(layers + 1, rows + 1, columns + 1, 0.);
  w_liq = GeoTensor<double>(layers + 1, rows + 1, columns + 1, 0.);
  w_ice = GeoTensor<double>(layers + 1, rows + 1, columns + 1, 0.);
  T = GeoTensor<double>(layers + 1, rows + 1, columns + 1, novalue);
}

Statevar1D::Statevar1D(double novalue, size_t layers)
{
  Dzl = GeoVector<double>(layers + 1, 0.0);
  T = GeoVector<double>(layers + 1, novalue);
  w_ice = GeoVector<double>(layers + 1, 0.0);
  w_liq = GeoVector<double>(layers + 1, 0.0);
}
