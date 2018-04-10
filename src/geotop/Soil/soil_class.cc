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
 * @brief Soil Data implementation
 */

#include "soil_class.h"

Soil::Soil(double novalue,
           size_t layers,
           size_t nrows,
           size_t ncols,
           size_t total_pixel):
  T_av_tensor{layers + 1, total_pixel + 1, 0.},
  thw_av_tensor{layers + 1, total_pixel + 1, 0.},
  thi_av_tensor{layers + 1, total_pixel + 1, 0.},
  Ptot{layers + 1, total_pixel + 1, novalue},
  th{layers + 1, total_pixel + 1, novalue},
  ET{layers + 1, nrows + 1, ncols + 1, 0.}
{}

// FIXME: Horrible hack needed to cope with legacy code structure
void Soil::allocate_data(double novalue,
                         size_t layers,
                         size_t nrows,
                         size_t ncols,
                         size_t total_pixel)
{
  T_av_tensor.resize(layers + 1, total_pixel + 1, 0.);
  thw_av_tensor.resize(layers + 1, total_pixel + 1, 0.);
  thi_av_tensor.resize(layers + 1, total_pixel + 1, 0.);
  Ptot.resize(layers + 1, total_pixel + 1, novalue);
  th.resize(layers + 1, total_pixel + 1, novalue);
  ET.resize(layers + 1, nrows + 1, ncols + 1, 0.);
}
