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

/*
 * @brief Soil Data definition
 */

#ifndef SOIL_CLASS_H
#define SOIL_CLASS_H

#include "soilstatevar.h"
#include "vegstatevar.h"

class Soil
{
public:
  SoilState *SS;
  StateVeg *VS;

  GeoMatrix<double> T_av_tensor;
  GeoMatrix<double> thw_av_tensor;
  GeoMatrix<double> thi_av_tensor;
  GeoMatrix<double> Ptot;
  GeoMatrix<double> th;
  GeoTensor<double> ET;

  // Computational Variables
  GeoMatrix<long> type;
  GeoVector<double> init_water_table_depth;
  GeoTensor<double> pa;
  GeoTensor<double> pa_bed;

  // Special plot e cumulated not used
  GeoMatrix<double> Tzplot;
  GeoMatrix<double> Tzavplot;
  GeoMatrix<double> Ptotzplot;
  GeoMatrix<double> Pzplot;
  GeoMatrix<double> thzplot;
  GeoMatrix<double> thzavplot;
  GeoMatrix<double> thizplot;
  GeoMatrix<double> thizavplot;
  GeoVector<double> Pnetcum;  // TODO mattiu
  GeoVector<double> ETcum;

  Soil() {}

  Soil(double novalue,
       size_t layers,
       size_t nrows,
       size_t ncols,
       size_t total_pixel);
  void allocate_data(double novalue,
                     size_t layers,
                     size_t nrows,
                     size_t ncols,
                     size_t total_pixel);
};

#endif
