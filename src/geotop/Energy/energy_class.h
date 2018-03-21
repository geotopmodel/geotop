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
 * @brief Energy Data definition
 */

#ifndef ENERGY_CLASS_H
#define ENERGY_CLASS_H

#include "../statevar.h"

class Energy
{
public:
  // TODO Variable_mean have to be deleted. Now they are declared only to
  // compile correctly
  GeoVector<double> Rn_mean;
  GeoVector<double> LWin_mean;
  GeoVector<double> LW_mean;
  GeoVector<double> SW_mean;
  GeoVector<double> ET_mean;
  GeoVector<double> H_mean;
  GeoVector<double> SEB_mean;
  GeoVector<double>
  Ts_mean; /*averaged surface Temperature(on nDt_output_basin Dt time
                intervals)*/
  GeoVector<double> Rswdown_mean;
  GeoVector<double> Rswbeam_mean;

  GeoVector<long> nDt_shadow;
  GeoVector<long> nDt_sun;
  GeoVector<double> Rn;
  GeoVector<double> LWin;
  GeoVector<double> LW;
  GeoVector<double> SW;
  GeoVector<double> LE;
  GeoVector<double> H;
  GeoVector<double> G;
  GeoVector<double> Ts;
  GeoVector<double> SWin;
  GeoVector<double> SWinb;
  GeoVector<short> shad;

  GeoVector<double> Hgplot;
  GeoVector<double> LEgplot;
  GeoVector<double> Hvplot;
  GeoVector<double> LEvplot;
  GeoVector<double> SWinplot;
  GeoVector<double> SWgplot;
  GeoVector<double> SWvplot;
  GeoVector<double> LWinplot;
  GeoVector<double> LWgplot;
  GeoVector<double> LWvplot;
  GeoVector<double> Tgplot;
  GeoVector<double> Tsplot;
  GeoVector<double> Tvplot;
  GeoVector<double> Hgp;
  GeoVector<double> LEgp;
  GeoVector<double> Hvp;
  GeoVector<double> LEvp;
  GeoVector<double> SWinp;
  GeoVector<double> SWgp;
  GeoVector<double> SWvp;
  GeoVector<double> LWinp;
  GeoVector<double> LWgp;
  GeoVector<double> LWvp;
  GeoVector<double> Tgp;
  GeoVector<double> Tsp;

  double *sun;
  double hsun;
  double sinhsun;
  double dsun;

  GeoVector<double> Dlayer;
  GeoVector<double> liq;
  GeoVector<double> ice;
  GeoVector<double> Temp;
  GeoVector<double> deltaw;
  GeoVector<double> SWlayer;
  GeoVector<double> soil_transp_layer;
  GeoVector<double> dFenergy;
  GeoVector<double> udFenergy;
  GeoVector<double> Kth0;
  GeoVector<double> Kth1;
  GeoVector<double> Fenergy;
  GeoVector<double> Newton_dir;
  GeoVector<double> T0;
  GeoVector<double> T1;
  GeoVector<double> Tstar;
  GeoVector<double> THETA;
  GeoVector<double> soil_evap_layer_bare;
  GeoVector<double> soil_evap_layer_veg;
  GeoMatrix<double> Tgskin_surr;
  GeoMatrix<double> SWrefl_surr;
  Energy() {};
  Energy(double novalue, size_t total_pixel);
  void allocate_data(double novalue, size_t total_pixel);
};

#endif
