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

/**
 * @brief Meteo Stations Variables
 * @date September 2014
 */

#ifndef METEOSTATIONS_H
#define METEOSTATIONS_H

#include "../datastructs.h"

class MeteoStations
{
public:
  GeoVector<double> E;
  GeoVector<double> N;
  GeoVector<double> lat;
  GeoVector<double> lon;
  GeoVector<double> Z;
  GeoVector<double> sky;
  GeoVector<double> ST;
  GeoVector<double> Vheight;
  GeoVector<double> Theight;
  GeoVector<double> tau_cloud_av_meteoST;  // vector containing the tau_cloud_av
  // at each meteo stations measuring
  // SW radiation
  GeoVector<double> tau_cloud_meteoST;     // vector containing the tau_cloud at
  // each meteo stations measuring SW
  // radiation
  GeoVector<short> tau_cloud_av_yes_meteoST;  // flag indicating whether the
  // tau_cloud_av at each meteo
  // stations is available
  GeoVector<short> tau_cloud_yes_meteoST;     // flag indicating whether the
  // tau_cloud at each meteo stations
  // is available
  GeoVector<short> flag_SW_meteoST;  // flag vector saying whether a meteo
  // station accounts for SW radiation (0: no
  // SW, 1: SW available)

  MeteoStations() {};
  MeteoStations(size_t nmeteo_stations, double novalue);
};

#endif
