/*
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 - 31 December 2016
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1 is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to GEOtop Foundation and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model.
 Any feedback will be highly appreciated.
 
 */

/**
 * @brief Meteo Stations Variables implementation
 * @date September 2014
 */

#include "meteostations.h"

MeteoStations::MeteoStations(size_t nmeteo_stations, double novalue)
{
    E = GeoVector<double>(nmeteo_stations);
    N = GeoVector<double>(nmeteo_stations);
    lat = GeoVector<double>(nmeteo_stations);
    lon = GeoVector<double>(nmeteo_stations);
    Z = GeoVector<double>(nmeteo_stations);
    sky = GeoVector<double>(nmeteo_stations);
    ST = GeoVector<double>(nmeteo_stations);
    Vheight = GeoVector<double>(nmeteo_stations);
    Theight = GeoVector<double>(nmeteo_stations);

    tau_cloud_av_meteoST = GeoVector<double>(nmeteo_stations, novalue);
    tau_cloud_meteoST = GeoVector<double>(nmeteo_stations, novalue);
    tau_cloud_av_yes_meteoST = GeoVector<short>(nmeteo_stations, novalue);
    tau_cloud_yes_meteoST = GeoVector<short>(nmeteo_stations, novalue);
    flag_SW_meteoST = GeoVector<short>(nmeteo_stations, novalue);
}
