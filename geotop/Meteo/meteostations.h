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
    GeoVector<double> tau_cloud_av_meteoST;   // vector containing the tau_cloud_av at each meteo stations measuring SW radiation
    GeoVector<double> tau_cloud_meteoST;      // vector containing the tau_cloud at each meteo stations measuring SW radiation
    GeoVector<short> tau_cloud_av_yes_meteoST;// flag indicating whether the tau_cloud_av at each meteo stations is available
    GeoVector<short> tau_cloud_yes_meteoST;   // flag indicating whether the tau_cloud at each meteo stations is available
    GeoVector<short> flag_SW_meteoST;         // flag vector saying whether a meteo station accounts for SW radiation (0: no SW, 1: SW available)

    MeteoStations() {};
    MeteoStations(size_t nmeteo_stations, double novalue);
};


#endif
