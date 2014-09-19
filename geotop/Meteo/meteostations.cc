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
