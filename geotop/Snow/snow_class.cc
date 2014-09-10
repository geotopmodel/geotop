/*
 * @brief Snow Data implementation
 */

#include "snow_class.h"

Snow::Snow(double novalue, size_t total_pixel)
{
    yes = GeoVector<short>(total_pixel + 1, (short)0);
    t_snow = GeoVector<double>(total_pixel + 1, 0.);
    subl = GeoVector<double>(total_pixel + 1, 0.);
    melted = GeoVector<double>(total_pixel + 1, novalue);
    HNcum = GeoVector<double>(total_pixel + 1, 0.);
    age = GeoVector<double>(total_pixel + 1, novalue);
}

//FIXME: Horrible hack needed to cope with legacy code structure
void Snow::allocate_data(double novalue, size_t total_pixel)
{
    yes.resize(total_pixel + 1, (short)0);
    t_snow.resize(total_pixel + 1, 0.);
    subl.resize(total_pixel + 1, 0.);
    melted.resize(total_pixel + 1, novalue);
    HNcum.resize(total_pixel + 1, 0.);
    age.resize(total_pixel + 1, novalue);
}

