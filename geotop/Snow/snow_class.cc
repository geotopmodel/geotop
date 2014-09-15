/*
 * @brief Snow Data implementation
 */

#include "snow_class.h"

Snow::Snow(double novalue, size_t total_pixel)
{
    age = GeoVector<double>(total_pixel + 1, novalue);
    melted = GeoVector<double>(total_pixel + 1, novalue);
    HNcum = GeoVector<double>(total_pixel + 1, 0.);
    subl = GeoVector<double>(total_pixel + 1, 0.);
    t_snow = GeoVector<double>(total_pixel + 1, 0.);
    yes = GeoVector<short>(total_pixel + 1, (short)0);
    Wsubl_plot = GeoMatrix<double>(total_pixel + 1, 0.);
    Wtrans_plot = GeoMatrix<double>(total_pixel + 1, 0.);
    Dplot = GeoVector<double>(total_pixel + 1, 0.);
}

//FIXME: Horrible hack needed to cope with legacy code structure
void Snow::allocate_data(double novalue, size_t total_pixel)
{
    age.resize(total_pixel + 1, novalue);
    melted.resize(total_pixel + 1, novalue);
    HNcum.resize(total_pixel + 1, 0.);
    subl.resize(total_pixel + 1, 0.);
    t_snow.resize(total_pixel + 1, 0.);
    yes.resize(total_pixel + 1, (short)0);
    Wsubl_plot.resize(total_pixel + 1, 0.);
    Wtrans_plot.resize(total_pixel + 1, 0.);
    Dplot.resize(total_pixel + 1, 0.);
}

