/**
 * @brief Vegetation State Variables implementation
 * @date September 2014
 */

#include "vegstatevar.h"

StateVeg::StateVeg(size_t total_pixel)
{
    Tv = GeoVector<double>(total_pixel, 0.);
    wrain = GeoVector<double>(total_pixel, 0.);                  /*intercepted precipitation in mm*/
    wsnow = GeoVector<double>(total_pixel, 0.);                  /*intercepted precipitation in mm*/
}
