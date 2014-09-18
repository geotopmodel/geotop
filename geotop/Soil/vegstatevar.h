/**
 * @brief Vegetation State Variables
 * @date September 2014
 */

#ifndef VEGSTATEVAR_H
#define VEGSTATEVAR_H

#include "../datastructs.h"

class StateVeg
{
public:
    GeoVector<double> Tv;
    GeoVector<double> wrain;                  /*intercepted precipitation in mm*/
    GeoVector<double> wsnow;                  /*intercepted precipitation in mm*/

    StateVeg() {};
    StateVeg(size_t total_pixel);
};


#endif
