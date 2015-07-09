/**
 * @brief Soil State Variables
 * @date September 2014
 */

#ifndef SOILSTATEVAR_H
#define SOILSTATEVAR_H

#include "../datastructs.h"

class SoilState
{
public:
    GeoMatrix<double> P;
    GeoMatrix<double> thi;
    GeoMatrix<double> T;

    SoilState() {};
    SoilState(size_t total_pixel, size_t layers);
};


#endif
