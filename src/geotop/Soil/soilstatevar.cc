/**
 * @brief Soil State Variables implementation
 * @date September 2014
 */

#include "soilstatevar.h"

SoilState::SoilState(size_t total_pixel, size_t layers)
{
    P = GeoMatrix<double>(layers + 1, total_pixel + 1, 0.);
    thi = GeoMatrix<double>(layers + 1, total_pixel + 1, 0.);
    T = GeoMatrix<double>(layers + 1, total_pixel + 1);
}
