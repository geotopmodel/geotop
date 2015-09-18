/*
 * @brief Soil Data implementation
 */

#include "soil_class.h"

Soil::Soil(double novalue, size_t layers, size_t nrows, size_t ncols, size_t total_pixel)
{
    T_av_tensor = GeoMatrix<double>(layers + 1, total_pixel + 1, 0.);
    thw_av_tensor = GeoMatrix<double>(layers + 1, total_pixel + 1, 0.);
    thi_av_tensor = GeoMatrix<double>(layers + 1, total_pixel + 1, 0.);
    Ptot = GeoMatrix<double>(layers + 1, total_pixel + 1, novalue);
    th = GeoMatrix<double>(layers + 1, total_pixel + 1, novalue);
    ET = GeoTensor<double>(layers + 1, nrows + 1, ncols + 1, 0.);
}

//FIXME: Horrible hack needed to cope with legacy code structure
void Soil::allocate_data(double novalue, size_t layers, size_t nrows, size_t ncols, size_t total_pixel)
{
    T_av_tensor.resize(layers + 1, total_pixel + 1, 0.);
    thw_av_tensor.resize(layers + 1, total_pixel + 1, 0.);
    thi_av_tensor.resize(layers + 1, total_pixel + 1, 0.);
    Ptot.resize(layers + 1, total_pixel + 1, novalue);
    th.resize(layers + 1, total_pixel + 1, novalue);
    ET.resize(layers + 1, nrows + 1, ncols + 1, 0.);
}

