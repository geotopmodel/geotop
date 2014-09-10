/**
 * @brief Snow and Glacier State Variables
 * @date September 2014
 */

#ifndef STATEVAR_H
#define STATEVAR_H

#include "datastructs.h"

class Statevar3D
{
public:
    GeoMatrix<short> type;
    GeoMatrix<long> lnum;
    GeoTensor<double> Dzl;
    GeoTensor<double> w_liq;
    GeoTensor<double> w_ice;
    GeoTensor<double> T;
    Statevar3D() {};
    Statevar3D(double novalue, size_t layers, size_t rows, size_t columns);
};

class Statevar1D
{
public:
    short type;
    long lnum; //TODO: check this and see if it can be changed to private
    GeoVector<double> Dzl;
    GeoVector<double> w_liq;
    GeoVector<double> w_ice;
    GeoVector<double> T;
    Statevar1D() {};
    Statevar1D(double novalue, size_t layers);
};


#endif


