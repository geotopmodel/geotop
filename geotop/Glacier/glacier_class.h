/*
 * @brief Glacier Data definition
 */

#ifndef GLACIER_CLASS_H
#define GLACIER_CLASS_H

#include "../statevar.h"

class Glacier
{
    public:
        Statevar3D *G;
        GeoVector<double> MELTED;
        GeoVector<double> melted;
        GeoVector<double> SUBL;
        GeoVector<double> subl;
        Glacier() {};
        Glacier(double novalue, size_t total_pixel);
        void allocate_data(double novalue, size_t total_pixel);
};

#endif

