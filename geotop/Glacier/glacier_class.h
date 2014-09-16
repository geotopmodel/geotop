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

        GeoVector<double> melted;
        GeoVector<double> subl;

        GeoVector<double> MELTED; // this is the cumulative value of the melted instant variable, unuseful with the new output
        GeoVector<double> SUBL; // this is the cumulative value of the subl instant variable, unuseful with the new output
        Glacier() {};
        Glacier(double novalue, size_t total_pixel);
        void allocate_data(double novalue, size_t total_pixel);
};

#endif

