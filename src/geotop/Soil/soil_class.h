/*
 * @brief Soil Data definition
 */

#ifndef SOIL_CLASS_H
#define SOIL_CLASS_H

#include "soilstatevar.h"
#include "vegstatevar.h"

class Soil
{
    public:
        SoilState *SS;
        StateVeg *VS;

        GeoMatrix<double> T_av_tensor;
        GeoMatrix<double> thw_av_tensor;
        GeoMatrix<double> thi_av_tensor;
        GeoMatrix<double> Ptot;
        GeoMatrix<double> th;
        GeoTensor<double> ET;

        // Computational Variables
        GeoMatrix<long> type;
        GeoVector<double> init_water_table_depth;
        GeoTensor<double> pa;
        GeoTensor<double> pa_bed;

        // Special plot e cumulated not used
        GeoMatrix<double> Tzplot;
        GeoMatrix<double> Tzavplot;
        GeoMatrix<double> Ptotzplot;
        GeoMatrix<double> Pzplot;
        GeoMatrix<double> thzplot;
        GeoMatrix<double> thzavplot;
        GeoMatrix<double> thizplot;
        GeoMatrix<double> thizavplot;
        GeoVector<double> Pnetcum;                //TODO mattiu
        GeoVector<double> ETcum;

        Soil() {};
        Soil(double novalue, size_t layers,  size_t nrows, size_t ncols, size_t total_pixel);
        void allocate_data(double novalue, size_t layers,  size_t nrows, size_t ncols, size_t total_pixel);
};

#endif

