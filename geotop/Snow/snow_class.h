/*
 * @brief Snow Data definition
 */

#ifndef SNOW_CLASS_H
#define SNOW_CLASS_H

#include "../statevar.h"

class Snow
{
    public:
        Statevar3D *S;
        Statevar1D *S_for_BS;

        GeoVector<double> age;
        GeoVector<double> melted;
        GeoVector<double> HNcum;                  // TODO mattiu
        GeoVector<double> subl;
        GeoVector<double> t_snow;
        GeoVector<short> yes;
        GeoMatrix<double> Wsubl_plot;
        GeoMatrix<double> Wtrans_plot;
        GeoVector<double> Dplot;

        // computational variables
        GeoMatrix<double> Qsub;
        GeoMatrix<double> Qsub_x;
        GeoMatrix<double> Qsub_y;
        GeoMatrix<double> Nabla2_Qtrans;
        GeoMatrix<double> Qtrans;
        GeoMatrix<double> Qsalt;
        GeoMatrix<double> Qtrans_x;
        GeoMatrix<double> Qtrans_y;
        GeoVector<long> change_dir_wind;


        GeoVector<double> MELTED; // this is the cumulative value of the melted instant variable, unuseful with the new output
        GeoVector<double> SUBL; // this is the cumulative value of the subl instant variable, unuseful with the new output

        Snow() {};
        Snow(double novalue, size_t total_pixel);
        void allocate_data(double novalue, size_t total_pixel);
};

#endif

