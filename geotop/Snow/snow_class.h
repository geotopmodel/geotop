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
        GeoVector<double> MELTED;
        GeoVector<double> melted;
        GeoVector<double> HNcum;                  // TODO mattiu
        GeoVector<double> SUBL;
        GeoVector<double> subl;
        GeoVector<double> t_snow;
        GeoVector<short> yes;
        GeoMatrix<double> Qsub;
        GeoMatrix<double> Qsub_x;
        GeoMatrix<double> Qsub_y;
        GeoMatrix<double> Nabla2_Qtrans;
        GeoMatrix<double> Qtrans;
        GeoMatrix<double> Qsalt;
        GeoMatrix<double> Qtrans_x;
        GeoMatrix<double> Qtrans_y;
        GeoMatrix<double> Wsubl_plot;
        GeoMatrix<double> Wtrans_plot;
        GeoVector<double> Dplot;
        GeoVector<long> change_dir_wind;
};

#endif

