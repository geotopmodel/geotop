/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */


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

