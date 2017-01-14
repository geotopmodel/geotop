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
 * @brief Water Data definition
 */

#ifndef WATER_CLASS_H
#define WATER_CLASS_H

#include "../datastructs.h"

class Water                                       /*nstations=number of all the rain-stations,number_of_pixels=number of all the pixels of the basin R*C,*/
{
    public:                                       /* R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
        GeoMatrix<double> PrecTot;                /*total(snow+rain) precipitation in mm (in a Dt)*/
        GeoMatrix<double> Pnet;                  /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                                                   of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                                                   same subroutine and in "water.balance.c" module*/
        GeoMatrix<double> HN;                     // map of new snow TODO mattiu
        GeoVector<double> PrTOT_mean;             /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
        GeoVector<double> PrSNW_mean;
        GeoVector<double> Pt;
        GeoVector<double> Ps;
        GeoVector<double> h_sup;


        GeoMatrix<double> error;
        GeoVector<double> Lx;
        GeoVector<double> H0;
        GeoVector<double> H1;
        GeoVector<double> dH;
        GeoVector<double> B;
        GeoVector<double> f;
        GeoVector<double> df;
        GeoMatrix<double> Klat;
        GeoMatrix<double> Kbottom;
        double Voutlandsub;
        double Voutlandsup;
        double Voutbottom;
        Water() {};
        Water(double novalue, size_t nrows, size_t ncols, size_t total_pixel);
        void allocate_data(double novalue, size_t nrows, size_t ncols, size_t total_pixel);
		
};

#endif

