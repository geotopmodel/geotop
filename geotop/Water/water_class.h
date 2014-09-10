/*
 * @brief Water Data definition
 */

#ifndef WATER_CLASS_H
#define WATER_CLASS_H

// #include "../statevar.h"
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

