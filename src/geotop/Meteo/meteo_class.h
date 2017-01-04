/*
 * @brief Meteo Data definition
 */

#ifndef METEO_CLASS_H
#define METEO_CLASS_H

#include "meteostations.h"

class Meteo
{
public:
    // Meteo() : st(NULL), data(NULL), numlines(NULL), horizonlines(NULL),
    //           var(NULL), line_interp_WEB(NULL), line_interp_Bsnow(NULL),
    //           line_interp_WEB_LR(0), line_interp_Bsnow_LR(0), tau_cloud(0.),
    //           tau_cloud_av(0.), tau_cloud_yes(0), tau_cloud_av_yes(0) {}

    MeteoStations *st;

    double ***data;
    long *numlines;
    double ***horizon;
    long *horizonlines;
    double **var;
    long *line_interp_WEB;
    long *line_interp_Bsnow;
    long line_interp_WEB_LR;
    long line_interp_Bsnow_LR;
    double **LRs;                             //matrix read from the external value
    long LRsnr;                               //number of lines of the matrix
    double *LRv;                              //vector of interpolated values
    double **LRc;                             //cyclic values from the parameter file (one vector for each LR variable)
    long *LRcnc;                              //number of components of the vector (for each component)
    double *LRd;                              //vector of default values

    double **qins;
    double *qinv;
    long qinsnr;
    long qinline;


//    MB 3.1.2017
//    double tau_cloud;                         // tau_cloud for the chosen meteo station used to derive cloud
//    double tau_cloud_av;                      // tau_cloud for the chosen meteo station used to derive cloud
//    short tau_cloud_yes;
//    short tau_cloud_av_yes;
    GeoMatrix<double> tau_cloud;                         // tau_cloud used for shortwave
    GeoMatrix<double> tau_cloud_av;                      // tau_cloud used for longwave (averaged in a wider interval)
    GeoMatrix<short> tau_cloud_yes;                      // it is read (1) or used default (0)
    GeoMatrix<short> tau_cloud_av_yes;
    
    GeoMatrix<double> Tgrid;
    GeoMatrix<double> Pgrid;
    GeoMatrix<double> Vgrid;
    GeoMatrix<double> Vdir;
    GeoMatrix<double> RHgrid;
    GeoMatrix<double> ILWRgrid;

    GeoVector<double> Tamean;
    GeoVector<double> Vspdmean;
    GeoVector<double> Vdirmean;
    GeoVector<double> RHmean;

    GeoVector<double> Taplot;
    GeoVector<double> Vxplot;
    GeoVector<double> Vyplot;
    GeoVector<double> RHplot;

    double V;

    long nstcloud;                            // meteo station ID (1...n) to use for the cloudiness
    long numstcloud;                          // number of meteo stations measuring cloudiness
    long nstsrad;
    long nstlrad;
    long nstTs;

    GeoVector<long> imeteo_stations;
    Meteo() {};
    Meteo(double novalue, size_t nrows, size_t ncols);
    Meteo(size_t nrows, size_t ncols, double Vmin, double Pa);
    void allocate_data(double novalue, size_t nrows, size_t ncols, size_t Z0nrows, size_t Z0ncols, double Vmin, double Pa);
};

#endif

