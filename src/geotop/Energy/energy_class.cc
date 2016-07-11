/*
 * @brief Energy Data implementation
 */

#include "energy_class.h"

Energy::Energy(double novalue, size_t total_pixel)
{
    Rn_mean = GeoVector<double>(total_pixel + 1, 0.);
    LWin_mean = GeoVector<double>(total_pixel + 1, 0.);
    LW_mean = GeoVector<double>(total_pixel + 1, 0.);
    SW_mean = GeoVector<double>(total_pixel + 1, 0.);
    ET_mean = GeoVector<double>(total_pixel + 1, 0.);
    H_mean = GeoVector<double>(total_pixel + 1, 0.);
    SEB_mean = GeoVector<double>(total_pixel + 1, 0.);
    Ts_mean = GeoVector<double>(total_pixel + 1, 0.);
    Rswdown_mean = GeoVector<double>(total_pixel + 1, 0.);
    Rswbeam_mean = GeoVector<double>(total_pixel + 1, 0.);


    nDt_shadow = GeoVector<long>(total_pixel + 1, 0.);
    nDt_sun = GeoVector<long>(total_pixel + 1, 0.);
    Rn = GeoVector<double>(total_pixel + 1);
    LWin = GeoVector<double>(total_pixel + 1);
    LW = GeoVector<double>(total_pixel + 1);
    SW = GeoVector<double>(total_pixel + 1);
    LE = GeoVector<double>(total_pixel + 1);
    H = GeoVector<double>(total_pixel + 1);
    G = GeoVector<double>(total_pixel + 1);
    Ts = GeoVector<double>(total_pixel + 1);
    SWin = GeoVector<double>(total_pixel + 1);
    SWinb = GeoVector<double>(total_pixel + 1);
    shad = GeoVector<short>(total_pixel + 1, 0.);


    GeoVector<double> Hgplot;
    GeoVector<double> LEgplot;
    GeoVector<double> Hvplot;
    GeoVector<double> LEvplot;
    GeoVector<double> SWinplot;
    GeoVector<double> SWgplot;
    GeoVector<double> SWvplot;
    GeoVector<double> LWinplot;
    GeoVector<double> LWgplot;
    GeoVector<double> LWvplot;
    GeoVector<double> Tgplot;
    GeoVector<double> Tsplot;
    GeoVector<double> Tvplot;
    GeoVector<double> Hgp;
    GeoVector<double> LEgp;
    GeoVector<double> Hvp;
    GeoVector<double> LEvp;
    GeoVector<double> SWinp;
    GeoVector<double> SWgp;
    GeoVector<double> SWvp;
    GeoVector<double> LWinp;
    GeoVector<double> LWgp;
    GeoVector<double> LWvp;
    GeoVector<double> Tgp;
    GeoVector<double> Tsp;

	//those variables are not used inside of this code
    //double *sun;
    //double hsun;
	//double sinhsun;
    //double dsun;

    GeoVector<double> Dlayer;
    GeoVector<double> liq;
    GeoVector<double> ice;
    GeoVector<double> Temp;
    GeoVector<double> deltaw;
    GeoVector<double> SWlayer;
    GeoVector<double> soil_transp_layer;
    GeoVector<double> dFenergy;
    GeoVector<double> udFenergy;
    GeoVector<double> Kth0;
    GeoVector<double> Kth1;
    GeoVector<double> Fenergy;
    GeoVector<double> Newton_dir;
    GeoVector<double> T0;
    GeoVector<double> T1;
    GeoVector<double> Tstar;
    GeoVector<double> THETA;
    GeoVector<double> soil_evap_layer_bare;
    GeoVector<double> soil_evap_layer_veg;
    GeoMatrix<double> Tgskin_surr;
    GeoMatrix<double> SWrefl_surr;
}

//FIXME: Horrible hack needed to cope with legacy code structure
void Energy::allocate_data(double novalue, size_t total_pixel)
{
    Rn_mean.resize(total_pixel + 1, 0.);
    LWin_mean.resize(total_pixel + 1, 0.);
    LW_mean.resize(total_pixel + 1, 0.);
    SW_mean.resize(total_pixel + 1, 0.);
    ET_mean.resize(total_pixel + 1, 0.);
    H_mean.resize(total_pixel + 1, 0.);
    SEB_mean.resize(total_pixel + 1, 0.);
    Ts_mean.resize(total_pixel + 1, 0.);
    Rswdown_mean.resize(total_pixel + 1, 0.);
    Rswbeam_mean.resize(total_pixel + 1, 0.);
    

    nDt_shadow.resize(total_pixel + 1, 0.);
    nDt_sun.resize(total_pixel + 1, 0.);
    Rn.resize(total_pixel + 1);
    LWin.resize(total_pixel + 1);
    LW.resize(total_pixel + 1);
    SW.resize(total_pixel + 1);
    LE.resize(total_pixel + 1);
    H.resize(total_pixel + 1);
    G.resize(total_pixel + 1);
    Ts.resize(total_pixel + 1);
    SWin.resize(total_pixel + 1);
    SWinb.resize(total_pixel + 1);
    shad.resize(total_pixel + 1, 0.);


    GeoVector<double> Hgplot;
    GeoVector<double> LEgplot;
    GeoVector<double> Hvplot;
    GeoVector<double> LEvplot;
    GeoVector<double> SWinplot;
    GeoVector<double> SWgplot;
    GeoVector<double> SWvplot;
    GeoVector<double> LWinplot;
    GeoVector<double> LWgplot;
    GeoVector<double> LWvplot;
    GeoVector<double> Tgplot;
    GeoVector<double> Tsplot;
    GeoVector<double> Tvplot;
    GeoVector<double> Hgp;
    GeoVector<double> LEgp;
    GeoVector<double> Hvp;
    GeoVector<double> LEvp;
    GeoVector<double> SWinp;
    GeoVector<double> SWgp;
    GeoVector<double> SWvp;
    GeoVector<double> LWinp;
    GeoVector<double> LWgp;
    GeoVector<double> LWvp;
    GeoVector<double> Tgp;
    GeoVector<double> Tsp;

	//those variables are not used inside of this code
    //double *sun;
    //double hsun;
    //double sinhsun;
    //double dsun;

    GeoVector<double> Dlayer;
    GeoVector<double> liq;
    GeoVector<double> ice;
    GeoVector<double> Temp;
    GeoVector<double> deltaw;
    GeoVector<double> SWlayer;
    GeoVector<double> soil_transp_layer;
    GeoVector<double> dFenergy;
    GeoVector<double> udFenergy;
    GeoVector<double> Kth0;
    GeoVector<double> Kth1;
    GeoVector<double> Fenergy;
    GeoVector<double> Newton_dir;
    GeoVector<double> T0;
    GeoVector<double> T1;
    GeoVector<double> Tstar;
    GeoVector<double> THETA;
    GeoVector<double> soil_evap_layer_bare;
    GeoVector<double> soil_evap_layer_veg;
    GeoMatrix<double> Tgskin_surr;
    GeoMatrix<double> SWrefl_surr;
}

