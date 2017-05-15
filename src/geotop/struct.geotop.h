
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
#ifndef STRUCT_GEOTOP_H
#define STRUCT_GEOTOP_H
#include "datastructs.h"

#include <vector>

/*---------------------------------------------------------------------------*/

class Energy
{
    public:
        GeoVector<double> Rn_mean;
        GeoVector<double> LWin_mean;
        GeoVector<double> LW_mean;
        GeoVector<double> SW_mean;
        GeoVector<double> ET_mean;
        GeoVector<double> H_mean;
        GeoVector<double> SEB_mean;
        GeoVector<double> Ts_mean;                /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
        GeoVector<double> Rswdown_mean;
        GeoVector<double> Rswbeam_mean;
        GeoVector<long> nDt_shadow;
        GeoVector<long> nDt_sun;
        GeoVector<double> Rn;
        GeoVector<double> LWin;
        GeoVector<double> LW;
        GeoVector<double> SW;
        GeoVector<double> LE;
        GeoVector<double> H;
        GeoVector<double> G;
        GeoVector<double> Ts;
        GeoVector<double> SWin;
        GeoVector<double> SWinb;
        GeoVector<short> shad;
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

        double *sun;
        double hsun;
        double sinhsun;
        double dsun;

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
};

/*---------------------------------------------------------------------------*/

class SoilState
{
    public:
        GeoMatrix<double> P;
        GeoMatrix<double> thi;
        GeoMatrix<double> T;
};

/*---------------------------------------------------------------------------*/

class StateVeg
{
    public:
        GeoVector<double> Tv;
        GeoVector<double> wrain;                  /*intercepted precipitation in mm*/
        GeoVector<double> wsnow;                  /*intercepted precipitation in mm*/
};

/*---------------------------------------------------------------------------*/

class Soil
{
    public:
        GeoMatrix<long> type;
        GeoVector<double> init_water_table_depth;
        GeoTensor<double> pa;
        GeoTensor<double> pa_bed;
        GeoMatrix<double> T_av_tensor;
        GeoMatrix<double> thw_av_tensor;
        GeoMatrix<double> thi_av_tensor;
        GeoMatrix<double> Ptot;
        GeoMatrix<double> th;
        GeoTensor<double> ET;
        GeoMatrix<double> Tzplot;
        GeoMatrix<double> Tzavplot;
        GeoMatrix<double> Ptotzplot;
        GeoMatrix<double> Pzplot;
        GeoMatrix<double> thzplot;
        GeoMatrix<double> thzavplot;
        GeoMatrix<double> thizplot;
        GeoMatrix<double> thizavplot;
        SoilState *SS;
        StateVeg *VS;
        GeoVector<double> Pnetcum;                //TODO mattiu
        GeoVector<double> ETcum;
};

/*---------------------------------------------------------------------------*/

class Topo
{
    public:
        GeoMatrix<double> Z0;                     //elevetions of each pixel (DEM)
        GeoTensor<double> Z;
        GeoMatrix<double> sky;                    //view factor (of the sky) for each pixel
        GeoMatrix<short> pixel_type;

        GeoMatrix<double> aspect;                 /*aspect; ex: matr_ev->azimuth*/
        GeoMatrix<double>  slope;                 /*slope of the pixels; ex: matr_ev->slope*/

        GeoMatrix<double> curvature1;
        GeoMatrix<double> curvature2;
        GeoMatrix<double> curvature3;
        GeoMatrix<double> curvature4;

        double ***horizon_height;
        long *horizon_numlines;
        GeoMatrix<long> horizon_point;

        long num_horizon_point;

        long ***i_cont;
        GeoMatrix<long> lrc_cont;

        long **j_cont;
        GeoMatrix<long> rc_cont;

        GeoVector<long> Lp;
        GeoVector<long> Li;

        GeoMatrix<long> Jdown;
        GeoMatrix<double> Qdown;

        GeoMatrix<short> is_on_border;

        GeoMatrix<double> East;
        GeoMatrix<double> North;

        GeoMatrix<long> BC_counter;

        GeoVector<double> BC_DepthFreeSurface;

        GeoMatrix<double> dzdE;
        GeoMatrix<double> dzdN;

        GeoMatrix<double> latitude;
        GeoMatrix<double> longitude;
};

/*---------------------------------------------------------------------------*/
class Land
{
    public:
        GeoMatrix<double> LC;                     //land cover (wood,lake,town,...) for each pixel*/
        GeoMatrix<double> delay;
        GeoMatrix<short> shadow;                  //=1 if shadow, =0 if not*/
        GeoMatrix<double> ty;

        double ***vegpars;
        double **vegparv;
        GeoVector<double> vegpar;

        long *NumlinesVegTimeDepData;

        GeoMatrix<double> root_fraction;
};                                                /*all this data are calculated on the basis of land use data and some other par*/

/*---------------------------------------------------------------------------*/

class Channel
{
    public:
        /*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                      R=number of rows of the basin,C=number of columns in the basin*/

        GeoVector<long> r;                        /*array of rows of the channel-pixels; dimension=nch*/
        GeoVector<long> c;                        /*array of columns of the channel-pixels; dimension=nch*/
        GeoMatrix<long> ch;
        GeoVector<long> ch_down;
        GeoVector<double> Vsup;
        GeoVector<double> Vsub;
        GeoVector<double> h_sup;
        GeoVector<double> length;

        double Vout;
        long **ch3;
        GeoMatrix<long> lch;
        GeoVector<long> soil_type;

        GeoMatrix<double> th;
        GeoMatrix<double> ET;

        GeoVector<double> Kbottom;

        SoilState *SS;
};

/*---------------------------------------------------------------------------*/

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
};

/*---------------------------------------------------------------------------*/
class Times
{
    public:
        GeoVector<double> JD_plots;
        double time;                              /*time=current time from the begin of simulation [s]*/
        long iplot;

};

/*---------------------------------------------------------------------------*/
class Par
{
    public:
        double Dt;                                /*Dt=the integration time interval [s]*/
        double ST;
        short print;                              /*1 IF YOU WANT TO PRINT MATRICES WITH INTERMEDIATE RESULTS, 0 OTHERWISE*/
        short monin_obukhov;
        double gamma_m;                           /*Exponent of the law of uniform motion on the surface*/
        double T_rain;                            /*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
        double T_snow;                            /*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
        double aep;                               /*ALBEDO EXTINCTION PARAMETER [m]*/
        double avo;                               /*NEW SNOW VISIBLE BAND REFLECTANCE*/
        double airo;                              /*NEW NEAR INFRARED BAND REFLECTANCE*/
        double Sr;                                /*WATER FRACTION RETAINED BY CAPILLARY FORCES IN SNOW*/
        double rho_ice;                           /*Ice density [kg/mc]*/
        long total_pixel;                         /*The number of the valid pixel of the whole basin*/
        long total_channel;
        double total_area;

        double max_weq_snow;
        long max_snow_layers;
        GeoVector<long> inf_snow_layers;

        double max_weq_glac;
        long max_glac_layers;
        GeoVector<long> inf_glac_layers;

        double Sr_glac;

        short state_turb;
        short state_lwrad;

        double imp;
        double f_bound_Richards;

        double epsilon_snow;

        GeoVector<double> output_soil;
        GeoVector<double> output_snow;
        GeoVector<double> output_glac;
        GeoVector<double> output_surfenergy;
        GeoVector<double> output_vegetation;
        GeoVector<double> output_meteo;

        short output_soil_bin;
        short output_snow_bin;
        short output_glac_bin;
        short output_surfenergy_bin;
        short output_meteo_bin;

        GeoMatrix<double> chkpt;
        GeoMatrix<long> rc;
        GeoVector<long> jplot;

        short recover;

        double Vmin;

        double snowcorrfact;
        double raincorrfact;

        double RHmin;

        short format_out;
        short sky;

        GeoVector<double> saving_points;
        double ContRecovery;
        long n_ContRecovery;

        short point_sim;

        double snow_maxpor;
        double snow_density_cutoff;
        double drysnowdef_rate;
        double wetsnowdef_rate;
        double snow_viscosity;

        double latitude;
        double longitude;

        double z0_snow;
        size_t n_landuses;

        short blowing_snow;

        GeoVector<long> r_points;
        GeoVector<long> c_points;

        double psimin;
        double stmin;

        short wat_balance;
        short en_balance;

        long nLC;

        double fetch_up;
        double fetch_down;

        short iobsint;
        double dn;
        double slopewt;
        double curvewt;
        double slopewtI;
        double curvewtI;
        double slopewtD;
        double curvewtD;

        short LRflag;

        GeoVector<short> vegflag;

        long MaxiterTol;
        double MaxErrWb;
        double TolVWb;

        double incr_V_factor;

        double alpha_snow;
        double tol_energy;
        long maxiter_energy;
        long maxiter_canopy;
        long maxiter_Ts;
        long maxiter_Loc;
        long maxiter_Businger;
        short stabcorr_incanopy;

        double TolCG;
        long MaxiterCorr;
        short UpdateK;

        short bedrock;

        double thres_hsup_1;
        double thres_hsup_2;

        double thres_hchannel;

        double w_dx;

        double RelTolVWb;

        double snow_smin;
        double snow_smax;
        double snow_curv;

        double Zboundary;
        double Tboundary;
        double Fboundary;

        double Ks_channel;
        double depr_channel;

        GeoVector<short> plot_discharge_with_Dt_integration;
        GeoVector<short> plot_point_with_Dt_integration;
        GeoVector<short> plot_basin_with_Dt_integration;

        GeoVector<double> Dtplot_point;
        GeoVector<double> Dtplot_basin;
        double Dtplot_discharge;

        short state_pixel;
        short state_discharge;
        short state_basin;

        double Dt_PBSM;

        long lowpass;
        long lowpass_curvatures;

        short dew;

        double init_date;
        double end_date;
        long run_times;

        double delay_day_recover;

        short all_point;
        short all_basin;
        short all_snow;
        short all_glac;
        short all_soil;

        double Wice_PBSM;

        GeoMatrix<double> maxSWE;

        long soil_type_land_default;
        long soil_type_chan_default;
        long soil_type_bedr_default;

        double MinIncrFactWithElev;
        double MaxIncrFactWithElev;

        long nsurface;

        double max_courant_land;
        double max_courant_channel;
        double max_courant_land_channel;
        double min_hsup_land;
        double min_hsup_channel;
        double min_dhsup_land_channel_in;
        double min_dhsup_land_channel_out;
        double dtmin_sup;

        long nsoiltypes;

        GeoVector<long> IDpoint;

        double min_lambda_en;
        long max_times_min_lambda_en;
        short exit_lambda_min_en;

        double min_lambda_wat;
        long max_times_min_lambda_wat;
        short exit_lambda_min_wat;

        double free_drainage_bottom;
        double free_drainage_lateral;

        short surroundings;

        GeoVector<double> soil_plot_depths;
        GeoVector<double> snow_plot_depths;
        GeoVector<double> glac_plot_depths;

        short ric_cloud;
        short vap_as_RH;
        short vap_as_Td;
        long ndivdaycloud;
        short cast_shadow;
        short wind_as_dir;
        short wind_as_xy;

        double snow_aging_vis;
        double snow_aging_nir;

        double DepthFreeSurface;

        short prec_as_intensity;

        GeoVector<short> linear_interpolation_meteo;

        short output_vertical_distances;

        short upwindblowingsnow;

        double Wmin_BS;
        double SWE_top;
        double SWE_bottom;
        double GWE_top;
        double GWE_bottom;

        double min_Dt;
        double dem_rotation;
        short usemeteoio;                         // flag indicating whether MeteoIO library is used
        bool use_meteoio_cloud;
        short use_meteoio_meteodata;
	bool use_ilwr_wrf;  // this flag enable ILWR from WRF plugin

        short qin;
        short flag1D;

        double k_to_ksat;
        short RunIfAnOldRunIsPresent;

        GeoVector<long> Nl_spinup ;

        short newperiodinit;

        short Tzrun;
        short wzrun;
        short dUzrun;
        short SWErun;

        short Tzmaxrun;
        short wzmaxrun;
        short Tzminrun;
        short wzminrun;

        double k1;
        double k2;
        double Lozone;
        double alpha_iqbal;
        double beta_iqbal;

        short albedoSWin;
        short micro;

        double EB;
        double Cair;
        double Tsup;

        double Tair_default;
        double RH_default;
        double V_default;
        double Vdir_default;
        double IPrec_default;

        double simulation_hours;

        double minP_torestore_A;
        short snow_conductivity;
        short snow_wind_compaction_1D;
        short snow_plot;

        short DDchannel;
        short DDland;

        double cum_prec;//cumulated precipitation since the beginning of an event [mm]
        double cum_da_up;//cumulated precipitation since the reset of the albedo [mm]
        int time_wo_prec;//time without precipitation after the previous event [sec]
        short evento;//states whether at the previous time stamp the event was ongoing
        short up_albedo;//states whether the albedo has already been reset for the ongoing event
        double tres_wo_prec;//max time interval allowed to differentiate between two contiguous events

        double VegVpdStess;
        double TvegMin;
        double TvegMax;
        double TvegRes;
        long VegRswStress;
        long VegVPDStress;
        long VegTempStress;
        long VegWaterStress;

};

class Statevar3D
{
    public:
        GeoMatrix<short> type;
        GeoMatrix<long> lnum;
        GeoTensor<double> Dzl;
        GeoTensor<double> w_liq;
        GeoTensor<double> w_ice;
        GeoTensor<double> T;
};

class Statevar1D
{
    public:
        short type;
        long lnum;
        GeoVector<double> Dzl;
        GeoVector<double>  w_liq;
        GeoVector<double> w_ice;
        GeoVector<double> T;
};

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

class Glacier
{
    public:
        Statevar3D *G;
        GeoVector<double> MELTED;
        GeoVector<double> melted;
        GeoVector<double> SUBL;
        GeoVector<double> subl;
};

class MeteoStations
{
    public:
        GeoVector<double> E;
        GeoVector<double> N;
        GeoVector<double> lat;
        GeoVector<double> lon;
        GeoVector<double> Z;

        GeoVector<double> sky;
        GeoVector<double> ST;
        GeoVector<double> Vheight;
        GeoVector<double> Theight;
        GeoVector<double> tau_cloud_av_meteoST;   // vector containing the tau_cloud_av at each meteo stations measuring SW radiation
        GeoVector<double> tau_cloud_meteoST;      // vector containing the tau_cloud at each meteo stations measuring SW radiation
        GeoVector<short> tau_cloud_av_yes_meteoST;// flag indicating whether the tau_cloud_av at each meteo stations is available
        GeoVector<short> tau_cloud_yes_meteoST;   // flag indicating whether the tau_cloud at each meteo stations is available
        GeoVector<short> flag_SW_meteoST;         // flag vector saying whether a meteo station accounts for SW radiation (0: no SW, 1: SW available)
};

class Meteo
{
    public:
        Meteo() : st(NULL), data(NULL), numlines(NULL), horizonlines(NULL),
            var(NULL), line_interp_WEB(NULL), line_interp_Bsnow(NULL), line_interp_WEB_LR(0), line_interp_Bsnow_LR(0){}

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

        GeoMatrix<double> tau_cloud;                         // tau_cloud used for shortwave
        GeoMatrix<double> tau_cloud_av;                      // tau_cloud used for longwave (averaged in a wider interval)
        GeoMatrix<short> tau_cloud_yes;                      // it is read (1) or used default (0)
        GeoMatrix<short> tau_cloud_av_yes;
 
    
//  MB 3.1.2017
//        GeoMatrix<double> tau_cl_map;             // matrix containing the tau_cloud for each grid point
//        GeoMatrix<double> tau_cl_av_map;          // matrix containing the tau_cloud_av for each grid point
//        GeoMatrix<short> tau_cl_map_yes;          // boolean matrix saying whether the grid point has tau_cl value
//        GeoMatrix<short> tau_cl_av_map_yes;       // boolean matrix saying whether the grid point has tau_cl_av value

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
};


class AllData
{
    public:
        Soil *S;
        Water *W;
        Land *L;
        Par *P;
        Topo *T;
        Channel *C;
        Energy *E;
        Snow *N;
        Glacier *G;
        Meteo *M;
        Times *I;

};
#endif
