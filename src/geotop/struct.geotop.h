
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#include "turtle.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "tensor3D.h"
#include <memory>

/*---------------------------------------------------------------------------*/
typedef struct
{

  std::unique_ptr<Vector<double>> Rn_mean;
  std::unique_ptr<Vector<double>> LWin_mean;
  std::unique_ptr<Vector<double>> LW_mean;
  std::unique_ptr<Vector<double>> SW_mean;
  std::unique_ptr<Vector<double>> ET_mean;
  std::unique_ptr<Vector<double>> H_mean;
  std::unique_ptr<Vector<double>> SEB_mean;
  DOUBLEVECTOR
  *Ts_mean;  /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
  std::unique_ptr<Vector<double>> Rswdown_mean;
  std::unique_ptr<Vector<double>> Rswbeam_mean;
  LONGVECTOR *nDt_shadow;
  LONGVECTOR *nDt_sun;

  std::unique_ptr<Vector<double>> Rn;
  std::unique_ptr<Vector<double>> LWin;
  std::unique_ptr<Vector<double>> LW;
  std::unique_ptr<Vector<double>> SW;
  std::unique_ptr<Vector<double>> LE;
  std::unique_ptr<Vector<double>> H;
  std::unique_ptr<Vector<double>> G;
  std::unique_ptr<Vector<double>> Ts;
  std::unique_ptr<Vector<double>> SWin;
  std::unique_ptr<Vector<double>> SWinb;
  SHORTVECTOR *shad;

  std::unique_ptr<Vector<double>> Hgplot;
  std::unique_ptr<Vector<double>> LEgplot;
  std::unique_ptr<Vector<double>> Hvplot;
  std::unique_ptr<Vector<double>> LEvplot;
  std::unique_ptr<Vector<double>> SWinplot;
  std::unique_ptr<Vector<double>> SWgplot;
  std::unique_ptr<Vector<double>> SWvplot;
  std::unique_ptr<Vector<double>> LWinplot;
  std::unique_ptr<Vector<double>> LWgplot;
  std::unique_ptr<Vector<double>> LWvplot;
  std::unique_ptr<Vector<double>> Tgplot;
  std::unique_ptr<Vector<double>> Tsplot;
  std::unique_ptr<Vector<double>> Tvplot;

  std::unique_ptr<Vector<double>> Hgp;
  std::unique_ptr<Vector<double>> LEgp;
  std::unique_ptr<Vector<double>> Hvp;
  std::unique_ptr<Vector<double>> LEvp;
  std::unique_ptr<Vector<double>> SWinp;
  std::unique_ptr<Vector<double>> SWgp;
  std::unique_ptr<Vector<double>> SWvp;
  std::unique_ptr<Vector<double>> LWinp;
  std::unique_ptr<Vector<double>> LWgp;
  std::unique_ptr<Vector<double>> LWvp;
  std::unique_ptr<Vector<double>> Tgp;
  std::unique_ptr<Vector<double>> Tsp;

  double *sun;
  double hsun;
  double sinhsun;
  double dsun;

  std::unique_ptr<Vector<double>> Dlayer;
  std::unique_ptr<Vector<double>> liq;
  std::unique_ptr<Vector<double>> ice;
  std::unique_ptr<Vector<double>> Temp;
  std::unique_ptr<Vector<double>> deltaw;
  std::unique_ptr<Vector<double>> SWlayer;
  std::unique_ptr<Vector<double>> soil_transp_layer;
  std::unique_ptr<Vector<double>> dFenergy;
  std::unique_ptr<Vector<double>> udFenergy;
  std::unique_ptr<Vector<double>> Kth0;
  std::unique_ptr<Vector<double>> Kth1;
  std::unique_ptr<Vector<double>> Fenergy;
  std::unique_ptr<Vector<double>> Newton_dir;
  std::unique_ptr<Vector<double>> T0;
  std::unique_ptr<Vector<double>> T1;
  std::unique_ptr<Vector<double>> Tstar;
  std::unique_ptr<Vector<double>> THETA;
  std::unique_ptr<Vector<double>> soil_evap_layer_bare;
  std::unique_ptr<Vector<double>> soil_evap_layer_veg;

  DOUBLEMATRIX *Tgskin_surr;
  DOUBLEMATRIX *SWrefl_surr;

} ENERGY;

/*---------------------------------------------------------------------------*/

typedef struct
{
  DOUBLEMATRIX *P;
  DOUBLEMATRIX *thi;
  DOUBLEMATRIX *T;
} SOIL_STATE;

/*---------------------------------------------------------------------------*/

typedef struct
{

  std::unique_ptr<Vector<double>> Tv;
  std::unique_ptr<Vector<double>> wrain;       /*intercepted precipitation in mm*/
  std::unique_ptr<Vector<double>> wsnow;       /*intercepted precipitation in mm*/

} STATE_VEG;

/*---------------------------------------------------------------------------*/

typedef struct
{
  LONGMATRIX *type;
  DOUBLETENSOR *pa;
  DOUBLEMATRIX *T_av_tensor;
  DOUBLEMATRIX *thw_av_tensor;
  DOUBLEMATRIX *thi_av_tensor;
  DOUBLEMATRIX *Ptot;
  DOUBLEMATRIX *th;
  DOUBLETENSOR *ET;
  DOUBLEMATRIX *Tzplot;
  DOUBLEMATRIX *Tzavplot;
  DOUBLEMATRIX *Ptotzplot;
  DOUBLEMATRIX *Pzplot;
  DOUBLEMATRIX *thzplot;
  DOUBLEMATRIX *thzavplot;
  DOUBLEMATRIX *thizplot;
  DOUBLEMATRIX *thizavplot;
  DOUBLEMATRIX *satratio;
  SOIL_STATE *SS;
  STATE_VEG *VS;

  DOUBLEMATRIX *Tzrun;
  DOUBLEMATRIX *wzrun;
  DOUBLEMATRIX *dUzrun;
  DOUBLEMATRIX *SWErun;

  DOUBLEMATRIX *Tzmaxrun;
  DOUBLEMATRIX *wzmaxrun;
  DOUBLEMATRIX *Tzminrun;
  DOUBLEMATRIX *wzminrun;

  std::unique_ptr<Vector<double>> Pnetcum;
  std::unique_ptr<Vector<double>> ETcum;

} SOIL;



/*---------------------------------------------------------------------------*/
typedef struct
{
  DOUBLEMATRIX *Z0;         //elevation of each pixel (DEM)
  DOUBLETENSOR *Z;

  DOUBLEMATRIX *sky;        //view factor (of the sky) for each pixel
  SHORTMATRIX *pixel_type;

  //SHORTMATRIX *DD;      //Drainage Directions for each pixel; ex matr_ev->slope*/
  //LONGMATRIX *DDup;
  //LONGVECTOR *DDdown;
  //DOUBLEMATRIX *i_DD;       /*slope along Drainage Direction for each pixel*/

  DOUBLEMATRIX *aspect;     /*aspect; ex: matr_ev->azimuth*/
  DOUBLEMATRIX *slope;     /*slope of the pixels; ex: matr_ev->slope*/

  DOUBLEMATRIX *curvature1;
  DOUBLEMATRIX *curvature2;
  DOUBLEMATRIX *curvature3;
  DOUBLEMATRIX *curvature4;

  double ***horizon_height;
  long *horizon_numlines;
  LONGMATRIX *horizon_point;
  long num_horizon_point;

  long ***i_cont;
  LONGMATRIX *lrc_cont;

  long **j_cont;
  LONGMATRIX *rc_cont;

  LONGVECTOR *Lp;
  LONGVECTOR *Li;
  //LONGVECTOR *Up;
  //LONGVECTOR *Ui;

  LONGMATRIX *Jdown;
  DOUBLEMATRIX *Qdown;

  SHORTMATRIX *is_on_border;

  DOUBLEMATRIX *East;
  DOUBLEMATRIX *North;

  LONGMATRIX *BC_counter;
  std::unique_ptr<Vector<double>> BC_DepthFreeSurface;

  DOUBLEMATRIX *dzdE;
  DOUBLEMATRIX *dzdN;

  DOUBLEMATRIX *latitude;
  DOUBLEMATRIX *longitude;

} TOPO;


/*---------------------------------------------------------------------------*/
typedef struct
{
  DOUBLEMATRIX *LC;            //land cover for each pixel
  DOUBLEMATRIX *delay;
  SHORTMATRIX *shadow;      //=1 if shadow, =0 if not
  DOUBLEMATRIX *ty;

  double ***vegpars;
  double **vegparv;
  std::unique_ptr<Vector<double>> vegpar;
  long *NumlinesVegTimeDepData;

  DOUBLEMATRIX *root_fraction;

} LAND;/*all this data are calculated on the basis of land use data and some other par*/


/*---------------------------------------------------------------------------*/


typedef struct
{
  /*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                    R=number of rows of the basin,C=number of columns in the basin*/
  LONGVECTOR *r;          /*array of rows of the channel-pixels; dimension=nch*/
  LONGVECTOR
  *c;          /*array of columns of the channel-pixels; dimension=nch*/
  LONGMATRIX *ch;
  LONGVECTOR *ch_down;
  std::unique_ptr<Vector<double>> Vsup;
  std::unique_ptr<Vector<double>> Vsub;
  std::unique_ptr<Vector<double>> h_sup;
  std::unique_ptr<Vector<double>> length;
  double Vout;
  long **ch3;
  LONGMATRIX *lch;
  LONGVECTOR *soil_type;
  DOUBLEMATRIX *th;
  DOUBLEMATRIX *ET;
  std::unique_ptr<Vector<double>> Kbottom;
  SOIL_STATE *SS;
} CHANNEL;




/*---------------------------------------------------------------------------*/

typedef struct
{
  /*nstations=number of all the rain-stations,number_of_pixels=number of all the pixels of the basin R*C,
                     R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
  DOUBLEMATRIX *PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/
  DOUBLEMATRIX
  *Pnet;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                              same subroutine and in "water.balance.c" module*/

  DOUBLEVECTOR
  *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
  std::unique_ptr<Vector<double>> PrSNW_mean;
  std::unique_ptr<Vector<double>> Pt;
  std::unique_ptr<Vector<double>> Ps;

  std::unique_ptr<Vector<double>> h_sup;

  DOUBLEMATRIX *error;

  std::unique_ptr<Vector<double>> Lx;
  std::unique_ptr<Vector<double>> Ux;

  std::unique_ptr<Vector<double>> H0;
  std::unique_ptr<Vector<double>> H1;
  std::unique_ptr<Vector<double>> dH;
  std::unique_ptr<Vector<double>> B;
  std::unique_ptr<Vector<double>> f;
  std::unique_ptr<Vector<double>> df;
  DOUBLEMATRIX *Klat;
  DOUBLEMATRIX *Kbottom;

  double Voutlandsub;
  double Voutlandsup;
  double Voutbottom;

} WATER;


/*---------------------------------------------------------------------------*/
typedef struct
{

  std::unique_ptr<Vector<double>> JD_plots;
  double time;    /*time=current time from the begin of simulation [s]*/
  long iplot;
  double **Dt_matrix;
  long numlinesDt_matrix;
  double *Dt_vector;

} TIMES;


/*---------------------------------------------------------------------------*/
typedef struct
{
  double Dt;      /*Dt=the integration time interval [s]*/
  double ST;
  short print;         /*1 IF YOU WANT TO PRINT MATRICES WITH INTERMEDIATE RESULTS, 0 OTHERWISE*/
  short monin_obukhov;
  double gamma_m;   /*Exponent of the law of uniform motion on the surface*/
  double T_rain;    /*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
  double T_snow;    /*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
  double aep;       /*ALBEDO EXTINCTION PARAMETER [m]*/
  double avo;       /*NEW SNOW VISIBLE BAND REFLECTANCE*/
  double airo;      /*NEW NEAR INFRARED BAND REFLECTANCE*/
  double Sr;      /*WATER FRACTION RETAINED BY CAPILLARY FORCES IN SNOW*/
  double rho_ice;     /*Ice density [kg/mc]*/
  long total_pixel;    /*The number of the valid pixel of the whole basin*/
  long total_channel;
  double total_area;

  double max_weq_snow;
  long max_snow_layers;
  LONGVECTOR *inf_snow_layers;

  double max_weq_glac;
  long max_glac_layers;
  LONGVECTOR *inf_glac_layers;

  double Sr_glac;

  short state_turb;
  short state_lwrad;

  double imp;
  double f_bound_Richards;

  double epsilon_snow;

  std::unique_ptr<Vector<double>> output_soil;
  std::unique_ptr<Vector<double>> output_snow;
  std::unique_ptr<Vector<double>> output_glac;
  std::unique_ptr<Vector<double>> output_surfenergy;
  std::unique_ptr<Vector<double>> output_vegetation;
  std::unique_ptr<Vector<double>> output_meteo;

  short output_soil_bin;
  short output_snow_bin;
  short output_glac_bin;
  short output_surfenergy_bin;
  short output_meteo_bin;

  DOUBLEMATRIX *chkpt;
  LONGMATRIX *rc;
  LONGVECTOR *jplot;

  short recover;

  double Vmin;

  double snowcorrfact;
  double raincorrfact;

  double RHmin;

  short format_out;
  short sky;

  std::unique_ptr<Vector<double>> saving_points;
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
  long n_landuses;

  short blowing_snow;

  LONGVECTOR *r_points;
  LONGVECTOR *c_points;

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

  SHORTVECTOR *vegflag;

  //short harm_or_arit_mean_normal;
  //short harm_or_arit_mean_parallel;

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

  short tsteps_from_file;

  SHORTVECTOR *plot_discharge_with_Dt_integration;
  SHORTVECTOR *plot_point_with_Dt_integration;
  SHORTVECTOR *plot_basin_with_Dt_integration;

  std::unique_ptr<Vector<double>> Dtplot_point;
  std::unique_ptr<Vector<double>> Dtplot_basin;
  std::unique_ptr<Vector<double>> Dtplot_discharge;

  short state_pixel;
  short state_discharge;
  short state_basin;

  double Dt_PBSM;

  long lowpass;
  long lowpass_curvatures;

  short dew;

  std::unique_ptr<Vector<double>> init_date;
  std::unique_ptr<Vector<double>> end_date;
  LONGVECTOR *run_times;

  double delay_day_recover;

  short all_point;
  short all_basin;
  short all_snow;
  short all_glac;
  short all_soil;

  double Wice_PBSM;

  DOUBLEMATRIX *maxSWE;

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

  LONGVECTOR *IDpoint;

  double min_lambda_en;
  long max_times_min_lambda_en;
  short exit_lambda_min_en;

  double min_lambda_wat;
  long max_times_min_lambda_wat;
  short exit_lambda_min_wat;

  double free_drainage_bottom;
  double free_drainage_lateral;

  short surroundings;

  std::unique_ptr<Vector<double>> soil_plot_depths;
  std::unique_ptr<Vector<double>> snow_plot_depths;
  std::unique_ptr<Vector<double>> glac_plot_depths;

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

  SHORTVECTOR *linear_interpolation_meteo;

  short output_vertical_distances;

  short upwindblowingsnow;

  double Wmin_BS;
  double SWE_top;
  double SWE_bottom;
  double GWE_top;
  double GWE_bottom;

  double min_Dt;
  double dem_rotation;

  short qin;
  short flag1D;

  double k_to_ksat;
  short RunIfAnOldRunIsPresent;

  LONGVECTOR *Nl_spinup;

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
  double Tbottom;

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

} PAR;



/*---------------------------------------------------------------------------*/
typedef struct
{
  SHORTMATRIX  *type;
  LONGMATRIX   *lnum;
  DOUBLETENSOR *Dzl;
  DOUBLETENSOR *w_liq;
  DOUBLETENSOR *w_ice;
  DOUBLETENSOR *T;
} STATEVAR_3D;


typedef struct
{
  short type;
  long lnum;
  std::unique_ptr<Vector<double>> Dzl;
  std::unique_ptr<Vector<double>> w_liq;
  std::unique_ptr<Vector<double>> w_ice;
  std::unique_ptr<Vector<double>> T;
} STATEVAR_1D;

typedef struct
{
  STATEVAR_3D *S;
  STATEVAR_1D *S_for_BS;
  std::unique_ptr<Vector<double>> age;
  std::unique_ptr<Vector<double>> MELTED;
  std::unique_ptr<Vector<double>> melted;
  std::unique_ptr<Vector<double>> SUBL;
  std::unique_ptr<Vector<double>> subl;
  std::unique_ptr<Vector<double>> t_snow;
  SHORTVECTOR *yes;
  DOUBLEMATRIX *Qsub;
  DOUBLEMATRIX *Qsub_x;
  DOUBLEMATRIX *Qsub_y;
  DOUBLEMATRIX *Nabla2_Qtrans;
  DOUBLEMATRIX *Qtrans;
  DOUBLEMATRIX *Qsalt;
  DOUBLEMATRIX *Qtrans_x;
  DOUBLEMATRIX *Qtrans_y;
  DOUBLEMATRIX *Wsubl_plot;
  DOUBLEMATRIX *Wtrans_plot;
  std::unique_ptr<Vector<double>> Dplot;
  LONGVECTOR *change_dir_wind;
} SNOW;

typedef struct
{
  STATEVAR_3D *G;
  std::unique_ptr<Vector<double>> MELTED;
  std::unique_ptr<Vector<double>> melted;
  std::unique_ptr<Vector<double>> SUBL;
  std::unique_ptr<Vector<double>> subl;
} GLACIER;

typedef struct
{
  std::unique_ptr<Vector<double>> E;
  std::unique_ptr<Vector<double>> N;
  std::unique_ptr<Vector<double>> lat;
  std::unique_ptr<Vector<double>> lon;
  std::unique_ptr<Vector<double>> Z;
  std::unique_ptr<Vector<double>> sky;
  std::unique_ptr<Vector<double>> ST;
  std::unique_ptr<Vector<double>> Vheight;
  std::unique_ptr<Vector<double>> Theight;
} METEO_STATIONS;


typedef struct
{
  METEO_STATIONS *st;

  double ***data;
  long *numlines;
  double ***horizon;
  long *horizonlines;
  double **var;
  long *line_interp_WEB;
  long *line_interp_Bsnow;
  long line_interp_WEB_LR;
  long line_interp_Bsnow_LR;

  double **LRs; //matrix read from the external value
  long LRsnr;   //number of lines of the matrix
  double *LRv;  //vector of interpolatedvalues
  double **LRc; //cyclic values from the parameter file (one vector for each LR variable)
  long *LRcnc;  //number of components of the vector (for each component)
  double *LRd;  //vector of default values

  double **qins;
  double *qinv;
  long qinsnr;
  long qinline;

  double tau_cloud;
  double tau_cloud_av;
  short tau_cloud_yes;
  short tau_cloud_av_yes;

  DOUBLEMATRIX *Tgrid;
  DOUBLEMATRIX *Pgrid;
  DOUBLEMATRIX *Vgrid;
  DOUBLEMATRIX *Vdir;
  DOUBLEMATRIX *RHgrid;

  std::unique_ptr<Vector<double>> Tamean;
  std::unique_ptr<Vector<double>> Vspdmean;
  std::unique_ptr<Vector<double>> Vdirmean;
  std::unique_ptr<Vector<double>> RHmean;

  std::unique_ptr<Vector<double>> Taplot;
  std::unique_ptr<Vector<double>> Vxplot;
  std::unique_ptr<Vector<double>> Vyplot;
  std::unique_ptr<Vector<double>> RHplot;

  double V;

  DOUBLEMATRIX *Tday;
  DOUBLEMATRIX *Tvar;

  long nstsrad;
  long nstlrad;
  long nstcloud;
  long nstTs;
  long nstTbottom;

  LONGVECTOR *imeteo_stations;

} METEO;


class ALLDATA
{
  std::unique_ptr<SOIL> _S;
  std::unique_ptr<WATER> _W;
  std::unique_ptr<LAND> _L;
  std::unique_ptr<PAR> _P;
  std::unique_ptr<TOPO> _T;
  std::unique_ptr<CHANNEL> _C;
  std::unique_ptr<ENERGY> _E;
  std::unique_ptr<SNOW> _N;
  std::unique_ptr<GLACIER> _G;
  std::unique_ptr<METEO> _M;
  std::unique_ptr<TIMES> _I;
public:
  SOIL       *S;
  WATER      *W;
  LAND       *L;
  PAR        *P;
  TOPO       *T;
  CHANNEL    *C;
  ENERGY     *E;
  SNOW       *N;
  GLACIER    *G;
  METEO      *M;
  TIMES      *I;
  ALLDATA():
    _S{new SOIL},
  _W{new WATER},
  _L{new LAND},
  _P{new PAR},
  _T{new TOPO},
  _C{new CHANNEL},
  _E{new ENERGY},
  _N{new SNOW},
  _G{new GLACIER},
  _M{new METEO},
  _I{new TIMES},
  S{_S.get()},
  W{_W.get()},
  L{_L.get()},
  P{_P.get()},
  T{_T.get()},
  C{_C.get()},
  E{_E.get()},
  N{_N.get()},
  G{_G.get()},
  M{_M.get()},
  I{_I.get()} {}
};
