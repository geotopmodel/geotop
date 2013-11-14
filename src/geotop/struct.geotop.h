
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 1.225-15
 
 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
 #include "turtle.h"
 #include "t_datamanipulation.h"
 #include "t_utilities.h"
 #include "tensor3D.h"

 
/*---------------------------------------------------------------------------*/
typedef struct {

    DOUBLEVECTOR *Rn_mean; 
	DOUBLEVECTOR *LWin_mean;
	DOUBLEVECTOR *LW_mean;
	DOUBLEVECTOR *SW_mean;	
    DOUBLEVECTOR *ET_mean; 
    DOUBLEVECTOR *H_mean; 
    DOUBLEVECTOR *SEB_mean; 
    DOUBLEVECTOR *Ts_mean;  /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
	DOUBLEVECTOR *Rswdown_mean;
	DOUBLEVECTOR *Rswbeam_mean;
	LONGVECTOR *nDt_shadow;
	LONGVECTOR *nDt_sun;
	
	DOUBLEVECTOR *Rn;
	DOUBLEVECTOR *LWin;
	DOUBLEVECTOR *LW;
	DOUBLEVECTOR *SW;
	DOUBLEVECTOR *LE;
	DOUBLEVECTOR *H;
	DOUBLEVECTOR *G;
	DOUBLEVECTOR *Ts;
	DOUBLEVECTOR *SWin;
	DOUBLEVECTOR *SWinb;
	SHORTVECTOR *shad;
	
	DOUBLEVECTOR *Hgplot;
	DOUBLEVECTOR *LEgplot;
	DOUBLEVECTOR *Hvplot;
	DOUBLEVECTOR *LEvplot;	
	DOUBLEVECTOR *SWinplot;
	DOUBLEVECTOR *SWgplot;
	DOUBLEVECTOR *SWvplot;	
	DOUBLEVECTOR *LWinplot;
	DOUBLEVECTOR *LWgplot;
	DOUBLEVECTOR *LWvplot;
	DOUBLEVECTOR *Tgplot;
	DOUBLEVECTOR *Tsplot;	
	DOUBLEVECTOR *Tvplot;

	DOUBLEVECTOR *Hgp;
	DOUBLEVECTOR *LEgp;
	DOUBLEVECTOR *Hvp;
	DOUBLEVECTOR *LEvp;	
	DOUBLEVECTOR *SWinp;
	DOUBLEVECTOR *SWgp;
	DOUBLEVECTOR *SWvp;	
	DOUBLEVECTOR *LWinp;
	DOUBLEVECTOR *LWgp;
	DOUBLEVECTOR *LWvp;
	DOUBLEVECTOR *Tgp;
	DOUBLEVECTOR *Tsp;	
	
	double *sun;
	double hsun;
	double sinhsun;
	double dsun;

	DOUBLEVECTOR *Dlayer;
	DOUBLEVECTOR *liq;
	DOUBLEVECTOR *ice;
	DOUBLEVECTOR *Temp; 
	DOUBLEVECTOR *deltaw;
	DOUBLEVECTOR *SWlayer;
	DOUBLEVECTOR *soil_transp_layer;
	DOUBLEVECTOR *dFenergy;
	DOUBLEVECTOR *udFenergy;
	DOUBLEVECTOR *Kth0;
	DOUBLEVECTOR *Kth1;
	DOUBLEVECTOR *Fenergy;
	DOUBLEVECTOR *Newton_dir;
	DOUBLEVECTOR *T0;
	DOUBLEVECTOR *T1;
	DOUBLEVECTOR *Tstar;
	DOUBLEVECTOR *THETA;
	DOUBLEVECTOR *soil_evap_layer_bare;
	DOUBLEVECTOR *soil_evap_layer_veg;
	
	DOUBLEMATRIX *Tgskin_surr;
	DOUBLEMATRIX *SWrefl_surr;
	
} ENERGY;

/*---------------------------------------------------------------------------*/

typedef struct {
	DOUBLEMATRIX *P;
	DOUBLEMATRIX *thi;
	DOUBLEMATRIX *T;
	DOUBLEMATRIX *C;
} SOIL_STATE;

/*---------------------------------------------------------------------------*/

typedef struct {
	
	DOUBLEVECTOR *Tv;
	DOUBLEVECTOR *wrain;       /*intercepted precipitation in mm*/
    DOUBLEVECTOR *wsnow;       /*intercepted precipitation in mm*/
	
} STATE_VEG;

/*---------------------------------------------------------------------------*/

typedef struct {
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
	
	DOUBLEVECTOR *Pnetcum;
	DOUBLEVECTOR *ETcum;
		
} SOIL;


	
/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *Z0;         //elevation of each pixel (DEM)
	DOUBLETENSOR *Z;

    DOUBLEMATRIX *sky;        //view factor (of the sky) for each pixel
    SHORTMATRIX *pixel_type; 
	
	//SHORTMATRIX *DD;		  //Drainage Directions for each pixel; ex matr_ev->slope*/
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
	DOUBLEVECTOR *BC_DepthFreeSurface;
		
	DOUBLEMATRIX *dzdE;
	DOUBLEMATRIX *dzdN;
	
	DOUBLEMATRIX *latitude;
	DOUBLEMATRIX *longitude;
			
} TOPO;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *LC;            //land cover for each pixel
	DOUBLEMATRIX *delay;
	SHORTMATRIX *shadow;		  //=1 if shadow, =0 if not
	DOUBLEMATRIX *ty;
	
	double ***vegpars;
	double **vegparv;
	DOUBLEVECTOR *vegpar;
	long *NumlinesVegTimeDepData;
	
	DOUBLEMATRIX *root_fraction;
	
} LAND;/*all this data are calculated on the basis of land use data and some other par*/


/*---------------------------------------------------------------------------*/


typedef struct {/*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                  R=number of rows of the basin,C=number of columns in the basin*/
    LONGVECTOR *r;          /*array of rows of the channel-pixels; dimension=nch*/
    LONGVECTOR *c;          /*array of columns of the channel-pixels; dimension=nch*/
	LONGMATRIX *ch;
	LONGVECTOR *ch_down;
	DOUBLEVECTOR *Vsup;
	DOUBLEVECTOR *Vsub;	
	DOUBLEVECTOR *h_sup;
	DOUBLEVECTOR *length;
	double Vout;
	long **ch3;
	LONGMATRIX *lch;
	LONGVECTOR *soil_type;
	DOUBLEMATRIX *th;
	DOUBLEMATRIX *ET;
	DOUBLEVECTOR *Kbottom;
	SOIL_STATE *SS;
} CHANNEL;




/*---------------------------------------------------------------------------*/

typedef struct { /*nstations=number of all the rain-stations,number_of_pixels=number of all the pixels of the basin R*C,
                   R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
    DOUBLEMATRIX *PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/
    DOUBLEMATRIX *Pnet;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                              same subroutine and in "water.balance.c" module*/

    DOUBLEVECTOR *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
    DOUBLEVECTOR *PrSNW_mean;
	DOUBLEVECTOR *Pt;
	DOUBLEVECTOR *Ps;

	DOUBLEVECTOR *h_sup;
		
	DOUBLEMATRIX *error;
	
	DOUBLEVECTOR *Lx;
	DOUBLEVECTOR *Ux;

	DOUBLEVECTOR *H0;
	DOUBLEVECTOR *H1;
	DOUBLEVECTOR *dH;
	DOUBLEVECTOR *B;
	DOUBLEVECTOR *f;
	DOUBLEVECTOR *df;
	DOUBLEMATRIX *Klat;
	DOUBLEMATRIX *Kbottom;
	
	double Voutlandsub;
	double Voutlandsup;
	double Voutbottom;
	
} WATER;

/*---------------------------------------------------------------------------*/
typedef struct {
	
	DOUBLEVECTOR *JD_plots;
	double time;    /*time=current time from the begin of simulation [s]*/
	long iplot;
	double **Dt_matrix;
	long numlinesDt_matrix;
	double *Dt_vector;

} TIMES;


/*---------------------------------------------------------------------------*/
typedef struct {
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
    double Sr;		  /*WATER FRACTION RETAINED BY CAPILLARY FORCES IN SNOW*/
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
	
	DOUBLEVECTOR *output_soil;
	DOUBLEVECTOR *output_snow;
	DOUBLEVECTOR *output_glac;
	DOUBLEVECTOR *output_surfenergy;
	DOUBLEVECTOR *output_vegetation;
	DOUBLEVECTOR *output_meteo;
	
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
	
	DOUBLEVECTOR *saving_points;
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
	short transport_model;	//by Flo
				
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
	
    DOUBLEVECTOR *Dtplot_point;  
	DOUBLEVECTOR *Dtplot_basin;
	DOUBLEVECTOR *Dtplot_discharge;
	
	short state_pixel;
	short state_discharge;
	short state_basin;
	
	double Dt_PBSM;
	
	long lowpass;
	long lowpass_curvatures;
	
	short dew;
	
	DOUBLEVECTOR *init_date;
	DOUBLEVECTOR *end_date;
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
	
	DOUBLEVECTOR *soil_plot_depths;
	DOUBLEVECTOR *snow_plot_depths;
	DOUBLEVECTOR *glac_plot_depths;
	
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
typedef struct {
	SHORTMATRIX  *type;
	LONGMATRIX	 *lnum;
	DOUBLETENSOR *Dzl;
	DOUBLETENSOR *w_liq;
	DOUBLETENSOR *w_ice;
	DOUBLETENSOR *T;
} STATEVAR_3D;


typedef struct {
	short type;
	long lnum;
	DOUBLEVECTOR *Dzl;
	DOUBLEVECTOR *w_liq;
	DOUBLEVECTOR *w_ice;
	DOUBLEVECTOR *T;
} STATEVAR_1D;

typedef struct {
	STATEVAR_3D *S;
	STATEVAR_1D *S_for_BS;
	DOUBLEVECTOR *age;	
	DOUBLEVECTOR *MELTED;
	DOUBLEVECTOR *melted;
	DOUBLEVECTOR *SUBL;
	DOUBLEVECTOR *subl;
	DOUBLEVECTOR *t_snow;
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
	DOUBLEVECTOR *Dplot;	
	LONGVECTOR *change_dir_wind;
} SNOW;

typedef struct { 
	STATEVAR_3D *G;
	DOUBLEVECTOR *MELTED;
	DOUBLEVECTOR *melted;
	DOUBLEVECTOR *SUBL;
	DOUBLEVECTOR *subl;
} GLACIER;

typedef struct{
	DOUBLEVECTOR *E;
	DOUBLEVECTOR *N;
	DOUBLEVECTOR *lat;
	DOUBLEVECTOR *lon;
	DOUBLEVECTOR *Z;
	DOUBLEVECTOR *sky;
	DOUBLEVECTOR *ST;
	DOUBLEVECTOR *Vheight;
	DOUBLEVECTOR *Theight;
} METEO_STATIONS;

typedef struct {
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
	
	double **LRs;	//matrix read from the external value
	long LRsnr;		//number of lines of the matrix
	double *LRv;	//vector of interpolatedvalues
	double **LRc;	//cyclic values from the parameter file (one vector for each LR variable)
	long *LRcnc;	//number of components of the vector (for each component)
	double *LRd;	//vector of default values
	
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
	
	DOUBLEVECTOR *Tamean;
	DOUBLEVECTOR *Vspdmean;
	DOUBLEVECTOR *Vdirmean;
	DOUBLEVECTOR *RHmean;
	
	DOUBLEVECTOR *Taplot;
	DOUBLEVECTOR *Vxplot;
	DOUBLEVECTOR *Vyplot;
	DOUBLEVECTOR *RHplot;
		
	double V;
		
	DOUBLEMATRIX *Tday;	
	DOUBLEMATRIX *Tvar;
		
	long nstsrad;
	long nstlrad;
	long nstcloud;
	long nstTs;
	
	LONGVECTOR *imeteo_stations;
	
} METEO;


typedef struct {
	SOIL *S;
	WATER *W;
	TRANSPORT *Tr;	//by Flo
	LAND *L;
	PAR *P;
	TOPO *T;
	CHANNEL *C;
	ENERGY *E;
	SNOW *N;
	GLACIER *G;
	METEO *M;
	TIMES *I;	
} ALLDATA;
