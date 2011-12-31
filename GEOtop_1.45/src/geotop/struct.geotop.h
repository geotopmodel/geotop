
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
 #include "../libraries/fluidturtle/turtle.h"
 #include "../libraries/fluidturtle/t_datamanipulation.h"
 #include "../libraries/fluidturtle/t_utilities.h"
 #include "../libraries/fluidturtle/tensor3D.h"
 //#include "turtle2umfpack.h"

 
/*---------------------------------------------------------------------------*/
typedef struct {

    DOUBLEMATRIX *Rn_mean; 
	DOUBLEMATRIX *Rn_max;
	DOUBLEMATRIX *Rn_min;   
	DOUBLEMATRIX *LWin_mean;
	DOUBLEMATRIX *LW_mean;
	DOUBLEMATRIX *LW_max;
	DOUBLEMATRIX *LW_min;	
	DOUBLEMATRIX *SW_mean;	
	DOUBLEMATRIX *SW_max;			
    DOUBLEMATRIX *ET_mean; 
	DOUBLEMATRIX *ET_max;
	DOUBLEMATRIX *ET_min;
    DOUBLEMATRIX *H_mean; 
	DOUBLEMATRIX *H_max;
	DOUBLEMATRIX *H_min;    
    DOUBLEMATRIX *SEB_mean; 
	DOUBLEMATRIX *G_max;
	DOUBLEMATRIX *G_min;
	DOUBLEMATRIX *G_snowsoil;	
    DOUBLEMATRIX *Ts_mean;  /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
	DOUBLEMATRIX *Ts_max;
	DOUBLEMATRIX *Ts_min;
	DOUBLEMATRIX *Rswdown_mean;
	DOUBLEMATRIX *Rswdown_max;
	
	DOUBLEMATRIX *Rswbeam_mean;
	LONGMATRIX *nDt_shadow;
	LONGMATRIX *nDt_sun;
	
	DOUBLEMATRIX *Hgplot;
	DOUBLEMATRIX *LEgplot;
	DOUBLEMATRIX *Hvplot;
	DOUBLEMATRIX *LEvplot;	
	
	DOUBLEMATRIX *SWinplot;
	DOUBLEMATRIX *SWgplot;
	DOUBLEMATRIX *SWvplot;
	
	DOUBLEMATRIX *LWinplot;
	DOUBLEMATRIX *LWgplot;
	DOUBLEMATRIX *LWvplot;
	
	DOUBLEMATRIX *Tgplot;
	DOUBLEMATRIX *Tvplot;
	DOUBLEMATRIX *Tsplot;	
	
	double *sun;
	DOUBLEVECTOR *hsun;
	DOUBLEVECTOR *sinhsun;
	DOUBLEVECTOR *dsun;

	DOUBLEVECTOR *Dlay;
	DOUBLEVECTOR *wliq;
	DOUBLEVECTOR *wice;
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
	
	DOUBLEMATRIX *Tgskin;

	DOUBLEMATRIX *Tgskinsurr;
	DOUBLEMATRIX *Asurr;

} ENERGY;


/*---------------------------------------------------------------------------*/
typedef struct {

	LONGMATRIX *type;
	DOUBLEVECTOR *init_water_table_height;
	DOUBLETENSOR *pa;
	DOUBLETENSOR *pa_bed;
	DOUBLETENSOR *P;
	DOUBLETENSOR *Ptot;
	DOUBLETENSOR *T;
	DOUBLETENSOR *T_av_tensor;
	DOUBLETENSOR *thice;
	DOUBLETENSOR *th;
	DOUBLEMATRIX *Jinf;
	DOUBLEMATRIX *Tv;
	SHORTMATRIX *bc;
	DOUBLETENSOR *ET;
	
	DOUBLEMATRIX *Tzplot;
	DOUBLEMATRIX *Tzavplot;
	DOUBLEMATRIX *Ptotzplot;
	DOUBLEMATRIX *Pzplot;
	DOUBLEMATRIX *thzplot;
	DOUBLEMATRIX *thzavplot;
	DOUBLEMATRIX *thicezplot;
	DOUBLEMATRIX *thicezavplot;
	
	
} SOIL;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *Z0;         //elevetions of each pixel (DEM)
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
	
	LONGMATRIX *Rdown;
	LONGMATRIX *Cdown;
	
	SHORTMATRIX *is_on_border;
	
	DOUBLEMATRIX *East;
	DOUBLEMATRIX *North;
	
	LONGMATRIX *BC_counter;
	DOUBLEVECTOR *BC_LatDistance;
	DOUBLEVECTOR *BC_DepthFreeSurface;
	
	DOUBLEMATRIX *dzdE;
	DOUBLEMATRIX *dzdN;
	
	DOUBLEMATRIX *latitude;
	DOUBLEMATRIX *longitude;
			
} TOPO;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *LC;            //land cover (wood,lake,town,...) for each pixel*/
	SHORTMATRIX *shadow;		  //=1 if shadow, =0 if not*/
	DOUBLEMATRIX *ty;
	
	double ***vegpars;
	double **vegparv;
	DOUBLEVECTOR *vegpar;
	long *NumlinesVegTimeDepData;
	
	DOUBLEMATRIX *root_fraction;
	
} LANDCOVER;/*all this data are calculated on the basis of land use data and some other par*/


/*---------------------------------------------------------------------------*/
typedef struct {/*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                  R=number of rows of the basin,C=number of columns in the basin*/
    LONGVECTOR *r;          /*array of rows of the channel-pixels; dimension=nch*/
    LONGVECTOR *c;          /*array of columns of the channel-pixels; dimension=nch*/
	LONGMATRIX *ch;
	
	LONGVECTOR *ch_down;
	
	DOUBLEVECTOR *Vsup;
	DOUBLEVECTOR *Vsub;	
	//DOUBLEVECTOR *Vsup_cum;
	//DOUBLEVECTOR *Vsub_cum;	
	
	DOUBLEVECTOR *h_sup;
	
	DOUBLEVECTOR *length;
	
	double Vout;
	
	long **ch3;
	LONGMATRIX *lch;
	
	LONGVECTOR *soil_type;
	DOUBLEMATRIX *P;
	DOUBLEMATRIX *th;
	DOUBLEMATRIX *thice;
	DOUBLEMATRIX *T;
	DOUBLEMATRIX *ET;
	DOUBLEVECTOR *Tgskin;
			
} CHANNEL;




/*---------------------------------------------------------------------------*/
typedef struct { /*nstations=number of all the rain-stations,number_of_pixels=number of all the pixels of the basin R*C,
                   R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
    DOUBLETENSOR *PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/
    DOUBLEMATRIX *Pnet;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                              same subroutine and in "water.balance.c" module*/
    DOUBLEMATRIX *wcan_rain;       /*intercepted precipitation in mm*/
    DOUBLEMATRIX *wcan_snow;       /*intercepted precipitation in mm*/

    DOUBLEMATRIX *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
    DOUBLEMATRIX *PrSNW_mean;

	DOUBLEMATRIX *h_sup;
		
	DOUBLEMATRIX *error;
		
	//UMFPACK_REAL_TRIPLET *Jtriplet;
	//UMFPACK_REAL_MATRIX *Jmatrix;	
	
	DOUBLEVECTOR *Lx;
	DOUBLEVECTOR *Ux;

	DOUBLEVECTOR *P0;
	DOUBLEVECTOR *H0;
	DOUBLEVECTOR *H1;
	DOUBLEVECTOR *dH;
	DOUBLEVECTOR *B;
	DOUBLEVECTOR *f;
	DOUBLEVECTOR *df;
	DOUBLEMATRIX *Klat;
	
	double Voutland;
	
		
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
 	long snowlayer_max;
	long snowlayer_inf;
	DOUBLEVECTOR *Dmin;
	DOUBLEVECTOR *Dmax;
	double Sr_glac;
	long glaclayer_max;
	long glaclayer_inf;
	DOUBLEVECTOR *Dmin_glac;
	DOUBLEVECTOR *Dmax_glac;

	short state_turb;
	short state_lwrad;	

	double imp;
	double f_bound_Richards;

	double epsilon_snow;
	
	double output_soil;
	double output_snow;
	double output_glac;
	double output_surfenergy;
	double output_vegetation;
	double output_meteo;
	
	DOUBLEMATRIX *chkpt;
	LONGMATRIX *rc;
	short state_px_coord;
		
	short recover;

	double Vmin;
		
	double snowcorrfact;
	double raincorrfact;
	
	double RHmin;
	
	short format_out;
	short sky;
	
	DOUBLEVECTOR *saving_points;
	
	//double Vis; //visibility in km (>5 km)
	//double Lozone; //thickness of the stratospheric ozone layer (in cm normal conditions)
	
	short point_sim;/* =0 distributed simulation, =1 point simulation (the parameter files are different in the two cases) */
			
	double snow_maxpor;
	double snow_density_cutoff;
	double drysnowdef_rate;
	double wetsnowdef_rate;
	double snow_viscosity;
	
	double latitude;
	double longitude;
	
	double z0_snow;
	long n_landuses;
		
	short meteodistr;
	short blowing_snow;
	
	LONGVECTOR *r_points;
	LONGVECTOR *c_points;
	
	double psimin;
	double stmin;
		    	
	short wat_balance;
	short en_balance;
	
	short distr_stat;
	
	double Dpsi;
	double dtmin;
	
	double PsiInf;
	double TolPsiInf;
	
	double q1;
	double q2;

	double ***transect;
	double **vtrans;
	
	LONGVECTOR *cont_trans;
	LONGVECTOR *ibeg;
	
	long nLC;	
	
	double fetch_up;
	double fetch_down;
	
	short iobsint;
	double dn;
	double slopewt;
	double curvewt;	
	
	short LRflag;
	
	SHORTVECTOR *vegflag;
	
	short harm_or_arit_mean_normal;
	short harm_or_arit_mean_parallel;
	
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
	
	short default_point;
	short default_basin;
	short default_snow;
	short default_glac;
	short default_soil;
	
	double Wice_PBSM;
	
	DOUBLEMATRIX *maxSWE;
		
	long soil_type_land_default;
	long soil_type_chan_default;
	
	double MinIncrFactWithElev;
	double MaxIncrFactWithElev;
	
	long nsurface;
	
	double max_courant_land;
	double max_courant_channel;
	double min_hsup_land;
	double min_hsup_channel;
	double dtmin_sup;
	
	long nsoiltypes;
	
	LONGVECTOR *IDpoint;
	
	double min_lambda_en;
	long max_times_min_lambda_en;
	short exit_lambda_min_en;
	long max_times_halving_time_step_en;
	
	double min_lambda_wat;
	long max_times_min_lambda_wat;
	short exit_lambda_min_wat;	
	long max_times_halving_time_step_wat;
		
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
	
	DOUBLEMATRIX *nondimens_age;
	DOUBLEMATRIX *dimens_age;	
	DOUBLEMATRIX *max;
	DOUBLEMATRIX *average;	
	DOUBLEMATRIX *MELTED;
	DOUBLEMATRIX *SUBL;
	DOUBLEMATRIX *t_snow;
		
	DOUBLEMATRIX *rho_newsnow;
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
	/*DOUBLEMATRIX *Qtrans_plot;
	DOUBLEMATRIX *Qtrans_eq_plot;
	DOUBLEMATRIX *Qsub_plot;
	DOUBLEMATRIX *Qsub_eq_plot;*/
	DOUBLEMATRIX *ListonSWE;
	DOUBLEMATRIX *softSWE;
	DOUBLEMATRIX *softSWE1;	
	DOUBLEMATRIX *Dplot;	
	LONGVECTOR *change_dir_wind;
	
} SNOW;

typedef struct { 
	
	STATEVAR_3D *G;
	
	DOUBLEMATRIX *MELTED;
	DOUBLEMATRIX *SUBL;	
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
	double *LRv;	//vector of interpolatedvalues (used by meteodistr)
	double **LRc;	//cyclic values from the parameter file (one vector for each LR variable)
	long *LRcnc;	//number of components of the vector (for each component)
	double *LRd;	//vector of default values
	
	DOUBLEVECTOR *tau_cloud;
	DOUBLEVECTOR *tau_cloud_av;
	SHORTVECTOR *tau_cloud_yes;
	SHORTVECTOR *tau_cloud_av_yes;
	
	DOUBLETENSOR *Tgrid;
	DOUBLETENSOR *Pgrid;
	DOUBLETENSOR *Vgrid;
	DOUBLETENSOR *Vdir;
	DOUBLETENSOR *RHgrid;
	
	DOUBLEMATRIX *Vspdmean;
	DOUBLEMATRIX *Vdirmean;
	DOUBLEMATRIX *RHmean;
	
	DOUBLEMATRIX *Ta_mean;
	DOUBLEMATRIX *Ta_max;
	DOUBLEMATRIX *Ta_min;	
	
	DOUBLEMATRIX *Taplot;
	DOUBLEMATRIX *Vspdplot;
	DOUBLEMATRIX *Vdirplot;
	DOUBLEMATRIX *RHplot;
	
	double V;
		
	DOUBLEMATRIX *Tday;	
	DOUBLEMATRIX *Tvar;
		
	long nstsrad;
	long nstlrad;
	long nstcloud;
	
	LONGVECTOR *imeteo_stations;
	
} METEO;

#ifdef USE_NETCDF

typedef struct {
	DOUBLETENSOR * soil_thw_cum; // cumulated version of S->th
	DOUBLETENSOR * soil_thi_cum;// cumulated version of S->thice
	DOUBLETENSOR * soil_T_cum;// cumulated version of S->T
	DOUBLETENSOR * soil_P_cum;// cumulated version of S->P
	DOUBLETENSOR * soil_Ptot_cum;// cumulated version of S->Ptot
	// snow
	DOUBLETENSOR* snowD_cum; // cumulated version of snow->S->Dzl
	DOUBLETENSOR* snowT_cum; // cumulated version of snow->S->T
	DOUBLETENSOR* snowI_cum; // cumulated version of snow->S->w_ice
	DOUBLETENSOR* snowW_cum; // cumulated version of snow->S->w_liq
	// glacier
	DOUBLETENSOR* glacD_cum; // cumulated version of glac->G->Dzl
	DOUBLETENSOR* glacT_cum; // cumulated version of glac->G->T
	DOUBLETENSOR* glacI_cum; // cumulated version of glac->G->w_ice
	DOUBLETENSOR* glacW_cum; // cumulated version of glac->G->w_liq
} OUTPUT_NCDATA;

#endif

typedef struct {
	SOIL *S;
	WATER *W;
	LANDCOVER *L;
	PAR *P;
	TOPO *T;
	CHANNEL *C;
	ENERGY *E;
	SNOW *N;
	GLACIER *G;
	METEO *M;
	TIMES *I;
#ifdef USE_NETCDF
	int ncid; // pointer to netCDF archive file
	long counter_snow; // counter for time print of snow
	long counter_surface_energy; //conter for surface energy maps
	OUTPUT_NCDATA *outnc;
#endif
}ALLDATA;

#ifdef USE_HPC
typedef struct GCSTRUCT	//Struct containing subdomain ghost-cells adjacency data for MPI SEND/RECV commands
{
	int rank;
	char calltype[4];
	int top;
	int left;	
	int bottom;	
	int right;	
	struct GCSTRUCT *next;
};
typedef struct GCSTRUCT;

struct WORKAREA	//Struct containing local subdomain coords
{
	int rank;
	int top;
	int left;	
	int bottom;	
	int right;
};
typedef struct WORKAREA;
#endif
