
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef STRUCT_GEOTOP_H
#define STRUCT_GEOTOP_H
//#include "../libraries/fluidturtle/turtle.h"
#include "../libraries/fluidturtle/tensors3D.h"
#include "datastructs.h"

#include <vector>
#include <meteoio/MeteoIO.h>
 
/*---------------------------------------------------------------------------*/

//typedef struct {
class Energy {
 public:
//	DOUBLEVECTOR *Rn_mean;
	GeoVector<double> Rn_mean;
//	DOUBLEVECTOR *LWin_mean;
	GeoVector<double> LWin_mean;
//	DOUBLEVECTOR *LW_mean;
	GeoVector<double> LW_mean;
//	DOUBLEVECTOR *SW_mean;
	GeoVector<double> SW_mean;
//	DOUBLEVECTOR *ET_mean;
	GeoVector<double> ET_mean;
//	DOUBLEVECTOR *H_mean;
	GeoVector<double> H_mean;
//	DOUBLEVECTOR *SEB_mean;
	GeoVector<double> SEB_mean;
//	DOUBLEVECTOR *Ts_mean;  /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
    GeoVector<double> Ts_mean;  /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
//	DOUBLEVECTOR *Rswdown_mean;
    GeoVector<double> Rswdown_mean;
//	DOUBLEVECTOR *Rswbeam_mean;
    GeoVector<double> Rswbeam_mean;
//	LONGVECTOR *nDt_shadow;
	GeoVector<long> nDt_shadow;
//	LONGVECTOR *nDt_sun;
	GeoVector<long> nDt_sun;
//	DOUBLEVECTOR *Rn;
	GeoVector<double> Rn;
//	DOUBLEVECTOR *LWin;
    GeoVector<double> LWin;
//	DOUBLEVECTOR *LW;
	GeoVector<double> LW;
//	DOUBLEVECTOR *SW;
	GeoVector<double> SW;
//	DOUBLEVECTOR *LE;
	GeoVector<double> LE;
//	DOUBLEVECTOR *H;
	GeoVector<double> H;
//	DOUBLEVECTOR *G;
	GeoVector<double> G;
//	DOUBLEVECTOR *Ts;
	GeoVector<double> Ts;
//	DOUBLEVECTOR *SWin;
	GeoVector<double> SWin;
//	DOUBLEVECTOR *SWinb;
	GeoVector<double> SWinb;
//	SHORTVECTOR *shad;
	GeoVector<short> shad;
//	DOUBLEVECTOR *Hgplot;
	GeoVector<double> Hgplot;
//	DOUBLEVECTOR *LEgplot;
	GeoVector<double> LEgplot;
//	DOUBLEVECTOR *Hvplot;
	GeoVector<double> Hvplot;
//	DOUBLEVECTOR *LEvplot;
	GeoVector<double> LEvplot;
//	DOUBLEVECTOR *SWinplot;
	GeoVector<double> SWinplot;
//	DOUBLEVECTOR *SWgplot;
	GeoVector<double> SWgplot;
//	DOUBLEVECTOR *SWvplot;
	GeoVector<double> SWvplot;
//	DOUBLEVECTOR *LWinplot;
	GeoVector<double> LWinplot;
//	DOUBLEVECTOR *LWgplot;
	GeoVector<double> LWgplot;
//	DOUBLEVECTOR *LWvplot;
	GeoVector<double> LWvplot;
//	DOUBLEVECTOR *Tgplot;
	GeoVector<double> Tgplot;
//	DOUBLEVECTOR *Tsplot;
	GeoVector<double> Tsplot;
//	DOUBLEVECTOR *Tvplot;
	GeoVector<double> Tvplot;
//	DOUBLEVECTOR *Hgp;
	GeoVector<double> Hgp;
//	DOUBLEVECTOR *LEgp;
	GeoVector<double> LEgp;
//	DOUBLEVECTOR *Hvp;
	GeoVector<double> Hvp;
//	DOUBLEVECTOR *LEvp;
	GeoVector<double> LEvp;
//	DOUBLEVECTOR *SWinp;
	GeoVector<double> SWinp;
//	DOUBLEVECTOR *SWgp;
	GeoVector<double> SWgp;
//	DOUBLEVECTOR *SWvp;
	GeoVector<double> SWvp;
//	DOUBLEVECTOR *LWinp;
	GeoVector<double> LWinp;
//	DOUBLEVECTOR *LWgp;
	GeoVector<double> LWgp;
//	DOUBLEVECTOR *LWvp;
	GeoVector<double> LWvp;
//	DOUBLEVECTOR *Tgp;
	GeoVector<double> Tgp;
//	DOUBLEVECTOR *Tsp;
	GeoVector<double> Tsp;
	
	double *sun;
	double hsun;
	double sinhsun;
	double dsun;

//	DOUBLEVECTOR *Dlayer;
	GeoVector<double> Dlayer;
//	DOUBLEVECTOR *liq;
	GeoVector<double> liq;
//	DOUBLEVECTOR *ice;
	GeoVector<double> ice;
//	DOUBLEVECTOR *Temp;
	GeoVector<double> Temp;
//	DOUBLEVECTOR *deltaw;
	GeoVector<double> deltaw;
//	DOUBLEVECTOR *SWlayer;
	GeoVector<double> SWlayer;
//	DOUBLEVECTOR *soil_transp_layer;
	GeoVector<double> soil_transp_layer;
//	DOUBLEVECTOR *dFenergy;
	GeoVector<double> dFenergy;
//	DOUBLEVECTOR *udFenergy;
	GeoVector<double> udFenergy;
//	DOUBLEVECTOR *Kth0;
	GeoVector<double> Kth0;
//	DOUBLEVECTOR *Kth1;
	GeoVector<double> Kth1;
//	DOUBLEVECTOR *Fenergy;
	GeoVector<double> Fenergy;
//	DOUBLEVECTOR *Newton_dir;
	GeoVector<double> Newton_dir;
//	DOUBLEVECTOR *T0;
	GeoVector<double> T0;
//	DOUBLEVECTOR *T1;
	GeoVector<double> T1;
//	DOUBLEVECTOR *Tstar;
	GeoVector<double> Tstar;
//	DOUBLEVECTOR *THETA;
	GeoVector<double> THETA;
//	DOUBLEVECTOR *soil_evap_layer_bare;
	GeoVector<double> soil_evap_layer_bare;
//	DOUBLEVECTOR *soil_evap_layer_veg;
	GeoVector<double> soil_evap_layer_veg;
//	DOUBLEMATRIX *Tgskin_surr;
	GeoMatrix<double> Tgskin_surr;
//	DOUBLEMATRIX *SWrefl_surr;
	GeoMatrix<double> SWrefl_surr;
	};
//	} ENERGY;

/*---------------------------------------------------------------------------*/

//typedef struct {
class SoilState{
public:
//	DOUBLEMATRIX *P;
	GeoMatrix<double> P;
//	DOUBLEMATRIX *thi;
	GeoMatrix<double> thi;
//	DOUBLEMATRIX *T;
	GeoMatrix<double> T;

//	GeoMatrix<double> *P1;
	};

//	} SOIL_STATE;

/*---------------------------------------------------------------------------*/

//typedef struct {
class StateVeg{
public:
//	DOUBLEVECTOR *Tv;
	GeoVector<double> Tv;
//	DOUBLEVECTOR *wrain;       /*intercepted precipitation in mm*/
	GeoVector<double> wrain;   /*intercepted precipitation in mm*/
//	DOUBLEVECTOR *wsnow;       /*intercepted precipitation in mm*/
	GeoVector<double> wsnow;   /*intercepted precipitation in mm*/
};
//	} STATE_VEG;

/*---------------------------------------------------------------------------*/

//typedef struct {
class Soil{
public:
//	LONGMATRIX *type;
	GeoMatrix<long> type;
//	DOUBLEVECTOR *init_water_table_depth;
	GeoVector<double> init_water_table_depth;
//	DOUBLETENSOR *pa;
	GeoTensor<double> pa;
	GeoTensor<double> pa_bed;
//	DOUBLEMATRIX *T_av_tensor;
	GeoMatrix<double> T_av_tensor;
//	DOUBLEMATRIX *thw_av_tensor;
	GeoMatrix<double> thw_av_tensor;
//	DOUBLEMATRIX *thi_av_tensor;
	GeoMatrix<double> thi_av_tensor;
//	DOUBLEMATRIX *Ptot;
	GeoMatrix<double> Ptot;
//	DOUBLEMATRIX *th;
	GeoMatrix<double> th;
//	DOUBLETENSOR *ET;
	GeoTensor<double> ET;
//	DOUBLEMATRIX *Tzplot;
	GeoMatrix<double> Tzplot;
//	DOUBLEMATRIX *Tzavplot;
	GeoMatrix<double> Tzavplot;
//	DOUBLEMATRIX *Ptotzplot;
	GeoMatrix<double> Ptotzplot;
//	DOUBLEMATRIX *Pzplot;
	GeoMatrix<double> Pzplot;
//	DOUBLEMATRIX *thzplot;
	GeoMatrix<double> thzplot;
//	DOUBLEMATRIX *thzavplot;
	GeoMatrix<double> thzavplot;
//	DOUBLEMATRIX *thizplot;
	GeoMatrix<double> thizplot;
//	DOUBLEMATRIX *thizavplot;
	GeoMatrix<double> thizavplot;
//	SOIL_STATE *SS;
	SoilState *SS;
//	STATE_VEG *VS;
	StateVeg *VS;
};
//	} SOIL;

	
/*---------------------------------------------------------------------------*/

//typedef struct {
class Topo{
	public:
//  DOUBLEMATRIX *Z0;         //elevetions of each pixel (DEM)
	GeoMatrix<double> Z0;	  //elevetions of each pixel (DEM)
//	DOUBLETENSOR *Z;
    GeoTensor<double> Z;
//  DOUBLEMATRIX *sky;        //view factor (of the sky) for each pixel
	GeoMatrix<double> sky;        //view factor (of the sky) for each pixel
//  SHORTMATRIX *pixel_type;
    GeoMatrix<short> pixel_type;

//  DOUBLEMATRIX *aspect;     	  /*aspect; ex: matr_ev->azimuth*/
    GeoMatrix<double> aspect;     /*aspect; ex: matr_ev->azimuth*/
//  DOUBLEMATRIX *slope;          /*slope of the pixels; ex: matr_ev->slope*/
	GeoMatrix<double>  slope;     /*slope of the pixels; ex: matr_ev->slope*/

//	DOUBLEMATRIX *curvature1;
    GeoMatrix<double> curvature1;
//	DOUBLEMATRIX *curvature2;
    GeoMatrix<double> curvature2;
//	DOUBLEMATRIX *curvature3;
    GeoMatrix<double> curvature3;
//	DOUBLEMATRIX *curvature4;
    GeoMatrix<double> curvature4;
	
	double ***horizon_height;
	long *horizon_numlines;
//	LONGMATRIX *horizon_point;
	GeoMatrix<long> horizon_point;

	long num_horizon_point;
		
	long ***i_cont;
//	LONGMATRIX *lrc_cont;
	GeoMatrix<long> lrc_cont;
	
	long **j_cont;
//	LONGMATRIX *rc_cont;
	GeoMatrix<long> rc_cont;
		
//	LONGVECTOR *Lp;
	GeoVector<long> Lp;
//	LONGVECTOR *Li;
	GeoVector<long> Li;
	
//	LONGMATRIX *Jdown;
	GeoMatrix<long> Jdown;
//	DOUBLEMATRIX *Qdown;
	GeoMatrix<double> Qdown;

//	SHORTMATRIX *is_on_border;
	GeoMatrix<short> is_on_border;

//	DOUBLEMATRIX *East;
	GeoMatrix<double> East;
//	DOUBLEMATRIX *North;
	GeoMatrix<double> North;
	
//	LONGMATRIX *BC_counter;
	GeoMatrix<long> BC_counter;

//	DOUBLEVECTOR *BC_DepthFreeSurface;
	GeoVector<double> BC_DepthFreeSurface;

//	DOUBLEMATRIX *dzdE;
	GeoMatrix<double> dzdE;
//	DOUBLEMATRIX *dzdN;
	GeoMatrix<double> dzdN;
	
//	DOUBLEMATRIX *latitude;
	GeoMatrix<double> latitude;
//	DOUBLEMATRIX *longitude;
	GeoMatrix<double> longitude;
};
//} TOPO;


/*---------------------------------------------------------------------------*/
//typedef struct {
class Land{
public:
//  DOUBLEMATRIX *LC;            //land cover (wood,lake,town,...) for each pixel*/
    GeoMatrix<double> LC;        //land cover (wood,lake,town,...) for each pixel*/
//	DOUBLEMATRIX *delay;
	GeoMatrix<double> delay;
//	SHORTMATRIX *shadow;		  //=1 if shadow, =0 if not*/
	GeoMatrix<short> shadow; 	  //=1 if shadow, =0 if not*/
//	DOUBLEMATRIX *ty;
	GeoMatrix<double> ty;
	
	double ***vegpars;
	double **vegparv;
//	DOUBLEVECTOR *vegpar;
	GeoVector<double> vegpar;

	long *NumlinesVegTimeDepData;
	
//	DOUBLEMATRIX *root_fraction;
	GeoMatrix<double> root_fraction;
}; /*all this data are calculated on the basis of land use data and some other par*/
//} LAND;/*all this data are calculated on the basis of land use data and some other par*/


/*---------------------------------------------------------------------------*/


//typedef struct {
class Channel{
public:
	/*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                  R=number of rows of the basin,C=number of columns in the basin*/

//	LONGVECTOR *r;          /*array of rows of the channel-pixels; dimension=nch*/
	GeoVector<long> r;    /*array of rows of the channel-pixels; dimension=nch*/
//  LONGVECTOR *c;          /*array of columns of the channel-pixels; dimension=nch*/
    GeoVector<long> c;    /*array of columns of the channel-pixels; dimension=nch*/
//	LONGMATRIX *ch;
	GeoMatrix<long> ch;
//	LONGVECTOR *ch_down;
	GeoVector<long> ch_down;
//	DOUBLEVECTOR *Vsup;
	GeoVector<double> Vsup;
//	DOUBLEVECTOR *Vsub;
	GeoVector<double> Vsub;
//	DOUBLEVECTOR *h_sup;
	GeoVector<double> h_sup;
//	DOUBLEVECTOR *length;
	GeoVector<double> length;

	double Vout;
	long **ch3;
//	LONGMATRIX *lch;
	GeoMatrix<long> lch;
//	LONGVECTOR *soil_type;
	GeoVector<long> soil_type;

//	DOUBLEMATRIX *th;
	GeoMatrix<double> th;
//	DOUBLEMATRIX *ET;
	GeoMatrix<double> ET;

//	DOUBLEVECTOR *Kbottom;
	GeoVector<double> Kbottom;

//	SOIL_STATE *SS;
	SoilState *SS;
};
//	} CHANNEL;




/*---------------------------------------------------------------------------*/

//typedef struct {

class Water { /*nstations=number of all the rain-stations,number_of_pixels=number of all the pixels of the basin R*C,*/
    public:              /* R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
//  DOUBLEMATRIX *PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/
	GeoMatrix<double> PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/

//  DOUBLEMATRIX *Pnet;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              /*of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                                same subroutine and in "water.balance.c" module*/
    GeoMatrix<double> Pnet;   /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                                of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                                same subroutine and in "water.balance.c" module*/

//	DOUBLEVECTOR *PrTOT_mean;  	 /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
    GeoVector<double> PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
//	DOUBLEVECTOR *PrSNW_mean;
    GeoVector<double> PrSNW_mean;
//	DOUBLEVECTOR *Pt;
	GeoVector<double> Pt;
//	DOUBLEVECTOR *Ps;
	GeoVector<double> Ps;

//	DOUBLEVECTOR *h_sup;
    GeoVector<double> h_sup;

//	DOUBLEMATRIX *error;
	GeoMatrix<double> error;

//	DOUBLEVECTOR *Lx;
	GeoVector<double> Lx;

//	DOUBLEVECTOR *Ux;

//	DOUBLEVECTOR *H0;
	GeoVector<double> H0;
//	DOUBLEVECTOR *H1;
	GeoVector<double> H1;
//	DOUBLEVECTOR *dH;
	GeoVector<double> dH;

//	DOUBLEVECTOR *B;
	GeoVector<double> B;
//	DOUBLEVECTOR *f;
	GeoVector<double> f;
//	DOUBLEVECTOR *df;
	GeoVector<double> df;

//	DOUBLEMATRIX *Klat;
	GeoMatrix<double> Klat;
//	DOUBLEMATRIX *Kbottom;
	GeoMatrix<double> Kbottom;
	
	double Voutlandsub;
	double Voutlandsup;
	double Voutbottom;
};
//} WATER;


/*---------------------------------------------------------------------------*/
//typedef struct {
class Times{
public:
//	DOUBLEVECTOR *JD_plots;
	GeoVector<double> JD_plots;
	double time;    /*time=current time from the begin of simulation [s]*/
	long iplot;
	double **Dt_matrix;
	long numlinesDt_matrix;
	double *Dt_vector;
};
//} TIMES;


/*---------------------------------------------------------------------------*/
//typedef struct {
class Par{
public:
    double Dt;      	/*Dt=the integration time interval [s]*/
	double ST;
    short print;        /*1 IF YOU WANT TO PRINT MATRICES WITH INTERMEDIATE RESULTS, 0 OTHERWISE*/
    short monin_obukhov;  
    double gamma_m;   	/*Exponent of the law of uniform motion on the surface*/
    double T_rain;    	/*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
    double T_snow;    	/*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
    double aep;       	/*ALBEDO EXTINCTION PARAMETER [m]*/
    double avo;       	/*NEW SNOW VISIBLE BAND REFLECTANCE*/
    double airo;      	/*NEW NEAR INFRARED BAND REFLECTANCE*/
    double Sr;		  	/*WATER FRACTION RETAINED BY CAPILLARY FORCES IN SNOW*/
    double rho_ice;     /*Ice density [kg/mc]*/
    long total_pixel;   /*The number of the valid pixel of the whole basin*/
	long total_channel;
	double total_area;
	
	double max_weq_snow;
	long max_snow_layers;
//	LONGVECTOR *inf_snow_layers;
	GeoVector<long> inf_snow_layers;

	double max_weq_glac;
	long max_glac_layers;
//	LONGVECTOR *inf_glac_layers;
	GeoVector<long> inf_glac_layers;

	double Sr_glac;

	short state_turb;
	short state_lwrad;	

	double imp;
	double f_bound_Richards;

	double epsilon_snow;
	
//	DOUBLEVECTOR *output_soil;
	GeoVector<double> output_soil;
//	DOUBLEVECTOR *output_snow;
	GeoVector<double> output_snow;
//	DOUBLEVECTOR *output_glac;
	GeoVector<double> output_glac;
//	DOUBLEVECTOR *output_surfenergy;
	GeoVector<double> output_surfenergy;
//	DOUBLEVECTOR *output_vegetation;
	GeoVector<double> output_vegetation;
//	DOUBLEVECTOR *output_meteo;
	GeoVector<double> output_meteo;

	short output_soil_bin;
	short output_snow_bin;
	short output_glac_bin;
	short output_surfenergy_bin;
	short output_meteo_bin;
	
//	DOUBLEMATRIX *chkpt;
	GeoMatrix<double> chkpt;
//	LONGMATRIX *rc;
	GeoMatrix<long> rc;
//	LONGVECTOR *jplot;
	GeoVector<long> jplot;
		
	short recover;

	double Vmin;
		
	double snowcorrfact;
	double raincorrfact;
	
	double RHmin;
	
	short format_out;
	short sky;
	
//	DOUBLEVECTOR *saving_points;
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
	long n_landuses;
		
	short blowing_snow;
	
//	LONGVECTOR *r_points;
	GeoVector<long> r_points;
//	LONGVECTOR *c_points;
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
	
//	SHORTVECTOR *vegflag;
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
	
	short tsteps_from_file;
	
//	SHORTVECTOR *plot_discharge_with_Dt_integration;
	GeoVector<short> plot_discharge_with_Dt_integration;
//	SHORTVECTOR *plot_point_with_Dt_integration;
	GeoVector<short> plot_point_with_Dt_integration;
//	SHORTVECTOR *plot_basin_with_Dt_integration;
	GeoVector<short> plot_basin_with_Dt_integration;
	
//  DOUBLEVECTOR *Dtplot_point;
	GeoVector<double> Dtplot_point;
//	DOUBLEVECTOR *Dtplot_basin;
    GeoVector<double> Dtplot_basin;
//	DOUBLEVECTOR *Dtplot_discharge;
	GeoVector<double> Dtplot_discharge;

	short state_pixel;
	short state_discharge;
	short state_basin;
	
	double Dt_PBSM;
	
	long lowpass;
	long lowpass_curvatures;
	
	short dew;
	
//	DOUBLEVECTOR *init_date;
	GeoVector<double> init_date;
//	DOUBLEVECTOR *end_date;
	GeoVector<double> end_date;
//	LONGVECTOR *run_times;
	GeoVector<long> run_times;
	
	double delay_day_recover;
	
	short all_point;
	short all_basin;
	short all_snow;
	short all_glac;
	short all_soil;
	
	double Wice_PBSM;
	
//	DOUBLEMATRIX *maxSWE;
	GeoMatrix<double> maxSWE;
		
	long soil_type_land_default;
	long soil_type_chan_default;
	
	double MinIncrFactWithElev;
	double MaxIncrFactWithElev;
	
	long nsurface;
	
	double max_courant_land;
	double max_courant_channel;
	double min_hsup_land;
	double min_hsup_channel;
	double min_dhsup_land_channel;
	double dtmin_sup;
	
	long nsoiltypes;
	
//	LONGVECTOR *IDpoint;
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
	
//	SHORTVECTOR *linear_interpolation_meteo;
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
	short usemeteoio;// flag indicating whether MeteoIO library is used
	bool use_meteoio_cloud;
	short use_meteoio_meteodata;
	
	short qin;
	short flag1D;
	
	double k_to_ksat;
	short RunIfAnOldRunIsPresent;

	double Lozone;
	double alpha_iqbal;
	double beta_iqbal;
};
//} PAR;


//typedef struct {
class Statevar3D{
public:
//	SHORTMATRIX  *type;
	GeoMatrix<short> type;
//	LONGMATRIX	 *lnum;
	GeoMatrix<long> lnum;
//	DOUBLETENSOR *Dzl;
	GeoTensor<double> Dzl;
//	DOUBLETENSOR *w_liq;
	GeoTensor<double> w_liq;
//	DOUBLETENSOR *w_ice;
	GeoTensor<double> w_ice;
//	DOUBLETENSOR *T;
	GeoTensor<double> T;
	};
//	} STATEVAR_3D;


//typedef struct {
class Statevar1D{
public:
	short type;
	long lnum;
//	DOUBLEVECTOR *Dzl;
	GeoVector<double> Dzl;
//	DOUBLEVECTOR *w_liq;
	GeoVector<double>  w_liq;
//	DOUBLEVECTOR *w_ice;
	GeoVector<double> w_ice;
//	DOUBLEVECTOR *T;
	GeoVector<double> T;
    };
//	} STATEVAR_1D;

//typedef struct {
class Snow{
	 public:
//	STATEVAR_3D *S;
	Statevar3D *S;
//	STATEVAR_1D *S_for_BS;
	Statevar1D *S_for_BS;
//	DOUBLEVECTOR *age;
	GeoVector<double> age;
//	DOUBLEVECTOR *MELTED;
	GeoVector<double> MELTED;
//	DOUBLEVECTOR *melted;
	GeoVector<double> melted;
//	DOUBLEVECTOR *SUBL;
	GeoVector<double> SUBL;
//	DOUBLEVECTOR *subl;
	GeoVector<double> subl;
//	DOUBLEVECTOR *t_snow;
	GeoVector<double> t_snow;
//	SHORTVECTOR *yes;
	GeoVector<short> yes;
//	DOUBLEMATRIX *Qsub;
	GeoMatrix<double> Qsub;
//	DOUBLEMATRIX *Qsub_x;
	GeoMatrix<double> Qsub_x;
//	DOUBLEMATRIX *Qsub_y;
	GeoMatrix<double> Qsub_y;
//	DOUBLEMATRIX *Nabla2_Qtrans;
	GeoMatrix<double> Nabla2_Qtrans;
//	DOUBLEMATRIX *Qtrans;
	GeoMatrix<double> Qtrans;
//	DOUBLEMATRIX *Qsalt;
	GeoMatrix<double> Qsalt;
//	DOUBLEMATRIX *Qtrans_x;
	GeoMatrix<double> Qtrans_x;
//	DOUBLEMATRIX *Qtrans_y;
	GeoMatrix<double> Qtrans_y;
//	DOUBLEMATRIX *Wsubl_plot;
	GeoMatrix<double> Wsubl_plot;
//	DOUBLEMATRIX *Wtrans_plot;
	GeoMatrix<double> Wtrans_plot;
//	DOUBLEVECTOR *Dplot;
	GeoVector<double> Dplot;
//	LONGVECTOR *change_dir_wind;
	GeoVector<long> change_dir_wind;
	};
//	} SNOW;

//typedef struct {
class Glacier{
 public:
//	STATEVAR_3D *G;
	Statevar3D *G;
//	DOUBLEVECTOR *MELTED;
	GeoVector<double> MELTED;
//	DOUBLEVECTOR *melted;
	GeoVector<double> melted;
//	DOUBLEVECTOR *SUBL;
	GeoVector<double> SUBL;
//	DOUBLEVECTOR *subl;
	GeoVector<double> subl;
};
// } GLACIER;

//typedef struct{
class MeteoStations {
 public:
//	DOUBLEVECTOR *E;
	GeoVector<double> E;
//	DOUBLEVECTOR *N;
	GeoVector<double> N;
//	DOUBLEVECTOR *lat;
	GeoVector<double> lat;
//	DOUBLEVECTOR *lon;
	GeoVector<double> lon;
//	DOUBLEVECTOR *Z;
	GeoVector<double> Z;

//	DOUBLEVECTOR *sky;
	GeoVector<double> sky;
//	DOUBLEVECTOR *ST;
	GeoVector<double> ST;
//	DOUBLEVECTOR *Vheight;
	GeoVector<double> Vheight;
//	DOUBLEVECTOR *Theight;
	GeoVector<double> Theight;
//	DOUBLEVECTOR *tau_cloud_av_meteoST;// vector containing the tau_cloud_av at each meteo stations measuring SW radiation
	GeoVector<double> tau_cloud_av_meteoST; // vector containing the tau_cloud_av at each meteo stations measuring SW radiation
//	DOUBLEVECTOR *tau_cloud_meteoST;// vector containing the tau_cloud at each meteo stations measuring SW radiation
	GeoVector<double> tau_cloud_meteoST;// vector containing the tau_cloud at each meteo stations measuring SW radiation
//	SHORTVECTOR *tau_cloud_av_yes_meteoST;// flag indicating whether the tau_cloud_av at each meteo stations is available
	GeoVector<short> tau_cloud_av_yes_meteoST;// flag indicating whether the tau_cloud_av at each meteo stations is available
//	SHORTVECTOR *tau_cloud_yes_meteoST;// flag indicating whether the tau_cloud at each meteo stations is available
	GeoVector<short> tau_cloud_yes_meteoST;// flag indicating whether the tau_cloud at each meteo stations is available
//	SHORTVECTOR *flag_SW_meteoST;// flag vector saying whether a meteo station accounts for SW radiation (0: no SW, 1: SW available)
	GeoVector<short> flag_SW_meteoST;// flag vector saying whether a meteo station accounts for SW radiation (0: no SW, 1: SW available)
};
//	}METEO_STATIONS;


//typedef struct {
class Meteo {
	public:
//	METEO_STATIONS *st;
	MeteoStations *st;
	
	double ***data;
	long *numlines;
	double ***horizon;
	long *horizonlines;

//#ifdef USE_INTERNAL_METEODISTR
	double **var;
	long *line_interp_WEB;
	long *line_interp_Bsnow;
	long line_interp_WEB_LR;
	long line_interp_Bsnow_LR;
//#endif
	double **LRs;	//matrix read from the external value
	long LRsnr;		//number of lines of the matrix
	double *LRv;	//vector of interpolated values
	double **LRc;	//cyclic values from the parameter file (one vector for each LR variable)
	long *LRcnc;	//number of components of the vector (for each component)
	double *LRd;	//vector of default values
	
	double **qins;
	double *qinv;
	long qinsnr;
	long qinline;
	
	double tau_cloud;// tau_cloud for the chosen meteo station used to derive cloud
	double tau_cloud_av;// tau_cloud for the chosen meteo station used to derive cloud
	short tau_cloud_yes;
	short tau_cloud_av_yes;
//	DOUBLEMATRIX* tau_cl_map;// matrix containing the tau_cloud for each grid point
	GeoMatrix<double> tau_cl_map;// matrix containing the tau_cloud for each grid point
//	DOUBLEMATRIX* tau_cl_av_map;// matrix containing the tau_cloud_av for each grid point
	GeoMatrix<double> tau_cl_av_map;// matrix containing the tau_cloud_av for each grid point
//	SHORTMATRIX* tau_cl_map_yes;// boolean matrix saying whether the grid point has tau_cl value
	GeoMatrix<short> tau_cl_map_yes;// boolean matrix saying whether the grid point has tau_cl value
//	SHORTMATRIX* tau_cl_av_map_yes;// boolean matrix saying whether the grid point has tau_cl_av value
	GeoMatrix<short> tau_cl_av_map_yes;// boolean matrix saying whether the grid point has tau_cl_av value

//	DOUBLEMATRIX *Tgrid;
	GeoMatrix<double> Tgrid;
//	DOUBLEMATRIX *Pgrid;
	GeoMatrix<double> Pgrid;
//	DOUBLEMATRIX *Vgrid;
	GeoMatrix<double> Vgrid;
//	DOUBLEMATRIX *Vdir;
	GeoMatrix<double> Vdir;
//	DOUBLEMATRIX *RHgrid;
	GeoMatrix<double> RHgrid;
	
//	DOUBLEVECTOR *Tamean;
	GeoVector<double> Tamean;
//	DOUBLEVECTOR *Vspdmean;
	GeoVector<double> Vspdmean;
//	DOUBLEVECTOR *Vdirmean;
	GeoVector<double> Vdirmean;
//	DOUBLEVECTOR *RHmean;
	GeoVector<double> RHmean;
	
//	DOUBLEVECTOR *Taplot;
	GeoVector<double> Taplot;
//	DOUBLEVECTOR *Vxplot;
	GeoVector<double> Vxplot;
//	DOUBLEVECTOR *Vyplot;
	GeoVector<double> Vyplot;
//	DOUBLEVECTOR *RHplot;
	GeoVector<double> RHplot;

	double V;
		
//	DOUBLEMATRIX *Tday;
//	DOUBLEMATRIX *Tvar;
		

	long nstcloud;		// meteo station ID (1...n) to use for the cloudiness
	long numstcloud;	// number of meteo stations measuring cloudiness
	long nstsrad;
	long nstlrad;
	long nstTs;
	
//	LONGVECTOR *imeteo_stations;
	GeoVector<long> imeteo_stations;
	};
//	} METEO;


#ifdef USE_NETCDF
//typedef struct {
class OutputNCData {
	public:
//	DOUBLEMATRIX*	soil_thw_cum; // cumulated version of S->th
	GeoMatrix<double>  *soil_thw_cum; // cumulated version of S->th
//	DOUBLEMATRIX*	soil_thi_cum;// cumulated version of S->thice
	GeoMatrix<double>  *soil_thi_cum;// cumulated version of S->thice
//	DOUBLEMATRIX*	soil_T_cum;// cumulated version of S->T
	GeoMatrix<double>  *soil_T_cum;// cumulated version of S->T
//	DOUBLEMATRIX*	soil_P_cum;// cumulated version of S->P
	GeoMatrix<double>  *soil_P_cum;// cumulated version of S->P
//	DOUBLEMATRIX*	soil_Ptot_cum;// cumulated version of S->Ptot
	GeoMatrix<double>  *soil_Ptot_cum;// cumulated version of S->Ptot
	// snow
	
//	DOUBLETENSOR* snowD_cum; // cumulated version of snow->S->Dzl
	GeoTensor<double> *snowD_cum; // cumulated version of snow->S->Dzl
//	DOUBLETENSOR* snowT_cum; // cumulated version of snow->S->T
	GeoTensor<double> *snowT_cum; // cumulated version of snow->S->T
//	DOUBLETENSOR* snowI_cum; // cumulated version of snow->S->w_ice
	GeoTensor<double> *snowI_cum; // cumulated version of snow->S->w_ice
//	DOUBLETENSOR* snowW_cum; // cumulated version of snow->S->w_liq
	GeoTensor<double> *snowW_cum; // cumulated version of snow->S->w_liq
	// glacier
//	DOUBLETENSOR* glacD_cum; // cumulated version of glac->G->Dzl
	GeoTensor<double> *glacD_cum; // cumulated version of glac->G->Dzl
//	DOUBLETENSOR* glacT_cum; // cumulated version of glac->G->T
	GeoTensor<double> *glacT_cum; // cumulated version of glac->G->T
//	DOUBLETENSOR* glacI_cum; // cumulated version of glac->G->w_ice
	GeoTensor<double> *glacI_cum; // cumulated version of glac->G->w_ice
//	DOUBLETENSOR* glacW_cum; // cumulated version of glac->G->w_liq
	GeoTensor<double> *glacW_cum; // cumulated version of glac->G->w_liq
	// surface variables
//	DOUBLEMATRIX* egy_Rn_mean_cum;
//	DOUBLEMATRIX* egy_LWin_mean_cum;
//	DOUBLEMATRIX* egy_LW_mean_cum;
//	DOUBLEMATRIX* egy_SW_mean_cum;
//	DOUBLEMATRIX* egy_Rswdown_mean_cum;
//	DOUBLEMATRIX* egy_Rswbeam_mean_cum;
//	DOUBLEMATRIX* egy_nDt_shadow_cum;
//	DOUBLEMATRIX* egy_nDt_sun_cum;
//	DOUBLEMATRIX* egy_Rn_max_cum;
//	DOUBLEMATRIX* egy_Rn_min_cum;
//	DOUBLEMATRIX* egy_LW_max_cum;
//	DOUBLEMATRIX* egy_LW_min_cum;
//	DOUBLEMATRIX* egy_SW_max_cum;
//	DOUBLEMATRIX* egy_Rswdown_max_cum;
//	DOUBLEMATRIX* egy_SEB_mean_cum;
//	DOUBLEMATRIX* egy_G_max_cum;
//	DOUBLEMATRIX* egy_G_min_cum;
//	DOUBLEMATRIX* egy_H_mean_cum;
//	DOUBLEMATRIX* egy_H_max_cum;
//	DOUBLEMATRIX* egy_H_min_cum;
//	DOUBLEMATRIX* egy_ET_mean_cum;
//	DOUBLEMATRIX* egy_ET_max_cum;
//	DOUBLEMATRIX* egy_ET_min_cum;
//	DOUBLEMATRIX* egy_Ts_mean_cum;
//	DOUBLEMATRIX* egy_Ts_max_cum;
//	DOUBLEMATRIX* egy_Ts_min_cum;
};// OUTPUT_NCDATA;
#endif

//	typedef struct {
class AllData {
		public:
//	SOIL *S;
	Soil *S;
//	WATER *W;
	Water *W;
//	LAND *L;
	Land *L;
//	PAR *P;
	Par *P;
//	TOPO *T;
	Topo *T;
//	CHANNEL *C;
	Channel *C;
	Energy *E;
//	SNOW *N;
	Snow *N;
//	GLACIER *G;
	Glacier *G;
//	METEO *M;
	Meteo *M;
//	TIMES *I;
	Times *I;
#ifdef USE_NETCDF
	int ncid; // pointer to netCDF archive file
	long counter_snow; // counter for time print of snow
	long counter_surface_energy; //counter for surface energy maps
	long counter_soil;//counter for printing soil properties
	long counter_glac;//counter for printing glacier properties
	long counter_point; // counter for printing point data
	int z_point_var_type;// variable type for z_point (e.g. temperature in depth in a point)
	int point_var_type;// variable type (e.g. radiation in a point)
	int unstruct_z_point_var_type;// variable type for z_point (e.g. temperature in depth in a point) in unstructured grid
	int unstruct_point_var_type; // variable type (e.g. radiation in a point) in unstructured grid
//	OUTPUT_NCDATA *outnc;
	OutputNCData *outnc;
#endif
};
//	}ALLDATA;


#endif
