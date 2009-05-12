
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/



 #include "turtle.h"

/*---------------------------------------------------------------------------*/
typedef struct {

    DOUBLEMATRIX *Rn_mean;/*Map of averaged Net Radiation (on nDt_output_basin Dt time intervals). Created if(par->output_Rn>0) */
	DOUBLEMATRIX *Rn_max;/*Map of Max Net Radiation (on nDt_output_basin Dt time intervals) Created if (par->output_Rn>0 && par->distr_stat==1) */
	DOUBLEMATRIX *Rn_min;/*Map of Min Net Radiation (on nDt_output_basin Dt time intervals) Created if (par->output_Rn>0 && par->distr_stat==1) */
	DOUBLEMATRIX *LW_in;/*Map of Incoming LW Radiation (on nDt_output_basin Dt time intervals). Created if(par->output_Rn>0) */
	DOUBLEMATRIX *LW_out;/*Map of Outgoing LW Radiation (on nDt_output_basin Dt time intervals). Created if(par->output_Rn>0) */
	DOUBLEMATRIX *LW_max;/*Map of Max LW Radiation (on nDt_output_basin Dt time intervals) Created if (par->output_Rn>0 && par->distr_stat==1) */
	DOUBLEMATRIX *LW_min;/*Map of Min LW Radiation (on nDt_output_basin Dt time intervals) Created if (par->output_Rn>0 && par->distr_stat==1) */
	DOUBLEMATRIX *SW;/*Map of SW Radiation (on nDt_output_basin Dt time intervals). Created if(par->output_Rn>0) */
	DOUBLEMATRIX *SW_max;/*Map of Max SW Radiation (on nDt_output_basin Dt time intervals). Created if (par->output_Rn>0 && par->distr_stat==1) */
    DOUBLEMATRIX *ET_mean;/*Map of averaged EvapoTranspiration (on nDt_output_basin Dt time intervals). Created if(par->output_ET>0) */
	DOUBLEMATRIX *ET_max;/*Map of Max EvapoTranspiration (on nDt_output_basin Dt time intervals). Created if (par->output_ET>0 && par->distr_stat==1) */
	DOUBLEMATRIX *ET_min;/*Map of Min EvapoTranspiration (on nDt_output_basin Dt time intervals). Created if (par->output_ET>0 && par->distr_stat==1) */
    DOUBLEMATRIX *H_mean;/*Map of averaged Sensible Heat flux (on nDt_output_basin Dt time intervals). Created if(par->output_H>0) */
	DOUBLEMATRIX *H_max;/*Map of Max Sensible Heat flux (on nDt_output_basin Dt time intervals). Created if(par->output_H>0 && par->distr_stat==1) */
	DOUBLEMATRIX *H_min;/*Map of Min Sensible Heat flux (on nDt_output_basin Dt time intervals). Created if(par->output_H>0 && par->distr_stat==1) */
    DOUBLEMATRIX *G_mean;/*Map of Averaged Ground Heat flux (on nDt_output_basin Dt time intervals). Created if(par->output_G>0) */
	DOUBLEMATRIX *G_max;/*Map of Max Ground Heat flux (on nDt_output_basin Dt time intervals). Created if(par->output_G>0 && par->distr_stat==1) */
	DOUBLEMATRIX *G_min;/*Map of Min Ground Heat flux (on nDt_output_basin Dt time intervals). Created if(par->output_G>0 && par->distr_stat==1) */
	DOUBLEMATRIX *G_snowsoil;/*Map of Heat flux between snow and soil (on nDt_output_basin Dt time intervals). Created if(par->output_G>0) */
    DOUBLEMATRIX *Ts_mean;/*Map of averaged Surface Temperature (on nDt_output_basin Dt time intervals) Created if(par->output_Ts>0) */
	DOUBLEMATRIX *Ts_max;/*Map of Max Surface Temperature (on nDt_output_basin Dt time intervals). Created if(par->output_Ts>0 && par->distr_stat==1) */
	DOUBLEMATRIX *Ts_min;/*Map of Min Surface Temperature (on nDt_output_basin Dt time intervals). Created if(par->output_Ts>0 && par->distr_stat==1) */
	DOUBLEMATRIX *Rswdown_mean;/*Map of averaged SWin Radiation (on nDt_output_basin Dt time intervals) Created if(par->output_Rswdown>0) */
	DOUBLEMATRIX *Rswdown_max;/*Map of Max SWin Radiation (on nDt_output_basin Dt time intervals) Created if(par->output_Rswdown>0 && par->distr_stat==1) */
	DOUBLEMATRIX *Ta_mean;/*Map of averaged Air temperature (on nDt_output_basin Dt time intervals) Created if(par->output_meteo>0) */
	DOUBLEMATRIX *Ta_max;/*Map of Max Air temperature (on nDt_output_basin Dt time intervals) Created if(par->output_meteo>0 && par->distr_stat==1) */
	DOUBLEMATRIX *Ta_min;/*Map of Min Air temperature (on nDt_output_basin Dt time intervals) Created if(par->output_meteo>0 && par->distr_stat==1) */
	DOUBLEMATRIX *out1;/*matrix with pixel's values for the output*/
	DOUBLEVECTOR *out2;/*vector with basin's averages for the output*/
	DOUBLEMATRIX *out3;/* matrix with altimetric stripes values for the output*/
	DOUBLEMATRIX *Rswbeam;/*Map of averaged SWbeam Radiation (on nDt_output_basin Dt time intervals) Created if(par->output_Rswdown>0) */
	LONGMATRIX *nDt_shadow;/* Map of number of Dt in which the pixel is in shadow */
	LONGMATRIX *nDt_sun;/* Map of number of Dt in which the pixel is in sun. Created if(par->ES_num>0) */

	DOUBLEMATRIX *Hgplot;/* sensible heat flux: activated if par->JD_plots has more than one component */
	DOUBLEMATRIX *LEgplot;/* EvapoTranspiration flux: activated if par->JD_plots has more than one component */
	DOUBLEMATRIX *Hvplot;
	DOUBLEMATRIX *LEvplot;

	DOUBLEMATRIX *SWinplot;/* SWin flux: activated if par->JD_plots has more than one component */
	DOUBLEMATRIX *SWgplot;
	DOUBLEMATRIX *SWvplot;

	DOUBLEMATRIX *LWinplot;
	DOUBLEMATRIX *LWgplot;
	DOUBLEMATRIX *LWvplot;

	DOUBLEMATRIX *Tgplot;
	DOUBLEMATRIX *Tvplot;
	DOUBLEMATRIX *Tsplot;/* Surface Temperature flux: activated if par->JD_plots has more than one component */

	DOUBLEMATRIX *SWin;
	DOUBLEMATRIX *LWin;

	DOUBLEMATRIX *Hgrid;/* map of sensible heat flux [W/m2] */
	DOUBLEMATRIX *Tsgrid;/* map of surface temperature [ûC] */

	DOUBLEVECTOR *VSFA;
	DOUBLEVECTOR *HSFA;

	double hsun;/* solar elevation angle (radiants) */
	double dsun;/* solar azimuth angle (from N clockwise, radiants) */

} ENERGY;


/*---------------------------------------------------------------------------*/
typedef struct {

	SHORTMATRIX *type;/* matrix of soil type */
	DOUBLETENSOR *pa; /* doubletensor of soil parameters */
	DOUBLETENSOR *P; /* soil pressure */
	DOUBLETENSOR *T; /* soil temperature */
	DOUBLETENSOR *thice; /* soil theta_ice */
	DOUBLEMATRIX *Jinf; /* water infiltration in the soil surface */
	DOUBLETENSOR *J; /*  water flux outgoing from one layer to the lower one */
	DOUBLEMATRIX *Tv;
	DOUBLETENSOR *Tav;
	DOUBLETENSOR *thwav;
	DOUBLETENSOR *thiav;

} SOIL;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *Z0;         /*elevetions of each pixel (DEM)*/
	DOUBLEMATRIX *Z0dp;		  /*DEM depitted*/
	DOUBLEMATRIX *Z0ext;      /*DEM extended (to avoid curvature problems)*/
    DOUBLEMATRIX *sky;        /*view factor (of the sky) for each pixel*/
    SHORTMATRIX *pixel_type; /*0=land,9=novalue,10=channel,11=lake,12=sea*/
    SHORTMATRIX *DD;         /*Drainage Directions for each pixel; ex matr_ev->slopes*/
    DOUBLEMATRIX *i_DD;       /*slope along Drainage Direction for each pixel*/
    DOUBLEMATRIX *dz_dx;      /*slope along x direction for each pixel*/
    DOUBLEMATRIX *dz_dy;      /*slope along y direction for each pixel*/
    DOUBLEMATRIX *top_index; /*topographic index for water content distribution*/
    SHORTMATRIX *curv;       /*topology of curvature (1 and 0) ex matr_ev->curv*/
    DOUBLEMATRIX *area;       /*area of a pixel considering the slopeex matr_ev->area*/
    DOUBLEMATRIX *aspect;     /*aspect; ex: matr_ev->azimuth*/
    DOUBLEMATRIX *slopes;     /*slope of the pixels; ex: matr_ev->slopes*/
    DOUBLEMATRIX *i_ch;       /*slope to calculate the surface velocity of the channel incoming flow*/
	DOUBLEMATRIX *pixel_distance;
	double Zmin;
	double Zmax;
	LONGVECTOR *ES_pixel;
	DOUBLEVECTOR *ES_aspect;
	DOUBLEVECTOR *ES_slope;
	double ****horizon_height;
} TOPO;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *LC;            //land cover (wood,lake,town,...) for each pixel*/
	SHORTMATRIX *LC2;			  //FURTHER LAND CLASSIFICATION
    DOUBLEMATRIX *albedo;         //albedo calculated for each pixel*/
	SHORTMATRIX *shadow;		  //=1 if shadow, =0 if not*/
	LONGVECTOR *clax;
	LONGMATRIX *cont;
	DOUBLEMATRIX *ty; /* land type (better would be land cover, i.e. pasture, rock, forest, paved_road...) */
	DOUBLEVECTOR *LAI; /* leaf area index */

} LAND;/*all this data are calculated on the basis of land use data and some other par*/


/*---------------------------------------------------------------------------*/
typedef struct {/*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                  R=number of rows of the basin,C=number of columns in the basin*/
    LONGVECTOR *r;          /*array of rows of the channel-pixels; dimension=nch*/
    LONGVECTOR *c;          /*array of columns of the channel-pixels; dimension=nch*/
	DOUBLEVECTOR *Q;
    DOUBLEVECTOR *s0;       /*distance of each channel-pixel from outlet; dimension=nch*/
    DOUBLEMATRIX *fraction_spread;/*fraction of flow for each s(=virtual distance) of all channel-pixel;
                                    dimension=nch*ns */
    /*cnet-flow in virtual stretches of channel:*/
    DOUBLEVECTOR *Q_sup_s;   /*derived from q_sup[mc/s]; dimension=ns*/
    DOUBLEVECTOR *Q_sub_s;   /*derived from q_sub[mc/s]; dimension=ns*/
	DOUBLEVECTOR *Q_sup_spread;
    DOUBLEVECTOR *Q_sub_spread;
} CHANNEL;




/*---------------------------------------------------------------------------*/
typedef struct { /*nstations=number of all the rain-stations,npixel=number of all the pixels of the basin R*C,
                   R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
    DOUBLEMATRIX *weights_Kriging; /*dimension=npixel*nstations */
    DOUBLETENSOR *q_sub;            /*ground flow (subsurface) in [mm/s]; positive if go out from the cell; dimension=(number of layers)*R*C*/
	DOUBLEMATRIX *q_sup;/*runoff discharge in [mm/s] */
    DOUBLEMATRIX *h_sup;    /*height of water over the sl-surface not infiltrated in mm; dimension=R*C*/
    DOUBLEMATRIX *total;    /*total(snow+rain) precipitation in mm (in a Dt)*/
    DOUBLEMATRIX *Pn;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                              same subroutine and in "water.balance.c" module*/
    DOUBLEMATRIX *wt;       /*intercepted precipitation in mm*/
    DOUBLEMATRIX *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
    DOUBLEMATRIX *PrSNW_mean;/* average of precipitation fallen as snow */
	DOUBLEMATRIX *Psnow;/* precipitation calculated as snow (as result of the air temperature) */

	//output variables
    DOUBLEMATRIX *out1;/*matrix with pixel's values for the output*/
	DOUBLEVECTOR *out2;/*matrix with basin's values for the output*/

	DOUBLEMATRIX *hsupav;/* average of water height on a pixel in the output time step [mm]. Created if(par->output_h_sup>0) */

	DOUBLEMATRIX *outfluxes;

} WATER;


/*---------------------------------------------------------------------------*/
typedef struct {
    short iter;    /*current iteration of a time interval for egy-mass balance*/
    short n_iter;  /*n_iter=number of iterations to do for each time-step*/
    double TH;     /*TH=last of all the simulation in hours*/
    long i_pixel;  /*counter for the output of a pixel*/
    long n_pixel;  /*nDt_output_pixel=number of Dt after which the output of a pixel are printed*/
	long i_basin;
	long n_basin;
	long i_plot;
	long n_plot;/* hour interval after which the special output are printed */
	long nt_plot;
	long d_plot;
    double JD;      /*day=current Julian day during the simulation*/
    long DD;       /*current day of the month during the simulation*/
    long MM;       /*current month of the year during the simulation*/
    long AAAA;     /*current year*/
    long hh;       /*current hour of day during the simulation*/
    long mm;       /*current minute of hour during the simulation*/
    double Dt;      /*Dt=the integration time interval [s]*/
    double time;    /*time=current time from the begin of simulation [s]*/
    double egy;  /*the time of egy subroutine [s]*/
    double vert_wb; /*the time of water-vertical-distribution subroutine [s]*/
    double horiz_wb;/*the time of water-horizontal-distribution subroutine [s]*/
    double writeout;/*the time of write-output subroutine [s]*/
} TIMES;


/*---------------------------------------------------------------------------*/
typedef struct {
    double Dt; /*Dt=the integration time interval [s]*/
	double JD0;/* Decimal julian day of the beginning of simulation (0.0 - 365.99)*/
	long year0;/* Year of the beginning of the simulation*/
	double ST;/* Standard time to which all the output data are referred (difference respect UMT, in hour) */
    short print; /*1 IF YOU WANT TO PRINT MATRICES WITH INTERMEDIATE RESULTS, 0 OTHERWISE*/
    short monin_obukhov;
    double gamma_m;/*Exponent of the law of uniform motion on the surface*/
    double T_rain; /*TEMPERATURE ABOVE WICH ALL PRECIPITAION IS RAIN [C]*/
    double T_snow; /*TEMPERATURE BELOW WICH ALL PRECIPITAION IS SNOW [C]*/
    double aep; /*ALBEDO EXTINCTION PARAMETER [m]*/
    double avo; /*NEW SNOW VISIBLE BAND REFLECTANCE*/
    double airo; /*NEW NEAR INFRARED BAND REFLECTANCE*/
    double Sr; /*WATER FRACTION RETAINED BY CAPILLARY FORCES IN SNOW*/
    double rho_ice; /*Ice density [kg/mc]*/
    long total_pixel;/*The number of the valid pixel of the whole basin*/
    long n_error; /*Current number of error of the simulation*/
    long max_error; /*Maximum number of error for the simulation*/
	long snowlayer_max;/* max number of snow layers*/
	long snowlayer_inf;/* layer of unlimited thickness (beginning from below) */
	DOUBLEVECTOR *Dmin;/* vector of Min snow-layer thickness [mm] */
	DOUBLEVECTOR *Dmax;/* vector of Max snow-layer thickness [mm] */
	double Sr_glac;
	long glaclayer_max;/* maximum number of glacier layers, it must be at least 1 (if it is 0 the glacier module is switched off)*/
	DOUBLEVECTOR *Dmin_glac;/* Min thickness of glacier layers [mm] */
	DOUBLEVECTOR *Dmax_glac;/* Max thickness of glacier layers [mm] */

    short state_snow;/* 1 if you have a file with initial snow depth (in days), 0 otherwise */
	short state_glac;/* 1 if you have a file with glacier depth, 0 otherwise */
	short state_turb;
	short state_snow_age;/* 1 if you have a file with snow age (in days), 0 otherwise */
	short state_lwrad;

	double imp;
	double f_bound_Richards;

	double epsilon_snow;

	double output_Txy; /* 1 if you want to display output Txy distributed MAPS, 0 otherwise */
	double output_TETAxy;/* 1 if you want to display output TETAxy distributed MAPS, 0 otherwise */
	double output_TETAICExy; /* 1 if you want to display output TETA_ICExy distributed MAPS, 0 otherwise */
	double output_PSIxy; /* 1 if you want to display output PSIxy distributed MAPS, 0 otherwise */
	double output_snow; /* 1 if you want to display output SNOW distributed MAPS, 0 otherwise */
	double output_glac;/* 1 if you want to display output GLACIER distributed MAPS, 0 otherwise */
	double output_h_sup;/* 1 if you want to display output h_sup distributed MAPS, 0 otherwise */
	double output_albedo;/* 1 if you want to display output albedo distributed MAPS, 0 otherwise */
    double output_Rn;/* 1 if you want to display output radiation distributed MAPS, 0 otherwise */
	double output_G;/* 1 if you want to display output G distributed MAPS, 0 otherwise */
	double output_H;/* 1 if you want to display output H distributed MAPS, 0 otherwise */
	double output_ET;/* 1 if you want to display output ET distributed MAPS, 0 otherwise */
	double output_Ts;/* 1 if you want to display output Ts distributed MAPS, 0 otherwise */
	double output_P;/* 1 if you want to display output Precipitation distributed MAPS, 0 otherwise */
	double output_Wr;/* 1 if you want to display output Wr (water stored in vegetation) distributed MAPS, 0 otherwise */
	double output_balancesn;/* 1 if you want to display output snow melting and sublimation distributed MAPS, 0 otherwise */
	double output_balancegl;/* 1 if you want to display output glacier melting and sublimation distributed MAPS, 0 otherwise */
	double output_Rswdown;/* 1 if you want to display output SWin MAPS, 0 otherwise */
	double output_meteo;/* 1 if you want to display output meteo MAPS, 0 otherwise */

	DOUBLEMATRIX *chkpt;/* Matrix (n X 3) where n is the number of points that have to be plotted and 3 is the number of parameters: if state_px_coord==1 (E, N, layer) ---  if state_px_coord==0 (row,col,layer) (max 9999 points)*/
	LONGMATRIX *rc;/* Matrix (n X 2) where n is the number of points that have to be plotted. The first represents the row of the point in the raster, the second the column */
	short ES_num; /* see __control file: >=1 (or<=1) how many altimetric stripes you want to consider (up to 99) (negative value means that this is done only for glacier pixels) */
	short state_px_coord;/* 1 if all coordinates are in (East-North) format, 0 if in (row, columns) format */

	double integr_scale_rain;/* range of the variogram */
	double variance_rain;/* sill of the variogram */

	short recover;/* =1 if you want to recover a simulation, 0 otherwise */

	double Vmin;/* minimum wind speed [m/s]*/

	double snowcorrfact;/* correction factor for snow precipitation */
	double raincorrfact;/* correction factor for rain precipitation */

	double RHmin; /* minimum relative humidity of the air [%] */

	short format_out;
	short nsky;
	double channel_thres;

	DOUBLEVECTOR *saving_points;

	double Vis; //visibility in km (>5 km)
	double Lozone; //thickness of the stratospheric ozone layer (in cm normal conditions)

	short point_sim;/* =0 distributed simulation, =1 point simulation (the parameter files are different in the two cases) */
	double snow_maxpor;/* Max allowd snow porosity [-]*/
	double snow_density_cutoff;/* Snow density cutoff [Kg/m3] to change snow deformation rate */
	double drysnowdef_rate;/* Snow compaction (% per hour) due to destructive metamorphism for snow density<snow_density_cutoff and dry snow */
	double wetsnowdef_rate;/* Enhancement factor in presence of wet snow */
	double snow_viscosity;/* snow viscosity coefficient [Kg s /m2] at T=0ûC and snow_density=0 */

	double latitude;
	double longitude;

	double z0_snow;
	long n_landuses;/* the maximum digit in the first row of block 1 in parameters file */

	LONGVECTOR *JD_plots; /* vector of Julian Days in which energy balance and meteo data are plotted with a very short time step
	see block 5 _options file */
	short micromet1;/* Use Micromet for wind,T,RH,Prec,Pres (=1), otherwise (=0)*/
	short micromet2;/* Use Micromet for SWin (=1), otherwise (=0)*/
	short micromet3;/* Use Micromet for LWin (=1), otherwise (=0)*/
	short snowtrans;/* Use SnowTrans3D (=1), otherwise (=0)*/

	LONGVECTOR *r_points;
	LONGVECTOR *c_points;

	double psimin;/* Absolute minimum admitted suction potential */
	double Esoil;/* Soil comprimibility if oversaturated water is added*/
	double stmin;

	double glac_thr;
	short state_pixel;

	short wat_balance;/* 1 calculates water_balance, 0 do not calculate */
	short en_balance;/* 1 calculates energy_balance, 0 do not calculate */

	long nDt_water;
	long MaxiterVWB;
	double TolVWb;
	double MaxerrVWB;
	double dtminVWB;

	short distr_stat;/* 1 to plot distribute statistics (min, max, var), 0 otherwise */

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

} PAR;



/*---------------------------------------------------------------------------*/
typedef struct {
	SHORTMATRIX  *type;/* snow type (0,1 or 2) depending on snow height */
	LONGMATRIX	 *lnum;/* number of snow layer */
	DOUBLETENSOR *Dzl;/* snow layer depth [mm]*/
	DOUBLETENSOR *w_liq;/* mass of liquid water present in the snow per unit of surface [kg/m2] */
	DOUBLETENSOR *w_ice;/* equivalent mass of liquid water stored in the snow per unit of surface [kg/m2] */
	DOUBLETENSOR *T;/* temperature of the snow layer [ûC]*/
	DOUBLEMATRIX *age;/* matrix of snow age (in days)*/
	double evap_basin;
	double subl_basin;
	double melted_basin;
	DOUBLEVECTOR *evap;
	DOUBLEVECTOR *subl;
	DOUBLEVECTOR *melted;
	DOUBLEMATRIX *max;/* Max snow depth. Created if(par->output_snow>0)*/
	DOUBLEMATRIX *average;/* average snow depth in the print output interval. Created if(par->output_snow>0)*/
	DOUBLEMATRIX *MELTED;
	DOUBLEMATRIX *SUBL;
	DOUBLEMATRIX *t_snow;
	DOUBLEMATRIX *totav_snow;

	DOUBLEMATRIX *DDF;
	DOUBLEMATRIX *DDF1;
	DOUBLEMATRIX *DDFvar;
	LONGMATRIX *DDFcont;
	DOUBLEMATRIX *DDFTmin;
	DOUBLEMATRIX *DDFmeltTL0;
	DOUBLEMATRIX *DDFmelt;
	//DOUBLETENSOR *MC;

	DOUBLEMATRIX *rho_newsnow;/* density of the new snow [Kg/m3] */
	DOUBLEMATRIX *Wsubl;
	DOUBLEMATRIX *Wsusp;
	DOUBLEMATRIX *Wsalt;
	DOUBLEMATRIX *Wsubgrid;
	DOUBLEMATRIX *Wtot_cum;
	DOUBLEMATRIX *Wsubl_cum;
	DOUBLEMATRIX *Wsusp_cum;
	DOUBLEMATRIX *Wsalt_cum;
	DOUBLEMATRIX *Wsubgrid_cum;
	DOUBLEMATRIX *out_bs;
	DOUBLEMATRIX *ListonSWE;
	DOUBLEMATRIX *softSWE;
	DOUBLEMATRIX *softSWE1;
	DOUBLEMATRIX *Dplot;/* snow depth: activated if par->JD_plots has more than one component */

	DOUBLEVECTOR *CR1;/* destructive metamorphism in snow compaction rate */
	DOUBLEVECTOR *CR2;/* overburden load in snow compaction rate*/
	DOUBLEVECTOR *CR3;
	DOUBLEVECTOR *CR1m;/* average destructive metamorphism snow compaction rate */
	DOUBLEVECTOR *CR2m;/* overburden load in snow compaction rate*/
	DOUBLEVECTOR *CR3m;

} SNOW;

typedef struct {
	LONGMATRIX	 *lnum;/* number of glacier layers */
	DOUBLETENSOR *Dzl;/* layer depth [mm] */
	DOUBLETENSOR *w_liq;/* mass of liquid water present in the ice per unit of surface [kg/m2] */
	DOUBLETENSOR *w_ice;/* equivalent mass of liquid water stored in the ice per unit of surface [kg/m2] */
	DOUBLETENSOR *T;/* temperature of the glacier layer [ûC]*/
	double evap_basin;
	double subl_basin;
	double melted_basin;
	DOUBLEVECTOR *evap;
	DOUBLEVECTOR *subl;
	DOUBLEVECTOR *melted;
	DOUBLEMATRIX *MELTED;
	DOUBLEMATRIX *SUBL;
	DOUBLEMATRIX *DDF;
	DOUBLEMATRIX *DDF1;
	DOUBLEMATRIX *DDFvar;
	LONGMATRIX *DDFcont;
	DOUBLEMATRIX *DDFTmin;
	DOUBLEMATRIX *DDFmeltTL0;
	DOUBLEMATRIX *DDFmelt;
	//DOUBLETENSOR *MC;
} GLACIER;

typedef struct{
	DOUBLEVECTOR *E;/* East coordinate [m] of the meteo station */
	DOUBLEVECTOR *N;/* North coordinate [m] of the meteo station */
	DOUBLEVECTOR *lat;/* Latitude [rad] of the meteo station */
	DOUBLEVECTOR *lon; /* Longitude [rad] of the meteo station */
	DOUBLEVECTOR *Z;/* Elevation [m] of the meteo station */
	DOUBLEVECTOR *sky;/* Sky-view-factor [-] of the meteo station */
	DOUBLEVECTOR *ST;/* Standard time minus UTM [hours] of the meteo station */
	DOUBLEVECTOR *Vheight;/* Wind velocity measurement height [m] (a.g.l.)  */
	DOUBLEVECTOR *Theight;/* Air temperature measurement height [m] (a.g.l.)  */
	DOUBLEVECTOR *JD0;/* Decimal Julian Day of the first data */
	LONGVECTOR *Y0;/* Year of the first data */
	DOUBLEVECTOR *Dt;/* Dt of sampling of the data [sec]*/
	LONGVECTOR *offset;/* offset column */
} METEO_STATIONS;


typedef struct {
	METEO_STATIONS *st;
	double ***data;/* tensor (#meteo_stat X #rows_meteo_data X #actual_meteo_variables) of ALL the data variables for all the stations */
	long **column;/* matrix (#meteo_stat X #allowed_meteo_variables) indicating in which column lies the correspondent meteo variable */
	double ***horizon;
	double **var;/* matrix (#meteo_stat X #allowed_meteo_variables) of the interpolated value of the meteo data for the current time */

	DOUBLEMATRIX *Tgrid;/* Map of air temperature in each point */
	DOUBLEMATRIX *Pgrid;/* Map of air pressure in each point */
	DOUBLEMATRIX *Vgrid;/* Map of wind velocity in each point */
	DOUBLEMATRIX *Vdir;/* Map of wind direction in each point */
	DOUBLEMATRIX *RHgrid;/* Map of relative humidity in each point */
	DOUBLEMATRIX *CFgrid;

	DOUBLEMATRIX *Vspdmean;
	DOUBLEMATRIX *Vdirmean;
	DOUBLEMATRIX *RHmean;

	DOUBLEMATRIX *Taplot;
	DOUBLEMATRIX *Vspdplot;
	DOUBLEMATRIX *Vdirplot;
	DOUBLEMATRIX *RHplot;

	double V;/* wind velocity [m/s] */
	double RH;/* relative humidity [%]*/
	double LapseRate;

	DOUBLEMATRIX *Tday;
	DOUBLEMATRIX *Tvar;

	//Liston
	float *LT;
	float *Lrh;
	float *Lws;
	float *Lwd;
	float *LP;

	long nstsrad; /* number of the station from which one gets the data of SW radiation */
	long nstlrad; /* number of the station from which one gets the data of LW radiation */
	long nstcloud; /* number of the station from which one gets the data of cloudiness */

} METEO;


