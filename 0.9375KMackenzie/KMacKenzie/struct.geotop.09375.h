
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

    DOUBLEMATRIX *Rn_mean;
	DOUBLEMATRIX *Rn_max;
	DOUBLEMATRIX *Rn_min;
	DOUBLEMATRIX *LW_in;
	DOUBLEMATRIX *LW_out;
	DOUBLEMATRIX *LW_max;
	DOUBLEMATRIX *LW_min;
	DOUBLEMATRIX *SW;
	DOUBLEMATRIX *SW_max;
    DOUBLEMATRIX *ET_mean;
	DOUBLEMATRIX *ET_max;
	DOUBLEMATRIX *ET_min;
    DOUBLEMATRIX *H_mean;
	DOUBLEMATRIX *H_max;
	DOUBLEMATRIX *H_min;
    DOUBLEMATRIX *G_mean;
	DOUBLEMATRIX *G_max;
	DOUBLEMATRIX *G_min;
	DOUBLEMATRIX *G_snowsoil;
    DOUBLEMATRIX *Ts_mean;  /*averaged surface Temperature(on nDt_output_basin Dt time intervals)*/
	DOUBLEMATRIX *Ts_max;
	DOUBLEMATRIX *Ts_min;
	DOUBLEMATRIX *Rswdown_mean;
	DOUBLEMATRIX *Rswdown_max;
	DOUBLEMATRIX *Ta_mean;
	DOUBLEMATRIX *Ta_max;
	DOUBLEMATRIX *Ta_min;
	DOUBLEMATRIX *out1;			/*matrix with pixel's values for the output*/
	DOUBLEVECTOR *out2;			/*vector with basin's means for the output*/
	DOUBLEMATRIX *out3;
	DOUBLEMATRIX *Rswbeam;
	LONGMATRIX *nDt_shadow;
	LONGMATRIX *nDt_sun;

	DOUBLEMATRIX *Hplot;
	DOUBLEMATRIX *LEplot;
	DOUBLEMATRIX *SWinplot;
	DOUBLEMATRIX *SWoutplot;
	DOUBLEMATRIX *LWinplot;
	DOUBLEMATRIX *LWoutplot;
	DOUBLEMATRIX *Tsplot;

	DOUBLEMATRIX *SWin;
	DOUBLEMATRIX *LWin;

	DOUBLEMATRIX *Hgrid;
	DOUBLEMATRIX *Tsgrid;

	double VSFA;
	double HSFA;

	double hsun;
	double dsun;

} ENERGY;


/*---------------------------------------------------------------------------*/
typedef struct {

	SHORTMATRIX *type;
	DOUBLETENSOR *pa;
	DOUBLETENSOR *P;
	DOUBLETENSOR *T;
	DOUBLETENSOR *thice;
	DOUBLEMATRIX *Jinf;
	DOUBLETENSOR *J;

} SOIL;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *Z0;        /*elevetions of each pixel (DEM)*/
	DOUBLEMATRIX *Z0dp;
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
    SHORTMATRIX *use;             /*land use (wood,lake,town,...) for each pixel*/
    DOUBLEMATRIX *albedo;         /*albedo calculated for each pixel*/
	SHORTMATRIX *shadow;		  /*=1 if shadow, =0 if not*/
	LONGVECTOR *clax;
	LONGMATRIX *cont;
	DOUBLEMATRIX *ty;
	DOUBLEVECTOR *LAI;

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
    DOUBLETENSOR *q_sub;           /*subsuperficial flows in mm/s; positive if go out from the cell; dimension=(number of layers)*R*C*/
	DOUBLEMATRIX *q_sup;
    DOUBLEMATRIX *h_sup;    /*height of water over the sl-surface not infiltrated in mm; dimension=R*C*/
    DOUBLEMATRIX *total;    /*total(snow+rain) precipitation in mm (in a Dt)*/
    DOUBLEMATRIX *Pn;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                              same subroutine and in "water.balance.c" module*/
    DOUBLEMATRIX *wt;       /*intercepted precipitation in mm*/
    DOUBLEMATRIX *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
    DOUBLEMATRIX *PrSNW_mean;
	DOUBLEMATRIX *Psnow;

	//output variables
    DOUBLEMATRIX *out1;
	DOUBLEVECTOR *out2;

	DOUBLEMATRIX *hsupav;

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
	long n_plot;
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
    double Dt;      /*Dt=the integration time interval [s]*/
	double JD0;
	long year0;
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
    long n_error;        /*Current number of error of the simulation*/
    long max_error;      /*Maximum number of error for the simulation*/
	long snowlayer_max;
	long snowlayer_inf;
	DOUBLEVECTOR *Dmin;
	DOUBLEVECTOR *Dmax;
	double Sr_glac;
	long glaclayer_max;
	DOUBLEVECTOR *Dmin_glac;
	DOUBLEVECTOR *Dmax_glac;

    short state_snow;    /*0 IF YOU DO NOT HAVE A FILE WITH SNOW HIGTH, 1 OTHERWISE*/
	short state_glac;
	short state_turb;
	short state_snow_age;
	short state_lwrad;

	double imp;
	double f_bound_Richards;

	double epsilon_snow;

	double output_Txy;
	double output_TETAxy;
	double output_TETAICExy;
	double output_PSIxy;
	double output_snow;
	double output_glac;
	double output_h_sup;
	double output_albedo;
    double output_Rn;
	double output_G;
	double output_H;
	double output_ET;
	double output_Ts;
	double output_P;
	double output_Wr;
	double output_balancesn;
	double output_balancegl;
	double output_Rswdown;
	double output_meteo;

	DOUBLEMATRIX *chkpt;
	LONGMATRIX *rc;
	short ES_num;
	short state_px_coord;

	double integr_scale_rain;
	double variance_rain;

	short recover;

	double Vmin;

	double snowcorrfact;
	double raincorrfact;

	double RHmin;

	short format_out;
	short nsky;
	double channel_thres;

	DOUBLEVECTOR *saving_points;

	double Vis; //visibility in km (>5 km)
	double Lozone; //thickness of the stratospheric ozone layer (in cm normal conditions)

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

	LONGVECTOR *JD_plots;

	short micromet1;
	short micromet2;
	short micromet3;
	short snowtrans;

	LONGVECTOR *r_points;
	LONGVECTOR *c_points;

	double psimin;
	double Esoil;
	double stmin;

	double glac_thr;
	short state_pixel;

	short wat_balance;
	short en_balance;

	long nDt_water;
	long MaxiterVWB;
	double TolVWb;
	double MaxerrVWB;
	double dtminVWB;

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

} PAR;



/*---------------------------------------------------------------------------*/
typedef struct {
	SHORTMATRIX  *type;
	LONGMATRIX	 *lnum;
	DOUBLETENSOR *Dzl;
	DOUBLETENSOR *w_liq;
	DOUBLETENSOR *w_ice;
	DOUBLETENSOR *T;
	DOUBLEMATRIX *age;
	double evap_basin;
	double subl_basin;
	double melted_basin;
	DOUBLEVECTOR *evap;
	DOUBLEVECTOR *subl;
	DOUBLEVECTOR *melted;
	DOUBLEMATRIX *max;
	DOUBLEMATRIX *average;
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

	DOUBLEMATRIX *rho_newsnow;
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
	DOUBLEMATRIX *Dplot;

	DOUBLEVECTOR *CR1;
	DOUBLEVECTOR *CR2;
	DOUBLEVECTOR *CR3;
	DOUBLEVECTOR *CR1m;
	DOUBLEVECTOR *CR2m;
	DOUBLEVECTOR *CR3m;

} SNOW;

typedef struct {
	LONGMATRIX	 *lnum;
	DOUBLETENSOR *Dzl;
	DOUBLETENSOR *w_liq;
	DOUBLETENSOR *w_ice;
	DOUBLETENSOR *T;
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
	DOUBLEVECTOR *E;
	DOUBLEVECTOR *N;
	DOUBLEVECTOR *lat;
	DOUBLEVECTOR *lon;
	DOUBLEVECTOR *Z;
	DOUBLEVECTOR *sky;
	DOUBLEVECTOR *ST;
	DOUBLEVECTOR *Vheight;
	DOUBLEVECTOR *Theight;
	DOUBLEVECTOR *JD0;
	LONGVECTOR *Y0;
	DOUBLEVECTOR *Dt;
	LONGVECTOR *offset;
} METEO_STATIONS;


typedef struct {
	METEO_STATIONS *st;
	double ***data;
	long **column;
	double ***horizon;
	double **var;

	DOUBLEMATRIX *Tgrid;
	DOUBLEMATRIX *Pgrid;
	DOUBLEMATRIX *Vgrid;
	DOUBLEMATRIX *Vdir;
	DOUBLEMATRIX *RHgrid;

	DOUBLEMATRIX *Vspdmean;
	DOUBLEMATRIX *Vdirmean;
	DOUBLEMATRIX *RHmean;

	DOUBLEMATRIX *Taplot;
	DOUBLEMATRIX *Vspdplot;
	DOUBLEMATRIX *Vdirplot;
	DOUBLEMATRIX *RHplot;

	double V;
	double RH;
	double LapseRate;

	DOUBLEMATRIX *Tday;
	DOUBLEMATRIX *Tvar;

	//Liston
	float *LT;
	float *Lrh;
	float *Lws;
	float *Lwd;
	float *LP;

	long nstsrad;
	long nstlrad;
	long nstcloud;

} METEO;


