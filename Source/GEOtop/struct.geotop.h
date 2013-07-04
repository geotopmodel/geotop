
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
    
 #include "turtle.h"
 #include "t_datamanipulation.h"
 #include "t_utilities.h"
 #include "tensor3D.h"
 //#include "turtle2umfpack.h"

 
/*---------------------------------------------------------------------------*/
typedef struct {

    DOUBLEMATRIX *Rn_mean; 
	DOUBLEMATRIX *Rn_max;
	DOUBLEMATRIX *Rn_min;   
	DOUBLEMATRIX *LW_in;
	DOUBLEMATRIX *LW;
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
    DOUBLEMATRIX *SEB_mean; 
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
	
	DOUBLEMATRIX *Rswbeam;
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

	DOUBLEMATRIX *SWin;
	DOUBLEMATRIX *LWin;
	
	double hsun;
	double dsun;

	DOUBLEVECTOR *Dlay;
	DOUBLEVECTOR *wliq;
	DOUBLEVECTOR *wice;
	DOUBLEVECTOR *Temp; 
	DOUBLEVECTOR *deltaw;
	DOUBLEVECTOR *SWlayer;
	DOUBLEVECTOR *soil_transp_layer;
	DOUBLEVECTOR *dFenergy;
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
	
} ENERGY;


/*---------------------------------------------------------------------------*/
typedef struct {

	LONGMATRIX *type;
	DOUBLETENSOR *pa;
	DOUBLETENSOR *P;
	DOUBLETENSOR *Ptot;
	DOUBLETENSOR *T;
	DOUBLETENSOR *thice;
	DOUBLETENSOR *th;
	DOUBLEMATRIX *Jinf;
	DOUBLEMATRIX *Tv;
	SHORTMATRIX *bc;
	DOUBLETENSOR *ET;
	
	DOUBLEMATRIX *th_av;
	DOUBLEMATRIX *thice_av;
	DOUBLEMATRIX *T_av;

} SOIL;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *Z0;         /*elevetions of each pixel (DEM)*/
	DOUBLEMATRIX *Z1;
	DOUBLEMATRIX *Z0dp;		  /*DEM depitted*/
	DOUBLEMATRIX *Z0ext;      /*DEM extended (to avoid curvature problems)*/
	DOUBLETENSOR *Z;
    DOUBLEMATRIX *sky;        /*view factor (of the sky) for each pixel*/
    SHORTMATRIX *pixel_type; 
    
	SHORTMATRIX *DD;         /*Drainage Directions for each pixel; ex matr_ev->slopes*/
	LONGMATRIX *DDup;
	LONGVECTOR *DDdown;
	
    DOUBLEMATRIX *i_DD;       /*slope along Drainage Direction for each pixel*/
    DOUBLEMATRIX *top_index; /*topographic index for water content distribution*/
    DOUBLEMATRIX *area;       /*area of a pixel considering the slopeex matr_ev->area*/
    DOUBLEMATRIX *aspect;     /*aspect; ex: matr_ev->azimuth*/
    DOUBLEMATRIX *slopes;     /*slope of the pixels; ex: matr_ev->slopes*/
	double ****horizon_height;
	
	//for micromet
	DOUBLEMATRIX *Zm;
	DOUBLEMATRIX *curv_m;		
	DOUBLEMATRIX *slope_m;
	DOUBLEMATRIX *slopeaz_m;
	
	long ***i_cont;
	LONGMATRIX *lrc_cont;
	
	long **j_cont;
	LONGMATRIX *rc_cont;
		
	LONGVECTOR *Lp;
	LONGVECTOR *Li;
	LONGVECTOR *Up;
	LONGVECTOR *Ui;
			
} TOPO;


/*---------------------------------------------------------------------------*/
typedef struct {
    DOUBLEMATRIX *LC;            //land cover (wood,lake,town,...) for each pixel*/
	//SHORTMATRIX *LC2;			  //FURTHER LAND CLASSIFICATION 	
    DOUBLEMATRIX *albedo;         //albedo calculated for each pixel*/
	SHORTMATRIX *shadow;		  //=1 if shadow, =0 if not*/
	LONGMATRIX *cont;
	DOUBLEMATRIX *ty;
	
	LONGMATRIX *vegparp;
	double ***vegpars;
	double **vegparv;
	DOUBLEVECTOR *vegpar;
	
	DOUBLEMATRIX *root_fraction;
	
} LAND;/*all this data are calculated on the basis of land use data and some other par*/


/*---------------------------------------------------------------------------*/
typedef struct {/*nch=number of channel-pixel,ns=number of virtual stretches of channel,L=number of layers,
                  R=number of rows of the basin,C=number of columns in the basin*/
    LONGVECTOR *r;          /*array of rows of the channel-pixels; dimension=nch*/
    LONGVECTOR *c;          /*array of columns of the channel-pixels; dimension=nch*/
	LONGMATRIX *ch;
	DOUBLEVECTOR *Qsup;
	DOUBLEVECTOR *Qsub;	
    //DOUBLEVECTOR *s0;       /*distance of each channel-pixel from outlet; dimension=nch*/
    //DOUBLEMATRIX *fraction_spread;/*fraction of flow for each s(=virtual distance) of all channel-pixel*/
    //DOUBLEVECTOR *Q_sup_s;   /*derived from q_sup[mc/s]; dimension=ns*/
    //DOUBLEVECTOR *Q_sub_s;   /*derived from q_sub[mc/s]; dimension=ns*/
	//DOUBLEVECTOR *Qsup_spread;
    //DOUBLEVECTOR *Qsub_spread;
	
	DOUBLEVECTOR *h_sup;
	DOUBLEVECTOR *dh_sup;
	
	DOUBLEMATRIX *hsupav;
	
	DOUBLEVECTOR *length;
	
	double Q_out;
			
} CHANNEL;




/*---------------------------------------------------------------------------*/
typedef struct { /*nstations=number of all the rain-stations,npixel=number of all the pixels of the basin R*C,
                   R=number of rows,C=number of columns,nt=number of time-step of the whole similation*/
    DOUBLEMATRIX *weights_Kriging; /*dimension=npixel*nstations */
    DOUBLEMATRIX *PrecTot;    /*total(snow+rain) precipitation in mm (in a Dt)*/
    DOUBLEMATRIX *Pnet;       /*liquid precipitation which reaches the sl surface in mm in a Dt as input
                              of "punctual_energy" subroutine, rain intensity in mm/s as output of the
                              same subroutine and in "water.balance.c" module*/
    DOUBLEMATRIX *wcan_rain;       /*intercepted precipitation in mm*/
    DOUBLEMATRIX *wcan_snow;       /*intercepted precipitation in mm*/

    DOUBLEMATRIX *PrTOT_mean;  /*Total precipitation [mm](on nDt_output_basin Dt time intervals)*/
    DOUBLEMATRIX *PrSNW_mean;

	DOUBLEMATRIX *hsupav;
	DOUBLEMATRIX *dh_sup;
		
	DOUBLEMATRIX *error;
		
	//UMFPACK_REAL_TRIPLET *Jtriplet;
	//UMFPACK_REAL_MATRIX *Jmatrix;	
	
	DOUBLEVECTOR *Lx;
	DOUBLEVECTOR *Ux;
	
	DOUBLEVECTOR *H0;
	DOUBLEVECTOR *H1;
	DOUBLEVECTOR *dH;
	DOUBLEVECTOR *B;
	DOUBLEVECTOR *f;
	DOUBLEVECTOR *df;
	
		
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
	long i_discharge;
	long n_discharge;
	long i_plot;
	long n_plot;
	long nt_plot;
	long d_plot;
    double JD;      /*day=current Julian day during the simulation*/
    long day;       /*current day of the month during the simulation*/
    long month;       /*current month of the year during the simulation*/
    long year;     /*current year*/
    long hour;       /*current hour of day during the simulation*/
    long min;       /*current minute of hour during the simulation*/
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
	long total_channel;
 	long snowlayer_max;
	long snowlayer_inf;
	DOUBLEVECTOR *Dmin;
	DOUBLEVECTOR *Dmax;
	double Sr_glac;
	long glaclayer_max;
	DOUBLEVECTOR *Dmin_glac;
	DOUBLEVECTOR *Dmax_glac;

	short state_turb;
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
	
	DOUBLEVECTOR *JD_plots;
	
	short micromet;
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
	
	double snow_fetch;
	
	short ifill;
	short iobsint;
	double dn;
	double curve_len_scale;
	double slopewt;
	double curvewt;	
	short topoflag;
	
	short LRflag;
	
	SHORTVECTOR *vegflag;
	
	short harm_or_arit_mean;
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
	
	short state_pixel;
	short state_basin;
	short state_discharge;
			
	double DtminWb;
	long nredDtWb;
	
	double TolCG;
	long MaxiterCorr;
	short UpdateK;
	
	short channel_network;
	short bedrock;
	
	double thres_hsup;
	double thres_hchannel;
	
	//double cm_hsup_0;
	//double max_Courant_sup_sub;
	
	double Kch_b;
	double w_dx;
	
	double RelTolVWb;
	
	double snow_smin;
	double snow_smax;
	double snow_curv;
	
	double Zboundary;
	double Tboundary;
	
	double Ks_channel;
	double depr_channel;
	
} PAR;



/*---------------------------------------------------------------------------*/
typedef struct { 
	SHORTMATRIX  *type;
	LONGMATRIX	 *lnum;
	DOUBLETENSOR *Dzl;
	DOUBLETENSOR *w_liq;
	DOUBLETENSOR *w_ice;
	DOUBLETENSOR *T;
	DOUBLEMATRIX *nondimens_age;
	DOUBLEMATRIX *dimens_age;	
	DOUBLEMATRIX *max;
	DOUBLEMATRIX *average;	
	DOUBLEMATRIX *MELTED;
	DOUBLEMATRIX *SUBL;
	DOUBLEMATRIX *t_snow;
	
	DOUBLEMATRIX *DDF;
	DOUBLEMATRIX *DDF1;
	DOUBLEMATRIX *DDFvar;	
	LONGMATRIX *DDFcont;
	DOUBLEMATRIX *DDFTmin;
	DOUBLEMATRIX *DDFmeltTL0;
	DOUBLEMATRIX *DDFmelt;
	
	DOUBLEMATRIX *rho_newsnow;
	DOUBLEMATRIX *Qsub;
	DOUBLEMATRIX *Wtrans;
	DOUBLEMATRIX *Qtrans;
	DOUBLEMATRIX *Qtrans_x;
	DOUBLEMATRIX *Qtrans_y;	
	DOUBLEMATRIX *Wtot;
	DOUBLEMATRIX *Wsubl_cum;
	DOUBLEMATRIX *Wsusp_cum;
	DOUBLEMATRIX *Wtrans_cum;
	DOUBLEMATRIX *Wsubgrid_cum;
	DOUBLEMATRIX *out_bs;
	DOUBLEMATRIX *ListonSWE;
	DOUBLEMATRIX *softSWE;
	DOUBLEMATRIX *softSWE1;	
	DOUBLEMATRIX *Dplot;	
	
	DOUBLEVECTOR *CR1;
	DOUBLEVECTOR *CR2;
	DOUBLEVECTOR *CR3;
	
	LONGVECTOR *change_dir_wind;
	DOUBLEMATRIX *Psnow;	
	
} SNOW;

typedef struct { 
	LONGMATRIX	 *lnum;
	DOUBLETENSOR *Dzl;
	DOUBLETENSOR *w_liq;
	DOUBLETENSOR *w_ice;
	DOUBLETENSOR *T;	
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
	
	double **LRs;
	double *LRv;
	LONGVECTOR *LRp;
	
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
		
	DOUBLEMATRIX *Tday;	
	DOUBLEMATRIX *Tvar;
		
	long nstsrad;
	long nstlrad;
	long nstcloud;
	
} METEO;


typedef struct {
	SOIL *S;
	WATER *W;
	LAND *L;
	PAR *P;
	TOPO *T;
	CHANNEL *C;
	ENERGY *E;
	SNOW *N;
	GLACIER *G;
	METEO *M;
	TIMES *I;	
}ALLDATA;
