
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie 

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Matteo Dall'Amico

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
    
    
//Structure including all the par, matrices, vectors and tensors handled by MicroMet and SnowTrans3D
typedef struct {
	int nx;
	int ny;
	float dt;	
	float deltax;
	float deltay;	
	float undef;
	float topoflag;
	int nvegtypes;	//not in call subroutine
	
	float *topo_land;
	float *topo;
	float *xlat_grid;
	float *vegtype;
	float *snow_d_init;
	float *soft_snow_d;
	float *snowdepth;
	float *vegsnowdepth;
	float *veg_z0;
	float *forest_LAI;
	float *forest_LAIw;
	float *forest_LAIs;	
	float *sum_sprec;
	float *gap_frac;
	
	double xmn;
	double ymn;
	int iyear_init;
	int imonth_init;
	int iday_init;
	float xhour_init;
	int ifill;
	int iobsint;
	float dn;
	float curve_len_scale;
	float slopewt;
	float curvewt;
	int i_tair_flag;
	int i_rh_flag;
	int i_wind_flag;
	int i_solar_flag;
	int i_prec_flag;	
	int i_longwave_flag;
	int lapse_rate_user_flag;
	int iprecip_lapse_rate_user_flag; 	
	float use_shortwave_obs;
	float use_longwave_obs;
	float use_sfc_pressure_obs;
	float calc_subcanopy_met;
	float cloud_frac_factor;
	float barnes_lg_domain;
	float windspd_min;
	float windspd_flag;
	float winddir_flag;
	int nstat;	
	double *met_stE;
	double *met_stN;
	float *met_stZ;
	
	float *curvature;
	float *terrain_slope;
	float *slope_az;
	float *tair_grid;
	float *rh_grid;
	float *uwind_grid;
	float *vwind_grid;
	float *Qsi_grid;
	float *prec_grid;
	float *windspd_grid;
	float *winddir_grid;
	float *sprec;
	float *Qli_grid;
	float *sfc_pressure;	
	int *k_stn;
	
	float bc_flag;
	float bs_flag;
	float C_z;
	float dz_susp;
	float fall_vel;
	float fetch;
	float gravity;
	float h_const;
	float ht_rhobs;
	float ht_windobs;
	float ro_air;
	float ro_snow;
	float ro_water;
	float snow_z0;
	float subgrid_flag;
	float tp_scale;
	float Up_const;
	float Ur_const;
	float vonKarman;
	float wind_min;
	float xmu;
	float ztop_susp;
	float erosion_dist;
	float tabler_dir;
	float slope_adjust;
	float Utau_t_const;
	float Utau_t_flag;
	
	float *conc_salt;
	float *dh_salt;
	float *dh_salt_u;
	float *dh_salt_v;
	float *dh_susp;
	float *dh_susp_u;
	float *dh_susp_v;
	float *h_star;
	float *Qsalt;
	float *Qsalt_u;
	float *Qsalt_v;
	float *Qsalt_max;
	float *Qsalt_maxu;
	float *Qsalt_maxv;
	float *Qsubl;
	float *Qsusp;
	float *Qsusp_u;
	float *Qsusp_v;
	float *wbal_salt;
	float *wbal_susp;
	float *wbal_qsubl;
	float *wbal_subgrid;
	float *Utau;
	float *Utau_t;
	float *z_0;
	float *ro_soft_snow_old;
	float *ro_soft_snow;	
	float *sprec_geotop;
	float *ro_nsnow;
	
 } LISTON;

//initialize Liston variables (to be called every time step)
void initialize_liston(DOUBLEMATRIX *landtype, SHORTMATRIX *LU, DOUBLEMATRIX *Z, PAR *par, METEO_STATIONS *st, DOUBLEMATRIX *snowD_init, double rhosnow0, LISTON *liston);
void deallocate_liston(LISTON *liston);

void call_MicroMet(double t, double dt, DOUBLEMATRIX *Z, METEO *met, ENERGY *egy, SNOW *snow, DOUBLEMATRIX *prec, DOUBLEVECTOR *LAI, LISTON *liston, PAR *par);
void call_SnowTrans(double t, double dt, DOUBLEMATRIX *Z, SNOW *snow, DOUBLEMATRIX *Psnow, LISTON *liston);

DOUBLEMATRIX *liston2matrix(float *vector, long nr, long nc);
void liston2matrix_2(DOUBLEMATRIX *M, float *vector);
float *matrix2liston(DOUBLEMATRIX *M);
void matrix2liston_2(float *vector, DOUBLEMATRIX *M);
float *vector2liston(DOUBLEVECTOR *V);
double *vector2listondouble(DOUBLEVECTOR *V);
void liston2doublevector(DOUBLEVECTOR *V, float *vector);

void extend_topography(DOUBLEMATRIX *M, double novalue);
void find_the_nearest(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);
short no_novalue(long r, long c, DOUBLEMATRIX *M, double novalue, long *rr, long *cc);
void m_dmatrix(DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double f, double novalue);

void set_to_point_simulation(DOUBLEMATRIX *Mbigger, DOUBLEMATRIX *Msmall, LONGVECTOR *r, LONGVECTOR *c);
