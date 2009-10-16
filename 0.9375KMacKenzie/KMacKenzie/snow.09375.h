
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
    
    
    
//Author: Stefano Endrizzi
//Date: 13 November 2005
//Contents: Snow subroutines

/*----------------------------------------------------------------------------------------------------------*/
double rho_newlyfallensnow(double u, double Tatm, double Tfreez);
double kidrmax_snow(double rho);

/*----------------------------------------------------------------------------------------------------------*/
void snow_compactation(long r, long c, long l, SNOW *snow, double slope, PAR *par, double *CR1, double *CR2);


/*----------------------------------------------------------------------------------------------------------*/
void snow_layer_combination(long r, long c, SNOW *snow, double Ta, long linf, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, double time);


/*----------------------------------------------------------------------------------------------------------*/
void glac_layer_combination(long r, long c, GLACIER *glac, double Ta, long max, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, double time);


/*----------------------------------------------------------------------------------------------------------*/
double DEPTH(long r, long c, LONGMATRIX *n, DOUBLETENSOR *Dz);
double get_SWE(long r, long c, LONGMATRIX *n, DOUBLETENSOR *w1, DOUBLETENSOR *w2);

/*----------------------------------------------------------------------------------------------------------*/
void snowlayer_merging(long r, long c, SNOW *snow, long l1, long l2, long lres);
double internal_energy(double w_ice, double w_liq, double T);
void from_internal_energy(long r, long c, double h, double *w_ice, double *w_liq, double *T);

/*----------------------------------------------------------------------------------------------------------*/
//void get_softsnow(SNOW *snow, double **Z);

/*----------------------------------------------------------------------------------------------------------*/
void write_snow(long r, long c, long l, SNOW *snow);
void write_snow_all(long r, long c, SNOW *snow);

/*----------------------------------------------------------------------------------------------------------*/
void set_snow(long r, long c, double *wliq, double *wice, double *T, double *D, long l1, long l2, double Dlim);
void split_layers(long r, long c, SNOW *snow, long l1);
void merge_layers(long r, long c, SNOW *snow, long l1);
void min_max_layer(long n, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, DOUBLEVECTOR *Dmin2, DOUBLEVECTOR *Dmax2, long linf);
void initialize_snow(long r, long c, long l, SNOW *snow);

/*----------------------------------------------------------------------------------------------------------*/
void show_Dminmax(long r, long c, double *Dmin, double *Dmax, long n);

/*----------------------------------------------------------------------------------------------------------*/
void update_snow_age(double Psnow, double Ts, double Dt, double *tsnow_dim, double *tsnow_nondim);

/*----------------------------------------------------------------------------------------------------------*/
double snow_albedo(double ground_alb, double snowD, double AEP, double freshsnow_alb, double C, double tsnow, double cosinc, double ( *F)(double x));
double Fzen(double cosinc);

/*----------------------------------------------------------------------------------------------------------*/
void set_shallow_snowpack(long r, long c, double Dt, SNOW *snow, double *SW, double *Mr, long *n);

/*----------------------------------------------------------------------------------------------------------*/
double k_thermal_snow_Sturm(double density);

double k_thermal_snow_Yen(double density);

void non_dimensionalize_snowage(double *snowage, double Ta);

void snow_properties(long r, long c, long beg, long end, double *D, double *wl, double *wi, double *k, double *T, double ***tdz, double ***twl, double ***twi, 
	double ***tT, double (* kfunct)(double rho));
					 
void glac2snow(long r, long c, SNOW *snow, GLACIER *glac, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax);

void liqWBsnow(long r, long c, SNOW *snow, double *Mr, double *PonS, PAR *par, double slope, double P, double *wi, double *wl, double *T, double Edt);

void iceWBsnow(long r, long c, SNOW *snow, double P, double Ta);

void glacier_init_t0(long r, long c, double Ta, GLACIER *glac, SNOW *snow, PAR *par, double time);

void WBglacier(long ns, long r, long c, GLACIER *glac, double *Mr, PAR *par, double *wi, double *wl, double *T, double Edt);

void output_snow(SNOW *snow, double **Z, PAR *par);

void find_SCA(short **LC2, SNOW *snow, PAR *par, double **Z, double t);

void snow_fluxes_H(short **LC2, SNOW *snow, DOUBLEMATRIX *Ts, DOUBLEMATRIX *H, DOUBLEMATRIX *Z, DOUBLEMATRIX *Ta, TIMES *times, PAR *par, double *VSFA, double *Hv);

double adv_efficiency(double SFA);

double theta_snow(double T);
double dtheta_snow(double T);
double max_dtheta_snow();

