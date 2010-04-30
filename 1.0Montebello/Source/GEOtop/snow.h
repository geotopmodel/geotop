
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
    
    
    
//Author: Stefano Endrizzi
//Date: 13 November 2005
//Contents: Snow subroutines

/*----------------------------------------------------------------------------------------------------------*/
double rho_newlyfallensnow(double u, double Tatm, double Tfreez);
double kidrmax_snow(double rho);

/*----------------------------------------------------------------------------------------------------------*/
void snow_compactation(long r, long c, long l, SNOW *snow, double slope, PAR *par, double *CR1, double *CR2);


/*----------------------------------------------------------------------------------------------------------*/
void snow_layer_combination(double a, long r, long c, SNOW *snow, double Ta, long linf, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, double time);


/*----------------------------------------------------------------------------------------------------------*/
void glac_layer_combination(long r, long c, GLACIER *glac, double Ta, long max, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, double time);


/*----------------------------------------------------------------------------------------------------------*/
double DEPTH(long r, long c, LONGMATRIX *n, DOUBLETENSOR *Dz);
double get_SWE(long r, long c, LONGMATRIX *n, DOUBLETENSOR *w1, DOUBLETENSOR *w2);

/*----------------------------------------------------------------------------------------------------------*/
void snowlayer_merging(double a, long r, long c, SNOW *snow, long l1, long l2, long lres);
double internal_energy(double w_ice, double w_liq, double T);
void from_internal_energy(double a, long r, long c, double h, double *w_ice, double *w_liq, double *T);

/*----------------------------------------------------------------------------------------------------------*/
//void get_softsnow(SNOW *snow, double **Z);

/*----------------------------------------------------------------------------------------------------------*/
void write_snow(long r, long c, long l, SNOW *snow);
void write_snow_all(long r, long c, SNOW *snow);

/*----------------------------------------------------------------------------------------------------------*/
void set_snow(double a, long r, long c, double *wliq, double *wice, double *T, double *D, long l1, long l2, double Dlim);
void split_layers(long r, long c, SNOW *snow, long l1);
void merge_layers(double a, long r, long c, SNOW *snow, long l1);
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
					 
void glac2snow(long r, long c, SNOW *snow, GLACIER *glac, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax);

void liqWBsnow(long r, long c, SNOW *snow, double *Mr, double *PonS, PAR *par, double slope, double P, double *wi, double *wl, double *T, double Edt);

void iceWBsnow(double a, long r, long c, SNOW *snow, double P, double Ta);

void glacier_init_t0(long r, long c, double Ta, GLACIER *glac, SNOW *snow, PAR *par, double time);

void WBglacier(long ns, long r, long c, GLACIER *glac, double *Mr, PAR *par, double *wi, double *wl, double *T, double Edt);

void output_snow(SNOW *snow, double **Z, PAR *par);

void find_SCA(SNOW *snow, PAR *par, double **Z, double t);

double theta_snow(double a, double T);
double dtheta_snow(double a, double T);
double max_dtheta_snow(double a);

void remove_snow(long r, long c, SNOW *snow, long n);
void add_snow(long r, long c, double Dsnow, double rho_snow, double Tsnow, SNOW *snow, PAR *par);
