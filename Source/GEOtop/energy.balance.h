
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



void energy_balance(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, SNOW *snow, GLACIER *glac);

void PointEnergyBalance(long r, long c, ENERGY *egy, LAND *land, SOIL *sl, PAR *par, long ns, long ng, double zmu, double zmT, 
	double z0s, double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, 
	double P, double LR, double eps, double fc, double LSAI, double decaycoeff0, double *Wcrn, double Wcrnmax, double *Wcsn, 
	double Wcsnmax, double SWin, double LWin, double SWv,double *LW, double *H, double *E, double *LWv, double *Hv, double *LEv,
	double *Etrans, double *Ts, double *Qs, double Hadd, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Qv, 
	double *Qg, double *Lobukhov, double *rh, double *rv, double *rb, double *rc, double *r_uc, double *u_top, 
	double *decay, double *Locc, double *LWup_above_v);

void update_F_energy(long n, DOUBLEVECTOR *F, double w, DOUBLEVECTOR *K, double *T);

void update_diag_dF_energy(long n, DOUBLEVECTOR *dF, double w, DOUBLEVECTOR *K);

double calc_C(long l, long r, long c, long nsng, double a, double *wi, double *wl, double *dw, double *D, double *T, SOIL *sl);

double calc_k(long l, long r, long c, long ns, long ng, double *wi, double *wl, double *dw, double *T, double *D, 
	double (* kfunct_snow)(double rho), double (* kfunct_glac)(double rho), SOIL *sl, PAR *par);	

void EnergyFluxes(double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
	double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR,
	double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax, 
	double *dWcrn, double *dWcsn, double *theta, double **soil, double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer, 
	double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv,
	double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, 
	double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb, double *r_uc, double *rh_g, 
	double *rv_g, double *Qg, double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T, 
	DOUBLEVECTOR *soil_evap_layer_bare, DOUBLEVECTOR *soil_evap_layer_veg );

void EnergyFluxes_no_rec_turbulence(double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
	double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR,
	double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax, 
	double *dWcrn, double *dWcsn, double *theta, double **soil, double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer, 
	double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv,
	double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, 
	double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb, double *r_uc, double *rh_g, 
	double *rv_g, double *Qg, double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T, 
	DOUBLEVECTOR *soil_evap_layer_bare, DOUBLEVECTOR *soil_evap_layer_veg );

double surface(long r, long c, long ns, long ng, DOUBLETENSOR *snow, DOUBLETENSOR *ice, DOUBLETENSOR *sl);

double k_thermal_soil(double th_liq, double th_ice, double th_sat, double T, double k_solid);

double flux(long i, long icol, long **col, double **met, double k, double est);

void output_map_plots(long r, long c, PAR *par, double t, long n, ENERGY *egy, METEO *met, SNOW *snow, double Hg, double LEg, double Hv, double LEv, double SWin, 
	double SWg, double SWv, double LWin, double LWg, double LWv, double Ts, double Tg, double Tv);
	
void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e, double *T, SHORTVECTOR *mf);

double soil_red_evap(double psi, double T);

double red_evap(long n, double psi, double T);

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD, double thres, double z0snow, double *z0_ris, double *d0_ris, double *z0_z0t_ris);

void check_continuity(long n, double *dw, double *wi, double *wl);

void prepare_output(double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, double prec_rain, double prec_snow_atm, 
	ENERGY *egy, WATER *wat, SNOW *snow, GLACIER *glac, LAND *land, TOPO *top, SOIL *sl, METEO *met, TIMES *times, PAR *par, long r, long c, double A, 
	double LE, double surfEB, double H, double G, double Ts, double SWin, double SW, double SWbeam, double eps, double LWin, double LW, double cosinc, double Precpoint);

void output_pixel(long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double Sr_soil, double Er_soil, 
	double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, double Evt, double LE, double H, 
	double surfEB, double G, double Tg, double A, double eps, double LSAI, double LWin, double SWin, double LW, double SW, double epsa, 
	double epsa_min, double epsa_max, double SWbeam, double SWdiff, double Tdew, long n, double Lobukhov, PAR *par, WATER *wat, ENERGY *egy, 
	TOPO *top, METEO *met, SNOW *snow, GLACIER *glac, LAND *land, SOIL *sl, double Tpoint, double Ppoint, double Vpoint, double RHpoint, 
	double prec_snow_atm, double prec_rain_atm, double z0soil, double z0, double d0, double SWv, double LWv, double Tv, double Ts, double Hg0, 
	double Hg1, double Eg0, double Eg1, double fc, double rh, double rv, double rb, double rc, double r_uc, double Hv, double LEv, 
	double Qv, double Qg, double Qa, double Qs, double u_top, double decay, double Locc, double SWup_above_v, double LWup_above_v);
		
void output_basin(long n, double npixel, double prec_rain, double prec_snow, double prec_rain_atm, double prec_snow_atm, double Ta, double Tg, 
	double Tv, double Eg, double Evt, double LE, double H, double SW, double LW, double LEv, double Hv, double SWv, double LWv, double SWin, double LWin,
	double Dt);
	