
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



//Authors: Stefano Endrizzi and Giacomo Bertoldi
//Date: 13 November 2005
//Contents: Energy balance (and also mass balance for snow and glacier)


void energy_balance(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, SNOW *snow, GLACIER *glac);

void PointEnergyBalance(long r, long c, long ns, long ng, double zmu, double zmT, double z0s, double d0s, double rz0s, double z0v, double d0v,
	double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR, double eps, double fc, double LAI, double *Wcrn, double Wcrnmax,
	double *Wcsn, double Wcsnmax, double *theta, double *land, double *root, DOUBLEVECTOR *ftcl, DOUBLEVECTOR *turb_rep, double SWin, double LWin,
	double SWv, double *LW, double *H, double *E, double *LWv, double *Hv, double *LEv, double *Etrans, double *Ts, double *Qs, double Hadd,
	double *SW, double *k, double *D, double *wi, double *wl, double *dw, double *T, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Qv,
	double *Qg, double *rh, double *rv, double *rb, double *rc, double *rh_ic, double *rv_ic, double *u_top, SOIL *sl, PAR *par, double *decay,
	double *Locc, double TsupN, double TsupNp1);

void coeff(long r, long c, long ns, long ng, long nlim, double *ad, double *ads, double *adi, double *b, double *k0, double *D, double *T0,
	double *T, double *wl, double *wi, double *dw, double EB0, double EB1, double dEB, double *SW, double Zb, double Tb, double Dt,
	double (* kfunct_snow)(double rho), double (* kfunct_glac)(double rho), SOIL *sl, PAR *par, double TsupNp1);

double calc_C(long l, long r, long c, long nsng, double *wi, double *wl, double *dw, double *D, SOIL *sl);

double calc_C0(long l, long r, long c, long nsng, double *wi, double *wl, double *dw, double *D, SOIL *sl);

double calc_k(long l, long r, long c, long ns, long ng, double *wi, double *wl, double *dw, double *T, double *D,
	double (* kfunct_snow)(double rho), double (* kfunct_glac)(double rho), SOIL *sl, PAR *par);

void EnergyFluxes(double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, double d0s, double rz0s,
	double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR, double psi, double sat, double e,
	double fc, double LAI, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn, double *dWcsn, double *theta, double **soil,
	double *land, double *root, PAR *par, DOUBLEVECTOR *ftcl, DOUBLEVECTOR *rep, double SWin, double LWin, double SWv, double *LW, double *H,
	double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv, double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs,
	double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *rh, double *rv, double *rc, double *rb, double *rh_ic, double *rv_ic, double *Qg,
	double *u_top, double *decay, double *Locc);

double surface(long r, long c, long ns, long ng, DOUBLETENSOR *snow, DOUBLETENSOR *ice, DOUBLETENSOR *sl);

void soil_properties(long r, long c, long beg, long end, double *th_l, double ***th_i, double *D, double *wl, double *wi, double *k, double *T, SOIL *sl);

double k_thermal_soil(double th_liq, double th_ice, double th_sat, double T, double k_solid);

void k_interface(long n, double *k, double *D, double *ki);

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
	double LE, double surfEB, double H, double G, double Ts, double SWin, double SWout, double SWbeam, double eps, double LWin, double LWout, double cosinc, double Precpoint);

void output_pixel(long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double Sr_soil, double Er_soil,
	double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, double Evt, double LE, double H, double surfEB,
	double G, double Tg, double A, double eps, double LAI, double LWin, double SWin, double LWout, double SWout, double epsa, double epsa_min,
	double epsa_max, double SWbeam, double SWdiff, double DTcorr, double Tdew, long n, DOUBLEVECTOR *turbulence, PAR *par, WATER *wat, ENERGY *egy,
	TOPO *top, METEO *met, SNOW *snow, GLACIER *glac, LAND *land, double Tpoint, double Ppoint, double Vpoint, double RHpoint, double prec_snow_atm,
	double prec_rain_atm, double z0soil, double z0, double d0, double SWv, double LWv, double Hv, double LEv, double Tv, double Ts, double Hg0,
	double Hg1, double Eg0, double Eg1, double fc, double rh, double rv, double rb, double rc, double rh_ic, double rv_ic, double Hv1, double LEv1,
	double Qv, double Qg, double Qa, double Qs, double u_top, double decay, double Locc);

void output_basin(long n, double prec_rain, double prec_snow, double Ta, double Tg, double Tv, double Eg, double Evt, double Eve, double LE, double H, double SW, double LW,
	double LEv, double Hv, double SWv, double LWv, double SWin, double *wat, double *en, double Dt);

void output_altrank(long ES, double prec_rain, double prec_snow, double Ts, double Er_soil, double Sr_soil, double Evt, double Eve, double LE, double H, double surfEB, double SWin, double SWout,
	double LWin, double eps, double V, double Pn, double runoff, double Er_snow, double Sr_snow, double Mr_snow, double Er_glac, double Sr_glac, double Mr_glac, double Dt, long n, double Z, double Zmin,
	double Zmax, double glacD, double glacDmin, double **out1);

void update_soil(long l, long r, long c, double evap, double *theta, double *theta_ice, double *T, SOIL *sl, PAR *par);
void energy_balance_superfast(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, SNOW *snow, GLACIER *glac);
