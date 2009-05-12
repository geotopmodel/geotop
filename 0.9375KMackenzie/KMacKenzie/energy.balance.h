
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

/*----------------------------------------------------------------------------------------------------------*/
void energy_balance(TIMES *times, PAR *par, LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, SNOW *snow, GLACIER *glac, LISTON *liston);
					
/*----------------------------------------------------------------------------------------------------------*/
void soil_freezing(long l, long r, long c, DOUBLEVECTOR *e0, DOUBLEVECTOR *e1, DOUBLEVECTOR *dw, DOUBLEVECTOR *T, DOUBLEVECTOR *C, SOIL *sl, double psimin);

/*----------------------------------------------------------------------------------------------------------*/
double k_thermal_snow_Sturm(double density);
double k_thermal_snow_Yen(double density);

double SB(double T);
double dSB_dT(double T);

double surface(long r, long c, long ns, long ng, DOUBLETENSOR *snow, DOUBLETENSOR *ice, DOUBLETENSOR *sl);

/*----------------------------------------------------------------------------------------------------------*/
void non_dimensionalize_snowage(double *snowage, double Ta);

void glacier_init_t0(long r, long c, double Ta, GLACIER *glac, SNOW *snow, PAR *par, double time);
			   
/*----------------------------------------------------------------------------------------------------------*/
void SolveEB(long n1, long n2, double dt, double Tg, double h0, double h1, double dhdT, double *SW, long nlim, double Gb, double *k, double *C, double *D, double *T, 
		double *dw, double *wi, double *wl, short snowtype, double Ta, SOIL *sl, long r, long c, double t, PAR *par);

double k_thermal_soil(double th_liq, double th_ice, double th_sat, double T, double k_solid);

void soil_properties(long r, long c, long beg, long end, double *th_l, double ***th_i, double *D, double *wl, double *wi, double *k, double *C, double *T, SOIL *sl);
void snow_properties(long r, long c, long beg, long end, double *D, double *wl, double *wi, double *k, double *C, double *T, double ***tdz, double ***twl, double ***twi, 
		double ***tT, double (* kfunct)(double r));
					 
void k_interface(long n, double *k, double *D, double *ki);

/*----------------------------------------------------------------------------------------------------------*/
void prepare_output(double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, double prec_rain, double prec_snow_atm, 
		ENERGY *egy, WATER *wat, SNOW *snow, GLACIER *glac, LAND *land, TOPO *top, SOIL *sl, METEO *met, TIMES *times, PAR *par, long r, long c, double A, 
		double LE, double surfEB, double H, double G, double Ts, double SWin, double SWout, double SWbeam, double eps, double LWin, double cosinc);

void output_pixel(long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double Sr_soil, double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, 
			double Mr_glac, double Er_glac, double Sr_glac, double Evt, double LE, double H, double surfEB, double G, double Tg, double A, double eps, double LAI,
			double LWin, double SWin, double LWout, double SWout, double epsa, double epsa_min, double epsa_max, double SWbeam, double SWdiff, double DTcorr, double Tdew, 
			long n, DOUBLEVECTOR *turbulence, PAR *par, WATER *wat, ENERGY *egy, TOPO *top, METEO *met, SNOW *snow, GLACIER *glac, LAND *land, double Vpoint, double RHpoint, 
			double prec_snow_atm, double prec_rain_atm, double maxstorage_c, double evap_c, double drip_c, double z0, double d0, double SWv, double LWv, double Hv, double LEv, 
			double Tv, double Ts, double Hg0, double Hg1, double Eg0, double Eg1, double fc, double rh, double rv, double rb, double rc, double rh_ic, double rv_ic,
			double Hv1, double LEv1, double Qv, double Qg, double Qa, double Qs);
									
void output_basin(long n, double prec_rain, double prec_snow, double Ta, double Tg, double Tv, double Eg, double Evt, double Eve, double LE, double H, double SW, double LW, 
		double LEv, double Hv, double SWv, double LWv, double SWin, double *wat, double *en, double Dt);
				
void output_altrank(long ES, double prec_rain, double prec_snow, double Ts, double Er_soil, double Sr_soil, double Evt, double Eve, double LE, double H, double surfEB, double SWin, double SWout, 
		double LWin, double eps, double V, double Pn, double runoff, double Er_snow, double Sr_snow, double Mr_snow, double Er_glac, double Sr_glac, double Mr_glac, double Dt, long n, double Z, double Zmin,
		double Zmax, double glacD, double glacDmin, double **out1);

/*----------------------------------------------------------------------------------------------------------*/	
double flux(long i, long icol, long **col, double **met, double k, double est);

/*----------------------------------------------------------------------------------------------------------*/	
void evaporation(long r, long c, long n, double E, double *wi, double *wl, double *D, double dt, double *Ss, double *Es, double *Eg);

/*----------------------------------------------------------------------------------------------------------*/	
void update_soil(long ntot, long r, long c, SOIL *sl, double *T, double *dw, double *wi, double *th, double *ft, double Evt, double E, PAR *par);

/*----------------------------------------------------------------------------------------------------------*/	
void liqWBsnow(long r, long c, SNOW *snow, double *Mr, double *PonS, PAR *par, double slope, double P, double *wi, double *wl, double *T);
void iceWBsnow(long r, long c, SNOW *snow, double P, double Ta);

/*----------------------------------------------------------------------------------------------------------*/	
void WBglacier(long ns, long r, long c, GLACIER *glac, double *Mr, PAR *par, double *wi, double *wl, double *T);
void glac2snow(long r, long c, SNOW *snow, GLACIER *glac, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax);

/*----------------------------------------------------------------------------------------------------------*/	
void output_snow(SNOW *snow, double **Z, PAR *par);

/*----------------------------------------------------------------------------------------------------------*/	
void output_map_plots(long r, long c, PAR *par, double t, long n, ENERGY *egy, METEO *met, SNOW *snow, double Hg, double LEg, double Hv, double LEv, double SWin, 
	double SWg, double SWv, double LWin, double LWg, double LWv, double Ts, double Tg, double Tv);
						  		
/*----------------------------------------------------------------------------------------------------------*/	
void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e, double *T, SHORTVECTOR *mf);

/*----------------------------------------------------------------------------------------------------------*/	
void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, SNOW *snow);


/*----------------------------------------------------------------------------------------------------------*/	
void root(double d, double *D, double *root_fraction);

/*----------------------------------------------------------------------------------------------------------*/	
double soil_red_evap(double psi, double T);
double red_evap(long n, double psi, double T);

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD, double thres, double z0snow, double *z0_ris, double *d0_ris, double *z0_z0t_ris);
		
void EnergyFluxes(double Tg, long r, long c, long n, double zmu, double zmT, double z0s, double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta,
	double Qa, double P, double LR, double psi, double ice, double e, double fc, double LAI, double Wc, double fwet, double fsnow, double *theta, double **soil, 
	double *land, double *root, PAR *par, DOUBLEVECTOR *ftcl, DOUBLEVECTOR *rep, double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E, 
	double *dE_dT, double *LWv, double *Hv, double *Ev, double *Evt, double *Tv, double *Qv, double *Ts, double *Qs, long cont, double *Hg0, double *Hg1, double *Eg0, 
	double *Eg1, double *rh, double *rv, double *rc, double *rb, double *rh_ic, double *rv_ic, double *Qg);
			
void PointEnergyBalance(long r, long c, long ns, long ng, double zmu, double zmT, double z0s, double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, 
	double Ta, double Qa, double P, double LR, double eps, double fc, double LAI, double Wc, double fwet, double fsnow, double *theta, double *land, double *root, 
	PAR *par, DOUBLEVECTOR *ftcl, DOUBLEVECTOR *turb_rep, double SWin, double LWin, double SWv, double *LW, double *H, double *E, double *LWv, double *Hv, double *Ev, 
	double *Evt, double *Tv, double *Ts, double *Qs, short sy, SOIL *sl, double Hadd, double *SW, double t, double Dt, double *k, double *C, double *D, double *wi, 
	double *wl, short sntype, double *dw, double *T, double *Hg0, double *Hg1, double *Eg0, double *Eg1);	

void coeff(long n, long nlim, double *ad, double *ads, double *adi, double *b, double *k, double *C, double *D, double *T, double hup, double *h, double Zboundary, double Tboundary, double Dt);

void update_coeff(long r, long c, long n1, long n2, short *occurs, double *AD, double *ADS, double *ADI, double *B, short *mf, double *ad, double *ads, double *adi, 
	double *b, double *e0, double *e, double *wi, double *wl, double *C, DOUBLEVECTOR *e0g, DOUBLEVECTOR *e1g, DOUBLEVECTOR *Cg, DOUBLEVECTOR *Tg, DOUBLEVECTOR *dwg, 
	SOIL *sl, double *SFflag, PAR *par, double Dt, double *T, double *dw);
	
void assign_values(long n1, short *mf, double *e, double *wi, double *wl, double *dw, double *T);

void check_continuity(long n, double *dw, double *wi, double *wl);

void find_tau_cloud_station(long i, METEO *met, PAR *par, double alpha, double E0, double A, double sky, double *tau, double *sa);
	