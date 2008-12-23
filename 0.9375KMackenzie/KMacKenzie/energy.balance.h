
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion MacLavagna 

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 MacLavagna. 
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
void energy_balance(TIMES *times, PAR *par, LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, 
					SNOW *snow, GLACIER *glac, LISTON *liston);
					
/*----------------------------------------------------------------------------------------------------------*/
void soil_freezing(long l, long r, long c, DOUBLEVECTOR *e0, DOUBLEVECTOR *e1, DOUBLEVECTOR *dw, DOUBLEVECTOR *T, DOUBLEVECTOR *C, SOIL *sl, double psimin);

/*----------------------------------------------------------------------------------------------------------*/
void cloudinessSWglob(double alpha, double E0, double P, double RH, double T, double SW, double A, double sky, double Vis, double L03, double *fcloud, double *tau_atm, double *kd, double *sa, double *tau_cloud);
void cloudinessSWdiff(double alpha, double E0, double P, double RH, double T, double SWd, double SWb, double A, double sky, double Vis, double L03, double *fcloud, double *tau_atm, double *kd, double *sa, double *tau_cloud);

void cloudiness(long i, METEO *met, double A, double alpha, double direction, double E0, double *tau_atm, short *need_cloud_next_day, 
		double *cloud_p_v, double *cloud_p_t, double *cloud_n_v, double *cloud_n_t, short *check_rain_night, TIMES *times, PAR *par, double *fcloud);

void check_clouds_nextday(double t, long i, METEO_STATIONS *stm, long **col, double ***data_meteo, double *cloud_n_t, double *cloud_n_v, short *check_rain_night, double A, double **hor_height, 
						  double tol_mount, double tol_flat, PAR *par);

double interpolate_function(double x1, double y1, double x2, double y2, double x);

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
void SolveEB(long n1, long n2, double dt, double h, double dhdT, DOUBLEVECTOR *abs, long nlim, double Gb, double *k, double *C, double *D, double *T, double *dw, double *wi, double *wl, 
			 short snowtype, double Ta, SOIL *sl, long r, long c, double t, PAR *par);

double k_thermal_soil(double th_liq, double th_ice, double th_sat, double T, double k_solid);

void soil_properties(long r, long c, long beg, long end, double *th_l, double ***th_i, double *D, double *wl, double *wi, double *k, double *C, double *T, SOIL *sl);
void snow_properties(long r, long c, long beg, long end, double *D, double *wl, double *wi, double *k, double *C, double *T, double ***tdz, double ***twl, double ***twi, 
					 double ***tT, double (* kfunct)(double r));
					 
void k_interface(long n, double *k, double *D, double *ki);

/*----------------------------------------------------------------------------------------------------------*/
void prepare_output(double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, double prec_rain, double prec_snow_atm, 
			ENERGY *egy, WATER *wat, SNOW *snow, GLACIER *glac, LAND *land, TOPO *top, SOIL *sl, METEO *met, TIMES *times, PAR *par, long r, long c, double tausn, double albedo_value, 
			double LE, double surfEB, double H, double G, double Ts, double SWin, double SWbeam, double eps, double LWin, double alpha, double cosinc);

void output_pixel(long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double Sr_soil, double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, 
			double Mr_glac, double Er_glac, double Sr_glac, double Epc, double fec, double ftc, double LE, double H, double surfEB, double G, double Ts, double albedo_value, double eps, 
			double LWin, double SWin, double LWout, double SWout, double epsa, double epsa_min, double epsa_max, double SWbeam, double SWdiff, double DTcorr, double Tdew, 
			long n, DOUBLEVECTOR *turbulence, PAR *par, WATER *wat, ENERGY *egy, TOPO *top, METEO *met, SNOW *snow, GLACIER *glac, LAND *land, double Vpoint, double RHpoint, 
			double prec_snow_atm, double prec_rain_atm, double maxstorage_c, double evap_c, double drip_c, double z0, double d0);
			
void output_basin(double prec_rain, double prec_snow, double Ts, double Er_soil, double Sr_soil, double Epc, double fec, double ftc, double LE, double H, double surfEB, 
		double Eimm, double SWin, double LWin, double albedo, double eps, double V, double Dt, long n, double *wat, double *en);
				
void output_altrank(long ES, double prec_rain, double prec_snow, double Ts, double Er_soil, double Sr_soil, double Epc, double fec, double ftc, double LE, double H, double surfEB, 
		double Eimm, double SWin, double LWin, double albedo, double eps, double V, double Pn, double runoff, double Er_snow, double Sr_snow, double Mr_snow, double Er_glac,
		double Sr_glac, double Mr_glac, double Dt, long n, double Z, double Zmin, double Zmax, double glacD, double glacDmin, double **out);

/*----------------------------------------------------------------------------------------------------------*/	
double flux(long i, long icol, long **col, double **met, double k, double est);

/*----------------------------------------------------------------------------------------------------------*/	
void evaporation(long r, long c, long n, double E, double Rg, double *wi, double *wl, double *D, double dt, double *Ss, double *Es, double *Eg);

/*----------------------------------------------------------------------------------------------------------*/	
void update_soil(long ntot, long r, long c, SOIL *sl, double *T, double *dw, double *wi, double *th, double fc, double *ft, double Ep, double Eg, PAR *par);

/*----------------------------------------------------------------------------------------------------------*/	
void liqWBsnow(long r, long c, SNOW *snow, double *Mr, double *PonS, PAR *par, double slope, double P, double *wi, double *wl, double *T);
void iceWBsnow(long r, long c, SNOW *snow, double P, double Ta);

/*----------------------------------------------------------------------------------------------------------*/	
void WBglacier(long ns, long r, long c, GLACIER *glac, double *Mr, PAR *par, double *wi, double *wl, double *T);
void glac2snow(long r, long c, SNOW *snow, GLACIER *glac, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax);

/*----------------------------------------------------------------------------------------------------------*/	
void output_snow(SNOW *snow, double **Z, PAR *par);

/*----------------------------------------------------------------------------------------------------------*/	
void output_map_plots(long r, long c, PAR *par, double t, long n, ENERGY *egy, METEO *met, SNOW *snow, double H, double LE, double SWin, double SWout, double LWin, double LWout, double Ts);
					  
/*----------------------------------------------------------------------------------------------------------*/	
void evaporation_soil(long nsnow, double E, double dEdT, double theta, double res, double sat, double Ratm, double Ts, double *RG, double *LE_s, double *dLEdT_s, double *Evc, 
    double *Dc, double *Msc, double *Epc, double *dEpcdT, double *fec, double *ftc, DOUBLEVECTOR *ftcl, double *LE_c, double *dLEdT_c);
	
void evaporation_canopy(long r, long c, double Ptot, double *theta, double Ratm, double rho, double Tv, double Ta, double Qa, double Pa, double SW, double LAI,
	double *land, double *root, double *Wc, double *LE, double *dLEdT, double *Evc, double *Dc, double *Msc, double *Epc, double *dEpcdT, double *fec, double *ftc, 
	double *ftcl, PAR *par);
		
/*----------------------------------------------------------------------------------------------------------*/	
void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e, double *T, SHORTVECTOR *mf);

/*----------------------------------------------------------------------------------------------------------*/	
void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, SNOW *snow);

/*----------------------------------------------------------------------------------------------------------*/	
void snow_fluxes_H(SNOW *snow, DOUBLEMATRIX *Ts, DOUBLEMATRIX *H, DOUBLEMATRIX *Z, DOUBLEMATRIX *Ta, TIMES *times, PAR *par, double *VSFA, double *Hv);
double adv_efficiency(double SFA);

/*----------------------------------------------------------------------------------------------------------*/	
void root(double d, double *D, double *root_fraction);

/*----------------------------------------------------------------------------------------------------------*/	
