
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
//Contents: Meteorological subrotines (included turbulent transfer)

/*----------------------------------------------------------------------------------------------------------*/
void meteo_distr(METEO *met, LISTON *liston, ENERGY *egy, WATER *wat, LAND *land, TOPO *top, SNOW *snow, double time, PAR *par);


/*----------------------------------------------------------------------------------------------------------*/
void meteo_vert_distr(long i, DOUBLEMATRIX *Z, METEO *met, PAR *par);


/*----------------------------------------------------------------------------------------------------------*/
void shortwave_radiation(long r, long c, double alpha, double direction, double E0, short shadow, double sky, double tau_cloud, double sa, double slope, double aspect, 
	double tau_atm, double *met_data, long *met_col, double sky_st, double A, double *SWbeam, double *SWdiff, double *cosinc, LONGMATRIX *nDt_shadow, 
	LONGMATRIX *nDt_sun);						 
												
/*----------------------------------------------------------------------------------------------------------*/						
void shadow_n(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow);
		

/*----------------------------------------------------------------------------------------------------------*/
void kriging_distr(double t, METEO_STATIONS *met_st, float *data, float novalue, DOUBLEMATRIX *Z, double int_scale, double variance, DOUBLEMATRIX *out);
void vert_distr(DOUBLEMATRIX *V, DOUBLEMATRIX *Z0, double Z_st, double V_st, double gamma, double (*f)(double a, double b, double c));
double pressure(double Dz, double P0, double gamma);
double temperature(double Dz, double T0, double gamma);

/*----------------------------------------------------------------------------------------------------------*/
void meteo_interp(double **data, double Dt, double t, double *out);

/*----------------------------------------------------------------------------------------------------------*/
void part_snow(double prec_total, double *prec_rain, double *prec_snow, double temperature, double t_rain, double t_snow);


/*----------------------------------------------------------------------------------------------------------*/
void sat_vap_pressure(double *p, double *dp_dT, double T, double P);
void sat_vap_pressure_inv(double *T, double p, double P);
double spec_humidity(double p, double P);
void sat_spec_humidity(double *Q, double *dQ_dT, double RH, double T, double P);
double air_density(double T, double Q, double P);
double air_cp(double T);

/*----------------------------------------------------------------------------------------------------------*/					  
void aero_resistance(double zmu, double zmt, double z0, double d0, double z0_z0t, double v, double Ta, double T, double Qa, double Q, double P, double gmT, 
		DOUBLEVECTOR *rep, double *rm, double *rh, double *rv, short state_turb, short MO);
			
void turbulent_fluxes(double rh, double rv, double P, double Ta, double T, double Qa, double Q, double dQdT, double *H, double *dHdT, double *E, double *dEdT);
							  					  
double Psim(double z);

double Psih(double z);		

double Zero(double z);

double PsiStab(double z);

double diff2glob(double a);

double atm_transmittance(double alpha, double P, double RH, double T, double A, double Vis, double Lozone);

void longwave_radiation(short state, double pvap, double RH, double T, double taucloud, double *eps, double *eps_max, double *eps_min);

void Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, double v, double *rm, double *rh, double *rv, DOUBLEVECTOR *w);

double cz(double zmeas, double z0, double d0, double L, double (* unstab)(double z), double (* stab)(double z));

double CZ(short state, double zmeas, double z0, double d0, double L, double (*Psi)(double z));

void Star(short a, double zmeas, double z0, double d0, double L, double u, double delta, double M, double N, double R, double *var, double *c, double *z0v, double (*Psi)(double z), double (*roughness)(double x, double y, double z) );

double roughT(double M, double N, double R);

double roughQ(double M, double N, double R);

void Businger(short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w);
void Businger2(long r, long c, short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w);

double Levap(double T);

double latent(double Ts, double Le);

short shadows_point(double **hor_height, double alpha, double azimuth, double tol_mount, double tol_flat);

void sun(double hour, double JD, double *alpha, double *direction, double *E0, double latitude, double longitude, double standard_time);

void data_meteo_for_liston(METEO *met, double novalue);

float select_data(double *meteo_t, long *col, long cod, double novalue);
