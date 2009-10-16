
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
void meteo_distr(METEO *met, ENERGY *egy, WATER *wat, TOPO *top, SNOW *snow, double time, PAR *par);


/*----------------------------------------------------------------------------------------------------------*/
void meteo_vert_distr(long i, DOUBLEMATRIX *Z, METEO *met, PAR *par);


/*----------------------------------------------------------------------------------------------------------*/
void kriging_distr(double t, METEO_STATIONS *met_st, float *data, float novalue, DOUBLEMATRIX *Z, double int_scale, double variance, DOUBLEMATRIX *out);
void vert_distr(DOUBLEMATRIX *V, DOUBLEMATRIX *Z0, double Z_st, double V_st, double gamma, double (*f)(double a, double b, double c));
double pressure(double Dz, double P0, double gamma);
double temperature(double Dz, double T0, double gamma);

/*----------------------------------------------------------------------------------------------------------*/
void meteo_interp(double **data, double Dt, double t, double *out);
void meteo_interp2(double *out, LONGVECTOR *col_data, long col_JD, long col_y, double **data, double t, PAR *par);

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
