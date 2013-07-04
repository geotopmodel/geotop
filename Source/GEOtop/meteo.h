
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
void meteo_interp2(long *ibeg, double *out, long *col_data, long n_data, long col_JD, long col_y, double **data, double t, double JD0, long y0, double novalue);

/*----------------------------------------------------------------------------------------------------------*/
void part_snow(double prec_total, double *prec_rain, double *prec_snow, double temperature, double t_rain, double t_snow);


/*----------------------------------------------------------------------------------------------------------*/
void sat_vap_pressure(double *p, double *dp_dT, double T, double P);
void sat_vap_pressure_inv(double *T, double p, double P);
void sat_vap_pressure_2(double *p, double T, double P);
double spec_humidity(double p, double P);
void sat_spec_humidity(double *Q, double *dQ_dT, double RH, double T, double P);
void sat_spec_humidity_2(double *Q, double RH, double T, double P);
double air_density(double T, double Q, double P);
double air_cp(double T);

/*----------------------------------------------------------------------------------------------------------*/					  
