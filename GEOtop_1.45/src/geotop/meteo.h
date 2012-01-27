
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

/*----------------------------------------------------------------------------------------------------------*/
void meteo_distr(short update, long k, long *line, long lineLR, METEO *met, WATER *wat, TOPO *top, PAR *par, double JDbeg, double JDend);

/*----------------------------------------------------------------------------------------------------------*/
double pressure(double Z, double Z0, double P0);

/*----------------------------------------------------------------------------------------------------------*/
double temperature(double Z, double Z0, double T0, double gamma);

/*----------------------------------------------------------------------------------------------------------*/
void part_snow(double prec_total, double *prec_rain, double *prec_snow, double temperature, double t_rain, double t_snow);

/*----------------------------------------------------------------------------------------------------------*/
double SatVapPressure(double T, double P);
void SatVapPressure_2(double *e, double *de_dT, double T, double P);
double TfromSatVapPressure(double e, double P);
double SpecHumidity(double e, double P);
void SpecHumidity_2(double *Q, double *dQ_dT, double RH, double T, double P);
double VapPressurefromSpecHumidity(double Q, double P);
double Tdew(double T, double RH, double Z);
double RHfromTdew(double T, double Tdew, double Z);
double air_density(double T, double Q, double P);
double air_cp(double T);

/*----------------------------------------------------------------------------------------------------------*/




