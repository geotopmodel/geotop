
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#ifndef METEO_H
#define METEO_H
#include "constants.h"
#include "struct.geotop.h"
#include "meteodata.h"
#include "meteodistr.h"
#include "../libraries/ascii/tabs.h"
#include "times.h"
#include "snow.h"
#include "../meteoio_plugin/meteoioplugin.h"

extern long number_novalue;
//number_absent;
//extern T_INIT *UV;
extern TInit *UV;
//extern char *WORKING_DIRECTORY;
//extern char *logfile;
extern std::string logfile;

extern long  Nr, Nc;//Nl
//extern long i_sim;
/*----------------------------------------------------------------------------------------------------------*/
//void meteo_distr(long *line, long lineLR, METEO *met, WATER *wat, TOPO *top, PAR *par, double JD0, double JDbeg, double JDend);
void meteo_distr(long *line, long lineLR, Meteo *met, Water *wat, Topo *top, Par *par, double JD0, double JDbeg, double JDend);

/*----------------------------------------------------------------------------------------------------------*/
double pressure(double Z);

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

#endif


