
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225-9 'Moab' - 24 Aug 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225-9 'Moab'
 
 GEOtop 1.225-9 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225-9 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef RADIATION_H
#define RADIATION_H

#include "constants.h"
#include "struct.geotop.h"
#include "meteo.h"
#include "../libraries/ascii/tabs.h"
#include "times.h"
#include "../libraries/math/util_math.h"
#include "../libraries/geomorphology/geomorphology.h"

extern long number_novalue, number_absent;

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern char *FailedRunFile;
extern long Nl, Nr, Nc;

void sun(double JDfrom0, double *E0, double *Et, double *Delta);

double SolarHeight(double JD, double latitude, double Delta, double dh);

double SolarHeight_(double JD, double *others);

double SolarHeight__(double JD, void *others);

double SolarAzimuth(double JD, double latitude, double Delta, double dh);

double SolarAzimuth_(double JD, double *others);

double SolarAzimuth__(double JD, void *others);

double TauatmCosinc(double JD, double *others);

double TauatmCosinc_(double JD, void *others);

double TauatmSinalpha(double JD, double *others);

double TauatmSinalpha_(double JD, void *others);

double Cosinc(double JD, double *others);

double Cosinc_(double JD, void *others);

double Sinalpha(double JD, double *others);

double Sinalpha_(double JD, void *others);

double Tauatm(double JD, double *others);

double Tauatm_(double JD, void *others);

void shortwave_radiation(double JDbeg, double JDend, double *others, double sin_alpha, double E0, double sky, double SWrefl_surr, 
						 double tau_cloud, short shadow, double *SWb, double *SWd, double *cos_inc_bd, 
						 double *tau_atm_sin_alpha, short *SWb_yes);
	
double diff2glob(double a);

double atm_transmittance(double X, double P, double RH, double T, double Lozone, double a, double b, double rho_g);

void longwave_radiation(short state, double pvap, double RH, double T, double k1, double k2, double taucloud, double *eps, double *eps_max, double *eps_min);

double SB(double T);

double dSB_dT(double T);
	
void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, STATEVAR_3D *snow);

double cloud_transmittance(double JDbeg, double JDend, double lat, double Delta, double dh, double RH, double T,
						   double P, double SWd, double SWb, double SW, double E0, double sky, double SWrefl_surr,
						   double Lozone, double alpha, double beta, double albedo);

double find_tau_cloud_station(double JDbeg, double JDend, long i, METEO *met, double Delta, double E0, double Et, 
							  double ST, double SWrefl_surr, double Lozone, double alpha, double beta, double albedo);

short shadows_point(double **hor_height, long hor_lines, double alpha, double azimuth, double tol_mount, double tol_flat);

void shadow_haiden(DOUBLEMATRIX *Z, double alpha, double direction, SHORTMATRIX *SH);

double find_albedo(double dry_albedo, double sat_albedo, double wat_content, double residual_wc, double saturated_wc);

void find_actual_cloudiness(double *tau_cloud, double *tau_cloud_av, short *tau_cloud_yes, short *tau_cloud_av_yes, int meteo_stat_num,
							METEO *met, double JDb, double JDe, double Delta, double E0, double Et, double ST, double A,
							double Lozone, double alpha, double beta, double albedo);
#endif
