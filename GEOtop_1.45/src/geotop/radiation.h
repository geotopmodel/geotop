
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

void shortwave_radiation(double JDbeg, double JDend, double *others, double sin_alpha, double E0, double sky, double A, 
						 double tau_cloud, short shadow, double *SWb, double *SWd, double *cos_inc_bd, 
						 double *tau_atm_sin_alpha, short *SWb_yes);
	
double diff2glob(double a);

double atm_transmittance(double X, double P, double RH, double T);

void longwave_radiation(short state, double pvap, double RH, double T, double taucloud, double *eps, double *eps_max, double *eps_min);

double SB(double T);

double dSB_dT(double T);
	
void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, STATEVAR_3D *snow);

double cloud_transmittance(double JDbeg, double JDend, double lat, double Delta, double dh, double RH, double T,
						   double P, double SWd, double SWb, double SW, double E0, double sky, double A);

double find_tau_cloud_station(double JDbeg, double JDend, long i, METEO *met, double Delta, double E0, double Et, 
							  double ST, double A);

short shadows_point(double **hor_height, long hor_lines, double alpha, double azimuth, double tol_mount, double tol_flat);

void shadow_haiden(TOPO *top, double alpha, double direction, SHORTMATRIX *shadow);

double find_albedo(double dry_albedo, double sat_albedo, double wat_content, double residual_wc, double saturated_wc);

void find_actual_cloudiness(double *tau_cloud, double *tau_cloud_av, short *tau_cloud_yes, short *tau_cloud_av_yes, 
							METEO *met, double JDb, double JDe, double Delta, double E0, double Et, double ST, double A);
