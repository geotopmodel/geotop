
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

void shortwave_radiation(long r, long c, double alpha, double direction, double E0, short shadow, double sky, double tau_cloud, double sa, double slope, double aspect, 
	double tau_atm, double *met_data, long *met_col, double sky_st, double A, double *SWbeam, double *SWdiff, double *cosinc, LONGMATRIX *nDt_shadow, 
	LONGMATRIX *nDt_sun);
	
double diff2glob(double a);

double atm_transmittance(double X, double P, double RH, double T);

void longwave_radiation(short state, double pvap, double RH, double T, double taucloud, double *eps, double *eps_max, double *eps_min);

double SB(double T);

double dSB_dT(double T);

short shadows_point(double **hor_height, double alpha, double azimuth, double tol_mount, double tol_flat);

void sun(double JD, double *alpha, double *direction, double *E0, double latitude, double longitude, double standard_time);

void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, SNOW *snow);

void find_tau_cloud_station(long i, METEO *met, PAR *par, double alpha, double E0, double sky, double A, double *tau, double *sa);

void shadow_haiden(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow);
