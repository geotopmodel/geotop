
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 2.0.0
 
 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#ifndef SNOW_H
#define SNOW_H
#include "constants.h"
#include "struct.geotop.h"
#include "../libraries/ascii/rw_maps.h"
#include "times.h"
#include "output.h"

#include "PBSM.h"


//Author: Stefano Endrizzi
//Date: 13 November 2005
//Contents: Snow subroutines

/*----------------------------------------------------------------------------------------------------------*/
double rho_newlyfallensnow(double u, double Tatm);

/*----------------------------------------------------------------------------------------------------------*/
void snow_layer_combination(double a, long r, long c, Statevar3D *snow, double Ta, GeoVector<long>& inf, double SWEmax_layer, double SWEmax_tot);

/*----------------------------------------------------------------------------------------------------------*/
double DEPTH(long r, long c, GeoMatrix<long>& n, GeoTensor<double>& Dz);

/*----------------------------------------------------------------------------------------------------------*/
double internal_energy(double w_ice, double w_liq, double T);

/*----------------------------------------------------------------------------------------------------------*/
void from_internal_energy(double a, double h, double *w_ice, double *w_liq, double *T);

/*----------------------------------------------------------------------------------------------------------*/
void write_snow_all(long r, long c, Statevar3D *snow);

/*----------------------------------------------------------------------------------------------------------*/
short set_snowice_min(Statevar1D *snow, long l1, long l2, double wicemin);

/*----------------------------------------------------------------------------------------------------------*/
void update_snow_age(double Psnow, double Ts, double Dt, double Prestore, double *tsnow_nondim) ;

/*----------------------------------------------------------------------------------------------------------*/
double snow_albedo(double ground_alb, double snowD, double AEP, double freshsnow_alb, double C, double tsnow, double cosinc, double ( *F)(const double& x));

/*----------------------------------------------------------------------------------------------------------*/
double Fzen(const double& cosinc);

/*----------------------------------------------------------------------------------------------------------*/
void non_dimensionalize_snowage(double *snowage, double Ta);

/*----------------------------------------------------------------------------------------------------------*/
void WBsnow(double Dt, long ns, long r, long c, Statevar3D *snow, double *Melt, double *RainOnSnow, Par *par, double slope, double Rain, Energy *E, double Evap);

/*----------------------------------------------------------------------------------------------------------*/
void new_snow(double a, long r, long c, Statevar3D *snow, double P, double Dz, double T);

/*----------------------------------------------------------------------------------------------------------*/
void WBglacier(long ns, long ng, long r, long c, Statevar3D *glac, double *Melt, Par *par, Energy *E, double Evap);

/*----------------------------------------------------------------------------------------------------------*/
void find_SCA(Statevar3D *snow, Par *par, GeoMatrix<double>& Z, double t);

/*----------------------------------------------------------------------------------------------------------*/
double theta_snow(double a, double b, double T);

/*----------------------------------------------------------------------------------------------------------*/
double dtheta_snow(double a, double b, double T);

/*----------------------------------------------------------------------------------------------------------*/
void allocate_and_initialize_statevar_3D(Statevar3D *V, double nan, long nl, long nr, long nc);

/*----------------------------------------------------------------------------------------------------------*/
void deallocate_statevar_3D(Statevar3D *V);

/*----------------------------------------------------------------------------------------------------------*/
void allocate_and_initialize_statevar_1D(Statevar1D *V, double nan, long nl);

/*----------------------------------------------------------------------------------------------------------*/
void deallocate_statevar_1D(Statevar1D *V);

/*----------------------------------------------------------------------------------------------------------*/
short copy_statevar_from3D_to1D(long r, long c, Statevar3D *origin, Statevar1D *destination);

/*----------------------------------------------------------------------------------------------------------*/
double interpolate_snow(long r, long c, double h, long max, GeoTensor<double>& Dz, GeoTensor<double>& Q, short k);

/*----------------------------------------------------------------------------------------------------------*/
#endif
