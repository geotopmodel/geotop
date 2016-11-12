
/* STATEMENT:

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1 
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef ENERGY_BALANCE_H
#define ENERGY_BALANCE_H

#define __MATHOPTIM_H__
#include <meteoio/MeteoIO.h>
#include "constants.h"
#include "struct.geotop.h"
#include "meteo.h"
#include "snow.h"
#include "pedo.funct.h"
#include "vegetation.h"
#include "radiation.h"
#include "turbulence.h"
#include "times.h"
#include "tables.h"
#include "blowingsnow.h"
#include "meteodata.h"

#define MM 1
#define ni_en 1.E-4
#define num_iter_after_which_surfenergy_balance_not_recalculated 50
#define Tmin_surface_below_which_surfenergy_balance_recalculated -50
#define num_iter_after_which_only_neutrality 10
#define Csnow_at_T_greater_than_0 1.E20
#define ratio_max_storage_RAIN_over_canopy_to_LSAI 0.1
#define ratio_max_storage_SNOW_over_canopy_to_LSAI 5.0
#define Tmax 100.0
#define Tmin -90.0


//short EnergyBalance(double Dt, double JD0, double JDb, double JDe, SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G, STATE_VEG *V, DOUBLEVECTOR *snowage, ALLDATA *A, double *W, std::vector<mio::MeteoData>& vecmeteo);
short EnergyBalance(double Dt, double JD0, double JDb, double JDe, SoilState *L, SoilState *C, Statevar3D *S, Statevar3D *G, StateVeg *V, GeoVector<double>& snowage, AllData *A, double *W, std::vector<mio::MeteoData>& vecmeteo);

//short PointEnergyBalance(long i, long r, long c, double Dt, double JDb, double JDe, SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G, STATE_VEG *V,
//						 DOUBLEVECTOR *snowage, ALLDATA *A, double E0, double Et, double Dtplot, double W, FILE *f, double *SWupabove_v, double *Tgskin);

//short PointEnergyBalance(long i, long r, long c, double Dt, double JDb, double JDe, SOIL_STATE *L, SOIL_STATE *C, Statevar3D *S, Statevar3D *G, STATE_VEG *V,
//						 DOUBLEVECTOR *snowage, ALLDATA *A, double E0, double Et, double Dtplot, double W, FILE *f, double *SWupabove_v, double *Tgskin);

short PointEnergyBalance(long i, long r, long c, double Dt, double JDb, double JDe, SoilState *L, SoilState *C, Statevar3D *S, Statevar3D *G, StateVeg *V,
		 GeoVector<double>& snowage, AllData *A, double E0, double Et, double Dtplot, double W, FILE *f, double *SWupabove_v, double *Tgskin);



short SolvePointEnergyBalance(const short surfacemelting, double Tgd, double EBd, double Convd, short surfacebalance, const double t, const double Dt, const long i, const long j, const long r, const long c, SoilState * SL, SoilState *SC,
						StateVeg *V, Energy *egy, Land *land, Soil *sl, Channel *cnet, Topo *topog, Par *par, const long ns, const long ng, const double zmu, const double zmT,
						const double z0s, const double d0s, const double rz0s,const double z0v, const double d0v, const double rz0v, const double hveg, const double v,
						const double Ta, const double Qa, const double P, const double LR, const double eps, const double fc, const double LSAI,
						const double decaycoeff0, double *Wcrn, const double Wcrnmax, double *Wcsn, const double Wcsnmax, const double SWin, const double LWin,
						const double SWv, double *LW, double *H, double *E, double *LWv, double *Hv, double *LEv, double *Etrans, double *Ts, double *Qs, 
						const double Eadd, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Qv, double *Qg, double *Lob, double *rh, double *rv,
						double *rb, double *rc, double *ruc, double *u_top, double *decay, double *Locc, double *LWup_ab_v, long *lpb,  double *dUsl);


void update_soil_land(long nsurf, long n, long i, long r, long c, double fc, double Dt, Energy *egy, GeoTensor<double>& pa, long sy, SoilState *S, GeoTensor<double>& ET, GeoMatrix<double>& th);

//void update_soil_channel(long nsurf, long n, long ch, double fc, double Dt, Energy *egy, double **pa, SOIL_STATE *S, DOUBLEMATRIX *ET, DOUBLEMATRIX *th);
  void update_soil_channel(long nsurf, long n, long ch, double fc, double Dt, Energy *egy, GeoTensor<double>& pa, long sy, SoilState *S, GeoMatrix<double>& ET, GeoMatrix<double>& th);


//void update_F_energy(long nbeg, long nend, DOUBLEVECTOR *F, double w, DOUBLEVECTOR *K, double *T);
void update_F_energy(long nbeg, long nend, GeoVector<double>& F, double w, const GeoVector<double>& K, const double *T);

//void update_diag_dF_energy(long nbeg, long nend, DOUBLEVECTOR *dF, double w, DOUBLEVECTOR *K);
void update_diag_dF_energy(long nbeg, long nend, GeoVector<double>& dF, double w, const GeoVector<double>& K);

//double calc_C(long l, long nsng, double a, double *wi, double *wl, double *dw, double *D, double **pa);
double calc_C(long l, long nsng, double a, double *wi, double *wl, double *dw, double *D, const GeoTensor<double>& pa, long sy);
	
//void EnergyFluxes(double t, double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s,
//				  double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR,
//				  double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax,
//				  double *dWcrn, double *dWcsn, double *theta, double **soil, double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer,
//				  double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv,
//				  double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0,
//				  double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb, double *ruc, double *rh_g,
//				  double *rv_g, double *Qg, double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T,
//				  DOUBLEVECTOR *soil_evap_layer_bare, DOUBLEVECTOR *soil_evap_layer_veg, double point_elev, double point_slope, double point_aspect, double point_sky, int point_lc, long point_sy, double snowD);


void EnergyFluxes(double t, double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
				  double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR,
				  double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax, 
			   double *dWcrn, double *dWcsn, double *theta, GeoTensor<double>& soil, long sy, const GeoMatrix<double>& land, long lu, const GeoMatrix<double>& root, Par *par, GeoVector<double>& soil_transp_layer,
				  double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv,
				  double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, 
				  double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb, double *ruc, double *rh_g, 
				  double *rv_g, double *Qg, double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T, 
				  GeoVector<double>& soil_evap_layer_bare, GeoVector<double>& soil_evap_layer_veg, double point_elev, double point_slope, double point_aspect, double point_sky, int point_lc, long point_sy, double snowD);


//void EnergyFluxes_no_rec_turbulence(double t, double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s,
//									double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR,
//									double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax,
//									double *dWcrn, double *dWcsn, double *theta, double **soil, double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer,
//									double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv,
//									double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0,
//									double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb, double *ruc, double *rh_g,
//									double *rv_g, double *Qg, double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T,
//									DOUBLEVECTOR *soil_evap_layer_bare, DOUBLEVECTOR *soil_evap_layer_veg, short flagTmin, long cont);

void EnergyFluxes_no_rec_turbulence(double t, double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
									double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR,
									double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax, 
									double *dWcrn, double *dWcsn, double *theta, GeoTensor<double>& soil, long sy, const GeoMatrix<double>& land, long lu,const GeoMatrix<double>& root, Par *par, GeoVector<double>& soil_transp_layer,
									double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv,
									double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, 
									double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb, double *ruc, double *rh_g, 
									double *rv_g, double *Qg, double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T, 
									GeoVector<double>& soil_evap_layer_bare, GeoVector<double>& soil_evap_layer_veg, double point_sky, short flagTmin, long cont);

double k_thermal(short snow, short a, double th_liq, double th_ice, double th_sat, double k_solid);

double flux(long i, long icol, double **met, double k, double est);

//void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e, double *T, SHORTVECTOR *mf);
void check_errors(long r, long c, long n, GeoVector<double>& adi, GeoVector<double>& ad, GeoVector<double>& ads, GeoVector<double>& b, GeoVector<double>& e,
				  double *T, GeoVector<short>& mf);


double soil_red_evap(double psi, double T);

double red_evap(long n, double psi, double T);

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD, double thres, double z0snow, double *z0_ris, double *d0_ris, double *z0_z0t_ris);

void check_continuity(long n, double *dw, double *wi, double *wl);

//void merge(double a, DOUBLEVECTOR *ice, DOUBLEVECTOR *liq, DOUBLEVECTOR *Temp, DOUBLEVECTOR *D, long l, long lup, long ldw, long tot);
  void merge(double a, GeoVector<double>& ice, GeoVector<double>& liq, GeoVector<double>& Temp, GeoVector<double>& D, long l, long lup, long ldw, long tot);

void sux_minus6_condition(double ic, double wa, double rho, double D1, Energy *E);
#endif
