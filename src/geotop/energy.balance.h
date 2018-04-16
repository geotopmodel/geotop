
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

short EnergyBalance(double Dt, double JD0, double JDb, double JDe,
                    SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G, STATE_VEG *V,
                    DOUBLEVECTOR *snowage,
                    ALLDATA *A, double *W);

short PointEnergyBalance(long i, long r, long c, double Dt, double JDb,
                         double JDe, SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G,
                         STATE_VEG *V,
                         DOUBLEVECTOR *snowage, ALLDATA *A, double E0, double Et, double Dtplot,
                         double W, FILE *f, double *SWupabove_v, double *Tgskin);

short SolvePointEnergyBalance(short surfacemelting, double Tgd,
                              double Tbottom, double EBd, double Convd, short surfacebalance, double t,
                              double Dt, long i, long j, long r, long c, SOIL_STATE *SL,
                              SOIL_STATE *SC, STATE_VEG *V, ENERGY *egy, LAND *land, SOIL *sl,
                              CHANNEL *cnet, TOPO *top, PAR *par, long ns, long ng, double zmu, double zmT,
                              double z0s, double d0s,
                              double rz0s, double z0v, double d0v, double rz0v, double hveg, double v,
                              double Ta, double Qa, double P, double LR, double eps, double fc, double LSAI,
                              double decaycoeff0, double *Wcrn, double Wcrnmax, double *Wcsn,
                              double Wcsnmax, double SWin, double LWin, double SWv, double *LW, double *H,
                              double *E,
                              double *LWv, double *Hv, double *LEv, double *Etrans, double *Ts, double *Qs,
                              double Eadd, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Qv,
                              double *Qg,
                              double *Lob, double *rh, double *rv, double *rb, double *rc, double *ruc,
                              double *u_top, double *decay, double *Locc, double *LWup_ab_v, long *lpb,
                              double *dUsl);

void update_soil_land(long nsurf, long n, long i, long r, long c, double fc,
                      double Dt, ENERGY *egy, double **pa, SOIL_STATE *S, DOUBLETENSOR *ET,
                      DOUBLEMATRIX *th);

void update_soil_channel(long nsurf, long n, long ch, double fc, double Dt,
                         ENERGY *egy, double **pa, SOIL_STATE *S, DOUBLEMATRIX *ET, DOUBLEMATRIX *th);

void update_F_energy(long nbeg, long nend, DOUBLEVECTOR *F, double w,
                     DOUBLEVECTOR *K, double *T);

void update_diag_dF_energy(long nbeg, long nend, DOUBLEVECTOR *dF, double w,
                           DOUBLEVECTOR *K);

double calc_C(long l, long nsng, double a, double *wi, double *wl, double *dw,
              double *D, double **pa);

void EnergyFluxes(double t, double Tg, long r, long c, long n, double Tg0,
                  double Qg0, double Tv0, double zmu, double zmT, double z0s,
                  double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg,
                  double v, double Ta, double Qa, double P, double LR,
                  double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn,
                  double Wcrnmax, double Wcsn, double Wcsnmax,
                  double *dWcrn, double *dWcsn, double *theta, double **soil, double *land,
                  double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer,
                  double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT,
                  double *E, double *dE_dT, double *LWv, double *Hv,
                  double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs,
                  double *Hg0, double *Hg1, double *Eg0,
                  double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb,
                  double *ruc, double *rh_g,
                  double *rv_g, double *Qg, double *u_top, double *decay, double *Locc,
                  double *LWup_above_v, double *T,
                  DOUBLEVECTOR *soil_evap_layer_bare, DOUBLEVECTOR *soil_evap_layer_veg,
                  double sky);

void EnergyFluxes_no_rec_turbulence(double t, double Tg, long r, long c,
                                    long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT,
                                    double z0s,
                                    double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg,
                                    double v, double Ta, double Qa, double P, double LR,
                                    double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn,
                                    double Wcrnmax, double Wcsn, double Wcsnmax,
                                    double *dWcrn, double *dWcsn, double *theta, double **soil, double *land,
                                    double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer,
                                    double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT,
                                    double *E, double *dE_dT, double *LWv, double *Hv,
                                    double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs,
                                    double *Hg0, double *Hg1, double *Eg0,
                                    double *Eg1, double *Lobukhov, double *rh, double *rv, double *rc, double *rb,
                                    double *ruc, double *rh_g,
                                    double *rv_g, double *Qg, double *u_top, double *decay, double *Locc,
                                    double *LWup_above_v, double *T,
                                    DOUBLEVECTOR *soil_evap_layer_bare, DOUBLEVECTOR *soil_evap_layer_veg,
                                    double sky, short flagTmin, long cont);

double k_thermal(short snow, short a, double th_liq, double th_ice,
                 double th_sat, double k_solid);

double flux(long i, long icol, double **met, double k, double add,
            double est);

void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR *ad,
                  DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e, double *T,
                  SHORTVECTOR *mf);

double soil_red_evap(double psi, double T);

double red_evap(long n, double psi, double T);

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD,
                           double thres, double z0snow, double *z0_ris, double *d0_ris,
                           double *z0_z0t_ris);

void check_continuity(long n, double *dw, double *wi, double *wl);

void merge(double a, DOUBLEVECTOR *ice, DOUBLEVECTOR *liq, DOUBLEVECTOR *Temp,
           DOUBLEVECTOR *D, long l, long lup, long ldw, long tot);

void sux_minus6_condition(double ic, double wa, double rho, double D1,
                          ENERGY *E);
