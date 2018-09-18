
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



void aero_resistance(double zmu, double zmt, double z0, double d0,
                     double z0_z0t, double v, double Ta, double T, double Qa, double Q, double P,
                     double gmT,
                     double *Lobukhov, double *rm, double *rh, double *rv, short state_turb,
                     short MO, long maxiter);

void turbulent_fluxes(double rh, double rv, double P, double Ta, double T,
                      double Qa, double Q, double dQdT, double *H, double *dHdT, double *E,
                      double *dEdT);

double Psim(double z);

double Psih(double z);

double Zero(double z);

double PsiStab(double z);

void Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t,
           double Ta, double Ts, double v, double *rm, double *rh, double *rv,
           Vector<double>* w);

double cz(double zmeas, double z0, double d0, double L,
          double (* unstab)(double z), double (* stab)(double z));

double CZ(short state, double zmeas, double z0, double d0, double L,
          double (*Psi)(double z));

void Star(short a, double zmeas, double z0, double d0, double L, double u,
          double delta, double M, double N, double R, double *var, double *c,
          double *z0v,
          double (*Psi)(double z), double (*roughness)(double x, double y, double z) );

double roughT(double M, double N, double R);

double roughQ(double M, double N, double R);

void Businger(short a, double zmu, double zmt, double d0, double z0, double v,
              double T, double DT, double DQ, double z0_z0t, double *rm, double *rh,
              double *rv, double *Lobukhov, long maxiter);

double Levap(double T);

double latent(double Ts, double Le);

void find_actual_evaporation_parameters(long R, long C, double *alpha,
                                        double *beta, Vector<double> *evap_layer, double *theta,
                                        MatrixView<double> &&soil, double *T, double psi, double P, double rv,
                                        double Ta,
                                        double Qa, double Qgsat, long nsnow);
