
/* STATEMENT:

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)

 Copyright (c), 2016 - GEOtop Foundation

 This file is part of GEOtop 2.1

 GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

 GEOtop 2.1  is distributed as a free software in the hope to create and support
 a community of developers and users that constructively interact. If you just
 use the code, please give feedback to the authors and the community. Any way
 you use the model, may be the most trivial one, is significantly helpful for
 the future development of the GEOtop model. Any feedback will be highly
 appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */
#ifndef PEDO_FUNCT_H
#define PEDO_FUNCT_H

#include "struct.geotop.h"
#include "constants.h"
#include <vector>
#include <string>

double psi_teta(double w,
                double i,
                double s,
                double r,
                double a,
                double n,
                double m,
                double pmin,
                double st);
double teta_psi(double psi,
                double i,
                double s,
                double r,
                double a,
                double n,
                double m,
                double pmin,
                double st);
double dteta_dpsi(double psi,
                  double i,
                  double s,
                  double r,
                  double a,
                  double n,
                  double m,
                  double pmin,
                  double st);

double k_hydr_soil(double psi,
                   double ksat,
                   double imp,
                   double i,
                   double s,
                   double r,
                   double a,
                   double n,
                   double m,
                   double v,
                   double T,
                   double ratio);

double psi_saturation(double i,
                      double s,
                      double r,
                      double a,
                      double n,
                      double m);

double Harmonic_Mean(double D1, double D2, double K1, double K2);
double Arithmetic_Mean(double D1, double D2, double K1, double K2);
double Mean(short a, double D1, double D2, double K1, double K2);

double Psif(double T);

// double theta_from_psi(double psi, double ice, long l, double **pa, double
// pmin);
double theta_from_psi(double psi,
                      double ice,
                      long l,
                      const GeoTensor<double> &pa,
                      const long &sy,
                      double pmin);
double psi_from_theta(double th,
                      double ice,
                      long l,
                      const GeoTensor<double> &pa,
                      const long &sy,
                      double pmin);
// double dtheta_dpsi_from_psi(double psi, double ice, long l, double **pa,
// double pmin);
double dtheta_dpsi_from_psi(double psi,
                            double ice,
                            long l,
                            GeoTensor<double> &pa,
                            const long &sy,
                            double pmin);
// double k_from_psi(long jK, double psi, double ice, double T, long l, double
// **pa, double imp, double ratio);
double k_from_psi(long jK,
                  double psi,
                  double ice,
                  double T,
                  long l,
                  GeoTensor<double> &pa,
                  const long &sy,
                  double imp,
                  double ratio);
// double psisat_from(double ice, long l, double **pa);
double psisat_from(double ice, long l, GeoTensor<double> &pa, const long &sy);

#endif
