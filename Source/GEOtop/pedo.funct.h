
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

double psi_teta(double w, double i, double s, double r, double a, double n, double m, double pmin, double st);
double teta_psi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double st);
double dteta_dpsi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double st);

double T_max_dteta(double a, double n, double m);
double P_max_dteta(double a, double n, double m);
double d2theta(double P, double n, double m);
double d3theta(double P, double n, double m);

double K(double psi, double K_sat, double imp, double i, double s, double r, double a, double n, double m, double v, double pmin, double T);

double psi_saturation(double i, double s, double r, double a, double n, double m);

double Harmonic_Mean(double D1, double D2, double K1, double K2);
double Arithmetic_Mean(double D1, double D2, double K1, double K2);
double Mean(short a, double D1, double D2, double K1, double K2);

double Psif(double T);

double theta_from_psi(double psi, long l, long r, long c, SOIL *sl, double pmin);
double psi_from_theta(double th, long l, long r, long c, SOIL *sl, double pmin);
double dtheta_dpsi_from_psi(double psi, long l, long r, long c, SOIL *sl, double pmin);
double k_from_psi(long jK, double psi, long l, long r, long c, SOIL *sl, double imp);
double psisat_from(long l, long r, long c, SOIL *sl);

