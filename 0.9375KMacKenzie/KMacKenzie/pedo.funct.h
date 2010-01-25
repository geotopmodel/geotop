
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie 

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie. 
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
    
    

double psi_teta(double w, double i, double s, double r, double a, double n, double m, double pmin, double E);
double psi_teta2(double w, double i, double s, double r, double a, double n, double m, double pmin, double E);

double teta_psi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double E);

double dteta_dpsi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double E);

double T_max_dteta(double a, double n, double m);

double d2theta(double P, double n, double m);

double d3theta(double P, double n, double m);

double K(double psi, double K_sat, double imp, double i, double s, double r, double a, double n, double m, double v, double pmin, double T);

double psi_saturation(double i, double s, double r, double a, double n, double m);

double Harmonic_Mean(double D1, double D2, double K1, double K2);

double Arithmetic_Mean(double D1, double D2, double K1, double K2);

double Mean(short a, double D1, double D2, double K1, double K2);

double Psif(double T);

double funct_T(double T, double W, double h, double **soil, long l, double psimin);

double dfunct_T(double T, double W, double h, double **soil, long l, double psimin);

double theta_from_psi(double psi, long l, long r, long c, SOIL *sl, double Esoil);
double psi_from_theta(double th, long l, long r, long c, SOIL *sl, double Esoil);
double dtheta_dpsi_from_psi(double psi, long l, long r, long c, SOIL *sl, double Esoil);
double k_from_psi(long jK, double psi, long l, long r, long c, SOIL *sl, double imp);
double k_from_psi2(long jK, double psi, long l, long r, long c, SOIL *sl, double imp);

double psisat_from(long l, long r, long c, SOIL *sl);


