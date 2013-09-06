
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef TURBULENCE_H
#define TURBULENCE_H
#include "constants.h"
#include "struct.geotop.h"
#include "meteo.h"
#include <vector>

//extern T_INIT *UV;
extern TInit *UV;
//extern char *WORKING_DIRECTORY;
extern std::string WORKING_DIRECTORY;

//extern char *logfile;
extern std::string logfile;

extern long Nl, Nr, Nc;
//extern char *FailedRunFile;
extern std::string FailedRunFile;

class Turbulence {

	public:
		static void aero_resistance(double zmu, double zmt, double z0, double d0, double z0_z0t, double v, double Ta, 
						 double T, double Qa, double Q, double P, double gmT, double *Lobukhov, double *rm, 
						 double *rh, double *rv, short state_turb, short MO, long maxiter);
	
		static void turbulent_fluxes(const double& rh, const double& rv, const double& P, const double& Ta, const double& T, const double& Qa,
							    const double& Q, const double& dQdT, double& H, double& dHdT, double& E, double& dEdT);

		static double Psim(const double& z);

		static double Psih(const double& z);

		static double Zero(const double&);

		static double PsiStab(const double& z);

		static void Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, 
				 double v, double *rm, double *rh, double *rv, GeoVector<double>& w);

		static double cz(double zmeas, double z0, double d0, double L, double (* unstab)(const double& z), double (* stab)(const double& z));

		static double CZ(short state, double zmeas, double z0, double d0, double L, double (*Psi)(const double& z));

		static void Businger(short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, 
						 double DQ, double z0_z0t, double *rm, double *rh, double *rv, double *Lobukhov, long maxiter);

		static double Levap(const double& T);

		static double latent(const double& Ts, const double& Le);

		static void find_actual_evaporation_parameters(long R, long C, double *alpha, double *beta, GeoVector<double>& evap_layer, double *theta,
											  //double **soil, double *T, double psi, double P, double rv, double Ta, double Qa, 
											  GeoTensor<double>& soil, long sy, double *T, double psi, double P, double rv, double Ta, double Qa, 
											  double Qgsat, long nsnow);

	private:
		static void Star(short a, double zmeas, double z0, double d0, double L, double u, double delta, 
					  double M, double N, double R, double *var, double *c, double *z0v, 
					  double (*Psi)(const double& z), double (*roughness)(const double& x, const double& y, const double& z));
	
		static double roughT(const double& M, const double& N, const double& R);

		static double roughQ(const double& M, const double& N, const double& R);
};
#endif
