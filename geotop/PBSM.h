/*
 
 PRAIRIE BLOWING SNOW MODEL CODE
 
 Code written by Stefano Endrizzi by translating and adapting the idea behind the PBSM Fortran Code 
 by John Pomeroy. The author does not guarantee the perfect conformance of this code with the Fortran one. 
 However, he asks to give credit to Pomeroy, 1993, when using it with satisfaction. 
 
 Reference:
 John Pomeroy
 The Prairie Blowing Snow Model: characteristics, validation, operation
 Journal of Hydrology, 144 (1993) 165-192
 */
#ifndef PBSM_H
#define PBSM_H

#include "constants.h"
#include "struct.geotop.h"
//#include "../libraries/math/util_math.h"

//extern T_INIT *UV;
extern TInit *UV;
extern long Nl, Nr, Nc;
//extern char *WORKING_DIRECTORY;
extern std::string WORKING_DIRECTORY;

//global variables
extern double T, Ustar, Z0, Zr;
extern double k_atm, Diff, B, SvDens, Sigma;
extern float CC1, CC2, CC3, CC4, CC5;

void Pbsm (long r, long c, double Fetch, double N, double dv, double Hv, double rho_sn, double zmeas, double V, double Ta, double RH, 
		   double *Trans, double *Subl, double *Salt, double Dsnow, double slope, FILE *flog);

double suspension(double Z);

double sublimation(double Z);

#endif
