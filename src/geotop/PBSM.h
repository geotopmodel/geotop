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

void Pbsm (long r, long c, double Fetch, double N, double dv, double Hv,
           double rho_sn, double zmeas, double V, double Ta, double RH,
           double *Trans, double *Subl, double *Salt, double Dsnow, double slope,
           FILE *flog);

double suspension(double Z);

double sublimation(double Z);

