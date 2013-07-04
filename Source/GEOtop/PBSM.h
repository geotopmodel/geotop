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

void BlowingSnow(long r, long c, float FetchDist, SNOW *snow, double T, double RH, double Uten, double StubHt, double *Transport, double *Sublimation, double *LatentH);

void Pbsm (float Meht, float Fht, float Fetch, float Uthr, float Uten, float Temp, float Rel_H, float *Trans, float *Subl, float *LatentH);

void Sum(float TQsalt, float TQsusp, float SBsum, float SBsalt, float *TransFlux, float *SubFlux, float *HeatFlux);

void Probability_Threshold(double Wliq, double Wice, int Snow_Age, double Psnow, double Temperature, float Uten, float *Probability, float *Threshold);

void set_windtrans_snow(SNOW *snow, METEO *met, LAND *land, PAR *par, double t);

void print_windtrans_snow(long r, long c, SNOW *snow, PAR *par);

void extend_topography(DOUBLEMATRIX *M, double novalue);

void extend_topography_row(DOUBLEMATRIX *M, double novalue);

void extend_topography_column(DOUBLEMATRIX *M, double novalue);

void find_the_nearest(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

void find_the_nearest_row(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

void find_the_nearest_column(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc);

short no_novalue(long r, long c, DOUBLEMATRIX *M, double novalue, long *rr, long *cc);

void set_no_value(DOUBLEMATRIX *M, DOUBLEMATRIX *N, double undef);