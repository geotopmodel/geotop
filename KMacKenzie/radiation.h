void shortwave_radiation(long r, long c, double alpha, double direction, double E0, short shadow, double sky, double tau_cloud, double sa, double slope, double aspect,
	double tau_atm, double *met_data, long *met_col, double sky_st, double A, double *SWbeam, double *SWdiff, double *cosinc, LONGMATRIX *nDt_shadow,
	LONGMATRIX *nDt_sun);

void shadow_n(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow);

double diff2glob(double a);

double atm_transmittance(double alpha, double P, double RH, double T, double A, double Vis, double Lozone);

void longwave_radiation(short state, double pvap, double RH, double T, double taucloud, double *eps, double *eps_max, double *eps_min);

double SB(double T);

double dSB_dT(double T);

short shadows_point(double **hor_height, double alpha, double azimuth, double tol_mount, double tol_flat);

void sun(double hour, double JD, double *alpha, double *direction, double *E0, double latitude, double longitude, double standard_time);

void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, SNOW *snow);

void find_tau_cloud_station(long i, METEO *met, PAR *par, double alpha, double E0, double A, double sky, double *tau, double *sa);

void shadow_haiden(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow);

