void aero_resistance(double zmu, double zmt, double z0, double d0, double z0_z0t, double v, double Ta, double T, double Qa, double Q, double P, double gmT, DOUBLEVECTOR *rep, double *rm, double *rh, 
	double *rv, short state_turb, short MO);
	
void aero_resistance2(double zmu, double zmt, double z0, double d0, double z0_z0t, double hveg, double v, double Ta, double T, double Qa, double Q, double P, 
	double gmT, double LAI, DOUBLEVECTOR *rep, double *rm, double *rh, double *rv, double *u_top, double *Lmo, short state_turb, short MO);
	
void turbulent_fluxes(double rh, double rv, double P, double Ta, double T, double Qa, double Q, double dQdT, double *H, double *dHdT, double *E, double *dEdT);

double Psim(double z);

double Psih(double z);

double Zero(double z);

double PsiStab(double z);

void Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, double v, double *rm, double *rh, double *rv, DOUBLEVECTOR *w);

double cz(double zmeas, double z0, double d0, double L, double (* unstab)(double z), double (* stab)(double z));

double CZ(short state, double zmeas, double z0, double d0, double L, double (*Psi)(double z));

void Star(short a, double zmeas, double z0, double d0, double L, double u, double delta, double M, double N, double R, double *var, double *c, double *z0v,
	double (*Psi)(double z), double (*roughness)(double x, double y, double z) );
	
double roughT(double M, double N, double R);

double roughQ(double M, double N, double R);

void Businger(short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w);

void Businger2(long r, long c, short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w);

double Levap(double T);

double latent(double Ts, double Le);


	