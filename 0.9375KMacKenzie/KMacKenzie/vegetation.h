
void Tcanopy(long r, long c, double Tv0, double Tg, double Qg, double dQgdT, double Tg0, double Qg0, double Ta, double Qa, double zmu, double zmT, 
	double z0, double z0s, double d0, double z0r, double hveg, double v, double LR, double P, double SW, double SWv, double LW, double e, double LAI, 
	double *land, double Wcrn0, double Wcrnmax, double Wcsn0, double Wcsnmax, double *dWcrn, double *dWcsn, double *LWv, double *LWg, double *Hv, 
	double *Hg, double *dHgdT, double *LEv, double *Eg, double *dEgdT, double *Ts, double *Qs, double *froot, double *theta, double *ftl, 
	DOUBLEVECTOR *rep, PAR *par, long n, double sat, double *rh, double *rv, double *rc, double *rb, double *rh_ic, double *rv_ic, double *u_top, 
	double *Etrans, double *Tv, double *Qv, double *decay, double *Locc);			 
					
void canopy_fluxes(long r, long c, double Tv, double Tg, double Ta, double Qg, double Qa, double zmu, double zmT, double z0, double z0s, 
	double d0, double z0r, double hveg, double v, double LR, double P, double SW, double LW, double e, double LAI, double *land, double Wcrn, 
	double Wcrnmax, double Wcsn, double Wcsnmax, double *Esubl, double *Etrans, double *LWv, double *LWg, double *H, double *LE, double *h, 
	double *dhdT, double *Ts, double *Qs, double *rh_ic, double *rv_ic, double *froot, double *theta, double *ftl, long chgsgn, double *Lobukhov, 
	DOUBLEVECTOR *rep, PAR *par, long n, double sat, double *rh, double *rv, double *rc, double *rb, double *u_top, double *decay, double *Locc);
				
void shortwave_vegetation(double Sd, double Sb, double x, double fwsn, double wsn, double Bsnd, double Bsnb, double Agd, double Agb, double C, double R, double T, double L, 
	double *Sv, double *Sg);
	
void longwave_vegetation(double Lin, double eg, double Tg, double Tv, double L, double *Lv, double *Lg, double *dLv_dTv);

void canopy_rain_interception(double rain_max_loading, double LAI, double Prain, double *max_storage, double *storage, double *drip);

void canopy_snow_interception(double snow_max_loading, double LAI, double Psnow, double Tc, double v, double Dt, double *max_storage, double *storage, double *drip);

void update_roughness_veg(double hc, double snowD, double zmu, double zmt, double *z0_ris, double *d0_ris, double *hc_ris);

void root(double d, double *D, double *root_fraction);

void canopy_evapotranspiration(double rbv, double Tv, double Qa, double Pa, double SWin, double *theta, double *land, double *root, double *f, double *fl);

void veg_transmittance(long r, long c, double rm, double v, double Ts, double Tg, double z0soil, double LAI, double z0veg, double d0veg, 
	double Hveg, double u_top, double Lo, double *rb, double *rh, double *ra, double *decay, double Loc);