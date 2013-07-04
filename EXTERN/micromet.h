
void Micromet(T_INIT *UV, DOUBLEMATRIX *topo, DOUBLEMATRIX *curvature, DOUBLEMATRIX *terrain_slope, DOUBLEMATRIX *slope_az, METEO *met, 
	double slopewt, double curvewt, double windspd_min, double dn, short ifill, short iobsint, long Tcode, long RHcode, long Vcode, long Vdircode, 
	long Pcode, DOUBLEMATRIX *Tair_grid, DOUBLEMATRIX *RH_grid, DOUBLEMATRIX *windspd_grid, DOUBLEMATRIX *winddir_grid, DOUBLEMATRIX *sfc_pressure,
	DOUBLEMATRIX *prec_grid, double T_lapse_rate, double Td_lapse_rate, double Prec_lapse_rate);
	
void get_temperature(T_INIT *UV, METEO *met, long Tcode, DOUBLEMATRIX *Tair_grid, double dn, DOUBLEMATRIX *topo, short ifill, short iobsint, 
	double T_lapse_rate);
	
void get_relative_humidity(T_INIT *UV, METEO *met, long RHcode, long Tcode, DOUBLEMATRIX *RH_grid, DOUBLEMATRIX *Tair_grid, double dn, 
	DOUBLEMATRIX *topo, short ifill, short iobsint, double Td_lapse_rate);
	
void topo_mod_winds(DOUBLEMATRIX *winddir_grid, DOUBLEMATRIX *windspd_grid, double slopewt, double curvewt, DOUBLEMATRIX *curvature, 
	DOUBLEMATRIX *slope_az, DOUBLEMATRIX *terrain_slope, DOUBLEMATRIX *topo, double undef);
	
void get_wind(T_INIT *UV, METEO *met, long Vcode, long Vdircode, DOUBLEMATRIX *windspd_grid, DOUBLEMATRIX *winddir_grid, DOUBLEMATRIX *curvature, 
	DOUBLEMATRIX *slope_az, DOUBLEMATRIX *terrain_slope, double slopewt, double curvewt, double windspd_min, double dn, DOUBLEMATRIX *topo,
	short ifill, short iobsint);
	
void get_precipitation(T_INIT *UV, METEO *met, long Pcode, DOUBLEMATRIX *prec_grid, double dn, DOUBLEMATRIX *topo, short ifill, short iobsint, 
	double Prec_lapse_rate);
	
void get_pressure(DOUBLEMATRIX *topo, DOUBLEMATRIX *sfc_pressure, double undef);

double find_cloudfactor(double Tair, double RH, double Z, double T_lapse_rate, double Td_lapse_rate);

void interpolate_meteo(T_INIT *UV, METEO_STATIONS *allmstn, double **value, long **metcol, long metcod, DOUBLEMATRIX *grid, double dn0, short ifill, 
	short iobsint);
	
void get_dn(long nc, long nr, double deltax, double deltay, long nstns, double *dn);

void barnes_oi(long nc, long nr, double deltax, double deltay, double xmn, double ymn, long nstns, DOUBLEVECTOR *xstn, DOUBLEVECTOR *ystn, 
	DOUBLEVECTOR *var, double dn, DOUBLEMATRIX *grid, double undef, short ifill);

double Tdew(double T, double RH, double Z);

double RHfromTdew(double T, double Tdew, double Z);

void topo_data(double deltax, double deltay, DOUBLEMATRIX *topo, DOUBLEMATRIX *curvature, DOUBLEMATRIX *terrain_slope, DOUBLEMATRIX *slope_az, 
	double curve_len_scale, double undef);

