/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 _________________________

 Note on meteodistr.c - meteodistr.h

 Basic ideas of the routines distributing wind-precipitation-temperature-relative humidity are derived from the Micromet Fortran Code by Liston and Elder.

 Reference:
 Liston, Glen E.; Elder, Kelly
 A meteorological distribution system for high-resolution terrestrial modeling (MicroMet)
 Journal of Hydrometeorology. 7(April): 217-234.

 However this code is significantly different from the above mentioned-code.

 */


void Meteodistr(double dE, double dN, Matrix<double> *E, Matrix<double> *N, Matrix<double> *topo, Matrix<double> *curvature1,
                Matrix<double> *curvature2, Matrix<double> *curvature3, Matrix<double> *curvature4,
                Matrix<double> *terrain_slope, Matrix<double> *slope_az, METEO *met, double slopewtD, double curvewtD,
                double slopewtI, double curvewtI, double windspd_min, double RH_min, double dn, short iobsint,
                long Tcode, long Tdcode, long Vxcode, long Vycode, long VScode, long Pcode, Matrix<double> *Tair_grid,
                Matrix<double> *RH_grid, Matrix<double> *windspd_grid, Matrix<double> *winddir_grid, Matrix<double> *sfc_pressure,
                Matrix<double> *prec_grid, double T_lapse_rate, double Td_lapse_rate, double Prec_lapse_rate,
                double maxfactorP, double minfactorP, short dew, double Train, double Tsnow, double snow_corr_factor,
                double rain_corr_factor);

short get_temperature(double dE, double dN, Matrix<double> *E, Matrix<double> *N, METEO *met, long Tcode,
                      Matrix<double> *Tair_grid, double dn, Matrix<double> *topo, short iobsint, double lapse_rate);

short get_relative_humidity(double dE, double dN, Matrix<double> *E, Matrix<double> *N, METEO *met, long Tdcode,
                            Matrix<double> *RH_grid, Matrix<double> *Tair_grid, double RH_min, double dn,
                            Matrix<double> *topo,
                            short iobsint, double lapse_rate);

void topo_mod_winds(Matrix<double> *winddir_grid, Matrix<double> *windspd_grid,
                    double slopewtD, double curvewtD, double slopewtI, double curvewtI,
                    Matrix<double> *curvature1, Matrix<double> *curvature2, Matrix<double> *curvature3,
                    Matrix<double> *curvature4, Matrix<double> *slope_az,
                    Matrix<double> *terrain_slope, Matrix<double> *topo, double undef);

short get_wind(double dE, double dN, Matrix<double> *E, Matrix<double> *N, METEO *met, long ucode, long vcode,
               long Vscode, Matrix<double> *windspd_grid, Matrix<double> *winddir_grid, Matrix<double> *curvature1,
               Matrix<double> *curvature2, Matrix<double> *curvature3, Matrix<double> *curvature4,
               Matrix<double> *slope_az, Matrix<double> *terrain_slope, double slopewtD, double curvewtD,
               double slopewtI, double curvewtI, double windspd_min, double dn, Matrix<double> *topo, short iobsint);

short get_precipitation(double dE, double dN, Matrix<double> *E,
                        Matrix<double> *N, METEO *met, long Pcode, long Tcode,
                        long Tdcode, Matrix<double> *prec_grid, double dn, Matrix<double> *topo, short iobsint,
                        double lapse_rate,
                        double max, double min, short dew, double Train, double Tsnow,
                        double snow_corr_factor, double rain_corr_factor);

void get_pressure(Matrix<double> *topo, Matrix<double> *sfc_pressure, double undef);

double find_cloudfactor(double Tair, double RH, double Z, double T_lapse_rate,
                        double Td_lapse_rate);

short interpolate_meteo(short flag, double dX, double dY, Matrix<double> *Xpoint, Matrix<double> *Ypoint, // 5
                        Vector<double> *Xst, Vector<double> *Yst, double **value, long metcod,
                        Matrix<double> *grid, // 5
                        double dn0, short iobsint); // 2

void get_dn(long nc, long nr, double deltax, double deltay, long nstns,
            double *dn);

void barnes_oi(short flag, Matrix<double> *xpoint, Matrix<double> *ypoint, // 3
               Vector<double> *xstnall, Vector<double> *ystnall, Vector<double> *xstn, // 3
               Vector<double> *ystn, Vector<double> *var, double dn, // 3
               double undef, Matrix<double> *grid, double **value_station, long metcode); // 4 => TOT = 13 parameters

