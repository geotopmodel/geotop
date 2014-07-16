/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 20 Jun 2013
 
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

#ifndef METEODISTR_H
#define METEODISTR_H
#include "constants.h"
#include "struct.geotop.h"
#include "radiation.h"
#include "../libraries/ascii/rw_maps.h"
#include "meteo.h"

void Meteodistr(double dE, double dN, GeoMatrix<double>& E, GeoMatrix<double>& N, GeoMatrix<double>& topo, GeoMatrix<double>& curvature1, GeoMatrix<double>& curvature2,
                GeoMatrix<double>& curvature3, GeoMatrix<double>& curvature4, GeoMatrix<double>& terrain_slope, GeoMatrix<double>& slope_az, Meteo *met,
                double slopewtD, double curvewtD, double slopewtI, double curvewtI, double windspd_min, double RH_min, double dn, short iobsint, 
                long Tcode, long Tdcode, long Vxcode, long Vycode, long VScode, long Pcode, GeoMatrix<double>& Tair_grid, GeoMatrix<double>& RH_grid,
                GeoMatrix<double>& windspd_grid, GeoMatrix<double>& winddir_grid, GeoMatrix<double>& sfc_pressure,GeoMatrix<double>& prec_grid,
                double T_lapse_rate, double Td_lapse_rate, double Prec_lapse_rate, double maxfactorP, double minfactorP, 
                short dew, double Train, double Tsnow, double snow_corr_factor, double rain_corr_factor);

double find_cloudfactor(double Tair, double RH, double Z, double T_lapse_rate, double Td_lapse_rate);

#endif

