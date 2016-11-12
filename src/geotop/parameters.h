
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1 
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "struct.geotop.h"
//#include "input.h"
#include "constants.h"
#include "../libraries/ascii/tabs.h"
//#include "../libraries/ascii/extensions.h"
#include "../libraries/ascii/rw_maps.h"
#include "times.h"
#include "pedo.funct.h"
#include "input.h"
#define __MATHOPTIM_H__
#include <meteoio/MeteoIO.h>
#include <string>

short read_inpts_par(Par *par, Land *land, Times *times, Soil *sl, Meteo *met, InitTools *itools, FILE *flog) ;

short read_soil_parameters(std::string name, InitTools *IT, Soil *sl, long bed, FILE *flog);

short read_point_file(std::string name, std::vector<std::string> key_header, Par *par, FILE *flog);

short read_meteostations_file(const GeoVector<long>& i, MeteoStations *S, std::string name, std::vector<std::string> key_header, FILE *flog);

#endif
