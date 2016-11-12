
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

#ifndef TABLES_H
#define TABLES_H

#include "struct.geotop.h"

//double find_activelayerdepth_up(long i, long ty, SOIL *sl);
double find_activelayerdepth_up(long i, long ty, Soil *sl);

//double find_activelayerdepth_dw(long i, long ty, SOIL *sl);
double find_activelayerdepth_dw(long i, long ty, Soil *sl);

//double find_watertabledepth_up(long i, long ty, SOIL *sl);
double find_watertabledepth_up(long i, long ty, Soil *sl);

//double find_watertabledepth_dw(long i, long ty, SOIL *sl);
double find_watertabledepth_dw(long i, long ty, Soil *sl);

#endif
