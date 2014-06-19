
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.0.0 - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 2.0.0 
 
 GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef BLOWINGSNOW_H
#define BLOWINGSNOW_H

#include "input.h"
#include "constants.h"
#include "struct.geotop.h"
#include "snow.h"
#include "PBSM.h"
#include "meteo.h"
#include "vegetation.h"
#include "energy.balance.h"
#include "meteodata.h"
#include "geotop_common.h"

void windtrans_snow(Snow *snow, Meteo *met, Land *land, Topo *top, Par *par, double t0);

void set_inhomogeneous_fetch(Snow *snow, Meteo *met, Land *land, Par *par, Topo *top, short *yes);

#endif
