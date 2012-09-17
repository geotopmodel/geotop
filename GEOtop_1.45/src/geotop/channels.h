
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225-9 'Moab' - 24 Aug 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225-9 'Moab'
 
 GEOtop 1.225-9 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225-9 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef CHANNELS_H
#define CHANNELS_H
#include "struct.geotop.h"
#include "constants.h"
extern T_INIT *UV;


void enumerate_channels(CHANNEL *cnet, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, DOUBLEMATRIX *Z, DOUBLEMATRIX *slope, long novalue);
void next_down_channel_pixel( long r, long c, double **Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long novalue, long *R, long *C);
void find_max_constraint( double **Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long novalue, long *R, long *C);
short neighboring_down_channel_pixel( long r, long c, long ir, long ic, double **Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long novalue);


#endif
