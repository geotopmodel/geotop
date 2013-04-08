/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.223 'Wallis' - 26 Jul 2011
 
 Copyright (c), 2011 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.223 'Wallis'
 
 GEOtop 1.223 'Wallis' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.223 'Wallis' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
#ifndef IMPORT_ASCII_H
#define IMPORT_ASCII_H
#include "../fluidturtle/turtle.h"
#include "extensions.h"
#include "../fluidturtle/t_utilities.h"

#define max_figures 30

double *read_grassascii(double *header, double novalue, char *name);

double *read_esriascii(double *header, double novalue, char *name);

void error_message(short format, long n, long n1, long n2, long n3, char *name);
#endif
