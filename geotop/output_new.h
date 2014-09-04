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

/**
 *@file output_new.h
 *@author Gianfranco Gallizia
 *@date 07/2014
 *@brief Definitions for the new output post-processing
 */

#ifndef OUTPUT_NEW_H
#define OUTPUT_NEW_H

#include "struct.geotop.h"

#ifdef METEOIO_OUTPUT
void output_file_preproc(AllData* A, mio::Config& mioConfig);
#else
void output_file_preproc(AllData* A);
#endif

void write_output_new(AllData* A);

void deallocate_output_new();

#endif

