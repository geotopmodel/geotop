/**
 * @file output_new.h
 * @author Gianfranco Gallizia
 * @copyright (C) 2014 eXact lab srl
 * @date 07/2014
 * @brief Definitions for the new output post-processing
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

