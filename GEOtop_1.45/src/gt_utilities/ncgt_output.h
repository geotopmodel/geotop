#ifdef USE_NETCDF

#ifndef NCGT_OUTPUT_H
#define NCGT_OUTPUT_H

#include "../libraries/fluidturtle/turtle.h"
#include "../libraries/fluidturtle/t_utilities.h"
#include "../libraries/fluidturtle/t_datamanipulation.h"
#include "../libraries/ascii/init.h"
#include "../libraries/fluidturtle/tensor3D.h"
#include "ncgt_turtle2netcdf.h"
#include "gt_symbols.h"


//long ncgt_add_output_var(int ncid, DOUBLEMATRIX *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,const char* dimension_y, long counter, short upgrade, short rotate_y, double number_novalue,DOUBLEMATRIX *rc);
long ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,const char* dimension_y, long counter, short update, short rotate_y, double number_novalue, LONGMATRIX *rc);

int ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue);

long ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
		const char* dimension_y, long counter, short reinitialize, short update, short rotate_y, double number_novalue, LONGMATRIX* rc);

int ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue);



void* ncgt_new_output_var(void * m0, short nlimdim, double novalue, char* suffix, double print_flag);

#endif
#endif
