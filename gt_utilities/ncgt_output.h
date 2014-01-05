#ifdef USE_NETCDF

#ifndef NCGT_OUTPUT_H
#define NCGT_OUTPUT_H

#include "../libraries/fluidturtle/turtle.h"
//#include "t_utilities.h"
//#include "t_datamanipulation.h"
//#include "init.h"
#include "../libraries/fluidturtle/t_alloc.h"
#include "../libraries/fluidturtle/t_io.h"
#include "../libraries/fluidturtle/tensors3D.h"
#include "ncgt_turtle2netcdf.h"
#include "gt_symbols.h"

long ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
		const char* dimension_y, long counter, short update, short rotate_y, double geotop::input::gDoubleNoValue, GeoMatrix<long>* rc, long **j, long total_pixel);

int ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue);

//long ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
//		const char* dimension_y, long counter, short reinitialize, short update, short rotate_y, double geotop::input::gDoubleNoValue, LONGMATRIX *rc, long **j_cont, long total_pixel);
long ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
		const char* dimension_y, long counter, short reinitialize, short update, short rotate_y, double geotop::input::gDoubleNoValue, GeoMatrix<long>* rc, long **j_cont, long total_pixel);


int ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue);

void* ncgt_new_output_var(const void * m0, const short& nlimdim, const double& novalue, const std::string& suffix, const double& print_flag);

#endif
#endif
