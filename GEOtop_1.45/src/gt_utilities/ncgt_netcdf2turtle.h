/*
 * gt_netcdf2turtle.h
 *
 *  Created on: Nov 30, 2011
 *      Author: matteo
 */

/*
 * TO BE
 *
 */



#ifdef USE_NETCDF

#ifndef NCGT_NETCDF2TURTLE_H
#define NCGT_NETCDF2TURTLE_H

#define ERRCODE 2
#define ERROR_MESSAGE(e,n_function,n_ncfunction) {printf("Error in %s() function: %s",n_function,n_ncfunction); printf("\nError: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define GLOBAL_ATTRIBUTE "global_attribute"
#define INIT_VALUE -9998

#include "../libraries/fluidturtle/turtle.h"
//#include <netcdf.h>
#include "gt_utilities.h"

DOUBLEVECTOR *ncgt_new_doublevector(int ncid,const char *varname, const char *dimension);

DOUBLEMATRIX *ncgt_new_doublematrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y);
FLOATVECTOR *ncgt_new_floatvector(int ncid,const char *varname, const char *dimension);

INTVECTOR *ncgt_new_intvector(int ncid,const char *varname, const char *dimension);

LONGVECTOR *ncgt_new_longvector(int ncid,const char *varname, const char *dimension);

SHORTVECTOR *ncgt_new_shortvector(int ncid,const char *varname, const char *dimension);

FLOATMATRIX *ncgt_new_floatmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y);

SHORTMATRIX *ncgt_new_shortmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y);

INTMATRIX *ncgt_new_intmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y);

LONGMATRIX *ncgt_new_longmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y);

DOUBLETENSOR *ncgt_new_doubletensor(int ncid,const char *varname, const char *dimension_x,const char *dimension_y,const char *dimension_z);

int ncgt_get_doublematrix_level_l(DOUBLEMATRIX *M, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_l);

//221009_s
int ncgt_get_doublevector_level_l(DOUBLEVECTOR *v, long l,int ncid,const char *dimension_x, const char *dimension_l);

int ncgt_get_doubletensor_level_l(DOUBLETENSOR *t, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_z, const char *dimension_l);

signed char ncgt_get_byte(int ncid,const char *varname);


// NC_INT
int nc_get_int(int ncid,const char *varname);


// NC_SHORT
short nc_get_short(int ncid,const char *varname);
// NC_FLOAT
float nc_get_float(int ncid,const char *varname);

// NC_DOUBLE
double nc_get_double(int ncid,const char *varname);

// NC_LONG
long nc_get_long(int ncid,const char *varname);

DOUBLETENSOR *ncgt_new_doubletensor_rotate180y(int ncid,const char *varname, const char *dimension_x,const char *dimension_y,const char *dimension_z);
DOUBLEMATRIX *ncgt_new_doublematrix_rotate180y(int ncid,const char *varname, const char *dimension_x,const char *dimension_y);
int ncgt_get_doubletensor_rotate180y_level_l(DOUBLETENSOR *t, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_z, const char *dimension_l);
int ncgt_get_doublematrix_rotate180y_level_l(DOUBLEMATRIX *M, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_l);
void nc_get_global_attr_lat_lon_min_max(int ncid,double *long_min,double *long_max,double *lat_min,double *lat_max);

void nc_get_global_attr_missing_value(int ncid,double *missing_value);

void nc_get_global_attr_resolution(int ncid,double *resolution);

void nc_get_global_value(int ncid,char *attr_name,double *attr_value);
#ifdef USE_HPC


#endif
#endif
#endif
