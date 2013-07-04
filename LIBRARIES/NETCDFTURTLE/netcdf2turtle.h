
/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file netcdf2turtle.h

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico and Riccardo Rigon

This file is part of numerioc_solver.
	 Turtle_NetCdf is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

     Turtle_NetCdf is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 *
 * \file netcdf2turtle.h
 * \data 15 September 2009
 * \author Emanuele Cordano
 */


DOUBLEVECTOR *nc_t_new_doublevector(const char *filename,const char *varname, const char *dimension);
DOUBLEMATRIX *nc_t_new_doublematrix(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y);
FLOATVECTOR *nc_t_new_floatvector(const char *filename,const char *varname, const char *dimension);
INTVECTOR *nc_t_new_intvector(const char *filename,const char *varname, const char *dimension);
LONGVECTOR *nc_t_new_longvector(const char *filename,const char *varname, const char *dimension);
SHORTVECTOR *nc_t_new_shortvector(const char *filename,const char *varname, const char *dimension);
//
FLOATMATRIX *nc_t_new_floatmatrix(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y);
SHORTMATRIX *nc_t_new_shortmatrix(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y);
INTMATRIX *nc_t_new_intmatrix(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y);
LONGMATRIX *nc_t_new_longmatrix(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y);
//
DOUBLETENSOR *nc_t_new_doubletensor(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y,const char *dimension_z);
//
int nc_t_get_doublematrix_level_l(DOUBLEMATRIX *M, long l,const char *filename,const char *dimension_x, const char *dimension_y, const char *dimension_l);
//SCALAR TYPES
signed char nc_get_byte(const char *filename,const char *varname);
double nc_get_double(const char *filename,const char *varname);
float nc_get_float(const char *filename,const char *varname);
short nc_get_short(const char *filename,const char *varname);
int nc_get_int(const char *filename,const char *varname);
long nc_get_long(const char *filename,const char *varname);
//
int nc_t_get_doublevector_level_l(DOUBLEVECTOR *v, long l,const char *filename,const char *dimension_x, const char *dimension_l);
int nc_t_get_doubletensor_level_l(DOUBLETENSOR *t, long l,const char *filename,const char *dimension_x, const char *dimension_y, const char *dimension_z, const char *dimension_l);
//131109_s
DOUBLETENSOR *nc_t_new_doubletensor_rotate180y(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y,const char *dimension_z);
DOUBLEMATRIX *nc_t_new_doublematrix_rotate180y(const char *filename,const char *varname, const char *dimension_x,const char *dimension_y);
int nc_t_get_doubletensor_rotate180y_level_l(DOUBLETENSOR *t, long l,const char *filename,const char *dimension_x, const char *dimension_y, const char *dimension_z, const char *dimension_l);
int nc_t_get_doublematrix_rotate180y_level_l(DOUBLEMATRIX *M, long l,const char *filename,const char *dimension_x, const char *dimension_y, const char *dimension_l);
//131109_e
void nc_get_global_attr_lat_lon_min_max(const char *filename,double *long_min,double *long_max,double *lat_min,double *lat_max);
void nc_get_global_attr_resolution(const char *filename,double *resolution);
void nc_get_global_attr_missing_value(const char *filename,double *missing_value);
void nc_get_global_value(const char *filename,char *attr_name,double *attr_value);
int nc_t_exist_doublevector(const char *filename,const char *varname, const char *dimension);
