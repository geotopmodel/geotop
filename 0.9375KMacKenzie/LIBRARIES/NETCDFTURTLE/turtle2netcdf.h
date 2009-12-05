
/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file turtle2netcdf.h

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
 * \file turtle2netcdf.h
 * \data 11 September 2009
 * \author Emanuele Cordano
 */

int t_nc_newemptyfile(const char *filename);
int t_nc_put_doublevector(DOUBLEVECTOR *v,const char *filename,const char *dimension);
int t_nc_put_var_textattributes(const char *filename,const char *varname, const char *attribute_name, const char *attribute_text);
int t_nc_put_doublematrix(DOUBLEMATRIX *m, const char *filename, const char *dimension_x, const char *dimension_y);
int t_nc_put_floatvector(FLOATVECTOR *v, const char *filename, const char *dimension);
int t_nc_put_intvector(INTVECTOR *v, const char *filename, const char *dimension);
int t_nc_put_longvector(LONGVECTOR *v, const char *filename, const char *dimension);
int t_nc_put_shortvector(SHORTVECTOR *v, const char *filename, const char *dimension);
//
int t_nc_put_floatmatrix(FLOATMATRIX *m, const char *filename, const char *dimension_x, const char *dimension_y);
int t_nc_put_shortmatrix(SHORTMATRIX *m, const char *filename, const char *dimension_x, const char *dimension_y);
int t_nc_put_intmatrix(INTMATRIX *m, const char *filename, const char *dimension_x, const char *dimension_y);
int t_nc_put_longmatrix(LONGMATRIX *m, const char *filename, const char *dimension_x, const char *dimension_y);
//
int t_nc_put_doubletensor(DOUBLETENSOR *dt, const char *filename, const char *dimension_x, const char *dimension_y, const char *dimension_z);
//
int t_nc_put_doublematrix_vs_time(DOUBLEMATRIX *m, long k, const char *filename, const char *dimension_t,  const char *dimension_x, const char *dimension_y);
//
int nc_put_byte(signed char bVal, const char *filename,const char *varname, const char *units, const char *description,const char *standard_name,const char *long_name);
int nc_put_double(double dval, const char *filename,const char *varname, const char *units, const char *description,const char *standard_name,const char *long_name);
int nc_put_float(float fval, const char *filename,const char *varname, const char *units, const char *description,const char *standard_name,const char *long_name);
int nc_put_short(short sval, const char *filename,const char *varname, const char *units, const char *description,const char *standard_name,const char *long_name);
int nc_put_int(int ival, const char *filename,const char *varname, const char *units, const char *description,const char *standard_name,const char *long_name);
int nc_put_long(long lval, const char *filename,const char *varname, const char *units, const char *description,const char *standard_name,const char *long_name);
//
int t_nc_put_doublevector_vs_time(DOUBLEVECTOR *v, long k, const char *filename, const char *dimension_t,  const char *dimension_x);
int t_nc_put_doubletensor_vs_time(DOUBLETENSOR *t, long k, const char *filename, const char *dimension_t,  const char *dimension_x, const char *dimension_y, const char *dimension_z);
int t_nc_put_rotate180_y_doublematrix(DOUBLEMATRIX *m, const char *filename, const char *dimension_x, const char *dimension_y);
int t_nc_put_rotate180_y_doubletensor(DOUBLETENSOR *m, const char *filename, const char *dimension_x, const char *dimension_y, const char *dimension_z);
int t_nc_put_rotate180_y_doubletensor_vs_time(DOUBLETENSOR *t, long k, const char *filename, const char *dimension_t,  const char *dimension_x, const char *dimension_y, const char *dimension_z);
//int create_empty_netcdf4_file(const char *filename);
int t_nc_put_rotate180_y_doublematrix_vs_time(DOUBLEMATRIX *m, long k, const char *filename, const char *dimension_t,  const char *dimension_x, const char *dimension_y);
void nc_add_global_attr_lat_lon_min_max(const char *filename,double long_min,double long_max,double lat_min,double lat_max);
void nc_add_variable_attr_missing_value(const char *filename,const char *varname,double missing_value);
void nc_add_global_attr_missing_value(const char *filename,double missing_value);
void nc_add_global_attr_resolution(const char *filename,double map_resolution);
void nc_add_global_attr_double(const char *filename,char *attr_name,double attr_value);
