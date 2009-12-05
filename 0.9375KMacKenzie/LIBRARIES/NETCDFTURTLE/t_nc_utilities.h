
/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file t_nc_utilities.h

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
 * \file t_nc_utilities.h
 * \data 15 September 2009
 * \author Emanuele Cordano
 */

#define GLOB_ATTR_LONG_MIN "longitude_min"
#define GLOB_ATTR_LONG_MAX "longitude_max"
#define GLOB_ATTR_LAT_MIN "latitude_min"
#define GLOB_ATTR_LAT_MAX "latitude_max"
#define GLOB_ATTR_MISSING_VALUE "missing_value"
#define ATTR_MISSING_VALUE GLOB_ATTR_MISSING_VALUE
#define NaN 9.969209968386869e+36
#define GLOB_ATTR_MAP_RESOLUTION "GridDepth"

char *copy_stringnames(const char *origin);
//061009_s
void nc_add_cf_convention_header(const char *filename,const char *title, const char *institution,const char *source,const char *history,const char *references,const char *comment,const char *conventions);
void nc_add_cf_convention_var_attributes(const char *filename,const char *varname, const char *units,const char *standard_name,const char *long_name,const char *ancillary_variables,const char *fillvalue,
										 const char *missing_value);//,const char *valid_max,const char *valid_min,const char *valid_range);
//061009_e

int rotate180_y_doublematrix(DOUBLEMATRIX *M); // added by Emanuele Cordano on 22 Oct 2009
int rotate180_y_doubletensor(DOUBLETENSOR *M); // added by Emanuele Cordano on 22 Oct 2009


int rotate180_y_shortmatrix(SHORTMATRIX *M); // added by Emanuele Cordano on 26 Oct 2009
int rotate180_y_floatmatrix(FLOATMATRIX *M); // added by Emanuele Cordano on 26 Oct 2009
int rotate180_y_longmatrix(LONGMATRIX *M); // added by Emanuele Cordano on 26 Oct 2009
int rotate180_y_intmatrix(INTMATRIX *M); // added by Emanuele Cordano on 26 Oct 2009
int invert_order_doublevector(DOUBLEVECTOR *v); // added by Emanuele Cordano on 26 Oct 2009
int invert_order_longvector(LONGVECTOR *v);  // added by Emanuele Cordano on 26 Oct 2009
