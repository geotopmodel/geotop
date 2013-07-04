/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file inpts2netcdf.h

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

    This files contains the function to load from the inpts geotop  file the details about
    a specific netcdf variable:

*/
/*
 * inpts2netcdf.h
 *
 *  Created on: 17-nov-2009
 *      Author: everri
 */

#ifndef INPTS2NETCDF_H_
#define INPTS2NETCDF_H_

#define EXTRACT_MODE_SUFFIX 0
#define EXTRACT_MODE_NETCDF_VAR_INFO 1

char * inp_get_nc_var_name(char *inputstring);
int inp_get_nc_attr_num(char *inputstring);
char * inp_get_nc_attr_name(char *inputstring,int attrindex);
char * inp_get_nc_attr_value_str(char *inputstring,int attrindex);
char * inp_get_nc_variable_suffix(char *inputstring,int extractmode);
#endif /* INPTS2NETCDF_H_ */
