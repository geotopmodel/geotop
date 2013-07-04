/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file inpts2netcdf.c

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

    This files contains the function to load from the inpts.cf geotop input file, the details about
    a specific netcdf variable:

*/
/*
 * inpts2netcdf.c
 *
 *  Created on: 17-nov-2009
 *      Author: everri
 */
/*i.e sample .inpts.cf line identifier to use in the netcdf version of .inpts file :
structure is
#VAR_NAME#NUM_ATTR#NAME_ATTR1#VAL_ATTR2#....#NAME_ATTR_N#VAL_ATTR_N#
sample :
#height_of_terrain_above_reference_ellipsoid#2##axis#x#units#meters#
*/
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "inpts2netcdf.h"

#define END_VAL_SEPARATOR "#" 	//this separator (token) is used at the start and end of any values
#define SUFFIX_VARIABLE_NAME_SEPARATOR "$" 	//this separator (token) is used at the start of a suffix variable name
//type variable is managed at runtime
#define VAR_NAME_FIELD 1  //the field number containing the variable name (in the example is "height_of_terrain_above_reference_ellipsoid")
#define ATTR_NUM_FIELD 2 //the field number containing the number of attributes  (in the example is "2")
#define ATTR_NAME_FIELD 3 //the field number containing the attribute name (in the example are "axis" and "units")
#define ATTR_VAL_FIELD 4 //the field number containing the attribute value (in the example are "x" and "meters")
#define ATTR_FIELD_USED 2 //use only name and value for attributes

char * inp_get_nc_generic_token(char *inputstring,int token_number,int attrindex){
	/*!
	 *\param inputstring - (char *) info line from .inpts file
	 *
	 *\brief This private function find in the inputstring the token number and return corresponding value
	 *
	 * \author Enrico Verri
	 * \date November 2009
	 *
	 * \return variable id or NULL if search fails.

	 *
	 */
	//const char *function_name="inp_get_nc_generic_token";
	char *res;
	char *cp;
	int i=0;

	if (inputstring == NULL)
		return NULL;
	else
	{
		cp = strdup (inputstring); /* Make writable copy.  */
		if (token_number>=VAR_NAME_FIELD){
			res = strtok (cp,END_VAL_SEPARATOR);
			for(i=ATTR_NUM_FIELD;i<=token_number+attrindex;i++){
				res = strtok (NULL,END_VAL_SEPARATOR);
			}
			return res;
		}
		else
			return NULL;
	}
}



char * inp_get_nc_var_name(char *inputstring){
	/*!
	 *\param inputstring - (char *) info line from .inpts file
	 *
	 *\brief This function find in the inputstring the VAR_NAME_FIELD and return corresponding value
	 *
	 * \author Enrico Verri
	 * \date November 2009
	 *
	 * \return variable name or FATAL_ERROR if search fails (i.e if no variable name is specified)
	 * 		   the variable name is mandatory for netcdf standard

	 *
	 */
	//const char *function_name="inp_get_nc_var_name";
	char *res=inp_get_nc_generic_token(inputstring,VAR_NAME_FIELD,0);
	if (res==NULL)
		t_error("Missing netcdf variable name in the *.inpts.cf file   ");
	return res;
}


int inp_get_nc_attr_num(char *inputstring){
	/*!
	 *\param inputstring - (char *) info line from .inpts file
	 *
	 *\brief This function find in the inputstring the ATTR_NUM_FIELD and return corresponding numeric value
	 *
	 * \author Enrico Verri
	 * \date November 2009
	 *
	 * \return number of attribute for specific variable or -1 if search fails.

	 *
	 */
	//const char *function_name="inp_get_nc_attr_num";
	char *res= inp_get_nc_generic_token(inputstring,ATTR_NUM_FIELD,0);
	if (res == NULL ){
		return 0;
		}
	else
		return atoi(res);
}


char * inp_get_nc_attr_name(char *inputstring,int attrindex){
	/*!
	 *\param inputstring - (char *) info line from .inpts file
	 *\param attrindex : use 1 for extract first attribute name, etc ...
	 *
	 *\brief This function find in the inputstring the ATTR_NAME_FIELD and return corresponding numeric value
	 *
	 * \author Enrico Verri
	 * \date November 2009
	 *
	 * \return name of attribute for specific variable or -1 if search fails.

	 *
	 */
	//const char *function_name="inp_get_nc_attr_name";
	return inp_get_nc_generic_token(inputstring,ATTR_NAME_FIELD,(ATTR_FIELD_USED * (attrindex-1)));

}

char * inp_get_nc_attr_value_str(char *inputstring,int attrindex){
	/*!
	 *\param inputstring - (char *) info line from .inpts file
	 *\param attrindex : use 1 for extract first attribute value string , etc ...
	 *
	 *\brief This function find in the inputstring the ATTR_VAL_FIELD and return corresponding string	 *
	 * \author Enrico Verri
	 * \date November 2009
	 *
	 * \return value of attribute for specific variable or NULL if search fails.

	 *
	 */
	//const char *function_name="inp_get_nc_attr_value_str";
	return inp_get_nc_generic_token(inputstring,ATTR_VAL_FIELD,(ATTR_FIELD_USED * (attrindex-1)));
}

char * inp_get_nc_variable_suffix(char *inputstring,int extractmode){
	/*!
	 *\param inputstring - (char *) info line from .inpts file
	 *       extractmode
	 *\brief This  function find in the inputstring the token $ and return corresponding variable suffix value
	 *
	 * \author Enrico Verri
	 * \date November 2009
	 *
	 * \return variable suffix if extractmode = 0
	 * \       netcdf varname + attributes if extractmode =1
	 * \       or NULL if search fails.

	 *
	 */
	char *res;
	char *cp;

	if (inputstring == NULL)
		return NULL;
	else
	{
		cp = strdup (inputstring); /* Make writable copy.  */
		res = strtok (cp,SUFFIX_VARIABLE_NAME_SEPARATOR); //this is variable name + attributes standard line
		if (extractmode == EXTRACT_MODE_NETCDF_VAR_INFO)
			return res;

		res = strtok (NULL,END_VAL_SEPARATOR); //return suffix, if any
		return res;

	}
}

