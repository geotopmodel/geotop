
/* GEOTRIVIALUtilities is an under-construction set of  C functions which supports i/o interface
   and other utilities for models like GEOtop
GEOTRIVIALUtilities Version 1.0

file read_command_line_netcdf.c

Copyright (c), 2011 Emanuele Cordano

This file is part of GEOTRIVIALUtilities.
 GEOTRIVIALUtilities is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOTRIVIALUtilities is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License v. 3.0 for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "geo_trivial_utilities.h"
#include "geo_trivial_symbols.h"

#define SUCCESS 1
#define NO_SUCCESS 0




int  nc_open_from_option_string(int argc,char *argv[], char *option_f,short define_mode,short print)
{
	/*!
	 * \author Emanuele Cordano
	 * \date October 2011
	 *
	 *\param (int) - argc
	 *\param (char *) - argv
	 *\param (char *) - option_f - a string related to the flag which can appear in the command
	 *\param (short) -  print
	 *
	 *\return Returns the pointer (int)  to the netCF archive whose name is the string followed by the string option_f in the command line
	 * \brief The netCDF archive is created or opened in the DEFINE MODE if define_fleg is set to 1 (NC_GEOTOP_DEFINE)
	 *
	*/
	int i,s;
	int status;
	char *function_name="nc_open_open_option_string";
	char *filename=read_option_string(argc,argv,option_f,NC_GEOTOP_NULL_EXIT,GEOT_VERBOSE);
	int ncid=NC_GEOTOP_MISSING;

	if (strcmp(filename,NC_GEOTOP_NULL_EXIT)) {


		status=nc_open(filename,NC_WRITE|NC_SHARE,&ncid);
		if(status==NC_NOERR) {
			status=nc_redef(ncid);
			if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_redef");
		} else {
			status=nc_create(filename,NC_GEOTOP_NEW_EMPTY_FILE, &ncid);
			if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_create");
		}

		if (define_mode==NC_GEOTOP_DEFINE) {
			status=nc_enddef(ncid);
			if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_enddef");
		}
	}


	return ncid;
}


int nc_close_geotop_archive(int ncid)  {

	/*!
	 *
	 *\param ncid netCDF ID
	 *
	 * \date Octobr 2011
	 *
	 * \author Emanuele Cordano
	 *
	 *
	 *\brief It closes a netCDF archive file.
	 *
	 */

	int status;
	char *function_name="nc_close_geotop_archive";

	if (ncid!=NC_GEOTOP_MISSING) {
		status=nc_close(ncid);
		if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_close");

	}

	return ncid;
}




