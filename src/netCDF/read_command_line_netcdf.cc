
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

#include "config.h"

#ifdef USE_NETCDF
#include "read_command_line_netcdf.h"
using namespace std;



int  ncgt_open_from_option_string(int argc,char *argv[], char *option_f,short define_mode,short print) // ncge_open_from_option_string
//int  ncgt_open_from_option_string(int argc,char* argv[],const std::string & option_f,short define_mode,short print)
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
	//int i,s;
	int status;
	char *function_name="ncgt_open_from_option_string";
	//string function_name="ncgt_open_from_option_string";
	char *filename=read_option_string(argc,argv,option_f,NC_GEOTOP_NULL_EXIT,GEOT_VERBOSE);
	//string filename=read_option_string(argc,argv,option_f, NC_GEOTOP_NULL_EXIT,GEOT_VERBOSE);
				  //read_option_string(int argc,std::string argv[], std::string& option_f,std::string& no_option_argument,short print)

	int ncid=NC_GEOTOP_MISSING;

	if (strcmp(filename,NC_GEOTOP_NULL_EXIT)) {
	//if (strcmp(filename.c_str(),NC_GEOTOP_NULL_EXIT)) {


		status=nc_open(filename,NC_WRITE|NC_SHARE,&ncid);
		//status=nc_open(filename.c_str(),NC_WRITE|NC_SHARE,&ncid);
		if(status==NC_NOERR) {
			status=nc_redef(ncid);
			if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_redef");
		} else {
			status=nc_create(filename,NC_GEOTOP_NEW_EMPTY_FILE, &ncid);
			//status=nc_create(filename.c_str(),NC_GEOTOP_NEW_EMPTY_FILE, &ncid);
			if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_create");


/* This line are to test the first function of libCf just installed:
 * this functions add the convention CF 1.0 used as global attribute
 */
//			status=nccf_def_convention(ncid);
//			if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nccf_def_convention");

		}

		if (define_mode==NC_GEOTOP_DEFINE) {
			status=nc_enddef(ncid);
			if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_enddef");
		}
	}


	return ncid;
}


int ncgt_close_geotop_archive(int ncid)  {

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
	char *function_name="ncgt_close_geotop_archive";
	//string function_name="ncgt_close_geotop_archive";

	if (ncid!=NC_GEOTOP_MISSING) {
		status=nc_close(ncid);
		if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_close");

	}

	return ncid;
}

#endif


