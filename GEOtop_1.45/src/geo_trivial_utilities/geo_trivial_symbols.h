
/* GEOTRIVIALUtilities is an under-construction set of  C functions which supports i/o interface
   and other utilities for models like GEOtop
GEOTRIVIALUtilities Version 1.0

file geo_trivial_utilities.h

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





#define GEOT_VERBOSE 0
#define GEOT_DEAFAULT -1234

#ifdef USE_NETCDF

	#define NC_GEOTOP_NULL_EXIT "null_exit"
	#define NC_GEOTOP_1D_OUTPUT_OPTION "-nc-1D-output"
	#define NC_GEOTOP_ARCHIVE_OPTION "-nc-archive"

	#define NC_GEOTOP_MISSING 12340000


	#define ERRCODE 2
	#define NC_GEOTOP_ERROR_MESSAGE(e,n_function,n_ncfunction) {printf("Error in %s() function: %s",n_function,n_ncfunction); printf("\nError: %s\n", nc_strerror(e)); exit(ERRCODE);}
	#define NC_GEOTOP_ERROR_MESSAGE_VERBOSE(filename,varname,e,n_function,n_ncfunction) {printf("Error in %s() function: %s - netcdf_filename: %s - Variable: %s",n_function,n_ncfunction,filename,varname); printf("\nError: %s\n", nc_strerror(e)); exit(ERRCODE);}
	#define NC_GEOTOP_GLOBAL_ATTRIBUTE "global_attribute"
	#define NC_GEOTOP_DEFINE 1
	#define NC_GEOTOP_NODEFINE 0

	#define NC_GEOTOP_NEW_EMPTY_FILE NC_CLOBBER|NC_NETCDF4

	/* 26.11.2011 space-time dimension */
	#define NC_GEOTOP_XLAT "XLAT" // Latitude coordinate (m)
	#define NC_GEOTOP_YLON "YLON" // Longitude coordinate (m)
	#define NC_GEOTOP_TIME_GENERIC "time"
	#define NC_GEOTOP_TIMEPREFIX "TIME_FOR" // prefix for the temporal axis ()
#endif
