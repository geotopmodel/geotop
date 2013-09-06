
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




#ifdef USE_NETCDF
#ifndef GT_SYMBOLS_H
#define GT_SYMBOSL_H

#define GEOT_VERBOSE 0
#define GEOT_DEAFAULT -1234


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
	#define NC_GEOTOP_YLAT "YLAT" // North coordinate (m)
	#define NC_GEOTOP_XLON "XLON" // East coordinate (m)
	#define NC_GEOTOP_Z_GENERIC "SOIL_DEPTH"
	#define NC_GEOTOP_Z_SOIL "SOIL_DEPTH"
	#define NC_GEOTOP_Z_SNOW "SNOW_LAYER"
	#define NC_GEOTOP_TIME_GENERIC "time"
	#define NC_GEOTOP_TIME_FOR_SOIL "time_for_soil"
	#define NC_GEOTOP_TIME_FOR_SNOW "time_for_snow"
	#define NC_GEOTOP_TIME_FOR_GLAC "time_for_glac"
	#define NC_GEOTOP_TIME_FOR_POINT_DATA "time_for_points"
	#define NC_GEOTOP_TIME_FOR_SURFACE_ENERGY "time_for_surface_energy"
	#define NC_GEOTOP_POINT_DIM_GENERIC "ID_control_point"
	#define NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT "_in_control_points"
	#define NC_GEOTOP_MISSING_DIMENSION "missing_dimension"
	#define NC_GEOTOP_TIMEPREFIX "TIME_FOR_" // prefix for the temporal axis ()
	#define NC_GEOTOP_MISSING_VALUE_ATTRIBUTE "missing_value"
	#define NC_GEOTOP_UPDATE_COUNTER_TIME 1
	#define NC_GEOTOP_NOUPDATE_COUNTER_TIME 0
	#define NC_GEOTOP_REINITIALIZE_VARIABLE 1
	#define NC_GEOTOP_NOREINITIALIZE_VARIABLE 0
	#define NC_GEOTOP_NOVALUE -9999
	#define NC_GEOTOP_ROTATE_Y 1
	#define NC_GEOTOP_NOROTATE_Y 0
	#define NC_GEOTOP_2D_MAP 2
	#define NC_GEOTOP_3D_MAP 3
	#define NC_GEOTOP_2D_MAP_IN_CONTROL_POINT 52
	#define NC_GEOTOP_3D_MAP_IN_CONTROL_POINT 53
	#define NC_GEOTOP_POINT_VAR 1
	#define NC_GEOTOP_Z_POINT_VAR 22
	#define NC_GEOTOP_0DIM_VAR 0
	#define NC_GEOTOP_UNSTRUCT_MAP 42
	#define NC_GEOTOP_Z_UNSTRUCT_MAP 43
	#define NC_GEOTOP_UNSTRUCT_MAP_IN_CONTROL_POINT 72
	#define NC_GEOTOP_Z_UNSTRUCT_MAP_IN_CONTROL_POINT 73
#endif
#endif
