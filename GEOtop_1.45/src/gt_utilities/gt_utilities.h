
/* GEOTRIVIALUtilities is an under-construction set of  C functions which supports i/o interface
   and other utilities for models like GEOtop
GEOTRIVIALUtilities Version 1.0

file geo_trivial_symbols.h

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




#include "read_command_line.h"




#ifdef USE_NETCDF




//	#include "netcdf.h"
	#include "libcf_src.h"

//	#include <assert.h>
//	#include <assert.h>
	#include "netcdf.h"

//#include "nccf_grid.h"
//#include "nccf_global.h"
//#include "nccf_data.h"
//#include "nccf_mosaic.h"
//#include "nccf_host.h"
//#include "nccf_utility_functions.h"
//#include "nccf_handle_error.h"




//	#include "geo_trivial_netcdf_utilities.h"
	#include "read_command_line_netcdf.h"
	#include "gt_turtle2netcdf.h"
	#include "gt_netcdf2turtle.h"

#endif
