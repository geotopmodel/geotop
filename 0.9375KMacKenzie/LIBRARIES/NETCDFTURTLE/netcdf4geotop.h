/*
 * netcdf4geotop.h
 *
 *  Created on: 24-nov-2009
 *      Author: Enrico Verri - CGS
 */

#ifndef NETCDF4GEOTOP_H_
#define NETCDF4GEOTOP_H_


#endif /* NETCDF4GEOTOP_H_ */

#define LATITUDE_VAR_NAME "lat";
#define LONGITUDE_VAR_NAME "lon";
#define TIME_VAR_NAME "time";
#define TIME_UNLIMITED_PREFIX "time_for_"
#define VERTICAL_DIMENSION "vertical_downward_depth_from_terrain_surface"
#define STANDARD_NO_VALUE -9999.0
#define GEOTOP_NETCDF_MAPS_FILENAME "geotop_maps.nc" //placed in WORKING_DIRECTORY
#define VAL_NETCDF_FORMAT 4 //output map parameter
#define ATTR_X0MAP_NAME "X0map"
#define ATTR_Y0MAP_NAME "Y0map"

DOUBLEMATRIX *read_netcdf(char *filename, DOUBLEMATRIX *Mref, double *Dxmap,double *Dymap, double *X0map, double *Y0map);
void write_doubletensor_in_netcdf(long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV,double time_in_sec);
void write_vertical_dimension_in_netcdf(DOUBLETENSOR *sl_pa,long jdz_val,double par_jd0,long int par_year0);
void write_netcdf(long i, char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV,double time_in_sec);
