/* netcdf4geotop.c CONTAINS FUNCTIONS TO read/write netcdf map in/to geotop Version 0.9375 KMackenzie
to use NETCDFTURLE library in Geotop you must:
1)in Gcccompiler.settings.symbols.defined symbol add USE_NETCDF_MAP
2)add NETCDFTURTLE directory in project.GCC C compiler.settings.directories
****** TO USE A NETCDF4 CONFIGURATION (handle map file > 2GB  + 1 or more variable with unlimted dimension allowed )********
3) netcdf4 lib and include directory (for *.h)
   with all the libraries, netCDF, HDF5, zlib and respective *.h
   The libraries, which must be listed in the CORRECT ORDER: netcdf hdf5_hl hdf5 z
   On some systems (i.e. LINUX) you must also include m for the math library.
   Finally in Gcccompiler.settings.symbols.defined symbol add  USE_NETCDF4
4) if error "NC_NETCDF4  undeclared" appears you  haven't correctly set the lib or path for netcdf4
******IN ALTERNATIVE YOU CAN USE A NETCDF3 CONFIGURATION (handle map file < 2GB  + only one variable with unlimited dimension allowed )********
5) netcdf4 lib and include directory (for *.h)
   Finally in Gcccompiler.settings remove the variable USE_NETCDF4 (if already defined)

file netcdf2turtle.c

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico and Riccardo Rigon

This file is part of NETCDFTURTLE.
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
/*
 * netcdf4geotop.c
 *
 *  Created on: 24-nov-2009
 *      Author: Enrico Verri - CGS
 */
#include "turtle.h"
#include "t_io.h"
#include <netcdf.h>
#include "turtle2netcdf.h"
#include "netcdf2turtle.h"
#include "t_nc_utilities.h"
#include "inpts2netcdf.h"
#include "netcdf4geotop.h"
#include "time.h"

#define NETCDF_CF 1
/*----------    Global GEOTOP variables  ------------*/
T_INIT *UV;
char *MAP_WORKING_DIRECTORY='\0';
/*----------    Global member variables  ------------*/
STRINGBIN *esri_ascii_files;
STRINGBIN *netcdf_variables_and_atttributes;

extern char *WORKING_DIRECTORY;

double simulation_start_time_in_sec = 0.0;
double map_resolution,lat_min,lat_max,long_min,long_max;
long nrow,ncol;
const char * filename_wr=GEOTOP_NETCDF_MAPS_FILENAME;
char * time_units = "seconds since 1970-01-01 00:00:00 UTC";

DOUBLEMATRIX *read_netcdf(char *filename, DOUBLEMATRIX *Mref, double *Dxmap,double *Dymap, double *X0map, double *Y0map){
	char *netcdf_filename=join_strings(WORKING_DIRECTORY,GEOTOP_NETCDF_MAPS_FILENAME);
	DOUBLEMATRIX *dm_currentmap;

	char *current_var_name;
	char *current_var_attr_name;
	char *current_var_attr_val;
	int number_of_attributes,i;
	char *lat_var_name, *long_var_name;
	double lat_min,lat_max,long_min,long_max,map_resolution;

	lat_var_name = LATITUDE_VAR_NAME;
	long_var_name = LONGITUDE_VAR_NAME;

	 if (strncmp(filename,"#",1) == 0){
		 //if current filename indicates a netcdf variable (string start with #)
		 //load map from netcdf variable
		//extract current variable name
		current_var_name=inp_get_nc_var_name(filename);
		number_of_attributes=inp_get_nc_attr_num(filename);
		//printf("ATTR_NAME= %sATTR NUM = %d\n",current_var_name,number_of_attributes);
		//get attributes names+values
		for(i=1;i<=number_of_attributes;i++){
			//extract from current entry
			current_var_attr_name=inp_get_nc_attr_name(filename,i);
			current_var_attr_val=inp_get_nc_attr_value_str(filename,i);
			//printf("ATTR NAME %d:%s\n",i,current_var_attr_name);
			//printf("ATTR_VAL %d:%s\n",i,current_var_attr_val);
		}
		dm_currentmap=nc_t_new_doublematrix_rotate180y(netcdf_filename,current_var_name,long_var_name,lat_var_name);

		nc_get_global_attr_resolution(netcdf_filename,&map_resolution);
		nc_get_global_attr_lat_lon_min_max(netcdf_filename,&long_min,&long_max,&lat_min,&lat_max);

		*Dxmap=map_resolution;
		*Dymap=map_resolution;
		nc_get_global_value(netcdf_filename,ATTR_X0MAP_NAME,X0map);
		nc_get_global_value(netcdf_filename,ATTR_Y0MAP_NAME,Y0map);
	 }
	 return dm_currentmap;
}

void write_doubletensor_in_netcdf(long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV,double time_in_sec){
	char *current_var_name;
	char *current_var_attr_name;
	char *current_var_attr_val;
	char *time_var_name;
	char *lat_var_name, *long_var_name, *layer_var_name;
	int number_of_attributes,idx_of_attribute;
	long idx_time_stamp;
	char *netcdf_filename=join_strings(WORKING_DIRECTORY,GEOTOP_NETCDF_MAPS_FILENAME);
	char *current_var_name_suffix;
	char *cur_filename;

	//if current filename indicates a netcdf variable (string start with #)
	 if (strncmp(filename,"#",1) == 0){
		//extract current variable name
		current_var_name=inp_get_nc_var_name(filename);
		//printf("%s\n",filename);
		//check if filename token expression contains a suffix for variable name
		current_var_name_suffix = inp_get_nc_variable_suffix(filename,EXTRACT_MODE_SUFFIX);
		//printf("%s\n",current_var_name_suffix);
		if (current_var_name_suffix != NULL){
			//ADD SUFFIX
			//printf("SUFFIX %s\n",current_var_name_suffix);
			current_var_name=join_strings(current_var_name,current_var_name_suffix);
			//printf("NEW VAR NAME %s\n",current_var_name);
		}

		cur_filename=join_strings(filename,"");
		//printf("NEW fileNAME %s\n",cur_filename);

		//append TIME_UNLIMITED_PREFIX at variable name to make unlimited time variable name
		time_var_name=join_strings(TIME_UNLIMITED_PREFIX,current_var_name);
		lat_var_name = LATITUDE_VAR_NAME;
		long_var_name = LONGITUDE_VAR_NAME;
		layer_var_name = VERTICAL_DIMENSION;


		//Set here tensor name
		T->name =current_var_name;
		t_nc_put_rotate180_y_doubletensor_vs_time (T, i-1, netcdf_filename,time_var_name,long_var_name,lat_var_name,layer_var_name);

		//insert attributes + missing value for the variable
		number_of_attributes=inp_get_nc_attr_num(cur_filename);
		//printf("ATTR_NAME= %sATTR NUM = %d\n",current_var_name,number_of_attributes);
		//get attributes names+values
		for(idx_of_attribute=1;idx_of_attribute<=number_of_attributes;idx_of_attribute++){
			//extract from current entry and eventually remove suffix
			if (current_var_name_suffix != NULL){
					cur_filename= inp_get_nc_variable_suffix(filename,EXTRACT_MODE_NETCDF_VAR_INFO);
			}
			current_var_attr_name=inp_get_nc_attr_name(cur_filename,idx_of_attribute);
			current_var_attr_val=inp_get_nc_attr_value_str(cur_filename,idx_of_attribute);
			//printf("ATTR NAME %d:%s\n",idx_of_attribute,current_var_attr_name);
			//printf("ATTR_VAL %d:%s\n",idx_of_attribute,current_var_attr_val);
			t_nc_put_var_textattributes(netcdf_filename,current_var_name, current_var_attr_name,current_var_attr_val);

		}
		nc_add_variable_attr_missing_value(netcdf_filename,current_var_name,(double)STANDARD_NO_VALUE);
		//WRITE time variable
		//IMPORTANT! put time vector values after time unlimited definition
		DOUBLEVECTOR *dv_time,*dv_time_new;
		//load dv_time from netcdf
		int time_var_already_exist;
		time_var_already_exist= nc_t_exist_doublevector(netcdf_filename,time_var_name,time_var_name);

		if (time_var_already_exist == 0) {
			//variable time exist read old value and add new value
			dv_time=nc_t_new_doublevector(netcdf_filename,time_var_name,time_var_name);
			dv_time_new=new_doublevector(dv_time->nh +1);
			dv_time_new->name=time_var_name;
			//insert new time stamp
			for(idx_time_stamp=dv_time->nl;idx_time_stamp < dv_time->nh;idx_time_stamp++){//091209
				dv_time_new->element[idx_time_stamp]=dv_time->element[idx_time_stamp];
			}
			dv_time_new->element[idx_time_stamp]=simulation_start_time_in_sec + time_in_sec;
			//printf("WRITE %s in %s->%s\n",dv_time_new->name,filename_wr,time_var_name);
			t_nc_put_doublevector(dv_time_new,netcdf_filename,time_var_name);
			//t_nc_put_var_textattributes(netcdf_filename,dv_time_new->name, "units",time_units);//091209
			free_doublevector(dv_time);
			free_doublevector(dv_time_new);

		}
		else{
			dv_time_new=new_doublevector(1);
			dv_time_new->name=time_var_name;
			dv_time_new->element[1]=simulation_start_time_in_sec + time_in_sec;
			//printf("WRITE %s in %s->%s\n",dv_time_new->name,filename_wr,time_var_name);
			t_nc_put_doublevector(dv_time_new,netcdf_filename,time_var_name);
			t_nc_put_var_textattributes(netcdf_filename,dv_time_new->name, "units",time_units);
			t_nc_put_var_textattributes(netcdf_filename,dv_time_new->name, "axis","t");//091209
			free_doublevector(dv_time_new);
		}


	 }
	 free(cur_filename);
}

void write_netcdf(long i, char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV,double time_in_sec){
	char *current_var_name;
	char *current_var_attr_name;
	char *current_var_attr_val;
	char *time_var_name;
	char *lat_var_name, *long_var_name;
	int number_of_attributes,idx_of_attribute;
	long idx_time_stamp;
	char *netcdf_filename=join_strings(WORKING_DIRECTORY,GEOTOP_NETCDF_MAPS_FILENAME);
	char *current_var_name_suffix;
	char *cur_filename;

	//if current filename indicates a netcdf variable (string start with #)
	 if (strncmp(filename,"#",1) == 0){
		//extract current variable name
		current_var_name=inp_get_nc_var_name(filename);
		//printf("%s\n",filename);
		//check if filename token expression contains a suffix for variable name
		current_var_name_suffix = inp_get_nc_variable_suffix(filename,EXTRACT_MODE_SUFFIX);
		if (current_var_name_suffix != NULL){
			//ADD SUFFIX
			//printf("SUFFIX %s\n",current_var_name_suffix);
			current_var_name=join_strings(current_var_name,current_var_name_suffix);
			//printf("NEW VAR NAME %s\n",current_var_name);
		}

		cur_filename=join_strings(filename,"");
		//printf("NEW fileNAME %s\n",cur_filename);

		//append TIME_UNLIMITED_PREFIX at variable name to make unlimited time variable name
		time_var_name=join_strings(TIME_UNLIMITED_PREFIX,current_var_name);
		lat_var_name = LATITUDE_VAR_NAME;
		long_var_name = LONGITUDE_VAR_NAME;




		//Set here tensor name
		M->name =current_var_name;
		t_nc_put_rotate180_y_doublematrix_vs_time (M, i-1, netcdf_filename,time_var_name,long_var_name,lat_var_name);

		//insert attributes + missing value for the variable
		number_of_attributes=inp_get_nc_attr_num(cur_filename);
		//printf("ATTR_NAME= %sATTR NUM = %d\n",current_var_name,number_of_attributes);
		//get attributes names+values
		for(idx_of_attribute=1;idx_of_attribute<=number_of_attributes;idx_of_attribute++){
			//eventually remove suffix
			if (current_var_name_suffix != NULL){
					cur_filename= inp_get_nc_variable_suffix(filename,EXTRACT_MODE_NETCDF_VAR_INFO);
			}
			//extract from current entry
			current_var_attr_name=inp_get_nc_attr_name(cur_filename,idx_of_attribute);
			current_var_attr_val=inp_get_nc_attr_value_str(cur_filename,idx_of_attribute);
			//printf("ATTR NAME %d:%s\n",idx_of_attribute,current_var_attr_name);
			//printf("ATTR_VAL %d:%s\n",idx_of_attribute,current_var_attr_val);
			t_nc_put_var_textattributes(netcdf_filename,current_var_name, current_var_attr_name,current_var_attr_val);

		}
		nc_add_variable_attr_missing_value(netcdf_filename,current_var_name,(double)STANDARD_NO_VALUE);
		//WRITE time variable
		//IMPORTANT! put time vector values after time unlimited definition
		DOUBLEVECTOR *dv_time,*dv_time_new;
		//load dv_time from netcdf
		int time_var_already_exist;
		time_var_already_exist= nc_t_exist_doublevector(netcdf_filename,time_var_name,time_var_name);

		if (time_var_already_exist == 0) {
			//variable time exist read old value and add new value
			dv_time=nc_t_new_doublevector(netcdf_filename,time_var_name,time_var_name);
			dv_time_new=new_doublevector(dv_time->nh +1);
			dv_time_new->name=time_var_name;
			//insert new time stamp
			for(idx_time_stamp=dv_time->nl;idx_time_stamp < dv_time->nh;idx_time_stamp++){ //091209
				dv_time_new->element[idx_time_stamp]=dv_time->element[idx_time_stamp];
			}
			dv_time_new->element[idx_time_stamp]=simulation_start_time_in_sec + time_in_sec;
			//printf("WRITE %s in %s->%s\n",dv_time_new->name,filename_wr,time_var_name);
			t_nc_put_doublevector(dv_time_new,netcdf_filename,time_var_name);
			//t_nc_put_var_textattributes(netcdf_filename,dv_time_new->name, "units",time_units); //091209
			free_doublevector(dv_time);
			free_doublevector(dv_time_new);

		}
		else{
			dv_time_new=new_doublevector(1);
			dv_time_new->name=time_var_name;
			dv_time_new->element[1]=simulation_start_time_in_sec + time_in_sec;
			//printf("WRITE %s in %s->%s\n",dv_time_new->name,filename_wr,time_var_name);
			t_nc_put_doublevector(dv_time_new,netcdf_filename,time_var_name);
			t_nc_put_var_textattributes(netcdf_filename,dv_time_new->name, "units",time_units);
			t_nc_put_var_textattributes(netcdf_filename,dv_time_new->name, "axis","t");//091209
			free_doublevector(dv_time_new);
		}

	 }
	 free(cur_filename);
}

double convert_gregorian_into_julian(long int Y, long int M, long int D){
	//Convert Gregorian calendar date to Julian Day Number
	//from http://en.wikipedia.org/wiki/Julian_day#cite_note-10
	return  (1461 * (Y + 4800 + (M - 14)/12))/4 +(367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075 ;
}

char * make_time_units(double par_jd0,long int par_year0){
	//return string containing measure units for time
	//calculate simulation start time in sec
	double jd0_1970 ; //julian day corresponding to 1970-01-01 00:00:00
	double diff_jd,jd_startsim;
	long int M,D,Y;
	//
	M=1;
	D=1;
	Y= 1970;
	//julian day corresponding to 1970-01-01 00:00:00
	jd0_1970 = convert_gregorian_into_julian(Y,M,D);

	Y=par_year0;
	//calc jd_startsim referring to year0 M D
	jd_startsim =  convert_gregorian_into_julian(Y,M,D);
	diff_jd=jd_startsim - jd0_1970;

	//convert jd into secods
	simulation_start_time_in_sec=(double)(diff_jd * 100 * 864);

	return "seconds since 1970-01-01 00:00:00 UTC";
}

void write_vertical_dimension_in_netcdf(DOUBLETENSOR *sl_pa,long jdz_val,double par_jd0,long int par_year0){
	//make here the time_units + start time with parameter file values
	time_units = make_time_units(par_jd0,par_year0);

	DOUBLEVECTOR *z_layer_v;
	char *netcdf_filename=join_strings(WORKING_DIRECTORY,GEOTOP_NETCDF_MAPS_FILENAME);

	long f_cnt;

	z_layer_v=new_doublevector(sl_pa->nch);
	double z_cnt=0.0;

	for(f_cnt=z_layer_v->nl;f_cnt<=z_layer_v->nh;f_cnt++) {
		z_layer_v->element[f_cnt]=z_cnt+sl_pa->element[1][jdz_val][f_cnt]/2.0 ; //WARNING : All soil classes must have the same layer thicknesses
		z_cnt=z_cnt+sl_pa->element[1][jdz_val][f_cnt];
	}
	z_layer_v->name=VERTICAL_DIMENSION;
	t_nc_put_doublevector(z_layer_v,netcdf_filename, z_layer_v->name);
	t_nc_put_var_textattributes(netcdf_filename,z_layer_v->name,"units","millimeters");
	t_nc_put_var_textattributes(netcdf_filename,z_layer_v->name, "axis","z");//091209
	free_doublevector(z_layer_v);

}
