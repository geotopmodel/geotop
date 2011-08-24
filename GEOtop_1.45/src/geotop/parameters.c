
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "struct.geotop.h"
#include "input.h"
#include "parameters.h"
#include "constants.h"
#include "../libraries/ascii/tabs.h"
#include "../libraries/ascii/extensions.h"
#include "../libraries/ascii/rw_maps.h"
#include "times.h"

extern long number_novalue, number_absent;
extern char *string_novalue;

extern char *WORKING_DIRECTORY;

extern char **files;

extern long *outputpoint, noutputpoint, *outputbasin, noutputbasin, *outputsnow, noutputsnow;
extern long *outputglac, noutputglac, *outputsoil, noutputsoil;
extern char **headerpoint, **headerbasin, **headersnow, **headerglac, **headersoil;
extern char *keywords_num[num_par_number] , *keywords_char[num_par_char];

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		

short read_inpts_par(PAR *par, LAND *land, TIMES *times, SOIL *sl, METEO *met, INIT_TOOLS *itools, char *filename, FILE *flog){
	
	//variables
	FILE *f;
	
	short res;
	long i, j, inum, istr, n;
	
	long keylength, stringlength, numberlength;
	long *key, *string;
	double *number;
	
	char **keywords_num_read, **keywords_str_read;
	long **string_read;
	long *string_length_read;
	double **number_read;
	long *number_comp_read;
	
	long *num_param_components;
	double **num_param;
	
	char **string_param;
	
	short endoffile, ok;
	
	char *temp;
	
	char **keywords_num_lower_case, **keywords_char_lower_case;
	
	long beg=0, end=0;
	
	//convert keyword listed on top of the file in lower case
	n = (long)num_par_number;
	keywords_num_lower_case = (char**)malloc(n*sizeof(char*));
	for (i=0; i<n; i++) {
		//printf("i:%ld    %s\n",i,keywords_num[i]);
		keywords_num_lower_case[i] = assign_string(keywords_num[i]);
		convert_string_in_lower_case(keywords_num_lower_case[i]);
	}
	
	n = (long)num_par_char;
	keywords_char_lower_case = (char**)malloc(n*sizeof(char*));
	for (i=0; i<n; i++) {
		//printf("i:%ld    %s\n",i,keywords_char[i]);
		keywords_char_lower_case[i] = assign_string(keywords_char[i]);
		convert_string_in_lower_case(keywords_char_lower_case[i]);
	}
	
	//Allocation	
	n = (long)max_charstring;
	key = (long*)malloc(n*sizeof(long));
	string = (long*)malloc(n*sizeof(long));
	
	n = (long)max_numvect;
	number = (double*)malloc(n*sizeof(double));
	
	//read how many lines there are
	inum=0;
	istr=0;
	f = t_fopen(filename, "r");
	do{
		res=readline_par(f, 33, 61, 44, max_charstring, max_numvect, key, &keylength, string, &stringlength, number, &numberlength, &endoffile);
		if(res==1) inum++;
		if(res==2) istr++;
	}while (endoffile==0);
	t_fclose(f);
	
	string_length_read=(long*)malloc(istr*sizeof(long));
	number_comp_read=(long*)malloc(inum*sizeof(long));
	
	keywords_str_read=(char**)malloc(istr*sizeof(char*));
	keywords_num_read=(char**)malloc(inum*sizeof(char*));
	
	string_read=(long**)malloc(istr*sizeof(long*));
	number_read=(double**)malloc(inum*sizeof(double*));
	
	//read single lines
	f = fopen(filename, "r");
	inum=0;
	istr=0;
	do{
		res=readline_par(f, 33, 61, 44, max_charstring, max_numvect, key, &keylength, string, &stringlength, number, &numberlength, &endoffile);
		if(res==1){
			inum++;
			keywords_num_read[inum-1]=find_string(key, keylength);
			convert_string_in_lower_case(keywords_num_read[inum-1]);
			number_read[inum-1]=find_number_vector(number, numberlength);
			number_comp_read[inum-1]=numberlength;
		}else if(res==2){
			istr++;
			keywords_str_read[istr-1]=find_string(key, keylength);
			convert_string_in_lower_case(keywords_str_read[istr-1]);
			string_read[istr-1]=find_string_int(string, stringlength);
			string_length_read[istr-1]=stringlength;
		}
	}while (endoffile==0);
	fclose(f);
	free(key);
	free(string);
	free(number);
	
	//compare keywords number	
	n = (long)num_par_number;
	num_param_components = (long*)malloc(n*sizeof(long));
	num_param = (double**)malloc(n*sizeof(double*));
	for(i=0; i<num_par_number; i++){
		ok=0;
		for(j=0; j<inum; j++){
			if( strcmp (keywords_num_lower_case[i], keywords_num_read[j]) == 0){
				ok=1;
				num_param_components[i] = number_comp_read[j];
				num_param[i] = find_number_vector(number_read[j], number_comp_read[j]);
			}
		}
		if(ok==0){ //parameter not read
			num_param_components[i] = 1;
			num_param[i] = (double*)malloc(sizeof(double));
			num_param[i][0] = (double)number_novalue;
		}
	}
	
	//deallocate read arrays
	for(j=0; j<inum; j++){
		free(number_read[j]);
		free(keywords_num_read[j]);
	}
	free(number_read);
	free(keywords_num_read);
	free(number_comp_read);
	
	//assign parameter
	assign_numeric_parameters(par, land, times, sl, met, itools, num_param, num_param_components, keywords_num, flog);
	
	//deallocate keyword arrays
	for(i=0; i<num_par_number; i++){
		free(num_param[i]);
	}
	free(num_param);
	free(num_param_components);
	
	//compare keywords string
	n = (long)num_par_char;
	string_param = (char**)malloc(n*sizeof(char*));
	for(i=0; i<num_par_char; i++){
		ok=0;
		for(j=0; j<istr; j++){
			if( strcmp (keywords_char_lower_case[i], keywords_str_read[j] ) == 0){
				ok=1;
				string_param[i] = find_string(string_read[j], string_length_read[j]);
			}
		}
		if (ok==0) {
			n = strlen(string_novalue)+1;
			string_param[i] = (char*)malloc(n*sizeof(char));
			string_param[i] = strcpy(string_param[i], string_novalue);
		}
	}
	
	//deallocate read arrays
	for(j=0; j<istr; j++){
		free(string_read[j]);
		free(keywords_str_read[j]);
	}
	free(string_read);
	free(keywords_str_read);
	free(string_length_read);
	
	//assign parameter
	end += nmet;
	itools->met_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += nsoilprop;
	itools->soil_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 2;
	itools->horizon_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char); 
	beg = end;
	end += nfiles;
	files = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += otot;
	headerpoint = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += ootot;
	headerbasin = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 12;
	headersnow = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 12;	
	headerglac = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 6;	
	headersoil = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += ptTOT;		
	itools->point_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += nlstot;
	itools->lapserates_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 8;
	itools->meteostations_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	
	//add working path to the file name
	for (i=0; i<nfiles; i++) {
		if (strcmp(files[i], string_novalue) != 0) {
			temp = assign_string(files[i]);
			free(files[i]);
			files[i] = join_strings(WORKING_DIRECTORY, temp);
			free(temp);
		}
	}
	
	//deallocate keyword arrays	
	for(i=0; i<num_par_char; i++){
		free(string_param[i]);
	}
	free(string_param);
	
	n = (long)num_par_number;
	for (i=0; i<n; i++) {
		free(keywords_num_lower_case[i]);
	}
	free(keywords_num_lower_case);
	
	n = (long)num_par_char;
	for (i=0; i<n; i++) {
		free(keywords_char_lower_case[i]);
	}
	free(keywords_char_lower_case);
	
	
	//replace none value with some default values
	
	//Horizon
	for (j=0; j<2; j++) {
		if(strcmp(itools->horizon_col_names[j], string_novalue) == 0){
			free(itools->horizon_col_names[j]);
			if (j==0) { itools->horizon_col_names[j] = assign_string("AngleFromNorthClockwise");
			}else if (j==1) {  itools->horizon_col_names[j] = assign_string("HorizonHeight");
			}
		}
	}
	
	//HeaderPointFile
	for (j=0; j<otot; j++) {
		if(strcmp(headerpoint[j], string_novalue) == 0){
			free(headerpoint[j]);
			if(j==odate12){ headerpoint[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==oJDfrom0){ headerpoint[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==odaysfromstart){ headerpoint[j] = assign_string("TimeFromStart[days]"); 
			}else if(j==operiod){ headerpoint[j] = assign_string("Simulation_Period");
			}else if(j==orun){ headerpoint[j] = assign_string("Run");	
			}else if(j==opoint){ headerpoint[j] = assign_string("IDpoint"); 
			}else if(j==osnowover){ headerpoint[j] = assign_string("Psnow_over_canopy[mm]");     	
			}else if(j==orainover){ headerpoint[j] = assign_string("Prain_over_canopy[mm]"); 		
			}else if(j==oprecsnow){ headerpoint[j] = assign_string("Psnow_under_canopy[mm]"); 
			}else if(j==oprecrain){ headerpoint[j] = assign_string("Prain_under_canopy[mm]"); 		
			}else if(j==orainonsnow){ headerpoint[j] = assign_string("Prain_rain_on_snow[mm]");
			}else if(j==oV){ headerpoint[j] = assign_string("Wind_speed[m/s]");      	   
			}else if(j==oVdir){ headerpoint[j] = assign_string("Wind_direction[deg]");  
			}else if(j==oRH){ headerpoint[j] = assign_string("Relative_Humidity[-]");    
			}else if(j==oPa){ headerpoint[j] = assign_string("Pressure[mbar]");    
			}else if(j==oTa){ headerpoint[j] = assign_string("Tair[C]");    
			}else if(j==oTdew){ headerpoint[j] = assign_string("Tdew[C]");  
			}else if(j==oTg){ headerpoint[j] = assign_string("Tsurface[C]");    
			}else if(j==oTv){ headerpoint[j] = assign_string("Tvegetation[C]");    
			}else if(j==oTs){ headerpoint[j] = assign_string("Tcanopyair[C]");    
			}else if(j==oEB){ headerpoint[j] = assign_string("Surface_Energy_balance[W/m2]");    
			}else if(j==oG){ headerpoint[j] = assign_string("Soil_heat_flux[W/m2]");     
			}else if(j==oSWin){ headerpoint[j] = assign_string("SWin[W/m2]");  
			}else if(j==oSWb){ headerpoint[j] = assign_string("SWbeam[W/m2]");   
			}else if(j==oSWd){ headerpoint[j] = assign_string("SWdiff[W/m2]");  
			}else if(j==oLWin){ headerpoint[j] = assign_string("LWin[W/m2]"); 
			}else if(j==ominLWin){ headerpoint[j] = assign_string("LWin_min[W/m2]"); 
			}else if(j==omaxLWin){ headerpoint[j] = assign_string("LWin_max[W/m2]");
			}else if(j==oSW){ headerpoint[j] = assign_string("SWnet[W/m2]");     
			}else if(j==oLW){ headerpoint[j] = assign_string("LWnet[W/m2]");     
			}else if(j==oH){ headerpoint[j] = assign_string("H[W/m2]");      
			}else if(j==oLE){ headerpoint[j] = assign_string("LE[W/m2]");     
			}else if(j==ofc){ headerpoint[j] = assign_string("Canopy_fraction[-]");     
			}else if(j==oLSAI){ headerpoint[j] = assign_string("LSAI[m2/m2]");   
			}else if(j==oz0v){ headerpoint[j] = assign_string("z0veg[m]");    
			}else if(j==od0v){ headerpoint[j] = assign_string("d0veg[m]");    
			}else if(j==oEcan){ headerpoint[j] = assign_string("Estored_canopy[W/m2]");   
			}else if(j==oSWv){ headerpoint[j] = assign_string("SWv[W/m2]");    
			}else if(j==oLWv){ headerpoint[j] = assign_string("LWv[W/m2]");    
			}else if(j==oHv){ headerpoint[j] = assign_string("Hv[W/m2]");     
			}else if(j==oLEv){ headerpoint[j] = assign_string("LEv[W/m2]");    
			}else if(j==oHg0){ headerpoint[j] = assign_string("Hg_unveg[W/m2]");    
			}else if(j==oLEg0){ headerpoint[j] = assign_string("LEg_unveg[W/m2]");   
			}else if(j==oHg1){ headerpoint[j] = assign_string("Hg_veg[W/m2]");    
			}else if(j==oLEg1){ headerpoint[j] = assign_string("LEg_veg[W/m2]");   
			}else if(j==oevapsur){ headerpoint[j] = assign_string("Evap_surface[mm]");  
			}else if(j==otrasp){ headerpoint[j] = assign_string("Trasp_canopy[mm]");    
			}else if(j==owcan_rain){ headerpoint[j] = assign_string("Water_on_canopy[mm]");
			}else if(j==owcan_snow){ headerpoint[j] = assign_string("Snow_on_canopy[mm]");
			}else if(j==oQv){ headerpoint[j] = assign_string("Qvegetation[-]");   
			}else if(j==oQg){ headerpoint[j] = assign_string("Qsurface[-]");   
			}else if(j==oQa){ headerpoint[j] = assign_string("Qair[-]");   
			}else if(j==oQs){ headerpoint[j] = assign_string("Qcanopyair[-]");   
			}else if(j==oLobuk){ headerpoint[j] = assign_string("LObukhov[m]");
			}else if(j==oLobukcan){ headerpoint[j] = assign_string("LObukhovcanopy[m]");
			}else if(j==outop){ headerpoint[j] = assign_string("Wind_speed_top_canopy[m/s]");    
			}else if(j==odecay){ headerpoint[j] = assign_string("Decay_of_K_in_canopy[-]");   
			}else if(j==oSWup){ headerpoint[j] = assign_string("SWup[W/m2]");   
			}else if(j==oLWup){ headerpoint[j] = assign_string("LWup[W/m2]");   
			}else if(j==oHup){ headerpoint[j] = assign_string("Hup[W/m2]");    
			}else if(j==oLEup){ headerpoint[j] = assign_string("LEup[W/m2]");   
			}else if(j==osnowdepth){ headerpoint[j] = assign_string("snow_depth[mm]"); 
			}else if(j==oSWE){ headerpoint[j] = assign_string("snow_water_equivalent[mm]"); 
			}else if(j==osnowdens){ headerpoint[j] = assign_string("snow_density[kg/m3]"); 
			}else if(j==osnowT){ headerpoint[j] = assign_string("snow_temperature[C]"); 	
			}else if(j==omrsnow){ headerpoint[j] = assign_string("snow_melted[mm]"); 
			}else if(j==osrsnow){ headerpoint[j] = assign_string("snow_subl[mm]"); 
			}else if(j==oblowingsnowtrans){ headerpoint[j] = assign_string("snow_blown_away[mm]"); 
			}else if(j==oblowingsnowsubl){ headerpoint[j] = assign_string("snow_subl_while_blown[mm]");			
			}else if(j==oglacdepth){ headerpoint[j] = assign_string("glac_depth[mm]"); 
			}else if(j==oGWE){ headerpoint[j] = assign_string("glac_water_equivalent[mm]"); 
			}else if(j==oglacdens){ headerpoint[j] = assign_string("glac_density[kg/m3]"); 
			}else if(j==oglacT){ headerpoint[j] = assign_string("glac_temperature[C]"); 
			}else if(j==omrglac){ headerpoint[j] = assign_string("glac_melted[mm]"); 
			}else if(j==osrglac){ headerpoint[j] = assign_string("glac_subl[mm]"); 
			}else if(j==othawed){ headerpoint[j] = assign_string("thawed_soil_depth[mm]"); 
			}else if(j==owtable){ headerpoint[j] = assign_string("water_table_depth[mm]"); 
			}
		}
	}
	
	//HeaderBasinFile
	for (j=0; j<ootot; j++) {
		if(strcmp(headerbasin[j], string_novalue) == 0){
			free(headerbasin[j]);
			if(j==oodate12){ headerbasin[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==ooJDfrom0){ headerbasin[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==oodaysfromstart){ headerbasin[j] = assign_string("TimeFromStart[days]");    
			}else if(j==ooperiod){ headerbasin[j] = assign_string("Simulation_Period");
			}else if(j==oorun){ headerbasin[j] = assign_string("Run");	
			}else if(j==ooprecrain){ headerbasin[j] = assign_string("Prain_below_canopy[mm]");     	
			}else if(j==ooprecsnow){ headerbasin[j] = assign_string("Psnow_below_canopy[mm]");     	
			}else if(j==oorainover){ headerbasin[j] = assign_string("Prain_above_canopy[mm]");     	
			}else if(j==oosnowover){ headerbasin[j] = assign_string("Prain_above_canopy[mm]");     	
			}else if(j==ooTa){ headerbasin[j] = assign_string("Tair[C]");     	
			}else if(j==ooTg){ headerbasin[j] = assign_string("Tsurface[C]");     	
			}else if(j==ooTv){ headerbasin[j] = assign_string("Tvegetation[C]");     	
			}else if(j==ooevapsur){ headerbasin[j] = assign_string("Evap_surface[mm]");     	
			}else if(j==ootrasp){ headerbasin[j] = assign_string("Transpiration_canopy[mm]");     	
			}else if(j==ooLE){ headerbasin[j] = assign_string("LE[W/m2]");     	
			}else if(j==ooH){ headerbasin[j] = assign_string("H[W/m2]");     	
			}else if(j==ooSW){ headerbasin[j] = assign_string("SW[W/m2]");     	
			}else if(j==ooLW){ headerbasin[j] = assign_string("LW[W/m2]");     	
			}else if(j==ooLEv){ headerbasin[j] = assign_string("LEv[W/m2]");     	
			}else if(j==ooHv){ headerbasin[j] = assign_string("Hv[W/m2]");     	
			}else if(j==ooSWv){ headerbasin[j] = assign_string("SWv[W/m2]");     	
			}else if(j==ooLWv){ headerbasin[j] = assign_string("LWv[W/m2]");     	
			}else if(j==ooSWin){ headerbasin[j] = assign_string("SWin[W/m2]");     	
			}else if(j==ooLWin){ headerbasin[j] = assign_string("LWin[W/m2]");     	
			}else if(j==oomasserror){ headerbasin[j] = assign_string("Mass_balance_error[mm]");     	
			}
		}
	}
	
	//HeaderSnowFile
	for (j=0; j<12; j++) {
		if(strcmp(headersnow[j], string_novalue) == 0){
			free(headersnow[j]);
			if(j==0){ headersnow[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==1){ headersnow[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==2){ headersnow[j] = assign_string("TimeFromStart[days]");  
			}else if(j==3){ headersnow[j] = assign_string("Simulation_Period");
			}else if(j==4){ headersnow[j] = assign_string("Run");					
			}else if(j==5){ headersnow[j] = assign_string("IDpoint"); 
			}else if(j==6){ headersnow[j] = assign_string("WaterEquivalent[mm]");     	
			}else if(j==7){ headersnow[j] = assign_string("Depth[mm]");     	
			}else if(j==8){ headersnow[j] = assign_string("Density[kg/m3]");     	
			}else if(j==9){ headersnow[j] = assign_string("Temperature[C]");     	
			}else if(j==10){ headersnow[j] = assign_string("IceContent[-]");     	
			}else if(j==11){ headersnow[j] = assign_string("LiqWaterContent[-]");     	
			}
		}
	}
	
	//HeaderGlacierFile
	for (j=0; j<12; j++) {
		if(strcmp(headerglac[j], string_novalue) == 0){
			
			free(headerglac[j]);
			if(j==0){ headerglac[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==1){ headerglac[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==2){ headerglac[j] = assign_string("TimeFromStart[days]");   
			}else if(j==3){ headerglac[j] = assign_string("Simulation_Period");
			}else if(j==4){ headerglac[j] = assign_string("Run");					
			}else if(j==5){ headerglac[j] = assign_string("IDpoint"); 
			}else if(j==6){ headerglac[j] = assign_string("WaterEquivalent[mm]");     	
			}else if(j==7){ headerglac[j] = assign_string("Depth[mm]");     	
			}else if(j==8){ headerglac[j] = assign_string("Density[kg/m3]");     	
			}else if(j==9){ headerglac[j] = assign_string("Temperature[C]");     	
			}else if(j==10){ headerglac[j] = assign_string("IceContent[-]");     	
			}else if(j==11){ headerglac[j] = assign_string("LiqWaterContent[-]");     	
			}
		}
	}				
	
	
	//HeaderGlacierFile
	for (j=0; j<6; j++) {
		if(strcmp(headersoil[j], string_novalue) == 0){
			free(headersoil[j]);
			if(j==0){ headersoil[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==1){ headersoil[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==2){ headersoil[j] = assign_string("TimeFromStart[days]"); 
			}else if(j==3){ headersoil[j] = assign_string("Simulation_Period");
			}else if(j==4){ headersoil[j] = assign_string("Run");									
			}else if(j==5){ headersoil[j] = assign_string("IDpoint"); 
			}
		}	
	}		
	
	return 1;
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

void assign_numeric_parameters(PAR *par, LAND *land, TIMES *times, SOIL *sl, METEO *met, INIT_TOOLS *itools, double **num_param, long *num_param_components, char **keyword, FILE *flog){
	
	short occurring;
	long cod, codn, i, j, k, n, nsoillayers, nmeteo_stations, npoints;
	double a, minDt=1.E99;
	
	fprintf(flog,"\n");
	
	par->print=0;
	
	//find components of times->Dt_vector
	cod = 0;
	n = (long)max_cols_time_steps_file + 1;
	times->Dt_vector=(double *)malloc(n*sizeof(double));
	times->Dt_vector[0] = 0.;//it is the space for the date in case of time variable time step
	times->Dt_vector[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 1);
	for (i=2; i<n; i++) {
		if (i <= num_param_components[cod]) {
			times->Dt_vector[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
		}else {
			times->Dt_vector[i] = (double)number_novalue;
		}
	}
	for (i=1; i<n; i++) {
		if((long)times->Dt_vector[i] != number_novalue){
			if (times->Dt_vector[i] < minDt) minDt = times->Dt_vector[i];
		}
	}
	
	//init date
	cod = 1;
	par->init_date = new_doublevector(num_param_components[cod]);
	par->init_date->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 1);
	par->init_date->co[1] = convert_dateeur12_JDfrom0(par->init_date->co[1]);
	for (i=2; i<=par->init_date->nh; i++) {
		par->init_date->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, (double)number_novalue, 1);
		par->init_date->co[i] = convert_dateeur12_JDfrom0(par->init_date->co[i]);
	}
	
	//end date
	cod = 2;
	par->end_date = new_doublevector(num_param_components[cod]);
	if (par->end_date->nh != par->init_date->nh) t_error("Error:: End date has a number of components different from Init Date");
	par->end_date->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 1);
	par->end_date->co[1] = convert_dateeur12_JDfrom0(par->end_date->co[1]);
	for (i=2; i<=par->end_date->nh; i++) {
		par->end_date->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, (double)number_novalue, 1);
		par->end_date->co[i] = convert_dateeur12_JDfrom0(par->end_date->co[i]);
	}
	
	//run times
	cod = 3;
	par->run_times = new_longvector(par->init_date->nh);
	par->run_times->co[1] = (long)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 1., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->run_times->co[i] = (long)assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, (double)par->run_times->co[i-1], 0);
	}
			
	par->ST = assignation_number(flog, 4, 0, keyword, num_param, num_param_components, 0., 0);
	
	cod = 5;
	par->Dtplot_discharge = new_doublevector(par->init_date->nh);
	par->Dtplot_discharge->co[1] = 3600.*assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->Dtplot_discharge->co[i] = 3600.*assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_discharge->co[i-1], 0);
	}
	par->plot_discharge_with_Dt_integration = new_shortvector(par->init_date->nh);
	par->state_discharge = 0;
	for (i=1; i<=par->init_date->nh; i++) {
		if(par->Dtplot_discharge->co[i] > 1.E-5 && par->Dtplot_discharge->co[i] <= minDt){
			par->plot_discharge_with_Dt_integration->co[i]=1;
		}else{
			par->plot_discharge_with_Dt_integration->co[i]=0;
		}
		if(par->Dtplot_discharge->co[i] > 1.E-5) par->state_discharge = 1;
	}
	
	cod = 6;
	par->Dtplot_point = new_doublevector(par->init_date->nh);
	par->Dtplot_point->co[1] = 3600.*assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->Dtplot_point->co[i] = 3600.*assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_point->co[i-1], 0);
	}
	par->plot_point_with_Dt_integration = new_shortvector(par->init_date->nh);
	par->state_pixel = 0;
	for (i=1; i<=par->init_date->nh; i++) {
		if(par->Dtplot_point->co[i] > 1.E-5 && par->Dtplot_point->co[i] <= minDt){
			par->plot_point_with_Dt_integration->co[i]=1;
		}else{
			par->plot_point_with_Dt_integration->co[i]=0;
		}
		if(par->Dtplot_point->co[i] > 1.E-5) par->state_pixel = 1;
	}
	
	cod = 7;
	par->Dtplot_basin = new_doublevector(par->init_date->nh);
	par->Dtplot_basin->co[1] = 3600.*assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->Dtplot_basin->co[i] = 3600.*assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_basin->co[i-1], 0);
	}
	par->plot_basin_with_Dt_integration = new_shortvector(par->init_date->nh);
	par->state_basin = 0;
	for (i=1; i<=par->init_date->nh; i++) {
		if(par->Dtplot_basin->co[i] > 1.E-5 && par->Dtplot_basin->co[i] <= minDt){
			par->plot_basin_with_Dt_integration->co[i]=1;
		}else{
			par->plot_basin_with_Dt_integration->co[i]=0;
		}
		if(par->Dtplot_basin->co[i] > 1.E-5) par->state_basin = 1;
	}
	
	
	par->lowpass = (long)assignation_number(flog, 8, 0, keyword, num_param, num_param_components, 0., 0);
	par->lowpass_curvatures = (long)assignation_number(flog, 9, 0, keyword, num_param, num_param_components, 0., 0);
	par->sky = (short)assignation_number(flog, 10, 0, keyword, num_param, num_param_components, 0., 0);
	par->format_out = (short)assignation_number(flog, 11, 0, keyword, num_param, num_param_components, 3., 0);
	par->point_sim = (short)assignation_number(flog, 12, 0, keyword, num_param, num_param_components, 0., 0);
	par->recover = (short)assignation_number(flog, 13, 0, keyword, num_param, num_param_components, 0., 0);
	
	
	//land cover types
	par->n_landuses = (long)assignation_number(flog, 14, 0, keyword, num_param, num_param_components, 1., 0);
	
	land->ty = new_doublematrix(par->n_landuses, nlandprop);
	
	land->ty->co[1][jz0] = assignation_number(flog, 15, 0, keyword, num_param, num_param_components, 10., 0);
	land->ty->co[1][jz0thressoil] = assignation_number(flog, 16, 0, keyword, num_param, num_param_components, land->ty->co[1][jz0], 0);
	land->ty->co[1][jHveg] = assignation_number(flog, 17, 0, keyword, num_param, num_param_components, 1000., 0);
	land->ty->co[1][jz0thresveg] = assignation_number(flog, 18, 0, keyword, num_param, num_param_components, land->ty->co[1][jHveg], 0);
	land->ty->co[1][jz0thresveg2] = assignation_number(flog, 19, 0, keyword, num_param, num_param_components, land->ty->co[1][jz0thresveg], 0);
	land->ty->co[1][jLSAI] = assignation_number(flog, 20, 0, keyword, num_param, num_param_components, 1., 0);
	land->ty->co[1][jcf] = assignation_number(flog, 21, 0, keyword, num_param, num_param_components, 0., 0);
	land->ty->co[1][jdecay0] = assignation_number(flog, 22, 0, keyword, num_param, num_param_components, 2.5, 0);
	land->ty->co[1][jexpveg] = assignation_number(flog, 23, 0, keyword, num_param, num_param_components, 1., 0);
	land->ty->co[1][jroot] = assignation_number(flog, 24, 0, keyword, num_param, num_param_components, 300., 0);
	land->ty->co[1][jrs] = assignation_number(flog, 25, 0, keyword, num_param, num_param_components, 60., 0);
	land->ty->co[1][jvR_vis] = assignation_number(flog, 26, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty->co[1][jvR_nir] = assignation_number(flog, 27, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty->co[1][jvT_vis] = assignation_number(flog, 28, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty->co[1][jvT_nir] = assignation_number(flog, 29, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty->co[1][jvCh] = assignation_number(flog, 30, 0, keyword, num_param, num_param_components, 0., 0);
	land->ty->co[1][jcd] = assignation_number(flog, 31, 0, keyword, num_param, num_param_components, 2., 0);
	land->ty->co[1][ja_vis_dry] = assignation_number(flog, 32, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty->co[1][ja_nir_dry] = assignation_number(flog, 33, 0, keyword, num_param, num_param_components, land->ty->co[1][ja_vis_dry], 0);
	land->ty->co[1][ja_vis_sat] = assignation_number(flog, 34, 0, keyword, num_param, num_param_components, land->ty->co[1][ja_vis_dry], 0);
	land->ty->co[1][ja_nir_sat] = assignation_number(flog, 35, 0, keyword, num_param, num_param_components, land->ty->co[1][ja_nir_dry], 0);
	land->ty->co[1][jemg] = assignation_number(flog, 36, 0, keyword, num_param, num_param_components, 0.96, 0);
	land->ty->co[1][jcm] = assignation_number(flog, 37, 0, keyword, num_param, num_param_components, 0.5, 0);
	land->ty->co[1][jN] = assignation_number(flog, 38, 0, keyword, num_param, num_param_components, 0., 0);
	land->ty->co[1][jdv] = assignation_number(flog, 39, 0, keyword, num_param, num_param_components, 50., 0);
	
	for (i=2; i<=par->n_landuses; i++) {
		for (j=1; j<=nlandprop; j++) {
			land->ty->co[i][j] = assignation_number(flog, 15+j-1, i-1, keyword, num_param, num_param_components, land->ty->co[i-1][j], 0);
		}
	}
	
	//former block 2
	par->imp = assignation_number(flog, 40, 0, keyword, num_param, num_param_components, 2., 0);
	par->TolVWb = assignation_number(flog, 41, 0, keyword, num_param, num_param_components, 1.E-8, 0);
	par->RelTolVWb = RelativeErrorRichards;
	par->MaxErrWb = 1.E99;
	par->MaxiterTol = (long)assignation_number(flog, 42, 0, keyword, num_param, num_param_components, 100., 0);
	par->TolCG = assignation_number(flog, 43, 0, keyword, num_param, num_param_components, 0.01, 0);
	par->UpdateK = 0;
	par->min_lambda_wat = assignation_number(flog, 44, 0, keyword, num_param, num_param_components, 1.E-10, 0);
	par->max_times_min_lambda_wat = (long)assignation_number(flog, 45, 0, keyword, num_param, num_param_components, 10.0, 0);
	par->exit_lambda_min_wat = (short)assignation_number(flog, 46, 0, keyword, num_param, num_param_components, 1., 0);
	par->max_times_halving_time_step_wat = (long)assignation_number(flog, 47, 0, keyword, num_param, num_param_components, 10.0, 0);	
	par->harm_or_arit_mean_normal = 1;	
	par->harm_or_arit_mean_parallel = 0;	
	par->gamma_m = assignation_number(flog, 48, 0, keyword, num_param, num_param_components, 2./3., 0);
	par->thres_hsup_1 = assignation_number(flog, 49, 0, keyword, num_param, num_param_components, 0., 0);
	par->thres_hsup_2 = assignation_number(flog, 50, 0, keyword, num_param, num_param_components, 50., 0);
	par->Ks_channel = assignation_number(flog, 51, 0, keyword, num_param, num_param_components, 20., 0);
	par->thres_hchannel = assignation_number(flog, 52, 0, keyword, num_param, num_param_components, 50., 0);
	par->w_dx = assignation_number(flog, 53, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->depr_channel = assignation_number(flog, 54, 0, keyword, num_param, num_param_components, 500., 0);
	par->max_courant_land = assignation_number(flog, 55, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->max_courant_channel = assignation_number(flog, 56, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->min_hsup_land = assignation_number(flog, 57, 0, keyword, num_param, num_param_components, 1.0, 0);
	par->min_hsup_channel = assignation_number(flog, 58, 0, keyword, num_param, num_param_components, 1.0, 0);
	par->dtmin_sup = assignation_number(flog, 59, 0, keyword, num_param, num_param_components, 0.1, 0);
	//former block 3
	par->latitude = assignation_number(flog, 60, 0, keyword, num_param, num_param_components, 45., 0);
	par->longitude = assignation_number(flog, 61, 0, keyword, num_param, num_param_components, 0., 0);
	par->Vmin = assignation_number(flog, 62, 0, keyword, num_param, num_param_components, 0.5, 0);
	par->RHmin = assignation_number(flog, 63, 0, keyword, num_param, num_param_components, 10., 0);
	par->alpha_snow = assignation_number(flog, 64, 0, keyword, num_param, num_param_components, 1.E2, 0);
	par->nsurface = (long)assignation_number(flog, 65, 0, keyword, num_param, num_param_components, 0., 0);
	par->tol_energy = assignation_number(flog, 66, 0, keyword, num_param, num_param_components, 1.E-4, 0);
	par->maxiter_energy = (long)assignation_number(flog, 67, 0, keyword, num_param, num_param_components, 100., 0);
	par->min_lambda_en = assignation_number(flog, 68, 0, keyword, num_param, num_param_components, 1.E-5, 0);
	par->max_times_min_lambda_en = (long)assignation_number(flog, 69, 0, keyword, num_param, num_param_components, 0.0, 0);
	par->exit_lambda_min_en = (short)assignation_number(flog, 70, 0, keyword, num_param, num_param_components, 0., 0);
	par->max_times_halving_time_step_en = (long)assignation_number(flog, 71, 0, keyword, num_param, num_param_components, 5.0, 0);	
	par->maxiter_canopy = (long)assignation_number(flog, 72, 0, keyword, num_param, num_param_components, 3., 0);
	par->maxiter_Businger = (long)assignation_number(flog, 73, 0, keyword, num_param, num_param_components, 5., 0);
	par->maxiter_Ts = (long)assignation_number(flog, 74, 0, keyword, num_param, num_param_components, 2., 0);
	par->maxiter_Loc = (long)assignation_number(flog, 75, 0, keyword, num_param, num_param_components, 3., 0);
	par->stabcorr_incanopy = (short)assignation_number(flog, 76, 0, keyword, num_param, num_param_components, 1., 0);
	par->iobsint = (short)assignation_number(flog, 77, 0, keyword, num_param, num_param_components, 1., 0);
	par->dn = assignation_number(flog, 78, 0, keyword, num_param, num_param_components, 1., 0);
	par->slopewt = assignation_number(flog, 79, 0, keyword, num_param, num_param_components, 0., 0);
	par->curvewt = assignation_number(flog, 80, 0, keyword, num_param, num_param_components, 0., 0);
	par->Zboundary = assignation_number(flog, 81, 0, keyword, num_param, num_param_components, 1.E20, 0);
	par->Tboundary = assignation_number(flog, 82, 0, keyword, num_param, num_param_components, 20., 0);
	par->Fboundary = assignation_number(flog, 83, 0, keyword, num_param, num_param_components, 0., 0);
	//former block 4
	itools->snow0 = assignation_number(flog, 84, 0, keyword, num_param, num_param_components, 0., 0);
	itools->rhosnow0 = assignation_number(flog, 85, 0, keyword, num_param, num_param_components, 200., 0);
	itools->Tsnow0 = assignation_number(flog, 86, 0, keyword, num_param, num_param_components, -3., 0);
	itools->agesnow0 = assignation_number(flog, 87, 0, keyword, num_param, num_param_components, 0., 0);
	par->T_rain = assignation_number(flog, 88, 0, keyword, num_param, num_param_components, 3., 0);
	par->T_snow = assignation_number(flog, 89, 0, keyword, num_param, num_param_components, -1., 0);
	par->dew = (short)assignation_number(flog, 90, 0, keyword, num_param, num_param_components, 0., 0);
	par->aep = assignation_number(flog, 91, 0, keyword, num_param, num_param_components, 10., 0);
	par->avo = assignation_number(flog, 92, 0, keyword, num_param, num_param_components, 0.9, 0);
	par->airo = assignation_number(flog, 93, 0, keyword, num_param, num_param_components, 0.65, 0); 
	par->Sr = assignation_number(flog, 94, 0, keyword, num_param, num_param_components, 0.02, 0);  
	par->epsilon_snow = assignation_number(flog, 95, 0, keyword, num_param, num_param_components, 0.98, 0); 
	par->z0_snow = 0.001*assignation_number(flog, 96, 0, keyword, num_param, num_param_components, 0.1, 0); 
	par->snowcorrfact = assignation_number(flog, 97, 0, keyword, num_param, num_param_components, 1., 0); 
	par->raincorrfact = assignation_number(flog, 98, 0, keyword, num_param, num_param_components, 1., 0); 
	par->snowlayer_max = (long)assignation_number(flog, 99, 0, keyword, num_param, num_param_components, 3., 0);
	par->snowlayer_inf = (long)assignation_number(flog, 100, 0, keyword, num_param, num_param_components, 2., 0);
	par->snow_maxpor = assignation_number(flog, 101, 0, keyword, num_param, num_param_components, 0.7, 0);
	par->drysnowdef_rate = assignation_number(flog, 102, 0, keyword, num_param, num_param_components, 1., 0);
	par->snow_density_cutoff = assignation_number(flog, 103, 0, keyword, num_param, num_param_components, 100., 0);
	par->wetsnowdef_rate = assignation_number(flog, 104, 0, keyword, num_param, num_param_components, 1.5, 0);
	par->snow_viscosity = assignation_number(flog, 105, 0, keyword, num_param, num_param_components, 1.E6, 0);
	par->fetch_up = assignation_number(flog, 106, 0, keyword, num_param, num_param_components, 1000., 0);
	par->fetch_down = assignation_number(flog, 107, 0, keyword, num_param, num_param_components, 100., 0);
	par->Wice_PBSM = assignation_number(flog, 108, 0, keyword, num_param, num_param_components, 0., 0);
	par->Dt_PBSM = assignation_number(flog, 109, 0, keyword, num_param, num_param_components, times->Dt_vector[1], 0);
	par->snow_smin = assignation_number(flog, 110, 0, keyword, num_param, num_param_components, 30., 0);
	par->snow_smax = assignation_number(flog, 111, 0, keyword, num_param, num_param_components, 80., 0);
	par->snow_curv = assignation_number(flog, 112, 0, keyword, num_param, num_param_components, -200., 0);
	
	//former blocks 5/6
	par->Dmin = new_doublevector(par->snowlayer_max);
	par->Dmax = new_doublevector(par->snowlayer_max);
	for (i=1; i<=par->snowlayer_max; i++) {
		par->Dmin->co[i] = assignation_number(flog, 113, i-1, keyword, num_param, num_param_components, 10., 0);
		par->Dmax->co[i] = assignation_number(flog, 114, i-1, keyword, num_param, num_param_components, 100., 0);	
	}
	par->Dmax->co[par->snowlayer_inf] = 1.E10;
	for (i=1; i<=par->snowlayer_max-1; i++) {
		if(par->Dmin->co[i] < par->Dmin->co[par->snowlayer_max]){
			printf("Error:: The last snow layer must be the shallowest, you set Dmin[%ld]=%f < Dmin[%ld]=%f. This is not possible\n",i,par->Dmin->co[i],par->snowlayer_max,par->Dmin->co[par->snowlayer_max]);
			stop_execution();
			t_error("Not Possible to Continue");
		}
	}
	
	//former block 7
	itools->Dglac0 = assignation_number(flog, 115, 0, keyword, num_param, num_param_components, 0., 0);
	itools->rhoglac0 = assignation_number(flog, 116, 0, keyword, num_param, num_param_components, 800., 0);
	itools->Tglac0 = assignation_number(flog, 117, 0, keyword, num_param, num_param_components, -3., 0);
	par->Sr_glac = assignation_number(flog, 118, 0, keyword, num_param, num_param_components, 0.02, 0);
	par->glaclayer_max = (long)assignation_number(flog, 119, 0, keyword, num_param, num_param_components, 0, 0);
	par->glaclayer_inf = 1;
	
	//former block 8
	if(par->glaclayer_max > 0){
		par->Dmin_glac = new_doublevector(par->glaclayer_max);
		par->Dmax_glac = new_doublevector(par->glaclayer_max);		
		for (i=1; i<=par->glaclayer_max; i++) {
			par->Dmin_glac->co[i] = assignation_number(flog, 120, i-1, keyword, num_param, num_param_components, 10., 0);
			par->Dmax_glac->co[i] = assignation_number(flog, 121, i-1, keyword, num_param, num_param_components, 100., 0);		
		}
		par->Dmax_glac->co[par->glaclayer_max] = 1.E10;
	}
	
	//former block 9
	par->state_turb = 1;
	par->state_lwrad = (short)assignation_number(flog, 122, 0, keyword, num_param, num_param_components, 9., 0);
	par->monin_obukhov = (short)assignation_number(flog, 123, 0, keyword, num_param, num_param_components, 1., 0);
		
	//distributed option file
	par->wat_balance = (short)assignation_number(flog, 124, 0, keyword, num_param, num_param_components, 0., 0);
	par->en_balance = (short)assignation_number(flog, 125, 0, keyword, num_param, num_param_components, 0., 0);
	par->blowing_snow = (short)assignation_number(flog, 126, 0, keyword, num_param, num_param_components, 0., 0);
	par->state_px_coord = (short)assignation_number(flog, 127, 0, keyword, num_param, num_param_components, 1., 0);
	
	cod = 128;
	npoints = 0;
	for (j=1; j<=16; j++) {
		if (npoints < num_param_components[cod + j-1]) npoints = num_param_components[cod + j-1];
	}
	
	if (par->point_sim == 1) {
		par->chkpt = new_doublematrix(npoints, ptTOT);
	}else {
		par->chkpt = new_doublematrix(npoints, 3);
	}
	
	for (i=1; i<=par->chkpt->nrh; i++) {
		for (j=1; j<=par->chkpt->nch; j++) {
			par->chkpt->co[i][j] = assignation_number(flog, cod + j-1, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
		}
	}
		
	cod = 147;
	par->saving_points = new_doublevector(num_param_components[cod]);
	for (i=1; i<=par->saving_points->nh; i++) {
		par->saving_points->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 0);
	}
	
	par->distr_stat = 0;
	par->output_soil = assignation_number(flog, 148, 0, keyword, num_param, num_param_components, 0., 0);
	par->output_snow = assignation_number(flog, 149, 0, keyword, num_param, num_param_components, 0., 0);
	par->output_glac = assignation_number(flog, 150, 0, keyword, num_param, num_param_components, 0., 0);
	par->output_surfenergy = assignation_number(flog, 151, 0, keyword, num_param, num_param_components, 0., 0);
	par->output_vegetation = assignation_number(flog, 152, 0, keyword, num_param, num_param_components, 0., 0);
	par->output_meteo = assignation_number(flog, 153, 0, keyword, num_param, num_param_components, 0., 0);
	
	cod = 154;
	codn = 155;
	if(num_param_components[cod] != num_param_components[codn]){
		printf("ERROR: Number of components of parameters %s and %s must be equal\n",keyword[cod],keyword[codn]);
		stop_execution();
		t_error("Not Possible To Continue (6) ");
	}
	times->JD_plots = new_doublevector(num_param_components[cod] + num_param_components[codn]);
	for (i=1; i<=(long)(times->JD_plots->nh/2.); i++) {
		times->JD_plots->co[2*i-1] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 0);
		times->JD_plots->co[2*i  ] = assignation_number(flog, codn, i-1, keyword, num_param, num_param_components, 0., 0);				
	}
	if (times->JD_plots->nh == 2 && times->JD_plots->co[1] < 1.E-5 && times->JD_plots->co[2] < 1.E-5) {
		free_doublevector(times->JD_plots);
		times->JD_plots = new_doublevector(1);
		initialize_doublevector(times->JD_plots, 0.);
	}	
	if (times->JD_plots->nh > 1) {
		for (i=1; i<=times->JD_plots->nh; i++) {
			times->JD_plots->co[i] = convert_dateeur12_JDfrom0(times->JD_plots->co[i]);
		}
	}
	
	//initial condition on the water pressure
	par->nsoiltypes = (long)assignation_number(flog, 156, 0, keyword, num_param, num_param_components, 1., 0);
	if (par->nsoiltypes < 1) par->nsoiltypes = 1;
	
	cod = 157;
	sl->init_water_table_height = new_doublevector(par->nsoiltypes);
	sl->init_water_table_height->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (k=2; k<=par->nsoiltypes; k++) {
		sl->init_water_table_height->co[k] = assignation_number(flog, cod, k-1, keyword, num_param, num_param_components, sl->init_water_table_height->co[k-1], 0);
	}
	
	//soil properties and discretization
	par->soil_type_land_default = (long)assignation_number(flog, 158, 0, keyword, num_param, num_param_components, 1., 0);
	if(par->soil_type_land_default<1 || par->soil_type_land_default>par->nsoiltypes) t_error("soil_type_land_default lower than 0 or higher than soil types numbers");
	par->soil_type_chan_default = (long)assignation_number(flog, 159, 0, keyword, num_param, num_param_components, 1., 0);
	if(par->soil_type_chan_default<1 || par->soil_type_chan_default>par->nsoiltypes) t_error("soil_type_chan_default lower than 0 or higher than soil types numbers");
	
	cod = 160;
	a = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);
	
	//there is a specific layer discretization
	if ( (long)a != number_novalue && num_param_components[cod] > 1) {
		
		nsoillayers = num_param_components[cod];
		sl->pa = new_doubletensor(par->nsoiltypes, nsoilprop, nsoillayers);
		
		sl->pa->co[1][jdz][1] = a;
		for (i=2; i<=sl->pa->nch; i++) {
			sl->pa->co[1][jdz][i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, sl->pa->co[1][jdz][i-1], 0);
		}
				
	}else {
		
		if((long)a==number_novalue) a=100.;
		nsoillayers = (long)assignation_number(flog, 161, 0, keyword, num_param, num_param_components, 5., 0);
		
		sl->pa = new_doubletensor(par->nsoiltypes, nsoilprop, nsoillayers);
		
		for (i=1; i<=sl->pa->nch; i++) {
			sl->pa->co[1][jdz][i] = a;
		}

	}
	
	//first layer
	i = 1;
	sl->pa->co[1][jpsi][i] = assignation_number(flog, 162, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	sl->pa->co[1][jT][i] = assignation_number(flog, 163, i-1, keyword, num_param, num_param_components, 5., 0);
	sl->pa->co[1][jKn][i] = assignation_number(flog, 164, i-1, keyword, num_param, num_param_components, 1.E-4, 0); 
	sl->pa->co[1][jKl][i] = assignation_number(flog, 165, i-1, keyword, num_param, num_param_components, 1.E-4, 0); 
	sl->pa->co[1][jres][i] = assignation_number(flog, 166, i-1, keyword, num_param, num_param_components, 0.05, 0);  
	sl->pa->co[1][jwp][i] = assignation_number(flog, 167, i-1, keyword, num_param, num_param_components, 0.15, 0);  
	sl->pa->co[1][jfc][i] = assignation_number(flog, 168, i-1, keyword, num_param, num_param_components, 0.25, 0);  
	sl->pa->co[1][jsat][i] = assignation_number(flog, 169, i-1, keyword, num_param, num_param_components, 0.5, 0);  
	sl->pa->co[1][ja][i] = assignation_number(flog, 170, i-1, keyword, num_param, num_param_components, 0.004, 0);  
	sl->pa->co[1][jns][i] = assignation_number(flog, 171, i-1, keyword, num_param, num_param_components, 1.3, 0);  
	sl->pa->co[1][jv][i] = assignation_number(flog, 172, i-1, keyword, num_param, num_param_components, 0.5, 0);  
	sl->pa->co[1][jkt][i] = assignation_number(flog, 173, i-1, keyword, num_param, num_param_components, 2.5, 0);  
	sl->pa->co[1][jct][i] = assignation_number(flog, 174, i-1, keyword, num_param, num_param_components, 1.E6, 0);  
	sl->pa->co[1][jss][i] = assignation_number(flog, 175, i-1, keyword, num_param, num_param_components, 1.E-7, 0);  
	
	//other layers
	for (i=2; i<=sl->pa->nch; i++) {
		for (j=2; j<=sl->pa->nrh; j++) {
			sl->pa->co[1][j][i] = assignation_number(flog, 162 + j-2, i-1, keyword, num_param, num_param_components, sl->pa->co[1][j][i-1], 0);
		}
	}
	
	//other soil types
	for (k=2; k<=par->nsoiltypes; k++) {
		for (i=1; i<=sl->pa->nch; i++) {
			for (j=1; j<=sl->pa->nrh; j++) {
				sl->pa->co[k][j][i] = sl->pa->co[1][j][i];
			}
		}	
	}
	
	//use water table for water pressure
	for (k=1; k<=par->nsoiltypes; k++) {
		occurring = 0;//check if psi initial has at least one novalue
		for (i=1; i<=sl->pa->nch; i++) {
			if ( (long)sl->pa->co[k][jpsi][i] == number_novalue) occurring = 1;
		}
		if (occurring == 0) sl->init_water_table_height->co[k] = (double)number_novalue;
	}
	
	cod = 176;
	sl->pa_bed = new_doubletensor(1, nsoilprop, nsoillayers);
	for (i=1; i<=nsoillayers; i++) {
		sl->pa_bed->co[1][jdz][i] = sl->pa->co[1][jdz][i];
	}
	for (j=1; j<=nsoilprop; j++) {
		if(j != jdz){
			sl->pa_bed->co[1][j][1] = assignation_number(flog, cod+j-2, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);
		}
	}
	for (i=2; i<=nsoillayers; i++) {
		for (j=1; j<=nsoilprop; j++) {
			if(j != jdz) sl->pa_bed->co[1][j][i] = assignation_number(flog, cod+j-2, i-1, keyword, num_param, num_param_components, sl->pa_bed->co[1][j][i-1], 0);
		}
	}	
	
	//meteo stations
	cod = 190;
	met->imeteo_stations = new_longvector(num_param_components[cod]);
	met->imeteo_stations->co[1] = (long)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);
	if ( met->imeteo_stations->co[1] != number_novalue ) {
		for (i=2; i<=num_param_components[cod]; i++) {
			met->imeteo_stations->co[i] = (long)assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 1);
		}
		nmeteo_stations = num_param_components[cod];
	}else {
		nmeteo_stations = (long)assignation_number(flog, 191, 0, keyword, num_param, num_param_components, 1., 0);
	}

	met->st=(METEO_STATIONS *)malloc(sizeof(METEO_STATIONS));	
	if(!met->st) t_error("meteo_stations was not allocated"); 
	met->st->E=new_doublevector(nmeteo_stations);
	met->st->N=new_doublevector(nmeteo_stations);
	met->st->lat=new_doublevector(nmeteo_stations);
	met->st->lon=new_doublevector(nmeteo_stations);
	met->st->Z=new_doublevector(nmeteo_stations);
	met->st->sky=new_doublevector(nmeteo_stations);
	met->st->ST=new_doublevector(nmeteo_stations);
	met->st->Vheight=new_doublevector(nmeteo_stations);
	met->st->Theight=new_doublevector(nmeteo_stations);
	
	i=1;
	met->st->E->co[i] = assignation_number(flog, 192, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	met->st->N->co[i] = assignation_number(flog, 193, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	met->st->lat->co[i] = assignation_number(flog, 194, i-1, keyword, num_param, num_param_components, par->latitude, 0);
	met->st->lon->co[i] = assignation_number(flog, 195, i-1, keyword, num_param, num_param_components, par->longitude, 0);
	met->st->Z->co[i] = assignation_number(flog, 196, i-1, keyword, num_param, num_param_components, 0., 0);
	met->st->sky->co[i] = assignation_number(flog, 197, i-1, keyword, num_param, num_param_components, 1., 0);
	met->st->ST->co[i] = assignation_number(flog, 198, i-1, keyword, num_param, num_param_components, par->ST, 0);
	
	for (i=2; i<=nmeteo_stations; i++) {
		met->st->E->co[i] = assignation_number(flog, 192, i-1, keyword, num_param, num_param_components, met->st->E->co[i-1], 0);
		met->st->N->co[i] = assignation_number(flog, 193, i-1, keyword, num_param, num_param_components, met->st->N->co[i-1], 0);
		met->st->lat->co[i] = assignation_number(flog, 194, i-1, keyword, num_param, num_param_components, met->st->lat->co[i-1], 0);
		met->st->lon->co[i] = assignation_number(flog, 195, i-1, keyword, num_param, num_param_components, met->st->lon->co[i-1], 0);
		met->st->Z->co[i] = assignation_number(flog, 196, i-1, keyword, num_param, num_param_components, met->st->Z->co[i-1], 0);
		met->st->sky->co[i] = assignation_number(flog, 197, i-1, keyword, num_param, num_param_components, met->st->sky->co[i-1], 0);
		met->st->ST->co[i] = assignation_number(flog, 198, i-1, keyword, num_param, num_param_components, met->st->ST->co[i-1], 0);
 	}
		
	a = assignation_number(flog, 199, 0, keyword, num_param, num_param_components, 10., 0);
	for (i=1; i<=nmeteo_stations; i++) {
		met->st->Vheight->co[i] = a;
	}
	
	a = assignation_number(flog, 200, 0, keyword, num_param, num_param_components, 2., 0);
	for (i=1; i<=nmeteo_stations; i++) {
		met->st->Theight->co[i] = a;
	}
	
	//lapse rates (cyclic)
	n = (long)nlstot;
	met->LRcnc = (long*)malloc(n*sizeof(long));
	met->LRcnc[ilsDate12] = 1;
	met->LRcnc[ilsTa] = num_param_components[201];
	met->LRcnc[ilsTdew] = num_param_components[202];
	met->LRcnc[ilsPrec] = num_param_components[203];
	met->LRc = (double**)malloc(n*sizeof(double*));
	for (i=0; i<nlstot; i++) {
		met->LRc[i] = (double*)malloc(met->LRcnc[i]*sizeof(double));
		for (j=0; j<met->LRcnc[i]; j++) {
			if(i==ilsDate12) met->LRc[i][j] = 0.;
			if(i==ilsTa) met->LRc[i][j] = assignation_number(flog, 201, j, keyword, num_param, num_param_components, (double)number_novalue, 0);
			if(i==ilsTdew) met->LRc[i][j] = assignation_number(flog, 202, j, keyword, num_param, num_param_components, (double)number_novalue, 0);
			if(i==ilsPrec) met->LRc[i][j] = assignation_number(flog, 203, j, keyword, num_param, num_param_components, (double)number_novalue, 0);
		}
	}
	
	par->MinIncrFactWithElev = assignation_number(flog, 204, j, keyword, num_param, num_param_components, 0.1, 0);
	par->MaxIncrFactWithElev = assignation_number(flog, 205, j, keyword, num_param, num_param_components, 4.4, 0);
		
	//output point column
	n = (long)otot;
	outputpoint = (long*)malloc(n*sizeof(long));
	
	par->default_point = (short)assignation_number(flog, 283, 0, keyword, num_param, num_param_components, 1., 0);
	
	if (par->default_point == 1) {
		
		for (i=0; i<n; i++) {
			outputpoint[i] = i;
		}
		
		noutputpoint = n;
		
	}else {
		
		for (i=0; i<n; i++) {
			outputpoint[i] = -1;
		}
		for (i=0; i<n; i++) {
			j = (long)assignation_number(flog, 206+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) outputpoint[j-1] = i;
		}
		noutputpoint = 0;
		for (i=0; i<n; i++) {
			if(outputpoint[i] > 0) noutputpoint = i+1;
		}
	}
	
	//output basin column
	n = (long)ootot;
	outputbasin = (long*)malloc(n*sizeof(long));
	
	par->default_basin = (short)assignation_number(flog, 309, 0, keyword, num_param, num_param_components, 1., 0);
	
	if (par->default_basin == 1) {
		
		for (i=0; i<n; i++) {
			outputbasin[i] = i;
		}
		
		noutputbasin = (long)ootot;
		
	}else {
		
		for (i=0; i<n; i++) {
			outputbasin[i] = -1;
		}
		for (i=0; i<n; i++) {
			j = (long)assignation_number(flog, 284+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) outputbasin[j-1] = i;
		}
		noutputbasin = 0;
		for (i=0; i<n; i++) {
			if(outputbasin[i] > 0) noutputbasin = i+1;
		}
	}
	
	//output snow column
	n = 6 + 6*par->snowlayer_max;
	outputsnow = (long*)malloc(n*sizeof(long));
	
	par->default_snow = (short)assignation_number(flog, 322, 0, keyword, num_param, num_param_components, 1., 0);
	
	if (par->default_snow == 1) {
		
		for (i=0; i<n; i++) {
			outputsnow[i] = i;
		}
		
		noutputsnow = n;
		
	}else {
		
		cod = 310;
		
		for (i=0; i<n; i++) {
			outputsnow[i] = -1;
		}
		for (i=0; i<6; i++) {
			j = (long)assignation_number(flog, cod+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) outputsnow[j-1] = i;
		}
		for (i=6; i<12; i++) {
			for (k=0; k<par->snowlayer_max; k++) {
				j = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (j>=1 && j<=n) outputsnow[j-1] = (i-4)*par->snowlayer_max+k+4;
			}
		}
		noutputsnow = 0;
		for (i=0; i<n; i++) {
			if(outputsnow[i] > 0) noutputsnow = i+1;
		}
	}
	
	//output glacier column
	n = 6 + 6*par->glaclayer_max;
	outputglac = (long*)malloc(n*sizeof(long));
	
	par->default_glac = (short)assignation_number(flog, 335, 0, keyword, num_param, num_param_components, 1., 0);
	
	if (par->default_glac == 1) {
		
		for (i=0; i<n; i++) {
			outputglac[i] = i;
		}
		
		noutputglac = n;
		
	}else {
		
		cod = 323;
		
		for (i=0; i<n; i++) {
			outputglac[i] = -1;
		}
		for (i=0; i<6; i++) {
			j = (long)assignation_number(flog, cod+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) outputglac[j-1] = i;
		}
		for (i=6; i<10; i++) {
			for (k=0; k<par->glaclayer_max; k++) {
				j = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (j>=1 && j<=n) outputglac[j-1] = (i-4)*par->glaclayer_max+k+4;
			}
		}
		noutputglac = 0;
		for (i=0; i<n; i++) {
			if(outputglac[i] > 0) noutputglac = i+1;
		}
	}
	
	//output soil column
	n = 6;
	outputsoil = (long*)malloc(n*sizeof(long));
	
	par->default_soil = (short)assignation_number(flog, 342, 0, keyword, num_param, num_param_components, 1., 0);
	
	if (par->default_soil == 1) {
		
		for (i=0; i<n; i++) {
			outputsoil[i] = i;
		}
		
		noutputsoil = n;
		
	}else {
		
		for (i=0; i<n; i++) {
			outputsoil[i] = -1;
		}
		for (i=0; i<n; i++) {
			j = (long)assignation_number(flog, 336+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) outputsoil[j-1] = i;
		}
		noutputsoil = 0;
		for (i=0; i<n; i++) {
			if(outputsoil[i] > 0) noutputsoil = i+1;
		}
	}
	
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

char **assign_string_parameter(FILE *f, long beg, long end, char **string_param, char **keyword){
	
	long i;
	char **a;
	
	a = (char**)malloc((end-beg)*sizeof(char*));
	
	for (i=0; i<end-beg; i++) {
		a[i] = assignation_string(f, i+beg, keyword, string_param);
	}
	
	return(a);
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

double assignation_number(FILE *f, long i, long j, char **keyword, double **num_param, long *num_param_components, double default_value, short code_error){
	
	double a;
	
	if (j < num_param_components[i]) {
		a = num_param[i][j];
	}else {
		a = (double)number_novalue;
	}
	
	if((long)a == number_novalue){		
		if (code_error==0) {
			a = default_value;
			fprintf(f,"WARNING: The parameter %s component %ld assigned at the default value (%e) i.e. (%f)\n", keyword[i], j+1, a, a);
		}else{
			printf("ERROR: The parameter %s component %ld not assigned\n", keyword[i], j+1);
			stop_execution();
			t_error("ERROR");
		}
	}else{
		fprintf(f,"The parameter %s component %ld assigned at %e i.e. %f\n", keyword[i], j+1, a, a);
	}
	
	return a;
}	

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

char *assignation_string(FILE *f, long i, char **keyword, char **string_param){
	
	char *a;
	long j, dimstring = strlen(string_param[i]);
	
	a = (char*)malloc((dimstring+1)*sizeof(char));
	
	for (j=0; j<dimstring; j++) {
		a[j] = string_param[i][j];
	}
	a[dimstring] = 0;
	
	fprintf(f,"%s assigned at %s\n", keyword[i], a);
	
	return(a);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_soil_parameters(char *name, char **key_header, SOIL *sl, FILE *flog){
	
	short ok;
	long i, j, k, n, nlines, nlinesprev;
	char *temp;
	double **soildata;
	DOUBLETENSOR *old_sl_par;
	
	//look if there is at least 1 soil file
	i = 0;
	ok = 0;
	nlinesprev = -1;
	
	do{
		
		temp = namefile_i_we2(name, i+1);
		
		if(existing_file_text(temp)==1){
			free(temp);
			ok = 1;
			temp = namefile_i(name, i+1);
			nlines = count_lines(temp, 33, 44);
			if (nlinesprev >= 0 && nlines != nlinesprev){
				printf("Error:: The file %s with soil paramaters has a number of layers %ld, which different from the numbers %ld of the other soil parameter files\n",temp,nlines,nlinesprev);
				printf("In GEOtop it is only possible to have the same number of layers in any soil parameter files\n");
				stop_execution();
				t_error("It is not possible to continue");
			}
			nlinesprev = nlines;
		}
			
		free(temp);
		i++;
		
	}while (ok == 0 && i < sl->pa->ndh);
	
	if (ok == 1){
		
		//save sl->pa in a new doubletensor and deallocate
		old_sl_par = new_doubletensor(sl->pa->ndh, sl->pa->nrh, sl->pa->nch);
		for (i=1; i<=sl->pa->ndh; i++) {
			for (n=1; n<=sl->pa->nrh; n++) {
				for (j=1; j<=sl->pa->nch; j++) {
					old_sl_par->co[i][n][j] = sl->pa->co[i][n][j];
				}
			}
		}
		free_doubletensor(sl->pa);
		
		//reallocate
		sl->pa = new_doubletensor(old_sl_par->ndh, old_sl_par->nrh, nlines);
		
		for (i=1; i<=sl->pa->ndh; i++) {
			
			//read files
			temp = namefile_i_we2(name, i);
			
			if(existing_file_text(temp)==1){
				free(temp);
				temp = namefile_i(name, i);
				soildata = read_txt_matrix(temp, 33, 44, key_header, nsoilprop, &nlines, flog);
			}else {
				soildata = (double**)malloc(nlines*sizeof(double*));
				for (j=0; j<nlines; j++) {
					k = (long)nsoilprop;
					soildata[j] = (double*)malloc(k*sizeof(double));
					for (n=0; n<k; n++) {
						soildata[j][n] = (double)number_absent;
					}
				}
			}

			free(temp);
			
			//assign soildata to soil->pa
			for (n=1; n<=nsoilprop; n++) {
				for (j=1; j<=sl->pa->nch; j++) { //j is the layer index
					sl->pa->co[i][n][j] = soildata[j-1][n-1];
				}
			}
			 
			//deallocate soildata
			for (j=0; j<nlines; j++) {
				free(soildata[j]);
			}
			free(soildata);
			
			//fix layer thickness
			n = jdz;
			for (j=1; j<=sl->pa->nch; j++) { //j is the layer index
				if ((long)sl->pa->co[i][n][j] != number_novalue && (long)sl->pa->co[i][n][j] != number_absent) {
					if ( i > 1 && fabs( sl->pa->co[i][n][j] - sl->pa->co[i-1][n][j] ) > 1.E-5 )  {
						printf("For soil type %ld it has been given a set of soil layer thicknesses different from the other ones\n",i);
						printf("In GEOtop it is only possible to have the soil layer discretization in any soil parameter files\n");
						stop_execution();
						t_error("It is not possible to continue");
					}
				}else if (i == 1) {
					if (j <= old_sl_par->nch) {
						sl->pa->co[i][n][j] = old_sl_par->co[i][n][j];
					}else {
						sl->pa->co[i][n][j] = sl->pa->co[i][n][j-1];
					}
				}else {
					sl->pa->co[i][n][j] = sl->pa->co[i-1][n][j];
				}
			}
			
			//all other variables
			for (n=1; n<=nsoilprop; n++) {
				if (n != jdz) {
					for (j=1; j<=sl->pa->nch; j++) { //j is the layer index
						if ((long)sl->pa->co[i][n][j] == number_novalue || (long)sl->pa->co[i][n][j] == number_absent) {
							if (j <= old_sl_par->nch) {
								sl->pa->co[i][n][j] = old_sl_par->co[i][n][j];
							}else {
								sl->pa->co[i][n][j] = sl->pa->co[i][n][j-1];
							}
						}
					}
				}
			}						
		}
		
		free_doubletensor(old_sl_par);
			
	}
	
	//write on the screen the soil paramater
	fprintf(flog,"\n");
	k = (long)nmet;
	fprintf(flog,"Soil Layers: %ld\n",sl->pa->nch);
	for (i=1; i<=sl->pa->ndh; i++) {
		fprintf(flog,"-> Soil Type: %ld\n",i);
		for (n=1; n<=nsoilprop; n++) {
			fprintf(flog,"%s: ",keywords_char[k+n-1]);
			for (j=1; j<=sl->pa->nch; j++) {
				fprintf(flog,"%f(%.2e)",sl->pa->co[i][n][j],sl->pa->co[i][n][j]);
				if(j<sl->pa->nch)fprintf(flog,", ");
			}
			fprintf(flog,"\n");
		}
	}
	
	//bedrock
	old_sl_par = new_doubletensor(1, sl->pa_bed->nrh, sl->pa_bed->nch);
	for (n=1; n<=sl->pa_bed->nrh; n++) {
		for (j=1; j<=sl->pa_bed->nch; j++) {
			old_sl_par->co[1][n][j] = sl->pa_bed->co[1][n][j];
		}
	}
	free_doubletensor(sl->pa_bed);
	
	sl->pa_bed = new_doubletensor(sl->pa->ndh, sl->pa->nrh, sl->pa->nch);
	for (i=1; i<=sl->pa_bed->ndh; i++) {
		for (n=1; n<=sl->pa_bed->nrh; n++) {
			if (n == jdz) {
				for (j=1; j<=sl->pa_bed->nch; j++) {
					sl->pa_bed->co[i][n][j] = sl->pa->co[1][n][j];
				}
			}else {
				for (j=1; j<=sl->pa_bed->nch; j++) {
					if (j <= old_sl_par->nch) {
						sl->pa_bed->co[i][n][j] = old_sl_par->co[1][n][j];
					}else {
						sl->pa_bed->co[i][n][j] = sl->pa_bed->co[i][n][j-1];
					}
				}
				for (j=1; j<=sl->pa_bed->nch; j++) {
					if ( (long)sl->pa_bed->co[i][n][j] == number_novalue ) sl->pa_bed->co[i][n][j] = sl->pa->co[i][n][j];
				}
			}
		}
	}
	free_doubletensor(old_sl_par);
	
	fprintf(flog,"\n");
	k = (long)nmet;
	fprintf(flog,"Soil Bedrock Layers: %ld\n",sl->pa->nch);
	for (i=1; i<=sl->pa_bed->ndh; i++) {
		fprintf(flog,"-> Soil Type: %ld\n",i);
		for (n=1; n<=nsoilprop; n++) {
			fprintf(flog,"%s: ",keywords_char[k+n-1]);
			for (j=1; j<=sl->pa->nch; j++) {
				fprintf(flog,"%f(%.2e)",sl->pa_bed->co[i][n][j],sl->pa_bed->co[i][n][j]);
				if(j<sl->pa->nch)fprintf(flog,", ");
			}
			fprintf(flog,"\n");
		}
	}
	fprintf(flog,"\n");
	
	return 1;
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		

short read_point_file(char *name, char **key_header, PAR *par, FILE *flog){
	
	DOUBLEMATRIX *chkpt2;
	double **points;
	long nlines, n, j;
	char *temp;
		
	if(existing_file_text(name)==1){	
		temp = join_strings(name, textfile);
		points = read_txt_matrix(temp, 33, 44, key_header, par->chkpt->nch, &nlines, flog);
		free(temp);
				
		chkpt2 = new_doublematrix(par->chkpt->nrh, par->chkpt->nch);
		copy_doublematrix(par->chkpt, chkpt2);
		free_doublematrix(par->chkpt);
		
		par->chkpt = new_doublematrix(nlines, chkpt2->nch);
		for (n=1; n<=nlines; n++) {
			for (j=1; j<=chkpt2->nch; j++) {
				par->chkpt->co[n][j] = points[n-1][j-1];
				if ( (long)par->chkpt->co[n][j] == number_novalue || (long)par->chkpt->co[n][j] == number_absent ) {
					if ( n <= chkpt2->nrh ) {
						par->chkpt->co[n][j] = chkpt2->co[n][j];
					}else {
						par->chkpt->co[n][j] = chkpt2->co[chkpt2->nrh][j];
					}
				}
			}
			
			if ( par->point_sim != 1 ){
				if( (long)par->chkpt->co[n][ptX] == number_novalue || (long)par->chkpt->co[n][ptY] == number_novalue ) {
					par->state_pixel = 0;
				}
			}
			
			free(points[n-1]);
		}
		free_doublematrix(chkpt2);
		free(points);
	}
	
	if (par->point_sim != 1){
		for (n=1; n<=par->chkpt->nrh; n++) {
			if ( (long)par->chkpt->co[n][ptX] == number_novalue || (long)par->chkpt->co[n][ptY] == number_novalue) {
				fprintf(flog,"Warning: The points to plot specific results are not completely specified\n");
				fprintf(flog,"Output for single point output is deactivated.\n");
				printf("Warning: The points to plot specific results are not completely specified\n");
				printf("Output for single point output is deactivated.\n");
				par->state_pixel = 0;
			}
		}
	}
	
	
	return 1;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		
	
short read_meteostations_file(LONGVECTOR *i, METEO_STATIONS *S, char *name, char **key_header, FILE *flog){
	
	double **M;
	long nlines, n, j, k;
	char *temp;
		
	if(existing_file_text(name)==1){	
		temp = join_strings(name, textfile);
		M = read_txt_matrix(temp, 33, 44, key_header, 8, &nlines, flog);
		free(temp);
				
		for (j=1; j<=i->nh; j++) {
			for (n=1; n<=nlines; n++) {
				if ((long)M[n-1][0] == i->co[j]) {
					for (k=1; k<8; k++) {
						if ((long)M[n-1][k] != number_novalue && (long)M[n-1][k] != number_absent) {
							if (k==1) {
								S->E->co[j] = M[n-1][k];
							}else if (k==2) {
								S->N->co[j] = M[n-1][k];
							}else if (k==3) {
								S->lat->co[j] = M[n-1][k];
							}else if (k==4) {
								S->lon->co[j] = M[n-1][k];
							}else if (k==5) {
								S->Z->co[j] = M[n-1][k];
							}else if (k==6) {
								S->sky->co[j] = M[n-1][k];
							}else if (k==7) {
								S->ST->co[j] = M[n-1][k];
							}
						}
					}
				}					
			}
		}
		
		for (n=1; n<=nlines; n++) {
			free(M[n-1]);
		}
		
		free(M);
		
		return 1;
	}else {
		return 0;
	}
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		
