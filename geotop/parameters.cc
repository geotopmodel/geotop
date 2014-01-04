/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 1.225-15
 
 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "parameters.h"
#include "constants.h"
#include <iomanip>
#include <inputKeywords.h>

using namespace std;
using namespace boost::assign ;

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

/** @brief get a double parameter from the configuration file and return it, applying default value.
  * @param[in] pConfigStore a pointer to a previously initialized configuration store
  * @param[in] pName the keyword of the parameter to read
  * @param[in] pDefaultValue a default value to be applied if the value of the parameter is geotop::input::gDoubleNoValue
  * @param[in] pAllowNoValue if true the returning of no value is permitted, by default is false
  * @return the parameter value. If the keyword does not exists or if the return value is
  *     geotop::input::gDoubleNoValue (and pAllowNoValue == false)
  *     the function will cause an exit from the program by calling the t_error function
  */
static double getDoubleValueWithDefault(const boost::shared_ptr<geotop::input::ConfigStore> pConfigStore, const std::string pName, const double pDefaultValue, const bool pAllowNoValue = false){
    
	double lValue;
    bool lGetResult = pConfigStore->get(pName, lValue) ;
    
    if(not lGetResult) {
        t_error(std::string("Fatal Error: unable to get parameter: ") + pName);
    }
    
    if(lValue == geotop::input::gDoubleNoValue) {
        lValue = pDefaultValue ;
    }
    
    if(lValue == geotop::input::gDoubleNoValue && not pAllowNoValue) {
        t_error(std::string("Fatal Error: value not assigned: ") + pName);
    }
    
    return lValue ;
}

/** @brief get a double vector parameter from the configuration file and return a modified version (extend the size, apply default values, etc.)
  * @param[in] pConfigStore a pointer to a previously initialized configuration store
  * @param[in] pName the keyword of the parameter to read
  * @param[in] pDefaultValue a default value to be applied if the value of the parameter is geotop::input::gDoubleNoValue
  * @param[in] pUsePrevElement if true, use the element repeat the last element instead of the default value, use the default value only
  *     if the latest value is geotop::input::gDoubleNoValue. If false use the default value for each element.
  * @param[in] pLength the expected length of the output vector
  * @param[in] pAllowNoValue if true the returning of a vector containing no value is permitted, by default is false
  * @return the parameter with the array of double. If the keyword does not exists or if the return value is
  *     geotop::input::gDoubleNoValue (and pAllowNoValue == false)
  *     the function will cause an exit from the program by calling the t_error function
  */
static std::vector<double> getDoubleVectorValueWithDefault(const boost::shared_ptr<geotop::input::ConfigStore> pConfigStore, const std::string pName, const double pDefaultValue, const bool pUsePrevElement, const size_t pLength, const bool pAllowNoValue = false){

    std::vector<double> lValue;
    bool lGetResult = pConfigStore->get(pName, lValue) ;
    
    size_t lDesiredLength = pLength ;
    if ( pLength == 0 ) {
        lDesiredLength = lValue.size() ;
    }
    
    if(not lGetResult) {
        t_error(std::string("Fatal Error: unable to get parameter: ") + pName);
    }

    size_t lLength = lValue.size() ;
    
    for (size_t i = 0; i < lDesiredLength ; i++) {
        double lElementValue = geotop::input::gDoubleNoValue ;
        
        if(i < lLength) {
            lElementValue = lValue[i] ;
        } else {
            lValue.push_back( geotop::input::gDoubleNoValue ) ;
        }

        if ( lElementValue == geotop::input::gDoubleNoValue ) {
            if( i > 0 ) {
                if(pUsePrevElement)
                    lElementValue = lValue[i-1] ;
                else
                    lElementValue = pDefaultValue ;
            } else {
                lElementValue = pDefaultValue ;
            }
        }

        if ( not pAllowNoValue && lElementValue == geotop::input::gDoubleNoValue ) {
            t_error(std::string("Fatal Error: array value not assigned: ") + pName);
        } else {
            lValue[i] = lElementValue ;
        }
    }
    
    return lValue ;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_inpts_par(Par *par, Land *land, Times *times, Soil *sl, Meteo *met, InitTools *itools, std::string filename, FILE *flog){

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
	
    std::string temp;
    std::string path_rec_files;
	
	char **keywords_num_lower_case, **keywords_char_lower_case;
	
	long beg=0, end=0;

    boost::shared_ptr<geotop::input::ConfigStore> lConfigStore = geotop::input::ConfigStoreSingletonFactory::getInstance() ;

#ifdef STAGED_FOR_REMOVING
	//convert keyword listed on top of the file in lower case
	n = (long)num_par_number;
	keywords_num_lower_case = (char**)malloc(n*sizeof(char*));
	for (i=0; i<n; i++) {
		//printf("%ld,%s\n",i,keywords_num[i]);
		keywords_num_lower_case[i] = assign_string(keywords_num[i].c_str());
		convert_string_in_lower_case(keywords_num_lower_case[i]);
	}
#endif

	n = (long)num_par_char;
	keywords_char_lower_case = (char**)malloc(n*sizeof(char*));
	for (i=0; i<n; i++) {
		//printf("%ld,%s\n",i,keywords_char[i]);
		keywords_char_lower_case[i] = assign_string(keywords_char[i].c_str());
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
	f = t_fopen(filename.c_str(), "r");
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
	f = fopen(filename.c_str(), "r");
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
			num_param[i][0] = geotop::input::gDoubleNoValue;
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
#ifdef STAGED_FOR_REMOVING
	assign_numeric_parameters(par, land, times, sl, met, itools, num_param, num_param_components, keywords_num, flog);
#else
	assign_numeric_parameters(par, land, times, sl, met, itools, num_param, num_param_components, flog);    
#endif

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

	//assign string parameter
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
	files = assign_string_parameter_v(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += otot;
	hpnt = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += ootot;
	hbsn = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 10;
	hsnw = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 10;	
	hglc = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 6;	
	hsl = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += ptTOT;		
	itools->point_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += nlstot;
	itools->lapserates_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	beg = end;
	end += 8;
	itools->meteostations_col_names = assign_string_parameter(flog, beg, end, string_param, keywords_char);
	
	beg = end;
	end += 1;
	temp = assignation_string(flog, beg, keywords_char, string_param);
	if (strcmp(string_novalue,temp.c_str()) == 0) {
		temp = "_SUCCESSFUL_RUN";
	}
	SuccessfulRunFile = WORKING_DIRECTORY + temp ;
			
	beg = end;
	end += 1;
	temp = assignation_string(flog, beg, keywords_char, string_param);
	if (strcmp(string_novalue,temp.c_str()) == 0) {
		temp = "_FAILED_RUN";
	}	
	FailedRunFile = WORKING_DIRECTORY + temp;
	
	beg = end;
	end += 1;
	path_rec_files = assignation_string(flog, beg, keywords_char, string_param);//path of recovery files
			
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
		if(strcmp(hpnt[j], string_novalue) == 0){
			free(hpnt[j]);
			if(j==odate12){ hpnt[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==oJDfrom0){ hpnt[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==odaysfromstart){ hpnt[j] = assign_string("TimeFromStart[days]"); 
			}else if(j==operiod){ hpnt[j] = assign_string("Simulation_Period");
			}else if(j==orun){ hpnt[j] = assign_string("Run");	
			}else if(j==opoint){ hpnt[j] = assign_string("IDpoint"); 
			}else if(j==osnowover){ hpnt[j] = assign_string("Psnow_over_canopy[mm]");     	
			}else if(j==orainover){ hpnt[j] = assign_string("Prain_over_canopy[mm]"); 		
			}else if(j==oprecsnow){ hpnt[j] = assign_string("Psnow_under_canopy[mm]"); 
			}else if(j==oprecrain){ hpnt[j] = assign_string("Prain_under_canopy[mm]"); 		
			}else if(j==orainonsnow){ hpnt[j] = assign_string("Prain_rain_on_snow[mm]");
			}else if(j==oV){ hpnt[j] = assign_string("Wind_speed[m/s]");      	   
			}else if(j==oVdir){ hpnt[j] = assign_string("Wind_direction[deg]");  
			}else if(j==oRH){ hpnt[j] = assign_string("Relative_Humidity[-]");    
			}else if(j==oPa){ hpnt[j] = assign_string("Pressure[mbar]");    
			}else if(j==oTa){ hpnt[j] = assign_string("Tair[C]");    
			}else if(j==oTdew){ hpnt[j] = assign_string("Tdew[C]");  
			}else if(j==oTg){ hpnt[j] = assign_string("Tsurface[C]");    
			}else if(j==oTv){ hpnt[j] = assign_string("Tvegetation[C]");    
			}else if(j==oTs){ hpnt[j] = assign_string("Tcanopyair[C]");    
			}else if(j==oEB){ hpnt[j] = assign_string("Surface_Energy_balance[W/m2]");    
			}else if(j==oG){ hpnt[j] = assign_string("Soil_heat_flux[W/m2]");     
			}else if(j==oSWin){ hpnt[j] = assign_string("SWin[W/m2]");  
			}else if(j==oSWb){ hpnt[j] = assign_string("SWbeam[W/m2]");   
			}else if(j==oSWd){ hpnt[j] = assign_string("SWdiff[W/m2]");  
			}else if(j==oLWin){ hpnt[j] = assign_string("LWin[W/m2]"); 
			}else if(j==ominLWin){ hpnt[j] = assign_string("LWin_min[W/m2]"); 
			}else if(j==omaxLWin){ hpnt[j] = assign_string("LWin_max[W/m2]");
			}else if(j==oSW){ hpnt[j] = assign_string("SWnet[W/m2]");     
			}else if(j==oLW){ hpnt[j] = assign_string("LWnet[W/m2]");     
			}else if(j==oH){ hpnt[j] = assign_string("H[W/m2]");      
			}else if(j==oLE){ hpnt[j] = assign_string("LE[W/m2]");     
			}else if(j==ofc){ hpnt[j] = assign_string("Canopy_fraction[-]");     
			}else if(j==oLSAI){ hpnt[j] = assign_string("LSAI[m2/m2]");   
			}else if(j==oz0v){ hpnt[j] = assign_string("z0veg[m]");    
			}else if(j==od0v){ hpnt[j] = assign_string("d0veg[m]");    
			}else if(j==oEcan){ hpnt[j] = assign_string("Estored_canopy[W/m2]");   
			}else if(j==oSWv){ hpnt[j] = assign_string("SWv[W/m2]");    
			}else if(j==oLWv){ hpnt[j] = assign_string("LWv[W/m2]");    
			}else if(j==oHv){ hpnt[j] = assign_string("Hv[W/m2]");     
			}else if(j==oLEv){ hpnt[j] = assign_string("LEv[W/m2]");    
			}else if(j==oHg0){ hpnt[j] = assign_string("Hg_unveg[W/m2]");    
			}else if(j==oLEg0){ hpnt[j] = assign_string("LEg_unveg[W/m2]");   
			}else if(j==oHg1){ hpnt[j] = assign_string("Hg_veg[W/m2]");    
			}else if(j==oLEg1){ hpnt[j] = assign_string("LEg_veg[W/m2]");   
			}else if(j==oevapsur){ hpnt[j] = assign_string("Evap_surface[mm]");  
			}else if(j==otrasp){ hpnt[j] = assign_string("Trasp_canopy[mm]");    
			}else if(j==owcan_rain){ hpnt[j] = assign_string("Water_on_canopy[mm]");
			}else if(j==owcan_snow){ hpnt[j] = assign_string("Snow_on_canopy[mm]");
			}else if(j==oQv){ hpnt[j] = assign_string("Qvegetation[-]");   
			}else if(j==oQg){ hpnt[j] = assign_string("Qsurface[-]");   
			}else if(j==oQa){ hpnt[j] = assign_string("Qair[-]");   
			}else if(j==oQs){ hpnt[j] = assign_string("Qcanopyair[-]");   
			}else if(j==oLobuk){ hpnt[j] = assign_string("LObukhov[m]");
			}else if(j==oLobukcan){ hpnt[j] = assign_string("LObukhovcanopy[m]");
			}else if(j==outop){ hpnt[j] = assign_string("Wind_speed_top_canopy[m/s]");    
			}else if(j==odecay){ hpnt[j] = assign_string("Decay_of_K_in_canopy[-]");   
			}else if(j==oSWup){ hpnt[j] = assign_string("SWup[W/m2]");   
			}else if(j==oLWup){ hpnt[j] = assign_string("LWup[W/m2]");   
			}else if(j==oHup){ hpnt[j] = assign_string("Hup[W/m2]");    
			}else if(j==oLEup){ hpnt[j] = assign_string("LEup[W/m2]");   
			}else if(j==osnowdepth){ hpnt[j] = assign_string("snow_depth[mm]"); 
			}else if(j==oSWE){ hpnt[j] = assign_string("snow_water_equivalent[mm]"); 
			}else if(j==osnowdens){ hpnt[j] = assign_string("snow_density[kg/m3]"); 
			}else if(j==osnowT){ hpnt[j] = assign_string("snow_temperature[C]"); 	
			}else if(j==omrsnow){ hpnt[j] = assign_string("snow_melted[mm]"); 
			}else if(j==osrsnow){ hpnt[j] = assign_string("snow_subl[mm]"); 
			}else if(j==oblowingsnowtrans){ hpnt[j] = assign_string("snow_blown_away[mm]"); 
			}else if(j==oblowingsnowsubl){ hpnt[j] = assign_string("snow_subl_while_blown[mm]");			
			}else if(j==oglacdepth){ hpnt[j] = assign_string("glac_depth[mm]"); 
			}else if(j==oGWE){ hpnt[j] = assign_string("glac_water_equivalent[mm]"); 
			}else if(j==oglacdens){ hpnt[j] = assign_string("glac_density[kg/m3]"); 
			}else if(j==oglacT){ hpnt[j] = assign_string("glac_temperature[C]"); 
			}else if(j==omrglac){ hpnt[j] = assign_string("glac_melted[mm]"); 
			}else if(j==osrglac){ hpnt[j] = assign_string("glac_subl[mm]"); 
			}else if(j==othawedup){ hpnt[j] = assign_string("lowest_thawed_soil_depth[mm]"); 
			}else if(j==othaweddw){ hpnt[j] = assign_string("highest_thawed_soil_depth[mm]"); 
			}else if(j==owtableup){ hpnt[j] = assign_string("lowest_water_table_depth[mm]"); 
			}else if(j==owtabledw){ hpnt[j] = assign_string("highest_water_table_depth[mm]"); 
			}
		}
	}
	
	//HeaderBasinFile
	for (j=0; j<ootot; j++) {
		if(strcmp(hbsn[j], string_novalue) == 0){
			free(hbsn[j]);
			if(j==oodate12){ hbsn[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==ooJDfrom0){ hbsn[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==oodaysfromstart){ hbsn[j] = assign_string("TimeFromStart[days]");    
			}else if(j==ooperiod){ hbsn[j] = assign_string("Simulation_Period");
			}else if(j==oorun){ hbsn[j] = assign_string("Run");	
			}else if(j==ooprecrain){ hbsn[j] = assign_string("Prain_below_canopy[mm]");     	
			}else if(j==ooprecsnow){ hbsn[j] = assign_string("Psnow_below_canopy[mm]");     	
			}else if(j==oorainover){ hbsn[j] = assign_string("Prain_above_canopy[mm]");     	
			}else if(j==oosnowover){ hbsn[j] = assign_string("Prain_above_canopy[mm]"); 
			}else if(j==oopnet){ hbsn[j] = assign_string("Pnet[mm]"); 			
			}else if(j==ooTa){ hbsn[j] = assign_string("Tair[C]");     	
			}else if(j==ooTg){ hbsn[j] = assign_string("Tsurface[C]");     	
			}else if(j==ooTv){ hbsn[j] = assign_string("Tvegetation[C]");     	
			}else if(j==ooevapsur){ hbsn[j] = assign_string("Evap_surface[mm]");     	
			}else if(j==ootrasp){ hbsn[j] = assign_string("Transpiration_canopy[mm]");     	
			}else if(j==ooLE){ hbsn[j] = assign_string("LE[W/m2]");     	
			}else if(j==ooH){ hbsn[j] = assign_string("H[W/m2]");     	
			}else if(j==ooSW){ hbsn[j] = assign_string("SW[W/m2]");     	
			}else if(j==ooLW){ hbsn[j] = assign_string("LW[W/m2]");     	
			}else if(j==ooLEv){ hbsn[j] = assign_string("LEv[W/m2]");     	
			}else if(j==ooHv){ hbsn[j] = assign_string("Hv[W/m2]");     	
			}else if(j==ooSWv){ hbsn[j] = assign_string("SWv[W/m2]");     	
			}else if(j==ooLWv){ hbsn[j] = assign_string("LWv[W/m2]");     	
			}else if(j==ooSWin){ hbsn[j] = assign_string("SWin[W/m2]");     	
			}else if(j==ooLWin){ hbsn[j] = assign_string("LWin[W/m2]");     	
			}else if(j==oomasserror){ hbsn[j] = assign_string("Mass_balance_error[mm]");     	
			}else if(j==ootimestep){ hbsn[j] = assign_string("Mean_Time_Step[s]");     					
			}
		}
	}
	
	//HeaderSnowFile
	for (j=0; j<10; j++) {
		if(strcmp(hsnw[j], string_novalue) == 0){
			free(hsnw[j]);
			if(j==0){ hsnw[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==1){ hsnw[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==2){ hsnw[j] = assign_string("TimeFromStart[days]");  
			}else if(j==3){ hsnw[j] = assign_string("Simulation_Period");
			}else if(j==4){ hsnw[j] = assign_string("Run");					
			}else if(j==5){ hsnw[j] = assign_string("IDpoint"); 
			}
		}
	}
	
	//HeaderGlacierFile
	for (j=0; j<10; j++) {
		if(strcmp(hglc[j], string_novalue) == 0){
			
			free(hglc[j]);
			if(j==0){ hglc[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==1){ hglc[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==2){ hglc[j] = assign_string("TimeFromStart[days]");   
			}else if(j==3){ hglc[j] = assign_string("Simulation_Period");
			}else if(j==4){ hglc[j] = assign_string("Run");					
			}else if(j==5){ hglc[j] = assign_string("IDpoint"); 
			}else if(j==6){ hglc[j] = assign_string("Temperature[C]");     	
			}else if(j==7){ hglc[j] = assign_string("wice[kg/m2]");     	
			}else if(j==8){ hglc[j] = assign_string("wliq[kg/m2]");     	
			}else if(j==9){ hglc[j] = assign_string("Dz[mm]");
			}
		}
	}
	
	
	//HeaderGlacierFile
	for (j=0; j<6; j++) {
		if(strcmp(hsl[j], string_novalue) == 0){
			free(hsl[j]);
			if(j==0){ hsl[j] = assign_string("Date12[DDMMYYYYhhmm]");    		
			}else if(j==1){ hsl[j] = assign_string("JulianDayFromYear0[days]");   		
			}else if(j==2){ hsl[j] = assign_string("TimeFromStart[days]"); 
			}else if(j==3){ hsl[j] = assign_string("Simulation_Period");
			}else if(j==4){ hsl[j] = assign_string("Run");									
			}else if(j==5){ hsl[j] = assign_string("IDpoint"); 
			}
		}	
	}	
	
	//Recovery Files
	for (j=rpsi; j<=rsux; j++) {
		if(files[j] == string_novalue){
			if (j==rpsi) { files[j] = assign_string("SoilPressure");
			}else if (j==riceg) { files[j] = assign_string("SoilIceContent");
			}else if (j==rTg) { files[j] = assign_string("SoilTemperature");
			}else if (j==rDzs) { files[j] = assign_string("SnowThickness");
			}else if (j==rwls) { files[j] = assign_string("SnowLiqWaterContent");
			}else if (j==rwis) { files[j] = assign_string("SnowIceContent");
			}else if (j==rTs) { files[j] = assign_string("SnowTemperature");
			}else if (j==rDzi) { files[j] = assign_string("GlacThickness");
			}else if (j==rwli) { files[j] = assign_string("GlacLiqWaterContent");
			}else if (j==rwii) { files[j] = assign_string("GlacIceContent");
			}else if (j==rTi) { files[j] = assign_string("GlacTemperature");
			}else if (j==rns) { files[j] = assign_string("SnowLayersNumber");
			}else if (j==rni) { files[j] = assign_string("GlacLayersNumber");
			}else if (j==rsnag) { files[j] = assign_string("SnowAge");
			}else if (j==rwcrn) { files[j] = assign_string("RainOnCanopy");
			}else if (j==rwcsn) { files[j] = assign_string("SnowOnCanopy");
			}else if (j==rTv) { files[j] = assign_string("VegTemperature");
			}else if (j==rpsich) { files[j] = assign_string("SoilChannelPressure");
			}else if (j==ricegch) { files[j] = assign_string("SoilChannelIceContent");
			}else if (j==rTgch) { files[j] = assign_string("SoilChannelTemperature");
			}else if (j==rTrun) { files[j] = assign_string("RunMeanSoilTemperature");
			}else if (j==rwrun) { files[j] = assign_string("RunMeanSoilTotWater");
			}else if (j==rdUrun) { files[j] = assign_string("RunSoilInternalEnergy");
			}else if (j==rSWErun) { files[j] = assign_string("RunMeanSWE");
			}else if (j==rTmaxrun) { files[j] = assign_string("RunMaxSoilTemperature");
			}else if (j==rTminrun) { files[j] = assign_string("RunMinSoilTemperature");
			}else if (j==rwmaxrun) { files[j] = assign_string("RunMaxSoilTotWater");
			}else if (j==rwminrun) { files[j] = assign_string("RunMinSoilTotWater");
			}else if (j==rtime) { files[j] = assign_string("RecoveryTime");		
			}else if (j==rsux) { files[j] = assign_string("SuccessfulRecovery");				
			}
		}
	}

	//add path to recovery files
	if(path_rec_files != string_novalue){
		temp = path_rec_files;
		path_rec_files = temp + std::string("/");
		for (i=rpsi; i<=rsux; i++) {
			temp = files[i];
			files[i] = path_rec_files + temp;
		}
	}

	//add working path to the file name
	for (i=0; i<nfiles; i++) {
		if (files[i] != string_novalue) {
			temp = files[i];
			files[i] = WORKING_DIRECTORY + std::string(temp);
		}
	}

	
	return 1;
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
#ifdef STAGED_FOR_REMOVING
void assign_numeric_parameters(Par *par, Land *land, Times *times, Soil *sl, Meteo *met, InitTools *itools, double **num_param, long *num_param_components, std::string keyword[], FILE *flog)
#else
void assign_numeric_parameters(Par *par, Land *land, Times *times, Soil *sl, Meteo *met, InitTools *itools, double **num_param, long *num_param_components, FILE *flog)
#endif
{
	short occurring;
#ifdef STAGED_FOR_REMOVING
	size_t cod, codn ;
#endif
    size_t k, n, m, nsoillayers, nmeteo_stations, npoints;
	double a;
    double minDt=1.E99;

    std::vector<double> lDoubleTempVector ;
    double lDoubleTempValue ;
    bool lConfParamGetResult ;

	fprintf(flog,"\n");
	
	par->print=0;

    boost::shared_ptr<geotop::input::ConfigStore> lConfigStore = geotop::input::ConfigStoreSingletonFactory::getInstance() ;

#ifdef STAGED_FOR_REMOVING
	//find components of times->Dt_vector
	cod = 0;
	n = (long)GTConst::max_cols_time_steps_file + 1;
	times->Dt_vector=(double *)malloc(n*sizeof(double));
	times->Dt_vector[0] = 0.;//it is the space for the date in case of time variable time step
	times->Dt_vector[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 1);
	for (size_t i=2; i<n; i++) {
		if (i <= num_param_components[cod]) {
			times->Dt_vector[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
		}else {
			times->Dt_vector[i] = geotop::input::gDoubleNoValue;
		}
	}
	for (size_t i=1; i<n; i++) {
		if((long)times->Dt_vector[i] != geotop::input::gDoubleNoValue){
			if (times->Dt_vector[i] < minDt) minDt = times->Dt_vector[i];
		}
	}
#else
    n = (long)GTConst::max_cols_time_steps_file + 1;
    std::vector<double> lTimeStepEnergyAndWater = getDoubleVectorValueWithDefault(lConfigStore, "TimeStepEnergyAndWater", geotop::input::gDoubleNoValue, false, (long)GTConst::max_cols_time_steps_file, true);
    times->Dt_vector=(double *)malloc(n*sizeof(double));
    times->Dt_vector[0] = 0.;//it is the space for the date in case of time variable time step
    for (size_t i=1; i<n; i++) {
        times->Dt_vector[i] = lTimeStepEnergyAndWater[i-1] ;
        if((long)times->Dt_vector[i] != geotop::input::gDoubleNoValue){
            if (times->Dt_vector[i] < minDt) minDt = times->Dt_vector[i];
        }
    }
#endif

#ifdef STAGED_FOR_REMOVING
	//init date
	cod = 1;
	par->init_date.resize(num_param_components[cod] + 1, 0);
	for (size_t i=1; i<par->init_date.size(); i++) {
		par->init_date[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 010119000000., 0);
		par->init_date[i] = convert_dateeur12_JDfrom0(par->init_date[i]);
	}
	
	//simulation time
	cod = 398;
	par->simulation_hours = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 1., 0);
#else
	//init date
        std::vector<double> lInitDateDDMMYYYYhhmm = getDoubleVectorValueWithDefault(lConfigStore, "InitDateDDMMYYYYhhmm", geotop::input::gDoubleNoValue, false, 0, false);
    par->init_date.resize(lInitDateDDMMYYYYhhmm.size() + 1 , 0);
	for (size_t i=1; i<par->init_date.size(); i++) {
		par->init_date[i] = lInitDateDDMMYYYYhhmm[i-1];
		par->init_date[i] = convert_dateeur12_JDfrom0(par->init_date[i]);
	}
	
	//simulation time
    par->simulation_hours = getDoubleValueWithDefault(lConfigStore, "SimulationHours", geotop::input::gDoubleNoValue, false) ;
#endif

#ifdef STAGED_FOR_REMOVING
	//end date
	cod = 2;
//	par->end_date = new_doublevector(num_param_components[cod]);
    par->end_date.resize(num_param_components[cod] + 1, 0) ;
	if (par->end_date.size() != par->init_date.size()){
		fprintf(flog, "Error:: End date has a number of components different from Init Date");
		printf("Error:: End date has a number of components different from Init Date");
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}
		
	for (size_t i=1; i<par->end_date.size(); i++) {
		par->end_date[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
		if ((long)par->end_date[i] == geotop::input::gDoubleNoValue){
			par->end_date[i] = par->init_date[i] + par->simulation_hours/24.;
		}else {
			par->end_date[i] = convert_dateeur12_JDfrom0(par->end_date[i]);
		}		
	}
#else
	//end date
    std::vector<double> lEndDateDDMMYYYYhhmm = getDoubleVectorValueWithDefault(lConfigStore, "EndDateDDMMYYYYhhmm", geotop::input::gDoubleNoValue, false, 0, false);
    
    //	par->end_date = new_doublevector(num_param_components[cod]);
    par->end_date.resize(lEndDateDDMMYYYYhhmm.size() + 1, 0) ;
	if (par->end_date.size() != par->init_date.size()){
		fprintf(flog, "Error:: End date has a number of components different from Init Date");
		printf("Error:: End date has a number of components different from Init Date");
		t_error("Fatal Error! Geotop is closed. See failing report.");
	}
    
	for (size_t i=1; i<par->end_date.size(); i++) {
		par->end_date[i] = lEndDateDDMMYYYYhhmm[i-1];
		if ((long)par->end_date[i] == geotop::input::gDoubleNoValue){
			par->end_date[i] = par->init_date[i] + par->simulation_hours/24.;
		}else {
			par->end_date[i] = convert_dateeur12_JDfrom0(par->end_date[i]);
		}
	}
#endif

    par->run_times.resize(par->init_date.size() + 1, 0);

#ifdef STAGED_FOR_REMOVING
	//run times
	cod = 3;
	par->run_times[1] = (long)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 1., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->run_times[i] = (long)assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, (double)par->run_times[i-1], 0);
	}
#else
    std::vector<double> lNumSimulationTimes = getDoubleVectorValueWithDefault(lConfigStore, "NumSimulationTimes", 1., true, par->init_date.size(), false);
    size_t lNumSimulationTimesSize = lNumSimulationTimes.size() ;
	for (size_t i=1; i<par->run_times.size(); i++) {
		par->run_times[i] = lNumSimulationTimes[i-1] ;
	}
#endif

#ifdef STAGED_FOR_REMOVING
	par->ST = assignation_number(flog, 4, 0, keyword, num_param, num_param_components, 0., 0);
#else
	par->ST = getDoubleValueWithDefault(lConfigStore, "StandardTimeSimulation", 0., false) ;
#endif

	par->Dtplot_discharge.resize(par->init_date.size() + 1, 0);

#ifdef STAGED_FOR_REMOVING
	cod = 5;
	par->Dtplot_discharge[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->Dtplot_discharge[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_discharge[i-1], 0);
	}
#else
    std::vector<double> lDtPlotDischarge = getDoubleVectorValueWithDefault(lConfigStore, "TimeStepEnergyAndWater", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->Dtplot_discharge.size(); i++) {
		par->Dtplot_discharge[i] = lDtPlotDischarge[i-1] ;
	}
#endif
    
	par->plot_discharge_with_Dt_integration.resize(par->init_date.size() + 1, 0);
	par->state_discharge = 0;
	for (size_t i=1; i<par->init_date.size(); i++) {
		par->Dtplot_discharge[i] *= 3600.;
		if(par->Dtplot_discharge[i] > 1.E-5 && par->Dtplot_discharge[i] <= minDt){
			par->plot_discharge_with_Dt_integration[i]=1;
		}else{
			par->plot_discharge_with_Dt_integration[i]=0;
		}
		if(par->Dtplot_discharge[i] > 1.E-5) par->state_discharge = 1;
	}
	
	par->Dtplot_point.resize(par->init_date.size() + 1, 0);

#ifdef STAGED_FOR_REMOVING
	cod = 6;
	par->Dtplot_point[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->Dtplot_point[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_point[i-1], 0);
	}
#else
    std::vector<double> lDtPlotPoint = getDoubleVectorValueWithDefault(lConfigStore, "DtPlotPoint", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->Dtplot_point.size(); i++) {
		par->Dtplot_point[i] = lDtPlotPoint[i-1];
	}
#endif

	par->plot_point_with_Dt_integration.resize(par->init_date.size() + 1, 0);
	par->state_pixel = 0;
	for (size_t i=1; i<par->init_date.size(); i++) {
		par->Dtplot_point[i] *= 3600.;
		if(par->Dtplot_point[i] > 1.E-5 && par->Dtplot_point[i] <= minDt){
			par->plot_point_with_Dt_integration[i]=1;
		}else{
			par->plot_point_with_Dt_integration[i]=0;
		}
		if(par->Dtplot_point[i] > 1.E-5) par->state_pixel = 1;
	}

    par->Dtplot_basin.resize(par->init_date.size() + 1, 0);

#ifdef STAGED_FOR_REMOVING
	cod = 7;
	par->Dtplot_basin[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->Dtplot_basin[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_basin[i-1], 0);
	}
#else
    std::vector<double> lDtPlotBasin = getDoubleVectorValueWithDefault(lConfigStore, "DtPlotBasin", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->Dtplot_basin.size(); i++) {
		par->Dtplot_basin[i] = lDtPlotBasin[i-1];
	}
#endif
    
	par->plot_basin_with_Dt_integration.resize(par->init_date.size() + 1, 0);
	par->state_basin = 0;
	for (size_t i=1; i<par->init_date.size(); i++) {
		par->Dtplot_basin[i] *= 3600.;
		if(par->Dtplot_basin[i] > 1.E-5 && par->Dtplot_basin[i] <= minDt){
			par->plot_basin_with_Dt_integration[i]=1;
		}else{
			par->plot_basin_with_Dt_integration[i]=0;
		}
		if(par->Dtplot_basin[i] > 1.E-5) par->state_basin = 1;
	}

#ifdef STAGED_FOR_REMOVING
	par->lowpass = (long)assignation_number(flog, 8, 0, keyword, num_param, num_param_components, 0., 0);
	par->lowpass_curvatures = (long)assignation_number(flog, 9, 0, keyword, num_param, num_param_components, 0., 0);
	par->sky = (short)assignation_number(flog, 10, 0, keyword, num_param, num_param_components, 0., 0);
	par->format_out = (short)assignation_number(flog, 11, 0, keyword, num_param, num_param_components, 3., 0);
	par->point_sim = (short)assignation_number(flog, 12, 0, keyword, num_param, num_param_components, 0., 0);
	par->recover = (short)assignation_number(flog, 13, 0, keyword, num_param, num_param_components, 0., 0);

    //land cover types
	par->n_landuses = (long)assignation_number(flog, 14, 0, keyword, num_param, num_param_components, 1., 0);
#else
    par->lowpass = (long)getDoubleValueWithDefault(lConfigStore, "NumLowPassFilterOnDemForAll", 0., false) ;
	par->lowpass_curvatures = (long)getDoubleValueWithDefault(lConfigStore, "NumLowPassFilterOnDemForCurv", 0., false) ;
	par->sky = (short)getDoubleValueWithDefault(lConfigStore, "FlagSkyViewFactor", 0., false) ;
	par->format_out = (short)getDoubleValueWithDefault(lConfigStore, "FormatOutputMaps", 3., false) ;
	par->point_sim = (short)getDoubleValueWithDefault(lConfigStore, "PointSim", 0., false) ;
	par->recover = (short)getDoubleValueWithDefault(lConfigStore, "RecoverSim", 0., false) ;
    
    //land cover types
	par->n_landuses = (long)getDoubleValueWithDefault(lConfigStore, "NumLandCoverTypes", 0., false) ;
#endif
	
	land->ty.resize(par->n_landuses + 1, nlandprop + 1, 0);

#ifdef STAGED_FOR_REMOVING
	land->ty[1][jz0] = assignation_number(flog, 15, 0, keyword, num_param, num_param_components, 10., 0);
	land->ty[1][jz0thressoil] = assignation_number(flog, 16, 0, keyword, num_param, num_param_components, land->ty[1][jz0], 0);
	land->ty[1][jHveg] = assignation_number(flog, 17, 0, keyword, num_param, num_param_components, 1000., 0);
	land->ty[1][jz0thresveg] = assignation_number(flog, 18, 0, keyword, num_param, num_param_components, land->ty[1][jHveg], 0);
	land->ty[1][jz0thresveg2] = assignation_number(flog, 19, 0, keyword, num_param, num_param_components, land->ty[1][jz0thresveg], 0);
	land->ty[1][jLSAI] = assignation_number(flog, 20, 0, keyword, num_param, num_param_components, 1., 0);
	land->ty[1][jcf] = assignation_number(flog, 21, 0, keyword, num_param, num_param_components, 0., 0);
	land->ty[1][jdecay0] = assignation_number(flog, 22, 0, keyword, num_param, num_param_components, 2.5, 0);
	land->ty[1][jexpveg] = assignation_number(flog, 23, 0, keyword, num_param, num_param_components, 1., 0);
	land->ty[1][jroot] = assignation_number(flog, 24, 0, keyword, num_param, num_param_components, 300., 0);
	land->ty[1][jrs] = assignation_number(flog, 25, 0, keyword, num_param, num_param_components, 60., 0);
	land->ty[1][jvR_vis] = assignation_number(flog, 26, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty[1][jvR_nir] = assignation_number(flog, 27, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty[1][jvT_vis] = assignation_number(flog, 28, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty[1][jvT_nir] = assignation_number(flog, 29, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty[1][jvCh] = assignation_number(flog, 30, 0, keyword, num_param, num_param_components, 0., 0);
	land->ty[1][jcd] = assignation_number(flog, 31, 0, keyword, num_param, num_param_components, 2., 0);
	land->ty[1][ja_vis_dry] = assignation_number(flog, 32, 0, keyword, num_param, num_param_components, 0.2, 0);
	land->ty[1][ja_nir_dry] = assignation_number(flog, 33, 0, keyword, num_param, num_param_components, land->ty[1][ja_vis_dry], 0);
	land->ty[1][ja_vis_sat] = assignation_number(flog, 34, 0, keyword, num_param, num_param_components, land->ty[1][ja_vis_dry], 0);
	land->ty[1][ja_nir_sat] = assignation_number(flog, 35, 0, keyword, num_param, num_param_components, land->ty[1][ja_nir_dry], 0);
	land->ty[1][jemg] = assignation_number(flog, 36, 0, keyword, num_param, num_param_components, 0.96, 0);
	land->ty[1][jcm] = assignation_number(flog, 37, 0, keyword, num_param, num_param_components, 0.5, 0);
	land->ty[1][jN] = assignation_number(flog, 38, 0, keyword, num_param, num_param_components, 0., 0);
	land->ty[1][jdv] = assignation_number(flog, 39, 0, keyword, num_param, num_param_components, 50., 0);
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilRoughness", 10., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jz0] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "ThresSnowSoilRough", land->ty[1][jz0], true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jz0thressoil] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "VegHeight", 1000., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jHveg] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "ThresSnowVegUp", land->ty[1][jHveg], true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jz0thresveg] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "ThresSnowVegDown", land->ty[1][jz0thresveg], true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jz0thresveg2] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "LSAI", 1., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jLSAI] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "CanopyFraction", 0., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jcf] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "DecayCoeffCanopy", 2.5, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jdecay0] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "VegSnowBurying", 1., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jexpveg] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "RootDepth", 300., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jroot] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MinStomatalRes", 60., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jrs] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "VegReflectVis", 0.2, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jvR_vis] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "VegReflNIR", 0.2, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jvR_nir] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "VegTransVis", 0.2, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jvT_vis] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "VegTransNIR", 0.2, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jvT_nir] = lDoubleTempVector[i-1];
    }
    
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "LeafAngles", 0., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jvCh] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "CanDensSurface", 2., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jcd] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilAlbVisDry", 0.2, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][ja_vis_dry] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilAlbNIRDry", land->ty[1][ja_vis_dry], true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][ja_nir_dry] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilAlbVisWet", land->ty[1][ja_vis_dry], true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][ja_vis_sat] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilAlbNIRWet", land->ty[1][ja_nir_dry], true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][ja_nir_sat] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilEmissiv", 0.96, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jemg] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SurFlowResLand", 0.5, true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jcm] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "RoughElemXUnitArea", 0., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jN] = lDoubleTempVector[i-1];
    }

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "RoughElemDiam", 50., true, par->n_landuses, false) ;
	for (size_t i=1; i<land->ty.getNx(); i++) {
        land->ty[i][jdv] = lDoubleTempVector[i-1];
    }
#endif

#ifdef STAGED_FOR_REMOVING
	for (size_t i=2; i<=par->n_landuses; i++) {
		for (size_t j=1; j<=nlandprop; j++) {
			land->ty[i][j] = assignation_number(flog, 15+j-1, i-1, keyword, num_param, num_param_components, land->ty[i-1][j], 0);
		}
	}
#endif

#ifdef STAGED_FOR_REMOVING
	//former block 2
	par->imp = assignation_number(flog, 40, 0, keyword, num_param, num_param_components, 7., 0);
	par->free_drainage_bottom = assignation_number(flog, 41, 0, keyword, num_param, num_param_components, 0., 0);
	par->free_drainage_lateral = assignation_number(flog, 42, 0, keyword, num_param, num_param_components, 1., 0);
	par->TolVWb = assignation_number(flog, 43, 0, keyword, num_param, num_param_components, 1.E-6, 0);
	par->RelTolVWb = GTConst::RelativeErrorRichards;
	par->MaxErrWb = 1.E99;
	par->MaxiterTol = (long)assignation_number(flog, 44, 0, keyword, num_param, num_param_components, 100., 0);
	par->TolCG = assignation_number(flog, 45, 0, keyword, num_param, num_param_components, 0.01, 0);
	par->min_lambda_wat = assignation_number(flog, 46, 0, keyword, num_param, num_param_components, 1.E-7, 0);
	par->max_times_min_lambda_wat = (long)assignation_number(flog, 47, 0, keyword, num_param, num_param_components, 0.0, 0);
	par->exit_lambda_min_wat = (short)assignation_number(flog, 48, 0, keyword, num_param, num_param_components, 1., 0);
	par->min_Dt = assignation_number(flog, 49, 0, keyword, num_param, num_param_components, 10., 0);
	par->gamma_m = assignation_number(flog, 50, 0, keyword, num_param, num_param_components, 2./3., 0);
	par->thres_hsup_1 = assignation_number(flog, 51, 0, keyword, num_param, num_param_components, 0., 0);
	par->thres_hsup_2 = assignation_number(flog, 52, 0, keyword, num_param, num_param_components, 0., 0);
	par->Ks_channel = assignation_number(flog, 53, 0, keyword, num_param, num_param_components, 20., 0);
	par->thres_hchannel = assignation_number(flog, 54, 0, keyword, num_param, num_param_components, 50., 0);
	par->w_dx = assignation_number(flog, 55, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->depr_channel = assignation_number(flog, 56, 0, keyword, num_param, num_param_components, 500., 0);
	par->max_courant_land = assignation_number(flog, 57, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->max_courant_channel = assignation_number(flog, 58, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->min_hsup_land = assignation_number(flog, 59, 0, keyword, num_param, num_param_components, 1., 0);
	par->min_hsup_channel = assignation_number(flog, 60, 0, keyword, num_param, num_param_components, 1., 0);
	par->min_dhsup_land_channel_in = assignation_number(flog, 61, 0, keyword, num_param, num_param_components, 1., 0);
	par->dtmin_sup = assignation_number(flog, 62, 0, keyword, num_param, num_param_components, 0.01, 0);
#else
	//former block 2
    
	par->imp = getDoubleValueWithDefault(lConfigStore, "FrozenSoilHydrCondReduction", geotop::input::gDoubleNoValue, false) ;
    par->free_drainage_bottom = getDoubleValueWithDefault(lConfigStore, "FreeDrainageAtBottom", geotop::input::gDoubleNoValue, false) ;
    par->free_drainage_lateral = getDoubleValueWithDefault(lConfigStore, "FreeDrainageAtLateralBorder", geotop::input::gDoubleNoValue, false) ;
    par->TolVWb = getDoubleValueWithDefault(lConfigStore, "RichardTol", geotop::input::gDoubleNoValue, false) ;

	par->RelTolVWb = GTConst::RelativeErrorRichards;
	par->MaxErrWb = 1.E99;

    par->MaxiterTol = (long) getDoubleValueWithDefault(lConfigStore, "RichardMaxIter", geotop::input::gDoubleNoValue, false) ;
	par->TolCG =  getDoubleValueWithDefault(lConfigStore, "RichardInitForc", geotop::input::gDoubleNoValue, false) ;
    par->min_lambda_wat =  getDoubleValueWithDefault(lConfigStore, "MinLambdaWater", geotop::input::gDoubleNoValue, false) ;
    par->max_times_min_lambda_wat = (long)getDoubleValueWithDefault(lConfigStore, "MaxTimesMinLambdaWater", geotop::input::gDoubleNoValue, false) ;
    par->exit_lambda_min_wat = (short)getDoubleValueWithDefault(lConfigStore, "ExitMinLambdaWater", geotop::input::gDoubleNoValue, false) ;
    par->min_Dt = getDoubleValueWithDefault(lConfigStore, "MinTimeStep", geotop::input::gDoubleNoValue, false) ;
    par->gamma_m = getDoubleValueWithDefault(lConfigStore, "SurFlowResExp", geotop::input::gDoubleNoValue, false) ;
    par->thres_hsup_1 = getDoubleValueWithDefault(lConfigStore, "ThresWaterDepthLandInf", geotop::input::gDoubleNoValue, false) ;
    par->thres_hsup_2 = getDoubleValueWithDefault(lConfigStore, "ThresWaterDepthLandSup", geotop::input::gDoubleNoValue, false) ;
    par->Ks_channel = getDoubleValueWithDefault(lConfigStore, "SurFlowResChannel", geotop::input::gDoubleNoValue, false) ;
    par->thres_hchannel = getDoubleValueWithDefault(lConfigStore, "ThresWaterDepthChannel", geotop::input::gDoubleNoValue, false) ;
    par->w_dx = getDoubleValueWithDefault(lConfigStore, "RatioChannelWidthPixelWidth", geotop::input::gDoubleNoValue, false) ;
    par->depr_channel = getDoubleValueWithDefault(lConfigStore, "ChannelDepression", geotop::input::gDoubleNoValue, false) ;
    par->max_courant_land = getDoubleValueWithDefault(lConfigStore, "MaxCourantSupFlowLand", geotop::input::gDoubleNoValue, false) ;
    par->max_courant_channel = getDoubleValueWithDefault(lConfigStore, "MaxCourantSupFlowChannel", geotop::input::gDoubleNoValue, false) ;
    par->min_hsup_land = getDoubleValueWithDefault(lConfigStore, "MinSupWaterDepthLand", geotop::input::gDoubleNoValue, false) ;
    par->min_hsup_channel = getDoubleValueWithDefault(lConfigStore, "MinSupWaterDepthChannel", geotop::input::gDoubleNoValue, false) ;
    par->min_dhsup_land_channel_in = getDoubleValueWithDefault(lConfigStore, "MinDiffSupWaterDepthLandChannel", geotop::input::gDoubleNoValue, false) ;
	par->dtmin_sup = getDoubleValueWithDefault(lConfigStore, "MinTimeStepSupFlow", geotop::input::gDoubleNoValue, false) ;
    
#endif

#ifdef STAGED_FOR_REMOVING
	//former block 3
	par->latitude = assignation_number(flog, 63, 0, keyword, num_param, num_param_components, 45., 0);
	par->longitude = assignation_number(flog, 64, 0, keyword, num_param, num_param_components, 0., 0);
	par->Vmin = assignation_number(flog, 65, 0, keyword, num_param, num_param_components, 0.5, 0);
	par->RHmin = assignation_number(flog, 66, 0, keyword, num_param, num_param_components, 10., 0);
	par->alpha_snow = assignation_number(flog, 67, 0, keyword, num_param, num_param_components, 1.E5, 0);
	par->nsurface = (long)assignation_number(flog, 68, 0, keyword, num_param, num_param_components, 0., 0);
	par->tol_energy = assignation_number(flog, 69, 0, keyword, num_param, num_param_components, 1.E-4, 0);
	par->maxiter_energy = (long)assignation_number(flog, 70, 0, keyword, num_param, num_param_components, 500., 0);
	par->min_lambda_en = assignation_number(flog, 71, 0, keyword, num_param, num_param_components, 1.E-5, 0);
	par->max_times_min_lambda_en = (long)assignation_number(flog, 72, 0, keyword, num_param, num_param_components, 0.0, 0);
	par->exit_lambda_min_en = (short)assignation_number(flog, 73, 0, keyword, num_param, num_param_components, 0., 0);
	par->dem_rotation = assignation_number(flog, 74, 0, keyword, num_param, num_param_components, 0., 0);
	par->maxiter_canopy = (long)assignation_number(flog, 75, 0, keyword, num_param, num_param_components, 3., 0);
	par->maxiter_Businger = (long)assignation_number(flog, 76, 0, keyword, num_param, num_param_components, 5., 0);
	par->maxiter_Ts = (long)assignation_number(flog, 77, 0, keyword, num_param, num_param_components, 2., 0);
	par->maxiter_Loc = (long)assignation_number(flog, 78, 0, keyword, num_param, num_param_components, 3., 0);
	par->stabcorr_incanopy = (short)assignation_number(flog, 79, 0, keyword, num_param, num_param_components, 1., 0);
	par->iobsint = (short)assignation_number(flog, 80, 0, keyword, num_param, num_param_components, 1., 0);
	par->dn = assignation_number(flog, 81, 0, keyword, num_param, num_param_components, 1., 0);
	par->slopewt = assignation_number(flog, 82, 0, keyword, num_param, num_param_components, 0., 0);
	par->curvewt = assignation_number(flog, 83, 0, keyword, num_param, num_param_components, 0., 0);
	par->slopewtD = assignation_number(flog, 84, 0, keyword, num_param, num_param_components, par->slopewt, 0);
	par->curvewtD = assignation_number(flog, 85, 0, keyword, num_param, num_param_components, par->curvewt, 0);
	par->slopewtI = assignation_number(flog, 86, 0, keyword, num_param, num_param_components, par->slopewt, 0);
	par->curvewtI = assignation_number(flog, 87, 0, keyword, num_param, num_param_components, par->curvewt, 0);
	par->Zboundary = assignation_number(flog, 88, 0, keyword, num_param, num_param_components, 1.E20, 0);
	par->Tboundary = assignation_number(flog, 89, 0, keyword, num_param, num_param_components, 20., 0);
	par->Fboundary = assignation_number(flog, 90, 0, keyword, num_param, num_param_components, 0., 0);
#else
    //former block 3
    par->latitude = getDoubleValueWithDefault(lConfigStore, "Latitude", geotop::input::gDoubleNoValue, false) ;
    par->longitude = getDoubleValueWithDefault(lConfigStore, "Longitude", geotop::input::gDoubleNoValue, false) ;
    par->Vmin = getDoubleValueWithDefault(lConfigStore, "Vmin", geotop::input::gDoubleNoValue, false) ;
    par->RHmin = getDoubleValueWithDefault(lConfigStore, "RHmin", geotop::input::gDoubleNoValue, false) ;
    par->alpha_snow = getDoubleValueWithDefault(lConfigStore, "AlphaSnow", geotop::input::gDoubleNoValue, false) ;
    par->nsurface = (long)getDoubleValueWithDefault(lConfigStore, "HighestNodeCorrespondsToLayer", geotop::input::gDoubleNoValue, false) ;
    par->tol_energy = getDoubleValueWithDefault(lConfigStore, "HeatEqTol", geotop::input::gDoubleNoValue, false) ;
    par->maxiter_energy = (long)getDoubleValueWithDefault(lConfigStore, "HeatEqMaxIter", geotop::input::gDoubleNoValue, false) ;
    par->min_lambda_en = getDoubleValueWithDefault(lConfigStore, "MinLambdaEnergy", geotop::input::gDoubleNoValue, false) ;
    par->max_times_min_lambda_en = (long)getDoubleValueWithDefault(lConfigStore, "MaxTimesMinLambdaEnergy", geotop::input::gDoubleNoValue, false) ;
    par->exit_lambda_min_en = (short)getDoubleValueWithDefault(lConfigStore, "ExitMinLambdaEnergy", geotop::input::gDoubleNoValue, false) ;
    par->dem_rotation = getDoubleValueWithDefault(lConfigStore, "DEMRotationAngle", geotop::input::gDoubleNoValue, false) ;
    par->maxiter_canopy = (long)getDoubleValueWithDefault(lConfigStore, "CanopyMaxIter", geotop::input::gDoubleNoValue, false) ;
    par->maxiter_Businger = (long)getDoubleValueWithDefault(lConfigStore, "BusingerMaxIter", geotop::input::gDoubleNoValue, false) ;
    par->maxiter_Ts = (long)getDoubleValueWithDefault(lConfigStore, "TsMaxIter", geotop::input::gDoubleNoValue, false) ;
    par->maxiter_Loc = (long)getDoubleValueWithDefault(lConfigStore, "LocMaxIter", geotop::input::gDoubleNoValue, false) ;
    par->stabcorr_incanopy = (short)getDoubleValueWithDefault(lConfigStore, "CanopyStabCorrection", geotop::input::gDoubleNoValue, false) ;
    par->iobsint = (short)getDoubleValueWithDefault(lConfigStore, "Iobsint", geotop::input::gDoubleNoValue, false) ;
    par->dn = getDoubleValueWithDefault(lConfigStore, "Dn", geotop::input::gDoubleNoValue, false) ;
    par->slopewt = getDoubleValueWithDefault(lConfigStore, "SlopeWeight", geotop::input::gDoubleNoValue, false) ;
    par->curvewt = getDoubleValueWithDefault(lConfigStore, "CurvatureWeight", geotop::input::gDoubleNoValue, false) ;
    par->slopewtD = getDoubleValueWithDefault(lConfigStore, "SlopeWeightD", geotop::input::gDoubleNoValue, false) ;
    par->curvewtD = getDoubleValueWithDefault(lConfigStore, "CurvatureWeightD", geotop::input::gDoubleNoValue, false) ;
    par->slopewtI = getDoubleValueWithDefault(lConfigStore, "SlopeWeightI", geotop::input::gDoubleNoValue, false) ;
    par->curvewtI = getDoubleValueWithDefault(lConfigStore, "CurvatureWeightI", geotop::input::gDoubleNoValue, false) ;
    par->Zboundary = getDoubleValueWithDefault(lConfigStore, "ZeroTempAmplitDepth", geotop::input::gDoubleNoValue, false) ;
    par->Tboundary = getDoubleValueWithDefault(lConfigStore, "ZeroTempAmplitTemp", geotop::input::gDoubleNoValue, false) ;
	par->Fboundary = getDoubleValueWithDefault(lConfigStore, "BottomBoundaryHeatFlux", geotop::input::gDoubleNoValue, false) ;
#endif

#ifdef STAGED_FOR_REMOVING
	//former block 4
	itools->swe0 = assignation_number(flog, 91, 0, keyword, num_param, num_param_components, 0., 0);
	itools->rhosnow0 = assignation_number(flog, 92, 0, keyword, num_param, num_param_components, 200., 0);
	itools->Tsnow0 = assignation_number(flog, 93, 0, keyword, num_param, num_param_components, -3., 0);
	itools->agesnow0 = assignation_number(flog, 94, 0, keyword, num_param, num_param_components, 0., 0);
	par->T_rain = assignation_number(flog, 95, 0, keyword, num_param, num_param_components, 3., 0);
	par->T_snow = assignation_number(flog, 96, 0, keyword, num_param, num_param_components, -1., 0);
	par->dew = (short)assignation_number(flog, 97, 0, keyword, num_param, num_param_components, 0., 0);
	par->aep = assignation_number(flog, 98, 0, keyword, num_param, num_param_components, 10., 0);
	par->avo = assignation_number(flog, 99, 0, keyword, num_param, num_param_components, 0.9, 0);
	par->airo = assignation_number(flog, 100, 0, keyword, num_param, num_param_components, 0.65, 0); 
	par->Sr = assignation_number(flog, 101, 0, keyword, num_param, num_param_components, 0.02, 0);  
	par->epsilon_snow = assignation_number(flog, 102, 0, keyword, num_param, num_param_components, 0.98, 0); 
	par->z0_snow = 0.001*assignation_number(flog, 103, 0, keyword, num_param, num_param_components, 0.1, 0); 
	par->snowcorrfact = assignation_number(flog, 104, 0, keyword, num_param, num_param_components, 1., 0); 
	par->raincorrfact = assignation_number(flog, 105, 0, keyword, num_param, num_param_components, 1., 0); 
	par->snow_maxpor = assignation_number(flog, 106, 0, keyword, num_param, num_param_components, 0.7, 0);
	par->drysnowdef_rate = assignation_number(flog, 107, 0, keyword, num_param, num_param_components, 1., 0);
	par->snow_density_cutoff = assignation_number(flog, 108, 0, keyword, num_param, num_param_components, 100., 0);
	par->wetsnowdef_rate = assignation_number(flog, 109, 0, keyword, num_param, num_param_components, 1.5, 0);
	par->snow_viscosity = assignation_number(flog, 110, 0, keyword, num_param, num_param_components, 1.E6, 0);
	par->fetch_up = assignation_number(flog, 111, 0, keyword, num_param, num_param_components, 1000., 0);
	par->fetch_down = assignation_number(flog, 112, 0, keyword, num_param, num_param_components, 100., 0);
	par->Wice_PBSM = assignation_number(flog, 113, 0, keyword, num_param, num_param_components, 0., 0);
	par->Dt_PBSM = assignation_number(flog, 114, 0, keyword, num_param, num_param_components, times->Dt_vector[1], 0);
	par->snow_smin = assignation_number(flog, 115, 0, keyword, num_param, num_param_components, 30., 0);
	par->snow_smax = assignation_number(flog, 116, 0, keyword, num_param, num_param_components, 80., 0);
	par->snow_curv = assignation_number(flog, 117, 0, keyword, num_param, num_param_components, -200., 0);
	
	//former blocks 5/6
	par->max_weq_snow = assignation_number(flog, 118, 0, keyword, num_param, num_param_components, 5., 0);
	par->max_snow_layers = (long)assignation_number(flog, 119, 0, keyword, num_param, num_param_components, 10., 0);

	cod = 120;
#else
    //former block 4
    itools->swe0 = getDoubleValueWithDefault(lConfigStore, "InitSWE", geotop::input::gDoubleNoValue, false) ; ;
    itools->rhosnow0 = getDoubleValueWithDefault(lConfigStore, "InitSnowDensity", geotop::input::gDoubleNoValue, false) ; ;
    itools->Tsnow0 = getDoubleValueWithDefault(lConfigStore, "InitSnowTemp", geotop::input::gDoubleNoValue, false) ; ;
    itools->agesnow0 = getDoubleValueWithDefault(lConfigStore, "InitSnowAge", geotop::input::gDoubleNoValue, false) ; ;
    par->T_rain = getDoubleValueWithDefault(lConfigStore, "ThresTempRain", geotop::input::gDoubleNoValue, false) ; ;
    par->T_snow = getDoubleValueWithDefault(lConfigStore, "ThresTempSnow", geotop::input::gDoubleNoValue, false) ; ;
    par->dew = (short)getDoubleValueWithDefault(lConfigStore, "DewTempOrNormTemp", geotop::input::gDoubleNoValue, false) ; ;
    par->aep = getDoubleValueWithDefault(lConfigStore, "AlbExtParSnow", geotop::input::gDoubleNoValue, false) ; ;
    par->avo = getDoubleValueWithDefault(lConfigStore, "FreshSnowReflVis", geotop::input::gDoubleNoValue, false) ; ;
    par->airo = getDoubleValueWithDefault(lConfigStore, "FreshSnowReflNIR", geotop::input::gDoubleNoValue, false) ; ;
    par->Sr = getDoubleValueWithDefault(lConfigStore, "IrriducibleWatSatSnow", geotop::input::gDoubleNoValue, false) ; ;
    par->epsilon_snow = getDoubleValueWithDefault(lConfigStore, "SnowEmissiv", geotop::input::gDoubleNoValue, false) ; ;
    par->z0_snow = 0.001*getDoubleValueWithDefault(lConfigStore, "SnowRoughness", geotop::input::gDoubleNoValue, false) ; ;
    par->snowcorrfact = getDoubleValueWithDefault(lConfigStore, "SnowCorrFactor", geotop::input::gDoubleNoValue, false) ; ;
    par->raincorrfact = getDoubleValueWithDefault(lConfigStore, "RainCorrFactor", geotop::input::gDoubleNoValue, false) ; ;
    par->snow_maxpor = getDoubleValueWithDefault(lConfigStore, "MaxSnowPorosity", geotop::input::gDoubleNoValue, false) ; ;
	par->drysnowdef_rate = getDoubleValueWithDefault(lConfigStore, "DrySnowDefRate", geotop::input::gDoubleNoValue, false) ; ;
    par->snow_density_cutoff = getDoubleValueWithDefault(lConfigStore, "SnowDensityCutoff", geotop::input::gDoubleNoValue, false) ; ;
    par->wetsnowdef_rate = getDoubleValueWithDefault(lConfigStore, "WetSnowDefRate", geotop::input::gDoubleNoValue, false) ; ;
    par->snow_viscosity = getDoubleValueWithDefault(lConfigStore, "SnowViscosity", geotop::input::gDoubleNoValue, false) ; ;
    par->fetch_up = getDoubleValueWithDefault(lConfigStore, "FetchUp", geotop::input::gDoubleNoValue, false) ; ;
    par->fetch_down = getDoubleValueWithDefault(lConfigStore, "FetchDown", geotop::input::gDoubleNoValue, false) ; ;
    par->Wice_PBSM = getDoubleValueWithDefault(lConfigStore, "BlowingSnowSoftLayerIceContent", geotop::input::gDoubleNoValue, false) ; ;
    par->Dt_PBSM = getDoubleValueWithDefault(lConfigStore, "TimeStepBlowingSnow", geotop::input::gDoubleNoValue, false) ; ;
    par->snow_smin = getDoubleValueWithDefault(lConfigStore, "SnowSMIN", geotop::input::gDoubleNoValue, false) ; ;
    par->snow_smax = getDoubleValueWithDefault(lConfigStore, "SnowSMAX", geotop::input::gDoubleNoValue, false) ; ;
	par->snow_curv = getDoubleValueWithDefault(lConfigStore, "SnowCURV", geotop::input::gDoubleNoValue, false) ; ;
	
	//former blocks 5/6
	par->max_weq_snow = getDoubleValueWithDefault(lConfigStore, "MaxWaterEqSnowLayerContent", geotop::input::gDoubleNoValue, false) ; ;
	par->max_snow_layers = (long)getDoubleValueWithDefault(lConfigStore, "MaxSnowLayersMiddle", geotop::input::gDoubleNoValue, false) ; ;
#endif

#ifdef STAGED_FOR_REMOVING
	n = par->max_snow_layers ;
	if(n < 1){
		fprintf(flog,"Error:: %s must be 1 or larger\n",keyword[119].c_str());
		printf("Error:: %s must be 1 or larger\n",keyword[119].c_str());
		t_error("Fatal Error! Geotop is closed.");
	}
#else
	n = par->max_snow_layers ;
	if(n < 1){
		fprintf(flog,"Error:: MaxSnowLayersMiddle must be 1 or larger\n");
		printf("Error:: MaxSnowLayersMiddle must be 1 or larger\n");
		t_error("Fatal Error! Geotop is closed.");
	}
#endif

#ifdef STAGED_FOR_REMOVING
	par->SWE_bottom = assignation_number(flog, 120, 0, keyword, num_param, num_param_components, 20., 0);
	par->SWE_top = assignation_number(flog, 121, 0, keyword, num_param, num_param_components, 20., 0);
#else
	par->SWE_bottom = getDoubleValueWithDefault(lConfigStore, "SWEbottom", geotop::input::gDoubleNoValue, false) ;
	par->SWE_top = getDoubleValueWithDefault(lConfigStore, "SWEtop", geotop::input::gDoubleNoValue, false) ;
#endif
    
	par->max_snow_layers = (long)floor(par->SWE_bottom/par->max_weq_snow) + (long)floor(par->SWE_top/par->max_weq_snow) + n;
	par->inf_snow_layers.resize(n + 1, 0);
	fprintf(flog,"Max snow layer number: %ld, of which %.0f at the bottom, %ld in the middle, and %.0f at the top.\n",par->max_snow_layers,floor(par->SWE_bottom/par->max_weq_snow),n,floor(par->SWE_top/par->max_weq_snow));
	fprintf(flog,"Infinite Snow layer numbers are numbers: ");
	for (size_t i=1; i<=n; i++) {
		par->inf_snow_layers[i] = (long)floor(par->SWE_bottom/par->max_weq_snow) + i;
		fprintf(flog, "%ld ",par->inf_snow_layers[i]);
	}
	fprintf(flog,"\n");
	
#ifdef STAGED_FOR_REMOVING
	//former block 7
	itools->Dglac0 = assignation_number(flog, 122, 0, keyword, num_param, num_param_components, 0., 0);
	itools->rhoglac0 = assignation_number(flog, 123, 0, keyword, num_param, num_param_components, 800., 0);
	itools->Tglac0 = assignation_number(flog, 124, 0, keyword, num_param, num_param_components, -3., 0);
	par->Sr_glac = assignation_number(flog, 125, 0, keyword, num_param, num_param_components, 0.02, 0);

	//former block 8
	par->max_weq_glac = assignation_number(flog, 126, 0, keyword, num_param, num_param_components, 5., 0);
	n = (long)assignation_number(flog, 127, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->GWE_bottom = assignation_number(flog, 128, 0, keyword, num_param, num_param_components, 0., 0);
	par->GWE_top = assignation_number(flog, 129, 0, keyword, num_param, num_param_components, 0., 0);
#else
	//former block 7
    itools->Dglac0 = getDoubleValueWithDefault(lConfigStore, "InitGlacierDepth", geotop::input::gDoubleNoValue, false) ;
    itools->rhoglac0 = getDoubleValueWithDefault(lConfigStore, "InitGlacierDensity", geotop::input::gDoubleNoValue, false) ;
    itools->Tglac0 = getDoubleValueWithDefault(lConfigStore, "InitGlacierTemp", geotop::input::gDoubleNoValue, false) ;
    par->Sr_glac = getDoubleValueWithDefault(lConfigStore, "IrriducibleWatSatGlacier", geotop::input::gDoubleNoValue, false) ;

    //former block 8
    par->max_weq_glac = getDoubleValueWithDefault(lConfigStore, "MaxWaterEqGlacLayerContent", geotop::input::gDoubleNoValue, false) ;
    n = (long)getDoubleValueWithDefault(lConfigStore, "MaxGlacLayersMiddle", geotop::input::gDoubleNoValue, false) ;

    par->GWE_bottom = getDoubleValueWithDefault(lConfigStore, "GWEbottom", geotop::input::gDoubleNoValue, false) ;
	par->GWE_top = getDoubleValueWithDefault(lConfigStore, "GWEtop", geotop::input::gDoubleNoValue, false) ;
#endif
    
#ifdef STAGED_FOR_REMOVING
	if(n < 1 && (par->GWE_bottom > 0 || par->GWE_top > 0)){
		fprintf(flog,"Error:: %s must be 1 or larger\n",keyword[127].c_str());
		printf("Error:: %s must be 1 or larger\n",keyword[127].c_str());
		t_error("Fatal Error! Geotop is closed.");
	}
#else
	if(n < 1 && (par->GWE_bottom > 0 || par->GWE_top > 0)){
		fprintf(flog,"Error:: MaxGlacLayersMiddle must be 1 or larger\n");
		printf("Error:: MaxGlacLayersMiddle must be 1 or larger\n");
		t_error("Fatal Error! Geotop is closed.");	
	}
#endif

	par->max_glac_layers = (long)floor(par->GWE_bottom/par->max_weq_glac) + (long)floor(par->GWE_top/par->max_weq_glac) + n;
	par->inf_glac_layers.resize(n + 1, 0);
	fprintf(flog,"Max glac layer number: %ld, of which %.0f at the bottom, %ld in the middle, and %.0f at the top.\n",par->max_glac_layers,floor(par->GWE_bottom/par->max_weq_glac),n,floor(par->GWE_top/par->max_weq_glac));
	fprintf(flog,"Infinite Glac layer numbers are numbers: ");
	for (size_t i=1; i<=n; i++) {
		par->inf_glac_layers[i] = (long)floor(par->GWE_bottom/par->max_weq_glac) + i;
		fprintf(flog, "%ld ",par->inf_glac_layers[i]);
	}
	fprintf(flog,"\n");

    par->state_turb = 1;
    
#ifdef STAGED_FOR_REMOVING
	//former block 9
	par->state_lwrad = (short)assignation_number(flog, 130, 0, keyword, num_param, num_param_components, 9., 0);
	par->monin_obukhov = (short)assignation_number(flog, 131, 0, keyword, num_param, num_param_components, 1., 0);
	par->surroundings = (short)assignation_number(flog, 132, 0, keyword, num_param, num_param_components, 0., 0);
		
	//distributed option file
	par->wat_balance = (short)assignation_number(flog, 133, 0, keyword, num_param, num_param_components, 0., 0);
	par->en_balance = (short)assignation_number(flog, 134, 0, keyword, num_param, num_param_components, 0., 0);
	par->blowing_snow = (short)assignation_number(flog, 135, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->Wmin_BS = assignation_number(flog, 136, 0, keyword, num_param, num_param_components, 8., 0);
#else
	//former block 9
    par->state_lwrad = (short)getDoubleValueWithDefault(lConfigStore, "LWinParameterization", geotop::input::gDoubleNoValue, false) ;
    par->monin_obukhov = (short)getDoubleValueWithDefault(lConfigStore, "MoninObukhov", geotop::input::gDoubleNoValue, false) ;
    par->surroundings = (short)getDoubleValueWithDefault(lConfigStore, "Surroundings", geotop::input::gDoubleNoValue, false) ;

    //distributed option file
    par->wat_balance = (short)getDoubleValueWithDefault(lConfigStore, "WaterBalance", geotop::input::gDoubleNoValue, false) ;
    par->en_balance = (short)getDoubleValueWithDefault(lConfigStore, "EnergyBalance", geotop::input::gDoubleNoValue, false) ;
    par->blowing_snow = (short)getDoubleValueWithDefault(lConfigStore, "BlowingSnow", geotop::input::gDoubleNoValue, false) ;
	par->Wmin_BS = getDoubleValueWithDefault(lConfigStore, "MinIceContentForBlowingSnow", geotop::input::gDoubleNoValue, false) ;
#endif
    
#ifdef STAGED_FOR_REMOVING
	cod = 137;
	npoints = 0;
	for (size_t j=1; j<=16; j++) {
		if (npoints < num_param_components[cod + j-1]) {
            npoints = num_param_components[cod + j-1];
        }
	}
#endif

#ifdef STAGED_FOR_REMOVING
	if (par->point_sim == 1) {
		par->chkpt.resize(npoints + 1, ptTOT + 1, 0);
	}else {
		par->chkpt.resize(npoints + 1, 3 + 1, 0);
	}
#else
    std::vector<std::string> lKeywordString ;
    lKeywordString += "PointID", "CoordinatePointX", "CoordinatePointY" ;
	if (par->point_sim == 1) {
        lKeywordString += "PointElevation", "PointLandCoverType", "PointSoilType",
        "PointSlope", "PointAspect", "PointSkyViewFactor", "PointCurvatureNorthSouthDirection", "PointCurvatureWestEastDirection",
        "PointCurvatureNorthwestSoutheastDirection", "PointCurvatureNortheastSouthwestDirection",
        "PointDepthFreeSurface", "PointHorizon", "PointMaxSWE", "PointLatitude", "PointLongitude" ;
	}

    npoints = 0 ;
    for (size_t j=1; j<par->chkpt.getCols(); j++) {
        lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, lKeywordString[j-1], geotop::input::gDoubleNoValue, false, 0, true) ;
		if (npoints < lDoubleTempVector.size()) {
            npoints = lDoubleTempVector.size();
        }
        
	}
    par->chkpt.resize(npoints + 1, lKeywordString.size() + 1, 0);
#endif
    
#ifdef STAGED_FOR_REMOVING
	for (size_t i=1; i<par->chkpt.getRows(); i++) {
		for (size_t j=1; j<par->chkpt.getCols(); j++) {
			par->chkpt[i][j] = assignation_number(flog, cod + j-1, i-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
		}
	}
#else
    for (size_t j=1; j<par->chkpt.getCols(); j++) {
        lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, lKeywordString[j-1], geotop::input::gDoubleNoValue, false, 0, true) ;
        for (size_t i=1; i<par->chkpt.getRows(); i++) {
			par->chkpt[i][j] = lDoubleTempVector[i-1] ;
		}
	}
#endif

#ifdef STAGED_FOR_REMOVING
	cod = 156;
	par->saving_points.resize(num_param_components[cod] + 1,0);
	for (size_t i=1; i<par->saving_points.size(); i++) {
		par->saving_points[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SavingPoints", geotop::input::gDoubleNoValue, false, 0, false) ;
	par->saving_points.resize(lDoubleTempVector.size() + 1, 0);
	for (size_t i=1; i<par->saving_points.size(); i++) {
		par->saving_points[i] = lDoubleTempVector[i-1];
	}
#endif
			
	par->output_soil.resize(par->init_date.size() + 1, 0);
	par->output_snow.resize(par->init_date.size() + 1, 0);
	par->output_glac.resize(par->init_date.size() + 1, 0);
	par->output_surfenergy.resize(par->init_date.size() + 1, 0);
	par->output_vegetation.resize(par->init_date.size() + 1, 0);
	par->output_meteo.resize(par->init_date.size() + 1, 0);

#ifdef STAGED_FOR_REMOVING
	cod = 157;
	par->output_soil[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->output_soil[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_soil[i-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "OutputSoilMaps", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->output_soil.size(); i++) {
		par->output_soil[i] = lDoubleTempVector[i-1] ;
	}
#endif

#ifdef STAGED_FOR_REMOVING
	cod = 158;
	par->output_snow[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->output_snow[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_snow[i-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "OutputSnowMaps", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->output_snow.size(); i++) {
		par->output_snow[i] = lDoubleTempVector[i-1] ;
	}
#endif
    
#ifdef STAGED_FOR_REMOVING
	cod = 159;
	par->output_glac[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->output_glac[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_glac[i-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "OutputGlacierMaps", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->output_glac.size(); i++) {
		par->output_glac[i] = lDoubleTempVector[i-1] ;
	}
#endif
    

#ifdef STAGED_FOR_REMOVING
	cod = 160;
	par->output_surfenergy[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->output_surfenergy[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_surfenergy[i-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "OutputSurfEBALMaps", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->output_surfenergy.size(); i++) {
		par->output_surfenergy[i] = lDoubleTempVector[i-1] ;
	}
#endif
	
#ifdef STAGED_FOR_REMOVING
	cod = 161;
	par->output_vegetation[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->output_vegetation[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_vegetation[i-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "OutputVegetationMaps", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->output_vegetation.size(); i++) {
		par->output_vegetation[i] = lDoubleTempVector[i-1] ;
	}
#endif

#ifdef STAGED_FOR_REMOVING
	cod = 162;
	par->output_meteo[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<par->init_date.size(); i++) {
		par->output_meteo[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_meteo[i-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "OutputMeteoMaps", 0., true, par->init_date.size(), false) ;
	for (size_t i=1; i<par->output_meteo.size(); i++) {
		par->output_meteo[i] = lDoubleTempVector[i-1] ;
	}
#endif
    
	par->output_soil_bin = 0;
	par->output_snow_bin = 0;
	par->output_glac_bin = 0;
	par->output_surfenergy_bin = 0;
	par->output_meteo_bin = 0;

	for (size_t i=1; i<par->init_date.size(); i++) {
		if (par->output_soil[i] > 0) par->output_soil_bin = 1;		
		if (par->output_snow[i] > 0) par->output_snow_bin = 1;
		if (par->output_glac[i] > 0) par->output_glac_bin = 1;
		if (par->output_surfenergy[i] > 0) par->output_surfenergy_bin = 1;
		if (par->output_meteo[i] > 0) par->output_meteo_bin = 1;
	}

#ifdef STAGED_FOR_REMOVING
	cod = 163;
	codn = 164;
	
	if(num_param_components[cod] != num_param_components[codn]){
		fprintf(flog, "Error:: Number of components of parameters %s and %s must be equal\n",keyword[cod].c_str(),keyword[codn].c_str());
		printf("Error:: Number of components of parameters %s and %s must be equal\n",keyword[cod].c_str(),keyword[codn].c_str());
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}

	times->JD_plots.resize(num_param_components[cod] + num_param_components[codn] + 1, 0) ;
	for (size_t i=1; i<(long)(times->JD_plots.size()/2.); i++) {
		times->JD_plots[2*i-1] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 0);
		times->JD_plots[2*i  ] = assignation_number(flog, codn, i-1, keyword, num_param, num_param_components, 0., 0);				
	}
	if (times->JD_plots.size() == 3 && times->JD_plots[1] < 1.E-5 && times->JD_plots[2] < 1.E-5) {
//		free_doublevector(times->JD_plots);
		times->JD_plots.resize(1 + 1, 0) ;
	}
	if (times->JD_plots.size() > 2) {
		for (size_t i=1; i<times->JD_plots.size(); i++) {

			times->JD_plots[i] = convert_dateeur12_JDfrom0(times->JD_plots[i]);
		}
	}
#else
    std::vector<double> lSpecialPlotBegin = getDoubleVectorValueWithDefault(lConfigStore, "SpecialPlotBegin", geotop::input::gDoubleNoValue, false, 0, false) ;
    std::vector<double> lSpecialPlotEnd = getDoubleVectorValueWithDefault(lConfigStore, "SpecialPlotEnd", geotop::input::gDoubleNoValue, false, 0, false) ;

    if(lSpecialPlotBegin.size() != lSpecialPlotEnd.size()){
		fprintf(flog, "Error:: Number of components of parameters SpecialPlotBegin and SpecialPlotEnd must be equal\n");
		printf("Error:: Number of components of parameters SpecialPlotBegin and SpecialPlotEnd must be equal\n");
		t_error("Fatal Error! Geotop is closed. See failing report.");
	}

    times->JD_plots.resize(lSpecialPlotBegin.size() + lSpecialPlotEnd.size() + 1, 0) ;

    for (size_t i=1; i<(size_t)(times->JD_plots.size()/2.); i++) {
		times->JD_plots[2*i-1] = lSpecialPlotBegin[i-1];
		times->JD_plots[2*i  ] = lSpecialPlotEnd[i-1];
	}
	if (times->JD_plots.size() == 3 && times->JD_plots[1] < 1.E-5 && times->JD_plots[2] < 1.E-5) {
        //		free_doublevector(times->JD_plots);
		times->JD_plots.resize(1 + 1, 0) ;
	}
	if (times->JD_plots.size() > 2) {
		for (size_t i=1; i<times->JD_plots.size(); i++) {
            
			times->JD_plots[i] = convert_dateeur12_JDfrom0(times->JD_plots[i]);
		}
	}
#endif

	//initial condition on the water pressure
#ifdef STAGED_FOR_REMOVING
	par->nsoiltypes = (long)assignation_number(flog, 165, 0, keyword, num_param, num_param_components, 1., 0);
#else
	par->nsoiltypes = (long)getDoubleValueWithDefault(lConfigStore, "SoilLayerTypes", geotop::input::gDoubleNoValue, false) ; ;
#endif
    if (par->nsoiltypes < 1)
        par->nsoiltypes = 1;

	itools->init_water_table_depth.resize(par->nsoiltypes + 1, 0);

#ifdef STAGED_FOR_REMOVING
	cod = 166;
	itools->init_water_table_depth[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 5000., 0);
	for (k=2; k<=par->nsoiltypes; k++) {
		itools->init_water_table_depth[k] = assignation_number(flog, cod, k-1, keyword, num_param, num_param_components, itools->init_water_table_depth[k-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "InitWaterTableDepth", 0., true, par->nsoiltypes, false) ;
	for (size_t i=1; i< itools->init_water_table_depth.size(); i++) {
		itools->init_water_table_depth[i] = lDoubleTempVector[i-1] ;
	}
#endif

    //soil properties and discretization
#ifdef STAGED_FOR_REMOVING
	par->soil_type_land_default = (long)assignation_number(flog, 167, 0, keyword, num_param, num_param_components, 1., 0);
#else
	par->soil_type_land_default = (long)getDoubleValueWithDefault(lConfigStore, "DefaultSoilTypeLand", geotop::input::gDoubleNoValue, false) ;
#endif

	if(par->soil_type_land_default<1 || par->soil_type_land_default>par->nsoiltypes){
		fprintf(flog, "Error:: Soil_type_land_default lower than 0 or higher than soil types numbers");
		printf("Error:: Soil_type_land_default lower than 0 or higher than soil types numbers");
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}

#ifdef STAGED_FOR_REMOVING
	par->soil_type_chan_default = (long)assignation_number(flog, 168, 0, keyword, num_param, num_param_components, 1., 0);
#else
	par->soil_type_chan_default = (long)getDoubleValueWithDefault(lConfigStore, "DefaultSoilTypeChannel", geotop::input::gDoubleNoValue, false) ;
#endif

	if(par->soil_type_chan_default<1 || par->soil_type_chan_default>par->nsoiltypes){
		fprintf(flog, "Error:  Soil_type_chan_default lower than 0 or higher than soil types numbers");
		printf("Error:  Soil_type_chan_default lower than 0 or higher than soil types numbers");
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}

#ifdef STAGED_FOR_REMOVING
	cod = 169;
	a = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
    
	//there is a specific layer discretization
	if ( (long)a != geotop::input::gDoubleNoValue && num_param_components[cod] > 1) {
		
		nsoillayers = num_param_components[cod];
		sl->pa.resize(par->nsoiltypes, nsoilprop, nsoillayers) ;
		
		sl->pa[1][jdz][1] = a;
		for (size_t i=2; i<=sl->pa.getCh(); i++) {
			sl->pa[1][jdz][i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, sl->pa[1][jdz][i-1], 0);
		}
				
	}else {
		
		if((long)a==geotop::input::gDoubleNoValue) a=100.;
		nsoillayers = (long)assignation_number(flog, 170, 0, keyword, num_param, num_param_components, 5., 0);
		
        sl->pa.resize(par->nsoiltypes + 1, nsoilprop + 1, nsoillayers + 1);
        
		for (size_t i=1; i<sl->pa.getCh(); i++) {
			sl->pa[1][jdz][i] = a;
		}

	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilLayerThicknesses", geotop::input::gDoubleNoValue, true, 0, true) ;
    a = lDoubleTempVector[0] ;
    
	//there is a specific layer discretization
	if ( (long)a != geotop::input::gDoubleNoValue && lDoubleTempVector.size() > 1) {
		
		nsoillayers = lDoubleTempVector.size();
		sl->pa.resize(par->nsoiltypes + 1, nsoilprop + 1, nsoillayers + 1) ;
		
        lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilLayerThicknesses", a, true, nsoillayers, true) ;
		for (size_t i=1; i<sl->pa.getCh(); i++) {
			sl->pa[1][jdz][i] = lDoubleTempVector[i-1] ;
		}
        
	}else {
		
		if((long)a==geotop::input::gDoubleNoValue) a=100.;
		nsoillayers = (long)getDoubleValueWithDefault(lConfigStore, "SoilLayerNumber", geotop::input::gDoubleNoValue, false) ;
		
        sl->pa.resize(par->nsoiltypes + 1, nsoilprop + 1, nsoillayers + 1);
        
		for (size_t i=1; i<sl->pa.getCh(); i++) {
			sl->pa[1][jdz][i] = a;
		}
        
	}
#endif
	
	//first layer
	size_t lStartIndex = 1;
#ifdef STAGED_FOR_REMOVING
	sl->pa[1][jpsi][lStartIndex] = assignation_number(flog, 171, lStartIndex-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	sl->pa[1][jT][lStartIndex] = assignation_number(flog, 172, lStartIndex-1, keyword, num_param, num_param_components, 5., 0);
	sl->pa[1][jKn][lStartIndex] = assignation_number(flog, 173, lStartIndex-1, keyword, num_param, num_param_components, 1.E-4, 0);
	sl->pa[1][jKl][lStartIndex] = assignation_number(flog, 174, lStartIndex-1, keyword, num_param, num_param_components, 1.E-4, 0);
	sl->pa[1][jres][lStartIndex] = assignation_number(flog, 175, lStartIndex-1, keyword, num_param, num_param_components, 0.05, 0);
	sl->pa[1][jwp][lStartIndex] = assignation_number(flog, 176, lStartIndex-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	sl->pa[1][jfc][lStartIndex] = assignation_number(flog, 177, lStartIndex-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	sl->pa[1][jsat][lStartIndex] = assignation_number(flog, 178, lStartIndex-1, keyword, num_param, num_param_components, 0.5, 0);
	sl->pa[1][ja][lStartIndex] = assignation_number(flog, 179, lStartIndex-1, keyword, num_param, num_param_components, 0.004, 0);
	sl->pa[1][jns][lStartIndex] = assignation_number(flog, 180, lStartIndex-1, keyword, num_param, num_param_components, 1.3, 0);
	sl->pa[1][jv][lStartIndex] = assignation_number(flog, 181, lStartIndex-1, keyword, num_param, num_param_components, 0.5, 0);
	sl->pa[1][jkt][lStartIndex] = assignation_number(flog, 182, lStartIndex-1, keyword, num_param, num_param_components, 2.5, 0);
	sl->pa[1][jct][lStartIndex] = assignation_number(flog, 183, lStartIndex-1, keyword, num_param, num_param_components, 1.E6, 0);
	sl->pa[1][jss][lStartIndex] = assignation_number(flog, 184, lStartIndex-1, keyword, num_param, num_param_components, 1.E-7, 0);
#else
    sl->pa[1][jpsi][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "InitSoilPressure", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jT][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "InitSoilTemp", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jKn][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "NormalHydrConductivity", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jKl][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "LateralHydrConductivity", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jres][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "ThetaRes", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jwp][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "WiltingPoint", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jfc][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "FieldCapacity", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jsat][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "ThetaSat", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][ja][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "AlphaVanGenuchten", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jns][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "NVanGenuchten", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jv][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "VMualem", geotop::input::gDoubleNoValue, false, 0, false)[lStartIndex-1] ;
	sl->pa[1][jkt][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "ThermalConductivitySoilSolids", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jct][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "ThermalCapacitySoilSolids", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
	sl->pa[1][jss][lStartIndex] = getDoubleVectorValueWithDefault(lConfigStore, "SpecificStorativity", geotop::input::gDoubleNoValue, false, 0, true)[lStartIndex-1] ;
#endif

#ifdef STAGED_FOR_REMOVING
	//other layers
	for (size_t i=2; i<sl->pa.getCh(); i++) {
		for (size_t j=2; j<sl->pa.getRh(); j++) {
            sl->pa(1,j,i) = assignation_number(flog, 171 + j-2, i-1, keyword, num_param, num_param_components, sl->pa(1,j,i-1), 0);
		}
	}
#else
    //other layers
    lKeywordString.clear() ;
    lKeywordString += "InitSoilPressure", "InitSoilTemp", "NormalHydrConductivity", "LateralHydrConductivity",
        "ThetaRes", "WiltingPoint", "FieldCapacity", "ThetaSat", "AlphaVanGenuchten", "NVanGenuchten",
        "VMualem", "ThermalConductivitySoilSolids", "ThermalCapacitySoilSolids", "SpecificStorativity" ;

    for (size_t j=2; j<sl->pa.getRh(); j++) {
        lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, lKeywordString[j-2], sl->pa(1,j,1), true, sl->pa.getCh(), true) ;
        for (size_t i=2; i<sl->pa.getCh(); i++) {
            sl->pa(1,j,i) = lDoubleTempVector[i-1];
		}
	}
#endif

	//field capacity (-0.333 bar) and wilting point (-15 bar)
	for (size_t i=1; i<sl->pa.getCh(); i++){
		if( (long)sl->pa(1,jfc,i) == geotop::input::gDoubleNoValue){
			sl->pa[1][jfc][i] = teta_psi( (-1./3.)*1.E5/GTConst::GRAVITY, 0., sl->pa(1,jsat,i), sl->pa(1,jres,i), sl->pa(1,ja,i),
										 sl->pa(1,jns,i), 1.-1./sl->pa(1,jns,i), GTConst::PsiMin, sl->pa(1,jss,i));
		}
		
		if( (long)sl->pa(1,jwp,i) == geotop::input::gDoubleNoValue){
			sl->pa(1,jwp,i) = teta_psi( -15.*1.E5/GTConst::GRAVITY, 0., sl->pa(1,jsat,i), sl->pa(1,jres,i), sl->pa(1,ja,i),
                                       sl->pa(1,jns,i), 1.-1./sl->pa(1,jns,i), GTConst::PsiMin, sl->pa(1,jss,i));
		}
	}
	
	//other soil types
    for (k=2; k<=par->nsoiltypes; k++) {
        for (size_t i=1; i<sl->pa.getCh(); i++) {
            for (size_t j=1; j<sl->pa.getRh(); j++) {
                sl->pa(k,j,i) = sl->pa(1,j,i);
            }
        }
    }
	
	//use water table for water pressure
    for (k=1; k<=par->nsoiltypes; k++) {
        occurring = 0;//check if psi initial has at least one novalue
        /*for (i=1; i<=sl->pa->nch; i++) {
         if ( (long)sl->pa->co[k][jpsi][i] == geotop::input::gDoubleNoValue) occurring = 1;*/
        for (size_t i=1; i<sl->pa.getCh(); i++) {
            if ( (long)sl->pa(k,jpsi,i) == geotop::input::gDoubleNoValue) occurring = 1;
        }
        //      if (occurring == 0) sl->init_water_table_depth->co[k] = geotop::input::gDoubleNoValue;
        if (occurring == 0) sl->init_water_table_depth[k] = geotop::input::gDoubleNoValue;
    }

#ifdef STAGED_FOR_REMOVING
	cod = 185;
	itools->pa_bed.resize(1 + 1, nsoilprop + 1, nsoillayers + 1);
	for (size_t i=1; i<=nsoillayers; i++) {
		itools->pa_bed(1,jdz,i) = sl->pa(1,jdz,i);
	}
	for (size_t j=1; j<=nsoilprop; j++) {
		if(j != jdz){
			itools->pa_bed(1,j,1) = assignation_number(flog, cod+j-2, 0, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
		}
	}
	for (size_t i=2; i<=nsoillayers; i++) {
		for (size_t j=1; j<=nsoilprop; j++) {
			if(j != jdz)
                itools->pa_bed(1,j,i) = assignation_number(flog, cod+j-2, i-1, keyword, num_param, num_param_components, itools->pa_bed(1,j,i-1), 0);
		}
	}
#else
	itools->pa_bed.resize(1 + 1, nsoilprop + 1, nsoillayers + 1);
	for (size_t i=1; i<=nsoillayers; i++) {
		itools->pa_bed(1,jdz,i) = sl->pa(1,jdz,i);
	}
    
    //other layers
    lKeywordString.clear() ;
    lKeywordString += "InitSoilPressureBedrock", "InitSoilTempBedrock", "NormalHydrConductivityBedrock",
        "LateralHydrConductivityBedrock", "ThetaResBedrock", "WiltingPointBedrock",
        "FieldCapacityBedrock", "ThetaSatBedrock", "AlphaVanGenuchtenBedrock", "NVanGenuchtenBedrock",
        "VMualemBedrock", "ThermalConductivitySoilSolidsBedrock", "ThermalCapacitySoilSolidsBedrock",
        "SpecificStorativityBedrock" ;

    for (size_t j=1; j<itools->pa_bed.getCh(); j++) {
        lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, lKeywordString[j-1], geotop::input::gDoubleNoValue, true, nsoillayers, true) ;
        for (size_t i=1; i<itools->pa_bed.getCh(); i++) {
			if(j != jdz)
                itools->pa_bed(1,j,i) = lDoubleTempVector[i-1] ;
		}
	}
#endif
    
	//field capacity (-0.333 bar) and wilting point (-15 bar)
	for (size_t i=1; i<sl->pa.getCh(); i++){
		if( (long)itools->pa_bed(1,jsat,i) != geotop::input::gDoubleNoValue && (long)itools->pa_bed(1,jres,i) != geotop::input::gDoubleNoValue &&
			(long)itools->pa_bed(1,ja,i) != geotop::input::gDoubleNoValue && (long)itools->pa_bed(1,jns,i) != geotop::input::gDoubleNoValue &&
			(long)itools->pa_bed(1,jss,i) ) {
			if ( (long)itools->pa_bed(1,jfc,i) == geotop::input::gDoubleNoValue) {
				itools->pa_bed(1,jfc,i) = teta_psi( (-1./3.)*1.E5/GTConst::GRAVITY, 0., itools->pa_bed[1][jsat][i], itools->pa_bed(1,jres,i), itools->pa_bed(1,ja,i),
										 itools->pa_bed(1,jns,i), 1.-1./itools->pa_bed(1,jns,i), GTConst::PsiMin, itools->pa_bed(1,jss,i));
			}
			if ( (long)itools->pa_bed(1,jwp,i) == geotop::input::gDoubleNoValue) {
				itools->pa_bed(1,jwp,i) = teta_psi( -15.*1.E5/GTConst::GRAVITY, 0., itools->pa_bed(1,jsat,i), itools->pa_bed(1,jres,i), itools->pa_bed(1,ja,i),
										 itools->pa_bed(1,jns,i), 1.-1./itools->pa_bed(1,jns,i), GTConst::PsiMin, itools->pa_bed(1,jss,i));
			}
		}	
	}

	//meteo stations
#ifdef STAGED_FOR_REMOVING
	cod = 199;
	met->imeteo_stations.resize(num_param_components[cod] + 1);
	met->imeteo_stations[1] = (long)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	if ( met->imeteo_stations[1] != geotop::input::gDoubleNoValue ) {
		for (size_t i=2; i<=num_param_components[cod]; i++) {
			met->imeteo_stations[i] = (long)assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 1);
		}
		nmeteo_stations = num_param_components[cod];
	}else {
		nmeteo_stations = (long)assignation_number(flog, 200, 0, keyword, num_param, num_param_components, 1., 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationsID", geotop::input::gDoubleNoValue, false, 0, true) ;
	met->imeteo_stations.resize(lDoubleTempVector.size() + 1);
    met->imeteo_stations[1] = lDoubleTempVector[0] ;
	if ( met->imeteo_stations[1] != geotop::input::gDoubleNoValue ) {
		for (size_t i=2; i < met->imeteo_stations.size() ; i++) {
            double lValue = lDoubleTempVector[i-1];
            if (lValue == geotop::input::gDoubleNoValue)
                lValue = 0. ;
			met->imeteo_stations[i] = (long)lValue;
		}
		nmeteo_stations = lDoubleTempVector.size();
	}else {
		nmeteo_stations = (long) getDoubleValueWithDefault(lConfigStore, "NumberOfMeteoStations", geotop::input::gDoubleNoValue, false);
	}
#endif

    //	met->st=(METEO_STATIONS *)malloc(sizeof(METEO_STATIONS));
    met->st = new MeteoStations();
    size_t lMeteoStationContainerSize = nmeteo_stations+1 ;
    if(!met->st) t_error("meteo_stations was not allocated");
    met->st->E.resize(lMeteoStationContainerSize);
    met->st->N.resize(lMeteoStationContainerSize);
    met->st->lat.resize(lMeteoStationContainerSize);
    met->st->lon.resize(lMeteoStationContainerSize);
    met->st->Z.resize(lMeteoStationContainerSize);
    met->st->sky.resize(lMeteoStationContainerSize);
    met->st->ST.resize(lMeteoStationContainerSize);
    met->st->Vheight.resize(lMeteoStationContainerSize);
    met->st->Theight.resize(lMeteoStationContainerSize);
    
#ifdef STAGED_FOR_REMOVING
	lStartIndex = 1;
	met->st->E[lStartIndex] = assignation_number(flog, 201, lStartIndex-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	met->st->N[lStartIndex] = assignation_number(flog, 202, lStartIndex-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	met->st->lat[lStartIndex] = assignation_number(flog, 203, lStartIndex-1, keyword, num_param, num_param_components, par->latitude, 0);
	met->st->lon[lStartIndex] = assignation_number(flog, 204, lStartIndex-1, keyword, num_param, num_param_components, par->longitude, 0);
	met->st->Z[lStartIndex] = assignation_number(flog, 205, lStartIndex-1, keyword, num_param, num_param_components, 0., 0);
	met->st->sky[lStartIndex] = assignation_number(flog, 206, lStartIndex-1, keyword, num_param, num_param_components, 1., 0);
	met->st->ST[lStartIndex] = assignation_number(flog, 207, lStartIndex-1, keyword, num_param, num_param_components, par->ST, 0);
	
	for (size_t i=2; i<=nmeteo_stations; i++) {
		met->st->E[i] = assignation_number(flog, 201, i-1, keyword, num_param, num_param_components, met->st->E[i-1], 0);
		met->st->N[i] = assignation_number(flog, 202, i-1, keyword, num_param, num_param_components, met->st->N[i-1], 0);
		met->st->lat[i] = assignation_number(flog, 203, i-1, keyword, num_param, num_param_components, met->st->lat[i-1], 0);
		met->st->lon[i] = assignation_number(flog, 204, i-1, keyword, num_param, num_param_components, met->st->lon[i-1], 0);
		met->st->Z[i] = assignation_number(flog, 205, i-1, keyword, num_param, num_param_components, met->st->Z[i-1], 0);
		met->st->sky[i] = assignation_number(flog, 206, i-1, keyword, num_param, num_param_components, met->st->sky[i-1], 0);
		met->st->ST[i] = assignation_number(flog, 207, i-1, keyword, num_param, num_param_components, met->st->ST[i-1], 0);
 	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationCoordinateX", geotop::input::gDoubleNoValue, true, nmeteo_stations, false);
    for(size_t i=1 ; i < lMeteoStationContainerSize ; i++) {
        met->st->E[i] = lDoubleTempVector[i-1];
    }
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationCoordinateY", geotop::input::gDoubleNoValue, true, nmeteo_stations, false);
    for(size_t i=1 ; i < lMeteoStationContainerSize ; i++) {
        met->st->N[i] = lDoubleTempVector[i-1];
    }
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationLatitude", par->latitude, true, nmeteo_stations, false);
    for(size_t i=1 ; i < lMeteoStationContainerSize ; i++) {
        met->st->lat[i] = lDoubleTempVector[i-1];
    }
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationLongitude", par->longitude, true, nmeteo_stations, false);
    for(size_t i=1 ; i < lMeteoStationContainerSize ; i++) {
        met->st->lon[i] = lDoubleTempVector[i-1];
    }
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationElevation", 0., true, nmeteo_stations, false);
    for(size_t i=1 ; i < lMeteoStationContainerSize ; i++) {
        met->st->Z[i] = lDoubleTempVector[i-1];
    }
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationSkyViewFactor", 1., true, nmeteo_stations, false);
    for(size_t i=1 ; i < lMeteoStationContainerSize ; i++) {
        met->st->sky[i] = lDoubleTempVector[i-1];
    }
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationStandardTime", par->ST, true, nmeteo_stations, false);
    for(size_t i=1 ; i < lMeteoStationContainerSize ; i++) {
        met->st->ST[i] = lDoubleTempVector[i-1];
    }
#endif

#ifdef STAGED_FOR_REMOVING
	a = assignation_number(flog, 208, 0, keyword, num_param, num_param_components, 10., 0);
	for (size_t i=1; i<=nmeteo_stations; i++) {
		met->st->Vheight[i] = a;
	}
	
	a = assignation_number(flog, 209, 0, keyword, num_param, num_param_components, 2., 0);
	for (size_t i=1; i<=nmeteo_stations; i++) {
		met->st->Theight[i] = a;
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationWindVelocitySensorHeight", geotop::input::gDoubleNoValue, true, nmeteo_stations, false) ;
	for (size_t i=1; i<lMeteoStationContainerSize; i++) {
		met->st->Vheight[i] = lDoubleTempVector[i-1];
	}

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "MeteoStationTemperatureSensorHeight", geotop::input::gDoubleNoValue, true, nmeteo_stations, false) ;
	for (size_t i=1; i<lMeteoStationContainerSize; i++) {
		met->st->Theight[i] = lDoubleTempVector[i-1];
	}
#endif

	//lapse rates (cyclic)
	n = (long)nlstot;
	met->LRcnc = (long*)malloc(n*sizeof(long));
	met->LRcnc[ilsDate12] = 1;
#ifdef STAGED_FOR_REMOVING
	met->LRcnc[ilsTa] = num_param_components[210];
	met->LRcnc[ilsTdew] = num_param_components[211];
	met->LRcnc[ilsPrec] = num_param_components[212];
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "LapseRateTemp", geotop::input::gDoubleNoValue, false, 0, false) ;
	met->LRcnc[ilsTa] = lDoubleTempVector.size() ;
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "LapseRateDewTemp", geotop::input::gDoubleNoValue, false, 0, true) ;
	met->LRcnc[ilsTdew] = lDoubleTempVector.size() ;
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "LapseRatePrec", geotop::input::gDoubleNoValue, false, 0, false) ;
	met->LRcnc[ilsPrec] = lDoubleTempVector.size() ;
#endif
	met->LRc = (double**)malloc(n*sizeof(double*));

#ifdef STAGED_FOR_REMOVING
    size_t lColumnIndexJ ;
	for (size_t i=0; i<nlstot; i++) {
		met->LRc[i] = (double*)malloc(met->LRcnc[i]*sizeof(double));
		for (lColumnIndexJ=0; lColumnIndexJ<met->LRcnc[i]; lColumnIndexJ++) {
			if(i==ilsDate12)
                met->LRc[i][lColumnIndexJ] = 0.;
			if(i==ilsTa) {
                met->LRc[i][lColumnIndexJ] = assignation_number(flog, 210, lColumnIndexJ, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
            }
			if(i==ilsTdew) {
                met->LRc[i][lColumnIndexJ] = assignation_number(flog, 211, lColumnIndexJ, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
            }
			if(i==ilsPrec) {
                met->LRc[i][lColumnIndexJ] = assignation_number(flog, 212, lColumnIndexJ, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
            }
		}
	}
	par->MinIncrFactWithElev = assignation_number(flog, 213, lColumnIndexJ, keyword, num_param, num_param_components, 0.1, 0);
	par->MaxIncrFactWithElev = assignation_number(flog, 214, lColumnIndexJ, keyword, num_param, num_param_components, 4.4, 0);
#else
    size_t lColumnIndexJ ;
    std::vector<double> lLapseRateTemp = getDoubleVectorValueWithDefault(lConfigStore, "LapseRateTemp", geotop::input::gDoubleNoValue, false, 0, true) ;
    std::vector<double> lLapseRateDewTemp = getDoubleVectorValueWithDefault(lConfigStore, "LapseRateDewTemp", geotop::input::gDoubleNoValue, false, 0, true) ;
    std::vector<double> lLapseRatePrec = getDoubleVectorValueWithDefault(lConfigStore, "LapseRatePrec", geotop::input::gDoubleNoValue, false, 0, true) ;
	for (size_t i=0; i<nlstot; i++) {
		met->LRc[i] = (double*)malloc(met->LRcnc[i]*sizeof(double));
		for (lColumnIndexJ=0; lColumnIndexJ<met->LRcnc[i]; lColumnIndexJ++) {
			if(i==ilsDate12)
                met->LRc[i][lColumnIndexJ] = 0.;
			if(i==ilsTa) {
                if (lColumnIndexJ <= lLapseRateTemp.size())
                    met->LRc[i][lColumnIndexJ] = lLapseRateTemp[lColumnIndexJ];
                else
                    met->LRc[i][lColumnIndexJ] = geotop::input::gDoubleNoValue ;
            }
			if(i==ilsTdew) {
                if (lColumnIndexJ <= lLapseRateDewTemp.size())
                    met->LRc[i][lColumnIndexJ] = lLapseRateDewTemp[lColumnIndexJ];
                else
                    met->LRc[i][lColumnIndexJ] = geotop::input::gDoubleNoValue ;
            }
			if(i==ilsPrec) {
                if (lColumnIndexJ <= lLapseRatePrec.size())
                    met->LRc[i][lColumnIndexJ] = lLapseRatePrec[lColumnIndexJ];
                else
                    met->LRc[i][lColumnIndexJ] = geotop::input::gDoubleNoValue ;
            }
		}
	}

	par->MinIncrFactWithElev =  getDoubleValueWithDefault(lConfigStore, "MinPrecIncreaseFactorWithElev", geotop::input::gDoubleNoValue, false) ;
	par->MaxIncrFactWithElev =  getDoubleValueWithDefault(lConfigStore, "MaxPrecDecreaseFactorWithElev", geotop::input::gDoubleNoValue, false) ;
#endif
	
	//output point column
	n = (long)otot;
	opnt = (long*)malloc(n*sizeof(long));
	ipnt = (short*)malloc(n*sizeof(short));

#ifdef STAGED_FOR_REMOVING
	par->all_point = (short)assignation_number(flog, 294, 0, keyword, num_param, num_param_components, 0., 0);
#else
    par->all_point = (short)getDoubleValueWithDefault(lConfigStore, "PointAll", geotop::input::gDoubleNoValue, false) ;
#endif

	if (par->all_point == 1) {
		
		for (size_t i=0; i<n; i++) {
			ipnt[i] = 1;
			opnt[i] = i;
		}
		
		nopnt = n;
		
	}else {
		
		for (size_t i=0; i<n; i++) {
			ipnt[i] = 0;
			opnt[i] = -1;
		}

#ifdef STAGED_FOR_REMOVING
		for (size_t i=0; i<n; i++) {
			lColumnIndexJ = (long)assignation_number(flog, 215+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n){
				opnt[lColumnIndexJ-1] = i;
				ipnt[i] = 1;
			}
		}
#else
        lKeywordString.clear() ;
        lKeywordString += "DatePoint", "JulianDayFromYear0Point", "TimeFromStartPoint",
        "PeriodPoint", "RunPoint", "IDPointPoint", "PsnowPoint",
        "PrainPoint", "PsnowNetPoint", "PrainNetPoint", "PrainOnSnowPoint",
        "WindSpeedPoint", "WindDirPoint", "RHPoint", "AirPressPoint",
        "AirTempPoint", "TDewPoint", "TsurfPoint", "TvegPoint",
        "TCanopyAirPoint", "SurfaceEBPoint", "SoilHeatFluxPoint", "SWinPoint",
        "SWbeamPoint", "SWdiffPoint", "LWinPoint", "LWinMinPoint",
        "LWinMaxPoint", "SWNetPoint", "LWNetPoint", "HPoint",
        "LEPoint", "CanopyFractionPoint", "LSAIPoint", "z0vegPoint",
        "d0vegPoint", "EstoredCanopyPoint", "SWvPoint", "LWvPoint",
        "HvPoint", "LEvPoint", "HgUnvegPoint", "LEgUnvegPoint",
        "HgVegPoint", "LEgVegPoint", "EvapSurfacePoint", "TraspCanopyPoint",
        "WaterOnCanopyPoint", "SnowOnCanopyPoint", "QVegPoint", "QSurfPoint",
        "QAirPoint", "QCanopyAirPoint", "LObukhovPoint", "LObukhovCanopyPoint",
        "WindSpeedTopCanopyPoint", "DecayKCanopyPoint", "SWupPoint", "LWupPoint",
        "HupPoint", "LEupPoint", "SnowDepthPoint", "SWEPoint",
        "SnowDensityPoint", "SnowTempPoint", "SnowMeltedPoint", "SnowSublPoint",
        "SWEBlownPoint", "SWESublBlownPoint", "GlacDepthPoint", "GWEPoint",
        "GlacDensityPoint", "GlacTempPoint", "GlacMeltedPoint", "GlacSublPoint",
        "LowestThawedSoilDepthPoint", "HighestThawedSoilDepthPoint", "LowestWaterTableDepthPoint", "HighestWaterTableDepthPoint" ;

		for (size_t i=0; i<n; i++) {
			lColumnIndexJ = (long)getDoubleValueWithDefault(lConfigStore, lKeywordString[i], geotop::input::gDoubleNoValue, false) ;
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n){
				opnt[lColumnIndexJ-1] = i;
				ipnt[i] = 1;
			}
		}
#endif
		nopnt = 0;
		for (size_t i=0; i<n; i++) {
			if(opnt[i] > 0) nopnt = i+1;
		}
	}
	
	//output basin column
	n = (long)ootot;
	obsn = (long*)malloc(n*sizeof(long));
	ibsn = (short*)malloc(n*sizeof(short));
	
#ifdef STAGED_FOR_REMOVING
	par->all_basin = (short)assignation_number(flog, 322, 0, keyword, num_param, num_param_components, 0., 0);
#else
    par->all_basin = (short)getDoubleValueWithDefault(lConfigStore, "BasinAll", geotop::input::gDoubleNoValue, false) ;
#endif

	if (par->all_basin == 1) {
		
		for (size_t i=0; i<n; i++) {
			ibsn[i] = 1;
			obsn[i] = i;
		}
		
		nobsn = (long)ootot;
		
	}else {
		
		for (size_t i=0; i<n; i++) {
			ibsn[i] = 0;
			obsn[i] = -1;
		}
#ifdef STAGED_FOR_REMOVING
		for (size_t i=0; i<n; i++) {
			lColumnIndexJ = (long)assignation_number(flog, 295+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n){
				obsn[lColumnIndexJ-1] = i;
				ibsn[i] = 1;
			}
		}
#else
        lKeywordString.clear() ;
        lKeywordString += "DateBasin", "JulianDayFromYear0Basin", "TimeFromStartBasin",
            "PeriodBasin", "RunBasin", "PRainNetBasin", "PSnowNetBasin", "PRainBasin",
            "PSnowBasin", "PNetBasin", "AirTempBasin", "TSurfBasin", "TvegBasin",
            "EvapSurfaceBasin", "TraspCanopyBasin", "LEBasin", "HBasin", "SWNetBasin",
            "LWNetBasin", "LEvBasin", "HvBasin", "SWvBasin", "LWvBasin",
            "SWinBasin", "LWinBasin", "MassErrorBasin", "MeanTimeStep" ;

		for (size_t i=0; i<n; i++) {
			lColumnIndexJ = (long)getDoubleValueWithDefault(lConfigStore, lKeywordString[i], geotop::input::gDoubleNoValue, false) ;
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n){
				obsn[lColumnIndexJ-1] = i;
				ibsn[i] = 1;
			}
		}
#endif
		nobsn = 0;
		for (size_t i=0; i<n; i++) {
			if(obsn[i] > 0) nobsn = i+1;
		}
	}
	
#ifdef STAGED_FOR_REMOVING
	cod = 356;
	n = 0;
	do{
		a = assignation_number(flog, cod, n, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
		if ((long)a != geotop::input::gDoubleNoValue) n++;
	}while ((long)a != geotop::input::gDoubleNoValue && n<=1000000);	
	if (n==0) n=1;
	par->soil_plot_depths.resize(n + 1);
	for (n=1; n<par->soil_plot_depths.size(); n++) {
		par->soil_plot_depths[n] = assignation_number(flog, cod, n-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	}	
	
	cod = 357;
	n = 0;
	do{
		a = assignation_number(flog, cod, n, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
		if ((long)a != geotop::input::gDoubleNoValue) n++;
	}while ((long)a != geotop::input::gDoubleNoValue && n<=1000000);	
	if (n==0) n=1;
	par->snow_plot_depths.resize(n + 1);
	for (n=1; n<par->snow_plot_depths.size(); n++) {
		par->snow_plot_depths[n] = assignation_number(flog, cod, n-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	}	
	
	cod = 358;
	n = 0;
	do{
		a = assignation_number(flog, cod, n, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
		if ((long)a != geotop::input::gDoubleNoValue) n++;
	}while ((long)a != geotop::input::gDoubleNoValue && n<=1000000);	
	if (n==0) n=1;
	par->glac_plot_depths.resize(n + 1);
	for (n=1; n<par->glac_plot_depths.size(); n++) {
		par->glac_plot_depths[n] = assignation_number(flog, cod, n-1, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	}	
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SoilPlotDepths", geotop::input::gDoubleNoValue, false, 0, true) ;
    size_t lSoilPlotDepthsSize = lDoubleTempVector.size() + 1;
    par->soil_plot_depths.resize(lSoilPlotDepthsSize);
	for (size_t i=1; i<lSoilPlotDepthsSize; i++) {
		par->soil_plot_depths[i] = lDoubleTempVector[i-1] ;
	}

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SnowPlotDepths", geotop::input::gDoubleNoValue, false, 0, true) ;
    lSoilPlotDepthsSize = lDoubleTempVector.size() + 1;
    par->snow_plot_depths.resize(lSoilPlotDepthsSize);
	for (size_t i=1; i<lSoilPlotDepthsSize; i++) {
		par->snow_plot_depths[i] = lDoubleTempVector[i-1] ;
	}

    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "GlacPlotDepths", geotop::input::gDoubleNoValue, false, 0, true) ;
    lSoilPlotDepthsSize = lDoubleTempVector.size() + 1;
    par->glac_plot_depths.resize(lSoilPlotDepthsSize);
	for (size_t i=1; i<lSoilPlotDepthsSize; i++) {
		par->glac_plot_depths[i] = lDoubleTempVector[i-1] ;
	}
#endif

	//output snow column
	if ((long)par->snow_plot_depths[1] != geotop::input::gDoubleNoValue) {
		m = par->snow_plot_depths.size();
	}else {
		m = par->max_snow_layers;
	}
	n = 6;
	
	osnw = (long*)malloc(n*sizeof(long));

#ifdef STAGED_FOR_REMOVING
	par->all_snow = (short)assignation_number(flog, 335, 0, keyword, num_param, num_param_components, 0., 0);
#else
    par->all_snow = (short)getDoubleValueWithDefault(lConfigStore, "SnowAll", geotop::input::gDoubleNoValue, false) ;
#endif

	if (par->all_snow == 1) {
		
		for (size_t i=0; i<n; i++) {
			osnw[i] = i;
		}
		
		nosnw = n;
		
	}else {

#ifdef STAGED_FOR_REMOVING
		cod = 323;
#endif
		
		for (size_t i=0; i<n; i++) {
			osnw[i] = -1;
		}
#ifdef STAGED_FOR_REMOVING
		for (size_t i=0; i<6; i++) {
			lColumnIndexJ = (long)assignation_number(flog, cod+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n)
                osnw[lColumnIndexJ-1] = i;
		}
#else
        lKeywordString.clear() ;
        lKeywordString += "DateSnow", "JulianDayFromYear0Snow", "TimeFromStartSnow",
            "PeriodSnow", "RunSnow", "IDPointSnow" ;

		for (size_t i=0; i<6; i++) {
			lColumnIndexJ = (long)getDoubleValueWithDefault(lConfigStore, lKeywordString[i], geotop::input::gDoubleNoValue, false) ;
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n)
                osnw[lColumnIndexJ-1] = i;
		}
#endif
		nosnw = 0;
		for (size_t i=0; i<n; i++) {
			if(osnw[i] > 0) nosnw = i+1;
		}
	}
		
	
	//output glacier column
	if ((long)par->glac_plot_depths[1] != geotop::input::gDoubleNoValue) {
		m = par->glac_plot_depths.size();
	}else {
		m = par->max_glac_layers;
	}
	n = 6 + 3*m + 1*par->max_glac_layers;
	oglc = (long*)malloc(n*sizeof(long));

#ifdef STAGED_FOR_REMOVING
	par->all_glac = (short)assignation_number(flog, 348, 0, keyword, num_param, num_param_components, 0., 0);
    cod = 335;
#else
    par->all_glac = (short)getDoubleValueWithDefault(lConfigStore, "GlacAll", geotop::input::gDoubleNoValue, false) ;
#endif

	if (par->all_glac == 1) {
		
		for (size_t i=0; i<n; i++) {
			oglc[i] = i;
		}
		
		noglc = n;
		
	}else {
				
		for (size_t i=0; i<n; i++) {
			oglc[i] = -1;
		}
#ifdef STAGED_FOR_REMOVING
		for (size_t i=0; i<6; i++) {
			lColumnIndexJ = (long)assignation_number(flog, cod+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n) oglc[lColumnIndexJ-1] = i;
		}
		for (size_t i=6; i<9; i++) {
			for (k=0; k<m; k++) {
				lColumnIndexJ = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (lColumnIndexJ>=1 && lColumnIndexJ<=n) oglc[lColumnIndexJ-1] = (i-6)*m + k + 6;
			}
		}
		for (size_t i=9; i<10; i++) {
			for (k=0; k<par->max_glac_layers; k++) {
				lColumnIndexJ = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (lColumnIndexJ>=1 && lColumnIndexJ<=n) oglc[lColumnIndexJ-1] = (i-9)*par->max_glac_layers + k + 6 + 3*m;
			}
		}
#else
        lKeywordString.clear() ;
        lKeywordString += "SnowAll", "DateGlac", "JulianDayFromYear0Glac", "TimeFromStartGlac", "PeriodGlac", "RunGlac" ;
        for(size_t i = 0 ; i < lKeywordString.size() ; i++)
        {
            lColumnIndexJ = getDoubleValueWithDefault(lConfigStore, "SnowAll", geotop::input::gDoubleNoValue, false) ;
            if (lColumnIndexJ>=1 && lColumnIndexJ<=n) oglc[lColumnIndexJ-1] = i;
        }

        lKeywordString.clear() ;
        lKeywordString += "IDPointGlac", "TempGlac", "IceContentGlac" ;
        for(size_t i = 0 ; i < lKeywordString.size() ; i++)
        {
            lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, lKeywordString[i], -1, false, 0, false) ;
			for (k=0; k<m; k++) {
				lColumnIndexJ = (size_t)lDoubleTempVector[k];
				if (lColumnIndexJ>=1 && lColumnIndexJ<=n) oglc[lColumnIndexJ-1] = (i)*m + k + 6;
			}
        }

        lKeywordString.clear() ;
        lKeywordString += "WatContentGlac" ;
        for(size_t i = 0 ; i < lKeywordString.size() ; i++)
        {
            lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, lKeywordString[i], -1, false, 0, false) ;
			for (k=0; k<par->max_glac_layers; k++) {
				lColumnIndexJ = (size_t)lDoubleTempVector[k];
				if (lColumnIndexJ>=1 && lColumnIndexJ<=n) oglc[lColumnIndexJ-1] = (i-9)*par->max_glac_layers + k + 6 + 3*m;
			}
        }
#endif
		noglc = 0;
		for (size_t i=0; i<n; i++) {
			if(oglc[i] > 0) noglc = i+1;
		}
	}
	
	//output soil column
	n = 6;
	osl = (long*)malloc(n*sizeof(long));

#ifdef STAGED_FOR_REMOVING
	par->all_soil = (short)assignation_number(flog, 355, 0, keyword, num_param, num_param_components, 0., 0);
#else
    par->all_soil = (short)getDoubleValueWithDefault(lConfigStore, "SoilAll", geotop::input::gDoubleNoValue, false) ;
#endif

	if (par->all_soil == 1) {
		
		for (size_t i=0; i<n; i++) {
			osl[i] = i;
		}
		
		nosl = n;
		
	}else {
		
		for (size_t i=0; i<n; i++) {
			osl[i] = -1;
		}
#ifdef STAGED_FOR_REMOVING
		for (size_t i=0; i<n; i++) {
			lColumnIndexJ = (long)assignation_number(flog, 349+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n) osl[lColumnIndexJ-1] = i;
		}
#else
        lKeywordString.clear() ;
        lKeywordString += "DateSoil", "JulianDayFromYear0Soil", "TimeFromStartSoil",
            "PeriodSoil", "RunSoil", "IDPointSoil" ;

		for (size_t i=0; i<n; i++) {
			lColumnIndexJ = (long)getDoubleValueWithDefault(lConfigStore, lKeywordString[i], geotop::input::gDoubleNoValue, false) ;
			if (lColumnIndexJ>=1 && lColumnIndexJ<=n) osl[lColumnIndexJ-1] = i;
		}
#endif
		nosl = 0;
		for (size_t i=0; i<n; i++) {
			if(osl[i] > 0) nosl = i+1;
		}
	}

#ifdef STAGED_FOR_REMOVING
	par->ric_cloud = (short)assignation_number(flog, 359, 0, keyword, num_param, num_param_components, 0., 0);
	par->vap_as_RH = (short)assignation_number(flog, 360, 0, keyword, num_param, num_param_components, 1., 0);
	par->vap_as_Td = (short)assignation_number(flog, 361, 0, keyword, num_param, num_param_components, 0., 0);
	par->ndivdaycloud = (long)assignation_number(flog, 362, 0, keyword, num_param, num_param_components, 3., 0);
	par->cast_shadow = (short)assignation_number(flog, 363, 0, keyword, num_param, num_param_components, 1., 0);
	par->wind_as_dir = (short)assignation_number(flog, 364, 0, keyword, num_param, num_param_components, 1., 0);
	par->wind_as_xy = (short)assignation_number(flog, 365, 0, keyword, num_param, num_param_components, 0., 0);
	par->snow_aging_vis = assignation_number(flog, 366, 0, keyword, num_param, num_param_components, 0.2, 0);
	par->snow_aging_nir = assignation_number(flog, 367, 0, keyword, num_param, num_param_components, 0.5, 0);
	par->DepthFreeSurface = assignation_number(flog, 368, 0, keyword, num_param, num_param_components, 0., 0);
	par->prec_as_intensity = assignation_number(flog, 369, 0, keyword, num_param, num_param_components, 0., 0);
#else
	par->ric_cloud = (short)getDoubleValueWithDefault(lConfigStore, "RicalculateCloudiness", geotop::input::gDoubleNoValue, false) ;
	par->vap_as_RH = (short)getDoubleValueWithDefault(lConfigStore, "DewTemperatureAsRH", geotop::input::gDoubleNoValue, false) ;
	par->vap_as_Td = (short)getDoubleValueWithDefault(lConfigStore, "RHAsDewTemperature", geotop::input::gDoubleNoValue, false) ;
	par->ndivdaycloud = (long)getDoubleValueWithDefault(lConfigStore, "NumberDayIntervalsToCalculateCloudiness", geotop::input::gDoubleNoValue, false) ;
	par->cast_shadow = (short)getDoubleValueWithDefault(lConfigStore, "CalculateCastShadow", geotop::input::gDoubleNoValue, false) ;
	par->wind_as_dir = (short)getDoubleValueWithDefault(lConfigStore, "WindAsSpeedAndDirection", geotop::input::gDoubleNoValue, false) ;
	par->wind_as_xy = (short)getDoubleValueWithDefault(lConfigStore, "WindAsWindXAndWindY", geotop::input::gDoubleNoValue, false) ;
	par->snow_aging_vis = getDoubleValueWithDefault(lConfigStore, "SnowAgingCoeffVis", geotop::input::gDoubleNoValue, false) ;
	par->snow_aging_nir = getDoubleValueWithDefault(lConfigStore, "SnowAgingCoeffNIR", geotop::input::gDoubleNoValue, false) ;
	par->DepthFreeSurface = getDoubleValueWithDefault(lConfigStore, "DepthFreeSurfaceAtTheBoundary", geotop::input::gDoubleNoValue, false) ;
	par->prec_as_intensity = getDoubleValueWithDefault(lConfigStore, "PrecAsIntensity", geotop::input::gDoubleNoValue, false) ;
#endif

    par->linear_interpolation_meteo.resize(nmeteo_stations + 1);
    
#ifdef STAGED_FOR_REMOVING
	cod = 370;
	par->linear_interpolation_meteo[1] = (short)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (size_t i=2; i<=nmeteo_stations; i++) {
		par->linear_interpolation_meteo[i] = (short)assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->linear_interpolation_meteo[i-1], 0);
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "LinearInterpolation", 0., true, nmeteo_stations, false) ;
	for (size_t i=1; i<par->linear_interpolation_meteo.size(); i++) {
		par->linear_interpolation_meteo[i] = lDoubleTempVector[i-1] ;
	}
#endif

#ifdef STAGED_FOR_REMOVING
	cod = 371;
    par->output_vertical_distances = (short)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	if(par->point_sim != 1){
		if(par->output_vertical_distances == 1){
			printf("Only for point simulations the parameter %s can be assigned to 1, layers are defined vertically\n",keyword[cod].c_str());
			fprintf(flog, "Only for point simulations the parameter %s can be assigned to 1, layers are defined vertically\n",keyword[cod].c_str());
		}
	}
#else
    par->output_vertical_distances = (short)getDoubleValueWithDefault(lConfigStore, "OutputDepthsVertical", geotop::input::gDoubleNoValue, false) ;
	if(par->point_sim != 1){
		if(par->output_vertical_distances == 1){
			printf("Only for point simulations the parameter OutputDepthsVertical can be assigned to 1, layers are defined vertically\n");
			fprintf(flog, "Only for point simulations the parameter OutputDepthsVertical can be assigned to 1, layers are defined vertically\n");
		}
	}
#endif

#ifdef STAGED_FOR_REMOVING
	par->upwindblowingsnow = (short)assignation_number(flog, 372, 0, keyword, num_param, num_param_components, 0., 0);

	par->UpdateK = (short)assignation_number(flog, 373, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->ContRecovery = assignation_number(flog, 374, 0, keyword, num_param, num_param_components, 0., 0);
	par->flag1D = (short)assignation_number(flog, 375, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->k_to_ksat = assignation_number(flog, 376, 0, keyword, num_param, num_param_components, 0., 0);
	par->RunIfAnOldRunIsPresent = (short)assignation_number(flog, 377, 0, keyword, num_param, num_param_components, 1., 0);
	
	par->max_courant_land_channel = assignation_number(flog, 378, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->min_dhsup_land_channel_out = assignation_number(flog, 379, 0, keyword, num_param, num_param_components, 1., 0);
#else
	par->upwindblowingsnow = (short)getDoubleValueWithDefault(lConfigStore, "UpwindBorderBlowingSnow", geotop::input::gDoubleNoValue, false) ;
    
	par->UpdateK = (short)getDoubleValueWithDefault(lConfigStore, "UpdateHydraulicConductivity", geotop::input::gDoubleNoValue, false) ;

	par->ContRecovery = getDoubleValueWithDefault(lConfigStore, "ContinuousRecovery", geotop::input::gDoubleNoValue, false) ;
	par->flag1D = (short)getDoubleValueWithDefault(lConfigStore, "ActualOrProjectedArea", geotop::input::gDoubleNoValue, false) ;
	
	par->k_to_ksat = getDoubleValueWithDefault(lConfigStore, "MinRatioKactualToKSat", geotop::input::gDoubleNoValue, false) ;
	par->RunIfAnOldRunIsPresent = (short)getDoubleValueWithDefault(lConfigStore, "RunIfAnOldRunIsPresent", geotop::input::gDoubleNoValue, false) ;
	
	par->max_courant_land_channel = getDoubleValueWithDefault(lConfigStore, "MaxCourantSupFlowChannelLand", geotop::input::gDoubleNoValue, false) ;
	par->min_dhsup_land_channel_out = getDoubleValueWithDefault(lConfigStore, "MinDiffSupWaterDepthChannelLand", geotop::input::gDoubleNoValue, false) ;
#endif

	par->Nl_spinup.resize(par->end_date.size() + 1);

#ifdef STAGED_FOR_REMOVING
	cod = 381;
	par->Nl_spinup[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 10000., 0);
	for (size_t i=2; i<par->end_date.size(); i++) {
		par->Nl_spinup[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Nl_spinup[i-1], 0);
	}
	if(par->Nl_spinup[1]<10000. && par->point_sim!=1){
		printf("You can use %s only if %s is set to 1\n",keyword[cod].c_str(),keyword[12].c_str());
		fprintf(flog,"You can use %s only if %s is set to 1\n",keyword[cod].c_str(),keyword[12].c_str());
		t_error("Not possible to continue");
	}
#else
    lDoubleTempVector = getDoubleVectorValueWithDefault(lConfigStore, "SpinUpLayerBottom", 0., true, par->end_date.size(), false) ;
	for (size_t i=1; i<par->Nl_spinup.size(); i++) {
		par->Nl_spinup[i] = lDoubleTempVector[i-1] ;
	}
	if(par->Nl_spinup[1]<10000. && par->point_sim!=1){
		printf("You can use SpinUpLayerBottom only if PointSim is set to 1\n");
		fprintf(flog,"You can use SpinUpLayerBottom only if PointSim is set to 1\n");
		t_error("Not possible to continue");
	}
#endif

#ifdef STAGED_FOR_REMOVING
	par->newperiodinit = (short)assignation_number(flog, 382, 0, keyword, num_param, num_param_components, 0., 0);
#else
	par->newperiodinit = (short)getDoubleValueWithDefault(lConfigStore, "InitInNewPeriods", geotop::input::gDoubleNoValue, false) ;
#endif
	if(par->newperiodinit != 0 && par->point_sim != 1){
		printf("You can use InitInNewPeriods only if PointSim is set to 1\n");
		fprintf(flog,"You can use InitInNewPeriods only if PointSim is set to 1\n");
		t_error("Not possible to continue");
	}

#ifdef STAGED_FOR_REMOVING
	par->k1 = assignation_number(flog, 383, 0, keyword, num_param, num_param_components, 0.484, 0);
	par->k2 = assignation_number(flog, 384, 0, keyword, num_param, num_param_components, 8., 0);
	par->Lozone = assignation_number(flog, 385, 0, keyword, num_param, num_param_components, 0.3, 0);
	par->alpha_iqbal = assignation_number(flog, 386, 0, keyword, num_param, num_param_components, 1.3, 0);
	par->beta_iqbal = assignation_number(flog, 387, 0, keyword, num_param, num_param_components, 0.1, 0);
		
	par->albedoSWin = (short)assignation_number(flog, 388, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->micro = (short)assignation_number(flog, 389, 0, keyword, num_param, num_param_components, 1., 0);
	par->EB = assignation_number(flog, 390, 0, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	par->Cair = assignation_number(flog, 391, 0, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);
	par->Tsup = assignation_number(flog, 392, 0, keyword, num_param, num_param_components, geotop::input::gDoubleNoValue, 0);

	par->Tair_default = assignation_number(flog, 393, 0, keyword, num_param, num_param_components, 5., 0);
	par->RH_default = assignation_number(flog, 394, 0, keyword, num_param, num_param_components, 70., 0)/100.;
	par->V_default = assignation_number(flog, 395, 0, keyword, num_param, num_param_components, par->Vmin, 0);
	par->Vdir_default = assignation_number(flog, 396, 0, keyword, num_param, num_param_components, 0., 0);
	par->IPrec_default = assignation_number(flog, 397, 0, keyword, num_param, num_param_components, 0., 0);
#else
	par->k1 = getDoubleValueWithDefault(lConfigStore, "KonzelmannA", geotop::input::gDoubleNoValue, false) ;
	par->k2 = getDoubleValueWithDefault(lConfigStore, "KonzelmannB", geotop::input::gDoubleNoValue, false) ;
	par->Lozone = getDoubleValueWithDefault(lConfigStore, "Lozone", geotop::input::gDoubleNoValue, false) ;
	par->alpha_iqbal = getDoubleValueWithDefault(lConfigStore, "AngstromAlpha", geotop::input::gDoubleNoValue, false) ;
	par->beta_iqbal = getDoubleValueWithDefault(lConfigStore, "AngstromBeta", geotop::input::gDoubleNoValue, false) ;
    
	par->albedoSWin = (short)getDoubleValueWithDefault(lConfigStore, "ConsiderAlbedoInSWin", geotop::input::gDoubleNoValue, false) ;
	
	par->micro = (short)getDoubleValueWithDefault(lConfigStore, "ConsiderMicrometeorology", geotop::input::gDoubleNoValue, false) ;
	par->EB = getDoubleValueWithDefault(lConfigStore, "SurfaceEnergyFlux", geotop::input::gDoubleNoValue, true) ;
	par->Cair = getDoubleValueWithDefault(lConfigStore, "ConvectiveHeatTransferCoefficient", geotop::input::gDoubleNoValue, true) ;
	par->Tsup = getDoubleValueWithDefault(lConfigStore, "SurfaceTemperature", geotop::input::gDoubleNoValue, true) ;
    
	par->Tair_default = getDoubleValueWithDefault(lConfigStore, "BaseAirTemperature", geotop::input::gDoubleNoValue, false) ;
	par->RH_default = getDoubleValueWithDefault(lConfigStore, "BaseRelativeHumidity", geotop::input::gDoubleNoValue, false)/100.;
	par->V_default = getDoubleValueWithDefault(lConfigStore, "BaseWindSpeed", geotop::input::gDoubleNoValue, false) ;
	par->Vdir_default = getDoubleValueWithDefault(lConfigStore, "BaseWindDirection", geotop::input::gDoubleNoValue, false) ;
	par->IPrec_default = getDoubleValueWithDefault(lConfigStore, "BaseIPrec", geotop::input::gDoubleNoValue, false) ;
#endif

#ifdef STAGED_FOR_REMOVING
	par->soil_type_bedr_default = (long)assignation_number(flog, 399, 0, keyword, num_param, num_param_components, 1., 0);
#else
	par->soil_type_bedr_default = (long)getDoubleValueWithDefault(lConfigStore, "DefaultSoilTypeBedrock", geotop::input::gDoubleNoValue, false) ;
#endif

    if(par->soil_type_bedr_default<1 || par->soil_type_bedr_default>par->nsoiltypes){
		fprintf(flog, "Error:  soil_type_bedr_default lower than 0 or higher than soil types numbers");
		printf("Error:  soil_type_bedr_default lower than 0 or higher than soil types numbers");
		t_error("Fatal Error! Geotop is closed. See failing report.");
	}

#ifdef STAGED_FOR_REMOVING
	par->minP_torestore_A = assignation_number(flog, 400, 0, keyword, num_param, num_param_components, 10., 0);
	par->snow_conductivity = (short)assignation_number(flog, 401, 0, keyword, num_param, num_param_components, 1., 0);
	par->snow_wind_compaction_1D = (short)assignation_number(flog, 402, 0, keyword, num_param, num_param_components, 0., 0);
#else
	par->minP_torestore_A = getDoubleValueWithDefault(lConfigStore, "MinPrecToRestoreFreshSnowAlbedo", geotop::input::gDoubleNoValue, false) ;
	par->snow_conductivity = (short)getDoubleValueWithDefault(lConfigStore, "SnowThermalConductivityPar", geotop::input::gDoubleNoValue, false) ;
	par->snow_wind_compaction_1D = (short)getDoubleValueWithDefault(lConfigStore, "WindCompaction1D", geotop::input::gDoubleNoValue, false) ;
#endif

	if(par->snow_wind_compaction_1D == 1) par->blowing_snow = 1;

#ifdef STAGED_FOR_REMOVING
	par->DDchannel = (short)assignation_number(flog, 403, 0, keyword, num_param, num_param_components, 1., 0);
	par->DDland = (short)assignation_number(flog, 404, 0, keyword, num_param, num_param_components, 1., 0);
#else
	par->DDchannel = (short)getDoubleValueWithDefault(lConfigStore, "DDChannel", geotop::input::gDoubleNoValue, false) ;
	par->DDland = (short)getDoubleValueWithDefault(lConfigStore, "DDLand", geotop::input::gDoubleNoValue, false) ;
#endif
    
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
std::vector<std::string> assign_string_parameter_v(FILE *f, long beg, long end, char **string_param, string keyword[]){
    
	long i;
	std::vector<std::string> a;

	for (i=0; i<end-beg; i++) {
        std::string lString = assignation_string(f, i+beg, keyword, string_param) ;
        a.push_back(lString);
	}
    
	return(a);
}

char **assign_string_parameter(FILE *f, long beg, long end, char **string_param, string keyword[]){

	long i;
	char **a;

	a = (char**)malloc((end-beg)*sizeof(char*));

	for (i=0; i<end-beg; i++) {
        std::string lString = assignation_string(f, i+beg, keyword, string_param);
        a[i] = strdup(lString.c_str()) ;
	}

	return(a);

}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
#ifdef STAGED_FOR_REMOVING
// overloaded function
double assignation_number(FILE *f, long i, long j, string keyword[], double **num_param, long *num_param_components, double default_value, short code_error){

	double a;
    
	if (j < num_param_components[i]) {
		a = num_param[i][j];
	}else {
		a = geotop::input::gDoubleNoValue;
	}

	if((long)a == geotop::input::gDoubleNoValue){
		if (code_error==0) {
			a = default_value;
			fprintf(f,"%s[%ld] = %e (default) \n", keyword[i].c_str(), j+1, a);
		}else{
			f = fopen(FailedRunFile.c_str(), "w");
			fprintf(f, "%s[%ld] not assigned\n", keyword[i].c_str(), j+1);
			fclose(f);
			t_error("Fatal Error!");
		}
	}else{
		fprintf(f,"%s[%ld] = %e \n", keyword[i].c_str(), j+1, a);
	}

    //TODO: removeme
    //std::cout << "DEBUG_KEYWORDS assignation_number: name(" << keyword[i] << "),Value(" << std::setprecision(12) << a << "),i(" << i << "),j(" << j << "),num_param_components(" << *num_param_components << "),default(" << std::setprecision(12) << default_value << ")" << std::endl ;

	return a;
}
#endif

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

std::string assignation_string(FILE *f, long i, string keyword[], char **string_param){

    std::string a;
	long j, dimstring = strlen(string_param[i]);

	for (j=0; j<dimstring; j++) {
		a.push_back(string_param[i][j]);
	}

    //TODO: removeme
    //std::cout << "DEBUG_KEYWORDS assignation_string: name(" << keyword[i] << "),Value("<< a << "),i(" << i << "),string_param[i](" << string_param[i] << ")" << std::endl ;
	fprintf(f,"%s = %s\n", keyword[i].c_str(), a.c_str());
	return(a);
}


/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_soil_parameters(std::string name, InitTools *IT, Soil *sl, long bed, FILE *flog){
	
	short ok;
	long i, j, k, n, nlines, nlinesprev;
    std::string temp;
	double **soildata;
    GeoTensor<double> old_sl_par;
    FILE *f ;
	
	//look if there is at least 1 soil file
i = 0;
	ok = 0;
	nlinesprev = -1;
	
	do{
		
		temp = namefile_i_we2(name, i+1);
		
		if (mio::IOUtils::fileExists(string(temp) + string(textfile))) {
			ok = 1;
			temp = namefile_i(name, i+1);
			nlines = count_lines(temp, 33, 44);
            
			if (nlinesprev >= 0 && nlines != nlinesprev){
				f = fopen(FailedRunFile.c_str(), "w");
				fprintf(f,"Error:: The file %s with soil paramaters has a number of layers %ld, which different from the numbers %ld of the other soil parameter files\n",temp.c_str(), nlines, nlinesprev);
				fprintf(f,"In GEOtop it is only possible to have the same number of layers in any soil parameter files\n");
				fclose(f);
				t_error("Fatal Error! GEOtop is closed. See failing report.");
			}
			nlinesprev = nlines;
		}else {
			if (i==0 && name != string_novalue) {
				f = fopen(FailedRunFile.c_str(), "w");
				fprintf(f,"Error:: Soil file %s not existing.\n", name.c_str());
				fclose(f);
				t_error("Fatal Error! GEOtop is closed. See failing report.");	
			}
		}
			
		i++;
		
	} while(ok == 0 && i < sl->pa.getDh()-1);
	
	if (ok == 1){
		
		//save sl->pa in a new doubletensor and deallocate
		old_sl_par.resize(sl->pa.getDh() + 1, sl->pa.getRh() + 1, sl->pa.getCh() + 1);
		for (i=1; i<sl->pa.getDh(); i++) {
			for (n=1; n<sl->pa.getRh(); n++) {
				for (j=1; j<sl->pa.getCh(); j++) {
					old_sl_par[i][n][j] = sl->pa(i,n,j);
				}
			}
		}
	//	free_doubletensor(sl->pa);

		//reallocate
		sl->pa.resize(old_sl_par.getDh(), old_sl_par.getRh(), nlines+1);

		for (i=1; i<sl->pa.getDh(); i++) {
			
			//read files
			temp = namefile_i_we2(name, i);
			
			if (mio::IOUtils::fileExists(string(temp) + string(textfile))) {
				temp = namefile_i(name, i);
				soildata = read_txt_matrix(temp, 33, 44, IT->soil_col_names, nsoilprop, &nlines, flog);
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
			
			//assign soildata to soil->pa
			for (n=1; n<=nsoilprop; n++) {
				for (j=1; j<sl->pa.getCh(); j++) { //j is the layer index
					sl->pa[i][n][j] = soildata[j-1][n-1];
				}
			}
			 
			//deallocate soildata
			for (j=0; j<nlines; j++) {
				free(soildata[j]);
			}
			free(soildata);
			
			//fix layer thickness
			n = jdz;
			for (j=1; j<sl->pa.getCh(); j++) { //j is the layer index
				if ((long)sl->pa[i][n][j] != geotop::input::gDoubleNoValue && (long)sl->pa[i][n][j] != number_absent) {
					if ( i > 1 && fabs( sl->pa[i][n][j] - sl->pa[i-1][n][j] ) > 1.E-5 )  {
						f = fopen(FailedRunFile.c_str(), "w");
						fprintf(f,"Error:: For soil type %ld it has been given a set of soil layer thicknesses different from the other ones.\n",i);
						fprintf(f,"In Geotop it is only possible to have the soil layer discretization in any soil parameter files.\n");
						fclose(f);
						t_error("Fatal Error! Geotop is closed. See failing report.");	
					}
				} else if (i == 1) {
					if (j < old_sl_par.getCh()) {
							sl->pa(i,n,j) = old_sl_par(i,n,j);
					} else {
							sl->pa(i,n,j) = sl->pa(i,n,j-1);
					}
				} else {
					sl->pa(i,n,j) = sl->pa(i-1,n,j);
				}
			}
			
			//all other variables
            for (n=1; n<=nsoilprop; n++) {
                if (n != jdz) {
                    //for (j=1; j<=sl->pa->nch; j++) { //j is the layer index
                    for (j=1; j<sl->pa.getCh(); j++) { //j is the layer index
                        if ((long)sl->pa[i][n][j] == geotop::input::gDoubleNoValue || (long)sl->pa[i][n][j] == number_absent) {
                            //if (j <= old_sl_par->nch) {
                            if (j < old_sl_par.getCh()) {
                                sl->pa[i][n][j] = old_sl_par[i][n][j];
                            }else {
                                sl->pa[i][n][j] = sl->pa[i][n][j-1];
                            }
                        }
                    }
                }
            }
			
			//field capacity and wilting point
			for (j=1; j<sl->pa.getCh(); j++){
				if( (long)sl->pa[i][jfc][j] == geotop::input::gDoubleNoValue){
					sl->pa[i][jfc][j] = teta_psi( (-1./3.)*1.E5/GTConst::GRAVITY, 0., sl->pa[i][jsat][j], sl->pa[i][jres][j], sl->pa[i][ja][j],
													 sl->pa[i][jns][j], 1.-1./sl->pa[i][jns][j], GTConst::PsiMin, sl->pa[i][jss][j]);
				}
				
				if( (long)sl->pa[i][jwp][j] == geotop::input::gDoubleNoValue){
					sl->pa[i][jwp][j] = teta_psi( -15.*1.E5/GTConst::GRAVITY, 0., sl->pa[i][jsat][j], sl->pa[i][jres][j], sl->pa[i][ja][j],
													 sl->pa[i][jns][j], 1.-1./sl->pa[i][jns][j], GTConst::PsiMin, sl->pa[i][jss][j]);
				}
			}
			
			//pressure
            ok = 1;
            for (j=1; j<sl->pa.getCh(); j++){
                if( (long)sl->pa[i][jpsi][j] == geotop::input::gDoubleNoValue) ok = 0;
            }


            if (ok == 1 )
            {
                if (IT->init_water_table_depth.size() <= i ) {

			IT->init_water_table_depth.resize(IT->init_water_table_depth.size() + 1);
                }
                IT->init_water_table_depth[i] = geotop::input::gDoubleNoValue;
            }
        }

		//free_doubletensor(old_sl_par);
			
	}
	
	//write on the screen the soil paramater
	fprintf(flog,"\n");
	k = (long)nmet;
	fprintf(flog,"Soil Layers: %ud\n",sl->pa.getCh()-1);
	for (i=1; i<sl->pa.getDh(); i++) {
		fprintf(flog,"-> Soil Type: %ld\n",i);
		for (n=1; n<=nsoilprop; n++) {
			fprintf(flog,"%s: ",keywords_char[k+n-1].c_str());
			for (j=1; j<sl->pa.getCh(); j++) {
				fprintf(flog,"%f(%.2e)",sl->pa[i][n][j],sl->pa[i][n][j]);
				if (j<sl->pa.getCh()-1) fprintf(flog,", ");
			}
			fprintf(flog,"\n");
		}
	}
	
	//bedrock
	old_sl_par.resize(1 + 1, IT->pa_bed.getRh() + 1, IT->pa_bed.getCh() + 1);
	for (n=1; n<sl->pa_bed.getRh(); n++) {
		for (j=1; j<sl->pa_bed.getCh(); j++) {
			old_sl_par[1][n][j] = IT->pa_bed[1][n][j];
		}
	}
	//free_doubletensor(sl->pa_bed);
	
	IT->pa_bed.resize(sl->pa.getDh(), sl->pa.getRh(), sl->pa.getCh());
	for (i=1; i<IT->pa_bed.getDh(); i++) {
		for (n=1; n<IT->pa_bed.getRh(); n++) {
			if (n == jdz) {
				for (j=1; j<IT->pa_bed.getCh(); j++) {
					IT->pa_bed(i,n,j) = sl->pa(1,n,j);
				}
			}else {
				for (j=1; j<IT->pa_bed.getCh(); j++) {
					if (j < old_sl_par.getCh()) {
						IT->pa_bed(i,n,j) = old_sl_par(1,n,j);
					}else {
						IT->pa_bed(i,n,j) = sl->pa_bed(i,n,j-1);
					}
				}
				for (j=1; j<IT->pa_bed.getCh(); j++) {
					if ( (long)IT->pa_bed(i,n,j) == geotop::input::gDoubleNoValue ) IT->pa_bed(i,n,j) = sl->pa(bed,n,j);
				}
			}
		}		
	}
	//free_doubletensor(old_sl_par);
	
	fprintf(flog,"\n");
	k = (long)nmet;
	fprintf(flog,"Soil Bedrock Layers: %ud\n",sl->pa.getCh()-1);
	for (i=1; i<IT->pa_bed.getDh(); i++) {
		fprintf(flog,"-> Soil Type: %ld\n",i);
		for (n=1; n<=nsoilprop; n++) {
			fprintf(flog,"%s: ",keywords_char[k+n-1].c_str());
			for (j=1; j<sl->pa.getCh(); j++) {
				fprintf(flog,"%f(%.2e)",IT->pa_bed[i][n][j],IT->pa_bed[i][n][j]);
				if(j < sl->pa.getCh()-1) fprintf(flog,", ");
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

short read_point_file(std::string name, char **key_header, Par *par, FILE *flog){

	GeoMatrix<double>  chkpt2;
	double **points;
	long nlines, n, j;
      std::string temp;

	if (mio::IOUtils::fileExists(string(name) + string(textfile))) {
		temp = name + std::string(textfile);
		printf("%s\n",temp.c_str());
		points = read_txt_matrix(temp, 34, 44, key_header, par->chkpt.getCols()-1, &nlines, flog);
				
		chkpt2.resize(par->chkpt.getRows() + 1, par->chkpt.getCols() + 1);

		chkpt2=par->chkpt;
	//	free_doublematrix(par->chkpt);
		
		par->chkpt.resize(nlines+1, chkpt2.getCols());
		for (n=1; n<=nlines; n++) {
			for (j=1; j< chkpt2.getCols(); j++) {
				par->chkpt[n][j] = points[n-1][j-1];
				if ( (long)par->chkpt[n][j] == geotop::input::gDoubleNoValue || (long)par->chkpt[n][j] == number_absent ) {
					if ( n < chkpt2.getRows() ) {
						par->chkpt[n][j] = chkpt2[n][j];
					}else {
						par->chkpt[n][j] = chkpt2[chkpt2.getRows()-1][j];
					}
				}
			}
			
			if ( par->point_sim != 1 ){
				if( (long)par->chkpt[n][ptX] == geotop::input::gDoubleNoValue || (long)par->chkpt[n][ptY] == geotop::input::gDoubleNoValue ) {
					par->state_pixel = 0;
				}
			}
			
			free(points[n-1]);
		}
				
	//	free_doublematrix(chkpt2);
		free(points);
	}
	
	if (par->point_sim != 1){
		for (n=1; n< par->chkpt.getRows(); n++) {
			if ( (long)par->chkpt[n][ptX] == geotop::input::gDoubleNoValue || (long)par->chkpt[n][ptY] == geotop::input::gDoubleNoValue) {
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
	
//short read_meteostations_file(LONGVECTOR *i, METEO_STATIONS *S, char *name, char **key_header, FILE *flog){
short read_meteostations_file(const GeoVector<long>& i, MeteoStations *S, std::string name, char **key_header, FILE *flog){
	double **M;
	long nlines, n, j, k;
    std::string temp;
		
	if (mio::IOUtils::fileExists(name + string(textfile))) {
		temp = name + textfile ;
		M = read_txt_matrix(temp, 33, 44, key_header, 8, &nlines, flog);
				
		for (j=1; j<i.size(); j++) {
			for (n=1; n<=nlines; n++) {
				if ((long)M[n-1][0] == i[j]) {
					for (k=1; k<8; k++) {
						if ((long)M[n-1][k] != geotop::input::gDoubleNoValue && (long)M[n-1][k] != number_absent) {
							if (k==1) {
								S->E[j] = M[n-1][k];
							}else if (k==2) {
								S->N[j] = M[n-1][k];
							}else if (k==3) {
								S->lat[j] = M[n-1][k];
							}else if (k==4) {
								S->lon[j] = M[n-1][k];
							}else if (k==5) {
								S->Z[j] = M[n-1][k];
							}else if (k==6) {
								S->sky[j] = M[n-1][k];
							}else if (k==7) {
								S->ST[j] = M[n-1][k];
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
