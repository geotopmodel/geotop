
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225-9 'Moab' - 24 Aug 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225-9 'Moab'
 
 GEOtop 1.225-9 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
  
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#include "input.h"
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//! Subroutine which reads input data, performs  geomporphological analisys and allocates data
#ifdef USE_HPC
void get_all_input(long argc, char *argv[], TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet,
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times, WORKAREA *rankArea)
#else
void get_all_input(long argc, char *argv[], TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet, 
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times)
#endif
 
{
	
	FILE *flog, *f;
	DOUBLEMATRIX *M;
	INIT_TOOLS *IT;

	short a, success, added_JDfrom0=0, added_wind_xy=0, added_wind_dir=0, added_cloud=0, added_Tdew=0, added_RH=0, added_Pint=0, old=0;
	long l, r, c, i, ist, j, n, sy, num_cols, num_lines, day, month, year, hour, minute;
	double z, th_oversat, JD, k_snowred, maxSWE, SWE, D, cosslope, **matrix;
	char *temp, *name, **temp2;	
	char rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	
	IT=(INIT_TOOLS *)malloc(sizeof(INIT_TOOLS));
	
	if(!argv[1]){
		WORKING_DIRECTORY=get_workingdirectory();
	} else if (argc==2) {
    // modified by Emanuele Cordano on Aug 2011
		WORKING_DIRECTORY=assign_string(argv[1]);
	}	else {
	// modified by Emanuele Cordano on Aug 2011
	//	WORKING_DIRECTORY=assign_string(read_option_string(argc,argv,"-wpath",".",0)); // assign_string(argv[1]); // MODIFY HERE EC
		WORKING_DIRECTORY=get_workingdirectory();
	}
	
	//add "/" if it is missing
	if (WORKING_DIRECTORY[strlen(WORKING_DIRECTORY)-1] != 47) {
		temp = assign_string(WORKING_DIRECTORY);
		free(WORKING_DIRECTORY);
		WORKING_DIRECTORY = join_strings(temp, "/");
		free(temp);
	}	
	
	logfile = join_strings(WORKING_DIRECTORY, logfile_name);
	flog = fopen(logfile, "w");
	
	printf("STATEMENT:\n");
	printf("\n");	 
	printf("GEOtop 1.225-9 'Moab' - 24 Aug 2012\n\n");	 
	//printf("Copyright (c), 2012 - Stefano Endrizzi \n\n");	 
	printf("GEOtop 1.225-9 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
	printf("WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
 	//printf("Riccardo Rigon is acknowledged as he coded the Fluidturtle routines (GPL Licenced), which are used in GEOtop 1.225-9 'Moab'.\n");
	//printf("Riccardo Rigon is also acknowledged as first founder of the GEOtop model in 1997.\n");
	//printf("Riccardo Rigon and his research group are acknowledged as GEOtop 1.225-9 'Moab' uses most of their modelling achievements.\n");
	//printf("John Pomeroy is acknowledged as he freely provided the Prairie Blowing Snow Model Code.\n");
	//printf("Glen Liston and Kelly Elder are acknowledged as they freely provided their Micromet code, from which the routines that distribute wind-air temperature-relative humidity-precipitation in GEOtop 1.225-9 'Moab' are derived.\n");
	//printf("However, the routine that distributes the meteo data in this GEOtop version is named Meteodistr and it significantly differs from Micromet.\n\n");
	//printf("If you have satisfactorily used the code, please acknowledge the authors.\n");	
	printf("\nWORKING DIRECTORY: %s\n",WORKING_DIRECTORY);
	printf("\nLOGFILE: %s\n",logfile);
	
	fprintf(flog,"STATEMENT:\n");
	fprintf(flog,"\n");	 
	fprintf(flog,"GEOtop 1.225-9 'Moab' - 24 Aug 2012\n\n");	 
	//fprintf(flog,"Copyright (c), 2012 - Stefano Endrizzi \n\n");	 
	fprintf(flog,"GEOtop 1.225-9 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
	fprintf(flog,"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
	//fprintf(flog,"Riccardo Rigon is acknowledged as he coded the Fluidturtle routines (GPL Licenced), which are used in GEOtop 1.225-9 'Moab'.\n");
	//fprintf(flog,"Riccardo Rigon is also acknowledged as first founder of the GEOtop model in 1997.\n");
	//fprintf(flog,"Riccardo Rigon and his research group are acknowledged as GEOtop 1.225-9 'Moab' uses most of their modelling achievements.\n");
	//fprintf(flog,"John Pomeroy is acknowledged as he freely provided the Prairie Blowing Snow Model Code.\n");
	//fprintf(flog,"Glen Liston and Kelly Elder are acknowledged as they freely provided their Micromet code, from which the routines that distribute wind-air temperature-relative humidity-precipitation in GEOtop 1.225-9 'Moab' are derived.\n");
	//fprintf(flog,"However, the routine that distributes the meteo data in this GEOtop version is named Meteodistr and it significantly differs from Micromet.\n\n");
	//fprintf(flog,"If you have satisfactorily used the code, please acknowledge the authors.\n");	 	
	fprintf(flog,"\nWORKING DIRECTORY: %s\n",WORKING_DIRECTORY);
	
	//reads the parameters in __control_parameters
	temp = join_strings(WORKING_DIRECTORY, program_name);
	success = read_inpts_par(par, land, times, sl, met, IT, temp, flog); 
	free(temp);
	
	//correct state pixel
	par->Tzrun = 0;
	par->wzrun = 0;
	par->Tzmaxrun = 0;
	par->Tzminrun = 0;
	par->wzmaxrun = 0;
	par->wzminrun = 0;	
	par->dUzrun = 0;
	par->SWErun = 0;	
	
	if(strcmp(files[fTrun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->Tzrun = 1;
	}
	if(strcmp(files[fwrun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->wzrun = 1;
	}
	if(strcmp(files[fTmaxrun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->Tzmaxrun = 1;

	}
	if(strcmp(files[fwmaxrun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->wzmaxrun = 1;
	}
	if(strcmp(files[fTminrun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->Tzminrun = 1;
	}
	if(strcmp(files[fwminrun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->wzminrun = 1;
	}
	if(strcmp(files[fdUrun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->dUzrun = 1;
	}
	if(strcmp(files[fSWErun] , string_novalue) != 0){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->SWErun = 1;
	}
	if (par->newperiodinit == 2 && (par->Tzrun == 0 || par->wzrun == 0)){
		f = fopen(FailedRunFile, "w");
		fprintf(f, "Error: You have to assign a name to the Tzrun and wzrun files\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");		
	}	
	
	//soil parameters
	success = read_soil_parameters(files[fspar], IT, sl, par->soil_type_bedr_default, flog);	
	Nl=sl->pa->nch;	
	
	//pointlist files
	success = read_point_file(files[fpointlist], IT->point_col_names, par, flog);
	
	//max time that the simulation is supposed to model
	max_time = 0.;
	for (i=1; i<=par->init_date->nh; i++) {
		max_time += par->run_times->co[i]*(par->end_date->co[i] - par->init_date->co[i])*86400.;//seconds
	}
	
	//Recovering
	par->delay_day_recover = 0.0;
	par->n_ContRecovery = 0;
	
	if(par->recover > 0){
		if(par->saving_points->nh < par->recover){
			f = fopen(FailedRunFile, "w");
			fprintf(f, "Error: recover index higher than the length of the saving points vector");
			fclose(f);
			t_error("Fatal Error! Geotop is closed. See failing report (1).");
		}
		par->delay_day_recover = par->saving_points->co[par->recover];
	}
	
	//Continuous Recovering
	if (par->RunIfAnOldRunIsPresent != 1) {
		if (existing_file_wext(SuccessfulRunFile,"") == 1){
			temp = join_strings(SuccessfulRunFile, ".old");
			rename(SuccessfulRunFile, temp);
			free(temp);
			temp = join_strings(FailedRunFile, ".old");
			rename(FailedRunFile, temp);
			free(temp);	
			f = fopen(FailedRunFile, "w");
			fprintf(f, "This simulation has successfully reached the end, you cannot recover it.\n");
			fclose(f);
			t_error("Fatal Error! Geotop is closed. See failing report.");
		}
		if (existing_file_wext(FailedRunFile,"") == 1){
			temp = join_strings(SuccessfulRunFile, ".old");
			rename(SuccessfulRunFile, temp);
			free(temp);
			temp = join_strings(FailedRunFile, ".old");
			rename(FailedRunFile, temp);
			free(temp);	
			f = fopen(FailedRunFile, "w");
			fprintf(f, "This simulation has failed, you cannot recover it.\n");
			fclose(f);
			t_error("Fatal Error! Geotop is closed. See failing report.");
		}
	}
	
	temp = join_strings(SuccessfulRunFile, ".old");
	rename(SuccessfulRunFile, temp);
	free(temp);
	temp = join_strings(FailedRunFile, ".old");
	rename(FailedRunFile, temp);
	free(temp);	
	
	if (par->ContRecovery > 0 && par->recover == 0) {
		
		if (existing_file_woext(files[rsux]) != 1){
			temp = join_strings(files[rsux], ".old");
			if (existing_file_woext(temp) == 1) old = 1;
			free(temp);
		}

		if (old == 1) {
			name = join_strings(files[rtime], ".old");
		}else {
			name = assign_string(files[rtime]);
		}
		
		if (existing_file_wext(name, textfile) == 1){
			temp = join_strings(name, textfile);
			matrix = read_txt_matrix_2(temp, 33, 44, 8, &num_lines);
			par->delay_day_recover = matrix[0][1];
			i_run0 = (long)matrix[0][3];
			i_sim0 = (long)matrix[0][4];
			cum_time = matrix[0][5];
			elapsed_time_start = matrix[0][6];
			par->n_ContRecovery = (long)matrix[0][7]+1;
			for (i=0; i<num_lines; i++) {
				free(matrix[i]);
			}
			free(matrix);
			free(temp);
		}
		
		free(name);
	}
	
	
	
	//Time indices	
	convert_JDfrom0_JDandYear(par->init_date->co[i_sim0], &JD, &year);
	convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);
		
	/****************************************************************************************************/
	/*! Reading of the Input files:                                                                     */
	/****************************************************************************************************/
	
	if(par->point_sim!=1){  //distributed simulation
		read_inputmaps(top, land, sl, par, IT, flog);
	}else{
		read_optionsfile_point(par, top, land, sl, times, IT, flog);
	}
	
	Nr=top->Z0->nrh;
	Nc=top->Z0->nch;
	
	par->total_pixel=0;
	par->total_area=0.;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if ((long)land->LC->co[r][c]!=number_novalue){
				par->total_pixel++;
				if(par->point_sim != 1){
					par->total_area += (UV->U->co[1]*UV->U->co[2])/cos(top->slope->co[r][c]*Pi/180.);
				}else {
					par->total_area += (UV->U->co[1]*UV->U->co[2]);
				}
			}
		}
	}
		
	top->Z = find_Z_of_any_layer(top->Z0, top->slope, land->LC, sl, par->point_sim);
	
	top->Jdown = new_longmatrix(par->total_pixel, 4);
	top->Qdown = new_doublematrix(par->total_pixel, 4);
	for (i=1; i<=par->total_pixel; i++) {
		for (j=1; j<=4; j++) {
			top->Jdown->co[i][j] = i;
			top->Qdown->co[i][j] = 0.;
		}
	}
	
	fprintf(flog,"Valid pixels: %ld\n",par->total_pixel);
	fprintf(flog,"Number of nodes: %ld\n",(Nl+1)*par->total_pixel);
	fprintf(flog,"Novalue pixels: %ld\n",(Nr*Nc-par->total_pixel));
	fprintf(flog,"Basin area: %f km2\n",(double)par->total_pixel*UV->U->co[1]*UV->U->co[2]/1.E6);
	
	printf("\nValid pixels: %ld\n",par->total_pixel);
	printf("Number of nodes: %ld\n",(Nl+1)*par->total_pixel);
	printf("Novalue pixels: %ld\n",(Nr*Nc-par->total_pixel));
	printf("Basin area: %f km2\n\n",(double)par->total_pixel*UV->U->co[1]*UV->U->co[2]/1.E6);
	
	/****************************************************************************************************/
	//Reading of RAIN data file,	METEO data file, 	and CLOUD data file
	
	num_cols = (long)nmet;
	
	//meteo data
	met->data=(double***)malloc(met->st->E->nh*sizeof(double**));
	//number of line of meteo data
	met->numlines=(long*)malloc(met->st->E->nh*sizeof(long));
	//meteo variables for the current instant
	met->var=(double**)malloc(met->st->E->nh*sizeof(double*));
	//horizon for meteo stations
	met->horizon=(double***)malloc(met->st->E->nh*sizeof(double**));
	//number of line in the horizon file
	met->horizonlines=(long*)malloc(met->st->E->nh*sizeof(long));
	//line of met->data used (stored in memory to avoid from searching from the first line)
	met->line_interp_WEB=(long*)malloc(met->st->E->nh*sizeof(long));
	met->line_interp_Bsnow=(long*)malloc(met->st->E->nh*sizeof(long));
	met->line_interp_WEB_LR=0;
	met->line_interp_Bsnow_LR=0;	
	
	success = read_meteostations_file(met->imeteo_stations, met->st, files[fmetstlist], IT->meteostations_col_names, flog);
	par->usemeteoio=0;
	par->use_meteoio_meteodata=0;
	par->use_meteoio_cloud=0;
//#ifdef USE_METEOIO
	par->usemeteoio=1;
	par->use_meteoio_meteodata=1;
	par->use_meteoio_cloud=1;
	if(par->usemeteoio==1) meteoio_init();
//#endif
	for(i=1;i<=met->st->E->nh;i++){
				
		if (met->imeteo_stations->co[1] != number_novalue) {
			ist = met->imeteo_stations->co[i];
		}else {
			ist = i;
		}

		//initialize 
		met->line_interp_WEB[i-1] = 0;
		met->line_interp_Bsnow[i-1] = 0;

		//allocate var
		met->var[i-1] = (double*)malloc(num_cols*sizeof(double));
		
		//read horizon
		met->horizon[i-1] = read_horizon(1, ist, files[fhormet], IT->horizon_col_names, &num_lines, flog);
		met->horizonlines[i-1] = num_lines;
		
		//filename
		if (strcmp(files[fmet], string_novalue) != 0){
						
			//read matrix
			temp=namefile_i(files[fmet], ist);
			met->data[i-1] = read_txt_matrix(temp, 33, 44, IT->met_col_names, nmet, &num_lines, flog);
			
			if ((long)met->data[i-1][0][iDate12] == number_absent && (long)met->data[i-1][0][iJDfrom0] == number_absent) {
				f = fopen(FailedRunFile, "w");
				fprintf(f, "Error:: Date Column missing in file %s\n",temp);
				fclose(f);
				t_error("Fatal Error! Geotop is closed. See failing report (2).");
			}
			met->numlines[i-1] = num_lines;
			
			//fixing dates: converting times in the same standard time set for the simulation and fill JDfrom0
			added_JDfrom0 = fixing_dates(ist, met->data[i-1], par->ST, met->st->ST->co[i], met->numlines[i-1], iDate12, iJDfrom0);
			
			check_times(ist, met->data[i-1], met->numlines[i-1], iJDfrom0);
			
			//find clouds
			if(strcmp(IT->met_col_names[itauC], string_novalue) != 0){
				if((long)met->data[i-1][0][itauC] == number_absent || par->ric_cloud == 1){
					added_cloud = fill_meteo_data_with_cloudiness(met->data[i-1], met->numlines[i-1], met->horizon[i-1], met->horizonlines[i-1], 
						met->st->lat->co[i], met->st->lon->co[i], par->ST, met->st->Z->co[i], met->st->sky->co[i], 0.0, par->ndivdaycloud, par->dem_rotation,
						par->Lozone, par->alpha_iqbal, par->beta_iqbal, 0.);
				}
			}
						
			//calcululate Wx and Wy if Wspeed and direction are given
			if (par->wind_as_xy == 1) {
				added_wind_xy = fill_wind_xy(met->data[i-1], met->numlines[i-1], iWs, iWdir, iWsx, iWsy, IT->met_col_names[iWsx], IT->met_col_names[iWsy]);
			}

			//calcululate Wspeed and direction if Wx and Wy are given
			if (par->wind_as_dir == 1) {
				added_wind_dir = fill_wind_dir(met->data[i-1], met->numlines[i-1], iWs, iWdir, iWsx, iWsy, IT->met_col_names[iWs], IT->met_col_names[iWdir]);
			}
			
			//find Tdew
			if(par->vap_as_Td == 1){
			   added_Tdew = fill_Tdew(i, met->st->Z, met->data[i-1], met->numlines[i-1], iRh, iT, iTdew, IT->met_col_names[iTdew], par->RHmin);
			}
			
			//find RH
			if(par->vap_as_RH == 1){
				added_RH = fill_RH(i,  met->st->Z, met->data[i-1], met->numlines[i-1], iRh, iT, iTdew, IT->met_col_names[iRh]);
			}
			
			//find Prec Intensity
			if ( par->linear_interpolation_meteo->co[i] == 1 && (long)met->data[i-1][0][iPrec] != number_absent) {
				f = fopen(FailedRunFile, "w");
				fprintf(f,"Meteo data for station %ld contain precipitation as volume, but Linear Interpolation is set. This is not possible, the precipitation data are removed.\n",i);
				fprintf(f,"If you want to use precipitation as volume, you cannot set keyword LinearInterpolation at 1.\n");
				fclose(f);
				t_error("Fatal Error! Geotop is closed. See failing report (3).");
			}
					
			if(par->prec_as_intensity == 1){
				added_Pint = fill_Pint(i, met->data[i-1], met->numlines[i-1], iPrec, iPrecInt, iJDfrom0, IT->met_col_names[iPrecInt]);
			}
				
			//rewrite completed files 
			rewrite_meteo_files(met->data[i-1], met->numlines[i-1], IT->met_col_names, temp, added_JDfrom0, added_wind_xy, added_wind_dir, added_cloud, added_Tdew, added_RH, added_Pint);
				
			//calcululate Wx and Wy if Wspeed and direction are given
			if (par->wind_as_xy != 1) {
				added_wind_xy = fill_wind_xy(met->data[i-1], met->numlines[i-1], iWs, iWdir, iWsx, iWsy, IT->met_col_names[iWsx], IT->met_col_names[iWsy]);
			}
			
			//find Prec Intensity
			if(par->prec_as_intensity != 1){
				added_Pint = fill_Pint(i, met->data[i-1], met->numlines[i-1], iPrec, iPrecInt, iJDfrom0, IT->met_col_names[iPrecInt]);
			}
			
			//find Tdew
			if(par->vap_as_Td != 1){
				added_Tdew = fill_Tdew(i, met->st->Z, met->data[i-1], met->numlines[i-1], iRh, iT, iTdew, IT->met_col_names[iTdew], par->RHmin);
			}
			
			free(temp);
			
		}else {
			
			fprintf(flog, "Warning: File meteo not in the list, meteo data not read, used default values\n");
			printf("Warning: File meteo not in the list, meteo data not read, used default values\n");

			met->data[i-1] = (double**)malloc(2*sizeof(double*));
			
			for (n=1; n<=2; n++) {
				met->data[i-1][n-1] = (double*)malloc(num_cols*sizeof(double));
				for (j=1; j<=nmet; j++) {
					met->data[i-1][n-1][j-1] = (double)number_absent;
				}
			}
			
			met->data[i-1][0][iJDfrom0] = 0.;
			met->data[i-1][1][iJDfrom0] = 1.E10;
			
		}
	}
	
		
	//read LAPSE RATES FILE  
	
	if(strcmp(files[fLRs] , string_novalue) != 0){   //s stands for string
		
		if(existing_file_wext(files[fLRs], textfile)==0) printf("Lapse rate file unavailable. Check input files. If you do not have a lapse rate file, remove its name and keywords from input file\n");
		temp = join_strings(files[fLRs], textfile);	
		met->LRs = read_txt_matrix(temp, 33, 44, IT->lapserates_col_names, nlstot, &num_lines, flog);
		free(temp);
		met->LRsnr = num_lines;
		par->LRflag=1;
		printf("\nLapse rate file read\n");
		
	}else{		
		
		par->LRflag=0;
		
	}
	
	n = (long)nlstot;
	met->LRv = (double*)malloc(n*sizeof(double));
	met->LRd = (double*)malloc(n*sizeof(double));
	for( i=0; i<nlstot; i++){
		met->LRv[i] = (double)number_novalue;
		if (i == ilsTa) {
			met->LRd[i] = LapseRateTair;
		}else if (i == ilsTdew) {
			met->LRd[i] = LapseRateTdew;
		}else if (i == ilsPrec) {
			met->LRd[i] = LapseRatePrec;
		}else {
			met->LRd[i] = 0.0;
		}
	}
	
	
	//FIND A STATION WITH SHORTWAVE RADIATION DATA
	met->tau_cl_map=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	initialize_doublematrix(met->tau_cl_map, (double)number_novalue);

	met->tau_cl_av_map=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	initialize_doublematrix(met->tau_cl_av_map, (double)number_novalue);

	met->tau_cl_map_yes=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
	initialize_shortmatrix(met->tau_cl_map_yes, (short)number_novalue);

	met->tau_cl_av_map_yes=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
	initialize_shortmatrix(met->tau_cl_av_map_yes, (short)number_novalue);

	// vector defining which meteo station has the SW radiation information
	met->st->flag_SW_meteoST=new_shortvector(met->st->Z->nh);
	initialize_shortvector(met->st->flag_SW_meteoST,number_novalue);

	met->st->tau_cloud_av_yes_meteoST=new_shortvector(met->st->Z->nh);
	initialize_shortvector(met->st->tau_cloud_av_yes_meteoST,number_novalue);

	met->st->tau_cloud_yes_meteoST=new_shortvector(met->st->Z->nh);
	initialize_shortvector(met->st->tau_cloud_yes_meteoST,number_novalue);

	// find number of meteo stations measuring SW radiation
	met->nstsrad=0;// meteo station ID (1...n) to use for the SW
	met->numstsrad=0;// counter of meteo stations containing SW info
	a=0;// counter
	do{
		a++;
		if( (long)met->data[a-1][0][iSW]!=number_absent || ((long)met->data[a-1][0][iSWb]!=number_absent && (long)met->data[a-1][0][iSWd]!=number_absent ) ) {
			if(met->nstsrad==0) met->nstsrad=a;// records the first station that has cloud info
			met->numstsrad++;
			met->st->flag_SW_meteoST->co[a]=1;
		}
	}while(a<met->st->Z->nh);//}while(met->nstsrad<met->st->Z->nh && a==0);
	if(met->numstsrad==0){
		printf("WARNING: NO shortwave radiation measurements available\n");
		fprintf(flog,"WARNING: NO shortwave radiation measurements available\n");
		met->nstsrad=1;
	}else{
		printf("Shortwave radiation measurements available from %ld stations: will be used data from station %ld\n",met->numstsrad,met->nstsrad);
		fprintf(flog,"Shortwave radiation measurements available from %ld stations: will be used data from station %ld\n",met->numstsrad,met->nstsrad);
	}
	
	//FIND A STATION WITH CLOUD DATA
	met->nstcloud=0; // meteo station ID (1...n) to use for the cloudiness
	met->numstcloud=0;// counter of meteo stations containing cloud info
	a=0;// counter
	do{
		a++;
		if( (long)met->data[a-1][0][iC]!=number_absent || (long)met->data[a-1][0][itauC]!=number_absent ) {
			if(met->nstcloud==0) met->nstcloud=a;// records the first station that has cloud info
			met->numstcloud++;// increase the counter of meteo stations measuring cloud
		}
	}while(a<met->st->Z->nh);//while(met->nstcloud<met->st->Z->nh && a==0);
	if(met->numstcloud==0){
		printf("WARNING: NO cloudiness measurements available\n");
		fprintf(flog,"WARNING: NO cloudiness measurements available\n");
		met->nstcloud=1;
	}else{
		printf("Cloudiness measurements from %ld stations: will be used data from station %ld\n",met->numstcloud,met->nstcloud);
		fprintf(flog,"Cloudiness measurements from station %ld\n",met->nstcloud);
	}
	
	met->st->tau_cloud_av_meteoST=new_doublevector(met->st->Z->nh);
	initialize_doublevector(met->st->tau_cloud_av_meteoST,(double)number_novalue);
	met->st->tau_cloud_meteoST=new_doublevector(met->st->Z->nh);
	initialize_doublevector(met->st->tau_cloud_meteoST,(double)number_novalue);


	//FIND A STATION WITH LONGWAVE RADIATION DATA
	met->nstlrad=0;
	do{
		met->nstlrad++;
		a=0;
		if( (long)met->data[met->nstlrad-1][0][iLWi]!=number_absent) a=1;
	}while(met->nstlrad<met->st->Z->nh && a==0);
	if(a==0){
		printf("WARNING: NO longwave radiation measurements available\n");
		fprintf(flog,"WARNING: NO longwave radiation measurements available\n");
	}else{
		printf("Longwave radiation measurements from station %ld\n",met->nstlrad);
		fprintf(flog,"Longwave radiation measurements from station %ld\n",met->nstlrad);
	}
	//i.show details on checkpoints
	fprintf(flog,"\nCHECKPOINTS:\n");
	fprintf(flog,"ID,r,c,Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-]\n");
	for(i=1;i<=par->rc->nrh;i++){
		r=par->rc->co[i][1];
		c=par->rc->co[i][2];
		fprintf(flog,"%ld,%ld,%ld,%f,%d,%ld,%f,%f,%f\n",i,r,c,top->Z0->co[r][c],(short)land->LC->co[r][c],sl->type->co[r][c],top->slope->co[r][c],top->aspect->co[r][c],top->sky->co[r][c]);
	}

	//i.show meteo stations
	fprintf(flog,"\nMETEO STATIONS:\n");
	fprintf(flog,"ID,East[m],North[m],Lat[deg],Lon[deg],Elev[m a.s.l.],Sky[deg],Stand_Time[h],WindSensHeight[m],TempSensHeight[m]\n");
	for(r=1;r<=met->st->E->nh;r++){
		fprintf(flog,"%ld,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",r,met->st->E->co[r],met->st->N->co[r],met->st->lat->co[r], met->st->lon->co[r],
				met->st->Z->co[r],met->st->sky->co[r],met->st->ST->co[r],met->st->Vheight->co[r],met->st->Theight->co[r]);
	}
	
	//FIND A STATION WITH SURFACE TEMPERATURE
	met->nstTs=0;
	do{
		met->nstTs++;
		a=0;
		if( (long)met->data[met->nstTs-1][0][iTs]!=number_absent) a=1;
	}while(met->nstTs<met->st->Z->nh && a==0);
	if(a==0){
		printf("WARNING: NO Surface temperature measurements available\n");
		fprintf(flog,"WARNING: NO Surface temperature measurements available\n");
	}else{
		printf("Surface temperature measurements from station %ld\n",met->nstTs);
		fprintf(flog,"Surface temperature measurements from station %ld\n",met->nstTs);
	}
	
	/****************************************************************************************************/
	//read INCOMING DISCHARGE

	met->qinv = (double*)malloc(2*sizeof(double));
	met->qinv[0] = 0.;
	met->qinv[1] = 0.;
	
	if(strcmp(files[fqin] , string_novalue) != 0){   
		temp = join_strings(files[fqin], textfile);	
		temp2 = (char **)malloc(2*sizeof(char*));
		temp2[0] = assign_string("Time");
		temp2[1] = assign_string("Qx");
		met->qins = read_txt_matrix(temp, 33, 44, temp2, 2, &num_lines, flog);
		free(temp);
		free(temp2[0]);
		free(temp2[1]);
		free(temp2);
		met->qinsnr = num_lines;
		par->qin = 1;
		met->qinline = 0;
		printf("\nIncoming discharge file read\n");
	}else{		
		par->qin = 0;
	}
	
	
	/****************************************************************************************************/
	/*! Completing the several time-indipendent input variables with the data of input-files:           */
	/****************************************************************************************************/
	/****************************************************************************************************/
	// Completing of "land" (of the type LAND):  
	
	//Initialize matrix of shadow
	land->shadow=new_shortmatrix(Nr,Nc);
	initialize_shortmatrix(land->shadow,0);/* initialized as if it was always NOT in shadow*/
	
	//Check that there aren't cell with an undefined land use value
	z = 0.;
	l = 0;
	do{	
		l++;
		z += sl->pa->co[1][jdz][l];
	}while(l<Nl && z < z_transp);
	land->root_fraction=new_doublematrix(par->n_landuses, l);
	initialize_doublematrix(land->root_fraction, 0.0);
	
	//check vegetation variable consistency
	if(jHveg!=jdHveg+jHveg-1) t_error("Vegetation variables not consistent");
	if(jz0thresveg!=jdz0thresveg+jHveg-1) t_error("Vegetation variables not consistent");
	if(jz0thresveg2!=jdz0thresveg2+jHveg-1) t_error("Vegetation variables not consistent");
	if(jLSAI!=jdLSAI+jHveg-1) t_error("Vegetation variables not consistent");
	if(jcf!=jdcf+jHveg-1) t_error("Vegetation variables not consistent");
	if(jdecay0!=jddecay0+jHveg-1) t_error("Vegetation variables not consistent");
	if(jexpveg!=jdexpveg+jHveg-1) t_error("Vegetation variables not consistent");
	if(jroot!=jdroot+jHveg-1) t_error("Vegetation variables not consistent");
	if(jrs!=jdrs+jHveg-1) t_error("Vegetation variables not consistent");
	
	//variables used to assign vegetation properties that change with time	
	num_cols = jdvegprop + 1;
	land->vegpars=(double ***)malloc(par->n_landuses*sizeof(double**));
	land->vegparv=(double **)malloc(par->n_landuses*sizeof(double*));
	land->NumlinesVegTimeDepData=(long*)malloc(par->n_landuses*sizeof(long));
	
	land->vegpar=new_doublevector(jdvegprop);
	
	par->vegflag=new_shortvector(par->n_landuses);	
	initialize_shortvector(par->vegflag, 0);
	
	//time dependent vegetation parameters	
	for(i=1;i<=par->n_landuses;i++){
				
		if(strcmp(files[fvegpar] , string_novalue) != 0){   //s stands for string
			
			temp = namefile_i_we2(files[fvegpar], i);
						
			if(existing_file_wext(temp, textfile)==1){ 
				
				printf("There is a specific vegetation parameter file for land cover type = %ld\n",i);
				free(temp);
				temp = namefile_i(files[fvegpar], i);
				land->vegpars[i-1] = read_txt_matrix_2(temp, 33, 44, num_cols, &num_lines);
				free(temp);
				land->NumlinesVegTimeDepData[i-1] = num_lines;
				par->vegflag->co[i]=1;
				
			}else{
				
				free(temp);
				printf("There is NOT a specific vegetation parameter file for land cover type = %ld\n",i);
				
			}
			
		}
		
		land->vegparv[i-1]=(double*)malloc(num_cols*sizeof(double));
		for(j=0; j<num_cols; j++){
			land->vegparv[i-1][j] = (double)number_novalue;
		}
		
		//z0 (convert in m)
		land->ty->co[i][jz0]*=0.001;
		
		//find root fraction
		root(land->root_fraction->nch, land->ty->co[i][jroot], 0.0, sl->pa->co[1][jdz], land->root_fraction->co[i]);
		
		//error messages
		for(l=1;l<=met->st->Z->nh;l++){
			if(0.001*land->ty->co[i][jHveg]>met->st->Vheight->co[l] || 0.001*land->ty->co[i][jHveg]>met->st->Theight->co[l]){
				f = fopen(FailedRunFile, "w");
				fprintf(f, "hc:%f m, zmu:%f m, zmt:%f m - hc must be lower than measurement height - land cover %ld, meteo station %ld\n",
							0.001*land->ty->co[i][jHveg],met->st->Vheight->co[l],met->st->Theight->co[l],i,l);
				fclose(f);
				t_error("Fatal Error! Geotop is closed. See failing report (5).");
			}
		}
	}
	
	//******************************************	
	//file with time steps		
	
	if(strcmp(files[ftsteps] , string_novalue) != 0){   
		
		temp=join_strings(files[ftsteps],textfile);
		times->Dt_matrix = read_txt_matrix_2(temp, 33, 44, max_cols_time_steps_file+1, &num_lines);
		free(temp);
		times->numlinesDt_matrix = num_lines;
		par->tsteps_from_file=1;
		
	}else{
		
		par->tsteps_from_file=0;
		
	}
	
	/****************************************************************************************************/
	/*! Filling up of the struct "channel" (of the type CHANNEL):                                        */
	
	/*The number of channel-pixel are counted:*/
	i=0;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if (top->pixel_type->co[r][c]>=10) i++;
		}
	}
	fprintf(flog,"Channel pixels: %ld\n",i);
	par->total_channel = i;
	
	//allocate channel vectors/matrixes
	if(i==0) i=1;	
	
	cnet->Vout = 0.;
	
	cnet->r=new_longvector(i);
	initialize_longvector(cnet->r, 0);
	
	cnet->c=new_longvector(i);
	initialize_longvector(cnet->c, 0);
	
	cnet->ch=new_longmatrix(Nr,Nc);
	initialize_longmatrix(cnet->ch, 0);
	
	cnet->ch_down=new_longvector(i);
	initialize_longvector(cnet->ch_down, 0);
	
	cnet->length=new_doublevector(i);
	initialize_doublevector(cnet->length, 0.);
	
	cnet->Vsup=new_doublevector(i);
	initialize_doublevector(cnet->Vsup, 0.);
	
	cnet->Vsub=new_doublevector(i);
	initialize_doublevector(cnet->Vsub, 0.);
		
	cnet->h_sup=new_doublevector(i);
	initialize_doublevector(cnet->h_sup, 0.);
	
	cnet->soil_type = new_longvector(i);
	initialize_longvector(cnet->soil_type, par->soil_type_chan_default);
	
	if(par->total_channel>1) enumerate_channels(cnet, land->LC, top->pixel_type, top->Z0, top->slope, number_novalue);
				
	cnet->ch3 = (long**)malloc((Nl+1)*sizeof(long*));
	for (l=0; l<=Nl; l++) {
		cnet->ch3[l] = (long*)malloc((i+1)*sizeof(long));
	}
	
	cnet->lch = new_longmatrix( (Nl+1)*i, 2);
	initialize_longmatrix(cnet->lch, 0);
	
	lch3_cont(cnet->ch3, cnet->lch, Nl, par->total_channel);
	
	
	/****************************************************************************************************/	
	//Cont for Richards 3D
	
	n = Fminlong(par->Nl_spinup->co[i_sim0],Nl);
	
	//3D
	top->i_cont=(long ***)malloc((n+1)*sizeof(long**));
	for(l=0;l<=n;l++){
		top->i_cont[l]=(long **)malloc((Nr+1)*sizeof(long*));
		for(r=1;r<=Nr;r++){
			top->i_cont[l][r]=(long *)malloc((Nc+1)*sizeof(long));
		}
	}
	
	top->lrc_cont=new_longmatrix( (n+1)*par->total_pixel, 3);
	initialize_longmatrix(top->lrc_cont, 0);
	
	i_lrc_cont(land->LC, top->i_cont, top->lrc_cont, n, Nr, Nc);
	
	//2D
	top->j_cont=(long **)malloc((Nr+1)*sizeof(long*));
	for (r=1; r<=Nr; r++) {
		top->j_cont[r]=(long *)malloc((Nc+1)*sizeof(long));
		for (c=1; c<=Nc; c++) {
			top->j_cont[r][c] = 0;
		}
	}
	
	top->rc_cont=new_longmatrix(par->total_pixel, 2);
	initialize_longmatrix(top->rc_cont, 0);
	
	j_rc_cont(land->LC, top->j_cont, top->rc_cont, Nr, Nc);

	//plotted points
	if(par->state_pixel == 1){
		par->jplot = new_longvector(par->total_pixel);
		initialize_longvector(par->jplot, 0);
		for (i=1; i<=par->total_pixel; i++) {
			for (j=1; j<=par->rc->nrh; j++) {
				if (top->rc_cont->co[i][1] == par->rc->co[j][1] && top->rc_cont->co[i][2] == par->rc->co[j][2]) {
					par->jplot->co[i] = j;
				}
			}
		}
	}	
	
	//BEDROCK (adjusting soil properties)
	set_bedrock(IT, sl, cnet, par, top, land->LC, flog);	
		
	/****************************************************************************************************/
	/*! Completing of the initialization of SOIL structure                               */
	/****************************************************************************************************/
	
	sl->SS = (SOIL_STATE *)malloc(sizeof(SOIL_STATE));
	initialize_soil_state(sl->SS, par->total_pixel, Nl);
	
	sl->VS = (STATE_VEG *)malloc(sizeof(STATE_VEG));
	initialize_veg_state(sl->VS, par->total_pixel);
		
	sl->th=new_doublematrix(Nl,par->total_pixel);
	initialize_doublematrix(sl->th,(double)number_novalue);
	
	sl->Ptot=new_doublematrix(Nl,par->total_pixel);
	initialize_doublematrix(sl->Ptot,(double)number_novalue);

	sl->ET=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(sl->ET,0.);
	
	if(par->output_soil_bin == 1){
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0){
			sl->T_av_tensor=new_doublematrix(Nl,par->total_pixel);
			initialize_doublematrix(sl->T_av_tensor,0.);
		}
		
		if(strcmp(files[ficeav] , string_novalue) != 0){
			sl->thi_av_tensor=new_doublematrix(Nl,par->total_pixel);
			initialize_doublematrix(sl->thi_av_tensor,0.);
		}
		
		if(strcmp(files[fliqav] , string_novalue) != 0){
			sl->thw_av_tensor=new_doublematrix(Nl,par->total_pixel);
			initialize_doublematrix(sl->thw_av_tensor,0.);
		}
		
		if(strcmp(files[fpnet] , string_novalue) != 0){
			sl->Pnetcum=new_doublevector(par->total_pixel);
			initialize_doublevector(sl->Pnetcum,0.);
		}	
		
		if(strcmp(files[fevap] , string_novalue) != 0){
			sl->ETcum=new_doublevector(par->total_pixel);
			initialize_doublevector(sl->ETcum,0.);
		}			
		
	}
						
	if(existing_file(files[fwt0]) == 0){
		
		for (i=1; i<=par->total_pixel; i++) {

			r = top->rc_cont->co[i][1];
			c = top->rc_cont->co[i][2];
			
			sy=sl->type->co[r][c];
						
			if ((long)IT->init_water_table_depth->co[sy] != number_novalue) {
				z = 0.;
				sl->SS->P->co[0][i] = -IT->init_water_table_depth->co[sy]*cos(top->slope->co[r][c]*Pi/180.);
				for(l=1;l<=Nl;l++){
					z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
					sl->SS->P->co[l][i] = sl->SS->P->co[0][i] + z;
					z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
				}
			}else {
				for(l=1;l<=Nl;l++){
					sl->SS->P->co[l][i] = sl->pa->co[sy][jpsi][l];
				}
			}
		}
		
	}else {
		
		M = read_map(2, files[fwt0], land->LC, UV, (double)number_novalue);
		
		for (i=1; i<=par->total_pixel; i++) {
			
			r = top->rc_cont->co[i][1];
			c = top->rc_cont->co[i][2];

			sy=sl->type->co[r][c];
										
			z = 0.;
			sl->SS->P->co[0][i] = -M->co[r][c]*cos(top->slope->co[r][c]*Pi/180.);
			for(l=1;l<=Nl;l++){
				z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
				sl->SS->P->co[l][i] = sl->SS->P->co[0][i] + z;
				z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
			}
		}

		free_doublematrix(M);
	}
			
	for (i=1; i<=par->total_pixel; i++) {
		
		r = top->rc_cont->co[i][1];
		c = top->rc_cont->co[i][2];
		
		sy=sl->type->co[r][c];

		for(l=1;l<=Nl;l++){
			
			sl->SS->T->co[l][i]=sl->pa->co[sy][jT][l];
					
			sl->Ptot->co[l][i] = sl->SS->P->co[l][i];
			sl->th->co[l][i] = teta_psi(sl->SS->P->co[l][i], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
										   sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
					
			th_oversat = Fmax( sl->SS->P->co[l][i] , 0.0 ) * sl->pa->co[sy][jss][l];
			sl->th->co[l][i] -= th_oversat;
											
			if(sl->SS->T->co[l][i]<=Tfreezing){
				
				//Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
				sl->SS->thi->co[l][i] = sl->th->co[l][i] - teta_psi(Psif(sl->SS->T->co[l][i]), 0.0, sl->pa->co[sy][jsat][l], 
																		sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 
																		1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
										
				//if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
				if(sl->SS->thi->co[l][i]<0) sl->SS->thi->co[l][i]=0.0;
				
				//Psi is updated taking into account the freezing
				sl->th->co[l][i] -= sl->SS->thi->co[l][i];

				sl->SS->P->co[l][i] = psi_teta(sl->th->co[l][i] + th_oversat, sl->SS->thi->co[l][i], sl->pa->co[sy][jsat][l],
											   sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 
											   1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
			}
		}
	}
	
				
	if(par->state_pixel == 1){
		if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot = new_doublematrix_0(par->rc->nrh, Nl);
		if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thizplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot = new_doublematrix(par->rc->nrh, Nl);
		
		for (i=1; i<=par->rc->nrh; i++) {
			r = top->rc_cont->co[i][1];
			c = top->rc_cont->co[i][2];
			j = top->j_cont[r][c];
			for(l=1;l<=Nl;l++){
				if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot->co[i][l] = sl->SS->T->co[l][j];
				if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] = sl->SS->T->co[l][j];
				if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot->co[i][l] = sl->th->co[l][j];
				if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] = sl->th->co[l][j];
				if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thizplot->co[i][l] = sl->SS->thi->co[l][j];
				if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot->co[i][l] = sl->SS->thi->co[l][j];
				if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot->co[i][l] = sl->Ptot->co[l][j];
			}
			for(l=0;l<=Nl;l++) {
				if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot->co[i][l] = sl->SS->P->co[l][j];
			}
		}

	}
	
	if(par->delay_day_recover > 0){
		
		assign_recovered_tensor_vector(old, par->recover, files[riceg], sl->SS->thi, top->rc_cont, par, land->LC);
		assign_recovered_tensor_vector(old, par->recover, files[rTg], sl->SS->T, top->rc_cont, par, land->LC);
		assign_recovered_tensor_vector(old, par->recover, files[rpsi], sl->SS->P, top->rc_cont, par, land->LC);
		
		assign_recovered_map_vector(old, par->recover, files[rwcrn], sl->VS->wrain, top->rc_cont, par, land->LC);
		assign_recovered_map_vector(old, par->recover, files[rwcsn], sl->VS->wsnow, top->rc_cont, par, land->LC);
		assign_recovered_map_vector(old, par->recover, files[rTv], sl->VS->Tv, top->rc_cont, par, land->LC);
		
	}

	
	//channel soil
	cnet->SS = (SOIL_STATE *)malloc(sizeof(SOIL_STATE));
	initialize_soil_state(cnet->SS, cnet->r->nh, Nl);
	
	cnet->th = new_doublematrix(Nl, cnet->r->nh);
	
	cnet->ET = new_doublematrix(Nl, cnet->r->nh);	
	initialize_doublematrix(cnet->ET, 0.0);
	
	cnet->Kbottom = new_doublevector(cnet->r->nh);
	initialize_doublevector(cnet->Kbottom, 0.0);
	
	for(j=1;j<=par->total_channel;j++){
				
		sy=cnet->soil_type->co[j];
		r=cnet->r->co[j];
		c=cnet->c->co[j];
		
		cnet->SS->P->co[0][j] = sl->SS->P->co[0][top->j_cont[r][c]] + par->depr_channel;
		
		for(l=1;l<=Nl;l++){
			cnet->SS->P->co[l][j] = sl->Ptot->co[l][top->j_cont[r][c]] + par->depr_channel;
		}
		
		for(l=1;l<=Nl;l++){
			
			cnet->SS->T->co[l][j]=sl->pa->co[sy][jT][l];
						
			cnet->th->co[l][j] = teta_psi(cnet->SS->P->co[l][j], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], 
										  sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], 
										  PsiMin, sl->pa->co[sy][jss][l]);
			
			th_oversat = Fmax( cnet->SS->P->co[l][j] , 0.0 ) * sl->pa->co[sy][jss][l];
			cnet->th->co[l][j] -= th_oversat;
			
			if(cnet->SS->T->co[l][j]<=Tfreezing){
				//Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
				cnet->SS->thi->co[l][j] = cnet->th->co[l][j] - teta_psi(Psif(cnet->SS->T->co[l][j]), 0.0, sl->pa->co[sy][jsat][l], 
																		sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 
																		1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
				
				//if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
				if(cnet->SS->thi->co[l][j]<0) cnet->SS->thi->co[l][j]=0.0;
				
				//Psi is updated taking into account the freezing
				cnet->th->co[l][j] -= cnet->SS->thi->co[l][j];
				cnet->SS->P->co[l][j] = psi_teta(cnet->th->co[l][j] + th_oversat, cnet->SS->thi->co[l][j], 
												 sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
												 sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
			}
		}
	}
	
	if(par->delay_day_recover > 0 && par->total_channel > 0){
		assign_recovered_tensor_channel(old, par->recover, files[rpsich], cnet->SS->P, cnet->r, cnet->c, top->Z0);
		assign_recovered_tensor_channel(old, par->recover, files[ricegch], cnet->SS->thi, cnet->r, cnet->c, top->Z0);
		assign_recovered_tensor_channel(old, par->recover, files[rTgch], cnet->SS->T, cnet->r, cnet->c, top->Z0);
		
		for(i=1; i<=par->total_channel; i++){
			for (l=1; l<=Nl; l++) {
				sy = cnet->soil_type->co[i];
				cnet->th->co[l][i] = teta_psi(Fmin(cnet->SS->P->co[l][i], psi_saturation(cnet->SS->thi->co[l][i], sl->pa->co[sy][jsat][l], 
						sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l])), 
						cnet->SS->thi->co[l][i], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
						sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);				
			}
		}
		
	}
	
	//WRITE INITIAL CONDITION
	write_output_headers(met->st->Z->nh, times, wat, par, top, land, sl, egy, snow, glac);
	
	if(par->state_pixel == 1){
		for(j=1;j<=par->rc->nrh;j++){
			
			if(par->output_vertical_distances == 1){
				r = par->rc->co[j][1];
				c = par->rc->co[j][2];				
				cosslope = cos( Fmin(max_slope, top->slope->co[r][c]) * Pi/180. );
			}else {
				cosslope = 1.;
			}

			write_soil_output(j, par->IDpoint->co[j], par->init_date->co[1], par->init_date->co[1], JD, day, month, year, hour, minute, par->soil_plot_depths, sl, par, PsiMin, cosslope);
		}
	}
	
	//z boundary condition
	for(l=1;l<=Nl;l++){
		par->Zboundary -= sl->pa->co[1][jdz][l];
	}
	
	if(par->Zboundary < 0){		
		f = fopen(FailedRunFile, "w");
		fprintf(f, "Z at which 0 annual temperature takes place is not lower than the soil column\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report (6).");		
	}
	
	par->Zboundary *= 1.E-3;	//convert in [m]
	
	/****************************************************************************************************/
	/*! Initialization of the struct "egy" (of the type ENERGY):*/
	
	if(par->output_surfenergy_bin == 1){
		if(strcmp(files[fradnet] , string_novalue) != 0){
			egy->Rn_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Rn_mean, 0.);
			egy->Rn = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fradLWin] , string_novalue) != 0){
			egy->LWin_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->LWin_mean, 0.);
			egy->LWin = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fradLW] , string_novalue) != 0){
			egy->LW_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->LW_mean, 0.);
		}
		if(strcmp(files[fradLW] , string_novalue) != 0 || strcmp(files[fradnet] , string_novalue) != 0){
			egy->LW = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fradSW] , string_novalue) != 0){
			egy->SW_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->SW_mean, 0.);
		}
		if(strcmp(files[fradSW] , string_novalue) != 0 || strcmp(files[fradnet] , string_novalue) != 0){
			egy->SW = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fLE] , string_novalue) != 0){
			egy->ET_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->ET_mean, 0.);
			egy->LE = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fH] , string_novalue) != 0){
			egy->H_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->H_mean, 0.);
			egy->H = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fG] , string_novalue) != 0){
			egy->SEB_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->SEB_mean, 0.);
			egy->G = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fTs] , string_novalue) != 0){
			egy->Ts_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Ts_mean, 0.);
			egy->Ts = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fradSWin] , string_novalue) != 0){
			egy->Rswdown_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Rswdown_mean, 0.);
			egy->SWin = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fradSWinbeam] , string_novalue) != 0){
			egy->Rswbeam_mean = new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Rswbeam_mean, 0.);
			egy->SWinb = new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fshadow] , string_novalue) != 0){
			egy->nDt_shadow=new_longvector(par->total_pixel);
			egy->nDt_sun=new_longvector(par->total_pixel);
			initialize_longvector(egy->nDt_shadow,0);
			initialize_longvector(egy->nDt_sun,0); 	
			egy->shad=new_shortvector(par->total_pixel);
		}
	}
	
	egy->sun = (double*)malloc(12*sizeof(double));
	
	if(times->JD_plots->nh > 1){
		if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			egy->Hgplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Hgplot, 0.);
			egy->Hgp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0){
			egy->Hvplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Hvplot, 0.);
			egy->Hvp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			egy->LEgplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->LEgplot, 0.);
			egy->LEgp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0){
			egy->LEvplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->LEvplot, 0.);
			egy->LEvp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pSWin] , string_novalue) != 0){
			egy->SWinplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->SWinplot, 0.);
			egy->SWinp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			egy->SWgplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->SWgplot, 0.);
			egy->SWgp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pSWv] , string_novalue) != 0){
			egy->SWvplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->SWvplot, 0.);
			egy->SWvp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pLWin] , string_novalue) != 0){
			egy->LWinplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->LWinplot, 0.);
			egy->LWinp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			egy->LWgplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->LWgplot, 0.);
			egy->LWgp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pLWv] , string_novalue) != 0){
			egy->LWvplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->LWvplot, 0.);
			egy->LWvp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pTs] , string_novalue) != 0){
			egy->Tsplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Tsplot, 0.);
			egy->Tsp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pTg] , string_novalue) != 0){
			egy->Tgplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Tgplot, 0.);
			egy->Tgp=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[pTv] , string_novalue) != 0){
			egy->Tvplot=new_doublevector(par->total_pixel);
			initialize_doublevector(egy->Tvplot, 0.);
		}
	}
		
	
	//vectors used in energy_balance()	
	egy->Tgskin_surr = new_doublematrix(Nr, Nc);
	initialize_doublematrix(egy->Tgskin_surr, 0.);
	
	egy->SWrefl_surr = new_doublematrix(Nr, Nc);
	initialize_doublematrix(egy->SWrefl_surr, 0.);
	
	egy->Dlayer = new_doublevector( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->liq = new_doublevector( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->ice = new_doublevector( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->Temp = new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->deltaw = new_doublevector( Nl + par->max_snow_layers + par->max_glac_layers );
	
	egy->SWlayer = new_doublevector0( par->max_snow_layers + 1 );
	
	egy->soil_transp_layer = new_doublevector(land->root_fraction->nch);
	initialize_doublevector(egy->soil_transp_layer, 0.);
	
	egy->dFenergy = new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->udFenergy = new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers - 1);
	egy->Kth0=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers - 1 );
	egy->Kth1=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers - 1);
	egy->Fenergy=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->Newton_dir=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->T0=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->T1=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
	egy->Tstar=new_doublevector(Nl); //soil temperature at which freezing begins
	egy->THETA=new_doublevector(Nl);	//water content (updated in the iterations)
	
	//allocate vector	of soil layer contributions to evaporation (up to z_evap)
	z = 0.;
	l = 0;
	do{	
		l++;
		z += sl->pa->co[1][jdz][l];
	}while(l<Nl && z < z_evap);
	egy->soil_evap_layer_bare = new_doublevector(l);
	egy->soil_evap_layer_veg = new_doublevector(l);
	initialize_doublevector(egy->soil_evap_layer_bare, 0.);
	initialize_doublevector(egy->soil_evap_layer_veg, 0.);
	
	printf("Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare->nh);
	printf("Soil water transpires from the first %ld layers\n",egy->soil_transp_layer->nh);
	fprintf(flog,"Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare->nh);
	fprintf(flog,"Soil water transpires from the first %ld layers\n",egy->soil_transp_layer->nh);
	
	/****************************************************************************************************/
	/*! Completing of the struct "water" (of the type WATER) */
	
	wat->Voutlandsub = 0.;
	wat->Voutlandsup = 0.;
	wat->Voutbottom = 0.;
	
	/* Initialization of wat->Pnet (liquid precipitation that reaches the sl surface in mm):*/
	wat->Pnet=new_doublematrix(Nr,Nc);
	initialize_doublematrix(wat->Pnet,0.0);

	/* Initialization of wat->HN (new snow precipitation in mm):*/
	wat->HN=new_doublematrix(Nr,Nc);
	initialize_doublematrix(wat->HN,0.0);
		
	/* Initialization of wat->PrecTot (total precipitation (rain+snow) precipitation):*/
	wat->PrecTot=new_doublematrix(Nr,Nc);
	initialize_doublematrix(wat->PrecTot,par->IPrec_default);
	
	/* Initialization of the matrices with the output of total precipitation and interception:*/
	if (par->output_meteo_bin == 1 && strcmp(files[fprec], string_novalue) != 0){
		wat->PrTOT_mean=new_doublevector(par->total_pixel);
		initialize_doublevector(wat->PrTOT_mean, 0.);
		wat->PrSNW_mean=new_doublevector(par->total_pixel);
		initialize_doublevector(wat->PrSNW_mean, 0.);
		wat->Pt=new_doublevector(par->total_pixel);
		wat->Ps=new_doublevector(par->total_pixel);
	}
	
	wat->h_sup=new_doublevector(par->total_pixel);
	initialize_doublevector(wat->h_sup, 0.);
		
	/****************************************************************************************************/
	/*! Initialization of the struct "snow" (of the type SNOW):*/
	
	/***************************************************************************************************/
	snow->S=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
	allocate_and_initialize_statevar_3D(snow->S, (double)number_novalue, par->max_snow_layers, Nr, Nc);

	//initial snow depth 
	if( strcmp(files[fsn0] , string_novalue) != 0 && strcmp(files[fswe0] , string_novalue) != 0 ){
		printf("Initial condition on snow depth from file %s\n",files[fsn0]);
		M=read_map(2, files[fsn0], land->LC, UV, (double)number_novalue);
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				snow->S->Dzl->co[1][r][c] = M->co[r][c];
			}
		}
		free_doublematrix(M);
		
		printf("Initial condition on snow water equivalent from file %s\n",files[fswe0]);
		M=read_map(2, files[fswe0], land->LC, UV, (double)number_novalue);
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				snow->S->w_ice->co[1][r][c] = M->co[r][c];
			}
		}
		free_doublematrix(M);
		
	}else if( strcmp(files[fsn0] , string_novalue) != 0 ){
		printf("Initial condition on snow depth from file %s\n",files[fsn0]);
		M=read_map(2, files[fsn0], land->LC, UV, (double)number_novalue);
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				snow->S->Dzl->co[1][r][c] = M->co[r][c];
			}
		}
		free_doublematrix(M);
		
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				if ((long)land->LC->co[r][c] != number_novalue) snow->S->w_ice->co[1][r][c] = snow->S->Dzl->co[1][r][c]*IT->rhosnow0/rho_w;
			}
		}
		
	}else if( strcmp(files[fswe0] , string_novalue) != 0 ){
		printf("Initial condition on snow water equivalent from file %s\n",files[fswe0]);
		M=read_map(2, files[fswe0], land->LC, UV, (double)number_novalue);
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				snow->S->w_ice->co[1][r][c] = M->co[r][c];
			}
		}
		free_doublematrix(M);
		
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				if ((long)land->LC->co[r][c] != number_novalue) snow->S->Dzl->co[1][r][c] = snow->S->w_ice->co[1][r][c]*rho_w/IT->rhosnow0;
			}
		}
		
	}else {	
		
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				if ((long)land->LC->co[r][c] != number_novalue){
					snow->S->w_ice->co[1][r][c] = IT->swe0;
					snow->S->Dzl->co[1][r][c] = IT->swe0*rho_w/IT->rhosnow0;
				}
			}
		}
	}

	
	//Optional reading of snow age in the whole basin
	if( strcmp(files[fsnag0] , string_novalue) != 0 ){
		printf("Snow age initial condition from file %s\n",files[fsnag0]+1);
		snow->age = read_map_vector(2, files[fsnag0], land->LC, UV, (double)number_novalue, top->rc_cont);
	}else{
		snow->age = new_doublevector(par->total_pixel);
		initialize_doublevector(snow->age, IT->agesnow0);
	}
	

	if(times->JD_plots->nh > 1){
		if(strcmp(files[pD] , string_novalue) != 0) {
			snow->Dplot=new_doublevector(par->total_pixel);
			initialize_doublevector(snow->Dplot, 0.);
		}
	}

	if(par->blowing_snow==1){
		
		snow->S_for_BS=(STATEVAR_1D *)malloc(sizeof(STATEVAR_1D));
		allocate_and_initialize_statevar_1D(snow->S_for_BS, (double)number_novalue, par->max_snow_layers);
				
		snow->change_dir_wind=new_longvector(Fmaxlong(Nr,Nc));
		
		snow->Qtrans=new_doublematrix(Nr,Nc);
		snow->Qsub=new_doublematrix(Nr,Nc);
		initialize_doublematrix(snow->Qtrans, 0.0);
		initialize_doublematrix(snow->Qsub, 0.0);
		
		snow->Qsalt=new_doublematrix(Nr,Nc);
		
		snow->Nabla2_Qtrans=new_doublematrix(Nr,Nc);
		snow->Qsub_x=new_doublematrix(Nr,Nc);
		snow->Qsub_y=new_doublematrix(Nr,Nc);
		snow->Qtrans_x=new_doublematrix(Nr,Nc);
		snow->Qtrans_y=new_doublematrix(Nr,Nc);
				
		if(par->output_snow_bin == 1){
			snow->Wtrans_plot=new_doublematrix(Nr,Nc);
			snow->Wsubl_plot=new_doublematrix(Nr,Nc);
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if((long)land->LC->co[r][c]==number_novalue){
						snow->Wtrans_plot->co[r][c]=(double)number_novalue;		
						snow->Wsubl_plot->co[r][c]=(double)number_novalue;	
						
					}else{
						snow->Wtrans_plot->co[r][c]=0.0;					
						snow->Wsubl_plot->co[r][c]=0.0;			
						
					}
				}
			}
		}
	}
	
	if(par->output_snow_bin == 1){
		if(strcmp(files[fsnowmelt] , string_novalue) != 0){
			snow->MELTED=new_doublevector(par->total_pixel);
			initialize_doublevector(snow->MELTED,0.);
			snow->melted=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fsnowsubl] , string_novalue) != 0){
			snow->SUBL=new_doublevector(par->total_pixel);
			initialize_doublevector(snow->SUBL,0.);
			snow->subl=new_doublevector(par->total_pixel);
		}
		if(strcmp(files[fsndur] , string_novalue) != 0){
			snow->t_snow=new_doublevector(par->total_pixel);
			initialize_doublevector(snow->t_snow,0.);
			snow->yes=new_shortvector(par->total_pixel);
		}
		if(strcmp(files[fHN] , string_novalue) != 0){//TODO (mattiu)
			snow->HNcum=new_doublevector(par->total_pixel);
			initialize_doublevector(snow->HNcum,0.);
			snow->yes=new_shortvector(par->total_pixel);
		}
	}
		
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			
			if( (long)land->LC->co[r][c]!=number_novalue){
				
				//Adjusting snow init depth in case of steep slope (contribution by Stephan Gruber)
				if (par->snow_curv > 0 && top->slope->co[r][c] > par->snow_smin){
					if (top->slope->co[r][c] <= par->snow_smax){
						k_snowred = ( exp(-pow(top->slope->co[r][c] - par->snow_smin, 2.)/par->snow_curv) -
									 exp(-pow(par->snow_smax, 2.)/par->snow_curv) );
					}else{
						k_snowred = 0.0;
					}
					snow->S->Dzl->co[1][r][c] *= k_snowred;
					snow->S->w_ice->co[1][r][c] *= k_snowred;
				}
				
				D = snow->S->Dzl->co[1][r][c];
				SWE = snow->S->w_ice->co[1][r][c];
								
				if(D<0 || SWE<0){
					f = fopen(FailedRunFile, "w");
					fprintf(f, "Error: negative initial snow depth %e or snow water equivalent %e\n",D,SWE);
					fclose(f);
					t_error("Fatal Error! Geotop is closed. See failing report (7).");		
					
				}else if(D<1.E-5 && SWE>1.E-5){
					f = fopen(FailedRunFile, "w");
					fprintf(f, "Error: Initial snow water equivalent %e > 0 and initial snow depth %e\n",SWE,D);
					fclose(f);
					t_error("Fatal Error! Geotop is closed. See failing report (8).");		
					
				}else if(D>1.E-5 && SWE<1.E-5){
					f = fopen(FailedRunFile, "w");
					fprintf(f, "Error: Initial snow depth %e > 0 and initial snow water equivalent %e\n",D,SWE);
					fclose(f);
					t_error("Fatal Error! Geotop is closed. See failing report (9).");		
					
				}else if(D>1.E-5 || SWE>1.E-5){
															
					snow->age->co[top->j_cont[r][c]]*=86400.0;	//now in [s]
							
					if (SWE <= par->max_weq_snow * par->max_snow_layers ) {

						i = floor( SWE/par->max_weq_snow );
						
						if (i>0) {
							
							for (n=1; n<=i; n++) {
								snow->S->w_ice->co[n][r][c] = par->max_weq_snow;
							}
							
							if (SWE - i * par->max_weq_snow > 0.1 * par->max_weq_snow) {
								snow->S->w_ice->co[i+1][r][c] = SWE - i * par->max_weq_snow;
								snow->S->lnum->co[r][c] = i+1;
							}else {
								snow->S->w_ice->co[i][r][c] += (SWE - i * par->max_weq_snow);
								snow->S->lnum->co[r][c] = i;
							}
						
						}else {
							
							snow->S->w_ice->co[1][r][c] = SWE;
							snow->S->lnum->co[r][c] = 1;
						
						}
						
					}else {
												
						snow->S->lnum->co[r][c] = par->max_snow_layers;

						for (n=1; n<=par->max_snow_layers; n++) {
							
							a = 0;
							
							for (i=1; i<=par->inf_snow_layers->nh; i++) {
								if (n == abs(par->inf_snow_layers->co[i])) a = 1;
							}
							
							if (a == 0) {
								snow->S->w_ice->co[n][r][c] = par->max_weq_snow;
							}else {
								snow->S->w_ice->co[n][r][c] = ( SWE - par->max_weq_snow * ( par->max_snow_layers - par->inf_snow_layers->nh ) ) / par->inf_snow_layers->nh;
							}
						}
					}
					
					for (n=1; n<=snow->S->lnum->co[r][c]; n++) {
						snow->S->Dzl->co[n][r][c] = D * (snow->S->w_ice->co[n][r][c] / SWE);
						snow->S->T->co[n][r][c] = IT->Tsnow0;
					}
					
				}
				
				non_dimensionalize_snowage(&(snow->age->co[top->j_cont[r][c]]), IT->Tsnow0);		
				
				if (par->point_sim == 1) {
					maxSWE = par->maxSWE->co[r][c];
				}else {
					maxSWE = 1.E10;
				}
								
				f = fopen(logfile, "a");
				snow_layer_combination(par->alpha_snow, r, c, snow->S, -0.1, par->inf_snow_layers, par->max_weq_snow, maxSWE, f);
				fclose(f);

			}
		}
	}

	if(par->delay_day_recover > 0){
		initialize_shortmatrix(snow->S->type, 2);
		assign_recovered_map_long(old, par->recover, files[rns], snow->S->lnum, par, land->LC);
		assign_recovered_map_vector(old, par->recover, files[rsnag], snow->age, top->rc_cont, par, land->LC);
		assign_recovered_tensor(old, par->recover, files[rDzs], snow->S->Dzl, par, land->LC);
		assign_recovered_tensor(old, par->recover, files[rwls], snow->S->w_liq, par, land->LC);
		assign_recovered_tensor(old, par->recover, files[rwis], snow->S->w_ice, par, land->LC);
		assign_recovered_tensor(old, par->recover, files[rTs], snow->S->T, par, land->LC);				
	}
	
	
				
	/****************************************************************************************************/
	/*! Initialization of the struct "glac" (of the type GLACIER):*/
	/***************************************************************************************************/ 
	/*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
	if( par->point_sim!=1 && strcmp(files[fgl0] , string_novalue) != 0 && par->max_glac_layers==0){
		printf("Warning: Glacier map present, but glacier represented with 0 layers\n");
		fprintf(flog,"Warning: Glacier map present, but glacier represented with 0 layers\n");
	}
	
	if(par->max_glac_layers==0 && IT->Dglac0>0){
		printf("\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
		fprintf(flog,"\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
	}
	
	//If the max number of glacier layers is greater than 1, the matrices (or tensors) lnum, Dzl. w_liq, w_ice, T and print matrices are defined, according to the respective flags
	if(par->max_glac_layers>0){
		
		if( par->point_sim!=1 && strcmp(files[fgl0] , string_novalue) != 0 ){
			printf("Glacier initial condition from file %s\n",files[fgl0]+1);
			fprintf(flog,"Glacier initial condition from file %s\n",files[fgl0]+1);
			M=read_map(2, files[fgl0], land->LC, UV, (double)number_novalue);
		}else{
			M=copydoublematrix_const(IT->Dglac0, land->LC, (double)number_novalue);
		}		
		
		glac->G=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
		allocate_and_initialize_statevar_3D(glac->G, (double)number_novalue, par->max_glac_layers, Nr, Nc);
		
		if(par->output_glac_bin == 1){
			if(strcmp(files[fglacmelt] , string_novalue) != 0){
				glac->MELTED=new_doublevector(par->total_pixel);
				initialize_doublevector(glac->MELTED,0.);
				glac->melted=new_doublevector(par->total_pixel);
			}
			if(strcmp(files[fglacsubl] , string_novalue) != 0){
				glac->SUBL=new_doublevector(par->total_pixel);
				initialize_doublevector(glac->SUBL, 0.);
				glac->subl=new_doublevector(par->total_pixel);
			}
		}
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if( (long)land->LC->co[r][c]!=number_novalue){
					
					if(M->co[r][c]<0){
						f = fopen(FailedRunFile, "w");
						fprintf(f, "Error: negative glacier data\n");
						fclose(f);
						t_error("Fatal Error! Geotop is closed. See failing report (10).");		
						
					}else if(M->co[r][c]>1.E-5){
														
						if (IT->rhoglac0 * M->co[r][c] / rho_w < par->max_weq_glac * par->max_glac_layers ) {
								
							n = 0;
							z = 0.;
								
							do{
									
								n++;
									
								if (IT->rhoglac0 * M->co[r][c] / rho_w < par->max_weq_glac * n) {
									glac->G->w_ice->co[n][r][c] = par->max_weq_glac;
								}else {
									glac->G->w_ice->co[n][r][c] = IT->rhoglac0 * M->co[r][c] / 1000. - z;
								}
									
								glac->G->Dzl->co[n][r][c] = rho_w * glac->G->w_ice->co[n][r][c] / IT->rhoglac0;
								glac->G->T->co[n][r][c] = IT->Tglac0;
									
								z += glac->G->w_ice->co[n][r][c];
									
							}while (fabs(z - IT->rhoglac0 * M->co[r][c] / rho_w) < 1.E-6);
								
							glac->G->lnum->co[r][c] = n;
								
						}else {
								
							glac->G->lnum->co[r][c] = par->max_glac_layers;
								
							for (n=1; n<=par->max_glac_layers; n++) {
									
								a = 0;
								
								for (i=1; i<=par->inf_glac_layers->nh; i++) {
									if (n == i) a = 1;
								}
									
								if (a == 0) {
									glac->G->w_ice->co[n][r][c] = par->max_weq_glac;
								}else {
									glac->G->w_ice->co[n][r][c] = ( IT->rhoglac0 * M->co[r][c] / rho_w - par->max_weq_glac * ( par->max_glac_layers - par->inf_glac_layers->nh ) ) / par->inf_glac_layers->nh;
								}
									
								glac->G->Dzl->co[n][r][c] = rho_w * glac->G->w_ice->co[n][r][c] / IT->rhoglac0;
								glac->G->T->co[n][r][c] = IT->Tglac0;
									
							}
						}
					}
										
					f = fopen(logfile, "a");
					snow_layer_combination(par->alpha_snow, r, c, glac->G, -0.1, par->inf_glac_layers, par->max_weq_glac, 1.E10, f);							
					fclose(f);
									
				}
			}
		}
		
		free_doublematrix(M);
		
		if(par->delay_day_recover > 0){
			assign_recovered_map_long(old, par->recover, files[rni], glac->G->lnum, par, land->LC);
			assign_recovered_tensor(old, par->recover, files[rDzi], glac->G->Dzl, par, land->LC);
			assign_recovered_tensor(old, par->recover, files[rwli], glac->G->w_liq, par, land->LC);
			assign_recovered_tensor(old, par->recover, files[rwii], glac->G->w_ice, par, land->LC);
			assign_recovered_tensor(old, par->recover, files[rTi], glac->G->T, par, land->LC);
		}
	}																				
			
	
	
	//***************************************************************************************************
	// Filling up of the struct "met" (of the type METEO):
		
	met->Tgrid=new_doublematrix(Nr,Nc);
	initialize_doublematrix(met->Tgrid, par->Tair_default);
	met->Pgrid=new_doublematrix(Nr,Nc);
	initialize_doublematrix(met->Pgrid, Pa0);
	met->RHgrid=new_doublematrix(Nr,Nc);
	initialize_doublematrix(met->RHgrid, par->RH_default);
	met->Vgrid=new_doublematrix(Nr,Nc);
	initialize_doublematrix(met->Vgrid, par->V_default);
	met->Vdir=new_doublematrix(Nr,Nc);
	initialize_doublematrix(met->Vdir, par->Vdir_default);
		
	if (par->output_meteo_bin == 1){
		if(strcmp(files[fTa] , string_novalue) != 0){
			met->Tamean=new_doublevector(par->total_pixel);
			initialize_doublevector(met->Tamean, 0.);
		}
		if(strcmp(files[fwspd] , string_novalue) != 0){
			met->Vspdmean=new_doublevector(par->total_pixel);
			initialize_doublevector(met->Vspdmean, 0.);
		}
		if(strcmp(files[fwdir] , string_novalue) != 0){
			met->Vdirmean=new_doublevector(par->total_pixel);
			initialize_doublevector(met->Vdirmean, 0.);
		}
		if(strcmp(files[frh] , string_novalue) != 0){
			met->RHmean=new_doublevector(par->total_pixel);
			initialize_doublevector(met->RHmean, 0.);
		}
	}	
	
	
	//plot output
	if(times->JD_plots->nh > 1){
		if(strcmp(files[pTa] , string_novalue) != 0){
			met->Taplot=new_doublevector(par->total_pixel);
			initialize_doublevector(met->Taplot, 0.);
		}
		if(strcmp(files[pRH] , string_novalue) != 0){
			met->RHplot=new_doublevector(par->total_pixel);
			initialize_doublevector(met->RHplot, 0.);
		}
		if(strcmp(files[pVspd] , string_novalue) != 0 || strcmp(files[pVdir] , string_novalue) != 0){
			met->Vxplot=new_doublevector(par->total_pixel);
			met->Vyplot=new_doublevector(par->total_pixel);
			initialize_doublevector(met->Vxplot, 0.);
			initialize_doublevector(met->Vyplot, 0.);
		}		
	}
	
	/****************************************************************************************************/			
	//SpinUp variables
	
	if(par->recover > 0) write_suffix(rec, par->recover, 4);
	if(par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);
	
	if(par->Tzrun == 1){
		
		sl->Tzrun = new_doublematrix(par->rc->nrh, Nl);
		initialize_doublematrix(sl->Tzrun, 0.);
		
		if (par->recover>0) {
			temp = join_strings(files[fTrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fTrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fTrun], textfile);
		}
		
		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point");
		for (l=1; l<=Nl; l++) {
			fprintf(f, ",l[%ld]",l);
		}
		fprintf(f, "\n");
		free(name);
		fclose(f);
	}
	
	if(par->wzrun == 1){
		sl->wzrun = new_doublematrix(par->rc->nrh, Nl);
		initialize_doublematrix(sl->wzrun, 0.);
		
		if (par->recover>0) {
			temp = join_strings(files[fwrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fwrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fwrun], textfile);
		}
		
		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point");
		for (l=1; l<=Nl; l++) {
			fprintf(f, ",l[%ld]",l);
		}
		fprintf(f, "\n");
		free(name);
		fclose(f);
	}
	
	if(par->Tzmaxrun == 1){
		sl->Tzmaxrun = new_doublematrix(par->rc->nrh, Nl);
		initialize_doublematrix(sl->Tzmaxrun, -1.E99);
		
		if (par->recover>0) {
			temp = join_strings(files[fTmaxrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fTmaxrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fTmaxrun], textfile);
		}
		
		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point");
		for (l=1; l<=Nl; l++) {
			fprintf(f, ",l[%ld]",l);
		}
		fprintf(f, "\n");
		free(name);
		fclose(f);
	}
	
	if(par->wzmaxrun == 1){
		sl->wzmaxrun = new_doublematrix(par->rc->nrh, Nl);
		initialize_doublematrix(sl->wzmaxrun, -1.E99);
		
		if (par->recover>0) {
			temp = join_strings(files[fwmaxrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fwmaxrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fwmaxrun], textfile);
		}

		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point");
		for (l=1; l<=Nl; l++) {
			fprintf(f, ",l[%ld]",l);
		}
		fprintf(f, "\n");
		free(name);
		fclose(f);
	}
	
	if(par->Tzminrun == 1){
		sl->Tzminrun = new_doublematrix(par->rc->nrh, Nl);
		initialize_doublematrix(sl->Tzminrun, 1.E99);
		
		if (par->recover>0) {
			temp = join_strings(files[fTminrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fTminrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fTminrun], textfile);
		}
		
		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point");
		for (l=1; l<=Nl; l++) {
			fprintf(f, ",l[%ld]",l);
		}
		fprintf(f, "\n");
		free(name);
		fclose(f);
	}
	
	if(par->wzminrun == 1){
		sl->wzminrun = new_doublematrix(par->rc->nrh, Nl);
		initialize_doublematrix(sl->wzminrun, 1.E99);
		
		if (par->recover>0) {
			temp = join_strings(files[fwminrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fwminrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fwminrun], textfile);
		}
		
		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point");
		for (l=1; l<=Nl; l++) {
			fprintf(f, ",l[%ld]",l);
		}
		fprintf(f, "\n");
		free(name);
		fclose(f);	
	}
	
	if(par->dUzrun == 1){
		sl->dUzrun = new_doublematrix(par->rc->nrh, Nl);
		initialize_doublematrix(sl->dUzrun, 0.);
		
		if (par->recover>0) {
			temp = join_strings(files[fdUrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fdUrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fdUrun], textfile);
		}

		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point");
		for (l=1; l<=Nl; l++) {
			fprintf(f, ",l[%ld]",l);
		}
		fprintf(f, "\n");
		free(name);
		fclose(f);
	}		
	
	if(par->SWErun == 1){
		sl->SWErun = new_doublematrix(par->rc->nrh, 3);//mean,max,min
		for (l=1; l<=par->rc->nrh; l++) {
			sl->SWErun->co[l][1] = 0.;
			sl->SWErun->co[l][2] = -1.E99;
			sl->SWErun->co[l][3] = 1.E99;				
		}
		
		if (par->recover>0) {
			temp = join_strings(files[fSWErun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fSWErun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fSWErun], textfile);
		}
		
		f = fopen(name, "w");
		fprintf(f, "Period,Run,Point,SWEmean[mm],SWEmax[mm],SWEmin[mm]\n");
		free(name);
		fclose(f);
	}				

	if(par->delay_day_recover > 0){
		if(par->Tzrun == 1) recover_run_averages(old, sl->Tzrun, files[rTrun], land->LC, top->rc_cont, par, Nl);
		if(par->wzrun == 1) recover_run_averages(old, sl->wzrun, files[rwrun], land->LC, top->rc_cont, par, Nl);
		if(par->Tzminrun == 1) recover_run_averages(old, sl->Tzminrun, files[rTminrun], land->LC, top->rc_cont, par, Nl);
		if(par->wzminrun == 1) recover_run_averages(old, sl->wzminrun, files[rwminrun], land->LC, top->rc_cont, par, Nl);
		if(par->Tzmaxrun == 1) recover_run_averages(old, sl->Tzmaxrun, files[rTmaxrun], land->LC, top->rc_cont, par, Nl);
		if(par->wzmaxrun == 1) recover_run_averages(old, sl->wzmaxrun, files[rwmaxrun], land->LC, top->rc_cont, par, Nl);
		if(par->dUzrun == 1) recover_run_averages(old, sl->dUzrun, files[rdUrun], land->LC, top->rc_cont, par, Nl);
		if(par->SWErun == 1) recover_run_averages(old, sl->SWErun, files[rSWErun], land->LC, top->rc_cont, par, 3);
	}

	/****************************************************************************************************/			
	
	//Free the struct allocated in this subroutine:
	//free_doublematrix(par->chkpt);
	
	free_doublevector(IT->init_water_table_depth);
	free_doublematrix(IT->bed);
	free_doubletensor(IT->pa_bed);
	
	for (i=0; i<nmet; i++) {
		free(IT->met_col_names[i]);
	}
	free(IT->met_col_names);
	
	for (i=0; i<nsoilprop; i++) {
		free(IT->soil_col_names[i]);
	}	
	free(IT->soil_col_names);
	
	for (i=0; i<2; i++) {
		free(IT->horizon_col_names[i]);
	}
	free(IT->horizon_col_names);
	
	for (i=0; i<ptTOT; i++) {
		free(IT->point_col_names[i]);
	}
	free(IT->point_col_names);

	for (i=0; i<nlstot; i++) {
		free(IT->lapserates_col_names[i]);
	}
	free(IT->lapserates_col_names);

	for (i=0; i<8; i++) {
		free(IT->meteostations_col_names[i]);
	}	
	free(IT->meteostations_col_names);
	
	free(IT);
	
	n = Fminlong(par->Nl_spinup->co[i_sim0],Nl);

	if (par->point_sim != 1) {
		cont_nonzero_values_matrix2(&i, &j, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel, par->total_channel, n);
		top->Li = new_longvector(i);
		top->Lp = new_longvector(j);
		wat->Lx = new_doublevector(i);	
		cont_nonzero_values_matrix3(top->Lp, top->Li, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel, par->total_channel, n);
	}else {
		i = n;
		j = n+1;
		top->Li = new_longvector(i);
		top->Lp = new_longvector(j);
		wat->Lx = new_doublevector(i);	
		for (l=1; l<=n; l++) {
			top->Li->co[l] = l+1;
			top->Lp->co[l] = l;
		}
		top->Lp->co[j] = i;
	}

	wat->H0 = new_doublevector(j);
	wat->H1 = new_doublevector(j);
	wat->dH = new_doublevector(j);
	wat->B = new_doublevector(j);
	wat->f = new_doublevector(j);
	wat->df = new_doublevector(j);
	
	wat->Kbottom = new_doublematrix(Nr, Nc);
	initialize_doublematrix(wat->Kbottom, 0.);
	
	wat->Klat = new_doublematrix(top->BC_DepthFreeSurface->nh, Nl);
	initialize_doublematrix(wat->Klat, 0.);

	fclose(flog);
	
			
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void read_inputmaps(TOPO *top, LAND *land, SOIL *sl, PAR *par, INIT_TOOLS *IT, FILE *flog){

	long r, c, i, cont;					
	DOUBLEMATRIX *M;
	SHORTMATRIX *curv;
	short flag;
	char *temp;
	double min, max;
	FILE *f;
	
	//reading TOPOGRAPHY
	flag = file_exists(fdem, flog);	
	if(flag == 1){
		M=new_doublematrix(1,1);
		top->Z0=read_map(0, files[fdem], M, UV, (double)number_novalue); //topography
		free_doublematrix(M);
		write_map(files[fdem], 0, par->format_out, top->Z0, UV, number_novalue);
					
		//filtering
		M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
		multipass_topofilter(par->lowpass, top->Z0, M, (double)number_novalue, 1);
		copy_doublematrix(M, top->Z0);
		free_doublematrix(M);
		//write_map(files[fdem], 0, par->format_out, top->Z0, UV, number_novalue);
		
		//calculate East and North
		top->East = new_doublematrix(top->Z0->nrh, top->Z0->nch);
		top->North = new_doublematrix(top->Z0->nrh, top->Z0->nch);
		for (r=1; r<=top->Z0->nrh; r++) {
			for (c=1; c<=top->Z0->nch; c++) {
				top->East->co[r][c] = UV->U->co[4] + (c-0.5)*UV->U->co[2];
				top->North->co[r][c] = UV->U->co[3] + (top->Z0->nrh-(r-0.5))*UV->U->co[1];
			}
		}
						
	}else{
		
		f = fopen(FailedRunFile, "w");
		fprintf(f, "Error: It is impossible to proceed without giving the digital elevation model\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report (11).");		
		
	}

	//reading LAND COVER TYPE
	flag = file_exists(flu, flog);
	if(flag == 1){
		land->LC=read_map(1, files[flu], top->Z0, UV, (double)number_novalue);
		//Check borders
		for(r=1;r<=land->LC->nrh;r++){
			land->LC->co[r][1]=(double)number_novalue;
			land->LC->co[r][land->LC->nch]=(double)number_novalue;
		}
		for(c=1;c<=land->LC->nch;c++){
			land->LC->co[1][c]=(double)number_novalue;
			land->LC->co[land->LC->nrh][c]=(double)number_novalue;
		}		
		for(r=1;r<=land->LC->nrh;r++){
			for(c=1;c<=land->LC->nch;c++){
				if ((long)land->LC->co[r][c] != number_novalue) {
					if ((long)land->LC->co[r][c] < 1 || (long)land->LC->co[r][c] > par->n_landuses){
						f = fopen(FailedRunFile, "w");
						fprintf(f, "Error: It is not possible to assign Value < 1 or > n_landuses to the land cover type\n");
						fclose(f);
						t_error("Fatal Error! Geotop is closed. See failing report (12).");		
					}
				}
			}
		}
		
		//Land use is the official mask
		for(r=1;r<=land->LC->nrh;r++){
			for(c=1;c<=land->LC->nch;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					if((long)top->Z0->co[r][c]==number_novalue){
						printf("ERROR Land use mask include DTM novalue pixels");
						printf("\nr:%ld c:%ld Z:%f landuse:%f\n",r,c,top->Z0->co[r][c],land->LC->co[r][c]);
						land->LC->co[r][c]=(double)number_novalue;
						printf("LANDUSE set at novalue where DTM is not available\n"); 
					}
				}
			}
		}

	}else{

		//Write land->LC (land cover)
		printf("Land cover type assumed to be always 1\n");
		land->LC=copydoublematrix_const(1.0, top->Z0, (double)number_novalue);
		
		for(r=1;r<=land->LC->nrh;r++){
			land->LC->co[r][1]=(double)number_novalue;
			land->LC->co[r][land->LC->nch]=(double)number_novalue;
		}
		for(c=1;c<=land->LC->nch;c++){
			land->LC->co[1][c]=(double)number_novalue;
			land->LC->co[land->LC->nrh][c]=(double)number_novalue;
		}		
	}
	if(flag >= 0) write_map(files[flu], 1, par->format_out, land->LC, UV, number_novalue);
	
	if(par->state_pixel == 1){
		par->rc=new_longmatrix(par->chkpt->nrh,2);
		par->IDpoint=new_longvector(par->chkpt->nrh);
		for(i=1;i<=par->chkpt->nrh;i++){	

			par->rc->co[i][1]=row(par->chkpt->co[i][ptY], top->Z0->nrh, UV, number_novalue);
			par->rc->co[i][2]=col(par->chkpt->co[i][ptX], top->Z0->nch, UV, number_novalue);
						
			if (par->rc->co[i][1] == number_novalue || par->rc->co[i][2] == number_novalue) {
				printf("Point #%4ld is out of the domain",i);
				
				fprintf(flog, "Point #%4ld is out of the domain",i);
				fclose(flog);

				f = fopen(FailedRunFile, "w");
				fprintf(f, "Point #%4ld is out of the domain",i);				
				fclose(f);
				t_error("Fatal Error! Geotop is closed. See failing report.");						
			}
			
			if((long)land->LC->co[par->rc->co[i][1]][par->rc->co[i][2]]==number_novalue){
				printf("Point #%4ld corresponds to NOVALUE pixel",i);

				fprintf(flog, "Point #%4ld corresponds to NOVALUE pixel",i);
				fclose(flog);
				
				f = fopen(FailedRunFile, "w");
				fprintf(f, "Point #%4ld corresponds to NOVALUE pixel",i);
				fclose(f);
				t_error("Fatal Error! Geotop is closed. See failing report.");										
			}
			
			if ((long)par->chkpt->co[i][ptID]!=number_novalue) {
				par->IDpoint->co[i]=(long)par->chkpt->co[i][ptID];
			}else {
				par->IDpoint->co[i]=i;
			}

		}
	}
	
	
/****************************************************************************************************/
	
	//reading SKY VIEW FACTOR
	flag = file_exists(fsky, flog);
	if(flag == 1){
		top->sky=read_map(2, files[fsky], land->LC, UV, (double)number_novalue);
	}else{/*The sky view factor file "top->sky" must be calculated*/
		top->sky = new_doublematrix(top->Z0->nrh,top->Z0->nch);
		if (par->sky == 0) {
			initialize_doublematrix(top->sky, 1.);
		}else {
			curv = new_shortmatrix(top->Z0->nrh,top->Z0->nch);
			nablaquadro_mask(top->Z0, curv, UV->U, UV->V);
			sky_view_factor(top->sky, 36, UV, top->Z0, curv, number_novalue);
			free_shortmatrix(curv);
		}
	}
	if(flag >= 0) write_map(files[fsky], 0, par->format_out, top->sky, UV, number_novalue);

/****************************************************************************************************/
	
	//reading DELAY
	flag = file_exists(fdelay, flog);
	if(flag == 1){
		land->delay = read_map(2, files[fdelay], land->LC, UV, (double)number_novalue);
	}else{
		land->delay = new_doublematrix(top->Z0->nrh,top->Z0->nch);
		initialize_doublematrix(land->delay, 0);
	}
	if(flag >= 0) write_map(files[fdelay], 0, par->format_out, land->delay, UV, number_novalue);
	
/****************************************************************************************************/
	//reading SOIL MAP
	flag = file_exists(fsoil, flog);
	if(flag == 1){
		M=read_map(2, files[fsoil], land->LC, UV, (double)number_novalue);
		sl->type=copylong_doublematrix(M);
		for(r=1;r<=land->LC->nrh;r++){
			for(c=1;c<=land->LC->nch;c++){
				if ((long)land->LC->co[r][c] != number_novalue) {
					if (sl->type->co[r][c] < 1 || sl->type->co[r][c] > par->nsoiltypes){
						f = fopen(FailedRunFile, "w");
						fprintf(f, "Error: It is not possible to assign Value < 1 or > nsoiltypes to the soil type map");
						fclose(f);
						t_error("Fatal Error! GEOtop is closed. See failing report (13).");
					}	
				}
			}
		}
	}else{
		M=copydoublematrix_const(par->soil_type_land_default, land->LC, (double)number_novalue);
		sl->type=copylong_doublematrix(M);
	}
	if(flag >= 0) write_map(files[fsoil], 1, par->format_out, M, UV, number_novalue);
	free_doublematrix(M);
	

/****************************************************************************************************/
	//SLOPE	
	top->dzdE = new_doublematrix(land->LC->nrh, land->LC->nch);
	top->dzdN = new_doublematrix(land->LC->nrh, land->LC->nch);
	find_slope(UV->U->co[1], UV->U->co[2], top->Z0, top->dzdE, top->dzdN, (double)number_novalue);
	
	flag = file_exists(fslp, flog);
	if(flag == 1){
		top->slope=read_map(2, files[fslp], land->LC, UV, (double)number_novalue);		//reads in degrees
	}else{
		top->slope=find_max_slope(top->Z0, top->dzdE, top->dzdN, (double)number_novalue);
	}
	if(flag >= 0) write_map(files[fslp], 0, par->format_out, top->slope, UV, number_novalue);
	
	find_min_max(top->slope, (double)number_novalue, &max, &min);
	printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min,tan(max*Pi/180.),max);
	fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min,tan(max*Pi/180.),max);
	
/****************************************************************************************************/
	//ASPECT
	flag = file_exists(fasp, flog);
	if(flag == 1){
		top->aspect=read_map(2, files[fasp], land->LC, UV, (double)number_novalue);
	}else{
		top->aspect=find_aspect(top->Z0, top->dzdE, top->dzdN, (double)number_novalue);
	}
	if(flag >= 0) write_map(files[fasp], 0, par->format_out, top->aspect, UV, number_novalue);


/****************************************************************************************************/
	//curvature
	top->curvature1=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	top->curvature2=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	top->curvature3=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	top->curvature4=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	
	//filtering
	M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
	multipass_topofilter(par->lowpass_curvatures, top->Z0, M, (double)number_novalue, 1);
	curvature(UV->U->co[1], UV->U->co[2], M, top->curvature1, top->curvature2, top->curvature3, top->curvature4, (double)number_novalue);
	free_doublematrix(M);

	if(strcmp(files[fcurv] , string_novalue) != 0){

		temp = join_strings(files[fcurv], "N-S");
		write_map(temp, 0, par->format_out, top->curvature1, UV, number_novalue);
		free(temp);
		
		temp = join_strings(files[fcurv], "W-E");
		write_map(temp, 0, par->format_out, top->curvature2, UV, number_novalue);
		free(temp);
		
		temp = join_strings(files[fcurv], "NW-SE");
		write_map(temp, 0, par->format_out, top->curvature3, UV, number_novalue);
		free(temp);		
		
		temp = join_strings(files[fcurv], "NE-SW");
		write_map(temp, 0, par->format_out, top->curvature4, UV, number_novalue);
		free(temp);							
	
	}
	
	find_min_max(top->curvature1, (double)number_novalue, &max, &min);
	printf("Curvature N-S Min:%f  Max:%f \n",min,max);
	fprintf(flog,"Curvature N-S Min:%f  Max:%f \n",min,max);
	
	find_min_max(top->curvature2, (double)number_novalue, &max, &min);
	printf("Curvature W-E Min:%f  Max:%f \n",min,max);
	fprintf(flog,"Curvature W-E Min:%f  Max:%f \n",min,max);

	find_min_max(top->curvature3, (double)number_novalue, &max, &min);
	printf("Curvature NW-SE Min:%f  Max:%f \n",min,max);
	fprintf(flog,"Curvature NW-SE Min:%f  Max:%f \n",min,max);

	find_min_max(top->curvature4, (double)number_novalue, &max, &min);
	printf("Curvature NE-SW Min:%f  Max:%f \n",min,max);
	fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
			
/****************************************************************************************************/	
//Channel network (in top->pixel_type)
	
	//pixel type = 0 land pixel (if it is on the border, the border is impermeable, water is free only on the surface)
	//pixel type = 1 or 2 land pixel (it it is on the border, the border is permeable above an user-defined elevation in the saturated part, weir-wise)
	//pixel type = 10 channel pixel (if it is on the border, the border is impermeable, water is free only on the surface)
	//pixel type = -1 land pixel where an incoming discharge from outside is considered (as rain)
	
	flag = file_exists(fnet, flog);
	if(flag == 1){
		M=read_map(2, files[fnet], land->LC, UV, (double)number_novalue);
		top->pixel_type=copyshort_doublematrix(M);

		cont = 0;
		for(r=1;r<=top->Z0->nrh;r++){
			for(c=1;c<=top->Z0->nch;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					if(top->pixel_type->co[r][c]!=0 && top->pixel_type->co[r][c]!=1 && top->pixel_type->co[r][c]!=2 && top->pixel_type->co[r][c]!=10){
						f = fopen(FailedRunFile, "w");
						fprintf(f, "Error: Only the following values are admitted in the network map: 0, 1, 10\n");
						fclose(f);
						t_error("Fatal Error! GEOtop is closed. See failing report (14).");
					}
					if(top->pixel_type->co[r][c]==10) cont++;
				}
			}
		}
	
		printf("Channel networks has %ld pixels set to channel\n",cont);
		
		if(flag >= 0) write_map(files[fnet], 1, par->format_out, M, UV, number_novalue);
		free_doublematrix(M);
	
	}else{
		
		top->pixel_type=new_shortmatrix(land->LC->nrh, land->LC->nch);
		initialize_shortmatrix(top->pixel_type, 0);
				
	}
	
	/****************************************************************************************************/	
	
	//border
	top->is_on_border=new_shortmatrix(land->LC->nrh, land->LC->nch);
	for(r=1;r<=land->LC->nrh;r++){
		for(c=1;c<=land->LC->nch;c++){
			if ( (long)land->LC->co[r][c]!=number_novalue){
				top->is_on_border->co[r][c] = is_boundary(r, c, land->LC, number_novalue);
			}else{
				top->is_on_border->co[r][c] = -1;
			}
		}
	}
	
	//count the pixels having pixel_type = 1, 2 or -1
	cont = 0;
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			if(top->is_on_border->co[r][c]==1){
				if (top->pixel_type->co[r][c] == -1 || top->pixel_type->co[r][c] == 1 || top->pixel_type->co[r][c] == 2) cont ++;
			}
		}
	}
	
	top->BC_counter = new_longmatrix(top->Z0->nrh, top->Z0->nch);
	initialize_longmatrix(top->BC_counter, 0);
	
	if (cont > 0) {
		top->BC_DepthFreeSurface = new_doublevector(cont);
		cont = 0;
		for(r=1;r<=top->Z0->nrh;r++){
			for(c=1;c<=top->Z0->nch;c++){
				if(top->is_on_border->co[r][c]==1){
					if (top->pixel_type->co[r][c] == -1 || top->pixel_type->co[r][c] == 1 || top->pixel_type->co[r][c] == 2){
						cont ++;
						top->BC_counter->co[r][c] = cont;
						top->BC_DepthFreeSurface->co[cont] = par->DepthFreeSurface; //[mm]
					}
				}
			}
		}
	}else {
		top->BC_DepthFreeSurface = new_doublevector(1);	
		initialize_doublevector(top->BC_DepthFreeSurface, (double)number_novalue);
	}
	
	//bedrock
	flag = file_exists(fbed, flog);
	if(flag == 1){
		IT->bed=read_map(2, files[fbed], land->LC, UV, (double)number_novalue);
	}else{
		IT->bed=new_doublematrix(top->Z0->nrh, top->Z0->nch);
		initialize_doublematrix(IT->bed, 1.E99);
	}
	if(flag>=0) write_map(files[fbed], 0, par->format_out, top->aspect, UV, number_novalue);


}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void read_optionsfile_point(PAR *par, TOPO *top, LAND *land, SOIL *sl, TIMES *times, INIT_TOOLS *IT, FILE *flog){

	long i, r, c, num_lines;
	DOUBLEMATRIX *Q=NULL, *P=NULL, *R=NULL, *S=NULL, *T=NULL, *Z=NULL, *LU=NULL; // ec 2012 08 22
	SHORTMATRIX *curv;
	short read_dem, read_lu, read_soil, read_sl, read_as, read_sk, read_bed, read_curv, flag, coordinates;
	char *temp;
	double min, max;
	FILE *f;
	
	//4. CALCULATE TOPOGRAPHIC PROPERTIES
	//check if there are point coordinates
	coordinates = 1;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if ( (long)par->chkpt->co[i][ptX]==number_novalue || (long)par->chkpt->co[i][ptY]==number_novalue ) coordinates = 0;
	}
	/*if (coordinates == 0 && par->recover>0){
		printf("Warning: Not possible to recover the simulation because at least one point has no coordinates\n");
		printf("Starting from normal initial condition\n");
		fprintf(flog,"Warning: Not possible to recover the simulation because at least one point has no coordinates\n");
		fprintf(flog,"Starting from normal initial condition\n");
		par->recover = 0;		
	}*/
	
	//a. read dem	
	read_dem=0;
	//if(par->recover>0) read_dem=1;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptLC]==number_novalue || (long)par->chkpt->co[i][ptSY]==number_novalue || 
		   (long)par->chkpt->co[i][ptS]==number_novalue || (long)par->chkpt->co[i][ptA]==number_novalue ||
		   (long)par->chkpt->co[i][ptCNS]==number_novalue || (long)par->chkpt->co[i][ptCWE]==number_novalue ||
		   (long)par->chkpt->co[i][ptCNwSe]==number_novalue || (long)par->chkpt->co[i][ptCNeSw]==number_novalue){
			read_dem=1;
		}
	}
	if(read_dem == 1 && coordinates == 0){
		printf("Warning: Not possible to read from dem because at least one point has no coordinates\n");
		fprintf(flog,"Warning: Not possible to read from dem because at least one point has no coordinates\n");
		read_dem = 0;
	}	
	if(read_dem==1){
		flag = file_exists(fdem, flog);	
		if(flag == 1){			
			printf("Warning: Dem file %s present\n",files[fdem]+1);
			fprintf(flog,"Warning: Dem file %s present\n",files[fdem]+1);
			
			Q=new_doublematrix(1,1);
			Z=read_map(0, files[fdem], Q, UV, (double)number_novalue); //topography
			free_doublematrix(Q);
			
			Q=new_doublematrix(Z->nrh,Z->nch);
			multipass_topofilter(par->lowpass, Z, Q, (double)number_novalue, 1);
			copy_doublematrix(Q, Z);
			free_doublematrix(Q);
			
		}else{
			
			read_dem=0;
			printf("Warning: Dem file not present\n");
			fprintf(flog,"Warning: Dem file not present\n");
			
			/*if(par->recover>0){
				printf("Warning: Not possible to recover the simulation because there is no dem\n");
				printf("Starting from normal initial condition\n");
				fprintf(flog,"Warning: Not possible to recover the simulation because there is no dem\n");
				fprintf(flog,"Starting from normal initial condition\n");
				par->recover = 0;
			}*/	
		}
	}
	
	
	if(read_dem==1){ 
		par->r_points=new_longvector(par->chkpt->nrh);
		par->c_points=new_longvector(par->chkpt->nrh);		
		for(i=1;i<=par->chkpt->nrh;i++){ 
			par->r_points->co[i]=row(par->chkpt->co[i][ptY], Z->nrh, UV, number_novalue);
			par->c_points->co[i]=col(par->chkpt->co[i][ptX], Z->nch, UV, number_novalue);
			if((long)par->chkpt->co[i][ptZ]==number_novalue) par->chkpt->co[i][ptZ]=Z->co[par->r_points->co[i]][par->c_points->co[i]];
			printf("ok");
		}
	}

	//b. read land use
	read_lu=0;
	//if(par->recover>0) read_lu=1;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptLC]==number_novalue) read_lu=1;
	}
	if(read_lu==1 && coordinates==0) read_lu=0;
	if(read_lu==1){
		flag = file_exists(flu, flog);	
		if(flag == 1){		
			if(read_dem==0){
				Q=new_doublematrix(1,1);
				LU=read_map(0, files[flu], Q, UV, (double)number_novalue);
				free_doublematrix(Q);
			}else{
				LU=read_map(1, files[flu], Z, UV, (double)number_novalue);
			}
			
		}else{
			printf("Warning: Landuse file not present, uniform cover considered\n");
			fprintf(flog,"Warning: Landuse file not present, uniform cover considered\n");
			if(read_dem==1){
				LU=copydoublematrix_const(1.0, Z, (double)number_novalue);
			}else{
				read_lu=0;
				/*if(par->recover>0){
					printf("Warning: Not possible to recover the simulation because there is no dem\n");
					printf("Starting from normal initial condition\n");
					fprintf(flog,"Warning: Not possible to recover the simulation because there is no dem\n");
					fprintf(flog,"Starting from normal initial condition\n");
					par->recover = 0;
				}*/						
			}
		}
	}
	
	if(read_lu==1){
		for(i=1;i<=par->chkpt->nrh;i++){ 
			if((long)par->chkpt->co[i][ptLC]==number_novalue){
				r=row(par->chkpt->co[i][ptY], LU->nrh, UV, number_novalue);
				c=col(par->chkpt->co[i][ptX], LU->nch, UV, number_novalue);
				par->chkpt->co[i][ptLC]=LU->co[r][c];
			}
		}
	}
	
	//c. read soil type
	read_soil=0;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptSY]==number_novalue) read_soil=1; 
	}
	if(read_soil==1 && coordinates==0) read_soil=0;
	if(read_soil==1){
		flag = file_exists(fsoil, flog);	
		if(flag == 1){		
			if(read_dem==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files[fsoil], Q, UV, (double)number_novalue);
				free_doublematrix(Q);
			}else{
				P=read_map(1, files[fsoil], Z, UV, (double)number_novalue);
			}
			
		}else{
			printf("Warning: Soiltype file not present\n");
			fprintf(flog,"Warning: Soiltype file not present\n");
			read_soil=0;
		}
	}
	if(read_soil==1){
		for(i=1;i<=par->chkpt->nrh;i++){ 
			if((long)par->chkpt->co[i][ptSY]==number_novalue){
				r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
				c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				par->chkpt->co[i][ptSY]=P->co[r][c];
			}
		}
		free_doublematrix(P);
	}
		
	//d. read slope		
	read_sl=0;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptS]==number_novalue) read_sl=1; 
	}
	if(read_sl==1 && coordinates==0) read_sl=0;
	if(read_sl==1){
		flag = file_exists(fslp, flog);	
		if(flag == 1){				
			if(read_dem==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files[fslp], Q, UV, (double)number_novalue);
				free_doublematrix(Q);
			}else{
				P=read_map(1, files[fslp], Z, UV, (double)number_novalue);
			}
			
		}else{
			if(read_dem==0){
				printf("Warning: Slopes file not present\n");
				fprintf(flog,"Warning: Slopes file not present\n");
				read_sl=0;
			}else{
				Q=new_doublematrix(Z->nrh,Z->nch);
				R=new_doublematrix(Z->nrh,Z->nch);
				find_slope(UV->U->co[1], UV->U->co[2], Z, Q, R, (double)number_novalue);
				P=find_max_slope(Z, Q, R, (double)number_novalue);
				free_doublematrix(Q);
				free_doublematrix(R);
				if(flag==0) write_map(files[fslp], 0, par->format_out, P, UV, number_novalue);
			}
		}
	}

	if(read_sl==1){
		find_min_max(P, (double)number_novalue, &max, &min);
		printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min,tan(max*Pi/180.),max);
		fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min,tan(max*Pi/180.),max);
		
		for(i=1;i<=par->chkpt->nrh;i++){ 

			if((long)par->chkpt->co[i][ptS]==number_novalue){
				r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
				c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				par->chkpt->co[i][ptS]=P->co[r][c];			
			}
		}
		free_doublematrix(P);
	}

	//e. read aspect		
	read_as=0;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptA]==number_novalue) read_as=1; 
	}
	if(read_as==1 && coordinates==0) read_as=0;
	if(read_as==1){
		flag = file_exists(fasp, flog);	
		if(flag == 1){		
			if(read_dem==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files[fasp], Q, UV, (double)number_novalue);
				free_doublematrix(Q);
			}else{
				P=read_map(1, files[fasp], Z, UV, (double)number_novalue);
			}
		}else{
			if(read_dem==0){
				printf("Warning: Aspect file not present\n");
				fprintf(flog,"Warning: Aspect file not present\n");
				read_as=0;
			}else{
				Q=new_doublematrix(Z->nrh,Z->nch);
				R=new_doublematrix(Z->nrh,Z->nch);
				find_slope(UV->U->co[1], UV->U->co[2], Z, Q, R, (double)number_novalue);
				P=find_aspect(Z, Q, R, (double)number_novalue);
				free_doublematrix(Q);
				free_doublematrix(R);
				if(flag==0) write_map(files[fasp], 0, par->format_out, P, UV, number_novalue);
			}
		}
	}

	if(read_as==1){
		for(i=1;i<=par->chkpt->nrh;i++){ 
			if((long)par->chkpt->co[i][ptA]==number_novalue){
				r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
				c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				par->chkpt->co[i][ptA]=P->co[r][c];		
			}
		}
		free_doublematrix(P);
	}
	
	//f. sky view factor file
	read_sk=0;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptSKY]==number_novalue) read_sk=1; 
	}
	if(read_sk==1 && coordinates==0) read_sk=0;
	if(read_sk==1){
		flag = file_exists(fsky, flog);	
		if(flag == 1){				
			if(read_dem==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files[fsky], Q, UV, (double)number_novalue);
				free_doublematrix(Q);
			}else{
				P=read_map(1, files[fsky], Z, UV, (double)number_novalue);
			}
		}else{
			if(read_dem==0){
				printf("Warning: Sky view factor file not present\n");
				fprintf(flog,"Warning: Sky view factor file not present\n");
				read_sk=0;
			}else{
				P=new_doublematrix(Z->nrh,Z->nch);
				curv=new_shortmatrix(Z->nrh,Z->nch);
				nablaquadro_mask(Z, curv, UV->U, UV->V);
				sky_view_factor(P, 36, UV, Z, curv, number_novalue);
				free_shortmatrix(curv);
				if(flag==0) write_map(files[fsky], 0, par->format_out, P, UV, number_novalue);				
			}
		}
	}
	
	if(read_sk==1){
		for(i=1;i<=par->chkpt->nrh;i++){ 
			if((long)par->chkpt->co[i][ptSKY]==number_novalue){				
				r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
				c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				par->chkpt->co[i][ptSKY]=P->co[r][c];
			}
		}
		free_doublematrix(P);
	}	
	
	//f2. bedrock file
	read_bed=0;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptBED]==number_novalue) read_bed=1; 
	}
	if(read_bed==1 && coordinates==0) read_bed=0;
	if(read_bed==1){
		flag = file_exists(fbed, flog);	
		if(flag == 1){				
			if(read_dem==0){
				Q=new_doublematrix(1,1);
				P=read_map(0, files[fbed], Q, UV, (double)number_novalue);
				free_doublematrix(Q);
			}else{
				P=read_map(1, files[fbed], Z, UV, (double)number_novalue);
			}
		}else{
			printf("Warning: Bedrock depth file not present\n");
			fprintf(flog,"Warning: Sky view factor file not present\n");
			read_bed=0;
		}
	}
	
	if(read_bed==1){
		for(i=1;i<=par->chkpt->nrh;i++){ 
			if((long)par->chkpt->co[i][ptBED]==number_novalue){				
				r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
				c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				par->chkpt->co[i][ptBED]=P->co[r][c];
			}
		}
		free_doublematrix(P);
	}	
	
	//g.curvature
	read_curv=0;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if( (long)par->chkpt->co[i][ptCNS]==number_novalue || (long)par->chkpt->co[i][ptCWE]==number_novalue || 
		    (long)par->chkpt->co[i][ptCNwSe]==number_novalue || (long)par->chkpt->co[i][ptCNeSw]==number_novalue  ) read_curv=1; 
	}
	if(read_curv==1 && coordinates==0) read_curv=0;
	if(read_curv==1){
		if(read_dem==0){
			printf("Warning: Dem file is not present, and therefore it is not possible to calculate curvature\n");
			fprintf(flog,"Warning: Dem file is not present, and therefore it is not possible to calculate curvature\n");
			read_curv=0;
		}else{
			Q=new_doublematrix(Z->nrh,Z->nch);
			
			P=new_doublematrix(Z->nrh,Z->nch);
			R=new_doublematrix(Z->nrh,Z->nch);
			S=new_doublematrix(Z->nrh,Z->nch);
			T=new_doublematrix(Z->nrh,Z->nch);	
			
			multipass_topofilter(par->lowpass_curvatures, Z, Q, (double)number_novalue, 1);
			curvature(UV->U->co[1], UV->U->co[2], Q, P, R, S, T, (double)number_novalue);
			free_doublematrix(Q);
			
			if(strcmp(files[fcurv],string_novalue) != 0){
				temp = join_strings(files[fcurv], "N-S");
				write_map(temp, 0, par->format_out, P, UV, number_novalue);
				free(temp);
				temp = join_strings(files[fcurv], "W-E");
				write_map(temp, 0, par->format_out, R, UV, number_novalue);
				free(temp);
				temp = join_strings(files[fcurv], "NW-SE");
				write_map(temp, 0, par->format_out, S, UV, number_novalue);
				free(temp);		
				temp = join_strings(files[fcurv], "NE-SW");
				write_map(temp, 0, par->format_out, T, UV, number_novalue);
				free(temp);							
			}
			
			find_min_max(P, (double)number_novalue, &max, &min);
			printf("Curvature N-S Min:%f  Max:%f \n",min,max);
			fprintf(flog,"Curvature N-S Min:%f  Max:%f \n",min,max);
			
			find_min_max(R, (double)number_novalue, &max, &min);
			printf("Curvature W-E Min:%f  Max:%f \n",min,max);
			fprintf(flog,"Curvature W-E Min:%f  Max:%f \n",min,max);
			
			find_min_max(S, (double)number_novalue, &max, &min);
			printf("Curvature NW-SE Min:%f  Max:%f \n",min,max);
			fprintf(flog,"Curvature NW-SE Min:%f  Max:%f \n",min,max);
			
			find_min_max(T, (double)number_novalue, &max, &min);
			printf("Curvature NE-SW Min:%f  Max:%f \n",min,max);
			fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
		}
	}
	if(read_curv==1){
		for(i=1;i<=par->chkpt->nrh;i++){ 
			r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
			c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
			if((long)par->chkpt->co[i][ptCNS]==number_novalue) par->chkpt->co[i][ptCNS]=P->co[r][c];
			if((long)par->chkpt->co[i][ptCWE]==number_novalue) par->chkpt->co[i][ptCWE]=R->co[r][c];
			if((long)par->chkpt->co[i][ptCNwSe]==number_novalue) par->chkpt->co[i][ptCNwSe]=S->co[r][c];
			if((long)par->chkpt->co[i][ptCNeSw]==number_novalue) par->chkpt->co[i][ptCNeSw]=T->co[r][c];
		}
		free_doublematrix(P);
		free_doublematrix(R);
		free_doublematrix(S);
		free_doublematrix(T);
	}
	
	//h. no value check
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if((long)par->chkpt->co[i][ptZ]==number_novalue) par->chkpt->co[i][ptZ]=0.0;
		if((long)par->chkpt->co[i][ptLC]==number_novalue) par->chkpt->co[i][ptLC]=1.0;
		if((long)par->chkpt->co[i][ptSY]==number_novalue) par->chkpt->co[i][ptSY]=1.0;
		if((long)par->chkpt->co[i][ptS]==number_novalue) par->chkpt->co[i][ptS]=0.0;
		if((long)par->chkpt->co[i][ptA]==number_novalue) par->chkpt->co[i][ptA]=0.0;
		if((long)par->chkpt->co[i][ptSKY]==number_novalue) par->chkpt->co[i][ptSKY]=1.0;
		if((long)par->chkpt->co[i][ptCNS]==number_novalue) par->chkpt->co[i][ptCNS]=0.0;
		if((long)par->chkpt->co[i][ptCWE]==number_novalue) par->chkpt->co[i][ptCWE]=0.0;
		if((long)par->chkpt->co[i][ptCNwSe]==number_novalue) par->chkpt->co[i][ptCNwSe]=0.0;
		if((long)par->chkpt->co[i][ptCNeSw]==number_novalue) par->chkpt->co[i][ptCNeSw]=0.0;
		if((long)par->chkpt->co[i][ptDrDEPTH]==number_novalue) par->chkpt->co[i][ptDrDEPTH]=par->DepthFreeSurface;//[mm]		
		if((long)par->chkpt->co[i][ptMAXSWE]==number_novalue) par->chkpt->co[i][ptMAXSWE]=1.E10;//[mm]		
		if((long)par->chkpt->co[i][ptLAT]==number_novalue) par->chkpt->co[i][ptLAT]=par->latitude;		
		if((long)par->chkpt->co[i][ptLON]==number_novalue) par->chkpt->co[i][ptLON]=par->longitude;		
		if((long)par->chkpt->co[i][ptID]==number_novalue) par->chkpt->co[i][ptID]=(double)i;		
		if((long)par->chkpt->co[i][ptHOR]==number_novalue) par->chkpt->co[i][ptHOR]=par->chkpt->co[i][ptID];	
	}
	
	//i.show results
	fprintf(flog,"\nPOINTS:\n");
	fprintf(flog,"ID,East[m],North[m],Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-],CurvatureN-S[1/m],CurvatureW-E[1/m],CurvatureNW-SE[1/m],CurvatureNE-SW[1/m],DepthFreeSurface[mm],Hor,maxSWE[mm],Lat[deg],Long[deg]\n");
	for(r=1;r<=par->chkpt->nrh;r++){
		for(c=1;c<=ptTOT;c++){
			fprintf(flog,"%f",par->chkpt->co[r][c]);
			if (c<ptTOT) fprintf(flog, ",");
		}
		fprintf(flog,"\n");
	}	
	
	//l. set UV
	if(read_dem==0 && read_lu==0 && read_soil==0 && read_sl==0 && read_as==0 && read_sk==0){
		UV->U=new_doublevector(4);
		UV->V=new_doublevector(2);
	}
	UV->U->co[2]=1.0;
	UV->U->co[1]=1.0;
	UV->U->co[4]=0.0;
	UV->U->co[3]=0.0;
	UV->V->co[2]=(double)number_novalue;
	if(UV->V->co[2]<0){
		UV->V->co[1] = -1.;
	}else{
		UV->V->co[1] = 1.;
	}
	
	//n. deallocation
	if(read_dem==1) free_doublematrix(Z);
	if(read_lu==1) free_doublematrix(LU);
	//if(par->recover==0 && read_dem==1){
	if(read_dem==1){
		free_longvector(par->r_points);
		free_longvector(par->c_points);
	}
			
	//5. SET CHECKPOINT
	if(par->state_pixel == 1){
		par->rc=new_longmatrix(par->chkpt->nrh,2);
		for(i=1;i<=par->chkpt->nrh;i++){
			par->rc->co[i][1]=1;
			par->rc->co[i][2]=i;
		}
	}
	
	//6. SET PROPERTIES
	top->East=new_doublematrix(1,par->chkpt->nrh);
	top->North=new_doublematrix(1,par->chkpt->nrh);
	top->Z0=new_doublematrix(1,par->chkpt->nrh);
	land->LC=new_doublematrix(1,par->chkpt->nrh);	
	land->delay=new_doublematrix(1,par->chkpt->nrh);
	sl->type=new_longmatrix(1,par->chkpt->nrh);	
	top->slope=new_doublematrix(1,par->chkpt->nrh);
	top->aspect=new_doublematrix(1,par->chkpt->nrh);	
	top->curvature1=new_doublematrix(1,par->chkpt->nrh);	
	top->curvature2=new_doublematrix(1,par->chkpt->nrh);	
	top->curvature3=new_doublematrix(1,par->chkpt->nrh);	
	top->curvature4=new_doublematrix(1,par->chkpt->nrh);	
	top->sky=new_doublematrix(1,par->chkpt->nrh);
	top->pixel_type=new_shortmatrix(1,par->chkpt->nrh);
	top->BC_counter=new_longmatrix(1,par->chkpt->nrh);
	top->BC_DepthFreeSurface=new_doublevector(par->chkpt->nrh);
	par->maxSWE=new_doublematrix(1,par->chkpt->nrh);
	top->horizon_point=new_longmatrix(1,par->chkpt->nrh);
	top->dzdE=new_doublematrix(1,par->chkpt->nrh);
	top->dzdN=new_doublematrix(1,par->chkpt->nrh);
	top->latitude=new_doublematrix(1,par->chkpt->nrh);
	top->longitude=new_doublematrix(1,par->chkpt->nrh);
	par->IDpoint=new_longvector(par->chkpt->nrh);
	IT->bed=new_doublematrix(1,par->chkpt->nrh);
	for(i=1;i<=par->chkpt->nrh;i++){
		top->East->co[1][i]=par->chkpt->co[i][ptX];
		top->North->co[1][i]=par->chkpt->co[i][ptY];
		top->Z0->co[1][i]=par->chkpt->co[i][ptZ];
		land->LC->co[1][i]=par->chkpt->co[i][ptLC];
		
		if((long)land->LC->co[1][i] <= 0){
			f = fopen(FailedRunFile, "w");
			fprintf(f, "Error:: Point %ld has land cover type <= 0. This is not admitted.\n",i);
			fclose(f);
			t_error("Fatal Error! Geotop is closed. See failing report (15).");	
		}
		
		sl->type->co[1][i]=(long)par->chkpt->co[i][ptSY];
		
		if(sl->type->co[1][i] <= 0){
			f = fopen(FailedRunFile, "w");
			fprintf(f, "Error:: Point %ld has soil type <= 0. This is not admitted.\n",i);
			fclose(f);
			t_error("Fatal Error! Geotop is closed. See failing report (16).");	
		}
		
		top->slope->co[1][i]=par->chkpt->co[i][ptS];
		top->aspect->co[1][i]=par->chkpt->co[i][ptA];
		top->sky->co[1][i]=par->chkpt->co[i][ptSKY];
		top->curvature1->co[1][i]=par->chkpt->co[i][ptCNS];
		top->curvature2->co[1][i]=par->chkpt->co[i][ptCWE];
		top->curvature3->co[1][i]=par->chkpt->co[i][ptCNwSe];
		top->curvature4->co[1][i]=par->chkpt->co[i][ptCNeSw];
		top->pixel_type->co[1][i]=1;
		top->BC_counter->co[1][i]=i;
		top->BC_DepthFreeSurface->co[i]=par->chkpt->co[i][ptDrDEPTH];
		top->horizon_point->co[1][i]=(long)par->chkpt->co[i][ptHOR];
		top->dzdE->co[1][i]=0.;
		top->dzdN->co[1][i]=0.;		
		land->delay->co[1][i]=0.;
		
		if(sl->type->co[1][i] <= 0){
			f = fopen(FailedRunFile, "w");
			fprintf(f, "Error:: Point %ld has horizon type <= 0. This is not admitted.\n",i);
			fclose(f);
			t_error("Fatal Error! Geotop is closed. See failing report (17).");	
		}
		
		par->maxSWE->co[1][i]=par->chkpt->co[i][ptMAXSWE];
		top->latitude->co[1][i]=par->chkpt->co[i][ptLAT];
		top->longitude->co[1][i]=par->chkpt->co[i][ptLON];
		par->IDpoint->co[i]=(long)par->chkpt->co[i][ptID];
		
		IT->bed->co[1][i]=par->chkpt->co[i][ptBED];
		if( (long)IT->bed->co[1][i] == number_novalue ) IT->bed->co[1][i] = 1.E99;
	}
		
	//7. SET PAR
	for (i=1; i<=par->init_date->nh; i++) {
		par->output_soil->co[i]=0.;
		par->output_snow->co[i]=0.;
		par->output_glac->co[i]=0.;
		par->output_surfenergy->co[i]=0.;
		par->output_vegetation->co[i]=0.;
		par->output_meteo->co[i]=0.;
	}
	
	par->output_soil_bin = 0;
	par->output_snow_bin = 0;
	par->output_glac_bin = 0;
	par->output_surfenergy_bin = 0;
	par->output_meteo_bin = 0;
	
	//8. READ HORIZONS
	//find max top->horizon_point
	top->num_horizon_point=0;
	for(r=1;r<=top->horizon_point->nrh;r++){
		for(c=1;c<=top->horizon_point->nch;c++){
			if (top->horizon_point->co[r][c] > top->num_horizon_point) top->num_horizon_point = top->horizon_point->co[r][c];
		}
	}
	top->horizon_height=(double ***)malloc(top->num_horizon_point*sizeof(double**));	
	top->horizon_numlines=(long *)malloc(top->num_horizon_point*sizeof(long));
	for(i=1; i<=top->num_horizon_point; i++){
		
		c=0;
		do{
			flag = 0;
			if (c < par->chkpt->nrh) {
				if (top->horizon_point->co[1][c+1] != i) c++;
			}
			if (c < par->chkpt->nrh) {
				if (top->horizon_point->co[1][c+1] != i) flag=1;
			}
			
		}while (flag == 1 && c < par->chkpt->nrh);
				
		if (c < par->chkpt->nrh) {
			top->horizon_height[i-1] = read_horizon(0, i, files[fhorpoint], IT->horizon_col_names, &num_lines, flog);
			top->horizon_numlines[i-1] = num_lines;
		}else {
			top->horizon_height[i-1] = (double**)malloc(sizeof(double*));
			top->horizon_height[i-1][0] = (double*)malloc(2.*sizeof(double));
			top->horizon_height[i-1][0][0] = 0.;
			top->horizon_height[i-1][0][1] = 0.;
			top->horizon_numlines[i-1] = 1;
		}
	}
		
}
	
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void set_bedrock(INIT_TOOLS *IT, SOIL *sl, CHANNEL *cnet, PAR *par, TOPO *top, DOUBLEMATRIX *LC, FILE *flog){
	
	DOUBLETENSOR *T;
	DOUBLEVECTOR *WT;
	long i, j, l, r, c, sy, synew;
	double zlim, z;
	short yes=0;
	FILE *f;
	
	//check if bedrock depth is above soil lower border, otherwise we do not need to calculate anything
	z = 0.;
	for (l=1; l<=Nl; l++) {
		z += sl->pa->co[1][jdz][l];
	}
	for (i=1; i<=par->total_pixel; i++){
		r = top->rc_cont->co[i][1];
		c = top->rc_cont->co[i][2];
		if (IT->bed->co[r][c] < z) yes = 1;
	}
	
	if (yes == 1){
	
		//consistency check
		if (IT->init_water_table_depth->nh != sl->pa->ndh){
			f = fopen(FailedRunFile, "w");
			fprintf(f, "Error:: Error in bedrock calculations");
			fclose(f);
			t_error("Fatal Error! Geotop is closed. See failing report (19).");	
		}
		
		//rewrite soil type
		T=new_doubletensor(sl->pa->ndh, nsoilprop, Nl);
		for(i=1;i<=sl->pa->ndh;i++){
			for(j=1;j<=nsoilprop;j++){
				for(l=1;l<=Nl;l++){		
					T->co[i][j][l]=sl->pa->co[i][j][l];
				}
			}
		}	
		free_doubletensor(sl->pa);
		sl->pa=new_doubletensor(par->total_pixel+par->total_channel, nsoilprop, Nl);
		
		//rewrite initial water table depth
		WT=new_doublevector(IT->init_water_table_depth->nh);
		for(i=1;i<=IT->init_water_table_depth->nh;i++) {
			WT->co[i]=IT->init_water_table_depth->co[i];
		}
		free_doublevector(IT->init_water_table_depth);
		IT->init_water_table_depth=new_doublevector(par->total_pixel+par->total_channel);
		
		//assign jdz (is needed later)
		for(i=1;i<=sl->pa->ndh;i++){
			for(l=1;l<=Nl;l++){		
				sl->pa->co[i][jdz][l]=T->co[1][jdz][l];
			}
		}	
		
		for (i=1; i<=par->total_pixel+par->total_channel; i++) {
			
			if (i<=par->total_pixel) {
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				sy = sl->type->co[r][c];
				synew = i;
				sl->type->co[r][c] = synew;
				z = 0.;
			}else {
				r = cnet->r->co[i-par->total_pixel];
				c = cnet->c->co[i-par->total_pixel];
				sy = cnet->soil_type->co[i-par->total_pixel];
				synew = i;
				cnet->soil_type->co[i-par->total_pixel] = synew;
				z = par->depr_channel;
			}
			
			IT->init_water_table_depth->co[synew] = WT->co[sy];
			
			zlim = IT->bed->co[r][c];
			
			for(l=1;l<=Nl;l++){			
				
				z += 0.5*sl->pa->co[synew][jdz][l];
				
				if(z <= zlim){
					
					for(j=1;j<=nsoilprop;j++){
						sl->pa->co[synew][j][l] = T->co[sy][j][l];
					}
					
				}else{
					
					for(j=1;j<=nsoilprop;j++){
						sl->pa->co[synew][j][l] = IT->pa_bed->co[sy][j][l] ;
					}				
					
				}
				
				z += 0.5*sl->pa->co[synew][jdz][l];
				
			}
		}
		
		free_doubletensor(T);
		free_doublevector(WT);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLETENSOR *find_Z_of_any_layer(DOUBLEMATRIX *Zsurface, DOUBLEMATRIX *slope, DOUBLEMATRIX *LC, SOIL *sl, short point){
	
	DOUBLETENSOR *Z;
	double Zaverage=0., z, cosine;
	long l, r, c, n, sy;
	
	if(point!=1){
		Zaverage=0.;
		n=0;
		for(r=1;r<=Zsurface->nrh;r++){
			for(c=1;c<=Zsurface->nch;c++){
				if((long)LC->co[r][c]!=number_novalue){
					n++;
					Zaverage += Zsurface->co[r][c];
				}
			}
		}
		Zaverage/=(double)n;
	}

	
	Z=new_doubletensor0(sl->pa->nch, Zsurface->nrh, Zsurface->nch);
	initialize_doubletensor(Z, (double)number_novalue);

	for(r=1;r<=Zsurface->nrh;r++){
		for(c=1;c<=Zsurface->nch;c++){
			if((long)LC->co[r][c]!=number_novalue){
				
				cosine = cos(slope->co[r][c]*Pi/180.);
						
				sy=sl->type->co[r][c];
				if (point!=1){
					z=1.E3*(Zsurface->co[r][c]-Zaverage);//[mm]
				}else {
					z=0.;
				}
			
				l=0;
				Z->co[l][r][c]=z;
			
				do{
					l++;
					z -= 0.5*sl->pa->co[sy][jdz][l]*cosine;
					Z->co[l][r][c]=z;
					z -= 0.5*sl->pa->co[sy][jdz][l]*cosine;
				}while(l<sl->pa->nch);
			}
		}
	}
	
	return Z;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short file_exists(short key, FILE *flog){
	
	//no keyword -> -1
	//keyword but does not exist -> 0
	//keyword and exists -> 1
	
	printf("Attempting to read '%s' in the file '%s': ",keywords_char[key+nmet+nsoilprop+2],files[key]);
	fprintf(flog,"Attempting to read '%s' in the file '%s': ",keywords_char[key+nmet+nsoilprop+2],files[key]);
	
	if(strcmp(files[key] , string_novalue) == 0){
		printf("not present in file list\n");
		fprintf(flog,"not present in file list\n");		
		return (-1);
	}else{
		if(existing_file(files[key] )>0){ 
			printf("EXISTING in format %d\n",existing_file(files[key]));
			fprintf(flog, "EXISTING in format %d\n",existing_file(files[key]));
			return (1);
		}else{
			printf("not existing\n");
			fprintf(flog, "not existing\n");
			return (0);
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double peat_thickness(double dist_from_channel){
	
	double D;
	
	if(dist_from_channel<45.23){
		D = 10.*(47.383 - 0.928*dist_from_channel + 0.010*pow(dist_from_channel,2.));
	}else{
		D = 10.*26.406;
	}
	
	return(D);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initialize_soil_state(SOIL_STATE *S, long n, long nl){

	S->T = new_doublematrix(nl, n);
	initialize_doublematrix(S->T, 0.);
	S->P = new_doublematrix0_(nl, n);
	initialize_doublematrix(S->P, 0.);
	S->thi = new_doublematrix(nl, n);
	initialize_doublematrix(S->thi, 0.);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void copy_soil_state(SOIL_STATE *from, SOIL_STATE *to){
	
	long l,i;
	long nl=from->T->nrh,n=from->T->nch;
	
	for (i=1; i<=n; i++) {
		to->P->co[0][i] = from->P->co[0][i];
		for (l=1; l<=nl; l++) {
			to->P->co[l][i] = from->P->co[l][i];
			to->T->co[l][i] = from->T->co[l][i];
			to->thi->co[l][i] = from->thi->co[l][i];
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initialize_veg_state(STATE_VEG *V, long n){
	
	V->Tv = new_doublevector(n);
	initialize_doublevector(V->Tv, 0.);
	V->wsnow = new_doublevector(n);
	initialize_doublevector(V->wsnow, 0.);
	V->wrain = new_doublevector(n);
	initialize_doublevector(V->wrain, 0.);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void copy_veg_state(STATE_VEG *from, STATE_VEG *to){
	
	long i, n=from->Tv->nh;
	for (i=1; i<=n; i++) {
		to->Tv->co[i] = from->Tv->co[i];
		to->wrain->co[i] = from->wrain->co[i];
		to->wsnow->co[i] = from->wsnow->co[i];
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

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
	char *path_rec_files;
	
	char **keywords_num_lower_case, **keywords_char_lower_case;
	
	long beg=0, end=0;
	
	//convert keyword listed on top of the file in lower case
	n = (long)num_par_number;
	keywords_num_lower_case = (char**)malloc(n*sizeof(char*));
	for (i=0; i<n; i++) {
		//printf("%ld,%s\n",i,keywords_num[i]);
		keywords_num_lower_case[i] = assign_string(keywords_num[i]);
		convert_string_in_lower_case(keywords_num_lower_case[i]);
	}
	
	n = (long)num_par_char;
	keywords_char_lower_case = (char**)malloc(n*sizeof(char*));
	for (i=0; i<n; i++) {
		//printf("%ld,%s\n",i,keywords_char[i]);
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
	if (strcmp(string_novalue,temp) == 0) {
		free(temp);
		temp = assign_string("_SUCCESSFUL_RUN");
	}
	SuccessfulRunFile = join_strings(WORKING_DIRECTORY, temp);
	free(temp);
			
	beg = end;
	end += 1;
	temp = assignation_string(flog, beg, keywords_char, string_param);
	if (strcmp(string_novalue,temp) == 0) {
		free(temp);
		temp = assign_string("_FAILED_RUN");
	}	
	FailedRunFile = join_strings(WORKING_DIRECTORY, temp);
	free(temp);
	
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
			}else if(j==6){ hsnw[j] = assign_string("Temperature[C]");     	
			}else if(j==7){ hsnw[j] = assign_string("wice[kg/m2]");     	
			}else if(j==8){ hsnw[j] = assign_string("wliq[kg/m2]");     	
			}else if(j==9){ hsnw[j] = assign_string("Dz[mm]");
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
		if(strcmp(files[j], string_novalue) == 0){
			free(files[j]);
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
	if(strcmp(path_rec_files, string_novalue) != 0){
		temp = assign_string(path_rec_files);
		free(path_rec_files);
		path_rec_files = join_strings(temp, "/");
		free(temp);
		for (i=rpsi; i<=rsux; i++) {
			temp = assign_string(files[i]);
			free(files[i]);
			files[i] = join_strings(path_rec_files, temp);
			free(temp);
		}
	}
	free(path_rec_files);

	//add working path to the file name
	for (i=0; i<nfiles; i++) {
		if (strcmp(files[i], string_novalue) != 0) {
			temp = assign_string(files[i]);
			free(files[i]);
			files[i] = join_strings(WORKING_DIRECTORY, temp);
			free(temp);
		}
	}

	
	return 1;
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

void assign_numeric_parameters(PAR *par, LAND *land, TIMES *times, SOIL *sl, METEO *met, INIT_TOOLS *itools, double **num_param, long *num_param_components, char **keyword, FILE *flog){
	
	short occurring;
	long cod, codn, i, j, k, n, m, nsoillayers, nmeteo_stations, npoints;
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
	for (i=1; i<=par->init_date->nh; i++) {
		par->init_date->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 010119000000., 0);
		par->init_date->co[i] = convert_dateeur12_JDfrom0(par->init_date->co[i]);
	}
	
	//simulation time
	cod = 398;
	par->simulation_hours = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 1., 0);
	
	//end date
	cod = 2;
	par->end_date = new_doublevector(num_param_components[cod]);
	
	if (par->end_date->nh != par->init_date->nh){
		fprintf(flog, "Error:: End date has a number of components different from Init Date");
		printf("Error:: End date has a number of components different from Init Date");
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}
		
	for (i=1; i<=par->end_date->nh; i++) {
		par->end_date->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
		if ((long)par->end_date->co[i] == number_novalue){
			par->end_date->co[i] = par->init_date->co[i] + par->simulation_hours/24.;
		}else {
			par->end_date->co[i] = convert_dateeur12_JDfrom0(par->end_date->co[i]);
		}		
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
	par->Dtplot_discharge->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->Dtplot_discharge->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_discharge->co[i-1], 0);
	}
	par->plot_discharge_with_Dt_integration = new_shortvector(par->init_date->nh);
	par->state_discharge = 0;
	for (i=1; i<=par->init_date->nh; i++) {
		par->Dtplot_discharge->co[i] *= 3600.;
		if(par->Dtplot_discharge->co[i] > 1.E-5 && par->Dtplot_discharge->co[i] <= minDt){
			par->plot_discharge_with_Dt_integration->co[i]=1;
		}else{
			par->plot_discharge_with_Dt_integration->co[i]=0;
		}
		if(par->Dtplot_discharge->co[i] > 1.E-5) par->state_discharge = 1;
	}
	
	cod = 6;
	par->Dtplot_point = new_doublevector(par->init_date->nh);
	par->Dtplot_point->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->Dtplot_point->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_point->co[i-1], 0);
	}
	par->plot_point_with_Dt_integration = new_shortvector(par->init_date->nh);
	par->state_pixel = 0;
	for (i=1; i<=par->init_date->nh; i++) {
		par->Dtplot_point->co[i] *= 3600.;
		if(par->Dtplot_point->co[i] > 1.E-5 && par->Dtplot_point->co[i] <= minDt){
			par->plot_point_with_Dt_integration->co[i]=1;
		}else{
			par->plot_point_with_Dt_integration->co[i]=0;
		}
		if(par->Dtplot_point->co[i] > 1.E-5) par->state_pixel = 1;
	}
	
	cod = 7;
	par->Dtplot_basin = new_doublevector(par->init_date->nh);
	par->Dtplot_basin->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->Dtplot_basin->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Dtplot_basin->co[i-1], 0);
	}
	par->plot_basin_with_Dt_integration = new_shortvector(par->init_date->nh);
	par->state_basin = 0;
	for (i=1; i<=par->init_date->nh; i++) {
		par->Dtplot_basin->co[i] *= 3600.;
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
	par->imp = assignation_number(flog, 40, 0, keyword, num_param, num_param_components, 7., 0);
	par->free_drainage_bottom = assignation_number(flog, 41, 0, keyword, num_param, num_param_components, 0., 0);
	par->free_drainage_lateral = assignation_number(flog, 42, 0, keyword, num_param, num_param_components, 1., 0);
	par->TolVWb = assignation_number(flog, 43, 0, keyword, num_param, num_param_components, 1.E-6, 0);
	par->RelTolVWb = RelativeErrorRichards;
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
	par->max_weq_snow = assignation_number(flog, 118, 0, keyword, num_param, num_param_components, 5., 0);
	n = (long)assignation_number(flog, 119, 0, keyword, num_param, num_param_components, 2., 0);
	if(n < 1){
		fprintf(flog,"Error:: %s must be 1 or larger\n",keyword[119]);
		printf("Error:: %s must be 1 or larger\n",keyword[119]);
		t_error("Fatal Error! Geotop is closed.");	
	}	
	par->SWE_bottom = assignation_number(flog, 120, 0, keyword, num_param, num_param_components, 20., 0);
	par->SWE_top = assignation_number(flog, 121, 0, keyword, num_param, num_param_components, 20., 0);
	par->max_snow_layers = (long)floor(par->SWE_bottom/par->max_weq_snow) + (long)floor(par->SWE_top/par->max_weq_snow) + n;
	par->inf_snow_layers = new_longvector(n);
	fprintf(flog,"Max snow layer number: %ld, of which %.0f at the bottom, %ld in the middle, and %.0f at the top.\n",par->max_snow_layers,floor(par->SWE_bottom/par->max_weq_snow),n,floor(par->SWE_top/par->max_weq_snow));
	fprintf(flog,"Infinite Snow layer numbers are numbers: ");
	for (i=1; i<=n; i++) {
		par->inf_snow_layers->co[i] = (long)floor(par->SWE_bottom/par->max_weq_snow) + i;
		fprintf(flog, "%ld ",par->inf_snow_layers->co[i]);
	}
	fprintf(flog,"\n");
	
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
	
	if(n < 1 && (par->GWE_bottom > 0 || par->GWE_top > 0)){
		fprintf(flog,"Error:: %s must be 1 or larger\n",keyword[127]);
		printf("Error:: %s must be 1 or larger\n",keyword[127]);
		t_error("Fatal Error! GEOtop is closed.");	
	}
	
	par->max_glac_layers = (long)floor(par->GWE_bottom/par->max_weq_glac) + (long)floor(par->GWE_top/par->max_weq_glac) + n;
	par->inf_glac_layers = new_longvector(n);
	fprintf(flog,"Max glac layer number: %ld, of which %.0f at the bottom, %ld in the middle, and %.0f at the top.\n",par->max_glac_layers,floor(par->GWE_bottom/par->max_weq_glac),n,floor(par->GWE_top/par->max_weq_glac));
	fprintf(flog,"Infinite Glac layer numbers are numbers: ");
	for (i=1; i<=n; i++) {
		par->inf_glac_layers->co[i] = (long)floor(par->GWE_bottom/par->max_weq_glac) + i;
		fprintf(flog, "%ld ",par->inf_glac_layers->co[i]);
	}
	fprintf(flog,"\n");
	
	//former block 9
	par->state_turb = 1;
	par->state_lwrad = (short)assignation_number(flog, 130, 0, keyword, num_param, num_param_components, 9., 0);
	par->monin_obukhov = (short)assignation_number(flog, 131, 0, keyword, num_param, num_param_components, 1., 0);
	par->surroundings = (short)assignation_number(flog, 132, 0, keyword, num_param, num_param_components, 0., 0);
		
	//distributed option file
	par->wat_balance = (short)assignation_number(flog, 133, 0, keyword, num_param, num_param_components, 0., 0);
	par->en_balance = (short)assignation_number(flog, 134, 0, keyword, num_param, num_param_components, 0., 0);
	par->blowing_snow = (short)assignation_number(flog, 135, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->Wmin_BS = assignation_number(flog, 136, 0, keyword, num_param, num_param_components, 8., 0);

	cod = 137;
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
		
	cod = 156;
	par->saving_points = new_doublevector(num_param_components[cod]);
	for (i=1; i<=par->saving_points->nh; i++) {
		par->saving_points->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 0);
	}
			
	par->output_soil = new_doublevector(par->init_date->nh);
	par->output_snow = new_doublevector(par->init_date->nh);
	par->output_glac = new_doublevector(par->init_date->nh);
	par->output_surfenergy = new_doublevector(par->init_date->nh);
	par->output_vegetation = new_doublevector(par->init_date->nh);
	par->output_meteo = new_doublevector(par->init_date->nh);
	
	cod = 157;
	par->output_soil->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->output_soil->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_soil->co[i-1], 0);
	}
	
	cod = 158;
	par->output_snow->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->output_snow->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_snow->co[i-1], 0);
	}
	
	cod = 159;
	par->output_glac->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->output_glac->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_glac->co[i-1], 0);
	}	
	
	cod = 160;
	par->output_surfenergy->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->output_surfenergy->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_surfenergy->co[i-1], 0);
	}	
	
	cod = 161;
	par->output_vegetation->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->output_vegetation->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_vegetation->co[i-1], 0);
	}
	
	cod = 162;
	par->output_meteo->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=par->init_date->nh; i++) {
		par->output_meteo->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->output_meteo->co[i-1], 0);
	}
	
	par->output_soil_bin = 0;
	par->output_snow_bin = 0;
	par->output_glac_bin = 0;
	par->output_surfenergy_bin = 0;
	par->output_meteo_bin = 0;

	for (i=1; i<=par->init_date->nh; i++) {
		if (par->output_soil->co[i] > 0) par->output_soil_bin = 1;		
		if (par->output_snow->co[i] > 0) par->output_snow_bin = 1;
		if (par->output_glac->co[i] > 0) par->output_glac_bin = 1;
		if (par->output_surfenergy->co[i] > 0) par->output_surfenergy_bin = 1;
		if (par->output_meteo->co[i] > 0) par->output_meteo_bin = 1;
	}

	cod = 163;
	codn = 164;
	
	if(num_param_components[cod] != num_param_components[codn]){
		fprintf(flog, "Error:: Number of components of parameters %s and %s must be equal\n",keyword[cod],keyword[codn]);
		printf("Error:: Number of components of parameters %s and %s must be equal\n",keyword[cod],keyword[codn]);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
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
	par->nsoiltypes = (long)assignation_number(flog, 165, 0, keyword, num_param, num_param_components, 1., 0);
	if (par->nsoiltypes < 1) par->nsoiltypes = 1;
	
	cod = 166;
	itools->init_water_table_depth = new_doublevector(par->nsoiltypes);
	itools->init_water_table_depth->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 5000., 0);
	for (k=2; k<=par->nsoiltypes; k++) {
		itools->init_water_table_depth->co[k] = assignation_number(flog, cod, k-1, keyword, num_param, num_param_components, itools->init_water_table_depth->co[k-1], 0);
	}
	
	//soil properties and discretization
	par->soil_type_land_default = (long)assignation_number(flog, 167, 0, keyword, num_param, num_param_components, 1., 0);
	
	if(par->soil_type_land_default<1 || par->soil_type_land_default>par->nsoiltypes){
		fprintf(flog, "Error:: Soil_type_land_default lower than 0 or higher than soil types numbers");
		printf("Error:: Soil_type_land_default lower than 0 or higher than soil types numbers");
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}
	
	par->soil_type_chan_default = (long)assignation_number(flog, 168, 0, keyword, num_param, num_param_components, 1., 0);
	
	if(par->soil_type_chan_default<1 || par->soil_type_chan_default>par->nsoiltypes){
		fprintf(flog, "Error:  Soil_type_chan_default lower than 0 or higher than soil types numbers");
		printf("Error:  Soil_type_chan_default lower than 0 or higher than soil types numbers");
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}

	cod = 169;
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
		nsoillayers = (long)assignation_number(flog, 170, 0, keyword, num_param, num_param_components, 5., 0);
		
		sl->pa = new_doubletensor(par->nsoiltypes, nsoilprop, nsoillayers);
		
		for (i=1; i<=sl->pa->nch; i++) {
			sl->pa->co[1][jdz][i] = a;
		}

	}
	
	//first layer
	i = 1;
	sl->pa->co[1][jpsi][i] = assignation_number(flog, 171, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	sl->pa->co[1][jT][i] = assignation_number(flog, 172, i-1, keyword, num_param, num_param_components, 5., 0);
	sl->pa->co[1][jKn][i] = assignation_number(flog, 173, i-1, keyword, num_param, num_param_components, 1.E-4, 0); 
	sl->pa->co[1][jKl][i] = assignation_number(flog, 174, i-1, keyword, num_param, num_param_components, 1.E-4, 0); 
	sl->pa->co[1][jres][i] = assignation_number(flog, 175, i-1, keyword, num_param, num_param_components, 0.05, 0);  
	sl->pa->co[1][jwp][i] = assignation_number(flog, 176, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);  
	sl->pa->co[1][jfc][i] = assignation_number(flog, 177, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	sl->pa->co[1][jsat][i] = assignation_number(flog, 178, i-1, keyword, num_param, num_param_components, 0.5, 0);  
	sl->pa->co[1][ja][i] = assignation_number(flog, 179, i-1, keyword, num_param, num_param_components, 0.004, 0);  
	sl->pa->co[1][jns][i] = assignation_number(flog, 180, i-1, keyword, num_param, num_param_components, 1.3, 0);  
	sl->pa->co[1][jv][i] = assignation_number(flog, 181, i-1, keyword, num_param, num_param_components, 0.5, 0);  
	sl->pa->co[1][jkt][i] = assignation_number(flog, 182, i-1, keyword, num_param, num_param_components, 2.5, 0);  
	sl->pa->co[1][jct][i] = assignation_number(flog, 183, i-1, keyword, num_param, num_param_components, 1.E6, 0);  
	sl->pa->co[1][jss][i] = assignation_number(flog, 184, i-1, keyword, num_param, num_param_components, 1.E-7, 0);  

	//other layers
	for (i=2; i<=sl->pa->nch; i++) {
		for (j=2; j<=sl->pa->nrh; j++) {
			sl->pa->co[1][j][i] = assignation_number(flog, 171 + j-2, i-1, keyword, num_param, num_param_components, sl->pa->co[1][j][i-1], 0);
		}
	}
	
	//field capacity (-0.333 bar) and wilting point (-15 bar)
	for (i=1; i<=sl->pa->nch; i++){
		if( (long)sl->pa->co[1][jfc][i] == number_novalue){
			sl->pa->co[1][jfc][i] = teta_psi( (-1./3.)*1.E5/g, 0., sl->pa->co[1][jsat][i], sl->pa->co[1][jres][i], sl->pa->co[1][ja][i], 
										 sl->pa->co[1][jns][i], 1.-1./sl->pa->co[1][jns][i], PsiMin, sl->pa->co[1][jss][i]);
		}
		
		if( (long)sl->pa->co[1][jwp][i] == number_novalue){
			sl->pa->co[1][jwp][i] = teta_psi( -15.*1.E5/g, 0., sl->pa->co[1][jsat][i], sl->pa->co[1][jres][i], sl->pa->co[1][ja][i], 
										 sl->pa->co[1][jns][i], 1.-1./sl->pa->co[1][jns][i], PsiMin, sl->pa->co[1][jss][i]);		
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
		if (occurring == 0) itools->init_water_table_depth->co[k] = (double)number_novalue;
	}
	
	cod = 185;
	itools->pa_bed = new_doubletensor(1, nsoilprop, nsoillayers);
	for (i=1; i<=nsoillayers; i++) {
		itools->pa_bed->co[1][jdz][i] = sl->pa->co[1][jdz][i];
	}
	for (j=1; j<=nsoilprop; j++) {
		if(j != jdz){
			itools->pa_bed->co[1][j][1] = assignation_number(flog, cod+j-2, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);
		}
	}
	for (i=2; i<=nsoillayers; i++) {
		for (j=1; j<=nsoilprop; j++) {
			if(j != jdz) itools->pa_bed->co[1][j][i] = assignation_number(flog, cod+j-2, i-1, keyword, num_param, num_param_components, itools->pa_bed->co[1][j][i-1], 0);
		}
	}	
	//field capacity (-0.333 bar) and wilting point (-15 bar)
	for (i=1; i<=sl->pa->nch; i++){
		if( (long)itools->pa_bed->co[1][jsat][i] != number_novalue && (long)itools->pa_bed->co[1][jres][i] != number_novalue && 
			(long)itools->pa_bed->co[1][ja][i] != number_novalue && (long)itools->pa_bed->co[1][jns][i] != number_novalue &&
			(long)itools->pa_bed->co[1][jss][i] ) {
			if ( (long)itools->pa_bed->co[1][jfc][i] == number_novalue) {
				itools->pa_bed->co[1][jfc][i] = teta_psi( (-1./3.)*1.E5/g, 0., itools->pa_bed->co[1][jsat][i], itools->pa_bed->co[1][jres][i], itools->pa_bed->co[1][ja][i], 
										 itools->pa_bed->co[1][jns][i], 1.-1./itools->pa_bed->co[1][jns][i], PsiMin, itools->pa_bed->co[1][jss][i]);
			}
			if ( (long)itools->pa_bed->co[1][jwp][i] == number_novalue) {
				itools->pa_bed->co[1][jwp][i] = teta_psi( -15.*1.E5/g, 0., itools->pa_bed->co[1][jsat][i], itools->pa_bed->co[1][jres][i], itools->pa_bed->co[1][ja][i], 
										 itools->pa_bed->co[1][jns][i], 1.-1./itools->pa_bed->co[1][jns][i], PsiMin, itools->pa_bed->co[1][jss][i]);
			}
		}	
	}
	
	//meteo stations
	cod = 199;
	met->imeteo_stations = new_longvector(num_param_components[cod]);
	met->imeteo_stations->co[1] = (long)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);
	if ( met->imeteo_stations->co[1] != number_novalue ) {
		for (i=2; i<=num_param_components[cod]; i++) {
			met->imeteo_stations->co[i] = (long)assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, 0., 1);
		}
		nmeteo_stations = num_param_components[cod];
	}else {
		nmeteo_stations = (long)assignation_number(flog, 200, 0, keyword, num_param, num_param_components, 1., 0);
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
	met->st->E->co[i] = assignation_number(flog, 201, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	met->st->N->co[i] = assignation_number(flog, 202, i-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	met->st->lat->co[i] = assignation_number(flog, 203, i-1, keyword, num_param, num_param_components, par->latitude, 0);
	met->st->lon->co[i] = assignation_number(flog, 204, i-1, keyword, num_param, num_param_components, par->longitude, 0);
	met->st->Z->co[i] = assignation_number(flog, 205, i-1, keyword, num_param, num_param_components, 0., 0);
	met->st->sky->co[i] = assignation_number(flog, 206, i-1, keyword, num_param, num_param_components, 1., 0);
	met->st->ST->co[i] = assignation_number(flog, 207, i-1, keyword, num_param, num_param_components, par->ST, 0);
	
	for (i=2; i<=nmeteo_stations; i++) {
		met->st->E->co[i] = assignation_number(flog, 201, i-1, keyword, num_param, num_param_components, met->st->E->co[i-1], 0);
		met->st->N->co[i] = assignation_number(flog, 202, i-1, keyword, num_param, num_param_components, met->st->N->co[i-1], 0);
		met->st->lat->co[i] = assignation_number(flog, 203, i-1, keyword, num_param, num_param_components, met->st->lat->co[i-1], 0);
		met->st->lon->co[i] = assignation_number(flog, 204, i-1, keyword, num_param, num_param_components, met->st->lon->co[i-1], 0);
		met->st->Z->co[i] = assignation_number(flog, 205, i-1, keyword, num_param, num_param_components, met->st->Z->co[i-1], 0);
		met->st->sky->co[i] = assignation_number(flog, 206, i-1, keyword, num_param, num_param_components, met->st->sky->co[i-1], 0);
		met->st->ST->co[i] = assignation_number(flog, 207, i-1, keyword, num_param, num_param_components, met->st->ST->co[i-1], 0);
 	}
		
	a = assignation_number(flog, 208, 0, keyword, num_param, num_param_components, 10., 0);
	for (i=1; i<=nmeteo_stations; i++) {
		met->st->Vheight->co[i] = a;
	}
	
	a = assignation_number(flog, 209, 0, keyword, num_param, num_param_components, 2., 0);
	for (i=1; i<=nmeteo_stations; i++) {
		met->st->Theight->co[i] = a;
	}
	
	//lapse rates (cyclic)
	n = (long)nlstot;
	met->LRcnc = (long*)malloc(n*sizeof(long));
	met->LRcnc[ilsDate12] = 1;
	met->LRcnc[ilsTa] = num_param_components[210];
	met->LRcnc[ilsTdew] = num_param_components[211];
	met->LRcnc[ilsPrec] = num_param_components[212];
	met->LRc = (double**)malloc(n*sizeof(double*));
	for (i=0; i<nlstot; i++) {
		met->LRc[i] = (double*)malloc(met->LRcnc[i]*sizeof(double));
		for (j=0; j<met->LRcnc[i]; j++) {
			if(i==ilsDate12) met->LRc[i][j] = 0.;
			if(i==ilsTa) met->LRc[i][j] = assignation_number(flog, 210, j, keyword, num_param, num_param_components, (double)number_novalue, 0);
			if(i==ilsTdew) met->LRc[i][j] = assignation_number(flog, 211, j, keyword, num_param, num_param_components, (double)number_novalue, 0);
			if(i==ilsPrec) met->LRc[i][j] = assignation_number(flog, 212, j, keyword, num_param, num_param_components, (double)number_novalue, 0);
		}
	}
	
	par->MinIncrFactWithElev = assignation_number(flog, 213, j, keyword, num_param, num_param_components, 0.1, 0);
	par->MaxIncrFactWithElev = assignation_number(flog, 214, j, keyword, num_param, num_param_components, 4.4, 0);
	
	
	//output point column
	n = (long)otot;
	opnt = (long*)malloc(n*sizeof(long));
	ipnt = (short*)malloc(n*sizeof(short));
	
	par->all_point = (short)assignation_number(flog, 294, 0, keyword, num_param, num_param_components, 0., 0);
	
	if (par->all_point == 1) {
		
		for (i=0; i<n; i++) {
			ipnt[i] = 1;
			opnt[i] = i;
		}
		
		nopnt = n;
		
	}else {
		
		for (i=0; i<n; i++) {
			ipnt[i] = 0;
			opnt[i] = -1;
		}
		for (i=0; i<n; i++) {
			j = (long)assignation_number(flog, 215+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n){
				opnt[j-1] = i;
				ipnt[i] = 1;
			}
		}
		nopnt = 0;
		for (i=0; i<n; i++) {
			if(opnt[i] > 0) nopnt = i+1;
		}
	}
	
	//output basin column
	n = (long)ootot;
	obsn = (long*)malloc(n*sizeof(long));
	ibsn = (short*)malloc(n*sizeof(short));
	
	par->all_basin = (short)assignation_number(flog, 322, 0, keyword, num_param, num_param_components, 0., 0);
	
	if (par->all_basin == 1) {
		
		for (i=0; i<n; i++) {
			ibsn[i] = 1;
			obsn[i] = i;
		}
		
		nobsn = (long)ootot;
		
	}else {
		
		for (i=0; i<n; i++) {
			ibsn[i] = 0;
			obsn[i] = -1;
		}
		for (i=0; i<n; i++) {
			j = (long)assignation_number(flog, 295+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n){
				obsn[j-1] = i;
				ibsn[i] = 1;
			}
		}
		nobsn = 0;
		for (i=0; i<n; i++) {
			if(obsn[i] > 0) nobsn = i+1;
		}
	}
	
	cod = 356;
	n = 0;
	do{
		a = assignation_number(flog, cod, n, keyword, num_param, num_param_components, (double)number_novalue, 0);
		if ((long)a != number_novalue) n++;
	}while ((long)a != number_novalue && n<=100);	
	if (n==0) n=1;
	par->soil_plot_depths = new_doublevector(n);
	for (n=1; n<=par->soil_plot_depths->nh; n++) {
		par->soil_plot_depths->co[n] = assignation_number(flog, cod, n-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	}	
	
	cod = 357;
	n = 0;
	do{
		a = assignation_number(flog, cod, n, keyword, num_param, num_param_components, (double)number_novalue, 0);
		if ((long)a != number_novalue) n++;
	}while ((long)a != number_novalue && n<=100);	
	if (n==0) n=1;
	par->snow_plot_depths = new_doublevector(n);
	for (n=1; n<=par->snow_plot_depths->nh; n++) {
		par->snow_plot_depths->co[n] = assignation_number(flog, cod, n-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	}	
	
	cod = 358;
	n = 0;
	do{
		a = assignation_number(flog, cod, n, keyword, num_param, num_param_components, (double)number_novalue, 0);
		if ((long)a != number_novalue) n++;
	}while ((long)a != number_novalue && n<=100);	
	if (n==0) n=1;
	par->glac_plot_depths = new_doublevector(n);
	for (n=1; n<=par->glac_plot_depths->nh; n++) {
		par->glac_plot_depths->co[n] = assignation_number(flog, cod, n-1, keyword, num_param, num_param_components, (double)number_novalue, 0);
	}	
	
	//output snow column
	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
		m = par->snow_plot_depths->nh;
	}else {
		m = par->max_snow_layers;
	}
	n = 6 + 3*m + 1*par->max_snow_layers;
	
	osnw = (long*)malloc(n*sizeof(long));
	
	par->all_snow = (short)assignation_number(flog, 335, 0, keyword, num_param, num_param_components, 0., 0);
	
	if (par->all_snow == 1) {
		
		for (i=0; i<n; i++) {
			osnw[i] = i;
		}
		
		nosnw = n;
		
	}else {
		
		cod = 323;
		
		for (i=0; i<n; i++) {
			osnw[i] = -1;
		}
		for (i=0; i<6; i++) {
			j = (long)assignation_number(flog, cod+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) osnw[j-1] = i;
		}
		for (i=6; i<9; i++) {
			for (k=0; k<m; k++) {
				j = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (j>=1 && j<=n) osnw[j-1] = (i-6)*m + k + 6;
			}
		}
		for (i=9; i<10; i++) {
			for (k=0; k<par->max_snow_layers; k++) {
				j = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (j>=1 && j<=n) osnw[j-1] = (i-9)*par->max_snow_layers + k + 6 + 3*m;
			}
		}
		nosnw = 0;
		for (i=0; i<n; i++) {
			if(osnw[i] > 0) nosnw = i+1;
		}
	}
		
	//output glacier column
	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
		m = par->glac_plot_depths->nh;
	}else {
		m = par->max_glac_layers;
	}
	n = 6 + 3*m + 1*par->max_glac_layers;
	oglc = (long*)malloc(n*sizeof(long));
	
	par->all_glac = (short)assignation_number(flog, 348, 0, keyword, num_param, num_param_components, 0., 0);
	
	if (par->all_glac == 1) {
		
		for (i=0; i<n; i++) {
			oglc[i] = i;
		}
		
		noglc = n;
		
	}else {
		
		cod = 335;
		
		for (i=0; i<n; i++) {
			oglc[i] = -1;
		}
		for (i=0; i<6; i++) {
			j = (long)assignation_number(flog, cod+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) oglc[j-1] = i;
		}
		for (i=6; i<9; i++) {
			for (k=0; k<m; k++) {
				j = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (j>=1 && j<=n) oglc[j-1] = (i-6)*m + k + 6;
			}
		}
		for (i=9; i<10; i++) {
			for (k=0; k<par->max_glac_layers; k++) {
				j = (long)assignation_number(flog, cod+i, k, keyword, num_param, num_param_components, -1., 0);
				if (j>=1 && j<=n) oglc[j-1] = (i-9)*par->max_glac_layers + k + 6 + 3*m;
			}
		}
		noglc = 0;
		for (i=0; i<n; i++) {
			if(oglc[i] > 0) noglc = i+1;
		}
	}
	
	//output soil column
	n = 6;
	osl = (long*)malloc(n*sizeof(long));
	
	par->all_soil = (short)assignation_number(flog, 355, 0, keyword, num_param, num_param_components, 0., 0);
	
	if (par->all_soil == 1) {
		
		for (i=0; i<n; i++) {
			osl[i] = i;
		}
		
		nosl = n;
		
	}else {
		
		for (i=0; i<n; i++) {
			osl[i] = -1;
		}
		for (i=0; i<n; i++) {
			j = (long)assignation_number(flog, 349+i, 0, keyword, num_param, num_param_components, -1., 0);
			if (j>=1 && j<=n) osl[j-1] = i;
		}
		nosl = 0;
		for (i=0; i<n; i++) {
			if(osl[i] > 0) nosl = i+1;
		}
	}
	
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
	
	cod = 370;
	par->linear_interpolation_meteo = new_shortvector(nmeteo_stations);
	par->linear_interpolation_meteo->co[1] = (short)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
	for (i=2; i<=nmeteo_stations; i++) {
		par->linear_interpolation_meteo->co[i] = (short)assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->linear_interpolation_meteo->co[i-1], 0);
	}
	
	cod = 371;
	if(par->point_sim != 1){
		par->output_vertical_distances = (short)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 0., 0);
		if(par->output_vertical_distances == 1){
			printf("Only for point simulations the parameter %s can be assigned to 1, layers are defined vertically\n",keyword[cod]);
			fprintf(flog, "Only for point simulations the parameter %s can be assigned to 1, layers are defined vertically\n",keyword[cod]);
		}
	}else {
		par->output_vertical_distances = (short)assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 1., 0);
	}
	
	par->upwindblowingsnow = (short)assignation_number(flog, 372, 0, keyword, num_param, num_param_components, 0., 0);

	par->UpdateK = (short)assignation_number(flog, 373, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->ContRecovery = assignation_number(flog, 374, 0, keyword, num_param, num_param_components, 0., 0);
	par->flag1D = (short)assignation_number(flog, 375, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->k_to_ksat = assignation_number(flog, 376, 0, keyword, num_param, num_param_components, 0., 0);
	par->RunIfAnOldRunIsPresent = (short)assignation_number(flog, 377, 0, keyword, num_param, num_param_components, 1., 0);
	
	par->max_courant_land_channel = assignation_number(flog, 378, 0, keyword, num_param, num_param_components, 0.1, 0);
	par->min_dhsup_land_channel_out = assignation_number(flog, 379, 0, keyword, num_param, num_param_components, 1., 0);
	

	cod = 381;
	par->Nl_spinup = new_longvector(par->end_date->nh);
	par->Nl_spinup->co[1] = assignation_number(flog, cod, 0, keyword, num_param, num_param_components, 10000., 0);
	for (i=2; i<=par->end_date->nh; i++) {
		par->Nl_spinup->co[i] = assignation_number(flog, cod, i-1, keyword, num_param, num_param_components, par->Nl_spinup->co[i-1], 0);
	}
	if(par->Nl_spinup->co[1]<10000. && par->point_sim!=1){
		printf("You can use %s only if %s is set to 1\n",keyword[cod],keyword[12]);
		fprintf(flog,"You can use %s only if %s is set to 1\n",keyword[cod],keyword[12]);
		t_error("Not possible to continue");
	}

	par->newperiodinit = (short)assignation_number(flog, 382, 0, keyword, num_param, num_param_components, 0., 0);
	if(par->newperiodinit != 0 && par->point_sim != 1){
		printf("You can use %s only if %s is set to 1\n",keyword[382],keyword[12]);
		fprintf(flog,"You can use %s only if %s is set to 1\n",keyword[382],keyword[12]);
		t_error("Not possible to continue");
	}
	
	par->k1 = assignation_number(flog, 383, 0, keyword, num_param, num_param_components, 0.484, 0);
	par->k2 = assignation_number(flog, 384, 0, keyword, num_param, num_param_components, 8., 0);
	par->Lozone = assignation_number(flog, 385, 0, keyword, num_param, num_param_components, 0.3, 0);
	par->alpha_iqbal = assignation_number(flog, 386, 0, keyword, num_param, num_param_components, 1.3, 0);
	par->beta_iqbal = assignation_number(flog, 387, 0, keyword, num_param, num_param_components, 0.1, 0);
		
	par->albedoSWin = (short)assignation_number(flog, 388, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->micro = (short)assignation_number(flog, 389, 0, keyword, num_param, num_param_components, 1., 0);
	par->EB = assignation_number(flog, 390, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);
	par->Cair = assignation_number(flog, 391, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);
	par->Tsup = assignation_number(flog, 392, 0, keyword, num_param, num_param_components, (double)number_novalue, 0);

	par->Tair_default = assignation_number(flog, 393, 0, keyword, num_param, num_param_components, 5., 0);
	par->RH_default = assignation_number(flog, 394, 0, keyword, num_param, num_param_components, 70., 0)/100.;
	par->V_default = assignation_number(flog, 395, 0, keyword, num_param, num_param_components, par->Vmin, 0);
	par->Vdir_default = assignation_number(flog, 396, 0, keyword, num_param, num_param_components, 0., 0);
	par->IPrec_default = assignation_number(flog, 397, 0, keyword, num_param, num_param_components, 0., 0);
	
	par->soil_type_bedr_default = (long)assignation_number(flog, 399, 0, keyword, num_param, num_param_components, 1., 0);
	if(par->soil_type_bedr_default<1 || par->soil_type_bedr_default>par->nsoiltypes){
		fprintf(flog, "Error:  soil_type_bedr_default lower than 0 or higher than soil types numbers");
		printf("Error:  soil_type_bedr_default lower than 0 or higher than soil types numbers");
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}

}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

short read_soil_parameters(char *name, INIT_TOOLS *IT, SOIL *sl, long bed, FILE *flog){

	short ok;
	long i, j, k, n, nlines, nlinesprev;
	char *temp;
	double **soildata;
	DOUBLETENSOR *old_sl_par;
	FILE *f;

	//look if there is at least 1 soil file
	i = 0;
	ok = 0;
	nlinesprev = -1;

	do{

		temp = namefile_i_we2(name, i+1);

		if(existing_file_wext(temp, textfile)==1){
			free(temp);
			ok = 1;
			temp = namefile_i(name, i+1);
			nlines = count_lines(temp, 33, 44);
			if (nlinesprev >= 0 && nlines != nlinesprev){
				f = fopen(FailedRunFile, "w");
				fprintf(f,"Error:: The file %s with soil paramaters has a number of layers %ld, which different from the numbers %ld of the other soil parameter files\n",temp,nlines,nlinesprev);
				fprintf(f,"In GEOtop it is only possible to have the same number of layers in any soil parameter files\n");
				fclose(f);
				t_error("Fatal Error! Geotop is closed. See failing report.");
			}
			nlinesprev = nlines;
		}else {
			if (i==0 && strcmp(name, string_novalue) != 0) {
				f = fopen(FailedRunFile, "w");
				fprintf(f,"Error:: Soil file %s not existing.\n",name);
				fclose(f);
				t_error("Fatal Error! Geotop is closed. See failing report.");
			}
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

			if(existing_file_wext(temp, textfile)==1){
				free(temp);
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
						f = fopen(FailedRunFile, "w");
						fprintf(f,"Error:: For soil type %ld it has been given a set of soil layer thicknesses different from the other ones.\n",i);
						fprintf(f,"In GEOtop it is only possible to have the soil layer discretization in any soil parameter files.\n");
						fclose(f);
						t_error("Fatal Error! Geotop is closed. See failing report.");
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

			//field capacity and wilting point
			for (j=1; j<=sl->pa->nch; j++){
				if( (long)sl->pa->co[i][jfc][j] == number_novalue){
					sl->pa->co[i][jfc][j] = teta_psi( (-1./3.)*1.E5/g, 0., sl->pa->co[i][jsat][j], sl->pa->co[i][jres][j], sl->pa->co[i][ja][j],
													 sl->pa->co[i][jns][j], 1.-1./sl->pa->co[i][jns][j], PsiMin, sl->pa->co[i][jss][j]);
				}

				if( (long)sl->pa->co[i][jwp][j] == number_novalue){
					sl->pa->co[i][jwp][j] = teta_psi( -15.*1.E5/g, 0., sl->pa->co[i][jsat][j], sl->pa->co[i][jres][j], sl->pa->co[i][ja][j],
													 sl->pa->co[i][jns][j], 1.-1./sl->pa->co[i][jns][j], PsiMin, sl->pa->co[i][jss][j]);
				}
			}

			//pressure
			ok = 1;
			for (j=1; j<=sl->pa->nch; j++){
				if( (long)sl->pa->co[i][jpsi][j] == number_novalue) ok = 0;
			}
			if (ok == 1) IT->init_water_table_depth->co[i] = (double)number_novalue;
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
	old_sl_par = new_doubletensor(1, IT->pa_bed->nrh, IT->pa_bed->nch);
	for (n=1; n<=IT->pa_bed->nrh; n++) {
		for (j=1; j<=IT->pa_bed->nch; j++) {
			old_sl_par->co[1][n][j] = IT->pa_bed->co[1][n][j];
		}
	}
	free_doubletensor(IT->pa_bed);

	IT->pa_bed = new_doubletensor(sl->pa->ndh, sl->pa->nrh, sl->pa->nch);
	for (i=1; i<=IT->pa_bed->ndh; i++) {
		for (n=1; n<=IT->pa_bed->nrh; n++) {
			if (n == jdz) {
				for (j=1; j<=IT->pa_bed->nch; j++) {
					IT->pa_bed->co[i][n][j] = sl->pa->co[1][n][j];
				}
			}else {
				for (j=1; j<=IT->pa_bed->nch; j++) {
					if (j <= old_sl_par->nch) {
						IT->pa_bed->co[i][n][j] = old_sl_par->co[1][n][j];
					}else {
						IT->pa_bed->co[i][n][j] = IT->pa_bed->co[i][n][j-1];
					}
				}
				for (j=1; j<=IT->pa_bed->nch; j++) {
					if ( (long)IT->pa_bed->co[i][n][j] == number_novalue ) IT->pa_bed->co[i][n][j] = sl->pa->co[bed][n][j];
				}
			}
		}
	}
	free_doubletensor(old_sl_par);

	fprintf(flog,"\n");
	k = (long)nmet;
	fprintf(flog,"Soil Bedrock Layers: %ld\n",sl->pa->nch);
	for (i=1; i<=IT->pa_bed->ndh; i++) {
		fprintf(flog,"-> Soil Type: %ld\n",i);
		for (n=1; n<=nsoilprop; n++) {
			fprintf(flog,"%s: ",keywords_char[k+n-1]);
			for (j=1; j<=sl->pa->nch; j++) {
				fprintf(flog,"%f(%.2e)",IT->pa_bed->co[i][n][j],IT->pa_bed->co[i][n][j]);
				if(j<sl->pa->nch)fprintf(flog,", ");
			}
			fprintf(flog,"\n");
		}
	}
	fprintf(flog,"\n");

	return 1;

}
