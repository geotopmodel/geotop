
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
#include "../libraries/geomorphology/geomorphology.0875.h"
#include "../libraries/geomorphology/geomorphology.h"
#include "pedo.funct.h"
#include "../libraries/geomorphology/networks.h"
#include "constants.h"
#include "../libraries/geomorphology/dtm_resolution.h"
#include "../libraries/ascii/rw_maps.h"
#include "../libraries/ascii/extensions.h"
#include "../libraries/ascii/tabs.h"
#include "snow.h"
#include "meteodistr.h"
#include "vegetation.h"
#include "output.h"
#include "meteodistr.h"
#include "times.h"
#include "clouds.h"
#include "meteo.h"
#include "meteodata.h"
#include "channels.h"
#include "indices.h"
#include "recovering.h"


#include "../gt_utilities/gt_utilities.h"

extern long number_novalue, number_absent;
extern char *string_novalue;

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern char **files, *logfile;
extern long Nl, Nr, Nc;
extern long *outputpoint, noutputpoint, *outputbasin, noutputbasin, *outputsnow, noutputsnow;
extern long *outputglac, noutputglac, *outputsoil, noutputsoil;
extern char **headerpoint, **headerbasin, **headersnow, **headerglac, **headersoil;
extern char *keywords_num[num_par_number] , *keywords_char[num_par_char];


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//! Subroutine which reads input data, performs  geomporphological analisys and allocates data
void get_all_input(long argc, char *argv[], TOPO *top, SOIL *sl, LANDCOVER *land, METEO *met, WATER *wat, CHANNEL *cnet, 
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times)

{
	
	FILE *flog;
	DOUBLEMATRIX *M;
	INIT_TOOLS *IT;

	short a, success, added_JDfrom0, added_wind, added_cloud, added_Tdew;
	long l, r, c, i, ist, j, n, sy, num_cols, num_lines, day, month, year, hour, minute;
	double z, th_oversat, JD, k_snowred, maxSWE, ave;
	char *temp;			
	
	IT=(INIT_TOOLS *)malloc(sizeof(INIT_TOOLS));

	if(!argv[1]){

		WORKING_DIRECTORY=get_workingdirectory();
	} else if (argc==2) {
    // modified by Emanuele Cordano on Aug 2011
		WORKING_DIRECTORY=assign_string(argv[1]);
	}	else {
	// modified by Emanuele Cordano on Aug 2011
		WORKING_DIRECTORY=assign_string(read_option_string(argc,argv,"-wpath",".",0)); // assign_string(argv[1]); // MODIFY HERE EC

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
	printf("GEOtop 1.145 'Montebello' - 8 Nov 2010\n\n");	 
	printf("Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch \n\n");	 
	printf("GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
	printf("WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
 	printf("Riccardo Rigon is acknowledged as he coded the Fluidturtle routines (GPL Licenced), which are used in GEOtop 1.145 'Montebello'.\n");
	printf("Riccardo Rigon is also acknowledged as first founder of the GEOtop model in 1997.\n");
	printf("Riccardo Rigon and his research group are acknowledged as GEOtop 1.145 'Montebello' uses most of their modelling achievements.\n");
	printf("John Pomeroy is acknowledged as he freely provided the Prairie Blowing Snow Model Code.\n");
	printf("Glen Liston and Kelly Elder are acknowledged as they freely provided their Micromet code, from which the routines that distribute wind-air temperature-relative humidity-precipitation in GEOtop 1.145 'Montebello' are derived.\n");
	printf("However, the routine that distributes the meteo data in this GEOtop version is named Meteodistr and it significantly differs from Micromet.\n\n");
	printf("If you have satisfactorily used the code, please acknowledge the authors.\n");	
	printf("\nWORKING DIRECTORY: %s\n",WORKING_DIRECTORY);
	printf("\nLOGFILE: %s\n",logfile);
	
	fprintf(flog,"STATEMENT:\n");
	fprintf(flog,"\n");	 
	fprintf(flog,"GEOtop 1.145 'Montebello' - 8 Nov 2010\n\n");	 
	fprintf(flog,"Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch \n\n");	 
	fprintf(flog,"GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
	fprintf(flog,"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
	fprintf(flog,"Riccardo Rigon is acknowledged as he coded the Fluidturtle routines (GPL Licenced), which are used in GEOtop 1.145 'Montebello'.\n");
	fprintf(flog,"Riccardo Rigon is also acknowledged as first founder of the GEOtop model in 1997.\n");
	fprintf(flog,"Riccardo Rigon and his research group are acknowledged as GEOtop 1.145 'Montebello' uses most of their modelling achievements.\n");
	fprintf(flog,"John Pomeroy is acknowledged as he freely provided the Prairie Blowing Snow Model Code.\n");
	fprintf(flog,"Glen Liston and Kelly Elder are acknowledged as they freely provided their Micromet code, from which the routines that distribute wind-air temperature-relative humidity-precipitation in GEOtop 1.145 'Montebello' are derived.\n");
	fprintf(flog,"However, the routine that distributes the meteo data in this GEOtop version is named Meteodistr and it significantly differs from Micromet.\n\n");
	fprintf(flog,"If you have satisfactorily used the code, please acknowledge the authors.\n");	 	
	fprintf(flog,"\nWORKING DIRECTORY: %s\n",WORKING_DIRECTORY);
	
	//reads the parameters in __control_parameters
	temp = join_strings(WORKING_DIRECTORY, program_name);
	success = read_inpts_par(par, land, times, sl, met, IT, temp, flog); 
	free(temp);
	success = read_soil_parameters(files[fspar], IT->soil_col_names, sl, flog);
	
	Nl=sl->pa->nch;	
	
	success = read_point_file(files[fpointlist], IT->point_col_names, par, flog);
	
	//Time indices	
	par->delay_day_recover = 0.0;
	if(par->recover > 0){
		if(par->saving_points->nh < par->recover) t_error("Error: recover index higher than the length of the saving points vector");
		par->delay_day_recover = par->saving_points->co[par->recover];
	}
	par->init_date->co[1] += par->delay_day_recover;	
	
	convert_JDfrom0_JDandYear(par->init_date->co[1], &JD, &year);
	convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);
		
	/****************************************************************************************************/
	/*! Reading of the Input files:                                                                     */
	/****************************************************************************************************/
	
	if(par->point_sim!=1){  //distributed simulation
		read_inputmaps(top, land, sl, par, flog);
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
				par->total_area += (UV->U->co[1]*UV->U->co[2])/cos(top->slope->co[r][c]*Pi/180.);
			}
		}
	}
	
	
	top->Z = find_Z_of_any_layer(top->Z0, top->slope, land->LC, sl);
	
	top->Rdown = new_longmatrix(Nr, Nc);
	top->Cdown = new_longmatrix(Nr, Nc);	
	
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
	if (par->usemeteoio==1){
		
	}else{// classic GEOtop parsing
		for(i=1;i<=met->st->E->nh;i++){// for all the meteo stations
			
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
					printf("Error:: Date Column missing in file %s\n",temp);
					stop_execution();
					t_error("Not Possible To Continue (1) ");
				}
				met->numlines[i-1] = num_lines;

				//fixing dates: converting times in the same standard time set for the simulation and fill JDfrom0
				added_JDfrom0 = fixing_dates(ist, met->data[i-1], par->ST, met->st->ST->co[i], met->numlines[i-1], iDate12, iJDfrom0);

				//calcululate Wx and Wy if Wspeed and direction are given
				added_wind = fill_wind_speed(met->data[i-1], met->numlines[i-1], iWs, iWdir, iWsx, iWsy, IT->met_col_names[iWsx], IT->met_col_names[iWsy]);

				//find clouds
				added_cloud = fill_meteo_data_with_cloudiness(met->data[i-1], met->numlines[i-1], met->horizon[i-1], met->horizonlines[i-1],
												met->st->lat->co[i], met->st->lon->co[i], par->ST, met->st->Z->co[i], met->st->sky->co[i], 0.0);

				//find Tdew
				added_Tdew = fill_Tdew(i, met->st->Z, met->data[i-1], met->numlines[i-1], iRh, iT, iTdew, IT->met_col_names[iTdew], par->RHmin);

				//rewrite completed files
				rewrite_meteo_files(met->data[i-1], met->numlines[i-1], IT->met_col_names, temp, added_JDfrom0, added_wind, added_cloud, added_Tdew);

				free(temp);

			}else {

				t_error("Error: File meteo not in the list, meteo data not read, not possible to continue\n");

			}
		}

	}
	//read LAPSE RATES FILE  
	
	if(strcmp(files[fLRs] , string_novalue) != 0){   //s stands for string
		
		if(existing_file_text(files[fLRs])==0) printf("Lapse rate file unavailable. Check input files. If you do not have a lapse rate file, remove its name and keywords from input file\n");
		
		temp = join_strings(files[fLRs],textfile);	
		
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
	met->nstsrad=0;
	do{
		met->nstsrad++;
		a=0;
		if( (long)met->data[met->nstsrad-1][0][iSW]!=number_absent || ((long)met->data[met->nstsrad-1][0][iSWb]!=number_absent && (long)met->data[met->nstsrad-1][0][iSWd]!=number_absent ) ) a=1;
	}while(met->nstsrad<met->st->Z->nh && a==0);
	if(a==0){
		printf("WARNING: NO shortwave radiation measurements available\n");
		fprintf(flog,"WARNING: NO shortwave radiation measurements available\n");
	}else{
		printf("Shortwave radiation measurements from station %ld\n",met->nstsrad);
		fprintf(flog,"Shortwave radiation measurements from station %ld\n",met->nstsrad);
	}
	
	//FIND A STATION WITH CLOUD DATA
	met->nstcloud=0;
	do{
		met->nstcloud++;
		a=0;
		if( (long)met->data[met->nstcloud-1][0][iC]!=number_absent || (long)met->data[met->nstcloud-1][0][itauC]!=number_absent ) a=1;
	}while(met->nstcloud<met->st->Z->nh && a==0);
	if(a==0){
		printf("WARNING: NO cloudiness measurements available\n");
		fprintf(flog,"WARNING: NO cloudiness measurements available\n");
	}else{
		printf("Cloudiness measurements from station %ld\n",met->nstcloud);
		fprintf(flog,"Cloudiness measurements from station %ld\n",met->nstcloud);
	}
	
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
	
	
	/****************************************************************************************************/
	/*! Completing the several time-indipendent input variables with the data of input-files:           */
	/****************************************************************************************************/
	/****************************************************************************************************/
	// Completing of "land" (of the type LANDCOVER):  
	
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
			
			temp = namefile_i_we(files[fvegpar], i);
			
			if(existing_file_text(temp)==1){ 
				
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
				printf("hc:%f m, zmu:%f m, zmt:%f m - set hc lower than measurement height - land cover %ld, meteo station %ld\n",
					   0.001*land->ty->co[i][jHveg],met->st->Vheight->co[l],met->st->Theight->co[l],i,l);
				t_error("ERROR 1");
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
	cnet->Vsub=new_doublevector(i);
	
	//cnet->Vsup_cum=new_doublevector(i);
	//initialize_doublevector(cnet->Vsup_cum, 0.0);

	//cnet->Vsub_cum=new_doublevector(i);
	//initialize_doublevector(cnet->Vsub_cum, 0.0);
	
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
	
	lch3_cont(cnet->r->nh, cnet->ch3, cnet->lch);
	
	
	/****************************************************************************************************/
	//Calculates distance from the main channel
	
	/*DOUBLEMATRIX *M;
	M=new_doublematrix(land->LC->nrh, land->LC->nch);
	distance_from_channel2(M, top->pixel_type, cnet->r, cnet->c);
	
	long R=Nr;	
	for (i=1; i<=cnet->r->nh; i++) {
		if (cnet->r->co[i]<=R) R=cnet->r->co[i];
	}
	long C;
	for(r=1;r<=Nr;r++){
		C=Nc;
		for (i=1; i<=cnet->r->nh; i++) {
			if (r==cnet->r->co[i]) C=cnet->c->co[i];
		}
		for(c=1;c<=Nc;c++){
			if((long)land->LC->co[r][c]!=number_novalue){
				M->co[r][c]*=UV->U->co[1];
				if(M->co[r][c]>45) M->co[r][c]=(double)number_novalue;
				if(c>=C) M->co[r][c]=(double)number_novalue;
				if(r<=R) M->co[r][c]=(double)number_novalue;
			}
		}
	}
	
	temp=join_strings(WORKING_DIRECTORY, "dist_from_channel");
	write_map(temp, 0, par->format_out, M, UV, number_novalue);
	free(temp);
	free_doublematrix(M);*/
	
	//Cont for Richards 3D
	//There are not channels
	top->i_cont=(long ***)malloc((Nl+1)*sizeof(long**));
	for(l=0;l<=Nl;l++){
		top->i_cont[l]=(long **)malloc((Nr+1)*sizeof(long*));
		for(r=1;r<=Nr;r++){
			top->i_cont[l][r]=(long *)malloc((Nc+1)*sizeof(long));
		}
	}
	
	top->lrc_cont=new_longmatrix( (Nl+1)*par->total_pixel, 3);
	initialize_longmatrix(top->lrc_cont, 0);
	
	i_lrc_cont(land->LC, top->i_cont, top->lrc_cont);
	
	top->j_cont=(long **)malloc((Nr+1)*sizeof(long*));
	for (r=1; r<=Nr; r++) {
		top->j_cont[r]=(long *)malloc((Nc+1)*sizeof(long));
	}
	
	top->rc_cont=new_longmatrix(par->total_pixel, 2);
	initialize_longmatrix(top->rc_cont, 0);
	
	j_rc_cont(land->LC, top->j_cont, top->rc_cont);
	
	//BEDROCK (adjusting soil properties)
	par->bedrock = 0;	
	if( strcmp(files[fbed] , string_novalue) != 0 ) set_bedrock(sl, cnet, par, top, land->LC, flog);	
	free_doubletensor(sl->pa_bed);	
	
	/****************************************************************************************************/
	/*! Completing of the initialization of SOIL structure                               */
	/****************************************************************************************************/
		
	sl->th=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(sl->th,0.0);
	
	sl->thice=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(sl->thice,0.0);
	
	sl->ET=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(sl->ET,0.0);
	
	sl->P=new_doubletensor0(Nl,Nr,Nc);
	initialize_doubletensor(sl->P,0.0);
	
	sl->Ptot=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(sl->Ptot,0.0);
	
	sl->T=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(sl->T,0.0);
	
	sl->T_av_tensor=new_doubletensor(Nl,Nr,Nc);
	initialize_doubletensor(sl->T_av_tensor,0.0);
	
	sl->Tv=new_doublematrix(Nr,Nc);
	initialize_doublematrix(sl->Tv,0.0);
		
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if((long)land->LC->co[r][c]!=number_novalue){
				
				for(i=1;i<=sl->pa->ndh;i++){
					
					if(sl->type->co[r][c]==i){
						
						sy=sl->type->co[r][c];
						
						//sl->Tv->co[r][c]=sl->pa->co[sy][jT][1];
						
						if ((long)sl->init_water_table_height->co[sy] != number_novalue) {
							z = 0.;
							sl->P->co[0][r][c] = sl->init_water_table_height->co[sy];
							for(l=1;l<=Nl;l++){
								z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
								sl->P->co[l][r][c] = sl->init_water_table_height->co[sy] + z;
								z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
							}
						}else {
							for(l=1;l<=Nl;l++){
								sl->P->co[l][r][c] = sl->pa->co[sy][jpsi][l];
							}
						}
						
						for(l=1;l<=Nl;l++){

							sl->Ptot->co[l][r][c] = sl->P->co[l][r][c];
							
							sl->T->co[l][r][c]=sl->pa->co[sy][jT][l];
							sl->T_av_tensor->co[l][r][c]=0.0;
							
							sl->th->co[l][r][c] = teta_psi(sl->P->co[l][r][c], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
														   sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
							
							th_oversat = Fmax( sl->P->co[l][r][c] , 0.0 ) * sl->pa->co[sy][jss][l];
							sl->th->co[l][r][c] -= th_oversat;
							
							if(sl->T->co[l][r][c]<=Tfreezing){
								//Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
								sl->thice->co[l][r][c] = sl->th->co[l][r][c] - teta_psi(Psif(sl->T->co[l][r][c]), 0.0, sl->pa->co[sy][jsat][l], 
												sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 
												1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
								
								//if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
								if(sl->thice->co[l][r][c]<0) sl->thice->co[l][r][c]=0.0;
								
								//Psi is updated taking into account the freezing
								sl->th->co[l][r][c] -= sl->thice->co[l][r][c];
								sl->P->co[l][r][c] = psi_teta(sl->th->co[l][r][c] + th_oversat, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l],
															  sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 
															  1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
							}
						}
					}
				}
				
			}else{
				
				sl->P->co[0][r][c]=(double)number_novalue;
				
				for(l=1;l<=Nl;l++){
					sl->P->co[l][r][c]=(double)number_novalue;
					sl->T->co[l][r][c]=(double)number_novalue;
					sl->T_av_tensor->co[l][r][c]=(double)number_novalue;
					sl->thice->co[l][r][c]=(double)number_novalue;
					sl->th->co[l][r][c]=(double)number_novalue;	
					sl->Ptot->co[l][r][c]=(double)number_novalue;
					sl->Tv->co[r][c]=(double)number_novalue;
				}
				
			}
		}	
	}
				
	if(par->state_pixel == 1){
		if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thicezplot = new_doublematrix(par->rc->nrh, Nl);
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thicezavplot = new_doublematrix(par->rc->nrh, Nl);
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					for(i=1;i<=par->rc->nrh;i++){		
						if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){
							for(l=1;l<=Nl;l++){
								if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot->co[i][l] = sl->T->co[l][r][c];
								if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] = sl->T->co[l][r][c];
								if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot->co[i][l] = sl->Ptot->co[l][r][c];
								if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot->co[i][l] = sl->P->co[l][r][c];
								if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot->co[i][l] = sl->th->co[l][r][c];
								if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] = sl->th->co[l][r][c];
								if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thicezplot->co[i][l] = sl->thice->co[l][r][c];
								if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thicezavplot->co[i][l] = sl->thice->co[l][r][c];
							}
						}
					}
				}
			}
		}
	}
	
	if(par->recover > 0){
		assign_recovered_map(par->recover, files[rTv], sl->Tv, par, land->LC, IT->LU);
		assign_recovered_tensor(par->recover, files[rpsi], sl->P, par, land->LC, IT->LU);
		assign_recovered_tensor(par->recover, files[riceg], sl->thice, par, land->LC, IT->LU);
		assign_recovered_tensor(par->recover, files[rTg], sl->T, par, land->LC, IT->LU);
		
		for(r=1; r<=Nr; r++){
			for(c=1; c<=Nc; c++){
				if( (long)land->LC->co[r][c]!=number_novalue){
					for (l=1; l<=Nl; l++) {
						sl->th->co[l][r][c] = theta_from_psi( Fmin ( sl->P->co[l][r][c], psisat_from(l, r, c, sl) ), l, r, c, sl, PsiMin );
					}
				}
			}
		}
	}
	
	//channel soil
	cnet->th = new_doublematrix(Nl, cnet->r->nh);
	cnet->thice = new_doublematrix(Nl, cnet->r->nh);
	cnet->T = new_doublematrix(Nl, cnet->r->nh);
	cnet->P = new_doublematrix0(Nl, cnet->r->nh);
	
	cnet->ET = new_doublematrix(Nl, cnet->r->nh);	
	initialize_doublematrix(cnet->ET, 0.0);
	
	cnet->Tgskin = new_doublevector(cnet->r->nh);
	initialize_doublevector(cnet->Tgskin, 0.0);
	
	for(j=1;j<=par->total_channel;j++){
		for(i=1;i<=sl->pa->ndh;i++){
			if(cnet->soil_type->co[j]==i){
				
				sy=cnet->soil_type->co[j];
				r=cnet->r->co[j];
				c=cnet->c->co[j];
				
				if ((long)sl->init_water_table_height->co[sy] != number_novalue) {
					z = 0.;
					cnet->P->co[0][j] = sl->init_water_table_height->co[sy];
				
					for(l=1;l<=Nl;l++){
						z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
						cnet->P->co[l][j] = sl->init_water_table_height->co[sy] + z;
						z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
					}
				}else {
					for(l=1;l<=Nl;l++){
						cnet->P->co[l][j]= sl->pa->co[sy][jpsi][l];
					}
				}

				for(l=1;l<=Nl;l++){
					
					cnet->T->co[l][j]=sl->pa->co[sy][jT][l];
					cnet->th->co[l][j] = teta_psi(cnet->P->co[l][j], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], 
												  sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], 
												  PsiMin, sl->pa->co[sy][jss][l]);
					
					th_oversat = Fmax( cnet->P->co[l][j] , 0.0 ) * sl->pa->co[sy][jss][l];
					cnet->th->co[l][j] -= th_oversat;
					
					if(cnet->T->co[l][j]<=Tfreezing){
						//Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
						cnet->thice->co[l][j] = cnet->th->co[l][j] - teta_psi(Psif(cnet->T->co[l][j]), 0.0, sl->pa->co[sy][jsat][l], 
																			  sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 
																			  1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
						
						//if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
						if(cnet->thice->co[l][j]<0) cnet->thice->co[l][j]=0.0;
						
						//Psi is updated taking into account the freezing
						cnet->th->co[l][j] -= cnet->thice->co[l][j];
						cnet->P->co[l][j] = psi_teta(cnet->th->co[l][j] + th_oversat, cnet->thice->co[l][j], 
													 sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
													 sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
					}
				}
			}
		}	
	}
	
	if(par->recover > 0 && par->total_channel > 0){
		assign_recovered_tensor_channel(par->recover, files[rpsich], cnet->P, cnet->r, cnet->c, top->Z0);
		assign_recovered_tensor_channel(par->recover, files[ricegch], cnet->thice, cnet->r, cnet->c, top->Z0);
		assign_recovered_tensor_channel(par->recover, files[rTgch], cnet->T, cnet->r, cnet->c, top->Z0);
		
		for(i=1; i<=par->total_channel; i++){
			for (l=1; l<=Nl; l++) {
				sy = cnet->soil_type->co[i];
				cnet->th->co[l][i] = teta_psi(Fmin(cnet->P->co[l][i], psi_saturation(cnet->thice->co[l][i], sl->pa->co[sy][jsat][l], 
						sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l])), 
						cnet->thice->co[l][i], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], 
						sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);				
			}
		}
		
	}
	
	//WRITE INITIAL CONDITION
	write_output_headers(met->st->Z->nh, times, wat, par, top, land, sl, egy, snow, glac);
	
	if(par->state_pixel == 1){
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]!=number_novalue){//if the pixel is not a novalue				
					for(j=1;j<=par->rc->nrh;j++){
						if(r==par->rc->co[j][1] && c==par->rc->co[j][2]) write_soil_output(j, par->IDpoint->co[j], par->init_date->co[1], par->init_date->co[1], JD, day, month, year, hour, minute, par->rc, sl, PsiMin);
					}
				}
			}
		}
	}
	
	//z boundary condition
	for(l=1;l<=Nl;l++){
		par->Zboundary -= sl->pa->co[1][jdz][l];
	}
	if(par->Zboundary < sl->pa->co[1][jdz][Nl]/2.){
		printf("Z at which 0 annual temperature takes place is not lower than the soil column\n");
		t_error("Please correct");
	}
	par->Zboundary *= 1.E-3;	//convert in [m]
	
	/****************************************************************************************************/
	/*! Initialization of the struct "egy" (of the type ENERGY):*/
	
	egy->hsun = new_doublevector((long)pow(2.,par->max_times_halving_time_step_en+1.)-1);
	egy->dsun = new_doublevector((long)pow(2.,par->max_times_halving_time_step_en+1.)-1);
	egy->sinhsun = new_doublevector((long)pow(2.,par->max_times_halving_time_step_en+1.)-1);
	
	if(par->output_surfenergy>0){
		egy->Rn_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->Rn_min=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->Rn_max=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->LW_max=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->LW_min=new_doublematrix(Nr,Nc);
		egy->LWin_mean=new_doublematrix(Nr,Nc);
		egy->LW_mean=new_doublematrix(Nr,Nc);
		egy->SW_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->SW_max=new_doublematrix(Nr,Nc);
		
		egy->ET_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->ET_max=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->ET_min=new_doublematrix(Nr,Nc);
		
		egy->H_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->H_max=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->H_min=new_doublematrix(Nr,Nc); 
		
		egy->SEB_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->G_max=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->G_min=new_doublematrix(Nr,Nc);
		egy->G_snowsoil=new_doublematrix(Nr,Nc);
		
		egy->Ts_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->Ts_max=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->Ts_min=new_doublematrix(Nr,Nc);
		
		egy->Rswdown_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)egy->Rswdown_max=new_doublematrix(Nr,Nc);
		egy->Rswbeam_mean=new_doublematrix(Nr,Nc);
	}
	
	egy->nDt_shadow=new_longmatrix(Nr,Nc);
	egy->nDt_sun=new_longmatrix(Nr,Nc); 
	initialize_longmatrix(egy->nDt_shadow,0);
	initialize_longmatrix(egy->nDt_sun,0); 
	
	egy->sun = (double*)malloc(8*sizeof(double));
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if((long)land->LC->co[r][c]!=number_novalue){
				
				if(par->output_surfenergy>0){
					egy->Rn_mean->co[r][c]=0.0;
					if(par->distr_stat==1)egy->Rn_max->co[r][c]=-9.0E99;
					if(par->distr_stat==1)egy->Rn_min->co[r][c]=+9.0E99;
					if(par->distr_stat==1)egy->LW_min->co[r][c]=+9.0E99;
					if(par->distr_stat==1)egy->LW_max->co[r][c]=-9.0E99;
					egy->LWin_mean->co[r][c]=0.0;
					egy->LW_mean->co[r][c]=0.0;
					egy->SW_mean->co[r][c]=0.0;
					if(par->distr_stat==1)egy->SW_max->co[r][c]=-9.0E99;
					
					egy->ET_mean->co[r][c]=0.0;
					if(par->distr_stat==1)egy->ET_max->co[r][c]=-9.0E99;
					if(par->distr_stat==1)egy->ET_min->co[r][c]=+9.0E99;
					
					egy->H_mean->co[r][c]=0.0;
					if(par->distr_stat==1)egy->H_max->co[r][c]=-9.0E99;
					if(par->distr_stat==1)egy->H_min->co[r][c]=+9.0E99;
					
					egy->SEB_mean->co[r][c]=0.0;
					if(par->distr_stat==1)egy->G_max->co[r][c]=-9.0E99;
					if(par->distr_stat==1)egy->G_min->co[r][c]=+9.0E99;
					egy->G_snowsoil->co[r][c]=0.0;
					
					egy->Ts_mean->co[r][c]=0.0;
					if(par->distr_stat==1)egy->Ts_max->co[r][c]=-9.0E99;
					if(par->distr_stat==1)egy->Ts_min->co[r][c]=+9.0E99;
					
					egy->Rswdown_mean->co[r][c]=0.0;
					if(par->distr_stat==1)egy->Rswdown_max->co[r][c]=-9.0E99;
					egy->Rswbeam_mean->co[r][c]=0.0;
					
				}
				
			}else{
				
				if(par->output_surfenergy>0){
					egy->Rn_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->Rn_max->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->Rn_min->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->LW_min->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->LW_max->co[r][c]=(double)number_novalue;
					egy->LWin_mean->co[r][c]=(double)number_novalue;
					egy->LW_mean->co[r][c]=(double)number_novalue;
					egy->SW_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->SW_max->co[r][c]=(double)number_novalue;			
					
					egy->ET_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->ET_max->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->ET_min->co[r][c]=(double)number_novalue;
					
					egy->H_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->H_max->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->H_min->co[r][c]=(double)number_novalue;
					
					egy->SEB_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->G_max->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->G_min->co[r][c]=(double)number_novalue;
					egy->G_snowsoil->co[r][c]=(double)number_novalue;
					
					egy->Ts_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->Ts_max->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->Ts_min->co[r][c]=(double)number_novalue;
					
					egy->Rswdown_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)egy->Rswdown_max->co[r][c]=(double)number_novalue;
					egy->Rswbeam_mean->co[r][c]=(double)number_novalue;
					
				}
			}
		}
	}
	
	if(times->JD_plots->nh > 1){
		egy->Hgplot=new_doublematrix(Nr,Nc);
		egy->LEgplot=new_doublematrix(Nr,Nc);
		egy->Hvplot=new_doublematrix(Nr,Nc);
		egy->LEvplot=new_doublematrix(Nr,Nc);
		
		egy->SWinplot=new_doublematrix(Nr,Nc);
		egy->SWgplot=new_doublematrix(Nr,Nc);
		egy->SWvplot=new_doublematrix(Nr,Nc);
		
		egy->LWinplot=new_doublematrix(Nr,Nc);
		egy->LWgplot=new_doublematrix(Nr,Nc);
		egy->LWvplot=new_doublematrix(Nr,Nc);
		
		egy->Tsplot=new_doublematrix(Nr,Nc);
		egy->Tgplot=new_doublematrix(Nr,Nc);
		egy->Tvplot=new_doublematrix(Nr,Nc);
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]==number_novalue){
					egy->Hgplot->co[r][c]=(double)number_novalue;
					egy->LEgplot->co[r][c]=(double)number_novalue;
					egy->Hvplot->co[r][c]=(double)number_novalue;
					egy->LEvplot->co[r][c]=(double)number_novalue;
					egy->SWinplot->co[r][c]=(double)number_novalue;
					egy->SWgplot->co[r][c]=(double)number_novalue;
					egy->SWvplot->co[r][c]=(double)number_novalue;
					egy->LWinplot->co[r][c]=(double)number_novalue;
					egy->LWgplot->co[r][c]=(double)number_novalue;
					egy->LWvplot->co[r][c]=(double)number_novalue;
					egy->Tsplot->co[r][c]=(double)number_novalue;
					egy->Tgplot->co[r][c]=(double)number_novalue;
					egy->Tvplot->co[r][c]=(double)number_novalue;				
				}else{
					egy->Hgplot->co[r][c]=0.0;
					egy->LEgplot->co[r][c]=0.0;
					egy->Hvplot->co[r][c]=0.0;
					egy->LEvplot->co[r][c]=0.0;
					egy->SWinplot->co[r][c]=0.0;
					egy->SWgplot->co[r][c]=0.0;
					egy->SWvplot->co[r][c]=0.0;
					egy->LWinplot->co[r][c]=0.0;
					egy->LWgplot->co[r][c]=0.0;
					egy->LWvplot->co[r][c]=0.0;
					egy->Tsplot->co[r][c]=0.0;
					egy->Tgplot->co[r][c]=0.0;
					egy->Tvplot->co[r][c]=0.0;				
				}
			}
		}
	}
	
	//vectors used in energy_balance()
	egy->Tgskin = new_doublematrix(Nr, Nc);
	initialize_doublematrix(egy->Tgskin, 0.);
	
	egy->Tgskinsurr = new_doublematrix(Nr, Nc);
	initialize_doublematrix(egy->Tgskinsurr, 0.);
	
	egy->Asurr = new_doublematrix(Nr, Nc);
	initialize_doublematrix(egy->Asurr, 0.);
	
	egy->Dlay = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->wliq = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->wice = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->Temp = new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->deltaw = new_doublevector( Nl + par->snowlayer_max + par->glaclayer_max );
	
	egy->SWlayer = new_doublevector0( par->snowlayer_max + 1 );
	
	egy->soil_transp_layer = new_doublevector(land->root_fraction->nch);
	initialize_doublevector(egy->soil_transp_layer, 0.);
	
	egy->dFenergy = new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->udFenergy = new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max - 1);
	egy->Kth0=new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max - 1 );
	egy->Kth1=new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max - 1);
	egy->Fenergy=new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->Newton_dir=new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->T0=new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max );
	egy->T1=new_doublevector0( Nl + par->snowlayer_max + par->glaclayer_max );
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
	/*! Completing of the struct "water" (of the type WATER) with the initializations of the remanent
	 matrices (wat->rain, wat->Pnet, wat->wcan_rain):        */
	
	/* Initialization of wat->Pnet (liquid precipitation that reaches the sl surface in mm):*/
	wat->Pnet=new_doublematrix(Nr,Nc);
	initialize_doublematrix(wat->Pnet,0.0);
	
	/* Initialization of wat->wcan_rain: (liquid precipitation intercepted by vegetation in mm):*/
	wat->wcan_rain=new_doublematrix(Nr,Nc);
	wat->wcan_snow=new_doublematrix(Nr,Nc);
	
	/* Initialization of wat->PrecTot (total precipitation (rain+snow) precipitation):*/
	wat->PrecTot=new_doubletensor((long)pow(2.,par->max_times_halving_time_step_en+1.)-1,Nr,Nc);
	initialize_doubletensor(wat->PrecTot,0.0);
	
	/* Initialization of the matrices with the output of total precipitation and interception:*/
	if(par->output_meteo>0){
		wat->PrTOT_mean=new_doublematrix(Nr,Nc);
		wat->PrSNW_mean=new_doublematrix(Nr,Nc);
	}
	
	wat->h_sup=new_doublematrix(Nr, Nc);
	
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if((long)land->LC->co[r][c]!=number_novalue){
				wat->wcan_rain->co[r][c]=0.0;
				wat->wcan_snow->co[r][c]=0.0;
				if(par->output_meteo>0){
					wat->PrTOT_mean->co[r][c]=0.0;
					wat->PrSNW_mean->co[r][c]=0.0;
				}
			}else{
				wat->wcan_rain->co[r][c]=(double)number_novalue;
				wat->wcan_snow->co[r][c]=(double)number_novalue;
				if(par->output_meteo>0){
					wat->PrTOT_mean->co[r][c]=(double)number_novalue;
					wat->PrSNW_mean->co[r][c]=(double)number_novalue;
				}
			}
		}
	}
	
	if(par->recover > 0){
		assign_recovered_map(par->recover, files[rwcrn], wat->wcan_rain, par, land->LC, IT->LU);
		assign_recovered_map(par->recover, files[rwcsn], wat->wcan_rain, par, land->LC, IT->LU);
	}
	
	
	
	/****************************************************************************************************/
	/*! Initialization of the struct "snow" (of the type SNOW):*/
	
	/***************************************************************************************************/
	/*! Optional reading of initial real snow thickness map in the whole basin ("SNOW0"):    */
	if( strcmp(files[fsn0] , string_novalue) != 0 ){
		printf("Snow initial condition from file %s\n",files[fsn0]+1);
		M=read_map(2, files[fsn0], land->LC, UV, (double)number_novalue);
	}else{
		M=copydoublematrix_const(IT->snow0*rho_w/IT->rhosnow0, land->LC, (double)number_novalue);
	}
	
	/*! Optional reading of snow age in the whole basin  */
	if( strcmp(files[fsnag0] , string_novalue) != 0 ){
		printf("Snow age initial condition from file %s\n",files[fsnag0]+1);
		snow->dimens_age=read_map(2, files[fsnag0], land->LC, UV, (double)number_novalue); 
	}else{
		snow->dimens_age=copydoublematrix_const(IT->agesnow0, land->LC, (double)number_novalue);
	}
	
	if(times->JD_plots->nh > 1){
		snow->Dplot=new_doublematrix(Nr,Nc);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]==number_novalue){
					snow->Dplot->co[r][c]=(double)number_novalue;
				}else{
					snow->Dplot->co[r][c]=0.0;
				}
			}
		}
	}
	
	snow->S=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
	allocate_and_initialize_statevar_3D(snow->S, (double)number_novalue, par->snowlayer_max, Nr, Nc);
	
	snow->nondimens_age=new_doublematrix(Nr,Nc);
	
	if(par->blowing_snow==1){
		
		snow->S_for_BS=(STATEVAR_1D *)malloc(sizeof(STATEVAR_1D));
		allocate_and_initialize_statevar_1D(snow->S_for_BS, (double)number_novalue, Nl);
				
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
				
		if(par->output_snow>0){
			snow->Wtrans_plot=new_doublematrix(Nr,Nc);
			snow->Wsubl_plot=new_doublematrix(Nr,Nc);
			/*snow->Qtrans_eq_plot=new_doublematrix(Nr,Nc);
			 snow->Qtrans_plot=new_doublematrix(Nr,Nc);
			 snow->Qsub_eq_plot=new_doublematrix(Nr,Nc);
			 snow->Qsub_plot=new_doublematrix(Nr,Nc);*/
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if((long)land->LC->co[r][c]==number_novalue){
						snow->Wtrans_plot->co[r][c]=(double)number_novalue;		
						snow->Wsubl_plot->co[r][c]=(double)number_novalue;	
						/*snow->Qtrans_eq_plot->co[r][c]=(double)number_novalue;
						 snow->Qtrans_plot->co[r][c]=(double)number_novalue;
						 snow->Qsub_eq_plot->co[r][c]=(double)number_novalue;
						 snow->Qsub_plot->co[r][c]=(double)number_novalue;*/
						
					}else{
						snow->Wtrans_plot->co[r][c]=0.0;					
						snow->Wsubl_plot->co[r][c]=0.0;			
						/*snow->Qtrans_eq_plot->co[r][c]=0.0;
						 snow->Qtrans_plot->co[r][c]=0.0;
						 snow->Qsub_eq_plot->co[r][c]=0.0;
						 snow->Qsub_plot->co[r][c]=0.0;*/
						
					}
				}
			}
		}
	}
	
	if(par->output_snow>0){
		snow->MELTED=new_doublematrix(Nr,Nc);
		initialize_doublematrix(snow->MELTED,number_novalue);
		snow->SUBL=new_doublematrix(Nr,Nc);
		initialize_doublematrix(snow->SUBL,number_novalue);
		snow->t_snow=new_doublematrix(Nr,Nc);
		initialize_doublematrix(snow->t_snow,number_novalue);
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
					M->co[r][c] *= k_snowred;
				}
				
				
				if(M->co[r][c]<0){
					t_error("Error: negative snow data");
				}else if(M->co[r][c]>1.E-5){
					snow->dimens_age->co[r][c]*=86400.0;	//now in [s]
					
					snow->S->lnum->co[r][c]=1;
					snow->S->Dzl->co[1][r][c]=M->co[r][c];
					snow->S->w_ice->co[1][r][c]=IT->rhosnow0*0.001*snow->S->Dzl->co[1][r][c];
					snow->S->T->co[1][r][c]=IT->Tsnow0;
				}
												
				snow->nondimens_age->co[r][c]=snow->dimens_age->co[r][c];
				non_dimensionalize_snowage(&(snow->nondimens_age->co[r][c]), IT->Tsnow0);		
				
				if (par->point_sim == 1) {
					maxSWE = par->maxSWE->co[r][c];
				}else {
					maxSWE = 1.E10;
				}
				snow_layer_combination(par->alpha_snow, r, c, snow->S, -1., par->snowlayer_inf, par->Dmin, par->Dmax, maxSWE, 0., flog);	
				
				if(par->output_snow>0){
					//snow->max->co[r][c]=0.0;
					//snow->average->co[r][c]=0.0;
					snow->MELTED->co[r][c]=0.0;
					snow->SUBL->co[r][c]=0.0;
					snow->t_snow->co[r][c]=0.0;
				}
				
			}
		}
	}
	
	free_doublematrix(M);
		
	if(par->recover > 0){
		initialize_shortmatrix(snow->S->type, 2);
		assign_recovered_map_long(par->recover, files[rns], snow->S->lnum, par, land->LC, IT->LU);
		assign_recovered_map(par->recover, files[rsnag_nondim], snow->nondimens_age, par, land->LC, IT->LU);
		assign_recovered_map(par->recover, files[rsnag_dim], snow->dimens_age, par, land->LC, IT->LU);
		assign_recovered_tensor(par->recover, files[rDzs], snow->S->Dzl, par, land->LC, IT->LU);
		assign_recovered_tensor(par->recover, files[rwls], snow->S->w_liq, par, land->LC, IT->LU);
		assign_recovered_tensor(par->recover, files[rwis], snow->S->w_ice, par, land->LC, IT->LU);
		assign_recovered_tensor(par->recover, files[rTs], snow->S->T, par, land->LC, IT->LU);
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if( (long)land->LC->co[r][c]!=number_novalue){
					snow_layer_combination(par->alpha_snow, r, c, snow->S, -1., par->snowlayer_inf, par->Dmin, par->Dmax, maxSWE, 0., flog);	
				}
			}
		}
				
	}
	
	
	/****************************************************************************************************/
	/*! Initialization of the struct "glac" (of the type GLACIER):*/
	
	/***************************************************************************************************/ 
	/*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
	if( par->point_sim!=1 && strcmp(files[fgl0] , string_novalue) != 0 && par->glaclayer_max==0){
		printf("Warning: Glacier map present, but glacier represented with 0 layers\n");
		fprintf(flog,"Warning: Glacier map present, but glacier represented with 0 layers\n");
	}
	
	if(par->glaclayer_max==0 && IT->Dglac0>0){
		printf("\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
		fprintf(flog,"\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
	}
	
	//If the max number of glacier layers is greater than 1, the matrices (or tensors) lnum, Dzl. w_liq, w_ice, T and print matrices are defined, according to the respective flags
	if(par->glaclayer_max>0){
		
		if( par->point_sim!=1 && strcmp(files[fgl0] , string_novalue) != 0 ){
			printf("Glacier initial condition from file %s\n",files[fgl0]+1);
			fprintf(flog,"Glacier initial condition from file %s\n",files[fgl0]+1);
			M=read_map(2, files[fgl0], land->LC, UV, (double)number_novalue);
		}else{
			M=copydoublematrix_const(IT->Dglac0, land->LC, (double)number_novalue);
		}		
		
		glac->G=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
		allocate_and_initialize_statevar_3D(glac->G, (double)number_novalue, par->glaclayer_max, Nr, Nc);
		
		if(par->output_glac>0){
			glac->MELTED=new_doublematrix(Nr,Nc);
			initialize_doublematrix(glac->MELTED,number_novalue);
			glac->SUBL=new_doublematrix(Nr,Nc);
			initialize_doublematrix(glac->SUBL,number_novalue);
		}
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if( (long)land->LC->co[r][c]!=number_novalue){
					
					if(M->co[r][c]<0){
						t_error("Error: negative glacier data");
					}else if(M->co[r][c]>1.E-5){
						glac->G->lnum->co[r][c]=1;
						glac->G->Dzl->co[1][r][c]=M->co[r][c];
						glac->G->w_ice->co[1][r][c]=IT->rhoglac0*0.001*glac->G->Dzl->co[1][r][c];
						glac->G->T->co[1][r][c]=IT->Tglac0;
					}
					
					if(par->output_glac>0){
						glac->MELTED->co[r][c]=0.0;
						glac->SUBL->co[r][c]=0.0;
					}
										
					snow_layer_combination(par->alpha_snow, r, c, glac->G, -1., par->glaclayer_inf, par->Dmin_glac, par->Dmax_glac, 1.E10, 0., flog);	
					
				}
			}
		}
		
		if(par->recover > 0){
			assign_recovered_map_long(par->recover, files[rni], glac->G->lnum, par, land->LC, IT->LU);
			assign_recovered_tensor(par->recover, files[rDzi], glac->G->Dzl, par, land->LC, IT->LU);
			assign_recovered_tensor(par->recover, files[rwli], glac->G->w_liq, par, land->LC, IT->LU);
			assign_recovered_tensor(par->recover, files[rwii], glac->G->w_ice, par, land->LC, IT->LU);
			assign_recovered_tensor(par->recover, files[rTi], glac->G->T, par, land->LC, IT->LU);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						snow_layer_combination(par->alpha_snow, r, c, glac->G, -1., par->glaclayer_inf, par->Dmin_glac, par->Dmax_glac, 1.E10, 0., flog);
					}
				}
			}					
		}
	}																					
	
	
	
	//***************************************************************************************************
	// Filling up of the struct "met" (of the type METEO):
	
	met->tau_cloud = new_doublevector((long)pow(2.,par->max_times_halving_time_step_en+1.)-1);
	met->tau_cloud_av = new_doublevector((long)pow(2.,par->max_times_halving_time_step_en+1.)-1);
	met->tau_cloud_yes = new_shortvector((long)pow(2.,par->max_times_halving_time_step_en+1.)-1);
	met->tau_cloud_av_yes = new_shortvector((long)pow(2.,par->max_times_halving_time_step_en+1.)-1); 
	
	met->Tgrid=new_doubletensor((long)pow(2.,par->max_times_halving_time_step_en+1.)-1,Nr,Nc);
	initialize_doubletensor(met->Tgrid, 5.);
	met->Pgrid=new_doubletensor((long)pow(2.,par->max_times_halving_time_step_en+1.)-1,Nr,Nc);
	initialize_doubletensor(met->Pgrid, Pa0);
	met->RHgrid=new_doubletensor((long)pow(2.,par->max_times_halving_time_step_en+1.)-1,Nr,Nc);
	initialize_doubletensor(met->RHgrid, 0.7);
	
	met->Vgrid=new_doubletensor((long)pow(2.,par->max_times_halving_time_step_en+1.)-1,Nr,Nc);
	initialize_doubletensor(met->Vgrid, par->Vmin);
	met->Vdir=new_doubletensor((long)pow(2.,par->max_times_halving_time_step_en+1.)-1,Nr,Nc);
	initialize_doubletensor(met->Vdir, 0.0);
		
	if(par->output_meteo>0){
		met->Ta_mean=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)met->Ta_max=new_doublematrix(Nr,Nc);
		if(par->distr_stat==1)met->Ta_min=new_doublematrix(Nr,Nc);

		met->Vspdmean=new_doublematrix(Nr,Nc);
		met->Vdirmean=new_doublematrix(Nr,Nc);
		met->RHmean=new_doublematrix(Nr,Nc);
		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if( (long)land->LC->co[r][c]!=number_novalue){
					met->Ta_mean->co[r][c]=0.0;
					if(par->distr_stat==1)met->Ta_max->co[r][c]=-9.0E99;
					if(par->distr_stat==1)met->Ta_min->co[r][c]=+9.0E99;
					met->Vspdmean->co[r][c]=0.0;
					met->Vdirmean->co[r][c]=0.0;
					met->RHmean->co[r][c]=0.0;		
					
				}else{
					met->Ta_mean->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)met->Ta_max->co[r][c]=(double)number_novalue;
					if(par->distr_stat==1)met->Ta_min->co[r][c]=(double)number_novalue;
					met->Vspdmean->co[r][c]=(double)number_novalue;
					met->Vdirmean->co[r][c]=(double)number_novalue;
					met->RHmean->co[r][c]=(double)number_novalue;
					
				}
			}
		}	
	}	
	
	
	//plot output
	if(times->JD_plots->nh > 1){
		met->Taplot=new_doublematrix(Nr,Nc);
		met->Vspdplot=new_doublematrix(Nr,Nc);
		met->Vdirplot=new_doublematrix(Nr,Nc);
		met->RHplot=new_doublematrix(Nr,Nc);		
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]==number_novalue){
					met->Taplot->co[r][c]=(double)number_novalue;
					met->Vspdplot->co[r][c]=(double)number_novalue;
					met->Vdirplot->co[r][c]=(double)number_novalue;
					met->RHplot->co[r][c]=(double)number_novalue;					
				}else{
					met->Taplot->co[r][c]=0.0;
					met->Vspdplot->co[r][c]=0.0;
					met->Vdirplot->co[r][c]=0.0;
					met->RHplot->co[r][c]=0.0;					
				}
			}
		}
	}
	
	//UPDATE Tskin
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if((long)land->LC->co[r][c]!=number_novalue){
				if (snow->S->lnum->co[r][c] > 0) {
					egy->Tgskin->co[r][c] = snow->S->T->co[1][r][c];
					if (cnet->ch->co[r][c] > 0) cnet->Tgskin->co[cnet->ch->co[r][c]] = egy->Tgskin->co[r][c];
				}else if (par->glaclayer_max > 0) {
					if (glac->G->lnum->co[r][c] > 0) {
						egy->Tgskin->co[r][c] = glac->G->T->co[1][r][c];
						if (cnet->ch->co[r][c] > 0) cnet->Tgskin->co[cnet->ch->co[r][c]] = egy->Tgskin->co[r][c];
					}else {
						egy->Tgskin->co[r][c] = sl->T->co[1][r][c];
						if (cnet->ch->co[r][c] > 0) cnet->Tgskin->co[cnet->ch->co[r][c]] = cnet->T->co[1][cnet->ch->co[r][c]];
					}
				}else {
					egy->Tgskin->co[r][c] = sl->T->co[1][r][c];
					if (cnet->ch->co[r][c] > 0) cnet->Tgskin->co[cnet->ch->co[r][c]] = cnet->T->co[1][cnet->ch->co[r][c]];
				}
			}
		}
	}
	
	//FIND SURROUDINGS ALBEDO AND SKIN TEMPERATURE
	if (par->point_sim == 1) {
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					egy->Tgskinsurr->co[r][c] = egy->Tgskin->co[r][c];
					egy->Asurr->co[r][c] = 0.;
				}
			}
		}
	}else {
		ave = 0.;
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					ave += egy->Tgskin->co[r][c] / (double)par->total_pixel;
				}
			}
		}
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					egy->Tgskinsurr->co[r][c] = ave;
					egy->Asurr->co[r][c] = 0.;
				}
			}
		}
	}

	/****************************************************************************************************/			
	
	/*Free the struct allocated in this subroutine:*/
	
	free_doublematrix(par->chkpt);
	free_doublevector(sl->init_water_table_height);

	if(par->point_sim==1 && par->recover>0){
		free_doublematrix(IT->LU);
		free_longvector(par->r_points);
		free_longvector(par->c_points);
	}

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
		
	cont_nonzero_values_matrix2(&i, &j, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel, par->point_sim);
	top->Li = new_longvector(i);
	top->Lp = new_longvector(j);
	wat->Lx = new_doublevector(i);	
	cont_nonzero_values_matrix3(top->Lp, top->Li, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel, par->point_sim);

	wat->P0 = new_doublevector(j);
	wat->H0 = new_doublevector(j);
	wat->H1 = new_doublevector(j);
	wat->dH = new_doublevector(j);
	wat->B = new_doublevector(j);
	wat->f = new_doublevector(j);
	wat->df = new_doublevector(j);
	
	wat->Klat = new_doublematrix(top->BC_LatDistance->nh,Nl);
	
	fclose(flog);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void read_inputmaps(TOPO *top, LANDCOVER *land, SOIL *sl, PAR *par, FILE *flog){

	long r, c, i, cont;					
	DOUBLEMATRIX *M;
	SHORTMATRIX *curv;
	short flag;
	char *temp;
	double min, max;
	
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
		t_error("It is impossible to proceed without giving the digital elevation model");
	}

	//reading LAND COVER TYPE
	flag = file_exists(flu, flog);
	if(flag == 1){
		land->LC=read_map(1, files[flu], top->Z0, UV, (double)number_novalue);
		//Check borders
		for(r=1;r<=land->LC->nrh;r++){
			for(i=1;i<=2;i++){
				land->LC->co[r][i]=(double)number_novalue;
				land->LC->co[r][land->LC->nch+1-i]=(double)number_novalue;
			}
		}
		for(c=1;c<=land->LC->nch;c++){
			for(i=1;i<=1;i++){
				land->LC->co[i][c]=(double)number_novalue;
				land->LC->co[land->LC->nrh+1-i][c]=(double)number_novalue;
			}
		}		
		for(r=1;r<=land->LC->nrh;r++){
			for(c=1;c<=land->LC->nch;c++){
				if ((long)land->LC->co[r][c] != number_novalue) {
					if ((long)land->LC->co[r][c] < 1 || (long)land->LC->co[r][c] > par->n_landuses) t_error("Error: It is not possible to assign Value < 1 or > n_landuses to the land cover type");
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
			
			if(par->state_px_coord==1){
				par->rc->co[i][1]=row(par->chkpt->co[i][ptY], top->Z0->nrh, UV, number_novalue);
				par->rc->co[i][2]=col(par->chkpt->co[i][ptX], top->Z0->nch, UV, number_novalue);
			}else{
				par->rc->co[i][1]=(long)par->chkpt->co[i][ptX];
				par->rc->co[i][2]=(long)par->chkpt->co[i][ptY];
			}
			
			if((long)land->LC->co[par->rc->co[i][1]][par->rc->co[i][2]]==number_novalue){
				printf("\nWARNING: Point #%4ld corresponds to NOVALUE pixel",i);
				stop_execution();
				t_error("Not Possible To Continue (2) ");
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
	//reading SOIL MAP
	flag = file_exists(fsoil, flog);
	if(flag == 1){
		M=read_map(2, files[fsoil], land->LC, UV, (double)number_novalue);
		sl->type=copylong_doublematrix(M);
		for(r=1;r<=land->LC->nrh;r++){
			for(c=1;c<=land->LC->nch;c++){
				if ((long)land->LC->co[r][c] != number_novalue) {
					if (sl->type->co[r][c] < 1 || sl->type->co[r][c] > par->nsoiltypes) t_error("Error: It is not possible to assign Value < 1 or > nsoiltypes to the soil type map");
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
		top->slope=read_map(2, files[fslp], top->Z0, UV, (double)number_novalue);		//reads in degrees
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
		top->aspect=read_map(2, files[fasp], top->Z0, UV, (double)number_novalue);
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
	//pixel type = 1 land pixel (it it is on the border, the border is permeable above an user-defined elevation in the saturated part (weir-wise)
	//pixel type = 10 channel pixel (if it is on the border, the border is impermeable, water is free only on the surface)
	
	flag = file_exists(fnet, flog);
	if(flag == 1){
		M=read_map(2, files[fnet], land->LC, UV, (double)number_novalue);
		top->pixel_type=copyshort_doublematrix(M);

		cont = 0;
		for(r=1;r<=top->Z0->nrh;r++){
			for(c=1;c<=top->Z0->nch;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					if(top->pixel_type->co[r][c]!=0 && top->pixel_type->co[r][c]!=10){
						t_error("Only values of 0 and 10 in the network map are admitted");
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
	
	//count the pixels having pixel_type = 1
	cont = 0;
	for(r=1;r<=top->Z0->nrh;r++){
		for(c=1;c<=top->Z0->nch;c++){
			if((long)land->LC->co[r][c]!=number_novalue){
				if (top->pixel_type->co[r][c] == 1) cont ++;
			}
		}
	}
	
	top->BC_counter = new_longmatrix(top->Z0->nrh, top->Z0->nch);
	initialize_longmatrix(top->BC_counter, 0);
	
	if (cont > 0) {
		top->BC_LatDistance = new_doublevector(cont);
		top->BC_DepthFreeSurface = new_doublevector(cont);
		cont = 0;
		for(r=1;r<=top->Z0->nrh;r++){
			for(c=1;c<=top->Z0->nch;c++){
				if((long)land->LC->co[r][c]!=number_novalue){
					if (top->pixel_type->co[r][c] == 1){
						cont ++;
						top->BC_counter->co[r][c] = cont;
						top->BC_LatDistance->co[cont] = 0.5*UV->U->co[1];	//[m]
						top->BC_DepthFreeSurface->co[cont] = 500.0;			//[mm]
					}
				}
			}
		}
	}else {
		top->BC_LatDistance = new_doublevector(1);
		initialize_doublevector(top->BC_LatDistance, (double)number_novalue);
		top->BC_DepthFreeSurface = new_doublevector(1);	
		initialize_doublevector(top->BC_DepthFreeSurface, (double)number_novalue);
	}


	
	/****************************************************************************************************/	
	
	top->is_on_border=new_shortmatrix(land->LC->nrh, land->LC->nch);
	for(r=1;r<=land->LC->nrh;r++){
		for(c=1;c<=land->LC->nch;c++){
			if ( (long)land->LC->co[r][c]!=number_novalue && par->point_sim!=1){
				top->is_on_border->co[r][c] = is_boundary(r, c, land->LC, number_novalue);
			}else{
				top->is_on_border->co[r][c] = -1;
			}
		}
	}
	
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void read_optionsfile_point(PAR *par, TOPO *top, LANDCOVER *land, SOIL *sl, TIMES *times, INIT_TOOLS *IT, FILE *flog){

	long i, r, c, num_lines;
	DOUBLEMATRIX *Q, *P, *R, *S, *T, *Z, *LU;
	SHORTMATRIX *curv;
	short read_dem, read_lu, read_soil, read_sl, read_as, read_sk, read_curv, flag, coordinates;
	char *temp;
	double min, max;
	
	//4. CALCULATE TOPOGRAPHIC PROPERTIES
	//check if there are point coordinates
	coordinates = 1;
	for(i=1;i<=par->chkpt->nrh;i++){ 
		if ( (long)par->chkpt->co[i][ptX]==number_novalue || (long)par->chkpt->co[i][ptY]==number_novalue ) coordinates = 0;
	}
	if (coordinates == 0 && par->recover>0){
		printf("Warning: Not possible to recover the simulation because at least one point has no coordinates\n");
		printf("Starting from normal initial condition\n");
		fprintf(flog,"Warning: Not possible to recover the simulation because at least one point has no coordinates\n");
		fprintf(flog,"Starting from normal initial condition\n");
		par->recover = 0;		
	}
	
	//a. read dem	
	read_dem=0;
	if(par->recover>0) read_dem=1;
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
			
			if(par->recover>0){
				printf("Warning: Not possible to recover the simulation because there is no dem\n");
				printf("Starting from normal initial condition\n");
				fprintf(flog,"Warning: Not possible to recover the simulation because there is no dem\n");
				fprintf(flog,"Starting from normal initial condition\n");
				par->recover = 0;
			}	
		}
	}
	
	
	if(read_dem==1){ 
		par->r_points=new_longvector(par->chkpt->nrh);
		par->c_points=new_longvector(par->chkpt->nrh);		
		for(i=1;i<=par->chkpt->nrh;i++){ 
			if(par->state_px_coord==0){
				par->r_points->co[i]=(long)par->chkpt->co[i][ptX];
				par->c_points->co[i]=(long)par->chkpt->co[i][ptY];
			}else{
				par->r_points->co[i]=row(par->chkpt->co[i][ptY], Z->nrh, UV, number_novalue);
				par->c_points->co[i]=col(par->chkpt->co[i][ptX], Z->nch, UV, number_novalue);
			}
			if((long)par->chkpt->co[i][ptZ]==number_novalue) par->chkpt->co[i][ptZ]=Z->co[par->r_points->co[i]][par->c_points->co[i]];
		}
	}

	//b. read land use
	read_lu=0;
	if(par->recover>0) read_lu=1;
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
				if(par->recover>0){
					printf("Warning: Not possible to recover the simulation because there is no dem\n");
					printf("Starting from normal initial condition\n");
					fprintf(flog,"Warning: Not possible to recover the simulation because there is no dem\n");
					fprintf(flog,"Starting from normal initial condition\n");
					par->recover = 0;
				}						
			}
		}
	}
	
	if(read_lu==1){
		for(i=1;i<=par->chkpt->nrh;i++){ 
			if((long)par->chkpt->co[i][ptLC]==number_novalue){
				if(par->state_px_coord==0){
					r=(long)par->chkpt->co[i][ptX];
					c=(long)par->chkpt->co[i][ptY];
				}else{				
					r=row(par->chkpt->co[i][ptY], LU->nrh, UV, number_novalue);
					c=col(par->chkpt->co[i][ptX], LU->nch, UV, number_novalue);
				}
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
				if(par->state_px_coord==0){
					r=(long)par->chkpt->co[i][ptX];
					c=(long)par->chkpt->co[i][ptY];
				}else{					
					r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
					c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				}
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
				if(par->state_px_coord==0){
					r=(long)par->chkpt->co[i][ptX];
					c=(long)par->chkpt->co[i][ptY];
				}else{					
					r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
					c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				}
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
				if(par->state_px_coord==0){
					r=(long)par->chkpt->co[i][ptX];
					c=(long)par->chkpt->co[i][ptY];
				}else{					
					r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
					c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				}
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
				if(par->state_px_coord==0){
					r=(long)par->chkpt->co[i][ptX];
					c=(long)par->chkpt->co[i][ptY];
				}else{					
					r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
					c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
				}
				par->chkpt->co[i][ptSKY]=P->co[r][c];
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
			if(par->state_px_coord==0){
				r=(long)par->chkpt->co[i][ptX];
				c=(long)par->chkpt->co[i][ptY];
			}else{					
				r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
				c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
			}
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
		if((long)par->chkpt->co[i][ptDrDIST]==number_novalue) par->chkpt->co[i][ptDrDIST]=1.E10;//[m]
		if((long)par->chkpt->co[i][ptDrDEPTH]==number_novalue) par->chkpt->co[i][ptDrDEPTH]=500.;//[mm]		
		if((long)par->chkpt->co[i][ptMAXSWE]==number_novalue) par->chkpt->co[i][ptMAXSWE]=1.E10;//[mm]		
		if((long)par->chkpt->co[i][ptLAT]==number_novalue) par->chkpt->co[i][ptLAT]=par->latitude;		
		if((long)par->chkpt->co[i][ptLON]==number_novalue) par->chkpt->co[i][ptLON]=par->longitude;		
		if((long)par->chkpt->co[i][ptID]==number_novalue) par->chkpt->co[i][ptID]=(double)i;		
		if((long)par->chkpt->co[i][ptHOR]==number_novalue) par->chkpt->co[i][ptHOR]=par->chkpt->co[i][ptID];	
	}
	
	
	//i.show results
	fprintf(flog,"\nPOINTS:\n");
	fprintf(flog,"ID,East[m],North[m],Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-],CurvatureN-S[1/m],CurvatureW-E[1/m],CurvatureNW-SE[1/m],CurvatureNE-SW[1/m],LatDrainageDist[m],DepthFreeSurface[mm],Hor,maxSWE[mm],Lat[deg],Long[deg]\n");
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
	}
	
	//m. set IT->LU
	if(par->recover>0){
		IT->LU=new_doublematrix(Z->nrh, Z->nch);
		copy_doublematrix(LU, IT->LU);
	}
	
	//n. deallocation
	if(read_dem==1) free_doublematrix(Z);
	if(read_lu==1) free_doublematrix(LU);
	if(par->recover==0 && read_dem==1){
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
	top->BC_LatDistance=new_doublevector(par->chkpt->nrh);
	top->BC_DepthFreeSurface=new_doublevector(par->chkpt->nrh);
	par->maxSWE=new_doublematrix(1,par->chkpt->nrh);
	top->horizon_point=new_longmatrix(1,par->chkpt->nrh);
	top->dzdE=new_doublematrix(1,par->chkpt->nrh);
	top->dzdN=new_doublematrix(1,par->chkpt->nrh);
	top->latitude=new_doublematrix(1,par->chkpt->nrh);
	top->longitude=new_doublematrix(1,par->chkpt->nrh);
	par->IDpoint=new_longvector(par->chkpt->nrh);
	for(i=1;i<=par->chkpt->nrh;i++){
		if(par->state_px_coord==0){
			r=(long)par->chkpt->co[i][ptX];
			c=(long)par->chkpt->co[i][ptY];
			top->East->co[1][i] = UV->U->co[4] + (c-0.5)*UV->U->co[2];
			top->North->co[1][i] = UV->U->co[3] + (top->Z0->nrh-(r-0.5))*UV->U->co[1];
		}else{				
			top->East->co[1][i]=par->chkpt->co[i][ptX];
			top->North->co[1][i]=par->chkpt->co[i][ptY];
		}
		top->Z0->co[1][i]=par->chkpt->co[i][ptZ];
		land->LC->co[1][i]=par->chkpt->co[i][ptLC];
		
		if((long)land->LC->co[1][i] <= 0){
			printf("Error:: Point %ld has land cover type <= 0. This is not admitted.\n",i);
			stop_execution();
			t_error("Not Possible to Continue");
		}
		
		sl->type->co[1][i]=(long)par->chkpt->co[i][ptSY];
		
		if(sl->type->co[1][i] <= 0){
			printf("Error:: Point %ld has soil type <= 0. This is not admitted.\n",i);
			stop_execution();
			t_error("Not Possible to Continue");
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
		top->BC_LatDistance->co[i]=par->chkpt->co[i][ptDrDIST];
		top->BC_DepthFreeSurface->co[i]=par->chkpt->co[i][ptDrDEPTH];
		top->horizon_point->co[1][i]=(long)par->chkpt->co[i][ptHOR];
		top->dzdE->co[1][i]=0.;
		top->dzdN->co[1][i]=0.;		
		if(sl->type->co[1][i] <= 0){
			printf("Error:: Point %ld has horizon type <= 0. This is not admitted.\n",i);
			stop_execution();
			t_error("Not Possible to Continue");
		}
		
		par->maxSWE->co[1][i]=par->chkpt->co[i][ptMAXSWE];
		top->latitude->co[1][i]=par->chkpt->co[i][ptLAT];
		top->longitude->co[1][i]=par->chkpt->co[i][ptLON];
		par->IDpoint->co[i]=(long)par->chkpt->co[i][ptID];
	}
		
	//7. SET PAR
	par->output_soil=0.;
	par->output_snow=0.;
	par->output_glac=0.;
	par->output_surfenergy=0.;
	par->output_vegetation=0.;
	par->output_meteo=0.;

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
			if (c < par->chkpt->nrh) {
				if (top->horizon_point->co[1][c+1] != i) c++;
			}
		}while (top->horizon_point->co[1][c+1] != i && c < par->chkpt->nrh);
		
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

void set_bedrock(SOIL *sl, CHANNEL *cnet, PAR *par, TOPO *top, DOUBLEMATRIX *LC, FILE *flog){
	
	DOUBLEMATRIX *B;
	DOUBLETENSOR *T;
	DOUBLEVECTOR *WT;
	long i, j, l, r, c, sy, synew;
	double zlim, z;
	
	if(existing_file(files[fbed])==0){
		printf("File %s is missing. Please check if you have a bedrock topography map. If it is not available, remove the file name and keyword from input file\n",files[fbed]+1);
		t_error("Check input files");
	}

	printf("A bedrock depth map has been assigned and read from %s\n\n",files[fbed]);
	fprintf(flog,"A bedrock depth map has been assigned and read from %s\n\n",files[fbed]);

	par->bedrock = 1;
	B = read_map(2, files[fbed], LC, UV, (double)number_novalue);
	
	if (sl->init_water_table_height->nh != sl->pa->ndh) t_error("Error in bedrock calculations");

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
	WT=new_doublevector(sl->init_water_table_height->nh);
	for(i=1;i<=sl->init_water_table_height->nh;i++) {
		WT->co[i]=sl->init_water_table_height->co[i];
	}
	free_doublevector(sl->init_water_table_height);
	sl->init_water_table_height=new_doublevector(par->total_pixel+par->total_channel);

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
			z = 0.0;
		}else {
			r = cnet->r->co[i-par->total_pixel];
			c = cnet->c->co[i-par->total_pixel];
			sy = cnet->soil_type->co[i-par->total_pixel];
			synew = i;
			cnet->soil_type->co[i-par->total_pixel] = synew;
			z = par->depr_channel;
		}

		sl->init_water_table_height->co[synew] = WT->co[sy];
		
		zlim = B->co[r][c];

		for(l=1;l<=Nl;l++){			
				
			z += 0.5*sl->pa->co[synew][jdz][l];
				
			if(z <= zlim){
					
				for(j=1;j<=nsoilprop;j++){
					sl->pa->co[synew][j][l] = T->co[sy][j][l];
				}
					
			}else{
					
				for(j=1;j<=nsoilprop;j++){
					sl->pa->co[synew][j][l] = sl->pa_bed->co[sy][j][l] ;
				}				
					
			}
				
			z += 0.5*sl->pa->co[synew][jdz][l];
				
		}
	}
			
	free_doubletensor(T);
	free_doublematrix(B);
	free_doublevector(WT);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLETENSOR *find_Z_of_any_layer(DOUBLEMATRIX *Zsurface, DOUBLEMATRIX *slope, DOUBLEMATRIX *LC, SOIL *sl){
	
	DOUBLETENSOR *Z;
	double Zaverage, z;
	long l, r, c, n, sy;
	
	Zaverage=0.0;
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
	
	Z=new_doubletensor0(sl->pa->nch, Zsurface->nrh, Zsurface->nch);
	initialize_doubletensor(Z, (double)number_novalue);

	for(r=1;r<=Zsurface->nrh;r++){
		for(c=1;c<=Zsurface->nch;c++){
			if((long)LC->co[r][c]!=number_novalue){
						
				sy=sl->type->co[r][c];
				z=1.E3*(Zsurface->co[r][c]-Zaverage);//[mm]
			
				l=0;
				Z->co[l][r][c]=z;
			
				do{
					l++;
					z -= 0.5*sl->pa->co[sy][jdz][l]*cos(slope->co[r][c]*Pi/180.);
					Z->co[l][r][c]=z;
					z -= 0.5*sl->pa->co[sy][jdz][l]*cos(slope->co[r][c]*Pi/180.);
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
