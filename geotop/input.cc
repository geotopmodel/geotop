
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 1.225-15
 
 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
  
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
 
#include "config.h"
#include "input.h"
#include "parameters.h"
#include <unistd.h>
#include <inputKeywords.h>
#include "geotop_common.h"

using namespace std ;

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//! Subroutine which reads input data, performs  geomporphological analisys and allocates data

//void get_all_input(long argc, char *argv[], TOPO *top, SOIL *sl, LAN D *land, METEO *met, WATER *wat, CHANNEL *cnet,
//			    PAR *par, Energy *egy, SNOW *snow, GLACIER *glac, TIMES *times, mio::IOManager& iomanager)

void get_all_input(long argc, char *argv[], Topo *top, Soil *sl, Land *land, Meteo *met, Water *wat, Channel *cnet,
                   Par *par, Energy *egy, Snow *snow, Glacier *glac, Times *times, mio::IOManager& iomanager)

{

    FILE *flog, *f;
    //	DOUBLEMATRIX *M;
    GeoMatrix<double> M;
    //	INIT_TOOLS *IT;
    InitTools *IT;

    size_t a;
    short success, added_JDfrom0=0, added_wind_xy=0, added_wind_dir=0, added_cloud=0, added_Tdew=0, added_RH=0, added_Pint=0;
    long l, r, c, i, ist, j, n, sy, num_cols, num_lines, day, month, year, hour, minute;
    double z, th_oversat, JD, k_snowred, maxSWE, SWE, D, cosslope, **matrix;
    //	char *temp, **temp2;
    std::vector<std::string> temp2 ;
//    char  **temp2;
    string temp;

    //	IT=(INIT_TOOLS *)malloc(sizeof(INIT_TOOLS));
    IT = new InitTools();

    if (geotop::common::Variables::WORKING_DIRECTORY != "")
    {
        
    } else if(!argv[1]){
        geotop::common::Variables::WORKING_DIRECTORY=get_workingdirectory();
    }else if (argc==2){
        // modified by Emanuele Cordano on Aug 2011
        geotop::common::Variables::WORKING_DIRECTORY = argv[1] ;
    } else {
        // modified by Emanuele Cordano on Aug 2011
#ifdef USE_NETCDF
        //	WORKING_DIRECTORY=assign_string(read_option_string(argc,argv,"-wpath",".",0)); // assign_string(argv[1]); // MODIFY HERE EC
        WORKING_DIRECTORY=read_option_string(argc,argv,"-wpath",".",0); // assign_string(argv[1]); // MODIFY HERE EC
        //printf("WORKING_DIRECTORY=%s",WORKING_DIRECTORY);stop_execution();
#endif	
    }

    //add "/" if it is missing
    if (geotop::common::Variables::WORKING_DIRECTORY[strlen(geotop::common::Variables::WORKING_DIRECTORY.c_str())-1] != 47) {
        //	temp = assign_string(WORKING_DIRECTORY);
        temp = geotop::common::Variables::WORKING_DIRECTORY;
        //	free(WORKING_DIRECTORY);
        //	WORKING_DIRECTORY = join_strings(temp, "/");
        geotop::common::Variables::WORKING_DIRECTORY = temp + "/";
        //	free(temp);
    }

    //	logfile = join_strings(WORKING_DIRECTORY, logfile_name);
    geotop::common::Variables::logfile = geotop::common::Variables::WORKING_DIRECTORY + logfile_name;
    flog = fopen(geotop::common::Variables::logfile.c_str(), "w");

    printf("STATEMENT:\n");
    printf("\n");
    printf("GEOtop 1.225 'MOAB' - 9 Mar 2012\n\n");
    printf("Copyright (c), 2012 - Stefano Endrizzi \n\n");
    printf("TN -EXACT version (tmp)");
    printf("GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
    printf("WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
    printf("\nWORKING DIRECTORY: %s\n",geotop::common::Variables::WORKING_DIRECTORY.c_str());
    printf("\nLOGFILE: %s\n",geotop::common::Variables::logfile.c_str());

    fprintf(flog,"STATEMENT:\n");
    fprintf(flog,"\n");
    fprintf(flog,"GEOtop 1.225 'Moab' - 9 Mar 2012\n\n");
    fprintf(flog,"Copyright (c), 2012 - Stefano Endrizzi \n\n");
    fprintf(flog,"GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
    fprintf(flog,"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
    fprintf(flog,"\nWORKING DIRECTORY: %s\n",geotop::common::Variables::WORKING_DIRECTORY.c_str());

    //reads the parameters in __control_parameters
    //	temp = join_strings(WORKING_DIRECTORY, program_name);
    temp = geotop::common::Variables::WORKING_DIRECTORY + program_name;
    
    boost::shared_ptr<geotop::input::ConfigStore> lConfigStore = geotop::input::ConfigStoreSingletonFactory::getInstance() ;
    const std::string lFilePath (temp) ;
    bool lParsingRes = lConfigStore->parse(lFilePath) ;
    if(not lParsingRes)
    {
        t_error("Fatal Error! Unable to parse configuration file: " + lFilePath);
    }

    success = read_inpts_par(par, land, times, sl, met, IT, temp, flog);
    //	free(temp);

	//correct state pixel
	par->Tzrun = 0;
	par->wzrun = 0;
	par->Tzmaxrun = 0;
	par->Tzminrun = 0;
	par->wzmaxrun = 0;
	par->wzminrun = 0;	
	par->dUzrun = 0;
	par->SWErun = 0;	
	
	if(geotop::common::Variables::files[fTrun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->Tzrun = 1;
	}
	if(geotop::common::Variables::files[fwrun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->wzrun = 1;
	}
	if(geotop::common::Variables::files[fTmaxrun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->Tzmaxrun = 1;

	}
	if(geotop::common::Variables::files[fwmaxrun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->wzmaxrun = 1;
	}
	if(geotop::common::Variables::files[fTminrun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->Tzminrun = 1;
	}
	if(geotop::common::Variables::files[fwminrun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->wzminrun = 1;
	}
	if(geotop::common::Variables::files[fdUrun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->dUzrun = 1;
	}
	if(geotop::common::Variables::files[fSWErun] != geotop::input::gStringNoValue){
		if(par->point_sim == 1) par->state_pixel = 1;
		if(par->state_pixel == 1) par->SWErun = 1;
	}
	if (par->newperiodinit == 2 && (par->Tzrun == 0 || par->wzrun == 0)){
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f, "Error: You have to assign a name to the Tzrun and wzrun files\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");		
	}	

    std::cout << "SPAR: " << fspar << " : " << geotop::common::Variables::files[fspar] << std::endl ;
    success = read_soil_parameters(geotop::common::Variables::files[fspar], IT, sl, par->soil_type_bedr_default, flog);

    //Nl=sl->pa->nch;
    geotop::common::Variables::Nl = sl->pa.getCh() - 1;

    success = read_point_file(geotop::common::Variables::files[fpointlist], IT->point_col_names, par, flog);

    geotop::common::Variables::max_time = 0.;
    //	for (i=1; i<=par->init_date->nh; i++) {
    for (size_t i=1; i<par->init_date.size(); i++) {
        //	max_time += par->run_times->co[i]*(par->end_date->co[i] - par->init_date->co[i])*86400.;//seconds
        geotop::common::Variables::max_time += par->run_times[i]*(par->end_date[i] - par->init_date[i])*86400.;//seconds
    }

    //Recovering
    par->delay_day_recover = 0.0;
    par->n_ContRecovery = 0;

    if(par->recover > 0){
        //	if(par->saving_points->nh < par->recover){
        if(par->saving_points.size() -1< par->recover){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error: recover index higher than the length of the saving points vector");
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (1).");
        }
        //	par->delay_day_recover = par->saving_points->co[par->recover];
        par->delay_day_recover = par->saving_points[par->recover];
    }

    //Continuous Recovering
    if (par->RunIfAnOldRunIsPresent != 1) {

        if (mio::IOUtils::fileExists(geotop::common::Variables::SuccessfulRunFile)) {
            //	temp = join_strings(SuccessfulRunFile, ".old");
            temp = geotop::common::Variables::SuccessfulRunFile + ".old";
            rename(geotop::common::Variables::SuccessfulRunFile.c_str(), temp.c_str());
            //	free(temp);
            //	temp = join_strings(FailedRunFile, ".old");
            temp = geotop::common::Variables::FailedRunFile + ".old";
            rename(geotop::common::Variables::FailedRunFile.c_str(), temp.c_str());
            //	free(temp);
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "This simulation has successfully reached the end, you cannot recover it.\n");
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report.");
        }
        if (mio::IOUtils::fileExists(string(geotop::common::Variables::FailedRunFile))) {
            //	temp = join_strings(SuccessfulRunFile, ".old");
            temp = geotop::common::Variables::SuccessfulRunFile
            + ".old";
            rename(geotop::common::Variables::SuccessfulRunFile.c_str(), temp.c_str());
            //	free(temp);
            //	temp = join_strings(FailedRunFile, ".old");
            temp = geotop::common::Variables::FailedRunFile
            + ".old";
            rename(geotop::common::Variables::FailedRunFile.c_str(), temp.c_str());
            //	free(temp);
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "This simulation has failed, you cannot recover it.\n");
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report.");
        }
    }

    //	temp = join_strings(SuccessfulRunFile, ".old");
    temp = geotop::common::Variables::SuccessfulRunFile + ".old";
    rename(geotop::common::Variables::SuccessfulRunFile.c_str(), temp.c_str());
    //	free(temp);
    //	temp = join_strings(FailedRunFile, ".old");
    temp = geotop::common::Variables::FailedRunFile + ".old";
    rename(geotop::common::Variables::FailedRunFile.c_str(), temp.c_str());
    //	free(temp);

    if (par->ContRecovery > 0 && par->recover == 0) {
        if (mio::IOUtils::fileExists(string(geotop::common::Variables::files[rtime]) + string(textfile))) {
            //	temp = join_strings(geotop::common::Variables::files[rtime], textfile);
            temp = geotop::common::Variables::files[rtime] + string(textfile);
            matrix = read_txt_matrix_2(temp, 33, 44, 7, &num_lines);
            par->delay_day_recover = matrix[0][1];
            par->n_ContRecovery = (long)matrix[0][2];
            geotop::common::Variables::i_run0 = (long)matrix[0][3];
            geotop::common::Variables::i_sim0 = (long)matrix[0][4];
            geotop::common::Variables::cum_time = matrix[0][5];
            geotop::common::Variables::elapsed_time_start = matrix[0][6];
            for (i=0; i<num_lines; i++) {
                free(matrix[i]);
            }
            free(matrix);
            //	free(temp);
        }
    }

    //Time indices
    //	par->init_date->co[i_sim0] += par->delay_day_recover;
    par->init_date[geotop::common::Variables::i_sim0] += par->delay_day_recover;
    //	convert_JDfrom0_JDandYear(par->init_date->co[i_sim0], &JD, &year);
    convert_JDfrom0_JDandYear(par->init_date[geotop::common::Variables::i_sim0], &JD, &year);
    convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);

	geotop::common::Variables::i_run = geotop::common::Variables::i_run0;//Run index
		
    /****************************************************************************************************/
    /*! Reading of the Input files:                                                                     */
    /****************************************************************************************************/

    // ##################################################################################################################################
    // ##################################################################################################################################
    par->use_meteoio_cloud = false;
#ifndef USE_INTERNAL_METEODISTR
    //par->use_meteoio_meteodata=1;
    par->use_meteoio_cloud = true;
#endif

    meteoio_init(iomanager);
    // ##################################################################################################################################
    // ##################################################################################################################################

    if(par->point_sim!=1){  //distributed simulation
        read_inputmaps(top, land, sl, par, flog, iomanager);
    }else{
        read_optionsfile_point(par, top, land, sl, times, IT, flog);
    }

    //	Nr=top->Z0->nrh;
    geotop::common::Variables::Nr=top->Z0.getRows()-1;
    //	Nc=top->Z0->nch;
    geotop::common::Variables::Nc=top->Z0.getCols()-1;

    par->total_pixel=0;
    par->total_area=0.;
    for(r=1;r<=geotop::common::Variables::Nr;r++){
        for(c=1;c<=geotop::common::Variables::Nc;c++){
            //	if ((long)land->LC->co[r][c]!=geotop::input::gDoubleNoValue){
            if ((long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                par->total_pixel++;
                if(par->point_sim != 1){
                    //	par->total_area += (UV->U->co[1]*UV->U->co[2])/cos(top->slope->co[r][c]*GTConst::Pi/180.);
                    par->total_area += (geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2])/cos(top->slope[r][c]*GTConst::Pi/180.);
                }else {
                    par->total_area += (geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
                }
            }
        }
    }

    top->Z = find_Z_of_any_layer(top->Z0, top->slope, land->LC, sl, par->point_sim);

    //DOUBLEMATRIX *testm= new_doublematrix(10000, 3000);  // just to check if its possible to allocate -noori

    //	top->Jdown = new_longmatrix(par->total_pixel, 4);
    top->Jdown.resize(par->total_pixel+1, 4+1);
    //	top->Qdown = new_doublematrix(par->total_pixel, 4);
    top->Qdown.resize(par->total_pixel+1, 4+1);
    for (i=1; i<=par->total_pixel; i++) {
        for (j=1; j<=4; j++) {
            //	top->Jdown->co[i][j] = i;
            top->Jdown[i][j] = i;
            //	top->Qdown->co[i][j] = 0.;
            top->Qdown[i][j] = 0.;
        }
    }

    fprintf(flog,"Valid pixels: %ld\n",par->total_pixel);
    fprintf(flog,"Number of nodes: %ld\n",(geotop::common::Variables::Nl+1)*par->total_pixel);
    fprintf(flog,"Novalue pixels: %ld\n",(geotop::common::Variables::Nr*geotop::common::Variables::Nc-par->total_pixel));

    //	fprintf(flog,"Basin area: %f km2\n",(double)par->total_pixel*UV->U->co[1]*UV->U->co[2]/1.E6);

#ifdef USE_DOUBLE_PRECISION_OUTPUT
    fprintf(flog,"Basin area: %12g km2\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#else
    fprintf(flog,"Basin area: %f km2\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#endif

    printf("\nValid pixels: %ld\n",par->total_pixel);
    printf("Number of nodes: %ld\n",(geotop::common::Variables::Nl+1)*par->total_pixel);
    printf("Novalue pixels: %ld\n",(geotop::common::Variables::Nr*geotop::common::Variables::Nc-par->total_pixel));

    //	printf("Basin area: %f km2\n\n",(double)par->total_pixel*UV->U->co[1]*UV->U->co[2]/1.E6);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Basin area: %12g km2\n\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#else
    printf("Basin area: %f km2\n\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#endif
    /****************************************************************************************************/
    //Reading of RAIN data file,	METEO data file, 	and CLOUD data file

    num_cols = (long)nmet;

    //	meteo data
    //	met->data=(double***)malloc(met->st->E->nh*sizeof(double**));
    met->data=(double***)malloc(met->st->E.size()*sizeof(double**));
    //	number of line of meteo data
    //	met->numlines=(long*)malloc(met->st->E->nh*sizeof(long));
    met->numlines=(long*)malloc(met->st->E.size()*sizeof(long));

    //	horizon for meteo stations
    //	met->horizon=(double***)malloc(met->st->E->nh*sizeof(double**));
    met->horizon=(double***)malloc(met->st->E.size()*sizeof(double**));
    //	number of line in the horizon file
    //	met->horizonlines=(long*)malloc(met->st->E->nh*sizeof(long));
    met->horizonlines=(long*)malloc(met->st->E.size()*sizeof(long));
    //	line of met->data used (stored in memory to avoid from searching from the first line)

    success = read_meteostations_file(met->imeteo_stations, met->st, geotop::common::Variables::files[fmetstlist], IT->meteostations_col_names, flog);

    //	for(i=1;i<=met->st->E->nh;i++){
#ifdef USE_INTERNAL_METEODISTR
    met->var=(double**)malloc((met->st->E.size()-1)*sizeof(double*));
    met->line_interp_WEB=(long*)malloc(met->st->E.size()*sizeof(long));
    met->line_interp_Bsnow=(long*)malloc(met->st->E.size()*sizeof(long));
    met->line_interp_WEB_LR=0;
    met->line_interp_Bsnow_LR=0;
#endif
    long num_met_stat=met->st->E.size()-1;
    for(size_t i=1; i <= num_met_stat; i++){
        //	if (met->imeteo_stations->co[1] != geotop::input::gDoubleNoValue) {
        if (met->imeteo_stations[1] != geotop::input::gDoubleNoValue) {
            //	ist = met->imeteo_stations->co[i];
            ist = met->imeteo_stations[i];
        }else {
            ist = i;
        }
        
        //	read horizon
        met->horizon[i-1] = read_horizon(1, ist, geotop::common::Variables::files[fhormet], IT->horizon_col_names, &num_lines, flog);
        met->horizonlines[i-1] = num_lines;

        // ##################################################################################################################################
        // #####################Probably the next lines should not be necessary if we use meteoIO   #####################################
        // ##################################################################################################################################
#ifdef USE_INTERNAL_METEODISTR
        //initialize
        met->line_interp_WEB[i-1] = 0;
        met->line_interp_Bsnow[i-1] = 0;

        //allocate var
        met->var[i-1] = (double*)malloc(num_cols*sizeof(double));

        //filename
        if (geotop::common::Variables::files[fmet] != geotop::input::gStringNoValue){

            //read matrix
            temp=namefile_i(geotop::common::Variables::files[fmet], ist);

            met->data[i-1] = read_txt_matrix(temp, 33, 44, IT->met_col_names, nmet, &num_lines, flog);

            if ((long)met->data[i-1][0][iDate12] == geotop::input::gDoubleAbsent && (long)met->data[i-1][0][iJDfrom0] == geotop::input::gDoubleAbsent) {
                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                fprintf(f, "Error:: Date Column missing in file %s\n",temp.c_str());
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report (2).");
            }
            met->numlines[i-1] = num_lines;

            //fixing dates: converting times in the same standard time set for the simulation and fill JDfrom0
            short added_JDfrom0 = fixing_dates(ist, met->data[i-1], par->ST, met->st->ST[i], met->numlines[i-1], iDate12, iJDfrom0);

            check_times(ist, met->data[i-1], met->numlines[i-1], iJDfrom0);

            //find clouds
            if(IT->met_col_names[itauC] != geotop::input::gStringNoValue){
                if((long)met->data[i-1][0][itauC] == geotop::input::gDoubleAbsent || par->ric_cloud == 1){
			added_cloud = fill_meteo_data_with_cloudiness(met->data[i-1], met->numlines[i-1], met->horizon[i-1], met->horizonlines[i-1], 
			    met->st->lat[i], met->st->lon[i], par->ST, met->st->Z[i], met->st->sky[i], 0.0, par->ndivdaycloud, par->dem_rotation,
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

            // find Tdew
            if(par->vap_as_Td == 1){
                added_Tdew = fill_Tdew(i, met->st->Z, met->data[i-1], met->numlines[i-1], iRh, iT, iTdew, IT->met_col_names[iTdew], par->RHmin);
            }

            // find RH
            if(par->vap_as_RH == 1){
                added_RH = fill_RH(i,  met->st->Z, met->data[i-1], met->numlines[i-1], iRh, iT, iTdew, IT->met_col_names[iRh]);
            }

            //find Prec Intensity
            if ( par->linear_interpolation_meteo[i] == 1 && (long)met->data[i-1][0][iPrec] != geotop::input::gDoubleAbsent) {
                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                fprintf(f,"Meteo data for station %ld contain precipitation as volume, but Linear Interpolation is set. This is not possible, the precipitation data are removed.\n",i);
                fprintf(f,"If you want to use precipitation as volume, you cannot set keyword LinearInterpolation at 1.\n");
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report (3).");
            }

            if(par->prec_as_intensity == 1){
                added_Pint = fill_Pint(i, met->data[i-1], met->numlines[i-1], iPrec, iPrecInt, iJDfrom0, IT->met_col_names[iPrecInt]);
            }

            //rewrite completed files
            //rewrite_meteo_files(met->data[i-1], met->numlines[i-1], IT->met_col_names, temp.c_str(), added_JDfrom0, added_wind_xy, added_wind_dir, added_cloud, added_Tdew, added_RH, added_Pint);

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

            //free(temp);

        }else {

            fprintf(flog, "Warning: File meteo not in the list, meteo data not read, used default values\n");
            printf("Warning: File meteo not in the list, meteo data not read, used default values\n");

            met->data[i-1] = (double**)malloc(2*sizeof(double*));

            for (n=1; n<=2; n++) {
                met->data[i-1][n-1] = (double*)malloc(num_cols*sizeof(double));
                for (j=1; j<=nmet; j++) {
                    met->data[i-1][n-1][j-1] = (double)geotop::input::gDoubleAbsent;
                }
            }

            met->data[i-1][0][iJDfrom0] = 0.;
            met->data[i-1][1][iJDfrom0] = 1.E10;

        }
#endif
    }

    //read LAPSE RATES FILE

    if(geotop::common::Variables::files[fLRs] != geotop::input::gStringNoValue){   //s stands for string

        if (!mio::IOUtils::fileExists(string(geotop::common::Variables::files[fLRs]) + string(textfile)))
            printf("Lapse rate file unavailable. Check input files. If you do not have a lapse rate file, remove its name and keywords from input file\n");
        //	temp = join_strings(geotop::common::Variables::files[fLRs], textfile);
        temp = geotop::common::Variables::files[fLRs] + string(textfile);
        met->LRs = read_txt_matrix(temp, 33, 44, IT->lapserates_col_names, nlstot, &num_lines, flog);
        //	free(temp);
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
        met->LRv[i] = geotop::input::gDoubleNoValue;
        if (i == ilsTa) {
            met->LRd[i] = GTConst::LapseRateTair;
        }else if (i == ilsTdew) {
            met->LRd[i] = GTConst::LapseRateTdew;
        }else if (i == ilsPrec) {
            met->LRd[i] = GTConst::LapseRatePrec;
        }else {
            met->LRd[i] = 0.0;
        }
    }


    //	FIND A STATION WITH SHORTWAVE RADIATION DATA
#ifdef USE_INTERNAL_METEODISTR
    met->nstsrad=0;
    do{
        met->nstsrad++;
        a=0;
//			printf("met->nstsrad-1=%ld",met->nstsrad-1);
//			printf("met->data[met->nstsrad-1][0][iSW]=%ld",(long)met->data[met->nstsrad-1][0][iSW]);
//			printf("met->data[met->nstsrad-1][0][iSWb]=%ld",(long)met->data[met->nstsrad-1][0][iSWb]);
//			printf("met->data[met->nstsrad-1][0][iSWd]=%ld",(long)met->data[met->nstsrad-1][0][iSWd]);
        if( (long)met->data[met->nstsrad-1][0][iSW]!=geotop::input::gDoubleAbsent || ((long)met->data[met->nstsrad-1][0][iSWb]!=geotop::input::gDoubleAbsent && (long)met->data[met->nstsrad-1][0][iSWd]!=geotop::input::gDoubleAbsent ) ) a=1;
		}while(met->nstsrad< num_met_stat && a==0);
    if(a==0){
        printf("WARNING: NO shortwave radiation measurements available\n");
        fprintf(flog,"WARNING: NO shortwave radiation measurements available\n");
    }else{
        printf("Shortwave radiation measurements from station %ld\n",met->nstsrad);
        fprintf(flog,"Shortwave radiation measurements from station %ld\n",met->nstsrad);
    }
#endif

#ifdef USE_INTERNAL_METEODISTR
	//FIND A STATION WITH CLOUD DATA
	met->nstcloud=0;
	do{
		met->nstcloud++;
		a=0;
		if( (long)met->data[met->nstcloud-1][0][iC]!=geotop::input::gDoubleAbsent || (long)met->data[met->nstcloud-1][0][itauC]!=geotop::input::gDoubleAbsent ) a=1;
	}while(met->nstcloud<met->st->Z.size()-1 && a==0);
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
		if( (long)met->data[met->nstlrad-1][0][iLWi]!=geotop::input::gDoubleAbsent) a=1;
	}while(met->nstlrad<met->st->Z.size()-1 && a==0);
	if(a==0){
		printf("WARNING: NO longwave radiation measurements available\n");
		fprintf(flog,"WARNING: NO longwave radiation measurements available\n");
	}else{
		printf("Longwave radiation measurements from station %ld\n",met->nstlrad);
		fprintf(flog,"Longwave radiation measurements from station %ld\n",met->nstlrad);
	}

#endif

    //	met->tau_cl_map=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    //	initialize_doublematrix(met->tau_cl_map, geotop::input::gDoubleNoValue);
    met->tau_cl_map.resize(top->Z0.getRows(),top->Z0.getCols(),geotop::input::gDoubleNoValue);

    //	met->tau_cl_av_map=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    //	initialize_doublematrix(met->tau_cl_av_map, geotop::input::gDoubleNoValue);
    met->tau_cl_av_map.resize(top->Z0.getRows(),top->Z0.getCols(), geotop::input::gDoubleNoValue);

    //	met->tau_cl_map_yes=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
    //	initialize_shortmatrix(met->tau_cl_map_yes, (short)geotop::input::gDoubleNoValue);
    met->tau_cl_map_yes.resize(top->Z0.getRows(),top->Z0.getCols(), (short)geotop::input::gDoubleNoValue);

    //	met->tau_cl_av_map_yes=new_shortmatrix(top->Z0->nrh,top->Z0->nch);
    //	initialize_shortmatrix(met->tau_cl_av_map_yes, (short)geotop::input::gDoubleNoValue);
    met->tau_cl_av_map_yes.resize(top->Z0.getRows(),top->Z0.getCols(), (short)geotop::input::gDoubleNoValue);

    //	vector defining which meteo station has the SW radiation information
    //	met->st->flag_SW_meteoST=new_shortvector(met->st->Z->nh);
    //	initialize_shortvector(met->st->flag_SW_meteoST,geotop::input::gDoubleNoValue);
    met->st->flag_SW_meteoST.resize(met->st->Z.size(),geotop::input::gDoubleNoValue);

    //	met->st->tau_cloud_av_yes_meteoST=new_shortvector(met->st->Z->nh);
    //	initialize_shortvector(met->st->tau_cloud_av_yes_meteoST,geotop::input::gDoubleNoValue);
    met->st->tau_cloud_av_yes_meteoST.resize(met->st->Z.size(), geotop::input::gDoubleNoValue);

    //	met->st->tau_cloud_yes_meteoST=new_shortvector(met->st->Z->nh);
    //	initialize_shortvector(met->st->tau_cloud_yes_meteoST,geotop::input::gDoubleNoValue);
    met->st->tau_cloud_yes_meteoST.resize(met->st->Z.size(),geotop::input::gDoubleNoValue);

    //	met->st->tau_cloud_av_meteoST=new_doublevector(met->st->Z->nh);
    //	initialize_doublevector(met->st->tau_cloud_av_meteoST,geotop::input::gDoubleNoValue);
    met->st->tau_cloud_av_meteoST.resize(met->st->Z.size(), geotop::input::gDoubleNoValue);
    //	met->st->tau_cloud_meteoST=new_doublevector(met->st->Z->nh);
    //	initialize_doublevector(met->st->tau_cloud_meteoST,geotop::input::gDoubleNoValue);
    met->st->tau_cloud_meteoST.resize(met->st->Z.size(),geotop::input::gDoubleNoValue);

    //i.show details on checkpoints
    fprintf(flog,"\nCHECKPOINTS:\n");
    fprintf(flog,"ID,r,c,Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-]\n");
    //	for(i=1;i<=par->rc->nrh;i++){
    for(i=1;i<par->rc.getRows();i++){
        //	r=par->rc->co[i][1];
        r=par->rc[i][1];
        //	c=par->rc->co[i][2];
        c=par->rc[i][2];
        //	fprintf(flog,"%ld,%ld,%ld,%f,%d,%ld,%f,%f,%f\n",i,r,c,top->Z0->co[r][c],(short)land->LC->co[r][c],sl->type->co[r][c],top->slope->co[r][c],top->aspect->co[r][c],top->sky->co[r][c]);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
        fprintf(flog,"%ld,%ld,%ld,%12g,%d,%ld,%12g,%12g,%12g\n",i,r,c,top->Z0[r][c],(short)land->LC[r][c],sl->type[r][c],top->slope[r][c],top->aspect[r][c],top->sky[r][c]);
#else
        fprintf(flog,"%ld,%ld,%ld,%f,%d,%ld,%f,%f,%f\n",i,r,c,top->Z0[r][c],(short)land->LC[r][c],sl->type[r][c],top->slope[r][c],top->aspect[r][c],top->sky[r][c]);
#endif
    }

    //i.show meteo stations
    fprintf(flog,"\nMETEO STATIONS:\n");
    fprintf(flog,"ID,East[m],North[m],Lat[deg],Lon[deg],Elev[m a.s.l.],Sky[deg],Stand_Time[h],WindSensHeight[m],TempSensHeight[m]\n");
    //	for(r=1;r<=met->st->E->nh;r++){
    for(size_t r=1;r<met->st->E.size();r++){
        //	fprintf(flog,"%ld,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",r,met->st->E->co[r],met->st->N->co[r],met->st->lat->co[r], met->st->lon->co[r],
        //			met->st->Z->co[r],met->st->sky->co[r],met->st->ST->co[r],met->st->Vheight->co[r],met->st->Theight->co[r]);

#ifdef USE_DOUBLE_PRECISION_OUTPUT
        fprintf(flog,"%ld,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g\n",r,met->st->E[r],met->st->N[r],met->st->lat[r], met->st->lon[r],
                met->st->Z[r],met->st->sky[r],met->st->ST[r],met->st->Vheight[r],met->st->Theight[r]);
#else
        fprintf(flog,"%ld,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",r,met->st->E[r],met->st->N[r],met->st->lat[r], met->st->lon[r],
                met->st->Z[r],met->st->sky[r],met->st->ST[r],met->st->Vheight[r],met->st->Theight[r]);
#endif
    }

    /****************************************************************************************************/
    //read INCOMING DISCHARGE

    met->qinv = (double*)malloc(2*sizeof(double));
    met->qinv[0] = 0.;
    met->qinv[1] = 0.;

    if(geotop::common::Variables::files[fqin] != geotop::input::gStringNoValue){
        //	temp = join_strings(geotop::common::Variables::files[fqin], textfile);
        temp =geotop::common::Variables::files[fqin] + string(textfile);
        temp2.push_back("Time") ;
        temp2.push_back("Qx") ;
        met->qins = read_txt_matrix(temp, 33, 44, temp2, 2, &num_lines, flog);
        //	free(temp);
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

    //	Initialize matrix of shadow
    //	land->shadow=new_shortmatrix(Nr,Nc);
    //	initialize_shortmatrix(land->shadow,0);/* initialized as if it was always NOT in shadow*/
    land->shadow.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0);

    //	Check that there aren't cell with an undefined land use value
    z = 0.;
    l = 0;
    do{
        l++;
        //	z += sl->pa->co[1][jdz][l];
        z += sl->pa(1,jdz,l);
    }while(l<geotop::common::Variables::Nl && z < GTConst::z_transp);
    //	land->root_fraction=new_doublematrix(par->n_landuses, l);
    //	initialize_doublematrix(land->root_fraction, 0.0);
    land->root_fraction.resize(par->n_landuses+1, l+1 ,0.0);

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

    //	variables used to assign vegetation properties that change with time
    num_cols = jdvegprop + 1;
    land->vegpars=(double ***)malloc(par->n_landuses*sizeof(double**));
    land->vegparv=(double **)malloc(par->n_landuses*sizeof(double*));
    land->NumlinesVegTimeDepData=(long*)malloc(par->n_landuses*sizeof(long));

    //	land->vegpar=new_doublevector(jdvegprop);
    land->vegpar.resize(jdvegprop+1);

    //	par->vegflag=new_shortvector(par->n_landuses);
    //	initialize_shortvector(par->vegflag, 0);
    par->vegflag.resize(par->n_landuses+1,0);

    //	time dependent vegetation parameters
    for(i=1;i<=par->n_landuses;i++){

        if (geotop::common::Variables::files[fvegpar] != geotop::input::gStringNoValue) {   //s stands for string

            temp = namefile_i_we2(geotop::common::Variables::files[fvegpar], i);

            if (mio::IOUtils::fileExists(string(temp) + string(textfile))) {
                printf("There is a specific vegetation parameter file for land cover type = %ld\n",i);
                //	free(temp);
                temp = namefile_i(geotop::common::Variables::files[fvegpar], i);
                land->vegpars[i-1] = read_txt_matrix_2(temp, 33, 44, num_cols, &num_lines);
                //	free(temp);
                land->NumlinesVegTimeDepData[i-1] = num_lines;
                //	par->vegflag->co[i]=1;
                par->vegflag[i]=1;
            } else {
                //	free(temp);
                printf("There is NOT a specific vegetation parameter file for land cover type = %ld\n",i);
            }
        }

        land->vegparv[i-1]=(double*)malloc(num_cols*sizeof(double));
        for(j=0; j<num_cols; j++){
            land->vegparv[i-1][j] = geotop::input::gDoubleNoValue;
        }

        //	z0 (convert in m)
        //	land->ty->co[i][jz0]*=0.001;
        land->ty[i][jz0]*=0.001;

        //	find root fraction
        //	root(land->root_fraction->nch, land->ty->co[i][jroot], 0.0, sl->pa->co[1][jdz], land->root_fraction->co[i]);
        //	root(land->root_fraction.getCols(), land->ty[i][jroot], 0.0, sl->pa->co[1][jdz], land->root_fraction, i);
        root(land->root_fraction.getCols(), land->ty[i][jroot], 0.0, sl->pa, land->root_fraction, i);

        //	error messages
        //	for(l=1;l<=met->st->Z->nh;l++){
        for(size_t l=1;l<met->st->Z.size();l++){
            //	if(0.001*land->ty->co[i][jHveg]>met->st->Vheight->co[l] || 0.001*land->ty->co[i][jHveg]>met->st->Theight->co[l]){
            if(0.001*land->ty[i][jHveg]>met->st->Vheight[l] || 0.001*land->ty[i][jHveg]>met->st->Theight[l]){
                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                //	fprintf(f, "hc:%f m, zmu:%f m, zmt:%f m - hc must be lower than measurement height - land cover %ld, meteo station %ld\n",
                //								0.001*land->ty->co[i][jHveg],met->st->Vheight->co[l],met->st->Theight->co[l],i,l);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
                fprintf(f, "hc:%12g m, zmu:%12g m, zmt:%12g m - hc must be lower than measurement height - land cover %ld, meteo station %ld\n",
                        0.001*land->ty[i][jHveg],met->st->Vheight[l],met->st->Theight[l],i,l);
#else
                fprintf(f, "hc:%f m, zmu:%f m, zmt:%f m - hc must be lower than measurement height - land cover %ld, meteo station %ld\n",
                        0.001*land->ty[i][jHveg],met->st->Vheight[l],met->st->Theight[l],i,l);
#endif

                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report (5).");
            }
        }
    }

    //******************************************
    //file with time steps

    if(geotop::common::Variables::files[ftsteps] != geotop::input::gStringNoValue){

        //	temp=join_strings(geotop::common::Variables::files[ftsteps],textfile);
        temp=geotop::common::Variables::files[ftsteps] + string(textfile);
        times->Dt_matrix = read_txt_matrix_2(temp, 33, 44, GTConst::max_cols_time_steps_file+1, &num_lines);
        //	free(temp);
        times->numlinesDt_matrix = num_lines;
        par->tsteps_from_file=1;

    }else{

        par->tsteps_from_file=0;

    }

    /****************************************************************************************************/
    /*! Filling up of the struct "channel" (of the type CHANNEL):                                        */

    /*The number of channel-pixel are counted:*/
    i=0;
    for(r=1;r<=geotop::common::Variables::Nr;r++){
        for(c=1;c<=geotop::common::Variables::Nc;c++){
            //	if (top->pixel_type->co[r][c]>=10) i++;
            if (top->pixel_type[r][c]>=10) i++;
        }
    }
    fprintf(flog,"Channel pixels: %ld\n",i);
    par->total_channel = i;

    //allocate channel vectors/matrixes
    if(i==0) i=1;

    cnet->Vout = 0.;

    //	cnet->r=new_longvector(i);
    //	initialize_longvector(cnet->r, 0);
    cnet->r.resize(i+1,0);

    //	cnet->c=new_longvector(i);
    //	initialize_longvector(cnet->c, 0);
    cnet->c.resize(i+1,0);

    //	cnet->ch=new_longmatrix(Nr,Nc);
    //	initialize_longmatrix(cnet->ch, 0);
    cnet->ch.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0);

    //	cnet->ch_down=new_longvector(i);
    //	initialize_longvector(cnet->ch_down, 0);
    cnet->ch_down.resize(i+1,0);

    //	cnet->length=new_doublevector(i);
    //	initialize_doublevector(cnet->length, 0.);
    cnet->length.resize(i+1,0.0);

    //	cnet->Vsup=new_doublevector(i);
    //	initialize_doublevector(cnet->Vsup, 0.);
    cnet->Vsup.resize(i+1,0.);
    //	cnet->Vsub=new_doublevector(i);
    //	initialize_doublevector(cnet->Vsub, 0.);
    cnet->Vsub.resize(i+1,0.);

    //	cnet->h_sup=new_doublevector(i);
    //	initialize_doublevector(cnet->h_sup, 0.);
    cnet->h_sup.resize(i+1,0.0);

    //	cnet->soil_type = new_longvector(i);
    //	initialize_longvector(cnet->soil_type, par->soil_type_chan_default);
    cnet->soil_type.resize(i+1, par->soil_type_chan_default);

    if(par->total_channel>1) enumerate_channels(cnet, land->LC, top->pixel_type, top->Z0, top->slope, geotop::input::gDoubleNoValue);

    cnet->ch3 = (long**)malloc((geotop::common::Variables::Nl+1)*sizeof(long*));
    for (l=0; l<=geotop::common::Variables::Nl; l++) {
        cnet->ch3[l] = (long*)malloc((i+1)*sizeof(long));
    }

    //	cnet->lch = new_longmatrix( (geotop::common::Variables::Nl+1)*i, 2);
    //	initialize_longmatrix(cnet->lch, 0);
    cnet->lch.resize( ((geotop::common::Variables::Nl+1)*i)+1, 2+1 , 0);

    //	lch3_cont(cnet->r->nh, cnet->ch3, cnet->lch);
    lch3_cont(cnet->ch3, cnet->lch,geotop::common::Variables::Nl,par->total_channel);


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
    for(r=1;r<=geotop::common::Variables::Nr;r++){
        C=Nc;
        for (i=1; i<=cnet->r->nh; i++) {
            if (r==cnet->r->co[i]) C=cnet->c->co[i];
        }
        for(c=1;c<=geotop::common::Variables::Nc;c++){
            if((long)land->LC->co[r][c]!=geotop::input::gDoubleNoValue){
                M->co[r][c]*=UV->U->co[1];
                if(M->co[r][c]>45) M->co[r][c]=geotop::input::gDoubleNoValue;
                if(c>=C) M->co[r][c]=geotop::input::gDoubleNoValue;
                if(r<=R) M->co[r][c]=geotop::input::gDoubleNoValue;
            }
        }
    }

    temp=join_strings(WORKING_DIRECTORY, "dist_from_channel");
    write_map(temp, 0, par->format_out, M, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
    free(temp);
    free_doublematrix(M);*/

    //Cont for Richards 3D
    //There are not channels
    top->i_cont=(long ***)malloc((geotop::common::Variables::Nl+1)*sizeof(long**));
    for(l=0;l<=geotop::common::Variables::Nl;l++){
        top->i_cont[l]=(long **)malloc((geotop::common::Variables::Nr+1)*sizeof(long*));
        for(r=1;r<=geotop::common::Variables::Nr;r++){
            top->i_cont[l][r]=(long *)malloc((geotop::common::Variables::Nc+1)*sizeof(long));
        }
    }

    //	top->lrc_cont=new_longmatrix( (geotop::common::Variables::Nl+1)*par->total_pixel, 3);
    //	initialize_longmatrix(top->lrc_cont, 0);
    top->lrc_cont.resize( (geotop::common::Variables::Nl+1)*par->total_pixel+1, 3+1, 0);

    i_lrc_cont(land->LC, top->i_cont, top->lrc_cont);

    top->j_cont=(long **)malloc((geotop::common::Variables::Nr+1)*sizeof(long*));
    for (r=1; r<=geotop::common::Variables::Nr; r++) {
        top->j_cont[r]=(long *)malloc((geotop::common::Variables::Nc+1)*sizeof(long));
        for (c=1; c<=geotop::common::Variables::Nc; c++) {
            top->j_cont[r][c] = 0;
        }
    }

    //	top->rc_cont=new_longmatrix(par->total_pixel, 2);
    //	initialize_longmatrix(top->rc_cont, 0);
    top->rc_cont.resize(par->total_pixel+1, 2+1);

    j_rc_cont(land->LC, top->j_cont, top->rc_cont);

    if(par->state_pixel == 1){
        //	par->jplot = new_longvector(par->total_pixel);
        //	initialize_longvector(par->jplot, 0);
        par->jplot.resize(par->total_pixel+1, 0);

        for (i=1; i<=par->total_pixel; i++) {
            //	for (j=1; j<=par->rc->nrh; j++) {
            for (j=1; j<par->rc.getRows(); j++) {
                //	if (top->rc_cont->co[i][1] == par->rc->co[j][1] && top->rc_cont->co[i][2] == par->rc->co[j][2]) {
                if (top->rc_cont[i][1] == par->rc[j][1] && top->rc_cont[i][2] == par->rc[j][2]) {
                    //	par->jplot->co[i] = j;
                    par->jplot[i] = j;
                }
            }
        }
    }

    //BEDROCK (adjusting soil properties)
    par->bedrock = 0;
    if(geotop::common::Variables::files[fbed] != geotop::input::gStringNoValue) set_bedrock(sl, cnet, par, top, land->LC, flog);
    //free_doubletensor(sl->pa_bed);

    /****************************************************************************************************/
    /*! Completing of the initialization of SOIL structure                               */
    /****************************************************************************************************/

    //	sl->SS = (SOIL_STATE *)malloc(sizeof(SOIL_STATE));
    sl->SS =new SoilState();
    initialize_soil_state(sl->SS, par->total_pixel, geotop::common::Variables::Nl);

    //	sl->VS = (STATE_VEG *)malloc(sizeof(STATE_VEG));
    sl->VS =new StateVeg();
    initialize_veg_state(sl->VS, par->total_pixel);

    //	sl->th=new_doublematrix(geotop::common::Variables::Nl,par->total_pixel);
    //	initialize_doublematrix(sl->th,geotop::input::gDoubleNoValue);
    sl->th.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,geotop::input::gDoubleNoValue);

    //	sl->Ptot=new_doublematrix(geotop::common::Variables::Nl,par->total_pixel);
    //	initialize_doublematrix(sl->Ptot,geotop::input::gDoubleNoValue);
    sl->Ptot.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,geotop::input::gDoubleNoValue);

    if(geotop::common::Variables::files[fTav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTavsup] != geotop::input::gStringNoValue){
        //	sl->T_av_tensor=new_doublematrix(geotop::common::Variables::Nl,par->total_pixel);
        //	initialize_doublematrix(sl->T_av_tensor,0.);
        sl->T_av_tensor.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,0.0);
    }

    if(geotop::common::Variables::files[ficeav] != geotop::input::gStringNoValue){
        //	sl->thi_av_tensor=new_doublematrix(geotop::common::Variables::Nl,par->total_pixel);
        //	initialize_doublematrix(sl->thi_av_tensor,0.);
        sl->thi_av_tensor.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,0.0);
    }

    if(geotop::common::Variables::files[fliqav] != geotop::input::gStringNoValue){
        //		sl->thw_av_tensor=new_doublematrix(geotop::common::Variables::Nl,par->total_pixel);
        //		initialize_doublematrix(sl->thw_av_tensor,0.);
        sl->thw_av_tensor.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,0.0);
    }

    //	sl->ET=new_doubletensor(geotop::common::Variables::Nl,Nr,Nc);
    //	initialize_doubletensor(sl->ET,0.);
    sl->ET.resize(geotop::common::Variables::Nl+1,geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.);

    if (!mio::IOUtils::fileExists(string(geotop::common::Variables::files[fwt0]) + string(ascii_esri))){

        for (i=1; i<=par->total_pixel; i++) {

            //	r = top->rc_cont->co[i][1];
            r = top->rc_cont[i][1];
            //	c = top->rc_cont->co[i][2];
            c = top->rc_cont[i][2];

            //	sy=sl->type->co[r][c];
            sy=sl->type[r][c];

            //	if ((long)sl->init_water_table_depth->co[sy] != geotop::input::gDoubleNoValue) {
            if ((long)IT->init_water_table_depth[sy] != geotop::input::gDoubleNoValue) {
                z = 0.;
                //	sl->SS->P->co[0][i] = -sl->init_water_table_depth->co[sy]*cos(top->slope->co[r][c]*GTConst::Pi/180.);
                sl->SS->P[0][i] = -IT->init_water_table_depth[sy]*cos(top->slope[r][c]*GTConst::Pi/180.);

                for(l=1;l<=geotop::common::Variables::Nl;l++){
                    //	z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*GTConst::Pi/180.);
                    z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
                    //	sl->SS->P->co[l][i] = sl->SS->P->co[0][i] + z;
                    sl->SS->P[l][i] = sl->SS->P[0][i] + z;
                    //	z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*GTConst::Pi/180.);
                    z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
                }
            }else {
                for(l=1;l<=geotop::common::Variables::Nl;l++){
                    //	sl->SS->P->co[l][i] = sl->pa->co[sy][jpsi][l];
                    sl->SS->P[l][i] = sl->pa(sy,jpsi,l);
                }
            }
        }

    }else {

        //	M = read_map(2,geotop::common::Variables::files[fwt0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fwt0]), M);

        for (i=1; i<=par->total_pixel; i++) {
            //	r = top->rc_cont->co[i][1];
            r = top->rc_cont[i][1];
            //	c = top->rc_cont->co[i][2];
            c = top->rc_cont[i][2];

            //	sy=sl->type->co[r][c];
            sy=sl->type[r][c];

            z = 0.;
            //	sl->SS->P->co[0][i] = -M->co[r][c]*cos(top->slope->co[r][c]*GTConst::Pi/180.);
            sl->SS->P[0][i] = -M[r][c]*cos(top->slope[r][c]*GTConst::Pi/180.);
            for(l=1;l<=geotop::common::Variables::Nl;l++){
                //	z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*GTConst::Pi/180.);
                z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
                //	sl->SS->P->co[l][i] = sl->SS->P->co[0][i] + z;
                sl->SS->P[l][i] = sl->SS->P[0][i] + z;
                //	z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*GTConst::Pi/180.);
                z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
            }
        }

        //	free_doublematrix(M);
    }

    for (i=1; i<=par->total_pixel; i++) {

        //	r = top->rc_cont->co[i][1];
        r = top->rc_cont[i][1];
        //	c = top->rc_cont->co[i][2];
        c = top->rc_cont[i][2];

        //	sy=sl->type->co[r][c];
        sy=sl->type[r][c];

        for(l=1;l<=geotop::common::Variables::Nl;l++){

            //	sl->SS->T->co[l][i]=sl->pa->co[sy][jT][l];
            sl->SS->T[l][i]=sl->pa(sy,jT,l);

            //	sl->Ptot->co[l][i] = sl->SS->P->co[l][i];
            sl->Ptot[l][i] = sl->SS->P[l][i];
            //	sl->th->co[l][i] = teta_psi(sl->SS->P->co[l][i], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
            //								   sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], GTConst::PsiMin, sl->pa->co[sy][jss][l]);

            sl->th[l][i] = teta_psi(sl->SS->P[l][i], 0.0, sl->pa(sy,jsat,l), sl->pa(sy,jres,l), sl->pa(sy,ja,l),
                                    sl->pa(sy,jns,l), 1-1/sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));

            //	th_oversat = Fmax( sl->SS->P->co[l][i] , 0.0 ) * sl->pa->co[sy][jss][l];
            th_oversat = Fmax( sl->SS->P[l][i] , 0.0 ) * sl->pa(sy,jss,l);
            //	sl->th->co[l][i] -= th_oversat;
            sl->th[l][i] -= th_oversat;

            //	if(sl->SS->T->co[l][i]<=GTConst::Tfreezing){
            if(sl->SS->T[l][i]<=GTConst::Tfreezing){

                //	Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
                //	sl->SS->thi->co[l][i] = sl->th->co[l][i] - teta_psi(Psif(sl->SS->T->co[l][i]), 0.0, sl->pa->co[sy][jsat][l],
                //															sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
                //															1.-1./sl->pa->co[sy][jns][l], GTConst::PsiMin, sl->pa->co[sy][jss][l]);

                sl->SS->thi[l][i] = sl->th[l][i] - teta_psi(Psif(sl->SS->T[l][i]), 0.0, sl->pa(sy,jsat,l),
                                                            sl->pa(sy,jres,l), sl->pa(sy,ja,l), sl->pa(sy,jns,l), 1.-1./sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));

                //	if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
                if(sl->SS->thi[l][i]<0) sl->SS->thi[l][i]=0.0;

                //	Psi is updated taking into account the freezing
                //	sl->th->co[l][i] -= sl->SS->thi->co[l][i];
                sl->th[l][i] -= sl->SS->thi[l][i];

                //	sl->SS->P->co[l][i] = psi_teta(sl->th->co[l][i] + th_oversat, sl->SS->thi->co[l][i], sl->pa->co[sy][jsat][l],
                //								   sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
                //								   1-1/sl->pa->co[sy][jns][l], GTConst::PsiMin, sl->pa->co[sy][jss][l]);
                sl->SS->P[l][i] = psi_teta(sl->th[l][i] + th_oversat, sl->SS->thi[l][i], sl->pa(sy,jsat,l),
                                           sl->pa(sy,jres,l), sl->pa(sy,ja,l), sl->pa(sy,jns,l),
                                           1-1/sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));
            }
        }
    }


    if(par->state_pixel == 1){
        if(geotop::common::Variables::files[fTz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTzwriteend] != geotop::input::gStringNoValue) sl->Tzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fTzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTzavwriteend] != geotop::input::gStringNoValue) sl->Tzavplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fpsiztot] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fpsiztotwriteend] != geotop::input::gStringNoValue) sl->Ptotzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fpsiz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fpsizwriteend] != geotop::input::gStringNoValue) sl->Pzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fliqz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fliqzwriteend] != geotop::input::gStringNoValue) sl->thzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fliqzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fliqzavwriteend] != geotop::input::gStringNoValue) sl->thzavplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[ficez] != geotop::input::gStringNoValue ||geotop::common::Variables::files[ficezwriteend] != geotop::input::gStringNoValue) sl->thizplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[ficezav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[ficezavwriteend] != geotop::input::gStringNoValue) sl->thizavplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);

        //	for (i=1; i<=par->rc->nrh; i++) {
        for (i=1; i<par->rc.getRows(); i++) {
            //	r = top->rc_cont->co[i][1];
            r = top->rc_cont[i][1];
            //	c = top->rc_cont->co[i][2];
            c = top->rc_cont[i][2];
            j = top->j_cont[r][c];
            for(l=1;l<=geotop::common::Variables::Nl;l++){
                if(geotop::common::Variables::files[fTz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTzwriteend] != geotop::input::gStringNoValue) sl->Tzplot[i][l] = sl->SS->T[l][j];
                if(geotop::common::Variables::files[fTzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTzavwriteend] != geotop::input::gStringNoValue) sl->Tzavplot[i][l] = sl->SS->T[l][j];
                if(geotop::common::Variables::files[fliqzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fliqzavwriteend] != geotop::input::gStringNoValue) sl->thzavplot[i][l] = sl->th[l][j];
                if(geotop::common::Variables::files[ficez] != geotop::input::gStringNoValue ||geotop::common::Variables::files[ficezwriteend] != geotop::input::gStringNoValue) sl->thizplot[i][l] = sl->SS->thi[l][j];
                if(geotop::common::Variables::files[ficezav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[ficezavwriteend] != geotop::input::gStringNoValue) sl->thizavplot[i][l] = sl->SS->thi[l][j];
                if(geotop::common::Variables::files[fpsiztot] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fpsiztotwriteend] != geotop::input::gStringNoValue) sl->Ptotzplot[i][l] = sl->Ptot[l][j];
            }
            for(l=0;l<=geotop::common::Variables::Nl;l++) {
                if(geotop::common::Variables::files[fpsiz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fpsizwriteend] != geotop::input::gStringNoValue) sl->Pzplot[i][l] = sl->SS->P[l][j];
            }
        }

    }

	printf("%f\n",par->delay_day_recover);
	
    if(par->delay_day_recover > 0){

        assign_recovered_tensor_vector(1, par->recover, geotop::common::Variables::files[riceg], sl->SS->thi, top->rc_cont, par, land->LC, IT->LU);
        assign_recovered_tensor_vector(1, par->recover, geotop::common::Variables::files[rTg], sl->SS->T, top->rc_cont, par, land->LC, IT->LU);
        assign_recovered_tensor_vector(0, par->recover, geotop::common::Variables::files[rpsi], sl->SS->P, top->rc_cont, par, land->LC, IT->LU);

        assign_recovered_map_vector(par->recover, geotop::common::Variables::files[rwcrn], sl->VS->wrain, top->rc_cont, par, land->LC, IT->LU);
        assign_recovered_map_vector(par->recover, geotop::common::Variables::files[rwcsn], sl->VS->wsnow, top->rc_cont, par, land->LC, IT->LU);
        assign_recovered_map_vector(par->recover, geotop::common::Variables::files[rTv], sl->VS->Tv, top->rc_cont, par, land->LC, IT->LU);
    }


    //	channel soil
    //	cnet->SS = (SOIL_STATE *)malloc(sizeof(SOIL_STATE));
    cnet->SS = new SoilState();
    //	initialize_soil_state(cnet->SS, cnet->r->nh, Nl);
    initialize_soil_state(cnet->SS, cnet->r.size(), geotop::common::Variables::Nl);

    //	cnet->th = new_doublematrix(geotop::common::Variables::Nl, cnet->r->nh);
    cnet->th.resize(geotop::common::Variables::Nl+1, cnet->r.size());

    //	cnet->ET = new_doublematrix(geotop::common::Variables::Nl, cnet->r->nh);
    //	initialize_doublematrix(cnet->ET, 0.0);
    cnet->ET.resize(geotop::common::Variables::Nl+1, cnet->r.size(),0.0);

    //	cnet->Kbottom = new_doublevector(cnet->r->nh);
    //	initialize_doublevector(cnet->Kbottom, 0.0);
    cnet->Kbottom.resize(cnet->r.size(),0.0);

    for(j=1;j<=par->total_channel;j++){

        //	sy=cnet->soil_type->co[j];
        sy=cnet->soil_type[j];
        //	r=cnet->r->co[j];
        r=cnet->r[j];
        //	c=cnet->c->co[j];
        c=cnet->c[j];

        //	cnet->SS->P->co[0][j] = sl->SS->P->co[0][top->j_cont[r][c]];
        cnet->SS->P[0][j] = sl->SS->P[0][top->j_cont[r][c]];
        for(l=1;l<=geotop::common::Variables::Nl;l++){
            //	cnet->SS->P->co[l][j] = sl->Ptot->co[l][top->j_cont[r][c]];
            cnet->SS->P[l][j] = sl->Ptot[l][top->j_cont[r][c]];
        }

        for(l=1;l<=geotop::common::Variables::Nl;l++){

            //	cnet->SS->T->co[l][j]=sl->pa->co[sy][jT][l];
            cnet->SS->T[l][j]=sl->pa(sy,jT,l);

            //	cnet->th->co[l][j] = teta_psi(cnet->SS->P->co[l][j], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
            //								  sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l],
            //								  GTConst::PsiMin, sl->pa->co[sy][jss][l]);
            cnet->th[l][j] = teta_psi(cnet->SS->P[l][j], 0.0, sl->pa(sy,jsat,l), sl->pa(sy,jres,l),
                                      sl->pa(sy,ja,l), sl->pa(sy,jns,l), 1.-1./sl->pa(sy,jns,l),
                                      GTConst::PsiMin, sl->pa(sy,jss,l));

            //	th_oversat = Fmax( cnet->SS->P->co[l][j] , 0.0 ) * sl->pa->co[sy][jss][l];
            th_oversat = Fmax( cnet->SS->P[l][j] , 0.0 ) * sl->pa(sy,jss,l);
            //	cnet->th->co[l][j] -= th_oversat;
            cnet->th[l][j] -= th_oversat;

            //	if(cnet->SS->T->co[l][j]<=GTConst::Tfreezing){
            if(cnet->SS->T[l][j]<=GTConst::Tfreezing){
                //	Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
                //	cnet->SS->thi->co[l][j] = cnet->th->co[l][j] - teta_psi(Psif(cnet->SS->T->co[l][j]), 0.0, sl->pa->co[sy][jsat][l],
                //															sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
                //															1.-1./sl->pa->co[sy][jns][l], GTConst::PsiMin, sl->pa->co[sy][jss][l]);
                cnet->SS->thi[l][j] = cnet->th[l][j] - teta_psi(Psif(cnet->SS->T[l][j]), 0.0, sl->pa(sy,jsat,l),
                                                                sl->pa(sy,jres,l), sl->pa(sy,ja,l), sl->pa(sy,jns,l),
                                                                1.-1./sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));

                //	if Theta(without freezing)<Theta_unfrozen(in equilibrium with T) Theta_ice is set at 0
                //	if(cnet->SS->thi->co[l][j]<0) cnet->SS->thi->co[l][j]=0.0;
                if(cnet->SS->thi[l][j]<0) cnet->SS->thi[l][j]=0.0;

                //Psi is updated taking into account the freezing
                //	cnet->th->co[l][j] -= cnet->SS->thi->co[l][j];
                cnet->th[l][j] -= cnet->SS->thi[l][j];
                //	cnet->SS->P->co[l][j] = psi_teta(cnet->th->co[l][j] + th_oversat, cnet->SS->thi->co[l][j],
                //									 sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
                //									 sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], GTConst::PsiMin, sl->pa->co[sy][jss][l]);
                cnet->SS->P[l][j] = psi_teta(cnet->th[l][j] + th_oversat, cnet->SS->thi[l][j],
                                             sl->pa(sy,jsat,l), sl->pa(sy,jres,l), sl->pa(sy,ja,l),
                                             sl->pa(sy,jns,l), 1.-1./sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));
            }
        }
    }

    if(par->delay_day_recover > 0 && par->total_channel > 0){
        assign_recovered_tensor_channel(0, par->recover, geotop::common::Variables::files[rpsich], cnet->SS->P, cnet->r, cnet->c, top->Z0);
        assign_recovered_tensor_channel(1, par->recover, geotop::common::Variables::files[ricegch], cnet->SS->thi, cnet->r, cnet->c, top->Z0);
        assign_recovered_tensor_channel(1, par->recover, geotop::common::Variables::files[rTgch], cnet->SS->T, cnet->r, cnet->c, top->Z0);

        for(i=1; i<=par->total_channel; i++){
            for (l=1; l<=geotop::common::Variables::Nl; l++) {
                //	sy = cnet->soil_type->co[i];
                sy = cnet->soil_type[i];
                //	cnet->th->co[l][i] = teta_psi(Fmin(cnet->SS->P->co[l][i], psi_saturation(cnet->SS->thi->co[l][i], sl->pa->co[sy][jsat][l],
                //			sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l])),
                //			cnet->SS->thi->co[l][i], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
                //			sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], GTConst::PsiMin, sl->pa->co[sy][jss][l]);
                cnet->th[l][i] = teta_psi(Fmin(cnet->SS->P[l][i], psi_saturation(cnet->SS->thi[l][i], sl->pa(sy,jsat,l),
                                                                                 sl->pa(sy,jres,l), sl->pa(sy,ja,l), sl->pa(sy,jns,l), 1.-1./sl->pa(sy,jns,l))),
                                          cnet->SS->thi[l][i], sl->pa(sy,jsat,l), sl->pa(sy,jres,l), sl->pa(sy,ja,l),
                                          sl->pa(sy,jns,l), 1.-1./sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));

            }
        }
    }

    //	WRITE INITIAL CONDITION
    //	write_output_headers(met->st->Z->nh, times, wat, par, top, land, sl, egy, snow, glac);
    write_output_headers(met->st->Z.size(), times, wat, par, top, land, sl, egy, snow, glac);

    if(par->state_pixel == 1){
        //	for(j=1;j<=par->rc->nrh;j++){
        for(j=1;j<par->rc.getRows();j++){
            if(par->output_vertical_distances == 1){
                //	r = par->rc->co[j][1];
                r = par->rc[j][1];
                //	c = par->rc->co[j][2];
                c = par->rc[j][2];
                //	cosslope = cos( Fmin(GTConst::max_slope, top->slope->co[r][c]) * GTConst::Pi/180.0 );
                cosslope = cos( Fmin(GTConst::max_slope, top->slope[r][c]) * GTConst::Pi/180.0 );
            }else {
                cosslope = 1.;
            }

            //write_soil_output(j, par->IDpoint->co[j], par->init_date->co[1], par->init_date->co[1], JD, day, month, year, hour, minute, par->soil_plot_depths, sl, par, GTConst::PsiMin, cosslope);
            write_soil_output(j, par->IDpoint[j], par->init_date[1], par->init_date[1], JD, day, month, year, hour, minute, par->soil_plot_depths, sl, par, GTConst::PsiMin, cosslope);
        }
    }

    //	z boundary condition
    for(l=1;l<=geotop::common::Variables::Nl;l++){
        //par->Zboundary -= sl->pa->co[1][jdz][l];
        par->Zboundary -= sl->pa(1,jdz,l);
    }

    if(par->Zboundary < 0){
        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
        fprintf(f, "Z at which 0 annual temperature takes place is not lower than the soil column\n");
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report (6).");
    }

    par->Zboundary *= 1.E-3;	//convert in [m]

    /****************************************************************************************************/
    /*! Initialization of the struct "egy" (of the type ENERGY):*/
    // revision performed on 24.12.2013// 
	
    if(par->output_surfenergy_bin == 1){
        if(geotop::common::Variables::files[fradnet] != geotop::input::gStringNoValue){
		    egy->Rn_mean.resize(par->total_pixel+1,0.0);
            egy->Rn.resize(par->total_pixel+1,0.0);
        }
        if(geotop::common::Variables::files[fradLWin] != geotop::input::gStringNoValue){
            egy->LWin_mean.resize(par->total_pixel+1, 0.0);
            egy->LWin.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fradLW] != geotop::input::gStringNoValue){
            egy->LW_mean.resize(par->total_pixel+1,0.0);
            egy->LW.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fradSW] != geotop::input::gStringNoValue){
            egy->SW_mean.resize(par->total_pixel+1,0.0);
            egy->SW.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fLE] != geotop::input::gStringNoValue){
            egy->ET_mean.resize(par->total_pixel+1,0.0);
            egy->LE.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fH] != geotop::input::gStringNoValue){
            egy->H_mean.resize(par->total_pixel+1,0.0);
            egy->H.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fG] != geotop::input::gStringNoValue){
            egy->SEB_mean.resize(par->total_pixel+1,0.0);
            egy->G.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fTs] != geotop::input::gStringNoValue){
            egy->Ts_mean.resize(par->total_pixel+1,0.0);
            egy->Ts.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fradSWin] != geotop::input::gStringNoValue){
            egy->Rswdown_mean.resize(par->total_pixel+1,0.0);
            egy->SWin.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fradSWinbeam] != geotop::input::gStringNoValue){
            egy->Rswbeam_mean.resize(par->total_pixel+1,0.0);
            egy->SWinb.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fshadow] != geotop::input::gStringNoValue){
            egy->nDt_shadow.resize(par->total_pixel+1,0.0);
            egy->nDt_sun.resize(par->total_pixel+1,0.0);
            egy->shad.resize(par->total_pixel+1,0.0);
        }
    }

    egy->sun = (double*)malloc(12*sizeof(double));


    //	if(times->JD_plots->nh > 1){
    if(times->JD_plots.size() > 1){
        if(geotop::common::Variables::files[pH] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pHg] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pG] != geotop::input::gStringNoValue){
            //	egy->Hgplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->Hgplot, 0.);
            egy->Hgplot.resize(par->total_pixel+1,0.0);
            //	egy->Hgp=new_doublevector(par->total_pixel);
            egy->Hgp.resize(par->total_pixel+1);

        }
        if(geotop::common::Variables::files[pH] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pHv] != geotop::input::gStringNoValue){
            //	egy->Hvplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->Hvplot, 0.);
            egy->Hvplot.resize(par->total_pixel+1,0.0);
            //	egy->Hvp=new_doublevector(par->total_pixel);
            egy->Hvp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pLE] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pLEg] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pG] != geotop::input::gStringNoValue){
            //	egy->LEgplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->LEgplot, 0.);
            egy->LEgplot.resize(par->total_pixel+1,0.0);
            //	egy->LEgp=new_doublevector(par->total_pixel);
            egy->LEgp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pLE] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pLEv] != geotop::input::gStringNoValue){
            //	egy->LEvplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->LEvplot, 0.);
            egy->LEvplot.resize(par->total_pixel+1,0.0);
            //	egy->LEvp=new_doublevector(par->total_pixel);
            egy->LEvp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pSWin] != geotop::input::gStringNoValue){
            //	egy->SWinplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->SWinplot, 0.);
            egy->SWinplot.resize(par->total_pixel+1,0.0);
            //	egy->SWinp=new_doublevector(par->total_pixel);
            egy->SWinp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pSWg] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pG] != geotop::input::gStringNoValue){
            //	egy->SWgplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->SWgplot, 0.);
            //	egy->SWgp=new_doublevector(par->total_pixel);
            egy->SWgplot.resize(par->total_pixel+1,0.0);
            egy->SWgp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pSWv] != geotop::input::gStringNoValue){
            //	egy->SWvplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->SWvplot, 0.);
            egy->SWvplot.resize(par->total_pixel+1,0.0);
            //	egy->SWvp=new_doublevector(par->total_pixel);
            egy->SWvp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pLWin] != geotop::input::gStringNoValue){
            //	egy->LWinplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->LWinplot, 0.);
            egy->LWinplot.resize(par->total_pixel+1,0.0);
            //	egy->LWinp=new_doublevector(par->total_pixel);
            egy->LWinp.resize(par->total_pixel);
        }
        if(geotop::common::Variables::files[pLWg] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pG] != geotop::input::gStringNoValue){
            //	egy->LWgplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->LWgplot, 0.);
            egy->LWgplot.resize(par->total_pixel+1,0.0);
            //	egy->LWgp=new_doublevector(par->total_pixel);
            egy->LWgp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pLWv] != geotop::input::gStringNoValue){
            //	egy->LWvplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->LWvplot, 0.);
            egy->LWvplot.resize(par->total_pixel+1,0.0);
            //	egy->LWvp=new_doublevector(par->total_pixel);
            egy->LWvp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pTs] != geotop::input::gStringNoValue){
            //	egy->Tsplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->Tsplot, 0.);
            egy->Tsplot.resize(par->total_pixel+1,0.0);
            //	egy->Tsp=new_doublevector(par->total_pixel);
            egy->Tsp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pTg] != geotop::input::gStringNoValue){
            //	egy->Tgplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->Tgplot, 0.);
            egy->Tgplot.resize(par->total_pixel+1,0.0);
            //	egy->Tgp=new_doublevector(par->total_pixel);
            egy->Tgp.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[pTv] != geotop::input::gStringNoValue){
            //	egy->Tvplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(egy->Tvplot, 0.);
            egy->Tvplot.resize(par->total_pixel+1,0.0);
        }
    }


    //vectors used in energy_balance()
	
    egy->Tgskin_surr.resize(geotop::common::Variables::Nr+1, geotop::common::Variables::Nc+1, 0.0);
    egy->SWrefl_surr.resize(geotop::common::Variables::Nr+1, geotop::common::Variables::Nc+1, 0.0);
    egy->Dlayer.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers +1);
    egy->liq.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers +1);
    egy->ice.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers +1);
    egy->Temp.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers +1);
    egy->deltaw.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers +1);
    egy->SWlayer.resize(par->max_snow_layers+1+1);

    //  tolto +1 dalla linea qua sotto 24.12.2013 S.C.&S.E. 
    egy->soil_transp_layer.resize(land->root_fraction.getCols());

    //	egy->dFenergy = new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
    egy->dFenergy.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers+1);
    //	egy->udFenergy = new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers - 1);
    egy->udFenergy.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers +1);
    //	egy->Kth0=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers - 1 );
    egy->Kth0.resize(geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers+1);
    //	egy->Kth1=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers - 1);
    egy->Kth1.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers+1);
    //	egy->Fenergy=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
    egy->Fenergy.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers +1);
    //	egy->Newton_dir=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
    egy->Newton_dir.resize(geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers+1);
    //	egy->T0=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
    egy->T0.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers+1);
    //	egy->T1=new_doublevector0( Nl + par->max_snow_layers + par->max_glac_layers );
    egy->T1.resize( geotop::common::Variables::Nl + par->max_snow_layers + par->max_glac_layers+1);
    //	egy->Tstar=new_doublevector(geotop::common::Variables::Nl); //soil temperature at which freezing begins
    egy->Tstar.resize(geotop::common::Variables::Nl+1);
    //	egy->THETA=new_doublevector(geotop::common::Variables::Nl);	//water content (updated in the iterations)
    egy->THETA.resize(geotop::common::Variables::Nl+1);

    //allocate vector	of soil layer contributions to evaporation (up to z_evap)
    z = 0.;
    l = 0;
    do{
        l++;
        //z += sl->pa->co[1][jdz][l];
        z += sl->pa(1,jdz,l);
    }while(l<geotop::common::Variables::Nl && z < GTConst::z_evap);
	
    //	egy->soil_evap_layer_bare = new_doublevector(l);
    //	initialize_doublevector(egy->soil_evap_layer_bare, 0.);
    egy->soil_evap_layer_bare.resize(l+1,0.0);
    //	egy->soil_evap_layer_veg = new_doublevector(l);
    //	initialize_doublevector(egy->soil_evap_layer_veg, 0.);
    egy->soil_evap_layer_veg.resize(l+1,0.0);

    //	printf("Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare->nh);
    printf("Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare.size()-1);
    //	printf("Soil water transpires from the first %ld layers\n",egy->soil_transp_layer->nh);
    printf("Soil water transpires from the first %ld layers\n",egy->soil_transp_layer.size()-1);
    //	fprintf(flog,"Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare->nh);
    fprintf(flog,"Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare.size()-1);
    //	fprintf(flog,"Soil water transpires from the first %ld layers\n",egy->soil_transp_layer->nh);
    fprintf(flog,"Soil water transpires from the first %ld layers\n",egy->soil_transp_layer.size()-1);

    /****************************************************************************************************/
    /*! Completing of the struct "water" (of the type WATER) */

    wat->Voutlandsub = 0.;
    wat->Voutlandsup = 0.;
    wat->Voutbottom = 0.;

    /* Initialization of wat->Pnet (liquid precipitation that reaches the sl surface in mm):*/
    //	wat->Pnet=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(wat->Pnet,0.0);
    wat->Pnet.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.0);

    /* Initialization of wat->PrecTot (total precipitation (rain+snow) precipitation):*/
    //	wat->PrecTot=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(wat->PrecTot,0.0);
    wat->PrecTot.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.0);

    /* Initialization of the matrices with the output of total precipitation and interception:*/
    if (par->output_meteo_bin == 1 &&geotop::common::Variables::files[fprec] != geotop::input::gStringNoValue){
        //	wat->PrTOT_mean=new_doublevector(par->total_pixel);
        //	initialize_doublevector(wat->PrTOT_mean, 0.);
        wat->PrTOT_mean.resize(par->total_pixel+1,0.0);

        //	wat->PrSNW_mean=new_doublevector(par->total_pixel);
        //	initialize_doublevector(wat->PrSNW_mean, 0.);
        wat->PrSNW_mean.resize(par->total_pixel+1,0.0);

        //	wat->Pt=new_doublevector(par->total_pixel);
        wat->Pt.resize(par->total_pixel+1);
        //	wat->Ps=new_doublevector(par->total_pixel);
        wat->Ps.resize(par->total_pixel+1);
    }

    //	wat->h_sup=new_doublevector(par->total_pixel);
    //	initialize_doublevector(wat->h_sup, 0.);
    wat->h_sup.resize(par->total_pixel+1,0.0);

    /****************************************************************************************************/
    /*! Initialization of the struct "snow" (of the type SNOW):*/

    /***************************************************************************************************/
    //snow->S=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
    snow->S=new Statevar3D();
    allocate_and_initialize_statevar_3D(snow->S, geotop::input::gDoubleNoValue, par->max_snow_layers, geotop::common::Variables::Nr, geotop::common::Variables::Nc);

    //initial snow depth
    if(geotop::common::Variables::files[fsn0] != geotop::input::gStringNoValue &&geotop::common::Variables::files[fswe0] != geotop::input::gStringNoValue ){
        printf("Initial condition on snow depth from file %s\n",geotop::common::Variables::files[fsn0].c_str());
        //	M=read_map(2,geotop::common::Variables::files[fsn0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fsn0]), M);

        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	snow->S->Dzl->co[1][r][c] = M->co[r][c];
                snow->S->Dzl[1][r][c] = M[r][c];
            }
        }
        //	free_doublematrix(M);

        printf("Initial condition on snow water equivalent from file %s\n",geotop::common::Variables::files[fswe0].c_str());
        //	M=read_map(2,geotop::common::Variables::files[fswe0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fswe0]), M);

        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	snow->S->w_ice->co[1][r][c] = M->co[r][c];
                snow->S->w_ice[1][r][c] = M[r][c];
            }
        }
        //	free_doublematrix(M);

    }else if(geotop::common::Variables::files[fsn0] != geotop::input::gStringNoValue ){
        printf("Initial condition on snow depth from file %s\n",geotop::common::Variables::files[fsn0].c_str());
        //	M=read_map(2,geotop::common::Variables::files[fsn0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fsn0]), M);

        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	snow->S->Dzl->co[1][r][c] = M->co[r][c];
                snow->S->Dzl[1][r][c] = M[r][c];
            }
        }
        //	free_doublematrix(M);

        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	if ((long)land->LC->co[r][c] != geotop::input::gDoubleNoValue) snow->S->w_ice->co[1][r][c] = snow->S->Dzl->co[1][r][c]*IT->rhosnow0/GTConst::rho_w;
                if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue) snow->S->w_ice[1][r][c] = snow->S->Dzl[1][r][c]*IT->rhosnow0/GTConst::rho_w;
            }
        }

    }else if(geotop::common::Variables::files[fswe0] != geotop::input::gStringNoValue ){
        printf("Initial condition on snow water equivalent from file %s\n",geotop::common::Variables::files[fswe0].c_str());
        //	M=read_map(2,geotop::common::Variables::files[fswe0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fswe0]), M);
        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	snow->S->w_ice->co[1][r][c] = M->co[r][c];
                snow->S->w_ice[1][r][c] = M[r][c];
            }
        }
        //	free_doublematrix(M);

        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	if ((long)land->LC->co[r][c] != geotop::input::gDoubleNoValue) snow->S->Dzl->co[1][r][c] = snow->S->w_ice->co[1][r][c]*GTConst::rho_w/IT->rhosnow0;
                if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue) snow->S->Dzl[1][r][c] = snow->S->w_ice[1][r][c]*GTConst::rho_w/IT->rhosnow0;
            }
        }

    }else {

        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	if ((long)land->LC->co[r][c] != geotop::input::gDoubleNoValue){
                if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue){
                    //	snow->S->w_ice->co[1][r][c] = IT->swe0;
                    snow->S->w_ice[1][r][c] = IT->swe0;
                    snow->S->Dzl[1][r][c] = IT->swe0*GTConst::rho_w/IT->rhosnow0;
                }
            }
        }
    }


    //Optional reading of snow age in the whole basin
    if(geotop::common::Variables::files[fsnag0] != geotop::input::gStringNoValue ){
        printf("Snow age initial condition from file %s\n",geotop::common::Variables::files[fsnag0 + 1].c_str());
        snow->age = read_map_vector(2,geotop::common::Variables::files[fsnag0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue, top->rc_cont);
    }else{
        //	snow->age = new_doublevector(par->total_pixel);
        //	initialize_doublevector(snow->age, IT->agesnow0);
        snow->age.resize(par->total_pixel+1,IT->agesnow0);
    }


    //	if(times->JD_plots->nh > 1){
    if(times->JD_plots.size() > 1){
        if(geotop::common::Variables::files[pD] != geotop::input::gStringNoValue) {
            //	snow->Dplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(snow->Dplot, 0.);
            snow->Dplot.resize(par->total_pixel+1,0.0);
        }
    }

    if(par->blowing_snow==1){

        //	snow->S_for_BS=(STATEVAR_1D *)malloc(sizeof(STATEVAR_1D));
        snow->S_for_BS=new Statevar1D();
        allocate_and_initialize_statevar_1D(snow->S_for_BS, geotop::input::gDoubleNoValue, par->max_snow_layers);

        //	snow->change_dir_wind=new_longvector(Fmaxlong(Nr,Nc));
        snow->change_dir_wind.resize(Fmaxlong(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1));

        //	snow->Qtrans=new_doublematrix(Nr,Nc);
        //	initialize_doublematrix(snow->Qtrans, 0.0);
        snow->Qtrans.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.0);

        //	snow->Qsub=new_doublematrix(Nr,Nc);
        //	initialize_doublematrix(snow->Qsub, 0.0);
        snow->Qsub.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.0);

        //	snow->Qsalt=new_doublematrix(Nr,Nc);
        snow->Qsalt.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);

        //	snow->Nabla2_Qtrans=new_doublematrix(Nr,Nc);
        snow->Nabla2_Qtrans.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);

        //	snow->Qsub_x=new_doublematrix(Nr,Nc);
        snow->Qsub_x.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);

        //	snow->Qsub_y=new_doublematrix(Nr,Nc);
        snow->Qsub_y.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);

        //	snow->Qtrans_x=new_doublematrix(Nr,Nc);
        snow->Qtrans_x.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);

        //	snow->Qtrans_y=new_doublematrix(Nr,Nc);
        snow->Qtrans_y.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);

        if(par->output_snow_bin == 1){
            //	snow->Wtrans_plot=new_doublematrix(Nr,Nc);
            snow->Wtrans_plot.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);
            //	snow->Wsubl_plot=new_doublematrix(Nr,Nc);
            snow->Wsubl_plot.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1);

            for(r=1;r<=geotop::common::Variables::Nr;r++){
                for(c=1;c<=geotop::common::Variables::Nc;c++){
                    //	if((long)land->LC->co[r][c]==geotop::input::gDoubleNoValue){
                    if((long)land->LC[r][c]==geotop::input::gDoubleNoValue){
                        //	snow->Wtrans_plot->co[r][c]=geotop::input::gDoubleNoValue;
                        snow->Wtrans_plot[r][c]=geotop::input::gDoubleNoValue;
                        //	snow->Wsubl_plot->co[r][c]=geotop::input::gDoubleNoValue;
                        snow->Wsubl_plot[r][c]=geotop::input::gDoubleNoValue;

                    }else{
                        //	snow->Wtrans_plot->co[r][c]=0.0;
                        snow->Wtrans_plot[r][c]=0.0;
                        //	snow->Wsubl_plot->co[r][c]=0.0;
                        snow->Wsubl_plot[r][c]=0.0;

                    }
                }
            }
        }
    }

    if(par->output_snow_bin == 1){
        if(geotop::common::Variables::files[fsnowmelt] != geotop::input::gStringNoValue){
            //	snow->MELTED=new_doublevector(par->total_pixel);
            //	initialize_doublevector(snow->MELTED,0.);
            snow->MELTED.resize(par->total_pixel+1,0.0);
            //	snow->melted=new_doublevector(par->total_pixel);
            snow->melted.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fsnowsubl] != geotop::input::gStringNoValue){
            //	snow->SUBL=new_doublevector(par->total_pixel);
            //	initialize_doublevector(snow->SUBL,0.);
            snow->SUBL.resize(par->total_pixel+1,0.0);
            //	snow->subl=new_doublevector(par->total_pixel);
            snow->subl.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fsndur] != geotop::input::gStringNoValue){
            //	snow->t_snow=new_doublevector(par->total_pixel);
            //	initialize_doublevector(snow->t_snow,0.);
            snow->t_snow.resize(par->total_pixel+1,0.0);
            //	snow->yes=new_shortvector(par->total_pixel);
            snow->yes.resize(par->total_pixel+1,0.0);
        }
    }

    for(r=1;r<=geotop::common::Variables::Nr;r++){
        for(c=1;c<=geotop::common::Variables::Nc;c++){

            //	if( (long)land->LC->co[r][c]!=geotop::input::gDoubleNoValue){
            if( (long)land->LC[r][c]!=geotop::input::gDoubleNoValue){

                //Adjusting snow init depth in case of steep slope (contribution by Stephan Gruber)
                //	if (par->snow_curv > 0 && top->slope->co[r][c] > par->snow_smin){
                if (par->snow_curv > 0 && top->slope[r][c] > par->snow_smin){
                    //	if (top->slope->co[r][c] <= par->snow_smax){
                    if (top->slope[r][c] <= par->snow_smax){
                        //	k_snowred = ( exp(-pow(top->slope->co[r][c] - par->snow_smin, 2.)/par->snow_curv) -
                        //				 exp(-pow(par->snow_smax, 2.)/par->snow_curv) );
                        k_snowred = ( exp(-pow(top->slope[r][c] - par->snow_smin, 2.)/par->snow_curv) -
                                      exp(-pow(par->snow_smax, 2.)/par->snow_curv) );
                    }else{
                        k_snowred = 0.0;
                    }
                    //	snow->S->Dzl->co[1][r][c] *= k_snowred;
                    snow->S->Dzl[1][r][c] *= k_snowred;
                    //	snow->S->w_ice->co[1][r][c] *= k_snowred;
                    snow->S->w_ice[1][r][c] *= k_snowred;
                }

                //	D = snow->S->Dzl->co[1][r][c];
                D = snow->S->Dzl[1][r][c];
                //	SWE = snow->S->w_ice->co[1][r][c];
                SWE = snow->S->w_ice[1][r][c];

                if(D<0 || SWE<0){
                    f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
#ifdef USE_DOUBLE_PRECISION_OUTPUT
                    fprintf(f, "Error: negative initial snow depth %12g or snow water equivalent %12g\n",D,SWE);
#else
                    fprintf(f, "Error: negative initial snow depth %e or snow water equivalent %e\n",D,SWE);
#endif

                    fclose(f);
                    t_error("Fatal Error! Geotop is closed. See failing report (7).");

                }else if(D<1.E-5 && SWE>1.E-5){
                    f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
#ifdef USE_DOUBLE_PRECISION_OUTPUT
                    fprintf(f, "Error: Initial snow water equivalent %12g > 0 and initial snow depth %12g\n",SWE,D);
#else
                    fprintf(f, "Error: Initial snow water equivalent %e > 0 and initial snow depth %e\n",SWE,D);
#endif

                    fclose(f);
                    t_error("Fatal Error! Geotop is closed. See failing report (8).");

                }else if(D>1.E-5 && SWE<1.E-5){
                    f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
#ifdef USE_DOUBLE_PRECISION_OUTPUT
                    fprintf(f, "Error: Initial snow depth %12g > 0 and initial snow water equivalent %12g\n",D,SWE);
#else
                    fprintf(f, "Error: Initial snow depth %e > 0 and initial snow water equivalent %e\n",D,SWE);
#endif
                    fclose(f);
                    t_error("Fatal Error! Geotop is closed. See failing report (9).");

                }else if(D>1.E-5 || SWE>1.E-5){

                    //	snow->age->co[top->j_cont[r][c]]*=86400.0;	//now in [s]
                    snow->age[top->j_cont[r][c]]*=86400.0;	//now in [s]
                    if (SWE <= par->max_weq_snow * par->max_snow_layers ) {

                        i = floor( SWE/par->max_weq_snow );

                        if (i>0) {

                            for (n=1; n<=i; n++) {
                                //	snow->S->w_ice->co[n][r][c] = par->max_weq_snow;
                                snow->S->w_ice[n][r][c] = par->max_weq_snow;
                            }

                            if (SWE - i * par->max_weq_snow > 0.1 * par->max_weq_snow) {
                                //	snow->S->w_ice->co[i+1][r][c] = SWE - i * par->max_weq_snow;
                                snow->S->w_ice[i+1][r][c] = SWE - i * par->max_weq_snow;
                                //	snow->S->lnum->co[r][c] = i+1;
                                snow->S->lnum[r][c] = i+1;
                            }else {
                                //	snow->S->w_ice->co[i][r][c] += (SWE - i * par->max_weq_snow);
                                snow->S->w_ice[i][r][c] += (SWE - i * par->max_weq_snow);
                                //	snow->S->lnum->co[r][c] = i;
                                snow->S->lnum[r][c] = i;
                            }

                        }else {

                            //	snow->S->w_ice->co[1][r][c] = SWE;
                            snow->S->w_ice[1][r][c] = SWE;
                            //	snow->S->lnum->co[r][c] = 1;
                            snow->S->lnum[r][c] = 1;

                        }

                    }else {

                        //	snow->S->lnum->co[r][c] = par->max_snow_layers;
                        snow->S->lnum[r][c] = par->max_snow_layers;

                        for (n=1; n<=par->max_snow_layers; n++) {

                            a = 0;

                            //	for (i=1; i<=par->inf_snow_layers->nh; i++) {
                            for (size_t i=1; i<par->inf_snow_layers.size(); i++) {
                                //	if (n == abs(par->inf_snow_layers->co[i])) a = 1;
                                if (n == abs(par->inf_snow_layers[i])) a = 1;
                            }

                            if (a == 0) {
                                //	snow->S->w_ice->co[n][r][c] = par->max_weq_snow;
                                snow->S->w_ice[n][r][c] = par->max_weq_snow;
                            }else {
                                //	snow->S->w_ice->co[n][r][c] = ( SWE - par->max_weq_snow * ( par->max_snow_layers - par->inf_snow_layers->nh ) ) / par->inf_snow_layers->nh;
                                snow->S->w_ice[n][r][c] = ( SWE - par->max_weq_snow * ( par->max_snow_layers - par->inf_snow_layers.size() ) ) / par->inf_snow_layers.size();
                            }
                        }
                    }

                    //	for (n=1; n<=snow->S->lnum->co[r][c]; n++) {
                    for (n=1; n<=snow->S->lnum[r][c]; n++) {
                        //	snow->S->Dzl->co[n][r][c] = D * (snow->S->w_ice->co[n][r][c] / SWE);
                        snow->S->Dzl[n][r][c] = D * (snow->S->w_ice[n][r][c] / SWE);
                        //	snow->S->T->co[n][r][c] = IT->Tsnow0;
                        snow->S->T[n][r][c] = IT->Tsnow0;
                    }

                }

                //				non_dimensionalize_snowage(&(snow->age->co[top->j_cont[r][c]]), IT->Tsnow0);
                non_dimensionalize_snowage(&(snow->age[top->j_cont[r][c]]), IT->Tsnow0);

                if (par->point_sim == 1) {
                    //	maxSWE = par->maxSWE->co[r][c];
                    maxSWE = par->maxSWE[r][c];
                }else {
                    maxSWE = 1.E10;
                }

                f = fopen(geotop::common::Variables::logfile.c_str(), "a");
                snow_layer_combination(par->alpha_snow, r, c, snow->S, -0.1, par->inf_snow_layers, par->max_weq_snow, maxSWE, f);
                fclose(f);

            }
        }
    }

    if(par->delay_day_recover > 0){
        //	initialize_shortmatrix(snow->S->type, 2);
        snow->S->type.resize(snow->S->type.getRows(),snow->S->type.getCols(),2);

        assign_recovered_map_long(par->recover, geotop::common::Variables::files[rns], snow->S->lnum, par, land->LC, IT->LU);
        assign_recovered_map_vector(par->recover, geotop::common::Variables::files[rsnag], snow->age, top->rc_cont, par, land->LC, IT->LU);
        assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rDzs], snow->S->Dzl, par, land->LC, IT->LU);
        assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rwls], snow->S->w_liq, par, land->LC, IT->LU);
        assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rwis], snow->S->w_ice, par, land->LC, IT->LU);
        assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rTs], snow->S->T, par, land->LC, IT->LU);
    }



    /****************************************************************************************************/
    /*! Initialization of the struct "glac" (of the type GLACIER):*/

    /***************************************************************************************************/
    /*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
    if( par->point_sim!=1 &&geotop::common::Variables::files[fgl0] != geotop::input::gStringNoValue && par->max_glac_layers==0){
        printf("Warning: Glacier map present, but glacier represented with 0 layers\n");
        fprintf(flog,"Warning: Glacier map present, but glacier represented with 0 layers\n");
    }

    if(par->max_glac_layers==0 && IT->Dglac0>0){
        printf("\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
        fprintf(flog,"\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
    }

    //If the max number of glacier layers is greater than 1, the matrices (or tensors) lnum, Dzl. w_liq, w_ice, T and print matrices are defined, according to the respective flags
    if(par->max_glac_layers>0){

        if( par->point_sim!=1 &&geotop::common::Variables::files[fgl0] != geotop::input::gStringNoValue ){
            printf("Glacier initial condition from file %s\n",geotop::common::Variables::files[fgl0+1].c_str());
            fprintf(flog,"Glacier initial condition from file %s\n",geotop::common::Variables::files[fgl0+1].c_str());
            //M=read_map(2,geotop::common::Variables::files[fgl0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            meteoio_readMap(string(geotop::common::Variables::files[fgl0]), M);
        }else{
            //	M=copydoublematrix_const(IT->Dglac0, land->LC, geotop::input::gDoubleNoValue);
            copydoublematrix_const(IT->Dglac0, land->LC, M, geotop::input::gDoubleNoValue);
        }

        //	glac->G=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
        glac->G=new Statevar3D();
        allocate_and_initialize_statevar_3D(glac->G, geotop::input::gDoubleNoValue, par->max_glac_layers, geotop::common::Variables::Nr, geotop::common::Variables::Nc);

        if(par->output_glac_bin == 1){
            if(geotop::common::Variables::files[fglacmelt] != geotop::input::gStringNoValue){
                //	glac->MELTED=new_doublevector(par->total_pixel);
                //	initialize_doublevector(glac->MELTED,0.);
                glac->MELTED.resize(par->total_pixel+1,0.0);
                //	glac->melted=new_doublevector(par->total_pixel);
                glac->melted.resize(par->total_pixel+1);
            }
            if(geotop::common::Variables::files[fglacsubl] != geotop::input::gStringNoValue){
                //	glac->SUBL=new_doublevector(par->total_pixel);
                //	initialize_doublevector(glac->SUBL, 0.);
                glac->SUBL.resize(par->total_pixel+1,0.0);
                //	glac->subl=new_doublevector(par->total_pixel);
                glac->subl.resize(par->total_pixel+1);
            }
        }

        for(r=1;r<=geotop::common::Variables::Nr;r++){
            for(c=1;c<=geotop::common::Variables::Nc;c++){
                //	if( (long)land->LC->co[r][c]!=geotop::input::gDoubleNoValue){
                if( (long)land->LC[r][c]!=geotop::input::gDoubleNoValue){

                    //	if(M->co[r][c]<0){
                    if(M[r][c]<0){
                        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                        fprintf(f, "Error: negative glacier data\n");
                        fclose(f);
                        t_error("Fatal Error! Geotop is closed. See failing report (10).");

                        //	}else if(M->co[r][c]>1.E-5){
                    }else if(M[r][c]>1.E-5){

                        //	if (IT->rhoglac0 * M->co[r][c] / GTConst::rho_w < par->max_weq_glac * par->max_glac_layers ) {
                        if (IT->rhoglac0 * M[r][c] / GTConst::rho_w < par->max_weq_glac * par->max_glac_layers ) {

                            n = 0;
                            z = 0.;

                            do{

                                n++;

                                //	if (IT->rhoglac0 * M->co[r][c] / GTConst::rho_w < par->max_weq_glac * n) {
                                if (IT->rhoglac0 * M[r][c] / GTConst::rho_w < par->max_weq_glac * n) {
                                    //	glac->G->w_ice->co[n][r][c] = par->max_weq_glac;
                                    glac->G->w_ice[n][r][c] = par->max_weq_glac;
                                }else {
                                    //	glac->G->w_ice->co[n][r][c] = IT->rhoglac0 * M->co[r][c] / 1000. - z;
                                    glac->G->w_ice[n][r][c] = IT->rhoglac0 * M[r][c] / 1000. - z;
                                }

                                //	glac->G->Dzl->co[n][r][c] = GTConst::rho_w * glac->G->w_ice->co[n][r][c] / IT->rhoglac0;
                                glac->G->Dzl[n][r][c] = GTConst::rho_w * glac->G->w_ice[n][r][c] / IT->rhoglac0;
                                //	glac->G->T->co[n][r][c] = IT->Tglac0;
                                glac->G->T[n][r][c] = IT->Tglac0;

                                //	z += glac->G->w_ice->co[n][r][c];
                                z += glac->G->w_ice[n][r][c];

                                //	}while (fabs(z - IT->rhoglac0 * M->co[r][c] / GTConst::rho_w) < 1.E-6);
                            }while (fabs(z - IT->rhoglac0 * M[r][c] / GTConst::rho_w) < 1.E-6);

                            //	glac->G->lnum->co[r][c] = n;
                            glac->G->lnum[r][c] = n;

                        }else {

                            //	glac->G->lnum->co[r][c] = par->max_glac_layers;
                            glac->G->lnum[r][c] = par->max_glac_layers;

                            for (n=1; n<=par->max_glac_layers; n++) {

                                a = 0;

                                //	for (i=1; i<=par->inf_glac_layers->nh; i++) {
                                for (i=1; i< par->inf_glac_layers.size(); i++) {
                                    if (n == i) a = 1;
                                }

                                if (a == 0) {
                                    //	glac->G->w_ice->co[n][r][c] = par->max_weq_glac;
                                    glac->G->w_ice[n][r][c] = par->max_weq_glac;
                                }else {
                                    //	glac->G->w_ice->co[n][r][c] = ( IT->rhoglac0 * M->co[r][c] / GTConst::rho_w - par->max_weq_glac * ( par->max_glac_layers - par->inf_glac_layers->nh ) ) / par->inf_glac_layers->nh;
                                    glac->G->w_ice[n][r][c] = ( IT->rhoglac0 * M[r][c] / GTConst::rho_w - par->max_weq_glac * ( par->max_glac_layers - par->inf_glac_layers.size() ) ) / par->inf_glac_layers.size();
                                }

                                //	glac->G->Dzl->co[n][r][c] = GTConst::rho_w * glac->G->w_ice->co[n][r][c] / IT->rhoglac0;
                                glac->G->Dzl[n][r][c] = GTConst::rho_w * glac->G->w_ice[n][r][c] / IT->rhoglac0;
                                //	glac->G->T->co[n][r][c] = IT->Tglac0;
                                glac->G->T[n][r][c] = IT->Tglac0;

                            }
                        }
                    }

                    f = fopen(geotop::common::Variables::logfile.c_str(), "a");
                    snow_layer_combination(par->alpha_snow, r, c, glac->G, -0.1, par->inf_glac_layers, par->max_weq_glac, 1.E10, f);
                    fclose(f);

                }
            }
        }

        //	free_doublematrix(M);

        if(par->delay_day_recover > 0){
            assign_recovered_map_long(par->recover, geotop::common::Variables::files[rni], glac->G->lnum, par, land->LC, IT->LU);
            assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rDzi], glac->G->Dzl, par, land->LC, IT->LU);
            assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rwli], glac->G->w_liq, par, land->LC, IT->LU);
            assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rwii], glac->G->w_ice, par, land->LC, IT->LU);
            assign_recovered_tensor(1, par->recover, geotop::common::Variables::files[rTi], glac->G->T, par, land->LC, IT->LU);
        }
    }



    //***************************************************************************************************
    // Filling up of the struct "met" (of the type METEO):

    //	met->Tgrid=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(met->Tgrid, 5.);
    met->Tgrid.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,5.);

    //	met->Pgrid=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(met->Pgrid, GTConst::Pa0);
    met->Pgrid.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,GTConst::Pa0);

    //	met->RHgrid=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(met->RHgrid, 0.7);
    met->RHgrid.resize(geotop::common::Variables::Nr+1, geotop::common::Variables::Nc+1, 0.7);

    //	met->Vgrid=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(met->Vgrid, par->Vmin);
    met->Vgrid.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1, par->Vmin);

    //	met->Vdir=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(met->Vdir, 0.0);
    met->Vdir.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.0);

    if (par->output_meteo_bin == 1){
        if(geotop::common::Variables::files[fTa] != geotop::input::gStringNoValue){
            //	met->Tamean=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->Tamean, 0.);
            met->Tamean.resize(par->total_pixel+1);
        }
        if(geotop::common::Variables::files[fwspd] != geotop::input::gStringNoValue){
            //	met->Vspdmean=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->Vspdmean, 0.);
            met->Vspdmean.resize(par->total_pixel+1,0.0);
        }
        if(geotop::common::Variables::files[fwdir] != geotop::input::gStringNoValue){
            //	met->Vdirmean=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->Vdirmean, 0.);
            met->Vdirmean.resize(par->total_pixel+1,0.0);
        }
        if(geotop::common::Variables::files[frh] != geotop::input::gStringNoValue){
            //	met->RHmean=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->RHmean, 0.);
            met->RHmean.resize(par->total_pixel+1,0.0);
        }
    }


    //plot output
    //	if(times->JD_plots->nh > 1){
    if(times->JD_plots.size() > 1){
        if(geotop::common::Variables::files[pTa] != geotop::input::gStringNoValue){
            //	met->Taplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->Taplot, 0.);
            met->Taplot.resize(par->total_pixel+1,0.0);
        }
        if(geotop::common::Variables::files[pRH] != geotop::input::gStringNoValue){
            //	met->RHplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->RHplot, 0.);
            met->RHplot.resize(par->total_pixel+1,0.0);
        }
        if(geotop::common::Variables::files[pVspd] != geotop::input::gStringNoValue ||geotop::common::Variables::files[pVdir] != geotop::input::gStringNoValue){
            //	met->Vxplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->Vxplot, 0.);
            met->Vxplot.resize(par->total_pixel+1,0.0);
            //	met->Vyplot=new_doublevector(par->total_pixel);
            //	initialize_doublevector(met->Vyplot, 0.);
            met->Vyplot.resize(par->total_pixel+1,0.0);
        }
    }

    /****************************************************************************************************/

    /*Free the struct allocated in this subroutine:*/
    //free_doublematrix(par->chkpt);
    //free_doublevector(sl->init_water_table_depth);

    delete(IT);
    //free(IT);

    if (par->point_sim != 1) {

        cont_nonzero_values_matrix2(&i, &j, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel);
		
        //	top->Li = new_longvector(i);
        top->Li.resize(i+1);
        //	top->Lp = new_longvector(j);
        top->Lp.resize(j+1);

        //	wat->Lx = new_doublevector(i);
        wat->Lx.resize(i+1);
        cont_nonzero_values_matrix3(top->Lp, top->Li, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel);

    }else {

        i = geotop::common::Variables::Nl;
        j = geotop::common::Variables::Nl + 1;
        //	top->Li = new_longvector(i);
        top->Li.resize(i+1);
        //	top->Lp = new_longvector(j);
        top->Lp.resize(j+1);

        //	wat->Lx = new_doublevector(i);
        wat->Lx.resize(i+1);

        for (l=1; l<=geotop::common::Variables::Nl; l++) {
            //	top->Li->co[l] = l+1;
            top->Li[l] = l+1;
            //	top->Lp->co[l] = l;
            top->Lp[l] = l;
        }
        //	top->Lp->co[l] = i;
        top->Lp[l] = i;
    }

    //	wat->H0 = new_doublevector(j);
    wat->H0.resize(j+1);
    //	wat->H1 = new_doublevector(j);
    wat->H1.resize(j+1);
    //	wat->dH = new_doublevector(j);
    wat->dH.resize(j+1);
    //	wat->B = new_doublevector(j);
    wat->B.resize(j+1);
    //	wat->f = new_doublevector(j);
    wat->f.resize(j+1);
    //	wat->df = new_doublevector(j);
    wat->df.resize(j+1);
    //	wat->Kbottom = new_doublematrix(Nr, Nc);
    //	initialize_doublematrix(wat->Kbottom, 0.);
    wat->Kbottom.resize(geotop::common::Variables::Nr+1, geotop::common::Variables::Nc+1, 0.);

    //	wat->Klat = new_doublematrix(top->BC_DepthFreeSurface->nh, Nl);
    //	initialize_doublematrix(wat->Klat, 0.);
    wat->Klat.resize(top->BC_DepthFreeSurface.size(), geotop::common::Variables::Nl+1,0.0);

    fclose(flog);

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//void read_inputmaps(TOPO *top, LAND *land, SOIL *sl, PAR *par, FILE *flog, mio::IOManager& iomanager){
void read_inputmaps(Topo *top, Land *land, Soil *sl, Par *par, FILE *flog, mio::IOManager& iomanager){

    long r, c, i, cont;
    //	DOUBLEMATRIX *M;
    GeoMatrix<double> M;
    //	SHORTMATRIX *curv;
    GeoMatrix<short> curv;
    short flag;
    //	char *temp;
    string temp;
    double min, max;
    FILE *f;

    meteoio_readDEM(top->Z0);

    //	reading TOPOGRAPHY
    flag = file_exists(fdem, flog);
    if(flag == 1){
        //	M=new_doublematrix(1,1);
        //	top->Z0=read_map(0,geotop::common::Variables::files[fdem], M, geotop::common::Variables::UV, geotop::input::gDoubleNoValue); //topography
        //	free_doublematrix(M);

        //write_map(geotop::common::Variables::files[fdem], 0, par->format_out, top->Z0, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        //	meteoio_writeEsriasciiMap(geotop::common::Variables::files[fdem], pUV, top->Z0 ,geotop::input::gDoubleNoValue) ;

        //filtering

        //	M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
        M.resize(top->Z0.getRows(),top->Z0.getCols());
        multipass_topofilter(par->lowpass, top->Z0, M, geotop::input::gDoubleNoValue, 1);
        //	copy_doublematrix(M, top->Z0);
        top->Z0=M;
        //	free_doublematrix(M);
        //	write_map(geotop::common::Variables::files[fdem], 0, par->format_out, top->Z0, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

        //	calculate East and North
        //	top->East = new_doublematrix(top->Z0->nrh, top->Z0->nch);
        top->East.resize(top->Z0.getRows(), top->Z0.getCols());
        //	top->North = new_doublematrix(top->Z0->nrh, top->Z0->nch);
        top->North.resize(top->Z0.getRows(), top->Z0.getCols());
        //	for (r=1; r<=top->Z0->nrh; r++) {
        for (r=1; r<top->Z0.getRows(); r++) {
            //	for (c=1; c<=top->Z0->nch; c++) {
            for (c=1; c<top->Z0.getCols(); c++) {
                //	top->East->co[r][c] = UV->U->co[4] + (c-0.5)*UV->U->co[2];
                top->East[r][c] = geotop::common::Variables::UV->U[4] + (c-0.5)*geotop::common::Variables::UV->U[2];
                //	top->North->co[r][c] = UV->U->co[3] + (top->Z0->nrh-(r-0.5))*UV->U->co[1];
                top->North[r][c] = geotop::common::Variables::UV->U[3] + ((top->Z0.getRows()-1)-(r-0.5))*geotop::common::Variables::UV->U[1];
            }
        }

    }else{

        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
        fprintf(f, "Error: It is impossible to proceed without giving the digital elevation model\n");
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report (11).");

    }

    //	reading LAND COVER TYPE
    flag = file_exists(flu, flog);
    if(flag == 1){
        //	land->LC=read_map(1,geotop::common::Variables::files[flu], top->Z0, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[flu]), land->LC);

        //	Check borders
        //	for(r=1;r<=land->LC->nrh;r++){
        for(r=1;r<land->LC.getRows();r++){
            //	land->LC->co[r][1]=geotop::input::gDoubleNoValue;
            land->LC[r][1]=geotop::input::gDoubleNoValue;
            //	land->LC->co[r][land->LC->nch]=geotop::input::gDoubleNoValue;
            land->LC[r][land->LC.getCols()-1]=geotop::input::gDoubleNoValue;
        }
        //	for(c=1;c<=land->LC->nch;c++){
        for(c=1;c<land->LC.getCols();c++){
            //	land->LC->co[1][c]=geotop::input::gDoubleNoValue;
            land->LC[1][c]=geotop::input::gDoubleNoValue;
            //	land->LC->co[land->LC->nrh][c]=geotop::input::gDoubleNoValue;
            land->LC[land->LC.getRows()-1][c]=geotop::input::gDoubleNoValue;
        }
        //	for(r=1;r<=land->LC->nrh;r++){
        for(r=1;r<land->LC.getRows();r++){
            //	for(c=1;c<=land->LC->nch;c++){
            for(c=1;c<land->LC.getCols();c++){
                //	if ((long)land->LC->co[r][c] != geotop::input::gDoubleNoValue) {
                if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue) {
                    //	if ((long)land->LC->co[r][c] < 1 || (long)land->LC->co[r][c] > par->n_landuses){
                    if ((long)land->LC[r][c] < 1 || (long)land->LC[r][c] > par->n_landuses){
                        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                        fprintf(f, "Error: It is not possible to assign Value < 1 or > n_landuses to the land cover type\n");
                        fclose(f);
                        t_error("Fatal Error! Geotop is closed. See failing report (12).");
                    }
                }
            }
        }

        //Land use is the official mask
        //	for(r=1;r<=land->LC->nrh;r++){
        for(r=1;r< land->LC.getRows();r++){
            //	for(c=1;c<=land->LC->nch;c++){
            for(c=1;c<land->LC.getCols();c++){
                //	if((long)land->LC->co[r][c]!=geotop::input::gDoubleNoValue){
                if((long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                    //	if((long)top->Z0->co[r][c]==geotop::input::gDoubleNoValue){
                    if((long)top->Z0[r][c]==geotop::input::gDoubleNoValue){
                        printf("ERROR Land use mask include DTM novalue pixels");
                        //	printf("\nr:%ld c:%ld Z:%f landuse:%f\n",r,c,top->Z0->co[r][c],land->LC->co[r][c]);
                        printf("\nr:%ld c:%ld Z:%f landuse:%f\n",r,c,top->Z0[r][c],land->LC[r][c]);
                        //	land->LC->co[r][c]=geotop::input::gDoubleNoValue;
                        land->LC[r][c]=geotop::input::gDoubleNoValue;
                        printf("LANDUSE set at novalue where DTM is not available\n");
                    }
                }
            }
        }

    }else{

        //	Write land->LC (land cover)
        printf("Land cover type assumed to be always 1\n");
        //	land->LC=copydoublematrix_const(1.0, top->Z0, geotop::input::gDoubleNoValue);
        copydoublematrix_const(1.0, top->Z0,land->LC ,geotop::input::gDoubleNoValue);

        //	for(r=1;r<=land->LC->nrh;r++){
        for(r=1;r<land->LC.getRows();r++){
            //	land->LC->co[r][1]=geotop::input::gDoubleNoValue;
            land->LC[r][1]=geotop::input::gDoubleNoValue;
            //	land->LC->co[r][land->LC->nch]=geotop::input::gDoubleNoValue;
            land->LC[r][land->LC.getCols()-1]=geotop::input::gDoubleNoValue;
        }
        //	for(c=1;c<=land->LC->nch;c++){
        for(c=1;c<land->LC.getCols();c++){
            //	land->LC->co[1][c]=geotop::input::gDoubleNoValue;
            land->LC[1][c]=geotop::input::gDoubleNoValue;
            //	land->LC->co[land->LC->nrh][c]=geotop::input::gDoubleNoValue;
            land->LC[land->LC.getRows()-1][c]=geotop::input::gDoubleNoValue;
        }
    }
    //if(flag >= 0) write_map(geotop::common::Variables::files[flu], 1, par->format_out, land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

    if(par->state_pixel == 1){
        //	par->rc=new_longmatrix(par->chkpt->nrh,2);
        par->rc.resize(par->chkpt.getRows(),2+1);
        //	par->IDpoint=new_longvector(par->chkpt->nrh);
        par->IDpoint.resize(par->chkpt.getRows());
        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	par->rc->co[i][1]=row(par->chkpt->co[i][ptY], top->Z0->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            par->rc[i][1]=row(par->chkpt[i][ptY], top->Z0.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            //	par->rc->co[i][2]=col(par->chkpt->co[i][ptX], top->Z0->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            par->rc[i][2]=col(par->chkpt[i][ptX], top->Z0.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

            //	if (par->rc->co[i][1] == geotop::input::gDoubleNoValue || par->rc->co[i][2] == geotop::input::gDoubleNoValue) {
            if (par->rc[i][1] == geotop::input::gDoubleNoValue || par->rc[i][2] == geotop::input::gDoubleNoValue) {
                printf("Point #%4ld is out of the domain",i);

                fprintf(flog, "Point #%4ld is out of the domain",i);
                fclose(flog);

                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                fprintf(f, "Point #%4ld is out of the domain",i);
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report.");
            }

            //	if((long)land->LC->co[par->rc->co[i][1]][par->rc->co[i][2]]==geotop::input::gDoubleNoValue){
            if((long)land->LC[par->rc[i][1]][par->rc[i][2]]==geotop::input::gDoubleNoValue){
                printf("Point #%4ld corresponds to NOVALUE pixel",i);

                fprintf(flog, "Point #%4ld corresponds to NOVALUE pixel",i);
                fclose(flog);

                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                fprintf(f, "Point #%4ld corresponds to NOVALUE pixel",i);
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report.");
            }

            //	if ((long)par->chkpt->co[i][ptID]!=geotop::input::gDoubleNoValue) {
            if ((long)par->chkpt[i][ptID]!=geotop::input::gDoubleNoValue) {
                //	par->IDpoint->co[i]=(long)par->chkpt->co[i][ptID];
                par->IDpoint[i]=(long)par->chkpt[i][ptID];
            }else {
                //	par->IDpoint->co[i]=i;
                par->IDpoint[i]=i;
            }

        }
    }


    /****************************************************************************************************/

    //reading SKY VIEW FACTOR
    flag = file_exists(fsky, flog);
    if(flag == 1){
        //	top->sky=read_map(2,geotop::common::Variables::files[fsky], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fsky]), top->sky);
    }else{/*The sky view factor file "top->sky" must be calculated*/
        //	top->sky = new_doublematrix(top->Z0->nrh,top->Z0->nch);
        top->sky.resize(top->Z0.getRows(),top->Z0.getCols());
        if (par->sky == 0) {
            //	initialize_doublematrix(top->sky, 1.);
            top->sky.resize(top->Z0.getRows(),top->Z0.getCols(),1.);
        }else {
            //	curv = new_shortmatrix(top->Z0->nrh,top->Z0->nch);
            curv.resize(top->Z0.getRows(),top->Z0.getCols());
            nablaquadro_mask(top->Z0, curv, geotop::common::Variables::UV->U, geotop::common::Variables::UV->V);
            sky_view_factor(top->sky, 36, geotop::common::Variables::UV, top->Z0, curv, geotop::input::gDoubleNoValue);
            //	free_shortmatrix(curv);
        }
    }
    //if(flag >= 0) write_map(geotop::common::Variables::files[fsky], 0, par->format_out, top->sky, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

    /****************************************************************************************************/

    //reading DELAY
    flag = file_exists(fdelay, flog);
    if(flag == 1){
        //	land->delay = read_map(2,geotop::common::Variables::files[fdelay], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fdelay]), land->delay);
    }else{
        //	land->delay = new_doublematrix(top->Z0->nrh,top->Z0->nch);
        //	initialize_doublematrix(land->delay, 0);
        land->delay.resize(top->Z0.getRows(),top->Z0.getCols(),0);
    }
    //if(flag >= 0) write_map(geotop::common::Variables::files[fdelay], 0, par->format_out, land->delay, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

    /****************************************************************************************************/

    //reading SOIL MAP
    flag = file_exists(fsoil, flog);
    if(flag == 1){
        //	M=read_map(2,geotop::common::Variables::files[fsoil], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fsoil]), M);

        //	sl->type=copylong_doublematrix(M);
        copylong_doublematrix(sl->type, M);
        //	for(r=1;r<=land->LC->nrh;r++){
        for(r=1;r<land->LC.getRows();r++){
            //	for(c=1;c<=land->LC->nch;c++){
            for(c=1;c<land->LC.getCols();c++){
                //	if ((long)land->LC->co[r][c] != geotop::input::gDoubleNoValue) {
                if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue) {
                    //	if (sl->type->co[r][c] < 1 || sl->type->co[r][c] > par->nsoiltypes){
                    if (sl->type[r][c] < 1 || sl->type[r][c] > par->nsoiltypes){
                        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                        fprintf(f, "Error: It is not possible to assign Value < 1 or > nsoiltypes to the soil type map");
                        fclose(f);
                        t_error("Fatal Error! GEOtop is closed. See failing report (13).");
                    }
                }
            }
        }
    }else{
        //	M=copydoublematrix_const(par->soil_type_land_default, land->LC, geotop::input::gDoubleNoValue);
        copydoublematrix_const(par->soil_type_land_default, land->LC, M, geotop::input::gDoubleNoValue);

        //	sl->type=copylong_doublematrix(M);
        copylong_doublematrix(sl->type,M);
    }
    //if(flag >= 0) write_map(geotop::common::Variables::files[fsoil], 1, par->format_out, M, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
    //	free_doublematrix(M);


    /****************************************************************************************************/
    //SLOPE
    //	top->dzdE = new_doublematrix(land->LC->nrh, land->LC->nch);
    top->dzdE.resize(land->LC.getRows(), land->LC.getCols());
    //	top->dzdN = new_doublematrix(land->LC->nrh, land->LC->nch);
    top->dzdN.resize(land->LC.getRows(), land->LC.getCols());

    //	find_slope(UV->U->co[1], UV->U->co[2], top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue);
    find_slope(geotop::common::Variables::UV->U[1], geotop::common::Variables::UV->U[2], top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue);

    flag = file_exists(fslp, flog);
    if(flag == 1){
        //	top->slope=read_map(2,geotop::common::Variables::files[fslp], top->Z0, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);		//reads in degrees
        meteoio_readMap(string(geotop::common::Variables::files[fslp]), top->slope);
    }else{
        //	top->slope=find_max_slope(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue);
        find_max_slope(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue, top->slope);
    }
    //if(flag >= 0) write_map(geotop::common::Variables::files[fslp], 0, par->format_out, top->slope, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

    find_min_max(top->slope, geotop::input::gDoubleNoValue, &max, &min);

#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
    fprintf(flog,"Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#else
    printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
    fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#endif

    /****************************************************************************************************/
    //ASPECT
    flag = file_exists(fasp, flog);
    if(flag == 1){
        //	top->aspect=read_map(2,geotop::common::Variables::files[fasp], top->Z0, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fasp]), top->aspect);
    }else{
        //	top->aspect=find_aspect(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue);
        find_aspect(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue,top->aspect);
    }
	if(flag >= 0) write_map(geotop::common::Variables::files[fasp], 0, par->format_out, top->aspect, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);


    /****************************************************************************************************/
    //curvature
    //	top->curvature1=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    top->curvature1.resize(top->Z0.getRows(),top->Z0.getCols());
    //	top->curvature2=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    top->curvature2.resize(top->Z0.getRows(),top->Z0.getCols());
    //	top->curvature3=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    top->curvature3.resize(top->Z0.getRows(),top->Z0.getCols());
    //	top->curvature4=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    top->curvature4.resize(top->Z0.getRows(),top->Z0.getCols());

    //filtering
    //	M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    M.resize(top->Z0.getRows(),top->Z0.getCols());
    //	multipass_topofilter(par->lowpass_curvatures, top->Z0, M, geotop::input::gDoubleNoValue, 1);
    multipass_topofilter(par->lowpass_curvatures, top->Z0, M, geotop::input::gDoubleNoValue, 1);
    //	curvature(UV->U->co[1], UV->U->co[2], M, top->curvature1, top->curvature2, top->curvature3, top->curvature4, geotop::input::gDoubleNoValue);
    curvature(geotop::common::Variables::UV->U[1], geotop::common::Variables::UV->U[2], M, top->curvature1, top->curvature2, top->curvature3, top->curvature4, geotop::input::gDoubleNoValue);
    //	free_doublematrix(M);

    if(geotop::common::Variables::files[fcurv] != geotop::input::gStringNoValue){
        //	temp = join_strings(geotop::common::Variables::files[fcurv], "N-S");
        temp =geotop::common::Variables::files[fcurv] + string("N-S");
        write_map(temp, 0, par->format_out, top->curvature1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        //	free(temp);

        //	temp = join_strings(geotop::common::Variables::files[fcurv], "W-E");
        temp =geotop::common::Variables::files[fcurv] + string("W-E");
        write_map(temp, 0, par->format_out, top->curvature2, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        //	free(temp);

        //	temp = join_strings(geotop::common::Variables::files[fcurv], "NW-SE");
        temp =geotop::common::Variables::files[fcurv] + string("NW-SE");
        write_map(temp, 0, par->format_out, top->curvature3, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        //	free(temp);

        //	temp = join_strings(geotop::common::Variables::files[fcurv], "NE-SW");
        temp =geotop::common::Variables::files[fcurv] + string("NE-SW");
        write_map(temp, 0, par->format_out, top->curvature4, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        //	free(temp);
    }

    find_min_max(top->curvature1, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature N-S Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature N-S Min:%12g  Max:%12g \n",min,max);
#else
    printf("Curvature N-S Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature N-S Min:%f  Max:%f \n",min,max);
#endif

    find_min_max(top->curvature2, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature W-E Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature W-E Min:%12g  Max:%12g \n",min,max);
#else
    printf("Curvature W-E Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature W-E Min:%f  Max:%f \n",min,max);
#endif

    find_min_max(top->curvature3, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
#else
    printf("Curvature NW-SE Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature NW-SE Min:%f  Max:%f \n",min,max);
#endif

    find_min_max(top->curvature4, geotop::input::gDoubleNoValue, &max, &min);

#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
#else
    printf("Curvature NE-SW Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
#endif

    /****************************************************************************************************/
    //Channel network (in top->pixel_type)

    //pixel type = 0 land pixel (if it is on the border, the border is impermeable, water is free only on the surface)
    //pixel type = 1 or 2 land pixel (it it is on the border, the border is permeable above an user-defined elevation in the saturated part, weir-wise)
    //pixel type = 10 channel pixel (if it is on the border, the border is impermeable, water is free only on the surface)
    //pixel type = -1 land pixel where an incoming discharge from outside is considered (as rain)

    flag = file_exists(fnet, flog);
    if(flag == 1){
        //M=read_map(2,geotop::common::Variables::files[fnet], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fnet]), M);

        //	top->pixel_type=copyshort_doublematrix(M);
        copyshort_doublematrix(top->pixel_type, M);

        cont = 0;
        //	for(r=1;r<=top->Z0->nrh;r++){
        for(r=1;r<top->Z0.getRows();r++){
            //	for(c=1;c<=top->Z0->nch;c++){
            for(c=1;c<top->Z0.getCols();c++){
                //	if((long)land->LC->co[r][c]!=geotop::input::gDoubleNoValue){
                if((long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                    //	if(top->pixel_type->co[r][c]!=0 && top->pixel_type->co[r][c]!=1 && top->pixel_type->co[r][c]!=2 && top->pixel_type->co[r][c]!=10){
                    if(top->pixel_type[r][c]!=0 && top->pixel_type[r][c]!=1 && top->pixel_type[r][c]!=2 && top->pixel_type[r][c]!=10){
                        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                        fprintf(f, "Error: Only the following values are admitted in the network map: 0, 1, 10\n");
                        fclose(f);
                        t_error("Fatal Error! GEOtop is closed. See failing report (14).");
                    }
                    //	if(top->pixel_type->co[r][c]==10) cont++;
                    if(top->pixel_type[r][c]==10) cont++;
                }
            }
        }

        printf("Channel networks has %ld pixels set to channel\n",cont);

        if(flag >= 0) write_map(geotop::common::Variables::files[fnet], 1, par->format_out, M, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        //	free_doublematrix(M);

    }else{

        //	top->pixel_type=new_shortmatrix(land->LC->nrh, land->LC->nch);
        //	initialize_shortmatrix(top->pixel_type, 0);
        top->pixel_type.resize(land->LC.getRows(), land->LC.getCols(), 0);

    }

    /****************************************************************************************************/

    //border
    //	top->is_on_border=new_shortmatrix(land->LC->nrh, land->LC->nch);
    top->is_on_border.resize(land->LC.getRows(), land->LC.getCols());
    //	for(r=1;r<=land->LC->nrh;r++){
    for(r=1;r<land->LC.getRows();r++){
        //	for(c=1;c<=land->LC->nch;c++){
        for(c=1;c<land->LC.getCols();c++){
            //	if ( (long)land->LC->co[r][c]!=geotop::input::gDoubleNoValue){
            if ( (long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                //	top->is_on_border->co[r][c] = is_boundary(r, c, land->LC, geotop::input::gDoubleNoValue);
                top->is_on_border[r][c] = is_boundary(r, c, land->LC, geotop::input::gDoubleNoValue);
            }else{
                //	top->is_on_border->co[r][c] = -1;
                top->is_on_border[r][c] = -1;
            }
        }
    }

    //count the pixels having pixel_type = 1, 2 or -1
    cont = 0;
    //	for(r=1;r<=top->Z0->nrh;r++){
    for(r=1;r<top->Z0.getRows();r++){
        //	for(c=1;c<=top->Z0->nch;c++){
        for(c=1;c<top->Z0.getCols();c++){
            //	if(top->is_on_border->co[r][c]==1){
            if(top->is_on_border[r][c]==1){
                //	if (top->pixel_type->co[r][c] == -1 || top->pixel_type->co[r][c] == 1 || top->pixel_type->co[r][c] == 2) cont ++;
                if (top->pixel_type[r][c] == -1 || top->pixel_type[r][c] == 1 || top->pixel_type[r][c] == 2) cont ++;
            }
        }
    }

    //	top->BC_counter = new_longmatrix(top->Z0->nrh, top->Z0->nch);
    //	initialize_longmatrix(top->BC_counter, 0);
    top->BC_counter.resize(top->Z0.getRows(), top->Z0.getCols(), 0);


    if (cont > 0) {
        //	top->BC_DepthFreeSurface = new_doublevector(cont);
        top->BC_DepthFreeSurface.resize(cont+1);

        cont = 0;
        //	for(r=1;r<=top->Z0->nrh;r++){
        for(r=1;r<top->Z0.getRows();r++){
            //	for(c=1;c<=top->Z0->nch;c++){
            for(c=1;c<top->Z0.getCols();c++){
                //	if(top->is_on_border->co[r][c]==1){
                if(top->is_on_border[r][c]==1){
                    //	if (top->pixel_type->co[r][c] == -1 || top->pixel_type->co[r][c] == 1 || top->pixel_type->co[r][c] == 2){
                    if (top->pixel_type[r][c] == -1 || top->pixel_type[r][c] == 1 || top->pixel_type[r][c] == 2){
                        cont ++;
                        //	top->BC_counter->co[r][c] = cont;
                        top->BC_counter[r][c] = cont;
                        //	top->BC_DepthFreeSurface->co[cont] = par->DepthFreeSurface; //[mm]
                        top->BC_DepthFreeSurface[cont] = par->DepthFreeSurface; //[mm]
                    }
                }
            }
        }
    }else {
        //	top->BC_DepthFreeSurface = new_doublevector(1);
        //	initialize_doublevector(top->BC_DepthFreeSurface, geotop::input::gDoubleNoValue);
        top->BC_DepthFreeSurface.resize(2,geotop::input::gDoubleNoValue);
    }

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//void read_optionsfile_point(PAR *par, TOPO *top, LAND *land, SOIL *sl, TIMES *times, INIT_TOOLS *IT, FILE *flog){
void read_optionsfile_point(Par *par, Topo *top, Land *land, Soil *sl, Times *times, InitTools *IT, FILE *flog){

    long i, r, c, num_lines;
    //	DOUBLEMATRIX *Q, *P, *R, *S, *T, *Z, *LU;
    GeoMatrix<double> Q, P, R, S, T, Z, LU;
    //	SHORTMATRIX *curv;
    GeoMatrix<short> curv;
    short read_dem, read_lu, read_soil, read_sl, read_as, read_sk, read_curv, flag, coordinates;
    //	char *temp;
    string temp;
    double min, max;
    FILE *f;

    //4. CALCULATE TOPOGRAPHIC PROPERTIES
    //check if there are point coordinates
    coordinates = 1;
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i<par->chkpt.getRows();i++){
        //	if ( (long)par->chkpt->co[i][ptX]==geotop::input::gDoubleNoValue || (long)par->chkpt->co[i][ptY]==geotop::input::gDoubleNoValue ) coordinates = 0;
        if ( (long)par->chkpt[i][ptX]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptY]==geotop::input::gDoubleNoValue ) coordinates = 0;
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
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i< par->chkpt.getRows();i++){
        //		if((long)par->chkpt->co[i][ptLC]==geotop::input::gDoubleNoValue || (long)par->chkpt->co[i][ptSY]==geotop::input::gDoubleNoValue ||
        //		   (long)par->chkpt->co[i][ptS]==geotop::input::gDoubleNoValue || (long)par->chkpt->co[i][ptA]==geotop::input::gDoubleNoValue ||
        //		   (long)par->chkpt->co[i][ptCNS]==geotop::input::gDoubleNoValue || (long)par->chkpt->co[i][ptCWE]==geotop::input::gDoubleNoValue ||
        //		   (long)par->chkpt->co[i][ptCNwSe]==geotop::input::gDoubleNoValue || (long)par->chkpt->co[i][ptCNeSw]==geotop::input::gDoubleNoValue){
        if((long)par->chkpt[i][ptLC]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptSY]==geotop::input::gDoubleNoValue ||
                (long)par->chkpt[i][ptS]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptA]==geotop::input::gDoubleNoValue ||
                (long)par->chkpt[i][ptCNS]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptCWE]==geotop::input::gDoubleNoValue ||
                (long)par->chkpt[i][ptCNwSe]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptCNeSw]==geotop::input::gDoubleNoValue){
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
            printf("Warning: Dem file %s present\n",geotop::common::Variables::files[fdem+1].c_str());
            fprintf(flog,"Warning: Dem file %s present\n",geotop::common::Variables::files[fdem+1].c_str());

            //	Q=new_doublematrix(1,1);
            //	Z=read_map(0,geotop::common::Variables::files[fdem], Q, geotop::common::Variables::UV, geotop::input::gDoubleNoValue); //topography
            meteoio_readDEM(Z);
            //	free_doublematrix(Q);

            //	Q=new_doublematrix(Z->nrh,Z->nch);
            Q.resize(Z.getRows(),Z.getCols());
            multipass_topofilter(par->lowpass, Z, Q, geotop::input::gDoubleNoValue, 1);
            //	copy_doublematrix(Q, Z);
            Z = Q;
            //	free_doublematrix(Q);

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
        //	par->r_points=new_longvector(par->chkpt->nrh);
        par->r_points.resize(par->chkpt.getRows());
        //	par->c_points=new_longvector(par->chkpt->nrh);
        par->c_points.resize(par->chkpt.getRows());

        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	par->r_points->co[i]=row(par->chkpt->co[i][ptY], Z->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            par->r_points[i]=row(par->chkpt[i][ptY], Z.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            //	par->c_points->co[i]=col(par->chkpt->co[i][ptX], Z->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            par->c_points[i]=col(par->chkpt[i][ptX], Z.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            //	if((long)par->chkpt->co[i][ptZ]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptZ]=Z->co[par->r_points->co[i]][par->c_points->co[i]];
            if((long)par->chkpt[i][ptZ]==geotop::input::gDoubleNoValue) par->chkpt[i][ptZ]=Z[par->r_points[i]][par->c_points[i]];
        }
    }

    //b. read land use
    read_lu=0;
    //if(par->recover>0) read_lu=1;
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i<par->chkpt.getRows();i++){
        //	if((long)par->chkpt->co[i][ptLC]==geotop::input::gDoubleNoValue) read_lu=1;
        if((long)par->chkpt[i][ptLC]==geotop::input::gDoubleNoValue) read_lu=1;
    }
    if(read_lu==1 && coordinates==0) read_lu=0;
    if(read_lu==1){
        flag = file_exists(flu, flog);
        if(flag == 1){
            meteoio_readMap(string(geotop::common::Variables::files[flu]), LU); //HACK: add consitency check in meteoioplugin
            /*if(read_dem==0){
                Q=new_doublematrix(1,1);
                LU=read_map(0,geotop::common::Variables::files[flu], Q, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                free_doublematrix(Q);
            }else{
                LU=read_map(1,geotop::common::Variables::files[flu], Z, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                }*/

        }else{
            printf("Warning: Landuse file not present, uniform cover considered\n");
            fprintf(flog,"Warning: Landuse file not present, uniform cover considered\n");
            if(read_dem==1){
                //	LU=copydoublematrix_const(1.0, Z, geotop::input::gDoubleNoValue);
                copydoublematrix_const(1.0, Z, LU,geotop::input::gDoubleNoValue);
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
        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	if((long)par->chkpt->co[i][ptLC]==geotop::input::gDoubleNoValue){
            if((long)par->chkpt[i][ptLC]==geotop::input::gDoubleNoValue){
                //	r=row(par->chkpt->co[i][ptY], LU->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                r=row(par->chkpt[i][ptY], LU.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	c=col(par->chkpt->co[i][ptX], LU->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                c=col(par->chkpt[i][ptX], LU.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	par->chkpt->co[i][ptLC]=LU->co[r][c];
                par->chkpt[i][ptLC]=LU[r][c];
            }
        }
    }

    //c. read soil type
    read_soil=0;
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i<par->chkpt.getRows();i++){
        //	if((long)par->chkpt->co[i][ptSY]==geotop::input::gDoubleNoValue) read_soil=1;
        if((long)par->chkpt[i][ptSY]==geotop::input::gDoubleNoValue) read_soil=1;
    }
    if(read_soil==1 && coordinates==0) read_soil=0;
    if(read_soil==1){
        flag = file_exists(fsoil, flog);
        if(flag == 1){
            meteoio_readMap(string(geotop::common::Variables::files[fsoil]), P); //HACK: add consitency check in meteoioplugin
            /*
            if(read_dem==0){
                Q=new_doublematrix(1,1);
                P=read_map(0,geotop::common::Variables::files[fsoil], Q, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                free_doublematrix(Q);
            }else{
                P=read_map(1,geotop::common::Variables::files[fsoil], Z, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }
            */
        }else{
            printf("Warning: Soiltype file not present\n");
            fprintf(flog,"Warning: Soiltype file not present\n");
            read_soil=0;
        }
    }
    if(read_soil==1){
        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	if((long)par->chkpt->co[i][ptSY]==geotop::input::gDoubleNoValue){
            if((long)par->chkpt[i][ptSY]==geotop::input::gDoubleNoValue){
                //	r=row(par->chkpt->co[i][ptY], P->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                r=row(par->chkpt[i][ptY], P.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	c=col(par->chkpt->co[i][ptX], P->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                c=col(par->chkpt[i][ptX], P.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	par->chkpt->co[i][ptSY]=P->co[r][c];
                par->chkpt[i][ptSY]=P[r][c];
            }
        }
        //	free_doublematrix(P);
    }

    //d. read slope
    read_sl=0;
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i< par->chkpt.getRows();i++){
        //	if((long)par->chkpt->co[i][ptS]==geotop::input::gDoubleNoValue) read_sl=1;
        if((long)par->chkpt[i][ptS]==geotop::input::gDoubleNoValue) read_sl=1;
    }
    if(read_sl==1 && coordinates==0) read_sl=0;
    if(read_sl==1){
        flag = file_exists(fslp, flog);
        if(flag == 1){
            meteoio_readMap(string(geotop::common::Variables::files[fslp]), P); //HACK: add consitency check in meteoioplugin
            /*
            if(read_dem==0){
                Q=new_doublematrix(1,1);
                P=read_map(0,geotop::common::Variables::files[fslp], Q, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                free_doublematrix(Q);
            }else{
                P=read_map(1,geotop::common::Variables::files[fslp], Z, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }
            */
        }else{
            if(read_dem==0){
                printf("Warning: Slopes file not present\n");
                fprintf(flog,"Warning: Slopes file not present\n");
                read_sl=0;
            }else{
                //	Q=new_doublematrix(Z->nrh,Z->nch);
                Q.resize(Z.getRows(),Z.getCols());
                //	R=new_doublematrix(Z->nrh,Z->nch);
                R.resize(Z.getRows(),Z.getCols());
                //	find_slope(UV->U->co[1], UV->U->co[2], Z, Q, R, geotop::input::gDoubleNoValue);
                find_slope(geotop::common::Variables::UV->U[1], geotop::common::Variables::UV->U[2], Z, Q, R, geotop::input::gDoubleNoValue);
                //	P=find_max_slope(Z, Q, R, geotop::input::gDoubleNoValue);
                find_max_slope(Z, Q, R, geotop::input::gDoubleNoValue, P);
                //	free_doublematrix(Q);
                //	free_doublematrix(R);
                if(flag==0) write_map(geotop::common::Variables::files[fslp], 0, par->format_out, P, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }
        }
    }

    if(read_sl==1){
        find_min_max(P, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
        printf("Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
        fprintf(flog,"Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#else
        printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
        fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#endif

        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	if((long)par->chkpt->co[i][ptS]==geotop::input::gDoubleNoValue){
            if((long)par->chkpt[i][ptS]==geotop::input::gDoubleNoValue){
                //	r=row(par->chkpt->co[i][ptY], P->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                r=row(par->chkpt[i][ptY], P.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	c=col(par->chkpt->co[i][ptX], P->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                c=col(par->chkpt[i][ptX], P.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	par->chkpt->co[i][ptS]=P->co[r][c];
                par->chkpt[i][ptS]=P[r][c];
            }
        }
        //	free_doublematrix(P);
    }

    //e. read aspect
    read_as=0;
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i< par->chkpt.getRows();i++){
        //	if((long)par->chkpt->co[i][ptA]==geotop::input::gDoubleNoValue) read_as=1;
        if((long)par->chkpt[i][ptA]==geotop::input::gDoubleNoValue) read_as=1;
    }
    if(read_as==1 && coordinates==0) read_as=0;
    if(read_as==1){
        flag = file_exists(fasp, flog);
        if(flag == 1){
            meteoio_readMap(string(geotop::common::Variables::files[fasp]), P); //HACK: add consitency check in meteoioplugin
            /*
            if(read_dem==0){
                Q=new_doublematrix(1,1);
                P=read_map(0,geotop::common::Variables::files[fasp], Q, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                free_doublematrix(Q);
            }else{
                P=read_map(1,geotop::common::Variables::files[fasp], Z, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }
            */
        }else{
            if(read_dem==0){
                printf("Warning: Aspect file not present\n");
                fprintf(flog,"Warning: Aspect file not present\n");
                read_as=0;
            }else{
                //	Q=new_doublematrix(Z->nrh,Z->nch);
                Q.resize(Z.getRows(),Z.getCols());
                //	R=new_doublematrix(Z->nrh,Z->nch);
                R.resize(Z.getRows(),Z.getCols());
                //	find_slope(UV->U->co[1], UV->U->co[2], Z, Q, R, geotop::input::gDoubleNoValue);
                find_slope(geotop::common::Variables::UV->U[1], geotop::common::Variables::UV->U[2], Z, Q, R, geotop::input::gDoubleNoValue);
                //	P=find_aspect(Z, Q, R, geotop::input::gDoubleNoValue);
                find_aspect(Z, Q, R, geotop::input::gDoubleNoValue, P);

                //	free_doublematrix(Q);
                //	free_doublematrix(R);
                if(flag==0) write_map(geotop::common::Variables::files[fasp], 0, par->format_out, P, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }
        }
    }

    if(read_as==1){
        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	if((long)par->chkpt->co[i][ptA]==geotop::input::gDoubleNoValue){
            if((long)par->chkpt[i][ptA]==geotop::input::gDoubleNoValue){
                //	r=row(par->chkpt->co[i][ptY], P->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                r=row(par->chkpt[i][ptY], P.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	c=col(par->chkpt->co[i][ptX], P->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                c=col(par->chkpt[i][ptX], P.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	par->chkpt->co[i][ptA]=P->co[r][c];
                par->chkpt[i][ptA]=P[r][c];
            }
        }
        //	free_doublematrix(P);
    }

    //f. sky view factor file
    read_sk=0;
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i< par->chkpt.getRows();i++){
        //	if((long)par->chkpt->co[i][ptSKY]==geotop::input::gDoubleNoValue) read_sk=1;
        if((long)par->chkpt[i][ptSKY]==geotop::input::gDoubleNoValue) read_sk=1;
    }
    if(read_sk==1 && coordinates==0) read_sk=0;
    if(read_sk==1){
        flag = file_exists(fsky, flog);
        if(flag == 1){
            meteoio_readMap(string(geotop::common::Variables::files[fsky]), P); //HACK: add consitency check in meteoioplugin
            /*
            if(read_dem==0){
                Q=new_doublematrix(1,1);
                P=read_map(0,geotop::common::Variables::files[fsky], Q, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                free_doublematrix(Q);
            }else{
                P=read_map(1,geotop::common::Variables::files[fsky], Z, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }
            */
        }else{
            if(read_dem==0){
                printf("Warning: Sky view factor file not present\n");
                fprintf(flog,"Warning: Sky view factor file not present\n");
                read_sk=0;
            }else{
                //	P=new_doublematrix(Z->nrh,Z->nch);
                P.resize(Z.getRows(),Z.getCols());
                //	curv=new_shortmatrix(Z->nrh,Z->nch);
                curv.resize(Z.getRows(),Z.getCols());
                nablaquadro_mask(Z, curv, geotop::common::Variables::UV->U, geotop::common::Variables::UV->V);
                sky_view_factor(P, 36, geotop::common::Variables::UV, Z, curv, geotop::input::gDoubleNoValue);
                //	free_shortmatrix(curv);
                if(flag==0) write_map(geotop::common::Variables::files[fsky], 0, par->format_out, P, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }
        }
    }

    if(read_sk==1){
        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	if((long)par->chkpt->co[i][ptSKY]==geotop::input::gDoubleNoValue){
            if((long)par->chkpt[i][ptSKY]==geotop::input::gDoubleNoValue){
                //	r=row(par->chkpt->co[i][ptY], P->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                r=row(par->chkpt[i][ptY], P.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	c=col(par->chkpt->co[i][ptX], P->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                c=col(par->chkpt[i][ptX], P.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	par->chkpt->co[i][ptSKY]=P->co[r][c];
                par->chkpt[i][ptSKY]=P[r][c];
            }
        }
        //	free_doublematrix(P);
    }

    //g.curvature
    read_curv=0;
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i< par->chkpt.getRows();i++){
        //	if( (long)par->chkpt->co[i][ptCNS]==geotop::input::gDoubleNoValue || (long)par->chkpt->co[i][ptCWE]==geotop::input::gDoubleNoValue ||
        //	    (long)par->chkpt->co[i][ptCNwSe]==geotop::input::gDoubleNoValue || (long)par->chkpt->co[i][ptCNeSw]==geotop::input::gDoubleNoValue  ) read_curv=1;
        if( (long)par->chkpt[i][ptCNS]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptCWE]==geotop::input::gDoubleNoValue ||
                (long)par->chkpt[i][ptCNwSe]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptCNeSw]==geotop::input::gDoubleNoValue  ) read_curv=1;
    }
    if(read_curv==1 && coordinates==0) read_curv=0;
    if(read_curv==1){
        if(read_dem==0){
            printf("Warning: Dem file is not present, and therefore it is not possible to calculate curvature\n");
            fprintf(flog,"Warning: Dem file is not present, and therefore it is not possible to calculate curvature\n");
            read_curv=0;
        }else{
            //	Q=new_doublematrix(Z->nrh,Z->nch);
            Q.resize(Z.getRows(),Z.getCols());

            //	P=new_doublematrix(Z->nrh,Z->nch);
            P.resize(Z.getRows(),Z.getCols());

            //	R=new_doublematrix(Z->nrh,Z->nch);
            R.resize(Z.getRows(),Z.getCols());

            //	S=new_doublematrix(Z->nrh,Z->nch);
            S.resize(Z.getRows(),Z.getCols());
            //	T=new_doublematrix(Z->nrh,Z->nch);
            T.resize(Z.getRows(),Z.getCols());

            multipass_topofilter(par->lowpass_curvatures, Z, Q, geotop::input::gDoubleNoValue, 1);
            //	curvature(UV->U->co[1], UV->U->co[2], Q, P, R, S, T, geotop::input::gDoubleNoValue);
            curvature(geotop::common::Variables::UV->U[1], geotop::common::Variables::UV->U[2], Q, P, R, S, T, geotop::input::gDoubleNoValue);
            //	free_doublematrix(Q);

            if(geotop::common::Variables::files[fcurv] != geotop::input::gStringNoValue){
                //	temp = join_strings(geotop::common::Variables::files[fcurv], "N-S");
                temp =geotop::common::Variables::files[fcurv] + string("N-S");
                write_map(temp, 0, par->format_out, P, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	free(temp);
                //	temp = join_strings(geotop::common::Variables::files[fcurv], "W-E");
                temp =geotop::common::Variables::files[fcurv] + string("W-E");
                write_map(temp, 0, par->format_out, R, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	free(temp);
                //	temp = join_strings(geotop::common::Variables::files[fcurv], "NW-SE");
                temp =geotop::common::Variables::files[fcurv] + string("NW-SE");
                write_map(temp, 0, par->format_out, S, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	free(temp);
                //	temp = join_strings(geotop::common::Variables::files[fcurv], "NE-SW");
                temp =geotop::common::Variables::files[fcurv] + string("NE-SW");
                write_map(temp, 0, par->format_out, T, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
                //	free(temp);
            }

            find_min_max(P, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature N-S Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature N-S Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature N-S Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature N-S Min:%f  Max:%f \n",min,max);
#endif

            find_min_max(R, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature W-E Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature W-E Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature W-E Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature W-E Min:%f  Max:%f \n",min,max);
#endif

            find_min_max(S, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature NW-SE Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature NW-SE Min:%f  Max:%f \n",min,max);
#endif

            find_min_max(T, geotop::input::gDoubleNoValue, &max, &min);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature NE-SW Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
#endif
        }
    }
    if(read_curv==1){
        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	r=row(par->chkpt->co[i][ptY], P->nrh, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            r=row(par->chkpt[i][ptY], P.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            //	c=col(par->chkpt->co[i][ptX], P->nch, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            c=col(par->chkpt[i][ptX], P.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            //	if((long)par->chkpt->co[i][ptCNS]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCNS]=P->co[r][c];
            if((long)par->chkpt[i][ptCNS]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCNS]=P[r][c];
            //	if((long)par->chkpt->co[i][ptCWE]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCWE]=R->co[r][c];
            if((long)par->chkpt[i][ptCWE]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCWE]=R[r][c];
            //	if((long)par->chkpt->co[i][ptCNwSe]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCNwSe]=S->co[r][c];
            if((long)par->chkpt[i][ptCNwSe]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCNwSe]=S[r][c];
            //	if((long)par->chkpt->co[i][ptCNeSw]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCNeSw]=T->co[r][c];
            if((long)par->chkpt[i][ptCNeSw]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCNeSw]=T[r][c];
        }
        //	free_doublematrix(P);
        //	free_doublematrix(R);
        //	free_doublematrix(S);
        //	free_doublematrix(T);
    }

    //h. no value check
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i< par->chkpt.getRows();i++){
        //	if((long)par->chkpt->co[i][ptZ]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptZ]=0.0;
        if((long)par->chkpt[i][ptZ]==geotop::input::gDoubleNoValue) par->chkpt[i][ptZ]=0.0;
        //	if((long)par->chkpt->co[i][ptLC]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptLC]=1.0;
        if((long)par->chkpt[i][ptLC]==geotop::input::gDoubleNoValue) par->chkpt[i][ptLC]=1.0;
        //	if((long)par->chkpt->co[i][ptSY]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptSY]=1.0;
        if((long)par->chkpt[i][ptSY]==geotop::input::gDoubleNoValue) par->chkpt[i][ptSY]=1.0;
        //	if((long)par->chkpt->co[i][ptS]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptS]=0.0;
        if((long)par->chkpt[i][ptS]==geotop::input::gDoubleNoValue) par->chkpt[i][ptS]=0.0;
        //	if((long)par->chkpt->co[i][ptA]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptA]=0.0;
        if((long)par->chkpt[i][ptA]==geotop::input::gDoubleNoValue) par->chkpt[i][ptA]=0.0;
        //	if((long)par->chkpt->co[i][ptSKY]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptSKY]=1.0;
        if((long)par->chkpt[i][ptSKY]==geotop::input::gDoubleNoValue) par->chkpt[i][ptSKY]=1.0;
        //	if((long)par->chkpt->co[i][ptCNS]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCNS]=0.0;
        if((long)par->chkpt[i][ptCNS]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCNS]=0.0;
        //	if((long)par->chkpt->co[i][ptCWE]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCWE]=0.0;
        if((long)par->chkpt[i][ptCWE]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCWE]=0.0;
        //	if((long)par->chkpt->co[i][ptCNwSe]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCNwSe]=0.0;
        if((long)par->chkpt[i][ptCNwSe]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCNwSe]=0.0;
        //	if((long)par->chkpt->co[i][ptCNeSw]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptCNeSw]=0.0;
        if((long)par->chkpt[i][ptCNeSw]==geotop::input::gDoubleNoValue) par->chkpt[i][ptCNeSw]=0.0;
        //	if((long)par->chkpt->co[i][ptDrDEPTH]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptDrDEPTH]=par->DepthFreeSurface;//[mm]
        if((long)par->chkpt[i][ptDrDEPTH]==geotop::input::gDoubleNoValue) par->chkpt[i][ptDrDEPTH]=par->DepthFreeSurface;//[mm]
        //	if((long)par->chkpt->co[i][ptMAXSWE]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptMAXSWE]=1.E10;//[mm]
        if((long)par->chkpt[i][ptMAXSWE]==geotop::input::gDoubleNoValue) par->chkpt[i][ptMAXSWE]=1.E10;//[mm]
        //	if((long)par->chkpt->co[i][ptLAT]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptLAT]=par->latitude;
        if((long)par->chkpt[i][ptLAT]==geotop::input::gDoubleNoValue) par->chkpt[i][ptLAT]=par->latitude;
        //	if((long)par->chkpt->co[i][ptLON]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptLON]=par->longitude;
        if((long)par->chkpt[i][ptLON]==geotop::input::gDoubleNoValue) par->chkpt[i][ptLON]=par->longitude;
        //	if((long)par->chkpt->co[i][ptID]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptID]=(double)i;
        if((long)par->chkpt[i][ptID]==geotop::input::gDoubleNoValue) par->chkpt[i][ptID]=(double)i;
        //	if((long)par->chkpt->co[i][ptHOR]==geotop::input::gDoubleNoValue) par->chkpt->co[i][ptHOR]=par->chkpt->co[i][ptID];
        if((long)par->chkpt[i][ptHOR]==geotop::input::gDoubleNoValue) par->chkpt[i][ptHOR]=par->chkpt[i][ptID];
    }

    //i.show results
    fprintf(flog,"\nPOINTS:\n");
    fprintf(flog,"ID,East[m],North[m],Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-],CurvatureN-S[1/m],CurvatureW-E[1/m],CurvatureNW-SE[1/m],CurvatureNE-SW[1/m],DepthFreeSurface[mm],Hor,maxSWE[mm],Lat[deg],Long[deg]\n");
    //	for(r=1;r<=par->chkpt->nrh;r++){
    for(r=1;r< par->chkpt.getRows();r++){
        for(c=1;c<ptTOT;c++){
            //	fprintf(flog,"%f",par->chkpt->co[r][c]);
            fprintf(flog,"%f",par->chkpt[r][c]);
            if (c<ptTOT) fprintf(flog, ",");
        }
        fprintf(flog,"\n");
    }

    //l. set UV
    if(read_dem==0 && read_lu==0 && read_soil==0 && read_sl==0 && read_as==0 && read_sk==0){
        geotop::common::Variables::UV->U.resize(4+1);
        geotop::common::Variables::UV->V.resize(2+1);
    }
    geotop::common::Variables::UV->U[2]=1.0;
    geotop::common::Variables::UV->U[1]=1.0;
    geotop::common::Variables::UV->U[4]=0.0;
    geotop::common::Variables::UV->U[3]=0.0;
    geotop::common::Variables::UV->V[2]=geotop::input::gDoubleNoValue;
    if(geotop::common::Variables::UV->V[2]<0){
        geotop::common::Variables::UV->V[1] = -1.;
    }else{
        geotop::common::Variables::UV->V[1] = 1.;
    }

    //m. set IT->LU
    /*if(par->recover>0){
        IT->LU=new_doublematrix(Z->nrh, Z->nch);
        copy_doublematrix(LU, IT->LU);
    }*/

    //n. deallocation
    //	if(read_dem==1) free_doublematrix(Z);
    //	if(read_lu==1) free_doublematrix(LU);
    //	if(par->recover==0 && read_dem==1){
    if(read_dem==1){
        //	free_longvector(par->r_points);
        //	free_longvector(par->c_points);
    }

    //5. SET CHECKPOINT
    if(par->state_pixel == 1){
        //	par->rc=new_longmatrix(par->chkpt->nrh,2);
        par->rc.resize(par->chkpt.getRows(),2+1);
        //	for(i=1;i<=par->chkpt->nrh;i++){
        for(i=1;i< par->chkpt.getRows();i++){
            //	par->rc->co[i][1]=1;
            par->rc[i][1]=1;
            //	par->rc->co[i][2]=i;
            par->rc[i][2]=i;
        }
    }

    //6. SET PROPERTIES
    //	top->East=new_doublematrix(1,par->chkpt->nrh);
    top->East.resize(1+1,par->chkpt.getRows());
    //	top->North=new_doublematrix(1,par->chkpt->nrh);
    top->North.resize(1+1,par->chkpt.getRows());
    //	top->Z0=new_doublematrix(1,par->chkpt->nrh);
    top->Z0.resize(1+1,par->chkpt.getRows());
    //	land->LC=new_doublematrix(1,par->chkpt->nrh);
    land->LC.resize(1+1,par->chkpt.getRows());
    //	land->delay=new_doublematrix(1,par->chkpt->nrh);
    land->delay.resize(1+1,par->chkpt.getRows());
    //	sl->type=new_longmatrix(1,par->chkpt->nrh);
    sl->type.resize(1+1,par->chkpt.getRows());
    //	top->slope=new_doublematrix(1,par->chkpt->nrh);
    top->slope.resize(1+1,par->chkpt.getRows());
    //	top->aspect=new_doublematrix(1,par->chkpt->nrh);
    top->aspect.resize(1+1,par->chkpt.getRows());
    //	top->curvature1=new_doublematrix(1,par->chkpt->nrh);
    top->curvature1.resize(1+1,par->chkpt.getRows());
    //	top->curvature2=new_doublematrix(1,par->chkpt->nrh);
    top->curvature2.resize(1+1,par->chkpt.getRows());
    //	top->curvature3=new_doublematrix(1,par->chkpt->nrh);
    top->curvature3.resize(1+1,par->chkpt.getRows());
    //	top->curvature4=new_doublematrix(1,par->chkpt->nrh);
    top->curvature4.resize(1+1,par->chkpt.getRows());
    //	top->sky=new_doublematrix(1,par->chkpt->nrh);
    top->sky.resize(1+1,par->chkpt.getRows());
    //	top->pixel_type=new_shortmatrix(1,par->chkpt->nrh);
    top->pixel_type.resize(1+1,par->chkpt.getRows());
    //	top->BC_counter=new_longmatrix(1,par->chkpt->nrh);
    top->BC_counter.resize(1+1,par->chkpt.getRows());
    //	top->BC_DepthFreeSurface=new_doublevector(par->chkpt->nrh);
    top->BC_DepthFreeSurface.resize(par->chkpt.getRows());
    //	par->maxSWE=new_doublematrix(1,par->chkpt->nrh);
    par->maxSWE.resize(2,par->chkpt.getRows());
    //	top->horizon_point=new_longmatrix(1,par->chkpt->nrh);
    top->horizon_point.resize(1+1,par->chkpt.getRows());
    //	top->dzdE=new_doublematrix(1,par->chkpt->nrh);
    top->dzdE.resize(1+1,par->chkpt.getRows());
    //	top->dzdN=new_doublematrix(1,par->chkpt->nrh);
    top->dzdN.resize(1+1,par->chkpt.getRows());
    //	top->latitude=new_doublematrix(1,par->chkpt->nrh);
    top->latitude.resize(2,par->chkpt.getRows());
    //	top->longitude=new_doublematrix(1,par->chkpt->nrh);
    top->longitude.resize(2,par->chkpt.getRows());
    //	par->IDpoint=new_longvector(par->chkpt->nrh);
    par->IDpoint.resize(par->chkpt.getRows());
    //	for(i=1;i<=par->chkpt->nrh;i++){
    for(i=1;i< par->chkpt.getRows();i++){
        //	top->East->co[1][i]=par->chkpt->co[i][ptX];
        top->East[1][i]=par->chkpt[i][ptX];
        //	top->North->co[1][i]=par->chkpt->co[i][ptY];
        top->North[1][i]=par->chkpt[i][ptY];
        //	top->Z0->co[1][i]=par->chkpt->co[i][ptZ];
        top->Z0[1][i]=par->chkpt[i][ptZ];
        //	land->LC->co[1][i]=par->chkpt->co[i][ptLC];
        land->LC[1][i]=par->chkpt[i][ptLC];

        //	if((long)land->LC->co[1][i] <= 0){
        if((long)land->LC[1][i] <= 0){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error:: Point %ld has land cover type <= 0. This is not admitted.\n",i);
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (15).");
        }

        //	sl->type->co[1][i]=(long)par->chkpt->co[i][ptSY];
        sl->type[1][i]=(long)par->chkpt[i][ptSY];

        //	if(sl->type->co[1][i] <= 0){
        if(sl->type[1][i] <= 0){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error:: Point %ld has soil type <= 0. This is not admitted.\n",i);
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (16).");
        }

        //	top->slope->co[1][i]=par->chkpt->co[i][ptS];
        top->slope[1][i]=par->chkpt[i][ptS];
        //	top->aspect->co[1][i]=par->chkpt->co[i][ptA];
        top->aspect[1][i]=par->chkpt[i][ptA];
        //	top->sky->co[1][i]=par->chkpt->co[i][ptSKY];
        top->sky[1][i]=par->chkpt[i][ptSKY];
        //	top->curvature1->co[1][i]=par->chkpt->co[i][ptCNS];
        top->curvature1[1][i]=par->chkpt[i][ptCNS];
        //	top->curvature2->co[1][i]=par->chkpt->co[i][ptCWE];
        top->curvature2[1][i]=par->chkpt[i][ptCWE];
        //	top->curvature3->co[1][i]=par->chkpt->co[i][ptCNwSe];
        top->curvature3[1][i]=par->chkpt[i][ptCNwSe];
        //	top->curvature4->co[1][i]=par->chkpt->co[i][ptCNeSw];
        top->curvature4[1][i]=par->chkpt[i][ptCNeSw];
        //	top->pixel_type->co[1][i]=1;
        top->pixel_type[1][i]=1;
        //	top->BC_counter->co[1][i]=i;
        top->BC_counter[1][i]=i;
        //	top->BC_DepthFreeSurface->co[i]=par->chkpt->co[i][ptDrDEPTH];
        top->BC_DepthFreeSurface[i]=par->chkpt[i][ptDrDEPTH];
        //	top->horizon_point->co[1][i]=(long)par->chkpt->co[i][ptHOR];
        top->horizon_point[1][i]=(long)par->chkpt[i][ptHOR];
        //	top->dzdE->co[1][i]=0.;
        top->dzdE[1][i]=0.;
        //	top->dzdN->co[1][i]=0.;
        top->dzdN[1][i]=0.;
        //	land->delay->co[1][i]=0.;
        land->delay[1][i]=0.;

        //	if(sl->type->co[1][i] <= 0){
        if(sl->type[1][i] <= 0){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error:: Point %ld has horizon type <= 0. This is not admitted.\n",i);
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (17).");
        }

        //	par->maxSWE->co[1][i]=par->chkpt->co[i][ptMAXSWE];
        par->maxSWE[1][i]=par->chkpt[i][ptMAXSWE];
        //	top->latitude->co[1][i]=par->chkpt->co[i][ptLAT];
        top->latitude[1][i]=par->chkpt[i][ptLAT];
        //	top->longitude->co[1][i]=par->chkpt->co[i][ptLON];
        top->longitude[1][i]=par->chkpt[i][ptLON];
        //	par->IDpoint->co[i]=(long)par->chkpt->co[i][ptID];
        par->IDpoint[i]=(long)par->chkpt[i][ptID];
    }

    //7. SET PAR
    //	for (i=1; i<=par->init_date->nh; i++) {
    for (size_t i=1; i< par->init_date.size(); i++) {
        //	par->output_soil->co[i]=0.;
        par->output_soil[i]=0.;
        //	par->output_snow->co[i]=0.;
        par->output_snow[i]=0.;
        //	par->output_glac->co[i]=0.;
        par->output_glac[i]=0.;
        //	par->output_surfenergy->co[i]=0.;
        par->output_surfenergy[i]=0.;
        //	par->output_vegetation->co[i]=0.;
        par->output_vegetation[i]=0.;
        //	par->output_meteo->co[i]=0.;
        par->output_meteo[i]=0.;
    }

    par->output_soil_bin = 0;
    par->output_snow_bin = 0;
    par->output_glac_bin = 0;
    par->output_surfenergy_bin = 0;
    par->output_meteo_bin = 0;

    //8. READ HORIZONS
    //find max top->horizon_point
    top->num_horizon_point=0;
    //	for(r=1;r<=top->horizon_point->nrh;r++){
    for(r=1;r<=top->horizon_point.getRows()-1;r++){
        //	for(c=1;c<=top->horizon_point->nch;c++){
        for(c=1;c<=top->horizon_point.getCols()-1;c++){
            //	if (top->horizon_point->co[r][c] > top->num_horizon_point) top->num_horizon_point = top->horizon_point->co[r][c];
            if (top->horizon_point[r][c] > top->num_horizon_point) top->num_horizon_point = top->horizon_point[r][c];
        }
    }
    top->horizon_height=(double ***)malloc(top->num_horizon_point*sizeof(double**));
    top->horizon_numlines=(long *)malloc(top->num_horizon_point*sizeof(long));
    for(i=1; i<=top->num_horizon_point; i++){

        c=0;
        do{
            flag = 0;
            //	if (c < par->chkpt->nrh) {
            if (c < par->chkpt.getRows()-1) {
                //	if (top->horizon_point->co[1][c+1] != i) c++;
                if (top->horizon_point[1][c+1] != i) c++;
            }
            //	if (c < par->chkpt->nrh) {
            if (c < par->chkpt.getRows()-1) {
                //	if (top->horizon_point->co[1][c+1] != i) flag=1;
                if (top->horizon_point[1][c+1] != i) flag=1;
            }

            //	}while (flag == 1 && c < par->chkpt->nrh);
        }while (flag == 1 && c < par->chkpt.getRows()-1);

        //	if (c < par->chkpt->nrh) {
        if (c < par->chkpt.getRows()-1) {
            top->horizon_height[i-1] = read_horizon(0, i,geotop::common::Variables::files[fhorpoint], IT->horizon_col_names, &num_lines, flog);
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

//void set_bedrock(SOIL *sl, CHANNEL *cnet, PAR *par, TOPO *top, DOUBLEMATRIX *LC, FILE *flog){
void set_bedrock(Soil *sl, Channel *cnet, Par *par, Topo *top, GeoMatrix<double>& LC, FILE *flog){

    //	DOUBLEMATRIX *B;
    GeoMatrix<double> B;
    //	DOUBLETENSOR *T;
    //	DOUBLEVECTOR *WT;
    GeoVector<double> WT;

    long i, j, l, r, c, sy, synew;
    double zlim, z;
    FILE *f;

    if (!mio::IOUtils::fileExists(string(geotop::common::Variables::files[fbed]) + string(ascii_esri))) {
        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
        fprintf(f, "Error:: File %s is missing. Please check if you have a bedrock topography map. If it is not available, remove the file name and keyword from input file\n",geotop::common::Variables::files[fbed+1].c_str());
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report (18).");
    }

    printf("A bedrock depth map has been assigned and read from %s\n\n",geotop::common::Variables::files[fbed].c_str());
    fprintf(flog,"A bedrock depth map has been assigned and read from %s\n\n",geotop::common::Variables::files[fbed].c_str());

    par->bedrock = 1;
    //	B = read_map(2,geotop::common::Variables::files[fbed], LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
    meteoio_readMap(string(geotop::common::Variables::files[fbed]), B);

    //	if (sl->init_water_table_depth->nh != sl->pa->ndh){
    if (sl->init_water_table_depth.size() != sl->pa.getDh()){
        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
        fprintf(f, "Error:: Error in bedrock calculations");
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report (19).");
    }

    //	rewrite soil type
    //	T=new_doubletensor(sl->pa->ndh, nsoilprop, Nl);
    GeoTensor<double> T;
    T.resize(sl->pa.getDh(), nsoilprop+1, geotop::common::Variables::Nl+1);
    //	for(i=1;i<=sl->pa->ndh;i++){
    for(i=1;i<sl->pa.getDh();i++){
        for(j=1;j<=nsoilprop;j++){
            for(l=1;l<=geotop::common::Variables::Nl;l++){
                //	T->co[i][j][l]=sl->pa->co[i][j][l];
                T(i,j,l) = sl->pa(i,j,l);
            }
        }
    }
    //	free_doubletensor(sl->pa);
    //	sl->pa=new_doubletensor(par->total_pixel+par->total_channel, nsoilprop, Nl);
    sl->pa.resize(par->total_pixel+par->total_channel+1, nsoilprop+1, geotop::common::Variables::Nl+1);

    //	rewrite initial water table depth
    //	WT=new_doublevector(sl->init_water_table_depth->nh);
    //	for(i=1;i<=sl->init_water_table_depth->nh;i++) {
    for(i=1;i<sl->init_water_table_depth.size();i++) {
        //	WT.push_back(sl->init_water_table_depth->co[i]);
        WT.data.push_back(sl->init_water_table_depth[i]);
        //	WT->co[i]=sl->init_water_table_depth->co[i];
    }
    //	free_doublevector(sl->init_water_table_depth);
    //	sl->init_water_table_depth=new_doublevector(par->total_pixel+par->total_channel);
    sl->init_water_table_depth.resize(par->total_pixel+par->total_channel+1);


    //	assign jdz (is needed later)
    //	for(i=1;i<=sl->pa->ndh;i++){
    for(i=1;i<sl->pa.getDh();i++){
        for(l=1;l<=geotop::common::Variables::Nl;l++){
            //	sl->pa->co[i][jdz][l]=T->co[1][jdz][l];
            sl->pa(i,jdz,l) = T(1,jdz,l);
        }
    }

    for (i=1; i<=par->total_pixel+par->total_channel; i++) {

        if (i<=par->total_pixel) {
            //	r = top->rc_cont->co[i][1];
            r = top->rc_cont[i][1];
            //	c = top->rc_cont->co[i][2];
            c = top->rc_cont[i][2];
            //	sy = sl->type->co[r][c];
            sy = sl->type[r][c];
            synew = i;
            //	sl->type->co[r][c] = synew;
            sl->type[r][c] = synew;
            z = 0.0;
        }else {
            //	r = cnet->r->co[i-par->total_pixel];
            r = cnet->r[i-par->total_pixel];
            //	c = cnet->c->co[i-par->total_pixel];
            c = cnet->c[i-par->total_pixel];
            //	sy = cnet->soil_type->co[i-par->total_pixel];
            sy = cnet->soil_type[i-par->total_pixel];
            synew = i;
            //	cnet->soil_type->co[i-par->total_pixel] = synew;
            cnet->soil_type[i-par->total_pixel] = synew;
            z = par->depr_channel;
        }

        //	sl->init_water_table_depth->co[synew] = WT[sy-1]; //WT->co[sy];
        sl->init_water_table_depth[synew] = WT[sy-1]; //WT->co[sy];

        //	zlim = B->co[r][c];
        zlim = B[r][c];

        for(l=1;l<=geotop::common::Variables::Nl;l++){

            //	z += 0.5*sl->pa->co[synew][jdz][l];
            z += 0.5*sl->pa(synew,jdz,l);

            if(z <= zlim){

                for(j=1;j<=nsoilprop;j++){
                    //	sl->pa->co[synew][j][l] = T->co[sy][j][l];
                    sl->pa(synew,j,l) = T(sy,j,l);
                }

            }else{

                for(j=1;j<=nsoilprop;j++){
                    //	sl->pa->co[synew][j][l] = sl->pa_bed(sy,j,l);
                    sl->pa(synew,j,l) = sl->pa_bed(sy,j,l);
                }

            }

            //	z += 0.5*sl->pa->co[synew][jdz][l];
            z += 0.5*sl->pa(synew,jdz,l);

        }
    }

    //	free_doubletensor(T);
    //	free_doublematrix(B);
    //	free_doublevector(WT);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//DOUBLETENSOR *find_Z_of_any_layer(DOUBLEMATRIX *Zsurface, DOUBLEMATRIX *slope, DOUBLEMATRIX *LC, SOIL *sl, short point){
GeoTensor<double> find_Z_of_any_layer(GeoMatrix<double>& Zsurface, GeoMatrix<double>& slope, GeoMatrix<double>& LC, Soil *sl, short point){

    //	DOUBLETENSOR *Z;
    GeoTensor<double> Z;
    double Zaverage, z, cosine;
    long l, r, c, n, sy;

    if(point!=1){
        Zaverage=0.;
        n=0;
        //	for(r=1;r<=Zsurface->nrh;r++){
        for(r=1;r<Zsurface.getRows();r++){
            //	for(c=1;c<=Zsurface->nch;c++){
            for(c=1;c<Zsurface.getCols();c++){
                //	if((long)LC->co[r][c]!=geotop::input::gDoubleNoValue){
                if((long)LC[r][c]!=geotop::input::gDoubleNoValue){
                    n++;
                    //	Zaverage += Zsurface->co[r][c];
                    Zaverage += Zsurface[r][c];
                }
            }
        }
        Zaverage/=(double)n;
    }


    //	Z=new_doubletensor0(sl->pa->nch, Zsurface->nrh, Zsurface->nch);
    //	initialize_doubletensor(Z, geotop::input::gDoubleNoValue);
    Z.resize(sl->pa.getCh(), Zsurface.getRows(), Zsurface.getCols(), geotop::input::gDoubleNoValue);


    //	for(r=1;r<=Zsurface->nrh;r++){
    for(r=1;r<Zsurface.getRows();r++){
        //	for(c=1;c<=Zsurface->nch;c++){
        for(c=1;c<Zsurface.getCols();c++){
            //	if((long)LC->co[r][c]!=geotop::input::gDoubleNoValue){
            if((long)LC[r][c]!=geotop::input::gDoubleNoValue){

                //	cosine = cos(slope->co[r][c]*GTConst::Pi/180.);
                cosine = cos(slope[r][c]*GTConst::Pi/180.);

                //	sy=sl->type->co[r][c];
                sy=sl->type[r][c];
                if (point!=1){
                    //	z=1.E3*(Zsurface->co[r][c]-Zaverage);//[mm]
                    z=1.E3*(Zsurface[r][c]-Zaverage);//[mm]
                }else {
                    z=0.;
                }

                l=0;
                //	Z->co[l][r][c]=z;
                Z[l][r][c]=z;

                do{
                    l++;
                    //z -= 0.5*sl->pa->co[sy][jdz][l]*cosine;
                    z -= 0.5 * sl->pa(sy,jdz,l) * cosine;
                    //	Z->co[l][r][c]=z;
                    Z[l][r][c]=z;
                    //z -= 0.5*sl->pa->co[sy][jdz][l]*cosine;
                    z -= 0.5 * sl->pa(sy,jdz,l) * cosine;
                    //}while(l<sl->pa->nch);
                } while(l < sl->pa.getCh()-1);
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
/*
    printf("Attempting to read '%s' in the file '%s': ",keywords_char[key+nmet+nsoilprop+2].c_str(),geotop::common::Variables::files[key].c_str());
    fprintf(flog,"Attempting to read '%s' in the file '%s': ",keywords_char[key+nmet+nsoilprop+2].c_str(),geotop::common::Variables::files[key].c_str());
*/

    if(geotop::common::Variables::files[key] == geotop::input::gStringNoValue){
        printf("not present in file list\n");
        fprintf(flog,"not present in file list\n");
        return (-1);
    }else{
        bool is_present = mio::IOUtils::fileExists(string(geotop::common::Variables::files[key]) + string(ascii_esri));

        if (is_present) {
            printf("EXISTING in format %d\n", 3);
            fprintf(flog, "EXISTING in format %d\n", 3);
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

void initialize_soil_state(SoilState *S, long n, long nl){

    //	S->T = new_doublematrix(nl, n);
    //	initialize_doublematrix(S->T, 0.);
    S->T.resize(nl+1, n+1, 0.);

    //	S->P = new_doublematrix0_(nl, n);
    //	initialize_doublematrix(S->P, 0.);
    //	TODO: notice that "new_doublematrix0_()"
    S->P.resize(nl+1, n+1, 0.);

    //	S->thi = new_doublematrix(nl, n);
    //	initialize_doublematrix(S->thi, 0.);
    S->thi.resize(nl+1, n+1,0.0);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void copy_soil_state(SoilState *from, SoilState *to){

    long l,i;
    long nl=from->T.getRows(), n=from->T.getCols();

    for (i=1; i<n; i++) {
        //	to->P->co[0][i] = from->P->co[0][i];
        to->P[0][i] = from->P[0][i];
        for (l=1; l<nl; l++) {
            to->P[l][i] = from->P[l][i];
            to->T[l][i] = from->T[l][i];
            to->thi[l][i] = from->thi[l][i];
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initialize_veg_state(StateVeg *V, long n){
    //	V->Tv = new_doublevector(n);
    //	initialize_doublevector(V->Tv, 0.);
    V->Tv.resize(n+1,0.0);

    //	V->wsnow = new_doublevector(n);
    //	initialize_doublevector(V->wsnow, 0.);
    V->wsnow.resize(n+1,0.0);

    //	V->wrain = new_doublevector(n);
    //	initialize_doublevector(V->wrain, 0.);
    V->wrain.resize(n+1,0.);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void copy_veg_state(StateVeg *from, StateVeg *to){

    long i, n=from->Tv.size();
    for (i=1; i<n; i++) {
        to->Tv[i] = from->Tv[i];
        to->wrain[i] = from->wrain[i];
        to->wsnow[i] = from->wsnow[i];
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
