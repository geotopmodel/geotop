
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 2.0.0
 
 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
  
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
 
#include "config.h"
#include "input.h"
#include "parameters.h"
#include <unistd.h>
#include <inputKeywords.h>
#include "geotop_common.h"

#ifdef WITH_LOGGER
#include <iostream>
#include "global_logger.h"
#endif

using namespace std ;

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//! Subroutine which reads input data, performs  geomporphological analisys and allocates data

void get_all_input(long argc, char *argv[], Topo *top, Soil *sl, Land *land, Meteo *met, Water *wat, Channel *cnet,
                   Par *par, Energy *egy, Snow *snow, Glacier *glac, Times *times, mio::IOManager& iomanager)

{

    FILE *flog, *f;
    GeoMatrix<double> M;
    InitTools *IT;

    size_t a;
    short success, added_JDfrom0=0, added_wind_xy=0, added_wind_dir=0, added_cloud=0, added_Tdew=0, added_RH=0, added_Pint=0;
    long l, r, c, i, ist, j, n, sy, num_cols, num_lines, day, month, year, hour, minute;
    double z, th_oversat, JD, k_snowred, maxSWE, SWE, D, cosslope, **matrix;
	short count_file_missing=0;  // needed to check if some output map files are present or not 
    std::vector<std::string> temp2 ;
    string temp;


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

        WORKING_DIRECTORY=read_option_string(argc,argv,"-wpath",".",0); // assign_string(argv[1]); // MODIFY HERE EC
#endif	
    }
#ifdef WITH_LOGGER
    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();

    lg->writeAll("STATEMENT:\n");
    lg->writeAll("\n");
    lg->writeAll("GEOtop 2.0.0 'MOAB' - 9 Mar 2012\n\n");
    lg->writeAll("Copyright (c), 2012 - Stefano Endrizzi \n\n");
    lg->writeAll("TN -EXACT version (tmp)\n\n");
    lg->writeAll("GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
    lg->writeAll("WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
    lg->writefAll("\nWORKING DIRECTORY: %s\n",geotop::common::Variables::WORKING_DIRECTORY.c_str());
    lg->log("Using Experimental Logger");

    //TODO: remove these lines AFTER correncting all the functions that
    //use flog to log
    //8<==================================
    if (geotop::common::Variables::WORKING_DIRECTORY[strlen(geotop::common::Variables::WORKING_DIRECTORY.c_str())-1] != 47) {
        temp = geotop::common::Variables::WORKING_DIRECTORY;
        geotop::common::Variables::WORKING_DIRECTORY = temp + "/";
    }

    geotop::common::Variables::logfile = geotop::common::Variables::WORKING_DIRECTORY + logfile_name;
    flog = fopen(geotop::common::Variables::logfile.c_str(), "w");
    fprintf(flog,"STATEMENT:\n");
    fprintf(flog,"\n");
    fprintf(flog,"GEOtop 2.0.0 - 9 Mar 2012\n\n");
    fprintf(flog,"Copyright (c), 2012 - Stefano Endrizzi \n\n");
 	fprintf(flog,"TN -EXACT version (tmp)\n\n");
    fprintf(flog,"GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
    fprintf(flog,"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
    fprintf(flog,"\nWORKING DIRECTORY: %s\n",geotop::common::Variables::WORKING_DIRECTORY.c_str());
    //8<==================================

    std::clog << "\nLOGFILE: " << lg->getLogFilePath() << std::endl ;

#else
    //add "/" if it is missing
    if (geotop::common::Variables::WORKING_DIRECTORY[strlen(geotop::common::Variables::WORKING_DIRECTORY.c_str())-1] != 47) {
        temp = geotop::common::Variables::WORKING_DIRECTORY;
        geotop::common::Variables::WORKING_DIRECTORY = temp + "/";
    }

    geotop::common::Variables::logfile = geotop::common::Variables::WORKING_DIRECTORY + logfile_name;
    flog = fopen(geotop::common::Variables::logfile.c_str(), "w");

    printf("STATEMENT:\n");
    printf("\n");
    printf("GEOtop 2.0.0 'MOAB' - 9 Mar 2012\n\n");
    printf("Copyright (c), 2012 - Stefano Endrizzi \n\n");
    printf("TN -EXACT version (tmp)\n\n");
    printf("GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
    printf("WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
    printf("\nWORKING DIRECTORY: %s\n",geotop::common::Variables::WORKING_DIRECTORY.c_str());
    printf("\nLOGFILE: %s\n",geotop::common::Variables::logfile.c_str());

    fprintf(flog,"STATEMENT:\n");
    fprintf(flog,"\n");
    fprintf(flog,"GEOtop 2.0.0 - 9 Mar 2012\n\n");
    fprintf(flog,"Copyright (c), 2012 - Stefano Endrizzi \n\n");
 	fprintf(flog,"TN -EXACT version  \n\n");
    fprintf(flog,"GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>\n");
    fprintf(flog,"WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
    fprintf(flog,"\nWORKING DIRECTORY: %s\n",geotop::common::Variables::WORKING_DIRECTORY.c_str());
#endif

    //reads the parameters in __control_parameters
 
    temp = geotop::common::Variables::WORKING_DIRECTORY + program_name;
    
    boost::shared_ptr<geotop::input::ConfigStore> lConfigStore = geotop::input::ConfigStoreSingletonFactory::getInstance() ;
    const std::string lFilePath (temp) ;
    bool lParsingRes = lConfigStore->parse(lFilePath) ;
    if(not lParsingRes)
    {
#ifdef WITH_LOGGER
        lg->log("Unable to parse configuration file: " + lFilePath, geotop::logger::CRITICAL);
        exit(1);
#else
        t_error("Fatal Error! Unable to parse configuration file: " + lFilePath);
#endif
    }

    //TODO: correct this BEFORE the flog variable removal
    success = read_inpts_par(par, land, times, sl, met, IT, temp, flog);

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
#ifdef WITH_LOGGER
        lg->log("You have to assign a name to the Tzrun and wzrun files",
                geotop::logger::ERROR);
        lg->log("Geotop failed. See failing report.",
                geotop::logger::CRITICAL);
        exit(1);
#else
		t_error("Fatal Error! Geotop is closed. See failing report.");
#endif
	}	

    std::cout << "SPAR: " << fspar << " : " << geotop::common::Variables::files[fspar] << std::endl ;
    success = read_soil_parameters(geotop::common::Variables::files[fspar], IT, sl, par->soil_type_bedr_default, flog);

    geotop::common::Variables::Nl = sl->pa.getCh() - 1;

    success = read_point_file(geotop::common::Variables::files[fpointlist], IT->point_col_names, par, flog);

    geotop::common::Variables::max_time = 0.;
    geotop::common::Variables::max_time = (par->end_date - par->init_date)*86400.;//seconds

    //Recovering
    par->delay_day_recover = 0.0;
    par->n_ContRecovery = 0;
// this below should be removed (S.C. 22.06.2014)
    if(par->recover > 0){
        if(par->saving_points.size() -1< par->recover){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error: recover index higher than the length of the saving points vector");
            fclose(f);
#ifdef WITH_LOGGER
            lg->log("Recover index higher than the length of the saving points vector",
                    geotop::logger::ERROR);
            lg->log("Geotop failed. See failing report (1).",
                    geotop::logger::CRITICAL);
#else
            t_error("Fatal Error! Geotop is closed. See failing report (1).");
#endif
        }
        par->delay_day_recover = par->saving_points[par->recover];
    }

    //Continuous Recovering
/*    if (par->RunIfAnOldRunIsPresent != 1) {

        if (mio::IOUtils::fileExists(geotop::common::Variables::SuccessfulRunFile)) {

            temp = geotop::common::Variables::SuccessfulRunFile + ".old";
            rename(geotop::common::Variables::SuccessfulRunFile.c_str(), temp.c_str());
            temp = geotop::common::Variables::FailedRunFile + ".old";
            rename(geotop::common::Variables::FailedRunFile.c_str(), temp.c_str());
            
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "This simulation has successfully reached the end, you cannot recover it.\n");
            fclose(f);

#ifdef WITH_LOGGER
            lg->log("This simulation has successfully reached the end, you cannot recover it.",
                    geotop::logger::ERROR);
            lg->log("Geotop failed. See failing report.",
                    geotop::logger::CRITICAL);
            exit(1);

#else
            t_error("Fatal Error! Geotop is closed. See failing report.");
#endif
        }
        if (mio::IOUtils::fileExists(string(geotop::common::Variables::FailedRunFile))) {
            temp = geotop::common::Variables::SuccessfulRunFile + ".old";
            rename(geotop::common::Variables::SuccessfulRunFile.c_str(), temp.c_str());
            temp = geotop::common::Variables::FailedRunFile + ".old";
            rename(geotop::common::Variables::FailedRunFile.c_str(), temp.c_str());
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "This simulation has failed, you cannot recover it.\n");
            fclose(f);
#ifdef WITH_LOGGER
            lg->log("This simulation has failed, you cannot recover it.",
                    geotop::logger::ERROR);
            lg->log("Geotop failed. See failing report.",
                    geotop::logger::CRITICAL);
            exit(1);
#else
            t_error("Fatal Error! Geotop is closed. See failing report.");
#endif
        }
    }*/

    temp = geotop::common::Variables::SuccessfulRunFile + ".old";
    rename(geotop::common::Variables::SuccessfulRunFile.c_str(), temp.c_str());
    temp = geotop::common::Variables::FailedRunFile + ".old";
    rename(geotop::common::Variables::FailedRunFile.c_str(), temp.c_str());

    if (par->ContRecovery > 0 && par->recover == 0) {
        if (mio::IOUtils::fileExists(string(geotop::common::Variables::files[rtime]) + string(textfile))) {
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
        }
    }

    //Time indices
	
    par->init_date += par->delay_day_recover;
    convert_JDfrom0_JDandYear(par->init_date, &JD, &year);
    convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);

	geotop::common::Variables::i_run = geotop::common::Variables::i_run0;//Run index
		
    /****************************************************************************************************/
    /*! Reading of the Input files:                                                                     */
    /****************************************************************************************************/

    // ##################################################################################################################################
    // ##################################################################################################################################
    par->use_meteoio_cloud = false;
#ifndef USE_INTERNAL_METEODISTR
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

    geotop::common::Variables::Nr=top->Z0.getRows()-1;
    geotop::common::Variables::Nc=top->Z0.getCols()-1;

    par->total_pixel=0;
    par->total_area=0.;
    for(r=1;r<=geotop::common::Variables::Nr;r++){
        for(c=1;c<=geotop::common::Variables::Nc;c++){
            if ((long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                par->total_pixel++;
                if(par->point_sim != 1){
                    par->total_area += (geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2])/cos(top->slope[r][c]*GTConst::Pi/180.);
                }else {
                    par->total_area += (geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]);
                }
            }
        }
    }

    top->Z = find_Z_of_any_layer(top->Z0, top->slope, land->LC, sl, par->point_sim);


    top->Jdown.resize(par->total_pixel+1, 4+1);
    top->Qdown.resize(par->total_pixel+1, 4+1);
    for (i=1; i<=par->total_pixel; i++) {
        for (j=1; j<=4; j++) {
            top->Jdown[i][j] = i;
            top->Qdown[i][j] = 0.;
        }
    }

#ifdef WITH_LOGGER
    lg->logf("Valid pixels: %ld", par->total_pixel);
    lg->logf("Number of nodes: %ld",
            (geotop::common::Variables::Nl+1) * par->total_pixel);
    lg->logf("Novalue pixels: %ld",
            (geotop::common::Variables::Nr * geotop::common::Variables::Nc-par->total_pixel)
            );
#else
    fprintf(flog,"Valid pixels: %ld\n",par->total_pixel);
    fprintf(flog,"Number of nodes: %ld\n",(geotop::common::Variables::Nl+1)*par->total_pixel);
    fprintf(flog,"Novalue pixels: %ld\n",(geotop::common::Variables::Nr*geotop::common::Variables::Nc-par->total_pixel));
#endif


#ifdef WITH_LOGGER
    lg->logf("Basin area: %12g km2",
            (double)par->total_pixel*geotop::common::Variables::UV->U[1] *
            geotop::common::Variables::UV->U[2] / 1.E6);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    fprintf(flog,"Basin area: %12g km2\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#else
    fprintf(flog,"Basin area: %f km2\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#endif  //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

#ifdef WITH_LOGGER
    lg->logf("Valid pixels: %ld", par->total_pixel);
    lg->logf("Number of nodes: %ld",
            (geotop::common::Variables::Nl+1) * par->total_pixel);
    lg->logf("Novalue pixels: %ld",
            (geotop::common::Variables::Nr *
             geotop::common::Variables::Nc-par->total_pixel));
#else
    printf("\nValid pixels: %ld\n",par->total_pixel);
    printf("Number of nodes: %ld\n",(geotop::common::Variables::Nl+1)*par->total_pixel);
    printf("Novalue pixels: %ld\n",(geotop::common::Variables::Nr*geotop::common::Variables::Nc-par->total_pixel));
#endif

#ifdef WITH_LOGGER
    lg->logf("Basin area: %12g km2",
            (double)par->total_pixel * geotop::common::Variables::UV->U[1] *
            geotop::common::Variables::UV->U[2] / 1.E6);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Basin area: %12g km2\n\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#else
    printf("Basin area: %f km2\n\n",(double)par->total_pixel*geotop::common::Variables::UV->U[1]*geotop::common::Variables::UV->U[2]/1.E6);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

    /****************************************************************************************************/
    //Reading of RAIN data file,	METEO data file, 	and CLOUD data file

    num_cols = (long)nmet;

    //	meteo data
    met->data=(double***)malloc(met->st->E.size()*sizeof(double**));
    //	number of line of meteo data
    met->numlines=(long*)malloc(met->st->E.size()*sizeof(long));

    //	horizon for meteo stations
    met->horizon=(double***)malloc(met->st->E.size()*sizeof(double**));
    //	number of line in the horizon file
    met->horizonlines=(long*)malloc(met->st->E.size()*sizeof(long));
    //	line of met->data used (stored in memory to avoid from searching from the first line)

    success = read_meteostations_file(met->imeteo_stations, met->st, geotop::common::Variables::files[fmetstlist], IT->meteostations_col_names, flog);


#ifdef USE_INTERNAL_METEODISTR
    met->var=(double**)malloc((met->st->E.size()-1)*sizeof(double*));
    met->line_interp_WEB=(long*)malloc(met->st->E.size()*sizeof(long));
    met->line_interp_Bsnow=(long*)malloc(met->st->E.size()*sizeof(long));
    met->line_interp_WEB_LR=0;
    met->line_interp_Bsnow_LR=0;
#endif
    long num_met_stat=met->st->E.size()-1;
    for(size_t i=1; i <= num_met_stat; i++){
        if (met->imeteo_stations[1] != geotop::input::gDoubleNoValue) {
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

#ifdef WITH_LOGGER
                lg->log("Date Column missing in file " + temp,
                        geotop::logger::ERROR);
                lg->log("Geotop failed. See failing report (2).",
                       geotop::logger::CRITICAL);
                exit(1);
#else
                t_error("Fatal Error! Geotop is closed. See failing report (2).");
#endif
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
#ifdef WITH_LOGGER
                lg->logsf(geotop::logger::ERROR,
                        "Meteo data for station %ld contain precipitation as volume, but Linear Interpolation is set. This is not possible, the precipitation data are removed.",
                        i);
                lg ->log("Geotop failed. See failing report (3).",
                        geotop::logger::CRITICAL);
                exit(1);
#else
                t_error("Fatal Error! Geotop is closed. See failing report (3).");
#endif
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

#ifdef WITH_LOGGER
            lg->log("File meteo not in the list, meteo data not read, used default values",
                    geotop::logger::WARNING);
#else
            fprintf(flog, "Warning: File meteo not in the list, meteo data not read, used default values\n");
            printf("Warning: File meteo not in the list, meteo data not read, used default values\n");
#endif

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
        {
#ifdef WITH_LOGGER
            lg->log("Lapse rate file unavailable. Check input files. If you do not have a lapse rate file, remove its name and keywords from input file",
                    geotop::logger::WARNING);
#else
            printf("Lapse rate file unavailable. Check input files. If you do not have a lapse rate file, remove its name and keywords from input file\n");
#endif
        }
        temp = geotop::common::Variables::files[fLRs] + string(textfile);
        met->LRs = read_txt_matrix(temp, 33, 44, IT->lapserates_col_names, nlstot, &num_lines, flog);
        met->LRsnr = num_lines;
        par->LRflag=1;
#ifdef WITH_LOGGER
        lg->log("Lapse rate file read");
#else
        printf("\nLapse rate file read\n");
#endif

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

#ifdef USE_INTERNAL_METEODISTR

    //Find the first station with shortwave radiation data
    met->nstsrad=0;

    do{
        met->nstsrad++;
        a = 0;
        if ((long)met->data[met->nstsrad-1][0][iSW]  != geotop::input::gDoubleAbsent ||
           ((long)met->data[met->nstsrad-1][0][iSWb] != geotop::input::gDoubleAbsent &&
           (long)met->data[met->nstsrad-1][0][iSWd]  != geotop::input::gDoubleAbsent ))
        {
            a = 1;
        }
    } while (met->nstsrad < num_met_stat && a == 0);

#ifdef WITH_LOGGER
    if (a == 0) {
        lg->log("NO shortwave radiation measurements available",
                geotop::logger::WARNING);
    } else {
        lg->logf("Shortwave radiation measurements from station %ld\n",
                met->nstsrad);
    }
#else
    if(a==0){
        printf("WARNING: NO shortwave radiation measurements available\n");
        fprintf(flog,"WARNING: NO shortwave radiation measurements available\n");
    }else{
        printf("Shortwave radiation measurements from station %ld\n",met->nstsrad);
        fprintf(flog,"Shortwave radiation measurements from station %ld\n",met->nstsrad);
    }
#endif //WITH_LOGGER

	//Find the first station with cloud data
	met->nstcloud=0;
	do{
		met->nstcloud++;
		a=0;
        if( (long)met->data[met->nstcloud-1][0][iC]    != geotop::input::gDoubleAbsent ||
            (long)met->data[met->nstcloud-1][0][itauC] != geotop::input::gDoubleAbsent )
        {
            a = 1;
        }
	} while (met->nstcloud < met->st->Z.size() - 1 && a == 0);

#ifdef WITH_LOGGER
    if (a == 0) {
        lg->log("No cloudiness measurements available",
                geotop::logger::WARNING);
    } else {
        lg->logf("Cloudiness measurements from station %ld",
                met->nstcloud);
    }
#else
	if(a==0){
		printf("WARNING: NO cloudiness measurements available\n");
		fprintf(flog,"WARNING: NO cloudiness measurements available\n");
	}else{
		printf("Cloudiness measurements from station %ld\n",met->nstcloud);
		fprintf(flog,"Cloudiness measurements from station %ld\n",met->nstcloud);
	}
#endif //WITH_LOGGER
	
	//Find the first station with longwave radiation data
	met->nstlrad=0;
	do{
		met->nstlrad++;
		a=0;
		if( (long)met->data[met->nstlrad-1][0][iLWi]!=geotop::input::gDoubleAbsent) a=1;
	} while (met->nstlrad < met->st->Z.size() - 1 && a == 0);

#ifdef WITH_LOGGER
    if (a == 0) {
        lg->log("No longwave radiation measurements available",
                geotop::logger::WARNING);
    } else {
        lg->logf("Longwave radiation measurements from station %ld\n",
                met->nstlrad);
    }
#else
	if(a==0){
		printf("WARNING: NO longwave radiation measurements available\n");
		fprintf(flog,"WARNING: NO longwave radiation measurements available\n");
	}else{
		printf("Longwave radiation measurements from station %ld\n",met->nstlrad);
		fprintf(flog,"Longwave radiation measurements from station %ld\n",met->nstlrad);
	}
#endif //WITH_LOGGER

#endif //USE_INTERNAL_METEODISTR

    met->tau_cl_map.resize(top->Z0.getRows(),top->Z0.getCols(),geotop::input::gDoubleNoValue);

    met->tau_cl_av_map.resize(top->Z0.getRows(),top->Z0.getCols(), geotop::input::gDoubleNoValue);

    met->tau_cl_map_yes.resize(top->Z0.getRows(),top->Z0.getCols(), (short)geotop::input::gDoubleNoValue);

    met->tau_cl_av_map_yes.resize(top->Z0.getRows(),top->Z0.getCols(), (short)geotop::input::gDoubleNoValue);

    //	vector defining which meteo station has the SW radiation information

    met->st->flag_SW_meteoST.resize(met->st->Z.size(),geotop::input::gDoubleNoValue);
	
	met->st->tau_cloud_av_yes_meteoST.resize(met->st->Z.size(), geotop::input::gDoubleNoValue);

    met->st->tau_cloud_yes_meteoST.resize(met->st->Z.size(),geotop::input::gDoubleNoValue);

    met->st->tau_cloud_av_meteoST.resize(met->st->Z.size(), geotop::input::gDoubleNoValue);

    met->st->tau_cloud_meteoST.resize(met->st->Z.size(),geotop::input::gDoubleNoValue);

    //i.show details on checkpoints
#ifdef WITH_LOGGER
    lg->writeAll("\nCHECKPOINTS:\n");
    lg->writeAll("ID,r,c,Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-]\n");
#else
    fprintf(flog,"\nCHECKPOINTS:\n");
    fprintf(flog,"ID,r,c,Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-]\n");
#endif

    for(i=1;i<par->rc.getRows();i++){
        r=par->rc[i][1];
        c=par->rc[i][2];
#ifdef WITH_LOGGER
        lg->writefAll("%ld,%ld,%ld,%12g,%d,%ld,%12g,%12g,%12g\n", i, r, c,
                top->Z0[r][c], (short)land->LC[r][c], sl->type[r][c],
                top->slope[r][c], top->aspect[r][c], top->sky[r][c]);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
        fprintf(flog,"%ld,%ld,%ld,%12g,%d,%ld,%12g,%12g,%12g\n",i,r,c,top->Z0[r][c],(short)land->LC[r][c],sl->type[r][c],top->slope[r][c],top->aspect[r][c],top->sky[r][c]);
#else
        fprintf(flog,"%ld,%ld,%ld,%f,%d,%ld,%f,%f,%f\n",i,r,c,top->Z0[r][c],(short)land->LC[r][c],sl->type[r][c],top->slope[r][c],top->aspect[r][c],top->sky[r][c]);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER
    }

    //i.show meteo stations
#ifdef WITH_LOGGER
    lg->writeAll("\nMETEO STATIONS:\n");
    lg->writeAll("ID,East[m],North[m],Lat[deg],Lon[deg],Elev[m a.s.l.],Sky[deg],Stand_Time[h],WindSensHeight[m],TempSensHeight[m]\n");
#else
    fprintf(flog,"\nMETEO STATIONS:\n");
    fprintf(flog,"ID,East[m],North[m],Lat[deg],Lon[deg],Elev[m a.s.l.],Sky[deg],Stand_Time[h],WindSensHeight[m],TempSensHeight[m]\n");
#endif

    for(size_t r=1;r<met->st->E.size();r++){

#ifdef WITH_LOGGER
        lg->writefAll("%ld,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g\n",
                r, met->st->E[r], met->st->N[r], met->st->lat[r],
                met->st->lon[r], met->st->Z[r], met->st->sky[r],
                met->st->ST[r], met->st->Vheight[r], met->st->Theight[r]);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
        fprintf(flog,"%ld,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g\n",r,met->st->E[r],met->st->N[r],met->st->lat[r], met->st->lon[r],
                met->st->Z[r],met->st->sky[r],met->st->ST[r],met->st->Vheight[r],met->st->Theight[r]);
#else
        fprintf(flog,"%ld,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",r,met->st->E[r],met->st->N[r],met->st->lat[r], met->st->lon[r],
                met->st->Z[r],met->st->sky[r],met->st->ST[r],met->st->Vheight[r],met->st->Theight[r]);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER
    }

    /****************************************************************************************************/
    //read INCOMING DISCHARGE

    met->qinv = (double*)malloc(2*sizeof(double));
    met->qinv[0] = 0.;
    met->qinv[1] = 0.;

    if(geotop::common::Variables::files[fqin] != geotop::input::gStringNoValue){
        temp =geotop::common::Variables::files[fqin] + string(textfile);
        temp2.push_back("Time") ;
        temp2.push_back("Qx") ;
        met->qins = read_txt_matrix(temp, 33, 44, temp2, 2, &num_lines, flog);
        met->qinsnr = num_lines;
        par->qin = 1;
        met->qinline = 0;
#ifdef WITH_LOGGER
        lg->log("Incoming discharge file read");
#else
        printf("\nIncoming discharge file read\n");
#endif
    }else{
        par->qin = 0;
    }

    /****************************************************************************************************/
    /*! Completing the several time-indipendent input variables with the data of input-files:           */
    /****************************************************************************************************/
    /****************************************************************************************************/
    // Completing of "land" (of the type LAND):

    //	Initialize matrix of shadow

    land->shadow.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0);

    //	Check that there aren't cell with an undefined land use value
    z = 0.;
    l = 0;
    do{
        l++;
        z += sl->pa(1,jdz,l);
    } while (l < geotop::common::Variables::Nl && z < GTConst::z_transp);
    
	land->root_fraction.resize(par->n_landuses + 1, l + 1 ,0.0);

#ifdef WITH_LOGGER
    //check vegetation variable consistency
    if( (jHveg != jdHveg + jHveg - 1)               ||
        (jz0thresveg != jdz0thresveg + jHveg - 1)   ||
        (jz0thresveg2 != jdz0thresveg2 + jHveg - 1) ||
        (jLSAI != jdLSAI + jHveg - 1)               ||
        (jcf != jdcf + jHveg - 1)                   ||
        (jdecay0 != jddecay0 + jHveg - 1)           ||
        (jexpveg != jdexpveg + jHveg - 1)           ||
        (jroot != jdroot + jHveg - 1)               ||
        (jrs != jdrs + jHveg - 1))
    {
        lg->log("Vegetation variables not consistent",
                geotop::logger::CRITICAL);
        exit(1);
    }
#else
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
#endif

    //	variables used to assign vegetation properties that change with time
    num_cols = jdvegprop + 1;
    land->vegpars=(double ***)malloc(par->n_landuses*sizeof(double**));
    land->vegparv=(double **)malloc(par->n_landuses*sizeof(double*));
    land->NumlinesVegTimeDepData=(long*)malloc(par->n_landuses*sizeof(long));

    land->vegpar.resize(jdvegprop+1);

    par->vegflag.resize(par->n_landuses+1,0);

    for(i=1;i<=par->n_landuses;i++){

        if (geotop::common::Variables::files[fvegpar] != geotop::input::gStringNoValue) {   //s stands for string

            temp = namefile_i_we2(geotop::common::Variables::files[fvegpar], i);

            if (mio::IOUtils::fileExists(string(temp) + string(textfile))) {
#ifdef WITH_LOGGER
                lg->logf("There is a specific vegetation parameter file for land cover type = %ld",
                        i);
#else
                printf("There is a specific vegetation parameter file for land cover type = %ld\n",i);
#endif
                temp = namefile_i(geotop::common::Variables::files[fvegpar], i);
                land->vegpars[i-1] = read_txt_matrix_2(temp, 33, 44, num_cols, &num_lines);
                land->NumlinesVegTimeDepData[i-1] = num_lines;
                par->vegflag[i]=1;
            } else {
#ifdef WITH_LOGGER
                lg->logsf(geotop::logger::WARNING,
                        "There is NOT a specific vegetation parameter file for land cover type = %ld",
                        i);
#else
                printf("There is NOT a specific vegetation parameter file for land cover type = %ld\n",i);
#endif
            }
        }

        land->vegparv[i-1]=(double*)malloc(num_cols*sizeof(double));
        for(j=0; j<num_cols; j++){
            land->vegparv[i-1][j] = geotop::input::gDoubleNoValue;
        }

        //	z0 (convert in m)
        land->ty[i][jz0]*=0.001;

        //	find root fraction
        root(land->root_fraction.getCols(), land->ty[i][jroot], 0.0, sl->pa, land->root_fraction, i);

        //	error messages
        for(size_t l=1;l<met->st->Z.size();l++){
            if(0.001*land->ty[i][jHveg]>met->st->Vheight[l] || 0.001*land->ty[i][jHveg]>met->st->Theight[l]){
                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
#ifdef USE_DOUBLE_PRECISION_OUTPUT
                fprintf(f, "hc:%12g m, zmu:%12g m, zmt:%12g m - hc must be lower than measurement height - land cover %ld, meteo station %ld\n",
                        0.001*land->ty[i][jHveg],met->st->Vheight[l],met->st->Theight[l],i,l);
#else
                fprintf(f, "hc:%f m, zmu:%f m, zmt:%f m - hc must be lower than measurement height - land cover %ld, meteo station %ld\n",
                        0.001*land->ty[i][jHveg],met->st->Vheight[l],met->st->Theight[l],i,l);
#endif

                fclose(f);

#ifdef WITH_LOGGER
                lg->logsf(geotop::logger::ERROR,
                        "hc:%12g m, zmu:%12g m, zmt:%12g m - hc must be lower than measurement height - land cover %ld, meteo station %ld",
                        0.001*land->ty[i][jHveg],
                        met->st->Vheight[l],
                        met->st->Theight[l],
                        i,
                        l);
                lg->log("Geotop failed. See failing report (5).",
                        geotop::logger::CRITICAL);
                exit(1);
#else
                t_error("Fatal Error! Geotop is closed. See failing report (5).");
#endif
            }
        }
    }


    /****************************************************************************************************/
    /*! Filling up of the struct "channel" (of the type CHANNEL):                                        */

    /*The number of channel-pixel are counted:*/
    i=0;
    for(r=1;r<=geotop::common::Variables::Nr;r++){
        for(c=1;c<=geotop::common::Variables::Nc;c++){
            if (top->pixel_type[r][c]>=10) i++;
        }
    }

#ifdef WITH_LOGGER
    lg->logf("Channel pixels: %ld", i);
#else
    fprintf(flog,"Channel pixels: %ld\n",i);
#endif
    par->total_channel = i;

    //allocate channel vectors/matrixes
    if(i==0) i=1;

    cnet->Vout = 0.;

    cnet->r.resize(i+1,0);

    cnet->c.resize(i+1,0);

    cnet->ch.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0);

    cnet->ch_down.resize(i+1,0);

    cnet->length.resize(i+1,0.0);

    cnet->Vsup.resize(i+1,0.);
    cnet->Vsub.resize(i+1,0.);

    cnet->h_sup.resize(i+1,0.0);

    cnet->soil_type.resize(i+1, par->soil_type_chan_default);

    if(par->total_channel>1) enumerate_channels(cnet, land->LC, top->pixel_type, top->Z0, top->slope, geotop::input::gDoubleNoValue);

    cnet->ch3 = (long**)malloc((geotop::common::Variables::Nl+1)*sizeof(long*));
    for (l=0; l<=geotop::common::Variables::Nl; l++) {
        cnet->ch3[l] = (long*)malloc((i+1)*sizeof(long));
    }

    cnet->lch.resize( ((geotop::common::Variables::Nl+1)*i)+1, 2+1 , 0);

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

    top->lrc_cont.resize( (geotop::common::Variables::Nl+1)*par->total_pixel+1, 3+1, 0);

    i_lrc_cont(land->LC, top->i_cont, top->lrc_cont);

    top->j_cont=(long **)malloc((geotop::common::Variables::Nr+1)*sizeof(long*));
    for (r=1; r<=geotop::common::Variables::Nr; r++) {
        top->j_cont[r]=(long *)malloc((geotop::common::Variables::Nc+1)*sizeof(long));
        for (c=1; c<=geotop::common::Variables::Nc; c++) {
            top->j_cont[r][c] = 0;
        }
    }

    top->rc_cont.resize(par->total_pixel+1, 2+1);

    j_rc_cont(land->LC, top->j_cont, top->rc_cont);

    if(par->state_pixel == 1){
        par->jplot.resize(par->total_pixel+1, 0);

        for (i=1; i<=par->total_pixel; i++) {
            for (j=1; j<par->rc.getRows(); j++) {
                if (top->rc_cont[i][1] == par->rc[j][1] && top->rc_cont[i][2] == par->rc[j][2]) {
                    par->jplot[i] = j;
                }
            }
        }
    }

    //BEDROCK (adjusting soil properties)
    par->bedrock = 0;
    if(geotop::common::Variables::files[fbed] != geotop::input::gStringNoValue) set_bedrock(IT, sl, cnet, par, top, land->LC, flog);
    
    /****************************************************************************************************/
    /*! Completing of the initialization of SOIL structure                               */
    /****************************************************************************************************/

    sl->SS =new SoilState();
    initialize_soil_state(sl->SS, par->total_pixel, geotop::common::Variables::Nl);

    sl->VS =new StateVeg();
    initialize_veg_state(sl->VS, par->total_pixel);
    
	sl->th.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,geotop::input::gDoubleNoValue);

	sl->Ptot.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,geotop::input::gDoubleNoValue);

    if(geotop::common::Variables::files[fTav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTavsup] != geotop::input::gStringNoValue){
 
		sl->T_av_tensor.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,0.0);
    }

    if(geotop::common::Variables::files[ficeav] != geotop::input::gStringNoValue){
  
		sl->thi_av_tensor.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,0.0);
    }

    if(geotop::common::Variables::files[fliqav] != geotop::input::gStringNoValue){

        sl->thw_av_tensor.resize(geotop::common::Variables::Nl+1,par->total_pixel+1,0.0);
    }

    if(geotop::common::Variables::files[fpnet] != geotop::input::gStringNoValue){//TODO mattiu
    	sl->Pnetcum.resize(par->total_pixel+1,0.0);
    }
    if(geotop::common::Variables::files[fevap] != geotop::input::gStringNoValue){
    	sl->ETcum.resize(par->total_pixel+1,0.0);
    }//end mattiu

    sl->ET.resize(geotop::common::Variables::Nl+1,geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.);

    if (!mio::IOUtils::fileExists(string(geotop::common::Variables::files[fwt0]) + string(ascii_esri))){

        for (i=1; i<=par->total_pixel; i++) {

            r = top->rc_cont[i][1];
		    c = top->rc_cont[i][2];

            sy=sl->type[r][c];

            if ((long)IT->init_water_table_depth[sy] != geotop::input::gDoubleNoValue) {
                z = 0.;
			    sl->SS->P[0][i] = -IT->init_water_table_depth[sy]*cos(top->slope[r][c]*GTConst::Pi/180.);

                for(l=1;l<=geotop::common::Variables::Nl;l++){
                    z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
                    sl->SS->P[l][i] = sl->SS->P[0][i] + z;
                    z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
                }
            }else {
                for(l=1;l<=geotop::common::Variables::Nl;l++){
                    sl->SS->P[l][i] = sl->pa(sy,jpsi,l);
                }
            }
        }

    }else {

        meteoio_readMap(string(geotop::common::Variables::files[fwt0]), M);

        for (i=1; i<=par->total_pixel; i++) {
            r = top->rc_cont[i][1];
            c = top->rc_cont[i][2];

            sy=sl->type[r][c];

            z = 0.;
            sl->SS->P[0][i] = -M[r][c]*cos(top->slope[r][c]*GTConst::Pi/180.);
            for(l=1;l<=geotop::common::Variables::Nl;l++){
                z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
                sl->SS->P[l][i] = sl->SS->P[0][i] + z;
                z += 0.5*sl->pa(sy,jdz,l)*cos(top->slope[r][c]*GTConst::Pi/180.);
            }
        }

    }

    for (i=1; i<=par->total_pixel; i++) {

        r = top->rc_cont[i][1];
        c = top->rc_cont[i][2];

        sy=sl->type[r][c];

        for(l=1;l<=geotop::common::Variables::Nl;l++){

            sl->SS->T[l][i]=sl->pa(sy,jT,l);

            sl->Ptot[l][i] = sl->SS->P[l][i];

            sl->th[l][i] = teta_psi(sl->SS->P[l][i], 0.0, sl->pa(sy,jsat,l), sl->pa(sy,jres,l), sl->pa(sy,ja,l),
                                    sl->pa(sy,jns,l), 1-1/sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));

            th_oversat = Fmax( sl->SS->P[l][i] , 0.0 ) * sl->pa(sy,jss,l);
            sl->th[l][i] -= th_oversat;

            if(sl->SS->T[l][i]<=GTConst::Tfreezing){

                //	Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
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

// we should insert here a check to see if all the files requested are in place // 
        
		if(geotop::common::Variables::files[fTz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTzwriteend] != geotop::input::gStringNoValue) sl->Tzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fTzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTzavwriteend] != geotop::input::gStringNoValue) sl->Tzavplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fpsiztot] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fpsiztotwriteend] != geotop::input::gStringNoValue) sl->Ptotzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fpsiz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fpsizwriteend] != geotop::input::gStringNoValue) sl->Pzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fliqz] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fliqzwriteend] != geotop::input::gStringNoValue) sl->thzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[fliqzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fliqzavwriteend] != geotop::input::gStringNoValue) sl->thzavplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[ficez] != geotop::input::gStringNoValue ||geotop::common::Variables::files[ficezwriteend] != geotop::input::gStringNoValue) sl->thizplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
        if(geotop::common::Variables::files[ficezav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[ficezavwriteend] != geotop::input::gStringNoValue) sl->thizavplot.resize(par->rc.getRows(), geotop::common::Variables::Nl+1);
		
        for (i=1; i<par->rc.getRows(); i++) {
            
            r = top->rc_cont[i][1];
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

#ifdef WITH_LOGGER
	lg->logf("INPUT: printing delay_day_recover %f\n",par->delay_day_recover);
#else
	printf(" INPUT: printing delay_day_recover %f\n",par->delay_day_recover);
#endif
	
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
                sy = cnet->soil_type[i];
 
			    cnet->th[l][i] = teta_psi(Fmin(cnet->SS->P[l][i], psi_saturation(cnet->SS->thi[l][i], sl->pa(sy,jsat,l),
                                                                                 sl->pa(sy,jres,l), sl->pa(sy,ja,l), sl->pa(sy,jns,l), 1.-1./sl->pa(sy,jns,l))),
                                          cnet->SS->thi[l][i], sl->pa(sy,jsat,l), sl->pa(sy,jres,l), sl->pa(sy,ja,l),
                                          sl->pa(sy,jns,l), 1.-1./sl->pa(sy,jns,l), GTConst::PsiMin, sl->pa(sy,jss,l));

            }
        }
    }

    //	WRITE INITIAL CONDITION
    write_output_headers(met->st->Z.size(), times, wat, par, top, land, sl, egy, snow, glac);

    if(par->state_pixel == 1){
        for(j=1;j<par->rc.getRows();j++){
            if(par->output_vertical_distances == 1){
                r = par->rc[j][1];
                c = par->rc[j][2];
                cosslope = cos( Fmin(GTConst::max_slope, top->slope[r][c]) * GTConst::Pi/180.0 );
            }else {
                cosslope = 1.;
            }

            write_soil_output(j, par->IDpoint[j], par->init_date, par->init_date, JD, day, month, year, hour, minute, par->soil_plot_depths, sl, par, GTConst::PsiMin, cosslope);
        }
    }

    //	z boundary condition
    for(l=1;l<=geotop::common::Variables::Nl;l++){
         par->Zboundary -= sl->pa(1,jdz,l);
    }

    if(par->Zboundary < 0){
        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
        fprintf(f, "Z at which 0 annual temperature takes place is not lower than the soil column\n");
        fclose(f);
#ifdef WITH_LOGGER
        lg->log("Z at which 0 annual temperature takes place is not lower than the soil column",
               geotop::logger::ERROR);
        lg->log("Geotop failed. See failing report (6).",
                geotop::logger::CRITICAL);
        exit(1);
#else
        t_error("Fatal Error! Geotop is closed. See failing report (6).");
#endif
    }

    par->Zboundary *= 1.E-3;	//convert in [m]

    /****************************************************************************************************/
    /*! Initialization of the struct "egy" (of the type ENERGY):*/
    // revision performed on 24.12.2013// 

#ifdef WITH_LOGGER
    lg->log("Checking for map files undefined...");
#else
    fprintf(flog,"Checking for map files undefined...");
#endif
	
    if(par->output_surfenergy_bin == 1){
				
        if(geotop::common::Variables::files[fradnet] != geotop::input::gStringNoValue){
		    egy->Rn_mean.resize(par->total_pixel+1,0.0);
            egy->Rn.resize(par->total_pixel+1,0.0);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File NetRadiationMapFile [usually defined in output_maps/RadNet] NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File NetRadiationMapFile [usually defined in output_maps/RadNet] NOT DEFINED\n");			
#endif
		}
		
		}
        if(geotop::common::Variables::files[fradLWin] != geotop::input::gStringNoValue){
            egy->LWin_mean.resize(par->total_pixel+1, 0.0);
            egy->LWin.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File InLongwaveRadiationMapFile [usually defined in output_maps/LWin] NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File InLongwaveRadiationMapFile [usually defined in output_maps/LWin] NOT DEFINED\n");
#endif
		}		
        if((geotop::common::Variables::files[fradLW] != geotop::input::gStringNoValue)||(geotop::common::Variables::files[fradnet] != geotop::input::gStringNoValue)) {
            egy->LW_mean.resize(par->total_pixel+1,0.0);
            egy->LW.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File InLongwaveRadiationMapFile[usually defined in output_maps/LWin] NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File InLongwaveRadiationMapFile[usually defined in output_maps/LWin] NOT DEFINED\n");			
#endif
		}
	
        if((geotop::common::Variables::files[fradSW] != geotop::input::gStringNoValue)||(geotop::common::Variables::files[fradnet] != geotop::input::gStringNoValue)){
            egy->SW_mean.resize(par->total_pixel+1,0.0);
            egy->SW.resize(par->total_pixel+1);
        }else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File  with fradSWin identifier NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File  with fradSWin identifier NOT DEFINED\n");
#endif
		}
	
        if(geotop::common::Variables::files[fLE] != geotop::input::gStringNoValue){
            egy->ET_mean.resize(par->total_pixel+1,0.0);
            egy->LE.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File  SurfaceLatentHeatFluxMapFile [= maps/LE] NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File  SurfaceLatentHeatFluxMapFile [= maps/LE] NOT DEFINED\n"); 
#endif
		}
        if(geotop::common::Variables::files[fH] != geotop::input::gStringNoValue){
            egy->H_mean.resize(par->total_pixel+1,0.0);
            egy->H.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File SurfaceSensibleHeatFluxMapFile [= maps/H] NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File SurfaceSensibleHeatFluxMapFile [= maps/H] NOT DEFINED\n");
#endif
		}
        if(geotop::common::Variables::files[fG] != geotop::input::gStringNoValue){
            egy->SEB_mean.resize(par->total_pixel+1,0.0);
            egy->G.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File  with fG identifier  NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File  with fG identifier  NOT DEFINED\n");
#endif
		}	
        if(geotop::common::Variables::files[fTs] != geotop::input::gStringNoValue){
            egy->Ts_mean.resize(par->total_pixel+1,0.0);
            egy->Ts.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File  with fTs identifier  NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog,"File  with fTs identifier  NOT DEFINED\n");
#endif
			
		}	
        if(geotop::common::Variables::files[fradSWin] != geotop::input::gStringNoValue){
            egy->Rswdown_mean.resize(par->total_pixel+1,0.0);
            egy->SWin.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File  with fradSWin identifier  NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog, "File  with fradSWin identifier  NOT DEFINED\n");
#endif
		}	
        if(geotop::common::Variables::files[fradSWinbeam] != geotop::input::gStringNoValue){
			
            egy->Rswbeam_mean.resize(par->total_pixel+1,0.0);
            egy->SWinb.resize(par->total_pixel+1);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File  with fradSwinbeam identifier  NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog, "File  with fradSwinbeam identifier  NOT DEFINED\n");
#endif
		}
	
        if(geotop::common::Variables::files[fshadow] != geotop::input::gStringNoValue){
            egy->nDt_shadow.resize(par->total_pixel+1,0.0);
            egy->nDt_sun.resize(par->total_pixel+1,0.0);
            egy->shad.resize(par->total_pixel+1,0.0);
        }
		else{
			count_file_missing++;
#ifdef WITH_LOGGER
			lg->log("File ShadowFractionTimeMapFile [= maps/Shadow???] NOT DEFINED",
                    geotop::logger::WARNING);
#else
			fprintf(flog, "File ShadowFractionTimeMapFile [= maps/Shadow???] NOT DEFINED\n");
#endif

		}
	
	if (count_file_missing >0){
#ifdef WITH_LOGGER
		lg->logsf(geotop::logger::WARNING,
                "%d mapfiles undefined: see above for names of missing files",
                count_file_missing);
#else
		fprintf(flog, "Warning: %d mapfiles undefined: see above for names of missing files\n",count_file_missing);
#endif
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


#ifdef WITH_LOGGER
    lg->logf("Soil water evaporates from the first %ld layers",egy->soil_evap_layer_bare.size() - 1);
    lg->logf("Soil water transpires from the first %ld layers",egy->soil_transp_layer.size() - 1);
#else
    printf("Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare.size()-1);
    printf("Soil water transpires from the first %ld layers\n",egy->soil_transp_layer.size()-1);
    fprintf(flog,"Soil water evaporates from the first %ld layers\n",egy->soil_evap_layer_bare.size()-1);
    fprintf(flog,"Soil water transpires from the first %ld layers\n",egy->soil_transp_layer.size()-1);
#endif

    /****************************************************************************************************/
    /*! Completing of the struct "water" (of the type WATER) */

    wat->Voutlandsub = 0.;
    wat->Voutlandsup = 0.;
    wat->Voutbottom = 0.;

    /* Initialization of wat->Pnet (liquid precipitation that reaches the sl surface in mm):*/
    //	wat->Pnet=new_doublematrix(Nr,Nc);
    //	initialize_doublematrix(wat->Pnet,0.0);
    wat->Pnet.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.0);

    wat->HN.resize(geotop::common::Variables::Nr+1,geotop::common::Variables::Nc+1,0.0);//TODO mattiu
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
#ifdef WITH_LOGGER
        lg->logf("Initial condition on snow depth from file %s",
                geotop::common::Variables::files[fsn0].c_str());
#else
        printf("Initial condition on snow depth from file %s\n",geotop::common::Variables::files[fsn0].c_str());
#endif
        //	M=read_map(2,geotop::common::Variables::files[fsn0], land->LC, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
        meteoio_readMap(string(geotop::common::Variables::files[fsn0]), M);

        for (r=1; r<=geotop::common::Variables::Nr; r++) {
            for (c=1; c<=geotop::common::Variables::Nc; c++) {
                //	snow->S->Dzl->co[1][r][c] = M->co[r][c];
                snow->S->Dzl[1][r][c] = M[r][c];
            }
        }
        //	free_doublematrix(M);

#ifdef WITH_LOGGER
        lg->logf("Initial condition on snow water equivalent from file %s",
                geotop::common::Variables::files[fswe0].c_str());
#else
        printf("Initial condition on snow water equivalent from file %s\n",geotop::common::Variables::files[fswe0].c_str());
#endif
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
#ifdef WITH_LOGGER
        lg->logf("Initial condition on snow depth from file %s",
                geotop::common::Variables::files[fsn0].c_str());
#else
        printf("Initial condition on snow depth from file %s\n",geotop::common::Variables::files[fsn0].c_str());
#endif
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
#ifdef WITH_LOGGER
        lg->logf("Initial condition on snow water equivalent from file %s",
                geotop::common::Variables::files[fswe0].c_str());
#else
        printf("Initial condition on snow water equivalent from file %s\n",geotop::common::Variables::files[fswe0].c_str());
#endif
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
#ifdef WITH_LOGGER
        lg->logf("Snow age initial condition from file %s",
                geotop::common::Variables::files[fsnag0 + 1].c_str());
#else
        printf("Snow age initial condition from file %s\n",geotop::common::Variables::files[fsnag0 + 1].c_str());
#endif
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

        if(geotop::common::Variables::files[fHN] != geotop::input::gStringNoValue){//TODO mattiu
			snow->HNcum.resize(par->total_pixel+1,0.0);
			//snow->yes.resize(par->total_pixel+1,0.0);//boh
		}//end mattiu
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
#ifdef WITH_LOGGER
                    lg->logsf(geotop::logger::ERROR,
                            "Error: negative initial snow depth %12g or snow water equivalent %12g",
                            D, SWE);
                    lg->log("Geotop failed. See failing report (7).",
                            geotop::logger::CRITICAL);
                    exit(1);
#else
                    t_error("Fatal Error! Geotop is closed. See failing report (7).");
#endif

                }else if(D<1.E-5 && SWE>1.E-5){
                    f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
#ifdef USE_DOUBLE_PRECISION_OUTPUT
                    fprintf(f, "Error: Initial snow water equivalent %12g > 0 and initial snow depth %12g\n",SWE,D);
#else
                    fprintf(f, "Error: Initial snow water equivalent %e > 0 and initial snow depth %e\n",SWE,D);
#endif

                    fclose(f);
#ifdef WITH_LOGGER
                    lg->logsf(geotop::logger::ERROR,
                            "Initial snow water equivalent %12g > 0 and initial snow depth %12g",
                            SWE, D);
                    lg->log("Geotop failed. See failing report (8).",
                            geotop::logger::CRITICAL);
                    exit(1);
#else
                    t_error("Fatal Error! Geotop is closed. See failing report (8).");
#endif

                }else if(D>1.E-5 && SWE<1.E-5){
                    f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
#ifdef USE_DOUBLE_PRECISION_OUTPUT
                    fprintf(f, "Error: Initial snow depth %12g > 0 and initial snow water equivalent %12g\n",D,SWE);
#else
                    fprintf(f, "Error: Initial snow depth %e > 0 and initial snow water equivalent %e\n",D,SWE);
#endif
                    fclose(f);
#ifdef WITH_LOGGER
                    lg->logsf(geotop::logger::ERROR,
                            "Initial snow depth %12g > 0 and initial snow water equivalent %12g\n",
                            D,SWE);
                    lg->log("Geotop failed. See failing report (9).",
                            geotop::logger::CRITICAL);
                    exit(1);
#else
                    t_error("Fatal Error! Geotop is closed. See failing report (9).");
#endif

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

    f = fopen(geotop::common::Variables::logfile.c_str(), "a");
	for(r=1;r<=geotop::common::Variables::Nr;r++){
    	for(c=1;c<=geotop::common::Variables::Nc;c++){
    		if( (long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
    			snow_layer_combination(par->alpha_snow, r, c, snow->S, 0., par->inf_snow_layers, par->max_weq_snow, maxSWE, f);
    		}
    	}
    }
    fclose(f);

    /****************************************************************************************************/
    /*! Initialization of the struct "glac" (of the type GLACIER):*/

    /***************************************************************************************************/
    /*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
    if( par->point_sim!=1 &&geotop::common::Variables::files[fgl0] != geotop::input::gStringNoValue && par->max_glac_layers==0){
#ifdef WITH_LOGGER
        lg->log("Glacier map present, but glacier represented with 0 layers",
                geotop::logger::WARNING);
#else
        printf("Warning: Glacier map present, but glacier represented with 0 layers\n");
        fprintf(flog,"Warning: Glacier map present, but glacier represented with 0 layers\n");
#endif
    }

    if(par->max_glac_layers==0 && IT->Dglac0>0){
#ifdef WITH_LOGGER
        lg->log("You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.",
                geotop::logger::WARNING);
#else
        printf("\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
        fprintf(flog,"\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, but you assigned a value of the glacier depth. The latter will be ignored.\n");
#endif
    }

    //If the max number of glacier layers is greater than 1, the matrices (or tensors) lnum, Dzl. w_liq, w_ice, T and print matrices are defined, according to the respective flags
    if(par->max_glac_layers>0){

        if( par->point_sim!=1 &&geotop::common::Variables::files[fgl0] != geotop::input::gStringNoValue ){
#ifdef WITH_LOGGER
        lg->logf("Glacier initial condition from file %s",
                geotop::common::Variables::files[fgl0+1].c_str());
#else
            printf("Glacier initial condition from file %s\n",geotop::common::Variables::files[fgl0+1].c_str());
            fprintf(flog,"Glacier initial condition from file %s\n",geotop::common::Variables::files[fgl0+1].c_str());
#endif
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
#ifdef WITH_LOGGER
                        lg->log("Negative glacier data",
                                geotop::logger::ERROR);
                        lg->log("Geotop failed. See failing report (10).",
                                geotop::logger::CRITICAL);
                        exit(1);
#else
                        t_error("Fatal Error! Geotop is closed. See failing report (10).");
#endif

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
#ifdef WITH_LOGGER
    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();
#endif

    meteoio_readDEM(top->Z0);

    //	reading TOPOGRAPHY
#ifdef WITH_LOGGER
    flag = file_exists(fdem);
#else
    flag = file_exists(fdem, flog);
#endif
    if(flag == 1){

        //filtering
        M.resize(top->Z0.getRows(),top->Z0.getCols());
        multipass_topofilter(par->lowpass, top->Z0, M, geotop::input::gDoubleNoValue, 1);
        top->Z0=M;

        //	calculate East and North
        top->East.resize(top->Z0.getRows(), top->Z0.getCols());
        top->North.resize(top->Z0.getRows(), top->Z0.getCols());
        for (r=1; r<top->Z0.getRows(); r++) {
            for (c=1; c<top->Z0.getCols(); c++) {
                top->East[r][c] = geotop::common::Variables::UV->U[4] + (c-0.5)*geotop::common::Variables::UV->U[2];
                top->North[r][c] = geotop::common::Variables::UV->U[3] + ((top->Z0.getRows()-1)-(r-0.5))*geotop::common::Variables::UV->U[1];
            }
        }

    }else{

        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
        fprintf(f, "Error: It is impossible to proceed without giving the digital elevation model\n");
        fclose(f);
#ifdef WITH_LOGGER
        lg->log("It is impossible to proceed without giving the digital elevation model",
                geotop::logger::ERROR);
        lg->log("Geotop failed. See failing report (11).",
                geotop::logger::CRITICAL);
        exit(1);
#else
        t_error("Fatal Error! Geotop is closed. See failing report (11).");
#endif

    }

    //	reading LAND COVER TYPE
#ifdef WITH_LOGGER
    flag = file_exists(flu);
#else
    flag = file_exists(flu, flog);
#endif
    if(flag == 1){
        meteoio_readMap(string(geotop::common::Variables::files[flu]), land->LC);

        //	Check borders
        for(r=1;r<land->LC.getRows();r++){
            land->LC[r][1]=geotop::input::gDoubleNoValue;
            land->LC[r][land->LC.getCols()-1]=geotop::input::gDoubleNoValue;
        }
        for(c=1;c<land->LC.getCols();c++){
            land->LC[1][c]=geotop::input::gDoubleNoValue;
            land->LC[land->LC.getRows()-1][c]=geotop::input::gDoubleNoValue;
        }
        for(r=1;r<land->LC.getRows();r++){
            for(c=1;c<land->LC.getCols();c++){
                if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue) {
                    if ((long)land->LC[r][c] < 1 || (long)land->LC[r][c] > par->n_landuses){
                        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                        fprintf(f, "Error: It is not possible to assign Value < 1 or > n_landuses to the land cover type\n");
                        fclose(f);
#ifdef WITH_LOGGER
                        lg->log("It is not possible to assign Value < 1 or > n_landuses to the land cover type",
                                geotop::logger::ERROR);
                        lg->log("Geotop failed. See failing report (12).",
                                geotop::logger::CRITICAL);
                        exit(1);
#else
                        t_error("Fatal Error! Geotop is closed. See failing report (12).");
#endif
                    }
                }
            }
        }

        //Land use is the official mask
        for(r=1;r< land->LC.getRows();r++){
            for(c=1;c<land->LC.getCols();c++){
                if((long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                    if((long)top->Z0[r][c]==geotop::input::gDoubleNoValue){
#ifdef WITH_LOGGER
                        lg->log("Land use mask include DTM novalue pixels",
                                geotop::logger::WARNING);
#else
                        printf("ERROR Land use mask include DTM novalue pixels");
                        printf("\nr:%ld c:%ld Z:%f landuse:%f\n",r,c,top->Z0[r][c],land->LC[r][c]);
                        land->LC[r][c]=geotop::input::gDoubleNoValue;
                        printf("LANDUSE set at novalue where DTM is not available\n");
#endif
                    }
                }
            }
        }

    }else{

        //	Write land->LC (land cover)
#ifdef WITH_LOGGER
        lg->log("Land cover type assumed to be always 1");
#else
        printf("Land cover type assumed to be always 1\n");
#endif
        copydoublematrix_const(1.0, top->Z0,land->LC ,geotop::input::gDoubleNoValue);

        for(r=1;r<land->LC.getRows();r++){
            land->LC[r][1]=geotop::input::gDoubleNoValue;
            land->LC[r][land->LC.getCols()-1]=geotop::input::gDoubleNoValue;
        }
        for(c=1;c<land->LC.getCols();c++){
            land->LC[1][c]=geotop::input::gDoubleNoValue;
            land->LC[land->LC.getRows()-1][c]=geotop::input::gDoubleNoValue;
        }
    }

#ifdef WITH_LOGGER
	lg->logf("par->state_pixel=%ld\n",
            par->state_pixel);
#else
	printf("par->state_pixel=%ld\n",par->state_pixel);
#endif
	
    if(par->state_pixel == 1){
        par->rc.resize(par->chkpt.getRows(),2+1);
        par->IDpoint.resize(par->chkpt.getRows());
		
        for(i=1;i< par->chkpt.getRows();i++){
            par->rc[i][1]=row(par->chkpt[i][ptY], top->Z0.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            par->rc[i][2]=col(par->chkpt[i][ptX], top->Z0.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

            if (par->rc[i][1] == geotop::input::gDoubleNoValue || par->rc[i][2] == geotop::input::gDoubleNoValue) {
#ifdef WITH_LOGGER
                lg->logsf(geotop::logger::ERROR,
                        "Point #%4ld is out of the domain",
                        i);
#else
                printf("Point #%4ld is out of the domain",i);
                fprintf(flog, "Point #%4ld is out of the domain",i);
                fclose(flog);
#endif

                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                fprintf(f, "Point #%4ld is out of the domain",i);
                fclose(f);

#ifdef WITH_LOGGER
                lg->log("Geotop failed. See failing report.",
                        geotop::logger::CRITICAL);
                exit(1);
#else
                t_error("Fatal Error! Geotop is closed. See failing report.");
#endif
            }

            if((long)land->LC[par->rc[i][1]][par->rc[i][2]]==geotop::input::gDoubleNoValue){
#ifdef WITH_LOGGER
                lg->logsf(geotop::logger::ERROR,
                        "Point #%4ld corresponds to NOVALUE pixel",
                        i);
#else
                printf("Point #%4ld corresponds to NOVALUE pixel",i);
                fprintf(flog, "Point #%4ld corresponds to NOVALUE pixel",i);
                fclose(flog);
#endif

                f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                fprintf(f, "Point #%4ld corresponds to NOVALUE pixel",i);
                fclose(f);
#ifdef WITH_LOGGER
                lg->log("Geotop failed. See failing report.",
                        geotop::logger::CRITICAL);
                exit(1);
#else
                t_error("Fatal Error! Geotop is closed. See failing report.");
#endif
            }
			
#ifdef WITH_LOGGER
			lg->writefAll("i:%ld %f\n",i,par->chkpt[i][ptID]);
#else
			printf("i:%ld %f\n",i,par->chkpt[i][ptID]);
#endif
			
            if ((long)par->chkpt[i][ptID]!=geotop::input::gDoubleNoValue) {
                par->IDpoint[i]=(long)par->chkpt[i][ptID];
#ifdef WITH_LOGGER
				lg->writefAll("A i:%ld %ld\n",i,par->IDpoint[i]);
#else
				printf("A i:%ld %ld\n",i,par->IDpoint[i]);
#endif
            }else {
                par->IDpoint[i]=i;
#ifdef WITH_LOGGER
				lg->writefAll("B i:%ld %ld\n",i,par->IDpoint[i]);
#else
				printf("B i:%ld %ld\n",i,par->IDpoint[i]);
#endif
            }

        }
    }


    /****************************************************************************************************/

    //reading SKY VIEW FACTOR
#ifdef WITH_LOGGER
    flag = file_exists(fsky);
#else
    flag = file_exists(fsky, flog);
#endif
    if(flag == 1){
        meteoio_readMap(string(geotop::common::Variables::files[fsky]), top->sky);
    }else{/*The sky view factor file "top->sky" must be calculated*/
        top->sky.resize(top->Z0.getRows(),top->Z0.getCols());
        if (par->sky == 0) {
            //	initialize_doublematrix(top->sky, 1.);
            top->sky.resize(top->Z0.getRows(),top->Z0.getCols(),1.);
        }else {
            curv.resize(top->Z0.getRows(),top->Z0.getCols());
            nablaquadro_mask(top->Z0, curv, geotop::common::Variables::UV->U, geotop::common::Variables::UV->V);
            sky_view_factor(top->sky, 36, geotop::common::Variables::UV, top->Z0, curv, geotop::input::gDoubleNoValue);
        }
    }

    /****************************************************************************************************/

    //reading DELAY
#ifdef WITH_LOGGER
    flag = file_exists(fdelay);
#else
    flag = file_exists(fdelay, flog);
#endif
    if(flag == 1){
        meteoio_readMap(string(geotop::common::Variables::files[fdelay]), land->delay);
    }else{
        land->delay.resize(top->Z0.getRows(),top->Z0.getCols(),0);
    }

    /****************************************************************************************************/

    //reading SOIL MAP
#ifdef WITH_LOGGER
    flag = file_exists(fsoil);
#else
    flag = file_exists(fsoil, flog);
#endif
    if(flag == 1){
        meteoio_readMap(string(geotop::common::Variables::files[fsoil]), M);

        copylong_doublematrix(sl->type, M);
        for(r=1;r<land->LC.getRows();r++){
            for(c=1;c<land->LC.getCols();c++){
                if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue) {
                    if (sl->type[r][c] < 1 || sl->type[r][c] > par->nsoiltypes){
                        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                        fprintf(f, "Error: It is not possible to assign Value < 1 or > nsoiltypes to the soil type map");
                        fclose(f);
#ifdef WITH_LOGGER
                        lg->log("It is not possible to assign Value < 1 or > nsoiltypes to the soil type map",
                                geotop::logger::ERROR);
                        lg->log("GEOtop failed. See failing report (13).",
                                geotop::logger::CRITICAL);
                        exit(1);
#else
                        t_error("Fatal Error! GEOtop is closed. See failing report (13).");
#endif
                    }
                }
            }
        }
    }else{
        copydoublematrix_const(par->soil_type_land_default, land->LC, M, geotop::input::gDoubleNoValue);

        copylong_doublematrix(sl->type,M);
    }

    /****************************************************************************************************/
    //SLOPE
    top->dzdE.resize(land->LC.getRows(), land->LC.getCols());
    top->dzdN.resize(land->LC.getRows(), land->LC.getCols());

    find_slope(geotop::common::Variables::UV->U[1], geotop::common::Variables::UV->U[2], top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue);

#ifdef WITH_LOGGER
    flag = file_exists(fslp);
#else
    flag = file_exists(fslp, flog);
#endif
    if(flag == 1){
        meteoio_readMap(string(geotop::common::Variables::files[fslp]), top->slope);
    }else{
        find_max_slope(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue, top->slope);
    }

    find_min_max(top->slope, geotop::input::gDoubleNoValue, &max, &min);

#ifdef WITH_LOGGER
    lg->logf("Slope Min:%12g (%12g deg) Max:%12g (%12g deg)",
            tan(min*GTConst::Pi/180.), min,
            tan(max*GTConst::Pi/180.), max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
    fprintf(flog,"Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#else
    printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
    fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

    /****************************************************************************************************/
    //ASPECT

#ifdef WITH_LOGGER
    flag = file_exists(fasp);
#else
    flag = file_exists(fasp, flog);
#endif
    if(flag == 1){
        meteoio_readMap(string(geotop::common::Variables::files[fasp]), top->aspect);
    }else{
        find_aspect(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue,top->aspect);
    }
	if(flag >= 0) write_map(geotop::common::Variables::files[fasp], 0, par->format_out, top->aspect, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);


    /****************************************************************************************************/
    //curvature

    top->curvature1.resize(top->Z0.getRows(),top->Z0.getCols());
    top->curvature2.resize(top->Z0.getRows(),top->Z0.getCols());
    top->curvature3.resize(top->Z0.getRows(),top->Z0.getCols());
    top->curvature4.resize(top->Z0.getRows(),top->Z0.getCols());

    //filtering
    M.resize(top->Z0.getRows(),top->Z0.getCols());
    multipass_topofilter(par->lowpass_curvatures, top->Z0, M, geotop::input::gDoubleNoValue, 1);
    curvature(geotop::common::Variables::UV->U[1], geotop::common::Variables::UV->U[2], M, top->curvature1, top->curvature2, top->curvature3, top->curvature4, geotop::input::gDoubleNoValue);

    if(geotop::common::Variables::files[fcurv] != geotop::input::gStringNoValue){
        temp =geotop::common::Variables::files[fcurv] + string("N-S");
        write_map(temp, 0, par->format_out, top->curvature1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

        temp =geotop::common::Variables::files[fcurv] + string("W-E");
        write_map(temp, 0, par->format_out, top->curvature2, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

        temp =geotop::common::Variables::files[fcurv] + string("NW-SE");
        write_map(temp, 0, par->format_out, top->curvature3, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

        temp =geotop::common::Variables::files[fcurv] + string("NE-SW");
        write_map(temp, 0, par->format_out, top->curvature4, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
    }

    find_min_max(top->curvature1, geotop::input::gDoubleNoValue, &max, &min);
#ifdef WITH_LOGGER
    lg->logf("Curvature N-S Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature N-S Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature N-S Min:%12g  Max:%12g \n",min,max);
#else
    printf("Curvature N-S Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature N-S Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

    find_min_max(top->curvature2, geotop::input::gDoubleNoValue, &max, &min);
#ifdef WITH_LOGGER
    lg->logf("Curvature W-E Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature W-E Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature W-E Min:%12g  Max:%12g \n",min,max);
#else
    printf("Curvature W-E Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature W-E Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

    find_min_max(top->curvature3, geotop::input::gDoubleNoValue, &max, &min);
#ifdef WITH_LOGGER
    lg->logf("Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
#else
    printf("Curvature NW-SE Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature NW-SE Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

    find_min_max(top->curvature4, geotop::input::gDoubleNoValue, &max, &min);

#ifdef WITH_LOGGER
    lg->logf("Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
    printf("Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
    fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
#else
    printf("Curvature NE-SW Min:%f  Max:%f \n",min,max);
    fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

    /****************************************************************************************************/
    //Channel network (in top->pixel_type)

    //pixel type = 0 land pixel (if it is on the border, the border is impermeable, water is free only on the surface)
    //pixel type = 1 or 2 land pixel (it it is on the border, the border is permeable above an user-defined elevation in the saturated part, weir-wise)
    //pixel type = 10 channel pixel (if it is on the border, the border is impermeable, water is free only on the surface)
    //pixel type = -1 land pixel where an incoming discharge from outside is considered (as rain)

#ifdef WITH_LOGGER
    flag = file_exists(fnet);
#else
    flag = file_exists(fnet, flog);
#endif
    if(flag == 1){
        meteoio_readMap(string(geotop::common::Variables::files[fnet]), M);

        copyshort_doublematrix(top->pixel_type, M);

        cont = 0;
        for(r=1;r<top->Z0.getRows();r++){
            for(c=1;c<top->Z0.getCols();c++){
                if((long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                    if(top->pixel_type[r][c]!=0 && top->pixel_type[r][c]!=1 && top->pixel_type[r][c]!=2 && top->pixel_type[r][c]!=10){
                        f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                        fprintf(f, "Error: Only the following values are admitted in the network map: 0, 1, 10\n");
                        fclose(f);
#ifdef WITH_LOGGER
                        lg->log("Only the following values are admitted in the network map: 0, 1, 10",
                                geotop::logger::ERROR);
                        lg->log("Geotop failed. See failing report (14).",
                                geotop::logger::CRITICAL);
                        exit(1);
#else
                        t_error("Fatal Error! Geotop is closed. See failing report (14).");
#endif
                    }
                    if(top->pixel_type[r][c]==10) cont++;
                }
            }
        }

#ifdef WITH_LOGGER
        lg->logf("Channel networks has %ld pixels set to channel",cont);
#else
        printf("Channel networks has %ld pixels set to channel\n",cont);
#endif

        if(flag >= 0) write_map(geotop::common::Variables::files[fnet], 1, par->format_out, M, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

    }else{

        top->pixel_type.resize(land->LC.getRows(), land->LC.getCols(), 0);

    }

    /****************************************************************************************************/

    //border
    top->is_on_border.resize(land->LC.getRows(), land->LC.getCols());
    for(r=1;r<land->LC.getRows();r++){
        for(c=1;c<land->LC.getCols();c++){
            if ( (long)land->LC[r][c]!=geotop::input::gDoubleNoValue){
                top->is_on_border[r][c] = is_boundary(r, c, land->LC, geotop::input::gDoubleNoValue);
            }else{
                top->is_on_border[r][c] = -1;
            }
        }
    }

    //count the pixels having pixel_type = 1, 2 or -1
    cont = 0;

    for(r=1;r<top->Z0.getRows();r++){
        for(c=1;c<top->Z0.getCols();c++){
            if(top->is_on_border[r][c]==1){
                if (top->pixel_type[r][c] == -1 || top->pixel_type[r][c] == 1 || top->pixel_type[r][c] == 2) cont ++;
            }
        }
    }

    top->BC_counter.resize(top->Z0.getRows(), top->Z0.getCols(), 0);

    if (cont > 0) {
        top->BC_DepthFreeSurface.resize(cont+1);

        cont = 0;
        for(r=1;r<top->Z0.getRows();r++){
            for(c=1;c<top->Z0.getCols();c++){
                if(top->is_on_border[r][c]==1){
                    if (top->pixel_type[r][c] == -1 || top->pixel_type[r][c] == 1 || top->pixel_type[r][c] == 2){
                        cont ++;
                        top->BC_counter[r][c] = cont;
                        top->BC_DepthFreeSurface[cont] = par->DepthFreeSurface; //[mm]
                    }
                }
            }
        }
    }else {
        top->BC_DepthFreeSurface.resize(2,geotop::input::gDoubleNoValue);
    }

    //bedrock: the map is read inside the set_bedrock function
//    flag = file_exists(fbed, flog);
//    if(flag == 1){
//    	GeoMatrix<double> B;
//    	//IT->bed=read_map(2, files[fbed], land->LC, UV, (double)number_novalue);
//    	 meteoio_readMap(string(geotop::common::Variables::files[fbed]), B);
//    }

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void read_optionsfile_point(Par *par, Topo *top, Land *land, Soil *sl, Times *times, InitTools *IT, FILE *flog){

    long i, r, c, num_lines;
    GeoMatrix<double> Q, P, R, S, T, Z, LU;
    GeoMatrix<short> curv;
    short read_dem, read_lu, read_soil, read_sl, read_as, read_sk, read_bed, read_curv, flag, coordinates;
    string temp;
    double min, max;
    FILE *f;
#ifdef WITH_LOGGER
    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();
#endif

    //4. CALCULATE TOPOGRAPHIC PROPERTIES
    //check if there are point coordinates
    coordinates = 1;
    for(i=1;i<par->chkpt.getRows();i++){
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
    for(i=1;i< par->chkpt.getRows();i++){
        if((long)par->chkpt[i][ptLC]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptSY]==geotop::input::gDoubleNoValue ||
                (long)par->chkpt[i][ptS]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptA]==geotop::input::gDoubleNoValue ||
                (long)par->chkpt[i][ptCNS]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptCWE]==geotop::input::gDoubleNoValue ||
                (long)par->chkpt[i][ptCNwSe]==geotop::input::gDoubleNoValue || (long)par->chkpt[i][ptCNeSw]==geotop::input::gDoubleNoValue){
            read_dem=1;
        }
    }
    if(read_dem == 1 && coordinates == 0){
#ifdef WITH_LOGGER
        lg->log("Not possible to read from dem because at least one point has no coordinates",
                geotop::logger::WARNING);
#else
        printf("Warning: Not possible to read from dem because at least one point has no coordinates\n");
        fprintf(flog,"Warning: Not possible to read from dem because at least one point has no coordinates\n");
#endif
        read_dem = 0;
    }
    if(read_dem==1){
#ifdef WITH_LOGGER
        flag = file_exists(fdem);
#else
        flag = file_exists(fdem, flog);
#endif
        if(flag == 1){
#ifdef WITH_LOGGER
        lg->logsf(geotop::logger::WARNING,
                "Dem file %s present",
                geotop::common::Variables::files[fdem+1].c_str());
#else
            printf("Warning: Dem file %s present\n",geotop::common::Variables::files[fdem+1].c_str());
            fprintf(flog,"Warning: Dem file %s present\n",geotop::common::Variables::files[fdem+1].c_str());
#endif

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
#ifdef WITH_LOGGER
        lg->log("Dem file not present",
                geotop::logger::WARNING);
#else
            printf("Warning: Dem file not present\n");
            fprintf(flog,"Warning: Dem file not present\n");
#endif

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
#ifdef WITH_LOGGER
        flag = file_exists(flu);
#else
        flag = file_exists(flu, flog);
#endif
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
#ifdef WITH_LOGGER
        lg->log("Landuse file not present, uniform cover considered",
                geotop::logger::WARNING);
#else
            printf("Warning: Landuse file not present, uniform cover considered\n");
            fprintf(flog,"Warning: Landuse file not present, uniform cover considered\n");
#endif
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
#ifdef WITH_LOGGER
        flag = file_exists(fsoil);
#else
        flag = file_exists(fsoil, flog);
#endif
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
#ifdef WITH_LOGGER
        lg->log("Soiltype file not present",
                geotop::logger::WARNING);
#else
            printf("Warning: Soiltype file not present\n");
            fprintf(flog,"Warning: Soiltype file not present\n");
#endif
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
#ifdef WITH_LOGGER
        flag = file_exists(fslp);
#else
        flag = file_exists(fslp, flog);
#endif
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
#ifdef WITH_LOGGER
        lg->log("Slopes file not present",
                geotop::logger::WARNING);
#else
                printf("Warning: Slopes file not present\n");
                fprintf(flog,"Warning: Slopes file not present\n");
#endif
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
#ifdef WITH_LOGGER
        lg->writefAll("Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
        printf("Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
        fprintf(flog,"Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#else
        printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
        fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*GTConst::Pi/180.),min,tan(max*GTConst::Pi/180.),max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

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
#ifdef WITH_LOGGER
        flag = file_exists(fasp);
#else
        flag = file_exists(fasp, flog);
#endif
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
#ifdef WITH_LOGGER
        lg->log("Aspect file not present",
                geotop::logger::WARNING);
#else
                printf("Warning: Aspect file not present\n");
                fprintf(flog,"Warning: Aspect file not present\n");
#endif
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
#ifdef WITH_LOGGER
        flag = file_exists(fsky);
#else
        flag = file_exists(fsky, flog);
#endif
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
#ifdef WITH_LOGGER
        lg->log("Sky view factor file not present",
                geotop::logger::WARNING);
#else
                printf("Warning: Sky view factor file not present\n");
                fprintf(flog,"Warning: Sky view factor file not present\n");
#endif
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

    //f2. bedrock file
	read_bed=0;cout << "par->chkpt.getRows()=" << par->chkpt.getRows() << " par->chkpt[i][ptBED]=" << par->chkpt[1][ptBED] <<  endl;
	for(i=1;i<par->chkpt.getRows();i++){
		if((long)par->chkpt[i][ptBED]==geotop::input::gDoubleNoValue) read_bed=1;
	}
	if(read_bed==1 && coordinates==0) read_bed=0;
	if(read_bed==1){
#ifdef WITH_LOGGER
		flag = file_exists(fbed);
#else
		flag = file_exists(fbed, flog);
#endif
		if(flag == 1){
				meteoio_readMap(string(geotop::common::Variables::files[fbed]), P);
		}else{
#ifdef WITH_LOGGER
        lg->log("Bedrock depth file not present",
                geotop::logger::WARNING);
#else
			printf("Warning: Bedrock depth file not present\n");
			fprintf(flog,"Warning: Bedrock view factor file not present\n");
#endif
			read_bed=0;
		}
	}

	if(read_bed==1){
		for(i=1;i< par->chkpt.getRows();i++){
			if((long)par->chkpt[i][ptBED]==geotop::input::gDoubleNoValue){
				r=row(par->chkpt[i][ptY], P.getRows()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
				c=col(par->chkpt[i][ptX], P.getCols()-1, geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
				par->chkpt[i][ptBED]=P[r][c];
			}
		}
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
#ifdef WITH_LOGGER
        lg->log("Dem file is not present, and therefore it is not possible to calculate curvature",
                geotop::logger::WARNING);
#else
            printf("Warning: Dem file is not present, and therefore it is not possible to calculate curvature\n");
            fprintf(flog,"Warning: Dem file is not present, and therefore it is not possible to calculate curvature\n");
#endif
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
#ifdef WITH_LOGGER
            lg->writefAll("Curvature N-S Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature N-S Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature N-S Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature N-S Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature N-S Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

            find_min_max(R, geotop::input::gDoubleNoValue, &max, &min);
#ifdef WITH_LOGGER
            lg->writefAll("Curvature W-E Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature W-E Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature W-E Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature W-E Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature W-E Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

            find_min_max(S, geotop::input::gDoubleNoValue, &max, &min);
#ifdef WITH_LOGGER
            lg->writefAll("Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature NW-SE Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature NW-SE Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature NW-SE Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER

            find_min_max(T, geotop::input::gDoubleNoValue, &max, &min);
#ifdef WITH_LOGGER
            lg->writefAll("Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
#else
#ifdef USE_DOUBLE_PRECISION_OUTPUT
            printf("Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
            fprintf(flog,"Curvature NE-SW Min:%12g  Max:%12g \n",min,max);
#else
            printf("Curvature NE-SW Min:%f  Max:%f \n",min,max);
            fprintf(flog,"Curvature NE-SW Min:%f  Max:%f \n",min,max);
#endif //USE_DOUBLE_PRECISION_OUTPUT
#endif //WITH_LOGGER
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
#ifdef WITH_LOGGER
    lg->writefAll("\nPOINTS:\n");
    lg->writefAll("ID,East[m],North[m],Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-],CurvatureN-S[1/m],CurvatureW-E[1/m],CurvatureNW-SE[1/m],CurvatureNE-SW[1/m],DepthFreeSurface[mm],Hor,maxSWE[mm],Lat[deg],Long[deg]\n");
#else
    fprintf(flog,"\nPOINTS:\n");
    fprintf(flog,"ID,East[m],North[m],Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],SkyViewFactor[-],CurvatureN-S[1/m],CurvatureW-E[1/m],CurvatureNW-SE[1/m],CurvatureNE-SW[1/m],DepthFreeSurface[mm],Hor,maxSWE[mm],Lat[deg],Long[deg]\n");
#endif

    for(r=1;r< par->chkpt.getRows();r++){
        for(c=1;c<ptTOT;c++){
#ifdef WITH_LOGGER
            lg->writefAll("%f",par->chkpt[r][c]);
#else
            fprintf(flog,"%f",par->chkpt[r][c]);
#endif
            if (c<ptTOT) {
#ifdef WITH_LOGGER
                lg->writefAll(",");
#else
                fprintf(flog, ",");
#endif
            }
        }
#ifdef WITH_LOGGER
        lg->writefAll("\n");
#else
        fprintf(flog,"\n");
#endif
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
        par->rc.resize(par->chkpt.getRows(),2+1);
        for(i=1;i< par->chkpt.getRows();i++){
            par->rc[i][1]=1;
            par->rc[i][2]=i;
        }
    }

    //6. SET PROPERTIES
    top->East.resize(1+1,par->chkpt.getRows());
    top->North.resize(1+1,par->chkpt.getRows());
    top->Z0.resize(1+1,par->chkpt.getRows());
    land->LC.resize(1+1,par->chkpt.getRows());
    land->delay.resize(1+1,par->chkpt.getRows());
    sl->type.resize(1+1,par->chkpt.getRows());
    top->slope.resize(1+1,par->chkpt.getRows());
    top->aspect.resize(1+1,par->chkpt.getRows());
    top->curvature1.resize(1+1,par->chkpt.getRows());
    top->curvature2.resize(1+1,par->chkpt.getRows());
    top->curvature3.resize(1+1,par->chkpt.getRows());
    top->curvature4.resize(1+1,par->chkpt.getRows());
    top->sky.resize(1+1,par->chkpt.getRows());
    top->pixel_type.resize(1+1,par->chkpt.getRows());
    top->BC_counter.resize(1+1,par->chkpt.getRows());
    top->BC_DepthFreeSurface.resize(par->chkpt.getRows());
    par->maxSWE.resize(2,par->chkpt.getRows());
    top->horizon_point.resize(1+1,par->chkpt.getRows());
    top->dzdE.resize(1+1,par->chkpt.getRows());
    top->dzdN.resize(1+1,par->chkpt.getRows());
    top->latitude.resize(2,par->chkpt.getRows());
    top->longitude.resize(2,par->chkpt.getRows());
    par->IDpoint.resize(par->chkpt.getRows());
	
	
    for(i=1;i< par->chkpt.getRows();i++){
        top->East[1][i]=par->chkpt[i][ptX];
        top->North[1][i]=par->chkpt[i][ptY];
        top->Z0[1][i]=par->chkpt[i][ptZ];
        land->LC[1][i]=par->chkpt[i][ptLC];

        if((long)land->LC[1][i] <= 0){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error:: Point %ld has land cover type <= 0. This is not admitted.\n",i);
            fclose(f);
#ifdef WITH_LOGGER
            lg->logsf(geotop::logger::ERROR,
            "Point %ld has land cover type <= 0. This is not admitted.",
            i);
            lg->log("Geotop failed. See failing report (15).",
            geotop::logger::CRITICAL);
            exit(1);
#else
            t_error("Fatal Error! Geotop is closed. See failing report (15).");
#endif
        }

        sl->type[1][i]=(long)par->chkpt[i][ptSY];

        if(sl->type[1][i] <= 0){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error:: Point %ld has soil type <= 0. This is not admitted.\n",i);
            fclose(f);
#ifdef WITH_LOGGER
            lg->logsf(geotop::logger::ERROR,
            "Point %ld has soil type <= 0. This is not admitted.",
            i);
            lg->log("Geotop failed. See failing report (16).",
            geotop::logger::CRITICAL);
            exit(1);
#else
            t_error("Fatal Error! Geotop is closed. See failing report (16).");
#endif
        }

        top->slope[1][i]=par->chkpt[i][ptS];
        top->aspect[1][i]=par->chkpt[i][ptA];
        top->sky[1][i]=par->chkpt[i][ptSKY];
        top->curvature1[1][i]=par->chkpt[i][ptCNS];
        top->curvature2[1][i]=par->chkpt[i][ptCWE];
        top->curvature3[1][i]=par->chkpt[i][ptCNwSe];
        top->curvature4[1][i]=par->chkpt[i][ptCNeSw];
        top->pixel_type[1][i]=1;
        top->BC_counter[1][i]=i;
        top->BC_DepthFreeSurface[i]=par->chkpt[i][ptDrDEPTH];
        top->horizon_point[1][i]=(long)par->chkpt[i][ptHOR];
        top->dzdE[1][i]=0.;
        top->dzdN[1][i]=0.;
        land->delay[1][i]=0.;

        if(sl->type[1][i] <= 0){
            f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
            fprintf(f, "Error:: Point %ld has horizon type <= 0. This is not admitted.\n",i);
            fclose(f);
#ifdef WITH_LOGGER
            lg->logsf(geotop::logger::ERROR,
            "Point %ld has horizon type <= 0. This is not admitted.",
            i);
            lg->log("Geotop failed. See failing report (17).",
            geotop::logger::CRITICAL);
            exit(1);
#else
            t_error("Fatal Error! Geotop is closed. See failing report (17).");
#endif
        }

        par->maxSWE[1][i]=par->chkpt[i][ptMAXSWE];
        top->latitude[1][i]=par->chkpt[i][ptLAT];
        top->longitude[1][i]=par->chkpt[i][ptLON];
        par->IDpoint[i]=(long)par->chkpt[i][ptID];

//        IT->bed[1][i]=par->chkpt[i][ptBED];
//        if( (long)IT->bed[1][i] == geotop::input::gDoubleNoValue ) IT->bed[1][i] = 1.E99;

    }

    //7. SET PAR these are all scalar from now on.. 
   // for (size_t i=1; i< par->init_date.size(); i++) {
        par->output_soil[1]=0.;
        par->output_snow[1]=0.;
        par->output_glac[1]=0.;
        par->output_surfenergy[1]=0.;
        par->output_vegetation[1]=0.;
        par->output_meteo[1]=0.;
    //}

    par->output_soil_bin = 0;
    par->output_snow_bin = 0;
    par->output_glac_bin = 0;
    par->output_surfenergy_bin = 0;
    par->output_meteo_bin = 0;

    //8. READ HORIZONS
    //find max top->horizon_point
    top->num_horizon_point=0;
    for(r=1;r<=top->horizon_point.getRows()-1;r++){
        for(c=1;c<=top->horizon_point.getCols()-1;c++){
            if (top->horizon_point[r][c] > top->num_horizon_point) top->num_horizon_point = top->horizon_point[r][c];
        }
    }
    top->horizon_height=(double ***)malloc(top->num_horizon_point*sizeof(double**));
    top->horizon_numlines=(long *)malloc(top->num_horizon_point*sizeof(long));
    for(i=1; i<=top->num_horizon_point; i++){

        c=0;
        do{
            flag = 0;
            if (c < par->chkpt.getRows()-1) {
                if (top->horizon_point[1][c+1] != i) c++;
            }
            if (c < par->chkpt.getRows()-1) {
                if (top->horizon_point[1][c+1] != i) flag=1;
            }

        }while (flag == 1 && c < par->chkpt.getRows()-1);

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


void set_bedrock(InitTools *IT, Soil *sl, Channel *cnet, Par *par, Topo *top, GeoMatrix<double>& LC, FILE *flog){

	GeoMatrix<double> B;
	GeoTensor<double> T;
	GeoVector<double> WT;
	long i, j, l, r, c, sy, synew;
	double zlim, z;
	short yes=0;
	FILE *f;

#ifdef WITH_LOGGER
    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();
#endif

	if (!mio::IOUtils::fileExists(string(geotop::common::Variables::files[fbed]) + string(ascii_esri))) {
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f, "Error:: File %s is missing. Please check if you have a bedrock topography map. If it is not available, remove the file name and keyword from input file\n",geotop::common::Variables::files[fbed+1].c_str());
		fclose(f);
#ifdef WITH_LOGGER
            lg->logsf(geotop::logger::ERROR,
            "File %s is missing. Please check if you have a bedrock topography map. If it is not available, remove the file name and keyword from input file.",
            geotop::common::Variables::files[fbed+1].c_str());
            lg->log("Geotop failed. See failing report (18).",
            geotop::logger::CRITICAL);
            exit(1);
#else
		t_error("Fatal Error! Geotop is closed. See failing report (18).");
#endif
	}

#ifdef WITH_LOGGER
    lg->logf("A bedrock depth map has been assigned and read from %s",
            geotop::common::Variables::files[fbed].c_str());
#else
	printf("A bedrock depth map has been assigned and read from %s\n\n",geotop::common::Variables::files[fbed].c_str());
	fprintf(flog,"A bedrock depth map has been assigned and read from %s\n\n",geotop::common::Variables::files[fbed].c_str());
#endif

	par->bedrock = 1;
	meteoio_readMap(string(geotop::common::Variables::files[fbed]), B);

	//check if bedrock depth is above soil lower border, otherwise we do not need to calculate anything
	z = 0.;
	for (l=1; l<geotop::common::Variables::Nl+1; l++) {
		z += sl->pa[1][jdz][l];
	}
	for (i=1; i<=par->total_pixel; i++){
		r = top->rc_cont[i][1];
		c = top->rc_cont[i][2];
		if(B[r][c] < z) yes=1;//if (IT->bed[r][c] < z) yes = 1;
	}

	if (yes == 1){

		//consistency check
		if (IT->init_water_table_depth.size() != sl->pa.getDh()){
			f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
			fprintf(f, "Error:: Error in bedrock calculations");
			fclose(f);
#ifdef WITH_LOGGER
            lg->logsf(geotop::logger::DEBUG,
                    "IT->init_water_table_depth.size()=%u sl->pa.getDh()=%u",
                    IT->init_water_table_depth.size(),
                    sl->pa.getDh());
            lg->log("Error in bedrock calculations.",
                    geotop::logger::ERROR);
			lg->log(" Geotop failed. See failing report (19).",
                    geotop::logger::CRITICAL);
            exit(1);
#else
			cout << "IT->init_water_table_depth.size()" << IT->init_water_table_depth.size() << " sl->pa.getDh()" << sl->pa.getDh() << endl;
			t_error("Fatal Error! Geotop is closed. See failing report (19).");
#endif
		}

		//	rewrite soil type
		T.resize(sl->pa.getDh(), nsoilprop+1, geotop::common::Variables::Nl+1);
		for(i=1;i<sl->pa.getDh();i++){
			for(j=1;j<=nsoilprop;j++){
				for(l=1;l<=geotop::common::Variables::Nl;l++){
					T(i,j,l) = sl->pa(i,j,l);
				}
			}
		}
		sl->pa.resize(par->total_pixel+par->total_channel+1, nsoilprop+1, geotop::common::Variables::Nl+1);

		//	rewrite initial water table depth
		WT.resize(IT->init_water_table_depth.size(),geotop::input::gDoubleNoValue);
		for(size_t k=1;k<IT->init_water_table_depth.size();k++) {
			WT.data.push_back(IT->init_water_table_depth[k]);
		}
		IT->init_water_table_depth.resize(par->total_pixel+par->total_channel+1);


		//assign jdz (is needed later)
		for(i=1;i<sl->pa.getDh();i++){
			for(l=1;l<=geotop::common::Variables::Nl;l++){
				sl->pa(i,jdz,l) = T(1,jdz,l);
			}
		}


		for (i=1; i<=par->total_pixel+par->total_channel; i++) {

			if (i<=par->total_pixel) {
				r = top->rc_cont[i][1];
				c = top->rc_cont[i][2];
				sy = sl->type[r][c];
				synew = i;
				sl->type[r][c] = synew;
				z = 0.0;
			}else {
				r = cnet->r[i-par->total_pixel];
				c = cnet->c[i-par->total_pixel];
				sy = cnet->soil_type[i-par->total_pixel];
				synew = i;
				cnet->soil_type[i-par->total_pixel] = synew;
				z = par->depr_channel;
			}

			IT->init_water_table_depth[synew] = WT[sy-1]; //WT->co[sy];

			//zlim = IT->bed->co[r][c];
			zlim = B[r][c];

			for(l=1;l<=geotop::common::Variables::Nl;l++){

				z += 0.5*sl->pa(synew,jdz,l);

				if(z <= zlim){

					for(j=1;j<=nsoilprop;j++){
						sl->pa(synew,j,l) = T(sy,j,l);
					}

				}else{

					for(j=1;j<=nsoilprop;j++){
						sl->pa(synew,j,l) = IT->pa_bed(sy,j,l);//IT->pa_bed->co[sy][j][l] ;
					}
				}

				z += 0.5*sl->pa(synew,jdz,l);

			}
		}

	}
}





/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


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

#ifdef WITH_LOGGER
//TODO: Add meaningful information about the files checked
short file_exists(short key){
    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();
#else
short file_exists(short key, FILE *flog){
#endif

    //no keyword -> -1
    //keyword but does not exist -> 0
    //keyword and exists -> 1

    if(geotop::common::Variables::files[key] == geotop::input::gStringNoValue){
#ifdef WITH_LOGGER
        lg->log("not present in file list");
#else
        printf("not present in file list\n");
        fprintf(flog,"not present in file list\n");
#endif
        return (-1);
    }else{
        bool is_present = mio::IOUtils::fileExists(string(geotop::common::Variables::files[key]) + string(ascii_esri));

        if (is_present) {
#ifdef WITH_LOGGER
            lg->logf("EXISTING in format %d", 3);
#else
            printf("EXISTING in format %d\n", 3);
            fprintf(flog, "EXISTING in format %d\n", 3);
#endif
            return (1);
        }else{
#ifdef WITH_LOGGER
            lg->log("File" + geotop::common::Variables::files[key] + " not existing", geotop::logger::WARNING);
#else
            printf("not existing\n");
            fprintf(flog, "not existing\n");
#endif
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
