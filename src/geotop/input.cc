
/** STATEMENT:

    Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
    Geotop 2.0.0 - 31 Oct 2013

    Copyright (c), 2013 - Stefano Endrizzi

    This file is part of Geotop 2.0.0

    Geotop 2.0.0  is a free software and is distributed under GNU General Public License 
    v. 3.0 <http://www.gnu.org/licenses/>
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
    PARTICULAR PURPOSE

    If you have satisfactorily used the code, please acknowledge the authors.

*/

#include "struct.geotop.h"
#include "input.h"
#include "parameters.h"
#include "geomorphology.0875.h"
#include "geomorphology.h"
#include "pedo.funct.h"
#include "networks.h"
#include "constants.h"
#include "dtm_resolution.h"
#include "rw_maps.h"
#include "extensions.h"
#include "tabs.h"
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
#include "logger.h"
#include <iostream>

extern long number_novalue, number_absent;
extern char *string_novalue;

extern T_INIT *UV;
extern const char *WORKING_DIRECTORY;
extern char **files, *logfile;
extern long Nl, Nr, Nc;
extern char *keywords_char[num_par_char];
extern char *SuccessfulRunFile, *FailedRunFile;
extern long i_sim0, i_run0, i_run;
extern double elapsed_time_start, cum_time, max_time;

//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************
//****************************************************************************************************

/** Subroutine which reads input data, performs  geomporphological analisys and allocates data */
void get_all_input(long argc, char *argv[], TOPO *top, SOIL *sl, LAND *land,
                   METEO *met, WATER *wat, CHANNEL *cnet,
                   PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times)

{
    GEOLOG_PREFIX(__func__);

    FILE *flog; /** log file */
    FILE *f; /** failed run file*/
    DOUBLEMATRIX *M;
    std::unique_ptr<INIT_TOOLS> IT;

    short a, success, added_JDfrom0=0, added_wind_xy=0, added_wind_dir=0,
            added_cloud=0, added_Tdew=0, added_RH=0, added_Pint=0, old=0, recovered=0;
    long l, r, c, i, ist, j, n, sy, num_cols, num_lines, day, month, year, hour,
            minute;
    double z, th_oversat, JD, k_snowred, maxSWE, SWE, D, cosslope, **matrix;
    char *temp, *name, **temp2;
    char rec[ ]= {"_recNNNN"},crec[ ]= {"_crecNNNN"};

    IT.reset(new INIT_TOOLS{});



    // logfile =  nullptr;  //join_strings(WORKING_DIRECTORY, logfile_name);
    flog = fopen(logfile, "w");


    // reads the parameters in __control_parameters
    temp = join_strings(WORKING_DIRECTORY, program_name);
    success = read_inpts_par(par, land, times, sl, met, IT.get(), temp);
    free(temp);

    // correct state pixel
    par->Tzrun = 0;
    par->wzrun = 0;
    par->Tzmaxrun = 0;
    par->Tzminrun = 0;
    par->wzmaxrun = 0;
    par->wzminrun = 0;
    par->dUzrun = 0;
    par->SWErun = 0;

    if (strcmp(files[fTrun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->Tzrun = 1;
    }
    if (strcmp(files[fwrun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->wzrun = 1;
    }
    if (strcmp(files[fTmaxrun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->Tzmaxrun = 1;

    }
    if (strcmp(files[fwmaxrun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->wzmaxrun = 1;
    }
    if (strcmp(files[fTminrun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->Tzminrun = 1;
    }
    if (strcmp(files[fwminrun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->wzminrun = 1;
    }
    if (strcmp(files[fdUrun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->dUzrun = 1;
    }
    if (strcmp(files[fSWErun], string_novalue) != 0)
    {
        if (par->point_sim == 1) par->state_pixel = 1;
        if (par->state_pixel == 1) par->SWErun = 1;
    }
    if (par->newperiodinit == 2 && (par->Tzrun == 0 || par->wzrun == 0))
    {
        f = fopen(FailedRunFile, "w");
        fprintf(f, "Error: You have to assign a name to the Tzrun and wzrun files\n");
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report.");
    }

    // soil parameters
    success = read_soil_parameters(files[fspar], IT.get(), sl, par->soil_type_bedr_default);
    Nl=sl->pa->nch;

    // pointlist files
    success = read_point_file(files[fpointlist], IT->point_col_names, par);

    // max time that the simulation is supposed to model
    max_time = 0.;
    for (i=1; i<=par->init_date->nh; i++)
    {
        max_time += (*par->run_times)(i)*(par->end_date->co[i] -
                                           par->init_date->co[i])*86400.;//seconds
    }

    // recovering
    par->delay_day_recover = 0.0;
    par->n_ContRecovery = 0;

    if (par->recover > 0)
    {
        if (par->saving_points->nh < par->recover)
        {
            f = fopen(FailedRunFile, "w");
            fprintf(f,
                    "Error: recover index higher than the length of the saving points vector");
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (1).");
        }
        par->delay_day_recover = par->saving_points->co[par->recover];
        recovered = 1;
    }

    // continuous recovering
    if (par->RunIfAnOldRunIsPresent != 1)
    {
        if (existing_file_wext(SuccessfulRunFile,"") == 1)
        {
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
        if (existing_file_wext(FailedRunFile,"") == 1)
        {
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

    if (par->ContRecovery > 0 && par->recover == 0)
    {

        if (existing_file_woext(files[rsux]) != 1)
        {
            temp = join_strings(files[rsux], ".old");
            if (existing_file_woext(temp) == 1) old = 1;
            free(temp);
        }

        if (old == 1)
        {
            name = join_strings(files[rtime], ".old");
        }
        else
        {
            name = assign_string(files[rtime]);
        }

        if (existing_file_wext(name, textfile) == 1)
        {
            temp = join_strings(name, textfile);
            matrix = read_txt_matrix_2(temp, 33, 44, 8, &num_lines);
            par->delay_day_recover = matrix[0][0]/secinday;
            i_run0 = (long)matrix[0][3];
            i_sim0 = (long)matrix[0][4];
            cum_time = matrix[0][5];
            elapsed_time_start = matrix[0][6];
            par->n_ContRecovery = (long)matrix[0][7]+1;

            for (i=0; i<num_lines; i++)
            {
                free(matrix[i]);
            }
            free(matrix);
            free(temp);
            recovered = 1;
        }
        free(name);
    }

    // time indices
    convert_JDfrom0_JDandYear(par->init_date->co[i_sim0], &JD, &year);
    convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);

    i_run = i_run0;//Run index

    /**************************************************************************************************/
    /*! Reading of the Input files:                                                                   */
    /**************************************************************************************************/

    if (par->point_sim!=1)  //distributed simulation
    {
        read_inputmaps(top, land, sl, par, IT.get(), flog);
    }
    else
    {
        read_optionsfile_point(par, top, land, sl, times, IT.get(), flog);
    }

    Nr=top->Z0->nrh;
    Nc=top->Z0->nch;

    par->total_pixel=0;
    par->total_area=0.;
    for (r=1; r<=Nr; r++)
    {
        for (c=1; c<=Nc; c++)
        {
            if ((long)land->LC->co[r][c]!=number_novalue)
            {
                par->total_pixel++;
                if (par->point_sim != 1)
                {
                    par->total_area += (UV->U->co[1]*UV->U->co[2])/cos(top->slope->co[r][c]*Pi/180.);
                }
                else
                {
                    par->total_area += (UV->U->co[1]*UV->U->co[2]);
                }
            }
        }
    }

    top->Z = find_Z_of_any_layer(top->Z0, top->slope, land->LC, sl, par->point_sim);

    top->Jdown = new_longmatrix(par->total_pixel, 4);
    top->Qdown = new_doublematrix(par->total_pixel, 4);
    for (i=1; i<=par->total_pixel; i++)
    {
        for (j=1; j<=4; j++)
        {
            top->Jdown->co[i][j] = i;
            top->Qdown->co[i][j] = 0.;
        }
    }

    geolog << "Valid pixels: " << par->total_pixel << std::endl;
    geolog << "Number of nodes: " << (Nl+1)*par->total_pixel << std::endl;
    geolog << "Novalue pixels: " << (Nr*Nc-par->total_pixel) << std::endl;
    geolog << "Basin area: " << (double)par->total_pixel*UV->U->co[1]*UV->U->co[2]/1.E6
           << " km2" << std::endl;


    /**************************************************************************************************/
    // Reading  RAIN data file, METEO data file and CLOUD data file

    num_cols = (long)nmet;

    // meteo data
    met->data=(double ***)malloc(met->st->E->nh*sizeof(double **));
    // number of line of meteo data
    met->numlines=(long *)malloc(met->st->E->nh*sizeof(long));
    // meteo variables for the current instant
    met->var=(double **)malloc(met->st->E->nh*sizeof(double *));
    // horizon for meteo stations
    met->horizon=(double ***)malloc(met->st->E->nh*sizeof(double **));
    // number of line in the horizon file
    met->horizonlines=(long *)malloc(met->st->E->nh*sizeof(long));
    // line of met->data used (stored in memory to avoid from searching from the first line)
    met->line_interp_WEB=(long *)malloc(met->st->E->nh*sizeof(long));
    met->line_interp_Bsnow=(long *)malloc(met->st->E->nh*sizeof(long));
    met->line_interp_WEB_LR=0;
    met->line_interp_Bsnow_LR=0;

    success = read_meteostations_file(met->imeteo_stations.get(), met->st.get(), files[fmetstlist],
                                      IT->meteostations_col_names);

    for (i=1; i<=met->st->E->nh; i++)
    {

        if ((*met->imeteo_stations)(1) != number_novalue)
        {
            ist = (*met->imeteo_stations)(i);
        }
        else
        {
            ist = i;
        }

        // initialize
        met->line_interp_WEB[i-1] = 0;
        met->line_interp_Bsnow[i-1] = 0;

        // allocate var
        met->var[i-1] = (double *)malloc(num_cols*sizeof(double));

        // read horizon
        met->horizon[i-1] = read_horizon(1, ist, files[fhormet],
                                         IT->horizon_col_names, &num_lines, flog);
        met->horizonlines[i-1] = num_lines;

        // filename
        if (strcmp(files[fmet], string_novalue) != 0)
        {

            // read matrix
            temp=namefile_i(files[fmet], ist);
            met->data[i-1] = read_txt_matrix(temp, 33, 44, IT->met_col_names, nmet, &num_lines);

            if ((long)met->data[i-1][0][iDate12] == number_absent
                && (long)met->data[i-1][0][iJDfrom0] == number_absent)
            {
                f = fopen(FailedRunFile, "w");
                fprintf(f, "Error:: Date Column missing in file %s\n",temp);
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report (2).");
            }
            met->numlines[i-1] = num_lines;

            /* fixing dates: converting times in the same standard time set for the simulation and fill
           JDfrom0 */
            added_JDfrom0 = fixing_dates(ist, met->data[i-1], par->ST, met->st->ST->co[i],
                                         met->numlines[i-1], iDate12, iJDfrom0);

            check_times(ist, met->data[i-1], met->numlines[i-1], iJDfrom0);

            // find clouds
            if (strcmp(IT->met_col_names[itauC], string_novalue) != 0)
            {
                if ((long)met->data[i-1][0][itauC] == number_absent || par->ric_cloud == 1)
                {
                    added_cloud = fill_meteo_data_with_cloudiness(met->data[i-1],
                                                                  met->numlines[i-1],
                                                                  met->horizon[i-1],
                                                                  met->horizonlines[i-1],
                                                                  met->st->lat->co[i],
                                                                  met->st->lon->co[i],
                                                                  par->ST,
                                                                  met->st->Z->co[i],
                                                                  met->st->sky->co[i],
                                                                  0.0, par->ndivdaycloud,
                                                                  par->dem_rotation,
                                                                  par->Lozone,
                                                                  par->alpha_iqbal,
                                                                  par->beta_iqbal,
                                                                  0.);
                }
            }

            // calculate Wx and Wy if Wspeed and direction are given
            if (par->wind_as_xy == 1)
            {
                added_wind_xy = fill_wind_xy(met->data[i-1], met->numlines[i-1], iWs, iWdir,
                                             iWsx, iWsy, IT->met_col_names[iWsx],
                                             IT->met_col_names[iWsy]);
            }

            // calcululate Wspeed and direction if Wx and Wy are given
            if (par->wind_as_dir == 1)
            {
                added_wind_dir = fill_wind_dir(met->data[i-1], met->numlines[i-1], iWs, iWdir,
                                               iWsx, iWsy, IT->met_col_names[iWs],
                                               IT->met_col_names[iWdir]);
            }

            // find Tdew
            if (par->vap_as_Td == 1)
            {
                added_Tdew = fill_Tdew(i, met->st->Z.get(), met->data[i-1], met->numlines[i-1], iRh,
                                       iT, iTdew, IT->met_col_names[iTdew], par->RHmin);
            }

            // find RH
            if (par->vap_as_RH == 1)
            {
                added_RH = fill_RH(i,  met->st->Z.get(), met->data[i-1], met->numlines[i-1], iRh,
                                   iT, iTdew, IT->met_col_names[iRh]);
            }

            // find Prec Intensity
            if ( (*par->linear_interpolation_meteo)(i) == 1
                 && (long)met->data[i-1][0][iPrec] != number_absent)
            {
                f = fopen(FailedRunFile, "w");
                fprintf(f, "Meteo data for station %ld contain precipitation as volume, \
but Linear Interpolation is set. This is not possible, the precipitation data are removed.\n", i);
                fprintf(f, "If you want to use precipitation as volume, you cannot set \
keyword LinearInterpolation at 1.\n");
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report (3).");
            }

            if (par->prec_as_intensity == 1)
            {
                added_Pint = fill_Pint(i, met->data[i-1], met->numlines[i-1], iPrec, iPrecInt,
                                       iJDfrom0, IT->met_col_names[iPrecInt]);
            }

            // rewrite completed files
            rewrite_meteo_files(met->data[i-1], met->numlines[i-1], IT->met_col_names,
                                temp, added_JDfrom0, added_wind_xy, added_wind_dir, added_cloud,
                                added_Tdew, added_RH, added_Pint);

            // calculate Wx and Wy if Wspeed and direction are given
            if (par->wind_as_xy != 1)
            {
                added_wind_xy = fill_wind_xy(met->data[i-1], met->numlines[i-1], iWs, iWdir,
                                             iWsx, iWsy, IT->met_col_names[iWsx],
                                             IT->met_col_names[iWsy]);
            }

            // find Prec Intensity
            if (par->prec_as_intensity != 1)
            {
                added_Pint = fill_Pint(i, met->data[i-1], met->numlines[i-1], iPrec, iPrecInt,
                                       iJDfrom0, IT->met_col_names[iPrecInt]);
            }

            // find Tdew
            if (par->vap_as_Td != 1)
            {
                added_Tdew = fill_Tdew(i, met->st->Z.get(), met->data[i-1], met->numlines[i-1], iRh,
                                       iT, iTdew, IT->met_col_names[iTdew], par->RHmin);
            }

            free(temp);
        }
        else
        {

            geolog << "Warning: File meteo not in the list, meteo data not read, \
used default values" << std::endl;

            met->data[i-1] = (double **)malloc(2*sizeof(double *));

            for (n=1; n<=2; n++)
            {
                met->data[i-1][n-1] = (double *)malloc(num_cols*sizeof(double));
                for (j=1; j<=nmet; j++)
                {
                    met->data[i-1][n-1][j-1] = (double)number_absent;
                }
            }

            met->data[i-1][0][iJDfrom0] = 0.;
            met->data[i-1][1][iJDfrom0] = 1.E10;

        }
    }


    // read LAPSE RATES FILE

    if (strcmp(files[fLRs], string_novalue) != 0)    //s stands for string
    {
        if (existing_file_wext(files[fLRs],
                               textfile)==0)
            printf("Lapse rate file unavailable. Check input files. If you do not have a lapse rate file,\
 remove its name and keywords from input file\n");
        temp = join_strings(files[fLRs], textfile);
        met->LRs = read_txt_matrix(temp, 33, 44, IT->lapserates_col_names, nlstot, &num_lines);
        free(temp);
        met->LRsnr = num_lines;
        par->LRflag=1;
        printf("\nLapse rate file read\n");
    }
    else
    {
        par->LRflag=0;
    }

    n = (long)nlstot;
    met->LRv = (double *)malloc(n*sizeof(double));
    met->LRd = (double *)malloc(n*sizeof(double));
    for ( i=0; i<nlstot; i++)
    {
        met->LRv[i] = (double)number_novalue;
        if (i == ilsTa)
        {
            met->LRd[i] = LapseRateTair;
        }
        else if (i == ilsTdew)
        {
            met->LRd[i] = LapseRateTdew;
        }
        else if (i == ilsPrec)
        {
            met->LRd[i] = LapseRatePrec;
        }
        else
        {
            met->LRd[i] = 0.0;
        }
    }


    // FIND A STATION WITH SHORTWAVE RADIATION DATA
    met->nstsrad=0;
    do
    {
        met->nstsrad++;
        a=0;
        if ( (long)met->data[met->nstsrad-1][0][iSW]!=number_absent
             || ((long)met->data[met->nstsrad-1][0][iSWb]!=number_absent
                 && (long)met->data[met->nstsrad-1][0][iSWd]!=number_absent ) )
            a=1;
    }
    while (met->nstsrad<met->st->Z->nh && a==0);
    if (a==0)
    {
        printf("WARNING: NO shortwave radiation measurements available\n");
        fprintf(flog,"WARNING: NO shortwave radiation measurements available\n");
    }
    else
    {
        printf("Shortwave radiation measurements from station %ld\n",met->nstsrad);
        fprintf(flog,"Shortwave radiation measurements from station %ld\n", met->nstsrad);
    }

    // FIND A STATION WITH CLOUD DATA
    met->nstcloud=0;
    do
    {
        met->nstcloud++;
        a=0;
        if ( (long)met->data[met->nstcloud-1][0][iC]!=number_absent
             || (long)met->data[met->nstcloud-1][0][itauC]!=number_absent )
            a=1;
    }
    while (met->nstcloud<met->st->Z->nh && a==0);
    if (a==0)
    {
        printf("WARNING: NO cloudiness measurements available\n");
        fprintf(flog,"WARNING: NO cloudiness measurements available\n");
    }
    else
    {
        printf("Cloudiness measurements from station %ld\n",met->nstcloud);
        fprintf(flog,"Cloudiness measurements from station %ld\n",met->nstcloud);
    }

    // FIND A STATION WITH LONGWAVE RADIATION DATA
    met->nstlrad=0;
    do
    {
        met->nstlrad++;
        a=0;
        if ( (long)met->data[met->nstlrad-1][0][iLWi]!=number_absent)
            a=1;
    }
    while (met->nstlrad<met->st->Z->nh && a==0);
    if (a==0)
    {
        printf("WARNING: NO longwave radiation measurements available\n");
        fprintf(flog,"WARNING: NO longwave radiation measurements available\n");
    }
    else
    {
        printf("Longwave radiation measurements from station %ld\n",met->nstlrad);
        fprintf(flog,"Longwave radiation measurements from station %ld\n",
                met->nstlrad);
    }

    // FIND A STATION WITH SURFACE TEMPERATURE ABOVE
    met->nstTs=0;
    do
    {
        met->nstTs++;
        a=0;
        if ( (long)met->data[met->nstTs-1][0][iTs]!=number_absent)
            a=1;
    }
    while (met->nstTs<met->st->Z->nh && a==0);
    if (a==0)
    {
        printf("WARNING: NO Surface temperature measurements available\n");
        fprintf(flog,"WARNING: NO Surface temperature measurements available\n");
    }
    else
    {
        printf("Surface temperature measurements from station %ld\n",met->nstTs);
        fprintf(flog,"Surface temperature  measurements from station %ld\n",
                met->nstTs);
    }

    //FIND A STATION WITH BOTTOM TEMPERATURE ABOVE
    met->nstTbottom=0;
    do
    {
        met->nstTbottom++;
        a=0;
        if ((long)met->data[met->nstTbottom-1][0][iTbottom]!=number_absent)
            a=1;
    }
    while (met->nstTbottom<met->st->Z->nh && a==0);
    if (a==0)
    {
        printf("WARNING: NO Bottom temperature measurements available\n");
        fprintf(flog,"WARNING: NO Bottom temperature measurements available\n");
    }
    else
    {
        printf("Bottom temperature measurements from station %ld\n", met->nstTbottom);
        fprintf(flog,"Bottom temperature measurements from station %ld\n", met->nstTbottom);
    }

    /**************************************************************************************************/
    // read INCOMING DISCHARGE

    met->qinv = (double *)malloc(2*sizeof(double));
    met->qinv[0] = 0.;
    met->qinv[1] = 0.;

    if (strcmp(files[fqin], string_novalue) != 0)
    {
        temp = join_strings(files[fqin], textfile);
        temp2 = (char **)malloc(2*sizeof(char *));
        temp2[0] = assign_string("Time");
        temp2[1] = assign_string("Qx");
        met->qins = read_txt_matrix(temp, 33, 44, temp2, 2, &num_lines);
        free(temp);
        free(temp2[0]);
        free(temp2[1]);
        free(temp2);
        met->qinsnr = num_lines;
        par->qin = 1;
        met->qinline = 0;
        printf("\nIncoming discharge file read\n");
    }
    else
    {
        par->qin = 0;
    }


    /**************************************************************************************************/
    /*! Completing the several time-indipendent input variables with the data of input-files:         */
    /**************************************************************************************************/
    /**************************************************************************************************/
    // Completing of "land" (of the type LAND):

    // Initialize matrix of shadow
    land->shadow=new_shortmatrix(Nr,Nc);
    initialize_shortmatrix(land->shadow, 0); /* initialized as if it was always NOT in shadow */

    // Check that there aren't cell with an undefined land use value
    z = 0.;
    l = 0;
    do
    {
        l++;
        z += sl->pa->co[1][jdz][l];
    }
    while (l<Nl && z < z_transp);
    land->root_fraction=new_doublematrix(par->n_landuses, l);
    initialize_doublematrix(land->root_fraction, 0.0);

    // check vegetation variable consistency
    if (jHveg!=jdHveg+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jz0thresveg!=jdz0thresveg+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jz0thresveg2!=jdz0thresveg2+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jLSAI!=jdLSAI+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jcf!=jdcf+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jdecay0!=jddecay0+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jexpveg!=jdexpveg+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jroot!=jdroot+jHveg-1)
        t_error("Vegetation variables not consistent");
    if (jrs!=jdrs+jHveg-1)
        t_error("Vegetation variables not consistent");

    // variables used to assign vegetation properties that change with time
    num_cols = jdvegprop + 1;
    land->vegpars=(double ***)malloc(par->n_landuses*sizeof(double **));
    land->vegparv=(double **)malloc(par->n_landuses*sizeof(double *));
    land->NumlinesVegTimeDepData=(long *)malloc(par->n_landuses*sizeof(long));

    land->vegpar.reset(new Vector<double>{jdvegprop});

    par->vegflag.reset(new Vector<short>{par->n_landuses});

    // time dependent vegetation parameters
    for (i=1; i<=par->n_landuses; i++)
    {
        if (strcmp(files[fvegpar], string_novalue) != 0) // s stands for string
        {
            temp = namefile_i_we2(files[fvegpar], i);

            if (existing_file_wext(temp, textfile)==1)
            {
                printf("There is a specific vegetation parameter file for land cover type = %ld\n", i);
                free(temp);
                temp = namefile_i(files[fvegpar], i);
                land->vegpars[i-1] = read_txt_matrix_2(temp, 33, 44, num_cols, &num_lines);
                free(temp);
                land->NumlinesVegTimeDepData[i-1] = num_lines;
                par->vegflag->co[i]=1;
            }
            else
            {
                free(temp);
                printf("There is NOT a specific vegetation parameter file for land cover type = %ld\n",
                       i);
            }
        }

        land->vegparv[i-1]=(double *)malloc(num_cols*sizeof(double));
        for (j=0; j<num_cols; j++)
        {
            land->vegparv[i-1][j] = (double)number_novalue;
        }

        // z0 (convert in m)
        land->ty->co[i][jz0]*=0.001;

        // find root fraction
        root(land->root_fraction->nch, land->ty->co[i][jroot], 0.0,
             sl->pa->co[1][jdz], land->root_fraction->co[i]);

        // error messages
        for (l=1; l<=met->st->Z->nh; l++)
        {
            if (0.001*land->ty->co[i][jHveg]>met->st->Vheight->co[l]
                || 0.001*land->ty->co[i][jHveg]>met->st->Theight->co[l])
            {
                f = fopen(FailedRunFile, "w");
                fprintf(f, "hc:%f m, zmu:%f m, zmt:%f m - hc must be lower than measurement height - \
land cover %ld, meteo station %ld\n",
                        0.001*land->ty->co[i][jHveg], met->st->Vheight->co[l],
                        met->st->Theight->co[l], i, l);
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report (5).");
            }
        }
    }

    //******************************************
    // file with time steps

    if (strcmp(files[ftsteps], string_novalue) != 0)
    {

        temp=join_strings(files[ftsteps],textfile);
        times->Dt_matrix = read_txt_matrix_2(temp, 33, 44, max_cols_time_steps_file+1, &num_lines);
        free(temp);
        times->numlinesDt_matrix = num_lines;
        par->tsteps_from_file=1;

    }
    else
    {
        par->tsteps_from_file=0;
    }

    /**************************************************************************************************/
    /*! Filling up of the struct "channel" (of the type CHANNEL):                                     */

    /* The number of channel-pixel are counted: */
    i=0;
    for (r=1; r<=Nr; r++)
    {
        for (c=1; c<=Nc; c++)
        {
            if (top->pixel_type->co[r][c]>=10) i++;
        }
    }
    fprintf(flog,"Channel pixels: %ld\n",i);
    par->total_channel = i;

    // allocate channel vectors/matrixes
    if (i==0)
        i=1;

    cnet->Vout = 0.;

    cnet->r.reset(new Vector<long>{i});
    cnet->c.reset(new Vector<long>{i});

    cnet->ch=new_longmatrix(Nr,Nc);
    initialize_longmatrix(cnet->ch, 0);

    cnet->ch_down.reset(new Vector<long>{i});

    cnet->length.reset(new Vector<double>{i});

    cnet->Vsup.reset(new Vector<double>{i});

    cnet->Vsub.reset(new Vector<double>{i});

    cnet->h_sup.reset(new Vector<double>{i});

    cnet->soil_type.reset(new Vector<long>{i});
    (*cnet->soil_type) = par->soil_type_chan_default;

    if (par->total_channel>1)
        enumerate_channels(cnet, land->LC, top->pixel_type, top->Z0, top->slope, number_novalue);

    cnet->ch3 = (long **)malloc((Nl+1)*sizeof(long *));
    for (l=0; l<=Nl; l++)
    {
        cnet->ch3[l] = (long *)malloc((i+1)*sizeof(long));
    }

    cnet->lch = new_longmatrix( (Nl+1)*i, 2);
    initialize_longmatrix(cnet->lch, 0);

    lch3_cont(cnet->ch3, cnet->lch, Nl, par->total_channel);


    /**************************************************************************************************/
    // Cont for Richards 3D
    n = Fminlong((*par->Nl_spinup)(i_sim0),Nl);

    // 3D
    top->i_cont=(long ***)malloc((n+1)*sizeof(long **));
    for (l=0; l<=n; l++)
    {
        top->i_cont[l]=(long **)malloc((Nr+1)*sizeof(long *));
        for (r=1; r<=Nr; r++)
        {
            top->i_cont[l][r]=(long *)malloc((Nc+1)*sizeof(long));
        }
    }

    top->lrc_cont=new_longmatrix( (n+1)*par->total_pixel, 3);
    initialize_longmatrix(top->lrc_cont, 0);

    i_lrc_cont(land->LC, top->i_cont, top->lrc_cont, n, Nr, Nc);

    // 2D
    top->j_cont=(long **)malloc((Nr+1)*sizeof(long *));
    for (r=1; r<=Nr; r++)
    {
        top->j_cont[r]=(long *)malloc((Nc+1)*sizeof(long));
        for (c=1; c<=Nc; c++)
        {
            top->j_cont[r][c] = 0;
        }
    }

    top->rc_cont=new_longmatrix(par->total_pixel, 2);
    initialize_longmatrix(top->rc_cont, 0);

    j_rc_cont(land->LC, top->j_cont, top->rc_cont, Nr, Nc);

    // plotted points
    if (par->state_pixel == 1)
    {
        par->jplot.reset(new Vector<long>{par->total_pixel});
        for (i=1; i<=par->total_pixel; i++)
        {
            for (j=1; j<=par->rc->nrh; j++)
            {
                if (top->rc_cont->co[i][1] == par->rc->co[j][1]
                    && top->rc_cont->co[i][2] == par->rc->co[j][2])
                {
                    (*par->jplot)(i) = j;
                }
            }
        }
    }

    // BEDROCK (adjusting soil properties)
    set_bedrock(IT.get(), sl, cnet, par, top, land->LC, flog);

    /**************************************************************************************************/
    /*! Completing of the initialization of SOIL structure                                            */
    /**************************************************************************************************/

    sl->SS = new SOIL_STATE {par->total_pixel, Nl};

    sl->VS.reset(new STATE_VEG{});
    initialize_veg_state(sl->VS.get(), par->total_pixel);

    sl->th=new_doublematrix(Nl,par->total_pixel);
    initialize_doublematrix(sl->th,(double)number_novalue);

    sl->Ptot=new_doublematrix(Nl,par->total_pixel);
    initialize_doublematrix(sl->Ptot,(double)number_novalue);

    sl->ET=new_doubletensor(Nl,Nr,Nc);
    initialize_doubletensor(sl->ET,0.);

    if (par->output_soil_bin == 1)
    {
        if (strcmp(files[fTav], string_novalue) != 0
            || strcmp(files[fTavsup], string_novalue) != 0)
        {
            sl->T_av_tensor=new_doublematrix(Nl,par->total_pixel);
            initialize_doublematrix(sl->T_av_tensor,0.);
        }

        if (strcmp(files[ficeav], string_novalue) != 0)
        {
            sl->thi_av_tensor=new_doublematrix(Nl,par->total_pixel);
            initialize_doublematrix(sl->thi_av_tensor,0.);
        }

        if (strcmp(files[fliqav], string_novalue) != 0)
        {
            sl->thw_av_tensor=new_doublematrix(Nl,par->total_pixel);
            initialize_doublematrix(sl->thw_av_tensor,0.);
        }

        if (strcmp(files[fpnet], string_novalue) != 0)
        {
            sl->Pnetcum.reset(new Vector<double>{par->total_pixel});
        }

        if (strcmp(files[fevap], string_novalue) != 0)
        {
            sl->ETcum.reset(new Vector<double>{par->total_pixel});
        }

    }

    if (existing_file(files[fwt0]) == 0)
    {

        for (i=1; i<=par->total_pixel; i++)
        {

            r = top->rc_cont->co[i][1];
            c = top->rc_cont->co[i][2];

            sy=sl->type->co[r][c];

            if ((long)IT->init_water_table_depth->co[sy] != number_novalue)
            {
                z = 0.;
                sl->SS->P->co[0][i] = -IT->init_water_table_depth->co[sy] *
                                      cos(top->slope->co[r][c]*Pi/180.);
                for (l=1; l<=Nl; l++)
                {
                    z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
                    sl->SS->P->co[l][i] = sl->SS->P->co[0][i] + z;
                    z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
                }
            }
            else
            {
                for (l=1; l<=Nl; l++)
                {
                    sl->SS->P->co[l][i] = sl->pa->co[sy][jpsi][l];
                }
            }
        }
    }
    else
    {

        M = read_map(2, files[fwt0], land->LC, UV, (double)number_novalue);

        for (i=1; i<=par->total_pixel; i++)
        {
            r = top->rc_cont->co[i][1];
            c = top->rc_cont->co[i][2];

            sy=sl->type->co[r][c];

            z = 0.;
            sl->SS->P->co[0][i] = -M->co[r][c]*cos(top->slope->co[r][c]*Pi/180.);
            for (l=1; l<=Nl; l++)
            {
                z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
                sl->SS->P->co[l][i] = sl->SS->P->co[0][i] + z;
                z += 0.5*sl->pa->co[sy][jdz][l]*cos(top->slope->co[r][c]*Pi/180.);
            }
        }
        free_doublematrix(M);
    }

    for (i=1; i<=par->total_pixel; i++)
    {
        r = top->rc_cont->co[i][1];
        c = top->rc_cont->co[i][2];

        sy=sl->type->co[r][c];

        for (l=1; l<=Nl; l++)
        {
            sl->SS->T->co[l][i]=sl->pa->co[sy][jT][l];

            sl->Ptot->co[l][i] = sl->SS->P->co[l][i];
            sl->th->co[l][i] = teta_psi(sl->SS->P->co[l][i], 0.0, sl->pa->co[sy][jsat][l],
                                        sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
                                        sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin,
                                        sl->pa->co[sy][jss][l]);

            th_oversat = Fmax( sl->SS->P->co[l][i], 0.0 ) * sl->pa->co[sy][jss][l];
            sl->th->co[l][i] -= th_oversat;

            if (sl->SS->T->co[l][i]<=Tfreezing)
            {

                // Theta_ice = Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
                sl->SS->thi->co[l][i] = sl->th->co[l][i] - teta_psi(Psif(sl->SS->T->co[l][i]),
                                                                    0.0,
                                                                    sl->pa->co[sy][jsat][l],
                                                                    sl->pa->co[sy][jres][l],
                                                                    sl->pa->co[sy][ja][l],
                                                                    sl->pa->co[sy][jns][l],
                                                                    1.-1./sl->pa->co[sy][jns][l],
                                                                    PsiMin,
                                                                    sl->pa->co[sy][jss][l]);

                // if Theta(without freezing) < Theta_unfrozen(in equilibrium with T):
                // Theta_ice is set at 0
                if (sl->SS->thi->co[l][i]<0)
                    sl->SS->thi->co[l][i]=0.0;

                // Psi is updated taking into account the freezing
                sl->th->co[l][i] -= sl->SS->thi->co[l][i];

                sl->SS->P->co[l][i] = psi_teta(sl->th->co[l][i] + th_oversat,
                                               sl->SS->thi->co[l][i], sl->pa->co[sy][jsat][l],
                                               sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
                                               sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l],
                                               PsiMin, sl->pa->co[sy][jss][l]);
            }
        }
    }


    if (par->state_pixel == 1)
    {
        if (strcmp(files[fTz], string_novalue) != 0
            || strcmp(files[fTzwriteend],
                      string_novalue) != 0) sl->Tzplot = new_doublematrix(par->rc->nrh, Nl);
        if (strcmp(files[fTzav], string_novalue) != 0
            || strcmp(files[fTzavwriteend],
                      string_novalue) != 0) sl->Tzavplot = new_doublematrix(par->rc->nrh, Nl);
        if (strcmp(files[fpsiztot], string_novalue) != 0
            || strcmp(files[fpsiztotwriteend],
                      string_novalue) != 0) sl->Ptotzplot = new_doublematrix(par->rc->nrh, Nl);
        if (strcmp(files[fpsiz], string_novalue) != 0
            || strcmp(files[fpsizwriteend],
                      string_novalue) != 0) sl->Pzplot = new_doublematrix_0(par->rc->nrh, Nl);
        if (strcmp(files[fliqz], string_novalue) != 0
            || strcmp(files[fliqzwriteend],
                      string_novalue) != 0) sl->thzplot = new_doublematrix(par->rc->nrh, Nl);
        if (strcmp(files[fliqzav], string_novalue) != 0
            || strcmp(files[fliqzavwriteend],
                      string_novalue) != 0) sl->thzavplot = new_doublematrix(par->rc->nrh, Nl);
        if (strcmp(files[ficez], string_novalue) != 0
            || strcmp(files[ficezwriteend],
                      string_novalue) != 0) sl->thizplot = new_doublematrix(par->rc->nrh, Nl);
        if (strcmp(files[ficezav], string_novalue) != 0
            || strcmp(files[ficezavwriteend],
                      string_novalue) != 0) sl->thizavplot = new_doublematrix(par->rc->nrh, Nl);
        if (strcmp(files[fsatz],
                   string_novalue) != 0) sl->satratio = new_doublematrix(par->rc->nrh, Nl);

        for (i=1; i<=par->rc->nrh; i++)
        {
            r = top->rc_cont->co[i][1];
            c = top->rc_cont->co[i][2];
            j = top->j_cont[r][c];
            sy = sl->type->co[r][c];
            for (l=1; l<=Nl; l++)
            {
                if (strcmp(files[fTz], string_novalue) != 0
                    || strcmp(files[fTzwriteend],
                              string_novalue) != 0) sl->Tzplot->co[i][l] = sl->SS->T->co[l][j];
                if (strcmp(files[fTzav], string_novalue) != 0
                    || strcmp(files[fTzavwriteend],
                              string_novalue) != 0) sl->Tzavplot->co[i][l] = sl->SS->T->co[l][j];
                if (strcmp(files[fliqz], string_novalue) != 0
                    || strcmp(files[fliqzwriteend],
                              string_novalue) != 0) sl->thzplot->co[i][l] = sl->th->co[l][j];
                if (strcmp(files[fliqzav], string_novalue) != 0
                    || strcmp(files[fliqzavwriteend],
                              string_novalue) != 0) sl->thzavplot->co[i][l] = sl->th->co[l][j];
                if (strcmp(files[ficez], string_novalue) != 0
                    || strcmp(files[ficezwriteend],
                              string_novalue) != 0) sl->thizplot->co[i][l] = sl->SS->thi->co[l][j];
                if (strcmp(files[ficezav], string_novalue) != 0
                    || strcmp(files[ficezavwriteend],
                              string_novalue) != 0) sl->thizavplot->co[i][l] = sl->SS->thi->co[l][j];
                if (strcmp(files[fpsiztot], string_novalue) != 0
                    || strcmp(files[fpsiztotwriteend],
                              string_novalue) != 0) sl->Ptotzplot->co[i][l] = sl->Ptot->co[l][j];
                if (strcmp(files[fsatz],
                           string_novalue) != 0) sl->satratio->co[i][l] = (sl->SS->thi->co[l][j] +
                                                                           sl->th->co[l][j] -
                                                                           sl->pa->co[sy][jres][l]) /
                                                                          (sl->pa->co[sy][jsat][l] - sl->pa->co[sy][jres][l]);
            }
            for (l=0; l<=Nl; l++)
            {
                if (strcmp(files[fpsiz], string_novalue) != 0
                    || strcmp(files[fpsizwriteend],
                              string_novalue) != 0) sl->Pzplot->co[i][l] = sl->SS->P->co[l][j];
            }
        }
    }

    if (recovered > 0)
    {
        assign_recovered_tensor_vector(old, par->recover, files[riceg], sl->SS->thi,
                                       top->rc_cont, par, land->LC);
        assign_recovered_tensor_vector(old, par->recover, files[rTg], sl->SS->T,
                                       top->rc_cont, par, land->LC);
        assign_recovered_tensor_vector(old, par->recover, files[rpsi], sl->SS->P,
                                       top->rc_cont, par, land->LC);

        assign_recovered_map_vector(old, par->recover, files[rwcrn], sl->VS->wrain.get(),
                                    top->rc_cont, par, land->LC);
        assign_recovered_map_vector(old, par->recover, files[rwcsn], sl->VS->wsnow.get(),
                                    top->rc_cont, par, land->LC);
        assign_recovered_map_vector(old, par->recover, files[rTv], sl->VS->Tv.get(),
                                    top->rc_cont, par, land->LC);
    }


    // channel soil
    cnet->SS = new SOIL_STATE {cnet->r->nh, Nl};

    cnet->th = new_doublematrix(Nl, cnet->r->nh);

    cnet->ET = new_doublematrix(Nl, cnet->r->nh);
    initialize_doublematrix(cnet->ET, 0.0);

    cnet->Kbottom.reset(new Vector<double>{cnet->r->nh});

    for (j=1; j<=par->total_channel; j++)
    {
        sy=(*cnet->soil_type)(j);
        r=(*cnet->r)(j);
        c=(*cnet->c)(j);

        cnet->SS->P->co[0][j] = sl->SS->P->co[0][top->j_cont[r][c]] +
                                par->depr_channel;

        for (l=1; l<=Nl; l++)
        {
            cnet->SS->P->co[l][j] = sl->Ptot->co[l][top->j_cont[r][c]] +
                                    par->depr_channel;
        }

        for (l=1; l<=Nl; l++)
        {
            cnet->SS->T->co[l][j]=sl->pa->co[sy][jT][l];

            cnet->th->co[l][j] = teta_psi(cnet->SS->P->co[l][j], 0.0,
                                          sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l],
                                          sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l],
                                          1.-1./sl->pa->co[sy][jns][l],
                                          PsiMin, sl->pa->co[sy][jss][l]);

            th_oversat = Fmax( cnet->SS->P->co[l][j], 0.0 ) * sl->pa->co[sy][jss][l];
            cnet->th->co[l][j] -= th_oversat;

            if (cnet->SS->T->co[l][j]<=Tfreezing)
            {
                // Theta_ice = Theta(without freezing) - Theta_unfrozen(in equilibrium with T)
                cnet->SS->thi->co[l][j] = cnet->th->co[l][j] - teta_psi(Psif( cnet->SS->T->co[l][j]),
                                                                        0.0,
                                                                        sl->pa->co[sy][jsat][l],
                                                                        sl->pa->co[sy][jres][l],
                                                                        sl->pa->co[sy][ja][l],
                                                                        sl->pa->co[sy][jns][l],
                                                                        1.-1./sl->pa->co[sy][jns][l],
                                                                        PsiMin,
                                                                        sl->pa->co[sy][jss][l]);

                // if Theta(without freezing)< Theta_unfrozen(in equilibrium with T)
                // Theta_ice is set at 0
                if (cnet->SS->thi->co[l][j]<0) cnet->SS->thi->co[l][j]=0.0;

                // Psi is updated taking into account the freezing
                cnet->th->co[l][j] -= cnet->SS->thi->co[l][j];
                cnet->SS->P->co[l][j] = psi_teta(cnet->th->co[l][j] + th_oversat,
                                                 cnet->SS->thi->co[l][j],
                                                 sl->pa->co[sy][jsat][l],
                                                 sl->pa->co[sy][jres][l],
                                                 sl->pa->co[sy][ja][l],
                                                 sl->pa->co[sy][jns][l],
                                                 1.-1./sl->pa->co[sy][jns][l],
                                                 PsiMin,
                                                 sl->pa->co[sy][jss][l]);
            }
        }
    }

    if (recovered > 0 && par->total_channel > 0)
    {
        assign_recovered_tensor_channel(old, par->recover, files[rpsich], cnet->SS->P,
                                        cnet->r.get(), cnet->c.get(), top->Z0);
        assign_recovered_tensor_channel(old, par->recover, files[ricegch],
                                        cnet->SS->thi, cnet->r.get(), cnet->c.get(), top->Z0);
        assign_recovered_tensor_channel(old, par->recover, files[rTgch], cnet->SS->T,
                                        cnet->r.get(), cnet->c.get(), top->Z0);

        for (i=1; i<=par->total_channel; i++)
        {
            for (l=1; l<=Nl; l++)
            {
                sy = (*cnet->soil_type)(i);
                cnet->th->co[l][i] = teta_psi(Fmin(cnet->SS->P->co[l][i],
                                                   psi_saturation(cnet->SS->thi->co[l][i],
                                                                  sl->pa->co[sy][jsat][l],
                                                                  sl->pa->co[sy][jres][l],
                                                                  sl->pa->co[sy][ja][l],
                                                                  sl->pa->co[sy][jns][l],
                                                                  1.-1./sl->pa->co[sy][jns][l])),
                                              cnet->SS->thi->co[l][i],
                                              sl->pa->co[sy][jsat][l],
                                              sl->pa->co[sy][jres][l],
                                              sl->pa->co[sy][ja][l],
                                              sl->pa->co[sy][jns][l],
                                              1.-1./sl->pa->co[sy][jns][l],
                                              PsiMin,
                                              sl->pa->co[sy][jss][l]);
            }
        }
    }

    // z boundary condition
    for (l=1; l<=Nl; l++)
    {
        par->Zboundary -= sl->pa->co[1][jdz][l];
    }

    if (par->Zboundary < 0)
    {
        f = fopen(FailedRunFile, "w");
        fprintf(f, "Z at which 0 annual temperature takes place is not lower than the soil column\n");
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report (6).");
    }

    par->Zboundary *= 1.E-3;  // convert in [m]

    /**************************************************************************************************/
    /*! Initialization of the struct "egy" (of the type ENERGY):*/

    if (par->output_surfenergy_bin == 1)
    {
        if (strcmp(files[fradnet], string_novalue) != 0)
        {
            egy->Rn_mean.reset(new Vector<double>{par->total_pixel});
            egy->Rn.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fradLWin], string_novalue) != 0)
        {
            egy->LWin_mean.reset(new Vector<double>{par->total_pixel});
            egy->LWin.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fradLW], string_novalue) != 0)
        {
            egy->LW_mean.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fradLW], string_novalue) != 0
            || strcmp(files[fradnet], string_novalue) != 0)
        {
            egy->LW.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fradSW], string_novalue) != 0)
        {
            egy->SW_mean.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fradSW], string_novalue) != 0
            || strcmp(files[fradnet], string_novalue) != 0)
        {
            egy->SW.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fLE], string_novalue) != 0)
        {
            egy->ET_mean.reset(new Vector<double>{par->total_pixel});
            egy->LE.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fH], string_novalue) != 0)
        {
            egy->H_mean.reset(new Vector<double>{par->total_pixel});
            egy->H.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fG], string_novalue) != 0)
        {
            egy->SEB_mean.reset(new Vector<double>{par->total_pixel});
            egy->G.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fTs], string_novalue) != 0)
        {
            egy->Ts_mean.reset(new Vector<double>{par->total_pixel});
            egy->Ts.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fradSWin], string_novalue) != 0)
        {
            egy->Rswdown_mean.reset(new Vector<double>{par->total_pixel});
            egy->SWin.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fradSWinbeam], string_novalue) != 0)
        {
            egy->Rswbeam_mean.reset(new Vector<double>{par->total_pixel});
            egy->SWinb.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fshadow], string_novalue) != 0)
        {
            egy->nDt_shadow.reset(new Vector<long>{par->total_pixel});
            egy->nDt_sun.reset(new Vector<long>{par->total_pixel});
            (*egy->nDt_shadow) = 0;
            (*egy->nDt_sun) = 0;
            egy->shad.reset(new Vector<short>{par->total_pixel});
        }
    }

    egy->sun = (double *)malloc(12*sizeof(double));

    if (times->JD_plots->nh > 1)
    {
        if (strcmp(files[pH], string_novalue) != 0
            || strcmp(files[pHg], string_novalue) != 0
            || strcmp(files[pG], string_novalue) != 0)
        {
            egy->Hgplot.reset(new Vector<double>{par->total_pixel});
            egy->Hgp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pH], string_novalue) != 0
            || strcmp(files[pHv], string_novalue) != 0)
        {
            egy->Hvplot.reset(new Vector<double>{par->total_pixel});
            egy->Hvp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pLE], string_novalue) != 0
            || strcmp(files[pLEg], string_novalue) != 0
            || strcmp(files[pG], string_novalue) != 0)
        {
            egy->LEgplot.reset(new Vector<double>{par->total_pixel});
            egy->LEgp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pLE], string_novalue) != 0
            || strcmp(files[pLEv], string_novalue) != 0)
        {
            egy->LEvplot.reset(new Vector<double>{par->total_pixel});
            egy->LEvp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pSWin], string_novalue) != 0)
        {
            egy->SWinplot.reset(new Vector<double>{par->total_pixel});
            egy->SWinp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pSWg], string_novalue) != 0
            || strcmp(files[pG], string_novalue) != 0)
        {
            egy->SWgplot.reset(new Vector<double>{par->total_pixel});
            egy->SWgp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pSWv], string_novalue) != 0)
        {
            egy->SWvplot.reset(new Vector<double>{par->total_pixel});
            egy->SWvp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pLWin], string_novalue) != 0)
        {
            egy->LWinplot.reset(new Vector<double>{par->total_pixel});
            egy->LWinp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pLWg], string_novalue) != 0
            || strcmp(files[pG], string_novalue) != 0)
        {
            egy->LWgplot.reset(new Vector<double>{par->total_pixel});
            egy->LWgp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pLWv], string_novalue) != 0)
        {
            egy->LWvplot.reset(new Vector<double>{par->total_pixel});
            egy->LWvp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pTs], string_novalue) != 0)
        {
            egy->Tsplot.reset(new Vector<double>{par->total_pixel});
            egy->Tsp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pTg], string_novalue) != 0)
        {
            egy->Tgplot.reset(new Vector<double>{par->total_pixel});
            egy->Tgp.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pTv], string_novalue) != 0)
        {
            egy->Tvplot.reset(new Vector<double>{par->total_pixel});
        }
    }


    // vectors used in energy_balance()
    egy->Tgskin_surr = new_doublematrix(Nr, Nc);
    initialize_doublematrix(egy->Tgskin_surr, 0.);

    egy->SWrefl_surr = new_doublematrix(Nr, Nc);
    initialize_doublematrix(egy->SWrefl_surr, 0.);

    egy->Dlayer.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers });
    egy->liq.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers });
    egy->ice.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers });
    egy->Temp.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers,0} );
    egy->deltaw.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers });

    egy->SWlayer.reset(new Vector<double>{ par->max_snow_layers + 1,0} );

    egy->soil_transp_layer.reset(new Vector<double>{land->root_fraction->nch});

    egy->dFenergy.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers, 0} );
    egy->udFenergy.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers - 1, 0});
    egy->Kth0.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers - 1, 0});
    egy->Kth1.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers - 1, 0});
    egy->Fenergy.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers, 0});
    egy->Newton_dir.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers, 0});
    egy->T0.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers, 0} );
    egy->T1.reset(new Vector<double>{ Nl + par->max_snow_layers + par->max_glac_layers, 0} );
    egy->Tstar.reset(new Vector<double>{Nl}); // soil temperature at which freezing begins
    egy->THETA.reset(new Vector<double>{Nl});  // water content (updated in the iterations)

    // allocate vector of soil layer contributions to evaporation (up to z_evap)
    z = 0.;
    l = 0;
    do
    {
        l++;
        z += sl->pa->co[1][jdz][l];
    }
    while (l<Nl && z < z_evap);
    egy->soil_evap_layer_bare.reset(new Vector<double> {l});
    egy->soil_evap_layer_veg.reset(new Vector<double> {l});

    printf("Soil water evaporates from the first %ld layers\n", egy->soil_evap_layer_bare->nh);
    printf("Soil water transpires from the first %ld layers\n", egy->soil_transp_layer->nh);
    fprintf(flog,"Soil water evaporates from the first %ld layers\n", egy->soil_evap_layer_bare->nh);
    fprintf(flog,"Soil water transpires from the first %ld layers\n", egy->soil_transp_layer->nh);

    /**************************************************************************************************/
    /*! Completing of the struct "water" (of the type WATER) */

    wat->Voutlandsub = 0.;
    wat->Voutlandsup = 0.;
    wat->Voutbottom = 0.;

    /* Initialization of wat->Pnet (liquid precipitation that reaches the sl surface in mm): */
    wat->Pnet=new_doublematrix(Nr,Nc);
    initialize_doublematrix(wat->Pnet,0.0);

    /* Initialization of wat->PrecTot (total precipitation (rain+snow) precipitation): */
    wat->PrecTot=new_doublematrix(Nr,Nc);
    initialize_doublematrix(wat->PrecTot,par->IPrec_default);

    /* Initialization of the matrices with the output of total precipitation and interception: */
    if (par->output_meteo_bin == 1 && strcmp(files[fprec], string_novalue) != 0)
    {
        wat->PrTOT_mean.reset(new Vector<double>{par->total_pixel});
        wat->PrSNW_mean.reset(new Vector<double>{par->total_pixel});
        wat->Pt.reset(new Vector<double>{par->total_pixel});
        wat->Ps.reset(new Vector<double>{par->total_pixel});
    }

    wat->h_sup.reset(new Vector<double>{par->total_pixel});

    /**************************************************************************************************/
    /*! Initialization of the struct "snow" (of the type SNOW): */
    /*************************************************************************************************/
    snow->S=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
    snow->S = new STATEVAR_3D{(double)number_novalue, par->max_snow_layers, Nr, Nc};

    // initial snow depth
    if ( strcmp(files[fsn0], string_novalue) != 0
         && strcmp(files[fswe0], string_novalue) != 0 )
    {
        printf("Initial condition on snow depth from file %s\n",files[fsn0]);
        M=read_map(2, files[fsn0], land->LC, UV, (double)number_novalue);
        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                snow->S->Dzl->co[1][r][c] = M->co[r][c];
            }
        }
        free_doublematrix(M);

        printf("Initial condition on snow water equivalent from file %s\n", files[fswe0]);
        M=read_map(2, files[fswe0], land->LC, UV, (double)number_novalue);
        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                snow->S->w_ice->co[1][r][c] = M->co[r][c];
            }
        }
        free_doublematrix(M);

    }
    else if ( strcmp(files[fsn0], string_novalue) != 0 )
    {
        printf("Initial condition on snow depth from file %s\n",files[fsn0]);
        M=read_map(2, files[fsn0], land->LC, UV, (double)number_novalue);
        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                snow->S->Dzl->co[1][r][c] = M->co[r][c];
            }
        }
        free_doublematrix(M);

        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                if ((long)land->LC->co[r][c] != number_novalue) snow->S->w_ice->co[1][r][c] =
                                                                        snow->S->Dzl->co[1][r][c] *
                                                                        IT->rhosnow0/rho_w;
            }
        }

    }
    else if ( strcmp(files[fswe0], string_novalue) != 0 )
    {
        printf("Initial condition on snow water equivalent from file %s\n", files[fswe0]);
        M=read_map(2, files[fswe0], land->LC, UV, (double)number_novalue);
        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                snow->S->w_ice->co[1][r][c] = M->co[r][c];
            }
        }
        free_doublematrix(M);

        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                if ((long)land->LC->co[r][c] != number_novalue) snow->S->Dzl->co[1][r][c] =
                                                                        snow->S->w_ice->co[1][r][c] *
                                                                        rho_w/IT->rhosnow0;
            }
        }
    }
    else
    {

        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                if ((long)land->LC->co[r][c] != number_novalue)
                {
                    snow->S->w_ice->co[1][r][c] = IT->swe0;
                    snow->S->Dzl->co[1][r][c] = IT->swe0*rho_w/IT->rhosnow0;
                }
            }
        }
    }


    // Optional reading of snow age in the whole basin
    if ( strcmp(files[fsnag0], string_novalue) != 0 )
    {
        printf("Snow age initial condition from file %s\n",files[fsnag0]+1);
        snow->age = read_map_vector(2, files[fsnag0], land->LC, UV,
                                    (double)number_novalue, top->rc_cont);
    }
    else
    {
        snow->age.reset(new Vector<double>{par->total_pixel});
        *(snow->age) = IT->agesnow0;
    }


    if (times->JD_plots->nh > 1)
    {
        if (strcmp(files[pD], string_novalue) != 0)
        {
            snow->Dplot.reset(new Vector<double>{par->total_pixel});
        }
    }

    if (par->blowing_snow==1)
    {

        snow->S_for_BS=(STATEVAR_1D *)malloc(sizeof(STATEVAR_1D));
        allocate_and_initialize_statevar_1D(snow->S_for_BS, (double)number_novalue,
                                            par->max_snow_layers);

        snow->change_dir_wind.reset(new Vector<long>{Fmaxlong(Nr,Nc)});

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

        if (par->output_snow_bin == 1)
        {
            snow->Wtrans_plot=new_doublematrix(Nr,Nc);
            snow->Wsubl_plot=new_doublematrix(Nr,Nc);

            for (r=1; r<=Nr; r++)
            {
                for (c=1; c<=Nc; c++)
                {
                    if ((long)land->LC->co[r][c]==number_novalue)
                    {
                        snow->Wtrans_plot->co[r][c]=(double)number_novalue;
                        snow->Wsubl_plot->co[r][c]=(double)number_novalue;

                    }
                    else
                    {
                        snow->Wtrans_plot->co[r][c]=0.0;
                        snow->Wsubl_plot->co[r][c]=0.0;

                    }
                }
            }
        }
    }

    if (par->output_snow_bin == 1)
    {
        if (strcmp(files[fsnowmelt], string_novalue) != 0)
        {
            snow->MELTED.reset(new Vector<double>{par->total_pixel});
            snow->melted.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fsnowsubl], string_novalue) != 0)
        {
            snow->SUBL.reset(new Vector<double>{par->total_pixel});
            snow->subl.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fsndur], string_novalue) != 0)
        {
            snow->t_snow.reset(new Vector<double>{par->total_pixel});
            snow->yes.reset(new Vector<short>{par->total_pixel});
        }
    }

    for (r=1; r<=Nr; r++)
    {
        for (c=1; c<=Nc; c++)
        {

            if ( (long)land->LC->co[r][c]!=number_novalue)
            {

                // Adjusting snow init depth in case of steep slope (contribution by Stephan Gruber)
                if (par->snow_curv > 0 && top->slope->co[r][c] > par->snow_smin)
                {
                    if (top->slope->co[r][c] <= par->snow_smax)
                    {
                        k_snowred = ( exp(-pow(top->slope->co[r][c] - par->snow_smin,
                                               2.)/par->snow_curv) -
                                      exp(-pow(par->snow_smax, 2.)/par->snow_curv) );
                    }
                    else
                    {
                        k_snowred = 0.0;
                    }
                    snow->S->Dzl->co[1][r][c] *= k_snowred;
                    snow->S->w_ice->co[1][r][c] *= k_snowred;
                }

                D = snow->S->Dzl->co[1][r][c];
                SWE = snow->S->w_ice->co[1][r][c];

                if (D<0 || SWE<0)
                {
                    f = fopen(FailedRunFile, "w");
                    fprintf(f, "Error: negative initial snow depth %e or snow water equivalent %e\n",
                            D,SWE);
                    fclose(f);
                    t_error("Fatal Error! Geotop is closed. See failing report (7).");

                }
                else if (D<1.E-5 && SWE>1.E-5)
                {
                    f = fopen(FailedRunFile, "w");
                    fprintf(f,"Error: Initial snow water equivalent %e > 0 and initial snow depth %e\n",
                            SWE,D);
                    fclose(f);
                    t_error("Fatal Error! Geotop is closed. See failing report (8).");

                }
                else if (D>1.E-5 && SWE<1.E-5)
                {
                    f = fopen(FailedRunFile, "w");
                    fprintf(f,"Error: Initial snow depth %e > 0 and initial snow water equivalent %e\n",
                            D, SWE);
                    fclose(f);
                    t_error("Fatal Error! Geotop is closed. See failing report (9).");

                }
                else if (D>1.E-5 || SWE>1.E-5)
                {

                    snow->age->co[top->j_cont[r][c]]*=86400.0;  // now in [s]

                    if (SWE <= par->max_weq_snow * par->max_snow_layers )
                    {
                        i = floor( SWE/par->max_weq_snow );

                        if (i>0)
                        {

                            for (n=1; n<=i; n++)
                            {
                                snow->S->w_ice->co[n][r][c] = par->max_weq_snow;
                            }

                            if (SWE - i * par->max_weq_snow > 0.1 * par->max_weq_snow)
                            {
                                snow->S->w_ice->co[i+1][r][c] = SWE - i * par->max_weq_snow;
                                snow->S->lnum->co[r][c] = i+1;
                            }
                            else
                            {
                                snow->S->w_ice->co[i][r][c] += (SWE - i * par->max_weq_snow);
                                snow->S->lnum->co[r][c] = i;
                            }
                        }
                        else
                        {
                            snow->S->w_ice->co[1][r][c] = SWE;
                            snow->S->lnum->co[r][c] = 1;
                        }
                    }
                    else
                    {
                        snow->S->lnum->co[r][c] = par->max_snow_layers;

                        for (n=1; n<=par->max_snow_layers; n++)
                        {
                            a = 0;

                            for (i=1; i<=par->inf_snow_layers->nh; i++)
                            {
                                if ( n == abs((*par->inf_snow_layers)(i)) )
                                    a = 1;
                            }

                            if (a == 0)
                            {
                                snow->S->w_ice->co[n][r][c] = par->max_weq_snow;
                            }
                            else
                            {
                                snow->S->w_ice->co[n][r][c] = ( SWE - par->max_weq_snow *
                                                                      ( par->max_snow_layers -
                                                                        par->inf_snow_layers->nh ) ) /
                                                              par->inf_snow_layers->nh;
                            }
                        }
                    }

                    for (n=1; n<=snow->S->lnum->co[r][c]; n++)
                    {
                        snow->S->Dzl->co[n][r][c] = D * (snow->S->w_ice->co[n][r][c] / SWE);
                        snow->S->T->co[n][r][c] = IT->Tsnow0;
                    }
                }

                non_dimensionalize_snowage(&(snow->age->co[top->j_cont[r][c]]), IT->Tsnow0);

                if (par->point_sim == 1)
                {
                    maxSWE = par->maxSWE->co[r][c];
                }
                else
                {
                    maxSWE = 1.E10;
                }

                f = fopen(logfile, "a");
                snow_layer_combination(par->alpha_snow, r, c, snow->S, -0.1,
                                       par->inf_snow_layers.get(), par->max_weq_snow, maxSWE, f);
                fclose(f);

            }
        }
    }

    if (recovered > 0)
    {
        initialize_shortmatrix(snow->S->type, 2);

        assign_recovered_map_long(old, par->recover, files[rns], snow->S->lnum, par, land->LC);
        // initialize_longmatrix(snow->S->lnum, snow->S->Dzl->ndh);

        assign_recovered_map_vector(old, par->recover, files[rsnag], snow->age.get(),
                                    top->rc_cont, par, land->LC);
        assign_recovered_tensor(old, par->recover, files[rDzs], snow->S->Dzl, par,
                                land->LC);
        assign_recovered_tensor(old, par->recover, files[rwls], snow->S->w_liq, par,
                                land->LC);
        assign_recovered_tensor(old, par->recover, files[rwis], snow->S->w_ice, par,
                                land->LC);
        assign_recovered_tensor(old, par->recover, files[rTs], snow->S->T, par,
                                land->LC);

        /* f = fopen(logfile, "a");
       for(r=1;r<=Nr;r++){
       for(c=1;c<=Nc;c++){
       if( (long)land->LC->co[r][c]!=number_novalue){
       snow_layer_combination(par->alpha_snow, r, c, snow->S, 0., par->inf_snow_layers,
       par->max_weq_snow, maxSWE, f);
       }
       }
       }
       fclose(f);*/
    }



    /**************************************************************************************************/
    /*! Initialization of the struct "glac" (of the type GLACIER):*/
    /**************************************************************************************************/
    /*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
    if ( par->point_sim!=1 && strcmp(files[fgl0], string_novalue) != 0
         && par->max_glac_layers==0)
    {
        printf("Warning: Glacier map present, but glacier represented with 0 layers\n");
        fprintf(flog,
                "Warning: Glacier map present, but glacier represented with 0 layers\n");
    }

    if (par->max_glac_layers==0 && IT->Dglac0>0)
    {
        printf("\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, \
but you assigned a value of the glacier depth. The latter will be ignored.\n");
        fprintf(flog, "\nWARNING: You have chosen 0 glacier layers in block 10 in the parameter file, \
but you assigned a value of the glacier depth. The latter will be ignored.\n");
    }

    // If the max number of glacier layers is greater than 1, the matrices (or tensors) lnum,
    // Dzl. w_liq, w_ice, T and print matrices are defined, according to the respective flags
    if (par->max_glac_layers>0)
    {

        if ( par->point_sim!=1 && strcmp(files[fgl0], string_novalue) != 0 )
        {
            printf("Glacier initial condition from file %s\n",files[fgl0]+1);
            fprintf(flog,"Glacier initial condition from file %s\n",files[fgl0]+1);
            M=read_map(2, files[fgl0], land->LC, UV, (double)number_novalue);
        }
        else
        {
            M=copydoublematrix_const(IT->Dglac0, land->LC, (double)number_novalue);
        }

        glac->G = new STATEVAR_3D {(double)number_novalue, par->max_glac_layers, Nr, Nc};

        if (par->output_glac_bin == 1)
        {
            if (strcmp(files[fglacmelt], string_novalue) != 0)
            {
                glac->MELTED.reset(new Vector<double>{par->total_pixel});
                glac->melted.reset(new Vector<double>{par->total_pixel});
            }
            if (strcmp(files[fglacsubl], string_novalue) != 0)
            {
                glac->SUBL.reset(new Vector<double>{par->total_pixel});
                glac->subl.reset(new Vector<double>{par->total_pixel});
            }
        }

        for (r=1; r<=Nr; r++)
        {
            for (c=1; c<=Nc; c++)
            {
                if ( (long)land->LC->co[r][c]!=number_novalue)
                {

                    if (M->co[r][c]<0)
                    {
                        f = fopen(FailedRunFile, "w");
                        fprintf(f, "Error: negative glacier data\n");
                        fclose(f);
                        t_error("Fatal Error! Geotop is closed. See failing report (10).");

                    }
                    else if (M->co[r][c]>1.E-5)
                    {

                        if (IT->rhoglac0 * M->co[r][c] / rho_w < par->max_weq_glac *
                                                                 par->max_glac_layers )
                        {

                            n = 0;
                            z = 0.;

                            do
                            {
                                n++;

                                if (IT->rhoglac0 * M->co[r][c] / rho_w < par->max_weq_glac * n)
                                {
                                    glac->G->w_ice->co[n][r][c] = par->max_weq_glac;
                                }
                                else
                                {
                                    glac->G->w_ice->co[n][r][c] = IT->rhoglac0 * M->co[r][c] /
                                                                  1000. - z;
                                }

                                glac->G->Dzl->co[n][r][c] = rho_w * glac->G->w_ice->co[n][r][c] /
                                                            IT->rhoglac0;
                                glac->G->T->co[n][r][c] = IT->Tglac0;

                                z += glac->G->w_ice->co[n][r][c];

                            }
                            while (fabs(z - IT->rhoglac0 * M->co[r][c] / rho_w) < 1.E-6);

                            glac->G->lnum->co[r][c] = n;

                        }
                        else
                        {

                            glac->G->lnum->co[r][c] = par->max_glac_layers;

                            for (n=1; n<=par->max_glac_layers; n++)
                            {
                                a = 0;

                                for (i=1; i<=par->inf_glac_layers->nh; i++)
                                {
                                    if (n == i)
                                        a = 1;
                                }

                                if (a == 0)
                                {
                                    glac->G->w_ice->co[n][r][c] = par->max_weq_glac;
                                }
                                else
                                {
                                    glac->G->w_ice->co[n][r][c] = ( IT->rhoglac0 * M->co[r][c] /
                                                                    rho_w -
                                                                    par->max_weq_glac *
                                                                    ( par->max_glac_layers -
                                                                      par->inf_glac_layers->nh ) ) /
                                                                  par->inf_glac_layers->nh;
                                }

                                glac->G->Dzl->co[n][r][c] = rho_w * glac->G->w_ice->co[n][r][c] /
                                                            IT->rhoglac0;
                                glac->G->T->co[n][r][c] = IT->Tglac0;

                            }
                        }
                    }

                    f = fopen(logfile, "a");
                    snow_layer_combination(par->alpha_snow, r, c, glac->G, -0.1,
                                           par->inf_glac_layers.get(), par->max_weq_glac, 1.E10, f);
                    fclose(f);

                }
            }
        }

        free_doublematrix(M);

        if (recovered > 0)
        {
            assign_recovered_map_long(old, par->recover, files[rni], glac->G->lnum, par, land->LC);
            assign_recovered_tensor(old, par->recover, files[rDzi], glac->G->Dzl, par, land->LC);
            assign_recovered_tensor(old, par->recover, files[rwli], glac->G->w_liq, par, land->LC);
            assign_recovered_tensor(old, par->recover, files[rwii], glac->G->w_ice, par, land->LC);
            assign_recovered_tensor(old, par->recover, files[rTi], glac->G->T, par, land->LC);
        }
    }



    //*************************************************************************************************
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

    if (par->output_meteo_bin == 1)
    {
        if (strcmp(files[fTa], string_novalue) != 0)
        {
            met->Tamean.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fwspd], string_novalue) != 0)
        {
            met->Vspdmean.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[fwdir], string_novalue) != 0)
        {
            met->Vdirmean.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[frh], string_novalue) != 0)
        {
            met->RHmean.reset(new Vector<double>{par->total_pixel});
        }
    }


    // plot output
    if (times->JD_plots->nh > 1)
    {
        if (strcmp(files[pTa], string_novalue) != 0)
        {
            met->Taplot.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pRH], string_novalue) != 0)
        {
            met->RHplot.reset(new Vector<double>{par->total_pixel});
        }
        if (strcmp(files[pVspd], string_novalue) != 0
            || strcmp(files[pVdir], string_novalue) != 0)
        {
            met->Vxplot.reset(new Vector<double>{par->total_pixel});
            met->Vyplot.reset(new Vector<double>{par->total_pixel});
        }
    }

    /**************************************************************************************************/
    // SpinUp variables

    if (par->recover > 0)
        write_suffix(rec, par->recover, 4);
    if (par->n_ContRecovery > 0)
        write_suffix(crec, par->n_ContRecovery, 5);

    if (par->Tzrun == 1)
    {
        sl->Tzrun = new_doublematrix(par->rc->nrh, Nl);
        initialize_doublematrix(sl->Tzrun, 0.);

        if (par->recover>0)
        {
            temp = join_strings(files[fTrun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fTrun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fTrun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point");
        for (l=1; l<=Nl; l++)
        {
            fprintf(f, ",l[%ld]",l);
        }
        fprintf(f, "\n");
        free(name);
        fclose(f);
    }

    if (par->wzrun == 1)
    {
        sl->wzrun = new_doublematrix(par->rc->nrh, Nl);
        initialize_doublematrix(sl->wzrun, 0.);

        if (par->recover>0)
        {
            temp = join_strings(files[fwrun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fwrun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fwrun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point");
        for (l=1; l<=Nl; l++)
        {
            fprintf(f, ",l[%ld]",l);
        }
        fprintf(f, "\n");
        free(name);
        fclose(f);
    }

    if (par->Tzmaxrun == 1)
    {
        sl->Tzmaxrun = new_doublematrix(par->rc->nrh, Nl);
        initialize_doublematrix(sl->Tzmaxrun, -1.E99);

        if (par->recover>0)
        {
            temp = join_strings(files[fTmaxrun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fTmaxrun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fTmaxrun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point");
        for (l=1; l<=Nl; l++)
        {
            fprintf(f, ",l[%ld]",l);
        }
        fprintf(f, "\n");
        free(name);
        fclose(f);
    }

    if (par->wzmaxrun == 1)
    {
        sl->wzmaxrun = new_doublematrix(par->rc->nrh, Nl);
        initialize_doublematrix(sl->wzmaxrun, -1.E99);

        if (par->recover>0)
        {
            temp = join_strings(files[fwmaxrun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fwmaxrun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fwmaxrun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point");
        for (l=1; l<=Nl; l++)
        {
            fprintf(f, ",l[%ld]",l);
        }
        fprintf(f, "\n");
        free(name);
        fclose(f);
    }

    if (par->Tzminrun == 1)
    {
        sl->Tzminrun = new_doublematrix(par->rc->nrh, Nl);
        initialize_doublematrix(sl->Tzminrun, 1.E99);

        if (par->recover>0)
        {
            temp = join_strings(files[fTminrun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fTminrun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fTminrun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point");
        for (l=1; l<=Nl; l++)
        {
            fprintf(f, ",l[%ld]",l);
        }
        fprintf(f, "\n");
        free(name);
        fclose(f);
    }

    if (par->wzminrun == 1)
    {
        sl->wzminrun = new_doublematrix(par->rc->nrh, Nl);
        initialize_doublematrix(sl->wzminrun, 1.E99);

        if (par->recover>0)
        {
            temp = join_strings(files[fwminrun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fwminrun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fwminrun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point");
        for (l=1; l<=Nl; l++)
        {
            fprintf(f, ",l[%ld]",l);
        }
        fprintf(f, "\n");
        free(name);
        fclose(f);
    }

    if (par->dUzrun == 1)
    {
        sl->dUzrun = new_doublematrix(par->rc->nrh, Nl);
        initialize_doublematrix(sl->dUzrun, 0.);

        if (par->recover>0)
        {
            temp = join_strings(files[fdUrun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fdUrun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fdUrun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point");
        for (l=1; l<=Nl; l++)
        {
            fprintf(f, ",l[%ld]",l);
        }
        fprintf(f, "\n");
        free(name);
        fclose(f);
    }

    if (par->SWErun == 1)
    {
        sl->SWErun = new_doublematrix(par->rc->nrh, 3); //mean,max,min
        for (l=1; l<=par->rc->nrh; l++)
        {
            sl->SWErun->co[l][1] = 0.;
            sl->SWErun->co[l][2] = -1.E99;
            sl->SWErun->co[l][3] = 1.E99;
        }

        if (par->recover>0)
        {
            temp = join_strings(files[fSWErun], rec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else if (par->n_ContRecovery>0)
        {
            temp = join_strings(files[fSWErun], crec);
            name = join_strings(temp, textfile);
            free(temp);
        }
        else
        {
            name = join_strings(files[fSWErun], textfile);
        }

        f = fopen(name, "w");
        fprintf(f, "Period,Run,Point,SWEmean[mm],SWEmax[mm],SWEmin[mm]\n");
        free(name);
        fclose(f);
    }

    if (recovered > 0)
    {
        if (par->Tzrun == 1)
            recover_run_averages(old, sl->Tzrun, files[rTrun], land->LC, top->rc_cont, par, Nl);

        if (par->wzrun == 1)
            recover_run_averages(old, sl->wzrun, files[rwrun], land->LC, top->rc_cont, par, Nl);

        if (par->Tzminrun == 1)
            recover_run_averages(old, sl->Tzminrun, files[rTminrun], land->LC, top->rc_cont, par, Nl);

        if (par->wzminrun == 1)
            recover_run_averages(old, sl->wzminrun, files[rwminrun], land->LC, top->rc_cont, par, Nl);

        if (par->Tzmaxrun == 1)
            recover_run_averages(old, sl->Tzmaxrun, files[rTmaxrun], land->LC, top->rc_cont, par, Nl);

        if (par->wzmaxrun == 1)
            recover_run_averages(old, sl->wzmaxrun, files[rwmaxrun], land->LC, top->rc_cont, par, Nl);

        if (par->dUzrun == 1)
            recover_run_averages(old, sl->dUzrun, files[rdUrun], land->LC, top->rc_cont, par, Nl);

        if (par->SWErun == 1)
            recover_run_averages(old, sl->SWErun, files[rSWErun], land->LC, top->rc_cont, par, 3);
    }

    // WRITE INITIAL CONDITION
    write_output_headers(met->st->Z->nh, times, wat, par, top, land, sl, egy, snow, glac);

    if (par->state_pixel == 1)
    {
        for (j=1; j<=par->rc->nrh; j++)
        {

            r = par->rc->co[j][1];
            c = par->rc->co[j][2];

            if (par->output_vertical_distances == 1)
            {
                cosslope = cos( Fmin(max_slope, top->slope->co[r][c]) * Pi/180. );
            }
            else
            {
                cosslope = 1.;
            }

            write_soil_output(j, (*par->IDpoint)(j), par->init_date->co[1],
                              par->end_date->co[1], par->init_date->co[1], JD, day, month, year, hour,
                              minute, par->soil_plot_depths.get(), sl, par, (double)PsiMin, cosslope);
            write_snow_output(j, (*par->IDpoint)(j), r, c, par->init_date->co[1],
                              par->end_date->co[1], par->init_date->co[1], JD, day, month, year, hour,
                              minute, par->snow_plot_depths.get(), snow->S, par, cosslope);
        }
    }

    /**************************************************************************************************/
    // Free the struct allocated in this subroutine:
    free_doublematrix(par->chkpt);
    free_doublematrix(IT->bed);
    free_doubletensor(IT->pa_bed);

    for (i=0; i<nmet; i++)
    {
        free(IT->met_col_names[i]);
    }
    free(IT->met_col_names);

    for (i=0; i<nsoilprop; i++)
    {
        free(IT->soil_col_names[i]);
    }
    free(IT->soil_col_names);

    for (i=0; i<2; i++)
    {
        free(IT->horizon_col_names[i]);
    }
    free(IT->horizon_col_names);

    for (i=0; i<ptTOT; i++)
    {
        free(IT->point_col_names[i]);
    }
    free(IT->point_col_names);

    for (i=0; i<nlstot; i++)
    {
        free(IT->lapserates_col_names[i]);
    }
    free(IT->lapserates_col_names);

    for (i=0; i<8; i++)
    {
        free(IT->meteostations_col_names[i]);
    }
    free(IT->meteostations_col_names);

    n = Fminlong((*par->Nl_spinup)(i_sim0),Nl);

    if (par->point_sim != 1)
    {
        cont_nonzero_values_matrix2(&i, &j, cnet, land->LC, top->lrc_cont,
                                    top->i_cont, par->total_pixel, par->total_channel, n);
        top->Li.reset(new Vector<long>{i});
        top->Lp.reset(new Vector<long>{j});
        wat->Lx.reset(new Vector<double> {i});
        cont_nonzero_values_matrix3(top->Lp.get(), top->Li.get(), cnet, land->LC, top->lrc_cont,
                                    top->i_cont, par->total_pixel, par->total_channel, n);
    }
    else
    {
        i = n;
        j = n+1;
        top->Li.reset(new Vector<long>{i});
        top->Lp.reset(new Vector<long>{j});
        wat->Lx.reset(new Vector<double>{i});
        for (l=1; l<=n; l++)
        {
            (*top->Li)(l) = l+1;
            (*top->Lp)(l) = l;
        }
        (*top->Lp)(j) = i;
    }

    wat->H0.reset(new Vector<double>{j});
    wat->H1.reset(new Vector<double>{j});
    wat->dH.reset(new Vector<double>{j});
    wat->B.reset(new Vector<double>{j});
    wat->f.reset(new Vector<double>{j});
    wat->df.reset(new Vector<double>{j});

    wat->Kbottom = new_doublematrix(Nr, Nc);
    initialize_doublematrix(wat->Kbottom, 0.);

    wat->Klat = new_doublematrix(top->BC_DepthFreeSurface->nh, Nl);
    initialize_doublematrix(wat->Klat, 0.);

    fclose(flog);

}

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

void read_inputmaps(TOPO *top, LAND *land, SOIL *sl, PAR *par, INIT_TOOLS *IT,
                    FILE *flog)
{
    GEOLOG_PREFIX(__func__);
    long r, c, i, cont;
    DOUBLEMATRIX *M;
    SHORTMATRIX *curv;
    short flag;
    char *temp;
    double min, max;
    FILE *f;

    // reading TOPOGRAPHY
    flag = file_exists(fdem);
    if (flag == 1)
    {
        M=new_doublematrix(1,1);
        top->Z0=read_map(0, files[fdem], M, UV, (double)number_novalue); // topography
        free_doublematrix(M);
        write_map(files[fdem], 0, par->format_out, top->Z0, UV, number_novalue);

        // filtering
        M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
        multipass_topofilter(par->lowpass, top->Z0, M, (double)number_novalue, 1);
        copy_doublematrix(M, top->Z0);
        free_doublematrix(M);
        // write_map(files[fdem], 0, par->format_out, top->Z0, UV, number_novalue);

        // calculate East and North
        top->East = new_doublematrix(top->Z0->nrh, top->Z0->nch);
        top->North = new_doublematrix(top->Z0->nrh, top->Z0->nch);
        for (r=1; r<=top->Z0->nrh; r++)
        {
            for (c=1; c<=top->Z0->nch; c++)
            {
                top->East->co[r][c] = UV->U->co[4] + (c-0.5)*UV->U->co[2];
                top->North->co[r][c] = UV->U->co[3] + (top->Z0->nrh-(r-0.5))*UV->U->co[1];
            }
        }

    }
    else
    {

        f = fopen(FailedRunFile, "w");
        fprintf(f, "Error: It is impossible to proceed without giving the digital elevation model\n");
        fclose(f);
        t_error("Fatal Error! Geotop is closed. See failing report (11).");

    }

    // reading LAND COVER TYPE
    flag = file_exists(flu);
    if (flag == 1)
    {
        land->LC=read_map(1, files[flu], top->Z0, UV, (double)number_novalue);
        // Check borders
        for (r=1; r<=land->LC->nrh; r++)
        {
            land->LC->co[r][1]=(double)number_novalue;
            land->LC->co[r][land->LC->nch]=(double)number_novalue;
        }
        for (c=1; c<=land->LC->nch; c++)
        {
            land->LC->co[1][c]=(double)number_novalue;
            land->LC->co[land->LC->nrh][c]=(double)number_novalue;
        }
        for (r=1; r<=land->LC->nrh; r++)
        {
            for (c=1; c<=land->LC->nch; c++)
            {
                if ((long)land->LC->co[r][c] != number_novalue)
                {
                    if ((long)land->LC->co[r][c] < 1
                        || (long)land->LC->co[r][c] > par->n_landuses)
                    {
                        f = fopen(FailedRunFile, "w");
                        fprintf(f, "Error: It is not possible to assign Value < 1 or > n_landuses \
to the land cover type\n");
                        fclose(f);
                        t_error("Fatal Error! Geotop is closed. See failing report (12).");
                    }
                }
            }
        }

        // Land use is the official mask
        for (r=1; r<=land->LC->nrh; r++)
        {
            for (c=1; c<=land->LC->nch; c++)
            {
                if ((long)land->LC->co[r][c]!=number_novalue)
                {
                    if ((long)top->Z0->co[r][c]==number_novalue)
                    {
                        printf("ERROR Land use mask include DTM novalue pixels");
                        printf("\nr:%ld c:%ld Z:%f landuse:%f\n",r,c,top->Z0->co[r][c],
                               land->LC->co[r][c]);
                        land->LC->co[r][c]=(double)number_novalue;
                        printf("LANDUSE set at novalue where DTM is not available\n");
                    }
                }
            }
        }

    }
    else
    {
        // Write land->LC (land cover)
        printf("Land cover type assumed to be always 1\n");
        land->LC=copydoublematrix_const(1.0, top->Z0, (double)number_novalue);

        for (r=1; r<=land->LC->nrh; r++)
        {
            land->LC->co[r][1]=(double)number_novalue;
            land->LC->co[r][land->LC->nch]=(double)number_novalue;
        }
        for (c=1; c<=land->LC->nch; c++)
        {
            land->LC->co[1][c]=(double)number_novalue;
            land->LC->co[land->LC->nrh][c]=(double)number_novalue;
        }
    }
    if (flag >= 0)
        write_map(files[flu], 1, par->format_out, land->LC, UV, number_novalue);

    if (par->state_pixel == 1)
    {
        par->rc=new_longmatrix(par->chkpt->nrh,2);
        par->IDpoint.reset(new Vector<long>{par->chkpt->nrh});
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            par->rc->co[i][1]=row(par->chkpt->co[i][ptY], top->Z0->nrh, UV, number_novalue);
            par->rc->co[i][2]=col(par->chkpt->co[i][ptX], top->Z0->nch, UV, number_novalue);

            if (par->rc->co[i][1] == number_novalue
                || par->rc->co[i][2] == number_novalue)
            {
                geolog << "Point #" << i << " is out of the domain";
                fprintf(flog, "Point #%4ld is out of the domain",i);
                fclose(flog);

                f = fopen(FailedRunFile, "w");
                fprintf(f, "Point #%4ld is out of the domain",i);
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report.");
            }

            if ((long)land->LC->co[par->rc->co[i][1]][par->rc->co[i][2]]==number_novalue)
            {
                printf("Point #%4ld corresponds to NOVALUE pixel",i);
                fprintf(flog, "Point #%4ld corresponds to NOVALUE pixel",i);
                fclose(flog);

                f = fopen(FailedRunFile, "w");
                fprintf(f, "Point #%4ld corresponds to NOVALUE pixel",i);
                fclose(f);
                t_error("Fatal Error! Geotop is closed. See failing report.");
            }

            if ((long)par->chkpt->co[i][ptID]!=number_novalue)
            {
                (*par->IDpoint)(i)=(long)par->chkpt->co[i][ptID];
            }
            else
            {
                (*par->IDpoint)(i)=i;
            }
        }
    }

    /*************************************************************************************************/
    // reading SKY VIEW FACTOR
    flag = file_exists(fsky);
    if (flag == 1)
    {
        top->sky=read_map(2, files[fsky], land->LC, UV, (double)number_novalue);
    }
    else  /* The sky view factor file "top->sky" must be calculated */
    {
        top->sky = new_doublematrix(top->Z0->nrh,top->Z0->nch);
        if (par->sky == 0)
        {
            initialize_doublematrix(top->sky, 1.);
        }
        else
        {
            curv = new_shortmatrix(top->Z0->nrh,top->Z0->nch);
            nablaquadro_mask(top->Z0, curv, UV->U.get(), UV->V.get());
            sky_view_factor(top->sky, 36, UV, top->Z0, curv, number_novalue);
            free_shortmatrix(curv);
        }
    }
    if (flag >= 0)
        write_map(files[fsky], 0, par->format_out, top->sky, UV, number_novalue);

    /**************************************************************************************************/
    // reading DELAY
    flag = file_exists(fdelay);
    if (flag == 1)
    {
        land->delay = read_map(2, files[fdelay], land->LC, UV,
                               (double)number_novalue);
    }
    else
    {
        land->delay = new_doublematrix(top->Z0->nrh,top->Z0->nch);
        initialize_doublematrix(land->delay, 0);
    }
    if (flag >= 0)
        write_map(files[fdelay], 0, par->format_out, land->delay, UV, number_novalue);

    /**************************************************************************************************/
    // reading SOIL MAP
    flag = file_exists(fsoil);
    if (flag == 1)
    {
        M=read_map(2, files[fsoil], land->LC, UV, (double)number_novalue);
        sl->type=copylong_doublematrix(M);
        for (r=1; r<=land->LC->nrh; r++)
        {
            for (c=1; c<=land->LC->nch; c++)
            {
                if ((long)land->LC->co[r][c] != number_novalue)
                {
                    if (sl->type->co[r][c] < 1 || sl->type->co[r][c] > par->nsoiltypes)
                    {
                        f = fopen(FailedRunFile, "w");
                        fprintf(f, "Error: It is not possible to assign Value < 1 or > nsoiltypes \
to the soil type map");
                        fclose(f);
                        t_error("Fatal Error! Geotop is closed. See failing report (13).");
                    }
                }
            }
        }
    }
    else
    {
        M=copydoublematrix_const(par->soil_type_land_default, land->LC, (double)number_novalue);
        sl->type=copylong_doublematrix(M);
    }
    if (flag >= 0)
        write_map(files[fsoil], 1, par->format_out, M, UV, number_novalue);
    free_doublematrix(M);

    /**************************************************************************************************/
    // SLOPE
    top->dzdE = new_doublematrix(land->LC->nrh, land->LC->nch);
    top->dzdN = new_doublematrix(land->LC->nrh, land->LC->nch);
    find_slope(UV->U->co[1], UV->U->co[2], top->Z0, top->dzdE, top->dzdN,
               (double)number_novalue);

    flag = file_exists(fslp);
    if (flag == 1)
    {
        top->slope=read_map(2, files[fslp], land->LC, UV, (double)number_novalue); // reads in degrees
    }
    else
    {
        top->slope=find_max_slope(top->Z0, top->dzdE, top->dzdN, (double)number_novalue);
    }
    if (flag >= 0)
        write_map(files[fslp], 0, par->format_out, top->slope, UV, number_novalue);

    find_min_max(top->slope, (double)number_novalue, &max, &min);
    printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min, tan(max*Pi/180.),max);
    fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min,tan(max*Pi/180.),max);

    /**************************************************************************************************/
    // ASPECT
    flag = file_exists(fasp);
    if (flag == 1)
    {
        top->aspect=read_map(2, files[fasp], land->LC, UV, (double)number_novalue);
    }
    else
    {
        top->aspect=find_aspect(top->Z0, top->dzdE, top->dzdN, (double)number_novalue);
    }
    if (flag >= 0)
        write_map(files[fasp], 0, par->format_out, top->aspect, UV, number_novalue);

    /**************************************************************************************************/
    // curvature
    top->curvature1=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    top->curvature2=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    top->curvature3=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    top->curvature4=new_doublematrix(top->Z0->nrh,top->Z0->nch);

    // filtering
    M=new_doublematrix(top->Z0->nrh,top->Z0->nch);
    multipass_topofilter(par->lowpass_curvatures, top->Z0, M, (double)number_novalue, 1);
    curvature(UV->U->co[1], UV->U->co[2], M, top->curvature1, top->curvature2,
              top->curvature3, top->curvature4, (double)number_novalue);
    free_doublematrix(M);

    if (strcmp(files[fcurv], string_novalue) != 0)
    {
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

    /**************************************************************************************************/
    /*
    //Channel network (in top->pixel_type)

    //pixel type = 0 land pixel (if it is on the border, the border is impermeable,
    water is free only on the surface)
    //pixel type = 1 or 2 land pixel (it it is on the border, the border is permeable above an
    user-defined elevation in the saturated part, weir-wise)
    //pixel type = 10 channel pixel (if it is on the border, the border is impermeable,
    water is free only on the surface)
    //pixel type = -1 land pixel where an incoming discharge from outside is considered (as rain)
    */

    flag = file_exists(fnet);
    if (flag == 1)
    {
        M=read_map(2, files[fnet], land->LC, UV, (double)number_novalue);
        top->pixel_type=copyshort_doublematrix(M);

        cont = 0;
        for (r=1; r<=top->Z0->nrh; r++)
        {
            for (c=1; c<=top->Z0->nch; c++)
            {
                if ((long)land->LC->co[r][c]!=number_novalue)
                {
                    if (top->pixel_type->co[r][c]!=0 && top->pixel_type->co[r][c]!=1
                        && top->pixel_type->co[r][c]!=2 && top->pixel_type->co[r][c]!=10
                        && top->pixel_type->co[r][c]!=11 && top->pixel_type->co[r][c]!=12
                        && top->pixel_type->co[r][c]!=-1)
                    {
                        f = fopen(FailedRunFile, "w");
                        fprintf(f, "Error: Only the following values are admitted in the network map: \
-1, 0, 1, 2, 10\n");
                        fclose(f);
                        t_error("Fatal Error! Geotop is closed. See failing report (14).");
                    }
                    if (top->pixel_type->co[r][c]==10) cont++;
                }
            }
        }

        printf("Channel networks has %ld pixels set to channel\n",cont);

        if (flag >= 0) write_map(files[fnet], 1, par->format_out, M, UV, number_novalue);
        free_doublematrix(M);
    }
    else
    {
        top->pixel_type=new_shortmatrix(land->LC->nrh, land->LC->nch);
        initialize_shortmatrix(top->pixel_type, 0);
    }

    /**************************************************************************************************/
    // border
    top->is_on_border=new_shortmatrix(land->LC->nrh, land->LC->nch);
    for (r=1; r<=land->LC->nrh; r++)
    {
        for (c=1; c<=land->LC->nch; c++)
        {
            if ( (long)land->LC->co[r][c]!=number_novalue)
            {
                top->is_on_border->co[r][c] = is_boundary(r, c, land->LC, number_novalue);
            }
            else
            {
                top->is_on_border->co[r][c] = -1;
            }
        }
    }

    // count the pixels having pixel_type = 1, 2 or -1
    cont = 0;
    for (r=1; r<=top->Z0->nrh; r++)
    {
        for (c=1; c<=top->Z0->nch; c++)
        {
            if (top->is_on_border->co[r][c]==1)
            {
                if (top->pixel_type->co[r][c] == -1 || top->pixel_type->co[r][c] == 1
                    || top->pixel_type->co[r][c] == 2 || top->pixel_type->co[r][c] == 11
                    || top->pixel_type->co[r][c] == 12) cont ++;
            }
        }
    }

    top->BC_counter = new_longmatrix(top->Z0->nrh, top->Z0->nch);
    initialize_longmatrix(top->BC_counter, 0);

    if (cont > 0)
    {
        top->BC_DepthFreeSurface.reset(new Vector<double>{cont});
        cont = 0;
        for (r=1; r<=top->Z0->nrh; r++)
        {
            for (c=1; c<=top->Z0->nch; c++)
            {
                if (top->is_on_border->co[r][c]==1)
                {
                    if (top->pixel_type->co[r][c] == -1 || top->pixel_type->co[r][c] == 1
                        || top->pixel_type->co[r][c] == 2 || top->pixel_type->co[r][c] == 11
                        || top->pixel_type->co[r][c] == 12)
                    {
                        cont ++;
                        top->BC_counter->co[r][c] = cont;
                        top->BC_DepthFreeSurface->co[cont] = par->DepthFreeSurface; //[mm]
                    }
                }
            }
        }
    }
    else
    {
        top->BC_DepthFreeSurface.reset(new Vector<double>{1});
        *(top->BC_DepthFreeSurface) = double(number_novalue);
    }

    // bedrock
    flag = file_exists(fbed);
    if (flag == 1)
    {
        IT->bed=read_map(2, files[fbed], land->LC, UV, (double)number_novalue);
    }
    else
    {
        IT->bed=new_doublematrix(top->Z0->nrh, top->Z0->nch);
        initialize_doublematrix(IT->bed, 1.E99);
    }
    if (flag>=0)
        write_map(files[fbed], 0, par->format_out, IT->bed, UV, number_novalue);
}

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

void read_optionsfile_point(PAR *par, TOPO *top, LAND *land, SOIL *sl,
                            TIMES *times, INIT_TOOLS *IT, FILE *flog)
{
    GEOLOG_PREFIX(__func__);
    long i, r, c, num_lines;
    DOUBLEMATRIX *Q=NULL, *P=NULL, *R=NULL, *S=NULL, *T=NULL, *Z=NULL, *LU=NULL; // ec 2012 08 22
    SHORTMATRIX *curv;
    short read_dem, read_lu, read_soil, read_sl, read_as, read_sk, read_bed,
            read_curv, flag, coordinates;
    char *temp;
    double min, max;
    FILE *f;

    // 4. CALCULATE TOPOGRAPHIC PROPERTIES
    // check if there are point coordinates
    coordinates = 1;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ( (long)par->chkpt->co[i][ptX]==number_novalue
             || (long)par->chkpt->co[i][ptY]==number_novalue ) coordinates = 0;
    }
    /* if (coordinates == 0 && par->recover>0){
       printf("Warning: Not possible to recover the simulation because at least one point has no \
       coordinates\n");
       printf("Starting from normal initial condition\n");
       fprintf(flog,"Warning: Not possible to recover the simulation because at least one point has no \
       coordinates\n");
       fprintf(flog,"Starting from normal initial condition\n");
       par->recover = 0;
       } */

    // ---------------------- (a) Read dem ----------------------
    read_dem=0;
    // if(par->recover>0) read_dem=1;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptLC]==number_novalue
            || (long)par->chkpt->co[i][ptSY]==number_novalue ||
            (long)par->chkpt->co[i][ptS]==number_novalue
            || (long)par->chkpt->co[i][ptA]==number_novalue ||
            (long)par->chkpt->co[i][ptCNS]==number_novalue
            || (long)par->chkpt->co[i][ptCWE]==number_novalue ||
            (long)par->chkpt->co[i][ptCNwSe]==number_novalue
            || (long)par->chkpt->co[i][ptCNeSw]==number_novalue)
        {
            read_dem=1;
        }
    }
    if (read_dem == 1 && coordinates == 0)
    {
        printf("Warning: \
Not possible to read from dem because at least one point has no coordinates\n");
        fprintf(flog, "Warning: \
Not possible to read from dem because at least one point has no coordinates\n");
        read_dem = 0;
    }
    if (read_dem==1)
    {
        flag = file_exists(fdem);
        if (flag == 1)
        {
            printf("Warning: Dem file %s present\n",files[fdem]+1);
            fprintf(flog,"Warning: Dem file %s present\n",files[fdem]+1);

            Q=new_doublematrix(1,1);
            Z=read_map(0, files[fdem], Q, UV, (double)number_novalue); //topography
            free_doublematrix(Q);

            Q=new_doublematrix(Z->nrh,Z->nch);
            multipass_topofilter(par->lowpass, Z, Q, (double)number_novalue, 1);
            copy_doublematrix(Q, Z);
            free_doublematrix(Q);

        }
        else
        {
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


    if (read_dem==1)
    {
        par->r_points.reset(new Vector<long>{par->chkpt->nrh});
        par->c_points.reset(new Vector<long>{par->chkpt->nrh});
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            (*par->r_points)(i)=row(par->chkpt->co[i][ptY], Z->nrh, UV, number_novalue);
            (*par->c_points)(i)=col(par->chkpt->co[i][ptX], Z->nch, UV, number_novalue);
            if ((long)par->chkpt->co[i][ptZ]==number_novalue)
                par->chkpt->co[i][ptZ] = Z->co[(*par->r_points)(i)][(*par->c_points)(i)];
            printf("ok");
        }
    }

    // ---------------------- (b) Read land use ----------------------
    read_lu=0;
    // if(par->recover>0) read_lu=1;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptLC]==number_novalue)
            read_lu=1;
    }
    if (read_lu==1 && coordinates==0)
        read_lu=0;
    if (read_lu==1)
    {
        flag = file_exists(flu);
        if (flag == 1)
        {
            if (read_dem==0)
            {
                Q=new_doublematrix(1,1);
                LU=read_map(0, files[flu], Q, UV, (double)number_novalue);
                free_doublematrix(Q);
            }
            else
            {
                LU=read_map(1, files[flu], Z, UV, (double)number_novalue);
            }

        }
        else
        {
            printf("Warning: Landuse file not present, uniform cover considered\n");
            fprintf(flog,"Warning: Landuse file not present, uniform cover considered\n");
            if (read_dem==1)
            {
                LU=copydoublematrix_const(1.0, Z, (double)number_novalue);
            }
            else
            {
                read_lu=0;
                /* if(par->recover>0){
           printf("Warning: Not possible to recover the simulation because there is no dem\n");
           printf("Starting from normal initial condition\n");
           fprintf(flog,"Warning: \
           Not possible to recover the simulation because there is no dem\n");
           fprintf(flog,"Starting from normal initial condition\n");
           par->recover = 0;
           }*/
            }
        }
    }

    if (read_lu==1)
    {
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            if ((long)par->chkpt->co[i][ptLC]==number_novalue)
            {
                r=row(par->chkpt->co[i][ptY], LU->nrh, UV, number_novalue);
                c=col(par->chkpt->co[i][ptX], LU->nch, UV, number_novalue);
                par->chkpt->co[i][ptLC]=LU->co[r][c];
            }
        }
    }

    // ---------------------- (c) Read soil type ----------------------
    read_soil=0;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptSY]==number_novalue)
            read_soil=1;
    }
    if (read_soil==1 && coordinates==0)
        read_soil=0;
    if (read_soil==1)
    {
        flag = file_exists(fsoil);
        if (flag == 1)
        {
            if (read_dem==0)
            {
                Q=new_doublematrix(1,1);
                P=read_map(0, files[fsoil], Q, UV, (double)number_novalue);
                free_doublematrix(Q);
            }
            else
            {
                P=read_map(1, files[fsoil], Z, UV, (double)number_novalue);
            }

        }
        else
        {
            printf("Warning: Soiltype file not present\n");
            fprintf(flog,"Warning: Soiltype file not present\n");
            read_soil=0;
        }
    }
    if (read_soil==1)
    {
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            if ((long)par->chkpt->co[i][ptSY]==number_novalue)
            {
                r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
                c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
                par->chkpt->co[i][ptSY]=P->co[r][c];
            }
        }
        free_doublematrix(P);
    }

    // ---------------------- (d) Read slope ----------------------
    read_sl=0;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptS]==number_novalue)
            read_sl=1;
    }
    if (read_sl==1 && coordinates==0)
        read_sl=0;
    if (read_sl==1)
    {
        flag = file_exists(fslp);
        if (flag == 1)
        {
            if (read_dem==0)
            {
                Q=new_doublematrix(1,1);
                P=read_map(0, files[fslp], Q, UV, (double)number_novalue);
                free_doublematrix(Q);
            }
            else
            {
                P=read_map(1, files[fslp], Z, UV, (double)number_novalue);
            }

        }
        else
        {
            if (read_dem==0)
            {
                printf("Warning: Slopes file not present\n");
                fprintf(flog,"Warning: Slopes file not present\n");
                read_sl=0;
            }
            else
            {
                Q=new_doublematrix(Z->nrh,Z->nch);
                R=new_doublematrix(Z->nrh,Z->nch);
                find_slope(UV->U->co[1], UV->U->co[2], Z, Q, R, (double)number_novalue);
                P=find_max_slope(Z, Q, R, (double)number_novalue);
                free_doublematrix(Q);
                free_doublematrix(R);
                if (flag==0)
                    write_map(files[fslp], 0, par->format_out, P, UV, number_novalue);
            }
        }
    }

    if (read_sl==1)
    {
        find_min_max(P, (double)number_novalue, &max, &min);
        printf("Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min,
               tan(max*Pi/180.),max);
        fprintf(flog,"Slope Min:%f (%f deg) Max:%f (%f deg) \n",tan(min*Pi/180.),min,
                tan(max*Pi/180.),max);

        for (i=1; i<=par->chkpt->nrh; i++)
        {
            if ((long)par->chkpt->co[i][ptS]==number_novalue)
            {
                r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
                c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
                par->chkpt->co[i][ptS]=P->co[r][c];
            }
        }
        free_doublematrix(P);
    }

    // ---------------------- (e) Read aspect ----------------------
    read_as=0;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptA]==number_novalue)
            read_as=1;
    }
    if (read_as==1 && coordinates==0)
        read_as=0;
    if (read_as==1)
    {
        flag = file_exists(fasp);
        if (flag == 1)
        {
            if (read_dem==0)
            {
                Q=new_doublematrix(1,1);
                P=read_map(0, files[fasp], Q, UV, (double)number_novalue);
                free_doublematrix(Q);
            }
            else
            {
                P=read_map(1, files[fasp], Z, UV, (double)number_novalue);
            }
        }
        else
        {
            if (read_dem==0)
            {
                printf("Warning: Aspect file not present\n");
                fprintf(flog,"Warning: Aspect file not present\n");
                read_as=0;
            }
            else
            {
                Q=new_doublematrix(Z->nrh,Z->nch);
                R=new_doublematrix(Z->nrh,Z->nch);
                find_slope(UV->U->co[1], UV->U->co[2], Z, Q, R, (double)number_novalue);
                P=find_aspect(Z, Q, R, (double)number_novalue);
                free_doublematrix(Q);
                free_doublematrix(R);
                if (flag==0)
                    write_map(files[fasp], 0, par->format_out, P, UV, number_novalue);
            }
        }
    }

    if (read_as==1)
    {
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            if ((long)par->chkpt->co[i][ptA]==number_novalue)
            {
                r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
                c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
                par->chkpt->co[i][ptA]=P->co[r][c];
            }
        }
        free_doublematrix(P);
    }

    // ---------------------- (f) Sky view factor file ----------------------
    read_sk=0;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptSKY]==number_novalue)
            read_sk=1;
    }
    if (read_sk==1 && coordinates==0)
        read_sk=0;
    if (read_sk==1)
    {
        flag = file_exists(fsky);
        if (flag == 1)
        {
            if (read_dem==0)
            {
                Q=new_doublematrix(1,1);
                P=read_map(0, files[fsky], Q, UV, (double)number_novalue);
                free_doublematrix(Q);
            }
            else
            {
                P=read_map(1, files[fsky], Z, UV, (double)number_novalue);
            }
        }
        else
        {
            if (read_dem==0)
            {
                printf("Warning: Sky view factor file not present\n");
                fprintf(flog,"Warning: Sky view factor file not present\n");
                read_sk=0;
            }
            else
            {
                P=new_doublematrix(Z->nrh,Z->nch);
                curv=new_shortmatrix(Z->nrh,Z->nch);
                nablaquadro_mask(Z, curv, UV->U.get(), UV->V.get());
                sky_view_factor(P, 36, UV, Z, curv, number_novalue);
                free_shortmatrix(curv);
                if (flag==0) write_map(files[fsky], 0, par->format_out, P, UV,
                                       number_novalue);
            }
        }
    }

    if (read_sk==1)
    {
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            if ((long)par->chkpt->co[i][ptSKY]==number_novalue)
            {
                r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
                c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
                par->chkpt->co[i][ptSKY]=P->co[r][c];
            }
        }
        free_doublematrix(P);
    }

    // ---------------------- (f2) Bedrock file ----------------------
    read_bed=0;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptBED]==number_novalue)
            read_bed=1;
    }
    if (read_bed==1 && coordinates==0)
        read_bed=0;
    if (read_bed==1)
    {
        flag = file_exists(fbed);
        if (flag == 1)
        {
            if (read_dem==0)
            {
                Q=new_doublematrix(1,1);
                P=read_map(0, files[fbed], Q, UV, (double)number_novalue);
                free_doublematrix(Q);
            }
            else
            {
                P=read_map(1, files[fbed], Z, UV, (double)number_novalue);
            }
        }
        else
        {
            printf("Warning: Bedrock depth file not present\n");
            fprintf(flog,"Warning: Sky view factor file not present\n");
            read_bed=0;
        }
    }

    if (read_bed==1)
    {
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            if ((long)par->chkpt->co[i][ptBED]==number_novalue)
            {
                r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
                c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);
                par->chkpt->co[i][ptBED]=P->co[r][c];
            }
        }
        free_doublematrix(P);
    }

    // ---------------------- (g) Curvature ----------------------
    read_curv=0;
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ( (long)par->chkpt->co[i][ptCNS]==number_novalue
             || (long)par->chkpt->co[i][ptCWE]==number_novalue ||
             (long)par->chkpt->co[i][ptCNwSe]==number_novalue
             || (long)par->chkpt->co[i][ptCNeSw]==number_novalue  ) read_curv=1;
    }
    if (read_curv==1 && coordinates==0)
        read_curv=0;
    if (read_curv==1)
    {
        if (read_dem==0)
        {
            printf("Warning: \
Dem file is not present, and therefore it is not possible to calculate curvature\n");
            fprintf(flog, "Warning: \
Dem file is not present, and therefore it is not possible to calculate curvature\n");
            read_curv=0;
        }
        else
        {
            Q=new_doublematrix(Z->nrh,Z->nch);
            P=new_doublematrix(Z->nrh,Z->nch);
            R=new_doublematrix(Z->nrh,Z->nch);
            S=new_doublematrix(Z->nrh,Z->nch);
            T=new_doublematrix(Z->nrh,Z->nch);

            multipass_topofilter(par->lowpass_curvatures, Z, Q, (double)number_novalue, 1);
            curvature(UV->U->co[1], UV->U->co[2], Q, P, R, S, T, (double)number_novalue);
            free_doublematrix(Q);

            if (strcmp(files[fcurv],string_novalue) != 0)
            {
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
    if (read_curv==1)
    {
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            r=row(par->chkpt->co[i][ptY], P->nrh, UV, number_novalue);
            c=col(par->chkpt->co[i][ptX], P->nch, UV, number_novalue);

            if ((long)par->chkpt->co[i][ptCNS]==number_novalue)
                par->chkpt->co[i][ptCNS]= P->co[r][c];

            if ((long)par->chkpt->co[i][ptCWE]==number_novalue)
                par->chkpt->co[i][ptCWE]= R->co[r][c];

            if ((long)par->chkpt->co[i][ptCNwSe]==number_novalue)
                par->chkpt->co[i][ptCNwSe]=S->co[r][c];

            if ((long)par->chkpt->co[i][ptCNeSw]==number_novalue)
                par->chkpt->co[i][ptCNeSw]=T->co[r][c];
        }
        free_doublematrix(P);
        free_doublematrix(R);
        free_doublematrix(S);
        free_doublematrix(T);
    }

    // ---------------------- (h) No value check ----------------------
    /** Assign default value if NO value is inserted */
    for (i=1; i<=par->chkpt->nrh; i++)
    {
        if ((long)par->chkpt->co[i][ptZ]==number_novalue)
            par->chkpt->co[i][ptZ] = 0.0;

        if ((long)par->chkpt->co[i][ptLC]==number_novalue)
            par->chkpt->co[i][ptLC] = 1.0;

        if ((long)par->chkpt->co[i][ptSY]==number_novalue)
            par->chkpt->co[i][ptSY] = (double)par->soil_type_land_default;

        if ((long)par->chkpt->co[i][ptS]==number_novalue)
            par->chkpt->co[i][ptS] = 0.0;

        if ((long)par->chkpt->co[i][ptA]==number_novalue)
            par->chkpt->co[i][ptA] = 0.0;

        if ((long)par->chkpt->co[i][ptSKY]==number_novalue)
            par->chkpt->co[i][ptSKY] = 1.0;

        if ((long)par->chkpt->co[i][ptCNS]==number_novalue)
            par->chkpt->co[i][ptCNS] = 0.0;

        if ((long)par->chkpt->co[i][ptCWE]==number_novalue)
            par->chkpt->co[i][ptCWE] = 0.0;

        if ((long)par->chkpt->co[i][ptCNwSe]==number_novalue)
            par->chkpt->co[i][ptCNwSe]=0.0;

        if ((long)par->chkpt->co[i][ptCNeSw]==number_novalue)
            par->chkpt->co[i][ptCNeSw]=0.0;

        if ((long)par->chkpt->co[i][ptDrDEPTH]==number_novalue)
            par->chkpt->co[i][ptDrDEPTH]=par->DepthFreeSurface; // [mm]

        if ((long)par->chkpt->co[i][ptMAXSWE]==number_novalue)
            par->chkpt->co[i][ptMAXSWE]=1.E10; // [mm]

        if ((long)par->chkpt->co[i][ptLAT]==number_novalue)
            par->chkpt->co[i][ptLAT]= par->latitude;

        if ((long)par->chkpt->co[i][ptLON]==number_novalue)
            par->chkpt->co[i][ptLON] = par->longitude;

        if ((long)par->chkpt->co[i][ptID]==number_novalue)
            par->chkpt->co[i][ptID] = (double)i;

        if ((long)par->chkpt->co[i][ptHOR]==number_novalue)
            par->chkpt->co[i][ptHOR] = par->chkpt->co[i][ptID];
    }

    // ---------------------- (i) Show results ----------------------
    fprintf(flog,"\nPOINTS:\n");
    fprintf(flog, "ID,East[m],North[m],Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],\
SkyViewFactor[-],CurvatureN-S[1/m],CurvatureW-E[1/m],CurvatureNW-SE[1/m],CurvatureNE-SW[1/m],\
DepthFreeSurface[mm],Hor,maxSWE[mm],Lat[deg],Long[deg]\n");

    for (r=1; r<=par->chkpt->nrh; r++)
    {
        for (c=1; c<=ptTOT; c++)
        {
            fprintf(flog,"%f",par->chkpt->co[r][c]);
            if (c<ptTOT) fprintf(flog, ",");
        }
        fprintf(flog,"\n");
    }

    // ---------------------- (l) Set UV ----------------------
    if (read_dem==0 && read_lu==0 && read_soil==0 && read_sl==0 && read_as==0 && read_sk==0)
    {
        UV->U.reset(new Vector<double>{4});
        UV->V.reset(new Vector<double>{2});
    }
    UV->U->co[2]=1.0;
    UV->U->co[1]=1.0;
    UV->U->co[4]=0.0;
    UV->U->co[3]=0.0;
    UV->V->co[2]=(double)number_novalue;
    if (UV->V->co[2]<0)
    {
        UV->V->co[1] = -1.;
    }
    else
    {
        UV->V->co[1] = 1.;
    }

    // ---------------------- (n) Deallocation ----------------------
    if (read_dem==1)
        free_doublematrix(Z);
    if (read_lu==1)
        free_doublematrix(LU);
    // if(par->recover==0 && read_dem==1){
    if (read_dem==1)
    {
    }

    // 5. SET CHECKPOINT
    if (par->state_pixel == 1)
    {
        par->rc=new_longmatrix(par->chkpt->nrh,2);
        for (i=1; i<=par->chkpt->nrh; i++)
        {
            par->rc->co[i][1]=1;
            par->rc->co[i][2]=i;
        }
    }

    // 6. SET PROPERTIES
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
    top->BC_DepthFreeSurface.reset(new Vector<double>{par->chkpt->nrh});
    par->maxSWE=new_doublematrix(1,par->chkpt->nrh);
    top->horizon_point=new_longmatrix(1,par->chkpt->nrh);
    top->dzdE=new_doublematrix(1,par->chkpt->nrh);
    top->dzdN=new_doublematrix(1,par->chkpt->nrh);
    top->latitude=new_doublematrix(1,par->chkpt->nrh);
    top->longitude=new_doublematrix(1,par->chkpt->nrh);
    par->IDpoint.reset(new Vector<long>{par->chkpt->nrh});
    IT->bed=new_doublematrix(1,par->chkpt->nrh);

    for (i=1; i<=par->chkpt->nrh; i++)
    {
        top->East->co[1][i]=par->chkpt->co[i][ptX];
        top->North->co[1][i]=par->chkpt->co[i][ptY];
        top->Z0->co[1][i]=par->chkpt->co[i][ptZ];
        land->LC->co[1][i]=par->chkpt->co[i][ptLC];

        if ((long)land->LC->co[1][i] <= 0)
        {
            f = fopen(FailedRunFile, "w");
            fprintf(f, "Error: Point %ld has land cover type <= 0. This is not admitted.\n",i);
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (15).");
        }

        sl->type->co[1][i]=(long)par->chkpt->co[i][ptSY];

        if (sl->type->co[1][i] <= 0)
        {
            f = fopen(FailedRunFile, "w");
            fprintf(f, "Error: Point %ld has soil type <= 0. This is not admitted.\n",i);
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

        if (sl->type->co[1][i] <= 0)
        {
            f = fopen(FailedRunFile, "w");
            fprintf(f, "Error: Point %ld has horizon type <= 0. This is not admitted.\n", i);
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (17).");
        }

        par->maxSWE->co[1][i]=par->chkpt->co[i][ptMAXSWE];
        top->latitude->co[1][i]=par->chkpt->co[i][ptLAT];
        top->longitude->co[1][i]=par->chkpt->co[i][ptLON];
        (*par->IDpoint)(i)=(long)par->chkpt->co[i][ptID];

        IT->bed->co[1][i]=par->chkpt->co[i][ptBED];
        if ( (long)IT->bed->co[1][i] == number_novalue )
            IT->bed->co[1][i] = 1.E99;
    }

    // 7. SET PAR
    for (i=1; i<=par->init_date->nh; i++)
    {
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

    // 8. READ HORIZONS
    // find max top->horizon_point
    top->num_horizon_point=0;
    for (r=1; r<=top->horizon_point->nrh; r++)
    {
        for (c=1; c<=top->horizon_point->nch; c++)
        {
            if (top->horizon_point->co[r][c] > top->num_horizon_point)
                top->num_horizon_point = top->horizon_point->co[r][c];
        }
    }
    top->horizon_height=(double ***)malloc(top->num_horizon_point*sizeof(double **));
    top->horizon_numlines=(long *)malloc(top->num_horizon_point*sizeof(long));
    for (i=1; i<=top->num_horizon_point; i++)
    {
        c=0;
        do
        {
            flag = 0;
            if (c < par->chkpt->nrh)
            {
                if (top->horizon_point->co[1][c+1] != i)
                    c++;
            }
            if (c < par->chkpt->nrh)
            {
                if (top->horizon_point->co[1][c+1] != i)
                    flag=1;
            }

        }
        while (flag == 1 && c < par->chkpt->nrh);

        if (c < par->chkpt->nrh)
        {
            top->horizon_height[i-1] = read_horizon(0, i, files[fhorpoint],
                                                    IT->horizon_col_names, &num_lines, flog);
            top->horizon_numlines[i-1] = num_lines;
        }
        else
        {
            top->horizon_height[i-1] = (double **)malloc(sizeof(double *));
            top->horizon_height[i-1][0] = (double *)malloc(2.*sizeof(double));
            top->horizon_height[i-1][0][0] = 0.;
            top->horizon_height[i-1][0][1] = 0.;
            top->horizon_numlines[i-1] = 1;
        }
    }
}

//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************
//***************************************************************************************************

void set_bedrock(INIT_TOOLS *IT, SOIL *sl, CHANNEL *cnet, PAR *par, TOPO *top,
                 DOUBLEMATRIX *LC, FILE *flog)
{
    GEOLOG_PREFIX(__func__);
    DOUBLETENSOR *T;
    std::unique_ptr<Vector<double>> WT;
    long i, j, l, r, c, sy, synew;
    double zlim, z;
    short yes=0;
    FILE *f;

    // check if bedrock depth is above soil lower border, otherwise we do not need to calculate anything
    z = 0.;
    for (l=1; l<=Nl; l++)
    {
        z += sl->pa->co[1][jdz][l];
    }
    for (i=1; i<=par->total_pixel; i++)
    {
        r = top->rc_cont->co[i][1];
        c = top->rc_cont->co[i][2];
        if (IT->bed->co[r][c] < z) yes = 1;
    }

    if (yes == 1)
    {
        // consistency check
        if (IT->init_water_table_depth->nh != sl->pa->ndh)
        {
            f = fopen(FailedRunFile, "w");
            fprintf(f, "Error: Error in bedrock calculations");
            fclose(f);
            t_error("Fatal Error! Geotop is closed. See failing report (19).");
        }

        // rewrite soil type
        T=new_doubletensor(sl->pa->ndh, nsoilprop, Nl);
        for (i=1; i<=sl->pa->ndh; i++)
        {
            for (j=1; j<=nsoilprop; j++)
            {
                for (l=1; l<=Nl; l++)
                {
                    T->co[i][j][l]=sl->pa->co[i][j][l];
                }
            }
        }
        free_doubletensor(sl->pa);
        sl->pa=new_doubletensor(par->total_pixel+par->total_channel, nsoilprop, Nl);

        // rewrite initial water table depth
        WT.reset(new Vector<double>{IT->init_water_table_depth->nh});
        for (i=1; i<=IT->init_water_table_depth->nh; i++)
        {
            WT->co[i]=IT->init_water_table_depth->co[i];
        }
        IT->init_water_table_depth.reset(new Vector<double>{par->total_pixel
                                                            +par->total_channel});

        // assign jdz (is needed later)
        for (i=1; i<=sl->pa->ndh; i++)
        {
            for (l=1; l<=Nl; l++)
            {
                sl->pa->co[i][jdz][l]=T->co[1][jdz][l];
            }
        }

        for (i=1; i<=par->total_pixel+par->total_channel; i++)
        {

            if (i<=par->total_pixel)
            {
                r = top->rc_cont->co[i][1];
                c = top->rc_cont->co[i][2];
                sy = sl->type->co[r][c];
                synew = i;
                sl->type->co[r][c] = synew;
                z = 0.;
            }
            else
            {
                r = (*cnet->r)(i-par->total_pixel);
                c = (*cnet->c)(i-par->total_pixel);
                sy = (*cnet->soil_type)(i-par->total_pixel);
                synew = i;
                (*cnet->soil_type)(i-par->total_pixel) = synew;
                z = par->depr_channel;
            }

            IT->init_water_table_depth->co[synew] = WT->co[sy];

            zlim = IT->bed->co[r][c];

            for (l=1; l<=Nl; l++)
            {
                z += 0.5*sl->pa->co[synew][jdz][l];

                if (z <= zlim)
                {

                    for (j=1; j<=nsoilprop; j++)
                    {
                        sl->pa->co[synew][j][l] = T->co[sy][j][l];
                    }
                }
                else
                {
                    for (j=1; j<=nsoilprop; j++)
                    {
                        sl->pa->co[synew][j][l] = IT->pa_bed->co[sy][j][l] ;
                    }
                }
                z += 0.5*sl->pa->co[synew][jdz][l];
            }
        }
        free_doubletensor(T);
    }
}

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

DOUBLETENSOR *find_Z_of_any_layer(DOUBLEMATRIX *Zsurface, DOUBLEMATRIX *slope,
                                  DOUBLEMATRIX *LC, SOIL *sl, short point)
{

    GEOLOG_PREFIX(__func__);
    DOUBLETENSOR *Z;
    double Zaverage=0., z, cosine;
    long l, r, c, n, sy;

    if (point!=1)
    {
        Zaverage=0.;
        n=0;
        for (r=1; r<=Zsurface->nrh; r++)
        {
            for (c=1; c<=Zsurface->nch; c++)
            {
                if ((long)LC->co[r][c]!=number_novalue)
                {
                    n++;
                    Zaverage += Zsurface->co[r][c];
                }
            }
        }
        Zaverage/=(double)n;
    }


    Z=new_doubletensor0(sl->pa->nch, Zsurface->nrh, Zsurface->nch);
    initialize_doubletensor(Z, (double)number_novalue);

    for (r=1; r<=Zsurface->nrh; r++)
    {
        for (c=1; c<=Zsurface->nch; c++)
        {
            if ((long)LC->co[r][c]!=number_novalue)
            {
                cosine = cos(slope->co[r][c]*Pi/180.);
                sy=sl->type->co[r][c];

                if (point!=1)
                {
                    z=1.E3*(Zsurface->co[r][c]-Zaverage);//[mm]
                }
                else
                {
                    z=0.;
                }

                l=0;
                Z->co[l][r][c]=z;

                do
                {
                    l++;
                    z -= 0.5*sl->pa->co[sy][jdz][l]*cosine;
                    Z->co[l][r][c]=z;
                    z -= 0.5*sl->pa->co[sy][jdz][l]*cosine;
                }
                while (l<sl->pa->nch);
            }
        }
    }
    return Z;
}

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

short file_exists(short key)
{
    GEOLOG_PREFIX(__func__);
    //no keyword -> -1
    //keyword but does not exist -> 0
    //keyword and exists -> 1

    geolog << "Attempting to read '"<< keywords_char[key+nmet+nsoilprop+2] << "' in the file '" << files[key]<< "': ";

    if (strcmp(files[key], string_novalue) == 0)
    {
        geolog << "not present in file list" << std::endl;
        return (-1);
    }
    else
    {
        if (existing_file(files[key])>0)
        {
            geolog << "EXISTING in format " << existing_file(files[key]) << std::endl;
            return (1);
        }
        else
        {
            geolog << "not existing" << std::endl;
            return (0);
        }
    }
}

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

void copy_soil_state(SOIL_STATE *from, SOIL_STATE *to)
{
    GEOLOG_PREFIX(__func__);
    long l,i;
    long nl=from->T->nrh,n=from->T->nch;

    for (i=1; i<=n; i++)
    {
        to->P->co[0][i] = from->P->co[0][i];
        for (l=1; l<=nl; l++)
        {
            to->P->co[l][i] = from->P->co[l][i];
            to->T->co[l][i] = from->T->co[l][i];
            to->thi->co[l][i] = from->thi->co[l][i];
        }
    }
}

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

void initialize_veg_state(STATE_VEG *V, long n)
{
    GEOLOG_PREFIX(__func__);
    V->Tv.reset(new Vector<double>{n});
    V->wsnow.reset(new Vector<double>{n});
    V->wrain.reset(new Vector<double>{n});

}

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

void copy_veg_state(STATE_VEG *from, STATE_VEG *to)
{
    GEOLOG_PREFIX(__func__);
    long i, n=from->Tv->nh;
    for (i=1; i<=n; i++)
    {
        to->Tv->co[i] = from->Tv->co[i];
        to->wrain->co[i] = from->wrain->co[i];
        to->wsnow->co[i] = from->wsnow->co[i];
    }
}

/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

