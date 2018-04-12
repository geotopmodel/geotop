/*

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 - 31 December 2016

 Copyright (c), 2016 - GEOtop Foundation

 This file is part of GEOtop 2.1

 GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

 GEOtop 2.1 is distributed as a free software in the hope to create and support
 a community of developers and users that constructively interact. If you just
 use the code, please give feedback to GEOtop Foundation and the community. Any
 way you use the model, may be the most trivial one, is significantly helpful
 for the future development of the GEOtop model. Any feedback will be highly
 appreciated.

 */

#include "input.h"
#include "parameters.h"
#include <unistd.h>
#include "inputKeywords.h"
#include "geotop_common.h"
#include <iostream>
#include "global_logger.h"

using namespace std;

//********************************************************************************************************
//********************************************************************************************************
//********************************************************************************************************
//********************************************************************************************************

//! Subroutine which reads input data, performs  geomporphological analisys and
//! allocates data

void get_all_input(long /*argc*/,
                   char ** /*argv*/,
                   Topo *top,
                   Soil *sl,
                   Land *land,
                   Meteo *met,
                   Water *wat,
                   Channel *cnet,
                   Par *par,
                   Energy *egy,
                   Snow *snow,
                   Glacier *glac,
                   Times *times,
                   mio::IOManager &iomanager)

{
  FILE *f;
  GeoMatrix<double> M;
  InitTools *IT;

  size_t a;
  short success, /*added_JDfrom0=0,*/ added_wind_xy = 0, added_wind_dir = 0,
                                      added_cloud = 0, added_Tdew = 0,
                                      added_RH = 0, added_Pint = 0;
  long l, r, c, ist, i, j, n, sy, num_cols, num_lines, day, month, year, hour,
       minute;
  double z, th_oversat, JD, k_snowred, maxSWE, SWE, D, cosslope, **matrix;

  // these variables are used wen METEOIO if OFF, so we cast them to void to
  // silence warnings when using the METEOIO library
  double initial_date_meteo, end_date_meteo;
  long d, m, y, mi, h;
  (void)initial_date_meteo;
  (void)end_date_meteo;
  (void)d;
  (void)m;
  (void)y;
  (void)mi;
  (void)h;

  short count_file_missing =
    0;  // needed to check if some output map files are present or not
  std::vector<std::string> temp2;
  string temp;

  // silence some warnings about unused vars because they appear only inside and
  // if branch
  (void)added_wind_xy;
  (void)added_wind_dir;
  (void)added_Tdew;
  (void)added_RH;
  (void)added_Pint;
  (void)added_cloud;
  maxSWE = 1.E10;  // silence warning about possibly uninitialized variable

  IT = new InitTools();

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  if (geotop::common::Variables::WORKING_DIRECTORY
      [strlen(geotop::common::Variables::WORKING_DIRECTORY.c_str()) - 1] !=
      47)
    {
      temp = geotop::common::Variables::WORKING_DIRECTORY;
      geotop::common::Variables::WORKING_DIRECTORY = temp + "/";
    }

  // reads the parameters in __control_parameters

  temp = geotop::common::Variables::WORKING_DIRECTORY + program_name;

  std::shared_ptr<geotop::input::ConfigStore> lConfigStore =
    geotop::input::ConfigStoreSingletonFactory::getInstance();
  const std::string lFilePath(temp);
  bool lParsingRes = lConfigStore->parse(lFilePath);
  if (not lParsingRes)
    {
      lg->log("Unable to parse configuration file: " + lFilePath,
              geotop::logger::CRITICAL);
      exit(1);
    }

  success = read_inpts_par(par, land, times, sl, met, IT);

  // correct state pixel
  par->Tzrun = 0;
  par->wzrun = 0;
  par->Tzmaxrun = 0;
  par->Tzminrun = 0;
  par->wzmaxrun = 0;
  par->wzminrun = 0;
  par->dUzrun = 0;
  par->SWErun = 0;

  if (geotop::common::Variables::files[fTrun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->Tzrun = 1;
    }
  if (geotop::common::Variables::files[fwrun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->wzrun = 1;
    }
  if (geotop::common::Variables::files[fTmaxrun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->Tzmaxrun = 1;
    }
  if (geotop::common::Variables::files[fwmaxrun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->wzmaxrun = 1;
    }
  if (geotop::common::Variables::files[fTminrun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->Tzminrun = 1;
    }
  if (geotop::common::Variables::files[fwminrun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->wzminrun = 1;
    }
  if (geotop::common::Variables::files[fdUrun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->dUzrun = 1;
    }
  if (geotop::common::Variables::files[fSWErun] !=
      geotop::input::gStringNoValue)
    {
      if (par->point_sim == 1) par->state_pixel = 1;
      if (par->state_pixel == 1) par->SWErun = 1;
    }
  if (par->newperiodinit == 2 && (par->Tzrun == 0 || par->wzrun == 0))
    {
      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
      fprintf(f,
              "Error: You have to assign a name to the Tzrun and wzrun files\n");
      fclose(f);
      lg->log("You have to assign a name to the Tzrun and wzrun files",
              geotop::logger::ERROR);
      lg->log("Geotop failed. See failing report.", geotop::logger::CRITICAL);
      exit(1);
    }

  lg->logf("file index SPAR (soil parameter) %d", fspar);
  lg->logf("file name assigned : %s \n",
           geotop::common::Variables::files[fspar].c_str());

  success = read_soil_parameters(geotop::common::Variables::files[fspar], IT,
                                 sl, par->soil_type_bedr_default);

  geotop::common::Variables::Nl = sl->pa.getCh() - 1;

  success = read_point_file(geotop::common::Variables::files[fpointlist],
                            IT->point_col_names, par);

  geotop::common::Variables::max_time = 0.;
  geotop::common::Variables::max_time =
    (par->end_date - par->init_date) * 86400.;  // seconds

  // Recovering
  par->delay_day_recover = 0.0;
  par->n_ContRecovery = 0;

  temp = geotop::common::Variables::SuccessfulRunFile + ".old";
  rename(geotop::common::Variables::SuccessfulRunFile.c_str(), temp.c_str());
  temp = geotop::common::Variables::FailedRunFile + ".old";
  rename(geotop::common::Variables::FailedRunFile.c_str(), temp.c_str());

  if (par->ContRecovery > 0)
    {
      if (mio::IOUtils::fileExists(
            string(geotop::common::Variables::files[rtime]) + string(textfile)))
        {
          temp = geotop::common::Variables::files[rtime] + string(textfile);
          matrix = read_txt_matrix_2(temp, 33, 44, 7, &num_lines);
          par->delay_day_recover = matrix[0][1];
          par->n_ContRecovery = (long)matrix[0][2];
          geotop::common::Variables::i_run0 = (long)matrix[0][3];
          geotop::common::Variables::i_sim0 = (long)matrix[0][4];
          geotop::common::Variables::cum_time = matrix[0][5];
          geotop::common::Variables::elapsed_time_start = matrix[0][6];
          for (i = 0; i < num_lines; i++)
            {
              free(matrix[i]);
            }
          free(matrix);
        }
    }

  // Time indices

  par->init_date += par->delay_day_recover;
  convert_JDfrom0_JDandYear(par->init_date, &JD, &year);
  convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);
  // TODO remove this, after having understood line above which contains i_run0

  geotop::common::Variables::i_run =
    geotop::common::Variables::i_run0;  // Run index

  /****************************************************************************************************/
  /*! Reading of the Input files: */
  /****************************************************************************************************/
  par->cum_prec = 0;
  par->cum_da_up = 0;
  par->time_wo_prec = 0;
  par->evento = 0;     // false
  par->up_albedo = 0;  // false
  par->tres_wo_prec = 12 * 3600;

  // ##########################################################################################################################
  // #######################################################################################################################
  par->use_meteoio_cloud = false;
#ifndef USE_INTERNAL_METEODISTR
  par->use_meteoio_cloud = true;
#endif

#ifdef ILWR_PRESENT
  par->use_ilwr_wrf = true;  // TODO: convert to cmake flag
#else
  par->use_ilwr_wrf = false;
#endif
  meteoio_init(iomanager);
  // ##################################################################################################################################
  // ##################################################################################################################################

  if (par->point_sim != 1)    // distributed simulation
    {
      read_inputmaps(top, land, sl, par, IT, iomanager);
    }
  else
    {
      read_optionsfile_point(par, top, land, sl, times, IT);
    }

  geotop::common::Variables::Nr = top->Z0.getRows() - 1;
  geotop::common::Variables::Nc = top->Z0.getCols() - 1;

  par->total_pixel = 0;
  par->total_area = 0.;
  for (r = 1; r <= geotop::common::Variables::Nr; r++)
    {
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
            {
              par->total_pixel++;
              if (par->point_sim != 1)
                {
                  par->total_area += (geotop::common::Variables::UV->U[1] *
                                      geotop::common::Variables::UV->U[2]) /
                                     cos(top->slope[r][c] * GTConst::Pi / 180.);
                }
              else
                {
                  par->total_area += (geotop::common::Variables::UV->U[1] *
                                      geotop::common::Variables::UV->U[2]);
                }
            }
        }
    }

  top->Z =
    find_Z_of_any_layer(top->Z0, top->slope, land->LC, sl, par->point_sim);

  top->Jdown.resize(par->total_pixel + 1, 4 + 1);
  top->Qdown.resize(par->total_pixel + 1, 4 + 1);
  for (i = 1; i <= long(par->total_pixel); i++)
    {
      for (j = 1; j <= 4; j++)
        {
          top->Jdown[i][j] = i;
          top->Qdown[i][j] = 0.;
        }
    }

  lg->logf("INFO about the simulated basin:");
  lg->logf("Valid pixels: %ld", par->total_pixel);
  lg->logf("Number of nodes: %ld",
           (geotop::common::Variables::Nl + 1) * par->total_pixel);
  lg->logf("Novalue pixels: %ld",
           geotop::common::Variables::Nr * geotop::common::Variables::Nc -
           par->total_pixel);

  lg->logf("Basin area: %12g km2",
           (double)par->total_pixel * geotop::common::Variables::UV->U[1] *
           geotop::common::Variables::UV->U[2] / 1.E6);

  /****************************************************************************************************/
  // Reading METEO data file
  lg->log("Start Reading METEO data files...");

  num_cols = (long)nmet;  // number of possible columns of meteodata

  //  meteo data
  met->data = (double ** *)malloc(met->st->E.size() * sizeof(double **));
  //  number of line of meteo data
  met->numlines = (long *)malloc(met->st->E.size() * sizeof(long));

  success = read_meteostations_file(
              met->imeteo_stations, met->st, geotop::common::Variables::files[fmetstlist],
              IT->meteostations_col_names);
  if (success == 0)
    {
      lg->logf(" File for keyword %s is not assigned",
               geotop::common::Variables::filenames[fmetstlist].c_str());
    }
  mio::Config cfg = iomanager.getConfig();
  std::string input_meteo_plugin = cfg.get("METEO", "Input");
  if (input_meteo_plugin != "GEOTOP")
    {
      success = fill_GTmeteostations_meta(par->init_date, iomanager, met);
    }

  //  horizon for meteo stations
  met->horizon = (double ** *)malloc(met->st->E.size() * sizeof(double **));
  //  number of line in the horizon file
  met->horizonlines = (long *)malloc(met->st->E.size() * sizeof(long));
  //  line of met->data used (stored in memory to avoid from searching from
  //the first line)

#ifdef USE_INTERNAL_METEODISTR
  met->var = (double **)malloc((met->st->E.size() - 1) * sizeof(double *));
  met->line_interp_WEB = (long *)malloc(met->st->E.size() * sizeof(long));
  met->line_interp_Bsnow = (long *)malloc(met->st->E.size() * sizeof(long));
  met->line_interp_WEB_LR = 0;
  met->line_interp_Bsnow_LR = 0;
#endif
  long num_met_stat = met->st->E.size() - 1;
  lg->logf("Numbers of Meteo Stations read: %d", num_met_stat);

  if (num_met_stat < 0) num_met_stat = 0;

  // loop over meteo files..

  for (size_t ii = 1; ii <= (size_t)num_met_stat; ii++)
    {
      // check if there is a list of meteostation
      if (met->imeteo_stations[1] != geotop::input::gDoubleNoValue)
        {
          ist = met->imeteo_stations[ii];
        }
      else
        {
          ist = ii;
        }

      lg->logf("Reading horizon for : %d meteo station", ii);
      //  read horizon
      met->horizon[ii - 1] =
        read_horizon(1, ist, geotop::common::Variables::files[fhormet],
                     IT->horizon_col_names, &num_lines);
      met->horizonlines[ii - 1] = num_lines;

      // ### Probably the next lines should not be necessary if we use meteoIO ###

#ifdef USE_INTERNAL_METEODISTR
      // initialize
      met->line_interp_WEB[ii - 1] = 0;
      met->line_interp_Bsnow[ii - 1] = 0;

      // allocate var
      met->var[ii - 1] = (double *)malloc(num_cols * sizeof(double));

      // filename
      if (geotop::common::Variables::files[fmet] !=
          geotop::input::gStringNoValue)
        {
          // read matrix
          temp = namefile_i(geotop::common::Variables::files[fmet], ist);

          met->data[ii - 1] =
            read_txt_matrix(temp, 33, 44, IT->met_col_names, nmet, &num_lines);

          if ((long)met->data[ii - 1][0][iDate12] == geotop::input::gDoubleAbsent &&
              (long)met->data[ii - 1][0][iJDfrom0] ==
              geotop::input::gDoubleAbsent)
            {
              f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
              fprintf(f, "Error:: Date Column missing in file %s\n", temp.c_str());
              fclose(f);

              lg->log("Date Column missing in file " + temp, geotop::logger::ERROR);
              lg->log("Geotop failed. See failing report (2).",
                      geotop::logger::CRITICAL);
              exit(1);
            }
          met->numlines[ii - 1] = num_lines;
          lg->logf(
            "Read meteo data for meteo station #: %d with number of lines=%d", ii,
            num_lines);

          // get initial and final dates and compares against simulation periods
          // let us make comparison in JDfrom0 format.

          initial_date_meteo =
            convert_dateeur12_JDfrom0(met->data[ii - 1][1][iDate12]);
          end_date_meteo =
            convert_dateeur12_JDfrom0(met->data[ii - 1][num_lines - 1][iDate12]);

          convert_dateeur12_daymonthyearhourmin(met->data[ii - 1][0][iDate12], &d,
                                                &m, &y, &h, &mi);
          lg->logf(
            "Initial date of meteo data : %02.0f/%02.0f/%04.0f %02.0f:%02.0f",
            (float)d, (float)m, (float)y, (float)h, (float)mi);

          convert_dateeur12_daymonthyearhourmin(
            met->data[ii - 1][num_lines - 1][iDate12], &d, &m, &y, &h, &mi);

          lg->logf(
            "Final date of meteo data   : %02.0f/%02.0f/%04.0f %02.0f:%02.0f",
            (float)d, (float)m, (float)y, (float)h, (float)mi);

          if (initial_date_meteo > par->init_date)
            {
              lg->log(
                "Initial date for meteo data greater that starting date of "
                "simulation",
                geotop::logger::WARNING);
            }
          if (end_date_meteo < par->end_date)
            {
              lg->log(
                "Ending date for meteo data smaller than ending date of simulation",
                geotop::logger::WARNING);
            }

          // fixing dates: converting times in the same standard time set for the
          // simulation and fill JDfrom0
          short added_JDfrom0 =
            fixing_dates(ist, met->data[ii - 1], par->ST, met->st->ST[ii],
                         met->numlines[ii - 1], iDate12, iJDfrom0);

          switch (added_JDfrom0)
            {
            case 1:
              // Conversion performed
              break;
            case 0:
              // No conversion has been performed
              break;
            case -1:
              // An error occurred
              exit(1);
              break;
            default:
              // An unknown error occurred
              exit(666);
              break;
            }

          check_times(ist, met->data[ii - 1], met->numlines[ii - 1], iJDfrom0);

          // find clouds
          if (IT->met_col_names[itauC] != geotop::input::gStringNoValue)
            {
              if ((long)met->data[ii - 1][0][itauC] == geotop::input::gDoubleAbsent ||
                  par->ric_cloud == 1)
                {
                  added_cloud = fill_meteo_data_with_cloudiness(
                                  met->data[ii - 1], met->numlines[ii - 1], met->horizon[ii - 1],
                                  met->horizonlines[ii - 1], met->st->lat[ii], met->st->lon[ii],
                                  par->ST, met->st->Z[ii], met->st->sky[ii], 0.0, par->ndivdaycloud,
                                  par->dem_rotation, par->Lozone, par->alpha_iqbal, par->beta_iqbal,
                                  0.);
                }
            }

          // calculate Wx and Wy if Wspeed and direction are given
          if (par->wind_as_xy == 1)
            {
              added_wind_xy = fill_wind_xy(
                                met->data[ii - 1], met->numlines[ii - 1], iWs, iWdir, iWsx, iWsy,
                                IT->met_col_names[iWsx], IT->met_col_names[iWsy]);
            }

          // calculate Wspeed and direction if Wx and Wy are given
          if (par->wind_as_dir == 1)
            {
              added_wind_dir = fill_wind_dir(
                                 met->data[ii - 1], met->numlines[ii - 1], iWs, iWdir, iWsx, iWsy,
                                 IT->met_col_names[iWs], IT->met_col_names[iWdir]);
            }

          // find Tdew
          if (par->vap_as_Td == 1)
            {
              added_Tdew =
                fill_Tdew(ii, met->st->Z, met->data[ii - 1], met->numlines[ii - 1],
                          iRh, iT, iTdew, IT->met_col_names[iTdew], par->RHmin);
            }

          // find RH
          if (par->vap_as_RH == 1)
            {
              added_RH =
                fill_RH(ii, met->st->Z, met->data[ii - 1], met->numlines[ii - 1], iRh,
                        iT, iTdew, IT->met_col_names[iRh]);
            }

          // find Prec Intensity
          if (par->linear_interpolation_meteo[ii] == 1 &&
              (long)met->data[ii - 1][0][iPrec] != geotop::input::gDoubleAbsent)
            {
              f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
              fprintf(f,
                      "Meteo data for station %ld contain precipitation as volume, "
                      "but Linear Interpolation is set. This is not possible, the "
                      "precipitation data are removed.\n",
                      ii);
              fprintf(f,
                      "If you want to use precipitation as volume, you cannot set "
                      "keyword LinearInterpolation at 1.\n");
              fclose(f);
              lg->logsf(geotop::logger::ERROR,
                        "Meteo data for station %ld contain precipitation as volume, "
                        "but Linear Interpolation is set. This is not possible, the "
                        "precipitation data are removed.",
                        ii);
              lg->log("Geotop failed. See failing report (3).",
                      geotop::logger::CRITICAL);
              exit(1);
            }

          if (par->prec_as_intensity == 1)
            {
              added_Pint = fill_Pint(met->data[ii - 1], met->numlines[ii - 1], iPrec,
                                     iPrecInt, iJDfrom0, IT->met_col_names[iPrecInt]);
            }

          // calculate Wx and Wy if Wspeed and direction are given
          if (par->wind_as_xy != 1)
            {
              added_wind_xy = fill_wind_xy(
                                met->data[ii - 1], met->numlines[ii - 1], iWs, iWdir, iWsx, iWsy,
                                IT->met_col_names[iWsx], IT->met_col_names[iWsy]);
            }

          // find Prec Intensity
          if (par->prec_as_intensity != 1)
            {
              added_Pint = fill_Pint(met->data[ii - 1], met->numlines[ii - 1], iPrec,
                                     iPrecInt, iJDfrom0, IT->met_col_names[iPrecInt]);
            }

          // find Tdew
          if (par->vap_as_Td != 1)
            {
              added_Tdew =
                fill_Tdew(ii, met->st->Z, met->data[ii - 1], met->numlines[ii - 1],
                          iRh, iT, iTdew, IT->met_col_names[iTdew], par->RHmin);
            }

          // free(temp);

        }
      else
        {
          lg->log(
            "File meteo not in the list, meteo data not read, used default values",
            geotop::logger::WARNING);

          met->data[ii - 1] = (double **)malloc(2 * sizeof(double *));

          for (n = 1; n <= 2; n++)
            {
              met->data[ii - 1][n - 1] = (double *)malloc(num_cols * sizeof(double));
              for (j = 1; j <= nmet; j++)
                {
                  met->data[ii - 1][n - 1][j - 1] =
                    (double)geotop::input::gDoubleAbsent;
                }
            }

          met->data[ii - 1][0][iJDfrom0] = 0.;
          met->data[ii - 1][1][iJDfrom0] = 1.E10;
        }
#endif  // USE_INTERNAL_METEODISTR
    }

  // read LAPSE RATES FILE

  if (geotop::common::Variables::files[fLRs] !=
      geotop::input::gStringNoValue)    // s stands for string
    {

      if (!mio::IOUtils::fileExists(
            string(geotop::common::Variables::files[fLRs]) + string(textfile)))
        {
          lg->log(
            "Lapse rate file unavailable. Check input files. If you do not have a "
            "lapse rate file, remove its name and keywords from input file",
            geotop::logger::WARNING);
        }
      temp = geotop::common::Variables::files[fLRs] + string(textfile);
      met->LRs = read_txt_matrix(temp, 33, 44, IT->lapserates_col_names, nlstot,
                                 &num_lines);
      met->LRsnr = num_lines;
      par->LRflag = 1;
      lg->log("Lapse rate file read");

    }
  else
    {
      par->LRflag = 0;
    }

  n = (long)nlstot;
  met->LRv = (double *)malloc(n * sizeof(double));
  met->LRd = (double *)malloc(n * sizeof(double));
  for (i = 0; i < nlstot; i++)
    {
      met->LRv[i] = geotop::input::gDoubleNoValue;
      if (i == ilsTa)
        {
          met->LRd[i] = GTConst::LapseRateTair;
        }
      else if (i == ilsTdew)
        {
          met->LRd[i] = GTConst::LapseRateTdew;
        }
      else if (i == ilsPrec)
        {
          met->LRd[i] = GTConst::LapseRatePrec;
        }
      else
        {
          met->LRd[i] = 0.0;
        }
    }

#ifdef USE_INTERNAL_METEODISTR

  // Find the first station with shortwave radiation data
  met->nstsrad = 0;

  do
    {
      met->nstsrad++;
      a = 0;
      if ((long)met->data[met->nstsrad - 1][0][iSW] !=
          geotop::input::gDoubleAbsent ||
          ((long)met->data[met->nstsrad - 1][0][iSWb] !=
           geotop::input::gDoubleAbsent &&
           (long)met->data[met->nstsrad - 1][0][iSWd] !=
           geotop::input::gDoubleAbsent))
        {
          a = 1;
        }
    }
  while (met->nstsrad < num_met_stat && a == 0);

  if (a == 0)
    {
      lg->log("NO shortwave radiation measurements available",
              geotop::logger::WARNING);
    }
  else
    {
      lg->logf("Shortwave radiation measurements from station %ld\n",
               met->nstsrad);
    }

  // Find the first station with cloud data
  met->nstcloud = 0;
  do
    {
      met->nstcloud++;
      a = 0;
      if ((long)met->data[met->nstcloud - 1][0][iC] !=
          geotop::input::gDoubleAbsent ||
          (long)met->data[met->nstcloud - 1][0][itauC] !=
          geotop::input::gDoubleAbsent)
        {
          a = 1;
        }
    }
  while ((size_t)(met->nstcloud) < met->st->Z.size() - 1 && a == 0);

  if (a == 0)
    {
      lg->log("No cloudiness measurements available", geotop::logger::WARNING);
    }
  else
    {
      lg->logf("Cloudiness measurements from station %ld", met->nstcloud);
    }

  // Find the first station with longwave radiation data
  met->nstlrad = 0;
  do
    {
      met->nstlrad++;
      a = 0;
      if ((long)met->data[met->nstlrad - 1][0][iLWi] !=
          geotop::input::gDoubleAbsent)
        a = 1;
    }
  while ((size_t)(met->nstlrad) < met->st->Z.size() - 1 && a == 0);

  if (a == 0)
    {
      lg->log("No longwave radiation measurements available",
              geotop::logger::WARNING);
    }
  else
    {
      lg->logf("Longwave radiation measurements from station %ld\n",
               met->nstlrad);
    }

#endif  // USE_INTERNAL_METEODISTR

  met->tau_cloud.resize(top->Z0.getRows(), top->Z0.getCols(),
                        geotop::input::gDoubleNoValue);

  met->tau_cloud_av.resize(top->Z0.getRows(), top->Z0.getCols(),
                           geotop::input::gDoubleNoValue);

  met->tau_cloud_yes.resize(top->Z0.getRows(), top->Z0.getCols(),
                            (short)geotop::input::gDoubleNoValue);

  met->tau_cloud_av_yes.resize(top->Z0.getRows(), top->Z0.getCols(),
                               (short)geotop::input::gDoubleNoValue);

  //  vector defining which meteo station has the SW radiation information

  met->st->flag_SW_meteoST.resize(met->st->Z.size(),
                                  geotop::input::gDoubleNoValue);

  met->st->tau_cloud_av_yes_meteoST.resize(met->st->Z.size(),
                                           geotop::input::gDoubleNoValue);

  met->st->tau_cloud_yes_meteoST.resize(met->st->Z.size(),
                                        geotop::input::gDoubleNoValue);

  met->st->tau_cloud_av_meteoST.resize(met->st->Z.size(),
                                       geotop::input::gDoubleNoValue);

  met->st->tau_cloud_meteoST.resize(met->st->Z.size(),
                                    geotop::input::gDoubleNoValue);

  // i.show details on checkpoints
  lg->writeAll("\nCHECKPOINTS:");
  lg->writeAll(
    "ID,r,c,Elevation[masl],LandCoverType,SoilType,Slope[deg],Aspect[deg],"
    "SkyViewFactor[-]\n");

  for (i = 1; i < long(par->rc.getRows()); i++)
    {
      r = par->rc[i][1];
      c = par->rc[i][2];
      lg->writefAll("%ld,%ld,%ld,%12g,%d,%ld,%12g,%12g,%12g\n", i, r, c,
                    top->Z0[r][c], (short)land->LC[r][c], sl->type[r][c],
                    top->slope[r][c], top->aspect[r][c], top->sky[r][c]);
    }

  // i.show meteo stations
  lg->writeAll("\nMETEO STATIONS:\n");
  lg->writeAll(
    "ID,East[m],North[m],Lat[deg],Lon[deg],Elev[m "
    "a.s.l.],Sky[deg],Stand_Time[h],WindSensHeight[m],TempSensHeight[m]\n");

  for (size_t rr = 1; rr < met->st->E.size(); rr++)
    {
      lg->writefAll("%ld,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g,%12g\n", rr,
                    met->st->E[rr], met->st->N[rr], met->st->lat[rr],
                    met->st->lon[rr], met->st->Z[rr], met->st->sky[rr],
                    met->st->ST[rr], met->st->Vheight[rr], met->st->Theight[rr]);
    }

  /****************************************************************************************************/
  // read INCOMING DISCHARGE

  met->qinv = (double *)malloc(2 * sizeof(double));
  met->qinv[0] = 0.;
  met->qinv[1] = 0.;

  if (geotop::common::Variables::files[fqin] != geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fqin] + string(textfile);
      temp2.push_back("Time");
      temp2.push_back("Qx");
      met->qins = read_txt_matrix(temp, 33, 44, temp2, 2, &num_lines);
      met->qinsnr = num_lines;
      par->qin = 1;
      met->qinline = 0;
      lg->log("Incoming discharge file read");
    }
  else
    {
      par->qin = 0;
    }

  /****************************************************************************************************/
  /*! Completing the several time-indipendent input variables with the data of
   * input-files:           */
  /****************************************************************************************************/
  /****************************************************************************************************/
  // Completing of "land" (of the type LAND):

  //  Initialize matrix of shadow

  land->shadow.resize(geotop::common::Variables::Nr + 1,
                      geotop::common::Variables::Nc + 1, 0);

  //  Check that there aren't cell with an undefined land use value
  z = 0.;
  l = 0;
  do
    {
      l++;
      z += sl->pa(1, jdz, l);
    }
  while (l < geotop::common::Variables::Nl && z < GTConst::z_transp);

  land->root_fraction.resize(par->n_landuses + 1, l + 1, 0.0);

  // check vegetation variable consistency
  if ((jHveg != jdHveg + jHveg - 1) ||
      (jz0thresveg != jdz0thresveg + jHveg - 1) ||
      (jz0thresveg2 != jdz0thresveg2 + jHveg - 1) ||
      (jLSAI != jdLSAI + jHveg - 1) || (jcf != jdcf + jHveg - 1) ||
      (jdecay0 != jddecay0 + jHveg - 1) || (jexpveg != jdexpveg + jHveg - 1) ||
      (jroot != jdroot + jHveg - 1) || (jrs != jdrs + jHveg - 1))
    {
      lg->log("Vegetation variables not consistent", geotop::logger::CRITICAL);
      exit(1);
    }

  //  variables used to assign vegetation properties that change with time
  num_cols = jdvegprop + 1;
  land->vegpars = (double ** *)malloc(par->n_landuses * sizeof(double **));
  land->vegparv = (double **)malloc(par->n_landuses * sizeof(double *));
  land->NumlinesVegTimeDepData = (long *)malloc(par->n_landuses * sizeof(long));

  land->vegpar.resize(jdvegprop + 1);

  par->vegflag.resize(par->n_landuses + 1, 0);

  for (i = 1; (size_t)i <= par->n_landuses; i++)
    {
      if (geotop::common::Variables::files[fvegpar] !=
          geotop::input::gStringNoValue)    // s stands for string
        {

          temp = namefile_i_we2(geotop::common::Variables::files[fvegpar], i);

          if (mio::IOUtils::fileExists(string(temp) + string(textfile)))
            {
              lg->logf(
                "There is a specific vegetation parameter file for land cover type = "
                "%ld",
                i);
              temp = namefile_i(geotop::common::Variables::files[fvegpar], i);
              land->vegpars[i - 1] =
                read_txt_matrix_2(temp, 33, 44, num_cols, &num_lines);
              land->NumlinesVegTimeDepData[i - 1] = num_lines;
              par->vegflag[i] = 1;
            }
          else
            {
              lg->logsf(geotop::logger::WARNING,
                        "There is NOT a specific vegetation parameter file for land "
                        "cover type = %ld",
                        i);
            }
        }

      land->vegparv[i - 1] = (double *)malloc(num_cols * sizeof(double));
      for (j = 0; j < num_cols; j++)
        {
          land->vegparv[i - 1][j] = geotop::input::gDoubleNoValue;
        }

      //  z0 (convert in m)
      land->ty[i][jz0] *= 0.001;

      //  find root fraction
      root(land->root_fraction.getCols(), land->ty[i][jroot], 0.0, sl->pa,
           land->root_fraction, i);

      //  error messages
      for (size_t ll = 1; ll < met->st->Z.size(); ll++)
        {
          if (0.001 * land->ty[i][jHveg] > met->st->Vheight[ll] ||
              0.001 * land->ty[i][jHveg] > met->st->Theight[ll])
            {
              f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
              fprintf(f,
                      "hc:%f m, zmu:%f m, zmt:%f m - hc must be lower than "
                      "measurement height - land cover %ld, meteo station %ld\n",
                      0.001 * land->ty[i][jHveg], met->st->Vheight[ll],
                      met->st->Theight[ll], i, ll);
              fclose(f);

              lg->logsf(geotop::logger::ERROR,
                        "hc:%12g m, zmu:%12g m, zmt:%12g m - hc must be lower than "
                        "measurement height - land cover %ld, meteo station %ld",
                        0.001 * land->ty[i][jHveg], met->st->Vheight[ll],
                        met->st->Theight[ll], i, ll);
              lg->log("Geotop failed. See failing report (5).",
                      geotop::logger::CRITICAL);
              exit(1);
            }
        }
    }

  /****************************************************************************************************/
  /*! Filling up of the struct "channel" (of the type CHANNEL): */

  /*The number of channel-pixel are counted:*/
  i = 0;
  for (r = 1; r <= geotop::common::Variables::Nr; r++)
    {
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          if (top->pixel_type[r][c] >= 10) i++;
        }
    }

  lg->logf("Channel pixels: %ld", i);
  par->total_channel = i;

  // allocate channel vectors/matrixes
  if (i == 0) i = 1;

  cnet->Vout = 0.;

  cnet->r.resize(i + 1, 0);

  cnet->c.resize(i + 1, 0);

  cnet->ch.resize(geotop::common::Variables::Nr + 1,
                  geotop::common::Variables::Nc + 1, 0);

  cnet->ch_down.resize(i + 1, 0);

  cnet->length.resize(i + 1, 0.0);

  cnet->Vsup.resize(i + 1, 0.);
  cnet->Vsub.resize(i + 1, 0.);

  cnet->h_sup.resize(i + 1, 0.0);

  cnet->soil_type.resize(i + 1, par->soil_type_chan_default);

  if (par->total_channel > 1)
    enumerate_channels(cnet, land->LC, top->pixel_type, top->Z0, top->slope,
                       geotop::input::gDoubleNoValue);

  cnet->ch3 =
    (long **)malloc((geotop::common::Variables::Nl + 1) * sizeof(long *));
  for (l = 0; l <= geotop::common::Variables::Nl; l++)
    {
      cnet->ch3[l] = (long *)malloc((i + 1) * sizeof(long));
    }

  cnet->lch.resize(((geotop::common::Variables::Nl + 1) * i) + 1, 2 + 1, 0);

  lch3_cont(cnet->ch3, cnet->lch, geotop::common::Variables::Nl,
            par->total_channel);

  /****************************************************************************************************/
  // Calculates distance from the main channel

  // Cont for Richards 3D
  // There are not channels

  top->i_cont =
    (long ** *)malloc((geotop::common::Variables::Nl + 1) * sizeof(long **));
  for (l = 0; l <= geotop::common::Variables::Nl; l++)
    {
      top->i_cont[l] =
        (long **)malloc((geotop::common::Variables::Nr + 1) * sizeof(long *));
      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          top->i_cont[l][r] =
            (long *)malloc((geotop::common::Variables::Nc + 1) * sizeof(long));
        }
    }

  top->lrc_cont.resize(
    (geotop::common::Variables::Nl + 1) * par->total_pixel + 1, 3 + 1, 0);

  i_lrc_cont(land->LC, top->i_cont, top->lrc_cont);

  top->j_cont =
    (long **)malloc((geotop::common::Variables::Nr + 1) * sizeof(long *));
  for (r = 1; r <= geotop::common::Variables::Nr; r++)
    {
      top->j_cont[r] =
        (long *)malloc((geotop::common::Variables::Nc + 1) * sizeof(long));
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          top->j_cont[r][c] = 0;
        }
    }

  top->rc_cont.resize(par->total_pixel + 1, 2 + 1);

  j_rc_cont(land->LC, top->j_cont, top->rc_cont);

  if (par->state_pixel == 1)
    {
      par->jplot.resize(par->total_pixel + 1, 0);

      for (i = 1; i <= long(par->total_pixel); i++)
        {
          for (j = 1; j < long(par->rc.getRows()); j++)
            {
              if (top->rc_cont[i][1] == par->rc[j][1] &&
                  top->rc_cont[i][2] == par->rc[j][2])
                {
                  par->jplot[i] = j;
                }
            }
        }
    }

  // BEDROCK (adjusting soil properties)
  par->bedrock = 0;
  if (geotop::common::Variables::files[fbed] != geotop::input::gStringNoValue)
    set_bedrock(IT, sl, cnet, par, top, land->LC);

  /****************************************************************************************************/
  /*! Completing of the initialization of SOIL structure */
  /****************************************************************************************************/

  sl->SS = new SoilState();
  initialize_soil_state(sl->SS, par->total_pixel,
                        geotop::common::Variables::Nl);

  sl->VS = new StateVeg();
  initialize_veg_state(sl->VS, par->total_pixel);

  sl->th.resize(geotop::common::Variables::Nl + 1, par->total_pixel + 1,
                geotop::input::gDoubleNoValue);

  sl->Ptot.resize(geotop::common::Variables::Nl + 1, par->total_pixel + 1,
                  geotop::input::gDoubleNoValue);

  if (geotop::common::Variables::files[fTav] != geotop::input::gStringNoValue ||
      geotop::common::Variables::files[fTavsup] !=
      geotop::input::gStringNoValue)
    {
      sl->T_av_tensor.resize(geotop::common::Variables::Nl + 1,
                             par->total_pixel + 1, 0.0);
    }

  if (geotop::common::Variables::files[ficeav] !=
      geotop::input::gStringNoValue)
    {
      sl->thi_av_tensor.resize(geotop::common::Variables::Nl + 1,
                               par->total_pixel + 1, 0.0);
    }

  if (geotop::common::Variables::files[fliqav] !=
      geotop::input::gStringNoValue)
    {
      sl->thw_av_tensor.resize(geotop::common::Variables::Nl + 1,
                               par->total_pixel + 1, 0.0);
    }

  if (geotop::common::Variables::files[fpnet] !=
      geotop::input::gStringNoValue)    // TODO mattiu
    {
      sl->Pnetcum.resize(par->total_pixel + 1, 0.0);
    }
  if (geotop::common::Variables::files[fevap] !=
      geotop::input::gStringNoValue)
    {
      sl->ETcum.resize(par->total_pixel + 1, 0.0);
    }  // end mattiu

  sl->ET.resize(geotop::common::Variables::Nl + 1,
                geotop::common::Variables::Nr + 1,
                geotop::common::Variables::Nc + 1, 0.);

  if (!mio::IOUtils::fileExists(string(geotop::common::Variables::files[fwt0]) +
                                string(ascii_esri)))
    {
      for (i = 1; i <= long(par->total_pixel); i++)
        {
          r = top->rc_cont[i][1];
          c = top->rc_cont[i][2];

          sy = sl->type[r][c];

          if ((long)IT->init_water_table_depth[sy] !=
              geotop::input::gDoubleNoValue)
            {
              z = 0.;
              sl->SS->P[0][i] = -IT->init_water_table_depth[sy] *
                                cos(top->slope[r][c] * GTConst::Pi / 180.);

              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  z += 0.5 * sl->pa(sy, jdz, l) *
                       cos(top->slope[r][c] * GTConst::Pi / 180.);
                  sl->SS->P[l][i] = sl->SS->P[0][i] + z;
                  z += 0.5 * sl->pa(sy, jdz, l) *
                       cos(top->slope[r][c] * GTConst::Pi / 180.);
                }
            }
          else
            {
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  sl->SS->P[l][i] = sl->pa(sy, jpsi, l);
                }
            }
        }

    }
  else
    {
      meteoio_readMap(string(geotop::common::Variables::files[fwt0]), M);

      for (i = 1; i <= long(par->total_pixel); i++)
        {
          r = top->rc_cont[i][1];
          c = top->rc_cont[i][2];

          sy = sl->type[r][c];

          z = 0.;
          sl->SS->P[0][i] = -M[r][c] * cos(top->slope[r][c] * GTConst::Pi / 180.);
          for (l = 1; l <= geotop::common::Variables::Nl; l++)
            {
              z +=
                0.5 * sl->pa(sy, jdz, l) * cos(top->slope[r][c] * GTConst::Pi / 180.);
              sl->SS->P[l][i] = sl->SS->P[0][i] + z;
              z +=
                0.5 * sl->pa(sy, jdz, l) * cos(top->slope[r][c] * GTConst::Pi / 180.);
            }
        }
    }

  for (i = 1; i <= long(par->total_pixel); i++)
    {
      r = top->rc_cont[i][1];
      c = top->rc_cont[i][2];

      sy = sl->type[r][c];

      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          sl->SS->T[l][i] = sl->pa(sy, jT, l);

          sl->Ptot[l][i] = sl->SS->P[l][i];

          sl->th[l][i] = teta_psi(sl->SS->P[l][i], 0.0, sl->pa(sy, jsat, l),
                                  sl->pa(sy, jres, l), sl->pa(sy, ja, l),
                                  sl->pa(sy, jns, l), 1 - 1 / sl->pa(sy, jns, l),
                                  GTConst::PsiMin, sl->pa(sy, jss, l));

          th_oversat = Fmax(sl->SS->P[l][i], 0.0) * sl->pa(sy, jss, l);
          sl->th[l][i] -= th_oversat;

          if (sl->SS->T[l][i] <= GTConst::Tfreezing)
            {
              //  Theta_ice=Theta(without freezing) - Theta_unfrozen(in equilibrium
              //with T)
              sl->SS->thi[l][i] =
                sl->th[l][i] - teta_psi(Psif(sl->SS->T[l][i]), 0.0,
                                        sl->pa(sy, jsat, l), sl->pa(sy, jres, l),
                                        sl->pa(sy, ja, l), sl->pa(sy, jns, l),
                                        1. - 1. / sl->pa(sy, jns, l), GTConst::PsiMin,
                                        sl->pa(sy, jss, l));

              //  if Theta(without freezing)<Theta_unfrozen(in equilibrium with T)
              //Theta_ice is set at 0
              if (sl->SS->thi[l][i] < 0) sl->SS->thi[l][i] = 0.0;

              //  Psi is updated taking into account the freezing
              sl->th[l][i] -= sl->SS->thi[l][i];

              sl->SS->P[l][i] = psi_teta(
                                  sl->th[l][i] + th_oversat, sl->SS->thi[l][i], sl->pa(sy, jsat, l),
                                  sl->pa(sy, jres, l), sl->pa(sy, ja, l), sl->pa(sy, jns, l),
                                  1 - 1 / sl->pa(sy, jns, l), GTConst::PsiMin, sl->pa(sy, jss, l));
            }
        }
    }

  if (par->state_pixel == 1)
    {
      // we should insert here a check to see if all the files requested are in
      // place //

      if (geotop::common::Variables::files[fTz] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fTzwriteend] !=
          geotop::input::gStringNoValue)
        sl->Tzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl + 1);
      if (geotop::common::Variables::files[fTzav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fTzavwriteend] !=
          geotop::input::gStringNoValue)
        sl->Tzavplot.resize(par->rc.getRows(), geotop::common::Variables::Nl + 1);
      if (geotop::common::Variables::files[fpsiztot] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fpsiztotwriteend] !=
          geotop::input::gStringNoValue)
        sl->Ptotzplot.resize(par->rc.getRows(),
                             geotop::common::Variables::Nl + 1);
      if (geotop::common::Variables::files[fpsiz] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fpsizwriteend] !=
          geotop::input::gStringNoValue)
        sl->Pzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl + 1);
      if (geotop::common::Variables::files[fliqz] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fliqzwriteend] !=
          geotop::input::gStringNoValue)
        sl->thzplot.resize(par->rc.getRows(), geotop::common::Variables::Nl + 1);
      if (geotop::common::Variables::files[fliqzav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fliqzavwriteend] !=
          geotop::input::gStringNoValue)
        sl->thzavplot.resize(par->rc.getRows(),
                             geotop::common::Variables::Nl + 1);
      if (geotop::common::Variables::files[ficez] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[ficezwriteend] !=
          geotop::input::gStringNoValue)
        sl->thizplot.resize(par->rc.getRows(), geotop::common::Variables::Nl + 1);
      if (geotop::common::Variables::files[ficezav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[ficezavwriteend] !=
          geotop::input::gStringNoValue)
        sl->thizavplot.resize(par->rc.getRows(),
                              geotop::common::Variables::Nl + 1);

      for (i = 1; i < long(par->rc.getRows()); i++)
        {
          r = top->rc_cont[i][1];
          c = top->rc_cont[i][2];
          j = top->j_cont[r][c];
          for (l = 1; l <= geotop::common::Variables::Nl; l++)
            {
              if (geotop::common::Variables::files[fTz] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[fTzwriteend] !=
                  geotop::input::gStringNoValue)
                sl->Tzplot[i][l] = sl->SS->T[l][j];
              if (geotop::common::Variables::files[fTzav] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[fTzavwriteend] !=
                  geotop::input::gStringNoValue)
                sl->Tzavplot[i][l] = sl->SS->T[l][j];
              if (geotop::common::Variables::files[fliqzav] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[fliqzavwriteend] !=
                  geotop::input::gStringNoValue)
                sl->thzavplot[i][l] = sl->th[l][j];
              if (geotop::common::Variables::files[ficez] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[ficezwriteend] !=
                  geotop::input::gStringNoValue)
                sl->thizplot[i][l] = sl->SS->thi[l][j];
              if (geotop::common::Variables::files[ficezav] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[ficezavwriteend] !=
                  geotop::input::gStringNoValue)
                sl->thizavplot[i][l] = sl->SS->thi[l][j];
              if (geotop::common::Variables::files[fpsiztot] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[fpsiztotwriteend] !=
                  geotop::input::gStringNoValue)
                sl->Ptotzplot[i][l] = sl->Ptot[l][j];
            }
          for (l = 0; l <= geotop::common::Variables::Nl; l++)
            {
              if (geotop::common::Variables::files[fpsiz] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[fpsizwriteend] !=
                  geotop::input::gStringNoValue)
                sl->Pzplot[i][l] = sl->SS->P[l][j];
            }
        }
    }

  lg->logf("INPUT: printing delay_day_recover %f\n", par->delay_day_recover);

  if (par->delay_day_recover > 0)
    {
      assign_recovered_tensor_vector(1, par->recover,
                                     geotop::common::Variables::files[riceg],
                                     sl->SS->thi, top->rc_cont);
      assign_recovered_tensor_vector(1, par->recover,
                                     geotop::common::Variables::files[rTg],
                                     sl->SS->T, top->rc_cont);
      assign_recovered_tensor_vector(0, par->recover,
                                     geotop::common::Variables::files[rpsi],
                                     sl->SS->P, top->rc_cont);

      assign_recovered_map_vector(par->recover,
                                  geotop::common::Variables::files[rwcrn],
                                  sl->VS->wrain, top->rc_cont);
      assign_recovered_map_vector(par->recover,
                                  geotop::common::Variables::files[rwcsn],
                                  sl->VS->wsnow, top->rc_cont);
      assign_recovered_map_vector(par->recover,
                                  geotop::common::Variables::files[rTv],
                                  sl->VS->Tv, top->rc_cont);
    }

  //  channel soil
  cnet->SS = new SoilState();
  initialize_soil_state(cnet->SS, cnet->r.size(),
                        geotop::common::Variables::Nl);

  cnet->th.resize(geotop::common::Variables::Nl + 1, cnet->r.size());

  cnet->ET.resize(geotop::common::Variables::Nl + 1, cnet->r.size(), 0.0);

  cnet->Kbottom.resize(cnet->r.size(), 0.0);

  for (j = 1; j <= long(par->total_channel); j++)
    {
      sy = cnet->soil_type[j];
      r = cnet->r[j];
      c = cnet->c[j];

      cnet->SS->P[0][j] = sl->SS->P[0][top->j_cont[r][c]];
      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          cnet->SS->P[l][j] = sl->Ptot[l][top->j_cont[r][c]];
        }

      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          cnet->SS->T[l][j] = sl->pa(sy, jT, l);

          cnet->th[l][j] = teta_psi(
                             cnet->SS->P[l][j], 0.0, sl->pa(sy, jsat, l), sl->pa(sy, jres, l),
                             sl->pa(sy, ja, l), sl->pa(sy, jns, l), 1. - 1. / sl->pa(sy, jns, l),
                             GTConst::PsiMin, sl->pa(sy, jss, l));

          th_oversat = Fmax(cnet->SS->P[l][j], 0.0) * sl->pa(sy, jss, l);
          cnet->th[l][j] -= th_oversat;

          if (cnet->SS->T[l][j] <= GTConst::Tfreezing)
            {
              cnet->SS->thi[l][j] =
                cnet->th[l][j] - teta_psi(Psif(cnet->SS->T[l][j]), 0.0,
                                          sl->pa(sy, jsat, l), sl->pa(sy, jres, l),
                                          sl->pa(sy, ja, l), sl->pa(sy, jns, l),
                                          1. - 1. / sl->pa(sy, jns, l),
                                          GTConst::PsiMin, sl->pa(sy, jss, l));

              if (cnet->SS->thi[l][j] < 0) cnet->SS->thi[l][j] = 0.0;

              // Psi is updated taking into account the freezing
              cnet->th[l][j] -= cnet->SS->thi[l][j];
              cnet->SS->P[l][j] = psi_teta(
                                    cnet->th[l][j] + th_oversat, cnet->SS->thi[l][j], sl->pa(sy, jsat, l),
                                    sl->pa(sy, jres, l), sl->pa(sy, ja, l), sl->pa(sy, jns, l),
                                    1. - 1. / sl->pa(sy, jns, l), GTConst::PsiMin, sl->pa(sy, jss, l));
            }
        }
    }

  if (par->delay_day_recover > 0 && par->total_channel > 0)
    {
      assign_recovered_tensor_channel(0, par->recover,
                                      geotop::common::Variables::files[rpsich],
                                      cnet->SS->P, cnet->r, cnet->c);
      assign_recovered_tensor_channel(1, par->recover,
                                      geotop::common::Variables::files[ricegch],
                                      cnet->SS->thi, cnet->r, cnet->c);
      assign_recovered_tensor_channel(1, par->recover,
                                      geotop::common::Variables::files[rTgch],
                                      cnet->SS->T, cnet->r, cnet->c);

      for (i = 1; i <= long(par->total_channel); i++)
        {
          for (l = 1; l <= geotop::common::Variables::Nl; l++)
            {
              sy = cnet->soil_type[i];

              cnet->th[l][i] = teta_psi(
                                 Fmin(
                                   cnet->SS->P[l][i],
                                   psi_saturation(cnet->SS->thi[l][i], sl->pa(sy, jsat, l),
                                                  sl->pa(sy, jres, l), sl->pa(sy, ja, l),
                                                  sl->pa(sy, jns, l), 1. - 1. / sl->pa(sy, jns, l))),
                                 cnet->SS->thi[l][i], sl->pa(sy, jsat, l), sl->pa(sy, jres, l),
                                 sl->pa(sy, ja, l), sl->pa(sy, jns, l), 1. - 1. / sl->pa(sy, jns, l),
                                 GTConst::PsiMin, sl->pa(sy, jss, l));
            }
        }
    }

  //  WRITE INITIAL CONDITION
  write_output_headers(met->st->Z.size(), times, wat, par, top, land, sl, egy,
                       snow, glac);

  if (par->state_pixel == 1)
    {
      for (j = 1; j < long(par->rc.getRows()); j++)
        {
          if (par->output_vertical_distances == 1)
            {
              r = par->rc[j][1];
              c = par->rc[j][2];
              cosslope =
                cos(Fmin(GTConst::max_slope, top->slope[r][c]) * GTConst::Pi / 180.0);
            }
          else
            {
              cosslope = 1.;
            }

          write_soil_output(j, par->IDpoint[j], par->init_date, par->init_date, JD,
                            day, month, year, hour, minute, par->soil_plot_depths,
                            sl, par, GTConst::PsiMin, cosslope);
        }
    }

  //  z boundary condition
  for (l = 1; l <= geotop::common::Variables::Nl; l++)
    {
      par->Zboundary -= sl->pa(1, jdz, l);
    }

  if (par->Zboundary < 0)
    {
      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
      fprintf(f,
              "Z at which 0 annual temperature takes place is not lower than the "
              "soil column\n");
      fclose(f);
      lg->log(
        "Z at which 0 annual temperature takes place is not lower than the soil "
        "column",
        geotop::logger::ERROR);
      lg->log("Geotop failed. See failing report (6).", geotop::logger::CRITICAL);
      exit(1);
    }

  par->Zboundary *= 1.E-3;  // convert in [m]

  /****************************************************************************************************/
  /*! Initialization of the struct "egy" (of the type ENERGY):*/
  // revision performed on 24.12.2013//

  lg->log("Checking for map files undefined...");

  if (par->output_surfenergy_bin == 1)
    {
      if (geotop::common::Variables::files[fradnet] !=
          geotop::input::gStringNoValue)
        {
          egy->Rn_mean.resize(par->total_pixel + 1, 0.0);
          egy->Rn.resize(par->total_pixel + 1, 0.0);
        }
      else
        {
          count_file_missing++;
          lg->log(
            "File NetRadiationMapFile [usually defined in output_maps/RadNet] NOT "
            "DEFINED",
            geotop::logger::WARNING);
        }
    }
  if (geotop::common::Variables::files[fradLWin] !=
      geotop::input::gStringNoValue)
    {
      egy->LWin_mean.resize(par->total_pixel + 1, 0.0);
      egy->LWin.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log(
        "File InLongwaveRadiationMapFile [usually defined in output_maps/LWin] "
        "NOT DEFINED",
        geotop::logger::WARNING);
    }
  if ((geotop::common::Variables::files[fradLW] !=
       geotop::input::gStringNoValue) ||
      (geotop::common::Variables::files[fradnet] !=
       geotop::input::gStringNoValue))
    {
      egy->LW_mean.resize(par->total_pixel + 1, 0.0);
      egy->LW.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log(
        "File InLongwaveRadiationMapFile[usually defined in output_maps/LWin] "
        "NOT DEFINED",
        geotop::logger::WARNING);
    }

  if ((geotop::common::Variables::files[fradSW] !=
       geotop::input::gStringNoValue) ||
      (geotop::common::Variables::files[fradnet] !=
       geotop::input::gStringNoValue))
    {
      egy->SW_mean.resize(par->total_pixel + 1, 0.0);
      egy->SW.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log("File  with fradSWin identifier NOT DEFINED",
              geotop::logger::WARNING);
    }

  if (geotop::common::Variables::files[fLE] != geotop::input::gStringNoValue)
    {
      egy->ET_mean.resize(par->total_pixel + 1, 0.0);
      egy->LE.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log("File  SurfaceLatentHeatFluxMapFile [= maps/LE] NOT DEFINED",
              geotop::logger::WARNING);
    }
  if (geotop::common::Variables::files[fH] != geotop::input::gStringNoValue)
    {
      egy->H_mean.resize(par->total_pixel + 1, 0.0);
      egy->H.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log("File SurfaceSensibleHeatFluxMapFile [= maps/H] NOT DEFINED",
              geotop::logger::WARNING);
    }
  if (geotop::common::Variables::files[fG] != geotop::input::gStringNoValue)
    {
      egy->SEB_mean.resize(par->total_pixel + 1, 0.0);
      egy->G.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log("File  with fG identifier  NOT DEFINED", geotop::logger::WARNING);
    }
  if (geotop::common::Variables::files[fTs] != geotop::input::gStringNoValue)
    {
      egy->Ts_mean.resize(par->total_pixel + 1, 0.0);
      egy->Ts.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log("File  with fTs identifier  NOT DEFINED", geotop::logger::WARNING);
    }
  if (geotop::common::Variables::files[fradSWin] !=
      geotop::input::gStringNoValue)
    {
      egy->Rswdown_mean.resize(par->total_pixel + 1, 0.0);
      egy->SWin.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log("File  with fradSWin identifier  NOT DEFINED",
              geotop::logger::WARNING);
    }
  if (geotop::common::Variables::files[fradSWinbeam] !=
      geotop::input::gStringNoValue)
    {
      egy->Rswbeam_mean.resize(par->total_pixel + 1, 0.0);
      egy->SWinb.resize(par->total_pixel + 1);
    }
  else
    {
      count_file_missing++;
      lg->log("File  with fradSwinbeam identifier  NOT DEFINED",
              geotop::logger::WARNING);
    }

  if (geotop::common::Variables::files[fshadow] !=
      geotop::input::gStringNoValue)
    {
      egy->nDt_shadow.resize(par->total_pixel + 1, 0.0);
      egy->nDt_sun.resize(par->total_pixel + 1, 0.0);
      egy->shad.resize(par->total_pixel + 1, 0.0);
    }
  else
    {
      count_file_missing++;
      lg->log("File ShadowFractionTimeMapFile [= maps/Shadow???] NOT DEFINED",
              geotop::logger::WARNING);
    }

  if (count_file_missing > 0)
    {
      lg->logsf(geotop::logger::WARNING,
                "%d mapfiles undefined: see above for names of missing files",
                count_file_missing);
    }

  egy->sun = (double *)malloc(12 * sizeof(double));

  if (times->JD_plots.size() > 1)
    {
      if (geotop::common::Variables::files[pH] != geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pHg] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
        {
          egy->Hgplot.resize(par->total_pixel + 1, 0.0);
          egy->Hgp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pH] != geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pHv] !=
          geotop::input::gStringNoValue)
        {
          egy->Hvplot.resize(par->total_pixel + 1, 0.0);
          egy->Hvp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pLE] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pLEg] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
        {
          egy->LEgplot.resize(par->total_pixel + 1, 0.0);
          egy->LEgp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pLE] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pLEv] !=
          geotop::input::gStringNoValue)
        {
          egy->LEvplot.resize(par->total_pixel + 1, 0.0);
          egy->LEvp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pSWin] !=
          geotop::input::gStringNoValue)
        {
          egy->SWinplot.resize(par->total_pixel + 1, 0.0);
          egy->SWinp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pSWg] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
        {
          egy->SWgplot.resize(par->total_pixel + 1, 0.0);
          egy->SWgp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pSWv] !=
          geotop::input::gStringNoValue)
        {
          egy->SWvplot.resize(par->total_pixel + 1, 0.0);
          egy->SWvp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pLWin] !=
          geotop::input::gStringNoValue)
        {
          egy->LWinplot.resize(par->total_pixel + 1, 0.0);
          egy->LWinp.resize(par->total_pixel);
        }
      if (geotop::common::Variables::files[pLWg] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
        {
          egy->LWgplot.resize(par->total_pixel + 1, 0.0);
          egy->LWgp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pLWv] !=
          geotop::input::gStringNoValue)
        {
          egy->LWvplot.resize(par->total_pixel + 1, 0.0);
          egy->LWvp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pTs] !=
          geotop::input::gStringNoValue)
        {
          egy->Tsplot.resize(par->total_pixel + 1, 0.0);
          egy->Tsp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pTg] !=
          geotop::input::gStringNoValue)
        {
          egy->Tgplot.resize(par->total_pixel + 1, 0.0);
          egy->Tgp.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[pTv] !=
          geotop::input::gStringNoValue)
        {
          egy->Tvplot.resize(par->total_pixel + 1, 0.0);
        }
    }

  // vectors used in energy_balance()

  egy->Tgskin_surr.resize(geotop::common::Variables::Nr + 1,
                          geotop::common::Variables::Nc + 1, 0.0);
  egy->SWrefl_surr.resize(geotop::common::Variables::Nr + 1,
                          geotop::common::Variables::Nc + 1, 0.0);
  egy->Dlayer.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                     par->max_glac_layers + 1);
  egy->liq.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                  par->max_glac_layers + 1);
  egy->ice.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                  par->max_glac_layers + 1);
  egy->Temp.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                   par->max_glac_layers + 1);
  egy->deltaw.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                     par->max_glac_layers + 1);
  egy->SWlayer.resize(par->max_snow_layers + 1 + 1);

  //  tolto +1 dalla linea qua sotto 24.12.2013 S.C.&S.E.
  egy->soil_transp_layer.resize(land->root_fraction.getCols());

  egy->dFenergy.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                       par->max_glac_layers + 1);
  egy->udFenergy.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                        par->max_glac_layers + 1);
  egy->Kth0.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                   par->max_glac_layers + 1);
  egy->Kth1.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                   par->max_glac_layers + 1);
  egy->Fenergy.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                      par->max_glac_layers + 1);
  egy->Newton_dir.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                         par->max_glac_layers + 1);
  egy->T0.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                 par->max_glac_layers + 1);
  egy->T1.resize(geotop::common::Variables::Nl + par->max_snow_layers +
                 par->max_glac_layers + 1);
  egy->Tstar.resize(geotop::common::Variables::Nl + 1);
  egy->THETA.resize(geotop::common::Variables::Nl + 1);

  // allocate vector  of soil layer contributions to evaporation (up to
  // z_evap)
  z = 0.;
  l = 0;
  do
    {
      l++;
      z += sl->pa(1, jdz, l);
    }
  while (l < geotop::common::Variables::Nl && z < GTConst::z_evap);

  egy->soil_evap_layer_bare.resize(l + 1, 0.0);
  egy->soil_evap_layer_veg.resize(l + 1, 0.0);

  lg->logf("Soil water evaporates from the first %ld layers",
           egy->soil_evap_layer_bare.size() - 1);
  lg->logf("Soil water transpires from the first %ld layers",
           egy->soil_transp_layer.size() - 1);

  /****************************************************************************************************/
  /*! Completing of the struct "water" (of the type WATER) */

  wat->Voutlandsub = 0.;
  wat->Voutlandsup = 0.;
  wat->Voutbottom = 0.;

  /* Initialization of wat->Pnet (liquid precipitation that reaches the sl
   * surface in mm):*/
  wat->Pnet.resize(geotop::common::Variables::Nr + 1,
                   geotop::common::Variables::Nc + 1, 0.0);

  wat->HN.resize(geotop::common::Variables::Nr + 1,
                 geotop::common::Variables::Nc + 1, 0.0);  // TODO mattiu
  /* Initialization of wat->PrecTot (total precipitation (rain+snow)
   * precipitation):*/
  wat->PrecTot.resize(geotop::common::Variables::Nr + 1,
                      geotop::common::Variables::Nc + 1, 0.0);

  /* Initialization of the matrices with the output of total precipitation and
   * interception:*/
  if (par->output_meteo_bin == 1 && geotop::common::Variables::files[fprec] !=
      geotop::input::gStringNoValue)
    {
      wat->PrTOT_mean.resize(par->total_pixel + 1, 0.0);

      wat->PrSNW_mean.resize(par->total_pixel + 1, 0.0);

      wat->Pt.resize(par->total_pixel + 1);
      wat->Ps.resize(par->total_pixel + 1);
    }

  wat->h_sup.resize(par->total_pixel + 1, 0.0);

  /****************************************************************************************************/
  /*! Initialization of the struct "snow" (of the type SNOW):*/

  /***************************************************************************************************/
  snow->S = new Statevar3D();
  allocate_and_initialize_statevar_3D(
    snow->S, geotop::input::gDoubleNoValue, par->max_snow_layers,
    geotop::common::Variables::Nr, geotop::common::Variables::Nc);

  // initial snow depth
  if (geotop::common::Variables::files[fsn0] != geotop::input::gStringNoValue &&
      geotop::common::Variables::files[fswe0] !=
      geotop::input::gStringNoValue)
    {
      lg->logf("Initial condition on snow depth from file %s",
               geotop::common::Variables::files[fsn0].c_str());

      meteoio_readMap(string(geotop::common::Variables::files[fsn0]), M);

      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              snow->S->Dzl[1][r][c] = M[r][c];
            }
        }

      lg->logf("Initial condition on snow water equivalent from file %s",
               geotop::common::Variables::files[fswe0].c_str());

      meteoio_readMap(string(geotop::common::Variables::files[fswe0]), M);

      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              snow->S->w_ice[1][r][c] = M[r][c];
            }
        }

    }
  else if (geotop::common::Variables::files[fsn0] !=
           geotop::input::gStringNoValue)
    {
      lg->logf("Initial condition on snow depth from file %s",
               geotop::common::Variables::files[fsn0].c_str());

      meteoio_readMap(string(geotop::common::Variables::files[fsn0]), M);

      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              snow->S->Dzl[1][r][c] = M[r][c];
            }
        }

      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                snow->S->w_ice[1][r][c] =
                  snow->S->Dzl[1][r][c] * IT->rhosnow0 / GTConst::rho_w;
            }
        }

    }
  else if (geotop::common::Variables::files[fswe0] !=
           geotop::input::gStringNoValue)
    {
      lg->logf("Initial condition on snow water equivalent from file %s",
               geotop::common::Variables::files[fswe0].c_str());

      meteoio_readMap(string(geotop::common::Variables::files[fswe0]), M);
      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              snow->S->w_ice[1][r][c] = M[r][c];
            }
        }

      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                snow->S->Dzl[1][r][c] =
                  snow->S->w_ice[1][r][c] * GTConst::rho_w / IT->rhosnow0;
            }
        }

    }
  else
    {
      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  snow->S->w_ice[1][r][c] = IT->swe0;
                  snow->S->Dzl[1][r][c] = IT->swe0 * GTConst::rho_w / IT->rhosnow0;
                }
            }
        }
    }

  // Optional reading of snow age in the whole basin
  if (geotop::common::Variables::files[fsnag0] !=
      geotop::input::gStringNoValue)
    {
      lg->logf("Snow age initial condition from file %s",
               geotop::common::Variables::files[fsnag0 + 1].c_str());
      snow->age =
        read_map_vector(geotop::common::Variables::files[fsnag0], top->rc_cont);
    }
  else
    {
      snow->age.resize(par->total_pixel + 1, IT->agesnow0);
    }

  if (times->JD_plots.size() > 1)
    {
      if (geotop::common::Variables::files[pD] != geotop::input::gStringNoValue)
        {
          snow->Dplot.resize(par->total_pixel + 1, 0.0);
        }
    }

  if (par->blowing_snow == 1)
    {
      snow->S_for_BS = new Statevar1D();
      allocate_and_initialize_statevar_1D(
        snow->S_for_BS, geotop::input::gDoubleNoValue, par->max_snow_layers);

      snow->change_dir_wind.resize(Fmaxlong(geotop::common::Variables::Nr + 1,
                                            geotop::common::Variables::Nc + 1));

      snow->Qtrans.resize(geotop::common::Variables::Nr + 1,
                          geotop::common::Variables::Nc + 1, 0.0);

      snow->Qsub.resize(geotop::common::Variables::Nr + 1,
                        geotop::common::Variables::Nc + 1, 0.0);

      snow->Qsalt.resize(geotop::common::Variables::Nr + 1,
                         geotop::common::Variables::Nc + 1);

      snow->Nabla2_Qtrans.resize(geotop::common::Variables::Nr + 1,
                                 geotop::common::Variables::Nc + 1);

      snow->Qsub_x.resize(geotop::common::Variables::Nr + 1,
                          geotop::common::Variables::Nc + 1);

      snow->Qsub_y.resize(geotop::common::Variables::Nr + 1,
                          geotop::common::Variables::Nc + 1);

      snow->Qtrans_x.resize(geotop::common::Variables::Nr + 1,
                            geotop::common::Variables::Nc + 1);

      snow->Qtrans_y.resize(geotop::common::Variables::Nr + 1,
                            geotop::common::Variables::Nc + 1);

      if (par->output_snow_bin == 1)
        {
          snow->Wtrans_plot.resize(geotop::common::Variables::Nr + 1,
                                   geotop::common::Variables::Nc + 1);
          snow->Wsubl_plot.resize(geotop::common::Variables::Nr + 1,
                                  geotop::common::Variables::Nc + 1);

          for (r = 1; r <= geotop::common::Variables::Nr; r++)
            {
              for (c = 1; c <= geotop::common::Variables::Nc; c++)
                {
                  if ((long)land->LC[r][c] == geotop::input::gDoubleNoValue)
                    {
                      snow->Wtrans_plot[r][c] = geotop::input::gDoubleNoValue;
                      snow->Wsubl_plot[r][c] = geotop::input::gDoubleNoValue;

                    }
                  else
                    {
                      snow->Wtrans_plot[r][c] = 0.0;
                      snow->Wsubl_plot[r][c] = 0.0;
                    }
                }
            }
        }
    }

  if (par->output_snow_bin == 1)
    {
      // if(geotop::common::Variables::files[fsnowmelt] !=
      // geotop::input::gStringNoValue){
      snow->MELTED.resize(par->total_pixel + 1, 0.0);
      snow->melted.resize(par->total_pixel + 1);
      // }
      // if(geotop::common::Variables::files[fsnowsubl] !=
      // geotop::input::gStringNoValue){
      snow->SUBL.resize(par->total_pixel + 1, 0.0);
      snow->subl.resize(par->total_pixel + 1);
      // }
      // if(geotop::common::Variables::files[fsndur] !=
      // geotop::input::gStringNoValue){
      snow->t_snow.resize(par->total_pixel + 1, 0.0);
      snow->yes.resize(par->total_pixel + 1, 0.0);
      // }

      // if(geotop::common::Variables::files[fHN] !=
      // geotop::input::gStringNoValue){//TODO mattiu
      snow->HNcum.resize(par->total_pixel + 1, 0.0);
      // }//end mattiu
    }

  for (r = 1; r <= geotop::common::Variables::Nr; r++)
    {
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
            {
              // Adjusting snow init depth in case of steep slope (contribution by
              // Stephan Gruber)
              if (par->snow_curv > 0 && top->slope[r][c] > par->snow_smin)
                {
                  if (top->slope[r][c] <= par->snow_smax)
                    {
                      k_snowred = (exp(-pow(top->slope[r][c] - par->snow_smin, 2.) /
                                       par->snow_curv) -
                                   exp(-pow(par->snow_smax, 2.) / par->snow_curv));
                    }
                  else
                    {
                      k_snowred = 0.0;
                    }
                  snow->S->Dzl[1][r][c] *= k_snowred;
                  snow->S->w_ice[1][r][c] *= k_snowred;
                }

              D = snow->S->Dzl[1][r][c];
              SWE = snow->S->w_ice[1][r][c];

              if (D < 0 || SWE < 0)
                {
                  f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                  fprintf(f,
                          "Error: negative initial snow depth %e or snow water "
                          "equivalent %e\n",
                          D, SWE);
                  fclose(f);
                  lg->logsf(geotop::logger::ERROR,
                            "Error: negative initial snow depth %12g or snow water "
                            "equivalent %12g",
                            D, SWE);
                  lg->log("Geotop failed. See failing report (7).",
                          geotop::logger::CRITICAL);
                  exit(1);

                }
              else if (D < 1.E-5 && SWE > 1.E-5)
                {
                  f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                  fprintf(f,
                          "Error: Initial snow water equivalent %e > 0 and initial "
                          "snow depth %e\n",
                          SWE, D);
                  fclose(f);
                  lg->logsf(geotop::logger::ERROR,
                            "Initial snow water equivalent %12g > 0 and initial snow "
                            "depth %12g",
                            SWE, D);
                  lg->log("Geotop failed. See failing report (8).",
                          geotop::logger::CRITICAL);
                  exit(1);

                }
              else if (D > 1.E-5 && SWE < 1.E-5)
                {
                  f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                  fprintf(f,
                          "Error: Initial snow depth %e > 0 and initial snow water "
                          "equivalent %e\n",
                          D, SWE);
                  fclose(f);
                  lg->logsf(geotop::logger::ERROR,
                            "Initial snow depth %12g > 0 and initial snow water "
                            "equivalent %12g\n",
                            D, SWE);
                  lg->log("Geotop failed. See failing report (9).",
                          geotop::logger::CRITICAL);
                  exit(1);

                }
              else if (D > 1.E-5 || SWE > 1.E-5)
                {
                  snow->age[top->j_cont[r][c]] *= 86400.0;  // now in [s]
                  if (SWE <= par->max_weq_snow * par->max_snow_layers)
                    {
                      i = floor(SWE / par->max_weq_snow);

                      if (i > 0)
                        {
                          for (n = 1; n <= i; n++)
                            {
                              snow->S->w_ice[n][r][c] = par->max_weq_snow;
                            }

                          if (SWE - i * par->max_weq_snow > 0.1 * par->max_weq_snow)
                            {
                              snow->S->w_ice[i + 1][r][c] = SWE - i * par->max_weq_snow;
                              snow->S->lnum[r][c] = i + 1;
                            }
                          else
                            {
                              snow->S->w_ice[i][r][c] += (SWE - i * par->max_weq_snow);
                              snow->S->lnum[r][c] = i;
                            }

                        }
                      else
                        {
                          snow->S->w_ice[1][r][c] = SWE;
                          snow->S->lnum[r][c] = 1;
                        }

                    }
                  else
                    {
                      snow->S->lnum[r][c] = par->max_snow_layers;

                      for (n = 1; n <= par->max_snow_layers; n++)
                        {
                          a = 0;

                          for (size_t ii = 1; ii < par->inf_snow_layers.size(); ii++)
                            {
                              if (n == abs(par->inf_snow_layers[ii])) a = 1;
                            }

                          if (a == 0)
                            {
                              snow->S->w_ice[n][r][c] = par->max_weq_snow;
                            }
                          else
                            {
                              snow->S->w_ice[n][r][c] =
                                (SWE -
                                 par->max_weq_snow * (par->max_snow_layers -
                                                      (par->inf_snow_layers.size() - 1))) /
                                (par->inf_snow_layers.size() - 1);
                            }
                        }
                    }

                  for (n = 1; n <= snow->S->lnum[r][c]; n++)
                    {
                      snow->S->Dzl[n][r][c] = D * (snow->S->w_ice[n][r][c] / SWE);
                      snow->S->T[n][r][c] = IT->Tsnow0;
                    }
                }

              non_dimensionalize_snowage(&(snow->age[top->j_cont[r][c]]), IT->Tsnow0);

              if (par->point_sim == 1)
                {
                  maxSWE = par->maxSWE[r][c];
                }
              else
                {
                  maxSWE = 1.E10;
                }

              //                f = fopen(geotop::common::Variables::logfile.c_str(),
              //                "a");
              snow_layer_combination(par->alpha_snow, r, c, snow->S, -0.1,
                                     par->inf_snow_layers, par->max_weq_snow, maxSWE);
              //                fclose(f);
            }
        }
    }

  if (par->delay_day_recover > 0)
    {
      snow->S->type.resize(snow->S->type.getRows(), snow->S->type.getCols(), 2);

      assign_recovered_map_long(
        par->recover, geotop::common::Variables::files[rns], snow->S->lnum);
      assign_recovered_map_vector(par->recover,
                                  geotop::common::Variables::files[rsnag],
                                  snow->age, top->rc_cont);
      assign_recovered_tensor(
        1, par->recover, geotop::common::Variables::files[rDzs], snow->S->Dzl);
      assign_recovered_tensor(
        1, par->recover, geotop::common::Variables::files[rwls], snow->S->w_liq);
      assign_recovered_tensor(
        1, par->recover, geotop::common::Variables::files[rwis], snow->S->w_ice);
      assign_recovered_tensor(1, par->recover,
                              geotop::common::Variables::files[rTs], snow->S->T);
    }

  f = fopen(geotop::common::Variables::logfile.c_str(), "a");
  for (r = 1; r <= geotop::common::Variables::Nr; r++)
    {
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
            {
              snow_layer_combination(par->alpha_snow, r, c, snow->S, 0.,
                                     par->inf_snow_layers, par->max_weq_snow, maxSWE);
            }
        }
    }
  //    fclose(f);

  /****************************************************************************************************/
  /*! Initialization of the struct "glac" (of the type GLACIER):*/

  /***************************************************************************************************/
  /*! Optional reading of glacier depth in the whole basin ("GLACIER0"):    */
  if (geotop::common::Variables::files[fgl0] != geotop::input::gStringNoValue &&
      par->max_glac_layers == 0)
    {
      lg->log("Glacier map present, but glacier represented with 0 layers",
              geotop::logger::WARNING);
    }

  if (par->max_glac_layers == 0 && IT->Dglac0 > 0)
    {
      lg->log(
        "You have chosen 0 glacier layers in block 10 in the parameter file, but "
        "you assigned a value of the glacier depth. The latter will be ignored.",
        geotop::logger::WARNING);
    }

  // If the max number of glacier layers is greater than 1, the matrices (or
  // tensors) lnum, Dzl. w_liq, w_ice, T and print matrices are defined,
  // according to the respective flags
  if (par->max_glac_layers > 0)
    {
      glac->G = new Statevar3D();
      allocate_and_initialize_statevar_3D(
        glac->G, geotop::input::gDoubleNoValue, par->max_glac_layers,
        geotop::common::Variables::Nr, geotop::common::Variables::Nc);

      if (geotop::common::Variables::files[fgl0] !=
          geotop::input::gStringNoValue)
        {
          lg->logf("Glacier initial condition from file %s",
                   geotop::common::Variables::files[fgl0].c_str());  // FIXED fgl0+1

          meteoio_readMap(string(geotop::common::Variables::files[fgl0]), M);

          for (r = 1; r <= geotop::common::Variables::Nr; r++)
            {
              for (c = 1; c <= geotop::common::Variables::Nc; c++)
                {
                  glac->G->Dzl[1][r][c] = M[r][c];
                }
            }

          for (r = 1; r <= geotop::common::Variables::Nr; r++)
            {
              for (c = 1; c <= geotop::common::Variables::Nc; c++)
                {
                  if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                    glac->G->w_ice[1][r][c] =
                      glac->G->Dzl[1][r][c] * IT->rhoglac0 / GTConst::rho_w;
                }
            }

        }
      else
        {
          for (r = 1; r <= geotop::common::Variables::Nr; r++)
            {
              for (c = 1; c <= geotop::common::Variables::Nc; c++)
                {
                  if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                    {
                      glac->G->Dzl[1][r][c] = IT->Dglac0;
                      glac->G->w_ice[1][r][c] =
                        IT->Dglac0 * IT->rhoglac0 / GTConst::rho_w;
                    }
                }
            }
        }

      if (par->output_glac_bin == 1)
        {
          if (geotop::common::Variables::files[fglacmelt] !=
              geotop::input::gStringNoValue)
            {
              glac->MELTED.resize(par->total_pixel + 1, 0.0);
              glac->melted.resize(par->total_pixel + 1);
            }
          if (geotop::common::Variables::files[fglacsubl] !=
              geotop::input::gStringNoValue)
            {
              glac->SUBL.resize(par->total_pixel + 1, 0.0);
              glac->subl.resize(par->total_pixel + 1);
            }
        }

      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  // SKIP Adjusting snow init depth in case of steep slope

                  D = glac->G->Dzl[1][r][c];
                  SWE = glac->G->w_ice[1][r][c];  // phisically GWE but named SWE

                  if (D < 0 || SWE < 0)
                    {
                      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                      fprintf(f,
                              "Error: negative initial glacier depth %e or glacier water "
                              "equivalent %e\n",
                              D, SWE);
                      fclose(f);
                      lg->logsf(geotop::logger::ERROR,
                                "Error: negative initial glacier depth %12g or snow "
                                "glacier equivalent %12g",
                                D, SWE);
                      lg->log("Geotop failed. See failing report (10).",
                              geotop::logger::CRITICAL);
                      exit(1);

                    }
                  else if (D < 1.E-5 && SWE > 1.E-5)
                    {
                      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                      fprintf(f,
                              "Error: Initial glacier water equivalent %e > 0 and "
                              "initial glacier depth %e\n",
                              SWE, D);
                      fclose(f);
                      lg->logsf(geotop::logger::ERROR,
                                "Initial glacier water equivalent %12g > 0 and initial "
                                "glacier depth %12g",
                                SWE, D);
                      lg->log("Geotop failed. See failing report (11).",
                              geotop::logger::CRITICAL);
                      exit(1);

                    }
                  else if (D > 1.E-5 && SWE < 1.E-5)
                    {
                      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                      fprintf(f,
                              "Error: Initial glacier depth %e > 0 and initial glacier "
                              "water equivalent %e\n",
                              D, SWE);
                      fclose(f);
                      lg->logsf(geotop::logger::ERROR,
                                "Initial glacier depth %12g > 0 and initial glacier "
                                "water equivalent %12g\n",
                                D, SWE);
                      lg->log("Geotop failed. See failing report (12).",
                              geotop::logger::CRITICAL);
                      exit(1);

                    }
                  else if (D > 1.E-5 || SWE > 1.E-5)
                    {
                      if (SWE < par->max_weq_glac * par->max_glac_layers)
                        {
                          i = floor(SWE / par->max_weq_glac);

                          if (i > 0)
                            {
                              for (n = 1; n <= i; n++)
                                {
                                  glac->G->w_ice[n][r][c] = par->max_weq_glac;
                                }

                              if (SWE - i * par->max_weq_glac > 0.1 * par->max_weq_glac)
                                {
                                  glac->G->w_ice[i + 1][r][c] = SWE - i * par->max_weq_glac;
                                  glac->G->lnum[r][c] = i + 1;
                                }
                              else
                                {
                                  glac->G->w_ice[i][r][c] += (SWE - i * par->max_weq_glac);
                                  glac->G->lnum[r][c] = i;
                                }

                            }
                          else
                            {
                              glac->G->w_ice[1][r][c] = SWE;
                              glac->G->lnum[r][c] = 1;
                            }

                        }
                      else
                        {
                          glac->G->lnum[r][c] = par->max_glac_layers;

                          for (n = 1; n <= par->max_glac_layers; n++)
                            {
                              a = 0;

                              for (size_t ii = 1; ii < par->inf_glac_layers.size(); ii++)
                                {
                                  if (n == abs(par->inf_glac_layers[ii])) a = 1;
                                }

                              if (a == 0)
                                {
                                  glac->G->w_ice[n][r][c] = par->max_weq_glac;
                                }
                              else
                                {
                                  glac->G->w_ice[n][r][c] =
                                    (SWE -
                                     par->max_weq_glac * (par->max_glac_layers -
                                                          (par->inf_glac_layers.size() - 1))) /
                                    (par->inf_glac_layers.size() - 1);
                                }
                            }
                        }

                      for (n = 1; n <= glac->G->lnum[r][c]; n++)
                        {
                          glac->G->Dzl[n][r][c] = D * (glac->G->w_ice[n][r][c] / SWE);
                          glac->G->T[n][r][c] = IT->Tglac0;
                        }
                    }

                  if (par->point_sim == 1)
                    {
                      maxSWE = 1.E10;  // for snow: par->maxSWE[r][c]; // TODO
                    }
                  else
                    {
                      maxSWE = 1.E10;
                    }

                  f = fopen(geotop::common::Variables::logfile.c_str(), "a");
                  snow_layer_combination(par->alpha_snow, r, c, glac->G, -0.1,
                                         par->inf_glac_layers, par->max_weq_glac,
                                         maxSWE);
                  //                    fclose(f);
                }
            }
        }

      if (par->delay_day_recover > 0)
        {
          glac->G->type.resize(glac->G->type.getRows(), glac->G->type.getCols(), 2);

          assign_recovered_map_long(
            par->recover, geotop::common::Variables::files[rni], glac->G->lnum);
          assign_recovered_tensor(
            1, par->recover, geotop::common::Variables::files[rDzi], glac->G->Dzl);
          assign_recovered_tensor(1, par->recover,
                                  geotop::common::Variables::files[rwli],
                                  glac->G->w_liq);
          assign_recovered_tensor(1, par->recover,
                                  geotop::common::Variables::files[rwii],
                                  glac->G->w_ice);
          assign_recovered_tensor(
            1, par->recover, geotop::common::Variables::files[rTi], glac->G->T);
        }
    }

  //***************************************************************************************************
  // Filling up of the struct "met" (of the type METEO):

  met->Tgrid.resize(geotop::common::Variables::Nr + 1,
                    geotop::common::Variables::Nc + 1, 5.);

  met->Pgrid.resize(geotop::common::Variables::Nr + 1,
                    geotop::common::Variables::Nc + 1, GTConst::Pa0);

  met->RHgrid.resize(geotop::common::Variables::Nr + 1,
                     geotop::common::Variables::Nc + 1, 0.7);

  met->Vgrid.resize(geotop::common::Variables::Nr + 1,
                    geotop::common::Variables::Nc + 1, par->Vmin);

  met->Vdir.resize(geotop::common::Variables::Nr + 1,
                   geotop::common::Variables::Nc + 1, 0.0);

  met->ILWRgrid.resize(geotop::common::Variables::Nr + 1,
                       geotop::common::Variables::Nc + 1, 0.0);

  if (par->output_meteo_bin == 1)
    {
      if (geotop::common::Variables::files[fTa] !=
          geotop::input::gStringNoValue)
        {
          met->Tamean.resize(par->total_pixel + 1);
        }
      if (geotop::common::Variables::files[fwspd] !=
          geotop::input::gStringNoValue)
        {
          met->Vspdmean.resize(par->total_pixel + 1, 0.0);
        }
      if (geotop::common::Variables::files[fwdir] !=
          geotop::input::gStringNoValue)
        {
          met->Vdirmean.resize(par->total_pixel + 1, 0.0);
        }
      if (geotop::common::Variables::files[frh] !=
          geotop::input::gStringNoValue)
        {
          met->RHmean.resize(par->total_pixel + 1, 0.0);
        }
    }

  // plot output
  if (times->JD_plots.size() > 1)
    {
      if (geotop::common::Variables::files[pTa] !=
          geotop::input::gStringNoValue)
        {
          met->Taplot.resize(par->total_pixel + 1, 0.0);
        }
      if (geotop::common::Variables::files[pRH] !=
          geotop::input::gStringNoValue)
        {
          met->RHplot.resize(par->total_pixel + 1, 0.0);
        }
      if (geotop::common::Variables::files[pVspd] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[pVdir] !=
          geotop::input::gStringNoValue)
        {
          met->Vxplot.resize(par->total_pixel + 1, 0.0);
          met->Vyplot.resize(par->total_pixel + 1, 0.0);
        }
    }

  /****************************************************************************************************/

  delete (IT);

  if (par->point_sim != 1)
    {
      cont_nonzero_values_matrix2(&i, &j, cnet, land->LC, top->lrc_cont,
                                  top->i_cont, par->total_pixel);

      top->Li.resize(i + 1);
      top->Lp.resize(j + 1);

      wat->Lx.resize(i + 1);
      cont_nonzero_values_matrix3(top->Lp, top->Li, cnet, land->LC, top->lrc_cont,
                                  top->i_cont, par->total_pixel);

    }
  else
    {
      i = geotop::common::Variables::Nl;
      j = geotop::common::Variables::Nl + 1;
      top->Li.resize(i + 1);
      top->Lp.resize(j + 1);

      wat->Lx.resize(i + 1);

      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          top->Li[l] = l + 1;
          top->Lp[l] = l;
        }
      top->Lp[l] = i;
    }

  wat->H0.resize(j + 1);
  wat->H1.resize(j + 1);
  wat->dH.resize(j + 1);
  wat->B.resize(j + 1);
  wat->f.resize(j + 1);
  wat->df.resize(j + 1);
  wat->Kbottom.resize(geotop::common::Variables::Nr + 1,
                      geotop::common::Variables::Nc + 1, 0.);

  wat->Klat.resize(top->BC_DepthFreeSurface.size(),
                   geotop::common::Variables::Nl + 1, 0.0);

  //    fclose(flog);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void read_inputmaps(Topo *top,
                    Land *land,
                    Soil *sl,
                    Par *par,
                    InitTools * /*IT*/,
                    mio::IOManager &iomanager)
{
  size_t r, c, i, cont;
  GeoMatrix<double> M;
  GeoMatrix<short> curv;
  short flag;
  string temp;
  double min, max;
  FILE *f;
  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  // we check here  if DEM file read by meteoio is the same specified in
  // geotop.inpts introduced here by S.C. 17.01.2017 eventually we should merge
  // geotop.intps ands io_it.ini and make this check useless
  //
  mio::Config cfg = iomanager.getConfig();
  std::string filename_dem_path = cfg.get("DEMFILE", "Input");
  std::string filename_dem = mio::IOUtils::getFilename(filename_dem_path);
  std::string filename_geotop =
    mio::IOUtils::getFilename(geotop::common::Variables::files[fdem] + ".asc");
  lg->logsf(geotop::logger::NOTICE,
            "DEM FILE used by METEO-IO (defined in io_it.ini)   : %s",
            filename_dem.c_str());
  lg->logsf(geotop::logger::NOTICE,
            "DEM FILE used by GEOTOP   (defined in geotop.inpts): %s",
            filename_geotop.c_str());
  if (filename_dem != filename_geotop)
    {
      lg->log(
        " You specified two different file for DEM in two different files.. "
        "please define just one",
        geotop::logger::ERROR);
      lg->log("Geotop failed.", geotop::logger::CRITICAL);
      exit(1);
    }

  meteoio_readDEM(top->Z0);

  //  reading TOPOGRAPHY
  flag = file_exists(fdem);

  if (flag == 1)
    {
      // filtering
      M.resize(top->Z0.getRows(), top->Z0.getCols());
      multipass_topofilter(par->lowpass, top->Z0, M,
                           geotop::input::gDoubleNoValue, 1);
      top->Z0 = M;

      //  calculate East and North
      top->East.resize(top->Z0.getRows(), top->Z0.getCols());
      top->North.resize(top->Z0.getRows(), top->Z0.getCols());

      // to check <= or just < ::TODO
      for (r = 1; r < top->Z0.getRows(); r++)
        {
          for (c = 1; c < top->Z0.getCols(); c++)
            {
              top->East[r][c] = geotop::common::Variables::UV->U[4] +
                                (c - 0.5) * geotop::common::Variables::UV->U[2];
              top->North[r][c] = geotop::common::Variables::UV->U[3] +
                                 ((top->Z0.getRows() - 1) - (r - 0.5)) *
                                 geotop::common::Variables::UV->U[1];
            }
        }

    }
  else
    {
      lg->log(
        "It is impossible to proceed without giving the digital elevation model",
        geotop::logger::ERROR);
      lg->log("Geotop failed. See failing report (11).",
              geotop::logger::CRITICAL);
      exit(1);
    }

  //  reading LAND COVER TYPE
  flag = file_exists(flu);
  if (flag == 1)
    {
      meteoio_readMap(string(geotop::common::Variables::files[flu]), land->LC);

      //  Check borders
      for (r = 1; r < land->LC.getRows(); r++)
        {
          land->LC[r][1] = geotop::input::gDoubleNoValue;
          land->LC[r][land->LC.getCols() - 1] = geotop::input::gDoubleNoValue;
        }
      for (c = 1; c < land->LC.getCols(); c++)
        {
          land->LC[1][c] = geotop::input::gDoubleNoValue;
          land->LC[land->LC.getRows() - 1][c] = geotop::input::gDoubleNoValue;
        }
      for (r = 1; r < land->LC.getRows(); r++)
        {
          for (c = 1; c < land->LC.getCols(); c++)
            {
              if (land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  if (land->LC[r][c] < 1. ||
                      land->LC[r][c] > (double)(par->n_landuses))
                    {
                      lg->log(
                        "It is not possible to assign Value < 1 or > n_landuses to the "
                        "land cover type",
                        geotop::logger::ERROR);
                      lg->log("Geotop failed. See failing report (12).",
                              geotop::logger::CRITICAL);
                      exit(1);
                    }
                }
            }
        }

      // Land use is the official mask
      for (r = 1; r < land->LC.getRows(); r++)
        {
          for (c = 1; c < land->LC.getCols(); c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  if ((long)top->Z0[r][c] == geotop::input::gDoubleNoValue)
                    {
                      lg->log("Land use mask include DTM novalue pixels",
                              geotop::logger::WARNING);
                    }
                }
            }
        }

    }
  else
    {
      //  Write land->LC (land cover)
      lg->log("Land cover type assumed to be always 1");
      copydoublematrix_const(1.0, top->Z0, land->LC,
                             geotop::input::gDoubleNoValue);

      for (r = 1; r < land->LC.getRows(); r++)
        {
          land->LC[r][1] = geotop::input::gDoubleNoValue;
          land->LC[r][land->LC.getCols() - 1] = geotop::input::gDoubleNoValue;
        }
      for (c = 1; c < land->LC.getCols(); c++)
        {
          land->LC[1][c] = geotop::input::gDoubleNoValue;
          land->LC[land->LC.getRows() - 1][c] = geotop::input::gDoubleNoValue;
        }
    }
  // missing write of file : To CHECK

  lg->logf("par->state_pixel=%ld\n", par->state_pixel);
  // plot results for a certain numbevr of pixel

  if (par->state_pixel == 1)
    {
      par->rc.resize(par->chkpt.getRows(), 2 + 1);
      par->IDpoint.resize(par->chkpt.getRows());

      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          par->rc[i][1] =
            row(par->chkpt[i][ptY], top->Z0.getRows() - 1,
                geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
          par->rc[i][2] =
            col(par->chkpt[i][ptX], top->Z0.getCols() - 1,
                geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

          if (par->rc[i][1] == geotop::input::gDoubleNoValue ||
              par->rc[i][2] == geotop::input::gDoubleNoValue)
            {
              lg->logsf(geotop::logger::ERROR, "Point #%4ld is out of the domain", i);
              lg->log("Geotop failed. See failing report.", geotop::logger::CRITICAL);
              exit(1);
            }

          if ((long)land->LC[par->rc[i][1]][par->rc[i][2]] ==
              geotop::input::gDoubleNoValue)
            {
              lg->logsf(geotop::logger::ERROR,
                        "Point #%4ld corresponds to NOVALUE pixel", i);
              lg->log("Geotop failed. See failing report.", geotop::logger::CRITICAL);
              exit(1);
            }

          lg->writefAll("i:%ld %f\n", i, par->chkpt[i][ptID]);

          if ((long)par->chkpt[i][ptID] != geotop::input::gDoubleNoValue)
            {
              par->IDpoint[i] = (long)par->chkpt[i][ptID];
              lg->writefAll("A i:%ld %ld\n", i, par->IDpoint[i]);
            }
          else
            {
              par->IDpoint[i] = i;
              lg->writefAll("B i:%ld %ld\n", i, par->IDpoint[i]);
            }
        }
    }

  /****************************************************************************************************/

  // reading SKY VIEW FACTOR
  flag = file_exists(fsky);
  if (flag == 1)
    {
      meteoio_readMap(string(geotop::common::Variables::files[fsky]), top->sky);
    }
  else     /*The sky view factor file "top->sky" must be calculated*/
    {
      top->sky.resize(top->Z0.getRows(), top->Z0.getCols());
      if (par->sky == 0)
        {
          top->sky.resize(top->Z0.getRows(), top->Z0.getCols(), 1.);
        }
      else
        {
          curv.resize(top->Z0.getRows(), top->Z0.getCols());
          nablaquadro_mask(top->Z0, curv, geotop::common::Variables::UV->U,
                           geotop::common::Variables::UV->V);
          sky_view_factor(top->sky, 36, geotop::common::Variables::UV, top->Z0,
                          curv, geotop::input::gDoubleNoValue);
        }
    }
  // missing writin of sky with respect to c se27xx version

  /****************************************************************************************************/

  // reading DELAY
  flag = file_exists(fdelay);
  if (flag == 1)
    {
      meteoio_readMap(string(geotop::common::Variables::files[fdelay]),
                      land->delay);
    }
  else
    {
      land->delay.resize(top->Z0.getRows(), top->Z0.getCols(), 0);
    }

  /****************************************************************************************************/

  // reading SOIL MAP
  flag = file_exists(fsoil);
  if (flag == 1)
    {
      meteoio_readMap(string(geotop::common::Variables::files[fsoil]), M);

      copylong_doublematrix(sl->type, M);
      for (r = 1; r < land->LC.getRows(); r++)
        {
          for (c = 1; c < land->LC.getCols(); c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  if (sl->type[r][c] < 1 || sl->type[r][c] > par->nsoiltypes)
                    {
                      lg->log(
                        "It is not possible to assign Value < 1 or > nsoiltypes to the "
                        "soil type map",
                        geotop::logger::ERROR);
                      lg->log("GEOtop failed. See failing report (13).",
                              geotop::logger::CRITICAL);
                      exit(1);
                    }
                }
            }
        }
    }
  else
    {
      copydoublematrix_const(par->soil_type_land_default, land->LC, M,
                             geotop::input::gDoubleNoValue);

      copylong_doublematrix(sl->type, M);
    }

  /****************************************************************************************************/
  // SLOPE
  top->dzdE.resize(land->LC.getRows(), land->LC.getCols());
  top->dzdN.resize(land->LC.getRows(), land->LC.getCols());

  find_slope(geotop::common::Variables::UV->U[1],
             geotop::common::Variables::UV->U[2], top->Z0, top->dzdE, top->dzdN,
             geotop::input::gDoubleNoValue);

  flag = file_exists(fslp);
  if (flag == 1)
    {
      meteoio_readMap(string(geotop::common::Variables::files[fslp]), top->slope);
    }
  else
    {
      find_max_slope(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue,
                     top->slope);
    }

  find_min_max(top->slope, geotop::input::gDoubleNoValue, &max, &min);

  lg->logf("Slope Min:%12g (%12g deg) Max:%12g (%12g deg)",
           tan(min * GTConst::Pi / 180.), min, tan(max * GTConst::Pi / 180.),
           max);

  /****************************************************************************************************/
  // ASPECT

  flag = file_exists(fasp);
  if (flag == 1)
    {
      meteoio_readMap(string(geotop::common::Variables::files[fasp]),
                      top->aspect);
    }
  else
    {
      find_aspect(top->Z0, top->dzdE, top->dzdN, geotop::input::gDoubleNoValue,
                  top->aspect);
    }
  // why this below ??? Why are we rewriting the ASPECT Map ???? WHY ????
  // S.C. 10.12.2016 checked of  0 after discussion with S.Endrizzi

  if (flag == 0)
    write_map(geotop::common::Variables::files[fasp], 0, par->format_out,
              top->aspect, geotop::common::Variables::UV,
              geotop::input::gDoubleNoValue);

  /****************************************************************************************************/
  // curvature

  top->curvature1.resize(top->Z0.getRows(), top->Z0.getCols());
  top->curvature2.resize(top->Z0.getRows(), top->Z0.getCols());
  top->curvature3.resize(top->Z0.getRows(), top->Z0.getCols());
  top->curvature4.resize(top->Z0.getRows(), top->Z0.getCols());

  // filtering
  M.resize(top->Z0.getRows(), top->Z0.getCols());
  multipass_topofilter(par->lowpass_curvatures, top->Z0, M,
                       geotop::input::gDoubleNoValue, 1);
  curvature(geotop::common::Variables::UV->U[1],
            geotop::common::Variables::UV->U[2], M, top->curvature1,
            top->curvature2, top->curvature3, top->curvature4,
            geotop::input::gDoubleNoValue);

  if (geotop::common::Variables::files[fcurv] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fcurv] + string("N-S");
      write_map(temp, 0, par->format_out, top->curvature1,
                geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

      temp = geotop::common::Variables::files[fcurv] + string("W-E");
      write_map(temp, 0, par->format_out, top->curvature2,
                geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

      temp = geotop::common::Variables::files[fcurv] + string("NW-SE");
      write_map(temp, 0, par->format_out, top->curvature3,
                geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

      temp = geotop::common::Variables::files[fcurv] + string("NE-SW");
      write_map(temp, 0, par->format_out, top->curvature4,
                geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
    }

  find_min_max(top->curvature1, geotop::input::gDoubleNoValue, &max, &min);
  lg->logf("Curvature N-S Min:%12g  Max:%12g", min, max);

  find_min_max(top->curvature2, geotop::input::gDoubleNoValue, &max, &min);
  lg->logf("Curvature W-E Min:%12g  Max:%12g", min, max);

  find_min_max(top->curvature3, geotop::input::gDoubleNoValue, &max, &min);
  lg->logf("Curvature NW-SE Min:%12g  Max:%12g", min, max);

  find_min_max(top->curvature4, geotop::input::gDoubleNoValue, &max, &min);

  lg->logf("Curvature NE-SW Min:%12g  Max:%12g", min, max);

  /****************************************************************************************************/
  // Channel network (in top->pixel_type)

  // pixel type = 0 land pixel (if it is on the border, the border is
  // impermeable, water is free only on the surface)  pixel type = 1 or 2 land
  // pixel (it it is on the border, the border is permeable above an user-defined
  // elevation in the saturated part, weir-wise)  pixel type = 10 channel pixel
  // (if it is on the border, the border is impermeable, water is free only on
  // the surface)  pixel type = -1 land pixel where an incoming discharge from
  // outside is considered (as rain) (TO REMOVE: no longer used  pixel type = 11
  // or 12 the sum of the previous two

  flag = file_exists(fnet);
  if (flag == 1)
    {
      meteoio_readMap(string(geotop::common::Variables::files[fnet]), M);

      copyshort_doublematrix(top->pixel_type, M);

      cont = 0;
      for (r = 1; r < top->Z0.getRows(); r++)
        {
          for (c = 1; c < top->Z0.getCols(); c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  if (top->pixel_type[r][c] != 0 && top->pixel_type[r][c] != 1 &&
                      top->pixel_type[r][c] != 2 && top->pixel_type[r][c] != 10 &&
                      top->pixel_type[r][c] != 11 && top->pixel_type[r][c] != 12)
                    {
                      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
                      fprintf(f,
                              "Error: Only the following values are admitted in the "
                              "network map: 0, 1,2,10,11,12\n");
                      fclose(f);
                      lg->log(
                        "Only the following values are admitted in the network map: "
                        "0,1,10,11,12",
                        geotop::logger::ERROR);
                      lg->log("Geotop failed. See failing report (14).",
                              geotop::logger::CRITICAL);
                      exit(1);
                    }
                  if (top->pixel_type[r][c] >= 10) cont++;
                }
            }
        }

      lg->logf("Channel networks has %ld pixels set to channel", cont);
      // this below WHY commented: why do we have to write back the files just
      // read ?
      // if(flag >= 0) write_map(geotop::common::Variables::files[fnet], 1,
      // par->format_out, M, geotop::common::Variables::UV,
      // geotop::input::gDoubleNoValue);

    }
  else
    {
      top->pixel_type.resize(land->LC.getRows(), land->LC.getCols(), 0);
    }

  /****************************************************************************************************/

  // border
  top->is_on_border.resize(land->LC.getRows(), land->LC.getCols());
  for (r = 1; r < land->LC.getRows(); r++)
    {
      for (c = 1; c < land->LC.getCols(); c++)
        {
          if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
            {
              top->is_on_border[r][c] =
                is_boundary(r, c, land->LC, geotop::input::gDoubleNoValue);
            }
          else
            {
              top->is_on_border[r][c] = -1;
            }
        }
    }

  // count the pixels having pixel_type = 1, 2 or 11 ans 12
  cont = 0;

  for (r = 1; r < top->Z0.getRows(); r++)
    {
      for (c = 1; c < top->Z0.getCols(); c++)
        {
          if (top->is_on_border[r][c] == 1)
            {
              if (top->pixel_type[r][c] == 1 || top->pixel_type[r][c] == 2 ||
                  top->pixel_type[r][c] == 11 || top->pixel_type[r][c] == 12)
                cont++;
            }
        }
    }

  top->BC_counter.resize(top->Z0.getRows(), top->Z0.getCols(), 0);

  if (cont > 0)
    {
      top->BC_DepthFreeSurface.resize(cont + 1);
      cont = 0;
      for (r = 1; r < top->Z0.getRows(); r++)
        {
          for (c = 1; c < top->Z0.getCols(); c++)
            {
              if (top->is_on_border[r][c] == 1)
                {
                  if (top->pixel_type[r][c] == 1 || top->pixel_type[r][c] == 2 ||
                      top->pixel_type[r][c] == 11 || top->pixel_type[r][c] == 12)
                    {
                      cont++;
                      top->BC_counter[r][c] = cont;
                      top->BC_DepthFreeSurface[cont] = par->DepthFreeSurface;  //[mm]
                    }
                }
            }
        }
    }
  else
    {
      top->BC_DepthFreeSurface.resize(2, geotop::input::gDoubleNoValue);
    }

  //    lg->logf("Channel networks has %ld pixels set to channel",cont);
  // commented out: TODO check why we need it and why the map is not correct..
  // SC23.12.2016
  // if(flag >= 0) write_map(geotop::common::Variables::files[fnet], 1,
  // par->format_out, M, geotop::common::Variables::UV,
  // geotop::input::gDoubleNoValue);

  // Bedrock NOT considered here... TODO: fix this mess..
  // bedrock
  flag = file_exists(fbed);
  //         if(flag == 1){
  //                 IT->bed=read_map(2, files[fbed], land->LC, UV,
  //                 (double)number_novalue);
  //         }else{
  //                 IT->bed=new_doublematrix(top->Z0->nrh, top->Z0->nch);
  //                 initialize_doublematrix(IT->bed, 1.E99);
  //         }
  //         if(flag>=0) write_map(files[fbed], 0, par->format_out, IT->bed, UV,
  //         number_novalue);

  //   if(flag == 1){
  //   meteoio_readMap(string(geotop::common::Variables::files[fbed]),IT->bed) ;
  //  }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void read_optionsfile_point(Par *par,
                            Topo *top,
                            Land *land,
                            Soil *sl,
                            Times * /*times*/,
                            InitTools *IT)
{
  size_t i, r, c;
  long num_lines;
  GeoMatrix<double> Q, P, R, S, T, Z, LU;
  GeoMatrix<short> curv;
  short read_dem, read_lu, read_soil, read_sl, read_as, read_sk, read_bed,
        read_curv, flag, coordinates;
  string temp;
  double min, max;
  FILE *f;
  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  // 4. CALCULATE TOPOGRAPHIC PROPERTIES
  // check if there are point coordinates
  coordinates = 1;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptX] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptY] == geotop::input::gDoubleNoValue)
        coordinates = 0;
    }

  // a. read dem
  read_dem = 0;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptLC] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptSY] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptS] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptA] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptCNS] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptCWE] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptCNwSe] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptCNeSw] == geotop::input::gDoubleNoValue)
        {
          read_dem = 1;
        }
    }
  if (read_dem == 1 && coordinates == 0)
    {
      lg->log(
        "Not possible to read from dem because at least one point has no "
        "coordinates",
        geotop::logger::WARNING);
      read_dem = 0;
    }
  if (read_dem == 1)
    {
      flag = file_exists(fdem);
      if (flag == 1)
        {
          lg->logsf(
            geotop::logger::WARNING, "Dem file %s present",
            geotop::common::Variables::files[fdem].c_str());  // FIXED fdem+1

          meteoio_readDEM(Z);

          Q.resize(Z.getRows(), Z.getCols());
          multipass_topofilter(par->lowpass, Z, Q, geotop::input::gDoubleNoValue,
                               1);
          Z = Q;

        }
      else
        {
          read_dem = 0;
          lg->log("Dem file not present", geotop::logger::WARNING);
        }
    }

  if (read_dem == 1)
    {
      par->r_points.resize(par->chkpt.getRows());
      par->c_points.resize(par->chkpt.getRows());

      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          par->r_points[i] =
            row(par->chkpt[i][ptY], Z.getRows() - 1, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue);
          par->c_points[i] =
            col(par->chkpt[i][ptX], Z.getCols() - 1, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue);
          if ((long)par->chkpt[i][ptZ] == geotop::input::gDoubleNoValue)
            par->chkpt[i][ptZ] = Z[par->r_points[i]][par->c_points[i]];
        }
    }

  // b. read land use
  read_lu = 0;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptLC] == geotop::input::gDoubleNoValue) read_lu = 1;
    }
  if (read_lu == 1 && coordinates == 0) read_lu = 0;
  if (read_lu == 1)
    {
      flag = file_exists(flu);
      if (flag == 1)
        {
          meteoio_readMap(string(geotop::common::Variables::files[flu]),
                          LU);  // HACK: add consitency check in meteoioplugin
        }
      else
        {
          lg->log("Landuse file not present, uniform cover considered",
                  geotop::logger::WARNING);
          if (read_dem == 1)
            {
              copydoublematrix_const(1.0, Z, LU, geotop::input::gDoubleNoValue);
            }
          else
            {
              read_lu = 0;
            }
        }
    }

  if (read_lu == 1)
    {
      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          if ((long)par->chkpt[i][ptLC] == geotop::input::gDoubleNoValue)
            {
              r = row(par->chkpt[i][ptY], LU.getRows() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              c = col(par->chkpt[i][ptX], LU.getCols() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              par->chkpt[i][ptLC] = LU[r][c];
            }
        }
    }

  // c. read soil type
  read_soil = 0;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptSY] == geotop::input::gDoubleNoValue)
        read_soil = 1;
    }
  if (read_soil == 1 && coordinates == 0) read_soil = 0;
  if (read_soil == 1)
    {
      flag = file_exists(fsoil);
      if (flag == 1)
        {
          meteoio_readMap(string(geotop::common::Variables::files[fsoil]),
                          P);  // HACK: add consitency check in meteoioplugin
        }
      else
        {
          lg->log("Soiltype file not present", geotop::logger::WARNING);
          read_soil = 0;
        }
    }
  if (read_soil == 1)
    {
      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          if ((long)par->chkpt[i][ptSY] == geotop::input::gDoubleNoValue)
            {
              r = row(par->chkpt[i][ptY], P.getRows() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              c = col(par->chkpt[i][ptX], P.getCols() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              par->chkpt[i][ptSY] = P[r][c];
            }
        }
    }

  // d. read slope
  read_sl = 0;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptS] == geotop::input::gDoubleNoValue) read_sl = 1;
    }
  if (read_sl == 1 && coordinates == 0) read_sl = 0;
  if (read_sl == 1)
    {
      flag = file_exists(fslp);
      if (flag == 1)
        {
          meteoio_readMap(string(geotop::common::Variables::files[fslp]),
                          P);  // HACK: add consitency check in meteoioplugin
        }
      else
        {
          if (read_dem == 0)
            {
              lg->log("Slopes file not present", geotop::logger::WARNING);
              read_sl = 0;
            }
          else
            {
              Q.resize(Z.getRows(), Z.getCols());
              R.resize(Z.getRows(), Z.getCols());
              find_slope(geotop::common::Variables::UV->U[1],
                         geotop::common::Variables::UV->U[2], Z, Q, R,
                         geotop::input::gDoubleNoValue);
              find_max_slope(Z, Q, R, geotop::input::gDoubleNoValue, P);
              if (flag == 0)
                write_map(geotop::common::Variables::files[fslp], 0, par->format_out,
                          P, geotop::common::Variables::UV,
                          geotop::input::gDoubleNoValue);
            }
        }
    }

  if (read_sl == 1)
    {
      find_min_max(P, geotop::input::gDoubleNoValue, &max, &min);
      lg->writefAll("Slope Min:%12g (%12g deg) Max:%12g (%12g deg) \n",
                    tan(min * GTConst::Pi / 180.), min,
                    tan(max * GTConst::Pi / 180.), max);

      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          if ((long)par->chkpt[i][ptS] == geotop::input::gDoubleNoValue)
            {
              r = row(par->chkpt[i][ptY], P.getRows() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              c = col(par->chkpt[i][ptX], P.getCols() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              par->chkpt[i][ptS] = P[r][c];
            }
        }
    }

  // e. read aspect
  read_as = 0;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptA] == geotop::input::gDoubleNoValue) read_as = 1;
    }
  if (read_as == 1 && coordinates == 0) read_as = 0;
  if (read_as == 1)
    {
      flag = file_exists(fasp);
      if (flag == 1)
        {
          meteoio_readMap(string(geotop::common::Variables::files[fasp]),
                          P);  // HACK: add consitency check in meteoioplugin
        }
      else
        {
          if (read_dem == 0)
            {
              lg->log("Aspect file not present", geotop::logger::WARNING);
              read_as = 0;
            }
          else
            {
              Q.resize(Z.getRows(), Z.getCols());
              R.resize(Z.getRows(), Z.getCols());
              find_slope(geotop::common::Variables::UV->U[1],
                         geotop::common::Variables::UV->U[2], Z, Q, R,
                         geotop::input::gDoubleNoValue);
              find_aspect(Z, Q, R, geotop::input::gDoubleNoValue, P);

              if (flag == 0)
                write_map(geotop::common::Variables::files[fasp], 0, par->format_out,
                          P, geotop::common::Variables::UV,
                          geotop::input::gDoubleNoValue);
            }
        }
    }

  if (read_as == 1)
    {
      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          if ((long)par->chkpt[i][ptA] == geotop::input::gDoubleNoValue)
            {
              r = row(par->chkpt[i][ptY], P.getRows() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              c = col(par->chkpt[i][ptX], P.getCols() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              par->chkpt[i][ptA] = P[r][c];
            }
        }
    }

  // f. sky view factor file
  read_sk = 0;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptSKY] == geotop::input::gDoubleNoValue)
        read_sk = 1;
    }
  if (read_sk == 1 && coordinates == 0) read_sk = 0;
  if (read_sk == 1)
    {
      flag = file_exists(fsky);
      if (flag == 1)
        {
          meteoio_readMap(string(geotop::common::Variables::files[fsky]),
                          P);  // HACK: add consitency check in meteoioplugin
        }
      else
        {
          if (read_dem == 0)
            {
              lg->log("Sky view factor file not present", geotop::logger::WARNING);
              read_sk = 0;
            }
          else
            {
              P.resize(Z.getRows(), Z.getCols());
              curv.resize(Z.getRows(), Z.getCols());
              nablaquadro_mask(Z, curv, geotop::common::Variables::UV->U,
                               geotop::common::Variables::UV->V);
              sky_view_factor(P, 36, geotop::common::Variables::UV, Z, curv,
                              geotop::input::gDoubleNoValue);
              if (flag == 0)
                write_map(geotop::common::Variables::files[fsky], 0, par->format_out,
                          P, geotop::common::Variables::UV,
                          geotop::input::gDoubleNoValue);
            }
        }
    }

  if (read_sk == 1)
    {
      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          if ((long)par->chkpt[i][ptSKY] == geotop::input::gDoubleNoValue)
            {
              r = row(par->chkpt[i][ptY], P.getRows() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              c = col(par->chkpt[i][ptX], P.getCols() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              par->chkpt[i][ptSKY] = P[r][c];
            }
        }
    }

  // f2. bedrock file
  read_bed = 1;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptBED] == geotop::input::gDoubleNoValue)
        read_bed = 1;
    }
  if (read_bed == 1 && coordinates == 0) read_bed = 0;
  if (read_bed == 1)
    {
      flag = file_exists(fbed);
      if (flag == 1)
        {
          meteoio_readMap(string(geotop::common::Variables::files[fbed]), P);
        }
      else
        {
          lg->log("Bedrock depth file not present", geotop::logger::WARNING);
          read_bed = 0;
        }
    }

  if (read_bed == 1)
    {
      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          if ((long)par->chkpt[i][ptBED] == geotop::input::gDoubleNoValue)
            {
              r = row(par->chkpt[i][ptY], P.getRows() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              c = col(par->chkpt[i][ptX], P.getCols() - 1,
                      geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              par->chkpt[i][ptBED] = P[r][c];
            }
        }
    }

  // g.curvature
  read_curv = 0;
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptCNS] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptCWE] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptCNwSe] == geotop::input::gDoubleNoValue ||
          (long)par->chkpt[i][ptCNeSw] == geotop::input::gDoubleNoValue)
        read_curv = 1;
    }
  if (read_curv == 1 && coordinates == 0) read_curv = 0;
  if (read_curv == 1)
    {
      if (read_dem == 0)
        {
          lg->log(
            "Dem file is not present, and therefore it is not possible to "
            "calculate curvature",
            geotop::logger::WARNING);
          read_curv = 0;
        }
      else
        {
          Q.resize(Z.getRows(), Z.getCols());
          P.resize(Z.getRows(), Z.getCols());
          R.resize(Z.getRows(), Z.getCols());
          S.resize(Z.getRows(), Z.getCols());
          T.resize(Z.getRows(), Z.getCols());

          multipass_topofilter(par->lowpass_curvatures, Z, Q,
                               geotop::input::gDoubleNoValue, 1);
          curvature(geotop::common::Variables::UV->U[1],
                    geotop::common::Variables::UV->U[2], Q, P, R, S, T,
                    geotop::input::gDoubleNoValue);

          if (geotop::common::Variables::files[fcurv] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fcurv] + string("N-S");
              write_map(temp, 0, par->format_out, P, geotop::common::Variables::UV,
                        geotop::input::gDoubleNoValue);
              temp = geotop::common::Variables::files[fcurv] + string("W-E");
              write_map(temp, 0, par->format_out, R, geotop::common::Variables::UV,
                        geotop::input::gDoubleNoValue);
              temp = geotop::common::Variables::files[fcurv] + string("NW-SE");
              write_map(temp, 0, par->format_out, S, geotop::common::Variables::UV,
                        geotop::input::gDoubleNoValue);
              temp = geotop::common::Variables::files[fcurv] + string("NE-SW");
              write_map(temp, 0, par->format_out, T, geotop::common::Variables::UV,
                        geotop::input::gDoubleNoValue);
            }

          find_min_max(P, geotop::input::gDoubleNoValue, &max, &min);
          lg->writefAll("Curvature N-S Min:%12g  Max:%12g \n", min, max);

          find_min_max(R, geotop::input::gDoubleNoValue, &max, &min);
          lg->writefAll("Curvature W-E Min:%12g  Max:%12g \n", min, max);

          find_min_max(S, geotop::input::gDoubleNoValue, &max, &min);
          lg->writefAll("Curvature NW-SE Min:%12g  Max:%12g \n", min, max);

          find_min_max(T, geotop::input::gDoubleNoValue, &max, &min);
          lg->writefAll("Curvature NE-SW Min:%12g  Max:%12g \n", min, max);
        }
    }
  if (read_curv == 1)
    {
      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          r = row(par->chkpt[i][ptY], P.getRows() - 1,
                  geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
          c = col(par->chkpt[i][ptX], P.getCols() - 1,
                  geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
          if ((long)par->chkpt[i][ptCNS] == geotop::input::gDoubleNoValue)
            par->chkpt[i][ptCNS] = P[r][c];
          if ((long)par->chkpt[i][ptCWE] == geotop::input::gDoubleNoValue)
            par->chkpt[i][ptCWE] = R[r][c];
          if ((long)par->chkpt[i][ptCNwSe] == geotop::input::gDoubleNoValue)
            par->chkpt[i][ptCNwSe] = S[r][c];
          if ((long)par->chkpt[i][ptCNeSw] == geotop::input::gDoubleNoValue)
            par->chkpt[i][ptCNeSw] = T[r][c];
        }
    }

  // h. no value check
  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      if ((long)par->chkpt[i][ptZ] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptZ] = 0.0;
      if ((long)par->chkpt[i][ptLC] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptLC] = 1.0;
      if ((long)par->chkpt[i][ptSY] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptSY] = 1.0;
      if ((long)par->chkpt[i][ptS] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptS] = 0.0;
      if ((long)par->chkpt[i][ptA] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptA] = 0.0;
      if ((long)par->chkpt[i][ptSKY] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptSKY] = 1.0;
      if ((long)par->chkpt[i][ptCNS] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptCNS] = 0.0;
      if ((long)par->chkpt[i][ptCWE] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptCWE] = 0.0;
      if ((long)par->chkpt[i][ptCNwSe] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptCNwSe] = 0.0;
      if ((long)par->chkpt[i][ptCNeSw] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptCNeSw] = 0.0;
      if ((long)par->chkpt[i][ptDrDEPTH] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptDrDEPTH] = par->DepthFreeSurface;  //[mm]
      if ((long)par->chkpt[i][ptMAXSWE] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptMAXSWE] = 1.E10;  //[mm]
      if ((long)par->chkpt[i][ptLAT] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptLAT] = par->latitude;
      if ((long)par->chkpt[i][ptLON] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptLON] = par->longitude;
      if ((long)par->chkpt[i][ptID] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptID] = (double)i;
      if ((long)par->chkpt[i][ptHOR] == geotop::input::gDoubleNoValue)
        par->chkpt[i][ptHOR] = par->chkpt[i][ptID];
    }

  // i.show results
  lg->writefAll("\nPOINTS:\n");
  lg->writefAll(
    "ID,East[m],North[m],Elevation[masl],LandCoverType,SoilType,Slope[deg],"
    "Aspect[deg],SkyViewFactor[-],CurvatureN-S[1/m],CurvatureW-E[1/"
    "m],CurvatureNW-SE[1/m],CurvatureNE-SW[1/"
    "m],DepthFreeSurface[mm],Hor,maxSWE[mm],Lat[deg],Long[deg]\n");

  for (r = 1; r < par->chkpt.getRows(); r++)
    {
      for (c = 1; c < ptTOT; c++)
        {
          lg->writefAll("%f", par->chkpt[r][c]);
          if (c < ptTOT) { lg->writefAll(","); }
        }
      lg->writefAll("\n");
    }

  // l. set UV
  if (read_dem == 0 && read_lu == 0 && read_soil == 0 && read_sl == 0 &&
      read_as == 0 && read_sk == 0)
    {
      geotop::common::Variables::UV->U.resize(4 + 1);
      geotop::common::Variables::UV->V.resize(2 + 1);
    }
  geotop::common::Variables::UV->U[2] = 1.0;
  geotop::common::Variables::UV->U[1] = 1.0;
  geotop::common::Variables::UV->U[4] = 0.0;
  geotop::common::Variables::UV->U[3] = 0.0;
  geotop::common::Variables::UV->V[2] = geotop::input::gDoubleNoValue;
  if (geotop::common::Variables::UV->V[2] < 0)
    {
      geotop::common::Variables::UV->V[1] = -1.;
    }
  else
    {
      geotop::common::Variables::UV->V[1] = 1.;
    }

  // 5. SET CHECKPOINT
  if (par->state_pixel == 1)
    {
      par->rc.resize(par->chkpt.getRows(), 2 + 1);
      for (i = 1; i < par->chkpt.getRows(); i++)
        {
          par->rc[i][1] = 1;
          par->rc[i][2] = i;
        }
    }

  // 6. SET PROPERTIES
  top->East.resize(1 + 1, par->chkpt.getRows());
  top->North.resize(1 + 1, par->chkpt.getRows());
  top->Z0.resize(1 + 1, par->chkpt.getRows());
  land->LC.resize(1 + 1, par->chkpt.getRows());
  land->delay.resize(1 + 1, par->chkpt.getRows());
  sl->type.resize(1 + 1, par->chkpt.getRows());
  top->slope.resize(1 + 1, par->chkpt.getRows());
  top->aspect.resize(1 + 1, par->chkpt.getRows());
  top->curvature1.resize(1 + 1, par->chkpt.getRows());
  top->curvature2.resize(1 + 1, par->chkpt.getRows());
  top->curvature3.resize(1 + 1, par->chkpt.getRows());
  top->curvature4.resize(1 + 1, par->chkpt.getRows());
  top->sky.resize(1 + 1, par->chkpt.getRows());
  top->pixel_type.resize(1 + 1, par->chkpt.getRows());
  top->BC_counter.resize(1 + 1, par->chkpt.getRows());
  top->BC_DepthFreeSurface.resize(par->chkpt.getRows());
  par->maxSWE.resize(2, par->chkpt.getRows());
  top->horizon_point.resize(1 + 1, par->chkpt.getRows());
  top->dzdE.resize(1 + 1, par->chkpt.getRows());
  top->dzdN.resize(1 + 1, par->chkpt.getRows());
  top->latitude.resize(2, par->chkpt.getRows());
  top->longitude.resize(2, par->chkpt.getRows());
  par->IDpoint.resize(par->chkpt.getRows());

  for (i = 1; i < par->chkpt.getRows(); i++)
    {
      top->East[1][i] = par->chkpt[i][ptX];
      top->North[1][i] = par->chkpt[i][ptY];
      top->Z0[1][i] = par->chkpt[i][ptZ];
      land->LC[1][i] = par->chkpt[i][ptLC];

      if ((long)land->LC[1][i] <= 0)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(
            f,
            "Error:: Point %ld has land cover type <= 0. This is not admitted.\n",
            i);
          fclose(f);
          lg->logsf(geotop::logger::ERROR,
                    "Point %ld has land cover type <= 0. This is not admitted.", i);
          lg->log("Geotop failed. See failing report (15).",
                  geotop::logger::CRITICAL);
          exit(1);
        }

      sl->type[1][i] = (long)par->chkpt[i][ptSY];

      if (sl->type[1][i] <= 0)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(
            f, "Error:: Point %ld has soil type <= 0. This is not admitted.\n", i);
          fclose(f);
          lg->logsf(geotop::logger::ERROR,
                    "Point %ld has soil type <= 0. This is not admitted.", i);
          lg->log("Geotop failed. See failing report (16).",
                  geotop::logger::CRITICAL);
          exit(1);
        }

      top->slope[1][i] = par->chkpt[i][ptS];
      top->aspect[1][i] = par->chkpt[i][ptA];
      top->sky[1][i] = par->chkpt[i][ptSKY];
      top->curvature1[1][i] = par->chkpt[i][ptCNS];
      top->curvature2[1][i] = par->chkpt[i][ptCWE];
      top->curvature3[1][i] = par->chkpt[i][ptCNwSe];
      top->curvature4[1][i] = par->chkpt[i][ptCNeSw];
      top->pixel_type[1][i] = 1;
      top->BC_counter[1][i] = i;
      top->BC_DepthFreeSurface[i] = par->chkpt[i][ptDrDEPTH];
      top->horizon_point[1][i] = (long)par->chkpt[i][ptHOR];
      top->dzdE[1][i] = 0.;
      top->dzdN[1][i] = 0.;
      land->delay[1][i] = 0.;

      if (sl->type[1][i] <= 0)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(
            f, "Error:: Point %ld has horizon type <= 0. This is not admitted.\n",
            i);
          fclose(f);
          lg->logsf(geotop::logger::ERROR,
                    "Point %ld has horizon type <= 0. This is not admitted.", i);
          lg->log("Geotop failed. See failing report (17).",
                  geotop::logger::CRITICAL);
          exit(1);
        }

      par->maxSWE[1][i] = par->chkpt[i][ptMAXSWE];
      top->latitude[1][i] = par->chkpt[i][ptLAT];
      top->longitude[1][i] = par->chkpt[i][ptLON];
      par->IDpoint[i] = (long)par->chkpt[i][ptID];
    }

  // 7. SET PAR these are all scalar from now on..
  par->output_soil[1] = 0.;
  par->output_snow[1] = 0.;
  par->output_glac[1] = 0.;
  par->output_surfenergy[1] = 0.;
  par->output_vegetation[1] = 0.;
  par->output_meteo[1] = 0.;

  par->output_soil_bin = 0;
  par->output_snow_bin = 0;
  par->output_glac_bin = 0;
  par->output_surfenergy_bin = 0;
  par->output_meteo_bin = 0;

  // 8. READ HORIZONS
  // find max top->horizon_point
  top->num_horizon_point = 0;
  for (r = 1; r <= top->horizon_point.getRows() - 1; r++)
    {
      for (c = 1; c <= top->horizon_point.getCols() - 1; c++)
        {
          if (top->horizon_point[r][c] > top->num_horizon_point)
            top->num_horizon_point = top->horizon_point[r][c];
        }
    }
  top->horizon_height =
    (double ** *)malloc(top->num_horizon_point * sizeof(double **));
  top->horizon_numlines = (long *)malloc(top->num_horizon_point * sizeof(long));
  for (i = 1; i <= size_t(top->num_horizon_point); i++)
    {
      c = 0;
      do
        {
          flag = 0;
          if (c < par->chkpt.getRows() - 1)
            {
              if (top->horizon_point[1][c + 1] != long(i)) c++;
            }
          if (c < par->chkpt.getRows() - 1)
            {
              if (top->horizon_point[1][c + 1] != long(i)) flag = 1;
            }

        }
      while (flag == 1 && c < par->chkpt.getRows() - 1);

      if (c < par->chkpt.getRows() - 1)
        {
          top->horizon_height[i - 1] =
            read_horizon(0, i, geotop::common::Variables::files[fhorpoint],
                         IT->horizon_col_names, &num_lines);
          top->horizon_numlines[i - 1] = num_lines;
        }
      else
        {
          top->horizon_height[i - 1] = (double **)malloc(sizeof(double *));
          top->horizon_height[i - 1][0] = (double *)malloc(2. * sizeof(double));
          top->horizon_height[i - 1][0][0] = 0.;
          top->horizon_height[i - 1][0][1] = 0.;
          top->horizon_numlines[i - 1] = 1;
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void set_bedrock(InitTools *IT,
                 Soil *sl,
                 Channel *cnet,
                 Par *par,
                 Topo *top,
                 GeoMatrix<double> & /*LC*/)
{
  GeoMatrix<double> B;
  GeoTensor<double> T;
  GeoVector<double> WT;
  long j, l, r, c, sy, synew;
  size_t i;
  double zlim, z;
  short yes = 0;
  FILE *f;

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  if (!mio::IOUtils::fileExists(string(geotop::common::Variables::files[fbed]) +
                                string(ascii_esri)))
    {
      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
      fprintf(f,
              "Error:: File %s is missing. Please check if you have a bedrock "
              "topography map. If it is not available, remove the file name and "
              "keyword from input file\n",
              geotop::common::Variables::files[fbed + 1].c_str());
      fclose(f);
      lg->logsf(geotop::logger::ERROR,
                "File %s is missing. Please check if you have a bedrock "
                "topography map. If it is not available, remove the file name "
                "and keyword from input file.",
                geotop::common::Variables::files[fbed + 1].c_str());
      lg->log("Geotop failed. See failing report (18).",
              geotop::logger::CRITICAL);
      exit(1);
    }

  lg->logf("A bedrock depth map has been assigned and read from %s",
           geotop::common::Variables::files[fbed].c_str());

  par->bedrock = 1;
  meteoio_readMap(string(geotop::common::Variables::files[fbed]), B);

  // check if bedrock depth is above soil lower border, otherwise we do not need
  // to calculate anything
  z = 0.;
  for (l = 1; l < geotop::common::Variables::Nl + 1; l++)
    {
      z += sl->pa[1][jdz][l];
    }
  for (i = 1; i <= par->total_pixel; i++)
    {
      r = top->rc_cont[i][1];
      c = top->rc_cont[i][2];
      if (B[r][c] < z) yes = 1;
    }

  if (yes == 1)
    {
      // consistency check
      if (IT->init_water_table_depth.size() != sl->pa.getDh())
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Error:: Error in bedrock calculations");
          fclose(f);
          lg->logsf(geotop::logger::DEBUG,
                    "IT->init_water_table_depth.size()=%u sl->pa.getDh()=%u",
                    IT->init_water_table_depth.size(), sl->pa.getDh());
          lg->log("Error in bedrock calculations.", geotop::logger::ERROR);
          lg->log(" Geotop failed. See failing report (19).",
                  geotop::logger::CRITICAL);
          exit(1);
        }

      //  rewrite soil type
      T.resize(sl->pa.getDh(), nsoilprop + 1, geotop::common::Variables::Nl + 1);
      for (i = 1; i < sl->pa.getDh(); i++)
        {
          for (j = 1; j <= nsoilprop; j++)
            {
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  T(i, j, l) = sl->pa(i, j, l);
                }
            }
        }
      sl->pa.resize(par->total_pixel + par->total_channel + 1, nsoilprop + 1,
                    geotop::common::Variables::Nl + 1);

      //  rewrite initial water table depth
      WT.resize(IT->init_water_table_depth.size(), geotop::input::gDoubleNoValue);
      for (size_t k = 1; k < IT->init_water_table_depth.size(); k++)
        {
          WT.data.push_back(IT->init_water_table_depth[k]);
        }
      IT->init_water_table_depth.resize(par->total_pixel + par->total_channel +
                                        1);

      // assign jdz (is needed later)
      for (i = 1; i < sl->pa.getDh(); i++)
        {
          for (l = 1; l <= geotop::common::Variables::Nl; l++)
            {
              sl->pa(i, jdz, l) = T(1, jdz, l);
            }
        }

      for (i = 1; i <= par->total_pixel + par->total_channel; i++)
        {
          if (i <= par->total_pixel)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              sy = sl->type[r][c];
              synew = i;
              sl->type[r][c] = synew;
              z = 0.0;
            }
          else
            {
              r = cnet->r[i - par->total_pixel];
              c = cnet->c[i - par->total_pixel];
              sy = cnet->soil_type[i - par->total_pixel];
              synew = i;
              cnet->soil_type[i - par->total_pixel] = synew;
              z = par->depr_channel;
            }

          IT->init_water_table_depth[synew] = WT[sy - 1];

          zlim = B[r][c];

          for (l = 1; l <= geotop::common::Variables::Nl; l++)
            {
              z += 0.5 * sl->pa(synew, jdz, l);

              if (z <= zlim)
                {
                  for (j = 1; j <= nsoilprop; j++)
                    {
                      sl->pa(synew, j, l) = T(sy, j, l);
                    }

                }
              else
                {
                  for (j = 1; j <= nsoilprop; j++)
                    {
                      sl->pa(synew, j, l) = IT->pa_bed(sy, j, l);
                    }
                }

              z += 0.5 * sl->pa(synew, jdz, l);
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

GeoTensor<double> find_Z_of_any_layer(GeoMatrix<double> &Zsurface,
                                      GeoMatrix<double> &slope,
                                      GeoMatrix<double> &LC,
                                      Soil *sl,
                                      short point)
{
  GeoTensor<double> Z;
  double Zaverage = 0.0, z, cosine;
  size_t l, r, c, n;
  long sy;

  if (point != 1)
    {
      Zaverage = 0.;
      n = 0;
      for (r = 1; r < Zsurface.getRows(); r++)
        {
          for (c = 1; c < Zsurface.getCols(); c++)
            {
              if ((long)LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  n++;
                  Zaverage += Zsurface[r][c];
                }
            }
        }
      Zaverage /= (double)n;
    }

  Z.resize(sl->pa.getCh(), Zsurface.getRows(), Zsurface.getCols(),
           geotop::input::gDoubleNoValue);

  for (r = 1; r < Zsurface.getRows(); r++)
    {
      for (c = 1; c < Zsurface.getCols(); c++)
        {
          if ((long)LC[r][c] != geotop::input::gDoubleNoValue)
            {
              cosine = cos(slope[r][c] * GTConst::Pi / 180.);

              sy = sl->type[r][c];
              if (point != 1)
                {
                  z = 1.E3 * (Zsurface[r][c] - Zaverage);  //[mm]
                }
              else
                {
                  z = 0.;
                }

              l = 0;
              Z[l][r][c] = z;

              do
                {
                  l++;
                  z -= 0.5 * sl->pa(sy, jdz, l) * cosine;
                  Z[l][r][c] = z;
                  z -= 0.5 * sl->pa(sy, jdz, l) * cosine;
                }
              while (l < sl->pa.getCh() - 1);
            }
        }
    }

  return Z;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

// TODO: Add meaningful information about the files checked
short file_exists(short key)
{
  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  // no keyword -> -1
  // keyword but does not exist -> 0
  // keyword and exists -> 1

  if (geotop::common::Variables::files[key] == geotop::input::gStringNoValue)
    {
      lg->logf("File for keyword : %s  not present in file list",
               geotop::common::Variables::filenames[key].c_str());
      return (-1);
    }
  else
    {
      bool is_present = mio::IOUtils::fileExists(
                          string(geotop::common::Variables::files[key]) + string(ascii_esri));

      if (is_present)
        {
          lg->log("File " + geotop::common::Variables::files[key] +
                  " existing in format 3 (.asc) ",
                  geotop::logger::WARNING);
          lg->logf("The File is present: %s\n",
                   geotop::common::Variables::files[key].c_str());
          return (1);
        }
      else
        {
          lg->log("File " + geotop::common::Variables::files[key] + " not existing",
                  geotop::logger::WARNING);
          return (0);
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double peat_thickness(double dist_from_channel)
{
  double D;

  if (dist_from_channel < 45.23)
    {
      D = 10. * (47.383 - 0.928 * dist_from_channel +
                 0.010 * pow(dist_from_channel, 2.));
    }
  else
    {
      D = 10. * 26.406;
    }

  return (D);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initialize_soil_state(SoilState *S, long n, long nl)
{
  S->T.resize(nl + 1, n + 1, 0.);

  S->P.resize(nl + 1, n + 1, 0.);

  S->thi.resize(nl + 1, n + 1, 0.0);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void copy_soil_state(SoilState *from, SoilState *to)
{
  long l, i;
  long nl = from->T.getRows(), n = from->T.getCols();

  for (i = 1; i < n; i++)
    {
      to->P[0][i] = from->P[0][i];
      for (l = 1; l < nl; l++)
        {
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

void initialize_veg_state(StateVeg *V, long n)
{
  V->Tv.resize(n + 1, 0.0);

  V->wsnow.resize(n + 1, 0.0);

  V->wrain.resize(n + 1, 0.);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void copy_veg_state(StateVeg *from, StateVeg *to)
{
  long i, n = from->Tv.size();
  for (i = 1; i < n; i++)
    {
      to->Tv[i] = from->Tv[i];
      to->wrain[i] = from->wrain[i];
      to->wsnow[i] = from->wsnow[i];
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
short fill_GTmeteostations_meta(const double &JDE,
                                mio::IOManager &iomanager,
                                Meteo *met)
{
  mio::Config cfg = iomanager.getConfig();
#ifdef WRF_PLUGIN
  std::string tz = cfg.get("SIM_TIME_ZONE", "Input");
#else
  std::string tz = cfg.get("TIME_ZONE", "Input");
#endif
  double d_tz = atof(tz.c_str());

  mio::Date d1(JDE + mio::Date::Matlab_offset, d_tz);
  std::vector<mio::MeteoData> vec_meteo;
  iomanager.getMeteoData(d1, vec_meteo);

  // mio::MeteoData& tmp = vec_meteo[1]; //varible not used

#ifdef WRF_PLUGIN
  size_t i_sky = tmp.getParameterIndex("SKY");
#endif

  size_t lMeteoStationContainerSize = vec_meteo.size() + 1;

  met->st->E.resize(lMeteoStationContainerSize);
  met->st->N.resize(lMeteoStationContainerSize);
  met->st->lat.resize(lMeteoStationContainerSize);
  met->st->lon.resize(lMeteoStationContainerSize);
  met->st->Z.resize(lMeteoStationContainerSize);
  met->st->sky.resize(lMeteoStationContainerSize);
  met->st->ST.resize(lMeteoStationContainerSize);
  met->st->Vheight.resize(lMeteoStationContainerSize);
  met->st->Theight.resize(lMeteoStationContainerSize);

  for (size_t i = 1; i < lMeteoStationContainerSize; i++)
    {
      mio::MeteoData &tmpmeteo = vec_meteo[i - 1];

      met->st->E[i] = tmpmeteo.meta.position.getEasting();
      met->st->N[i] = tmpmeteo.meta.position.getNorthing();
      met->st->lat[i] = tmpmeteo.meta.position.getLat();
      met->st->lon[i] = tmpmeteo.meta.position.getLon();
      met->st->Z[i] = tmpmeteo.meta.position.getAltitude();
#ifdef WRF_PLUGIN
      met->st->sky[i] = tmpmeteo(i_sky);
#endif
      // met->st->sky[i] = 0.97;
      met->st->ST[i] = d_tz;
      met->st->Vheight[i] = 2;
      met->st->Theight[i] = 5;
    }

  return 1;
}
