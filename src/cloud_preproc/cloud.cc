/*
 * cloud.cc
 *
 *  Created on: Feb 3, 2014
 *      Author: matteo
 */

//#include <iostream>
#include <fstream>
//#include <string>
//#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <inputKeywords.h>
//#include "clouds.h"
#include "input.h"
#include "parameters.h"
//#include "meteodata.h"
//#include "radiation.h"

using namespace std;
using namespace mio;  // The MeteoIO namespace is called mio
using namespace geotop;

/** @brief This function calculates the average cloudiness of a dataset of meteo
 * variables recorded in a meteo station where the shortwave radiation is
 * measured. The day is considered to be the composed by: daylight (time where
 * the sun is present) and night. The start of the daylight is given by moment
 * (meteo line) in which the sun rises. The end of the daylight corresponds to
 * the moment (meteo line) when the sun sets. The function determines the meteo
 * lines when the sun rises (n0) and when the sun sets (n1), then divides this
 * time frame into ndivday intervals. Then, for each meteo line inside this
 * "daylight" interval, calculates the cloud_transmittance and then makes the
 * average inside the interval. For the cloudiness at night, linearly
 * interpolates between the previous sunset and the sun rise.
 *
 * @param[in] **meteo a pointer to a previously initialized configuration store
 * @param[in] **horizon
 * @param[in] long horizonlines
 * @param[in] double lat
 * @param[in] double lon
 * @param[in] double ST
 * @param[in] double Z
 * @param[in] double sky
 * @param[in] double SWrefl_surr
 * @param[in] double *cloudtrans
 * @param[in] long ndivday: number of intervals in which the daylight is divided
 * @param[in] double rotation
 * @param[in] double Lozone
 * @param[in] double alpha
 * @param[in] double beta
 * @param[in] double albedo
 * @return a vector of average cloud_transmissivity (*cloudtrans)
 */
void cloudiness_MeteoIO(double **meteo,
                        long meteolines,
                        double **horizon,
                        long horizonlines,
                        double lat,
                        double lon,
                        double ST,
                        double Z,
                        double sky,
                        double SWrefl_surr,
                        double *cloudtrans,
                        long ndivday,
                        double rotation,
                        double Lozone,
                        double alpha,
                        double beta,
                        double albedo)
{
  double tc;   // tau cloud
  double tc0;  // tau cloud at the previous time step
  double E0;
  double Et;
  double Delta;
  double height_sun;  // solar elevation angle
  double direction;   // solar azimuth angle
  long n00;           // line at which the computation begins
  long n0;            // line of the meteo file at which the sun rises
  long n1;            // line of the meteo file at which the sun sets
  long k, n;          // indexes
  long *ndiv;
  FILE *f;
  //  char *temp;
  std::string temp;

  // file header
  temp = geotop::common::Variables::WORKING_DIRECTORY + filecloud;
  f = fopen(temp.c_str(), "w");
  fprintf(f,
          "Date,SolarHeight[deg],SolarAzimuth[deg],SinSolarHeight,SWinMeasured["
          "W/m2],SWinClearSky[W/m2],AtmTransmissivity,CloudTransmissivity\n");
  fclose(f);
  //  free(temp);

  // average cloudiness in ndivday daily intervals
  n = ndivday + 1;
  ndiv = (long *)malloc(n * sizeof(long));

  // initialization
  n00 = 0;
  tc = geotop::input::gDoubleNoValue;

  // loop
  do
    {
      tc0 = tc;
      find_sunset(n00, &n0, &n1, meteo, meteolines, horizon, horizonlines, lat,
                  lon, ST, rotation);

      ndiv[0] = n0;  // line of the meteo data at which the sun rises
      for (k = 1; k <= ndivday - 1; k++)
        {
          ndiv[k] = (long)(n0 + k * (n1 - n0) / (double)ndivday);
        }
      ndiv[ndivday] = n1;  // line of the meteo data at which the sun sets

      for (k = 1; k <= ndivday; k++)
        {
          // find the average cloudiness during the daylight (between the line
          // ndiv[k-1] and ndiv[k])
          tc = average_cloudiness(ndiv[k - 1], ndiv[k], meteo, meteolines, lat, lon,
                                  ST, Z, sky, SWrefl_surr, rotation, Lozone, alpha,
                                  beta, albedo);

          // cloudiness at night (from n00<=n<n0)
          if (k == 1)
            {
              for (n = n00; n < n0; n++)
                {
                  if ((long)tc0 != geotop::input::gDoubleNoValue &&
                      (long)tc != geotop::input::gDoubleNoValue)
                    {
                      cloudtrans[n] = tc0 + (tc - tc0) * (n - n00) / (n0 - n00);
                    }
                  else
                    {
                      cloudtrans[n] = geotop::input::gDoubleNoValue;
                    }
                  printf("n = %ld/%ld\n", n + 1, meteolines);
                }
            }

          // assigns the values to the vector
          for (n = ndiv[k - 1]; n < ndiv[k]; n++)
            {
              cloudtrans[n] = tc;
              printf("n = %ld/%ld\n", n + 1, meteolines);
            }
        }

      n00 = n1;  // set the computation line at the line where the sun sets

    }
  while (n00 < meteolines - 1);

  n = n00;
  sun(meteo[n][iJDfrom0], &E0, &Et, &Delta);
  height_sun =
    SolarHeight(meteo[n][iJDfrom0], lat, Delta,
                (lon - ST * GTConst::Pi / 12. + Et) / GTConst::omega);
  direction =
    SolarAzimuth(meteo[n][iJDfrom0], lat, Delta,
                 (lon - ST * GTConst::Pi / 12. + Et) / GTConst::omega) +
    rotation * GTConst::Pi / 180.;
  if (direction < 0) direction += 2 * GTConst::Pi;
  if (direction > 2 * GTConst::Pi) direction -= 2 * GTConst::Pi;

  if (shadows_point(horizon, horizonlines, height_sun * 180. / GTConst::Pi,
                    direction * 180. / GTConst::Pi, GTConst::Tol_h_mount,
                    GTConst::Tol_h_flat) == 0)
    {
      tc = find_cloudiness(n, meteo, meteolines, lat, lon, ST, Z, sky,
                           SWrefl_surr, rotation, Lozone, alpha, beta, albedo);
    }
  else
    {
      tc = geotop::input::gDoubleNoValue;
    }
  cloudtrans[n] = tc;
  printf("n = %ld/%ld\n", n + 1, meteolines);
}

short fill_meteo_data_with_cloudiness_MeteoIO(double **meteo,
                                              long meteolines,
                                              double **horizon,
                                              long horizonlines,
                                              double lat,
                                              double lon,
                                              double ST,
                                              double Z,
                                              double sky,
                                              double SWrefl_surr,
                                              long ndivday,
                                              double rotation,
                                              double Lozone,
                                              double alpha,
                                              double beta,
                                              double albedo)
{
  double *cloudtrans;
  long n;

  // if there are radiation data, and no cloudiness
  if ((long)meteo[0][iSW] != geotop::input::gDoubleAbsent ||
      ((long)meteo[0][iSWb] != geotop::input::gDoubleAbsent &&
       (long)meteo[0][iSWd] != geotop::input::gDoubleAbsent))
    {
      cloudtrans = (double *)malloc(meteolines * sizeof(double));
      cloudiness_MeteoIO(meteo, meteolines, horizon, horizonlines,
                         lat * GTConst::Pi / 180., lon * GTConst::Pi / 180., ST,
                         Z, sky, SWrefl_surr, cloudtrans, ndivday, rotation,
                         Lozone, alpha, beta, albedo);

      for (n = 0; n < meteolines; n++)
        {
          meteo[n][itauC] = cloudtrans[n];
        }

      free(cloudtrans);
      return 1;

    }
  else
    {
      return 0;
    }
}

//******************************************************************************
//******************************************************************************
//******************************************************************************
int main(int argc, char **argv)
{
  AllData *adt;
  InitTools *IT;
  string temp;
  adt = new AllData();

  adt->I = new Times();
  if (!(adt->I)) t_error("times was not allocated");

  adt->S = new Soil();
  if (!(adt->S)) t_error("sl was not allocated");

  adt->L = new Land();
  if (!(adt->L)) t_error("land was not allocated");

  adt->M = new Meteo();
  if (!(adt->M)) t_error("met was not allocated");

  adt->P = new Par();
  if (!(adt->P)) t_error("par was not allocated");

  IT = new InitTools();

  short success, added_cloud = 0, added_wind_xy = 0, added_wind_dir = 0,
                 added_Tdew = 0, added_RH = 0, added_Pint = 0;
  long ist, num_lines, num_cols, day, month, year, hour, minute;
  double JD;
  FILE *flog, *f;

  //  string cfgfile = "io_it.ini";
  //  std::string lDataPath ;
  //  if (argc >= 2)
  //  {
  //    lDataPath = argv[1] ;
  //  } else {
  //    lDataPath = get_workingdirectory() ;
  //  }
  //
  //  if(lDataPath == "" )
  //  {
  //    std::cerr << "Error: data path is empty" << std::endl ;
  //    exit (200) ;
  //  }
  //
  //  chdir(lDataPath.c_str());
  //  char lCWD[8192] ;
  //  char * lCwdStr = getcwd(lCWD, sizeof(lCWD));
  //    if (lCWD == NULL){
  //      std::cerr << "Error: unable to get the current path: " <<
  //    strerror(errno) << std::endl ;  exit (201);
  //    }
  //    std::string lFullPath(lCwdStr) ;
  //  geotop::common::Variables::WORKING_DIRECTORY = lFullPath ;
  //  cfgfile = geotop::common::Variables::WORKING_DIRECTORY + "/" + cfgfile;
  //
  //  mio::Config cfg(cfgfile);
  //  cfg.addKey("GRID2DPATH", "Input", "");
  //  mio::IOManager iomanager(cfg);

  if (geotop::common::Variables::WORKING_DIRECTORY != "")
    {
    }
  else if (!argv[1])
    {
      geotop::common::Variables::WORKING_DIRECTORY = get_workingdirectory();
    }
  else if (argc == 2)
    {
      // modified by Emanuele Cordano on Aug 2011
      geotop::common::Variables::WORKING_DIRECTORY = argv[1];
    }
  else
    {
      // modified by Emanuele Cordano on Aug 2011
    }

  // add "/" if it is missing
  if (geotop::common::Variables::WORKING_DIRECTORY
      [strlen(geotop::common::Variables::WORKING_DIRECTORY.c_str()) - 1] !=
      47)
    {
      temp = geotop::common::Variables::WORKING_DIRECTORY;
      geotop::common::Variables::WORKING_DIRECTORY = temp + "/";
    }

  geotop::common::Variables::logfile =
    geotop::common::Variables::WORKING_DIRECTORY + logfile_name;
  flog = fopen(geotop::common::Variables::logfile.c_str(), "w");
  temp = geotop::common::Variables::WORKING_DIRECTORY + program_name;

  std::shared_ptr<geotop::input::ConfigStore> lConfigStore =
    geotop::input::ConfigStoreSingletonFactory::getInstance();
  const std::string lFilePath(temp);
  bool lParsingRes = lConfigStore->parse(lFilePath);
  if (not lParsingRes)
    {
      t_error("Fatal Error! Unable to parse configuration file: " + lFilePath);
    }

  geotop::common::Variables::i_sim0 = 1;
  geotop::common::Variables::cum_time = 0.;
  geotop::common::Variables::elapsed_time_start = 0.;

  success = read_inpts_par(adt->P, adt->L, adt->I, adt->S, adt->M, IT, flog);

  // Time indices
  adt->P->init_date[geotop::common::Variables::i_sim0] +=
    adt->P->delay_day_recover;
  convert_JDfrom0_JDandYear(
    adt->P->init_date[geotop::common::Variables::i_sim0], &JD, &year);
  convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);

  num_cols = (long)nmet;
  adt->M->horizon = (double ** *)malloc(adt->M->st->E.size() * sizeof(
                                          double **));
  adt->M->horizonlines = (long *)malloc(adt->M->st->E.size() * sizeof(long));
  adt->M->data = (double ** *)malloc(adt->M->st->E.size() * sizeof(double **));
  adt->M->var = (double **)malloc((adt->M->st->E.size() - 1) * sizeof(
                                    double *));
  adt->M->numlines = (long *)malloc(adt->M->st->E.size() * sizeof(long));
  success =
    read_meteostations_file(adt->M->imeteo_stations, adt->M->st,
                            geotop::common::Variables::files[fmetstlist],
                            IT->meteostations_col_names, flog);
  long num_met_stat = adt->M->st->E.size() - 1;

  for (size_t i = 1; i <= num_met_stat; i++)
    {
      //  if (met->imeteo_stations->co[1] != geotop::input::gDoubleNoValue) {
      if (adt->M->imeteo_stations[1] != geotop::input::gDoubleNoValue)
        {
          //  ist = met->imeteo_stations->co[i];
          ist = adt->M->imeteo_stations[i];
        }
      else
        {
          ist = i;
        }

      //  read horizon
      adt->M->horizon[i - 1] =
        read_horizon(1, ist, geotop::common::Variables::files[fhormet],
                     IT->horizon_col_names, &num_lines, flog);
      adt->M->horizonlines[i - 1] = num_lines;

      adt->M->var[i - 1] = (double *)malloc(num_cols * sizeof(double));

      if (geotop::common::Variables::files[fmet] !=
          geotop::input::gStringNoValue)
        {
          // read matrix
          temp = namefile_i(geotop::common::Variables::files[fmet], ist);

          adt->M->data[i - 1] = read_txt_matrix(temp, 33, 44, IT->met_col_names,
                                                nmet, &num_lines, flog);

          if ((long)adt->M->data[i - 1][0][iDate12] ==
              geotop::input::gDoubleAbsent &&
              (long)adt->M->data[i - 1][0][iJDfrom0] ==
              geotop::input::gDoubleAbsent)
            {
              f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
              fprintf(f, "Error:: Date Column missing in file %s\n", temp.c_str());
              fclose(f);
              t_error("Fatal Error! Geotop is closed. See failing report (2).");
            }
          adt->M->numlines[i - 1] = num_lines;

          // fixing dates: converting times in the same standard time set for the
          // simulation and fill JDfrom0
          short added_JDfrom0 =
            fixing_dates(ist, adt->M->data[i - 1], adt->P->ST, adt->M->st->ST[i],
                         adt->M->numlines[i - 1], iDate12, iJDfrom0);

          check_times(ist, adt->M->data[i - 1], adt->M->numlines[i - 1], iJDfrom0);

          // find clouds
          if (IT->met_col_names[itauC] != geotop::input::gStringNoValue)
            {
              if ((long)adt->M->data[i - 1][0][itauC] ==
                  geotop::input::gDoubleAbsent ||
                  adt->P->ric_cloud == 1)
                {
                  added_cloud = fill_meteo_data_with_cloudiness(
                                  adt->M->data[i - 1], adt->M->numlines[i - 1],
                                  adt->M->horizon[i - 1], adt->M->horizonlines[i - 1],
                                  adt->M->st->lat[i], adt->M->st->lon[i], adt->P->ST,
                                  adt->M->st->Z[i], adt->M->st->sky[i], 0.0, adt->P->ndivdaycloud,
                                  adt->P->dem_rotation, adt->P->Lozone, adt->P->alpha_iqbal,
                                  adt->P->beta_iqbal, 0.);
                }
            }

          rewrite_meteo_files(adt->M->data[i - 1], adt->M->numlines[i - 1],
                              IT->met_col_names, temp.c_str(), added_JDfrom0,
                              added_wind_xy, added_wind_dir, added_cloud,
                              added_Tdew, added_RH, added_Pint);
          // IT->met_col_name;
        }
    }
  fclose(flog);

  return 0;
}
