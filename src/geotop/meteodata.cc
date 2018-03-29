
/* STATEMENT:

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)

 Copyright (c), 2016 - GEOtop Foundation

 This file is part of GEOtop 2.1

 GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

 GEOtop 2.1  is distributed as a free software in the hope to create and support
 a community of developers and users that constructively interact. If you just
 use the code, please give feedback to the authors and the community. Any way
 you use the model, may be the most trivial one, is significantly helpful for
 the future development of the GEOtop model. Any feedback will be highly
 appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#include "meteodata.h"
#include "geotop_common.h"
#include "inputKeywords.h"
#include "global_logger.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void time_interp_linear(double t0,
                        double tbeg,
                        double tend,
                        double *out,
                        double **data,
                        long nlines,
                        long ncols,
                        long col_date,
                        short flag,
                        long *istart)

{
  short a0 = 0, abeg = 0, aend = 0;
  long i0, ibeg, iend, i, j;
  double t, add;

  i0 = find_line_data(flag, t0, *istart, data, col_date, nlines, &a0);
  ibeg = find_line_data(flag, tbeg, i0, data, col_date, nlines, &abeg);
  iend = find_line_data(flag, tend, ibeg, data, col_date, nlines, &aend);

#ifdef VERY_VERBOSE
  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  lg->logsf(geotop::logger::TRACE,
            "time_interp_linear istart:%ld ibeg:%ld iend:%ld tbeg:%f tend:%f "
            "abeg:%ld aend:%ld",
            *istart, ibeg, iend, tbeg, tend, abeg, aend);
#endif

  for (i = 0; i < ncols; i++)
    {
      if ((long)data[0][i] == geotop::input::gDoubleAbsent)
        {
          out[i] = (double)geotop::input::gDoubleAbsent;

        }
      else if (abeg == 1 && aend == 1)
        {
          out[i] = 0.;

          if ((long)out[i] != geotop::input::gDoubleNoValue)
            {
              add =
                integrate_meas_linear_beh(flag, tbeg, ibeg + 1, data, i, col_date);
#ifdef VERY_VERBOSE
              lg->logsf(geotop::logger::TRACE, "0:tbeg:%f add:%f", tbeg, add);
#endif
              if ((long)add != geotop::input::gDoubleNoValue)
                {
                  out[i] += add;
                }
              else
                {
                  out[i] = geotop::input::gDoubleNoValue;
                }
            }

          if ((long)out[i] != geotop::input::gDoubleNoValue)
            {
              add =
                -integrate_meas_linear_beh(flag, tend, iend + 1, data, i, col_date);
#ifdef VERY_VERBOSE
              lg->logsf(geotop::logger::TRACE, "end:tend:%f add:%f", tend, add);
#endif
              if ((long)add != geotop::input::gDoubleNoValue)
                {
                  out[i] += add;
                }
              else
                {
                  out[i] = geotop::input::gDoubleNoValue;
                }
            }

          j = ibeg + 1;
          while ((long)out[i] != geotop::input::gDoubleNoValue && j <= iend)
            {
              t = time_in_JDfrom0(flag, j, col_date, data);
              add = integrate_meas_linear_beh(flag, t, j + 1, data, i, col_date);
#ifdef VERY_VERBOSE
              lg->logsf(geotop::logger::TRACE, "j:%ld:t:%f iadd:%f", j, t, add);
#endif
              j++;
              if ((long)add != geotop::input::gDoubleNoValue)
                {
                  out[i] += add;
                }
              else
                {
                  out[i] = geotop::input::gDoubleNoValue;
                }
            }

          if ((long)out[i] != geotop::input::gDoubleNoValue)
            {
              out[i] /= (tend - tbeg);
            }

        }
      else
        {
          out[i] = geotop::input::gDoubleNoValue;
        }
    }

  *istart = i0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void time_interp_constant(double t0,
                          double tbeg,
                          double tend,
                          double *out,
                          double **data,
                          long nlines,
                          long ncols,
                          long col_date,
                          short flag,
                          long *istart)

{
  short a0 = 0, abeg = 0, aend = 0;
  long i0, ibeg, iend, i, j;
  double t, add;

  i0 = find_line_data(flag, t0, *istart, data, col_date, nlines, &a0);
  ibeg = find_line_data(flag, tbeg, i0, data, col_date, nlines, &abeg);
  iend = find_line_data(flag, tend, ibeg, data, col_date, nlines, &aend);

#ifdef VERY_VERBOSE
  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();
  lg->logsf(geotop::logger::TRACE,
            "time_interp_constant istart:%ld ibeg:%ld iend:%ld tbeg:%f tend:%f "
            "abeg:%d aend:%d",
            *istart, ibeg, iend, tbeg, tend, abeg, aend);
#endif

  for (i = 0; i < ncols; i++)
    {
      if ((long)data[0][i] == geotop::input::gDoubleAbsent)
        {
          out[i] = (double)geotop::input::gDoubleAbsent;

        }
      else if (abeg == 1 && aend == 1)
        {
          out[i] = 0.;

          j = ibeg;

          while ((long)out[i] != geotop::input::gDoubleNoValue && j < iend)
            {
              t = time_in_JDfrom0(flag, j + 1, col_date, data);
              add = integrate_meas_constant_beh(flag, t, j + 1, data, i, col_date);
#ifdef VERY_VERBOSE
              lg->logsf(geotop::logger::TRACE, "j:%ld t:%f add:%f va:%f\n", j, t, add,
                        data[j + 1][i]);
#endif
              j++;

              if ((long)add != geotop::input::gDoubleNoValue)
                {
                  out[i] += add;
                }
              else
                {
                  out[i] = geotop::input::gDoubleNoValue;
                }
            }

          if ((long)out[i] != geotop::input::gDoubleNoValue)
            {
              add =
                -integrate_meas_constant_beh(flag, tbeg, ibeg + 1, data, i, col_date);
#ifdef VERY_VERBOSE
              lg->logsf(geotop::logger::TRACE, "ibeg:%ld t:beg%f add:%f", ibeg, tbeg,
                        add);
#endif

              if ((long)add != geotop::input::gDoubleNoValue)
                {
                  out[i] += add;
                }
              else
                {
                  out[i] = geotop::input::gDoubleNoValue;
                }
            }

          if ((long)out[i] != geotop::input::gDoubleNoValue)
            {
              add =
                integrate_meas_constant_beh(flag, tend, iend + 1, data, i, col_date);
#ifdef VERY_VERBOSE
              lg->logsf(geotop::logger::TRACE, ":iend:%ld tend:%f add:%f va:%f\n",
                        iend, tend, add, data[iend + 1][i]);
#endif

              if ((long)add != geotop::input::gDoubleNoValue)
                {
                  out[i] += add;
                }
              else
                {
                  out[i] = geotop::input::gDoubleNoValue;
                }
            }

          if ((long)out[i] != geotop::input::gDoubleNoValue)
            {
              out[i] /= (tend - tbeg);
            }

#ifdef VERY_VERBOSE
          lg->logsf(geotop::logger::TRACE, "i:%ld out:%f", i, out[i]);
#endif

        }
      else
        {
          out[i] = geotop::input::gDoubleNoValue;

#ifdef VERY_VERBOSE
          lg->logsf(geotop::logger::TRACE, "<-i:%ld out:%f", i, out[i]);
#endif
        }
    }

  *istart = i0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void time_no_interp(short flag,
                    long *istart,
                    double *out,
                    double **data,
                    long nlines,
                    long ncols,
                    long col_date,
                    double tbeg)

{
  short abeg;
  long ibeg, i;

  ibeg =
    find_line_data(flag, tbeg + 1.E-5, *istart, data, col_date, nlines, &abeg);

  if (abeg == 3)
    {
      ibeg = nlines - 1;
      abeg = 1;
    }

  for (i = 0; i < ncols; i++)
    {
      if ((long)data[0][i] == geotop::input::gDoubleAbsent)
        {
          out[i] = (double)geotop::input::gDoubleAbsent;

        }
      else if (abeg == 1 &&
               (long)data[ibeg][i] != geotop::input::gDoubleNoValue)
        {
          out[i] = data[ibeg][i];

        }
      else
        {
          out[i] = geotop::input::gDoubleNoValue;
        }
    }

  *istart = ibeg;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double integrate_meas_linear_beh(short flag,
                                 double t,
                                 long i,
                                 double **data,
                                 long col,
                                 long col_date)
{
  double t0, t1, value, res;
  FILE *f;

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  if ((long)data[i][col] != geotop::input::gDoubleNoValue &&
      (long)data[i][col] != geotop::input::gDoubleAbsent &&
      (long)data[i - 1][col] != geotop::input::gDoubleNoValue &&
      (long)data[i - 1][col] != geotop::input::gDoubleAbsent)
    {
      t0 = time_in_JDfrom0(flag, i - 1, col_date, data);
      t1 = time_in_JDfrom0(flag, i, col_date, data);

      if (fabs(t0 - t1) < 1.E-5)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f,
                  "Error:: There are 2 consecutive line in a meteo file with same "
                  "date: t0:%12g(%ld) t1:%12g(%ld) equal",
                  t0, i - 1, t1, i);
          fclose(f);
          lg->logsf(geotop::logger::CRITICAL,
                    "There are 2 consecutive line in a meteo file with same date: "
                    "t0:%12g(%ld) t1:%12g(%ld) equal",
                    t0, i - 1, t1, i);
          exit(1);
        }

      value = ((t - t0) * data[i][col] + (t1 - t) * data[i - 1][col]) / (t1 - t0);

#ifdef VERY_VERBOSE
      lg->logsf(geotop::logger::TRACE,
                "Integrate: t:%f t0:%f t1:%f v0:%f v1:%f v:%f\n", t, t0, t1,
                data[i - 1][col], data[i][col], value);
#endif

      // area of trapezium
      res = 0.5 * (value + data[i][col]) * (t1 - t);

    }
  else
    {
      res = geotop::input::gDoubleNoValue;
    }

  return (res);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double integrate_meas_constant_beh(short flag,
                                   double t,
                                   long i,
                                   double **data,
                                   long col,
                                   long col_date)
{
  double t0, value, res;

  if ((long)data[i][col] != geotop::input::gDoubleNoValue &&
      (long)data[i][col] != geotop::input::gDoubleAbsent)
    {
      t0 = time_in_JDfrom0(flag, i - 1, col_date, data);

      value = data[i][col];

      // area of rectangle
      res = value * (t - t0);

    }
  else
    {
      res = geotop::input::gDoubleNoValue;
    }

  return (res);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long find_line_data(short flag,
                    double t,
                    long ibeg,
                    double **data,
                    long col_date,
                    long nlines,
                    short *a)
{
  long i;
  double t0, t1;

  i = ibeg;

  do
    {
      *a = 0;

      if (i + 1 <= nlines - 1)
        {
          t0 = time_in_JDfrom0(flag, i, col_date, data);
          t1 = time_in_JDfrom0(flag, i + 1, col_date, data);

          if (t0 <= t && t1 >= t)
            {
              *a = 1;
            }
          else if (t0 > t)
            {
              *a = 2;
            }
        }

      if (i + 1 >= nlines - 1 && *a == 0) *a = 3;

      if (*a == 0) i++;

    }
  while (*a == 0);

  return (i);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

// All the times in the following subroutines are in Julian days from year 0

// flag 0: dates read in DDMMYYYYhhmm and converted in JDfrom0
// flag 1: dates in JDfrom0

double time_in_JDfrom0(short flag, long i, long col, double **data)
{
  double t;

  if (flag == 0)
    {
      t = convert_dateeur12_JDfrom0(data[i][col]);
    }
  else
    {
      t = data[i][col];
    }

  return (t);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long find_station(long metvar, long nstat, double **var)
{
  long i = 0;

  while ((long)var[i][metvar] == geotop::input::gDoubleAbsent &&
         i < nstat - 1)
    {
      i++;
    }

  return i;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double **read_horizon(short a,
                      long i,
                      std::string name,
                      std::vector<std::string> ColDescr,
                      long *num_lines)
{
  FILE *f;
  long j;
  double **hor = NULL;
  std::string temp;
  short fileyes;

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  // check is the file exists
  if (name != geotop::input::gStringNoValue)
    {
      fileyes = 1;
    }
  else
    {
      fileyes = -1;
    }
  if (fileyes == 1)
    {
      temp = namefile_i(name, i);
      f = fopen(temp.c_str(), "r");
      if (f == NULL)
        {
          fileyes = 0;
        }
      else
        {
          fclose(f);
        }
    }

  // different cases
  if (fileyes == -1)
    {
      if (a == 0)
        {
          lg->logsf(geotop::logger::WARNING,
                    "No horizon file found for point type #%ld. In this case the "
                    "horizon will be considered always not obscured, i.e. "
                    "shadow=FALSE",
                    i);
        }
      else if (a == 1)
        {
          lg->logsf(geotop::logger::WARNING,
                    "No horizon file found for meteo station #%ld. In this case "
                    "the horizon will be considered always not obscured, i.e. "
                    "shadow=FALSE",
                    i);
        }

      *num_lines = 4;
      hor = (double **)malloc((*num_lines) * sizeof(double *));
      for (j = 0; j < (*num_lines); j++)
        {
          hor[j] = (double *)malloc(2 * sizeof(double));
          hor[j][0] = 45.0 + j * 90.0;
          hor[j][1] = 0.0;
        }

    }
  else if (fileyes == 0)
    {
      if (a == 0)
        {
          lg->logsf(geotop::logger::WARNING,
                    "No horizon file found for point type #%ld. In this case the "
                    "horizon will be considered always not obscured, i.e. "
                    "shadow=FALSE",
                    i);
        }
      else if (a == 1)
        {
          lg->logsf(geotop::logger::WARNING,
                    "No horizon file found for meteo station #%ld. In this case "
                    "the horizon will be considered always not obscured, i.e. "
                    "shadow=FALSE",
                    i);
        }

      temp = namefile_i(name, i);

      f = fopen(temp.c_str(), "w");
      fprintf(f, "! Horizon file for met station or point #%ld \n", i);
      fprintf(f, "! All measures in degrees\n");
      fprintf(f, "\n");
      fprintf(f, "%s,%s\n", ColDescr[0].c_str(), ColDescr[1].c_str());

      *num_lines = 4;
      hor = (double **)malloc((*num_lines) * sizeof(double *));
      for (j = 0; j < (*num_lines); j++)
        {
          hor[j] = (double *)malloc(2 * sizeof(double));
          hor[j][0] = 45.0 + j * 90.0;
          hor[j][1] = 0.0;
          fprintf(f, "%f,%f\n", hor[j][0], hor[j][1]);
        }

      fclose(f);

    }
  else if (fileyes == 1)
    {
      if (a == 0)
        {
          lg->logsf(geotop::logger::WARNING,
                    "Horizon file FOUND for point type #%ld", i);
        }
      else if (a == 1)
        {
          lg->logsf(geotop::logger::WARNING,
                    "Horizon file FOUND for meteo station #%ld", i);
        }

      temp = namefile_i(name, i);
      hor = read_txt_matrix(temp, 33, 44, ColDescr, 2, num_lines);

      if ((long)hor[0][0] == geotop::input::gDoubleAbsent ||
          (long)hor[0][1] == geotop::input::gDoubleAbsent)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f,
                  "Error:: In the file %s the columns %s and/or %s are missing\n",
                  temp.c_str(), ColDescr[0].c_str(), ColDescr[1].c_str());
          fclose(f);
          lg->logsf(geotop::logger::CRITICAL,
                    "In the file %s the columns %s and/or %s are missing\n",
                    temp.c_str(), ColDescr[0].c_str(), ColDescr[1].c_str());
          exit(1);
        }
    }

  return (hor);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fixing_dates(long imeteo,
                   double **data,
                   double ST,
                   double STstat,
                   long nlines,
                   long date12col,
                   long JDfrom0col)
{
  long i;
  FILE *f;

  if ((long)data[0][JDfrom0col] == geotop::input::gDoubleAbsent &&
      (long)data[0][date12col] != geotop::input::gDoubleAbsent)
    {
      for (i = 0; i < nlines; i++)
        {
          // converting in JDfrom0
          data[i][JDfrom0col] = convert_dateeur12_JDfrom0(data[i][date12col]);
          // setting ST
          data[i][JDfrom0col] += (ST - STstat) / 24.;
        }

      return 1;

    }
  else if ((long)data[0][JDfrom0col] != geotop::input::gDoubleAbsent)
    {
      return 0;

    }
  else
    {
      f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
      fprintf(f, "Date and Time not available for meteo station %ld\n", imeteo);
      fclose(f);

      geotop::logger::GlobalLogger *lg =
        geotop::logger::GlobalLogger::getInstance();
      lg->logsf(geotop::logger::CRITICAL,
                "Date and Time not available for meteo station %ld\n", imeteo);

      return -1;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fill_wind_xy(double **data,
                   long nlines,
                   long Wspeed,
                   long Wdir,
                   long Wx,
                   long Wy,
                   std::string HeaderWx,
                   std::string HeaderWy)
{
  long i;

  // if the columns Wspeed and Wdir are present, and the columns Wx and Wy are
  // not present
  if ((long)data[0][Wspeed] != geotop::input::gDoubleAbsent &&
      (long)data[0][Wdir] != geotop::input::gDoubleAbsent &&
      ((long)data[0][Wx] == geotop::input::gDoubleAbsent ||
       (long)data[0][Wy] == geotop::input::gDoubleAbsent))
    {
      for (i = 0; i < nlines; i++)
        {
          if ((long)data[i][Wspeed] != geotop::input::gDoubleNoValue &&
              (long)data[i][Wdir] != geotop::input::gDoubleNoValue)
            {
              data[i][Wx] =
                -data[i][Wspeed] * sin(data[i][Wdir] * GTConst::Pi / 180.);
              data[i][Wy] =
                -data[i][Wspeed] * cos(data[i][Wdir] * GTConst::Pi / 180.);
            }
          else
            {
              data[i][Wx] = geotop::input::gDoubleNoValue;
              data[i][Wy] = geotop::input::gDoubleNoValue;
            }

          if (HeaderWx != geotop::input::gStringNoValue &&
              HeaderWy != geotop::input::gStringNoValue)
            {
              data[i][Wspeed] = (double)geotop::input::gDoubleAbsent;
              data[i][Wdir] = (double)geotop::input::gDoubleAbsent;
            }
        }
      if (HeaderWx != geotop::input::gStringNoValue &&
          HeaderWy != geotop::input::gStringNoValue)
        {
          return 1;
        }
      else
        {
          return 0;
        }
    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fill_wind_dir(double **data,
                    long nlines,
                    long Wspeed,
                    long Wdir,
                    long Wx,
                    long Wy,
                    std::string HeaderWSpeed,
                    std::string HeaderWdir)
{
  long i;
  double a;

  // if the columns Wspeed and Wdir are present, and the columns Wx and Wy are
  // not present
  if ((long)data[0][Wx] != geotop::input::gDoubleAbsent &&
      (long)data[0][Wy] != geotop::input::gDoubleAbsent &&
      ((long)data[0][Wspeed] == geotop::input::gDoubleAbsent ||
       (long)data[0][Wdir] == geotop::input::gDoubleAbsent))
    {
      for (i = 0; i < nlines; i++)
        {
          if ((long)data[i][Wx] != geotop::input::gDoubleNoValue &&
              (long)data[i][Wy] != geotop::input::gDoubleNoValue)
            {
              data[i][Wspeed] = sqrt(pow(data[i][Wx], 2.) + pow(data[i][Wy], 2.));

              if (fabs(data[i][Wy]) < 1.E-10)
                {
                  a = GTConst::Pi / 2.;
                }
              else
                {
                  a = atan(fabs(data[i][Wx] / data[i][Wy]));
                }

              if (data[i][Wx] <= 0 && data[i][Wy] <= 0)
                {
                  data[i][Wdir] = a * 180. / GTConst::Pi;
                }
              else if (data[i][Wx] <= 0 && data[i][Wy] >= 0)
                {
                  data[i][Wdir] = a * 180. / GTConst::Pi + 90.;
                }
              else if (data[i][Wx] >= 0 && data[i][Wy] >= 0)
                {
                  data[i][Wdir] = a * 180. / GTConst::Pi + 180.;
                }
              else
                {
                  data[i][Wdir] = a * 180. / GTConst::Pi + 270.;
                }

            }
          else
            {
              data[i][Wspeed] = geotop::input::gDoubleNoValue;
              data[i][Wdir] = geotop::input::gDoubleNoValue;
            }

          if (HeaderWSpeed != geotop::input::gStringNoValue &&
              HeaderWdir != geotop::input::gStringNoValue)
            {
              data[i][Wx] = (double)geotop::input::gDoubleAbsent;
              data[i][Wy] = (double)geotop::input::gDoubleAbsent;
            }
        }
      if (HeaderWSpeed != geotop::input::gStringNoValue &&
          HeaderWdir != geotop::input::gStringNoValue)
        {
          return 1;
        }
      else
        {
          return 0;
        }
    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fill_Tdew(long imeteo,
                GeoVector<double> &Z,
                double **data,
                long nlines,
                long RH,
                long Tair,
                long Tairdew,
                std::string HeaderTdew,
                double RHmin)
{
  long i;

  if ((long)data[0][RH] != geotop::input::gDoubleAbsent &&
      (long)data[0][Tair] != geotop::input::gDoubleAbsent &&
      (long)data[0][Tairdew] == geotop::input::gDoubleAbsent)
    {
      for (i = 0; i < nlines; i++)
        {
          if ((long)data[i][RH] != geotop::input::gDoubleNoValue &&
              (long)data[i][Tair] != geotop::input::gDoubleNoValue)
            {
              data[i][Tairdew] =
                Tdew(data[i][Tair], Fmax(RHmin, data[i][RH]) / 100., Z[imeteo]);
            }
          else
            {
              data[i][Tairdew] = geotop::input::gDoubleNoValue;
            }

          if (HeaderTdew != geotop::input::gStringNoValue)
            {
              data[i][RH] = (double)geotop::input::gDoubleAbsent;
            }
        }
      if (HeaderTdew != geotop::input::gStringNoValue)
        {
          return 1;
        }
      else
        {
          return 0;
        }
    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fill_RH(long imeteo,
              GeoVector<double> &Z,
              double **data,
              long nlines,
              long RH,
              long Tair,
              long Tairdew,
              std::string HeaderRH)
{
  long i;

  if ((long)data[0][RH] == geotop::input::gDoubleAbsent &&
      (long)data[0][Tair] != geotop::input::gDoubleAbsent &&
      (long)data[0][Tairdew] != geotop::input::gDoubleAbsent)
    {
      for (i = 0; i < nlines; i++)
        {
          if ((long)data[i][Tairdew] != geotop::input::gDoubleNoValue &&
              (long)data[i][Tair] != geotop::input::gDoubleNoValue)
            {
              data[i][RH] =
                100. * RHfromTdew(data[i][Tair], data[i][Tairdew], Z[imeteo]);
            }
          else
            {
              data[i][RH] = geotop::input::gDoubleNoValue;
            }

          if (HeaderRH != geotop::input::gStringNoValue)
            {
              data[i][Tairdew] = (double)geotop::input::gDoubleAbsent;
            }
        }
      if (HeaderRH != geotop::input::gStringNoValue)
        {
          return 1;
        }
      else
        {
          return 0;
        }
    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fill_Pint(double **data,
                long nlines,
                long Prec,
                long PrecInt,
                long JDfrom0,
                std::string HeaderPrecInt)
{
  long i;
#ifdef VERY_VERBOSE
  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  lg->log("fill_Pint", geotop::logger::TRACE);
#endif

  if ((long)data[0][Prec] != geotop::input::gDoubleAbsent &&
      (long)data[0][PrecInt] == geotop::input::gDoubleAbsent)
    {
      data[0][PrecInt] = geotop::input::gDoubleNoValue;

      for (i = 1; i < nlines; i++)
        {
          if ((long)data[i][Prec] != geotop::input::gDoubleNoValue)
            {
              data[i][PrecInt] =
                data[i][Prec] / (data[i][JDfrom0] - data[i - 1][JDfrom0]);  //[mm/d]
              data[i][PrecInt] /= 24.;                                      //[mm/h]
#ifdef VERY_VERBOSE
              lg->logsf(geotop::logger::TRACE, "%ld %f %f", i, data[i][PrecInt],
                        data[i][Prec]);
#endif
            }
          else
            {
              data[i][PrecInt] = geotop::input::gDoubleNoValue;
            }

          if (HeaderPrecInt != geotop::input::gStringNoValue)
            {
              data[i][Prec] = (double)geotop::input::gDoubleAbsent;
            }
        }

      if (HeaderPrecInt != geotop::input::gStringNoValue)
        {
          return 1;
        }
      else
        {
          return 0;
        }

    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void check_times(long imeteo, double **data, long nlines, long JDfrom0)
{
  long i;
  FILE *f;

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  for (i = 1; i < nlines; i++)
    {
      if (data[i][JDfrom0] <= data[i - 1][JDfrom0])
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f,
                  "Error:: Time %12g is before Time %12g of previous line in meteo "
                  "file %ld at line %ld.\n",
                  data[i][JDfrom0], data[i - 1][JDfrom0], imeteo, i);
          fclose(f);
          lg->logsf(geotop::logger::CRITICAL,
                    "Time %12g is before Time %12g of previous line in meteo file "
                    "%ld at line %ld.\n",
                    data[i][JDfrom0], data[i - 1][JDfrom0], imeteo, i);
          exit(1);
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void rewrite_meteo_files(double **meteo,
                         long meteolines,
                         std::vector<std::string> header,
                         std::string name,
                         short added_JD,
                         short added_wind_xy,
                         short added_wind_dir,
                         short added_cloudiness,
                         short added_Tdew,
                         short added_RH,
                         short added_Pint)
{
  std::string newname;
  short first_column, write;
  long i, n, d, m, y, h, mi;
  FILE *f;

  newname = name;
  newname += ".old";
  f = fopen(newname.c_str(), "r");
  if (f == NULL)
    {
      write = 1;
    }
  else
    {
      fclose(f);
      write = 0;
    }

  if (added_cloudiness != 1 && added_wind_xy != 1 && added_wind_dir != 1 &&
      added_JD != 1 && added_Tdew != 1 && added_RH != 1 && added_Pint != 1)
    write = 0;

  if (write == 1)
    {
      rename(name.c_str(), newname.c_str());

      f = fopen(name.c_str(), "w");

      first_column = 1;
      for (i = 0; i < nmet; i++)
        {
          if ((long)meteo[0][i] != geotop::input::gDoubleAbsent &&
              header[i].c_str() != geotop::input::gStringNoValue)
            {
              if (first_column == 0)
                {
                  fprintf(f, ",");
                }
              else
                {
                  first_column = 0;
                }
              fprintf(f, "%s", header[i].c_str());
            }
          else if (i == iDate12)
            {
              if (first_column == 0)
                {
                  fprintf(f, ",");
                }
              else
                {
                  first_column = 0;
                }
              fprintf(f, "%s", header[i].c_str());
            }
        }
      fprintf(f, "\n");

      for (n = 0; n < meteolines; n++)
        {
          first_column = 1;
          for (i = 0; i < nmet; i++)
            {
              if ((long)meteo[0][i] != geotop::input::gDoubleAbsent &&
                  header[i].c_str() != geotop::input::gStringNoValue)
                {
                  if (first_column == 0)
                    {
                      fprintf(f, ",");
                    }
                  else
                    {
                      first_column = 0;
                    }
                  if (i == iDate12)
                    {
                      convert_dateeur12_daymonthyearhourmin(meteo[n][i], &d, &m, &y, &h,
                                                            &mi);
                      fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)d, (float)m,
                              (float)y, (float)h, (float)mi);
                    }
                  else
                    {
                      fprintf(f, "%12g", meteo[n][i]);
                    }
                }
              else if (i == iDate12)
                {
                  if (first_column == 0)
                    {
                      fprintf(f, ",");
                    }
                  else
                    {
                      first_column = 0;
                    }
                  convert_dateeur12_daymonthyearhourmin(
                    convert_JDfrom0_dateeur12(meteo[n][iJDfrom0]), &d, &m, &y, &h, &mi);
                  fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)d, (float)m,
                          (float)y, (float)h, (float)mi);
                }
            }
          fprintf(f, "\n");
        }

      fclose(f);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
