/* STATEMENT:

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 - 31 December 2016

     Copyright (c), 2016 - GEOtop Foundation

    This file is part of GEOtop 2.1

   GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

   GEOtop 2.1 is distributed as a free software in the hope to create and
 support a community of developers and users that constructively interact.

 If you just use the code, please give feedback to the authors and the
 community. Any way you use the model, may be the most trivial one, is
 significantly helpful for the future development of the Geotop model. Any
 feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 _________________________

 Note on meteodistr.c - meteodistr.h

 Basic ideas of the routines distributing
 wind-precipitation-temperature-relative humidity are derived from the Micromet
 Fortran Code by Liston and Elder.

 Reference:
 Liston, Glen E.; Elder, Kelly
 A meteorological distribution system for high-resolution terrestrial modeling
 (MicroMet) Journal of Hydrometeorology. 7(April): 217-234.

 However this code is significantly different from the above mentioned-code.

 */

#include "meteodistr.h"
#include "geotop_common.h"
#include "inputKeywords.h"
#include "global_logger.h"

//***************************************************************************************************************
//***************************************************************************************************************
// Local Declarations

static void get_dn(long nc,
                   long nr,
                   double deltax,
                   double deltay,
                   long nstns,
                   double *dn);

static short get_precipitation(double dE,
                               double dN,
                               GeoMatrix<double> &E,
                               GeoMatrix<double> &N,
                               Meteo *met,
                               long Pcode,
                               long Tcode,
                               long Tdcode,
                               GeoMatrix<double> &prec_grid,
                               double dn,
                               GeoMatrix<double> &topo,
                               short iobsint,
                               double lapse_rate,
                               double max,
                               double min,
                               short dew,
                               double Train,
                               double Tsnow,
                               double snow_corr_factor,
                               double rain_corr_factor);

static void get_pressure(GeoMatrix<double> &topo,
                         GeoMatrix<double> &sfc_pressure,
                         double undef);

static short get_relative_humidity(double dE,
                                   double dN,
                                   GeoMatrix<double> &E,
                                   GeoMatrix<double> &N,
                                   Meteo *met,
                                   long Tdcode,
                                   GeoMatrix<double> &RH_grid,
                                   GeoMatrix<double> &Tair_grid,
                                   double RH_min,
                                   double dn,
                                   GeoMatrix<double> &topo,
                                   short iobsint,
                                   double lapse_rate);

static short get_temperature(double dE,
                             double dN,
                             GeoMatrix<double> &E,
                             GeoMatrix<double> &N,
                             Meteo *met,
                             long Tcode,
                             GeoMatrix<double> &Tair_grid,
                             double dn,
                             GeoMatrix<double> &topo,
                             short iobsint,
                             double lapse_rate);

static short get_wind(double dE,
                      double dN,
                      GeoMatrix<double> &E,
                      GeoMatrix<double> &N,
                      Meteo *met,
                      long ucode,
                      long vcode,
                      long Vscode,
                      GeoMatrix<double> &windspd_grid,
                      GeoMatrix<double> &winddir_grid,
                      GeoMatrix<double> &curvature1,
                      GeoMatrix<double> &curvature2,
                      GeoMatrix<double> &curvature3,
                      GeoMatrix<double> &curvature4,
                      GeoMatrix<double> &slope_az,
                      GeoMatrix<double> &terrain_slope,
                      double slopewtD,
                      double curvewtD,
                      double slopewtI,
                      double curvewtI,
                      double windspd_min,
                      double dn,
                      GeoMatrix<double> &topo,
                      short iobsint);

static short interpolate_meteo(short flag,
                               double dX,
                               double dY,
                               GeoMatrix<double> &Xpoint,
                               GeoMatrix<double> &Ypoint,
                               GeoVector<double> &Xst,
                               GeoVector<double> &Yst,
                               double **value,
                               long metcod,
                               GeoMatrix<double> &grid,
                               double dn0,
                               short iobsint);

static void topo_mod_winds(GeoMatrix<double> &winddir_grid,
                           GeoMatrix<double> &windspd_grid,
                           double slopewtD,
                           double curvewtD,
                           double slopewtI,
                           double curvewtI,
                           GeoMatrix<double> &curvature1,
                           GeoMatrix<double> &curvature2,
                           GeoMatrix<double> &curvature3,
                           GeoMatrix<double> &curvature4,
                           GeoMatrix<double> &slope_az,
                           GeoMatrix<double> &terrain_slope,
                           GeoMatrix<double> &topo,
                           double undef);

static void barnes_oi(short flag,
                      GeoMatrix<double> &xpoint,
                      GeoMatrix<double> &ypoint,
                      GeoVector<double> &xstnall,
                      GeoVector<double> &ystnall,
                      GeoVector<double> &xstn,
                      GeoVector<double> &ystn,
                      GeoVector<double> &var,
                      double dn,
                      GeoMatrix<double> &grid,
                      double **value_station,
                      long metcode);

//***************************************************************************************************************
//***************************************************************************************************************

void Meteodistr(double dE,
                double dN,
                GeoMatrix<double> &E,
                GeoMatrix<double> &N,
                GeoMatrix<double> &topo,
                GeoMatrix<double> &curvature1,
                GeoMatrix<double> &curvature2,
                GeoMatrix<double> &curvature3,
                GeoMatrix<double> &curvature4,
                GeoMatrix<double> &terrain_slope,
                GeoMatrix<double> &slope_az,
                Meteo *met,
                double slopewtD,
                double curvewtD,
                double slopewtI,
                double curvewtI,
                double windspd_min,
                double RH_min,
                double dn,
                short iobsint,
                long Tcode,
                long Tdcode,
                long Vxcode,
                long Vycode,
                long VScode,
                long Pcode,
                long /*cloudcode*/,
                GeoMatrix<double> &Tair_grid,
                GeoMatrix<double> &RH_grid,
                GeoMatrix<double> &windspd_grid,
                GeoMatrix<double> &winddir_grid,
                GeoMatrix<double> &sfc_pressure,
                GeoMatrix<double> &prec_grid,
                GeoMatrix<double> & /*tau_cloud_av*/,
                double T_lapse_rate,
                double Td_lapse_rate,
                double Prec_lapse_rate,
                double maxfactorP,
                double minfactorP,
                short dew,
                double Train,
                double Tsnow,
                double snow_corr_factor,
                double rain_corr_factor)
{
  short ok;

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  ok = get_temperature(dE, dN, E, N, met, Tcode, Tair_grid, dn, topo, iobsint,
                       T_lapse_rate);

  if (ok == 0)
    {
      lg->log(
        "No temperature measurements, used the value of the previous time step",
        geotop::logger::DEBUG);
    }

  ok = get_relative_humidity(dE, dN, E, N, met, Tdcode, RH_grid, Tair_grid,
                             RH_min, dn, topo, iobsint, Td_lapse_rate);

  if (ok == 0)
    {
      lg->log("No RH measurements, used the value of the previous time step",
              geotop::logger::DEBUG);
    }

  ok = get_wind(dE, dN, E, N, met, Vxcode, Vycode, VScode, windspd_grid,
                winddir_grid, curvature1, curvature2, curvature3, curvature4,
                slope_az, terrain_slope, slopewtD, curvewtD, slopewtI, curvewtI,
                windspd_min, dn, topo, iobsint);

  if (ok == 0)
    {
      lg->log("No wind measurements, used the value of the previous time step",
              geotop::logger::DEBUG);
    }
  else if (ok == 1)
    {
      lg->log(
        "No wind direction measurements, used the value of the previous time "
        "step",
        geotop::logger::DEBUG);
    }

  ok = get_precipitation(dE, dN, E, N, met, Pcode, Tcode, Tdcode, prec_grid, dn,
                         topo, iobsint, Prec_lapse_rate, maxfactorP, minfactorP,
                         dew, Train, Tsnow, snow_corr_factor, rain_corr_factor);
  if (ok == 0)
    {
      lg->log("No precipitation measurements, considered it 0.0",
              geotop::logger::DEBUG);
    }

  get_pressure(topo, sfc_pressure, geotop::input::gDoubleNoValue);

  // SE we should create a get_cloudiness, for the moment cloud is a scalar
  // ok = interpolate_meteo(0,  dE, dN, E, N, met->st->E, met->st->N,
  // met->var, cloudcode, tau_cloud_av, dn, iobsint);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static short get_temperature(double dE,
                             double dN,
                             GeoMatrix<double> &E,
                             GeoMatrix<double> &N,
                             Meteo *met,
                             long Tcode,
                             GeoMatrix<double> &Tair_grid,
                             double dn,
                             GeoMatrix<double> &topo,
                             short iobsint,
                             double lapse_rate)
{
  double topo_ref;
  long r, c;
  size_t n;
  short ok;
  long nr = topo.getRows() - 1;
  long nc = topo.getCols() - 1;
  size_t nstation = met->st->Z.size() - 1;
  // Define the topographic reference surface.
  topo_ref = 0.0;

  // Convert the station data to sea level values in [K].
  for (n = 1; n < met->st->Z.size(); n++)
    {
      if ((long)met->var[n - 1][Tcode] != geotop::input::gDoubleAbsent &&
          (long)met->var[n - 1][Tcode] != geotop::input::gDoubleNoValue)
        {
          met->var[n - 1][Tcode] = temperature(topo_ref, met->st->Z[n],
                                               met->var[n - 1][Tcode], lapse_rate) +
                                   GTConst::tk;
        }
    }

  // Use the barnes oi scheme to interpolate the station data to the grid.
  ok = interpolate_meteo(1, dE, dN, E, N, met->st->E, met->st->N, met->var,
                         Tcode, Tair_grid, dn, iobsint);

  if (ok == 0) return ok;

  // Convert these grid values back to the actual gridded elevations [C].
  for (r = 1; r <= nr; r++)    // top->Z0.getRows(),top->Z0.getCols()
    {
      for (c = 1; c <= nc; c++)
        {
          if ((long)topo[r][c] != geotop::input::gDoubleNoValue)
            {
              Tair_grid[r][c] =
                temperature(topo[r][c], topo_ref, Tair_grid[r][c], lapse_rate) -
                GTConst::tk;
            }
        }
    }

  // Convert back the station data [C].
  for (n = 1; n <= nstation; n++)
    {
      if ((long)met->var[n - 1][Tcode] != geotop::input::gDoubleAbsent &&
          (long)met->var[n - 1][Tcode] != geotop::input::gDoubleNoValue)
        {
          met->var[n - 1][Tcode] = temperature(met->st->Z[n], topo_ref,
                                               met->var[n - 1][Tcode], lapse_rate) -
                                   GTConst::tk;
        }
    }

  return 1;
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static short get_relative_humidity(double dE,
                                   double dN,
                                   GeoMatrix<double> &E,
                                   GeoMatrix<double> &N,
                                   Meteo *met,
                                   long Tdcode,
                                   GeoMatrix<double> &RH_grid,
                                   GeoMatrix<double> &Tair_grid,
                                   double RH_min,
                                   double dn,
                                   GeoMatrix<double> &topo,
                                   short iobsint,
                                   double lapse_rate)
{
  // First convert stn relative humidity to dew-point temperature.  Use
  //   the Td lapse rate to take the stn Td to sea level.  Interpolate
  //   the stn Td to the grid.  Use the Td lapse rate to take the sea
  //   level grid to the actual elevations.  Convert each Td to
  //   relative humidity.

  double topo_ref;
  long n, r, c;
  short ok;
  long nr = topo.getRows() - 1;
  long nc = topo.getCols() - 1;
  long nstation = met->st->Z.size() - 1;
  // Define the topographic reference surface.
  topo_ref = 0.0;

  // Convert the station data to sea level values in [K].
  for (n = 1; n <= nstation; n++)
    {
      if ((long)met->var[n - 1][Tdcode] != geotop::input::gDoubleAbsent &&
          (long)met->var[n - 1][Tdcode] != geotop::input::gDoubleNoValue)
        {
          met->var[n - 1][Tdcode] =
            temperature(topo_ref, met->st->Z[n], met->var[n - 1][Tdcode],
                        lapse_rate) +
            GTConst::tk;
        }
    }

  // Use the barnes oi scheme to interpolate the station data to the grid.
  ok = interpolate_meteo(1, dE, dN, E, N, met->st->E, met->st->N, met->var,
                         Tdcode, RH_grid, dn, iobsint);
  if (ok == 0) return 0;

  // Convert these grid values back to the actual gridded elevations, and
  // convert to RH
  for (r = 1; r <= nr; r++)
    {
      for (c = 1; c <= nc; c++)
        {
          if ((long)topo[r][c] != geotop::input::gDoubleNoValue)
            {
              RH_grid[r][c] =
                temperature(topo[r][c], topo_ref, RH_grid[r][c], lapse_rate) -
                GTConst::tk;
              RH_grid[r][c] = RHfromTdew(Tair_grid[r][c], RH_grid[r][c], topo[r][c]);
              if (RH_grid[r][c] < RH_min / 100.) RH_grid[r][c] = RH_min / 100.;
            }
        }
    }

  // Convert back the station data [C].
  for (n = 1; n <= nstation; n++)
    {
      if ((long)met->var[n - 1][Tdcode] != geotop::input::gDoubleAbsent &&
          (long)met->var[n - 1][Tdcode] != geotop::input::gDoubleNoValue)
        {
          met->var[n - 1][Tdcode] =
            temperature(met->st->Z[n], topo_ref, met->var[n - 1][Tdcode],
                        lapse_rate) -
            GTConst::tk;
        }
    }

  return 1;
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static void topo_mod_winds(GeoMatrix<double> &winddir_grid,
                           GeoMatrix<double> &windspd_grid,
                           double slopewtD,
                           double curvewtD,
                           double slopewtI,
                           double curvewtI,
                           GeoMatrix<double> &curvature1,
                           GeoMatrix<double> &curvature2,
                           GeoMatrix<double> &curvature3,
                           GeoMatrix<double> &curvature4,
                           GeoMatrix<double> &slope_az,
                           GeoMatrix<double> &terrain_slope,
                           GeoMatrix<double> &topo,
                           double undef)
{
  long r, c, nc = topo.getCols() - 1, nr = topo.getRows() - 1;
  double deg2rad = GTConst::Pi / 180.0, dirdiff, windwt;

  GeoMatrix<double> wind_slope, wind_curv;

  // Compute the wind modification factor which is a function of topography and
  // wind direction following Liston and Sturm (1998).

  // Compute the slope in the direction of the wind.
  //  wind_slope=new_doublematrix(topo->nrh, topo->nch);
  //  initialize_doublematrix(wind_slope, 0.);
  wind_slope.resize(topo.getRows() + 1, topo.getCols() + 1, 0);

  for (r = 1; r <= nr; r++)
    {
      for (c = 1; c <= nc; c++)
        {
          if (topo[r][c] != undef)
            {
              wind_slope[r][c] =
                tan(deg2rad * terrain_slope[r][c] *
                    cos(deg2rad * (winddir_grid[r][c] - slope_az[r][c])));
              if (wind_slope[r][c] > 1) wind_slope[r][c] = 1.;
              if (wind_slope[r][c] < -1) wind_slope[r][c] = -1.;
            }
        }
    }

  // curvature
  //  wind_curv=new_doublematrix(topo->nrh, topo->nch);
  //  initialize_doublematrix(wind_curv, 0.);
  wind_curv.resize(topo.getRows() + 1, topo.getCols() + 1, 0);
  for (r = 1; r <= nr; r++)
    {
      for (c = 1; c <= nc; c++)
        {
          if (topo[r][c] != undef)
            {
              if ((winddir_grid[r][c] > 360. - 22.5 || winddir_grid[r][c] <= 22.5) ||
                  (winddir_grid[r][c] > 180. - 22.5 &&
                   winddir_grid[r][c] <= 180. + 22.5))
                {
                  wind_curv[r][c] = -curvature1[r][c];

                }
              else if ((winddir_grid[r][c] > 90. - 22.5 &&
                        winddir_grid[r][c] <= 90. + 22.5) ||
                       (winddir_grid[r][c] > 270. - 22.5 &&
                        winddir_grid[r][c] <= 270. + 22.5))
                {
                  wind_curv[r][c] = -curvature2[r][c];

                }
              else if ((winddir_grid[r][c] > 135. - 22.5 &&
                        winddir_grid[r][c] <= 135. + 22.5) ||
                       (winddir_grid[r][c] > 315. - 22.5 &&
                        winddir_grid[r][c] <= 315. + 22.5))
                {
                  wind_curv[r][c] = -curvature3[r][c];

                }
              else if ((winddir_grid[r][c] > 45. - 22.5 &&
                        winddir_grid[r][c] <= 45. + 22.5) ||
                       (winddir_grid[r][c] > 225. - 22.5 &&
                        winddir_grid[r][c] <= 225. + 22.5))
                {
                  wind_curv[r][c] = -curvature4[r][c];
                }
            }
        }
    }

  // Calculate the wind speed and direction adjustments.
  // The weights should be chosen so that they limit the total wind weight to
  // between 0.5 and 1.5 (but this is not required).

  for (r = 1; r <= nr; r++)
    {
      for (c = 1; c <= nc; c++)
        {
          if (topo[r][c] != undef)
            {
              // Generate the terrain-modified wind speed.
              windwt = 1.0;

              if (wind_slope[r][c] < 0)
                {
                  windwt += slopewtD * wind_slope[r][c];
                }
              else
                {
                  windwt += slopewtI * wind_slope[r][c];
                }

              if (wind_curv[r][c] < 0)
                {
                  windwt += curvewtD * wind_curv[r][c];
                }
              else
                {
                  windwt += curvewtI * wind_curv[r][c];
                }

              windspd_grid[r][c] *= windwt;

              // Modify the wind direction according to Ryan (1977).
              dirdiff = slope_az[r][c] - winddir_grid[r][c];
              winddir_grid[r][c] =
                winddir_grid[r][c] - 22.5 * Fmin(fabs(wind_slope[r][c]), 1.) *
                sin(deg2rad * (2.0 * dirdiff));
              if (winddir_grid[r][c] > 360.0)
                {
                  winddir_grid[r][c] = winddir_grid[r][c] - 360.0;
                }
              else if (winddir_grid[r][c] < 0.0)
                {
                  winddir_grid[r][c] = winddir_grid[r][c] + 360.0;
                }
            }
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
static short get_wind(double dE,
                      double dN,
                      GeoMatrix<double> &E,
                      GeoMatrix<double> &N,
                      Meteo *met,
                      long ucode,
                      long vcode,
                      long Vscode,
                      GeoMatrix<double> &windspd_grid,
                      GeoMatrix<double> &winddir_grid,
                      GeoMatrix<double> &curvature1,
                      GeoMatrix<double> &curvature2,
                      GeoMatrix<double> &curvature3,
                      GeoMatrix<double> &curvature4,
                      GeoMatrix<double> &slope_az,
                      GeoMatrix<double> &terrain_slope,
                      double slopewtD,
                      double curvewtD,
                      double slopewtI,
                      double curvewtI,
                      double windspd_min,
                      double dn,
                      GeoMatrix<double> &topo,
                      short iobsint)
{
  // This program takes the station wind speed and direction, converts
  //   them to u and v components, interpolates u and v to a grid,
  //   converts the gridded values to speed and direction, and then
  //   runs a simple wind model that adjusts those speeds and
  //   directions according to topographic slope and curvature
  //   relationships.  The resulting speeds and directions are
  //   converted to u and v components and passed back to the main
  //   program to be written to a file.  (All of the conversion
  //   between u-v and speed-dir is done because of the problems
  //   with interpolating over the 360/0 direction line.)

  GeoMatrix<double> u_grid, v_grid;
  long r, c, nc = topo.getCols() - 1, nr = topo.getRows() - 1;
  double rad2deg = 180.0 / GTConst::Pi;
  short oku, okv, ok;

  // Use the barnes oi scheme to interpolate the station data to
  // the grid.
  // U component.
  u_grid.resize(topo.getRows() + 1, topo.getCols() + 1, 0);
  oku = interpolate_meteo(0, dE, dN, E, N, met->st->E, met->st->N, met->var,
                          ucode, u_grid, dn, iobsint);

  // V component.
  v_grid.resize(topo.getRows() + 1, topo.getCols() + 1, 0);
  okv = interpolate_meteo(0, dE, dN, E, N, met->st->E, met->st->N, met->var,
                          vcode, v_grid, dn, iobsint);

  if (oku == 0 || okv == 0)
    {
      // There is not availability of U and V, just consider wind speed as a
      // scalar
      ok = interpolate_meteo(0, dE, dN, E, N, met->st->E, met->st->N, met->var,
                             Vscode, windspd_grid, dn, iobsint);
      if (ok == 0) { return ok; }

    }
  else
    {
      ok = 2;
      // Convert these u and v components to speed and directions.
      // Inside this loop, the u_grid and/or v_grid values will blow.
      // The value of windsp_grid is fine before the loop, and it explodes after.
      for (r = 1; r <= nr; r++)
        {
          for (c = 1; c <= nc; c++)
            {
              if ((long)topo[r][c] != geotop::input::gDoubleNoValue)
                {
                  winddir_grid[r][c] =
                    270.0 - rad2deg * atan2(v_grid[r][c], u_grid[r][c]);
                  if (winddir_grid[r][c] >= 360.0) winddir_grid[r][c] -= 360.0;
                  windspd_grid[r][c] =
                    sqrt(u_grid[r][c] * u_grid[r][c] + v_grid[r][c] * v_grid[r][c]);
                }
            }
        }

      // Modify the wind speed and direction according to simple wind-topography
      // relationships.
      topo_mod_winds(winddir_grid, windspd_grid, slopewtD, curvewtD, slopewtI,
                     curvewtI, curvature1, curvature2, curvature3, curvature4,
                     slope_az, terrain_slope, topo,
                     geotop::input::gDoubleNoValue);
    }

  // Avoid problems of zero (low) winds (for example, turbulence
  // theory, log wind profile, etc., says that we must have some
  // wind.  Thus, some equations blow up when the wind speed gets
  // very small).
  for (r = 1; r <= nr; r++)
    {
      for (c = 1; c <= nc; c++)
        {
          if ((long)topo[r][c] != geotop::input::gDoubleNoValue)
            {
              if (windspd_grid[r][c] < windspd_min) windspd_grid[r][c] = windspd_min;
            }
        }
    }

  return ok;
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static short get_precipitation(double dE,
                               double dN,
                               GeoMatrix<double> &E,
                               GeoMatrix<double> &N,
                               Meteo *met,
                               long Pcode,
                               long Tcode,
                               long Tdcode,
                               GeoMatrix<double> &prec_grid,
                               double dn,
                               GeoMatrix<double> &topo,
                               short iobsint,
                               double lapse_rate,
                               double max,
                               double min,
                               short dew,
                               double Train,
                               double Tsnow,
                               double snow_corr_factor,
                               double rain_corr_factor)
{
  // Interpolate the observed precipitation values to the grid.  Also
  //   interpolate the station elevations to a reference surface.  Use
  //   a precipitation "lapse rate", or adjustment factor to define
  //   the precipitation on the actual elevation grid.  The reason the
  //   interpolated station elevations are used as the topographic
  //   reference surface (instead of something like sea level), is
  //   because the precipitation adjustment factor is a non-linear
  //   function of elevation difference.

  // The adjustment factor that is used comes from: Thornton, P. E.,
  //   S. W. Running, and M. A. White, 1997: Generating surfaces of
  //   daily meteorological variables over large regions of complex
  //   terrain.  J. Hydrology, 190, 214-251.

  long n, r, c, nc = topo.getCols() - 1, nr = topo.getRows() - 1, cnt = 0;
  long nstation = met->st->Z.size() - 1;
  double alfa, prec, f, rain, snow;
  GeoMatrix<double> topo_ref_grid;
  short ok;

  // Use the barnes oi scheme to interpolate the station elevation data
  //   to the grid, so that it can be used as a topographic reference
  //   surface.
  for (n = 1; n <= nstation; n++)
    {
      if ((long)met->var[n - 1][Pcode] != geotop::input::gDoubleNoValue &&
          (long)met->var[n - 1][Pcode] != geotop::input::gDoubleAbsent)
        {
          prec = met->var[n - 1][Pcode];
          met->var[n - 1][Pcode] = met->st->Z[n];
          met->st->Z[n] = prec;
          cnt++;
        }
    }

  if (cnt == 0)
    {
      return 0;

    }
  else
    {
      // topo_ref_grid=new_doublematrix(nr,nc);
      topo_ref_grid.resize(topo.getRows() + 1, topo.getCols() + 1, 0);
      ok = interpolate_meteo(0, dE, dN, E, N, met->st->E, met->st->N, met->var,
                             Pcode, topo_ref_grid, dn, iobsint);

      for (n = 1; n <= nstation; n++)
        {
          if ((long)met->var[n - 1][Pcode] != geotop::input::gDoubleNoValue &&
              (long)met->var[n - 1][Pcode] != geotop::input::gDoubleAbsent)
            {
              prec = met->st->Z[n];
              met->st->Z[n] = met->var[n - 1][Pcode];
              if (dew == 1)
                {
                  part_snow(prec, &rain, &snow, met->var[n - 1][Tdcode], Train, Tsnow);
                }
              else
                {
                  part_snow(prec, &rain, &snow, met->var[n - 1][Tcode], Train, Tsnow);
                }
              met->var[n - 1][Pcode] =
                rain_corr_factor * rain + snow_corr_factor * snow;
            }
        }

      // Use the barnes oi scheme to interpolate the station data to
      //   the grid.
      ok = interpolate_meteo(0, dE, dN, E, N, met->st->E, met->st->N, met->var,
                             Pcode, prec_grid, dn, iobsint);

      //  Convert the gridded station data to the actual gridded elevations.
      for (c = 1; c <= nc; c++)
        {
          for (r = 1; r <= nr; r++)
            {
              if (topo[r][c] != geotop::input::gDoubleNoValue)
                {
                  alfa = 1.E-3 * lapse_rate * (topo[r][c] - topo_ref_grid[r][c]);
                  if (alfa < -1. + 1.E-6) alfa = -1. + 1.E-6;
                  f = (1.0 - alfa) / (1.0 + alfa);
                  if (f > max) f = max;
                  if (f < min) f = min;
                  prec_grid[r][c] = prec_grid[r][c] * f;
                }
            }
        }

      return ok;
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static void get_pressure(GeoMatrix<double> &topo,
                         GeoMatrix<double> &sfc_pressure,
                         double undef)
{
  long r, c, nc = topo.getCols() - 1, nr = topo.getRows() - 1;

  for (c = 1; c <= nc; c++)
    {
      for (r = 1; r <= nr; r++)
        {
          if (topo[r][c] != undef) { sfc_pressure[r][c] = pressure(topo[r][c]); }
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static short interpolate_meteo(short flag,
                               double dX,
                               double dY,
                               GeoMatrix<double> &Xpoint,
                               GeoMatrix<double> &Ypoint,
                               GeoVector<double> &Xst,
                               GeoVector<double> &Yst,
                               double **value,
                               long metcod,
                               GeoMatrix<double> &grid,
                               double dn0,
                               short iobsint)
{
  long r, c, nstn;
  size_t n;
  double dn;
  long nr = Xpoint.getRows() - 1, nc = Xpoint.getCols() - 1;
  nstn = 0;
  for (n = 1; n < Xst.size(); n++)
    {
      if ((long)value[n - 1][metcod] != geotop::input::gDoubleAbsent &&
          (long)value[n - 1][metcod] != geotop::input::gDoubleNoValue)
        nstn++;
    }

  if (nstn == 0) return 0;

  GeoVector<double> var;
  var.resize(nstn + 1, 0);
  GeoVector<double> xst;
  xst.resize(nstn + 1, 0);
  GeoVector<double> yst;
  yst.resize(nstn + 1, 0);
  nstn = 0;
  for (n = 1; n < Xst.size(); n++)
    {
      if ((long)value[n - 1][metcod] != geotop::input::gDoubleAbsent &&
          (long)value[n - 1][metcod] != geotop::input::gDoubleNoValue)
        {
          nstn++;
          var[nstn] = value[n - 1][metcod];
          xst[nstn] = Xst[n];
          yst[nstn] = Yst[n];
        }
    }

  if (nstn > 1)
    {
      if (iobsint == 1)
        {
          dn = dn0;
        }
      else
        {
          get_dn(nc, nr, dX, dY, nstn, &dn);
        }
      barnes_oi(flag, Xpoint, Ypoint, Xst, Yst, xst, yst, var, dn, grid, value,
                metcod);
    }
  else
    {
      for (r = 1; r <= nr; r++)
        {
          for (c = 1; c <= nc; c++)
            {
              grid[r][c] = var[nstn];
            }
        }
    }

  return nstn;
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static void get_dn(long nc,
                   long nr,
                   double deltax,
                   double deltay,
                   long nstns,
                   double *dn)
{
  double dn_max, dn_min;

  dn_max = pow(deltax * nc * deltay * nr, 0.5) *
           ((1.0 + pow((double)nstns, 0.5)) / ((double)nstns - 1.0));
  dn_min = pow((deltax * nc * deltay * nr) / (double)nstns, 0.5);

  *dn = 0.5 * (dn_min + dn_max);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

static void barnes_oi(short flag,
                      GeoMatrix<double> &xpoint,
                      GeoMatrix<double> &ypoint,
                      GeoVector<double> &xstnall,
                      GeoVector<double> &ystnall,
                      GeoVector<double> &xstn,
                      GeoVector<double> &ystn,
                      GeoVector<double> &var,
                      double dn,
                      GeoMatrix<double> &grid,
                      double **value_station,
                      long metcode)
{
  long r, c, mm, nn;
  long nr = xpoint.getRows() - 1, nc = xpoint.getCols() - 1;
  long nstns = xstn.size() - 1;
  long nstnsall = xstnall.size() - 1;

  double gamma = 0.2;

  double xg, yg;        // temporary x and y coords of current cell
  double w1, w2;        // weights for Gauss-weighted average
  double wtot1, wtot2;  // sum of weights
  double ftot1, ftot2;  // accumulators for values, corrections
  double dsq;           // delx**2 + dely**2
  double xa, ya;        // x, y coords of current station
  double xb, yb;        // x, y coords of current station
  double *dvar;
  dvar = (double *)calloc(nstns + 1, sizeof(double));

  double xkappa_1;  // Gauss constant for first pass
  double xkappa_2;  // Gauss constant for second pass
  double rmax_1;    // maximum scanning radii, for first
  double rmax_2;    // and second passes
  double anum_1;    // numerator, beyond scanning radius,
  double anum_2;    // for first and second passes

  // dvar=new_doublevector(nstns);

  // Compute the first and second pass values of the scaling parameter
  //   and the maximum scanning radius used by the Barnes scheme.
  //   Values above this maximum will use a 1/r**2 weighting.  Here I
  //   have assumed a gamma value of 0.2.

  // First-round values, Eqn (13).
  xkappa_1 = 5.052 * pow(2.0 * dn / GTConst::Pi, 2.0);

  // Define the maximum scanning radius to have weight defined by
  //   wt = 1.0 x 10**(-30) = exp(-rmax_1/xkappa_1)
  // Also scale the 1/r**2 wts so that when dsq = rmax, the wts match.
  rmax_1 = xkappa_1 * 30.0 * log(10.0);
  anum_1 = 1.E-30 * rmax_1;

  // Second-round values, Eqn (4).
  xkappa_2 = gamma * xkappa_1;
  rmax_2 = rmax_1 * gamma;
  anum_2 = 1.E-30 * rmax_2;

  // Scan each input data point and construct estimated error, dvar, at that
  // point.
  for (nn = 1; nn <= nstns; nn++)    // 222
    {

      xa = xstn[nn];
      ya = ystn[nn];
      wtot1 = 0.0;
      ftot1 = 0.0;

      for (mm = 1; mm <= nstns; mm++)    // 111
        {

          xb = xstn[mm];
          yb = ystn[mm];
          dsq = pow(xb - xa, 2.0) + pow(yb - ya, 2.0);

          if (dsq <= rmax_1)
            {
              w1 = exp((-dsq) / xkappa_1);
            }
          else
            {
              // Assume a 1/r**2 weight.
              w1 = anum_1 / dsq;
            }

          wtot1 = wtot1 + w1;
          ftot1 = ftot1 + w1 * var[mm];
        }  // end 111

      if (wtot1 == 0.0) printf("stn wt totals zero\n");

      dvar[nn] = var[nn] - ftot1 / wtot1;

    }  // end 222

  // Grid-prediction loop.  Generate the estimate using first set of
  //   weights, and correct using error estimates, dvar, and second
  //   set of weights.

  for (c = 1; c <= nc; c++)      // 666
    {
      for (r = 1; r <= nr; r++)    // 555
        {

          xg = xpoint[r][c];
          yg = ypoint[r][c];

          // Scan each input data point.
          ftot1 = 0.0;
          wtot1 = 0.0;
          ftot2 = 0.0;
          wtot2 = 0.0;

          for (nn = 1; nn <= nstns; nn++)    // 333
            {

              xa = xstn[nn];
              ya = ystn[nn];
              dsq = pow(xg - xa, 2.0) + pow(yg - ya, 2.0);

              if (dsq <= rmax_2)
                {
                  w1 = exp((-dsq) / xkappa_1);
                  w2 = exp((-dsq) / xkappa_2);
                }
              else if (dsq <= rmax_1)
                {
                  w1 = exp((-dsq) / xkappa_1);
                  w2 = anum_2 / dsq;
                }
              else
                {
                  // Assume a 1/r**2 weight.
                  w1 = anum_1 / dsq;
                  // With anum_2/dsq.
                  w2 = gamma * w1;
                }

              wtot1 = wtot1 + w1;
              wtot2 = wtot2 + w2;
              ftot1 = ftot1 + w1 * var[nn];
              ftot2 = ftot2 + w2 * dvar[nn];
            }  // end 333

          if (wtot1 == 0.0 || wtot2 == 0.0) printf("wts total zero\n");

          grid[r][c] = ftot1 / wtot1 + ftot2 / wtot2;

        }  // end 555
    }    // end 666

  // station points
  if (flag == 1)
    {
      for (c = 1; c <= nstnsall; c++)    // 666
        {

          xg = xstnall[c];
          yg = ystnall[c];

          // Scan each input data point.
          ftot1 = 0.0;
          wtot1 = 0.0;
          ftot2 = 0.0;
          wtot2 = 0.0;

          for (nn = 1; nn <= nstns; nn++)    // 333
            {

              xa = xstn[nn];
              ya = ystn[nn];
              dsq = pow(xg - xa, 2.0) + pow(yg - ya, 2.0);

              if (dsq <= rmax_2)
                {
                  w1 = exp((-dsq) / xkappa_1);
                  w2 = exp((-dsq) / xkappa_2);
                }
              else if (dsq <= rmax_1)
                {
                  w1 = exp((-dsq) / xkappa_1);
                  w2 = anum_2 / dsq;
                }
              else
                {
                  // Assume a 1/r**2 weight.
                  w1 = anum_1 / dsq;
                  // With anum_2/dsq.
                  w2 = gamma * w1;
                }

              wtot1 = wtot1 + w1;
              wtot2 = wtot2 + w2;
              ftot1 = ftot1 + w1 * var[nn];
              ftot2 = ftot2 + w2 * dvar[nn];

            }  // end 333

          if (wtot1 == 0.0 || wtot2 == 0.0) printf("wts total zero\n");

          value_station[c - 1][metcode] = ftot1 / wtot1 + ftot2 / wtot2;

        }  // end 666
    }
  free(dvar);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
