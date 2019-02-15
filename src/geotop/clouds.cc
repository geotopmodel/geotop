
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#include "struct.geotop.h"
#include "clouds.h"
#include "constants.h"
#include "meteo.h"
#include "radiation.h"
#include "times.h"
#include "timer.h"

#define filecloud "clouds.txt"

extern long number_novalue, number_absent;
extern const char *WORKING_DIRECTORY;

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//lat and lon in [deg]
short fill_meteo_data_with_cloudiness(double **meteo, long meteolines,
                                      double **horizon, long horizonlines, double lat,
                                      double lon, double ST, double Z, double sky, double SWrefl_surr, long ndivday,
                                      double rotation,
                                      double Lozone, double alpha, double beta, double albedo)
{
  GEOTIMER_SECTION(__func__);

  double *cloudtrans;
  long n;

  //if there are radiation data, and no cloudiness
  if ( (long)meteo[0][iSW] != number_absent
       || ( (long)meteo[0][iSWb] != number_absent
            && (long)meteo[0][iSWd] != number_absent ) )
    {

      cloudtrans = (double *)malloc(meteolines*sizeof(double));
      cloudiness(meteo, meteolines, horizon, horizonlines, lat*GTConst::FromDegToRad, lon*GTConst::FromDegToRad,
                 ST, Z, sky, SWrefl_surr, cloudtrans, ndivday, rotation, Lozone, alpha, beta,
                 albedo);

      for (n=0; n<meteolines; n++)
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

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

void cloudiness(double **meteo, long meteolines, double **horizon,
                long horizonlines, double lat, double lon,
                double ST, double Z, double sky, double SWrefl_surr, double *cloudtrans,
                long ndivday, double rotation,
                double Lozone, double alpha, double beta, double albedo)
{

  double tc, tc0;
  double E0, Et, Delta, height_sun, dir_sun;
  long n00, n0, n1, k, n;
  long *ndiv;
  FILE *f;
  char *temp;

  //file header
  temp = join_strings(WORKING_DIRECTORY,filecloud);
  f = fopen(temp,"w");
  fprintf(f,
          "Date,SolarHeight[deg],SolarAzimuth[deg],SinSolarHeight,SWinMeasured[W/m2],SWinClearSky[W/m2],AtmTransmissivity,CloudTransmissivity\n");
  fclose(f);
  free(temp);

  //average cloudiness in ndivday daily intervals
  n = ndivday + 1;
  ndiv = (long *)malloc(n*sizeof(long));

  //initialization
  n00 = 0;
  tc = (double)number_novalue;

  //loop
  do
    {
      tc0=tc;
      find_sunset(n00, &n0, &n1, meteo, meteolines, horizon, horizonlines, lat, lon,
                  ST, rotation);

      ndiv[0]=n0;
      for (k=1; k<=ndivday-1; k++)
        {
          ndiv[k]=(long)(n0+k*(n1-n0)/(double)ndivday);
        }
      ndiv[ndivday]=n1;

      for (k=1; k<=ndivday; k++)
        {
          tc = average_cloudiness(ndiv[k-1], ndiv[k], meteo, meteolines, lat, lon, ST,
                                  Z, sky, SWrefl_surr, rotation, Lozone, alpha, beta, albedo);

          //cloudiness at night (from n00<=n<n0)
          if (k==1)
            {
              for (n=n00; n<n0; n++)
                {
                  if ( (long)tc0 != number_novalue && (long)tc != number_novalue )
                    {
                      cloudtrans[n] = tc0 + (tc-tc0) * (n-n00) / (n0-n00);
                    }
                  else
                    {
                      cloudtrans[n] = (double)number_novalue;
                    }
                  printf("n = %ld/%ld\n",n+1,meteolines);
                }
            }

          //cloudiness during the day
          for (n=ndiv[k-1]; n<ndiv[k]; n++)
            {
              cloudtrans[n] = tc;
              printf("n = %ld/%ld\n",n+1,meteolines);
            }
        }

      n00 = n1;

    }
  while (n00 < meteolines-1);

  n = n00;
  sun(meteo[n][iJDfrom0], &E0, &Et, &Delta);
  height_sun = SolarHeight(meteo[n][iJDfrom0], lat, Delta,
                           (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega);
  dir_sun = SolarAzimuth(meteo[n][iJDfrom0], lat, Delta,
                         (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega) + rotation*GTConst::FromDegToRad;
  if (dir_sun < 0) dir_sun += 2*GTConst::Pi;
  if (dir_sun > 2*GTConst::Pi) dir_sun -= 2*GTConst::Pi;

  if ( shadows_point(horizon, horizonlines, height_sun*180./GTConst::Pi, dir_sun*180./GTConst::Pi,
                     GTConst::Tol_h_mount, GTConst::Tol_h_flat) == 0)
    {
      tc = find_cloudiness(n, meteo, meteolines, lat, lon, ST, Z, sky, SWrefl_surr,
                           rotation, Lozone, alpha, beta, albedo);
    }
  else
    {
      tc = (double)number_novalue;
    }
  cloudtrans[n] = tc;
  printf("n = %ld/%ld\n",n+1,meteolines);

}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

double find_cloudiness(long n, double **meteo, long meteolines, double lat,
                       double lon, double ST, double Z, double sky, double SWrefl_surr,
                       double rotation, double Lozone, double alpha, double beta, double albedo)
{

  double tau_cloud;
  double E0, Et, Delta;
  double JD, JDbegin, JDend;
  double P, RH, T;
  double height_sun, tau_atm;
  long d, m, y, h, mi;
  FILE *f;
  char *temp;

  //initial and final JD of the time step
  if (n == 0)
    {
      JDbegin = meteo[n][iJDfrom0] - 0.5 * ( meteo[n+1][iJDfrom0] -
                                             meteo[n][iJDfrom0] );
    }
  else
    {
      JDbegin = 0.5 * ( meteo[n-1][iJDfrom0] + meteo[n][iJDfrom0] );
    }
  if (n == meteolines - 1)
    {
      JDend = meteo[n][iJDfrom0] + 0.5 * ( meteo[n][iJDfrom0] - meteo[n
                                           -1][iJDfrom0] );
    }
  else
    {
      JDend = 0.5 * ( meteo[n][iJDfrom0] + meteo[n+1][iJDfrom0] );
    }

  //sun variables
  sun(meteo[n][iJDfrom0], &E0, &Et, &Delta);

  //pressure [mbar]
  P=pressure(Z);

  //relative humidity [-]
  RH=meteo[n][iRh];
  if ((long)RH != number_novalue && (long)RH != number_absent)
    {
      RH/=100.;
    }
  else
    {
      if ( (long)meteo[n][iT] != number_absent
           && (long)meteo[n][iT] != number_novalue
           && (long)meteo[n][iTdew] != number_absent
           && (long)meteo[n][iTdew] != number_novalue)
        {
          RH=RHfromTdew(meteo[n][iT], meteo[n][iTdew], Z);
        }
      else
        {
          RH=0.4;
        }
    }
  if (RH<0.01) RH=0.01;

  //air temperature [C]
  T=meteo[n][iT];
  if ((long)T == number_novalue || (long)T == number_absent) T=0.0;

  //cloudiness transmissivity
  tau_cloud = cloud_transmittance(JDbegin, JDend, lat, Delta,
                                  (lon-ST*GTConst::Pi/12.+Et)/GTConst::omega, RH, T, P, meteo[n][iSWd],
                                  meteo[n][iSWb], meteo[n][iSW], E0, sky, SWrefl_surr, Lozone, alpha, beta,
                                  albedo);

  //plotting
  temp = join_strings(WORKING_DIRECTORY,filecloud);
  f = fopen(temp,"a");
  convert_JDfrom0_JDandYear(meteo[n][iJDfrom0], &JD, &y);
  convert_JDandYear_daymonthhourmin(JD, y, &d, &m, &h, &mi);
  height_sun = SolarHeight(meteo[n][iJDfrom0], lat, Delta,
                           (lon-ST*GTConst::Pi/12.+Et)/GTConst::omega);
  tau_atm = atm_transmittance(height_sun, P, RH, T, Lozone, alpha, beta,
                              albedo);
  fprintf(f,"%02.f/%02.f/%04.f %02.f:%02.f,%f,%f,%f,%f,%f,%f,%f\n",(float)d,
          (float)m, (float)y, (float)h,(float)mi,
          height_sun*180./GTConst::Pi, rotation + (SolarAzimuth(meteo[n][iJDfrom0], lat, Delta,
                                                       (lon-ST*GTConst::Pi/12.+Et)/GTConst::omega)) * 180./GTConst::Pi,
          std::max<double>(sin(height_sun), 0.05), meteo[n][iSW], GTConst::Isc*E0*std::max<double>(sin(height_sun),
              0.05)*tau_atm,tau_atm,tau_cloud);
  fclose(f);
  free(temp);
  return tau_cloud;

}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

double average_cloudiness(long n0, long n1, double **meteo, long meteolines,
                          double lat, double lon, double ST, double Z, double sky, double SWrefl_surr,
                          double rotation, double Lozone, double alpha, double beta, double albedo)
{

  long n;
  double tc, tc_av=0.0;
  short is_novalue=0;

  for (n=n0; n<n1; n++)
    {
      tc = find_cloudiness(n, meteo, meteolines, lat, lon, ST, Z, sky, SWrefl_surr,
                           rotation, Lozone, alpha, beta, albedo);
      if ( (long)tc == number_novalue)
        {
          is_novalue = 1;
        }
      else
        {
          tc_av += tc / ((double)(n1-n0));
        }
    }

  if (is_novalue==1 || n0==n1) tc_av = (double)number_novalue;

  return tc_av;

}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//nist line at which the computation begins
//n0 line at sunrise
//n1 line at sunset
//lat and lon in [rad]
void find_sunset(long nist, long *n0, long *n1, double **meteo,
                 long meteolines, double **horizon, long horizonlines,
                 double lat, double lon, double ST, double rotation)
{

  short shad, shad0, a=0;
  long n=nist;
  double alpha, direction, E0, Et, Delta;

  shad=1;//suppose it is in shadow
  *n0=-1;

  do
    {

      //find if it in shadow or not
      shad0=shad;

      sun( meteo[n][iJDfrom0], &E0, &Et, &Delta );
      alpha = SolarHeight(meteo[n][iJDfrom0], lat, Delta,
                          (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega);
      direction = SolarAzimuth(meteo[n][iJDfrom0], lat, Delta,
                               (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega) + rotation*GTConst::FromDegToRad;
      if (direction < 0) direction += 2*GTConst::Pi;
      if (direction > 2*GTConst::Pi) direction -= 2*GTConst::Pi;

      shad=shadows_point(horizon, horizonlines, alpha*180./GTConst::Pi, direction*180./GTConst::Pi,
                         GTConst::Tol_h_mount, GTConst::Tol_h_flat);

      //from shadow to non-shadow = sunrise
      if (shad0==1 && shad==0) *n0=n;

      //from non-shadow to shadow = sunset
      if (shad0==0 && shad==1) a=1;

      n++;

    }
  while (a==0 && n<meteolines);

  *n1=n-1;
  if (*n0==-1) *n0=(*n1);

}


//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
