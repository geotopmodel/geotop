
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
#include "geotop_common.h"
#include "meteo.h"
#include "inputKeywords.h"
#include "global_logger.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void meteo_distr(long *line,
                 long lineLR,
                 Meteo *met,
                 Water *wat,
                 Topo *top,
                 Par *par,
                 double JD0,
                 double JDbeg,
                 double JDend)
{
  long r, c, year, j;
  size_t i;
  double JD;

  // INTERPOLATION OF METEO VARIABLES
  for (i = 1; i < met->st->Z.size(); i++)
    {
      if (par->linear_interpolation_meteo[i] == 1)
        {
          time_interp_linear(JD0, JDbeg, JDend, met->var[i - 1], met->data[i - 1],
                             met->numlines[i - 1], nmet, iJDfrom0, 1,
                             &(line[i - 1]));
        }
      else
        {
          time_interp_constant(JD0, JDbeg, JDend, met->var[i - 1], met->data[i - 1],
                               met->numlines[i - 1], nmet, iJDfrom0, 1,
                               &(line[i - 1]));
        }
    }

  // LAPSE RATES
  if (par->LRflag == 1)
    {
      time_interp_linear(JD0, JDbeg, JDend, met->LRv, met->LRs, met->LRsnr,
                         nlstot, 0, 0, &lineLR);
    }
  else
    {
      convert_JDfrom0_JDandYear((JDbeg + JDend) / 2., &JD, &year);
      for (i = 0; i < nlstot; i++)
        {
          j = floor(JD * met->LRcnc[i] / (365. + is_leap(year)));
          met->LRv[i] = met->LRc[i][j];
        }
    }

  for (i = 0; i < nlstot; i++)
    {
      if ((long)met->LRv[i] == geotop::input::gDoubleNoValue)
        met->LRv[i] = met->LRd[i];
    }

  // DISTRIBUTION OF METEROLOGICAL VARIABLES FROM MEASUREMENTS IN SOME STATIONS
  met->tau_cloud_av.resize(
    top->Z0.getRows(), top->Z0.getCols(),
    geotop::input::
    gDoubleNoValue);  // initialize_doublematrix(met->tau_cl_av_map,
  // (double)number_novalue);
  Meteodistr(
    geotop::common::Variables::UV->U[2], geotop::common::Variables::UV->U[1],
    top->East, top->North, top->Z0, top->curvature1, top->curvature2,
    top->curvature3, top->curvature4, top->slope, top->aspect, met,
    par->slopewtD, par->curvewtD, par->slopewtI, par->curvewtI, par->Vmin,
    par->RHmin, par->dn, par->iobsint, iT, iTdew, iWsx, iWsy, iWs, iPrecInt,
    itauC, met->Tgrid, met->RHgrid, met->Vgrid, met->Vdir, met->Pgrid,
    wat->PrecTot, met->tau_cloud_av, met->LRv[ilsTa], met->LRv[ilsTdew],
    met->LRv[ilsPrec], par->MaxIncrFactWithElev, par->MinIncrFactWithElev,
    par->dew, par->T_rain, par->T_snow, par->snowcorrfact, par->raincorrfact);

  if (par->en_balance == 0)
    {
      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              wat->Pnet[r][c] =
                par->raincorrfact * wat->PrecTot[r][c] * ((JDend - JDbeg) * 24.) *
                cos(top->slope[r][c] * GTConst::Pi / 180.);  // from [mm/h] to [mm]
              if (par->point_sim == 1)
                wat->Pnet[r][c] *= cos(top->slope[r][c] * GTConst::Pi / 180.);
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double pressure(double Z)
{
  double scale_ht = 8500.0;
  return GTConst::Pa0 * exp(-Z / scale_ht);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double temperature(double Z, double Z0, double T0, double gamma)
{
  return T0 - (Z - Z0) * gamma * 1.E-3;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void part_snow(double prec_total,
               double *prec_rain,
               double *prec_snow,
               double temperature,
               double t_rain,
               double t_snow)

/*  Inputs: prec_total:   matrix of total precipitation
            temperature:  matrix of air temperature
            t_rain:     temperature above wich all precipitation is rain
            t_snow:     temperature below wich all precipitation is snow
    Outputs:prec_snow:    matrix of solid precipitation
            prec_rain:    matrix of liquid precipitation
*/

{
  if (temperature <= t_snow)
    {
      *prec_snow = prec_total;
      *prec_rain = 0.0;
    }
  else if (temperature >= t_rain)
    {
      *prec_snow = 0.0;
      *prec_rain = prec_total;
    }
  else
    {
      *prec_snow = prec_total * (t_rain - temperature) / (t_rain - t_snow);
      *prec_rain = prec_total * (temperature - t_snow) / (t_rain - t_snow);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/*
In the following subroutines:
e: vapour pressure [mbar]
P: atmospheric pressure [mbar]
T: air temperature in [C]
*/

// Saturated vapour pressure
double SatVapPressure(double T, double P)
{
  double A, b, c, e;
  A = 6.1121 * (1.0007 + 3.46E-6 * P);
  b = 17.502;
  c = 240.97;
  e = A * exp(b * T / (c + T));
  if (e < 0.5 * P)    // vapour pressure limited to 0.5Patm
    {
      return e;
    }
  else
    {
      return 0.5 * P;
    }
}

// Finds the saturated vapour pressure and its derivative with respect to
// temperature
void SatVapPressure_2(double *e, double *de_dT, double T, double P)
{
  double A, b, c;
  A = 6.1121 * (1.0007 + 3.46E-6 * P);
  b = 17.502;
  c = 240.97;
  *e = A * exp(b * T / (c + T));
  if (*e < 0.5 * P)    // vapour pressure limited to 0.5Patm
    {
      *de_dT = (*e) * (b / (c + T) - b * T / pow(c + T, 2.0));
    }
  else
    {
      *e = 0.5 * P;
      *de_dT = 0.;
    }
}

// Temperature from saturated water vapour
double TfromSatVapPressure(double e, double P)
{
  double A, b, c;
  A = 6.1121 * (1.0007 + 3.46E-6 * P);
  b = 17.502;
  c = 240.97;
  return c * log(e / A) / (b - log(e / A));
}

double SpecHumidity(double e, double P)
{
  // printf("e:%f P:%f Q:%f\n",e,P,0.622*e/(P-0.378*e));
  return 0.622 * e / (P - 0.378 * e);
}

void SpecHumidity_2(double *Q, double *dQ_dT, double RH, double T, double P)
{
  double e, de_dT, dQ_de;
  SatVapPressure_2(&e, &de_dT, T, P);
  *Q = SpecHumidity(RH * e, P);
  // printf("RH:%f e:%f Q:%f\n",RH,e,*Q);
  dQ_de = 0.622 / (P - 0.378 * e) + 0.235116 * e / pow(P - 0.378 * e, 2.0);
  *dQ_dT = dQ_de * de_dT;
}

double VapPressurefromSpecHumidity(double Q, double P)
{
  return Q * P / (0.378 * Q + 0.622);
}

double Tdew(double T, double RH, double Z)
{
  double e, P;
  P = pressure(Z);
  e = SatVapPressure(T, P);
  return TfromSatVapPressure(e * Fmin(RH, 1.0), P);
}

double RHfromTdew(double T, double Tdew, double Z)
{
  double e, es, P, RH;
  P = pressure(Z);
  e = SatVapPressure(Tdew, P);
  es = SatVapPressure(T, P);
  RH = e / es;
  if (RH > 1) RH = 1.0;
  if (RH < 0) RH = 0.0;
  return RH;
}

double air_density(double T, double Q, double P)    //[kg/m3]
{
  double rho;
  rho = P * 100 / (287.04 * (T + 273.15)) *
        (1 - (Q * P / (0.622 + 0.368 * Q)) / P * (1 - 0.622));
  return (rho);
}

double air_cp(double
              T)    // air specific heat at constant pressure [J/(kg K)]
{
  // (Garrat,1992)
  double cp;
  cp = 1005.00 + (T + 23.15) * (T + 23.15) / 3364.0;
  return (cp);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
