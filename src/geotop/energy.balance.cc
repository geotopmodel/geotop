
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

// @author: Stefano Endrizzi
// @brief Energy balance (and also mass balance for snow and glacier)
#include "energy.balance.h"
#include "geotop_common.h"
#include "inputKeywords.h"
#include <assert.h>
#include "global_logger.h"

using namespace std;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short EnergyBalance(double Dt,
                    double JD0,
                    double JDb,
                    double JDe,
                    SoilState *L,
                    SoilState *C,
                    Statevar3D *S,
                    Statevar3D *G,
                    StateVeg *V,
                    GeoVector<double> &snowage,
                    AllData *A,
                    double *W,
                    std::vector<mio::MeteoData> &vec_meteo)
{
  (void)vec_meteo;  // silence warning when we do not use METEIO
  short sux;
  size_t i = 0, r, c, cnt = 0;
  static long line_interp;
  double Dtplot = 1., Delta, E0, Et, SWup, Tgskin, SWrefl_surr_ave = 0.,
         Tgskin_surr_ave = 0.;
  double tau_cloud, tau_cloud_av;
  short tau_cloud_yes, tau_cloud_av_yes;
  FILE *f = NULL;

  //  calculation to be done before plotting maps
  *W = 0.;
  if (A->I->time == 0) A->I->iplot = 1;
  assert(A->I->iplot >= 0);
  if (A->I->JD_plots.size() > 1 &&
      (size_t)(A->I->iplot) <= A->I->JD_plots.size())
    {
      i = 2 * A->I->iplot - 1;

      if (A->P->init_date + A->I->time / GTConst::secinday + 1.E-5 >=
          A->I->JD_plots[i] &&
          A->P->init_date + (A->I->time + Dt) / GTConst::secinday - 1.E-5 <=
          A->I->JD_plots[i + 1])
        {
          Dtplot = (A->I->JD_plots[i + 1] - A->I->JD_plots[i]) * GTConst::secinday;
          *W = Dt / Dtplot;
          f = fopen(geotop::common::Variables::logfile.c_str(), "a");

          fprintf(f, "Saving plot number %ld Weight:%f \n", A->I->iplot, *W);
          printf("Saving plot number %ld Weight:%f \n", A->I->iplot, *W);

          fclose(f);
        }
    }

  //  vegetation

  if (A->I->time == 0) line_interp = 0;
  for (size_t j = 1; j <= A->P->n_landuses; j++)
    {
      if (A->P->vegflag[j] == 1)
        {
          time_interp_linear(
            JD0, JDb, JDe, A->L->vegparv[j - 1], A->L->vegpars[j - 1],
            A->L->NumlinesVegTimeDepData[j - 1], jdvegprop + 1, 0, 0, &line_interp);
        }
    }

  //  SUN and SHADOWING
  sun((JDb + JDe) / 2., &E0, &Et, &Delta);
  A->E->sun[1] = Delta;
  // added from  UZH version
  A->E->sun[8] = A->P->Lozone;
  A->E->sun[9] = A->P->alpha_iqbal;
  A->E->sun[10] = A->P->beta_iqbal;
  A->E->sun[11] = 0.;

  if (A->P->point_sim != 1)
    {
      A->E->sun[0] = A->P->latitude * GTConst::Pi / 180.;
      A->E->sun[2] = (A->P->longitude * GTConst::Pi / 180. -
                      A->P->ST * GTConst::Pi / 12. + Et) /
                     GTConst::omega;
      A->E->hsun =
        adaptiveSimpsons2(SolarHeight__, A->E->sun, JDb, JDe, 1.E-6, 20) /
        (JDe - JDb);
      A->E->dsun =
        adaptiveSimpsons2(SolarAzimuth__, A->E->sun, JDb, JDe, 1.E-6, 20) /
        (JDe - JDb) +
        A->P->dem_rotation * GTConst::Pi / 180.;
      if (A->E->dsun < 0) A->E->dsun += 2 * GTConst::Pi;
      if (A->E->dsun > 2 * GTConst::Pi) A->E->dsun -= 2 * GTConst::Pi;
      A->E->sinhsun =
        adaptiveSimpsons2(Sinalpha_, A->E->sun, JDb, JDe, 1.E-6, 20) /
        (JDe - JDb);
      if (A->P->cast_shadow == 1)
        shadow_haiden(A->T->Z0, A->E->hsun, A->E->dsun, A->L->shadow);
    }

  //  INITIALIZE BASIN AVERAGES
  if (A->P->Dtplot_basin[geotop::common::Variables::i_sim] > 0)
    {
      for (long j = 0; j < geotop::common::Variables::nobsn; j++)
        {
          geotop::common::Variables::odb[j] = 0.;
        }
    }

  // CLOUDINESS
  // was: find_actual_cloudiness(&(A->M->tau_cloud), &(A->M->tau_cloud_av),
  // &(A->M->tau_cloud_yes), &(A->M->tau_cloud_av_yes), A->M, JDb, JDe, Delta,
  // E0, Et, A->P->ST, 0.);

#ifndef USE_INTERNAL_METEODISTR

  for (unsigned long ii = 1; ii < A->M->st->Z.size(); ii++)
    {
      if (A->M->st->flag_SW_meteoST[ii] ==
          1)    // if that meteo station measures cloudiness
        {
          find_actual_cloudiness_meteoio(&(A->M->st->tau_cloud_meteoST[ii]),
                                         &(A->M->st->tau_cloud_av_meteoST[ii]),
                                         &(A->M->st->tau_cloud_yes_meteoST[ii]),
                                         &(A->M->st->tau_cloud_av_yes_meteoST[ii]),
                                         ii, A->M, vec_meteo, JDb, JDe, Delta, E0,
                                         Et, A->P->ST, 0., A->P->Lozone,
                                         A->P->alpha_iqbal, A->P->beta_iqbal, 0.);
        }
    }

  // call the function that interpolates the cloudiness
  meteoio_interpolate_cloudiness(A->P, JDe, A->M->tau_cloud,
                                 A->M->st->tau_cloud_meteoST);
  meteoio_interpolate_cloudiness(
    A->P, JDe, A->M->tau_cloud_av,
    A->M->st->tau_cloud_av_meteoST);  // Matteo: just added 17.4.2012

  // tau_cloud and tau_cloud_av are available as matrices

  /*if (A->M->st->tau_cloud_meteoST[A->M->nstcloud] !=
  geotop::input::gDoubleNoValue) tau_cloud =
  A->M->st->tau_cloud_meteoST[A->M->nstcloud];

  if (A->M->st->tau_cloud_av_meteoST[A->M->nstcloud] !=
  geotop::input::gDoubleNoValue) tau_cloud_av =
  A->M->st->tau_cloud_av_meteoST[A->M->nstcloud];

  if (A->M->st->tau_cloud_yes_meteoST[A->M->nstcloud] !=
  geotop::input::gDoubleNoValue) tau_cloud_yes =
  A->M->st->tau_cloud_yes_meteoST[A->M->nstcloud];

  if (A->M->st->tau_cloud_av_yes_meteoST[A->M->nstcloud] !=
  geotop::input::gDoubleNoValue) tau_cloud_av_yes =
  A->M->st->tau_cloud_av_yes_meteoST[A->M->nstcloud];*/

  // tau_cloud(s) are available from meteoio, therefore tau_cloud_yes is always
  // set at 1

#else

  find_actual_cloudiness(&tau_cloud, &tau_cloud_av, &tau_cloud_yes,
                         &tau_cloud_av_yes, i, A->M, JDb, JDe, Delta, E0, Et,
                         A->P->ST, 0., A->P->Lozone, A->P->alpha_iqbal,
                         A->P->beta_iqbal, 0.);
  // TODO use tau_cloud tau_cloud_av as matrices

#endif

#ifdef VERY_VERBOSE
  printf("t:%f cl:%f clav:%f cly:%d clavy:%d", A->I->time, tau_cloud,
         tau_cloud_av, tau_cloud_yes, tau_cloud_av_yes);
#endif

  //  POINT ENERGY BALANCE
  for (i = 1; i <= A->P->total_channel + A->P->total_pixel; i++)
    {
      // channel or land?
      if (i <= A->P->total_channel)    // CHANNEL
        {

          r = A->C->r[i];
          c = A->C->c[i];
        }
      else      // LAND
        {

          r = A->T->rc_cont[i - A->P->total_channel][1];
          c = A->T->rc_cont[i - A->P->total_channel][2];
        }

#ifndef USE_INTERNAL_METEODISTR
      // this temporary.. it sohuld be removed once we will work everywhere with
      // tau matrices
      tau_cloud = A->M->tau_cloud[r][c];
      tau_cloud_av = A->M->tau_cloud_av[r][c];
      if ((long)tau_cloud == geotop::input::gDoubleNoValue)
        {
          tau_cloud_yes = 0;
        }
      else
        {
          tau_cloud_yes = 1;
        }
      if ((long)tau_cloud_av == geotop::input::gDoubleNoValue)
        {
          tau_cloud_av_yes = 0;
        }
      else
        {
          tau_cloud_av_yes = 1;
        }

#endif

      if (A->L->delay[r][c] <= A->I->time / GTConst::secinday)
        {
          cnt++;
          sux = PointEnergyBalance(i, r, c, Dt, JDb, JDe, L, C, S, G, V, snowage, A,
                                   E0, Et, Dtplot, *W, f, &SWup, &Tgskin, tau_cloud,
                                   tau_cloud_av, tau_cloud_yes, tau_cloud_av_yes);

          if (sux == 1) return 1;
        }
      else
        {
          SWup = 0.;
          Tgskin = 0.;
        }

      // surroundings
      if (i > A->P->total_channel)
        {
          if (A->P->surroundings == 1)
            {
              SWrefl_surr_ave += SWup / (double)A->P->total_pixel;
              Tgskin_surr_ave += Tgskin / (double)A->P->total_pixel;
            }
          else
            {
              A->E->SWrefl_surr[r][c] = SWup;
              A->E->Tgskin_surr[r][c] = Tgskin;
            }
        }
    }

  if (A->P->surroundings == 1)
    {
      for (i = 1; i <= A->P->total_pixel; i++)
        {
          r = A->T->rc_cont[i][1];
          c = A->T->rc_cont[i][2];
          A->E->SWrefl_surr[r][c] = SWrefl_surr_ave;
          A->E->Tgskin_surr[r][c] = Tgskin_surr_ave;
        }
    }

  return 0;
}

// end of "EnergyBalance" subroutine

//******************************************************************************************************************************************/

short PointEnergyBalance(size_t i,
                         long r,
                         long c,
                         double Dt,
                         double JDb,
                         double JDe,
                         SoilState *L,
                         SoilState *C,
                         Statevar3D *S,
                         Statevar3D *G,
                         StateVeg *V,
                         GeoVector<double> &snowage,
                         AllData *A,
                         double E0,
                         double Et,
                         double /*Dtplot*/,
                         double W,
                         FILE * /*f*/,
                         double *SWupabove_v,
                         double *Tgskin,
                         double tau_cloud,
                         double tau_cloud_av,
                         short tau_cloud_yes,
                         short tau_cloud_av_yes)
{
  long l = 0, j = 0, ns = 0, ng = 0;
  double SWin, SW, SWbeam, SWdiff, SWv_vis, SWv_nir, SWg_vis, SWg_nir, cosinc,
         avis_b, avis_d, anir_b, anir_d;
  double SWupabove_v_vis, SWupabove_v_nir;
  double avis_ground = 0., anir_ground = 0.;
  short SWb_yes = 0;
  double tauatm_sinhsun;
  double LWin, LW, LWv, epsa, eps, epsa_min, epsa_max, LWupabove_v;
  double H, LE, E, Hv, LEv, Etrans, Hg0, Eg0, Hg1, Eg1, Hadv = 0.0;
  double surfEB = 0.0, GEF = 0.0, dUsl = 0.0;
  double snowD, RainOnSnow = 0., fsnow, fsnownew, fsnowcan, Melt_snow = 0.,
                Evap_snow = 0., Melt_glac = 0., Evap_glac = 0., Evap_soil = 0.;
  double k_snowred, maxSWE, k;
  double Prain_over, Psnow_over, Prain = 0.0, Psnow = 0.0, drip_rain, drip_snow,
                                 max_wcan_rain = 0., max_wcan_snow = 0.,
                                 theta_sup;
  double z0, z0veg = 0., d0, d0veg = 0., z0_z0t, hveg = 0.;
  double rh, rv, rc, rb, ruc, Lobukhov;
  double fc, fc0, decaycoeff, Locc, u_top, Qv;
  double ea, Tdew, Qa, RHpoint, Vpoint, Ppoint, Tpoint, Precpoint, zmeas_T,
         zmeas_u, Tdirichlet;
  double ilwr_point;  // Incoming Longwave Radiation in the point (r,c) as an
  // interrogation from the map ILWRgrid that is given as
  // input
  double Ts, Qs, Qg;
  long sy;
  short lu;
  short sux, surface;  //(1=yes, 0=no)
  double ic = 0., wa, rho = 0.;
  long lpb;

  //    geotop::logger::GlobalLogger* lg =
  //    geotop::logger::GlobalLogger::getInstance();

  // initialization of cumulated water volumes and set soil ancillary state vars
  if (i <= A->P->total_channel)
    {
      j = i;
      lu = (short)A->L->LC[r][c];
      sy = A->C->soil_type[j];
      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          A->C->ET[l][j] = 0.0;
          A->C->th[l][j] = theta_from_psi(C->P[l][j], C->thi[l][j], l, A->S->pa, sy,
                                          GTConst::PsiMin);
          A->C->th[l][j] =
            Fmin(A->C->th[l][j], A->S->pa[sy][jsat][l] - C->thi[l][j]);
        }

    }
  else
    {
      j = i - A->P->total_channel;
      A->W->Pnet[r][c] = 0.0;
      A->W->HN[r][c] = 0.0;  // TODO mattiu
      lu = (short)A->L->LC[r][c];
      sy = A->S->type[r][c];

      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          A->S->ET[l][r][c] = 0.0;
          A->S->th[l][j] = theta_from_psi(L->P[l][j], L->thi[l][j], l, A->S->pa, sy,
                                          GTConst::PsiMin);
          A->S->th[l][j] =
            Fmin(A->S->th[l][j], A->S->pa[sy][jsat][l] - L->thi[l][j]);
        }
    }

  // METEO

  Tpoint = A->M->Tgrid[r][c];
  Ppoint = A->M->Pgrid[r][c];
  RHpoint = A->M->RHgrid[r][c];
  Vpoint = A->M->Vgrid[r][c];
  ilwr_point = A->M->ILWRgrid[r][c];  // TODO: to be added once WRF has ilwr

  Precpoint = A->W->PrecTot[r][c];
  // define prec as normal (not vertical)

  Precpoint *= cos(A->T->slope[r][c] * GTConst::Pi / 180.);

  // another cosine correction applied for point simulations (due to area
  // proejection)

  if (A->P->point_sim == 1 && A->P->flag1D == 0)
    Precpoint *= cos(A->T->slope[r][c] * GTConst::Pi / 180.);

  // TODO: integrate Tdirichlet
  // from UZH version:
  // Tdirichlet=A->M->var[A->M->nstTs-1][iTs];  <-segfault on this: to better
  // understand..
  Tdirichlet = geotop::input::gDoubleNoValue;  // this as temporary hack to make
  // things work

  if ((long)Tdirichlet == geotop::input::gDoubleNoValue ||
      (long)Tdirichlet == geotop::input::gDoubleAbsent)
    Tdirichlet = A->P->Tsup;

  // SNOW

  snowD = DEPTH(r, c, S->lnum, S->Dzl);
  ns = S->lnum[r][c];

  // vegetation parameters

  if (i > A->P->total_channel)
    {
      for (l = 1; l <= jdvegprop; l++)
        {
          if ((long)A->L->vegparv[lu - 1][l] != geotop::input::gDoubleNoValue)
            {
              A->L->vegpar[l] = A->L->vegparv[lu - 1][l];
              if (l == jdroot)
                root(A->L->root_fraction.getCols(), A->L->vegpar[jdroot], 0.0,
                     A->S->pa, A->L->root_fraction, lu);
            }
          else
            {
              A->L->vegpar[l] = A->L->ty[lu][l + jHveg - 1];
            }
        }

      if (snowD > A->L->vegpar[jdz0thresveg])
        {
          fsnow = 1.0;
        }
      else if (snowD > A->L->vegpar[jdz0thresveg2])
        {
          fsnow = (snowD - A->L->vegpar[jdz0thresveg2]) /
                  (A->L->vegpar[jdz0thresveg] - A->L->vegpar[jdz0thresveg2]);
        }
      else
        {
          fsnow = 0.0;
        }

      fc = A->L->vegpar[jdcf] * pow(1.0 - fsnow, A->L->vegpar[jdexpveg]);

      // control if LSAI is too low (in order to prevent numerical problems)

      if (A->L->vegpar[jdLSAI] < GTConst::LSAIthres) fc = 0.0;

      // glacier
      if (A->P->max_glac_layers > 0)
        {
          ng = G->lnum[r][c];
          if (ng > 0)
            {
              A->L->vegpar[jdLSAI] = 0.0;
              fc = 0.0;
            }
        }

    }
  else
    {
      A->L->vegpar.resize(A->L->vegpar.size() + 1, 0.0);
      fc = 0.0;
    }

  //  METEOROLOGICAL COMPUTATIONS

  // temperature and wind velocity measurement heights

  zmeas_u = A->M->st->Vheight[1];
  zmeas_T = A->M->st->Theight[1];

  // RAIN AND SNOW PRECIPITATION [mm]
  // convert total precipitation to [mm]

  Precpoint *=
    (Dt / 3600.0);  // from [mm/h] to [mm] TODO: HACK To verify with MeteoIO

  ea = RHpoint * SatVapPressure(Tpoint, Ppoint);  // Vapour Pressure [mbar]
  Qa = SpecHumidity(ea, Ppoint);                  // Specific Humidity
  Tdew = TfromSatVapPressure(ea, Ppoint);         // Dew Temperature

  // distinguish between rain and snow

  if (A->P->dew == 1)
    {
      // on the base of the dew temperature of the air
      part_snow(Precpoint, &Prain_over, &Psnow_over, Tdew, A->P->T_rain,
                A->P->T_snow);
    }
  else
    {
      // on the base of the temperature of the air
      part_snow(Precpoint, &Prain_over, &Psnow_over, Tpoint, A->P->T_rain,
                A->P->T_snow);
    }

#ifndef USE_INTERNAL_METEODISTR
  // added snow and rain correction factor (meteodistr uses them internally)
  Psnow_over *= A->P->snowcorrfact;
  Prain_over *= A->P->raincorrfact;
#endif

  // Adjusting snow precipitation in case of steep slope (contribution by
  // Stephan Gruber)
  if (A->P->snow_curv > 0 && A->T->slope[r][c] > A->P->snow_smin)
    {
      if (A->T->slope[r][c] <= A->P->snow_smax)
        {
          k_snowred =
            (exp(-pow(A->T->slope[r][c] - A->P->snow_smin, 2.) / A->P->snow_curv) -
             exp(-pow(A->P->snow_smax, 2.) / A->P->snow_curv));
        }
      else
        {
          k_snowred = 0.0;
        }
      Psnow_over *= k_snowred;
    }

  // precipitation at the surface in the unvegetated fraction
  Prain = (1. - fc) * Prain_over;
  Psnow = (1. - fc) * Psnow_over;

  if (fc > 0)
    {
      canopy_rain_interception(ratio_max_storage_RAIN_over_canopy_to_LSAI,
                               A->L->vegpar[jdLSAI], Prain_over, &max_wcan_rain,
                               &(V->wrain[j]), &drip_rain);
      canopy_snow_interception(ratio_max_storage_SNOW_over_canopy_to_LSAI,
                               A->L->vegpar[jdLSAI], Psnow_over, V->Tv[j], Vpoint,
                               Dt, &max_wcan_snow, &(V->wsnow[j]), &drip_snow);

      // precipitation at the surface in the vegetated fraction
      Prain += fc * drip_rain;
      Psnow += fc * drip_snow;

#ifdef USE_DOUBLE_PRECISION_OUTPUT
      if (V->wsnow[j] < 0)
        printf("Error 1 wcansnow:%12g %ld %ld\n", V->wsnow[j], r, c);
#else
      if (V->wsnow[j] < 0)
        printf("Error 1 wcansnow:%f %ld %ld\n", V->wsnow[j], r, c);
#endif

      if (V->wsnow[j] < 0) V->wsnow[j] = 0.0;
    }

  // SHORTWAVE RADIATION
  // initialization of shortwave radiation absorbed by soil
  SW = 0.0;
  *SWupabove_v = 0.0;

#ifdef VERY_VERBOSE
  printf("tcav:%f tc:%f Tp:%f Rh:%f Z:%f lsTa:%f lsTdew:%f \n", tau_cloud_av,
         tau_cloud, Tpoint, RHpoint, A->T->Z0[r][c], A->M->LRv[ilsTa],
         A->M->LRv[ilsTdew]);
#endif

  // CURRENT CLOUDINESS
  //  if averaged cloud transmissivity is not available it is estimated
  //through this micromet subroutine
  if (tau_cloud_av_yes == 0)
    tau_cloud_av =
      1. - 0.71 * find_cloudfactor(Tpoint, RHpoint, A->T->Z0[r][c],
                                   A->M->LRv[ilsTa],
                                   A->M->LRv[ilsTdew]);  // Kimball(1928)

  // in case of shortwave data not available
  if (tau_cloud_yes == 0) tau_cloud = tau_cloud_av;

#ifdef VERY_VERBOSE
  printf("tcav:%f tc:%f Tp:%f Rh:%f Z:%f lsTa:%f lsTdew:%f \n", tau_cloud_av,
         tau_cloud, Tpoint, RHpoint, A->T->Z0[r][c], A->M->LRv[ilsTa],
         A->M->LRv[ilsTdew]);
#endif

  //    std::string filename;
  //        filename = "TauCloud.txt";
  //        FILE * f1;
  //        long year,day,month,hour,minute;
  //        double JDfrom0,JD;
  //        JDfrom0 = convert_tfromstart_JDfrom0(A->I->time+A->P->Dt,
  //        A->P->init_date[geotop::common::Variables::i_sim]);
  //        convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
  //        convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour,
  //        &minute);
  //
  //        //if(r==233 & c==192){
  //        f1=fopen(filename.c_str(),"a");
  //        fprintf(f1,"%02.0f/%02.0f/%04.0f
  //    %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);
  //        fprintf(f1,",%f,%f\n",A->M->tau_cl_av_map[r][c],A->M->tau_cl_map[r][c]);
  //        fclose(f1);

  // albedo
  // update ground albedo
  if (i <= A->P->total_channel)
    {
      theta_sup = A->C->th[1][j];
    }
  else
    {
      theta_sup = A->S->th[1][j];
    }
  avis_ground =
    find_albedo(A->L->ty[lu][ja_vis_dry], A->L->ty[lu][ja_vis_sat], theta_sup,
                A->S->pa[sy][jres][1], A->S->pa[sy][jsat][1]);
  anir_ground =
    find_albedo(A->L->ty[lu][ja_nir_dry], A->L->ty[lu][ja_nir_sat], theta_sup,
                A->S->pa[sy][jres][1], A->S->pa[sy][jsat][1]);

  if (snowD > 0)
    {
      /*
       * This piece of code may be interesting in debugging snow_age and
       snow_albedo calculation
       * To add the header, open the connection in input.cc
       *  std::string filename;
          filename = geotop::common::Variables::WORKING_DIRECTORY+
       "log_albedo.txt"; FILE *f1; long year,day,month,hour,minute; double
       JDfrom0,JD; JDfrom0 = convert_tfromstart_JDfrom0(A->I->time+A->P->Dt,
       A->P->init_date); convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
          convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour,
       &minute); f1=fopen(filename.c_str(),"a"); fprintf(f1,"%02.0f/%02.0f/%04.0f
       %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);
          fprintf(f1,",%f, %f, %f, %f, %f, %f, %f\n",avis_ground, snowD,
       A->P->aep, A->P->avo, A->P->snow_aging_vis, snowage[j], cosinc);
          fclose(f1);
      */

      // update snow albedo
      // if(i>A->P->total_channel) update_snow_age(Psnow_over, S->T[ns][r][c], Dt,
      // A->P->minP_torestore_A, &(snowage[j]));
      if (i > A->P->total_channel)
        update_snow_age_cumEvent(Psnow_over, Dt, &(A->P->cum_prec),
                                 &(A->P->cum_da_up), &(A->P->time_wo_prec),
                                 &(A->P->evento), &(A->P->up_albedo),
                                 A->P->tres_wo_prec, A->P->minP_torestore_A,
                                 S->T[ns][r][c], &(snowage[j]));

      // from cosinc to 0 as matteo suggested 06.11.2013 SC
      avis_d =
        snow_albedo(avis_ground, snowD, A->P->aep, A->P->avo,
                    A->P->snow_aging_vis, snowage[j], 0., (&Turbulence::Zero));
      avis_b = avis_d;  // approx
      anir_d =
        snow_albedo(anir_ground, snowD, A->P->aep, A->P->airo,
                    A->P->snow_aging_nir, snowage[j], 0., (&Turbulence::Zero));
      anir_b = anir_d;  // approx
    }
  else
    {
      avis_b = avis_ground;
      avis_d = avis_ground;
      anir_b = anir_ground;
      anir_d = anir_ground;
    }

  // calculation of SWin
  if (A->P->point_sim == 1)
    {
      A->E->sun[0] = A->T->latitude[r][c] * GTConst::Pi / 180.;
      A->E->sun[2] = (A->T->longitude[r][c] * GTConst::Pi / 180. -
                      A->P->ST * GTConst::Pi / 12. + Et) /
                     GTConst::omega;

      A->E->hsun =
        adaptiveSimpsons2(SolarHeight__, A->E->sun, JDb, JDe, 1.E-6, 20) /
        (JDe - JDb);
      A->E->dsun =
        adaptiveSimpsons2(SolarAzimuth__, A->E->sun, JDb, JDe, 1.E-6, 20) /
        (JDe - JDb) +
        A->P->dem_rotation * GTConst::Pi / 180.;

      if (A->E->dsun < 0) A->E->dsun += 2 * GTConst::Pi;
      if (A->E->dsun > 2 * GTConst::Pi) A->E->dsun -= 2 * GTConst::Pi;
      A->E->sinhsun =
        adaptiveSimpsons2(Sinalpha_, A->E->sun, JDb, JDe, 1.E-6, 20) /
        (JDe - JDb);

      if (A->P->cast_shadow == 1)
        A->L->shadow[r][c] =
          shadows_point(A->T->horizon_height[A->T->horizon_point[r][c] - 1],
                        A->T->horizon_numlines[A->T->horizon_point[r][c] - 1],
                        A->E->hsun * 180. / GTConst::Pi,
                        A->E->dsun * 180 / GTConst::Pi, 0., 0.);
    }

  A->E->sun[3] = RHpoint;
  A->E->sun[4] = Tpoint;
  A->E->sun[5] = Ppoint;
  A->E->sun[6] = A->T->slope[r][c] * GTConst::Pi / 180.;
  A->E->sun[7] = A->T->aspect[r][c] * GTConst::Pi / 180.;

  // TODO: merge it #if(A->P->albedoSWin != 0) A->E->sun[11] = (avis_b + avis_d
  // + anir_b + anir_d)/4.;
  shortwave_radiation(JDb, JDe, A->E->sun, A->E->sinhsun, E0, A->T->sky[r][c],
                      A->E->SWrefl_surr[r][c], tau_cloud, A->L->shadow[r][c],
                      &SWbeam, &SWdiff, &cosinc, &tauatm_sinhsun, &SWb_yes);

  // TODO: merge it #SWbeam=flux(A->M->nstsrad, iSWb, A->M->var, 1.0, 0.0,
  // SWbeam);
  // TODO: merge it #SWdiff=flux(A->M->nstsrad, iSWd, A->M->var, 1.0,
  // (1.-A->T->sky->co[r][c])*A->E->SWrefl_surr->co[r][c], SWdiff);

  SWin = SWbeam + SWdiff;

  // update snow albedo
  if (snowD > 0)
    {
      avis_b = snow_albedo(avis_ground, snowD, A->P->aep, A->P->avo,
                           A->P->snow_aging_vis, snowage[j], cosinc, (*Fzen));
      anir_b = snow_albedo(anir_ground, snowD, A->P->aep, A->P->airo,
                           A->P->snow_aging_nir, snowage[j], cosinc, (*Fzen));
    }

  // shortwave absorbed by soil (when vegetation is not present)
  SW += (1.0 - fc) * (SWdiff * (1.0 - 0.5 * avis_d - 0.5 * anir_d) +
                      SWbeam * (1.0 - 0.5 * avis_b - 0.5 * anir_b));
  *SWupabove_v += (1.0 - fc) * (SWdiff * (0.5 * avis_d + 0.5 * anir_d) +
                                SWbeam * (0.5 * avis_b + 0.5 * anir_b));

  // shortwave radiation absorbed by canopy
  if (fc > 0)
    {
      if (A->E->hsun > 0)
        {
          if (V->wsnow[j] < 0)
            printf("Error wcansnow:%f %ld %ld\n", V->wsnow[j], r, c);

          if (V->wsnow[j] < 0) V->wsnow[j] = 0.0;
          fsnowcan = pow(V->wsnow[j] / max_wcan_snow, 2. / 3.);
          shortwave_vegetation(
            0.5 * SWdiff, 0.5 * SWbeam, cosinc, fsnowcan, GTConst::wsn_vis,
            GTConst::Bsnd_vis, GTConst::Bsnb_vis, avis_d, avis_b,
            A->L->ty[lu][jvCh], A->L->ty[lu][jvR_vis], A->L->ty[lu][jvT_vis],
            A->L->vegpar[jdLSAI], &SWv_vis, &SWg_vis, &SWupabove_v_vis);
          shortwave_vegetation(
            0.5 * SWdiff, 0.5 * SWbeam, cosinc, fsnowcan, GTConst::wsn_nir,
            GTConst::Bsnd_nir, GTConst::Bsnb_nir, anir_d, anir_b,
            A->L->ty[lu][jvCh], A->L->ty[lu][jvR_nir], A->L->ty[lu][jvT_nir],
            A->L->vegpar[jdLSAI], &SWv_nir, &SWg_nir, &SWupabove_v_nir);

        }
      else
        {
          SWv_vis = 0.0;
          SWg_vis = 0.0;
          SWv_nir = 0.0;
          SWg_nir = 0.0;
          SWupabove_v_vis = 0.0;
          SWupabove_v_nir = 0.0;
        }

      SW += fc * (SWg_vis + SWg_nir);
      *SWupabove_v += fc * (SWupabove_v_vis + SWupabove_v_nir);

    }
  else
    {
      SWv_vis = 0.0;
      SWv_nir = 0.0;
      SWupabove_v_vis = 0.0;
      SWupabove_v_nir = 0.0;
    }

  // correct in case of reading data
  // printf("iSWn: %d    iLWi: %d\n", iSWn, iLWi);
  // HACK: EGGER: the next line is superfluous for now, we're not reading iSWn:
  // TODO: merge it # SW=flux(A->M->nstsrad, iSWn, A->M->var, 1.0, 0.0, SW);
  // SW=flux(1, iSWn, A->M->var, 1.0, SW);
  // Extinction coefficient for SW in the snow layers
  rad_snow_absorption(r, c, A->E->SWlayer, SW, S);

  // LONGWAVE RADIATION
  // soil-snow emissivity
  if (snowD > 10)
    {
      eps = A->P->epsilon_snow;
    }
  else
    {
      eps = A->L->ty[lu][jemg];
    }

  if (!A->P->use_ilwr_wrf)
    {
      longwave_radiation(A->P->state_lwrad, ea, RHpoint, Tpoint, A->P->k1,
                         A->P->k2, tau_cloud_av, &epsa, &epsa_max, &epsa_min);
      LWin = A->T->sky[r][c] * epsa * SB(Tpoint);
      // if LWin is in the meteo file, usi this value instead
      // double flux(long i, long icol, double **met, double k, double est){
#ifdef USE_INTERNAL_METEODISTR
      LWin = flux(A->M->nstlrad, iLWi, A->M->var, 1.0, LWin);
#endif
    }
  else
    {
      LWin = ilwr_point;
    }

  if (A->P->surroundings == 1)
    {
      /* LWin corrected with temperature of surrounding terrain,
       * calculated as averaged. */
      LWin += (1. - A->T->sky[r][c]) * eps * SB(A->E->Tgskin_surr[r][c]);
    }

  // if incoming LW data are available, they are used (priority)
  // HACK: EGGER: the next line is superfluous for now, we're not reading iLWi:
  // TODO: merge it #LWin=flux(A->M->nstlrad, iLWi, A->M->var, 1.0, LWin);

  // roughness lengths
  update_roughness_soil(A->L->ty[lu][jz0], 0.0, 0.0, snowD,
                        A->L->ty[lu][jz0thressoil], A->P->z0_snow, &z0, &d0,
                        &z0_z0t);
  if (fc > 0)
    update_roughness_veg(A->L->vegpar[jdHveg], snowD, zmeas_u, zmeas_T, &z0veg,
                         &d0veg, &hveg);

  // variables used in the point_surface_balance
  for (l = 1; l <= geotop::common::Variables::Nl + ns + ng; l++)
    {
      if (l <= ns)    // snow
        {
          A->E->Dlayer[l] = 1.E-3 * S->Dzl[ns + 1 - l][r][c];
          A->E->liq[l] = S->w_liq[ns + 1 - l][r][c];
          A->E->ice[l] = S->w_ice[ns + 1 - l][r][c];
          A->E->Temp[l] = S->T[ns + 1 - l][r][c];
        }
      else if (l <= ns + ng)      // glacier
        {
          A->E->Dlayer[l] = 1.E-3 * G->Dzl[ns + ng + 1 - l][r][c];
          A->E->liq[l] = G->w_liq[ns + ng + 1 - l][r][c];
          A->E->ice[l] = G->w_ice[ns + ng + 1 - l][r][c];
          A->E->Temp[l] = G->T[ns + ng + 1 - l][r][c];
        }
      else      // soil
        {
          if (i <= A->P->total_channel)
            {
              A->E->Dlayer[l] = 1.E-3 * A->S->pa[sy][jdz][l - ns - ng];
              A->E->liq[l] =
                A->C->th[l - ns - ng][j] * A->E->Dlayer[l] * GTConst::rho_w;
              A->E->ice[l] =
                C->thi[l - ns - ng][j] * A->E->Dlayer[l] * GTConst::rho_w;
              A->E->Temp[l] = C->T[l - ns - ng][j];
            }
          else
            {
              A->E->Dlayer[l] = 1.E-3 * A->S->pa[sy][jdz][l - ns - ng];
              A->E->liq[l] =
                A->S->th[l - ns - ng][j] * A->E->Dlayer[l] * GTConst::rho_w;
              A->E->ice[l] =
                L->thi[l - ns - ng][j] * A->E->Dlayer[l] * GTConst::rho_w;
              A->E->Temp[l] = L->T[l - ns - ng][j];
            }
        }
    }

  // SOIL AND SNOW PROPERTIES
  sux = 0;
  lpb = 0;  // layerproblematic

  do
    {
      // init
      surface = 0;
      ic = 0.;
      wa = 0.;

      // snow or glacier melting with T[0] >= 0
      if (sux == -1)
        {
          surface = 1;

          // it is not admitted that a snow/glacier layer melts completely
          // snow layer not alone
        }
      else if (sux == -2)
        {
          merge(A->P->alpha_snow, A->E->ice, A->E->liq, A->E->Temp, A->E->Dlayer,
                lpb, 1, ns, ns + ng + geotop::common::Variables::Nl);
          ns--;

          // glacier layer not alone
        }
      else if (sux == -3)
        {
          merge(A->P->alpha_snow, A->E->ice, A->E->liq, A->E->Temp, A->E->Dlayer,
                lpb, ns + 1, ns + ng, ns + ng + geotop::common::Variables::Nl);
          ng--;

          // just one glacier layer and snow layers >= 0
        }
      else if (sux == -4)
        {
          merge(A->P->alpha_snow, A->E->ice, A->E->liq, A->E->Temp, A->E->Dlayer,
                ns + ng, ns, ns + ng, ns + ng + geotop::common::Variables::Nl);
          ng--;

          // just one snow layer and glacier layers > 0
        }
      else if (sux == -5)
        {
          merge(A->P->alpha_snow, A->E->ice, A->E->liq, A->E->Temp, A->E->Dlayer,
                ns, ns, ns + 1, ns + ng + geotop::common::Variables::Nl);
          ns--;

          // just one snow layer and glacier layers == 0
        }
      else if (sux == -6)
        {
          ic = A->E->ice[1];
          wa = A->E->liq[1];
          rho = ic / A->E->Dlayer[1];
          for (l = 2; l <= 1 + geotop::common::Variables::Nl; l++)
            {
              A->E->ice[l - 1] = A->E->ice[l];
              A->E->liq[l - 1] = A->E->liq[l];
              A->E->Temp[l - 1] = A->E->Temp[l];
              A->E->Dlayer[l - 1] = A->E->Dlayer[l];
            }
          ns--;
          A->E->Dlayer[1] *= (1. + ic / (A->E->ice[1] + A->E->liq[1]));

          A->E->ice[1] += ic;
          if (i <= A->P->total_channel)
            {
              C->thi[1][j] = A->E->ice[1] / (A->E->Dlayer[1] * GTConst::rho_w);
              A->C->th[1][j] = A->E->liq[1] / (A->E->Dlayer[1] * GTConst::rho_w);
            }
          else
            {
              L->thi[1][j] = A->E->ice[1] / (A->E->Dlayer[1] * GTConst::rho_w);
              A->S->th[1][j] = A->E->liq[1] / (A->E->Dlayer[1] * GTConst::rho_w);
            }
        }

      // skin temperature
      A->E->Temp[0] = A->E->Temp[1];

      // ENERGY BALANCE
      sux = SolvePointEnergyBalance(
              surface, Tdirichlet, A->P->EB, A->P->Cair, A->P->micro,
              JDb - A->P->init_date, Dt, i, j, r, c, L, C, V, A->E, A->L, A->S, A->C,
              A->T, A->P, ns, ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg, d0veg, 1.0,
              hveg, Vpoint, Tpoint, Qa, Ppoint, A->M->LRv[ilsTa], eps, fc,
              A->L->vegpar[jdLSAI], A->L->vegpar[jddecay0], &(V->wrain[j]),
              max_wcan_rain, &(V->wsnow[j]), max_wcan_snow, SWin, LWin,
              SWv_vis + SWv_nir, &LW, &H, &E, &LWv, &Hv, &LEv, &Etrans, &Ts, &Qs, Hadv,
              &Hg0, &Hg1, &Eg0, &Eg1, &Qv, &Qg, &Lobukhov, &rh, &rv, &rb, &rc, &ruc,
              &u_top, &decaycoeff, &Locc, &LWupabove_v, &lpb, &dUsl);

    }
  while (sux < 0);

  if (sux == 0)
    {
      // skin temperature
      *Tgskin = A->E->Temp[A->P->nsurface];

      // snow melting issue

      if (ic > 0)
        {
          sux_minus6_condition(ic, wa, rho, 1.E-3 * A->S->pa[sy][jdz][1], A->E);
          ns++;
        }

      if (i <= A->P->total_channel)
        {
          update_soil_channel(A->P->nsurface, ns + ng, i, fc, Dt, A->E, A->S->pa,
                              sy, C, A->C->ET, A->C->th);

        }
      else
        {
          // surface energy balance
          LE =
            E * Turbulence::latent(A->E->Temp[1], Turbulence::Levap(A->E->Temp[1]));
          surfEB = SW + LW - H - LE;

          if (ns == 0)
            {
              GEF = surfEB;
            }
          else
            {
              k = k_thermal(
                    0, 1,
                    Arithmetic_Mean(
                      A->E->Dlayer[ns], A->E->Dlayer[ns + 1],
                      A->E->liq[ns] / (GTConst::rho_w * A->E->Dlayer[ns]),
                      A->E->liq[ns + 1] / (GTConst::rho_w * A->E->Dlayer[ns + 1])),
                    Arithmetic_Mean(
                      A->E->Dlayer[ns], A->E->Dlayer[ns + 1],
                      A->E->ice[ns] / (GTConst::rho_w * A->E->Dlayer[ns]),
                      A->E->ice[ns + 1] / (GTConst::rho_w * A->E->Dlayer[ns + 1])),
                    Arithmetic_Mean(A->E->Dlayer[ns], A->E->Dlayer[ns + 1], 1.,
                                    A->S->pa[sy][jsat][1]),
                    Arithmetic_Mean(A->E->Dlayer[ns], A->E->Dlayer[ns + 1], 0.,
                                    A->S->pa[sy][jkt][1]));
              GEF = k * (A->E->Temp[ns] - A->E->Temp[ns + 1]) /
                    (0.5 * A->E->Dlayer[ns] + 0.5 * A->E->Dlayer[ns + 1]);
            }

          // assign evaporation
          if (ns > 0)
            {
              Evap_snow += E * Dt;
            }
          else if (ng > 0)
            {
              Evap_glac += E * Dt;
            }
          else
            {
              Evap_soil += E * Dt;
            }

          update_soil_land(A->P->nsurface, ns + ng, j, r, c, fc, Dt, A->E, A->S->pa,
                           sy, L, A->S->ET, A->S->th);

          // glacier
          if (A->P->max_glac_layers > 0)
            {
              // glacier melting
              WBglacier(ns, ng, r, c, G, &Melt_glac, A->P, A->E, Evap_glac);
              // adjust glacier layers
              snow_layer_combination(A->P->alpha_snow, r, c, G, Tpoint,
                                     A->P->inf_glac_layers, A->P->max_weq_glac,
                                     1.E10);
            }

          // SNOW
          // snow melting, compactation and percolation

          WBsnow(Dt, ns, r, c, S, &Melt_snow, &RainOnSnow, A->P, A->T->slope[r][c],
                 Prain, A->E, Evap_snow);

          // adjust snow layers
          if (A->P->point_sim == 1)
            {
              //  maxSWE = A->P->maxSWE->co[r][c];
              maxSWE = A->P->maxSWE[r][c];
            }
          else
            {
              maxSWE = 1.E10;
            }
          snow_layer_combination(A->P->alpha_snow, r, c, S, Tpoint,
                                 A->P->inf_snow_layers, A->P->max_weq_snow, maxSWE);

          //  add new snow; TODO: choose if snow density is to be given according
          //to Valt or Jordan
          if (Psnow > 0)
            new_snow(A->P->alpha_snow, r, c, S, Psnow,
                     Psnow * GTConst::rho_w / rho_newlyfallensnow(Vpoint, Tpoint),
                     Tpoint);  // Jordan
          // if(Psnow>0) new_snow(A->P->alpha_snow, r, c, S, Psnow,
          // Psnow*GTConst::rho_w/rho_valt(Tpoint), Tpoint); // Valt

          //  NET PRECIPITATION
          A->W->Pnet[r][c] += (Melt_snow + Melt_glac + Prain);
          // A->W->HN[r][c] += Psnow*GTConst::rho_w/rho_newlyfallensnow(Vpoint,
          // Tpoint); // Jordan
          A->W->HN[r][c] += Psnow * GTConst::rho_w / rho_valt(Tpoint);  // Valt
          // VEGETATION

          if (A->L->vegpar[jdLSAI] >= GTConst::LSAIthres && ng == 0)
            {
              snowD = DEPTH(r, c, S->lnum, S->Dzl);

              fc0 = fc;

              if (snowD > A->L->vegpar[jdz0thresveg])
                {
                  fsnownew = 1.0;

                }
              else if (snowD > A->L->vegpar[jdz0thresveg2])
                {
                  fsnownew = (snowD - A->L->vegpar[jdz0thresveg2]) /
                             (A->L->vegpar[jdz0thresveg] - A->L->vegpar[jdz0thresveg2]);
                }
              else
                {
                  fsnownew = 0.0;
                }

              fc = A->L->vegpar[jdcf] * pow(1.0 - fsnownew, A->L->vegpar[jdexpveg]);

              if (fc > fc0)
                {
                  V->wrain[j] *= (fc0 / fc);

                  V->wsnow[j] *= (fc0 / fc);

                  // b) fc decreases
                }
              else if (fc < fc0)
                {
                  new_snow(A->P->alpha_snow, r, c, S, V->wsnow[j] * (fc0 - fc),
                           V->wsnow[j] * (fc0 - fc) * 1000. / 300, V->Tv[j]);

                  A->W->Pnet[r][c] += V->wrain[j] * (fc0 - fc);
                }

              if (fc < 1.E-6)
                {
                  V->wrain[j] = 0.0;

                  V->wsnow[j] = 0.0;
                }
            }

          // write output

          if (A->P->output_snow[geotop::common::Variables::i_sim] > 0)
            {
              // if(geotop::common::Variables::files[fsnowmelt] !=
              // geotop::input::gStringNoValue)
              A->N->melted[j] = Melt_snow;

              // if(geotop::common::Variables::files[fsnowsubl] !=
              // geotop::input::gStringNoValue)
              A->N->subl[j] = Evap_snow;
              // if(geotop::common::Variables::files[fsndur] !=
              // geotop::input::gStringNoValue){
              if (snowD > 0)
                {
                  A->N->yes[j] = 1;
                }
              else
                {
                  A->N->yes[j] = 0;
                }
              // }
            }

          if (A->P->max_glac_layers > 0 &&
              A->P->output_glac[geotop::common::Variables::i_sim] > 0)
            {
              if (geotop::common::Variables::files[fglacmelt] !=
                  geotop::input::gStringNoValue)
                A->G->melted[j] = Melt_glac;

              if (geotop::common::Variables::files[fglacsubl] !=
                  geotop::input::gStringNoValue)
                A->G->subl[j] = Evap_glac;
            }

          if (A->P->output_surfenergy[geotop::common::Variables::i_sim] > 0)
            {
              if (geotop::common::Variables::files[fradnet] !=
                  geotop::input::gStringNoValue)
                A->E->Rn[j] = (SW + LW);

              if (geotop::common::Variables::files[fradLWin] !=
                  geotop::input::gStringNoValue)
                A->E->LWin[j] = LWin;

              if (geotop::common::Variables::files[fradLW] !=
                  geotop::input::gStringNoValue)
                A->E->LW[j] = LW;

              if (geotop::common::Variables::files[fradSW] !=
                  geotop::input::gStringNoValue)
                A->E->SW[j] = SW;

              if (geotop::common::Variables::files[fradSWin] !=
                  geotop::input::gStringNoValue)
                A->E->SWin[j] = SWin;

              if (geotop::common::Variables::files[fradSWinbeam] !=
                  geotop::input::gStringNoValue)
                A->E->SWinb[j] = SWbeam;

              if (geotop::common::Variables::files[fshadow] !=
                  geotop::input::gStringNoValue)
                A->E->shad[j] = SWb_yes;

              if (geotop::common::Variables::files[fG] !=
                  geotop::input::gStringNoValue)
                A->E->G[j] = surfEB;

              if (geotop::common::Variables::files[fH] !=
                  geotop::input::gStringNoValue)
                A->E->H[j] = H;

              if (geotop::common::Variables::files[fLE] !=
                  geotop::input::gStringNoValue)
                A->E->LE[j] = LE;

              if (geotop::common::Variables::files[fTs] !=
                  geotop::input::gStringNoValue)
                A->E->Ts[j] = (*Tgskin);
            }

          if (A->P->output_meteo[geotop::common::Variables::i_sim] > 0)
            {
              if (geotop::common::Variables::files[fprec] !=
                  geotop::input::gStringNoValue)
                {
                  A->W->Pt[j] = Precpoint;

                  A->W->Ps[j] = Psnow_over;
                }
            }

          if (A->I->JD_plots.size() > 1 && W > 0)
            {
              if (geotop::common::Variables::files[pH] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pHg] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pG] !=
                  geotop::input::gStringNoValue)
                A->E->Hgp[j] = W * H;

              if (geotop::common::Variables::files[pH] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pHv] !=
                  geotop::input::gStringNoValue)
                A->E->Hvp[j] = W * fc * Hv;

              if (geotop::common::Variables::files[pLE] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pLEg] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pG] !=
                  geotop::input::gStringNoValue)
                A->E->LEgp[j] = W * LE;

              if (geotop::common::Variables::files[pLE] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pLEv] !=
                  geotop::input::gStringNoValue)
                A->E->LEvp[j] = W * fc * LEv;

              if (geotop::common::Variables::files[pSWin] !=
                  geotop::input::gStringNoValue)
                A->E->SWinp[j] = W * SWin;

              if (geotop::common::Variables::files[pSWg] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pG] !=
                  geotop::input::gStringNoValue)
                A->E->SWgp[j] = W * SW;

              if (geotop::common::Variables::files[pSWv] !=
                  geotop::input::gStringNoValue)
                A->E->SWvp[j] = W * fc * (SWv_vis + SWv_nir);

              if (geotop::common::Variables::files[pLWin] !=
                  geotop::input::gStringNoValue)
                A->E->LWinp[j] = W * LWin;

              if (geotop::common::Variables::files[pLWg] !=
                  geotop::input::gStringNoValue ||
                  geotop::common::Variables::files[pG] !=
                  geotop::input::gStringNoValue)
                A->E->LWgp[j] = W * LW;

              if (geotop::common::Variables::files[pLWv] !=
                  geotop::input::gStringNoValue)
                A->E->LWvp[j] = W * fc * LWv;

              if (geotop::common::Variables::files[pTs] !=
                  geotop::input::gStringNoValue)
                A->E->Tsp[j] = W * Ts;

              if (geotop::common::Variables::files[pTg] !=
                  geotop::input::gStringNoValue)
                A->E->Tgp[j] = W * (*Tgskin);
            }

          if (A->P->state_pixel == 1)
            {
              if (A->P->jplot[j] > 0 &&
                  A->P->Dtplot_point[geotop::common::Variables::i_sim] > 0)
                {
                  for (l = 0; l < geotop::common::Variables::nopnt; l++)
                    {
                      if (geotop::common::Variables::opnt[l] > 0)
                        {
                          if (geotop::common::Variables::opnt[l] == osnowover)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Psnow_over;

                            }
                          else if (geotop::common::Variables::opnt[l] == orainover)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Prain_over;

                            }
                          else if (geotop::common::Variables::opnt[l] == oprecsnow)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Psnow;

                            }
                          else if (geotop::common::Variables::opnt[l] == oprecrain)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Prain;

                            }
                          else if (geotop::common::Variables::opnt[l] == orainonsnow)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                RainOnSnow;

                            }
                          else if (geotop::common::Variables::opnt[l] == oV)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Vpoint * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oVdir)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                A->M->Vdir[r][c] * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oRH)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                RHpoint * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oPa)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Ppoint * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oTa)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Tpoint * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oTdew)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Tdew * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oTg)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                *Tgskin * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oTv)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                V->Tv[j] * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oTs)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Ts * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oEB)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                surfEB * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];
                            }
                          else if (geotop::common::Variables::opnt[l] == oG)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                GEF * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oSWin)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                SWin * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oSWb)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                SWbeam * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oSWd)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                SWdiff * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLWin)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                LWin * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == ominLWin)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                (A->T->sky[r][c] * epsa_min * SB(Tpoint) +
                                 (1. - A->T->sky[r][c]) * eps *
                                 SB(A->E->Tgskin_surr[r][c])) *
                                Dt / A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == omaxLWin)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                (A->T->sky[r][c] * epsa_max * SB(Tpoint) +
                                 (1. - A->T->sky[r][c]) * eps *
                                 SB(A->E->Tgskin_surr[r][c])) *
                                Dt / A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oSW)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                SW * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLW)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                LW * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oH)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                H * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLE)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                LE * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == ofc)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                fc * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLSAI)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                A->L->vegpar[jdLSAI] * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oz0v)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                z0veg * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == od0v)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                d0veg * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oEcan)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                (SWv_vis + SWv_nir + LWv - Hv - LEv) * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oSWv)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                (SWv_vis + SWv_nir) * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLWv)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                LWv * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oHv)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Hv * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLEv)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                LEv * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oHg0)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Hg0 * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLEg0)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Turbulence::Levap(*Tgskin) * Eg0 * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oHg1)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Hg1 * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLEg1)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Turbulence::Levap(*Tgskin) * Eg1 * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oevapsur)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Evap_soil;

                            }
                          else if (geotop::common::Variables::opnt[l] == otrasp)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Etrans * Dt;

                            }
                          else if (geotop::common::Variables::opnt[l] == oQv)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Qv * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oQg)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Qg * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oQa)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Qa * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oQs)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Qs * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLobuk)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Lobukhov * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLobukcan)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Locc * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == outop)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                u_top * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == odecay)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                decaycoeff * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oSWup)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                *SWupabove_v * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == oLWup)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                LWupabove_v * Dt /
                                A->P->Dtplot_point[geotop::common::Variables::i_sim];

                            }
                          else if (geotop::common::Variables::opnt[l] == omrsnow)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Melt_snow;

                            }
                          else if (geotop::common::Variables::opnt[l] == osrsnow)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Evap_snow;

                            }
                          else if (geotop::common::Variables::opnt[l] == omrglac)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Melt_glac;

                            }
                          else if (geotop::common::Variables::opnt[l] == osrglac)
                            {
                              geotop::common::Variables::odp
                              [geotop::common::Variables::opnt[l]][A->P->jplot[j] - 1] =
                                Evap_glac;
                            }
                        }
                    }
                }
            }

          if (A->P->state_basin == 1)
            {
              if (A->P->Dtplot_basin[geotop::common::Variables::i_sim] > 0)
                {
                  for (l = 0; l < geotop::common::Variables::nobsn; l++)
                    {
                      if (geotop::common::Variables::obsn[l] > 0)
                        {
                          if (geotop::common::Variables::obsn[l] == ooprecrain)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                Prain / (double)A->P->total_pixel;
                            }
                          else if (geotop::common::Variables::obsn[l] == ooprecsnow)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                Psnow / (double)A->P->total_pixel;
                            }
                          else if (geotop::common::Variables::obsn[l] == oorainover)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                Prain_over / (double)A->P->total_pixel;
                            }
                          else if (geotop::common::Variables::obsn[l] == oosnowover)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                Psnow_over / (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooTa)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                Tpoint *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooTg)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                *Tgskin *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooTv)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                V->Tv[j] *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;
                            }
                          else if (geotop::common::Variables::obsn[l] == ooevapsur)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                Evap_soil / (double)A->P->total_pixel;
                            }
                          else if (geotop::common::Variables::obsn[l] == ootrasp)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                (Etrans * Dt) / (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooLE)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                LE *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooH)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                H *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooSW)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                SW *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooLW)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                LW *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooLEv)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                LEv *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooHv)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                Hv *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooSWv)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                (SWv_vis + SWv_nir) *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooLWv)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                LWv *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooSWin)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                SWin *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;

                            }
                          else if (geotop::common::Variables::obsn[l] == ooLWin)
                            {
                              geotop::common::Variables::odb
                              [geotop::common::Variables::obsn[l]] +=
                                LWin *
                                (Dt / A->P->Dtplot_basin[geotop::common::Variables::i_sim]) /
                                (double)A->P->total_pixel;
                            }
                        }
                    }
                }
            }
        }

      return 0;

    }
  else
    {
#ifdef VERBOSE
      lg->logf(
        "PointEnergyBalance not converging. surfacemelting=%d, time=%f, Dt=%f, "
        "i=%ld, j=%ld, r=%ld, c=%ld, elev=%f, aspect=%f, slope=%f, soil type=%d, "
        "land cover=%d "
        "snowD=%f, ns=%ld, ng=%ld, zmu=%f, zmT=%f, z0s=%f, d0s=%f, "
        "rz0s=%f,z0v=%f, d0v=%f, rz0v=%f,"
        "hveg=%f, WindPoint=%f, AirTPoint=%f, Qa=%f, AirPPoint=%f, "
        "LapseRateAirT=%f, eps=%f, fc=%f, LSAI=%f, decaycoeff0=%f, Wcrn=%f,\n "
        "Wcrnmax=%f, Wcsn=%f, Wcsnmax=%f, SWin=%f, LWin=%f, SWv=%f, LW=%f, H=%f, "
        "E=%f, LWv=%f, Hv=%f,"
        "LEv=%f, Etrans=%f, Ts=%f, Qs=%f, Hadv=%f, Hg0=%f, Hg1=%f, Eg0=%f, "
        "Eg1=%f, Qv=%f, Qg=%f,"
        "Lobukhov=%f, rh=%f, rv=%f, rb=%f, rc=%f, ruc=%f, u_top=%f, decay=%f, "
        "Locc=%f, LWupabove_v=%f, lpb=%ld\n",
        surface, JDb - A->P->init_date, Dt, i, j, r, c, A->T->Z0[r][c],
        A->T->aspect[r][c], A->T->slope[r][c], (int)A->L->ty[lu][jz0],
        (int)A->L->LC[r][c], snowD, ns, ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg,
        d0veg, 1.0, hveg, Vpoint, Tpoint, Qa, Ppoint, A->M->LRv[ilsTa], eps, fc,
        A->L->vegpar[jdLSAI], A->L->vegpar[jddecay0], V->wrain[j], max_wcan_rain,
        V->wsnow[j], max_wcan_snow, SWin, LWin, SWv_vis + SWv_nir, LW, H, E, LWv,
        Hv, LEv, Etrans, Ts, Qs, Hadv, Hg0, Hg1, Eg0, Eg1, Qv, Qg, Lobukhov, rh,
        rv, rb, rc, ruc, u_top, decaycoeff, Locc, LWupabove_v, lpb);
#endif  // VERBOSE

      return 1;
    }
}

// end of "PointEnergyBalance" subroutine

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short SolvePointEnergyBalance(

  const short surfacemelting,
  double Tgd,
  double EBd,
  double Convd,
  short surfacebalance,
  const double t,
  const double Dt,
  const long i,
  const long j,
  const long r,
  const long c,
  SoilState *SL,
  SoilState *SC,
  StateVeg *V,
  Energy *egy,
  Land *land,
  Soil *sl,
  Channel *cnet,
  Topo *topog,
  Par *par,
  const long ns,
  const long ng,
  const double zmu,
  const double zmT,
  const double z0s,
  const double d0s,
  const double rz0s,
  const double z0v,
  const double d0v,
  const double rz0v,
  const double hveg,
  const double v,
  const double Ta,
  const double Qa,
  const double P,
  const double LR,
  const double eps,
  const double fc,
  const double LSAI,
  const double decaycoeff0,
  double *Wcrn,
  const double Wcrnmax,
  double *Wcsn,
  const double Wcsnmax,
  const double SWin,
  const double LWin,
  const double SWv,
  double *LW,
  double *H,
  double *E,
  double *LWv,
  double *Hv,
  double *LEv,
  double *Etrans,
  double *Ts,
  double *Qs,
  const double Eadd,
  double *Hg0,
  double *Hg1,
  double *Eg0,
  double *Eg1,
  double *Qv,
  double *Qg,
  double *Lob,
  double *rh,
  double *rv,
  double *rb,
  double *rc,
  double *ruc,
  double *u_top,
  double *decay,
  double *Locc,
  double *LWup_ab_v,
  long *lpb,
  double *dUsl)
{
  short iter_close, iter_close2, lu = land->LC[r][c], flagTmin = 0, sux;
  short dirichlet = 0, neumann = 0, micro_sempl = 0, micro = 0;
  long sur, sy, l, m, cont = 0, cont2,
                      n = geotop::common::Variables::Nl + ns + ng,
                      cont_lambda_min = 0;
  double dH_dT, dE_dT, EB = 0.0, dEB_dT = 0.0, EB0, Tg, Tg0, psim0, psi0,
                       Qg0 = 0.0, Tv0 = 0.0, dWcsn = 0.0, dWcrn = 0.0, rh_g,
                       rv_g;
  double res, res0[3], res_av, res_prev[MM], lambda[3], C0, C1, th0, th1, kbb0,
         kbb1, kub0 = 0.0, kub1 = 0, thi = 0.0, thin = 0.0, thw = 0.0, thwn = 0.0,
               sat = 0.0, satn = 0.0, kt = 0.0, ktn = 0.0;
  (void)kub0;  // silence a warning about unused variable, since kub0 appears
  // only in a if branch
  double snowD_tmp = 0.0;  // to make it compile..
  // SnowD is NOT used at all in this routine: it is taken from above and passed
  // to energyfluxes routines to be printed out in case of NaN is happening on
  // some variable...

  FILE *f;
  //    geotop::logger::GlobalLogger* lg =
  //    geotop::logger::GlobalLogger::getInstance();

  //

  res0[0] = 0.;
  res0[1] = 0.;
  res0[2] = 0.;
  lambda[0] = 0.;
  lambda[1] = 0.;
  lambda[2] = 0.;

  //
  // boundary condition

  if ((long)Tgd != geotop::input::gDoubleNoValue) dirichlet = 1;
  if ((long)EBd != geotop::input::gDoubleNoValue) neumann = 1;
  if ((long)Convd != geotop::input::gDoubleNoValue) micro_sempl = 1;
  if (neumann == 0 && dirichlet == 0 && micro_sempl == 0) micro = 1;

  // Soil layer; this is not working...
  // n = Fminlong(par->Nl_spinup->[geotop::common::Variables::i_sim],Nl) + ns +
  // ng;

  // Surface conditions
  sur = par->nsurface;
  if (surfacemelting == 1) sur = 1;
  if (neumann == 1) EB = EBd;
  if (micro_sempl == 1) dEB_dT = -Convd;

  if (dirichlet == 1)
    {
      sur = 1;
      egy->Temp[0] = Tgd;
    }

  // Initializes vectors
  for (l = sur; l <= n; l++)
    {
      egy->T0[l] = egy->Temp[l];
    }

  // Initialize soil properties
  if (i <= long(par->total_channel))
    {
      sy = cnet->soil_type[j];
      psi0 = SC->P[0][j];
      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          //  water content to be modified at each iteration
          egy->THETA[l] = cnet->th[l][j];

          // total pressure (=pressure of the water would have if it was all in
          // liquid state)
          //
          psim0 = psi_teta(cnet->th[l][j] + SC->thi[l][j], 0.0, sl->pa[sy][jsat][l],
                           sl->pa[sy][jres][l], sl->pa[sy][ja][l],
                           sl->pa[sy][jns][l], 1 - 1 / sl->pa[sy][jns][l],
                           GTConst::PsiMin, sl->pa[sy][jss][l]);
          // max temperature at which the first particle of ice comes up

          egy->Tstar[l] =
            Fmin(psim0 / (1000.0 * GTConst::Lf /
                          (GTConst::GRAVITY * (GTConst::Tfreezing + GTConst::tk))),
                 0.0);
        }
    }
  else
    {
      sy = sl->type[r][c];
      psi0 = SL->P[0][j];
      for (l = 1; l <= geotop::common::Variables::Nl; l++)
        {
          //  water content to be modified at each iteration

          egy->THETA[l] = sl->th[l][j];

          //  total pressure (=pressure of the water would have if it was all in
          //liquind state)

          psim0 = psi_teta(sl->th[l][j] + SL->thi[l][j], 0.0, sl->pa[sy][jsat][l],
                           sl->pa[sy][jres][l], sl->pa[sy][ja][l],
                           sl->pa[sy][jns][l], 1 - 1 / sl->pa[sy][jns][l],
                           GTConst::PsiMin, sl->pa[sy][jss][l]);

          //  max temperature at which the first particle of ice comes up

          egy->Tstar[l] =
            Fmin(psim0 / (1000.0 * GTConst::Lf /
                          (GTConst::GRAVITY * (GTConst::Tfreezing + GTConst::tk))),
                 0.0);
        }
    }

  // Calculate the surface energy balance only if you need it
  if (surfacebalance == 1)
    {
      //  Surface skin temperature used to compute the surface energy fluxes

      Tg = egy->Temp[sur];
      Tg0 = Tg;

      // Qg0 is the specific humidity at the soil surface at the beginning of the
      // time step (assumed saturated)

      Qg0 = SpecHumidity(SatVapPressure(Tg0, P), P);

      // canopy temperature at the beginning of the time step

      Tv0 = V->Tv[j];

      /*calculates all the surface energy fluxes, includes the calculation of
       vegetation temperature (solving the vegetation energy balance) Returns:

       LW: net longwave radiation at the surface
       H: sensible heat flux at the surface (positive upwards)
       E: evaporation from the soil surface (mm/s)

       LWv: net longwave in the canopy
       Hv: sensible heat flux from the canopy
       LEv: Turbulence::latent heat flux from the canopy
       Etrans: canopy transpiration (as water flux in mm/s)

       sl->Tv->co[r][c] : vegetation temperature
       Qv : specific humidity at the canopy (assumed saturated)
       Ts : temperature of canopy air
       Qs : specific humidity of canopy air*/

      // surface energy balance

      // in this routine there are many variable not used to actually computing
      // something out of them. We could get rid of them..

      EnergyFluxes(
        t, Tg, r, c, ns + ng, Tg0, Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v, d0v,
        rz0v, hveg, v, Ta, Qa, P, LR, psi0, eps, fc, LSAI, decaycoeff0, *Wcrn,
        Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, &(egy->THETA[0]), sl->pa, sy,
        land->ty, lu, land->root_fraction, par, egy->soil_transp_layer, SWin,
        LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans, &(V->Tv[j]),
        Qv, Ts, Qs, Hg0, Hg1, Eg0, Eg1, Lob, rh, rv, rc, rb, ruc, &rh_g, &rv_g,
        Qg, u_top, decay, Locc, LWup_ab_v, &(egy->Temp[0]),
        egy->soil_evap_layer_bare, egy->soil_evap_layer_veg, topog->Z0[r][c],
        topog->slope[r][c], topog->aspect[r][c], topog->sky[r][c], land->LC[r][c],
        sl->type[r][c], snowD_tmp);

      if (micro == 1)
        {
          EB = *LW - (*H) - Turbulence::latent(Tg, Turbulence::Levap(Tg)) * (*E) +
               Eadd + egy->SWlayer[0];
          dEB_dT = -eps * topog->sky[r][c] * dSB_dT(Tg) - dH_dT -
                   Turbulence::latent(Tg, Turbulence::Levap(Tg)) *
                   dE_dT;  // ADDED top->sky[r][c] term
        }
    }

  if (micro_sempl == 1) EB = Convd * (Ta - egy->Temp[sur]);

  EB0 = EB;

  // Calculate -F(x0), given that the system to solve is F'(x0) * (x-x0) =
  // -F(x0)  F(T) = diag(egy->Fenergy(T)) + K(T)*T

  for (l = sur; l <= n; l++)
    {
      egy->Fenergy[l] = 0.0;
      if (l > 0) egy->deltaw[l] = 0.0;
    }

  for (l = 1; l <= n; l++)
    {
      // conductive M-matrix (lower diagonal part of this M-matrix is stored in a
      // vector)
      thw = (egy->liq[l] + egy->deltaw[l]) / (GTConst::rho_w * egy->Dlayer[l]);
      thi = (egy->ice[l] - egy->deltaw[l]) / (GTConst::rho_i * egy->Dlayer[l]);

      if (l > ns + ng)
        {
          sat = sl->pa[sy][jsat][l - ns - ng];
          kt = sl->pa[sy][jkt][l - ns - ng];
        }
      else
        {
          sat = 1.;
          kt = 0.;
        }

      // from UZH: condition on l=1
      //
      if (l == 1)
        {
          if (l <= ns + ng)
            {
              kub1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
            }
          else
            {
              kub1 = k_thermal(0, 1, thw, thi, sat, kt);
            }
          kub0 = kub1;
        }

      if (l > 1)
        {
          thwn = (egy->liq[l - 1] + egy->deltaw[l - 1]) /
                 (GTConst::rho_w * egy->Dlayer[l - 1]);
          thin = (egy->ice[l - 1] - egy->deltaw[l - 1]) /
                 (GTConst::rho_i * egy->Dlayer[l - 1]);

          if (l - 1 > ns + ng)    // if the layer (l) is <= num_of_snow_layers(ns) +
            {
              // num_of_glacies(ng)
              satn = sl->pa[sy][jsat][l - 1 - ns - ng];
              ktn = sl->pa[sy][jkt][l - 1 - ns - ng];
            }
          else
            {
              satn = 1.;
              ktn = 0.;
            }

          thwn = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], thw, thwn);
          thin = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], thi, thin);
          satn = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], sat, satn);
          ktn = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], kt, ktn);

          if (l <= ns + ng)    // if the layer (l) is <= num_of_snow_layers(ns) +
            {
              // num_of_glacies(ng)
              egy->Kth1[l - 1] =
                -2. * k_thermal(1, par->snow_conductivity, thwn, thin, satn, ktn) /
                (egy->Dlayer[l - 1] + egy->Dlayer[l]);
            }
          else
            {
              egy->Kth1[l - 1] = -2. * k_thermal(0, 1, thwn, thin, satn, ktn) /
                                 (egy->Dlayer[l - 1] + egy->Dlayer[l]);
            }

        }
      else if (l > sur)
        {
          if (l <= ns + ng)
            {
              egy->Kth1[l - 1] =
                -k_thermal(1, par->snow_conductivity, thw, thi, sat, kt) /
                (egy->Dlayer[l] / 2.);
            }
          else
            {
              egy->Kth1[l - 1] =
                -k_thermal(0, 1, thw, thi, sat, kt) / (egy->Dlayer[l] / 2.);
            }
        }

      egy->Kth0[l - 1] = egy->Kth1[l - 1];

      // diagonal part of F due to heat capacity and boundary condition (except
      // conduction)
      C0 = calc_C(l, ns + ng, 0.0, &(egy->ice[0]), &(egy->liq[0]),
                  &(egy->deltaw[0]), &(egy->Dlayer[0]), sl->pa, sy);
      C1 = calc_C(l, ns + ng, 1.0, &(egy->ice[0]), &(egy->liq[0]),
                  &(egy->deltaw[0]), &(egy->Dlayer[0]), sl->pa, sy);
      if (l <= ns + ng && egy->ice[l] - egy->deltaw[l] < 1.E-7)
        C1 = Csnow_at_T_greater_than_0;
      egy->Fenergy[l] +=
        (GTConst::Lf * egy->deltaw[l] + C1 * egy->Dlayer[l] * egy->Temp[l] -
         C0 * egy->Dlayer[l] * egy->T0[l]) /
        Dt;

      // shortwave radiation penetrating under the surface
      if (micro == 1 && l <= ns + 1) egy->Fenergy[l] -= egy->SWlayer[l];
    }

#ifdef VERBOSE
  lg->logsf(geotop::logger::DEBUG, "energy-balance: NS(%ld)\n", ns);

  // top boundary condition

#endif
  egy->Fenergy[sur] -= ((1. - GTConst::KNe) * EB + GTConst::KNe * EB0);
  // TODO: removed 4/11/2013 egy->Fenergy[sur] -= egy->SWlayer[0];
  //
  // bottom boundary condition (treated as sink)
  if (n <= ns + ng)
    {
      kbb1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
    }
  else
    {
      kbb1 = k_thermal(0, 1, thw, thi, sat, kt);
    }
  kbb0 = kbb1;
  egy->Fenergy[n] -=
    ((par->Tboundary - egy->Temp[n]) * (1. - GTConst::KNe) * kbb1 +
     (par->Tboundary - egy->T0[n]) * GTConst::KNe * kbb0) /
    (egy->Dlayer[n] / 2. + par->Zboundary);
  egy->Fenergy[n] -= par->Fboundary;

  // include conduction in F(x0)
  update_F_energy(sur, n, egy->Fenergy, 1. - GTConst::KNe, egy->Kth1,
                  &(egy->Temp[0]));
  update_F_energy(sur, n, egy->Fenergy, GTConst::KNe, egy->Kth0, &(egy->T0[0]));

  // calculate the norm
  res = norm_2(egy->Fenergy, sur, n);

  do
    {
      cont++;

      for (l = sur; l <= n; l++)
        {
          egy->T1[l] = egy->Temp[l];
        }

      // Calculates F'(x) referred to as Jacobian matrix

      // First the diagonal part is calculated and put into the vector dFenergy

      for (l = sur; l <= n; l++)
        {
          egy->dFenergy[l] = 0.0;
        }

      // Heat capacity part
      for (l = 1; l <= n; l++)
        {
          // real thermal capacity
          C1 = calc_C(l, ns + ng, 1.0, &(egy->ice[0]), &(egy->liq[0]),
                      &(egy->deltaw[0]), &(egy->Dlayer[0]), sl->pa, sy);
          if (l <= ns + ng && egy->ice[l] - egy->deltaw[l] < 1.E-7)
            C1 = Csnow_at_T_greater_than_0;
          // adds apparent thermal conductivity (due to phase change) for both snow
          // and soil
          if (l <= ns + ng)    // snow
            {
              C1 += GTConst::Lf * (egy->ice[l] + egy->liq[l]) *
                    dtheta_snow(par->alpha_snow, 1., egy->Temp[l]) / egy->Dlayer[l];
            }
          else      // soil
            {
              m = l - ns - ng;

              if (egy->Temp[l] <= egy->Tstar[m])
                {
                  C1 += GTConst::rho_w * GTConst::Lf *
                        (GTConst::Lf / (GTConst::GRAVITY * GTConst::tk) * 1.E3) *
                        dteta_dpsi(Psif(egy->Temp[l]), 0.0, sl->pa[sy][jsat][m],
                                   sl->pa[sy][jres][m], sl->pa[sy][ja][m],
                                   sl->pa[sy][jns][m], 1 - 1 / sl->pa[sy][jns][m],
                                   GTConst::PsiMin, 0.0);
                }
            }
          egy->dFenergy[l] += C1 * egy->Dlayer[l] / Dt;
        }

      // Upper Boundary Condition: Neumann
      egy->dFenergy[sur] -= (1. - GTConst::KNe) * dEB_dT;

      // Upper Boundary Condition: Dirichlet
      if (dirichlet == 1)
        egy->dFenergy[sur] += (1. - GTConst::KNe) * kub1 / (egy->Dlayer[1] / 2.);

      // Bottom Boundary Condition
      egy->dFenergy[n] +=
        (1. - GTConst::KNe) * kbb1 / (egy->Dlayer[n] / 2. + par->Zboundary);

      // Update the part of Jacobian due to conduction is included in Kth1 and
      // Kth0
      update_diag_dF_energy(sur, n, egy->dFenergy, 1. - GTConst::KNe, egy->Kth1);

      // Extradiagonal part of the Jacobian
      for (l = sur; l <= n - 1; l++)
        {
          egy->udFenergy[l] = (1. - GTConst::KNe) * egy->Kth1[l];
        }
      // Solve tridiagonal system
      sux = tridiag2(sur, n, egy->udFenergy, egy->dFenergy, egy->udFenergy,
                     egy->Fenergy, egy->Newton_dir);

      if (sux == 1)
        {
          if (ns != 1)
            {
              return 1;
            }
          else if (ng > 0)
            {
              return -5;
            }
          else
            {
              return -6;
            }
        }

      cont2 = 0;
      iter_close = 0;
      res0[0] = res;

      // non-monotonic line search (it is monotonic if M==1)
      for (m = Fminlong(cont, MM); m > 1; m--)
        {
          res_prev[m - 1] = res_prev[m - 2];
        }
      res_prev[0] = res;

      res_av = 0.0;
      for (m = 1; m <= Fminlong(cont, MM); m++)
        {
          res_av = Fmax(res_prev[m - 1], res_av);
        }

      do
        {
          cont2++;
          iter_close2 = 0;

          if (cont2 == 1)
            {
              lambda[0] = 1.0;

            }
          else if (cont2 == 2)
            {
              lambda[1] = lambda[0];
              res0[1] = res;
              lambda[0] = GTConst::thmax;

            }
          else
            {
              lambda[2] = lambda[1];
              res0[2] = res0[1];
              lambda[1] = lambda[0];
              res0[1] = res;
              lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1],
                                                  lambda[2], res0[2]);
            }
          for (l = sur; l <= n; l++)
            {
              egy->Temp[l] = egy->T1[l] + lambda[0] * egy->Newton_dir[l];
              // soil
              if (l > ns + ng)
                {
                  m = l - ns - ng;

                  if (i <= long(par->total_channel))
                    {
                      th0 = cnet->th[m][j];

                      th1 = teta_psi(
                              Psif(Fmin(egy->Tstar[m], egy->Temp[l])), 0.0, sl->pa[sy][jsat][m],
                              sl->pa[sy][jres][m], sl->pa[sy][ja][m], sl->pa[sy][jns][m],
                              1 - 1 / sl->pa[sy][jns][m], GTConst::PsiMin, sl->pa[sy][jss][m]);

                      if (th1 > cnet->th[m][j] + SC->thi[m][j])
                        th1 = cnet->th[m][j] + SC->thi[m][j];
                    }
                  else
                    {
                      th0 = sl->th[m][j];

                      th1 = teta_psi(
                              Psif(Fmin(egy->Tstar[m], egy->Temp[l])), 0.0, sl->pa[sy][jsat][m],
                              sl->pa[sy][jres][m], sl->pa[sy][ja][m], sl->pa[sy][jns][m],
                              1 - 1 / sl->pa[sy][jns][m], GTConst::PsiMin, sl->pa[sy][jss][m]);

                      if (th1 > sl->th[m][j] + SL->thi[m][j])
                        th1 = sl->th[m][j] + SL->thi[m][j];
                    }

                  egy->deltaw[l] = (th1 - th0) * egy->Dlayer[l] * GTConst::rho_w;

                  // snow
                }
              else if (l > 0)
                {
                  th0 = egy->liq[l] / (egy->ice[l] + egy->liq[l]);

                  th1 = theta_snow(par->alpha_snow, 1., egy->Temp[l]);

                  egy->deltaw[l] = (th1 - th0) * (egy->ice[l] + egy->liq[l]);
                  //                    printf ("DEBUG snow:%ld,%f\n", l,
                  //                    egy->deltaw[l]);
                }
            }

          Tg = egy->Temp[sur];
          if (Tg < Tmin_surface_below_which_surfenergy_balance_recalculated)
            flagTmin++;

          // this condition from UZH
          if (surfacebalance == 1)
            {
              if (cont < num_iter_after_which_surfenergy_balance_not_recalculated ||
                  flagTmin > 0)
                {
                  // update egy->THETA taking into account evaporation (if there is not
                  // snow)
                  if (ns + ng == 0)
                    {
                      for (size_t ll = ns + ng + 1; ll <= (size_t)n; ll++)
                        {
                          m = ll - ns - ng;

                          if (i <= long(par->total_channel))
                            {
                              egy->THETA[m] =
                                cnet->th[m][j] +
                                egy->deltaw[ll] / (GTConst::rho_w * egy->Dlayer[ll]);
                            }
                          else
                            {
                              egy->THETA[m] =
                                sl->th[m][j] +
                                egy->deltaw[ll] / (GTConst::rho_w * egy->Dlayer[ll]);
                            }

                          // add canopy transpiration
                          if (egy->THETA[m] > sl->pa[sy][jres][1] + 1.E-3 &&
                              ll <= egy->soil_transp_layer.size())
                            {
                              egy->THETA[m] -= Fmax(Dt * fc * egy->soil_transp_layer[m] /
                                                    (GTConst::rho_w * egy->Dlayer[ll]),
                                                    0.0);
                              if (egy->THETA[m] < sl->pa[sy][jres][m] + 1.E-3)
                                egy->THETA[m] = sl->pa[sy][jres][m] + 1.E-3;
                            }

                          // add soil evaporation

                          if (egy->THETA[m] > sl->pa[sy][jres][1] + 1.E-3 &&
                              ll < egy->soil_evap_layer_bare.size())
                            {
                              egy->THETA[m] -=
                                Fmax(Dt * (1. - fc) * egy->soil_evap_layer_bare[m] /
                                     (GTConst::rho_w * egy->Dlayer[ll]),
                                     0.0);

                              egy->THETA[m] -= Fmax(Dt * fc * egy->soil_evap_layer_veg[m] /
                                                    (GTConst::rho_w * egy->Dlayer[ll]),
                                                    0.0);
                              if (egy->THETA[m] < sl->pa[sy][jres][m] + 1.E-3)
                                egy->THETA[m] = sl->pa[sy][jres][m] + 1.E-3;
                            }
                        }
                    }

                  // surface energy balance

                  EnergyFluxes_no_rec_turbulence(
                    t, Tg, r, c, ns + ng, egy->T0[1], Qg0, Tv0, zmu, zmT, z0s, d0s,
                    rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, P, LR, SL->P[0][j], eps, fc,
                    LSAI, decaycoeff0, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn,
                    &(egy->THETA[0]), sl->pa, sy, land->ty, lu, land->root_fraction,
                    par, egy->soil_transp_layer, SWin, LWin, SWv, LW, H, &dH_dT, E,
                    &dE_dT, LWv, Hv, LEv, Etrans, &(V->Tv[j]), Qv, Ts, Qs, Hg0, Hg1,
                    Eg0, Eg1, Lob, rh, rv, rc, rb, ruc, &rh_g, &rv_g, Qg, u_top, decay,
                    Locc, LWup_ab_v, &(egy->Temp[0]), egy->soil_evap_layer_bare,
                    egy->soil_evap_layer_veg, topog->sky[r][c], flagTmin, cont);

                  if (micro == 1)
                    {
                      EB = *LW - (*H) -
                           Turbulence::latent(Tg, Turbulence::Levap(Tg)) * (*E) + Eadd +
                           egy->SWlayer[0];
                      dEB_dT =
                        -eps * topog->sky[r][c] * dSB_dT(Tg) - dH_dT -
                        Turbulence::latent(Tg, Turbulence::Levap(Tg)) *
                        dE_dT;  // ADDED top->sky[r][c] term to match uzh exact results
                    }
                }
            }  // endif for surfacebalance==1

          if (micro_sempl == 1) EB = Convd * (Ta - egy->Temp[sur]);

          // F(T) = diag(egy->Fenergy(T)) + K(T)*T
          for (l = sur; l <= n; l++)
            {
              egy->Fenergy[l] = 0.0;
            }

          *dUsl = 0.0;

          for (l = 1; l <= n; l++)
            {
              // conductive M-matrix (lower diagonal part of this M-matrix is stored
              // in a vector)

              thw =
                (egy->liq[l] + egy->deltaw[l]) / (GTConst::rho_w * egy->Dlayer[l]);
              thi =
                (egy->ice[l] - egy->deltaw[l]) / (GTConst::rho_i * egy->Dlayer[l]);

              if (l > ns + ng)
                {
                  sat = sl->pa[sy][jsat][l - ns - ng];
                  kt = sl->pa[sy][jkt][l - ns - ng];
                }
              else
                {
                  sat = 1.;
                  kt = 0.;
                }

              // ADDEDED from UZH

              if (l <= ns + ng)
                {
                  kub1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
                }
              else
                {
                  kub1 = k_thermal(0, 1, thw, thi, sat, kt);
                }

              if (l > 1)
                {
                  thwn = (egy->liq[l - 1] + egy->deltaw[l - 1]) /
                         (GTConst::rho_w * egy->Dlayer[l - 1]);
                  thin = (egy->ice[l - 1] - egy->deltaw[l - 1]) /
                         (GTConst::rho_i * egy->Dlayer[l - 1]);

                  if (l - 1 > ns + ng)
                    {
                      satn = sl->pa[sy][jsat][l - 1 - ns - ng];
                      ktn = sl->pa[sy][jkt][l - 1 - ns - ng];
                    }
                  else
                    {
                      satn = 1.;
                      ktn = 0.;
                    }

                  thwn = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], thw, thwn);
                  thin = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], thi, thin);
                  satn = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], sat, satn);
                  ktn = Arithmetic_Mean(egy->Dlayer[l], egy->Dlayer[l - 1], kt, ktn);

                  if (l <= ns + ng)
                    {
                      egy->Kth1[l - 1] =
                        -2. *
                        k_thermal(1, par->snow_conductivity, thwn, thin, satn, ktn) /
                        (egy->Dlayer[l - 1] + egy->Dlayer[l]);
                    }
                  else
                    {
                      egy->Kth1[l - 1] = -2. * k_thermal(0, 1, thwn, thin, satn, ktn) /
                                         (egy->Dlayer[l - 1] + egy->Dlayer[l]);
                    }

                }
              else if (l > sur)
                {
                  if (l <= ns + ng)
                    {
                      egy->Kth1[l - 1] =
                        -k_thermal(1, par->snow_conductivity, thw, thi, sat, kt) /
                        (egy->Dlayer[l] / 2.);
                    }
                  else
                    {
                      egy->Kth1[l - 1] =
                        -k_thermal(0, 1, thw, thi, sat, kt) / (egy->Dlayer[l] / 2.);
                    }
                }

              // diagonal part of F due to heat capacity and boundary condition
              // (except conduction)

              C0 = calc_C(l, ns + ng, 0.0, &(egy->ice[0]), &(egy->liq[0]),
                          &(egy->deltaw[0]), &(egy->Dlayer[0]), sl->pa, sy);
              C1 = calc_C(l, ns + ng, 1.0, &(egy->ice[0]), &(egy->liq[0]),
                          &(egy->deltaw[0]), &(egy->Dlayer[0]), sl->pa, sy);

              if (l <= ns + ng && egy->ice[l] - egy->deltaw[l] < 1.E-7)
                C1 = Csnow_at_T_greater_than_0;
              egy->Fenergy[l] +=
                (GTConst::Lf * egy->deltaw[l] + C1 * egy->Dlayer[l] * egy->Temp[l] -
                 C0 * egy->Dlayer[l] * egy->T0[l]) /
                Dt;

              // shortwave radiation penetrating under the surface

              if (micro == 1 && l <= ns + 1) egy->Fenergy[l] -= egy->SWlayer[l];

              // store soil internal energy
              if (l > ns + ng)
                *dUsl = *dUsl + GTConst::Lf * egy->deltaw[l] +
                        C1 * egy->Dlayer[l] * egy->Temp[l] -
                        C0 * egy->Dlayer[l] * egy->T0[l];
            }

          // top boundary condition

          egy->Fenergy[sur] -= ((1. - GTConst::KNe) * EB + GTConst::KNe * EB0);

          //     if(dirichlet == 1) egy->Fenergy[sur] -= (
          //(Tgd-egy->Temp[sur])*(1.-KNe)*kub1 + (Tgd-egy->T0[sur])*KNe*kub0 ) /
          //(egy->Dlayer[1]/2.);

          // TODO: remove 04/11/2013: egy->Fenergy[sur] -= egy->SWlayer[0];

          // bottom boundary condition (treated as sink)
          if (n <= ns + ng)
            {
              kbb1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
            }
          else
            {
              kbb1 = k_thermal(0, 1, thw, thi, sat, kt);
            }
          egy->Fenergy[n] -=
            ((par->Tboundary - egy->Temp[n]) * (1. - GTConst::KNe) * kbb1 +
             (par->Tboundary - egy->T0[n]) * GTConst::KNe * kbb0) /
            (egy->Dlayer[n] / 2. + par->Zboundary);
          egy->Fenergy[n] -= par->Fboundary;

          update_F_energy(sur, n, egy->Fenergy, 1. - GTConst::KNe, egy->Kth1,
                          &(egy->Temp[0]));
          update_F_energy(sur, n, egy->Fenergy, GTConst::KNe, egy->Kth0,
                          &(egy->T0[0]));
          res = norm_2(egy->Fenergy, sur, n);

          if (res <= res_av * (1.0 - ni_en * lambda[0])) iter_close2 = 1;
          if (lambda[0] <= par->min_lambda_en) cont_lambda_min++;
          if (cont_lambda_min > par->max_times_min_lambda_en)
            {
              if (par->exit_lambda_min_en == 1)
                {
                  return 1;
                }
              else
                {
                  iter_close2 = 1;
                  cont_lambda_min = 0;
                }
            }

        }
      while (iter_close2 != 1);

      if (iter_close >= 0 && res <= par->tol_energy) iter_close = 1;
      if (cont >= par->maxiter_energy) iter_close = 1;

    }
  while (iter_close != 1);

  // if there is no convergence, go out of the loop
  if (res > par->tol_energy)
    {
      if (ns > 1 && surfacemelting == 0)
        {
          return -1;
        }
      else if (ns == 1 && ng == 0)
        {
          return -6;
        }
      else if (ns == 1 && ng > 0)
        {
          return -5;
        }
      else
        {
          return 1;
        }
    }

  //  skin temperature

  if (sur == 1 && par->nsurface == 0) egy->Temp[par->nsurface] = egy->Temp[sur];

  //  it is not admitted that the skin layer temperature goes > 0 when snow or
  //glacier are present

  if (sur == 0 && ns + ng > 0 && egy->Temp[sur] > 0) return -1;

  //  it is not admitted that a snow/glacier layer melts completely
  for (l = 1; l <= ns + ng; l++)
    {
      //  snow layer not alone

      if (ns > 1 && l >= 1 && l <= ns && egy->ice[l] - egy->deltaw[l] < 1.E-5)
        {
          *lpb = l;
          return -2;

          //  glacier layer not alone
          //  }else if (ng > 1 && l >= ns+1 && l<= ns+ng && egy->ice->co[l] -
          //egy->deltaw->co[l] < 1.E-5){
        }
      else if (ng > 1 && l >= ns + 1 && l <= ns + ng &&
               egy->ice[l] - egy->deltaw[l] < 1.E-5)
        {
          *lpb = l;
          return -3;

          //  just one glacier layer and snow layers >= 0

        }
      else if (ng == 1 && l == ns + ng &&
               egy->ice[l] - egy->deltaw[l] < 1.E-5)
        {
          return -4;

          //  just one snow layer and glacier layers > 0

        }
      else if (ns == 1 && ng > 0 && l == ns &&
               egy->ice[l] - egy->deltaw[l] < 1.E-5)
        {
          return -5;

          //  just one snow layer and glacier layers == 0

        }
      else if (ns == 1 && ng == 0 && l == ns &&
               egy->ice[l] - egy->deltaw[l] < 1.E-5)
        {
          return -6;
        }
    }

  // error messages
  for (l = sur; l <= geotop::common::Variables::Nl + ns + ng; l++)
    {
      // no value conditions: here T is NaN
      if (egy->Temp[l] != egy->Temp[l])
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t);

          fprintf(f,
                  "T no value, PointEnergyBalance, l:%ld r:%ld c:%ld elev=%f "
                  "slope=%f aspect=%f sky=%f LC=%d SoilType=%ld\n",
                  l, r, c, topog->Z0[r][c], topog->slope[r][c], topog->aspect[r][c],
                  topog->sky[r][c], (int)land->LC[r][c], sl->type[r][c]);

          fclose(f);
          t_error("Fatal Error! GEOtop is closed. See failing report.");
        }

      if (egy->Temp[l] < Tmin || egy->Temp[l] > Tmax)
        {
          printf(
            "Warning, T outside of range, PointEnergyBalance, l:%ld r:%ld c:%ld "
            "elev=%f slope=%f aspect=%f sky=%f LC=%d SoilType=%ld snowD:%f T:%f "
            "LW:%f H:%f LE:%f Ta:%f\n",
            l, r, c, topog->Z0[r][c], topog->slope[r][c], topog->aspect[r][c],
            topog->sky[r][c], (int)land->LC[r][c], sl->type[r][c], snowD_tmp,
            egy->Temp[l], *LW, *H,
            Turbulence::latent(Tg, Turbulence::Levap(Tg)) * (*E), Ta);
        }
    }

  for (l = 1; l <= geotop::common::Variables::Nl + ns + ng; l++)
    {
      //  if(egy->deltaw->co[l] != egy->deltaw->co[l]){
      if (egy->deltaw[l] != egy->deltaw[l])
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t);
          fprintf(f,
                  "dw no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",
                  l, r, c);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }
    }

  // update canopy water
  *Wcrn = *Wcrn + dWcrn;
  *Wcsn = *Wcsn + dWcsn;

  return 0;
}

// end of "SolvePointEnergyBalance" subroutine
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_soil_land(long /*nsurf*/,
                      long n,
                      long i,
                      long r,
                      long c,
                      double fc,
                      double Dt,
                      Energy *egy,
                      GeoTensor<double> &pa,
                      long sy,
                      SoilState *S,
                      GeoTensor<double> &ET,
                      GeoMatrix<double> &th)
{
  // long l;//this variable is not used
  double th_oversat, psisat;

  // soil variables
  assert(geotop::common::Variables::Nl >= 0);
  for (size_t l = 1; l <= (size_t)(geotop::common::Variables::Nl); l++)
    {
      //  canopy transpiration
      if (l < egy->soil_transp_layer.size())
        ET[l][r][c] += fc * egy->soil_transp_layer[l] * Dt;

      //  soil evaporation
      if (l < egy->soil_evap_layer_bare.size())
        {
          ET[l][r][c] += (1. - fc) * egy->soil_evap_layer_bare[l] * Dt;
          ET[l][r][c] += fc * egy->soil_evap_layer_veg[l] * Dt;
        }

      //  water pressure and contents
      psisat = psi_saturation(
                 Fmax(0., egy->ice[l + n]) / (GTConst::rho_w * egy->Dlayer[l + n]),
                 pa[sy][jsat][l], pa[sy][jres][l], pa[sy][ja][l], pa[sy][jns][l],
                 1. - 1. / pa[sy][jns][l]);
      th_oversat = Fmax(S->P[l][i] - psisat, 0.0) * pa[sy][jss][l];

      th[l][i] = Fmax(0., egy->liq[l + n] + egy->deltaw[l + n]) /
                 (GTConst::rho_w * egy->Dlayer[l + n]);
      S->thi[l][i] = Fmax(0., egy->ice[l + n] - egy->deltaw[l + n]) /
                     (GTConst::rho_w * egy->Dlayer[l + n]);

      S->P[l][i] =
        psi_teta(th[l][i] + th_oversat, S->thi[l][i], pa[sy][jsat][l],
                 pa[sy][jres][l], pa[sy][ja][l], pa[sy][jns][l],
                 1. - 1. / pa[sy][jns][l], GTConst::PsiMin, pa[sy][jss][l]);

      //  temperature
      S->T[l][i] = egy->Temp[l + n];
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_soil_channel(long /*nsurf*/,
                         long n,
                         long ch,
                         double fc,
                         double Dt,
                         Energy *egy,
                         GeoTensor<double> &pa,
                         long sy,
                         SoilState *S,
                         GeoMatrix<double> &ET,
                         GeoMatrix<double> &th)
{
  // long l;//this variable is not used
  double th_oversat, psisat;

  // soil variables
  for (size_t l = 1; l <= (size_t)(geotop::common::Variables::Nl); l++)
    {
      //  canopy transpiration
      if (l <= egy->soil_transp_layer.size())
        ET[l][ch] += fc * egy->soil_transp_layer[l] * Dt;

      //  soil evaporation
      //  if(l <= egy->soil_evap_layer_bare->nh){
      if (l < egy->soil_evap_layer_bare.size())
        {
          //  ET->co[l][ch] += (1.-fc)*egy->soil_evap_layer_bare->co[l]*Dt;
          ET[l][ch] += (1. - fc) * egy->soil_evap_layer_bare[l] * Dt;
          //  ET->co[l][ch] += fc*egy->soil_evap_layer_veg->co[l]*Dt;
          ET[l][ch] += fc * egy->soil_evap_layer_veg[l] * Dt;
        }

      // water pressure and contents
      psisat =
        psi_saturation(S->thi[l][ch], pa[sy][jsat][l], pa[sy][jres][l],
                       pa[sy][ja][l], pa[sy][jns][l], 1. - 1. / pa[sy][jns][l]);
      //  th_oversat = Fmax( S->P->co[l][ch] - psisat , 0.0 ) * pa[jss][l];
      th_oversat = Fmax(S->P[l][ch] - psisat, 0.0) * pa[sy][jss][l];

      //  th->co[l][ch] =
      //Fmax(0.,egy->liq->co[l+n]+egy->deltaw->co[l+n])/(GTConst::rho_w*egy->Dlayer->co[l+n]);
      th[l][ch] = Fmax(0., egy->liq[l + n] + egy->deltaw[l + n]) /
                  (GTConst::rho_w * egy->Dlayer[l + n]);
      //  S->thi->co[l][ch] =
      //Fmax(0.,egy->ice->co[l+n]-egy->deltaw->co[l+n])/(GTConst::rho_w*egy->Dlayer->co[l+n]);
      S->thi[l][ch] = Fmax(0., egy->ice[l + n] - egy->deltaw[l + n]) /
                      (GTConst::rho_w * egy->Dlayer[l + n]);

      //  S->P->co[l][ch] = psi_teta(th->co[l][ch]+th_oversat, S->thi->co[l][ch],
      //pa[jsat][l], pa[jres][l], pa[ja][l], pa[jns][l], 1.-1./pa[jns][l],
      //GTConst::PsiMin, pa[jss][l]);
      S->P[l][ch] =
        psi_teta(th[l][ch] + th_oversat, S->thi[l][ch], pa[sy][jsat][l],
                 pa[sy][jres][l], pa[sy][ja][l], pa[sy][jns][l],
                 1. - 1. / pa[sy][jns][l], GTConst::PsiMin, pa[sy][jss][l]);

      //  temperature
      //  S->T->co[l][ch] = egy->Temp->co[l+n];
      S->T[l][ch] = egy->Temp[l + n];
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_F_energy(long nbeg,
                     long nend,
                     GeoVector<double> &F,
                     double w,
                     const GeoVector<double> &K,
                     const double *T)
{
  long l;

  for (l = nbeg; l <= nend; l++)
    {
      if (l == nbeg)
        {
          F[l] += w * (-K[l] * T[l] + K[l] * T[l + 1]);
        }
      else if (l < nend)
        {
          F[l] +=
            w * (K[l - 1] * T[l - 1] - (K[l] + K[l - 1]) * T[l] + K[l] * T[l + 1]);
        }
      else
        {
          F[l] += w * (K[l - 1] * T[l - 1] - K[l - 1] * T[l]);
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_diag_dF_energy(long nbeg,
                           long nend,
                           GeoVector<double> &dF,
                           double w,
                           const GeoVector<double> &K)
{
  long l;

  for (l = nbeg; l <= nend; l++)
    {
      if (l == nbeg)
        {
          dF[l] -= w * K[l];
        }
      else if (l < nend)
        {
          dF[l] -= w * (K[l] + K[l - 1]);
        }
      else
        {
          dF[l] -= w * K[l - 1];
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double calc_C(long l,
              long nsng,
              double a,
              double *wi,
              double *wl,
              double *dw,
              double *D,
              const GeoTensor<double> &pa,
              long sy)
{
  double C;

  if (l <= nsng)    // snow
    {
      C = (GTConst::c_ice * (wi[l] - a * dw[l]) +
           GTConst::c_liq * (wl[l] + a * dw[l])) /
          D[l];
    }
  else      // soil
    {
      C = pa(sy, jct, l - nsng) * (1. - pa(sy, jsat, l - nsng)) +
          (GTConst::c_ice * (wi[l] - a * dw[l]) +
           GTConst::c_liq * (wl[l] + a * dw[l])) /
          D[l];
    }

  return (C);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void EnergyFluxes(double t,
                  double Tg,
                  long r,
                  long c,
                  long n,
                  double Tg0,
                  double Qg0,
                  double Tv0,
                  double zmu,
                  double zmT,
                  double z0s,
                  double d0s,
                  double rz0s,
                  double z0v,
                  double d0v,
                  double rz0v,
                  double hveg,
                  double v,
                  double Ta,
                  double Qa,
                  double P,
                  double LR,
                  double psi,
                  double e,
                  double fc,
                  double LSAI,
                  double decaycoeff0,
                  double Wcrn,
                  double Wcrnmax,
                  double Wcsn,
                  double Wcsnmax,
                  double *dWcrn,
                  double *dWcsn,
                  double *theta,
                  GeoTensor<double> &soil,
                  long sy,
                  const GeoMatrix<double> &land,
                  long lu,
                  const GeoMatrix<double> &root,
                  Par *par,
                  GeoVector<double> &soil_transp_layer,
                  double SWin,
                  double LWin,
                  double SWv,
                  double *LW,
                  double *H,
                  double *dH_dT,
                  double *E,
                  double *dE_dT,
                  double *LWv,
                  double *Hv,
                  double *LEv,
                  double *Etrans,
                  double *Tv,
                  double *Qv,
                  double *Ts,
                  double *Qs,
                  double *Hg0,
                  double *Hg1,
                  double *Eg0,
                  double *Eg1,
                  double *Lobukhov,
                  double *rh,
                  double *rv,
                  double *rc,
                  double *rb,
                  double *ruc,
                  double *rh_g,
                  double *rv_g,
                  double *Qg,
                  double *u_top,
                  double *decay,
                  double *Locc,
                  double *LWup_above_v,
                  double *T,
                  GeoVector<double> &soil_evap_layer_bare,
                  GeoVector<double> &soil_evap_layer_veg,
                  double point_elev,
                  double point_slope,
                  double point_aspect,
                  double point_sky,
                  int point_lc,
                  long point_sy,
                  double snowD)
{
  double Hg, dHg_dT, Eg, dEg_dT, LWg, LWup;
  double dQgdT, rm;
  double alpha, beta;  // from Ye and Pielke, 1993
  FILE *f;

  *Tv = Tg;

  *H = 0.0;
  *E = 0.0;

  *dH_dT = 0.0;
  *dE_dT = 0.0;

  *LW = 0.0;
  *LWup_above_v = 0.0;

  *LWv = 0.0;
  *Hv = 0.0;
  *LEv = 0.0;

  *Etrans = 0.0;
  *Hg0 = 0.0;
  *Eg0 = 0.0;
  *Hg1 = 0.0;
  *Eg1 = 0.0;

  *decay = 0.0;
  *Locc = 0.0;

  // thermodynamical calculations
  SpecHumidity_2(Qg, &dQgdT, 1.0, Tg, P);

  // initialization
  if (fc == 0)
    {
      *Qs = geotop::input::gDoubleNoValue;
      *Ts = geotop::input::gDoubleNoValue;
      *Qv = geotop::input::gDoubleNoValue;
      *u_top = geotop::input::gDoubleNoValue;
      *rh = 1.E99;
      *rv = 1.E99;
      *rc = 1.E99;
      *rb = 1.E99;
      *ruc = 1.E99;

    }
  else if (fc == 1)
    {
      *rh_g = 1.E99;
      *rv_g = 1.E99;
    }

  // bare soil
  if (fc < 1)
    {
      Turbulence::aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P,
                                  LR, Lobukhov, &rm, rh_g, rv_g, par->state_turb,
                                  par->monin_obukhov, par->maxiter_Businger);

      Turbulence::find_actual_evaporation_parameters(
        r, c, &alpha, &beta, soil_evap_layer_bare, theta, soil, sy, T, psi, P,
        *rv_g, Ta, Qa, *Qg, n);
      Turbulence::turbulent_fluxes(*rh_g, *rv_g / beta, P, Ta, Tg, Qa,
                                   *Qg * alpha, dQgdT * alpha, Hg, dHg_dT, Eg,
                                   dEg_dT);

      *H += (1.0 - fc) * Hg;
      *E += (1.0 - fc) * Eg;

      *dH_dT += (1.0 - fc) * dHg_dT;
      *dE_dT += (1.0 - fc) * dEg_dT;

      // modified LP 230915
      *LW += (1.0 - fc) * (LWin - e * point_sky * SB(Tg));
      *LWup_above_v += (1.0 - fc) * ((1.0 - e) * LWin + e * SB(Tg));

      *Hg0 = Hg;
      *Eg0 = Eg;

      if (Hg != Hg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t);
          fprintf(f,
                  "Hg no value bare soil EnergyFluxes r=%ld c=%ld elev=%f slope=%f "
                  "aspect=%f sky=%f LC=%d soil type=%ld snowD=%f rh_g:%e rv_g:%e "
                  "Hg:%f Eg:%f\n",
                  r, c, point_elev, point_slope, point_aspect, point_sky,
                  (int)point_lc, point_sy, snowD, *rh_g, *rv_g, Hg, Eg);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);

          fclose(f);
          t_error("Fatal Error! GEOtop is closed. See failing report.");
        }

      if (Eg != Eg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
#ifdef USE_DOUBLE_PRECISION_OUTPUT
          fprintf(f, "Number of days after start:%12g\n", t);
          fprintf(f,
                  "Eg no value bare soil r=%ld c=%ld elev=%12g slope=%12g "
                  "aspect=%12g sky=%12g LC=%d soil type=%ld snowD=%12g\n",
                  r, c, point_elev, point_slope, point_aspect, point_sky,
                  (int)point_lc, point_sy, snowD);
          fprintf(f, "Ta:%12g Tg:%12g\n", Ta, Tg);
          // zmu: height of sensor of wind speed
          // zmt: height of sensor of AirT
          // v: wind velocity
          // d0s: displacement height
          // Tg: temperature of groudn surface
          fprintf(f,
                  "zmu:%12g zmT:%12g z0s:%12g d0s:%12g rz0s:%12g v:%12g Ta:%12g "
                  "Tg:%12g Qa:%12g Qg:%e P:%12g rm:%12g psi:%12g \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);
#else
          fprintf(f, "Number of days after start:%f\n", t);
          fprintf(f,
                  "Eg no value bare soil r=%ld c=%ld elev=%f slope=%f aspect=%f "
                  "sky=%f LC=%d soil type=%ld snowD=%f\n",
                  r, c, point_elev, point_slope, point_aspect, point_sky,
                  (int)point_lc, point_sy, snowD);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          // zmu: height of sensor of wind speed
          // zmt: height of sensor of AirT
          // v: wind velocity
          // d0s: displacement height
          // Tg: temperature of groudn surface
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);
#endif
          fclose(f);
          t_error("Fatal Error! GEOtop is closed. See failing report.");
        }
    }

  // vegetation
  if (fc > 0)
    {
      Tcanopy(r, c, Tv0, Tg, *Qg, dQgdT, Tg0, Qg0, Ta, Qa, zmu, zmT, z0v, z0s,
              d0v, rz0v, hveg, v, LR, P, SWin, SWv, LWin, e, LSAI, decaycoeff0,
              land, lu, Wcrn, Wcrnmax, Wcsn, Wcsnmax, dWcrn, dWcsn, LWv, &LWg, Hv,
              &Hg, &dHg_dT, LEv, &Eg, &dEg_dT, Ts, Qs, root, theta,
              soil_transp_layer, Lobukhov, par, n, &rm, rh, rv, rc, rb, ruc,
              u_top, Etrans, Tv, Qv, decay, Locc, &LWup, psi, soil, sy, T,
              soil_evap_layer_veg);

      *H += fc * Hg;
      *E += fc * Eg;

      *dH_dT += fc * dHg_dT;
      *dE_dT += fc * dEg_dT;

      *LW += fc * LWg;
      *LWup_above_v += fc * LWup;

      *Hg1 = Hg;
      *Eg1 = Eg;

      if (Hg != Hg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t / 86400.);
          fprintf(f,
                  "Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e "
                  "Hg:%f Eg:%f\n",
                  r, c, *rh_g, *rv_g, Hg, Eg);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);

          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (Eg != Eg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t / 86400.);
          fprintf(f, "Eg no value vegetated soil %ld %ld \n", r, c);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

// void EnergyFluxes_no_rec_turbulence(double t, double Tg, long r, long c, long
// n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s,
//                  double d0s, double rz0s, double z0v, double d0v, double
//rz0v, double hveg, double v, double Ta, double Qa,
//                  double P, double LR, double psi, double e, double fc, double
//LSAI, double decaycoeff0, double Wcrn,
//                  double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn,
//double *dWcsn, double *theta, double **soil,
//                  double *land, double *root, PAR *par, DOUBLEVECTOR
//*soil_transp_layer, double SWin, double LWin, double SWv, double *LW,
//                  double *H, double *dH_dT, double *E, double *dE_dT, double
//*LWv, double *Hv, double *LEv, double *Etrans,
//                  double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0,
//double *Hg1, double *Eg0, double *Eg1, double *Lobukhov,
//                  double *rh, double *rv, double *rc, double *rb, double *ruc,
//double *rh_g, double *rv_g, double *Qg,
//                  double *u_top, double *decay, double *Locc, double
//*LWup_above_v, double *T, DOUBLEVECTOR *soil_evap_layer_bare,
//                  DOUBLEVECTOR *soil_evap_layer_veg, short flagTmin, long
//cont){

void EnergyFluxes_no_rec_turbulence(double t,
                                    double Tg,
                                    long r,
                                    long c,
                                    long n,
                                    double /*Tg0*/,
                                    double /*Qg0*/,
                                    double /*Tv0*/,
                                    double zmu,
                                    double zmT,
                                    double z0s,
                                    double d0s,
                                    double rz0s,
                                    double /*z0v*/,
                                    double /*d0v*/,
                                    double /*rz0v*/,
                                    double /*hveg*/,
                                    double v,
                                    double Ta,
                                    double Qa,
                                    double P,
                                    double LR,
                                    double psi,
                                    double e,
                                    double fc,
                                    double LSAI,
                                    double /*decaycoeff0*/,
                                    double /*Wcrn*/,
                                    double /*Wcrnmax*/,
                                    double /*Wcsn*/,
                                    double /*Wcsnmax*/,
                                    double * /*dWcrn*/,
                                    double * /*dWcsn*/,
                                    double *theta,
                                    GeoTensor<double> &soil,
                                    long sy,
                                    const GeoMatrix<double> & /*land*/,
                                    long /*lu*/,
                                    const GeoMatrix<double> & /*root*/,
                                    Par *par,
                                    GeoVector<double> & /*soil_transp_layer*/,
                                    double /*SWin*/,
                                    double LWin,
                                    double /*SWv*/,
                                    double *LW,
                                    double *H,
                                    double *dH_dT,
                                    double *E,
                                    double *dE_dT,
                                    double *LWv,
                                    double * /*Hv*/,
                                    double * /*LEv*/,
                                    double * /*Etrans*/,
                                    double *Tv,
                                    double *Qv,
                                    double *Ts,
                                    double *Qs,
                                    double *Hg0,
                                    double *Hg1,
                                    double *Eg0,
                                    double *Eg1,
                                    double *Lobukhov,
                                    double *rh,
                                    double *rv,
                                    double *rc,
                                    double *rb,
                                    double *ruc,
                                    double *rh_g,
                                    double *rv_g,
                                    double *Qg,
                                    double * /*u_top*/,
                                    double * /*decay*/,
                                    double * /*Locc*/,
                                    double *LWup_above_v,
                                    double *T,
                                    GeoVector<double> &soil_evap_layer_bare,
                                    GeoVector<double> &soil_evap_layer_veg,
                                    double point_sky,
                                    short flagTmin,
                                    long cont)
{
  double Hg, dHg_dT, Eg, dEg_dT, LWg, LWup;
  double dQgdT;
  double dLWvdT, LWv_old;
  double alpha, beta;  // from Ye and Pielke, 1993
  double rm = 1.E20;
  FILE *f;

  // initalization
  *H = 0.0;
  *E = 0.0;

  *dH_dT = 0.0;
  *dE_dT = 0.0;

  *LW = 0.0;
  *LWup_above_v = 0.0;

  // thermodynamical calculations
  SpecHumidity_2(Qg, &dQgdT, 1.0, Tg, P);

  if (fc < 1)
    {
      if (flagTmin > 0 && cont > num_iter_after_which_only_neutrality)
        Turbulence::aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg,
                                    P, LR, Lobukhov, &rm, rh_g, rv_g,
                                    par->state_turb, 2, par->maxiter_Businger);

      Turbulence::find_actual_evaporation_parameters(
        r, c, &alpha, &beta, soil_evap_layer_bare, theta, soil, sy, T, psi, P,
        *rv_g, Ta, Qa, *Qg, n);
      Turbulence::turbulent_fluxes(*rh_g, *rv_g / beta, P, Ta, Tg, Qa,
                                   *Qg * alpha, dQgdT * alpha, Hg, dHg_dT, Eg,
                                   dEg_dT);

      *H += (1.0 - fc) * Hg;
      *E += (1.0 - fc) * Eg;

      *dH_dT += (1.0 - fc) * dHg_dT;
      *dE_dT += (1.0 - fc) * dEg_dT;

      if (par->surroundings == 1)
        {
          /* LWin corrected with temperature of surrounding terrain,
           * calculated as averaged. Also cfr. PointEnergyBalance() function. */
          *LW += (1.0 - fc) * (e * (LWin - SB(Tg)));

          /* FIXME: This should be corrected as well */
          *LWup_above_v += (1.0 - fc) * ((1.0 - e) * LWin + e * SB(Tg));
        }
      else
        {
          *LW += (1.0 - fc) * (e * (LWin - point_sky * SB(Tg)));
          *LWup_above_v += (1.0 - fc) * ((1.0 - e) * LWin + e * SB(Tg));
        }

      *Hg0 = Hg;
      *Eg0 = Eg;

      // error messages
      if (Hg != Hg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t / 86400.);
          fprintf(f,
                  "Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e "
                  "Hg:%f Eg:%f\n",
                  r, c, *rh_g, *rv_g, Hg, Eg);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (Eg != Eg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t / 86400.);
          fprintf(f, "Eg no value bare soil %ld %ld \n", r, c);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }
    }

  // TURBULENT FLUXES FOR CANOPY CANOPY
  if (fc > 0)
    {
      Turbulence::find_actual_evaporation_parameters(
        r, c, &alpha, &beta, soil_evap_layer_veg, theta, soil, sy, T, psi, P,
        *ruc, Ta, Qa, *Qg, n);
      *Ts = (Ta / (*rh) + Tg / (*ruc) + (*Tv) / (*rb)) /
            (1. / (*rh) + 1. / (*ruc) + 1. / (*rb));
      *Qs = (Qa / (*rv) + (*Qg) * alpha * beta / (*ruc) + (*Qv) / (*rc)) /
            (1. / (*rv) + beta / (*ruc) + 1. / (*rc));
      Turbulence::turbulent_fluxes(*ruc, *ruc / beta, P, *Ts, Tg, *Qs,
                                   *Qg * alpha, dQgdT * alpha, Hg, dHg_dT, Eg,
                                   dEg_dT);

      LWv_old = *LWv;
      longwave_vegetation(LWin, e, Tg, *Tv, LSAI, LWv, &LWg, &dLWvdT, &LWup);
      *LWv = LWv_old;

      *H += fc * Hg;
      *dH_dT += fc * dHg_dT;
      *E += fc * Eg;
      *dE_dT += fc * dEg_dT;

      *LW += fc * LWg;
      *LWup_above_v += fc * LWup;

      *Hg1 = Hg;
      *Eg1 = Eg;

      if (Hg != Hg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t / 86400.);
          fprintf(f,
                  "Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e "
                  "Hg:%f Eg:%f\n",
                  r, c, *rh_g, *rv_g, Hg, Eg);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);

          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (Eg != Eg)
        {
          f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
          fprintf(f, "Simulation Period:%ld\n", geotop::common::Variables::i_sim);
          fprintf(f, "Number of days after start:%f\n", t / 86400.);
          fprintf(f, "Eg no value vegetated soil %ld %ld \n", r, c);
          fprintf(f, "Ta:%f Tg:%f\n", Ta, Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f "
                  "Qg:%e P:%f rm:%f psi:%f \n",
                  zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, rm, psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double k_thermal(short snow,
                 short a,
                 double th_liq,
                 double th_ice,
                 double th_sat,
                 double k_solid)
{
  /*
   for soil:
   Quadratic parallel based on the analogy with electric lattice
   (Cosenza et al., 2003, European Journal of Soil Science, September 2003, 54,
   581–587)

   for snow:
   if a==1 Cosenza
   if a==2 Sturm(1997)
   if a==3 Jordan(1991)
   */

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();
  double r, k = 0.;

  if (snow == 0 || (a != 2 && a != 3))
    {
      k = pow((1. - th_sat) * sqrt(k_solid) + th_liq * sqrt(GTConst::k_liq) +
              th_ice * sqrt(GTConst::k_ice) +
              (th_sat - th_liq - th_ice) * sqrt(GTConst::k_air),
              2.);
      if (k != k)
        {
          lg->logsf(geotop::logger::TRACE,
                    "th_sat(%.12g) th_liq(%.12g) th_ice(%.12g)\n", th_sat, th_liq,
                    th_ice);
        }
    }
  else
    {
      r = th_liq * GTConst::rho_w + th_ice * GTConst::rho_i;

      if (a == 2)
        {
          if (r < 156)
            {
              k = 0.023 + 0.234 * r * 1.E-3;
            }
          else
            {
              k = 0.138 - 1.01 * r * 1.E-3 + 3.233 * r * r * 1.E-6;
            }
        }
      else if (a == 3)
        {
          k = GTConst::k_air +
              (7.75E-5 * r + 1.105E-6 * r * r) * (GTConst::k_ice - GTConst::k_air);
        }
    }

  return (k);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double flux(long i, long icol, double **met, double k, double est)
{
  double F = k * met[i - 1][icol];

  if ((long)F == geotop::input::gDoubleAbsent ||
      (long)F == geotop::input::gDoubleNoValue)
    {
      F = est;
    }

  return (F);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
// FIXME: this rourine below never called ????
// FIXME: should we remove it ? SC.21.04.2014)
// void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR
// *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e,
//          double *T, SHORTVECTOR *mf){

void check_errors(long r,
                  long c,
                  long n,
                  GeoVector<double> &adi,
                  GeoVector<double> &ad,
                  GeoVector<double> &ads,
                  GeoVector<double> &b,
                  GeoVector<double> &e,
                  double *T,
                  GeoVector<short> &mf)
{
  long l;
  short occ = 0;

  for (l = 1; l <= n; l++)
    {
      //  if(e->co[l]!=e->co[l]) occ=1;
      if (e[l] != e[l]) occ = 1;
    }

  if (occ == 1)
    {
      printf("NOvalue in Energy Balance: r:%ld c:%ld\n", r, c);
      //    printf("l:%d adi:--- ad:%f ads:%f b:%f e:%f T:%f mf:%d\n",
      //         1,ad->co[1],ads->co[1],b->co[1],e->co[1],T[1],mf->co[1]);
      printf("l:%d adi:--- ad:%f ads:%f b:%f e:%f T:%f mf:%d\n", 1, ad[1], ads[1],
             b[1], e[1], T[1], mf[1]);
      for (l = 2; l <= n - 1; l++)
        {
          //  printf("l:%ld adi:%f ad:%f ads:%f b:%f e:%f T:%f mf:%d\n",
          //       l,adi->co[l-1],ad->co[l],ads->co[l],b->co[l],e->co[l],T[l],
          //mf->co[l]);
          printf("l:%ld adi:%f ad:%f ads:%f b:%f e:%f T:%f mf:%d\n", l, adi[l - 1],
                 ad[l], ads[l], b[l], e[l], T[l], mf[l]);
        }
      //  printf("l:%ld adi:%f ad:%f ads:--- b:%f e:%f T:%f mf:%d\n",
      //       n,adi->co[n-1],ad->co[n],b->co[n],e->co[n],T[n],mf->co[n]);
      printf("l:%ld adi:%f ad:%f ads:--- b:%f e:%f T:%f mf:%d\n", n, adi[n - 1],
             ad[n], b[n], e[n], T[n], mf[n]);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double soil_red_evap(double psi, double T)
{
  double alpha;
  alpha =
    exp(1.0E-3 * psi * GTConst::GRAVITY / (461.50464 * (T + GTConst::tk)));
  return (alpha);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double red_evap(long n, double psi, double T)
{
  if (n > 0)
    {
      return (1.0);
    }
  else
    {
      return (soil_red_evap(psi, T));
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_roughness_soil(double z0,
                           double d0,
                           double z0_z0t,
                           double snowD,
                           double thres,
                           double z0snow,
                           double *z0_ris,
                           double *d0_ris,
                           double *z0_z0t_ris)
{
  if (snowD > thres)
    {
      *z0_ris = z0snow;
      *d0_ris = 0.0;
      *z0_z0t_ris = 0.0;
    }
  else
    {
      *z0_ris = z0;
      *d0_ris = d0;
      *z0_z0t_ris = z0_z0t;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void check_continuity(long n, double *dw, double *wi, double *wl)
{
  long l;

  for (l = 1; l <= n; l++)
    {
      if (wl[l] + dw[l] < 0) dw[l] = -wl[l];
      if (wi[l] - dw[l] < 0) dw[l] = wi[l];
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

// void merge(double a, DOUBLEVECTOR *ice, DOUBLEVECTOR *liq, DOUBLEVECTOR
// *Temp, DOUBLEVECTOR *D, long l, long lup, long ldw, long tot){
void merge(double a,
           GeoVector<double> &ice,
           GeoVector<double> &liq,
           GeoVector<double> &Temp,
           GeoVector<double> &D,
           long l,
           long lup,
           long ldw,
           long tot)
{
  long i, n, m;
  double h;

  if (l == lup)
    {
      m = l + 1;
      n = l;
    }
  else if (l < ldw)
    {
      //  if (ice->co[l+1] < ice->co[l-1]) {
      if (ice[l + 1] < ice[l - 1])
        {
          m = l + 1;
          n = l;
        }
      else
        {
          m = l;
          n = l - 1;
        }
    }
  else
    {
      m = l;
      n = l - 1;
    }

  //  h = internal_energy(ice->co[m], liq->co[m], Temp->co[m]) +
  //internal_energy(ice->co[n], liq->co[n], Temp->co[n]);
  h = internal_energy(ice[m], liq[m], Temp[m]) +
      internal_energy(ice[n], liq[n], Temp[n]);

  //  ice->co[n] += ice->co[m];
  ice[n] += ice[m];
  //  liq->co[n] += liq->co[m];
  liq[n] += liq[m];
  //  D->co[n] += D->co[m];
  D[n] += D[m];

  from_internal_energy(a, h, &(ice[n]), &(liq[n]), &(Temp[n]));

  for (i = m + 1; i <= tot; i++)
    {
      //  ice->co[i-1] = ice->co[i];
      ice[i - 1] = ice[i];
      //  liq->co[i-1] = liq->co[i];
      liq[i - 1] = liq[i];
      //  Temp->co[i-1] = Temp->co[i];
      Temp[i - 1] = Temp[i];
      //  D->co[i-1] = D->co[i];
      D[i - 1] = D[i];
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void sux_minus6_condition(double ic,
                          double wa,
                          double rho,
                          double D1,
                          Energy *E)
{
  long l;

  E->Dlayer[1] = D1;

  for (l = geotop::common::Variables::Nl; l >= 1; l--)
    {
#ifdef VERBOSE
      if (E->ice[l] != E->ice[l]) { printf("PRINT_DEBUG: %.12g", E->ice[l]); }
      if (E->liq[l] != E->liq[l]) { printf("PRINT_DEBUG: %.12g", E->liq[l]); }
#endif
      E->ice[l + 1] = E->ice[l];
      E->liq[l + 1] = E->liq[l];
      E->deltaw[l + 1] = E->deltaw[l];
      E->Temp[l + 1] = E->Temp[l];
      E->Dlayer[l + 1] = E->Dlayer[l];
    }

  E->ice[1] = ic;
  E->liq[1] = wa;
  E->Temp[1] = Fmin(E->Temp[2], -0.1);
  E->Dlayer[1] = ic / rho;

  if (E->deltaw[2] > ic)
    {
      E->deltaw[1] = ic;
      E->deltaw[2] -= ic;
    }
  else if (E->deltaw[2] > 0)
    {
      E->deltaw[1] = E->deltaw[2];
      E->deltaw[2] = 0.;
    }

  E->ice[2] -= ic;

#ifdef VERBOSE
  for (l = 1; l <= geotop::common::Variables::Nl + 1; l++)
    {
      printf("SUX-6 l:%ld ice:%f liq:%f dw:%f D:%f T:%f D1:%f ic:%f wa:%f \n", l,
             E->ice[l], E->liq[l], E->deltaw[l], E->Dlayer[l], E->Temp[l], D1, ic,
             wa);
    }
#endif
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
