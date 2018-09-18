
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


//Author: Stefano Endrizzi
//Contents: Energy balance (and also mass balance for snow and glacier)

#include "constants.h"
#include "struct.geotop.h"
#include "energy.balance.h"
#include "meteo.h"
#include "snow.h"
#include "pedo.funct.h"
#include "vegetation.h"
#include "radiation.h"
#include "turbulence.h"
#include "util_math.h"
#include "meteodistr.h"
#include "times.h"
#include "tables.h"
#include "blowingsnow.h"
#include "meteodata.h"

extern long number_novalue, number_absent;
extern char *string_novalue;

extern char *logfile;
extern long Nl;
extern long *opnt, nopnt, *obsn, nobsn;
extern double **odp, *odb;
extern long nopnt, *opnt;
extern long i_sim, i_run;
extern char *FailedRunFile;
extern char **files;

#define MM 1
#define ni 1.E-4
#define num_iter_after_which_surfenergy_balance_not_recalculated 50
#define Tmin_surface_below_which_surfenergy_balance_recalculated -50
#define num_iter_after_which_only_neutrality 10
#define Csnow_at_T_greater_than_0 1.E20
#define ratio_max_storage_RAIN_over_canopy_to_LSAI 0.1
#define ratio_max_storage_SNOW_over_canopy_to_LSAI 5.0
#define Tmax 100.0
#define Tmin -90.0

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short EnergyBalance(double Dt, double JD0, double JDb, double JDe,
                    SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G, STATE_VEG *V,
                    Vector<double> *snowage, ALLDATA *A, double *W)
{

  short sux;
  long i, r, c, cnt=0;
  static long line_interp;
  double Dtplot=0., Delta, E0, Et, SWup, Tgskin, SWrefl_surr_ave=0.,
         Tgskin_surr_ave=0.;
  FILE *f=nullptr;

  //calculation to be done before plotting maps
  *W = 0.;
  if (A->I->time==0) A->I->iplot=1;
  if (A->I->JD_plots->nh > 1 && A->I->iplot<=A->I->JD_plots->nh)
    {
      i=2*A->I->iplot-1;
      if ( (*A->P->init_date)(i_sim) + A->I->time/secinday+1.E-5 >=
           (*A->I->JD_plots)(i) &&
           (*A->P->init_date)(i_sim) + (A->I->time+Dt)/secinday-1.E-5 <=
           (*A->I->JD_plots)(i+1) )
        {
          Dtplot=((*A->I->JD_plots)(i+1) - (*A->I->JD_plots)(i))*secinday;
          *W = Dt / Dtplot;
          f = fopen(logfile, "a");
          fprintf(f,"Saving plot number %ld Weight:%f \n",A->I->iplot,*W);
          printf("Saving plot number %ld Weight:%f \n",A->I->iplot,*W);
          fclose(f);
        }
    }

  //vegetation
  if (A->I->time==0) line_interp=0;
  for (i=1; i<=A->P->n_landuses; i++)
    {
      if ((*A->P->vegflag)(i)==1)
        {
          time_interp_linear(JD0, JDb, JDe, A->L->vegparv[i-1], A->L->vegpars[i-1],
                             A->L->NumlinesVegTimeDepData[i-1], jdvegprop+1, 0, 0, &line_interp);
        }
    }

  //SUN and SHADOWING
  sun( (JDb+JDe)/2., &E0, &Et, &Delta );
  A->E->sun[1] = Delta;
  A->E->sun[8] = A->P->Lozone;
  A->E->sun[9] = A->P->alpha_iqbal;
  A->E->sun[10] = A->P->beta_iqbal;
  A->E->sun[11] = 0.;

  if (A->P->point_sim != 1)
    {
      A->E->sun[0] = A->P->latitude*Pi/180.;
      A->E->sun[2] = (A->P->longitude*Pi/180. - A->P->ST*Pi/12. + Et)/omega;
      A->E->hsun = adaptiveSimpsons2(SolarHeight__, A->E->sun, JDb, JDe, 1.E-6,
                                     20) / (JDe - JDb);
      A->E->dsun = adaptiveSimpsons2(SolarAzimuth__, A->E->sun, JDb, JDe, 1.E-6,
                                     20) / (JDe - JDb) + A->P->dem_rotation*Pi/180.;
      if (A->E->dsun < 0) A->E->dsun += 2*Pi;
      if (A->E->dsun > 2*Pi) A->E->dsun -= 2*Pi;
      A->E->sinhsun = adaptiveSimpsons2(Sinalpha_, A->E->sun, JDb, JDe, 1.E-6,
                                        20) / (JDe - JDb);
      if (A->P->cast_shadow==1) 
          shadow_haiden(A->T->Z0.get(), A->E->hsun, A->E->dsun, A->L->shadow.get());
    }

  //INITIALIZE BASIN AVERAGES
  if ((*A->P->Dtplot_basin)(i_sim)>0)
    {
      for (i=0; i<nobsn; i++)
        {
          odb[i] = 0.;
        }
    }

  //CLOUDINESS
  find_actual_cloudiness(&(A->M->tau_cloud), &(A->M->tau_cloud_av),
                         &(A->M->tau_cloud_yes), &(A->M->tau_cloud_av_yes), A->M.get(), JDb, JDe, Delta, E0,
                         Et, A->P->ST, 0., A->P->Lozone, A->P->alpha_iqbal, A->P->beta_iqbal, 0.);

  //POINT ENERGY BALANCE
  for (i=1; i<=A->P->total_channel+A->P->total_pixel; i++)
    {

      //channel or land?
      if (i<=A->P->total_channel)  //CHANNEL
        {
          r = (*A->C->r)(i);
          c = (*A->C->c)(i);
        }
      else   //LAND
        {
          r = (*A->T->rc_cont)(i - A->P->total_channel,1);
          c = (*A->T->rc_cont)(i - A->P->total_channel,2);
        }

      if ((*A->L->delay)(r,c) <= A->I->time/secinday)
        {

          cnt++;

          sux = PointEnergyBalance(i, r, c, Dt, JDb, JDe, L, C, S, G, V, snowage, A, E0,
                                   Et, Dtplot, *W, f, &SWup, &Tgskin);

          if (sux==1) return 1;

        }
      else
        {

          SWup = 0.;
          Tgskin = 0.;

        }

      //surroundings
      if (i>A->P->total_channel)
        {
          if (A->P->surroundings == 1)
            {
              SWrefl_surr_ave += SWup/(double)A->P->total_pixel;
              Tgskin_surr_ave += Tgskin/(double)A->P->total_pixel;
            }
          else
            {
                (*A->E->SWrefl_surr)(r,c) = SWup;
                (*A->E->Tgskin_surr)(r,c) = Tgskin;
            }
        }

    }

  if (A->P->surroundings == 1)
    {
      for (i=1; i<=A->P->total_pixel; i++)
        {
          r = (*A->T->rc_cont)(i,1);
          c = (*A->T->rc_cont)(i,2);
            (*A->E->SWrefl_surr)(r,c) = SWrefl_surr_ave;
            (*A->E->Tgskin_surr)(r,c) = Tgskin_surr_ave;
        }
    }

  return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short PointEnergyBalance(long i, long r, long c, double Dt, double JDb,
                         double JDe, SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G,
                         STATE_VEG *V,
                         Vector<double> *snowage, ALLDATA *A, double E0, double Et, double Dtplot,
                         double W, FILE *f, double *SWupabove_v, double *Tgskin)
{

  long l, j, ns, ng=0;
  double SWin, SW, SWbeam, SWdiff, SWv_vis, SWv_nir, SWg_vis, SWg_nir, cosinc,
         avis_b, avis_d, anir_b, anir_d;
  double SWupabove_v_vis, SWupabove_v_nir;
  double avis_ground=0., anir_ground=0.;
  short SWb_yes;
  double tauatm_sinhsun;
  double LWin, LW, LWv, epsa, eps, epsa_min, epsa_max, LWupabove_v;
  double H, LE, E, Hv, LEv, Etrans, Hg0, Eg0, Hg1, Eg1, Hadv=0.0;
  double surfEB, GEF, dUsl;
  double snowD, RainOnSnow=0., fsnow, fsnownew, fsnowcan, Melt_snow=0.,
                Evap_snow=0., Melt_glac=0., Evap_glac=0., Evap_soil=0.;
  double k_snowred, maxSWE, k;
  double Prain_over, Psnow_over, Prain=0.0, Psnow=0.0, drip_rain, drip_snow,
                                 max_wcan_rain = 0., max_wcan_snow = 0., theta_sup;
  double z0, z0veg = 0., d0, d0veg = 0., z0_z0t, hveg = 0.;
  double rh, rv, rc, rb, ruc, Lobukhov;
  double fc, fc0, decaycoeff, Locc, u_top, Qv;
  double ea, Tdew, Qa, RHpoint, Vpoint, Ppoint, Tpoint, Precpoint, zmeas_T,
         zmeas_u, Tdirichlet, Tdirichlet_bottom;
  double Ts, Qs, Qg;
  long sy;
  short lu;
  short sux, surface; //(1=yes, 0=no)
  double ic=0., wa, rho=0.;
  long lpb;

  //initialization of cumulated water volumes and set soil ancillary state vars

  if (i <= A->P->total_channel)
    {
      j = i;
      lu = (short)(*A->L->LC)(r,c);
      sy = (*A->C->soil_type)(j);

      for (l=1; l<=Nl; l++)
        {
          (*A->C->ET)(l,j) = 0.0;
          (*A->C->th)(l,j) = theta_from_psi((*C->P)(l,j), (*C->thi)(l,j), l, A->S->pa->matrix(sy), PsiMin);
          (*A->C->th)(l,j) = Fmin( (*A->C->th)(l,j), (*A->S->pa)(sy,jsat,l)-(*C->thi)(l,j) );
        }

    }
  else
    {
      j = i - A->P->total_channel;
      (*A->W->Pnet)(r,c) = 0.0;
      lu = (short)(*A->L->LC)(r,c);
      sy = (*A->S->type)(r,c);

      for (l=1; l<=Nl; l++)
        {
          (*A->S->ET)(l,r,c) = 0.0;
          (*A->S->th)(l,j) = theta_from_psi((*L->P)(l,j), (*L->thi)(l,j), l, A->S->pa->matrix(sy), PsiMin);
          (*A->S->th)(l,j) = Fmin( (*A->S->th)(l,j),(*A->S->pa)(sy,jsat,l)-(*L->thi)(l,j) );
        }
    }


  //METEO
  Tpoint=(*A->M->Tgrid)(r,c);
  Ppoint=(*A->M->Pgrid)(r,c);
  RHpoint=(*A->M->RHgrid)(r,c);
  Vpoint=(*A->M->Vgrid)(r,c);

  Precpoint=(*A->W->PrecTot)(r,c);
  //define prec as normal (not vertical)
  Precpoint*=cos((*A->T->slope)(r,c)*Pi/180.);
  //another cosine correction applied for point simulations (due to area proejection)
  if (A->P->point_sim==1 && A->P->flag1D==0)
      Precpoint*=cos((*A->T->slope)(r,c)*Pi/180.);

  Tdirichlet=A->M->var[A->M->nstTs-1][iTs];
  if ((long)Tdirichlet == number_novalue
      || (long)Tdirichlet == number_absent) Tdirichlet=A->P->Tsup;

  Tdirichlet_bottom=A->M->var[A->M->nstTbottom-1][iTbottom];
  if ((long)Tdirichlet_bottom == number_novalue
      || (long)Tdirichlet_bottom == number_absent) Tdirichlet_bottom=A->P->Tbottom;

  //SNOW
  snowD=DEPTH(r, c, S->lnum.get(), S->Dzl.get());
  ns=(*S->lnum)(r,c);

  //vegetation parameters
  if (i>A->P->total_channel)
    {

      for (l=1; l<=jdvegprop; l++)
        {
          if ( (long)A->L->vegparv[lu-1][l] != number_novalue )
            {
              A->L->vegpar->co[l] = A->L->vegparv[lu-1][l];
              if (l==jdroot)
                root(A->L->root_fraction->nch, A->L->vegpar->co[jdroot], 0.0, A->S->pa->row(1,jdz), A->L->root_fraction->row(lu));
            }
          else
            {
              A->L->vegpar->co[l] = (*A->L->ty)(lu,l+jHveg-1);
            }
        }

      if (snowD > A->L->vegpar->co[jdz0thresveg])
        {
          fsnow=1.0;
        }
      else if (snowD > A->L->vegpar->co[jdz0thresveg2])
        {
          fsnow=(snowD-A->L->vegpar->co[jdz0thresveg2])/(A->L->vegpar->co[jdz0thresveg]
                                                         -A->L->vegpar->co[jdz0thresveg2]);
        }
      else
        {
          fsnow=0.0;
        }

      fc = A->L->vegpar->co[jdcf] * pow(1.0-fsnow,A->L->vegpar->co[jdexpveg]);

      //control if LSAI is too low (in order to prevent numerical problems)
      if (A->L->vegpar->co[jdLSAI]<LSAIthres) fc=0.0;

      //glacier
      if (A->P->max_glac_layers>0)
        {
          ng=(*G->lnum)(r,c);
          if (ng>0)
            {
              A->L->vegpar->co[jdLSAI]=0.0;
              fc=0.0;
            }
        }

    }
  else
    {

      fc = 0.0;

    }

  //METEOROLOGICAL COMPUTATIONS

  //temperature and wind velocity measurement heights
  zmeas_u = (*A->M->st->Vheight)(1);
  zmeas_T = (*A->M->st->Theight)(1);

  //RAIN AND SNOW PRECIPITATION [mm]
  //convert total precipitation to [mm]
  Precpoint*=(Dt/3600.0); //from [mm/h] to [mm]

  ea = RHpoint * SatVapPressure(Tpoint, Ppoint);//Vapour Pressure [mbar]
  Qa = SpecHumidity(ea, Ppoint);//Specific Humidity
  Tdew = TfromSatVapPressure(ea, Ppoint);//Dew Temperature

  //distinguish between rain and snow
  if (A->P->dew==1)
    {
      //on the base of the dew temperature of the air
      part_snow(Precpoint, &Prain_over, &Psnow_over, Tdew, A->P->T_rain,
                A->P->T_snow);
    }
  else
    {
      //on the base of the temperature of the air
      part_snow(Precpoint, &Prain_over, &Psnow_over, Tpoint, A->P->T_rain,
                A->P->T_snow);
    }

  //Adjusting snow precipitation in case of steep slope (contribution by Stephan Gruber)
  if (A->P->snow_curv > 0 && (*A->T->slope)(r,c) > A->P->snow_smin)
    {
      if ((*A->T->slope)(r,c) <= A->P->snow_smax)
        {
          k_snowred = ( exp(-pow((*A->T->slope)(r,c) - A->P->snow_smin,
                                 2.)/A->P->snow_curv) -
                        exp(-pow(A->P->snow_smax, 2.)/A->P->snow_curv) );
        }
      else
        {
          k_snowred = 0.0;
        }
      Psnow_over *= k_snowred;
    }

  //precipitation at the surface in the unvegetated fraction
  Prain=(1.-fc)*Prain_over;
  Psnow=(1.-fc)*Psnow_over;

  if (fc>0)
    {
      canopy_rain_interception(ratio_max_storage_RAIN_over_canopy_to_LSAI,
                               A->L->vegpar->co[jdLSAI], Prain_over, &max_wcan_rain,
                               &(V->wrain->co[j]), &drip_rain);
      canopy_snow_interception(ratio_max_storage_SNOW_over_canopy_to_LSAI,
                               A->L->vegpar->co[jdLSAI], Psnow_over, V->Tv->co[j],
                               Vpoint, Dt, &max_wcan_snow, &(V->wsnow->co[j]), &drip_snow);

      //precipitation at the surface in the vegetated fraction
      Prain+=fc*drip_rain;
      Psnow+=fc*drip_snow;

      if (V->wsnow->co[j]<0) printf("Error 1 wcansnow:%f %ld %ld\n",V->wsnow->co[j],
                                      r,c);
      if (V->wsnow->co[j]<0) V->wsnow->co[j]=0.0;
    }

  //SHORTWAVE RADIATION
  //initialization of shortwave radiation absorbed by soil
  SW=0.0;
  *SWupabove_v=0.0;

  //if averaged cloud transmissivity is not available it is estimated through this micromet subroutine (not very reliable)
  if (A->M->tau_cloud_av_yes==0) A->M->tau_cloud_av = 1. -
                                                        0.71*find_cloudfactor(Tpoint, RHpoint, (*A->T->Z0)(r,c), A->M->LRv[ilsTa],
                                                            A->M->LRv[ilsTdew]);//Kimball(1928)

  //in case of shortwave data not available
  if (A->M->tau_cloud_yes==0) A->M->tau_cloud = A->M->tau_cloud_av;

  //albedo
  if (snowD>0)
    {
      if (i>A->P->total_channel) update_snow_age(Psnow_over, (*S->T)(ns,r,c), Dt,
                                                   A->P->minP_torestore_A, &(snowage->co[j]));
      avis_d=snow_albedo(avis_ground, snowD, A->P->aep, A->P->avo,
                         A->P->snow_aging_vis, snowage->co[j], 0., (*Zero));
      avis_b=avis_d;//approx
      anir_d=snow_albedo(anir_ground, snowD, A->P->aep, A->P->airo,
                         A->P->snow_aging_nir, snowage->co[j], 0., (*Zero));
      anir_b=anir_d;//approx
    }
  else
    {
      if (i<=A->P->total_channel)
        {
          theta_sup = (*A->C->th)(1,j);
        }
      else
        {
          theta_sup = (*A->S->th)(1,j);
        }
      avis_ground = find_albedo((*A->L->ty)(lu,ja_vis_dry),
                                (*A->L->ty)(lu,ja_vis_sat), theta_sup, (*A->S->pa)(sy,jres,1),
                                (*A->S->pa)(sy,jsat,1));
      anir_ground = find_albedo((*A->L->ty)(lu,ja_nir_dry),
                                (*A->L->ty)(lu,ja_nir_sat), theta_sup, (*A->S->pa)(sy,jres,1),
                                (*A->S->pa)(sy,jsat,1));
      avis_b=avis_ground;
      avis_d=avis_ground;
      anir_b=anir_ground;
      anir_d=anir_ground;
    }

  //calculation of SWin
  if (A->P->point_sim==1)
    {
      A->E->sun[0] = (*A->T->latitude)(r,c)*Pi/180.;
      A->E->sun[2] = ((*A->T->longitude)(r,c)*Pi/180. - A->P->ST*Pi/12. +
                      Et)/omega;
      A->E->hsun = adaptiveSimpsons2(SolarHeight__, A->E->sun, JDb, JDe, 1.E-6,
                                     20) / (JDe - JDb);
      A->E->dsun = adaptiveSimpsons2(SolarAzimuth__, A->E->sun, JDb, JDe, 1.E-6,
                                     20) / (JDe - JDb) + A->P->dem_rotation*Pi/180.;
      if (A->E->dsun < 0) A->E->dsun += 2*Pi;
      if (A->E->dsun > 2*Pi) A->E->dsun -= 2*Pi;
      A->E->sinhsun = adaptiveSimpsons2(Sinalpha_, A->E->sun, JDb, JDe, 1.E-6,
                                        20) / (JDe - JDb);
      if (A->P->cast_shadow==1) (*A->L->shadow)(r,c)=shadows_point(
                                                           A->T->horizon_height[(*A->T->horizon_point)(r,c)-1],
                                                           A->T->horizon_numlines[(*A->T->horizon_point)(r,c)-1], A->E->hsun*180./Pi,
                                                           A->E->dsun*180/Pi, 0., 0.);
    }

  A->E->sun[3] = RHpoint;
  A->E->sun[4] = Tpoint;
  A->E->sun[5] = Ppoint;
  A->E->sun[6] = (*A->T->slope)(r,c)*Pi/180.;
  A->E->sun[7] = (*A->T->aspect)(r,c)*Pi/180.;
  if (A->P->albedoSWin != 0) A->E->sun[11] = (avis_b + avis_d + anir_b +
                                                anir_d)/4.;

  shortwave_radiation(JDb, JDe, A->E->sun, A->E->sinhsun, E0,
                      (*A->T->sky)(r,c), (*A->E->SWrefl_surr)(r,c),
                      A->M->tau_cloud, (*A->L->shadow)(r,c), &SWbeam, &SWdiff, &cosinc,
                      &tauatm_sinhsun, &SWb_yes);

  SWbeam=flux(A->M->nstsrad, iSWb, A->M->var, 1.0, 0.0, SWbeam);
  SWdiff=flux(A->M->nstsrad, iSWd, A->M->var, 1.0,
              (1.-(*A->T->sky)(r,c)) * (*A->E->SWrefl_surr)(r,c), SWdiff);
  SWin=SWbeam+SWdiff;

  //update snow albedo
  if (snowD>0)
    {
      avis_b=snow_albedo(avis_ground, snowD, A->P->aep, A->P->avo,
                         A->P->snow_aging_vis, snowage->co[j], cosinc, (*Fzen));
      anir_b=snow_albedo(anir_ground, snowD, A->P->aep, A->P->airo,
                         A->P->snow_aging_nir, snowage->co[j], cosinc, (*Fzen));
    }

  //shortwave absorbed by soil (when vegetation is not present)
  SW += (1.0-fc)*( SWdiff*(1.0-0.5*avis_d-0.5*anir_d) + SWbeam*
                   (1.0-0.5*avis_b-0.5*anir_b) );
  *SWupabove_v += (1.0-fc)*( SWdiff*(0.5*avis_d+0.5*anir_d) + SWbeam*
                             (0.5*avis_b+0.5*anir_b) );

  //shortwave radiation absorbed by canopy
  if (fc>0)
    {
      if (A->E->hsun>0)
        {
          if (V->wsnow->co[j]<0) printf("Error wcansnow:%f %ld %ld\n",V->wsnow->co[j],r,
                                          c);
          if (V->wsnow->co[j]<0) V->wsnow->co[j]=0.0;
          fsnowcan=pow(V->wsnow->co[j]/max_wcan_snow, 2./3.);
          shortwave_vegetation(0.5*SWdiff, 0.5*SWbeam, cosinc, fsnowcan, wsn_vis,
                               Bsnd_vis, Bsnb_vis, avis_d, avis_b, (*A->L->ty)(lu,jvCh),
                               (*A->L->ty)(lu,jvR_vis), (*A->L->ty)(lu,jvT_vis),
                               A->L->vegpar->co[jdLSAI], &SWv_vis, &SWg_vis, &SWupabove_v_vis);
          shortwave_vegetation(0.5*SWdiff, 0.5*SWbeam, cosinc, fsnowcan, wsn_nir,
                               Bsnd_nir, Bsnb_nir, anir_d, anir_b, (*A->L->ty)(lu,jvCh),
                               (*A->L->ty)(lu,jvR_nir), (*A->L->ty)(lu,jvT_nir),
                               A->L->vegpar->co[jdLSAI], &SWv_nir, &SWg_nir, &SWupabove_v_nir);

        }
      else
        {
          SWv_vis=0.0;
          SWg_vis=0.0;
          SWv_nir=0.0;
          SWg_nir=0.0;
          SWupabove_v_vis = 0.0;
          SWupabove_v_nir = 0.0;
        }

      SW+=fc*(SWg_vis+SWg_nir);
      *SWupabove_v+=fc*(SWupabove_v_vis+SWupabove_v_nir);

    }
  else
    {

      SWv_vis=0.0;
      SWv_nir=0.0;
      SWupabove_v_vis = 0.0;
      SWupabove_v_nir = 0.0;

    }

  //correct in case of reading data
  SW=flux(A->M->nstsrad, iSWn, A->M->var, 1.0, 0.0, SW);

  //Extinction coefficient for SW in the snow layers
  rad_snow_absorption(r, c, A->E->SWlayer.get(), SW, S);

  //LONGWAVE RADIATION
  //soil-snow emissivity
  if (snowD>10)
    {
      eps=A->P->epsilon_snow;
    }
  else
    {
      eps=(*A->L->ty)(lu,jemg);
    }

  longwave_radiation(A->P->state_lwrad, ea, RHpoint, Tpoint, A->P->k1, A->P->k2,
                     A->M->tau_cloud_av, &epsa, &epsa_max, &epsa_min);
  LWin=(*A->T->sky)(r,c)*epsa*SB(Tpoint);

  //if incoming LW data are available, they are used (priority)
  LWin=flux(A->M->nstlrad, iLWi, A->M->var, 1.0, 0.0, LWin);

  //roughness lengths
  update_roughness_soil((*A->L->ty)(lu,jz0), 0.0, 0.0, snowD,
                        (*A->L->ty)(lu,jz0thressoil), A->P->z0_snow, &z0, &d0, &z0_z0t);
  if (fc>0) update_roughness_veg(A->L->vegpar->co[jdHveg], snowD, zmeas_u,
                                   zmeas_T, &z0veg, &d0veg, &hveg);

  //variables used in the point_surface_balance
  for (l=1; l<=Nl+ns+ng; l++)
    {
      if (l<=ns)  //snow
        {
          A->E->Dlayer->co[l] = 1.E-3*(*S->Dzl)(ns+1-l,r,c);
          A->E->liq->co[l] = (*S->w_liq)(ns+1-l,r,c);
          A->E->ice->co[l] = (*S->w_ice)(ns+1-l,r,c);
          A->E->Temp->co[l] = (*S->T)(ns+1-l,r,c);
        }
      else if (l<=ns+ng)    //glacier
        {
          A->E->Dlayer->co[l] = 1.E-3*(*G->Dzl)(ns+ng+1-l,r,c);
          A->E->liq->co[l] =(*G->w_liq)(ns+ng+1-l,r,c);
          A->E->ice->co[l] = (*G->w_ice)(ns+ng+1-l,r,c);
          A->E->Temp->co[l] =  (*G->T)(ns+ng+1-l,r,c);
        }
      else    //soil
        {
          if (i<=A->P->total_channel)
            {
              A->E->Dlayer->co[l] = 1.E-3 * (*A->S->pa)(sy,jdz,l-ns-ng);
              A->E->liq->co[l] = (*A->C->th)(l-ns-ng,j)*A->E->Dlayer->co[l]*rho_w;
              A->E->ice->co[l] = (*C->thi)(l-ns-ng,j)*A->E->Dlayer->co[l]*rho_w;
              A->E->Temp->co[l] = (*C->T)(l-ns-ng,j);
            }
          else
            {
              A->E->Dlayer->co[l] = 1.E-3*(*A->S->pa)(sy,jdz,l-ns-ng);
              A->E->liq->co[l] = (*A->S->th)(l-ns-ng,j)*A->E->Dlayer->co[l]*rho_w;
              A->E->ice->co[l] = (*L->thi)(l-ns-ng,j)*A->E->Dlayer->co[l]*rho_w;
              A->E->Temp->co[l] = (*L->T)(l-ns-ng,j);
            }
        }
    }

  //SOIL AND SNOW PROPERTIES
  sux = 0;
  lpb = 0;

  do
    {

      //init
      surface = 0;
      ic = 0.;
      wa = 0.;

      //snow or glacier melting with T[0] >= 0
      if (sux == -1)
        {
          surface = 1;

          //it is not admitted that a snow/glacier layer melts completely
          //snow layer not alone
        }
      else if (sux == -2)
        {
          merge(A->P->alpha_snow, A->E->ice.get(), A->E->liq.get(), A->E->Temp.get(), A->E->Dlayer.get(), lpb,
                1, ns, ns+ng+Nl);
          ns--;

          //glacier layer not alone
        }
      else if (sux == -3)
        {
          merge(A->P->alpha_snow, A->E->ice.get(), A->E->liq.get(), A->E->Temp.get(), A->E->Dlayer.get(), lpb,
                ns+1, ns+ng, ns+ng+Nl);
          ng --;

          //just one glacier layer and snow layers >= 0
        }
      else if (sux == -4)
        {
          merge(A->P->alpha_snow, A->E->ice.get(), A->E->liq.get(), A->E->Temp.get(), A->E->Dlayer.get(), ns+ng,
                ns, ns+ng, ns+ng+Nl);
          ng --;

          //just one snow layer and glacier layers > 0
        }
      else if (sux == -5)
        {
          merge(A->P->alpha_snow, A->E->ice.get(), A->E->liq.get(), A->E->Temp.get(), A->E->Dlayer.get(), ns,
                ns, ns+1, ns+ng+Nl);
          ns --;

          //just one snow layer and glacier layers == 0
        }
      else if (sux == -6)
        {

          ic = A->E->ice->co[1];
          wa = A->E->liq->co[1];
          rho = ic / A->E->Dlayer->co[1];
          for (l=2; l<=1+Nl; l++)
            {
              A->E->ice->co[l-1] = A->E->ice->co[l];
              A->E->liq->co[l-1] = A->E->liq->co[l];
              A->E->Temp->co[l-1] = A->E->Temp->co[l];
              A->E->Dlayer->co[l-1] = A->E->Dlayer->co[l];
            }
          ns --;
          A->E->Dlayer->co[1] *= ( 1. + ic / (A->E->ice->co[1]+A->E->liq->co[1]) );
          A->E->ice->co[1] += ic;
          if (i<=A->P->total_channel)
            {
              (*C->thi)(1,j) = A->E->ice->co[1]/(A->E->Dlayer->co[1]*rho_w);
              (*A->C->th)(1,j) = A->E->liq->co[1]/(A->E->Dlayer->co[1]*rho_w);
            }
          else
            {
              (*L->thi)(1,j) = A->E->ice->co[1]/(A->E->Dlayer->co[1]*rho_w);
              (*A->S->th)(1,j) = A->E->liq->co[1]/(A->E->Dlayer->co[1]*rho_w);
            }

        }

      //skin temperature
      A->E->Temp->co[0] = A->E->Temp->co[1];

      //ENERGY BALANCE
      sux=SolvePointEnergyBalance(surface, Tdirichlet, Tdirichlet_bottom, A->P->EB,
                                  A->P->Cair, A->P->micro, JDb-(*A->P->init_date)(i_sim), Dt, i, j, r, c, L, C,
                                  V, A->E.get(), A->L.get(),
                                  A->S.get(), A->C.get(), A->T.get(), A->P.get(), ns, ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg, d0veg,
                                  1.0, hveg, Vpoint, Tpoint, Qa, Ppoint, A->M->LRv[ilsTa],
                                  eps, fc, A->L->vegpar->co[jdLSAI], A->L->vegpar->co[jddecay0],
                                  &(V->wrain->co[j]), max_wcan_rain, &(V->wsnow->co[j]), max_wcan_snow,
                                  SWin, LWin, SWv_vis+SWv_nir, &LW, &H, &E, &LWv, &Hv, &LEv, &Etrans, &Ts, &Qs,
                                  Hadv, &Hg0, &Hg1, &Eg0, &Eg1, &Qv, &Qg, &Lobukhov,
                                  &rh, &rv, &rb, &rc, &ruc, &u_top, &decaycoeff, &Locc, &LWupabove_v, &lpb,
                                  &dUsl);


    }
  while (sux < 0);

  if (sux == 0 )
    {

      //skin temperature
      *Tgskin = A->E->Temp->co[A->P->nsurface];

      //snow melting issue
      if (ic > 0)
        {
          sux_minus6_condition(ic, wa, rho, 1.E-3 * (*A->S->pa)(sy,jdz,1), A->E.get());
          ns ++;
        }

      if (i<=A->P->total_channel)
        {

          update_soil_channel(A->P->nsurface, ns+ng, i, fc, Dt, A->E.get(), A->S->pa->matrix(sy),
                              C, A->C->ET.get(), A->C->th.get());

        }
      else
        {

          //surface energy balance
          LE = E*latent(A->E->Temp->co[1],Levap(A->E->Temp->co[1]));
          surfEB = SW + LW - H - LE;

          if (ns==0)
            {
              GEF = surfEB;
            }
          else
            {
              k = k_thermal(  0, 1,
                              Arithmetic_Mean(A->E->Dlayer->co[ns], A->E->Dlayer->co[ns+1],
                                              A->E->liq->co[ns]/(rho_w*A->E->Dlayer->co[ns]),
                                              A->E->liq->co[ns+1]/(rho_w*A->E->Dlayer->co[ns+1])),
                              Arithmetic_Mean(A->E->Dlayer->co[ns], A->E->Dlayer->co[ns+1],
                                              A->E->ice->co[ns]/(rho_w*A->E->Dlayer->co[ns]),
                                              A->E->ice->co[ns+1]/(rho_w*A->E->Dlayer->co[ns+1])),
                              Arithmetic_Mean(A->E->Dlayer->co[ns], A->E->Dlayer->co[ns+1], 1.,
                                              (*A->S->pa)(sy,jsat,1)),
                              Arithmetic_Mean(A->E->Dlayer->co[ns], A->E->Dlayer->co[ns+1], 0.,
                                              (*A->S->pa)(sy,jkt,1)) );
              GEF = k*(A->E->Temp->co[ns]-A->E->Temp->co[ns+1])/(0.5*A->E->Dlayer->co[ns]
                                                                 +0.5*A->E->Dlayer->co[ns+1]);
            }

          //assign evaporation
          if (ns>0)
            {
              Evap_snow += E*Dt;
            }
          else if (ng>0)
            {
              Evap_glac += E*Dt;
            }
          else
            {
              Evap_soil += E*Dt;
            }

          //soil
          update_soil_land(A->P->nsurface, ns+ng, j, r, c, fc, Dt, A->E.get(),
                           A->S->pa->matrix(sy), L, A->S->ET.get(), A->S->th.get());

          //glacier
          if (A->P->max_glac_layers>0)
            {
              //glacier melting
              WBglacier(ns, ng, r, c, G, &Melt_glac, A->P.get(), A->E.get(), Evap_glac);
              //adjust glacier layers
                snow_layer_combination(A->P->alpha_snow, r, c, G, Tpoint, A->P->inf_glac_layers.get(),
                                       A->P->max_weq_glac, 1.E10);
            }

          //SNOW
          //snow melting, compactation and percolation
          WBsnow(Dt, ns, r, c, S, &Melt_snow, &RainOnSnow, A->P.get(), (*A->T->slope)(r,c),
                 Prain, A->E.get(), Evap_snow);

          //adjust snow layers
          if (A->P->point_sim == 1)
            {
              maxSWE = (*A->P->maxSWE)(r,c);
            }
          else
            {
              maxSWE = 1.E10;
            }
            snow_layer_combination(A->P->alpha_snow, r, c, S, Tpoint, A->P->inf_snow_layers.get(), A->P->max_weq_snow,
                                   maxSWE);

          //add new snow
          if (Psnow>0) new_snow(A->P->alpha_snow, r, c, S, Psnow,
                                  Psnow*rho_w/rho_newlyfallensnow(Vpoint, Tpoint, Tfreezing), Tpoint);

          //NET PRECIPITATION
          (*A->W->Pnet)(r,c) += (Melt_snow + Melt_glac + Prain);

          //VEGETATION
          if ( A->L->vegpar->co[jdLSAI]>=LSAIthres && ng==0 )
            {

              snowD = DEPTH(r, c, S->lnum.get(), S->Dzl.get());

              fc0=fc;

              if (snowD > A->L->vegpar->co[jdz0thresveg])
                {
                  fsnownew=1.0;
                }
              else if (snowD > A->L->vegpar->co[jdz0thresveg2])
                {
                  fsnownew=(snowD-A->L->vegpar->co[jdz0thresveg2])/
                           (A->L->vegpar->co[jdz0thresveg]-A->L->vegpar->co[jdz0thresveg2]);
                }
              else
                {
                  fsnownew=0.0;
                }

              fc = A->L->vegpar->co[jdcf]*pow(1.0-fsnownew,A->L->vegpar->co[jdexpveg]);

              //a) fc increases
              if (fc>fc0)
                {
                  V->wrain->co[j]*=(fc0/fc);
                  V->wsnow->co[j]*=(fc0/fc);

                  //b) fc decreases
                }
              else if (fc<fc0)
                {
                  new_snow(A->P->alpha_snow, r, c, S, V->wsnow->co[j]*(fc0-fc),
                           V->wsnow->co[j]*(fc0-fc)*1000./300, V->Tv->co[j]);
                  (*A->W->Pnet)(r,c) +=  V->wrain->co[j]*(fc0-fc);
                }

              if (fc<1.E-6)
                {
                  V->wrain->co[j]=0.0;
                  V->wsnow->co[j]=0.0;
                }
            }

          //write output
          if ((*A->P->output_snow)(i_sim)>0)
            {
              if (strcmp(files[fsnowmelt],
                         string_novalue) != 0) (*A->N->melted)(j) = Melt_snow;
              if (strcmp(files[fsnowsubl],
                         string_novalue) != 0) (*A->N->subl)(j) = Evap_snow;
              if (strcmp(files[fsndur], string_novalue) != 0)
                {
                  if (snowD>0)
                    {
                      (*A->N->yes)(j) = 1;
                    }
                  else
                    {
                      (*A->N->yes)(j) = 0;
                    }
                }
            }

          if (A->P->max_glac_layers>0 && (*A->P->output_glac)(i_sim)>0)
            {
              if (strcmp(files[fglacmelt],
                         string_novalue) != 0) (*A->G->melted)(j) = Melt_glac;
              if (strcmp(files[fglacsubl],
                         string_novalue) != 0) (*A->G->subl)(j) = Evap_glac;
            }

          if (A->P->output_surfenergy->co[i_sim]>0)
            {
              if (strcmp(files[fradnet], string_novalue) != 0) A->E->Rn->co[j] = (SW+LW);
              if (strcmp(files[fradLWin], string_novalue) != 0) A->E->LWin->co[j] = LWin;
              if (strcmp(files[fradLW], string_novalue) != 0) A->E->LW->co[j] = LW;
              if (strcmp(files[fradSW], string_novalue) != 0)  A->E->SW->co[j] = SW;
              if (strcmp(files[fradSWin], string_novalue) != 0) A->E->SWin->co[j] = SWin;
              if (strcmp(files[fradSWinbeam], string_novalue) != 0) A->E->SWinb->co[j] = SWbeam;
              if (strcmp(files[fshadow], string_novalue) != 0) (*A->E->shad)(j) = SWb_yes;
              if (strcmp(files[fG], string_novalue) != 0) A->E->G->co[j] = surfEB;
              if (strcmp(files[fH], string_novalue) != 0)  A->E->H->co[j] = H;
              if (strcmp(files[fLE], string_novalue) != 0)  A->E->LE->co[j] = LE;
              if (strcmp(files[fTs], string_novalue) != 0) A->E->Ts->co[j] = (*Tgskin);
            }

          if ((*A->P->output_meteo)(i_sim)>0)
            {
              if (strcmp(files[fprec], string_novalue) != 0)
                {
                 (*A->W->Pt)(j) = Precpoint;
                 (*A->W->Ps)(j) = Psnow_over;
                }
            }

          if (A->I->JD_plots->nh > 1 && W>0)
            {
              if (strcmp(files[pH], string_novalue) != 0
                  || strcmp(files[pHg], string_novalue) != 0
                  || strcmp(files[pG], string_novalue) != 0) A->E->Hgp->co[j] = W*H;
              if (strcmp(files[pH], string_novalue) != 0
                  || strcmp(files[pHv], string_novalue) != 0) A->E->Hvp->co[j] = W*fc*Hv;
              if (strcmp(files[pLE], string_novalue) != 0
                  || strcmp(files[pLEg], string_novalue) != 0
                  || strcmp(files[pG], string_novalue) != 0) A->E->LEgp->co[j] = W*LE;
              if (strcmp(files[pLE], string_novalue) != 0
                  || strcmp(files[pLEv], string_novalue) != 0) A->E->LEvp->co[j] = W*fc*LEv;
              if (strcmp(files[pSWin], string_novalue) != 0) A->E->SWinp->co[j] = W*SWin;
              if (strcmp(files[pSWg], string_novalue) != 0
                  || strcmp(files[pG], string_novalue) != 0) A->E->SWgp->co[j] = W*SW;
              if (strcmp(files[pSWv], string_novalue) != 0) A->E->SWvp->co[j] = W*fc*
                    (SWv_vis+SWv_nir);
              if (strcmp(files[pLWin], string_novalue) != 0) A->E->LWinp->co[j] = W*LWin;
              if (strcmp(files[pLWg], string_novalue) != 0
                  || strcmp(files[pG], string_novalue) != 0) A->E->LWgp->co[j] = W*LW;
              if (strcmp(files[pLWv], string_novalue) != 0) A->E->LWvp->co[j] = W*fc*LWv;
              if (strcmp(files[pTs], string_novalue) != 0) A->E->Tsp->co[j] = W*Ts;
              if (strcmp(files[pTg], string_novalue) != 0) A->E->Tgp->co[j] = W*(*Tgskin);
            }

          if (A->P->state_pixel == 1)
            {
              if ((*A->P->jplot)(j) > 0)
                {
                  if ((*A->P->Dtplot_point)(i_sim)>0)
                    {
                      for (l=0; l<nopnt; l++)
                        {
                          if (opnt[l] > 0)
                            {
                              if (opnt[l] == osnowover)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Psnow_over;
                                }
                              else if (opnt[l] == orainover)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Prain_over;
                                }
                              else if (opnt[l] == oprecsnow)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Psnow;
                                }
                              else if (opnt[l] == oprecrain)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Prain;
                                }
                              else if (opnt[l] == orainonsnow)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = RainOnSnow;
                                }
                              else if (opnt[l] == oV)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Vpoint*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oVdir)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    (*A->M->Vdir)(r,c)*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oRH)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = RHpoint*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oPa)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    Ppoint*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oTa)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Tpoint*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oTdew)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Tdew*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oTg)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = *Tgskin*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oTv)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    V->Tv->co[j]*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oTs)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Ts*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oEB)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = surfEB*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oG)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = GEF*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oSWin)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = SWin*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oSWb)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = SWbeam*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oSWd)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = SWdiff*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLWin)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = LWin*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == ominLWin)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = ( (*A->T->sky)(r,c) * epsa_min * SB(Tpoint) + (1.-(*A->T->sky)(r,c)) * eps * SB((*A->E->Tgskin_surr)(r,c)) )
                                                                        * Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == omaxLWin)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = ((*A->T->sky)(r,c)*epsa_max*SB(Tpoint) + (1.-(*A->T->sky)(r,c)) * eps * SB((*A->E->Tgskin_surr)(r,c)) )
                                                                      * Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oSW)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = SW*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLW)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = LW*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oH)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = H*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLE)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = LE*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == ofc)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = fc*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLSAI)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    A->L->vegpar->co[jdLSAI]*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oz0v)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = z0veg*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == od0v)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = d0veg*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oEcan)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = (SWv_vis+SWv_nir+LWv-Hv-LEv)
                                                                       *Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oSWv)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = (SWv_vis+SWv_nir)
                                                                       *Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLWv)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = LWv*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oHv)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Hv*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLEv)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = LEv*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oHg0)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Hg0*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLEg0)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Levap(*Tgskin)
                                                                       *Eg0*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oHg1)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Hg1*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLEg1)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Levap(*Tgskin)
                                                                       *Eg1*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oevapsur)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Evap_soil;
                                }
                              else if (opnt[l] == otrasp)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Etrans*Dt;
                                }
                              else if (opnt[l] == oQv)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Qv*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oQg)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Qg*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oQa)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Qa*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oQs)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Qs*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLobuk)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    Lobukhov*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLobukcan)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Locc*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == outop)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = u_top*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == odecay)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    decaycoeff*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oSWup)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    *SWupabove_v*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == oLWup)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] =
                                    LWupabove_v*Dt/(*A->P->Dtplot_point)(i_sim);
                                }
                              else if (opnt[l] == omrsnow)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Melt_snow;
                                }
                              else if (opnt[l] == osrsnow)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Evap_snow;
                                }
                              else if (opnt[l] == omrglac)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Melt_glac;
                                }
                              else if (opnt[l] == osrglac)
                                {
                                  odp[opnt[l]][(*A->P->jplot)(j)-1] = Evap_glac;
                                }
                            }
                        }
                    }
                }
            }

          if (A->P->state_basin == 1)
            {
              if ((*A->P->Dtplot_basin)(i_sim)>0)
                {
                  for (l=0; l<nobsn; l++)
                    {
                      if (obsn[l] > 0)
                        {
                          if (obsn[l] == ooprecrain)
                            {
                              odb[obsn[l]] += Prain/(double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooprecsnow)
                            {
                              odb[obsn[l]] += Psnow/(double)A->P->total_pixel;
                            }
                          else if (obsn[l] == oorainover)
                            {
                              odb[obsn[l]] += Prain_over/(double)A->P->total_pixel;
                            }
                          else if (obsn[l] == oosnowover)
                            {
                              odb[obsn[l]] += Psnow_over/(double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooTa)
                            {
                              odb[obsn[l]] += Tpoint*(Dt/(*A->P->Dtplot_basin)(i_sim))/
                                              (double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooTg)
                            {
                              odb[obsn[l]] += *Tgskin*(Dt/(*A->P->Dtplot_basin)(i_sim))/
                                              (double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooTv)
                            {
                              odb[obsn[l]] += V->Tv->co[j]*(Dt/(*A->P->Dtplot_basin)(i_sim))/
                                              (double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooevapsur)
                            {
                              odb[obsn[l]] += Evap_soil/(double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ootrasp)
                            {
                              odb[obsn[l]] += (Etrans*Dt)/(double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooLE)
                            {
                              odb[obsn[l]] += LE*(Dt/(*A->P->Dtplot_basin)(i_sim))/(double)
                                              A->P->total_pixel;
                            }
                          else if (obsn[l] == ooH)
                            {
                              odb[obsn[l]] += H*(Dt/(*A->P->Dtplot_basin)(i_sim))/(double)
                                              A->P->total_pixel;
                            }
                          else if (obsn[l] == ooSW)
                            {
                              odb[obsn[l]] += SW*(Dt/(*A->P->Dtplot_basin)(i_sim))/(double)
                                              A->P->total_pixel;
                            }
                          else if (obsn[l] == ooLW)
                            {
                              odb[obsn[l]] += LW*(Dt/(*A->P->Dtplot_basin)(i_sim))/(double)
                                              A->P->total_pixel;
                            }
                          else if (obsn[l] == ooLEv)
                            {
                              odb[obsn[l]] += LEv*(Dt/(*A->P->Dtplot_basin)(i_sim))/
                                              (double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooHv)
                            {
                              odb[obsn[l]] += Hv*(Dt/(*A->P->Dtplot_basin)(i_sim))/(double)
                                              A->P->total_pixel;
                            }
                          else if (obsn[l] == ooSWv)
                            {
                              odb[obsn[l]] += (SWv_vis+SWv_nir)*(Dt/(*A->P->Dtplot_basin)(i_sim))/
                                              (double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooLWv)
                            {
                              odb[obsn[l]] += LWv*(Dt/(*A->P->Dtplot_basin)(i_sim))/
                                              (double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooSWin)
                            {
                              odb[obsn[l]] += SWin*(Dt/(*A->P->Dtplot_basin)(i_sim))/
                                              (double)A->P->total_pixel;
                            }
                          else if (obsn[l] == ooLWin)
                            {
                              odb[obsn[l]] += LWin*(Dt/(*A->P->Dtplot_basin)(i_sim))/
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

      return 1;

    }




}

//end of "PointEnergyBalance" subroutine

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short SolvePointEnergyBalance(short surfacemelting, double Tgd,
                              double Tbottom, double EBd, double Convd, short surfacebalance, double t,
                              double Dt, long i, long j, long r, long c, SOIL_STATE *SL,
                              SOIL_STATE *SC, STATE_VEG *V, ENERGY *egy, LAND *land, SOIL *sl,
                              CHANNEL *cnet, TOPO *top, PAR *par, long ns, long ng, double zmu, double zmT,
                              double z0s, double d0s,
                              double rz0s, double z0v, double d0v, double rz0v, double hveg, double v,
                              double Ta, double Qa, double P, double LR, double eps, double fc, double LSAI,
                              double decaycoeff0, double *Wcrn, double Wcrnmax, double *Wcsn,
                              double Wcsnmax, double SWin, double LWin, double SWv, double *LW, double *H,
                              double *E,
                              double *LWv, double *Hv, double *LEv, double *Etrans, double *Ts, double *Qs,
                              double Eadd, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Qv,
                              double *Qg,
                              double *Lob, double *rh, double *rv, double *rb, double *rc, double *ruc,
                              double *u_top, double *decay, double *Locc, double *LWup_ab_v, long *lpb,
                              double *dUsl)
{


  short iter_close, iter_close2, lu=(*land->LC)(r,c), flagTmin = 0, sux,
                                 dirichlet = 0, neumann = 0, micro_sempl = 0, micro = 0, dirichlet_bottom = 0;
  long sur, sy, l, m, cont=0, cont2, n, cont_lambda_min=0;
  double EB=0., dEB_dT=0., dH_dT, dE_dT, EB0, Tg, Tg0, psim0, psi0, Qg0=0.,
         Tv0=0., dWcsn=0.0, dWcrn=0.0, rh_g, rv_g;
  double res, res0[3], res_av, res_prev[MM], lambda[3], C0, C1, th0, th1, kbb0,
         kbb1, kub0=0., kub1=0., thi=0., thin, thw=0., thwn, sat=0., satn, kt=0., ktn;
  FILE *f;

  res0[0]=0.;
  res0[1]=0.;
  res0[2]=0.;
  lambda[0]=0.;
  lambda[1]=0.;
  lambda[2]=0.;

  //boundary condition
  //up
  if ((long)Tgd != number_novalue) dirichlet = 1;
  if ((long)EBd != number_novalue) neumann = 1;
  if ((long)Convd != number_novalue) micro_sempl = 1;
  if (neumann == 0 && dirichlet == 0 && micro_sempl == 0) micro = 1;
  //down
  if ((long)Tbottom != number_novalue) dirichlet_bottom = 1;

  //Soil layer
  n = Fminlong((*par->Nl_spinup)(i_sim),Nl) + ns + ng;

  //Surface conditions
  sur = par->nsurface;
  if (surfacemelting == 1) sur = 1;
  if (neumann == 1) EB = EBd;
  if (micro_sempl == 1) dEB_dT = -Convd;

  if (dirichlet == 1)
    {
      sur = 1;
      egy->Temp->co[0] = Tgd;
    }

  //Initializes vectors
  for (l=sur; l<=n; l++)
    {
      egy->T0->co[l] = egy->Temp->co[l];
    }

  //Initialize soil properties
  if (i<=par->total_channel)
    {
      sy = (*cnet->soil_type)(j);
      psi0 = (*SC->P)(0,j);
      for (l=1; l<=Nl; l++)
        {
          //water content to be modified at each iteration
          egy->THETA->co[l] = (*cnet->th)(l,j);
          //total pressure (=pressure of the water would have if it was all in liquid state)
          psim0=psi_teta((*cnet->th)(l,j)+(*SC->thi)(l,j), 0.0,
                         (*sl->pa)(sy,jsat,l), (*sl->pa)(sy,jres,l),
                         (*sl->pa)(sy,ja,l), (*sl->pa)(sy,jns,l), 1-1/(*sl->pa)(sy,jns,l),
                         PsiMin, (*sl->pa)(sy,jss,l));
          //max temperature at which the first particle of ice comes up
          egy->Tstar->co[l]=Fmin(psim0/(1000.0*Lf/(g*(Tfreezing+tk))), 0.0);
        }
    }
  else
    {
      sy = (*sl->type)(r,c);
      psi0 = (*SL->P)(0,j);
      for (l=1; l<=Nl; l++)
        {
          //water content to be modified at each iteration
          egy->THETA->co[l] = (*sl->th)(l,j);
          //total pressure (=pressure of the water would have if it was all in liquid state)
          psim0=psi_teta((*sl->th)(l,j)+(*SL->thi)(l,j), 0.0,
                         (*sl->pa)(sy,jsat,l), (*sl->pa)(sy,jres,l),
                         (*sl->pa)(sy,ja,l), (*sl->pa)(sy,jns,l), 1-1/(*sl->pa)(sy,jns,l),
                         PsiMin, (*sl->pa)(sy,jss,l));
          //max temperature at which the first particle of ice comes up
          egy->Tstar->co[l]=Fmin(psim0/(1000.0*Lf/(g*(Tfreezing+tk))), 0.0);
        }

    }

  //Calculate the surface energy balance only if you need it
  if (surfacebalance == 1)
    {

      //Surface skin temperature used to compute the surface energy fluxes
      Tg = egy->Temp->co[sur];
      Tg0 = Tg;

      //Qg0 is the specific humidity at the soil surface at the beginning of the time step (assumed saturated)
      Qg0 = SpecHumidity(SatVapPressure(Tg0, P), P);

      //canopy temperature at the beginning of the time step
      Tv0 = V->Tv->co[j];

      /*calculates all the surface energy fluxes, includes the calculation of vegetation temperature (solving the vegetation energy balance)
       Returns:

       LW: net longwave radiation at the surface
       H: sensible heat flux at the surface (positive upwards)
       E: evaporation from the soil surface (mm/s)

       LWv: net longwave in the canopy
       Hv: sensible heat flux from the canopy
       LEv: latent heat flux from the canopy
       Etrans: canopy transpiration (as water flux in mm/s)

       sl->Tv->co[r][c] : vegetation temperature
       Qv : specific humidity at the canopy (assumed saturated)
       Ts : temperature of canopy air
       Qs : specific humidity of canopy air*/

      //surface energy balance
      EnergyFluxes(t, Tg, r, c, ns+ng, Tg0, Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v,
                   d0v, rz0v, hveg, v, Ta, Qa, P, LR, psi0, eps, fc, LSAI,
                   decaycoeff0, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, &egy->THETA->co[0],
                   sl->pa->co[sy], land->ty->row(lu),
                   land->root_fraction->row(lu), par, egy->soil_transp_layer.get(), SWin, LWin, SWv, LW,
                   H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans,
                   &(V->Tv->co[j]), Qv, Ts, Qs, Hg0, Hg1, Eg0, Eg1, Lob, rh, rv, rc, rb, ruc,
                   &rh_g, &rv_g, Qg, u_top, decay, Locc, LWup_ab_v,
                   &egy->Temp->co[0], egy->soil_evap_layer_bare.get(), egy->soil_evap_layer_veg.get(),
                   (*top->sky)(r,c));

      if (micro == 1)
        {
          EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Eadd + egy->SWlayer->co[0];
          dEB_dT = -eps*(*top->sky)(r,c)*dSB_dT(Tg) - dH_dT - latent(Tg,
                   Levap(Tg))*dE_dT;//changed
        }
    }

  if (micro_sempl == 1)  EB = Convd * (Ta - egy->Temp->co[sur]);

  EB0 = EB;

  //Calculate -F(x0), given that the system to solve is F'(x0) * (x-x0) = -F(x0)
  //F(T) = diag(egy->Fenergy(T)) + K(T)*T

  for (l=sur; l<=n; l++)
    {
      egy->Fenergy->co[l] = 0.0;
      if (l>0)egy->deltaw->co[l] = 0.0;
    }

  for (l=1; l<=n; l++)
    {

      //conductive M-matrix (lower diagonal part of this M-matrix is stored in a vector)
      thw = (egy->liq->co[l]+egy->deltaw->co[l])/(rho_w*egy->Dlayer->co[l]);
      thi = (egy->ice->co[l]-egy->deltaw->co[l])/(rho_i*egy->Dlayer->co[l]);

      if (l > ns+ng)
        {
          sat = (*sl->pa)(sy,jsat,l-ns-ng);
          kt = (*sl->pa)(sy,jkt,l-ns-ng);
        }
      else
        {
          sat = 1.;
          kt = 0.;
        }

      if (l==1)
        {
          if (l<=ns+ng)
            {
              kub1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
            }
          else
            {
              kub1 = k_thermal(0, 1, thw, thi, sat, kt);
            }
          kub0 = kub1;
        }

      if (l>1)
        {

          thwn = (egy->liq->co[l-1]+egy->deltaw->co[l-1])/(rho_w*egy->Dlayer->co[l-1]);
          thin = (egy->ice->co[l-1]-egy->deltaw->co[l-1])/(rho_i*egy->Dlayer->co[l-1]);

          if (l-1 > ns+ng)
            {
              satn = (*sl->pa)(sy,jsat,l-1-ns-ng);
              ktn = (*sl->pa)(sy,jkt,l-1-ns-ng);
            }
          else
            {
              satn = 1.;
              ktn = 0.;
            }

          thwn = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], thw, thwn);
          thin = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], thi, thin);
          satn = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], sat, satn);
          ktn = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], kt, ktn);

          if (l <= ns+ng)
            {
              egy->Kth1->co[l-1] = - 2. * k_thermal(1, par->snow_conductivity, thwn, thin,
                                                    satn, ktn) / ( egy->Dlayer->co[l-1] + egy->Dlayer->co[l] );
            }
          else
            {
              egy->Kth1->co[l-1] = - 2. * k_thermal(0, 1, thwn, thin, satn,
                                                    ktn) / ( egy->Dlayer->co[l-1] + egy->Dlayer->co[l] );
            }

        }
      else if (l>sur)
        {

          if (l <= ns+ng)
            {
              egy->Kth1->co[l-1] = - k_thermal(1, par->snow_conductivity, thw, thi, sat,
                                               kt) / ( egy->Dlayer->co[l]/2. );
            }
          else
            {
              egy->Kth1->co[l-1] = - k_thermal(0, 1, thw, thi, sat,
                                               kt) / ( egy->Dlayer->co[l]/2. );
            }

        }

      egy->Kth0->co[l-1] = egy->Kth1->co[l-1];

      //diagonal part of F due to heat capacity and boundary condition (except conduction)
      C0 = calc_C(l, ns+ng, 0.0, &egy->ice->co[0], &egy->liq->co[0], &egy->deltaw->co[0],
                  &egy->Dlayer->co[0], sl->pa->co[sy]);
      C1 = calc_C(l, ns+ng, 1.0, &egy->ice->co[0], &egy->liq->co[0], &egy->deltaw->co[0],
                  &egy->Dlayer->co[0], sl->pa->co[sy]);
      if (l<=ns+ng
          && egy->ice->co[l]-egy->deltaw->co[l]<1.E-7) C1 = Csnow_at_T_greater_than_0;
      egy->Fenergy->co[l] += ( Lf*egy->deltaw->co[l] +
                               C1*egy->Dlayer->co[l]*egy->Temp->co[l] -
                               C0*egy->Dlayer->co[l]*egy->T0->co[l] ) / Dt;

      //shortwave radiation penetrating under the surface
      if (micro == 1 && l<=ns+1) egy->Fenergy->co[l] -= egy->SWlayer->co[l];

    }

  //top boundary condition
  egy->Fenergy->co[sur] -= ( (1.-KNe)*EB + KNe*EB0 );
  if (dirichlet == 1) egy->Fenergy->co[sur] -= ( (Tgd-egy->Temp->co[sur])*
                                                   (1.-KNe)*kub1 + (Tgd-egy->T0->co[sur])*KNe*kub0 ) / (egy->Dlayer->co[1]/2.);

  //bottom boundary condition (treated as sink)
  egy->Fenergy->co[n] -= par->Fboundary;

  if (n <= ns+ng)
    {
      kbb1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
    }
  else
    {
      kbb1 = k_thermal(0, 1, thw, thi, sat, kt);
    }
  kbb0 = kbb1;

  if (dirichlet_bottom == 1)
    {
      egy->Fenergy->co[n] -= ( (Tbottom-egy->Temp->co[n])*(1.-KNe)*kbb1 +
                               (Tbottom-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlayer->co[n]/2.);
    }
  else
    {
      egy->Fenergy->co[n] -= ( (par->Tboundary-egy->Temp->co[n])*(1.-KNe)*kbb1 +
                               (par->Tboundary-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlayer->co[n]/2.
                                   +par->Zboundary);
    }

  //include conduction in F(x0)
  update_F_energy(sur, n, egy->Fenergy.get(), 1.-KNe, egy->Kth1.get(), &egy->Temp->co[0]);
  update_F_energy(sur, n, egy->Fenergy.get(), KNe, egy->Kth0.get(), &egy->T0->co[0]);

  //calculate the norm
  res = norm_2(egy->Fenergy.get(), sur, n);


  do
    {

      cont++;

      for (l=sur; l<=n; l++)
        {
          egy->T1->co[l] = egy->Temp->co[l];
        }

      //Calculates F'(x) referred to as Jacobian matrix

      //First the diagonal part is calculated and put into the vector dFenergy

      for (l=sur; l<=n; l++)
        {
          egy->dFenergy->co[l] = 0.0;
        }

      //Heat capacity part
      for (l=1; l<=n; l++)
        {

          //real thermal capacity
          C1 = calc_C(l, ns+ng, 1.0, &egy->ice->co[0], &egy->liq->co[0], &egy->deltaw->co[0],
                      &egy->Dlayer->co[0], sl->pa->co[sy]);
          if (l<=ns+ng
              && egy->ice->co[l]-egy->deltaw->co[l]<1.E-7) C1 = Csnow_at_T_greater_than_0;
          //adds apparent thermal conductivity (due to phase change) for both snow and soil
          if (l<=ns+ng) //snow
            {
              C1 += Lf*(egy->ice->co[l]+egy->liq->co[l])*dtheta_snow(par->alpha_snow, 1.,
                                                                     egy->Temp->co[l])/egy->Dlayer->co[l];
            }
          else    //soil
            {
              m=l-ns-ng;
              if (egy->Temp->co[l]<=egy->Tstar->co[m]) C1 += rho_w * Lf * (Lf/(g*tk)*1.E3) 
                      * dteta_dpsi(Psif(egy->Temp->co[l]), 0.0, (*sl->pa)(sy,jsat,m),
                                  (*sl->pa)(sy,jres,m), (*sl->pa)(sy,ja,m), (*sl->pa)(sy,jns,m),
                                  1-1/(*sl->pa)(sy,jns,m), PsiMin, 0.0);
            }

          egy->dFenergy->co[l] += C1*egy->Dlayer->co[l] / Dt;

        }

      //Upper Boundary Condition: Neumann
      egy->dFenergy->co[sur] -= (1.-KNe)*dEB_dT;

      //Upper Boundary Condition: Dirichlet
      if (dirichlet == 1) egy->dFenergy->co[sur] += (1.-KNe)*kub1 /
                                                      (egy->Dlayer->co[1]/2.);

      //Bottom Boundary Condition
      egy->dFenergy->co[n] += (1.-KNe)*kbb1 / (egy->Dlayer->co[n]/2.
                                               +par->Zboundary);

      //Update the part of Jacobian due to conduction is included in Kth1 and Kth0
      update_diag_dF_energy(sur, n, egy->dFenergy.get(), 1.-KNe, egy->Kth1.get());

      //Extradiagonal part of the Jacobian
      for (l=sur; l<=n-1; l++)
        {
          egy->udFenergy->co[l] = (1.-KNe)*egy->Kth1->co[l];
        }

      //Solve tridiagonal system
      sux = tridiag2(1, r, c, sur, n, egy->udFenergy.get(), egy->dFenergy.get(), egy->udFenergy.get(),
                     egy->Fenergy.get(), egy->Newton_dir.get());

      if (sux == 1)
        {
          if (ns != 1)
            {
              return 1;
            }
          else if (ng>0)
            {
              return -5;
            }
          else
            {
              return -6;
            }
        }

      cont2=0;
      iter_close=0;
      res0[0]=res;

      //non-monotonic line search (it is monotonic if M==1)
      for (m=Fminlong(cont,MM); m>1; m--)
        {
          res_prev[m-1]=res_prev[m-2];
        }
      res_prev[0]=res;

      res_av=0.0;
      for (m=1; m<=Fminlong(cont,MM); m++)
        {
          res_av=Fmax(res_prev[m-1],res_av);
        }

      do
        {

          cont2++;
          iter_close2=0;

          if (cont2 == 1)
            {
              lambda[0] = 1.0;

            }
          else if (cont2 == 2)
            {
              lambda[1] = lambda[0];
              res0[1] = res;
              lambda[0] = thmax;

            }
          else
            {
              lambda[2] = lambda[1];
              res0[2] = res0[1];
              lambda[1] = lambda[0];
              res0[1] = res;
              lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2],
                                                  res0[2]);

            }

          for (l=sur; l<=n; l++)
            {

              egy->Temp->co[l] = egy->T1->co[l] + lambda[0] * egy->Newton_dir->co[l];


              //soil
              if (l > ns+ng)
                {
                  m=l-ns-ng;

                  if (i<=par->total_channel)
                    {
                      th0 = (*cnet->th)(m,j);
                      th1 = teta_psi(Psif(Fmin(egy->Tstar->co[m],egy->Temp->co[l])), 0.0,
                                     (*sl->pa)(sy,jsat,m),
                                     (*sl->pa)(sy,jres,m), (*sl->pa)(sy,ja,m), (*sl->pa)(sy,jns,m),
                                     1-1/(*sl->pa)(sy,jns,m), PsiMin, (*sl->pa)(sy,jss,m));
                      if (th1 > (*cnet->th)(m,j) + (*SC->thi)(m,j)) th1 = (*cnet->th)(m,j) +
                                                                                (*SC->thi)(m,j);
                    }
                  else
                    {
                      th0 = (*sl->th)(m,j);
                      th1 = teta_psi(Psif(Fmin(egy->Tstar->co[m],egy->Temp->co[l])), 0.0,
                                     (*sl->pa)(sy,jsat,m),
                                     (*sl->pa)(sy,jres,m), (*sl->pa)(sy,ja,m), (*sl->pa)(sy,jns,m),
                                     1-1/(*sl->pa)(sy,jns,m), PsiMin, (*sl->pa)(sy,jss,m));
                      if (th1 > (*sl->th)(m,j) + (*SL->thi)(m,j)) 
                          th1 = (*sl->th)(m,j) + (*SL->thi)(m,j);
                    }

                  egy->deltaw->co[l]=(th1-th0)*egy->Dlayer->co[l]*rho_w;


                  //snow
                }
              else if (l>0)
                {

                  th0=egy->liq->co[l]/(egy->ice->co[l]+egy->liq->co[l]);
                  th1=theta_snow(par->alpha_snow, 1., egy->Temp->co[l]);

                  egy->deltaw->co[l]=(th1-th0)*(egy->ice->co[l]+egy->liq->co[l]);

                }

            }

          Tg=egy->Temp->co[sur];
          if (Tg < Tmin_surface_below_which_surfenergy_balance_recalculated) flagTmin++;

          if (surfacebalance == 1)
            {
              if (cont < num_iter_after_which_surfenergy_balance_not_recalculated
                  || flagTmin>0)
                {

                  //update egy->THETA taking into account evaporation (if there is not snow)
                  if (ns+ng == 0)
                    {
                      for (l=ns+ng+1; l<=n; l++)
                        {

                          m=l-ns-ng;

                          if (i<=par->total_channel)
                            {
                              egy->THETA->co[m] = (*cnet->th)(m,j) + egy->deltaw->co[l]/(rho_w*egy->Dlayer->co[l]);
                            }
                          else
                            {
                              egy->THETA->co[m] = (*sl->th)(m,j) + egy->deltaw->co[l]/(rho_w*egy->Dlayer->co[l]);
                            }

                          //add canopy transpiration
                          if (egy->THETA->co[m] > (*sl->pa)(sy,jres,1) + 1.E-3
                              && l <= egy->soil_transp_layer->nh )
                            {
                              egy->THETA->co[m] -= Fmax( Dt*fc*egy->soil_transp_layer->co[m]/
                                                         (rho_w*egy->Dlayer->co[l]), 0.0 );
                              if (egy->THETA->co[m] < (*sl->pa)(sy,jres,m)+1.E-3) egy->THETA->co[m] =
                                  (*sl->pa)(sy,jres,m)+1.E-3;
                            }

                          //add soil evaporation
                          if (egy->THETA->co[m] > (*sl->pa)(sy,jres,1) + 1.E-3
                              && l <= egy->soil_evap_layer_bare->nh )
                            {
                              egy->THETA->co[m] -= Fmax( Dt*(1.-fc)*egy->soil_evap_layer_bare->co[m]/
                                                         (rho_w*egy->Dlayer->co[l]), 0.0 );
                              egy->THETA->co[m] -= Fmax( Dt*fc*egy->soil_evap_layer_veg->co[m]/
                                                         (rho_w*egy->Dlayer->co[l]), 0.0 );
                              if (egy->THETA->co[m] < (*sl->pa)(sy,jres,m)+1.E-3) egy->THETA->co[m] =
                                  (*sl->pa)(sy,jres,m)+1.E-3;
                            }
                        }
                    }

                  //surface energy balance
                  EnergyFluxes_no_rec_turbulence(t, Tg, r, c, ns+ng, egy->T0->co[1], Qg0, Tv0, zmu, zmT, // 10 parameters
                                                 z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, // 10 parameters
                                                 P, LR, (*SL->P)(0,j), eps, fc, LSAI, decaycoeff0, *Wcrn, Wcrnmax, *Wcsn,  // 10 parameters
                                                 Wcsnmax, &dWcrn, &dWcsn, &egy->THETA->co[0], sl->pa->co[sy], // 5 parameters
                                                 land->ty->row(lu), land->root_fraction->row(lu), par, egy->soil_transp_layer.get(), SWin, // 5 parameters
                                                 LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, // 10 parameters
                                                 Etrans, &(V->Tv->co[j]), Qv, Ts, Qs, Hg0, Hg1, Eg0, Eg1, Lob, // 10 parameters
                                                 rh, rv, rc, rb, ruc, &rh_g, &rv_g, Qg, u_top, decay, // 10 parameters
                                                 Locc, LWup_ab_v, &egy->Temp->co[0], egy->soil_evap_layer_bare.get(), egy->soil_evap_layer_veg.get(), // 5 parameters
                                                 (*top->sky)(r,c), flagTmin, cont); // 3 parameters => TOT = 78 parameters

                  if (micro == 1)
                    {
                      EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Eadd + egy->SWlayer->co[0];
                      dEB_dT = -eps*(*top->sky)(r,c)*dSB_dT(Tg) - dH_dT - latent(Tg,
                               Levap(Tg))*dE_dT;//changed
                    }

                }
            }

          if (micro_sempl == 1)  EB = Convd * (Ta - egy->Temp->co[sur]);

          //F(T) = diag(egy->Fenergy(T)) + K(T)*T
          for (l=sur; l<=n; l++)
            {
              egy->Fenergy->co[l] = 0.0;
            }

          *dUsl = 0.;

          for (l=1; l<=n; l++)
            {

              //conductive M-matrix (lower diagonal part of this M-matrix is stored in a vector)

              thw = (egy->liq->co[l]+egy->deltaw->co[l])/(rho_w*egy->Dlayer->co[l]);
              thi = (egy->ice->co[l]-egy->deltaw->co[l])/(rho_i*egy->Dlayer->co[l]);

              if (l > ns+ng)
                {
                  sat = (*sl->pa)(sy,jsat,l-ns-ng);
                  kt = (*sl->pa)(sy,jkt,l-ns-ng);
                }
              else
                {
                  sat = 1.;
                  kt = 0.;
                }

              if (l<=ns+ng)
                {
                  kub1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
                }
              else
                {
                  kub1 = k_thermal(0, 1, thw, thi, sat, kt);
                }

              if (l>1)
                {

                  thwn = (egy->liq->co[l-1]+egy->deltaw->co[l-1])/(rho_w*egy->Dlayer->co[l-1]);
                  thin = (egy->ice->co[l-1]-egy->deltaw->co[l-1])/(rho_i*egy->Dlayer->co[l-1]);

                  if (l-1 > ns+ng)
                    {
                      satn = (*sl->pa)(sy,jsat,l-1-ns-ng);
                      ktn = (*sl->pa)(sy,jkt,l-1-ns-ng);
                    }
                  else
                    {
                      satn = 1.;
                      ktn = 0.;
                    }

                  thwn = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], thw, thwn);
                  thin = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], thi, thin);
                  satn = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], sat, satn);
                  ktn = Arithmetic_Mean(egy->Dlayer->co[l], egy->Dlayer->co[l-1], kt, ktn);

                  if (l <= ns+ng)
                    {
                      egy->Kth1->co[l-1] = - 2. * k_thermal(1, par->snow_conductivity, thwn, thin,
                                                            satn, ktn) / ( egy->Dlayer->co[l-1] + egy->Dlayer->co[l] );
                    }
                  else
                    {
                      egy->Kth1->co[l-1] = - 2. * k_thermal(0, 1, thwn, thin, satn,
                                                            ktn) / ( egy->Dlayer->co[l-1] + egy->Dlayer->co[l] );
                    }

                }
              else if (l>sur)
                {

                  if (l <= ns+ng)
                    {
                      egy->Kth1->co[l-1] = - k_thermal(1, par->snow_conductivity, thw, thi, sat,
                                                       kt) / ( egy->Dlayer->co[l]/2. );
                    }
                  else
                    {
                      egy->Kth1->co[l-1] = - k_thermal(0, 1, thw, thi, sat,
                                                       kt) / ( egy->Dlayer->co[l]/2. );
                    }

                }

              //diagonal part of F due to heat capacity and boundary condition (except conduction)
              C0 = calc_C(l, ns+ng, 0.0, &egy->ice->co[0], &egy->liq->co[0], &egy->deltaw->co[0],
                          &egy->Dlayer->co[0], sl->pa->co[sy]);
              C1 = calc_C(l, ns+ng, 1.0, &egy->ice->co[0], &egy->liq->co[0], &egy->deltaw->co[0],
                          &egy->Dlayer->co[0], sl->pa->co[sy]);
              if (l<=ns+ng
                  && egy->ice->co[l]-egy->deltaw->co[l]<1.E-7) C1 = Csnow_at_T_greater_than_0;
              egy->Fenergy->co[l] += ( Lf*egy->deltaw->co[l] +
                                       C1*egy->Dlayer->co[l]*egy->Temp->co[l] -
                                       C0*egy->Dlayer->co[l]*egy->T0->co[l] ) / Dt;

              //shortwave radiation penetrating under the surface
              if (micro == 1 && l<=ns+1) egy->Fenergy->co[l] -= egy->SWlayer->co[l];

              //store soil internal energy
              if (l>ns+ng) *dUsl = *dUsl + Lf*egy->deltaw->co[l] +
                                     C1*egy->Dlayer->co[l]*egy->Temp->co[l] - C0*egy->Dlayer->co[l]*egy->T0->co[l];

            }


          //top boundary condition
          egy->Fenergy->co[sur] -= ( (1.-KNe)*EB + KNe*EB0 );
          if (dirichlet == 1) egy->Fenergy->co[sur] -= ( (Tgd-egy->Temp->co[sur])*
                                                           (1.-KNe)*kub1 + (Tgd-egy->T0->co[sur])*KNe*kub0 ) / (egy->Dlayer->co[1]/2.);

          //bottom boundary condition (treated as sink)
          egy->Fenergy->co[n] -= par->Fboundary;

          if (n <= ns+ng)
            {
              kbb1 = k_thermal(1, par->snow_conductivity, thw, thi, sat, kt);
            }
          else
            {
              kbb1 = k_thermal(0, 1, thw, thi, sat, kt);
            }

          if (dirichlet_bottom == 1)
            {
              egy->Fenergy->co[n] -= ( (Tbottom-egy->Temp->co[n])*(1.-KNe)*kbb1 +
                                       (Tbottom-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlayer->co[n]/2.);
            }
          else
            {
              egy->Fenergy->co[n] -= ( (par->Tboundary-egy->Temp->co[n])*(1.-KNe)*kbb1 +
                                       (par->Tboundary-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlayer->co[n]/2.
                                           +par->Zboundary);
            }

          update_F_energy(sur, n, egy->Fenergy.get(), 1.-KNe, egy->Kth1.get(), &egy->Temp->co[0]);
          update_F_energy(sur, n, egy->Fenergy.get(), KNe, egy->Kth0.get(), &egy->T0->co[0]);
          res = norm_2(egy->Fenergy.get(), sur, n);

          //printf("Dt:%.2f res:%e cont:%ld cont2:%ld T0:%f T1:%f EB:%f SW:%f Ta:%f Tg:%f LW:%.2f H:%.2f ET:%.2f dw:%f liq:%f ice:%f \n",Dt,res,cont,cont2,egy->Temp->co[sur],egy->Temp->co[1],EB,egy->SWlayer->co[0],Ta,Tg,*LW,(*H),latent(Tg,Levap(Tg))*(*E),egy->deltaw->co[1],egy->liq->co[1]+egy->deltaw->co[1],egy->ice->co[1]-egy->deltaw->co[1]);
          //printf("res:%e sur:%ld\n",res,sur);

          if (res <= res_av*(1.0 - ni*lambda[0])) iter_close2=1;
          if (lambda[0] <= par->min_lambda_en) cont_lambda_min++;
          if (cont_lambda_min > par->max_times_min_lambda_en)
            {
              if (par->exit_lambda_min_en == 1)
                {
                  return 1;
                }
              else
                {
                  iter_close2=1;
                  cont_lambda_min=0;
                }
            }

        }
      while (iter_close2!=1);

      if (iter_close>=0 && res<=par->tol_energy) iter_close=1;
      if (cont>=par->maxiter_energy) iter_close=1;

    }
  while (iter_close!=1);

  //printf("res:%e sur:%ld\n",res,sur);

  //if there is no convergence, go out of the loop
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

  //skin temperature
  if (sur == 1
      && par->nsurface == 0) egy->Temp->co[par->nsurface] = egy->Temp->co[sur];

  //it is not admitted that the skin layer temperature goes > 0 when snow or glacier are present
  if (sur == 0 && ns+ng > 0 && egy->Temp->co[sur] > 0) return -1;

  //it is not admitted that a snow/glacier layer melts completely
  for (l=1; l<=ns+ng; l++)
    {

      //snow layer not alone
      if (ns > 1 && l >= 1 && l <= ns
          && egy->ice->co[l] - egy->deltaw->co[l] < 1.E-5)
        {
          *lpb = l;
          return -2;

          //glacier layer not alone
        }
      else if (ng > 1 && l >= ns+1 && l<= ns+ng
               && egy->ice->co[l] - egy->deltaw->co[l] < 1.E-5)
        {
          *lpb = l;
          return -3;

          //just one glacier layer and snow layers >= 0
        }
      else if (ng == 1 && l == ns+ng
               && egy->ice->co[l] - egy->deltaw->co[l] < 1.E-5)
        {
          return -4;

          //just one snow layer and glacier layers > 0
        }
      else if ( ns == 1 && ng > 0 && l == ns
                && egy->ice->co[l] - egy->deltaw->co[l] < 1.E-5)
        {
          return -5;

          //just one snow layer and glacier layers == 0
        }
      else if ( ns == 1 && ng == 0 && l == ns
                && egy->ice->co[l] - egy->deltaw->co[l] < 1.E-5)
        {
          return -6;
        }
    }

  //error messages
  for (l=sur; l<=Nl+ns+ng; l++)
    {

      if (egy->Temp->co[l] != egy->Temp->co[l])
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t);
          fprintf(f, "T no value, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (egy->Temp->co[l] < Tmin || egy->Temp->co[l] > Tmax)
        {
          printf("Warning, T outside of range, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f LW:%f H:%f LE:%f Ta:%f\n",
                 l,r,c,egy->Temp->co[l],*LW,*H,latent(Tg,Levap(Tg))*(*E),Ta);
        }
    }

  for (l=1; l<=Nl+ns+ng; l++)
    {

      if (egy->deltaw->co[l] != egy->deltaw->co[l])
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t);
          fprintf(f, "dw no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,
                  r,c);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

    }

  //update canopy water
  *Wcrn = *Wcrn + dWcrn;
  *Wcsn = *Wcsn + dWcsn;

  return 0;

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_soil_land(long nsurf, long n, long i, long r, long c, double fc,
                      double Dt, ENERGY *egy, MatrixView<double> &&pa, SOIL_STATE *S, Tensor<double> *ET,
                      Matrix<double> *th)
{

  long l;
  double th_oversat, psisat;

  //soil variables
  for (l=1; l<=Nl; l++)
    {

      //canopy transpiration
      if (l <= egy->soil_transp_layer->nh)
        (*ET)(l,r,c) +=  fc*egy->soil_transp_layer->co[l]*Dt;

      //soil evaporation
      if (l <= egy->soil_evap_layer_bare->nh)
        {
          (*ET)(l,r,c) += (1.-fc)*egy->soil_evap_layer_bare->co[l]*Dt;
          (*ET)(l,r,c) += fc*egy->soil_evap_layer_veg->co[l]*Dt;
        }

      //water pressure and contents
      psisat = psi_saturation(Fmax(0.,
                                   egy->ice->co[l+n])/(rho_w*egy->Dlayer->co[l+n]), pa(jsat,l), pa(jres,l),
                              pa[ja][l], pa[jns][l], 1.-1./pa[jns][l]);
      th_oversat = Fmax( (*S->P)(l,i) - psisat, 0.0 ) * pa[jss][l];

      (*th)(l,i) = Fmax(0.,egy->liq->co[l+n]+egy->deltaw->co[l+n])/
                     (rho_w*egy->Dlayer->co[l+n]);
      (*S->thi)(l,i) = Fmax(0.,
                              egy->ice->co[l+n]-egy->deltaw->co[l+n])/(rho_w*egy->Dlayer->co[l+n]);

      (*S->P)(l,i) = psi_teta((*th)(l,i)+th_oversat, (*S->thi)(l,i),
                                pa(jsat,l), pa(jres,l), pa[ja][l], pa[jns][l], 1.-1./pa[jns][l], PsiMin,
                                pa[jss][l]);

      //temperature
        (*S->T)(l,i) = egy->Temp->co[l+n];

    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_soil_channel(long nsurf, long n, long ch, double fc, double Dt,
                         ENERGY *egy, MatrixView<double> &&pa, SOIL_STATE *S, Matrix<double> *ET, Matrix<double> *th)
{

  long l;
  double th_oversat, psisat;

  //soil variables
  for (l=1; l<=Nl; l++)
    {

      //canopy transpiration
      if (l <= egy->soil_transp_layer->nh) (*ET)(l,ch) +=
          fc*egy->soil_transp_layer->co[l]*Dt;

      //soil evaporation
      if (l <= egy->soil_evap_layer_bare->nh)
        {
          (*ET)(l,ch) += (1.-fc)*egy->soil_evap_layer_bare->co[l]*Dt;
          (*ET)(l,ch) += fc*egy->soil_evap_layer_veg->co[l]*Dt;
        }

      //water pressure and contents
      psisat = psi_saturation((*S->thi)(l,ch), pa(jsat,l), pa(jres,l),
                              pa[ja][l], pa[jns][l], 1.-1./pa[jns][l]);
      th_oversat = Fmax( (*S->P)(l,ch) - psisat, 0.0 ) * pa[jss][l];

      (*th)(l,ch) = Fmax(0., egy->liq->co[l+n]+egy->deltaw->co[l+n])/(rho_w*egy->Dlayer->co[l+n]);
      (*S->thi)(l,ch) = Fmax(0., egy->ice->co[l+n]-egy->deltaw->co[l+n])/(rho_w*egy->Dlayer->co[l+n]);

      (*S->P)(l,ch) = psi_teta((*th)(l,ch)+th_oversat, (*S->thi)(l,ch),
                                 pa(jsat,l), pa(jres,l), pa[ja][l], pa[jns][l], 1.-1./pa[jns][l], PsiMin,
                                 pa[jss][l]);

      //temperature
        (*S->T)(l,ch) = egy->Temp->co[l+n];

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_F_energy(long nbeg, long nend, Vector<double> *F, double w,
                     Vector<double> *K, double *T)
{

  long l;

  for (l=nbeg; l<=nend; l++)
    {

      if (l==nbeg)
        {
          F->co[l] += w*(-K->co[l]*T[l] + K->co[l]*T[l+1]);
        }
      else if (l<nend)
        {
          F->co[l] += w*(K->co[l-1]*T[l-1] - (K->co[l]+K->co[l-1])*T[l] + K->co[l]*T[l+1]);
        }
      else
        {
          F->co[l] += w*(K->co[l-1]*T[l-1] - K->co[l-1]*T[l]);
        }

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_diag_dF_energy(long nbeg, long nend, Vector<double> *dF, double w,
                           Vector<double> *K)
{

  long l;

  for (l=nbeg; l<=nend; l++)
    {

      if (l==nbeg)
        {
          dF->co[l] -= w*K->co[l];
        }
      else if (l<nend)
        {
          dF->co[l] -= w*(K->co[l]+K->co[l-1]);
        }
      else
        {
          dF->co[l] -= w*K->co[l-1];
        }

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double calc_C(long l, long nsng, double a, double *wi, double *wl, double *dw,
              double *D, double **pa)
{

  double C;

  if (l<=nsng)  //snow
    {
      C = (c_ice*(wi[l]-a*dw[l]) + c_liq*(wl[l]+a*dw[l]))/D[l];
    }
  else    //soil
    {
      C = pa[jct][l-nsng]*(1.-pa[jsat][l-nsng]) + (c_ice*(wi[l]-a*dw[l]) + c_liq*
                                                   (wl[l]+a*dw[l]))/D[l];
    }

  return (C);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void EnergyFluxes(double t, double Tg, long r, long c, long n, // 5 parameters
                  double Tg0, double Qg0, double Tv0, double zmu, double zmT, // 5 parameters
                  double z0s, double d0s, double rz0s, double z0v, double d0v, // 5 parameters
                  double rz0v, double hveg, double v, double Ta, double Qa, // 5 parameters
                  double P, double LR, double psi, double e, double fc, // 5 parameters
                  double LSAI, double decaycoeff0, double Wcrn, double Wcrnmax, double Wcsn, // 5 parameters
                  double Wcsnmax, double *dWcrn, double *dWcsn, double *theta, double **soil, // 5 parameters
                  RowView<double> &&land, RowView<double> &&root, PAR *par, Vector<double> *soil_transp_layer,
                  double SWin, // 5 parameters
                  double LWin, double SWv, double *LW, double *H, double *dH_dT, // 5 parameters
                  double *E, double *dE_dT, double *LWv, double *Hv, double *LEv, // 5 parameters
                  double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs, // 5 parameters
                  double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Lobukhov, // 5 parameters
                  double *rh, double *rv, double *rc, double *rb, double *ruc, // 5 parameters
                  double *rh_g, double *rv_g, double *Qg, double *u_top, double *decay, // 5 parameters
                  double *Locc, double *LWup_above_v, double *T, Vector<double> *soil_evap_layer_bare,
                  Vector<double> *soil_evap_layer_veg, // 5 parameters
                  double sky) // 1 parameters
{


  double Hg, dHg_dT, Eg, dEg_dT, LWg, LWup;
  double dQgdT, rm;
  double alpha, beta; //from Ye and Pielke, 1993
  FILE *f;

  *Tv=Tg;

  *H=0.0;
  *E=0.0;

  *dH_dT=0.0;
  *dE_dT=0.0;

  *LW=0.0;
  *LWup_above_v=0.0;

  *LWv=0.0;
  *Hv=0.0;
  *LEv=0.0;

  *Etrans=0.0;
  *Hg0=0.0;
  *Eg0=0.0;
  *Hg1=0.0;
  *Eg1=0.0;

  *decay=0.0;
  *Locc=0.0;

  //thermodynamical calculations
  SpecHumidity_2(Qg, &dQgdT, 1.0, Tg, P);

  //initialization
  if (fc==0)
    {

      *Qs=(double)number_novalue;
      *Ts=(double)number_novalue;
      *Qv=(double)number_novalue;
      *u_top=(double)number_novalue;
      *rh=1.E99;
      *rv=1.E99;
      *rc=1.E99;
      *rb=1.E99;
      *ruc=1.E99;

    }
  else if (fc==1)
    {

      *rh_g=1.E99;
      *rv_g=1.E99;

    }

  //bare soil
  if (fc<1)
    {
      aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, LR, Lobukhov,
                      &rm, rh_g, rv_g, par->state_turb, par->monin_obukhov, par->maxiter_Businger);

      find_actual_evaporation_parameters(r,c,&alpha, &beta, soil_evap_layer_bare,
                                         theta, soil, T, psi, P, *rv_g, Ta, Qa, *Qg, n);
      turbulent_fluxes(*rh_g, *rv_g/beta, P, Ta, Tg, Qa, *Qg*alpha, dQgdT*alpha,
                       &Hg, &dHg_dT, &Eg, &dEg_dT);

      *H+=(1.0-fc)*Hg;
      *E+=(1.0-fc)*Eg;

      *dH_dT+=(1.0-fc)*dHg_dT;
      *dE_dT+=(1.0-fc)*dEg_dT;

      *LW += (1.0-fc)*( e*(LWin-sky*SB(Tg)) );//changed
      *LWup_above_v+=(1.0-fc)*( (1.0-e)*LWin+e*SB(Tg) );

      *Hg0=Hg;
      *Eg0=Eg;

      if (Hg!=Hg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f,
                  "Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,
                  c,*rh_g,*rv_g,Hg,Eg);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (Eg!=Eg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f, "Eg no value bare soil %ld %ld \n",r,c);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

    }

  //vegetation
  if (fc>0)
    {

      Tcanopy(r, c, Tv0, Tg, *Qg, dQgdT, Tg0, Qg0, Ta, Qa, // 10
              zmu, zmT, z0v, z0s, d0v, rz0v, hveg, v, LR, P, // 10
              SWin, SWv, LWin, e, LSAI, decaycoeff0, std::forward<RowView<double>>(land), Wcrn, Wcrnmax, Wcsn, // 10
              Wcsnmax, dWcrn, dWcsn, LWv, &LWg, Hv, &Hg, &dHg_dT, LEv, &Eg, // 10
              &dEg_dT, Ts, Qs, std::forward<RowView<double>>(root), theta, soil_transp_layer, Lobukhov, par, n, &rm, // 10
              rh, rv, rc, rb, ruc, u_top, Etrans, Tv, Qv, decay, // 10
              Locc, &LWup, psi, soil, T, soil_evap_layer_veg); // 6


      *H+=fc*Hg;
      *E+=fc*Eg;

      *dH_dT+=fc*dHg_dT;
      *dE_dT+=fc*dEg_dT;

      *LW+=fc*LWg;
      *LWup_above_v+=fc*LWup;

      *Hg1=Hg;
      *Eg1=Eg;

      if (Hg!=Hg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f,
                  "Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",
                  r,c,*rh_g,*rv_g,Hg,Eg);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (Eg!=Eg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f, "Eg no value vegetated soil %ld %ld \n",r,c);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void EnergyFluxes_no_rec_turbulence(double t, double Tg, long r, long c,
                                    long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT,
                                    double z0s,
                                    double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg,
                                    double v, double Ta, double Qa,
                                    double P, double LR, double psi, double e, double fc, double LSAI,
                                    double decaycoeff0, double Wcrn,
                                    double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn, double *dWcsn,
                                    double *theta, double **soil,
                                    RowView<double> &&land, RowView<double> &&root, PAR *par,
                                    Vector<double> *soil_transp_layer,
                                    double SWin, double LWin, double SWv, double *LW,
                                    double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv,
                                    double *LEv, double *Etrans,
                                    double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1,
                                    double *Eg0, double *Eg1, double *Lobukhov,
                                    double *rh, double *rv, double *rc, double *rb, double *ruc, double *rh_g,
                                    double *rv_g, double *Qg,
                                    double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T,
                                    Vector<double> *soil_evap_layer_bare,
                                    Vector<double> *soil_evap_layer_veg, double sky, short flagTmin, long cont)
{


  double Hg, dHg_dT, Eg, dEg_dT, LWg, LWup;
  double dQgdT;
  double dLWvdT,LWv_old;
  double alpha, beta; //from Ye and Pielke, 1993
  double rm=1.E20;
  FILE *f;

  //initalization
  *H=0.0;
  *E=0.0;

  *dH_dT=0.0;
  *dE_dT=0.0;

  *LW=0.0;
  *LWup_above_v=0.0;

  //thermodynamical calculations
  SpecHumidity_2(Qg, &dQgdT, 1.0, Tg, P);

  if (fc<1)
    {

      if (flagTmin>0
          && cont>num_iter_after_which_only_neutrality) aero_resistance(zmu, zmT, z0s,
                d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, LR, Lobukhov, &rm, rh_g, rv_g,
                par->state_turb, 2, par->maxiter_Businger);

      find_actual_evaporation_parameters(r,c,&alpha, &beta, soil_evap_layer_bare,
                                         theta, soil, T, psi, P, *rv_g, Ta, Qa, *Qg, n);
      turbulent_fluxes(*rh_g, *rv_g/beta, P, Ta, Tg, Qa, *Qg*alpha, dQgdT*alpha,
                       &Hg, &dHg_dT, &Eg, &dEg_dT);

      *H+=(1.0-fc)*Hg;
      *E+=(1.0-fc)*Eg;

      *dH_dT+=(1.0-fc)*dHg_dT;
      *dE_dT+=(1.0-fc)*dEg_dT;

      *LW += (1.0-fc)*( e*(LWin-sky*SB(Tg)) );//changed
      *LWup_above_v+=(1.0-fc)*( (1.0-e)*LWin+e*SB(Tg) );

      *Hg0=Hg;
      *Eg0=Eg;

      //error messages
      if (Hg!=Hg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f,
                  "Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,
                  c,*rh_g,*rv_g,Hg,Eg);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (Eg!=Eg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f, "Eg no value bare soil %ld %ld \n",r,c);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

    }

  //TURBULENT FLUXES FOR CANOPY CANOPY
  if (fc>0)
    {

      find_actual_evaporation_parameters(r,c,&alpha, &beta, soil_evap_layer_veg,
                                         theta, soil, T, psi, P, *ruc, Ta, Qa, *Qg, n);
      *Ts=(Ta/(*rh)+Tg/(*ruc)+(*Tv)/(*rb))/(1./(*rh)+1./(*ruc)+1./(*rb));
      *Qs=(Qa/(*rv)+(*Qg)*alpha*beta/(*ruc)+(*Qv)/(*rc))/(1./(*rv)+beta/(*ruc)+1./
                                                          (*rc));
      turbulent_fluxes(*ruc, *ruc/beta, P, *Ts, Tg, *Qs, *Qg*alpha, dQgdT*alpha,
                       &Hg, &dHg_dT, &Eg, &dEg_dT);

      LWv_old=*LWv;
      longwave_vegetation(LWin, e, Tg, *Tv, LSAI, LWv, &LWg, &dLWvdT, &LWup);
      *LWv=LWv_old;

      *H+=fc*Hg;
      *dH_dT+=fc*dHg_dT;
      *E+=fc*Eg;
      *dE_dT+=fc*dEg_dT;

      *LW+=fc*LWg;
      *LWup_above_v+=fc*LWup;

      *Hg1=Hg;
      *Eg1=Eg;

      if (Hg!=Hg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f,
                  "Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",
                  r,c,*rh_g,*rv_g,Hg,Eg);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

      if (Eg!=Eg)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f, "Simulation Period:%ld\n",i_sim);
          fprintf(f, "Run Time:%ld\n",i_run);
          fprintf(f, "Number of days after start:%f\n",t/86400.);
          fprintf(f, "Eg no value vegetated soil %ld %ld \n",r,c);
          fprintf(f, "Ta:%f Tg:%f\n",Ta,Tg);
          fprintf(f,
                  "zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
                  zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double k_thermal(short snow, short a, double th_liq, double th_ice,
                 double th_sat, double k_solid)
{

  /*
   for soil:
   Quadratic parallel based on the analogy with electric lattice
   (Cosenza et al., 2003, European Journal of Soil Science, September 2003, 54, 581–587)

   for snow:
   if a==1 Cosenza
   if a==2 Sturm(1997)
   if a==3 Jordan(1991)
   */

  double r, k=0.;

  if (snow==0 || (a!=2 && a!=3))
    {
      k = pow ( (1.-th_sat)*sqrt(k_solid) + th_liq*sqrt(k_liq) + th_ice*sqrt(
                  k_ice) + (th_sat-th_liq-th_ice)*sqrt(k_air), 2. ) ;

    }
  else
    {

      r = th_liq * rho_w + th_ice * rho_i ;

      if (a==2)
        {
          if (r < 156)
            {
              k = 0.023 + 0.234 * r*1.E-3 ;
            }
          else
            {
              k = 0.138 - 1.01 * r*1.E-3 + 3.233 * r*r*1.E-6 ;
            }
        }
      else if (a==3)
        {
          k = k_air + (7.75E-5 * r + 1.105E-6 * r*r) * (k_ice-k_air) ;
        }

    }

  return (k);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


double flux(long i, long icol, double **met, double k, double add, double est)
{

  double F=met[i-1][icol];

  if ((long)F==number_absent || (long)F==number_novalue)
    {
      F=est;
    }
  else
    {
      F*=k;
      F+=add;
    }


  return (F);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD,
                           double thres, double z0snow, double *z0_ris,
                           double *d0_ris, double *z0_z0t_ris)
{

  if (snowD>thres)
    {
      *z0_ris=z0snow;
      *d0_ris=0.0;
      *z0_z0t_ris=0.0;
    }
  else
    {
      *z0_ris=z0;
      *d0_ris=d0;
      *z0_z0t_ris=z0_z0t;
    }

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void merge(double a, Vector<double> *ice, Vector<double> *liq, Vector<double> *Temp,
           Vector<double> *D, long l, long lup, long ldw, long tot)
{

  long i, n, m;
  double h;

  if (l == lup)
    {
      m = l+1;
      n = l;
    }
  else if (l < ldw)
    {
      if (ice->co[l+1] < ice->co[l-1])
        {
          m = l+1;
          n = l;
        }
      else
        {
          m = l;
          n = l-1;
        }
    }
  else
    {
      m = l;
      n = l-1;
    }

  h = internal_energy(ice->co[m], liq->co[m],
                      Temp->co[m]) + internal_energy(ice->co[n], liq->co[n], Temp->co[n]);

  ice->co[n] += ice->co[m];
  liq->co[n] += liq->co[m];
  D->co[n] += D->co[m];

  from_internal_energy(a, 0, 0, h, &(ice->co[n]), &(liq->co[n]),
                       &(Temp->co[n]));

  for (i=m+1; i<=tot; i++)
    {
      ice->co[i-1] = ice->co[i];
      liq->co[i-1] = liq->co[i];
      Temp->co[i-1] = Temp->co[i];
      D->co[i-1] = D->co[i];
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void sux_minus6_condition(double ic, double wa, double rho, double D1,
                          ENERGY *E)
{

  long l;

  E->Dlayer->co[1] = D1;

  for (l=Nl; l>=1; l--)
    {
      E->ice->co[l+1] = E->ice->co[l];
      E->liq->co[l+1] = E->liq->co[l];
      E->deltaw->co[l+1] = E->deltaw->co[l];
      E->Temp->co[l+1] = E->Temp->co[l];
      E->Dlayer->co[l+1] = E->Dlayer->co[l];
    }

  E->ice->co[1] = ic;
  E->liq->co[1] = wa;
  E->Temp->co[1] = Fmin(E->Temp->co[2], -0.1);
  E->Dlayer->co[1] = ic/rho;

  if (E->deltaw->co[2] > ic)
    {
      E->deltaw->co[1] = ic;
      E->deltaw->co[2] -= ic;
    }
  else if (E->deltaw->co[2] > 0)
    {
      E->deltaw->co[1] = E->deltaw->co[2];
      E->deltaw->co[2] = 0.;
    }

  E->ice->co[2] -= ic;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
