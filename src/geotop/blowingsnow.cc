
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


#include "blowingsnow.h"
#include "inputKeywords.h"

// Private Declarations

static void print_windtrans_snow(double Dt, Snow *snow, Par *par, Topo *top);
static void set_windtrans_snow(double Dt,
                               double t,
                               Snow *snow,
                               Meteo *met,
                               Land *land,
                               Par *par,
                               FILE *f);

/**
 * @ROUTINE windtrans_snow
 * @Author S. Endrizzi
 * @date   2009
 * @brief  snow trasport by wind
 *
 *
 */

void windtrans_snow(Snow *snow,
                    Meteo *met,
                    Land *land,
                    Topo *top,
                    Par *par,
                    double t0)
{
  long r, c, j, l;
  static long line_interp;
  short lu, yes, sux;
  (void)sux;  // silence warning about unused variable since sux appears only in
  // a if-branch
  size_t lux;
  double t, Dt, Dt0;
  double zmeas, D, wice, fsnow, fc;
  double canopy_height_over_snow, rho_snow_surface;
  double PBSM_fetch = 1000.;
  double dErdt;
  long ns;
  FILE *f;

  f = fopen(geotop::common::Variables::logfile.c_str(), "a");

  if (t0 == 0) line_interp = 0;

  t = 0.;

  // t0 = time from beginning simulation
  // t = time from beginning of time step

  do
    {
      // find Dt
      Dt = par->Dt_PBSM;
      if (t + Dt > par->Dt) Dt = par->Dt - t;
      Dt0 = Dt;

      // vegetation
      for (lux = 1; lux <= par->n_landuses; lux++)
        {
          if (par->vegflag[lux] == 1)
            {
              time_interp_linear(par->init_date + t0 / GTConst::secinday,
                                 par->init_date + (t0 + t) / GTConst::secinday,
                                 par->init_date + (t0 + t + Dt) / GTConst::secinday,
                                 land->vegparv[lux - 1], land->vegpars[lux - 1],
                                 land->NumlinesVegTimeDepData[lux - 1], jdvegprop + 1,
                                 0, 0, &line_interp);
            }
        }

      // loop for every pixel
      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
                  D = DEPTH(r, c, snow->S->lnum, snow->S->Dzl);
                  wice = DEPTH(r, c, snow->S->lnum, snow->S->w_ice);
                  if (wice > par->Wmin_BS)
                    {
                      // find canopy_height_over_snow
                      lu = (short)land->LC[r][c];
                      canopy_height_over_snow = 0.0;

                      zmeas = met->st->Vheight[1];

                      for (j = 1; j <= jdvegprop; j++)
                        {
                          if ((long)land->vegparv[lu - 1][j - 1 + 1] !=
                              geotop::input::gDoubleNoValue)
                            {
                              land->vegpar[j] = land->vegparv[lu - 1][j - 1 + 1];
                            }
                          else
                            {
                              land->vegpar[j] = land->ty[lu][j + jHveg - 1];
                            }
                        }

                      if (D > land->vegpar[jdz0thresveg])
                        {
                          fsnow = 1.0;
                        }
                      else if (D > land->vegpar[jdz0thresveg2])
                        {
                          fsnow =
                            (D - land->vegpar[jdz0thresveg2]) /
                            (land->vegpar[jdz0thresveg] - land->vegpar[jdz0thresveg2]);
                        }
                      else
                        {
                          fsnow = 0.0;
                        }
                      // canopy fraction
                      fc = land->vegpar[jdcf] * pow(1.0 - fsnow, land->vegpar[jdexpveg]);
                      if (land->vegpar[jdLSAI] < GTConst::LSAIthres) fc = 0.0;
                      if (fc > 0)
                        canopy_height_over_snow +=
                          fc * Fmax(land->vegpar[jdHveg] - D, 0.) * 1.E-3;

                      // rearrange snow layers
                      ns = snow->S->lnum[r][c];
                      sux = copy_statevar_from3D_to1D(r, c, snow->S, snow->S_for_BS);

                      if (snow->S_for_BS->w_ice[ns] < par->Wice_PBSM)
                        {
                          l = ns;
                          do
                            {
                              l--;
                              sux = set_snowice_min(snow->S_for_BS, ns, l, par->Wice_PBSM);
                            }
                          while (l > 1 && snow->S_for_BS->w_ice[ns] < par->Wice_PBSM);
                        }

                      // find equilibrium snow trasport
                      rho_snow_surface =
                        (snow->S_for_BS->w_ice[ns] + snow->S_for_BS->w_liq[ns]) /
                        (1.E-3 * snow->S_for_BS->Dzl[ns]);

                      Pbsm(r, c, PBSM_fetch, land->ty[lu][jN], 1.E-3 * land->ty[lu][jdv],
                           canopy_height_over_snow, rho_snow_surface, zmeas,
                           met->Vgrid[r][c], met->Tgrid[r][c], met->RHgrid[r][c],
                           &(snow->Qtrans[r][c]), &(snow->Qsub[r][c]),
                           &(snow->Qsalt[r][c]), D, top->slope[r][c]);
                    }
                  else
                    {
                      snow->Qtrans[r][c] = 0.0;
                      snow->Qsub[r][c] = 0.0;
                    }
                }
            }
        }

      set_inhomogeneous_fetch(snow, met, land, par, top, &yes);

      if (yes == 1)
        {
          for (r = 1; r <= geotop::common::Variables::Nr; r++)
            {
              for (c = 1; c <= geotop::common::Variables::Nc; c++)
                {
                  if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                    {
                      wice = DEPTH(r, c, snow->S->lnum, snow->S->w_ice);
                      dErdt =
                        sqrt(pow(snow->Qsub_x[r][c], 2.) + pow(snow->Qsub_y[r][c], 2.)) -
                        snow->Nabla2_Qtrans[r][c];

                      if (snow->S->lnum[r][c] > 0)
                        {
                          if (dErdt * Dt > 0.99 * wice) Dt = 0.99 * wice / dErdt;
                        }
                    }
                }
            }

          Dt = Fmin(Dt, Dt0);
          set_windtrans_snow(Dt, t0 + t, snow, met, land, par, f);
          print_windtrans_snow(Dt, snow, par, top);
        }

      t += Dt;

    }
  while (t < par->Dt);

  fclose(f);
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

void set_inhomogeneous_fetch(Snow *snow,
                             Meteo *met,
                             Land *land,
                             Par *par,
                             Topo *top,
                             short *yes)
{
  long i, r, c, num_change, r0, c0;
  double dx, dy, F1, F2;
  double Qtrans = 0.0;
  double Qup, Qdown, Sup, Sdown, k_snowred;

  dx = geotop::common::Variables::UV->U[1];
  dy = geotop::common::Variables::UV->U[2];

  F1 = par->fetch_up / 3.;
  F2 = par->fetch_down / 3.;

  *yes = 0;

  snow->Nabla2_Qtrans.resize(snow->Nabla2_Qtrans.getRows(),
                             snow->Nabla2_Qtrans.getCols(), 0.0);

  // check if there is snow transport
  for (r = 1; r <= geotop::common::Variables::Nr; r++)
    {
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
            {
              Qtrans += snow->Qtrans[r][c] / (double)par->total_pixel;
            }
        }
    }

  // if there is blowing snow
  if (Qtrans > 0)
    {
      *yes = 1;

      // wind in direction west-east
      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          // find west-east component
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              snow->Qtrans_x[r][c] = fabs(
                                       snow->Qtrans[r][c] * (-sin(met->Vdir[r][c] * GTConst::Pi / 180.)));
              snow->Qsub_x[r][c] =
                fabs(snow->Qsub[r][c] * (-sin(met->Vdir[r][c] * GTConst::Pi / 180.)));
            }

          // find when there is a wind direction inversion
          snow->change_dir_wind.resize(snow->change_dir_wind.size(), 0.0);

          num_change = 0;
          c = 1;

          num_change++;
          c0 = c;
          snow->change_dir_wind[num_change] = c;

          do
            {
              c = c0;
              do
                {
                  c++;
                }
              while ((-sin(met->Vdir[r][c] * GTConst::Pi / 180.)) *
                     (-sin(met->Vdir[r][c0] * GTConst::Pi / 180.)) >
                     0 &&
                     c < geotop::common::Variables::Nc);

              num_change++;
              c0 = c;
              snow->change_dir_wind[num_change] = c;

            }
          while (c0 < geotop::common::Variables::Nc);

          for (i = 1; i < num_change; i++)
            {
              // east wind
              if ((-sin(met->Vdir[r][snow->change_dir_wind[i]] * GTConst::Pi /
                        180.)) > 0)
                {
                  if (par->upwindblowingsnow == 1 && snow->change_dir_wind[i] != 1)
                    snow->Qtrans_x[r][snow->change_dir_wind[i]] = 0.0;
                  snow->Qtrans_x[r][snow->change_dir_wind[i]] = 0.0;
                  for (c = snow->change_dir_wind[i] + 1;
                       c <= snow->change_dir_wind[i + 1]; c++)
                    {
                      if (snow->change_dir_wind[i + 1] == geotop::common::Variables::Nc ||
                          (snow->change_dir_wind[i + 1] !=
                           geotop::common::Variables::Nc &&
                           c < snow->change_dir_wind[i + 1]))
                        {
                          Qup = snow->Qtrans_x[r][c - 1];
                          Qdown = snow->Qtrans_x[r][c];
                          Sup = snow->Qsub_x[r][c - 1];
                          Sdown = snow->Qsub_x[r][c];
                          if (Qdown >= Qup)
                            {
                              snow->Qtrans_x[r][c] = (Qdown + Qup * F1 / dx) / (1. + F1 / dx);
                              snow->Qsub_x[r][c] = (Sdown + Sup * F1 / dx) / (1. + F1 / dx);
                            }
                          else
                            {
                              snow->Qtrans_x[r][c] = (Qdown + Qup * F2 / dx) / (1. + F2 / dx);
                              snow->Qsub_x[r][c] = (Sdown + Sup * F2 / dx) / (1. + F2 / dx);
                            }
                          Qdown = snow->Qtrans_x[r][c];
                          snow->Nabla2_Qtrans[r][c] += (Qup - Qdown) / dx;
                        }
                    }

                  // west wind
                }
              else
                {
                  if (par->upwindblowingsnow == 1 &&
                      snow->change_dir_wind[i + 1] != geotop::common::Variables::Nc)
                    snow->Qtrans_x[r][snow->change_dir_wind[i + 1] - 1] = 0.0;
                  snow->Qtrans_x[r][snow->change_dir_wind[i + 1] - 1] = 0.0;
                  for (c = snow->change_dir_wind[i + 1] - 1;
                       c >= snow->change_dir_wind[i]; c--)
                    {
                      if (snow->change_dir_wind[i + 1] == geotop::common::Variables::Nc ||
                          (snow->change_dir_wind[i + 1] !=
                           geotop::common::Variables::Nc &&
                           c < snow->change_dir_wind[i + 1] - 1))
                        {
                          Qup = snow->Qtrans_x[r][c + 1];
                          Qdown = snow->Qtrans_x[r][c];
                          Sup = snow->Qsub_x[r][c + 1];
                          Sdown = snow->Qsub_x[r][c];
                          if (Qdown >= Qup)
                            {
                              snow->Qtrans_x[r][c] = (Qdown + Qup * F1 / dx) / (1. + F1 / dx);
                              snow->Qsub_x[r][c] = (Sdown + Sup * F1 / dx) / (1. + F1 / dx);
                            }
                          else
                            {
                              snow->Qtrans_x[r][c] = (Qdown + Qup * F2 / dx) / (1. + F2 / dx);
                              snow->Qsub_x[r][c] = (Sdown + Sup * F2 / dx) / (1. + F2 / dx);
                            }
                          Qdown = snow->Qtrans_x[r][c];
                          snow->Nabla2_Qtrans[r][c] += (Qup - Qdown) / dx;
                        }
                    }
                }
            }
        }

      // wind in direction south-north
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          // find south-north component
          for (r = 1; r <= geotop::common::Variables::Nr; r++)
            {
              snow->Qtrans_y[r][c] = fabs(
                                       snow->Qtrans[r][c] * (-cos(met->Vdir[r][c] * GTConst::Pi / 180.)));
              snow->Qsub_y[r][c] =
                fabs(snow->Qsub[r][c] * (-cos(met->Vdir[r][c] * GTConst::Pi / 180.)));
            }

          // find when there is a wind direction inversion
          snow->change_dir_wind.resize(snow->change_dir_wind.size(), 0);
          num_change = 0;
          r = 1;

          num_change++;
          r0 = r;
          snow->change_dir_wind[num_change] = r;

          do
            {
              r = r0;
              do
                {
                  r++;
                }
              while ((-cos(met->Vdir[r][c] * GTConst::Pi / 180.)) *
                     (-cos(met->Vdir[r0][c] * GTConst::Pi / 180.)) >
                     0 &&
                     r < geotop::common::Variables::Nr);

              num_change++;
              r0 = r;
              snow->change_dir_wind[num_change] = r;

            }
          while (r0 < geotop::common::Variables::Nr);

          for (i = 1; i < num_change; i++)
            {
              // south wind
              if ((-cos(met->Vdir[snow->change_dir_wind[i]][c] * GTConst::Pi /
                        180.)) < 0)
                {
                  if (par->upwindblowingsnow == 1 && snow->change_dir_wind[i] != 1)
                    snow->Qtrans_y[snow->change_dir_wind[i]][c] = 0.0;
                  snow->Qtrans_y[snow->change_dir_wind[i]][c] = 0.0;
                  for (r = snow->change_dir_wind[i] + 1;
                       r <= snow->change_dir_wind[i + 1]; r++)
                    {
                      if (snow->change_dir_wind[i + 1] == geotop::common::Variables::Nr ||
                          (snow->change_dir_wind[i + 1] !=
                           geotop::common::Variables::Nr &&
                           r < snow->change_dir_wind[i + 1]))
                        {
                          Qup = snow->Qtrans_y[r - 1][c];
                          Qdown = snow->Qtrans_y[r][c];
                          Sup = snow->Qsub_y[r - 1][c];
                          Sdown = snow->Qsub_y[r][c];
                          if (Qdown >= Qup)
                            {
                              snow->Qtrans_y[r][c] = (Qdown + Qup * F1 / dy) / (1. + F1 / dy);
                              snow->Qsub_y[r][c] = (Sdown + Sup * F1 / dy) / (1. + F1 / dy);
                            }
                          else
                            {
                              snow->Qtrans_y[r][c] = (Qdown + Qup * F2 / dy) / (1. + F2 / dy);
                              snow->Qsub_y[r][c] = (Sdown + Sup * F2 / dy) / (1. + F2 / dy);
                            }
                          Qdown = snow->Qtrans_y[r][c];
                          snow->Nabla2_Qtrans[r][c] += (Qup - Qdown) / dy;
                        }
                    }

                  // north wind
                }
              else
                {
                  if (par->upwindblowingsnow == 1 &&
                      snow->change_dir_wind[i + 1] != geotop::common::Variables::Nr)
                    snow->Qtrans_y[snow->change_dir_wind[i + 1] - 1][c] = 0.0;
                  snow->Qtrans_y[snow->change_dir_wind[i + 1] - 1][c] = 0.0;
                  for (r = snow->change_dir_wind[i + 1] - 1;
                       r >= snow->change_dir_wind[i]; r--)
                    {
                      if (snow->change_dir_wind[i + 1] == geotop::common::Variables::Nr ||
                          (snow->change_dir_wind[i + 1] !=
                           geotop::common::Variables::Nr &&
                           r < snow->change_dir_wind[i + 1] - 1))
                        {
                          Qup = snow->Qtrans_y[r + 1][c];
                          Qdown = snow->Qtrans_y[r][c];
                          Sup = snow->Qsub_y[r + 1][c];
                          Sdown = snow->Qsub_y[r][c];
                          if (Qdown >= Qup)
                            {
                              snow->Qtrans_y[r][c] = (Qdown + Qup * F1 / dy) / (1. + F1 / dy);
                              snow->Qsub_y[r][c] = (Sdown + Sup * F1 / dy) / (1. + F1 / dy);
                            }
                          else
                            {
                              snow->Qtrans_y[r][c] = (Qdown + Qup * F2 / dy) / (1. + F2 / dy);
                              snow->Qsub_y[r][c] = (Sdown + Sup * F2 / dy) / (1. + F2 / dy);
                            }
                          Qdown = snow->Qtrans_y[r][c];
                          snow->Nabla2_Qtrans[r][c] += (Qup - Qdown) / dy;
                        }
                    }
                }
            }
        }

      // Adjusting snow init depth in case of steep slope (contribution by Stephan
      // Gruber)
      for (r = 1; r <= geotop::common::Variables::Nr; r++)
        {
          for (c = 1; c <= geotop::common::Variables::Nc; c++)
            {
              if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
                {
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
                      if (snow->Nabla2_Qtrans[r][c] > 0)
                        snow->Nabla2_Qtrans[r][c] *= k_snowred;
                    }
                }
            }
        }
    }
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

static void set_windtrans_snow(double Dt,
                               double t,
                               Snow *snow,
                               Meteo *met,
                               Land *land,
                               Par *par,
                               FILE *f)
{
  long i, l, r, c, ns;
  double Qsub, DW, DWl, Wtrans_tot = 0.0, Wsubl_tot = 0.0, Utot = 0.0;
  short ok = 1;
  (void)ok;  // silence warning unused var

  // density of wind transported snow
  float rho_wind_transported_snow = 400.0;

  // snow compaction for effect of wind
  double CR;
  double overburden;
  float A4 = 7.81E-3;  // m kg^-1
  float D4 = 0.0884;   // m2 N-1

  // update snow depth
  for (r = 1; r <= geotop::common::Variables::Nr; r++)
    {
      for (c = 1; c <= geotop::common::Variables::Nc; c++)
        {
          if ((long)land->LC[r][c] != geotop::input::gDoubleNoValue)
            {
              ns = snow->S->lnum[r][c];

              Qsub = sqrt(pow(snow->Qsub_x[r][c], 2.) + pow(snow->Qsub_y[r][c], 2.));

              Wsubl_tot += Dt * Qsub / (double)par->total_pixel;
              Wtrans_tot += Dt * snow->Nabla2_Qtrans[r][c] / (double)par->total_pixel;
              Utot += met->Vgrid[r][c] / (double)par->total_pixel;

              DW = Dt * (snow->Nabla2_Qtrans[r][c] - Qsub);

              if (ns > 0)  // snow on the soil
                {
                  if (DW < 0)  // snow eroded
                    {
                      i = ns;
                      DWl = 0.0;

                      do
                        {
                          if (i < ns)
                            {
                              ok = 0;
                              if (snow->S->w_ice[i + 1][r][c] < 0)
                                {
                                  DW = snow->S->w_ice[i + 1][r][c];
                                  DWl = snow->S->w_liq[i + 1][r][c];
                                  snow->S->w_ice[i + 1][r][c] = 0.0;
                                  snow->S->w_liq[i + 1][r][c] = 0.0;
                                  snow->S->Dzl[i + 1][r][c] = 0.0;
                                  snow->S->lnum[r][c] -= 1;
                                }
                            }

                          snow->S->Dzl[i][r][c] *=
                            (snow->S->w_ice[i][r][c] + DW) / snow->S->w_ice[i][r][c];
                          // kg/m2
                          snow->S->w_ice[i][r][c] += DW;
                          // kg/m2
                          snow->S->w_liq[i][r][c] += DWl;

                          i--;

                        }
                      while (snow->S->w_ice[i + 1][r][c] < 0 && i > 0);

                      if (i == 0 && snow->S->w_ice[i + 1][r][c] < 0)
                        {
                          // kg/m2
                          snow->S->w_ice[i + 1][r][c] = 0.0;
                          // mm
                          snow->S->Dzl[i + 1][r][c] = 0.0;
                          snow->S->lnum[r][c] = 0;
                        }

                    }  // snow drifted
                  else
                    {
                      i = snow->S->lnum[r][c];
                      snow->S->w_ice[i][r][c] += DW;
                      snow->S->Dzl[snow->S->lnum[r][c]][r][c] +=
                        1.0E+3 * DW / rho_wind_transported_snow;
                    }

                }  // snot not on the soil
              else
                {
                  if (DW > 0)
                    {
                      snow->S->w_ice[1][r][c] += DW;
                      snow->S->Dzl[1][r][c] += 1.0E+3 * DW / rho_wind_transported_snow;
                      snow->S->T[1][r][c] = Fmin(-1., met->Tgrid[r][c]);
                    }
                }

              // from Liston(2007) - wind packing factor
              if (snow->Qtrans[r][c] > 1.E-10)
                {
                  overburden = 0.0;
                  for (l = snow->S->lnum[r][c]; l >= 1; l--)
                    {
                      overburden +=
                        (snow->S->w_ice[l][r][c] + snow->S->w_liq[l][r][c]) / 2.;

                      // compactation at the surface (10%/hour if U8=8m/s U8t=4m/s =>
                      // Qsalt=3.555342e-03 kg/m/s)
                      CR = -A4 * snow->Qsalt[r][c];
                      // decrease due to oberburden
                      CR *= exp(-D4 * GTConst::GRAVITY * overburden);

                      snow->S->Dzl[l][r][c] *= exp(CR * Dt);

                      if (snow->S->w_ice[l][r][c] /
                          (GTConst::rho_w * snow->S->Dzl[l][r][c] * 1.E-3) >
                          par->snow_maxpor)
                        {
                          snow->S->Dzl[l][r][c] = 1.E3 * snow->S->w_ice[l][r][c] /
                                                  (GTConst::rho_w * par->snow_maxpor);
                        }

                      overburden +=
                        (snow->S->w_ice[l][r][c] + snow->S->w_liq[l][r][c]) / 2.;
                    }
                }

              snow_layer_combination(par->alpha_snow, r, c, snow->S, met->Tgrid[r][c],
                                     par->inf_snow_layers, par->max_weq_snow, 1.E10);
            }
        }
    }

  fprintf(f, "BLOWING SNOW: t:%f Wsubl:%e Wtrans:%e U:%f Dt:%f \n", t,
          Wsubl_tot, Wtrans_tot, Utot, Dt);
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

static void print_windtrans_snow(double Dt, Snow *snow, Par *par, Topo *top)
{
  size_t i, r, c;
  double Qsub;

  for (i = 1; i <= par->total_pixel; i++)
    {
      r = top->rc_cont[i][1];
      c = top->rc_cont[i][2];

      Qsub = sqrt(pow(snow->Qsub_x[r][c], 2.) + pow(snow->Qsub_y[r][c], 2.));

      if (par->output_snow[geotop::common::Variables::i_sim] > 0)
        {
          snow->Wtrans_plot[r][c] += Dt * snow->Nabla2_Qtrans[r][c];
          snow->Wsubl_plot[r][c] += Dt * Qsub;
        }

      if (par->Dtplot_point[geotop::common::Variables::i_sim] > 1.E-5 &&
          par->state_pixel == 1 && par->jplot[i] > 0)
        {
          geotop::common::Variables::odp[oblowingsnowtrans][par->jplot[i] - 1] -=
            Dt * snow->Nabla2_Qtrans[r][c];
          geotop::common::Variables::odp[oblowingsnowsubl][par->jplot[i] - 1] +=
            Dt * Qsub;
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
