/* STATEMENT:

   GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
   GEOtop 2.1 release candidate  (release date: 31 december 2016)

   Copyright (c), 2016 - GEOtop Foundation

   This file is part of GEOtop 2.1

   GEOtop 2.1  is a free software and is distributed under GNU General Public
   License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
   even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE

   GEOtop 2.1  is distributed as a free software in the hope to create and
   support a community of developers and users that constructively interact. If
   you just use the code, please give feedback to the authors and the community.
   Any way you use the model, may be the most trivial one, is significantly
   helpful for the future development of the GEOtop model. Any feedback will be
   highly appreciated.

   If you have satisfactorily used the code, please acknowledge the authors.

*/

#include <string>
#include "output.h"
#include "constants.h"
#include "geotop_common.h"
#include "inputKeywords.h"

#include "global_logger.h"

using namespace std;

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_output(Times *times,
                  Water *wat,
                  Channel *cnet,
                  Par *par,
                  Topo *top,
                  Land *land,
                  Soil *sl,
                  Energy *egy,
                  Snow *snow,
                  Glacier *glac,
                  Meteo *met)

{
  /*internal auxiliary variables:*/
  size_t i;
  long j, l, r = 0L, c = 0L, m; /*counters*/
  long n_file; /*number of file of the type "TETAxySSSlZZ"(i.e. number of the
                  basin-time-step)*/
  std::string NNNNN = "NNNNN";
  std::string RRRRR = "RRRRR"; /*TODO: remove this*/
  std::string SSSSS = "SSSSS"; /*TODO: remove this*/
  std::string NNNN = "NNNN";
  std::string rec = "_recNNNN";
  std::string crec = "_crecNNNN";

  string name, temp1, temp2, s1, s2;
  FILE *f = NULL;

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  //  time variables
  time_t stop_time;
  double percent_done = 0.0, remaining_time, total_time;
  short first_column;
  double JD, JDfrom0;
  long day, month, year, hour, minute;

  //  static double Qsub_ch, Qsup_ch;
  // static long isavings;
  static double mass_error_tot;
  static double t_discharge, t_basin, t_point, t_rec;
  double Vchannel, Vsub, Vsup;

  //  other variables
  GeoVector<double> V;

  GeoMatrix<double> M;

  double D, cosslope;

  // initialize static variables
  if (times->time < 1.E-5)
    {
      mass_error_tot = 0.;
      t_discharge = 0.;
      t_basin = 0.;
      t_point = 0.;
      t_rec = 0.;
    }

  write_suffix(SSSSS, geotop::common::Variables::i_sim,
               1);  // TODO: remove this!
  write_suffix(RRRRR, geotop::common::Variables::i_run,
               1);  // TODO: remove this!

  // Time indices
  JDfrom0 = convert_tfromstart_JDfrom0(times->time + par->Dt, par->init_date);
  convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
  convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);

  // DISCHARGE
  //****************************************************************************************************
  //*************************************************************************************
  //**************

  if (par->state_discharge == 1 && par->Dtplot_discharge > 1.E-5 &&
      geotop::common::Variables::files[fQ] != geotop::input::gStringNoValue)
    {
      t_discharge += par->Dt;

      if (fabs(t_discharge - par->Dtplot_discharge) < 1.E-5)
        {
          Vchannel = 0.;
          Vsub = 0.;
          Vsup = 0.;
          for (l = 1; l <= long(par->total_channel); l++)
            {
              r = cnet->r[l];
              c = cnet->c[l];
              Vchannel += 1.E-3 * Fmax(cnet->SS->P[0][l], 0.) /
                          cos(top->slope[r][c] * GTConst::Pi / 180.) *
                          geotop::common::Variables::UV->U[1] * par->w_dx *
                          cnet->length[l];
              Vsub += cnet->Vsub[l];
              Vsup += cnet->Vsup[l];
            }

          if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

          if (par->n_ContRecovery > 0)
            {
              temp1 = geotop::common::Variables::files[fQ] + string(crec);
              name = temp1 + textfile;
            }
          else
            {
              name = geotop::common::Variables::files[fQ] + string(textfile);
            }

          f = fopen(name.c_str(), "a");
          fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)day, (float)month,
                  (float)year, (float)hour, (float)minute);

          fprintf(f, ",%f,%f,%f", (times->time + par->Dt) / GTConst::secinday,
                  JDfrom0, JD);
          fprintf(f, ",%e,%e,%e,%e,%e,%e,%e\n",
                  cnet->Vout / (double)par->Dtplot_discharge,
                  Vsup / (double)par->Dtplot_discharge,
                  Vsub / (double)par->Dtplot_discharge, Vchannel,
                  wat->Voutlandsup / (double)par->Dtplot_discharge,
                  wat->Voutlandsub / (double)par->Dtplot_discharge,
                  wat->Voutbottom / (double)par->Dtplot_discharge);
          fclose(f);

          t_discharge = 0.0;

          cnet->Vsub.resize(cnet->Vsub.size(), 0.);
          cnet->Vsup.resize(cnet->Vsup.size(), 0.);

          cnet->Vout = 0.;
          wat->Voutbottom = 0.;
          wat->Voutlandsub = 0.;
          wat->Voutlandsup = 0.;
        }
    }

  // DATA POINT
  //********************************************************************************************************
  //*********************************************************************************************************

  if (par->Dtplot_point[geotop::common::Variables::i_sim] > 1.E-5)
    {
      t_point += par->Dt;
      geotop::common::Variables::cum_time += par->Dt;

      if (par->state_pixel == 1)
        {
          for (i = 1; i < par->rc.getRows(); i++)
            {
              r = par->rc[i][1];
              c = par->rc[i][2];
              j = top->j_cont[r][c];

              for (long ll = 1; ll <= geotop::common::Variables::Nl; ll++)
                {
                  if (geotop::common::Variables::files[fTzav] !=
                      geotop::input::gStringNoValue ||
                      geotop::common::Variables::files[fTzavwriteend] !=
                      geotop::input::gStringNoValue)
                    sl->Tzavplot[i][ll] +=
                      sl->SS->T[ll][j] *
                      (par->Dt / par->Dtplot_point[geotop::common::Variables::i_sim]);
                  if (geotop::common::Variables::files[fliqzav] !=
                      geotop::input::gStringNoValue ||
                      geotop::common::Variables::files[fliqzavwriteend] !=
                      geotop::input::gStringNoValue)
                    sl->thzavplot[i][ll] +=
                      sl->th[ll][j] *
                      (par->Dt / par->Dtplot_point[geotop::common::Variables::i_sim]);
                  if (geotop::common::Variables::files[ficezav] !=
                      geotop::input::gStringNoValue ||
                      geotop::common::Variables::files[ficezavwriteend] !=
                      geotop::input::gStringNoValue)
                    sl->thizavplot[i][ll] +=
                      sl->SS->thi[ll][j] *
                      (par->Dt / par->Dtplot_point[geotop::common::Variables::i_sim]);
                }

              D = find_activelayerdepth_up(j, sl->type[r][c], sl);
              geotop::common::Variables::odpnt[othawedup][i - 1] +=
                D * (par->Dt / par->Dtplot_point[geotop::common::Variables::i_sim]);

              D = find_activelayerdepth_dw(j, sl->type[r][c], sl);
              geotop::common::Variables::odpnt[othaweddw][i - 1] +=
                D * (par->Dt / par->Dtplot_point[geotop::common::Variables::i_sim]);

              D = find_watertabledepth_up(j, sl->type[r][c], sl);
              if (D < 1.E-5 && sl->SS->P[0][j] > 0)
                D = -sl->SS->P[0][j] / cos(top->slope[r][c] * GTConst::Pi / 180.);
              geotop::common::Variables::odpnt[owtableup][i - 1] +=
                D * (par->Dt / par->Dtplot_point[geotop::common::Variables::i_sim]);

              if (sl->SS->P[0][j] > 0)
                {
                  D = -sl->SS->P[0][i] / cos(top->slope[r][c] * GTConst::Pi / 180.);
                }
              else
                {
                  D = find_watertabledepth_dw(i, sl->type[r][c], sl);
                }
              geotop::common::Variables::odpnt[owtabledw][i - 1] +=
                D * (par->Dt / par->Dtplot_point[geotop::common::Variables::i_sim]);
            }
        }

      //  Print of pixel-output every times->n_pixel time step
      if (fabs(t_point - par->Dtplot_point[geotop::common::Variables::i_sim]) <
          1.E-5)
        {
          if (par->state_pixel == 1)
            {
              if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

              for (i = 1; i < par->rc.getRows(); i++)
                {
                  write_suffix(NNNN, par->IDpoint[i], 0);
                  r = par->rc[i][1];
                  c = par->rc[i][2];
                  j = top->j_cont[r][c];

                  // sy = sl->type->co[r][c];

                  if (par->output_vertical_distances == 1)
                    {
                      cosslope = cos(Fmin(GTConst::max_slope, top->slope[r][c]) *
                                     GTConst::Pi / 180.);
                    }
                  else
                    {
                      cosslope = 1.;
                    }

                  //  soil data
                  for (l = 1; l <= geotop::common::Variables::Nl; l++)
                    {
                      if (geotop::common::Variables::files[fTz] !=
                          geotop::input::gStringNoValue ||
                          geotop::common::Variables::files[fTzwriteend] !=
                          geotop::input::gStringNoValue)
                        sl->Tzplot[i][l] = sl->SS->T[l][j];
                      if (geotop::common::Variables::files[fpsiztot] !=
                          geotop::input::gStringNoValue ||
                          geotop::common::Variables::files[fpsiztotwriteend] !=
                          geotop::input::gStringNoValue)
                        sl->Ptotzplot[i][l] = sl->Ptot[l][j];
                      if (geotop::common::Variables::files[fliqz] !=
                          geotop::input::gStringNoValue ||
                          geotop::common::Variables::files[fliqzwriteend] !=
                          geotop::input::gStringNoValue)
                        sl->thzplot[i][l] = sl->th[l][j];
                      if (geotop::common::Variables::files[ficez] !=
                          geotop::input::gStringNoValue ||
                          geotop::common::Variables::files[ficezwriteend] !=
                          geotop::input::gStringNoValue)
                        sl->thizplot[i][l] = sl->SS->thi[l][j];
                    }
                  for (l = 0; l <= geotop::common::Variables::Nl; l++)
                    {
                      if (geotop::common::Variables::files[fpsiz] !=
                          geotop::input::gStringNoValue ||
                          geotop::common::Variables::files[fpsizwriteend] !=
                          geotop::input::gStringNoValue)
                        sl->Pzplot[i][l] = sl->SS->P[l][j];
                    }

                  //  snow data

                  if (snow->S->lnum[r][c] > 0)
                    {
                      geotop::common::Variables::odpnt[osnowdepth][i - 1] = 0.0;
                      geotop::common::Variables::odpnt[oSWE][i - 1] = 0.0;
                      geotop::common::Variables::odpnt[osnowT][i - 1] = 0.0;
                      for (l = 1; l <= snow->S->lnum[r][c]; l++)
                        {
                          geotop::common::Variables::odpnt[osnowdepth][i - 1] +=
                            snow->S->Dzl[l][r][c];
                          geotop::common::Variables::odpnt[oSWE][i - 1] +=
                            1.0E+3 * (snow->S->w_liq[l][r][c] + snow->S->w_ice[l][r][c]) /
                            GTConst::rho_w;
                          geotop::common::Variables::odpnt[osnowT][i - 1] +=
                            snow->S->T[l][r][c] * snow->S->Dzl[l][r][c];
                        }
                      geotop::common::Variables::odpnt[osnowdens][i - 1] =
                        geotop::common::Variables::odpnt[oSWE][i - 1] * GTConst::rho_w /
                        geotop::common::Variables::odpnt[osnowdepth][i - 1];
                      geotop::common::Variables::odpnt[osnowT][i - 1] /=
                        geotop::common::Variables::odpnt[osnowdepth][i - 1];
                    }
                  else
                    {
                      geotop::common::Variables::odpnt[osnowdepth][i - 1] = 0.0;
                      geotop::common::Variables::odpnt[oSWE][i - 1] = 0.0;
                      geotop::common::Variables::odpnt[osnowdens][i - 1] =
                        geotop::input::gDoubleNoValue;
                      geotop::common::Variables::odpnt[osnowT][i - 1] =
                        geotop::input::gDoubleNoValue;
                    }

                  // glacier data
                  if (par->max_glac_layers > 0)
                    {
                      if (glac->G->lnum[r][c] > 0)
                        {
                          geotop::common::Variables::odpnt[oglacdepth][i - 1] = 0.0;
                          geotop::common::Variables::odpnt[oGWE][i - 1] = 0.0;
                          geotop::common::Variables::odpnt[oglacT][i - 1] = 0.0;
                          for (l = 1; l <= glac->G->lnum[r][c]; l++)
                            {
                              geotop::common::Variables::odpnt[oglacdepth][i - 1] +=
                                glac->G->Dzl[l][r][c];

                              geotop::common::Variables::odpnt[oGWE][i - 1] +=
                                1.0E+3 * (glac->G->w_liq[l][r][c] + glac->G->w_ice[l][r][c]) /
                                GTConst::rho_w;
                              geotop::common::Variables::odpnt[oglacT][i - 1] +=
                                glac->G->T[l][r][c] * glac->G->Dzl[l][r][c];
                            }
                          geotop::common::Variables::odpnt[oglacdens][i - 1] =
                            geotop::common::Variables::odpnt[oGWE][i - 1] * GTConst::rho_w /
                            geotop::common::Variables::odpnt[oglacdepth][i - 1];
                          geotop::common::Variables::odpnt[oglacT][i - 1] /=
                            geotop::common::Variables::odpnt[oglacdepth][i - 1];
                        }
                      else
                        {
                          geotop::common::Variables::odpnt[oglacdepth][i - 1] = 0.0;
                          geotop::common::Variables::odpnt[oGWE][i - 1] = 0.0;
                          geotop::common::Variables::odpnt[oglacdens][i - 1] =
                            geotop::input::gDoubleNoValue;
                          geotop::common::Variables::odpnt[oglacT][i - 1] =
                            geotop::input::gDoubleNoValue;
                        }
                    }

                  // Point data

                  if (geotop::common::Variables::files[fpoint] !=
                      geotop::input::gStringNoValue)
                    {
                      temp1 = geotop::common::Variables::files[fpoint] + string(NNNN);

                      if (par->n_ContRecovery > 0)
                        {
                          temp2 = temp1 + crec;
                          name = temp2 + textfile;
                        }
                      else
                        {
                          name = temp1 + textfile;
                        }

                      f = fopen(name.c_str(), "a");
                      first_column = 1;
                      for (j = 0; j < geotop::common::Variables::nopnt; j++)
                        {
                          if (first_column == 0)
                            {
                              fprintf(f, ",");
                            }
                          else
                            {
                              first_column = 0;
                            }
                          if (geotop::common::Variables::opnt[j] >= 0)
                            {
                              if (geotop::common::Variables::opnt[j] == odate12)
                                {
                                  fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)day,
                                          (float)month, (float)year, (float)hour,
                                          (float)minute);
                                }
                              else if (geotop::common::Variables::opnt[j] == oJDfrom0)
                                {
                                  fprintf(f, "%f", JDfrom0);
                                }
                              else if (geotop::common::Variables::opnt[j] ==
                                       odaysfromstart)
                                {
                                  fprintf(f, "%f", JDfrom0 - par->init_date);
                                }
                              else if (geotop::common::Variables::opnt[j] == operiod)
                                {
                                  fprintf(f, "%ld", geotop::common::Variables::i_sim);
                                }
                              else if (geotop::common::Variables::opnt[j] == opoint)
                                {
                                  fprintf(f, "%ld", par->IDpoint[i]);
                                }
                              else
                                {
                                  fprintf(f, "%f",
                                          geotop::common::Variables::odpnt
                                          [geotop::common::Variables::opnt[j]][i - 1]);
                                }
                            }
                          else
                            {
                              fprintf(f, "%ld", (long)geotop::input::gDoubleNoValue);
                            }
                        }
                      fprintf(f, "\n");
                      fclose(f);
                    }

                  if (geotop::common::Variables::files[fpointwriteend] !=
                      geotop::input::gStringNoValue)
                    {
                      first_column = 1;
                      for (j = 0; j < geotop::common::Variables::nopnt; j++)
                        {
                          if (first_column == 0)
                            {
                              fprintf(geotop::common::Variables::ffpoint, ",");
                            }
                          else
                            {
                              first_column = 0;
                            }
                          if (geotop::common::Variables::opnt[j] >= 0)
                            {
                              if (geotop::common::Variables::opnt[j] == odate12)
                                {
                                  fprintf(geotop::common::Variables::ffpoint,
                                          "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)day,
                                          (float)month, (float)year, (float)hour,
                                          (float)minute);
                                }
                              else if (geotop::common::Variables::opnt[j] == oJDfrom0)
                                {
                                  fprintf(geotop::common::Variables::ffpoint, "%f", JDfrom0);
                                }
                              else if (geotop::common::Variables::opnt[j] ==
                                       odaysfromstart)
                                {
                                  fprintf(geotop::common::Variables::ffpoint, "%f",
                                          JDfrom0 - par->init_date);
                                }
                              else if (geotop::common::Variables::opnt[j] == operiod)
                                {
                                  fprintf(geotop::common::Variables::ffpoint, "%ld",
                                          geotop::common::Variables::i_sim);
                                }
                              else if (geotop::common::Variables::opnt[j] == opoint)
                                {
                                  fprintf(geotop::common::Variables::ffpoint, "%ld",
                                          par->IDpoint[i]);
                                }
                              else
                                {
                                  fprintf(geotop::common::Variables::ffpoint, "%12g",
                                          geotop::common::Variables::odpnt
                                          [geotop::common::Variables::opnt[j]][i - 1]);
                                  fprintf(geotop::common::Variables::ffpoint, "%f",
                                          geotop::common::Variables::odpnt
                                          [geotop::common::Variables::opnt[j]][i - 1]);
                                }
                            }
                          else
                            {
                              fprintf(geotop::common::Variables::ffpoint, "%g",
                                      geotop::input::gDoubleNoValue);
                            }
                        }
                      fprintf(geotop::common::Variables::ffpoint, "\n");
                    }

                  // Glacier
                  if (par->max_glac_layers > 0)
                    {
                      if (geotop::common::Variables::files[fglz] !=
                          geotop::input::gStringNoValue)
                        {
                          temp1 = geotop::common::Variables::files[fglz] + string(NNNN);

                          if (par->n_ContRecovery > 0)
                            {
                              temp2 = temp1 + crec;
                              name = temp2 + textfile;

                            }
                          else
                            {
                              name = temp1 + textfile;
                            }

                          f = fopen(name.c_str(), "a");

                          if ((long)par->glac_plot_depths[1] !=
                              geotop::input::gDoubleNoValue)
                            {
                              m = par->glac_plot_depths.size();
                            }
                          else
                            {
                              m = par->max_glac_layers;
                            }

                          first_column = 1;
                          for (j = 0; j < geotop::common::Variables::noglc; j++)
                            {
                              if (first_column == 0)
                                {
                                  fprintf(f, ",");
                                }
                              else
                                {
                                  first_column = 0;
                                }
                              if (geotop::common::Variables::oglc[j] >= 0)
                                {
                                  if (geotop::common::Variables::oglc[j] == 0)
                                    {
                                      fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)day,
                                              (float)month, (float)year, (float)hour,
                                              (float)minute);
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 1)
                                    {
                                      fprintf(f, "%f", JDfrom0);
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 2)
                                    {
                                      fprintf(f, "%f", JDfrom0 - par->init_date);
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 3)
                                    {
                                      fprintf(f, "%ld", geotop::common::Variables::i_sim);
                                    }
                                  else if (geotop::common::Variables::oglc[j] ==
                                           4)  // Needed to balance the braces
                                    {
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 5)
                                    {
                                      fprintf(f, "%ld", par->IDpoint[i]);
                                    }
                                  else if (geotop::common::Variables::oglc[j] <= 5 + 1 * m)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 0 * m;
                                      if ((long)par->glac_plot_depths[1] !=
                                          geotop::input::gDoubleNoValue)
                                        {
                                          fprintf(
                                            f, "%f",
                                            interpolate_snow(
                                              r, c, par->glac_plot_depths[l] * cosslope,
                                              glac->G->lnum[r][c], glac->G->Dzl, glac->G->T, 0.));
                                        }
                                      else
                                        {
                                          fprintf(f, "%f", glac->G->T[l][r][c]);
                                        }
                                    }
                                  else if (geotop::common::Variables::oglc[j] <= 5 + 2 * m)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 1 * m;
                                      if ((long)par->glac_plot_depths[1] !=
                                          geotop::input::gDoubleNoValue)
                                        {
                                          fprintf(f, "%f",
                                                  interpolate_snow(
                                                    r, c, par->glac_plot_depths[l] * cosslope,
                                                    glac->G->lnum[r][c], glac->G->Dzl,
                                                    glac->G->w_ice, 0.));
                                        }
                                      else
                                        {
                                          fprintf(f, "%f", glac->G->w_ice[l][r][c]);
                                        }
                                    }
                                  else if (geotop::common::Variables::oglc[j] <= 5 + 3 * m)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 2 * m;
                                      if ((long)par->glac_plot_depths[1] !=
                                          geotop::input::gDoubleNoValue)
                                        {
                                          fprintf(f, "%f",
                                                  interpolate_snow(
                                                    r, c, par->glac_plot_depths[l] * cosslope,
                                                    glac->G->lnum[r][c], glac->G->Dzl,
                                                    glac->G->w_liq, 0.));

                                        }
                                      else
                                        {
                                          fprintf(f, "%f", glac->G->w_liq[l][r][c]);
                                        }
                                    }
                                  else if (geotop::common::Variables::oglc[j] <=
                                           5 + 3 * m + par->max_glac_layers)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 3 * m;
                                      fprintf(f, "%f", glac->G->Dzl[l][r][c]);
                                    }
                                }
                              else
                                {
                                  fprintf(f, "%f", geotop::input::gDoubleNoValue);
                                }
                            }
                          fprintf(f, "\n");
                          fclose(f);
                        }

                      if (geotop::common::Variables::files[fglzwriteend] !=
                          geotop::input::gStringNoValue)
                        {
                          if ((long)par->glac_plot_depths[1] !=
                              geotop::input::gDoubleNoValue)
                            {
                              m = par->glac_plot_depths.size();
                            }
                          else
                            {
                              m = par->max_glac_layers;
                            }
                          first_column = 1;
                          for (j = 0; j < geotop::common::Variables::noglc; j++)
                            {
                              if (first_column == 0)
                                {
                                  fprintf(f, ",");
                                }
                              else
                                {
                                  first_column = 0;
                                }
                              if (geotop::common::Variables::oglc[j] >= 0)
                                {
                                  if (geotop::common::Variables::oglc[j] == 0)
                                    {
                                      fprintf(geotop::common::Variables::ffglac,
                                              "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)day,
                                              (float)month, (float)year, (float)hour,
                                              (float)minute);
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 1)
                                    {
                                      fprintf(geotop::common::Variables::ffglac, "%f", JDfrom0);
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 2)
                                    {
                                      fprintf(geotop::common::Variables::ffglac, "%f",
                                              JDfrom0 - par->init_date);
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 3)
                                    {
                                      fprintf(geotop::common::Variables::ffglac, "%ld",
                                              geotop::common::Variables::i_sim);
                                    }
                                  else if (geotop::common::Variables::oglc[j] ==
                                           4)  // Needed to balance the braces
                                    {
                                    }
                                  else if (geotop::common::Variables::oglc[j] == 5)
                                    {
                                      fprintf(geotop::common::Variables::ffglac, "%ld",
                                              par->IDpoint[i]);
                                    }
                                  else if (geotop::common::Variables::oglc[j] <= 5 + 1 * m)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 0 * m;
                                      if ((long)par->glac_plot_depths[1] !=
                                          geotop::input::gDoubleNoValue)
                                        {
                                          fprintf(
                                            geotop::common::Variables::ffglac, "%f",
                                            interpolate_snow(
                                              r, c, par->glac_plot_depths[l] * cosslope,
                                              glac->G->lnum[r][c], glac->G->Dzl, glac->G->T, 0.));

                                        }
                                      else
                                        {
                                          fprintf(geotop::common::Variables::ffglac, "%f",
                                                  glac->G->T[l][r][c]);
                                        }
                                    }
                                  else if (geotop::common::Variables::oglc[j] <= 5 + 2 * m)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 1 * m;
                                      if ((long)par->glac_plot_depths[1] !=
                                          geotop::input::gDoubleNoValue)
                                        {
                                          fprintf(geotop::common::Variables::ffglac, "%f",
                                                  interpolate_snow(
                                                    r, c, par->glac_plot_depths[l] * cosslope,
                                                    glac->G->lnum[r][c], glac->G->Dzl,
                                                    glac->G->w_ice, 0.));
                                        }
                                      else
                                        {
                                          fprintf(geotop::common::Variables::ffglac, "%f",
                                                  glac->G->w_ice[l][r][c]);
                                        }
                                    }
                                  else if (geotop::common::Variables::oglc[j] <= 5 + 3 * m)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 2 * m;
                                      if ((long)par->glac_plot_depths[1] !=
                                          geotop::input::gDoubleNoValue)
                                        {
                                          fprintf(geotop::common::Variables::ffglac, "%f",
                                                  interpolate_snow(
                                                    r, c, par->glac_plot_depths[l] * cosslope,
                                                    glac->G->lnum[r][c], glac->G->Dzl,
                                                    glac->G->w_liq, 0.));
                                        }
                                      else
                                        {
                                          fprintf(geotop::common::Variables::ffglac, "%f",
                                                  glac->G->w_liq[l][r][c]);
                                        }
                                    }
                                  else if (geotop::common::Variables::oglc[j] <=
                                           5 + 3 * m + par->max_glac_layers)
                                    {
                                      l = geotop::common::Variables::oglc[j] - 5 - 3 * m;

                                      fprintf(geotop::common::Variables::ffglac, "%f",
                                              glac->G->Dzl[l][r][c]);
                                    }
                                }
                              else
                                {
                                  fprintf(geotop::common::Variables::ffglac, "%f",
                                          geotop::input::gDoubleNoValue);
                                }
                            }
                          fprintf(geotop::common::Variables::ffglac, "\n");
                        }
                    }

                  // sl output

                  write_soil_output(i, par->IDpoint[i], par->init_date, JDfrom0, JD,
                                    day, month, year, hour, minute,
                                    par->soil_plot_depths, sl, par, GTConst::PsiMin,
                                    cosslope);

                  // snow output : this to enable at later stage (SC: 26.12.2013) //
                  // done by Matteo and Leonardo, 13.08.2015
                  write_snow_output(i, par->IDpoint[i], r, c, par->init_date, JDfrom0,
                                    day, month, year, hour, minute,
                                    par->snow_plot_depths, snow->S, par, cosslope);

                  //  initialize
                  for (j = 0; j < otot; j++)
                    {
                      geotop::common::Variables::odpnt[j][i - 1] = 0.0;
                    }
                }
            }

          percent_done = 100. * geotop::common::Variables::cum_time /
                         geotop::common::Variables::max_time;

          time(&stop_time);
          geotop::common::Variables::elapsed_time =
            difftime(stop_time, geotop::common::Variables::start_time) +
            geotop::common::Variables::elapsed_time_start;

          if (percent_done > 1.0e-6)
            {
              total_time =
                geotop::common::Variables::elapsed_time * 100.0 / percent_done;
            }
          else
            {
              total_time = 1.E5;
            }

          remaining_time = (total_time - geotop::common::Variables::elapsed_time);

          // logging on file..
          lg->logsf(
            geotop::logger::NOTICE,
            "%ld/%ld/%ld %ld:%02.0f %.2f%% - Time elapsed (h:m:s) "
            "%2.0f:%02.0f:%02.0f Time remaining (h:m) %2.0f:%02.0f",
            day, month, year, hour, (float)minute, percent_done,
            floor(geotop::common::Variables::elapsed_time / 3600.0),
            floor(((geotop::common::Variables::elapsed_time / 3600) -
                   floor(geotop::common::Variables::elapsed_time / 3600.0)) *
                  60.),
            floor((((geotop::common::Variables::elapsed_time / 3600) -
                    floor(geotop::common::Variables::elapsed_time / 3600.0)) *
                   60. -
                   floor(((geotop::common::Variables::elapsed_time / 3600) -
                          floor(geotop::common::Variables::elapsed_time / 3600.0)) *
                         60.)) *
                  60.),
            floor(remaining_time / 3600.0),
            floor(((remaining_time / 3600) - floor(remaining_time / 3600.0)) *
                  60.));

          t_point = 0.0;
        }
    }

  // BASIN DATA
  //****************************************************************************************************************
  //****************************************************************************************************************

  if (par->Dtplot_basin[geotop::common::Variables::i_sim] > 1.E-5 &&
      par->state_basin == 1)
    {
      t_basin += par->Dt;

      if (fabs(t_basin - par->Dtplot_basin[geotop::common::Variables::i_sim]) <
          1.E-5)
        {
          if (geotop::common::Variables::files[fbas] !=
              geotop::input::gStringNoValue)
            {
              if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

              if (par->n_ContRecovery > 0)
                {
                  temp1 = geotop::common::Variables::files[fbas] + string(crec);
                  name = temp1 + textfile;
                }
              else
                {
                  name = geotop::common::Variables::files[fbas] + string(textfile);
                }

              f = fopen(name.c_str(), "a");
              first_column = 1;
              for (j = 0; j < geotop::common::Variables::nobsn; j++)
                {
                  if (first_column == 0)
                    {
                      fprintf(f, ",");
                    }
                  else
                    {
                      first_column = 0;
                    }
                  if (geotop::common::Variables::obsn[j] >= 0)
                    {
                      if (geotop::common::Variables::obsn[j] == oodate12)
                        {
                          fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)day,
                                  (float)month, (float)year, (float)hour, (float)minute);
                        }
                      else if (geotop::common::Variables::obsn[j] == ooJDfrom0)
                        {
                          fprintf(f, "%f", JDfrom0);

                        }
                      else if (geotop::common::Variables::obsn[j] == oodaysfromstart)
                        {
                          fprintf(f, "%f", JDfrom0 - par->init_date);

                        }
                      else
                        {
                          fprintf(f, "%f",
                                  geotop::common::Variables::odbsn
                                  [geotop::common::Variables::obsn[j]]);
                        }
                    }
                  else
                    {
                      fprintf(f, "%f", geotop::input::gDoubleNoValue);
                    }
                }
              fprintf(f, "\n");
              fclose(f);
            }

          if (geotop::common::Variables::files[fbaswriteend] !=
              geotop::input::gStringNoValue)
            {
              first_column = 1;
              for (j = 0; j < geotop::common::Variables::nobsn; j++)
                {
                  if (first_column == 0)
                    {
                      fprintf(geotop::common::Variables::ffbas, ",");
                    }
                  else
                    {
                      first_column = 0;
                    }
                  if (geotop::common::Variables::obsn[j] >= 0)
                    {
                      if (geotop::common::Variables::obsn[j] == oodate12)
                        {
                          fprintf(geotop::common::Variables::ffbas,
                                  "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)day,
                                  (float)month, (float)year, (float)hour, (float)minute);
                        }
                      else if (geotop::common::Variables::obsn[j] == ooJDfrom0)
                        {
                          fprintf(geotop::common::Variables::ffbas, "%f", JDfrom0);

                        }
                      else if (geotop::common::Variables::obsn[j] == oodaysfromstart)
                        {
                          fprintf(geotop::common::Variables::ffbas, "%f",
                                  JDfrom0 - par->init_date);
                        }
                      else
                        {
                          fprintf(geotop::common::Variables::ffbas, "%f",
                                  geotop::common::Variables::odbsn
                                  [geotop::common::Variables::obsn[j]]);
                        }
                    }
                  else
                    {
                      fprintf(geotop::common::Variables::ffbas, "%f",
                              geotop::input::gDoubleNoValue);
                    }
                }
              fprintf(geotop::common::Variables::ffbas, "\n");
            }

          mass_error_tot += geotop::common::Variables::odbsn[oomasserror];

          // logging facility

          lg->logsf(
            geotop::logger::NOTICE,
            "%ld/%ld/%ld %ld:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed!",
            day, month, year, hour, (float)minute, JD,
            (long)(floor(times->time / 86400)) + 1, percent_done);

          lg->logsf(
            geotop::logger::NOTICE,
            " t_meteo:%6.2f s,t_energy:%6.2f s, t_blowingsnow:%6.2f s, "
            "t_water:%6.2f s, t_sub:%6.2f s, t_sup:%6.2f s, t_out:%6.2f s",
            geotop::common::Variables::t_meteo, geotop::common::Variables::t_energy,
            geotop::common::Variables::t_blowingsnow,
            geotop::common::Variables::t_water, geotop::common::Variables::t_sub,
            geotop::common::Variables::t_sup, geotop::common::Variables::t_out);
          lg->logsf(geotop::logger::NOTICE,
                    "SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 "
                    "Pvert=%6.2f mm Prain=%6.2f mm  Psnow=%6.2f mm",
                    geotop::common::Variables::odbsn[ooSW],
                    geotop::common::Variables::odbsn[ooLW],
                    geotop::common::Variables::odbsn[ooH],
                    geotop::common::Variables::odbsn[ooLE],
                    geotop::common::Variables::odbsn[oopnet],
                    geotop::common::Variables::odbsn[oorainover],
                    geotop::common::Variables::odbsn[oosnowover]);
          lg->logsf(geotop::logger::NOTICE, "Max Error Richards=%e mm/h",
                    geotop::common::Variables::odbsn[oomasserror] * 3600.0 /
                    par->Dtplot_basin[geotop::common::Variables::i_sim]);
          lg->logsf(geotop::logger::NOTICE,
                    "Tot Error Richards=%e mm Mean Time Step=%f s", mass_error_tot,
                    geotop::common::Variables::odbsn[ootimestep]);

          //
          for (j = 0; j < ootot; j++)
            {
              geotop::common::Variables::odbsn[j] = 0.0;
            }
          t_basin = 0.0;
        }
    }

  // DISTRIBUTED OUTPUTS
  //****************************************************************************************************************
  //****************************************************************************************************************
  // averaging properties
  if (par->output_meteo[geotop::common::Variables::i_sim] > 0)
    {
      if (geotop::common::Variables::files[fTa] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              met->Tamean[i] +=
                met->Tgrid[top->rc_cont[i][1]][top->rc_cont[i][2]] /
                ((par->output_meteo[geotop::common::Variables::i_sim] * 3600.0) /
                 (par->Dt));
            }
        }
      if (geotop::common::Variables::files[fwspd] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              met->Vspdmean[i] +=
                met->Vgrid[top->rc_cont[i][1]][top->rc_cont[i][2]] /
                ((par->output_meteo[geotop::common::Variables::i_sim] * 3600.0) /
                 (par->Dt));
            }
        }
      if (geotop::common::Variables::files[fwdir] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              met->Vdirmean[i] +=
                met->Vdir[top->rc_cont[i][1]][top->rc_cont[i][2]] /
                ((par->output_meteo[geotop::common::Variables::i_sim] * 3600.0) /
                 (par->Dt));
            }
        }
      if (geotop::common::Variables::files[frh] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              met->RHmean[i] +=
                met->RHgrid[top->rc_cont[i][1]][top->rc_cont[i][2]] /
                ((par->output_meteo[geotop::common::Variables::i_sim] * 3600.0) /
                 (par->Dt));
            }
        }
    }

  if (par->output_soil[geotop::common::Variables::i_sim] > 0)
    {
      if (geotop::common::Variables::files[fTav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fTavsup] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  sl->T_av_tensor[l][i] +=
                    sl->SS->T[l][i] /
                    ((par->output_soil[geotop::common::Variables::i_sim] * 3600.0) /
                     (par->Dt));
                }
            }
        }
      if (geotop::common::Variables::files[fliqav] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  sl->thw_av_tensor[l][i] +=
                    sl->th[l][i] /
                    ((par->output_soil[geotop::common::Variables::i_sim] * 3600.0) /
                     (par->Dt));
                }
            }
        }
      if (geotop::common::Variables::files[ficeav] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  sl->thi_av_tensor[l][i] +=
                    sl->SS->thi[l][i] /
                    ((par->output_soil[geotop::common::Variables::i_sim] * 3600.0) /
                     (par->Dt));
                }
            }
        }
    }

  V.resize(par->total_pixel + 1, geotop::input::gDoubleNoValue);

  //  soil properties
  if (par->output_soil[geotop::common::Variables::i_sim] > 0 &&
      fmod(times->time + par->Dt + par->delay_day_recover * 86400.,
           par->output_soil[geotop::common::Variables::i_sim] * 3600.0) <
      1.E-5)
    {
      n_file =
        (long)((times->time + par->Dt + par->delay_day_recover * 86400.) /
               (par->output_soil[geotop::common::Variables::i_sim] * 3600.0));

      write_suffix(NNNNN, n_file, 1);
      s1 = NNNNN + string("");
      s2 = s1 + "";

      // theta liq tensor
      if (geotop::common::Variables::files[fliq] !=
          geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[fliq], 0, par->format_out,
                sl->th, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa,
                top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[fliq], 0, par->format_out,
                sl->th, geotop::common::Variables::UV, geotop::input::gDoubleNoValue,
                top->j_cont, geotop::common::Variables::Nr,
                geotop::common::Variables::Nc);
            }
        }

      // theta liq surface
      if (geotop::common::Variables::files[fliqsup] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              V[i] = sl->th[1][i];
            }

          temp1 = geotop::common::Variables::files[fliqsup] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      // write thw_av tensor
      if (geotop::common::Variables::files[fliqav] !=
          geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[fliqav], 0, par->format_out,
                sl->thw_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont,
                sl->pa, top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[fliqav], 0, par->format_out,
                sl->thw_av_tensor, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }
        }

      // initialize thw_av_tensor
      if (geotop::common::Variables::files[fliqav] !=
          geotop::input::gStringNoValue)
        sl->thw_av_tensor.resize(sl->thw_av_tensor.getRows(),
                                 sl->thw_av_tensor.getCols(), 0.0);

      // write T tensor
      if (geotop::common::Variables::files[fT] != geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[fT], 0, par->format_out,
                sl->SS->T, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa,
                top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[fT], 0, par->format_out,
                sl->SS->T, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }
        }

      // theta T surface
      if (geotop::common::Variables::files[fTsup] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              V[i] = sl->SS->T[1][i];
            }

          temp1 = geotop::common::Variables::files[fTsup] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      // write Tav tensor
      if (geotop::common::Variables::files[fTav] !=
          geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[fTav], 0, par->format_out,
                sl->T_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont,
                sl->pa, top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[fTav], 0, par->format_out,
                sl->T_av_tensor, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }
        }

      // theta Tav surface
      if (geotop::common::Variables::files[fTavsup] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              V[i] = sl->T_av_tensor[1][i];
            }

          temp1 = geotop::common::Variables::files[fTavsup] + s2;

          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      // initialize T_av_tensor
      if (geotop::common::Variables::files[fTav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fTavsup] !=
          geotop::input::gStringNoValue)
        sl->T_av_tensor.resize(sl->T_av_tensor.getRows(),
                               sl->T_av_tensor.getCols(), 0.0);

      // theta_ice tensor
      if (geotop::common::Variables::files[fice] !=
          geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[fice], 0, par->format_out,
                sl->SS->thi, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa,
                top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[fice], 0, par->format_out,
                sl->SS->thi, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }
        }

      // theta_ice surface
      if (geotop::common::Variables::files[ficesup] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              V[i] = sl->SS->thi[1][i];
            }

          temp1 = geotop::common::Variables::files[ficesup] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      // write thi_av tensor
      if (geotop::common::Variables::files[ficeav] !=
          geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[ficeav], 0, par->format_out,
                sl->thi_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont,
                sl->pa, top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[ficeav], 0, par->format_out,
                sl->thi_av_tensor, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }
        }

      // initialize thi_av_tensor
      if (geotop::common::Variables::files[ficeav] !=
          geotop::input::gStringNoValue)
        sl->thi_av_tensor.resize(sl->thi_av_tensor.getRows(),
                                 sl->thi_av_tensor.getCols(), 0.0);

      // write psi tensors
      if (geotop::common::Variables::files[fpsitot] !=
          geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[fpsitot], 0, par->format_out,
                sl->Ptot, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa,
                top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[fpsitot], 0, par->format_out,
                sl->Ptot, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }
        }

      if (geotop::common::Variables::files[fpsiliq] !=
          geotop::input::gStringNoValue)
        {
          if ((long)par->soil_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              write_tensorseries_soil(
                1, s2, geotop::common::Variables::files[fpsiliq], 0, par->format_out,
                sl->SS->P, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa,
                top->slope, par->output_vertical_distances);
            }
          else
            {
              write_tensorseries3_vector(
                s2, geotop::common::Variables::files[fpsiliq], 0, par->format_out,
                sl->SS->P, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }
        }

      // calculate saturation front depth
      if (geotop::common::Variables::files[fwtable_up] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              // normal
              V[i] = find_watertabledepth_up(i, sl->type[r][c], sl);
              if (V[i] < 1.E-5 && sl->SS->P[0][i] > 0)
                V[i] = -sl->SS->P[0][i] / cos(top->slope[r][c] * GTConst::Pi / 180.);
            }
          temp1 = geotop::common::Variables::files[fwtable_up] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      if (geotop::common::Variables::files[fwtable_dw] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              if (sl->SS->P[0][i] > 0)
                {
                  V[i] = -sl->SS->P[0][i] / cos(top->slope[r][c] * GTConst::Pi / 180.);
                }
              else
                {
                  // normal
                  V[i] = find_watertabledepth_dw(i, sl->type[r][c], sl);
                }
            }
          temp1 = geotop::common::Variables::files[fwtable_dw] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      // calculate active layer depth
      if (geotop::common::Variables::files[fthawed_up] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              // normal
              V[i] = find_activelayerdepth_up(i, sl->type[r][c], sl);
            }
          temp1 = geotop::common::Variables::files[fthawed_up] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      if (geotop::common::Variables::files[fthawed_dw] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              // normal
              V[i] = find_activelayerdepth_dw(i, sl->type[r][c], sl);
            }
          temp1 = geotop::common::Variables::files[fthawed_dw] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      // WATER OVER THE SURFACE
      if (geotop::common::Variables::files[fhsupland] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              V[i] =
                Fmax(0, sl->SS->P[0][i]) / cos(top->slope[r][c] * GTConst::Pi / 180.);
            }

          temp1 = geotop::common::Variables::files[fhsupland] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      if (geotop::common::Variables::files[fhsupch] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              if (cnet->ch[r][c] != 0)
                {
                  V[i] = cnet->SS->P[0][cnet->ch[r][c]] /
                         cos(top->slope[r][c] * GTConst::Pi / 180.);
                }
              else
                {
                  V[i] = geotop::input::gDoubleNoValue;
                }
            }

          temp1 = geotop::common::Variables::files[fhsupch] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }
      // TODO: mattiu
      if (geotop::common::Variables::files[fpnet] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fpnet] + s2;
          write_map_vector(
            temp1, 0, par->format_out, sl->Pnetcum, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          sl->Pnetcum.resize(sl->Pnetcum.size(), 0.0);
        }
      if (geotop::common::Variables::files[fevap] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fevap] + s2;
          write_map_vector(
            temp1, 0, par->format_out, sl->ETcum, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          sl->ETcum.resize(sl->ETcum.size(), 0.0);
        }  // end mattiu
    }

  // snow properties
  if (par->output_snow[geotop::common::Variables::i_sim] > 0 &&
      fmod(times->time + par->Dt + par->delay_day_recover * 86400.,
           par->output_snow[geotop::common::Variables::i_sim] * 3600.0) <
      1.E-5)
    {
      n_file =
        (long)((times->time + par->Dt + par->delay_day_recover * 86400.) /
               (par->output_snow[geotop::common::Variables::i_sim] * 3600.0));
      write_suffix(NNNNN, n_file, 1);
      s1 = NNNNN + string("");
      s2 = s1 + "";

      if (geotop::common::Variables::files[fsnowdepth] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              V[i] = 0.;
              for (l = 1; l <= snow->S->lnum[r][c]; l++)
                {
                  V[i] += snow->S->Dzl[l][r][c];
                }
            }

          temp1 = geotop::common::Variables::files[fsnowdepth] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      if (geotop::common::Variables::files[fHN] !=
          geotop::input::gStringNoValue)  // TODO mattiu
        {
          temp1 = geotop::common::Variables::files[fHN] + s2;
          write_map_vector(
            temp1, 0, par->format_out, snow->HNcum, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          snow->HNcum.resize(snow->HNcum.size(), 0.0);

        }  // end mattiu

      if (geotop::common::Variables::files[fsnowmelt] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fsnowmelt] + s2;
          write_map_vector(
            temp1, 0, par->format_out, snow->MELTED, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          snow->MELTED.resize(snow->MELTED.size(), 0.0);
        }

      if (geotop::common::Variables::files[fsnowsubl] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fsnowsubl] + s2;
          write_map_vector(
            temp1, 0, par->format_out, snow->SUBL, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          snow->SUBL.resize(snow->SUBL.size(), 0.0);
        }

      if (geotop::common::Variables::files[fswe] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              V[i] = 0.;
              for (l = 1; l <= snow->S->lnum[r][c]; l++)
                {
                  V[i] += (snow->S->w_liq[l][r][c] + snow->S->w_ice[l][r][c]);
                }
            }
          temp1 = geotop::common::Variables::files[fswe] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);

          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              V[i] = 0.;
              D = 0.;
              for (l = 1; l <= snow->S->lnum[r][c]; l++)
                {
                  V[i] += (snow->S->w_liq[l][r][c] + snow->S->w_ice[l][r][c]);
                  D += snow->S->Dzl[l][r][c];
                }
              // temporary fix to avoid division by zero: to be better understood...
              if (D != 0.) { V[i] /= (1.E-3 * D); }
            }
          temp1 = geotop::common::Variables::files[fswe] + string("DENSITY");
          temp2 = temp1 + s2;
          write_map_vector(
            temp2, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);

          if (par->blowing_snow == 1)
            {
              temp1 = geotop::common::Variables::files[fswe] + string("WindTrans");
              temp2 = temp1 + s2;
              write_map(temp2, 0, par->format_out, snow->Wtrans_plot,
                        geotop::common::Variables::UV, geotop::input::gDoubleNoValue);

              temp1 = geotop::common::Variables::files[fswe] + string("WindSubl");
              temp2 = temp1 + s2;
              write_map(temp2, 0, par->format_out, snow->Wsubl_plot,
                        geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
              initmatrix(0.0, snow->Wsubl_plot, land->LC,
                         geotop::input::gDoubleNoValue);
            }
        }

      if (geotop::common::Variables::files[fsndur] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fsndur] + s2;
          write_map_vector(
            temp1, 0, par->format_out, snow->t_snow, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          // snow->t_snow.resize(snow->t_snow.size(), 0.0);
        }
    }

  // glacier properties
  if (par->max_glac_layers > 0 &&
      par->output_glac[geotop::common::Variables::i_sim] > 0 &&
      fmod(times->time + par->Dt + par->delay_day_recover * 86400.,
           par->output_glac[geotop::common::Variables::i_sim] * 3600.0) <
      1.E-5)
    {
      n_file =
        (long)((times->time + par->Dt + par->delay_day_recover * 86400.) /
               (par->output_glac[geotop::common::Variables::i_sim] * 3600.0));
      write_suffix(NNNNN, n_file, 1);
      s1 = NNNNN + string("");
      s2 = s1 + "";

      if (geotop::common::Variables::files[fglacdepth] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              V[i] = 0.;
              for (l = 1; l <= glac->G->lnum[r][c]; l++)
                {
                  V[i] += glac->G->Dzl[l][r][c];  // LP like fsnowdepth
                  // V[i] += (glac->G->w_liq[l][r][c] + glac->G->w_ice[l][r][c]);
                }
            }
          temp1 = geotop::common::Variables::files[fglacdepth] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }

      if (geotop::common::Variables::files[fglacmelt] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fglacmelt] + s2;
          write_map_vector(
            temp1, 0, par->format_out, glac->MELTED, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          glac->MELTED.resize(glac->MELTED.size(), 0.);
        }

      if (geotop::common::Variables::files[fglacsubl] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fglacsubl] + s2;
          write_map_vector(
            temp1, 0, par->format_out, glac->SUBL, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          glac->SUBL.resize(glac->SUBL.size(), 0.);
        }

      if (geotop::common::Variables::files[fgwe] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              r = top->rc_cont[i][1];
              c = top->rc_cont[i][2];
              V[i] = 0.;
              for (l = 1; l <= glac->G->lnum[r][c]; l++)
                {
                  V[i] += (glac->G->w_liq[l][r][c] + glac->G->w_ice[l][r][c]);
                }
            }
          temp1 = geotop::common::Variables::files[fgwe] + s2;
          write_map_vector(
            temp1, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }
    }

  // SURFACE ENERGY BALANCE

  // RADIATION
  if (par->output_surfenergy[geotop::common::Variables::i_sim] > 0 &&
      fmod(times->time + par->Dt + par->delay_day_recover * 86400.,
           par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.0) <
      1.E-5)
    {
      n_file = (long)((times->time + par->Dt + par->delay_day_recover * 86400.) /
                      (par->output_surfenergy[geotop::common::Variables::i_sim] *
                       3600.0));
      write_suffix(NNNNN, n_file, 1);
      s1 = NNNNN + string("");
      s2 = s1 + "";
      if (geotop::common::Variables::files[fradnet] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fradnet] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->Rn_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->Rn_mean.resize(egy->Rn_mean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fradLWin] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fradLWin] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->LWin_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->LWin_mean.resize(egy->LWin_mean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fradLW] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fradLW] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->LW_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->LW_mean.resize(egy->LW_mean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fradSW] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fradSW] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->SW_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->SW_mean.resize(egy->SW_mean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fradSWin] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fradSWin] + s2;
          write_map_vector(name, 0, par->format_out, egy->Rswdown_mean,
                           geotop::common::Variables::UV,
                           geotop::input::gDoubleNoValue, top->j_cont,
                           geotop::common::Variables::Nr,
                           geotop::common::Variables::Nc);
          egy->Rswdown_mean.resize(egy->Rswdown_mean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fradSWinbeam] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fradSWinbeam] + s2;
          write_map_vector(name, 0, par->format_out, egy->Rswbeam_mean,
                           geotop::common::Variables::UV,
                           geotop::input::gDoubleNoValue, top->j_cont,
                           geotop::common::Variables::Nr,
                           geotop::common::Variables::Nc);
          egy->Rswbeam_mean.resize(egy->Rswbeam_mean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fshadow] !=
          geotop::input::gStringNoValue)
        {
          for (i = 1; i <= par->total_pixel; i++)
            {
              if (egy->nDt_sun[i] > 0)
                {
                  V[i] = egy->nDt_shadow[i] / (double)(egy->nDt_sun[i]);
                }
              else
                {
                  V[i] = -1.;
                }
            }

          name = geotop::common::Variables::files[fshadow] + s2;
          write_map_vector(
            name, 0, par->format_out, V, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->nDt_shadow.resize(egy->nDt_shadow.size(), 0.0);
          egy->nDt_sun.resize(egy->nDt_sun.size(), 0.0);
        }

      // GROUND HEAT FLUX
      if (geotop::common::Variables::files[fG] != geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fG] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->SEB_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->SEB_mean.resize(egy->SEB_mean.size(), 0.0);
        }

      // SENSIBLE HEAT FLUX
      if (geotop::common::Variables::files[fH] != geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fH] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->H_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->H_mean.resize(egy->H_mean.size(), 0.0);
        }

      // LATENT HEAT FLUX
      if (geotop::common::Variables::files[fLE] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fLE] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->ET_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->ET_mean.resize(egy->ET_mean.size(), 0.0);
        }

      // SURFACE TEMPERATURE
      if (geotop::common::Variables::files[fTs] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fTs] + s2;
          write_map_vector(
            name, 0, par->format_out, egy->Ts_mean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          egy->Ts_mean.resize(egy->Ts_mean.size(), 0.0);
        }
    }

  // vegetation variables
  if (par->output_vegetation[geotop::common::Variables::i_sim] > 0 &&
      fmod(times->time + par->Dt + par->delay_day_recover * 86400.,
           par->output_vegetation[geotop::common::Variables::i_sim] * 3600.0) <
      1.E-5)
    {
      n_file = (long)((times->time + par->Dt + par->delay_day_recover * 86400.) /
                      (par->output_vegetation[geotop::common::Variables::i_sim] *
                       3600.0));
      write_suffix(NNNNN, n_file, 1);
      s1 = NNNNN + string("");
      s1 = NNNNN + string(RRRRR);
      s2 = s1 + "";

      // INTERCEPTED PRECIPITATION
      if (geotop::common::Variables::files[fcint] !=
          geotop::input::gStringNoValue)
        {
          temp1 = geotop::common::Variables::files[fcint] + string("water");
          temp2 = temp1 + s2;
          write_map_vector(
            temp2, 0, par->format_out, sl->VS->wrain, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);

          temp1 = geotop::common::Variables::files[fcint] + string("snow");
          temp2 = temp1 + s2;
          write_map_vector(
            temp2, 0, par->format_out, sl->VS->wsnow, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
        }
    }

  // METEO
  if (par->output_meteo[geotop::common::Variables::i_sim] > 0 &&
      fmod(times->time + par->Dt + par->delay_day_recover * 86400.,
           par->output_meteo[geotop::common::Variables::i_sim] * 3600.0) <
      1.E-5)
    {
      n_file =
        (long)((times->time + par->Dt + par->delay_day_recover * 86400.) /
               (par->output_meteo[geotop::common::Variables::i_sim] * 3600.0));

      write_suffix(NNNNN, n_file, 1);
      s1 = NNNNN + string("");
      s2 = s1 + "";
      // s2 = s1 + SSSSS ;

      // AIR TEMPERATURE
      if (geotop::common::Variables::files[fTa] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fTa] + s2;
          write_map_vector(
            name, 0, par->format_out, met->Tamean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          met->Tamean.resize(met->Tamean.size(), 0.0);
        }

      // PRECIPITATION
      if (geotop::common::Variables::files[fprec] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fprec] + string("TOTAL");
          temp1 = name + s2;
          write_map_vector(temp1, 0, par->format_out, wat->PrTOT_mean,
                           geotop::common::Variables::UV,
                           geotop::input::gDoubleNoValue, top->j_cont,
                           geotop::common::Variables::Nr,
                           geotop::common::Variables::Nc);
          wat->PrTOT_mean.resize(wat->PrTOT_mean.size(), 0.0);

          name = geotop::common::Variables::files[fprec] + string("SNOW");
          temp1 = name + s2;
          write_map_vector(temp1, 0, par->format_out, wat->PrSNW_mean,
                           geotop::common::Variables::UV,
                           geotop::input::gDoubleNoValue, top->j_cont,
                           geotop::common::Variables::Nr,
                           geotop::common::Variables::Nc);
          wat->PrSNW_mean.resize(wat->PrSNW_mean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fwspd] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fwspd] + s2;
          write_map_vector(
            name, 0, par->format_out, met->Vspdmean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          met->Vspdmean.resize(met->Vspdmean.size(), 0.0);
        }

      if (geotop::common::Variables::files[fwdir] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[fwdir] + s2;
          write_map_vector(
            name, 0, par->format_out, met->Vdirmean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          met->Vdirmean.resize(met->Vdirmean.size(), 0.0);
        }

      if (geotop::common::Variables::files[frh] !=
          geotop::input::gStringNoValue)
        {
          name = geotop::common::Variables::files[frh] + s2;
          write_map_vector(
            name, 0, par->format_out, met->RHmean, geotop::common::Variables::UV,
            geotop::input::gDoubleNoValue, top->j_cont,
            geotop::common::Variables::Nr, geotop::common::Variables::Nc);
          met->RHmean.resize(met->RHmean.size(), 0.0);
        }
    }

  /**********************************************************************************************************/
  /**********************************************************************************************************/
  // SPECIAL PLOTS AT SOME DAYS
  /**********************************************************************************************************/
  /**********************************************************************************************************/

  if (times->JD_plots.size() > 1 &&
      (size_t)(times->iplot) < times->JD_plots.size())
    {
      i = times->iplot;
      j = 2 * i - 1;
      if (fabs(par->init_date + (times->time + par->Dt) / 86400. -
               times->JD_plots[j + 1]) < 1.E-5)
        {
          lg->logf("Printing plot number %ld \n", i);

          V.resize(par->total_pixel + 1);

          if (geotop::common::Variables::files[pH] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= par->total_pixel; i++)
                {
                  V[i] = egy->Hgplot[i] + egy->Hvplot[i];
                }

              plot(geotop::common::Variables::files[pH], i, V, par->format_out,
                   top->j_cont);
            }

          if (geotop::common::Variables::files[pLE] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= par->total_pixel; i++)
                {
                  V[i] = egy->LEgplot[i] + egy->LEvplot[i];
                }

              plot(geotop::common::Variables::files[pLE], i, V, par->format_out,
                   top->j_cont);
            }

          if (geotop::common::Variables::files[pG] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= par->total_pixel; i++)
                {
                  V[i] = egy->SWgplot[i] + egy->LWgplot[i] - egy->Hgplot[i] -
                         egy->LEgplot[i];
                }
              plot(geotop::common::Variables::files[pG], i, V, par->format_out,
                   top->j_cont);
            }

          if (geotop::common::Variables::files[pth] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= par->total_pixel; i++)
                {
                  V[i] = sl->th[1][i];
                }
              plot(geotop::common::Variables::files[pth], i, V, par->format_out,
                   top->j_cont);
            }

          if (geotop::common::Variables::files[pth] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= par->total_pixel; i++)
                {
                  V[i] = sl->th[1][i];
                }
              plot(geotop::common::Variables::files[pth], i, V, par->format_out,
                   top->j_cont);
            }

          if (geotop::common::Variables::files[pVspd] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= par->total_pixel; i++)
                {
                  V[i] = sqrt(pow(met->Vxplot[i], 2.0) + pow(met->Vyplot[i], 2.0));
                }
              plot(geotop::common::Variables::files[pVspd], i, V, par->format_out,
                   top->j_cont);
            }

          if (geotop::common::Variables::files[pVdir] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= par->total_pixel; i++)
                {
                  V[i] = 270.0 -
                         (180. / GTConst::Pi) * atan2(met->Vyplot[i], met->Vxplot[i]);
                  if (V[i] >= 360.0) V[i] -= 360.0;
                }
              plot(geotop::common::Variables::files[pVdir], i, V, par->format_out,
                   top->j_cont);
            }

          if (geotop::common::Variables::files[pHg] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pHg], i, egy->Hgplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pLEg] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pLEg], i, egy->LEgplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pHv] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pHv], i, egy->Hvplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pLEv] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pLEv], i, egy->LEvplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pSWin] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pSWin], i, egy->SWinplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pSWg] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pSWg], i, egy->SWgplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pSWv] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pSWv], i, egy->SWvplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pLWin] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pLWin], i, egy->LWinplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pLWg] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pLWg], i, egy->LWgplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pLWv] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pLWv], i, egy->LWvplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pTs] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pTs], i, egy->Tsplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pTg] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pTg], i, egy->Tgplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pTv] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pTv], i, egy->Tvplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pTa] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pTa], i, met->Taplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pD] != geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pD], i, snow->Dplot,
                 par->format_out, top->j_cont);
          if (geotop::common::Variables::files[pRH] !=
              geotop::input::gStringNoValue)
            plot(geotop::common::Variables::files[pRH], i, met->RHplot,
                 par->format_out, top->j_cont);

          if (geotop::common::Variables::files[pH] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pHg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->Hgplot.resize(egy->Hgplot.size(), 0.0);
          if (geotop::common::Variables::files[pLE] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pLEg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->LEgplot.resize(egy->LEgplot.size(), 0.0);
          if (geotop::common::Variables::files[pH] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pHv] !=
              geotop::input::gStringNoValue)
            egy->Hvplot.resize(egy->Hvplot.size(), 0.);
          if (geotop::common::Variables::files[pLE] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pLEv] !=
              geotop::input::gStringNoValue)
            egy->LEvplot.resize(egy->LEvplot.size(), 0.);
          if (geotop::common::Variables::files[pSWin] !=
              geotop::input::gStringNoValue)
            egy->SWinplot.resize(egy->SWinplot.size(), 0.);
          if (geotop::common::Variables::files[pSWg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->SWgplot.resize(egy->SWgplot.size(), 0.);
          if (geotop::common::Variables::files[pSWv] !=
              geotop::input::gStringNoValue)
            egy->SWvplot.resize(egy->SWvplot.size(), 0.);
          if (geotop::common::Variables::files[pLWin] !=
              geotop::input::gStringNoValue)
            egy->LWinplot.resize(egy->LWinplot.size(), 0.);
          if (geotop::common::Variables::files[pLWg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->LWgplot.resize(egy->LWgplot.size(), 0.);
          if (geotop::common::Variables::files[pLWv] !=
              geotop::input::gStringNoValue)
            egy->LWvplot.resize(egy->LWvplot.size(), 0.);
          if (geotop::common::Variables::files[pTs] !=
              geotop::input::gStringNoValue)
            egy->Tsplot.resize(egy->Tsplot.size(), 0.);
          if (geotop::common::Variables::files[pTg] !=
              geotop::input::gStringNoValue)
            egy->Tgplot.resize(egy->Tgplot.size(), 0.);
          if (geotop::common::Variables::files[pTv] !=
              geotop::input::gStringNoValue)
            egy->Tvplot.resize(egy->Tvplot.size(), 0.);
          if (geotop::common::Variables::files[pD] != geotop::input::gStringNoValue)
            snow->Dplot.resize(snow->Dplot.size(), 0.);
          if (geotop::common::Variables::files[pTa] !=
              geotop::input::gStringNoValue)
            met->Taplot.resize(met->Taplot.size(), 0.);
          if (geotop::common::Variables::files[pVspd] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pVdir] !=
              geotop::input::gStringNoValue)
            {
              met->Vxplot.resize(met->Vxplot.size(), 0.);
              met->Vyplot.resize(met->Vyplot.size(), 0.);
            }
          if (geotop::common::Variables::files[pRH] !=
              geotop::input::gStringNoValue)
            met->RHplot.resize(met->RHplot.size(), 0.);

          times->iplot++;
        }
    }

  /**********************************************************************************************************/
  /**********************************************************************************************************/
  // SAVING POINTS (No longer supported since June 3rd 2014)
  /**********************************************************************************************************/
  /**********************************************************************************************************/

  if (par->ContRecovery > 0)
    {
      t_rec += par->Dt;
      // used to be ContRecovery*secinday, replaced with timestep
      if (fabs(t_rec - par->ContRecovery * GTConst::secinday) < 1.E-5)
        {
          lg->logsf(geotop::logger::NOTICE,
                    "Writing continuous-recovering files at time (in day): %f "
                    "par->Dt:%f  ContRecovery:%f\n",
                    ((times->time + par->Dt) / GTConst::secinday), par->Dt,
                    par->ContRecovery);
          t_rec = 0.;
          if (geotop::common::Variables::files[rtime] !=
              geotop::input::gStringNoValue)
            {
              name = geotop::common::Variables::files[rtime] + string(textfile);
              f = fopen(name.c_str(), "w");
              fprintf(f,
                      "Time[s],Time[d],n,i_run,i_sim,cum_time[s],elapsed_time[s]\n");
              fprintf(
                f, "%f,%f,%ld,%ld,%ld,%f,%f",
                times->time + par->Dt + par->delay_day_recover * GTConst::secinday,
                (times->time + par->Dt) / GTConst::secinday + par->delay_day_recover,
                (long)(((times->time + par->Dt) / GTConst::secinday +
                        par->delay_day_recover) /
                       par->ContRecovery),
                geotop::common::Variables::i_run, geotop::common::Variables::i_sim,
                geotop::common::Variables::cum_time,
                geotop::common::Variables::elapsed_time);
              fclose(f);
            }
          // TODO: NNNN should be removed in recovery
          write_suffix(NNNN, 0, 0);

          for (l = 0; l <= geotop::common::Variables::Nl; l++)
            {
              if (geotop::common::Variables::files[rpsi] !=
                  geotop::input::gStringNoValue)
                write_tensorseries_vector(
                  1, l, 0, geotop::common::Variables::files[rpsi], -1,
                  par->format_out, sl->SS->P, geotop::common::Variables::UV,
                  geotop::input::gDoubleNoValue, top->j_cont,
                  geotop::common::Variables::Nr, geotop::common::Variables::Nc);
              if (l > 0)
                {
                  if (geotop::common::Variables::files[riceg] !=
                      geotop::input::gStringNoValue)
                    write_tensorseries_vector(
                      1, l, 0, geotop::common::Variables::files[riceg], -1,
                      par->format_out, sl->SS->thi, geotop::common::Variables::UV,
                      geotop::input::gDoubleNoValue, top->j_cont,
                      geotop::common::Variables::Nr, geotop::common::Variables::Nc);
                  if (geotop::common::Variables::files[rTg] !=
                      geotop::input::gStringNoValue)
                    write_tensorseries_vector(
                      1, l, 0, geotop::common::Variables::files[rTg], -1,
                      par->format_out, sl->SS->T, geotop::common::Variables::UV,
                      geotop::input::gDoubleNoValue, top->j_cont,
                      geotop::common::Variables::Nr, geotop::common::Variables::Nc);
                }
            }

          if (geotop::common::Variables::files[rwcrn] !=
              geotop::input::gStringNoValue)
            {
              name = geotop::common::Variables::files[rwcrn] + string(NNNN);
              write_map_vector(name, -1, par->format_out, sl->VS->wrain,
                               geotop::common::Variables::UV,
                               geotop::input::gDoubleNoValue, top->j_cont,
                               geotop::common::Variables::Nr,
                               geotop::common::Variables::Nc);
            }

          if (geotop::common::Variables::files[rwcsn] !=
              geotop::input::gStringNoValue)
            {
              name = geotop::common::Variables::files[rwcsn] + string(NNNN);
              write_map_vector(name, -1, par->format_out, sl->VS->wsnow,
                               geotop::common::Variables::UV,
                               geotop::input::gDoubleNoValue, top->j_cont,
                               geotop::common::Variables::Nr,
                               geotop::common::Variables::Nc);
            }

          for (i = 1; i <= par->total_pixel; i++)
            {
              if ((long)sl->VS->Tv[i] == geotop::input::gDoubleNoValue)
                sl->VS->Tv[i] = 0.;
            }

          if (geotop::common::Variables::files[rTv] !=
              geotop::input::gStringNoValue)
            {
              name = geotop::common::Variables::files[rTv] + string(NNNN);
              write_map_vector(
                name, -1, par->format_out, sl->VS->Tv, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }

          for (l = 1; l <= par->max_snow_layers; l++)
            {
              if (geotop::common::Variables::files[rDzs] !=
                  geotop::input::gStringNoValue)
                write_tensorseries(1, l, 0, geotop::common::Variables::files[rDzs],
                                   -1, par->format_out, snow->S->Dzl,
                                   geotop::common::Variables::UV,
                                   geotop::input::gDoubleNoValue);
              if (geotop::common::Variables::files[rwls] !=
                  geotop::input::gStringNoValue)
                write_tensorseries(1, l, 0, geotop::common::Variables::files[rwls],
                                   -1, par->format_out, snow->S->w_liq,
                                   geotop::common::Variables::UV,
                                   geotop::input::gDoubleNoValue);
              if (geotop::common::Variables::files[rwis] !=
                  geotop::input::gStringNoValue)
                write_tensorseries(1, l, 0, geotop::common::Variables::files[rwis],
                                   -1, par->format_out, snow->S->w_ice,
                                   geotop::common::Variables::UV,
                                   geotop::input::gDoubleNoValue);
              if (geotop::common::Variables::files[rTs] !=
                  geotop::input::gStringNoValue)
                write_tensorseries(1, l, 0, geotop::common::Variables::files[rTs], -1,
                                   par->format_out, snow->S->T,
                                   geotop::common::Variables::UV,
                                   geotop::input::gDoubleNoValue);
            }

          if (geotop::common::Variables::files[rsnag] !=
              geotop::input::gStringNoValue)
            {
              name = geotop::common::Variables::files[rsnag] + string(NNNN);
              write_map_vector(
                name, -1, par->format_out, snow->age, geotop::common::Variables::UV,
                geotop::input::gDoubleNoValue, top->j_cont,
                geotop::common::Variables::Nr, geotop::common::Variables::Nc);
            }

          if (geotop::common::Variables::files[rns] !=
              geotop::input::gStringNoValue)
            {
              name = geotop::common::Variables::files[rns] + string(NNNN);
              write_map(name, 1, par->format_out, snow->S->lnum,
                        geotop::common::Variables::UV, geotop::input::gDoubleNoValue);
            }

          if (par->max_glac_layers > 0)
            {
              for (l = 1; l <= par->max_glac_layers; l++)
                {
                  if (geotop::common::Variables::files[rDzi] !=
                      geotop::input::gStringNoValue)
                    write_tensorseries(1, l, 0, geotop::common::Variables::files[rDzi],
                                       -1, par->format_out, glac->G->Dzl,
                                       geotop::common::Variables::UV,
                                       geotop::input::gDoubleNoValue);
                  if (geotop::common::Variables::files[rwli] !=
                      geotop::input::gStringNoValue)
                    write_tensorseries(1, l, 0, geotop::common::Variables::files[rwli],
                                       -1, par->format_out, glac->G->w_liq,
                                       geotop::common::Variables::UV,
                                       geotop::input::gDoubleNoValue);
                  if (geotop::common::Variables::files[rwii] !=
                      geotop::input::gStringNoValue)
                    write_tensorseries(1, l, 0, geotop::common::Variables::files[rwii],
                                       -1, par->format_out, glac->G->w_ice,
                                       geotop::common::Variables::UV,
                                       geotop::input::gDoubleNoValue);
                  if (geotop::common::Variables::files[rTi] !=
                      geotop::input::gStringNoValue)
                    write_tensorseries(1, l, 0, geotop::common::Variables::files[rTi],
                                       -1, par->format_out, glac->G->T,
                                       geotop::common::Variables::UV,
                                       geotop::input::gDoubleNoValue);
                }

              if (geotop::common::Variables::files[rni] !=
                  geotop::input::gStringNoValue)
                {
                  name = geotop::common::Variables::files[rni] + string(NNNN);
                  write_map(name, 1, par->format_out, glac->G->lnum,
                            geotop::common::Variables::UV,
                            geotop::input::gDoubleNoValue);
                }
            }

          if (geotop::common::Variables::files[rpsich] !=
              geotop::input::gStringNoValue)
            {
              M.resize(geotop::common::Variables::Nl + 1, par->total_pixel + 1, 0.0);
              for (l = 0; l <= geotop::common::Variables::Nl; l++)
                {
                  for (i = 1; i <= par->total_channel; i++)
                    {
                      r = cnet->r[i];
                      c = cnet->c[i];
                      M[l][top->j_cont[r][c]] = cnet->SS->P[l][i];
                    }
                  write_tensorseries_vector(
                    1, l, 0, geotop::common::Variables::files[rpsich], -1,
                    par->format_out, M, geotop::common::Variables::UV,
                    geotop::input::gDoubleNoValue, top->j_cont,
                    geotop::common::Variables::Nr, geotop::common::Variables::Nc);
                }
            }

          if (geotop::common::Variables::files[rTgch] !=
              geotop::input::gStringNoValue)
            {
              M.resize(geotop::common::Variables::Nl + 1, par->total_pixel + 1, 0.0);
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  for (i = 1; i <= par->total_channel; i++)
                    {
                      r = cnet->r[i];
                      c = cnet->c[i];
                      M[l][top->j_cont[r][c]] = cnet->SS->T[l][i];
                    }
                  write_tensorseries_vector(
                    1, l, 0, geotop::common::Variables::files[rTgch], -1,
                    par->format_out, M, geotop::common::Variables::UV,
                    geotop::input::gDoubleNoValue, top->j_cont,
                    geotop::common::Variables::Nr, geotop::common::Variables::Nc);
                }
            }

          if (geotop::common::Variables::files[ricegch] !=
              geotop::input::gStringNoValue)
            {
              M.resize(geotop::common::Variables::Nl + 1, par->total_pixel + 1, 0.0);
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  for (i = 1; i <= par->total_channel; i++)
                    {
                      M[l][top->j_cont[r][c]] = cnet->SS->thi[l][i];
                    }
                  write_tensorseries_vector(
                    1, l, 0, geotop::common::Variables::files[ricegch], -1,
                    par->format_out, M, geotop::common::Variables::UV,
                    geotop::input::gDoubleNoValue, top->j_cont,
                    geotop::common::Variables::Nr, geotop::common::Variables::Nc);
                }
            }
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_output_headers(long n,
                          Times * /*times*/,
                          Water * /*wat*/,
                          Par *par,
                          Topo *top,
                          Land *land,
                          Soil *sl,
                          Energy * /*egy*/,
                          Snow * /*snow*/,
                          Glacier * /*glac*/)
{
  /*internal auxiliary variables:*/
  size_t i, m, r, c;
  long l, j;

  std::string NNNN = "NNNN";
  std::string rec = "_recNNNN";
  std::string crec = "_crecNNNN";
  string name, temp, temp2;
  long sy;
  short lu, first_column;
  GeoVector<double> root_fraction;
  FILE *f;

  geotop::logger::GlobalLogger *lg =
    geotop::logger::GlobalLogger::getInstance();

  if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

  // DISCHARGE
  if (par->state_discharge == 1 &&
      geotop::common::Variables::files[fQ] != geotop::input::gStringNoValue)
    {
      if (par->n_ContRecovery > 0)
        {
          temp = geotop::common::Variables::files[fQ] + string(crec);
          name = temp + textfile;
        }
      else
        {
          name = geotop::common::Variables::files[fQ] + string(textfile);
        }

      f = fopen(name.c_str(), "w");
      if (f == NULL)
        {
          lg->logsf(geotop::logger::CRITICAL, "Error opening file: %s\n",
                    name.c_str());
          t_error("Error opening file...");
        }
      fprintf(f,
              "DATE[day/month/year "
              "hour:min],t[days],JDfrom0,JD,Qtot[m3/s],Vsup/Dt[m3/s],Vsub/Dt[m3/"
              "s],Vchannel[m3],Qoutlandsup[m3/s],Qoutlandsub[m3/s],Qoutbottom[m3/"
              "s]\n");
      fclose(f);
    }

  if (par->state_pixel == 1)
    {
      // output matrix and vectors
      m = (long)otot;
      geotop::common::Variables::odpnt = (double **)malloc(m * sizeof(double *));
      geotop::common::Variables::odp = (double **)malloc(m * sizeof(double *));
      for (i = 0; i < otot; i++)
        {
          geotop::common::Variables::odpnt[i] =
            (double *)malloc(par->rc.getRows() * sizeof(double));
          geotop::common::Variables::odp[i] =
            (double *)malloc(par->rc.getRows() * sizeof(double));
          // to check: is this -1 below needed ? SC26.12.2013
          for (j = 0; j < long(par->rc.getRows() - 1); j++)
            {
              geotop::common::Variables::odpnt[i][j] = 0.;
              geotop::common::Variables::odp[i][j] = 0.;
            }
        }

      if (geotop::common::Variables::files[fpointwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fpointwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              ;
              name =
                geotop::common::Variables::files[fpointwriteend] + string(textfile);
            }

          geotop::common::Variables::ffpoint = fopen(name.c_str(), "w");
          first_column = 1;
          for (j = 0; j < geotop::common::Variables::nopnt; j++)
            {
              if (first_column == 0)
                {
                  fprintf(geotop::common::Variables::ffpoint, ",");
                }
              else
                {
                  first_column = 0;
                }
              if (geotop::common::Variables::opnt[j] >= 0)
                {
                  fprintf(
                    geotop::common::Variables::ffpoint, "%s",
                    geotop::common::Variables::hpnt[geotop::common::Variables::opnt[j]]
                    .c_str());
                }
              else
                {
                  fprintf(geotop::common::Variables::ffpoint, "None");
                }
            }
          fprintf(geotop::common::Variables::ffpoint, "\n");
        }

      if (geotop::common::Variables::files[fsnTzwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fsnTzwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[fsnTzwriteend] + string(textfile);
            }

          geotop::common::Variables::ffsnow = fopen(name.c_str(), "w");

          if ((long)par->snow_plot_depths[1] != geotop::input::gDoubleNoValue)
            {
              m = par->snow_plot_depths.size();
            }
          else
            {
              m = par->max_snow_layers;
            }

          first_column = 1;
          for (j = 0; j < geotop::common::Variables::nosnw; j++)
            {
              if (first_column == 0)
                {
                  fprintf(geotop::common::Variables::ffsnow, ",");
                }
              else
                {
                  first_column = 0;
                }
              if (geotop::common::Variables::osnw[j] >= 0 &&
                  geotop::common::Variables::osnw[j] <= 5)
                {
                  fprintf(
                    geotop::common::Variables::ffsnow, "%s",
                    geotop::common::Variables::hsnw[geotop::common::Variables::osnw[j]]
                    .c_str());
                }
              else if (geotop::common::Variables::osnw[j] >= 6 &&
                       geotop::common::Variables::osnw[j] < long(6 + 3 * m))
                {
                  l = (long)fmod((double)geotop::common::Variables::osnw[j] - 6.,
                                 (double)m) +
                      1;
                  n = floor(((double)geotop::common::Variables::osnw[j] - 6.) /
                            (double)m) +
                      6;
                  if ((long)par->snow_plot_depths[1] != geotop::input::gDoubleNoValue)
                    {
                      fprintf(geotop::common::Variables::ffsnow, "%s(%f)",
                              geotop::common::Variables::hsnw[n].c_str(),
                              par->snow_plot_depths[l]);
                    }
                  else
                    {
                      fprintf(geotop::common::Variables::ffsnow, "%s(%ld)",
                              geotop::common::Variables::hsnw[n].c_str(), l);
                    }
                }
              else if (geotop::common::Variables::osnw[j] >= long(6 + 3 * m))
                {
                  l = (long)fmod(
                        (double)geotop::common::Variables::osnw[j] - 6. - 3 * (double)m,
                        (double)par->max_snow_layers) +
                      1;
                  n = floor(((double)geotop::common::Variables::osnw[j] - 6. -
                             3. * (double)m) /
                            (double)par->max_snow_layers) +
                      6 + 3;
                  fprintf(geotop::common::Variables::ffsnow, "%s(%ld)",
                          geotop::common::Variables::hsnw[n].c_str(), l);
                }
              else
                {
                  fprintf(geotop::common::Variables::ffsnow, "None");
                }
            }
          fprintf(geotop::common::Variables::ffsnow, "\n");
        }

      if (par->max_glac_layers > 0)
        {
          if (geotop::common::Variables::files[fglzwriteend] !=
              geotop::input::gStringNoValue)
            {
              if (par->n_ContRecovery > 0)
                {
                  temp =
                    geotop::common::Variables::files[fpointwriteend] + string(crec);
                  name = temp + textfile;
                }
              else
                {
                  name =
                    geotop::common::Variables::files[fpointwriteend] + string(textfile);
                }

              geotop::common::Variables::ffglac = fopen(name.c_str(), "w");

              if ((long)par->glac_plot_depths[1] != geotop::input::gDoubleNoValue)
                {
                  m = par->glac_plot_depths.size();
                }
              else
                {
                  m = par->max_glac_layers;
                }
              first_column = 1;
              for (j = 0; j < geotop::common::Variables::noglc; j++)
                {
                  if (first_column == 0)
                    {
                      fprintf(geotop::common::Variables::ffglac, ",");
                    }
                  else
                    {
                      first_column = 0;
                    }
                  if (geotop::common::Variables::oglc[j] >= 0 &&
                      geotop::common::Variables::oglc[j] <= 5)
                    {
                      fprintf(geotop::common::Variables::ffglac, "%s",
                              geotop::common::Variables::hglc
                              [geotop::common::Variables::oglc[j]]
                              .c_str());
                    }
                  else if (geotop::common::Variables::oglc[j] >= 6 &&
                           geotop::common::Variables::oglc[j] < long(6 + 3 * m))
                    {
                      l = (long)fmod((double)geotop::common::Variables::oglc[j] - 6.,
                                     (double)m) +
                          1;
                      n = floor(((double)geotop::common::Variables::oglc[j] - 6.) /
                                (double)m) +
                          6;
                      if ((long)par->glac_plot_depths[1] !=
                          geotop::input::gDoubleNoValue)
                        {
                          fprintf(geotop::common::Variables::ffglac, "%s(%f)",
                                  geotop::common::Variables::hglc[n].c_str(),
                                  par->glac_plot_depths[l]);
                        }
                      else
                        {
                          fprintf(geotop::common::Variables::ffglac, "%s(%ld)",
                                  geotop::common::Variables::hglc[n].c_str(), l);
                        }
                    }
                  else if (geotop::common::Variables::oglc[j] >= long(6 + 3 * m))
                    {
                      l = (long)fmod((double)geotop::common::Variables::oglc[j] - 6. -
                                     3 * (double)m,
                                     (double)par->max_glac_layers) +
                          1;
                      n = floor(((double)geotop::common::Variables::oglc[j] - 6. -
                                 3. * (double)m) /
                                (double)par->max_glac_layers) +
                          6 + 3;
                      fprintf(geotop::common::Variables::ffglac, "%s(%ld)",
                              geotop::common::Variables::hglc[n].c_str(), l);
                    }
                  else
                    {
                      fprintf(geotop::common::Variables::ffglac, "None");
                    }
                }
              fprintf(geotop::common::Variables::ffglac, "\n");
            }
        }

      if (geotop::common::Variables::files[fTzwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fTzwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name = geotop::common::Variables::files[fTzwriteend] + string(textfile);
            }

          geotop::common::Variables::ffT = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::ffT, par->soil_plot_depths,
                            sl->pa);
        }

      if (geotop::common::Variables::files[fTzavwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fTzavwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[fTzavwriteend] + string(textfile);
            }

          geotop::common::Variables::ffTav = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::ffTav, par->soil_plot_depths,
                            sl->pa);
        }

      if (geotop::common::Variables::files[fpsiztotwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp =
                geotop::common::Variables::files[fpsiztotwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[fpsiztotwriteend] + string(textfile);
            }

          geotop::common::Variables::ffpsitot = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::ffpsitot,
                            par->soil_plot_depths, sl->pa);
        }

      if (geotop::common::Variables::files[fpsizwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fpsizwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[fpsizwriteend] + string(textfile);
            }

          geotop::common::Variables::ffpsi = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::ffpsi, par->soil_plot_depths,
                            sl->pa);
        }

      if (geotop::common::Variables::files[fliqzwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fliqzwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[fliqzwriteend] + string(textfile);
            }

          geotop::common::Variables::ffliq = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::ffliq, par->soil_plot_depths,
                            sl->pa);
        }

      if (geotop::common::Variables::files[fliqzavwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fliqzavwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[fliqzavwriteend] + string(textfile);
            }

          geotop::common::Variables::ffliqav = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::ffliqav,
                            par->soil_plot_depths, sl->pa);
        }

      if (geotop::common::Variables::files[ficezwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[ficezwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[ficezwriteend] + string(textfile);
            }

          geotop::common::Variables::ffice = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::ffice, par->soil_plot_depths,
                            sl->pa);
        }

      if (geotop::common::Variables::files[ficezavwriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[ficezavwriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[ficezavwriteend] + string(textfile);
            }

          geotop::common::Variables::fficeav = fopen(name.c_str(), "w");
          write_soil_header(geotop::common::Variables::fficeav,
                            par->soil_plot_depths, sl->pa);
        }

      // is this ok ? Nl +1 ??

      root_fraction.resize(geotop::common::Variables::Nl + 1);

      // DATA POINTS

      for (i = 1; i < par->rc.getRows(); i++)
        {
          write_suffix(NNNN, par->IDpoint[i], 0);
          r = par->rc[i][1];
          c = par->rc[i][2];
          sy = sl->type[r][c];
          lu = (short)land->LC[r][c];

          if (geotop::common::Variables::files[fpoint] !=
              geotop::input::gStringNoValue &&
              par->point_sim != 1)
            {
              name = geotop::common::Variables::files[fpoint] + string("_info_");
              temp = name + NNNN;
              temp2 = temp + textfile;
              f = fopen(temp2.c_str(), "w");

              fprintf(f,
                      " The main properties of the pixel E=%15.3f N=%15.3f, row=%4ld "
                      "col=%4ld are:\n",
                      par->chkpt[i][ptX], par->chkpt[i][ptY], r, c);
              fprintf(f, " Elevation above sea level: %10.3f m\n", top->Z0[r][c]);
              fprintf(f, " Gauckler-Strickler [m^1/3/s]: %f\n", land->ty[lu][jcm]);
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " Residual water content[-] of the layer %ld: %f\n", l,
                          sl->pa[sy][jres][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " Saturated water content[-] of the layer %ld: %f\n", l,
                          sl->pa[sy][jsat][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " Alpha of van Genuchten[mm^-1] of the layer %ld: %f\n", l,
                          sl->pa[sy][ja][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " n of van Genuchten[-] of the layer %ld: %f\n", l,
                          sl->pa[sy][jns][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " m of van Genuchten[-] of the layer %ld: %f\n", l,
                          1 - 1 / sl->pa[sy][jns][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " v of van Genuchten[-] of the layer %ld: %f\n", l,
                          sl->pa[sy][jv][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f,
                          " Water content of wilting point [-] of the layer %ld: %f\n",
                          l, sl->pa[sy][jwp][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f,
                          " Water content of field capacity [-] of the layer %ld: %f\n",
                          l, sl->pa[sy][jfc][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " Kv_sat of layer %ld [mm/s]: %f\n", l,
                          sl->pa[sy][jKn][l]);
                }
              for (l = 1; l <= geotop::common::Variables::Nl; l++)
                {
                  fprintf(f, " Kh_sat of layer %ld [mm/s]: %f\n", l,
                          sl->pa[sy][jKl][l]);
                }

              fprintf(f, " Terrain elevation [m]: %f\n", top->Z0[r][c]);
              fprintf(f, " Sky view factor [-]: %f\n", top->sky[r][c]);
              fprintf(f, " The pixel-type is %d \n", top->pixel_type[r][c]);
              fprintf(f, " Aspect [deg] [0=Nord, clockwise]: %f \n",
                      top->aspect[r][c]);
              fprintf(f, " Mean slope of the pixel [deg]: %f \n", top->slope[r][c]);
              fprintf(f, " Land use number is %d \n", (short)land->LC[r][c]);

              for (l = 1; l < long(land->root_fraction.getCols()); l++)
                {
                  fprintf(f, " The root fraction [-] of layer %ld: %f\n", l,
                          land->root_fraction[lu][l]);
                }

              fprintf(f, " Surface fraction of land covered by vegetation [-]: %f \n",
                      land->ty[lu][jcf]);
              fprintf(f, " Leaf and Stem Area Index [-]: %f \n", land->ty[lu][jLSAI]);
              fprintf(f, " Momentum roughness length z0soil [m]: %f \n",
                      land->ty[lu][jz0]);
              fprintf(f, " Vegetation height [m]: %f \n", land->ty[lu][jHveg]);

              fprintf(f, " \n");
              fclose(f);
            }

          if (geotop::common::Variables::files[fpoint] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fpoint] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              if (f == NULL)
                {
                  lg->logsf(geotop::logger::CRITICAL, "Error opening file: %s\n",
                            name.c_str());
                  t_error("Error opening file...");
                }

              first_column = 1;
              for (j = 0; j < geotop::common::Variables::nopnt; j++)
                {
                  if (first_column == 0)
                    {
                      fprintf(f, ",");
                    }
                  else
                    {
                      first_column = 0;
                    }
                  if (geotop::common::Variables::opnt[j] >= 0)
                    {
                      fprintf(f, "%s",
                              geotop::common::Variables::hpnt
                              [geotop::common::Variables::opnt[j]]
                              .c_str());
                    }
                  else
                    {
                      fprintf(f, "None");
                    }
                }
              fprintf(f, "\n");
              fclose(f);
            }

          if (par->max_glac_layers > 0)
            {
              if (geotop::common::Variables::files[fglz] !=
                  geotop::input::gStringNoValue)
                {
                  temp = geotop::common::Variables::files[fglz] + string(NNNN);

                  if (par->n_ContRecovery > 0)
                    {
                      temp2 = temp + crec;
                      name = temp2 + textfile;
                    }
                  else
                    {
                      name = temp + textfile;
                    }

                  if ((long)par->glac_plot_depths[1] != geotop::input::gDoubleNoValue)
                    {
                      m = par->glac_plot_depths.size();
                    }
                  else
                    {
                      m = par->max_glac_layers;
                    }

                  f = fopen(name.c_str(), "w");
                  first_column = 1;
                  for (j = 0; j < geotop::common::Variables::noglc; j++)
                    {
                      if (first_column == 0)
                        {
                          fprintf(f, ",");
                        }
                      else
                        {
                          first_column = 0;
                        }
                      if (geotop::common::Variables::oglc[j] >= 0 &&
                          geotop::common::Variables::oglc[j] <= 5)
                        {
                          fprintf(f, "%s",
                                  geotop::common::Variables::hglc
                                  [geotop::common::Variables::oglc[j]]
                                  .c_str());
                        }
                      else if (geotop::common::Variables::oglc[j] >= 6 &&
                               geotop::common::Variables::oglc[j] < long(6 + 3 * m))
                        {
                          l = (long)fmod((double)geotop::common::Variables::oglc[j] - 6.,
                                         (double)m) +
                              1;
                          n = floor(((double)geotop::common::Variables::oglc[j] - 6.) /
                                    (double)m) +
                              6;
                          if ((long)par->glac_plot_depths[1] !=
                              geotop::input::gDoubleNoValue)
                            {
                              fprintf(f, "%s(%f)", geotop::common::Variables::hglc[n].c_str(),
                                      par->glac_plot_depths[l]);
                            }
                          else
                            {
                              fprintf(f, "%s(%ld)",
                                      geotop::common::Variables::hglc[n].c_str(), l);
                            }
                        }
                      else if (geotop::common::Variables::oglc[j] >= long(6 + 3 * m))
                        {
                          l = (long)fmod((double)geotop::common::Variables::oglc[j] - 6. -
                                         3 * (double)m,
                                         (double)par->max_glac_layers) +
                              1;
                          n = floor(((double)geotop::common::Variables::oglc[j] - 6. -
                                     3. * (double)m) /
                                    (double)par->max_glac_layers) +
                              6 + 3;
                          fprintf(f, "%s(%ld)", geotop::common::Variables::hglc[n].c_str(),
                                  l);
                        }
                      else
                        {
                          fprintf(f, "None");
                        }
                    }
                  fprintf(f, "\n");
                  fclose(f);
                }
            }

          if (geotop::common::Variables::files[fTz] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fTz] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          if (geotop::common::Variables::files[fTzav] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fTzav] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          if (geotop::common::Variables::files[fpsiztot] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fpsiztot] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          if (geotop::common::Variables::files[fpsiz] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fpsiz] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          if (geotop::common::Variables::files[fliqz] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fliqz] + string(NNNN);
              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          if (geotop::common::Variables::files[fliqzav] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fliqzav] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          if (geotop::common::Variables::files[ficez] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[ficez] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          if (geotop::common::Variables::files[ficezav] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[ficezav] + string(NNNN);

              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }

              f = fopen(name.c_str(), "w");
              write_soil_header(f, par->soil_plot_depths, sl->pa);
              fclose(f);
            }

          // snow headers
          //          "SnowIceContentProfileFile",//fsniz
          if (geotop::common::Variables::files[fsniz] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fsniz] + string(NNNN);
              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }
              if ((long)par->snow_plot_depths[1] != geotop::input::gDoubleNoValue)
                {
                  m = par->snow_plot_depths.size();
                }
              else
                {
                  m = par->max_snow_layers;
                }
              f = fopen(name.c_str(), "w");
              write_snow_header(f, par->snow_plot_depths, par->max_snow_layers);
              fclose(f);
            }

          //"SnowTempProfileFile",//fsnTz
          if (geotop::common::Variables::files[fsnTz] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fsnTz] + string(NNNN);
              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }
              if ((long)par->snow_plot_depths[1] != geotop::input::gDoubleNoValue)
                {
                  m = par->snow_plot_depths.size();
                }
              else
                {
                  m = par->max_snow_layers;
                }
              f = fopen(name.c_str(), "w");
              write_snow_header(f, par->snow_plot_depths, par->max_snow_layers);
              fclose(f);
            }
          //"SnowLiqContentProfileFile",//fsnlz
          if (geotop::common::Variables::files[fsnlz] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fsnlz] + string(NNNN);
              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }
              if ((long)par->snow_plot_depths[1] != geotop::input::gDoubleNoValue)
                {
                  m = par->snow_plot_depths.size();
                }
              else
                {
                  m = par->max_snow_layers;
                }
              f = fopen(name.c_str(), "w");
              write_snow_header(f, par->snow_plot_depths, par->max_snow_layers);
              fclose(f);
            }

          //"SnowDepthLayersFile",//fsndz
          if (geotop::common::Variables::files[fsndz] !=
              geotop::input::gStringNoValue)
            {
              temp = geotop::common::Variables::files[fsndz] + string(NNNN);
              if (par->n_ContRecovery > 0)
                {
                  temp2 = temp + crec;
                  name = temp2 + textfile;
                }
              else
                {
                  name = temp + textfile;
                }
              if ((long)par->snow_plot_depths[1] != geotop::input::gDoubleNoValue)
                {
                  m = par->snow_plot_depths.size();
                }
              else
                {
                  m = par->max_snow_layers;
                }
              f = fopen(name.c_str(), "w");
              write_snow_header(f, par->snow_plot_depths, par->max_snow_layers);
              fclose(f);
            }
        }
    }

  m = (long)ootot;
  geotop::common::Variables::odbsn = (double *)malloc(m * sizeof(double));
  geotop::common::Variables::odb = (double *)malloc(m * sizeof(double));
  for (i = 0; i < ootot; i++)
    {
      geotop::common::Variables::odbsn[i] = 0.;
      geotop::common::Variables::odb[i] = 0.;
    }

  if (par->state_basin == 1)
    {
      // DATA BASIN
      if (geotop::common::Variables::files[fbaswriteend] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fbaswriteend] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name =
                geotop::common::Variables::files[fbaswriteend] + string(textfile);
            }

          geotop::common::Variables::ffbas = fopen(name.c_str(), "w");

          first_column = 1;
          for (j = 0; j < geotop::common::Variables::nobsn; j++)
            {
              if (first_column == 0)
                {
                  fprintf(geotop::common::Variables::ffbas, ",");
                }
              else
                {
                  first_column = 0;
                }
              if (geotop::common::Variables::obsn[j] >= 0)
                {
                  fprintf(
                    geotop::common::Variables::ffbas, "%s",
                    geotop::common::Variables::hbsn[geotop::common::Variables::obsn[j]]
                    .c_str());
                }
              else
                {
                  fprintf(geotop::common::Variables::ffbas, "None");
                }
            }
          fprintf(geotop::common::Variables::ffbas, "\n");
        }

      if (geotop::common::Variables::files[fbas] !=
          geotop::input::gStringNoValue)
        {
          if (par->n_ContRecovery > 0)
            {
              temp = geotop::common::Variables::files[fbas] + string(crec);
              name = temp + textfile;
            }
          else
            {
              name = geotop::common::Variables::files[fbas] + string(textfile);
            }

          f = fopen(name.c_str(), "w");

          first_column = 1;
          for (j = 0; j < geotop::common::Variables::nobsn; j++)
            {
              if (first_column == 0)
                {
                  fprintf(f, ",");
                }
              else
                {
                  first_column = 0;
                }
              if (geotop::common::Variables::obsn[j] >= 0)
                {
                  fprintf(
                    f, "%s",
                    geotop::common::Variables::hbsn[geotop::common::Variables::obsn[j]]
                    .c_str());
                }
              else
                {
                  fprintf(f, "None");
                }
            }
          fprintf(f, "\n");

          fclose(f);
        }
    }

  // SNOW COVERED AREA STATISTICS
  if (par->point_sim != 1 &&
      geotop::common::Variables::files[fSCA] != geotop::input::gStringNoValue)
    {
      if (par->n_ContRecovery > 0)
        {
          temp = geotop::common::Variables::files[fSCA] + string(crec);
          name = temp + textfile;
        }
      else
        {
          name = geotop::common::Variables::files[fSCA] + string(textfile);
        }

      f = fopen(name.c_str(), "w");
      fprintf(f,
              "DATE[day/month/year "
              "hour:min],t[days],JDfrom0,JD,snowDav,SWEav,Tav,Tsav,perc.SFA,perc."
              "SCA\n");
      fclose(f);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_output(long i,
                       long iname,
                       double init_date,
                       double JDfrom0,
                       double /*JD*/,
                       long day,
                       long month,
                       long year,
                       long hour,
                       long minute,
                       const GeoVector<double> &n,
                       Soil *sl,
                       Par *par,
                       double /*psimin*/,
                       double cosslope)
{
  std::string NNNN = "NNNN";
  string name, temp, temp2;

  std::string rec = "_recNNNN";
  std::string crec = "_crecNNNN";
  long l;
  FILE *f;

  write_suffix(NNNN, iname, 0);

  if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

  if (geotop::common::Variables::files[fTz] != geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fTz] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->Tzplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[fTzwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(1, iname, geotop::common::Variables::ffT, day, month, year,
                      hour, minute, JDfrom0, init_date, sl->Tzplot, i, n, sl->pa,
                      cosslope);
    }

  if (geotop::common::Variables::files[fTzav] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fTzav] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->Tzavplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[fTzavwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(1, iname, geotop::common::Variables::ffTav, day, month,
                      year, hour, minute, JDfrom0, init_date, sl->Tzavplot, i, n,
                      sl->pa, cosslope);
    }

  if (geotop::common::Variables::files[fpsiztot] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fpsiztot] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->Ptotzplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[fpsiztotwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(1, iname, geotop::common::Variables::ffpsitot, day, month,
                      year, hour, minute, JDfrom0, init_date, sl->Ptotzplot, i, n,
                      sl->pa, cosslope);
    }

  if (geotop::common::Variables::files[fpsiz] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fpsiz] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(0, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->Pzplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[fpsizwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(0, iname, geotop::common::Variables::ffpsi, day, month,
                      year, hour, minute, JDfrom0, init_date, sl->Pzplot, i, n,
                      sl->pa, cosslope);
    }

  if (geotop::common::Variables::files[fliqz] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fliqz] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->thzplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[fliqzwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(1, iname, geotop::common::Variables::ffliq, day, month,
                      year, hour, minute, JDfrom0, init_date, sl->thzplot, i, n,
                      sl->pa, cosslope);
    }

  if (geotop::common::Variables::files[fliqzav] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fliqzav] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->thzavplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[fliqzavwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(1, iname, geotop::common::Variables::ffliqav, day, month,
                      year, hour, minute, JDfrom0, init_date, sl->thzavplot, i, n,
                      sl->pa, cosslope);
    }

  if (geotop::common::Variables::files[ficez] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[ficez] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->thizplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[ficezwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(1, iname, geotop::common::Variables::ffice, day, month,
                      year, hour, minute, JDfrom0, init_date, sl->thizplot, i, n,
                      sl->pa, cosslope);
    }

  if (geotop::common::Variables::files[ficezav] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[ficezav] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0,
                      init_date, sl->thizavplot, i, n, sl->pa, cosslope);
      fclose(f);
    }

  if (geotop::common::Variables::files[ficezavwriteend] !=
      geotop::input::gStringNoValue)
    {
      write_soil_file(1, iname, geotop::common::Variables::fficeav, day, month,
                      year, hour, minute, JDfrom0, init_date, sl->thizavplot, i,
                      n, sl->pa, cosslope);
    }

  for (l = 1; l <= geotop::common::Variables::Nl; l++)
    {
      if (geotop::common::Variables::files[fTzav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fTzavwriteend] !=
          geotop::input::gStringNoValue)
        sl->Tzavplot[i][l] = 0.0;
      if (geotop::common::Variables::files[fliqzav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[fliqzavwriteend] !=
          geotop::input::gStringNoValue)
        sl->thzavplot[i][l] = 0.0;
      if (geotop::common::Variables::files[ficezav] !=
          geotop::input::gStringNoValue ||
          geotop::common::Variables::files[ficezavwriteend] !=
          geotop::input::gStringNoValue)
        sl->thizavplot[i][l] = 0.0;
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_file(long lmin,
                     long i,
                     FILE *f,
                     long d,
                     long m,
                     long y,
                     long h,
                     long mi,
                     double JDfrom0,
                     double JDfrom0init,
                     const GeoMatrix<double> &var,
                     long row,
                     const GeoVector<double> &n,
                     const GeoTensor<double> &dz,
                     double cosslope)
{
  short first_column = 1;
  long j;
  size_t l;

  for (j = 0; j < geotop::common::Variables::nosl; j++)
    {
      if (first_column == 0)
        {
          fprintf(f, ",");
        }
      else
        {
          first_column = 0;
        }
      if (geotop::common::Variables::osl[j] >= 0)
        {
          if (geotop::common::Variables::osl[j] == 0)
            {
              fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)d, (float)m,
                      (float)y, (float)h, (float)mi);
            }
          else if (geotop::common::Variables::osl[j] == 1)
            {
              fprintf(f, "%f", JDfrom0);
            }
          else if (geotop::common::Variables::osl[j] == 2)
            {
              fprintf(f, "%f", JDfrom0 - JDfrom0init);
            }
          else if (geotop::common::Variables::osl[j] == 3)
            {
              fprintf(f, "%ld", geotop::common::Variables::i_sim);
            }
          else if (geotop::common::Variables::osl[j] == 4)
            {
              fprintf(f, "%ld", geotop::common::Variables::i_run);
            }
          else if (geotop::common::Variables::osl[j] == 5)
            {
              fprintf(f, "%ld", i);
            }
        }
      else
        {
          fprintf(f, "%f", geotop::input::gDoubleNoValue);
        }
    }

  if ((long)n[1] != geotop::input::gDoubleNoValue)
    {
      for (l = 1; l < n.size(); l++)
        {
          fprintf(f, ",%f",
                  interpolate_soil(lmin, n[l] * cosslope,
                                   geotop::common::Variables::Nl, dz, var, row));
        }
    }
  else
    {
      for (l = 1; l <= (size_t)(geotop::common::Variables::Nl); l++)
        {
          fprintf(f, ",%f", var(row, l));
        }
    }

  fprintf(f, " \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_header(FILE *f,
                       const GeoVector<double> &n,
                       const GeoTensor<double> &dz)
{
  short first_column = 1;
  long j;
  size_t l;
  double z = 0.0;

  for (j = 0; j < geotop::common::Variables::nosl; j++)
    {
      if (first_column == 0)
        {
          fprintf(f, ",");
        }
      else
        {
          first_column = 0;
        }
      if (geotop::common::Variables::osl[j] >= 0 &&
          geotop::common::Variables::osl[j] <= 5)
        {
          fprintf(f, "%s",
                  geotop::common::Variables::hsl[geotop::common::Variables::osl[j]]
                  .c_str());
        }
      else
        {
          fprintf(f, "none");
        }
    }

  if ((long)n[1] != geotop::input::gDoubleNoValue)
    {
      for (l = 1; l < n.size(); l++)
        {
          fprintf(f, ",%f", n[l]);
        }
    }
  else
    {
      for (l = 1; l <= (size_t)(geotop::common::Variables::Nl); l++)
        {
          z += dz(1, jdz, l);
          fprintf(f, ",%f ", z - 0.5 * dz(1, jdz, l));
        }
    }

  fprintf(f, " \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void plot(std::string name,
          long i_plot,
          const GeoVector<double> &V,
          short format,
          long **J)
{
  std::string ADS = "iiii";
  string temp;

  write_suffix(ADS, i_plot, 0);
  temp = name + std::string(ADS);

  write_map_vector(temp, 0, format, V, geotop::common::Variables::UV,
                   geotop::input::gDoubleNoValue, J,
                   geotop::common::Variables::Nr,
                   geotop::common::Variables::Nc);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

// to check this function...

double interpolate_soil3(long lmin, double h, long max, double *Dz,
                         double *Q)
{
  double q, z, z0 = 0.;
  long l;

  l = lmin;
  q = geotop::input::gDoubleNoValue;

  do
    {
      if (l == lmin)
        {
          z = z0;
          if (l > 0) z += Dz[l] / 2.;
        }
      else if (l <= max)
        {
          z = z0 + Dz[l] / 2;
          if (l > 1) z += Dz[l - 1] / 2.;
        }
      else
        {
          z = z0 + Dz[max] / 2.;
        }

      if (h < z && h >= z0)
        {
          if (l == lmin)
            {
              q = Q[lmin];
            }
          else if (l <= max)
            {
              q = (Q[l - 1] * (z - h) + Q[l] * (h - z0)) / (z - z0);
            }
          else
            {
              q = Q[max];
            }
        }

      z0 = z;

      l++;

    }
  while ((long)q == geotop::input::gDoubleNoValue && l <= max + 1);

  return q;
}

double interpolate_soil(long lmin,
                        double h,
                        long max,
                        const GeoTensor<double> &Dz,
                        const GeoMatrix<double> &Q,
                        const long &row)
{
  double q, z, z0 = 0.;
  long l;

  l = lmin;
  q = geotop::input::gDoubleNoValue;

  do
    {
      if (l == lmin)
        {
          z = z0;
          if (l > 0) z += Dz(1, jdz, l) / 2.;
        }
      else if (l <= max)
        {
          z = z0 + Dz(1, jdz, l) / 2;
          if (l > 1) z += Dz(1, jdz, l - 1) / 2.;
        }
      else
        {
          z = z0 + Dz(1, jdz, max) / 2.;
        }

      if (h < z && h >= z0)
        {
          if (l == lmin)
            {
              q = Q(row, lmin);
            }
          else if (l <= max)
            {
              q = (Q(row, l - 1) * (z - h) + Q(row, l) * (h - z0)) / (z - z0);
            }
          else
            {
              q = Q(row, max);
            }
        }

      z0 = z;

      l++;

    }
  while ((long)q == geotop::input::gDoubleNoValue && l <= max + 1);

  return q;
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

double interpolate_soil2(long lmin,
                         double h,
                         long max,
                         const GeoTensor<double> &Dz,
                         GeoMatrix<double> &Q,
                         long i)
{
  double q, z, z0 = 0.;
  long l;

  l = lmin;
  q = geotop::input::gDoubleNoValue;

  do
    {
      if (l == lmin)
        {
          z = z0;
          if (l > 0) z += Dz(1, jdz, l) / 2.;
        }
      else if (l <= max)
        {
          z = z0 + Dz(1, jdz, l) / 2.;
          if (l > 1) z += Dz(1, jdz, l - 1) / 2.;
        }
      else
        {
          z = z0 + Dz(1, jdz, max) / 2.;
        }

      if (h < z && h >= z0)
        {
          if (l == lmin)
            {
              q = Q[l][i];
            }
          else if (l <= max)
            {
              q = (Q[l - 1][i] * (z - h) + Q[l][i] * (h - z0)) / (z - z0);
            }
          else
            {
              q = Q[max][i];
            }
        }

      z0 = z;

      l++;

    }
  while ((long)q == geotop::input::gDoubleNoValue && l <= max + 1);

  return q;
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_tensorseries_soil(long lmin,
                             std::string suf,
                             std::string filename,
                             short type,
                             short format,
                             GeoMatrix<double> &T,
                             const GeoVector<double> &n,
                             long **J,
                             GeoMatrix<long> &RC,
                             GeoTensor<double> &dz,
                             GeoMatrix<double> &slope,
                             short vertical)
{
  char LLLLL[] = {"LLLLL"};
  string temp1, temp2;
  long npoints = T.getCols();
  double cosslope = 1.;
  GeoVector<double> V(npoints);

  for (size_t l = 1; l < n.size(); l++)
    {
      temp1 = LLLLL + string(suf);
      write_suffix(temp1, l, 1);

      for (long i = 1; i < npoints; i++)
        {
          if (vertical == 1)
            cosslope = cos(Fmin(GTConst::max_slope, slope[RC[i][1]][RC[i][2]]) *
                           GTConst::Pi / 180.);
          V[i] = interpolate_soil2(lmin, n[l] * cosslope,
                                   geotop::common::Variables::Nl, dz, T, i);
        }

      temp2 = filename + temp1;
      write_map_vector(temp2, type, format, V, geotop::common::Variables::UV,
                       geotop::input::gDoubleNoValue, J, slope.getRows() - 1,
                       slope.getCols() - 1);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void fill_output_vectors(double Dt,
                         double W,
                         Energy *egy,
                         Snow *snow,
                         Glacier *glac,
                         Water *wat,
                         Meteo *met,
                         Par *par,
                         Times *time,
                         Topo *top,
                         Soil *sl)
{
  long i, j, r = 0, c = 0;

  for (j = 1; j <= long(par->total_pixel); j++)
    {
      // TODO mattiu
      if (par->output_soil[geotop::common::Variables::i_sim] > 0)
        {
          r = top->rc_cont[j][1];
          c = top->rc_cont[j][2];
          if (geotop::common::Variables::files[fpnet] !=
              geotop::input::gStringNoValue)
            sl->Pnetcum[j] += wat->Pnet[r][c];
          if (geotop::common::Variables::files[fevap] !=
              geotop::input::gStringNoValue)
            {
              for (i = 1; i <= geotop::common::Variables::Nl; i++)
                {
                  sl->ETcum[j] += sl->ET[i][r][c];
                }
            }
        }
      if (par->output_snow[geotop::common::Variables::i_sim] > 0)
        {
          r = top->rc_cont[j][1];
          c = top->rc_cont[j][2];
          if (geotop::common::Variables::files[fHN] !=
              geotop::input::gStringNoValue)
            snow->HNcum[j] += wat->HN[r][c];
        }  // end mattiu

      if (par->output_snow[geotop::common::Variables::i_sim] > 0)
        {
          if (geotop::common::Variables::files[fsnowmelt] !=
              geotop::input::gStringNoValue)
            snow->MELTED[j] += snow->melted[j];
          if (geotop::common::Variables::files[fsnowsubl] !=
              geotop::input::gStringNoValue)
            snow->SUBL[j] += snow->subl[j];
          if (geotop::common::Variables::files[fsndur] !=
              geotop::input::gStringNoValue)
            {
              if (snow->yes[j] == 1) snow->t_snow[j] += Dt / GTConst::secinday;
            }
        }

      if (par->max_glac_layers > 0 &&
          par->output_glac[geotop::common::Variables::i_sim] > 0)
        {
          if (geotop::common::Variables::files[fglacmelt] !=
              geotop::input::gStringNoValue)
            glac->MELTED[j] += glac->melted[j];
          if (geotop::common::Variables::files[fglacsubl] !=
              geotop::input::gStringNoValue)
            glac->SUBL[j] += glac->subl[j];
        }

      if (par->output_surfenergy[geotop::common::Variables::i_sim] > 0)
        {
          if (geotop::common::Variables::files[fradnet] !=
              geotop::input::gStringNoValue)
            egy->Rn_mean[j] +=
              (egy->SW[j] + egy->LW[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fradLWin] !=
              geotop::input::gStringNoValue)
            egy->LWin_mean[j] +=
              (egy->LWin[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fradLW] !=
              geotop::input::gStringNoValue)
            egy->LW_mean[j] +=
              (egy->LW[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fradSW] !=
              geotop::input::gStringNoValue)
            egy->SW_mean[j] +=
              (egy->SW[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fradSWin] !=
              geotop::input::gStringNoValue)
            egy->Rswdown_mean[j] +=
              (egy->SWin[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fradSWinbeam] !=
              geotop::input::gStringNoValue)
            egy->Rswbeam_mean[j] +=
              (egy->SWinb[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fG] != geotop::input::gStringNoValue)
            egy->SEB_mean[j] +=
              (egy->G[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fH] != geotop::input::gStringNoValue)
            egy->H_mean[j] +=
              (egy->H[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fLE] !=
              geotop::input::gStringNoValue)
            egy->ET_mean[j] +=
              (egy->LE[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);
          if (geotop::common::Variables::files[fTs] !=
              geotop::input::gStringNoValue)
            egy->Ts_mean[j] +=
              (egy->Ts[j]) * Dt /
              (par->output_surfenergy[geotop::common::Variables::i_sim] * 3600.);

          if (geotop::common::Variables::files[fshadow] !=
              geotop::input::gStringNoValue)
            {
              if (egy->shad[j] >= 0) egy->nDt_sun[j]++;
              if (egy->shad[j] == 0) egy->nDt_shadow[j]++;
            }
        }
      if (par->output_meteo[geotop::common::Variables::i_sim] > 0)
        {
          if (geotop::common::Variables::files[fprec] !=
              geotop::input::gStringNoValue)
            {
              wat->PrTOT_mean[j] += wat->Pt[j];
              wat->PrSNW_mean[j] += wat->Ps[j];
            }
        }

      if (time->JD_plots.size() > 1 && W > 0)
        {
          if (geotop::common::Variables::files[pH] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pHg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->Hgplot[j] += egy->Hgp[j];
          if (geotop::common::Variables::files[pH] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pHv] !=
              geotop::input::gStringNoValue)
            egy->Hvplot[j] += egy->Hvp[j];
          if (geotop::common::Variables::files[pLE] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pLEg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->LEgplot[j] += egy->LEgp[j];
          if (geotop::common::Variables::files[pLE] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pLEv] !=
              geotop::input::gStringNoValue)
            egy->LEvplot[j] += egy->LEvp[j];
          if (geotop::common::Variables::files[pSWin] !=
              geotop::input::gStringNoValue)
            egy->SWinplot[j] += egy->SWinp[j];
          if (geotop::common::Variables::files[pSWg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->SWgplot[j] += egy->SWgp[j];
          if (geotop::common::Variables::files[pSWv] !=
              geotop::input::gStringNoValue)
            egy->SWvplot[j] += egy->SWvp[j];
          if (geotop::common::Variables::files[pLWin] !=
              geotop::input::gStringNoValue)
            egy->LWinplot[j] += egy->LWinp[j];
          if (geotop::common::Variables::files[pLWg] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pG] != geotop::input::gStringNoValue)
            egy->LWgplot[j] += egy->LWgp[j];
          if (geotop::common::Variables::files[pLWv] !=
              geotop::input::gStringNoValue)
            egy->LWvplot[j] += egy->LWvp[j];
          if (geotop::common::Variables::files[pTs] !=
              geotop::input::gStringNoValue)
            egy->Tsplot[j] += egy->Tsp[j];
          if (geotop::common::Variables::files[pTg] !=
              geotop::input::gStringNoValue)
            egy->Tgplot[j] += egy->Tgp[j];
          if (geotop::common::Variables::files[pD] != geotop::input::gStringNoValue)
            snow->Dplot[j] += W * DEPTH(top->rc_cont[j][1], top->rc_cont[j][2],
                                        snow->S->lnum, snow->S->Dzl);
          if (geotop::common::Variables::files[pTa] !=
              geotop::input::gStringNoValue)
            met->Taplot[j] +=
              W * met->Tgrid[top->rc_cont[j][1]][top->rc_cont[j][2]];
          if (geotop::common::Variables::files[pRH] !=
              geotop::input::gStringNoValue)
            met->RHplot[j] +=
              W * met->RHgrid[top->rc_cont[j][1]][top->rc_cont[j][2]];
          if (geotop::common::Variables::files[pVspd] !=
              geotop::input::gStringNoValue ||
              geotop::common::Variables::files[pVdir] !=
              geotop::input::gStringNoValue)
            {
              met->Vxplot[j] -=
                W * met->Vgrid[top->rc_cont[j][1]][top->rc_cont[j][2]] *
                sin(met->Vdir[top->rc_cont[j][1]][top->rc_cont[j][2]] * GTConst::Pi /
                    180.);
              met->Vyplot[j] -=
                W * met->Vgrid[top->rc_cont[j][1]][top->rc_cont[j][2]] *
                cos(met->Vdir[top->rc_cont[j][1]][top->rc_cont[j][2]] * GTConst::Pi /
                    180.);
            }
        }
      if (par->state_pixel == 1)
        {
          if (par->jplot[j] > 0 &&
              par->Dtplot_point[geotop::common::Variables::i_sim] > 0)
            {
              for (i = 0; i < otot; i++)
                {
                  geotop::common::Variables::odpnt[i][par->jplot[j] - 1] +=
                    geotop::common::Variables::odp[i][par->jplot[j] - 1];
                }
            }
        }
    }

  if (par->state_basin == 1)
    {
      if (par->Dtplot_basin[geotop::common::Variables::i_sim] > 0)
        {
          for (i = 0; i < ootot; i++)
            {
              geotop::common::Variables::odbsn[i] +=
                geotop::common::Variables::odb[i];
            }
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_snow_output(long IDpoint,
                       long iname,
                       long r,
                       long c,
                       double init_date,
                       double JDfrom0,
                       long day,
                       long month,
                       long year,
                       long hour,
                       long minute,
                       const GeoVector<double> &plot_depth,
                       Statevar3D *Snow,
                       Par *par,
                       double cosslope)
{
  std::string NNNN = "NNNN";
  string name, temp, temp2;

  std::string rec = "_recNNNN";
  std::string crec = "_crecNNNN";
  FILE *f;

  write_suffix(NNNN, iname, 0);

  if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

  //    "SnowTempProfileFile",//fsnTz
  if (geotop::common::Variables::files[fsnTz] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fsnTz] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      long actual_snow_layer_numb = Snow->lnum(r, c);
      long max_snow_layer_numb = par->max_snow_layers;
      write_snow_file(0, IDpoint, r, c, actual_snow_layer_numb,
                      max_snow_layer_numb, f, day, month, year, hour, minute,
                      JDfrom0, init_date, plot_depth, Snow->Dzl, Snow->T,
                      cosslope);
      fclose(f);
    }
  //    //    "SnowTempProfileFileWriteEnd",//fsnTzwriteend
  //    if (geotop::common::Variables::files[fsnTzwriteend] !=
  //    geotop::input::gStringNoValue){
  //    long actual_snow_layer_numb=Snow->lnum(r,c);
  //    long max_snow_layer_numb=par->max_snow_layers;
  //    write_snow_file(1, par->IDpoint[i], r, c,actual_snow_layer_numb,
  //max_snow_layer_numb, ffsnowT, day, month, year, hour, minute, JDfrom0,
  //init_date,
  //                    par->snow_plot_depths, Snow->Dzl,
  //Snow->T, cosslope);
  //        }
  //    Hack: writeEnd possibility not implemented. In case, the above is
  // reported a possible implementation. Consider to open the file connection
  // elsewhere

  //    "SnowLiqContentProfileFile",//fsnlz
  if (geotop::common::Variables::files[fsnlz] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fsnlz] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      long actual_snow_layer_numb = Snow->lnum(r, c);
      long max_snow_layer_numb = par->max_snow_layers;
      write_snow_file(1, IDpoint, r, c, actual_snow_layer_numb,
                      max_snow_layer_numb, f, day, month, year, hour, minute,
                      JDfrom0, init_date, plot_depth, Snow->Dzl, Snow->w_liq,
                      cosslope);
      fclose(f);
    }
  //    "SnowLiqContentProfileFileWriteEnd",//fsnlzwriteend NOT IMPLEMENTED

  //    "SnowIceContentProfileFile",//fsniz
  if (geotop::common::Variables::files[fsniz] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fsniz] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      long actual_snow_layer_numb = Snow->lnum(r, c);
      long max_snow_layer_numb = par->max_snow_layers;
      write_snow_file(1, IDpoint, r, c, actual_snow_layer_numb,
                      max_snow_layer_numb, f, day, month, year, hour, minute,
                      JDfrom0, init_date, plot_depth, Snow->Dzl, Snow->w_ice,
                      cosslope);
      fclose(f);
    }
  //    "SnowIceContentProfileFileWriteEnd",//fsnizwriteend NOT IMPLEMENTED

  //    "SnowDepthLayersFile",//fsndz
  if (geotop::common::Variables::files[fsndz] !=
      geotop::input::gStringNoValue)
    {
      temp = geotop::common::Variables::files[fsndz] + string(NNNN);

      if (par->n_ContRecovery > 0)
        {
          temp2 = temp + crec;
          name = temp2 + textfile;
        }
      else
        {
          name = temp + textfile;
        }

      f = fopen(name.c_str(), "a");
      long actual_snow_layer_numb = Snow->lnum(r, c);
      long max_snow_layer_numb = par->max_snow_layers;
      write_snow_file(2, IDpoint, r, c, actual_snow_layer_numb,
                      max_snow_layer_numb, f, day, month, year, hour, minute,
                      JDfrom0, init_date, plot_depth, Snow->Dzl, Snow->Dzl,
                      cosslope);
      fclose(f);
    }
  //    "SnowDepthLayersFileWriteEnd",//fsndzwriteend NOT IMPLEMENTED
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_snow_file(short choice,
                     long IDpoint,
                     long r,
                     long c,
                     long actual_snow_layer_numb,
                     long max_snow_layer_numb,
                     FILE *f,
                     long d,
                     long m,
                     long y,
                     long h,
                     long mi,
                     double JDfrom0,
                     double JDfrom0init,
                     const GeoVector<double> &plot_depth,
                     const GeoTensor<double> &dz,
                     const GeoTensor<double> &var_to_print,
                     double cosslope)
// choice=0 snow(according to snow_depth_plot) and var used
// choice=1 snow(according to snow_depth_plot) and var/dz used
// choice=2 snow(according to layers) and var used
// choice=3 snow(according to layers) and var/dz used
{
  short first_column = 1;
  long j;

  for (j = 0; j < geotop::common::Variables::nosnw; j++)
    {
      if (first_column == 0)
        {
          fprintf(f, ",");
        }
      else
        {
          first_column = 0;
        }
      if (geotop::common::Variables::osnw[j] >= 0)
        {
          if (geotop::common::Variables::osnw[j] == 0)
            {
              fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f", (float)d, (float)m,
                      (float)y, (float)h, (float)mi);
            }
          else if (geotop::common::Variables::osnw[j] == 1)
            {
              fprintf(f, "%f", JDfrom0);
            }
          else if (geotop::common::Variables::osnw[j] == 2)
            {
              fprintf(f, "%f", JDfrom0 - JDfrom0init);
            }
          else if (geotop::common::Variables::osnw[j] == 3)
            {
              fprintf(f, "%ld", geotop::common::Variables::i_sim);
            }
          else if (geotop::common::Variables::osnw[j] == 4)
            {
              // need to balance the braces
            }
          else if (geotop::common::Variables::osnw[j] == 5)
            {
              fprintf(f, "%ld", IDpoint);
            }
        }
      else
        {
          fprintf(f, "%f", geotop::input::gDoubleNoValue);
        }
    }

  if ((choice == 0 || choice == 1) &&
      (long)plot_depth[1] != geotop::input::gDoubleNoValue)
    {
      // prints output for selected snow depth
      for (size_t l = 1; l < plot_depth.size(); l++)
        {
          fprintf(
            f, ",%f",
            interpolate_snow(r, c, plot_depth[l] * cosslope, actual_snow_layer_numb,
                             dz, var_to_print, choice));
        }

    }
  else
    {
      // prints output for all snow layers
      if (choice == 2 || choice == 0)
        {
          // var_to_print is NOT to be divided by layer depth
          for (long l = 1; l <= (long)(max_snow_layer_numb); l++)
            {
              if (dz(l, r, c) > 0)
                {
                  fprintf(f, ",%f", var_to_print(l, r, c));
                }
              else
                {
                  fprintf(f, ",%f", (double)geotop::input::gDoubleNoValue);
                }
            }

        }
      else if (choice == 3 || choice == 1)
        {
          // var_to_print is to be divided by layer depth
          for (long l = 1; l <= (long)(max_snow_layer_numb); l++)
            {
              if (dz(l, r, c) > 0)
                {
                  fprintf(f, ",%f", var_to_print(l, r, c) / dz(l, r, c));
                }
              else
                {
                  fprintf(f, ",%f", (double)geotop::input::gDoubleNoValue);
                }
            }
        }
    }

  fprintf(f, "\n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_snow_header(FILE *f,
                       const GeoVector<double> &plot_depth,
                       size_t max_snow_layer)
{
  short first_column = 1;
  long j;
  size_t l;

  for (j = 0; j < geotop::common::Variables::nosnw; j++)
    {
      if (first_column == 0)
        {
          fprintf(f, ",");
        }
      else
        {
          first_column = 0;
        }
      if (geotop::common::Variables::osnw[j] >= 0 &&
          geotop::common::Variables::osnw[j] <= 5)
        {
          fprintf(
            f, "%s",
            geotop::common::Variables::hsnw[geotop::common::Variables::osnw[j]]
            .c_str());
        }
      else
        {
          fprintf(f, "none");
        }
    }

  if ((long)plot_depth[1] != geotop::input::gDoubleNoValue)
    {
      for (l = 1; l < plot_depth.size(); l++)
        {
          fprintf(f, ",%f", plot_depth[l]);
        }
    }
  else
    {
      for (l = 1; l <= max_snow_layer; l++)
        {
          fprintf(f, ",L%ld", l);
        }
    }
  fprintf(f, "\n");
}
