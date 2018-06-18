
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
#include "constants.h"
#include "snow.h"
#include "deallocate.h"
#include "init.h"

extern char **files;
extern FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq,
       *ffliqav, *ffice, *fficeav, *ffsnowT, *ffsnowl, *ffsnowi, *ffsnowd, *ffglac;
extern double **odpnt, * *odp, *odbsn, *odb;
extern long *opnt, nopnt, *obsn, nobsn, *osnw, nosnw;
extern long *oglc, noglc, *osl, nosl;
extern short *ipnt, *ibsn;
extern char **hpnt, * *hbsn, * *hsnw, * *hglc, * *hsl;
extern const char *WORKING_DIRECTORY;
extern char *string_novalue;
extern long number_novalue;
extern long Nl, Nr, Nc;
extern T_INIT *UV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,
                 PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times)
{

  long i,j,r,l,n;

  printf("Close files\n");
  if (strcmp(files[fpointwriteend], string_novalue) != 0) fclose(ffpoint);
  if (strcmp(files[fsnTzwriteend], string_novalue) != 0) fclose(ffsnowT);
  if (strcmp(files[fsnlzwriteend], string_novalue) != 0) fclose(ffsnowl);
  if (strcmp(files[fsnizwriteend], string_novalue) != 0) fclose(ffsnowi);
  if (strcmp(files[fsndzwriteend], string_novalue) != 0) fclose(ffsnowd);
  if (strcmp(files[fglzwriteend], string_novalue) != 0) fclose(ffglac);
  if (strcmp(files[fbaswriteend], string_novalue) != 0) fclose(ffbas);
  if (strcmp(files[fTzwriteend], string_novalue) != 0) fclose(ffT);
  if (strcmp(files[fTzavwriteend], string_novalue) != 0) fclose(ffTav);
  if (strcmp(files[fpsizwriteend], string_novalue) != 0) fclose(ffpsi);
  if (strcmp(files[fpsiztotwriteend], string_novalue) != 0) fclose(ffpsitot);
  if (strcmp(files[fliqzwriteend], string_novalue) != 0) fclose(ffliq);
  if (strcmp(files[ficezwriteend], string_novalue) != 0) fclose(ffice);

  printf("Deallocating global variables\n");
  if (par->state_pixel == 1)
    {
      for (i=0; i<otot; i++)
        {
          free(odpnt[i]);
          free(odp[i]);
        }
      free(odpnt);
      free(odp);
    }
  for (i=0; i<otot; i++)
    {
      free(hpnt[i]);
    }
  free(hpnt);
  free(opnt);
  free(ipnt);

  for (i=0; i<ootot; i++)
    {
      free(hbsn[i]);
    }
  free(odbsn);
  free(odb);
  free(hbsn);
  free(obsn);
  free(ibsn);

  for (j=0; j<6; j++)
    {
      free(hsnw[j]);
    }

  for (j=0; j<10; j++)
    {
      free(hglc[j]);
    }
  free(hsnw);
  free(hglc);

  free(osnw);
  free(oglc);

  for (j=0; j<6; j++)
    {
      free(hsl[j]);
    }
  free(hsl);
  free(osl);

  /* Deallocation of struct SOIL "sl": */
  printf("Deallocating soil\n");
  free_doublematrix(sl->Ptot);
  if (par->output_soil_bin == 1)
    {
      if (strcmp(files[fTav], string_novalue) != 0
          || strcmp(files[fTavsup],
                    string_novalue) != 0) free_doublematrix(sl->T_av_tensor);
      if (strcmp(files[fliqav],
                 string_novalue) != 0) free_doublematrix(sl->thw_av_tensor);
      if (strcmp(files[ficeav],
                 string_novalue) != 0) free_doublematrix(sl->thi_av_tensor);
    }
  free_doublematrix(sl->th);
  free_longmatrix(sl->type);
  free_doubletensor(sl->pa);
  free_doubletensor(sl->ET);
  delete sl->SS;

  if (par->state_pixel == 1)
    {
      if (strcmp(files[fTz], string_novalue) != 0
          || strcmp(files[fTzwriteend],
                    string_novalue) != 0) free_doublematrix(sl->Tzplot);
      if (strcmp(files[fTzav], string_novalue) != 0
          || strcmp(files[fTzavwriteend],
                    string_novalue) != 0) free_doublematrix(sl->Tzavplot);
      if (strcmp(files[fpsiztot], string_novalue) != 0
          || strcmp(files[fpsiztotwriteend],
                    string_novalue) != 0) free_doublematrix(sl->Ptotzplot);
      if (strcmp(files[fpsiz], string_novalue) != 0
          || strcmp(files[fpsizwriteend],
                    string_novalue) != 0) free_doublematrix(sl->Pzplot);
      if (strcmp(files[fliqz], string_novalue) != 0
          || strcmp(files[fliqzwriteend],
                    string_novalue) != 0) free_doublematrix(sl->thzplot);
      if (strcmp(files[fliqzav], string_novalue) != 0
          || strcmp(files[fliqzavwriteend],
                    string_novalue) != 0) free_doublematrix(sl->thzavplot);
      if (strcmp(files[ficez], string_novalue) != 0
          || strcmp(files[ficezwriteend],
                    string_novalue) != 0) free_doublematrix(sl->thizplot);
      if (strcmp(files[ficezav], string_novalue) != 0
          || strcmp(files[ficezavwriteend],
                    string_novalue) != 0) free_doublematrix(sl->thizavplot);
      if (strcmp(files[fsatz],
                 string_novalue) != 0) free_doublematrix(sl->satratio);

      if (strcmp(files[fTrun], string_novalue) != 0) free_doublematrix(sl->Tzrun);
      if (strcmp(files[fwrun], string_novalue) != 0) free_doublematrix(sl->wzrun);
      if (strcmp(files[fdUrun], string_novalue) != 0) free_doublematrix(sl->dUzrun);
      if (strcmp(files[fSWErun],
                 string_novalue) != 0) free_doublematrix(sl->SWErun);
      if (strcmp(files[fTmaxrun],
                 string_novalue) != 0) free_doublematrix(sl->Tzmaxrun);
      if (strcmp(files[fTminrun],
                 string_novalue) != 0) free_doublematrix(sl->Tzminrun);
      if (strcmp(files[fwmaxrun],
                 string_novalue) != 0) free_doublematrix(sl->wzmaxrun);
      if (strcmp(files[fwminrun],
                 string_novalue) != 0) free_doublematrix(sl->wzminrun);

    }

  //  free(sl);

  /* Deallocation of struct TOPO "top": */
  printf("Deallocating top\n");

  if (par->point_sim==1)
    {
      for (i=1; i<=top->num_horizon_point; i++)
        {
          for (j=0; j<top->horizon_numlines[i-1]; j++)
            {
              free(top->horizon_height[i-1][j]);
            }
          free(top->horizon_height[i-1]);
        }
      free(top->horizon_height);
      free(top->horizon_numlines);

      free_longmatrix(top->horizon_point);
    }

  free_doublematrix(top->sky);
  free_doublematrix(top->Z0);
  free_doublematrix(top->East);
  free_doublematrix(top->North);
  free_shortmatrix(top->pixel_type);
  free_doublematrix(top->aspect);
  free_doublematrix(top->slope);
  free_doublematrix(top->dzdE);
  free_doublematrix(top->dzdN);
  free_doublematrix(top->curvature1);
  free_doublematrix(top->curvature2);
  free_doublematrix(top->curvature3);
  free_doublematrix(top->curvature4);
  free_longmatrix(top->Jdown);
  free_doublematrix(top->Qdown);
  if (par->point_sim==0) free_shortmatrix(top->is_on_border);

  free_longmatrix(top->lrc_cont);

  n = Fminlong(par->Nl_spinup->co[par->init_date->nh],Nl);
  for (l=0; l<=n; l++)
    {
      for (r=1; r<=Nr; r++)
        {
          free(top->i_cont[l][r]);
        }
      free(top->i_cont[l]);
    }
  free(top->i_cont);

  free_longmatrix(top->rc_cont);
  for (r=1; r<=Nr; r++)
    {
      free(top->j_cont[r]);
    }
  free(top->j_cont);

  free_doubletensor(top->Z);
  
  free_longmatrix(top->BC_counter);

  if (par->point_sim==1)
    {
      free_doublematrix(top->latitude);
      free_doublematrix(top->longitude);
    }

  //  free(top);


  /* Deallocation of struct LAND "land": */
  printf("Deallocating land\n");
  free_doublematrix(land->LC);
  free_doublematrix(land->delay);
  free_shortmatrix(land->shadow);
  free_doublematrix(land->ty);
  free_doublematrix(land->root_fraction);

  for (i=0; i<par->n_landuses; i++)
    {
      free(land->vegparv[i]);

      if (par->vegflag->co[i+1]==1)
        {
          for (j=0; j<land->NumlinesVegTimeDepData[i]; j++)
            {
              free(land->vegpars[i][j]);
            }
          free(land->vegpars[i]);
        }
    }

  free(land->vegpars);
  free(land->vegparv);
  free(land->NumlinesVegTimeDepData);

  //  free(land);

  /* Deallocation of struct WATER "water": */
  printf("Deallocating water\n");
  free_doublematrix(wat->PrecTot);
  free_doublematrix(wat->Pnet);

  //  free(wat);

  /* Deallocation of struct CHANNEL "channel": */
  printf("Deallocating channel network\n");
  free_longvector(cnet->r);
  free_longvector(cnet->c);
  free_longmatrix(cnet->ch);
  free_longvector(cnet->ch_down);
  for (l=0; l<=Nl; l++)
    {
      free(cnet->ch3[l]);
    }
  free(cnet->ch3);
  free_longmatrix(cnet->lch);
  free_longvector(cnet->soil_type);
  free_doublematrix(cnet->th);
  free_doublematrix(cnet->ET);
  delete cnet->SS;
  //  free(cnet);

  /* Deallocation of struct T_INIT "UV": */
  printf("Deallocating UV\n");

  /* Deallocation of struct ENERGY "egy": */
  printf("Deallocating egy\n");
  if (par->output_surfenergy_bin == 1)
    {
      if (strcmp(files[fshadow], string_novalue) != 0)
        {
          free_shortvector(egy->shad);
        }
    }

  free(egy->sun);

  if (times->JD_plots->nh > 1)
    {
      if (strcmp(files[pH], string_novalue) != 0
          || strcmp(files[pHg], string_novalue) != 0
          || strcmp(files[pG], string_novalue) != 0)
        {
        }
      if (strcmp(files[pLE], string_novalue) != 0
          || strcmp(files[pLEg], string_novalue) != 0
          || strcmp(files[pG], string_novalue) != 0)
        {
        }
      if (strcmp(files[pH], string_novalue) != 0
          || strcmp(files[pHv], string_novalue) != 0)
        {
        }
      if (strcmp(files[pLE], string_novalue) != 0
          || strcmp(files[pLEv], string_novalue) != 0)
        {
        }
      if (strcmp(files[pSWin], string_novalue) != 0)
        {
        }
      if (strcmp(files[pSWg], string_novalue) != 0
          || strcmp(files[pG], string_novalue) != 0)
        {
        }
      if (strcmp(files[pSWv], string_novalue) != 0)
        {
        }
      if (strcmp(files[pLWin], string_novalue) != 0)
        {
        }
      if (strcmp(files[pLWg], string_novalue) != 0
          || strcmp(files[pG], string_novalue) != 0)
        {
        }
      if (strcmp(files[pLWv], string_novalue) != 0)
        {
        }
      if (strcmp(files[pTs], string_novalue) != 0)
        {
        }
      if (strcmp(files[pTg], string_novalue) != 0)
        {
        }
    }

  free_doublematrix(egy->Tgskin_surr);
  free_doublematrix(egy->SWrefl_surr);

  //  free(egy);


  /* Deallocation of struct SNOW "snow": */
  printf("Deallocating snow\n");

  if (times->JD_plots->nh > 1)
    {
    }

  delete snow->S;

  if (par->blowing_snow==1)
    {
      deallocate_statevar_1D(snow->S_for_BS);
      free_doublematrix(snow->Nabla2_Qtrans);
      free_doublematrix(snow->Qsub);
      free_doublematrix(snow->Qsub_x);
      free_doublematrix(snow->Qsub_y);
      free_doublematrix(snow->Qtrans);
      free_doublematrix(snow->Qtrans_x);
      free_doublematrix(snow->Qtrans_y);
      free_doublematrix(snow->Qsalt);

      if (par->output_snow_bin == 1)
        {
          free_doublematrix(snow->Wtrans_plot);
          free_doublematrix(snow->Wsubl_plot);
        }
    }

  if (par->output_snow_bin == 1)
    {
      if (strcmp(files[fsndur], string_novalue) != 0)
        {
          free_shortvector(snow->yes);
        }
      if (strcmp(files[fsnowmelt], string_novalue) != 0)
        {
        }
      if (strcmp(files[fsnowsubl], string_novalue) != 0)
        {
        }
    }

  if (par->blowing_snow==1) free_longvector(snow->change_dir_wind);
  //  free(snow);

  printf("Deallocating glacier\n");
  if (par->max_glac_layers>0)
    {
      delete glac->G;
      if (par->output_glac_bin == 1)
        {
          if (strcmp(files[fglacmelt], string_novalue) != 0)
            {
            }
          if (strcmp(files[fglacsubl], string_novalue) != 0)
            {
            }
        }
    }
  //  free(glac);

  printf("Deallocating met\n");

  free_doublematrix(met->Tgrid);
  free_doublematrix(met->Pgrid);
  free_doublematrix(met->Vgrid);
  free_doublematrix(met->Vdir);
  free_doublematrix(met->RHgrid);

  if (times->JD_plots->nh > 1)
    {
      if (strcmp(files[pVspd], string_novalue) != 0
          || strcmp(files[pVdir], string_novalue) != 0)
        {
        }
    }

  for (i=0; i<met->st->Z->nh; i++)
    {

      for (j=0; j<met->numlines[i]; j++)
        {
          free(met->data[i][j]);
        }
      free(met->data[i]);

      free(met->var[i]);

      for (j=0; j<met->horizonlines[i]; j++)
        {
          free(met->horizon[i][j]);
        }
      free(met->horizon[i]);
    }

  free(met->data);
  free(met->numlines);
  free(met->var);
  free(met->line_interp_WEB);
  free(met->line_interp_Bsnow);
  free(met->horizon);
  free(met->horizonlines);

  if (par->LRflag==1)
    {
      for (i=0; i<met->LRsnr; i++)
        {
          free(met->LRs[i]);
        }
      free(met->LRs);
    }
  free(met->LRv);

  for (i=0; i<nlstot; i++)
    {
      free(met->LRc[i]);
    }
  free(met->LRc);
  free(met->LRcnc);
  free(met->LRd);

  free_longvector(met->imeteo_stations);
//  dealloc_meteostations(met->st);

  free(met->qinv);
  if (par->qin == 1)
    {
      for (i=0; i<met->qinsnr; i++)
        {
          free(met->qins[i]);
        }
    }

  //  free(met);

  printf("Deallocating times\n");
  free(times->Dt_vector);
  if (par->tsteps_from_file==1)
    {
      for (j=0; j<times->numlinesDt_matrix; j++)
        {
          free(times->Dt_matrix[j]);
        }
      free(times->Dt_matrix);
    }
  //  free(times);

  /* Deallocation of struct PAR "par": */
  printf("Deallocating par\n");
  if (par->state_pixel == 1)
    {
      free_longvector(par->jplot);
      free_longmatrix(par->rc);
      free_longvector(par->IDpoint);
    }


  free_longvector(par->run_times);

  if (par->point_sim == 1) free_doublematrix(par->maxSWE);

  free_shortvector(par->plot_discharge_with_Dt_integration);
  free_shortvector(par->plot_point_with_Dt_integration);
  free_shortvector(par->plot_basin_with_Dt_integration);




  free_shortvector(par->linear_interpolation_meteo);

  free_longvector(par->inf_snow_layers);
  free_longvector(par->inf_glac_layers);

  free_longvector(par->Nl_spinup);

  //  free(par);

  /* Deallocation of struct FILENAMES "filenames": */
  printf("Deallocating files\n");
  for (i=0; i<nfiles; i++)
    {
      free(files[i]);
    }
  free(files);

  printf("Deallocating novalues\n");
  free(string_novalue);


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

void reset_to_zero(PAR *par, SOIL *sl, LAND *land, SNOW *snow, GLACIER *glac,
                   ENERGY *egy, METEO *met, WATER *wat)
{

  long i, j;

  if (par->state_pixel == 1)
    {
      if (strcmp(files[fTzav], string_novalue) != 0
          || strcmp(files[fTzavwriteend],
                    string_novalue) != 0) initialize_doublematrix(sl->Tzavplot,0.);
      if (strcmp(files[fliqzav], string_novalue) != 0
          || strcmp(files[fliqzavwriteend],
                    string_novalue) != 0) initialize_doublematrix(sl->thzavplot,0.);
      if (strcmp(files[ficezav], string_novalue) != 0
          || strcmp(files[ficezavwriteend],
                    string_novalue) != 0) initialize_doublematrix(sl->thizavplot,0.);

      for (i=1; i<=par->rc->nrh; i++)
        {
          for (j=0; j<otot; j++)
            {
              odpnt[j][i-1]=0.0;
            }
        }
    }

  if (par->state_basin== 1)
    {
      for (j=0; j<ootot; j++)
        {
          odbsn[j]=0.0;
        }
    }

  if (par->output_soil_bin == 1)
    {
      if (strcmp(files[fTav], string_novalue) != 0
          || strcmp(files[fTavsup],
                    string_novalue) != 0) initialize_doublematrix(sl->T_av_tensor, 0.);
      if (strcmp(files[ficeav],
                 string_novalue) != 0) initialize_doublematrix(sl->thi_av_tensor, 0.);
      if (strcmp(files[fliqav],
                 string_novalue) != 0) initialize_doublematrix(sl->thw_av_tensor, 0.);
    }

  if (par->output_snow_bin == 1)
    {

      if (strcmp(files[fswe], string_novalue) != 0 && par->blowing_snow==1)
        {
          initmatrix(0.0, snow->Wtrans_plot, land->LC, number_novalue);
          initmatrix(0.0, snow->Wsubl_plot, land->LC, number_novalue);
        }
    }

  if (par->output_surfenergy_bin == 1)
    {

      if (strcmp(files[fshadow], string_novalue) != 0)
        {
          (*egy->nDt_shadow) = 0;
          (*egy->nDt_sun) = 0;
        }

    }


}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

