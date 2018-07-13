
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
#include <logger.h>

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
                 PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times) {

    long i, j, r, l, n;

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

    geolog << "Deallocating global variables" << std::endl;
    if (par->state_pixel == 1) {
        for (i = 0; i < otot; i++) {
            free(odpnt[i]);
            free(odp[i]);
        }
        free(odpnt);
        free(odp);
    }
    for (i = 0; i < otot; i++) {
        free(hpnt[i]);
    }
    free(hpnt);
    free(opnt);
    free(ipnt);

    for (i = 0; i < ootot; i++) {
        free(hbsn[i]);
    }
    free(odbsn);
    free(odb);
    free(hbsn);
    free(obsn);
    free(ibsn);

    for (j = 0; j < 6; j++) {
        free(hsnw[j]);
    }

    for (j = 0; j < 10; j++) {
        free(hglc[j]);
    }
    free(hsnw);
    free(hglc);

    free(osnw);
    free(oglc);

    for (j = 0; j < 6; j++) {
        free(hsl[j]);
    }
    free(hsl);
    free(osl);

    /* Deallocation of struct SOIL "sl": */
    geolog << "Deallocating soil" << std::endl;

    free_doubletensor(sl->pa);
    free_doubletensor(sl->ET);
    delete sl->SS;

    if (par->state_pixel == 1)
    {
    }

    //  free(sl);

    /* Deallocation of struct TOPO "top": */
    geolog << "Deallocating top" << std::endl;

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


    free_shortmatrix(top->pixel_type);
    free_longmatrix(top->Jdown);
    if (par->point_sim==0) free_shortmatrix(top->is_on_border);


    n = Fminlong((*par->Nl_spinup)(par->init_date->nh),Nl);
    for (l=0; l<=n; l++)
    {
        for (r=1; r<=Nr; r++)
        {
            free(top->i_cont[l][r]);
        }
        free(top->i_cont[l]);
    }
    free(top->i_cont);

    for (r=1; r<=Nr; r++)
    {
        free(top->j_cont[r]);
    }
    free(top->j_cont);

    free_doubletensor(top->Z);


    if (par->point_sim==1)
    {
    }

    //  free(top);


    /* Deallocation of struct LAND "land": */
    geolog << "Deallocating land" << std::endl;
    free_shortmatrix(land->shadow);

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
    geolog << "Deallocating water" << std::endl;

    //  free(wat);

    /* Deallocation of struct CHANNEL "channel": */
    geolog << "Deallocating channel network" << std::endl;
    free_longmatrix(cnet->ch);
    for (l=0; l<=Nl; l++)
    {
        free(cnet->ch3[l]);
    }
    free(cnet->ch3);
    free_longmatrix(cnet->lch);
    delete cnet->SS;
    //  free(cnet);

    /* Deallocation of struct T_INIT "UV": */
    geolog << "Deallocating UV" << std::endl;

    /* Deallocation of struct ENERGY "egy": */
    geolog << "Deallocating egy" << std::endl;

    if (par->output_surfenergy_bin == 1)
    {
        if (strcmp(files[fshadow], string_novalue) != 0)
        {
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

    //  free(egy);

    /* Deallocation of struct SNOW "snow": */
    geolog << "Deallocating snow" << std::endl;

    if (times->JD_plots->nh > 1)
    {
    }

    delete snow->S;

    if (par->blowing_snow==1)
    {
        deallocate_statevar_1D(snow->S_for_BS);

        if (par->output_snow_bin == 1)
        {
        }
    }

    if (par->output_snow_bin == 1)
    {
        if (strcmp(files[fsndur], string_novalue) != 0)
        {
        }
        if (strcmp(files[fsnowmelt], string_novalue) != 0)
        {
        }
        if (strcmp(files[fsnowsubl], string_novalue) != 0)
        {
        }
    }

    //  free(snow);

    geolog << "Deallocating glacier" << std::endl;
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

    geolog << "Deallocating met" << std::endl;

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


    geolog << "Deallocating times" << std::endl;
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
    geolog << "Deallocating par" << std::endl;
    if (par->state_pixel == 1)
    {
        free_longmatrix(par->rc);
    }

    if (par->point_sim == 1)
    {
    }

    //  free(par);

    /* Deallocation of struct FILENAMES "filenames": */
    geolog << "Deallocating files" << std::endl;
    for (i=0; i<nfiles; i++)
    {
        free(files[i]);
    }
    free(files);

    geolog << "Deallocating novalues" << std::endl;
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


void reset_to_zero(PAR *par, SOIL *sl, LAND *land, SNOW *snow, GLACIER *glac, ENERGY *egy, METEO *met, WATER *wat){

    long i, j;

    if(par->state_pixel == 1){
        if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0)
            (*sl->Tzavplot) = 0.;
        if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0)
            (*sl->thzavplot) = 0.;
        if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0)
            (*sl->thizavplot) = 0.;

        for(i=1;i<=par->rc->nrh;i++){
            for(j=0;j<otot;j++) {
                odpnt[j][i-1]=0.0;
            }
        }
    }

    if(par->state_basin== 1){
        for(j=0;j<ootot;j++){
            odbsn[j]=0.0;
        }
    }

    if(par->output_soil_bin == 1){
        if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0)
        (*sl->T_av_tensor) = 0.;

        if(strcmp(files[ficeav] , string_novalue) != 0)
            (*sl->thi_av_tensor) = 0.;
        if(strcmp(files[fliqav] , string_novalue) != 0)
            (*sl->thw_av_tensor) = 0.;
        if(strcmp(files[fpnet] , string_novalue) != 0)
            *(sl->Pnetcum) = 0;
        if(strcmp(files[fevap] , string_novalue) != 0)
            *(sl->ETcum) = 0;
    }

    if(par->output_snow_bin == 1){
        if(strcmp(files[fsnowmelt] , string_novalue) != 0) *(snow->MELTED)= 0.;
        if(strcmp(files[fsnowsubl] , string_novalue) != 0) *(snow->SUBL) = 0;
        if(strcmp(files[fswe] , string_novalue) != 0 && par->blowing_snow==1){
            initmatrix(0.0, snow->Wtrans_plot.get(), land->LC.get(), number_novalue);
            initmatrix(0.0, snow->Wsubl_plot.get(), land->LC.get(), number_novalue);
        }
        if(strcmp(files[fsndur] , string_novalue) != 0) *(snow->t_snow) = 0;
    }

    if(par->max_glac_layers>0 && par->output_glac_bin==1){
        if(strcmp(files[fglacmelt] , string_novalue) != 0) *(glac->MELTED) = 0;
        if(strcmp(files[fglacsubl] , string_novalue) != 0) *(glac->SUBL) = 0;
    }

    if(par->output_surfenergy_bin == 1){

        if(strcmp(files[fradnet] , string_novalue) != 0) *(egy->Rn_mean) = 0;
        if(strcmp(files[fradLWin] , string_novalue) != 0) *(egy->LWin_mean) = 0;
        if(strcmp(files[fradLW] , string_novalue) != 0) *(egy->LW_mean) = 0;
        if(strcmp(files[fradSW] , string_novalue) != 0) *(egy->SW_mean) = 0;
        if(strcmp(files[fradSWin] , string_novalue) != 0) *(egy->Rswdown_mean) = 0;
        if(strcmp(files[fradSWinbeam] , string_novalue) != 0) *(egy->Rswbeam_mean) = 0;
        if(strcmp(files[fshadow] , string_novalue) != 0){
            *(egy->nDt_shadow) = 0;
            *(egy->nDt_sun) = 0;
        }

        if(strcmp(files[fG] , string_novalue) != 0) *(egy->SEB_mean) = 0;
        if(strcmp(files[fH] , string_novalue) != 0) *(egy->H_mean) = 0;
        if(strcmp(files[fLE] , string_novalue) != 0) *(egy->ET_mean) = 0;
        if(strcmp(files[fTs] , string_novalue) != 0) *(egy->Ts_mean) = 0;
    }

    if (par->output_meteo_bin == 1){
        if(strcmp(files[fTa] , string_novalue) != 0) *(met->Tamean) = 0;
        if(strcmp(files[fwspd] , string_novalue) != 0) *(met->Vspdmean) = 0;
        if(strcmp(files[fwdir] , string_novalue) != 0) *(met->Vdirmean) = 0;
        if(strcmp(files[frh] , string_novalue) != 0) *(met->RHmean) = 0;
        if(strcmp(files[fprec] , string_novalue) != 0){
            *(wat->PrTOT_mean) = 0;
            *(wat->PrSNW_mean) = 0;
        }
    }

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

