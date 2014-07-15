
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.0.0 - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 2.0.0 
 
 GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#include "deallocate.h"
#include "geotop_common.h"
#include "inputKeywords.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

  void dealloc_all(Topo *top,Soil *sl,Land *land,Water *wat,Channel *cnet,Par *par,Energy *egy,Snow *snow, Glacier *glac, Meteo *met, Times *times){

	long i,j,r,l;
	
	printf("Close files\n");
    if(geotop::common::Variables::files[fpointwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffpoint)
        fclose(geotop::common::Variables::ffpoint);
    if(geotop::common::Variables::files[fsnTzwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffsnowT)
        fclose(geotop::common::Variables::ffsnowT);
    if(geotop::common::Variables::files[fsnlzwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffsnowl)
        fclose(geotop::common::Variables::ffsnowl);
    if(geotop::common::Variables::files[fsnizwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffsnowi)
        fclose(geotop::common::Variables::ffsnowi);
    if(geotop::common::Variables::files[fsndzwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffsnowd)
        fclose(geotop::common::Variables::ffsnowd);
    if(geotop::common::Variables::files[fglzwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffglac)
        fclose(geotop::common::Variables::ffglac);
    if(geotop::common::Variables::files[fbaswriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffbas)
        fclose(geotop::common::Variables::ffbas);
    if(geotop::common::Variables::files[fTzwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffT)
        fclose(geotop::common::Variables::ffT);
    if(geotop::common::Variables::files[fTzavwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffTav)
        fclose(geotop::common::Variables::ffTav);
    if(geotop::common::Variables::files[fpsizwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffpsi)
        fclose(geotop::common::Variables::ffpsi);
    if(geotop::common::Variables::files[fpsiztotwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffpsitot)
        fclose(geotop::common::Variables::ffpsitot);
    if(geotop::common::Variables::files[fliqzwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffliq)
        fclose(geotop::common::Variables::ffliq);
    if(geotop::common::Variables::files[ficezwriteend].compare(geotop::input::gStringNoValue) != 0
        && geotop::common::Variables::ffice)
        fclose(geotop::common::Variables::ffice);
	
	printf("Deallocating global variables\n"); 
	if(par->state_pixel==1){
		for (i=0; i<otot; i++) {
			free(geotop::common::Variables::odpnt[i]);
			free(geotop::common::Variables::odp[i]);
		}
		free(geotop::common::Variables::odpnt);
		free(geotop::common::Variables::odp);
	}

    free(geotop::common::Variables::opnt);
	free(geotop::common::Variables::ipnt);
	
	free(geotop::common::Variables::odbsn);
	free(geotop::common::Variables::odb);
	free(geotop::common::Variables::obsn);
	free(geotop::common::Variables::ibsn);
 	
	/* Deallocation of struct SOIL "sl": */
	printf("Deallocating soil\n");
	
	/* Deallocation of struct TOPO "top": */
	printf("Deallocating top\n");
	
	if(par->point_sim==1){
		for (i=1; i<=top->num_horizon_point; i++) {
			for (j=0; j<top->horizon_numlines[i-1]; j++) {
				free(top->horizon_height[i-1][j]);
			}
			free(top->horizon_height[i-1]);
		}
		free(top->horizon_height);
		free(top->horizon_numlines);
		
	}
	
	for(l=0;l<=geotop::common::Variables::Nl;l++){
		for(r=1;r<=geotop::common::Variables::Nr;r++){
			free(top->i_cont[l][r]);
		}
		free(top->i_cont[l]);
	}
	free(top->i_cont);
	
	for(r=1;r<=geotop::common::Variables::Nr;r++){
		free(top->j_cont[r]);
	}
	free(top->j_cont);
		
	if (par->point_sim==1) {
	}
	
	/* Deallocation of struct LAND "land": */
	printf("Deallocating land\n");
	
	for(size_t i=0;i<par->n_landuses;i++){
		free(land->vegparv[i]);
		
		if(par->vegflag[i+1]==1){
			for(j=0;j<land->NumlinesVegTimeDepData[i];j++){
				free(land->vegpars[i][j]);
			}
			free(land->vegpars[i]);
		}
	}
	
	free(land->vegpars);
	free(land->vegparv);
	free(land->NumlinesVegTimeDepData);
	
	/* Deallocation of struct WATER "water": */
	printf("Deallocating water\n"); 

	free(wat);
	
	/* Deallocation of struct CHANNEL "channel": */
	printf("Deallocating channel network\n"); 
	for (l=0; l<=geotop::common::Variables::Nl; l++) {
		free(cnet->ch3[l]);
	}
	free(cnet->ch3);
	
	/* Deallocation of struct T_INIT "UV": */
	printf("Deallocating UV\n"); 
	
	/* Deallocation of struct ENERGY "egy": */
	printf("Deallocating egy\n");  	
	
	free(egy->sun);
	
	/* Deallocation of struct SNOW "snow": */
	printf("Deallocating snow\n"); 

	if(par->blowing_snow==1){
		deallocate_statevar_1D(snow->S_for_BS);
	}
	
	if(par->blowing_snow==1)
        free(snow);
	
	printf("Deallocating glacier\n");
	if(par->max_glac_layers>0){ deallocate_statevar_3D(glac->G);
	}
	
	printf("Deallocating met\n"); 
	
    for(i=0;i<long(met->st->Z.size()-1);i++){

        for(j=0;j<met->horizonlines[i];j++){
            free(met->horizon[i][j]);
        }
        free(met->horizon[i]);
#ifdef USE_INTERNAL_METEODISTR
        free(met->var[i]);
#endif
    }
#ifdef USE_INTERNAL_METEODISTR
	free(met->var);
#endif
	free(met->data);
	free(met->numlines);
	free(met->horizon);
	free(met->horizonlines);
	
	if(par->LRflag==1){
		for (i=0; i<met->LRsnr; i++) {
			free(met->LRs[i]);
		}
		free(met->LRs);
	}
	free(met->LRv);
	
	for (i=0; i<nlstot; i++) {
		free(met->LRc[i]);
	}
	free(met->LRc);
	free(met->LRcnc);
	free(met->LRd);
	
	dealloc_meteostations(met->st); 
	
	free(met->qinv);
	if (par->qin == 1) {
		for (i=0; i<met->qinsnr; i++) {
			free(met->qins[i]);
		}
	}
	
	printf("Deallocating times\n");
	free(times);
	
	/* Deallocation of struct PAR "par": */
	printf("Deallocating par\n"); 
	
	/* Deallocation of struct FILENAMES "filenames": */
	
	printf("Deallocating novalues\n"); 
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

  void dealloc_meteostations(MeteoStations *st){
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_soil_state(SoilState *S){
	
	free(S);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_veg_state(StateVeg *V){

	free(V);
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

  void reset_to_zero(Par *par, Soil *sl, Land *land, Snow *snow, Glacier *glac, Energy *egy, Meteo *met, Water *wat){

	long i, j;
	
	if(par->state_pixel == 1){
		if(geotop::common::Variables::files[fTzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTzavwriteend] != geotop::input::gStringNoValue) {
			sl->Tzavplot.resize(sl->Tzavplot.getRows(), sl->Tzavplot.getCols(),0.0);
		}

		if(geotop::common::Variables::files[fliqzav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fliqzavwriteend] != geotop::input::gStringNoValue) {
			sl->thzavplot.resize(sl->thzavplot.getRows(),sl->thzavplot.getCols(),0.0);
		}
		if(geotop::common::Variables::files[ficezav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[ficezavwriteend] != geotop::input::gStringNoValue) {
			sl->thizavplot.resize(sl->thizavplot.getRows(),sl->thizavplot.getCols(),0.0);
		}
		
		for(i=1;i<par->rc.getRows();i++){
			for(j=0;j<otot;j++) { 
				geotop::common::Variables::odpnt[j][i-1]=0.0;
			}
		}	
	}
	
	if(par->state_basin== 1){
		for(j=0;j<ootot;j++){ 
			geotop::common::Variables::odbsn[j]=0.0; 
		}
	}
	
	if(par->output_soil_bin == 1){
		if(geotop::common::Variables::files[fTav] != geotop::input::gStringNoValue ||geotop::common::Variables::files[fTavsup] != geotop::input::gStringNoValue) sl->T_av_tensor.resize(sl->T_av_tensor.getRows(),sl->T_av_tensor.getCols(), 0.);
		if(geotop::common::Variables::files[ficeav] != geotop::input::gStringNoValue) sl->thi_av_tensor.resize(sl->thi_av_tensor.getRows(),sl->thi_av_tensor.getCols(), 0.0);
		if(geotop::common::Variables::files[fliqav] != geotop::input::gStringNoValue) sl->thw_av_tensor.resize(sl->thw_av_tensor.getRows(),sl->thw_av_tensor.getCols(), 0.0);
	}
	
	if(par->output_snow_bin == 1){
		if(geotop::common::Variables::files[fsnowmelt] != geotop::input::gStringNoValue) snow->MELTED.resize(snow->MELTED.size(), 0.);
		if(geotop::common::Variables::files[fsnowsubl] != geotop::input::gStringNoValue) snow->SUBL.resize(snow->SUBL.size(), 0.0);
		if(geotop::common::Variables::files[fswe] != geotop::input::gStringNoValue && par->blowing_snow==1){
			initmatrix(0.0, snow->Wtrans_plot, land->LC, geotop::input::gDoubleNoValue);
			initmatrix(0.0, snow->Wsubl_plot, land->LC, geotop::input::gDoubleNoValue);
		}
		if(geotop::common::Variables::files[fsndur] != geotop::input::gStringNoValue)
			{
			snow->t_snow.resize(snow->t_snow.size(),0.0);
			}
	}
	
	if(par->max_glac_layers>0 && par->output_glac_bin==1){
		if(geotop::common::Variables::files[fglacmelt] != geotop::input::gStringNoValue) glac->MELTED.resize(glac->MELTED.size(), 0.);
		if(geotop::common::Variables::files[fglacsubl] != geotop::input::gStringNoValue) glac->SUBL.resize(glac->SUBL.size(), 0.);
	}
	
	if(par->output_surfenergy_bin == 1){
		
		if(geotop::common::Variables::files[fradnet] != geotop::input::gStringNoValue) egy->Rn_mean.resize(egy->Rn_mean.size(), 0.0);
		if(geotop::common::Variables::files[fradLWin] != geotop::input::gStringNoValue) egy->LWin_mean.resize(egy->LWin_mean.size(), 0.0);
		if(geotop::common::Variables::files[fradLW] != geotop::input::gStringNoValue) egy->LW_mean.resize(egy->LW_mean.size(), 0.0);
		if(geotop::common::Variables::files[fradSW] != geotop::input::gStringNoValue) egy->SW_mean.resize(egy->SW_mean.size(), 0.0);
		if(geotop::common::Variables::files[fradSWin] != geotop::input::gStringNoValue) egy->Rswdown_mean.resize(egy->Rswdown_mean.size(), 0.0);
		if(geotop::common::Variables::files[fradSWinbeam] != geotop::input::gStringNoValue) egy->Rswbeam_mean.resize(egy->Rswbeam_mean.size(), 0.0);
		if(geotop::common::Variables::files[fshadow] != geotop::input::gStringNoValue){
			egy->nDt_shadow.resize(egy->nDt_shadow.size(),0.0);
			egy->nDt_sun.resize(egy->nDt_sun.size(),0.0);
		}
		
		if(geotop::common::Variables::files[fG] != geotop::input::gStringNoValue) egy->SEB_mean.resize(egy->SEB_mean.size(), 0.);
		if(geotop::common::Variables::files[fH] != geotop::input::gStringNoValue) egy->H_mean.resize(egy->H_mean.size(), 0.0);
		if(geotop::common::Variables::files[fLE] != geotop::input::gStringNoValue) egy->ET_mean.resize(egy->ET_mean.size(), 0.0);
		if(geotop::common::Variables::files[fTs] != geotop::input::gStringNoValue) egy->Ts_mean.resize(egy->Ts_mean.size(), 0.0);
	}
	
	if (par->output_meteo_bin == 1){
		if(geotop::common::Variables::files[fTa] != geotop::input::gStringNoValue) met->Tamean.resize(met->Tamean.size(), 0.0);
		if(geotop::common::Variables::files[fwspd] != geotop::input::gStringNoValue) met->Vspdmean.resize(met->Vspdmean.size(), 0.0);
		if(geotop::common::Variables::files[fwdir] != geotop::input::gStringNoValue) met->Vdirmean.resize(met->Vdirmean.size(),0.0);

		if(geotop::common::Variables::files[frh] != geotop::input::gStringNoValue) met->RHmean.resize(met->RHmean.size(), 0.0);
		if(geotop::common::Variables::files[fprec] != geotop::input::gStringNoValue){
			wat->PrTOT_mean.resize(wat->PrTOT_mean.size(),0.0);
			wat->PrSNW_mean.resize(wat->PrSNW_mean.size(),0.0);
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

