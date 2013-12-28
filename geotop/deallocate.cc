
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#include "deallocate.h"
#include "geotop_common.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,Energy *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times){
  void dealloc_all(Topo *top,Soil *sl,Land *land,Water *wat,Channel *cnet,Par *par,Energy *egy,Snow *snow, Glacier *glac, Meteo *met, Times *times){

	long i,j,r,l;
	
	printf("Close files\n");
	if(files[fpointwriteend] != string_novalue) fclose(ffpoint);
	if(files[fsnTzwriteend] != string_novalue) fclose(ffsnowT);
	if(files[fsnlzwriteend] != string_novalue) fclose(ffsnowl);
	if(files[fsnizwriteend] != string_novalue) fclose(ffsnowi);
	if(files[fsndzwriteend] != string_novalue) fclose(ffsnowd);
	if(files[fglzwriteend] != string_novalue) fclose(ffglac);
	if(files[fbaswriteend] != string_novalue) fclose(ffbas);
	if(files[fTzwriteend] != string_novalue) fclose(ffT);
	if(files[fTzavwriteend] != string_novalue) fclose(ffTav);
	if(files[fpsizwriteend] != string_novalue) fclose(ffpsi);
	if(files[fpsiztotwriteend] != string_novalue) fclose(ffpsitot);
	if(files[fliqzwriteend] != string_novalue) fclose(ffliq);
	if(files[ficezwriteend] != string_novalue) fclose(ffice);
	
	printf("Deallocating global variables\n"); 
	if(par->state_pixel==1){
		for (i=0; i<otot; i++) {
			free(odpnt[i]);
			free(odp[i]);
		}
		free(odpnt);
		free(odp);
	}
	for (i=0; i<otot; i++) {
		free(hpnt[i]);
	}
	free(hpnt);
	free(opnt);
	free(ipnt);
	
	for (i=0; i<ootot; i++) {
		free(hbsn[i]);
	}
	free(odbsn);
	free(odb);
	free(hbsn);
	free(obsn);
	free(ibsn);

/*
	for (j=0; j<10; j++) {
		free(hsnw[j]);
		free(hglc[j]);
	}
	free(hsnw);
	free(hglc);
	
	free(osnw);
	free(oglc);
		
	for (j=0; j<6; j++) {
		free(hsl[j]);
	}
	free(hsl);
	free(osl);
*/
 
//	free(WORKING_DIRECTORY);
	
	/* Deallocation of struct SOIL "sl": */
	printf("Deallocating soil\n");
//	free_doublematrix(sl->Ptot);

    //commenting the following statements, since they have no body
    /*
	if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0);// free_doublematrix(sl->T_av_tensor);
	if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0);// free_doublematrix(sl->T_av_tensor);
	if(strcmp(files[fliqav] , string_novalue) != 0);// free_doublematrix(sl->thw_av_tensor);
	if(strcmp(files[ficeav] , string_novalue) != 0);// free_doublematrix(sl->thi_av_tensor);
    */

//	free_doublematrix(sl->T_av_tensor);
//	free_doublematrix(sl->th);
//	free_longmatrix(sl->type);
//	free_doubletensor(sl->pa);
//	free_doubletensor(sl->ET);
//	deallocate_soil_state(sl->SS);
//	deallocate_veg_state(sl->VS);
	
    //commenting the following statements, since they have no body
    /*
	if(par->state_pixel == 1){
		if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) ;// free_doublematrix(sl->Tzplot);
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) ;//free_doublematrix(sl->Tzavplot);
		if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) ;// free_doublematrix(sl->Ptotzplot);
		if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) ;// free_doublematrix(sl->Pzplot);
		if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) ;// free_doublematrix(sl->thzplot);
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) ;// free_doublematrix(sl->thzavplot);
		if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) ;// free_doublematrix(sl->thizplot);
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) ;//free_doublematrix(sl->thizavplot);
	}
    */

	//free(sl);
	
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
		
	//	free_longmatrix(top->horizon_point);
	}
//	free_doublematrix(par->chkpt);
//	free_doublematrix(top->sky);
//	free_doublematrix(top->Z0);
//	free_doublematrix(top->East);
//	free_doublematrix(top->North);
//	free_shortmatrix(top->pixel_type);
//	free_doublematrix(top->aspect);
//	free_doublematrix(top->slope);
//	free_doublematrix(top->dzdE);
//	free_doublematrix(top->dzdN);
//	free_doublematrix(top->curvature1);
//	free_doublematrix(top->curvature2);
//	free_doublematrix(top->curvature3);
//	free_doublematrix(top->curvature4);
//	free_longmatrix(top->Jdown);
//	free_doublematrix(top->Qdown);
//	if(par->point_sim==0); free_shortmatrix(top->is_on_border);
	
//	free_longmatrix(top->lrc_cont);
	
	for(l=0;l<=Nl;l++){
		for(r=1;r<=Nr;r++){
			free(top->i_cont[l][r]);
		}
		free(top->i_cont[l]);
	}
	free(top->i_cont);
	
//	free_longmatrix(top->rc_cont);
	for(r=1;r<=Nr;r++){
		free(top->j_cont[r]);
	}
	free(top->j_cont);
		
//	free_doubletensor(top->Z);
	
//	free_longvector(top->Lp);
//	free_longvector(top->Li);
	
//	free_longmatrix(top->BC_counter);
//	free_doublevector(top->BC_DepthFreeSurface);
		
	if (par->point_sim==1) {
	//	free_doublematrix(top->latitude);
	//	free_doublematrix(top->longitude);
	}
	
	//free(top);
	
	/* Deallocation of struct LAND "land": */
	printf("Deallocating land\n");
//	free_doublematrix(land->LC);
//	free_doublematrix(land->delay);
//	free_shortmatrix(land->shadow);
//	free_doublematrix(land->ty);
//	free_doublematrix(land->root_fraction);
//	free_doublevector(land->vegpar);
	
	for(i=0;i<par->n_landuses;i++){
		free(land->vegparv[i]);
		
	//	if(par->vegflag->co[i+1]==1){
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
	
	//free(land);
	
	/* Deallocation of struct WATER "water": */
	printf("Deallocating water\n"); 
//	free_doublematrix(wat->PrecTot);
//	free_doublematrix(wat->Pnet);

	if (par->output_meteo_bin == 1 && files[fprec] != string_novalue){
		//free_doublevector(wat->PrTOT_mean);
		//free_doublevector(wat->PrSNW_mean);
		//free_doublevector(wat->Pt);
		//free_doublevector(wat->Ps);
	}
	
	//free_doublevector(wat->h_sup);
	//free_doublevector(wat->Lx);
	//free_doublevector(wat->H0);
	//free_doublevector(wat->H1);
	//free_doublevector(wat->dH);
	//free_doublevector(wat->B);
	//free_doublevector(wat->f);
	//free_doublevector(wat->df);
	//free_doublematrix(wat->Klat);
	//free_doublematrix(wat->Kbottom);
	free(wat);
	
	/* Deallocation of struct CHANNEL "channel": */
	printf("Deallocating channel network\n"); 
	//free_longvector(cnet->r);
	//free_longvector(cnet->c);
	//free_longmatrix(cnet->ch);
	//free_longvector(cnet->ch_down);
	//free_doublevector(cnet->Vsup);
	//free_doublevector(cnet->Vsub);
	//free_doublevector(cnet->Vsup_cum);
	//free_doublevector(cnet->Vsub_cum);
	//free_doublevector(cnet->h_sup);
	//free_doublevector(cnet->length);
	for (l=0; l<=Nl; l++) {
		free(cnet->ch3[l]);
	}
	free(cnet->ch3);
//	free_longmatrix(cnet->lch);
//	free_longvector(cnet->soil_type);
//	free_doublematrix(cnet->th);
//	free_doublematrix(cnet->ET);
//	free_doublevector(cnet->Kbottom);
	//deallocate_soil_state(cnet->SS);
//	free(cnet);
	
	/* Deallocation of struct T_INIT "UV": */
	printf("Deallocating UV\n"); 
//	free_doublevector(UV->U);
//	free_doublevector(UV->V);
//	free(UV);
	
	/* Deallocation of struct ENERGY "egy": */
	printf("Deallocating egy\n");  	
	if(par->output_surfenergy_bin == 1){
		if(files[fradnet] != string_novalue){
			//free_doublevector(egy->Rn_mean);
			//free_doublevector(egy->Rn);
		}
		if(files[fradLWin] != string_novalue){
			//free_doublevector(egy->LWin_mean);
			//free_doublevector(egy->LWin);
		}	
		if(files[fradLW] != string_novalue){
			//free_doublevector(egy->LW_mean);
			//free_doublevector(egy->LW);
		}	
		if(files[fradSW] != string_novalue){
			//free_doublevector(egy->SW_mean);
			//free_doublevector(egy->SW);
		}
		if(files[fradSWin] != string_novalue){
			//free_doublevector(egy->Rswdown_mean);
			//free_doublevector(egy->SWin);
		}
		if(files[fradSWinbeam] != string_novalue){
			//free_doublevector(egy->Rswbeam_mean);
			//free_doublevector(egy->SWinb);
		}
		if(files[fLE] != string_novalue){
			//free_doublevector(egy->ET_mean);
			//free_doublevector(egy->LE);
		}
		if(files[fG] != string_novalue){
			//free_doublevector(egy->SEB_mean);
			//free_doublevector(egy->G);
		}
		if(files[fH] != string_novalue){
			//free_doublevector(egy->H_mean);
			//free_doublevector(egy->H);
		}
		if(files[fTs] != string_novalue){
			//free_doublevector(egy->Ts_mean);
			//free_doublevector(egy->Ts);
		}
		if(files[fshadow] != string_novalue){
			//free_longvector(egy->nDt_shadow);
			//free_longvector(egy->nDt_sun);
			//free_shortvector(egy->shad);
		}
	}
	
	free(egy->sun);
	
//	if(times->JD_plots->nh > 1){
	if(times->JD_plots.size() > 1){
		if(files[pH] != string_novalue || files[pHg] != string_novalue || files[pG] != string_novalue){
			//free_doublevector(egy->Hgplot);
			//free_doublevector(egy->Hgp);
		}
		if(files[pLE] != string_novalue || files[pLEg] != string_novalue || files[pG] != string_novalue){
			//free_doublevector(egy->LEgplot);
			//free_doublevector(egy->LEgp);
		}
		if(files[pH] != string_novalue || files[pHv] != string_novalue){
			//free_doublevector(egy->Hvplot);
			//free_doublevector(egy->Hvp);
		}
		if(files[pLE] != string_novalue || files[pLEv] != string_novalue){
			//free_doublevector(egy->LEvplot);
			//free_doublevector(egy->LEvp);
		}
		if(files[pSWin] != string_novalue){
			//free_doublevector(egy->SWinplot);
			//free_doublevector(egy->SWinp);
		}
		if(files[pSWg] != string_novalue || files[pG] != string_novalue){
			//free_doublevector(egy->SWgplot);
			//free_doublevector(egy->SWgp);
		}
		if(files[pSWv] != string_novalue){
			//free_doublevector(egy->SWvplot);
			//free_doublevector(egy->SWvp);
		}
		if(files[pLWin] != string_novalue){
			//free_doublevector(egy->LWinplot);
			//free_doublevector(egy->LWinp);
		}
		if(files[pLWg] != string_novalue || files[pG] != string_novalue){
			//free_doublevector(egy->LWgplot);
			//free_doublevector(egy->LWgp);
		}
		if(files[pLWv] != string_novalue){
			//free_doublevector(egy->LWvplot);
			//free_doublevector(egy->LWvp);
		}
		if(files[pTs] != string_novalue){
			//free_doublevector(egy->Tsplot);
			//free_doublevector(egy->Tsp);
		}
		if(files[pTg] != string_novalue){
			//free_doublevector(egy->Tgplot);
			//free_doublevector(egy->Tgp);
		}
		if(files[pTv] != string_novalue){
			//free_doublevector(egy->Tvplot);
		}
	}
	
	//free_doublevector(egy->Dlayer);
	//free_doublevector(egy->liq);
	//free_doublevector(egy->ice);
	//free_doublevector(egy->Temp);
	//free_doublevector(egy->deltaw);
	//free_doublevector(egy->SWlayer);
	//free_doublevector(egy->soil_transp_layer);
	//free_doublevector(egy->dFenergy);
	//free_doublevector(egy->udFenergy);
	//free_doublevector(egy->Kth0);
	//free_doublevector(egy->Kth1);
	//free_doublevector(egy->Fenergy);
	//free_doublevector(egy->Newton_dir);
	//free_doublevector(egy->T0);
	//free_doublevector(egy->T1);
	//free_doublevector(egy->Tstar);
	//free_doublevector(egy->THETA);
	//free_doublevector(egy->soil_evap_layer_bare);
	//free_doublevector(egy->soil_evap_layer_veg);
	//free_doublematrix(egy->Tgskin_surr);
	//free_doublematrix(egy->SWrefl_surr);
	
	//free(egy);
	
	
	/* Deallocation of struct SNOW "snow": */
	printf("Deallocating snow\n"); 

//	if(times->JD_plots->nh > 1){
	if(times->JD_plots.size() > 1){
		if(files[pD] != string_novalue){
		//	free_doublevector(snow->Dplot);
		}
	}
			
//	deallocate_statevar_3D(snow->S);
	
	if(par->blowing_snow==1){
		deallocate_statevar_1D(snow->S_for_BS);
	//	free_doublematrix(snow->Nabla2_Qtrans);
	//	free_doublematrix(snow->Qsub);
	//	free_doublematrix(snow->Qsub_x);
	//	free_doublematrix(snow->Qsub_y);
	//	free_doublematrix(snow->Qtrans);
	//	free_doublematrix(snow->Qtrans_x);
	//	free_doublematrix(snow->Qtrans_y);
	//	free_doublematrix(snow->Qsalt);
		
		if(par->output_snow_bin == 1){
		//	free_doublematrix(snow->Wtrans_plot);
		//	free_doublematrix(snow->Wsubl_plot);
		}
	}
	
	//free_doublevector(snow->age);
	if(par->output_snow_bin == 1){
		if(files[fsndur] != string_novalue){
			//free_doublevector(snow->t_snow);
			//free_shortvector(snow->yes);
		}
		if(files[fsnowmelt] != string_novalue){
			//free_doublevector(snow->MELTED);
			//free_doublevector(snow->melted);
		}
		if(files[fsnowsubl] != string_novalue){
			//free_doublevector(snow->SUBL);
			//free_doublevector(snow->subl);
		}
	}
	
	// if(par->blowing_snow==1) free_longvector(snow->change_dir_wind);
	if(par->blowing_snow==1) // free_longvector(snow->change_dir_wind);
	free(snow);
	
	printf("Deallocating glacier\n");
	if(par->max_glac_layers>0){ deallocate_statevar_3D(glac->G);
		if(par->output_glac_bin == 1){
			if(files[fglacmelt] != string_novalue){
				//free_doublevector(glac->MELTED);
				//free_doublevector(glac->melted);
			}
			if(files[fglacsubl] != string_novalue){
				//free_doublevector(glac->SUBL);
				//free_doublevector(glac->subl);
			}
		}
	}
	//free(glac);
	
	printf("Deallocating met\n"); 
	
//	free_doublematrix(met->Tgrid);
//	free_doublematrix(met->Pgrid);
//	free_doublematrix(met->Vgrid);
//	free_doublematrix(met->Vdir);
//	free_doublematrix(met->RHgrid);
//	free_doublematrix(met->tau_cl_map);
//	free_doublematrix(met->tau_cl_av_map);
//	free_shortmatrix(met->tau_cl_map_yes);
//	free_shortmatrix(met->tau_cl_av_map_yes);
	if (par->output_meteo_bin == 1){
		//if(strcmp(files[fTa] , string_novalue) != 0) free_doublevector(met->Tamean);
		//if(strcmp(files[fwspd] , string_novalue) != 0) free_doublevector(met->Vspdmean);
		//if(strcmp(files[fwdir] , string_novalue) != 0) free_doublevector(met->Vdirmean);
		//if(strcmp(files[frh] , string_novalue) != 0) free_doublevector(met->RHmean);
	}		
//	if(times->JD_plots->nh > 1){
	if(times->JD_plots.size() > 1){
		if(files[pTa] != string_novalue) {
			//free_doublevector(met->Taplot);
		}
		if(files[pVspd] != string_novalue || files[pVdir] != string_novalue){
			//free_doublevector(met->Vxplot);
			//free_doublevector(met->Vyplot);
		}
		if(files[pRH] != string_novalue)
		{
			//free_doublevector(met->RHplot);
		}
	}		
	
//	for(i=0;i<met->st->Z->nh;i++){
	for(i=0;i<long(met->st->Z.size()-1);i++){
//		for(j=0;j<met->numlines[i];j++){
//			free(met->data[i][j]);
//		}
//		free(met->data[i]);
		
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
	
	//free_longvector(met->imeteo_stations);
	dealloc_meteostations(met->st); 
	
	free(met->qinv);
	if (par->qin == 1) {
		for (i=0; i<met->qinsnr; i++) {
			free(met->qins[i]);
		}
	}
	
	//free(met);
	
	printf("Deallocating times\n");
	//free_doublevector(times->JD_plots);
	free(times->Dt_vector);
	if(par->tsteps_from_file==1){
		for(j=0;j<times->numlinesDt_matrix;j++){
			free(times->Dt_matrix[j]);
		}
		free(times->Dt_matrix);
	}
	free(times);
	
	/* Deallocation of struct PAR "par": */
	printf("Deallocating par\n"); 
	//free_shortvector(par->vegflag);
	if(par->state_pixel == 1){
		//free_longvector(par->jplot);
		//free_longmatrix(par->rc);
		//free_longvector(par->IDpoint);
	}

	//free_doublevector(par->saving_points);
	
	//free_doublevector(par->init_date);
	//free_doublevector(par->end_date);
	//free_longvector(par->run_times);
	
	//if (par->point_sim == 1) //free_doublematrix(par->maxSWE);
	
	//free_shortvector(par->plot_discharge_with_Dt_integration);
	//free_shortvector(par->plot_point_with_Dt_integration);
	//free_shortvector(par->plot_basin_with_Dt_integration);
	
    //free_doublevector(par->Dtplot_point);
	//free_doublevector(par->Dtplot_basin);
	//free_doublevector(par->Dtplot_discharge);
	
	//free_doublevector(par->output_soil);
	//free_doublevector(par->output_snow);
	//free_doublevector(par->output_glac);
	//free_doublevector(par->output_surfenergy);
	//free_doublevector(par->output_vegetation);
	//free_doublevector(par->output_meteo);
	
	//free_doublevector(par->soil_plot_depths);
	//free_doublevector(par->snow_plot_depths);
	//free_doublevector(par->glac_plot_depths);
	
	//free_shortvector(par->linear_interpolation_meteo);
	
	//free_longvector(par->inf_snow_layers);
	//free_longvector(par->inf_glac_layers);
	
	//free(par);
	
	/* Deallocation of struct FILENAMES "filenames": */
	
	printf("Deallocating novalues\n"); 
	free(string_novalue);
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void dealloc_meteostations(METEO_STATIONS *st){
  void dealloc_meteostations(MeteoStations *st){
	
	//free_doublevector(st->E);
	//free_doublevector(st->N);
	//free_doublevector(st->lat);
	//free_doublevector(st->lon);
	//free_doublevector(st->Z);
	//free_doublevector(st->sky);
	//free_doublevector(st->ST);
	//free_doublevector(st->Vheight);
	//free_doublevector(st->Theight);
	//free_doublevector(st->tau_cloud_meteoST);
	//free_doublevector(st->tau_cloud_av_meteoST);
	//free_shortvector(st->flag_SW_meteoST);
	//free_shortvector(st->tau_cloud_av_yes_meteoST);
	//free_shortvector(st->tau_cloud_yes_meteoST);
	//free(st);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void deallocate_soil_state(SOIL_STATE *S){
void deallocate_soil_state(SoilState *S){
	
//	free_doublematrix(S->T);
//	free_doublematrix(S->P);
//	free_doublematrix(S->thi);
	free(S);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void deallocate_veg_state(STATE_VEG *V){
void deallocate_veg_state(StateVeg *V){

	//free_doublevector(V->Tv);
	//free_doublevector(V->wsnow);
	//free_doublevector(V->wrain);
	free(V);
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//void reset_to_zero(PAR *par, SOIL *sl, LAND *land, SNOW *snow, GLACIER *glac, Energy *egy, METEO *met, WATER *wat){
  void reset_to_zero(Par *par, Soil *sl, Land *land, Snow *snow, Glacier *glac, Energy *egy, Meteo *met, Water *wat){

	long i, j;
	
	if(par->state_pixel == 1){
		if(files[fTzav] != string_novalue || files[fTzavwriteend] != string_novalue) {
			//initialize_doublematrix(sl->Tzavplot,0.);
			sl->Tzavplot.resize(sl->Tzavplot.getRows(), sl->Tzavplot.getCols(),0.0);
		}

		if(files[fliqzav] != string_novalue || files[fliqzavwriteend] != string_novalue) {
		//	initialize_doublematrix(sl->thzavplot,0.);
			sl->thzavplot.resize(sl->thzavplot.getRows(),sl->thzavplot.getCols(),0.0);
		}
		if(files[ficezav] != string_novalue || files[ficezavwriteend] != string_novalue) {
		//	initialize_doublematrix(sl->thizavplot,0.);
			sl->thizavplot.resize(sl->thizavplot.getRows(),sl->thizavplot.getCols(),0.0);
		}
		
	//	for(i=1;i<=par->rc->nrh;i++){
		for(i=1;i<par->rc.getRows();i++){
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
		if(files[fTav] != string_novalue || files[fTavsup] != string_novalue) sl->T_av_tensor.resize(sl->T_av_tensor.getRows(),sl->T_av_tensor.getCols(), 0.);
		if(files[ficeav] != string_novalue) sl->thi_av_tensor.resize(sl->thi_av_tensor.getRows(),sl->thi_av_tensor.getCols(), 0.0);
		if(files[fliqav] != string_novalue) sl->thw_av_tensor.resize(sl->thw_av_tensor.getRows(),sl->thw_av_tensor.getCols(), 0.0);
	}
	
	if(par->output_snow_bin == 1){
		if(files[fsnowmelt] != string_novalue) snow->MELTED.resize(snow->MELTED.size(), 0.);
		if(files[fsnowsubl] != string_novalue) snow->SUBL.resize(snow->SUBL.size(), 0.0);
		if(files[fswe] != string_novalue && par->blowing_snow==1){
			initmatrix(0.0, snow->Wtrans_plot, land->LC, number_novalue);
			initmatrix(0.0, snow->Wsubl_plot, land->LC, number_novalue);
		}
		if(files[fsndur] != string_novalue)
			{
		//	initialize_doublevector(snow->t_snow, 0.);
			snow->t_snow.resize(snow->t_snow.size(),0.0);
			}
	}
	
	if(par->max_glac_layers>0 && par->output_glac_bin==1){
		if(files[fglacmelt] != string_novalue) glac->MELTED.resize(glac->MELTED.size(), 0.);
		if(files[fglacsubl] != string_novalue) glac->SUBL.resize(glac->SUBL.size(), 0.);
	}
	
	if(par->output_surfenergy_bin == 1){
		
		if(files[fradnet] != string_novalue) egy->Rn_mean.resize(egy->Rn_mean.size(), 0.0);
		if(files[fradLWin] != string_novalue) egy->LWin_mean.resize(egy->LWin_mean.size(), 0.0);
		if(files[fradLW] != string_novalue) egy->LW_mean.resize(egy->LW_mean.size(), 0.0);
		if(files[fradSW] != string_novalue) egy->SW_mean.resize(egy->SW_mean.size(), 0.0);
		if(files[fradSWin] != string_novalue) egy->Rswdown_mean.resize(egy->Rswdown_mean.size(), 0.0);
		if(files[fradSWinbeam] != string_novalue) egy->Rswbeam_mean.resize(egy->Rswbeam_mean.size(), 0.0);
		if(files[fshadow] != string_novalue){
		//	initialize_longvector(egy->nDt_shadow, 0);
			egy->nDt_shadow.resize(egy->nDt_shadow.size(),0.0);
		//	initialize_longvector(egy->nDt_sun, 0);
			egy->nDt_sun.resize(egy->nDt_sun.size(),0.0);
		}
		
		if(files[fG] != string_novalue) egy->SEB_mean.resize(egy->SEB_mean.size(), 0.);
		if(files[fH] != string_novalue) egy->H_mean.resize(egy->H_mean.size(), 0.0);
		if(files[fLE] != string_novalue) egy->ET_mean.resize(egy->ET_mean.size(), 0.0);
		if(files[fTs] != string_novalue) egy->Ts_mean.resize(egy->Ts_mean.size(), 0.0);
	}
	
	if (par->output_meteo_bin == 1){
		if(files[fTa] != string_novalue) met->Tamean.resize(met->Tamean.size(), 0.0);
		if(files[fwspd] != string_novalue) met->Vspdmean.resize(met->Vspdmean.size(), 0.0);
		if(files[fwdir] != string_novalue) met->Vdirmean.resize(met->Vdirmean.size(),0.0);

		if(files[frh] != string_novalue) met->RHmean.resize(met->RHmean.size(), 0.0);
		if(files[fprec] != string_novalue){
		//	initialize_doublevector(wat->PrTOT_mean, 0.);
			wat->PrTOT_mean.resize(wat->PrTOT_mean.size(),0.0);
		//	initialize_doublevector(wat->PrSNW_mean, 0.);
			wat->PrSNW_mean.resize(wat->PrSNW_mean.size(),0.0);
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

