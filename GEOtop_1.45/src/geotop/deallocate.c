
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
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times){
	
	long i,j,r,l;
	
	printf("Close files\n");
	if(strcmp(files[fpointwriteend] , string_novalue) != 0) fclose(ffpoint);
	if(strcmp(files[fsnzwriteend] , string_novalue) != 0) fclose(ffsnow);
	if(strcmp(files[fglzwriteend] , string_novalue) != 0) fclose(ffglac);
	if(strcmp(files[fbaswriteend] , string_novalue) != 0) fclose(ffbas);
	if(strcmp(files[fTzwriteend] , string_novalue) != 0) fclose(ffT);
	if(strcmp(files[fTzavwriteend] , string_novalue) != 0) fclose(ffTav);
	if(strcmp(files[fpsizwriteend] , string_novalue) != 0) fclose(ffpsi);
	if(strcmp(files[fpsiztotwriteend] , string_novalue) != 0) fclose(ffpsitot);
	if(strcmp(files[fliqzwriteend] , string_novalue) != 0) fclose(ffliq);
	if(strcmp(files[ficezwriteend] , string_novalue) != 0) fclose(ffice);
	
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
	
	free(WORKING_DIRECTORY);
	
	/* Deallocation of struct SOIL "sl": */
	printf("Deallocating soil\n");
	free_doublematrix(sl->Ptot);
	if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0) free_doublematrix(sl->T_av_tensor);
	if(strcmp(files[fliqav] , string_novalue) != 0) free_doublematrix(sl->thw_av_tensor);
	if(strcmp(files[ficeav] , string_novalue) != 0) free_doublematrix(sl->thi_av_tensor);
	//TODO: Hack
	if(strcmp(files[fpnet] , string_novalue) != 0) free_doublevector(sl->Pnetcum);
	if(strcmp(files[fevap] , string_novalue) != 0) free_doublevector(sl->ETcum);
	//free_doublematrix(sl->T_av_tensor);
	free_doublematrix(sl->th); 
	free_longmatrix(sl->type);
	free_doubletensor(sl->pa);
	free_doubletensor(sl->ET);	
	deallocate_soil_state(sl->SS);
	deallocate_veg_state(sl->VS);
	
	if(par->state_pixel == 1){
		if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) free_doublematrix(sl->Tzplot);
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) free_doublematrix(sl->Tzavplot);
		if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) free_doublematrix(sl->Ptotzplot);
		if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) free_doublematrix(sl->Pzplot);
		if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) free_doublematrix(sl->thzplot);
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) free_doublematrix(sl->thzavplot);
		if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) free_doublematrix(sl->thizplot);
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) free_doublematrix(sl->thizavplot);
	}
	free(sl);
	
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
		
		free_longmatrix(top->horizon_point);
	}
	free_doublematrix(par->chkpt);
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
	if(par->point_sim==0) free_shortmatrix(top->is_on_border);
	
	free_longmatrix(top->lrc_cont);
	
	for(l=0;l<=Nl;l++){
		for(r=1;r<=Nr;r++){
			free(top->i_cont[l][r]);
		}
		free(top->i_cont[l]);
	}
	free(top->i_cont);
	
	free_longmatrix(top->rc_cont);
	for(r=1;r<=Nr;r++){
		free(top->j_cont[r]);
	}
	free(top->j_cont);
		
	free_doubletensor(top->Z);
	
	free_longvector(top->Lp);
	free_longvector(top->Li);
	
	free_longmatrix(top->BC_counter);
	free_doublevector(top->BC_DepthFreeSurface);
		
	if (par->point_sim==1) {
		free_doublematrix(top->latitude);
		free_doublematrix(top->longitude);
	}
	
	free(top);
	
	
	/* Deallocation of struct LAND "land": */
	printf("Deallocating land\n");
	free_doublematrix(land->LC);
	free_doublematrix(land->delay);
	free_shortmatrix(land->shadow);
	free_doublematrix(land->ty);
	free_doublematrix(land->root_fraction);
	free_doublevector(land->vegpar);
	
	for(i=0;i<par->n_landuses;i++){
		free(land->vegparv[i]);
		
		if(par->vegflag->co[i+1]==1){ 
			for(j=0;j<land->NumlinesVegTimeDepData[i];j++){
				free(land->vegpars[i][j]);
			}
			free(land->vegpars[i]);
		}
	}
	
	free(land->vegpars);
	free(land->vegparv);
	free(land->NumlinesVegTimeDepData);
	
	free(land);
	
	/* Deallocation of struct WATER "water": */
	printf("Deallocating water\n"); 
	free_doublematrix(wat->PrecTot);
	free_doublematrix(wat->Pnet);

	if (par->output_meteo_bin == 1 && strcmp(files[fprec] , string_novalue) != 0){
		free_doublevector(wat->PrTOT_mean);
		free_doublevector(wat->PrSNW_mean);
		free_doublevector(wat->Pt);
		free_doublevector(wat->Ps);
	}
	
	free_doublevector(wat->h_sup);
	free_doublevector(wat->Lx);
	free_doublevector(wat->H0);
	free_doublevector(wat->H1);
	free_doublevector(wat->dH);
	free_doublevector(wat->B);
	free_doublevector(wat->f);
	free_doublevector(wat->df);
	free_doublematrix(wat->Klat);
	free_doublematrix(wat->Kbottom);
	free(wat);
	
	/* Deallocation of struct CHANNEL "channel": */
	printf("Deallocating channel network\n"); 
	free_longvector(cnet->r);
	free_longvector(cnet->c);
	free_longmatrix(cnet->ch);
	free_longvector(cnet->ch_down);
	free_doublevector(cnet->Vsup);
	free_doublevector(cnet->Vsub);
	//free_doublevector(cnet->Vsup_cum);
	//free_doublevector(cnet->Vsub_cum);
	free_doublevector(cnet->h_sup);
	free_doublevector(cnet->length);
	for (l=0; l<=Nl; l++) {
		free(cnet->ch3[l]);
	}
	free(cnet->ch3);
	free_longmatrix(cnet->lch);
	free_longvector(cnet->soil_type);
	free_doublematrix(cnet->th);
	free_doublematrix(cnet->ET);
	free_doublevector(cnet->Kbottom);
	deallocate_soil_state(cnet->SS);
	free(cnet);
	
	/* Deallocation of struct T_INIT "UV": */
	printf("Deallocating UV\n"); 
	free_doublevector(UV->U);
	free_doublevector(UV->V);
	free(UV);
	
	/* Deallocation of struct ENERGY "egy": */
	printf("Deallocating egy\n");  	
	if(par->output_surfenergy_bin == 1){
		if(strcmp(files[fradnet] , string_novalue) != 0){
			free_doublevector(egy->Rn_mean);
			free_doublevector(egy->Rn);
		}
		if(strcmp(files[fradLWin] , string_novalue) != 0){
			free_doublevector(egy->LWin_mean);
			free_doublevector(egy->LWin);
		}	
		if(strcmp(files[fradLW] , string_novalue) != 0){
			free_doublevector(egy->LW_mean);
			free_doublevector(egy->LW);
		}	
		if(strcmp(files[fradSW] , string_novalue) != 0){
			free_doublevector(egy->SW_mean);
			free_doublevector(egy->SW);
		}
		if(strcmp(files[fradSWin] , string_novalue) != 0){
			free_doublevector(egy->Rswdown_mean);
			free_doublevector(egy->SWin);
		}
		if(strcmp(files[fradSWinbeam] , string_novalue) != 0){
			free_doublevector(egy->Rswbeam_mean);
			free_doublevector(egy->SWinb);
		}
		if(strcmp(files[fLE] , string_novalue) != 0){
			free_doublevector(egy->ET_mean);
			free_doublevector(egy->LE);
		}
		if(strcmp(files[fG] , string_novalue) != 0){
			free_doublevector(egy->SEB_mean);
			free_doublevector(egy->G);
		}
		if(strcmp(files[fH] , string_novalue) != 0){
			free_doublevector(egy->H_mean);
			free_doublevector(egy->H);
		}
		if(strcmp(files[fTs] , string_novalue) != 0){
			free_doublevector(egy->Ts_mean);
			free_doublevector(egy->Ts);
		}
		if(strcmp(files[fshadow] , string_novalue) != 0){
			free_longvector(egy->nDt_shadow);
			free_longvector(egy->nDt_sun); 
			free_shortvector(egy->shad);
		}
	}
	
	free(egy->sun);
	
	if(times->JD_plots->nh > 1){
		if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			free_doublevector(egy->Hgplot);
			free_doublevector(egy->Hgp);
		}
		if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			free_doublevector(egy->LEgplot);
			free_doublevector(egy->LEgp);
		}
		if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0){
			free_doublevector(egy->Hvplot);
			free_doublevector(egy->Hvp);
		}
		if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0){
			free_doublevector(egy->LEvplot);
			free_doublevector(egy->LEvp);
		}
		if(strcmp(files[pSWin] , string_novalue) != 0){
			free_doublevector(egy->SWinplot);
			free_doublevector(egy->SWinp);
		}
		if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			free_doublevector(egy->SWgplot);
			free_doublevector(egy->SWgp);
		}
		if(strcmp(files[pSWv] , string_novalue) != 0){
			free_doublevector(egy->SWvplot);
			free_doublevector(egy->SWvp);
		}
		if(strcmp(files[pLWin] , string_novalue) != 0){
			free_doublevector(egy->LWinplot);
			free_doublevector(egy->LWinp);
		}
		if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0){
			free_doublevector(egy->LWgplot);
			free_doublevector(egy->LWgp);
		}
		if(strcmp(files[pLWv] , string_novalue) != 0){
			free_doublevector(egy->LWvplot);
			free_doublevector(egy->LWvp);
		}
		if(strcmp(files[pTs] , string_novalue) != 0){
			free_doublevector(egy->Tsplot);
			free_doublevector(egy->Tsp);
		}
		if(strcmp(files[pTg] , string_novalue) != 0){
			free_doublevector(egy->Tgplot);
			free_doublevector(egy->Tgp);
		}
		if(strcmp(files[pTv] , string_novalue) != 0) free_doublevector(egy->Tvplot);
	}
	
	free_doublevector(egy->Dlayer);
	free_doublevector(egy->liq);
	free_doublevector(egy->ice);
	free_doublevector(egy->Temp); 
	free_doublevector(egy->deltaw);
	free_doublevector(egy->SWlayer);
	free_doublevector(egy->soil_transp_layer);
	free_doublevector(egy->dFenergy);
	free_doublevector(egy->udFenergy);
	free_doublevector(egy->Kth0);
	free_doublevector(egy->Kth1);
	free_doublevector(egy->Fenergy);
	free_doublevector(egy->Newton_dir);
	free_doublevector(egy->T0);
	free_doublevector(egy->T1);
	free_doublevector(egy->Tstar);
	free_doublevector(egy->THETA);
	free_doublevector(egy->soil_evap_layer_bare);
	free_doublevector(egy->soil_evap_layer_veg);
	free_doublematrix(egy->Tgskin_surr);
	free_doublematrix(egy->SWrefl_surr);	
	
	free(egy);
	
	
	/* Deallocation of struct SNOW "snow": */
	printf("Deallocating snow\n"); 

	if(times->JD_plots->nh > 1){
		if(strcmp(files[pD] , string_novalue) != 0) free_doublevector(snow->Dplot);
	}
			
	deallocate_statevar_3D(snow->S);
	
	if(par->blowing_snow==1){
		deallocate_statevar_1D(snow->S_for_BS);
		free_doublematrix(snow->Nabla2_Qtrans);
		free_doublematrix(snow->Qsub);
		free_doublematrix(snow->Qsub_x);
		free_doublematrix(snow->Qsub_y);
		free_doublematrix(snow->Qtrans); 
		free_doublematrix(snow->Qtrans_x); 
		free_doublematrix(snow->Qtrans_y); 
		free_doublematrix(snow->Qsalt);
		
		if(par->output_snow_bin == 1){
			free_doublematrix(snow->Wtrans_plot);
			free_doublematrix(snow->Wsubl_plot);
		}
	}
	
	free_doublevector(snow->age); 
	if(par->output_snow_bin == 1){
		if(strcmp(files[fsndur] , string_novalue) != 0){
			free_doublevector(snow->t_snow);
			free_shortvector(snow->yes);
		}
		if(strcmp(files[fsnowmelt] , string_novalue) != 0){
			free_doublevector(snow->MELTED);
			free_doublevector(snow->melted);
		}
		if(strcmp(files[fsnowsubl] , string_novalue) != 0){
			free_doublevector(snow->SUBL);
			free_doublevector(snow->subl);
		}
	}
	
	if(par->blowing_snow==1) free_longvector(snow->change_dir_wind);
	free(snow);
	
	printf("Deallocating glacier\n");
	if(par->max_glac_layers>0){
		deallocate_statevar_3D(glac->G);
		if(par->output_glac_bin == 1){
			if(strcmp(files[fglacmelt] , string_novalue) != 0){
				free_doublevector(glac->MELTED);
				free_doublevector(glac->melted);
			}
			if(strcmp(files[fglacsubl] , string_novalue) != 0){
				free_doublevector(glac->SUBL);
				free_doublevector(glac->subl);
			}
		}
	}
	free(glac);
	
	printf("Deallocating met\n"); 
	
	free_doublematrix(met->Tgrid);
	free_doublematrix(met->Pgrid); 
	free_doublematrix(met->Vgrid);
	free_doublematrix(met->Vdir);
	free_doublematrix(met->RHgrid);
	free_doublematrix(met->tau_cl_map);
	free_doublematrix(met->tau_cl_av_map);
	free_shortmatrix(met->tau_cl_map_yes);
	free_shortmatrix(met->tau_cl_av_map_yes);
	if (par->output_meteo_bin == 1){
		if(strcmp(files[fTa] , string_novalue) != 0) free_doublevector(met->Tamean);
		if(strcmp(files[fwspd] , string_novalue) != 0) free_doublevector(met->Vspdmean);
		if(strcmp(files[fwdir] , string_novalue) != 0) free_doublevector(met->Vdirmean);
		if(strcmp(files[frh] , string_novalue) != 0) free_doublevector(met->RHmean);
	}		
	if(times->JD_plots->nh > 1){
		if(strcmp(files[pTa] , string_novalue) != 0) free_doublevector(met->Taplot);
		if(strcmp(files[pVspd] , string_novalue) != 0 || strcmp(files[pVdir] , string_novalue) != 0){
			free_doublevector(met->Vxplot);
			free_doublevector(met->Vyplot);
		}
		if(strcmp(files[pRH] , string_novalue) != 0) free_doublevector(met->RHplot);
	}		
	
	for(i=0;i<met->st->Z->nh;i++){
		
		for(j=0;j<met->numlines[i];j++){
			free(met->data[i][j]);
		}
		free(met->data[i]);
		
		free(met->var[i]);
		
		for(j=0;j<met->horizonlines[i];j++){
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
	
	free_longvector(met->imeteo_stations);
	dealloc_meteostations(met->st); 
	
	free(met->qinv);
	if (par->qin == 1) {
		for (i=0; i<met->qinsnr; i++) {
			free(met->qins[i]);
		}
	}
	
	free(met);
	
	printf("Deallocating times\n");
	free_doublevector(times->JD_plots);
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
	free_shortvector(par->vegflag);
	if(par->state_pixel == 1){
		free_longvector(par->jplot);
		free_longmatrix(par->rc);
		free_longvector(par->IDpoint);
	}

	free_doublevector(par->saving_points);
	
	free_doublevector(par->init_date);
	free_doublevector(par->end_date);
	free_longvector(par->run_times);
	
	if (par->point_sim == 1) free_doublematrix(par->maxSWE);
	
	free_shortvector(par->plot_discharge_with_Dt_integration);
	free_shortvector(par->plot_point_with_Dt_integration);
	free_shortvector(par->plot_basin_with_Dt_integration);
	
    free_doublevector(par->Dtplot_point);  
	free_doublevector(par->Dtplot_basin);
	free_doublevector(par->Dtplot_discharge);
	
	free_doublevector(par->output_soil);
	free_doublevector(par->output_snow);
	free_doublevector(par->output_glac);
	free_doublevector(par->output_surfenergy);
	free_doublevector(par->output_vegetation);
	free_doublevector(par->output_meteo);
	
	free_doublevector(par->soil_plot_depths);
	free_doublevector(par->snow_plot_depths);
	free_doublevector(par->glac_plot_depths);	
	
	free_shortvector(par->linear_interpolation_meteo);
	
	free_longvector(par->inf_snow_layers);
	free_longvector(par->inf_glac_layers);
	
	free(par);
	
	/* Deallocation of struct FILENAMES "filenames": */
	printf("Deallocating files\n"); 
	for (i=0; i<nfiles; i++) {
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

void dealloc_meteostations(METEO_STATIONS *st){
	
	free_doublevector(st->E);
	free_doublevector(st->N);
	free_doublevector(st->lat);
	free_doublevector(st->lon);
	free_doublevector(st->Z);
	free_doublevector(st->sky);
	free_doublevector(st->ST);
	free_doublevector(st->Vheight);
	free_doublevector(st->Theight);
	free_doublevector(st->tau_cloud_meteoST);
	free_doublevector(st->tau_cloud_av_meteoST);
	free_shortvector(st->flag_SW_meteoST);
	free_shortvector(st->tau_cloud_av_yes_meteoST);
	free_shortvector(st->tau_cloud_yes_meteoST);
	free(st);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_soil_state(SOIL_STATE *S){
	
	free_doublematrix(S->T);
	free_doublematrix(S->P);
	free_doublematrix(S->thi);
	free(S);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_veg_state(STATE_VEG *V){
	
	free_doublevector(V->Tv);
	free_doublevector(V->wsnow);
	free_doublevector(V->wrain);
	free(V);
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void reset_to_zero(PAR *par, SOIL *sl, LAND *land, SNOW *snow, GLACIER *glac, ENERGY *egy, METEO *met, WATER *wat){
	
	long i, j;
	
	if(par->state_pixel == 1){
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) initialize_doublematrix(sl->Tzavplot,0.); 
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) initialize_doublematrix(sl->thzavplot,0.); 
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) initialize_doublematrix(sl->thizavplot,0.); 
		
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
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0) initialize_doublematrix(sl->T_av_tensor, 0.);
		if(strcmp(files[ficeav] , string_novalue) != 0) initialize_doublematrix(sl->thi_av_tensor, 0.);
		if(strcmp(files[fliqav] , string_novalue) != 0) initialize_doublematrix(sl->thw_av_tensor, 0.);
	}
	
	if(par->output_snow_bin == 1){
		if(strcmp(files[fsnowmelt] , string_novalue) != 0) initialize_doublevector(snow->MELTED, 0.);	
		if(strcmp(files[fsnowsubl] , string_novalue) != 0) initialize_doublevector(snow->SUBL, 0.);
		if(strcmp(files[fswe] , string_novalue) != 0 && par->blowing_snow==1){
			initmatrix(0.0, snow->Wtrans_plot, land->LC, number_novalue);
			initmatrix(0.0, snow->Wsubl_plot, land->LC, number_novalue);
		}
		if(strcmp(files[fsndur] , string_novalue) != 0) initialize_doublevector(snow->t_snow, 0.);
	}
	
	if(par->max_glac_layers>0 && par->output_glac_bin==1){
		if(strcmp(files[fglacmelt] , string_novalue) != 0) initialize_doublevector(glac->MELTED, 0.);
		if(strcmp(files[fglacsubl] , string_novalue) != 0) initialize_doublevector(glac->SUBL, 0.);
	}
	
	if(par->output_surfenergy_bin == 1){
		
		if(strcmp(files[fradnet] , string_novalue) != 0) initialize_doublevector(egy->Rn_mean, 0.);
		if(strcmp(files[fradLWin] , string_novalue) != 0) initialize_doublevector(egy->LWin_mean, 0.);
		if(strcmp(files[fradLW] , string_novalue) != 0) initialize_doublevector(egy->LW_mean, 0.);
		if(strcmp(files[fradSW] , string_novalue) != 0) initialize_doublevector(egy->SW_mean, 0.);
		if(strcmp(files[fradSWin] , string_novalue) != 0) initialize_doublevector(egy->Rswdown_mean, 0.);
		if(strcmp(files[fradSWinbeam] , string_novalue) != 0) initialize_doublevector(egy->Rswbeam_mean, 0.);			
		if(strcmp(files[fshadow] , string_novalue) != 0){
			initialize_longvector(egy->nDt_shadow, 0);
			initialize_longvector(egy->nDt_sun, 0);
		}
		
		if(strcmp(files[fG] , string_novalue) != 0) initialize_doublevector(egy->SEB_mean, 0.);
		if(strcmp(files[fH] , string_novalue) != 0) initialize_doublevector(egy->H_mean, 0.);
		if(strcmp(files[fLE] , string_novalue) != 0) initialize_doublevector(egy->ET_mean, 0.);	
		if(strcmp(files[fTs] , string_novalue) != 0) initialize_doublevector(egy->Ts_mean, 0.);
	}
	
	if (par->output_meteo_bin == 1){
		if(strcmp(files[fTa] , string_novalue) != 0) initialize_doublevector(met->Tamean, 0.);
		if(strcmp(files[fwspd] , string_novalue) != 0) initialize_doublevector(met->Vspdmean, 0.);
		if(strcmp(files[fwdir] , string_novalue) != 0) initialize_doublevector(met->Vdirmean, 0.);
		if(strcmp(files[frh] , string_novalue) != 0) initialize_doublevector(met->RHmean, 0.);
		if(strcmp(files[fprec] , string_novalue) != 0){
			initialize_doublevector(wat->PrTOT_mean, 0.);
			initialize_doublevector(wat->PrSNW_mean, 0.);
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

