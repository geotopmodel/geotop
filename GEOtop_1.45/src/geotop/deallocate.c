
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 "Montebello" - 21 Jun 2011
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 "Montebello"
 
 GEOtop 1.145 "Montebello" is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 "Montebello" is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "struct.geotop.h"
#include "constants.h"
#include "snow.h"
#include "deallocate.h"
#include "../libraries/ascii/init.h"

extern char **files;
extern FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav, *ffice, *fficeav, *ffsnow, *ffglac;
extern double **outdata_point, *outdata_basin;
extern long *outputpoint, noutputpoint, *outputbasin, noutputbasin, *outputsnow, noutputsnow;
extern long *outputglac, noutputglac, *outputsoil, noutputsoil;
extern char **headerpoint, **headerbasin, **headersnow, **headerglac, **headersoil;
extern char *WORKING_DIRECTORY;
extern char *string_novalue;
extern long number_novalue;
extern long Nl, Nr, Nc;
extern T_INIT *UV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void dealloc_all(TOPO *top,SOIL *sl,LANDCOVER *land,WATER *wat,CHANNEL *cnet,PAR *par,ENERGY *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times){
	
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
	for (i=0; i<otot; i++) {
		if(par->state_pixel==1) free(outdata_point[i]);
		free(headerpoint[i]);
	}
	if(par->state_pixel==1) free(outdata_point);
	free(headerpoint);
	free(outputpoint);
	
	for (i=0; i<ootot; i++) {
		free(headerbasin[i]);
	}
	free(outdata_basin);
	free(headerbasin);
	free(outputbasin);
	
	for (j=0; j<12; j++) {
		free(headersnow[j]);
		free(headerglac[j]);
	}
	free(headersnow);
	free(headerglac);
	
	free(outputsnow);
	free(outputglac);
	
	for (j=0; j<6; j++) {
		free(headersoil[j]);
	}
	free(headersoil);
	
	free(outputsoil);
	
	free(WORKING_DIRECTORY);
	
	/* Deallocation of struct SOIL "sl": */
	printf("Deallocating soil\n");
	free_doubletensor(sl->P);
	free_doubletensor(sl->Ptot);
	free_doubletensor(sl->T);
	free_doubletensor(sl->T_av_tensor);
	free_doublematrix(sl->Tv);
	free_doubletensor(sl->thice); 
	free_doubletensor(sl->th); 
	free_longmatrix(sl->type);
	free_doubletensor(sl->pa);
	free_doubletensor(sl->ET);	
	
	if(par->state_pixel == 1){
		if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) free_doublematrix(sl->Tzplot);
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) free_doublematrix(sl->Tzavplot);
		if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) free_doublematrix(sl->Ptotzplot);
		if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) free_doublematrix(sl->Pzplot);
		if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) free_doublematrix(sl->thzplot);
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) free_doublematrix(sl->thzavplot);
		if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) free_doublematrix(sl->thicezplot);
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) free_doublematrix(sl->thicezavplot);
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
	free_longmatrix(top->Rdown);
	free_longmatrix(top->Cdown);
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
	free_doublevector(top->BC_LatDistance);
	free_doublevector(top->BC_DepthFreeSurface);
	
	if (par->point_sim==1) {
		free_doublematrix(top->latitude);
		free_doublematrix(top->longitude);
	}
	
	free(top);
	
	
	/* Deallocation of struct LANDCOVER "land": */
	printf("Deallocating land\n");
	free_doublematrix(land->LC);
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
	free_doubletensor(wat->PrecTot);
	free_doublematrix(wat->Pnet);
	free_doublematrix(wat->wcan_rain);
	free_doublematrix(wat->wcan_snow);
	
	if (par->output_meteo>0){
		free_doublematrix(wat->PrTOT_mean);
		free_doublematrix(wat->PrSNW_mean);
	}
	
	free_doublematrix(wat->h_sup);
	free_doublevector(wat->Lx);
	free_doublevector(wat->P0);
	free_doublevector(wat->H0);
	free_doublevector(wat->H1);
	free_doublevector(wat->dH);
	free_doublevector(wat->B);
	free_doublevector(wat->f);
	free_doublevector(wat->df);
	free_doublematrix(wat->Klat);
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
	free_doublematrix(cnet->P);
	free_doublematrix(cnet->T);
	free_doublematrix(cnet->th);
	free_doublematrix(cnet->thice);
	free_doublematrix(cnet->ET);
	free_doublevector(cnet->Tgskin);
	free(cnet);
	
	/* Deallocation of struct FILENAMES "filenames": */
	printf("Deallocating files\n"); 
	for (i=0; i<nfiles; i++) {
		free(files[i]);
	}
	free(files);
	
	/* Deallocation of struct T_INIT "UV": */
	printf("Deallocating UV\n"); 
	free_doublevector(UV->U);
	free_doublevector(UV->V);
	free(UV);
	
	/* Deallocation of struct ENERGY "egy": */
	printf("Deallocating egy\n");  
	free_doublevector(egy->hsun);
	free_doublevector(egy->dsun);
	free_doublevector(egy->sinhsun);
	if(par->output_surfenergy>0){
		free_doublematrix(egy->Rn_mean);
		if(par->distr_stat==1)free_doublematrix(egy->Rn_max);	
		if(par->distr_stat==1)free_doublematrix(egy->Rn_min);
		if(par->distr_stat==1)free_doublematrix(egy->LW_max);
		if(par->distr_stat==1)free_doublematrix(egy->LW_min);
		free_doublematrix(egy->LWin_mean);
		free_doublematrix(egy->LW_mean);
		free_doublematrix(egy->SW_mean);
		if(par->distr_stat==1)free_doublematrix(egy->SW_max);	
	}
	
	if(par->output_surfenergy>0){
		free_doublematrix(egy->ET_mean);
		if(par->distr_stat==1)free_doublematrix(egy->ET_max);
		if(par->distr_stat==1)free_doublematrix(egy->ET_min);
	}
	
	if(par->output_surfenergy>0){
		free_doublematrix(egy->H_mean);
		if(par->distr_stat==1)free_doublematrix(egy->H_max);
		if(par->distr_stat==1)free_doublematrix(egy->H_min);
	}
	
	if(par->output_surfenergy>0){
		free_doublematrix(egy->SEB_mean);
		if(par->distr_stat==1)free_doublematrix(egy->G_max);
		if(par->distr_stat==1)free_doublematrix(egy->G_min);
		free_doublematrix(egy->G_snowsoil);
	}
	
	if(par->output_surfenergy>0){
		free_doublematrix(egy->Ts_mean);
		if(par->distr_stat==1)free_doublematrix(egy->Ts_max);
		if(par->distr_stat==1)free_doublematrix(egy->Ts_min);
	}
	
	if(par->output_surfenergy>0){
		free_doublematrix(egy->Rswdown_mean);
		if(par->distr_stat==1)free_doublematrix(egy->Rswdown_max);
		free_doublematrix(egy->Rswbeam_mean);
	}
	
	free(egy->sun);
	free_longmatrix(egy->nDt_shadow);
	free_longmatrix(egy->nDt_sun); 
	
	if(times->JD_plots->nh>1){
		free_doublematrix(egy->Hgplot);
		free_doublematrix(egy->LEgplot);
		free_doublematrix(egy->Hvplot);
		free_doublematrix(egy->LEvplot);
		free_doublematrix(egy->SWinplot);
		free_doublematrix(egy->SWgplot);
		free_doublematrix(egy->SWvplot);
		free_doublematrix(egy->LWinplot);
		free_doublematrix(egy->LWgplot);
		free_doublematrix(egy->LWvplot);
		free_doublematrix(egy->Tsplot);
		free_doublematrix(egy->Tgplot);
		free_doublematrix(egy->Tvplot);
	}
	
	free_doublevector(egy->Dlay);
	free_doublevector(egy->wliq);
	free_doublevector(egy->wice);
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
	free_doublematrix(egy->Tgskin);
	free_doublematrix(egy->Tgskinsurr);
	free_doublematrix(egy->Asurr);	
	
	free(egy);
	
	
	/* Deallocation of struct SNOW "snow": */
	printf("Deallocating snow\n"); 
	if(times->JD_plots->nh > 1) free_doublematrix(snow->Dplot);
	
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
		
		if(par->output_snow>0){
			free_doublematrix(snow->Wtrans_plot);
			free_doublematrix(snow->Wsubl_plot);
			/*free_doublematrix(snow->Qsub_plot);
			 free_doublematrix(snow->Qsub_eq_plot);
			 free_doublematrix(snow->Qtrans_plot);
			 free_doublematrix(snow->Qtrans_eq_plot);*/
		}
	}
	free_doublematrix(snow->nondimens_age);
	free_doublematrix(snow->dimens_age); 
	
	if(par->output_snow>0){
		free_doublematrix(snow->t_snow);
		free_doublematrix(snow->MELTED);
		free_doublematrix(snow->SUBL);
	}
	
	if(par->blowing_snow==1) free_longvector(snow->change_dir_wind);
	free(snow);
	
	printf("Deallocating glacier\n");
	if(par->glaclayer_max>0){
		deallocate_statevar_3D(glac->G);
		if(par->output_glac>0) free_doublematrix(glac->MELTED);
		if(par->output_glac>0)	free_doublematrix(glac->SUBL);
	}
	free(glac);
	
	printf("Deallocating met\n"); 
	free_doublevector(met->tau_cloud);
	free_doublevector(met->tau_cloud_av);
	free_shortvector(met->tau_cloud_yes);
	free_shortvector(met->tau_cloud_av_yes);	
	free_doubletensor(met->Tgrid);
	free_doubletensor(met->Pgrid); 
	free_doubletensor(met->Vgrid);
	free_doubletensor(met->Vdir);
	free_doubletensor(met->RHgrid);
	if(par->output_meteo>0){
		free_doublematrix(met->Ta_mean);
		if(par->distr_stat==1)free_doublematrix(met->Ta_max);
		if(par->distr_stat==1)free_doublematrix(met->Ta_min);
		free_doublematrix(met->Vspdmean);
		free_doublematrix(met->Vdirmean);
		free_doublematrix(met->RHmean);
	}		
	if(times->JD_plots->nh > 1){
		free_doublematrix(met->Taplot);
		free_doublematrix(met->Vspdplot);
		free_doublematrix(met->Vdirplot);
		free_doublematrix(met->RHplot);
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
	free_doublevector(par->Dmin);
	free_doublevector(par->Dmax);
	if(par->glaclayer_max>0){
		free_doublevector(par->Dmin_glac);
		free_doublevector(par->Dmax_glac); 
	}
	
	if(par->state_pixel == 1){
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
		
	free(par);
	
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
	free(st);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void reset_to_zero(PAR *par, SOIL *sl, LANDCOVER *land, SNOW *snow, GLACIER *glac, ENERGY *egy, METEO *met, WATER *wat){
	
	long r, c, l, i, j;
	
	if(par->state_pixel == 1){
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) initialize_doublematrix(sl->Tzavplot,0.); 
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) initialize_doublematrix(sl->thzavplot,0.); 
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) initialize_doublematrix(sl->thicezavplot,0.); 
		
		for(i=1;i<=par->rc->nrh;i++){
			for(j=0;j<otot;j++) { 
				outdata_point[j][i-1]=0.0; 
			}
		}	
	}
	
	if(par->state_basin== 1){
		for(j=0;j<ootot;j++){ 
			outdata_basin[j]=0.0; 
		}
	}
	
	if(par->output_soil>0){
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0){
			for(r=1; r<=Nr; r++){
				for (c=1; c<=Nc; c++) {
					if( (long)land->LC->co[r][c] != number_novalue){
						for (l=1; l<=Nl; l++) {
							sl->T_av_tensor->co[l][r][c] = 0.0;
						}
					}
				}
			}
		}
	}
	
	if(par->output_snow>0){
		if(strcmp(files[fsnowmelt] , string_novalue) != 0) initmatrix(0.0, snow->MELTED, land->LC, number_novalue);	
		if(strcmp(files[fsnowsubl] , string_novalue) != 0) initmatrix(0.0, snow->SUBL, land->LC, number_novalue);	
		if(strcmp(files[fswe] , string_novalue) != 0 && par->blowing_snow==1){
			initmatrix(0.0, snow->Wtrans_plot, land->LC, number_novalue);
			initmatrix(0.0, snow->Wsubl_plot, land->LC, number_novalue);
		}
		if(strcmp(files[fsndur] , string_novalue) != 0) initmatrix(0.0, snow->t_snow, land->LC, number_novalue);
	}
	
	if(par->glaclayer_max>0 && par->output_glac>0){
		if(strcmp(files[fglacmelt] , string_novalue) != 0) initmatrix(0.0, glac->MELTED, land->LC, number_novalue);
		if(strcmp(files[fglacsubl] , string_novalue) != 0) initmatrix(0.0, glac->SUBL, land->LC, number_novalue);
	}
	
	if(par->output_surfenergy>0){
		if(strcmp(files[frad] , string_novalue) != 0){
			initmatrix(0.0, egy->Rn_mean, land->LC, number_novalue);
			initmatrix(0.0, egy->LWin_mean, land->LC, number_novalue);
			initmatrix(0.0, egy->LW_mean, land->LC, number_novalue);
			initmatrix(0.0, egy->SW_mean, land->LC, number_novalue);
			initmatrix(0.0, egy->Rswdown_mean, land->LC, number_novalue);
			initmatrix(0.0, egy->Rswbeam_mean, land->LC, number_novalue);			
			initlongmatrix(0, egy->nDt_shadow, land->LC, number_novalue);	
			initlongmatrix(0, egy->nDt_sun, land->LC, number_novalue);		
			if(par->distr_stat==1){
				initmatrix(-1.0E+9, egy->Rn_max, land->LC, number_novalue);
				initmatrix(1.0E+9, egy->Rn_min, land->LC, number_novalue);
				initmatrix(-1.0E+9, egy->LW_max, land->LC, number_novalue);
				initmatrix(1.0E+9, egy->LW_min, land->LC, number_novalue);
				initmatrix(-1.0E+9, egy->SW_max, land->LC, number_novalue);	
				initmatrix(-1.0E+9, egy->Rswdown_max, land->LC, number_novalue);	
			}
		}
		if(strcmp(files[fG] , string_novalue) != 0){
			initmatrix(0.0, egy->SEB_mean, land->LC, number_novalue);
			if(par->distr_stat==1){
				initmatrix(-1.0E+9, egy->G_max, land->LC, number_novalue);
				initmatrix(1.0E+9, egy->G_min, land->LC, number_novalue);
			}
		}
		if(strcmp(files[fH] , string_novalue) != 0){
			initmatrix(0.0, egy->H_mean, land->LC, number_novalue);	
			if(par->distr_stat==1){
				initmatrix(-1.0E+9, egy->H_max, land->LC, number_novalue);
				initmatrix(1.0E+9, egy->H_min, land->LC, number_novalue);
			}
		}
		if(strcmp(files[fLE] , string_novalue) != 0){
			initmatrix(0.0, egy->ET_mean, land->LC, number_novalue);	
			if(par->distr_stat==1){			
				initmatrix(-1.0E+9, egy->ET_max, land->LC, number_novalue);			
				initmatrix(1.0E+9, egy->ET_min, land->LC, number_novalue);	
			}
		}
		if(strcmp(files[fTs] , string_novalue) != 0){
			initmatrix(0.0, egy->Ts_mean, land->LC, number_novalue);
			if(par->distr_stat==1){
				initmatrix(-1.0E+9, egy->Ts_max, land->LC, number_novalue);
				initmatrix(1.0E+9, egy->Ts_min, land->LC, number_novalue);
			}
		}
	}
	
	if(par->output_meteo>0){
		if(strcmp(files[fTa] , string_novalue) != 0){				
			initmatrix(0.0, met->Ta_mean, land->LC, number_novalue);	
			if(par->distr_stat==1){
				initmatrix(-1.0E+9, met->Ta_max, land->LC, number_novalue);
				initmatrix(1.0E+9, met->Ta_min, land->LC, number_novalue);	
			}
		}
		if(strcmp(files[fprec] , string_novalue) != 0){
			initmatrix(0.0, wat->PrTOT_mean, land->LC, number_novalue);
			initmatrix(0.0, wat->PrSNW_mean, land->LC, number_novalue);
		}
		if(strcmp(files[fwspd] , string_novalue) != 0) initmatrix(0.0, met->Vspdmean, land->LC, number_novalue);	
		if(strcmp(files[fwdir] , string_novalue) != 0) initmatrix(0.0, met->Vdirmean, land->LC, number_novalue);
		if(strcmp(files[frh] , string_novalue) != 0) initmatrix(0.0, met->RHmean, land->LC, number_novalue);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

