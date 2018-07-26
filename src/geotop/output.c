
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
#include "output.h"
#include "pedo.funct.h"
#include "networks.h"
#include "rw_maps.h"
#include "constants.h"
#include "extensions.h"
#include "times.h"
#include "energy.balance.h"
#include "input.h"
#include "tabs.h"
#include "vegetation.h"
#include "tables.h"
#include "snow.h"
#include "init.h"
#include "water.balance.h"
#include "indices.h"
#include "recovering.h"


#include <time.h>
 
extern long number_novalue;
extern char *string_novalue;

extern T_INIT *UV;
extern char **files, *logfile;
extern long Nl, Nr, Nc;

extern double t_meteo, t_energy, t_water, t_sub, t_sup, t_out, t_blowingsnow;

extern double **odpnt, **odp, *odbsn, *odb;
extern long *opnt, nopnt, *obsn, nobsn, *osnw, nosnw;
extern long *oglc, noglc, *osl, nosl;
extern char **hpnt, **hbsn, **hsnw, **hglc, **hsl;

extern char *FailedRunFile;

extern FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav, *ffice, *fficeav, *ffsnowT, *ffsnowl, *ffsnowi, *ffsnowd, *ffglac;

extern long i_sim, i_run;

extern time_t start_time;
extern double elapsed_time, elapsed_time_start, cum_time, max_time;

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_output(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac, METEO *met)
				  
{
	/*internal auxiliary variables:*/
	long i,j,r=0,c=0,l,m,sy; /*counters*/
	long n_file;      /*number of file of the type "TETAxySSSlZZ"(i.e. number of the basin-time-step)*/
	char NNNNN[ ]={"NNNNN"};
	char RRRRR[ ]={"RRRRR"};
	char SSSSS[ ]={"SSSSS"};
	char NNNN[ ]={"NNNN"};
	char rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	char *name, *temp1, *temp2, *s1, *s2;  
	FILE *f=NULL, *flog;
	
	//time variables
	time_t stop_time;
	double percent_done=0., remaining_time, total_time;
	short first_column;
	double JD, JDfrom0;
	long day, month, year, hour, minute;
	
	// static double Qsub_ch, Qsup_ch;
	static long isavings;
	static double mass_error_tot;
	static double t_discharge, t_basin, t_point, t_rec;
	double Vchannel, Vsub, Vsup;
		
	//other variables	
	DOUBLEVECTOR *V;
	DOUBLEMATRIX *M;
	double D, Dthaw, cosslope;
	
	
	//initialize static variables
	if(times->time < 1.E-5){
		mass_error_tot=0.;
		t_discharge=0.;
		t_basin=0.;
		t_point=0.;
		t_rec=0.;
	}
	
	write_suffix(SSSSS, i_sim, 1);
	write_suffix(RRRRR, i_run, 1);
	 	
	//Time indices
	JDfrom0 = convert_tfromstart_JDfrom0(times->time+par->Dt, par->init_date->co[i_sim]);
	convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
	convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);
	cum_time += par->Dt;
	
	//DISCHARGE
	//****************************************************************************************************************
	//****************************************************************************************************************
	
	if(par->state_discharge == 1 && par->Dtplot_discharge->co[i_sim] > 1.E-5 && strcmp(files[fQ] , string_novalue) != 0){

		t_discharge += par->Dt;
				
		if (fabs(t_discharge - par->Dtplot_discharge->co[i_sim]) < 1.E-5){
			
			Vchannel = 0.;
			Vsub = 0.;
			Vsup = 0.;
			for(l=1;l<=par->total_channel;l++){
				r = cnet->r->co[l];
				c = cnet->c->co[l];
				Vchannel += 1.E-3 * Fmax(cnet->SS->P->co[0][l], 0.) / cos(top->slope->co[r][c]*Pi/180.) * UV->U->co[1] * par->w_dx * cnet->length->co[l];
				Vsub += cnet->Vsub->co[l];
				Vsup += cnet->Vsup->co[l];
			}
				
			if (par->recover > 0) write_suffix(rec, par->recover, 4);
			if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

			if (par->recover>0) {
				temp1 = join_strings(files[fQ], rec);
				name = join_strings(temp1, textfile);
				free(temp1);
			}else if (par->n_ContRecovery>0) {
				temp1 = join_strings(files[fQ], crec);
				name = join_strings(temp1, textfile);
				free(temp1);
			}else {
				name = join_strings(files[fQ], textfile);
			}
			
			f=fopen(name,"a");

			fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
			fprintf(f,",%f,%f,%f",(times->time+par->Dt)/secinday+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]),JDfrom0,JD);   
			fprintf(f,",%e,%e,%e,%e,%e,%e,%e\n",cnet->Vout/(double)par->Dtplot_discharge->co[i_sim],Vsup/(double)par->Dtplot_discharge->co[i_sim],
					Vsub/(double)par->Dtplot_discharge->co[i_sim],Vchannel,wat->Voutlandsup/(double)par->Dtplot_discharge->co[i_sim],
					wat->Voutlandsub/(double)par->Dtplot_discharge->co[i_sim],wat->Voutbottom/(double)par->Dtplot_discharge->co[i_sim]);
			fclose(f);
			free(name);
			
			t_discharge = 0.0;
			
			initialize_doublevector(cnet->Vsub, 0.);
			initialize_doublevector(cnet->Vsup, 0.);
			
			cnet->Vout = 0.;
			wat->Voutbottom = 0.;
			wat->Voutlandsub = 0.;
			wat->Voutlandsup = 0.;
			
		}
	}
	

	//****************************************************************************************************************
	//****************************************************************************************************************

	//DATA POINT
	if(par->Dtplot_point->co[i_sim] > 1.E-5){
		
		t_point += par->Dt;
		
		if(par->state_pixel == 1){
			
			for(i=1;i<=par->rc->nrh;i++){
				
				r=par->rc->co[i][1];
				c=par->rc->co[i][2];	
				j=top->j_cont[r][c];
				
				for(l=1;l<=Nl;l++){
					if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] += sl->SS->T->co[l][j]*(par->Dt/par->Dtplot_point->co[i_sim]);
					if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] += sl->th->co[l][j]*(par->Dt/par->Dtplot_point->co[i_sim]);
					if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot->co[i][l] += sl->SS->thi->co[l][j]*(par->Dt/par->Dtplot_point->co[i_sim]);
				}
				
				D = find_activelayerdepth_up(j, sl->type->co[r][c], sl);
				odpnt[othawedup][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);
				Dthaw = D;
				
				D = find_activelayerdepth_dw(j, sl->type->co[r][c], sl);
				odpnt[othaweddw][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);

				D = find_watertabledepth_up(Dthaw, j, sl->type->co[r][c], sl);
				odpnt[owtableup][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);
				
				D = find_watertabledepth_dw(Dthaw, j, sl->type->co[r][c], sl);
				odpnt[owtabledw][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);
			}
		}
			
		//Print of pixel-output every times->n_pixel time step
		if (fabs(t_point - par->Dtplot_point->co[i_sim])<1.E-5){
			
			if(par->state_pixel == 1){
				
				if (par->recover > 0) write_suffix(rec, par->recover, 4);
				if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);
			
				for(i=1;i<=par->rc->nrh;i++){
					
					write_suffix(NNNN, par->IDpoint->co[i], 0);
					r=par->rc->co[i][1];
					c=par->rc->co[i][2];
					j=top->j_cont[r][c];
					sy = sl->type->co[r][c];
					
					if (par->output_vertical_distances == 1) {
						cosslope = cos( Fmin(max_slope, top->slope->co[r][c]) * Pi/180. );
					}else {
						cosslope = 1.;
					}
										
					//soil data
					for(l=1;l<=Nl;l++){
						if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot->co[i][l] = sl->SS->T->co[l][j];
						if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot->co[i][l] = sl->Ptot->co[l][j];
						if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot->co[i][l] = sl->th->co[l][j];
						if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thizplot->co[i][l] = sl->SS->thi->co[l][j];
						if(strcmp(files[fsatz], string_novalue) != 0) sl->satratio->co[i][l] = (sl->SS->thi->co[l][j] + sl->th->co[l][j] - sl->pa->co[sy][jres][l]) / (sl->pa->co[sy][jsat][l] - sl->pa->co[sy][jres][l]);

						
					}			
					for(l=0;l<=Nl;l++){
						if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot->co[i][l] = sl->SS->P->co[l][j];
					}
					
					//snow data
					if(snow->S->lnum->co[r][c]>0){
						odpnt[osnowdepth][i-1] = 0.0;
						odpnt[oSWE][i-1] = 0.0;
						odpnt[osnowT][i-1] = 0.0;
						for(l=1;l<=snow->S->lnum->co[r][c];l++){
							odpnt[osnowdepth][i-1] += snow->S->Dzl->co[l][r][c];			   
							odpnt[oSWE][i-1] += 1.0E+3*(snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c])/rho_w;
							odpnt[osnowT][i-1] += snow->S->T->co[l][r][c]*snow->S->Dzl->co[l][r][c];
						}
						odpnt[osnowdens][i-1] = odpnt[oSWE][i-1]*rho_w/odpnt[osnowdepth][i-1];		
						odpnt[osnowT][i-1] /= odpnt[osnowdepth][i-1];   
					}else{
						odpnt[osnowdepth][i-1] = 0.0;
						odpnt[oSWE][i-1] = 0.0;
						odpnt[osnowdens][i-1] = (double)number_novalue;
						odpnt[osnowT][i-1] = (double)number_novalue;
					}
					
					//glacier data
					if(par->max_glac_layers>0){
						if(glac->G->lnum->co[r][c]>0){
							odpnt[oglacdepth][i-1] = 0.0;
							odpnt[oGWE][i-1] = 0.0;
							odpnt[oglacT][i-1] = 0.0;
							for(l=1;l<=glac->G->lnum->co[r][c];l++){
								odpnt[oglacdepth][i-1] += glac->G->Dzl->co[l][r][c];			   
								odpnt[oGWE][i-1] += 1.0E+3*(glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c])/rho_w;
								odpnt[oglacT][i-1] += glac->G->T->co[l][r][c]*glac->G->Dzl->co[l][r][c];
							}
							odpnt[oglacdens][i-1] = odpnt[oGWE][i-1]*rho_w/odpnt[oglacdepth][i-1];		
							odpnt[oglacT][i-1] /= odpnt[oglacdepth][i-1];   
						}else{
							odpnt[oglacdepth][i-1] = 0.0;
							odpnt[oGWE][i-1] = 0.0;
							odpnt[oglacdens][i-1] = (double)number_novalue;
							odpnt[oglacT][i-1] = (double)number_novalue;
						}
					}
					
					//Point data	
					
					if(strcmp(files[fpoint] , string_novalue) != 0){
						
						temp1=join_strings(files[fpoint],NNNN);
						
						if (par->recover>0) {
							temp2 = join_strings(temp1, rec);
							name = join_strings(temp2, textfile);
							free(temp2);
						}else if (par->n_ContRecovery>0) {
							temp2 = join_strings(temp1, crec);
							name = join_strings(temp2, textfile);
							free(temp2);
						}else {
							name = join_strings(temp1, textfile);
						}
						
						free(temp1);
						
						f=fopen(name,"a");				
						first_column=1;
						for (j=0; j<nopnt; j++) {
							if(first_column==0){
								fprintf(f,",");
							}else {
								first_column = 0;
							}
							if (opnt[j] >= 0) {
								if (opnt[j] == odate12) {
									fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (opnt[j] == oJDfrom0) {
									fprintf(f, "%f",JDfrom0);
								}else if (opnt[j] == odaysfromstart) {
									fprintf(f, "%f",(JDfrom0-par->init_date->co[i_sim])+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]));
								}else if (opnt[j] == operiod) {
									fprintf(f, "%ld",i_sim);
								}else if (opnt[j] == orun) {
									fprintf(f, "%ld",i_run);
								}else if (opnt[j] == opoint) {
									fprintf(f, "%ld",par->IDpoint->co[i]);
								}else {
									fprintf(f,"%f",odpnt[opnt[j]][i-1]);
								}
							}else {
								fprintf(f,"%ld",number_novalue);
							}
						}
						fprintf(f,"\n");
						fclose(f);
						free(name);			
					}
					
					if(strcmp(files[fpointwriteend] , string_novalue) != 0){
						
						first_column=1;
						for (j=0; j<nopnt; j++) {
							if(first_column==0){
								fprintf(ffpoint,",");
							}else {
								first_column = 0;
							}
							if (opnt[j] >= 0) {
								if (opnt[j] == odate12) {
									fprintf(ffpoint,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (opnt[j] == oJDfrom0) {
									fprintf(ffpoint, "%f",JDfrom0);
								}else if (opnt[j] == odaysfromstart) {
									fprintf(ffpoint, "%f",(JDfrom0-par->init_date->co[i_sim])+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]));
								}else if (opnt[j] == operiod) {
									fprintf(ffpoint, "%ld",i_sim);
								}else if (opnt[j] == orun) {
									fprintf(ffpoint, "%ld",i_run);		
								}else if (opnt[j] == opoint) {
									fprintf(ffpoint, "%ld",par->IDpoint->co[i]);
								}else {
									fprintf(ffpoint,"%f",odpnt[opnt[j]][i-1]);
								}
							}else {
								fprintf(ffpoint,"%ld",number_novalue);
							}
						}
						fprintf(ffpoint,"\n");
					}
															
					//Glacier
					if(par->max_glac_layers>0){
						
						if(strcmp(files[fglz] , string_novalue) != 0){
							
							temp1=join_strings(files[fglz],NNNN);
							
							if (par->recover>0) {
								temp2 = join_strings(temp1, rec);
								name = join_strings(temp2, textfile);
								free(temp2);
							}else if (par->n_ContRecovery>0) {
								temp2 = join_strings(temp1, crec);
								name = join_strings(temp2, textfile);
								free(temp2);
							}else {
								name = join_strings(temp1, textfile);
							}
							
							free(temp1);
							
							f=fopen(name,"a");
							
							if ((long)par->glac_plot_depths->co[1] != number_novalue) {
								m = par->glac_plot_depths->nh;
							}else {
								m = par->max_glac_layers;
							}										
						
							first_column=1;
							for (j=0; j<noglc; j++) {
								if(first_column==0){
									fprintf(f,",");
								}else {
									first_column = 0;
								}
								if (oglc[j] >= 0) {
									if (oglc[j] == 0) {
										fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
									}else if (oglc[j] == 1) {
										fprintf(f, "%f",JDfrom0);
									}else if (oglc[j] == 2) {
										fprintf(f, "%f",(JDfrom0-par->init_date->co[i_sim])+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]));
									}else if (oglc[j] == 3) {
										fprintf(f, "%ld",i_sim);
									}else if (oglc[j] == 4) {
										fprintf(f, "%ld",i_run);																
									}else if (oglc[j] == 5) {
										fprintf(f, "%ld",par->IDpoint->co[i]);
									}else if (oglc[j] <= 5 + 1*m) {
										l = oglc[j] - 5 - 0*m;
										if ((long)par->glac_plot_depths->co[1] != number_novalue) {
											fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->T, 0.));
										}else {
											fprintf(f, "%f", glac->G->T->co[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 2*m) {
										l = oglc[j] - 5 - 1*m;
										if ((long)par->glac_plot_depths->co[1] != number_novalue) {
											fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_ice, 0.));
										}else {
											fprintf(f, "%f",glac->G->w_ice->co[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m) {
										l = oglc[j] - 5 - 2*m;
										if ((long)par->glac_plot_depths->co[1] != number_novalue) {
											fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_liq, 0.));
										}else {
											fprintf(f, "%f",glac->G->w_liq->co[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m + par->max_glac_layers) {
										l = oglc[j] - 5 - 3*m;
										fprintf(f, "%f",glac->G->Dzl->co[l][r][c]);
									}
								}else {
									fprintf(f,"%f",(double)number_novalue);
								}
							}
							fprintf(f,"\n");
							fclose(f);
							free(name);
						}

						if(strcmp(files[fglzwriteend] , string_novalue) != 0){
							if ((long)par->glac_plot_depths->co[1] != number_novalue) {
								m = par->glac_plot_depths->nh;
							}else {
								m = par->max_glac_layers;
							}			
							first_column=1;
							for (j=0; j<noglc; j++) {
								if(first_column==0){
									fprintf(f,",");
								}else {
									first_column = 0;
								}
								if (oglc[j] >= 0) {
									if (oglc[j] == 0) {
										fprintf(ffglac,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
									}else if (oglc[j] == 1) {
										fprintf(ffglac, "%f",JDfrom0);
									}else if (oglc[j] == 2) {
										fprintf(f, "%f",(JDfrom0-par->init_date->co[i_sim])+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]));
									}else if (oglc[j] == 3) {
										fprintf(ffglac, "%ld",i_sim);
									}else if (oglc[j] == 4) {
										fprintf(ffglac, "%ld",i_run);																
									}else if (oglc[j] == 5) {
										fprintf(ffglac, "%ld",par->IDpoint->co[i]);
									}else if (oglc[j] <= 5 + 1*m) {
										l = oglc[j] - 5 - 0*m;
										if ((long)par->glac_plot_depths->co[1] != number_novalue) {
											fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->T, 0.));
										}else {
											fprintf(ffglac, "%f", glac->G->T->co[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 2*m) {
										l = oglc[j] - 5 - 1*m;
										if ((long)par->glac_plot_depths->co[1] != number_novalue) {
											fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_ice, 0.));
										}else {
											fprintf(ffglac, "%f",glac->G->w_ice->co[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m) {
										l = oglc[j] - 5 - 2*m;
										if ((long)par->glac_plot_depths->co[1] != number_novalue) {
											fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_liq, 0.));
										}else {
											fprintf(ffglac, "%f",glac->G->w_liq->co[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m + par->max_glac_layers) {
										l = oglc[j] - 5 - 3*m;
										fprintf(ffglac, "%f",glac->G->Dzl->co[l][r][c]);
									}
								}else {
									fprintf(ffglac,"%f",(double)number_novalue);
								}
							}
							fprintf(ffglac,"\n");
						}
						
					}
					
					//sl output
					write_soil_output(i, par->IDpoint->co[i], par->init_date->co[i_sim], par->end_date->co[i_sim], JDfrom0, JD, day, month, year, hour, minute, par->soil_plot_depths, sl, par, (double)PsiMin, cosslope);
				
					//snow output
					write_snow_output(i, par->IDpoint->co[i], r, c, par->init_date->co[i_sim], par->end_date->co[i_sim], JDfrom0, JD, day, month, year, hour, minute, par->snow_plot_depths, snow->S, par, cosslope);					
					
					//initialize
					for(j=0;j<otot;j++) { odpnt[j][i-1]=0.0; }

				}
				
			}
			
			if(strcmp(files[fSCA] , string_novalue) != 0) find_SCA(snow->S, par, land->LC->co, times->time + par->Dt);			
			
			percent_done = 100.*cum_time/max_time;
									
			time( &stop_time );
			elapsed_time = difftime( stop_time, start_time ) + elapsed_time_start;
			
			if( percent_done > 1.0e-6 ){
				total_time = elapsed_time * 100.0 / percent_done;
			}else{
				total_time = 1.E5;
			}
			
			remaining_time = (total_time - elapsed_time);
			
			printf("%ld/%ld/%ld %ld:%02.0f %.2f%% - Times: Elapsed (h:m:s) %2.0f:%02.0f:%02.0f Remaining (h:m) %2.0f:%02.0f  \n",
				   day,month,year,hour,(float)minute,percent_done, 
				   floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
				   floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
				   floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
			
			flog = fopen(logfile, "a");
			fprintf(flog,"%ld/%ld/%ld %ld:%02.0f %.2f%% - Time elapsed (h:m:s) %2.0f:%02.0f:%02.0f Time remaining (h:m) %2.0f:%02.0f  \n",
					day,month,year,hour,(float)minute,percent_done,
					floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
					floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
					floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
			fclose(flog);
			
			t_point = 0.0;
			
		}
	}
	

	
	//BASIN DATA
	//****************************************************************************************************************
	//****************************************************************************************************************
	if(par->Dtplot_basin->co[i_sim] > 1.E-5 && par->state_basin == 1){
		
		t_basin += par->Dt;
		
		if (fabs(t_basin - par->Dtplot_basin->co[i_sim])<1.E-5){
			
			if(strcmp(files[fbas] , string_novalue) != 0){
				
				if (par->recover > 0) write_suffix(rec, par->recover, 4);
				if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

				if (par->recover>0) {
					temp1 = join_strings(files[fbas], rec);
					name = join_strings(temp1, textfile);
					free(temp1);
				}else if (par->n_ContRecovery>0) {
					temp1 = join_strings(files[fbas], crec);
					name = join_strings(temp1, textfile);
					free(temp1);
				}else {
					name = join_strings(files[fbas], textfile);
				}
				
				f=fopen(name,"a");
				first_column=1;
				for (j=0; j<nobsn; j++) {
					if(first_column==0){
						fprintf(f,",");
					}else {
						first_column = 0;
					}
					if (obsn[j] >= 0) {
						if (obsn[j] == oodate12) {
							fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
						}else if (obsn[j] == ooJDfrom0) {
							fprintf(f, "%f",JDfrom0);
						}else if (obsn[j] == oodaysfromstart) {
							fprintf(f, "%f",(JDfrom0-par->init_date->co[i_sim])+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]));
						}else {
							fprintf(f,"%f",odbsn[obsn[j]]);
						}
					}else {
						fprintf(f,"%f",(double)number_novalue);
					}
				}
				fprintf(f,"\n");
				fclose(f);
				free(name);
			}
			
			if(strcmp(files[fbaswriteend] , string_novalue) != 0){
				first_column=1;
				for (j=0; j<nobsn; j++) {
					if(first_column==0){
						fprintf(ffbas,",");
					}else {
						first_column = 0;
					}
					if (obsn[j] >= 0) {
						if (obsn[j] == oodate12) {
							fprintf(ffbas,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
						}else if (obsn[j] == ooJDfrom0) {
							fprintf(ffbas, "%f",JDfrom0);
						}else if (obsn[j] == oodaysfromstart) {
							fprintf(f, "%f",(JDfrom0-par->init_date->co[i_sim])+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]));
						}else {
							fprintf(ffbas,"%f",odbsn[obsn[j]]);
						}
					}else {
						fprintf(ffbas,"%f",(double)number_novalue);
					}
				}
				fprintf(ffbas,"\n");
			}
			
			mass_error_tot += odbsn[oomasserror];

			printf("\n%ld/%ld/%ld %ld:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
				   day,month,year,hour,(float)minute,JD,(long)(floor(times->time/86400))+1,
				   percent_done);
			
			printf(" t_meteo:%6.2f s, t_energy:%6.2f s, t_blowingsnow:%6.2f s, t_water:%6.2f s, t_sub:%6.2f s, t_sup:%6.2f s, t_out:%6.2f s\n",t_meteo,t_energy,t_blowingsnow,t_water,t_sub,t_sup,t_out);
			printf(" SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Pvert=%6.2f mm Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm Mean Time Step=%f s\n\n",
				   odbsn[ooSW],odbsn[ooLW],odbsn[ooH],odbsn[ooLE],odbsn[oopnet],odbsn[oorainover],
				   odbsn[oosnowover],odbsn[oomasserror]*3600.0/par->Dtplot_basin->co[i_sim],mass_error_tot,odbsn[ootimestep]);
			
			flog = fopen(logfile, "a");
			fprintf(flog,"\n%ld/%ld/%ld %ld:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
					day,month,year,hour,(float)minute,JD,(long)(floor(times->time/86400))+1,
					percent_done);			
			fprintf(flog," t_meteo:%6.2f s, t_energy:%6.2f s, t_blowingsnow:%6.2f s, t_water:%6.2f s, t_sub:%6.2f s, t_sup:%6.2f s, t_out:%6.2f s\n",t_meteo,t_energy,t_blowingsnow,t_water,t_sub,t_sup,t_out);
			fprintf(flog," SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Pvert=%6.2f mm Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm Mean Time Step=%f s\n\n",
					odbsn[ooSW],odbsn[ooLW],odbsn[ooH],odbsn[ooLE],odbsn[oopnet],odbsn[oorainover],
					odbsn[oosnowover],odbsn[oomasserror]*3600.0/par->Dtplot_basin->co[i_sim],mass_error_tot,odbsn[ootimestep]);
			
			fclose(flog);
			
			
			for(j=0;j<ootot;j++){ 
				odbsn[j]=0.0; 
			}
			
			t_basin = 0.0;
			
		}
	}
	
		
	//DISTRIBUTED OUTPUTS
	//****************************************************************************************************************
	//****************************************************************************************************************
	//averaging properties
	if(par->output_meteo->co[i_sim]>0){
		if(strcmp(files[fTa] , string_novalue) != 0) {
			for(i=1;i<=par->total_pixel;i++){
				met->Tamean->co[i]+=met->Tgrid->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
			}
		}			
		if(strcmp(files[fwspd] , string_novalue) != 0){
			for(i=1;i<=par->total_pixel;i++){
				met->Vspdmean->co[i]+=met->Vgrid->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
			}
		}			
		if(strcmp(files[fwdir] , string_novalue) != 0){
			for(i=1;i<=par->total_pixel;i++){
				met->Vdirmean->co[i]+=met->Vdir->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
			}
		}			
		if(strcmp(files[frh] , string_novalue) != 0){
			for(i=1;i<=par->total_pixel;i++){
				met->RHmean->co[i]+=met->RHgrid->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
			}
		}			
	}			
		
	if (par->output_soil->co[i_sim]>0) {
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0){
			for (i=1; i<=par->total_pixel; i++) {
				for (l=1; l<=Nl; l++) {
					sl->T_av_tensor->co[l][i] += sl->SS->T->co[l][i]/((par->output_soil->co[i_sim]*3600.0)/(par->Dt));
				}
			}
		}
		if(strcmp(files[fliqav] , string_novalue) != 0){
			for (i=1; i<=par->total_pixel; i++) {
				for (l=1; l<=Nl; l++) {
					sl->thw_av_tensor->co[l][i] += sl->th->co[l][i]/((par->output_soil->co[i_sim]*3600.0)/(par->Dt));
				}
			}
		}
		if(strcmp(files[ficeav] , string_novalue) != 0){
			for (i=1; i<=par->total_pixel; i++) {
				for (l=1; l<=Nl; l++) {
					sl->thi_av_tensor->co[l][i] += sl->SS->thi->co[l][i]/((par->output_soil->co[i_sim]*3600.0)/(par->Dt));
				}
			}
		}
	}	
	
	
	V=new_doublevector(par->total_pixel);
	initialize_doublevector(V, (double)number_novalue);
		
	//soil properties
	if(par->output_soil->co[i_sim]>0 && fmod(times->time+par->Dt,par->output_soil->co[i_sim]*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt)/(par->output_soil->co[i_sim]*3600.0));
		write_suffix(NNNNN, n_file, 1);
		if (par->run_times->co[i_sim] == 1) {
			s1 = join_strings(NNNNN, "");
		}else {
			s1 = join_strings(NNNNN, RRRRR);
		}
		if (par->init_date->nh == 1) {
			s2 = join_strings(s1, "");
		}else {
			s2 = join_strings(s1, SSSSS);
		}
		free(s1);

		//theta liq tensor
		if(strcmp(files[fliq] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fliq], 0, par->format_out, sl->th, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fliq], 0, par->format_out, sl->th, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}

		//theta liq surface
		if(strcmp(files[fliqsup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
				V->co[i] = sl->th->co[1][i];
			}

			temp1=join_strings(files[fliqsup],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		//write thw_av tensor
		if(strcmp(files[fliqav] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fliqav], 0, par->format_out, sl->thw_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fliqav], 0, par->format_out, sl->thw_av_tensor, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//initialize thw_av_tensor
		if(strcmp(files[fliqav] , string_novalue) != 0) initialize_doublematrix(sl->thw_av_tensor, 0.);
		
		//write T tensor
		if(strcmp(files[fT] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fT], 0, par->format_out, sl->SS->T, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fT], 0, par->format_out, sl->SS->T, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//theta T surface
		if(strcmp(files[fTsup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
				V->co[i] = sl->SS->T->co[1][i];
			}
			
			temp1=join_strings(files[fTsup],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		//write Tav tensor
		if(strcmp(files[fTav] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fTav], 0, par->format_out, sl->T_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fTav], 0, par->format_out, sl->T_av_tensor, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//theta Tav surface
		if(strcmp(files[fTavsup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
				V->co[i] = sl->T_av_tensor->co[1][i];
			}

			temp1=join_strings(files[fTavsup],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		//initialize T_av_tensor
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0) initialize_doublematrix(sl->T_av_tensor, 0.);
		
		//theta_ice tensor
		if(strcmp(files[fice] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fice], 0, par->format_out, sl->SS->thi, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fice], 0, par->format_out, sl->SS->thi, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//theta_ice surface
		if(strcmp(files[ficesup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
				V->co[i] = sl->SS->thi->co[1][i];
			}

			temp1=join_strings(files[ficesup],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		//write thi_av tensor
		if(strcmp(files[ficeav] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[ficeav], 0, par->format_out, sl->thi_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[ficeav], 0, par->format_out, sl->thi_av_tensor, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//initialize thi_av_tensor
		if(strcmp(files[ficeav] , string_novalue) != 0) initialize_doublematrix(sl->thi_av_tensor, 0.);
		
		//write psi tensors
		if(strcmp(files[fpsitot] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fpsitot], 0, par->format_out, sl->Ptot, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fpsitot], 0, par->format_out, sl->Ptot, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		if(strcmp(files[fpsiliq] , string_novalue) != 0){
			if ((long)par->soil_plot_depths->co[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fpsiliq], 0, par->format_out, sl->SS->P, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa->co[1][jdz], top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fpsiliq], 0, par->format_out, sl->SS->P, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//calculate saturation front depth
		if( strcmp(files[fwtable_up] , string_novalue) != 0 ){

			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = find_watertabledepth_up(find_activelayerdepth_up(i, sl->type->co[r][c], sl), i, sl->type->co[r][c], sl);//normal
			}
			temp1=join_strings(files[fwtable_up],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}

		if( strcmp(files[fwtable_dw] , string_novalue) != 0 ){
			
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = find_watertabledepth_dw(find_activelayerdepth_up(i, sl->type->co[r][c], sl), i, sl->type->co[r][c], sl);//normal
			}
			temp1=join_strings(files[fwtable_dw],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		//calculate active layer depth		   
		if( strcmp(files[fthawed_up] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = find_activelayerdepth_up(i, sl->type->co[r][c], sl);   //normal
			}
			temp1=join_strings(files[fthawed_up],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}

		if( strcmp(files[fthawed_dw] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = find_activelayerdepth_dw(i, sl->type->co[r][c], sl);   //normal
			}
			temp1=join_strings(files[fthawed_dw],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		//WATER OVER THE SURFACE
		if( strcmp(files[fhsupland] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = Fmax(0, sl->SS->P->co[0][i]) / cos(top->slope->co[r][c] * Pi/180.);
			}
			
			temp1 = join_strings(files[fhsupland], s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
			
		if( strcmp(files[fhsupch] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				if (cnet->ch->co[r][c]!=0) {
					V->co[i] = cnet->SS->P->co[0][cnet->ch->co[r][c]] / cos(top->slope->co[r][c] * Pi/180.);
				}else {
					V->co[i] = (double)number_novalue;
				}
			}
			
			temp1 = join_strings(files[fhsupch], s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}

		if(strcmp(files[fpnet] , string_novalue) != 0){
			temp1=join_strings(files[fpnet], s2);
			write_map_vector(temp1, 0, par->format_out, sl->Pnetcum, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(sl->Pnetcum, 0.);
			free(temp1);
		}

		if(strcmp(files[fevap] , string_novalue) != 0){
			temp1=join_strings(files[fevap], s2);
			write_map_vector(temp1, 0, par->format_out, sl->ETcum, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(sl->ETcum, 0.);
			free(temp1);
		}
		
		free(s2);
	}
	
	//snow properties
	if(par->output_snow->co[i_sim]>0 && fmod(times->time+par->Dt,par->output_snow->co[i_sim]*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt)/(par->output_snow->co[i_sim]*3600.0));	
		write_suffix(NNNNN, n_file, 1);
		if (par->run_times->co[i_sim] == 1) {
			s1 = join_strings(NNNNN, "");
		}else {
			s1 = join_strings(NNNNN, RRRRR);
		}
		if (par->init_date->nh == 1) {
			s2 = join_strings(s1, "");
		}else {
			s2 = join_strings(s1, SSSSS);
		}
		free(s1);
		
		if(strcmp(files[fsnowdepth] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = 0.;
				for(l=1;l<=snow->S->lnum->co[r][c];l++){
					V->co[i] += snow->S->Dzl->co[l][r][c];
				}
			}
			
			temp1=join_strings(files[fsnowdepth],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		if(strcmp(files[fsnowmelt] , string_novalue) != 0){
			temp1=join_strings(files[fsnowmelt],s2);
			write_map_vector(temp1, 0, par->format_out, snow->MELTED, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(snow->MELTED, 0.);
			free(temp1);
		}
		
		if(strcmp(files[fsnowsubl] , string_novalue) != 0){
			temp1=join_strings(files[fsnowsubl],s2);
			write_map_vector(temp1, 0, par->format_out, snow->SUBL, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(snow->SUBL, 0.);
			free(temp1);
		}
		
		if(strcmp(files[fswe] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = 0.;
				for(l=1;l<=snow->S->lnum->co[r][c];l++){
					V->co[i] += (snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
				}
			}
			temp1=join_strings(files[fswe],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
			
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = 0.;
				D = 0.;
				for(l=1;l<=snow->S->lnum->co[r][c];l++){
					V->co[i] += (snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
					D += snow->S->Dzl->co[l][r][c];
				}
				V->co[i] /= (1.E-3*D);
			}
			temp1=join_strings(files[fswe], "DENSITY");
			temp2=join_strings(temp1,s2);
			write_map_vector(temp2, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp2);
			free(temp1);
			
			if(par->blowing_snow==1){
				temp1=join_strings(files[fswe],"WindTrans");
				temp2=join_strings(temp1, s2);
				write_map(temp2, 0, par->format_out, snow->Wtrans_plot, UV, number_novalue);
				initmatrix(0.0, snow->Wtrans_plot, land->LC, number_novalue);
				free(temp2);
				free(temp1);
				
				temp1=join_strings(files[fswe],"WindSubl");
				temp2=join_strings(temp1, s2);
				write_map(temp2, 0, par->format_out, snow->Wsubl_plot, UV, number_novalue);
				initmatrix(0.0, snow->Wsubl_plot, land->LC, number_novalue);
				free(temp2);
				free(temp1);
				
			}
		}
		
		if(strcmp(files[fsndur] , string_novalue) != 0){
			temp1=join_strings(files[fsndur],s2);
			write_map_vector(temp1, 0, par->format_out, snow->t_snow, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(snow->t_snow, 0.);
			free(temp1);
		}
		
		free(s2);
	}
	
	//glacier properties
	if(par->max_glac_layers>0 && par->output_glac->co[i_sim]>0 && fmod(times->time+par->Dt,par->output_glac->co[i_sim]*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt)/(par->output_glac->co[i_sim]*3600.0));	
		write_suffix(NNNNN, n_file, 1);
		if (par->run_times->co[i_sim] == 1) {
			s1 = join_strings(NNNNN, "");
		}else {
			s1 = join_strings(NNNNN, RRRRR);
		}
		if (par->init_date->nh == 1) {
			s2 = join_strings(s1, "");
		}else {
			s2 = join_strings(s1, SSSSS);
		}
		free(s1);
		
		if(strcmp(files[fglacdepth] , string_novalue) != 0){
			
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = 0.;
				for(l=1;l<=glac->G->lnum->co[r][c];l++){
					V->co[i] += (glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
				}
			}
			temp1=join_strings(files[fglacdepth],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}
		
		if(strcmp(files[fglacmelt] , string_novalue) != 0){
			temp1=join_strings(files[fglacmelt],s2);
			write_map_vector(temp1, 0, par->format_out, glac->MELTED, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(glac->MELTED, 0.);
			free(temp1);
		}
		
		if(strcmp(files[fglacsubl] , string_novalue) != 0){
			temp1=join_strings(files[fglacsubl],s2);
			write_map_vector(temp1, 0, par->format_out, glac->SUBL, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(glac->SUBL, 0.);
			free(temp1);
		}
		
		if(strcmp(files[fgwe] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
				r = top->rc_cont->co[i][1];
				c = top->rc_cont->co[i][2];
				V->co[i] = 0.;
				for(l=1;l<=glac->G->lnum->co[r][c];l++){
					V->co[i] += (glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
				}
			}
			temp1=join_strings(files[fgwe],s2);
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
		}

		free(s2);

	}
	
	//SURFACE ENERGY BALANCE	
	
	//RADIATION
	if(par->output_surfenergy->co[i_sim]>0 && fmod(times->time+par->Dt,par->output_surfenergy->co[i_sim]*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt)/(par->output_surfenergy->co[i_sim]*3600.0));	
		write_suffix(NNNNN, n_file, 1);
		if (par->run_times->co[i_sim] == 1) {
			s1 = join_strings(NNNNN, "");
		}else {
			s1 = join_strings(NNNNN, RRRRR);
		}
		if (par->init_date->nh == 1) {
			s2 = join_strings(s1, "");
		}else {
			s2 = join_strings(s1, SSSSS);
		}
		free(s1);
		
		if(strcmp(files[fradnet] , string_novalue) != 0){			
			name=join_strings(files[fradnet], s2);
			write_map_vector(name, 0, par->format_out, egy->Rn_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->Rn_mean, 0.);
			free(name);
		}
			
		if(strcmp(files[fradLWin] , string_novalue) != 0){			
			name=join_strings(files[fradLWin],s2);
			write_map_vector(name, 0, par->format_out, egy->LWin_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->LWin_mean, 0.);
			free(name);
		}
		
		if(strcmp(files[fradLW] , string_novalue) != 0){			
			name=join_strings(files[fradLW],s2);
			write_map_vector(name, 0, par->format_out, egy->LW_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->LW_mean, 0.);
			free(name);
		}
			
		if(strcmp(files[fradSW] , string_novalue) != 0){			
			name=join_strings(files[fradSW],s2);
			write_map_vector(name, 0, par->format_out, egy->SW_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->SW_mean, 0.);
			free(name);
		}
		
		if(strcmp(files[fradSWin] , string_novalue) != 0){			
			name=join_strings(files[fradSWin],s2);
			write_map_vector(name, 0, par->format_out, egy->Rswdown_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->Rswdown_mean, 0.);
			free(name);
		}
		
		if(strcmp(files[fradSWinbeam] , string_novalue) != 0){			
			name=join_strings(files[fradSWinbeam],s2);
			write_map_vector(name, 0, par->format_out, egy->Rswbeam_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->Rswbeam_mean, 0.);
			free(name);
		}
			
		if(strcmp(files[fshadow] , string_novalue) != 0){			

			for(i=1; i<=par->total_pixel; i++){
				if(egy->nDt_sun->co[i]>0){
					V->co[i] = egy->nDt_shadow->co[i]/(double)(egy->nDt_sun->co[i]);
				}else{
					V->co[i] = -1.;
				}
			}
			
			name=join_strings(files[fshadow],s2);
			write_map_vector(name, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_longvector(egy->nDt_shadow, 0.);
			initialize_longvector(egy->nDt_sun, 0.);
			free(name);
		}
					
		//GROUND HEAT FLUX
		if(strcmp(files[fG] , string_novalue) != 0){
			name=join_strings(files[fG],s2);
			write_map_vector(name, 0, par->format_out, egy->SEB_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->SEB_mean, 0.);
			free(name);
		}
				
		//SENSIBLE HEAT FLUX
		if(strcmp(files[fH] , string_novalue) != 0){
			name=join_strings(files[fH],s2);
			write_map_vector(name, 0, par->format_out, egy->H_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->H_mean, 0.);
			free(name);
		}
		
		
		//LATENT HEAT FLUX
		if(strcmp(files[fLE] , string_novalue) != 0){
			name=join_strings(files[fLE],s2);
			write_map_vector(name, 0, par->format_out, egy->ET_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->ET_mean, 0.);
			free(name);
		}
		
		//SURFACE TEMPERATURE
		if(strcmp(files[fTs] , string_novalue) != 0){
			name=join_strings(files[fTs],s2);
			write_map_vector(name, 0, par->format_out, egy->Ts_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(egy->Ts_mean, 0.);
			free(name);
		}
		
		free(s2);

	}
	
	//vegetation variables	
	if(par->output_vegetation->co[i_sim]>0 && fmod(times->time+par->Dt,par->output_vegetation->co[i_sim]*3600.0)<1.E-5){
		n_file=(long)((times->time+par->Dt)/(par->output_vegetation->co[i_sim]*3600.0));	
		write_suffix(NNNNN, n_file, 1);
		if (par->run_times->co[i_sim] == 1) {
			s1 = join_strings(NNNNN, "");
		}else {
			s1 = join_strings(NNNNN, RRRRR);
		}
		if (par->init_date->nh == 1) {
			s2 = join_strings(s1, "");
		}else {
			s2 = join_strings(s1, SSSSS);
		}
		free(s1);
		
		//INTERCEPTED PRECIPITATION
		if(strcmp(files[fcint] , string_novalue) != 0){
			
			temp1=join_strings(files[fcint],"water");
			temp2=join_strings(temp1,s2);
			write_map_vector(temp2, 0, par->format_out, sl->VS->wrain, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
			free(temp2);
			
			temp1=join_strings(files[fcint],"snow");
			temp2=join_strings(temp1,s2);
			write_map_vector(temp2, 0, par->format_out, sl->VS->wsnow, UV, number_novalue, top->j_cont, Nr, Nc);
			free(temp1);
			free(temp2);
		}
		
		free(s2);

	}
		
	//METEO
	if(par->output_meteo->co[i_sim]>0 && fmod(times->time+par->Dt,par->output_meteo->co[i_sim]*3600.0)<1.E-5){
		n_file=(long)((times->time+par->Dt)/(par->output_meteo->co[i_sim]*3600.0));	
		
		write_suffix(NNNNN, n_file, 1);
		if (par->run_times->co[i_sim] == 1) {
			s1 = join_strings(NNNNN, "");
		}else {
			s1 = join_strings(NNNNN, RRRRR);
		}
		if (par->init_date->nh == 1) {
			s2 = join_strings(s1, "");
		}else {
			s2 = join_strings(s1, SSSSS);
		}
		free(s1);
		
		//AIR TEMPERATURE
		if(strcmp(files[fTa] , string_novalue) != 0){
			name=join_strings(files[fTa],s2);
			write_map_vector(name, 0, par->format_out, met->Tamean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(met->Tamean, 0.);
			free(name);
		}
		
		//PRECIPITATION
		if(strcmp(files[fprec] , string_novalue) != 0){
			
			name=join_strings(files[fprec],"TOTAL");
			temp1=join_strings(name, s2);
			write_map_vector(temp1, 0, par->format_out, wat->PrTOT_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(wat->PrTOT_mean, 0.);
			free(temp1);
			free(name);
			
			name=join_strings(files[fprec],"SNOW");	
			temp1=join_strings(name, s2);
			write_map_vector(temp1, 0, par->format_out, wat->PrSNW_mean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(wat->PrSNW_mean, 0.);
			free(temp1);
			free(name);
		}
		
		if(strcmp(files[fwspd] , string_novalue) != 0){
			name=join_strings(files[fwspd],s2);
			write_map_vector(name, 0, par->format_out, met->Vspdmean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(met->Vspdmean, 0.);
			free(name);
		}
		
		if(strcmp(files[fwdir] , string_novalue) != 0){
			name=join_strings(files[fwdir],s2);
			write_map_vector(name, 0, par->format_out, met->Vdirmean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(met->Vdirmean, 0.);
			free(name);
		}
		
		if(strcmp(files[frh] , string_novalue) != 0){
			name=join_strings(files[frh],s2);
			write_map_vector(name, 0, par->format_out, met->RHmean, UV, number_novalue, top->j_cont, Nr, Nc);
			initialize_doublevector(met->RHmean, 0.);
			free(name);
		}
		
		free(s2);
	}	
	
	free_doublevector(V);
	
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	//SPECIAL PLOTS AT SOME DAYS
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	
	if(times->JD_plots->nh > 1 && times->iplot<=times->JD_plots->nh){
		
		i=times->iplot;
		j=2*i-1;
		if( fabs(par->init_date->co[i_sim]+(times->time+par->Dt)/86400. - times->JD_plots->co[j+1]) < 1.E-5 ){
			
			flog = fopen(logfile, "a");
			fprintf(flog,"Printing plot number %ld \n",i);
			fclose(flog);
			
			printf("Printing plot number %ld \n",i);
			
			V = new_doublevector(par->total_pixel);
			
			if(strcmp(files[pH] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
					V->co[i] = egy->Hgplot->co[i] + egy->Hvplot->co[i];
				}
				
				plot(files[pH], i, V, par->format_out, top->j_cont);
			}
			
			if(strcmp(files[pLE] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
					V->co[i] = egy->LEgplot->co[i] + egy->LEvplot->co[i];
				}
				
				plot(files[pLE], i, V, par->format_out, top->j_cont);
			}

			if(strcmp(files[pG] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
					V->co[i] = egy->SWgplot->co[i]+egy->LWgplot->co[i]-egy->Hgplot->co[i]-egy->LEgplot->co[i];
				}
				plot(files[pG], i, V, par->format_out, top->j_cont);
			}
			
			if(strcmp(files[pth] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
					V->co[i] = sl->th->co[1][i];
				}
				plot(files[pth], i, V, par->format_out, top->j_cont);
			}

			if(strcmp(files[pth] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
					V->co[i] = sl->th->co[1][i];
				}
				plot(files[pth], i, V, par->format_out, top->j_cont);
			}
			
			if(strcmp(files[pVspd] , string_novalue) != 0 ){
			   for (i=1; i<=par->total_pixel; i++) {
				   V->co[i] = sqrt(pow(met->Vxplot->co[i], 2.0) + pow(met->Vyplot->co[i], 2.0));
			   }
			   plot(files[pVspd], i, V, par->format_out, top->j_cont);
			}
			   
			if(strcmp(files[pVdir] , string_novalue) != 0){
				for (i=1; i<=par->total_pixel; i++) {
					V->co[i] = 270.0 - (180./Pi)*atan2(met->Vyplot->co[i],met->Vxplot->co[i]);
					if (V->co[i] >= 360.0) V->co[i] -= 360.0;
				}
				plot(files[pVdir], i, V, par->format_out, top->j_cont);
			}
			
			free_doublevector(V);
				
			if(strcmp(files[pHg] , string_novalue) != 0) plot(files[pHg], i, egy->Hgplot, par->format_out, top->j_cont);
			if(strcmp(files[pLEg] , string_novalue) != 0) plot(files[pLEg], i, egy->LEgplot, par->format_out, top->j_cont);
			if(strcmp(files[pHv] , string_novalue) != 0) plot(files[pHv], i, egy->Hvplot, par->format_out, top->j_cont);
			if(strcmp(files[pLEv] , string_novalue) != 0) plot(files[pLEv], i, egy->LEvplot, par->format_out, top->j_cont);
			if(strcmp(files[pSWin] , string_novalue) != 0) plot(files[pSWin], i, egy->SWinplot, par->format_out, top->j_cont);
			if(strcmp(files[pSWg] , string_novalue) != 0) plot(files[pSWg], i, egy->SWgplot, par->format_out, top->j_cont);
			if(strcmp(files[pSWv] , string_novalue) != 0) plot(files[pSWv], i, egy->SWvplot, par->format_out, top->j_cont);
			if(strcmp(files[pLWin] , string_novalue) != 0) plot(files[pLWin], i, egy->LWinplot, par->format_out, top->j_cont);
			if(strcmp(files[pLWg] , string_novalue) != 0) plot(files[pLWg], i, egy->LWgplot, par->format_out, top->j_cont);
			if(strcmp(files[pLWv] , string_novalue) != 0) plot(files[pLWv], i, egy->LWvplot, par->format_out, top->j_cont);
			if(strcmp(files[pTs] , string_novalue) != 0) plot(files[pTs], i, egy->Tsplot, par->format_out, top->j_cont);
			if(strcmp(files[pTg] , string_novalue) != 0) plot(files[pTg], i, egy->Tgplot, par->format_out, top->j_cont);
			if(strcmp(files[pTv] , string_novalue) != 0) plot(files[pTv], i, egy->Tvplot, par->format_out, top->j_cont);
			if(strcmp(files[pTa] , string_novalue) != 0) plot(files[pTa], i, met->Taplot, par->format_out, top->j_cont);
			if(strcmp(files[pD] , string_novalue) != 0) plot(files[pD], i, snow->Dplot, par->format_out, top->j_cont);
			if(strcmp(files[pRH] , string_novalue) != 0) plot(files[pRH], i, met->RHplot, par->format_out, top->j_cont);
						
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->Hgplot, 0.);
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->LEgplot, 0.);
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0) initialize_doublevector(egy->Hvplot, 0.);
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0) initialize_doublevector(egy->LEvplot, 0.);
			if(strcmp(files[pSWin] , string_novalue) != 0) initialize_doublevector(egy->SWinplot, 0.);
			if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->SWgplot, 0.);
			if(strcmp(files[pSWv] , string_novalue) != 0) initialize_doublevector(egy->SWvplot, 0.);
			if(strcmp(files[pLWin] , string_novalue) != 0) initialize_doublevector(egy->LWinplot, 0.);
			if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->LWgplot, 0.);
			if(strcmp(files[pLWv] , string_novalue) != 0) initialize_doublevector(egy->LWvplot, 0.);
			if(strcmp(files[pTs] , string_novalue) != 0) initialize_doublevector(egy->Tsplot, 0.);
			if(strcmp(files[pTg] , string_novalue) != 0) initialize_doublevector(egy->Tgplot, 0.);
			if(strcmp(files[pTv] , string_novalue) != 0) initialize_doublevector(egy->Tvplot, 0.);
			if(strcmp(files[pD] , string_novalue) != 0) initialize_doublevector(snow->Dplot, 0.);
			if(strcmp(files[pTa] , string_novalue) != 0) initialize_doublevector(met->Taplot, 0.);
			if(strcmp(files[pVspd] , string_novalue) != 0 || strcmp(files[pVdir] , string_novalue) != 0){
				initialize_doublevector(met->Vxplot, 0.);
				initialize_doublevector(met->Vyplot, 0.);
			}
			if(strcmp(files[pRH] , string_novalue) != 0) initialize_doublevector(met->RHplot, 0.);
			
			times->iplot ++;
		}
	}
	
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	//SAVING POINTS
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	
	if(times->time==0) isavings=par->recover;
	
	if(isavings < par->saving_points->nh){
		
		if(par->saving_points->nh==1 && par->saving_points->co[1]==0.0){
			
			isavings=1;
			
		}else{
						
			if(times->time+par->Dt >= par->saving_points->co[isavings+1]*86400.){
				isavings+=1;
				
				printf("Writing recovering files, saving point number %ld\n",isavings);
				
				flog = fopen(logfile, "a");
				fprintf(flog, "Writing recovering files, saving point number %ld\n",isavings);
				fclose(flog);
				
				write_suffix(NNNN, isavings, 0);
				
				for(l=0;l<=Nl;l++){
					if(strcmp(files[rpsi] , string_novalue) != 0) write_tensorseries_vector(1, l, isavings, files[rpsi], 0, par->format_out, sl->SS->P, UV, number_novalue, top->j_cont, Nr, Nc);
					if(l>0){
						if(strcmp(files[riceg] , string_novalue) != 0) write_tensorseries_vector(1, l, isavings, files[riceg], 0, par->format_out, sl->SS->thi, UV, number_novalue, top->j_cont, Nr, Nc);
						if(strcmp(files[rTg] , string_novalue) != 0) write_tensorseries_vector(1, l, isavings, files[rTg], 0, par->format_out, sl->SS->T, UV, number_novalue, top->j_cont, Nr, Nc);
					}
				}
				
				if(strcmp(files[rwcrn] , string_novalue) != 0){
					name = join_strings(files[rwcrn],NNNN);
					write_map_vector(name, 0, par->format_out, sl->VS->wrain, UV, number_novalue, top->j_cont, Nr, Nc);
					free(name);
				}
				
				if(strcmp(files[rwcsn] , string_novalue) != 0){
					name = join_strings(files[rwcsn],NNNN);
					write_map_vector(name, 0, par->format_out, sl->VS->wsnow, UV, number_novalue, top->j_cont, Nr, Nc);
					free(name);
				}
									
				for (i=1; i<=par->total_pixel; i++) {
					if ((long)sl->VS->Tv->co[i] == number_novalue) sl->VS->Tv->co[i] = 0.;
				}
				
				if(strcmp(files[rTv] , string_novalue) != 0){
					name = join_strings(files[rTv],NNNN);
					write_map_vector(name, 0, par->format_out, sl->VS->Tv, UV, number_novalue, top->j_cont, Nr, Nc);
					free(name);
				}	
									
				for(l=1;l<=par->max_snow_layers;l++){
					if(strcmp(files[rDzs] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rDzs], 0, par->format_out, snow->S->Dzl, UV, number_novalue);
					if(strcmp(files[rwls] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwls], 0, par->format_out, snow->S->w_liq, UV, number_novalue);
					if(strcmp(files[rwis] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwis], 0, par->format_out, snow->S->w_ice, UV, number_novalue);
					if(strcmp(files[rTs] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rTs], 0, par->format_out, snow->S->T, UV, number_novalue);
				}
				
				if(strcmp(files[rsnag] , string_novalue) != 0){
					name = join_strings(files[rsnag],NNNN);
					write_map_vector(name, 0, par->format_out, snow->age, UV, number_novalue, top->j_cont, Nr, Nc);
					free(name);
				}
				
				if(strcmp(files[rns] , string_novalue) != 0){
					M=copydouble_longmatrix(snow->S->lnum);
					name = join_strings(files[rns], NNNN);
					write_map(name, 1, par->format_out, M, UV, number_novalue);
					free_doublematrix(M);
					free(name);
				}
				
				if(par->max_glac_layers>0){
					for(l=1;l<=par->max_glac_layers;l++){
						if(strcmp(files[rDzi] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rDzi], 0, par->format_out, glac->G->Dzl, UV, number_novalue);
						if(strcmp(files[rwli] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwli], 0, par->format_out, glac->G->w_liq, UV, number_novalue);
						if(strcmp(files[rwii] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwii], 0, par->format_out, glac->G->w_ice, UV, number_novalue);
						if(strcmp(files[rTi] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rTi], 0, par->format_out, glac->G->T, UV, number_novalue);
					}
					
					if(strcmp(files[rni] , string_novalue) != 0){
						M=copydouble_longmatrix(glac->G->lnum);
						name = join_strings(files[rni], NNNN);
						write_map(name, 1, par->format_out, M, UV, number_novalue);
						free_doublematrix(M);
						free(name);
					}
				}		
				

				if(strcmp(files[rpsich] , string_novalue) != 0){
					M=new_doublematrix0_(Nl, par->total_pixel);
					for (l=0; l<=Nl; l++) {
						for (i=1; i<=par->total_pixel; i++) {
							M->co[l][i] = (double)number_novalue;
						}
						for (i=1; i<=par->total_channel; i++){
							r = cnet->r->co[i];
							c = cnet->c->co[i];
							M->co[l][top->j_cont[r][c]] = cnet->SS->P->co[l][i];
						}
						write_tensorseries_vector(1, l, isavings, files[rpsich], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
					}
					free_doublematrix(M);
				}

				if(strcmp(files[rTgch] , string_novalue) != 0){
					M=new_doublematrix(Nl, par->total_pixel);
					initialize_doublematrix(M, 0.0);
					for (l=1; l<=Nl; l++) {
						for (i=1; i<=par->total_pixel; i++) {
							M->co[l][i] = (double)number_novalue;
						}
						for (i=1; i<=par->total_channel; i++){
							r = cnet->r->co[i];
							c = cnet->c->co[i];
							M->co[l][top->j_cont[r][c]] = cnet->SS->T->co[l][i];
						}
						write_tensorseries_vector(1, l, isavings, files[rTgch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
					}
					free_doublematrix(M);
				}
				
				if(strcmp(files[ricegch] , string_novalue) != 0){
					M=new_doublematrix(Nl, par->total_pixel);
					initialize_doublematrix(M, 0.0);
					for (l=1; l<=Nl; l++) {
						for (i=1; i<=par->total_pixel; i++) {
							M->co[l][i] = (double)number_novalue;
						}
						for (i=1; i<=par->total_channel; i++){
							r = cnet->r->co[i];
							c = cnet->c->co[i];
							M->co[l][top->j_cont[r][c]] = cnet->SS->thi->co[l][i];
						}
						write_tensorseries_vector(1, l, isavings, files[ricegch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
					}
					free_doublematrix(M);
				}
				
			}
		}
		
	}
	
	//Continuous Recovery Option
	if (par->ContRecovery > 0) {
		t_rec += par->Dt;
		
		if (fabs(t_rec - par->ContRecovery*secinday)<1.E-5){
			t_rec = 0.;
			
			printf("Writing continuous-recovering files\n");
			
			flog = fopen(logfile, "a");
			fprintf(flog, "Writing continuous-recovering files\n");
			fclose(flog);
			
			if(strcmp(files[rtime] , string_novalue) != 0){
				name = join_strings(files[rtime], textfile);
				if(existing_file_woext(name) == 1){
					temp1 = join_strings(files[rtime], ".old");
					temp2 = join_strings(temp1, textfile);
					rename(name, temp2);
					free(temp1);
					free(temp2);
				}
				f = fopen(name, "w");
				fprintf(f,"Time[s],Time[d],n,i_run,i_sim,cum_time[s],elapsed_time[s],last_recover\n");
				fprintf(f,"%f,%f,%ld,%ld,%ld,%f,%f,%ld",times->time+par->Dt,(times->time+par->Dt)/secinday+(i_run-1)*(par->end_date->co[i_sim]-par->init_date->co[i_sim]),
						(long)(((times->time+par->Dt)/secinday)/par->ContRecovery),i_run,i_sim,cum_time,elapsed_time,par->n_ContRecovery);
				fclose(f);
				free(name);
			}
						
			write_suffix(NNNN, 0, 0);
			
			for(l=0;l<=Nl;l++){
				if(strcmp(files[rpsi] , string_novalue) != 0){
					rename_tensorseries(1, l, 0, files[rpsi]);
					write_tensorseries_vector(1, l, 0, files[rpsi], 0, par->format_out, sl->SS->P, UV, number_novalue, top->j_cont, Nr, Nc);
				}
				if(l>0){
					if(strcmp(files[riceg] , string_novalue) != 0){
						rename_tensorseries(1, l, 0, files[riceg]);
						write_tensorseries_vector(1, l, 0, files[riceg], 0, par->format_out, sl->SS->thi, UV, number_novalue, top->j_cont, Nr, Nc);
					}
					if(strcmp(files[rTg] , string_novalue) != 0){
						rename_tensorseries(1, l, 0, files[rTg]);
						write_tensorseries_vector(1, l, 0, files[rTg], 0, par->format_out, sl->SS->T, UV, number_novalue, top->j_cont, Nr, Nc);
					}
				}
			}
			
			if(strcmp(files[rwcrn] , string_novalue) != 0){
				name = join_strings(files[rwcrn],NNNN);
				rename_map(name);
				write_map_vector(name, 0, par->format_out, sl->VS->wrain, UV, number_novalue, top->j_cont, Nr, Nc);
				free(name);
			}
			
			if(strcmp(files[rwcsn] , string_novalue) != 0){
				name = join_strings(files[rwcsn],NNNN);
				rename_map(name);
				write_map_vector(name, 0, par->format_out, sl->VS->wsnow, UV, number_novalue, top->j_cont, Nr, Nc);
				free(name);
			}
			
			for (i=1; i<=par->total_pixel; i++) {
				if ((long)sl->VS->Tv->co[i] == number_novalue) sl->VS->Tv->co[i] = 0.;
			}
			
			if(strcmp(files[rTv] , string_novalue) != 0){
				name = join_strings(files[rTv],NNNN);				
				rename_map(name);
				write_map_vector(name, 0, par->format_out, sl->VS->Tv, UV, number_novalue, top->j_cont, Nr, Nc);
				free(name);
			}	
			
			for(l=1;l<=par->max_snow_layers;l++){
				if(strcmp(files[rDzs] , string_novalue) != 0){
					rename_tensorseries(1, l, 0, files[rDzs]);
					write_tensorseries(1, l, 0, files[rDzs], 0, par->format_out, snow->S->Dzl, UV, number_novalue);
				}
				if(strcmp(files[rwls] , string_novalue) != 0){
					rename_tensorseries(1, l, 0, files[rwls]);
					write_tensorseries(1, l, 0, files[rwls], 0, par->format_out, snow->S->w_liq, UV, number_novalue);
				}
				if(strcmp(files[rwis] , string_novalue) != 0){
					rename_tensorseries(1, l, 0, files[rwis]);
					write_tensorseries(1, l, 0, files[rwis], 0, par->format_out, snow->S->w_ice, UV, number_novalue);
				}
				if(strcmp(files[rTs] , string_novalue) != 0){
					rename_tensorseries(1, l, 0, files[rTs]);
					write_tensorseries(1, l, 0, files[rTs], 0, par->format_out, snow->S->T, UV, number_novalue);
				}
			}
			
			if(strcmp(files[rsnag] , string_novalue) != 0){
				name = join_strings(files[rsnag],NNNN);
				rename_map(name);
				write_map_vector(name, 0, par->format_out, snow->age, UV, number_novalue, top->j_cont, Nr, Nc);
				free(name);
			}
			
			if(strcmp(files[rns] , string_novalue) != 0){
				name = join_strings(files[rns], NNNN);
				rename_map(name);
				M=copydouble_longmatrix(snow->S->lnum);
				write_map(name, 1, par->format_out, M, UV, number_novalue);
				free_doublematrix(M);
				free(name);
			}
			
			if(par->max_glac_layers>0){
				for(l=1;l<=par->max_glac_layers;l++){
					if(strcmp(files[rDzi] , string_novalue) != 0){
						rename_tensorseries(1, l, 0, files[rDzi]);
						write_tensorseries(1, l, 0, files[rDzi], 0, par->format_out, glac->G->Dzl, UV, number_novalue);
					}
					if(strcmp(files[rwli] , string_novalue) != 0){
						rename_tensorseries(1, l, 0, files[rwli]);
						write_tensorseries(1, l, 0, files[rwli], 0, par->format_out, glac->G->w_liq, UV, number_novalue);
					}
					if(strcmp(files[rwii] , string_novalue) != 0){
						rename_tensorseries(1, l, 0, files[rwii]);
						write_tensorseries(1, l, 0, files[rwii], 0, par->format_out, glac->G->w_ice, UV, number_novalue);
					}
					if(strcmp(files[rTi] , string_novalue) != 0){
						rename_tensorseries(1, l, 0, files[rTi]);
						write_tensorseries(1, l, 0, files[rTi], 0, par->format_out, glac->G->T, UV, number_novalue);
					}
				}
				
				if(strcmp(files[rni] , string_novalue) != 0){
					name = join_strings(files[rni], NNNN);
					rename_map(name);				
					M=copydouble_longmatrix(glac->G->lnum);
					write_map(name, 1, par->format_out, M, UV, number_novalue);
					free_doublematrix(M);
					free(name);
				}
			}		
			
			if(strcmp(files[rpsich] , string_novalue) != 0){
				M=new_doublematrix0_(Nl, par->total_pixel);
				for (l=0; l<=Nl; l++) {
					rename_tensorseries(1, l, 0, files[rpsich]);
					for(i=1;i<=par->total_pixel;i++){
						M->co[l][i] = (double)number_novalue;
					}
					for (i=1; i<=par->total_channel; i++){
						r = cnet->r->co[i];
						c = cnet->c->co[i];
						M->co[l][top->j_cont[r][c]] = cnet->SS->P->co[l][i];
					}
					write_tensorseries_vector(1, l, 0, files[rpsich], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
				}
				free_doublematrix(M);
			}
			
			if(strcmp(files[rTgch] , string_novalue) != 0){
				M=new_doublematrix(Nl, par->total_pixel);
				for (l=1; l<=Nl; l++) {
					rename_tensorseries(1, l, 0, files[rTgch]);
					for(i=1;i<=par->total_pixel;i++){
						M->co[l][i] = (double)number_novalue;
					}
					for (i=1; i<=par->total_channel; i++){
						r = cnet->r->co[i];
						c = cnet->c->co[i];
						M->co[l][top->j_cont[r][c]] = cnet->SS->T->co[l][i];
					}
					write_tensorseries_vector(1, l, 0, files[rTgch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
				}
				free_doublematrix(M);
			}
			
			if(strcmp(files[ricegch] , string_novalue) != 0){
				M=new_doublematrix(Nl, par->total_pixel);
				for (l=1; l<=Nl; l++) {
					rename_tensorseries(1, l, 0, files[ricegch]);
					for(i=1;i<=par->total_pixel;i++){
						M->co[l][i] = (double)number_novalue;
					}
					for (i=1; i<=par->total_channel; i++){
						r = cnet->r->co[i];
						c = cnet->c->co[i];
						M->co[l][top->j_cont[r][c]] = cnet->SS->thi->co[l][i];
					}
					write_tensorseries_vector(1, l, 0, files[ricegch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
				}
				free_doublematrix(M);
			}
			
			if(par->Tzrun == 1 && strcmp(files[rTrun] , string_novalue) != 0) print_run_averages_for_recover(sl->Tzrun, files[rTrun], top->j_cont, par, Nl, Nr, Nc);
			if(par->wzrun == 1 && strcmp(files[rTrun] , string_novalue) != 0) print_run_averages_for_recover(sl->wzrun, files[rwrun], top->j_cont, par, Nl, Nr, Nc);
			if(par->Tzmaxrun == 1 && strcmp(files[rTmaxrun] , string_novalue) != 0) print_run_averages_for_recover(sl->Tzmaxrun, files[rTmaxrun], top->j_cont, par, Nl, Nr, Nc);
			if(par->wzmaxrun == 1 && strcmp(files[rwmaxrun] , string_novalue) != 0) print_run_averages_for_recover(sl->wzmaxrun, files[rwmaxrun], top->j_cont, par, Nl, Nr, Nc);
			if(par->Tzminrun == 1 && strcmp(files[rTminrun] , string_novalue) != 0) print_run_averages_for_recover(sl->Tzminrun, files[rTminrun], top->j_cont, par, Nl, Nr, Nc);
			if(par->wzminrun == 1 && strcmp(files[rwminrun] , string_novalue) != 0) print_run_averages_for_recover(sl->wzminrun, files[rwminrun], top->j_cont, par, Nl, Nr, Nc);
			if(par->dUzrun == 1 && strcmp(files[rdUrun] , string_novalue) != 0) print_run_averages_for_recover(sl->dUzrun, files[rdUrun], top->j_cont, par, Nl, Nr, Nc);
			if(par->SWErun == 1 && strcmp(files[rSWErun] , string_novalue) != 0) print_run_averages_for_recover(sl->SWErun, files[rSWErun], top->j_cont, par, 3, Nr, Nc);

			if(strcmp(files[rsux] , string_novalue) != 0){
				if (existing_file_woext(files[rsux]) == 1){
					temp1 = join_strings(files[rsux], ".old");
					rename(files[rsux], temp1);
					free(temp1);
				}
				f = fopen(files[rsux], "w");
				fclose(f);
			}

		}
	}
		
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac){

	/*internal auxiliary variables:*/
	long i,l,m,j,r,c;
	char *name,*temp,*temp2,NNNN[ ]={"NNNN"},rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	long sy;
	short lu, first_column;
	DOUBLEVECTOR *root_fraction;
	FILE *f; 
		
	if (par->recover > 0) write_suffix(rec, par->recover, 4);
	if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);
			
	//DISCHARGE
	if (par->state_discharge == 1 && strcmp(files[fQ] , string_novalue) != 0){
		if (par->recover>0) {
			temp = join_strings(files[fQ], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fQ], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fQ], textfile);
		}
		
		f=t_fopen(name,"w");
		fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,Qtot[m3/s],Vsup/Dt[m3/s],Vsub/Dt[m3/s],Vchannel[m3],Qoutlandsup[m3/s],Qoutlandsub[m3/s],Qoutbottom[m3/s]\n");
		t_fclose(f);
		free(name);
	}
	
	
	if(par->state_pixel == 1){
		
		//output matrix and vectors
		m=(long)otot;
		odpnt=(double**)malloc(m*sizeof(double*));
		odp=(double**)malloc(m*sizeof(double*));
		for (i=0; i<otot; i++) {
			odpnt[i]=(double*)malloc(par->rc->nrh*sizeof(double));
			odp[i]=(double*)malloc(par->rc->nrh*sizeof(double));
			for (j=0; j<par->rc->nrh; j++) {
				odpnt[i][j] = 0.;
				odp[i][j] = 0.;
			}
		}
						
		if(strcmp(files[fpointwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
				temp = join_strings(files[fpointwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fpointwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fpointwriteend], textfile);
			}
			
			ffpoint=fopen(name,"w");
			first_column=1;
			for(j=0;j<nopnt;j++){
				if(first_column==0){
					fprintf(ffpoint,",");
				}else {
					first_column = 0;
				}					
				if (opnt[j] >= 0) {
					fprintf(ffpoint,"%s",hpnt[opnt[j]]);					
				}else {
					fprintf(ffpoint, "None");
				}
			}
			fprintf(ffpoint,"\n");
			free(name);
		}
		
		if(par->max_glac_layers>0){
			if(strcmp(files[fglzwriteend] , string_novalue) != 0){
				
				if (par->recover>0) {
					temp = join_strings(files[fpointwriteend], rec);
					name = join_strings(temp, textfile);
					free(temp);
				}else if (par->n_ContRecovery>0) {
					temp = join_strings(files[fpointwriteend], crec);
					name = join_strings(temp, textfile);
					free(temp);
				}else {
					name = join_strings(files[fpointwriteend], textfile);
				}
				
				ffglac=t_fopen(name,"w");

				if ((long)par->glac_plot_depths->co[1] != number_novalue) {
					m = par->glac_plot_depths->nh;
				}else {
					m = par->max_glac_layers;
				}
				first_column=1;
				for(j=0;j<noglc;j++){
					if(first_column==0){
						fprintf(ffglac,",");
					}else {
						first_column = 0;
					}			
					if (oglc[j] >= 0 && oglc[j]<=5) {
						fprintf(ffglac,"%s",hglc[oglc[j]]);					
					}else if (oglc[j] >= 6 && oglc[j] < 6+3*m) {
						l = (long)fmod( (double)oglc[j]-6., (double)m ) + 1;
						n = floor( ( (double)oglc[j]-6.) / (double)m ) + 6;
						if ((long)par->glac_plot_depths->co[1] != number_novalue) {
							fprintf(ffglac, "%s(%f)",hglc[n],par->glac_plot_depths->co[l]);
						}else {
							fprintf(ffglac, "%s(%ld)",hglc[n],l);
						}
					}else if (oglc[j] >= 6+3*m) {
						l = (long)fmod( (double)oglc[j]-6.-3*(double)m, (double)par->max_glac_layers ) + 1;
						n = floor( ( (double)oglc[j]-6.-3.*(double)m) / (double)par->max_glac_layers ) + 6 + 3;
						fprintf(ffglac, "%s(%ld)",hglc[n],l);
					}else{
						fprintf(ffglac, "None");
					}
				}
				fprintf(ffglac,"\n");
				free(name);
			}
		}
		
		if(strcmp(files[fTzwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
				temp = join_strings(files[fTzwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fTzwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fTzwriteend], textfile);
			}
			
			ffT=fopen(name,"w");
			write_soil_header(ffT, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fTzavwriteend] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[fTzavwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fTzavwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fTzavwriteend], textfile);
			}
			
			ffTav=fopen(name,"w");		
			write_soil_header(ffTav, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fpsiztotwriteend] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[fpsiztotwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fpsiztotwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fpsiztotwriteend], textfile);
			}
			
			ffpsitot=fopen(name,"w");		
			write_soil_header(ffpsitot, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fpsizwriteend] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[fpsizwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fpsizwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fpsizwriteend], textfile);
			}

			ffpsi=fopen(name,"w");		
			write_soil_header(ffpsi, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fliqzwriteend] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[fliqzwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fliqzwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fliqzwriteend], textfile);
			}
			
			ffliq=fopen(name,"w");		
			write_soil_header(ffliq, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}

		if(strcmp(files[fliqzavwriteend] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[fliqzavwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fliqzavwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fliqzavwriteend], textfile);
			}

			ffliqav=fopen(name,"w");		
			write_soil_header(ffliqav, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}

		if(strcmp(files[ficezwriteend] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[ficezwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[ficezwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[ficezwriteend], textfile);
			}

			ffice=fopen(name,"w");				
			write_soil_header(ffice, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[ficezavwriteend] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[ficezavwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[ficezavwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[ficezavwriteend], textfile);
			}
			
			fficeav=fopen(name,"w");		
			write_soil_header(fficeav, par->soil_plot_depths, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fsnTzwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
				temp = join_strings(files[fsnTzwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fsnTzwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fsnTzwriteend], textfile);
			}
			
			ffsnowT=fopen(name,"w");	
			write_snow_header(0, 1, 1, ffsnowT, par->snow_plot_depths, snow->S->Dzl);
			free(name);
		}

		if(strcmp(files[fsnlzwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
				temp = join_strings(files[fsnlzwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fsnlzwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fsnlzwriteend], textfile);
			}
			
			ffsnowl=fopen(name,"w");	
			write_snow_header(0, 1, 1, ffsnowl, par->snow_plot_depths, snow->S->Dzl);
			free(name);
		}

		if(strcmp(files[fsnizwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
				temp = join_strings(files[fsnizwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fsnizwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fsnizwriteend], textfile);
			}
			
			ffsnowT=fopen(name,"w");	
			write_snow_header(0, 1, 1, ffsnowi, par->snow_plot_depths, snow->S->Dzl);
			free(name);
		}

		if(strcmp(files[fsndzwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
				temp = join_strings(files[fsndzwriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fsndzwriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fsndzwriteend], textfile);
			}
			
			ffsnowT=fopen(name,"w");	
			write_snow_header(2, 1, 1, ffsnowd, par->snow_plot_depths, snow->S->Dzl);
			free(name);
		}
		
		
		
		root_fraction=new_doublevector(Nl);
   
		
		//DATA POINTS
		for(i=1;i<=par->rc->nrh;i++){
						
			write_suffix(NNNN, par->IDpoint->co[i], 0);
			r=par->rc->co[i][1];
			c=par->rc->co[i][2];				
			sy=sl->type->co[r][c];
			lu=(short)land->LC->co[r][c];
						
			if(strcmp(files[fpoint] , string_novalue) != 0 && par->point_sim != 1){			
				name=join_strings(files[fpoint],"_info_");
				temp=join_strings(name,NNNN);
				temp2=join_strings(temp,textfile);
				f=t_fopen(temp2,"w");
				
				fprintf(f," The main properties of the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld are:\n",par->chkpt->co[i][ptX],par->chkpt->co[i][ptY],r,c);
				fprintf(f," Elevation above sea level: %10.3f m\n",top->Z0->co[r][c]);
				fprintf(f," Gauckler-Strickler [m^1/3/s]: %f\n",land->ty->co[lu][jcm]);
				for(l=1;l<=Nl;l++){
					fprintf(f," Residual water content[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jres][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Saturated water content[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jsat][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Alpha of van Genuchten[mm^-1] of the layer %ld: %f\n",l,sl->pa->co[sy][ja][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," n of van Genuchten[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jns][l]);
				}	
				for(l=1;l<=Nl;l++){
					fprintf(f," m of van Genuchten[-] of the layer %ld: %f\n",l,1-1/sl->pa->co[sy][jns][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," v of van Genuchten[-] of the layer %ld: %f\n",l,sl->pa->co[sy][jv][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Water content of wilting point [-] of the layer %ld: %f\n",l,sl->pa->co[sy][jwp][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Water content of field capacity [-] of the layer %ld: %f\n",l,sl->pa->co[sy][jfc][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Kv_sat of layer %ld [mm/s]: %f\n",l,sl->pa->co[sy][jKn][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Kh_sat of layer %ld [mm/s]: %f\n",l,sl->pa->co[sy][jKl][l]);
				}

				fprintf(f," Terrain elevation [m]: %f\n",top->Z0->co[r][c]);
				fprintf(f," Sky view factor [-]: %f\n",top->sky->co[r][c]);
				fprintf(f," The pixel-type is %d \n",top->pixel_type->co[r][c]);
				fprintf(f," Aspect [deg] [0=Nord, clockwise]: %f \n",top->aspect->co[r][c]);
				fprintf(f," Mean slope of the pixel [deg]: %f \n",top->slope->co[r][c]);
				fprintf(f," Land use number is %d \n",(short)land->LC->co[r][c]);

				int lmax = land->root_fraction->nch - land->root_fraction->ncl + 1;

				for(l=1;l<=lmax;l++){
					fprintf(f," The root fraction [-] of layer %ld: %f\n",l,land->root_fraction->co[lu][l]);
				}
		
				fprintf(f," Surface fraction of land covered by vegetation [-]: %f \n",land->ty->co[lu][jcf]);
				fprintf(f," Leaf and Stem Area Index [-]: %f \n",land->ty->co[lu][jLSAI]);
				fprintf(f," Momentum roughness length z0soil [m]: %f \n",land->ty->co[lu][jz0]);
				fprintf(f," Vegetation height [m]: %f \n",land->ty->co[lu][jHveg]);

				fprintf(f," \n");
				t_fclose(f);
				free(temp2);
				free(temp);
				free(name);
			}
			
			if(strcmp(files[fpoint] , string_novalue) != 0){
				
				temp=join_strings(files[fpoint],NNNN);

				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");
				first_column=1;
				for(j=0;j<nopnt;j++){
					if(first_column==0){
						fprintf(f,",");
					}else {
						first_column = 0;
					}					
					if (opnt[j] >= 0) {
						fprintf(f,"%s",hpnt[opnt[j]]);					
					}else {
						fprintf(f, "None");
					}
				}
				fprintf(f,"\n");
				t_fclose(f);
				free(name);
			}
			
			if(par->max_glac_layers>0){
				if(strcmp(files[fglz] , string_novalue) != 0){
					
					temp=join_strings(files[fglz],NNNN);
					
					if (par->recover>0) {
						temp2 = join_strings(temp, rec);
						name = join_strings(temp2, textfile);
						free(temp2);
					}else if (par->n_ContRecovery>0) {
						temp2 = join_strings(temp, crec);
						name = join_strings(temp2, textfile);
						free(temp2);
					}else {
						name = join_strings(temp, textfile);
					}
					
					free(temp);
					
					if ((long)par->glac_plot_depths->co[1] != number_novalue) {
						m = par->glac_plot_depths->nh;
					}else {
						m = par->max_glac_layers;
					}

					f=t_fopen(name,"w");
					first_column=1;
					for(j=0;j<noglc;j++){
						if(first_column==0){
							fprintf(f,",");
						}else {
							first_column = 0;
						}			
						if (oglc[j] >= 0 && oglc[j]<=5) {
							fprintf(f,"%s",hglc[oglc[j]]);					
						}else if (oglc[j] >= 6 && oglc[j] < 6+3*m) {
							l = (long)fmod( (double)oglc[j]-6., (double)m ) + 1;
							n = floor( ( (double)oglc[j]-6.) / (double)m ) + 6;
							if ((long)par->glac_plot_depths->co[1] != number_novalue) {
								fprintf(f, "%s(%f)",hglc[n],par->glac_plot_depths->co[l]);
							}else {
								fprintf(f, "%s(%ld)",hglc[n],l);
							}
						}else if (oglc[j] >= 6+3*m) {
							l = (long)fmod( (double)oglc[j]-6.-3*(double)m, (double)par->max_glac_layers ) + 1;
							n = floor( ( (double)oglc[j]-6.-3.*(double)m) / (double)par->max_glac_layers ) + 6 + 3;
							fprintf(f, "%s(%ld)",hglc[n],l);
						}else{
							fprintf(f, "None");
						}
					}
					fprintf(f,"\n");
					t_fclose(f);
					free(name);
				}
			}
				
			if(strcmp(files[fTz] , string_novalue) != 0){
				
				temp=join_strings(files[fTz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");			
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
			}
			
			if(strcmp(files[fTzav] , string_novalue) != 0){

				temp=join_strings(files[fTzav],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");		
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
			}
		
			if(strcmp(files[fpsiztot] , string_novalue) != 0){

				temp=join_strings(files[fpsiztot],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);

				f=t_fopen(name,"w");		
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);
				free(name);
			}
			
			if(strcmp(files[fpsiz] , string_novalue) != 0){

				temp=join_strings(files[fpsiz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);

				f=t_fopen(name,"w");		
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);
				free(name);
			}
		
			if(strcmp(files[fliqz] , string_novalue) != 0){

				temp=join_strings(files[fliqz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);

				f=t_fopen(name,"w");		
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);
				free(name);
			}
			
			if(strcmp(files[fliqzav] , string_novalue) != 0){

				temp=join_strings(files[fliqzav],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);

				f=t_fopen(name,"w");		
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
			}
		
			if(strcmp(files[ficez] , string_novalue) != 0){

				temp=join_strings(files[ficez],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);

				f=t_fopen(name,"w");				
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
			}
							
			if(strcmp(files[ficezav] , string_novalue) != 0){

				temp=join_strings(files[ficezav],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);

				f=t_fopen(name,"w");		
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
			}

			if(strcmp(files[fsatz] , string_novalue) != 0){
				
				temp=join_strings(files[fsatz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");				
				write_soil_header(f, par->soil_plot_depths, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
			}
			
			if(strcmp(files[fsnTz] , string_novalue) != 0){
				
				temp=join_strings(files[fsnTz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");
				write_snow_header(0, 1, 1, f, par->snow_plot_depths, snow->S->Dzl);
				t_fclose(f);	
				free(name);
			}
			
			if(strcmp(files[fsnlz] , string_novalue) != 0){
				
				temp=join_strings(files[fsnlz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");
				write_snow_header(0, 1, 1, f, par->snow_plot_depths, snow->S->Dzl);
				t_fclose(f);	
				free(name);
			}

			if(strcmp(files[fsniz] , string_novalue) != 0){
				
				temp=join_strings(files[fsniz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");
				write_snow_header(0, 1, 1, f, par->snow_plot_depths, snow->S->Dzl);
				t_fclose(f);	
				free(name);
			}

			if(strcmp(files[fsndz] , string_novalue) != 0){
				
				temp=join_strings(files[fsndz],NNNN);
				
				if (par->recover>0) {
					temp2 = join_strings(temp, rec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else if (par->n_ContRecovery>0) {
					temp2 = join_strings(temp, crec);
					name = join_strings(temp2, textfile);
					free(temp2);
				}else {
					name = join_strings(temp, textfile);
				}
				
				free(temp);
				
				f=t_fopen(name,"w");
				write_snow_header(2, 1, 1, f, par->snow_plot_depths, snow->S->Dzl);
				t_fclose(f);	
				free(name);
			}
						
		}
		
		free_doublevector(root_fraction);
		
	}
	
	m=(long)ootot;
	odbsn=(double*)malloc(m*sizeof(double));
	odb=(double*)malloc(m*sizeof(double));
	for (i=0; i<ootot; i++) {
		odbsn[i]=0.;
		odb[i]=0.;
	}		

	if(par->state_basin == 1){
				
		//DATA BASIN
		if(strcmp(files[fbaswriteend] , string_novalue) != 0){
						
			if (par->recover>0) {
				temp = join_strings(files[fbaswriteend], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fbaswriteend], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fbaswriteend], textfile);
			}
			
			ffbas=fopen(name,"w");
			
			first_column=1;
			for(j=0;j<nobsn;j++){
				if(first_column==0){
					fprintf(ffbas,",");
				}else {
					first_column = 0;
				}					
				if (obsn[j] >= 0) {
					fprintf(ffbas,"%s",hbsn[obsn[j]]);					
				}else {
					fprintf(ffbas, "None");
				}
			}
			fprintf(ffbas,"\n");
			free(name);
		}
		
		if(strcmp(files[fbas] , string_novalue) != 0){

			if (par->recover>0) {
				temp = join_strings(files[fbas], rec);
				name = join_strings(temp, textfile);
				free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fbas], crec);
				name = join_strings(temp, textfile);
				free(temp);
			}else {
				name = join_strings(files[fbas], textfile);
			}

			f=t_fopen(name,"w");
			
			first_column=1;
			for(j=0;j<nobsn;j++){
				if(first_column==0){
					fprintf(f,",");
				}else {
					first_column = 0;
				}					
				if (obsn[j] >= 0) {
					fprintf(f,"%s",hbsn[obsn[j]]);					
				}else {
					fprintf(f, "None");
				}
			}
			fprintf(f,"\n");
			
			t_fclose(f);
			free(name);
		}
				
	}
	
	//SNOW COVERED AREA STATISTICS
	if(par->point_sim!=1 && strcmp(files[fSCA] , string_novalue) != 0){
		
		if (par->recover>0) {
			temp = join_strings(files[fSCA], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fSCA], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fSCA], textfile);
		}
		
		f=t_fopen(name,"w");
		fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,snowDav,SWEav,Tav,Tsav,perc.SFA,perc.SCA\n");
		t_fclose(f);
		free(name);
	}
	
	

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_output(long i, long iname, double init_date, double end_date, double JDfrom0, double JD, long day, long month, long year, long hour, long minute, DOUBLEVECTOR *n, SOIL *sl, PAR *par, double psimin, double cosslope){

	char *name,*temp,*temp2,NNNN[ ]={"NNNN"};
	char rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	long l;
	FILE *f;
	
	write_suffix(NNNN, iname, 0);	
	
	if (par->recover > 0) write_suffix(rec, par->recover, 4);
	if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

	if(strcmp(files[fTz] , string_novalue) != 0){
		
		temp=join_strings(files[fTz],NNNN);

		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}
		
		f=fopen(name,"a");	
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Tzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fTzwriteend] , string_novalue) != 0){
		write_soil_file(1, iname, ffT, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Tzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}
	
	if(strcmp(files[fTzav] , string_novalue) != 0){
		temp=join_strings(files[fTzav],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");	
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Tzavplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fTzavwriteend] , string_novalue) != 0){
		write_soil_file(1, iname, ffTav, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Tzavplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}
	
	if(strcmp(files[fpsiztot] , string_novalue) != 0){
		temp=join_strings(files[fpsiztot],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");		
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Ptotzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fpsiztotwriteend] , string_novalue) != 0){
		write_soil_file(1, iname, ffpsitot, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Ptotzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}
	
	if(strcmp(files[fpsiz] , string_novalue) != 0){
		temp=join_strings(files[fpsiz],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");		
		write_soil_file(0, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Pzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fpsizwriteend] , string_novalue) != 0){
		write_soil_file(0, iname, ffpsi, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->Pzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}
	
	if(strcmp(files[fliqz] , string_novalue) != 0){
		temp=join_strings(files[fliqz],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");		
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fliqzwriteend] , string_novalue) != 0){
		write_soil_file(1, iname, ffliq, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thzplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}

	if(strcmp(files[fliqzav] , string_novalue) != 0){
		temp=join_strings(files[fliqzav],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");		
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thzavplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fliqzavwriteend] , string_novalue) != 0){
		write_soil_file(1, iname, ffliqav, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thzavplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}
	
	if(strcmp(files[ficez] , string_novalue) != 0){
		temp=join_strings(files[ficez],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");		
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thizplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[ficezwriteend] , string_novalue) != 0){
		write_soil_file(1, iname, ffice, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thizplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}

	if(strcmp(files[ficezav] , string_novalue) != 0){
		temp=join_strings(files[ficezav],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");		
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thizavplot->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[ficezavwriteend] , string_novalue) != 0){
		write_soil_file(1, iname, fficeav, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->thizavplot->co[i], n, sl->pa->co[1][jdz], cosslope);
	}
	
	if(strcmp(files[fsatz] , string_novalue) != 0){
		temp=join_strings(files[fsatz],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}		
		
		f=fopen(name,"a");		
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, end_date, sl->satratio->co[i], n, sl->pa->co[1][jdz], cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	for(l=1;l<=Nl;l++){
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] = 0.0;
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] = 0.0;
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot->co[i][l] = 0.0;
	}
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_snow_output(long i, long iname, long r, long c, double init_date, double end_date, double JDfrom0, double JD, 
					   long day, long month, long year, long hour, long minute, DOUBLEVECTOR *n, STATEVAR_3D *snow, PAR *par, double cosslope){
	
	char *name,*temp,*temp2,NNNN[ ]={"NNNN"};
	char rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	FILE *f;
	
	write_suffix(NNNN, iname, 0);	
	
	if (par->recover > 0) write_suffix(rec, par->recover, 4);
	if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);
	
	if(strcmp(files[fsnTz] , string_novalue) != 0){
		
		temp=join_strings(files[fsnTz],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}
		
		f=fopen(name,"a");	
		write_snow_file(0, iname, r, c, snow->lnum->co[r][c], f, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->T, cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fsnTzwriteend] , string_novalue) != 0){
		write_snow_file(0, iname, r, c, snow->lnum->co[r][c], ffsnowT, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->T, cosslope);
	}
	
	if(strcmp(files[fsnlz] , string_novalue) != 0){
		
		temp=join_strings(files[fsnlz],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}
		
		f=fopen(name,"a");	
		write_snow_file(1, iname, r, c, snow->lnum->co[r][c], f, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->w_liq, cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fsnlzwriteend] , string_novalue) != 0){
		write_snow_file(1, iname, r, c, snow->lnum->co[r][c], ffsnowl, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->w_liq, cosslope);
	}
	
	if(strcmp(files[fsniz] , string_novalue) != 0){
		
		temp=join_strings(files[fsniz],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}
		
		f=fopen(name,"a");	
		write_snow_file(1, iname, r, c, snow->lnum->co[r][c], f, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->w_ice, cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fsnizwriteend] , string_novalue) != 0){
		write_snow_file(1, iname, r, c, snow->lnum->co[r][c], ffsnowi, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->w_ice, cosslope);
	}
	
	if(strcmp(files[fsndz] , string_novalue) != 0){
		
		temp=join_strings(files[fsndz],NNNN);
		
		if (par->recover>0) {
			temp2 = join_strings(temp, rec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else if (par->n_ContRecovery>0) {
			temp2 = join_strings(temp, crec);
			name = join_strings(temp2, textfile);
			free(temp2);
		}else {
			name = join_strings(temp, textfile);
		}
		
		f=fopen(name,"a");	
		write_snow_file(2, iname, r, c, snow->lnum->co[r][c], f, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->Dzl, cosslope);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fsndzwriteend] , string_novalue) != 0){
		write_snow_file(2, iname, r, c, snow->lnum->co[r][c], ffsnowd, day, month, year, hour, minute, JDfrom0, init_date, end_date, n, snow->Dzl, snow->Dzl, cosslope);
	}
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_file(long lmin, long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, 
					 double JDfrom0end, double *var, DOUBLEVECTOR *n, double *dz, double cosslope){

	short first_column=1;
	long j, l;
	
	for (j=0; j<nosl; j++) {
		if(first_column==0){
			fprintf(f,",");
		}else {
			first_column = 0;
		}
		if (osl[j] >= 0) {
			if (osl[j] == 0) {
				fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)d,(float)m,(float)y,(float)h,(float)mi);		
			}else if (osl[j] == 1) {
				fprintf(f, "%f",JDfrom0);
			}else if (osl[j] == 2) {
				fprintf(f, "%f",(JDfrom0-JDfrom0init)+(i_run-1)*(JDfrom0end-JDfrom0init));
			}else if (osl[j] == 3) {
				fprintf(f, "%ld",i_sim);
			}else if (osl[j] == 4) {
				fprintf(f, "%ld",i_run);																				
			}else if (osl[j] == 5) {
				fprintf(f, "%ld",i);
			}
		}else {
			fprintf(f,"%f",(double)number_novalue);
		}
	}

	if ((long)n->co[1] != number_novalue) {
		for (l=1; l<=n->nh; l++) {
			fprintf(f, ",%f",interpolate_soil(lmin, n->co[l]*cosslope, Nl, dz, var));
		}
	}else{
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f",var[l]);
		}
	}
	
	fprintf(f," \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//a=0 snow(according to snow_depth_plot) and var used
//a=1 snow(according to snow_depth_plot) and var/dz used
//a=2 snow(according to layers) and var used
//a=3 snow(according to layers) and var/dz used

void write_snow_file(short a, long i, long r, long c, long lmax, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, 
					 double JDfrom0init, double JDfrom0end, DOUBLEVECTOR *n, DOUBLETENSOR *snowDz, DOUBLETENSOR *var, double cosslope){
	
	short first_column=1;
	long j, l;
	
	for (j=0; j<nosnw; j++) {
		if(first_column==0){
			fprintf(f,",");
		}else {
			first_column = 0;
		}
		if (osnw[j] >= 0) {
			if (osnw[j] == 0) {
				fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)d,(float)m,(float)y,(float)h,(float)mi);		
			}else if (osnw[j] == 1) {
				fprintf(f, "%f",JDfrom0);
			}else if (osnw[j] == 2) {
				fprintf(f, "%f",(JDfrom0-JDfrom0init)+(i_run-1)*(JDfrom0end-JDfrom0init));
			}else if (osnw[j] == 3) {
				fprintf(f, "%ld",i_sim);
			}else if (osnw[j] == 4) {
				fprintf(f, "%ld",i_run);																				
			}else if (osnw[j] == 5) {
				fprintf(f, "%ld",i);
			}
		}else {
			fprintf(f,"%f",(double)number_novalue);
		}
	}
	
	if ( (a==0 || a==1) && (long)n->co[1] != number_novalue ) {
		for (l=1; l<=n->nh; l++) {
			fprintf(f, ",%f",interpolate_snow(r, c, n->co[l]*cosslope, lmax, snowDz, var, a));
		}
	}else{
		if (a==2 || a==0) {
			for(l=1;l<=snowDz->ndh;l++){
				if (snowDz->co[l][r][c] > 0) {
					fprintf(f,",%f",var->co[l][r][c]);
				}else {
					fprintf(f,",%f",(double)number_novalue);
				}
			}
		}else if (a==3 || a==1){
			for(l=1;l<=snowDz->ndh;l++){
				if (snowDz->co[l][r][c] > 0) {
					fprintf(f,",%f",var->co[l][r][c]/snowDz->co[l][r][c]);
				}else {
					fprintf(f,",%f",(double)number_novalue);
				}
			}
		}
	}
	
	fprintf(f," \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_header(FILE *f, DOUBLEVECTOR *n, double *dz){
	
	short first_column=1;
	long j, l;
	double z=0.0;
	
	for(j=0;j<nosl;j++){
		if(first_column==0){
			fprintf(f,",");
		}else {
			first_column = 0;
		}					
		if (osl[j] >= 0 && osl[j]<=5) {
			fprintf(f,"%s",hsl[osl[j]]);					
		}else{
			fprintf(f, "none");
		}
	}

	if ((long)n->co[1] != number_novalue ) {
		for (l=1; l<=n->nh; l++) {
			fprintf(f, ",%f",n->co[l]);
		}
	}else{
		for(l=1;l<=Nl;l++){
			z += dz[l];
			fprintf(f,",%f ",z-0.5*dz[l]);
		}
	}
	
	fprintf(f," \n");
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************


//a=0 snow(according to snow_depth_plot) and var used
//a=1 snow(according to snow_depth_plot) and var/dz used
//a=2 snow(according to layers) and var used
//a=3 snow(according to layers) and var/dz used

void write_snow_header(short a, long r, long c, FILE *f, DOUBLEVECTOR *n, DOUBLETENSOR *Dz){
	
	short first_column=1;
	long j, l;
	
	for(j=0;j<nosnw;j++){
		if(first_column==0){
			fprintf(f,",");
		}else {
			first_column = 0;
		}				
		if (osnw[j] >= 0&& osnw[j]<=5) {
			fprintf(f,"%s",hsnw[osnw[j]]);					
		}else{
			fprintf(f, "none");
		}
	}
	
	if ((a==0 || a==1) && (long)n->co[1] != number_novalue ) {
		for (l=1; l<=n->nh; l++) {
			fprintf(f, ",%f",n->co[l]);
		}
	}else{
		for(l=1;l<=Dz->ndh;l++){
			fprintf(f,",L%ld ",l);
		}
	}
	
	fprintf(f," \n");
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void plot(char *name, long i_plot, DOUBLEVECTOR *V, short format, long **J){
	
	char ADS[ ]={"iiii"};
	char *temp;
	
	write_suffix(ADS, i_plot, 0);
	temp=join_strings(name,ADS);
	write_map_vector(temp, 0, format, V, UV, number_novalue, J, Nr, Nc);
	free(temp);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

double interpolate_soil(long lmin, double h, long max, double *Dz, double *Q){
	
	double q, z, z0=0.;
	long l;
		
	l = lmin;	
	q = (double)number_novalue;
	
	do{
		
		if (l == lmin){
			z = z0;
			if (l>0) z += Dz[l]/2.;
		}else if (l <= max) {
			z = z0 + Dz[l]/2;
			if (l>1) z += Dz[l-1]/2.;
		}else {
			z = z0 + Dz[max]/2.;
		}
		
		if(h < z && h >= z0){
			if (l == lmin) {
				q = Q[lmin];
			}else if (l <= max) {
				q = ( Q[l-1] * (z-h) + Q[l] * (h-z0) ) / (z - z0);
			}else {
				q = Q[max];
			}
		}
		
		z0 = z;
		
		l ++;
		
	}while ( (long)q == number_novalue && l <= max+1 );
	
	return q;
	
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

double interpolate_soil2(long lmin, double h, long max, double *Dz, DOUBLEMATRIX *Q, long i){
	
	double q, z, z0=0.;
	long l;
	
	l = lmin;
	q = (double)number_novalue;
	
	do{
		
		if (l == lmin){
			z = z0;
			if (l>0) z += Dz[l]/2.;
		}else if (l <= max) {
			z = z0 + Dz[l]/2.;
			if (l>1) z += Dz[l-1]/2.;
		}else {
			z = z0 + Dz[max]/2.;
		}
		
		if(h < z && h >= z0){
			if (l == lmin) {
				q = Q->co[l][i];
			}else if (l <= max) {
				q = ( Q->co[l-1][i] * (z-h) + Q->co[l][i] * (h-z0) ) / (z - z0);
			}else {
				q = Q->co[max][i];
			}
		}
		
		z0 = z;
		
		l ++;
		
	}while ( (long)q == number_novalue && l <= max+1 );
	
	return q;
	
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_tensorseries_soil(long lmin, char *suf, char *filename, short type, short format, DOUBLEMATRIX *T, DOUBLEVECTOR *n, long **J, LONGMATRIX *RC, double *dz, DOUBLEMATRIX *slope, short vertical){
	
	char LLLLL[ ]={"LLLLL"};
	char *temp1, *temp2;
	long i, l, npoints=T->nch;
	double cosslope=1.;
	DOUBLEVECTOR *V;
	
	for(l=1; l<=n->nh; l++){
		temp1 = join_strings(LLLLL, suf);
		write_suffix(temp1, l, 1);	
	
		V = new_doublevector(npoints);
	
		for(i=1; i<=npoints; i++){
			if (vertical == 1) cosslope = cos( Fmin(max_slope, slope->co[RC->co[i][1]][RC->co[i][2]]) * Pi/180. );
			V->co[i] = interpolate_soil2(lmin, n->co[l]*cosslope, Nl, dz, T, i);
		}
	
		temp2 = join_strings(filename, temp1);
		write_map_vector(temp2, type, format, V, UV, number_novalue, J, slope->nrh, slope->nch);
	
		free_doublevector(V);
		free(temp1);
		free(temp2);
	}
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void fill_output_vectors(double Dt, double W, ENERGY *egy, SNOW *snow, GLACIER *glac, WATER *wat, METEO *met, PAR *par, TIMES *time, TOPO *top, SOIL *sl){

	long i, j, r=0, c=0;
	double w;
	
	for (j=1; j<=par->total_pixel; j++) {

		if(par->output_soil->co[i_sim]>0){
			r = top->rc_cont->co[j][1];
			c = top->rc_cont->co[j][2];
			if(strcmp(files[fpnet] , string_novalue) != 0) sl->Pnetcum->co[j] += wat->Pnet->co[r][c];
			if(strcmp(files[fevap] , string_novalue) != 0){
				for (i=1; i<=Nl; i++) {
					sl->ETcum->co[j] += sl->ET->co[i][r][c];
				}
			}
		}
			
		if(par->output_snow->co[i_sim]>0){
			if(strcmp(files[fsnowmelt] , string_novalue) != 0) snow->MELTED->co[j] += snow->melted->co[j];			
			if(strcmp(files[fsnowsubl] , string_novalue) != 0) snow->SUBL->co[j] += snow->subl->co[j];
			if(strcmp(files[fsndur] , string_novalue) != 0){ 
				if(snow->yes->co[j] == 1) snow->t_snow->co[j] += Dt/secinday; 
			}
		}
		
		if(par->max_glac_layers>0 && par->output_glac->co[i_sim]>0){
			if(strcmp(files[fglacmelt] , string_novalue) != 0) glac->MELTED->co[j] += glac->melted->co[j];
			if(strcmp(files[fglacsubl] , string_novalue) != 0) glac->SUBL->co[j] += glac->subl->co[j];
		}
		
		if(par->output_surfenergy->co[i_sim]>0){
			if(strcmp(files[fradnet] , string_novalue) != 0) egy->Rn_mean->co[j] += (egy->SW->co[j]+egy->LW->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradLWin] , string_novalue) != 0) egy->LWin_mean->co[j] += (egy->LWin->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradLW] , string_novalue) != 0) egy->LW_mean->co[j] += (egy->LW->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradSW] , string_novalue) != 0)  egy->SW_mean->co[j] += (egy->SW->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradSWin] , string_novalue) != 0) egy->Rswdown_mean->co[j] += (egy->SWin->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradSWinbeam] , string_novalue) != 0) egy->Rswbeam_mean->co[j] += (egy->SWinb->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fG] , string_novalue) != 0) egy->SEB_mean->co[j] += (egy->G->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fH] , string_novalue) != 0) egy->H_mean->co[j] += (egy->H->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fLE] , string_novalue) != 0) egy->ET_mean->co[j] += (egy->LE->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fTs] , string_novalue) != 0) egy->Ts_mean->co[j] += (egy->Ts->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			
			if(strcmp(files[fshadow] , string_novalue) != 0){
				if (egy->shad->co[j] >= 0) egy->nDt_sun->co[j] ++;
				if (egy->shad->co[j] == 0) egy->nDt_shadow->co[j] ++;
			}
		}
		if(par->output_meteo->co[i_sim]>0){
			if(strcmp(files[fprec] , string_novalue) != 0){
				wat->PrTOT_mean->co[j] = wat->Pt->co[j];
				wat->PrSNW_mean->co[j] = wat->Ps->co[j];
			}
		}	
		
		if(time->JD_plots->nh > 1 && W>0){
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->Hgplot->co[j] += egy->Hgp->co[j];
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0) egy->Hvplot->co[j] += egy->Hvp->co[j];
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LEgplot->co[j] += egy->LEgp->co[j];
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0) egy->LEvplot->co[j] += egy->LEvp->co[j];
			if(strcmp(files[pSWin] , string_novalue) != 0) egy->SWinplot->co[j] += egy->SWinp->co[j];
			if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->SWgplot->co[j] += egy->SWgp->co[j];
			if(strcmp(files[pSWv] , string_novalue) != 0) egy->SWvplot->co[j] += egy->SWvp->co[j];
			if(strcmp(files[pLWin] , string_novalue) != 0) egy->LWinplot->co[j] += egy->LWinp->co[j];
			if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LWgplot->co[j] += egy->LWgp->co[j];
			if(strcmp(files[pLWv] , string_novalue) != 0) egy->LWvplot->co[j] += egy->LWvp->co[j];
			if(strcmp(files[pTs] , string_novalue) != 0) egy->Tsplot->co[j] += egy->Tsp->co[j];
			if(strcmp(files[pTg] , string_novalue) != 0) egy->Tgplot->co[j] += egy->Tgp->co[j];
			if(strcmp(files[pD] , string_novalue) != 0) snow->Dplot->co[j] += W * DEPTH(top->rc_cont->co[j][1], top->rc_cont->co[j][2], snow->S->lnum, snow->S->Dzl);
			if(strcmp(files[pTa] , string_novalue) != 0) met->Taplot->co[j] += W * met->Tgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]];
			if(strcmp(files[pRH] , string_novalue) != 0) met->RHplot->co[j] += W * met->RHgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]];
			if(strcmp(files[pVspd] , string_novalue) != 0 || strcmp(files[pVdir] , string_novalue) != 0){
				met->Vxplot->co[j] -= W * met->Vgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]] * sin(met->Vdir->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]]*Pi/180.);
				met->Vyplot->co[j] -= W * met->Vgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]] * cos(met->Vdir->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]]*Pi/180.);
			}
			
		}
		if (par->state_pixel==1){
			if (par->jplot->co[j] > 0 && par->Dtplot_point->co[i_sim]>0){
				for(i=0; i<otot; i++){
					odpnt[i][par->jplot->co[j]-1] += odp[i][par->jplot->co[j]-1];
				}
			}
			if (par->jplot->co[j] > 0){
				for(i=1;i<=Nl;i++){
					r = top->rc_cont->co[j][1];
					c = top->rc_cont->co[j][2];
					if(par->Tzrun == 1) sl->Tzrun->co[par->jplot->co[j]][i] += sl->SS->T->co[i][j] * Dt / ((par->end_date->co[i_sim] - par->init_date->co[i_sim])*86400.);
					if(par->Tzmaxrun == 1){if(sl->Tzmaxrun->co[par->jplot->co[j]][i] < sl->SS->T->co[i][j]) sl->Tzmaxrun->co[par->jplot->co[j]][i] = sl->SS->T->co[i][j];}
					if(par->Tzminrun == 1){if(sl->Tzminrun->co[par->jplot->co[j]][i] > sl->SS->T->co[i][j]) sl->Tzminrun->co[par->jplot->co[j]][i] = sl->SS->T->co[i][j];}
					if(par->wzrun == 1 || par->wzmaxrun == 1 || par->wzminrun == 1){
						w = (sl->SS->thi->co[i][j] + sl->th->co[i][j]) * sl->pa->co[sl->type->co[r][c]][jdz][i];
						if(par->wzrun == 1) sl->wzrun->co[par->jplot->co[j]][i] += w * Dt / ((par->end_date->co[i_sim] - par->init_date->co[i_sim])*86400.);
						if(par->wzmaxrun == 1){if(sl->wzmaxrun->co[par->jplot->co[j]][i] < w) sl->wzmaxrun->co[par->jplot->co[j]][i] = w;}
						if(par->wzminrun == 1){if(sl->wzminrun->co[par->jplot->co[j]][i] > w) sl->wzminrun->co[par->jplot->co[j]][i] = w;}
					}
				}
				if(par->SWErun == 1) {
					w = get_SWE(r, c, snow->S->lnum, snow->S->w_ice, snow->S->w_liq);
					sl->SWErun->co[par->jplot->co[j]][1] += w * Dt / ((par->end_date->co[i_sim] - par->init_date->co[i_sim])*86400.);
					if (sl->SWErun->co[par->jplot->co[j]][2]<w) sl->SWErun->co[par->jplot->co[j]][2]=w;
					if (sl->SWErun->co[par->jplot->co[j]][3]>w) sl->SWErun->co[par->jplot->co[j]][3]=w;
				}

			}			
		}
	}
	
	if (par->state_basin==1){
		if (par->Dtplot_basin->co[i_sim]>0){
			for(i=0; i<ootot; i++){
				odbsn[i] += odb[i];
			}
		}
	}
	
}
	
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void print_run_average(SOIL *sl, TOPO *top, PAR *par){
	
	long r,c,j,l,n;
	FILE *f;
	char *temp, *name;
	char rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};

	n = Fminlong(par->Nl_spinup->co[i_sim],Nl);
	
	if(par->recover > 0) write_suffix(rec, par->recover, 4);
	if(par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

	if(par->Tzrun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fTrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fTrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fTrun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			fprintf(f, "%ld,%ld,%ld",i_sim,i_run,par->IDpoint->co[j]);
			for (l=1; l<=n; l++){
				fprintf(f, ",%f",sl->Tzrun->co[j][l]);
			}
			for (l=n+1; l<=Nl; l++) {
				r = par->rc->co[j][1];
				c = par->rc->co[j][2];
				fprintf(f, ",%f",sl->SS->T->co[l][top->j_cont[r][c]]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
	
	
	if(par->wzrun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fwrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fwrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fwrun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			for (l=n+1; l<=Nl; l++) {
				r = par->rc->co[j][1];
				c = par->rc->co[j][2];
			}
		}
		fclose(f);
	}
		
	
	if(par->Tzmaxrun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fTmaxrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fTmaxrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fTmaxrun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			fprintf(f, "%ld,%ld,%ld",i_sim,i_run,par->IDpoint->co[j]);
			for (l=1; l<=n; l++){
				fprintf(f, ",%f",sl->Tzmaxrun->co[j][l]);
			}
			for (l=n+1; l<=Nl; l++) {
				r = par->rc->co[j][1];
				c = par->rc->co[j][2];
				fprintf(f, ",%f",sl->SS->T->co[l][top->j_cont[r][c]]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
	
	
	if(par->wzmaxrun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fwmaxrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fwmaxrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fwmaxrun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			fprintf(f, "%ld,%ld,%ld",i_sim,i_run,par->IDpoint->co[j]);
			for (l=1; l<=n; l++){
				fprintf(f, ",%f",sl->wzmaxrun->co[j][l]);
			}
			for (l=n+1; l<=Nl; l++) {
				r = par->rc->co[j][1];
				c = par->rc->co[j][2];
				fprintf(f, ",%f",(sl->SS->thi->co[l][top->j_cont[r][c]]+sl->th->co[l][top->j_cont[r][c]])*sl->pa->co[sl->type->co[r][c]][jdz][l]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

	
	if(par->Tzminrun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fTminrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fTminrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fTminrun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			fprintf(f, "%ld,%ld,%ld",i_sim,i_run,par->IDpoint->co[j]);
			for (l=1; l<=n; l++){
				fprintf(f, ",%f",sl->Tzminrun->co[j][l]);
			}
			for (l=n+1; l<=Nl; l++) {
				r = par->rc->co[j][1];
				c = par->rc->co[j][2];
				fprintf(f, ",%f",sl->SS->T->co[l][top->j_cont[r][c]]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

	
	if(par->wzminrun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fwminrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fwminrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fwminrun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			fprintf(f, "%ld,%ld,%ld",i_sim,i_run,par->IDpoint->co[j]);
			for (l=1; l<=n; l++){
				fprintf(f, ",%f",sl->wzminrun->co[j][l]);
			}
			for (l=n+1; l<=Nl; l++) {
				r = par->rc->co[j][1];
				c = par->rc->co[j][2];
				fprintf(f, ",%f",(sl->SS->thi->co[l][top->j_cont[r][c]]+sl->th->co[l][top->j_cont[r][c]])*sl->pa->co[sl->type->co[r][c]][jdz][l]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
	

	if(par->dUzrun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fdUrun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fdUrun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fdUrun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			fprintf(f, "%ld,%ld,%ld",i_sim,i_run,par->IDpoint->co[j]);
			for (l=1; l<=Nl; l++){
				fprintf(f, ",%f",sl->dUzrun->co[j][l]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}
		
	
 	if(par->SWErun == 1){
		if (par->recover>0) {
			temp = join_strings(files[fSWErun], rec);
			name = join_strings(temp, textfile);
			free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fSWErun], crec);
			name = join_strings(temp, textfile);
			free(temp);
		}else {
			name = join_strings(files[fSWErun], textfile);
		}
		f = fopen(name, "a");
		for (j=1; j<=par->rc->nrh; j++) {
			fprintf(f, "%ld,%ld,%ld",i_sim,i_run,par->IDpoint->co[j]);
			for (l=1; l<=3; l++){
				fprintf(f, ",%f",sl->SWErun->co[j][l]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}	

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void init_run(SOIL *sl, PAR *par){
	
	long j, l;

	if (par->state_pixel == 1) {
								
		for (j=1; j<=par->rc->nrh; j++) {
			for (l=1; l<=Nl; l++){
				if(par->Tzrun == 1) sl->Tzrun->co[j][l] = 0.;
				if(par->wzrun == 1) sl->wzrun->co[j][l] = 0.;
				if(par->dUzrun == 1) sl->dUzrun->co[j][l] = 0.;
				if(par->Tzmaxrun == 1) sl->Tzmaxrun->co[j][l] = -1.E99;
				if(par->Tzminrun == 1) sl->Tzminrun->co[j][l] = 1.E99;
				if(par->wzmaxrun == 1) sl->wzmaxrun->co[j][l] = -1.E99;
				if(par->wzminrun == 1) sl->wzminrun->co[j][l] = 1.E99;
			}
			if(par->SWErun == 1){
				sl->SWErun->co[j][1] = 0.;
				sl->SWErun->co[j][2] = -1.E99;				
				sl->SWErun->co[j][3] = 1.E99;
			}
		}
	}
}
				
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void end_period_1D(SOIL *sl, TOPO *top, PAR *par){
	
	long j, l, sy, n, m;
	double k, kn, thwn, thin, psin, Tn, T0n, z, T;
	double Ptlow=0., thwlow=0., thilow=0., Tlow=0.;
	FILE *f;
	
	n = Fminlong(par->Nl_spinup->co[i_sim],Nl);
	
	for (j=1; j<=par->total_pixel; j++) {
		sy = sl->type->co[1][j];
		
		//spinned up soil portion (above) -> == 1 (instanatenous values) == 2 (averaged values)
		if (par->newperiodinit == 2) {
			for (l=1; l<=n; l++) {
				sl->SS->T->co[l][j] = sl->Tzrun->co[j][l];
				sl->Ptot->co[l][j] = psi_from_theta(sl->wzrun->co[j][l]/sl->pa->co[sy][jdz][l], 0., l, sl->pa->co[sy], PsiMin);
				sl->SS->P->co[l][j] = Fmin(Psif(sl->Tzrun->co[j][l]),sl->Ptot->co[l][j]);
				sl->th->co[l][j] = theta_from_psi(sl->SS->P->co[l][j], 0., l, sl->pa->co[sy], PsiMin);
				sl->SS->thi->co[l][j] = sl->wzrun->co[j][l]/sl->pa->co[sy][jdz][l] - sl->th->co[l][j];
			}
			Tlow = sl->SS->T->co[n][j];
			Ptlow = sl->Ptot->co[n][j];
			thwlow = sl->th->co[n][j];
			thilow = sl->SS->thi->co[n][j];
		}else if (par->newperiodinit == 1) {
			Tlow = sl->Tzrun->co[j][n];
			Ptlow = psi_from_theta(sl->wzrun->co[j][n]/sl->pa->co[sy][jdz][n], 0., n, sl->pa->co[sy], PsiMin);
			thwlow = theta_from_psi(Fmin(Psif(Tlow),Ptlow), 0., n, sl->pa->co[sy], PsiMin);
			thilow = sl->wzrun->co[j][n]/sl->pa->co[sy][jdz][n] - thwlow;
		}

		z = 0.;
		for (l=1; l<=n; l++) {
			z += sl->pa->co[sy][jdz][l];
		}
		
		//non-spinned up soil portion (below) -> interpolation (== 1 or 2)
		for (l=n+1; l<=Nl; l++) {
			z += sl->pa->co[sy][jdz][l];
			
			//above the weir set equal to the value of the deepest node of the spinned up soil portion, below the weir assumed saturation
			if(z<top->BC_DepthFreeSurface->co[j]){
				sl->Ptot->co[l][j] = Ptlow;
			}else {
				sl->Ptot->co[l][j] = z - top->BC_DepthFreeSurface->co[j];
			}
	
			//thermal conductivity of layer above
			if (l-1 == n) {
				k = k_thermal(0, 1, thwlow, thilow, sl->pa->co[sy][jsat][l-1], sl->pa->co[sy][jkt][l-1]);
				T = Tlow;
			}else {
				k = k_thermal(0, 1, sl->th->co[l-1][j], sl->SS->thi->co[l-1][j], sl->pa->co[sy][jsat][l-1], sl->pa->co[sy][jkt][l-1]);
				T = sl->SS->T->co[l-1][j];
			}
			
			//iterative loop to calculate thermal conductivity of the layer l
			m = 0;
			Tn = T;//first guess
			
			do{
				psin = Fmin(Psif(Tn), sl->Ptot->co[l][j]);
				thwn = theta_from_psi(psin , 0., l, sl->pa->co[sy], PsiMin);
				thin = theta_from_psi(sl->Ptot->co[l][j], 0., l, sl->pa->co[sy], PsiMin) - thwn;
				kn = k_thermal(0, 1, thwn, thin, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jkt][l]);
				
				T0n = Tn;
				
				Tn = T + par->Fboundary * 1.E-3 * (sl->pa->co[sy][jdz][l-1]/k + sl->pa->co[sy][jdz][l]/kn);
				m++;
				
			}while(fabs(Tn-T0n)>0.0001 && m<100);
			
			if (m==100) {
				f = fopen(FailedRunFile, "w");
				fprintf(f, "No convergence in the initial condition in the new simulation period\n");
				fclose(f);		
				t_error("Fatal Error! Geotop is closed. See failing report.");
			}
						
			sl->SS->T->co[l][j] = Tn;
			sl->SS->P->co[l][j] = psin;
			sl->th->co[l][j] = thwn;
			sl->SS->thi->co[l][j] = thin;
				
		}
		
		//print	
		/*for (l=1; l<=Nl; l++) {
			printf("j:%ld l:%ld Pt:%f T:%f\n",j,l,sl->Ptot->co[l][j],sl->SS->T->co[l][j]);
		}*/
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void change_grid(long previous_sim, long next_sim, PAR *par, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet){
	
	long n_previous, n_next, l, r, i, j;
	
	n_previous = Fminlong(par->Nl_spinup->co[previous_sim],Nl);
	n_next = Fminlong(par->Nl_spinup->co[next_sim],Nl);
		
	if (n_previous != n_next){
		
		//deallocate vectors with n_previous component
		free_longmatrix(top->lrc_cont);
		for(l=0;l<=n_previous;l++){
			for(r=1;r<=Nr;r++){
				free(top->i_cont[l][r]);
			}
			free(top->i_cont[l]);
		}
		free(top->i_cont);
		
		free_doublevector(wat->Lx);
		free_longvector(top->Li);
		free_longvector(top->Lp);
		free_doublevector(wat->H0);
		free_doublevector(wat->H1);
		free_doublevector(wat->dH);
		free_doublevector(wat->B);
		free_doublevector(wat->df);				  
		
		//reallocate vectors with n_next components
		top->i_cont=(long ***)malloc((n_next+1)*sizeof(long**));
		for(l=0;l<=n_next;l++){
			top->i_cont[l]=(long **)malloc((Nr+1)*sizeof(long*));
			for(r=1;r<=Nr;r++){
				top->i_cont[l][r]=(long *)malloc((Nc+1)*sizeof(long));
			}
		}
		top->lrc_cont=new_longmatrix( (n_next+1)*par->total_pixel, 3);
		initialize_longmatrix(top->lrc_cont, 0);
		i_lrc_cont(land->LC, top->i_cont, top->lrc_cont, n_next, Nr, Nc);
		
		if (par->point_sim != 1) {
			cont_nonzero_values_matrix2(&i, &j, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel, par->total_channel, n_next);
			top->Li = new_longvector(i);
			top->Lp = new_longvector(j);
			wat->Lx = new_doublevector(i);	
			cont_nonzero_values_matrix3(top->Lp, top->Li, cnet, land->LC, top->lrc_cont, top->i_cont, par->total_pixel, par->total_channel, n_next);
		}else {
			i = n_next;
			j = n_next+1;
			top->Li = new_longvector(i);
			top->Lp = new_longvector(j);
			wat->Lx = new_doublevector(i);	
			for (l=1; l<=n_next; l++) {
				top->Li->co[l] = l+1;
				top->Lp->co[l] = l;
			}
			top->Lp->co[j] = i;
		}
		
		wat->H0 = new_doublevector(j);
		wat->H1 = new_doublevector(j);
		wat->dH = new_doublevector(j);
		wat->B = new_doublevector(j);
		wat->f = new_doublevector(j);
		wat->df = new_doublevector(j);
		
	}
		
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
