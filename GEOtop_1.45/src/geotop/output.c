
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
    
#include "struct.geotop.h"
#include "output.h"
#include "pedo.funct.h"
#include "../libraries/geomorphology/networks.h"
#include "../libraries/ascii/rw_maps.h"
#include "constants.h"
#include "../libraries/ascii/extensions.h"
#include "times.h"
#include "energy.balance.h"
#include "input.h"
#include "../libraries/ascii/tabs.h"
#include "vegetation.h"
#include "tables.h"
#include "snow.h"
#include "../libraries/ascii/init.h"
#include "water.balance.h"


#include <time.h>
 
extern long number_novalue, number_absent;
extern char *string_novalue;

extern char *WORKING_DIRECTORY;

extern T_INIT *UV;
extern char **files, *logfile;
extern long Nl, Nr, Nc;

extern double t_meteo, t_energy, t_water, t_sub, t_sup, t_out, t_blowingsnow;

extern double **outdata_point, *outdata_basin;
extern long *outputpoint, noutputpoint, *outputbasin, noutputbasin, *outputsnow, noutputsnow;
extern long *outputglac, noutputglac, *outputsoil, noutputsoil;
extern char **headerpoint, **headerbasin, **headersnow, **headerglac, **headersoil;

extern char *keywords_num[num_par_number] , *keywords_char[num_par_char];

extern FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav, *ffice, *fficeav, *ffsnow, *ffglac;

extern long i_sim, i_run;

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_output(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac, METEO *met)
				  
{
	/*internal auxiliary variables:*/
	long i,j,r,c,l; /*counters*/
	long n_file;      /*number of file of the type "TETAxySSSlZZ"(i.e. number of the basin-time-step)*/
	char SSSS[ ]={"SSSS"};
	char *name;	
	char *temp1,*temp2;  
	FILE *f, *flog;
	
	//time variables
	static time_t start_time, stop_time;
	double percent_done, elapsed_time, remaining_time, total_time;
	short first_column;
	double JD, JDfrom0;
	long day, month, year, hour, minute;
	
	// static double Qsub_ch, Qsup_ch;
	static double Vtot; /*averaged output flows*/
	static double Vsub, Vsup, Voutland;
	static long isavings;
	static double mass_error_tot;
	static double t_discharge, t_basin, t_point;
	double Vchannel;
	
	/*internal variables to memorize input par:*/
	double time_max;      
	
	//other variables	
	DOUBLEMATRIX *M, *K;
	DOUBLETENSOR *T;
	double D;
	
	flog = fopen(logfile, "a");
	
	//initialize static variables
	if(times->time < 1.E-5){
		time( &start_time );
		mass_error_tot=0.;
		Vtot=0.;
		Vsub=0.;
		Vsup=0.;
		Voutland=0.;
		t_discharge=0.;
		t_basin=0.;
		t_point=0.;
	}
	
	//Assignment to some internal variables of some input par
	time_max=(par->end_date->co[i_sim] - par->init_date->co[i_sim])*86400.;//seconds
 	
	//Time indices
	JDfrom0 = convert_tfromstart_JDfrom0(times->time+par->Dt, par->init_date->co[i_sim]);
	convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
	convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);
	
	
	//DISCHARGE
	//****************************************************************************************************************
	//****************************************************************************************************************
	
	if(par->state_discharge == 1 && par->Dtplot_discharge->co[i_sim] > 1.E-5 && strcmp(files[fQ] , string_novalue) != 0){

		t_discharge += par->Dt;

		Vtot += cnet->Vout;

		for(l=1;l<=par->total_channel;l++){
			Vsub += cnet->Vsub->co[l];
			Vsup += cnet->Vsup->co[l];
		}	

		Voutland += wat->Voutland;
		

		if (fabs(t_discharge - par->Dtplot_discharge->co[i_sim]) < 1.E-5){/*Print the outlet flows*/
			
			Vchannel = 0.;
			for(l=1;l<=par->total_channel;l++){
				r = cnet->r->co[l];
				c = cnet->c->co[l];
				Vchannel +=1.E-3*find_hsup(cnet->P->co[0][l], top->slope->co[r][c]) * UV->U->co[1] * par->w_dx * cnet->length->co[l];
			}
			
			name=join_strings(files[fQ],textfile);
			f=fopen(name,"a");
			fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
			fprintf(f,",%f,%f,%f",(times->time+par->Dt)/secinday,JDfrom0,JD);   
			fprintf(f,",%e,%e,%e,%e,%e\n",Vtot/(double)par->Dtplot_discharge->co[i_sim],Vsup/(double)par->Dtplot_discharge->co[i_sim],
					Vsub/(double)par->Dtplot_discharge->co[i_sim],Vchannel,Voutland/(double)par->Dtplot_discharge->co[i_sim]);
			fclose(f);
			free(name);
			
			t_discharge = 0.0;
			
			Vtot = 0.0;
			Vsub = 0.0;
			Vsup = 0.0;
			Voutland = 0.0;
			
			
		}
	}
	
	//DATA POINT
	//****************************************************************************************************************
	//****************************************************************************************************************
	
	
	if(par->Dtplot_point->co[i_sim] > 1.E-5){
		
		t_point += par->Dt;
		
		if(par->state_pixel == 1){
			for(i=1;i<=par->rc->nrh;i++){
				
				r=par->rc->co[i][1];
				c=par->rc->co[i][2];	
				
				for(l=1;l<=Nl;l++){
					if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] += sl->T->co[l][r][c]*(par->Dt/par->Dtplot_point->co[i_sim]);
					if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] += sl->th->co[l][r][c]*(par->Dt/par->Dtplot_point->co[i_sim]);
					if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thicezavplot->co[i][l] += sl->thice->co[l][r][c]*(par->Dt/par->Dtplot_point->co[i_sim]);
				}
				
				outdata_point[othawed][i-1]+=find_activelayerdepth(r, c, sl)*(par->Dt/par->Dtplot_point->co[i_sim]);				
				outdata_point[owtable][i-1]+=find_watertabledepth(r, c, sl)*(par->Dt/par->Dtplot_point->co[i_sim]);	
				
			}
		}
			
		//Print of pixel-output every times->n_pixel time step
		if (fabs(t_point - par->Dtplot_point->co[i_sim])<1.E-5){
			
			if(par->state_pixel == 1){
			
				for(i=1;i<=par->rc->nrh;i++){
					
					write_suffix(SSSS, par->IDpoint->co[i], 0);
					r=par->rc->co[i][1];
					c=par->rc->co[i][2];
										
					//soil data
					for(l=1;l<=Nl;l++){
						if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot->co[i][l] = sl->T->co[l][r][c];
						if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot->co[i][l] = sl->Ptot->co[l][r][c];
						if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot->co[i][l] = sl->P->co[l][r][c];
						if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot->co[i][l] = sl->th->co[l][r][c];
						if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thicezplot->co[i][l] = sl->thice->co[l][r][c];
					}				
					
					//snow data
					if(snow->S->lnum->co[r][c]>0){
						outdata_point[osnowdepth][i-1] = 0.0;
						outdata_point[oSWE][i-1] = 0.0;
						outdata_point[osnowT][i-1] = 0.0;
						for(l=1;l<=snow->S->lnum->co[r][c];l++){
							outdata_point[osnowdepth][i-1] += snow->S->Dzl->co[l][r][c];			   
							outdata_point[oSWE][i-1] += 1.0E+3*(snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c])/rho_w;
							outdata_point[osnowT][i-1] += snow->S->T->co[l][r][c]*snow->S->Dzl->co[l][r][c];
						}
						outdata_point[osnowdens][i-1] = outdata_point[oSWE][i-1]*rho_w/outdata_point[osnowdepth][i-1];		
						outdata_point[osnowT][i-1] /= outdata_point[osnowdepth][i-1];   
					}else{
						outdata_point[osnowdepth][i-1] = 0.0;
						outdata_point[oSWE][i-1] = 0.0;
						outdata_point[osnowdens][i-1] = (double)number_novalue;
						outdata_point[osnowT][i-1] = (double)number_novalue;
					}
					
					//glacier data
					if(par->glaclayer_max>0){
						if(glac->G->lnum->co[r][c]>0){
							outdata_point[oglacdepth][i-1] = 0.0;
							outdata_point[oGWE][i-1] = 0.0;
							outdata_point[oglacT][i-1] = 0.0;
							for(l=1;l<=glac->G->lnum->co[r][c];l++){
								outdata_point[oglacdepth][i-1] += glac->G->Dzl->co[l][r][c];			   
								outdata_point[oGWE][i-1] += 1.0E+3*(glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c])/rho_w;
								outdata_point[oglacT][i-1] += glac->G->T->co[l][r][c]*glac->G->Dzl->co[l][r][c];
							}
							outdata_point[oglacdens][i-1] = outdata_point[oGWE][i-1]*rho_w/outdata_point[oglacdepth][i-1];		
							outdata_point[oglacT][i-1] /= outdata_point[oglacdepth][i-1];   
						}else{
							outdata_point[oglacdepth][i-1] = 0.0;
							outdata_point[oGWE][i-1] = 0.0;
							outdata_point[oglacdens][i-1] = (double)number_novalue;
							outdata_point[oglacT][i-1] = (double)number_novalue;
						}
					}
					
					//Point data			
					if(strcmp(files[fpoint] , string_novalue) != 0){
						temp1=join_strings(files[fpoint],SSSS);
						name=join_strings(temp1,textfile);
						f=fopen(name,"a");				
						first_column=1;
						
						for (j=0; j<noutputpoint; j++) {
							if(first_column==0){
								fprintf(f,",");
							}else {
								first_column = 0;
							}
							if (outputpoint[j] >= 0) {
								if (outputpoint[j] == odate12) {
									fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (outputpoint[j] == oJDfrom0) {
									fprintf(f, "%f",JDfrom0);
								}else if (outputpoint[j] == odaysfromstart) {
									fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
								}else if (outputpoint[j] == operiod) {
									fprintf(f, "%ld",i_sim);
								}else if (outputpoint[j] == orun) {
									fprintf(f, "%ld",i_run);
								}else if (outputpoint[j] == opoint) {
									fprintf(f, "%ld",par->IDpoint->co[i]);
								}else {
									fprintf(f,"%f",outdata_point[outputpoint[j]][i-1]);
								}
							}else {
								fprintf(f,"%ld",number_novalue);
							}
						}
						fprintf(f,"\n");
						fclose(f);
						free(temp1);
						free(name);			
					}
					
					if(strcmp(files[fpointwriteend] , string_novalue) != 0){
						
						first_column=1;
						for (j=0; j<noutputpoint; j++) {
							if(first_column==0){
								fprintf(ffpoint,",");
							}else {
								first_column = 0;
							}
							if (outputpoint[j] >= 0) {
								if (outputpoint[j] == odate12) {
									fprintf(ffpoint,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (outputpoint[j] == oJDfrom0) {
									fprintf(ffpoint, "%f",JDfrom0);
								}else if (outputpoint[j] == odaysfromstart) {
									fprintf(ffpoint, "%f",JDfrom0-par->init_date->co[i_sim]);
								}else if (outputpoint[j] == operiod) {
									fprintf(ffpoint, "%ld",i_sim);
								}else if (outputpoint[j] == orun) {
									fprintf(ffpoint, "%ld",i_run);		
								}else if (outputpoint[j] == opoint) {
									fprintf(ffpoint, "%ld",par->IDpoint->co[i]);
								}else {
									fprintf(ffpoint,"%f",outdata_point[outputpoint[j]][i-1]);
								}
							}else {
								fprintf(ffpoint,"%ld",number_novalue);
							}
						}
						fprintf(ffpoint,"\n");
					}
					
					
					//Snow
					if(strcmp(files[fsnz] , string_novalue) != 0){
						temp1=join_strings(files[fsnz],SSSS);
						name=join_strings(temp1,textfile);
						f=fopen(name,"a");
						
						first_column=1;
						for (j=0; j<noutputsnow; j++) {
							if(first_column==0){
								fprintf(f,",");
							}else {
								first_column = 0;
							}
							if (outputsnow[j] >= 0) {
								if (outputsnow[j] == 0) {
									fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (outputsnow[j] == 1) {
									fprintf(f, "%f",JDfrom0);
								}else if (outputsnow[j] == 2) {
									fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
								}else if (outputsnow[j] == 3) {
									fprintf(f, "%ld",i_sim);
								}else if (outputsnow[j] == 4) {
									fprintf(f, "%ld",i_run);								
								}else if (outputsnow[j] == 5) {
									fprintf(f, "%ld",par->IDpoint->co[i]);
								}else if (outputsnow[j] <= 5 + 1*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 0*par->snowlayer_max;
									fprintf(f, "%f",snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
								}else if (outputsnow[j] <= 5 + 2*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 1*par->snowlayer_max;
									fprintf(f, "%f",snow->S->Dzl->co[l][r][c]);
								}else if (outputsnow[j] <= 5 + 3*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 2*par->snowlayer_max;
									if (snow->S->Dzl->co[l][r][c] > 0) {
										fprintf(f, "%f",(snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c])/(1.E-3*snow->S->Dzl->co[l][r][c]));
									}else {
										fprintf(f, "%f",(double)number_novalue);
									}
								}else if (outputsnow[j] <= 5 + 4*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 3*par->snowlayer_max;
									fprintf(f, "%f",snow->S->T->co[l][r][c]);
								}else if (outputsnow[j] <= 5 + 5*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 4*par->snowlayer_max;
									if (snow->S->Dzl->co[l][r][c] > 0) {
										fprintf(f, "%f",snow->S->w_ice->co[l][r][c]/(1.0E-3*snow->S->Dzl->co[l][r][c]*rho_i));
									}else {
										fprintf(f, "%f",(double)number_novalue);
									}
								}else if (outputsnow[j] <= 5 + 6*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 5*par->snowlayer_max;
									if (snow->S->Dzl->co[l][r][c] > 0) {
										fprintf(f, "%f",snow->S->w_liq->co[l][r][c]/(1.0E-3*snow->S->Dzl->co[l][r][c]*rho_w));
									}else {
										fprintf(f, "%f",(double)number_novalue);
									}
								}
							}else {
								fprintf(f,"%f",(double)number_novalue);
							}
						}
						fprintf(f,"\n");
						fclose(f);
						free(temp1);
						free(name);
					}
					
					if(strcmp(files[fsnzwriteend] , string_novalue) != 0){
						first_column=1;
						for (j=0; j<noutputsnow; j++) {
							if(first_column==0){
								fprintf(ffsnow,",");
							}else {
								first_column = 0;
							}
							if (outputsnow[j] >= 0) {
								if (outputsnow[j] == 0) {
									fprintf(ffsnow,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (outputsnow[j] == 1) {
									fprintf(ffsnow, "%f",JDfrom0);
								}else if (outputsnow[j] == 2) {
									fprintf(ffsnow, "%f",JDfrom0-par->init_date->co[i_sim]);
								}else if (outputsnow[j] == 3) {
									fprintf(ffsnow, "%ld",i_sim);
								}else if (outputsnow[j] == 4) {
									fprintf(ffsnow, "%ld",i_run);								
								}else if (outputsnow[j] == 5) {
									fprintf(ffsnow, "%ld",par->IDpoint->co[i]);
								}else if (outputsnow[j] <= 5 + 1*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 0*par->snowlayer_max;
									fprintf(ffsnow, "%f",snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
								}else if (outputsnow[j] <= 5 + 2*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 1*par->snowlayer_max;
									fprintf(ffsnow, "%f",snow->S->Dzl->co[l][r][c]);
								}else if (outputsnow[j] <= 5 + 3*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 2*par->snowlayer_max;
									if (snow->S->Dzl->co[l][r][c] > 0) {
										fprintf(ffsnow, "%f",(snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c])/(1.E-3*snow->S->Dzl->co[l][r][c]));
									}else {
										fprintf(ffsnow, "%f",(double)number_novalue);
									}
								}else if (outputsnow[j] <= 5 + 4*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 3*par->snowlayer_max;
									fprintf(ffsnow, "%f",snow->S->T->co[l][r][c]);
								}else if (outputsnow[j] <= 5 + 5*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 4*par->snowlayer_max;
									if (snow->S->Dzl->co[l][r][c] > 0) {
										fprintf(ffsnow, "%f",snow->S->w_ice->co[l][r][c]/(1.0E-3*snow->S->Dzl->co[l][r][c]*rho_i));
									}else {
										fprintf(ffsnow, "%f",(double)number_novalue);
									}
								}else if (outputsnow[j] <= 5 + 6*par->snowlayer_max) {
									l = outputsnow[j] - 5 - 5*par->snowlayer_max;
									if (snow->S->Dzl->co[l][r][c] > 0) {
										fprintf(ffsnow, "%f",snow->S->w_liq->co[l][r][c]/(1.0E-3*snow->S->Dzl->co[l][r][c]*rho_w));
									}else {
										fprintf(ffsnow, "%f",(double)number_novalue);
									}
								}
							}else {
								fprintf(ffsnow,"%f",(double)number_novalue);
							}
						}
						fprintf(ffsnow,"\n");
					}
					
					
					//Glacier
					if(par->glaclayer_max>0  && strcmp(files[fglz] , string_novalue) != 0){
						temp1=join_strings(files[fglz],SSSS);
						name=join_strings(temp1,textfile);
						f=fopen(name,"a");
						
						first_column=1;
						for (j=0; j<noutputglac; j++) {
							if(first_column==0){
								fprintf(f,",");
							}else {
								first_column = 0;
							}
							if (outputglac[j] >= 0) {
								if (outputglac[j] == 0) {
									fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (outputglac[j] == 1) {
									fprintf(f, "%f",JDfrom0);
								}else if (outputglac[j] == 2) {
									fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
								}else if (outputglac[j] == 3) {
									fprintf(f, "%ld",i_sim);
								}else if (outputglac[j] == 4) {
									fprintf(f, "%ld",i_run);																
								}else if (outputglac[j] == 5) {
									fprintf(f, "%ld",par->IDpoint->co[i]);
								}else if (outputglac[j] <= 5 + 1*par->glaclayer_max) {
									l = outputglac[j] - 5 - 0*par->glaclayer_max;
									fprintf(f, "%f",glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
								}else if (outputglac[j] <= 5 + 2*par->glaclayer_max) {
									l = outputglac[j] - 5 - 1*par->glaclayer_max;
									fprintf(f, "%f",glac->G->Dzl->co[l][r][c]);
								}else if (outputglac[j] <= 5 + 3*par->glaclayer_max) {
									l = outputglac[j] - 5 - 2*par->glaclayer_max;
									if (glac->G->Dzl->co[l][r][c] > 0) {
										fprintf(f, "%f",(glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c])/(1.E-3*glac->G->Dzl->co[l][r][c]));
									}else {
										fprintf(f, "%f",(double)number_novalue);
									}
								}else if (outputglac[j] <= 5 + 4*par->glaclayer_max) {
									l = outputglac[j] - 5 - 3*par->glaclayer_max;
									fprintf(f, "%f",glac->G->T->co[l][r][c]);
								}else if (outputglac[j] <= 5 + 5*par->glaclayer_max) {
									l = outputglac[j] - 5 - 4*par->glaclayer_max;
									if (glac->G->Dzl->co[l][r][c] > 0) {
										fprintf(f, "%f",glac->G->w_ice->co[l][r][c]/(1.0E-3*glac->G->Dzl->co[l][r][c]*rho_i));
									}else {
										fprintf(f, "%f",(double)number_novalue);
									}
								}else if (outputglac[j] <= 5 + 6*par->glaclayer_max) {
									l = outputglac[j] - 5 - 5*par->glaclayer_max;
									if (glac->G->Dzl->co[l][r][c] > 0) {
										fprintf(f, "%f",glac->G->w_liq->co[l][r][c]/(1.0E-3*glac->G->Dzl->co[l][r][c]*rho_w));
									}else {
										fprintf(f, "%f",(double)number_novalue);
									}
								}
							}else {
								fprintf(f,"%f",(double)number_novalue);
							}
						}
						fprintf(f,"\n");
						fclose(f);
						free(temp1);
						free(name);
					}
					
					if(par->glaclayer_max>0  && strcmp(files[fglzwriteend] , string_novalue) != 0){
						first_column=1;
						for (j=0; j<noutputglac; j++) {
							if(first_column==0){
								fprintf(ffglac,",");
							}else {
								first_column = 0;
							}
							if (outputglac[j] >= 0) {
								if (outputglac[j] == 0) {
									fprintf(ffglac,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (outputglac[j] == 1) {
									fprintf(ffglac, "%f",JDfrom0);
								}else if (outputglac[j] == 2) {
									fprintf(ffglac, "%f",JDfrom0-par->init_date->co[i_sim]);
								}else if (outputglac[j] == 3) {
									fprintf(ffglac, "%ld",i_sim);
								}else if (outputglac[j] == 4) {
									fprintf(ffglac, "%ld",i_run);																
								}else if (outputglac[j] == 5) {
									fprintf(ffglac, "%ld",par->IDpoint->co[i]);
								}else if (outputglac[j] <= 5 + 1*par->glaclayer_max) {
									l = outputglac[j] - 5 - 0*par->glaclayer_max;
									fprintf(ffglac, "%f",glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
								}else if (outputglac[j] <= 5 + 2*par->glaclayer_max) {
									l = outputglac[j] - 5 - 1*par->glaclayer_max;
									fprintf(ffglac, "%f",glac->G->Dzl->co[l][r][c]);
								}else if (outputglac[j] <= 5 + 3*par->glaclayer_max) {
									l = outputglac[j] - 5 - 2*par->glaclayer_max;
									if (glac->G->Dzl->co[l][r][c] > 0) {
										fprintf(ffglac, "%f",(glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c])/(1.E-3*glac->G->Dzl->co[l][r][c]));
									}else {
										fprintf(ffglac, "%f",(double)number_novalue);
									}
								}else if (outputglac[j] <= 5 + 4*par->glaclayer_max) {
									l = outputglac[j] - 5 - 3*par->glaclayer_max;
									fprintf(ffglac, "%f",glac->G->T->co[l][r][c]);
								}else if (outputglac[j] <= 5 + 5*par->glaclayer_max) {
									l = outputglac[j] - 5 - 4*par->glaclayer_max;
									if (glac->G->Dzl->co[l][r][c] > 0) {
										fprintf(ffglac, "%f",glac->G->w_ice->co[l][r][c]/(1.0E-3*glac->G->Dzl->co[l][r][c]*rho_i));
									}else {
										fprintf(ffglac, "%f",(double)number_novalue);
									}
								}else if (outputglac[j] <= 5 + 6*par->glaclayer_max) {
									l = outputglac[j] - 5 - 5*par->glaclayer_max;
									if (glac->G->Dzl->co[l][r][c] > 0) {
										fprintf(ffglac, "%f",glac->G->w_liq->co[l][r][c]/(1.0E-3*glac->G->Dzl->co[l][r][c]*rho_w));
									}else {
										fprintf(ffglac, "%f",(double)number_novalue);
									}
								}
							}else {
								fprintf(ffglac,"%f",(double)number_novalue);
							}
						}
						fprintf(ffglac,"\n");
					}
					
					//sl output
					write_soil_output(i, par->IDpoint->co[i], par->init_date->co[i_sim], JDfrom0, JD, day, month, year, hour, minute, par->rc, sl, PsiMin);
				
					//initialize
					for(j=0;j<otot;j++) { outdata_point[j][i-1]=0.0; }

				}
				
			}
								
			percent_done = 100.0*(double)times->time/time_max;
			time( &stop_time );
			elapsed_time = difftime( stop_time, start_time );
			if( percent_done > 1.0e-6 ){
				total_time = elapsed_time * 100.0 / percent_done;
			}else{
				total_time = -999.9;
			}
			remaining_time = (total_time - elapsed_time);
			
			
			printf("%ld/%ld/%ld %ld:%02.0f %.2f%% - Elapsed (h:m:s) %2.0f:%02.0f:%02.0f - Remaining (h:m) %2.0f:%02.0f  \n",
				   day,month,year,hour,(float)minute,(100.0*(times->time+par->Dt)/time_max), 
				   floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
				   floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
				   floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
			
			
			fprintf(flog,"%ld/%ld/%ld %ld:%02.0f %.2f%% - Elapsed (h:m:s) %2.0f:%02.0f:%02.0f - Remaining (h:m) %2.0f:%02.0f  \n",
					day,month,year,hour,(float)minute,(100.0*(times->time+par->Dt)/time_max), 
					floor(elapsed_time / 3600.0), floor(((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60.),
					floor((((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. - floor( ((elapsed_time/3600)-floor(elapsed_time / 3600.0))*60. ))*60.),
					floor(remaining_time / 3600.0), floor(((remaining_time/3600)-floor(remaining_time / 3600.0))*60.) );
			
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
				name=join_strings(files[fbas],textfile);
				f=fopen(name,"a");
				first_column=1;
				for (j=0; j<noutputbasin; j++) {
					if(first_column==0){
						fprintf(f,",");
					}else {
						first_column = 0;
					}
					if (outputbasin[j] >= 0) {
						if (outputbasin[j] == oodate12) {
							fprintf(f, "%012.0f",convert_daymonthyearhourmin_dateeur12(day, month, year, hour, minute));
						}else if (outputbasin[j] == ooJDfrom0) {
							fprintf(f, "%f",JDfrom0);
						}else if (outputbasin[j] == oodaysfromstart) {
							fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
						}else {
							fprintf(f,"%f",outdata_basin[outputbasin[j]]);
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
				for (j=0; j<noutputbasin; j++) {
					if(first_column==0){
						fprintf(ffbas,",");
					}else {
						first_column = 0;
					}
					if (outputbasin[j] >= 0) {
						if (outputbasin[j] == oodate12) {
							fprintf(ffbas, "%012.0f",convert_daymonthyearhourmin_dateeur12(day, month, year, hour, minute));
						}else if (outputbasin[j] == ooJDfrom0) {
							fprintf(ffbas, "%f",JDfrom0);
						}else if (outputbasin[j] == oodaysfromstart) {
							fprintf(ffbas, "%f",JDfrom0-par->init_date->co[i_sim]);
						}else {
							fprintf(ffbas,"%f",outdata_basin[outputbasin[j]]);
						}
					}else {
						fprintf(ffbas,"%f",(double)number_novalue);
					}
				}
				fprintf(ffbas,"\n");
			}
			
			printf("\n%ld/%ld/%ld %ld:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
				   day,month,year,hour,(float)minute,JD,(long)(floor(times->time/86400))+1,
				   (100.0*(double)(times->time+par->Dt)/time_max));
			
			fprintf(flog,"\n%ld/%ld/%ld %ld:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
					day,month,year,hour,(float)minute,JD,(long)(floor(times->time/86400))+1,
					(100.0*(double)(times->time+par->Dt)/time_max));
			
			
			mass_error_tot += outdata_basin[oomasserror];
			printf(" t_meteo:%6.2f s, t_energy:%6.2f s, t_blowingsnow:%6.2f s, t_water:%6.2f s, t_sub:%6.2f s, t_sup:%6.2f s, t_out:%6.2f s\n",t_meteo,t_energy,t_blowingsnow,t_water,t_sub,t_sup,t_out);
			printf(" SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm\n\n",
				   outdata_basin[ooSW],outdata_basin[ooLW],outdata_basin[ooH],outdata_basin[ooLE],outdata_basin[oorainover],
				   outdata_basin[oosnowover],outdata_basin[oomasserror]*3600.0/par->Dtplot_basin->co[i_sim],mass_error_tot);
			
			fprintf(flog," t_meteo:%6.2f s, t_energy:%6.2f s, t_blowingsnow:%6.2f s, t_water:%6.2f s, t_sub:%6.2f s, t_sup:%6.2f s, t_out:%6.2f s\n",t_meteo,t_energy,t_blowingsnow,t_water,t_sub,t_sup,t_out);
			fprintf(flog," SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm\n\n",
					outdata_basin[ooSW],outdata_basin[ooLW],outdata_basin[ooH],outdata_basin[ooLE],outdata_basin[oorainover],
					outdata_basin[oosnowover],outdata_basin[oomasserror]*3600.0/par->Dtplot_basin->co[i_sim],mass_error_tot);
			
			
			
			for(j=0;j<ootot;j++){ 
				outdata_basin[j]=0.0; 
			}
			
			t_basin = 0.0;
			
		}
	}
	
	
	//DISTRIBUTED OUTPUTS
	//****************************************************************************************************************
	//****************************************************************************************************************
	//averaging properties
	if(par->output_meteo>0){
		
		for(r=1; r<=Nr; r++){
			for (c=1; c<=Nc; c++) {
				if( (long)land->LC->co[r][c] != number_novalue){
					
					met->Ta_mean->co[r][c]+=met->Tgrid->co[1][r][c]/((par->output_meteo*3600.0)/(par->Dt));
					met->Vspdmean->co[r][c]+=met->Vgrid->co[1][r][c]/((par->output_meteo*3600.0)/(par->Dt));
					met->Vdirmean->co[r][c]+=met->Vdir->co[1][r][c]/((par->output_meteo*3600.0)/(par->Dt));
					met->RHmean->co[r][c]+=met->RHgrid->co[1][r][c]/((par->output_meteo*3600.0)/(par->Dt));
					
					if(par->distr_stat==1){
						if(met->Ta_max->co[r][c]<met->Tgrid->co[1][r][c]) met->Ta_max->co[r][c]=met->Tgrid->co[1][r][c];
						if(met->Ta_min->co[r][c]>met->Tgrid->co[1][r][c]) met->Ta_min->co[r][c]=met->Tgrid->co[1][r][c];
					}
				}
			}
		}
	}
	
	if (par->output_soil>0) {
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0){
			for(r=1; r<=Nr; r++){
				for (c=1; c<=Nc; c++) {
					if( (long)land->LC->co[r][c] != number_novalue){
						for (l=1; l<=Nl; l++) {
							sl->T_av_tensor->co[l][r][c] += sl->T->co[l][r][c]/((par->output_soil*3600.0)/(par->Dt));
						}
					}
				}
			}
		}
	}	
	
	
	M=new_doublematrix(Nr, Nc);
	initialize_doublematrix(M, (double)number_novalue);
	
	//soil properties
	if(par->output_soil>0 && fmod(times->time+par->Dt,par->output_soil*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt+par->delay_day_recover)/(par->output_soil*3600.0));
		write_suffix(SSSS, n_file, 0);
		
		//theta liq tensor
		if(strcmp(files[fliq] , string_novalue) != 0) write_tensorseries2(n_file, files[fliq], 0, par->format_out, sl->th, UV, number_novalue);	
		
		//theta liq surface
		if(strcmp(files[fliqsup] , string_novalue) != 0){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c] = sl->th->co[1][r][c];
					}
				}
			}
			temp1=join_strings(files[fliqsup],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);	
			free(temp1);
		}
		
		//write T tensor
		if(strcmp(files[fT] , string_novalue) != 0) write_tensorseries2(n_file, files[fT], 0, par->format_out, sl->T, UV, number_novalue);
		
		//theta T surface
		if(strcmp(files[fTsup] , string_novalue) != 0){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c] = sl->T->co[1][r][c];
					}
				}
			}
			temp1=join_strings(files[fTsup],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);	
			free(temp1);
		}
		
		//write Tav tensor
		if(strcmp(files[fTav] , string_novalue) != 0) write_tensorseries2(n_file, files[fTav], 0, par->format_out, sl->T_av_tensor, UV, number_novalue);
		
		//theta Tav surface
		if(strcmp(files[fTavsup] , string_novalue) != 0){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c] = sl->T_av_tensor->co[1][r][c];
					}
				}
			}
			temp1=join_strings(files[fTavsup],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);	
			free(temp1);
		}
		
		//initialize T_av_tensor
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
		
		//theta_ice tensor
		if(strcmp(files[fice] , string_novalue) != 0) write_tensorseries2(n_file, files[fice], 0, par->format_out, sl->thice, UV, number_novalue);
		
		//theta_ice surface
		if(strcmp(files[ficesup] , string_novalue) != 0){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c] = sl->thice->co[1][r][c];
					}
				}
			}
			temp1=join_strings(files[ficesup],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);	
			free(temp1);
		}
		
		//write psi tensors
		if(strcmp(files[fpsi] , string_novalue) != 0){
			temp1=join_strings(files[fpsi],"TOT");
			write_tensorseries2(n_file, temp1, 0, par->format_out, sl->Ptot, UV, number_novalue);
			free(temp1);
			
			temp1=join_strings(files[fpsi],"LIQ");
			write_tensorseries2(n_file, temp1, 0, par->format_out, sl->P, UV, number_novalue);
			free(temp1);		
		}
		
		//calculate saturation front depth
		if( strcmp(files[fwtable] , string_novalue) != 0 ){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=find_watertabledepth(r, c, sl);				
					}
				}
			}
			temp1=join_strings(files[fwtable],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);	
			free(temp1);
		}
		
		//calculate active layer depth		   
		if( strcmp(files[fthawed] , string_novalue) != 0 ){
			
			K=new_doublematrix(Nr,Nc);
			initialize_doublematrix(K,(double)number_novalue);	
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						K->co[r][c]=find_activelayerdepth(r, c, sl);					
					}
				}
			}
			temp1=join_strings(files[fthawed],SSSS);
			write_map(temp1, 0, par->format_out, K, UV, number_novalue);	
			free(temp1);
			
			//calculate frost table depth
			if( strcmp(files[fwtable] , string_novalue) != 0 && strcmp(files[fftable] , string_novalue) != 0){
				for(r=1;r<=Nr;r++){
					for(c=1;c<=Nc;c++){
						if( (long)land->LC->co[r][c]!=number_novalue){
							//the deeper between the previous 2		
							M->co[r][c]=Fmax(K->co[r][c], M->co[r][c]);
						}
					}
				}
				temp1=join_strings(files[fftable],SSSS);
				write_map(temp1, 0, par->format_out, M, UV, number_novalue);	
				free(temp1);
			}
			
			free_doublematrix(K);
		}
		
		//WATER OVER THE SURFACE
		if( strcmp(files[fhsup] , string_novalue) != 0 ){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c] = Fmax(0, sl->P->co[0][r][c]);
					}
				}
			}
			
			temp1 = join_strings(files[fhsup],"LAND");
			temp2 = join_strings(temp1, SSSS);
			write_map(temp2, 0, par->format_out, M, UV, number_novalue);
			free(temp1);
			free(temp2);
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						if (cnet->ch->co[r][c]!=0) {
							M->co[r][c] = cnet->h_sup->co[cnet->ch->co[r][c]];
						}else {
							M->co[r][c] = number_novalue;
						}
					}
				}
			}
			
			temp1 = join_strings(files[fhsup],"CHANNEL");
			temp2 = join_strings(temp1, SSSS);
			write_map(temp2, 0, par->format_out, M, UV, number_novalue);
			free(temp1);
			free(temp2);		
			
			
			//contribution in volume to the channel by surface flow
			/*for(r=1;r<=Nr;r++){
			 for(c=1;c<=Nc;c++){
			 if( (long)land->LC->co[r][c]!=number_novalue){
			 if (cnet->ch->co[r][c]!=0) {
			 M->co[r][c] = cnet->Vsup_cum->co[cnet->ch->co[r][c]];
			 }else {
			 M->co[r][c] = number_novalue;
			 }
			 }
			 }
			 }
			 
			 temp1 = join_strings(files[fhsup],"Vsup_to_channel");
			 temp2 = join_strings(temp1, SSSS);
			 write_map(temp2, 0, par->format_out, M, UV, number_novalue);
			 free(temp1);
			 free(temp2);		
			 initialize_doublevector(cnet->Vsup_cum, 0.0);
			 
			 
			 //contribution in volume to the channel by subsurface flow
			 for(r=1;r<=Nr;r++){
			 for(c=1;c<=Nc;c++){
			 if( (long)land->LC->co[r][c]!=number_novalue){
			 if (cnet->ch->co[r][c]!=0) {
			 M->co[r][c] = cnet->Vsub_cum->co[cnet->ch->co[r][c]];
			 }else {
			 M->co[r][c] = number_novalue;
			 }
			 }
			 }
			 }
			 
			 temp1 = join_strings(files[fhsup],"Vsub_to_channel");
			 temp2 = join_strings(temp1, SSSS);
			 write_map(temp2, 0, par->format_out, M, UV, number_novalue);
			 free(temp1);
			 free(temp2);		
			 initialize_doublevector(cnet->Vsub_cum, 0.0);*/
		}		
	}
	
	//snow properties
	if(par->output_snow>0 && fmod(times->time+par->Dt,par->output_snow*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt+par->delay_day_recover)/(par->output_snow*3600.0));	
		write_suffix(SSSS, n_file, 0);
		
		if(strcmp(files[fsnowdepth] , string_novalue) != 0){
			for(r=1;r<=snow->S->Dzl->nrh;r++){
				for(c=1;c<=snow->S->Dzl->nch;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=0.0;
						for(l=1;l<=snow->S->lnum->co[r][c];l++){
							M->co[r][c]+=snow->S->Dzl->co[l][r][c];
						}
					}
				}
			}
			
			temp1=join_strings(files[fsnowdepth],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);
			free(temp1);
		}
		
		if(strcmp(files[fsnowmelt] , string_novalue) != 0){
			temp1=join_strings(files[fsnowmelt],SSSS);
			write_map(temp1, 0, par->format_out, snow->MELTED, UV, number_novalue);		
			initmatrix(0.0, snow->MELTED, land->LC, number_novalue);	
			free(temp1);
		}
		
		if(strcmp(files[fsnowsubl] , string_novalue) != 0){
			temp1=join_strings(files[fsnowsubl],SSSS);
			write_map(temp1, 0, par->format_out, snow->SUBL, UV, number_novalue);	
			initmatrix(0.0, snow->SUBL, land->LC, number_novalue);	
			free(temp1);
		}
		
		if(strcmp(files[fswe] , string_novalue) != 0){
			for(r=1;r<=snow->S->Dzl->nrh;r++){
				for(c=1;c<=snow->S->Dzl->nch;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=0.0;
						for(l=1;l<=snow->S->lnum->co[r][c];l++){
							M->co[r][c]+=(snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
						}
					}
				}
			}
			temp1=join_strings(files[fswe],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);
			free(temp1);
			
			for(r=1;r<=snow->S->Dzl->nrh;r++){
				for(c=1;c<=snow->S->Dzl->nch;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=0.0;
						D=0.0;
						for(l=1;l<=snow->S->lnum->co[r][c];l++){
							M->co[r][c]+=(snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
							D+=snow->S->Dzl->co[l][r][c];
						}
						M->co[r][c]/=(1.E-3*D);
					}
				}
			}
			temp1=join_strings(files[fswe], "DENSITY");
			temp2=join_strings(temp1,SSSS);
			write_map(temp2, 0, par->format_out, M, UV, number_novalue);
			free(temp2);
			free(temp1);
			
			if(par->blowing_snow==1){
				temp1=join_strings(files[fswe],"WindTrans");
				temp2=join_strings(temp1, SSSS);
				write_map(temp2, 0, par->format_out, snow->Wtrans_plot, UV, number_novalue);
				initmatrix(0.0, snow->Wtrans_plot, land->LC, number_novalue);
				free(temp2);
				free(temp1);
				
				temp1=join_strings(files[fswe],"WindSubl");
				temp2=join_strings(temp1, SSSS);
				write_map(temp2, 0, par->format_out, snow->Wsubl_plot, UV, number_novalue);
				initmatrix(0.0, snow->Wsubl_plot, land->LC, number_novalue);
				free(temp2);
				free(temp1);
				
				/*temp1=join_strings(files[fswe],"Qtrans_eq");
				 temp2=join_strings(temp1, SSSS);
				 write_map(temp2, 0, par->format_out, snow->Qtrans_eq_plot, UV, number_novalue);
				 initmatrix(0.0, snow->Qtrans_eq_plot, land->LC, number_novalue);
				 free(temp2);
				 free(temp1);
				 
				 temp1=join_strings(files[fswe],"Qtrans");
				 temp2=join_strings(temp1, SSSS);
				 write_map(temp2, 0, par->format_out, snow->Qtrans_plot, UV, number_novalue);
				 initmatrix(0.0, snow->Qtrans_plot, land->LC, number_novalue);
				 free(temp2);
				 free(temp1);
				 
				 temp1=join_strings(files[fswe],"Qsubl_eq");
				 temp2=join_strings(temp1, SSSS);
				 write_map(temp2, 0, par->format_out, snow->Qsub_eq_plot, UV, number_novalue);
				 initmatrix(0.0, snow->Qsub_eq_plot, land->LC, number_novalue);
				 free(temp2);
				 free(temp1);
				 
				 temp1=join_strings(files[fswe],"Qsubl");
				 temp2=join_strings(temp1, SSSS);
				 write_map(temp2, 0, par->format_out, snow->Qsub_plot, UV, number_novalue);
				 initmatrix(0.0, snow->Qsub_plot, land->LC, number_novalue);
				 free(temp2);
				 free(temp1);*/
			}
		}
		
		if(strcmp(files[fsndur] , string_novalue) != 0){
			temp1=join_strings(files[fsndur],SSSS);
			write_map(temp1, 0, par->format_out, snow->t_snow, UV, number_novalue);	
			initmatrix(0.0, snow->t_snow, land->LC, number_novalue);			
			free(temp1);
		}
		
	}
	
	//glacier properties
	if(par->glaclayer_max>0 && par->output_glac>0 && fmod(times->time+par->Dt,par->output_glac*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt+par->delay_day_recover)/(par->output_glac*3600.0));	
		write_suffix(SSSS, n_file, 0);
		
		if(strcmp(files[fglacdepth] , string_novalue) != 0){
			
			for(r=1;r<=glac->G->Dzl->nrh;r++){
				for(c=1;c<=glac->G->Dzl->nch;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=0.0;
						for(l=1;l<=glac->G->lnum->co[r][c];l++){
							M->co[r][c]+=(glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
						}
					}
				}
			}
			temp1=join_strings(files[fglacdepth],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);
			free(temp1);
		}
		
		if(strcmp(files[fglacmelt] , string_novalue) != 0){
			temp1=join_strings(files[fglacmelt],SSSS);
			write_map(temp1, 0, par->format_out, glac->MELTED, UV, number_novalue);		
			initmatrix(0.0, glac->MELTED, land->LC, number_novalue);	
			free(temp1);
		}
		
		if(strcmp(files[fglacsubl] , string_novalue) != 0){
			temp1=join_strings(files[fglacsubl],SSSS);
			write_map(temp1, 0, par->format_out, glac->SUBL, UV, number_novalue);
			initmatrix(0.0, glac->SUBL, land->LC, number_novalue);				
			free(temp1);
		}
		
		if(strcmp(files[fgwe] , string_novalue) != 0){
			for(r=1;r<=glac->G->Dzl->nrh;r++){
				for(c=1;c<=glac->G->Dzl->nch;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=0.0;
						//D=0.0;
						for(l=1;l<=glac->G->lnum->co[r][c];l++){
							M->co[r][c]+=(glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
							//D+=0.001*glac->Dzl->co[l][r][c];
						}
						//if(D>0) M->co[r][c]/=D;
					}else{
						M->co[r][c]=(double)number_novalue;
					}
				}
			}
			temp1=join_strings(files[fgwe],SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);	
			free(temp1);
		}
	}
	
	//SURFACE ENERGY BALANCE	
	//RADIATION
	if(par->output_surfenergy>0 && fmod(times->time+par->Dt,par->output_surfenergy*3600.0)<1.E-5){
		
		n_file=(long)((times->time+par->Dt+par->delay_day_recover)/(par->output_surfenergy*3600.0));	
		write_suffix(SSSS, n_file, 0);
		
		if(strcmp(files[frad] , string_novalue) != 0){
			
			name=join_strings(files[frad],"net");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->Rn_mean, UV, number_novalue);
			initmatrix(0.0, egy->Rn_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			name=join_strings(files[frad],"LWin");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->LWin_mean, UV, number_novalue);
			initmatrix(0.0, egy->LWin_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			name=join_strings(files[frad],"LW");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->LW_mean, UV, number_novalue);
			initmatrix(0.0, egy->LW_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			name=join_strings(files[frad],"SW");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->SW_mean, UV, number_novalue);
			initmatrix(0.0, egy->SW_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			name=join_strings(files[frad],"SWin");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, egy->Rswdown_mean, UV, number_novalue);
			initmatrix(0.0, egy->Rswdown_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			name=join_strings(files[frad],"SWin_beam");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, egy->Rswbeam_mean, UV, number_novalue);	
			initmatrix(0.0, egy->Rswbeam_mean, land->LC, number_novalue);			
			free(temp1);
			free(name);
			
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						if(egy->nDt_sun->co[r][c]>0){
							M->co[r][c]=egy->nDt_shadow->co[r][c]/(double)(egy->nDt_sun->co[r][c]);
						}else{
							M->co[r][c]=0.0;
						}
					}
				}
			}
			name=join_strings(files[frad],"shadowfractime");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, M, UV, number_novalue);
			initlongmatrix(0, egy->nDt_shadow, land->LC, number_novalue);	
			initlongmatrix(0, egy->nDt_sun, land->LC, number_novalue);		
			free(temp1);
			free(name);
			
			
			if(par->distr_stat==1){
				name=join_strings(files[frad],"nmax");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->Rn_max, UV, number_novalue);
				initmatrix(-1.0E+9, egy->Rn_max, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[frad],"nmin");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->Rn_min, UV, number_novalue);
				initmatrix(1.0E+9, egy->Rn_min, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[frad],"LWmax");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->LW_max, UV, number_novalue);
				initmatrix(-1.0E+9, egy->LW_max, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[frad],"LWmin");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->LW_min, UV, number_novalue);
				initmatrix(1.0E+9, egy->LW_min, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[frad],"SWmax");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->SW_max, UV, number_novalue);
				initmatrix(-1.0E+9, egy->SW_max, land->LC, number_novalue);	
				free(temp1);
				free(name);
				
				name=join_strings(files[frad],"SWinmax");
				temp1=join_strings(name,SSSS);
				write_map(temp1, 0, par->format_out, egy->Rswdown_max, UV, number_novalue);
				initmatrix(-1.0E+9, egy->Rswdown_max, land->LC, number_novalue);	
				free(temp1);
				free(name);
			}
		}
		
		
		//GROUND HEAT FLUX
		if(strcmp(files[fG] , string_novalue) != 0){
			
			name=join_strings(files[fG],"mean");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->SEB_mean, UV, number_novalue);
			initmatrix(0.0, egy->SEB_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			/*name=join_strings(files[fG],"snowsoil");
			 temp1=join_strings(name, SSSS);
			 write_map(temp1, 0, par->format_out, egy->G_snowsoil, UV, number_novalue);
			 initmatrix(0.0, egy->G_snowsoil, land->LC, number_novalue);
			 free(temp1);
			 free(name);*/
			
			if(par->distr_stat==1){
				name=join_strings(files[fG],"max");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->G_max, UV, number_novalue);
				initmatrix(-1.0E+9, egy->G_max, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[fG],"min");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->G_min, UV, number_novalue);
				initmatrix(1.0E+9, egy->G_min, land->LC, number_novalue);	
				free(temp1);
				free(name);
			}
		}
		
		
		
		//SENSIBLE HEAT FLUX
		if(strcmp(files[fH] , string_novalue) != 0){
			
			name=join_strings(files[fH],"mean");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->H_mean, UV, number_novalue);
			initmatrix(0.0, egy->H_mean, land->LC, number_novalue);	
			free(temp1);
			free(name);
			
			if(par->distr_stat==1){
				name=join_strings(files[fH],"max");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->H_max, UV, number_novalue);
				initmatrix(-1.0E+9, egy->H_max, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[fH],"min");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->H_min, UV, number_novalue);
				initmatrix(1.0E+9, egy->H_min, land->LC, number_novalue);
				free(temp1);
				free(name);
			}
		}
		
		
		//LATENT HEAT FLUX
		if(strcmp(files[fLE] , string_novalue) != 0){
			
			name=join_strings(files[fLE],"mean");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->ET_mean, UV, number_novalue);
			initmatrix(0.0, egy->ET_mean, land->LC, number_novalue);	
			free(temp1);
			free(name);
			
			if(par->distr_stat==1){			
				name=join_strings(files[fLE],"max");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->ET_max, UV, number_novalue);
				initmatrix(-1.0E+9, egy->ET_max, land->LC, number_novalue);			
				free(temp1);
				free(name);
				
				name=join_strings(files[fLE],"min");	
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->ET_min, UV, number_novalue);	
				initmatrix(1.0E+9, egy->ET_min, land->LC, number_novalue);	
				free(temp1);
				free(name);
			}
		}
		
		//SURFACE TEMPERATURE
		if(strcmp(files[fTs] , string_novalue) != 0){
			name=join_strings(files[fTs],"mean");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, egy->Ts_mean, UV, number_novalue);		
			initmatrix(0.0, egy->Ts_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			if(par->distr_stat==1){
				name=join_strings(files[fTs],"max");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->Ts_max, UV, number_novalue);	
				initmatrix(-1.0E+9, egy->Ts_max, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[fTs],"min");
				temp1=join_strings(name, SSSS);
				write_map(temp1, 0, par->format_out, egy->Ts_min, UV, number_novalue);
				initmatrix(1.0E+9, egy->Ts_min, land->LC, number_novalue);
				free(temp1);
				free(name);
			}
		}
		
	}
	
	//vegetation variables	
	if(par->output_vegetation>0 && fmod(times->time+par->Dt,par->output_vegetation*3600.0)<1.E-5){
		n_file=(long)((times->time+par->Dt+par->delay_day_recover)/(par->output_vegetation*3600.0));	
		write_suffix(SSSS, n_file, 0);
		
		//INTERCEPTED PRECIPITATION
		if(strcmp(files[fcint] , string_novalue) != 0){
			
			temp1=join_strings(files[fcint],"water");
			temp2=join_strings(temp1,SSSS);
			write_map(temp2, 0, par->format_out, wat->wcan_rain, UV, number_novalue);
			free(temp1);
			free(temp2);
			
			temp1=join_strings(files[fcint],"snow");
			temp2=join_strings(temp1,SSSS);
			write_map(temp2, 0, par->format_out, wat->wcan_snow, UV, number_novalue);
			free(temp1);
			free(temp2);
		}
	}
		
	//METEO
	if(par->output_meteo>0 && fmod(times->time+par->Dt,par->output_meteo*3600.0)<1.E-5){
		n_file=(long)((times->time+par->Dt+par->delay_day_recover)/(par->output_meteo*3600.0));	
		write_suffix(SSSS, n_file, 0);
		
		//AIR TEMPERATURE
		if(strcmp(files[fTa] , string_novalue) != 0){
			name=join_strings(files[fTa],"mean");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, met->Ta_mean, UV, number_novalue);
			initmatrix(0.0, met->Ta_mean, land->LC, number_novalue);	
			free(temp1);
			free(name);
			
			if(par->distr_stat==1){
				name=join_strings(files[fTa],"max");
				temp1=join_strings(name,SSSS);
				write_map(temp1, 0, par->format_out, met->Ta_max, UV, number_novalue);
				initmatrix(-1.0E+9, met->Ta_max, land->LC, number_novalue);
				free(temp1);
				free(name);
				
				name=join_strings(files[fTa],"min");
				temp1=join_strings(name,SSSS);
				write_map(temp1, 0, par->format_out, met->Ta_min, UV, number_novalue);	
				initmatrix(1.0E+9, met->Ta_min, land->LC, number_novalue);	
				free(temp1);
				free(name);
			}
		}
		
		//PRECIPITATION
		if(strcmp(files[fprec] , string_novalue) != 0){
			
			name=join_strings(files[fprec],"TOTAL");
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, wat->PrTOT_mean, UV, number_novalue);	
			initmatrix(0.0, wat->PrTOT_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
			
			name=join_strings(files[fprec],"SNOW");	
			temp1=join_strings(name, SSSS);
			write_map(temp1, 0, par->format_out, wat->PrSNW_mean, UV, number_novalue);	
			initmatrix(0.0, wat->PrSNW_mean, land->LC, number_novalue);
			free(temp1);
			free(name);
		}
		
		if(strcmp(files[fwspd] , string_novalue) != 0){
			name=join_strings(files[fwspd],"mean");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, met->Vspdmean, UV, number_novalue);
			initmatrix(0.0, met->Vspdmean, land->LC, number_novalue);	
			free(temp1);
			free(name);
		}
		
		if(strcmp(files[fwdir] , string_novalue) != 0){
			name=join_strings(files[fwdir],"mean");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, met->Vdirmean, UV, number_novalue);	
			initmatrix(0.0, met->Vdirmean, land->LC, number_novalue);	
			free(temp1);
			free(name);
		}
		
		if(strcmp(files[frh] , string_novalue) != 0){
			name=join_strings(files[frh],"mean");
			temp1=join_strings(name,SSSS);
			write_map(temp1, 0, par->format_out, met->RHmean, UV, number_novalue);	
			initmatrix(0.0, met->RHmean, land->LC, number_novalue);
			free(temp1);
			free(name);
		}
	}	
	
	free_doublematrix(M);
	
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	//SPECIAL PLOTS AT SOME DAYS
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	
	if(times->JD_plots->nh > 1 && times->iplot<=times->JD_plots->nh){
		
		i=times->iplot;
		j=2*i-1;
		if( fabs(par->init_date->co[i_sim]+(times->time+par->Dt)/86400. - times->JD_plots->co[j+1]) < 1.E-5 ){
			
			fprintf(flog,"Printing plot number %ld \n",i);
			printf("Printing plot number %ld \n",i);
			
			M=new_doublematrix(Nr,Nc);
			
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=egy->Hgplot->co[r][c] + egy->Hvplot->co[r][c];
					}
				}
			}
			if(strcmp(files[pH] , string_novalue) != 0) plot(files[pH], i, land->LC, M, par->format_out);
			
			
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=egy->LEgplot->co[r][c] + egy->LEvplot->co[r][c];
					}
				}
			}
			if(strcmp(files[pLE] , string_novalue) != 0) plot(files[pLE], i, land->LC, M, par->format_out);
			
			
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						l=1;
						M->co[r][c]=sl->th->co[l][r][c];
					}
				}
			}
			if(strcmp(files[pth] , string_novalue) != 0) plot(files[pth], i, land->LC, M, par->format_out);
			
			
			
			if(strcmp(files[pHg] , string_novalue) != 0) plot(files[pHg], i, land->LC, egy->Hgplot, par->format_out);
			if(strcmp(files[pLEg] , string_novalue) != 0) plot(files[pLEg], i, land->LC, egy->LEgplot, par->format_out);
			if(strcmp(files[pHv] , string_novalue) != 0) plot(files[pHv], i, land->LC, egy->Hvplot, par->format_out);
			if(strcmp(files[pLEv] , string_novalue) != 0) plot(files[pLEv], i, land->LC, egy->LEvplot, par->format_out);
			if(strcmp(files[pSWin] , string_novalue) != 0) plot(files[pSWin], i, land->LC, egy->SWinplot, par->format_out);
			if(strcmp(files[pSWg] , string_novalue) != 0) plot(files[pSWg], i, land->LC, egy->SWgplot, par->format_out);
			if(strcmp(files[pSWv] , string_novalue) != 0) plot(files[pSWv], i, land->LC, egy->SWvplot, par->format_out);
			if(strcmp(files[pLWin] , string_novalue) != 0) plot(files[pLWin], i, land->LC, egy->LWinplot, par->format_out);
			if(strcmp(files[pLWg] , string_novalue) != 0) plot(files[pLWg], i, land->LC, egy->LWgplot, par->format_out);
			if(strcmp(files[pLWv] , string_novalue) != 0) plot(files[pLWv], i, land->LC, egy->LWvplot, par->format_out);
			if(strcmp(files[pTs] , string_novalue) != 0) plot(files[pTs], i, land->LC, egy->Tsplot, par->format_out);
			if(strcmp(files[pTg] , string_novalue) != 0) plot(files[pTg], i, land->LC, egy->Tgplot, par->format_out);
			if(strcmp(files[pTv] , string_novalue) != 0) plot(files[pTv], i, land->LC, egy->Tvplot, par->format_out);
			if(strcmp(files[pTa] , string_novalue) != 0) plot(files[pTa], i, land->LC, met->Taplot, par->format_out);
			if(strcmp(files[pD] , string_novalue) != 0) plot(files[pD], i, land->LC, snow->Dplot, par->format_out);
			if(strcmp(files[pVspd] , string_novalue) != 0) plot(files[pVspd], i, land->LC, met->Vspdplot, par->format_out);
			if(strcmp(files[pVdir] , string_novalue) != 0) plot(files[pVdir], i, land->LC, met->Vdirplot, par->format_out);
			if(strcmp(files[pRH] , string_novalue) != 0) plot(files[pRH], i, land->LC, met->RHplot, par->format_out);
						
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=egy->SWgplot->co[r][c]+egy->LWgplot->co[r][c]-egy->Hgplot->co[r][c]-egy->LEgplot->co[r][c];
					}
				}
			}
			if(strcmp(files[pG] , string_novalue) != 0) plot(files[pG], i, land->LC, M, par->format_out);
			
			
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=find_watertabledepth(r, c, sl);				
					}
				}
			}
			temp1=join_strings(WORKING_DIRECTORY,"plots/watertable");
			plot(temp1, i, land->LC, M, par->format_out);
			free(temp1);
			
			
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=find_activelayerdepth(r, c, sl);				
					}
				}
			}
			temp1=join_strings(WORKING_DIRECTORY,"plots/thaweddepth");
			plot(temp1, i, land->LC, M, par->format_out);
			free(temp1);
			
			
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						M->co[r][c]=egy->Hgplot->co[r][c] + egy->LEgplot->co[r][c];				
					}
				}
			}
			temp1=join_strings(WORKING_DIRECTORY,"plots/LHE");
			plot(temp1, i, land->LC, M, par->format_out);
			free(temp1);
			
			
			initialize_doublematrix(M,(double)number_novalue);
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						if(egy->SWgplot->co[r][c]>0){
							M->co[r][c]=(egy->Hgplot->co[r][c] + egy->LEgplot->co[r][c])/egy->SWgplot->co[r][c];
						}else{
							M->co[r][c]=0.;
						}
					}
				}
			}
			temp1=join_strings(WORKING_DIRECTORY,"plots/LHE_to_SW");
			plot(temp1, i, land->LC, M, par->format_out);
			free(temp1);
			
			
			
			
			
			free_doublematrix(M);
			
			
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if( (long)land->LC->co[r][c]!=number_novalue){
						egy->Hgplot->co[r][c]=0.0;
						egy->LEgplot->co[r][c]=0.0;
						egy->Hvplot->co[r][c]=0.0;
						egy->LEvplot->co[r][c]=0.0;
						egy->SWinplot->co[r][c]=0.0;
						egy->SWgplot->co[r][c]=0.0;
						egy->SWvplot->co[r][c]=0.0;
						egy->LWinplot->co[r][c]=0.0;
						egy->LWgplot->co[r][c]=0.0;
						egy->LWvplot->co[r][c]=0.0;
						egy->Tsplot->co[r][c]=0.0;
						egy->Tgplot->co[r][c]=0.0;
						egy->Tvplot->co[r][c]=0.0;	
						
						snow->Dplot->co[r][c]=0.0;
						met->Taplot->co[r][c]=0.0;
						met->Vspdplot->co[r][c]=0.0;
						met->Vdirplot->co[r][c]=0.0;
						met->RHplot->co[r][c]=0.0;
					}
				}
			}
			
			times->iplot ++;
		}
	}
	
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	//SAVING POINTS
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	
	if(times->time==0) isavings=0;
	
	if(isavings < par->saving_points->nh){
		
		if(par->saving_points->nh==1 && par->saving_points->co[1]==0.0){
			
			isavings=1;
			
		}else{
			
			if(times->time+par->Dt >= 86400.0*par->saving_points->co[isavings+1]){
				isavings+=1;
				
				printf("Writing recovering files, saving point number %ld\n",isavings);
				fprintf(flog, "Writing recovering files, saving point number %ld\n",isavings);
				
				write_suffix(SSSS, isavings+par->recover, 0);
				
				for(l=0;l<=Nl;l++){
					if(strcmp(files[rpsi] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rpsi], 0, par->format_out, sl->P, UV, number_novalue);
					if(l>0){
						if(strcmp(files[riceg] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[riceg], 0, par->format_out, sl->thice, UV, number_novalue);
						if(strcmp(files[rTg] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rTg], 0, par->format_out, sl->T, UV, number_novalue);
					}
				}
				
				if(strcmp(files[rwcrn] , string_novalue) != 0) write_map(join_strings(files[rwcrn],SSSS), 0, par->format_out, wat->wcan_rain, UV, number_novalue);			
				if(strcmp(files[rwcsn] , string_novalue) != 0) write_map(join_strings(files[rwcsn],SSSS), 0, par->format_out, wat->wcan_snow, UV, number_novalue);			
				
				for (r=1; r<=Nr; r++) {
					for (c=1; c<=Nc; c++) {
						if ((long)land->LC->co[r][c] != number_novalue) {
							if ((long)sl->Tv->co[r][c] == number_novalue) sl->Tv->co[r][c] = 0.;
						}
					}
				}
				
				if(strcmp(files[rTv] , string_novalue) != 0) write_map(join_strings(files[rTv],SSSS), 0, par->format_out, sl->Tv, UV, number_novalue);			
				
				for(l=1;l<=par->snowlayer_max;l++){
					if(strcmp(files[rDzs] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rDzs], 0, par->format_out, snow->S->Dzl, UV, number_novalue);
					if(strcmp(files[rwls] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rwls], 0, par->format_out, snow->S->w_liq, UV, number_novalue);
					if(strcmp(files[rwis] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rwis], 0, par->format_out, snow->S->w_ice, UV, number_novalue);
					if(strcmp(files[rTs] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rTs], 0, par->format_out, snow->S->T, UV, number_novalue);
				}
				if(strcmp(files[rsnag_nondim] , string_novalue) != 0) write_map(join_strings(files[rsnag_nondim],SSSS), 0, par->format_out, snow->nondimens_age, UV, number_novalue);			
				if(strcmp(files[rsnag_dim] , string_novalue) != 0) write_map(join_strings(files[rsnag_dim],SSSS), 0, par->format_out, snow->dimens_age, UV, number_novalue);			
				
				if(strcmp(files[rns] , string_novalue) != 0){
					M=copydouble_longmatrix(snow->S->lnum);
					write_map(join_strings(files[rns], SSSS), 1, par->format_out, M, UV, number_novalue);
					free_doublematrix(M);
				}
				
				if(par->glaclayer_max>0){
					for(l=1;l<=par->glaclayer_max;l++){
						if(strcmp(files[rDzi] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rDzi], 0, par->format_out, glac->G->Dzl, UV, number_novalue);
						if(strcmp(files[rwli] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rwli], 0, par->format_out, glac->G->w_liq, UV, number_novalue);
						if(strcmp(files[rwii] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rwii], 0, par->format_out, glac->G->w_ice, UV, number_novalue);
						if(strcmp(files[rTi] , string_novalue) != 0) write_tensorseries(1, l, isavings+par->recover, files[rTi], 0, par->format_out, glac->G->T, UV, number_novalue);
					}
					
					if(strcmp(files[rni] , string_novalue) != 0){
						M=copydouble_longmatrix(glac->G->lnum);
						write_map(join_strings(files[rni],SSSS), 1, par->format_out, M, UV, number_novalue);
						free_doublematrix(M);				
					}
				}		
				

				if(strcmp(files[rpsich] , string_novalue) != 0){
					T=new_doubletensor0(Nl, Nr, Nc);
					initialize_doubletensor(T, 0.0);
					for (l=0; l<=Nl; l++) {
						for (i=1; i<=par->total_channel; i++){
							T->co[l][cnet->r->co[i]][cnet->c->co[i]] = cnet->P->co[l][i];
						}
						write_tensorseries(1, l, isavings+par->recover, files[rpsich], 0, par->format_out, T, UV, number_novalue);
					}
					free_doubletensor(T);
				}

				if(strcmp(files[rTgch] , string_novalue) != 0){
					T=new_doubletensor(Nl, Nr, Nc);
					initialize_doubletensor(T, 0.0);
					for (l=1; l<=Nl; l++) {
						for (i=1; i<=par->total_channel; i++){
							T->co[l][cnet->r->co[i]][cnet->c->co[i]] = cnet->T->co[l][i];
						}
						write_tensorseries(1, l, isavings+par->recover, files[rTgch], 0, par->format_out, T, UV, number_novalue);
					}
					free_doubletensor(T);
				}
				
				if(strcmp(files[ricegch] , string_novalue) != 0){
					T=new_doubletensor(Nl, Nr, Nc);
					initialize_doubletensor(T, 0.0);
					for (l=1; l<=Nl; l++) {
						for (i=1; i<=par->total_channel; i++){
							T->co[l][cnet->r->co[i]][cnet->c->co[i]] = cnet->thice->co[l][i];
						}
						write_tensorseries(1, l, isavings+par->recover, files[ricegch], 0, par->format_out, T, UV, number_novalue);
					}
					free_doubletensor(T);
				}
				
			}
		}
		
	}
	
	
	fclose(flog);
	
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac){

	/*internal auxiliary variables:*/
	long i,l,m,j,r,c;
	char *name,*temp,*temp2,SSSS[ ]={"SSSS"};
	long sy;
	short lu, first_column;
	DOUBLEVECTOR *root_fraction;
	FILE *f; 
			
	//DISCHARGE
	if (par->state_discharge == 1 && strcmp(files[fQ] , string_novalue) != 0){
		temp=join_strings(files[fQ],textfile);
		f=t_fopen(temp,"w");
		//fprintf(f,"DATE[day/month/year hour:min],JDfrom0,JD,Q_tot[mc/s],Qsub_ch[mc/s],Qsup_ch[mc/s]\n");
		fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,Qtot[m3/s],Vsup/Dt[m3/s],Vsub/Dt[m3/s],Vchannel[m3],Qoutland[m3/s]\n");
		t_fclose(f);
		free(temp);
	}
	
	
	if(par->state_pixel == 1){
		
		//output matrix and vectors
		m = (long)otot;
		outdata_point=(double**)malloc(m*sizeof(double*));
		for (i=0; i<m; i++) {
			outdata_point[i]=(double*)malloc(par->rc->nrh*sizeof(double));
			for (j=0; j<par->rc->nrh; j++) {
				outdata_point[i][j] = 0.;
			}
		}
				
		if(strcmp(files[fpointwriteend] , string_novalue) != 0){
			name=join_strings(files[fpointwriteend],textfile);
			ffpoint=fopen(name,"w");
			first_column=1;
			for(j=0;j<noutputpoint;j++){
				if(first_column==0){
					fprintf(ffpoint,",");
				}else {
					first_column = 0;
				}					
				if (outputpoint[j] >= 0) {
					fprintf(ffpoint,"%s",headerpoint[outputpoint[j]]);					
				}else {
					fprintf(ffpoint, "None");
				}
			}
			fprintf(ffpoint,"\n");
			free(name);
		}
		
		if(strcmp(files[fsnzwriteend] , string_novalue) != 0){
			name=join_strings(files[fsnzwriteend],textfile);
			ffsnow=fopen(name,"w");
			first_column=1;
			for(j=0;j<noutputsnow;j++){
				if(first_column==0){
					fprintf(ffsnow,",");
				}else {
					first_column = 0;
				}			
				if (outputsnow[j] >= 0 && outputsnow[j]<=5) {
					fprintf(ffsnow,"%s",headersnow[outputsnow[j]]);					
				}else if (outputsnow[j] >= 6) {
					l = (long)fmod( (double)outputsnow[j]-6., (double)par->snowlayer_max ) + 1;
					n = floor( ( (double)outputsnow[j]-6.) / (double)par->snowlayer_max ) + 6;
					fprintf(ffsnow, "%s(%ld)",headersnow[n],l);
				}else{
					fprintf(ffsnow, "None");
				}
			}
			fprintf(ffsnow,"\n");
			free(name);
		}
		
		if(par->glaclayer_max>0 && strcmp(files[fglzwriteend] , string_novalue) != 0){
			name=join_strings(files[fglzwriteend],textfile);
			ffglac=fopen(name,"w");
			first_column=1;
			for(j=0;j<noutputglac;j++){
				if(first_column==0){
					fprintf(ffglac,",");
				}else {
					first_column = 0;
				}					
				if (outputglac[j] >= 0 && outputglac[j]<=5) {
					fprintf(ffglac,"%s",headerglac[outputglac[j]]);					
				}else if (outputglac[j] >= 6) {
					l = (long)fmod( (double)outputglac[j]-6., (double)par->glaclayer_max ) + 1;
					n = floor( ( (double)outputglac[j]-6.) / (double)par->glaclayer_max ) + 6;
					fprintf(ffglac, "%s(%ld)",headerglac[n],l);
				}else
					fprintf(ffglac, "None");
			}
			fprintf(ffglac,"\n");
			free(name);
		}
		
		if(strcmp(files[fTzwriteend] , string_novalue) != 0){
			name=join_strings(files[fTzwriteend],textfile);
			ffT=fopen(name,"w");
			write_soil_header(ffT, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fTzavwriteend] , string_novalue) != 0){
			name=join_strings(files[fTzavwriteend],textfile);
			ffTav=fopen(name,"w");		
			write_soil_header(ffTav, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fpsiztotwriteend] , string_novalue) != 0){
			name=join_strings(files[fpsiztotwriteend],textfile);
			ffpsitot=fopen(name,"w");		
			write_soil_header(ffpsitot, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fpsizwriteend] , string_novalue) != 0){
			name=join_strings(files[fpsizwriteend],textfile);
			ffpsi=fopen(name,"w");		
			write_soil_header(ffpsi, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[fliqzwriteend] , string_novalue) != 0){
			name=join_strings(files[fliqzwriteend],textfile);
			ffliq=fopen(name,"w");		
			write_soil_header(ffliq, sl->pa->co[1][jdz]);
			free(name);
		}

		if(strcmp(files[fliqzavwriteend] , string_novalue) != 0){
			name=join_strings(files[fliqzavwriteend],textfile);
			ffliqav=fopen(name,"w");		
			write_soil_header(ffliqav, sl->pa->co[1][jdz]);
			free(name);
		}

		if(strcmp(files[ficezwriteend] , string_novalue) != 0){
			name=join_strings(files[ficezwriteend],textfile);
			ffice=fopen(name,"w");				
			write_soil_header(ffice, sl->pa->co[1][jdz]);
			free(name);
		}
		
		if(strcmp(files[ficezavwriteend] , string_novalue) != 0){
			name=join_strings(files[ficezavwriteend],textfile);
			fficeav=fopen(name,"w");		
			write_soil_header(fficeav, sl->pa->co[1][jdz]);
			free(name);
		}
		
		root_fraction=new_doublevector(Nl);
   
		
		//DATA POINTS
		for(i=1;i<=par->rc->nrh;i++){
						
			write_suffix(SSSS, par->IDpoint->co[i], 0);
			r=par->rc->co[i][1];
			c=par->rc->co[i][2];				
			sy=sl->type->co[r][c];
			lu=(short)land->LC->co[r][c];
			
			if(strcmp(files[fpoint] , string_novalue) != 0 && par->point_sim != 1){			
				name=join_strings(files[fpoint],"_info_");
				temp=join_strings(name,SSSS);
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

				for(l=1;l<=Nl;l++){
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
				temp=join_strings(files[fpoint],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");
				first_column=1;
				for(j=0;j<noutputpoint;j++){
					if(first_column==0){
						fprintf(f,",");
					}else {
						first_column = 0;
					}					
					if (outputpoint[j] >= 0) {
						fprintf(f,"%s",headerpoint[outputpoint[j]]);					
					}else {
						fprintf(f, "None");
					}
				}
				fprintf(f,"\n");
				t_fclose(f);
				free(name);
				free(temp);
			}
			
			if(strcmp(files[fsnz] , string_novalue) != 0){
				temp=join_strings(files[fsnz],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");
				first_column=1;
				for(j=0;j<noutputsnow;j++){
					if(first_column==0){
						fprintf(f,",");
					}else {
						first_column = 0;
					}			
					if (outputsnow[j] >= 0 && outputsnow[j]<=5) {
						fprintf(f,"%s",headersnow[outputsnow[j]]);					
					}else if (outputsnow[j] >= 6) {
						l = (long)fmod( (double)outputsnow[j]-6., (double)par->snowlayer_max ) + 1;
						n = floor( ( (double)outputsnow[j]-6.) / (double)par->snowlayer_max ) + 6;
						fprintf(f, "%s(%ld)",headersnow[n],l);
					}else{
						fprintf(f, "None");
					}
				}
				fprintf(f,"\n");
				t_fclose(f);
				free(name);
				free(temp);
			}
			
			if(par->glaclayer_max>0 && strcmp(files[fglz] , string_novalue) != 0){
				temp=join_strings(files[fglz],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");
				first_column=1;
				for(j=0;j<noutputglac;j++){
					if(first_column==0){
						fprintf(f,",");
					}else {
						first_column = 0;
					}					
					if (outputglac[j] >= 0 && outputglac[j]<=5) {
						fprintf(f,"%s",headerglac[outputglac[j]]);					
					}else if (outputglac[j] >= 6) {
						l = (long)fmod( (double)outputglac[j]-6., (double)par->glaclayer_max ) + 1;
						n = floor( ( (double)outputglac[j]-6.) / (double)par->glaclayer_max ) + 6;
						fprintf(f, "%s(%ld)",headerglac[n],l);
					}else
						fprintf(f, "None");
				}
				fprintf(f,"\n");
				t_fclose(f);
				free(name);
				free(temp);				
			}
		
			if(strcmp(files[fTz] , string_novalue) != 0){
				temp=join_strings(files[fTz],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");			
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
				free(temp);
			}
			
			if(strcmp(files[fTzav] , string_novalue) != 0){
				temp=join_strings(files[fTzav],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
				free(temp);
			}
		
			if(strcmp(files[fpsiztot] , string_novalue) != 0){
				temp=join_strings(files[fpsiztot],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);
				free(name);
				free(temp);
			}
			
			if(strcmp(files[fpsiz] , string_novalue) != 0){
				temp=join_strings(files[fpsiz],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);
				free(name);
				free(temp);
			}
		
			if(strcmp(files[fliqz] , string_novalue) != 0){
				temp=join_strings(files[fliqz],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);
				free(name);
				free(temp);	
			}
			
			if(strcmp(files[fliqzav] , string_novalue) != 0){
				temp=join_strings(files[fliqzav],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
				free(temp);
			}
		
			if(strcmp(files[ficez] , string_novalue) != 0){
				temp=join_strings(files[ficez],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");				
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
				free(temp);
			}
							
			if(strcmp(files[ficezav] , string_novalue) != 0){
				temp=join_strings(files[ficezav],SSSS);
				name=join_strings(temp,textfile);
				f=t_fopen(name,"w");		
				write_soil_header(f, sl->pa->co[1][jdz]);
				t_fclose(f);	
				free(name);
				free(temp);
			}
						
		}
		
		free_doublevector(root_fraction);
		
	}
	
	m = (long)ootot;
	outdata_basin=(double*)malloc(m*sizeof(double));
	for (i=0; i<m; i++) {
		outdata_basin[i]=0.;
	}		

	if(par->state_basin == 1){
				
		//DATA BASIN
		if(strcmp(files[fbaswriteend] , string_novalue) != 0){
			temp=join_strings(files[fbaswriteend],textfile);
			ffbas=fopen(temp,"w");
			first_column=1;
			for(j=0;j<noutputbasin;j++){
				if(first_column==0){
					fprintf(ffbas,",");
				}else {
					first_column = 0;
				}					
				if (outputbasin[j] >= 0) {
					fprintf(ffbas,"%s",headerbasin[outputbasin[j]]);					
				}else {
					fprintf(ffbas, "None");
				}
			}
			fprintf(ffbas,"\n");
			free(temp);
		}
		
		if(strcmp(files[fbas] , string_novalue) != 0){
			temp=join_strings(files[fbas],textfile);
			f=t_fopen(temp,"w");
			first_column=1;
			for(j=0;j<noutputbasin;j++){
				if(first_column==0){
					fprintf(f,",");
				}else {
					first_column = 0;
				}					
				if (outputbasin[j] >= 0) {
					fprintf(f,"%s",headerbasin[outputbasin[j]]);					
				}else {
					fprintf(f, "None");
				}
			}
			fprintf(f,"\n");
			
			t_fclose(f);
			free(temp);
		}
				
	}
	
	//SNOW COVERED AREA STATISTICS
	if(par->point_sim!=1 && strcmp(files[fSCA] , string_novalue) != 0){
		temp=join_strings(files[fSCA],textfile);
		f=t_fopen(temp,"w");
		fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,snowDav,SWEav,Tav,Tsav,perc.SFA,perc.SCA\n");
		t_fclose(f);
		free(temp);
	}
	
	

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_output(long i, long iname, double init_date, double JDfrom0, double JD, long day, long month, long year, long hour, long minute, LONGMATRIX *rc, SOIL *sl, double psimin){

	char *name,*temp,SSSS[ ]={"SSSS"};
	long l;
	FILE *f;
	
	write_suffix(SSSS, iname, 0);	
	
	if(strcmp(files[fTz] , string_novalue) != 0){
		temp=join_strings(files[fTz],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");	
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fTzwriteend] , string_novalue) != 0){
		write_soil_file(iname, ffT, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzplot->co[i]);
	}
	
	if(strcmp(files[fTzav] , string_novalue) != 0){
		temp=join_strings(files[fTzav],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");	
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzavplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fTzavwriteend] , string_novalue) != 0){
		write_soil_file(iname, ffTav, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzavplot->co[i]);
	}
	
	if(strcmp(files[fpsiztot] , string_novalue) != 0){
		temp=join_strings(files[fpsiztot],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Ptotzplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fpsiztotwriteend] , string_novalue) != 0){
		write_soil_file(iname, ffpsitot, day, month, year, hour, minute, JDfrom0, init_date, sl->Ptotzplot->co[i]);
	}
	
	if(strcmp(files[fpsiz] , string_novalue) != 0){
		temp=join_strings(files[fpsiz],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Pzplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fpsizwriteend] , string_novalue) != 0){
		write_soil_file(iname, ffpsi, day, month, year, hour, minute, JDfrom0, init_date, sl->Pzplot->co[i]);
	}
	
	if(strcmp(files[fliqz] , string_novalue) != 0){
		temp=join_strings(files[fliqz],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thzplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[fliqzwriteend] , string_novalue) != 0){
		write_soil_file(iname, ffliq, day, month, year, hour, minute, JDfrom0, init_date, sl->thzplot->co[i]);
	}

	if(strcmp(files[fliqzav] , string_novalue) != 0){
		temp=join_strings(files[fliqzav],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thzavplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[fliqzavwriteend] , string_novalue) != 0){
		write_soil_file(iname, ffliqav, day, month, year, hour, minute, JDfrom0, init_date, sl->thzavplot->co[i]);
	}
	
	if(strcmp(files[ficez] , string_novalue) != 0){
		temp=join_strings(files[ficez],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thicezplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}

	if(strcmp(files[ficezwriteend] , string_novalue) != 0){
		write_soil_file(iname, ffice, day, month, year, hour, minute, JDfrom0, init_date, sl->thicezplot->co[i]);
	}

	if(strcmp(files[ficezav] , string_novalue) != 0){
		temp=join_strings(files[ficezav],SSSS);
		name=join_strings(temp,textfile);
		f=fopen(name,"a");		
		write_soil_file(iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thicezavplot->co[i]);
		fclose(f);
		free(name);
		free(temp);
	}
	
	if(strcmp(files[ficezavwriteend] , string_novalue) != 0){
		write_soil_file(iname, fficeav, day, month, year, hour, minute, JDfrom0, init_date, sl->thicezavplot->co[i]);
	}
	
	for(l=1;l<=Nl;l++){
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] = 0.0;
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] = 0.0;
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thicezavplot->co[i][l] = 0.0;
	}
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_file(long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, double *var){

	short first_column=1;
	long j, l;
	
	for (j=0; j<noutputsoil; j++) {
		if(first_column==0){
			fprintf(f,",");
		}else {
			first_column = 0;
		}
		if (outputsoil[j] >= 0) {
			if (outputsoil[j] == 0) {
				fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)d,(float)m,(float)y,(float)h,(float)mi);		
			}else if (outputsoil[j] == 1) {
				fprintf(f, "%f",JDfrom0);
			}else if (outputsoil[j] == 2) {
				fprintf(f, "%f",JDfrom0-JDfrom0init);
			}else if (outputsoil[j] == 3) {
				fprintf(f, "%ld",i_sim);
			}else if (outputsoil[j] == 4) {
				fprintf(f, "%ld",i_run);																				
			}else if (outputsoil[j] == 5) {
				fprintf(f, "%ld",i);
			}
		}else {
			fprintf(f,"%f",(double)number_novalue);
		}
	}
	
	for(l=1;l<=Nl;l++){
		fprintf(f,",%f",var[l]);
	}
	
	fprintf(f," \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void write_soil_header(FILE *f, double *dz){
	
	short first_column=1;
	long j, l;
	double z=0.0;
	
	for(j=0;j<noutputsoil;j++){
		if(first_column==0){
			fprintf(f,",");
		}else {
			first_column = 0;
		}					
		if (outputsoil[j] >= 0 && outputsoil[j]<=5) {
			fprintf(f,"%s",headersoil[outputsoil[j]]);					
		}else{
			fprintf(f, "none");
		}
	}

	for(l=1;l<=Nl;l++){
		z += dz[l];
		fprintf(f,",%.0f ",z-0.5*dz[l]);
	}
	fprintf(f," \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

void plot(char *name, long i_plot, DOUBLEMATRIX *LC, DOUBLEMATRIX *M, short format){
	
	char ADS[ ]={"iiii"};
	char *temp;
	
	write_suffix(ADS, i_plot, 0);
	
	temp=join_strings(name,ADS);
	write_map(temp, 0, format, M, UV, number_novalue);
	free(temp);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

