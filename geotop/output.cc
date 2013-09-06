
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.223 'Wallis' - 26 Jul 2011
 
 Copyright (c), 2011 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.223 'Wallis'
 
 GEOtop 1.223 'Wallis' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.223 'Wallis' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.

 */
    
#include <string>
#include "output.h"
using namespace std;

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//void write_output(TIMES *times, WATER *wat, CHANNEL *cnet, PAR *par, TOPO *top, LAND *land, SOIL *sl, Energy *egy, SNOW *snow, GLACIER *glac, METEO *met)
  void write_output(Times *times, Water *wat, Channel *cnet, Par *par, Topo *top, Land *land, Soil *sl, Energy *egy, Snow *snow, Glacier *glac, Meteo *met)

{
	/*internal auxiliary variables:*/
	long i,j,r,c,l,m; /*counters*/
	long n_file;      /*number of file of the type "TETAxySSSlZZ"(i.e. number of the basin-time-step)*/
	char NNNNN[ ]={"NNNNN"};
	char RRRRR[ ]={"RRRRR"};
	char SSSSS[ ]={"SSSSS"};
	char NNNN[ ]={"NNNN"};
	char rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
//	char *name, *temp1, *temp2, *s1, *s2;

	string name, temp1, temp2 , s1, s2;
	FILE *f, *flog;
	
//	time variables
	time_t stop_time;
	double percent_done, remaining_time, total_time;
	short first_column;
	double JD, JDfrom0;
	long day, month, year, hour, minute;
	
// 	static double Qsub_ch, Qsup_ch;
	static long isavings;
	static double mass_error_tot;
	static double t_discharge, t_basin, t_point, t_rec;
	double Vchannel, Vsub, Vsup;
		
//	other variables
//	DOUBLEVECTOR *V;
    GeoVector<double> V;

//	DOUBLEMATRIX *M;
    GeoMatrix<double> M;

	double D, cosslope;
	
	
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
//	JDfrom0 = convert_tfromstart_JDfrom0(times->time+par->Dt, par->init_date->co[i_sim]);
	JDfrom0 = convert_tfromstart_JDfrom0(times->time+par->Dt, par->init_date[i_sim]);
	convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
	convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);
	
	
	//DISCHARGE
	//****************************************************************************************************************
	//****************************************************************************************************************
	
	//if(par->state_discharge == 1 && par->Dtplot_discharge->co[i_sim] > 1.E-5 && strcmp(files[fQ] , string_novalue) != 0){
	if(par->state_discharge == 1 && par->Dtplot_discharge[i_sim] > 1.E-5 && strcmp(files[fQ] , string_novalue) != 0){
		t_discharge += par->Dt;
		
	//	if (fabs(t_discharge - par->Dtplot_discharge->co[i_sim]) < 1.E-5){
		if (fabs(t_discharge - par->Dtplot_discharge[i_sim]) < 1.E-5){
			Vchannel = 0.;
			Vsub = 0.;
			Vsup = 0.;
			for(l=1;l<=par->total_channel;l++){
			//	r = cnet->r->co[l];
				r = cnet->r[l];
			//	c = cnet->c->co[l];
				c = cnet->c[l];
			//	Vchannel += 1.E-3 * Fmax(cnet->SS->P->co[0][l], 0.) / cos(top->slope->co[r][c]*GTConst::Pi/180.) * UV->U->co[1] * par->w_dx * cnet->length->co[l];
				Vchannel += 1.E-3 * Fmax(cnet->SS->P[0][l], 0.) / cos(top->slope[r][c]*GTConst::Pi/180.) * UV->U[1] * par->w_dx * cnet->length[l];
			//	Vsub += cnet->Vsub->co[l];
				Vsub += cnet->Vsub[l];
			//	Vsup += cnet->Vsup->co[l];
				Vsup += cnet->Vsup[l];
			}
				
			if (par->recover > 0) write_suffix(rec, par->recover, 4);
			if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

			if (par->recover>0) {
				temp1 = join_strings(files[fQ], rec);
			//	name = join_strings(temp1, textfile);
				name = temp1+ textfile;
			//	free(temp1);
			}else if (par->n_ContRecovery>0) {
				temp1 = join_strings(files[fQ], crec);
			//	name = join_strings(temp1, textfile);
				name = temp1+ textfile;
			//	free(temp1);
			}else {
			//	name = join_strings(files[fQ], textfile);
				name = files[fQ] + string(textfile) ;
			}
			
			f=fopen(name.c_str(),"a");

			fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
			fprintf(f,",%f,%f,%f",(times->time+par->Dt)/GTConst::secinday,JDfrom0,JD);
//			fprintf(f,",%e,%e,%e,%e,%e,%e,%e\n",cnet->Vout/(double)par->Dtplot_discharge->co[i_sim],Vsup/(double)par->Dtplot_discharge->co[i_sim],
//					Vsub/(double)par->Dtplot_discharge->co[i_sim],Vchannel,wat->Voutlandsup/(double)par->Dtplot_discharge->co[i_sim],
//					wat->Voutlandsub/(double)par->Dtplot_discharge->co[i_sim],wat->Voutbottom/(double)par->Dtplot_discharge->co[i_sim]);

			fprintf(f,",%e,%e,%e,%e,%e,%e,%e\n",cnet->Vout/(double)par->Dtplot_discharge[i_sim],Vsup/(double)par->Dtplot_discharge[i_sim],
					Vsub/(double)par->Dtplot_discharge[i_sim],Vchannel,wat->Voutlandsup/(double)par->Dtplot_discharge[i_sim],
					wat->Voutlandsub/(double)par->Dtplot_discharge[i_sim],wat->Voutbottom/(double)par->Dtplot_discharge[i_sim]);


			fclose(f);
		//	free(name);
			
			t_discharge = 0.0;
			
		//	initialize_doublevector(cnet->Vsub, 0.);
			cnet->Vsub.resize(cnet->Vsub.size(),0.);
		//	initialize_doublevector(cnet->Vsup, 0.);
			cnet->Vsup.resize(cnet->Vsup.size(),0.);
			
			cnet->Vout = 0.;
			wat->Voutbottom = 0.;
			wat->Voutlandsub = 0.;
			wat->Voutlandsup = 0.;
			
		}
	}
	
	//DATA POINT
	//****************************************************************************************************************
	//****************************************************************************************************************	
	
//	if(par->Dtplot_point->co[i_sim] > 1.E-5){
	if(par->Dtplot_point[i_sim] > 1.E-5){
		
		t_point += par->Dt;
		cum_time += par->Dt;
		
		if(par->state_pixel == 1){
			
		//	for(i=1;i<=par->rc->nrh;i++){
			for(i=1;i<par->rc.getRows();i++){
				
			//	r=par->rc->co[i][1];
				r=par->rc[i][1];
			//	c=par->rc->co[i][2];
				c=par->rc[i][2];
				j=top->j_cont[r][c];
				
				for(l=1;l<=Nl;l++){
				//	if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] += sl->SS->T->co[l][j]*(par->Dt/par->Dtplot_point->co[i_sim]);
					if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot[i][l] += sl->SS->T[l][j]*(par->Dt/par->Dtplot_point[i_sim]);
				//	if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] += sl->th->co[l][j]*(par->Dt/par->Dtplot_point->co[i_sim]);
					if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot[i][l] += sl->th[l][j]*(par->Dt/par->Dtplot_point[i_sim]);
				//	if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot->co[i][l] += sl->SS->thi->co[l][j]*(par->Dt/par->Dtplot_point->co[i_sim]);
					if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot[i][l] += sl->SS->thi[l][j]*(par->Dt/par->Dtplot_point[i_sim]);
				}
				
			//	D = find_activelayerdepth_up(j, sl->type->co[r][c], sl);
				D = find_activelayerdepth_up(j, sl->type[r][c], sl);
			//	odpnt[othawedup][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);
				odpnt[othawedup][i-1] += D * (par->Dt/par->Dtplot_point[i_sim]);
				
			//	D = find_activelayerdepth_dw(j, sl->type->co[r][c], sl);
				D = find_activelayerdepth_dw(j, sl->type[r][c], sl);
			//	odpnt[othaweddw][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);
				odpnt[othaweddw][i-1] += D * (par->Dt/par->Dtplot_point[i_sim]);
				
			//	D = find_watertabledepth_up(j, sl->type->co[r][c], sl);
				D = find_watertabledepth_up(j, sl->type[r][c], sl);
			//	if (D < 1.E-5 && sl->SS->P->co[0][j] > 0) D = -sl->SS->P->co[0][j] / cos(top->slope->co[r][c] * GTConst::Pi/180.);
				if (D < 1.E-5 && sl->SS->P[0][j] > 0) D = -sl->SS->P[0][j] / cos(top->slope[r][c] * GTConst::Pi/180.);
			//	odpnt[owtableup][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);
				odpnt[owtableup][i-1] += D * (par->Dt/par->Dtplot_point[i_sim]);
				
			//	if (sl->SS->P->co[0][j] > 0){
				if (sl->SS->P[0][j] > 0){
				//	D = -sl->SS->P->co[0][i] / cos(top->slope->co[r][c] * GTConst::Pi/180.);
					D = -sl->SS->P[0][i] / cos(top->slope[r][c] * GTConst::Pi/180.);
				}else {
				//	D = find_watertabledepth_dw(i, sl->type->co[r][c], sl);
					D = find_watertabledepth_dw(i, sl->type[r][c], sl);
				}
			//	odpnt[owtabledw][i-1] += D * (par->Dt/par->Dtplot_point->co[i_sim]);
				odpnt[owtabledw][i-1] += D * (par->Dt/par->Dtplot_point[i_sim]);
			}
		}
			
	//	Print of pixel-output every times->n_pixel time step
	//	if (fabs(t_point - par->Dtplot_point->co[i_sim])<1.E-5){
		if (fabs(t_point - par->Dtplot_point[i_sim])<1.E-5){
			
			if(par->state_pixel == 1){
				
				if (par->recover > 0) write_suffix(rec, par->recover, 4);
				if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);
			
			//	for(i=1;i<=par->rc->nrh;i++){
				for(i=1;i<par->rc.getRows();i++){
					
				//	write_suffix(NNNN, par->IDpoint->co[i], 0);
					write_suffix(NNNN, par->IDpoint[i], 0);
				//	r=par->rc->co[i][1];
					r=par->rc[i][1];
				//	c=par->rc->co[i][2];
					c=par->rc[i][2];
					j=top->j_cont[r][c];
					
					if (par->output_vertical_distances == 1) {
					//	cosslope = cos( Fmin(GTConst::max_slope, top->slope->co[r][c]) * GTConst::Pi/180. );
						cosslope = cos( Fmin(GTConst::max_slope, top->slope[r][c]) * GTConst::Pi/180. );
					}else {
						cosslope = 1.;
					}
										
				//	soil data
					for(l=1;l<=Nl;l++){
					//	if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot->co[i][l] = sl->SS->T->co[l][j];
						if(strcmp(files[fTz] , string_novalue) != 0 || strcmp(files[fTzwriteend] , string_novalue) != 0) sl->Tzplot[i][l] = sl->SS->T[l][j];
					//	if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot->co[i][l] = sl->Ptot->co[l][j];
						if(strcmp(files[fpsiztot] , string_novalue) != 0 || strcmp(files[fpsiztotwriteend] , string_novalue) != 0) sl->Ptotzplot[i][l] = sl->Ptot[l][j];
					//	if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot->co[i][l] = sl->th->co[l][j];
						if(strcmp(files[fliqz] , string_novalue) != 0 || strcmp(files[fliqzwriteend] , string_novalue) != 0) sl->thzplot[i][l] = sl->th[l][j];
					//	if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thizplot->co[i][l] = sl->SS->thi->co[l][j];
						if(strcmp(files[ficez] , string_novalue) != 0 || strcmp(files[ficezwriteend] , string_novalue) != 0) sl->thizplot[i][l] = sl->SS->thi[l][j];
					}			
					for(l=0;l<=Nl;l++){
					//	if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot->co[i][l] = sl->SS->P->co[l][j];
						if(strcmp(files[fpsiz] , string_novalue) != 0 || strcmp(files[fpsizwriteend] , string_novalue) != 0) sl->Pzplot[i][l] = sl->SS->P[l][j];
					}
					
				//	snow data
				//	if(snow->S->lnum->co[r][c]>0){
					if(snow->S->lnum[r][c]>0){
						odpnt[osnowdepth][i-1] = 0.0;
						odpnt[oSWE][i-1] = 0.0;
						odpnt[osnowT][i-1] = 0.0;
					//	for(l=1;l<=snow->S->lnum->co[r][c];l++){
						for(l=1;l<=snow->S->lnum[r][c];l++){
						//	odpnt[osnowdepth][i-1] += snow->S->Dzl->co[l][r][c];
							odpnt[osnowdepth][i-1] += snow->S->Dzl[l][r][c];
						//	odpnt[oSWE][i-1] += 1.0E+3*(snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c])/GTConst::rho_w;
							odpnt[oSWE][i-1] += 1.0E+3*(snow->S->w_liq[l][r][c]+snow->S->w_ice[l][r][c])/GTConst::rho_w;
						//	odpnt[osnowT][i-1] += snow->S->T->co[l][r][c]*snow->S->Dzl->co[l][r][c];
							odpnt[osnowT][i-1] += snow->S->T[l][r][c]*snow->S->Dzl[l][r][c];
						}
						odpnt[osnowdens][i-1] = odpnt[oSWE][i-1]*GTConst::rho_w/odpnt[osnowdepth][i-1];
						odpnt[osnowT][i-1] /= odpnt[osnowdepth][i-1];   
					}else{
						odpnt[osnowdepth][i-1] = 0.0;
						odpnt[oSWE][i-1] = 0.0;
						odpnt[osnowdens][i-1] = (double)number_novalue;
						odpnt[osnowT][i-1] = (double)number_novalue;
					}
					
					//glacier data
					if(par->max_glac_layers>0){
					//	if(glac->G->lnum->co[r][c]>0){
						if(glac->G->lnum[r][c]>0){
							odpnt[oglacdepth][i-1] = 0.0;
							odpnt[oGWE][i-1] = 0.0;
							odpnt[oglacT][i-1] = 0.0;
						//	for(l=1;l<=glac->G->lnum->co[r][c];l++){
							for(l=1;l<=glac->G->lnum[r][c];l++){
							//	odpnt[oglacdepth][i-1] += glac->G->Dzl->co[l][r][c];
								odpnt[oglacdepth][i-1] += glac->G->Dzl[l][r][c];
							//	odpnt[oGWE][i-1] += 1.0E+3*(glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c])/GTConst::rho_w;
								odpnt[oGWE][i-1] += 1.0E+3*(glac->G->w_liq[l][r][c]+glac->G->w_ice[l][r][c])/GTConst::rho_w;
							//	odpnt[oglacT][i-1] += glac->G->T->co[l][r][c]*glac->G->Dzl->co[l][r][c];
								odpnt[oglacT][i-1] += glac->G->T[l][r][c]*glac->G->Dzl[l][r][c];
							}
							odpnt[oglacdens][i-1] = odpnt[oGWE][i-1]*GTConst::rho_w/odpnt[oglacdepth][i-1];
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
						//	temp2 = join_strings(temp1, rec);
							temp2 = temp1+ rec;
						//	name = join_strings(temp2, textfile);
							name = temp2+ textfile;
						//	free(temp2);
						}else if (par->n_ContRecovery>0) {
						//	temp2 = join_strings(temp1, crec);
							temp2 = temp1+ crec;
						//	name = join_strings(temp2, textfile);
							name = temp2+ textfile;
						//	free(temp2);
						}else {
						//	name = join_strings(temp1, textfile);
							name = temp1+ textfile;
						}
						
					//	free(temp1);
						
						f=fopen(name.c_str() ,"a");
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
								//	fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
									fprintf(f, "%f",JDfrom0-par->init_date[i_sim]);
								}else if (opnt[j] == operiod) {
									fprintf(f, "%ld",i_sim);
								}else if (opnt[j] == orun) {
									fprintf(f, "%ld",i_run);
								}else if (opnt[j] == opoint) {
								//	fprintf(f, "%ld",par->IDpoint->co[i]);
									fprintf(f, "%ld",par->IDpoint[i]);
								}else {
									fprintf(f,"%f",odpnt[opnt[j]][i-1]);
								}
							}else {
								fprintf(f,"%ld",number_novalue);
							}
						}
						fprintf(f,"\n");
						fclose(f);
					//	free(name);
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
								//	fprintf(ffpoint, "%f",JDfrom0-par->init_date->co[i_sim]);
									fprintf(ffpoint, "%f",JDfrom0-par->init_date[i_sim]);
								}else if (opnt[j] == operiod) {
									fprintf(ffpoint, "%ld",i_sim);
								}else if (opnt[j] == orun) {
									fprintf(ffpoint, "%ld",i_run);		
								}else if (opnt[j] == opoint) {
									//fprintf(ffpoint, "%ld",par->IDpoint->co[i]);
									fprintf(ffpoint, "%ld",par->IDpoint[i]);
								}else {
									fprintf(ffpoint,"%f",odpnt[opnt[j]][i-1]);
								}
							}else {
								fprintf(ffpoint,"%ld",number_novalue);
							}
						}
						fprintf(ffpoint,"\n");
					}
										
					//Snow
					if(strcmp(files[fsnz] , string_novalue) != 0){
						
						temp1=join_strings(files[fsnz],NNNN);
						
						if (par->recover>0) {
						//	temp2 = join_strings(temp1, rec);
							temp2 = temp1 + rec;
						//	name = join_strings(temp2, textfile);
							name = temp2 + textfile;
						//	free(temp2);
						}else if (par->n_ContRecovery>0) {
						//	temp2 = join_strings(temp1, crec);
							temp2 = temp1 + crec;
						//	name = join_strings(temp2, textfile);
							name = temp2 + textfile;
						//	free(temp2);
						}else {
						//	name = join_strings(temp1, textfile);
							name = temp1 + textfile;
						}
						
					//	free(temp1);
						
						f=fopen(name.c_str() ,"a");

					//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
						if ((long)par->snow_plot_depths[1] != number_novalue) {
						//	m = par->snow_plot_depths->nh;
							m = par->snow_plot_depths.size();
						}else {
							m = par->max_snow_layers;
						}
												
						first_column=1;
						for (j=0; j<nosnw; j++) {
							if(first_column==0){
								fprintf(f,",");
							}else {
								first_column = 0;
							}
							if (osnw[j] >= 0) {
								if (osnw[j] == 0) {
									fprintf(f,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (osnw[j] == 1) {
									fprintf(f, "%f",JDfrom0);
								}else if (osnw[j] == 2) {
								//	fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
									fprintf(f, "%f",JDfrom0-par->init_date[i_sim]);
								}else if (osnw[j] == 3) {
									fprintf(f, "%ld",i_sim);
								}else if (osnw[j] == 4) {
									fprintf(f, "%ld",i_run);								
								}else if (osnw[j] == 5) {
								//	fprintf(f, "%ld",par->IDpoint->co[i]);
									fprintf(f, "%ld",par->IDpoint[i]);
								}else if (osnw[j] <= 5 + 1*m) {
									l = osnw[j] - 5 - 0*m;
								//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
									if ((long)par->snow_plot_depths[1] != number_novalue) {
									//	fprintf(f, "%f",interpolate_snow(r, c, par->snow_plot_depths->co[l]*cosslope, snow->S->lnum->co[r][c], snow->S->Dzl, snow->S->T));
										fprintf(f, "%f",interpolate_snow(r, c, par->snow_plot_depths[l]*cosslope, snow->S->lnum[r][c], snow->S->Dzl, snow->S->T));
									}else {
									//	fprintf(f, "%f",snow->S->T->co[l][r][c]);
										fprintf(f, "%f",snow->S->T[l][r][c]);
									}
								}else if (osnw[j] <= 5 + 2*m) {
									l = osnw[j] - 5 - 1*m;
								//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
									if ((long)par->snow_plot_depths[1] != number_novalue) {
									//	fprintf(f, "%f",interpolate_snow(r, c, par->snow_plot_depths->co[l]*cosslope, snow->S->lnum->co[r][c], snow->S->Dzl, snow->S->w_ice));
										fprintf(f, "%f",interpolate_snow(r, c, par->snow_plot_depths[l]*cosslope, snow->S->lnum[r][c], snow->S->Dzl, snow->S->w_ice));
									}else {
									//	fprintf(f, "%f",snow->S->w_ice->co[l][r][c]);
										fprintf(f, "%f",snow->S->w_ice[l][r][c]);
									}
								}else if (osnw[j] <= 5 + 3*m) {
									l = osnw[j] - 5 - 2*m;
								//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
									if ((long)par->snow_plot_depths[1] != number_novalue) {
									//	fprintf(f, "%f",interpolate_snow(r, c, par->snow_plot_depths->co[l]*cosslope, snow->S->lnum->co[r][c], snow->S->Dzl, snow->S->w_liq));
										fprintf(f, "%f",interpolate_snow(r, c, par->snow_plot_depths[l]*cosslope, snow->S->lnum[r][c], snow->S->Dzl, snow->S->w_liq));
									}else {
									//	fprintf(f, "%f",snow->S->w_liq->co[l][r][c]);
										fprintf(f, "%f",snow->S->w_liq[l][r][c]);
									}
								}else if (osnw[j] <= 5 + 3*m + par->max_snow_layers) {
									l = osnw[j] - 5 - 3*m;
								//	fprintf(f, "%f",snow->S->Dzl->co[l][r][c]);
									fprintf(f, "%f",snow->S->Dzl[l][r][c]);
								}
							}else {
								fprintf(f,"%f",(double)number_novalue);
							}
						}
						fprintf(f,"\n");
						fclose(f);
					//	free(name);
					}
					
					if(strcmp(files[fsnzwriteend] , string_novalue) != 0){
					//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
						if ((long)par->snow_plot_depths[1] != number_novalue) {
						//	m = par->snow_plot_depths->nh;
							m = par->snow_plot_depths.size();
						}else {
							m = par->max_snow_layers;
						}
						first_column=1;
						for (j=0; j<nosnw; j++) {
							if(first_column==0){
								fprintf(ffsnow,",");
							}else {
								first_column = 0;
							}
							if (osnw[j] >= 0) {
								if (osnw[j] == 0) {
									fprintf(ffsnow,"%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)day,(float)month,(float)year,(float)hour,(float)minute);		
								}else if (osnw[j] == 1) {
									fprintf(ffsnow, "%f",JDfrom0);
								}else if (osnw[j] == 2) {
								//	fprintf(ffsnow, "%f",JDfrom0-par->init_date->co[i_sim]);
									fprintf(ffsnow, "%f",JDfrom0-par->init_date[i_sim]);
								}else if (osnw[j] == 3) {
									fprintf(ffsnow, "%ld",i_sim);
								}else if (osnw[j] == 4) {
									fprintf(ffsnow, "%ld",i_run);								
								}else if (osnw[j] == 5) {
								//	fprintf(ffsnow, "%ld",par->IDpoint->co[i]);
									fprintf(ffsnow, "%ld",par->IDpoint[i]);
								}else if (osnw[j] <= 5 + 1*m) {
									l = osnw[j] - 5 - 0*m;
								//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
									if ((long)par->snow_plot_depths[1] != number_novalue) {
									//	fprintf(ffsnow, "%f",interpolate_snow(r, c, par->snow_plot_depths->co[l]*cosslope, snow->S->lnum->co[r][c], snow->S->Dzl, snow->S->T));
										fprintf(ffsnow, "%f",interpolate_snow(r, c, par->snow_plot_depths[l]*cosslope, snow->S->lnum[r][c], snow->S->Dzl, snow->S->T));
									}else {
									//	fprintf(ffsnow, "%f",snow->S->T->co[l][r][c]);
										fprintf(ffsnow, "%f",snow->S->T[l][r][c]);
									}
								}else if (osnw[j] <= 5 + 2*m) {
									l = osnw[j] - 5 - 1*m;
								//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
									if ((long)par->snow_plot_depths[1] != number_novalue) {
									//	fprintf(ffsnow, "%f",interpolate_snow(r, c, par->snow_plot_depths->co[l]*cosslope, snow->S->lnum->co[r][c], snow->S->Dzl, snow->S->w_ice));
										fprintf(ffsnow, "%f",interpolate_snow(r, c, par->snow_plot_depths[l]*cosslope, snow->S->lnum[r][c], snow->S->Dzl, snow->S->w_ice));
									}else {
									//	fprintf(ffsnow, "%f",snow->S->w_ice->co[l][r][c]);
										fprintf(ffsnow, "%f",snow->S->w_ice[l][r][c]);
									}
								}else if (osnw[j] <= 5 + 3*m) {
									l = osnw[j] - 5 - 2*m;
								//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
									if ((long)par->snow_plot_depths[1] != number_novalue) {
									//	fprintf(ffsnow, "%f",interpolate_snow(r, c, par->snow_plot_depths->co[l]*cosslope, snow->S->lnum->co[r][c], snow->S->Dzl, snow->S->w_liq));
										fprintf(ffsnow, "%f",interpolate_snow(r, c, par->snow_plot_depths[l]*cosslope, snow->S->lnum[r][c], snow->S->Dzl, snow->S->w_liq));
									}else {
									//	fprintf(ffsnow, "%f",snow->S->w_liq->co[l][r][c]);
										fprintf(ffsnow, "%f",snow->S->w_liq[l][r][c]);
									}
								}else if (osnw[j] <= 5 + 3*m + par->max_snow_layers) {
									l = osnw[j] - 5 - 3*m;
								//	fprintf(ffsnow, "%f",snow->S->Dzl->co[l][r][c]);
									fprintf(ffsnow, "%f",snow->S->Dzl[l][r][c]);
								}
							}else {
								fprintf(ffsnow,"%f",(double)number_novalue);
							}
						}
						fprintf(ffsnow,"\n");
					}					
					
					//Glacier
					if(par->max_glac_layers>0){
						
						if(strcmp(files[fglz] , string_novalue) != 0){
							
							temp1=join_strings(files[fglz],NNNN);
							
							if (par->recover>0) {
							//	temp2 = join_strings(temp1, rec);
								temp2 = temp1 + rec;
							//	name = join_strings(temp2, textfile);
								name = temp2 + textfile;
							//	free(temp2);
							}else if (par->n_ContRecovery>0) {
							//	temp2 = join_strings(temp1, crec);
								temp2 = temp1 + crec;
							//	name = join_strings(temp2, textfile);
								name = temp2 + textfile;
							//	free(temp2);
							}else {
							//	name = join_strings(temp1, textfile);
								name = temp1 + textfile;
							}
							
						//	free(temp1);
							
							f=fopen(name.c_str() ,"a");
							
						//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
							if ((long)par->glac_plot_depths[1] != number_novalue) {
							//	m = par->glac_plot_depths->nh;
								m = par->glac_plot_depths.size();
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
									//	fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
										fprintf(f, "%f",JDfrom0-par->init_date[i_sim]);
									}else if (oglc[j] == 3) {
										fprintf(f, "%ld",i_sim);
									}else if (oglc[j] == 4) {
										fprintf(f, "%ld",i_run);																
									}else if (oglc[j] == 5) {
									//	fprintf(f, "%ld",par->IDpoint->co[i]);
										fprintf(f, "%ld",par->IDpoint[i]);
									}else if (oglc[j] <= 5 + 1*m) {
										l = oglc[j] - 5 - 0*m;
									//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
										if ((long)par->glac_plot_depths[1] != number_novalue) {
										//	fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->T));
											fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths[l]*cosslope, glac->G->lnum[r][c], glac->G->Dzl, glac->G->T));
										}else {
										//	fprintf(f, "%f", glac->G->T->co[l][r][c]);
											fprintf(f, "%f", glac->G->T[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 2*m) {
										l = oglc[j] - 5 - 1*m;
									//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
										if ((long)par->glac_plot_depths[1] != number_novalue) {
										//	fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_ice));
											fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths[l]*cosslope, glac->G->lnum[r][c], glac->G->Dzl, glac->G->w_ice));
										}else {
										//	fprintf(f, "%f",glac->G->w_ice->co[l][r][c]);
											fprintf(f, "%f",glac->G->w_ice[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m) {
										l = oglc[j] - 5 - 2*m;
									//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
										if ((long)par->glac_plot_depths[1] != number_novalue) {
										//	fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_liq));
											fprintf(f, "%f",interpolate_snow(r, c, par->glac_plot_depths[l]*cosslope, glac->G->lnum[r][c], glac->G->Dzl, glac->G->w_liq));
										}else {
										//	fprintf(f, "%f",glac->G->w_liq->co[l][r][c]);
											fprintf(f, "%f",glac->G->w_liq[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m + par->max_glac_layers) {
										l = oglc[j] - 5 - 3*m;
									//	fprintf(f, "%f",glac->G->Dzl->co[l][r][c]);
										fprintf(f, "%f",glac->G->Dzl[l][r][c]);
									}
								}else {
									fprintf(f,"%f",(double)number_novalue);
								}
							}
							fprintf(f,"\n");
							fclose(f);
						//	free(name);
						}

						if(strcmp(files[fglzwriteend] , string_novalue) != 0){
						//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
							if ((long)par->glac_plot_depths[1] != number_novalue) {
							//	m = par->glac_plot_depths->nh;
								m = par->glac_plot_depths.size();
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
									//	fprintf(ffglac, "%f",JDfrom0-par->init_date->co[i_sim]);
										fprintf(ffglac, "%f",JDfrom0-par->init_date[i_sim]);
									}else if (oglc[j] == 3) {
										fprintf(ffglac, "%ld",i_sim);
									}else if (oglc[j] == 4) {
										fprintf(ffglac, "%ld",i_run);																
									}else if (oglc[j] == 5) {
									//	fprintf(ffglac, "%ld",par->IDpoint->co[i]);
										fprintf(ffglac, "%ld",par->IDpoint[i]);
									}else if (oglc[j] <= 5 + 1*m) {
										l = oglc[j] - 5 - 0*m;
									//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
										if ((long)par->glac_plot_depths[1] != number_novalue) {
										//	fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->T));
											fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths[l]*cosslope, glac->G->lnum[r][c], glac->G->Dzl, glac->G->T));
										}else {
										//	fprintf(ffglac, "%f", glac->G->T->co[l][r][c]);
											fprintf(ffglac, "%f", glac->G->T[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 2*m) {
										l = oglc[j] - 5 - 1*m;
									//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
										if ((long)par->glac_plot_depths[1] != number_novalue) {
										//	fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_ice));
											fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths[l]*cosslope, glac->G->lnum[r][c], glac->G->Dzl, glac->G->w_ice));
										}else {
										//	fprintf(ffglac, "%f",glac->G->w_ice->co[l][r][c]);
											fprintf(ffglac, "%f",glac->G->w_ice[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m) {
										l = oglc[j] - 5 - 2*m;
									//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
										if ((long)par->glac_plot_depths[1] != number_novalue) {
										//	fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths->co[l]*cosslope, glac->G->lnum->co[r][c], glac->G->Dzl, glac->G->w_liq));
											fprintf(ffglac, "%f",interpolate_snow(r, c, par->glac_plot_depths[l]*cosslope, glac->G->lnum[r][c], glac->G->Dzl, glac->G->w_liq));
										}else {
										//	fprintf(ffglac, "%f",glac->G->w_liq->co[l][r][c]);
											fprintf(ffglac, "%f",glac->G->w_liq[l][r][c]);
										}
									}else if (oglc[j] <= 5 + 3*m + par->max_glac_layers) {
										l = oglc[j] - 5 - 3*m;
									//	fprintf(ffglac, "%f",glac->G->Dzl->co[l][r][c]);
										fprintf(ffglac, "%f",glac->G->Dzl[l][r][c]);
									}
								}else {
									fprintf(ffglac,"%f",(double)number_novalue);
								}
							}
							fprintf(ffglac,"\n");
						}
						
					}
					
					//sl output
				//	write_soil_output(i, par->IDpoint->co[i], par->init_date->co[i_sim], JDfrom0, JD, day, month, year, hour, minute, par->soil_plot_depths, sl, par, GTConst::PsiMin, cosslope);
					write_soil_output(i, par->IDpoint[i], par->init_date[i_sim], JDfrom0, JD, day, month, year, hour, minute, par->soil_plot_depths, sl, par, GTConst::PsiMin, cosslope);
				//	initialize
					for(j=0;j<otot;j++) { odpnt[j][i-1]=0.0; }

				}
				
			}
			
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
			
			flog = fopen(logfile.c_str(), "a");
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
//	if(par->Dtplot_basin->co[i_sim] > 1.E-5 && par->state_basin == 1){
	if(par->Dtplot_basin[i_sim] > 1.E-5 && par->state_basin == 1){
		
		t_basin += par->Dt;
		
	//	if (fabs(t_basin - par->Dtplot_basin->co[i_sim])<1.E-5){
		if (fabs(t_basin - par->Dtplot_basin[i_sim])<1.E-5){
			
			if(strcmp(files[fbas] , string_novalue) != 0){
				
				if (par->recover > 0) write_suffix(rec, par->recover, 4);
				if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

				if (par->recover>0) {
					temp1 = join_strings(files[fbas], rec);
				//	name = join_strings(temp1, textfile);
					name = temp1 + textfile;
				//	free(temp1);
				}else if (par->n_ContRecovery>0) {
					temp1 = join_strings(files[fbas], crec);
				//	name = join_strings(temp1, textfile);
					name = temp1 + textfile;
				//	free(temp1);
				}else {
				//	name = join_strings(files[fbas], textfile);
					name = files[fbas] + string(textfile);
				}
				
				f=fopen(name.c_str() ,"a");
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
						//	fprintf(f, "%f",JDfrom0-par->init_date->co[i_sim]);
							fprintf(f, "%f",JDfrom0-par->init_date[i_sim]);
						}else {
							fprintf(f,"%f",odbsn[obsn[j]]);
						}
					}else {
						fprintf(f,"%f",(double)number_novalue);
					}
				}
				fprintf(f,"\n");
				fclose(f);
			//	free(name);
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
						//	fprintf(ffbas, "%f",JDfrom0-par->init_date->co[i_sim]);
							fprintf(ffbas, "%f",JDfrom0-par->init_date[i_sim]);
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
		//	printf(" SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Pvert=%6.2f mm Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm Mean Time Step=%f s\n\n",
		//		   odbsn[ooSW],odbsn[ooLW],odbsn[ooH],odbsn[ooLE],odbsn[oopnet],odbsn[oorainover],
		//		   odbsn[oosnowover],odbsn[oomasserror]*3600.0/par->Dtplot_basin->co[i_sim],mass_error_tot,odbsn[ootimestep]);
			printf(" SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Pvert=%6.2f mm Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm Mean Time Step=%f s\n\n",
					odbsn[ooSW],odbsn[ooLW],odbsn[ooH],odbsn[ooLE],odbsn[oopnet],odbsn[oorainover],
					odbsn[oosnowover],odbsn[oomasserror]*3600.0/par->Dtplot_basin[i_sim],mass_error_tot,odbsn[ootimestep]);
			
			flog = fopen(logfile.c_str(), "a");
			fprintf(flog,"\n%ld/%ld/%ld %ld:%02.0f JD:%f (%ld^ simulation day) %5.2f%% completed! \n",
					day,month,year,hour,(float)minute,JD,(long)(floor(times->time/86400))+1,
					percent_done);			
			fprintf(flog," t_meteo:%6.2f s, t_energy:%6.2f s, t_blowingsnow:%6.2f s, t_water:%6.2f s, t_sub:%6.2f s, t_sup:%6.2f s, t_out:%6.2f s\n",t_meteo,t_energy,t_blowingsnow,t_water,t_sub,t_sup,t_out);
		//	fprintf(flog," SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Pvert=%6.2f mm Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm Mean Time Step=%f s\n\n",
		//			odbsn[ooSW],odbsn[ooLW],odbsn[ooH],odbsn[ooLE],odbsn[oopnet],odbsn[oorainover],
		//			odbsn[oosnowover],odbsn[oomasserror]*3600.0/par->Dtplot_basin->co[i_sim],mass_error_tot,odbsn[ootimestep]);

			fprintf(flog," SW=%6.2f W/m2  LW:%6.2f W/m2  H=%6.2f W/m2  LE=%6.2f W/m2 \n Pvert=%6.2f mm Prain=%6.2f mm  Psnow=%6.2f mm  \n Max Error Richards=%e mm/h \n Tot Error Richards=%e mm Mean Time Step=%f s\n\n",
							odbsn[ooSW],odbsn[ooLW],odbsn[ooH],odbsn[ooLE],odbsn[oopnet],odbsn[oorainover],
							odbsn[oosnowover],odbsn[oomasserror]*3600.0/par->Dtplot_basin[i_sim],mass_error_tot,odbsn[ootimestep]);
			
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
//	if(par->output_meteo->co[i_sim]>0){
	if(par->output_meteo[i_sim]>0){
		if(strcmp(files[fTa] , string_novalue) != 0) {
			for(i=1;i<=par->total_pixel;i++){
			//	met->Tamean->co[i]+=met->Tgrid->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
				met->Tamean[i]+=met->Tgrid[top->rc_cont[i][1]][top->rc_cont[i][2]]/((par->output_meteo[i_sim]*3600.0)/(par->Dt));
			}
		}			
		if(strcmp(files[fwspd] , string_novalue) != 0){
			for(i=1;i<=par->total_pixel;i++){
			//	met->Vspdmean->co[i]+=met->Vgrid->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
				met->Vspdmean[i]+=met->Vgrid[top->rc_cont[i][1]][top->rc_cont[i][2]]/((par->output_meteo[i_sim]*3600.0)/(par->Dt));
			}
		}			
		if(strcmp(files[fwdir] , string_novalue) != 0){
			for(i=1;i<=par->total_pixel;i++){
			//	met->Vdirmean->co[i]+=met->Vdir->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
				met->Vdirmean[i]+=met->Vdir[top->rc_cont[i][1]][top->rc_cont[i][2]]/((par->output_meteo[i_sim]*3600.0)/(par->Dt));
			}
		}			
		if(strcmp(files[frh] , string_novalue) != 0){
			for(i=1;i<=par->total_pixel;i++){
			//	met->RHmean->co[i]+=met->RHgrid->co[top->rc_cont->co[i][1]][top->rc_cont->co[i][2]]/((par->output_meteo->co[i_sim]*3600.0)/(par->Dt));
				met->RHmean[i]+=met->RHgrid[top->rc_cont[i][1]][top->rc_cont[i][2]]/((par->output_meteo[i_sim]*3600.0)/(par->Dt));
			}
		}			
	}			
		
//	if (par->output_soil->co[i_sim]>0) {
	if (par->output_soil[i_sim]>0) {
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0){
			for (i=1; i<=par->total_pixel; i++) {
				for (l=1; l<=Nl; l++) {
				//	sl->T_av_tensor->co[l][i] += sl->SS->T->co[l][i]/((par->output_soil->co[i_sim]*3600.0)/(par->Dt));
					sl->T_av_tensor[l][i] += sl->SS->T[l][i]/((par->output_soil[i_sim]*3600.0)/(par->Dt));
				}
			}
		}
		if(strcmp(files[fliqav] , string_novalue) != 0){
			for (i=1; i<=par->total_pixel; i++) {
				for (l=1; l<=Nl; l++) {
				//	sl->thw_av_tensor->co[l][i] += sl->th->co[l][i]/((par->output_soil->co[i_sim]*3600.0)/(par->Dt));
					sl->thw_av_tensor[l][i] += sl->th[l][i]/((par->output_soil[i_sim]*3600.0)/(par->Dt));
				}
			}
		}
		if(strcmp(files[ficeav] , string_novalue) != 0){
			for (i=1; i<=par->total_pixel; i++) {
				for (l=1; l<=Nl; l++) {
				//	sl->thi_av_tensor->co[l][i] += sl->SS->thi->co[l][i]/((par->output_soil->co[i_sim]*3600.0)/(par->Dt));
					sl->thi_av_tensor[l][i] += sl->SS->thi[l][i]/((par->output_soil[i_sim]*3600.0)/(par->Dt));
				}
			}
		}
	}	
	
	
//	V=new_doublevector(par->total_pixel);
//	initialize_doublevector(V, (double)number_novalue);
	V.resize(par->total_pixel+1,(double)number_novalue);

//	soil properties
//	if(par->output_soil->co[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_soil->co[i_sim]*3600.0)<1.E-5){
	if(par->output_soil[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_soil[i_sim]*3600.0)<1.E-5){
	//	n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_soil->co[i_sim]*3600.0));
		n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_soil[i_sim]*3600.0));

		write_suffix(NNNNN, n_file, 1);
	//	if (par->run_times->co[i_sim] == 1) {
		if (par->run_times[i_sim] == 1) {
		//	s1 = join_strings(NNNNN, "");
			s1 = NNNNN + string("");
		}else {
		//	s1 = join_strings(NNNNN, RRRRR);
			s1 = string(NNNNN) + string(RRRRR) ;
		}
	//	if (par->init_date->nh == 1) {
		if (par->init_date.size() == 1) {
		//	s2 = join_strings(s1, "");
			s2 = s1 + "";
		}else {
		//	s2 = join_strings(s1, SSSSS);
			s2 = s1 + SSSSS;
		}
	//	free(s1);

		//theta liq tensor
		if(strcmp(files[fliq] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fliq], 0, par->format_out, sl->th, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fliq], 0, par->format_out, sl->th, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}

		//theta liq surface
		if(strcmp(files[fliqsup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
			//	V->co[i] = sl->th->co[1][i];
				V[i] = sl->th[1][i];
			}

		//	temp1=join_strings(files[fliqsup],s2);
			temp1= files[fliqsup] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
		
		//write thw_av tensor
		if(strcmp(files[fliqav] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fliqav], 0, par->format_out, sl->thw_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fliqav], 0, par->format_out, sl->thw_av_tensor, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//initialize thw_av_tensor
	//	if(strcmp(files[fliqav] , string_novalue) != 0) initialize_doublematrix(sl->thw_av_tensor, 0.);
		if(strcmp(files[fliqav] , string_novalue) != 0) sl->thw_av_tensor.resize(sl->thw_av_tensor.getRows(),sl->thw_av_tensor.getCols(), 0.0);
		
		//write T tensor
		if(strcmp(files[fT] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fT], 0, par->format_out, sl->SS->T, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fT], 0, par->format_out, sl->SS->T, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//theta T surface
		if(strcmp(files[fTsup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
			//	V->co[i] = sl->SS->T->co[1][i];
				V[i] = sl->SS->T[1][i];
			}
			
		//	temp1=join_strings(files[fTsup],s2);
			temp1= files[fTsup] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
		
		//write Tav tensor
		if(strcmp(files[fTav] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fTav], 0, par->format_out, sl->T_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fTav], 0, par->format_out, sl->T_av_tensor, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//theta Tav surface
		if(strcmp(files[fTavsup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
			//	V->co[i] = sl->T_av_tensor->co[1][i];
				V[i] = sl->T_av_tensor[1][i];
			}

		//	temp1=join_strings(files[fTavsup],s2);
			temp1= files[fTavsup] + s2 ;

			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
		
		//initialize T_av_tensor
	//	if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0) initialize_doublematrix(sl->T_av_tensor, 0.);
		if(strcmp(files[fTav] , string_novalue) != 0 || strcmp(files[fTavsup] , string_novalue) != 0) sl->T_av_tensor.resize(sl->T_av_tensor.getRows(),sl->T_av_tensor.getCols(), 0.0);
		
		//theta_ice tensor
		if(strcmp(files[fice] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fice], 0, par->format_out, sl->SS->thi, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fice], 0, par->format_out, sl->SS->thi, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//theta_ice surface
		if(strcmp(files[ficesup] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
			//	V->co[i] = sl->SS->thi->co[1][i];
				V[i] = sl->SS->thi[1][i];
			}

		//	temp1=join_strings(files[ficesup],s2);
			temp1= files[ficesup] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
		
		//write thi_av tensor
		if(strcmp(files[ficeav] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[ficeav], 0, par->format_out, sl->thi_av_tensor, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[ficeav], 0, par->format_out, sl->thi_av_tensor, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//initialize thi_av_tensor
	//	if(strcmp(files[ficeav] , string_novalue) != 0) initialize_doublematrix(sl->thi_av_tensor, 0.);
		if(strcmp(files[ficeav] , string_novalue) != 0) sl->thi_av_tensor.resize(sl->thi_av_tensor.getRows(),sl->thi_av_tensor.getCols() , 0.0);
		
		//write psi tensors
		if(strcmp(files[fpsitot] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fpsitot], 0, par->format_out, sl->Ptot, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fpsitot], 0, par->format_out, sl->Ptot, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		if(strcmp(files[fpsiliq] , string_novalue) != 0){
		//	if ((long)par->soil_plot_depths->co[1] != number_novalue) {
			if ((long)par->soil_plot_depths[1] != number_novalue) {
				write_tensorseries_soil(1, s2, files[fpsiliq], 0, par->format_out, sl->SS->P, par->soil_plot_depths, top->j_cont, top->rc_cont, sl->pa, top->slope, par->output_vertical_distances);
			}else {
				write_tensorseries3_vector(s2,files[fpsiliq], 0, par->format_out, sl->SS->P, UV, number_novalue, top->j_cont, Nr, Nc);
			}
		}
		
		//calculate saturation front depth
		if( strcmp(files[fwtable_up] , string_novalue) != 0 ){
			
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = find_watertabledepth_up(i, sl->type->co[r][c], sl);//normal
				V[i] = find_watertabledepth_up(i, sl->type[r][c], sl);//normal
			//	if (V->co[i] < 1.E-5 && sl->SS->P->co[0][i] > 0) V->co[i] = -sl->SS->P->co[0][i] / cos(top->slope->co[r][c] * GTConst::Pi/180.);
				if (V[i] < 1.E-5 && sl->SS->P[0][i] > 0) V[i] = -sl->SS->P[0][i] / cos(top->slope[r][c] * GTConst::Pi/180.);
			}
		//	temp1=join_strings(files[fwtable_up],s2);
			temp1= files[fwtable_up] + s2;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}

		if( strcmp(files[fwtable_dw] , string_novalue) != 0 ){
			
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	if (sl->SS->P->co[0][i] > 0){
				if (sl->SS->P[0][i] > 0){
				//	V->co[i] = -sl->SS->P->co[0][i] / cos(top->slope->co[r][c] * GTConst::Pi/180.);
					V[i] = -sl->SS->P[0][i] / cos(top->slope[r][c] * GTConst::Pi/180.);
				}else {
				//	V->co[i] = find_watertabledepth_dw(i, sl->type->co[r][c], sl);//normal
					V[i] = find_watertabledepth_dw(i, sl->type[r][c], sl);//normal
				}
			}
		//	temp1=join_strings(files[fwtable_dw],s2);
			temp1= files[fwtable_dw] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
		
		//calculate active layer depth		   
		if( strcmp(files[fthawed_up] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = find_activelayerdepth_up(i, sl->type->co[r][c], sl);   //normal
				V[i] = find_activelayerdepth_up(i, sl->type[r][c], sl);   //normal
			}
		//	temp1=join_strings(files[fthawed_up],s2);
			temp1= files[fthawed_up] + s2;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}

		if( strcmp(files[fthawed_dw] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = find_activelayerdepth_dw(i, sl->type->co[r][c], sl);   //normal
				V[i] = find_activelayerdepth_dw(i, sl->type[r][c], sl);   //normal
			}
		//	temp1=join_strings(files[fthawed_dw],s2);
			temp1= files[fthawed_dw] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
		
		//WATER OVER THE SURFACE
		if( strcmp(files[fhsupland] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1]
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = Fmax(0, sl->SS->P->co[0][i]) / cos(top->slope->co[r][c] * GTConst::Pi/180.);
				V[i] = Fmax(0, sl->SS->P[0][i]) / cos(top->slope[r][c] * GTConst::Pi/180.);
			}
			
		//	temp1 = join_strings(files[fhsupland], s2);
			temp1 = files[fhsupland] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
			
		if( strcmp(files[fhsupch] , string_novalue) != 0 ){
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	if (cnet->ch->co[r][c]!=0) {
				if (cnet->ch[r][c]!=0) {
				//	V->co[i] = cnet->SS->P->co[0][cnet->ch->co[r][c]] / cos(top->slope->co[r][c] * GTConst::Pi/180.);
					V[i] = cnet->SS->P[0][cnet->ch[r][c]] / cos(top->slope[r][c] * GTConst::Pi/180.);
				}else {
				//	V->co[i] = (double)number_novalue;
					V[i] = (double)number_novalue;
				}
			}
			
		//	temp1 = join_strings(files[fhsupch], s2);
			temp1 = files[fhsupch] + s2;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}

	//	free(s2);
	}
	
	//snow properties
//	if(par->output_snow->co[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_snow->co[i_sim]*3600.0)<1.E-5){
	if(par->output_snow[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_snow[i_sim]*3600.0)<1.E-5){
	//	n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_snow->co[i_sim]*3600.0));
		n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_snow[i_sim]*3600.0));
		write_suffix(NNNNN, n_file, 1);
	//	if (par->run_times->co[i_sim] == 1) {
		if (par->run_times[i_sim] == 1) {
		//	s1 = join_strings(NNNNN, "");
			s1 = NNNNN + string("");
		}else {
		//	s1 = join_strings(NNNNN, RRRRR);
			s1 = NNNNN + string(RRRRR);
		}
	//	if (par->init_date->nh == 1) {
		if (par->init_date.size()-1 == 1) {
		//	s2 = join_strings(s1, "");
			s2 = s1 + "" ;
		}else {
		//	s2 = join_strings(s1, SSSSS);
			s2 = s1 + SSSSS ;
		}
	//	free(s1);
		
		if(strcmp(files[fsnowdepth] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = 0.;
				V[i] = 0.;
			//	for(l=1;l<=snow->S->lnum->co[r][c];l++){
				for(l=1;l<=snow->S->lnum[r][c];l++){
				//	V->co[i] += snow->S->Dzl->co[l][r][c];
					V[i] += snow->S->Dzl[l][r][c];
				}
			}
			
		//	temp1=join_strings(files[fsnowdepth],s2);
			temp1= files[fsnowdepth] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	meteoio_writeEsriasciiVector(s2, 0 , V, top->j_cont, Nr, Nc, UV, number_novalue);

		//	free(temp1);
		}
		
		if(strcmp(files[fsnowmelt] , string_novalue) != 0){
		//	temp1=join_strings(files[fsnowmelt],s2);
			temp1= files[fsnowmelt] + s2;
			write_map_vector(temp1, 0, par->format_out, snow->MELTED, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(snow->MELTED, 0.);
			snow->MELTED.resize(snow->MELTED.size(),0.0);
		//	free(temp1);
		}
		
		if(strcmp(files[fsnowsubl] , string_novalue) != 0){
		//	temp1=join_strings(files[fsnowsubl],s2);
			temp1= files[fsnowsubl] + s2 ;
			write_map_vector(temp1, 0, par->format_out, snow->SUBL, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(snow->SUBL, 0.);
			snow->SUBL.resize(snow->SUBL.size(),0.0);
		//	free(temp1);
		}
		
		if(strcmp(files[fswe] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = 0.;
				V[i] = 0.;
			//	for(l=1;l<=snow->S->lnum->co[r][c];l++){
				for(l=1;l<=snow->S->lnum[r][c];l++){
				//	V->co[i] += (snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
					V[i] += (snow->S->w_liq[l][r][c]+snow->S->w_ice[l][r][c]);
				}
			}
		//	temp1=join_strings(files[fswe],s2);
			temp1= files[fswe] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
			
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = 0.;
				V[i] = 0.;
				D = 0.;
			//	for(l=1;l<=snow->S->lnum->co[r][c];l++){
				for(l=1;l<=snow->S->lnum[r][c];l++){
				//	V->co[i] += (snow->S->w_liq->co[l][r][c]+snow->S->w_ice->co[l][r][c]);
					V[i] += (snow->S->w_liq[l][r][c]+snow->S->w_ice[l][r][c]);
					D += snow->S->Dzl[l][r][c];
				}
			//	V->co[i] /= (1.E-3*D);
				V[i] /= (1.E-3*D);
			}
		//	temp1=join_strings(files[fswe], "DENSITY");
			temp1= files[fswe] + string("DENSITY");
		//	temp2=join_strings(temp1,s2);
			temp2= temp1 + s2;
			write_map_vector(temp2, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp2);
		//	free(temp1);
			
			if(par->blowing_snow==1){
			//	temp1=join_strings(files[fswe],"WindTrans");
				temp1= files[fswe] + string("WindTrans");
			//	temp2=join_strings(temp1, s2);
				temp2= temp1 + s2;
				write_map(temp2, 0, par->format_out, snow->Wtrans_plot, UV, number_novalue);

				initmatrix(0.0, snow->Wtrans_plot, land->LC, number_novalue);
			//	free(temp2);
			//	free(temp1);
				
			//	temp1=join_strings(files[fswe],"WindSubl");
				temp1= files[fswe]+ string("WindSubl");
			//	temp2=join_strings(temp1, s2);
				temp2= temp1 + s2;
				write_map(temp2, 0, par->format_out, snow->Wsubl_plot, UV, number_novalue);
				initmatrix(0.0, snow->Wsubl_plot, land->LC, number_novalue);
			//	free(temp2);
			//	free(temp1);
				
			}
		}
		
		if(strcmp(files[fsndur] , string_novalue) != 0){
		//	temp1=join_strings(files[fsndur],s2);
			temp1= files[fsndur] + s2 ;
			write_map_vector(temp1, 0, par->format_out, snow->t_snow, UV, number_novalue, top->j_cont, Nr, Nc);
		//  initialize_doublevector(snow->t_snow, 0.);
			snow->t_snow.resize(snow->t_snow.size(), 0.0);
		//	free(temp1);
		}
		
	//	free(s2);
	}
	
	//glacier properties
//	if(par->max_glac_layers>0 && par->output_glac->co[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_glac->co[i_sim]*3600.0)<1.E-5){
	if(par->max_glac_layers>0 && par->output_glac[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_glac[i_sim]*3600.0)<1.E-5){

	//	n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_glac->co[i_sim]*3600.0));
		n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_glac[i_sim]*3600.0));
		write_suffix(NNNNN, n_file, 1);
	//	if (par->run_times->co[i_sim] == 1) {
		if (par->run_times[i_sim] == 1) {
		//	s1 = join_strings(NNNNN, "");
			s1 = NNNNN + string("");
		}else {
		//	s1 = join_strings(NNNNN, RRRRR);
			s1 = NNNNN + string(RRRRR);
		}
	//	if (par->init_date->nh == 1) {
		if (par->init_date.size() == 1) {
		//	s2 = join_strings(s1, "");
			s2 = s1 + "";
		}else {
		//	s2 = join_strings(s1, SSSSS);
			s2 = s1 + SSSSS;
		}
	//	free(s1);
		
		if(strcmp(files[fglacdepth] , string_novalue) != 0){
			
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = 0.;
				V[i] = 0.;
			//	for(l=1;l<=glac->G->lnum->co[r][c];l++){
				for(l=1;l<=glac->G->lnum[r][c];l++){
				//	V->co[i] += (glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
					V[i] += (glac->G->w_liq[l][r][c]+glac->G->w_ice[l][r][c]);
				}
			}
		//	temp1=join_strings(files[fglacdepth],s2);
			temp1= files[fglacdepth] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}
		
		if(strcmp(files[fglacmelt] , string_novalue) != 0){
		//	temp1=join_strings(files[fglacmelt],s2);
			temp1= files[fglacmelt] + s2 ;
			write_map_vector(temp1, 0, par->format_out, glac->MELTED, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(glac->MELTED, 0.);
			glac->MELTED.resize(glac->MELTED.size(), 0.);
		//	free(temp1);
		}
		
		if(strcmp(files[fglacsubl] , string_novalue) != 0){
		//	temp1=join_strings(files[fglacsubl],s2);
			temp1= files[fglacsubl] + s2 ;
			write_map_vector(temp1, 0, par->format_out, glac->SUBL, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(glac->SUBL, 0.);
			glac->SUBL.resize(glac->SUBL.size(), 0.);
		//	free(temp1);
		}
		
		if(strcmp(files[fgwe] , string_novalue) != 0){
			for(i=1; i<=par->total_pixel; i++){
			//	r = top->rc_cont->co[i][1];
				r = top->rc_cont[i][1];
			//	c = top->rc_cont->co[i][2];
				c = top->rc_cont[i][2];
			//	V->co[i] = 0.;
				V[i] = 0.;
			//	for(l=1;l<=glac->G->lnum->co[r][c];l++){
				for(l=1;l<=glac->G->lnum[r][c];l++){
				//	V->co[i] += (glac->G->w_liq->co[l][r][c]+glac->G->w_ice->co[l][r][c]);
					V[i] += (glac->G->w_liq[l][r][c]+glac->G->w_ice[l][r][c]);
				}
			}
		//	temp1=join_strings(files[fgwe],s2);
			temp1= files[fgwe] + s2 ;
			write_map_vector(temp1, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		}

	//	free(s2);

	}
	
	//SURFACE ENERGY BALANCE	
	
	//RADIATION
//	if(par->output_surfenergy->co[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_surfenergy->co[i_sim]*3600.0)<1.E-5){
	if(par->output_surfenergy[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_surfenergy[i_sim]*3600.0)<1.E-5){
		
	//	n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_surfenergy->co[i_sim]*3600.0));
		n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_surfenergy[i_sim]*3600.0));
		write_suffix(NNNNN, n_file, 1);
	//	if (par->run_times->co[i_sim] == 1) {
		if (par->run_times[i_sim] == 1) {
		//	s1 = join_strings(NNNNN, "");
			s1 = NNNNN + string("");
		}else {
		//	s1 = join_strings(NNNNN, RRRRR);
			s1 = NNNNN + string(RRRRR);
		}
	//	if (par->init_date->nh == 1) {
		if (par->init_date.size()-1 == 1) {
		//	s2 = join_strings(s1, "");
			s2 = s1 + "";
		}else {
		//	s2 = join_strings(s1, SSSSS);
			s2 = s1 + SSSSS;
		}
	//	free(s1);
		
		if(strcmp(files[fradnet] , string_novalue) != 0){			
		//	name=join_strings(files[fradnet], s2);
			name= files[fradnet] + s2 ;
			write_map_vector(name, 0, par->format_out, egy->Rn_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->Rn_mean, 0.);
			egy->Rn_mean.resize(egy->Rn_mean.size(),0.0);
		//	free(name);
		}
			
		if(strcmp(files[fradLWin] , string_novalue) != 0){			
		//	name=join_strings(files[fradLWin],s2);
			name= files[fradLWin] + s2 ;
			write_map_vector(name, 0, par->format_out, egy->LWin_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->LWin_mean, 0.);
			egy->LWin_mean.resize(egy->LWin_mean.size(),0.0);
		//	free(name);
		}
		
		if(strcmp(files[fradLW] , string_novalue) != 0){			
		//	name=join_strings(files[fradLW],s2);
			name= files[fradLW] + s2 ;
			write_map_vector(name , 0, par->format_out, egy->LW_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->LW_mean, 0.);
			egy->LW_mean.resize(egy->LW_mean.size(),0.0);
		//	free(name);
		}
			
		if(strcmp(files[fradSW] , string_novalue) != 0){			
		//	name=join_strings(files[fradSW],s2);
			name= files[fradSW] + s2 ;
			write_map_vector(name , 0, par->format_out, egy->SW_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->SW_mean, 0.);
			egy->SW_mean.resize(egy->SW_mean.size(),0.0);
		//	free(name);
		}
		
		if(strcmp(files[fradSWin] , string_novalue) != 0){			
		//	name=join_strings(files[fradSWin],s2);
			name= files[fradSWin] + s2;
			write_map_vector(name, 0, par->format_out, egy->Rswdown_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->Rswdown_mean, 0.);
			egy->Rswdown_mean.resize(egy->Rswdown_mean.size(),0.0);
		//	free(name);
		}
		
		if(strcmp(files[fradSWinbeam] , string_novalue) != 0){			
		//	name=join_strings(files[fradSWinbeam],s2);
			name= files[fradSWinbeam] + s2;
			write_map_vector(name, 0, par->format_out, egy->Rswbeam_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->Rswbeam_mean, 0.);
			egy->Rswbeam_mean.resize(egy->Rswbeam_mean.size(),0.0);
		//	free(name);
		}
			
		if(strcmp(files[fshadow] , string_novalue) != 0){			

			for(i=1; i<=par->total_pixel; i++){
			//	if(egy->nDt_sun->co[i]>0){
				if(egy->nDt_sun[i]>0){
				//	V->co[i] = egy->nDt_shadow->co[i]/(double)(egy->nDt_sun->co[i]);
					V[i] = egy->nDt_shadow[i]/(double)(egy->nDt_sun[i]);
				}else{
				//	V->co[i] = -1.;
					V[i] = -1.;
				}
			}
			
		//	name=join_strings(files[fshadow],s2);
			name= files[fshadow] + s2 ;
			write_map_vector(name, 0, par->format_out, V, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_longvector(egy->nDt_shadow, 0.);
			egy->nDt_shadow.resize(egy->nDt_shadow.size(),0.0);
		//	initialize_longvector(egy->nDt_sun, 0.);
			egy->nDt_sun.resize(egy->nDt_sun.size(),0.0);
		//	free(name);
		}
					
		//GROUND HEAT FLUX
		if(strcmp(files[fG] , string_novalue) != 0){
		//	name=join_strings(files[fG],s2);
			name= files[fG] + s2 ;
			write_map_vector(name, 0, par->format_out, egy->SEB_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->SEB_mean, 0.);
			egy->SEB_mean.resize(egy->SEB_mean.size(),0.0);
		//	free(name);
		}
				
		//SENSIBLE HEAT FLUX
		if(strcmp(files[fH] , string_novalue) != 0){
		//	name=join_strings(files[fH],s2);
			name= files[fH] + s2 ;
			write_map_vector(name, 0, par->format_out, egy->H_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->H_mean, 0.);
			egy->H_mean.resize(egy->H_mean.size(),0.0);
		//	free(name);
		}
		
		
		//LATENT HEAT FLUX
		if(strcmp(files[fLE] , string_novalue) != 0){
		//	name=join_strings(files[fLE],s2);
			name= files[fLE] + s2 ;
			write_map_vector(name, 0, par->format_out, egy->ET_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->ET_mean, 0.);
			egy->ET_mean.resize(egy->ET_mean.size(),0.0);
		//	free(name);
		}
		
		//SURFACE TEMPERATURE
		if(strcmp(files[fTs] , string_novalue) != 0){
		//	name=join_strings(files[fTs],s2);
			name= files[fTs] + s2;
			write_map_vector(name, 0, par->format_out, egy->Ts_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(egy->Ts_mean, 0.);
			egy->Ts_mean.resize(egy->Ts_mean.size(),0.0);
		//	free(name);
		}
		
	//	free(s2);

	}
	
	//vegetation variables	
//	if(par->output_vegetation->co[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_vegetation->co[i_sim]*3600.0)<1.E-5){
	if(par->output_vegetation[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_vegetation[i_sim]*3600.0)<1.E-5){
	//	n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_vegetation->co[i_sim]*3600.0));
		n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_vegetation[i_sim]*3600.0));
		write_suffix(NNNNN, n_file, 1);
	//	if (par->run_times->co[i_sim] == 1) {
		if (par->run_times[i_sim] == 1) {
		//	s1 = join_strings(NNNNN, "");
			s1 = NNNNN + string("");
		}else {
		//	s1 = join_strings(NNNNN, RRRRR);
			s1 = NNNNN + string(RRRRR);
		}
	//	if (par->init_date->nh == 1) {
		if (par->init_date.size() == 1) {
		//	s2 = join_strings(s1, "");
			s2 = s1 + "";
		}else {
		//	s2 = join_strings(s1, SSSSS);
			s2 = s1 + SSSSS;
		}
	//	free(s1);
		
		//INTERCEPTED PRECIPITATION
		if(strcmp(files[fcint] , string_novalue) != 0){
			
		//	temp1=join_strings(files[fcint],"water");
			temp1= files[fcint] + string("water");
		//	temp2=join_strings(temp1,s2);
			temp2= temp1 + s2;
			write_map_vector(temp2, 0, par->format_out, sl->VS->wrain, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		//	free(temp2);
			
		//	temp1=join_strings(files[fcint],"snow");
			temp1= files[fcint] + string("snow");
		//	temp2=join_strings(temp1,s2);
			temp2= temp1 + s2 ;
			write_map_vector(temp2, 0, par->format_out, sl->VS->wsnow, UV, number_novalue, top->j_cont, Nr, Nc);
		//	free(temp1);
		//	free(temp2);
		}
		
	//	free(s2);

	}
		
	//METEO
//	if(par->output_meteo->co[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_meteo->co[i_sim]*3600.0)<1.E-5){
	if(par->output_meteo[i_sim]>0 && fmod(times->time+par->Dt+par->delay_day_recover*86400.,par->output_meteo[i_sim]*3600.0)<1.E-5){
	//	n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_meteo->co[i_sim]*3600.0));
		n_file=(long)((times->time+par->Dt+par->delay_day_recover*86400.)/(par->output_meteo[i_sim]*3600.0));

		write_suffix(NNNNN, n_file, 1);
	//	if (par->run_times->co[i_sim] == 1) {
		if (par->run_times[i_sim] == 1) {
		//	s1 = join_strings(NNNNN, "");
			s1 = NNNNN + string("");
		}else {
		//	s1 = join_strings(NNNNN, RRRRR);
			s1 = NNNNN + string(RRRRR);
		}
	//	if (par->init_date->nh == 1) {
		if (par->init_date.size()-1 == 1) {
		//	s2 = join_strings(s1, "");
			s2 = s1 + "";
		}else {
		//	s2 = join_strings(s1, SSSSS);
			s2 = s1 + SSSSS ;
		}
	//	free(s1);
		
		//AIR TEMPERATURE
		if(strcmp(files[fTa] , string_novalue) != 0){
		//	name=join_strings(files[fTa],s2);
			name= files[fTa] + s2;
			write_map_vector(name, 0, par->format_out, met->Tamean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(met->Tamean, 0.);
			met->Tamean.resize(met->Tamean.size(),0.0);
		//	free(name);
		}
		
		//PRECIPITATION
		if(strcmp(files[fprec] , string_novalue) != 0){
			
		//	name=join_strings(files[fprec],"TOTAL");
			name= files[fprec]+string("TOTAL");
		//	temp1=join_strings(name , s2);
			temp1=name + s2;
			write_map_vector(temp1, 0, par->format_out, wat->PrTOT_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(wat->PrTOT_mean, 0.);
			wat->PrTOT_mean.resize(wat->PrTOT_mean.size(),0.0);

		//	free(temp1);
		//	free(name);
			
		//	name=join_strings(files[fprec],"SNOW");
			name= files[fprec]+ string("SNOW");
		//	temp1=join_strings(name, s2);
			temp1= name + s2;
			write_map_vector(temp1, 0, par->format_out, wat->PrSNW_mean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(wat->PrSNW_mean, 0.);
			wat->PrSNW_mean.resize(wat->PrSNW_mean.size(),0.0);
		//	free(temp1);
		//	free(name);
		}
		
		if(strcmp(files[fwspd] , string_novalue) != 0){
		//	name=join_strings(files[fwspd],s2);
			name= files[fwspd] + s2 ;
			write_map_vector(name, 0, par->format_out, met->Vspdmean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(met->Vspdmean, 0.);
			met->Vspdmean.resize(met->Vspdmean.size(), 0.0);
		//	free(name);
		}
		
		if(strcmp(files[fwdir] , string_novalue) != 0){
		//	name=join_strings(files[fwdir],s2);
			name= files[fwdir] + s2;
			write_map_vector(name, 0, par->format_out, met->Vdirmean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(met->Vdirmean, 0.);
			met->Vdirmean.resize(met->Vdirmean.size(), 0.0);
		//	free(name);
		}
		
		if(strcmp(files[frh] , string_novalue) != 0){
		//	name=join_strings(files[frh],s2);
			name= files[frh] + s2;
			write_map_vector(name, 0, par->format_out, met->RHmean, UV, number_novalue, top->j_cont, Nr, Nc);
		//	initialize_doublevector(met->RHmean, 0.);
 			met->RHmean.resize(met->RHmean.size(),0.0);
		//	free(name);
		}
		
	//	free(s2);
	}	
	
	//free_doublevector(V);
	
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	//SPECIAL PLOTS AT SOME DAYS
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	
//	if(times->JD_plots->nh > 1 && times->iplot<=times->JD_plots->nh){
	if(times->JD_plots.size() > 1 && times->iplot<=times->JD_plots.size()){
		
		i=times->iplot;
		j=2*i-1;
	//	if( fabs(par->init_date->co[i_sim]+(times->time+par->Dt)/86400. - times->JD_plots->co[j+1]) < 1.E-5 ){
	//	if( fabs(par->init_date[i_sim]+(times->time+par->Dt)/86400. - times->JD_plots->co[j+1]) < 1.E-5 ){
		if( fabs(par->init_date[i_sim]+(times->time+par->Dt)/86400. - times->JD_plots[j+1]) < 1.E-5 ){
			
			flog = fopen(logfile.c_str(), "a");
			fprintf(flog,"Printing plot number %ld \n",i);
			fclose(flog);
			
			printf("Printing plot number %ld \n",i);
			
		//	V = new_doublevector(par->total_pixel);
			V.resize(par->total_pixel+1);
			
			if(strcmp(files[pH] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
				//	V->co[i] = egy->Hgplot->co[i] + egy->Hvplot->co[i];
					V[i] = egy->Hgplot[i] + egy->Hvplot[i];
				}
				
				plot(files[pH], i, V, par->format_out, top->j_cont);
			}
			
			if(strcmp(files[pLE] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
				//	V->co[i] = egy->LEgplot->co[i] + egy->LEvplot->co[i];
					V[i] = egy->LEgplot[i] + egy->LEvplot[i];
				}
				
				plot(files[pLE], i, V, par->format_out, top->j_cont);
			}

			if(strcmp(files[pG] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
				//	V->co[i] = egy->SWgplot->co[i]+egy->LWgplot->co[i]-egy->Hgplot->co[i]-egy->LEgplot->co[i];
					V[i] = egy->SWgplot[i]+egy->LWgplot[i]-egy->Hgplot[i]-egy->LEgplot[i];
				}
				plot(files[pG], i, V, par->format_out, top->j_cont);
			}
			
			if(strcmp(files[pth] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
				//	V->co[i] = sl->th->co[1][i];
					V[i] = sl->th[1][i];
				}
				plot(files[pth], i, V, par->format_out, top->j_cont);
			}

			if(strcmp(files[pth] , string_novalue) != 0) {
				for (i=1; i<=par->total_pixel; i++) {
				//	V->co[i] = sl->th->co[1][i];
					V[i] = sl->th[1][i];
				}
				plot(files[pth], i, V, par->format_out, top->j_cont);
			}
			
			if(strcmp(files[pVspd] , string_novalue) != 0 ){
			   for (i=1; i<=par->total_pixel; i++) {
			  //   V->co[i] = sqrt(pow(met->Vxplot->co[i], 2.0) + pow(met->Vyplot->co[i], 2.0));
				   V[i] = sqrt(pow(met->Vxplot[i], 2.0) + pow(met->Vyplot[i], 2.0));
			   }
			   plot(files[pVspd], i, V, par->format_out, top->j_cont);
			}
			   
			if(strcmp(files[pVdir] , string_novalue) != 0){
				for (i=1; i<=par->total_pixel; i++) {
				//	V->co[i] = 270.0 - (180./GTConst::Pi)*atan2(met->Vyplot->co[i],met->Vxplot->co[i]);
					V[i] = 270.0 - (180./GTConst::Pi)*atan2(met->Vyplot[i],met->Vxplot[i]);
				//	if (V->co[i] >= 360.0) V->co[i] -= 360.0;
					if (V[i] >= 360.0) V[i] -= 360.0;
				}
				plot(files[pVdir], i, V, par->format_out, top->j_cont);
			}
			
			//free_doublevector(V);

			if(strcmp(files[pHg] , string_novalue) != 0)  plot(files[pHg], i, egy->Hgplot, par->format_out, top->j_cont);
			if(strcmp(files[pLEg] , string_novalue) != 0) plot(files[pLEg], i, egy->LEgplot, par->format_out, top->j_cont);
			if(strcmp(files[pHv] , string_novalue) != 0) plot(files[pHv], i, egy->Hvplot, par->format_out, top->j_cont);
			if(strcmp(files[pLEv] , string_novalue) != 0) plot(files[pLEv], i, egy->LEvplot, par->format_out, top->j_cont);
			if(strcmp(files[pSWin] , string_novalue) != 0)  plot(files[pSWin], i, egy->SWinplot, par->format_out, top->j_cont);
			if(strcmp(files[pSWg] , string_novalue) != 0)  plot(files[pSWg], i, egy->SWgplot, par->format_out, top->j_cont);
			if(strcmp(files[pSWv] , string_novalue) != 0)  plot(files[pSWv], i, egy->SWvplot, par->format_out, top->j_cont);
			if(strcmp(files[pLWin] , string_novalue) != 0) plot(files[pLWin], i, egy->LWinplot, par->format_out, top->j_cont);
			if(strcmp(files[pLWg] , string_novalue) != 0) plot(files[pLWg], i, egy->LWgplot, par->format_out, top->j_cont);
			if(strcmp(files[pLWv] , string_novalue) != 0)  plot(files[pLWv], i, egy->LWvplot, par->format_out, top->j_cont);
			if(strcmp(files[pTs] , string_novalue) != 0)  plot(files[pTs], i, egy->Tsplot, par->format_out, top->j_cont);
			if(strcmp(files[pTg] , string_novalue) != 0) plot(files[pTg], i, egy->Tgplot, par->format_out, top->j_cont);
			if(strcmp(files[pTv] , string_novalue) != 0)  plot(files[pTv], i, egy->Tvplot, par->format_out, top->j_cont);
			if(strcmp(files[pTa] , string_novalue) != 0) plot(files[pTa], i, met->Taplot, par->format_out, top->j_cont);
			if(strcmp(files[pD] , string_novalue) != 0) plot(files[pD], i, snow->Dplot, par->format_out, top->j_cont);
			if(strcmp(files[pRH] , string_novalue) != 0)  plot(files[pRH], i, met->RHplot, par->format_out, top->j_cont);
						
		//	if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->Hgplot, 0.);
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->Hgplot.resize(egy->Hgplot.size(), 0.0);
		//	if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->LEgplot, 0.);
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LEgplot.resize(egy->LEgplot.size(), 0.0);
		//	if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0) initialize_doublevector(egy->Hvplot, 0.);
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0) egy->Hvplot.resize(egy->Hvplot.size(), 0.);
		//	if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0) initialize_doublevector(egy->LEvplot, 0.);
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0) egy->LEvplot.resize(egy->LEvplot.size(), 0.);
		//	if(strcmp(files[pSWin] , string_novalue) != 0) initialize_doublevector(egy->SWinplot, 0.);
			if(strcmp(files[pSWin] , string_novalue) != 0) egy->SWinplot.resize(egy->SWinplot.size(), 0.);
		//	if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->SWgplot, 0.);
			if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->SWgplot.resize(egy->SWgplot.size(), 0.);
		//	if(strcmp(files[pSWv] , string_novalue) != 0) initialize_doublevector(egy->SWvplot, 0.);
			if(strcmp(files[pSWv] , string_novalue) != 0) egy->SWvplot.resize(egy->SWvplot.size(), 0.);
		//	if(strcmp(files[pLWin] , string_novalue) != 0) initialize_doublevector(egy->LWinplot, 0.);
			if(strcmp(files[pLWin] , string_novalue) != 0) egy->LWinplot.resize(egy->LWinplot.size(), 0.);
		//	if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) initialize_doublevector(egy->LWgplot, 0.);
			if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LWgplot.resize(egy->LWgplot.size(), 0.);
		//	if(strcmp(files[pLWv] , string_novalue) != 0) initialize_doublevector(egy->LWvplot, 0.);
			if(strcmp(files[pLWv] , string_novalue) != 0) egy->LWvplot.resize(egy->LWvplot.size(), 0.);
		//	if(strcmp(files[pTs] , string_novalue) != 0) initialize_doublevector(egy->Tsplot, 0.);
			if(strcmp(files[pTs] , string_novalue) != 0) egy->Tsplot.resize(egy->Tsplot.size(), 0.);
		//	if(strcmp(files[pTg] , string_novalue) != 0) initialize_doublevector(egy->Tgplot, 0.);
			if(strcmp(files[pTg] , string_novalue) != 0) egy->Tgplot.resize(egy->Tgplot.size(), 0.);
		//	if(strcmp(files[pTv] , string_novalue) != 0) initialize_doublevector(egy->Tvplot, 0.);
			if(strcmp(files[pTv] , string_novalue) != 0) egy->Tvplot.resize(egy->Tvplot.size(), 0.);
		//	if(strcmp(files[pD] , string_novalue) != 0) initialize_doublevector(snow->Dplot, 0.);
			if(strcmp(files[pD] , string_novalue) != 0) snow->Dplot.resize(snow->Dplot.size(), 0.);
		//	if(strcmp(files[pTa] , string_novalue) != 0) initialize_doublevector(met->Taplot, 0.);
			if(strcmp(files[pTa] , string_novalue) != 0) met->Taplot.resize(met->Taplot.size(), 0.);
			if(strcmp(files[pVspd] , string_novalue) != 0 || strcmp(files[pVdir] , string_novalue) != 0){
			//	initialize_doublevector(met->Vxplot, 0.);
				met->Vxplot.resize(met->Vxplot.size(), 0.);
			//	initialize_doublevector(met->Vyplot, 0.);
				met->Vyplot.resize(met->Vyplot.size(), 0.);
			}
		//	if(strcmp(files[pRH] , string_novalue) != 0) initialize_doublevector(met->RHplot, 0.);
			if(strcmp(files[pRH] , string_novalue) != 0) met->RHplot.resize(met->RHplot.size(), 0.);
			
			times->iplot ++;
		}
	}
	
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	//SAVING POINTS
	/**********************************************************************************************************/
	/**********************************************************************************************************/
	//char hold; cout << "point_sim=" << par->point_sim << " isavings=" << isavings << " par->saving_points.size()=" << par->saving_points.size() << " par->delay_day_recover=" << par->delay_day_recover; cin.get(hold);
	if(times->time==0) isavings=par->recover;
	
//	if(isavings < par->saving_points->nh){
	if(isavings < par->saving_points.size()-1){
	//	if(par->saving_points->nh==1 && par->saving_points->co[1]==0.0){
		if(par->saving_points.size()-1==1 && par->saving_points[1]==0.0){

			isavings=1;
			
		}else{
						
		//	if(times->time+par->Dt+par->delay_day_recover*86400. >= 86400.*par->saving_points->co[isavings+1]){
			if(times->time+par->Dt+par->delay_day_recover*86400. >= 86400.*par->saving_points[isavings+1]){
				isavings+=1;
				
				printf("Writing recovering files, saving point number %ld\n",isavings);
				
				flog = fopen(logfile.c_str(), "a");
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
				//	free(name);
				}
				
				if(strcmp(files[rwcsn] , string_novalue) != 0){
				//	name = join_strings(files[rwcsn],NNNN);
					name = files[rwcsn]+ string(NNNN);
					write_map_vector(name, 0, par->format_out, sl->VS->wsnow, UV, number_novalue, top->j_cont, Nr, Nc);
				//	free(name);
				}
									
				for (i=1; i<=par->total_pixel; i++) {
				//	if ((long)sl->VS->Tv->co[i] == number_novalue) sl->VS->Tv->co[i] = 0.;
					if ((long)sl->VS->Tv[i] == number_novalue) sl->VS->Tv[i] = 0.;
				}
				
				if(strcmp(files[rTv] , string_novalue) != 0){
				//	name = join_strings(files[rTv],NNNN);
					name = files[rTv] + string(NNNN);
					write_map_vector(name, 0, par->format_out, sl->VS->Tv, UV, number_novalue, top->j_cont, Nr, Nc);
				//	free(name);
				}	
									
				for(l=1;l<=par->max_snow_layers;l++){
					if(strcmp(files[rDzs] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rDzs], 0, par->format_out, snow->S->Dzl, UV, number_novalue);
					if(strcmp(files[rwls] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwls], 0, par->format_out, snow->S->w_liq, UV, number_novalue);
					if(strcmp(files[rwis] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwis], 0, par->format_out, snow->S->w_ice, UV, number_novalue);
					if(strcmp(files[rTs] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rTs], 0, par->format_out, snow->S->T, UV, number_novalue);
				}
				
				if(strcmp(files[rsnag] , string_novalue) != 0){
				//	name = join_strings(files[rsnag],NNNN);
					name = files[rsnag] + string(NNNN);
					write_map_vector(name, 0, par->format_out, snow->age, UV, number_novalue, top->j_cont, Nr, Nc);
				//	free(name);
				}
				
				if(strcmp(files[rns] , string_novalue) != 0){
				//	M=copydouble_longmatrix(snow->S->lnum);
				//	name = join_strings(files[rns], NNNN);
					name = files[rns] + string(NNNN);
				//	write_map(name, 1, par->format_out, M, UV, number_novalue);
					write_map(name, 1, par->format_out, snow->S->lnum, UV, number_novalue);
				//	free_doublematrix(M);
				//	free(name);
				}
				
				if(par->max_glac_layers>0){
					for(l=1;l<=par->max_glac_layers;l++){
						if(strcmp(files[rDzi] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rDzi], 0, par->format_out, glac->G->Dzl, UV, number_novalue);
						if(strcmp(files[rwli] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwli], 0, par->format_out, glac->G->w_liq, UV, number_novalue);
						if(strcmp(files[rwii] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rwii], 0, par->format_out, glac->G->w_ice, UV, number_novalue);
						if(strcmp(files[rTi] , string_novalue) != 0) write_tensorseries(1, l, isavings, files[rTi], 0, par->format_out, glac->G->T, UV, number_novalue);
					}
					
					if(strcmp(files[rni] , string_novalue) != 0){
					//	M=copydouble_longmatrix(glac->G->lnum);
					//	name = join_strings(files[rni], NNNN);
						name = files[rni]+ string(NNNN);
					//	write_map(name, 1, par->format_out, M, UV, number_novalue);
						write_map(name, 1, par->format_out, glac->G->lnum, UV, number_novalue);
					//	free_doublematrix(M);
					//	free(name);
					}
				}		
				
				if(strcmp(files[rpsich] , string_novalue) != 0){
				//	M=new_doublematrix0_(Nl, par->total_pixel);
				//	initialize_doublematrix(M, 0.0);
					M.resize(Nl+1, par->total_pixel+1, 0.0);

					for (l=0; l<=Nl; l++) {
						for (i=1; i<=par->total_channel; i++){
						//	r = cnet->r->co[i];
							r = cnet->r[i];
						//	c = cnet->c->co[i];
							c = cnet->c[i];
						//	M->co[l][top->j_cont[r][c]] = cnet->SS->P->co[l][i];
							M[l][top->j_cont[r][c]] = cnet->SS->P[l][i];
						}
						//M = cnet->SS->P;

						write_tensorseries_vector(1, l, isavings, files[rpsich], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
					}
				//	free_doublematrix(M);
				}

				if(strcmp(files[rTgch] , string_novalue) != 0){
				//	M=new_doublematrix(Nl, par->total_pixel);
				//	initialize_doublematrix(M, 0.0);
					M.resize(Nl+1, par->total_pixel+1,0.0);
					for (l=1; l<=Nl; l++) {
						for (i=1; i<=par->total_channel; i++){
						//	r = cnet->r->co[i];
							r = cnet->r[i];
						//	c = cnet->c->co[i];
							c = cnet->c[i];
						//	M->co[l][top->j_cont[r][c]] = cnet->SS->T->co[l][i];
							M[l][top->j_cont[r][c]] = cnet->SS->T[l][i];
						}
					//M=cnet->SS->T;
						write_tensorseries_vector(1, l, isavings, files[rTgch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
					}
				//	free_doublematrix(M);
				}
				
				if(strcmp(files[ricegch] , string_novalue) != 0){
				//	M=new_doublematrix(Nl, par->total_pixel);
				//	initialize_doublematrix(M, 0.0);
					M.resize(Nl+1, par->total_pixel+1, 0.0);
					for (l=1; l<=Nl; l++) {
						for (i=1; i<=par->total_channel; i++){
							//	r = cnet->r->co[i];
							r = cnet->r[i];
							//	c = cnet->c->co[i];
							c = cnet->c[i];
						//	M->co[l][top->j_cont[r][c]] = cnet->SS->thi->co[l][i];
							M[l][top->j_cont[r][c]] = cnet->SS->thi[l][i];
						}
					  //M =cnet->SS->thi;
						write_tensorseries_vector(1, l, isavings, files[ricegch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
					}
				//	free_doublematrix(M);
				}
			}
		}
		
	}
	
	if (par->ContRecovery > 0) {
		t_rec += par->Dt;
		printf("t_rec: %f    par->Dt:%f\n", t_rec, par->Dt);

		if (fabs(t_rec - par->ContRecovery*par->Dt)<1.E-5){ //used to be ContRecovery*secinday, replaced with timestep
			t_rec = 0.;
			
			printf("Writing continuous-recovering files\n");
			
			flog = fopen(logfile.c_str(), "a");
			fprintf(flog, "Writing continuous-recovering files\n");
			fclose(flog);
			
			if(strcmp(files[rtime] , string_novalue) != 0){
			//	name = join_strings(files[rtime], textfile);
				name = files[rtime] + string(textfile);
				f = fopen(name.c_str() , "w");
				fprintf(f,"Time[s],Time[d],n,i_run,i_sim,cum_time[s],elapsed_time[s]\n");
				fprintf(f,"%f,%f,%ld,%ld,%ld,%f,%f",times->time+par->Dt+par->delay_day_recover*GTConst::secinday,
									   (times->time+par->Dt)/GTConst::secinday+par->delay_day_recover,
									   (long)(((times->time+par->Dt)/GTConst::secinday+par->delay_day_recover)/par->ContRecovery),
										i_run,i_sim,cum_time,elapsed_time);
				fclose(f);
			//	free(name);
			}
						
			write_suffix(NNNN, 0, 0);
			
			for(l=0;l<=Nl;l++){
				if(strcmp(files[rpsi] , string_novalue) != 0) write_tensorseries_vector(1, l, 0, files[rpsi], 0, par->format_out, sl->SS->P, UV, number_novalue, top->j_cont, Nr, Nc);
				if(l>0){
					if(strcmp(files[riceg] , string_novalue) != 0) write_tensorseries_vector(1, l, 0, files[riceg], 0, par->format_out, sl->SS->thi, UV, number_novalue, top->j_cont, Nr, Nc);
					if(strcmp(files[rTg] , string_novalue) != 0) write_tensorseries_vector(1, l, 0, files[rTg], 0, par->format_out, sl->SS->T, UV, number_novalue, top->j_cont, Nr, Nc);
				}
			}
			
			if(strcmp(files[rwcrn] , string_novalue) != 0){
			//	name = join_strings(files[rwcrn],NNNN);
				name = files[rwcrn] + string(NNNN);
				write_map_vector(name, 0, par->format_out, sl->VS->wrain, UV, number_novalue, top->j_cont, Nr, Nc);
			//	free(name);
			}
			
			if(strcmp(files[rwcsn] , string_novalue) != 0){
			//	name = join_strings(files[rwcsn],NNNN);
				name = files[rwcsn]+ string(NNNN);
				write_map_vector(name, 0, par->format_out, sl->VS->wsnow, UV, number_novalue, top->j_cont, Nr, Nc);
			//	free(name);
			}
			
			for (i=1; i<=par->total_pixel; i++) {
			//	if ((long)sl->VS->Tv->co[i] == number_novalue) sl->VS->Tv->co[i] = 0.;
				if ((long)sl->VS->Tv[i] == number_novalue) sl->VS->Tv[i] = 0.;
			}
			
			if(strcmp(files[rTv] , string_novalue) != 0){
			//	name = join_strings(files[rTv],NNNN);
				name = files[rTv] + string(NNNN);
				write_map_vector(name, 0, par->format_out, sl->VS->Tv, UV, number_novalue, top->j_cont, Nr, Nc);
			//	free(name);
			}	
			
			for(l=1;l<=par->max_snow_layers;l++){
				if(strcmp(files[rDzs] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rDzs], 0, par->format_out, snow->S->Dzl, UV, number_novalue);
				if(strcmp(files[rwls] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rwls], 0, par->format_out, snow->S->w_liq, UV, number_novalue);
				if(strcmp(files[rwis] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rwis], 0, par->format_out, snow->S->w_ice, UV, number_novalue);
				if(strcmp(files[rTs] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rTs], 0, par->format_out, snow->S->T, UV, number_novalue);
			}
			
			if(strcmp(files[rsnag] , string_novalue) != 0){
			//	name = join_strings(files[rsnag],NNNN);
				name = files[rsnag]+ string(NNNN);
				write_map_vector(name, 0, par->format_out, snow->age, UV, number_novalue, top->j_cont, Nr, Nc);
			//	free(name);
			}
			
			if(strcmp(files[rns] , string_novalue) != 0){
			//	M=copydouble_longmatrix(snow->S->lnum);
			//	name = join_strings(files[rns], NNNN);
				name = files[rns]+ string(NNNN);
			//	write_map(name, 1, par->format_out, M, UV, number_novalue);
				write_map(name, 1, par->format_out, snow->S->lnum, UV, number_novalue);
			//	free_doublematrix(M);
			//	free(name);
			}
			
			if(par->max_glac_layers>0){
				for(l=1;l<=par->max_glac_layers;l++){
					if(strcmp(files[rDzi] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rDzi], 0, par->format_out, glac->G->Dzl, UV, number_novalue);
					if(strcmp(files[rwli] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rwli], 0, par->format_out, glac->G->w_liq, UV, number_novalue);
					if(strcmp(files[rwii] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rwii], 0, par->format_out, glac->G->w_ice, UV, number_novalue);
					if(strcmp(files[rTi] , string_novalue) != 0) write_tensorseries(1, l, 0, files[rTi], 0, par->format_out, glac->G->T, UV, number_novalue);
				}
				
				if(strcmp(files[rni] , string_novalue) != 0){
				//	M=copydouble_longmatrix(glac->G->lnum);
				//	name = join_strings(files[rni], NNNN);
					name = files[rni] + string(NNNN);
				//	write_map(name, 1, par->format_out, M, UV, number_novalue);
					write_map(name, 1, par->format_out, glac->G->lnum, UV, number_novalue);
				//	free_doublematrix(M);
				//	free(name);
				}
			}		
			
			
			if(strcmp(files[rpsich] , string_novalue) != 0){
			//	M=new_doublematrix0_(Nl, par->total_pixel);
			//	initialize_doublematrix(M, 0.0);
				M.resize(Nl+1, par->total_pixel+1, 0.0);
				for (l=0; l<=Nl; l++) {
					for (i=1; i<=par->total_channel; i++){
					//	r = cnet->r->co[i];
						r = cnet->r[i];
					//	c = cnet->c->co[i];
						c = cnet->c[i];
					//	M->co[l][top->j_cont[r][c]] = cnet->SS->P->co[l][i];
						M[l][top->j_cont[r][c]] = cnet->SS->P[l][i];
					}
				//M=cnet->SS->P;
					write_tensorseries_vector(1, l, 0, files[rpsich], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
				}
			//	free_doublematrix(M);
			}
			
			if(strcmp(files[rTgch] , string_novalue) != 0){
			//	M=new_doublematrix(Nl, par->total_pixel);
			//	initialize_doublematrix(M, 0.0);
				M.resize(Nl+1, par->total_pixel+1, 0.0);
				for (l=1; l<=Nl; l++) {
					for (i=1; i<=par->total_channel; i++){
					//	r = cnet->r->co[i];
						r = cnet->r[i];
					//	c = cnet->c->co[i];
						c = cnet->c[i];
					//	M->co[l][top->j_cont[r][c]] = cnet->SS->T->co[l][i];
						M[l][top->j_cont[r][c]] = cnet->SS->T[l][i];
					}
				//M= cnet->SS->T;
					write_tensorseries_vector(1, l, 0, files[rTgch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
				}
			//	free_doublematrix(M);
			}
			
			if(strcmp(files[ricegch] , string_novalue) != 0){
			//	M=new_doublematrix(Nl, par->total_pixel);
			//	initialize_doublematrix(M, 0.0);
				M.resize(Nl+1, par->total_pixel+1,0.0);
				for (l=1; l<=Nl; l++) {
					for (i=1; i<=par->total_channel; i++){
					//	M->co[l][top->j_cont[r][c]] = cnet->SS->thi->co[l][i];
						M[l][top->j_cont[r][c]] = cnet->SS->thi[l][i];
					}
					write_tensorseries_vector(1, l, 0, files[ricegch], 0, par->format_out, M, UV, number_novalue, top->j_cont, Nr, Nc);
				}
			//	free_doublematrix(M);
			}
		}
	}
		
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, Energy *egy, SNOW *snow, GLACIER *glac){
  void write_output_headers(long n, Times *times, Water *wat, Par *par, Topo *top, Land *land, Soil *sl, Energy *egy, Snow *snow, Glacier *glac){

	/*internal auxiliary variables:*/
	long i,l,m,j,r,c;
//	char *name,*temp,*temp2,NNNN[ ]={"NNNN"},rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	char NNNN[ ]={"NNNN"},rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	string name, temp, temp2;

	long sy;
	short lu, first_column;
//	DOUBLEVECTOR *root_fraction;
	GeoVector<double> root_fraction;
	FILE *f; 
		
	if (par->recover > 0) write_suffix(rec, par->recover, 4);
	if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);
			
	//DISCHARGE
	if (par->state_discharge == 1 && strcmp(files[fQ] , string_novalue) != 0){
		if (par->recover>0) {
		//	temp = join_strings(files[fQ], rec);
			temp = files[fQ] + string(rec);
		//	name = join_strings(temp, textfile);
			name = temp + textfile;
		//	free(temp);
		}else if (par->n_ContRecovery>0) {
			temp = join_strings(files[fQ], crec);
		//	name = join_strings(temp, textfile);
			name = temp + textfile;
		//	free(temp);
		}else {
		//	name = join_strings(files[fQ], textfile);
			name = files[fQ] + string(textfile) ;
		}
		
		f=t_fopen(name.c_str(),"w");
		fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,Qtot[m3/s],Vsup/Dt[m3/s],Vsub/Dt[m3/s],Vchannel[m3],Qoutlandsup[m3/s],Qoutlandsub[m3/s],Qoutbottom[m3/s]\n");
		t_fclose(f);
	//	free(name);
	}
	
	
	if(par->state_pixel == 1){
		
		//output matrix and vectors
		m=(long)otot;
		odpnt=(double**)malloc(m*sizeof(double*));
		odp=(double**)malloc(m*sizeof(double*));
		for (i=0; i<otot; i++) {
		//	odpnt[i]=(double*)malloc(par->rc->nrh*sizeof(double));
			odpnt[i]=(double*)malloc(par->rc.getRows()*sizeof(double));
		//	odp[i]=(double*)malloc(par->rc->nrh*sizeof(double));
			odp[i]=(double*)malloc(par->rc.getRows()*sizeof(double));
		//	for (j=0; j<par->rc->nrh; j++) {
			for (j=0; j<par->rc.getRows()-1; j++) {
				odpnt[i][j] = 0.;
				odp[i][j] = 0.;
			}
		}
						
		if(strcmp(files[fpointwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
			//	temp = join_strings(files[fpointwriteend], rec);
				temp = files[fpointwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fpointwriteend], crec);
				temp = files[fpointwriteend]+ string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fpointwriteend], textfile);
				name = files[fpointwriteend] + string(textfile);
			}
			
			ffpoint=fopen(name.c_str(),"w");
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
		//	free(name);
		}
		
		if(strcmp(files[fsnzwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
				temp = join_strings(files[fsnzwriteend], rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
				temp = join_strings(files[fsnzwriteend], crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fsnzwriteend], textfile);
				name = files[fsnzwriteend] + string(textfile) ;
			}
			
			ffsnow=fopen(name.c_str(),"w");
			
		//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
			if ((long)par->snow_plot_depths[1] != number_novalue) {
			//	m = par->snow_plot_depths->nh;
				m = par->snow_plot_depths.size();
			}else {
				m = par->max_snow_layers;
			}

			first_column=1;
			for(j=0;j<nosnw;j++){
				if(first_column==0){
					fprintf(ffsnow,",");
				}else {
					first_column = 0;
				}			
				if (osnw[j] >= 0 && osnw[j]<=5) {
					fprintf(ffsnow,"%s",hsnw[osnw[j]]);					
				}else if (osnw[j] >= 6 && osnw[j] < 6 + 3*m) {
					l = (long)fmod( (double)osnw[j]-6., (double)m ) + 1;
					n = floor( ( (double)osnw[j]-6.) / (double)m ) + 6;
				//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
					if ((long)par->snow_plot_depths[1] != number_novalue) {
					//	fprintf(ffsnow, "%s(%f)",hsnw[n],par->snow_plot_depths->co[l]);
						fprintf(ffsnow, "%s(%f)",hsnw[n],par->snow_plot_depths[l]);
					}else {
						fprintf(ffsnow, "%s(%ld)",hsnw[n],l);
					}
				}else if (osnw[j] >= 6 + 3*m) {
					l = (long)fmod( (double)osnw[j]-6.-3*(double)m, (double)par->max_snow_layers ) + 1;
					n = floor( ( (double)osnw[j]-6.-3.*(double)m) / (double)par->max_snow_layers ) + 6 + 3;
					fprintf(ffsnow, "%s(%ld)",hsnw[n],l);
				}else{
					fprintf(ffsnow, "None");
				}
			}
			fprintf(ffsnow,"\n");
		//	free(name);
		}
		
		if(par->max_glac_layers>0){
			if(strcmp(files[fglzwriteend] , string_novalue) != 0){
				
				if (par->recover>0) {
				//	temp = join_strings(files[fpointwriteend], rec);
					temp = files[fpointwriteend] + string(rec);
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				//	free(temp);
				}else if (par->n_ContRecovery>0) {
				//	temp = join_strings(files[fpointwriteend], crec);
					temp = files[fpointwriteend] + string(crec);
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				//	free(temp);
				}else {
				//	name = join_strings(files[fpointwriteend], textfile);
					name = files[fpointwriteend] + string(textfile);
				}
				
				ffglac=t_fopen(name.c_str(),"w");

			//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
				if ((long)par->glac_plot_depths[1] != number_novalue) {
				//	m = par->glac_plot_depths->nh;
					m = par->glac_plot_depths.size();
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
					//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
						if ((long)par->glac_plot_depths[1] != number_novalue) {
						//	fprintf(ffglac, "%s(%f)",hglc[n],par->glac_plot_depths->co[l]);
							fprintf(ffglac, "%s(%f)",hglc[n],par->glac_plot_depths[l]);
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
			//	free(name);
			}
		}
		
		if(strcmp(files[fTzwriteend] , string_novalue) != 0){
			
			if (par->recover>0) {
			//	temp = join_strings(files[fTzwriteend], rec);
				temp = files[fTzwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fTzwriteend], crec);
				temp = files[fTzwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fTzwriteend], textfile);
				name = files[fTzwriteend] + string(textfile);

			}
			
			ffT=fopen(name.c_str(),"w");
			write_soil_header(ffT, par->soil_plot_depths, sl->pa);
		//	free(name);
		}
		
		if(strcmp(files[fTzavwriteend] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[fTzavwriteend], rec);
				temp = files[fTzavwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fTzavwriteend], crec);
				temp = files[fTzavwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fTzavwriteend], textfile);
				name = files[fTzavwriteend]+ string(textfile);
			}
			
			ffTav=fopen(name.c_str(),"w");
			write_soil_header(ffTav, par->soil_plot_depths, sl->pa);
		//	free(name);
		}
		
		if(strcmp(files[fpsiztotwriteend] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[fpsiztotwriteend], rec);
				temp = files[fpsiztotwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fpsiztotwriteend], crec);
				temp = files[fpsiztotwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fpsiztotwriteend], textfile);
				name = files[fpsiztotwriteend]+ string(textfile);
			}
			
			ffpsitot=fopen(name.c_str(),"w");
			write_soil_header(ffpsitot, par->soil_plot_depths, sl->pa);
		//	free(name);
		}
		
		if(strcmp(files[fpsizwriteend] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[fpsizwriteend], rec);
				temp = files[fpsizwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fpsizwriteend], crec);
				temp = files[fpsizwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fpsizwriteend], textfile);
				name = files[fpsizwriteend]+ string(textfile);
			}

			ffpsi=fopen(name.c_str(),"w");
			write_soil_header(ffpsi, par->soil_plot_depths, sl->pa);
		//	free(name);
		}
		
		if(strcmp(files[fliqzwriteend] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[fliqzwriteend], rec);
				temp = files[fliqzwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fliqzwriteend], crec);
				temp = files[fliqzwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fliqzwriteend], textfile);
				name = files[fliqzwriteend] + string(textfile);
			}
			
			ffliq=fopen(name.c_str(),"w");
			write_soil_header(ffliq, par->soil_plot_depths, sl->pa);
		//	free(name);
		}

		if(strcmp(files[fliqzavwriteend] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[fliqzavwriteend], rec);
				temp = files[fliqzavwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fliqzavwriteend], crec);
				temp = files[fliqzavwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fliqzavwriteend], textfile);
				name = files[fliqzavwriteend]+ string(textfile);
			}

			ffliqav=fopen(name.c_str(),"w");
			write_soil_header(ffliqav, par->soil_plot_depths, sl->pa);
		//	free(name);
		}

		if(strcmp(files[ficezwriteend] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[ficezwriteend], rec);
				temp = files[ficezwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[ficezwriteend], crec);
				temp = files[ficezwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[ficezwriteend], textfile);
				name = files[ficezwriteend] + string(textfile);
			}

			ffice=fopen(name.c_str(),"w");
			write_soil_header(ffice, par->soil_plot_depths, sl->pa);
		//	free(name);
		}
		
		if(strcmp(files[ficezavwriteend] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[ficezavwriteend], rec);
				temp = files[ficezavwriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[ficezavwriteend], crec);
				temp = files[ficezavwriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[ficezavwriteend], textfile);
				name = files[ficezavwriteend]+ string(textfile);
			}
			
			fficeav=fopen(name.c_str(),"w");
			write_soil_header(fficeav, par->soil_plot_depths, sl->pa);
		//	free(name);
		}
		
	//	root_fraction=new_doublevector(Nl);
		root_fraction.resize(Nl+1);
   
		
		//DATA POINTS
	//	for(i=1;i<=par->rc->nrh;i++){
		for(i=1;i<par->rc.getRows();i++){
						
		//	write_suffix(NNNN, par->IDpoint->co[i], 0);
			write_suffix(NNNN, par->IDpoint[i], 0);
		//	r=par->rc->co[i][1];
			r=par->rc[i][1];
		//	c=par->rc->co[i][2];
			c=par->rc[i][2];
		//	sy=sl->type->co[r][c];
			sy=sl->type[r][c];
		//	lu=(short)land->LC->co[r][c];
			lu=(short)land->LC[r][c];
						
			if(strcmp(files[fpoint] , string_novalue) != 0 && par->point_sim != 1){			
			//	name=join_strings(files[fpoint],"_info_");
				name=files[fpoint]+string("_info_");
			//	temp=join_strings(name,NNNN);
				temp = name + NNNN;
			//	temp2 = join_strings(temp,textfile);
				temp2 = temp + textfile;
				f=t_fopen(temp2.c_str(),"w");
				
			//	fprintf(f," The main properties of the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld are:\n",par->chkpt->co[i][ptX],par->chkpt->co[i][ptY],r,c);
				fprintf(f," The main properties of the pixel E=%15.3f N=%15.3f, row=%4ld col=%4ld are:\n",par->chkpt[i][ptX],par->chkpt[i][ptY],r,c);
			//	fprintf(f," Elevation above sea level: %10.3f m\n",top->Z0->co[r][c]);
				fprintf(f," Elevation above sea level: %10.3f m\n",top->Z0[r][c]);
			//	fprintf(f," Gauckler-Strickler [m^1/3/s]: %f\n",land->ty->co[lu][jcm]);
				fprintf(f," Gauckler-Strickler [m^1/3/s]: %f\n",land->ty[lu][jcm]);
				for(l=1;l<=Nl;l++){
					fprintf(f," Residual water content[-] of the layer %ld: %f\n",l,sl->pa[sy][jres][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Saturated water content[-] of the layer %ld: %f\n",l,sl->pa[sy][jsat][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Alpha of van Genuchten[mm^-1] of the layer %ld: %f\n",l,sl->pa[sy][ja][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," n of van Genuchten[-] of the layer %ld: %f\n",l,sl->pa[sy][jns][l]);
				}	
				for(l=1;l<=Nl;l++){
					fprintf(f," m of van Genuchten[-] of the layer %ld: %f\n",l,1-1/sl->pa[sy][jns][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," v of van Genuchten[-] of the layer %ld: %f\n",l,sl->pa[sy][jv][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Water content of wilting point [-] of the layer %ld: %f\n",l,sl->pa[sy][jwp][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Water content of field capacity [-] of the layer %ld: %f\n",l,sl->pa[sy][jfc][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Kv_sat of layer %ld [mm/s]: %f\n",l,sl->pa[sy][jKn][l]);
				}
				for(l=1;l<=Nl;l++){
					fprintf(f," Kh_sat of layer %ld [mm/s]: %f\n",l,sl->pa[sy][jKl][l]);
				}

			//	fprintf(f," Terrain elevation [m]: %f\n",top->Z0->co[r][c]);
				fprintf(f," Terrain elevation [m]: %f\n",top->Z0[r][c]);
			//	fprintf(f," Sky view factor [-]: %f\n",top->sky->co[r][c]);
				fprintf(f," Sky view factor [-]: %f\n",top->sky[r][c]);
			//	fprintf(f," The pixel-type is %d \n",top->pixel_type->co[r][c]);
				fprintf(f," The pixel-type is %d \n",top->pixel_type[r][c]);
			//	fprintf(f," Aspect [deg] [0=Nord, clockwise]: %f \n",top->aspect->co[r][c]);
				fprintf(f," Aspect [deg] [0=Nord, clockwise]: %f \n",top->aspect[r][c]);
			//	fprintf(f," Mean slope of the pixel [deg]: %f \n",top->slope->co[r][c]);
				fprintf(f," Mean slope of the pixel [deg]: %f \n",top->slope[r][c]);
			//	fprintf(f," Land use number is %d \n",(short)land->LC->co[r][c]);
				fprintf(f," Land use number is %d \n",(short)land->LC[r][c]);

				for(l=1;l<=Nl;l++){
				//	fprintf(f," The root fraction [-] of layer %ld: %f\n",l,land->root_fraction->co[lu][l]);
					fprintf(f," The root fraction [-] of layer %ld: %f\n",l,land->root_fraction[lu][l]);
				}
		
			//	fprintf(f," Surface fraction of land covered by vegetation [-]: %f \n",land->ty->co[lu][jcf]);
				fprintf(f," Surface fraction of land covered by vegetation [-]: %f \n",land->ty[lu][jcf]);
			//	fprintf(f," Leaf and Stem Area Index [-]: %f \n",land->ty->co[lu][jLSAI]);
				fprintf(f," Leaf and Stem Area Index [-]: %f \n",land->ty[lu][jLSAI]);
			//	fprintf(f," Momentum roughness length z0soil [m]: %f \n",land->ty->co[lu][jz0]);
				fprintf(f," Momentum roughness length z0soil [m]: %f \n",land->ty[lu][jz0]);
			//	fprintf(f," Vegetation height [m]: %f \n",land->ty->co[lu][jHveg]);
				fprintf(f," Vegetation height [m]: %f \n",land->ty[lu][jHveg]);

				fprintf(f," \n");
				t_fclose(f);
			//	free(temp2);
			//	free(temp);
			//	free(name);
			}
			
			if(strcmp(files[fpoint] , string_novalue) != 0){
				
			//	temp=join_strings(files[fpoint],NNNN);
				temp=files[fpoint] + string(NNNN);

				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				}
				
			//	free(temp);
				
				f=t_fopen(name.c_str(),"w");
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
			//	free(name);
			}
			
			if(strcmp(files[fsnz] , string_novalue) != 0){
				
			//	temp=join_strings(files[fsnz],NNNN);
				temp= files[fsnz] + string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				}
				
			//	free(temp);
				
			//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
				if ((long)par->snow_plot_depths[1] != number_novalue) {
				//	m = par->snow_plot_depths->nh;
					m = par->snow_plot_depths.size();
				}else {
					m = par->max_snow_layers;
				}

				f=t_fopen(name.c_str(),"w");
				first_column=1;
				for(j=0;j<nosnw;j++){
					if(first_column==0){
						fprintf(f,",");
					}else {
						first_column = 0;
					}			
					if (osnw[j] >= 0 && osnw[j]<=5) {
						fprintf(f,"%s",hsnw[osnw[j]]);					
					}else if (osnw[j] >= 6 && osnw[j] < 6+3*m) {
						l = (long)fmod( (double)osnw[j]-6., (double)m ) + 1;
						n = floor( ( (double)osnw[j]-6.) / (double)m ) + 6;
					//	if ((long)par->snow_plot_depths->co[1] != number_novalue) {
						if ((long)par->snow_plot_depths[1] != number_novalue) {
						//	fprintf(f, "%s(%f)",hsnw[n],par->snow_plot_depths->co[l]);
							fprintf(f, "%s(%f)",hsnw[n],par->snow_plot_depths[l]);
						}else {
							fprintf(f, "%s(%ld)",hsnw[n],l);
						}
					}else if (osnw[j] >=6+3*m) {
						l = (long)fmod( (double)osnw[j]-6.-3*(double)m, (double)par->max_snow_layers ) + 1;
						n = floor( ( (double)osnw[j]-6.-3.*(double)m) / (double)par->max_snow_layers ) + 6 + 3;
						fprintf(f, "%s(%ld)",hsnw[n],l);
					}else{
						fprintf(f, "None");
					}
				}
				fprintf(f,"\n");
				t_fclose(f);
			//	free(name);
			}
			
			if(par->max_glac_layers>0){
				if(strcmp(files[fglz] , string_novalue) != 0){
					
				//	temp=join_strings(files[fglz],NNNN);
					temp= files[fglz]+string(NNNN);
					
					if (par->recover>0) {
					//	temp2 = join_strings(temp, rec);
						temp2 = temp + rec;
					//	name = join_strings(temp2, textfile);
						name = temp2 + textfile;
					//	free(temp2);
					}else if (par->n_ContRecovery>0) {
					//	temp2 = join_strings(temp, crec);
						temp2 = temp + crec;
					//	name = join_strings(temp2, textfile);
						name = temp2 + textfile ;
					//	free(temp2);
					}else {
					//	name = join_strings(temp, textfile);
						name = temp + textfile ;
					}
					
				//	free(temp);
					
				//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
					if ((long)par->glac_plot_depths[1] != number_novalue) {
					//	m = par->glac_plot_depths->nh;
						m = par->glac_plot_depths.size();
					}else {
						m = par->max_glac_layers;
					}

					f=t_fopen(name.c_str(),"w");
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
						//	if ((long)par->glac_plot_depths->co[1] != number_novalue) {
							if ((long)par->glac_plot_depths[1] != number_novalue) {
							//	fprintf(f, "%s(%f)",hglc[n],par->glac_plot_depths->co[l]);
								fprintf(f, "%s(%f)",hglc[n],par->glac_plot_depths[l]);
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
				//	free(name);
				}
			}
				
			if(strcmp(files[fTz] , string_novalue) != 0){
				
			//	temp=join_strings(files[fTz],NNNN);
				temp=files[fTz] +string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				}
				
			//	free(temp);
				
				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);	
			//	free(name);
			}
			
			if(strcmp(files[fTzav] , string_novalue) != 0){

			//	temp=join_strings(files[fTzav],NNNN);
				temp= files[fTzav] +string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile ;
				}
				
			//	free(temp);
				
				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);	
			//	free(name);
			}
		
			if(strcmp(files[fpsiztot] , string_novalue) != 0){

			//	temp=join_strings(files[fpsiztot],NNNN);
				temp= files[fpsiztot] +string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				}
				
			//	free(temp);

				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);
			//	free(name);
			}
			
			if(strcmp(files[fpsiz] , string_novalue) != 0){

			//	temp=join_strings(files[fpsiz],NNNN);
				temp= files[fpsiz] + string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile ;
				}
				
			//	free(temp);

				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);
			//	free(name);
			}
		
			if(strcmp(files[fliqz] , string_novalue) != 0){

			//	temp=join_strings(files[fliqz],NNNN);
				temp= files[fliqz] +string(NNNN);
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				}
				
			//	free(temp);

				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);
			//	free(name);
			}
			
			if(strcmp(files[fliqzav] , string_novalue) != 0){

			//	temp=join_strings(files[fliqzav],NNNN);
				temp= files[fliqzav] +string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile ;
				}
				
			//	free(temp);

				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);	
			//	free(name);
			}
		
			if(strcmp(files[ficez] , string_novalue) != 0){

			//	temp=join_strings(files[ficez],NNNN);
				temp= files[ficez] + string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile ;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile ;
				}
				
			//	free(temp);

				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);	
			//	free(name);
			}
							
			if(strcmp(files[ficezav] , string_novalue) != 0){

			//	temp=join_strings(files[ficezav],NNNN);
				temp= files[ficezav] + string(NNNN);
				
				if (par->recover>0) {
				//	temp2 = join_strings(temp, rec);
					temp2 = temp + rec ;
				//	name = join_strings(temp2, textfile);
					name = temp2 + textfile;
				//	free(temp2);
				}else if (par->n_ContRecovery>0) {
				//	temp2 = join_strings(temp, crec);
					temp2 = temp + crec;
				//	name = join_strings(temp2, textfile);
					name = temp2+ textfile;
				//	free(temp2);
				}else {
				//	name = join_strings(temp, textfile);
					name = temp + textfile;
				}
				
			//	free(temp);

				f=t_fopen(name.c_str(),"w");
				write_soil_header(f, par->soil_plot_depths, sl->pa);
				t_fclose(f);	
			//	free(name);
			}
						
		}
		
	//	free_doublevector(root_fraction);
		
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
			//	temp = join_strings(files[fbaswriteend], rec);
				temp = files[fbaswriteend] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile ;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fbaswriteend], crec);
				temp = files[fbaswriteend] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile;
			//	free(temp);
			}else {
			//	name = join_strings(files[fbaswriteend], textfile);
				name = files[fbaswriteend]+ string(textfile);
			}
			
			ffbas=fopen(name.c_str(),"w");
			
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
		//	free(name);
		}
		
		if(strcmp(files[fbas] , string_novalue) != 0){

			if (par->recover>0) {
			//	temp = join_strings(files[fbas], rec);
				temp = files[fbas] + string(rec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile ;
			//	free(temp);
			}else if (par->n_ContRecovery>0) {
			//	temp = join_strings(files[fbas], crec);
				temp = files[fbas] + string(crec);
			//	name = join_strings(temp, textfile);
				name = temp + textfile ;
			//	free(temp);
			}else {
			//	name = join_strings(files[fbas], textfile);
				name = files[fbas] + string(textfile);
			}

			f=t_fopen(name.c_str(),"w");
			
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
		//	free(name);
		}
				
	}
	
	//SNOW COVERED AREA STATISTICS
	if(par->point_sim!=1 && strcmp(files[fSCA] , string_novalue) != 0){
		
		if (par->recover>0) {
		//	temp = join_strings(files[fSCA], rec);
			temp = files[fSCA] + string(rec);
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		//	free(temp);
		}else if (par->n_ContRecovery>0) {
		//	temp = join_strings(files[fSCA], crec);
			temp = files[fSCA] + string(crec);
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		//	free(temp);
		}else {
		//	name = join_strings(files[fSCA], textfile);
			name = files[fSCA] + string(textfile);
		}
		
		f=t_fopen(name.c_str(),"w");
		fprintf(f,"DATE[day/month/year hour:min],t[days],JDfrom0,JD,snowDav,SWEav,Tav,Tsav,perc.SFA,perc.SCA\n");
		t_fclose(f);
	//	free(name);
	}
	
	

}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//	void write_soil_output(long i, long iname, double init_date, double JDfrom0, double JD, long day, long month, long year, long hour, long minute, DOUBLEVECTOR *n, Soil *sl, Par *par, double psimin, double cosslope){
  	void write_soil_output(long i, long iname, double init_date, double JDfrom0, double JD, long day, long month, long year, long hour, long minute, const GeoVector<double>& n, Soil *sl, Par *par, double psimin, double cosslope){

//  char *name,*temp,*temp2,NNNN[ ]={"NNNN"};
  	char NNNN[ ]={"NNNN"};
  	string name, temp, temp2;

	char rec[ ]={"_recNNNN"},crec[ ]={"_crecNNNN"};
	long l;
	FILE *f;
	
	write_suffix(NNNN, iname, 0);	
	
	if (par->recover > 0) write_suffix(rec, par->recover, 4);
	if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

	if(strcmp(files[fTz] , string_novalue) != 0){
		
	//	temp=join_strings(files[fTz],NNNN);
		temp= files[fTz] + string(NNNN);

		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		}
		
		f=fopen(name.c_str(),"a");
	//	write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}

	if(strcmp(files[fTzwriteend] , string_novalue) != 0){
	//	write_soil_file(1, iname, ffT, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, ffT, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzplot, i, n, sl->pa, cosslope);
	}
	
	if(strcmp(files[fTzav] , string_novalue) != 0){
	//	temp=join_strings(files[fTzav],NNNN);
		temp=files[fTzav] + string(NNNN);
		
		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		}		
		
		f=fopen(name.c_str() ,"a");
	//	write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzavplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzavplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}

	if(strcmp(files[fTzavwriteend] , string_novalue) != 0){
	//	write_soil_file(1, iname, ffTav, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzavplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, ffTav, day, month, year, hour, minute, JDfrom0, init_date, sl->Tzavplot, i, n, sl->pa, cosslope);
	}
	
	if(strcmp(files[fpsiztot] , string_novalue) != 0){
	//	temp=join_strings(files[fpsiztot],NNNN);
		temp= files[fpsiztot] + string(NNNN);
		
		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile;
		}		
		
		f=fopen(name.c_str(),"a");
	//	write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Ptotzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Ptotzplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}
	
	if(strcmp(files[fpsiztotwriteend] , string_novalue) != 0){
	//	write_soil_file(1, iname, ffpsitot, day, month, year, hour, minute, JDfrom0, init_date, sl->Ptotzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, ffpsitot, day, month, year, hour, minute, JDfrom0, init_date, sl->Ptotzplot, i, n, sl->pa, cosslope);
	}
	
	if(strcmp(files[fpsiz] , string_novalue) != 0){
	//	temp=join_strings(files[fpsiz],NNNN);
		temp= files[fpsiz] + string(NNNN);
		
		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		}		
		
		f=fopen(name.c_str(),"a");
	//	write_soil_file(0, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Pzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(0, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->Pzplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}

	if(strcmp(files[fpsizwriteend] , string_novalue) != 0){
	//	write_soil_file(0, iname, ffpsi, day, month, year, hour, minute, JDfrom0, init_date, sl->Pzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(0, iname, ffpsi, day, month, year, hour, minute, JDfrom0, init_date, sl->Pzplot, i, n, sl->pa, cosslope);
	}
	
	if(strcmp(files[fliqz] , string_novalue) != 0){
	//	temp=join_strings(files[fliqz],NNNN);
		temp=files[fliqz] + string(NNNN);
		
		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		}		
		
		f=fopen(name.c_str(),"a");
	//	write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thzplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}

	if(strcmp(files[fliqzwriteend] , string_novalue) != 0){
	//	write_soil_file(1, iname, ffliq, day, month, year, hour, minute, JDfrom0, init_date, sl->thzplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, ffliq, day, month, year, hour, minute, JDfrom0, init_date, sl->thzplot, i, n, sl->pa, cosslope);
	}

	if(strcmp(files[fliqzav] , string_novalue) != 0){
	//	temp=join_strings(files[fliqzav],NNNN);
		temp=files[fliqzav] + string(NNNN);
		
		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		}		
		
		f=fopen(name.c_str(),"a");
	//	write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thzavplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thzavplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}
	
	if(strcmp(files[fliqzavwriteend] , string_novalue) != 0){
	//	write_soil_file(1, iname, ffliqav, day, month, year, hour, minute, JDfrom0, init_date, sl->thzavplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, ffliqav, day, month, year, hour, minute, JDfrom0, init_date, sl->thzavplot, i, n, sl->pa, cosslope);
	}
	
	if(strcmp(files[ficez] , string_novalue) != 0){
	//	temp=join_strings(files[ficez],NNNN);
		temp=files[ficez] + string(NNNN);
		
		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		}		
		
		f=fopen(name.c_str(),"a");
	//	write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thizplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thizplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}

	if(strcmp(files[ficezwriteend] , string_novalue) != 0){
	//	write_soil_file(1, iname, ffice, day, month, year, hour, minute, JDfrom0, init_date, sl->thizplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, ffice, day, month, year, hour, minute, JDfrom0, init_date, sl->thizplot, i, n, sl->pa, cosslope);
	}

	if(strcmp(files[ficezav] , string_novalue) != 0){
	//	temp=join_strings(files[ficezav],NNNN);
		temp=files[ficezav] + string(NNNN);
		
		if (par->recover>0) {
		//	temp2 = join_strings(temp, rec);
			temp2 = temp + rec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else if (par->n_ContRecovery>0) {
		//	temp2 = join_strings(temp, crec);
			temp2 = temp + crec ;
		//	name = join_strings(temp2, textfile);
			name = temp2 + textfile ;
		//	free(temp2);
		}else {
		//	name = join_strings(temp, textfile);
			name = temp + textfile ;
		}		
		
		f=fopen(name.c_str(),"a");
	//	write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thizavplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, f, day, month, year, hour, minute, JDfrom0, init_date, sl->thizavplot, i, n, sl->pa, cosslope);
		fclose(f);
	//	free(name);
	//	free(temp);
	}
	
	if(strcmp(files[ficezavwriteend] , string_novalue) != 0){
	//	write_soil_file(1, iname, fficeav, day, month, year, hour, minute, JDfrom0, init_date, sl->thizavplot->co[i], n, sl->pa[1][jdz], cosslope);
		write_soil_file(1, iname, fficeav, day, month, year, hour, minute, JDfrom0, init_date, sl->thizavplot, i, n, sl->pa, cosslope);
	}
	
	for(l=1;l<=Nl;l++){
	//	if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot->co[i][l] = 0.0;
		if(strcmp(files[fTzav] , string_novalue) != 0 || strcmp(files[fTzavwriteend] , string_novalue) != 0) sl->Tzavplot[i][l] = 0.0;
	//	if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot->co[i][l] = 0.0;
		if(strcmp(files[fliqzav] , string_novalue) != 0 || strcmp(files[fliqzavwriteend] , string_novalue) != 0) sl->thzavplot[i][l] = 0.0;
	//	if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot->co[i][l] = 0.0;
		if(strcmp(files[ficezav] , string_novalue) != 0 || strcmp(files[ficezavwriteend] , string_novalue) != 0) sl->thizavplot[i][l] = 0.0;
	}
}


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************


//	void write_soil_file(long lmin, long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, double *var, DOUBLEVECTOR *n, double *dz, double cosslope){
void write_soil_file(long lmin, long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, const GeoMatrix<double>& var, long row, const GeoVector<double>& n, const GeoTensor<double>& dz, double cosslope){

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
				fprintf(f, "%f",JDfrom0-JDfrom0init);
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

//	if ((long)n->co[1] != number_novalue) {
	if ((long)n[1] != number_novalue) {
	//	for (l=1; l<=n->nh; l++) {
		for (l=1; l<=n.size(); l++) {
		//	fprintf(f, ",%f",interpolate_soil(lmin, n->co[l]*cosslope, Nl, dz, var));
			fprintf(f, ",%f",interpolate_soil(lmin, n[l]*cosslope, Nl, dz, var, row));
		}
	}else{
		for(l=1;l<=Nl;l++){
			fprintf(f,",%f", var(row,l));
		}
	}
	
	fprintf(f," \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

// void write_soil_header(FILE *f, DOUBLEVECTOR *n, double *dz){
void write_soil_header(FILE *f, const GeoVector<double>& n, const GeoTensor<double>& dz){
	
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

//	if ((long)n->co[1] != number_novalue) {
	if ((long)n[1] != number_novalue) {
	//	for (l=1; l<=n->nh; l++) {
		for (l=1; l<=n.size(); l++) {
		//	fprintf(f, ",%f",n->co[l]);
			fprintf(f, ",%f",n[l]);
		}
	}else{
		for(l=1;l<=Nl;l++){
			//z += dz[l];
			z += dz(1,jdz,l);
			//fprintf(f,",%f ",z-0.5*dz[l]);
			fprintf(f,",%f ",z-0.5*dz(1,jdz,l));
		}
	}
	
	fprintf(f," \n");
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//void plot(char *name, long i_plot, DOUBLEVECTOR *V, short format, long **J){
void plot(char *name, long i_plot, const GeoVector<double>& V, short format, long **J){
	
	char ADS[ ]={"iiii"};
//	char *temp;
	string  temp;
	
	write_suffix(ADS, i_plot, 0);
	temp=join_strings(name,ADS);

	write_map_vector(temp, 0, format, V, UV, number_novalue, J, Nr, Nc);
//	free(temp);
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

double interpolate_soil(long lmin, double h, long max, const GeoTensor<double>& Dz, const GeoMatrix<double>& Q, const long& row){
	
	double q, z, z0=0.;
	long l;
		
	l = lmin;	
	q = (double)number_novalue;
	
	do{
		
		if (l == lmin){
			z = z0;
			//if (l>0) z += Dz[l]/2.;
			if (l>0) z += Dz(1,jdz,l)/2.;
		}else if (l <= max) {
			//z = z0 + Dz[l]/2;
			//if (l>1) z += Dz[l-1]/2.;
			z = z0 + Dz(1,jdz,l)/2;
			if (l>1) z += Dz(1,jdz,l-1)/2.;
		}else {
			//z = z0 + Dz[max]/2.;
			z = z0 + Dz(1,jdz,max)/2.;
		}
		
		if(h < z && h >= z0){
			if (l == lmin) {
				q = Q(row, lmin);
			}else if (l <= max) {
				q = ( Q(row,l-1) * (z-h) + Q(row, l) * (h-z0) ) / (z - z0);
			}else {
				q = Q(row, max);
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

//double interpolate_soil2(long lmin, double h, long max, double *Dz, DOUBLEMATRIX *Q, long i){
double interpolate_soil2(long lmin, double h, long max, const GeoTensor<double>& Dz, GeoMatrix<double>& Q, long i){

	double q, z, z0=0.;
	long l;

	l = lmin;
	q = (double)number_novalue;

	do{

		if (l == lmin){
			z = z0;
			//if (l>0) z += Dz[1][jdz][l]/2.;
			if (l>0) z += Dz(1,jdz,l)/2.;
		}else if (l <= max) {
			//z = z0 + Dz[l]/2.;
			z = z0 + Dz(1,jdz,l)/2.;
			//if (l>1) z += Dz[1][jdz][l-1]/2.;
			if (l>1) z += Dz(1,jdz,l-1)/2.;
		}else {
			//z = z0 + Dz[1][jdz][max]/2.;
			z = z0 + Dz(1,jdz,max)/2.;
		}

		if(h < z && h >= z0){
			if (l == lmin) {
			//	q = Q->co[l][i];
				q = Q[l][i];
			}else if (l <= max) {
			//	q = ( Q->co[l-1][i] * (z-h) + Q->co[l][i] * (h-z0) ) / (z - z0);
				q = ( Q[l-1][i] * (z-h) + Q[l][i] * (h-z0) ) / (z - z0);
			}else {
			//	q = Q->co[max][i];
				q = Q[max][i];
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

//void write_tensorseries_soil(long lmin, char *suf, char *filename, short type, short format, DOUBLEMATRIX *T, DOUBLEVECTOR *n, long **J, LONGMATRIX *RC, double *dz, DOUBLEMATRIX *slope, short vertical){
void write_tensorseries_soil(long lmin, char *suf, char *filename, short type, short format, GeoMatrix<double>& T, const GeoVector<double>& n, long **J, GeoMatrix<long>& RC, GeoTensor<double>& dz, GeoMatrix<double>& slope, short vertical){

 	char LLLLL[ ]={"LLLLL"};
 //	char *temp1, *temp2;
 	string temp1,  temp2;
 //	long i, l, npoints=T->nch;
	long npoints=T.getCols();
 	double cosslope=1.;
 //	DOUBLEVECTOR *V;
 	GeoVector<double> V(npoints);

 //	for(l=1; l<=n->nh; l++){
 	for(size_t l=1; l<n.size(); l++){
 		temp1 = join_strings(LLLLL, suf);
 		write_suffix(temp1, l, 1);

		//for(i=1; i<=npoints; i++){
 		for(long i=1; i<npoints; i++){
 		//	if (vertical == 1) cosslope = cos( Fmin(GTConst::max_slope, slope->co[RC->co[i][1]][RC->co[i][2]]) * GTConst::Pi/180. );
 			if (vertical == 1) cosslope = cos( Fmin(GTConst::max_slope, slope[RC[i][1]][RC[i][2]]) * GTConst::Pi/180. );
 		//	V->co[i] = interpolate_soil2(lmin, n->co[l]*cosslope, Nl, dz, T, i);
 			V[i] = interpolate_soil2(lmin, n[l]*cosslope, Nl, dz, T, i);
 		}

 	//	temp2 = join_strings(filename, temp1);
 		temp2 = filename + temp1;
 	//	write_map_vector(temp2, type, format, V, UV, number_novalue, J, slope->nrh, slope->nch);
 		write_map_vector(temp2, type, format, V, UV, number_novalue, J, slope.getRows()-1, slope.getCols()-1);

 	//	free_doublevector(V);
 	//	free(temp1);
 	//	free(temp2);
 	}
 }
//overloaded function

void write_tensorseries_soil(long lmin, string suf, char *filename, short type, short format, GeoMatrix<double>& T, const GeoVector<double>& n, long **J, GeoMatrix<long>& RC, GeoTensor<double>& dz, GeoMatrix<double>& slope, short vertical){

 	char LLLLL[ ]={"LLLLL"};
 //	char *temp1, *temp2;
	string temp1, temp2;
 //	long i, l, npoints=T->nch;
	long i, l, npoints=T.getCols();
 	double cosslope=1.;
 //	DOUBLEVECTOR *V;
 	GeoVector<double> V(npoints);

 //	for(l=1; l<=n->nh; l++){
 	for(l=1; l<n.size(); l++){
 	//	temp1 = join_strings(LLLLL, suf);
 		temp1 = LLLLL + suf;
 		write_suffix(temp1, l, 1);

		//for(i=1; i<=npoints; i++){
 		for(i=1; i<npoints; i++){
 		//	if (vertical == 1) cosslope = cos( Fmin(GTConst::max_slope, slope->co[RC->co[i][1]][RC->co[i][2]]) * GTConst::Pi/180. );
 			if (vertical == 1) cosslope = cos( Fmin(GTConst::max_slope, slope[RC[i][1]][RC[i][2]]) * GTConst::Pi/180. );
 		//	V->co[i] = interpolate_soil2(lmin, n->co[l]*cosslope, Nl, dz, T, i);
 			V[i] = interpolate_soil2(lmin, n[l]*cosslope, Nl, dz, T, i);
 		}

 	//	temp2 = join_strings(filename, temp1);
 		temp2 = filename + temp1;
 	//	write_map_vector(temp2, type, format, V, UV, number_novalue, J, slope->nrh, slope->nch);
 		write_map_vector(temp2, type, format, V, UV, number_novalue, J, slope.getRows()-1, slope.getCols()-1);

 	//	free_doublevector(V);
 	//	free(temp1);
 	//	free(temp2);
 	}
 }


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

//void fill_output_vectors(double Dt, double W, Energy *egy, SNOW *snow, GLACIER *glac, WATER *wat, METEO *met, PAR *par, TIMES *time, TOPO *top){
  void fill_output_vectors(double Dt, double W, Energy *egy, Snow *snow, Glacier *glac, Water *wat, Meteo *met, Par *par, Times *time, Topo *top){

	long i, j;
	
	for (j=1; j<=par->total_pixel; j++) {

	//	if(par->output_snow->co[i_sim]>0){
		if(par->output_snow[i_sim]>0){
		//	if(strcmp(files[fsnowmelt] , string_novalue) != 0) snow->MELTED->co[j] += snow->melted->co[j];
			if(strcmp(files[fsnowmelt] , string_novalue) != 0) snow->MELTED[j] += snow->melted[j];
		//	if(strcmp(files[fsnowsubl] , string_novalue) != 0) snow->SUBL->co[j] += snow->subl->co[j];
			if(strcmp(files[fsnowsubl] , string_novalue) != 0) snow->SUBL[j] += snow->subl[j];
			if(strcmp(files[fsndur] , string_novalue) != 0){ 
			//	if(snow->yes->co[j] == 1) snow->t_snow->co[j] += Dt/GTConst::secinday;
				if(snow->yes[j] == 1) snow->t_snow[j] += Dt/GTConst::secinday;
			}
		}
		
	//	if(par->max_glac_layers>0 && par->output_glac->co[i_sim]>0){
		if(par->max_glac_layers>0 && par->output_glac[i_sim]>0){
		//	if(strcmp(files[fglacmelt] , string_novalue) != 0) glac->MELTED->co[j] += glac->melted->co[j];
			if(strcmp(files[fglacmelt] , string_novalue) != 0) glac->MELTED[j] += glac->melted[j];
		//	if(strcmp(files[fglacsubl] , string_novalue) != 0) glac->SUBL->co[j] += glac->subl->co[j];
			if(strcmp(files[fglacsubl] , string_novalue) != 0) glac->SUBL[j] += glac->subl[j];
		}
		
	//	if(par->output_surfenergy->co[i_sim]>0){
		if(par->output_surfenergy[i_sim]>0){
		//	if(strcmp(files[fradnet] , string_novalue) != 0) egy->Rn_mean->co[j] += (egy->SW->co[j]+egy->LW->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			//cout << "files[fradnet]=" << files[fradnet] << " par->output_surfenergy[i_sim]=" << par->output_surfenergy[i_sim] << endl;cout << " egy->Rn_mean[j]=" << egy->Rn_mean[j] << " egy->SW[j]=" << egy->SW[j] << " egy->LW[j]=" << egy->LW[j] << endl;
			if(strcmp(files[fradnet] , string_novalue) != 0) egy->Rn_mean[j] += (egy->SW[j]+egy->LW[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fradLWin] , string_novalue) != 0) egy->LWin_mean->co[j] += (egy->LWin->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradLWin] , string_novalue) != 0) egy->LWin_mean[j] += (egy->LWin[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fradLW] , string_novalue) != 0) egy->LW_mean->co[j] += (egy->LW->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradLW] , string_novalue) != 0) egy->LW_mean[j] += (egy->LW[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fradSW] , string_novalue) != 0)  egy->SW_mean->co[j] += (egy->SW->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradSW] , string_novalue) != 0)  egy->SW_mean[j] += (egy->SW[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fradSWin] , string_novalue) != 0) egy->Rswdown_mean->co[j] += (egy->SWin->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradSWin] , string_novalue) != 0) egy->Rswdown_mean[j] += (egy->SWin[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fradSWinbeam] , string_novalue) != 0) egy->Rswbeam_mean->co[j] += (egy->SWinb->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fradSWinbeam] , string_novalue) != 0) egy->Rswbeam_mean[j] += (egy->SWinb[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fG] , string_novalue) != 0) egy->SEB_mean->co[j] += (egy->G->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fG] , string_novalue) != 0) egy->SEB_mean[j] += (egy->G[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fH] , string_novalue) != 0) egy->H_mean->co[j] += (egy->H->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fH] , string_novalue) != 0) egy->H_mean[j] += (egy->H[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fLE] , string_novalue) != 0) egy->ET_mean->co[j] += (egy->LE->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fLE] , string_novalue) != 0) egy->ET_mean[j] += (egy->LE[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
		//	if(strcmp(files[fTs] , string_novalue) != 0) egy->Ts_mean->co[j] += (egy->Ts->co[j])*Dt/(par->output_surfenergy->co[i_sim]*3600.);
			if(strcmp(files[fTs] , string_novalue) != 0) egy->Ts_mean[j] += (egy->Ts[j])*Dt/(par->output_surfenergy[i_sim]*3600.);
			
			if(strcmp(files[fshadow] , string_novalue) != 0){
			//	if (egy->shad->co[j] >= 0) egy->nDt_sun->co[j] ++;
				if (egy->shad[j] >= 0) egy->nDt_sun[j] ++;
			//	if (egy->shad->co[j] == 0) egy->nDt_shadow->co[j] ++;
				if (egy->shad[j] == 0) egy->nDt_shadow[j] ++;
			}
		}
	//	if(par->output_meteo->co[i_sim]>0){
		if(par->output_meteo[i_sim]>0){
			if(strcmp(files[fprec] , string_novalue) != 0){
			//	wat->PrTOT_mean->co[j] = wat->Pt->co[j];
				wat->PrTOT_mean[j] = wat->Pt[j];
			//	wat->PrSNW_mean->co[j] = wat->Ps->co[j];
				wat->PrSNW_mean[j] = wat->Ps[j];
			}
		}	
		
	//	if(time->JD_plots->nh > 1 && W>0){
		if(time->JD_plots.size() > 1 && W>0){
		//	if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->Hgplot->co[j] += egy->Hgp->co[j];
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->Hgplot[j] += egy->Hgp[j];
		//	if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0) egy->Hvplot->co[j] += egy->Hvp->co[j];
			if(strcmp(files[pH] , string_novalue) != 0 || strcmp(files[pHv] , string_novalue) != 0) egy->Hvplot[j] += egy->Hvp[j];
		//	if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LEgplot->co[j] += egy->LEgp->co[j];
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LEgplot[j] += egy->LEgp[j];
		//	if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0) egy->LEvplot->co[j] += egy->LEvp->co[j];
			if(strcmp(files[pLE] , string_novalue) != 0 || strcmp(files[pLEv] , string_novalue) != 0) egy->LEvplot[j] += egy->LEvp[j];
		//	if(strcmp(files[pSWin] , string_novalue) != 0) egy->SWinplot->co[j] += egy->SWinp->co[j];
			if(strcmp(files[pSWin] , string_novalue) != 0) egy->SWinplot[j] += egy->SWinp[j];
		//	if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->SWgplot->co[j] += egy->SWgp->co[j];
			if(strcmp(files[pSWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->SWgplot[j] += egy->SWgp[j];
		//	if(strcmp(files[pSWv] , string_novalue) != 0) egy->SWvplot->co[j] += egy->SWvp->co[j];
			if(strcmp(files[pSWv] , string_novalue) != 0) egy->SWvplot[j] += egy->SWvp[j];
		//	if(strcmp(files[pLWin] , string_novalue) != 0) egy->LWinplot->co[j] += egy->LWinp->co[j];
			if(strcmp(files[pLWin] , string_novalue) != 0) egy->LWinplot[j] += egy->LWinp[j];
		//	if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LWgplot->co[j] += egy->LWgp->co[j];
			if(strcmp(files[pLWg] , string_novalue) != 0 || strcmp(files[pG] , string_novalue) != 0) egy->LWgplot[j] += egy->LWgp[j];
		//	if(strcmp(files[pLWv] , string_novalue) != 0) egy->LWvplot->co[j] += egy->LWvp->co[j];
			if(strcmp(files[pLWv] , string_novalue) != 0) egy->LWvplot[j] += egy->LWvp[j];
		//	if(strcmp(files[pTs] , string_novalue) != 0) egy->Tsplot->co[j] += egy->Tsp->co[j];
			if(strcmp(files[pTs] , string_novalue) != 0) egy->Tsplot[j] += egy->Tsp[j];
		//	if(strcmp(files[pTg] , string_novalue) != 0) egy->Tgplot->co[j] += egy->Tgp->co[j];
			if(strcmp(files[pTg] , string_novalue) != 0) egy->Tgplot[j] += egy->Tgp[j];
		//	if(strcmp(files[pD] , string_novalue) != 0) snow->Dplot->co[j] += W * DEPTH(top->rc_cont->co[j][1], top->rc_cont->co[j][2], snow->S->lnum, snow->S->Dzl);
			if(strcmp(files[pD] , string_novalue) != 0) snow->Dplot[j] += W * DEPTH(top->rc_cont[j][1], top->rc_cont[j][2], snow->S->lnum, snow->S->Dzl);
		//	if(strcmp(files[pTa] , string_novalue) != 0) met->Taplot->co[j] += W * met->Tgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]];
			if(strcmp(files[pTa] , string_novalue) != 0) met->Taplot[j] += W * met->Tgrid[top->rc_cont[j][1]][top->rc_cont[j][2]];
		//	if(strcmp(files[pRH] , string_novalue) != 0) met->RHplot->co[j] += W * met->RHgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]];
			if(strcmp(files[pRH] , string_novalue) != 0) met->RHplot[j] += W * met->RHgrid[top->rc_cont[j][1]][top->rc_cont[j][2]];
			if(strcmp(files[pVspd] , string_novalue) != 0 || strcmp(files[pVdir] , string_novalue) != 0){
			//	met->Vxplot->co[j] -= W * met->Vgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]] * sin(met->Vdir->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]]*GTConst::Pi/180.);
				met->Vxplot[j] -= W * met->Vgrid[top->rc_cont[j][1]][top->rc_cont[j][2]] * sin(met->Vdir[top->rc_cont[j][1]][top->rc_cont[j][2]]*GTConst::Pi/180.);
			//	met->Vyplot->co[j] -= W * met->Vgrid->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]] * cos(met->Vdir->co[top->rc_cont->co[j][1]][top->rc_cont->co[j][2]]*GTConst::Pi/180.);
				met->Vyplot[j] -= W * met->Vgrid[top->rc_cont[j][1]][top->rc_cont[j][2]] * cos(met->Vdir[top->rc_cont[j][1]][top->rc_cont[j][2]]*GTConst::Pi/180.);
			}
		}
		if (par->state_pixel==1){
		//	if (par->jplot->co[j] > 0 && par->Dtplot_point->co[i_sim]>0){
			if (par->jplot[j] > 0 && par->Dtplot_point[i_sim]>0){
				for(i=0; i<otot; i++){
				//	odpnt[i][par->jplot->co[j]-1] += odp[i][par->jplot->co[j]-1];
					odpnt[i][par->jplot[j]-1] += odp[i][par->jplot[j]-1];
				}
			}
		}
	}
	
	if (par->state_basin==1){
	//	if (par->Dtplot_basin->co[i_sim]>0){
		if (par->Dtplot_basin[i_sim]>0){
			for(i=0; i<ootot; i++){
				odbsn[i] += odb[i];
			}
		}
	}
}
