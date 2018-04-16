
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

#include "constants.h"
#include "struct.geotop.h"
#include "snow.h"
#include "PBSM.h"
#include "blowingsnow.h"
#include "meteo.h"
#include "vegetation.h"
#include "energy.balance.h"
#include "meteodata.h"

extern long number_novalue, number_absent;

extern T_INIT *UV;
extern char *logfile;
extern long Nl, Nr, Nc;
extern char *WORKING_DIRECTORY;
extern double **odp;

extern long i_sim;

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

void windtrans_snow(SNOW *snow, METEO *met, WATER *wat, LAND *land, TOPO *top, PAR *par, double t0){
	
	long r, c, j, l;
	static long line_interp;
	short lu, lux, yes, sux;
	double t, Dt, Dt0;
	double zmeas, D, wice, fsnow, fc;
	double canopy_height_over_snow, rho_snow_surface;
	double PBSM_fetch = 1000.;
	double dErdt;
	long ns;
	FILE *f;
	
	//double U=0;
	
	f=fopen(logfile,"a");	
		
	if(t0 == 0) line_interp = 0;
	
	t = 0.;	
	
	//t0 = time from beginning simulation
	//t = time from beginning of time step
		
	do{
		
		//find Dt
		Dt = par->Dt_PBSM;
		if( t+Dt > par->Dt ) Dt = par->Dt - t;
		Dt0 = Dt;
				
		//meteo distribution
		meteo_distr(met->line_interp_Bsnow, met->line_interp_Bsnow_LR, met, wat, top, par, par->init_date->co[i_sim]+t0/secinday, 
					par->init_date->co[i_sim]+(t0+t)/secinday, par->init_date->co[i_sim]+(t0+t+Dt)/secinday);
		

		//vegetation
		for(lux=1; lux<=par->n_landuses; lux++) {
			if(par->vegflag->co[lux]==1){
				time_interp_linear(par->init_date->co[i_sim]+t0/secinday, par->init_date->co[i_sim]+(t0+t)/secinday, par->init_date->co[i_sim]+(t0+t+Dt)/secinday,
								   land->vegparv[lux-1], land->vegpars[lux-1], land->NumlinesVegTimeDepData[lux-1], jdvegprop+1, 0, 0, &line_interp);
			}
		}
		
		//loop for every pixel
		for (r=1; r<=Nr; r++) {
			for (c=1; c<=Nc; c++) {
				if ( (long)land->LC->co[r][c]!=number_novalue) {
						

					D = DEPTH(r, c, snow->S->lnum, snow->S->Dzl);
					wice = DEPTH(r, c, snow->S->lnum, snow->S->w_ice);
					
					//U += met->Vgrid->co[r][c]/par->total_pixel;

					if(wice > par->Wmin_BS){
																														
						//find canopy_height_over_snow
						lu = (short)land->LC->co[r][c];
						canopy_height_over_snow = 0.0;						
						//zmeas = Fmax(0.1, met->st->Vheight->co[1]-1.E-3*D);
						zmeas = met->st->Vheight->co[1];	
						
						for(j=1;j<=jdvegprop;j++){
							if( (long)land->vegparv[lu-1][j-1+1] != number_novalue ){
								land->vegpar->co[j] = land->vegparv[lu-1][j-1+1];
							}else{
								land->vegpar->co[j] = land->ty->co[lu][j+jHveg-1];
							}
						}
						
						if(D > land->vegpar->co[jdz0thresveg]){
							fsnow=1.0;
						}else if(D > land->vegpar->co[jdz0thresveg2]){
							fsnow=(D-land->vegpar->co[jdz0thresveg2])/(land->vegpar->co[jdz0thresveg]-land->vegpar->co[jdz0thresveg2]);
						}else{
							fsnow=0.0;
						}						
						fc = land->vegpar->co[jdcf] * pow(1.0-fsnow,land->vegpar->co[jdexpveg]);//canopy fraction
						if(land->vegpar->co[jdLSAI]<LSAIthres) fc=0.0; 
						if(fc>0) canopy_height_over_snow += fc*Fmax(land->vegpar->co[jdHveg]-D, 0.)*1.E-3;
											
						//rearrange snow layers						
						ns = snow->S->lnum->co[r][c];
						sux = copy_statevar_from3D_to1D(r, c, snow->S, snow->S_for_BS);
												
						if ( snow->S_for_BS->w_ice->co[ns] < par->Wice_PBSM ) {
							l=ns;
							do{
								l--;
								sux = set_snowice_min(par->alpha_snow, r, c, snow->S_for_BS, ns, l, par->Wice_PBSM);
							}while ( l > 1 && snow->S_for_BS->w_ice->co[ns] < par->Wice_PBSM );
						}
						
						//find equilibrium snow trasport
						rho_snow_surface = (snow->S_for_BS->w_ice->co[ns] + snow->S_for_BS->w_liq->co[ns]) / (1.E-3*snow->S_for_BS->Dzl->co[ns]);
												
						Pbsm (r, c, PBSM_fetch, land->ty->co[lu][jN], 1.E-3*land->ty->co[lu][jdv], canopy_height_over_snow, rho_snow_surface, zmeas,
							  met->Vgrid->co[r][c], met->Tgrid->co[r][c], met->RHgrid->co[r][c], &(snow->Qtrans->co[r][c]), &(snow->Qsub->co[r][c]),
							  &(snow->Qsalt->co[r][c]), D, top->slope->co[r][c], f);
				
	
					}else{
						
						snow->Qtrans->co[r][c] = 0.0;
						snow->Qsub->co[r][c] = 0.0;
						
					}
											
				}
			}
		}
		
		if (par->point_sim == 0) {
			
			set_inhomogeneous_fetch(snow, met, land, par, top, &yes);
			if(yes==1){
				for (r=1; r<=Nr; r++) {
					for (c=1; c<=Nc; c++) {
						if ( (long)land->LC->co[r][c]!=number_novalue) {
							wice = DEPTH(r, c, snow->S->lnum, snow->S->w_ice);
							dErdt = sqrt(pow(snow->Qsub_x->co[r][c], 2.)+pow(snow->Qsub_y->co[r][c], 2.)) - snow->Nabla2_Qtrans->co[r][c];
							
							if( snow->S->lnum->co[r][c] > 0 ){
								if(dErdt*Dt > 0.99*wice) Dt=0.99*wice/dErdt;
							}
							
						}
					}
				}
				
				Dt=Fmin(Dt,Dt0);
				set_windtrans_snow(Dt, t0+t, snow, met, land, par, f);
				print_windtrans_snow(Dt, snow, par, top, met, land->LC);
			}
						
		}else if (par->snow_wind_compaction_1D == 1){
			
			for (r=1; r<=Nr; r++) {
				for (c=1; c<=Nc; c++) {
					if ( (long)land->LC->co[r][c]!=number_novalue) {
						wind_packing(snow, par, r, c, Dt);
					}
				}
			}
		
		}

		t += Dt;	
		
	}while (t<par->Dt);
	
	
	fclose(f);
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

void set_inhomogeneous_fetch(SNOW *snow, METEO *met, LAND *land, PAR *par, TOPO *top, short *yes){

	long i, r, c, num_change, r0, c0;
	double dx, dy, F1, F2;
	double Qtrans=0.0;
	double Qup, Qdown, Sup, Sdown, k_snowred;
	
	dx=UV->U->co[1];                                     
	dy=UV->U->co[2];      
	
	F1=par->fetch_up/3.;
	F2=par->fetch_down/3.;
	
	*yes = 0;
	
	//initialize Nabla2_Qtrans(snow deposition if >0, snow erosion if <0)
	initialize_doublematrix(snow->Nabla2_Qtrans, 0.0);
		 			
	//check if there is snow transport
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if( (long)land->LC->co[r][c]!=number_novalue){				
				Qtrans += snow->Qtrans->co[r][c]/(double)par->total_pixel;
			}
		}
	}
	
				
	//if there is blowing snow 
	if( Qtrans>0 ){
		
		*yes = 1;
		
		//wind in direction west-east
		/*set_no_value(snow->Qtrans, land->LC, number_novalue);
		extend_topography_row(snow->Qtrans, number_novalue);
		
		set_no_value(snow->Qsub, land->LC, number_novalue);
		extend_topography_row(snow->Qsub, number_novalue);*/

		for(r=1;r<=Nr;r++){
			
			//find west-east component
			for(c=1;c<=Nc;c++){
				snow->Qtrans_x->co[r][c]=fabs(snow->Qtrans->co[r][c]*(-sin(met->Vdir->co[r][c]*Pi/180.)));
				snow->Qsub_x->co[r][c]=fabs(snow->Qsub->co[r][c]*(-sin(met->Vdir->co[r][c]*Pi/180.)));
			}
			
			//find when there is a wind direction inversion
			initialize_longvector(snow->change_dir_wind,0);
			num_change=0;
			c=1;
			
			num_change++;				
			c0=c;
			snow->change_dir_wind->co[num_change]=c;
									
			do{
				c=c0;
				do{
					c++;
				}while( (-sin(met->Vdir->co[r][c]*Pi/180.))*(-sin(met->Vdir->co[r][c0]*Pi/180.))>0 && c<Nc );
				
				num_change++;				
				c0=c;				
				snow->change_dir_wind->co[num_change]=c;
												
			}while(c0<Nc);
			
			for(i=1;i<num_change;i++){
				//east wind
				if( (-sin(met->Vdir->co[r][snow->change_dir_wind->co[i]]*Pi/180.)) > 0 ){	
					if(par->upwindblowingsnow==1 && snow->change_dir_wind->co[i]!=1) snow->Qtrans_x->co[r][snow->change_dir_wind->co[i]]=0.0;
					snow->Qtrans_x->co[r][snow->change_dir_wind->co[i]]=0.0;
					for(c=snow->change_dir_wind->co[i]+1;c<=snow->change_dir_wind->co[i+1];c++){
						if(snow->change_dir_wind->co[i+1]==Nc || (snow->change_dir_wind->co[i+1]!=Nc && c<snow->change_dir_wind->co[i+1])){						
							Qup = snow->Qtrans_x->co[r][c-1];
							Qdown = snow->Qtrans_x->co[r][c];
							Sup = snow->Qsub_x->co[r][c-1];
							Sdown = snow->Qsub_x->co[r][c];
							if(Qdown>=Qup){
								snow->Qtrans_x->co[r][c] = (Qdown + Qup*F1/dx)/(1.+F1/dx);
								snow->Qsub_x->co[r][c] = (Sdown + Sup*F1/dx)/(1.+F1/dx);
							}else{
								snow->Qtrans_x->co[r][c] = (Qdown + Qup*F2/dx)/(1.+F2/dx);
								snow->Qsub_x->co[r][c] = (Sdown + Sup*F2/dx)/(1.+F2/dx);
							}	
							Qdown = snow->Qtrans_x->co[r][c];
							snow->Nabla2_Qtrans->co[r][c] += ( Qup - Qdown )/dx;	
						}
					}
					
				//west wind	
				}else{	
					if(par->upwindblowingsnow==1 && snow->change_dir_wind->co[i+1]!=Nc) snow->Qtrans_x->co[r][snow->change_dir_wind->co[i+1]-1]=0.0;	
					snow->Qtrans_x->co[r][snow->change_dir_wind->co[i+1]-1]=0.0;	
					for(c=snow->change_dir_wind->co[i+1]-1;c>=snow->change_dir_wind->co[i];c--){
						if(snow->change_dir_wind->co[i+1]==Nc || (snow->change_dir_wind->co[i+1]!=Nc && c<snow->change_dir_wind->co[i+1]-1)){
							Qup = snow->Qtrans_x->co[r][c+1];
							Qdown = snow->Qtrans_x->co[r][c];
							Sup = snow->Qsub_x->co[r][c+1];
							Sdown = snow->Qsub_x->co[r][c];
							if(Qdown>=Qup){
								snow->Qtrans_x->co[r][c] = (Qdown + Qup*F1/dx)/(1.+F1/dx);
								snow->Qsub_x->co[r][c] = (Sdown + Sup*F1/dx)/(1.+F1/dx);
							}else{
								snow->Qtrans_x->co[r][c] = (Qdown + Qup*F2/dx)/(1.+F2/dx);
								snow->Qsub_x->co[r][c] = (Sdown + Sup*F2/dx)/(1.+F2/dx);
							}	
							Qdown = snow->Qtrans_x->co[r][c];
							snow->Nabla2_Qtrans->co[r][c] += ( Qup - Qdown )/dx;	
						}
					}
				}
			}
		}	
		
		
		//wind in direction south-north
		/*set_no_value(snow->Qtrans, land->LC, number_novalue);	
		extend_topography_column(snow->Qtrans, number_novalue);

		set_no_value(snow->Qsub, land->LC, number_novalue);
		extend_topography_column(snow->Qsub, number_novalue);*/
		
		for(c=1;c<=Nc;c++){
			
			//find south-north component
			for(r=1;r<=Nr;r++){
				snow->Qtrans_y->co[r][c]=fabs(snow->Qtrans->co[r][c]*(-cos(met->Vdir->co[r][c]*Pi/180.)));
				snow->Qsub_y->co[r][c]=fabs(snow->Qsub->co[r][c]*(-cos(met->Vdir->co[r][c]*Pi/180.)));
			}
			
			//find when there is a wind direction inversion
			initialize_longvector(snow->change_dir_wind,0);
			num_change=0;
			r=1;
			
			num_change++;				
			r0=r;
			snow->change_dir_wind->co[num_change]=r;
						
			do{
				r=r0;
				do{
					r++;
				}while( (-cos(met->Vdir->co[r][c]*Pi/180.))*(-cos(met->Vdir->co[r0][c]*Pi/180.))>0 && r<Nr );
				
				num_change++;				
				r0=r;				
				snow->change_dir_wind->co[num_change]=r;
								
			}while(r0<Nr);
			
			for(i=1;i<num_change;i++){
				//south wind
				if( (-cos(met->Vdir->co[snow->change_dir_wind->co[i]][c]*Pi/180.)) < 0 ){
					if(par->upwindblowingsnow==1 && snow->change_dir_wind->co[i]!=1) snow->Qtrans_y->co[snow->change_dir_wind->co[i]][c]=0.0;					
					snow->Qtrans_y->co[snow->change_dir_wind->co[i]][c]=0.0;					
					for(r=snow->change_dir_wind->co[i]+1;r<=snow->change_dir_wind->co[i+1];r++){
						if(snow->change_dir_wind->co[i+1]==Nr || (snow->change_dir_wind->co[i+1]!=Nr && r<snow->change_dir_wind->co[i+1])){
							Qup = snow->Qtrans_y->co[r-1][c];
							Qdown = snow->Qtrans_y->co[r][c];
							Sup = snow->Qsub_y->co[r-1][c];
							Sdown = snow->Qsub_y->co[r][c];
							if(Qdown>=Qup){
								snow->Qtrans_y->co[r][c] = (Qdown + Qup*F1/dy)/(1.+F1/dy);
								snow->Qsub_y->co[r][c] = (Sdown + Sup*F1/dy)/(1.+F1/dy);
							}else{
								snow->Qtrans_y->co[r][c] = (Qdown + Qup*F2/dy)/(1.+F2/dy);
								snow->Qsub_y->co[r][c] = (Sdown + Sup*F2/dy)/(1.+F2/dy);
							}	
							Qdown = snow->Qtrans_y->co[r][c];
							snow->Nabla2_Qtrans->co[r][c] += ( Qup - Qdown )/dy;	
						}
					}
					
				//north wind	
				}else{
					if(par->upwindblowingsnow==1 && snow->change_dir_wind->co[i+1]!=Nr) snow->Qtrans_y->co[snow->change_dir_wind->co[i+1]-1][c]=0.0;		
					snow->Qtrans_y->co[snow->change_dir_wind->co[i+1]-1][c]=0.0;		
					for(r=snow->change_dir_wind->co[i+1]-1;r>=snow->change_dir_wind->co[i];r--){
						if(snow->change_dir_wind->co[i+1]==Nr || (snow->change_dir_wind->co[i+1]!=Nr && r<snow->change_dir_wind->co[i+1]-1)){
							Qup = snow->Qtrans_y->co[r+1][c];
							Qdown = snow->Qtrans_y->co[r][c];
							Sup = snow->Qsub_y->co[r+1][c];
							Sdown = snow->Qsub_y->co[r][c];
							if(Qdown>=Qup){
								snow->Qtrans_y->co[r][c] = (Qdown + Qup*F1/dy)/(1.+F1/dy);
								snow->Qsub_y->co[r][c] = (Sdown + Sup*F1/dy)/(1.+F1/dy);
							}else{
								snow->Qtrans_y->co[r][c] = (Qdown + Qup*F2/dy)/(1.+F2/dy);
								snow->Qsub_y->co[r][c] = (Sdown + Sup*F2/dy)/(1.+F2/dy);
							}	
							Qdown = snow->Qtrans_y->co[r][c];
							snow->Nabla2_Qtrans->co[r][c] += ( Qup - Qdown )/dy;	
						}
					}
				}
			}
		}
		
		//Adjusting snow init depth in case of steep slope (contribution by Stephan Gruber)
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if( (long)land->LC->co[r][c]!=number_novalue ){	
					if (par->snow_curv > 0 && top->slope->co[r][c] > par->snow_smin){
						if (top->slope->co[r][c] <= par->snow_smax){
							k_snowred = ( exp(-pow(top->slope->co[r][c] - par->snow_smin, 2.)/par->snow_curv) -
										 exp(-pow(par->snow_smax, 2.)/par->snow_curv) );
						}else{
							k_snowred = 0.0;
						}
						if( snow->Nabla2_Qtrans->co[r][c] > 0 ) snow->Nabla2_Qtrans->co[r][c] *= k_snowred;
					}
				}
			}
		}		
	}
}
				
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

void set_windtrans_snow(double Dt, double t, SNOW *snow, METEO *met, LAND *land, PAR *par, FILE *f){
	
	long i, r, c, ns;
	double Qsub, DW, DWl, Wtrans_tot=0.0, Wsubl_tot=0.0, Utot=0.0;
	short ok=1;

	//density of wind transported snow
	float rho_wind_transported_snow = 400.0;
			
	//update snow depth
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if((long)land->LC->co[r][c]!=number_novalue){
					
				ns = snow->S->lnum->co[r][c];
				
				Qsub = sqrt(pow(snow->Qsub_x->co[r][c], 2.)+pow(snow->Qsub_y->co[r][c], 2.));
				
				Wsubl_tot += Dt*Qsub/(double)par->total_pixel;
				Wtrans_tot += Dt*snow->Nabla2_Qtrans->co[r][c]/(double)par->total_pixel;
				Utot += met->Vgrid->co[r][c]/(double)par->total_pixel;
					
				DW = Dt*(snow->Nabla2_Qtrans->co[r][c] - Qsub);	

				if(ns>0){	//snow on the soil
					
					if(DW<0){	//snow eroded
							
						i=ns;
						DWl=0.0;
						
						do{					
							if(i<ns){ 
								ok=0;
								if(snow->S->w_ice->co[i+1][r][c]<0){
									DW=snow->S->w_ice->co[i+1][r][c];
									DWl=snow->S->w_liq->co[i+1][r][c];
									snow->S->w_ice->co[i+1][r][c]=0.0;
									snow->S->w_liq->co[i+1][r][c]=0.0;
									snow->S->Dzl->co[i+1][r][c]=0.0;
									snow->S->lnum->co[r][c]-=1;
								}
							}
							
							snow->S->Dzl->co[i][r][c]*=(snow->S->w_ice->co[i][r][c]+DW)/snow->S->w_ice->co[i][r][c];
							snow->S->w_ice->co[i][r][c]+=DW;				//kg/m2
							snow->S->w_liq->co[i][r][c]+=DWl;				//kg/m2
															
							i--;
								
						}while(snow->S->w_ice->co[i+1][r][c]<0 && i>0);
							
						if(i==0 && snow->S->w_ice->co[i+1][r][c]<0){
							snow->S->w_ice->co[i+1][r][c]=0.0;				//kg/m2
							snow->S->Dzl->co[i+1][r][c]=0.0;	//mm	
							snow->S->lnum->co[r][c]=0;
						}
							
					}else{	//snow drifted
							
						i = snow->S->lnum->co[r][c];
						snow->S->w_ice->co[i][r][c]+=DW;
						snow->S->Dzl->co[snow->S->lnum->co[r][c]][r][c] += 1.0E+3*DW/rho_wind_transported_snow;
							
					}
						
				}else{	//snot not on the soil
						
					if(DW>0){
							
						snow->S->w_ice->co[1][r][c]+=DW;
						snow->S->Dzl->co[1][r][c]+=1.0E+3*DW/rho_wind_transported_snow;
						snow->S->T->co[1][r][c]=Fmin(-1.,met->Tgrid->co[r][c]);
						
					}
					
				}

				wind_packing(snow, par, r, c, Dt);

				snow_layer_combination(par->alpha_snow, r, c, snow->S, met->Tgrid->co[r][c], par->inf_snow_layers, par->max_weq_snow, 1.E10, f);

			}
		}
	}

	fprintf(f,"BLOWING SNOW: t:%f Wsubl:%e Wtrans:%e U:%f Dt:%f \n",t,Wsubl_tot,Wtrans_tot,Utot, Dt);
	
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

void print_windtrans_snow(double Dt, SNOW *snow, PAR *par, TOPO *top, METEO *met, DOUBLEMATRIX *LC){

	long i, r, c;
	double Qsub;
	
	for(i=1; i<=par->total_pixel; i++){
		
		r = top->rc_cont->co[i][1];
		c = top->rc_cont->co[i][2];
		
		Qsub = sqrt(pow(snow->Qsub_x->co[r][c], 2.)+pow(snow->Qsub_y->co[r][c], 2.));
								
		if(par->output_snow->co[i_sim]>0){
			snow->Wtrans_plot->co[r][c] += Dt*snow->Nabla2_Qtrans->co[r][c];			
			snow->Wsubl_plot->co[r][c] += Dt*Qsub;
		}
		
		if(par->Dtplot_point->co[i_sim] > 1.E-5 && par->state_pixel == 1 && par->jplot->co[i] > 0){
			odp[oblowingsnowtrans][par->jplot->co[i]-1] -= Dt*snow->Nabla2_Qtrans->co[r][c];
			odp[oblowingsnowsubl][par->jplot->co[i]-1] += Dt*Qsub;
		}
		
	}
}


//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
void extend_topography(DOUBLEMATRIX *M, double novalue){
	
	long r,c,rr,cc;
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(M->co[r][c]==novalue){
				find_the_nearest(r, c, novalue, M, &rr, &cc);
				M->co[r][c]=M->co[rr][cc];
			}
		}
	}	
}

void extend_topography_row(DOUBLEMATRIX *M, double novalue){
	
	long r,c,rr,cc;
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(M->co[r][c]==novalue){
				find_the_nearest_row(r, c, novalue, M, &rr, &cc);
				M->co[r][c]=M->co[rr][cc];
			}
		}
	}	
}

void extend_topography_column(DOUBLEMATRIX *M, double novalue){
	
	long r,c,rr,cc;
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(M->co[r][c]==novalue){
				find_the_nearest_column(r, c, novalue, M, &rr, &cc);
				M->co[r][c]=M->co[rr][cc];
			}
		}
	}	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/				

void find_the_nearest(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc){
	
	long i=0;
	short k;
		
	do{
		i++;
		k=0;
		if(k==0) k=no_novalue(r-i, c,   M, novalue, rr, cc);
		if(k==0) k=no_novalue(r-i, c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r  , c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r+i, c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r+i, c,   M, novalue, rr, cc);
		if(k==0) k=no_novalue(r+i, c-i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r,   c-i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r-i, c-i, M, novalue, rr, cc);
	}while(k==0);

}

void find_the_nearest_row(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc){
	
	long i=0;
	short k;
	
	*rr = r;
	*cc = c;
	
	do{
		i++;
		k=0;
		//if(k==0) k=no_novalue(r-i, c,   M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r-i, c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r  , c+i, M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r+i, c+i, M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r+i, c,   M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r+i, c-i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r,   c-i, M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r-i, c-i, M, novalue, rr, cc);
	}while(k==0 && i<M->nch);
	
	
	i=0;
	if(k==0){
		do{
			i++;
			k=0;
			if(k==0) k=no_novalue(r-i, c,   M, novalue, rr, cc);
			if(k==0) k=no_novalue(r-i, c+i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r  , c+i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r+i, c+i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r+i, c,   M, novalue, rr, cc);
			if(k==0) k=no_novalue(r+i, c-i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r,   c-i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r-i, c-i, M, novalue, rr, cc);

		}while(k==0);
	}
	
	
}

void find_the_nearest_column(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc){
	
	long i=0;
	short k;
	
	do{
		i++;
		k=0;
		if(k==0) k=no_novalue(r-i, c,   M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r-i, c+i, M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r  , c+i, M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r+i, c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r+i, c,   M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r+i, c-i, M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r,   c-i, M, novalue, rr, cc);
		//if(k==0) k=no_novalue(r-i, c-i, M, novalue, rr, cc);
	}while(k==0 && i<M->nrh);
	
	i=0;
	if(k==0){
		do{
			i++;
			k=0;
			if(k==0) k=no_novalue(r-i, c,   M, novalue, rr, cc);
			if(k==0) k=no_novalue(r-i, c+i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r  , c+i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r+i, c+i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r+i, c,   M, novalue, rr, cc);
			if(k==0) k=no_novalue(r+i, c-i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r,   c-i, M, novalue, rr, cc);
			if(k==0) k=no_novalue(r-i, c-i, M, novalue, rr, cc);
		}while(k==0);
	}
	
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/				

short no_novalue(long r, long c, DOUBLEMATRIX *M, double novalue, long *rr, long *cc){
	
	short k=0;
		
	if((r>=1 && r<=M->nrh) && (c>=1 && c<=M->nch)){
		if(M->co[r][c]!=novalue){
			k=1;
			*rr=r;
			*cc=c;
		}
	}
	
	return(k);
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/				
void set_no_value(DOUBLEMATRIX *M, DOUBLEMATRIX *N, double undef){
	
	long r, c;
	
	for(r=1; r<=M->nrh; r++){
		for(c=1; c<=M->nch; c++){
			if(N->co[r][c]==undef) M->co[r][c]=undef;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/				

void wind_packing(SNOW *snow, PAR *par, long r, long c, double Dt){
	
	long l;
	double CR;
	double overburden=0.0;
	float A4 = 7.81E-3;//m kg^-1
	float D4 = 0.0884;//m2 N-1	
	
	//from Liston(2007) - wind packing factor
	if(snow->Qtrans->co[r][c]>1.E-10){
			
		for (l=snow->S->lnum->co[r][c]; l>=1; l--) {
			
			overburden += (snow->S->w_ice->co[l][r][c]+snow->S->w_liq->co[l][r][c])/2.;
			
			//compactation at the surface (10%/hour if U8=8m/s U8t=4m/s => Qsalt=3.555342e-03 kg/m/s)
			CR = -A4 * snow->Qsalt->co[r][c];
			//decrease due to oberburden
			CR *= exp( -D4*g*overburden );		
			
			snow->S->Dzl->co[l][r][c] *= exp(CR*Dt);
			
			if (snow->S->w_ice->co[l][r][c]/(rho_w*snow->S->Dzl->co[l][r][c]*1.E-3) > par->snow_maxpor) {
				snow->S->Dzl->co[l][r][c] = 1.E3*snow->S->w_ice->co[l][r][c]/(rho_w*par->snow_maxpor);
			}
			
			overburden += (snow->S->w_ice->co[l][r][c]+snow->S->w_liq->co[l][r][c])/2.;
		}
	}
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/				
