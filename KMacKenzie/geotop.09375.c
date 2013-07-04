
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Matteo Dall'Amico, Riccardo Rigon, Emanuele Cordano

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


    /*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/
#include <sys/stat.h>
#include <fenv.h>
#include "struct.geotop.09375.h"
#include "input.09375.h"
#include "output.09375.h"
#include "times.h"
#include "constant.h"
#include "meteo.09375.h"
#include "energy.balance.h"
#include "water.balance_1D.h"
#include "water.balance_3D.h"
#include "pedo.funct.h"
#include "recovery.h"

void time_loop(ALLDATA *all);
void time_loop_superfast(ALLDATA *all);
void init_structs(TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet,
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times);

/*----------   1. Global variables  ------------*/
T_INIT *UV;
STRINGBIN *files;
char *error_file_name; /* name of the file error added by Emanuele Cordano */
long Nl; // total number of soil layers (constant in the whole basin)
long Nr; // total number of rows (of the map)
long Nc;// total number of columns (of the map)
double NoV;
long j,r,c;

/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/
int main(int argc,char *argv[]){
	//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); //for halting the process on arithmetic exceptions
	/*structs' declarations:*/
	ALLDATA *adt;

	/*dinamic allocations:*/
	UV=(T_INIT *)malloc(sizeof(T_INIT));
	if(!UV) t_error("UV was not allocated");

	adt=(ALLDATA *)malloc(sizeof(ALLDATA));
	if(!adt) t_error("adt was not allocated");

//---------------------------------------------------
//---------------------------------------------------

	adt->T=(TOPO *)malloc(sizeof(TOPO));
	if(!(adt->T)) t_error("top was not allocated");

	adt->S=(SOIL *)malloc(sizeof(SOIL));
	if(!(adt->S)) t_error("sl was not allocated");

	adt->L=(LAND *)malloc(sizeof(LAND));
	if(!(adt->L)) t_error("land was not allocated");

	adt->W=(WATER *)malloc(sizeof(WATER));
	if(!(adt->W)) t_error("water was not allocated");

	adt->P=(PAR *)malloc(sizeof(PAR));
	if(!(adt->P)) t_error("par was not allocated");

	adt->C=(CHANNEL *)malloc(sizeof(CHANNEL));
	if(!(adt->C)) t_error("channel was not allocated");

	adt->E=(ENERGY *)malloc(sizeof(ENERGY));
	if(!(adt->E)) t_error("egy was not allocated");

	adt->N=(SNOW *)malloc(sizeof(SNOW));
	if(!(adt->N)) t_error("snow was not allocated");

	adt->G=(GLACIER *)malloc(sizeof(GLACIER));
	if(!(adt->G)) t_error("glac was not allocated");

	adt->M=(METEO *)malloc(sizeof(METEO));
	if(!(adt->M)) t_error("met was not allocated");

	adt->I=(TIMES *)malloc(sizeof(TIMES));
	if(!(adt->I)) t_error("times was not allocated");


init_structs(adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I);

/*------------------    3.  Acquisition of input data and initialization    --------------------*/
get_all_input(argc, argv, adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I);

if(adt->P->superfast!=1) {
	write_init_condit(adt->M->st->Z->nh, adt->I, adt->W, adt->P, adt->T, adt->L, adt->S, adt->E, adt->N, adt->G);
	if(adt->P->en_balance!=1){// write the initial condition if energy.balance is switched off
		for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(adt->L->LC->co[r][c]!=NoV){//if the pixel is not a novalue
						for(j=1;j<=adt->P->rc->nrh;j++){
							if(r==adt->P->rc->co[j][1] && c==adt->P->rc->co[j][2])
								write_soil_output(0, j, adt->I->time, 0.0, adt->P->year0, adt->P->JD0, adt->P->rc, adt->S, PSImin, adt->P->Esoil);
						}
					}
				}
			}
		}
	}


/*----------------- Try to recover simulation if wished for -----------------*/
if (adt->P->recover == 1) {//recover from saved data
	double JD;
	int n=0;
	long day, month, year, hour, minute;
	char timestamp[15]; 
	double TH = adt->I->TH;
	
	date_time(adt->I->time, adt->P->year0, adt->P->JD0, 0.0, &JD, &day, &month, &year, &hour, &minute);
	//printf("Reading recovery files for date: %04ld%02ld%02ldT%02ld%02ld   simulation time:%f\n\n", 
	//	   y2,mo2,d2,h2,mi2, times->time+par->Dt);
	
	n = sprintf (timestamp, "%04ld%02ld%02ldT%02ld%02ld", year, month, day, hour, minute);
	if (n != 13){
		fprintf(stderr, "%s (%d): Not able to construct recovery timestamp", __FILE__, __LINE__);
		exit(1);
	} 

	recover_simulation(timestamp, adt->I, adt->W, adt->C, adt->S, adt->E, adt->N, adt->G);

	adt->I->TH    = TH;
	adt->I->time += (long)adt->P->Dt;
}


/*-----------------   4. Time-loop for the balances of water-mass and egy   -----------------*/
if(adt->P->superfast!=1) {// regular version
	time_loop(adt);
}else{// superfast
	time_loop_superfast(adt);
}

/*--------------------   5.Completion of the output files and deallocations  --------------------*/
dealloc_all(adt->T, adt->S, adt->L, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->M, adt->I);
free(adt);

printf("End of simulation!\n");
//stop_execution();
return 0;
}


/*----------------   6. The most important subroutine of the main: "time_loop"   ---------------*/
void time_loop(ALLDATA *all)
{


 /*begin of time-cycle:*/
 do{
	//printf("par->point_sim=%d, par->superfast=%ld",all->P->point_sim,all->P->superfast); stop_execution();
    updates_times(all->I, all->P);

	meteo_distr(all->M, all->E, all->W, all->T, all->N, all->I->time, all->P);

	//printf("\n ENERGY START %10.2f\n",all->I->time);
	//stop_execution();

    if(all->P->en_balance==1) {
    	energy_balance(all->I, all->P, all->L, all->T, all->S, all->M, all->W, all->E, all->N, all->G);
		}

	//printf("\n ENERGY END %10.2f\n",all->I->time);
	//stop_execution();

	//printf("\n MASS START %10.2f\n",all->I->time);
	//stop_execution();
	if(all->P->wat_balance==1){
		water_balance_1D(all->T, all->S, all->L, all->W, all->C, all->P, all->I->time);
	}else if(all->P->wat_balance==3){
		water_balance_3D(all);
	}

	//printf("\n MASS END %10.2f\n",all->I->time);
	//stop_execution();

	//printf("\n WRITE OUTPUT b %10.2f\n",all->I->time);
	//stop_execution();

	 write_output(all->I, all->W, all->C, all->P, all->T, all->L, all->S, all->E, all->N, all->G, all->M);

	//printf("\n WRITE OUTPUT e %10.2f\n",all->I->time);
	//stop_execution();

	//Increase TIME!
	all->I->time+=(long)all->P->Dt;

 }while(all->I->time<=(all->I->TH*3600.0));/*end of time-cycle*/

 free(all->I);

}

/*--------------------------------------------------------------------------------------------*/
void time_loop_superfast(ALLDATA *all)
{

	double t_station;
	long l,r,c,i;
	short sy=0;

 /*begin of time-cycle:*/
 do{

    updates_times(all->I, all->P);

	//meteo_distr(all->M, all->E, all->W, all->T, all->N, all->I->time, all->P);// was like this

    //INTERPOLATION OF METEO VARIABLES
	time_conversion(all->P->JD0, all->P->year0, all->I->time+all->P->Dt, all->M->st->JD0->co[1], all->M->st->Y0->co[1], &t_station);
	t_station+=(all->M->st->ST->co[1]-all->P->ST)*3600.0;
	meteo_interp(all->M->data[0], all->M->st->Dt->co[1], t_station, all->M->var[0]);
	//printf("\n ENERGY START %10.2f\n",all->I->time);
	//stop_execution();

	// STORE INITIAL CONDITION
	if(all->I->time==0.0){
		for(i=1;i<=all->P->chkpt->nrh;i++){
			r=all->P->rc->co[i][1];
			c=all->P->rc->co[i][2];
			sy=all->S->type->co[r][c];
			for(l=1;l<=Nl;l++){
				all->S->output[0][0][l-1]=all->S->T->co[l][r][c];// Tmin
				all->S->output[0][1][l-1]=all->S->T->co[l][r][c];// Tmax
				all->S->output[0][2][l-1]=all->S->T->co[l][r][c];// Tmean
				//printf("time=%f, T[%ld]=%f",times->time, l, all->S->T->co[l][r][c]); stop_execution();
				all->S->output[0][3][l-1]=teta_psi(all->S->P->co[l][r][c],all->S->thice->co[l][r][c],
						all->S->pa->co[sy][jsat][l],all->S->pa->co[sy][jres][l],all->S->pa->co[sy][ja][l],all->S->pa->co[sy][jns][l],1-1/all->S->pa->co[sy][jns][l],
						PSImin,all->P->Esoil);//Theta_w min
				all->S->output[0][4][l-1]=teta_psi(all->S->P->co[l][r][c],all->S->thice->co[l][r][c],
						all->S->pa->co[sy][jsat][l],all->S->pa->co[sy][jres][l],all->S->pa->co[sy][ja][l],all->S->pa->co[sy][jns][l],1-1/all->S->pa->co[sy][jns][l],
						PSImin,all->P->Esoil);//Theta_w max
				all->S->output[0][5][l-1]=teta_psi(all->S->P->co[l][r][c],all->S->thice->co[l][r][c],
						all->S->pa->co[sy][jsat][l],all->S->pa->co[sy][jres][l],all->S->pa->co[sy][ja][l],all->S->pa->co[sy][jns][l],1-1/all->S->pa->co[sy][jns][l],
						PSImin,all->P->Esoil);//Theta_w mean

				all->S->output[0][6][l-1]=all->S->thice->co[l][r][c];//Theta_i min
				all->S->output[0][7][l-1]=all->S->thice->co[l][r][c];//Theta_i max
				all->S->output[0][8][l-1]=all->S->thice->co[l][r][c];//Theta_i mean
				}
			}
		}

	energy_balance_superfast(all->I, all->P, all->L, all->T, all->S, all->M, all->W, all->E, all->N, all->G);

	//printf("\n ENERGY END %10.2f\n",all->I->time);
	//stop_execution();

	//printf("\n MASS START %10.2f\n",all->I->time);
	//stop_execution();
	water_balance_1D(all->T, all->S, all->L, all->W, all->C, all->P, all->I->time);

	//printf("\n MASS END %10.2f\n",all->I->time);
	//stop_execution();

	//printf("\n WRITE OUTPUT b %10.2f\n",all->I->time);
	//stop_execution();

	write_output_superfast(all->I, all->W, all->C, all->P, all->T, all->L, all->S, all->E, all->N, all->G, all->M);

	//printf("\n WRITE OUTPUT e %10.2f\n",all->I->time);
	//stop_execution();

	//Increase TIME!
	all->I->time+=(long)all->P->Dt;

 }while(all->I->time<=(all->I->TH*3600.0));/*end of time-cycle*/

 free(all->I);

}


void init_structs(TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet,
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times)
{
	//ENERGY struct
	egy->Rn_mean = NULL;
	egy->Rn_max = NULL;
	egy->Rn_min = NULL;
	egy->LW_in = NULL;
	egy->LW_out = NULL;
	egy->LW_max = NULL;
	egy->LW_min = NULL;
	egy->SW = NULL;
	egy->SW_max = NULL;
	egy->ET_mean = NULL;
	egy->ET_max = NULL;
	egy->ET_min = NULL;
	egy->H_mean = NULL;
	egy->H_max = NULL;
	egy->H_min = NULL;
	egy->G_mean = NULL;
	egy->G_max = NULL;
	egy->G_min = NULL;
	egy->G_snowsoil = NULL;
	egy->Ts_mean = NULL;
	egy->Ts_max = NULL;
	egy->Ts_min = NULL;
	egy->Rswdown_mean = NULL;
	egy->Rswdown_max = NULL;
	egy->Ta_mean = NULL;
	egy->Ta_max = NULL;
	egy->Ta_min = NULL;
	egy->out1 = NULL;
	egy->out2 = NULL;
	egy->out3 = NULL;
	egy->Rswbeam = NULL;
	egy->nDt_shadow = NULL;
	egy->nDt_sun = NULL;
	egy->Hgplot = NULL;
	egy->LEgplot = NULL;
	egy->Hvplot = NULL;
	egy->LEvplot = NULL;
	egy->SWinplot = NULL;
	egy->SWgplot = NULL;
	egy->SWvplot = NULL;
	egy->LWinplot = NULL;
	egy->LWgplot = NULL;
	egy->LWvplot = NULL;
	egy->Tgplot = NULL;
	egy->Tvplot = NULL;
	egy->Tsplot = NULL;
	egy->SWin = NULL;
	egy->LWin = NULL;
	egy->Hgrid = NULL;
	egy->Tsgrid = NULL;
	egy->VSFA = NULL;
	egy->HSFA = NULL;

	//SOIL struct
	sl->type = NULL;
	sl->pa = NULL;
	sl->P = NULL;
	sl->T = NULL;
	sl->thice = NULL;
	sl->Jinf = NULL;
	sl->J = NULL;
	sl->Tv = NULL;
	sl->Tav = NULL;
	sl->thwav = NULL;
	sl->thiav = NULL;
	sl->Tmean = NULL;
	sl->thetaw_mean = NULL;
	sl->thetai_mean = NULL;
	sl->psi_mean = NULL;
	sl->Tmax = NULL;
	sl->thetaw_max = NULL;
	sl->thetai_max = NULL;
	sl->Tmin = NULL;
	sl->thetaw_min = NULL;
	sl->thetai_min = NULL;
	sl->output = NULL;
	sl->bc = NULL;
	sl->ET = NULL;

	//TOPO struct
	top->Z0 = NULL;         
	top->Z1 = NULL;
	top->Z0dp = NULL;		  
	top->Z0ext = NULL;      
	top->sky = NULL;        
	top->pixel_type = NULL; 
	top->DD = NULL;         
	top->i_DD = NULL;       
	top->dz_dx = NULL;      
	top->dz_dy = NULL;      
	top->top_index = NULL; 
	top->curv = NULL;       
	top->area = NULL;       
	top->aspect = NULL;     
	top->slopes = NULL;     
	top->i_ch = NULL;       
	top->pixel_distance = NULL;
	top->ES_pixel = NULL;
	top->ES_aspect = NULL;
	top->ES_slope = NULL;
	top->horizon_height = NULL;
	top->Zm = NULL;
	top->curv_m = NULL;
	top->slope_m = NULL;
	top->slopeaz_m = NULL;
	top->i_cont = NULL;
	top->lrc_cont = NULL;
	top->Z = NULL;
	top->slope_H = NULL;

	//LAND struct
	land->LC = NULL;
	land->LC2 = NULL;
	land->albedo = NULL;
	land->shadow = NULL;
	land->clax = NULL;
	land->cont = NULL;
	land->ty = NULL;

	//CHANNEL struct
	cnet->r = NULL;
	cnet->c = NULL;
	cnet->ch = NULL;
	cnet->Q = NULL;
	cnet->s0 = NULL;
	cnet->fraction_spread = NULL;
	cnet->Q_sup_s = NULL;
	cnet->Q_sub_s = NULL;
	cnet->Qsup_spread = NULL;
	cnet->Qsub_spread = NULL;
	cnet->Qsup = NULL;
	cnet->Qsub = NULL;

	//WATER struct
	wat->weights_Kriging = NULL;
	wat->q_sub = NULL;          
	wat->q_sup = NULL;
	wat->h_sup = NULL;
	wat->total = NULL;
	wat->Pn = NULL;
	wat->wcan_rain = NULL;
	wat->wcan_snow = NULL;
	wat->PrTOT_mean = NULL;
	wat->PrSNW_mean = NULL;
	wat->out1 = NULL;
	wat->out2 = NULL;
	wat->hsupav = NULL;
	wat->outfluxes = NULL;
	wat->error = NULL;

	//PAR struct
	par->Dmin = NULL;
	par->Dmax = NULL;
	par->Dmin_glac = NULL;
	par->Dmax_glac = NULL;
	par->chkpt = NULL;
	par->rc = NULL;
     par->saving_points = NULL;
     par->JD_plots = NULL; 
	par->r_points = NULL;
	par->c_points = NULL;
	par->cont_trans = NULL;
	par->ibeg = NULL;
	par->transect = NULL;
	par->vtrans = NULL;

	//SNOW struct
	snow->type = NULL;
	snow->lnum = NULL;
	snow->Dzl = NULL;
	snow->w_liq = NULL;
	snow->w_ice = NULL;
	snow->T = NULL;
	snow->nondimens_age = NULL;
	snow->dimens_age = NULL;
	snow->evap = NULL;
	snow->subl = NULL;
	snow->melted = NULL;
	snow->max = NULL;
	snow->average = NULL;
	snow->MELTED = NULL;
	snow->SUBL = NULL;
	snow->t_snow = NULL;
	snow->totav_snow = NULL;
	snow->DDF = NULL;
	snow->DDF1 = NULL;
	snow->DDFvar = NULL;
	snow->DDFcont = NULL;
	snow->DDFTmin = NULL;
	snow->DDFmeltTL0 = NULL;
	snow->DDFmelt = NULL;
	snow->rho_newsnow = NULL;
	snow->Qsub = NULL;
	snow->Wtrans = NULL;
	snow->Qtrans = NULL;
	snow->Qtrans_x = NULL;
	snow->Qtrans_y = NULL;
	snow->Wtot = NULL;
	snow->Wsubl_cum = NULL;
	snow->Wsusp_cum = NULL;
	snow->Wsalt_cum = NULL;
	snow->Wsubgrid_cum = NULL;
	snow->out_bs = NULL;
	snow->ListonSWE = NULL;
	snow->softSWE = NULL;
	snow->softSWE1 = NULL;
	snow->Dplot = NULL;
	snow->CR1 = NULL;
	snow->CR2 = NULL;
	snow->CR3 = NULL;
	snow->CR1m = NULL;
	snow->CR2m = NULL;
	snow->CR3m = NULL;
	snow->change_dir_wind = NULL;
	snow->Psnow = NULL;
	snow->rhoSOFT = NULL;

	//GALCIER struct
	glac->lnum = NULL;
	glac->Dzl = NULL;
	glac->w_liq = NULL;
	glac->w_ice = NULL;
	glac->T = NULL;
	glac->evap = NULL;
	glac->subl = NULL;
	glac->melted = NULL;
	glac->MELTED = NULL;
	glac->SUBL = NULL;
	glac->DDF = NULL;
	glac->DDF1 = NULL;
	glac->DDFvar = NULL;
	glac->DDFcont = NULL;
	glac->DDFTmin = NULL;
	glac->DDFmeltTL0 = NULL;
	glac->DDFmelt = NULL;

	//METEO struct
	met->st = NULL;
	met->data = NULL;
	met->column = NULL;
	met->horizon = NULL;
	met->var = NULL;
	met->LRs = NULL;
	met->LRv = NULL;
	met->LRp = NULL;
	met->Tgrid = NULL;
	met->Pgrid = NULL;
	met->Vgrid = NULL;
	met->Vdir = NULL;
	met->RHgrid = NULL;
	met->Vspdmean = NULL;
	met->Vdirmean = NULL;
	met->RHmean = NULL;
	met->Taplot = NULL;
	met->Vspdplot = NULL;
	met->Vdirplot = NULL;
	met->RHplot = NULL;
	met->Tday = NULL;
	met->Tvar = NULL;
}
