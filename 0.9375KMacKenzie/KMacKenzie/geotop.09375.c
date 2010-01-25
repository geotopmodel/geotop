
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

void time_loop(ALLDATA *all);
void time_loop_superfast(ALLDATA *all);


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
