
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
    
    
    /*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/

#include <sys/stat.h>
#include "struct.geotop.h"
#include "input.h"
#include "output.h"
#include "times.h"
#include "constants.h"
#include "energy.balance.h"
#include "meteo.h"
#include "water.balance.h"
#include "snow.h"
#include "blowingsnow.h"
#include "tabs.h"
#include "deallocate.h"
#include "pedo.funct.h"

#include <time.h>

void time_loop(ALLDATA *A);

               
/*----------   1. Global variables  ------------*/

#include "keywords.h"	//contains the definition of char** keywords_num and char** keywords_char

long number_novalue;
long number_absent;
char *string_novalue;

T_INIT *UV;

char *logfile;
char **files;

long Nl,Nr,Nc;
double t_meteo, t_energy, t_water, t_sub, t_sup, t_blowingsnow, t_out;

double **odpnt, **odp;
long *opnt, nopnt;
short *ipnt, *ibsn;
char **hpnt;

double *odbsn, *odb;
long *obsn, nobsn;
char **hbsn;

long *osnw, nosnw;
char **hsnw;

long *oglc, noglc;
char **hglc;

long *osl, nosl;
char **hsl;

FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav, *ffice, *fficeav, *ffsnowT, *ffsnowl, *ffsnowi, *ffsnowd, *ffglac;

long i_sim=0, i_run, i_sim0, i_run0;

char *SuccessfulRunFile, *FailedRunFile;

time_t start_time; 
double elapsed_time, elapsed_time_start, cum_time, max_time;

double MM1, MM2, MMR, MMo, MS1, MS2;

/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/

int main(int argc,char *argv[]){
	
	ALLDATA *adt;
	FILE *f;
	
	//assign novalues
	number_novalue = -9999;
	number_absent = -9998;
	string_novalue = assign_string("none");
	i_sim0 = 1;
	i_run0 = 1;
	cum_time = 0.;
	elapsed_time_start = 0.;
	
	MM1=0.;
	MM2=0.;
	MMR=0.;
	MMo=0.;
   
	/*dinamic allocations:*/
	UV=(T_INIT *)malloc(sizeof(T_INIT));
	if(!UV) t_error("UV was not allocated");
 
	adt=(ALLDATA *)malloc(sizeof(ALLDATA));
	if(!adt){
		t_error("adt was not allocated");
	}else {
		
		adt->I=(TIMES *)malloc(sizeof(TIMES));
		if(!(adt->I)) t_error("times was not allocated");	
		
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
		
		t_meteo=0.;
		t_energy=0.; 
		t_water=0.;
		t_sub=0.;
		t_sup=0.;
		t_out=0.;
		t_blowingsnow=0.;
		
		
		/*------------------    3.  Acquisition of input data and initialisation    --------------------*/
		get_all_input(argc, argv, adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I);
		
		/*-----------------   4. Time-loop for the balances of water-mass and egy   -----------------*/
		time_loop(adt);
		
		/*--------------------   5.Completion of the output files and deallocaions  --------------------*/
		dealloc_all(adt->T, adt->S, adt->L, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->M, adt->I);
		free(adt);

	}
		
	printf("End of simulation!\n");
	
	f = fopen(SuccessfulRunFile, "w");
	fclose(f);
	
	return 0;
}


/*----------------   6. The most important subroutine of the main: "time_loop"   ---------------*/
void time_loop(ALLDATA *A){ 

	clock_t tstart, tend;
	short en=0, wt=0, out;
	long i, sy, r, c, j, l;
	double t, Dt, JD0, JDb, JDe, W, th, th0;
	double Vout, Voutsub, Voutsup, Vbottom, C0, C1;
	FILE *f;
	
	//double mean;
	
	STATEVAR_3D *S=NULL, *G=NULL;
	SOIL_STATE *L, *C;
	STATE_VEG *V;
	DOUBLEVECTOR *a, *Vsup_ch, *Vsub_ch;
	
	
	S=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
	allocate_and_initialize_statevar_3D(S, (double)number_novalue, A->P->max_snow_layers, Nr, Nc);
	if(A->P->max_glac_layers>0){
		G=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
		allocate_and_initialize_statevar_3D(G, (double)number_novalue, A->P->max_glac_layers, Nr, Nc);
	}
	L=(SOIL_STATE *)malloc(sizeof(SOIL_STATE));
	initialize_soil_state(L, A->P->total_pixel, Nl);
	C=(SOIL_STATE *)malloc(sizeof(SOIL_STATE));
	initialize_soil_state(C, A->C->r->nh, Nl);	
	V=(STATE_VEG *)malloc(sizeof(STATE_VEG));
	initialize_veg_state(V, A->P->total_pixel);
	a=new_doublevector(A->P->total_pixel);
	Vsub_ch=new_doublevector(A->C->r->nh);
	Vsup_ch=new_doublevector(A->C->r->nh);
	
	time( &start_time );

	//periods
	i_sim = i_sim0;
	
	do{
	
		//runs
		A->I->time = A->P->delay_day_recover*86400.;//Initialize time	
		A->P->delay_day_recover = 0.;
		
		do{
			
			if( A->I->time > (A->P->end_date->co[i_sim] - A->P->init_date->co[i_sim])*86400. - 1.E-5){
				printf("Number of times the simulation #%ld has been run: %ld\n",i_sim,i_run);
				f=fopen(logfile, "a");
				fprintf(f,"Number of times the simulation #%ld has been run: %ld\n",i_sim,i_run);
				fclose(f);
				
				print_run_average(A->S, A->T, A->P);
				
				i_run++;
				A->I->time = 0.0;//Initialize time
				
				A->M->line_interp_WEB_LR = 0;
				A->M->line_interp_Bsnow_LR = 0;
				for (i=1; i<=A->M->st->Z->nh; i++) {
					A->M->line_interp_WEB[i-1] = 0;
					A->M->line_interp_Bsnow[i-1] = 0;
				}
				
				if(i_run <= A->P->run_times->co[i_sim]){
					reset_to_zero(A->P, A->S, A->L, A->N, A->G, A->E, A->M, A->W);
					init_run(A->S, A->P);
				}
				
			}else {
				
				//find time step from file or inpts
				set_time_step(A->P, A->I);
				
				//time at the beginning of the time step
				JD0 = A->P->init_date->co[i_sim]+A->I->time/secinday;			
				
				//time step variables
				t = 0.;
				Dt = A->P->Dt;
				
				//time step subdivisions
				do{
					
					JDb = A->P->init_date->co[i_sim]+(A->I->time+t)/secinday;
					
					if (t + Dt > A->P->Dt) Dt = A->P->Dt - t;
					
					//iterations
					do{
						
						JDe = A->P->init_date->co[i_sim]+(A->I->time+t+Dt)/secinday;
						
						//copy state variables on 
						copy_snowvar3D(A->N->S, S);
						copy_doublevector(A->N->age, a);
						if (A->P->max_glac_layers>0) copy_snowvar3D(A->G->G, G);
						copy_soil_state(A->S->SS, L);
						copy_soil_state(A->C->SS, C);
						copy_veg_state(A->S->VS, V);	
						
						/*for (j=1; j<=A->W->H1->nh; j++) {
							l=A->T->lrc_cont->co[j][1];
							r=A->T->lrc_cont->co[j][2];
							c=A->T->lrc_cont->co[j][3];
							printf("START %ld %ld %ld %e\n",l,r,c,A->S->SS->P->co[l][A->T->j_cont[r][c]]);
						}*/
						
						//init
						initialize_doublevector(Vsub_ch, 0.);
						initialize_doublevector(Vsup_ch, 0.);			
						Vout = 0.;
						Voutsub = 0.;
						Voutsup = 0.;
						Vbottom = 0.;
						
						//meteo
						tstart=clock();
						meteo_distr(A->M->line_interp_WEB, A->M->line_interp_WEB_LR, A->M, A->W, A->T, A->P, JD0, JDb, JDe);
						tend=clock();
						t_meteo+=(tend-tstart)/(double)CLOCKS_PER_SEC;
						
						if(A->P->en_balance == 1){
							tstart=clock();
							en = EnergyBalance(Dt, JD0, JDb, JDe, L, C, S, G, V, a, A, &W);
							tend=clock();
							t_energy+=(tend-tstart)/(double)CLOCKS_PER_SEC;
						}
						
						if(A->P->wat_balance == 1 && en == 0){
							tstart=clock();
							wt = water_balance(Dt, JD0, JDb, JDe, L, C, A, Vsub_ch, Vsup_ch, &Vout, &Voutsub, &Voutsup, &Vbottom);
							tend=clock();
							t_water+=(tend-tstart)/(double)CLOCKS_PER_SEC;
						}


						if (en != 0 || wt != 0) {
							
							if(Dt > A->P->min_Dt) Dt *= 0.5;
							out = 0;
							
							f = fopen(logfile, "a");
							if (en != 0) {
								fprintf(f,"Energy balance not converging\n");
							}else {
								fprintf(f,"Water balance not converging\n");
							}
							fprintf(f,"Reducing time step to %f s, t:%f s\n",Dt,t);
							fclose(f);
							
						}else {
							out = 1;
						}
						
						//printf("Dt:%f min:%f\n",Dt,A->P->min_Dt);
						
					}while( out == 0 && Dt > A->P->min_Dt ); 
					
					/*if (en != 0 || wt != 0) {
						f = fopen(FailedRunFile, "w");
						fprintf(f, "Simulation Period:%ld\n",i_sim);
						fprintf(f, "Run Time:%ld\n",i_run);
						fprintf(f, "Number of days after start:%f\n",A->I->time/86400.);	
						
						if (en != 0 && wt == 0) {
							fprintf(f, "ERROR: Energy balance does not converge, Dt:%f\n",Dt);
						}else if (en == 0 && wt != 0) {
							fprintf(f, "ERROR: Water balance does not converge, Dt:%f\n",Dt);
						}else {
							fprintf(f, "ERROR: Water and energy balance do not converge, Dt:%f\n",Dt);
						}
						
						fclose(f);
						t_error("Fatal Error! Geotop is closed. See failing report.");	
					}*/
					
					if (en != 0 || wt != 0) {
						//f = fopen(FailedRunFile, "w");
						
						f = fopen(logfile, "a");
						//fprintf(f, "Simulation Period:%ld\n",i_sim);
						//fprintf(f, "Run Time:%ld\n",i_run);
						//fprintf(f, "Number of days after start:%f\n",A->I->time/86400.);	
						
						if (en != 0 && wt == 0) {
							fprintf(f, "WARNING: Energy balance does not converge, Dt:%f\n",Dt);
						}else if (en == 0 && wt != 0) {
							fprintf(f, "WARNING: Water balance does not converge, Dt:%f\n",Dt);
						}else {
							fprintf(f, "WARNING: Water and energy balance do not converge, Dt:%f\n",Dt);
						}
						
						fclose(f);
						//t_error("Fatal Error! Geotop is closed. See failing report.");	
					}
					
					t += Dt;
					
					if (A->P->state_pixel == 1 && A->P->dUzrun == 1) {
						for (j=1; j<=A->P->rc->nrh; j++) {
							for (l=1; l<=Nl; l++){
								r = A->P->rc->co[j][1];
								c = A->P->rc->co[j][2];
								sy = A->S->type->co[r][c];
								
								th = theta_from_psi(A->S->SS->P->co[l][A->T->j_cont[r][c]], A->S->SS->thi->co[l][A->T->j_cont[r][c]], l, A->S->pa->co[sy], PsiMin);
								if(th > A->S->pa->co[sy][jsat][l]-A->S->SS->thi->co[l][A->T->j_cont[r][c]]) th = A->S->pa->co[sy][jsat][l]-A->S->SS->thi->co[l][A->T->j_cont[r][c]];
								C0 = A->S->pa->co[sy][jct][l]*(1.-A->S->pa->co[sy][jsat][l])*A->S->pa->co[sy][jdz][l] + c_ice*A->S->SS->thi->co[l][A->T->j_cont[r][c]] + c_liq*th;
								th0 = th;
								
								th = theta_from_psi(L->P->co[l][A->T->j_cont[r][c]], L->thi->co[l][A->T->j_cont[r][c]], l, A->S->pa->co[sy], PsiMin);
								if(th > A->S->pa->co[sy][jsat][l]-L->thi->co[l][A->T->j_cont[r][c]]) th = A->S->pa->co[sy][jsat][l]-L->thi->co[l][A->T->j_cont[r][c]];
								C1 = A->S->pa->co[sy][jct][l]*(1.-A->S->pa->co[sy][jsat][l])*A->S->pa->co[sy][jdz][l] + c_ice*L->thi->co[l][A->T->j_cont[r][c]] + c_liq*th;
								
								A->S->dUzrun->co[j][l] += 1.E-6*( 0.5*(C0+C1)*(L->T->co[l][A->T->j_cont[r][c]] - A->S->SS->T->co[l][A->T->j_cont[r][c]]) + Lf*(th-th0)*A->S->pa->co[sy][jdz][l] );
							}
						}
					}
					
					//write state variables
					copy_snowvar3D(S, A->N->S);
					copy_doublevector(a, A->N->age);
					if (A->P->max_glac_layers>0) copy_snowvar3D(G, A->G->G);
					copy_soil_state(L, A->S->SS);
					copy_soil_state(C, A->C->SS);
					copy_veg_state(V, A->S->VS);
					add_doublevector(Vsub_ch, A->C->Vsub);
					add_doublevector(Vsup_ch, A->C->Vsup);	
					A->C->Vout += Vout;
					A->W->Voutbottom += Vbottom;
					A->W->Voutlandsub += Voutsub;
					A->W->Voutlandsup += Voutsup;
					
					//printf("%f\n",A->I->time);
					
					//record time step
					odb[ootimestep] = Dt * (Dt/A->P->Dtplot_basin->co[i_sim]);
					
					//write output variables
					fill_output_vectors(Dt, W, A->E, A->N, A->G, A->W, A->M, A->P, A->I, A->T, A->S);
					
					//reset Dt
					if (Dt < A->P->Dt) Dt *= 2.;
					
				}while(t < A->P->Dt);
				
				if(A->P->blowing_snow==1){
					tstart=clock();
					windtrans_snow(A->N, A->M, A->W, A->L, A->T, A->P, A->I->time);
					tend=clock();
					t_blowingsnow+=(tend-tstart)/(double)CLOCKS_PER_SEC;
				}

				tstart=clock();
				write_output(A->I, A->W, A->C, A->P, A->T, A->L, A->S, A->E, A->N, A->G, A->M);
				tend=clock();
				t_out+=(tend-tstart)/(double)CLOCKS_PER_SEC;
				
				A->I->time += A->P->Dt;//Increase TIME
								
			}

		}while(i_run <= A->P->run_times->co[i_sim]);//end of time-cycle
						
		if (A->P->newperiodinit != 0) end_period_1D(A->S, A->T, A->P);
		if (i_sim < A->P->init_date->nh) change_grid(i_sim, i_sim+1, A->P, A->T, A->L, A->W, A->C);
		
		reset_to_zero(A->P, A->S, A->L, A->N, A->G, A->E, A->M, A->W);
		init_run(A->S, A->P);
		
		i_sim++;
		i_run0 = 1;
		i_run = i_run0;
				
	}while (i_sim <= A->P->init_date->nh);
	
	deallocate_statevar_3D(S);
	if(A->P->max_glac_layers>0) deallocate_statevar_3D(G);
	deallocate_soil_state(L);
	deallocate_soil_state(C);
	deallocate_veg_state(V);
	free_doublevector(a);
	free_doublevector(Vsub_ch);
	free_doublevector(Vsup_ch);

}

/*--------------------------------------------------------------------------------------------*/

