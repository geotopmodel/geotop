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

/*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/

#ifdef USE_HPC

#include <mpi.h>

#endif

#include <sys/stat.h>
#include "struct.geotop.h"
#include "hpc.geotop.h"
#include "input.h"
#include "output.h"
#include "times.h"
#include "constants.h"
#include "energy.balance.h"
#include "meteo.h"
#include "water.balance.h"
#include "snow.h"
#include "blowingsnow.h"
#include "../libraries/ascii/tabs.h"
#include "deallocate.h"

#include "../gt_utilities/gt_utilities.h"
#include "../gt_utilities/gt_symbols.h"
#include "../gt_utilities/ncgt_output.h"

#include <time.h>
#include "output_nc.h"
void time_loop(ALLDATA *all);

/*----------   1. Global variables  ------------*/

#include "keywords.h"	//contains the definition of char** keywords_num and char** keywords_char
long number_novalue;
long number_absent;
char *string_novalue;

T_INIT *UV;

char *logfile;
char **files;

long Nl, Nr, Nc;
double t_meteo, t_energy, t_water, t_sub, t_sup, t_blowingsnow, t_out;

double **outdata_point;
long *outputpoint, noutputpoint;
char **headerpoint;

double *outdata_basin;
long *outputbasin, noutputbasin;
char **headerbasin;

long *outputsnow, noutputsnow;
char **headersnow;

long *outputglac, noutputglac;
char **headerglac;

long *outputsoil, noutputsoil;
char **headersoil;

FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav,
		*ffice, *fficeav, *ffsnow, *ffglac;

long i_sim, i_run;

/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/

int main(int argc, char *argv[]) {

	//assign novalues
	number_novalue = -9999;
	number_absent = -9998;
	string_novalue = assign_string("none");
	i_sim = 1;
	i_run = 1;

	/*structs' declarations:*/
	ALLDATA *adt;

#ifdef USE_HPC

	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_comm comm;
	MPI_info info;

	std::cout << "Process " << getpid() << " is " << myrank << " of " << nprocs << " processes" << std::endl;

	// linked list for parallel partition structure record
	GCSTRUCT *start = new GCSTRUCT; // first element will NOT handle any info
	start->next = NULL;
	WORKAREA* rankArea = new WORKAREA;

	getpartitions(start, rankarea);

	// for debug only
	if (cprocs != nprocs) {
		std::cout << "The number of process (" << nprocs << ") is different from the one required (" << cprocs << ")" << std::endl;
	} else {
		std::cout << "The number of active processes (" << nprocs << ") is equal to the one calculated (" << cprocs << ")" << std::endl;
	}

#endif

	/*dinamic allocations:*/
	UV = (T_INIT *) malloc(sizeof(T_INIT));
	if (!UV)
		t_error("UV was not allocated");

	adt = (ALLDATA *) malloc(sizeof(ALLDATA));
	if (!adt) {
		t_error("adt was not allocated");
	} else {

		adt->I = (TIMES *) malloc(sizeof(TIMES));
		if (!(adt->I))
			t_error("times was not allocated");

		adt->T = (TOPO *) malloc(sizeof(TOPO));
		if (!(adt->T))
			t_error("top was not allocated");

		adt->S = (SOIL *) malloc(sizeof(SOIL));
		if (!(adt->S))
			t_error("sl was not allocated");

		adt->L = (LANDCOVER *) malloc(sizeof(LANDCOVER));
		if (!(adt->L))
			t_error("land was not allocated");

		adt->W = (WATER *) malloc(sizeof(WATER));
		if (!(adt->W))
			t_error("water was not allocated");

		adt->P = (PAR *) malloc(sizeof(PAR));
		if (!(adt->P))
			t_error("par was not allocated");

		adt->C = (CHANNEL *) malloc(sizeof(CHANNEL));
		if (!(adt->C))
			t_error("channel was not allocated");

		adt->E = (ENERGY *) malloc(sizeof(ENERGY));
		if (!(adt->E))
			t_error("egy was not allocated");

		adt->N = (SNOW *) malloc(sizeof(SNOW));
		if (!(adt->N))
			t_error("snow was not allocated");

		adt->G = (GLACIER *) malloc(sizeof(GLACIER));
		if (!(adt->G))
			t_error("glac was not allocated");

		adt->M = (METEO *) malloc(sizeof(METEO));
		if (!(adt->M))
			t_error("met was not allocated");

		t_meteo = 0.;
		t_energy = 0.;
		t_water = 0.;
		t_sub = 0.;
		t_sup = 0.;
		t_out = 0.;
		t_blowingsnow = 0.;

#ifdef USE_NETCDF

		int ncid=ncgt_open_from_option_string(argc,argv,NC_GEOTOP_ARCHIVE_OPTION,NC_GEOTOP_NODEFINE,GEOT_VERBOSE);
		adt->ncid=ncid;

#endif

		/*------------------    3.  Acquisition of input data and initialization    --------------------*/
#ifdef USE_HPC
		get_all_input(argc, argv, adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I, rankArea);
#else
		get_all_input(argc, argv, adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I);
#endif

#ifdef USE_NETCDF
#ifdef USE_HPC
		set_output_nc(adt, rankArea);
#else
		set_output_nc(adt);
#endif
#endif
		/*-----------------   4. Time-loop for the balances of water-mass and egy   -----------------*/
#ifdef USE_HPC
		time_loop(adt, start, rankArea);
#else
		time_loop(adt);
#endif

#ifdef USE_NETCDF

		ncgt_close_geotop_archive(ncid);
		deallocate_output_nc(adt->outnc);
#endif

		/*--------------------   5.Completion of the output files and deallocations  --------------------*/
		dealloc_all(adt->T, adt->S, adt->L, adt->W, adt->C, adt->P, adt->E,
				adt->N, adt->G, adt->M, adt->I);

		free(adt);

	}

	printf("End of simulation!\n");

	return 0;
}

/*----------------   6. The most important subroutine of the main: "time_loop"   ---------------*/

#ifdef USE_HPC
void time_loop(ALLDATA *all, GCSTRUCT *start, WORKAREA *rankArea)
#else
void time_loop(ALLDATA *all)
#endif
{

	clock_t start, end;

	long i;

	FILE *f;

#ifdef USE_HPC
	//Nr=top->Z0->nrh;
	//Nc=top->Z0->nch;

	// define start and count values to control dimension_x and dimension_y for parallel write
	// better to place here ouside the loop
	Nr = abs(rankArea->top - rankArea->bottom);
	Nc = abs(rankArea->left - rankArea->right);
	if (rankArea->top == 1) {
		Nr = Nr + 1;
		offsetNr = rankArea->top;
	} else if (rankArea->bottom == top->Z0->nrh) {
		Nr = Nr + 1;
		offsetNr = rankArea->top - 1;
	} else {
		offsetNr = rankArea->top - 1;
		Nr = Nr + 2;
	}
	if (rankArea->left == 1) {
		Nc = Nc + 1;
		offsetNc = rankArea->left;
	} else if (rankArea->right == top->Z0->nch) {
		Nc = Nc + 1;
		offsetNc = rankArea->left - 1;
	} else {
		Nc = Nc + 2;
		offsetNc = rankArea->left;
	}
#endif


#ifdef USE_NETCDF
	long nsim=1;
	// TODO: check number of simulation issue
#else
	long nsim = all->P->init_date->nh;
#endif
	for (i_sim = 1; i_sim <= nsim; i_sim++) {

#ifdef USE_NETCDF
		long nrun=1;
		// TODO: check number run issue
#else
		long nrun = all->P->run_times->co[i_sim];
#endif
		i_run = 1;//Run index

		all->I->time = 0.0;//Initialize time	
#ifdef USE_NETCDF
		// printing time counter initialization
		//TODO please revisit into details the initialization procedure for counters according to time coordinates
		all->counter_surface_energy=0;
		all->counter_snow=0;
		all->counter_soil=0;
		all->counter_glac=0;
		all->counter_point=0;
		if(all->P->point_sim!=1) { //distributed simulation
			all->point_var_type=NC_GEOTOP_2D_MAP_IN_CONTROL_POINT;
			all->z_point_var_type=NC_GEOTOP_3D_MAP_IN_CONTROL_POINT;
		} else {// point simulation
			all->z_point_var_type=NC_GEOTOP_Z_POINT_VAR;
			all->point_var_type=NC_GEOTOP_POINT_VAR;
		}
#endif
		do {

			start = clock();
			set_time_step(all->P, all->I);

			meteo_distr(1, 1, all->M->line_interp_WEB,
					all->M->line_interp_WEB_LR, all->M, all->W, all->T, all->P,
					all->P->init_date->co[i_sim] + all->I->time / 86400.,
					all->P->init_date->co[i_sim] + (all->I->time + all->P->Dt)
							/ 86400.);
			end = clock();
			t_meteo += (end - start) / (double) CLOCKS_PER_SEC;

			if (all->P->en_balance == 1) {
				start = clock();
				energy_balance(all->I, all->P, all->L, all->T, all->S, all->M,
						all->W, all->E, all->N, all->G, all->C);
				end = clock();
				t_energy += (end - start) / (double) CLOCKS_PER_SEC;
			}

			if (all->P->wat_balance == 1) {
				start = clock();
				water_balance(all);
				end = clock();
				t_water += (end - start) / (double) CLOCKS_PER_SEC;
			}

			if (all->P->blowing_snow == 1) {
				start = clock();
				windtrans_snow(all->N, all->M, all->W, all->L, all->T, all->P,
						all->I->time);
				if (strcmp(files[fSCA], string_novalue) != 0)
					find_SCA(all->N->S, all->P, all->L->LC->co, all->I->time
							+ all->P->Dt);
				end = clock();
				t_blowingsnow += (end - start) / (double) CLOCKS_PER_SEC;

			}

			start = clock();
#ifdef USE_NETCDF
			if(all->ncid==NC_GEOTOP_MISSING) {// there is no GEOtop netCDF archive, then use ascii modality
				write_output(all->I, all->W, all->C, all->P, all->T, all->L, all->S, all->E, all->N, all->G, all->M);
			}
			else {
#ifdef USE_HPC
				write_output_nc(all, rankArea);
				updateGhostcells(all, start);
#else
				write_output_nc(all);
#endif
			}
#else
			write_output(all->I, all->W, all->C, all->P, all->T, all->L, all->S, all->E, all->N, all->G, all->M);
#endif
			end = clock();
			t_out += (end - start) / (double) CLOCKS_PER_SEC;

			all->I->time += all->P->Dt;//Increase TIME	

			if (all->I->time > (all->P->end_date->co[i_sim]
					- all->P->init_date->co[i_sim]) * 86400. - 1.E-5) {
				printf(
						"Number of times the simulation #%ld has been run: %ld\n",
						i_sim, i_run);
				f = fopen(logfile, "a");
				fprintf(
						f,
						"Number of times the simulation #%ld has been run: %ld\n",
						i_sim, i_run);
				fclose(f);

				i_run++;
				all->I->time = 0.0;//Initialize time

				all->M->line_interp_WEB_LR = 0;
				all->M->line_interp_Bsnow_LR = 0;
				for (i = 1; i <= all->M->st->Z->nh; i++) {
					all->M->line_interp_WEB[i - 1] = 0;
					all->M->line_interp_Bsnow[i - 1] = 0;
				}
			}

		} while (i_run <= nrun);//end of time-cycle

		reset_to_zero(all->P, all->S, all->L, all->N, all->G, all->E, all->M, all->W);

	}
}

/*--------------------------------------------------------------------------------------------*/

