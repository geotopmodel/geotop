
/* STATEMENT:
 water_balance
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
    
    
    /*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/

#include "config.h"

#include <sys/stat.h>
#include "struct.geotop.h"
#include "input.h"
#include "output.h"
#include "times.h"
#include "constants.h"
#include "energy.balance.h"
#include "water.balance.h"
#include "meteo.h"
#include "snow.h"
#include "blowingsnow.h"
//#include "../libraries/ascii/tabs.h"
#include "deallocate.h"
#include <meteoio/MeteoIO.h>
#include <string>
#ifdef USE_NETCDF
//#include "../gt_utilities/gt_utilities.h"
//#include "../gt_utilities/gt_symbols.h"
//#include "../gt_utilities/ncgt_output.h"
#include "../netCDF/netcdfIO.h"
#include "../netCDF/read_command_line_netcdf.h"
#include "../netCDF/gt_symbols.h"
#include "output_nc.h"
#endif
#include <time.h>

using namespace std;

//void time_loop(ALLDATA *A, mio::IOManager& iomanager);
  void time_loop(AllData *A, mio::IOManager& iomanager);

               
/*----------   1. Global variables  ------------*/
// defined and declared into geotop.cc

/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/

int main(int argc,char *argv[]){

	clock_t start, end;
	double elapsed;
	start = clock();
	
//	ALLDATA *adt;
	AllData *adt;
	FILE *f;

	string cfgfile = "io_it.ini";
//	printf("argc=%d, argv[0]=%s argv[1]=%s argv[2]=%s argv[3]=%s argv[4]=%s",argc,argv[0],argv[1],argv[2],argv[3],argv[4]);stop_execution();
//	if (argc > 1) cfgfile = string(argv[1]) + "/" + cfgfile; //if a working directory is given, we prepend it to io_it.ini
    /* ANTONIO: This is just crazy. We have *two* incompatible way of
       calling geotop, and part of the command line parsing is done
       later on, in ncgt_open_from_option_string(), but only if
       USE_NETCDF is defined, so if its not defined some of the
       command line options are just silently ignored !

       We should use getopt() to parse command line arguments, store
       ignored command line options and pass them to
       ncgt_open_from_option_string().
    */
	if (argc == 2) cfgfile = string(argv[1]) + "/" + cfgfile; //if a working directory is given, we prepend it to io_it.ini
	else if (argc > 2) cfgfile = string(argv[2]) + "/" + cfgfile;
	mio::Config cfg(cfgfile);
	cfg.addKey("GRID2DPATH", "Input", ".");
	mio::IOManager iomanager(cfg);
	
	//assign novalues
	number_novalue = -9999;
	number_absent = -9998;
	string_novalue = assign_string("none");
	i_sim0 = 1;
	i_run0 = 1;
	cum_time = 0.;
	elapsed_time_start = 0.;
   
	/*dinamic allocations:*/
//	UV=(T_INIT *)malloc(sizeof(T_INIT));
	UV=new TInit();
	if(!UV) t_error("UV was not allocated");
 
//	adt=(ALLDATA *)malloc(sizeof(ALLDATA));
	adt=new AllData();

	if(!adt){
		t_error("adt was not allocated");
	}else {
		
	//	adt->I=(TIMES *)malloc(sizeof(TIMES));
		adt->I=new Times();
		if(!(adt->I)) t_error("times was not allocated");	
		
	//	adt->T=(TOPO *)malloc(sizeof(TOPO));
		adt->T=new Topo();
		if(!(adt->T)) t_error("top was not allocated");

	//	adt->S=(SOIL *)malloc(sizeof(SOIL));
		adt->S=new Soil();
		if(!(adt->S)) t_error("sl was not allocated");
		
	//	adt->L=(LAND *)malloc(sizeof(LAND));
		adt->L=new Land();
		if(!(adt->L)) t_error("land was not allocated");
		
	//	adt->W=(WATER *)malloc(sizeof(WATER));
		adt->W=new Water();
		if(!(adt->W)) t_error("water was not allocated");
		
	//	adt->P=(PAR *)malloc(sizeof(PAR));
		adt->P=new Par();
		if(!(adt->P)) t_error("par was not allocated");
		
	//	adt->C=(CHANNEL *)malloc(sizeof(CHANNEL));
		adt->C=new Channel();
		if(!(adt->C)) t_error("channel was not allocated");
		
	//	adt->E=(ENERGY *)malloc(sizeof(ENERGY));
		adt->E = new Energy();
		if(!(adt->E)) t_error("egy was not allocated");
		
	//	adt->N=(SNOW *)malloc(sizeof(SNOW));
		adt->N=new Snow();
		if(!(adt->N)) t_error("snow was not allocated");	
		
	//	adt->G=(GLACIER *)malloc(sizeof(GLACIER));
		adt->G = new Glacier();
		if(!(adt->G)) t_error("glac was not allocated"); 
		
	//	adt->M=(METEO *)malloc(sizeof(METEO));
		adt->M=new Meteo();
		if(!(adt->M)) t_error("met was not allocated"); 
		
		t_meteo=0.;
		t_energy=0.; 
		t_water=0.;
		t_sub=0.;
		t_sup=0.;
		t_out=0.;
		t_blowingsnow=0.;

#ifdef USE_NETCDF
		int ncid=ncgt_open_from_option_string(argc,argv,NC_GEOTOP_ARCHIVE_OPTION,NC_GEOTOP_DEFINE,GEOT_VERBOSE);
		adt->ncid=ncid;

#endif	
		/*------------------    3.  Acquisition of input data and initialisation    --------------------*/
		get_all_input(argc, argv, adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I, iomanager);
		
		/*-----------------   4. Time-loop for the balances of water-mass and egy   -----------------*/
#ifdef USE_NETCDF
		set_output_nc(adt);
#endif
		time_loop(adt, iomanager);
#ifdef USE_NETCDF

		ncgt_close_geotop_archive(ncid);
		deallocate_output_nc(adt->outnc);
#endif	
		end = clock();
		elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;

		printf("Duration of simulation: %f\n", elapsed);
	

		/*--------------------   5.Completion of the output files and deallocaions  --------------------*/
		dealloc_all(adt->T, adt->S, adt->L, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->M, adt->I);
		free(adt);

	}
		
	printf("End of simulation!\n");
	
	f = fopen(SuccessfulRunFile.c_str(), "w");
	fclose(f);
	
	return 0;
}


/*----------------   6. The most important subroutine of the main: "time_loop"   ---------------*/
void time_loop(AllData *A, mio::IOManager& iomanager){

	clock_t tstart, tend;
	short en=0, wt=0, out;
//	long i;
	double t, Dt, JD0, JDb, JDe, W;
	double Vout, Voutsub, Voutsup, Vbottom;
	FILE *f;
	
	//double mean;
	
//	STATEVAR_3D *S, *G;
	Statevar3D *S, *G;
//	SOIL_STATE *L, *C;
	SoilState *L, *C;
//	STATE_VEG *V;
	StateVeg *V;
//	DOUBLEVECTOR *a, *Vsup_ch, *Vsub_ch;
	GeoVector<double> a, Vsup_ch, Vsub_ch ;
	

//	S=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
	S=new Statevar3D();
	allocate_and_initialize_statevar_3D(S, (double)number_novalue, A->P->max_snow_layers, Nr, Nc);
	if(A->P->max_glac_layers>0){
	//	G=(STATEVAR_3D *)malloc(sizeof(STATEVAR_3D));
		G=new Statevar3D();
		allocate_and_initialize_statevar_3D(G, (double)number_novalue, A->P->max_glac_layers, Nr, Nc);
	}
//	L=(SOIL_STATE *)malloc(sizeof(SOIL_STATE));
	L= new SoilState();
	initialize_soil_state(L, A->P->total_pixel, Nl);
//	C=(SOIL_STATE *)malloc(sizeof(SOIL_STATE));
	C= new SoilState();
//	initialize_soil_state(C, A->C->r->nh, Nl);
	initialize_soil_state(C, A->C->r.size(), Nl);
//	V=(STATE_VEG *)malloc(sizeof(STATE_VEG));
	V=new StateVeg();
	initialize_veg_state(V, A->P->total_pixel);

//	a=new_doublevector(A->P->total_pixel);
	a.resize(A->P->total_pixel+1);
//	Vsub_ch=new_doublevector(A->C->r->nh);
	Vsub_ch.resize(A->C->r.size());
//	Vsup_ch=new_doublevector(A->C->r->nh);
	Vsup_ch.resize(A->C->r.size());
	
	time( &start_time );

#ifdef USE_NETCDF
	long nsim=1;
	// TODO: check number of simulation issue
#else
	//long nsim = A->P->init_date->nh;
	long nsim = A->P->init_date.size();
#endif

#ifdef USE_NETCDF
		long nrun=1;
		// TODO: check number run issue
#endif
	//periods
	i_sim = i_sim0;
#ifdef USE_NETCDF
		// printing time counter initialization
		//TODO please revisit into details the initialization procedure for counters according to time coordinates
		A->counter_surface_energy=0;
		A->counter_snow=0;
		A->counter_soil=0;
		A->counter_glac=0;
		A->counter_point=0;
		if(A->P->point_sim!=1) { //distributed simulation
			A->point_var_type=NC_GEOTOP_2D_MAP_IN_CONTROL_POINT;
			A->z_point_var_type=NC_GEOTOP_3D_MAP_IN_CONTROL_POINT;
			A->unstruct_point_var_type=NC_GEOTOP_UNSTRUCT_MAP_IN_CONTROL_POINT;
			A->unstruct_z_point_var_type=NC_GEOTOP_Z_UNSTRUCT_MAP_IN_CONTROL_POINT;
		} else {// point simulation
			A->z_point_var_type=NC_GEOTOP_Z_POINT_VAR;
			A->point_var_type=NC_GEOTOP_POINT_VAR;
		}
#endif
	do {
	
		//runs
		i_run = i_run0;//Run index
		A->I->time = 0.0;//Initialize time	
		
		do {
						
			//find time step from file or inpts
			set_time_step(A->P, A->I);
				
			//time at the beginning of the time step
		//	JD0 = A->P->init_date->co[i_sim]+A->I->time/GTConst::secinday;
			JD0 = A->P->init_date[i_sim]+A->I->time/GTConst::secinday;
			
			//time step variables
			t = 0.;
			Dt = A->P->Dt;
			
			//time step subdivisions
			do{

			//	JDb = A->P->init_date->co[i_sim]+(A->I->time+t)/GTConst::secinday;
				JDb = A->P->init_date[i_sim]+(A->I->time+t)/GTConst::secinday;

				if (t + Dt > A->P->Dt) Dt = A->P->Dt - t;

				//iterations
				do{
					
				//	JDe = A->P->init_date->co[i_sim]+(A->I->time+t+Dt)/GTConst::secinday;
					JDe = A->P->init_date[i_sim]+(A->I->time+t+Dt)/GTConst::secinday;

				//	copy state variables on
				//	copy_snowvar3D(A->N->S, S);
					S=A->N->S;
				//	copy_doublevector(A->N->age, a);
					a = A->N->age;

					if (A->P->max_glac_layers>0)
				//	copy_snowvar3D(A->G->G, G);
					G=A->G->G;
				//	copy_soil_state(A->S->SS, L);
					L=A->S->SS;
				//	copy_soil_state(A->C->SS, C);
					C=A->C->SS;
				//	copy_veg_state(A->S->VS, V);
					V=A->S->VS;
					
				//	init
				//	initialize_doublevector(Vsub_ch, 0.);
					Vsub_ch.resize(Vsub_ch.size(),0.);
				//	initialize_doublevector(Vsup_ch, 0.);
					Vsup_ch.resize(Vsup_ch.size(),0.);
					Vout = 0.;
					Voutsub = 0.;
					Voutsup = 0.;
					Vbottom = 0.;
					
				//	meteo
					tstart=clock();


					mio::Date d1;
					d1.setMatlabDate(JDb, TZ); // GEOtop use matlab offset of julian date

                    std::vector<mio::MeteoData> vec_meteo;
					if(A->M->st->Z.size()>2)
                        iomanager.getMeteoData(d1, vec_meteo);

#ifndef USE_INTERNAL_METEODISTR
					A->P->use_meteoio_cloud = iswr_present(vec_meteo, (JDb == A->P->init_date[i_sim]), A);
					//	printf("met->Tgrid->ndh=%ld, met->Tgrid->nrh=%ld, met->Tgrid->nch=%ld", met->Tgrid->ndh, met->Tgrid->nrh, met->Tgrid->nch);stop_execution();

					if(A->P->point_sim == 1){
						printf("Calling pointwise MeteoIO ");
						meteoio_interpolate_pointwise( A->P, JDb, A->M, A->W);
						//	f = fopen("mio_hnw_1d_log.txt", "a");
					}else{
						printf("Calling 2D grid MeteoIO ");
						meteoio_interpolate(A->P, JDb, A->M, A->W);
						//	f = fopen("mio_hnw_2d_log.txt", "a");
						}
#else
					meteo_distr(A->M->line_interp_WEB, A->M->line_interp_WEB_LR, A->M, A->W, A->T, A->P, JD0, JDb, JDe);
					long ii,jj;
                    // for(ii=1; ii<=Nr; ii++){for(jj=1; jj<=Nc; jj++){printf("Tgrid[%ld][%ld]=%f ",ii,jj,A->M->Tgrid[jj][jj]);}printf("\n");};
					// printf("\n");for(ii=1; ii<=Nr; ii++){for(jj=1; jj<=Nc; jj++){printf("RHgrid[%ld][%ld]=%f ",ii,jj,A->M->RHgrid[jj][jj]);}printf("\n");};
					// printf("\n");for(ii=1; ii<=Nr; ii++){for(jj=1; jj<=Nc; jj++){printf("Pgrid[%ld][%ld]=%f ",ii,jj,A->M->Pgrid[jj][jj]);}printf("\n");};
					// printf("\n");for(ii=1; ii<=Nr; ii++){for(jj=1; jj<=Nc; jj++){printf("Vgrid[%ld][%ld]=%f ",ii,jj,A->M->Vgrid[jj][jj]);}printf("\n");};

#endif
					tend=clock();
					t_meteo+=(tend-tstart)/(double)CLOCKS_PER_SEC;

					if(A->P->en_balance == 1){
						tstart=clock();

						en = EnergyBalance(Dt, JD0, JDb, JDe, L, C, S, G, V, a, A, &W, vec_meteo);

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
						
						Dt *= 0.5;
						out = 0;
						
						f = fopen(logfile.c_str(), "a");
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
					
				}while( out == 0 && Dt > A->P->min_Dt ); 
				
				if (en != 0 || wt != 0) {
					f = fopen(FailedRunFile.c_str(), "w");
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
				}
						
				t += Dt;
							
			//	write state variables
			//	copy_snowvar3D(S, A->N->S);
				A->N->S=S;

			//	copy_doublevector(a, A->N->age);
				A->N->age = a;

				if (A->P->max_glac_layers>0)
				//	copy_snowvar3D(G, A->G->G);
					A->G->G=G;
			//	copy_soil_state(L, A->S->SS);
				A->S->SS=L;
			//	copy_soil_state(C, A->C->SS);
				A->C->SS=C;
			//	copy_veg_state(V, A->S->VS);
				A->S->VS=V;

			//	add_doublevector(Vsub_ch, A->C->Vsub);
			// 	TODO:Noori- we use std::transform  instead of  add_doublevector() method, need to check if they behave the same
				std::transform(Vsub_ch.data.begin(), Vsub_ch.data.end(), A->C->Vsub.data.begin(), A->C->Vsub.data.begin(),plus<double>());
			//	add_doublevector(Vsup_ch, A->C->Vsup);
			// 	TODO:Noori- we use std::transform  instead of  add_doublevector() method, need to check if they behave the same
				std::transform(Vsup_ch.data.begin(), Vsup_ch.data.end(), A->C->Vsup.data.begin(), A->C->Vsup.data.begin(),plus<double>());
				A->C->Vout += Vout;
				A->W->Voutbottom += Vbottom;
				A->W->Voutlandsub += Voutsub;
				A->W->Voutlandsup += Voutsup;

			//	record time step
			//	odb[ootimestep] = Dt * (Dt/A->P->Dtplot_basin->co[i_sim]);
				odb[ootimestep] = Dt * (Dt/A->P->Dtplot_basin[i_sim]);
				
				//write output variables
				fill_output_vectors(Dt, W, A->E, A->N, A->G, A->W, A->M, A->P, A->I, A->T);
				
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

#ifdef USE_NETCDF
			if(A->ncid==NC_GEOTOP_MISSING) {// there is no GEOtop netCDF archive, then use ascii modality
				write_output(A->I, A->W, A->C, A->P, A->T, A->L, A->S, A->E, A->N, A->G, A->M);
			//	if(strcmp(files[fSCA] , string_novalue) != 0) find_SCA(A->N->S, A->P, A->L->LC->co, A->I->time+A->P->Dt);
				if(strcmp(files[fSCA] , string_novalue) != 0) find_SCA(A->N->S, A->P, A->L->LC, A->I->time+A->P->Dt);
			}
			else {
				write_output_nc(A);
			}
#else
			write_output(A->I, A->W, A->C, A->P, A->T, A->L, A->S, A->E, A->N, A->G, A->M);
		//	if(strcmp(files[fSCA] , string_novalue) != 0) find_SCA(A->N->S, A->P, A->L->LC->co, A->I->time+A->P->Dt);
			if(strcmp(files[fSCA] , string_novalue) != 0) find_SCA(A->N->S, A->P, A->L->LC, A->I->time+A->P->Dt);
#endif
			tend=clock();
			t_out+=(tend-tstart)/(double)CLOCKS_PER_SEC;
						
			A->I->time += A->P->Dt;//Increase TIME	
			
		//	if( A->I->time > (A->P->end_date->co[i_sim] - A->P->init_date->co[i_sim])*86400. - 1.E-5){
			if( A->I->time > (A->P->end_date[i_sim] - A->P->init_date[i_sim])*86400. - 1.E-5){
				printf("Number of times the simulation #%ld has been run: %ld\n",i_sim,i_run);
				f=fopen(logfile.c_str(), "a");
				fprintf(f,"Number of times the simulation #%ld has been run: %ld\n",i_sim,i_run);
				fclose(f);
				
				i_run++;
				A->I->time = 0.0;//Initialize time
				
			//	A->P->init_date->co[i_sim0] -= A->P->delay_day_recover;
				A->P->init_date[i_sim0] -= A->P->delay_day_recover;
				A->P->delay_day_recover = 0.0;

//				A->M->line_interp_WEB_LR = 0;
//				A->M->line_interp_Bsnow_LR = 0;
//				for (i=1; i<=A->M->st->Z->nh; i++) {
//					A->M->line_interp_WEB[i-1] = 0;
//					A->M->line_interp_Bsnow[i-1] = 0;
//				}

				reset_to_zero(A->P, A->S, A->L, A->N, A->G, A->E, A->M, A->W);

			}
						
	//	}while(i_run <= A->P->run_times->co[i_sim]);//end of time-cycle
		}while(i_run <= A->P->run_times[i_sim]);//end of time-cycle
		
		reset_to_zero(A->P, A->S, A->L, A->N, A->G, A->E, A->M, A->W);
		
		i_sim++;
		i_run0 = 1;
		
  //}while (i_sim <= A->P->init_date->nh);
//	}while (i_sim < A->P->init_date.size());
	}while (i_sim < nsim);

	
//	deallocate_statevar_3D(S);
	if(A->P->max_glac_layers>0) deallocate_statevar_3D(G);
//	deallocate_soil_state(L);
//	deallocate_soil_state(C);
//	deallocate_veg_state(V);
//	free_doublevector(a);
//	free_doublevector(Vsub_ch);
//	free_doublevector(Vsup_ch);

}

/*--------------------------------------------------------------------------------------------*/

