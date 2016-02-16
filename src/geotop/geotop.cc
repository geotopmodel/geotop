/* STATEMENT:

   GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
   GEOtop 2.0.0 - 9 Mar 2012

   Copyright (c), 2012 - Stefano Endrizzi 

   This file is part of GEOtop 2.0.0 

   GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

   GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
   If you just use the code, please give feedback to the authors and the community.
   Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.

   If you have satisfactorily used the code, please acknowledge the authors.

*/


/*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/

#include "../config.h"

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
#define __MATHOPTIM_H__
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
#include <unistd.h>
#include <errno.h>

#include "inputKeywords.h"
#include "output_new.h"

#include "global_logger.h"

using namespace std;

void time_loop(AllData *A, mio::IOManager& iomanager);


/*----------   1. Global variables  ------------*/
// defined and declared into geotop.cc

/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/

int main(int argc,char *argv[]){

	clock_t start, end;
	double elapsed;
	start = clock();

	AllData *adt;
	FILE *f;

    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();

	string cfgfile = "io_it.ini";

	/* ANTONIO: This is just crazy. We have *two* incompatible way of
	   calling geotop, and part of the command line parsing is done
	   later on, in ncgt_open_from_option_string(), but only if
	   USE_NETCDF is defined, so if its not defined some of the
	   command line options are just silently ignored !

	   We should use getopt() to parse command line arguments, store
	   ignored command line options and pass them to
	   ncgt_open_from_option_string().
	   */

	std::string lDataPath ;
	if (argc >= 2)
	{
		lDataPath = argv[1] ;
	} else {
		lDataPath = get_workingdirectory() ;
	}

	if(lDataPath == "" )
	{
		std::cerr << "Error: data path is empty" << std::endl ;
		exit (200) ;
	}

	chdir(lDataPath.c_str());
	char lCWD[8192] ;
	char * lCwdStr = getcwd(lCWD, sizeof(lCWD));
	if (lCWD == NULL)
	{
		std::cerr << "Error: unable to get the current path: " << strerror(errno) << std::endl ;
		exit (201) ;
	}
	std::string lFullPath(lCwdStr) ;
	geotop::common::Variables::WORKING_DIRECTORY = lFullPath ;
	cfgfile = geotop::common::Variables::WORKING_DIRECTORY + "/" + cfgfile;

	mio::Config cfg(cfgfile);
	cfg.addKey("GRID2DPATH", "Input", "");
	mio::IOManager iomanager(cfg);
#ifdef VERBOSE 
	//MeteoIO post-initialization output
	std::cout << cfg.toString() << endl;
#endif

	geotop::common::Variables::cum_time = 0.;
	geotop::common::Variables::elapsed_time_start = 0.;

	/*dinamic allocations:*/

	geotop::common::Variables::UV=new TInit();
	if(!geotop::common::Variables::UV) t_error("UV was not allocated");

	adt=new AllData();

	if(!adt){
		t_error("adt was not allocated");
	}else {


		adt->I=new Times();
		if(!(adt->I)) t_error("times was not allocated");	

		adt->T=new Topo();
		if(!(adt->T)) t_error("top was not allocated");

		adt->S=new Soil();
		if(!(adt->S)) t_error("sl was not allocated");

		adt->L=new Land();
		if(!(adt->L)) t_error("land was not allocated");

		adt->W=new Water();
		if(!(adt->W)) t_error("water was not allocated");

		adt->P=new Par();
		if(!(adt->P)) t_error("par was not allocated");

		adt->C=new Channel();
		if(!(adt->C)) t_error("channel was not allocated");

		adt->E = new Energy();
		if(!(adt->E)) t_error("egy was not allocated");

		adt->N=new Snow();
		if(!(adt->N)) t_error("snow was not allocated");	

		adt->G = new Glacier();
		if(!(adt->G)) t_error("glac was not allocated"); 

		adt->M=new Meteo();
		if(!(adt->M)) t_error("met was not allocated"); 

		geotop::common::Variables::t_meteo=0.;
		geotop::common::Variables::t_energy=0.;
		geotop::common::Variables::t_water=0.;
		geotop::common::Variables::t_sub=0.;
		geotop::common::Variables::t_sup=0.;
		geotop::common::Variables::t_out=0.;
		geotop::common::Variables::t_blowingsnow=0.;

#ifdef USE_NETCDF
		int ncid=ncgt_open_from_option_string(argc,argv,NC_GEOTOP_ARCHIVE_OPTION,NC_GEOTOP_DEFINE,GEOT_VERBOSE);
		adt->ncid=ncid;

#endif	
		/*------------------    3.  Acquisition of input data and initialisation    --------------------*/

		get_all_input(argc, argv, adt->T, adt->S, adt->L, adt->M, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->I, iomanager);
#ifdef METEOIO_OUTPUT
        output_file_preproc(adt, cfg);
#else
        output_file_preproc(adt);
#endif

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

		lg->logsf(geotop::logger::NOTICE, "Duration of simulation: %f\n", elapsed);


		/*--------------------   5.Completion of the output files and deallocaions  --------------------*/

        deallocate_output_new();
		dealloc_all(adt->T, adt->S, adt->L, adt->W, adt->C, adt->P, adt->E, adt->N, adt->G, adt->M, adt->I);
		free(adt);

	}

	lg->log("End of simulation!\n");

	f = fopen(geotop::common::Variables::SuccessfulRunFile.c_str(), "w");
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

    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();

	//double mean;

	Statevar3D *S, *G = NULL; //TODO: check how EnergyBalance handles G 19/6/2014 Gianfranco
	SoilState *L, *C;
	StateVeg *V;

	GeoVector<double> a, Vsup_ch, Vsub_ch ;



	S=new Statevar3D();
	allocate_and_initialize_statevar_3D(S, geotop::input::gDoubleNoValue, A->P->max_snow_layers, geotop::common::Variables::Nr, geotop::common::Variables::Nc);
	if(A->P->max_glac_layers>0){
		G=new Statevar3D();
		allocate_and_initialize_statevar_3D(G, geotop::input::gDoubleNoValue, A->P->max_glac_layers, geotop::common::Variables::Nr, geotop::common::Variables::Nc);
	}

	L= new SoilState();
	initialize_soil_state(L, A->P->total_pixel, geotop::common::Variables::Nl);

	// this should be for channel (SC 08.12.2013)

	C= new SoilState();
	initialize_soil_state(C, A->C->r.size(), geotop::common::Variables::Nl);


	V=new StateVeg();
	initialize_veg_state(V, A->P->total_pixel);

	// to carefully check sizes.. 	 (SC 08.12.2013)

	a.resize(A->P->total_pixel+1);
	Vsub_ch.resize(A->C->r.size());
	Vsup_ch.resize(A->C->r.size());

	time( &geotop::common::Variables::start_time );

	long i_steps=0;

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

//////////////////////////////////////////////////////////////////////////////////
//Still to DO:
//   -consider times in seconds and not in days to avoid round-off error.. 
//   -remove all the division by GTConst::secinday
//   -better organization of timing 
//   -introduce logger.. 
	A->I->time = 0.0;//Initialize time	
	// main loop over time-steps.. 
	do {

		//time at the beginning of the time step
		JD0 = A->P->init_date + A->I->time/GTConst::secinday;

		//time step variables
		t = 0.;
		Dt = A->P->Dt;

		//time step subdivisions
		do{
			JDb = A->P->init_date + (A->I->time+t)/GTConst::secinday;
			if (t + Dt > A->P->Dt) Dt = A->P->Dt - t;
			//iterations
			do{
				JDe = A->P->init_date+(A->I->time+t+Dt)/GTConst::secinday;
#ifdef VERBOSE
				lg->logsf(geotop::logger::DEBUG,
                          "TIME-LOOP:Dtdefault:%f, JD0:%f,JDb:%f,JDe:%f\n",
                          A->P->Dt,JD0,JDb,JDe); 
#endif 					
				//	copy state variables on
				S=A->N->S;
				a = A->N->age;

				if (A->P->max_glac_layers>0) G=A->G->G;
				L=A->S->SS;
				C=A->C->SS;
				V=A->S->VS;

				//	init
				Vsub_ch.resize(Vsub_ch.size(),0.);
				Vsup_ch.resize(Vsup_ch.size(),0.);
				Vout = 0.;
				Voutsub = 0.;
				Voutsup = 0.;
				Vbottom = 0.;

				//	meteo

				tstart=clock();

				mio::Date d1;
				d1.setMatlabDate(JDe, geotop::common::Variables::TZ); // GEOtop uses matlab offset of julian date

				std::vector<mio::MeteoData> vec_meteo; //Defined here because it is passed to EnergyBalance even when USE_INTERNAL_METEODISTR is defined 13/06/2014 GG

#ifndef USE_INTERNAL_METEODISTR
				iomanager.getMeteoData(d1, vec_meteo);

				A->P->use_meteoio_cloud = iswr_present(vec_meteo, (JDb == A->P->init_date), A);

				if(A->P->point_sim == 1){
					lg->log("Calling pointwise MeteoIO", geotop::logger::DEBUG);
					meteoio_interpolate_pointwise( A->P, JDe, A->M, A->W);
				}else{
					lg->log("Calling 2D grid MeteoIO", geotop::logger::DEBUG);
					meteoio_interpolate(A->P, JDe, A->M, A->W);
					}
#else
				meteo_distr(A->M->line_interp_WEB, A->M->line_interp_WEB_LR, A->M, A->W, A->T, A->P, JD0, JDb, JDe);
#endif
				tend=clock();
				geotop::common::Variables::t_meteo+=(tend-tstart)/(double)CLOCKS_PER_SEC;

				if(A->P->en_balance == 1){
					tstart=clock();

					en = EnergyBalance(Dt, JD0, JDb, JDe, L, C, S, G, V, a, A, &W, vec_meteo);

					tend=clock();
					geotop::common::Variables::t_energy+=(tend-tstart)/(double)CLOCKS_PER_SEC;
					}

				if(A->P->wat_balance == 1 && en == 0){
					tstart=clock();
					wt = water_balance(Dt, JD0, JDb, JDe, L, C, A, Vsub_ch, Vsup_ch, &Vout, &Voutsub, &Voutsup, &Vbottom);
					tend=clock();
					geotop::common::Variables::t_water+=(tend-tstart)/(double)CLOCKS_PER_SEC;
					}

				if (en != 0 || wt != 0) {

					if(Dt > A->P->min_Dt) Dt *= 0.5;
					out = 0;

					if (en != 0) {
						lg->log("Energy balance not converging", geotop::logger::WARNING);
					}else {
						lg->log("Water balance not converging", geotop::logger::WARNING);
					}
					lg->logsf(geotop::logger::WARNING,
                              "Reducing time step to %f s, t:%f s\n",Dt,t);

				}else {
					out = 1;
				}

			}while( out == 0 && Dt > A->P->min_Dt );

            //TODO: sort out this mess about _FAILED_RUN files (G.G.)
			if (en != 0 || wt != 0) {
				f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
				fprintf(f, "Simulation Period:%ld\n",geotop::common::Variables::i_sim);
				fprintf(f, "Run Time:%ld\n",geotop::common::Variables::i_run);
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
			A->N->S=S;

			A->N->age = a;

			if (A->P->max_glac_layers>0) A->G->G=G;
			A->S->SS=L;
			A->C->SS=C;
			A->S->VS=V;

			std::transform(Vsub_ch.data.begin(), Vsub_ch.data.end(), A->C->Vsub.data.begin(), Vsub_ch.data.begin(),plus<double>());
			std::transform(Vsup_ch.data.begin(), Vsup_ch.data.end(), A->C->Vsup.data.begin(), Vsup_ch.data.begin(),plus<double>());
			A->C->Vout += Vout;
			A->W->Voutbottom += Vbottom;
			A->W->Voutlandsub += Voutsub;
			A->W->Voutlandsup += Voutsup;

			//	record time step : To further check: do we need the line below ? 
			geotop::common::Variables::odb[ootimestep] = Dt * (Dt/A->P->Dtplot_basin[geotop::common::Variables::i_sim]);

			//write output variables
			fill_output_vectors(Dt, W, A->E, A->N, A->G, A->W, A->M, A->P, A->I, A->T, A->S);

			//reset Dt
			if (Dt < A->P->Dt) Dt *= 2.;

		}while(t < A->P->Dt); // convergence loop

        //TODO: check this
		if(A->P->blowing_snow==1){
			tstart=clock();
			windtrans_snow(A->N, A->M, A->L, A->T, A->P, A->I->time);
			tend=clock();
			geotop::common::Variables::t_blowingsnow+=(tend-tstart)/(double)CLOCKS_PER_SEC;
			}

		tstart=clock();

#ifdef USE_NETCDF
		if(A->ncid==NC_GEOTOP_MISSING) {// there is no GEOtop netCDF archive, then use ascii modality
			write_output(A->I, A->W, A->C, A->P, A->T, A->L, A->S, A->E, A->N, A->G, A->M);
			if(strcmp(geotop::common::Variables::files[fSCA] , geotop::input::gStringNoValue) != 0) find_SCA(A->N->S, A->P, A->L->LC, A->I->time+A->P->Dt);
			}
		else {
			write_output_nc(A);
			}
#else
		write_output(A->I, A->W, A->C, A->P, A->T, A->L, A->S, A->E, A->N, A->G, A->M);
		// line below to move in output..   
		if(geotop::common::Variables::files[fSCA] != geotop::input::gStringNoValue) find_SCA(A->N->S, A->P, A->L->LC, A->I->time+A->P->Dt);
#endif
		tend=clock();
		geotop::common::Variables::t_out+=(tend-tstart)/(double)CLOCKS_PER_SEC;

			
		A->I->time += A->P->Dt; //Increase TIME
        write_output_new(A);

		// counter...// 
		i_steps++;
		lg->logsf(geotop::logger::TRACE,
                  "time-loop: time:%f steps:%ld enddate:%f initdate:%f diff:%f\n",
                  A->I->time,i_steps,A->P->end_date,A->P->init_date,(A->P->end_date - A->P->init_date)*86400.); 
			

	}while( (long)((A->P->end_date - A->P->init_date)*86400.) > (long)A->I->time );

///////////////////////////////////////////////////// end of time-loop... 


	if(A->P->max_glac_layers>0) deallocate_statevar_3D(G);

}

/*--------------------------------------------------------------------------------------------*/

