
/* STATEMENT:
 
 
 */

    
/*#include "struct.geotop.h"
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

extern long i_sim, i_run;*/
#ifdef USE_NETCDF

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
#include "../libraries/ascii/tabs.h"
#include "deallocate.h"

#include "../gt_utilities/gt_utilities.h"
#include "../gt_utilities/gt_symbols.h"
#include "../gt_utilities/ncgt_output.h"

#include <time.h>

extern long i_sim,i_run;
extern long number_novalue;
//***************************************************************************************************************
//***************************************************************************************************************
void set_output_nc(ALLDATA *all){
/*
 * !
 */
	//all->E->Ts_mean->name="mean_surface_temperature_in_DtSurfPrint"; // to be deleted just to test
	//all->S->Tzavplot->name="averaged_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
	//all->S->Tzplot->name="instaneous_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
	// soil
	all->S->th->name="water_content_in_soil_depth";
	all->S->T->name="temperature_in_soil_depth";
	all->S->thice->name="ice_content_in_soil_depth";
	all->S->P->name="liquid_water_pressure_in_soil_depth";
	all->S->Ptot->name="total_liquid_and_ice_water_pressure_in_soil_depth";
	// snow
	all->N->S->type->name="type_of_snow_discretization";
	all->N->S->lnum->name ="number_of_snow_layers";
	all->N->S->Dzl->name ="snow_layer_depth";
	all->N->S->T->name ="snow_layer_temperature";
	all->N->S->w_ice->name ="ice_content_in_snow_layer";
	all->N->S->w_liq->name= "water_content_in_snow_layer";
	// glacier
	/*all->G->G->type->name="type_of_glacier_discretization";
	all->G->G->lnum->name ="number_of_glacier_layers";
	all->G->G->Dzl->name ="glacier_layer_depth";
	all->G->G->T->name ="glacier_layer_temperature";
	all->G->G->w_ice->name ="ice_content_in_glacier_layer";
	all->G->G->w_liq->name= "water_content_in_glacier_layer";
	*/
	// soil
	all->outnc->soil_P_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->S->P, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->soil_T_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->S->T, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->soil_thw_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->S->th, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->soil_thi_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->S->thice, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->soil_Ptot_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->S->Ptot, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	// snow
	all->outnc->snowD_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->N->S->Dzl, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->snowT_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->N->S->T, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->snowI_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->N->S->w_ice, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->snowW_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->N->S->w_liq, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	//glacier
	/*all->outnc->glacD_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->G->G->Dzl, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->glacT_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->G->G->T, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->glacI_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->G->G->w_ice, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	all->outnc->glacW_cum=(DOUBLETENSOR *)ncgt_new_output_var((void *)all->G->G->w_liq, NC_GEOTOP_3D_MAP, (double) number_novalue, "_cum_from_previous_timestep");
	*/
	// vegetation??? TO DO
	//all->S->Tv=new_doublematrix(Nr,Nc);
}

//***************************************************************************************************************
//***************************************************************************************************************
void deallocate_output_nc(OUTPUT_NCDATA* outnc){
/*
 *
 */
	// soil
	free_doubletensor(outnc->soil_P_cum);
	free_doubletensor(outnc->soil_T_cum);
	free_doubletensor(outnc->soil_thw_cum);
	free_doubletensor(outnc->soil_thi_cum);
	free_doubletensor(outnc->soil_Ptot_cum);
	// snow
	free_doubletensor(outnc->snowD_cum);
	free_doubletensor(outnc->snowT_cum);
	free_doubletensor(outnc->snowI_cum);
	free_doubletensor(outnc->snowW_cum);
	//glacier
	free_doubletensor(outnc->glacD_cum);
	free_doubletensor(outnc->glacT_cum);
	free_doubletensor(outnc->glacI_cum);
	free_doubletensor(outnc->glacW_cum);
	free(outnc);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
void write_output_nc(ALLDATA* all){
/*
 * !
 * ! ciao
 */

	//DATA POINT
	//****************************************************************************************************************
	//****************************************************************************************************************
	if(all->P->Dtplot_point->co[i_sim] > 1.E-5){
		//if(par->Dtplot_point->co[i] > 1.E-5) par->state_pixel = 1;
	}
	/* function to write in netCDF modality */
	//all->E->Ts_mean->name="mean_surface_temperature_in_DtSurfPrint"; // to be deleted just to test
	//all->S->T->name="instantaneous_temperature_in_soil_depth"; // to be deleted just to test
	//all->S->Tzavplot->name="averaged_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
	//all->S->Tzplot->name="instaneous_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
	if(all->P->output_surfenergy>0 && fmod(all->I->time+all->P->Dt,all->P->output_surfenergy*3600.0)<1.E-5){
#ifdef USE_NETCDF_ONGOING
		printf("\n sto stampando all->counter_surface_energy=%ld",all->counter_surface_energy);
		//2D map
		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->E->Ts_mean, all->I->time+all->P->Dt, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_GENERIC,
				NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
				NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE);
		//3D map
		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->S->T, all->I->time+all->P->Dt, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_GENERIC,
				NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
				NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE);
		//point Z variable
		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->S->Tzavplot, all->I->time+all->P->Dt, NC_GEOTOP_Z_POINT_VAR, NC_GEOTOP_TIME_GENERIC,
				NC_GEOTOP_Z_GENERIC,NC_GEOTOP_POINT_DIM_GENERIC,NC_GEOTOP_MISSING_DIMENSION, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
				NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE);
		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->S->Tzplot, all->I->time+all->P->Dt, NC_GEOTOP_Z_POINT_VAR, NC_GEOTOP_TIME_GENERIC,
									NC_GEOTOP_Z_GENERIC,NC_GEOTOP_POINT_DIM_GENERIC,NC_GEOTOP_MISSING_DIMENSION, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
									NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE);
#endif
	}
}

#endif
