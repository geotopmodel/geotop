
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
#include "keywords.h"


#include <time.h>
 
*/

#include "config.h"
#include "geotop_common.h"

#ifdef USE_NETCDF

#include "output_nc.h"

using namespace std;
NetCDFIO ncio;
//***************************************************************************************************************
//***************************************************************************************************************
void set_output_nc(AllData *all){
/*
 * !
 */
	char hold;

	// write East and North variables
	double Yll=UV->U[3];
	double Xll=UV->U[4];
	std::vector<double> Xcoord;
	std::vector<double> Ycoord;
	std::vector<double> Z_soil_coord;
	std::vector<double> Z_snow_lay_num;
	std::vector<double> ID_point;
	std::vector<double> Z0_point;
	std::vector<double> asp_point;
	std::vector<double> slp_point;
	std::vector<double> lu_point;
	std::vector<double> sy_point;
	std::vector<double> East_point;
	std::vector<double> North_point;
	std::vector<double> sky_point;
	unsigned int l,r,c;

	//cout << endl<< "soil depth";
	Z_soil_coord.resize(all->S->pa.getCh()-1,-9999);

	l=0;
	Z_soil_coord[0]= 0.5*all->S->pa(1,jdz,l+1);
	for(l=1; l<all->S->pa.getCh()-1; l++){
		Z_soil_coord[l]= Z_soil_coord[l-1]+0.5*all->S->pa(1,jdz,l)+0.5*all->S->pa(1,jdz,l+1);
		//Z_soil_coord[l]=all->S->pa(1,jdz,l+1);
		//cout << "Z["<<l<<"]=" << Z_soil_coord[l] << " ";
	}
	//cin.get(hold);

	Z_snow_lay_num.resize(all->P->max_snow_layers,-9999);
	for(l=0; l<all->P->max_snow_layers; l++){
		Z_snow_lay_num[l]=l+1;
	}

	if(all->P->rc.getRows()-1 >= 1){// either 1D simulation or 2D simulations with check points
		East_point.resize(all->P->rc.getRows()-1,-9999);
		North_point.resize(all->P->rc.getRows()-1,-9999);
		ID_point.resize(all->P->rc.getRows()-1,-9999);
		Z0_point.resize(all->P->rc.getRows()-1,-9999);
		asp_point.resize(all->P->rc.getRows()-1,-9999);
		slp_point.resize(all->P->rc.getRows()-1,-9999);
		lu_point.resize(all->P->rc.getRows()-1,-9999);
		sy_point.resize(all->P->rc.getRows()-1,-9999);
		sky_point.resize(all->P->rc.getRows()-1,-9999);
		for(l=0; l<all->P->rc.getRows()-1; l++){
			r=all->P->rc[l+1][1];
			c=all->P->rc[l+1][2];
			East_point[l]=all->P->chkpt[l+1][ptX];
			North_point[l]=all->P->chkpt[l+1][ptY];
			ID_point[l]=l+1;//cout << "ID_point["<<l<<"]=" << ID_point[l] << " ";
			Z0_point[l]=all->T->Z0[r][c];
			asp_point[l]=all->T->aspect[r][c];
			slp_point[l]=all->T->slope[r][c];;
			lu_point[l]=all->T->pixel_type[r][c];
			sy_point[l]=(short)all->L->LC[r][c];
			sky_point[l]=all->T->sky[r][c];
			}
		//cin.get(hold);
		ncio.ncgt_put_doublevector(ID_point, all->ncid, NC_GEOTOP_POINT_DIM_GENERIC, NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(East_point, all->ncid, "EastCoord_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(North_point, all->ncid, "NorthCoord_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(Z0_point, all->ncid, "Elevation_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(asp_point, all->ncid, "Aspect_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(slp_point, all->ncid, "Slope_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(lu_point, all->ncid, "LandCover_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(sy_point, all->ncid, "SoilType_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
		ncio.ncgt_put_doublevector(sky_point, all->ncid, "SkyViewFactor_in_control_point", NC_GEOTOP_POINT_DIM_GENERIC);
	}

	if (all->P->point_sim!=1){// 2D simulations
		ncio.get_coordinates_from_Array2D(all->T->Z0, Xll, Yll, all->T->Z0.getRows(), all->T->Z0.getCols(), UV->U[1], Xcoord, Ycoord);
		ncio.ncgt_put_doublevector(Xcoord, all->ncid, NC_GEOTOP_XLON, NC_GEOTOP_XLON);
		ncio.ncgt_put_doublevector(Ycoord, all->ncid, NC_GEOTOP_YLAT, NC_GEOTOP_YLAT);
	}

	ncio.ncgt_put_doublevector(Z_soil_coord, all->ncid, NC_GEOTOP_Z_SOIL, NC_GEOTOP_Z_SOIL);
	ncio.ncgt_put_doublevector(Z_snow_lay_num, all->ncid, NC_GEOTOP_Z_SNOW, NC_GEOTOP_Z_SNOW);

	//all->E->Ts_mean->name="mean_surface_temperature_in_DtSurfPrint"; // to be deleted just to test
	//all->S->Tzavplot->name="averaged_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
	//all->S->Tzplot->name="instaneous_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
	// soil

	//all->S->th.name="water_content_in_soil_depth";
	all->S->th.setMetaData("water_content_in_soil_depth", "mm", 3);
	all->S->SS->T.name="temperature_in_soil_depth";
	all->S->SS->thi.name="ice_content_in_soil_depth";
	all->S->SS->P.name="liquid_water_pressure_in_soil_depth";
	all->S->Ptot.name="total_liquid_and_ice_water_pressure_in_soil_depth";
	// snow
	all->N->S->type.name="type_of_snow_discretization";
	all->N->S->lnum.name ="number_of_snow_layers";
	all->N->S->Dzl.name ="snow_layer_depth";
	all->N->S->T.name ="snow_layer_temperature";
	all->N->S->w_ice.name ="ice_content_in_snow_layer";
	all->N->S->w_liq.name= "water_content_in_snow_layer";
// glacier: commented as in the Duron template there is no glacier to test
/*	all->G->G->type.name="type_of_glacier_discretization";
	all->G->G->lnum.name ="number_of_glacier_layers";
	all->G->G->Dzl.name ="glacier_layer_depth";
	all->G->G->T.name ="glacier_layer_temperature";
	all->G->G->w_ice.name ="ice_content_in_glacier_layer";
	all->G->G->w_liq.name= "water_content_in_glacier_layer";
	*/
	// surface energy variables
	if(all->P->output_surfenergy_bin == 1){
		all->E->Rn_mean.name="mean_net_radiation_in_print_time_step";
		all->E->LWin_mean.name="mean_LWin_radiation_in_print_time_step";
		all->E->LW_mean.name="mean_LW_radiation_in_print_time_step";
		all->E->SW_mean.name="mean_SW_radiation_in_print_time_step";
		all->E->Rswdown_mean.name="mean_Rsdown_in_print_time_step";
		all->E->Rswbeam_mean.name="mean_Rswbeam_in_print_time_step";
		all->E->nDt_shadow.name="shadowed_points_in_print_time_step";
		all->E->nDt_sun.name="sun_points_in_print_time_step";
	//	all->E->Rn_max.name="max_net_radiation_in_print_time_step";
	//	all->E->Rn_min.name="min_net_radiation_in_print_time_step";
	//	all->E->LW_max.name="max_LW_radiation_in_print_time_step";
	//	all->E->LW_min.name="min_LW_radiation_in_print_time_step";
	//	all->E->SW_max.name="max_SW_radiation_in_print_time_step";
	//	all->E->Rswdown_max.name="min_SW_radiation_in_print_time_step";
		all->E->SEB_mean.name="mean_surfaceEB_in_print_time_step";
	//	all->E->G_max.name="max_ground_heat_flux_in_print_time_step";
	//	all->E->G_min.name="min_ground_heat_flux_in_print_time_step";
		all->E->H_mean.name="mean_sensible_heat_flux_in_print_time_step";
	//	all->E->H_max.name="max_sensible_heat_flux_in_print_time_step";
	//	all->E->H_min.name="min_sensible_heat_flux_in_print_time_step";
		all->E->ET_mean.name="mean_evapo_transpiration_flux_in_print_time_step";
	//	all->E->ET_max.name="max_evapo_transpiration_flux_in_print_time_step";
	//	all->E->ET_min.name="min_evapo_transpiration_flux_in_print_time_step";
		all->E->Ts_mean.name="mean_surface_temperature_in_print_time_step";
	//	all->E->Ts_max.name="max_surface_temperature_in_print_time_step";
	//	all->E->Ts_min.name="min_surface_temperature_flux_in_print_time_step";
	}
//	all->outnc=(OUTPUT_NCDATA*) malloc(sizeof(OUTPUT_NCDATA));
	all->outnc=new OutputNCData();

	if (all->P->point_sim!=1){// distributed simulation
		// soil
		//TODO: for netCDF only 1 simulation allowed!!!
	//	all->S->SS->P has also the component [0] at the interface soil-atmosphere
	//	all->outnc->soil_P_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->SS->P, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep", all->P->output_soil[1]);
		all->outnc->soil_P_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->SS->P, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep", all->P->output_soil[1]);
	//	all->outnc->soil_T_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->SS->T, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
		all->outnc->soil_T_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->SS->T, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
	//	all->outnc->soil_thw_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->th, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
		all->outnc->soil_thw_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->th, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
	//	all->outnc->soil_thi_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->SS->thi, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
		all->outnc->soil_thi_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->SS->thi, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
	//	all->outnc->soil_Ptot_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->Ptot, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
		all->outnc->soil_Ptot_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->Ptot, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_soil[1]);
		// snow
//		all->outnc->snowD_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->Dzl, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
		all->outnc->snowD_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->Dzl, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
//		all->outnc->snowT_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->T, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
		all->outnc->snowT_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->T, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
//		all->outnc->snowI_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->w_ice, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
		all->outnc->snowI_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->w_ice, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
//		all->outnc->snowW_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->w_liq, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
		all->outnc->snowW_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->w_liq, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_snow[1]);
		//glacier
	  /*all->outnc->glacD_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->Dzl, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_glac);
		all->outnc->glacT_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->T, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_glac);
		all->outnc->glacI_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->w_ice, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_glac);
		all->outnc->glacW_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->w_liq, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->output_glac);
		*/
//i_cont
		// vegetation??? TODO
		//all->S->Tv=new_doublematrix(Nr,Nc);
	} else{// point simulation
		// soil
	//	all->outnc->soil_P_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->SS->P, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep", all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->soil_P_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->SS->P, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep", all->P->Dtplot_point[geotop::common::Variables::i_sim]);
	//	all->outnc->soil_T_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->SS->T, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->soil_T_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->SS->T, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
	//	all->outnc->soil_thw_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->th, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->soil_thw_cum=(GeoMatrix<double> *)ncio.ncgt_new_output_var((void *)&all->S->th, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
	//	all->outnc->soil_thi_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->SS->thi, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->soil_thi_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->SS->thi, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
	//	all->outnc->soil_Ptot_cum=(DOUBLEMATRIX *)ncio.ncgt_new_output_var((void *)all->S->Ptot, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->soil_Ptot_cum=(GeoMatrix<double>*)ncio.ncgt_new_output_var((void *)&all->S->Ptot, NC_GEOTOP_Z_UNSTRUCT_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		// snow
//		all->outnc->snowD_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->Dzl, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->snowD_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->Dzl, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
//		all->outnc->snowT_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->T, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->snowT_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->T, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
//		all->outnc->snowI_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->w_ice, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->snowI_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->w_ice, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
//		all->outnc->snowW_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)&all->N->S->w_liq, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		all->outnc->snowW_cum=(GeoTensor<double>*)ncio.ncgt_new_output_var((void *)&all->N->S->w_liq, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point[geotop::common::Variables::i_sim]);
		//glacier
		/*all->outnc->glacD_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->Dzl, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point->co[geotop::common::Variables::i_sim]);
		all->outnc->glacT_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->T, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point->co[geotop::common::Variables::i_sim]);
		all->outnc->glacI_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->w_ice, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point->co[geotop::common::Variables::i_sim]);
		all->outnc->glacW_cum=(DOUBLETENSOR *)ncio.ncgt_new_output_var((void *)all->G->G->w_liq, NC_GEOTOP_3D_MAP, (double) geotop::input::gDoubleNoValue, "_cum_from_previous_timestep",all->P->Dtplot_point->co[geotop::common::Variables::i_sim]);
		*/
	}
}

//***************************************************************************************************************
//***************************************************************************************************************
void deallocate_output_nc(OutputNCData* outnc){
/*
 *
 */
	// soil
	delete(outnc->soil_P_cum);
	delete(outnc->soil_T_cum);
	delete(outnc->soil_thw_cum);
	delete(outnc->soil_thi_cum);
	delete(outnc->soil_Ptot_cum);
//	free_doublematrix(outnc->soil_P_cum);
//	free_doublematrix(outnc->soil_T_cum);
//	free_doublematrix(outnc->soil_thw_cum);
//	free_doublematrix(outnc->soil_thi_cum);
//	free_doublematrix(outnc->soil_Ptot_cum);
	// snow
	delete(outnc->snowD_cum);
	delete(outnc->snowT_cum);
	delete(outnc->snowI_cum);
	delete(outnc->snowW_cum);
	//glacier
	/*free_doubletensor(outnc->glacD_cum);
	free_doubletensor(outnc->glacT_cum);
	free_doubletensor(outnc->glacI_cum);
	free_doubletensor(outnc->glacW_cum);*/
	free(outnc);
}

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//void write_output_nc_scratch(ALLDATA* all){
///*
// * !
// * ! ciao
// */
//
//	//DATA POINT
//	//****************************************************************************************************************
//	//****************************************************************************************************************
//	if(all->P->Dtplot_point->co[geotop::common::Variables::i_sim] > 1.E-5){
//		//if(par->Dtplot_point->co[i] > 1.E-5) par->state_pixel = 1;
//	}
//	/* function to write in netCDF modality */
//	all->E->Ts_mean->name="mean_surface_temperature_in_DtSurfPrint"; // to be deleted just to test
//	all->S->T->name="instantaneous_temperature_in_soil_depth"; // to be deleted just to test
//	all->S->Tzavplot->name="averaged_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
//	all->S->Tzplot->name="instaneous_temperature_in_soil_depth_in_check_points_in_Dt_output"; // to be deleted just to test
//	if(all->P->output_surfenergy>0 && fmod(all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy,all->P->output_surfenergy*3600.0)<1.E-5){
//#ifdef USE_NETCDF_ONGOING
//		printf("\n sto stampando all->counter_surface_energy=%ld",all->counter_surface_energy);
//		//2D map
//		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->E->Ts_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_GENERIC,
//				NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
//				NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//		//3D map
//		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->S->T, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_GENERIC,
//				NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
//				NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//		//point Z variable
//		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->S->Tzavplot, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy, NC_GEOTOP_Z_POINT_VAR, NC_GEOTOP_TIME_GENERIC,
//				NC_GEOTOP_Z_GENERIC,NC_GEOTOP_POINT_DIM_GENERIC,NC_GEOTOP_MISSING_DIMENSION, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
//				NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//		all->counter_surface_energy=ncgt_add_output_var(all->ncid, (void *)all->S->Tzplot, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy, NC_GEOTOP_Z_POINT_VAR, NC_GEOTOP_TIME_GENERIC,
//									NC_GEOTOP_Z_GENERIC,NC_GEOTOP_POINT_DIM_GENERIC,NC_GEOTOP_MISSING_DIMENSION, all->counter_surface_energy, NC_GEOTOP_REINITIALIZE_VARIABLE,
//									NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//#endif
//	}
//}





//***************************************************************************************************************
//***************************************************************************************************************
void write_output_nc(AllData* all){
/*
 * !
 * ! ciao
 */
	//DISCHARGE
	//****************************************************************************************************************
	//if(par->state_discharge == 1 && par->Dtplot_discharge->co[geotop::common::Variables::i_sim] > 1.E-5 && strcmp(geotop::common::Variables::files[fQ] , geotop::input::gStringNoValue) != 0){}


	//DATA POINT
	//****************************************************************************************************************
	//if(all->P->Dtplot_point->co[geotop::common::Variables::i_sim] > 1.E-5){
	//if(par->Dtplot_point->co[i] > 1.E-5) par->state_pixel = 1;
	//if(par->Dtplot_point->co[geotop::common::Variables::i_sim] > 1.E-5){
		//Print of pixel-output every times->n_pixel time step
		//if (fabs(t_point - par->Dtplot_point->co[geotop::common::Variables::i_sim])<1.E-5){
	//soil water pressure head (Psi)
//	soil

//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_P_cum, (void *)all->S->SS->P, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
//	printf("all->outnc->soil_P_cum->nrl=%ld, all->outnc->soil_P_cum->nrh=%ld, all->outnc->soil_P_cum->ncl=%ld, all->outnc->soil_P_cum->nch=%ld",all->outnc->soil_P_cum->nrl, all->outnc->soil_P_cum->nrh, all->outnc->soil_P_cum->ncl, all->outnc->soil_P_cum->nch); stop_execution();

	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_P_cum, (void *)&all->S->SS->P, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
										  NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, 
										  &all->P->rc, 
										  all->T->j_cont, 
										  all->P->total_pixel);

//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_T_cum,(void *)all->S->SS->T, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_T_cum,(void *)&all->S->SS->T, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);


//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thw_cum,(void *)all->S->th, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thw_cum,(void *)&all->S->th, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);


//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thi_cum,(void *)all->S->SS->thi, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thi_cum,(void *)&all->S->SS->thi, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);


//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_Ptot_cum,(void *)all->S->Ptot, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_Ptot_cum,(void *)&all->S->Ptot, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->unstruct_z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);

//TODO:HACK: COMMENTATO LA NEVE DA QUA
// 	snow
//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowD_cum,(void *)all->N->S->Dzl, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	char hold;



	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowD_cum,(void *)&all->N->S->Dzl, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);


//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowT_cum,(void *)all->N->S->T, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowT_cum,(void *)&all->N->S->T, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);


//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowI_cum,(void *)all->N->S->w_ice, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowI_cum,(void *)&all->N->S->w_ice, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);


//	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowW_cum,(void *)all->N->S->w_liq, all->I->time+all->P->Dt,all->P->Dt,
//		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
//		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowW_cum,(void *)&all->N->S->w_liq, all->I->time+all->P->Dt,all->P->Dt,
		all->P->Dtplot_point[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SNOW,NC_GEOTOP_POINT_DIM_GENERIC,"", all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);


	// glacier
	/*all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->glacD_cum,(void *)all->G->G->Dzl, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy,all->P->Dt,
		all->P->Dtplot_point->co[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->glacT_cum,(void *)all->G->G->T, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy,all->P->Dt,
		all->P->Dtplot_point->co[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->glacI_cum,(void *)all->G->G->w_ice, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy,all->P->Dt,
		all->P->Dtplot_point->co[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc);
	all->counter_point=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->glacW_cum,(void *)all->G->G->w_liq, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy,all->P->Dt,
		all->P->Dtplot_point->co[geotop::common::Variables::i_sim], all->z_point_var_type, NC_GEOTOP_TIME_FOR_POINT_DATA,NC_GEOTOP_Z_SOIL,NC_GEOTOP_POINT_DIM_GENERIC,NULL, all->counter_point,
		NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc);
				*/
	//TODO: HACK MATTEO
//		for(i=1;i<=par->rc->nrh;i++){
//			for (j=0; j<nopnt; j++) {
//				int i, j;
//					i=1;
//					j=1;
//					if (opnt[j] >= 0) {
//					//fprintf(f,"%f",odpnt[opnt[j]][i-1]);
//					all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid, NULL, (void *)odpnt[opnt[j]][i-1], all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600,
//							all->unstruct_point_var_type, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY, NC_GEOTOP_Z_GENERIC,NC_GEOTOP_POINT_DIM_GENERIC,NULL, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy,
//							NC_GEOTOP_UPDATE_COUNTER_TIME,NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
//
//					// IL PROBLEMA E' CHE odpnt dovrebbe essere una doublematrix e non un doppio puntatore...
//					// bisognerebbe creare:
//					//ncgt_put_vector_from_matrix_vs_time
//				}
//			}
//		}
	//BASIN DATA
	//****************************************************************************************************************
	//if(par->Dtplot_basin->co[geotop::common::Variables::i_sim] > 1.E-5 && par->state_basin == 1){
		//t_basin += par->Dt;
		//if (fabs(t_basin - par->Dtplot_basin->co[geotop::common::Variables::i_sim])<1.E-5){

	// DISTRIBUTED
	/**********************************************************************************************************/
	if (all->P->point_sim!=1){
	//soil properties
	//soil water pressure head (Psi)
//	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_P_cum,(void *)all->S->SS->P, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
//			NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_P_cum,(void *)&all->S->SS->P, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
			NC_GEOTOP_REINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

	// soil Temperature
//	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_T_cum,(void *)all->S->SS->T, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_T_cum,(void *)&all->S->SS->T, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

	// soil water content (theta_w)
//	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thw_cum,(void *)all->S->th, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thw_cum,(void *)&all->S->th, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


	//soil ice content (theta_i)
//	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thi_cum,(void *)all->S->SS->thi, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_thi_cum,(void *)&all->S->SS->thi, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


	//soil water + ice pressure head (Psi_tot)
//	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_Ptot_cum,(void *)all->S->Ptot, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_soil=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->soil_Ptot_cum,(void *)&all->S->Ptot, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_soil[1]*3600, NC_GEOTOP_Z_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SOIL,NC_GEOTOP_Z_SOIL,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_soil,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

//TODO: HACK: COMMENTATO LA NEVE DA QUI
	//snow properties
	// Snow depth
//	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowD_cum,(void *)all->N->S->Dzl, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowD_cum,(void *)&all->N->S->Dzl, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

	// snow temperature
//	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowT_cum,(void *)all->N->S->T, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowT_cum,(void *)&all->N->S->T, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


	// snow Ice content
//	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowI_cum,(void *)all->N->S->w_ice, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowI_cum,(void *)&all->N->S->w_ice, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


	// snow Water content
//	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowW_cum,(void *)all->N->S->w_liq, all->I->time+all->P->Dt,all->P->Dt,
//			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
//			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_snow=ncio.ncgt_add_output_var_cumtime(all->ncid, (void *)all->outnc->snowW_cum,(void *)&all->N->S->w_liq, all->I->time+all->P->Dt,all->P->Dt,
			all->P->output_snow[1]*3600, NC_GEOTOP_3D_MAP, NC_GEOTOP_TIME_FOR_SNOW,NC_GEOTOP_Z_SNOW,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, all->counter_snow,
			NC_GEOTOP_NOREINITIALIZE_VARIABLE,NC_GEOTOP_UPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


	//if(par->output_snow>0 && fmod(times->time+par->Dt,par->output_snow*3600.0)<1.E-5){// print condition

	//glacier properties
	//if(par->glaclayer_max>0 && par->output_glac>0 && fmod(times->time+par->Dt,par->output_glac*3600.0)<1.E-5){

	//SURFACE ENERGY BALANCE - RADIATION
	/**********************************************************************************************************/
	//if(par->output_surfenergy>0 && fmod(times->time+par->Dt,par->output_surfenergy*3600.0)<1.E-5){
//	269: long ncio.ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
	//2D surface energy map


//	long ncio.ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
//			const char* dimension_y, long counter, short reinitialize, short update, short rotate_y, double geotop::input::gDoubleNoValue, LONGMATRIX *rc){
	if(all->P->output_surfenergy_bin == 1){
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Ts_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT,NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y,NC_GEOTOP_NOVALUE,NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->Ts_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT,NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y,NC_GEOTOP_NOVALUE,NULL, all->T->j_cont, all->P->total_pixel);


//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Ts_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Ts_min, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Rn_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->Rn_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->LWin_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->LWin_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->LW_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->LW_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->SW_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->SW_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Rswdown_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->Rswdown_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Rswbeam_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->Rswbeam_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->nDt_shadow, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->nDt_shadow, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
			NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
			NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->nDt_sun, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->nDt_sun, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Rn_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Rn_min, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->LW_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->LW_min, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->SW_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->Rswdown_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->SEB_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->SEB_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->G_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->G_min, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->H_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->H_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->H_min, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);

//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->ET_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)&all->E->ET_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600, NC_GEOTOP_UNSTRUCT_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL, all->T->j_cont, all->P->total_pixel);


//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->ET_max, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid,NULL, (void *)all->E->ET_min, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy->co[1]*3600, NC_GEOTOP_2D_MAP, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY,
//					NC_GEOTOP_Z_GENERIC,NC_GEOTOP_XLON,NC_GEOTOP_YLAT, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy, NC_GEOTOP_NOUPDATE_COUNTER_TIME,
//					NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, NULL);


// in case one wants to plot in a particular point the surface energy variables
//	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid, NULL, (void *)all->E->Ts_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600,
//			all->unstruct_point_var_type, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY, NC_GEOTOP_Z_GENERIC,NC_GEOTOP_POINT_DIM_GENERIC,NULL, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy,
//			NC_GEOTOP_UPDATE_COUNTER_TIME,NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, all->P->rc, all->T->j_cont, all->P->total_pixel);
	all->counter_surface_energy=ncio.ncgt_add_output_var_cumtime(all->ncid, NULL, (void *)&all->E->Ts_mean, all->I->time+all->P->Dt,all->P->Dt,all->P->output_surfenergy[1]*3600,
			all->unstruct_point_var_type, NC_GEOTOP_TIME_FOR_SURFACE_ENERGY, NC_GEOTOP_Z_GENERIC,NC_GEOTOP_POINT_DIM_GENERIC,NULL, NC_GEOTOP_NOREINITIALIZE_VARIABLE, all->counter_surface_energy,
			NC_GEOTOP_UPDATE_COUNTER_TIME,NC_GEOTOP_NOROTATE_Y, NC_GEOTOP_NOVALUE, &all->P->rc, all->T->j_cont, all->P->total_pixel);

	}
//

	}
	//VEGETATION VARIABLES
	/**********************************************************************************************************/
	//if(par->output_vegetation>0 && fmod(times->time+par->Dt,par->output_vegetation*3600.0)<1.E-5){

	//METEO
	/**********************************************************************************************************/
	//if(par->output_meteo>0 && fmod(times->time+par->Dt,par->output_meteo*3600.0)<1.E-5){

	//SPECIAL PLOTS AT SOME DAYS
	/**********************************************************************************************************/
	//if(times->JD_plots->nh > 1 && times->iplot<=times->JD_plots->nh){

	//SAVING POINTS
	/**********************************************************************************************************/

}





#endif
