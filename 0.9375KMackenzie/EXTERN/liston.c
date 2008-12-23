
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie 

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Matteo Dall'Amico

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
    
    
#include "constant.h"
#include "struct.geotop.09375.h"
#include "liston.h"
#include "fortranlibrary.h"
#include "times.h"
#include "rw_maps.h"

extern T_INIT *UV;
	
extern void micromet_code_(int *nx,int *ny, double *xmn, double *ymn, float *deltax, float *deltay, int *iyear_init, int *imonth_init, int *iday_init, float *xhour_init, 
				float *dt, float *undef, int *ifill, int *iobsint, float *dn, int *iter, float *curve_len_scale, float *slopewt,float *curvewt, float *topo, float *curvature, 
				float *terrain_slope, float *slope_az, float *topoflag, float *snowdepth, float *tair_grid, float *rh_grid, float *uwind_grid, float *vwind_grid, float *Qsi_grid, float *prec_grid, 
				int *i_tair_flag, int *i_rh_flag, int *i_wind_flag, int *i_solar_flag, int *i_prec_flag, float *windspd_grid, float *winddir_grid, float *windspd_flag, 
				float *winddir_flag, float *sprec, float *windspd_min, float *Qli_grid, int *i_longwave_flag, float *vegtype, float *forest_LAI, int *iyear, int *imonth, 
				int *iday, float *xhour, int *lapse_rate_user_flag, int *iprecip_lapse_rate_user_flag,float *use_shortwave_obs, float *use_longwave_obs, float *use_sfc_pressure_obs, 
				float *sfc_pressure, float *calc_subcanopy_met, float *vegsnowdepth, float *gap_frac, float *cloud_frac_factor, float *barnes_lg_domain, int *k_stn, float *xlat_grid, 
				int *nstat, double *xstn_orig, double *ystn_orig, float *elev_orig, float *Tair_orig, float *rh_orig, float *windspd_orig, float *winddir_orig, float *prec_orig, 
				int *nveg_geotop, float *LAIwinter, float *LAIsummer);
				
/*extern void snowtran_code_(float *bc_flag, float *bs_flag, float *C_z, float *conc_salt, float *deltax, float *deltay, float *dh_salt, float *dh_salt_u, float *dh_salt_v, 
				float *dh_susp, float *dh_susp_u, float *dh_susp_v, float *dt, float *dz_susp, float *fall_vel, float *fetch, float *gravity, float *h_const, float *h_star, 
				float *ht_rhobs, float *ht_windobs, int *iter, int *nx, int *ny, float *Qsalt, float *Qsalt_max, float *Qsalt_maxu, float *Qsalt_maxv, float *Qsalt_u, 
				float *Qsalt_v, float *Qsubl, float *Qsusp, float *Qsusp_u, float *Qsusp_v, float *rh_grid, float *ro_air, float *ro_snow, float *ro_water, float *swe_depth,
				float *snow_d_init, float *snow_z0, float *soft_snow_d, float *sprec, float *subgrid_flag, float *wbal_salt, float *wbal_susp, float *wbal_qsubl, float *sum_sprec,
				float *tair_grid, float *topo, float *topo_land, float *topoflag, float *tp_scale, float *Up_const, float *Ur_const, float *Utau, float *Utau_t, 
				float *uwind_grid, float *veg_z0, float *vegsnowdepth, float *vegtype, float *vonKarman, float *vwind_grid, float *wind_min, float *winddir_flag, 
				float *winddir_grid, float *windspd_flag, float *windspd_grid, float *xmu, float *z_0, float *ztop_susp, float *erosion_dist, float *wbal_subgrid, 
				float *tabler_dir, float *slope_adjust, float *Utau_t_const, float *Utau_t_flag, float *ro_soft_snow_old, float *ro_soft_snow, float *ro_nsnow); */ 
	
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
void initialize_liston(DOUBLEMATRIX *landtype, SHORTMATRIX *LU, DOUBLEMATRIX *Z, PAR *par, METEO_STATIONS *st, DOUBLEMATRIX *snowD_init, double rhosnow0, LISTON *liston){ 

	double JD0;
	long day0, month0, year0, hour0, min0, r, c, nr=Z->nrh, nc=Z->nch;
	int i;
	DOUBLEMATRIX *M;
	
	//initial time and date
	date_time(0.0, par->year0, par->JD0, 0.0, &JD0, &day0, &month0, &year0, &hour0, &min0);
	
//1.0 liston general par
	liston->nx=(int)nc;
	liston->ny=(int)nr;
	liston->dt=(float)(par->Dt);	
	liston->deltax=(float)UV->U->co[2];
	liston->deltay=(float)UV->U->co[1];
	liston->undef=(float)(UV->V->co[2]);
	liston->nvegtypes=(int)par->n_landuses;

//2.0 Matrices input 
	M=new_doublematrix(nr,nc);
	
	//a. topography
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			M->co[r][c]=Z->co[r][c];
		}
	}
	extend_topography(M, UV->V->co[2]);
	liston->topo_land=matrix2liston(M);

	//b. initial snow depth
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(Z->co[r][c]!=UV->V->co[2]){
				snowD_init->co[r][c]*=0.001;	//from [mm] to [m]
			}
		}
	}
	extend_topography(snowD_init, UV->V->co[2]);		
	liston->snow_d_init=matrix2liston(snowD_init);
		
	//c. topography taking into account snow depth
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(Z->co[r][c]!=UV->V->co[2]){
				M->co[r][c]=Z->co[r][c];
			}else{
				M->co[r][c]=UV->V->co[2];
			}
		}
	}
	extend_topography(M, UV->V->co[2]);	
	liston->topo=matrix2liston(M);
	
	//d. latitude
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			M->co[r][c]=par->latitude*180/Pi;	
		}
	}
	liston->xlat_grid=matrix2liston(M);	
		
	//e. vegetation 
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(Z->co[r][c]!=UV->V->co[2]){
				M->co[r][c]=(double)LU->co[r][c];
			}else{
				M->co[r][c]=UV->V->co[2];
			}
		}
	}
	extend_topography(M, UV->V->co[2]);		
	liston->vegtype=matrix2liston(M);
		
	//f.snow depth and snow water equivalent matrices
	liston->snowdepth=allocatem_f(nr,nc);
	
	//g. vegetation type
	//vegetation snow depth set at a minimum value (0.01) for every type of vegetation (not definitive) 
	liston->vegsnowdepth=allocatev_f(liston->nvegtypes);
	liston->veg_z0=allocatev_f(liston->nvegtypes);
	liston->gap_frac=allocatev_f(liston->nvegtypes);

	for(i=0;i<liston->nvegtypes;i++){
		liston->vegsnowdepth[i]=(float)(0.001*landtype->co[i+1][jvholdsn]);	//[m]
		liston->veg_z0[i]=(float)(landtype->co[i+1][jz0]);					//[m]
		liston->gap_frac[i]=(float)(1.0-landtype->co[i+1][jfc]);	//1-fraction of canopy
	}

	liston->forest_LAIw=allocatev_f(liston->nvegtypes);	//MicroMet input
	liston->forest_LAIs=allocatev_f(liston->nvegtypes);	//MicroMet input
	for(i=0;i<liston->nvegtypes;i++){
		liston->forest_LAIw[i]=(float)(landtype->co[i+1][jLAIw]);
		liston->forest_LAIs[i]=(float)(landtype->co[i+1][jLAIs]);
	}
    liston->forest_LAI=allocatev_f(liston->nvegtypes);	//MicroMet output
	
	//h. cumulate snow precipitation initialization
	liston->sum_sprec=allocatem_f(nr,nc);	
	for(i=1;i<=nr*nc;i++){
		liston->sum_sprec[i-1]=0.0;
	}
									
	//end new matrices
	free_doublematrix(M);
			
//3.0 micromet par
	liston->xmn=(double)UV->U->co[4];
	liston->ymn=(double)UV->U->co[3];
	liston->iyear_init=(int)year0;
	liston->imonth_init=(int)month0;
	liston->iday_init=(int)day0;
	liston->xhour_init=(float)(hour0+min0/60.0);
	liston->i_tair_flag = 1;
	liston->i_rh_flag = 1;
	liston->i_wind_flag = 1;
	liston->i_solar_flag = 1;
	liston->i_prec_flag = 1;
	liston->i_longwave_flag = 1;
	liston->lapse_rate_user_flag=0;
	liston->iprecip_lapse_rate_user_flag=0;
	liston->use_shortwave_obs=0;
	liston->use_longwave_obs=0;
	liston->use_sfc_pressure_obs=0;
	liston->barnes_lg_domain = 0.0;
	liston->windspd_min=(float)par->Vmin;
	liston->met_stE=vector2listondouble(st->E);
	liston->met_stN=vector2listondouble(st->N);
	liston->met_stZ=vector2liston(st->Z);
	liston->nstat=(int)st->Z->nh;	
	
//4.0 matrices and tensors required by micromet
	liston->curvature=allocatem_f(nr,nc);
	liston->terrain_slope=allocatem_f(nr,nc);
	liston->slope_az=allocatem_f(nr,nc);
	liston->tair_grid=allocatem_f(nr,nc);
	liston->rh_grid=allocatem_f(nr,nc);
	liston->uwind_grid=allocatem_f(nr,nc);
	liston->vwind_grid=allocatem_f(nr,nc);
	liston->Qsi_grid=allocatem_f(nr,nc);
	liston->prec_grid=allocatem_f(nr,nc);
	liston->windspd_grid=allocatem_f(nr,nc);
	liston->winddir_grid=allocatem_f(nr,nc);
	liston->sprec=allocatem_f(nr,nc);
	liston->Qli_grid=allocatem_f(nr,nc);
	liston->sfc_pressure=allocatem_f(nr,nc);
	liston->k_stn=allocatet_int(nr,nc,liston->nstat);

//5.0 snowtrans3D
	liston->gravity=g;
	liston->ht_rhobs=(float)(st->Theight->co[1]);
	liston->ht_windobs=(float)(st->Vheight->co[1]);
	liston->ro_air=1.275;
	liston->ro_snow=rhosnow0;
	liston->ro_water=1000.0;
	liston->snow_z0=(float)(par->z0_snow);
	liston->vonKarman=ka;
			
	liston->conc_salt=allocatem_f(nr,nc);
	liston->dh_salt=allocatem_f(nr,nc);
	liston->dh_salt_u=allocatem_f(nr,nc);
	liston->dh_salt_v=allocatem_f(nr,nc);
	liston->dh_susp=allocatem_f(nr,nc);
	liston->dh_susp_u=allocatem_f(nr,nc);
	liston->dh_susp_v=allocatem_f(nr,nc);
	liston->h_star=allocatem_f(nr,nc);
	liston->Qsalt=allocatem_f(nr,nc);
	liston->Qsalt_u=allocatem_f(nr,nc);
	liston->Qsalt_v=allocatem_f(nr,nc);
	liston->Qsalt_max=allocatem_f(nr,nc);
	liston->Qsalt_maxu=allocatem_f(nr,nc);
	liston->Qsalt_maxv=allocatem_f(nr,nc);
	liston->Qsubl=allocatem_f(nr,nc);
	liston->Qsusp=allocatem_f(nr,nc);
	liston->Qsusp_u=allocatem_f(nr,nc);
	liston->Qsusp_v=allocatem_f(nr,nc);
	liston->wbal_salt=allocatem_f(nr,nc);
	liston->wbal_susp=allocatem_f(nr,nc);
	liston->wbal_qsubl=allocatem_f(nr,nc);
	liston->wbal_subgrid=allocatem_f(nr,nc);
	liston->Utau=allocatem_f(nr,nc);
	liston->Utau_t=allocatem_f(nr,nc);
	liston->z_0=allocatem_f(nr,nc);
	liston->ro_soft_snow=allocatem_f(nr,nc);		
	liston->ro_soft_snow_old=allocatem_f(nr,nc);
	liston->soft_snow_d=allocatem_f(nr,nc);
	liston->sprec_geotop=allocatem_f(nr,nc);
	liston->ro_nsnow=allocatem_f(nr,nc);
	
	//Define the initial density of the soft snow layer.
	for(r=1;r<=nr*nc;r++){
		liston->ro_soft_snow_old[r-1]=50.0;
	}
		
}

//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
void deallocate_liston(LISTON *liston){

	free(liston->topo_land);	
	free(liston->topo);
	free(liston->xlat_grid);
	free(liston->vegtype);
	free(liston->snow_d_init);
	free(liston->soft_snow_d);	
	free(liston->snowdepth);
	free(liston->vegsnowdepth);
	free(liston->veg_z0);
    free(liston->forest_LAI);
	free(liston->forest_LAIw);
	free(liston->forest_LAIs);	
	free(liston->sum_sprec);
	free(liston->gap_frac);

	free(liston->met_stE);
	free(liston->met_stN);
	free(liston->met_stZ);
	
	free(liston->curvature);
	free(liston->terrain_slope);
	free(liston->slope_az);
	free(liston->tair_grid);
	free(liston->rh_grid);
	free(liston->uwind_grid);
	free(liston->vwind_grid);
	free(liston->Qsi_grid);
	free(liston->prec_grid);
	free(liston->windspd_grid);
	free(liston->winddir_grid);
	free(liston->sprec);
	free(liston->Qli_grid);
	free(liston->sfc_pressure);
	free(liston->k_stn);
	
	free(liston->conc_salt);
	free(liston->dh_salt);
	free(liston->dh_salt_u);
	free(liston->dh_salt_v);
	free(liston->dh_susp);
	free(liston->dh_susp_u);
	free(liston->dh_susp_v);
	free(liston->h_star);
	free(liston->Qsalt);
	free(liston->Qsalt_u);
	free(liston->Qsalt_v);
	free(liston->Qsalt_max);
	free(liston->Qsalt_maxu);
	free(liston->Qsalt_maxv);
	free(liston->Qsubl);
	free(liston->Qsusp);
	free(liston->Qsusp_u);
	free(liston->Qsusp_v);
	free(liston->wbal_salt);
	free(liston->wbal_susp);
	free(liston->wbal_qsubl);
	free(liston->wbal_subgrid);
	free(liston->Utau);
	free(liston->Utau_t);
	free(liston->z_0);
	free(liston->ro_soft_snow_old);
	free(liston->ro_soft_snow);	
	free(liston->sprec_geotop);
	free(liston->ro_nsnow);
}
				
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
void call_MicroMet(double t, double dt, DOUBLEMATRIX *Z, METEO *met, ENERGY *egy, SNOW *snow, DOUBLEMATRIX *prec, DOUBLEVECTOR *LAI, LISTON *liston, PAR *par){
	
	int iter, iyear, imonth, iday;
	float xhour;	
	long l, r, c, nr, nc;
	DOUBLEMATRIX *M;
	
	iter=(int)(1+t/dt);
	nr=(long)liston->ny;
	nc=(long)liston->nx;
	M=new_doublematrix(nr,nc);
	
//	snow depth(m)
	initialize_doublematrix(M,0.0);
	if(par->point_sim==0){
		for(r=1;r<=nr;r++){
			for(c=1;c<=nc;c++){
				if(Z->co[r][c]!=UV->V->co[2]){
					for(l=1;l<=snow->lnum->co[r][c];l++){
						M->co[r][c]+=snow->Dzl->co[l][r][c]/1000.0;	//SWE in metres
					}
				}else{
					M->co[r][c]=UV->V->co[2];
				}
			}
		}
	}
	extend_topography(M, UV->V->co[2]);
	matrix2liston_2(liston->snowdepth, M);
		
	
//	run micromet
	micromet_code_(&(liston->nx), &(liston->ny), &(liston->xmn), &(liston->ymn), &(liston->deltax), &(liston->deltay), &(liston->iyear_init), &(liston->imonth_init),
		&(liston->iday_init), &(liston->xhour_init), &(liston->dt), &(liston->undef), &(liston->ifill), &(liston->iobsint), &(liston->dn), &iter, 
		&(liston->curve_len_scale), &(liston->slopewt), &(liston->curvewt), liston->topo, liston->curvature, liston->terrain_slope, liston->slope_az,
		&(liston->topoflag), liston->snowdepth, liston->tair_grid, liston->rh_grid, liston->uwind_grid, liston->vwind_grid, liston->Qsi_grid, liston->prec_grid,
		&(liston->i_tair_flag), &(liston->i_rh_flag), &(liston->i_wind_flag), &(liston->i_solar_flag), &(liston->i_prec_flag), liston->windspd_grid, 
		liston->winddir_grid, &(liston->windspd_flag), &(liston->winddir_flag), liston->sprec, &(liston->windspd_min), liston->Qli_grid, &(liston->i_longwave_flag), 
		liston->vegtype, liston->forest_LAI, &iyear, &imonth, &iday, &xhour, &(liston->lapse_rate_user_flag), &(liston->iprecip_lapse_rate_user_flag), 
		&(liston->use_shortwave_obs), &(liston->use_longwave_obs), &(liston->use_sfc_pressure_obs), liston->sfc_pressure, &(liston->calc_subcanopy_met), 
		liston->vegsnowdepth, liston->gap_frac, &(liston->cloud_frac_factor), &(liston->barnes_lg_domain), liston->k_stn, liston->xlat_grid, &(liston->nstat), 
		liston->met_stE, liston->met_stN, liston->met_stZ, met->LT, met->Lrh, met->Lws, met->Lwd, met->LP, &(liston->nvegtypes), liston->forest_LAIw, liston->forest_LAIs); 
	
//	recorder the outputs	
	if(par->point_sim!=1){
		
		if(par->micromet1==1){
			liston2matrix_2(met->Tgrid,liston->tair_grid);
			liston2matrix_2(met->Pgrid,liston->sfc_pressure);
			liston2matrix_2(met->Vgrid,liston->windspd_grid);
			liston2matrix_2(met->Vdir,liston->winddir_grid);
			liston2matrix_2(met->RHgrid,liston->rh_grid);
			liston2matrix_2(prec,liston->prec_grid);
			liston2doublevector(LAI, liston->forest_LAI);
			for(r=1;r<=nr;r++){
				for(c=1;c<=nc;c++){
					met->Tgrid->co[r][c]-=tk;			//MicroMet[K] -> GEOtop[C]
					met->Pgrid->co[r][c]/=100.0;		//MicroMet[Pa] -> GEOtop[mbar]
					met->RHgrid->co[r][c]/=100.0;		//MicroMet[%] -> GEOtop[-]
					prec->co[r][c]*=1000.0;				//MicroMet[m/h] -> GEOtop[mm/h]
				}
			}
			met->LapseRate=0.006509;	//normal lapse rate
		}
		if(par->micromet2==1) liston2matrix_2(egy->SWin,liston->Qsi_grid);
		if(par->micromet3==1) liston2matrix_2(egy->LWin,liston->Qli_grid);		

	}else{
		
		if(par->micromet1==1){
			
			liston2matrix_2(M,liston->tair_grid);
			set_to_point_simulation(M, met->Tgrid, par->r_points, par->c_points);
						
			liston2matrix_2(M,liston->sfc_pressure);
			set_to_point_simulation(M, met->Pgrid, par->r_points, par->c_points);
			
			liston2matrix_2(M,liston->windspd_grid);
			set_to_point_simulation(M, met->Vgrid, par->r_points, par->c_points);

			liston2matrix_2(M,liston->winddir_grid);
			set_to_point_simulation(M, met->Vdir, par->r_points, par->c_points);

			liston2matrix_2(M,liston->rh_grid);
			set_to_point_simulation(M, met->RHgrid, par->r_points, par->c_points);
			
			liston2matrix_2(M,liston->prec_grid);
			set_to_point_simulation(M, prec, par->r_points, par->c_points);
						
			liston2doublevector(LAI, liston->forest_LAI);
	
			for(c=1;c<=par->r_points->nh;c++){
				met->Tgrid->co[1][c]-=tk;			//MicroMet[K] -> GEOtop[C]
				met->Pgrid->co[1][c]/=100.0;			//MicroMet[Pa] -> GEOtop[mbar]
				met->RHgrid->co[1][c]/=100.0;		//MicroMet[%] -> GEOtop[-]
				prec->co[1][c]*=1000.0;				//MicroMet[m/h] -> GEOtop[mm/h]
			}
			met->LapseRate=0.006509;	//normal lapse rate
		
		}
		
		if(par->micromet2==1){
			liston2matrix_2(M,liston->Qsi_grid);
			set_to_point_simulation(M, egy->SWin, par->r_points, par->c_points);
		}
			
		if(par->micromet3==1){
			liston2matrix_2(M,liston->Qli_grid);
			set_to_point_simulation(M, egy->LWin, par->r_points, par->c_points);
			
		}
	}
	
	free_doublematrix(M);
	
}

//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
void call_SnowTrans(double t, double dt, DOUBLEMATRIX *Z, SNOW *snow, DOUBLEMATRIX *Psnow, LISTON *liston){

	int iter;
	long l,r,c,nr,nc;
	DOUBLEMATRIX *M;
	
	iter=(int)(1+t/dt);

	nr=snow->lnum->nrh;
	nc=snow->lnum->nch;
	M=new_doublematrix(nr,nc);

//	SWE(m)
	initialize_doublematrix(M,0.0);
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(Z->co[r][c]!=UV->V->co[2]){
				for(l=1;l<=snow->lnum->co[r][c];l++){
					M->co[r][c]+=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c])/1000.0;	//SWE in metres
				}
			}else{
				M->co[r][c]=UV->V->co[2];
			}
		}
	}
	extend_topography(M, UV->V->co[2]);
	matrix2liston_2(liston->snowdepth, M);

//	Psnow(m)
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(Z->co[r][c]!=UV->V->co[2]){			
				M->co[r][c]=Psnow->co[r][c]/1000.0;	//from mm (GEOtop) to m (SnowTrans3D)
			}else{
				M->co[r][c]=UV->V->co[2];
			}
		}
	}
	extend_topography(M, UV->V->co[2]);	
	matrix2liston_2(liston->sprec_geotop, M);
		
//	Density of the new snow that may fall (kg/m3)
	extend_topography(snow->rho_newsnow, UV->V->co[2]);	
	matrix2liston_2(liston->ro_nsnow, snow->rho_newsnow);

//	Soft snow depth
	fmultiplydoublematrix(M, snow->softSWE, 1.0/liston->ro_snow, UV->V->co[2]);	//M [m] - softSWE [kg/m2]
	extend_topography(M, UV->V->co[2]);
	matrix2liston_2(liston->soft_snow_d,M);
	
	/*snowtran_code_(&(liston->bc_flag), &(liston->bs_flag), &(liston->C_z), liston->conc_salt, &(liston->deltax), &(liston->deltay), liston->dh_salt, liston->dh_salt_u,
				liston->dh_salt_v, liston->dh_susp, liston->dh_susp_u, liston->dh_susp_v, &(liston->dt), &(liston->dz_susp), &(liston->fall_vel), &(liston->fetch), 
				&(liston->gravity), &(liston->h_const), liston->h_star, &(liston->ht_rhobs), &(liston->ht_windobs), &iter, &(liston->nx), &(liston->ny), liston->Qsalt, 
				liston->Qsalt_max, liston->Qsalt_maxu, liston->Qsalt_maxv, liston->Qsalt_u, liston->Qsalt_v, liston->Qsubl, liston->Qsusp, liston->Qsusp_u, liston->Qsusp_v,
				liston->rh_grid, &(liston->ro_air), &(liston->ro_snow), &(liston->ro_water), liston->snowdepth, liston->snow_d_init, &(liston->snow_z0), liston->soft_snow_d, 
				liston->sprec_geotop, &(liston->subgrid_flag), liston->wbal_salt, liston->wbal_susp, liston->wbal_qsubl, liston->sum_sprec, liston->tair_grid, 
				liston->topo, liston->topo_land, &(liston->topoflag), &(liston->tp_scale), &(liston->Up_const), &(liston->Ur_const), liston->Utau, liston->Utau_t, 
				liston->uwind_grid, liston->veg_z0, liston->vegsnowdepth, liston->vegtype, &(liston->vonKarman), liston->vwind_grid, &(liston->wind_min), 
				&(liston->winddir_flag), liston->winddir_grid, &(liston->windspd_flag), liston->windspd_grid, &(liston->xmu), liston->z_0, &(liston->ztop_susp),
				&(liston->erosion_dist), liston->wbal_subgrid, &(liston->tabler_dir), &(liston->slope_adjust), &(liston->Utau_t_const), &(liston->Utau_t_flag),
				liston->ro_soft_snow_old, liston->ro_soft_snow, liston->ro_nsnow);*/

	liston2matrix_2(snow->Wsalt,liston->wbal_salt);
	liston2matrix_2(snow->Wsusp,liston->wbal_susp);
	liston2matrix_2(snow->Wsubl,liston->wbal_qsubl);
	liston2matrix_2(snow->Wsubgrid,liston->wbal_subgrid);
	liston2matrix_2(snow->ListonSWE,liston->snowdepth);
		
	liston2matrix_2(M,liston->soft_snow_d);
	//assignnovalue(M, Z, UV->V->co[2]);
	fmultiplydoublematrix(snow->softSWE1, M, liston->ro_snow, UV->V->co[2]);	//M [m] - softSWE [kg/m2]
	
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			snow->Wsalt->co[r][c]*=1000.0;		//from m to mm (or kg/m2)
			snow->Wsusp->co[r][c]*=1000.0;
			snow->Wsubl->co[r][c]*=1000.0;
			snow->Wsubgrid->co[r][c]*=1000.0;
			snow->ListonSWE->co[r][c]*=1000.0;
		}
	}
	
//	deallocate
	free_doublematrix(M);
	
}

//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
DOUBLEMATRIX *liston2matrix(float *vector, long nr, long nc){
	
	DOUBLEMATRIX *M;
	long r,c;
	
	M=new_doublematrix(nr,nc);
	
	for(c=1;c<=nc;c++){
		for(r=1;r<=nr;r++){
			M->co[nr-r+1][c]=(double)(vector[im(c,r,nc,nr)-1]);
		}
	}
	
	return(M);
	//free_doublematrix(M);
}

//*******************************************************************************************************************************************************	
void liston2matrix_2(DOUBLEMATRIX *M, float *vector){
	
	long r,c,nr=M->nrh,nc=M->nch;

	for(c=1;c<=nc;c++){
		for(r=1;r<=nr;r++){
			M->co[nr-r+1][c]=(double)(vector[im(c,r,nc,nr)-1]);
		}
	}	
}
	

//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	

float *matrix2liston(DOUBLEMATRIX *M){
	
	float *vector;
	int r,c;
	
	vector=(float *)malloc(M->nrh*M->nch*sizeof(float));

	for(c=1;c<=M->nch;c++){
		for(r=1;r<=M->nrh;r++){
			vector[im(c,r,(int)(M->nch),(int)(M->nrh))-1]=(float)(M->co[M->nrh-r+1][c]);
		}
	}
	
	return(vector);
	//free(vector);
}

//*******************************************************************************************************************************************************	
void matrix2liston_2(float *vector, DOUBLEMATRIX *M){
	
	int r,c;

	for(c=1;c<=M->nch;c++){
		for(r=1;r<=M->nrh;r++){
			vector[im(c,r,(int)(M->nch),(int)(M->nrh))-1]=(float)(M->co[M->nrh-r+1][c]);
		}
	}

}

//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
float *vector2liston(DOUBLEVECTOR *V){

	float *vector;
	int i;
	
	vector=(float *)malloc(V->nh*sizeof(float));
	
	for(i=1;i<=V->nh;i++){
		vector[i-1]=(float)V->co[i];
	}
	
	return(vector);
	//free(vector);
}

//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
double *vector2listondouble(DOUBLEVECTOR *V){

	double *vector;
	int i;
	
	vector=(double *)malloc(V->nh*sizeof(double));
	
	for(i=1;i<=V->nh;i++){
		vector[i-1]=V->co[i];
	}
	
	return(vector);
	//free(vector);
}
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
void liston2doublevector(DOUBLEVECTOR *V, float *vector){

	int i;
	
	for(i=1;i<=V->nh;i++){
		V->co[i]=(double)vector[i-1];
	}

}

//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
void extend_topography(DOUBLEMATRIX *M, double novalue){
	
	long r,c,rr,cc;
	DOUBLEMATRIX *Q;
	
	Q=new_doublematrix(M->nrh,M->nch);
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(M->co[r][c]==novalue){
				//printf("r:%ld c:%ld rr:%ld cc:%ld\n",r,c,0,0);
				find_the_nearest(r, c, novalue, M, &rr, &cc);
				//printf("r:%ld c:%ld rr:%ld cc:%ld %f\n",r,c,rr,cc,M->co[rr][cc]);
				//stop_execution();
				Q->co[r][c]=M->co[rr][cc];
			}else{
				Q->co[r][c]=M->co[r][c];
			}
		}
	}
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			M->co[r][c]=Q->co[r][c];
		}
	}
	
	free_doublematrix(Q);
		
}
				
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************	
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
			
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************		
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
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************		
void set_to_point_simulation(DOUBLEMATRIX *Mbigger, DOUBLEMATRIX *Msmall, LONGVECTOR *r, LONGVECTOR *c){

	long i;
	
	for(i=1;i<=r->nh;i++){
		Msmall->co[1][i]=Mbigger->co[r->co[i]][c->co[i]];
	}
}
//*******************************************************************************************************************************************************
//*******************************************************************************************************************************************************		
	