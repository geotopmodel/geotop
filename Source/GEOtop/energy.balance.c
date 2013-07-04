
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
    
    
//Author: Stefano Endrizzi
//Contents: Energy balance (and also mass balance for snow and glacier)

#include "constant.h"
#include "keywords_file.h"
#include "struct.geotop.h"
#include "energy.balance.h"
#include "meteo.h"
#include "snow.h"
#include "pedo.funct.h"
#include "vegetation.h"
#include "radiation.h"
#include "turbulence.h"
#include "util_math.h"
#include "micromet.h"
#include "times.h"
#include "tables.h"

#include "PBSM.h"
//#include "SnowTran.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;
extern DOUBLEMATRIX *outdata_point;
extern DOUBLEVECTOR *outdata_basin;
extern char *MISSING_FILE;

//extern double Tgmean;

//constants used in the Newton method to solve energy balance
#define M 1
#define min_lambda 1.E-7
#define ni 1.E-4
#define num_iter_rec_energy_balance 20
#define Csnow_at_T_greater_than_0 1.E10
#define ratio_max_storage_RAIN_over_canopy_to_LSAI 0.1
#define ratio_max_storage_SNOW_over_canopy_to_LSAI 5.0
#define Tmax 100.0
#define Tmin -70.0

#define dLW_yes 1.0

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void energy_balance(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, 
					SNOW *snow, GLACIER *glac)

{

	/*r = row index
	 c = column index
	 l = layer index
	 ns = snow layer number
	 ng = glacier layer number
	 */ 
	long r, c, j, l, ns, ng=0;
	
	/*Radiation variables
	 SWin = incoming shortwave radiation
	 SW = net shortwave at the surface
	 SWbeam = direct incoming shortwave
	 SWdiff = diffuse incoming shortwave
	 SWv_vis = net shortwave absorbed by vegetation in the visible subspectrum
	 SWv_nir = net shortwave absorbed by vegetation in the near infrared subspectrum
	 SWg_vis = net shortwave absorbed by the surface in the visible subspectrum
	 SWg_nir = net shortwave absorbed by the surface in the near infrared subspectrum
	 cosinc = cosine of the sun incidence angle
	 avis_b, avis_d, anir_b, anir_d = albedo of the surface for direct or diffuse radiation (b or d) in the visible or near infrared (vis or nir)
	 E0 = sun-earth distance correction factor
	 sa = Fmax ( a , sin(alpha) ), where a is 0.05 to 0.1
	 SWdata = flag -> 0= no radiation data measured, 1=beam and diff measured, 2=global measured
	 SWupabove_v = shortwave above vegetation directed upwards
	 */
	double SWin, SW, SWbeam, SWdiff, SWv_vis, SWv_nir, SWg_vis, SWg_nir, cosinc, avis_b, avis_d, anir_b, anir_d, E0, sa=0.;
	double SWupabove_v_vis, SWupabove_v_nir, SWup_above_v;
	short SWdata;
	
	/*LWin = incoming longwave
	 LW = net longwave at the surface
	 LWv = longwave absorbed by vegetation
	 epsa = emissivity of the atmosphere
	 eps = emissivity of the surface
	 Tsurr = temperature of the surroundings
	 epsa_min = clear sky atmosphere emissivity
	 epsa_max = overcast sky atmosphere emissivity
	 LWupabove_v = shortwave above vegetation directed upwards
	 */
	double LWin, LW, LWv, epsa, eps, Tsurr, epsa_min, epsa_max, LWupabove_v;	
	
	/*TURBULENT FLUXES: [W/m2]
	 H = sensible heat flux at the surface
	 LE = latent heat flux at the surface
	 E = evaporation at the surface
	 Hv = sensible heat flux at the canopy
	 LEv = latent heat flux at the canopy
	 Etrans = canopy transpiration
	 Hg0 = sensible heat flux at the surface in the unvegetated fraction
	 Eg0 = evaporation at the surface in the unvegetated fraction
	 Hg1 = sensible heat flux at the surface in the vegetated fraction
	 Eg1 = evaporation at the surface in the vegetated fraction
	 Hadv = sensible heat carried by lateral advection, not used
	 */
	double H, LE, E, Hv, LEv, Etrans, Hg0, Eg0, Hg1, Eg1, Hadv=0.0;
	
	/*surfEB = total surface energy balance [W/m2]
	 G = heat flux going into the soil*/
	double surfEB, G;
	
	/*snow and glacier variables:
	 snowD = snow depth [mm]
	 glacD = glacier depth [mm]
	 Prain_on_snow [mm]
	 fsnow = measure from 0 to 1 of the part of canopy covered by 0, canopy fraction depends on fsnow
	 fsnownew = value of fsnow updated at the end of the subroutine
	 fsnowcan = fraction of the water on the canopy present as snow (has implication on canopy albedo)
	 Mr_snow = snow melting rate
	 Er_snow = snow evaporation rate
	 Sr_rate = snow sublimation rate
	 Mr_glac = glacier melting rate
	 Er_glac = glacier evaporation rate
	 Sr_glac = glacier sublimation rate*/
	double snowD, glacD=0.0, Prain_on_snow, fsnow, fsnownew, fsnowcan, Mr_snow, Er_snow, Sr_snow, Mr_glac, Er_glac, Sr_glac;
	
	 /*Precipitation variables
	  Prain_over, Psnow_over: rain and snow precipitation ABOVE the canopy
	  Prain, Psnow: rain and snow precipitation BELOW the canopy
	  drip_rain, drip_snow: rain and snow precipitation not intercepted or unloaded by canopy
	  max_wcan_rain, max_wcan_snow: max storage of rain and snow on the canopy*/
	double Prain_over, Psnow_over, Prain=0.0, Psnow=0.0, drip_rain, drip_snow, max_wcan_rain = 0., max_wcan_snow = 0.;
	
	/*Cloud related variables
	 tau_cloud: cloud SW radiation transmissivity from a small number to 1
	 tau_cloud_yes: 1 available / 0 not available
	 tau_cloud_av: daily average of tau_cloud (used for longwave)
	 tau_cloud_av_yes: 1 available / 0 not available
	 tau_atm = SW radiation transmissivity of clear sky*/
	double tau_cloud, tau_cloud_av;
	short tau_cloud_yes, tau_cloud_av_yes;
	double tau_atm;
	
	/*Roughness lengths, without suffix means of the surface, veg means vegetation, _ means ratio*/
	double z0, z0veg = 0., d0, d0veg = 0., z0_z0t, hveg = 0.;
	
	/*Aerodyamic resistances, h=heat v=water vapor, c=canopy for water vapour, b=boundary layer canopy i.e. canopy for heat, 
	 uc=undercanopy=surface below the canopy*/
	double rh, rv, rc, rb, r_uc, Lobukhov;
	
	/*Vegetation characteristics
	 fc = canopy fraction
	 Locc = below canopy Monin-Obukhov length
	 utop = wind speed at the canopy height
	 Qv = saturated specific humidity at the canopy*/
	 double fc, fc0, decaycoeff, Locc, u_top, Qv;
	
	/*soil characteristics
	 Er_soil = evaporation rate from soil surface
	 Sr_soil = sublimation rate from soil surface*/
	double Er_soil, Sr_soil;
	
	/*meteo*/
	double ee, dee, Tdew, Qa, RHpoint, Vpoint, Ppoint, Tpoint, Precpoint, zmeas_T, zmeas_u;	
	 
	/*others
	 Ts, Qs temperature and specific humidity of the canopy air (the CLM uses the term "surface" for that)
	 Qg specific humidity at the soil/snow surface*/
	double Ts, Qs, Qg;

	//soil type
	long sy;
	
	//land use
	short lu;	

	//parameter interpolation
	static long line_interp;
	if(times->time==0) line_interp=0;
	
	//if(times->time==0) Tgmean=-7.;
	
	
	//SHADOW & CLOUDINESS
	sun(times->JD, &(egy->hsun), &(egy->dsun), &E0, par->latitude, par->longitude, par->ST);
	shadow_haiden(par->point_sim, top, egy->hsun, egy->dsun, land->shadow);
	
	//availability of SW data to calculate tau_cloud
	if(met->column[met->nstsrad-1][iSWb]!=-1 && met->column[met->nstsrad-1][iSWd]!=-1){
		if(met->var[met->nstsrad-1][met->column[met->nstsrad-1][iSWb]]!=NoV && met->var[met->nstsrad-1][met->column[met->nstsrad-1][iSWd]]!=NoV){
			SWdata=2;
		}else{
			SWdata=0;
		}
	}else if(met->column[met->nstsrad-1][iSW]!=-1){
		if(met->var[met->nstsrad-1][met->column[met->nstsrad-1][iSW]]!=NoV){
			SWdata=1;
		}else{
			SWdata=0;
		}
	}else{
		SWdata=0;
	}
	
	if(SWdata>0 && egy->hsun>0){
		find_tau_cloud_station(met->nstsrad, met, par, egy->hsun, E0, met->st->sky->co[met->nstsrad], Asurr, &tau_cloud, &sa);
		tau_cloud_yes = 1;
	}else{
		tau_cloud_yes = 0;
	}
	
	if( met->column[met->nstcloud-1][itauC]!=-1 ){
	
		if(times->time==0) printf("\nCloudiness transmissity of solar radiation AVAILABLE\n\n");
		tau_cloud_av=met->var[met->nstcloud-1][met->column[met->nstcloud-1][itauC]];
		if(tau_cloud_av!=NoV){
			if(tau_cloud_av_yes>1) tau_cloud_av_yes=1.;
			if(tau_cloud_av_yes<0) tau_cloud_av_yes=0.;	
			tau_cloud_av_yes = 1;
		}else{
			tau_cloud_av_yes = 0;
		}
		
	}else if( met->column[met->nstcloud-1][iC]!=-1 ){
			
		if(times->time==0) printf("\nCloudiness factor AVAILABLE\n\n");
		tau_cloud_av=met->var[met->nstcloud-1][met->column[met->nstcloud-1][iC]];
		if(tau_cloud_av!=NoV){
			if(tau_cloud_av_yes>1) tau_cloud_av_yes=1.;
			if(tau_cloud_av_yes<0) tau_cloud_av_yes=0.;				
			tau_cloud_av = 1. - 0.71*tau_cloud_av;//from fraction of sky covered by clouds to cloud transmissivity after Kimball (1928)
			tau_cloud_av_yes = 1;
		}else{
			tau_cloud_av_yes = 0;
		}
	}else{

		if(times->time==0) printf("\nCloudiness data NOT AVAILABLE\n\n");
		tau_cloud_av_yes = 0;
	
	}
	
	//COMPUTATIONS FOR EACH CELL
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			
			if(land->LC->co[r][c]!=NoV){//if the pixel is not a novalue
														
				//METEO
				//RH and V				
				if(par->micromet==1){
					if(par->point_sim==1){
						RHpoint=met->RHgrid->co[par->r_points->co[c]][par->c_points->co[c]];
						Vpoint=met->Vgrid->co[par->r_points->co[c]][par->c_points->co[c]];
					}else{
						RHpoint=met->RHgrid->co[r][c];
						Vpoint=met->Vgrid->co[r][c];
					}
				}else{
					RHpoint=met->RH;
					Vpoint=met->V;
				}
						
				if(par->point_sim==1 && par->micromet==1){
					Tpoint=met->Tgrid->co[par->r_points->co[c]][par->c_points->co[c]];
					Precpoint=wat->PrecTot->co[par->r_points->co[c]][par->c_points->co[c]];
					Ppoint=met->Pgrid->co[par->r_points->co[c]][par->c_points->co[c]];					
				}else{
					Tpoint=met->Tgrid->co[r][c];
					Precpoint=wat->PrecTot->co[r][c];
					Ppoint=met->Pgrid->co[r][c];			
				}
																								
				snow->rho_newsnow->co[r][c]=rho_newlyfallensnow(Vpoint, Tpoint, Tfreezing);
			
				//value used to calculate soil characteristics
				sy=sl->type->co[r][c];
				lu=(short)land->LC->co[r][c];
						
				//INITIALIZATION
				//Sr=sublimation rate Er=evaporation rate Mr=melting rate [mm/s]
				Sr_snow=0.0; 
				Er_snow=0.0; 
				Mr_snow=0.0; 
				Sr_glac=0.0;	
				Er_glac=0.0; 
				Mr_glac=0.0; 
				Er_soil=0.0; 
				Sr_soil=0.0;
				
				//initial condition
				if(times->time==0.0){
					if(par->recover!=1){ 
						snow_layer_combination(par->alpha_snow, r, c, snow, Tpoint, par->snowlayer_inf, par->Dmin, par->Dmax, times->time);
						if(par->glaclayer_max>0) glacier_init_t0(r, c, Tpoint, glac, snow, par, times->time);
					}
				}

				//snow
				snowD=DEPTH(r, c, snow->lnum, snow->Dzl);			
				ns=snow->lnum->co[r][c];
			
				//vegetation parameters
				if(par->vegflag->co[lu]==1) meteo_interp2(&line_interp, land->vegparv[lu], land->vegparp->co[lu], jdvegprop, 1, 0, land->vegpars[lu], times->time+0.5*par->Dt, par->JD0, par->year0, NoV);

				for(j=1;j<=jdvegprop;j++){
					if( land->vegparv[lu][land->vegparp->co[lu][j]] != NoV ){
						land->vegpar->co[j] = land->vegparv[lu][land->vegparp->co[lu][j]];
						if(j==jdroot) root(land->root_fraction->nch, land->vegpar->co[jdroot], 0.0, sl->pa->co[1][jdz], land->root_fraction->co[lu]);
						/*if(j==jdz0thresveg2){
							if(land->vegpar->co[j]>=land->vegpar->co[jdz0thresveg]){
								printf("The thresholds on snow depth regarding vegetation snow burial are not set in the right way in land cover type %d : thresveg:%f thresveg2:%f\n",
									   lu,land->vegpar->co[jdz0thresveg],land->vegpar->co[j]);
								t_error("It has to be thresveg>thresveg2");
							}
						}*/
					}else{
						land->vegpar->co[j] = land->ty->co[lu][j+jHveg-1];
					}
				}
			
				if(snowD > land->vegpar->co[jdz0thresveg]){
					fsnow=1.0;
				}else if(snowD > land->vegpar->co[jdz0thresveg2]){
					fsnow=(snowD-land->vegpar->co[jdz0thresveg2])/(land->vegpar->co[jdz0thresveg]-land->vegpar->co[jdz0thresveg2]);
				}else{
					fsnow=0.0;
				}
			
				fc = land->vegpar->co[jdcf] * pow(1.0-fsnow,land->vegpar->co[jdexpveg]);
										
				//control if LSAI is too low (in order to prevent numerical problems)
				if(land->vegpar->co[jdLSAI]<LSAIthres) fc=0.0; 
															
				//glacier
				if(par->glaclayer_max>0){
					ng=glac->lnum->co[r][c];
					glacD=DEPTH(r, c, glac->lnum, glac->Dzl);
					if(ng>0){ land->vegpar->co[jdLSAI]=0.0; fc=0.0; }
				}
			

				//METEOROLOGICAL COMPUTATIONS
			
				//temperature and wind velocity measurement heights
				zmeas_u=met->st->Vheight->co[1]-1.0E-3*snowD;
				zmeas_T=met->st->Theight->co[1]-1.0E-3*snowD;
				//zmeas_u=met->st->Vheight->co[1];
				//zmeas_T=met->st->Theight->co[1];
				if(zmeas_u<0.1) zmeas_u=0.1;
				if(zmeas_T<0.1) zmeas_T=0.1;	

				//RAIN AND SNOW PRECIPITATION [mm]
				
				//calculate dew temperature (otherwise replace Tdew with met->Tgrid) to distinguish between rain and snow
				sat_vap_pressure(&ee,&dee,Tpoint,Ppoint);	
				Qa=RHpoint*spec_humidity(ee, Ppoint);			
				ee=Qa*Ppoint/(0.622+Qa*0.378);
				sat_vap_pressure_inv(&Tdew,ee,Ppoint);
				
				//convert total precipitation to [mm]
				Precpoint*=(par->Dt/3600.0);	//from [mm/h] to [mm]
				wat->PrecTot->co[r][c]=Precpoint;
						
				//distinguish between rain and snow
				part_snow(Precpoint, &Prain_over, &Psnow_over, Tdew, par->T_rain, par->T_snow);
				
				//modify rain and snow using correction factors
				Prain_over *= par->raincorrfact;
				Psnow_over *= par->snowcorrfact;
				
				//correction of snow precipitation in case of steep slopes (contribution by Stephan Gruber)
				if (par->snow_curv > 0){
					if (top->slopes->co[r][c]*180./Pi > par->snow_smin && top->slopes->co[r][c]*180./Pi <= par->snow_smax){
						//k<-exp(-((slope-smin)^2)/curv)-exp(-(smax^2)/curv)
						Psnow_over *= ( exp(-pow(top->slopes->co[r][c]*180./Pi - par->snow_smin, 2.)/par->snow_curv) -
									   exp(-pow(par->snow_smax, 2.)/par->snow_curv) );
					}else if(top->slopes->co[r][c]*180./Pi > par->snow_smin){
						Psnow_over = 0.;
					}
				}
				
				//precipitation at the surface in the unvegetated fraction
				Prain=(1.-fc)*Prain_over;
				Psnow=(1.-fc)*Psnow_over;
			
				if(fc>0){
					canopy_rain_interception(ratio_max_storage_RAIN_over_canopy_to_LSAI, land->vegpar->co[jdLSAI], Prain_over, &max_wcan_rain, 
											 &(wat->wcan_rain->co[r][c]), &drip_rain);
					canopy_snow_interception(ratio_max_storage_SNOW_over_canopy_to_LSAI, land->vegpar->co[jdLSAI], Psnow_over, sl->Tv->co[r][c], 
											 Vpoint, par->Dt, &max_wcan_snow, &(wat->wcan_snow->co[r][c]), &drip_snow);
					//precipitation at the surface in the vegetated fraction
					Prain+=fc*drip_rain;
					Psnow+=fc*drip_snow;
					if(wat->wcan_snow->co[r][c]<0) printf("Error 1 wcansnow:%f %ld %ld\n",wat->wcan_snow->co[r][c],r,c);
					if(wat->wcan_snow->co[r][c]<0) wat->wcan_snow->co[r][c]=0.0;				
				}else{
					wat->wcan_rain->co[r][c]=0.0;
					wat->wcan_snow->co[r][c]=0.0;
				}
						
				snow->Psnow->co[r][c]=Psnow;

				//SHORTWAVE RADIATION
				//initialization of shortwave radiation absorbed by soil
				SW=0.0;
				SWup_above_v=0.0;

				//if averaged cloud transmissivity is not available it is estimated through this micromet subroutine (not very reliable)
				if(tau_cloud_av_yes==0) tau_cloud_av = 1. - 0.71*find_cloudfactor(Tpoint, RHpoint, top->Z0->co[r][c], met->LRv[2], met->LRv[3]);//Kimball(1928)

				//in case of shortwave data not available
				if(tau_cloud_yes==0) tau_cloud=tau_cloud_av;	
			
				//calculation of SWin
				tau_atm=atm_transmittance(egy->hsun, Ppoint, RHpoint, Tpoint);	
				shortwave_radiation(r, c, egy->hsun, egy->dsun, E0, land->shadow->co[r][c], top->sky->co[r][c], tau_cloud, sa, top->slopes->co[r][c], 
									top->aspect->co[r][c], tau_atm, met->var[met->nstsrad-1], met->column[met->nstsrad-1], met->st->sky->co[met->nstsrad],  Asurr, &SWbeam, &SWdiff, 
									&cosinc, egy->nDt_shadow, egy->nDt_sun);
				SWin=SWbeam+SWdiff;
			
				//albedo
				if(snowD>0){
					update_snow_age(Psnow_over, snow->T->co[ns][r][c], par->Dt, &(snow->dimens_age->co[r][c]), &(snow->nondimens_age->co[r][c]));
					avis_b=snow_albedo(land->ty->co[lu][ja_vis], snowD, par->aep, par->avo, 0.2, snow->nondimens_age->co[r][c], cosinc, (*Fzen));
					avis_d=snow_albedo(land->ty->co[lu][ja_vis], snowD, par->aep, par->avo, 0.2, snow->nondimens_age->co[r][c], cosinc, (*Zero));
					anir_b=snow_albedo(land->ty->co[lu][ja_nir], snowD, par->aep, par->airo, 0.5, snow->nondimens_age->co[r][c], cosinc, (*Fzen));
					anir_d=snow_albedo(land->ty->co[lu][ja_nir], snowD, par->aep, par->airo, 0.5, snow->nondimens_age->co[r][c], cosinc, (*Zero));
				}else{
					avis_b=land->ty->co[lu][ja_vis];
					avis_d=land->ty->co[lu][ja_vis];
					anir_b=land->ty->co[lu][ja_nir];
					anir_d=land->ty->co[lu][ja_nir];
				}
				//shortwave absorbed by soil (when vegetation is not present)
				SW += (1.0-fc)*( SWdiff*(1.0-0.5*avis_d-0.5*anir_d) + SWbeam*(1.0-0.5*avis_b-0.5*anir_b) );
				SWup_above_v += (1.0-fc)*( SWdiff*(0.5*avis_d+0.5*anir_d) + SWbeam*(0.5*avis_b+0.5*anir_b) );

				//shortwave radiation absorbed by canopy 
				if(fc>0){
					if(egy->hsun>0){
						if(wat->wcan_snow->co[r][c]<0) printf("Error wcansnow:%f %ld %ld\n",wat->wcan_snow->co[r][c],r,c);
						if(wat->wcan_snow->co[r][c]<0) wat->wcan_snow->co[r][c]=0.0;
						fsnowcan=pow(wat->wcan_snow->co[r][c]/max_wcan_snow, 2./3.);
						shortwave_vegetation(0.5*SWdiff, 0.5*SWbeam, cosinc, fsnowcan, wsn_vis, Bsnd_vis, Bsnb_vis, avis_d, avis_b, land->ty->co[lu][jvCh],
											 land->ty->co[lu][jvR_vis], land->ty->co[lu][jvT_vis], land->vegpar->co[jdLSAI], &SWv_vis, &SWg_vis, &SWupabove_v_vis);				
						shortwave_vegetation(0.5*SWdiff, 0.5*SWbeam, cosinc, fsnowcan, wsn_nir, Bsnd_nir, Bsnb_nir, anir_d, anir_b, land->ty->co[lu][jvCh],
											land->ty->co[lu][jvR_nir], land->ty->co[lu][jvT_nir], land->vegpar->co[jdLSAI], &SWv_nir, &SWg_nir, &SWupabove_v_nir);
					}else{
						SWv_vis=0.0; 
						SWg_vis=0.0; 
						SWv_nir=0.0; 
						SWg_nir=0.0;
						SWupabove_v_vis = 0.0;
						SWupabove_v_nir = 0.0;
					}
					SW+=fc*(SWg_vis+SWg_nir);
					SWup_above_v+=fc*(SWupabove_v_vis+SWupabove_v_nir);
				}else{
					SWv_vis=0.0; 
					SWv_nir=0.0;
					SWupabove_v_vis = 0.0;
					SWupabove_v_nir = 0.0;
				}
			
				//correct in case of reading data
				SW=flux(1, iSWn, met->column, met->var, 1.0, SW);
						
				//Extinction coefficient for SW in the snow layers
				rad_snow_absorption(r, c, egy->SWlayer, SW, snow);			
				set_shallow_snowpack(r, c, par->Dt, snow, egy->SWlayer->co, &Mr_snow, &ns);
							
				//LONGWAVE RADIATION
				//soil-snow emissivity
				if(snowD>10){
					eps=par->epsilon_snow;
				}else{
					eps=land->ty->co[lu][jemg];
				}			
				longwave_radiation(par->state_lwrad, ee, RHpoint, Tpoint, tau_cloud_av, &epsa, &epsa_max, &epsa_min);
				Tsurr=surface(r, c, ns, ng, snow->T, glac->T, sl->T);	//Temperature of surrounding surfaces
				LWin=epsa*SB(Tpoint);

				//if incoming LW data are available, they are used (priority)
				LWin=flux(met->nstlrad, iLWi, met->column, met->var, 1.0, LWin);
						
				//roughness lengths
				update_roughness_soil(land->ty->co[lu][jz0], 0.0, 0.0, snowD, land->ty->co[lu][jz0thressoil], par->z0_snow, &z0, &d0, &z0_z0t);
				if(fc>0) update_roughness_veg(land->vegpar->co[jdHveg], snowD, zmeas_u, zmeas_T, &z0veg, &d0veg, &hveg);							


				//SOIL AND SNOW PROPERTIES
				for(l=1;l<=Nl+ns+ng;l++){
					
					if(l<=ns){	//snow
						egy->Dlay->co[l] = 1.E-3*snow->Dzl->co[ns+1-l][r][c];
						egy->wliq->co[l] = snow->w_liq->co[ns+1-l][r][c];
						egy->wice->co[l] = snow->w_ice->co[ns+1-l][r][c];
						egy->Temp->co[l] = snow->T->co[ns+1-l][r][c];
					}else if(l<=ns+ng){   //glacier
						egy->Dlay->co[l] = 1.E-3*glac->Dzl->co[ns+ng+1-l][r][c];
						egy->wliq->co[l] = glac->w_liq->co[ns+ng+1-l][r][c];
						egy->wice->co[l] = glac->w_ice->co[ns+ng+1-l][r][c];
						egy->Temp->co[l] = glac->T->co[ns+ng+1-l][r][c];
					}else{	//soil
						egy->Dlay->co[l] = 1.E-3*sl->pa->co[sy][jdz][l-ns-ng];
						egy->wliq->co[l] = sl->th->co[l-ns-ng][r][c]*egy->Dlay->co[l]*rho_w;
						egy->wice->co[l] = sl->thice->co[l-ns-ng][r][c]*egy->Dlay->co[l]*rho_w;
						egy->Temp->co[l] = sl->T->co[l-ns-ng][r][c];
					}
					
				}
				

				//ENERGY BALANCE		
				PointEnergyBalance(r, c, egy, land, sl, par, ns, ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg, d0veg, 1.0, hveg, Vpoint, Tpoint, Qa, Ppoint, 
								   met->LRv[2], eps, fc, land->vegpar->co[jdLSAI], land->vegpar->co[jddecay0], &(wat->wcan_rain->co[r][c]), max_wcan_rain, 
								   &(wat->wcan_snow->co[r][c]), max_wcan_snow, SWin, LWin, SWv_vis+SWv_nir, &LW, &H, &E, &LWv, &Hv, &LEv, &Etrans, &Ts, &Qs, 
								   Hadv, &Hg0, &Hg1, &Eg0, &Eg1, &Qv, &Qg, &Lobukhov, &rh, &rv, &rb, &rc, &r_uc, &u_top, &decaycoeff, &Locc, &LWupabove_v);		

			
				if(wat->wcan_snow->co[r][c]<0) printf("Error 2 wcansnow:%f %ld %ld\n",wat->wcan_snow->co[r][c],r,c);					
				if(wat->wcan_snow->co[r][c]<0) wat->wcan_snow->co[r][c]=0.0;
	
				LE=E*latent(egy->Temp->co[1],Levap(egy->Temp->co[1]));
				surfEB = SW + LW - H - LE;
					
				if(ns==0){
					G=surfEB;
				}else{
					G=egy->Kth1->co[ns]*(egy->Temp->co[ns]-egy->Temp->co[ns+1])/(0.5*egy->Dlay->co[ns]+0.5*egy->Dlay->co[ns+1]);				
				}
																														
				liqWBsnow(r, c, snow, &Mr_snow, &Prain_on_snow, par, top->slopes->co[r][c], Prain, egy->wice->co, egy->wliq->co, egy->Temp->co, E*par->Dt);
				iceWBsnow(par->alpha_snow, r, c, snow, Psnow, Tpoint);
				snowD=DEPTH(r, c, snow->lnum, snow->Dzl);
						
				snow_layer_combination(par->alpha_snow, r, c, snow, Tpoint, par->snowlayer_inf, par->Dmin, par->Dmax, times->time);

				if(par->glaclayer_max>0){
					WBglacier(ns, r, c, glac, &Mr_glac, par, egy->wice->co, egy->wliq->co, egy->Temp->co, E*par->Dt);
					glac2snow(r, c, snow, glac, par->Dmin, par->Dmax);
					if(par->glaclayer_max>1)glac_layer_combination(r, c, glac, Tpoint, par->glaclayer_max, par->Dmin_glac, par->Dmax_glac, times->time);
					ng=glac->lnum->co[r][c];
				}
			
				//NET PRECIPITATION
				wat->Pnet->co[r][c] = Mr_snow;

				if( land->vegpar->co[jdLSAI]>=LSAIthres && ng==0 ){
				
					fc0=fc;	
				
					if(snowD > land->vegpar->co[jdz0thresveg]){
						fsnownew=1.0;
					}else if(snowD > land->vegpar->co[jdz0thresveg2]){
						fsnownew=(snowD-land->vegpar->co[jdz0thresveg2])/(land->vegpar->co[jdz0thresveg]-land->vegpar->co[jdz0thresveg2]);
					}else{
						fsnownew=0.0;
					}
				
					fc = land->vegpar->co[jdcf]*pow(1.0-fsnownew,land->vegpar->co[jdexpveg]);	
				
					//a) fc increases
					if(fc>fc0){
						wat->wcan_rain->co[r][c]*=(fc0/fc);
						wat->wcan_snow->co[r][c]*=(fc0/fc);					
					
						//b) fc decreases
					}else if(fc<fc0){	
						iceWBsnow(par->alpha_snow, r, c, snow, wat->wcan_snow->co[r][c]*(fc0-fc), sl->Tv->co[r][c]);
						wat->Pnet->co[r][c] +=  wat->wcan_rain->co[r][c]*(fc0-fc)/par->Dt;
					}
				
				}
			
				wat->Pnet->co[r][c] += (Mr_glac + Prain/par->Dt);	
			
			
				if(wat->wcan_snow->co[r][c]<0) printf("Error 3 wcansnow:%f %ld %ld\n",wat->wcan_snow->co[r][c],r,c);
				if(wat->wcan_snow->co[r][c]<0) wat->wcan_snow->co[r][c]=0.0;

																													
				//OUTPUT DATA
				if(ns>0){
					Sr_snow=E;
				}else if(ng>0){
					Sr_glac=E;
				}else{
					Er_soil=E;
				}

				if(par->point_sim!=1) prepare_output(Er_soil, Mr_snow, Er_snow, Sr_snow, Mr_glac, Er_glac, Sr_glac, Prain, Psnow_over, 
						egy, wat, snow, glac, land, top, sl, met, times, par, r, c, 0.25*(avis_d+anir_d+avis_b+anir_b), LE, surfEB, 
						H, G, egy->Temp->co[1], SWin, SW, SWbeam, eps, LWin, LW, cosinc, Precpoint);
						
				output_pixel(r, c, Psnow, Prain-Prain_on_snow, Prain_on_snow, Sr_soil, Er_soil, Mr_snow, Er_snow, Sr_snow, Mr_glac, 
						Er_glac, Sr_glac, Etrans, LE, H, surfEB, G, egy->Temp->co[1], 0.25*(avis_d+anir_d+avis_b+anir_b), eps, land->vegpar->co[jdLSAI], 
						LWin, SWin, LW, SW, epsa, epsa_min, epsa_max, SWbeam, SWdiff, Tdew, times->n_pixel, Lobukhov, 
						par, wat, egy, top, met, snow, glac, land, sl, Tpoint, Ppoint, Vpoint, RHpoint, Psnow_over, Prain_over, z0, 
						z0veg, d0veg, (SWv_vis+SWv_nir), LWv, sl->Tv->co[r][c], Ts, Hg0, Hg1, Eg0, Eg1, fc, rh, 
						rv, rb, rc, r_uc, Hv, LEv, Qv, Qg, Qa, Qs, u_top, decaycoeff, Locc, SWup_above_v, LWupabove_v);
										
				output_basin(times->n_basin, par->total_pixel, Prain, Psnow, Prain_over, Psnow_over, Tpoint, egy->Temp->co[1], sl->Tv->co[r][c], 
						Er_soil, Etrans, LE, H, SW, LW, fc*LEv, fc*Hv, fc*(SWv_vis+SWv_nir), fc*LWv, SWin, LWin, par->Dt);

				output_map_plots(r, c, par, times->time, times->n_plot, egy, met, snow, H, LE, fc*Hv, fc*LEv, SWin, SW, 
						fc*(SWv_vis+SWv_nir), LWin, LW, fc*LWv, Ts, egy->Temp->co[1], sl->Tv->co[r][c]);							
			
			}
		}
	}
	
	//ACCOUNT FOR SNOW WIND TRANSPORT

	//PBSM
	if(par->blowing_snow==1) set_windtrans_snow(snow, met, land, par, times->time);

	//SnowTran3D
	//if(par->blowing_snow==1) set_windtrans_snow2(snow, met, land, top, par, times->time);

	//CALCULATES SNOW COVERED AREA
	if(par->point_sim!=1 && strcmp(files->co[fSCA]+1 , MISSING_FILE) != 0) find_SCA(snow, par, land->LC->co, times->time);

	//PREPARE SNOW OUTPUT
	output_snow(snow, land->LC->co, par);
	
	/*Tgmean=0.;
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(land->LC->co[r][c]!=NoV){
				Tgmean+=egy->Temp->co[1]/(double)par->total_pixel;
			}
		}
	}
	printf("Tgmean:%f ",Tgmean);
	Tgmean=Tgmean*Fmin(Fmax(1., (Tgmean-2.)),2.);
	printf("Tgmean:%f\n",Tgmean);*/
}

//end of "energy_balance" subroutine

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/************************************************************ƒ******************************************************************************/

void PointEnergyBalance(long r, long c, ENERGY *egy, LAND *land, SOIL *sl, PAR *par, long ns, long ng, double zmu, double zmT, 
			double z0s, double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, 
			double P, double LR, double eps, double fc, double LSAI, double decaycoeff0, double *Wcrn, double Wcrnmax, double *Wcsn, 
			double Wcsnmax, double SWin, double LWin, double SWv,double *LW, double *H, double *E, double *LWv, double *Hv, double *LEv,
			double *Etrans, double *Ts, double *Qs, double Hadd, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Qv, 
			double *Qg, double *Lobukhov, double *rh, double *rv, double *rb, double *rc, double *r_uc, double *u_top, 
			double *decay, double *Locc, double *LWup_above_v){

		
	short iter_close, iter_close2, lu=(short)land->LC->co[r][c];

	long sy=sl->type->co[r][c];
	long l, m, cont=0, cont2, n=Nl+ns+ng;
	
	double dH_dT, dE_dT, EB, dEB_dT, EB0, Tg, psim;
	double res, res0[3], res_av, res_prev[M], lambda[3], C0, C1, th0, th1, kbb0, kbb1, k, kn;
	double Qg0, dQ, Tv0, dWcsn=0.0, dWcrn=0.0, rh_g, rv_g, th_oversat;
				
	//Assign vectors		
	for(l=1;l<=n;l++){ 
		egy->T0->co[l] = egy->Temp->co[l];
		egy->deltaw->co[l] = 0.0; //mass that melts (if>0) or freezes (if<0)
	}
	
	for(l=1;l<=Nl;l++){
		egy->THETA->co[l] = sl->th->co[l][r][c];	//water content to be modified at each iteration	
	}
	
	//surface skin temperature used to compute the surface energy fluxes
	Tg = egy->Temp->co[1]; 
	//Tg=Tgmean;
		
	//canopy temperature at the beginning of the time step
	Tv0=sl->Tv->co[r][c];
	
	//Qg0 is the specific humidity at the soil surface at the beginning of the time step
	sat_spec_humidity(&Qg0, &dQ, 1., Tg, P);	
	
	
	/*calculates all the surface energy fluxes, includes the calculation of vegetation temperature (solving the vegetation energy balance)
	Returns:
	 
	 LW: net longwave radiation at the surface
	 H: sensible heat flux at the surface (positive upwards)
	 E: evaporation from the soil surface (mm/s)
	 
	 LWv: net longwave in the canopy
	 Hv: sensible heat flux from the canopy
	 LEv: latent heat flux from the canopy
	 Etrans: canopy transpiration (as water flux in mm/s)
	 
	 sl->Tv->co[r][c] : vegetation temperature
	 Qv : specific humidity at the canopy (assumed saturated)
	 Ts : temperature of canopy air
	 Qs : specific humidity of canopy air*/
	
	//surface energy balance
	EnergyFluxes(Tg, r, c, ns+ng, egy->T0->co[1], Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, 
				 P, LR, sl->P->co[0][r][c], eps, fc, LSAI, decaycoeff0, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, 
				 egy->THETA->co, sl->pa->co[sy], land->ty->co[lu], land->root_fraction->co[lu], par, egy->soil_transp_layer, 
				 SWin, LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans, &(sl->Tv->co[r][c]), Qv, Ts, Qs, Hg0, Hg1, 
				 Eg0, Eg1, Lobukhov, rh, rv, rc, rb, r_uc, &rh_g, &rv_g, Qg, u_top, decay, Locc, LWup_above_v,
				 egy->Temp->co, egy->soil_evap_layer_bare, egy->soil_evap_layer_bare);
	
	
	EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
	dEB_dT = dLW_yes*(-eps*dSB_dT(Tg)) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;
					

	EB0 = EB;
	

	//soil freezing/thawing
	for(l=1;l<=Nl;l++){
		//total pressure (=pressure of the water would have if it was all in liquind state)
		psim=psi_teta(sl->th->co[l][r][c]+sl->thice->co[l][r][c], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], 
					sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
		//max temperature at which the first particle of ice comes up
		egy->Tstar->co[l]=Fmin(psim/(1000.0*Lf/(g*(Tfreezing+tk))), 0.0);
	}
	
	//F(T) = diag(egy->Fenergy(T)) + K(T)*T
	for(l=1;l<=n;l++){
		
		if(l>1) kn = k;
		k = calc_k(l, r, c, ns, ng, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Temp->co, egy->Dlay->co, (*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par);		
		if(l>1){
			egy->Kth0->co[l-1] = -k * kn / ( k * 0.5 * egy->Dlay->co[l-1] + kn * 0.5 * egy->Dlay->co[l] );	// (-k/dz)
			egy->Kth1->co[l-1] = egy->Kth0->co[l-1];
		}
		
		C0 = calc_C(l, r, c, ns+ng, 0.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
		C1 = calc_C(l, r, c, ns+ng, 1.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
		if(l<=ns+ng && fabs(egy->deltaw->co[l]-egy->wice->co[l])<1.E-7) C1 = Csnow_at_T_greater_than_0;
		
		egy->Fenergy->co[l] = ( Lf*egy->deltaw->co[l] + C1*egy->Dlay->co[l]*egy->Temp->co[l] - C0*egy->Dlay->co[l]*egy->T0->co[l] ) / par->Dt;
		
		if(l<=ns+1) egy->Fenergy->co[l] -= egy->SWlayer->co[l];
		
	}
	
	//top boundary condition (treated as sink)
	egy->Fenergy->co[1] -= ( (1.-KNe)*EB + KNe*EB0 );
		
	//bottom boundary condition (treated as sink)
	kbb0 = k;
	kbb1 = k;
	egy->Fenergy->co[n] -= ( (par->Tboundary-egy->Temp->co[n])*(1.-KNe)*kbb1 + (par->Tboundary-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlay->co[n]/2.+par->Zboundary);

	
	update_F_energy(n, egy->Fenergy, 1.-KNe, egy->Kth1, egy->Temp->co);
	update_F_energy(n, egy->Fenergy, KNe, egy->Kth0, egy->T0->co);	
	res = norm_2(egy->Fenergy, n);
	
	do{
		
		cont++;
		
		for(l=1;l<=n;l++){
			egy->T1->co[l] = egy->Temp->co[l];
		}
		
		//calculates J(T)
		for(l=1;l<=n;l++){
				
			//thermal conductivity
			C1 = calc_C(l, r, c, ns+ng, 1.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
			if(l<=ns+ng && fabs(egy->deltaw->co[l]-egy->wice->co[l])<1.E-7) C1 = Csnow_at_T_greater_than_0;
			
			//apparent thermal conductivities (phase change)		
			if(l<=ns+ng){//snow
				C1 += Lf*(egy->wice->co[l]+egy->wliq->co[l])*dtheta_snow(par->alpha_snow, egy->Temp->co[l])/egy->Dlay->co[l];
			}else{	//soil
				m=l-ns-ng;	
				if(egy->Temp->co[l]<=egy->Tstar->co[m]) C1 += rho_w*Lf*(Lf/(g*tk)*1.E3)*dteta_dpsi(Psif(egy->Temp->co[l]), 0.0, sl->pa->co[sy][jsat][m], 
											sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 
											1-1/sl->pa->co[sy][jns][m], PsiMin, 0.0);
			}
			
						
			egy->dFenergy->co[l] = C1*egy->Dlay->co[l]/par->Dt;

		}
		
		egy->dFenergy->co[1] -= (1.-KNe)*dEB_dT;
		egy->dFenergy->co[n] += (1.-KNe)*kbb1 / (egy->Dlay->co[n]/2.+par->Zboundary);
			

		update_diag_dF_energy(n, egy->dFenergy, 1.-KNe, egy->Kth1);
		tridiag2(1, r, c, n, 1.-KNe, egy->Kth1, 1., egy->dFenergy, 1.-KNe, egy->Kth1, egy->Fenergy, egy->Newton_dir);
		
		cont2=0;
		iter_close=0;
		res0[0]=res;
		
		//non-monotonic line search (it is monotonic if M==1)
		for(m=Fminlong(cont,M);m>1;m--){
			res_prev[m-1]=res_prev[m-2];
		}
		res_prev[0]=res;
		
		res_av=0.0;
		for(m=1;m<=Fminlong(cont,M);m++){
			res_av=Fmax(res_prev[m-1],res_av);
		}
		
		do{

			cont2++;
			iter_close2=0;
			
			if(cont2 == 1){
				lambda[0] = 1.0;
				
			}else if(cont2 == 2){
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = thmax;
								
			}else{
				lambda[2] = lambda[1];
				res0[2] = res0[1];
				lambda[1] = lambda[0];
				res0[1] = res;
				lambda[0] = minimize_merit_function(res0[0], lambda[1], res0[1], lambda[2], res0[2]);

			}
									
			for(l=1;l<=n;l++){ 
	
				egy->Temp->co[l] = egy->T1->co[l] + lambda[0] * egy->Newton_dir->co[l];
				if(egy->Temp->co[l] != egy->Temp->co[l]) printf("T nan l:%ld r:%ld c:%ld cont:%ld\n",l,r,c,cont); 

				//soil
				if(l > ns+ng){
					m=l-ns-ng;
										
					th0=sl->th->co[m][r][c];
					th1=Fmin( teta_psi(Psif(Fmin(egy->Tstar->co[m],egy->Temp->co[l])), 0.0, sl->pa->co[sy][jsat][m], sl->pa->co[sy][jres][m], 
									  sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 1-1/sl->pa->co[sy][jns][m], PsiMin, 
									  sl->pa->co[sy][jss][m]), 
							  sl->th->co[m][r][c]+sl->thice->co[m][r][c] );
										
					egy->deltaw->co[l]=(th1-th0)*egy->Dlay->co[l]*rho_w;
									
				//snow	
				}else{
					
					th0=egy->wliq->co[l]/(egy->wice->co[l]+egy->wliq->co[l]);
					th1=theta_snow(par->alpha_snow, egy->Temp->co[l]);
					
					egy->deltaw->co[l]=(th1-th0)*(egy->wice->co[l]+egy->wliq->co[l]);

				}
												
			}
			
			if(cont < num_iter_rec_energy_balance){
				Tg=egy->Temp->co[1];
				//Tg=Tgmean;
			
				//update egy->THETA taking into account evaporation (if there is not snow)
				if(ns+ng == 0){
					for(l=ns+ng+1;l<=n;l++){
						m=l-ns-ng;
						egy->THETA->co[m] = sl->th->co[m][r][c] + egy->deltaw->co[l]/(rho_w*egy->Dlay->co[l]);
			
						//add canopy transpiration
						if(egy->THETA->co[m] > sl->pa->co[sy][jres][1] + 1.E-3 && l <= egy->soil_transp_layer->nh ){
							egy->THETA->co[m] -= Fmax( par->Dt*fc*egy->soil_transp_layer->co[m]/(rho_w*egy->Dlay->co[l]), 0.0 );
							if(egy->THETA->co[m] < sl->pa->co[sy][jres][m]+1.E-3) egy->THETA->co[m] = sl->pa->co[sy][jres][m]+1.E-3;
						}
				
						//add soil evaporation
						if(egy->THETA->co[m] > sl->pa->co[sy][jres][1] + 1.E-3 && l <= egy->soil_evap_layer_bare->nh ){
							egy->THETA->co[m] -= Fmax( par->Dt*(1.-fc)*egy->soil_evap_layer_bare->co[m]/(rho_w*egy->Dlay->co[l]), 0.0 );
							egy->THETA->co[m] -= Fmax( par->Dt*fc*egy->soil_evap_layer_veg->co[m]/(rho_w*egy->Dlay->co[l]), 0.0 );
							if(egy->THETA->co[m] < sl->pa->co[sy][jres][m]+1.E-3) egy->THETA->co[m] = sl->pa->co[sy][jres][m]+1.E-3;
						}				
					}
				}
			
				//surface energy balance
				EnergyFluxes_no_rec_turbulence(Tg, r, c, ns+ng, egy->T0->co[1], Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, 
						 P, LR, sl->P->co[0][r][c], eps, fc, LSAI, decaycoeff0, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, egy->THETA->co,
						 sl->pa->co[sy], land->ty->co[lu], land->root_fraction->co[lu], par, egy->soil_transp_layer, SWin,
						 LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans, &(sl->Tv->co[r][c]), Qv, Ts, Qs, Hg0, Hg1, 
						 Eg0, Eg1, Lobukhov, rh, rv, rc, rb, r_uc, &rh_g, &rv_g, Qg, u_top, decay, Locc, LWup_above_v, egy->Temp->co, 
						 egy->soil_evap_layer_bare, egy->soil_evap_layer_bare);
				
				EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
				dEB_dT = dLW_yes*(-eps*dSB_dT(Tg)) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;
			}
			
			//F(T) = diag(egy->Fenergy(T)) + K(T)*T
			for(l=1;l<=n;l++){
				
				if(l>1) kn = k;
				k = calc_k(l, r, c, ns, ng, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Temp->co, egy->Dlay->co, (*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par);		
				if(l>1){
					egy->Kth1->co[l-1] = -k * kn / ( k * 0.5 * egy->Dlay->co[l-1] + kn * 0.5 * egy->Dlay->co[l] );	// (-k/dz)
				}
				
				C0 = calc_C(l, r, c, ns+ng, 0.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
				C1 = calc_C(l, r, c, ns+ng, 1.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
				if(l<=ns+ng && fabs(egy->deltaw->co[l]-egy->wice->co[l])<1.E-7) C1 = Csnow_at_T_greater_than_0;

				egy->Fenergy->co[l] = ( Lf*egy->deltaw->co[l] + C1*egy->Dlay->co[l]*egy->Temp->co[l] - C0*egy->Dlay->co[l]*egy->T0->co[l] ) / par->Dt;
				
				if(l<=ns+1) egy->Fenergy->co[l] -= egy->SWlayer->co[l];
				
			}
			
			//top boundary condition (treated as sink)
			egy->Fenergy->co[1] -= ( (1.-KNe)*EB + KNe*EB0 );
			
			//bottom boundary condition (treated as sink)
			kbb1 = k;
			egy->Fenergy->co[n] -= ( (par->Tboundary-egy->Temp->co[n])*(1.-KNe)*kbb1 + (par->Tboundary-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlay->co[n]/2.+par->Zboundary);
			
			
			update_F_energy(n, egy->Fenergy, 1.-KNe, egy->Kth1, egy->Temp->co);
			update_F_energy(n, egy->Fenergy, KNe, egy->Kth0, egy->T0->co);	
			res = norm_2(egy->Fenergy, n);
						
			if(res <= res_av*(1.0 - ni*lambda[0])) iter_close2=1;
			if(lambda[0] <= min_lambda) iter_close2=1;
																		
		}while(iter_close2!=1);
				
		if(iter_close>=0 && res<=par->tol_energy) iter_close=1;
		if(cont>=par->maxiter_energy) iter_close=1;

	}while(iter_close!=1);
	
	if(res>par->tol_energy) printf("ERROR: ENERGY BALANCE DOES NOT CONVERGE res:%e T:%f r:%ld c:%ld cont:%ld cont2:%ld\n",res,egy->Temp->co[1],r,c,cont,cont2);
			
	for(l=1;l<=Nl+ns+ng;l++){ 	

		//checks
		if(egy->wice->co[l]-egy->deltaw->co[l]<0) egy->deltaw->co[l]=egy->wice->co[l];
		if(egy->wliq->co[l]+egy->deltaw->co[l]<0) egy->deltaw->co[l]=-egy->wliq->co[l];

		//update
		egy->wliq->co[l]+=egy->deltaw->co[l];  
		egy->wice->co[l]-=egy->deltaw->co[l];
		if(egy->Temp->co[l]>0) egy->wice->co[l]=0.0;

		//error messages
		if(egy->Temp->co[l]!=egy->Temp->co[l]) printf("T no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
		if(egy->deltaw->co[l]!=egy->deltaw->co[l]) printf("dw no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
		if(egy->Temp->co[l]<Tmin || egy->Temp->co[l]>Tmax) printf("T outside of range, error 2, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f SW:%f LW:%f H:%f LE:%f Ta:%f\n",l,r,c,egy->Temp->co[l],egy->SWlayer->co[1],*LW,*H,latent(Tg,Levap(Tg))*(*E),Ta);
			
		
		//update soil
		if(l>ns+ng){
			m=l-ns-ng;
			
			sl->ET->co[m][r][c] = 0.;
			
			//canopy transpiration
			if(l <= egy->soil_transp_layer->nh) sl->ET->co[m][r][c] += fc*egy->soil_transp_layer->co[m];
			
			//soil evaporation
			if(l <= egy->soil_evap_layer_bare->nh){
				sl->ET->co[m][r][c] += (1.-fc)*egy->soil_evap_layer_bare->co[m];
				sl->ET->co[m][r][c] += fc*egy->soil_evap_layer_veg->co[m];
			}

			th_oversat = Fmax( sl->P->co[m][r][c] - psisat_from(m, r, c, sl) , 0.0 ) * sl->pa->co[sy][jss][m];
			
			sl->th->co[m][r][c] =  egy->wliq->co[l]/(rho_w*egy->Dlay->co[l]);
			sl->thice->co[m][r][c] = egy->wice->co[l]/(rho_w*egy->Dlay->co[l]);
			sl->T->co[m][r][c] = egy->Temp->co[l];

			sl->P->co[m][r][c] = psi_from_theta(sl->th->co[m][r][c] + th_oversat, m, r, c, sl, PsiMin);		

		}
		
	} 	
		
	*Wcrn = *Wcrn + dWcrn;
	*Wcsn = *Wcsn + dWcsn;
		
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void update_F_energy(long n, DOUBLEVECTOR *F, double w, DOUBLEVECTOR *K, double *T){
	
	long l;
	
	for(l=1;l<=n;l++){
		
		if(l==1){
			F->co[l] += w*(-K->co[l]*T[l] + K->co[l]*T[l+1]); 
		}else if(l<n){
			F->co[l] += w*(K->co[l-1]*T[l-1] - (K->co[l]+K->co[l-1])*T[l] + K->co[l]*T[l+1]);
		}else{
			F->co[l] += w*(K->co[l-1]*T[l-1] - K->co[l-1]*T[l]);
		}
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_diag_dF_energy(long n, DOUBLEVECTOR *dF, double w, DOUBLEVECTOR *K){
	
	long l;
	
	for(l=1;l<=n;l++){
		
		if(l==1){
			dF->co[l] -= w*K->co[l]; 
		}else if(l<n){
			dF->co[l] -= w*(K->co[l]+K->co[l-1]);
		}else{
			dF->co[l] -= w*K->co[l-1];
		}
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double calc_C(long l, long r, long c, long nsng, double a, double *wi, double *wl, double *dw, double *D, double *T, SOIL *sl){

	double C;
	long sy = sl->type->co[r][c];
	
	if(l<=nsng){	//snow
		C = (c_ice*(wi[l]-a*dw[l]) + c_liq*(wl[l]+a*dw[l]))/D[l];
	}else{	//soil
		C = sl->pa->co[sy][jct][l-nsng]*(1.-sl->pa->co[sy][jsat][l-nsng]) + (c_ice*(wi[l]-a*dw[l]) + c_liq*(wl[l]+a*dw[l]))/D[l];
	}
		
	return(C);
	
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
	
double calc_k(long l, long r, long c, long ns, long ng, double *wi, double *wl, double *dw, double *T, double *D, 
			  double (* kfunct_snow)(double rho), double (* kfunct_glac)(double rho), SOIL *sl, PAR *par){	

	double k;
	long sy = sl->type->co[r][c];
		
	if(l<=ns){
		k = (*kfunct_snow)((wl[l]+wi[l])/D[l]);
	}else if(l<=ng+ns){
		k = (*kfunct_glac)((wl[l]+wi[l])/D[l]);
	}else{
		k = k_thermal_soil((wl[l]+dw[l])/(rho_w*D[l]), (wi[l]-dw[l])/(rho_w*D[l]), sl->pa->co[sy][jsat][l-ns-ng], 
						   sl->T->co[l-ns-ng][r][c], sl->pa->co[sy][jkt][l-ns-ng]);
	}
	
	return(k);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void EnergyFluxes(double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
		double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, 
		double P, double LR, double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, 
		double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn, double *dWcsn, double *theta, double **soil, 
		double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer, double SWin, double LWin, double SWv, double *LW, 
		double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv, double *LEv, double *Etrans, 
		double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Lobukhov, 
		double *rh, double *rv, double *rc, double *rb, double *r_uc, double *rh_g, double *rv_g, double *Qg, 
		double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T, DOUBLEVECTOR *soil_evap_layer_bare,
		DOUBLEVECTOR *soil_evap_layer_veg){

	
	double Hg, dHg_dT, Eg, dEg_dT, LWg, LWup;
	double dQgdT, rm;
	double alpha, beta;	//from Ye and Pielke, 1993
		
	*H=0.0;	
	*E=0.0;	
	
	*dH_dT=0.0;	
	*dE_dT=0.0;	
	
	*LW=0.0; 
	*LWup_above_v=0.0;
	
	*LWv=0.0;
	*Hv=0.0; 
	*LEv=0.0; 
	
	*Etrans=0.0; 																																													
	*Hg0=0.0; 
	*Eg0=0.0;
	*Hg1=0.0; 
	*Eg1=0.0;
	
	*decay=0.0; 
	*Locc=0.0;
	
	//thermodynamical calculations
	sat_spec_humidity(Qg, &dQgdT, 1.0, Tg, P);
	
	//initialization
	if(fc==0){ 
		
		*Qs=*Qg;  
		*Ts=Tg; 
		*Qv=*Qg; 
		*Tv=Tg; 
		*rh=1.E99; 
		*rv=1.E99; 
		*rc=1.E99; 
		*rb=1.E99; 
		*r_uc=1.E99;
		*u_top=0.0; 
		
	}else if(fc==1){ 
		
		*rh_g=1.E99; 
		*rv_g=1.E99; 
		
	}
		
	//bare soil
	if(fc<1){
		aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, LR, Lobukhov, &rm, rh_g, rv_g, par->state_turb, par->monin_obukhov, par->maxiter_Businger);
		
		find_actual_evaporation_parameters(r,c,&alpha, &beta, soil_evap_layer_bare, theta, soil, T, psi, P, *rv_g, Ta, Qa, *Qg, n);
		turbulent_fluxes(*rh_g, *rv_g/beta, P, Ta, Tg, Qa, *Qg*alpha, dQgdT*alpha, &Hg, &dHg_dT, &Eg, &dEg_dT);
								
		*H+=(1.0-fc)*Hg;	
		*E+=(1.0-fc)*Eg;	

		*dH_dT+=(1.0-fc)*dHg_dT;		
		*dE_dT+=(1.0-fc)*dEg_dT;		
				
		*LW+=(1.0-fc)*( e*(LWin-SB(Tg)) );
		*LWup_above_v+=(1.0-fc)*( (1.0-e)*LWin+e*SB(Tg) );
								
		*Hg0=Hg;
		*Eg0=Eg;		
		
		if(Hg!=Hg){
			printf("Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
				   zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			stop_execution();
		}
		if(Eg!=Eg){
			printf("Eg no value bare soil %ld %ld \n",r,c);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
				   zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			stop_execution();
		}

	}
	
	//vegetation
	if(fc>0){
										
		Tcanopy(r, c, Tv0, Tg, *Qg, dQgdT, Tg0, Qg0, Ta, Qa, zmu, zmT, z0v, z0s, d0v, rz0v, hveg, v, LR, P, SWin, SWv, LWin, 
				e, LSAI, decaycoeff0, land, Wcrn, Wcrnmax, Wcsn, Wcsnmax, dWcrn, dWcsn, LWv, &LWg, Hv, &Hg, &dHg_dT, LEv, &Eg, 
				&dEg_dT, Ts, Qs, root, theta, soil_transp_layer, Lobukhov, par, n, rh, rv, rc, rb, r_uc, u_top, Etrans, Tv, Qv, 
				decay, Locc, &LWup, psi, soil, T, soil_evap_layer_veg);
		
		*H+=fc*Hg;	
		*E+=fc*Eg;
		
		*dH_dT+=fc*dHg_dT;
		*dE_dT+=fc*dEg_dT;		
				
		*LW+=fc*LWg;
		*LWup_above_v+=fc*LWup;
				
		*Hg1=Hg;
		*Eg1=Eg;
		
		if(Hg!=Hg){
			printf("Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
				   zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			stop_execution();
		}		
		if(Eg!=Eg){
			printf("Eg no value bare soil %ld %ld \n",r,c);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",
				   zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			stop_execution();
		}		
	}	
		
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void EnergyFluxes_no_rec_turbulence(double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
		double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, 
		double P, double LR, double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, 
		double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn, double *dWcsn, double *theta, double **soil, 
		double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer, double SWin, double LWin, double SWv, double *LW, 
		double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv, double *LEv, double *Etrans, 
		double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Lobukhov, 
		double *rh, double *rv, double *rc, double *rb, double *r_uc, double *rh_g, double *rv_g, double *Qg, 
		double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T, DOUBLEVECTOR *soil_evap_layer_bare,
		DOUBLEVECTOR *soil_evap_layer_veg){

	
	double Hg, dHg_dT, Eg, dEg_dT, LWg, LWup;
	double dQgdT;
	double dLWvdT,LWv_old;
	double alpha, beta;	//from Ye and Pielke, 1993

	//initalization		
	*H=0.0;	
	*E=0.0;	
	
	*dH_dT=0.0;	
	*dE_dT=0.0;	
	
	*LW=0.0;
	*LWup_above_v=0.0;
			
	//thermodynamical calculations
	sat_spec_humidity(Qg, &dQgdT, 1.0, Tg, P);
		
	if(fc<1){		
		find_actual_evaporation_parameters(r,c,&alpha, &beta, soil_evap_layer_bare, theta, soil, T, psi, P, *rv_g, Ta, Qa, *Qg, n);
		turbulent_fluxes(*rh_g, *rv_g/beta, P, Ta, Tg, Qa, *Qg*alpha, dQgdT*alpha, &Hg, &dHg_dT, &Eg, &dEg_dT);
		
		*H+=(1.0-fc)*Hg;	
		*E+=(1.0-fc)*Eg;	

		*dH_dT+=(1.0-fc)*dHg_dT;
		*dE_dT+=(1.0-fc)*dEg_dT;		
				
		*LW+=(1.0-fc)*( e*LWin-SB(Tg) );
		*LWup_above_v+=(1.0-fc)*( (1.0-e)*LWin+e*SB(Tg) );
						
		*Hg0=Hg;
		*Eg0=Eg;		
		
		//error messages
		if(Hg!=Hg){
			printf("Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f psi:%f \n",
				   zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,psi);
			stop_execution();
		}
		if(Eg!=Eg){
			printf("Eg no value bare soil %ld %ld \n",r,c);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f psi:%f \n",
				   zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,psi);
			stop_execution();
		}
		
	}
	
	//TURBULENT FLUXES FOR CANOPY CANOPY
	if(fc>0){

		find_actual_evaporation_parameters(r,c,&alpha, &beta, soil_evap_layer_veg, theta, soil, T, psi, P, *r_uc, Ta, Qa, *Qg, n);
		*Ts=(Ta/(*rh)+Tg/(*r_uc)+(*Tv)/(*rb))/(1./(*rh)+1./(*r_uc)+1./(*rb));
		*Qs=(Qa/(*rv)+(*Qg)*alpha*beta/(*r_uc)+(*Qv)/(*rc))/(1./(*rv)+beta/(*r_uc)+1./(*rc));
		
		turbulent_fluxes(*r_uc, *r_uc/beta, P, *Ts, Tg, *Qs, *Qg*alpha, dQgdT*alpha, &Hg, &dHg_dT, &Eg, &dEg_dT);
		
		LWv_old=*LWv;
		longwave_vegetation(LWin, e, Tg, *Tv, LSAI, LWv, &LWg, &dLWvdT, &LWup);
		*LWv=LWv_old;
		
		*H+=fc*Hg;	*dH_dT+=fc*dHg_dT;
		*E+=fc*Eg;	*dE_dT+=fc*dEg_dT;		
				
		*LW+=fc*LWg;
		*LWup_above_v+=fc*LWup;
				
		*Hg1=Hg;
		*Eg1=Eg;
		
		//error messages
		if(Hg!=Hg) printf("Hg no value canopy %ld %ld \n",r,c);
		if(Eg!=Eg) printf("Eg no value canopy %ld %ld \n",r,c);
		
	}	
		
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


double surface(long r, long c, long ns, long ng, DOUBLETENSOR *snow, DOUBLETENSOR *ice, DOUBLETENSOR *sl){
	double S;
	if(ns>0){
		S=snow->co[ns][r][c];
	}else if(ng>0){
		S=ice->co[ng][r][c];
	}else{
		S=sl->co[1][r][c];
	}
	return(S);
}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double k_thermal_soil(double th_liq, double th_ice, double th_sat, double T, double k_solid){

	/*double rho_dry, sr, ke, k_dry, k_sat, k;

	//dry sl density [kg/mc]
	rho_dry=2700.0*(1.0-th_sat);

	//saturation degree
	sr=(th_liq+th_ice)/th_sat;

	//Kersten number
	if(T>=Tfreezing){
		//ke=Fmax( 0. , log10(sr)+1.0 );
		if(th_sat>0.8){
			ke=sr;
		}else{
			ke=Fmax(0. , log10(sr)+1.0);
		}
	}else{
		//ke=sr;
		if(th_sat>0.8){
			ke=sr*sr;
		}else{
			ke=sr;
		}
	}

	//dry soil thermal conductivity [W m^-1 K^-1]
	k_dry=(0.135*rho_dry+64.7)/(2700.0-0.947*rho_dry);

	//soil thermal conductivity [W m^-1 K^-1] Farouki (1981)
	if(sr>1.0E-7){	
		//saturated sl thermal conductivity [W m^-1 K^-1]
		if(T>Tfreezing){
			k_sat=pow(k_solid,1.0-th_sat)*pow(k_liq,th_sat);
		}else{
			k_sat=pow(k_solid,1.0-th_sat)*pow(k_liq,th_liq)*pow(k_ice,th_sat-th_liq);
		}
		//soil thermal conductivity [W m^-1 K^-1]
		k=ke*k_sat + (1.0-ke)*k_dry;
	}else{
		k=k_dry;
	}
	
	return(k);*/
	
	/*Quadratic parallel based on the analogy with electric lattice 
	 (Cosenza et al., 2003, European Journal of Soil Science, September 2003, 54, 581–587)*/
	
	return ( pow ( (1.-th_sat)*sqrt(k_solid) + th_liq*sqrt(k_liq) + th_ice*sqrt(k_ice) + (th_sat-th_liq-th_ice)*sqrt(k_air) , 2. ) );
	
	/*double liq=0.5*th_sat;
	double ice=th_sat-liq;
	return ( pow ( (1.-th_sat)*sqrt(k_solid) + liq*sqrt(k_liq) + ice*sqrt(k_ice) + (th_sat-liq-ice)*sqrt(k_air) , 2. ) );*/
	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


double flux(long i, long icol, long **col, double **met, double k, double est){
	
	double F;
								
	if(col[i-1][icol]!=-1){
		F=k*met[i-1][col[i-1][icol]];
		if(F==NoV) F=est;
	}else{
		F=est;
	}
		
	return(F);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_map_plots(long r, long c, PAR *par, double t, long n, ENERGY *egy, METEO *met, SNOW *snow, double Hg, 
					  double LEg, double Hv, double LEv, double SWin, double SWg, double SWv, double LWin, double LWg, 
					  double LWv, double Ts, double Tg, double Tv){
	
	long j, d, mon, y, h, m;
	double tmin, tmax, JD, N;

	if(par->JD_plots->co[1]>=0){
	
		date_time(t, par->year0, par->JD0, 0.0, &JD, &d, &mon, &y, &h, &m);

		for(j=1;j<=(long)(par->JD_plots->nh/2.);j++){
			
			get_time( &tmin, par->JD_plots->co[2*j-1], y, par->JD0, par->year0 );
			get_time( &tmax, par->JD_plots->co[2*j]  , y, par->JD0, par->year0 );

			if( floor((tmax-tmin)/(n*par->Dt)) != (tmax-tmin)/(n*par->Dt) ){
				N=floor(((tmax-tmin)/(n*par->Dt)));
				tmax=tmin+N*n*par->Dt;
			}
					
			if(t>=tmin && t<tmax){	
				egy->Hgplot->co[r][c]+=Hg/(double)n;

				egy->LEgplot->co[r][c]+=LEg/(double)n;
				egy->Hvplot->co[r][c]+=Hv/(double)n;
				egy->LEvplot->co[r][c]+=LEv/(double)n;
								
				egy->SWinplot->co[r][c]+=SWin/(double)n;
				egy->SWgplot->co[r][c]+=SWg/(double)n;
				egy->SWvplot->co[r][c]+=SWv/(double)n;				
				
				egy->LWinplot->co[r][c]+=LWin/(double)n;
				egy->LWgplot->co[r][c]+=LWg/(double)n;	
				egy->LWvplot->co[r][c]+=LWv/(double)n;	
							
				egy->Tsplot->co[r][c]+=Ts/(double)n;
				egy->Tgplot->co[r][c]+=Tg/(double)n;
				egy->Tvplot->co[r][c]+=Tv/(double)n;
				
				snow->Dplot->co[r][c]+=DEPTH(r,c,snow->lnum,snow->Dzl)/(double)n;
				met->Taplot->co[r][c]+=met->Tgrid->co[r][c]/(double)n;
				if(par->micromet==1){
					met->Vspdplot->co[r][c]+=met->Vgrid->co[r][c]/(double)n;
					met->Vdirplot->co[r][c]+=met->Vdir->co[r][c]/(double)n;
					met->RHplot->co[r][c]+=met->RHgrid->co[r][c]/(double)n;
				}
			}
		}
	}
}
				

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e, 
				  double *T, SHORTVECTOR *mf){

	long l;
	short occ=0;
	
	for(l=1;l<=n;l++){
		if(e->co[l]!=e->co[l]) occ=1;
	}
	
	if(occ==1 ){
		printf("NOvalue in Energy Balance: r:%ld c:%ld\n",r,c);
		printf("l:%d adi:--- ad:%f ads:%f b:%f e:%f T:%f mf:%d\n",
			   1,ad->co[1],ads->co[1],b->co[1],e->co[1],T[1],mf->co[1]);
		for(l=2;l<=n-1;l++){
			printf("l:%ld adi:%f ad:%f ads:%f b:%f e:%f T:%f mf:%d\n",
				   l,adi->co[l-1],ad->co[l],ads->co[l],b->co[l],e->co[l],T[l], mf->co[l]);
		}
		printf("l:%ld adi:%f ad:%f ads:--- b:%f e:%f T:%f mf:%d\n",
			   n,adi->co[n-1],ad->co[n],b->co[n],e->co[n],T[n],mf->co[n]);
		stop_execution();
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
	
double soil_red_evap(double psi, double T){
	double alpha;
	alpha=exp(1.0E-3*psi*g/(461.50464*(T+tk)));
	return(alpha);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double red_evap(long n, double psi, double T){
	if(n>0){
		return(1.0);
	}else{
		return(soil_red_evap(psi,T));
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD, double thres, double z0snow, double *z0_ris, 
						   double *d0_ris, double *z0_z0t_ris){ 

	if(snowD>thres){
		*z0_ris=z0snow;
		*d0_ris=0.0;
		*z0_z0t_ris=0.0;
	}else{
		*z0_ris=z0;
		*d0_ris=d0;
		*z0_z0t_ris=z0_z0t;
	}
		
}
	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void check_continuity(long n, double *dw, double *wi, double *wl){

	long l;
	
	for(l=1;l<=n;l++){
		if(wl[l]+dw[l]<0) dw[l]=-wl[l];			
		if(wi[l]-dw[l]<0) dw[l]=wi[l];
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void prepare_output(double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac,
				double prec_rain, double prec_snow_atm, ENERGY *egy, WATER *wat, SNOW *snow, GLACIER *glac, LAND *land, 
				TOPO *top, SOIL *sl, METEO *met, TIMES *times, PAR *par, long r, long c, double A, double LE, double surfEB, 
				double H, double G, double Ts, double SWin, double SW, double SWbeam, double eps, double LWin, double LW, 
				double cosinc, double Precpoint){
					
	FILE *f;

	//OUTPUT DISTRIBUITI OPZIONALI
	if(par->output_balancesn>0){
		snow->MELTED->co[r][c]+=Mr_snow*par->Dt;
		snow->SUBL->co[r][c]+=(Sr_snow+Er_snow)*par->Dt;
		if(snow->type->co[r][c]==2){
			snow->t_snow->co[r][c]+=par->Dt/3600.0;
		}
	}
	
	if(par->output_balancegl>0 && par->glaclayer_max>0) glac->MELTED->co[r][c]+=Mr_glac*par->Dt;
	if(par->output_balancegl>0 && par->glaclayer_max>0) glac->SUBL->co[r][c]+=(Sr_glac+Er_glac)*par->Dt;


	if(par->output_P>0){
		wat->PrTOT_mean->co[r][c]+=wat->PrecTot->co[r][c];	/*[mm]*/
		wat->PrSNW_mean->co[r][c]+=prec_snow_atm;								/*[mm]*/
	}

	if(par->output_ET>0){
		egy->ET_mean->co[r][c]+=LE/((par->output_ET*3600.0)/(par->Dt)); //[W/m^2]
		if(par->distr_stat==1){
			if(egy->ET_max->co[r][c]<LE) egy->ET_max->co[r][c]=LE;
			if(egy->ET_min->co[r][c]>LE) egy->ET_min->co[r][c]=LE;
		}
	}
	
	if(par->output_G>0){
		egy->SEB_mean->co[r][c]+=surfEB/((par->output_G*3600.0)/(par->Dt)); //[W/m^2]
		if(par->distr_stat==1){
			if(egy->G_max->co[r][c]<surfEB) egy->G_max->co[r][c]=surfEB;
			if(egy->G_min->co[r][c]>surfEB) egy->G_min->co[r][c]=surfEB;
		}
		egy->G_snowsoil->co[r][c]+=G/((par->output_G*3600.0)/(par->Dt)); //[W/m^2]
	}
	
	
	if(par->output_H>0){
		egy->H_mean->co[r][c]+=H/((par->output_H*3600.0)/(par->Dt)); //[W/m^2]
		if(par->distr_stat==1){
			if(egy->H_max->co[r][c]<H) egy->H_max->co[r][c]=H;
			if(egy->H_min->co[r][c]>H) egy->H_min->co[r][c]=H;
		}
	}
	
	if(par->output_Ts>0){
		egy->Ts_mean->co[r][c]+=Ts/((par->output_Ts*3600.0)/(par->Dt)); /*update of surface temperature [Celsius]*/
		if(par->distr_stat==1){
			if(egy->Ts_max->co[r][c]<Ts) egy->Ts_max->co[r][c]=Ts;
			if(egy->Ts_min->co[r][c]>Ts) egy->Ts_min->co[r][c]=Ts;
		}
	}
	
	if(par->output_Rswdown>0){
		egy->Rswdown_mean->co[r][c]+=SWin/((par->output_Rswdown*3600.0)/(par->Dt));
		egy->Rswbeam->co[r][c]+=SWbeam/((par->output_Rswdown*3600.0)/(par->Dt));
		if(par->distr_stat==1){
			if(egy->Rswdown_max->co[r][c]<SWin) egy->Rswdown_max->co[r][c]=SWin;
		}
	}
	

	if(par->output_meteo>0){
		egy->Ta_mean->co[r][c]+=met->Tgrid->co[r][c]/((par->output_meteo*3600.0)/(par->Dt));
		if(par->distr_stat==1){
			if(egy->Ta_max->co[r][c]<met->Tgrid->co[r][c]) egy->Ta_max->co[r][c]=met->Tgrid->co[r][c];
			if(egy->Ta_min->co[r][c]>met->Tgrid->co[r][c]) egy->Ta_min->co[r][c]=met->Tgrid->co[r][c];
		}
		if(par->micromet==1){
			met->Vspdmean->co[r][c]+=met->Vgrid->co[r][c]/((par->output_meteo*3600.0)/(par->Dt));
			met->Vdirmean->co[r][c]+=met->Vdir->co[r][c]/((par->output_meteo*3600.0)/(par->Dt));
			met->RHmean->co[r][c]+=met->RHgrid->co[r][c]/((par->output_meteo*3600.0)/(par->Dt));
		}
	}
	
	if(par->output_Rn>0){
		egy->Rn_mean->co[r][c]+=(SW + LW) / ((par->output_Rn*3600.0)/(par->Dt)); //[W/m^2]
		
		if(fabs(SW + LW)>1500){
			f=fopen(files->co[ferr]+1,"a");
			fprintf(f,"\ntime=%10.1f r=%4ld c=%4ld Rn=%10.3f Rsw=%10.3f Rlwdiff=%10.3f albedo=%10.8f eps=%10.8fTa=%10.5f Ts=%10.5f Rsw_meas=%f sin(alpha)=%f cos(inc)=%f\n",
				times->time,r,c,SW+LW,SWin,LWin,A,
				eps,met->Tgrid->co[r][c],Ts,met->var[0][iSW],sin(egy->hsun),cosinc);
			fprintf(f,"\nH:%f LE:%f\n",H,LE);
			fclose(f);
		}
				
		egy->LW_in->co[r][c]+=LWin/((par->output_Rn*3600.0)/(par->Dt));
		egy->LW->co[r][c]+=LW/((par->output_Rn*3600.0)/(par->Dt));
		egy->SW->co[r][c]+=SW/((par->output_Rn*3600.0)/(par->Dt));
		
			
		if(par->distr_stat==1){
			if(egy->LW_max->co[r][c]<LW) egy->LW_max->co[r][c]=LW;
			if(egy->LW_min->co[r][c]>LW) egy->LW_min->co[r][c]=LW;
			if(egy->Rn_max->co[r][c]<SW + LW) egy->Rn_max->co[r][c]=(SW + LW);
			if(egy->Rn_min->co[r][c]>SW + LW) egy->Rn_min->co[r][c]=(SW + LW);
			if(egy->SW_max->co[r][c]<SW) egy->SW_max->co[r][c]=SW;
		}
	}
}
		
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_pixel(long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double Sr_soil, 
	double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, 
	double Evt, double LE, double H, double surfEB, double G, double Tg, double A, double eps, double LSAI, double LWin, 
	double SWin, double LW, double SW, double epsa, double epsa_min, double epsa_max, double SWbeam, double SWdiff,
	double Tdew, long n, double Lobukhov, PAR *par, WATER *wat, ENERGY *egy, TOPO *top, METEO *met, SNOW *snow, 
	GLACIER *glac, LAND *land, SOIL *sl, double Tpoint, double Ppoint, double Vpoint, double RHpoint, double prec_snow_atm, 
	double prec_rain_atm, double z0soil, double z0, double d0, double SWv, double LWv, double Tv, double Ts, double Hg0, 
	double Hg1, double Eg0, double Eg1, double fc, double rh, double rv, double rb, double rc, double r_uc,  
	double Hv, double LEv, double Qv, double Qg, double Qa, double Qs, double u_top, double decay, double Locc, 
	double SWup_above_v, double LWup_above_v){

	long i;
	
	if(par->state_pixel==1){
		for(i=1;i<=par->chkpt->nrh;i++){
			if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){
			
				outdata_point->co[oprecsnow][i]+=prec_snow;
				outdata_point->co[oprecrain][i]+=(prec_rain_on_soil+prec_rain_on_snow);
				outdata_point->co[orainonsnow][i]+=prec_rain_on_snow;			
				outdata_point->co[osnowover][i]+=prec_snow_atm;
				outdata_point->co[orainover][i]+=prec_rain_atm;
				
				outdata_point->co[oevapsur][i]+=Er_soil*par->Dt;	//Eg[mm]
				outdata_point->co[otrasp][i]+=Evt*(1000.0/rho_w)*par->Dt;	//Etc[mm]
				
				outdata_point->co[oLE][i]+=LE/(double)n; //ET[W/m^2]
				outdata_point->co[oH][i]+=H/(double)n; //H[W/m^2]
				outdata_point->co[oEB][i]+=surfEB/(double)n; 
				outdata_point->co[oG][i]+=G/(double)n; 
				outdata_point->co[oTg][i]+=Tg/(double)n; //Ts[C]
				
				outdata_point->co[oSWin][i]+=SWin/(double)n;
				outdata_point->co[oLWin][i]+=LWin/(double)n;
				
				outdata_point->co[oSW][i]+=SW/(double)n;
				outdata_point->co[oLW][i]+=LW/(double)n;

				outdata_point->co[oV][i]+=Vpoint/(double)n;
				outdata_point->co[oRH][i]+=RHpoint/(double)n;
				outdata_point->co[oPa][i]+=Ppoint/(double)n;
				outdata_point->co[oTa][i]+=Tpoint/(double)n;
				outdata_point->co[oTdew][i]+=Tdew/(double)n;

				outdata_point->co[oLobuk][i]+=Lobukhov/(double)n;
				//outdata_point->co[oiterLo][i]+=turbulence->co[1]/(double)n;
				
				if(epsa_min>0){
					outdata_point->co[ominLWin][i]+=(epsa_min*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
				}else{
					outdata_point->co[ominLWin][i]=UV->V->co[2];
				}
				
				if(epsa_max>0){
					outdata_point->co[omaxLWin][i]+=(epsa_max*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
				}else{
					outdata_point->co[omaxLWin][i]=UV->V->co[2];
				}
				
				outdata_point->co[oSWb][i]+=(SWbeam/(double)n);
				outdata_point->co[oSWd][i]+=(SWdiff/(double)n);
				
				if(par->micromet==1){
					outdata_point->co[oVdir][i]+=(met->Vdir->co[r][c])/(double)n;
				}else{
					outdata_point->co[oVdir][i]=UV->V->co[2];
				}
				
				outdata_point->co[oLSAI][i]+=LSAI/(double)n;
				
				outdata_point->co[oz0v][i]+=z0/(double)n;
				
				outdata_point->co[od0v][i]+=d0/(double)n;				

				outdata_point->co[omrsnow][i]+=Mr_snow*par->Dt;	//[mm]
				outdata_point->co[oersnow][i]+=Er_snow*par->Dt;	//[mm]
				outdata_point->co[osrsnow][i]+=Sr_snow*par->Dt;	//[mm]

				outdata_point->co[omrglac][i]+=Mr_glac*par->Dt;	//[mm]
				outdata_point->co[oerglac][i]+=Er_glac*par->Dt;	//[mm]
				outdata_point->co[osrglac][i]+=Sr_glac*par->Dt;	//[mm]
				
				outdata_point->co[oSWv][i]+=SWv/(double)n;
				outdata_point->co[oLWv][i]+=LWv/(double)n;
				outdata_point->co[oHv][i]+=Hv/(double)n;
				outdata_point->co[oLEv][i]+=LEv/(double)n;
				outdata_point->co[oTv][i]+=Tv/(double)n;
				outdata_point->co[oTs][i]+=Ts/(double)n;
				
				outdata_point->co[oHg0][i]+=Hg0/(double)n;
				outdata_point->co[oLEg0][i]+=Levap(Tg)*Eg0/(double)n;
				outdata_point->co[oHg1][i]+=Hg1/(double)n;
				outdata_point->co[oLEg1][i]+=Levap(Tg)*Eg1/(double)n;
				outdata_point->co[ofc][i]+=fc/(double)n;

				outdata_point->co[oQv][i]+=(Qv)/(double)n;
				outdata_point->co[oQg][i]+=(Qg)/(double)n;
				outdata_point->co[oQa][i]+=(Qa)/(double)n;
				outdata_point->co[oQs][i]+=(Qs)/(double)n;

				outdata_point->co[outop][i]+=(u_top)/(double)n;
				
				outdata_point->co[odecay][i]+=(decay)/(double)n;
				outdata_point->co[oLobukcan][i]+=(Locc)/(double)n;
				
				outdata_point->co[oSWup][i]+=SWup_above_v/(double)n;
				outdata_point->co[oLWup][i]+=LWup_above_v/(double)n;
				outdata_point->co[oEcan][i]+=(SWv+LWv-Hv-LEv)/(double)n;
				
				outdata_point->co[oHup][i]+=(H+fc*Hv)/(double)n;
				outdata_point->co[oLEup][i]+=(LE+fc*LEv)/(double)n;
				
				outdata_point->co[owcan_rain][i]+=wat->wcan_rain->co[r][c]/(double)n;
				outdata_point->co[owcan_snow][i]+=wat->wcan_snow->co[r][c]/(double)n;			
				
			}
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_basin(long n, double npixel, double prec_rain, double prec_snow, double prec_rain_atm, double prec_snow_atm, 
				  double Ta, double Tg, double Tv, double Eg, double Evt, double LE, double H, double SW, double LW, double LEv, 
				  double Hv, double SWv, double LWv, double SWin, double LWin, double Dt){

	
	outdata_basin->co[ooprecrain]+=prec_rain/(double)npixel;										
	outdata_basin->co[ooprecsnow]+=prec_snow/(double)npixel;			
	outdata_basin->co[oorainover]+=prec_rain_atm/(double)npixel;													
	outdata_basin->co[oosnowover]+=prec_snow_atm/(double)npixel;		
	outdata_basin->co[ooTa]+=Ta/(double)(n*npixel);															
	outdata_basin->co[ooTg]+=Tg/(double)(n*npixel);	
	outdata_basin->co[ooTv]+=Tv/(double)(n*npixel);			
	outdata_basin->co[ooevapsur]+=Eg*Dt/(double)npixel;															
	outdata_basin->co[ootrasp]+=Evt*Dt/(double)npixel;														
	outdata_basin->co[ooLE]+=LE/(double)(n*npixel);																
	outdata_basin->co[ooH]+=H/(double)(n*npixel);																
	outdata_basin->co[ooSW]+=SW/(double)(n*npixel);													
	outdata_basin->co[ooLW]+=LW/(double)(n*npixel);	
	outdata_basin->co[ooLEv]+=LEv/(double)(n*npixel);	
	outdata_basin->co[ooHv]+=Hv/(double)(n*npixel);	
	outdata_basin->co[ooSWv]+=SWv/(double)(n*npixel);	
	outdata_basin->co[ooLWv]+=LWv/(double)(n*npixel);	
	outdata_basin->co[ooSWin]+=SWin/(double)(n*npixel);	
	outdata_basin->co[ooLWin]+=LWin/(double)(n*npixel);		
	
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

				 

	