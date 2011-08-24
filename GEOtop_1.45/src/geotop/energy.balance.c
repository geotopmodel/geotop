
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

    
//Author: Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch
//Contents: Energy balance (and also mass balance for snow and glacier)

#include "constants.h"
#include "struct.geotop.h"
#include "energy.balance.h"
#include "meteo.h"
#include "snow.h"
#include "pedo.funct.h"
#include "vegetation.h"
#include "radiation.h"
#include "turbulence.h"
#include "../libraries/math/util_math.h"
#include "meteodistr.h"
#include "times.h"
#include "tables.h"
#include "blowingsnow.h"
#include "meteodata.h"

extern long number_novalue, number_absent;
extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern char *logfile;
extern long Nl, Nr, Nc;
extern double **outdata_point, *outdata_basin;
extern long i_sim;

#define M 1
#define ni 1.E-4
#define num_iter_after_which_surfenergy_balance_not_recalculated 50
#define Tmin_surface_below_which_surfenergy_balance_recalculated -50
#define num_iter_after_which_only_neutrality 10
#define Csnow_at_T_greater_than_0 1.E20
#define ratio_max_storage_RAIN_over_canopy_to_LSAI 0.1
#define ratio_max_storage_SNOW_over_canopy_to_LSAI 5.0
#define Tmax 100.0
#define Tmin -90.0

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void energy_balance(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, 
					SNOW *snow, GLACIER *glac, CHANNEL *cnet)

{

	/*r = row index
	 c = column index
	 l = layer index
	 ns = snow layer number
	 ng = glacier layer number
	 */ 
	long r, c, i, j, l, ns, ng=0;
	
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
	 SWupabove_v = shortwave above vegetation directed upwards
	 */
	double SWin, SW, SWbeam, SWdiff, SWv_vis, SWv_nir, SWg_vis, SWg_nir, cosinc, avis_b, avis_d, anir_b, anir_d, E0, Et, Delta;
	double SWupabove_v_vis, SWupabove_v_nir, SWup_above_v, Asurr_ave=0.0;
	double avis_ground, anir_ground;
	short SWb_yes;
	double tauatm_sinhsun;
	
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
	double LWin, LW, LWv, epsa, eps, epsa_min, epsa_max, LWupabove_v, Tsurr_ave=0.0;
	
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
	 Sr_glac = glacier sublimation rate
	 k_snowred = reducing factor of snow depth for steep surfaces*/
	double snowD, RainOnSnow, fsnow, fsnownew, fsnowcan, Melt_snow, Evap_snow, Melt_glac=0., Evap_glac=0.;
	double k_snowred, maxSWE, ksnow, ksoil;
	
	 /*Precipitation variables
	  Prain_over, Psnow_over: rain and snow precipitation ABOVE the canopy
	  Prain, Psnow: rain and snow precipitation BELOW the canopy
	  drip_rain, drip_snow: rain and snow precipitation not intercepted or unloaded by canopy
	  max_wcan_rain, max_wcan_snow: max storage of rain and snow on the canopy*/
	double Prain_over, Psnow_over, Prain=0.0, Psnow=0.0, drip_rain, drip_snow, max_wcan_rain = 0., max_wcan_snow = 0., theta_sup;
		
	/*Roughness lengths, without suffix means of the surface, veg means vegetation, _ means ratio*/
	double z0, z0veg = 0., d0, d0veg = 0., z0_z0t, hveg = 0.;
	
	/*Aerodyamic resistances, h=heat v=water vapor, c=canopy for water vapour, b=boundary layer canopy i.e. canopy for heat, 
	 uc=undercanopy=surface below the canopy*/
	double rh, rv, rc, rb, ruc, Lobukhov;
	
	/*Vegetation characteristics
	 fc = canopy fraction
	 Locc = below canopy Monin-Obukhov length
	 utop = wind speed at the canopy height
	 Qv = saturated specific humidity at the canopy*/
	 double fc, fc0, decaycoeff, Locc, u_top, Qv;
	
	/*soil characteristics
	 Evap_soil = evaporation from soil surface*/
	double Evap_soil;
	
	/*meteo
	 Tdew=dew temperature
	 Qa=specific humidity of air
	 zmeas_T=measurement height of T and RH [m]
	 zmeas_U=measurement height of wind[m]*/
	double ea, Tdew, Qa, RHpoint, Vpoint, Ppoint, Tpoint, Precpoint, zmeas_T, zmeas_u;	
	 
	/*others
	 Ts, Qs temperature and specific humidity of the canopy air (the CLM uses the term "surface" for that)
	 Qg specific humidity at the soil/snow surface*/
	double Ts, Qs, Qg;

	//soil type
	long sy;
	
	//land use
	short lu, lux;	
	
	//success of point energy balance 
	short sux;	//(1=yes, 0=no)
	double te, Dt;
	long k, kk, n, maxk=-1;
	
	//parameter interpolation
	static long line_interp;
	
	//plotting variables
	double W, Dtplot;
	FILE *f;
	
	//time
	double JDbegin, JDend;
		
	//open logfile
	f = fopen(logfile, "a");
		
	//calculation to be done before plotting maps
	if(times->time==0) times->iplot=1;
	if(times->JD_plots->nh > 1 && times->iplot<=times->JD_plots->nh){
		i=times->iplot;
		j=2*i-1;
		if( par->init_date->co[i_sim]+times->time/86400.+1.E-5>=times->JD_plots->co[j] && 
		    par->init_date->co[i_sim]+(times->time+par->Dt)/86400.-1.E-5<=times->JD_plots->co[j+1] ){
			
			Dtplot=(times->JD_plots->co[j+1]-times->JD_plots->co[j])*86400.;
			W = par->Dt / Dtplot;				
			fprintf(f,"Saving plot number %ld Weight:%f \n",i,W);
			printf("Saving plot number %ld Weight:%f \n",i,W);
		}else {
			W = 0.0;
		}

	}else {
		W = 0.0;
	}

	//initialization static variables
	if(times->time==0) line_interp=0;
			
	//SUN
	sun( convert_tfromstart_JDfrom0(times->time+0.5*par->Dt, par->init_date->co[i_sim]), &E0, &Et, &Delta );
	egy->sun[1] = Delta;
	
	//VEGETATION
	for(lux=1; lux<=par->n_landuses; lux++) {
		if(par->vegflag->co[lux]==1) meteo_interp(0, 1, &line_interp, land->vegparv[lux-1], land->vegpars[lux-1], land->NumlinesVegTimeDepData[lux-1], 
							jdvegprop+1, 0, par->init_date->co[i_sim]+times->time/86400., par->init_date->co[i_sim]+(times->time+par->Dt)/86400.);
	}
	
	//COMPUTATIONS FOR EACH CELL	
	for (i=1; i<=par->total_pixel+par->total_channel; i++) {
		
		if (i<=par->total_channel) {
			//CHANNEL
			r = cnet->r->co[i];
			c = cnet->c->co[i];
			sy = cnet->soil_type->co[i];
			lu=(short)land->LC->co[r][c];
			
		}else {
			//LAND
			r = top->rc_cont->co[i-par->total_channel][1];
			c = top->rc_cont->co[i-par->total_channel][2];
			sy = sl->type->co[r][c];		
			lu=(short)land->LC->co[r][c];
		}
						
		//initialization of cumulated water volumes
		if (i<=par->total_channel) {
			for (l=1; l<=Nl; l++) {
				cnet->ET->co[l][i] = 0.0;
			}
		}else{
			wat->Pnet->co[r][c] = 0.0;
			for (l=1; l<=Nl; l++) {
				sl->ET->co[l][r][c] = 0.0;
			}
		}
		
		//integration of heat equation with adjusting time step
		Dt = par->Dt;
		te = 0.0;
		k = 0;
		
		do{
					
			if (k > maxk) {
				
				maxk = k;				
				kk = (long)pow(2.,k);

				for (n=kk; n<=(long)pow(2.,k+1)-1; n++) {
										
					JDbegin = convert_tfromstart_JDfrom0(times->time+(n-kk)*Dt, par->init_date->co[i_sim]);
					JDend = convert_tfromstart_JDfrom0(times->time+(n-kk+1)*Dt, par->init_date->co[i_sim]);
					
					if (par->point_sim!=1) {
						//SUN AND SHADOWING	
						egy->sun[0] = par->latitude*Pi/180.;
						egy->sun[2] = (par->longitude*Pi/180. - par->ST*Pi/12. + Et)/omega;	
						egy->hsun->co[n] = adaptiveSimpsons2(SolarHeight__, egy->sun, JDbegin, JDend, 1.E-4, 10) / (JDend - JDbegin);
						//egy->hsun->co[n] = SolarHeight_( 0.5*(JDbegin+JDend), egy->sun );
						egy->dsun->co[n] = adaptiveSimpsons2(SolarAzimuth__, egy->sun, JDbegin, JDend, 1.E-4, 10) / (JDend - JDbegin);
						//egy->dsun->co[n] = SolarAzimuth_( 0.5*(JDbegin+JDend), egy->sun );
						egy->sinhsun->co[n] = adaptiveSimpsons2(Sinalpha_, egy->sun, JDbegin, JDend, 1.E-4, 10) / (JDend - JDbegin);
						//egy->sinhsun->co[n] = sin(Fmax(egy->hsun, 0.05));
						if(n==1) shadow_haiden(top, egy->hsun->co[n], egy->dsun->co[n], land->shadow);
					}
					
					//CLOUDINESS
					find_actual_cloudiness(&(met->tau_cloud->co[n]), &(met->tau_cloud_av->co[n]), &(met->tau_cloud_yes->co[n]), &(met->tau_cloud_av_yes->co[n]),
										   met, JDbegin, JDend, Delta, E0, Et, par->ST, egy->Asurr->co[r][c]);
					
					//METEO
					if(n>1) meteo_distr(0, n, met->line_interp_WEB, met->line_interp_WEB_LR, met, wat, top, par, JDbegin, JDend);
					
				}
			}
			
			JDbegin = convert_tfromstart_JDfrom0(times->time+te, par->init_date->co[i_sim]);
			JDend = convert_tfromstart_JDfrom0(times->time+te+Dt, par->init_date->co[i_sim]);

			n = (long)pow(2.,k) + floor(te/Dt);
			
			//METEO
			Tpoint=met->Tgrid->co[n][r][c];
			Ppoint=met->Pgrid->co[n][r][c];	
			RHpoint=met->RHgrid->co[n][r][c];
			Vpoint=met->Vgrid->co[n][r][c];
			
			Precpoint=wat->PrecTot->co[n][r][c]*cos(top->slope->co[r][c]*Pi/180.);
			
			//SNOW
			snowD=DEPTH(r, c, snow->S->lnum, snow->S->Dzl);			
			ns=snow->S->lnum->co[r][c];
			
			//vegetation parameters
			if (i>par->total_channel) {
				
				for(j=1;j<=jdvegprop;j++){
					if( (long)land->vegparv[lu-1][j-1+1] != number_novalue ){
						land->vegpar->co[j] = land->vegparv[lu-1][j-1+1];
						if(j==jdroot) root(land->root_fraction->nch, land->vegpar->co[jdroot], 0.0, sl->pa->co[1][jdz], land->root_fraction->co[lu]);
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
					ng=glac->G->lnum->co[r][c];
					if(ng>0){ 
						land->vegpar->co[jdLSAI]=0.0; 
						fc=0.0; 
					}
				}
				
			}else {
				
				initialize_doublevector(land->vegpar, 0.0);
				fc = 0.0;
				
			}
			
			//METEOROLOGICAL COMPUTATIONS
			
			//temperature and wind velocity measurement heights
			zmeas_u=met->st->Vheight->co[1];
			zmeas_T=met->st->Theight->co[1];
			
			//RAIN AND SNOW PRECIPITATION [mm]
			//convert total precipitation to [mm]
			Precpoint*=(Dt/3600.0);	//from [mm/h] to [mm]
			
			ea = RHpoint * SatVapPressure(Tpoint, Ppoint);//Vapour Pressure [mbar]
			Qa = SpecHumidity(ea, Ppoint);//Specific Humidity
			Tdew = TfromSatVapPressure(ea, Ppoint);//Dew Temperature
			
			//distinguish between rain and snow
			if(par->dew==1){
				//on the base of the dew temperature of the air
				part_snow(Precpoint, &Prain_over, &Psnow_over, Tdew, par->T_rain, par->T_snow);
			}else{
				//on the base of the temperature of the air
				part_snow(Precpoint, &Prain_over, &Psnow_over, Tpoint, par->T_rain, par->T_snow);
			}
						
			//Adjusting snow precipitation in case of steep slope (contribution by Stephan Gruber)
			if (par->snow_curv > 0 && top->slope->co[r][c] > par->snow_smin){
				if (top->slope->co[r][c] <= par->snow_smax){
					k_snowred = ( exp(-pow(top->slope->co[r][c] - par->snow_smin, 2.)/par->snow_curv) -
								 exp(-pow(par->snow_smax, 2.)/par->snow_curv) );
				}else{
					k_snowred = 0.0;
				}
				Psnow_over *= k_snowred;
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
			}
			
			//SHORTWAVE RADIATION
			//initialization of shortwave radiation absorbed by soil
			SW=0.0;
			SWup_above_v=0.0;
			
			//if averaged cloud transmissivity is not available it is estimated through this meteodistr subroutine (not very reliable)
			if(met->tau_cloud_av_yes->co[n]==0) met->tau_cloud_av->co[n] = 1. - 0.71*find_cloudfactor(Tpoint, RHpoint, top->Z0->co[r][c], met->LRv[ilsTa], met->LRv[ilsTdew]);//Kimball(1928)
			
			//in case of shortwave data not available
			if(met->tau_cloud_yes->co[n]==0) met->tau_cloud->co[n] = met->tau_cloud_av->co[n];	
			
			//calculation of SWin
			if(par->point_sim==1){
				egy->sun[0] = top->latitude->co[r][c]*Pi/180.;
				egy->sun[2] = (top->longitude->co[r][c]*Pi/180. - par->ST*Pi/12. + Et)/omega;	
				egy->hsun->co[n] = adaptiveSimpsons2(SolarHeight__, egy->sun, JDbegin, JDend, 1.E-4, 10) / (JDend - JDbegin);
				//egy->hsun->co[n] = SolarHeight_( 0.5*(JDbegin+JDend), egy->sun );
				egy->dsun->co[n] = adaptiveSimpsons2(SolarAzimuth__, egy->sun, JDbegin, JDend, 1.E-4, 10) / (JDend - JDbegin);
				//egy->dsun->co[n] = SolarAzimuth_( 0.5*(JDbegin+JDend), egy->sun );
				egy->sinhsun->co[n] = adaptiveSimpsons2(Sinalpha_, egy->sun, JDbegin, JDend, 1.E-4, 10) / (JDend - JDbegin);
				//egy->sinhsun->co[n] = sin(Fmax(egy->hsun, 0.05));
				land->shadow->co[r][c]=shadows_point(top->horizon_height[top->horizon_point->co[r][c]-1], top->horizon_numlines[top->horizon_point->co[r][c]-1], 
											   egy->hsun->co[n]*180./Pi, egy->dsun->co[n]*180/Pi, Tol_h_mount, Tol_h_flat);
			}
			egy->sun[3] = RHpoint;
			egy->sun[4] = Tpoint;
			egy->sun[5] = Ppoint;
			egy->sun[6] = top->slope->co[r][c]*Pi/180.;
			egy->sun[7] = top->aspect->co[r][c]*Pi/180.;	
			shortwave_radiation(JDbegin, JDend, egy->sun, egy->sinhsun->co[n], E0, top->sky->co[r][c], egy->Asurr->co[r][c], 
								met->tau_cloud->co[n], land->shadow->co[r][c], &SWbeam, &SWdiff, &cosinc, &tauatm_sinhsun, &SWb_yes);
			SWin=SWbeam+SWdiff;
			
			if (n==1 && i>par->total_channel) {
				if (SWb_yes == 1) {
					egy->nDt_sun->co[r][c] ++;
				}else if (SWb_yes == 0) {
					egy->nDt_shadow->co[r][c] ++;
				}
			}
			
			//albedo
			if (i<=par->total_channel) {
				theta_sup = cnet->th->co[1][i];
			}else {
				theta_sup = sl->th->co[1][r][c];
			}
			
			avis_ground = find_albedo(land->ty->co[lu][ja_vis_dry], land->ty->co[lu][ja_vis_sat], theta_sup, sl->pa->co[sy][jres][1],  sl->pa->co[sy][jsat][1]);
			anir_ground = find_albedo(land->ty->co[lu][ja_nir_dry], land->ty->co[lu][ja_nir_sat], theta_sup, sl->pa->co[sy][jres][1],  sl->pa->co[sy][jsat][1]);										 
			if(snowD>0){
				if(i>par->total_channel) update_snow_age(Psnow_over, snow->S->T->co[ns][r][c], par->Dt, &(snow->dimens_age->co[r][c]), &(snow->nondimens_age->co[r][c]));
				avis_b=snow_albedo(avis_ground, snowD, par->aep, par->avo, 0.2, snow->nondimens_age->co[r][c], cosinc, (*Fzen));
				avis_d=snow_albedo(avis_ground, snowD, par->aep, par->avo, 0.2, snow->nondimens_age->co[r][c], cosinc, (*Zero));
				anir_b=snow_albedo(anir_ground, snowD, par->aep, par->airo, 0.5, snow->nondimens_age->co[r][c], cosinc, (*Fzen));
				anir_d=snow_albedo(anir_ground, snowD, par->aep, par->airo, 0.5, snow->nondimens_age->co[r][c], cosinc, (*Zero));
			}else{
				avis_b=avis_ground;
				avis_d=avis_ground;
				anir_b=anir_ground;
				anir_d=anir_ground;
			}
			
			//shortwave absorbed by soil (when vegetation is not present)
			SW += (1.0-fc)*( SWdiff*(1.0-0.5*avis_d-0.5*anir_d) + SWbeam*(1.0-0.5*avis_b-0.5*anir_b) );
			SWup_above_v += (1.0-fc)*( SWdiff*(0.5*avis_d+0.5*anir_d) + SWbeam*(0.5*avis_b+0.5*anir_b) );
			
			//shortwave radiation absorbed by canopy 
			if(fc>0){
				if(egy->hsun->co[n]>0){
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
			SW=flux(1, iSWn, met->var, 1.0, SW);
			
			//mean albedo shortwave radiation (as proxy of albedo of surrounding surfaces)
			if(n==1 && i>par->total_channel && egy->hsun->co[n]>0 && met->tau_cloud->co[n]>0){
				if (par->point_sim == 1) {
					egy->Asurr->co[r][c] = SWup_above_v/(Isc*E0*met->tau_cloud->co[n]*tauatm_sinhsun);
				}else {
					Asurr_ave += (SWup_above_v/(Isc*E0*met->tau_cloud->co[n]*tauatm_sinhsun))/(double)par->total_pixel;
				}
			}
			
			//Extinction coefficient for SW in the snow layers
			rad_snow_absorption(r, c, egy->SWlayer, SW, snow->S);			
			
			//LONGWAVE RADIATION
			//soil-snow emissivity
			if(snowD>10){
				eps=par->epsilon_snow;
			}else{
				eps=land->ty->co[lu][jemg];
			}			
			longwave_radiation(par->state_lwrad, ea, RHpoint, Tpoint, met->tau_cloud_av->co[n], &epsa, &epsa_max, &epsa_min);
			LWin=top->sky->co[r][c]*epsa*SB(Tpoint) + (1.-top->sky->co[r][c])*eps*SB(egy->Tgskinsurr->co[r][c]);
			
			//if incoming LW data are available, they are used (priority)
			LWin=flux(met->nstlrad, iLWi, met->var, 1.0, LWin);
			
			//roughness lengths
			update_roughness_soil(land->ty->co[lu][jz0], 0.0, 0.0, snowD, land->ty->co[lu][jz0thressoil], par->z0_snow, &z0, &d0, &z0_z0t);
			if(fc>0) update_roughness_veg(land->vegpar->co[jdHveg], snowD, zmeas_u, zmeas_T, &z0veg, &d0veg, &hveg);							
			
			//SOIL AND SNOW PROPERTIES
			if (i<=par->total_channel) {
				egy->Temp->co[par->nsurface] = cnet->Tgskin->co[i];
			}else{
				egy->Temp->co[par->nsurface] = egy->Tgskin->co[r][c];
			}
						
			for(l=1;l<=Nl+ns+ng;l++){
				if(l<=ns){	//snow
					egy->Dlay->co[l] = 1.E-3*snow->S->Dzl->co[ns+1-l][r][c];
					egy->wliq->co[l] = snow->S->w_liq->co[ns+1-l][r][c];
					egy->wice->co[l] = snow->S->w_ice->co[ns+1-l][r][c];
					egy->Temp->co[l] = snow->S->T->co[ns+1-l][r][c];
				}else if(l<=ns+ng){   //glacier
					egy->Dlay->co[l] = 1.E-3*glac->G->Dzl->co[ns+ng+1-l][r][c];
					egy->wliq->co[l] = glac->G->w_liq->co[ns+ng+1-l][r][c];
					egy->wice->co[l] = glac->G->w_ice->co[ns+ng+1-l][r][c];
					egy->Temp->co[l] = glac->G->T->co[ns+ng+1-l][r][c];
				}else{	//soil
					if (i<=par->total_channel) {
						egy->Dlay->co[l] = 1.E-3*sl->pa->co[sy][jdz][l-ns-ng];
						egy->wliq->co[l] = cnet->th->co[l-ns-ng][i]*egy->Dlay->co[l]*rho_w;
						egy->wice->co[l] = cnet->thice->co[l-ns-ng][i]*egy->Dlay->co[l]*rho_w;
						egy->Temp->co[l] = cnet->T->co[l-ns-ng][i];
					}else {
						egy->Dlay->co[l] = 1.E-3*sl->pa->co[sy][jdz][l-ns-ng];					
						egy->wliq->co[l] = sl->th->co[l-ns-ng][r][c]*egy->Dlay->co[l]*rho_w;
						egy->wice->co[l] = sl->thice->co[l-ns-ng][r][c]*egy->Dlay->co[l]*rho_w;
						egy->Temp->co[l] = sl->T->co[l-ns-ng][r][c];
					}
				}					
			}			
						
			//ENERGY BALANCE	
			
			sux=PointEnergyBalance(f, Dt, i, r, c, sy, egy, land, sl, cnet, par, ns, ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg, d0veg, 1.0, hveg, 
								   Vpoint, Tpoint, Qa, Ppoint, met->LRv[ilsTa], eps, fc, land->vegpar->co[jdLSAI], 
								   land->vegpar->co[jddecay0], &(wat->wcan_rain->co[r][c]), max_wcan_rain, &(wat->wcan_snow->co[r][c]), 
								   max_wcan_snow, SWin, LWin, SWv_vis+SWv_nir, &LW, &H, &E, &LWv, &Hv, &LEv, &Etrans, &Ts, &Qs, Hadv, 
								   &Hg0, &Hg1, &Eg0, &Eg1, &Qv, &Qg, &Lobukhov, &rh, &rv, &rb, &rc, &ruc, &u_top, &decaycoeff, &Locc, 
								   &LWupabove_v, Prain);	
									
			if (sux == 0) {	
				
				k++;
				Dt /= 2.;
				
			}else{
					
				te += Dt;
										
				if(i<=par->total_channel){

					update_soil_channel(par->nsurface, ns+ng, i, r, c, sy, fc, Dt, egy, sl, cnet->T, cnet->P, cnet->th, cnet->thice, cnet->ET, cnet->Tgskin);
					
				}else {	
												
					LE = E*latent(egy->Temp->co[1],Levap(egy->Temp->co[1]));
					surfEB = SW + LW - H - LE;
						
					if(ns==0){
						G = surfEB;
					}else{
						ksnow = calc_k(ns, r, c, ns, ng, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Temp->co, egy->Dlay->co, (*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par);		
						ksoil = calc_k(ns+1, r, c, ns, ng, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Temp->co, egy->Dlay->co, (*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par);		
						G = Harmonic_Mean(egy->Dlay->co[ns], egy->Dlay->co[ns+1], ksnow, ksoil)*(egy->Temp->co[ns]-egy->Temp->co[ns+1])/(0.5*egy->Dlay->co[ns]+0.5*egy->Dlay->co[ns+1]);				
					}		
					
					//SOIL
					update_soil_land(par->nsurface, ns+ng, r, c, sy, fc, Dt, egy, sl, sl->T, sl->P, sl->th, sl->thice, sl->ET, egy->Tgskin);

					//assign evaporation to soil
					if (ng==0 && snow->S->lnum->co[r][c]==0){
						Evap_soil = E*Dt;
					}else {
						Evap_soil = 0.0;
					}
						
					//GLACIER
					if(par->glaclayer_max>0){
							
						//assign evaporation to glacier
						if(snow->S->lnum->co[r][c]==0 && glac->G->lnum->co[r][c]>0){
							Evap_glac = E*Dt;
						}else {
							Evap_glac = 0.0;
						}
							
						//glacier melting
						WBglacier(ns, r, c, glac->G, &Melt_glac, par, egy->wice->co, egy->wliq->co, egy->Temp->co, Evap_glac);
							
						//adjust glacier layers
						snow_layer_combination(par->alpha_snow, r, c, glac->G, Tpoint, par->glaclayer_inf, par->Dmin_glac, par->Dmax_glac, 1.E10, times->time, f);							

						//convert glacier to snow, when the former is too thin
						glac2snow(par->alpha_snow, r, c, snow->S, glac->G);
							
					}
						
					//SNOW
					//assign evaporation to snow
					if(snow->S->lnum->co[r][c]>0){
						Evap_snow = E*Dt;
					}else {
						Evap_snow = 0.0;
					}
						
					//snow melting, compactation and percolation
					WBsnow(Dt, r, c, snow->S, &Melt_snow, &RainOnSnow, par, top->slope->co[r][c], Prain, egy->wice->co, egy->wliq->co, egy->Temp->co, Evap_snow);
					
					//add new snow
					if(Psnow>0) new_snow(par->alpha_snow, r, c, snow->S, Psnow, Psnow*rho_w/rho_newlyfallensnow(Vpoint, Tpoint, Tfreezing), Tpoint);
					
					//adjust snow layers
					if (par->point_sim == 1) {
						maxSWE = par->maxSWE->co[r][c];
					}else {
						maxSWE = 1.E10;
					}		
					
					snow_layer_combination(par->alpha_snow, r, c, snow->S, Tpoint, par->snowlayer_inf, par->Dmin, par->Dmax, maxSWE, times->time, f);

					//NET PRECIPITATION
					wat->Pnet->co[r][c] += (Melt_snow + Melt_glac + Prain);
					
					//VEGETATION
					if( land->vegpar->co[jdLSAI]>=LSAIthres && ng==0 ){
						
						snowD = DEPTH(r, c, snow->S->lnum, snow->S->Dzl);
						
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
							new_snow(par->alpha_snow, r, c, snow->S, wat->wcan_snow->co[r][c]*(fc0-fc), wat->wcan_snow->co[r][c]*(fc0-fc)*1000./300, sl->Tv->co[r][c]);
							wat->Pnet->co[r][c] +=  wat->wcan_rain->co[r][c]*(fc0-fc);
						}
						
						if (fc<1.E-6) {
							wat->wcan_rain->co[r][c]=0.0;
							wat->wcan_snow->co[r][c]=0.0;
						}
					}
						
					//OUTPUT 
					if(par->point_sim!=1) prepare_output(f, Dt, Evap_snow, Melt_snow, Evap_glac, Melt_glac, Prain, Psnow, egy, wat, snow, glac, land, top, 
								sl, met, times, par, r, c, LE, surfEB, H, G, egy->Tgskin->co[r][c], SWin, SW, SWbeam, eps, LWin, LW, cosinc);
						
					output_pixel(n, Dt, par->Dtplot_point->co[i_sim], r, c, Psnow, Prain-RainOnSnow, RainOnSnow, Psnow_over, Prain_over, 
								 Evap_soil, Melt_snow, Evap_snow, Melt_glac, Evap_glac, Etrans*Dt, LE, H, surfEB, G, egy->Tgskin->co[r][c], eps, 
								 land->vegpar->co[jdLSAI], LWin, SWin, LW, SW, epsa, epsa_min, epsa_max, SWbeam, SWdiff, Tdew,  
								 Lobukhov, par, wat, egy, top, met, snow, glac, land, sl, Tpoint, Ppoint, Vpoint, RHpoint, z0, z0veg, d0veg, (SWv_vis+SWv_nir), 
								 LWv, sl->Tv->co[r][c], Ts, Hg0, Hg1, Eg0, Eg1, fc, rh, rv, rb, rc, ruc, Hv, LEv, Qv, Qg, Qa, Qs, u_top, decaycoeff, Locc, SWup_above_v, LWupabove_v);
						
					output_basin(Dt, par->Dtplot_basin->co[i_sim], par->total_pixel, Prain, Psnow, Prain_over, Psnow_over, Tpoint, egy->Tgskin->co[r][c], 
								 sl->Tv->co[r][c], Evap_soil, Etrans*Dt, LE, H, SW, LW, fc*LEv, fc*Hv, fc*(SWv_vis+SWv_nir), fc*LWv, SWin, LWin);
						
					output_map_plots(n, r, c, W*Dt/par->Dt, par, egy, met, snow, H, LE, fc*Hv, fc*LEv, SWin, SW, fc*(SWv_vis+SWv_nir), LWin, LW, fc*LWv, Ts, egy->Tgskin->co[r][c], sl->Tv->co[r][c]);							
					
				}
			}
										
		}while (te < par->Dt && k <= par->max_times_halving_time_step_en);
			
		if(sux==0){
			fprintf(f, "Energy balance not converging (r:%ld c:%ld Dt:%e). GEOtop is closed.",r,c,Dt);
			fclose(f);
			printf("Energy balance not converging (r:%ld c:%ld Dt:%e). GEOtop is closed.",r,c,Dt);
			stop_execution();
			t_error("Energy balance not converging");
		}
										
		//NET PRECIPITATION
		if (i>par->total_channel) {		
			if (par->point_sim == 1) {
				egy->Tgskinsurr->co[r][c] = egy->Tgskin->co[r][c];
			}else {
				Tsurr_ave += egy->Tgskin->co[r][c]/(double)par->total_pixel;				
			}

		}
	}
	
		
	//Properties of the surroundings
	if (par->point_sim != 1) {
		for (i=1; i<=par->total_pixel; i++) {
			r = top->rc_cont->co[i][1];
			c = top->rc_cont->co[i][2];
			egy->Tgskinsurr->co[r][c] = Tsurr_ave;
			if(egy->hsun->co[1]>0 && met->tau_cloud->co[1]>0) egy->Asurr->co[r][c] = Asurr_ave;
		}
	}
		
	//close logfile
	fclose(f);
	
	
}

//end of "energy_balance" subroutine

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short PointEnergyBalance(FILE *f, double Dt, long i, long r, long c, long sy, ENERGY *egy, LAND *land, SOIL *sl, CHANNEL *cnet, PAR *par, 
			long ns, long ng, double zmu, double zmT, double z0s, double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, 
			double v, double Ta, double Qa, double P, double LR, double eps, double fc, double LSAI, double decaycoeff0, double *Wcrn, 
			double Wcrnmax, double *Wcsn, double Wcsnmax, double SWin, double LWin, double SWv,double *LW, double *H, double *E, double *LWv, 
			double *Hv, double *LEv, double *Etrans, double *Ts, double *Qs, double Hadd, double *Hg0, double *Hg1, double *Eg0, double *Eg1, 
			double *Qv, double *Qg, double *Lobukhov, double *rh, double *rv, double *rb, double *rc, double *ruc, double *u_top, 
			double *decay, double *Locc, double *LWup_above_v, double Prain){

	
	short iter_close, iter_close2, lu=land->LC->co[r][c], flagTmin=0;

	long l, m, cont=0, cont2, n=Nl+ns+ng;
	long cont_lambda_min=0;
	double dH_dT, dE_dT, EB, dEB_dT, EB0, Tg, Tg0, psim, psi0;
	double res, res0[3], res_av, res_prev[M], lambda[3], C0, C1, th0, th1, kbb0, kbb1, k=0., kn;
	double Qg0, Tv0, dWcsn=0.0, dWcrn=0.0, rh_g, rv_g;
	
	
	//Initializes vectors	
	for(l=par->nsurface;l<=n;l++){ 
		egy->T0->co[l] = egy->Temp->co[l];
	}

	if (i<=par->total_channel) {
		psi0 = cnet->P->co[0][i];
		for(l=1;l<=Nl;l++){
			//mass that melts (if>0) or freezes (if<0) 
			egy->deltaw->co[l] = 0.0;
			//water content to be modified at each iteration
			egy->THETA->co[l] = cnet->th->co[l][i];		
			//total pressure (=pressure of the water would have if it was all in liquind state)
			psim=psi_teta(cnet->th->co[l][i]+cnet->thice->co[l][i], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], 
						  sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
			//max temperature at which the first particle of ice comes up
			egy->Tstar->co[l]=Fmin(psim/(1000.0*Lf/(g*(Tfreezing+tk))), 0.0);
		}		
	}else {
		psi0 = sl->P->co[0][r][c];
		for(l=1;l<=Nl;l++){
			//mass that melts (if>0) or freezes (if<0) 
			egy->deltaw->co[l] = 0.0;
			//water content to be modified at each iteration
			egy->THETA->co[l] = sl->th->co[l][r][c];	
			//total pressure (=pressure of the water would have if it was all in liquind state)
			psim=psi_teta(sl->th->co[l][r][c]+sl->thice->co[l][r][c], 0.0, sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], 
						  sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
			//max temperature at which the first particle of ice comes up
			egy->Tstar->co[l]=Fmin(psim/(1000.0*Lf/(g*(Tfreezing+tk))), 0.0);
		}
		
	}

	
	//Surface skin temperature used to compute the surface energy fluxes
	Tg = egy->Temp->co[par->nsurface];
	Tg0 = Tg;

	//Qg0 is the specific humidity at the soil surface at the beginning of the time step (assumed saturated)
	Qg0 = SpecHumidity(SatVapPressure(Tg, P), P);

	//canopy temperature at the beginning of the time step
	Tv0=sl->Tv->co[r][c];
					
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
	EnergyFluxes(f, Tg, r, c, ns+ng, Tg0, Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, 
				 P, LR, psi0, eps, fc, LSAI, decaycoeff0, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, 
				 egy->THETA->co, sl->pa->co[sy], land->ty->co[lu], land->root_fraction->co[lu], par, egy->soil_transp_layer, 
				 SWin, LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans, &(sl->Tv->co[r][c]), Qv, Ts, Qs, Hg0, Hg1, 
				 Eg0, Eg1, Lobukhov, rh, rv, rc, rb, ruc, &rh_g, &rv_g, Qg, u_top, decay, Locc, LWup_above_v,
				 egy->Temp->co, egy->soil_evap_layer_bare, egy->soil_evap_layer_bare);
	
	EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
	dEB_dT = -eps*dSB_dT(Tg) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;
	
	EB0 = EB;
					
	//Calculate -F(x0), given that the system to solve is F'(x0) * (x-x0) = -F(x0)
	
	//F(T) = diag(egy->Fenergy(T)) + K(T)*T
	
	for(l=par->nsurface;l<=n;l++){
		egy->Fenergy->co[l] = 0.0;
	}

	for(l=1;l<=n;l++){
		
		//conductive M-matrix (lower diagonal part of this M-matrix is stored in a vector)
		if(l>1) kn = k;
		k = calc_k(l, r, c, ns, ng, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Temp->co, egy->Dlay->co, (*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par);		
		if(l>1){
			egy->Kth0->co[l-1] = -k * kn / ( k * 0.5 * egy->Dlay->co[l-1] + kn * 0.5 * egy->Dlay->co[l] );	// (-k/dz)
			egy->Kth1->co[l-1] = egy->Kth0->co[l-1];
		}else if (l>par->nsurface) {
			egy->Kth0->co[l-1] = -k / ( 0.5 * egy->Dlay->co[l] );
			egy->Kth1->co[l-1] = egy->Kth0->co[l-1];
		}
		
		//diagonal part of F due to heat capacity and boundary condition (except conduction)
		C0 = calc_C(l, r, c, ns+ng, 0.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
		C1 = calc_C(l, r, c, ns+ng, 1.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
		if(l<=ns+ng && fabs(egy->deltaw->co[l]-egy->wice->co[l])<1.E-7) C1 = Csnow_at_T_greater_than_0;
		egy->Fenergy->co[l] += ( Lf*egy->deltaw->co[l] + C1*egy->Dlay->co[l]*egy->Temp->co[l] - C0*egy->Dlay->co[l]*egy->T0->co[l] ) / Dt;
		
		if(l<=ns+1) egy->Fenergy->co[l] -= egy->SWlayer->co[l];
	}
		
	//top boundary condition
	egy->Fenergy->co[par->nsurface] -= ( (1.-KNe)*EB + KNe*EB0 + egy->SWlayer->co[0] );
	//fictitious thermal capacity added to the skin layer (correspondent heat flux storage added to 1st layer as heat source)
	/*if(par->nsurface == 0 && ns+ng > 0){
		egy->Fenergy->co[par->nsurface] += (theta_snow(alpha_skin, 1.E-4, egy->Temp->co[par->nsurface]) - theta_snow(alpha_skin, 0., egy->T0->co[par->nsurface])) / Dt;
		egy->Fenergy->co[1] -= (theta_snow(alpha_skin, 1.E-4, egy->Temp->co[par->nsurface]) - theta_snow(alpha_skin, 0., egy->T0->co[par->nsurface])) / Dt;
	}*/
		
	//bottom boundary condition (treated as sink)
	kbb0 = k;
	kbb1 = k;
	egy->Fenergy->co[n] -= ( (par->Tboundary-egy->Temp->co[n])*(1.-KNe)*kbb1 + (par->Tboundary-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlay->co[n]/2.+par->Zboundary);
	egy->Fenergy->co[n] -= par->Fboundary;
	
	//include conduction in F(x0)
	update_F_energy(par->nsurface, n, egy->Fenergy, 1.-KNe, egy->Kth1, egy->Temp->co);
	update_F_energy(par->nsurface, n, egy->Fenergy, KNe, egy->Kth0, egy->T0->co);	
	
	//calculate the norm
	res = norm_2(egy->Fenergy, par->nsurface, n);
	//printf("res00:%e\n",res);

		
	do{
		
		cont++;
		
		for(l=par->nsurface;l<=n;l++){
			egy->T1->co[l] = egy->Temp->co[l];
		}
		
		//Calculates F'(x) referred to as Jacobian matrix
		
		//First the diagonal part is calculated and put into the vector dFenergy
		
		for(l=par->nsurface;l<=n;l++){
			egy->dFenergy->co[l] = 0.0;
		}
		
		//Heat capacity part
		for(l=1;l<=n;l++){
				
			//real thermal capacity
			C1 = calc_C(l, r, c, ns+ng, 1.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
			if(l<=ns+ng && fabs(egy->deltaw->co[l]-egy->wice->co[l])<1.E-7) C1 = Csnow_at_T_greater_than_0;
			//adds apparent thermal conductivity (due to phase change) for both snow and soil		
			if(l<=ns+ng){//snow
				C1 += Lf*(egy->wice->co[l]+egy->wliq->co[l])*dtheta_snow(par->alpha_snow, 1., egy->Temp->co[l])/egy->Dlay->co[l];
			}else{	//soil
				m=l-ns-ng;	
				if(egy->Temp->co[l]<=egy->Tstar->co[m]) C1 += rho_w*Lf*(Lf/(g*tk)*1.E3)*dteta_dpsi(Psif(egy->Temp->co[l]), 0.0, sl->pa->co[sy][jsat][m], 
											sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 
											1-1/sl->pa->co[sy][jns][m], PsiMin, 0.0);
			}
			
			egy->dFenergy->co[l] += C1*egy->Dlay->co[l] / Dt;
			
		}
		
		//Upper Boundary Condition
		egy->dFenergy->co[par->nsurface] -= (1.-KNe)*dEB_dT;
			
		//Bottom Boundary Condition
		egy->dFenergy->co[n] += (1.-KNe)*kbb1 / (egy->Dlay->co[n]/2.+par->Zboundary);
							
		//Update the part of Jacobian due to conduction is included in Kth1 and Kth0
		update_diag_dF_energy(par->nsurface, n, egy->dFenergy, 1.-KNe, egy->Kth1);
		
		//Extradiagonal part of the Jacobian
		for(l=par->nsurface;l<=n-1;l++){
			egy->udFenergy->co[l] = (1.-KNe)*egy->Kth1->co[l];
		}
		
		//fictitious thermal capacity added to the skin layer (correspondent heat flux storage added to 1st layer as heat source)
		/*if(par->nsurface == 0 && ns+ng > 0){
			egy->dFenergy->co[par->nsurface] += dtheta_snow(alpha_skin, 1.E-4, egy->Temp->co[par->nsurface]) / Dt;
			egy->udFenergy->co[1] -= dtheta_snow(alpha_skin, 1.E-4, egy->Temp->co[par->nsurface]) / Dt;
		}*/
		
		//Solve tridiagonal system
		tridiag2(1, r, c, par->nsurface, n, egy->udFenergy, egy->dFenergy, egy->udFenergy, egy->Fenergy, egy->Newton_dir);
						
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
												
			for(l=par->nsurface;l<=n;l++){ 
	
				egy->Temp->co[l] = egy->T1->co[l] + lambda[0] * egy->Newton_dir->co[l];
								
				if(egy->Temp->co[l] != egy->Temp->co[l]) printf("T no_value l:%ld r:%ld c:%ld cont:%ld\n",l,r,c,cont); 
				
				//soil
				if(l > ns+ng){
					m=l-ns-ng;
								
					if (i<=par->total_channel) {
						th0 = cnet->th->co[m][i];
						th1 = teta_psi(Psif(Fmin(egy->Tstar->co[m],egy->Temp->co[l])), 0.0, sl->pa->co[sy][jsat][m], 
									   sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 
									   1-1/sl->pa->co[sy][jns][m], PsiMin, sl->pa->co[sy][jss][m]);
						if (th1 > cnet->th->co[m][i] + cnet->thice->co[m][i]) th1 = cnet->th->co[m][i] + cnet->thice->co[m][i];
					}else {
						th0 = sl->th->co[m][r][c];
						th1 = teta_psi(Psif(Fmin(egy->Tstar->co[m],egy->Temp->co[l])), 0.0, sl->pa->co[sy][jsat][m], 
									   sl->pa->co[sy][jres][m], sl->pa->co[sy][ja][m], sl->pa->co[sy][jns][m], 
									   1-1/sl->pa->co[sy][jns][m], PsiMin, sl->pa->co[sy][jss][m]);
						if (th1 > sl->th->co[m][r][c] + sl->thice->co[m][r][c]) th1 = sl->th->co[m][r][c] + sl->thice->co[m][r][c];
					}

					egy->deltaw->co[l]=(th1-th0)*egy->Dlay->co[l]*rho_w;
									
				//snow	
				}else if(l>0) {
					
					th0=egy->wliq->co[l]/(egy->wice->co[l]+egy->wliq->co[l]);
					th1=theta_snow(par->alpha_snow, 1., egy->Temp->co[l]);
					
					egy->deltaw->co[l]=(th1-th0)*(egy->wice->co[l]+egy->wliq->co[l]);

				}
												
			}
							
			Tg=egy->Temp->co[par->nsurface];
			if(Tg<Tmin_surface_below_which_surfenergy_balance_recalculated) flagTmin++;

			if(cont < num_iter_after_which_surfenergy_balance_not_recalculated || flagTmin>0){
			
				//update egy->THETA taking into account evaporation (if there is not snow)
				if(ns+ng == 0){
					
					for(l=ns+ng+1;l<=n;l++){
						
						m=l-ns-ng;
						
						if (i<=par->total_channel) {
							egy->THETA->co[m] = cnet->th->co[m][i] + egy->deltaw->co[l]/(rho_w*egy->Dlay->co[l]);
						}else {
							egy->THETA->co[m] = sl->th->co[m][r][c] + egy->deltaw->co[l]/(rho_w*egy->Dlay->co[l]);
						}

						//add canopy transpiration
						if(egy->THETA->co[m] > sl->pa->co[sy][jres][1] + 1.E-3 && l <= egy->soil_transp_layer->nh ){
							egy->THETA->co[m] -= Fmax( Dt*fc*egy->soil_transp_layer->co[m]/(rho_w*egy->Dlay->co[l]), 0.0 );
							if(egy->THETA->co[m] < sl->pa->co[sy][jres][m]+1.E-3) egy->THETA->co[m] = sl->pa->co[sy][jres][m]+1.E-3;
						}
				
						//add soil evaporation
						if(egy->THETA->co[m] > sl->pa->co[sy][jres][1] + 1.E-3 && l <= egy->soil_evap_layer_bare->nh ){
							egy->THETA->co[m] -= Fmax( Dt*(1.-fc)*egy->soil_evap_layer_bare->co[m]/(rho_w*egy->Dlay->co[l]), 0.0 );
							egy->THETA->co[m] -= Fmax( Dt*fc*egy->soil_evap_layer_veg->co[m]/(rho_w*egy->Dlay->co[l]), 0.0 );
							if(egy->THETA->co[m] < sl->pa->co[sy][jres][m]+1.E-3) egy->THETA->co[m] = sl->pa->co[sy][jres][m]+1.E-3;
						}				
					}
				}
			
				//surface energy balance
				EnergyFluxes_no_rec_turbulence(f, Tg, r, c, ns+ng, egy->T0->co[1], Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, 
						 P, LR, sl->P->co[0][r][c], eps, fc, LSAI, decaycoeff0, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, egy->THETA->co,
						 sl->pa->co[sy], land->ty->co[lu], land->root_fraction->co[lu], par, egy->soil_transp_layer, SWin,
						 LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans, &(sl->Tv->co[r][c]), Qv, Ts, Qs, Hg0, Hg1, 
						 Eg0, Eg1, Lobukhov, rh, rv, rc, rb, ruc, &rh_g, &rv_g, Qg, u_top, decay, Locc, LWup_above_v, egy->Temp->co, 
						 egy->soil_evap_layer_bare, egy->soil_evap_layer_bare, flagTmin, cont);
				
				EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
				dEB_dT = -eps*dSB_dT(Tg) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;
																				
			}
			
			//F(T) = diag(egy->Fenergy(T)) + K(T)*T
			for(l=par->nsurface;l<=n;l++){
				egy->Fenergy->co[l] = 0.0;
			}
			
			for(l=1;l<=n;l++){
				
				//conductive M-matrix (lower diagonal part of this M-matrix is stored in a vector)
				if(l>1) kn = k;
				k = calc_k(l, r, c, ns, ng, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Temp->co, egy->Dlay->co, (*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par);		
				if(l>1){
					egy->Kth1->co[l-1] = -k * kn / ( k * 0.5 * egy->Dlay->co[l-1] + kn * 0.5 * egy->Dlay->co[l] );	// (-k/dz)
				}else if (l>par->nsurface) {
					egy->Kth1->co[l-1] = -k / ( 0.5 * egy->Dlay->co[l] );
				}
				
				//diagonal part of F due to heat capacity and boundary condition (except conduction)
				C0 = calc_C(l, r, c, ns+ng, 0.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
				C1 = calc_C(l, r, c, ns+ng, 1.0, egy->wice->co, egy->wliq->co, egy->deltaw->co, egy->Dlay->co, egy->Temp->co, sl);
				if(l<=ns+ng && fabs(egy->deltaw->co[l]-egy->wice->co[l])<1.E-7) C1 = Csnow_at_T_greater_than_0;
				egy->Fenergy->co[l] += ( Lf*egy->deltaw->co[l] + C1*egy->Dlay->co[l]*egy->Temp->co[l] - C0*egy->Dlay->co[l]*egy->T0->co[l] ) / Dt;
				if(l<=ns+1) egy->Fenergy->co[l] -= egy->SWlayer->co[l];
				
			}
						
			//top boundary condition
			egy->Fenergy->co[par->nsurface] -= ( (1.-KNe)*EB + KNe*EB0 + egy->SWlayer->co[0] );
			//fictitious thermal capacity added to the skin layer (correspondent heat flux storage added to 1st layer as heat source)
			/*if(par->nsurface == 0 && ns+ng > 0){
				egy->Fenergy->co[par->nsurface] += (theta_snow(alpha_skin, 1.E-4, egy->Temp->co[par->nsurface]) - theta_snow(alpha_skin, 0., egy->T0->co[par->nsurface])) / Dt;
				egy->Fenergy->co[1] -= (theta_snow(alpha_skin, 1.E-4, egy->Temp->co[par->nsurface]) - theta_snow(alpha_skin, 0., egy->T0->co[par->nsurface])) / Dt;
			}*/
			
			//bottom boundary condition (treated as sink)
			kbb1 = k;
			egy->Fenergy->co[n] -= ( (par->Tboundary-egy->Temp->co[n])*(1.-KNe)*kbb1 + (par->Tboundary-egy->T0->co[n])*KNe*kbb0 ) / (egy->Dlay->co[n]/2.+par->Zboundary);
			egy->Fenergy->co[n] -= par->Fboundary;

			update_F_energy(par->nsurface, n, egy->Fenergy, 1.-KNe, egy->Kth1, egy->Temp->co);
			update_F_energy(par->nsurface, n, egy->Fenergy, KNe, egy->Kth0, egy->T0->co);	
			res = norm_2(egy->Fenergy, par->nsurface, n);
			
			if(res <= res_av*(1.0 - ni*lambda[0])) iter_close2=1;
			if(lambda[0] <= par->min_lambda_en) cont_lambda_min++;
			if(cont_lambda_min > par->max_times_min_lambda_en){
				if(par->exit_lambda_min_en == 1){
					return 0;
				}else {
					iter_close2=1;
					cont_lambda_min=0;
				}
			}
			
			//printf("res:%e lamnda:%e cont:%ld Dt:%f\n",res,lambda[0],cont,Dt);
				   
		}while(iter_close2!=1);
						
		if(iter_close>=0 && res<=par->tol_energy) iter_close=1;
		if(cont>=par->maxiter_energy) iter_close=1;		
										
	}while(iter_close!=1);
		
	//if there is no convergence, go out of the loop
	if(res > par->tol_energy) return 0;	
	
	//checks
	l = par->nsurface;
	//error messages
	if(egy->Temp->co[l]!=egy->Temp->co[l]){
		printf("T no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
		fprintf(f,"T no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
	}
	if(egy->Temp->co[l]<Tmin || egy->Temp->co[l]>Tmax){
		printf("T outside of range, error 2, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f LW:%f H:%f LE:%f Ta:%f\n",l,r,c,egy->Temp->co[l],*LW,*H,latent(Tg,Levap(Tg))*(*E),Ta);
		fprintf(f,"T outside of range, error 2, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f LW:%f H:%f LE:%f Ta:%f\n",l,r,c,egy->Temp->co[l],*LW,*H,latent(Tg,Levap(Tg))*(*E),Ta);
	}
	
	for(l=par->nsurface;l<=Nl+ns+ng;l++){ 	

		//error messages
		if(egy->Temp->co[l]!=egy->Temp->co[l]){
			printf("T no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
			fprintf(f,"T no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
		}
		if(egy->Temp->co[l]<Tmin || egy->Temp->co[l]>Tmax){
			printf("T outside of range, error 2, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f LW:%f H:%f LE:%f Ta:%f\n",l,r,c,egy->Temp->co[l],*LW,*H,latent(Tg,Levap(Tg))*(*E),Ta);
			fprintf(f,"T outside of range, error 2, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f LW:%f H:%f LE:%f Ta:%f\n",l,r,c,egy->Temp->co[l],*LW,*H,latent(Tg,Levap(Tg))*(*E),Ta);
		}
		
		if(l>0){
			if(egy->deltaw->co[l]!=egy->deltaw->co[l]){
				printf("dw no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
				fprintf(f,"dw no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
			}
			
			//checks
			if(egy->wice->co[l]-egy->deltaw->co[l]<0) egy->deltaw->co[l]=egy->wice->co[l];
			if(egy->wliq->co[l]+egy->deltaw->co[l]<0) egy->deltaw->co[l]=-egy->wliq->co[l];
			
			//update
			egy->wliq->co[l]+=egy->deltaw->co[l];  
			egy->wice->co[l]-=egy->deltaw->co[l];
			if(egy->Temp->co[l]>0) egy->wice->co[l]=0.0;
						
		}
	}

	
	//update water on the canopy	
	*Wcrn = *Wcrn + dWcrn;
	*Wcsn = *Wcsn + dWcsn;
	
	return 1;
		
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_soil_land(long nsurf, long n, long r, long c, long sy, double fc, double Dt, ENERGY *egy, SOIL *sl, DOUBLETENSOR *T, DOUBLETENSOR *P, 
					  DOUBLETENSOR *thliq, DOUBLETENSOR *thice, DOUBLETENSOR *ET, DOUBLEMATRIX *Tskin){

	long l, j;
	double th_oversat, psisat;
	
	//skin temperature
	l = nsurf;
	Tskin->co[r][c] = egy->Temp->co[l];

	//soil variables
	for (l=1; l<=Nl; l++) {

		j = l+n;
		
		//canopy transpiration
		if(l <= egy->soil_transp_layer->nh) ET->co[l][r][c] += fc*egy->soil_transp_layer->co[l]*Dt;
		
		//soil evaporation
		if(l <= egy->soil_evap_layer_bare->nh){
			ET->co[l][r][c] += (1.-fc)*egy->soil_evap_layer_bare->co[l]*Dt;
			ET->co[l][r][c] += fc*egy->soil_evap_layer_veg->co[l]*Dt;
		}
		
		//water pressure and contents
		psisat = psi_saturation(thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l]);
		th_oversat = Fmax( P->co[l][r][c] - psisat , 0.0 ) * sl->pa->co[sy][jss][l];
		
		thliq->co[l][r][c] = egy->wliq->co[j]/(rho_w*egy->Dlay->co[j]);
		thice->co[l][r][c] = egy->wice->co[j]/(rho_w*egy->Dlay->co[j]);
		
		P->co[l][r][c] = psi_teta(thliq->co[l][r][c] + th_oversat, thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
								  sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
				
		//temperature
		T->co[l][r][c] = egy->Temp->co[j];

	}
} 	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_soil_channel(long nsurf, long n, long ch, long r, long c, long sy, double fc, double Dt, ENERGY *egy, SOIL *sl, DOUBLEMATRIX *T, DOUBLEMATRIX *P, 
					  DOUBLEMATRIX *thliq, DOUBLEMATRIX *thice, DOUBLEMATRIX *ET, DOUBLEVECTOR *Tskin){
	
	long l, j;
	double th_oversat, psisat;
	
	//skin temperature
	l = nsurf;
	Tskin->co[ch] = egy->Temp->co[l];
	
	//soil variables
	for (l=1; l<=Nl; l++) {
		
		j = l+n;
		
		//canopy transpiration
		if(l <= egy->soil_transp_layer->nh) ET->co[l][ch] += fc*egy->soil_transp_layer->co[l]*Dt;
		
		//soil evaporation
		if(l <= egy->soil_evap_layer_bare->nh){
			ET->co[l][ch] += (1.-fc)*egy->soil_evap_layer_bare->co[l]*Dt;
			ET->co[l][ch] += fc*egy->soil_evap_layer_veg->co[l]*Dt;
		}
		
		//water pressure and contents
		psisat = psi_saturation(thice->co[l][ch], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
								sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l]);
		th_oversat = Fmax( P->co[l][ch] - psisat , 0.0 ) * sl->pa->co[sy][jss][l];
		
		thliq->co[l][ch] = egy->wliq->co[j]/(rho_w*egy->Dlay->co[j]);
		thice->co[l][ch] = egy->wice->co[j]/(rho_w*egy->Dlay->co[j]);
		
		P->co[l][ch] = psi_teta(thliq->co[l][ch] + th_oversat, thice->co[l][ch], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
								  sl->pa->co[sy][jns][l], 1.-1./sl->pa->co[sy][jns][l], PsiMin, sl->pa->co[sy][jss][l]);
		
		//temperature
		T->co[l][ch] = egy->Temp->co[j];
		
	}
} 

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_F_energy(long nbeg, long nend, DOUBLEVECTOR *F, double w, DOUBLEVECTOR *K, double *T){
	
	long l;
	
	for(l=nbeg;l<=nend;l++){
		
		if(l==nbeg){
			F->co[l] += w*(-K->co[l]*T[l] + K->co[l]*T[l+1]); 
		}else if(l<nend){
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

void update_diag_dF_energy(long nbeg, long nend, DOUBLEVECTOR *dF, double w, DOUBLEVECTOR *K){
	
	long l;
	
	for(l=nbeg;l<=nend;l++){
		
		if(l==nbeg){
			dF->co[l] -= w*K->co[l]; 
		}else if(l<nend){
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

void EnergyFluxes(FILE *f, double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
		double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, 
		double P, double LR, double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, 
		double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn, double *dWcsn, double *theta, double **soil, 
		double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer, double SWin, double LWin, double SWv, double *LW, 
		double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv, double *LEv, double *Etrans, 
		double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Lobukhov, 
		double *rh, double *rv, double *rc, double *rb, double *ruc, double *rh_g, double *rv_g, double *Qg, 
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
	SpecHumidity_2(Qg, &dQgdT, 1.0, Tg, P);
	
	//initialization
	if(fc==0){ 
				
		*Qs=(double)number_novalue;  
		*Ts=(double)number_novalue;  
		*Qv=(double)number_novalue;  
		*Tv=(double)number_novalue;  
		*u_top=(double)number_novalue;  
		*rh=1.E99; 
		*rv=1.E99; 
		*rc=1.E99; 
		*rb=1.E99; 
		*ruc=1.E99;
		
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
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");
		}
				
		if(Eg!=Eg){
			printf("Eg no value bare soil %ld %ld \n",r,c);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Eg no value bare soil %ld %ld \n",r,c);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");			
		}

	}
	
	//vegetation
	if(fc>0){
				
		Tcanopy(r, c, Tv0, Tg, *Qg, dQgdT, Tg0, Qg0, Ta, Qa, zmu, zmT, z0v, z0s, d0v, rz0v, hveg, v, LR, P, SWin, SWv, LWin, 
				e, LSAI, decaycoeff0, land, Wcrn, Wcrnmax, Wcsn, Wcsnmax, dWcrn, dWcsn, LWv, &LWg, Hv, &Hg, &dHg_dT, LEv, &Eg, 
				&dEg_dT, Ts, Qs, root, theta, soil_transp_layer, Lobukhov, par, n, &rm, rh, rv, rc, rb, ruc, u_top, Etrans, Tv, Qv, 
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
			printf("Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");
		}
		
		if(Eg!=Eg){
			printf("Eg no value vegetated soil %ld %ld \n",r,c);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Eg no value vegetated soil %ld %ld \n",r,c);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");			
		}
	}	
			
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void EnergyFluxes_no_rec_turbulence(FILE *f, double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, 
		double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, 
		double P, double LR, double psi, double e, double fc, double LSAI, double decaycoeff0, double Wcrn, 
		double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn, double *dWcsn, double *theta, double **soil, 
		double *land, double *root, PAR *par, DOUBLEVECTOR *soil_transp_layer, double SWin, double LWin, double SWv, double *LW, 
		double *H, double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv, double *LEv, double *Etrans, 
		double *Tv, double *Qv, double *Ts, double *Qs, double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *Lobukhov, 
		double *rh, double *rv, double *rc, double *rb, double *ruc, double *rh_g, double *rv_g, double *Qg, 
		double *u_top, double *decay, double *Locc, double *LWup_above_v, double *T, DOUBLEVECTOR *soil_evap_layer_bare,
		DOUBLEVECTOR *soil_evap_layer_veg, short flagTmin, long cont){

	
	double Hg, dHg_dT, Eg, dEg_dT, LWg, LWup;
	double dQgdT;
	double dLWvdT,LWv_old;
	double alpha, beta;	//from Ye and Pielke, 1993
	double rm=1.E20;

	//initalization		
	*H=0.0;	
	*E=0.0;	
	
	*dH_dT=0.0;	
	*dE_dT=0.0;	
	
	*LW=0.0;
	*LWup_above_v=0.0;
			
	//thermodynamical calculations
	SpecHumidity_2(Qg, &dQgdT, 1.0, Tg, P);
		
	if(fc<1){		

		if(flagTmin>0 && cont>num_iter_after_which_only_neutrality) aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, LR, Lobukhov, &rm, rh_g, rv_g, par->state_turb, 2, par->maxiter_Businger);
		
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
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Hg no value bare soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");
		}
		
		if(Eg!=Eg){
			printf("Eg no value bare soil %ld %ld \n",r,c);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Eg no value bare soil %ld %ld \n",r,c);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");			
		}
		
	}
	
	//TURBULENT FLUXES FOR CANOPY CANOPY
	if(fc>0){

		find_actual_evaporation_parameters(r,c,&alpha, &beta, soil_evap_layer_veg, theta, soil, T, psi, P, *ruc, Ta, Qa, *Qg, n);
		*Ts=(Ta/(*rh)+Tg/(*ruc)+(*Tv)/(*rb))/(1./(*rh)+1./(*ruc)+1./(*rb));
		*Qs=(Qa/(*rv)+(*Qg)*alpha*beta/(*ruc)+(*Qv)/(*rc))/(1./(*rv)+beta/(*ruc)+1./(*rc));
		turbulent_fluxes(*ruc, *ruc/beta, P, *Ts, Tg, *Qs, *Qg*alpha, dQgdT*alpha, &Hg, &dHg_dT, &Eg, &dEg_dT);
				
		LWv_old=*LWv;
		longwave_vegetation(LWin, e, Tg, *Tv, LSAI, LWv, &LWg, &dLWvdT, &LWup);
		*LWv=LWv_old;
				
		*H+=fc*Hg;	*dH_dT+=fc*dHg_dT;
		*E+=fc*Eg;	*dE_dT+=fc*dEg_dT;		
				
		*LW+=fc*LWg;
		*LWup_above_v+=fc*LWup;
				
		*Hg1=Hg;
		*Eg1=Eg;
		
		if(Hg!=Hg){
			printf("Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Hg no value vegetated soil EnergyFluxes %ld %ld rh_g:%e rv_g:%e Hg:%f Eg:%f\n",r,c,*rh_g,*rv_g,Hg,Eg);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");
		}
		
		if(Eg!=Eg){
			printf("Eg no value vegetated soil %ld %ld \n",r,c);
			printf("Ta:%f Tg:%f\n",Ta,Tg);
			printf("zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fprintf(f,"Eg no value vegetated soil %ld %ld \n",r,c);
			fprintf(f,"Ta:%f Tg:%f\n",Ta,Tg);
			fprintf(f,"zmu:%f zmT:%f z0s:%f d0s:%f rz0s:%f v:%f Ta:%f Tg:%f Qa:%f Qg:%e P:%f rm:%f psi:%f \n",zmu,zmT,z0s,d0s,rz0s,v,Ta,Tg,Qa,*Qg,P,rm,psi);
			fclose(f);		
			t_error("Not possible to continue");			
		}
		
	}	
		
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

	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


double flux(long i, long icol, double **met, double k, double est){
	
	double F=k*met[i-1][icol];
								
	if((long)F==number_absent || (long)F==number_novalue) F=est;
		
	return(F);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_map_plots(long n, long r, long c, double W, PAR *par, ENERGY *egy, METEO *met, SNOW *snow, double Hg, double LEg, double Hv, double LEv,
					  double SWin, double SWg, double SWv, double LWin, double LWg, double LWv, double Ts, double Tg, double Tv){
		
	if(W>0){

		egy->Hgplot->co[r][c]+=Hg*W;
		
		egy->LEgplot->co[r][c]+=LEg*W;
		egy->Hvplot->co[r][c]+=Hv*W;
		egy->LEvplot->co[r][c]+=LEv*W;
		
		egy->SWinplot->co[r][c]+=SWin*W;
		egy->SWgplot->co[r][c]+=SWg*W;
		egy->SWvplot->co[r][c]+=SWv*W;				
		
		egy->LWinplot->co[r][c]+=LWin*W;
		egy->LWgplot->co[r][c]+=LWg*W;	
		egy->LWvplot->co[r][c]+=LWv*W;	
		
		egy->Tsplot->co[r][c]+=Ts*W;
		egy->Tgplot->co[r][c]+=Tg*W;
		egy->Tvplot->co[r][c]+=Tv*W;
		
		snow->Dplot->co[r][c]+=DEPTH(r,c,snow->S->lnum,snow->S->Dzl)*W;
		met->Taplot->co[r][c]+=met->Tgrid->co[n][r][c]*W;
		
		met->Vspdplot->co[r][c]+=met->Vgrid->co[n][r][c]*W;
		met->Vdirplot->co[r][c]+=met->Vdir->co[n][r][c]*W;
		met->RHplot->co[r][c]+=met->RHgrid->co[n][r][c]*W;
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

void prepare_output(FILE *f, double Dt, double evap_snow, double melt_snow, double evap_glac, double melt_glac, double prec_rain, double prec_snow,
					ENERGY *egy, WATER *wat, SNOW *snow, GLACIER *glac, LAND *land, TOPO *top, SOIL *sl, METEO *met, TIMES *times, PAR *par, 
					long r, long c, double LE, double surfEB, double H, double G, double Ts, double SWin, double SW, double SWbeam, double eps, 
					double LWin, double LW, double cosinc){
	
	
	//Volumes [mm]
	if(par->output_snow>0){
		snow->MELTED->co[r][c] += melt_snow;
		snow->SUBL->co[r][c] += evap_snow;
		if(snow->S->type->co[r][c]==2) snow->t_snow->co[r][c] += Dt/3600.0;
	}
	if(par->output_glac>0 && par->glaclayer_max>0){
		glac->MELTED->co[r][c] += melt_glac;
		glac->SUBL->co[r][c] += evap_glac;
	}
	if(par->output_meteo>0){
		wat->PrTOT_mean->co[r][c] += prec_rain;		
		wat->PrSNW_mean->co[r][c] += prec_snow;		
	}
	
	
	if(par->output_surfenergy>0){
		egy->ET_mean->co[r][c] += LE/(par->output_surfenergy*3600.0/Dt); //[W/m^2]
		egy->SEB_mean->co[r][c] += surfEB/(par->output_surfenergy*3600.0/Dt); //[W/m^2]
		egy->G_snowsoil->co[r][c] += G/(par->output_surfenergy*3600.0/Dt); //[W/m^2]		
		egy->H_mean->co[r][c] += H/(par->output_surfenergy*3600.0/Dt); //[W/m^2]
		egy->Ts_mean->co[r][c] += Ts/(par->output_surfenergy*3600.0/Dt); 
		egy->Rswdown_mean->co[r][c] += SWin/(par->output_surfenergy*3600.0/Dt);
		egy->Rswbeam_mean->co[r][c] += SWbeam/(par->output_surfenergy*3600.0/Dt);
		egy->Rn_mean->co[r][c] += (SW+LW)/(par->output_surfenergy*3600.0/Dt); //[W/m^2]
		egy->LWin_mean->co[r][c] += LWin/(par->output_surfenergy*3600.0/Dt);
		egy->LW_mean->co[r][c] += LW/(par->output_surfenergy*3600.0/Dt);
		egy->SW_mean->co[r][c] += SW/(par->output_surfenergy*3600.0/Dt);
		if(par->distr_stat==1){
			if(egy->ET_max->co[r][c]<LE) egy->ET_max->co[r][c]=LE;
			if(egy->ET_min->co[r][c]>LE) egy->ET_min->co[r][c]=LE;
			if(egy->G_max->co[r][c]<surfEB) egy->G_max->co[r][c]=surfEB;
			if(egy->G_min->co[r][c]>surfEB) egy->G_min->co[r][c]=surfEB;
			if(egy->H_max->co[r][c]<H) egy->H_max->co[r][c]=H;
			if(egy->H_min->co[r][c]>H) egy->H_min->co[r][c]=H;
			if(egy->Ts_max->co[r][c]<Ts) egy->Ts_max->co[r][c]=Ts;
			if(egy->Ts_min->co[r][c]>Ts) egy->Ts_min->co[r][c]=Ts;
			if(egy->Rswdown_max->co[r][c]<SWin) egy->Rswdown_max->co[r][c]=SWin;
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

void output_pixel(long n, double Dt, double Dt_plot, long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double prec_snow_atm, 
				  double prec_rain_atm, double evap_soil, double melt_snow, double evap_snow, double melt_glac, double evap_glac, double transpiration,
				  double LE, double H, double surfEB, double G, double Tg, double eps, double LSAI, double LWin, double SWin, double LW, double SW, 
				  double epsa, double epsa_min, double epsa_max, double SWbeam, double SWdiff, double Tdew, double Lobukhov, PAR *par, 
				  WATER *wat, ENERGY *egy, TOPO *top, METEO *met, SNOW *snow, GLACIER *glac, LAND *land, SOIL *sl, double Tpoint, double Ppoint, 
				  double Vpoint, double RHpoint, double z0soil, double z0, double d0, double SWv, double LWv, double Tv, double Ts, double Hg0, 
				  double Hg1, double Eg0, double Eg1, double fc, double rh, double rv, double rb, double rc, double ruc, double Hv, double LEv, 
				  double Qv, double Qg, double Qa, double Qs, double u_top, double decay, double Locc, double SWup_above_v, double LWup_above_v){

	long i;
	
	if(par->state_pixel ==1 && par->Dtplot_point->co[i_sim] > 1.E-5){
		for(i=0;i<par->rc->nrh;i++){
			if(r==par->rc->co[i+1][1] && c==par->rc->co[i+1][2]){
			
				//All volumes [mm]
				outdata_point[oprecsnow][i] += prec_snow;
				outdata_point[oprecrain][i] += (prec_rain_on_soil+prec_rain_on_snow);
				outdata_point[orainonsnow][i] += prec_rain_on_snow;	
				outdata_point[osnowover][i] += prec_snow_atm;
				outdata_point[orainover][i] += prec_rain_atm;
				outdata_point[oevapsur][i] += evap_soil;
				outdata_point[otrasp][i] += transpiration;
				outdata_point[omrsnow][i] += melt_snow;
				outdata_point[osrsnow][i] += evap_snow;
				outdata_point[omrglac][i] += melt_glac;
				outdata_point[osrglac][i] += evap_glac;
				
				outdata_point[oLE][i] += LE*(Dt/Dt_plot); //ET[W/m^2]
				outdata_point[oH][i] += H*(Dt/Dt_plot); //H[W/m^2]
				outdata_point[oEB][i] += surfEB*(Dt/Dt_plot); 
				outdata_point[oG][i] += G*(Dt/Dt_plot); 
				outdata_point[oTg][i] += Tg*(Dt/Dt_plot); //Ts[C]
				outdata_point[oSWin][i] += SWin*(Dt/Dt_plot);
				outdata_point[oLWin][i] += LWin*(Dt/Dt_plot);
				outdata_point[oSW][i] += SW*(Dt/Dt_plot);
				outdata_point[oLW][i] += LW*(Dt/Dt_plot);
				outdata_point[oV][i] += Vpoint*(Dt/Dt_plot);
				outdata_point[oRH][i] += RHpoint*(Dt/Dt_plot);
				outdata_point[oPa][i] += Ppoint*(Dt/Dt_plot);
				outdata_point[oTa][i] += Tpoint*(Dt/Dt_plot);
				outdata_point[oTdew][i] += Tdew*(Dt/Dt_plot);
				outdata_point[oLobuk][i] += Lobukhov*(Dt/Dt_plot);
				outdata_point[oSWb][i] += SWbeam*(Dt/Dt_plot);
				outdata_point[oSWd][i] += SWdiff*(Dt/Dt_plot);
				outdata_point[oSWv][i] += SWv*(Dt/Dt_plot);
				outdata_point[oLWv][i] += LWv*(Dt/Dt_plot);
				outdata_point[oHv][i] += Hv*(Dt/Dt_plot);
				outdata_point[oLEv][i] += LEv*(Dt/Dt_plot);
				outdata_point[oHg0][i] += Hg0*(Dt/Dt_plot);
				outdata_point[oLEg0][i] += Levap(Tg)*Eg0*(Dt/Dt_plot);
				outdata_point[oHg1][i] += Hg1*(Dt/Dt_plot);
				outdata_point[oLEg1][i] += Levap(Tg)*Eg1*(Dt/Dt_plot);
				outdata_point[oSWup][i] += SWup_above_v*(Dt/Dt_plot);
				outdata_point[oLWup][i] += LWup_above_v*(Dt/Dt_plot);
				outdata_point[oEcan][i] += (SWv+LWv-Hv-LEv)*(Dt/Dt_plot);
				outdata_point[oHup][i] += (H+fc*Hv)*(Dt/Dt_plot);
				outdata_point[oLEup][i] += (LE+fc*LEv)*(Dt/Dt_plot);
				
				if(epsa_min>0){
					outdata_point[ominLWin][i] += (epsa_min*5.67E-8*pow(Tpoint+tk,4.0))*(Dt/Dt_plot);
				}else{
					outdata_point[ominLWin][i] = (double)number_novalue;
				}
				
				if(epsa_max>0){
					outdata_point[omaxLWin][i] += (epsa_max*5.67E-8*pow(Tpoint+tk,4.0))*(Dt/Dt_plot);
				}else{
					outdata_point[omaxLWin][i] = (double)number_novalue;
				}
				
				outdata_point[oVdir][i] += (met->Vdir->co[n][r][c])*(Dt/Dt_plot);
				
				outdata_point[oLSAI][i] += LSAI*(Dt/Dt_plot);
				outdata_point[oz0v][i] += z0*(Dt/Dt_plot);
				outdata_point[od0v][i] += d0*(Dt/Dt_plot);				
				outdata_point[oTv][i] += Tv*(Dt/Dt_plot);
				outdata_point[oTs][i] += Ts*(Dt/Dt_plot);
				outdata_point[ofc][i] += fc*(Dt/Dt_plot);
				outdata_point[oQv][i] += (Qv)*(Dt/Dt_plot);
				outdata_point[oQg][i] += (Qg)*(Dt/Dt_plot);
				outdata_point[oQa][i] += (Qa)*(Dt/Dt_plot);
				outdata_point[oQs][i] += (Qs)*(Dt/Dt_plot);
				outdata_point[outop][i] += (u_top)*(Dt/Dt_plot);
				outdata_point[odecay][i] += (decay)*(Dt/Dt_plot);
				outdata_point[oLobukcan][i] += (Locc)*(Dt/Dt_plot);
				
				outdata_point[owcan_rain][i] += wat->wcan_rain->co[r][c]*(Dt/Dt_plot);
				outdata_point[owcan_snow][i] += wat->wcan_snow->co[r][c]*(Dt/Dt_plot);			
				
			}
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_basin(double Dt, double Dt_plot, double number_of_pixels, double prec_rain, double prec_snow, double prec_rain_atm, 
				  double prec_snow_atm, double Ta, double Tg, double Tv, double Eg, double Evt, double LE, double H, double SW, double LW, 
				  double LEv, double Hv, double SWv, double LWv, double SWin, double LWin){

	
	outdata_basin[ooprecrain] += prec_rain/(double)number_of_pixels;										
	outdata_basin[ooprecsnow] += prec_snow/(double)number_of_pixels;			
	outdata_basin[oorainover] += prec_rain_atm/(double)number_of_pixels;													
	outdata_basin[oosnowover] += prec_snow_atm/(double)number_of_pixels;
	outdata_basin[ooevapsur] += Eg/(double)number_of_pixels;															
	outdata_basin[ootrasp] += Evt/(double)number_of_pixels;														
	
	outdata_basin[ooTa] += Ta*(Dt/Dt_plot)/(double)number_of_pixels;															
	outdata_basin[ooTg] += Tg*(Dt/Dt_plot)/(double)number_of_pixels;	
	outdata_basin[ooTv] += Tv*(Dt/Dt_plot)/(double)number_of_pixels;			
	outdata_basin[ooLE] += LE*(Dt/Dt_plot)/(double)number_of_pixels;																
	outdata_basin[ooH] += H*(Dt/Dt_plot)/(double)number_of_pixels;																
	outdata_basin[ooSW] += SW*(Dt/Dt_plot)/(double)number_of_pixels;													
	outdata_basin[ooLW] += LW*(Dt/Dt_plot)/(double)number_of_pixels;	
	outdata_basin[ooLEv] += LEv*(Dt/Dt_plot)/(double)number_of_pixels;	
	outdata_basin[ooHv] += Hv*(Dt/Dt_plot)/(double)number_of_pixels;	
	outdata_basin[ooSWv] += SWv*(Dt/Dt_plot)/(double)number_of_pixels;	
	outdata_basin[ooLWv] += LWv*(Dt/Dt_plot)/(double)number_of_pixels;	
	outdata_basin[ooSWin] += SWin*(Dt/Dt_plot)/(double)number_of_pixels;	
	outdata_basin[ooLWin] += LWin*(Dt/Dt_plot)/(double)number_of_pixels;		
	
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

				 

	
