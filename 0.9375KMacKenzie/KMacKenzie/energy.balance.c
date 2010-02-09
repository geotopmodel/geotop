
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Matteo Dall'Amico, Emanuele Cordano

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


//Author: Stefano Endrizzi
//Date: 3 October 2008
//Contents: Energy balance (and also mass balance for snow and glacier)
#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "times.h"
#include "meteo.09375.h"
#include "snow.09375.h"
#include "pedo.funct.h"
#include "vegetation.h"
#include "radiation.h"
#include "turbulence.h"
#include "util_math.h"
#include "micromet.h"
#include "output.09375.h"
#include "energy.balance.h"
#include "PBSM.h"
//#include "SnowTran.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern char *error_file_name;
extern long Nl; // total number of soil layers (constant in the whole basin)
extern long Nr; // total number of rows (of the map)
extern long Nc;// total number of columns (of the map)
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void energy_balance(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, SNOW *snow, GLACIER *glac)

	{

	//1.DICHIARAZIONE DELLE VARIABILI
	//*******************************

	long r;//Counters of rows of the basin
	long c;//Counters of columns of the basin
	long l;//Counters of layers of the basin
	long j;//Counter
	long Ns;// Maximum number of snow layers
	long Ng;// Maximum number of glacier layers
	long ns;//Number of snow layer in a pixel,
	long ng=0;// number of glacier layer in a pixel
	double E0;//soil-earth correction [-]
	double Ts;//GROUND SURFACE TEMPERATURE
	double Qs;// Specific humidity of the surface
	double snowD;// Snow depth [mm] in a pixel
	double glacD=0.0;// glacier depth [mm] in a pixel
	double Prain_over; //Precipitation as rain [mm] over the canopy
	double Psnow_over;//Precipitation as snow [mm water equivalent] over the canopy
	//double Psoil=0.0; //precipitation that reaches the soil [mm] below canopy /* modified by Emanuele Cordano on 24/9/9 */
	double Prain=0.0;// precipitation that reaches the soil as rain [mm] below canopy
	double Prain_on_snow; // precipitation that reaches the snow as rain [mm] below canopy
	double Psnow=0.0;	//precipitation that reaches the soil as snow [mm] below canopy

	double SWin; //Short wave radiation incoming [W/mq]
	double	SW;// SW net to the ground [W/m2]
	double SWbeam; // Incoming Direct SW radition [W/m2]
	double SWdiff; // Diffuse SW radition [W/m2]
	double SWv_vis;// SW going to vegetation nel visibile
	double	SWv_nir;// SW going to vegetation nel near infrared
	double	SWg_vis;// SW che va nel ground nel visibile
	double	SWg_nir;// SW going to the ground near infrared
	double	cosinc=0.0; // cosine of the incidence angle [-] for a pixel
	double	avis_b;// albedo nel visibile beam
	double	avis_d;// albedo nel visibile diffuse
	double	anir_b;// albedo near infrared beam
	double	anir_d;// albedo near infrared diffuse

	double epsa=0.0; //Atmospheric absorbance for long wave radiation [-]
	double	LWin; // incoming Long Wave long wave radiation diffused towards land [W/m2] for a pixel
	double	LW;// LW net to the ground
	double	LWv;// Long wave of vegetation
	double	Tsurr;// Temperature of the sorrounding surfaces (not of the pixel) [ûC]

	double Hv; // Sensible heat flux of the vegetation [W/m2]
	double H;//Sensible heat flux of the ground [W/m2] for a pixel

	double LE;// Latent heat flux
	double	LEv;
	double	E;//Evaporation fluxes [kg/(s*m2)]
	double	Etrans;// transpiration from the vegetation
	//double	fwet;
	double	fsnow;
	double fsnownew;
	double	fsnowcan;// % of canopy covered by snow

	double surfEB;//heat flux on the surface (be it snow or soil) [W/m2]
	double G;// Heat flux going into the soil surface [W/m2]

	DOUBLEVECTOR *theta;//Adimensional water content

	//DOUBLEVECTOR *c_heat; //Thermal capacity [J m^-3 K^-1] (for a pixel)
	DOUBLEVECTOR *k_thermal; //, thermal conductivity [W m^-1 K^-1] for a snow/soil layer (for a pixel)
	DOUBLEVECTOR *k_thermal_interface; // thermal conductivity [W m^-1 K^-1] at the interface between 2 layers (for a pixel)
	//The 1st component is for the highest layer (be it of snow or soil), the (nl+ns)th component for the deepest soil layer

	DOUBLEVECTOR *D; //vector of column (snow+soil) layer thickness [m]
	DOUBLEVECTOR *wliq; // vector of column (snow+soil)liquid water content [kg/mq]
	DOUBLEVECTOR *wice; // vector of column (snow+soil) ice content [kg/mq]
	DOUBLEVECTOR *Temp; // vector of column (snow+soil) temperature [Celsius]
	DOUBLEVECTOR *deltaw; // vector of column (snow+soil) melting ice (if>0) or icing water (if<0) [kg/mq] in Dt
	//The 1st component is for the highest layer (be it of snow or soil), the (nl+ns)th component for the deepest soil layer

	//Other auxilairy variables
	double z0;// surface roughness [mm] i.e. the height above the ground where the wind velocity profile is 0, following the rule of thumb it's roughly 1/10 of the height of the object.
	double z0veg;// surface roughness [mm] for the vegetation
	double d0; //zero-plane displacement [mm]. It is usually 2/3*z0. It only plays a role when the object are so dense that the wind will hardly move though them, so it's like as the ground was displaced up by d0.
	double d0veg;//zero-plane displacement for vegetation
	double z0_z0t;//ratio z0/z0T
	double hveg;// vegetation height (not present in this version)
	double eps;// surface LW emissivity (be it soil or snow) [-]
	double ee=0.0;// saturated water vapour pressure
	double dee=0.0;// derivative of saturated water pressure
	double Tdew=0.0; // dew temperature (temp. di rugiada)
	double Qa;// Air specific humidity
	double drip_rain;
	double max_wcan_rain;
	double max_wcan_snow;
	//double evap_c;
	//double drip_c;
	double drip_snow;
	//double maxstorage_c;
	double LAI=0.5;// leaf area index
	double fc=0.0;// canopy fraction
	double fc0;
	//double drip_sn; //unused variable
	//double	expon;// unused variable
	double epsa_min;
	double	epsa_max;
	double	tau_atm;// atmospheri transmittance
	double	RHpoint;// Relative humidity of the point [%]
	double	Vpoint;// wind velocity of the point [m/s]
	double Ppoint;
	double Tpoint;
	double Precpoint;
	DOUBLEVECTOR *turbulence;
	DOUBLEVECTOR *ftcl;
	DOUBLEVECTOR *SWabsf;
	DOUBLEVECTOR *froot;
	DOUBLEVECTOR *DPsi;
	double DTcorr=0.0;
	double Mr_snow;// snow melting rate [mm/s]
	double	Er_snow;// snow evaporation rate [mm/s]
	double	Sr_snow;// snow sublimation rate [mm/s]
	double	Mr_glac;// glacier melting rate [mm/s]
	double	Er_glac;// glacier evaporation rate [mm/s]
	double	Sr_glac;// glacier sublimation rate [mm/s]
	double	Er_soil;// soil evaporation rate [mm/s]
	double	Sr_soil;// soil sublimation rate [mm/s]
	double zmeas_T; // [m] elevation of the air temperature sensor (with respect to the ground surface)
	double	zmeas_u;// [m] elevation of the wind temperature sensor (with respect to the ground surface)
	double fcloud;// cloud fraction of the sky (from 0 to 1)
	double	tau_cloud=0.0;// clear sky trasmissivity
	double	tau_cloud_av=0.0;//average clear sky trasmissivity /* modified by Emanuele Cordano on 24/9/9 */
	double	sa;// sin(alpha) where alpha is the solar altitude angle
	double	Hadv=0.0;
	double psisat;
	double u_top;
	short sy; // soil type
	short lu;// land use

	double Hg0,
		Hg1,
		Eg0,
		Eg1;
	//	dE_dT,
	//	dH_dT;
	//long cont;
	double rh,// resistenza del calore sensibile
		rv,// resistenza del calore latente
		rc,// resistenza della canopy per calore latente
		rb,// resistenza della canopy per calore sensibile
		rh_ic,// resistenza del calore sensibile undrcanopy
		rv_ic,// resistenza del calore latente undercanopy
		Qv,// umiditˆ specifica vicina alla canopy
		Qg,//umiditˆ specifica vicina alla superficie del suolo
		decaycoeff,
		Locc;

	double LAIthres=0.01;

	short SWdata;
	//FILE *f;// file pointer

	//ALLOCATION
	Ns=par->snowlayer_max; /*maximum number of snow layers*/
	Ng=par->glaclayer_max; /*maximum number of glacier layers*/

	//The 1st component is for the highest layer (be it of snow or glacier or soil), the (Nl+ns)th component for the deepest soil layer
	k_thermal=new_doublevector(Nl+Ns+Ng);
	k_thermal_interface=new_doublevector(Nl+Ns+Ng);
	D=new_doublevector(Nl+Ns+Ng);
	wliq=new_doublevector(Nl+Ns+Ng);
	wice=new_doublevector(Nl+Ns+Ng);
	Temp=new_doublevector(Nl+Ns+Ng);
	deltaw=new_doublevector(Nl+Ns+Ng);
	SWabsf=new_doublevector(Ns+1);
	theta=new_doublevector(Nl);
	ftcl=new_doublevector(Nl);
	froot=new_doublevector(Nl);
	DPsi=new_doublevector(Nl);
	turbulence=new_doublevector(11);
	// INITIALIZATION
	initialize_doublematrix(snow->rho_newsnow,UV->V->co[2]);
	snow->melted_basin=0.0; snow->evap_basin=0.0; snow->subl_basin=0.0; glac->melted_basin=0.0; glac->evap_basin=0.0; glac->subl_basin=0.0;
	initialize_doublevector(snow->CR1,0.0); initialize_doublevector(snow->CR2,0.0); initialize_doublevector(snow->CR3,0.0);

	//SHADOW & CLOUDINESS
	sun((double)(times->hh+times->mm/60.0), times->JD, &(egy->hsun), &(egy->dsun), &E0, par->latitude, par->longitude, par->ST);
	//shadow_n(par->point_sim, top, egy->hsun, egy->dsun, land->shadow);// old routine
	shadow_haiden(par->point_sim, top, egy->hsun, egy->dsun, land->shadow);

	//availability of SW data to calculate tau_cloud

	if(met->column[met->nstsrad-1][iSWb]!=-1 && met->column[met->nstsrad-1][iSWd]!=-1){// there is no station measuring neither SWbeam nor SWdiff
		if(met->var[met->nstsrad-1][met->column[met->nstsrad-1][iSWb]]!=NoV && met->var[met->nstsrad-1][met->column[met->nstsrad-1][iSWd]]!=NoV){// there is both SWbeam AND SWdiff
			SWdata=2;
		}else{
			SWdata=0;// no data of SW
		}
	}else if(met->column[met->nstsrad-1][iSW]!=-1){//there is SWglobal
		if(met->var[met->nstsrad-1][met->column[met->nstsrad-1][iSW]]!=NoV){
			SWdata=1;// just SWglobal
		}else{
			SWdata=0;
		}
	}else{
		SWdata=0;
	}

	if(SWdata>0 && egy->hsun>0){
		find_tau_cloud_station(met->nstsrad, met, par, egy->hsun, E0, Asurr, met->st->sky->co[met->nstsrad], &tau_cloud, &sa);
		//printf("SWdata:%ld tau:%f\n",SWdata,tau_cloud);
	}

	if( met->column[met->nstcloud-1][iC]!=-1 || met->column[met->nstsrad-1][itauC]!=-1){
		if(times->time==0) printf("\nCloudiness data AVAILABLE\n\n");

		if(met->column[met->nstcloud-1][itauC]!=-1){
			tau_cloud_av=met->var[met->nstcloud-1][met->column[met->nstcloud-1][itauC]];
			if(tau_cloud_av<0 || tau_cloud_av==NoV){
				fcloud=-1.0;
			}else if(tau_cloud_av<0.25){
				tau_cloud_av=0.25;
				fcloud=1.0;
			}else if(tau_cloud_av>1){
				tau_cloud_av=1.0;
				fcloud=0.0;
			}else{
				fcloud=pow( (1.0-tau_cloud_av)/0.75, 1.0/3.4 );
			}

		}else{/* no station has data on cloudiness or SW*/
			fcloud=met->var[met->nstcloud-1][met->column[met->nstcloud-1][iC]];
			if(fcloud>1){
				fcloud=1.0;
				tau_cloud_av=0.25;
			}else if(fcloud>0){
				tau_cloud_av=1.0-0.75*pow(fcloud,3.4);
			}else if(fcloud<0 || fcloud==NoV){
				fcloud=-1.0;
			}
		}

	}else{

		if(times->time==0) printf("\nCloudiness data NOT AVAILABLE\n\n");
		fcloud=-1.0;	//novalue initializer

	}

	//printf("SWdata:%ld\n",SWdata);

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
					Precpoint=wat->total->co[par->r_points->co[c]][par->c_points->co[c]];
					Ppoint=met->Pgrid->co[par->r_points->co[c]][par->c_points->co[c]];
				}else{
					Tpoint=met->Tgrid->co[r][c];
					Precpoint=wat->total->co[r][c];
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
				initialize_doublevector(deltaw,0.0);

				//initial condition
				if(times->time==0.0){
					snow_layer_combination(r, c, snow, Tpoint, par->snowlayer_inf, par->Dmin, par->Dmax, times->time);
					if(par->recover!=1){ if(par->glaclayer_max>0) glacier_init_t0(r, c, Tpoint, glac, snow, par, times->time);}

					for(j=1;j<=par->rc->nrh;j++){
						if(r==par->rc->co[j][1] && c==par->rc->co[j][2] && par->superfast!=1) write_soil_output(0, j, times->time, 0.0, par->year0, par->JD0, par->rc, sl, PSImin, par->Esoil);
					}
				}

				//snow
				ns=snow->lnum->co[r][c];
				snowD=DEPTH(r, c, snow->lnum, snow->Dzl);

				//accounting vegetation buried by snow
				fsnow=0.0;
				fc=land->ty->co[lu][jcf];
				LAI=land->ty->co[lu][jLAI];
				//condition on vegetation that can be buried : Hveg > snow_depth_above_which_z0_over_snow_occurs = buried_veg_height
				if(1.E3*land->ty->co[lu][jHveg] > land->ty->co[lu][jz0thresveg] || snowD>land->ty->co[lu][jz0thresveg]){
					fsnow=Fmin(1.0, snowD/land->ty->co[lu][jz0thresveg]);
					fc*=pow(1.0-fsnow,veg_jumping_exp);
					//printf("fc:%f fsnow:%f\n",fc,fsnow);
				}

				//control if LAI is too low (in order to prevent numerical problems)
				if(LAI<LAIthres) fc=0.0;

				//calculates theta from psi(state variable)
				for(l=1;l<=Nl;l++){
					//printf("\nr=%ld, c=%ld, l=%ld, sl->thice->co[l][r][c]=%f, top->Z0->co[r][c]=%f, sy=%d, sl->pa->co[sy][jsat][l]=%f, sl->pa->co[sy][jres][l]=%f, sl->pa->co[sy][ja][l]=%f, sl->pa->co[sy][jns][l]=%f, theta->co[l]=%f",r,c,l,sl->thice->co[l][r][c], top->Z0->co[r][c],sy,sl->pa->co[sy][jsat][l],sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],sl->pa->co[sy][jns][l],theta->co[l]);//stop_execution();
					psisat=psi_saturation(sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l]);
					DPsi->co[l]=Fmax(sl->P->co[l][r][c]-psisat, 0.0);
					sl->P->co[l][r][c]=Fmin(sl->P->co[l][r][c], psisat);
					theta->co[l]=teta_psi(sl->P->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PSImin, par->Esoil);
					if(theta->co[l]!=theta->co[l])printf("theta no value l:%ld teta:%f P:%f DPSi:%f T:%f\n",l,theta->co[l],sl->P->co[l][r][c],DPsi->co[l],sl->T->co[l][r][c]);
				}

				//glacier
				if(par->glaclayer_max>0){
					ng=glac->lnum->co[r][c];
					glacD=DEPTH(r, c, glac->lnum, glac->Dzl);
					if(ng>0){ LAI=0.0; fc=0.0; }
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
				wat->total->co[r][c]=Precpoint;
				// CORRECTION OF PRECIPITATION DUE TO TOPOGRAPHY (JUST FOR 1D SIMULATION)
				if (par->point_sim==1){
					wat->total->co[r][c]*=cos(top->slopes->element[r][c]);
					}
				//distinguish between rain and snow
				part_snow(Precpoint, &Prain_over, &Psnow_over, Tdew, par->T_rain, par->T_snow);
				//modify rain and snow using correction factors
				Prain_over*=par->raincorrfact;
				Psnow_over*=par->snowcorrfact;
				//if(times->time>3000000) Psnow_over=0.0;
				//printf("t:%f Prain:%f Psnow:%f\n",times->time,Prain_over,Psnow_over);

				Prain=(1.-fc)*Prain_over;
				Psnow=(1.-fc)*Psnow_over;

				if(fc>0){
					canopy_rain_interception(0.1, LAI, Prain_over, &max_wcan_rain, &(wat->wcan_rain->co[r][c]), &drip_rain);
					canopy_snow_interception(5.0, LAI, Psnow_over, sl->Tv->co[r][c], Vpoint, par->Dt, &max_wcan_snow, &(wat->wcan_snow->co[r][c]), &drip_snow);
					Prain+=fc*drip_rain;
					Psnow+=fc*drip_snow;
				}else{
					wat->wcan_rain->co[r][c]=0.0;
					wat->wcan_snow->co[r][c]=0.0;
				}

				if(wat->wcan_snow->co[r][c]<0) printf("Error 1 wcansnow:%f %ld %ld\n",wat->wcan_snow->co[r][c],r,c);
				if(wat->wcan_snow->co[r][c]<0) wat->wcan_snow->co[r][c]=0.0;

				snow->Psnow->co[r][c]=Psnow;

				//SHORTWAVE RADIATION
				//initialization of shortwave radiation absorbed by soil
				SW=0.0;

				if(fcloud<0){
					fcloud=find_cloudfactor(Tpoint, RHpoint, top->Z0->co[r][c], met->LRv[2], met->LRv[3]);
					tau_cloud_av=1.0-0.75*pow(fcloud,3.4);
				}
				//in case of shortwave data not available (tau_cloud was initialized to 0) and in case of values too low
				if(tau_cloud<0.1) tau_cloud=tau_cloud_av;

				//printf("taucloud:%f taucloudav:%f ",tau_cloud,tau_cloud_av);

				//calculation of SWin
				tau_atm=atm_transmittance(egy->hsun, Ppoint, RHpoint, Tpoint, Asurr, par->Vis, par->Lozone);
				shortwave_radiation(r, c, egy->hsun, egy->dsun, E0, land->shadow->co[r][c], top->sky->co[r][c], tau_cloud, sa, top->slopes->co[r][c],
					top->aspect->co[r][c], tau_atm, met->var[met->nstsrad-1], met->column[met->nstsrad-1], met->st->sky->co[met->nstsrad],  Asurr, &SWbeam, &SWdiff,
					&cosinc, egy->nDt_shadow, egy->nDt_sun);
				SWin=SWbeam+SWdiff;

				//printf("SWin:%f\n",SWin);

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
				SW+=(1.0-fc)*( SWdiff*(1.0-0.5*avis_d-0.5*anir_d) + SWbeam*(1.0-0.5*avis_b-0.5*anir_b) );

				//shortwave radiation absorbed by canopy
				if(fc>0){
					if(egy->hsun>0){
						if(wat->wcan_snow->co[r][c]<0) printf("Error wcansnow:%f %ld %ld\n",wat->wcan_snow->co[r][c],r,c);
						if(wat->wcan_snow->co[r][c]<0) wat->wcan_snow->co[r][c]=0.0;
						fsnowcan=pow(wat->wcan_snow->co[r][c]/max_wcan_snow, 2./3.);
						shortwave_vegetation(0.5*SWdiff, 0.5*SWbeam, cosinc, fsnowcan, wsn_vis, Bsnd_vis, Bsnb_vis, avis_d, avis_b, land->ty->co[lu][jvCh],
							land->ty->co[lu][jvR_vis], land->ty->co[lu][jvT_vis], LAI, &SWv_vis, &SWg_vis);
						shortwave_vegetation(0.5*SWdiff, 0.5*SWbeam, cosinc, fsnowcan, wsn_nir, Bsnd_nir, Bsnb_nir, anir_d, anir_b, land->ty->co[lu][jvCh],
							land->ty->co[lu][jvR_nir], land->ty->co[lu][jvT_nir], LAI, &SWv_nir, &SWg_nir);
					}else{
						SWv_vis=0.0; SWg_vis=0.0; SWv_nir=0.0; SWg_nir=0.0;
					}
					SW+=fc*(SWg_vis+SWg_nir);
				}else{
					SWv_vis=0.0;
					SWv_nir=0.0;
				}

				//correct in case of reading data
				SW=flux(1, iSWn, met->column, met->var, 1.0, SW);

				//Extinction coefficient for SW in the snow layers
				rad_snow_absorption(r, c, SWabsf, SW, snow);
				set_shallow_snowpack(r, c, par->Dt, snow, SWabsf->co, &Mr_snow, &ns);


				//LONGWAVE RADIATION
				//soil-snow emissivity
				if(snowD>10){
					eps=par->epsilon_snow;
				}else{
					eps=land->ty->co[lu][jemg];
				}
				longwave_radiation(par->state_lwrad, ee, RHpoint, Tpoint, tau_cloud_av, &epsa, &epsa_max, &epsa_min);
				Tsurr=surface(r, c, ns, ng, snow->T, glac->T, sl->T);	//Temperature of surrounding surfaces
				LWin=top->sky->co[r][c]*epsa*SB(Tpoint)+(1-top->sky->co[r][c])*eps*SB(Tsurr);

				//if incoming LW data are available, they are used (priority)
				LWin=flux(met->nstlrad, iLWi, met->column, met->var, 1.0, LWin);

				//roughness lengths
				update_roughness_soil(land->ty->co[lu][jz0], 0.0, 0.0, snowD, land->ty->co[lu][jz0thressoil], par->z0_snow, &z0, &d0, &z0_z0t);
				if(fc>0){
					update_roughness_veg(land->ty->co[lu][jHveg], snowD, zmeas_u, zmeas_T, &z0veg, &d0veg, &hveg);
					root(land->ty->co[lu][jroot], sl->pa->co[sy][jdz], froot->co);
				}


				//SOIL AND SNOW PROPERTIES
				soil_properties(r, c, ns+ng+1, ns+ng+Nl, theta->co, sl->thice->co, D->co, wliq->co, wice->co, k_thermal->co, Temp->co, sl);

				if(ns>0) snow_properties(r, c, 1, ns, D->co, wliq->co, wice->co, k_thermal->co, Temp->co, snow->Dzl->co, snow->w_liq->co,
								snow->w_ice->co, snow->T->co, (*k_thermal_snow_Sturm));

				if(ng>0) snow_properties(r, c, ns+1, ns+ng, D->co, wliq->co, wice->co, k_thermal->co, Temp->co, glac->Dzl->co, glac->w_liq->co,
								glac->w_ice->co, glac->T->co, (*k_thermal_snow_Yen));

				k_interface(ns+ng+Nl, k_thermal->co, D->co, k_thermal_interface->co);

				
				double TsupNp1=0;
				if (met->column[0][iTsup] != -1)
					TsupNp1 = met->var[0][met->column[0][iTsup]];// for Dirichlet Boundary condition at the top of the SOIL column
				if(times->time==0) sl->TsupN=Temp->co[1];
				//ADVECTION
				//if(par->point_sim==0){
				//	if(land->LC2->co[r][c]>0 && snow->Dzl->co[1][r][c]>0){
				//		if(egy->VSFA->co[land->LC2->co[r][c]]!=1.0) Hadv=egy->HSFA->co[land->LC2->co[r][c]]*egy->VSFA->co[land->LC2->co[r][c]]*
				//			adv_efficiency(egy->VSFA->co[land->LC2->co[r][c]])/(1.0-egy->VSFA->co[land->LC2->co[r][c]]);
				//		if(Hadv>egy->HSFA->co[land->LC2->co[r][c]]) Hadv=egy->HSFA->co[land->LC2->co[r][c]];
				//		if(Hadv<0) Hadv=0.0;
				//	}else{
				//		Hadv=0.0;
				//	}
				//}

				//ENERGY BALANCE
				PointEnergyBalance(r, c, ns, ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg, d0veg, 1.0, hveg, Vpoint, Tpoint, Qa, Ppoint, met->LRv[2],
					eps, fc, LAI, &(wat->wcan_rain->co[r][c]), max_wcan_rain, &(wat->wcan_snow->co[r][c]), max_wcan_snow, theta->co, land->ty->co[lu],
					froot->co, ftcl, turbulence, SWin, LWin, SWv_vis+SWv_nir, &LW, &H, &E, &LWv, &Hv, &LEv, &Etrans, &Ts, &Qs, Hadv, SWabsf->co,
					k_thermal_interface->co, D->co, wice->co, wliq->co, deltaw->co, Temp->co, &Hg0, &Hg1, &Eg0, &Eg1, &Qv, &Qg, &rh, &rv, &rb, &rc,
					&rh_ic, &rv_ic, &u_top, sl, par, &decaycoeff, &Locc, sl->TsupN,TsupNp1);

				sl->TsupN=TsupNp1;// memorize the top Dirichlet boundary condition
				if(wat->wcan_snow->co[r][c]<0) printf("Error 2 wcansnow:%f %ld %ld\n",wat->wcan_snow->co[r][c],r,c);
				if(wat->wcan_snow->co[r][c]<0) wat->wcan_snow->co[r][c]=0.0;

				LE=E*latent(Temp->co[1],Levap(Temp->co[1]));
				surfEB = SW + LW - H - LE;

				if(ns==0){
					G=surfEB;
				}else{
					G=k_thermal_interface->co[ns]*(Temp->co[ns]-Temp->co[ns+1])/(0.5*D->co[ns]+0.5*D->co[ns+1]);
				}

				liqWBsnow(r, c, snow, &Mr_snow, &Prain_on_snow, par, top->slopes->co[r][c], Prain, wice->co, wliq->co, Temp->co, E*par->Dt);
				iceWBsnow(r, c, snow, Psnow, Tpoint);
				snowD=DEPTH(r, c, snow->lnum, snow->Dzl);

				snow_layer_combination(r, c, snow, Tpoint, par->snowlayer_inf, par->Dmin, par->Dmax, times->time);

				if(par->glaclayer_max>0){
					WBglacier(ns, r, c, glac, &Mr_glac, par, wice->co, wliq->co, Temp->co, E*par->Dt);
					glac2snow(r, c, snow, glac, par->Dmin, par->Dmax);
					if(par->glaclayer_max>1)glac_layer_combination(r, c, glac, Tpoint, Ng, par->Dmin_glac, par->Dmax_glac, times->time);
					ng=glac->lnum->co[r][c];
				}

				//NET PRECIPITATION
				wat->Pn->co[r][c] = Mr_snow;

				if( LAI>=LAIthres && ng==0 && ( 1.E3*land->ty->co[lu][jHveg] > land->ty->co[lu][jz0thresveg] || (fsnow==0 && snowD>land->ty->co[lu][jz0thresveg]) ) ){

					fc0=fc;
					fsnownew=Fmin(1.0, snowD/land->ty->co[lu][jz0thresveg]);
					fc=pow(1.0-fsnownew,veg_jumping_exp)*land->ty->co[lu][jcf];

					//a) fc increases
					if(fc>fc0){
						wat->wcan_rain->co[r][c]*=(fc0/fc);
						wat->wcan_snow->co[r][c]*=(fc0/fc);

					//b) fc decreases
					}else{
						wat->Pn->co[r][c] +=  wat->wcan_rain->co[r][c]*(fc0-fc);
						iceWBsnow(r, c, snow, wat->wcan_snow->co[r][c]*(fc0-fc), sl->Tv->co[r][c]);

					}
				}

				wat->Pn->co[r][c] += (Mr_glac + Prain/par->Dt);

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

				if(par->point_sim!=1) prepare_output(Er_soil, Mr_snow, Er_snow, Sr_snow, Mr_glac, Er_glac, Sr_glac, Prain, Psnow_over, egy, wat, snow, glac, land, top, sl, met, times,
					par, r, c, 0.25*(avis_d+anir_d+avis_b+anir_b), LE, surfEB, H, G, Temp->co[1], SWin, SWin-SW, SWbeam, eps, LWin, LW-LWin, cosinc, Precpoint);

				output_pixel(r, c, Psnow, Prain-Prain_on_snow, Prain_on_snow, Sr_soil, Er_soil, Mr_snow, Er_snow, Sr_snow, Mr_glac, Er_glac, Sr_glac,
					Etrans, LE, H, surfEB, G, Temp->co[1], 0.25*(avis_d+anir_d+avis_b+anir_b), eps, LAI, LWin, SWin, LWin-LW, SWin-SW, epsa, epsa_min,
					epsa_max, SWbeam, SWdiff, DTcorr, Tdew, times->n_pixel, turbulence, par, wat, egy, top, met, snow, glac, land, Tpoint, Ppoint, Vpoint,
					RHpoint, Psnow_over, Prain_over, z0, z0veg, d0veg, (SWv_vis+SWv_nir), LWv, fc*Hv, fc*LEv, sl->Tv->co[r][c], Ts, Hg0, Hg1, Eg0, Eg1,
					fc, rh, rv, rb, rc, rh_ic, rv_ic, Hv, LEv, Qv, Qg, Qa, Qs, u_top, decaycoeff, Locc);

				output_basin(times->n_basin, Prain, Psnow, Tpoint, Temp->co[1], sl->Tv->co[r][c], Er_soil+Sr_soil, Etrans, 0.0, LE, H, SW, LW, fc*LEv,
					fc*Hv, fc*(SWv_vis+SWv_nir), fc*LWv, SWin, wat->out2->co, egy->out2->co, par->Dt);

				if(par->ES_num>0) output_altrank(par->ES_num, Prain, Psnow, Temp->co[1], Er_soil, Sr_soil, Etrans, 0.0, LE, H, surfEB, SWin, SWin-SW, LWin, eps,
					top->sky->co[r][c], wat->Pn->co[r][c], wat->Pn->co[r][c]-sl->Jinf->co[r][c], Er_snow, Sr_snow, Mr_snow, Er_glac, Sr_glac, Mr_glac, par->Dt, times->n_basin,\
					top->Z0->co[r][c], top->Zmin, top->Zmax, glacD, par->glac_thr, egy->out3->co);

				output_map_plots(r, c, par, times->time, times->n_plot, egy, met, snow, H, LE, fc*Hv, fc*LEv, SWin, SW, fc*(SWv_vis+SWv_nir), LWin, LW, fc*LWv, Ts,
					Temp->co[1], sl->Tv->co[r][c]);

				egy->Hgrid->co[r][c]=H+fc*Hv;
				egy->Tsgrid->co[r][c]=Temp->co[1];

				for(l=1;l<=Nl;l++){
					sl->P->co[l][r][c]+=DPsi->co[l];
					if(sl->P->co[l][r][c]!=sl->P->co[l][r][c])printf("Psi no value l:%ld teta:%f P:%f T:%f\n",l,theta->co[l],sl->P->co[l][r][c],sl->T->co[l][r][c]);
				}
			}
		}
	}


	//ACCOUNT FOR SNOW WIND TRANSPORT

	//PBSM
	if(par->blowing_snow==1) set_windtrans_snow(snow, met, land, par, times->time);

	//SnowTran3D
	//if(par->blowing_snow==1) set_windtrans_snow2(snow, met, land, top, par, times->time);

	//CALCULATE ADVECTION PARAMETERS FROM SNOW FREE TO SNOW COVERED AREA
	/*if(par->point_sim!=1){
		find_SCA(land->LC2->co, snow, par, land->LC->co, times->time);
		//snow_fluxes_H(land->LC2->co, snow, egy->Tsgrid, egy->Hgrid, land->LC, met->Tgrid, times, par, egy->VSFA->co, egy->HSFA->co);
	}*/

	//PREPARE SNOW OUTPUT
	output_snow(snow, land->LC->co, par);


	//DEALLOCATION
	free_doublevector(k_thermal);
	free_doublevector(k_thermal_interface);
	free_doublevector(D);
	free_doublevector(wliq);
	free_doublevector(wice);
	free_doublevector(Temp);
	free_doublevector(deltaw);
	free_doublevector(theta);
	free_doublevector(ftcl);
	free_doublevector(turbulence);
	free_doublevector(SWabsf);
	free_doublevector(froot);
	free_doublevector(DPsi);

}
//end of "energy_balance" subroutine

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
// C in [J m^-3 K^-1]



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double calc_C(long l, long r, long c, long nsng, double *wi, double *wl, double *dw, double *D, SOIL *sl){

	double C;
	short sy = sl->type->co[r][c];

	if(l<=nsng){	//snow
		C = (c_ice*(wi[l]-dw[l]) + c_liq*(wl[l]+dw[l]))/D[l];
	}else{	//soil
		C = sl->pa->co[sy][jct][l-nsng]*(1.-sl->pa->co[sy][jsat][l-nsng]) + (c_ice*(wi[l]-dw[l]) + c_liq*(wl[l]+dw[l]))/D[l];
	}

	return(C);

}

double calc_C0(long l, long r, long c, long nsng, double *wi, double *wl, double *dw, double *D, SOIL *sl){

	double C;
	short sy = sl->type->co[r][c];

	if(l<=nsng){	//snow
		C = (c_ice*(wi[l]-0.0*dw[l]) + c_liq*(wl[l]+0.0*dw[l]))/D[l];
	}else{	//soil
		C = sl->pa->co[sy][jct][l-nsng]*(1.-sl->pa->co[sy][jsat][l-nsng]) + (c_ice*(wi[l]-0.0*dw[l]) + c_liq*(wl[l]+0.0*dw[l]))/D[l];
	}

	return(C);

}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double calc_k(long l, long r, long c, long ns, long ng, double *wi, double *wl, double *dw, double *T, double *D, double (* kfunct_snow)(double rho),
	double (* kfunct_glac)(double rho), SOIL *sl, PAR *par){

	double k;
	short sy = sl->type->co[r][c];

	if(l<=ns){
		k = (*kfunct_snow)((wl[l]+wi[l])/D[l]);
	}else if(l<=ng+ns){
		k = (*kfunct_glac)((wl[l]+wi[l])/D[l]);
	}else{
		k = k_thermal_soil((wl[l]+dw[l])/(rho_w*D[l]), (wi[l]-dw[l])/(rho_w*D[l]), sl->pa->co[sy][jsat][l-ns-ng], sl->T->co[l-ns-ng][r][c], sl->pa->co[sy][jkt][l-ns-ng]);
	}

	return(k);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

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


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void soil_properties(long r, long c, long beg, long end, double *th_l, double ***th_i, double *D, double *wl, double *wi, double *k, double *T, SOIL *sl){

	long l, m;
	short sy=sl->type->co[r][c];

	for(l=beg;l<=end;l++){

		//sl index
		m=l-beg+1;

		//liquid and ice content
		D[l]=1.0E-3*sl->pa->co[sy][jdz][m];	//[m]
		wl[l]=rho_w*D[l]*th_l[m];			//[kg m^(-2)]
		wi[l]=rho_w*D[l]*th_i[m][r][c];		//[kg m^(-2)]

		//thermal conductivity [W m^-1 K^-1] - from Farouki
		k[l]=k_thermal_soil(th_l[m], th_i[m][r][c], sl->pa->co[sy][jsat][m], sl->T->co[m][r][c], sl->pa->co[sy][jkt][m]);

		//temperatura [Deg Celsius]
		T[l]=sl->T->co[m][r][c];
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double k_thermal_soil(double th_liq, double th_ice, double th_sat, double T, double k_solid){

	double rho_dry, sr, ke, k_dry, k_sat, ir, k;

	//dry sl density [kg/mc]
	rho_dry=2700.0*(1.0-th_sat);

	//saturation degree
	sr=(th_liq+th_ice)/th_sat;

	//Kersten number
	if(T>=Tfreezing){
		ke=log(sr)+1.0;
		if(ke<=0.0) ke=0.0;
	}else{
		ke=sr;
	}

	//dry soil thermal conductivity [W m^-1 K^-1]
	k_dry=(0.135*rho_dry+64.7)/(2700.0-0.947*rho_dry);

	//soil thermal conductivity [W m^-1 K^-1] Farouki (1981)
	if(sr>1.0E-7){
		//saturated sl thermal conductivity [W m^-1 K^-1]
		if(T>Tfreezing){
			k_sat=pow(k_solid,1.0-th_sat)*pow(k_liq,th_sat);
		}else{
			ir=th_ice/(th_liq+th_ice);
			k_sat=pow(k_solid,1.0-th_sat)*pow(k_liq,th_sat*(1-ir))*pow(k_ice,th_sat*ir);
		}
		//soil thermal conductivity [W m^-1 K^-1]
		k=ke*k_sat + (1.0-ke)*k_dry;
	}else{
		k=k_dry;
	}

	return(k);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void k_interface(long n, double *k, double *D, double *ki){

	long l;

	for(l=1;l<=n;l++){
		if(l<n){
			ki[l]=(k[l]*k[l+1]*0.5*(D[l]+D[l+1]))/(k[l]*0.5*D[l+1]+k[l+1]*0.5*D[l]);
		}else{
			ki[l]=k[l];
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void prepare_output(double Er_soil, double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, double prec_rain, double prec_snow_atm,
			ENERGY *egy, WATER *wat, SNOW *snow, GLACIER *glac, LAND *land, TOPO *top, SOIL *sl, METEO *met, TIMES *times, PAR *par, long r, long c, double A,
			double LE, double surfEB, double H, double G, double Ts, double SWin, double SWout, double SWbeam, double eps, double LWin, double LWout, double cosinc, double Precpoint){

	FILE *f;
	long l;

	//STRUCT SNOW
	snow->melted_basin+=Mr_snow*par->Dt;	//[mm]
	snow->evap_basin+=Er_snow*par->Dt;
	snow->subl_basin+=Sr_snow*par->Dt;

	//STRUCT GLAC
	glac->melted_basin+=Mr_glac*par->Dt;
	glac->evap_basin+=Er_glac*par->Dt;
	glac->subl_basin+=Sr_glac*par->Dt;

	//OUTPUT DISTRIBUITI OPZIONALI
	if(par->output_balancesn>0){
		snow->MELTED->co[r][c]+=Mr_snow*par->Dt;
		snow->SUBL->co[r][c]+=(Sr_snow+Er_snow)*par->Dt;
		if(snow->type->co[r][c]==2){
			snow->t_snow->co[r][c]+=par->Dt/86400.0;
			for(l=1;l<=snow->lnum->co[r][c];l++){
				snow->totav_snow->co[r][c]+=snow->Dzl->co[l][r][c];
			}
		}
	}

	if(par->output_balancegl>0 && par->glaclayer_max>0) glac->MELTED->co[r][c]+=Mr_glac*par->Dt;
	if(par->output_balancegl>0 && par->glaclayer_max>0) glac->SUBL->co[r][c]+=(Sr_glac+Er_glac)*par->Dt;


	if(par->output_P>0){
		wat->PrTOT_mean->co[r][c]+=wat->total->co[r][c];	/*[mm]*/
		wat->PrSNW_mean->co[r][c]+=prec_snow_atm;								/*[mm]*/
	}

	if(par->output_albedo>0){
		if(SWout>0){
			land->albedo->co[r][c]+=(SWout/SWin)/((par->output_albedo*3600.0)/(par->Dt));
		}else{
			land->albedo->co[r][c]+=A/((par->output_albedo*3600.0)/(par->Dt));
		}
	}

	if(par->output_ET>0){
		egy->ET_mean->co[r][c]+=LE/((par->output_ET*3600.0)/(par->Dt)); //[W/m^2]
		/*if(snow->totav_snow->co[r][c]>0 && LE>=0){
			egy->ET_mean->co[r][c]+=LE*par->Dt*1.E-6;
		}else{
			egy->ET_mean2->co[r][c]+=LE*par->Dt*1.E-6;
		}*/
		if(par->distr_stat==1){
			if(egy->ET_max->co[r][c]<LE) egy->ET_max->co[r][c]=LE;
			if(egy->ET_min->co[r][c]>LE) egy->ET_min->co[r][c]=LE;
		}
	}

	if(par->output_G>0){
		egy->G_mean->co[r][c]+=surfEB/((par->output_G*3600.0)/(par->Dt)); //[W/m^2]
		//if(snow->totav_snow->co[r][c]>0)egy->G_mean->co[r][c]+=surfEB*par->Dt*1.E-6;
		if(par->distr_stat==1){
			if(egy->G_max->co[r][c]<surfEB) egy->G_max->co[r][c]=surfEB;
			if(egy->G_min->co[r][c]>surfEB) egy->G_min->co[r][c]=surfEB;
		}
		egy->G_snowsoil->co[r][c]+=G/((par->output_G*3600.0)/(par->Dt)); //[W/m^2]
	}


	if(par->output_H>0){
		egy->H_mean->co[r][c]+=H/((par->output_H*3600.0)/(par->Dt)); //[W/m^2]
		/*if(snow->totav_snow->co[r][c]>0 && H>=0){
			egy->H_mean->co[r][c]+=H*par->Dt*1.E-6;
		}else{
			egy->H_mean2->co[r][c]+=H*par->Dt*1.E-6;
		}*/
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
		egy->Rn_mean->co[r][c]+=(SWin-SWout + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)))/((par->output_Rn*3600.0)/(par->Dt)); //[W/m^2]
		if(fabs(SWin-SWout + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)))>1500){
			f=fopen(error_file_name,"a");
			fprintf(f,"\ntime=%10.1f r=%4ld c=%4ld Rn=%10.3f Rsw=%10.3f Rlwdiff=%10.3f albedo=%10.8f eps=%10.8fTa=%10.5f Ts=%10.5f Rsw_meas=%f sin(alpha)=%f cos(inc)=%f\n",
				times->time,r,c,SWin-SWout + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)),SWin,LWin,A,
				eps,met->Tgrid->co[r][c],Ts,met->var[0][iSW],sin(egy->hsun),cosinc);
			fprintf(f,"\nH:%f LE:%f\n",H,LE);
			fclose(f);
		}

		//egy->LW_in->co[r][c]+=( top->sky->co[r][c]*5.67E-8*eps*pow((Ts+tk),4.0) )/((par->output_Rn*3600.0)/(par->Dt));
		//egy->LW_out->co[r][c]+=( top->sky->co[r][c]*eps*LWin )/((par->output_Rn*3600.0)/(par->Dt));
		//egy->SW->co[r][c]+=( SWin-SWout )/((par->output_Rn*3600.0)/(par->Dt));

		//if(snow->totav_snow->co[r][c]>0){
			egy->LW_in->co[r][c]+=LWin*par->Dt*1.E-6;
			egy->LW_out->co[r][c]+=LWout*par->Dt*1.E-6;
			egy->SW->co[r][c]+=(SWin-SWout)*par->Dt*1.E-6;
		//}

		if(par->distr_stat==1){
			if(egy->LW_max->co[r][c]<LWin-LWout) egy->LW_max->co[r][c]=LWin-LWout;
			if(egy->LW_min->co[r][c]>LWin-LWout) egy->LW_min->co[r][c]=LWin-LWout;
			if(egy->Rn_max->co[r][c]<SWin-SWout + LWin-LWout) egy->Rn_max->co[r][c]=(SWin-SWout + LWin-LWout);
			if(egy->Rn_min->co[r][c]>SWin-SWout + LWin-LWout) egy->Rn_min->co[r][c]=(SWin-SWout + LWin-LWout);
			if(egy->SW_max->co[r][c]<SWin-SWout) egy->SW_max->co[r][c]=SWin-SWout;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_pixel(long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double Sr_soil, double Er_soil,
	double Mr_snow, double Er_snow, double Sr_snow, double Mr_glac, double Er_glac, double Sr_glac, double Evt, double LE, double H, double surfEB,
	double G, double Tg, double A, double eps, double LAI, double LWin, double SWin, double LWout, double SWout, double epsa, double epsa_min,
	double epsa_max, double SWbeam, double SWdiff, double DTcorr, double Tdew, long n, DOUBLEVECTOR *turbulence, PAR *par, WATER *wat, ENERGY *egy,
	TOPO *top, METEO *met, SNOW *snow, GLACIER *glac, LAND *land, double Tpoint, double Ppoint, double Vpoint, double RHpoint, double prec_snow_atm,
	double prec_rain_atm, double z0soil, double z0, double d0, double SWv, double LWv, double Hv, double LEv, double Tv, double Ts, double Hg0,
	double Hg1, double Eg0, double Eg1, double fc, double rh, double rv, double rb, double rc, double rh_ic, double rv_ic, double Hv1, double LEv1,
	double Qv, double Qg, double Qa, double Qs, double u_top, double decay, double Locc){

	long i, l;
	long j=0; /* modified by Emanuele Cordano on 24/9/9 (j inizialization) */
	double ea, es, de_dT, Q;

	if(par->state_pixel==1){
		for(i=1;i<=par->chkpt->nrh;i++){
			if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){

				wat->out1->co[3][i]+=prec_snow;
				wat->out1->co[4][i]+=(prec_rain_on_soil+prec_rain_on_snow);

				wat->out1->co[16][i]+=prec_snow_atm;
				wat->out1->co[17][i]+=prec_rain_atm;
				//wat->out1->co[18][i]+=maxstorage_c/(double)n;
				//wat->out1->co[19][i]+=evap_c;
				//wat->out1->co[20][i]+=drip_c;

				wat->out1->co[21][i]+=prec_snow;
				wat->out1->co[22][i]+=prec_rain_on_soil;
				wat->out1->co[23][i]+=prec_snow_atm;
				wat->out1->co[24][i]+=prec_rain_atm;
				//wat->out1->co[25][i]+=evap_c;
				//wat->out1->co[26][i]+=drip_c;
				wat->out1->co[27][i]+=prec_rain_on_snow;

				egy->out1->co[1][i]+=Er_soil*par->Dt;	//Eg[mm]
				egy->out1->co[2][i]+=Sr_soil*par->Dt;	//Sg[mm]
				//egy->out1->co[3][i]+=Epc*fec*(1000.0/rho_w)*par->Dt;	//Evc[mm]
				egy->out1->co[4][i]+=Evt*(1000.0/rho_w)*par->Dt;	//Etc[mm]
				egy->out1->co[5][i]+=LE/(double)n; //ET[W/m^2]
				egy->out1->co[6][i]+=H/(double)n; //H[W/m^2]
				egy->out1->co[7][i]+=surfEB/(double)n;
				egy->out1->co[8][i]+=G/(double)n;
				egy->out1->co[9][i]+=((prec_rain_on_soil+prec_rain_on_snow)*Lf/par->Dt)/(double)n; //Qrain[W/m^2]
				egy->out1->co[10][i]+=Tg/(double)n; //Ts[C]
				//Rnet and Rnet_cumulated
				egy->out1->co[11][i]+=(SWin-SWout+eps*LWin-LWout)/(double)n;
				egy->out1->co[21][i]+=(SWin-SWout+eps*LWin-LWout)*par->Dt*1.0E-6;
				//Rlw_in and Rlw_in_cumulated
				egy->out1->co[12][i]+=LWin/(double)n;
				egy->out1->co[22][i]+=LWin*par->Dt*1.0E-6;
				//Rlw_out and Rlw_out_cumulated
				egy->out1->co[13][i]-=LWout/(double)n;
				egy->out1->co[23][i]-=LWout*par->Dt*1.0E-6;
				//Rsw, Rsw_incoming_cumulated, albedo, Rsw_outcoming_cumulated
				egy->out1->co[14][i]+=SWin/(double)n;
				egy->out1->co[24][i]+=SWin*par->Dt*1.0E-6;
				/*if(SWout>0){
					egy->out1->co[15][i]+=(SWout/SWin)/(double)n;
				}else{
					egy->out1->co[15][i]+=A/(double)n;
				}*/
				egy->out1->co[28][i]-=SWout*par->Dt*1.0E-6;
				egy->out1->co[45][i]-=SWout/(double)n;
				//atmosphere emissivity
				egy->out1->co[16][i]+=epsa/(double)n;
				//wind speed
				egy->out1->co[17][i]+=Vpoint/(double)n;
				//relative humidity
				egy->out1->co[18][i]+=RHpoint/(double)n;
				//atmospheric pressure
				egy->out1->co[19][i]+=Ppoint/(double)n;
				//air temperature
				egy->out1->co[20][i]+=Tpoint/(double)n;
				//ET cumulated
				egy->out1->co[25][i]+=LE*par->Dt*1.0E-6;
				//H cumulated
				egy->out1->co[26][i]+=H*par->Dt*1.0E-6;
				//G cumulated
				egy->out1->co[27][i]+=surfEB*par->Dt*1.0E-6;

				//L Obukhov
				egy->out1->co[29][i]+=turbulence->co[2]/(double)n;
				//number of iteration
				egy->out1->co[30][i]+=turbulence->co[1]/(double)n;

				//CH, CL (transfer coefficient for H and L)
				/*egy->out1->co[35][i]+=turbulence->co[10]/(double)n;
				egy->out1->co[36][i]+=turbulence->co[11]/(double)n;
				egy->out1->co[37][i]+=DTcorr/(double)n;
				egy->out1->co[38][i]+=turbulence->co[3]/(double)n;
				egy->out1->co[39][i]+=turbulence->co[4]/(double)n;
				egy->out1->co[40][i]+=turbulence->co[5]/(double)n;
				egy->out1->co[41][i]+=turbulence->co[6]/(double)n;
				egy->out1->co[42][i]+=turbulence->co[7]/(double)n;
				egy->out1->co[43][i]+=turbulence->co[8]/(double)n;
				egy->out1->co[44][i]+=turbulence->co[9]/(double)n;*/

				//ea(mbar) Q(-)
				egy->out1->co[46][i]+=Tdew/(double)n;
				sat_vap_pressure(&ea, &de_dT, Tpoint, Ppoint);
				Q=RHpoint*0.622*ea/(Ppoint-0.378*ea);
				ea=Q*Ppoint/(0.622+Q*0.378);
				egy->out1->co[47][i]+=ea/(double)n;
				egy->out1->co[48][i]+=Q/(double)n;
				sat_vap_pressure(&es, &de_dT, Ts, Ppoint);
				egy->out1->co[49][i]+=es/(double)n;
				egy->out1->co[50][i]+=0.622*es/(Ppoint-0.378*es)/(double)n;

				if(epsa_min>0){
					egy->out1->co[51][i]+=(epsa_min*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
				}else{
					egy->out1->co[51][i]=UV->V->co[2];
				}
				if(epsa_max>0){
					egy->out1->co[52][i]+=(epsa_max*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
				}else{
					egy->out1->co[52][i]=UV->V->co[2];
				}
				egy->out1->co[53][i]+=(SWbeam/(double)n);
				egy->out1->co[54][i]+=(SWdiff/(double)n);

				if(par->micromet==1){
					egy->out1->co[55][i]+=(met->Vdir->co[r][c])/(double)n;
				}else{
					egy->out1->co[55][i]=UV->V->co[2];
				}
				egy->out1->co[56][i]+=LAI/(double)n;

				egy->out1->co[57][i]+=z0/(double)n;
				egy->out1->co[58][i]+=d0/(double)n;

				snow->melted->co[i]+=Mr_snow*par->Dt;	//[mm]
				snow->evap->co[i]+=Er_snow*par->Dt;	//[mm]
				snow->subl->co[i]+=Sr_snow*par->Dt;	//[mm]

				glac->melted->co[i]+=Mr_glac*par->Dt;	//[mm]
				glac->evap->co[i]+=Er_glac*par->Dt;	//[mm]
				glac->subl->co[i]+=Sr_glac*par->Dt;	//[mm]

				for(l=1;l<=par->snowlayer_max;l++){
					snow->CR1m->co[l]+=(snow->CR1->co[l]*3600.0/(double)n);
					snow->CR2m->co[l]+=(snow->CR2->co[l]*3600.0/(double)n);
					snow->CR3m->co[l]+=(snow->CR3->co[l]*3600.0/(double)n);
				}

				egy->out1->co[59][i]+=SWv/(double)n;
				egy->out1->co[60][i]+=LWv/(double)n;
				egy->out1->co[61][i]+=Hv/(double)n;
				egy->out1->co[62][i]+=LEv/(double)n;
				egy->out1->co[63][i]+=Tv/(double)n;
				egy->out1->co[64][i]+=Ts/(double)n;

				egy->out1->co[65][i]+=Hg0/(double)n;
				egy->out1->co[66][i]+=Levap(Tg)*Eg0/(double)n;
				egy->out1->co[67][i]+=Hg1/(double)n;
				egy->out1->co[68][i]+=Levap(Tg)*Eg1/(double)n;
				egy->out1->co[69][i]+=fc/(double)n;

				egy->out1->co[70][i]+=(1./rh)/(double)n;
				egy->out1->co[71][i]+=(1./rv)/(double)n;
				egy->out1->co[72][i]+=(1./rb)/(double)n;
				egy->out1->co[73][i]+=(1./rc)/(double)n;
				egy->out1->co[74][i]+=(1./rh_ic)/(double)n;
				egy->out1->co[77][i]+=(1./rv_ic)/(double)n;

				egy->out1->co[75][i]+=(Hv1)/(double)n;
				egy->out1->co[76][i]+=(LEv1)/(double)n;

				egy->out1->co[78][i]+=(Qv)/(double)n;
				egy->out1->co[79][i]+=(Qg)/(double)n;
				egy->out1->co[80][i]+=(Qa)/(double)n;
				egy->out1->co[81][i]+=(Qs)/(double)n;

				egy->out1->co[82][i]+=(u_top)/(double)n;

				egy->out1->co[83][i]+=(decay)/(double)n;
				egy->out1->co[84][i]+=(Locc)/(double)n;


			}
		}
	}

	//OUTPUT precipitation+melting for each land use (split in snow covered and snow free)
	for(i=1;i<=land->clax->nh;i++){
		if((short)land->LC->co[r][c]==land->clax->co[i]) j=i;
	}
	if(snow->type->co[r][c]==2){
		wat->outfluxes->co[1][j]+=prec_rain_on_soil;
		wat->outfluxes->co[2][j]+=Mr_snow*par->Dt;
		wat->outfluxes->co[3][j]+=Mr_glac*par->Dt;
		land->cont->co[j][2]+=1;
	}else{
		wat->outfluxes->co[4][j]+=prec_rain_on_soil;
		wat->outfluxes->co[5][j]+=Mr_snow*par->Dt;
		wat->outfluxes->co[6][j]+=Mr_glac*par->Dt;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void output_basin(long n, double prec_rain, double prec_snow, double Ta, double Tg, double Tv, double Eg, double Evt, double Eve, double LE, double H, double SW, double LW,
	double LEv, double Hv, double SWv, double LWv, double SWin, double *wat, double *en, double Dt){

	//OUTPUT MEDI PER BACINO
	wat[1]+=prec_rain;																/*[mm]*/
	wat[2]+=prec_snow;
	en[1]+=Ta/(double)n;
	en[2]+=Tg/(double)n;															/*Ts[C]*/
	en[3]+=Eg*Dt;																	/*Eg[mm]*/
	en[4]+=Eve*Dt;																	/*Evc[mm]*/
	en[5]+=Evt*Dt;																	/*Etc[mm]*/
	en[6]+=LE/(double)n;															/*ET[W/m2]*/
	en[7]+=H/(double)n;																/*H[W/m2]*/
	en[8]+=SW/(double)n;
	en[9]+=LW/(double)n;
	en[10]+=LEv/(double)n;
	en[11]+=Hv/(double)n;
	en[12]+=SWv/(double)n;
	en[13]+=LWv/(double)n;
	en[14]+=SWin/(double)n;
	en[15]+=Tv/(double)n;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void output_altrank(long ES, double prec_rain, double prec_snow, double Ts, double Er_soil, double Sr_soil, double Evt, double Eve, double LE, double H, double surfEB, double SWin, double SWout,
	double LWin, double eps, double V, double Pn, double runoff, double Er_snow, double Sr_snow, double Mr_snow, double Er_glac, double Sr_glac, double Mr_glac, double Dt, long n, double Z, double Zmin,
	double Zmax, double glacD, double glacDmin, double **out1){

	long i;

	//OUTPUT MEDI PER FASCE ALTIMETRICHE
	for(i=1;i<=ES;i++){
		if( ( (Z>=Zmin+(i-1)*(Zmax-Zmin)/(double)ES && Z<Zmin+i*(Zmax-Zmin)/(double)ES) || (Z==Zmin+i*(Zmax-Zmin)/(double)ES && i==ES) ) && glacD>=glacDmin ){
			out1[1][i]+=Er_soil*Dt;					//Eg[mm]
			out1[2][i]+=Sr_soil*Dt;					//Sg[mm]
			out1[3][i]+=Eve*(1000.0/rho_w)*Dt;		//Evc[mm]
			out1[4][i]+=Evt*(1000.0/rho_w)*Dt;		//Etc[mm]
			out1[5][i]+=LE/(double)n;				//ET[W/m^2]
			out1[6][i]+=H/(double)n;					//H[W/m^2]
			out1[7][i]+=surfEB/(double)n;			//G[W/m^2]
			//out1[8][i]+=Eimm/(double)n;				//Ecanopy[W/m^2]
			out1[9][i]+=(prec_rain*Lf/Dt)/(double)n; //Qrain[W/m^2]
			out1[10][i]+=Ts/(double)n;				//Ts[C]
			out1[11][i]+=(SWin-SWout + V*(eps*LWin - eps*SB(Ts)))/(double)n; //[W/m^2]
			out1[12][i]+=prec_rain;					//[mm]
			out1[13][i]+=prec_snow;					//[mm]
			out1[14][i]+=Pn*Dt;						//[mm]
			out1[15][i]+=runoff*Dt;					//[mm]
			out1[16][i]+=Mr_snow*Dt;
			out1[17][i]+=Er_snow*Dt;
			out1[18][i]+=Sr_snow*Dt;
			out1[20][i]+=Mr_glac*Dt;
			out1[21][i]+=Er_glac*Dt;
			out1[22][i]+=Sr_glac*Dt;
		}
	}
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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void update_soil(long l, long r, long c, double evap, double *theta, double *theta_ice, double *T, SOIL *sl, PAR *par){

	short sy=sl->type->co[r][c];
	//double h;

	if(*theta <= 1.1*sl->pa->co[sy][jres][l]){
		evap = Fmin(0.0, evap);
	}else{
		if(*theta - evap/sl->pa->co[sy][jdz][l] <= 1.1*sl->pa->co[sy][jres][l]){
			evap = ( 1.1*sl->pa->co[sy][jres][l] - (*theta) )*sl->pa->co[sy][jdz][l];
		}
	}

	//h=internal_energy_soil(*theta, *theta_ice, *T, sl->pa->co[sy][jdz][l], sl->pa->co[sy][jct][l], sl->pa->co[sy][jsat][l]);

	*theta = (*theta) - evap/sl->pa->co[sy][jdz][l];

	//from_internal_soil_energy(r, c, l, h-Lf*evap, theta, theta_ice, T, sl->pa->co[sy], PSImin);

}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void output_map_plots(long r, long c, PAR *par, double t, long n, ENERGY *egy, METEO *met, SNOW *snow, double Hg, double LEg, double Hv, double LEv, double SWin,
	double SWg, double SWv, double LWin, double LWg, double LWv, double Ts, double Tg, double Tv){

	long j, d, M, y, h, m;
	double tmin, tmax, JD, N;

	if(par->JD_plots->co[1]>=0){

		date_time(t, par->year0, par->JD0, 0.0, &JD, &d, &M, &y, &h, &m);

		for(j=1;j<=(long)(par->JD_plots->nh/2.);j++){

			get_time( &tmin, par->JD_plots->co[2*j-1], y, par->JD0, par->year0 );
			get_time( &tmax, par->JD_plots->co[2*j]  , y, par->JD0, par->year0 );

			if( floor((tmax-tmin)/(n*par->Dt)) != (tmax-tmin)/(n*par->Dt) ){
				N=floor(((tmax-tmin)/(n*par->Dt)));
				tmax=tmin+N*n*par->Dt;
			}

			if(t>=tmin && t<tmax){
				egy->Hgplot->co[r][c]+=Hg/(double)n;
				if(r==2 && c==2)printf("t:%f tmin:%f tmax:%f Hg:%f Hgplot:%f n:%ld\n",t,tmin,tmax,Hg,egy->Hgplot->co[r][c],n);

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

void check_errors(long r, long c, long n, DOUBLEVECTOR *adi, DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, DOUBLEVECTOR *e, double *T, SHORTVECTOR *mf){
	/* warnings corrected by Emanuele Cordano on 24 Sept 2009 */
	long l;
	short occ=0;

	for(l=1;l<=n;l++){
		if(e->co[l]!=e->co[l]) occ=1;
	}

	if(occ==1 ){
		printf("NOvalue in Energy Balance: r:%ld c:%ld\n",r,c);
		printf("l:%d adi:--- ad:%f ads:%f b:%f e:%f T:%f mf:%d\n",1,ad->co[1],ads->co[1],b->co[1],e->co[1],T[1],mf->co[1]);
		for(l=2;l<=n-1;l++){
			printf("l:%ld adi:%f ad:%f ads:%f b:%f e:%f T:%f mf:%d\n",l,adi->co[l-1],ad->co[l],ads->co[l],b->co[l],e->co[l],T[l],mf->co[l]);
		}
		printf("l:%ld adi:%f ad:%f ads:--- b:%f e:%f T:%f mf:%d\n",n,adi->co[n-1],ad->co[n],b->co[n],e->co[n],T[n],mf->co[n]);
		stop_execution();
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/





/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/



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

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD, double thres, double z0snow, double *z0_ris, double *d0_ris, double *z0_z0t_ris){
	/* Author:  Stefano Endrizzi  Year:
	 * function that updates the soil surface roughness
	 * Input:	z0s: 	roughness length of the soil surface
	 * 			d0:		zero-plane displacement [mm]
	 * 			z0_z0t: ratio z0/z0t
	 * 			snowD:	snow depth [mm]
	 * 			thresh: threshold on snow depth to change roughness
	 * 			z0snow:	roughness length over snow [mm]
	 *
	 * Output:	z0_ris:	resulting roughness length
	 * 			d0_ris: resulting zero_plane displacement
	 * 			z0_z0t_ris:	resulting ratio z0/zot
	 * comment: Matteo Dall'Amico, May 2009 */
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
void EnergyFluxes(double Tg, long r, long c, long n, double Tg0, double Qg0, double Tv0, double zmu, double zmT, double z0s, double d0s, double rz0s,
	double z0v, double d0v, double rz0v, double hveg, double v, double Ta, double Qa, double P, double LR, double psi, double sat, double e,
	double fc, double LAI, double Wcrn, double Wcrnmax, double Wcsn, double Wcsnmax, double *dWcrn, double *dWcsn, double *theta, double **soil,
	double *land, double *root, PAR *par, DOUBLEVECTOR *ftcl, DOUBLEVECTOR *rep, double SWin, double LWin, double SWv, double *LW, double *H,
	double *dH_dT, double *E, double *dE_dT, double *LWv, double *Hv, double *LEv, double *Etrans, double *Tv, double *Qv, double *Ts, double *Qs,
	double *Hg0, double *Hg1, double *Eg0, double *Eg1, double *rh, double *rv, double *rc, double *rb, double *rh_ic, double *rv_ic, double *Qg,
	double *u_top, double *decay, double *Locc){
	/* Author:  Stefano Endrizzi  Year:2009
	* Input:	Tg: temperature of the first soil layer
	* 			r: row number
	* 			c: column number
	* 			n: number of snow+glacier layers
	* 			zmu: height of wind velocity sensor with respect to the ground surface [m]
	* 			zmT: height of temperature sensor with respect to the ground surface [m]
	* 			z0s: roughness length of the soil surface [mm]
	* 			d0s: zero-plane displacement for the soil [mm]
	* 			rz0s: as to do with roughness length of soil
	* 			z0v:roughness length for vegetation [mm]
	* 			d0v:zero-plane displacement for vegetation [mm]
	* 			rz0v:has to do with roughness length of vegetation
	* 			hveg: vegetation height
	* 			v:	wind velocity in the point
	*			Ta: air temperature
	*			Qa: air specific humidity
	*			P: air pressure
	*			LR: air lapse rate
	*			psi: psi of the first soil layer
	*			ice: ice content of the first soil layer
	*			e: surface LW emissivity (be it snow or soil)
	*			fc: canopy fraction
	*			LAI: leaf area index
	*			Wc: liquid precipitation intercepted by precipitation [mm]
	*			fwet: ??
	*			fsnow: fraction of canopy covered by snow vector of land cover properties for the specified land cover
	*			theta: vector of adimensional water content
	*			soil: matrix of soil parameters [#_of_soil_properties][#_of_soil_layers] for the given soil type
	*			land: vector of land cover properties for the specified land cover
	*			root: vector of root fraction?
	*			par: parameters structure
	*			ftcl: vector of?
	*			rep: vector of?
	*			SWin: incoming SW
	*			LWin: incoming LW
	*			SWv: SW to the vegetation
	*			cont3:
	* Output:
	*			LW: Net LW to the ground
	*			H: sensible heat flux to the ground
	*			E: Evaporation flux to the ground [Kg/m2]
	*			LWv: LW of vegetation
	*			Hv: sensible heat flux for vegetation
	*			Ev: Evaporation flux for vegetation
	*			Evt:
	*			Tv:
	*			Ts: surface temperature
	*			Qs: surface specific humidity
	*			Hg0: sensible heat flux in the ground at t=t0
	*			Hg1: sensible heat flux in the ground at t1=t0+Dt
	*			Eg0: evaporation heat flux in the ground at t=t0
	*			Eg1: evaporation heat flux in the ground at t1=t0+Dt
	*			rh:
	*			rv:
	*			rc:
	*			rb:
	*			rh_ic:
	*			rv_ic:
	*			Qg:
	* comment: Matteo Dall'Amico, May 2009 */

	double Hg, dHg_dT, Eg, dEg_dT, LWg;
	double dQgdT, rm;

	//initalization
	initialize_doublevector(ftcl,0.0);
	*H=0.0;	*dH_dT=0.0;
	*E=0.0;	*dE_dT=0.0;
	*LW=0.0; *LWv=0.0;
	*Hv=0.0; *LEv=0.0; *Etrans=0.0;
	*Hg0=0.0; *Eg0=0.0;
	*Hg1=0.0; *Eg1=0.0;
	*decay=0.0; *Locc=0.0;

	//thermodynamical calculations
	sat_spec_humidity(Qg, &dQgdT, red_evap(n, psi, Tg), Tg, P);

	if(fc==0){ *Qs=*Qg;  *Ts=Tg;  *Qv=0.0; *Tv=0.0; *rc=1.E99; *rb=1.E99; *rh_ic=1.E99; *rv_ic=1.E99; *u_top=0.0; }

	if(fc<1){
		aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, LR, rep, &rm, rh, rv, par->state_turb, par->monin_obukhov);
		if(*Qg>Qa && n==0) *rv=(*rv)+exp(8.206-4.255*sat);

		turbulent_fluxes(*rh, *rv, P, Ta, Tg, Qa, *Qg, dQgdT, &Hg, &dHg_dT, &Eg, &dEg_dT);

		*H+=(1.0-fc)*Hg;	*dH_dT+=(1.0-fc)*dHg_dT;
		*E+=(1.0-fc)*Eg;	*dE_dT+=(1.0-fc)*dEg_dT;

		*LW+=(1.0-fc)*( e*LWin-SB(Tg) );

		*Hg0=Hg;
		*Eg0=Eg;

		//error messages
		if(Hg!=Hg){
			printf("Hg no value bare soil %ld %ld \n",r,c);
			printf("rh:%e rv:%e P:%f Ta:%f Tg:%f Qa:%f Qg:%f Hg:%f Eg:%f\n",*rh,*rv,P,Ta,Tg,Qa,*Qg,Hg,Eg);
			stop_execution();
		}
		if(Eg!=Eg) printf("Eg no value bare soil %ld %ld \n",r,c);

	}

	//TURBULENT FLUXES FOR CANOPY CANOPY
	if(fc>0){

		Tcanopy(r, c, Tv0, Tg, *Qg, dQgdT, Tg0, Qg0, Ta, Qa, zmu, zmT, z0v, z0s, d0v, rz0v, hveg, v, LR, P, SWin, SWv, LWin, e, LAI, land, Wcrn,
			Wcrnmax, Wcsn, Wcsnmax, dWcrn, dWcsn, LWv, &LWg, Hv, &Hg, &dHg_dT, LEv, &Eg, &dEg_dT, Ts, Qs, root, theta, ftcl->co, rep, par, n,
			sat, rh, rv, rc, rb, rh_ic, rv_ic, u_top, Etrans, Tv, Qv, decay, Locc);

		*H+=fc*Hg;	*dH_dT+=fc*dHg_dT;
		*E+=fc*Eg;	*dE_dT+=fc*dEg_dT;

		*LW+=fc*LWg;

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
void PointEnergyBalance(long r,			long c,			long ns,		long ng,		double zmu,		double zmT,		double z0s,			double d0s,			double rz0s,								double z0v,
						double d0v,		double rz0v,	double hveg,	double v,		double Ta,		double Qa,		double P,					double LR,				double eps,									double fc,
						double LAI,		double *Wcrn,	double Wcrnmax,	double *Wcsn,	double Wcsnmax,	double *theta,	double *land,	double *root,	DOUBLEVECTOR *ftcl,	DOUBLEVECTOR *turb_rep,
						double SWin,	double LWin,	double SWv,		double *LW,		double *H,		double *E,		double *LWv,		double *Hv,			double *LEv,								double *Etrans,
						double *Ts,		double *Qs,		double Hadd,	double *SW,		double *k,		double *D,		double *wi,			double *wl,			double *dw,									double *T,
						double *Hg0,	double *Hg1,	double *Eg0,	double *Eg1,	double *Qv,		double *Qg,		double *rh,			double *rv,			double *rb,									double *rc,
						double *rh_ic,	double *rv_ic,	double *u_top,	SOIL *sl,		PAR *par,		double *decay,	double *Locc,	double TsupN,	double TsupNp1){

	short sy=sl->type->co[r][c];
	long l, cont=0, cont2, n=Nl+ns+ng;
	double dH_dT, dE_dT, EB, dEB_dT, C;
	DOUBLEVECTOR *ad, *ads, *adi, *b, *e, *T0, *T1, *DT, *Tstar, *theta0;
	double Zboundary=land[jzb], Tboundary=land[jtb];
	double error=1.E99, error0, nw=1.0, EB0, DTn, Tg, psim, psiliq;
	double Qg0, dQ, Tv0, dWcsn=0.0, dWcrn=0.0, dU, C0, C1, th0, th1, sink;
	double max_tol;// max tolerance
	short max_iter;// max number of iterations
	short iter_close;
	short Evap; // 1: evaporation is active, 0: evaporation excluded
	if (par->superfast!=1) {// regular version
		Evap=1;
		max_tol=1E-4;
		max_iter=10;
	}else{// when the superfast version is on, then evaporation is excluded by default
		Evap=0;
		max_tol=1E-3;
		max_iter=5;
	}

	//ALLOCATE
	ad=new_doublevector(n);
	ads=new_doublevector(n-1);
	adi=new_doublevector(n-1);
	b=new_doublevector(n);
	e=new_doublevector(n);
	T0=new_doublevector(n);
	T1=new_doublevector(n);
	DT=new_doublevector(n);
	Tstar=new_doublevector(Nl);
	theta0=new_doublevector(Nl);

	//snow+soil
	for(l=1;l<=n;l++){
		T0->co[l]=T[l];
		T1->co[l]=T[l];
		dw[l]=0.0;
		if(l>ns+ng){
			Zboundary-=D[l];
			theta0->co[l-ns-ng]=theta[l-ns-ng];
		}
	}

	Tg=T[1];
	sat_spec_humidity(&Qg0, &dQ, red_evap(ns+ng, sl->P->co[1][r][c], sl->T->co[1][r][c]), Tg, P);
	Tv0=sl->Tv->co[r][c];
	psiliq=sl->P->co[1][r][c];
	if(par->superfast!=1){
		EnergyFluxes(Tg, r, c, ns+ng, T0->co[1], Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, P, LR, psiliq,
			(sl->thice->co[1][r][c]+theta[1])/sl->pa->co[sy][jsat][1], eps, fc, LAI, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, theta,
			sl->pa->co[sy], land, root, par, ftcl, turb_rep, SWin, LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans, &(sl->Tv->co[r][c]),
			Qv, Ts, Qs, Hg0, Hg1, Eg0, Eg1, rh, rv, rc, rb, rh_ic, rv_ic, Qg, u_top, decay, Locc);
		EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
		dEB_dT = -eps*dSB_dT(Tg) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;
	}
	/* Decomment the following if you want to compare with Hansson et al (2004)
	 SW[1]=0.0;
	EB=-28.*(Tg+6);
	dEB_dT=-28.;*/

	EB0 = EB;

	if(Evap!=1){*Etrans=0.0;*E=0.0;}

	//Tstar only changed for soil
	for(l=ns+ng+1;l<=n;l++){
		psim=psi_teta(theta[l-ns-ng]+sl->thice->co[l-ns-ng][r][c], 0.0, sl->pa->co[sy][jsat][l-ns-ng], sl->pa->co[sy][jres][l-ns-ng], sl->pa->co[sy][ja][l-ns-ng], sl->pa->co[sy][jns][l-ns-ng], 1-1/sl->pa->co[sy][jns][l-ns-ng], PSImin, par->Esoil);
		Tstar->co[l-ns-ng]=Fmin(psim/(1000.0*Lf/(g*(Tfreezing+tk))), 0.0);
	}

	do{

		coeff(r, c, ns, ng, ns+1, ad->co, ads->co, adi->co, b->co, k, D, T0->co, T, wl, wi, dw, EB0, EB, dEB_dT, SW, Zboundary, Tboundary, par->Dt,
			(*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par, TsupNp1);

		for(l=1;l<=n;l++){
			if(l<=ns+ng){//snow
				C = Lf*(wi[l]+wl[l])*dtheta_snow(T[l])/D[l];
			}else{	//soil
				if(T[l]<=Tstar->co[l-ns-ng]){
					C = rho_w*Lf*(Lf/(g*tk)*1.E3)*dteta_dpsi(Psif(T[l]), 0.0, sl->pa->co[sy][jsat][l-ns-ng], sl->pa->co[sy][jres][l-ns-ng], sl->pa->co[sy][ja][l-ns-ng], sl->pa->co[sy][jns][l-ns-ng], 1-1/sl->pa->co[sy][jns][l-ns-ng], PSImin, par->Esoil);
				}else{
					C = 0.0;
				}
			}

			C0 = calc_C0(l, r, c, ns+ng, wi, wl, dw, D, sl);
			C1 = calc_C(l, r, c, ns+ng, wi, wl, dw, D, sl);

			dU = Lf*dw[l] + C1*D[l]*T[l] - C0*D[l]*T0->co[l];

			ad->co[l] += (C+C1)*D[l]/par->Dt;
			b->co[l] += ((C+C1)*D[l]*T[l] - dU)/par->Dt;

		}
		if(ad->co[1]==0.0){
				printf("r=%ld, c=%ld, ns+ng=%ld, wi[1]=%f, wl[1]=%f, dw[1]=%f, D[1]=%f C1=%f, C=%f\n",r, c, ns+ng, wi[1], wl[1], dw[1], D[1], C1, C);
				double lambdaT1,lambdaT2,lambdaTfin;
				lambdaT1=calc_k(1, r, c, ns, ng, wi, wl, dw, T, D, *k_thermal_snow_Sturm, *k_thermal_snow_Yen, sl, par);
				lambdaT2=calc_k(2, r, c, ns, ng, wi, wl, dw, T, D, *k_thermal_snow_Sturm, *k_thermal_snow_Yen, sl, par);
				lambdaTfin=(lambdaT1*lambdaT2*0.5*(D[1]+D[2]))/(lambdaT1*0.5*D[2]+lambdaT2*0.5*D[1]);
				printf("\nad[1]=%f, lambdaT1=%f, lambdaT2=%f",ad->co[1],lambdaT1,lambdaT2);
				t_error("Energy balance: ad[1]=0.0");
			}
		tridiag(1, r, c, n, adi, ad, ads, b, e);

		for(l=1;l<=n;l++){
			T[l]=e->co[l];
			DT->co[l]=T[l]-T1->co[l];
		}

		cont2=0;
		nw=1.0;
		error0=error;
		iter_close=0;
		do{

			for(l=1;l<=n;l++){
				T[l]=T1->co[l] + nw*DT->co[l];

				//soil
				if(l>ns+ng){
					if( cont==0 && ( (T0->co[l]<Tstar->co[l-ns-ng] && T[l]>=Tstar->co[l-ns-ng]) || (T0->co[l]>=Tstar->co[l-ns-ng] && T[l]<Tstar->co[l-ns-ng]) ) ){
						T[l]=Fmin(Tstar->co[l-ns-ng], T_max_dteta(sl->pa->co[sy][ja][l-ns-ng], sl->pa->co[sy][jns][l-ns-ng], 1-1/sl->pa->co[sy][jns][l-ns-ng]));
						iter_close=-1;
					}

					/*th0=teta_psi(Psif(Fmin(Tstar->co[l-ns-ng],T0->co[l])), 0.0, sl->pa->co[sy][jsat][l-ns-ng], sl->pa->co[sy][jres][l-ns-ng],
						sl->pa->co[sy][ja][l-ns-ng], sl->pa->co[sy][jns][l-ns-ng], 1-1/sl->pa->co[sy][jns][l-ns-ng], PSImin, par->Esoil);*/
					th0=theta0->co[l-ns-ng];
					th1=teta_psi(Psif(Fmin(Tstar->co[l-ns-ng],T[l])), 0.0, sl->pa->co[sy][jsat][l-ns-ng], sl->pa->co[sy][jres][l-ns-ng],
						sl->pa->co[sy][ja][l-ns-ng], sl->pa->co[sy][jns][l-ns-ng], 1-1/sl->pa->co[sy][jns][l-ns-ng], PSImin, par->Esoil);

					dw[l]=(th1-th0)*D[l]*rho_w;


				//snow
				}else{

					if(T0->co[l]<0 && T[l]>=0){
						T[l]=max_dtheta_snow();
						iter_close=-1;
					}

					//th0=theta_snow(T0->co[l]);
					th0=wl[l]/(wi[l]+wl[l]);
					th1=theta_snow(T[l]);

					dw[l]=(th1-th0)*(wi[l]+wl[l]);

				}

				//printf("l:%ld T0:%f T10:%f T1:%f wl:%f wi:%f dw:%f cont:%ld cont2:%ld\n",l,T0->co[l],T1->co[l],T[l],wl[l],wi[l],dw[l],cont,cont2);

			}

			Tg=T[1];
			for(l=ns+ng+1;l<=n;l++){
				theta[l-ns-ng] = theta0->co[l-ns-ng] + 0.5*dw[l]/(rho_w*D[l]);

				//admit evaporation and evapotranspiration from soil only if theta>=1.1*theta_res
				if(theta[l-ns-ng] >= 1.1*sl->pa->co[sy][jres][l-ns-ng]){
					theta[l-ns-ng] -= 0.5*Fmax( par->Dt*(*Etrans)*ftcl->co[l-ns-ng]/(rho_w*D[l]), 0.0 );
					if(theta[l-ns-ng] < 1.1*sl->pa->co[sy][jres][l-ns-ng]) theta[l-ns-ng] = 1.1*sl->pa->co[sy][jres][l-ns-ng];
				}
			}

			if(ns+ng==0){
				if(theta[1] >= 1.1*sl->pa->co[sy][jres][1]){
					theta[1] -= 0.5*Fmax( par->Dt*(*E)/(rho_w*D[1]), 0.0 );
					if(theta[1] < 1.1*sl->pa->co[sy][jres][1]) theta[1] = 1.1*sl->pa->co[sy][jres][1];
				}
				psiliq=teta_psi(theta[1], 0.0, sl->pa->co[sy][jsat][1], sl->pa->co[sy][jres][1], sl->pa->co[sy][ja][1],
					sl->pa->co[sy][jns][1], 1-1/sl->pa->co[sy][jns][1], PSImin, par->Esoil);
			}
			if(par->superfast!=1){
				EnergyFluxes(Tg, r, c, ns+ng, T0->co[1], Qg0, Tv0, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, P, LR, psiliq,
					(sl->thice->co[1][r][c]+theta[1])/sl->pa->co[sy][jsat][1], eps, fc, LAI, *Wcrn, Wcrnmax, *Wcsn, Wcsnmax, &dWcrn, &dWcsn, theta,
					sl->pa->co[sy], land, root, par, ftcl, turb_rep, SWin, LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, LEv, Etrans, &(sl->Tv->co[r][c]),
					Qv, Ts, Qs, Hg0, Hg1, Eg0, Eg1, rh, rv, rc, rb, rh_ic, rv_ic, Qg, u_top, decay, Locc);
				EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
				dEB_dT = -eps*dSB_dT(Tg) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;
			}
			/*SW[1]=0.0;
			EB=-28.*(Tg+6);
			dEB_dT=-28.;*/

			if(Evap!=1){*Etrans=0.0;*E=0.0;}

			coeff(r, c, ns, ng, ns+1, ad->co, ads->co, adi->co, b->co, k, D, T0->co, T, wl, wi, dw, EB0, EB, dEB_dT, SW, Zboundary, Tboundary, par->Dt,
				(*k_thermal_snow_Sturm), (*k_thermal_snow_Yen), sl, par, TsupNp1);

			for(l=1;l<=n;l++){
				C0 = calc_C0(l, r, c, ns+ng, wi, wl, dw, D, sl);
				C1 = calc_C(l, r, c, ns+ng, wi, wl, dw, D, sl);

				dU = Lf*dw[l] + C1*D[l]*T[l] - C0*D[l]*T0->co[l];

				b->co[l] -= dU/par->Dt;
			}

			error=0.0;
			for(l=1;l<=n;l++){
				if(l==1){
					error+=pow(ad->co[l]*T[l] + ads->co[l]*T[l+1] - b->co[l], 2.0);
				}else if(l>1 && l<n){
					error+=pow(adi->co[l-1]*T[l-1] + ad->co[l]*T[l] + ads->co[l]*T[l+1] - b->co[l], 2.0);
				}else{
					error+=pow(adi->co[l-1]*T[l-1] + ad->co[l]*T[l] - b->co[l], 2.0);
				}
			}
			error=pow(error,0.5);

			cont2++;
			nw/=3.0;

		}while(error>error0 && cont2<5);

		DTn=0.0;
		for(l=1;l<=n;l++){
			DTn+=pow(T[l]-T1->co[l], 2.0);
			T1->co[l]=T[l];
		}
		DTn=pow(DTn,0.5);

		cont++;

		if(iter_close>=0 && DTn<=max_tol) iter_close=1;
		if(cont>=max_iter) iter_close=1;

	}while(iter_close!=1);

	/*if(DTn>1.E-3){
		printf("PointEnergyBalance not converging DTn:%f l:%ld r:%ld c:%ld\n",DTn,l,r,c);
	}*/

	//Update
	for(l=1;l<=Nl+ns+ng;l++){

		//error messages

		if(T[l]!=T[l]) printf("T no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
		if(dw[l]!=dw[l]) printf("dw no value, error 1, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);
		if(T[l]<-80 || T[l]>50) printf("T outside of range, error 2, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f\n",l,r,c,T[l]);

		wl[l]+=dw[l];
		wi[l]-=dw[l];
		wi[l]=Fmax(0.0, wi[l]);
		wl[l]=Fmax(0.0, wl[l]);
		if(T[l]>0) wi[l]=0.0;

		if(l>ns+ng){
			theta[l-ns-ng] = wl[l]/(rho_w*D[l]);
			sl->thice->co[l-ns-ng][r][c] = wi[l]/(rho_w*D[l]);
			if(sl->thice->co[l-ns-ng][r][c]>sl->pa->co[sy][jsat][l-ns-ng]-0.00001) sl->thice->co[l-ns-ng][r][c]=sl->pa->co[sy][jsat][l-ns-ng]-0.00001;
			sl->T->co[l-ns-ng][r][c] = T[l];

			if(l==1 && ns+ng==0){
				sl->ET->co[l-ns-ng][r][c] = Fmax( (*Etrans)*ftcl->co[l-ns-ng], 0.0) + Fmax( (*E), 0.0);
			}else{
				sl->ET->co[l-ns-ng][r][c] = Fmax( (*Etrans)*ftcl->co[l-ns-ng], 0.0);
			}

			if(T[l]!=T[l]) printf("T no value, error 3, PointEnergyBalance, l:%ld r:%ld c:%ld T:%f\n",l,r,c,T[l]);
			if(T[l]<-80 || T[l]>50) printf("T outside of range, error 3, PointEnergyBalance, l:%ld r:%ld c:%ld\n",l,r,c);

		}
	}


	*Wcrn = *Wcrn + dWcrn;
	*Wcsn = *Wcsn + dWcsn;

	//Deallocate
	free_doublevector(ad);
	free_doublevector(ads);
	free_doublevector(adi);
	free_doublevector(b);
	free_doublevector(e);
	free_doublevector(T0);
	free_doublevector(T1);
	free_doublevector(DT);
	free_doublevector(Tstar);
	free_doublevector(theta0);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void coeff(long r, long c, long ns, long ng, long nlim, double *ad, double *ads, double *adi, double *b, double *k0, double *D, double *T0,
	double *T, double *wl, double *wi, double *dw, double EB0, double EB1, double dEB, double *SW, double Zb, double Tb, double Dt,
	double (* kfunct_snow)(double rho), double (* kfunct_glac)(double rho), SOIL *sl, PAR *par, double Ts){

	double Dcorr, k, kn, k00;
	long l, n=Nl+ns+ng;

	//upper boundary condition
	short Neumann;
	if (par->superfast!=1) {// regular version: Neumann condition
		Neumann=1;
	}else {// if the superfast version is on, then by default the top boundary condition is Dirichlet
		Neumann=0;
	}
	//double Ts=-5.0;
	//bottom boundary condition
	short Dirichlet=1;
	double Jbottom=0.0;

	//FIRST LAYER
	l=1;

	Dcorr=D[l];
	//Dcorr=0.5*(0.5*D[1]+CA*(D[1]+0.5*D[2]));

	k = calc_k(l, r, c, ns, ng, wi, wl, dw, T, D, kfunct_snow, kfunct_glac, sl, par);

	if(Neumann!=1){	//Dirichlet
		ad[l] = 2.0*(1.-KNe)*k/(D[l]);
		b[l] = 2.0*(1.-KNe)*k*Ts/(D[l]);
	}else{	//Neumann
		ad[l] = 0.0;
		b[l] = 0.0;
	}

	kn = calc_k(l+1, r, c, ns, ng, wi, wl, dw, T, D, kfunct_snow, kfunct_glac, sl, par);
	k = (k*kn*0.5*(D[l]+D[l+1]))/(k*0.5*D[l+1]+kn*0.5*D[l]);

	ad[l] += 2.0*(1.-KNe)*k/(D[l]+D[l+1]);
	ads[l] = -2.0*(1.-KNe)*k/(D[l]+D[l+1]);
	b[l] += 2.0*KNe*k0[l]*(T0[l+1]-T0[l])/(D[l]+D[l+1]);

	if(Neumann!=1){// Dirichlet
		k00 = calc_k(l, r, c, ns, ng, wi, wl, dw, T0, D, kfunct_snow, kfunct_glac, sl, par);
		b[l] += 2.0*KNe*k00*(Ts-T0[l])/(D[l]);
	}else{
		b[l] += (KNe*EB0 + (1.0-KNe)*(EB1-dEB*T[l]) + SW[l]);
		ad[l] -= (1.0-KNe)*dEB;
	}

	//INTERMEDIATE LAYERS
	for(l=2;l<=n-1;l++){
		adi[l-1] = ads[l-1];	//symmetrical matrix
		ad[l] = 2.0*(1.-KNe)*k/(D[l-1]+D[l]);

		k = calc_k(l, r, c, ns, ng, wi, wl, dw, T, D, kfunct_snow, kfunct_glac, sl, par);
		kn = calc_k(l+1, r, c, ns, ng, wi, wl, dw, T, D, kfunct_snow, kfunct_glac, sl, par);
		k = (k*kn*0.5*(D[l]+D[l+1]))/(k*0.5*D[l+1]+kn*0.5*D[l]);

		ad[l] += 2.0*(1.-KNe)*k/(D[l]+D[l+1]);
		ads[l] = -2.0*(1.-KNe)*k/(D[l]+D[l+1]);
		b[l] = - 2.0*KNe*k0[l-1]*(T0[l]-T0[l-1])/(D[l-1]+D[l]) + 2.0*KNe*k0[l]*(T0[l+1]-T0[l])/(D[l]+D[l+1]);

		if(l<=nlim) b[l]+=SW[l];
	}


	//LAST LAYER
	l=n;

	adi[l-1] = ads[l-1];
	ad[l] = 2.0*(1.-KNe)*k/(D[l-1]+D[l]);

	b[l] = - 2.0*KNe*k0[l-1]*(T0[l]-T0[l-1])/(D[l-1]+D[l]);

	if(Dirichlet==1){
		k = calc_k(l, r, c, ns, ng, wi, wl, dw, T, D, kfunct_snow, kfunct_glac, sl, par);
		if(Zb<0){
			k = 0.0;
			k0[l] = 0.0;
		}
		b[l] += (2.0*KNe*k0[l]*(Tb-T0[l])/(D[l]+2.0*Zb) + 2.0*(1.-KNe)*k*Tb/(D[l]+2.0*Zb));
		ad[l] += 2.0*(1.-KNe)*k/(D[l]+2.0*Zb);
	}else{
		b[l] += Jbottom;
	}

	if(l<=nlim) b[l]+=SW[l];
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void energy_balance_superfast(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, SNOW *snow, GLACIER *glac)

	{

	//1.DICHIARAZIONE DELLE VARIABILI
	//*******************************

	long r;//Counters of rows of the basin
	long c;//Counters of columns of the basin
	long l;//Counters of layers of the basin
	long Ns;// Maximum number of snow layers
	long Ng;// Maximum number of glacier layers
	long ns;//Number of snow layer in a pixel,
	long ng=0;// number of glacier layer in a pixel
	double Ts;//GROUND SURFACE TEMPERATURE
	double Qs;// Specific humidity of the surface
	double SWin; //Short wave radiation incoming [W/mq]
	double SWv_vis;// SW going to vegetation nel visibile
	double	SWv_nir;// SW going to vegetation nel near infrared
	double	LWin; // incoming Long Wave long wave radiation diffused towards land [W/m2] for a pixel
	double	LW;// LW net to the ground
	double	LWv;// Long wave of vegetation
	double Hv; // Sensible heat flux of the vegetation [W/m2]
	double H;//Sensible heat flux of the ground [W/m2] for a pixel
	double	LEv;
	double	E;//Evaporation fluxes [kg/(s*m2)]
	double	Etrans;// transpiration from the vegetation
	double z0;// surface roughness [mm] i.e. the height above the ground where the wind velocity profile is 0, following the rule of thumb it's roughly 1/10 of the height of the object.
	double z0veg;// surface roughness [mm] for the vegetation
	double d0veg;//zero-plane displacement for vegetation
	double hveg;// vegetation height (not present in this version)
	double eps;// surface LW emissivity (be it soil or snow) [-]
	double Qa;// Air specific humidity
	double max_wcan_rain;
	double max_wcan_snow;
	double TsupNp1;// Dirichlte Top boundary condition at the current time step
	double LAI;// leaf area index
	double fc;// canopy fraction
	double	Vpoint;// wind velocity of the point [m/s]
	double Ppoint;
	double Tpoint;
	double zmeas_T; // [m] elevation of the air temperature sensor (with respect to the ground surface)
	double	zmeas_u;// [m] elevation of the wind temperature sensor (with respect to the ground surface)
	double	Hadv;
	double psisat;
	double u_top;
	short sy; // soil type
	short lu;// land use
	double Hg0;
	double Hg1;
	double Eg0;
	double Eg1;
	double rh;// resistenza del calore sensibile
	double rv;// resistenza del calore latente
	double rc;// resistenza della canopy per calore latente
	double rb;// resistenza della canopy per calore sensibile
	double rh_ic;// resistenza del calore sensibile undrcanopy
	double rv_ic;// resistenza del calore latente undercanopy
	double Qv;// umiditˆ specifica vicina alla canopy
	double Qg;//umiditˆ specifica vicina alla superficie del suolo
	double decaycoeff;
	double Locc;
	DOUBLEVECTOR *theta;//Adimensional water content
	DOUBLEVECTOR *k_thermal; //, thermal conductivity [W m^-1 K^-1] for a snow/soil layer (for a pixel)
	DOUBLEVECTOR *k_thermal_interface; // thermal conductivity [W m^-1 K^-1] at the interface between 2 layers (for a pixel)
	DOUBLEVECTOR *D; //vector of column (snow+soil) layer thickness [m]
	DOUBLEVECTOR *wliq; // vector of column (snow+soil)liquid water content [kg/mq]
	DOUBLEVECTOR *wice; // vector of column (snow+soil) ice content [kg/mq]
	DOUBLEVECTOR *Temp; // vector of column (snow+soil) temperature [Celsius]
	DOUBLEVECTOR *deltaw; // vector of column (snow+soil) melting ice (if>0) or icing water (if<0) [kg/mq] in Dt. The 1st component is for the highest layer (be it of snow or soil), the (nl+ns)th component for the deepest soil layer
	DOUBLEVECTOR *turbulence;
	DOUBLEVECTOR *ftcl;
	DOUBLEVECTOR *SWabsf;
	DOUBLEVECTOR *froot;
	DOUBLEVECTOR *DPsi;
	//ALLOCATION
	Ns=par->snowlayer_max; /*maximum number of snow layers*/
	Ng=par->glaclayer_max; /*maximum number of glacier layers*/

	//The 1st component is for the highest layer (be it of snow or glacier or soil), the (Nl+ns)th component for the deepest soil layer
	k_thermal=new_doublevector(Nl+Ns+Ng);
	k_thermal_interface=new_doublevector(Nl+Ns+Ng);
	D=new_doublevector(Nl+Ns+Ng);
	wliq=new_doublevector(Nl+Ns+Ng);
	wice=new_doublevector(Nl+Ns+Ng);
	Temp=new_doublevector(Nl+Ns+Ng);
	deltaw=new_doublevector(Nl+Ns+Ng);
	SWabsf=new_doublevector(Ns+1);
	theta=new_doublevector(Nl);
	ftcl=new_doublevector(Nl);
	froot=new_doublevector(Nl);
	DPsi=new_doublevector(Nl);
	turbulence=new_doublevector(11);

	// INITIALIZATION
	snow->melted_basin=0.0;
	snow->evap_basin=0.0;
	snow->subl_basin=0.0;
	glac->melted_basin=0.0;
	glac->evap_basin=0.0;
	glac->subl_basin=0.0;
	initialize_doublevector(snow->CR1,0.0);
	initialize_doublevector(snow->CR2,0.0);
	initialize_doublevector(snow->CR3,0.0);

	//COMPUTATIONS FOR EACH CELL
	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){

			//printf("%ld %ld\n",r,c);
			if(times->time==0) sl->TsupN=sl->T->co[1][r][c];
			if(land->LC->co[r][c]!=NoV){//if the pixel is not a novalue

				//value used to calculate soil characteristics
				sy=sl->type->co[r][c];
				lu=(short)land->LC->co[r][c];

				//INITIALIZATION
				initialize_doublevector(deltaw,0.0);//ok
				initialize_doublevector(SWabsf,0.0);//ok
				initialize_doublevector(ftcl,0.0);//ok
				initialize_doublevector(froot,0.0);//ok
				initialize_doublevector(turbulence,0.0);//ok
				//snow
				snow->lnum->co[r][c]=0;

				//calculates theta from psi(state variable)
				for(l=1;l<=Nl;l++){
					psisat=psi_saturation(sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l]);
					DPsi->co[l]=Fmax(sl->P->co[l][r][c]-psisat, 0.0);
					sl->P->co[l][r][c]=Fmin(sl->P->co[l][r][c], psisat);
					theta->co[l]=teta_psi(sl->P->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l], sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], PSImin, par->Esoil);
					if(theta->co[l]!=theta->co[l])printf("theta no value l:%ld teta:%f P:%f DPSi:%f T:%f\n",l,theta->co[l],sl->P->co[l][r][c],DPsi->co[l],sl->T->co[l][r][c]);
				}

				ns=0;// ok snow layers
				ng=0;//ok glacier

				//RAIN AND SNOW PRECIPITATION [mm]
				wat->total->co[r][c]=0.0;
				wat->wcan_rain->co[r][c]=0.0;
				wat->wcan_snow->co[r][c]=0.0;
				snow->Psnow->co[r][c]=0.0;

				//SOIL AND SNOW PROPERTIES
				soil_properties(r, c, ns+ng+1, ns+ng+Nl, theta->co, sl->thice->co, D->co, wliq->co, wice->co, k_thermal->co, Temp->co, sl);

				k_interface(ns+ng+Nl, k_thermal->co, D->co, k_thermal_interface->co);

				SWv_vis=0.0;//shortwave radiation absorbed by canopy
				SWv_nir=0.0;//ok
				zmeas_u=1;//ok
				zmeas_T=1;//ok
				z0=1;//ok
				z0veg=1E-6;//ok
				d0veg=1E-6;//ok
				hveg=1;// ok
				TsupNp1=met->var[0][met->column[0][iTsup]];// for Dirichlet Boundary condition at the top of the SOIL column
				Vpoint=1;//ok
				Tpoint=TsupNp1; // ok
				Qa=1;// ok
				Ppoint=1;//ok
				met->LRv[2]=1;// Tair lapse rate
				eps=land->ty->co[lu][jemg];//ok
				fc=0;//ok
				LAI=0.5;//ok
				wat->wcan_rain->co[r][c]=0.0;//ok
				max_wcan_rain=0;//ok
				wat->wcan_snow->co[r][c]=0;//ok
				max_wcan_snow=0;//ok
				SWin=1.0;//ok
				LWin=1.0;//ok
				SWv_vis=1;//ok
				SWv_nir=1;//ok
				LW=1;//ok
				H=1;//ok
				E=1;//ok
				LWv=1;//ok
				Hv=1;//ok
				LEv=1;//ok
				Etrans=1;//ok
				Ts=TsupNp1;//ok
				Qs=1E-6;//ok
				Hadv=1E-6;//ok
				Hg0=1;//ok
				Hg1=1;//ok
				Eg0=1;//ok
				Eg1=1;//ok
				Qv=0;//ok
				Qg=0;//ok
				rh=1.0;//ok
				rv=1.0;//ok
				rb=1.0;//ok
				rc=1.0;//ok
				rh_ic=1.0;//ok
				rv_ic=1.0;//ok
				u_top=1.0;//ok
				decaycoeff=1.0;//ok
				Locc=1.0;//ok

				//ADVECTION
				//if(par->point_sim==0){
				//	if(land->LC2->co[r][c]>0 && snow->Dzl->co[1][r][c]>0){
				//		if(egy->VSFA->co[land->LC2->co[r][c]]!=1.0) Hadv=egy->HSFA->co[land->LC2->co[r][c]]*egy->VSFA->co[land->LC2->co[r][c]]*
				//			adv_efficiency(egy->VSFA->co[land->LC2->co[r][c]])/(1.0-egy->VSFA->co[land->LC2->co[r][c]]);
				//		if(Hadv>egy->HSFA->co[land->LC2->co[r][c]]) Hadv=egy->HSFA->co[land->LC2->co[r][c]];
				//		if(Hadv<0) Hadv=0.0;
				//	}else{
				//		Hadv=0.0;
				//	}
				//}

				//ENERGY BALANCE
				PointEnergyBalance(r,		c,								ns,					ng,								zmeas_u,																		zmeas_T,					z0,															0.0,									0.0,								z0veg,
								   d0veg,	1.0,							hveg,				Vpoint,							Tpoint,																			Qa,										Ppoint,											met->LRv[2],	eps,								fc,
								   LAI,		&(wat->wcan_rain->co[r][c]),	max_wcan_rain,		&(wat->wcan_snow->co[r][c]),	max_wcan_snow,												theta->co,			land->ty->co[lu],	froot->co,			ftcl,							turbulence,
								   SWin,	LWin, 							SWv_vis+SWv_nir,	&LW,							&H,																							&E,										&LWv,													&Hv,									&LEv,							&Etrans,
								   &Ts,		&Qs,							Hadv,				SWabsf->co,						k_thermal_interface->co,		D->co,							wice->co,									wliq->co,				deltaw->co, Temp->co,
								   &Hg0,	&Hg1,							&Eg0,				&Eg1,							&Qv,																						&Qg,									&rh,														&rv,									&rb,								&rc,
								   &rh_ic,	&rv_ic,							&u_top,				sl,								par,																						&decaycoeff,	&Locc,												sl->TsupN,			TsupNp1);

				sl->TsupN=TsupNp1;// memorize the top Dirichlet boundary condition

				//NET PRECIPITATION
				wat->Pn->co[r][c] = 0.0;
				wat->wcan_snow->co[r][c]=0.0;
				//Er_soil=0.0;
				for(l=1;l<=Nl;l++){
					sl->P->co[l][r][c]+=DPsi->co[l];
					if(sl->P->co[l][r][c]!=sl->P->co[l][r][c])printf("Psi no value l:%ld teta:%f P:%f T:%f\n",l,theta->co[l],sl->P->co[l][r][c],sl->T->co[l][r][c]);
				}

			}
		}
	}

	//DEALLOCATION
	free_doublevector(k_thermal);
	free_doublevector(k_thermal_interface);
	free_doublevector(D);
	free_doublevector(wliq);
	free_doublevector(wice);
	free_doublevector(Temp);
	free_doublevector(deltaw);
	free_doublevector(theta);
	free_doublevector(ftcl);
	free_doublevector(turbulence);
	free_doublevector(SWabsf);
	free_doublevector(froot);
	free_doublevector(DPsi);

}//end of "superfast energy_balance" subroutine

