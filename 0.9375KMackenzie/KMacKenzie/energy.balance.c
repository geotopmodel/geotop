
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

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
#include "liston.h"
#include "energy.balance.h"
#include "write_dem.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "meteo.09375.h"
#include "snow.09375.h"
#include "pedo.funct.h"
#include "util_math.h"
#include "rw_maps.h"
#include "times.h"
#include "tabs.h"
#include "output.09375.h"
#include "extensions.h"

void checkErrorSize(char *errfilepath);

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void energy_balance(TIMES *times, PAR *par,	LAND *land, TOPO *top, SOIL *sl, METEO *met, WATER *wat, ENERGY *egy, SNOW *snow, GLACIER *glac, LISTON *liston)

{

//1.DICHIARAZIONE DELLE VARIABILI
//*******************************

long r,c,l,j;//Counters of rows, columns, layers of the basin
long Ns,Ng;//Number of rows, columns,layers of the soil in the whole basin, Maximum number of snow layers, Maximum number of glacier layers
long ns,ng=0;//Number of snow layer in a pixel, number of glacier layer in a pixel
double E0;//soil-earth correction [-]
double Ts;//GROUND SURFACE TEMPERATURE
double Qs;// Specific humidity of the surface
double snowD, glacD=0.0;//Snow depth [mm] and glacier depth [mm]
double Prain_over; //Precipitation as rain [mm] over the canopy
double Psnow_over;//Precipitation as snow [mm water equivalent] over the canopy
double Psoil; //precipitation that reaches the soil [mm] below canopy
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
double	LWv;// Long wave of
double	Tsurr;// Temperature of the sorrounding surfaces (not of the pixel) [ûC]

double Hv; // Sensible heat flux of the vegetation [W/m2]
double H;//Sensible heat fluxe of the groudn [W/m2] for a pixel

double LE;// Latent heat flux
double	Ev; // evaporation fluxes for vegetation
double	E;//Evaporation fluxes [kg/(s*m2)]
double	Evt;// transpiration from the vegetation
double	fwet;
double	fsnow;
double	fsnowcan;// % of canopy covered by snow

double surfEB;//heat flux on the surface (be it snow or soil) [W/m2]
double G;// Heat flux going into the soil surface [W/m2]

DOUBLEVECTOR *theta;//Adimensional water content

DOUBLEVECTOR *c_heat; //Thermal capacity [J m^-3 K^-1] (for a pixel)
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
double z0,
	z0veg,
	d0,
	d0veg,
	z0_z0t,
	hveg,
	eps;
double	ee=0.0;// saturated water vapour pressure
double	dee=0.0;// derivative of saturated water pressure
double	Tdew=0.0; // dew temperature (temp. di rugiada)
double	Qa;// specific humidity
double evap_c;
double	drip_c;
double	maxstorage_c;
double	LAI=0.5;// leaf area index
double fc=0.0;// canopy fraction
double drip_sn; //unused variable
double	expon;// unused variable
double epsa_min;
double	epsa_max;
double	tau_atm;// atmospheri transmittance
double	RHpoint;// Relative humidity of the point [%]
double	Vpoint;// wind velocity of the point [m/s]
DOUBLEVECTOR *turbulence;
DOUBLEVECTOR *ftcl, *SWabsf, *froot;
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
double	tau_cloud;// clear sky trasmissivity
double	tau_cloud_av;//average clear sky trasmissivity
double	sa;// sin(alpha) where alpha is the solar altitude angle
double	Hadv=0.0;
short sy; // soil type
short lu;// land use

double Hg0,
	Hg1,
	Eg0,
	Eg1,
	dE_dT,
	dH_dT;
long cont;
double rh,
	rv,
	rc,
	rb,
	rh_ic,
	rv_ic,
	Qv,
	Qg;
FILE *f;// file pointer

//ALLOCATION
Ns=par->snowlayer_max; /*maximum number of snow layers*/
Ng=par->glaclayer_max; /*maximum number of glacier layers*/

//The 1st component is for the highest layer (be it of snow or glacier or soil), the (Nl+ns)th component for the deepest soil layer
c_heat=new_doublevector(Nl+Ns+Ng);
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
turbulence=new_doublevector(11);
// INITIALIZATION
initialize_doublematrix(snow->rho_newsnow,UV->V->co[2]);
snow->melted_basin=0.0;
snow->evap_basin=0.0;
snow->subl_basin=0.0;
glac->melted_basin=0.0;
glac->evap_basin=0.0;
glac->subl_basin=0.0;
initialize_doublevector(snow->CR1,0.0);
initialize_doublevector(snow->CR2,0.0);
initialize_doublevector(snow->CR3,0.0);

//SHADOW & CLOUDINESS
sun((double)(times->hh+times->mm/60.0), times->JD, &(egy->hsun), &(egy->dsun), &E0, par->latitude, par->longitude, par->ST);
shadow_n(par->point_sim, top, egy->hsun, egy->dsun, land->shadow);
//shadows calculated only when shortwave is not calculated with micromet
//if(par->micromet2!=1) shadow_n(par->point_sim, top, egy->hsun, egy->dsun, land->shadow);

if( par->micromet2==1 ){/* Use Micromet for SWin (=1) */

	if(times->time==0) printf("\nCloudiness data from MICROMET\n\n");

	if(times->time==0){
		f=t_fopen(join_strings(WORKING_DIRECTORY,"ii_clouds.txt"),"w");
		//fprintf(f,"Daily cloud cover inferred by shortwave radiation measurements\n");
		fprintf(f,"time;date;alpha;direction;sin(alpha);fcloud;tau_cloud\n");
		t_fclose(f);
	}

}else{ /* Micromet for SW switched off */

	if( met->column[met->nstcloud-1][iC]!=-1 || met->column[met->nstsrad-1][itauC]!=-1){ /* at least one station has data on cloudiness or SW*/

		if(times->time==0) printf("\nCloudiness data AVAILABLE\n\n");

		if(met->column[met->nstcloud-1][itauC]!=-1){/* cloudiness file available */
			tau_cloud_av=met->var[met->nstcloud-1][met->column[met->nstcloud-1][itauC]];
		}else{
			fcloud=met->var[met->nstcloud-1][met->column[met->nstcloud-1][iC]];
			if(fcloud>1) fcloud=1.0;
			tau_cloud_av=1.0-0.75*pow(fcloud,3.4);
		}

	}else{/* no station has data on cloudiness or SW*/

		fcloud=0.5;	//default value in case no data of cloudiness are available
		tau_cloud_av=1.0-0.75*pow(fcloud,3.4);

		if(times->time==0){
			printf("Warning: Cloudiness data are not provided. Cloud fraction always set at 0.5!!!!\n");
			f=fopen(files->co[ferr]+1,"a");
			fprintf(f,"Warning: Cloudiness data are not provided. Cloud fraction always set at 0.5!!!!\n");
			fclose(f);
		}

	}

	if( (met->column[met->nstsrad-1][iSW]!=-1 || (met->column[met->nstsrad-1][iSWb]!=-1 && met->column[met->nstsrad-1][iSWd]!=-1))  ){
		if(egy->hsun>0) find_tau_cloud_station(met->nstcloud, met, par, egy->hsun, E0, Asurr, met->st->sky->co[met->nstsrad], &tau_cloud, &sa);
	}else{
		tau_cloud=tau_cloud_av;
	}
}


//COMPUTATIONS FOR EACH CELL
for(r=1;r<=Nr;r++){
	for(c=1;c<=Nc;c++){

        if(land->LC->co[r][c]!=NoV){//if the pixel is not a novalue

			//value used to calculate soil characteristics
			sy=sl->type->co[r][c];
			lu=(short)land->LC->co[r][c];

			//INITIALIZATION
			//Sr=sublimation rate Er=evaporation rate Mr=melting rate [mm/s]
			Sr_snow=0.0; Er_snow=0.0; Mr_snow=0.0; Sr_glac=0.0;	Er_glac=0.0; Mr_glac=0.0; Er_soil=0.0; Sr_soil=0.0;
			initialize_doublevector(deltaw,0.0);

			//initial condition
			if(times->time==0.0){

				//snow layer resetting
				snow_layer_combination(r, c, snow, met->Tgrid->co[r][c], par->snowlayer_inf, par->Dmin, par->Dmax, times->time);

				if(par->recover!=1){
					//initial non-dimensional snow age
					non_dimensionalize_snowage(&(snow->age->co[r][c]), met->Tgrid->co[r][c]);

					//glacier
					if(par->glaclayer_max>0) glacier_init_t0(r, c, met->Tgrid->co[r][c], glac, snow, par, times->time);
				}

				for(j=1;j<=par->rc->nrh;j++){
					if(r==par->rc->co[j][1] && c==par->rc->co[j][2])
						write_soil_output(0, j, times->time, 0.0, par->year0, par->JD0, par->rc, sl, par->psimin, par->Esoil);
				}
			}

			//snow
			ns=snow->lnum->co[r][c];
			snowD=DEPTH(r, c, snow->lnum, snow->Dzl);
			if(snowD!=snowD){
				printf("Novalue in egy balance(a): r:%ld c:%ld SnowD:%f lnum:%ld\n",r,c,DEPTH(r, c, snow->lnum, snow->Dzl),snow->lnum->co[r][c]);
				write_snow_all(r, c, snow);
				stop_execution();
			}


			//calculates theta from psi(state variable)
			for(l=1;l<=Nl;l++){
				//printf("a) l:%ld P:%f ice:%f\n",l,sl->P->co[l][r][c],sl->thice->co[l][r][c]);
				theta->co[l]=teta_psi(sl->P->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
					sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);
				//printf("b) l:%ld P:%f th:%f \n",l,sl->P->co[l][r][c],theta->co[l]);
			}

			//glacier
			if(par->glaclayer_max>0){
				ng=glac->lnum->co[r][c];
				glacD=DEPTH(r, c, glac->lnum, glac->Dzl);
			}


			//METEOROLOGICAL COMPUTATIONS

			//temperature and wind velocity measurement heights
			zmeas_u=met->st->Vheight->co[1]-1.0E-3*snowD;
			zmeas_T=met->st->Theight->co[1]-1.0E-3*snowD;
			//zmeas_u=met->st->Vheight->co[1];
			//zmeas_T=met->st->Theight->co[1];
			if(zmeas_u<0.1) zmeas_u=0.1;
			if(zmeas_T<0.1) zmeas_T=0.1;

			//RH and V
			if(par->micromet1==1){
				RHpoint=met->RHgrid->co[r][c];//+0.15;
				if(RHpoint>1) RHpoint=1.0;
				Vpoint=met->Vgrid->co[r][c];
			}else{
				RHpoint=met->RH;
				Vpoint=met->V;
			}

			if(Vpoint<par->Vmin) Vpoint=par->Vmin;
			if(RHpoint<0.01*par->RHmin) RHpoint=0.01*par->RHmin;
			snow->rho_newsnow->co[r][c]=rho_newlyfallensnow(Vpoint, met->Tgrid->co[r][c], Tfreezing);


			//RAIN AND SNOW PRECIPITATION [mm]
			//calculate dew temperature (otherwise replace Tdew with met->Tgrid) to distinguish between rain and snow
			sat_vap_pressure(&ee,&dee,met->Tgrid->co[r][c],met->Pgrid->co[r][c]);
			Qa=RHpoint*spec_humidity(ee, met->Pgrid->co[r][c]);
			ee=Qa*met->Pgrid->co[r][c]/(0.622+Qa*0.378);
			sat_vap_pressure_inv(&Tdew,ee,met->Pgrid->co[r][c]);
			//convert total precipitation to [mm]
			wat->total->co[r][c]*=(par->Dt/3600.0);	//from [mm/h] to [mm]
			//distinguish between rain and snow
			part_snow(wat->total->co[r][c], &Prain_over, &Psnow_over, Tdew, par->T_rain, par->T_snow);
			//modify rain and snow using correction factors
			Prain_over*=par->raincorrfact;
			Psnow_over*=par->snowcorrfact;
			wat->total->co[r][c]=Prain_over+Psnow_over;
			//precipitation reaching soil below canopy
			Psoil=(1.0-fc)*wat->total->co[r][c];
			//rough partition of precipitation reaching the sl (below canopy) in rain and snow
			part_snow(Psoil, &Prain, &Psnow, Tdew, par->T_rain, par->T_snow);
			wat->Psnow->co[r][c]=Psnow;


			//SHORTWAVE RADIATION
			//initialization of shortwave radiation absorbed by soil
			SW=0.0;

			if(par->micromet2==1){/* Use Micromet for SWin */
				fcloud=met->CFgrid->co[r][c];
				tau_cloud_av=1.0-0.75*pow(fcloud,3.4);
				tau_cloud=tau_cloud_av;
				if(r==par->rc->co[1][1] && c==par->rc->co[1][2]){
					f=fopen(join_strings(WORKING_DIRECTORY,"ii_clouds.txt"),"a");
					if(times->mm>=10){
						fprintf(f,"%f;%ld/%ld/%ld %ld:%ld;%f;%f;%f;%f;%f\n",times->time,times->DD,times->MM,times->AAAA,times->hh,times->mm,egy->hsun*180.0/Pi,egy->dsun*180.0/Pi,
							sin(egy->hsun),fcloud,tau_cloud_av);
					}else{
						fprintf(f,"%f;%ld/%ld/%ld %ld:0%ld;%f;%f;%f;%f;%f\n",times->time,times->DD,times->MM,times->AAAA,times->hh,times->mm,egy->hsun*180.0/Pi,egy->dsun*180.0/Pi,
							sin(egy->hsun),fcloud,tau_cloud_av);
					}
					fclose(f);
				}
			}

			//calculation of SWin
			tau_atm=atm_transmittance(egy->hsun, met->Pgrid->co[r][c], RHpoint, met->Tgrid->co[r][c], Asurr, par->Vis, par->Lozone);
			shortwave_radiation(r, c, egy->hsun, egy->dsun, E0, land->shadow->co[r][c], top->sky->co[r][c], tau_cloud, sa, top->slopes->co[r][c],
				top->aspect->co[r][c], tau_atm, met->var[met->nstsrad-1], met->column[met->nstsrad-1], met->st->sky->co[met->nstsrad],  Asurr, &SWbeam, &SWdiff,
				&cosinc, egy->nDt_shadow, egy->nDt_sun);
			SWin=SWbeam+SWdiff;

			//albedo
			if(snowD>0){
				update_snow_age(Psnow_over, snow->T->co[ns][r][c], par->Dt, &(snow->age->co[r][c]));
				avis_b=snow_albedo(land->ty->co[lu][ja_vis], snowD, par->aep, par->avo, 0.2, snow->age->co[r][c], cosinc, (*Fzen));
				avis_d=snow_albedo(land->ty->co[lu][ja_vis], snowD, par->aep, par->avo, 0.2, snow->age->co[r][c], cosinc, (*Zero));
				anir_b=snow_albedo(land->ty->co[lu][ja_nir], snowD, par->aep, par->airo, 0.5, snow->age->co[r][c], cosinc, (*Fzen));
				anir_d=snow_albedo(land->ty->co[lu][ja_nir], snowD, par->aep, par->airo, 0.5, snow->age->co[r][c], cosinc, (*Zero));
			}else{
				avis_b=land->ty->co[lu][ja_vis];
				avis_d=land->ty->co[lu][ja_vis];
				anir_b=land->ty->co[lu][ja_nir];
				anir_d=land->ty->co[lu][ja_nir];
			}
			//shortwave absorbed by soil (when vegetation is not present)
			SW+=( SWdiff*(1.0-0.5*avis_d-0.5*anir_d) + SWbeam*(1.0-0.5*avis_b-0.5*anir_b) );

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
			if(par->micromet3==1){/* Use Micromet for LWin */
				epsa=UV->V->co[2];
				epsa_max=UV->V->co[2];
				epsa_min=UV->V->co[2];
				LWin=egy->LWin->co[r][c];
			}else{
				longwave_radiation(par->state_lwrad, ee, RHpoint, met->Tgrid->co[r][c], tau_cloud_av, &epsa, &epsa_max, &epsa_min);
				Tsurr=surface(r, c, ns, ng, snow->T, glac->T, sl->T);	//Temperature of surrounding surfaces
				LWin=top->sky->co[r][c]*epsa*SB(met->Tgrid->co[r][c])+(1-top->sky->co[r][c])*eps*SB(Tsurr);
			}
			//if incoming LW data are available, they are used (priority)
			LWin=flux(met->nstlrad, iLWi, met->column, met->var, 1.0, LWin);

			//roughness lengths
			update_roughness_soil(land->ty->co[lu][jz0soil], 0.0, 0.0, snowD, land->ty->co[lu][jz0thressoil], par->z0_snow, &z0, &d0, &z0_z0t);

			//SOIL AND SNOW PROPERTIES
			soil_properties(r, c, ns+ng+1, ns+ng+Nl, theta->co, sl->thice->co, D->co, wliq->co, wice->co, k_thermal->co, c_heat->co, Temp->co, sl);

			if(ns>0) snow_properties(r, c, 1, ns, D->co, wliq->co, wice->co, k_thermal->co, c_heat->co, Temp->co, snow->Dzl->co, snow->w_liq->co,
							snow->w_ice->co, snow->T->co, (*k_thermal_snow_Sturm));

			if(ng>0) snow_properties(r, c, ns+1, ns+ng, D->co, wliq->co, wice->co, k_thermal->co, c_heat->co, Temp->co, glac->Dzl->co, glac->w_liq->co,
							glac->w_ice->co, glac->T->co, (*k_thermal_snow_Yen));

			k_interface(ns+ng+Nl, k_thermal->co, D->co, k_thermal_interface->co);

			//ENERGY BALANCE
			PointEnergyBalance(r, c, ns, ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg, d0veg, 1.0, hveg, Vpoint, met->Tgrid->co[r][c], Qa, met->Pgrid->co[r][c], met->LapseRate,
				eps, fc, LAI, wat->wt->co[r][c], fwet, fsnowcan, theta->co, land->ty->co[lu], froot->co, par, ftcl, turbulence, SWin, LWin, SWv_vis+SWv_nir, &LW, &H, &E,
				&LWv, &Hv, &Ev, &Evt, &(sl->Tv->co[r][c]), &Ts, &Qs, sy, sl, Hadv, SWabsf->co, times->time, par->Dt, k_thermal->co, c_heat->co, D->co, wice->co, wliq->co,
				snow->type->co[r][c], deltaw->co, Temp->co, &Hg0, &Hg1, &Eg0, &Eg1);

			rb=0.0; rc=0.0; rh_ic=0.0; rv_ic=0.0;
			EnergyFluxes(Temp->co[1], r, c, ns+ng, zmeas_u, zmeas_T, z0, 0.0, 0.0, z0veg, d0veg, 1.0, hveg, Vpoint, met->Tgrid->co[r][c], Qa, met->Pgrid->co[r][c], met->LapseRate,
				sl->P->co[1][r][c], sl->thice->co[1][r][c], eps, fc, LAI, wat->wt->co[r][c], fwet, fsnowcan, theta->co, sl->pa->co[sy], land->ty->co[lu], froot->co, par, ftcl,
				turbulence, SWin, LWin, SWv_vis+SWv_nir, &LW, &H, &dH_dT, &E, &dE_dT, &LWv, &Hv, &Ev, &Evt, &(sl->Tv->co[r][c]), &Qv, &Ts, &Qs, cont, &Hg0, &Hg1, &Eg0, &Eg1, &rh,
				&rv, &rc, &rb, &rh_ic, &rv_ic, &Qg);

			LE=E*latent(Temp->co[1],Levap(Temp->co[1]));
			surfEB=SW + LW  - H - LE;

			if(ns==0){
				G=surfEB;
			}else{
				G=k_thermal_interface->co[ns]*(Temp->co[ns]-Temp->co[ns+1])/(0.5*D->co[ns]+0.5*D->co[ns+1]);
			}

			//UPDATE SNOW AND SOIL PROPERTIES

			for(l=1;l<=Nl+ns+ng;l++){
				if(wliq->co[l]!=wliq->co[l]) printf("%ld %ld wliq error 1\n",r,c);
				if(wice->co[l]!=wice->co[l]) printf("%ld %ld wice error 1\n",r,c);
			}

			if(E!=E) printf("E error 2\n");
			evaporation(r,c, ns+ng, E, wice->co, wliq->co, D->co, par->Dt, &Sr_snow, &Er_snow, &Er_soil);
			if(ng>0){
				Sr_glac=Sr_snow;
				Er_glac=Er_snow;
				Sr_snow=0.0;
				Er_snow=0.0;
			}

			for(l=1;l<=Nl+ns+ng;l++){
				if(wliq->co[l]!=wliq->co[l]) printf("%ld %ld wliq error 2\n");
				if(wice->co[l]!=wice->co[l]) printf("%ld %ld wice error 2\n");
			}

			update_soil(ng+ns+Nl, r, c, sl, Temp->co, deltaw->co, wice->co, theta->co, ftcl->co, Evt, Er_soil, par);

			liqWBsnow(r, c, snow, &Mr_snow, &Prain_on_snow, par, top->slopes->co[r][c], Prain, wice->co, wliq->co, Temp->co);
			iceWBsnow(r, c, snow, Psnow, met->Tgrid->co[r][c]);
			snow_layer_combination(r, c, snow, met->Tgrid->co[r][c], par->snowlayer_inf, par->Dmin, par->Dmax, times->time);

			if(par->glaclayer_max>0){
				WBglacier(snow->lnum->co[r][c], r, c, glac, &Mr_glac, par, wice->co, wliq->co, Temp->co);
				glac2snow(r, c, snow, glac, par->Dmin, par->Dmax);
				if(par->glaclayer_max>1)glac_layer_combination(r, c, glac, met->Tgrid->co[r][c], Ng, par->Dmin_glac, par->Dmax_glac, times->time);
				ng=glac->lnum->co[r][c];
			}

			snowD=DEPTH(r, c, snow->lnum, snow->Dzl);
			if(snowD!=snowD){
				printf("Novalue in egy balance: r:%ld c:%ld SnowD:%f lnum:%ld\n",r,c,DEPTH(r, c, snow->lnum, snow->Dzl),snow->lnum->co[r][c]);
				write_snow_all(r, c, snow);
				stop_execution();
			}


			//NET PRECIPITATION
			wat->Pn->co[r][c] = Mr_glac + Prain/par->Dt + Mr_snow;


			//OUTPUT DATA
			prepare_output(Er_soil, Mr_snow, Er_snow, Sr_snow, Mr_glac, Er_glac, Sr_glac, Prain, Psnow_over, egy, wat, snow, glac, land, top, sl, met, times,
				par, r, c, 0.25*(avis_d+anir_d+avis_b+anir_b), LE, surfEB, H, G, Temp->co[1], SWin, SWin-SW, SWbeam, eps, LWin, cosinc);

			output_pixel(r, c, Psnow, Prain-Prain_on_snow, Prain_on_snow, Sr_soil, Er_soil, Mr_snow, Er_snow, Sr_snow, Mr_glac, Er_glac, Sr_glac, Evt, LE, H,
				surfEB, G, Temp->co[1], 0.25*(avis_d+anir_d+avis_b+anir_b), eps, LAI, LWin, SWin, LWin-LW, SWin-SW, epsa, epsa_min, epsa_max, SWbeam, SWdiff, DTcorr, Tdew,
				times->n_pixel, turbulence, par, wat, egy, top, met, snow, glac, land, Vpoint, RHpoint, Psnow_over, Prain_over, maxstorage_c, evap_c, drip_c, z0, d0,
				(SWv_vis+SWv_nir), LWv, fc*Hv, fc*Levap(sl->Tv->co[r][c])*Ev, sl->Tv->co[r][c], Ts, Hg0, Hg1, Eg0, Eg1, fc, rh, rv, rb, rc, rh_ic, rv_ic, Hv,
				Levap(sl->Tv->co[r][c])*Ev, Qv, Qg, Qa, Qs);

			output_basin(times->n_basin, Prain, Psnow, met->Tgrid->co[r][c], Temp->co[1], sl->Tv->co[r][c], Er_soil+Sr_soil, Evt, Ev-Evt, LE, H, SW, LW, fc*Levap(sl->Tv->co[r][c])*Ev,
				fc*Hv, fc*(SWv_vis+SWv_nir), fc*LWv, SWin, wat->out2->co, egy->out2->co, times->Dt);

			if(par->ES_num>0) output_altrank(par->ES_num, Prain, Psnow, Temp->co[1], Er_soil, Sr_soil, Evt, Ev-Evt, LE, H, surfEB, SWin, SWin-SW, LWin, eps,
				top->sky->co[r][c], wat->Pn->co[r][c], wat->Pn->co[r][c]-sl->Jinf->co[r][c], Er_snow, Sr_snow, Mr_snow, Er_glac, Sr_glac, Mr_glac, par->Dt, times->n_basin,\
				top->Z0->co[r][c], top->Zmin, top->Zmax, glacD, par->glac_thr, egy->out3->co);

			output_map_plots(r, c, par, times->time, times->n_plot, egy, met, snow, H, LE, fc*Hv, fc*Levap(sl->Tv->co[r][c])*Ev, SWin, SW, fc*(SWv_vis+SWv_nir), LWin, LW, fc*LWv, Ts,
				Temp->co[1], sl->Tv->co[r][c]);

			egy->Hgrid->co[r][c]=H+fc*Hv;
			egy->Tsgrid->co[r][c]=Temp->co[1];

		}
	}
}

//PREPARE SNOW OUTPUT
output_snow(snow, land->LC->co, par);

//DEALLOCATION
free_doublevector(c_heat);
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

}
//end of "energy_balance" subroutine







/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
// C in [J m^-3 K^-1]

void soil_freezing(long l, long r, long c, DOUBLEVECTOR *e0, DOUBLEVECTOR *e1, DOUBLEVECTOR *dw, DOUBLEVECTOR *T, DOUBLEVECTOR *C, SOIL *sl, double psimin){

	double P,Q,th,th1,thi,thi1,theq,A,Tin,ct,Teq,d,Cs,sat,res,a,n,m,h;
	long ri=0, ci=0;
	short sy=sl->type->co[r][c];
	//FILE *f;

	d=sl->pa->co[sy][jdz][l];
	Cs=sl->pa->co[sy][jct][l];
	sat=sl->pa->co[sy][jsat][l];
	res=sl->pa->co[sy][jres][l];
	a=sl->pa->co[sy][ja][l];
	n=sl->pa->co[sy][jns][l];
	m=1-1/n;

	if( (e0->co[l]>Tfreezing && e1->co[l]<Tfreezing) || e0->co[l]<=Tfreezing){

		A=rho_w*0.001*d;		//kg m^-2
		ct=C->co[l]*0.001*d;	//J K-1 m-2
		Tin=fmin(e0->co[l]-Tfreezing, Tfreezing);
		thi=sl->thice->co[l][r][c];
		P=fmin(sl->P->co[l][r][c], psi_saturation(thi, sat, res, a, n, m));
	    theq=teta_psi(Psif2(Tin), 0.0, sat, res, a, n, m, psimin, 1.0);
		th=teta_psi(P, thi, sat, res, a, n, m, fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], psimin)), 1.0);

		if(r==ri && c==ci) printf("l:%ld r:%ld c:%ld thi:%f th:%f P:%f \n",l,r,c,thi,th,P);

		//heat gained(+) or lost(-)
		Q=(e1->co[l]-Tin)*ct;							//K * J * m^-2 * K^-1 = J * m^-2

		//if heat lost
		if(Q<=0){

			if(th>=theq){
				if(Lf*(theq-th)*A<Q){						//J kg^-1 kg m^-2 = j m^-2
					dw->co[l]=Q/Lf;							//J m^-2 J^-1 kg
					T->co[l]=e0->co[l];
					if(r==ri && c==ci) printf("l:%ld case1 Q:%f e0:%f e1:%f T:%f dw:%f wliq:%f wice:%f th:%f theq:%f \n",l,Q,e0->co[l],e1->co[l],T->co[l],dw->co[l],th*A,thi*A,th,theq);
				}else{
					//heat to reach the curve from above
					dw->co[l]=(theq-th)*A;	//<0
					if(r==ri && c==ci) printf("case2 :Q:%f e1:%f e0:%f ct:%f dw:%f th:%f theq:%f qadd:%f\n",Q,e1->co[l],e0->co[l],ct,dw->co[l],th,theq,Lf*dw->co[l]);
					Q-=(Lf*dw->co[l]);
					//move on the curve
					th1=th+dw->co[l]/A;
					thi1=thi-dw->co[l]/A;
					h=internal_energy_soil(th1, thi1, Tin, d, Cs, sat);
					if(r==ri && c==ci) printf("l:%ld case2a) theq:%e thi:%e th:%e h:%e Q:%e dw:%e T:%e\n",l,theq,thi1,th1,h,Q,dw->co[l],Tin);
					from_internal_soil_energy(r, c, l, h+Q, &th1, &thi1, &(T->co[l]), sl->pa->co[sy], psimin);
					dw->co[l]+=(th1-theq)*A;
					if(r==ri && c==ci) printf("l:%ld case2b) thi:%e th:%e T:%f dw:%e theq:%e\n",l,thi1,th1,T->co[l],dw->co[l],teta_psi(Psif2(T->co[l]), 0.0, sat, res, a, n, m, psimin, 1.0));
				}

			}else{
				Teq=P*(g*(tk+Tfreezing))/(1.0E3*Lf);
				if(e0->co[l]+Q/ct >= Teq){
					dw->co[l]=0.0;
					T->co[l]=e1->co[l];
					if(r==ri && c==ci) printf("l:%ld case3 Q:%f e0:%f e1:%f T:%f dw:%f wliq:%f wice:%f th:%f theq:%f \n",l,Q,e0->co[l],e1->co[l],T->co[l],dw->co[l],th*A,thi*A,th,theq);
				}else{
					Q-=ct*(Teq-e0->co[l]);
					//move on the curve
					th1=th;
					thi1=thi;
					h=internal_energy_soil(th1, thi1, Teq, d, Cs, sat);
					if(r==ri && c==ci) printf("l:%ld case4a) e0:%f e1:%f ct(e1-e0):%f theq:%f thi:%f th:%f h:%f Q:%f dw:%f T:%f Teq:%f\n",l,e0->co[l],e1->co[l],ct*(e1->co[l]-e0->co[l]),theq,thi1,th1,h,Q,dw->co[l],Tin,Teq);
					from_internal_soil_energy(r, c, l, h+Q, &th1, &thi1, &(T->co[l]), sl->pa->co[sy], psimin);
					dw->co[l]=(th1-th)*A;
					if(r==ri && c==ci) printf("l:%ld case4b) thi:%f th:%f T:%f dw:%f theq:%f\n",l,thi1,th1,T->co[l],dw->co[l],teta_psi(Psif2(T->co[l]), 0.0, sat, res, a, n, m, psimin, 1.0));
				}
			}

		//if heat gained
		}else{

			if(th<=theq){

				//there is not ice enough for th to reach theq
				if(thi<theq-th){
					if(Lf*thi*A<Q){	//all ice is melted
						dw->co[l]=thi*A;
						T->co[l]=e0->co[l]+(Q-Lf*dw->co[l])/ct;
						if(r==ri && c==ci) printf("l:%ld case5 Q:%f e0:%f e1:%f T:%f dw:%f wliq:%f wice:%f th:%f theq:%f \n",l,Q,e0->co[l],e1->co[l],T->co[l],dw->co[l],th*A,thi*A,th,theq);
					}else{	//part of ice is melted
						dw->co[l]=Q/Lf;
						T->co[l]=e0->co[l];
						if(r==ri && c==ci) printf("l:%ld case6 Q:%f e0:%f e1:%f T:%f dw:%f wliq:%f wice:%f th:%f theq:%f \n",l,Q,e0->co[l],e1->co[l],T->co[l],dw->co[l],th*A,thi*A,th,theq);
					}
				//there is ice enough for th to reach theq
				}else{
					if(Lf*(theq-th)*A>Q){	//theq can't be reached
						dw->co[l]=Q/Lf;
						T->co[l]=e0->co[l];
						if(r==ri && c==ci) printf("l:%ld case7 Q:%f e0:%f e1:%f T:%f dw:%f wliq:%f wice:%f th:%f theq:%f \n",l,Q,e0->co[l],e1->co[l],T->co[l],dw->co[l],th*A,thi*A,th,theq);
					}else{	//theq is reached
						dw->co[l]=(theq-th)*A;
						Q-=(Lf*dw->co[l]);
						//move on the curve
						th1=th+dw->co[l]/A;
						thi1=thi-dw->co[l]/A;
						h=internal_energy_soil(th1, thi1, Tin, d, Cs, sat);
						if(r==ri && c==ci) printf("l:%ld case8a) theq:%f thi:%f th:%f h:%f Q:%f dw:%f T:%f\n",l,theq,thi1,th1,h,Q,dw->co[l],Tin);
						from_internal_soil_energy(r, c, l, h+Q, &th1, &thi1, &(T->co[l]), sl->pa->co[sy], psimin);
						dw->co[l]+=(th1-theq)*A;
						if(r==ri && c==ci) printf("l:%ld case8b) thi:%f th:%f T:%f dw:%f theq:%f\n",l,thi1,th1,T->co[l],dw->co[l],teta_psi(Psif2(T->co[l]), 0.0, sat, res, a, n, m, psimin, 1.0));
					}
				}

			}else{
				Teq=P*(g*(tk+Tfreezing))/(1.0E3*Lf);
				if(e0->co[l]+Q/ct <= Teq){
					dw->co[l]=0.0;
					T->co[l]=e1->co[l];
					if(r==ri && c==ci)printf("l:%ld case9 Q:%f e0:%f Teq:%f e1:%f T:%f dw:%f wliq:%f wice:%f th:%f theq:%f \n",l,Q,e0->co[l],Teq,e1->co[l],T->co[l],dw->co[l],th*A,thi*A,th,theq);
				}else{
					Q-=ct*(Teq-e0->co[l]);
					//move on the curve
					th1=th;
					thi1=thi;
					h=internal_energy_soil(th1, thi1, Teq, d, Cs, sat);
					if(r==ri && c==ci) printf("l:%ld case10a) theq:%f thi:%f th:%f h:%f Q:%f dw:%f T:%f Teq:%f\n",l,theq,thi1,th1,h,Q,dw->co[l],Tin,Teq);
					from_internal_soil_energy(r, c, l, h+Q, &th1, &thi1, &(T->co[l]), sl->pa->co[sy], psimin);
					dw->co[l]=(th1-th)*A;
					if(r==ri && c==ci) printf("l:%ld case10b) thi:%f th:%f T:%f dw:%f theq:%f\n",l,thi1,th1,T->co[l],dw->co[l],teta_psi(Psif2(T->co[l]), 0.0, sat, res, a, n, m, psimin, 1.0));
				}
			}
		}

		//stop_execution();

		/*if(e1->co[l]<=e0->co[l]){
			if(T->co[l]<e1->co[l]-1.E-5 || T->co[l]>e0->co[l]+1.E-5){
				printf("ERROR 1 in freezing sl r:%ld c:%ld l:%ld e0:%f e1:%f T:%f dw:%f th:%f thi:%f\n",r,c,l,e0->co[l],e1->co[l],T->co[l],dw->co[l],th,thi);
				f=fopen(files->co[ferr]+1,"a");
				fprintf(f,"ERROR 1 in freezing sl r:%ld c:%ld l:%ld e0:%f e1:%f T:%f dw:%f th:%f thi:%f\n",r,c,l,e0->co[l],e1->co[l],T->co[l],dw->co[l],th,thi);
				fclose(f);
				//stop_execution();
			}
		}else{
			if(T->co[l]<e0->co[l]-1.E-5 || T->co[l]>e1->co[l]+1.E-5){
				printf("ERROR 2 in freezing sl r:%ld c:%ld l:%ld e0:%f e1:%f T:%f dw:%f th:%f thi:%f\n",r,c,l,e0->co[l],e1->co[l],T->co[l],dw->co[l],th,thi);
				f=fopen(files->co[ferr]+1,"a");
				fprintf(f,"ERROR 2 in freezing sl r:%ld c:%ld l:%ld e0:%f e1:%f T:%f dw:%f th:%f thi:%f\n",r,c,l,e0->co[l],e1->co[l],T->co[l],dw->co[l],th,thi);
				fclose(f);
				//stop_execution();
			}
		}*/
		if(r==ri && c==ci) printf("\n");

	}else{
		T->co[l]=e1->co[l];
		dw->co[l]=0.0;
	}

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double k_thermal_snow_Sturm(double density){	//W m^-1 K^-1 (Sturm, 1997)

	double kt;
	density*=0.001;
	if(density<0.156){
		kt=0.023+0.234*density;
	}else{
		kt=0.138-1.01*density+3.233*pow(density,2.0);
	}
	return(kt);
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double k_thermal_snow_Yen(double density){	//W m^-1 K^-1 (Yen, 1981)

	double kt;
	kt=k_ice*pow((density/rho_w),1.88);
	return(kt);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double SB(double T){	//Stefan-Boltzmann law
	double R;
	R=5.67E-8*pow(T+tk,4.0);
	return(R);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double dSB_dT(double T){
	double dR_dT;
	dR_dT=4.0*5.67E-8*pow(T+tk,3.0);
	return(dR_dT);
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
void non_dimensionalize_snowage(double *snowage, double Ta){
	/* Author:  Stefano Endrizzi  Year:
	* Input:	snowage: snow age (days) of the snow on the pixel
	* 			Ta: air temperature
	* Output: snowage: new non-dimensionalized snow age
	* see Stefano Endrizzi PhD thesis, pg.69
	* comment: Matteo Dall'Amico, May 2009 */
	double r1, r2, r3;

	r1=exp(5000.0*(1.0/273.16-1.0/(Ta+273.16)));//effect of grain growth due to vapour diffusion
	r2=pow(r1,10);// additional effect near and at freezing point due to melting and refreezing
	if(r2>1.0) r2=1.0;
	r3=0.3;

	*snowage*=((r1+r2+r3)*1.0E-6);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void glacier_init_t0(long r, long c, double Ta, GLACIER *glac, SNOW *snow, PAR *par, double time){
	/* Author:  Stefano Endrizzi  Year:
	* function that initializes the glacier column verifying that the first glacier layer is deeper than
	* the minimum admitted depth. If this is true, than all the mass and energy parameters are then given
	* to the snow otherwise see glac_layer_combination()
	* Input:	r: row of the pixel in the raster
	* 			c: column of the pixel in the raster
	* 			Ta: air temperature
	* 		snow:	snow structure
	* 		glac: glacier structure
	* 		par: parameter structure
	* 		time: current time of the simulation
	* Output: snow and glac modified
	* comment: Matteo Dall'Amico, May 2009 */
	double D;
	double h;// internal energy of the layer
	long l;

	D=DEPTH(r, c, glac->lnum, glac->Dzl);

	//when glacier depth is too low, glacier becomes a snow layer and all the mass and internal energy is given to the first glacier layer
	if(D>0 && D<par->Dmin_glac->co[1]){
		h=(c_ice*glac->w_ice->co[1][r][c] + c_liq*glac->w_liq->co[1][r][c])*(glac->T->co[1][r][c]-Tfreezing) + Lf*glac->w_liq->co[1][r][c];
		for(l=2;l<=glac->lnum->co[r][c];l++){
			glac->Dzl->co[1][r][c]+=glac->Dzl->co[l][r][c];
			glac->w_liq->co[1][r][c]+=glac->w_liq->co[l][r][c];
			glac->w_ice->co[1][r][c]+=glac->w_ice->co[l][r][c];
			//means that glac->T is not a novalue
			if(glac->T->co[l][r][c]>-98.999) h+=(c_ice*glac->w_ice->co[l][r][c] + c_liq*glac->w_liq->co[l][r][c])*(glac->T->co[l][r][c]-Tfreezing) + Lf*glac->w_liq->co[l][r][c];

			glac->Dzl->co[l][r][c]=0.0;
			glac->w_liq->co[l][r][c]=0.0;
			glac->w_ice->co[l][r][c]=0.0;
			glac->T->co[l][r][c]=-99.0;
		}
		// the mass and internal energy of the glacier layer is then added with that of the snow
		h+=(c_ice*snow->w_ice->co[1][r][c] + c_liq*snow->w_liq->co[1][r][c])*(snow->T->co[1][r][c]-Tfreezing) + Lf*snow->w_liq->co[1][r][c];
		snow->Dzl->co[1][r][c]+=glac->Dzl->co[1][r][c];
		snow->w_ice->co[1][r][c]+=glac->w_ice->co[1][r][c];
		snow->w_liq->co[1][r][c]+=glac->w_liq->co[1][r][c];
		if(h<0){
			snow->T->co[1][r][c]=Tfreezing + h/(c_ice*snow->w_ice->co[1][r][c]+c_liq*snow->w_liq->co[1][r][c]);
		}else{
			snow->T->co[1][r][c]=0.0;
		}

		glac->Dzl->co[1][r][c]=0.0;
		glac->w_liq->co[1][r][c]=0.0;
		glac->w_ice->co[1][r][c]=0.0;
		glac->T->co[1][r][c]=-99.0;
		glac->lnum->co[r][c]=0;

		snow_layer_combination(r, c, snow, Ta, par->snowlayer_inf, par->Dmin, par->Dmax, time);
	}

	//setting
	if(par->glaclayer_max>1)glac_layer_combination(r, c, glac, Ta, par->glaclayer_max, par->Dmin_glac, par->Dmax_glac, time);
}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void SolveEB(long n1, long n2, double dt, double Tg, double h0, double h1, double dhdT, double *SW, long nlim, double Gb, double *k, double *C, double *D, double *T,
			 double *dw, double *wi, double *wl, short snowtype, double Ta, SOIL *sl, long r, long c, double t, PAR *par){

	//n1: snow (and glacier)
	//n2: sl
	//n=n1+n2: total

	short occuring, sy=sl->type->co[r][c];
	long l, n=n1+n2;
	double Dcorr;
	SHORTVECTOR *mf;
	DOUBLEVECTOR *ad, *ads, *adi, *b, *ad1, *ads1, *adi1, *b1, *e;
	DOUBLEVECTOR *e0g,*e1g, *dwg, *Tsoil, *Cg;
	long ri=12,ci=13;
	//FILE *f;

	mf=new_shortvector(n);
	initialize_shortvector(mf,0);

	ad=new_doublevector(n);
	ads=new_doublevector(n-1);
	adi=new_doublevector(n-1);
	b=new_doublevector(n);
	ad1=new_doublevector(n);
	ads1=new_doublevector(n-1);
	adi1=new_doublevector(n-1);
	b1=new_doublevector(n);

	e=new_doublevector(n);

	e0g=new_doublevector(n2);
	e1g=new_doublevector(n2);
	dwg=new_doublevector(n2);
	Tsoil=new_doublevector(n2);
	Cg=new_doublevector(n2);

	Dcorr=D[1];
	//Dcorr=0.5*(0.5*D->co[1]+CA*(D->co[1]+0.5*D->co[2]));

	ad->co[1]=-(1.0-KNe)*dhdT + 2.0*(1-KNe)*k[1]/(D[1]+D[2]) + C[1]*Dcorr/dt;
	ads->co[1]=2.0*(KNe-1)*k[1]/(D[1]+D[2]);
	b->co[1]=C[1]*Dcorr*T[1]/dt + 2.0*KNe*k[1]*(T[2]-T[1])/(D[1]+D[2]) + SW[1] + KNe*h0 + (1.0-KNe)*(h1-dhdT*Tg);
	if(r==ri && c==ci) printf("B l:%ld b:%f 1:%f 2:%f 3:%f 4:%f\n", 1, b->co[1], C[1]*Dcorr*T[1]/dt, 2.0*KNe*k[1]*(T[2]-T[1])/(D[1]+D[2]), SW[1]+KNe*h0, (1.0-KNe)*(h1-dhdT*Tg));

	for(l=2;l<=n-1;l++){
		ad->co[l]=C[l]*D[l]/dt + 2.0*(1-KNe)*k[l-1]/(D[l-1]+D[l]) + 2.0*(1-KNe)*k[l]/(D[l]+D[l+1]);
		ads->co[l]=2.0*(KNe-1)*k[l]/(D[l]+D[l+1]);
		adi->co[l-1]=ads->co[l-1];	//la matrice e' simmetrica
		b->co[l]=C[l]*D[l]*T[l]/dt - 2.0*KNe*k[l-1]*(T[l]-T[l-1])/(D[l-1]+D[l]) + 2.0*KNe*k[l]*(T[l+1]-T[l])/(D[l]+D[l+1]);
		if(l<=nlim) b->co[l]+=SW[l];
		if(r==ri && c==ci) printf("B l:%ld b:%f T:%f 1:%f 2:%f 3:%f ",l,b->co[l],T[l],C[l]*D[l]*T[l]/dt,-2.0*KNe*k[l-1]*(T[l]-T[l-1])/(D[l-1]+D[l]),2.0*KNe*k[l]*(T[l+1]-T[l])/(D[l]+D[l+1]));
		if(l<=nlim && r==ri && c==ci) printf("4:%f",SW[l]);
		if(r==ri && c==ci) printf("\n");
	}

	ad->co[n]=C[n]*D[n]/dt + 2.0*(1-KNe)*k[n-1]/(D[n-1]+D[n]);
	adi->co[n-1]=ads->co[n-1];		//la matrice e' simmetrica
	b->co[n]=C[n]*D[n]*T[n]/dt - 2.0*KNe*k[n-1]*(T[n]-T[n-1])/(D[n-1]+D[n]) + Gb;
	if(r==ri && c==ci) printf("B l:%ld b:%f 1:%f 2:%f 3:%f\n",n,b->co[n],C[n]*D[n]*T[n]/dt,-2.0*KNe*k[n-1]*(T[n]-T[n-1])/(D[n-1]+D[n]),Gb);

	tridiag(1,r,c,n,adi,ad,ads,b,e);
	check_errors(r,c,n,adi,ad,ads,b,e,T,mf);

	/*for(l=1;l<=n;l++){
		if(r==ri && c==ci){
			printf("1. l:%ld h:%f abs:%f sum:%f ad:%f b:%f T:%f dw:%f wliq:%f wice:%f e:%f mf:%ld",l,h,abs->co[1],h+abs->co[1],ad->co[l],b->co[l],T[l],dw[l],wl[l],wi[l],e->co[l],mf->co[l]);
			if(l>=n1+1 && l<=n1+n2) printf(" wice2:%f ",sl->thice->co[l-n1][r][c]*sl->pa->co[sy][jdz][l-n1]);
			printf("\n");
		}
	}*/

	//nel caso neve semplificato (type=1), se non si verifica scioglimento, non si fa il bilancio di energia!!
	if(snowtype==1 && e->co[1]<Ta-20){
		mf->co[1]=99;		//Codice che memorizza questo passaggio
		ad->co[1]=1.0;
		ads->co[1]=0.0;
		if(Ta-20<Tfreezing){	//Si pone la temperatura della neve pari a quella dell'aria
			b->co[1]=Ta-20;
		}else{
			b->co[1]=Tfreezing;
		}
		tridiag(1,r,c,n,adi,ad,ads,b,e);
		check_errors(r,c,n,adi,ad,ads,b,e,T,mf);
		/*for(l=1;l<=n;l++){
			if(r==ri && c==ci) printf("1b. l:%ld h:%f sum:%f abs:%f ad:%f b:%f T:%f dw:%f wliq:%f wice:%f e:%f mf:%ld \n",l,h,abs->co[1],h+abs->co[1],ad->co[l],b->co[l],T[l],dw[l],wl[l],wi[l],e->co[l],mf->co[l]);
		}*/
	}

	//melting e freezing
	//riinizializzazione
	occuring=0;
	for(l=1;l<=n;l++){
		b1->co[l]=b->co[l];
		ad1->co[l]=ad->co[l];
		if(l!=n){
			ads1->co[l]=ads->co[l];
			adi1->co[l]=adi->co[l];
		}
	}

	//snow and glacier
	for(l=1;l<=n1;l++){
		if(mf->co[l]!=99){
			if(e->co[l]>Tfreezing && wi[l]>0.0){
				mf->co[l]=1;//melting
				occuring=1;
			}else if(e->co[l]<Tfreezing && wl[l]>0.0){
				mf->co[l]=2;//freezing
				occuring=1;
			}
		}
	}

	//sl freezing
	for(l=1;l<=n2;l++){
		e0g->co[l]=T[l+n1];
		e1g->co[l]=e->co[l+n1];
		Cg->co[l]=C[l+n1];
		if(sl->pa->co[sy][jsf][l]==1){
			occuring=1;
			soil_freezing(l, r, c, e0g, e1g, dwg, Tsoil, Cg, sl, par->psimin);
			if(r==ri && c==ci) printf("frez1. l:%ld e0g:%f e1g:%f Tsoil:%f dwg:%f\n",l,e0g->co[l],e1g->co[l],Tsoil->co[l],dwg->co[l]);
			if(Tsoil->co[l]!=Tsoil->co[l]) printf("NOvalue Tg1: r:%ld c:%ld l:%ld wi:%f wl:%f dw:%f mf:%ld T:%f\n",r,c,l,wi[l],wl[l],dw[l],mf->co[l],T[l]);
		}
	}

	if(occuring==0){	//non devo fare altri calcoli, scrivo i risultati

		for(l=1;l<=n;l++){
			T[l]=e->co[l];
			dw[l]=0.0;
		}

	}else{  //cioe' se accade melting o freezing

		for(l=1;l<=n;l++){

			/*se AVVIENE MELTING O FREEZING, l'equazione cambia*/

			//NEVE O GHIACCIO: l'incognita diventa l'acqua che si scioglie
			if(mf->co[l]==1 || mf->co[l]==2){
				ad1->co[l]=Lf/dt;
				b1->co[l]-=ad->co[l]*Tfreezing;
			}

			/*vengono corretti i termini di un layer se il layer soprastante e' soggetto a melt/freez*/
			if(l!=1){
				if(mf->co[l-1]==1 || mf->co[l-1]==2){
					if(mf->co[l]!=99){
						adi1->co[l-1]=0.0;
						b1->co[l]-=adi->co[l-1]*Tfreezing;
					}
				}
			}

			/*vengono corretti i termini di un layer se il layer sottostante e' soggetto a melt/freez*/
			if(l!=n){
				if(mf->co[l+1]==1 || mf->co[l+1]==2){
					if(mf->co[l]!=99){
						ads1->co[l]=0.0;
						b1->co[l]-=ads->co[l]*Tfreezing;
					}
				}
			}

			/*sl freezing*/
			if(l>n1){
				if(sl->pa->co[sy][jsf][l-n1]==1){
					if(l>1){
						b1->co[l-1]-=ads->co[l-1]*Tsoil->co[l-n1];
						ads1->co[l-1]=0.0;
						if(r==ri && c==ci) printf("#. l:%ld B:%f B1:%f Tgbelow:%f m:%ld Ads:%f\n",l-1,b->co[l-1],b1->co[l-1],Tsoil->co[l-n1],l-n1,ads->co[l-1]);
					}
					if(l<n){
						b1->co[l+1]-=adi->co[l]*Tsoil->co[l-n1];
						adi1->co[l]=0.0;
					}
				}
			}
		}

		tridiag(1,r,c,n,adi1,ad1,ads1,b1,e);
		check_errors(r,c,n,adi1,ad1,ads1,b1,e,T,mf);

		/*for(l=1;l<=n;l++){
			if(r==ri && c==ci) printf("2. l:%ld h:%f abs:%f sum:%f ad:%f b:%f T:%f dw:%f wliq:%f wice:%f e:%f mf:%ld \n",l,h,abs->co[1],h+abs->co[1],ad1->co[l],b1->co[l],T[l],dw[l],wl[l],wi[l],e->co[l],mf->co[l]);
		}*/

		//riinizializzazione
		occuring=0;
		for(l=1;l<=n;l++){
			b1->co[l]=b->co[l];
			ad1->co[l]=ad->co[l];
			if(l!=n){
				ads1->co[l]=ads->co[l];
				adi1->co[l]=adi->co[l];
			}
		}

		//snow and glacier
		for(l=1;l<=n1;l++){

			//MELTING
			if(mf->co[l]==1){
				//a. si verifica melting (mf=1), ma e e' negativo (FALSO MELTING)
				/*if(e->co[l]<0.0){
					mf->co[l]=100;
					occuring=1;
				//b. c'e' meno ghiaccio di quanto potrebbe sciogliersi
				}else*/ if(e->co[l]>wi[l]){
					mf->co[l]=10;
					occuring=1;
				}

			//FREEZING (snow and glacier)
			}else if(mf->co[l]==2){
				//a. si verifica freezing (mf=2), ma e e' positivo (FALSO FREEZING)
				/*if(e->co[l]>0.0){
					mf->co[l]=200;
					occuring=1;
				//c'e' meno acqua liquida di quanta ne potrebbe ghiacciare
				}else*/if(-e->co[l]>wl[l]){
					mf->co[l]=20;
					occuring=1;
				}
			}
		}

		//sl freezing
		for(l=1;l<=n2;l++){
			//e0g->co[l]=T[l+n1];
			e1g->co[l]=e->co[l+n1];
			//Cg->co[l]=C[l+n1];
			if(sl->pa->co[sy][jsf][l]==1){
				occuring=1;
				soil_freezing(l, r, c, e0g, e1g, dwg, Tsoil, Cg, sl, par->psimin);
				if(r==ri && c==ci) printf("frez2. l:%ld e0g:%f e1g:%f Tsoil:%f dwg:%f\n",l,e0g->co[l],e1g->co[l],Tsoil->co[l],dwg->co[l]);
			}
		}

		if(occuring==0){	//scrivo i risultati

			for(l=1;l<=n;l++){
				if(mf->co[l]==0 || mf->co[l]==99){
					T[l]=e->co[l];
					dw[l]=0.0;
				}else if(mf->co[l]==1 || mf->co[l]==2){	//melting or freezing
					T[l]=Tfreezing;
					dw[l]=e->co[l];
				}
			}

		}else{		//allora le equazioni cambiano

			for(l=1;l<=n;l++){

				if(mf->co[l]==1 || mf->co[l]==2){
					ad1->co[l]=Lf/dt;
					b1->co[l]-=ad->co[l]*Tfreezing;

				}else if(mf->co[l]==10){
					b1->co[l]-=Lf*wi[l]/dt;

				}else if(mf->co[l]==20){
					b1->co[l]+=Lf*wl[l]/dt;
				}

				/*vengono corretti i termini di un layer se il layer soprastante e' soggetto a melt/freez*/
				if(l!=1){
					if(mf->co[l-1]==1 || mf->co[l-1]==2){
						if(mf->co[l]!=99){
							adi1->co[l-1]=0.0;
							b1->co[l]-=adi->co[l-1]*Tfreezing;
						}
					}
				}

				/*vengono corretti i termini di un layer se il layer sottostante e' soggetto a melt/freez*/
				if(l!=n){
					if(mf->co[l+1]==1 || mf->co[l+1]==2){
						if(mf->co[l]!=99){
							ads1->co[l]=0.0;
							b1->co[l]-=ads->co[l]*Tfreezing;
						}
					}
				}

				/*sl freezing*/
				if(l>n1){
					if(sl->pa->co[sy][jsf][l-n1]==1){
						if(l>1){
							b1->co[l-1]-=ads->co[l-1]*Tsoil->co[l-n1];
							ads1->co[l-1]=0.0;
							if(r==ri && c==ci) printf("#. l:%ld B:%f B1:%f Tgbelow:%f m:%ld Ads:%f\n",l-1,b->co[l-1],b1->co[l-1],Tsoil->co[l-n1],l-n1,ads->co[l-1]);
						}
						if(l<n){
							b1->co[l+1]-=adi->co[l]*Tsoil->co[l-n1];
							adi1->co[l]=0.0;
						}
					}
				}
			}

			tridiag(1,r,c,n,adi1,ad1,ads1,b1,e);
			check_errors(r,c,n,adi1,ad1,ads1,b1,e,T,mf);

			for(l=1;l<=n;l++){
				if(r==ri && c==ci) printf("3. l:%ld ad:%f b:%f T:%f dw:%f wliq:%f wice:%f e:%f mf:%ld \n",l,ad1->co[l],b1->co[l],T[l],dw[l],wl[l],wi[l],e->co[l],mf->co[l]);
			}

			occuring=0;
			for(l=1;l<=n;l++){
				if(e->co[l]>Tfreezing && mf->co[l]==10){
					occuring=1;
					mf->co[l]=98;
					ad1->co[l]=1.0;
					b->co[l]=Tfreezing;
					if(l!=1){
						adi1->co[l-1]=0.0;
						if(mf->co[l-1]!=99){
							b1->co[l-1]-=ads->co[l-1]*Tfreezing;
							ads1->co[l-1]=0.0;
						}
					}
					if(l!=n){
						ads1->co[l]=0.0;
						if(mf->co[l+1]!=99){
							b1->co[l+1]-=adi->co[l]*Tfreezing;
							adi1->co[l]=0.0;
						}
					}
				}
			}

			if(occuring==1){
				tridiag(1,r,c,n,adi1,ad1,ads1,b,e);
				check_errors(r,c,n,adi1,ad1,ads1,b,e,T,mf);
				for(l=1;l<=n;l++){
					if(r==ri && c==ci)printf("3b. l:%ld ad:%f b:%f T:%f dw:%f wliq:%f wice:%f e:%f mf:%ld\n ",l,ad1->co[l],b1->co[l],T[l],dw[l],wl[l],wi[l],e->co[l],mf->co[l]);
				}
			}

			//scrivo i risultati

			//snow and glacier
			for(l=1;l<=n1;l++){
				if(mf->co[l]==0 || mf->co[l]==99){
					T[l]=e->co[l];
					dw[l]=0.0;
				}else if(mf->co[l]==1){
					T[l]=Tfreezing;
					dw[l]=e->co[l];
				}else if(mf->co[l]==2){
					T[l]=Tfreezing;
					dw[l]=e->co[l];
				}else if(mf->co[l]==10 || mf->co[l]==98){
					T[l]=e->co[l];
					dw[l]=wi[l];
				}else if(mf->co[l]==20){
					T[l]=e->co[l];
					dw[l]=-wl[l];
				}else if(mf->co[l]==100){	//FALSO MELTING
					T[l]=e->co[l];
					dw[l]=0.0;
				}else if(mf->co[l]==200){	//FALSO FREEZING
					T[l]=e->co[l];
					dw[l]=0.0;
				}
				if(wl[l]!=wl[l]) printf("NOvalue wl[i]snow: r:%ld c:%ld l:%ld wi:%f wl:%f dw:%f mf:%ld T:%f\n",r,c,l,wi[l],wl[l],dw[l],mf->co[l],T[l]);
			}

			//sl freezing
			for(l=1;l<=n2;l++){
				//e0g->co[l]=T[l+n1];
				e1g->co[l]=e->co[l+n1];
				//Cg->co[l]=C[l+n1];
				if(sl->pa->co[sy][jsf][l]==1){
					soil_freezing(l, r, c, e0g, e1g, dwg, Tsoil, Cg, sl, par->psimin);
					if(r==ri && c==ci) printf("frez3. l:%ld e0g:%f e1g:%f Tsoil:%f dwg:%f\n",l,e0g->co[l],e1g->co[l],Tsoil->co[l],dwg->co[l]);
				}
			}

			for(l=n1+1;l<=n1+n2;l++){
				if(sl->pa->co[sy][jsf][l-n1]==1){
					T[l]=Tsoil->co[l-n1];
					dw[l]=dwg->co[l-n1];
					if(wl[l]!=wl[l]) printf("NOvalue wl[i]soil: r:%ld c:%ld l:%ld wi:%f wl:%f dw:%f mf:%ld T:%f\n",r,c,l,wi[l],wl[l],dw[l],mf->co[l],T[l]);
					if(r==ri && c==ci) printf("l:%ld m:%ld T:%f dw:%f\n",l,l-n1,T[l],dw[l]);
				}else{
					T[l]=e->co[l];
					dw[l]=0.0;
				}
			}



		}
	}

	for(l=1;l<=n;l++){

		if(r==ri && c==ci) printf("END l:%ld ad:%f b:%f T:%f dw:%f wliq:%f wice:%f e:%f mf:%ld ",l,ad->co[l],b->co[l],T[l],dw[l],wl[l],wi[l],e->co[l],mf->co[l]);

		//Continuity Check
		//1
		if(wl[l]+dw[l]<0){
			//f=fopen(files->co[ferr]+1,"a");
			//fprintf(f,"\nLIQ CONTENT IN SNOW TOO LOW OR NEGATIVE, cell r=%4ld c=%4ld l=%4ld time=%10.3f, CORRECTED, mf=%4ld tetaliq=%f dth:%f\n",r,c,l,t,mf->co[l],wl[l]/(rho_w*D[l]),dw[l]/(rho_w*D[l]));
			//fclose(f);
			dw[l]=-wl[l];
		}
		//2
		if(wi[l]-dw[l]<0){
			//f=fopen(files->co[ferr]+1,"a");
			//fprintf(f,"\nICE CONTENT IN SNOW NEGATIVE, cell r=%4ld c=%4ld l=%4ld time=%10.3f, CORRECTED, mf=%4ld tetaice=%f dth:%f\n",r,c,l,t,mf->co[l],wi[l]/(rho_w*D[l]),dw[l]/(rho_w*D[l]));
			//fclose(f);
			dw[l]=wi[l];
		}

		//Assign
		wl[l]+=dw[l];
		wi[l]-=dw[l];

		if(r==ri && c==ci) printf(" dw:%f\n",dw[l]);

		if(wi[l]!=wi[l]){
			printf("NOvalue in egy balance4: r:%ld c:%ld l:%ld wi:%f wl:%f dw:%f mf:%ld T:%f\n",r,c,l,wi[l],wl[l],dw[l],mf->co[l],T[l]);
			stop_execution();
		}

		if(wl[l]!=wl[l]){
			printf("NOvalue in egy balance5: r:%ld c:%ld l:%ld wi:%f wl:%f dw:%f mf:%ld T:%f\n",r,c,l,wi[l],wl[l],dw[l],mf->co[l],T[l]);
			stop_execution();
		}

	}


	free_shortvector(mf);
	free_doublevector(ad);
	free_doublevector(ads);
	free_doublevector(adi);
	free_doublevector(b);
	free_doublevector(ad1);
	free_doublevector(ads1);
	free_doublevector(adi1);
	free_doublevector(b1);
	free_doublevector(e);
	free_doublevector(e0g);
	free_doublevector(e1g);
	free_doublevector(dwg);
	free_doublevector(Tsoil);
	free_doublevector(Cg);


}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void soil_properties(long r, long c, long beg, long end, double *th_l, double ***th_i, double *D, double *wl, double *wi, double *k, double *C, double *T, SOIL *sl){

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

		//sl heat capacity [J m^-3 K^-1]
		C[l]=sl->pa->co[sy][jct][m]*(1-sl->pa->co[sy][jsat][m]) + c_ice*th_i[m][r][c]*rho_w + c_liq*th_l[m]*rho_w;

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
void snow_properties(long r, long c, long beg, long end, double *D, double *wl, double *wi, double *k, double *C, double *T, double ***tdz, double ***twl, double ***twi,
					 double ***tT, double (* kfunct)(double r)){

	long l,m;
	double rho;

	for(l=beg;l<=end;l++){

		m=end-l+1;

		D[l]=1.0E-3*tdz[m][r][c];
		wl[l]=twl[m][r][c];
		wi[l]=twi[m][r][c];
		T[l]=tT[m][r][c];
		C[l]=(c_ice*wi[l]+c_liq*wl[l])/D[l];

		rho=(wl[l]+wi[l])/D[l];
		k[l]=(*kfunct)(rho);
	}
}

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
			double LE, double surfEB, double H, double G, double Ts, double SWin, double SWout, double SWbeam, double eps, double LWin, double cosinc){

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
		if(par->distr_stat==1){
			if(egy->ET_max->co[r][c]<LE) egy->ET_max->co[r][c]=LE;
			if(egy->ET_min->co[r][c]>LE) egy->ET_min->co[r][c]=LE;
		}
	}

	if(par->output_G>0){
		egy->G_mean->co[r][c]+=surfEB/((par->output_G*3600.0)/(par->Dt)); //[W/m^2]
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
		if(par->micromet1==1){
			met->Vspdmean->co[r][c]+=met->Vgrid->co[r][c]/((par->output_meteo*3600.0)/(par->Dt));
			met->Vdirmean->co[r][c]+=met->Vdir->co[r][c]/((par->output_meteo*3600.0)/(par->Dt));
			met->RHmean->co[r][c]+=met->RHgrid->co[r][c]/((par->output_meteo*3600.0)/(par->Dt));
		}
	}

	if(par->output_Rn>0){
		egy->Rn_mean->co[r][c]+=(SWin-SWout + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)))/((par->output_Rn*3600.0)/(par->Dt)); //[W/m^2]
		if(fabs(SWin-SWout + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)))>1500){
			f=fopen(files->co[ferr]+1,"a");
			fprintf(f,"\ntime=%10.1f r=%4ld c=%4ld Rn=%10.3f Rsw=%10.3f Rlwdiff=%10.3f albedo=%10.8f eps=%10.8fTa=%10.5f Ts=%10.5f Rsw_meas=%f sin(alpha)=%f cos(inc)=%f\n",
				times->time,r,c,SWin-SWout + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)),SWin,LWin,A,
				eps,met->Tgrid->co[r][c],Ts,met->var[0][iSW],sin(egy->hsun),cosinc);
			fprintf(f,"\nH:%f LE:%f\n",H,LE);
			fclose(f);
		}
		egy->LW_in->co[r][c]+=( top->sky->co[r][c]*5.67E-8*eps*pow((Ts+tk),4.0) )/((par->output_Rn*3600.0)/(par->Dt));
		egy->LW_out->co[r][c]+=( top->sky->co[r][c]*eps*LWin )/((par->output_Rn*3600.0)/(par->Dt));
		egy->SW->co[r][c]+=( SWin-SWout )/((par->output_Rn*3600.0)/(par->Dt));
		if(par->distr_stat==1){
			if(egy->LW_max->co[r][c]<top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)))
				egy->LW_max->co[r][c]=top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0));
			if(egy->LW_min->co[r][c]>top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)))
				egy->LW_min->co[r][c]=top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0));
			if(egy->Rn_max->co[r][c]<(SWin-SWout + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0))))
				egy->Rn_max->co[r][c]=(SWin-SWout  + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)));
			if(egy->Rn_min->co[r][c]>(SWin-SWout  + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0))))
				egy->Rn_min->co[r][c]=(SWin-SWout  + top->sky->co[r][c]*(eps*LWin - 5.67E-8*eps*pow((Ts+tk),4.0)));
			if(egy->SW_max->co[r][c]<(SWin-SWout)) egy->SW_max->co[r][c]=SWin-SWout;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_pixel(long r, long c, double prec_snow, double prec_rain_on_soil, double prec_rain_on_snow, double Sr_soil, double Er_soil, double Mr_snow, double Er_snow, double Sr_snow,
			double Mr_glac, double Er_glac, double Sr_glac, double Evt, double LE, double H, double surfEB, double G, double Tg, double A, double eps, double LAI,
			double LWin, double SWin, double LWout, double SWout, double epsa, double epsa_min, double epsa_max, double SWbeam, double SWdiff, double DTcorr, double Tdew,
			long n, DOUBLEVECTOR *turbulence, PAR *par, WATER *wat, ENERGY *egy, TOPO *top, METEO *met, SNOW *snow, GLACIER *glac, LAND *land, double Vpoint, double RHpoint,
			double prec_snow_atm, double prec_rain_atm, double maxstorage_c, double evap_c, double drip_c, double z0, double d0, double SWv, double LWv, double Hv, double LEv,
			double Tv, double Ts, double Hg0, double Hg1, double Eg0, double Eg1, double fc, double rh, double rv, double rb, double rc, double rh_ic, double rv_ic,
			double Hv1, double LEv1, double Qv, double Qg, double Qa, double Qs){


	long i, j, l;
	double ea, es, de_dT, Q;

	if(par->state_pixel==1){
		for(i=1;i<=par->chkpt->nrh;i++){
			if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){

				wat->out1->co[3][i]+=prec_snow;
				wat->out1->co[4][i]+=(prec_rain_on_soil+prec_rain_on_snow);

				wat->out1->co[16][i]+=prec_snow_atm;
				wat->out1->co[17][i]+=prec_rain_atm;
				wat->out1->co[18][i]+=maxstorage_c/(double)n;
				wat->out1->co[19][i]+=evap_c;
				wat->out1->co[20][i]+=drip_c;

				wat->out1->co[21][i]+=prec_snow;
				wat->out1->co[22][i]+=prec_rain_on_soil;
				wat->out1->co[23][i]+=prec_snow_atm;
				wat->out1->co[24][i]+=prec_rain_atm;
				wat->out1->co[25][i]+=evap_c;
				wat->out1->co[26][i]+=drip_c;
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
				egy->out1->co[12][i]+=eps*LWin/(double)n;
				egy->out1->co[22][i]+=eps*LWin*par->Dt*1.0E-6;
				//Rlw_out and Rlw_out_cumulated
				egy->out1->co[13][i]-=LWout/(double)n;
				egy->out1->co[23][i]-=LWout*par->Dt*1.0E-6;
				//Rsw, Rsw_incoming_cumulated, albedo, Rsw_outcoming_cumulated
				egy->out1->co[14][i]+=SWin/(double)n;
				egy->out1->co[24][i]+=SWin*par->Dt*1.0E-6;
				if(SWout>0){
					egy->out1->co[15][i]+=(SWout/SWin)/(double)n;
				}else{
					egy->out1->co[15][i]+=A/(double)n;
				}
				egy->out1->co[28][i]-=SWout*par->Dt*1.0E-6;
				egy->out1->co[45][i]-=SWout/(double)n;
				//atmosphere emissivity
				egy->out1->co[16][i]+=epsa/(double)n;
				//wind speed
				egy->out1->co[17][i]+=Vpoint/(double)n;
				//relative humidity
				egy->out1->co[18][i]+=RHpoint/(double)n;
				//atmospheric pressure
				egy->out1->co[19][i]+=met->Pgrid->co[r][c]/(double)n;
				//air temperature
				egy->out1->co[20][i]+=met->Tgrid->co[r][c]/(double)n;
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
				sat_vap_pressure(&ea, &de_dT, met->Tgrid->co[r][c], met->Pgrid->co[r][c]);
				Q=RHpoint*0.622*ea/(met->Pgrid->co[r][c]-0.378*ea);
				ea=Q*met->Pgrid->co[r][c]/(0.622+Q*0.378);
				egy->out1->co[47][i]+=ea/(double)n;
				egy->out1->co[48][i]+=Q/(double)n;
				sat_vap_pressure(&es, &de_dT, Ts, met->Pgrid->co[r][c]);
				egy->out1->co[49][i]+=es/(double)n;
				egy->out1->co[50][i]+=0.622*es/(met->Pgrid->co[r][c]-0.378*es)/(double)n;

				if(epsa_min>0){
					egy->out1->co[51][i]+=(eps*epsa_min*5.67E-8*pow(met->Tgrid->co[r][c]+tk,4.0))/(double)n;
				}else{
					egy->out1->co[51][i]=UV->V->co[2];
				}
				if(epsa_max>0){
					egy->out1->co[52][i]+=(eps*epsa_max*5.67E-8*pow(met->Tgrid->co[r][c]+tk,4.0))/(double)n;
				}else{
					egy->out1->co[52][i]=UV->V->co[2];
				}
				egy->out1->co[53][i]+=(SWbeam/(double)n);
				egy->out1->co[54][i]+=(SWdiff/(double)n);

				if(par->micromet1==1){
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
	}else{
		F=est;
	}

	//printf("%ld %ld %f %f\n",col[i-1][icol],icol,est,F);

	return(F);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void evaporation(long r, long c, long n, double E, double *wi, double *wl, double *D, double dt, double *Ss, double *Es, double *Eg){

	double thres=0.20;	//threshold on liquid content in the snow after which evaporation (not sublimation) occurs
	double th_liq;		//liquid content

	*Es=0.0;
	*Ss=0.0;
	*Eg=0.0;

	if(n>0){ //snow or glacier
		th_liq=wl[1]/(rho_w*D[1]);
		if(th_liq>thres){
			wl[1]-=E*dt;
			if(wl[1]<0){
				wi[1]+=wl[1];
				E+=wl[1]/dt;
				wl[1]=0.0;
				if(wi[1]<0){
					E+=wi[1]/dt;
					wi[1]=0.0;
				}
			}
			*Es=E;	//[mm/s] water equivalent
		}else{
			wi[1]-=E*dt;
			if(wi[1]<0){
				E+=wi[1]/dt;
				wi[1]=0.0;
			}
			*Ss=E;
		}
	}else{	//sl
		*Eg=E;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void update_soil(long ntot, long r, long c, SOIL *sl, double *T, double *dw, double *wi, double *th, double *ft, double Evt, double E, PAR *par){

	long l, m;
	double source,theta,A,dt=par->Dt;
	short sy=sl->type->co[r][c];

	for(l=1;l<=Nl;l++){
		m=l+ntot-Nl;
		A=sl->pa->co[sy][jdz][l]*0.001*rho_w;

		sl->T->co[l][r][c]=T[m];
		sl->thice->co[l][r][c]=wi[m]/A;		//[-]

	}

	for(l=1;l<=Nl;l++){
		m=l+ntot-Nl;
		A=sl->pa->co[sy][jdz][l]*0.001*rho_w;

		if(l==1){
			source=dw[m] - fmax( dt*(E+Evt*ft[1]), 0.0 );	//[kg/m2]=[mm water]
			//source=dw[m];
		}else{
			source=dw[m] - fmax( dt*Evt*ft[l], 0.0 );
			//source=dw[m];
		}

		theta=th[l]+source/A;
		if(th[l]>=sl->pa->co[sy][jsat][l]){
			if(theta>th[l]) theta=th[l];
		}else{
			if(theta>sl->pa->co[sy][jsat][l]) theta=sl->pa->co[sy][jsat][l];
		}

		A=sl->pa->co[sy][jdz][l]*0.001*rho_w;

		//printf("%ld psiante:%f theta:%f theta:%f\n",l,sl->P->co[l][r][c],teta_psi(sl->P->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		//		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil),th[l]);

		sl->P->co[l][r][c]=psi_teta(theta, sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
				sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil);

		//printf("psipost:%f theta:%f theta:%f source:%f dw:%f e:%f ice:%f\n",sl->P->co[l][r][c],teta_psi(sl->P->co[l][r][c], sl->thice->co[l][r][c], sl->pa->co[sy][jsat][l], sl->pa->co[sy][jres][l], sl->pa->co[sy][ja][l],
		//		sl->pa->co[sy][jns][l], 1-1/sl->pa->co[sy][jns][l], fmin(sl->pa->co[sy][jpsimin][l], Psif(sl->T->co[l][r][c], par->psimin)), par->Esoil),theta,
		//		source/A,dw[m]/A,(source-dw[m])/A,sl->thice->co[l][r][c]);

		//stop_execution();

	}
	//printf("\n\n");
}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void liqWBsnow(long r, long c, SNOW *snow, double *Mr, double *PonS, PAR *par, double slope, double P, double *wi, double *wl, double *T){

	long l, m;
	double Se, theta_liq, thice, Dw=P, wice_old;
	DOUBLEVECTOR *q;

	//rain on snow
	*PonS=0.0;	//initialization

	//a. ordinary case
	if(snow->type->co[r][c]==2){

		q=new_doublevector(snow->lnum->co[r][c]);

		for(l=snow->lnum->co[r][c];l>=1;l--){

			m=snow->lnum->co[r][c] - l + 1;	//increasing downwards (while l increases upwards)

			//ADD LIQUID WATER FROM ABOVE AND ACCOUNT FOR POSSIBLE REFREEZING
			if(l<snow->lnum->co[r][c]) Dw=q->co[l+1]*par->Dt;

			wl[m]+=Dw;

			//find new internal energy, taking into account the latent heat flux of the incoming water flux
			from_internal_energy(r, c, internal_energy(wi[m], wl[m], T[m]), &(wi[m]),&(wl[m]),&(T[m]));

			//assign snow state variables
			wice_old=snow->w_ice->co[l][r][c];
			snow->T->co[l][r][c]=T[m];
			snow->w_ice->co[l][r][c]=wi[m];
			snow->w_liq->co[l][r][c]=wl[m];

			//ACCOUNT FOR SNOW COMPACTION
			//a)destructive metamorphism and overburden
			snow_compactation(r, c, l, snow, slope, par, &(snow->CR1->co[l]), &(snow->CR2->co[l]));

			//b)melting: snow depth decreases maintaining the same density
			if(wi[m]/wice_old<1){
				snow->Dzl->co[l][r][c]*=(wi[m]/wice_old);
				snow->CR3->co[l]=(wi[m]-wice_old)/(wice_old*par->Dt);
			}else{
				snow->CR3->co[l]=0.0;
			}

			//check that snow porosity is not too large
			thice=snow->w_ice->co[l][r][c]/(1.0E-3*snow->Dzl->co[l][r][c]*rho_i);		//[-]
			if(thice>par->snow_maxpor){
				snow->Dzl->co[l][r][c]=1000.0*snow->w_ice->co[l][r][c]/(rho_i*par->snow_maxpor);
				thice=par->snow_maxpor;
			}

			//CALCULATE LIQUID WATER GOING BELOW
			if(snow->w_ice->co[l][r][c]<0.001){
				q->co[l]=snow->w_liq->co[l][r][c]/par->Dt;
				snow->w_liq->co[l][r][c]=0.0;
			}else{
				theta_liq=snow->w_liq->co[l][r][c]/(1.0E-3*snow->Dzl->co[l][r][c]*rho_w);		//[-]
				Se=(theta_liq - par->Sr*(1.0-thice))/(1.0-thice - par->Sr*(1.0-thice));
				if(Se<0) Se=0.0;
				if(Se>1) Se=1.0;
				q->co[l]=fabs(cos(slope))*kidrmax_snow((snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c])/(0.001*snow->Dzl->co[l][r][c]))*pow(Se,3.0);

				if(theta_liq>par->Sr*(1.0-thice)){		//when theta_liq is larger than the capillary retenction
					snow->w_liq->co[l][r][c]-=q->co[l]*par->Dt;
					theta_liq=snow->w_liq->co[l][r][c]/(1.0E-3*snow->Dzl->co[l][r][c]*rho_w);
					if(theta_liq<=par->Sr*(1.0-thice)){
						q->co[l]-=((par->Sr*(1.0-thice) - theta_liq)*1.0E-3*snow->Dzl->co[l][r][c]*rho_w)/par->Dt;
						snow->w_liq->co[l][r][c]=par->Sr*(1.0-thice)*1.0E-3*snow->Dzl->co[l][r][c]*rho_w;
					}else if(theta_liq > (1.0-thice) ){
						snow->w_liq->co[l][r][c]=(1.0-thice)*rho_w*1.0E-3*snow->Dzl->co[l][r][c];
						q->co[l]+=( theta_liq - (1.0-thice) )*(rho_w*1.0E-3*snow->Dzl->co[l][r][c])/par->Dt;  //liquid excess in snow goes downwards
					}
				}
			}

		}

		//melting rate
		*Mr=q->co[1]-P/par->Dt;	//[mm/s]

		//rain on snow
		*PonS=P;	//[mm]

		free_doublevector(q);



	//b. simplified case
	}else if(snow->type->co[r][c]==1){

		snow->T->co[1][r][c]=T[1];
		*Mr=1.0E+3*wl[1]/(rho_w*par->Dt);	//[mm/s]
		snow->w_liq->co[1][r][c]=0.0;
		snow->Dzl->co[1][r][c]*=(wi[1]/snow->w_ice->co[1][r][c]);
		snow->w_ice->co[1][r][c]=wi[1];

		snow->CR1->co[1]=0.0;
		snow->CR2->co[1]=0.0;
		snow->CR3->co[1]=(wi[1]-snow->w_ice->co[1][r][c])/(snow->w_ice->co[1][r][c]*par->Dt);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void iceWBsnow(long r, long c, SNOW *snow, double P, double Ta){

	long ns;
	double Dz, h;

	if(P>0){

		Dz=P*rho_w/snow->rho_newsnow->co[r][c];	//snow depth addition due to snow precipitation

		if(snow->type->co[r][c]==0){

			snow->Dzl->co[1][r][c]+=Dz;
			snow->w_ice->co[1][r][c]+=P;

		}else if(snow->type->co[r][c]==1){

			snow->Dzl->co[1][r][c]+=Dz;
			snow->w_ice->co[1][r][c]+=P;

		}else if(snow->type->co[r][c]==2){

			ns=snow->lnum->co[r][c];

			h=(c_ice*snow->w_ice->co[ns][r][c] + c_liq*snow->w_liq->co[ns][r][c])*(snow->T->co[ns][r][c] - Tfreezing) + Lf*snow->w_liq->co[ns][r][c];
			h+=(c_ice*P)*(Ta - Tfreezing);

			snow->Dzl->co[ns][r][c]+=Dz;
			snow->w_ice->co[ns][r][c]+=P;

			if(h<0.0){
				snow->T->co[ns][r][c]=Tfreezing + h/(c_ice*snow->w_ice->co[ns][r][c] + c_liq*snow->w_liq->co[ns][r][c]);
			}else{
				snow->T->co[ns][r][c]=Tfreezing;
			}

		}


	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void WBglacier(long ns, long r, long c, GLACIER *glac, double *Mr, PAR *par, double *wi, double *wl, double *T){

	long l, m;
	double Dz, thice, theta_liq;

	*Mr=0.0;

	for(l=glac->lnum->co[r][c];l>=1;l--){

		m=ns + glac->lnum->co[r][c] - l + 1;

		if(wi[m]/glac->w_ice->co[l][r][c]<1 && glac->w_ice->co[l][r][c]/(1.0E-3*glac->Dzl->co[l][r][c]*rho_i)<0.95 )
			glac->Dzl->co[l][r][c]*=(wi[m]/glac->w_ice->co[l][r][c]);

		glac->T->co[l][r][c]=T[m];
		glac->w_ice->co[l][r][c]=wi[m];
		glac->w_liq->co[l][r][c]=wl[m];

		Dz=1.0E-3*glac->Dzl->co[l][r][c];				//[m]
		thice=glac->w_ice->co[l][r][c]/(Dz*rho_i);	//[-]
		theta_liq=glac->w_liq->co[l][r][c]/(Dz*rho_w);	//[-]

		if(thice>0.95){
			glac->Dzl->co[l][r][c]=1000.0*glac->w_ice->co[l][r][c]/(rho_i*par->snow_maxpor);
			Dz=1.0E-3*glac->Dzl->co[l][r][c];
			thice=par->snow_maxpor;
			theta_liq=glac->w_liq->co[l][r][c]/(Dz*rho_w);
		}

		if(theta_liq>par->Sr_glac*(1.0-thice)){
			*Mr+=( rho_w*Dz*(theta_liq-par->Sr_glac*(1.0-thice))/par->Dt );	//in [mm/s]
			glac->w_liq->co[l][r][c]=rho_w*Dz*par->Sr_glac*(1.0-thice);
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void glac2snow(long r, long c, SNOW *snow, GLACIER *glac, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax){

	long l;
	double D, h;

	D=DEPTH(r, c, glac->lnum, glac->Dzl);

	//if ice is less thick than Dmin[1], ice is considered as snow
	if(D<Dmin->co[1] && D>0.1){
		h=(c_ice*glac->w_ice->co[1][r][c] + c_liq*glac->w_liq->co[1][r][c])*(glac->T->co[1][r][c]-Tfreezing) + Lf*glac->w_liq->co[1][r][c];

		for(l=2;l<=glac->lnum->co[r][c];l++){
			glac->Dzl->co[1][r][c]+=glac->Dzl->co[l][r][c];
			glac->w_liq->co[1][r][c]+=glac->w_liq->co[l][r][c];
			glac->w_ice->co[1][r][c]+=glac->w_ice->co[l][r][c];
			if(glac->T->co[l][r][c]>-98.999) h+=(c_ice*glac->w_ice->co[l][r][c] + c_liq*glac->w_liq->co[l][r][c])*(glac->T->co[l][r][c]-Tfreezing) + Lf*glac->w_liq->co[l][r][c];
			glac->Dzl->co[l][r][c]=0.0;
			glac->w_liq->co[l][r][c]=0.0;
			glac->w_ice->co[l][r][c]=0.0;
			glac->T->co[l][r][c]=-99.0;
		}

		if(snow->lnum->co[r][c]==0){
			snow->lnum->co[r][c]=1;
			snow->type->co[r][c]=2;
		}else{
			h+=(c_ice*snow->w_ice->co[1][r][c] + c_liq*snow->w_liq->co[1][r][c])*(snow->T->co[1][r][c]-Tfreezing) + Lf*snow->w_liq->co[2][r][c];
		}

		snow->Dzl->co[1][r][c]+=glac->Dzl->co[1][r][c];
		snow->w_ice->co[1][r][c]+=glac->w_ice->co[1][r][c];
		snow->w_liq->co[1][r][c]+=glac->w_liq->co[1][r][c];
		if(h<0){
			snow->T->co[1][r][c]=Tfreezing + h/(c_ice*snow->w_ice->co[1][r][c]+c_liq*snow->w_liq->co[1][r][c]);
		}else{
			snow->T->co[1][r][c]=0.0;
		}

		glac->Dzl->co[1][r][c]=0.0;
		glac->w_liq->co[1][r][c]=0.0;
		glac->w_ice->co[1][r][c]=0.0;
		glac->T->co[1][r][c]=-99.0;	//novalue
		glac->lnum->co[r][c]=0;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void output_snow(SNOW *snow, double **Z, PAR *par){

	long r, c, i, l, nr, nc;
	double D, DW, n=(par->output_snow*3600.0)/(par->Dt);

	nr=snow->T->nrh;
	nc=snow->T->nch;

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(Z[r][c]!=UV->V->co[2]){

				if(par->snowtrans==1){

					DW=snow->Wsalt->co[r][c]+snow->Wsalt->co[r][c]+snow->Wsubl->co[r][c]+snow->Wsubgrid->co[r][c];

					if(par->output_snow>0){
						snow->Wtot_cum->co[r][c]+=DW;
						snow->Wsalt_cum->co[r][c]+=snow->Wsalt->co[r][c];
						snow->Wsusp_cum->co[r][c]+=snow->Wsusp->co[r][c];
						snow->Wsubl_cum->co[r][c]+=snow->Wsubl->co[r][c];
						snow->Wsubgrid_cum->co[r][c]+=snow->Wsubgrid->co[r][c];
					}

					if(par->state_pixel==1){
						for(i=1;i<=par->chkpt->nrh;i++){
							if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){
								snow->out_bs->co[1][i]+=snow->Wsalt->co[r][c];
								snow->out_bs->co[2][i]+=snow->Wsalt->co[r][c];
								snow->out_bs->co[3][i]+=snow->Wsusp->co[r][c];
								snow->out_bs->co[4][i]+=snow->Wsusp->co[r][c];
								snow->out_bs->co[5][i]+=snow->Wsubl->co[r][c];
								snow->out_bs->co[6][i]+=snow->Wsubl->co[r][c];
								snow->out_bs->co[7][i]+=snow->Wsubgrid->co[r][c];
								snow->out_bs->co[8][i]+=snow->Wsubgrid->co[r][c];
								snow->out_bs->co[9][i]+=DW;
								snow->out_bs->co[10][i]+=DW;
							}
						}
					}
				}

				if(par->output_snow>0){
					D=0.0;
					for(l=1;l<=snow->lnum->co[r][c];l++){
						//D+=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c]);
						D+=snow->Dzl->co[l][r][c];
					}
					if(D>snow->max->co[r][c]) snow->max->co[r][c]=D;
					snow->average->co[r][c]+=D/n;
				}
			}
		}
	}
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void output_map_plots(long r, long c, PAR *par, double t, long n, ENERGY *egy, METEO *met, SNOW *snow, double Hg, double LEg, double Hv, double LEv, double SWin,
	double SWg, double SWv, double LWin, double LWg, double LWv, double Ts, double Tg, double Tv){

	long j, d, M, y, h, m, N;
	double tmin, tmax, JD;

	if(par->JD_plots->co[1]!=0){

		date_time(t, par->year0, par->JD0, 0.0, &JD, &d, &M, &y, &h, &m);

		for(j=1;j<=par->JD_plots->nh;j++){

			//tmin=get_time( (double)(par->JD_plots->co[j]-1), y, par->JD0, par->year0 );
			tmin=0.0;
			tmax=get_time( (double)(par->JD_plots->co[j]  ), y, par->JD0, par->year0 );
			if(fmod(tmax-tmin,n*par->Dt)!=0.0){
				N=floor(tmax-tmin/(n*par->Dt))+1;
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
				if(par->micromet1==1){
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

	long l;
	short occ=0;

	for(l=1;l<=n;l++){
		if(e->co[l]!=e->co[l] && e->co[1]>50) occ=1;
	}

	if(occ==1){
		printf("NOvalue in Energy Balance: r:%ld c:%ld\n",r,c);
		printf("l:%ld adi:--- ad:%f ads:%f b:%f e:%f T:%f mf:%ld\n",1,ad->co[1],ads->co[1],b->co[1],e->co[1],T[1],mf->co[1]);
		for(l=2;l<=n-1;l++){
			printf("l:%ld adi:%f ad:%f ads:%f b:%f e:%f T:%f mf:%ld\n",l,adi->co[l-1],ad->co[l],ads->co[l],b->co[l],e->co[l],T[l],mf->co[l]);
		}
		printf("l:%ld adi:%f ad:%f ads:--- b:%f e:%f T:%f mf:%ld\n",n,adi->co[n-1],ad->co[n],b->co[n],e->co[n],T[n],mf->co[n]);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, SNOW *snow){

	double z=0.0, res=R, rho, k;
	long l, m;

	for(l=snow->lnum->co[r][c];l>=1;l--){
		m=snow->lnum->co[r][c]-l+1;
		z+=0.001*snow->Dzl->co[l][r][c];
		rho=(snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c])/(0.001*snow->Dzl->co[l][r][c]);
		k=rho/3.0+50.0;
		//k*=1.E10;
		frac->co[m]=res-R*exp(-k*z);
		res=R*exp(-k*z);
		//printf("l:%ld res:%f R:%f k:%f %f %f %f\n",l,res,R,k,snow->w_ice->co[l][r][c],snow->w_liq->co[l][r][c],snow->Dzl->co[l][r][c]);
		if(l==1){
			frac->co[m]+=res;
			res=0.0;
		}
	}

	frac->co[snow->lnum->co[r][c]+1]=res;

}



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void root(double d, double *D, double *root_fraction){

	long l;
	double z=0.0;

	for(l=1;l<=Nl;l++){
		z+=D[l];
		if(d>z){
			root_fraction[l]=D[l]/d;
		}else{
			if(d>(z-D[l])){
				root_fraction[l]=(d-(z-D[l]))/d;
			}else{
				root_fraction[l]=0.0;
			}
		}
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

void update_roughness_soil(double z0, double d0, double z0_z0t, double snowD, double thres, double z0snow, double *z0_ris, double *d0_ris, double *z0_z0t_ris){

	if(snowD>thres){
		*z0_ris=z0snow;
		*d0_ris=0.0;
		*z0_z0t_ris=0.0;
	}else{
		*z0_ris=z0;
		*d0_ris=d0-0.001*snowD;
		*z0_z0t_ris=z0_z0t;
	}

}

void update_roughness_veg(double z0, double d0, double hc, double fsnow, double thres, double zmu, double zmt, double *z0_ris, double *d0_ris, double *hc_ris){

	*z0_ris = 0.9*(1.0-fsnow)*z0 + 0.1*z0;
	*d0_ris = 0.9*(1.0-fsnow)*d0 + 0.1*d0;
	*hc_ris = 0.9*(1.0-fsnow)*hc + 0.1*hc;

	if(*hc_ris>zmu) *hc_ris=0.99*zmu;
	if(*hc_ris>zmt) *hc_ris=0.99*zmt;
	if(*hc_ris-(*d0_ris)<(*z0_ris)) *hc_ris=1.01*(*z0_ris)+(*d0_ris);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void EnergyFluxes(double Tg, long r, long c, long n, double zmu, double zmT, double z0s, double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v, double Ta,
	double Qa, double P, double LR, double psi, double ice, double e, double fc, double LAI, double Wc, double fwet, double fsnow, double *theta, double **soil,
	double *land, double *root, PAR *par, DOUBLEVECTOR *ftcl, DOUBLEVECTOR *rep, double SWin, double LWin, double SWv, double *LW, double *H, double *dH_dT, double *E,
	double *dE_dT, double *LWv, double *Hv, double *Ev, double *Evt, double *Tv, double *Qv, double *Ts, double *Qs, long cont, double *Hg0, double *Hg1, double *Eg0,
	double *Eg1, double *rh, double *rv, double *rc, double *rb, double *rh_ic, double *rv_ic, double *Qg){

	double Hg, dHg_dT, Eg, dEg_dT, LWg;
	double dQgdT, rm;

	//initalization
	initialize_doublevector(ftcl,0.0);
	*H=0.0;	*dH_dT=0.0;
	*E=0.0;	*dE_dT=0.0;
	*LW=0.0;
	*Hv=0.0; *Ev=0.0; *Evt=0.0; *LWv=0.0;

	//thermodynamical calculations
	sat_spec_humidity(Qg, &dQgdT, red_evap(n, psi, Tg), Tg, P);

	//TURBULENT FLUXES IN CASE OF BARE SOIL
	*Hg0=0.0;
	*Eg0=0.0;
	if(fc<1){
		if(cont<=10){
			aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, LR, rep, &rm, rh, rv, par->state_turb, par->monin_obukhov);
		}else{
			aero_resistance(zmu, zmT, z0s, d0s, rz0s, v, Ta, Tg, Qa, *Qg, P, LR, rep, &rm, rh, rv, par->state_turb, 4);
		}
		if(*Qg>Qa && n==0) *rv=(*rv)+exp(8.206-4.255*(ice+theta[1])/soil[jsat][1]);
		turbulent_fluxes(*rh, *rv, P, Ta, Tg, Qa, *Qg, dQgdT, &Hg, &dHg_dT, &Eg, &dEg_dT);

		*H+=(1.0-fc)*Hg;	*dH_dT+=(1.0-fc)*dHg_dT;
		*E+=(1.0-fc)*Eg;	*dE_dT+=(1.0-fc)*dEg_dT;

		*LW+=(1.0-fc)*( e*LWin-SB(Tg) );

		*Ts=Tg;

		if(Eg!=Eg) printf("Error in bare soil %ld %ld\n",r,c);

		*Hg0=Hg;
		*Eg0=Eg;
	}


}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void PointEnergyBalance(long r, long c, long ns, long ng, double zmu, double zmT, double z0s, double d0s, double rz0s, double z0v, double d0v, double rz0v, double hveg, double v,
	double Ta, double Qa, double P, double LR, double eps, double fc, double LAI, double Wc, double fwet, double fsnow, double *theta, double *land, double *root,
	PAR *par, DOUBLEVECTOR *ftcl, DOUBLEVECTOR *turb_rep, double SWin, double LWin, double SWv, double *LW, double *H, double *E, double *LWv, double *Hv, double *Ev,
	double *Evt, double *Tv, double *Ts, double *Qs, short sy, SOIL *sl, double Hadd, double *SW, double t, double Dt, double *k, double *C, double *D, double *wi,
	double *wl, short sntype, double *dw, double *T, double *Hg0, double *Hg1, double *Eg0, double *Eg1){

	short occurs=0;
	long l, cont1, cont2=0, cont3=0, n=Nl+ns+ng;
	double Tg, dH_dT, dE_dT, EB, dEB_dT;
	double AD0, B0, L0;
	SHORTVECTOR *mf;
	DOUBLEVECTOR *ad0, *ads0, *adi0, *b0, *ad1, *ads1, *adi1, *b1, *e, *T0;
	DOUBLEVECTOR *e0g, *e1g, *Cg, *Tpg, *dwg;
	double rh, rv, rc, rb, rh_ic, rv_ic, Qv, Qg;
	double Zboundary=land[jzb], Tboundary=land[jtb];


	//TO CANCEL
	//********************************************************************************************************
	long ri=0, ci=0;
	//********************************************************************************************************

	//ALLOCATE
	mf=new_shortvector(n);

	ad0=new_doublevector(n);
	ads0=new_doublevector(n-1);
	adi0=new_doublevector(n-1);
	b0=new_doublevector(n);
	ad1=new_doublevector(n);
	ads1=new_doublevector(n-1);
	adi1=new_doublevector(n-1);
	b1=new_doublevector(n);
	e=new_doublevector(n);
	T0=new_doublevector(n);

	e0g=new_doublevector(Nl);
	e1g=new_doublevector(Nl);
	Cg=new_doublevector(Nl);
	Tpg=new_doublevector(Nl);
	dwg=new_doublevector(Nl);

	//TO CANCEL
	//********************************************************************************************************
	if(r==ri && c==ci){	for(l=1;l<=Nl+ns+ng;l++){	printf("\n.%ld T:%f",l,T[l]);	}	}
	//********************************************************************************************************

	Tg=T[1];
	for(l=1;l<=n;l++){ T0->co[l]=T[l]; }

	EnergyFluxes(Tg, r, c, ns+ng, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, P, LR, sl->P->co[1][r][c], sl->thice->co[1][r][c], eps, fc, LAI, Wc, fwet, fsnow,
			theta, sl->pa->co[sy], land, root, par, ftcl, turb_rep, SWin, LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, Ev, Evt, Tv, &Qv, Ts, Qs, cont3, Hg0, Hg1, Eg0,
			Eg1, &rh, &rv, &rc, &rb, &rh_ic, &rv_ic, &Qg);
	EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
	dEB_dT = -eps*dSB_dT(Tg) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;

	coeff(n, ns+1, ad0->co, ads0->co, adi0->co, b0->co, k, C, D, T, KNe*EB, SW, Zboundary, Tboundary, Dt);
	AD0=ad0->co[1];	B0=b0->co[1];

	L0=turb_rep->co[2];

	//TO CANCEL
	//********************************************************************************************************
	if(r==ri && c==ci){printf("\nTg:%f T1:%f Ta:%f Tv:%f H:%f E:%e dH_dT:%f EB:%f SW:%f LE:%f LW:%f \n",Tg,T[1],Ta,*Tv,*H,*E,dH_dT,EB+SW[1],SW[1],latent(Tg,Levap(Tg))*(*E),*LW);}
	//********************************************************************************************************

	do{

		Tg=T[1];
		ad0->co[1]=AD0;	b0->co[1]=B0;
		for(l=1;l<=n;l++){ mf->co[l]=0; }

		if(cont2>0){
			EnergyFluxes(Tg, r, c, ns+ng, zmu, zmT, z0s, d0s, rz0s, z0v, d0v, rz0v, hveg, v, Ta, Qa, P, LR, sl->P->co[1][r][c], sl->thice->co[1][r][c], eps, fc, LAI,
				Wc, fwet, fsnow, theta, sl->pa->co[sy], land, root, par, ftcl, turb_rep, SWin, LWin, SWv, LW, H, &dH_dT, E, &dE_dT, LWv, Hv, Ev, Evt, Tv, &Qv, Ts, Qs,
				cont3, Hg0, Hg1, Eg0, Eg1, &rh, &rv, &rc, &rb, &rh_ic, &rv_ic, &Qg);
			EB = *LW - (*H) - latent(Tg,Levap(Tg))*(*E) + Hadd;
			dEB_dT = -eps*dSB_dT(Tg) - dH_dT - latent(Tg,Levap(Tg))*dE_dT;
		}

		if(L0*turb_rep->co[2]<0) cont3++;
		L0=turb_rep->co[2];

		ad0->co[1]-=(1.0-KNe)*dEB_dT;
		b0->co[1]+=(1.0-KNe)*(EB-dEB_dT*Tg);

		tridiag(1, r, c, n, adi0, ad0, ads0, b0, e);
		check_errors(r, c, n, adi0, ad0, ads0, b0, e, T, mf);

		occurs=0;
		update_coeff(r, c, ns+ng, Nl, &occurs, ad1->co, ads1->co, adi1->co, b1->co, mf->co, ad0->co, ads0->co, adi0->co, b0->co, T0->co, e->co, wi, wl, C,
			e0g, e1g, Cg, Tpg, dwg, sl, sl->pa->co[sy][jsf], par, Dt, T, dw);

		//TO CANCEL
		//********************************************************************************************************
		if(r==ri && c==ci){
			for(l=1;l<=n;l++){
				if(l<n) printf("%ld ad:%f adi:%f ads:%f b:%f e:%f\n",l,ad0->co[l],adi0->co[l],ads0->co[l],b0->co[l],e->co[l]);
				if(l<n) printf("%ld ad:%f adi:%f ads:%f b:%f e:%f mf:%ld\n",l,ad1->co[l],adi1->co[l],ads1->co[l],b1->co[l],e->co[l],mf->co[l]);
			}
		}
		//********************************************************************************************************

		cont1=0;

		do{

			tridiag(1, r, c, n, adi1, ad1, ads1, b1, e);
			check_errors(r, c, n, adi1, ad1, ads1, b1, e, T, mf);

			//TO CANCEL
			//********************************************************************************************************
			if(r==ri && c==ci){
				for(l=1;l<=n;l++){
					if(l<n) printf("%ld ad:%f adi:%f ads:%f b:%f e:%f mf:%ld\n",l,ad1->co[l],adi1->co[l],ads1->co[l],b1->co[l],e->co[l],mf->co[l]);
				}
			}
			//********************************************************************************************************

			assign_values(ns+ng, mf->co, e->co, wi, wl, dw, T);

			occurs=0;
			update_coeff(r, c, ns+ng, Nl, &occurs, ad1->co, ads1->co, adi1->co, b1->co, mf->co, ad0->co, ads0->co, adi0->co, b0->co, T0->co, e->co, wi, wl, C,
				e0g, e1g, Cg, Tpg, dwg, sl, sl->pa->co[sy][jsf], par, Dt, T, dw);

			cont1++;

		}while(occurs!=0 && cont1<10);

		check_continuity(n, dw, wi, wl);

		cont2++;

		//TO CANCEL
		//********************************************************************************************************
		//if(cont1==20) printf("cont1:%ld r:%ld c:%ld\n",cont1,r,c);
		//if(cont1==20) stop_execution();
		if(r==ri && c==ci)for(l=1;l<=Nl+ns+ng;l++){ printf("\nb) cont1:%ld %ld T:%f dw:%f mf:%ld wi:%f wl:%f\n",cont1,l,T[l],dw[l],mf->co[l],wi[l],wl[l]); }
		if(r==ri && c==ci) printf("\nTg:%f T1:%f Ta:%f Tv:%f H:%f E:%e dH_dT:%f EB:%f SW:%f cont2:%ld cont3:%ld \n",Tg,T[1],Ta,*Tv,*H,*E,dH_dT,EB,SW[1],cont2,cont3);
		/*if((Tg>0 || T[1]>0) && ns>0){
			for(l=1;l<=Nl+ns+ng;l++){ printf("b) cont2:%ld %ld T:%f dw:%f mf:%ld wi:%f wl:%f \n",cont2,l,T[l],dw[l],mf->co[l],wi[l],wl[l]); }
			printf("r:%ld c:%ld Tg:%f T1:%f Ta:%f Tv:%f H:%f E:%e dH_dT:%f EB:%f SW:%f cont2:%ld cont3:%ld \n",r,c,Tg,T[1],Ta,*Tv,*H,*E,dH_dT,EB,SW[1],cont2,cont3);
			stop_execution();
		}*/
		//if(r==ri && c==ci) stop_execution();
		//********************************************************************************************************

	//only one iteration is temporarily done!!
	}while(cont2<0);

	//TO CANCEL
	//********************************************************************************************************
	//if(r==ri && c==ci) printf("CONVERGENCE\n");
	//if(r==ri && c==ci) stop_execution();
	for(l=1;l<=Nl+ns+ng;l++){
		if(T[l]!=T[l] || dw[l]!=dw[l]){
			printf("ERRROR IN ENERGY BALANCE r:%ld c:%ld l:%ld T:%f dw:%f wi:%f wl:%f\n",r,c,l,T[l],dw[l],wi[l],wl[l]);
			stop_execution();
		}
	}
	//********************************************************************************************************

	//Update
	*H=(*H)+dH_dT*(T[1]-Tg);
	*E=(*E)+dE_dT*(T[1]-Tg);
	*LW=(*LW)-eps*dSB_dT(Tg)*(T[1]-Tg);
	for(l=1;l<=Nl+ns+ng;l++){
		wl[l]+=dw[l];
		wi[l]-=dw[l];
	}


	//Deallocate
	free_doublevector(ad0);
	free_doublevector(ads0);
	free_doublevector(adi0);
	free_doublevector(b0);
	free_doublevector(ad1);
	free_doublevector(ads1);
	free_doublevector(adi1);
	free_doublevector(b1);
	free_doublevector(e);
	free_doublevector(T0);

	free_doublevector(e0g);
	free_doublevector(e1g);
	free_doublevector(Cg);
	free_doublevector(Tpg);
	free_doublevector(dwg);

	free_shortvector(mf);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void coeff(long n, long nlim, double *ad, double *ads, double *adi, double *b, double *k, double *C, double *D, double *T, double hup, double *h, double Zboundary, double Tboundary, double Dt){

	double Dcorr;
	long l;

	Dcorr=D[1];
	//Dcorr=0.5*(0.5*D[1]+CA*(D[1]+0.5*D[2]));

	ad[1] = 2.0*(1-KNe)*k[1]/(D[1]+D[2]) + C[1]*Dcorr/Dt;
	ads[1] = 2.0*(KNe-1)*k[1]/(D[1]+D[2]);
	b[1] = C[1]*Dcorr*T[1]/Dt + 2.0*KNe*k[1]*(T[2]-T[1])/(D[1]+D[2]) + hup + h[1];

	for(l=2;l<=n-1;l++){
		ad[l] = C[l]*D[l]/Dt + 2.0*(1-KNe)*k[l-1]/(D[l-1]+D[l]) + 2.0*(1-KNe)*k[l]/(D[l]+D[l+1]);
		ads[l] = 2.0*(KNe-1)*k[l]/(D[l]+D[l+1]);
		adi[l-1] = ads[l-1];	//symmetrical matrix
		b[l] = C[l]*D[l]*T[l]/Dt - 2.0*KNe*k[l-1]*(T[l]-T[l-1])/(D[l-1]+D[l]) + 2.0*KNe*k[l]*(T[l+1]-T[l])/(D[l]+D[l+1]);
		if(l<=nlim) b[l]+=h[l];
	}

	if(Zboundary<=0) k[n]=0.0;
	k[n]=0.0;
	ad[n] = C[n]*D[n]/Dt + 2.0*(1-KNe)*k[n-1]/(D[n-1]+D[n]) + 2.0*(1-KNe)*k[n]/(D[n]+2.0*Zboundary);
	adi[n-1] = ads[n-1];
	b[n] = C[n]*D[n]*T[n]/Dt - 2.0*KNe*k[n-1]*(T[n]-T[n-1])/(D[n-1]+D[n]) + 2.0*KNe*k[n]*(Tboundary-T[n])/(D[n]+2.0*Zboundary) + 2.0*(1-KNe)*k[n]*Tboundary/(D[n]+2.0*Zboundary);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void update_coeff(long r, long c, long n1, long n2, short *occurs, double *AD, double *ADS, double *ADI, double *B, short *mf, double *ad, double *ads, double *adi,
	double *b, double *e0, double *e, double *wi, double *wl, double *C, DOUBLEVECTOR *e0g, DOUBLEVECTOR *e1g, DOUBLEVECTOR *Cg, DOUBLEVECTOR *Tg, DOUBLEVECTOR *dwg,
	SOIL *sl, double *SFflag, PAR *par, double Dt, double *T, double *dw){

	long l, m;
	double tol=0.0;

	//for snow: changes coefficients and solves the equation
	//for soil: solves the equation and changes coefficients

	//INITIALIZE
	for(l=1;l<=n1+n2;l++){
		AD[l]=ad[l];
		B[l]=b[l];
		if(l!=n1+n2) ADS[l]=ads[l];
		if(l!=n1+n2) ADI[l]=adi[l];
	}

	//SNOW AND GLACIER
	//1.check cases
	for(l=1;l<=n1;l++){

		if(mf[l]==0){
			if(e[l]>Tfreezing && wi[l]>0.0){
				mf[l]=1;//melting
				*occurs=1;
			}else if(e[l]<Tfreezing && wl[l]>0.0){
				mf[l]=2;//freezing
				*occurs=1;
			}

		}else if(mf[l]==1){
			if(e[l]>wi[l]){
				mf[l]=10;
				*occurs=1;
			}

		}else if(mf[l]==2){
			if(-e[l]>wl[l]){
				mf[l]=20;
				*occurs=1;
			}
		}
	}

	//2.change coefficients
	for(l=1;l<=n1;l++){
		if(mf[l]==1 || mf[l]==2){
			AD[l]=Lf/Dt;
			B[l]-=ad[l]*Tfreezing;
		}else if(mf[l]==10){
			B[l]-=Lf*wi[l]/Dt;
		}else if(mf[l]==20){
			B[l]+=Lf*wl[l]/Dt;
		}
		if(l>1){
			if(mf[l-1]==1 || mf[l-1]==2){
				ADI[l-1]=0.0;
				B[l]-=adi[l-1]*Tfreezing;
			}
		}
		if(l<n1+n2){
			if(mf[l+1]==1 || mf[l+1]==2){
				ADS[l]=0.0;
				B[l]-=ads[l]*Tfreezing;
			}
		}
	}

	//SOIL FREEZING
	for(l=1;l<=n2;l++){
		e0g->co[l]=e0[l+n1];
		e1g->co[l]=e[l+n1];
		Cg->co[l]=C[l+n1];
		if(SFflag[l]==1){
			soil_freezing(l, r, c, e0g, e1g, dwg, Tg, Cg, sl, par->psimin);
			tol+=(fabs(T[l+n1]-Tg->co[l]));
			if(tol>0.1){
				*occurs=1;
				mf[l+n1]=-1;
			}
			T[l+n1]=Tg->co[l];
			dw[l+n1]=dwg->co[l];
		}else{
			T[l+n1]=e1g->co[l];
			dw[l+n1]=0.0;
		}
	}
	for(l=1;l<=n2;l++){
		if(SFflag[l]==1){
			m=l+n1;
			if(m>1){
				B[m-1]-=ads[m-1]*Tg->co[l];
				ADS[m-1]=0.0;
			}
			if(m<n1+n2){
				B[m+1]-=adi[m]*Tg->co[l];
				ADI[m]=0.0;
			}
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void assign_values(long n1, short *mf, double *e, double *wi, double *wl, double *dw, double *T){

	long l;

	//snow and glacier
	for(l=1;l<=n1;l++){
		if(mf[l]==0){
			T[l]=e[l];
			dw[l]=0.0;
		}else if(mf[l]==1 || mf[l]==2){
			T[l]=Tfreezing;
			dw[l]=e[l];
		}else if(mf[l]==10){
			T[l]=Tfreezing;
			dw[l]=wi[l];
		}else if(mf[l]==20){
			T[l]=e[l];
			dw[l]=-wl[l];
			if(T[l]>Tfreezing){
				dw[l]+=c_ice*(T[l]-Tfreezing)/Lf;
				T[l]=Tfreezing;
			}
		}
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
void find_tau_cloud_station(long i, METEO *met, PAR *par, double alpha, double E0, double A, double sky, double *tau, double *sa){
	/* Author: Stefano Endrizzi    Year:
		 * function that calculates the transmittance of the sky given measures of SW radiation
		 * Input:	i:	number of meteo station under consideration
		 * 			met: meteo data
		 * 			par: parameters
		 * 			alpha: solar altitude angle
		 * 			E0: earth-sun correction
		 * 			A:
		 * 			sky: sky view factor on the considered station
		 *
		 * Output:	tau: transmittance of the sky
		 * 			sa: sin(alpha)
		 * comment: Matteo Dall'Amico, April 2009 */
	double P, RH, T, tau_atm, kd, kd0;

	if(met->column[i-1][iPs]!=-1){/* there is a column of air pressure in the station meteo file */
		P=met->var[i-1][met->column[i-1][iPs]];
	}else{
		P=pressure(met->st->Z->co[i], Pa0, 0.0);
	}

	if(met->column[i-1][iRh]!=-1){/* there is a column of relative humidity in the station meteo file */
		RH=0.01*met->var[i-1][met->column[i-1][iRh]];
	}else{
		RH=0.5;
	}

	if(met->column[i-1][iT]!=-1){/* there is a column of air temperature in the station meteo file */
		T=met->var[i-1][met->column[i-1][iT]];
	}else{
		T=0.0;
	}

	tau_atm=atm_transmittance(alpha, P, RH, T, A, par->Vis, par->Lozone);

	if(met->column[i-1][iSWb]!=-1 && met->column[i-1][iSWd]!=-1){/* there is a column of SWglobal or SWdirect in the station meteo file */
		if(met->var[i-1][met->column[i-1][iSWb]]+met->var[i-1][met->column[i-1][iSWd]]>0){
			*sa=fmax(0.01,sin(alpha));
			kd=met->var[i-1][met->column[i-1][iSWd]]/(met->var[i-1][met->column[i-1][iSWb]]+met->var[i-1][met->column[i-1][iSWd]]);
			*tau=fmin(1.0,(met->var[i-1][met->column[i-1][iSWd]]/(Isc*E0*tau_atm*(*sa)))/(sky*kd+(1-sky)*A));
		}else{
			*tau=0.0;
		}

	}else if(met->column[i-1][iSW]!=-1){/* there is a column of SWin in the station meteo file */
		kd=0.2;
		*sa=sin(alpha);
		if(*sa<0.05) *sa=0.05;
		if(*sa<0.10) kd=(kd*(*sa-0.05)+1.0*(0.10-(*sa)))/(0.10-0.05);

		//printf("kd:%f\n",kd);

		do{
			kd0=kd;
			//SW = (1-kd(T))*Isc*T*sin + sky*Kd(T)*Isc*T*sin + (1-sky)*A*Isc*T*sin
			*tau=(met->var[i-1][met->column[i-1][iSW]]/(Isc*E0*tau_atm*(*sa)))/((1-kd)+sky*kd+(1-sky)*A);
			if(*tau>1) *tau=1.0;
			if(*tau<0) *tau=0.0;
			kd=diff2glob(*tau*tau_atm);
			if(*sa<0.10) kd=(kd*(*sa-0.05)+1.0*(0.10-(*sa)))/(0.10-0.05);
			//printf("kd:%f kd0:%f tau:%f ta:%f sa:%f alpha:%f SW:%f sky:%f\n",kd,kd0,*tau,tau_atm,*sa,alpha*180./Pi,met->var[i-1][met->column[i-1][iSW]],sky);
		}while(fabs(kd0-kd)>0.005);

		//printf("conv\n");
	}
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
