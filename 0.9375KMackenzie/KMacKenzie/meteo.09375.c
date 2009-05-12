
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



//Authors: Stefano Endrizzi and Giacomo Bertoldi
//Date: 13 November 2005
//Contents: Meteorological subroutines (included turbulent transfer)
#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "liston.h"
#include "meteo.09375.h"
#include "write_dem.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "shadows.h"
#include "tabs.h"
#include "times.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void meteo_distr(METEO *met, LISTON *liston, ENERGY *egy, WATER *wat, LAND *land, TOPO *top, SNOW *snow, double time, PAR *par){

	long i,r,c;
	double t_station;
	//INTERPOLATION OF METEO VARIABLES
	for(i=1;i<=met->st->Z->nh;i++){// for each meteo station
		time_conversion(par->JD0, par->year0, time+par->Dt, met->st->JD0->co[i], met->st->Y0->co[i], &t_station);
		t_station+=(met->st->ST->co[i]-par->ST)*3600.0;// accounts for the UTM standard time
		//printf("par->JDO=%f, par->year0=%ld, time=%f, Dt=%f, met->st->JD0=%f, met->st->Y0=%ld, meteo_distr: t_station=%f",par->JD0,par->year0,time, par->Dt,met->st->JD0->co[i], met->st->Y0->co[i],t_station);stop_execution();
		meteo_interp(met->data[i-1], met->st->Dt->co[i], t_station, met->var[i-1]);
	}

	//DISTRIBUTION OF METEROLOGICAL VARIABLES FROM MEASUREMENTS IN SOME STATIONS
	data_meteo_for_liston(met, NoV);
	if(par->micromet1==1 || par->micromet2==1 || par->micromet3==1){

		call_MicroMet(time, par->Dt, top->Z0, met, egy, snow, wat->total, land->LAI, liston, par);


	}else{

		kriging_distr(time, met->st, met->LP, -0.01, top->Z0, par->integr_scale_rain, par->variance_rain, wat->total);
		meteo_vert_distr(1, top->Z0, met, par);		//use the data of the first station

	}

	if(par->en_balance==0){
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				wat->Pn->co[r][c]=wat->total->co[r][c]/3600.0;	//from [mm/h] to [mm/s]
			}
		}
	}

}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void meteo_vert_distr(long i, DOUBLEMATRIX *Z, METEO *met, PAR *par){

	//1. Wind speed
	if(met->column[i-1][iWs]==-1){
		printf("WARNING: in met file no data of wind velocity, set at default values\n");
		met->V=0.0;
	}else{
		met->V=met->var[i-1][met->column[i-1][iWs]];
	}

	//2.Relative humidity
	if(met->column[i-1][iRh]==-1){
		printf("WARNING: in met file no data of relative humidity, set at default values\n");
		met->RH=0.0;
	}else{
		met->RH=met->var[i-1][met->column[i-1][iRh]]/100.0;
		if(met->RH>1.0) met->RH=1.0;
	}

	//3. Pressure
	if(met->column[i-1][iPs]==-1){
		vert_distr(met->Pgrid, Z, 0.0, Pa0, 0.0, (*pressure));
	}else{
		vert_distr(met->Pgrid, Z, met->st->Z->co[1], met->var[i-1][met->column[i-1][iPs]], 0.0, (*pressure));
	}

	//4. Temperature
	if(met->column[i-1][iTlr]!=-1){
		met->LapseRate=met->var[i-1][met->column[i-1][iTlr]];
	}else{
		met->LapseRate=0.006509;	//normal lapse rate
	}

	if(met->column[i-1][iT]==-1){
		printf("WARNING: in met file no data of temperature, set at default values\n");
	}else{
		vert_distr(met->Tgrid, Z, met->st->Z->co[1], met->var[i-1][met->column[i-1][iT]], met->LapseRate, (*temperature));
	}

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void kriging_distr(double t, METEO_STATIONS *met_st, float *data, float novalue, DOUBLEMATRIX *Z, double int_scale, double variance, DOUBLEMATRIX *out)

{

	long i, j=0, m, n=met_st->Z->nh, index_nodata=0, r, c;
	short *station_nodata;
	DOUBLEMATRIX *krigWeights;
	DOUBLEVECTOR *Nst_use, *Est_use;
	FILE *f;

	station_nodata=(short*)malloc(n*sizeof(short));

	//verifico se ci sono pioggie < 0 (nodata) e in caso rifaccio il kriging
	for(i=0;i<n;i++){
		station_nodata[i]=0;
		if(data[i]<=novalue) station_nodata[i]=1;
		index_nodata+=station_nodata[i];
	}

	initialize_doublematrix(out, 0.0);

	if(index_nodata==n){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"\nERROR: NO RAIN MEASUREMENTS for time=%10.3f s, RAIN set to 0 mm\n",t);
		fclose(f);

		for(r=1;r<=Z->nrh;r++){
			for(c=1;c<=Z->nch;c++){
				if(Z->co[r][c]!=UV->V->co[2]){
					out->co[r][c]=0.0;
				}else{
					out->co[r][c]=UV->V->co[2];
				}
			}
		}

	}else{

		krigWeights=new_doublematrix(Z->nrh*Z->nch,n-index_nodata);
		initialize_doublematrix(krigWeights,0.999999);
		Est_use=new_doublevector(n-index_nodata);
		Nst_use=new_doublevector(n-index_nodata);

		for(i=1;i<=n;i++){
			if(station_nodata[i-1]==0){
				j+=1;
				Est_use->co[j]=met_st->E->co[i];
				Nst_use->co[j]=met_st->N->co[i];
			}
		}

		ordi_kriging2(krigWeights, Est_use, Nst_use, Z, UV, int_scale, variance);

		for(r=1;r<=Z->nrh;r++){
			for(c=1;c<=Z->nch;c++){
				if(Z->co[r][c]!=UV->V->co[2]){
					// numero riga matrice di kriging (h*k)
					m=(r-1)*Z->nch+c;
					j=0;
					out->co[r][c]=0.0;
					for(i=1;i<=n;i++){
						if(station_nodata[i-1]==0){
							j++;
							out->co[r][c]+=krigWeights->co[m][j]*data[i-1];
						}
					}
				}else{
					out->co[r][c]=UV->V->co[2];
				}
			}
		}

		free_doublematrix(krigWeights);
		free_doublevector(Est_use);
		free_doublevector(Nst_use);

	}

	free(station_nodata);
}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void vert_distr(DOUBLEMATRIX *V, DOUBLEMATRIX *Z0, double Z_st, double V_st, double gamma, double (*f)(double a, double b, double c)){

	long r,c;

	for(r=1;r<=Z0->nrh;r++){
		for(c=1;c<=Z0->nch;c++){
			if(Z0->co[r][c]!=UV->V->co[2]){
				V->co[r][c]=(*f)(Z0->co[r][c]-Z_st, V_st, gamma);
			}else{
				V->co[r][c]=UV->V->co[2];
			}
		}
	}

}

/*==================================================================================================================*/
double pressure(double Dz, double P0, double gamma){

	double P;
	P=P0*exp(-Dz*0.00013);
	return(P);
}

/*==================================================================================================================*/
double temperature(double Dz, double T0, double gamma){

	double T;
	T=(T0+tk)*exp(-(gamma/(T0+tk))*Dz)-tk;
	return(T);
}




/*==================================================================================================================*/
/*==================================================================================================================*/
/*2. subroutine SHORTWAVE_RADIATION*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void shortwave_radiation(long r, long c, double alpha, double direction, double E0, short shadow, double sky, double tau_cloud, double sa, double slope, double aspect,
	double tau_atm, double *met_data, long *met_col, double sky_st, double A, double *SWbeam, double *SWdiff, double *cosinc, LONGMATRIX *nDt_shadow,
	LONGMATRIX *nDt_sun)

{

	double kd;

	if(alpha>0){ //DI GIORNO

		/* coseno dell'angolo di incidenza */
		*cosinc=cos(slope)*sin(alpha) + sin(slope)*cos(alpha)*cos(-aspect+direction);
		//*cosinc=sin(alpha);

		/* nel caso vi sia ombra propria (self shadow)*/
		if(*cosinc<=0.0) shadow=1;

		//direct and diffuse radiation
		kd=diff2glob(tau_cloud*tau_atm);
		*SWbeam=(1-kd)*Isc*E0*tau_cloud*tau_atm*(*cosinc);
		*SWdiff=sky*kd*Isc*E0*tau_cloud*tau_atm*sin(alpha) + (1-sky)*A*Isc*E0*sin(alpha);

		//shadows
		*SWbeam*=(1-shadow);

		if(shadow==1) nDt_shadow->co[r][c]+=1;
		nDt_sun->co[r][c]+=1;

	}else{ /* se di notte */

		*cosinc=0;

		*SWbeam=0.0;
		*SWdiff=0.0;

	}

}



/*==================================================================================================================*/
/*==================================================================================================================*/
/*4. subroutine SHADOW_N*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void shadow_n(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow)
/*	Author: Year:
 * 	Function that calculates if each pixel (matrix of shadows) is in shadow or at sun given
 *  a solar elevation angle and a solar azimuth.
 *  Inputs:	top 		elevation DTM (m)
 *  		alpha:		solar elevation angle (radiants)
 *  		direction:	solar azimuth angle (from N clockwise, radiants)
 *  		point: flag indicating whether the simulation is a point (=1) or a distributed simulation
 *  Outputs:shadow: 	shadows matrix (1 shadow, 0 sun)
*/
{
long quadrata;
double beta;
long i,j;
double r2d=180.0/Pi;	//from rad to degree

initialize_shortmatrix(shadow,0); /* initialized as if it was always NOT in shadow*/

/* ###################################################################   */
if (point==1){/* PUNCTUAL (1D) simulation */

	for(i=1;i<=Nr;i++){
		for(j=1;j<=Nc;j++){
			if(top->Z0->co[i][j]!=UV->V->co[2]) shadow->co[i][j]=shadows_point(top->horizon_height[i-1][j-1], alpha*r2d, direction*r2d, 0.0, 0.0);
		}
	}

/* ###################################################################   */
}else{ /* DISTRIBUTED (2D) simulation */


	/**======== CALCOLO DELLE OMBRE ===============================*/
	quadrata=2*(Nc+Nr);

 	/*  Chiama Orizzonte# in geotoplib.c
		Inputs:  	U->co[1]: dim. pixel (funziona solo per pixel quadrati)
 	    			quadrata: dimensione matrice
 	    			alpha: altezza solare
 	    			Z0: matrice elevazioni
 	    			curv: matrice curvature
 	     			beta: azimuth +#Pi/4
 	     			novalue: novalue per Z0
 	    Outputs:	shadow: matrice ombre (1 ombra 0 sole) */

	if(direction>=0. && direction<=Pi/4.){
		beta=direction;
		Orizzonte1(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi/4. && direction<=Pi/2.){
		beta=(Pi/2.-direction);
		Orizzonte2(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi/2. && direction<=Pi*3./4.){
		beta=(direction-Pi/2.);
		Orizzonte3(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*3./4. && direction<=Pi){
		beta=(Pi-direction);
		Orizzonte4(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi && direction<=Pi*5./4.){
		beta=(direction-Pi);
		Orizzonte5(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*5./4. && direction<=Pi*3./2.){
		beta=(Pi*3./2.-direction);
		Orizzonte6(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*3./2. && direction<=Pi*7./4.){
		beta=(direction-Pi*3./2.);
		Orizzonte7(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*7./4. && direction<2.*Pi){
		beta=(2.*Pi-direction);
		Orizzonte1(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);
		//error!!!
	}
} /* end of 2D simulation case*/


}










/*==================================================================================================================*/
/*==================================================================================================================*/
/*8. subroutine INTERPOLA_METEO*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void meteo_interp(double **data, double Dt, double t, double *out)
	/* Author:    Year:
	 * function that linearly interpolates the meteo variables given the current time
	 * Input:	t:		current time [second] based on a new origin given by the beginning of the first datum of the meteo station
	 * 			Dt:		time step of the meteo station
	 * 			data:	matrix of the meteo data
	 * Output:	out:	vector of interpolated value of the meteo data for the current time
	 * comment: Matteo Dall'Amico, April 2009 */
{

	long i,j,n;

	i=floor(t/Dt);	//previous instant
	n=dim2(data);	//total number of rows in the data matrix

	//printf("time:%f i0:%ld i1:%ld n:%ld Dt=%f\n",t,i,i+1,n,Dt);stop_execution();

	if(i<0){
		t_error("ERROR 1 in the met data!!");

	}else if(i>n-1){// exceeded the dimension of the data matrix
		t_error("ERROR 2 in the met data!!");

	}else if(i==n-1){// we are at the end of the data matrix => the last value is given
		for(j=0;j<dim1(data[j]);j++){
			out[j]=data[i-1+1][j];
		}

	}else{
		for(j=0;j<dim1(data[j]);j++){
			out[j]=data[i-1+1][j]+(data[i-1+2][j]-data[i-1+1][j])*(t-i*Dt)/Dt;
			//printf("j:%ld out:%f\n",j,out[j]);
		}
	}

}




/*==================================================================================================================*/
/*==================================================================================================================*/
/*9. subroutine PART_SNOW*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void part_snow(double prec_total, double *prec_rain, double *prec_snow, double temperature, double t_rain, double t_snow)

/*	Partiziona la precipitazione in pioggia e neve a seconda della T dell'aria (Tarboton)
	Inputs:	prec_total: 	matrix of total precipitation
			temperature:	matrix of air temperature
			t_rain:			temperature above wich all precipitation is rain
			t_snow:			temperature below wich all precipitation is snow
	Outputs:prec_snow:		matrix of solid precipitation
			prec_rain:		matrix of liquid precipitation
*/

{

if(temperature<t_snow){
	*prec_snow=prec_total;
	*prec_rain=0.0;
}else if(temperature>t_rain){
	*prec_snow=0.0;
	*prec_rain=prec_total;
}else{
	*prec_snow=prec_total*(t_rain-temperature)/(t_rain-t_snow);
	*prec_rain=prec_total*(temperature-t_snow)/(t_rain-t_snow);
}

}






/*==================================================================================================================*/
/*==================================================================================================================*/
/*10. subroutine SAT_VAP_PRESSURE*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void sat_vap_pressure(double *p, double *dp_dT, double T, double P)
	/* Author:  Stefano Endrizzi  Year:
	* function calculates the saturated water vapour pressure
	* Input:
	* 		P:	air pressure [mbar]
	* 		T: air temperature [ûC]
	* Output:
	* 		p: saturated vapour pressure
	* 	dp_dT: derivative of the vapour pressure with respect to Temperature
	* comment: Matteo Dall'Amico, May 2009 */
{	//water vapour pressure p(T,P)
//pressure in [mbar] - p vapour pressure - P atmospheric pressure, temperature in [C]
	double A, b, c;
	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;
	*p=A*exp(b*T/(c+T));
	*dp_dT=*p*(b/(c+T)-b*T/pow(c+T,2.0));
}

void sat_vap_pressure_inv(double *T, double p, double P){	//temperature(vapour pressure p,P)
	double A, b, c;
	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;
	*T=c*log(p/A)/(b-log(p/A));
}

double spec_humidity(double p, double P){
	double Q;
	Q=0.622*p/(P-0.378*p);
	return(Q);
}

void sat_spec_humidity(double *Q, double *dQ_dT, double RH, double T, double P){
	double p, dp_dT;
	sat_vap_pressure(&p, &dp_dT, T, P);
	*Q=RH*spec_humidity(p, P);
	*dQ_dT=RH*dp_dT*0.622*P/pow(P-0.378*p,2.0);
}


double air_density(double T, double Q, double P){		//[kg/m3]
	double rho;
	rho=P*100/(287.04*(T+273.15))*(1- (Q * P/(0.622+0.368*Q) ) / P*(1-0.622) );
	return(rho);
}

double air_cp(double T){	//air specific heat at constant pressure [J/(kg K)] (Garrat,1992)
	double cp;
	cp=1005.00+(T+23.15)*(T+23.15)/3364.0;
	return(cp);
}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*11. subroutine TURBULENT_FLUXES*/
/*==================================================================================================================*/
/*==================================================================================================================*/


void aero_resistance(double zmu, double zmt, double z0, double d0, double z0_z0t, double v, double Ta, double T, double Qa, double Q, double P, double gmT, DOUBLEVECTOR *rep, double *rm, double *rh,
	double *rv, short state_turb, short MO){

	//calculate resistences

	//double p=pow((1000.0/P),(0.286*(1-0.23*Qa)));
	double p=1.0;
	double Tpa=(Ta+tk)*p,Tp=(T+tk)*p;	//potential temperatures

	if(state_turb==0){
		Lewis(zmu, zmt, d0, z0, z0_z0t, Tpa, Tp, v, rm, rh, rv, rep);
	}else if(state_turb==1){
		Businger(MO, zmu, zmt, d0, z0, v, 0.5*T+0.5*Ta, T-Ta, Q-Qa, z0_z0t, rm, rh, rv, rep);
	}else if(state_turb==2){	//catabatic flows - OERLEMANS & GRISOGONO (2002)
		*rm=1.0/(v*ka*ka/(log((zmu-d0)/z0 )*log((zmu-d0)/z0 )));
		if(p*(-gmT+0.0098)>0.0015){
			*rh=1/( 0.0004*(Ta-T)*pow(g/(tk*p*(-gmT+0.0098)*5),0.5) );
		}else{
			*rh=1/( 0.0004*(Ta-T)*pow(g/(tk*(0.0015)*5),0.5) );
		}
		*rv=*rh;
	}

}



void turbulent_fluxes(double rh, double rv, double P, double Ta, double T, double Qa, double Q, double dQdT, double *H, double *dHdT, double *E, double *dEdT){

	double rho, cp, pot;

	rho=air_density(0.5*(Ta+T), 0.5*(Qa+Q), P);
	cp=air_cp(0.5*(Ta+T));
	//pot=pow((1000.0/P),(0.286*(1-0.23*Qa)));
	pot=1.0;

	//sensible heat flux [W/m2]
	*H=cp*rho*pot*(T-Ta)/rh;
	*dHdT=cp*rho*pot/rh;
	//evaporation [kg/(s*m2)]
	*E=rho*(Q-Qa)/rv;
	*dEdT=rho*dQdT/rv;

}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

//****PSIm
double Psim(double z)
{
double x,psi;
x=pow(1.0-15.0*z,0.25);
psi=2.0*log((1.0+x)/2.0)+log((1.0+x*x)/2.0)-2.0*atan(x)+0.5*Pi;
return(psi);
}

//****Psih
double Psih(double z)
{
double x,psi;
x=pow(1.0-15.0*z,0.25);
psi=2.0*log((1.0+x*x)/2.0);
return(psi);
}

//****Zero
double Zero(double z)
{
return(0.0);
}

//****PsiHolstag&deBruin
double PsiStab(double z)
{
double psi;
psi=10.71 + 0.7*z + 0.75*(z-14.28)*exp(-0.35*z);	//Holstag&De Bruin
/*if(z<=1){	//Brutsaert
	psi=5*z;
}else{
	psi=5*(1+log(z));
}*/
return(psi);
}


//****Ratio diffuse to global radiation - Erbs et al.(1982)
double diff2glob(double a)
{
double k;
if(a<0.22){
	k=1.0-0.09*a;
}else if(a<0.80){
	k=0.9511-0.1604*a+4.388*pow(a,2.0)-16.638*pow(a,3.0)+12.336*pow(a,4.0);
}else{
	k=0.165;
}
return(k);
}



/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
double atm_transmittance(double alpha, double P, double RH, double T, double A, double Vis, double Lozone)
/* Author:  Stefano Endrizzi  Year:
	* function that calculates the atmospheric transmittance
	* Input:
	* 		alpha:	solar altitude angle []
	* 		P:	air pressure
	* 		RH: air humidity [%]
	* 		T: air temperature [ûC]
	* 		A: Albedo of sorrounding terrain
	* 		Vis:	Visibility parameter (meteorological range) [Km]: gives an idea of the turbidity of the air 5<Vis<180Km
	* 		Lozone: Ozone level [cm(NTP)]: it is the height of gaseous ozone if all the ozone in a vertical column of unit area were brought
	* 				to normal temperature and surface pressure (NTP). This level changes according to latitude and seasons. Around the equator,
	* 				total ozone sverages about 0.24cm(NTP) and the amount increases with latitude, e.g. in the polar regions, it is up to 0.46cm(NTP).
	* 				It is maximum in Spring and minimum in fall (Iqbal 1983)
	* Output:
	* 		tau_atm: atmospheric transmittance
	* comment: Matteo Dall'Amico, May 2009 */
{
	double mr,//relative optical air mass
		 ma,//relative optical air mass for local condition (pressure corrected)
		 w, // pressure- and temperature-corrected precipitable water
		 beta, // Angstrom turbidity coefficient
		 g_beta,
		 alpha_wv, // absorbtance of direct irradiance by water vapor
		 alpha_g, // absorbtance of uniformly mixex gases
		 alpha_o, // ozone absorptance
		 tau_r, // trasmittance due to Rayleigh scattering
		 tau_as, // transmittance due to scattering aerosols
		 alpha_a, // absorbtance of direct irradiance by aerosols
		 tau_r_p, // trasmittance due to Rayleigh scattering when ma=1.66*P/P0
		 tau_as_p, // trasmittance due to scattering aerosols when ma=1.66*P/P0
		 tau_atm; // output

	//double p,dp;

	//transmissivity under cloudless sky (Iqbal par. 7.5)
	mr = 1.0/(sin(alpha)+0.15*(pow((3.885+alpha*180.0/Pi),-1.253)));// eq. 5.7.2 pg. 100 Iqbal
	ma = mr*P/1013.25;// eq. 5.7.2 pg. 100 Iqbal
	w = 0.493*RH*(exp(26.23-5416.0/(T+tk)))/(T+tk); //[cm] by Leckner: see eq. 5.4.6 pg. 94 Iqbal
	beta = 0.55*(3.912/Vis-0.01162)*(0.02472*(Vis-5)+1.132); // eq. 6.6.2 pg. 119 Iqbal (5Km<Vis<180 km)
	g_beta = -0.914 + 1.909267*exp(-0.667023*beta);// eq. 7.5.8 pg. 183 Iqbal
	alpha_wv = 0.110*(pow(w*mr+6.31E-4,0.3))-0.0121;// eq. 7.5.2 pg. 182 Iqbal
	alpha_g = 0.00235*pow(126*ma+0.0129,0.26)-7.5E-4+7.5E-3*pow(ma,0.875);// eq. 7.5.3 pg. 182 Iqbal
	alpha_o = 0.045*(pow(Lozone*mr+8.34E-4,0.38))-3.1E-3;// eq. 7.5.4 pg. 182 Iqbal
	tau_r = 0.615958+0.375566*exp(-0.221185*ma);// eq. 7.5.6 pg. 183 Iqbal
	tau_as = pow(g_beta,ma);// eq. 7.5.7 pg. 183 Iqbal
	alpha_a = 0.05*tau_as;// eq. 7.5.9 pg. 183 Iqbal
	tau_r_p = 0.615958+0.375566*exp(-0.221185*1.66*P/1013.25);// eq. 7.5.14a pg. 185 Iqbal
	tau_as_p = pow(g_beta,1.66*P/1013.25);// eq. 7.5.14b pg. 185 Iqbal
	tau_atm = (1 - alpha_wv - alpha_g - alpha_o - alpha_a)*( tau_r*tau_as + 0.5*(1-tau_r) + 0.75*(1-tau_as) );// eq. 7.5.1  eq.7.5.10 and eq.7.5.11 Iqbal
	tau_atm *= (1.0 + A*(1-alpha_wv-alpha_g-alpha_o-alpha_a)*(0.5*(1-tau_r_p) + 0.25*(1-tau_as_p)));

	/*sat_vap_pressure(&p, &dp, T, P);
	tau_atm=sin(alpha)/(1.2*sin(alpha)+RH*p*(1.0+sin(alpha))*1.E-3+0.0455);*/

	//tau_atm=0.47+0.47*sin(alpha);

	return(tau_atm);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void longwave_radiation(short state, double pvap, double RH, double T, double taucloud, double *eps, double *eps_max, double *eps_min){

	double eep, fcloud;

	if(taucloud>=1){
		fcloud=0.0;
	}else{
		fcloud=pow((1.0-taucloud)/0.75,1/3.4);
		if(fcloud>1) fcloud=1.0;
	}

	if(state==0){
		*eps=1.5689*pow(pvap/(T+tk),0.176)*(1.0+0.20*pow(fcloud,2.0));	//fitting of Zongo data
		*eps_min=1.5689*pow(pvap/(T+tk),0.176);
		*eps_max=*eps_min*(1.0+0.2);

	}else if(state==1){
		*eps=1.24*pow((pvap/(T+tk)),0.14285714)*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0); //Brutsaert, 1975 + Pirazzini
		//*eps=1.24*pow((pvap/(T+tk)),0.14285714)*(1.0+0.22*pow(fcloud,3.0)); //Brutsaert, 1975 + Pirazzini
		*eps_min=1.24*pow((pvap/(T+tk)),0.14285714);
		//*eps_max=1.24*pow((pvap/(T+tk)),0.14285714)*1.22;
		*eps_max=0.979;
		//*eps=*eps_min*(1.0+0.44*RH-0.18*taucloud);

	}else if(state==2){
		*eps=1.08*(1.0-exp(-pow(pvap,(T+tk)/2016.0)))*(1.0+0.40*pow(fcloud,2.0));	//Satterlund, 1979
		*eps_min=1.08*(1.0-exp(-pow(pvap,(T+tk)/2016.0)));
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==3){
		*eps=0.765*(1.0+0.40*pow(fcloud,2.0));	//Koenig-Langlo & Augstein, 1994
		*eps_min=0.765;
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==4){
		eep=(0.7+5.95*0.00001*pvap*exp(1500/(T+tk)));	//IDSO + HODGES + PIRAZZINI
		eep=-0.792+3.161*eep-1.573*eep*eep;
		*eps=eep*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0);	//Idso, 1981
		*eps_min=eep;
		*eps_max=0.979;

	}else if(state==5){
		*eps=(0.601+5.95*0.00001*pvap*exp(1500.0/(T+tk)))*(1.0+0.40*pow(fcloud,2.0));	//Andreas and Ackley, 1982
		*eps_min=(0.601+5.95*0.00001*pvap*exp(1500.0/(T+tk)));
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==6){
		*eps=(0.23+0.484*pow(pvap/(T+tk),0.125))*(1.0+0.40*pow(fcloud,2.0));	//Konzelmann (1994)
		*eps_min=(0.23+0.484*pow(pvap/(T+tk),0.125));
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==11){
		*eps=1.24*pow((pvap/(T+tk)),0.14285714)*(1.0+0.26*fcloud); //Brutsaert, 1975+Jacobs
		*eps_min=1.24*pow((pvap/(T+tk)),0.14285714);
		*eps_max=*eps_min*(1.0+0.26);

	}else if(state==14){
		eep=(0.7+5.95*0.00001*pvap*exp(1500/(T+tk)));	//IDSO + PIRAZZINI
		*eps=eep*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0);
		*eps_min=eep;
		*eps_max=0.979;

	}else{
		t_error("Incorrect value for longwave radiation formula");
	}

	//Pirazzini et al. (2000)
	if(*eps>0.979) *eps=0.979;
	if(*eps_max>0.979) *eps_max=0.979;

}



/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
//Loius' scheme (Kot & Song, 1998)
void Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, double v, double *rm, double *rh, double *rv, DOUBLEVECTOR *w){

	double z0t, f, Rib;
	double Chn, bh, c1h, c2h, c3h, c4h, Chx, ch, Fh;
	double Cmn, bm, c1m, c2m, c3m, c4m, Cmx, cm, Fm;

	//check
	//if(zmt-zmu>0.5 || zmt-zmu<-0.5) t_error("If you use Louis' scheme, wind and temperature should be measured approximately at the same elevation on the ground");

	//roughness
	if(z0_z0t==0.0){	//rigid surface
		z0t=z0/8.5;
	}else{				//bending surface
		z0t=z0/z0_z0t;
	}

	/* CH neutrale [m/s]*/
	Cmn=ka*ka/(log((zmu-d0)/z0 )*log((zmu-d0)/z0 ));
	Chn=ka*ka/(log((zmu-d0)/z0 )*log((zmt-d0)/z0t));

	/* calcola la funzione di stabilita' secondo la trattazione semplificata di Garrat, 1992 */
	f=(zmu/zmt-z0/zmt)/pow(zmt-z0t,0.5);

	/*calculation of the Richardson'number of Bulk */
	Rib=9.81*zmt*(Ta-Ts)/(0.5*(Ta+Ts)*v*v)*f*f;

	bm=8.0;
	bh=23.0;
	if (Rib < 0.0) {
		c1m=-0.9848;
		c2m=2.5398;
		c3m=-0.2325;
		c4m=14.1727;
		if (z0/z0t<100.0) {
			c1h= -1.1790;
			c2h= -1.9256;
			c3h= 0.1007;
			c4h= 16.6796;
		}else{
			c1h= -1.0487;
			c2h= -1.0689;
			c3h= 0.0952;
			c4h= 11.7828;
		}
	}else{
		if (z0/z0t<100.0) {
			c1m=-0.4738;
			c2m=-0.3268;
			c3m=0.0204;
			c4m=10.0715;
			c1h= -0.5128;
			c2h= -0.9448;
			c3h= 0.0643;
			c4h= 10.8925;
		}else{
			c1m=-0.4613;
			c2m=-0.2402;
			c3m=0.0146;
			c4m=8.9172;
			c1h= -0.3169;
			c2h= -0.3803;
			c3h= 0.0205;
			c4h= 7.5213;
		}
	}

	Cmx=c1m*log(zmu/z0)+c2m*log(z0/z0t)+c3m*log(zmt/z0t)+c4m;
	Chx=c1h*log(zmu/z0)+c2h*log(z0/z0t)+c3h*log(zmt/z0t)+c4h;

	cm=Cmx*Cmn*bm*f*pow(pow(zmu/z0 ,1.0/3.0)-1.0,1.5);
	ch=Chx*Chn*bh*f*pow(pow(zmt/z0t,1.0/3.0)-1.0,1.5);

	if (Rib < 0.0) {
		Fm=1.0-bm*Rib/(1.0+cm*pow((-Rib),0.5));
		Fh=1.0-bh*Rib/(1.0+ch*pow((-Rib),0.5));
	}else{
		Fm=1.0/(pow(1.0+Cmx*Rib,2.0));
		Fh=1.0/(pow(1.0+Chx*Rib,2.0));
	}

	/* calcola la resistenza aereodinamica */
	*rm=1.0/(Cmn*Fm*v);
	*rh=1.0/(Chn*Fh*v);
	*rv=*rh;

	//check variables
	w->co[1]=Rib;
	w->co[2]=Fh;
}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
double cz(double zmeas, double z0, double d0, double L, double (* unstab)(double z), double (* stab)(double z)){

	double c,zeta;

	zeta=(zmeas-d0)/L;
	//printf("zmeas:%f d0:%f L:%f zeta:%f\n",zmeas,d0,L,zeta);

	if(zeta<0){
		c=log((zmeas-d0)/z0) - (*unstab)(zeta) + (*unstab)(z0/L);
		//printf("d0:%f z0:%f z0/L:%e %f %f %f %f\n",d0,z0,z0/L,log((zmeas-d0)/z0),-(*unstab)(zeta),(*unstab)(z0/L),c);
	}else{
		c=log((zmeas-d0)/z0) + (*stab)(zeta) - (*stab)(z0/L);
		//printf("d0:%f z0:%f z0/L:%e %f %f %f %f\n",d0,z0,z0/L,log((zmeas-d0)/z0),(*stab)(zeta),-(*stab)(z0/L),c);
	}

	return(c);

}

/*==================================================================================================================*/
double CZ(short state, double zmeas, double z0, double d0, double L, double (*Psi)(double z)){

	double c;

	if(state==1){			//both instability and stability considered
		c=cz(zmeas,z0,d0,L,(Psi),(*PsiStab));
	}else if(state==2){		//instability considered & stability not considered
		c=cz(zmeas,z0,d0,L,(Psi),(*Zero));
	}else if(state==3){		//instability not considered & stability considered
		c=cz(zmeas,z0,d0,L,(*Zero),(*PsiStab));
	}else if(state==4){		//both instability and stability not considered
		c=cz(zmeas,z0,d0,L,(*Zero),(*Zero));
	}else{
		t_error("Value not admitted in CV");
	}

	return(c);
}

/*==================================================================================================================*/
void Star(short a, double zmeas, double z0, double d0, double L, double u, double delta, double M, double N, double R, double *var, double *c, double *z0v,
	double (*Psi)(double z), double (*roughness)(double x, double y, double z) ){

	*z0v=z0*(*roughness)(M, N, R);
	//if(*z0v<1.0E-5) *z0v=1.0E-5;
	*c=CZ(a,zmeas,*z0v,d0,L,(Psi));
	*var=delta*ka/(*c);
	//printf("delta:%f c:%f z0:%f z0v:%f zmeas:%f\n",delta,*c,z0,*z0v,zmeas);

}
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
double roughT(double M, double N, double R){

	double b0,b1,b2,fr;

	if(M<=0.135){
		b0=1.250;
		b1=0.0;
		b2=0.0;
	}else if(M<2.5){
		b0=0.149;
		b1=-0.550;
		b2=0.0;
	}else{
		b0=0.317;
		b1=-0.565;
		b2=-0.183;
	}
	fr=R+N*exp(b0+b1*log(M)+b2*pow(log(M),2.0));

	return(fr);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
double roughQ(double M, double N, double R){

	double b0,b1,b2,fr;

	if(M<=0.135){
		b0=1.610;
		b1=0.0;
		b2=0.0;
	}else if(M<2.5){
		b0=0.351;
		b1=-0.628;
		b2=0.0;
	}else{
		b0=0.396;
		b1=-0.512;
		b2=-0.180;
	}
	fr=R+N*exp(b0+b1*log(M)+b2*pow(log(M),2.0));
	//fr=R+N*exp(-ka*(7.3*pow(M,0.25)*pow(0.595,0.5)-5));
	return(fr);

}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void Businger(short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w){

	double L, L0, cm, u_star, ch, T_star, cv, Q_star;
	double z0v, z0t, z0q;
	long cont;
	double tol;
	FILE *f;

	//first guess of Obukhov length
	if(DT<0){
		L=1.E5;
	}else{
		L=-1.E5;
	}

	cont=0;
	tol=1.0E99;

	do{

		if(cont>0) tol=10*T_star+100*u_star+1000*Q_star;

		L0=L;

		//Conductances
		Star(a, zmu, z0, d0, L, 0.0, v, 1.0, 0.0, 1.0, &u_star, &cm, &z0v, (*Psim), (*roughT));	//momentum
		if(z0_z0t==0.0){	//rigid surface
			Star(a, zmt, z0, d0, L, u_star, DT, u_star*z0/1.4E-5, 1.0, 0.0, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, u_star*z0/1.4E-5, 1.0, 0.0, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
			//printf("rigid\n");
		}else{	//bending surface
			Star(a, zmt, z0, d0, L, u_star, DT, 1.0, 0.0, 1.0/z0_z0t, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, 1.0, 0.0, 1.0/z0_z0t, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
			//printf("bending\n");
		}

		//Obukhov length
		L=-u_star*u_star*(T+tk)/(ka*g*(T_star+0.61*Q_star*(T+tk)));
		//if(L*L0<0) L=L0;

		cont++;

	}while(fabs(10*T_star+100*u_star+1000*Q_star-tol)>0.01 && cont<=100);

	//CHECK
	/*if(fabs(L-L0)>0.001){
		if(par->n_error<par->max_error){
			par->n_error++;
			f=fopen(O_0ERRORSname,"a");
			fprintf(f,"Obukhov length iteration scheme doesn't converge after %4d steps (L-L0= %15.12f, L=%10.5f) in cell [%ld,%ld]\n",cont,L-L0,L,r,c);
			fclose(f);
		}
	}*/

	if(d0>zmu || d0>zmt){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"ERROR: Displacement height greater than measurement elevations");
		fclose(f);
	}
	if(zmu<=z0 || zmt<=z0t || zmt<=z0q){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"ERROR: Elevation of sensors lower than roughness length: zmu=%f zmt=%f z0=%f z0t=%f\n",zmu,zmt,z0,z0t);
		fclose(f);
	}

	*rm=cm*cm/(ka*ka*v);
	*rh=ch*cm/(ka*ka*v);
	*rv=cv*cm/(ka*ka*v);

	w->co[1]=(double)cont;
	w->co[2]=L;

}


void Businger2(long r, long c, short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w){

	double L, L0, cm, u_star, ch, T_star, cv, Q_star;
	double z0v, z0t, z0q;
	double tol;
	long cont;
	FILE *f;

	//first guess of Obukhov length
	if(DT<0){
		L=1.E5;
	}else{
		L=-1.E5;
	}

	cont=0;
	tol=1.0E99;

	do{

		if(cont>0) tol=10*T_star+100*u_star+1000*Q_star;

		L0=L;

		//Conductances
		Star(a, zmu, z0, d0, L, 0.0, v, 1.0, 0.0, 1.0, &u_star, &cm, &z0v, (*Psim), (*roughT));	//momentum
		if(z0_z0t==0.0){	//rigid surface
			Star(a, zmt, z0, d0, L, u_star, DT, u_star*z0/1.4E-5, 1.0, 0.0, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, u_star*z0/1.4E-5, 1.0, 0.0, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
		}else{	//bending surface
			Star(a, zmt, z0, d0, L, u_star, DT, 1.0, 0.0, 1.0/z0_z0t, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, 1.0, 0.0, 1.0/z0_z0t, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
		}

		//Obukhov length
		L=-u_star*u_star*(T+tk)/(ka*g*(T_star+0.61*Q_star*(T+tk)));
		//if(L*L0<0) L=L0;

		cont+=1;

	}while(fabs(10*T_star+100*u_star+1000*Q_star-tol)>0.01 && cont<=100);

	//CHECK
	/*if(fabs(L-L0)>0.001){
		if(par->n_error<par->max_error){
			par->n_error++;
			f=fopen(O_0ERRORSname,"a");
			fprintf(f,"Obukhov length iteration scheme doesn't converge after %4d steps (L-L0= %15.12f, L=%10.5f) in cell [%ld,%ld]\n",cont,L-L0,L,r,c);
			fclose(f);
		}
	}*/

	if(d0>zmu || d0>zmt){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"ERROR: Displacement height greater than measurement elevations");
		fclose(f);
	}
	if(zmu<=z0 || zmt<=z0t || zmt<=z0q){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"ERROR: Elevation of sensors lower than roughness length: zmu=%f zmt=%f z0=%f z0t=%f\n",zmu,zmt,z0,z0t);
		fclose(f);
	}

	*rm=cm*cm/(ka*ka*v);
	*rh=ch*cm/(ka*ka*v);
	*rv=cv*cm/(ka*ka*v);

	w->co[1]=(double)cont;
	w->co[2]=L;

}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
double Levap(double T){
	double Lv;
	if(T>0.0){
		Lv=2501000.0+(2406000.0-2501000.0)/40.0*T;
	}else{
		Lv=2501000.0;
	}
	return(Lv);
}

double latent(double Ts, double Le){
	double L;
	if(Ts<0){
		L=Le+Lf;	//sublimazione
	}else{
		L=Le;		//evaporazione
	}
	return(L);
}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
short shadows_point(double **hor_height, double alpha, double azimuth, double tol_mount, double tol_flat)

/*routine that tells you whether a point is in shadow or not, depending on the solar azimuth, elevation and horizon file at that point
 * Author: Matteo Dall'Amico, May 2008
Inputs: DOUBLEMATRIX* hor_height: matrix of horizon_height at the point
		double alpha: solar altitude (degree)
		double azimuth: solar azimuth (degree)
		double tol_mount: tolerance over a mountaneaus horizon to have a reliable cloud datum (degree)
		double tol_flat: tolerance over a mountaneaus horizon to have a reliable cloud datum (degree)
Output: shad: 1=the point is in shadow, 0 the point is in sun
*/
{
	double horiz_H;// horizon elevation at a defined solar azimuth
	long i,buf,n=dim2(hor_height);
	short shad=1; //  initialized as if it was in shadow

	if(azimuth>=hor_height[n-1][0] || azimuth<hor_height[0][0])
		buf=1;
	/* compare the current solar azimuth with the horizon matrix */
	for (i=1; i<=n-1; i++){
		if(azimuth>=hor_height[i-1][0] && azimuth<hor_height[i][0] )
			buf=i+1;
		}
	horiz_H=hor_height[buf-1][1]; // horizon elevation at a particular time

	if(alpha<tol_flat){
		shad=1;
	}else if(alpha<horiz_H+tol_mount){
		shad=1;
	}else{
		shad=0;
	}

	return (shad);
}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void sun(double hour, double JD, double *alpha, double *direction, double *E0, double latitude, double longitude, double standard_time){
	/* Author:
	* Year:
	* calculations based on Iqbal
	* function that calculates the elevation angle and the azimuth angle of the sun given the time (day of the year
	* and hour of the day) and given the geographic coordinate system (latitude, longitude)
	* Input:
	* JD: Julian day
	* hour: hour of the day (hour)
	* standard time: Standard time to which all the output data are referred (difference respect UMT, in hour)
	* longitude and latitude of the point
	* Output:
	* alpha: elevation angle in radiants
	* direction: azimuth angle in radiants
	* E0: Earth-Sun correction
	*  */
	double Gamma, Delta, long_standard, Et, h;

	Gamma=2.0*Pi*JD/365.;

	/* latitude e longitude Meteo station */
	long_standard=standard_time*Pi/12.0;

	/* Earth-Sun correction */
	*E0=1.00011+0.034221*cos(Gamma)+0.00128*sin(Gamma)+0.000719*cos(2*Gamma)+0.000077*sin(2*Gamma);
	/* correction of the given hour in UTC to obtain the local hour*/
	Et=0.000075 + 0.001868*cos(Gamma) - 0.032077*sin(Gamma) - 0.014615*cos(2*Gamma) - 0.04089*sin(2*Gamma);	/*Correction for sideral day (rad)*/
	// Solar Declination (Declinazione solare)
	Delta=0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.006758*cos(2*Gamma)+0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.00148*sin(3*Gamma);

	h=hour + (longitude-long_standard)/omega + Et/omega;	/*Iqbal: formula 1.4.2*/
	if(h>=24) h-=24.0;
	if(h<0) h+=24.0;

	/**========= CALCULATION OF THE SUN MOVEMENT ===============*/
	/*!!non funziona per latitudine=90gradi e -90gradi!!*/

	/* alpha: solar elevation angle in radiants */
	*alpha=asin( sin(latitude)*sin(Delta)+cos(latitude)*cos(Delta)*cos(omega*(12-h)) );
	//printf("lat:%f\n",latitude*180/3.1415927);

	/* direction: azimuth angle (from N, clockwise) in radiants */
	if(h<=12){
		if(*alpha==Pi/2.0){	/*sun at Zenith */
			*direction=Pi/2.0;
		}else{
			*direction=Pi - acos((sin(*alpha)*sin(latitude)-sin(Delta))/(cos(*alpha)*cos(latitude)));
		}
	}else{
		if(*alpha==Pi/2.0){ /*sun at Zenith */
			*direction=3*Pi/2.0;
		}else{
			*direction=Pi + acos((sin(*alpha)*sin(latitude)-sin(Delta))/(cos(*alpha)*cos(latitude)));
		}
	}

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void data_meteo_for_liston(METEO *met, double novalue){

	long i;

	for(i=1;i<=met->st->Z->nh;i++){

		met->LT[i-1]=select_data(met->var[i-1], met->column[i-1], iT, novalue);
		met->Lrh[i-1]=select_data(met->var[i-1], met->column[i-1], iRh, novalue);
		met->Lws[i-1]=select_data(met->var[i-1], met->column[i-1], iWs, novalue);
		met->Lwd[i-1]=select_data(met->var[i-1], met->column[i-1], iWd, novalue);
		met->LP[i-1]=select_data(met->var[i-1], met->column[i-1], iPt, novalue);
	}
}

/*==================================================================================================================*/

float select_data(double *meteo_t, long *col, long cod, double novalue){

	float a;
	if(col[cod]==-1){
		a=(float)novalue;
	}else{
		a=(float)meteo_t[col[cod]];
	}
	return(a);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

