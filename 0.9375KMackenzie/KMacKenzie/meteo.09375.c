
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion MacLavagna

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 MacLavagna.
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

#include "constant.h"
#include "keywords_file.h"
#include "struct.geotop.09375.h"
#include "liston.h"
#include "meteo.09375.h"
#include "write_dem.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "shadows.h"
#include "tabs.h"

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

	long i;
	double t_station;
	//FILE *f;

	//INTERPOLATION OF METEO VARIABLES
	for(i=1;i<=met->st->Z->nh;i++){
		time_conversion(par->JD0, par->year0, time+par->Dt, met->st->JD0->co[i], met->st->Y0->co[i], &t_station);
		t_station+=(met->st->ST->co[i]-par->ST)*3600.0;
		meteo_interp(met->data[i-1], met->st->Dt->co[i], t_station, met->var[i-1]);
	}

	//DISTRIBUTION OF METEROLOGICAL VARIABLES FROM MEASUREMENTS IN SOME STATIONS
	data_meteo_for_liston(met, NoV);
	if(par->micromet1==1 || par->micromet2==1 || par->micromet3==1){

		call_MicroMet(time, par->Dt, top->Z0, met, egy, snow, wat->total, land->LAI, liston, par);

		/*if(time<60){
			f=t_fopen(join_strings(WORKING_DIRECTORY,"_meteodata.txt"),"w");
			t_fclose(f);
			printf("%s\n",join_strings(WORKING_DIRECTORY,"_meteodata.txt"));
		}

		f=fopen(join_strings(WORKING_DIRECTORY,"_meteodata.txt"),"a");
		fprintf(f,"%f,%f,%f,%f,%f\n",met->Vgrid->co[1][1],met->RHgrid->co[1][1],met->Tgrid->co[1][1],met->var[0][iSW],met->var[0][iPt]);
		fclose(f);*/

	}else{

		kriging_distr(time, met->st, met->LP, -0.01, top->Z0, par->integr_scale_rain, par->variance_rain, wat->total);
		meteo_vert_distr(1, top->Z0, met, par);		//use the data of the first station

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
void shortwave_radiation(long r, long c, double alpha, double direction, double E0, short shadow, double sky, double fcloud, double slope,
						 double aspect, double tau_atm, double *met_data, long *met_col, double sky_st, double tau_atm_st, double A, double *SWbeam, double *SWdiff,
						 double *cosinc, LONGMATRIX *nDt_shadow, LONGMATRIX *nDt_sun)

{

	double tau_cloud, kd, kd0, sa;

	if(alpha>0){ //DI GIORNO

		/* attenuazione dovuta alla copertura nuvolosa */
		tau_cloud=1.0-0.75*pow(fcloud,3.4);

		/* coseno dell'angolo di incidenza */
		*cosinc=cos(slope)*sin(alpha)+sin(slope)*cos(alpha)*cos(-aspect+direction);
		//*cosinc=sin(alpha);

		/* nel caso vi sia ombra propria (self shadow)*/
		if(*cosinc<=0.0) shadow=1;

		//direct and diffuse radiation available
		if(met_col[iSWb]!=-1 && met_col[iSWd]!=-1){
			sa=sin(alpha);
			if(sa<0.01) sa=0.01;
			//SWdiff = sky*Kd(T)*Isc*T*sin + (1-sky)*A*Isc*T*sin
			if(met_data[met_col[iSWb]]+met_data[met_col[iSWd]]>0){
				kd=met_data[met_col[iSWd]]/(met_data[met_col[iSWb]]+met_data[met_col[iSWd]]);
				tau_cloud=(met_data[met_col[iSWd]]/(Isc*E0*tau_atm_st*sa))/(sky_st*kd+(1-sky_st)*A);
			}else{
				kd=0.0;
				tau_cloud=0.0;
			}


			*SWbeam=(1-kd)*Isc*E0*tau_cloud*tau_atm*(*cosinc);
			*SWdiff=Isc*E0*tau_cloud*tau_atm*sin(alpha)*( sky*kd + (1-sky)*A );

			if(*SWbeam+*SWdiff>1500){
				printf("%r:ld c:%ld SWmeas:%f B:%f D:%f kd=%f a2st=%f tau_cloud=%f tau_atm=%f sinalpha=%f cosinc=%f sin(alpha)=%f\n",r,c,met_data[met_col[iSW]],*SWbeam,*SWdiff,kd,tau_atm_st,tau_cloud,tau_atm,sin(alpha),*cosinc,sin(alpha));
				printf("..%f %f %f %f \n",met_data[met_col[iSW]]/((1-kd)*tau_atm_st*sin(alpha)+kd*tau_atm_st*sky_st*sin(alpha)),(*cosinc*tau_atm),(sin(alpha)*sky*tau_atm),sky);
				//stop_execution();
			}

		//global radiation  available
		}else if(met_col[iSW]!=-1){
			kd=diff2glob(tau_cloud*tau_atm_st);
			sa=sin(alpha);
			if(sa<0.05) sa=0.05;
			if(sa<0.10) kd=(kd*(sa-0.05)+1.0*(0.10-sa))/(0.10-0.05);

			do{
				kd0=kd;
				//SW = (1-kd(T))*Isc*T*sin + sky*Kd(T)*Isc*T*sin + (1-sky)*A*Isc*T*sin
				tau_cloud=(met_data[met_col[iSW]]/(Isc*E0*tau_atm_st))/((1-kd)*sa+sky_st*kd*sa+(1-sky_st)*A*sa);
				kd=diff2glob(tau_cloud*tau_atm_st);
				if(sa<0.10) kd=(kd*(sa-0.05)+1.0*(0.10-sa))/(0.10-0.05);
			}while(fabs(kd0-kd)>0.05);

			*SWbeam=(1-kd)*Isc*E0*tau_cloud*tau_atm*(*cosinc);
			*SWdiff=Isc*E0*tau_cloud*tau_atm*sin(alpha)*( sky*kd + (1-sky)*A );

			if(*SWbeam+*SWdiff>1500){
				printf("%r:ld c:%ld SWmeas:%f B:%f D:%f kd=%f a2st=%f tau_cloud=%f tau_atm=%f sinalpha=%f cosinc=%f sin(alpha)=%f\n",r,c,met_data[met_col[iSW]],*SWbeam,*SWdiff,kd,tau_atm_st,tau_cloud,tau_atm,sin(alpha),*cosinc,sin(alpha));
				printf("..%f %f %f %f \n",met_data[met_col[iSW]]/((1-kd)*tau_atm_st*sin(alpha)+kd*tau_atm_st*sky_st*sin(alpha)),(*cosinc*tau_atm),(sin(alpha)*sky*tau_atm),sky);
				//stop_execution();
			}

		//calculates direct and diffuse radiation
		}else{

			kd=diff2glob(tau_cloud*tau_atm);
			*SWbeam=(1-kd)*Isc*E0*tau_cloud*tau_atm*(*cosinc);
			*SWdiff=sky*kd*Isc*E0*tau_cloud*tau_atm*sin(alpha) + (1-sky)*A*Isc*E0*sin(alpha);

		}

		//shadows
		*SWbeam*=(1-shadow);

		if(shadow==1) nDt_shadow->co[r][c]+=1;
		nDt_sun->co[r][c]+=1;

	}else{ /* se di notte */

		*cosinc=0;

		if(met_col[iSW]==-1){
			*SWbeam=0.0;
			*SWdiff=0.0;

		}else{
			*SWbeam=0.0;
			*SWdiff=0.0;
		}

	}

}







/*==================================================================================================================*/
/*==================================================================================================================*/
/*3. subroutine ALBEDO*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void albedo(double T_surface, double coszen, double prec_snow, double albedo_constant, double snowD, double *albedo, double *tausn, short pixel_type, double alpha,
			double aep, double avo, double airo, double Dt, double fcloud)

/*	Computation of albedo as function of snow-age and incidence angle ( see Dickinson page 21):
	Inputs:	T_surface: 	    snow surface temperature
			coszen:			cosine incidence angle from normal [radiant]
			prec_snow:		matrix of solid precipitation [mm]
			albedo_constant:bare ground albedo [-]
			snowD:			snow depth [mm]
			alpha:          sun height [radiant]
			aep:			albedo extintion parameter [mm]
			avo:			new snow visible band reflectance [0.85]
			airo:			new snow near infared band reflectance [0.65]
	Outputs:albedo  		albedo considering snow and alpha [-]
			tausn			non-dimensional snow age [-]
    Author: Davide Tamanini
    Date:   15-05-2004  */

{

double b, cs, cn, fzen;
double r1, r2, r3, fage, avd, avis, aird, anir, rr;

/*1.If the pixel is a lake (=11) or sea (=12):*/
if(pixel_type==11 || pixel_type==12){
	*albedo=0.05;
	if(fcloud<0.5){
		if(alpha>0.1){
			*albedo=2.2*pow((alpha*180.0/Pi),-0.9);
		}else{
			*albedo=0.4;
		}
	}else{
		if(alpha>0.1){
			*albedo=0.95*pow((alpha*180.0/Pi),-0.75);
		}else{
			*albedo=0.26;
		}
	}

/*2.Otherwise:*/
}else{
	/*2.1.With snow:*/
	if(snowD>1){

		/**2.1.a.Parameters assignment:**/
		b=2.0;/*parameter albedo dipendence from solar angle*/
		cs=0.2;/*parameter albedo dipendence from snow age*/
		cn=0.5;/*parameter albedo dipendence from snow age*/

		/**2.1.b.Calculate snow age:**/
		/*effect snow surface temperature*/
		r1=exp(5000.0*(1.0/tk-1.0/(T_surface+tk)));
		/*effect melt and refreezing*/
		r2=pow(r1,10);
		if(r2>1.0) r2=1.0;
		/*effect of dirt*/
		r3=0.03;
		/*non-dimensional snow age: 10 mm of snow precipitation restore snow age Dt(s)*/
		*tausn=(*tausn+(r1+r2+r3)*Dt*1.0E-6)*(1.0-prec_snow/10.0);

		if(*tausn<0.0) *tausn=0.0;

		if((*tausn)!=(*tausn)) printf("Tausn no value - tausn:%f P:%f Ts:%f r1:%f r2:%f r3:%f\n",*tausn,prec_snow,T_surface,r1,r2,r3);

		/**2.1.c.Calculate snow albedo:**/
		/*dipendence from solar angle*/
		if(coszen<0.5){
			fzen=1.0/b*((b+1.0)/(1.0+2.0*b*coszen)-1.0);
		}else{
			fzen=0.0;
		}
		fzen=0.0;
		/*dipendence from snow age */
		fage=*tausn/(1.0+*tausn);
		//fage=1.0-1.0/(1.0+*tausn);
		/*diffuse visible albedo*/
		avd=(1.0-cs*fage)*avo;
		/*global visible albedo*/
		avis=avd+0.4*fzen*(1.0-avd);
		/*diffuse near infared albedo*/
		aird=(1.0-cn*fage)*airo;
		/*global near infared albedo*/
		anir=aird+0.4*fzen*(1.0-aird);
		/*albedo is taken as average*/
		*albedo=(avis+anir)/2.0;

		/**2.1.d.Linear transition from snow albedo to bare ground albedo:**/
		if(snowD<aep){
			rr=(1.0-snowD/aep)*exp(-snowD*0.5/aep);
			*albedo=rr*albedo_constant+(1.0-rr)*(*albedo);
		}

		if((*albedo)!=(*albedo)){
			printf("Novalue ALBEDO with snow: albedo:%f albedoconst:%f aep:%f snowD:%f rr:%f tausn:%f fage:%f avd:%f fzen:%f coszen:%f avis:%f aird:%f anir:%f P:%f Ts:%f r1:%f r2:%f r3:%f\n",
				*albedo,albedo_constant,aep,snowD,rr,*tausn,fage,avd,fzen,coszen,avis,aird,anir,prec_snow,T_surface,r1,r2,r3);
			stop_execution();
		}

	/*2.2.Without snow:*/
	}else{

		if(alpha>0.2){/*albedo depends on solar height*/
			*albedo=albedo_constant+(1-albedo_constant)*exp(-0.1*alpha*180.0/Pi)*(1-fcloud);
		}else{
			*albedo=albedo_constant+(1-albedo_constant)*0.31*(1-fcloud/8.0);
		}

		if((*albedo)!=(*albedo)){
			printf("Novalue ALBEDO with snow: albedo:%f albedoconst:%f alpha:%f fcloud:%f\n",*albedo,albedo_constant,alpha,fcloud);
			stop_execution();
		}
	}

	if(*albedo<albedo_constant) *albedo=albedo_constant;

	if(*albedo>0.9) *albedo=0.9;
}

}




/*==================================================================================================================*/
/*==================================================================================================================*/
/*4. subroutine SHADOW_N*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void shadow_n(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow)


/*
Inputs: Z0 		matrice delle elevazioni
		curv	matrice delle curvature
		day 	giorno di osservazione
		hour	ora del giorno in questione

Outputs:alpha altezza solare in radianti
		direction: azimuth sole (da N, orario) in radianti
		shadow: matrice con le ombre (1 ombra 0 sole)
		E0:		correzione distanza Terra-Sole
*/
{
long quadrata;
double beta;
long i,j;
double r2d=180.0/Pi;	//from rad to degree

initialize_shortmatrix(shadow,0); /* initialized as if it was always NOT in shadow*/

/* ###################################################################   */
if (point==1){/* PUNCTUAL (1D) simulation */

	for(i=1;i<=top->Z0->nrh;i++){
		for(j=1;j<=top->Z0->nch;j++){
			if(top->Z0->co[i][j]!=UV->V->co[2]) shadow->co[i][j]=shadows_point(top->horizon_height[i-1][j-1], alpha*r2d, direction*r2d, 0.0, 0.0);
		}
	}

/* ###################################################################   */
}else{ /* DISTRIBUTED (2D) simulation */


	/**======== CALCOLO DELLE OMBRE ===============================*/
	quadrata=2*(top->Z0->nch+top->Z0->nrh);

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

{

	long i,j,n;

	i=floor(t/Dt);	//previous instant
	n=dim2(data);	//number of time steps in the data matrix

	if(i<0){
		t_error("ERROR 1 in the met data!!");

	}else if(i>n-1){
		t_error("ERROR 2 in the met data!!");

	}else if(i==n-1){
		for(j=0;j<dim1(data[j]);j++){
			out[j]=data[i-1+1][j];
		}

	}else{
		for(j=0;j<dim1(data[j]);j++){
			out[j]=data[i-1+1][j]+(data[i-1+2][j]-data[i-1+1][j])*(t-i*Dt)/Dt;
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

//calcola la pressione di vapore [mbar] a saturazione in dipendenza dalla temperatura [gradi Celsius]

{
	double A, b, c;

	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;

	*p=A*exp(b*T/(c+T));
	*dp_dT=*p*(b/(c+T)-b*T/pow(c+T,2.0));
}

/*==================================================================================================================*/
void sat_vap_pressure_inv(double *T, double p, double P)

//Formula inversa della precedente

{
	double A, b, c;

	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;

	*T=c*log(p/A)/(b-log(p/A));
}





/*==================================================================================================================*/
/*==================================================================================================================*/
/*11. subroutine TURBULENT_FLUXES*/
/*==================================================================================================================*/
/*==================================================================================================================*/


void turbulent_fluxes(long r, long c, double zmu, double zmt, double z0, double d0, double z0_z0t, double v, double P, double RH, double gammaT, double Ta, double Ts,
					  double *H, double *E, double *dH_dT, double *dE_dT, double *r_v, double *rho, double *Qa, DOUBLEVECTOR *report, PAR *par){


	double ea, es, de_dT, cp, p, Tpa, Tps, r_h, Qs;

	//thermodynamical properties
	sat_vap_pressure(&ea, &de_dT, Ta, P);
	*Qa=0.622*ea/(P-0.378*ea);
	*Qa*=RH;
	sat_vap_pressure(&es, &de_dT, Ts, P);
	Qs=0.622*es/(P-0.378*es);
	*rho=P*100/(287.04*(Ta+273.15))*(1- (*Qa * P/(0.622+0.368*(*Qa)) ) / P*(1-0.622) );  //densita' aria [kg/m^3]
	cp=1005.00+(Ta+23.15)*(Ta+23.15)/3364;     //calore specifico a pressione costante dell'aria [J/(kg K)] (Garrat,1992)
	p=pow((1000/P),(0.286*(1-0.23*(*Qa))));		//potential temperatures

	//p=1;
	Tps=(Ts+tk)*p;
	Tpa=(Ta+tk)*p;

	// calculate resistences
	if(par->state_turb==0){
		Lewis(zmu, zmt, d0, z0, z0_z0t, Tpa, Tps, v, &r_h, r_v, report);
	}else if(par->state_turb==1){
		Businger(par->monin_obukhov, r, c, zmu, zmt, d0, z0, v, 0.5*Tps+0.5*Tpa, Tps-Tpa, Qs-(*Qa), z0_z0t, &r_h, r_v, report);
	}else if(par->state_turb==2){	//catabatic flows - OERLEMANS & GRISOGONO (2002)
		if(p*(-gammaT+0.0098)>0.0015){
			r_h=1/( 0.0004*(Ta-Ts)*pow(g/(tk*p*(-gammaT+0.0098)*5),0.5) );
		}else{
			r_h=1/( 0.0004*(Ta-Ts)*pow(g/(tk*(0.0015)*5),0.5) );
		}
		*r_v=r_h;
	}

	//fluxes and derivatives
	*H=(*rho)*cp*(Tps-Tpa)/r_h;	//[W/m2]
	*dH_dT=(*rho)*cp*p/r_h;
	*E=(*rho)*(Qs-(*Qa))/(*r_v);	//[kg/(s*m2)]
	*dE_dT=(*rho)*de_dT*(0.622*P/pow(P-0.378*es,2.0))/(*r_v);

	//*H=(*density_air)*(*cp)*(theta_sur-theta_air)/r_h + 1.00*(theta_sur-theta_air);
	//*dH_dT=(*density_air)*(*cp)*p/r_h + 1.00*p;
	//Epot=(*density_air)*(Q_sat_surface-Q)/r_v + (es-ea*RH)*2.0/(2501000.0+Lf);
	//dEpot_dT=(*density_air)*de_dT*(0.622*P/pow(P-0.378*es,2.0))/r_v + de_dT*2.0/(2501000.0+Lf);

	if(*H!=*H || *E!=*E){
		printf("\nr=%ld c=%ld\n",r,c);
		printf("\nrho=%f cp=%f t1=%f t2=%f rh=%f rv=%f ka=%f v=%f H=%f\n E=%f",*rho,cp,Tps,Tpa,r_h,*r_v,ka,v,*H,*E);
		printf("\nzmu=%f zmt=%f d0=%f z0=%f\n",zmu,zmt,d0,z0);
		stop_execution();
	}

	//report->co[3]=z0t;
	report->co[3]=0.0;
	//report->co[4]=z0q;
	report->co[4]=0.0;
	//report->co[5]=cm;
	report->co[5]=0.0;
	//report->co[6]=ch;
	report->co[6]=0.0;
	//report->co[7]=cv;
	report->co[7]=0.0;
	report->co[8]=zmu;
	report->co[9]=zmt;
	report->co[10]=1.0/r_h;
	report->co[11]=1.0/(*r_v);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
double Eg_Epot(double r_v, double theta, double res, double sat){

	double r_srf, C;

	r_srf=r_v*(1.0-(theta-res)/(sat-res))/((theta-res)/(sat-res));

	if (r_srf>1.0E9) r_srf=1.0E9;

	C=r_v/(r_v+r_srf);
	if(theta<res+0.001) C=0.0;

	if(C!=C) printf("EGEpot not a number: C:%f r_srf:%f rv:%f theta:%f res:%f sat:%f\n",C,r_srf,r_v,theta,res,sat);

	return(C);
}

/*==================================================================================================================*/
/*==================================================================================================================*/
void canopy_fluxes(double r_v, double rho, double Tv, double Ta, double Qa, double Pa, double SWin, double fW, double *theta, double LAI, double *land, double *root,
	double *Ecp, double *dEcp_dT, double *fv, double *ft, double *ftl){

	double ev, ea, de_dT, Qv, fS, fe, fTemp, Rsmin;
	long l;

	//potential evaporation from canopy (Ecp) in [kg/(s*m2)]
	sat_vap_pressure(&ev, &de_dT, Tv, Pa);
	Qv=0.622*ev/(Pa-0.378*ev);
	*Ecp=rho*(Qv-Qa)/r_v;
	//*dEcp_dT=rho*de_dT*(0.622*Pa/pow(Pa-0.378*ev,2.0))/r_v;
	*dEcp_dT=0.0; //assumes Tv=Ta

	//WET CANOPY EVAPORATION
	if(fW>=1){
		*fv=1;	//wet canopy fraction
	}else{
		*fv=pow(fW,0.66666666);
	}

	//CANOPY TRANSPIRATION (parameters from Best (1998))
	//solar radiation [Best, (1998); Dolman et al., 1991]
	fS=SWin/(SWin+250.0)*1.25;

	//pressure deficit [Best, (1998); Dickinson et al., 1991]*/
	ea=Qa*Pa/(0.378*Qa+0.622);
	fe=1.0+(ev-ea)/40.0;	//dipendenza dal deficit di pressione di vapore

	//temperature [Best, (1998); Dickinson et al., 1991]*/
	if(Ta<=0){
		fTemp=1E-12;
	}else if(Ta>=50.0){
		fTemp=1E-12;
	}else{
		fTemp=(Ta-0.0)*(50.0-Ta)/625.0;
	}

	*ft=0.0;
	for(l=1;l<=Nl;l++){

		//water content [Wigmosta et al., (1994); Feddes et al.(1978)]
		if (theta[l]>=land[jtfc]){
			ftl[l]=1.0;
		}else if(theta[l]>land[jtwp]){
			ftl[l]=(theta[l]-land[jtwp])/(land[jtfc]-land[jtwp]);
		}else{
			ftl[l]=0.0;
		}

		//stomata resistance for each layer
		if (fS*fe*fTemp*ftl[l]<6.0E-11){
			ftl[l]=1.0E12;
		}else{
			if(LAI<=1){
				Rsmin=land[jrs];
			}else{
				Rsmin=land[jrs]/LAI;
			}
			ftl[l]=Rsmin/(fS*fe*fTemp*ftl[l]);
		}

		//transpiration fraction for each layer
		ftl[l]=root[l]*(1-(*fv))*(r_v/(r_v+ftl[l]));

		//transpiration for all the column (as a fraction of Epc)
		*ft+=ftl[l];
	}

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
double psi;
psi=0.0;
return(psi);
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
double atm_transmittance(double alpha, double P, double RH, double T, double A, double Vis, double Lozone){

	double mr, ma, w, beta, g_beta, alpha_wv, alpha_g, alpha_o, tau_r, tau_as, alpha_a, tau_r_p, tau_as_p, tau_atm;

	//transmissivity under cloudless sky (Iqbal par. 7.5)
	mr = 1.0/(sin(alpha)+0.15*(pow((3.885+alpha*180.0/Pi),-1.253)));
	ma = mr*P/1013.25;
	w = 0.493*RH*(exp(26.23-5416.0/(T+tk)))/(T+tk); //cm
	beta = 0.55*(3.912/Vis-0.01162)*(0.02472*(Vis-5)+1.132); //Vis in km (>5km)
	g_beta = -0.914 + 1.909267*exp(-0.667023*beta);
	alpha_wv = 0.110*(pow(w*mr+6.31E-4,0.3))-0.0121;
	alpha_g = 0.00235*pow(126*ma+0.0129,0.26)-7.5E-4+7.5E-3*pow(ma,0.875);
	alpha_o = 0.045*(pow(Lozone*mr+8.34E-4,0.38))-3.1E-3;
	tau_r = 0.615958+0.375566*exp(-0.221185*ma);
	tau_as = pow(g_beta,ma);
	alpha_a = 0.05*tau_as;
	tau_r_p = 0.615958+0.375566*exp(-0.221185*1.66*P/1013.25);
	tau_as_p = pow(g_beta,1.66*P/1013.25);
	tau_atm = (1 - alpha_wv - alpha_g - alpha_o - alpha_a)*( tau_r*tau_as + 0.5*(1-tau_r) + 0.75*(1-tau_as) );
	tau_atm *= (1.0 + A*(1-alpha_wv-alpha_g-alpha_o-alpha_a)*(0.5*(1-tau_r_p) + 0.25*(1-tau_as_p)));
	return(tau_atm);

}

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void longwave_radiation(short state, double pvap, double T, double fcloud, double *eps, double *eps_max, double *eps_min){

	double eep;

	if(state==0){
		*eps=1.5689*pow(pvap/(T+tk),0.176)*(1.0+0.20*pow(fcloud,2.0));	//fitting of Zongo data
		*eps_min=1.5689*pow(pvap/(T+tk),0.176);
		*eps_max=*eps_min*(1.0+0.2);

	}else if(state==1){
		*eps=1.24*pow((pvap/(T+tk)),0.14285714)*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0); //Brutsaert, 1975 + Pirazzini
		*eps_min=1.24*pow((pvap/(T+tk)),0.14285714);
		*eps_max=0.979;

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
void Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, double v, double *rh, double *rv, DOUBLEVECTOR *w){

	double z0t;
	double CHn, f, Rib, bh, c1h, c2h, c3h, c4h, Chx, ch, FH;

	//check
	if(zmt-zmu>0.5 || zmt-zmu<-0.5) t_error("If you use Louis' scheme, wind and temperature should be measured approximately at the same elevation on the ground");

	//roughness
	if(z0_z0t==0.0){	//rigid surface
		z0t=z0/8.5;
	}else{				//bending surface
		z0t=z0/z0_z0t;
	}

	/* CH neutrale [m/s]*/
	CHn=ka*ka/(log((zmu-d0)/z0)*log((zmt-d0)/z0t));

	/* calcola la funzione di stabilita' secondo la trattazione semplificata di Garrat, 1992 */
	f=(zmu/zmt-z0/zmt)/pow(zmt-z0t,0.5);

	/*calculation of the Richardson'number of Bulk */
	Rib=9.81*zmt*(Ta-Ts)/(0.5*(Ta+Ts)*v*v)*f*f;

	bh=23.0;
	if (Rib < 0.0) {
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
			c1h= -0.5128;
			c2h= -0.9448;
			c3h= 0.0643;
			c4h= 10.8925;
		}else{
			c1h= -0.3169;
			c2h= -0.3803;
			c3h= 0.0205;
			c4h= 7.5213;
		}
	}

	Chx=c1h*log(zmu/z0)+c2h*log(z0/z0t)+c3h*log(zmt/z0t)+c4h;

	ch=Chx*CHn*bh*f*pow(pow(zmt/z0t,1.0/3.0)-1.0,1.5);

	if (Rib < 0.0) {
		FH=1.0-bh*Rib/(1.0+ch*pow((-Rib),0.5));
	}else{
		FH=1.0/(pow(1.0+Chx*Rib,2.0));
	}

	/* calcola la resistenza aereodinamica */
	*rh=1/(CHn*FH*v);
	*rv=*rh;

	//check variables
	w->co[1]=Rib;
	w->co[2]=FH;
}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
double cz(double zmeas, double z0, double d0, double L, double (* unstab)(double z), double (* stab)(double z)){

	double c,zeta;

	zeta=(zmeas-d0)/L;

	if(zeta<0){
		c=log((zmeas-d0)/z0) - (*unstab)(zeta) + (*unstab)(z0/L);
	}else{
		c=log((zmeas-d0)/z0) + (*stab)(zeta) - (*stab)(z0/L);
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
void Businger(short a, long r, long c, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rh, double *rv, DOUBLEVECTOR *w){

	double L, L0, cm, u_star, ch, T_star, cv, Q_star;
	double z0v, z0t, z0q;
	long cont;
	FILE *f;

	//first guess of Obukhov length
	if(DT<0){
		L=10.0;
	}else{
		L=-10.0;
	}

	cont=0;

	do{

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

		//printf(".a...%f %f ....\n",ch,cv);
		//stop_execution();

		//Obukhov length
		L=-u_star*u_star*T/(ka*g*(T_star+0.61*Q_star*T));
		if(L*L0<0) L=L0;
		cont+=1;

	}while(fabs(L-L0)>0.001 && cont<=100 && a<4);

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
		fprintf(f,"ERROR: Displacement height greater than measurement elevations in cell [%ld,%ld]\n",r,c);
		fclose(f);
	}
	if(zmu<=z0 || zmt<=z0t || zmt<=z0q){
		f=fopen(files->co[ferr]+1,"a");
		fprintf(f,"ERROR: Elevation of sensors lower than roughness length in cell [%ld,%ld]: zmu=%f zmt=%f z0=%f z0t=%f\n",r,c,zmu,zmt,z0,z0t);
		fclose(f);
	}

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
	double horiz_H,// horizon elevation at a defined solar azimuth
		   tol; // toleraance on alpha to get the sun (if alpha>horiz_H+tol => sun )
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
	if (horiz_H>0) {
		tol=tol_mount;// at that particular azimuth, there is a mountain on the horizon, i.e. tol_mount has to be used
	} else
		tol=tol_flat;

	if (alpha>=horiz_H+tol) { /* sun higher than horizon => shadow=FALSE */
		shad=0; // SUN=TRUE, shadow=FALSE
	}else{ /* sun lower than horizon => shadow=TRUE */
		shad=1; // SUN=FALSE
	}
	return (shad);
}


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void sun(double hour, double JD, double *alpha, double *direction, double *E0, double latitude, double longitude, double standard_time){

	double Gamma, Delta, long_standard, Et, h;

	Gamma=2.0*Pi*JD/365.;

	/* latitudine e longitudine stazione met */
	long_standard=standard_time*Pi/12.0;

	/* correzione distanza Terra-Sole */
	*E0=1.00011+0.034221*cos(Gamma)+0.00128*sin(Gamma)+0.000719*cos(2*Gamma)+0.000077*sin(2*Gamma);
	/* correggo l'ora data in UTC ed ottengo l'ora reale locale*/
	Et=0.000075 + 0.001868*cos(Gamma) - 0.032077*sin(Gamma) - 0.014615*cos(2*Gamma) - 0.04089*sin(2*Gamma);	/*Correction for sideral day (rad)*/
	// Solar Declination (Declinazione solare)
	Delta=0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.006758*cos(2*Gamma)+0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.00148*sin(3*Gamma);

	h=hour + (longitude-long_standard)/omega + Et/omega;	/*Iqbal: formula 1.4.2*/
	if(h>=24) h-=24.0;
	if(h<0) h+=24.0;

	/**========= CALCOLO IL MOTO DEL SOLE ===============*/
	/*!!non funziona per latitudine=90gradi e -90gradi!!*/

	/* alpha: altezza solare in radianti */
	*alpha=asin( sin(latitude)*sin(Delta)+cos(latitude)*cos(Delta)*cos(omega*(12-h)) );

	/* direction: azimuth sole (da N, orario) in radianti */
	if(h<=12){
		if(*alpha==Pi/2.0){	/*sole allo zenit*/
			*direction=Pi/2.0;
		}else{
			*direction=Pi - acos((sin(*alpha)*sin(latitude)-sin(Delta))/(cos(*alpha)*cos(latitude)));
		}
	}else{
		if(*alpha==Pi/2.0){ /*sole allo zenit*/
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

