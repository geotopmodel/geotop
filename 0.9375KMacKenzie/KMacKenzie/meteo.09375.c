
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



//Authors: Stefano Endrizzi
//Date: december 2008
//Contents: Meteorological subroutines (included turbulent transfer)
#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "geo_statistic.09375.h"
#include "times.h"
#include "snow.09375.h"
#include "meteo.09375.h"
#include "micromet.h"
#include "shadows.h"
#include "tabs.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern char *error_file_name;
extern long Nl, Nr, Nc;
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void meteo_distr(METEO *met, ENERGY *egy, WATER *wat, TOPO *top, SNOW *snow, double time, PAR *par){
	/* Author:  Stefano Endrizzi  Year:
		 * function that ...
		 * Input:	met:		current time [second] based on a new origin given by the beginning of the first datum of the meteo station
		 * 			egy:		time step of the meteo station
		 * 			wat:	matrix of the meteo data
		 * 			top:
		 * 			snow:
		 * 			time:
		 * 			par:
		 * comment: Matteo Dall'Amico, October 2009 */
	long i,r,c;
	double t_station;


	//INTERPOLATION OF METEO VARIABLES
	for(i=1;i<=met->st->Z->nh;i++){
		time_conversion(par->JD0, par->year0, time+par->Dt, met->st->JD0->co[i], met->st->Y0->co[i], &t_station);
		t_station+=(met->st->ST->co[i]-par->ST)*3600.0;
		meteo_interp(met->data[i-1], met->st->Dt->co[i], t_station, met->var[i-1]);
	}

	//LAPSE RATES
	if(par->LRflag==1) meteo_interp2(met->LRv, met->LRp, 1, 0, met->LRs, time+0.5*par->Dt, par);
	if(met->LRv[2]==NoV) met->LRv[2]=-6.5;	//Tair lapse rate
	if(met->LRv[3]==NoV) met->LRv[3]=-2.0;	//Tdew lapse rate
	if(met->LRv[4]==NoV) met->LRv[4]=0.0;	//Prec lapse rate

	//DISTRIBUTION OF METEROLOGICAL VARIABLES FROM MEASUREMENTS IN SOME STATIONS
	if(par->micromet==1){
		if(par->topoflag==1){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(top->Z0->co[r][c]!=UV->V->co[2]){
						top->Zm->co[r][c] = top->Z0->co[r][c] + DEPTH(r, c, snow->lnum, snow->Dzl);
						topo_data(UV->U->co[1], UV->U->co[2], top->Zm, top->curv_m, top->slope_m, top->slopeaz_m, par->curve_len_scale, NoV);
					}
				}
			}
		}

		Micromet(UV, top->Zm, top->curv_m, top->slope_m, top->slopeaz_m, met, par->slopewt, par->curvewt, par->Vmin, par->dn, par->ifill,
			par->iobsint, iT, iRh, iWs, iWd, iPt, met->Tgrid, met->RHgrid, met->Vgrid, met->Vdir, met->Pgrid, wat->total, met->LRv[2],
			met->LRv[3], met->LRv[4]);


	}else{

		met->V=5.0;
		met->RH=0.7;
		if(met->column[0][iWs]!=-1) met->V=Fmax(par->Vmin,met->var[0][met->column[0][iWs]]);
		if(met->column[0][iRh]!=-1) met->RH=Fmax(par->RHmin/100.,met->var[0][met->column[0][iRh]]/100.);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){

				//default values
				met->Tgrid->co[r][c]=0.0;
				met->Pgrid->co[r][c]=1.E3;
				wat->total->co[r][c]=0.0;

				//constant values
				if(met->column[0][iT]!=-1) met->Tgrid->co[r][c]=met->var[0][met->column[0][iT]];
				if(met->column[0][iPs]!=-1) met->Pgrid->co[r][c]=met->var[0][met->column[0][iPs]];
				//printf("\nmet->column[0][iWs]=%ld, met->column[0][iRh]=%ld,met->column[0][iT]=%ld, met->column[0][iPs]=%ld",met->column[0][iWs],met->column[0][iRh],met->column[0][iT],met->column[0][iPs]);stop_execution();
				//if(met->column[1][iPt]!=-1) wat->total->co[r][c]=met->var[1][met->column[1][iPt]];// was like this, Matteo 7/10/09
				if(met->column[0][iPt]!=-1) wat->total->co[r][c]=met->var[0][met->column[0][iPt]];// ho messo 0 al posto di 1, Matteo 7/10/09
			}
		}

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
	if(met->column[i-1][iT]==-1){
		printf("WARNING: in met file no data of temperature, set at default values\n");
	}else{
		vert_distr(met->Tgrid, Z, met->st->Z->co[1], met->var[i-1][met->column[i-1][iT]], met->LRv[2], (*temperature));
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

	//check if rain < 0 (nodata) and in that cae I do kriging another time
	for(i=0;i<n;i++){
		station_nodata[i]=0;
		if(data[i]<=novalue) station_nodata[i]=1;
		index_nodata+=station_nodata[i];
	}

	initialize_doublematrix(out, 0.0);

	if(index_nodata==n){
		f=fopen(error_file_name,"a");
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
					// number of rows (h*k)
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
	n=dim2(data);	//number of time rows in the data matrix

	if(i<0){
		//t_error("ERROR 1 in the met data!!");
		for(j=0;j<dim1(data[i]);j++){
			out[j]=NoV;
		}

	}else if(i>n-1){
		//t_error("ERROR 2 in the met data!!");
		for(j=0;j<dim1(data[i]);j++){
			out[j]=NoV;
		}

	}else if(i==n-1){// we are at the end of the data matrix => the last value is given
		for(j=0;j<dim1(data[i]);j++){
			out[j]=data[i-1+1][j];
			//printf("n-1: i:%ld j:%ld data:%f out:%f\n",i,j,data[i-1+1][j],out[j]);
		}

	}else{
		for(j=0; j<dim1(data[i]); j++){
			out[j]=data[i-1+1][j]+(data[i-1+2][j]-data[i-1+1][j])*(t-i*Dt)/Dt;
			//printf("  1: i:%ld j:%ld data:%f data:%f out:%f\n",i,j,data[i-1+2][j],data[i-1+1][j],out[j]);
		}
	}

}


void meteo_interp2(double *out, LONGVECTOR *col_data, long col_JD, long col_y, double **data, double t, PAR *par){
/* Matteo: va bene per i Lapse rates. a differenza di meteo_interp e' piu' evoluta perche' non ha bisogno del file meteo.txt che ti
 * dice qual e' il delta_t e l'istante di partenza, ma accetta dati irregolari
 * col_data le colonne dei dati (JD e anno):
 * out: vettore di output */
	long i,j,col;
	double t0, t10, t11;

	//t0 time in [s] from the origin of the data in the file
	time_conversion(par->JD0, par->year0, t, data[0][col_JD], (long)data[0][col_y], &t0);

	if(t0<0) t_error("Error 1 in subroutine interp2");

	i=0;
	do{
		if(i-1>=dim2(data)) t_error("Error 2 in subroutine interp2");
		get_time(&t10, data[i][col_JD], (long)data[i][col_y], data[0][col_JD], (long)data[0][col_y]);
		get_time(&t11, data[i+1][col_JD], (long)data[i+1][col_y], data[0][col_JD], (long)data[0][col_y]);
		i++;
	}while(t0<t10 || t0>t11);

	for(j=1;j<=col_data->nh;j++){
		col=col_data->co[j];
		out[col] = (data[i][col]*(t0-t10)+data[i-1][col]*(t11-t0))/(t11-t10);
	}
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

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


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/




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
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/





/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/



/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

/*==================================================================================================================*/
/*==================================================================================================================*/

/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/


/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/




/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/
/*==================================================================================================================*/

