
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


//Authors: Stefano Endrizzi
//Date: december 2010
//Contents: Meteorological subroutines (included turbulent transfer)

#include "constant.h"
#include "keywords_file.h"
#include "struct.geotop.h"
#include "meteo.h"
#include "micromet.h"
#include "tabs.h"
#include "times.h"
#include "snow.h"
#include "geo_statistic.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void meteo_distr(METEO *met, ENERGY *egy, WATER *wat, TOPO *top, SNOW *snow, double time, PAR *par){

	long i,r,c;
	double t_station;
	
	static long line_interp;
	if(time==0) line_interp=0;

	//INTERPOLATION OF METEO VARIABLES
	for(i=1;i<=met->st->Z->nh;i++){
		time_conversion(par->JD0, par->year0, time+par->Dt, met->st->JD0->co[i], met->st->Y0->co[i], &t_station);
		t_station+=(met->st->ST->co[i]-par->ST)*3600.0;
		meteo_interp(met->data[i-1], met->st->Dt->co[i], t_station, met->var[i-1]);
	}

	//LAPSE RATES
	if(par->LRflag==1) meteo_interp2(&line_interp, met->LRv, met->LRp->co, 3, 1, 0, met->LRs, time+0.5*par->Dt, par->JD0, par->year0, NoV);
	if(met->LRv[met->LRp->co[1]]==NoV) met->LRv[met->LRp->co[1]]=-6.5;	//Tair lapse rate
	if(met->LRv[met->LRp->co[2]]==NoV) met->LRv[met->LRp->co[2]]=-2.0;	//Tdew lapse rate
	if(met->LRv[met->LRp->co[3]]==NoV) met->LRv[met->LRp->co[3]]=0.0;	//Prec lapse rate
	//printf("%ld %f %ld %f %ld %f\n",met->LRp->co[1],met->LRv[met->LRp->co[1]],met->LRp->co[2],met->LRv[met->LRp->co[2]],met->LRp->co[3],met->LRv[met->LRp->co[3]]);

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
			par->iobsint, iT, iRh, iWs, iWd, iPt, met->Tgrid, met->RHgrid, met->Vgrid, met->Vdir, met->Pgrid, wat->PrecTot, met->LRv[2],
			met->LRv[3], met->LRv[4]);
		
	}else{

		met->V=5.0;
		met->RH=0.7;
		if(met->column[0][iWs]!=-1) met->V=Fmax(par->Vmin,met->var[0][met->column[0][iWs]]);
		if(met->column[0][iRh]!=-1) met->RH=Fmax(par->RHmin/100.,met->var[0][met->column[0][iRh]]/100.);

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				
				if(top->Z0->co[r][c]!=UV->V->co[2]){

					//default values
					met->Tgrid->co[r][c]=0.0;
					met->Pgrid->co[r][c]=pressure(top->Z0->co[r][c], Pa0, 0.0);
					wat->PrecTot->co[r][c]=0.0;

					//constant values
					if(met->column[0][iT]!=-1) met->Tgrid->co[r][c]=met->var[0][met->column[0][iT]];
					if(met->column[0][iPs]!=-1) met->Pgrid->co[r][c]=met->var[0][met->column[0][iPs]];
					if(met->column[0][iPt]!=-1) wat->PrecTot->co[r][c]=met->var[0][met->column[0][iPt]];

				}
			}
		}
	}

	if(par->en_balance==0){
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				wat->Pnet->co[r][c]=wat->PrecTot->co[r][c]/3600.0;	//from [mm/h] to [mm/s]
			}
		}
	}
	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

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

double pressure(double Dz, double P0, double gamma){

	double P;
	P=P0*exp(-Dz*0.00013);
	return(P);
}

double temperature(double Dz, double T0, double gamma){

	double T;
	T=(T0+tk)*exp(-(gamma/(T0+tk))*Dz)-tk;
	return(T);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void meteo_interp(double **data, double Dt, double t, double *out)

{

	long i,j,n;

	//printf("t:%f Dt:%f\n",t,Dt);

	i=floor(t/Dt);	//previous instant
	n=dim2(data);	//number of time steps in the data matrix

	//printf("i:%ld %ld\n",i,dim2(data));

	if(i<0){
		//t_error("ERROR 1 in the met data!!");
		for(j=0;j<dim1(out);j++){
			out[j]=NoV;
		}

	}else if(i>n-1){
		//t_error("ERROR 2 in the met data!!");
		for(j=0;j<dim1(out);j++){
			out[j]=NoV;
		}

	}else if(i==n-1){
		for(j=0;j<dim1(out);j++){
			out[j]=data[i-1+1][j];
			//printf("1: i:%ld j:%ld data:%f out:%f\n",i,j,data[i-1+1][j],out[j]);
		}

	}else{
		for(j=0;j<dim1(out);j++){
			if(data[i-1+2][j]!=NoV && data[i-1+1][j]!=NoV){
				out[j]=data[i-1+1][j]+(data[i-1+2][j]-data[i-1+1][j])*(t-i*Dt)/Dt;
				//printf("2: i:%ld j:%ld data:%f data:%f out:%f\n",i,j,data[i-1+2][j],data[i-1+1][j],out[j]);
			}else{
				out[j]=NoV;
			}
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void meteo_interp2(long *ibeg, double *out, long *col_data, long n_data, long col_JD, long col_y, double **data, double t, double JD0, long y0, double novalue){

	long i,j,col;
	double t0, t10, t11;
	
	//t0 time in [s] from the origin of the data in the file
	time_conversion(JD0, y0, t, data[0][col_JD], (long)data[0][col_y], &t0);
	//printf("JD0:%f y0:%ld t:%f %f %ld, %f\n",JD0, y0, t, data[0][col_JD], (long)data[0][col_y], t0);
	//stop_execution();

	if(t0<0) t_error("Error 1 in subroutine interp2");

	i=(*ibeg);
	
	do{
		if(i-1>=dim2(data)) t_error("Error 2 in subroutine interp2");
		get_time(&t10, data[i][col_JD], (long)data[i][col_y], data[0][col_JD], (long)data[0][col_y]);
		//printf("t10:%f\n",t10);
		get_time(&t11, data[i+1][col_JD], (long)data[i+1][col_y], data[0][col_JD], (long)data[0][col_y]);
		//printf("t10:%f\n",t11);
		//printf("i:%ld\n",i);
		i++;
	}while(t0<t10 || t0>t11);

	for(j=1;j<=n_data;j++){
		col=col_data[j];
		if(data[i][col]==novalue || data[i-1][col]==novalue){
			out[col]=novalue;
		}else{
			out[col] = (data[i][col]*(t0-t10)+data[i-1][col]*(t11-t0))/(t11-t10);
		}
	}
	
	*ibeg=i-1;
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

void sat_vap_pressure(double *p, double *dp_dT, double T, double P){	//water vapour pressure p(T,P)
//pressure in [mbar] - p vapour pressure - P atmospheric pressure, temperature in [C]
	double A, b, c;
	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;
	*p=A*exp(b*T/(c+T));
	*dp_dT=*p*(b/(c+T)-b*T/pow(c+T,2.0));
}

void sat_vap_pressure_2(double *p, double T, double P){	//water vapour pressure p(T,P)
	//pressure in [mbar] - p vapour pressure - P atmospheric pressure, temperature in [C]
	double A, b, c;
	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;
	*p=A*exp(b*T/(c+T));
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

void sat_spec_humidity_2(double *Q, double RH, double T, double P){
	double p;
	sat_vap_pressure_2(&p, T, P);
	*Q=RH*spec_humidity(p, P);
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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


