
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
//Date: 13 November 2005
//Contents: Snow subroutines
#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "snow.09375.h"
#include "write_dem.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;

/*==================================================================================================================*/
/*==================================================================================================================*/
/*1. subroutine RHO_NEWFALLENSNOW*/
/*==================================================================================================================*/
/*==================================================================================================================*/

double rho_newlyfallensnow(double u, double Tatm, double Tfreez)
	/* Author:  Stefano Endrizzi  Year:
	* function that calculates the density of the newly fallen snow based on Jordan et al, (1999)
	* Input:	u: wind velocity [m/s]
	* 			Tatm: Air temperature [ûC]
	* 			Tfreez: Freezing temperature
	* Output: snow density [Kg/m3]
	* comment: Matteo Dall'Amico, May 2009 */
{

double rho,T;

if(Tatm>260.15-tk){
	T=Tatm;
	if(T>=275.65-tk) T=275.65-tk;
	rho=500.0*(1.0-0.951*exp(-1.4*pow(278.15-tk-T,-1.15)-0.008*pow(u,1.7)));
}else{
	rho=500.0*(1.0-0.904*exp(-0.008*pow(u,1.7)));
}

/*if(Tatm>Tfreez+2.0){
	rho=50.0+1.7*pow(17,1.5);
}else if(Tatm<=Tfreez+2.0 && Tatm>Tfreez-15.0){
	rho=50.0+1.7*pow(Tatm-Tfreez+15.0,1.5);
}else{
	rho=50.0;
}*/

return(rho);

}



double kidrmax_snow(double rho)

{

double d,ks,K;

if(rho<50){
	d=2.0E-4;
}else if(rho>=50 && rho<100){
	d=2.0E-4*(100-rho)+4.0E-4*(rho-50);
	d/=50.0;
}else if(rho>=100 && rho<200){
	d=4.0E-4*(200-rho)+8.0E-4*(rho-100);
	d/=100.0;
}else if(rho>=200 && rho<300){
	d=8.0E-4*(300-rho)+1.4E-3*(rho-200);
	d/=100.0;
}else if(rho>=300 && rho<400){
	d=1.4E-3*(400-rho)+2.4E-3*(rho-300);
	d/=100.0;
}else if(rho>=400 && rho<500){
	d=2.4E-3*(500-rho)+3.5E-3*(rho-400);
	d/=100.0;
}else if(rho>=500 && rho<600){
	d=3.5E-3*(600-rho)+4.6E-3*(rho-500);
	d/=100.0;
}else if(rho>=600 && rho<700){
	d=4.6E-3*(700-rho)+4.6E-3*(rho-600);
	d/=100;
}else if(rho>=700){
	d=4.6E-3;
}

//Shimizu (1970)
ks=0.077*pow(d,2.0)*exp(-7.8*rho/rho_w);
K=ks*rho_w*g/mu_l;

return(K);

}




/*==================================================================================================================*/
/*==================================================================================================================*/
/*2. subroutine SNOW_COMPACTATION*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void snow_compactation(long r, long c, long l, SNOW *snow, double slope, PAR *par, double *CR1, double *CR2)

{

long m;
double theta_i,rho,c1,c2,c3,c4,c5,load,eta0,eta;

theta_i=snow->w_ice->co[l][r][c]/(0.001*snow->Dzl->co[l][r][c]*rho_i);	
rho=(snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c])/(0.001*snow->Dzl->co[l][r][c]);

if( theta_i < par->snow_maxpor ){

	//DESTRUCTIVE METAMORPHISM
	//gamma=theta*rho
	if(theta_i*rho_i<=par->snow_density_cutoff){
		c1=par->drysnowdef_rate;
	}else{
		c1=par->drysnowdef_rate*exp(-0.046*(rho_i*theta_i - par->snow_density_cutoff));
	}
	if(snow->w_liq->co[l][r][c]>0) c1*=par->wetsnowdef_rate;
	c2=2.777E-6; //[s^-1]
	c3=0.04;	 //[K^-1]
	*CR1=-c1*c2*exp(-c3*(Tfreezing-snow->T->co[l][r][c]));

	//OVERBURDEN
	eta0=par->snow_viscosity;  //[kg s m^-2]
	c4=0.08;	 //[K^-1]
	c5=0.023;    //[m^3 kg^-1]
	load=0.0;    //[kg m^-2]
	for(m=l;m<=snow->lnum->co[r][c];m++){
		load+=(snow->w_ice->co[m][r][c]+snow->w_liq->co[m][r][c]);
	}
	load*=fabs(cos(slope));
	eta=eta0*exp(c4*(Tfreezing-snow->T->co[l][r][c])+c5*(rho_i*theta_i));
	*CR2=-load/eta;

	//printf("l:%ld %e %e eta0:%f T:%f thi:%f Dzl:%f\n",l,*CR1,*CR2,eta0,snow->T->co[l][r][c],theta_i,snow->Dzl->co[l][r][c]);

	snow->Dzl->co[l][r][c]*=(1.0 + ((*CR1)+(*CR2))*par->Dt);

}else{

	*CR1=0.0;
	*CR2=0.0;

}

}



/*==================================================================================================================*/
/*==================================================================================================================*/
/*3. subroutine SNOW_LAYER_COMBINATION*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void snow_layer_combination(long r, long c, SNOW *snow, double Ta, long linf, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, double time)
	/* Author:  Stefano Endrizzi  Year:
	* function that finds the equilibrium number of snow layers depending on the given layer depth
	* and on the allowed minimum and maximum depth.
	* If snow[l]<Dmin[n-l+1] => a layer has less snow than allowed, therefore its mass and energy
	* has to be shared with the lower one.
	* If snow[l]>Dmax[n-l+1] => a layer has more snow than allowed, therefore a new snow layer has
	* to be initialized.
	* The scheme is symmetrical
	* The routine goes on until an equilibrium is found.
	* Input:	r: row of the pixel in the raster
	* 			c: column of the pixel in the raster
	* 			Dmin: minimum allowed snow depth
	* 			Dmax: maximum allowed snow depth
	* 			snow: snow structure
	* 			Ta: not used
	* 			linf: layer of unlimited thickness (beginning from below)
	* 			time: not used
	* Output: structure "snow" modified with correct number of snow layers
	* comment: Matteo Dall'Amico, May 2009 */
{

long l, cont=0;
long max=Dmin->nh;// max number of snow layers
short occuring;
double D0=0.0,D1=0.0;
double D=0.0;//D=snow depth(mm)
DOUBLEVECTOR *Dmin2,*Dmax2,*wice,*wliq,*Dz,*Temp;

for(l=1;l<=max;l++){
	D+=snow->Dzl->co[l][r][c];
}

//1. If the snow depth is too small, then it is reset to 0
if(D<=0.000001){

	snow->lnum->co[r][c]=0;
	snow->type->co[r][c]=0;
	for(l=1;l<=max;l++){
		snow->Dzl->co[l][r][c]=0.0;
		snow->w_ice->co[l][r][c]=0.0;
		snow->w_liq->co[l][r][c]=0.0;
		snow->T->co[l][r][c]=-99.0;	//Temperatura di inizializzazione
	}


//2. If D<Dmin(1), we are in the simplified case
}else if(D<Dmin->co[1] && snow->lnum->co[r][c]>0){

	snow->type->co[r][c]=1;
	if(snow->lnum->co[r][c]>1){

		for(l=snow->lnum->co[r][c];l>1;l--){
			snowlayer_merging(r, c, snow, l, l-1, l-1);
		}

		for(l=2;l<=snow->lnum->co[r][c];l++){
			snow->T->co[l][r][c]=-99.0;
			snow->Dzl->co[l][r][c]=0.0;
			snow->w_liq->co[l][r][c]=0.0;
			snow->w_ice->co[l][r][c]=0.0;
		}

		snow->lnum->co[r][c]=1;

	}
}


//3. se non c'e' ancora nessun layer di neve e 1mm<=z<Dmin(1), si inizializza un nuovo layer di neve (caso semplificato)
if(snow->lnum->co[r][c]==0 && snow->Dzl->co[1][r][c]>=0.000001 && snow->Dzl->co[1][r][c]<Dmin->co[1]){

	snow->lnum->co[r][c]=1;
	snow->type->co[r][c]=1;
	//si assegna alla neve la temperatura dell'aria
	if(Ta<Tfreezing){
		snow->T->co[1][r][c]=Ta;
	}else{
		snow->T->co[1][r][c]=Tfreezing;
	}
}


//4. se non c'e' ancora nessun layer di neve e z>=Dmin(1), si inizializza un nuovo layer di neve (caso semplificato)
if(snow->lnum->co[r][c]==0 && snow->Dzl->co[1][r][c]>=Dmin->co[1]){

	snow->lnum->co[r][c]=1;
	snow->type->co[r][c]=2;
	//si assegna alla neve la temperatura dell'aria
	if(Ta<Tfreezing){
		snow->T->co[1][r][c]=Ta;
	}else{
		snow->T->co[1][r][c]=Tfreezing;
	}
}

//5. se siamo nel caso semplificato e z>=Dmin(1), si entra nel caso ordinario
if(snow->type->co[r][c]==1 && snow->Dzl->co[1][r][c]>=Dmin->co[1]){

	snow->type->co[r][c]=2;

	//si assegna alla neve la temperatura dell'aria
	/*if(Ta<Tfreezing){
		snow->T->co[1][r][c]=Ta;
	}else{
		snow->T->co[1][r][c]=Tfreezing;
	}*/
}


//6. caso ordinario

// SIMMETRICAL PARAMETERIZATION SCHEME (new)

// if you want to use the old one, COMMENT FROM HERE

Dmax->co[linf]=1.0E10;	//snow layer of unlimited thickness

if(snow->type->co[r][c]==2){

	wliq=new_doublevector(max);
	wice=new_doublevector(max);
	Temp=new_doublevector(max);
	Dz=new_doublevector(max);

	Dmin2=new_doublevector(max);
	Dmax2=new_doublevector(max);
	min_max_layer(snow->lnum->co[r][c], Dmin, Dmax, Dmin2, Dmax2, linf);

	for(l=1;l<=max;l++){
		D0+=snow->Dzl->co[l][r][c];
	}

	do{

		//show_Dminmax(r, c, Dmin2->co, Dmax2->co, snow->lnum->co[r][c]);
		//write_snow_all(r, c, snow);

		cont+=1;

		if(cont>10){
			show_Dminmax(r, c, Dmin2->co, Dmax2->co, snow->lnum->co[r][c]);
			write_snow_all(r, c, snow);
			printf("Iteration to assign new thicknesses to the snow layers does not converge, r:%ld c:%ld\n",r,c);
			t_error("Error 1 in snow combination - Change values of blocks 5 and 6 parameter file");
		}

		occuring=0;

		//assign
		for(l=1;l<=max;l++){
			wice->co[l]=snow->w_ice->co[l][r][c];
			wliq->co[l]=snow->w_liq->co[l][r][c];
			Dz->co[l]=snow->Dzl->co[l][r][c];
			Temp->co[l]=snow->T->co[l][r][c];
		}

		//trying to adjust the layer thicknesses maintaining the same layer number
		for(l=snow->lnum->co[r][c];l>=linf;l--){
			if(Dz->co[l]>Dmax2->co[l] && l>1) set_snow(r, c, wliq->co, wice->co, Temp->co, Dz->co, l, l-1, Dmax2->co[l]);
			if(Dz->co[l]<Dmin2->co[l] && Dz->co[l-1]>=Dmin2->co[l]-Dz->co[l] && l>1) set_snow(r, c, wliq->co, wice->co, Temp->co, Dz->co, l, l-1, Dmin2->co[l]);
		}

		for(l=1;l<=linf;l++){
			if(Dz->co[l]>Dmax2->co[l] && l<snow->lnum->co[r][c]) set_snow(r, c, wliq->co, wice->co, Temp->co, Dz->co, l, l+1, Dmax2->co[l]);
			if(Dz->co[l]<Dmin2->co[l] && Dz->co[l+1]>=Dmin2->co[l]-Dz->co[l] && l<snow->lnum->co[r][c]) set_snow(r, c, wliq->co, wice->co, Temp->co, Dz->co, l, l+1, Dmin2->co[l]);
		}

		//checking if it is ok now
		for(l=1;l<=snow->lnum->co[r][c];l++){

			if(Dz->co[l]<Dmin2->co[l] || Dz->co[l]>Dmax2->co[l]) occuring=1;
		}

		//if it is ok, it is done; otherwise try to split or merge layers
		if(occuring==0){

			for(l=1;l<=snow->lnum->co[r][c];l++){
				snow->w_ice->co[l][r][c]=wice->co[l];
				snow->w_liq->co[l][r][c]=wliq->co[l];
				snow->Dzl->co[l][r][c]=Dz->co[l];
				snow->T->co[l][r][c]=Temp->co[l];
			}

		}else{
		
			for(l=1;l<=snow->lnum->co[r][c];l++){	
			
				//write_snow_all(r,c,snow);
				
				if(snow->Dzl->co[l][r][c]<Dmin2->co[l]){
					merge_layers(r, c, snow, l);	
					min_max_layer(snow->lnum->co[r][c], Dmin, Dmax, Dmin2, Dmax2, linf);	
				}
												
				if(snow->lnum->co[r][c]<max && snow->Dzl->co[l][r][c]>Dmax2->co[l]){
					split_layers(r, c, snow, l);
					min_max_layer(snow->lnum->co[r][c], Dmin, Dmax, Dmin2, Dmax2, linf);
				}
			}
		}

	}while(occuring==1);

	for(l=1;l<=max;l++){
		D1+=snow->Dzl->co[l][r][c];
	}
	if(fabs(D0-D1)>0.001){
		printf("r:%ld c:%ld Dold:%f Dnew:%f\n",r,c,D0,D1);
		t_error("Error 2 in snow combination - Change values of blocks 5 and 6 parameter file");
	}

	free_doublevector(wice);
	free_doublevector(wliq);
	free_doublevector(Temp);
	free_doublevector(Dz);

	free_doublevector(Dmin2);
	free_doublevector(Dmax2);
}

// TO HERE (comment if you want the old snow scheme)






// SNTHERM (CLASSICAL) SNOW PARAMETERIZATION SCHEME (new)

// if you want to use the new one, COMMENT FROM HERE
/*
if( Dmax->nh<max || Dmin->nh<max) t_error("Snow Dmin and Dmax given are not enough. Provide them for all the layers");
Dmax->co[max]=1.0E10;
long n, m;
double h;
if(snow->type->co[r][c]==2){
	n=snow->lnum->co[r][c];
	do{
		occuring=0;
		for(l=n;l>=1;l--){	//dall'alto in basso
			if(n==max && l==1){
				occuring=0;
			}else if(snow->Dzl->co[l][r][c]>Dmax->co[n-l+1]){
				if(n<max){//posso tranquillamente aggiungere un layer in piu'
					n=n+1;
					for(m=n;m>=l+2;m--){
						snow->Dzl->co[m][r][c]=snow->Dzl->co[m-1][r][c];
						snow->w_ice->co[m][r][c]=snow->w_ice->co[m-1][r][c];
						snow->w_liq->co[m][r][c]=snow->w_liq->co[m-1][r][c];
						snow->T->co[m][r][c]=snow->T->co[m-1][r][c];
					}
					snow->Dzl->co[l][r][c]*=0.5;
					snow->w_ice->co[l][r][c]*=0.5;
					snow->w_liq->co[l][r][c]*=0.5;
					snow->Dzl->co[l+1][r][c]=snow->Dzl->co[l][r][c];
					snow->w_ice->co[l+1][r][c]=snow->w_ice->co[l][r][c];
					snow->w_liq->co[l+1][r][c]=snow->w_liq->co[l][r][c];
					snow->T->co[l+1][r][c]=snow->T->co[l][r][c];
				}else{//devo trovare il posto per un layer in piu'
					if(n!=max) t_error("Error 1 in snow combination");	//controllo
					if(l>=3){	//combino i 2 layer piu' bassi
						h=(c_ice*snow->w_ice->co[1][r][c] + c_liq*snow->w_liq->co[1][r][c])*(snow->T->co[1][r][c]-Tfreezing) +
								Lf*snow->w_liq->co[1][r][c];
						h+=(c_ice*snow->w_ice->co[2][r][c] + c_liq*snow->w_liq->co[2][r][c])*(snow->T->co[2][r][c]-Tfreezing) +
								Lf*snow->w_liq->co[2][r][c];
						snow->Dzl->co[1][r][c]+=snow->Dzl->co[2][r][c];
						snow->w_ice->co[1][r][c]+=snow->w_ice->co[2][r][c];
						snow->w_liq->co[1][r][c]+=snow->w_liq->co[2][r][c];
						if(h<0.0){
							snow->T->co[1][r][c]=Tfreezing + h/(c_ice*snow->w_ice->co[1][r][c]+c_liq*snow->w_liq->co[1][r][c]);
						}else if(h>Lf*snow->w_liq->co[1][r][c]){
							snow->T->co[1][r][c]=Tfreezing + (h-Lf*snow->w_liq->co[1][r][c])/(c_ice*snow->w_ice->co[1][r][c]+
									c_liq*snow->w_liq->co[1][r][c]);
						}else{
							snow->T->co[1][r][c]=Tfreezing;
						}
						for(m=2;m<=n-1;m++){
							snow->Dzl->co[m][r][c]=snow->Dzl->co[m+1][r][c];
							snow->w_ice->co[m][r][c]=snow->w_ice->co[m+1][r][c];
							snow->w_liq->co[m][r][c]=snow->w_liq->co[m+1][r][c];
							snow->T->co[m][r][c]=snow->T->co[m+1][r][c];
						}
						for(m=n;m>=l+1;m--){
							snow->Dzl->co[m][r][c]=snow->Dzl->co[m-1][r][c];
							snow->w_ice->co[m][r][c]=snow->w_ice->co[m-1][r][c];
							snow->w_liq->co[m][r][c]=snow->w_liq->co[m-1][r][c];
							snow->T->co[m][r][c]=snow->T->co[m-1][r][c];
						}
						snow->Dzl->co[l-1][r][c]*=0.5;
						snow->w_ice->co[l-1][r][c]*=0.5;
						snow->w_liq->co[l-1][r][c]*=0.5;
						snow->Dzl->co[l][r][c]=snow->Dzl->co[l-1][r][c];
						snow->w_ice->co[l][r][c]=snow->w_ice->co[l-1][r][c];
						snow->w_liq->co[l][r][c]=snow->w_liq->co[l-1][r][c];
						snow->T->co[l][r][c]=snow->T->co[l-1][r][c];
					}else{
						if(l!=2) t_error("Error 2 in snow combination");	//controllo
						snow->Dzl->co[2][r][c]*=0.5;
						snow->w_ice->co[2][r][c]*=0.5;
						snow->w_liq->co[2][r][c]*=0.5;
						h=(c_ice*snow->w_ice->co[1][r][c] + c_liq*snow->w_liq->co[1][r][c])*(snow->T->co[1][r][c]-Tfreezing) +
								Lf*snow->w_liq->co[1][r][c];
						h+=(c_ice*snow->w_ice->co[2][r][c] + c_liq*snow->w_liq->co[2][r][c])*(snow->T->co[2][r][c]-Tfreezing) +
								Lf*snow->w_liq->co[2][r][c];
						snow->Dzl->co[1][r][c]+=snow->Dzl->co[2][r][c];
						snow->w_ice->co[1][r][c]+=snow->w_ice->co[2][r][c];
						snow->w_liq->co[1][r][c]+=snow->w_liq->co[2][r][c];
						if(h<0.0){
							snow->T->co[1][r][c]=Tfreezing + h/(c_ice*snow->w_ice->co[1][r][c]+c_liq*snow->w_liq->co[1][r][c]);
						}else if(h>Lf*snow->w_liq->co[1][r][c]){
							snow->T->co[1][r][c]=Tfreezing + (h-Lf*snow->w_liq->co[1][r][c])/(c_ice*snow->w_ice->co[1][r][c]+
								c_liq*snow->w_liq->co[1][r][c]);
						}else{
							snow->T->co[1][r][c]=Tfreezing;
						}
					}
				}
			}
		}

		for(l=1;l<=n;l++){	//dal basso in alto
			if(snow->Dzl->co[l][r][c]<Dmin->co[n-l+1]){
				if(l==1 && n==2){
					snow->Dzl->co[1][r][c]+=snow->Dzl->co[2][r][c];
					snow->w_ice->co[1][r][c]+=snow->w_ice->co[2][r][c];
					snow->w_liq->co[1][r][c]+=snow->w_liq->co[2][r][c];
					snow->T->co[1][r][c]=snow->T->co[2][r][c];
					n=n-1;
					snow->Dzl->co[2][r][c]=0.0;
					snow->w_ice->co[2][r][c]=0.0;
					snow->w_liq->co[2][r][c]=0.0;
					snow->T->co[2][r][c]=-99.0;
				}else if(l==1 && n>2){
					h=(c_ice*snow->w_ice->co[1][r][c] + c_liq*snow->w_liq->co[1][r][c])*(snow->T->co[1][r][c]-Tfreezing) +
						Lf*snow->w_liq->co[1][r][c];
					h+=(c_ice*snow->w_ice->co[2][r][c] + c_liq*snow->w_liq->co[2][r][c])*(snow->T->co[2][r][c]-Tfreezing) +
						Lf*snow->w_liq->co[2][r][c];
					snow->Dzl->co[1][r][c]+=snow->Dzl->co[2][r][c];
					snow->w_ice->co[1][r][c]+=snow->w_ice->co[2][r][c];
					snow->w_liq->co[1][r][c]+=snow->w_liq->co[2][r][c];
					if(h<0.0){
						snow->T->co[1][r][c]=Tfreezing + h/(c_ice*snow->w_ice->co[1][r][c]+c_liq*snow->w_liq->co[1][r][c]);
					}else if(h>Lf*snow->w_liq->co[1][r][c]){
						snow->T->co[1][r][c]=Tfreezing + (h-Lf*snow->w_liq->co[1][r][c])/(c_ice*snow->w_ice->co[1][r][c]+
							c_liq*snow->w_liq->co[1][r][c]);
					}else{
						snow->T->co[1][r][c]=Tfreezing;
					}
					n=n-1;
					for(m=2;m<=n;m++){
						snow->Dzl->co[m][r][c]=snow->Dzl->co[m+1][r][c];
						snow->w_ice->co[m][r][c]=snow->w_ice->co[m+1][r][c];
						snow->w_liq->co[m][r][c]=snow->w_liq->co[m+1][r][c];
						snow->T->co[m][r][c]=snow->T->co[m+1][r][c];
					}
					snow->Dzl->co[n+1][r][c]=0.0;
					snow->w_ice->co[n+1][r][c]=0.0;
					snow->w_liq->co[n+1][r][c]=0.0;
					snow->T->co[n+1][r][c]=-99.0;
				}else if(l!=n){
					h=(c_ice*snow->w_ice->co[l][r][c] + c_liq*snow->w_liq->co[l][r][c])*(snow->T->co[l][r][c]-Tfreezing) +
						Lf*snow->w_liq->co[l][r][c];
					h+=(c_ice*snow->w_ice->co[l-1][r][c] + c_liq*snow->w_liq->co[l-1][r][c])*(snow->T->co[l-1][r][c]-Tfreezing) +
						Lf*snow->w_liq->co[l-1][r][c];
					snow->Dzl->co[l-1][r][c]+=snow->Dzl->co[l][r][c];
					snow->w_ice->co[l-1][r][c]+=snow->w_ice->co[l][r][c];
					snow->w_liq->co[l-1][r][c]+=snow->w_liq->co[l][r][c];
					if(h<0.0){
						snow->T->co[l-1][r][c]=Tfreezing + h/(c_ice*snow->w_ice->co[l-1][r][c]+c_liq*snow->w_liq->co[l-1][r][c]);
					}else if(h>Lf*snow->w_liq->co[l-1][r][c]){
						snow->T->co[l-1][r][c]=Tfreezing + (h-Lf*snow->w_liq->co[l-1][r][c])/(c_ice*snow->w_ice->co[l-1][r][c]+
							c_liq*snow->w_liq->co[l-1][r][c]);
					}else{
						snow->T->co[l-1][r][c]=Tfreezing;
					}
					n=n-1;
					for(m=l;m<=n;m++){
						snow->Dzl->co[m][r][c]=snow->Dzl->co[m+1][r][c];
						snow->w_ice->co[m][r][c]=snow->w_ice->co[m+1][r][c];
						snow->w_liq->co[m][r][c]=snow->w_liq->co[m+1][r][c];
						snow->T->co[m][r][c]=snow->T->co[m+1][r][c];
					}
					snow->Dzl->co[n+1][r][c]=0.0;
					snow->w_ice->co[n+1][r][c]=0.0;
					snow->w_liq->co[n+1][r][c]=0.0;
					snow->T->co[n+1][r][c]=-99.0;
				}else{
					if(l!=n) t_error("Error 3 in snow combination");	//controllo
					snow->Dzl->co[l-1][r][c]+=snow->Dzl->co[l][r][c];
					snow->w_ice->co[l-1][r][c]+=snow->w_ice->co[l][r][c];
					snow->w_liq->co[l-1][r][c]+=snow->w_liq->co[l][r][c];
					snow->T->co[l-1][r][c]=snow->T->co[l][r][c];
					n=n-1;
					snow->Dzl->co[n+1][r][c]=0.0;
					snow->w_ice->co[n+1][r][c]=0.0;
					snow->w_liq->co[n+1][r][c]=0.0;
					snow->T->co[n+1][r][c]=-99.0;
				}
			}
		}

		for(l=1;l<=n;l++){
			if(snow->Dzl->co[l][r][c]<Dmin->co[n-l+1]) occuring=1;
		}

		for(l=1;l<=n;l++){
			if(n==max && l==1){
				occuring=0;
			}else if(snow->Dzl->co[l][r][c]>Dmax->co[n-l+1]){
				occuring=1;
			}
		}

	}while(occuring==1);

	snow->lnum->co[r][c]=n;

}*/

//TO HERE (comment if you want the new snow scheme)






}




/*==================================================================================================================*/
/*==================================================================================================================*/
/*4. subroutine GLAC_LAYER_COMBINATION*/
/*==================================================================================================================*/
/*==================================================================================================================*/

void glac_layer_combination(long r, long c, GLACIER *glac, double Ta, long max, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, double time)
	/* Author:  Stefano Endrizzi  Year:
	* function that finds the equilibrium number of glacier layers depending on the given layer depth
	* and on the allowed minimum and maximum depth.
	* If glac[l]<Dmin[n-l+1] => a layer has less ice than allowed, therefore its mass and energy
	* has to be shared with the lower one.
	* If glac[l]>Dmax[n-l+1] => a layer has more ice than allowed, therefore a new glacier layer has
	* to be initialized.
	* The routine goes on until an equilibrium is found.
	* Input:	r: row of the pixel in the raster
	* 			c: column of the pixel in the raster
	* 			Dmin: minimum allowed glacier depth
	* 			Dmax: maximum allowed glacier depth
	* 			glac: glacier structure
	* 			Ta: not used
	* 			max: maximum admitted glacier layer
	* 			time: not used
	* Output: glac modified with correct number of glacier layers
	* comment: Matteo Dall'Amico, May 2009 */
//Dmin, Dmax (1) il piu' superficiale
//Dmin, Dmax (5) il piu' profondo

{

long l;
long n;// number of glacier layers
long m;
short occuring;
double h;// internal energy of the glacier layer
double D=0.0;// glacier depth [mm]

for(l=1;l<=max;l++){
	D+=glac->Dzl->co[l][r][c];
}


if( (D==0.0 && glac->T->co[1][r][c]>-98.999) || (D>0 && D<=0.1) ){// no glacier
	for(l=2;l<=max;l++){
		glac->Dzl->co[l][r][c]=0.0;
		glac->w_liq->co[l][r][c]=0.0;
		glac->w_ice->co[l][r][c]=0.0;
		glac->T->co[l][r][c]=-99.0;
	}
	glac->Dzl->co[1][r][c]=0.0;
	glac->w_liq->co[1][r][c]=0.0;
	glac->w_ice->co[1][r][c]=0.0;
	glac->T->co[1][r][c]=-99.0;

	glac->lnum->co[r][c]=0;

}else if(D>0.1){
	n=glac->lnum->co[r][c];
	do{

		occuring=0;

		for(l=n;l>=1;l--){	//dall'alto in basso

			if(n==max && l==1){
				occuring=0;
			}else if(glac->Dzl->co[l][r][c]>Dmax->co[n-l+1]){
				if(n<max){//posso tranquillamente aggiungere un layer in piu'
					n=n+1;
					for(m=n;m>=l+2;m--){
						glac->Dzl->co[m][r][c]=glac->Dzl->co[m-1][r][c];
						glac->w_ice->co[m][r][c]=glac->w_ice->co[m-1][r][c];
						glac->w_liq->co[m][r][c]=glac->w_liq->co[m-1][r][c];
						glac->T->co[m][r][c]=glac->T->co[m-1][r][c];
					}
					glac->Dzl->co[l][r][c]*=0.5;
					glac->w_ice->co[l][r][c]*=0.5;
					glac->w_liq->co[l][r][c]*=0.5;
					glac->Dzl->co[l+1][r][c]=glac->Dzl->co[l][r][c];
					glac->w_ice->co[l+1][r][c]=glac->w_ice->co[l][r][c];
					glac->w_liq->co[l+1][r][c]=glac->w_liq->co[l][r][c];
					glac->T->co[l+1][r][c]=glac->T->co[l][r][c];
				}else{//devo trovare il posto per un layer in piu'
					if(n!=max) t_error("Error 1 in glac combination");	//controllo
					if(l>=3){	//combino i 2 layer piu' bassi
						h=(c_ice*glac->w_ice->co[1][r][c] + c_liq*glac->w_liq->co[1][r][c])*(glac->T->co[1][r][c]-Tfreezing) +
								Lf*glac->w_liq->co[1][r][c];
						h+=(c_ice*glac->w_ice->co[2][r][c] + c_liq*glac->w_liq->co[2][r][c])*(glac->T->co[2][r][c]-Tfreezing) +
								Lf*glac->w_liq->co[2][r][c];
						glac->Dzl->co[1][r][c]+=glac->Dzl->co[2][r][c];
						glac->w_ice->co[1][r][c]+=glac->w_ice->co[2][r][c];
						glac->w_liq->co[1][r][c]+=glac->w_liq->co[2][r][c];
						if(h<0.0){
							glac->T->co[1][r][c]=Tfreezing + h/(c_ice*glac->w_ice->co[1][r][c]+c_liq*glac->w_liq->co[1][r][c]);
						}else if(h>Lf*glac->w_liq->co[1][r][c]){
							glac->T->co[1][r][c]=Tfreezing + (h-Lf*glac->w_liq->co[1][r][c])/(c_ice*glac->w_ice->co[1][r][c]+
								c_liq*glac->w_liq->co[1][r][c]);
						}else{
							glac->T->co[1][r][c]=Tfreezing;
						}
						for(m=2;m<=n-1;m++){
							glac->Dzl->co[m][r][c]=glac->Dzl->co[m+1][r][c];
							glac->w_ice->co[m][r][c]=glac->w_ice->co[m+1][r][c];
							glac->w_liq->co[m][r][c]=glac->w_liq->co[m+1][r][c];
							glac->T->co[m][r][c]=glac->T->co[m+1][r][c];
						}
						for(m=n;m>=l+1;m--){
							glac->Dzl->co[m][r][c]=glac->Dzl->co[m-1][r][c];
							glac->w_ice->co[m][r][c]=glac->w_ice->co[m-1][r][c];
							glac->w_liq->co[m][r][c]=glac->w_liq->co[m-1][r][c];
							glac->T->co[m][r][c]=glac->T->co[m-1][r][c];
						}
						glac->Dzl->co[l-1][r][c]*=0.5;
						glac->w_ice->co[l-1][r][c]*=0.5;
						glac->w_liq->co[l-1][r][c]*=0.5;
						glac->Dzl->co[l][r][c]=glac->Dzl->co[l-1][r][c];
						glac->w_ice->co[l][r][c]=glac->w_ice->co[l-1][r][c];
						glac->w_liq->co[l][r][c]=glac->w_liq->co[l-1][r][c];
						glac->T->co[l][r][c]=glac->T->co[l-1][r][c];
					}else{
						if(l!=2) t_error("Error 2 in glac combination");	//controllo
						glac->Dzl->co[2][r][c]*=0.5;
						glac->w_ice->co[2][r][c]*=0.5;
						glac->w_liq->co[2][r][c]*=0.5;
						h=(c_ice*glac->w_ice->co[1][r][c] + c_liq*glac->w_liq->co[1][r][c])*(glac->T->co[1][r][c]-Tfreezing) +
							Lf*glac->w_liq->co[1][r][c];
						h+=(c_ice*glac->w_ice->co[2][r][c] + c_liq*glac->w_liq->co[2][r][c])*(glac->T->co[2][r][c]-Tfreezing) +
							Lf*glac->w_liq->co[2][r][c];
						glac->Dzl->co[1][r][c]+=glac->Dzl->co[2][r][c];
						glac->w_ice->co[1][r][c]+=glac->w_ice->co[2][r][c];
						glac->w_liq->co[1][r][c]+=glac->w_liq->co[2][r][c];
						if(h<0.0){
							glac->T->co[1][r][c]=Tfreezing + h/(c_ice*glac->w_ice->co[1][r][c]+c_liq*glac->w_liq->co[1][r][c]);
						}else if(h>Lf*glac->w_liq->co[1][r][c]){
							glac->T->co[1][r][c]=Tfreezing + (h-Lf*glac->w_liq->co[1][r][c])/(c_ice*glac->w_ice->co[1][r][c]+
								c_liq*glac->w_liq->co[1][r][c]);
						}else{
							glac->T->co[1][r][c]=Tfreezing;
						}
					}
				}
			}

		}
		for(l=1;l<=n;l++){	//dal basso in alto

			if(glac->Dzl->co[l][r][c]<Dmin->co[n-l+1]){
				if(l==1 && n==2){
					glac->Dzl->co[1][r][c]+=glac->Dzl->co[2][r][c];
					glac->w_ice->co[1][r][c]+=glac->w_ice->co[2][r][c];
					glac->w_liq->co[1][r][c]+=glac->w_liq->co[2][r][c];
					glac->T->co[1][r][c]=glac->T->co[2][r][c];
					n=n-1;
					glac->Dzl->co[2][r][c]=0.0;
					glac->w_ice->co[2][r][c]=0.0;
					glac->w_liq->co[2][r][c]=0.0;
					glac->T->co[2][r][c]=-99.0;
				}else if(l==1 && n>2){
					h=(c_ice*glac->w_ice->co[1][r][c] + c_liq*glac->w_liq->co[1][r][c])*(glac->T->co[1][r][c]-Tfreezing) +
						Lf*glac->w_liq->co[1][r][c];
					h+=(c_ice*glac->w_ice->co[2][r][c] + c_liq*glac->w_liq->co[2][r][c])*(glac->T->co[2][r][c]-Tfreezing) +
						Lf*glac->w_liq->co[2][r][c];
					glac->Dzl->co[1][r][c]+=glac->Dzl->co[2][r][c];
					glac->w_ice->co[1][r][c]+=glac->w_ice->co[2][r][c];
					glac->w_liq->co[1][r][c]+=glac->w_liq->co[2][r][c];
					if(h<0.0){
						glac->T->co[1][r][c]=Tfreezing + h/(c_ice*glac->w_ice->co[1][r][c]+c_liq*glac->w_liq->co[1][r][c]);
					}else if(h>Lf*glac->w_liq->co[1][r][c]){
						glac->T->co[1][r][c]=Tfreezing + (h-Lf*glac->w_liq->co[1][r][c])/(c_ice*glac->w_ice->co[1][r][c]+
							c_liq*glac->w_liq->co[1][r][c]);
					}else{
						glac->T->co[1][r][c]=Tfreezing;
					}
					n=n-1;
					for(m=2;m<=n;m++){
						glac->Dzl->co[m][r][c]=glac->Dzl->co[m+1][r][c];
						glac->w_ice->co[m][r][c]=glac->w_ice->co[m+1][r][c];
						glac->w_liq->co[m][r][c]=glac->w_liq->co[m+1][r][c];
						glac->T->co[m][r][c]=glac->T->co[m+1][r][c];
					}
					glac->Dzl->co[n+1][r][c]=0.0;
					glac->w_ice->co[n+1][r][c]=0.0;
					glac->w_liq->co[n+1][r][c]=0.0;
					glac->T->co[n+1][r][c]=-99.0;
				}else if(l!=n){
					h=(c_ice*glac->w_ice->co[l][r][c] + c_liq*glac->w_liq->co[l][r][c])*(glac->T->co[l][r][c]-Tfreezing) +
						Lf*glac->w_liq->co[l][r][c];
					h+=(c_ice*glac->w_ice->co[l-1][r][c] + c_liq*glac->w_liq->co[l-1][r][c])*(glac->T->co[l-1][r][c]-Tfreezing) +
						Lf*glac->w_liq->co[l-1][r][c];
					glac->Dzl->co[l-1][r][c]+=glac->Dzl->co[l][r][c];
					glac->w_ice->co[l-1][r][c]+=glac->w_ice->co[l][r][c];
					glac->w_liq->co[l-1][r][c]+=glac->w_liq->co[l][r][c];
					if(h<0.0){
						glac->T->co[l-1][r][c]=Tfreezing + h/(c_ice*glac->w_ice->co[l-1][r][c]+c_liq*glac->w_liq->co[l-1][r][c]);
					}else if(h>Lf*glac->w_liq->co[l-1][r][c]){
						glac->T->co[l-1][r][c]=Tfreezing + (h-Lf*glac->w_liq->co[l-1][r][c])/(c_ice*glac->w_ice->co[l-1][r][c]+
							c_liq*glac->w_liq->co[l-1][r][c]);
					}else{
						glac->T->co[l-1][r][c]=Tfreezing;
					}
					n=n-1;
					for(m=l;m<=n;m++){
						glac->Dzl->co[m][r][c]=glac->Dzl->co[m+1][r][c];
						glac->w_ice->co[m][r][c]=glac->w_ice->co[m+1][r][c];
						glac->w_liq->co[m][r][c]=glac->w_liq->co[m+1][r][c];
						glac->T->co[m][r][c]=glac->T->co[m+1][r][c];
					}
					glac->Dzl->co[n+1][r][c]=0.0;
					glac->w_ice->co[n+1][r][c]=0.0;
					glac->w_liq->co[n+1][r][c]=0.0;
					glac->T->co[n+1][r][c]=-99.0;
				}else{
					if(l!=n) t_error("Error 3 in glac combination");	//controllo
					glac->Dzl->co[l-1][r][c]+=glac->Dzl->co[l][r][c];
					glac->w_ice->co[l-1][r][c]+=glac->w_ice->co[l][r][c];
					glac->w_liq->co[l-1][r][c]+=glac->w_liq->co[l][r][c];
					glac->T->co[l-1][r][c]=glac->T->co[l][r][c];
					n=n-1;
					glac->Dzl->co[n+1][r][c]=0.0;
					glac->w_ice->co[n+1][r][c]=0.0;
					glac->w_liq->co[n+1][r][c]=0.0;
					glac->T->co[n+1][r][c]=-99.0;
				}
			}


		}
		for(l=1;l<=n;l++){
			if(glac->Dzl->co[l][r][c]<Dmin->co[n-l+1]) occuring=1;
		}

		for(l=1;l<=n;l++){
			if(n==max && l==1){
				occuring=0;
			}else if(glac->Dzl->co[l][r][c]>Dmax->co[n-l+1]){
				occuring=1;
			}
		}

	}while(occuring==1);

	glac->lnum->co[r][c]=n;

}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


double DEPTH(long r, long c, LONGMATRIX *n, DOUBLETENSOR *Dz){

	double d=0.0;
	long l;

	for(l=1;l<=n->co[r][c];l++){
		d+=Dz->co[l][r][c];
	}

	return(d);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double get_SWE(long r, long c, LONGMATRIX *n, DOUBLETENSOR *w1, DOUBLETENSOR *w2){

	double d=0.0;
	long l;

	for(l=1;l<=n->co[r][c];l++){
		d+=(w1->co[l][r][c]+w2->co[l][r][c]);
	}

	return(d);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void snowlayer_merging(long r, long c, SNOW *snow, long l1, long l2, long lres){
	/* Author:  Stefano Endrizzi  Year:
	* function that merges two snow layers in one based accounting for the mass and the energy of the system.
	* Input:	r: row of the pixel in the raster
	* 			c: column of the pixel in the raster
	* 			l1: first layer to merge
	* 			l2: second layer to merge
	* 			lres: resulting merged layer
	* Output: snow structure modified with correct merged layers
	* comment: Matteo Dall'Amico, May 2009 */
	double h;
/*printf("\nr=%ld, c=%ld, l1=%ld, l2=%ld ",r,c,l1,l2);
printf("snow->w_ice[l1]=%f, snow->w_liq[l1]=%f, snow->T[l1]=%f",snow->w_ice->co[l1][r][c],snow->w_liq->co[l1][r][c],snow->T->co[l1][r][c]);stop_execution();
printf("snow->T[l2]=%f",snow->T->co[l2][r][c]);stop_execution();
printf("snow->w_liq[l2]=%f",snow->w_liq->co[l2][r][c]);stop_execution();
printf("snow->w_ice[l2]=%f",snow->w_ice->co[l2][r][c]);stop_execution();*/
	h=internal_energy(snow->w_ice->co[l1][r][c],snow->w_liq->co[l1][r][c],snow->T->co[l1][r][c])+internal_energy(snow->w_ice->co[l2][r][c],snow->w_liq->co[l2][r][c],snow->T->co[l2][r][c]);
	snow->Dzl->co[lres][r][c]=snow->Dzl->co[l1][r][c]+snow->Dzl->co[l2][r][c];
	snow->w_ice->co[lres][r][c]=snow->w_ice->co[l1][r][c]+snow->w_ice->co[l2][r][c];
	snow->w_liq->co[lres][r][c]=snow->w_liq->co[l1][r][c]+snow->w_liq->co[l2][r][c];
	if(snow->Dzl->co[lres][r][c]<0 || snow->w_ice->co[lres][r][c]<0 || snow->w_liq->co[lres][r][c]<0){
		printf("ERROR 1 in snow layer merging r:%ld c:%ld l1:%ld l2:%ld lres:%ld\n",r,c,l1,l2,lres);
		write_snow_all(r, c, snow);
		stop_execution();
	}
	from_internal_energy(r+6000,c+6000,h,&(snow->w_ice->co[lres][r][c]),&(snow->w_liq->co[lres][r][c]),&(snow->T->co[lres][r][c]));
	if(snow->T->co[lres][r][c]>0){
		printf("ERROR 2 in snow layer merging r:%ld c:%ld l1:%ld l2:%ld lres:%ld\n",r,c,l1,l2,lres);
		write_snow_all(r, c, snow);
		stop_execution();
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double internal_energy(double w_ice, double w_liq, double T){

	double h;

	h=(c_ice*w_ice+c_liq*w_liq)*(T-Tfreezing) + Lf*w_liq;

	return(h);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void from_internal_energy(long r, long c, double h, double *w_ice, double *w_liq, double *T){
	/* Author:  Stefano Endrizzi  Year:
	* function that merges two snow layers in one based accounting for the mass and the energy of the system.
	* Input:	r: row of the pixel in the raster
	* 			c: column of the pixel in the raster
	* 			h: internal energy of the layer
	* Output:
	* 		wice : equivalent mass of liquid water stored in the snow per unit of surface [kg/m2]
	* 		wliq : mass of liquid water present in the snow per unit of surface [kg/m2]
	* 			T: temperature of the snow layer
	* comment: Matteo Dall'Amico, May 2009 */
	double SWE=(*w_ice)+(*w_liq);

	if(h<0){
		*w_ice=SWE;
		*w_liq=0.0;
		*T=Tfreezing + h/(c_ice*SWE);
		if(*w_ice<0 || *w_liq<0 || *T>Tfreezing) printf("Error 1 in H calculation r:%ld c:%ld : wice:%f wliq:%f T:%f h:%f SWE:%f Lf*SWE:%f\n",r,c,*w_ice,*w_liq,*T,h,SWE,Lf*SWE);
	}else if(0.0<=h && h<Lf*SWE){
		*T=Tfreezing;
		*w_liq=h/Lf;
		*w_ice=SWE-(*w_liq);
		if(*w_ice<0 || *w_liq<0 || *T>Tfreezing) printf("Error 2 in H calculation r:%ld c:%ld : wice:%f wliq:%f T:%f h:%f SWE:%f Lf*SWE:%f\n",r,c,*w_ice,*w_liq,*T,h,SWE,Lf*SWE);
	}else{
		*w_liq=SWE;
		*w_ice=0.0;
		*T=Tfreezing;
		if(*w_ice<0 || *w_liq<0 || *T>Tfreezing) printf("Error 3 in H calculation r:%ld c:%ld : wice:%f wliq:%f T:%f h:%f SWE:%f Lf*SWE:%f\n",r,c,*w_ice,*w_liq,*T,h,SWE,Lf*SWE);	
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/*void get_softsnow(SNOW *snow, double **Z){

	long l,r,c,nr=snow->T->nrh,nc=snow->T->nch;
	short wet;

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(Z[r][c]!=UV->V->co[2]){
				snow->softSWE->co[r][c]=0.0;
				for(l=snow->lnum->co[r][c];l>=1;l--){
					if(l==snow->lnum->co[r][c]) wet=0;
					if(snow->w_liq->co[l][r][c]>0) wet=1;
					if(wet==0) snow->softSWE->co[r][c]+=snow->w_ice->co[l][r][c];
				}
			}
		}
	}
}*/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void write_snow(long r, long c, long l, SNOW *snow){

	printf("r:%ld c:%ld wice(%ld/%ld):%f wliq(%ld):%f T(%ld):%f Dz(%ld):%f\n",r,c,l,snow->lnum->co[r][c],snow->w_ice->co[l][r][c],l,snow->w_liq->co[l][r][c],
			l,snow->T->co[l][r][c],l,snow->Dzl->co[l][r][c]);

}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void write_snow_all(long r, long c, SNOW *snow){
	long l;
	for(l=1;l<=snow->lnum->co[r][c];l++){
		write_snow(r,c,l,snow);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void set_snow(long r, long c, double *wliq, double *wice, double *T, double *D, long l1, long l2, double Dlim){

	double h, f, dwl, dwi;

	if(D[l1]<Dlim){	//l1 too shallow and takes mass from l2
		h=internal_energy(wice[l1],wliq[l1],T[l1]);
		f=(Dlim-D[l1])/D[l2];
		dwl=f*wliq[l2];
		wliq[l1]+=dwl;
		wliq[l2]-=dwl;
		dwi=f*wice[l2];
		wice[l1]+=dwi;
		wice[l2]-=dwi;
		D[l2]-=(Dlim-D[l1]);
		D[l1]=Dlim;
		h+=internal_energy(dwi, dwl, T[l2]);
		from_internal_energy(r+1000,c+1000,h,&(wice[l1]),&(wliq[l1]),&(T[l1]));

	}else if(D[l1]>Dlim){	//l1 too thick and gaves mass to l2
		h=internal_energy(wice[l2],wliq[l2],T[l2]);
		f=(D[l1]-Dlim)/D[l1];
		dwl=f*wliq[l1];
		wliq[l1]-=dwl;
		wliq[l2]+=dwl;
		dwi=f*wice[l1];
		wice[l1]-=dwi;
		wice[l2]+=dwi;
		D[l2]+=(D[l1]-Dlim);
		D[l1]=Dlim;
		h+=internal_energy(dwi, dwl, T[l1]);
		from_internal_energy(r+2000,c+2000,h,&(wice[l2]),&(wliq[l2]),&(T[l2]));
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void split_layers(long r, long c, SNOW *snow, long l1){

	long l;

	if(l1>snow->lnum->co[r][c]) t_error("Error 1 in split_layers");

	snow->w_ice->co[l1][r][c]*=0.5;
	snow->w_liq->co[l1][r][c]*=0.5;
	snow->Dzl->co[l1][r][c]*=0.5;

	for(l=snow->lnum->co[r][c];l>=l1;l--){
		snow->w_ice->co[l+1][r][c]=snow->w_ice->co[l][r][c];
		snow->w_liq->co[l+1][r][c]=snow->w_liq->co[l][r][c];
		snow->T->co[l+1][r][c]=snow->T->co[l][r][c];
		snow->Dzl->co[l+1][r][c]=snow->Dzl->co[l][r][c];
	}

	snow->lnum->co[r][c]+=1;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void merge_layers(long r, long c, SNOW *snow, long l1){

	long l;

	if(l1>snow->lnum->co[r][c]) t_error("Error 1 in merge_layers");

	if(l1==snow->lnum->co[r][c]){
		snowlayer_merging(r, c, snow, l1, l1-1, l1-1);
		initialize_snow(r, c, l1, snow);
	}else{
		snowlayer_merging(r, c, snow, l1, l1+1, l1);
		for(l=l1+1;l<snow->lnum->co[r][c];l++){
			snow->w_ice->co[l][r][c]=snow->w_ice->co[l+1][r][c];
			snow->w_liq->co[l][r][c]=snow->w_liq->co[l+1][r][c];
			snow->T->co[l][r][c]=snow->T->co[l+1][r][c];
			snow->Dzl->co[l][r][c]=snow->Dzl->co[l+1][r][c];
		}
		initialize_snow(r, c, snow->lnum->co[r][c], snow);
	}
	snow->lnum->co[r][c]-=1;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void min_max_layer(long n, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, DOUBLEVECTOR *Dmin2, DOUBLEVECTOR *Dmax2, long linf){
	/* Author:  Stefano Endrizzi  Year:
	* function that finds the range of minimum and maximum allowed layer thickness of a new snow layer
	* Input:
	* 		n:	number of snow layers
	* 	 Dmin:	vector of minimum snow layer thickness
	* 	 Dmax:	vector of maximum snow layer thickness
	* 	 linf:	layer of unlimited thickness (beginning from below)
	* Output:
	* 	 Dmin2:	vector of new minimum snow layer thickness
	* 	 Dmax2: vector of new maximum snow layer thickness
	* comment: Matteo Dall'Amico, May 2009 */
	long l,mup,mdw,N=Dmin->nh;

	if(n==N){// number of snow layers equals the maximum admitted snow layers

		for(l=1;l<=n;l++){
			Dmin2->co[l]=Dmin->co[l];
			Dmax2->co[l]=Dmax->co[l];
		}

	}else{

		mup=ceil(n/2.0);
		mdw=floor(n/2.0);

		//printf("a. n=%ld linf=%ld mup=%ld mdw=%ld\n",n,linf,mup,mdw);stop_execution();

		if(linf<N-mup && linf<=mdw){
			mdw=linf-1;
			mup=n-mdw;
		}else if(linf>mdw && linf>=N-mup){
			mup=N-linf;
			mdw=n-mup;
		}

		//printf("b. n=%ld linf=%ld mup=%ld mdw=%ld\n",n,linf,mup,mdw);

		for(l=1;l<=mdw;l++){
			Dmin2->co[l]=Dmin->co[l];
			Dmax2->co[l]=Dmax->co[l];
		}

		for(l=n;l>n-mup;l--){
			Dmin2->co[l]=Dmin->co[N+l-n];
			Dmax2->co[l]=Dmax->co[N+l-n];
		}
		
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initialize_snow(long r, long c, long l, SNOW *snow){

	snow->w_ice->co[l][r][c]=0.0;
	snow->w_liq->co[l][r][c]=0.0;
	snow->Dzl->co[l][r][c]=0.0;
	snow->T->co[l][r][c]=-99.0;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void show_Dminmax(long r, long c, double *Dmin, double *Dmax, long n){

	long l;

	for(l=1;l<=n;l++){
		printf("l:%ld ltot:%ld Dmin:%f Dmax:%f\n",l,n,Dmin[l],Dmax[l]);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_snow_age(double Psnow, double Ts, double Dt, double *tsnow){

	double r1, r2, r3;

	//effect of grain growth due to vapour diffusion
	r1=exp(5000.0*(1.0/tk-1.0/(Ts+tk)));

	//effect melt and refreezing*/
	r2=pow(r1,10);
	if(r2>1.0) r2=1.0;

	//effect of dirt
	r3=0.3;

	//non-dimensional snow age: 10 mm of snow precipitation restore snow age Dt(s)
	*tsnow=(*tsnow+(r1+r2+r3)*Dt*1.0E-6)*(1.0-Psnow/10.0);
	if(*tsnow<0.0) *tsnow=0.0;
	if((*tsnow)!=(*tsnow)) printf("tsnow no value - tausn:%f P:%f Ts:%f r1:%f r2:%f r3:%f\n",*tsnow,Psnow,Ts,r1,r2,r3);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double snow_albedo(double ground_alb, double snowD, double AEP, double freshsnow_alb, double C, double tsnow, double cosinc, double ( *F)(double x)){
	double A, Fage=1.0-1.0/(1.0+tsnow), w;
	A=freshsnow_alb*(1.0-C*Fage);
	A+=0.4*(1.0-A)*(*F)(cosinc);
	if(snowD<AEP){	//if snow is shallow (<AEP), interpolate between snow and ground albedo
		w=(1.0-snowD/AEP)*exp(-snowD*0.5/AEP);
		A=w*ground_alb+(1.0-w)*A;
	}
	return(A);
}

double Fzen(double cosinc){
	double f, b=2.0;
	if(cosinc<0.5){
		f=(1.0/b)*((b+1.0)/(1.0+2.0*b*cosinc)-1.0);
	}else{
		f=0.0;
	}
	return(f);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//Too shallow snow layers can give numerical instabilities. Just get rid of them at the beginning if there is enough energy

void set_shallow_snowpack(long r, long c, double Dt, SNOW *snow, double *SW, double *Mr, long *n){

	if(snow->type->co[r][c]==1){

		//printf("D:%f %f %f %f\n",snow->Dzl->co[1][r][c],SW[1],SW[2],Lf*snow->w_ice->co[1][r][c]/Dt);

		if(SW[1]+SW[2]>Lf*snow->w_ice->co[1][r][c]/Dt){

			//printf("%f %f %f\n",SW[1],SW[2],Lf*snow->w_ice->co[1][r][c]/Dt);

			*Mr=snow->w_ice->co[1][r][c]/Dt; //[mm/s]
			SW[1]+=(SW[2]-Lf*snow->w_ice->co[1][r][c]/Dt);
			SW[2]=0.0;

			//printf("%f %f\n",SW[1],SW[2]);
			//stop_execution();

			snow->w_ice->co[1][r][c]=0.0;
			snow->w_liq->co[1][r][c]=0.0;
			snow->T->co[1][r][c]=-99.0;
			snow->Dzl->co[1][r][c]=0.0;
			snow->lnum->co[r][c]=0;
			snow->type->co[r][c]=0;

			*n=0;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

