
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
#include "times.h"
#include "output.09375.h"
#include "snow.09375.h"
#include "extensions.h"
#include "rw_maps.h"

#include "PBSM.h"
//#include "SnowTran.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

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

double d=2.0E-4,ks,K; /* modified by Emanuele Cordano on 24/9/9 */

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
ks=0.077*pow(d,2.0)*exp(-0.0078*rho);
K=ks*rho_w*rho_w*g/mu_l;

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
double theta_i,theta_w,rho,c1,c2,c3,c4,c5,load,eta0,eta;

theta_i=snow->w_ice->co[l][r][c]/(0.001*snow->Dzl->co[l][r][c]*rho_i);
theta_w=snow->w_liq->co[l][r][c]/(0.001*snow->Dzl->co[l][r][c]*rho_w);

rho=(snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c])/(0.001*snow->Dzl->co[l][r][c]);

if( theta_i < par->snow_maxpor ){

	//DESTRUCTIVE METAMORPHISM
	//gamma=theta*rho
	if(theta_i*rho_i<=par->snow_density_cutoff){
		c1=par->drysnowdef_rate;
	}else{
		c1=par->drysnowdef_rate*exp(-0.046*(rho_i*theta_i - par->snow_density_cutoff));
	}
	if(theta_w>0.001) c1*=par->wetsnowdef_rate;
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

	//snow->Dzl->co[l][r][c]*=(1.0 + ((*CR1)+(*CR2))*par->Dt);// was like this
	snow->Dzl->co[l][r][c]*=Fmax(0.1, 1.0 + ((*CR1)+(*CR2))*par->Dt);// new by stefano on 18/11/09

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

long l, cont=0, max=Dmin->nh;
short occuring;
double D0=0.0,D1=0.0,D=0.0;
DOUBLEVECTOR *Dmin2,*Dmax2,*wice,*wliq,*Dz,*Temp;

//D=snow depth(mm)
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

		//if it is ok, it is done; otherwise try to spplit or merge layers
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

long l,n,m;
short occuring;
double h;// internal energy of the glacier layer
double D=0.0;// glacier depth [mm]

//D=glac depth(mm)
for(l=1;l<=max;l++){
	D+=glac->Dzl->co[l][r][c];
}


if( (D==0.0 && glac->T->co[1][r][c]>-98.999) || (D>0 && D<=0.1) ){
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
	double a=alpha_snow, T0;
	double A=-h/(c_ice*SWE), B=c_liq/(c_ice*a*a), C=(Lf*SWE-h)/(c_ice*SWE*a*a);
	long cont=0;

	if(SWE>0){
		*T=0.0;
		do{
			T0=*T;
			*T=T0-(pow(T0,3.0)+A*pow(T0,2.0)+B*T0+C)/(3.0*pow(T0,2.0)+2*A*T0+B);
			cont++;
			if(cont>100) printf("%ld %e %e %e %e\n",cont,SWE,h,T0,*T);
		}while(fabs(*T-T0)>1.E-10 && cont<100);

		*w_liq=theta_snow(*T)*SWE;
		*w_ice=SWE-(*w_liq);
	}else{
		*T=0.0;
		*w_liq=0.0;
		*w_ice=0.0;
	}

	/*if(h<0){
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
	}*/


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

		if(linf<N-mup && linf<=mdw){
			mdw=linf-1;
			mup=n-mdw;
		}else if(linf>mdw && linf>=N-mup){
			mup=N-linf;
			mdw=n-mup;
		}

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

void update_snow_age(double Psnow, double Ts, double Dt, double *tsnow_dim, double *tsnow_nondim){

	double r1, r2, r3;

	//increase snow age
	*tsnow_dim=Fmax( 0.0, (*tsnow_dim+Dt)*(1.0-Psnow/10.0) );

	//effect of grain growth due to vapour diffusion
	r1=exp(5000.0*(1.0/tk-1.0/(Ts+tk)));

	//effect melt and refreezing*/
	r2=pow(r1,10);
	if(r2>1.0) r2=1.0;

	//effect of dirt
	r3=0.3;

	//non-dimensional snow age: 10 mm of snow precipitation restore snow age Dt(s)
	*tsnow_nondim=Fmax( 0.0, (*tsnow_nondim+(r1+r2+r3)*Dt*1.0E-6)*(1.0-Psnow/10.0) );
	if((*tsnow_nondim)!=(*tsnow_nondim)) printf("tsnow no value - tausn:%f P:%f Ts:%f r1:%f r2:%f r3:%f\n",*tsnow_nondim,Psnow,Ts,r1,r2,r3);

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

		if(SW[1]+SW[2]>Lf*snow->w_ice->co[1][r][c]/Dt){

			*Mr=snow->w_ice->co[1][r][c]/Dt; //[mm/s]
			SW[1]+=(SW[2]-Lf*snow->w_ice->co[1][r][c]/Dt);
			SW[2]=0.0;

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


void snow_properties(long r, long c, long beg, long end, double *D, double *wl, double *wi, double *k, double *T, double ***tdz, double ***twl, double ***twi,
					 double ***tT, double (* kfunct)(double rho)){

	long l,m;
	double rho;

	for(l=beg;l<=end;l++){

		m=end-l+1;

		D[l]=1.0E-3*tdz[m][r][c];
		wl[l]=twl[m][r][c];
		wi[l]=twi[m][r][c];
		T[l]=tT[m][r][c];
		//C[l]=(c_ice*wi[l]+c_liq*wl[l])/D[l];

		//printf("l:%ld wi:%f wl:%f\n",l,wi[l],wl[l]);

		rho=(wl[l]+wi[l])/D[l];
		k[l]=(*kfunct)(rho);
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

void liqWBsnow(long r, long c, SNOW *snow, double *Mr, double *PonS, PAR *par, double slope, double P, double *wi, double *wl, double *T, double Edt){

	long l, m;
	double Se, theta_liq, thice, Dwl=P, Dwi, wice_old, rho;
	//double h;

	//rain on snow
	*PonS=0.0;	//initialization


	//a. ordinary case
	if(snow->type->co[r][c]==2){

		for(l=snow->lnum->co[r][c];l>=1;l--){

			m=snow->lnum->co[r][c] - l + 1;	//increasing downwards (while l increases upwards)

			//SUBTRACT SUBLIMATION in the upper layer
			if(m==1){
				Dwi=Fmax(0.0, wi[m]-Edt) - wi[m];	//negative
			}else{
				Dwi=0.0;
			}

			//find new internal energy, taking into account the latent heat flux of the incoming water flux
			//h = internal_energy(wi[m], wl[m], T[m]);

			wl[m]+=Dwl;
			wi[m]+=Dwi;

			//from_internal_energy(r, c, h + internal_energy(Dwi, 0.0, T[m]) + Lf*Dwl, &(wi[m]), &(wl[m]), &(T[m]));

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
				Dwl=snow->w_liq->co[l][r][c];
				snow->w_liq->co[l][r][c]=0.0;
			}else{
				theta_liq=snow->w_liq->co[l][r][c]/(1.0E-3*snow->Dzl->co[l][r][c]*rho_w);		//[-]
				Se=(theta_liq - par->Sr*(1.0-thice))/( (1.0-thice) - par->Sr*(1.0-thice));
				if(Se<0) Se=0.0;
				if(Se>1) Se=1.0;
				rho=snow->w_ice->co[l][r][c]/(0.001*snow->Dzl->co[l][r][c]);
				if(theta_liq>par->Sr*(1.0-thice)){
					Dwl=Fmin(5.0*pow(Se,3.0)*par->Dt , (theta_liq - par->Sr*(1.0-thice))*1.0E-3*snow->Dzl->co[l][r][c]*rho_w);
					//q=(theta_liq - par->Sr*(1.0-thice))*1.0E-3*snow->Dzl->co[l][r][c]*rho_w/par->Dt;
				}else{
					Dwl=0.0;
				}
				snow->w_liq->co[l][r][c]-=Dwl;
				//printf("l:%ld theta:%f thice:%f Dwl:%f\n",l,theta_liq,thice,Dwl);

				/*if(theta_liq>par->Sr*(1.0-thice)){
					q=(theta_liq - par->Sr*(1.0-thice))*1.0E-3*snow->Dzl->co[l][r][c]*rho_w/par->Dt;
					snow->w_liq->co[l][r][c]-=q*par->Dt;
					if(q!=q) stop_execution();
				}else{
					q=0.0;
				}*/

				/*if(theta_liq>par->Sr*(1.0-thice)){		//when theta_liq is larger than the capillary retenction
					snow->w_liq->co[l][r][c]-=q->co[l]*par->Dt;
					theta_liq=snow->w_liq->co[l][r][c]/(1.0E-3*snow->Dzl->co[l][r][c]*rho_w);
					if(theta_liq<=par->Sr*(1.0-thice)){
						q->co[l]-=((par->Sr*(1.0-thice) - theta_liq)*1.0E-3*snow->Dzl->co[l][r][c]*rho_w)/par->Dt;
						snow->w_liq->co[l][r][c]=par->Sr*(1.0-thice)*1.0E-3*snow->Dzl->co[l][r][c]*rho_w;
					}else if(theta_liq > (1.0-thice) ){
						snow->w_liq->co[l][r][c]=(1.0-thice)*rho_w*1.0E-3*snow->Dzl->co[l][r][c];
						q->co[l]+=( theta_liq - (1.0-thice) )*(rho_w*1.0E-3*snow->Dzl->co[l][r][c])/par->Dt;  //liquid excess in snow goes downwards
					}
				}*/
			}

		}

		//melting rate
		*Mr=(Dwl-P)/par->Dt;	//[mm/s]
		//printf("Mr:%f %f %f\n",*Mr,Dwl,P);

		//rain on snow
		*PonS=P;	//[mm]



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

	//stop_execution();
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void iceWBsnow(long r, long c, SNOW *snow, double P, double T){

	long ns;
	double Dz;
	//double h;

	if(P>0){

		Dz=P*rho_w/snow->rho_newsnow->co[r][c];	//snow depth addition due to snow precipitation

		if(snow->type->co[r][c]==0){

			snow->Dzl->co[1][r][c]+=Dz;
			snow->w_ice->co[1][r][c]+=P;

		}else{

			ns=snow->lnum->co[r][c];

			//h=internal_energy(snow->w_ice->co[ns][r][c], snow->w_liq->co[ns][r][c], snow->T->co[ns][r][c]);
			//h+=(c_ice*P)*(Fmin(T, -0.1) - Tfreezing);

			snow->Dzl->co[ns][r][c]+=Dz;
			snow->w_ice->co[ns][r][c]+=P;

			//from_internal_energy(r, c, h, &(snow->w_ice->co[ns][r][c]),&(snow->w_liq->co[ns][r][c]),&(snow->T->co[ns][r][c]));

		}


	}
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

void WBglacier(long ns, long r, long c, GLACIER *glac, double *Mr, PAR *par, double *wi, double *wl, double *T, double Edt){

	long l, m;
	double Dz, thice, theta_liq, Dwi;
	//double h;

	*Mr=0.0;

	for(l=glac->lnum->co[r][c];l>=1;l--){

		m=ns + glac->lnum->co[r][c] - l + 1;

		//SUBTRACT SUBLIMATION in the upper layer
		if(m==1){
			Dwi=Fmax(0.0, wi[m]-Edt) - wi[m];
		}else{
			Dwi=0.0;
		}

		//h = internal_energy(wi[m], wl[m], T[m]);

		wi[m]+=Dwi;

		//from_internal_energy(r, c, h + internal_energy(Dwi, 0.0, T[m]), &(wi[m]), &(wl[m]), &(T[m]));

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

void output_snow(SNOW *snow, double **Z, PAR *par){

	long r, c, l, nr, nc;
	double D, n=(par->output_snow*3600.0)/(par->Dt);

	nr=snow->T->nrh;
	nc=snow->T->nch;

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(Z[r][c]!=UV->V->co[2]){

				//PBSM
				if(par->blowing_snow==1) print_windtrans_snow(r, c, snow, par);

				//SnowTran3D
				//if(par->blowing_snow==1) print_windtrans_snow2(r, c, snow, par);

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

void find_SCA(short **LC2, SNOW *snow, PAR *par, double **Z, double t){

	long l, r, c, cont=0, conttot=0, d2, mo2, y2, h2, mi2;
	double T, D, SWE, Tmean=0.0, Tsmean=0.0, Dmean=0.0, SWEmean=0.0, SCA, JD;
	FILE *f;

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if(Z[r][c]!=UV->V->co[2] && LC2[r][c]>0){
				D=0.0; T=0.0; SWE=0.0;
				for(l=1;l<=snow->lnum->co[r][c];l++){
					D+=snow->Dzl->co[l][r][c];
					SWE+=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c]);
					T+=snow->T->co[l][r][c]*snow->Dzl->co[l][r][c];
				}
				if(D>0)T/=D;

				conttot++;
				if(D>10){
					cont++;
					Dmean+=D;
					SWEmean+=SWE;
					Tmean+=T;
					Tsmean+=snow->T->co[snow->lnum->co[r][c]][r][c];
				}

			}
		}
	}

	if(cont>0){
		Dmean/=(double)cont;
		SWEmean/=(double)cont;
		Tmean/=(double)cont;
		Tsmean/=(double)cont;
		SCA=cont/(double)conttot;
	}else{
		SCA=0.0;
	}

	date_time(t, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
	f=fopen(join_strings(files->co[fHpatch]+1,textfile),"a");
	write_date(f, d2, mo2, y2, h2, mi2);
	fprintf(f,",%f,%f,%f,%f,%f,%f,%f\n",JD+(double)(daysfrom0(y2)),Dmean,SWEmean,Tmean,Tsmean,(1.0-SCA)*100.0,SCA*100.0);
	fclose(f);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snow_fluxes_H(short **LC2, SNOW *snow, DOUBLEMATRIX *Ts, DOUBLEMATRIX *H, DOUBLEMATRIX *Z, DOUBLEMATRIX *Ta, TIMES *times, PAR *par, double *VSFA, double *Hv){

	long n,l,r,c,contLC,cont,contT,rmin,cmin,rmax,cmax;
	double Hsca,Hsfa,Hmax,Hmin,Hrsm,Tsmax,Tsmin,Tsmean,Tamean,SFA;
	double Hvmean,Hvmax,Hvmin,Hvrsm,Tsvmean;
	double Dmean, D;
	double JD;
	long d2,mo2,y2,h2,mi2;

	char *temp; /* added by Emanuele Cordano on 21 September 2009 */
	FILE *f;

	for(n=1;n<=par->nLC;n++){
		rmin=0;cmin=0;rmax=0;cmax=0;contLC=0;cont=0;contT=0;
		Hsca=0.0;Hsfa=0.0;Hmax=-10000.0;Hmin=10000.0;Hrsm=0.0;Tsmax=-300.0;Tsmin=300.0;Tsmean=0.0;Tamean=0.0;
		Hvmean=0.0;Hvmax=-10000.0;Hvmin=10000.0;Hvrsm=0.0;Tsvmean=0.0;Dmean=0.0;

		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(Z->co[r][c]!=UV->V->co[2]){
					if(LC2[r][c]==n){

						contLC++;
						Tamean+=Ta->co[r][c];

						D=0.0;
						for(l=1;l<=snow->lnum->co[r][c];l++){ D+=snow->Dzl->co[l][r][c];}

						if(D>10){
							cont++;
							Dmean+=D;
							Hsca+=H->co[r][c];
						}else{
							Hsfa+=H->co[r][c];
							Tsmean+=Ts->co[r][c];
							if(H->co[r][c]>Hmax) Hmax=H->co[r][c];
							if(H->co[r][c]<Hmin) Hmin=H->co[r][c];
							if(Ts->co[r][c]>Tsmax) Tsmax=Ts->co[r][c];
							if(Ts->co[r][c]<Tsmin) Tsmin=Ts->co[r][c];
						}

						if(D==0 && H->co[r][c]>0){
							contT++;
							Hvmean+=H->co[r][c];
							Tsvmean+=Ts->co[r][c];
							if(H->co[r][c]>Hvmax) Hvmax=H->co[r][c];
							if(H->co[r][c]<Hvmin) Hvmin=H->co[r][c];
						}
					}
				}
			}
		}

		if(cont > 0) Dmean/=(double)cont;
		//printf("n:%ld Dmean:%f cont:%ld\n",n,Dmean,cont);

		if(cont > 0) Hsca/=(double)cont;
		if(cont < contLC) Hsfa/=(double)(contLC-cont);
		if(cont < contLC) Tsmean/=(double)(contLC-cont);
		Tamean/=(double)contLC;
		SFA=(contLC-cont)/(double)contLC;

		if(contT>0){
			Hvmean/=(double)contT;
			Tsvmean/=(double)contT;
		}

		if(cont<contLC){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(Z->co[r][c]!=UV->V->co[2]){
						if(snow->Dzl->co[1][r][c]==0.0){
							Hrsm+=pow(H->co[r][c]-Hsfa,2.0)/(double)(contLC-cont);
							if(Ts->co[r][c]==Tsmin){
								rmin=r;
								cmin=c;
							}
							if(Ts->co[r][c]==Tsmax){
								rmax=r;
								cmax=c;
							}
						}
					}
				}
			}
			Hrsm=pow(Hrsm,0.5);
		}

		if(contT>0){
			for(r=1;r<=Nr;r++){
				for(c=1;c<=Nc;c++){
					if(Z->co[r][c]!=UV->V->co[2]){
						if(snow->Dzl->co[1][r][c]==0.0) Hvrsm+=pow(H->co[r][c]-Hvmean,2.0)/(double)contT;
					}
				}
			}
			Hvrsm=pow(Hvrsm,0.5);
		}

		date_time(times->time, par->year0, par->JD0, 0.0, &JD, &d2, &mo2, &y2, &h2, &mi2);
		temp=namefile_i(files->co[fHpatch]+1,n);
		f=fopen(temp,"a");
		free(temp);
		write_date(f, d2, mo2, y2, h2, mi2);
		fprintf(f,",%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,MIN:%ld %ld,MAX:%ld %ld\n",JD+(double)(daysfrom0(y2)),Dmean,SFA,1-SFA,Hsca,Hsfa,Hmax,Hmin,Hrsm,Tsmean,Tsmin,Tsmax,Tamean,contT/(double)contLC,Hvmean,Hvmin,Hvmax,Hvrsm,Tsvmean,rmin,cmin,rmax,cmax);
		fclose(f);

		Hv[n]=Hvmean;
		VSFA[n]=contT/(double)contLC;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double adv_efficiency(double SFA){

	double Y;

	Y=-0.043+1.0859*exp(-2.9125*SFA);
	if(SFA>0.85) Y=0.0;

	return(Y);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double max_dtheta_snow(){

	double a=alpha_snow;
	double maxT;

	maxT=-pow(1./(3*a*a), 0.5);
	return(maxT);
}

double theta_snow(double T){

	double a=alpha_snow;
	double th;

	th=1./(1.+pow(a*T, 2.));
	if(T>0) th=1.;
	return(th);
}

double dtheta_snow(double T){

	double a=alpha_snow;
	double dth;

	dth=-2.*a*a*T/pow(1.+pow(a*T, 2.), 2.);
	if(T>0) dth=-2.*a*a*max_dtheta_snow()/pow(1.+pow(a*max_dtheta_snow(), 2.), 2.);
	return(dth);
}

