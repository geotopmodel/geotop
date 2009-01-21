

/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion KMackenzie

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 KMackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOtop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


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

//Jordan et al., 1999

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
double theta_i,c1,c2,c3,c4,c5,load,eta0,eta;

theta_i=snow->w_ice->co[l][r][c]/(0.001*snow->Dzl->co[l][r][c]*rho_i);

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

	//printf("l:%ld %e %e\n",l,*CR1,*CR2);

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
//Dmin, Dmax (1) il piu' superficiale
//Dmin, Dmax (5) il piu' profondo

{

long l,n,m;
short occuring;
double h,D=0.0;

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

void set_windtrans_snow(SNOW *snow, double **Z, double **Ta, double **P, double t, PAR *par, METEO *met){

	double DW,DW0,SWE,D,rho;
	long i,r,c,nr,nc,ns;
	FILE *f;

	nr=snow->T->nrh;
	nc=snow->T->nch;

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(Z[r][c]!=UV->V->co[2]){

				//printf("%ld %ld T:%f V:%f\n",r,c,met->Tgrid->co[r][c],met->Vgrid->co[r][c]);
				//stop_execution();

				D=DEPTH(r,c,snow->lnum,snow->Dzl);
				if(D!=D){
					printf("Novalue in set_windtrans_snow(a): r:%ld c:%ld SnowD:%f lnum:%ld\n",r,c,DEPTH(r, c, snow->lnum, snow->Dzl),snow->lnum->co[r][c]);
					write_snow_all(r, c, snow);
				}
				SWE=get_SWE(r,c,snow->lnum,snow->w_ice,snow->w_liq);
				ns=snow->lnum->co[r][c];

				DW=snow->Wsalt->co[r][c] + snow->Wsubl->co[r][c] + snow->Wsusp->co[r][c] + snow->Wsubgrid->co[r][c];

				//if(r==par->rc->co[1][1] && c==par->rc->co[1][2]) printf("DW:%f\n",DW);

				DW0=DW;

				if(SWE+DW-snow->ListonSWE->co[r][c]>1.E-4){
					f=fopen(files->co[ferr]+1,"a");
					fprintf(f,"Not accordance with Liston: r:%ld c:%ld nl:%ld D:%f rho:%f SWE:%f Lis1:%f List2:%f DW:%f DWsalt:%f DWsubl:%f DWsusp:%f DWsubgrid:%f\n",
						r,c,ns,D,SWE*1000/D,SWE,snow->ListonSWE->co[r][c]-DW,snow->ListonSWE->co[r][c],DW,snow->Wsalt->co[r][c],snow->Wsubl->co[r][c],
						snow->Wsusp->co[r][c],snow->Wsubgrid->co[r][c]);
					fclose(f);
				}

				if(snow->lnum->co[r][c]>0){	//snow on the soil

					if(DW<0){	//snow eroded

						i=ns;
						do{
							if(i<ns){

								if(snow->w_ice->co[i+1][r][c]<0){
									DW=snow->w_ice->co[i+1][r][c];
									snow->w_ice->co[i+1][r][c]=0.0;
									snow->Dzl->co[i+1][r][c]=0.0;
									snow->lnum->co[r][c]-=1;
								}
							}

							if(fabs(DW)>0.05*snow->w_ice->co[i][r][c] && snow->w_liq->co[i][r][c]>0){
								printf("\nPositive snow water content with snow erosion by wind in cell r:%ld c:%ld snowlayer:%ld\n",r,c,i);
								f=fopen(files->co[ferr]+1,"a");
								fprintf(f,"\nPositive snow water content with snow erosion by wind in cell r:%ld c:%ld snowlayer:%ld\n",r,c,i);
								fprintf(f,"DW:%f DWsalt:%f DWsubl:%f DWsusp:%f DWsubgrid:%f i:%ld ns:%ld Wliq:%f Wice:%f Dlayer:%f T:%f Psnow:%f softSWEbeg:%f softSWEend:%f D:%f SWE:%f SWEliston:%f DW:%f check:%f\n",
									DW,snow->Wsalt->co[r][c],snow->Wsubl->co[r][c],snow->Wsusp->co[r][c],snow->Wsubgrid->co[r][c],i,ns,
									snow->w_liq->co[i][r][c],snow->w_ice->co[i][r][c],snow->Dzl->co[i][r][c],snow->T->co[i][r][c],P[r][c],
									snow->softSWE->co[r][c],snow->softSWE1->co[r][c],D,SWE,snow->ListonSWE->co[r][c]-DW0,DW0,
									snow->ListonSWE->co[r][c]-DW0-SWE);
								fclose(f);
								DW=snow->Wsubgrid->co[r][c];
								snow->Wsalt->co[r][c]=0.0;
								snow->Wsusp->co[r][c]=0.0;
								snow->Wsubl->co[r][c]=0.0;
							}

							rho=snow->w_ice->co[i][r][c]/(0.001*snow->Dzl->co[i][r][c]);	//rho snow [kg/m3]
							snow->w_ice->co[i][r][c]+=DW;				//kg/m2
							DW=0.0;
							snow->Dzl->co[i][r][c]+=1.0E3*(DW/rho);	//mm

							i--;

						}while(snow->w_ice->co[i+1][r][c]<0 && i>0);

						if(i==0 && snow->w_ice->co[i+1][r][c]<0){
							snow->w_ice->co[i+1][r][c]=0.0;				//kg/m2
							snow->Dzl->co[i+1][r][c]=0.0;	//mm
							snow->lnum->co[r][c]=0;

							f=fopen(files->co[ferr]+1,"a");
							fprintf(f,"\nSnow eroded more than snow present in cell r:%ld c:%ld\n",r,c);
							fprintf(f,"DW:%f DWsalt:%f DWsubl:%f DWsusp:%f DWsubgrid:%f i:%ld ns:%ld Wliq:%f Wice:%f Dlayer:%f T:%f Psnow:%f softSWE:%f D0:%f SWE0:%f SWEliston0:%f DW0:%f check:%f\n",
								DW,snow->Wsalt->co[r][c],snow->Wsubl->co[r][c],snow->Wsusp->co[r][c],snow->Wsubgrid->co[r][c],i,ns,
								snow->w_liq->co[i][r][c],snow->w_ice->co[i][r][c],snow->Dzl->co[i][r][c],snow->T->co[i][r][c],
								P[r][c],snow->softSWE->co[r][c],D,SWE,snow->ListonSWE->co[r][c]-DW0,DW0,snow->ListonSWE->co[r][c]-DW0-SWE);
							fclose(f);
						}

					}else{	//snow drifted

						snow->w_ice->co[snow->lnum->co[r][c]][r][c]+=DW;
						snow->Dzl->co[snow->lnum->co[r][c]][r][c]+=1.0E+3*DW/rho_newlyfallensnow(met->Vgrid->co[r][c], met->Tgrid->co[r][c], Tfreezing);

					}

				}else{	//snot not on the soil

					if(DW<0){

						f=fopen(files->co[ferr]+1,"a");
						fprintf(f,"Snow erosion >0 for snowD=0 in cell r:%ld c:%ld\n",r,c);
						fclose(f);

					}else{

						snow->w_ice->co[1][r][c]+=DW;
						snow->Dzl->co[1][r][c]+=1.0E+3*DW/rho_newlyfallensnow(met->Vgrid->co[r][c], met->Tgrid->co[r][c], Tfreezing);
						snow->T->co[1][r][c]=Ta[r][c];
						if(snow->T->co[1][r][c]>Tfreezing) snow->T->co[1][r][c]=Tfreezing;
					}
				}
				D=DEPTH(r,c,snow->lnum,snow->Dzl);
				if(D!=D){
					printf("Novalue in set_windtrans_snow(b): r:%ld c:%ld SnowD:%f lnum:%ld\n",r,c,DEPTH(r, c, snow->lnum, snow->Dzl),snow->lnum->co[r][c]);
					printf("Snow density :%f V:%f T:%f\n",rho_newlyfallensnow(met->Vgrid->co[r][c], met->Tgrid->co[r][c], Tfreezing),met->Vgrid->co[r][c], met->Tgrid->co[r][c]);
					write_snow_all(r, c, snow);
					printf("%f %f %f %f\n",DW,snow->Wsalt->co[r][c],snow->Wsubl->co[r][c],snow->Wsusp->co[r][c],snow->Wsubgrid->co[r][c]);
				}
				snow_layer_combination(r, c, snow, Ta[r][c], par->snowlayer_inf, par->Dmin, par->Dmax, t);
				D=DEPTH(r,c,snow->lnum,snow->Dzl);
				if(D!=D){
					printf("Novalue in set_windtrans_snow(c): r:%ld c:%ld SnowD:%f lnum:%ld\n",r,c,DEPTH(r, c, snow->lnum, snow->Dzl),snow->lnum->co[r][c]);
					write_snow_all(r, c, snow);
					printf("%f %f %f %f\n",DW,snow->Wsalt->co[r][c],snow->Wsubl->co[r][c],snow->Wsusp->co[r][c],snow->Wsubgrid->co[r][c]);
				}
			}
		}
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
void get_softsnow(SNOW *snow, double **Z){

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
}
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

	long l,mup,mdw,N=Dmin->nh;

	if(n==N){

		for(l=1;l<=n;l++){
			Dmin2->co[l]=Dmin->co[l];
			Dmax2->co[l]=Dmax->co[l];
		}

	}else{

		mup=ceil(n/2.0);
		mdw=floor(n/2.0);

		//printf("a. n=%ld linf=%ld mup=%ld mdw=%ld\n",n,linf,mup,mdw);

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
