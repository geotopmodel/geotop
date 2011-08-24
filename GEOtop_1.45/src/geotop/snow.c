
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
//Date: 13 November 2005
//Contents: Snow subroutines

#include "constants.h"
#include "struct.geotop.h"
#include "snow.h"
#include "../libraries/ascii/extensions.h"
#include "../libraries/ascii/rw_maps.h"
#include "times.h"
#include "output.h"

#include "PBSM.h"
//#include "SnowTran.h"

extern long number_novalue, number_absent;

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern char **files;
extern long Nl, Nr, Nc;
extern long i_sim;


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//Jordan et al., 1999
double rho_newlyfallensnow(double u, double Tatm, double Tfreez){

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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snow_compactation(double Dt, long r, long c, long l, STATEVAR_3D *snow, double slope, PAR *par){

	long m;
	double theta_i,theta_w,c1,c2,c3,c4,c5,load,eta0,eta,CR1,CR2;

	theta_i=snow->w_ice->co[l][r][c]/(0.001*snow->Dzl->co[l][r][c]*rho_i);	
	theta_w=snow->w_liq->co[l][r][c]/(0.001*snow->Dzl->co[l][r][c]*rho_w);	

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
		CR1=-c1*c2*exp(-c3*(Tfreezing-snow->T->co[l][r][c]));
	
		//OVERBURDEN
		eta0=par->snow_viscosity;  //[kg s m^-2]
		c4=0.08;	 //[K^-1]
		c5=0.023;    //[m^3 kg^-1]
		load=0.0;    //[kg m^-2]
		for(m=l;m<=snow->lnum->co[r][c];m++){
			load+=(snow->w_ice->co[m][r][c]+snow->w_liq->co[m][r][c]);
		}
		load*=fabs(cos(slope*Pi/180.));
		eta=eta0*exp(c4*(Tfreezing-snow->T->co[l][r][c])+c5*(rho_i*theta_i));
		CR2=-load/eta;
	
		snow->Dzl->co[l][r][c] *= exp( (CR1 + CR2) * Dt ); // Cuell

	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snow_layer_combination(double a, long r, long c, STATEVAR_3D *snow, double Ta, long linf, DOUBLEVECTOR *Dmin, DOUBLEVECTOR *Dmax, double SWEmax, double time, FILE *flog)

{
	
	long l, cont, max=Dmin->nh;
	short occuring, sux;
	double D0, D1, D, SWE, ice;
	DOUBLEVECTOR *Dmin2, *Dmax2;
	
	//check on SWEmax[kg/m2]
	SWE = 0.0;
	l=1;
	while (l <= snow->lnum->co[r][c]){
		SWE += snow->w_ice->co[l][r][c];
		ice = snow->w_ice->co[l][r][c];
		snow->w_ice->co[l][r][c] -= Fmax(0.0, SWE-SWEmax);
		if (snow->w_ice->co[l][r][c] < 0) snow->w_ice->co[l][r][c] = 0.0;
		if (ice > 0){
			snow->Dzl->co[l][r][c] *= snow->w_ice->co[l][r][c]/ice;
			snow->w_liq->co[l][r][c] *= snow->w_ice->co[l][r][c]/ice;
		}
		l++;
	}
	
	//D=snow depth(mm)
	D = 0.0;
	for(l=1;l<=max;l++){
		D += snow->Dzl->co[l][r][c];
	}
	
	//PREPROCESSING
	//1. If the snow depth is too small, then it is reset to 0
	if(D<=1.E-6){
		snow->lnum->co[r][c]=0;
		snow->type->co[r][c]=0;
		for(l=1;l<=max;l++){
			snow->Dzl->co[l][r][c]=0.0;
			snow->w_ice->co[l][r][c]=0.0;
			snow->w_liq->co[l][r][c]=0.0;
			snow->T->co[l][r][c]=0.0;	//Temperatura di inizializzazione
		}
		
		
	//2. If D<Dmin(max), we are in the simplified case
	}else if(snow->lnum->co[r][c]>0 && D<Dmin->co[max]){
		
		snow->type->co[r][c]=1;
		if(snow->lnum->co[r][c]>1){
			for(l=snow->lnum->co[r][c];l>1;l--){
				snowlayer_merging(a, r, c, snow, l, l-1, l-1);
			}
			for(l=2;l<=max;l++){
				snow->T->co[l][r][c]=0.0;
				snow->Dzl->co[l][r][c]=0.0;
				snow->w_liq->co[l][r][c]=0.0;
				snow->w_ice->co[l][r][c]=0.0;
			}
			snow->lnum->co[r][c]=1;
		}

	//3. if z>=Dmin(max), ordinary case
	}else if(snow->lnum->co[r][c]>0 && D>=Dmin->co[max]){
		
		snow->type->co[r][c]=2;	
			
	//4. if there is not yet a snow layer and 1mm<=D<Dmin(max), simplified case
	}else if(snow->lnum->co[r][c]==0 && D<Dmin->co[max]){
		
		snow->lnum->co[r][c]=1;
		snow->type->co[r][c]=1;
		snow->T->co[1][r][c]=Fmin(Ta,-0.1);

	
	//5. if there is not yet a snow layer and D>=Dmin(max), simplified case
	}else if(snow->lnum->co[r][c]==0 && D>=Dmin->co[max]){
		
		snow->lnum->co[r][c]=1;
		snow->type->co[r][c]=2;
		snow->T->co[1][r][c]=Fmin(Ta,-0.1);

	}
	
	// SIMMETRICAL PARAMETERIZATION SCHEME (new)
			
	if(snow->type->co[r][c]==2){
		
		Dmax->co[linf]=1.E10;	//snow layer of unlimited thickness 
				
		Dmin2=new_doublevector(max);
		Dmax2=new_doublevector(max);
		
		min_max_layer(snow->lnum->co[r][c], Dmin, Dmax, Dmin2, Dmax2, linf);
		
		//show_Dminmax(r, c, Dmin2->co, Dmax2->co, snow->lnum->co[r][c]);
		//write_snow_all(r, c, snow);	
		
		D0=0.0;
		for(l=1;l<=max;l++){
			D0+=snow->Dzl->co[l][r][c];
		}		
				
		cont=0;
		
		do{		
			
			cont++;
			
			if(cont>10){
				
				for(l=1;l<=snow->lnum->co[r][c];l++){
					if(snow->Dzl->co[l][r][c] < Dmin2->co[l]){
						printf("Iteration to assign new thicknesses to the snow layers does not converge, r:%ld c:%ld\n",r,c);
						printf("Forced change of MINIMUM thickness: ns:%ld l:%ld Dmin_old:%f Dmin_new:%f\n",snow->lnum->co[r][c],l,Dmin2->co[l],snow->Dzl->co[l][r][c]);
						fprintf(flog,"Iteration to assign new thicknesses to the snow layers does not converge, r:%ld c:%ld\n",r,c);
						fprintf(flog,"Forced change of MINIMUM thickness: ns:%ld l:%ld Dmin_old:%f Dmin_new:%f\n",snow->lnum->co[r][c],l,Dmin2->co[l],snow->Dzl->co[l][r][c]);
						Dmin2->co[l]=snow->Dzl->co[l][r][c];
					}
					if(snow->Dzl->co[l][r][c] > Dmax2->co[l]){
						printf("Iteration to assign new thicknesses to the snow layers does not converge, r:%ld c:%ld\n",r,c);
						printf("Forced change of MAXIMUM thickness: ns:%ld l:%ld Dmax_old:%f Dmax_new:%f\n",snow->lnum->co[r][c],l,Dmax2->co[l],snow->Dzl->co[l][r][c]);
						fprintf(flog,"Iteration to assign new thicknesses to the snow layers does not converge, r:%ld c:%ld\n",r,c);
						fprintf(flog,"Forced change of MAXIMUM thickness: ns:%ld l:%ld Dmax_old:%f Dmax_new:%f\n",snow->lnum->co[r][c],l,Dmax2->co[l],snow->Dzl->co[l][r][c]);
						Dmax2->co[l]=snow->Dzl->co[l][r][c];
					}				
				}
			}
						
			//trying to adjust the layer thicknesses maintaining the same layer number
			for(l=snow->lnum->co[r][c];l>=linf;l--){
				if (l>1) {
					sux = set_snow_max(a, r, c, snow, l, l-1, Dmax2->co[l]);
					sux = set_snow_min(a, r, c, snow, l, l-1, Dmin2->co[l]);
				}
			}
			
			for(l=1;l<=linf;l++){
				if (l<snow->lnum->co[r][c]) {
					sux = set_snow_max(a, r, c, snow, l, l+1, Dmax2->co[l]);
					sux = set_snow_min(a, r, c, snow, l, l+1, Dmin2->co[l]);
				}
			}		

			//checking if it is ok now
			occuring=0;
							
			for(l=1;l<=snow->lnum->co[r][c];l++){	
								
				if(Dmin2->co[l] - snow->Dzl->co[l][r][c] > 1.E-3){
					merge_layers(a, r, c, snow, l);	
					min_max_layer(snow->lnum->co[r][c], Dmin, Dmax, Dmin2, Dmax2, linf);
					occuring = 1;
				}
				
				if(Dmax2->co[l] - snow->Dzl->co[l][r][c] < -1.E-3){
					if(snow->lnum->co[r][c] < max){
						split_layers(r, c, snow, l);
						min_max_layer(snow->lnum->co[r][c], Dmin, Dmax, Dmin2, Dmax2, linf);
					}
					occuring = 1;
				}
			}
			
		}while(occuring==1);
				
		D1=0.0;
		for(l=1;l<=max;l++){
			D1+=snow->Dzl->co[l][r][c];
		}
		
		if(fabs(D0-D1)>0.001){
			printf("r:%ld c:%ld Dold:%f Dnew:%f\n",r,c,D0,D1);
			t_error("Error 2 in snow combination - Change values of blocks 5 and 6 parameter file");
		}
		
		//show_Dminmax(r, c, Dmin2->co, Dmax2->co, snow->lnum->co[r][c]);
		//write_snow_all(r, c, snow);	
				
		free_doublevector(Dmin2);
		free_doublevector(Dmax2);
	
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

void snowlayer_merging(double a, long r, long c, STATEVAR_3D *snow, long l1, long l2, long lres){

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
	from_internal_energy(a, r+6000,c+6000,h,&(snow->w_ice->co[lres][r][c]),&(snow->w_liq->co[lres][r][c]),&(snow->T->co[lres][r][c]));
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
		
	return (c_ice*w_ice+c_liq*w_liq)*(T-Tfreezing) + Lf*w_liq;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void from_internal_energy(double a, long r, long c, double h, double *w_ice, double *w_liq, double *T){

	double SWE=(*w_ice)+(*w_liq);
	double T0;
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
	
		*w_liq=theta_snow(a, 1., *T)*SWE;
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

void write_snow(long r, long c, long l, STATEVAR_3D *snow){

	printf("r:%ld c:%ld wice(%ld/%ld):%f wliq(%ld):%f T(%ld):%f Dz(%ld):%f\n",r,c,l,snow->lnum->co[r][c],snow->w_ice->co[l][r][c],l,snow->w_liq->co[l][r][c],
			l,snow->T->co[l][r][c],l,snow->Dzl->co[l][r][c]);
		
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_snow_all(long r, long c, STATEVAR_3D *snow){
	long l;
	for(l=1;l<=snow->lnum->co[r][c];l++){
		write_snow(r,c,l,snow);
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short set_snow_min(double a, long r, long c, STATEVAR_3D *snow, long l1, long l2, double Dmin){
	
	double h, f, dwl, dwi, dd;
	
	//printf("::l1:%ld Dz:%f min:%f \n",l1,snow->Dzl->co[l1][r][c],Dmin);
	if(snow->Dzl->co[l1][r][c] < Dmin && snow->Dzl->co[l2][r][c] > 1.E-5){	//l1 too shallow and takes mass from l2
		//printf("min.l1:%ld Dz:%f wice:%F wliq:%f l1:%ld Dz:%f wice:%F wliq:%f\n",l1,snow->Dzl->co[l1][r][c],snow->w_ice->co[l1][r][c],snow->w_liq->co[l1][r][c],l2,snow->Dzl->co[l2][r][c],snow->w_ice->co[l2][r][c],snow->w_liq->co[l2][r][c]);
		f = Fmin(Dmin - snow->Dzl->co[l1][r][c], snow->Dzl->co[l2][r][c])/snow->Dzl->co[l2][r][c];
		h = internal_energy(snow->w_ice->co[l1][r][c], snow->w_liq->co[l1][r][c], snow->T->co[l1][r][c]);
		dd = f*snow->Dzl->co[l2][r][c];
		snow->Dzl->co[l1][r][c] += dd;
		snow->Dzl->co[l2][r][c] -= dd;
		dwl = f*snow->w_liq->co[l2][r][c];
		snow->w_liq->co[l1][r][c] += dwl;
		snow->w_liq->co[l2][r][c] -= dwl;
		dwi = f*snow->w_ice->co[l2][r][c];
		snow->w_ice->co[l1][r][c] += dwi;
		snow->w_ice->co[l2][r][c] -= dwi;
		h += internal_energy(dwi, dwl, snow->T->co[l2][r][c]);
		from_internal_energy(a, r+1000, c+1000, h, &(snow->w_ice->co[l1][r][c]), &(snow->w_liq->co[l1][r][c]), &(snow->T->co[l1][r][c]));
		//printf(".l1:%ld Dz:%f wice:%F wliq:%f l1:%ld Dz:%f wice:%F wliq:%f\n",l1,snow->Dzl->co[l1][r][c],snow->w_ice->co[l1][r][c],snow->w_liq->co[l1][r][c],l2,snow->Dzl->co[l2][r][c],snow->w_ice->co[l2][r][c],snow->w_liq->co[l2][r][c]);
		return 1;
	}else {
		return 0;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
	
short set_snow_max(double a, long r, long c, STATEVAR_3D *snow, long l1, long l2, double Dmax){
	
	double h, f, dwl, dwi, dd;
		
	//printf("::l1:%ld Dz:%f max:%f\n",l1,snow->Dzl->co[l1][r][c],Dmax);
	if(snow->Dzl->co[l1][r][c] > Dmax){	//l1 too thick and gives mass to l2
		//printf("max.l1:%ld Dz:%f wice:%F wliq:%f l1:%ld Dz:%f wice:%F wliq:%f\n",l1,snow->Dzl->co[l1][r][c],snow->w_ice->co[l1][r][c],snow->w_liq->co[l1][r][c],l2,snow->Dzl->co[l2][r][c],snow->w_ice->co[l2][r][c],snow->w_liq->co[l2][r][c]);
		f = (snow->Dzl->co[l1][r][c] - Dmax)/snow->Dzl->co[l1][r][c];
		h = internal_energy(snow->w_ice->co[l2][r][c], snow->w_liq->co[l2][r][c], snow->T->co[l2][r][c]);
		dd = f*snow->Dzl->co[l1][r][c];
		snow->Dzl->co[l1][r][c] -= dd;
		snow->Dzl->co[l2][r][c] += dd;				
		dwl = f*snow->w_liq->co[l1][r][c];
		snow->w_liq->co[l1][r][c] -= dwl;
		snow->w_liq->co[l2][r][c] += dwl;
		dwi = f*snow->w_ice->co[l1][r][c];
		snow->w_ice->co[l1][r][c] -= dwi;
		snow->w_ice->co[l2][r][c] += dwi;
		h += internal_energy(dwi, dwl, snow->T->co[l1][r][c]);
		from_internal_energy(a, r+2000, c+2000, h, &(snow->w_ice->co[l2][r][c]), &(snow->w_liq->co[l2][r][c]), &(snow->T->co[l2][r][c]));
		//printf(".l1:%ld Dz:%f wice:%F wliq:%f l1:%ld Dz:%f wice:%F wliq:%f\n",l1,snow->Dzl->co[l1][r][c],snow->w_ice->co[l1][r][c],snow->w_liq->co[l1][r][c],l2,snow->Dzl->co[l2][r][c],snow->w_ice->co[l2][r][c],snow->w_liq->co[l2][r][c]);
		return 1;
	}else {
		return 0;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short set_snowice_min(double a, long r, long c, STATEVAR_1D *snow, long l1, long l2, double wicemin){
	
	//double h;
	double f, dwl, dwi, dd;
	
	if(snow->w_ice->co[l1] < wicemin && snow->w_ice->co[l2] > 1.E-6){	//l1 too shallow and takes mass from l2
		f = Fmin(wicemin - snow->w_ice->co[l1], snow->w_ice->co[l2])/snow->w_ice->co[l2];
		//h = internal_energy(snow->w_ice->co[l1], snow->w_liq->co[l1], snow->T->co[l1]);
		dd = f*snow->Dzl->co[l2];
		snow->Dzl->co[l1] += dd;
		snow->Dzl->co[l2] -= dd;
		dwl = f*snow->w_liq->co[l2];
		snow->w_liq->co[l1] += dwl;
		snow->w_liq->co[l2] -= dwl;
		dwi = f*snow->w_ice->co[l2];
		snow->w_ice->co[l1] += dwi;
		snow->w_ice->co[l2] -= dwi;
		//h += internal_energy(dwi, dwl, snow->T->co[l2]);
		//from_internal_energy(a, r+1000, c+1000, h, &(snow->w_ice->co[l1]), &(snow->w_liq->co[l1]), &(snow->T->co[l1]));
		return 1;
	}else {
		return 0;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void split_layers(long r, long c, STATEVAR_3D *snow, long l1){

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

void merge_layers(double a, long r, long c, STATEVAR_3D *snow, long l1){

	long l;
	
	if(l1>snow->lnum->co[r][c]) t_error("Error 1 in merge_layers");
	
	if(l1==snow->lnum->co[r][c]){
		snowlayer_merging(a, r, c, snow, l1, l1-1, l1-1);
		initialize_snow(r, c, l1, snow);
	}else{
		snowlayer_merging(a, r, c, snow, l1, l1+1, l1);
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

void initialize_snow(long r, long c, long l, STATEVAR_3D *snow){
	
	snow->w_ice->co[l][r][c]=0.0;
	snow->w_liq->co[l][r][c]=0.0;
	snow->Dzl->co[l][r][c]=0.0;
	snow->T->co[l][r][c]=0.0;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
	
void show_Dminmax(long r, long c, double *Dmin, double *Dmax, long n){

	long l;
	
	printf("n:%ld\n",n);
	
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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

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

	double r1, r2, r3;
					
	r1=exp(5000.0*(1.0/273.16-1.0/(Ta+273.16)));
	r2=pow(r1,10);
	if(r2>1.0) r2=1.0;
	r3=0.3;
	
	*snowage*=((r1+r2+r3)*1.0E-6);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void glac2snow(double a, long r, long c, STATEVAR_3D *snow, STATEVAR_3D *glac){

	double h;

	//if ice is less thick than Dmin, ice is considered as snow
	if(glac->type->co[r][c] == 1){
		
		h = internal_energy(glac->w_ice->co[1][r][c], glac->w_liq->co[1][r][c], glac->T->co[1][r][c]);
		h += internal_energy(snow->w_ice->co[1][r][c], snow->w_liq->co[1][r][c], snow->T->co[1][r][c]);		

		snow->Dzl->co[1][r][c] += glac->Dzl->co[1][r][c];
		snow->w_ice->co[1][r][c] += glac->w_ice->co[1][r][c];
		snow->w_liq->co[1][r][c] += glac->w_liq->co[1][r][c];
		from_internal_energy(a, r, c, h, &(snow->w_ice->co[1][r][c]), &(snow->w_ice->co[1][r][c]), &(snow->T->co[1][r][c]));
		
		glac->Dzl->co[1][r][c]=0.0;
		glac->w_liq->co[1][r][c]=0.0;
		glac->w_ice->co[1][r][c]=0.0;
		glac->T->co[1][r][c]=0.0;
		glac->lnum->co[r][c]=0;
		glac->type->co[r][c]=0;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void WBsnow(double Dt, long r, long c, STATEVAR_3D *snow, double *Melt, double *RainOnSnow, PAR *par, double slope, double Rain, double *wi, double *wl, double *T, double Edt){

	long l, m;
	double Se, thliq, thice, Dwl, Dwi, wice_old;
	
	Dwl = Rain;
	

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

			wl[m]+=Dwl;	
			wi[m]+=Dwi;	
						
			//assign snow state variables
			wice_old=snow->w_ice->co[l][r][c];
			snow->T->co[l][r][c]=T[m];
			snow->w_ice->co[l][r][c]=wi[m];
			snow->w_liq->co[l][r][c]=wl[m];
			
			//ACCOUNT FOR SNOW COMPACTION
			//a)destructive metamorphism and overburden
			snow_compactation(Dt, r, c, l, snow, slope, par);		
						
			//b)melting: snow depth decreases maintaining the same density
			if(wi[m]/wice_old<1) snow->Dzl->co[l][r][c]*=(wi[m]/wice_old);
					
			//check that ice volumetric content is not too large
			thice=snow->w_ice->co[l][r][c]/(1.E-3*snow->Dzl->co[l][r][c]*rho_i);		//[-]				
			if(thice>par->snow_maxpor){
				thice=par->snow_maxpor;
				snow->Dzl->co[l][r][c]=1.E3*snow->w_ice->co[l][r][c]/(rho_i*thice);
			}
				
			//CALCULATE LIQUID WATER GOING BELOW
			if(snow->w_ice->co[l][r][c]<0.001){
				Dwl=snow->w_liq->co[l][r][c];
				snow->w_liq->co[l][r][c]=0.0;
			}else{
				thliq=snow->w_liq->co[l][r][c]/(1.0E-3*snow->Dzl->co[l][r][c]*rho_w);		//[-]
				Se=(thliq - par->Sr*(1.0-thice))/( (1.0-thice) - par->Sr*(1.0-thice));
				if(Se<0) Se=0.0;
				if(Se>1) Se=1.0;
				if(thliq>par->Sr*(1.0-thice)){
					Dwl=Fmin(5.0*pow(Se,3.0)*Dt , (thliq - par->Sr*(1.0-thice)))*snow->Dzl->co[l][r][c]*1.E-3*rho_w;
				}else{
					Dwl=0.0;
				}
				snow->w_liq->co[l][r][c]-=Dwl;
			}
						
		}
		
		//melting rate
		*Melt = Dwl - Rain;			//[mm]
		
		//rain on snow
		*RainOnSnow = Rain;			//[mm]				
		
	//b. simplified case
	}else if(snow->type->co[r][c]==1){

		*Melt = wl[1];					//[mm]
		*RainOnSnow = 0.0;

		snow->T->co[1][r][c]=T[1];
		snow->w_liq->co[1][r][c]=0.0;
		snow->Dzl->co[1][r][c]*=(wi[1]/snow->w_ice->co[1][r][c]);
		snow->w_ice->co[1][r][c]=wi[1];
		
	}else {
		
		*Melt = 0.0;
		*RainOnSnow = 0.0;

	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void new_snow(double a, long r, long c, STATEVAR_3D *snow, double P, double Dz, double T){

	long ns;
	double h;
	
											
	if(snow->type->co[r][c]==0){

		snow->Dzl->co[1][r][c]+=Dz;
		snow->w_ice->co[1][r][c]+=P;

	}else{
					
		ns=snow->lnum->co[r][c];
			
		h=internal_energy(snow->w_ice->co[ns][r][c], snow->w_liq->co[ns][r][c], snow->T->co[ns][r][c]);
		h+=(c_ice*P)*(Fmin(T, -0.1) - Tfreezing);

		snow->Dzl->co[ns][r][c]+=Dz;
		snow->w_ice->co[ns][r][c]+=P;

		from_internal_energy(a, r, c, h, &(snow->w_ice->co[ns][r][c]),&(snow->w_liq->co[ns][r][c]),&(snow->T->co[ns][r][c]));
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void WBglacier(long ns, long r, long c, STATEVAR_3D *glac, double *Melt, PAR *par, double *wi, double *wl, double *T, double Edt){

	long l, m;
	double thice, thliq, Dwi;
								
	*Melt = 0.0;
	
	for(l=glac->lnum->co[r][c];l>=1;l--){

		m = ns + glac->lnum->co[r][c] - l + 1;
		
		//SUBTRACT SUBLIMATION in the upper layer
		if(m==1){
			Dwi=Fmax(0.0, wi[m]-Edt) - wi[m];
		}else{
			Dwi=0.0;
		}
																									
		wi[m]+=Dwi;	
							
		if(wi[m]/glac->w_ice->co[l][r][c]<1) glac->Dzl->co[l][r][c]*=(wi[m]/glac->w_ice->co[l][r][c]);
		glac->T->co[l][r][c]=T[m];
		glac->w_ice->co[l][r][c]=wi[m];
		glac->w_liq->co[l][r][c]=wl[m];

		//check that ice volumetric content is not too large
		thice=glac->w_ice->co[l][r][c]/(1.E-3*glac->Dzl->co[l][r][c]*rho_i);		//[-]				
		if(thice>0.95){
			thice=0.95;
			glac->Dzl->co[l][r][c]=1.E3*glac->w_ice->co[l][r][c]/(rho_i*thice);
		}
		
		//water flow
		thliq=glac->w_liq->co[l][r][c]/(1.E-3*glac->Dzl->co[l][r][c]*rho_w);		//[-]
		if(thliq>par->Sr_glac*(1.0-thice)){
			*Melt = *Melt + rho_w*1.E-3*glac->Dzl->co[l][r][c]*(thliq-par->Sr_glac*(1.0-thice));	//[mm]
			glac->w_liq->co[l][r][c] = rho_w*1.E-3*glac->Dzl->co[l][r][c]*par->Sr_glac*(1.0-thice);
		}
	}
	
}

		
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_SCA(STATEVAR_3D *snow, PAR *par, double **Z, double t){
	
	long l, r, c, cont=0, conttot=0;
	double T, D, SWE, Tmean=0.0, Tsmean=0.0, Dmean=0.0, SWEmean=0.0, SCA;
	
	double JD, JDfrom0;
	long day, month, year, hour, minute;
	
	FILE *f;
	
	JDfrom0 = convert_tfromstart_JDfrom0(t, par->init_date->co[i_sim]);
	convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
	convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute); 

	for(r=1;r<=Nr;r++){
		for(c=1;c<=Nc;c++){
			if((long)Z[r][c]!=number_novalue){	
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
	
	f=fopen(join_strings(files[fSCA],textfile),"a");
	fprintf(f,"%ld/%ld/%ld %ld:%02.0f",day,month,year,hour,(float)minute);
	fprintf(f,",%f,%f,%f",JDfrom0-par->init_date->co[i_sim],JDfrom0,JD);  
	fprintf(f,",%f,%f,%f,%f,%f,%f\n",Dmean,SWEmean,Tmean,Tsmean,(1.0-SCA)*100.0,SCA*100.0);
	fclose(f);
}
			

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double max_dtheta_snow(double a, double b){
	
	double maxT;
	
	maxT=-pow(b/(3*a*a), 0.5);
	return(maxT);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double theta_snow(double a, double b, double T){
	
	double th;
	
	th=1./(b+pow(a*T, 2.));
	if(T>0) th=1./(b+pow(a*0., 2.));
	return(th);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double dtheta_snow(double a, double b, double T){
	
	double dth;
	
	dth=-2.*a*a*T/pow(b+pow(a*T, 2.), 2.);
	if(T>0) dth=-2.*a*a*max_dtheta_snow(a,b)/pow(b+pow(a*max_dtheta_snow(a,b), 2.), 2.);
	return(dth);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void allocate_and_initialize_statevar_3D(STATEVAR_3D *V, double nan, long nl, long nr, long nc){
	
	V->type = new_shortmatrix(nr, nc);
	initialize_shortmatrix(V->type, 0);
	V->lnum = new_longmatrix(nr, nc);
	initialize_longmatrix(V->lnum, 0);
	V->Dzl = new_doubletensor(nl, nr, nc);
	initialize_doubletensor(V->Dzl, 0.);
	V->T = new_doubletensor(nl, nr, nc);
	initialize_doubletensor(V->T, nan);
	V->w_ice = new_doubletensor(nl, nr, nc);
	initialize_doubletensor(V->w_ice, 0.);
	V->w_liq = new_doubletensor(nl, nr, nc);
	initialize_doubletensor(V->w_liq, 0.);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_statevar_3D(STATEVAR_3D *V){
	
	free_shortmatrix(V->type);
	free_longmatrix(V->lnum);
	free_doubletensor(V->Dzl);
	free_doubletensor(V->T);
	free_doubletensor(V->w_ice);
	free_doubletensor(V->w_liq);
	free(V);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void allocate_and_initialize_statevar_1D(STATEVAR_1D *V, double nan, long nl){
	
	V->Dzl = new_doublevector(nl);
	initialize_doublevector(V->Dzl, 0.);
	V->T = new_doublevector(nl);
	initialize_doublevector(V->T, nan);
	V->w_ice = new_doublevector(nl);
	initialize_doublevector(V->w_ice, 0.);
	V->w_liq = new_doublevector(nl);
	initialize_doublevector(V->w_liq, 0.);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_statevar_1D(STATEVAR_1D *V){
	
	free_doublevector(V->Dzl);
	free_doublevector(V->T);
	free_doublevector(V->w_ice);
	free_doublevector(V->w_liq);
	free(V);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short copy_statevar_from3D_to1D(long r, long c, STATEVAR_3D *origin, STATEVAR_1D *destination){
	
	long nl, l;
	
	nl = origin->Dzl->ndh;

	if(r<1 || r>origin->type->nrh) return 0;
	if(c<1 || c>origin->type->nch) return 0;
	
	if(nl != destination->Dzl->nh) return 0;
	
	destination->type = origin->type->co[r][c];
	destination->lnum = origin->lnum->co[r][c];
	for (l=1; l<=nl; l++) {
		destination->Dzl->co[l] = origin->Dzl->co[l][r][c];
		destination->T->co[l] = origin->T->co[l][r][c];
		destination->w_ice->co[l] = origin->w_ice->co[l][r][c];
		destination->w_liq->co[l] = origin->w_liq->co[l][r][c];
	}
	
	return 1;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


