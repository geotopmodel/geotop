#include "constants.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "util_math.h"
#include "transport.model.h"
#include "meteodata.h"
#include <time.h>

extern long Nl, Nr, Nc; //max l, max r, max c
extern T_INIT *UV;	//dem information

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************


short transport_model(double Dt, ALLDATA *adt){
/*
	if (Dt>Dt_transport_model(Courant, M, C, H, K, D, adt){
		return 1;
	}else{
		calculate_new_concentrations(Dt, M, C, H, K, D, adt);
		return 0;
	}
*/

printf("Hallo ich bin das Transport Modell");
return 0;

}


//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

double Dt_transport_model(double Courant, DOUBLEVECTOR *M, DOUBLEVECTOR *C, DOUBLEVECTOR *H, DOUBLEVECTOR *K, DOUBLEVECTOR *D, ALLDATA *adt){
	
	long d; //direction index
	long i; //cell index
	long j; //neighbouring cell index
	long l, r, c;	//layer,row,column
	long ch;	//channel index
	long n=(Nl+1)*adt->P->total_pixel;	//number of land cells
	long cnt=0; //counter
	
	//directional vectors
	long il[7] = {0, 1,  0, 0,  0, 0, 0};
	long ir[7] = {0, 0, -1, 1,  0, 0, 0};
	long ic[7] = {0, 0,  0, 0, -1, 1, 0};
	
	double adv,diff;//advective and diffusive fluxes
	double advc;
	
	double Dt=1.E99;
	
	//i goes from 1 to n+nch, where n is the number of land cells and nch the number of channel cells
	for (i=1; i<=H->nh; i++) {
		
		if( i<=n){//land
			
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			ch=adt->C->ch->co[r][c];
			
		}else{//channel
			
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
			
		}
		
		/*we check the neighbouring cell in 5 different directions (direction up is excluded because we are looking
		 only to the lower diagonal part of the matrix*/
		
		for (d=1; d<=6; d++) {
			
			//the neighbouring cell must be within the control volume
			if (r+ir[d]>=1 && r+ir[d]<=Nr && c+ic[d]>=1 && c+ic[d]<=Nc && l+il[d]>=0 && l+il[d]<=Nl) {
				
				if(d<=5){//look for neighbouring land cell
					j = adt->T->i_cont[l+il[d]][r+ir[d]][c+ic[d]];
				}else {//look for neighbouring channel cell
					if (ch>0) {//there is actually a neighbouring channel cell
						j = n + adt->C->ch3[l][ch];
					}else {//there is NO neighbouring channel cell
						j = i;
					}
				}
				
				//we consider here only the lower diagonal
				if (j > i) {
					
					/*in the channel lateral transport in the soil is not considered
					 (i>n) it is a channel, then there is transport only in vertical direction
					 or (i<=n) it is not a channel, so there is transport in any direction (1 to 5)
					 if the pixel has a channel, we have to consider the exchange land-channel (d==6), which occurs for l>0*/
					if ( (i>n && d==1) || (i<=n && d<=5) || (i<=n && d==6 && l>0 && ch>0) ) {
						
						//counter is updated, counter is the correspondent position in the Kij matrix
						cnt++; 
						
						//if l==0 concentration is not defined, but we have to increase cnt, otherwise the matrixes are not consistent
						if(l>0){
							
							//ADVECTION
							//K is in [m2/s]
							//H is in [mm]
							//adv is in [m3/s]
							//M is in [kg]
							//C is in [kg/m3]
							adv = K->co[cnt] * 1.E-3*(H->co[i]-H->co[j]);
							
							if(adv >= 0){ //advection from j to i
								advc = Dt * (adv * C->co[j]);
							}else { //advection from i to j
								advc = Dt * (adv * C->co[i]);
							}
							
							//DIFFUSION-DISPERSION
							diff = D->co[cnt] * (C->co[j] - C->co[i]);
							
							//case 1
							//advection and diffusion in different direction
							if (adv >= 0 && diff<0){
								
								//case 1a
								if (fabs(advc) >= fabs(diff)) {
									if(Dt < Courant * M->co[j] / (fabs(advc)-fabs(diff))) Dt = Courant * M->co[j] / (fabs(advc)-fabs(diff));
								//case 1b
								}else {
									if(Dt < Courant * M->co[i] / (fabs(advc)-fabs(diff))) Dt = Courant * M->co[i] / (fabs(advc)-fabs(diff));
								}
								

						}
					}					
				}
			}
		}
	}
	
	return Dt;
}





//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

/*
 Dt = time step [s]
 M = tracer mass [kg], it is modified by this routine
 C = tracer concentration [kg/m3], input
 H = total head [mm], input
 K = matrix of (-area*k/d) [m2/s], input
 D = diffusion-dispersion matrix (area*diff/d) [m3/s], input
*/

void calculate_new_concentrations(double Dt, DOUBLEVECTOR *M, DOUBLEVECTOR *C, DOUBLEVECTOR *H, DOUBLEVECTOR *K, DOUBLEVECTOR *D, ALLDATA *adt){
	
	long d; //direction index
	long i; //cell index
	long j; //neighbouring cell index
	long l, r, c;	//layer,row,column
	long ch;	//channel index
	long n=(Nl+1)*adt->P->total_pixel;	//number of land cells
	long cnt=0; //counter
	
	//directional vectors
	long il[7] = {0, 1,  0, 0,  0, 0, 0};
	long ir[7] = {0, 0, -1, 1,  0, 0, 0};
	long ic[7] = {0, 0,  0, 0, -1, 1, 0};
	
	double adv,diff;//advective and diffusive fluxes
		
	//i goes from 1 to n+nch, where n is the number of land cells and nch the number of channel cells
	for (i=1; i<=H->nh; i++) {
				
		if( i<=n){//land
			
			l=adt->T->lrc_cont->co[i][1];
			r=adt->T->lrc_cont->co[i][2];
			c=adt->T->lrc_cont->co[i][3];
			ch=adt->C->ch->co[r][c];
						
		}else{//channel
			
			l=adt->C->lch->co[i-n][1];
			ch=adt->C->lch->co[i-n][2];
			r=adt->C->r->co[ch];
			c=adt->C->c->co[ch];
		
		}
		
		/*we check the neighbouring cell in 5 different directions (direction up is excluded because we are looking
		 only to the lower diagonal part of the matrix*/
		
		for (d=1; d<=6; d++) {
			
			//the neighbouring cell must be within the control volume
			if (r+ir[d]>=1 && r+ir[d]<=Nr && c+ic[d]>=1 && c+ic[d]<=Nc && l+il[d]>=0 && l+il[d]<=Nl) {
				
				if(d<=5){//look for neighbouring land cell
					j = adt->T->i_cont[l+il[d]][r+ir[d]][c+ic[d]];
				}else {//look for neighbouring channel cell
					if (ch>0) {//there is actually a neighbouring channel cell
						j = n + adt->C->ch3[l][ch];
					}else {//there is NO neighbouring channel cell
						j = i;
					}
				}

				//we consider here only the lower diagonal
				if (j > i) {
					
					/*in the channel lateral transport in the soil is not considered
					 (i>n) it is a channel, then there is transport only in vertical direction
					 or (i<=n) it is not a channel, so there is transport in any direction (1 to 5)
					 if the pixel has a channel, we have to consider the exchange land-channel (d==6), which occurs for l>0*/
					if ( (i>n && d==1) || (i<=n && d<=5) || (i<=n && d==6 && l>0 && ch>0) ) {
						
						//counter is updated, counter is the correspondent position in the Kij matrix
						cnt++; 
						
						//if l==0 concentration is not defined, but we have to increase cnt, otherwise the matrixes are not consistent
						if(l>0){
							
							//ADVECTION
							//K is in [m2/s]
							//H is in [mm]
							//adv is in [m3/s]
							//M is in [kg]
							//C is in [kg/m3]
							adv = K->co[cnt] * 1.E-3*(H->co[i]-H->co[j]);
							
							if(adv >= 0){ //advection from j to i
								M->co[i] += Dt * (adv * C->co[j]);
								M->co[j] -= Dt * (adv * C->co[j]);
							}else { //advection from i to j
								M->co[i] += Dt * (adv * C->co[i]);
								M->co[j] -= Dt * (adv * C->co[i]); 
							}
							
							//DIFFUSION-DISPERSION
							diff = D->co[cnt] * (C->co[j] - C->co[i]);
							
							M->co[i] += Dt * diff;
							M->co[j] -= Dt * diff;
						}
					}					
				}
			}
		}
	}
	
}


	
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

//Conversion from Mland, Mchannel -> V
void C_conversion_matrices_to_vector(DOUBLEVECTOR *V, DOUBLEMATRIX *Mland, DOUBLEMATRIX *Mchannel, ALLDATA *adt){
	
	long i; //cell index
	long j; //surface index
	long n=(Nl+1)*adt->P->total_pixel;	//number of land cells
	long l, r, c;	//layer,row,column
	
	for(i=1; i<=V->nh; i++){
		
		if (i<=n) { //land
			
			l = adt->T->lrc_cont->co[i][1];
			r = adt->T->lrc_cont->co[i][2];
			c = adt->T->lrc_cont->co[i][3];
			j = adt->T->j_cont[r][c];
			
			if (l==0) {
				V->co[i] = 0.;
			}else {
				V->co[i] = Mland->co[l][j];
			}
			
		}else {  //channel
			
			l = adt->C->lch->co[i-n][1];
			j = adt->C->lch->co[i-n][2];
			
			if (l==0) {
				V->co[i] = 0.;
			}else {
				V->co[i] = Mchannel->co[l][j];
			}
		}
	}
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

//Conversion from V -> Mland, Mchannel
void C_conversion_vector_to_matrices(DOUBLEVECTOR *V, DOUBLEMATRIX *Mland, DOUBLEMATRIX *Mchannel, ALLDATA *adt){
	
	long i; //cell index
	long j; //surface index
	long n=(Nl+1)*adt->P->total_pixel;	//number of land cells
	long l, r, c;	//layer,row,column
	
	for(i=1; i<=V->nh; i++){
		
		if (i<=n) { //land
			
			l = adt->T->lrc_cont->co[i][1];
			r = adt->T->lrc_cont->co[i][2];
			c = adt->T->lrc_cont->co[i][3];
			j = adt->T->j_cont[r][c];
			
			if (l>0) {Mland->co[l][j] = V->co[i];}

		}else {  //channel
			
			l = adt->C->lch->co[i-n][1];
			j = adt->C->lch->co[i-n][2];
			
			if (l>0) {Mchannel->co[l][j] = V->co[i];}

		}
	}
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

/*


 ^-^
(0 0)
  v
(   )
 V V
 m m


*/

