
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "channels.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void enumerate_channels(CHANNEL *cnet, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, DOUBLEMATRIX *Z, DOUBLEMATRIX *slope, long novalue){
	
	long r, c, rnext, cnext, i=0;
	
	initialize_doublevector(cnet->length, 0.);
	
	do{
		
		find_max_constraint( Z->co, LC, pixel_type, cnet->ch, novalue, &r, &c);
		
		if(r>0){
			
			i++;
			cnet->r->co[i]=r;
			cnet->c->co[i]=c;
			cnet->ch->co[r][c]=i;
			
			do{
				
				next_down_channel_pixel( r, c, Z->co, LC, pixel_type, cnet->ch, novalue, &rnext, &cnext);
				if(rnext>0){
					i++;
					if(fabs(rnext-r)==1 && fabs(cnext-c)==1){
						cnet->length->co[i-1]+=0.5*sqrt(2.);
						cnet->length->co[i]+=0.5*sqrt(2.);
					}else{
						cnet->length->co[i-1]+=0.5;
						cnet->length->co[i]+=0.5;				
					}
					
					cnet->r->co[i]=rnext;
					cnet->c->co[i]=cnext;
					cnet->ch->co[rnext][cnext]=i;
					r=rnext;
					c=cnext;
					
				}
				
			}while(rnext>0);
			
			if(rnext<0){
				if(fabs(-rnext-r)==1 && fabs(-cnext-c)==1){
					cnet->length->co[i]+=0.5*sqrt(2.);
				}else{
					cnet->length->co[i]+=0.5;				
				}
			}else{
				cnet->length->co[i]+=0.5;				
			}		
		}
		
	}while(r>0);
	
	for(i=1;i<=cnet->r->nh;i++){
		r = cnet->r->co[i];
		c = cnet->c->co[i];
		cnet->length->co[i] *= (UV->U->co[1]/cos(slope->co[r][c]*Pi/180.));
	}
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void next_down_channel_pixel( long r, long c, double **Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long novalue, long *R, long *C){
	
	*R=0;
	*C=0;
	
	if(neighboring_down_channel_pixel(r, c, 1, 1, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r+1;
		*C=c+1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, -1, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r+1;
		*C=c-1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 1, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r-1;
		*C=c+1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}																																																														
	
	if(neighboring_down_channel_pixel(r, c, -1, -1, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r-1;
		*C=c-1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, 0, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r+1;
		*C=c;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, 1, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r;
		*C=c+1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 0, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r-1;
		*C=c;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, -1, Z, LC, pixel_type, CH, novalue) == -1){
		*R=r;
		*C=c-1;
		*R=(*R)*(-1.);
		*C=(*C)*(-1.);
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, 0, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r+1;
		*C=c;
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, 1, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r;
		*C=c+1;
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 0, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r-1;
		*C=c;
	}
	
	if(neighboring_down_channel_pixel(r, c, 0, -1, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r;
		*C=c-1;
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, 1, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r+1;
		*C=c+1;
	}
	
	if(neighboring_down_channel_pixel(r, c, 1, -1, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r+1;
		*C=c-1;
	}
	
	if(neighboring_down_channel_pixel(r, c, -1, 1, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r-1;
		*C=c+1;
	}																																																														
	
	if(neighboring_down_channel_pixel(r, c, -1, -1, Z, LC, pixel_type, CH, novalue) == 1){
		*R=r-1;
		*C=c-1;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//find highest channel pixel that has not been enumerated yet	
void find_max_constraint( double **Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long novalue, long *R, long *C){
	
	long r, c;
	double z = -1.E99;
	
	*R=0;
	*C=0;
	
	for(r=1;r<=CH->nrh;r++){
		for(c=1;c<=CH->nch;c++){
			if((long)LC->co[r][c]!=novalue){
				if(pixel_type->co[r][c]>=10 && CH->co[r][c]==0){
					if(Z[r][c]>z){
						z=Z[r][c];
						*R=r;
						*C=c;
					}
				}
			}
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short neighboring_down_channel_pixel( long r, long c, long ir, long ic, double **Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long novalue){
	
	short yes=0;
	long R=r+ir, C=c+ic;
	
	if(R>=1 && R<=CH->nrh && C>=1 && C<=CH->nch){
		if((long)LC->co[R][C]!=novalue){
			//neighboring pixel is a channel
			if(Z[R][C]<=Z[r][c] && pixel_type->co[R][C]>=10) yes=-1;
			//neighboring pixel is a channels that has not been enumerated yet
			if(Z[R][C]<=Z[r][c] && pixel_type->co[R][C]>=10 && CH->co[R][C]==0) yes=1;
		}
	}
	
	return(yes);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
