
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1 
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "channels.h"
#include "geotop_common.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void enumerate_channels(Channel *cnet, GeoMatrix<double>& LC, GeoMatrix<short>& pixel_type, GeoMatrix<double>& Z, GeoMatrix<double>& slope, long novalue){
	
	long r, c, rnext, cnext, i=0;
	
	cnet->length.resize(cnet->length.size(),0.0);
	do{
		
		find_max_constraint( Z, LC, pixel_type, cnet->ch, novalue, &r, &c);
		
		if(r>0){
			
			i++;
			cnet->r[i]=r;
			cnet->c[i]=c;
			cnet->ch[r][c]=i;
			
			do{
				
				next_down_channel_pixel( r, c, Z, LC, pixel_type, cnet->ch, novalue, &rnext, &cnext);
				if(rnext>0){
					i++;
					if(fabs(rnext-r)==1 && fabs(cnext-c)==1){
						cnet->length[i-1]+=0.5*sqrt(2.);
						cnet->length[i]+=0.5*sqrt(2.);
					}else{
						cnet->length[i-1]+=0.5;
						cnet->length[i]+=0.5;
					}
					
					cnet->r[i]=rnext;
					cnet->c[i]=cnext;
					cnet->ch[rnext][cnext]=i;
					r=rnext;
					c=cnext;
					
				}
				
			}while(rnext>0);
			
			if(rnext<0){
				if(fabs(-rnext-r)==1 && fabs(-cnext-c)==1){
					cnet->length[i]+=0.5*sqrt(2.);
				}else{
					cnet->length[i]+=0.5;
				}
			}else{
				cnet->length[i]+=0.5;
			}		
		}
		
	}while(r>0);

    {
        size_t i; //Quiet the compiler
	
        for(i=1;i<cnet->r.size();i++){
            r = cnet->r[i];
            c = cnet->c[i];
            cnet->length[i] *= (geotop::common::Variables::UV->U[1]/cos(slope[r][c]*GTConst::Pi/180.));
        }

    }
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void next_down_channel_pixel( long r, long c, GeoMatrix<double>& Z, GeoMatrix<double>& LC, GeoMatrix<short>& pixel_type, GeoMatrix<long>& CH, long novalue, long *R, long *C){
	
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

/// ~brief find highest channel pixel that has not been enumerated yet
void find_max_constraint( GeoMatrix<double>& Z, GeoMatrix<double>& LC, GeoMatrix<short>& pixel_type, GeoMatrix<long>& CH, long novalue, long *R, long *C){
	
  size_t r, c;
	double z = -1.E99;
	
	*R=0;
	*C=0;
	
	for(r=1;r<CH.getRows();r++){
		for(c=1;c<CH.getCols();c++){
			if((long)LC[r][c]!=novalue){
				if(pixel_type[r][c]>=10 && CH[r][c]==0){
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

short neighboring_down_channel_pixel( long r, long c, long ir, long ic, GeoMatrix<double>& Z, GeoMatrix<double>& LC, GeoMatrix<short>& pixel_type, GeoMatrix<long>& CH, long novalue){
	
	short yes=0;
  size_t R=r+ir, C=c+ic;
	
	if(R>=1 && R<CH.getRows() && C>=1 && C<CH.getCols()){
		if((long)LC[R][C]!=novalue){
		//	neighboring pixel is a channel
			if(Z[R][C]<=Z[r][c] && pixel_type[R][C]>=10) yes=-1;
		//	neighboring pixel is a channels that has not been enumerated yet
			if(Z[R][C]<=Z[r][c] && pixel_type[R][C]>=10 && CH[R][C]==0) yes=1;
		}
	}
	
	return(yes);
}


