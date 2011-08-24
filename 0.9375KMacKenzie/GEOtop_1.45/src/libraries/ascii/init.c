
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

#include "../fluidturtle/turtle.h"

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initmatrix(double val, DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue){
	
	long r,c;
	for(r=1;r<=destination->nrh;r++){
		for(c=1;c<=destination->nch;c++){
			if(origin->co[r][c]!=novalue) destination->co[r][c]=val;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initlongmatrix(long val, LONGMATRIX *destination, DOUBLEMATRIX *origin, double novalue){
	
	long r,c;
	for(r=1;r<=destination->nrh;r++){
		for(c=1;c<=destination->nch;c++){
			if(origin->co[r][c]!=novalue) destination->co[r][c]=val;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void inittensor(double val, DOUBLETENSOR *destination, DOUBLEMATRIX *origin, double novalue){
	
	long l,r,c;
	for(r=1;r<=destination->nrh;r++){
		for(c=1;c<=destination->nch;c++){
			if(origin->co[r][c]!=novalue){
				for(l=1;l<=destination->ndh;l++){
					destination->co[l][r][c]=val;
				}
			}
		}
	}
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

