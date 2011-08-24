
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
#include "write_ascii.h"
#include "extensions.h"

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

void write_grassascii(char *name, short type, DOUBLEMATRIX *DTM, T_INIT *UV, long novalue){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;

	f=fopen(join_strings(name,ascii_grass),"w");
	
	fprintf(f,"north:%f\n",UV->U->co[3]+DTM->nrh*UV->U->co[1]);
	fprintf(f,"south:%f\n",UV->U->co[3]);
	fprintf(f,"east:%f\n",UV->U->co[4]+DTM->nch*UV->U->co[2]);
	fprintf(f,"west:%f\n",UV->U->co[4]);
	fprintf(f,"rows:%ld\n",DTM->nrh);
	fprintf(f,"cols:%ld\n",DTM->nch);
	for(r=1;r<=DTM->nrh;r++){
		for(c=1;c<=DTM->nch;c++){
			if((long)DTM->co[r][c]==(long)novalue){
				fprintf(f,"*");
			}else{
				if(type==1){
					fprintf(f,"%ld",(long)(DTM->co[r][c]));
				}else{
					fprintf(f,"%f",DTM->co[r][c]);
				}
			}
			if(c<DTM->nch) fprintf(f," ");
		}
		if(r<DTM->nrh) fprintf(f,"\n");
	}
	fclose(f);
		
}

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

void write_esriascii(char *name, short type, DOUBLEMATRIX *DTM, T_INIT *UV, long novalue){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;
	char *temp;

	if(UV->U->co[1]!=UV->U->co[2]){
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		t_error("Fatal error");
	}

	temp = join_strings(name,ascii_esri);
	f=fopen(temp,"w");
	
	fprintf(f,"ncols         %ld\n",DTM->nch);
	fprintf(f,"nrows         %ld\n",DTM->nrh);
	fprintf(f,"xllcorner     %f\n",UV->U->co[4]);
	fprintf(f,"yllcorner     %f\n",UV->U->co[3]);
	fprintf(f,"cellsize      %f\n",UV->U->co[1]);
	fprintf(f,"NODATA_value  %ld.0\n",novalue);
	for(r=1;r<=DTM->nrh;r++){
		for(c=1;c<=DTM->nch;c++){
			if((long)DTM->co[r][c]==novalue){
				fprintf(f,"%ld.0",novalue);
			}else{
				if(type==1){
					fprintf(f,"%ld",(long)(DTM->co[r][c]));
				}else{
					fprintf(f,"%f",DTM->co[r][c]);
				}
			}
			if(c<DTM->nch) fprintf(f," ");
		}
		if(r<DTM->nrh) fprintf(f,"\n");
	}
	fclose(f);
	free(temp);
}
