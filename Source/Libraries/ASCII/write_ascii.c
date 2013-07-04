
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
    
#include "turtle.h"
#include "write_ascii.h"
#include "extensions.h"


void write_fluidturtle(char *name, short type, DOUBLEMATRIX *DTM, T_INIT *UV){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;

	f=t_fopen(join_strings(name,ascii_ft),"w");

	fprintf(f,"/** MAP written by write_fluidturtle subroutine */\n");
	fprintf(f,"index{3}\n");
	fprintf(f,"1: double array pixels size  {%f,%f,%f,%f}\n",UV->U->co[1],UV->U->co[2],UV->U->co[3],UV->U->co[4]);
	fprintf(f,"2: double array novalues  {%ld.0,%ld.0}\n",(long)(UV->V->co[1]),(long)(UV->V->co[2]));
	fprintf(f,"3: double matrix state_variable  {%ld,%ld}\n",DTM->nrh,DTM->nch);
	for(r=1;r<=DTM->nrh;r++){
		for(c=1;c<=DTM->nch;c++){
			if(DTM->co[r][c]==UV->V->co[1]*fabs(UV->V->co[2])){
				fprintf(f,"%ld.0",(long)(DTM->co[r][c]));
			}else{
				if(type==1){
					fprintf(f,"%ld",(long)(DTM->co[r][c]));
				}else{
					fprintf(f,"%f",DTM->co[r][c]);
				}
			}
			if(c!=DTM->nch) fprintf(f," ");
		}
		if(r!=DTM->nrh) fprintf(f,"\n");
	}
	t_fclose(f);		
}

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

void write_grassascii(char *name, short type, DOUBLEMATRIX *DTM, T_INIT *UV){

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
			if(DTM->co[r][c]==UV->V->co[1]*fabs(UV->V->co[2])){
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

void write_esriascii(char *name, short type, DOUBLEMATRIX *DTM, T_INIT *UV){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;

	if(UV->U->co[1]!=UV->U->co[2]){
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		t_error("Fatal error");
	}

	f=fopen(join_strings(name,ascii_esri),"w");
	
	fprintf(f,"ncols         %ld\n",DTM->nch);
	fprintf(f,"nrows         %ld\n",DTM->nrh);
	fprintf(f,"xllcorner     %f\n",UV->U->co[4]);
	fprintf(f,"yllcorner     %f\n",UV->U->co[3]);
	fprintf(f,"cellsize      %f\n",UV->U->co[1]);
	fprintf(f,"NODATA_value  %ld.0\n",(long)(UV->V->co[1]*fabs(UV->V->co[2])));
	for(r=1;r<=DTM->nrh;r++){
		for(c=1;c<=DTM->nch;c++){
			if(DTM->co[r][c]==UV->V->co[1]*fabs(UV->V->co[2])){
				fprintf(f,"%ld.0",(long)(DTM->co[r][c]));
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
