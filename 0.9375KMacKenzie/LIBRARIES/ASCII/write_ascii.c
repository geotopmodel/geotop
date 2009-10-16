
/* STATEMENT:

ASCII-GEOtop LIBRARIES

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon

 LICENSE:

 This file is part of ASCII-GEOtop LIBRARIES
 ASCII-GEOtop LIBRARIES is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

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
	char *filename;

	filename=join_strings(name,ascii_grass);
	f=fopen(filename,"w");
	free(filename);
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
		if(r<=DTM->nrh) fprintf(f,"\n");
	}
	fclose(f);
}
