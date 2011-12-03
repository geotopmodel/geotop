
/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file ncgt_utilities.c

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico and Riccardo Rigon

This file is part of numerioc_solver.
	 Turtle_NetCdf is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

     Turtle_NetCdf is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 *
 * \file t_nc_utilities.c
 * \data 15 September 2009
 * \author Emanuele Cordano
 *
 */




#ifdef USE_NETCDF_ONGOING

#include "../libraries/fluidturtle/turtle.h"
//#include <netcdf.h>
//#include "gt_utilities.h"
#include "ncgt_utilities.h"

char *copy_stringnames(const char *origin){
/*!
 * \author Emanuele Cordano
 * \date April 2009
 *
 *
 */
char *dest;

//dest=origin;
dest=(char *)malloc((size_t)(strlen(origin)+4)*sizeof(char));
strcpy(dest,origin);
return dest;
}

//061009_s
//http://cf-pcmdi.llnl.gov/documents/cf-standard-names/standard-name-table/12/cf-standard-name-table.html
//void nc_add_cf_convention_header(const char *filename,const char *title, const char *institution,const char *source,const char *history,const char *references,const char *comment,const char *conventions){



int rotate180_y_doublematrix(DOUBLEMATRIX *M) {
/*
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLEMATRIX *M) matrix to rotate
 */
	long r,c;
	double val=0.0;


	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int rotate180_y_doubletensor(DOUBLETENSOR *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLETENSOR *M) matrix to rotate
 */
	long r,c,l;
	double val=0.0;

	long r_max=(M->nrh-M->nrl+1)/2;
//	printf("r_max=%ld %ld %ld \n",r_max,M->nrh,M->nrl);

	for(l=M->ndl;l<=M->ndh;l++) {
		for(r=M->nrl;r<=r_max;r++) {
			for (c=M->ncl;c<=M->nch;c++) {
				val=M->element[l][r][c];
				M->element[l][r][c]=M->element[l][M->nrh-r+1][c];
				M->element[l][M->nrh-r+1][c]=val;
			}
		}
	}
	return 0;
}

int rotate180_y_floatmatrix(FLOATMATRIX *M) {
/*
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLEMATRIX *M) matrix to rotate
 */
	long r,c;
	float val=0.0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}



int rotate180_y_longmatrix(LONGMATRIX *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (LONGMATRIX *M) matrix to rotate
 */
	long r,c;
	long val=0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int rotate180_y_intmatrix(INTMATRIX *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (INTMATRIX *M) matrix to rotate
 */
	long r,c;
	int val=0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int rotate180_y_shortmatrix(SHORTMATRIX *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (SHORTMATRIX *M) matrix to rotate
 */
	long r,c;
	short val=0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int invert_order_doublevector(DOUBLEVECTOR *v) {
/*!
 *
 * \author Emanuele Cordano
 * \date October 2009
 *
 *\param (DOUBLEVECTOR *) the vector which is inverted
 */

	long j;
	long j_max=(v->nh-v->nl+1)/2;
	double val=0.0;
	for (j=v->nl;j<=j_max;j++) {
		val=v->element[j];
		v->element[j]=v->element[v->nh-j+1];
		v->element[v->nh-j+1]=val;
	}

	return 0;
}

int invert_order_longvector(LONGVECTOR *v) {
/*!
 *
 * \author Emanuele Cordano
 * \date October 2009
 *
 *\param (LONGVECTOR *) the vector which is inverted
 */

	long j;
	long j_max=(v->nh-v->nl+1)/2;
	long val=0;
	for (j=v->nl;j<=j_max;j++) {
		val=v->element[j];
		v->element[j]=v->element[v->nh-j+1];
		v->element[v->nh-j+1]=val;
	}

	return 0;
}

#endif
