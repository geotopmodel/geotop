
/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file ncgt_utilities.c

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico 

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

#include "../../config.h"

#ifdef USE_NETCDF

#include "ncgt_utilities.h"

//char *copy_stringnames(const char *origin){
///*!
// * \author Emanuele Cordano
// * \date April 2009
// *
// *
// */
//char *dest;
//
////dest=origin;
//dest=(char *)malloc((size_t)(strlen(origin)+4)*sizeof(char));
//strcpy(dest,origin);
//return dest;
//}


//061009_s
//http://cf-pcmdi.llnl.gov/documents/cf-standard-names/standard-name-table/12/cf-standard-name-table.html
//void nc_add_cf_convention_header(const char *filename,const char *title, const char *institution,const char *source,const char *history,const char *references,const char *comment,const char *conventions){



int rotate180_y_doublematrix(GeoMatrix<double>& M) {
/*
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLEMATRIX *M) matrix to rotate
 */
	long r,c;
	double val=0.0;


	//long r_max=(M->nrh-M->nrl+1)/2;
	long r_max=(M.getRows()-1+1)/2;

	//for(r=M->nrl;r<=r_max;r++) {
	for(r=1;r<=r_max;r++) {
		//for (c=M->ncl;c<=M->nch;c++) {
		for (c=1;c<=M.getCols();c++) {
			//val=M->co[r][c];
			val=M(r,c);
			//M->co[r][c]=M->co[M->nrh-r+1][c];
			M(r,c)=M(M.getRows()-r+1,c);
			//M->co[M->nrh-r+1][c]=val;
			M(M.getRows()-r+1,c)=val;
		}
	}

	return 0;
}


int rotate180_y_doubletensor(GeoTensor<double>& M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLETENSOR *M) matrix to rotate
 */
	long r,c,l;
	double val=0.0;

	//long r_max=(M->nrh-M->nrl+1)/2;
	long r_max=(M.getRh()-1+1)/2;
//	printf("r_max=%ld %ld %ld \n",r_max,M->nrh,M->nrl);

	//for(l=M->ndl;l<=M->ndh;l++) {
	for(l=1;l<M.getDh();l++) {
		//for(r=M->nrl;r<=r_max;r++) {
		for(r=1;r<=r_max;r++) {
			//for (c=M->ncl;c<=M->nch;c++) {
			for (c=1;c<M.getCh();c++) {
				//val=M->co[l][r][c];
				val=M(l,r,c);
				//M->co[l][r][c]=M->co[l][M->nrh-r+1][c];
				M(l,r,c)=M(l,M.getRh()-r,c);
				//M->co[l][M->nrh-r+1][c]=val;
				M(l,M.getRh()-r,c)=val;
			}
		}
	}
	return 0;
}







#endif
