
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

#include "struct.geotop.h"
#include "../libraries/ascii/rw_maps.h"

extern T_INIT *UV;
extern long number_novalue;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map(long n, char *name, DOUBLEMATRIX *assign, PAR *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint){
	
	long r, c, i;
	DOUBLEMATRIX *M;
	char *temp;
	
	temp = namefile_i_we2(name, n);
	
	if(par->point_sim == 0){
		M = read_map(1, temp, Zdistr, UV, (double)number_novalue);
		for (r=1; r<=M->nrh; r++) {
			for (c=1; c<=M->nch; c++) {
				assign->co[r][c] = M->co[r][c];
			}
		}
		
	}else{
		M = read_map(1, temp, Zpoint, UV, (double)number_novalue);
		for(i=1; i<=par->r_points->nh; i++){
			r = par->r_points->co[i];
			c = par->c_points->co[i];
			assign->co[1][i] = M->co[r][c];
		}
	}
	
	free_doublematrix(M);
	free(temp);	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map_long(long n, char *name, LONGMATRIX *assign, PAR *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint){
	
	long r, c, i;
	DOUBLEMATRIX *M;
	char *temp;
	
	temp = namefile_i_we2(name, n);
	
	if(par->point_sim == 0){
		M = read_map(1, temp, Zdistr, UV, (double)number_novalue);
		for (r=1; r<=M->nrh; r++) {
			for (c=1; c<=M->nch; c++) {
				assign->co[r][c] = (long)M->co[r][c];
			}
		}
		
	}else{
		M = read_map(1, temp, Zpoint, UV, (double)number_novalue);
		for(i=1; i<=par->r_points->nh; i++){
			r = par->r_points->co[i];
			c = par->c_points->co[i];
			assign->co[1][i] = (long)M->co[r][c];
		}
	}
	
	free_doublematrix(M);
	free(temp);	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_tensor(long n, char *name, DOUBLETENSOR *assign, PAR *par, DOUBLEMATRIX *Zdistr, DOUBLEMATRIX *Zpoint){
	
	long r, c, i, l;
	DOUBLEMATRIX *M;
	char *temp1, *temp2;
	
	for (l=assign->ndl; l<=assign->ndh; l++) {
		
		temp1 = namefile_i_we2(name, n);
		temp2 = namefile_i_we(temp1, l);
		
		if(par->point_sim == 0){
			M = read_map(1, temp2, Zdistr, UV, (double)number_novalue);
			for (r=1; r<=M->nrh; r++) {
				for (c=1; c<=M->nch; c++) {
					assign->co[l][r][c] = M->co[r][c];
				}
			}
			
		}else{
			M = read_map(1, temp2, Zpoint, UV, (double)number_novalue);
			for(i=1; i<=par->r_points->nh; i++){
				r = par->r_points->co[i];
				c = par->c_points->co[i];
				assign->co[l][1][i] = M->co[r][c];
			}
		}
		
		free_doublematrix(M);
		free(temp2);	
		free(temp1);
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_tensor_channel(long n, char *name, DOUBLEMATRIX *assign, LONGVECTOR *r, LONGVECTOR *c, DOUBLEMATRIX *Zdistr){
	
	long ch, l;
	DOUBLEMATRIX *M;
	char *temp1, *temp2;
	
	for (l=assign->nrl; l<=assign->nrh; l++) {
		
		temp1 = namefile_i_we2(name, n);
		temp2 = namefile_i_we(temp1, l);
		
		M = read_map(1, temp2, Zdistr, UV, (double)number_novalue);
		
		for (ch=1; ch<=r->nh; ch++) {
			if(r->co[ch] > 0) assign->co[l][ch] = M->co[r->co[ch]][c->co[ch]];
		}
		
		free_doublematrix(M);
		free(temp2);	
		free(temp1);
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
