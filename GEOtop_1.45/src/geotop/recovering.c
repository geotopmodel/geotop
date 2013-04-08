
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225-9 'Moab' - 24 Aug 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225-9 'Moab'
 
 GEOtop 1.225-9 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225-9 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "recovering.h"
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map(short old, long n, char *name, DOUBLEMATRIX *assign, PAR *par, DOUBLEMATRIX *Zdistr){
	
	long r, c;
	//long i;
	DOUBLEMATRIX *M;
	char *temp, *temp2;
	
	temp = namefile_i_we2(name, n);
	
	if (old == 1) {
		temp2 = join_strings(temp, ".old");
		free(temp);
		temp = assign_string(temp2);
		free(temp2);
	}
		
	M = read_map(1, temp, Zdistr, UV, (double)number_novalue);
	for (r=1; r<=M->nrh; r++) {
		for (c=1; c<=M->nch; c++) {
			assign->co[r][c] = M->co[r][c];
		}
	}
			
	free_doublematrix(M);
	free(temp);	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map_vector(short old, long n, char *name, DOUBLEVECTOR *assign, LONGMATRIX *rc, PAR *par, DOUBLEMATRIX *Zdistr){
	
	long i, r, c;
	DOUBLEMATRIX *M;
	char *temp, *temp2;
	
	temp = namefile_i_we2(name, n);
	
	if (old == 1) {
		temp2 = join_strings(temp, ".old");
		free(temp);
		temp = assign_string(temp2);
		free(temp2);
	}
	
	M = read_map(1, temp, Zdistr, UV, (double)number_novalue);
	for (i=1; i<=rc->nrh; i++) {
		r = rc->co[i][1];
		c = rc->co[i][2];
		assign->co[i] = M->co[r][c];
	}
		
	free_doublematrix(M);
	free(temp);	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map_long(short old, long n, char *name, LONGMATRIX *assign, PAR *par, DOUBLEMATRIX *Zdistr){
	
	long r, c;
	DOUBLEMATRIX *M;
	char *temp, *temp2;
	
	temp = namefile_i_we2(name, n);
	
	if (old == 1) {
		temp2 = join_strings(temp, ".old");
		free(temp);
		temp = assign_string(temp2);
		free(temp2);
	}
	
	M = read_map(1, temp, Zdistr, UV, (double)number_novalue);
	for (r=1; r<=M->nrh; r++) {
		for (c=1; c<=M->nch; c++) {
			assign->co[r][c] = (long)M->co[r][c];
		}
	}
	
	free_doublematrix(M);
	free(temp);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_tensor(short old, long n, char *name, DOUBLETENSOR *assign, PAR *par, DOUBLEMATRIX *Zdistr){
	
	long r, c, l;
	//long i;
	DOUBLEMATRIX *M;
	char *temp1, *temp2, *temp3;
	
	for (l=assign->ndl; l<=assign->ndh; l++) {
		
		temp1 = namefile_i_we2(name, n);
		temp2 = namefile_i_we(temp1, l);
		
		if (old == 1) {
			temp3 = join_strings(temp2, ".old");
			free(temp2);
			temp2 = assign_string(temp3);
			free(temp3);
		}
		
		M = read_map(1, temp2, Zdistr, UV, (double)number_novalue);
		
		for (r=1; r<=M->nrh; r++) {
			for (c=1; c<=M->nch; c++) {
				assign->co[l][r][c] = M->co[r][c];
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

void assign_recovered_tensor_vector(short old, long n, char *name, DOUBLEMATRIX *assign, LONGMATRIX *rc, PAR *par, DOUBLEMATRIX *Zdistr){
	
	long r, c, i, l;
	DOUBLEMATRIX *M;
	char *temp1, *temp2, *temp3;
	
	for (l=assign->nrl; l<=assign->nrh; l++) {
				
		temp1 = namefile_i_we2(name, n);
		temp2 = namefile_i_we(temp1, l);
		
		if (old == 1) {
			temp3 = join_strings(temp2, ".old");
			free(temp2);
			temp2 = assign_string(temp3);
			free(temp3);
		}
		
		M = read_map(1, temp2, Zdistr, UV, (double)number_novalue);
		for (i=1; i<=rc->nrh; i++) {
			r = rc->co[i][1];
			c = rc->co[i][2];
			assign->co[l][i] = M->co[r][c];
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

void assign_recovered_tensor_channel(short old, long n, char *name, DOUBLEMATRIX *assign, LONGVECTOR *r, LONGVECTOR *c, DOUBLEMATRIX *Zdistr){
	
	long ch, l;
	DOUBLEMATRIX *M;
	char *temp1, *temp2, *temp3;
	
	for (l=assign->nrl; l<=assign->nrh; l++) {
				
		temp1 = namefile_i_we2(name, n);
		temp2 = namefile_i_we(temp1, l);
		
		if (old == 1) {
			temp3 = join_strings(temp2, ".old");
			free(temp2);
			temp2 = assign_string(temp3);
			free(temp3);
		}
		
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

void recover_run_averages(short old, DOUBLEMATRIX *A, char *name, DOUBLEMATRIX *LC, LONGMATRIX *rc, PAR *par, long n){

	DOUBLEMATRIX *M;
	long j, l;
	
	M = new_doublematrix(n, par->total_pixel);
	assign_recovered_tensor_vector(old, par->recover, name, M, rc, par, LC);
	for (j=1; j<=par->total_pixel; j++) {
		if (par->jplot->co[j] > 0) {
			for (l=1; l<=n; l++) {
				A->co[par->jplot->co[j]][l] = M->co[l][j];
			}
		}
	}
	free_doublematrix(M);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void print_run_averages_for_recover(DOUBLEMATRIX *A, char *name, long **j_cont, PAR *par, long n, long nr, long nc){

	DOUBLEMATRIX *M;
	long j, l;

	M = new_doublematrix(n, par->total_pixel);
	initialize_doublematrix(M, (double)number_novalue);
	for (j=1; j<=par->total_pixel; j++) {
		if (par->jplot->co[j] > 0) {
			for (l=1; l<=n; l++) {
				M->co[l][j] = A->co[par->jplot->co[j]][l];
			}
		}
	}
	for (l=1; l<=n; l++) {
		rename_tensorseries(1, l, 0, name);
		write_tensorseries_vector(1, l, 0, name, 0, par->format_out, M, UV, number_novalue, j_cont, nr, nc);
	}
	free_doublematrix(M);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
