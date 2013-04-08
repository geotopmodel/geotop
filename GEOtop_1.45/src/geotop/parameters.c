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

#include "parameters.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

char **assign_string_parameter(FILE *f, long beg, long end, char **string_param, char **keyword){
	
	long i;
	char **a;
	
	a = (char**)malloc((end-beg)*sizeof(char*));
		
	for (i=0; i<end-beg; i++) {
		a[i] = assignation_string(f, i+beg, keyword, string_param);
	}
	
	return(a);
	
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

double assignation_number(FILE *flog, long i, long j, char **keyword, double **num_param, long *num_param_components, double default_value, short code_error){
	
	double a;
		
	if (j < num_param_components[i]) {
		a = num_param[i][j];
	}else {
		a = (double)number_novalue;
	}
	
	if((long)a == number_novalue){		
		if (code_error==0) {
			a = default_value;
			fprintf(flog,"%s[%ld] = %e (default) \n", keyword[i], j+1, a);
		}else{
			fprintf(flog, "%s[%ld] not assigned\n", keyword[i], j+1);
			fclose(flog);
			t_error("Fatal Error, See geotop.log!");	
		}
	}else{
		fprintf(flog,"%s[%ld] = %e \n", keyword[i], j+1, a);
	}
	
	return a;
}	

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/

char *assignation_string(FILE *f, long i, char **keyword, char **string_param){
	
	char *a;
	long j, dimstring = strlen(string_param[i]);
		
	a = (char*)malloc((dimstring+1)*sizeof(char));
	
	for (j=0; j<dimstring; j++) {
		a[j] = string_param[i][j];
	}
	a[dimstring] = 0;
	
	fprintf(f,"%s = %s\n", keyword[i], a);
	
	return(a);
}



/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		

short read_point_file(char *name, char **key_header, PAR *par, FILE *flog){
	
	DOUBLEMATRIX *chkpt2;
	double **points;
	long nlines, n, j;
	char *temp;
		
	if(existing_file_wext(name, textfile)==1){	
		temp = join_strings(name, textfile);
		printf("%s\n",temp);
		points = read_txt_matrix(temp, 33, 44, key_header, par->chkpt->nch, &nlines, flog);
		free(temp);
				
		chkpt2 = new_doublematrix(par->chkpt->nrh, par->chkpt->nch);
		copy_doublematrix(par->chkpt, chkpt2);
		free_doublematrix(par->chkpt);
		
		par->chkpt = new_doublematrix(nlines, chkpt2->nch);
		
		for (n=1; n<=nlines; n++) {
			for (j=1; j<=chkpt2->nch; j++) {
				par->chkpt->co[n][j] = points[n-1][j-1];
				if ( (long)par->chkpt->co[n][j] == number_novalue || (long)par->chkpt->co[n][j] == number_absent ) {
					if ( n <= chkpt2->nrh ) {
						par->chkpt->co[n][j] = chkpt2->co[n][j];
					}else {
						par->chkpt->co[n][j] = chkpt2->co[chkpt2->nrh][j];
					}
				}
			}
			
			if ( par->point_sim != 1 ){
				if( (long)par->chkpt->co[n][ptX] == number_novalue || (long)par->chkpt->co[n][ptY] == number_novalue ) {
					par->state_pixel = 0;
				}
			}
			
			free(points[n-1]);
		}
				
		free_doublematrix(chkpt2);
		free(points);
	}
	
	if (par->point_sim != 1){
		for (n=1; n<=par->chkpt->nrh; n++) {
			if ( (long)par->chkpt->co[n][ptX] == number_novalue || (long)par->chkpt->co[n][ptY] == number_novalue) {
				fprintf(flog,"Warning: The points to plot specific results are not completely specified\n");
				fprintf(flog,"Output for single point output is deactivated.\n");
				printf("Warning: The points to plot specific results are not completely specified\n");
				printf("Output for single point output is deactivated.\n");
				par->state_pixel = 0;
			}
		}
	}
	
	
	return 1;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		
	
short read_meteostations_file(LONGVECTOR *i, METEO_STATIONS *S, char *name, char **key_header, FILE *flog){
	
	double **M;
	long nlines, n, j, k;
	char *temp;
		
	if(existing_file_wext(name, textfile)==1){	
		temp = join_strings(name, textfile);
		M = read_txt_matrix(temp, 33, 44, key_header, 8, &nlines, flog);
		free(temp);
				
		for (j=1; j<=i->nh; j++) {
			for (n=1; n<=nlines; n++) {
				if ((long)M[n-1][0] == i->co[j]) {
					for (k=1; k<8; k++) {
						if ((long)M[n-1][k] != number_novalue && (long)M[n-1][k] != number_absent) {
							if (k==1) {
								S->E->co[j] = M[n-1][k];
							}else if (k==2) {
								S->N->co[j] = M[n-1][k];
							}else if (k==3) {
								S->lat->co[j] = M[n-1][k];
							}else if (k==4) {
								S->lon->co[j] = M[n-1][k];
							}else if (k==5) {
								S->Z->co[j] = M[n-1][k];
							}else if (k==6) {
								S->sky->co[j] = M[n-1][k];
							}else if (k==7) {
								S->ST->co[j] = M[n-1][k];
							}
						}
					}
				}					
			}
		}
		
		for (n=1; n<=nlines; n++) {
			free(M[n-1]);
		}
		
		free(M);
		
		return 1;
	}else {
		return 0;
	}
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/***********************************************************/		
