/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
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

double assignation_number(FILE *f, long i, long j, char **keyword, double **num_param, long *num_param_components, double default_value, short code_error){
	
	double a;
		
	if (j < num_param_components[i]) {
		a = num_param[i][j];
	}else {
		a = (double)number_novalue;
	}
	
	if((long)a == number_novalue){		
		if (code_error==0) {
			a = default_value;
			fprintf(f,"%s[%ld] = %e (default) \n", keyword[i], j+1, a);
		}else{
			f = fopen(FailedRunFile, "w");
			fprintf(f, "%s[%ld] not assigned\n", keyword[i], j+1);
			fclose(f);
			t_error("Fatal Error!");	
		}
	}else{
		fprintf(f,"%s[%ld] = %e \n", keyword[i], j+1, a);
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

short read_soil_parameters(char *name, char **key_header, SOIL *sl, FILE *flog){
	
	short ok;
	long i, j, k, n, nlines, nlinesprev;
	char *temp;
	double **soildata;
	DOUBLETENSOR *old_sl_par;
	FILE *f;
	
	//look if there is at least 1 soil file
	i = 0;
	ok = 0;
	nlinesprev = -1;
	
	do{
		
		temp = namefile_i_we2(name, i+1);
		
		if(existing_file_wext(temp, textfile)==1){
			free(temp);
			ok = 1;
			temp = namefile_i(name, i+1);
			nlines = count_lines(temp, 33, 44);
			if (nlinesprev >= 0 && nlines != nlinesprev){
				f = fopen(FailedRunFile, "w");
				fprintf(f,"Error:: The file %s with soil paramaters has a number of layers %ld, which different from the numbers %ld of the other soil parameter files\n",temp,nlines,nlinesprev);
				fprintf(f,"In GEOtop it is only possible to have the same number of layers in any soil parameter files\n");
				fclose(f);
				t_error("Fatal Error! GEOtop is closed. See failing report.");
			}
			nlinesprev = nlines;
		}else {
			if (i==0 && strcmp(name, string_novalue) != 0) {
				f = fopen(FailedRunFile, "w");
				fprintf(f,"Error:: Soil file %s not existing.\n",name);
				fclose(f);
				t_error("Fatal Error! GEOtop is closed. See failing report.");	
			}
		}

			
		free(temp);
		i++;
		
	}while (ok == 0 && i < sl->pa->ndh);
	
	if (ok == 1){
		
		//save sl->pa in a new doubletensor and deallocate
		old_sl_par = new_doubletensor(sl->pa->ndh, sl->pa->nrh, sl->pa->nch);
		for (i=1; i<=sl->pa->ndh; i++) {
			for (n=1; n<=sl->pa->nrh; n++) {
				for (j=1; j<=sl->pa->nch; j++) {
					old_sl_par->co[i][n][j] = sl->pa->co[i][n][j];
				}
			}
		}
		free_doubletensor(sl->pa);
		
		//reallocate
		sl->pa = new_doubletensor(old_sl_par->ndh, old_sl_par->nrh, nlines);
		
		for (i=1; i<=sl->pa->ndh; i++) {
			
			//read files
			temp = namefile_i_we2(name, i);
			
			if(existing_file_wext(temp, textfile)==1){
				free(temp);
				temp = namefile_i(name, i);
				soildata = read_txt_matrix(temp, 33, 44, key_header, nsoilprop, &nlines, flog);
			}else {
				soildata = (double**)malloc(nlines*sizeof(double*));
				for (j=0; j<nlines; j++) {
					k = (long)nsoilprop;
					soildata[j] = (double*)malloc(k*sizeof(double));
					for (n=0; n<k; n++) {
						soildata[j][n] = (double)number_absent;
					}
				}
			}

			free(temp);
			
			//assign soildata to soil->pa
			for (n=1; n<=nsoilprop; n++) {
				for (j=1; j<=sl->pa->nch; j++) { //j is the layer index
					sl->pa->co[i][n][j] = soildata[j-1][n-1];
				}
			}
			 
			//deallocate soildata
			for (j=0; j<nlines; j++) {
				free(soildata[j]);
			}
			free(soildata);
			
			//fix layer thickness
			n = jdz;
			for (j=1; j<=sl->pa->nch; j++) { //j is the layer index
				if ((long)sl->pa->co[i][n][j] != number_novalue && (long)sl->pa->co[i][n][j] != number_absent) {
					if ( i > 1 && fabs( sl->pa->co[i][n][j] - sl->pa->co[i-1][n][j] ) > 1.E-5 )  {
						f = fopen(FailedRunFile, "w");
						fprintf(f,"Error:: For soil type %ld it has been given a set of soil layer thicknesses different from the other ones.\n",i);
						fprintf(f,"In GEOtop it is only possible to have the soil layer discretization in any soil parameter files.\n");
						fclose(f);
						t_error("Fatal Error! GEOtop is closed. See failing report.");	
					}
				}else if (i == 1) {
					if (j <= old_sl_par->nch) {
						sl->pa->co[i][n][j] = old_sl_par->co[i][n][j];
					}else {
						sl->pa->co[i][n][j] = sl->pa->co[i][n][j-1];
					}
				}else {
					sl->pa->co[i][n][j] = sl->pa->co[i-1][n][j];
				}
			}
			
			//all other variables
			for (n=1; n<=nsoilprop; n++) {
				if (n != jdz) {
					for (j=1; j<=sl->pa->nch; j++) { //j is the layer index
						if ((long)sl->pa->co[i][n][j] == number_novalue || (long)sl->pa->co[i][n][j] == number_absent) {
							if (j <= old_sl_par->nch) {
								sl->pa->co[i][n][j] = old_sl_par->co[i][n][j];
							}else {
								sl->pa->co[i][n][j] = sl->pa->co[i][n][j-1];
							}
						}
					}
				}
			}
			
			//field capacity and wilting point
			for (j=1; j<=sl->pa->nch; j++){
				if( (long)sl->pa->co[i][jfc][j] == number_novalue){
					sl->pa->co[i][jfc][j] = teta_psi( (-1./3.)*1.E5/g, 0., sl->pa->co[i][jsat][j], sl->pa->co[i][jres][j], sl->pa->co[i][ja][j], 
													 sl->pa->co[i][jns][j], 1.-1./sl->pa->co[i][jns][j], PsiMin, sl->pa->co[i][jss][j]);
				}
				
				if( (long)sl->pa->co[i][jwp][j] == number_novalue){
					sl->pa->co[i][jwp][j] = teta_psi( -15.*1.E5/g, 0., sl->pa->co[i][jsat][j], sl->pa->co[i][jres][j], sl->pa->co[i][ja][j], 
													 sl->pa->co[i][jns][j], 1.-1./sl->pa->co[i][jns][j], PsiMin, sl->pa->co[i][jss][j]);		
				}
			}
			
			//pressure
			ok = 1;
			for (j=1; j<=sl->pa->nch; j++){
				if( (long)sl->pa->co[i][jpsi][j] == number_novalue) ok = 0;
			}
			if (ok == 1) sl->init_water_table_depth->co[i] = (double)number_novalue;
		}
		
		free_doubletensor(old_sl_par);
			
	}
	
	//write on the screen the soil paramater
	fprintf(flog,"\n");
	k = (long)nmet;
	fprintf(flog,"Soil Layers: %ld\n",sl->pa->nch);
	for (i=1; i<=sl->pa->ndh; i++) {
		fprintf(flog,"-> Soil Type: %ld\n",i);
		for (n=1; n<=nsoilprop; n++) {
			fprintf(flog,"%s: ",keywords_char[k+n-1]);
			for (j=1; j<=sl->pa->nch; j++) {
				fprintf(flog,"%f(%.2e)",sl->pa->co[i][n][j],sl->pa->co[i][n][j]);
				if(j<sl->pa->nch)fprintf(flog,", ");
			}
			fprintf(flog,"\n");
		}
	}
	
	//bedrock
	old_sl_par = new_doubletensor(1, sl->pa_bed->nrh, sl->pa_bed->nch);
	for (n=1; n<=sl->pa_bed->nrh; n++) {
		for (j=1; j<=sl->pa_bed->nch; j++) {
			old_sl_par->co[1][n][j] = sl->pa_bed->co[1][n][j];
		}
	}
	free_doubletensor(sl->pa_bed);
	
	sl->pa_bed = new_doubletensor(sl->pa->ndh, sl->pa->nrh, sl->pa->nch);
	for (i=1; i<=sl->pa_bed->ndh; i++) {
		for (n=1; n<=sl->pa_bed->nrh; n++) {
			if (n == jdz) {
				for (j=1; j<=sl->pa_bed->nch; j++) {
					sl->pa_bed->co[i][n][j] = sl->pa->co[1][n][j];
				}
			}else {
				for (j=1; j<=sl->pa_bed->nch; j++) {
					if (j <= old_sl_par->nch) {
						sl->pa_bed->co[i][n][j] = old_sl_par->co[1][n][j];
					}else {
						sl->pa_bed->co[i][n][j] = sl->pa_bed->co[i][n][j-1];
					}
				}
				for (j=1; j<=sl->pa_bed->nch; j++) {
					if ( (long)sl->pa_bed->co[i][n][j] == number_novalue ) sl->pa_bed->co[i][n][j] = sl->pa->co[i][n][j];
				}
			}
		}		
	}
	free_doubletensor(old_sl_par);
	
	fprintf(flog,"\n");
	k = (long)nmet;
	fprintf(flog,"Soil Bedrock Layers: %ld\n",sl->pa->nch);
	for (i=1; i<=sl->pa_bed->ndh; i++) {
		fprintf(flog,"-> Soil Type: %ld\n",i);
		for (n=1; n<=nsoilprop; n++) {
			fprintf(flog,"%s: ",keywords_char[k+n-1]);
			for (j=1; j<=sl->pa->nch; j++) {
				fprintf(flog,"%f(%.2e)",sl->pa_bed->co[i][n][j],sl->pa_bed->co[i][n][j]);
				if(j<sl->pa->nch)fprintf(flog,", ");
			}
			fprintf(flog,"\n");
		}
	}
	fprintf(flog,"\n");
	
	return 1;
	
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
		points = read_txt_matrix(temp, 33, 44, key_header, par->chkpt->nch, &nlines, flog);
		free(temp);
				
		chkpt2 = new_doublematrix(par->chkpt->nrh, par->chkpt->nch);
		copy_doublematrix(par->chkpt, chkpt2);
		free_doublematrix(par->chkpt);
		
		par->chkpt = new_doublematrix(nlines, chkpt2->nch);
		for (n=1; n<=nlines; n++) {
			for (j=1; j<=chkpt2->nch; j++) {
				par->chkpt->co[n][j] = points[n-1][j-1];
				printf("n=%d j=%d points[n-1][j-1] %f \n",n,j, points[n-1][j-1]);
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
