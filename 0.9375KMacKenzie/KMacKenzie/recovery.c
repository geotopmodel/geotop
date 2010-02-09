/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "recovery.h"

short file_exists(char *filename)
{
	FILE *file;
	if ((file = fopen(filename, "r"))){
		fclose(file);
		return 1;
	}
	return 0;
}

short check_file_exists(char* filename)
{
	short retval = file_exists(filename);
	if (retval==1) filecounter++;
	
	return retval;
}

void compare_matrix(DOUBLEMATRIX *one, DOUBLEMATRIX *two)
{
	long ii, jj;

	printf("Comparing: %ld x %ld, %ld x %ld\n", one->nrh, one->nch, two->nrh, two->nch);

	for (ii=1; ii <= one->nrh; ii++){
		for (jj=1; jj <= one->nch; jj++){
			if (one->co[ii][jj] != two->co[ii][jj]){
				printf("UNGLEICH\n");
			}
		}
	}
}

void compare_lmatrix(LONGMATRIX *one, LONGMATRIX *two)
{
	long ii, jj;

	printf("Comparing: %ld x %ld, %ld x %ld\n", one->nrh, one->nch, two->nrh, two->nch);

	for (ii=1; ii <= one->nrh; ii++){
		for (jj=1; jj <= one->nch; jj++){
			if (one->co[ii][jj] != two->co[ii][jj]){
				printf("UNGLEICH\n");
			} 
		}
	}
}

void compare_tensor(DOUBLETENSOR *one, DOUBLETENSOR *two)
{
	long ii, jj, kk;

	//printf("Dimensions: %ld x %ld x %ld,   %ld x %ld x %ld\n", one->ndh, one->nrh, 
		  //	  one->nch, two->ndh, two->nrh, two->nch);

	for (kk=1; kk <= one->ndh; kk++){
		for (ii=1; ii <= one->nrh; ii++){
			for (jj=1; jj <= one->nch; jj++){
				if (one->co[kk][ii][jj] != two->co[kk][ii][jj]){
					//printf("UNGLEICH\n");
					printf("%.*f != %.*f\n", (int) LDBL_DIG, one->co[kk][ii][jj], 
						  (int) LDBL_DIG, two->co[kk][ii][jj]);
				} else {
					//printf("%f == %f\n", one->co[kk][ii][jj], two->co[kk][ii][jj]);
				}
			}
		}
	}
}

void read_tensor3D(char *_filename, char *timestamp, DOUBLETENSOR **tensor_out)
{
	//basic idea: read doublematrices and copy them into the tensor
	long layers=0, ii=0, jj=0, kk=0;
	char filename[256];
	DOUBLEMATRIX **mymatrices = NULL;


	if ((*tensor_out) != NULL)
		(*tensor_out) = free_doubletensor(*tensor_out);

	do {
		char LLLLL[ ]={"LLLLL"};
		DOUBLEMATRIX *mymatrix = NULL;
		write_suffix(LLLLL, (layers + 1), 1);
		sprintf(filename, "%s%s_%s", _filename, LLLLL, timestamp);

		//printf("Trying to read %s\n", filename);

		read_grass_map(filename, 0, (void*)&mymatrix);
		if (mymatrix == NULL){ //file does not exist
			break;
		} else {
			//realloc memory and add another matrix
			layers++;
			mymatrices = (DOUBLEMATRIX **) realloc(mymatrices, layers*sizeof(DOUBLEMATRIX*));
			mymatrices[layers-1] = mymatrix;
		}
		
	} while(layers < 999); // this is a hypothetical condition only

	//Now allocate tensor and copy data of matrices into tensor
	if (layers > 0)
		(*tensor_out) = new_doubletensor(layers, mymatrices[0]->nrh, mymatrices[0]->nch);
	else
		return;

	for (ii=1; ii<=(*tensor_out)->ndh; ii++){
		//printf("processing layer %ld\n", ii);
		for (jj=1; jj<=(*tensor_out)->nrh; jj++){
			for (kk=1; kk<=(*tensor_out)->nch; kk++){
				(*tensor_out)->co[ii][jj][kk] = mymatrices[ii-1]->co[jj][kk];
			}
		}
	}

	//Deallocate mymatrices
	for (ii=0; ii<layers; ii++) mymatrices[ii] = free_doublematrix(mymatrices[ii]);
	free(mymatrices);
}

void read_grass_map(char* _filename, short type, void **matrix_out)
{
	if (type > 2){
		fprintf(stderr, "Unrecognized type when invoking read_grass_map: type number %hd", type);
		exit(1);		
	}

	char filename[256];
	char line[2000];
	long nch, nrh, ii;
	char *token;

	FILE *file;
	DOUBLEMATRIX *doublematrix = NULL;
	LONGMATRIX   *longmatrix   = NULL;
	SHORTMATRIX  *shortmatrix  = NULL;

	sprintf(filename, "%s%s", _filename, ".grass"); //Construct correct filename

	if ((*matrix_out) != NULL){
		if (type == 0){
			doublematrix = (DOUBLEMATRIX*) (*matrix_out);
			doublematrix = free_doublematrix(doublematrix);
		} else if (type == 1){
			longmatrix = (LONGMATRIX*) (*matrix_out);
			longmatrix = free_longmatrix(longmatrix);
		} else if (type == 2){
			shortmatrix = (SHORTMATRIX*) (*matrix_out);
			shortmatrix = free_shortmatrix(shortmatrix);
		}		
	}

	if (check_file_exists(filename)){
		//open it, skip first 4 lines
		file = fopen(filename, "r");
		if ((fgets(line, 2000, file)) == NULL) t_error("Error reading grass file");
		if ((fgets(line, 2000, file)) == NULL) t_error("Error reading grass file");
		if ((fgets(line, 2000, file)) == NULL) t_error("Error reading grass file");
		if ((fgets(line, 2000, file)) == NULL) t_error("Error reading grass file");
		fscanf(file, "rows:%ld\n", &nrh);
		fscanf(file, "cols:%ld\n", &nch);
		//printf("Dimensions: %ld x %ld\n", nrh, nch);

		if (type == 0)	doublematrix = new_doublematrix(nrh, nch);
		else if (type == 1)	longmatrix = new_longmatrix(nrh, nch);
		else if (type == 2) shortmatrix = new_shortmatrix(nrh, nch);

		for (ii=1; ii<=nrh; ii++){
			fgets(line, 2000, file);
			size_t last = strlen(line) - 1;
			if (last >=0){
				if (line[last] == '\n') {
					line[last] = '\0';
				}
			}

			token = strtok(line, " ");
			long ccount=0;
			while (token != NULL){
				ccount++;
				if (strcmp(token, "*") == 0){
					if (type == 0) doublematrix->co[ii][ccount] = -9999.0;
					else if (type == 1) longmatrix->co[ii][ccount] = -9999;
					else if (type == 2) shortmatrix->co[ii][ccount] = -9999;
				} else {
					if (type == 0) doublematrix->co[ii][ccount] = atof(token);
					else if (type == 1) longmatrix->co[ii][ccount] = atol(token);
					else if (type == 2) shortmatrix->co[ii][ccount] = atoi(token);
				}

				//if (type == 0)
				//printf("Line %ld, Col %ld (token '%s'): %f\n", ii, ccount, token, doublematrix->co[ii][ccount]);

				token = strtok (NULL, " ");
			}
			if (ccount != nch){
				fprintf(stderr, "Error parsing GRASS file %s, line number: %ld\n", filename, ii);
				exit(1);
			}
		}

		if (type == 0) (*matrix_out) = doublematrix;
		else if (type == 1) (*matrix_out) = longmatrix;
		else if (type == 2) (*matrix_out) = shortmatrix;
		fclose(file);
	} else {
		return; //file does not exist
	}
}

void write_double_to_file(char* filename, double d)
{
	FILE *f;
	f=t_fopen(filename, "w");
	fprintf(f, "%.*f\n", (int)LDBL_DIG, d);
	t_fclose(f);
}

void write_long_to_file(char* filename, long l)
{
	FILE *f;
	f=t_fopen(filename, "w");
	fprintf(f, "%ld\n", l);
	t_fclose(f);
}

void write_short_to_file(char* filename, short s)
{
	FILE *f;
	f=t_fopen(filename, "w");
	fprintf(f, "%hd\n", s);
	t_fclose(f);
}

void read_double_from_file(char* filename, double* d)
{
	*d = 0.0;

	if (check_file_exists(filename)){
		FILE *file;
		file = fopen(filename, "r");

		//Now read one double value and close the file
		if ((fscanf(file, "%lf", d)) != 1){
			fclose(file);
			fprintf(stderr, "Error when trying to read file '%s', was expecting one double value", filename);
			exit(1);			
		}
		fclose(file);
	}
}

void read_long_from_file(char* filename, long* l)
{
	*l = 0;

	if (check_file_exists(filename)){
		FILE *file;
		file = fopen(filename, "r");

		//Now read one double value and close the file
		if ((fscanf(file, "%ld", l)) != 1){
			fclose(file);
			fprintf(stderr, "Error when trying to read file '%s', was expecting one long value", filename);
			exit(1);			
		}
		fclose(file);
	}
}

void read_short_from_file(char* filename, short* s)
{
	*s = 0;

	if (check_file_exists(filename)){
		FILE *file;
		file = fopen(filename, "r");

		//Now read one double value and close the file
		if ((fscanf(file, "%hd", s)) != 1){
			fclose(file);
			fprintf(stderr, "Error when trying to read file '%s', was expecting one short value", filename);
			exit(1);			
		}
		fclose(file);
	}
}

void read_vector_from_file(char* filename, short type, void **vec_out)
{
	//	type=0  floating point
	//	type=1  long
	//	type=2  short

	if (type > 2){
		fprintf(stderr, "Unrecognized type when invoking write_vector_to_file: type number %d", type);
		exit(1);		
	}

	DOUBLEVECTOR *doublevec = NULL;
	LONGVECTOR *longvec     = NULL;
	SHORTVECTOR *shortvec   = NULL;

	//printf("Reading vector: %s\n", filename);
	if ((*vec_out) != NULL){
		if (type == 0){
			doublevec = (DOUBLEVECTOR*) (*vec_out);
			doublevec = free_doublevector(doublevec);
		} else if (type == 1){
			longvec = (LONGVECTOR*) (*vec_out);
			longvec = free_longvector(longvec);
		} else if (type == 2){
			shortvec = (SHORTVECTOR*) (*vec_out);
			shortvec = free_shortvector(shortvec);
		}		
	}

	FILE *file;
	long lines = 0;
	char c;

	if (check_file_exists(filename)){
		file = fopen(filename, "r");
		
		//Count number of lines == number of elements
		do {
			c = fgetc(file);
			if (c == '\n') lines++;
		} while (c != EOF);
		//printf(" --> Number of lines: %ld\n", lines);
		fclose (file);
	} else {
		//printf(" --> File does not exist!\n");
		return; //file does not exist
	}

	//Now allocate and fill vector
	file = fopen(filename, "r");
	long counter = 1;

	if (type == 0){
		doublevec = new_doublevector(lines);
		while(feof(file)== 0){
			fscanf(file,"%lf\n", &doublevec->co[counter]);
			//printf("Read vector value[%ld] = %f\n", counter, doublevec->co[counter]);
			counter++;
		}
	} else if (type == 1){	
		longvec = new_longvector(lines);
		while(feof(file)== 0){
			fscanf(file,"%ld\n", &longvec->co[counter]);
			counter++;
		}
	} else if (type == 2){	
		shortvec = new_shortvector(lines);
		while(feof(file)== 0){
			fscanf(file,"%hd\n", &shortvec->co[counter]);
			counter++;
		}
	}
	fclose(file);

	if (type == 0) (*vec_out) = doublevec;
	else if (type == 1) (*vec_out) = longvec;
	else if (type == 2) (*vec_out) = shortvec;
}

void write_vector_to_file(char* filename, short type, void *vec_in)
{
	//	type=0  floating point
	//	type=1  long
	//	type=2  short

	if (vec_in == NULL)
		return;

	DOUBLEVECTOR *doublevec = NULL;
	LONGVECTOR *longvec     = NULL;
	SHORTVECTOR *shortvec     = NULL;
	long nh = 0;

	if (type == 0){
		doublevec = (DOUBLEVECTOR *)vec_in;
		nh = doublevec->nh;
	} else if (type == 1){
		longvec = (LONGVECTOR *)vec_in;
		nh = longvec->nh;
	} else if (type == 2){
		shortvec = (SHORTVECTOR *)vec_in;
		nh = shortvec->nh;
	} else {
		fprintf(stderr, "Unrecognized type when invoking write_vector_to_file: type number %d", type);
		exit(1);
	}

	FILE *f;
	long ii=0;

	f=t_fopen(filename, "w");	
	for (ii=1; ii<=nh; ii++){
		if (type == 0){
			fprintf(f, "%.*f\n", (int)LDBL_DIG, doublevec->co[ii]);
		} else if (type == 1){
			fprintf(f, "%ld\n", longvec->co[ii]);
		} else if (type == 2){
			fprintf(f, "%hd\n", shortvec->co[ii]);
		}
	}
	t_fclose(f);
}

void write_map2D(char *filename, char *timestamp, short type, void *M){
//	type=0  floating point
//	type=1  integer
//	type=2  short

	if (M == NULL)
		return;

	char tmp_fname[128];
	sprintf(tmp_fname, "%s_%s", filename, timestamp); 
	write_grassascii(tmp_fname, type, M, UV);
}

void read_map2D(char *filename, char *timestamp, short type, void **matrix_out)
{
	char tmp_fname[128];
	sprintf(tmp_fname, "%s_%s", filename, timestamp);
	read_grass_map(tmp_fname, type, matrix_out);
}

void read_vector1D(char *filename, char *timestamp, short type, void **vector_out)
{
	char tmp_fname[128];
	sprintf(tmp_fname, "%s_%s", filename, timestamp); 
	read_vector_from_file(tmp_fname, type, vector_out);
}

void read_double(char *filename, char *timestamp, double *d)
{
	char tmp_fname[128];
	sprintf(tmp_fname, "%s_%s", filename, timestamp); 
	read_double_from_file(tmp_fname, d);
}

void read_long(char *filename, char *timestamp, long *d)
{
	char tmp_fname[128];
	sprintf(tmp_fname, "%s_%s", filename, timestamp); 
	read_long_from_file(tmp_fname, d);
}

void read_short(char *filename, char *timestamp, short *d)
{
	char tmp_fname[128];
	sprintf(tmp_fname, "%s_%s", filename, timestamp); 
	read_short_from_file(tmp_fname, d);
}

void write_output_variables(char* timestamp)
{
	FILE *f;
	char tmp_fname[128];

	sprintf(tmp_fname, "%s_%s", "output/rec/output_variables", timestamp); 

	f=t_fopen(tmp_fname, "w");
	fprintf(f, "%.*f ", (int)LDBL_DIG, wt0_basin);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Ssup);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Ssub);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Rout);
	fprintf(f, "%.*f ", (int)LDBL_DIG, R_G);
	fprintf(f, "%.*f ", (int)LDBL_DIG, S_ch0);
	fprintf(f, "%.*f ", (int)LDBL_DIG, S_ch1);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Qsub_ch);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Qsup_ch);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Q_G);
	fprintf(f, "%.*f ", (int)LDBL_DIG, SWE_previous);
	fprintf(f, "%.*f ", (int)LDBL_DIG, GWE_previous);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Smelt);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Ssubl);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Sevap);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Gmelt);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Gsubl);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Gevap);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Smelt_previous);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Ssubl_previous);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Sevap_previous);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Gmelt_previous);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Gsubl_previous);
	fprintf(f, "%.*f ", (int)LDBL_DIG, Gevap_previous);
	fprintf(f, "%ld\n",  isavings);

	t_fclose(f);
}

void read_output_variables(char* timestamp)
{
	FILE *f;
	char tmp_fname[128];

	sprintf(tmp_fname, "%s_%s", "output/rec/output_variables", timestamp); 

	if (check_file_exists(tmp_fname)){
		f = t_fopen(tmp_fname, "r");
		fscanf(f, "%lf", &wt0_basin);
		fscanf(f, "%lf", &Ssup);
		fscanf(f, "%lf", &Ssub);
		fscanf(f, "%lf", &Rout);
		fscanf(f, "%lf", &R_G);
		fscanf(f, "%lf", &S_ch0);
		fscanf(f, "%lf", &S_ch1);
		fscanf(f, "%lf", &Qsub_ch);
		fscanf(f, "%lf", &Qsup_ch);
		fscanf(f, "%lf", &Q_G);
		fscanf(f, "%lf", &SWE_previous);
		fscanf(f, "%lf", &GWE_previous);
		fscanf(f, "%lf", &Smelt);
		fscanf(f, "%lf", &Ssubl);
		fscanf(f, "%lf", &Sevap);
		fscanf(f, "%lf", &Gmelt);
		fscanf(f, "%lf", &Gsubl);
		fscanf(f, "%lf", &Gevap);
		fscanf(f, "%lf", &Smelt_previous);
		fscanf(f, "%lf", &Ssubl_previous);
		fscanf(f, "%lf", &Sevap_previous);
		fscanf(f, "%lf", &Gmelt_previous);
		fscanf(f, "%lf", &Gsubl_previous);
		fscanf(f, "%lf", &Gevap_previous);
		fscanf(f, "%ld", &isavings);
		t_fclose(f);
	}
}

void recover_simulation(char *timestamp, TIMES *times, WATER *wat, CHANNEL *cnet, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac)
{
	printf("Attempting to recover simulation from data at timestamp %s\n", timestamp);
	filecounter=0;

	/*****************************************************************************/
	/* Actual reading of recovery files follows                                  */
	/*****************************************************************************/
	/* 1. READING SOIL DATA */
	read_map2D("output/rec/soil_type", timestamp, 2, (void**)&sl->type);

	// The following tensors have a depth of Nl (except sl->pa which has a depth of 1)
	read_tensor3D("output/rec/soil_pa", timestamp, &sl->pa);
	read_tensor3D("output/rec/soil_P", timestamp, &sl->P);
	read_tensor3D("output/rec/soil_T", timestamp, &sl->T);
	read_tensor3D("output/rec/soil_thice", timestamp, &sl->thice);
	read_tensor3D("output/rec/soil_J", timestamp, &sl->J);
	read_tensor3D("output/rec/soil_Tav", timestamp, &sl->Tav);
	read_tensor3D("output/rec/soil_thwav", timestamp, &sl->thwav);
	read_tensor3D("output/rec/soil_thiav", timestamp, &sl->thiav);

	read_map2D("output/rec/soil_Jinf", timestamp, 0, (void**)&sl->Jinf);
	read_map2D("output/rec/soil_Tv", timestamp, 0, (void**)&sl->Tv);
	read_map2D("output/rec/soil_Tmean", timestamp, 0, (void**)&sl->Tmean);
	read_map2D("output/rec/soil_thetaw_mean", timestamp, 0, (void**)&sl->thetaw_mean);
	read_map2D("output/rec/soil_thetai_mean", timestamp, 0, (void**)&sl->thetai_mean);
	read_map2D("output/rec/soil_psi_mean", timestamp, 0, (void**)&sl->psi_mean);
	read_map2D("output/rec/soil_Tmax", timestamp, 0, (void**)&sl->Tmax);
	read_map2D("output/rec/soil_thetaw_max", timestamp, 0, (void**)&sl->thetaw_max);
	read_map2D("output/rec/soil_thetai_max", timestamp, 0, (void**)&sl->thetai_max);
	read_map2D("output/rec/soil_Tmin", timestamp, 0, (void**)&sl->Tmin);
	read_map2D("output/rec/soil_thetaw_min", timestamp, 0, (void**)&sl->thetaw_min);
	read_map2D("output/rec/soil_thetai_min", timestamp, 0, (void**)&sl->thetai_min);
	read_double("output/rec/soil_TsupN", timestamp, &sl->TsupN);

	read_map2D("output/rec/soil_bc", timestamp, 2, (void**)&sl->bc);
	read_tensor3D("output/rec/soil_ET", timestamp, &sl->ET);
	/* FINISHED READING SOIL DATA */

	/* 2. READING WATER DATA */
	// The following tensor has a depth of Nl
	read_tensor3D("output/rec/water_q_sub", timestamp, &wat->q_sub);

	read_map2D("output/rec/water_weights_Kriging", timestamp, 0, (void**)&wat->weights_Kriging);
	read_map2D("output/rec/water_q_sup", timestamp, 0, (void**)&wat->q_sup);
	read_map2D("output/rec/water_h_sup", timestamp, 0, (void**)&wat->h_sup);
	read_map2D("output/rec/water_total", timestamp, 0, (void**)&wat->total);
	read_map2D("output/rec/water_Pn", timestamp, 0, (void**)&wat->Pn);
	read_map2D("output/rec/water_wcan_rain", timestamp, 0, (void**)&wat->wcan_rain);
	read_map2D("output/rec/water_wcan_snow", timestamp, 0, (void**)&wat->wcan_snow);
	read_map2D("output/rec/water_PrTOT_mean", timestamp, 0, (void**)&wat->PrTOT_mean);
	read_map2D("output/rec/water_PrSNW_mean", timestamp, 0, (void**)&wat->PrSNW_mean);
	read_map2D("output/rec/water_out1", timestamp, 0, (void**)&wat->out1);

	read_vector1D("output/rec/water_out2", timestamp, 0, (void**)&wat->out2);

	read_map2D("output/rec/water_hsupav", timestamp, 0, (void**)&wat->hsupav);
	read_map2D("output/rec/water_outfluxes", timestamp, 0, (void**)&wat->outfluxes);
	read_map2D("output/rec/water_error", timestamp, 0, (void**)&wat->error);
	/* FINISHED READING WATER DATA */

	/* 3. READING SNOW DATA */
	read_map2D("output/rec/snow_type", timestamp, 2, (void**)&snow->type);
	read_map2D("output/rec/snow_lnum", timestamp, 1, (void**)&snow->lnum);

	read_tensor3D("output/rec/snow_Dzl", timestamp, &snow->Dzl);
	read_tensor3D("output/rec/snow_w_liq", timestamp, &snow->w_liq);
	read_tensor3D("output/rec/snow_w_ice", timestamp, &snow->w_ice);
	read_tensor3D("output/rec/snow_T", timestamp, &snow->T);

	read_map2D("output/rec/snow_nondimens_age", timestamp, 0, (void**)&snow->nondimens_age);
	read_map2D("output/rec/snow_dimens_age", timestamp, 0, (void**)&snow->dimens_age);

	read_double("output/rec/snow_evap_basin", timestamp, &snow->evap_basin);
	read_double("output/rec/snow_subl_basin", timestamp, &snow->subl_basin);
	read_double("output/rec/snow_melted_basin", timestamp, &snow->melted_basin);

	read_vector1D("output/rec/snow_evap", timestamp, 0, (void**)&snow->evap);
	read_vector1D("output/rec/snow_subl", timestamp, 0, (void**)&snow->subl);
	read_vector1D("output/rec/snow_melted", timestamp, 0, (void**)&snow->melted);

	read_map2D("output/rec/snow_max", timestamp, 0, (void**)&snow->max);
	read_map2D("output/rec/snow_average", timestamp, 0, (void**)&snow->average);
	read_map2D("output/rec/snow_MELTED", timestamp, 0, (void**)&snow->MELTED);
	read_map2D("output/rec/snow_SUBL", timestamp, 0, (void**)&snow->SUBL);
	read_map2D("output/rec/snow_t_snow", timestamp, 0, (void**)&snow->t_snow);
	read_map2D("output/rec/snow_totav_snow", timestamp, 0, (void**)&snow->totav_snow);
	read_map2D("output/rec/snow_DDF", timestamp, 0, (void**)&snow->DDF);
	read_map2D("output/rec/snow_DDF1", timestamp, 0, (void**)&snow->DDF1);
	read_map2D("output/rec/snow_DDFvar", timestamp, 0, (void**)&snow->DDFvar);
	read_map2D("output/rec/snow_DDFcont", timestamp, 1, (void**)&snow->DDFcont);
	read_map2D("output/rec/snow_DDFTmin", timestamp, 0, (void**)&snow->DDFTmin);
	read_map2D("output/rec/snow_DDFmeltTL0", timestamp, 0, (void**)&snow->DDFmeltTL0);
	read_map2D("output/rec/snow_DDFmelt", timestamp, 0, (void**)&snow->DDFmelt);
	read_map2D("output/rec/snow_rho_newsnow", timestamp, 0, (void**)&snow->rho_newsnow);
	read_map2D("output/rec/snow_Qsub", timestamp, 0, (void**)&snow->Qsub);
	read_map2D("output/rec/snow_Wtrans", timestamp, 0, (void**)&snow->Wtrans);
	read_map2D("output/rec/snow_Qtrans", timestamp, 0, (void**)&snow->Qtrans);
	read_map2D("output/rec/snow_Qtrans_x", timestamp, 0, (void**)&snow->Qtrans_x);

	read_map2D("output/rec/snow_Qtrans_y", timestamp, 0, (void**)&snow->Qtrans_y);
	read_map2D("output/rec/snow_Wtot", timestamp, 0, (void**)&snow->Wtot);
	read_map2D("output/rec/snow_Wsubl_cum", timestamp, 0, (void**)&snow->Wsubl_cum);
	read_map2D("output/rec/snow_Wsusp_cum", timestamp, 0, (void**)&snow->Wsusp_cum);
	read_map2D("output/rec/snow_Wsalt_cum", timestamp, 0, (void**)&snow->Wsalt_cum);
	read_map2D("output/rec/snow_Wsubgrid_cum", timestamp, 0, (void**)&snow->Wsubgrid_cum);
	read_map2D("output/rec/snow_out_bs", timestamp, 0, (void**)&snow->out_bs);
	read_map2D("output/rec/snow_ListonSWE", timestamp, 0, (void**)&snow->ListonSWE);
	read_map2D("output/rec/snow_softSWE", timestamp, 0, (void**)&snow->softSWE);
	read_map2D("output/rec/snow_softSWE1", timestamp, 0, (void**)&snow->softSWE1);
	read_map2D("output/rec/snow_Dplot", timestamp, 0, (void**)&snow->Dplot);

	read_vector1D("output/rec/snow_CR1", timestamp, 0, (void**)&snow->CR1);
	read_vector1D("output/rec/snow_CR2", timestamp, 0, (void**)&snow->CR2);
	read_vector1D("output/rec/snow_CR3", timestamp, 0, (void**)&snow->CR3);
	read_vector1D("output/rec/snow_CR1m", timestamp, 0, (void**)&snow->CR1m);
	read_vector1D("output/rec/snow_CR2m", timestamp, 0, (void**)&snow->CR2m);
	read_vector1D("output/rec/snow_CR3m", timestamp, 0, (void**)&snow->CR3m);
	read_vector1D("output/rec/snow_change_dir_wind", timestamp, 1, (void**)&snow->change_dir_wind);

	read_map2D("output/rec/snow_Psnow", timestamp, 0, (void**)&snow->Psnow);
	read_map2D("output/rec/snow_rhoSOFT", timestamp, 0, (void**)&snow->rhoSOFT);
	/* FINISHED READING SNOW DATA */
	
	/* 4. READING GLACIER DATA */
	read_map2D("output/rec/glacier_lnum", timestamp, 1, (void**)&glac->lnum);

	read_tensor3D("output/rec/glacier_Dzl", timestamp, &glac->Dzl);
	read_tensor3D("output/rec/glacier_w_liq", timestamp, &glac->w_liq);
	read_tensor3D("output/rec/glacier_w_ice", timestamp, &glac->w_ice);
	read_tensor3D("output/rec/glacier_T", timestamp, &glac->T);

	read_double("output/rec/glacier_evap_basin", timestamp, &glac->evap_basin);
	read_double("output/rec/glacier_subl_basin", timestamp, &glac->subl_basin);
	read_double("output/rec/glacier_melted_basin", timestamp, &glac->melted_basin);

	read_vector1D("output/rec/glacier_evap", timestamp, 0, (void**)&glac->evap);
	read_vector1D("output/rec/glacier_subl", timestamp, 0, (void**)&glac->subl);
	read_vector1D("output/rec/glacier_melted", timestamp, 0, (void**)&glac->melted);

	read_map2D("output/rec/glacier_MELTED", timestamp, 0, (void**)&glac->MELTED);
	read_map2D("output/rec/glacier_SUBL", timestamp, 0, (void**)&glac->SUBL);
	read_map2D("output/rec/glacier_DDF", timestamp, 0, (void**)&glac->DDF);
	read_map2D("output/rec/glacier_DDF1", timestamp, 0, (void**)&glac->DDF1);
	read_map2D("output/rec/glacier_DDFvar", timestamp, 0, (void**)&glac->DDFvar);
	read_map2D("output/rec/glacier_DDFcont", timestamp, 1, (void**)&glac->DDFcont);
	read_map2D("output/rec/glacier_DDFTmin", timestamp, 0, (void**)&glac->DDFTmin);
	read_map2D("output/rec/glacier_DDFmeltTL0", timestamp, 0, (void**)&glac->DDFmeltTL0);
	read_map2D("output/rec/glacier_DDFmelt", timestamp, 0, (void**)&glac->DDFmelt);
	/* FINISHED READING GLACIER DATA */
	
	/* 5. READING CHANNEL DATA */
	read_vector1D("output/rec/channel_r", timestamp, 1, (void**)&cnet->r);
	read_vector1D("output/rec/channel_c", timestamp, 1, (void**)&cnet->c);
	read_map2D("output/rec/channel_ch", timestamp, 1, (void**)&cnet->ch);
	read_vector1D("output/rec/channel_Q", timestamp, 0, (void**)&cnet->Q);
	read_vector1D("output/rec/channel_s0", timestamp, 0, (void**)&cnet->s0);

	read_map2D("output/rec/channel_fraction_spread", timestamp, 0, (void**)&cnet->fraction_spread);

	read_vector1D("output/rec/channel_Q_sup_s", timestamp, 0, (void**)&cnet->Q_sup_s);
	read_vector1D("output/rec/channel_Q_sub_s", timestamp, 0, (void**)&cnet->Q_sub_s);
	read_vector1D("output/rec/channel_Qsup_spread", timestamp, 0, (void**)&cnet->Qsup_spread);
	read_vector1D("output/rec/channel_Qsub_spread", timestamp, 0, (void**)&cnet->Qsub_spread);
	read_vector1D("output/rec/channel_Qsup", timestamp, 0, (void**)&cnet->Qsup);
	read_vector1D("output/rec/channel_Qsub", timestamp, 0, (void**)&cnet->Qsub);
	/* FINISHED READING CHANNEL DATA */

	/* 5. READING ENERGY DATA */
	read_map2D("output/rec/energy_Rn_mean", timestamp, 0, (void**)&egy->Rn_mean);
	read_map2D("output/rec/energy_Rn_max", timestamp, 0, (void**)&egy->Rn_max);
	read_map2D("output/rec/energy_Rn_min", timestamp, 0, (void**)&egy->Rn_min);
	read_map2D("output/rec/energy_LW_in", timestamp, 0, (void**)&egy->LW_in);
	read_map2D("output/rec/energy_LW_out", timestamp, 0, (void**)&egy->LW_out);
	read_map2D("output/rec/energy_LW_max", timestamp, 0, (void**)&egy->LW_max);
	read_map2D("output/rec/energy_LW_min", timestamp, 0, (void**)&egy->LW_min);
	read_map2D("output/rec/energy_SW", timestamp, 0, (void**)&egy->SW);
	read_map2D("output/rec/energy_SW_max", timestamp, 0, (void**)&egy->SW_max);
	read_map2D("output/rec/energy_ET_mean", timestamp, 0, (void**)&egy->ET_mean);
	read_map2D("output/rec/energy_ET_max", timestamp, 0, (void**)&egy->ET_max);
	read_map2D("output/rec/energy_ET_min", timestamp, 0, (void**)&egy->ET_min);
	read_map2D("output/rec/energy_H_mean", timestamp, 0, (void**)&egy->H_mean);
	read_map2D("output/rec/energy_H_max", timestamp, 0, (void**)&egy->H_max);
	read_map2D("output/rec/energy_H_min", timestamp, 0, (void**)&egy->H_min);
	read_map2D("output/rec/energy_G_mean", timestamp, 0, (void**)&egy->G_mean);
	read_map2D("output/rec/energy_G_max", timestamp, 0, (void**)&egy->G_max);
	read_map2D("output/rec/energy_G_min", timestamp, 0, (void**)&egy->G_min);
	read_map2D("output/rec/energy_G_snowsoil", timestamp, 0, (void**)&egy->G_snowsoil);
	read_map2D("output/rec/energy_Ts_mean", timestamp, 0, (void**)&egy->Ts_mean);
	read_map2D("output/rec/energy_Ts_max", timestamp, 0, (void**)&egy->Ts_max);
	read_map2D("output/rec/energy_Ts_min", timestamp, 0, (void**)&egy->Ts_min);
	read_map2D("output/rec/energy_Rswdown_mean", timestamp, 0,(void**) &egy->Rswdown_mean);
	read_map2D("output/rec/energy_Rswdown_max", timestamp, 0, (void**)&egy->Rswdown_max);
	read_map2D("output/rec/energy_Ta_mean", timestamp, 0, (void**)&egy->Ta_mean);
	read_map2D("output/rec/energy_Ta_max", timestamp, 0, (void**)&egy->Ta_max);
	read_map2D("output/rec/energy_Ta_min", timestamp, 0, (void**)&egy->Ta_min);
	read_map2D("output/rec/energy_out1", timestamp, 0, (void**)&egy->out1);

	read_vector1D("output/rec/energy_out2", timestamp, 0, (void**)&egy->out2);

	read_map2D("output/rec/energy_out3", timestamp, 0, (void**)&egy->out3);
	read_map2D("output/rec/energy_Rswbeam", timestamp, 0, (void**)&egy->Rswbeam);
	read_map2D("output/rec/energy_nDt_shadow", timestamp, 1, (void**)&egy->nDt_shadow);
	read_map2D("output/rec/energy_nDt_sun", timestamp, 1, (void**)&egy->nDt_sun);
	read_map2D("output/rec/energy_Hgplot", timestamp, 0, (void**)&egy->Hgplot);
	read_map2D("output/rec/energy_LEgplot", timestamp, 0, (void**)&egy->LEgplot);
	read_map2D("output/rec/energy_Hvplot", timestamp, 0, (void**)&egy->Hvplot);
	read_map2D("output/rec/energy_LEvplot", timestamp, 0, (void**)&egy->LEvplot);
	read_map2D("output/rec/energy_SWinplot", timestamp, 0, (void**)&egy->SWinplot);
	read_map2D("output/rec/energy_SWgplot", timestamp, 0, (void**)&egy->SWgplot);
	read_map2D("output/rec/energy_SWvplot", timestamp, 0, (void**)&egy->SWvplot);
	read_map2D("output/rec/energy_LWinplot", timestamp, 0, (void**)&egy->LWinplot);
	read_map2D("output/rec/energy_LWgplot", timestamp, 0, (void**)&egy->LWgplot);
	read_map2D("output/rec/energy_LWvplot", timestamp, 0, (void**)&egy->LWvplot);
	read_map2D("output/rec/energy_Tgplot", timestamp, 0, (void**)&egy->Tgplot);
	read_map2D("output/rec/energy_Tvplot", timestamp, 0, (void**)&egy->Tvplot);
	read_map2D("output/rec/energy_Tsplot", timestamp, 0, (void**)&egy->Tsplot);
	read_map2D("output/rec/energy_SWin", timestamp, 0, (void**)&egy->SWin);
	read_map2D("output/rec/energy_Hgrid", timestamp, 0, (void**)&egy->Hgrid);
	read_map2D("output/rec/energy_Tsgrid", timestamp, 0, (void**)&egy->Tsgrid);

	read_vector1D("output/rec/energy_VSFA", timestamp, 0, (void**)&egy->VSFA);
	read_vector1D("output/rec/energy_HSFA", timestamp, 0, (void**)&egy->HSFA);

	read_double("output/rec/energy_hsun", timestamp, &egy->hsun);
	read_double("output/rec/energy_dsun", timestamp, &egy->dsun);
	/* FINISHED READING ENERGY DATA */

	/* 6. READ TIMES DATA */
	read_short("output/rec/times_iter", timestamp, &times->iter);
	read_short("output/rec/times_n_iter", timestamp, &times->n_iter);
	read_double("output/rec/times_TH", timestamp, &times->TH);

	read_long("output/rec/times_i_pixel", timestamp, &times->i_pixel);
	read_long("output/rec/times_n_pixel", timestamp, &times->n_pixel);
	read_long("output/rec/times_i_basin", timestamp, &times->i_basin);
	read_long("output/rec/times_n_basin", timestamp, &times->n_basin);
	read_long("output/rec/times_i_plot", timestamp, &times->i_plot);
	read_long("output/rec/times_n_plot", timestamp, &times->n_plot);
	read_long("output/rec/times_nt_plot", timestamp, &times->nt_plot);
	read_long("output/rec/times_d_plot", timestamp, &times->d_plot);

	read_double("output/rec/times_JD", timestamp, &times->JD);

	read_long("output/rec/times_DD", timestamp, &times->DD);
	read_long("output/rec/times_MM", timestamp, &times->MM);
	read_long("output/rec/times_AAAA", timestamp, &times->AAAA);
	read_long("output/rec/times_hh", timestamp, &times->hh);
	read_long("output/rec/times_mm", timestamp, &times->mm);

	read_double("output/rec/times_time", timestamp, &times->time);
	read_double("output/rec/times_egy", timestamp, &times->egy);
	read_double("output/rec/times_vert_wb", timestamp, &times->vert_wb);
	read_double("output/rec/times_horiz_wb", timestamp, &times->horiz_wb);
	read_double("output/rec/times_writeout", timestamp, &times->writeout);

	read_short("output/rec/times_count_super_printed", timestamp, &times->count_super_printed);
	/* FINISHED READING TIMES DATA */

	read_output_variables(timestamp);

	printf("Read %u files for timestamp %s\n\n", filecounter, timestamp);
	if (filecounter<10){
		fprintf(stderr, "\nNot enough data to recover simulation! Exiting\n\n");
		exit(1);
	}
}

void write_recovery_files(char *timestamp, TIMES *times, WATER *wat, CHANNEL *cnet, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac)
{
	char tmp_fname[128];

	/*****************************************************************************/
	/* Actual writing of recovery files (and constructing the filenames) follows */
	/*****************************************************************************/
	/* 1. WRITING SOIL DATA */
	write_map2D("output/rec/soil_type", timestamp, 2, sl->type);

	// The following tensors have a depth of Nl (except sl->pa which has a depth of 1)
	write_tensorseries(1, timestamp, "output/rec/soil_pa", 0, 2, sl->pa, UV);
	write_tensorseries(1, timestamp, "output/rec/soil_P", 0, 2, sl->P, UV);
	write_tensorseries(1, timestamp, "output/rec/soil_T", 0, 2, sl->T, UV);
	write_tensorseries(1, timestamp, "output/rec/soil_thice", 0, 2, sl->thice, UV);
	write_tensorseries(1, timestamp, "output/rec/soil_J", 0, 2, sl->J, UV);

	/* The following DOUBLETENSORs seem to be unused or NULL */
	write_tensorseries(1, timestamp, "output/rec/soil_Tav", 0, 2, sl->Tav, UV);
	write_tensorseries(1, timestamp, "output/rec/soil_thwav", 0, 2, sl->thwav, UV);
	write_tensorseries(1, timestamp, "output/rec/soil_thiav", 0, 2, sl->thiav, UV);
		 

	write_map2D("output/rec/soil_Jinf", timestamp, 0, sl->Jinf);
	write_map2D("output/rec/soil_Tv", timestamp, 0, sl->Tv);
	write_map2D("output/rec/soil_Tmean", timestamp, 0, sl->Tmean);
	write_map2D("output/rec/soil_thetaw_mean", timestamp, 0, sl->thetaw_mean);
	write_map2D("output/rec/soil_thetai_mean", timestamp, 0, sl->thetai_mean);
	write_map2D("output/rec/soil_psi_mean", timestamp, 0, sl->psi_mean);
	write_map2D("output/rec/soil_Tmax", timestamp, 0, sl->Tmax);
	write_map2D("output/rec/soil_thetaw_max", timestamp, 0, sl->thetaw_max);
	write_map2D("output/rec/soil_thetai_max", timestamp, 0, sl->thetai_max);
	write_map2D("output/rec/soil_Tmin", timestamp, 0, sl->Tmin);
	write_map2D("output/rec/soil_thetaw_min", timestamp, 0, sl->thetaw_min);
	write_map2D("output/rec/soil_thetai_min", timestamp, 0, sl->thetai_min);

	sprintf(tmp_fname, "%s_%s", "output/rec/soil_TsupN", timestamp); 
	write_double_to_file(tmp_fname, sl->TsupN);

	write_map2D("output/rec/soil_bc", timestamp, 2, sl->bc);

	write_tensorseries(1, timestamp, "output/rec/soil_ET", 0, 2, sl->ET, UV);
	/* FINISHED WRITING SOIL DATA */
	
	/* 2. WRITING WATER DATA */
	// The following tensor has a depth of Nl
	write_tensorseries(1, timestamp, "output/rec/water_q_sub", 0, 2, wat->q_sub, UV);
	
	write_map2D("output/rec/water_weights_Kriging", timestamp, 0, wat->weights_Kriging);
	write_map2D("output/rec/water_q_sup", timestamp, 0, wat->q_sup);
	write_map2D("output/rec/water_h_sup", timestamp, 0, wat->h_sup);
	write_map2D("output/rec/water_total", timestamp, 0, wat->total);
	write_map2D("output/rec/water_Pn", timestamp, 0, wat->Pn);
	write_map2D("output/rec/water_wcan_rain", timestamp, 0, wat->wcan_rain);
	write_map2D("output/rec/water_wcan_snow", timestamp, 0, wat->wcan_snow);
	write_map2D("output/rec/water_PrTOT_mean", timestamp, 0, wat->PrTOT_mean);
	write_map2D("output/rec/water_PrSNW_mean", timestamp, 0, wat->PrSNW_mean);

	write_map2D("output/rec/water_out1", timestamp, 0, wat->out1);
	sprintf(tmp_fname, "%s_%s", "output/rec/water_out2", timestamp); 
	write_vector_to_file(tmp_fname, 0, wat->out2);

	write_map2D("output/rec/water_hsupav", timestamp, 0, wat->hsupav);
	write_map2D("output/rec/water_outfluxes", timestamp, 0, wat->outfluxes);
	write_map2D("output/rec/water_error", timestamp, 0, wat->error);
	/* FINISHED WRITING WATER DATA */

	/* 3. WRITING SNOW DATA */
	write_map2D("output/rec/snow_type", timestamp, 2, snow->type);
	write_map2D("output/rec/snow_lnum", timestamp, 1, snow->lnum);
	
	// The following tensors have a depth of par->snowlayer_max		 
	write_tensorseries(1, timestamp, "output/rec/snow_Dzl", 0, 2, snow->Dzl, UV);
	write_tensorseries(1, timestamp, "output/rec/snow_w_liq", 0, 2, snow->w_liq, UV);
	write_tensorseries(1, timestamp, "output/rec/snow_w_ice", 0, 2, snow->w_ice, UV);
	write_tensorseries(1, timestamp, "output/rec/snow_T", 0, 2, snow->T, UV);

	write_map2D("output/rec/snow_nondimens_age", timestamp, 0, snow->nondimens_age);
	write_map2D("output/rec/snow_dimens_age", timestamp, 0, snow->dimens_age);

	sprintf(tmp_fname, "%s_%s", "output/rec/snow_evap_basin", timestamp); 
	write_double_to_file(tmp_fname, snow->evap_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_subl_basin", timestamp); 
	write_double_to_file(tmp_fname, snow->subl_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_melted_basin", timestamp); 
	write_double_to_file(tmp_fname, snow->melted_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_evap", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->evap);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_subl", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->subl);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_melted", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->melted);

	write_map2D("output/rec/snow_max", timestamp, 0, snow->max);
	write_map2D("output/rec/snow_average", timestamp, 0, snow->average);
	write_map2D("output/rec/snow_MELTED", timestamp, 0, snow->MELTED);
	write_map2D("output/rec/snow_SUBL", timestamp, 0, snow->SUBL);
	write_map2D("output/rec/snow_t_snow", timestamp, 0, snow->t_snow);
	write_map2D("output/rec/snow_totav_snow", timestamp, 0, snow->totav_snow);
	write_map2D("output/rec/snow_DDF", timestamp, 0, snow->DDF);
	write_map2D("output/rec/snow_DDF1", timestamp, 0, snow->DDF1);
	write_map2D("output/rec/snow_DDFvar", timestamp, 0, snow->DDFvar);
	write_map2D("output/rec/snow_DDFcont", timestamp, 1, snow->DDFcont);

	write_map2D("output/rec/snow_DDFTmin", timestamp, 0, snow->DDFTmin);
	write_map2D("output/rec/snow_DDFmeltTL0", timestamp, 0, snow->DDFmeltTL0);
	write_map2D("output/rec/snow_DDFmelt", timestamp, 0, snow->DDFmelt);
	write_map2D("output/rec/snow_rho_newsnow", timestamp, 0, snow->rho_newsnow);
	write_map2D("output/rec/snow_Qsub", timestamp, 0, snow->Qsub);
	write_map2D("output/rec/snow_Wtrans", timestamp, 0, snow->Wtrans);
	write_map2D("output/rec/snow_Qtrans", timestamp, 0, snow->Qtrans);
	write_map2D("output/rec/snow_Qtrans_x", timestamp, 0, snow->Qtrans_x);
	write_map2D("output/rec/snow_Qtrans_y", timestamp, 0, snow->Qtrans_y);
	write_map2D("output/rec/snow_Wtot", timestamp, 0, snow->Wtot);
	write_map2D("output/rec/snow_Wsubl_cum", timestamp, 0, snow->Wsubl_cum);
	write_map2D("output/rec/snow_Wsusp_cum", timestamp, 0, snow->Wsusp_cum);
	write_map2D("output/rec/snow_Wsalt_cum", timestamp, 0, snow->Wsalt_cum);
	write_map2D("output/rec/snow_Wsubgrid_cum", timestamp, 0, snow->Wsubgrid_cum);
	write_map2D("output/rec/snow_out_bs", timestamp, 0, snow->out_bs);
	write_map2D("output/rec/snow_ListonSWE", timestamp, 0, snow->ListonSWE);
	write_map2D("output/rec/snow_softSWE", timestamp, 0, snow->softSWE);
	write_map2D("output/rec/snow_softSWE1", timestamp, 0, snow->softSWE1);
	write_map2D("output/rec/snow_Dplot", timestamp, 0, snow->Dplot);

	sprintf(tmp_fname, "%s_%s", "output/rec/snow_CR1", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->CR1);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_CR2", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->CR2);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_CR3", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->CR3);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_CR1m", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->CR1m);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_CR2m", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->CR2m);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_CR3m", timestamp); 
	write_vector_to_file(tmp_fname, 0, snow->CR3m);
	sprintf(tmp_fname, "%s_%s", "output/rec/snow_change_dir_wind", timestamp); 
	write_vector_to_file(tmp_fname, 1, snow->change_dir_wind);

	write_map2D("output/rec/snow_Psnow", timestamp, 0, snow->Psnow);
	write_map2D("output/rec/snow_rhoSOFT", timestamp, 0, snow->rhoSOFT);
	/* FINISHED WRITING SNOW DATA */

	/* 4. WRITING GLACIER DATA */
	write_map2D("output/rec/glacier_lnum", timestamp, 1, glac->lnum);
	//CONTINUE HERE
	write_tensorseries(1, timestamp, "output/rec/glacier_Dzl", 0, 2, glac->Dzl, UV);
	write_tensorseries(1, timestamp, "output/rec/glacier_w_liq", 0, 2, glac->w_liq, UV);
	write_tensorseries(1, timestamp, "output/rec/glacier_w_ice", 0, 2, glac->w_ice, UV);
	write_tensorseries(1, timestamp, "output/rec/glacier_T", 0, 2, glac->T, UV);

	sprintf(tmp_fname, "%s_%s", "output/rec/glacier_evap_basin", timestamp); 
	write_double_to_file(tmp_fname, glac->evap_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/glacier_subl_basin", timestamp); 
	write_double_to_file(tmp_fname, glac->subl_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/glacier_melted_basin", timestamp); 
	write_double_to_file(tmp_fname, glac->melted_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/glacier_evap", timestamp); 
	write_vector_to_file(tmp_fname, 0, glac->evap);
	sprintf(tmp_fname, "%s_%s", "output/rec/glacier_subl", timestamp); 
	write_vector_to_file(tmp_fname, 0, glac->subl);
	sprintf(tmp_fname, "%s_%s", "output/rec/glacier_melted", timestamp); 
	write_vector_to_file(tmp_fname, 0, glac->melted);
		 
	write_map2D("output/rec/glacier_MELTED", timestamp, 0, glac->MELTED);
	write_map2D("output/rec/glacier_SUBL", timestamp, 0, glac->SUBL);
	write_map2D("output/rec/glacier_DDF", timestamp, 0, glac->DDF);
	write_map2D("output/rec/glacier_DDF1", timestamp, 0, glac->DDF1);
	write_map2D("output/rec/glacier_DDFvar", timestamp, 0, glac->DDFvar);
	write_map2D("output/rec/glacier_DDFcont", timestamp, 1, glac->DDFcont);
	write_map2D("output/rec/glacier_DDFTmin", timestamp, 0, glac->DDFTmin);
	write_map2D("output/rec/glacier_DDFmeltTL0", timestamp, 0, glac->DDFmeltTL0);
	write_map2D("output/rec/glacier_DDFmelt", timestamp, 0, glac->DDFmelt);
	/* FINISHED WRITING GLACIER DATA */

	/* 5. WRITING CHANNEL DATA */
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_r", timestamp);
	write_vector_to_file(tmp_fname, 1, cnet->r);
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_c", timestamp);
	write_vector_to_file(tmp_fname, 1, cnet->c);

	write_map2D("output/rec/channel_ch", timestamp, 1, cnet->ch);

	sprintf(tmp_fname, "%s_%s", "output/rec/channel_Q", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->Q);
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_s0", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->s0);

	write_map2D("output/rec/channel_fraction_spread", timestamp, 0, cnet->fraction_spread);

	sprintf(tmp_fname, "%s_%s", "output/rec/channel_Q_sup_s", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->Q_sup_s);
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_Q_sub_s", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->Q_sub_s);
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_Qsup_spread", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->Qsup_spread);
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_Qsub_spread", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->Qsub_spread);
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_Qsup", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->Qsup);
	sprintf(tmp_fname, "%s_%s", "output/rec/channel_Qsub", timestamp);
	write_vector_to_file(tmp_fname, 0, cnet->Qsub);
	/* FINISHED WRITING CHANNEL DATA */

	/* 5. WRITING ENERGY DATA */
	write_map2D("output/rec/energy_Rn_mean", timestamp, 0, egy->Rn_mean);
	write_map2D("output/rec/energy_Rn_max", timestamp, 0, egy->Rn_max);
	write_map2D("output/rec/energy_Rn_min", timestamp, 0, egy->Rn_min);
	write_map2D("output/rec/energy_LW_in", timestamp, 0, egy->LW_in);
	write_map2D("output/rec/energy_LW_out", timestamp, 0, egy->LW_out);
	write_map2D("output/rec/energy_LW_max", timestamp, 0, egy->LW_max);
	write_map2D("output/rec/energy_LW_min", timestamp, 0, egy->LW_min);
	write_map2D("output/rec/energy_SW", timestamp, 0, egy->SW);
	write_map2D("output/rec/energy_SW_max", timestamp, 0, egy->SW_max);
	write_map2D("output/rec/energy_ET_mean", timestamp, 0, egy->ET_mean);
	write_map2D("output/rec/energy_ET_max", timestamp, 0, egy->ET_max);
	write_map2D("output/rec/energy_ET_min", timestamp, 0, egy->ET_min);
	write_map2D("output/rec/energy_H_mean", timestamp, 0, egy->H_mean);
	write_map2D("output/rec/energy_H_max", timestamp, 0, egy->H_max);
	write_map2D("output/rec/energy_H_min", timestamp, 0, egy->H_min);

	write_map2D("output/rec/energy_G_mean", timestamp, 0, egy->G_mean);
	write_map2D("output/rec/energy_G_max", timestamp, 0, egy->G_max);
	write_map2D("output/rec/energy_G_min", timestamp, 0, egy->G_min);
	write_map2D("output/rec/energy_G_snowsoil", timestamp, 0, egy->G_snowsoil);
	write_map2D("output/rec/energy_Ts_mean", timestamp, 0, egy->Ts_mean);
	write_map2D("output/rec/energy_Ts_max", timestamp, 0, egy->Ts_max);
	write_map2D("output/rec/energy_Ts_min", timestamp, 0, egy->Ts_min);
	write_map2D("output/rec/energy_Rswdown_mean", timestamp, 0, egy->Rswdown_mean);
	write_map2D("output/rec/energy_Rswdown_max", timestamp, 0, egy->Rswdown_max);
	write_map2D("output/rec/energy_Ta_mean", timestamp, 0, egy->Ta_mean);
	write_map2D("output/rec/energy_Ta_max", timestamp, 0, egy->Ta_max);
	write_map2D("output/rec/energy_Ta_min", timestamp, 0, egy->Ta_min);
	write_map2D("output/rec/energy_out1", timestamp, 0, egy->out1);

	sprintf(tmp_fname, "%s_%s", "output/rec/energy_out2", timestamp);
	write_vector_to_file(tmp_fname, 0, egy->out2);
	write_map2D("output/rec/energy_out3", timestamp, 0, egy->out3);
	
	write_map2D("output/rec/energy_Rswbeam", timestamp, 0, egy->Rswbeam);
	write_map2D("output/rec/energy_nDt_shadow", timestamp, 1, egy->nDt_shadow);
	write_map2D("output/rec/energy_nDt_sun", timestamp, 1, egy->nDt_sun);
	write_map2D("output/rec/energy_Hgplot", timestamp, 0, egy->Hgplot);
	write_map2D("output/rec/energy_LEgplot", timestamp, 0, egy->LEgplot);
	write_map2D("output/rec/energy_Hvplot", timestamp, 0, egy->Hvplot);
	write_map2D("output/rec/energy_LEvplot", timestamp, 0, egy->LEvplot);
	write_map2D("output/rec/energy_SWinplot", timestamp, 0, egy->SWinplot);
	write_map2D("output/rec/energy_SWgplot", timestamp, 0, egy->SWgplot);
	write_map2D("output/rec/energy_SWvplot", timestamp, 0, egy->SWvplot);
	write_map2D("output/rec/energy_LWinplot", timestamp, 0, egy->LWinplot);
	write_map2D("output/rec/energy_LWgplot", timestamp, 0, egy->LWgplot);
	write_map2D("output/rec/energy_LWvplot", timestamp, 0, egy->LWvplot);
	write_map2D("output/rec/energy_Tgplot", timestamp, 0, egy->Tgplot);
	write_map2D("output/rec/energy_Tvplot", timestamp, 0, egy->Tvplot);
	write_map2D("output/rec/energy_Tsplot", timestamp, 0, egy->Tsplot);
	write_map2D("output/rec/energy_SWin", timestamp, 0, egy->SWin);
	write_map2D("output/rec/energy_Hgrid", timestamp, 0, egy->Hgrid);
	write_map2D("output/rec/energy_Tsgrid", timestamp, 0, egy->Tsgrid);

	sprintf(tmp_fname, "%s_%s", "output/rec/energy_VSFA", timestamp);
	write_vector_to_file(tmp_fname, 0, egy->VSFA);
	sprintf(tmp_fname, "%s_%s", "output/rec/energy_HSFA", timestamp);
	write_vector_to_file(tmp_fname, 0, egy->HSFA);
	
	sprintf(tmp_fname, "%s_%s", "output/rec/energy_hsun", timestamp); 
	write_double_to_file(tmp_fname, egy->hsun);
	sprintf(tmp_fname, "%s_%s", "output/rec/energy_dsun", timestamp); 
	write_double_to_file(tmp_fname, egy->dsun);
	/* FINISHED WRITING ENERGY DATA */

	/* 6. WRITING TIMES DATA */
	sprintf(tmp_fname, "%s_%s", "output/rec/times_iter", timestamp); 
	write_short_to_file(tmp_fname, times->iter);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_n_iter", timestamp); 
	write_short_to_file(tmp_fname, times->n_iter);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_TH", timestamp); 
	write_double_to_file(tmp_fname, times->TH);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_i_pixel", timestamp); 
	write_long_to_file(tmp_fname, times->i_pixel);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_n_pixel", timestamp); 
	write_long_to_file(tmp_fname, times->n_pixel);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_i_basin", timestamp); 
	write_long_to_file(tmp_fname, times->i_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_n_basin", timestamp); 
	write_long_to_file(tmp_fname, times->n_basin);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_i_plot", timestamp); 
	write_long_to_file(tmp_fname, times->i_plot);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_n_plot", timestamp); 
	write_long_to_file(tmp_fname, times->n_plot);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_nt_plot", timestamp); 
	write_long_to_file(tmp_fname, times->nt_plot);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_d_plot", timestamp); 
	write_long_to_file(tmp_fname, times->d_plot);

	sprintf(tmp_fname, "%s_%s", "output/rec/times_JD", timestamp); 
	write_double_to_file(tmp_fname, times->JD);

	sprintf(tmp_fname, "%s_%s", "output/rec/times_DD", timestamp); 
	write_long_to_file(tmp_fname, times->DD);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_MM", timestamp); 
	write_long_to_file(tmp_fname, times->MM);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_AAAA", timestamp); 
	write_long_to_file(tmp_fname, times->AAAA);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_hh", timestamp); 
	write_long_to_file(tmp_fname, times->hh);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_mm", timestamp); 
	write_long_to_file(tmp_fname, times->mm);

	sprintf(tmp_fname, "%s_%s", "output/rec/times_time", timestamp); 
	write_double_to_file(tmp_fname, times->time);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_egy", timestamp); 
	write_double_to_file(tmp_fname, times->egy);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_vert_wb", timestamp); 
	write_double_to_file(tmp_fname, times->vert_wb);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_horiz_wb", timestamp); 
	write_double_to_file(tmp_fname, times->horiz_wb);
	sprintf(tmp_fname, "%s_%s", "output/rec/times_writeout", timestamp); 
	write_double_to_file(tmp_fname, times->writeout);

	sprintf(tmp_fname, "%s_%s", "output/rec/times_count_super_printed", timestamp); 
	write_short_to_file(tmp_fname, times->count_super_printed);
	/* FINISHED WRITING TIMES DATA */

	write_output_variables(timestamp);
}
