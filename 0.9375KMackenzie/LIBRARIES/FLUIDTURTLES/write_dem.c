#include "turtle.h"

#include "write_dem.h"

#include "t_utilities.h"

#include "t_datamanipulation.h"





/**-----------------------------------------------------------------------*/

/** Writes a matrix of short int to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void shortmatrix_dem(SHORTMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print)

{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

	outputname=join_strings(WORKING_DIRECTORY,outputname);

	outputfile=t_fopen(outputname,"w");

	fprintf(outputfile,"/** File %s */\n",outputname);

	fprintf(outputfile,"/**  %s */\n",comment);

	fprintf(outputfile,"index{3,DEM}\n");

	fprintf(outputfile,"1: float array pixels size  ");

	write_floatarray_elements(outputfile,U,600);

	fprintf(outputfile,"2: float array novalues ");

	write_floatarray_elements(outputfile,V,600);

	fprintf(outputfile,"3: short matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

	write_shortmatrix_elements(outputfile,matrix,matrix->nch+1);

	t_fclose(outputfile);

} else {

	return;

}

}


/**-----------------------------------------------------------------------*/

/** Writes a matrix of short int to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void shortmatrix_dem3(SHORTMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,

	char *outputname, char *comment,short print)

{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

	outputname=join_strings(WORKING_DIRECTORY,outputname);

	outputfile=t_fopen(outputname,"w");

	fprintf(outputfile,"/** File %s */\n",outputname);

	fprintf(outputfile,"/**  %s */\n",comment);

	fprintf(outputfile,"index{3,DEM}\n");

	fprintf(outputfile,"1: double array pixels size  ");

	write_doublearray_elements(outputfile,U,600);

	fprintf(outputfile,"2: double array novalues ");

	write_doublearray_elements(outputfile,V,600);

	fprintf(outputfile,"3: short matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

	write_shortmatrix_elements(outputfile,matrix,matrix->nch+1);

	t_fclose(outputfile);

} else {

	return;

}

}




/**-----------------------------------------------------------------------*/

/** Writes a matrix of long int to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void longmatrix_dem(LONGMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

outputname=join_strings(WORKING_DIRECTORY,outputname);

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/** File %s */\n",outputname);

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3,DEM}\n");

fprintf(outputfile,"1: float array pixels size  ");

write_floatarray_elements(outputfile,U,600);

fprintf(outputfile,"2: float array novalues ");

write_floatarray_elements(outputfile,V,600);

fprintf(outputfile,"3: long matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_longmatrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of int to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void intmatrix_dem(INTMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

outputname=join_strings(WORKING_DIRECTORY,outputname);

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/** File %s */\n",outputname);

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3,DEM}\n");

fprintf(outputfile,"1: float array pixels size  ");

write_floatarray_elements(outputfile,U,600);

fprintf(outputfile,"2: float array novalues ");

write_floatarray_elements(outputfile,V,600);

fprintf(outputfile,"3: int matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_intmatrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of float to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void floatmatrix_dem(FLOATMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

//char *newoname;

FILE *outputfile;

extern char *WORKING_DIRECTORY;

//printf("i1qu\n");



/* Writes output file */

if(print==1){

outputname=join_strings(WORKING_DIRECTORY,outputname);

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/** File %s */\n",outputname);

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3,DEM}\n");

fprintf(outputfile,"1: float array pixels size  ");

write_floatarray_elements(outputfile,U,600);

fprintf(outputfile,"2: float array novalues ");

write_floatarray_elements(outputfile,V,600);

fprintf(outputfile,"3: float matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_floatmatrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of double to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void doublematrix_dem(DOUBLEMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

outputname=join_strings(WORKING_DIRECTORY,outputname);

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/** File %s */\n",outputname);

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3,DEM}\n");

fprintf(outputfile,"1: float array pixels size  ");

write_floatarray_elements(outputfile,U,600);

fprintf(outputfile,"2: float array novalues ");

write_floatarray_elements(outputfile,V,600);

fprintf(outputfile,"3: double matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_doublematrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}

/**-----------------------------------------------------------------------*/

/** Writes a matrix of double to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void doublematrix_dem3(DOUBLEMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

outputname=join_strings(WORKING_DIRECTORY,outputname);

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/** File %s */\n",outputname);

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3,DEM}\n");

fprintf(outputfile,"1: double array pixels size  ");

write_doublearray_elements(outputfile,U,600);

fprintf(outputfile,"2: double array novalues ");

write_doublearray_elements(outputfile,V,600);

fprintf(outputfile,"3: double matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_doublematrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of double to a turtle file */



/* Stampa una matrice in un file in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void doublematrix_control(DOUBLEMATRIX *matrix,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

	outputname=join_strings(WORKING_DIRECTORY,outputname);

	outputfile=t_fopen(outputname,"w");

	fprintf(outputfile,"/** File %s */\n",outputname);

	fprintf(outputfile,"/**  %s */\n",comment);

	fprintf(outputfile,"index{1}\n");

	fprintf(outputfile,"1: double matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

	write_doublematrix_elements(outputfile,matrix,matrix->nch+1);

	t_fclose(outputfile);

} else {

	return;

}

}



/**-----------------------------------------------------------------------*/

/** Writes a vector of float to a turtle file */



/* Stampa un vettore in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void floatvector_dem(FLOATVECTOR *vector_dem,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

	outputname=join_strings(WORKING_DIRECTORY,outputname);

	outputfile=t_fopen(outputname,"w");

	fprintf(outputfile,"/** File %s */\n",outputname);

	fprintf(outputfile,"/**  %s */\n",comment);

	fprintf(outputfile,"index{1}\n");

	fprintf(outputfile,"1: float array vector  ");

	write_floatarray_elements(outputfile,vector_dem,vector_dem->nh+1);

	t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a vector of double to a turtle file */



/* Stampa un vettore in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void doublevector_dem(DOUBLEVECTOR *vector_dem,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;

extern char *WORKING_DIRECTORY;



/* Writes output file */

if(print==1){

	outputname=join_strings(WORKING_DIRECTORY,outputname);

	outputfile=t_fopen(outputname,"w");

	fprintf(outputfile,"/** File %s */\n",outputname);

	fprintf(outputfile,"/**  %s */\n",comment);

	fprintf(outputfile,"index{1}\n");

	fprintf(outputfile,"1:double array vector  ");

	write_doublearray_elements(outputfile,vector_dem,vector_dem->nh+1);

	t_fclose(outputfile);

} else {

	return;

}

}







/**-------------------------------------------------------------------------*/

/** Writes a layer of a 3D double tensor (i.e. a matrix) to a turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	WORKING_DIRECTORY deve essere dichiarata come variabile esterna nel main()

	Inputs:	tensor		tensor

	        layer		numero del livello del tensore che vuoi stampare

			U,V 		vettori di header del DEM

			outputname 	nome del file di output senza estensione

			comment		eventuali commenti da scrivere nel file */



void doubletensor_dem(DOUBLETENSOR *tensor,long layer,DOUBLEVECTOR *U,

                      DOUBLEVECTOR *V,char *outputname,char *comment,short print)



{

/* Declaration variables */

long r,c;

FILE *outputfile;

extern char *WORKING_DIRECTORY;

DOUBLEMATRIX *matrix;



/* Allocation and initialisation of matrix */

matrix=new_doublematrix(tensor->nrh,tensor->nch);

for(r=1;r<=tensor->nrh;r++){

   for(c=1;c<=tensor->nch;c++){

      matrix->co[r][c]=tensor->co[layer][r][c];

   }

}



/* Writes output file */

if(print==1){

outputname=join_strings(WORKING_DIRECTORY,outputname);

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/** File %s */\n",outputname);

fprintf(outputfile,"/**  %s (layer %d)*/\n",comment,layer);

fprintf(outputfile,"index{3,DEM}\n");

fprintf(outputfile,"1: float array pixels size  ");

write_doublearray_elements(outputfile,U,600);

fprintf(outputfile,"2: float array novalues ");

write_doublearray_elements(outputfile,V,600);

fprintf(outputfile,"3: double matrix layer_%d {%d,%d}\n",layer,matrix->nrh,

        matrix->nch);

write_doublematrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

free_doublematrix(matrix);

} else {

	return;

}

}



































/**-----------------------------------------------------------------------*/

/** Writes a matrix of short int to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

   !!! non viene usata la WORKING_DIRECTORY !!!

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void shortmatrix_dem2(SHORTMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,

	char *outputname, char *comment,short print)

{

/* Declaration variables */

FILE *outputfile;



/* Writes output file */

if(print==1){

	outputfile=t_fopen(outputname,"w");

	fprintf(outputfile,"/**  %s */\n",comment);

	fprintf(outputfile,"index{3}\n");

	fprintf(outputfile,"1: double array pixels size  ");

	write_doublearray_elements(outputfile,U,600);

	fprintf(outputfile,"2: double array novalues ");

	write_doublearray_elements(outputfile,V,600);

	fprintf(outputfile,"3: short matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

	write_shortmatrix_elements(outputfile,matrix,matrix->nch+1);

	t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of long int to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

   !!! non viene usata la WORKING_DIRECTORY !!!

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void longmatrix_dem2(LONGMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;



/* Writes output file */

if(print==1){

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3}\n");

fprintf(outputfile,"1: double array pixels size  ");

write_doublearray_elements(outputfile,U,600);

fprintf(outputfile,"2: double array novalues ");

write_doublearray_elements(outputfile,V,600);

fprintf(outputfile,"3: long matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_longmatrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of int to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

   !!! non viene usata la WORKING_DIRECTORY !!!

   Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void intmatrix_dem2(INTMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;



/* Writes output file */

if(print==1){

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3}\n");

fprintf(outputfile,"1: float array pixels size  ");

write_floatarray_elements(outputfile,U,600);

fprintf(outputfile,"2: float array novalues ");

write_floatarray_elements(outputfile,V,600);

fprintf(outputfile,"3: int matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_intmatrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of float to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

   !!! non viene usata la WORKING_DIRECTORY !!!

   Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void floatmatrix_dem2(FLOATMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;



/* Writes output file */

if(print==1){

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3}\n");

fprintf(outputfile,"1: float array pixels size  ");



write_floatarray_elements(outputfile,U,600);



fprintf(outputfile,"2: float array novalues ");



write_floatarray_elements(outputfile,V,600);



fprintf(outputfile,"3: float matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_floatmatrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}





/**-----------------------------------------------------------------------*/

/** Writes a matrix of double to a DEM turtle file */



/* Stampa un file  DEM in formato Fluidturtle

	!!! non viene usata la WORKING_DIRECTORY !!!

	Inputs:	matrix		matrice

			U,V 		vettori di header del DEM

			otputname 	nome del file di output

			comment		eventuali commenti da scrivere nel file */



void doublematrix_dem2(DOUBLEMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,

	char *outputname, char *comment,short print)



{

/* Declaration variables */

FILE *outputfile;



/* Writes output file */

if(print==1){

outputfile=t_fopen(outputname,"w");

fprintf(outputfile,"/**  %s */\n",comment);

fprintf(outputfile,"index{3}\n");

fprintf(outputfile,"1: double array pixels size  ");

write_doublearray_elements(outputfile,U,600);

fprintf(outputfile,"2: double array novalues ");

write_doublearray_elements(outputfile,V,600);

fprintf(outputfile,"3: double matrix data {%d,%d}\n",matrix->nrh,matrix->nch);

write_doublematrix_elements(outputfile,matrix,matrix->nch+1);

t_fclose(outputfile);

} else {

	return;

}

}













