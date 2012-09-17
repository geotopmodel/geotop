#ifndef WRITE_DEM_H
#define WRITE_DEM_H

#include "turtle.h"
#include "t_utilities.h"
#include "t_datamanipulation.h"

/**
Name: elements_dem
Version: 1.0
Synopsis:	void shortmatrix_dem(SHORTMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,
			char *outputname, char *comment, short print);
Description: write an array of the specified type on a specified file in the turtle DEM format
Authors & Date: Giacomo Bertoldi, September 2000
FILE: LIBRARIES/GEOMORPHOLOGYLIB/write_dem.h, LIBRARIES/GEOMORPHOLOGYLIB/write_dem.c
Inputs:	matrix		pointer to the matrix or vector with the data

		U			pointer to the vector with the header (usually [dx, dy, x0, y0])

		V 			pointer to the vector with novalue (usually [X, NV], where X is -1 if NV is less than the mininmum value of the data,

					0 if NV is the 0, 1 if NV is greater than the maximum value of the data)

		outputname 	output files name (if is declared a string WORKING_DIRECTORY as extern char, the path not must be specificated)

		comment		a string with comments to write into the file

		print		if print=PRINT the command is executed, if print=NOPRINT the command is skipped
Return: void
See Also: write__
Keywords: turtle file
Examples: APPLICATIONS/HYDROLOGY/geotop/geotop.c
*/

void shortmatrix_dem(SHORTMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,
	char *outputname, char *comment,short print);
	

void shortmatrix_dem3(SHORTMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,
	char *outputname, char *comment,short print);


void longmatrix_dem(LONGMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,
	char *outputname, char *comment,short print);



void intmatrix_dem(INTMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,
	char *outputname, char *comment,short print);



void floatmatrix_dem(FLOATMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,
	char *outputname, char *comment,short print);



void doublematrix_dem(DOUBLEMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,
	char *outputname, char *comment,short print);
	

void doublematrix_dem3(DOUBLEMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,
	char *outputname, char *comment,short print);



void doublematrix_control(DOUBLEMATRIX *matrix,
	char *outputname, char *comment,short print);



void floatvector_dem(FLOATVECTOR *vector,
	char *outputname, char *comment,short print);



void doublevector_dem(DOUBLEVECTOR *vector,
	char *outputname, char *comment,short print);



void doubletensor_dem(DOUBLETENSOR *tensor,long layer,DOUBLEVECTOR *U,
                      DOUBLEVECTOR *V,char *outputname,char *comment,short print);



void shortmatrix_dem2(SHORTMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,
	char *outputname, char *comment,short print);



void longmatrix_dem2(LONGMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,
	char *outputname, char *comment,short print);



void intmatrix_dem2(INTMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,
	char *outputname, char *comment,short print);



void floatmatrix_dem2(FLOATMATRIX *matrix, FLOATVECTOR *U, FLOATVECTOR *V,

	char *outputname, char *comment,short print);

void doublematrix_dem2(DOUBLEMATRIX *matrix, DOUBLEVECTOR *U, DOUBLEVECTOR *V,
	char *outputname, char *comment,short print);


#endif








