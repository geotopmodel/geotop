/*
 * hpc.geotop.h
 *
 *  Created on: 02/mar/2012
 *      Author: Giuseppe Onorevoli
 */

#ifdef USE_HPC

#include <valarray>

struct GCSTRUCT	//Struct containing subdomain ghost-cells adjacency data for MPI SEND/RECV commands
{
	int rank;
	char calltype[4];
	int top;
	int left;
	int bottom;
	int right;
	struct GCSTRUCT *next;
};
typedef struct GCSTRUCT;

struct WORKAREA	//Struct containing local subdomain coords
{
	int rank;
	int top;
	int left;
	int bottom;
	int right;
	size_t start[];
	size_t count[];
	size_t stride[];
};
typedef struct WORKAREA;

#endif
