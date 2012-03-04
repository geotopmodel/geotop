/*
 * hpc.geotop.h
 *
 *  Created on: 02/mar/2012
 *      Author: giuseppe
 */

#ifdef USE_HPC

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
};
typedef struct WORKAREA;

#endif
