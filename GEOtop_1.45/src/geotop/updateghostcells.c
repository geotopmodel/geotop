/*
 * updateghostcells.c
 *
 *  Created on: 02/mar/2012
 *      Author: giuseppe
 */

#ifdef USE_HPC

// include
#include "hpc.geotop.h"

#include <iostream>
#include <string>
#include <fstream>
#include <mpi.h>

// namespace
using namespace std;
using namespace MPI;


//************************************************************************************************************************************
void updateGhostcells(ALLDATA *all, GCSTRUCT *start)
{
	int sendtag, recvtag, sendcount, recvcount, startValue, stopValue, sendRank, recvRank;
	MPI_Status status;
	float sendbuf, recvbuf;
	sendtag = 0;
	recvtag = 0;

	GCSTRUCT *ptr = start->next;
	while(ptr)
	{
		//char *str = new char[ptr->calltype.size()];
		char *str = new char[5];
		//strcpy(str, ptr->calltype.c_str());
		strcpy(str, ptr->calltype);
		if (strcmp(str, "SEND") == 0) {
			//predispone il buffer delle ghost-cells da inviare
			if (ptr->top == ptr->bottom){	// buffer lungo riga
				startValue = min(ptr->left, ptr->right);
				stopValue = max(ptr->left, ptr->right);
				for (int i = startValue; i <= stopValue; i++){
					// caricamento del buffer di SEND
				}
			}
			else {	// buffer lungo colonna
				startValue = min(ptr->top, ptr->bottom);
				stopValue = max(ptr->top, ptr->bottom);
				for (int i = startValue; i <= stopValue; i++){
					// caricamento del buffer di SEND
				}
			}
			sendcount = stopValue - startValue;
			sendRank = ptr->rank;
		} else {
			//predispone il buffer delle ghost-cells da ricevere
			if (ptr->top == ptr->bottom){
				// buffer lungo riga
				startValue = min(ptr->left, ptr->right);
				stopValue = max(ptr->left, ptr->right);
				for (int i = startValue; i <= stopValue; i++){
					// caricamento del buffer
				}
			}
			else {
				// buffer lungo colonna
				startValue = min(ptr->top, ptr->bottom);
				stopValue = max(ptr->top, ptr->bottom);
				for (int i = startValue; i <= stopValue; i++){
					// caricamento del buffer
				}
			}
			recvcount = stopValue - startValue;
			recvRank = ptr->rank;
			// aggiorna le ghost-cells
			MPI_Sendrecv(&sendbuf, sendcount, MPI_DOUBLE, sendRank, sendtag, &recvbuf, recvcount, MPI_DOUBLE, recvRank, recvtag, MPI_COMM_WORLD, &status);
			std::cout << "Sending to: " << sendRank << " and receiving from: " << recvRank << std::endl;
			if (ptr->top == ptr->bottom){	// buffer lungo riga
				for (int i = startValue; i <= stopValue; i++){
					// scaricamento del buffer di RECV
				}
			}
			else {	// buffer lungo colonna
				for (int i = startValue; i <= stopValue; i++){
					// scaricamento del buffer di RECV
				}
			}
		}
		ptr = ptr->next;
	}
	return;
}

#endif
