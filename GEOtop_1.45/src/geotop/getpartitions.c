/*
 * getpartitions.c
 *
 *  Created on: 04/mar/2012
 *      Author: giuseppe
 */

#ifdef USE_HPC

#include "hpc.geotop.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

void getpartitions(GCSTRUCT *start, WORKAREA *rankArea)
{

	int count, rank, myrank, nprocs, cprocs;
	int Np, Nrt, Nct, Nr, Nc, offsetNr, offsetNc;
	string buffer,line;
	line = "none";
	char *str = new char[line.size()];

	ifstream dataFile("partition.txt", ios::in);

	if (!dataFile) {	//file not open
		std::cout << "Error opening data file.\n";
		return;
	}

	dataFile.exceptions(ios::badbit | ios::failbit | ios::eofbit);

	try {

		//scorre il file in cerca della riga corrispondente a myrank
		count = 0;
		while(getline(dataFile,line)){
			std::cout << "Stream pointer currently in line " << count++ << ".\n";
			std::cout << "Confirmation: " << line << "\n";
			char *str = new char[line.size()];
			//strcpy(str, line.c_str());
			//str = strtok (str, " ");
			str = strtok (str, " ");
			if (count == 1) cprocs = atoi(str);
			if (atoi(str) == myrank){
				while (str != NULL)
				{
					str = strtok (NULL, " ,.-)(");
					if (str == NULL) break;
					if (strcmp(str, "DATA") == 0){
						// Memorizza i parametri dell'area di lavoro
						Np = rank;
						Nrt = atoi(strtok (NULL, " ,.-)("));
						Nct = atoi(strtok (NULL, " ,.-)("));
					} else if (strcmp(str, "SIZE") == 0){
						// Memorizza i parametri dell'area di lavoro
						rankArea->rank = rank;
						rankArea->top = atoi(strtok (NULL, " ,.-)("));
						rankArea->left = atoi(strtok (NULL, " ,.-)("));
						rankArea->bottom = atoi(strtok (NULL, " ,.-)("));
						rankArea->right = atoi(strtok (NULL, " ,.-)("));
						//nr = abs(top - bottom);
						//nc = abs(left - right);
					} else if ((strcmp(str, "SEND") == 0) || (strcmp(str, "RECV")) == 0)  {
						GCSTRUCT* newStruct = new GCSTRUCT;
						newStruct->rank = rank;
						//newStruct->calltype = str;
						strcpy(newStruct->calltype, str);
						newStruct->top = atoi(strtok (NULL, " ,.-)("));
						newStruct->left = atoi(strtok (NULL, " ,.-)("));
						newStruct->bottom = atoi(strtok (NULL, " ,.-)("));
						newStruct->right = atoi(strtok (NULL, " ,.-)("));
						if(start == NULL) {
							start = newStruct;	//if the first node (first link) is null, set the memory there
							start->next = NULL;
							return;
						}

						GCSTRUCT *ptr = start;
						while(ptr->next != NULL) // loop fino all'ultimo elemento della lista
							ptr = ptr->next;

						// siamo arrivati all'ultimo elemento della lista
						ptr->next = newStruct;
						ptr->next->next = NULL;

					} else {
						// il rank del nuovo pacchetto di dati
						rank = atoi(str);
					}
				}
			}

			// definisce le dimensioni del sotto-dominio e completa i parametri necessari a dimensionare i sotto-domini in lettura e scrittura
			Nr = abs(rankArea->top - rankArea->bottom);
			Nc = abs(rankArea->left - rankArea->right);
			// write block offsets
			rankArea->countw[0]=rankArea->top;
			rankArea->countw[1]=rankArea->right;
			// write block dimensions
			rankArea->countw[0]=Nr;
			rankArea->countw[1]=Nc;
			if (rankArea->top == 1) {
				Nr = Nr + 1;
				offsetNr = rankArea->top;
			} else if (rankArea->bottom == Nrt) {
				Nr = Nr + 1;
				offsetNr = rankArea->top - 1;
			} else {
				offsetNr = rankArea->top - 1;
				Nr = Nr + 2;
			}
			if (rankArea->left == 1) {
				Nc = Nc + 1;
				offsetNc = rankArea->left;
			} else if (rankArea->right == Nct) {
				Nc = Nc + 1;
				offsetNc = rankArea->left - 1;
			} else {
				Nc = Nc + 2;
				offsetNc = rankArea->left;
			}
			// read blocks offsets
			rankArea->startr[0]=Nc;
			rankArea->startr[1]=Nr;
			// read block dimensions
			rankArea->startw[0]=offsetNc;
			rankArea->startw[1]=offsetNr;

			// sub-sampling steps (not used at the moment)
			rankArea->strider[0]=0;
			rankArea->strider[1]=0;
			// sub-sampling steps (not used at the moment)
			rankArea->stridew[0]=0;
			rankArea->stridew[1]=0;

		}

	} catch(ios_base::failure exc) {
		std::cout << "Error reading data file.\n";
	}
	try {
		dataFile.close();
	} catch (ios_base::failure exc) {
		std::cout << "Error closing data file.";
	}
	return;
}

#endif
