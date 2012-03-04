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
					if (strcmp(str, "SIZE") == 0){
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
