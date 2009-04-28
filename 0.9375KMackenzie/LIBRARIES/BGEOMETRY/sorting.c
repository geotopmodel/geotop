/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
BGEOMETRY Version 0.9375 KMackenzie

file sorting.c

Copyright, 2009 Emanuele Cordano and Riccardo Rigon

This file is part of BGEOMETRY.
 BGEOMETRY is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BGEOMETRY is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "turtle.h"
#include "sorting.h"

/*
int swap(long a, long b) {

	 *
	 * \param long a,b two (long) integer numbers
	 *
	 * \date March 2009
	 * \author Emanuele Cordano;
	 *
	 * \return 0 if a>b and swaps a and b, 1 if a<=b

	long mk;

	if (a>b) {
		mk=a;
		a=b;
		b=mk;
		return 0;
	} else {
		return 1;
	}

	return 0;

}

*/
int bubble_sort(LONGVECTOR *v,short print) {
	/*
	 *\param v (LONGVECTOR *) - vector to be sorted
	 *
	 *\date March 2008
	 *\author Emanuele Cordano
	 *
	 *\brief it makes a bubble sort of the vector v (see the pseudo-code on http://en.wikipedia.org/wiki/Bubble_sort#Pseudocode_implementation )
	 */
	long i,mk;

	int swapped=1;

	long n=v->nh;
	do {
		swapped=0;
		n=n-1;
		for(i=v->nl;i<=n;i++){
			if (v->element[i]>v->element[i+1]){
				/* swap  v->element[i] and v->element[i+1]*/
				mk=v->element[i];
				v->element[i]=v->element[i+1];
				v->element[i+1]=mk;
				swapped=1;
			}
		}
		if (n<v->nl && swapped==1) {
			printf("Error in bubble_sort longvector was not correctly sorted!! ");
			print_longvector_elements(v,print);
			swapped=1;
		}
	} while (swapped==1);

	return 0;
}

LONGVECTOR *addresses_bubble_sort(LONGVECTOR *v,short print) {
	/*
	 *\param v (LONGVECTOR *) - vector to be sorted
	 *
	 *\date March 2008
	 *\author Emanuele Cordano
	 *
	 *\brief it makes a bubble sort of the vector v (see the pseudo-code on http://en.wikipedia.org/wiki/Bubble_sort#Pseudocode_implementation )
	 *
	 *\return the address of the "old" array to be sorted.
	 *
	 *
	 *
	 */

	long i,mk,lk;

	int swapped=1;

	long n=v->nh;

	LONGVECTOR *la;

	la=new_longvector(v->nh);

	for (i=la->nl;i<=la->nh;i++) {
		la->element[i]=i;
	}


	do {
		swapped=0;
		n=n-1;
		for(i=v->nl;i<=n;i++){
			if (v->element[i]>v->element[i+1]){
				/* swap  v->element[i] and v->element[i+1]*/
				mk=v->element[i];
				v->element[i]=v->element[i+1];
				v->element[i+1]=mk;
				/* swap  v->element[i] and v->element[i+1]*/
				lk=la->element[i];
				la->element[i]=la->element[i+1];
				la->element[i+1]=lk;
				swapped=1;
			}
		}
		if (n<v->nl && swapped==1) {
			printf("Error in bubble_sort longvector was not correctly sorted!! ");
			print_longvector_elements(v,print);
			swapped=1;
		}
	} while (swapped==1);

	return la;
}


int bubble_sort_matrix(LONGVECTOR *v, LONGMATRIX *m, short print) {
	/*
	 * \author Emanuele Cordano
	 * \date March 2009
	 *
	 * \param v (LONGVECTOR *) - vector to be sorted
	 * \param M (LONGMATRIX *) - matrix
	 * \param print (short)
	 *
	 *
	 * \brief It sorts v and then replace
	 */

	LONGVECTOR *iv;
	long i,r,c,mk,lk;

	if (v->nh!=m->nrh) printf("Error in bubble_sort_matrix the vector to b sorted and the matrix has different size: %ld and %ld (rows)/n",v->nh,m->nrh);


	iv=addresses_bubble_sort(v,print);
	for (r=m->nrl;r<=m->nrh;r++) {
		lk=iv->element[r];
		if (lk!=r) {
			for (c=m->ncl;c<=m->nch;c++){
				mk=m->element[r][c];
				m->element[r][c]=m->element[lk][c];
				m->element[lk][c]=mk;
			}
			iv->element[r]=r;
			iv->element[lk]=lk;
		}
	}




	if (print==1) printf("Function bubble_sort_matrix was successfully executed! \n");

	return 0;

}


//int sort_by_address(LONGVECTOR *all){
	/*
	 *
	 *\author Emanuele Cordano
	 * \date March 2009
	 *
	 *
	 */

//	LONGVECTOR *addresses_bubble_sort(LONGVECTOR *v,short print)
//}//
