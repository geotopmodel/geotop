
/* MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file conjugate_gradient.c

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico and Riccardo Rigon

This file is part of MATH2.
 MATH2 is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MATH2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include "turtle.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "linear_span.h"

#include "conjugate_gradient.h"

#define MAX_REITERATON 1

#define DELTA_MIN epsilon

/* DA SCRIVERE!!!!
 *
 * int (* temporary_function)(DOUBLEVECTOR *y,DOUBLEVECTOR *x);
 * DOUBLEVECTOR *diagonal_elements;
 *
 * DOUBLEVECTOR *diagonal(int (* funz)(DOUBLEVECTOR *y,DOUBLEVECTOR *x)) { da fare }
 *
 * int conditioned_function(DOUBLEVECTOR *y, DOUBLEVECTOR *x) {
 *
 * }
 *
 * fare rutines gradiente coniugato precondizionato
 *
int (*condtioning (int (* funz)(DOUBLEVECTOR *y,DOUBLEVECTOR *x),DOUBLEVECTOR *b, DOUBLEVECTOR *b_cond))(DOUBLEVECTOR *x,DOUBLEVECTOR *y) {

int (* funz)(DOUBLEVECTOR *y,DOUBLEVECTOR *x);

}
*/

long conjugate_gradient_search(long icnt, double epsilon,  DOUBLEVECTOR *x, DOUBLEVECTOR *b, int (* funz)(DOUBLEVECTOR *y,DOUBLEVECTOR *x)){

	/*
	 *\param icnt  - (long)
	 *\param epsilon - (double) required tollerance (2-order norm of the residuals)
	 *\param x     - (DOUBLEVECTRO *) vector of the unknowns x in Ax=b
	 *\param b     - (DOUBLEVECTOR *) vector of b in Ax=b
	 *\param (* funz)(DOUBLEVECTOR *y,DOUBLEVECTOR *x) - (int) pointer to the application A (x and y doublevector y=A(param)x ) it return 0 in case of success, -1 otherwise.
	 *
	 *\return the number of reitarations
	 *\brief algorithm proposed by Jonathan Richard Shewckuck in http://www.cs.cmu.edu/~jrs/jrspapers.html#cg
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 */


	double delta,delta_new,alpha,beta;
	DOUBLEVECTOR *r, *d,*q,*y;
	//,*x_new, *r_new, *d_new;
	int s;
	long icnt_max;

	r=new_doublevector(x->nh);
	d=new_doublevector(x->nh);
	q=new_doublevector(x->nh);
	y=new_doublevector(x->nh);

//	r_new=new_doublevector(x->nh);
//	d_new=new_doublevector(x->nh);
//	x_new=new_doublevector(x->nh);

	icnt=0;
	icnt_max=x->nh;



    s=(* funz)(y,x);
	linear_comb_doublevector(r,b,y,1.0,-1.0);

	copy_doublevector(r,d);
	delta_new=prodscal(r,r);
	printf("delta_init= %lf\n",delta_new);

	delta=delta_new;

	while ((icnt<=icnt_max) && (delta_new>epsilon)) {


		s=(* funz)(q,d);
		alpha=delta_new/prodscal(d,q);
		linear_comb_doublevector(x,x,d,1.0,alpha);

		if (icnt%MAX_REITERATON==0) {
			s=(* funz)(y,x);
			linear_comb_doublevector(r,b,y,1.0,-1.0);
		}else{
			linear_comb_doublevector(r,r,q,1.0,-alpha);
//		    copy_doublevector(r_new,r);
		}
		delta=delta_new;
		delta_new=prodscal(r,r);
		printf ("delta_new=%le icnt=%ld icnt_max=%ld deltamin=%le\n",delta_new,icnt,icnt_max,epsilon);

		beta=delta_new/delta;
		if (delta_new>delta) {
			printf("Problem does not converge delta=%lf delta_new=%lf beta=%lf icnt=%ld icnt_max=%ld epsilon=%lf \n ",delta,delta_new,beta,icnt,icnt_max,epsilon);
	//		stop_execution();
		}
		linear_comb_doublevector(d,r,d,1.0,beta);

		icnt++;


	}


	free_doublevector(r);
	free_doublevector(d);
	free_doublevector(q);
	free_doublevector(y);



	return icnt;

}



int linear_comb_doublevector(DOUBLEVECTOR *result,DOUBLEVECTOR *a, DOUBLEVECTOR *b, double ca, double cb) {
	/*
	 * \param result - (DOUBLEVECTOR *) result vector (r=a*ca+b*cb)
	 * \param a      - (DOUBLECTOR *) vector
	 * \param b      - (DOUBLECTOR *) vector
	 * \param ca      - coefficient for "a" vector
	 * \param cb      - coefficient for "b" vector
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 *\return if r is solved correctly or -1 otherwise
	 *
	 */
	long i;

	if((result->nh!=a->nh) || (result->nh!=b->nh)  || (b->nh!=a->nh) ) {
		printf("Error: in linear_comb vectors result and a and b do not have the same number of elements!!  /n");
		return -1;
	}

	for (i=a->nl;i<=a->nh;i++){
		result->element[i]=ca*a->element[i]+cb*b->element[i];

	}

	return 0;

}


