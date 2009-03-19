
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
	 *\note modified by Emanuele Cordano (March 2009)
	 */


	double delta,delta_new,alpha,beta,delta0;
	DOUBLEVECTOR *r, *d,*q,*y;

	int s,sl;
	long icnt_max;
	long j;
	double p;

	r=new_doublevector(x->nh);
	d=new_doublevector(x->nh);
	q=new_doublevector(x->nh);
	y=new_doublevector(x->nh);



	icnt=0;
	icnt_max=x->nh;



    s=(* funz)(y,x);
    delta_new=0.0;
    for (j=y->nl;j<=y->nh;j++) {
    	r->element[j]=b->element[j]-y->element[j];
    	d->element[j]=r->element[j];
    	delta_new+=r->element[j]*r->element[j];
    }
    delta0=1.0;
  //  printf("start delta_new epsilon :%le %le ",delta_new,epsilon);
 //   stop_execution();


	//while ((icnt<=icnt_max) && (delta_new>pow(epsilon,2.0)*delta0)) {
	while ((icnt<=icnt_max) && (max_doublevector(r)>epsilon)) {
		delta=delta_new;
		s=(* funz)(q,d);
		p=0.0;
		for(j=q->nl;j<=q->nh;j++) {
			p+=q->element[j]*d->element[j];
		/*	if (((j==y->nl) && sl==0 ) || (sl==1)) {
				    		printf("p =%le q[%ld]=%le d[%ld]=%le (j=%ld)  ",p,j,q->element[j],j,d->element[j],j);
				//    		if (delta_new==0.0) sl=1;
			}*/
		}
		alpha=delta_new/p;
		for(j=x->nl;j<=x->nh;j++) {
			x->element[j]=x->element[j]+alpha*d->element[j];
		}

	    s=(* funz)(y,x);
	    delta_new=0.0;
	    sl=0;
	    for (j=y->nl;j<=y->nh;j++) {
	   // 	r->element[j]=b->element[j]-y->element[j];
	    	r->element[j]=r->element[j]-alpha*q->element[j];
	    	delta_new+=r->element[j]*r->element[j];
/*    	if (((j==y->nl) && sl==0 ) || (sl==1)) {
	    		printf("delta_new =%le (j=%ld)  ",delta_new,j);
	//    		if (delta_new==0.0) sl=1;
	    	}*/
	    }
	    beta=delta_new/delta;
	   // double aa=1.0e-21;
	    printf("delta_new =%le p=%le alpha=%le beta=%le delta_max=%le\n",delta_new,p,alpha,beta,max_doublevector(r));
	    //stop_execution();

//		if (delta_new>delta) {
	//		printf("Problem does not converge delta=%lf delta_new=%lf beta=%lf icnt=%ld icnt_max=%ld epsilon=%lf \n ",delta,delta_new,beta,icnt,icnt_max,epsilon);
	//		stop_execution();
	//	}
	//	linear_comb_doublevector(d,r,d,1.0,beta);
		for (j=d->nl;j<=d->nh;j++) {
			 d->element[j]=r->element[j]+beta*d->element[j];
		}

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


long conjugate_gradient_search_LONG(long epsilon,LONGVECTOR *x,LONGVECTOR *b, int (* funz)(LONGVECTOR *y,LONGVECTOR *x)){

	/*
	 *\param icnt  - (long)
	 *\param epsilon - (double) required tollerance (2-order norm of the residuals)
	 *\param x     - (LONGVECTRO *) vector of the unknowns x in Ax=b
	 *\param b     - (LONGVECTOR *) vector of b in Ax=b
	 *\param (* funz)(LONGVECTOR *y,LONGVECTOR *x) - (int) pointer to the application A (x and y doublevector y=A(param)x ) it return 0 in case of success, -1 otherwise.
	 *
	 *\return the number of reitarations
	 *\brief algorithm proposed by Jonathan Richard Shewckuck in http://www.cs.cmu.edu/~jrs/jrspapers.html#cg . All the physical variables are expressed as long integer.
	 *
	 * \author Emanuele Cordano
	 * \date March 2009
	 *
	 *\note modified by Emanuele Cordano (March 2009)
	 */


	long delta,delta_new,alpha,beta,delta0;
	LONGVECTOR *r, *d,*q,*y;
	long icnt;
	int s,sl;
	long icnt_max;
	long j;
	long p;

	r=new_longvector(x->nh);
	d=new_longvector(x->nh);
	q=new_longvector(x->nh);
	y=new_longvector(x->nh);




	icnt_max=x->nh;
	icnt=0;


    s=(* funz)(y,x);
    delta_new=0.0;
    for (j=y->nl;j<=y->nh;j++) {
    	r->element[j]=b->element[j]-y->element[j];
    	d->element[j]=r->element[j];
    	delta_new+=r->element[j]*r->element[j];
    }
    printf("delta_new =%ld delta_max=%ld \n",delta_new,epsilon*epsilon);


//	while ((icnt<=icnt_max) && (delta_new>epsilon*epsilon))
	while ((icnt<=icnt_max) && (max_doublevector(r)<=epsilon)) {
		delta=delta_new;
		s=(* funz)(q,d);
		p=0.0;
		for(j=q->nl;j<=q->nh;j++) {
			p+=q->element[j]*d->element[j];

		}
		alpha=delta_new/p;
		for(j=x->nl;j<=x->nh;j++) {
			x->element[j]=x->element[j]+alpha*d->element[j];
		}

	    s=(* funz)(y,x);
	    delta_new=0.0;
	    sl=0;
	    for (j=y->nl;j<=y->nh;j++) {
	   // 	r->element[j]=b->element[j]-y->element[j];
	    	r->element[j]=r->element[j]-alpha*q->element[j];
	    	delta_new+=r->element[j]*r->element[j];
/*    	if (((j==y->nl) && sl==0 ) || (sl==1)) {
	    		printf("delta_new =%ld (j=%ld)  ",delta_new,j);
	//    		if (delta_new==0) sl=1;
	    	}*/
	    }
	    beta=delta_new/delta;
	   // double aa=1.0e-21;
	   printf("delta_new =%ld p=%ld alpha=%ld beta=%ld delta_max=%ld \n",delta_new,p,alpha,beta,epsilon*epsilon);
	    //stop_execution();

//		if (delta_new>delta) {
	//		printf("Problem does not converge delta=%lf delta_new=%lf beta=%lf icnt=%ld icnt_max=%ld epsilon=%lf \n ",delta,delta_new,beta,icnt,icnt_max,epsilon);
	//		stop_execution();
	//	}
	//	linear_comb_doublevector(d,r,d,1.0,beta);
		for (j=d->nl;j<=d->nh;j++) {
			 d->element[j]=r->element[j]+beta*d->element[j];
		}

		icnt++;


	}


	free_longvector(r);
	free_longvector(d);
	free_longvector(q);
	free_longvector(y);

	return icnt;

}


double max_doublevector(DOUBLEVECTOR *v) {
	/*
	 * \date March 2009
	 * \author Emanuele Cordano
	 *
	 * \brief maximum value in a doublevector
	 *
	 */
	 long j;
	 double MK=abs(v->element[v->nl]);

	 for (j=v->nl+1;j<=v->nh;j++){
		 MK=fmax(MK,abs(v->element[j]));
	 }

	 return MK;
}
