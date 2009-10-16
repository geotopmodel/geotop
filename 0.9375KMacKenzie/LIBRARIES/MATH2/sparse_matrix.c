
/*! MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file sparse_matrix.c

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
/*!
 * \file sparse_matrix.c
 *
 * \author Emanuele Cordano
 *
 *
 */
#include "turtle.h"
#include "sparse_matrix.h"
#include "t_utilities.h"
#include "util_math.h"

#define MAX_VALUE_DIAG 1e-8
#define MAX_ITERATIONS 50

int get_diagonal(DOUBLEVECTOR *diagonal, t_Matrix_element_with_voidp Matrix, void *data) {
	/*
	 *
	 * \author Emanuele Cordano
	 * \date May 2008
	 *
	 *\param diagonal (DOUBLEVECTOR *) - the diagonal doublevector of the matrix
	 *\param Matrix (t_Matrix) - a matrix from which the diagonal is extracted
	 *
	 *\brief it saved the square root of diagonal of Matrix in a doublevector
	 *
	 *\return 0 in case of success , -1 otherwise
	 *
	 */

	long i;
	DOUBLEVECTOR *x_v;

	x_v=new_doublevector(diagonal->nh);
	for (i=x_v->nl;i<=x_v->nh;i++) {
		x_v->co[i]=0.0;

	}
	for (i=x_v->nl;i<=x_v->nh;i++) {
			x_v->co[i]=1.0;
			diagonal->co[i]=(*Matrix)(i,x_v,data);
			//diagonal->co[i]=1.;
			x_v->co[i]=0.0;
	}

	free_doublevector(x_v);


	return 0;
}

int get_upper_diagonal(DOUBLEVECTOR *udiagonal, t_Matrix_element_with_voidp Matrix, void *data) {
	/*
	 *
	 * \author Emanuele Cordano,Stefano Endrizzi
	 * \date May August 2009
	 *
	 *\param diagonal (DOUBLEVECTOR *) - the upper diagonal doublevector of the matrix
	 *\param Matrix (t_Matrix) - a matrix from which the diagonal is extracted
	 *
	 *\brief it saved  the upper diagonal of Matrix in a doublevector
	 *
	 *\return 0 in case of success , -1 otherwise
	 *
	 */

	long i;
	DOUBLEVECTOR *x_v;

	x_v=new_doublevector(udiagonal->nh+1);
	for (i=x_v->nl;i<=x_v->nh;i++) {
		x_v->co[i]=0.0;

	}
	for (i=udiagonal->nl;i<=udiagonal->nh;i++) {
			x_v->co[i]=1.0;
			udiagonal->co[i]=(*Matrix)(i+1,x_v,data);
			//diagonal->co[i]=1.;
			x_v->co[i]=0.0;
	}

	free_doublevector(x_v);


	return 0;
}

long tridiag_preconditioned_conjugate_gradient_search(double epsilon,  DOUBLEVECTOR *x, DOUBLEVECTOR *b,
	t_Matrix_element_with_voidp function, void *data){

	/*!
	 *\param icnt  - (long) number of reiterations
	 *\param epsilon - (double) required tollerance (2-order norm of the residuals)
	 *\param x     - (DOUBLEVECTOR *) vector of the unknowns x in Ax=b
	 *\param b     - (DOUBLEVECTOR *) vector of b in Ax=b
	 *\param funz  - (t_Matrix_element_with_voidp) - (int) pointer to the application A (x and y doublevector y=A(param)x ) it return 0 in case of success, -1 otherwise.
	 *\param data - (void *) data and parameters related to the  argurment  t_Matrix_element_with_voidp funz
	 *
	 *
	 *\brief algorithm proposed by Jonathan Richard Shewckuck in http://www.cs.cmu.edu/~jrs/jrspapers.html#cg and http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
	 *
	 * \author Emanuele Cordano
	 * \date June 2009
	 *
	 *\return the number of reitarations
	 */


	double delta,delta_new,alpha,beta;
	DOUBLEVECTOR *r, *d,*q,*y,*sr,*diag,*udiag;

	int sl;
	long icnt_max;
	long icnt;
	long j;
	double p;

	r=new_doublevector(x->nh);
	d=new_doublevector(x->nh);
	q=new_doublevector(x->nh);
	y=new_doublevector(x->nh);
	sr=new_doublevector(x->nh);
	diag=new_doublevector(x->nh);
	udiag=new_doublevector(x->nh-1);

	icnt=0;
	icnt_max=x->nh;
	//icnt_max=100.0;

	for (j=x->nl;j<=x->nh;j++){
		y->co[j]=(*function)(j,x,data);
		//printf("j:%ld %f\n",j,y->co[j]);
	}
	//stop_execution();


    get_diagonal(diag,function,data);
    get_upper_diagonal(udiag,function,data);

//    print_doublevector_elements(diag,PRINT);
//    stop_execution();

    delta_new=0.0;

    for (j=y->nl;j<=y->nh;j++) {

    	r->co[j]=b->co[j]-y->co[j];
		//printf("j:%ld r:%le\n",j,r->co[j]);
    	if (diag->co[j]<0.0) {
    		diag->co[j]=1.0;
    		printf("\n Error in jacobi_preconditioned_conjugate_gradient_search function: diagonal of the matrix (%lf) is negative at %ld \n",diag->co[j],j);
    		stop_execution();
    	}
    	diag->co[j]=fmax(diag->co[j],MAX_VALUE_DIAG);

    }
    tridiag(0,0,0,x->nh,udiag,diag,udiag,r,d);

    for (j=y->nl;j<=y->nh;j++) {
    	//d->co[j]=r->co[j]/diag->co[j];
    	delta_new+=r->co[j]*d->co[j];
    }

	//printf("%f\n",max_doublevector(r));
	//stop_execution();

	while ((icnt<=icnt_max) && (max_doublevector(r)>epsilon)) {

		delta=delta_new;
	//	s=(* funz)(q,d);
		p=0.0;

		for(j=q->nl;j<=q->nh;j++) {
			q->co[j]=(*function)(j,d,data);
			p+=q->co[j]*d->co[j];

		}
		alpha=delta_new/p;
		for(j=x->nl;j<=x->nh;j++) {
			x->co[j]=x->co[j]+alpha*d->co[j];
		}


	    delta_new=0.0;
	    sl=0;
	    for (j=y->nl;j<=y->nh;j++) {
	    	if (icnt%MAX_ITERATIONS==0) {
					y->co[j]=(*function)(j,x,data);
					r->co[j]=b->co[j]-y->co[j];
	    	} else {
					r->co[j]=r->co[j]-alpha*q->co[j];
	    	}
	    }

	    tridiag(0,0,0,x->nh,udiag,diag,udiag,r,sr);

	    for (j=y->nl;j<=y->nh;j++) {
	    	delta_new+=sr->co[j]*r->co[j];

	    }
	    beta=delta_new/delta;
	   // double aa=1.0e-21;
	 // printf("alpha=%le beta=%le delta_max=%le\n",alpha,beta,max_doublevector(r));
	  //printf(" iter:%ld/%ld \n ",icnt,icnt_max);
	   //stop_execution();
		for (j=d->nl;j<=d->nh;j++) {
			 d->co[j]=sr->co[j]+beta*d->co[j];
		}

		icnt++;


	}

	printf("alpha=%le beta=%le delta_max=%le\n",alpha,beta,max_doublevector(r));
	printf("iter:%ld/%ld\n",icnt,icnt_max);

	free_doublevector(udiag);
	free_doublevector(diag);
	free_doublevector(sr);
	free_doublevector(r);
	free_doublevector(d);
	free_doublevector(q);
	free_doublevector(y);




	return icnt;

}


double max_doublevector(DOUBLEVECTOR *v) {
	/*!
	 * \date March 2009
	 * \author Emanuele Cordano
	 *
	 * \brief maximum value in a doublevector
	 *
	 */
	 long j;
	 double MK=fabs(v->co[v->nl]);

	 //printf("MK:%f\n",MK);

		//stop_execution();

	  for (j=v->nl+1;j<=v->nh;j++){
		 MK=fmax(MK,fabs(v->co[j]));
		// printf("j:%ld MK:%f abs:%f\n",j,MK,fabs(v->co[j]));
	 }

	 return MK;
}


/*
 * for r:= 1 step 1 until n-1 do
	d := 1/ arr
		for i := (r+1) step 1 until n do
			if (i,r)S then
				e := dai,r;
				ai,r := e ;
				for j := (r+1) step 1 until n do
					if ( (i,j)S ) and ( (r,j)S ) then
						ai,j := ai,j - e ar,j
					end if
				end (j-loop)
			end if
	end (i-loop)
end (r-loop)
 */

/*LONGMATRIX *extract_longmatrix(t_Matrix_element_with_voidp function, void *data,long i_nodes) {
	long i;
	LONGVECTOR *e;

	e=new_longvector(i_nodes);
	for (i=e->nl;i<=e->nh;i++) {
		e->element[i]=0.0;
	}

}*/


