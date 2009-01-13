
/* MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file linear_span.c

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
#include "tensor3D.h"
#include "networks.h"
#include "t_utilities.h"
#include "linear_span.h"
//#include "gridded.element.input.geotop.h"
#include "rw_maps.h"


int no_value_function(double x,DOUBLEVECTOR *V) {
	/*!<      int no_value_function(double x,FLOATVECTOR* V)
	 * \param x a generic value of a map
	 * \param V the no-value doublevector of a T_INIT struct)
	 *
	 * \brief Function which recognizes novalue in a map (doublematrix)
	 * \datails it folloes FLUIDTURLE FORMALISM
	 * \return 0 if x is not a no-value,  1 if x is a no-value
	 *
	 * \author Emanuele Cordano
	 * \date 7 January 2007
	 *
	 * \version 0.9375b
	 *
	 * */
	if ((V->element[1]>0 && x>=V->element[2]) || (V->element[1]<0 && x<=V->element[2])) {
		/*!<  printf("1");*/
		return 1;
	}else{
		/*!<   printf("0");*/
		return 0;
	}

}

DOUBLEMATRIX *extract_a_new_map(DOUBLETENSOR *xtensor, long l,T_INIT *UVref) {
		/*!<      DOUBLEMATRIX *extract_a_new_map(DOUBLETENSOR *xtensor, long l,T_INIT *UVref)
		 *
		 * \author  Emanuele Cordano
		 * \date November, 2007
		 *
		 * \brief it allocates and extracts a new map from a doubletensor xtensor related to a specific layers;
		 * \param DOUBLETENSOR tensor releted to variable x
		 * \param l (integer) number of layer at which map is requessted
		 * \param UVref : T_INIT struct containing no-value information
		 *
		 * \return the double matrix map of x at the requested layer
		 *
		 * \details The output corresponds to tha map at the \f$l\f$-th layer of the input doubletensor.
		 *
		 * \warning l must be equal or less than the number  of layers of input doubletensor x
		 *
		 * \relates FLUIDTURTLE_LIBRARIES,  no_value_function
		 */
		 long r,c;
		 DOUBLEMATRIX *xm;
		 char *error_message;
		 xm=new_doublematrix(xtensor->nrh,xtensor->nch);
		 if (l<xtensor->ndl && l>xtensor->ndh) {

					sprintf(error_message,"Number of layers does not correspond in extract_new_map function with doubletensor NAME: %s  l=%ld ndl=%ld ndh=%ld \n",xtensor->name,l,xtensor->ndl,xtensor->ndh);
					t_error(error_message);

			}else{

				for(r=xm->nrl;r<=xm->nrh;r++){
					for(c=xm->ncl;c<=xm->nch;c++){
					    if (no_value_function(xtensor->element[l][r][c],UVref->V)==1){
					    	xm->element[r][c]=UVref->V->element[2];
					    }else{
					    	xm->element[r][c]=xtensor->element[l][r][c];
					    }
					}

				}
			return xm;
			}

		return NULL;
		}



DOUBLEVECTOR *prod_doublematvet(DOUBLEMATRIX *m, DOUBLEVECTOR *v)
{


	/*!
	 *    DOUBLEVECTOR *prod_doublematvet(DOUBLEMATRIX *m, DOUBLEVECTOR *v)
	 *
	 * \author Matteo Dall'Amico
	 *
	 *
	 *
	 * \details (in Italian) riceve in input una matrice e un vettore double e ne ritorna il vettore prodotto */

	DOUBLEVECTOR *p;
	p=new_doublevector (m->nch);
	int i,j;
	if (m->nch != v->nh) t_error("Error in prod_doublematvet(): The matrix and the vector have not proper dimensions\n");
	for(i=m->nrl;i<=m->nrh;i++) {
		double buf=0.0;
		for (j=m->ncl;j<=m->nch;j++) {
				buf+=m->element[i][j]*v->element[j];
				p->element[i]=buf;
		}
	}
	return p;
}


double prodscal (DOUBLEVECTOR *a, DOUBLEVECTOR *b)
/*!<
 *
 * \author Matteo Dall'Amico, Emanuele Cordano
 * \date 2006 / 2008
 *
 *
 * \param a - (DOUBLEVECTOR *)
 * \param b - (DOUBLEVECTOR *)
 *
 * \return scalar product between a and b
 *
 * */
{
	double p=0.0;
	int i;
	if (a->nh != b->nh)
	t_error("Error in prodscal(): The two vectors have not equal dimensions\n");
	for(i=a->nl;i<=a->nh;i++) {
		p+=a->element[i]*b->element[i];
		}
	return p;
}

DOUBLEVECTOR *scalxvet (double a, DOUBLEVECTOR *b)
/*!   Funzione che da' in output un vettore che e' il prodotto tra lo scalare a e il vettore b */
{
	DOUBLEVECTOR *p;
	p=new_doublevector(b->nh);
	int i;
	for(i=b->nl;i<=b->nh;i++) {
			p->element[i]= a * b->element[i];
	}
	return p;
}

	/*!
	 *    DOUBLETENSOR *linear_span_doubletensor(double c1, double c2, DOUBLETENSOR *T1,DOUBLETENSOR *T2)
	 * \author Emanuele Cordano
	 *
	 * \date 9 January 2008
	 * \param c1 first double coefficient
	 * \param c2 second double coefficient
	 * \param T1 first doubletensor
	 * \param T2 second doubletensor
	 * \param V doublecetor with no_balue information
	 *
	 *
	 *
	 * \return TL doubletensor
	 *
	 *
	 *
	 * \details it solves the linear span between T1 and T2 : $\f TL=c1 * T1 + c2 * T2 $\f . It takes into account no-value data.
	 *
	 */

DOUBLETENSOR *linear_span_doubletensor(double c1, double c2, DOUBLETENSOR *T1,DOUBLETENSOR *T2, DOUBLEVECTOR *V){

	long r,c,l;
	DOUBLETENSOR* TL;

	if (T1->nrh!=T2->nrh) {
		 t_error("Error in function linear_span_doubletensor: doubltensors do  not have the same numbers of rows!!  ");
	}else if (T1->nch!=T2->nch){
		t_error("Error in function linear_span_doubletensor: doubletensors do not have the same numbers of columns!! ");
	}else if(T1->ndh!=T2->ndh) {
		t_error("linear_span_doubletensor: doubletensors do not have the same numbers of columns!!");
	} else{
		TL=new_doubletensor(T1->ndh,T1->nrh,T1->nch);
		for(r=T1->nrl;r<=T1->nrh;r++){
			for(c=T1->ncl;c<=T1->nch;c++){
				for(l=T1->ndl;l<=T1->ndh;l++){
					if (no_value_function(T1->element[l][r][c],V)==1 && no_value_function(T2->element[l][r][c],V)==1){
						TL->element[l][r][c]=V->element[2];
					}else {
						TL->element[l][r][c]=V->element[2]=c1*T1->element[l][r][c]+c2*T2->element[l][r][c];
					}
					}
				}

			}
		    return TL;
	}
return NULL;

}


DOUBLEMATRIX *linear_span_doublematrix(double c1, double c2, DOUBLEMATRIX *M1a,DOUBLEMATRIX *M2a, DOUBLEVECTOR *V){
	/*!
	 *    DOUBLETENSOR *linear_span_doublematrix(double c1, double c2, DOUBLEMATRIX *T1,DOUBLEMATRIX *T2)
	 * \author Emanuele Cordano
	 *
	 * \date 9 January 2008
	 * \param c1 first double coefficient
	 * \param c2 second double coefficient
	 * \param M1a first doublematrix
	 * \param M2a second doublematrix
	 * \param V floatvector with no_balue information
	 * \brief
	 *
	 *
	 * \return ML doublematrix
	 *
	 * \relates no_value_function
	 *
	 * \details it solves the linear span between M1 and M2 : $\f ML=c1 * M1 + c2 * M2 $\f . It takes into account no-value data.
	 *
	 */


	long r,c;
	DOUBLEMATRIX *ML;

	if (M1a->nrh!=M2a->nrh) {
		 t_error("Error in function linear_span_doublematrix: doublematrices do  not have the same numbers of rows!!  ");
	}else if (M1a->nch!=M2a->nch) t_error("Error in function linear_span_doublematrix: doublematrices do not have the same numbers of columns!! ");

	ML=new_doublematrix(M1a->nrh,M1a->nch);
	for(r=M1a->nrl;r<=M1a->nrh;r++){
		for(c=M1a->ncl;c<=M1a->nch;c++){
				if (no_value_function(M1a->element[r][c],V)==1 && no_value_function(M2a->element[r][c],V)==1){
					ML->element[r][c]=V->element[2];
		        }else {
					ML->element[r][c]=V->element[2]=c1*M1a->element[r][c]+c2*M2a->element[r][c];

				}
		}

	}
	return ML;
	}


DOUBLEMATRIX *transpose_doublematrix(DOUBLEMATRIX *M){
	/*!
	 *    DOUBLEMATRIX *transpose_doublematrix(DOUBLEMATRIX *M)
	 *
	 *    \param (DOULEMATRIX *) - M a matrix of double to be transposed
	 *
	 *
	 *    \return  It calculates the transpose matrix of M
	 *
	 *    \author Emanuele Cordano
	 *    \date August 2008
	 *
	 *    */
	DOUBLEMATRIX *MR;
	long r,c;

	MR=new_doublematrix(M->nch,M->nrh);
	for(r=MR->nrl;r<=MR->nrh;r++){
		for(r=MR->nrl;r<=MR->nrh;r++){
			MR->element[r][c]=M->element[c][r];
		}

	}

	return MR;

}

DOUBLEVECTOR *extract_a_column_from_doublematrix(long d,DOUBLEMATRIX *M) {
	/*!
		 *    DOUBLEVECTOR *extract_a_coluumn_from_doublematrix(long l,DOUBLEMATRIX *M)
		 *
		 *    \param (DOULEMATRIX *) - M a matrix of double from which the column will be estracted
		 *    \param (long)          - d clunm at
		 *
		 *
		 *    \return  DOUBLEVECTOR - a doublevector extracted and corresponding to the d-th column
		 *    \author Emanuele Cordano
		 *    \date August 2008
		 *

		 *    */

		DOUBLEVECTOR *VC;
		long r,c;

		VC=new_doublevector(M->nrh);
		if ((d>=M->ncl) && (d<=M->nch)) {
			for(r=M->nrl;r<=M->nrh;r++){
						VC->element[r]=M->element[r][d];
			}
		} else {
			printf ("Warning: d=%ld is out of range: %ld to %ld columns",d,M->ncl,M->nch);
		}

		return VC;


}


DOUBLEVECTOR *extract_a_row_from_doublematrix(long d,DOUBLEMATRIX *M) {
	/*!
		 *    DOUBLEVECTOR *extract_a_row_from_doublematrix(long l,DOUBLEMATRIX *M)
		 *
		 *    \param (DOULEMATRIX *) - M a matrix of double from which the row will be extracted
		 *    \param (long)          - d row
		 *
		 *
		 *    \return  DOUBLEVECTOR - a doublevector extracted and corresponding to the d-th row
		 *    \author Emanuele Cordano
		 *    \date August 2008
		 *

		 *    */

		DOUBLEVECTOR *VC;
		long c;

		VC=new_doublevector(M->nch);
		if ((d>=M->nrl) && (d<=M->nrh)) {
			for(c=M->nrl;c<=M->nrh;c++){
						VC->element[c]=M->element[d][c];
			}
		} else {
			printf ("Warning: d=%ld is out of range: %ld to %ld rows",d,M->nrl,M->nrh);
		}

		return VC;


}




DOUBLEVECTOR *extract_a_vertical_column_from_doubletensor(long r,long c, DOUBLETENSOR *T) {
	/*!
		 *    DOUBLEVECTOR *extract_a_vertical_column_from_doubletensor (long r,long c, DOUBLETENSOR *T)
		 *
		 *    \param (DOULEMATRIX *) - T a doubletensor of double from which the vertical column will be estracted
		 *    \param (long)          - r  row index of the vertical column;
		 *	  \param (long)          - c  column index of the verical column
		 *
		 *   \return  DOUBLEVECTOR - a doublevector extracted and corresponding to the list of all layers at the row r and the column c of the doubletensor T
		 *
		 *    \author Emanuele Cordano
		 *    \date August 2008
		 *
		 *    */

		DOUBLEVECTOR *VC;
		long i;

		VC=new_doublevector(T->ndh);
		if ((c>=T->ncl) && (c<=T->nch) && (r>=T->nrl) && (r<=T->nrh) ) {
			for(i=VC->nl;i<=VC->nh;i++){
				VC->element[i]=T->element[i][r][c];
			}
		} else {
			printf ("Warning: r=%ld  c=%ld is out of range: %ld to %ld columns and %ld to %ld rows ",r,c,T->nrl,T->nrh,T->ncl,T->nch);
		}

		return VC;



}





