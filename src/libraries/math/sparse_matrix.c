
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 1.225-15
 
 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "turtle.h"
#include "sparse_matrix.h"
#include "t_utilities.h"
#include "util_math.h"

#define MAX_VALUE_DIAG 1e-8
#define MAX_ITERATIONS 50

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int get_diagonal(DOUBLEVECTOR *diagonal, DOUBLEVECTOR *x0, double dt, t_Matrix_element_with_voidp Matrix, void *data) {

	long i;
	DOUBLEVECTOR *x_v;

	x_v=new_doublevector(diagonal->nh);
	for (i=x_v->nl;i<=x_v->nh;i++) {
		x_v->co[i]=0.0;
	}
	for (i=x_v->nl;i<=x_v->nh;i++) {
		x_v->co[i]=1.0;
		diagonal->co[i]=Fmax( (*Matrix)(i,x_v,x0,dt,data) , MAX_VALUE_DIAG );
		x_v->co[i]=0.0;
	}

	free_doublevector(x_v);


	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

int get_upper_diagonal(DOUBLEVECTOR *udiagonal, DOUBLEVECTOR *x0, double dt, t_Matrix_element_with_voidp Matrix, void *data) {
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
			udiagonal->co[i]=(*Matrix)(i+1,x_v,x0,dt,data);
			//diagonal->co[i]=1.;
			x_v->co[i]=0.0;
	}

	free_doublevector(x_v);


	return 0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long CG(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *x0, double dt, 
		DOUBLEVECTOR *b, t_Matrix_element_with_voidp function, void *data){

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

	short sux;
	double delta,alpha,beta,delta_new;
	DOUBLEVECTOR *r, *d,*q,*y,*sr,*diag,*udiag;

	long icnt_max;
	long icnt;
	long j;
	double p;
	double norm_r0;

	r=new_doublevector(x->nh);
	d=new_doublevector(x->nh);
	q=new_doublevector(x->nh);
	y=new_doublevector(x->nh);
	sr=new_doublevector(x->nh);
	diag=new_doublevector(x->nh);
	udiag=new_doublevector(x->nh-1);
	
	icnt=0;
	icnt_max=x->nh;
	//icnt_max=(long)(sqrt((double)x->nh));

	for (j=x->nl;j<=x->nh;j++){
		y->co[j]=(*function)(j,x,x0,dt,data);
	}

    get_diagonal(diag,x0,dt,function,data);
    get_upper_diagonal(udiag,x0,dt,function,data);
	
    delta_new=0.0;

    for (j=y->nl;j<=y->nh;j++) {

    	r->co[j]=b->co[j]-y->co[j];

    	if (diag->co[j]<0.0) {
    		diag->co[j]=1.0;
    		printf("\n Error in jacobi_preconditioned_conjugate_gradient_search function: diagonal of the matrix (%lf) is negative at %ld \n",diag->co[j],j);
    		t_error("Fatal Error");
    	}
		
    }
	
    sux=tridiag(0,0,0,x->nh,udiag,diag,udiag,r,d);
	if(sux==0) return(-1);

    for (j=y->nl;j<=y->nh;j++) {
    	//d->co[j]=r->co[j]/diag->co[j];
    	delta_new+=r->co[j]*d->co[j];
    }
	
	norm_r0 = norm_2(r, r->nl, r->nh);
					
	while ( icnt<=icnt_max && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {

		delta=delta_new;
		p=0.0;

		for(j=q->nl;j<=q->nh;j++) {
			q->co[j]=(*function)(j,d,x0,dt,data);
			p+=q->co[j]*d->co[j];

		}
		alpha=delta_new/p;
		for(j=x->nl;j<=x->nh;j++) {
			x->co[j]=x->co[j]+alpha*d->co[j];
		}


	    delta_new=0.0;
	    for (j=y->nl;j<=y->nh;j++) {
	    	if (icnt%MAX_ITERATIONS==0) {
					y->co[j]=(*function)(j,x,x0,dt,data);
					r->co[j]=b->co[j]-y->co[j];
	    	} else {
					r->co[j]=r->co[j]-alpha*q->co[j];
	    	}
	    }

	    sux=tridiag(0,0,0,x->nh,udiag,diag,udiag,r,sr);
		if(sux==0) return(-1);

	    for (j=y->nl;j<=y->nh;j++) {
	    	delta_new+=sr->co[j]*r->co[j];

	    }
	    beta=delta_new/delta;
		for (j=d->nl;j<=d->nh;j++) {
			 d->co[j]=sr->co[j]+beta*d->co[j];
		}

		icnt++;

	}

	free_doublevector(udiag);
	free_doublevector(diag);
	free_doublevector(sr);
	free_doublevector(r);
	free_doublevector(d);
	free_doublevector(q);
	free_doublevector(y);


	return icnt;

}

//******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long BiCGSTAB(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *x0, 
			  double dt, DOUBLEVECTOR *b, t_Matrix_element_with_voidp function, void *data){

	/* \author Stefano Endrizzi
	   \date March 2011	
	   from BI-CGSTAB: A FAST AND SMOOTHLY CONVERGING VARIANT OF BI-CG FOR THE SOLUTION OF NONSYMMETRIC LINEAR SYSTEMS
	   by H. A. VAN DER VORST - SIAM J. ScI. STAT. COMPUT. Vol. 13, No. 2, pp. 631-644, March 1992
	*/

	DOUBLEVECTOR *r0, *r, *p, *v, *s, *t, *diag, *udiag, *y, *z;
	//DOUBLEVECTOR *tt, *ss;
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j;
	short sux;
	
	r0 = new_doublevector(x->nh);
	r = new_doublevector(x->nh);
	p = new_doublevector(x->nh);
	v = new_doublevector(x->nh);
	s = new_doublevector(x->nh);
	t = new_doublevector(x->nh);
	diag = new_doublevector(x->nh);
	udiag = new_doublevector(x->nh-1);
	y = new_doublevector(x->nh);
	z = new_doublevector(x->nh);
	//tt = new_doublevector(x->nh);
	//ss = new_doublevector(x->nh);
	
	get_diagonal(diag,x0,dt,function,data);
    get_upper_diagonal(udiag,x0,dt,function,data);

    for (j=x->nl;j<=x->nh;j++ ) {		
    	r0->co[j] = b->co[j] - (*function)(j,x,x0,dt,data);
		r->co[j] = r0->co[j];
		p->co[j] = 0.;
		v->co[j] = 0.;
    }
	
	norm_r0 = norm_2(r0, r0->nl, r0->nh);
	
	rho = 1.;
	alpha = 1.;
	omeg = 1.;
	
	while ( i<=x->nh && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);
		
		beta = (rho1/rho)*(alpha/omeg);
		
		rho = rho1;
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			p->co[j] = r->co[j] + beta*(p->co[j] - omeg*v->co[j]);
		}
		
		sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, p, y);
		if(sux==0) return(-1);
		
		for (j=x->nl;j<=x->nh;j++ ) {	
			v->co[j] = (*function)(j,y,x0,dt,data);
		}
		
		alpha = rho/product(r0, v);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			s->co[j] = r->co[j] - alpha*v->co[j];
		}
		
		sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, s, z);	
		if(sux==0) return(-1);

		for (j=x->nl;j<=x->nh;j++ ) {	
			t->co[j] = (*function)(j,z,x0,dt,data);
		}
		
		/*tridiag(0, 0, 0, x->nh, udiag, diag, udiag, t, tt);
		tridiag(0, 0, 0, x->nh, udiag, diag, udiag, s, ss);
		omeg = product(tt, ss)/product(tt, tt);*/
		
		omeg = product(t, s)/product(t, t);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			x->co[j] += (alpha*y->co[j] + omeg*z->co[j]);
			r->co[j] = s->co[j] - omeg*t->co[j];
		}
		
		i++;
				
	}
	
	
	free_doublevector(r0);
	free_doublevector(r);
	free_doublevector(p);
	free_doublevector(v);
	free_doublevector(s);
	free_doublevector(t);
	free_doublevector(diag);
	free_doublevector(udiag);
	free_doublevector(y);
	free_doublevector(z);
	//free_doublevector(tt);
	//free_doublevector(ss);
	
	return i;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/*ALL THE FOLLOWING SUBROUTINES USE THE STORING METHOD OF Ai, Ap, Ax instead of using the function Jd giving the product of the Jacobian
by the Newton direction */

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void product_using_only_lower_diagonal_part(DOUBLEVECTOR *product, DOUBLEVECTOR *x, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){

	long c, r, i;
	
	for(i=1;i<=x->nh;i++){
		product->co[i] = 0.0;
	}
	
	c = 1;
	for(i=1;i<=Ai->nh;i++){
		r = Ai->co[i];
		
		product->co[c] += x->co[r] * Ax->co[i];
		if(r > c){
			product->co[r] += x->co[c] * Ax->co[i];
		}else if(r < c){
			printf("r:%ld c:%ld i:%ld Ap[c]:%ld tot:%ld %ld %ld\n",r,c,i,Ap->co[c],Ai->nh,Ap->co[x->nh],x->nh);
			t_error("matrix is not L");
		}
		
		if(i == Ap->co[c]) c ++;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void product_using_only_upper_diagonal_part(DOUBLEVECTOR *product, DOUBLEVECTOR *x, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){
	
	long c, r, i;
	
	for(i=1;i<=x->nh;i++){
		product->co[i] = 0.0;
	}
	
	c = 1;
	for(i=1;i<=Ai->nh;i++){
		r = Ai->co[i];
		
		product->co[c] += x->co[r] * Ax->co[i];
		if(r < c){
			product->co[r] += x->co[c] * Ax->co[i];
		}else if(r > c){
			printf("r:%ld c:%ld i:%ld Ap[c]:%ld tot:%ld %ld %ld\n",r,c,i,Ap->co[c],Ai->nh,Ap->co[x->nh],x->nh);
			t_error("matrix is not U");
		}
		
		if(i == Ap->co[c]) c ++;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void lower_matrix_product(DOUBLEVECTOR *product, DOUBLEVECTOR *x, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){
	
	long c, r, i;
	
	for(i=1;i<=x->nh;i++){
		product->co[i] = 0.0;
	}
	
	c = 1;
	for(i=1;i<=Ai->nh;i++){
		r = Ai->co[i];
		
		product->co[r] += x->co[c] * Ax->co[i];
		if(r < c){
			printf("r:%ld c:%ld i:%ld Ap[c]:%ld tot:%ld %ld %ld\n",r,c,i,Ap->co[c],Ai->nh,Ap->co[x->nh],x->nh);
			t_error("matrix is not L");
		}
		
		if(i == Ap->co[c]) c ++;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void upper_matrix_product(DOUBLEVECTOR *product, DOUBLEVECTOR *x, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){
	
	long c, r, i;
	
	for(i=1;i<=x->nh;i++){
		product->co[i] = 0.0;
	}
	
	c = 1;
	for(i=1;i<=Ai->nh;i++){
		r = Ai->co[i];
		
		product->co[r] += x->co[c] * Ax->co[i];
		if(r > c){
			printf("r:%ld c:%ld i:%ld Ap[c]:%ld tot:%ld %ld %ld\n",r,c,i,Ap->co[c],Ai->nh,Ap->co[x->nh],x->nh);
			t_error("matrix is not U");
		}
		
		if(i == Ap->co[c]) c ++;
	}
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long BiCGSTAB_unpreconditioned(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *b, 
							   LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){
	
	DOUBLEVECTOR *r0, *r, *p, *v, *s, *t;
	
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j;
	
	r0 = new_doublevector(x->nh);
	r = new_doublevector(x->nh);
	p = new_doublevector(x->nh);
	v = new_doublevector(x->nh);
	s = new_doublevector(x->nh);
	t = new_doublevector(x->nh);
	
	product_using_only_lower_diagonal_part(r, x, Li, Lp, Lx);
	
    for (j=x->nl;j<=x->nh;j++ ) {	
    	r->co[j] = b->co[j] - r->co[j];
		r0->co[j] = r->co[j];
		p->co[j] = 0.;
		v->co[j] = 0.;
    }
	
	norm_r0 = norm_2(r0, r0->nl, r0->nh);
	
	rho = 1.;
	alpha = 1.;
	omeg = 1.;
	
	while ( i<=x->nh && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);
		
		beta = (rho1/rho)*(alpha/omeg);
		
		rho = rho1;
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			p->co[j] = r->co[j] + beta*(p->co[j] - omeg*v->co[j]);
		}
		
		product_using_only_lower_diagonal_part(v, p, Li, Lp, Lx);
		
		alpha = rho/product(r0, v);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			s->co[j] = r->co[j] - alpha*v->co[j];
		}
		
		product_using_only_lower_diagonal_part(t, s, Li, Lp, Lx);
		
		omeg = product(t, s)/product(t, t);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			x->co[j] += (alpha*p->co[j] + omeg*s->co[j]);
			r->co[j] = s->co[j] - omeg*t->co[j];
		}
		
		i++;
				
	}
	
	
	free_doublevector(r0);
	free_doublevector(r);
	free_doublevector(p);
	free_doublevector(v);
	free_doublevector(s);
	free_doublevector(t);
	
	return i;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long BiCGSTAB_diag(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *b, 
				   LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){
	
	DOUBLEVECTOR *r0, *r, *p, *v, *s, *t, *y, *z, *d;
	//DOUBLEVECTOR *ss, *tt;
	
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j, jlim;
	
	r0 = new_doublevector(x->nh);
	r = new_doublevector(x->nh);
	p = new_doublevector(x->nh);
	v = new_doublevector(x->nh);
	s = new_doublevector(x->nh);
	t = new_doublevector(x->nh);
	y = new_doublevector(x->nh);
	z = new_doublevector(x->nh);
	d = new_doublevector(x->nh);
	//ss = new_doublevector(x->nh);
	//tt = new_doublevector(x->nh);
	
	for(j=1;j<=x->nh;j++){
		if(j>1){
			jlim = Lp->co[j-1]+1;
		}else{
			jlim = 1;
		}
		d->co[j] = Lx->co[jlim];
	}
	
	product_using_only_lower_diagonal_part(r, x, Li, Lp, Lx);
	
    for (j=x->nl;j<=x->nh;j++ ) {	
    	r->co[j] = b->co[j] - r->co[j];
		r0->co[j] = r->co[j];
		p->co[j] = 0.;
		v->co[j] = 0.;
    }
	
	norm_r0 = norm_2(r0, r0->nl, r0->nh);
	
	rho = 1.;
	alpha = 1.;
	omeg = 1.;
	
	while ( i<=x->nh && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);
		
		beta = (rho1/rho)*(alpha/omeg);
		
		rho = rho1;
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			p->co[j] = r->co[j] + beta*(p->co[j] - omeg*v->co[j]);
			y->co[j] = p->co[j]/d->co[j];
		}
		
		product_using_only_lower_diagonal_part(v, y, Li, Lp, Lx);
		
		alpha = rho/product(r0, v);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			s->co[j] = r->co[j] - alpha*v->co[j];
			z->co[j] = s->co[j]/d->co[j];
		}
		
		product_using_only_lower_diagonal_part(t, z, Li, Lp, Lx);
		
		omeg = product(t, s)/product(t, t);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			x->co[j] += (alpha*y->co[j] + omeg*z->co[j]);
			r->co[j] = s->co[j] - omeg*t->co[j];
		}
		
		i++;
				
	}
	
	
	free_doublevector(r0);
	free_doublevector(r);
	free_doublevector(p);
	free_doublevector(v);
	free_doublevector(s);
	free_doublevector(t);
	free_doublevector(y);
	free_doublevector(z);
	free_doublevector(d);
	//free_doublevector(ss);
	//free_doublevector(tt);
	
	return i;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void get_diag_lower_matrix(DOUBLEVECTOR *diag, DOUBLEVECTOR *udiag, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){

	long i, ilim;
	
	for(i=1;i<=diag->nh;i++){
		if(i>1){
			ilim = Ap->co[i-1]+1;
		}else{
			ilim = 1;
		}
		diag->co[i] = Ax->co[ilim];
		if(i<diag->nh){
			if(Ai->co[ilim+1] == i+1 && ilim+1 <= Ap->co[i]){
				udiag->co[i] = Ax->co[ilim+1];
			}else{
				udiag->co[i] = 0.0;
			}
		}	
	}
	
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void get_diag_upper_matrix(DOUBLEVECTOR *diag, DOUBLEVECTOR *udiag, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){
	
	long i, ilim;	
	
	for(i=1;i<=diag->nh;i++){
		if(i>1){
			ilim = Ap->co[i-1]+1;
		}else{
			ilim = 1;
		}
		diag->co[i] = Ax->co[Ap->co[i]];
		if(i>1){
			if(Ai->co[Ap->co[i]-1] == i-1 && Ap->co[i]-1 >= ilim){
				udiag->co[i-1] = Ax->co[Ap->co[i]-1];
			}else{
				udiag->co[i-1] = 0.0;
			}
		}	
	}
}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long BiCGSTAB_lower(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *b, 
			 LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){
	
	DOUBLEVECTOR *r0, *r, *p, *v, *s, *t, *diag, *udiag, *y, *z;
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j;
	short sux;
	
	r0 = new_doublevector(x->nh);
	r = new_doublevector(x->nh);
	p = new_doublevector(x->nh);
	v = new_doublevector(x->nh);
	s = new_doublevector(x->nh);
	t = new_doublevector(x->nh);
	diag = new_doublevector(x->nh);
	udiag = new_doublevector(x->nh-1);
	y = new_doublevector(x->nh);
	z = new_doublevector(x->nh);
	
	get_diag_lower_matrix(diag, udiag, Ai, Ap, Ax);
	
	product_using_only_lower_diagonal_part(r, x, Ai, Ap, Ax);
	
    for (j=x->nl;j<=x->nh;j++ ) {	
    	r->co[j] = b->co[j] - r->co[j];
		r0->co[j] = r->co[j];
		p->co[j] = 0.;
		v->co[j] = 0.;
    }
	
	norm_r0 = norm_2(r0, r0->nl, r0->nh);
	
	rho = 1.;
	alpha = 1.;
	omeg = 1.;
	
	while ( i<=x->nh && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);
		
		beta = (rho1/rho)*(alpha/omeg);
		
		rho = rho1;
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			p->co[j] = r->co[j] + beta*(p->co[j] - omeg*v->co[j]);
		}
		
		sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, p, y);
		if(sux==0) return(-1);
		
		product_using_only_lower_diagonal_part(v, y, Ai, Ap, Ax);
		
		alpha = rho/product(r0, v);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			s->co[j] = r->co[j] - alpha*v->co[j];
		}
		
		sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, s, z);	
		if(sux==0) return(-1);

		product_using_only_lower_diagonal_part(t, z, Ai, Ap, Ax);
		
		omeg = product(t, s)/product(t, t);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			x->co[j] += (alpha*y->co[j] + omeg*z->co[j]);
			r->co[j] = s->co[j] - omeg*t->co[j];
		}
		
		i++;
				
	}
	
	
	free_doublevector(r0);
	free_doublevector(r);
	free_doublevector(p);
	free_doublevector(v);
	free_doublevector(s);
	free_doublevector(t);
	free_doublevector(diag);
	free_doublevector(udiag);
	free_doublevector(y);
	free_doublevector(z);
	
	return i;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long BiCGSTAB_upper(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *b, 
					LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax){
	
	DOUBLEVECTOR *r0, *r, *p, *v, *s, *t, *diag, *udiag, *y, *z;
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j;
	short sux;
	
	r0 = new_doublevector(x->nh);
	r = new_doublevector(x->nh);
	p = new_doublevector(x->nh);
	v = new_doublevector(x->nh);
	s = new_doublevector(x->nh);
	t = new_doublevector(x->nh);
	diag = new_doublevector(x->nh);
	udiag = new_doublevector(x->nh-1);
	y = new_doublevector(x->nh);
	z = new_doublevector(x->nh);
	
	get_diag_upper_matrix(diag, udiag, Ai, Ap, Ax);
	
	product_using_only_upper_diagonal_part(r, x, Ai, Ap, Ax);
	
    for (j=x->nl;j<=x->nh;j++ ) {	
    	r->co[j] = b->co[j] - r->co[j];
		r0->co[j] = r->co[j];
		p->co[j] = 0.;
		v->co[j] = 0.;
    }
	
	norm_r0 = norm_2(r0, r0->nl, r0->nh);
	
	rho = 1.;
	alpha = 1.;
	omeg = 1.;
	
	while ( i<=x->nh && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);
		
		beta = (rho1/rho)*(alpha/omeg);
		
		rho = rho1;
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			p->co[j] = r->co[j] + beta*(p->co[j] - omeg*v->co[j]);
		}
		
		sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, p, y);
		if(sux==0) return(-1);
		
		product_using_only_upper_diagonal_part(v, y, Ai, Ap, Ax);
		
		alpha = rho/product(r0, v);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			s->co[j] = r->co[j] - alpha*v->co[j];
		}
		
		sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, s, z);	
		if(sux==0) return(-1);

		product_using_only_upper_diagonal_part(t, z, Ai, Ap, Ax);
		
		omeg = product(t, s)/product(t, t);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			x->co[j] += (alpha*y->co[j] + omeg*z->co[j]);
			r->co[j] = s->co[j] - omeg*t->co[j];
		}
		
		i++;
				
	}
	
	
	free_doublevector(r0);
	free_doublevector(r);
	free_doublevector(p);
	free_doublevector(v);
	free_doublevector(s);
	free_doublevector(t);
	free_doublevector(diag);
	free_doublevector(udiag);
	free_doublevector(y);
	free_doublevector(z);
	
	return i;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void solve_upper_diagonal_system(DOUBLEVECTOR *x, DOUBLEVECTOR *B, LONGVECTOR *Ait, LONGVECTOR *Apt, DOUBLEVECTOR *Axt){
	
	//Ait, Apt, Axt are the transposed matrix (lower diagonal)
	
	long i, c, ilim;
	
	for(c=x->nh;c>=1;c--){
		
		if(c>1){
			ilim = Apt->co[c-1]+1;
		}else{
			ilim = 1;
		}
		
		x->co[c] = B->co[c];
		
		for(i=Apt->co[c];i>ilim;i--){
			x->co[c] -= Axt->co[i] * x->co[Ait->co[i]];
		}
		
		x->co[c] /= Axt->co[ilim];
		
	}
			
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void solve_lower_diagonal_system(DOUBLEVECTOR *x, DOUBLEVECTOR *B, LONGVECTOR *Ait, LONGVECTOR *Apt, DOUBLEVECTOR *Axt){
	
	//Ait, Apt, Axt are the transposed matrix (upper diagonal)
	
	long i, r, ilim;
	
	for(r=1;r<=x->nh;r++){
		
		if(r>1){
			ilim = Apt->co[r-1]+1;
		}else{
			ilim = 1;
		}
		
		x->co[r] = B->co[r];
		
		for(i=ilim;i<Apt->co[r];i++){
			x->co[r] -= Axt->co[i] * x->co[Ait->co[i]];
		}
		
		x->co[r] /= Axt->co[Apt->co[r]];
		
	}
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void solve_SSOR_preconditioning(double w, DOUBLEVECTOR *x, DOUBLEVECTOR *B, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx, 
							    LONGVECTOR *Ui, LONGVECTOR *Up, DOUBLEVECTOR *Ux){
	
	//L D-1 U x  = B
	
	//L y = B
	
	//D-1 U x  = y
	
	//U x = D y
	
	long i;	
	DOUBLEVECTOR *y;
	
	y = new_doublevector(x->nh);

	Lx->co[1] /= w;
	Ux->co[1] /= w;
	for(i=2;i<=x->nh;i++){
		Lx->co[Lp->co[i-1]+1] /= w;
		Ux->co[Up->co[i]] /= w;
	}

	solve_lower_diagonal_system(y, B, Ui, Up, Ux);

	y->co[1] *= Lx->co[1]*(2.-w);
	for(i=2;i<=x->nh;i++){
		y->co[i] *= Lx->co[Lp->co[i-1]+1]*(2.-w);
	}
	
	solve_upper_diagonal_system(x, y, Li, Lp, Lx);

    Lx->co[1] *= w;
	Ux->co[1] *= w;
	for(i=2;i<=x->nh;i++){
		Lx->co[Lp->co[i-1]+1] *= w;
		Ux->co[Up->co[i]] *= w;
	}
					 
	free_doublevector(y);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long BiCGSTAB_LU_SSOR(double w, double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *b, 
					LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx, LONGVECTOR *Ui, LONGVECTOR *Up, DOUBLEVECTOR *Ux){
	
	DOUBLEVECTOR *r0, *r, *p, *v, *s, *t, *y, *z, *ss, *tt;
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j;
	
	r0 = new_doublevector(x->nh);
	r = new_doublevector(x->nh);
	p = new_doublevector(x->nh);
	v = new_doublevector(x->nh);
	s = new_doublevector(x->nh);
	t = new_doublevector(x->nh);
	y = new_doublevector(x->nh);
	z = new_doublevector(x->nh);
	ss = new_doublevector(x->nh);
	tt = new_doublevector(x->nh);
		
	product_using_only_upper_diagonal_part(r, x, Ui, Up, Ux);
	
    for (j=x->nl;j<=x->nh;j++ ) {	
    	r->co[j] = b->co[j] - r->co[j];
		r0->co[j] = r->co[j];
		p->co[j] = 0.;
		v->co[j] = 0.;
    }
	
	norm_r0 = norm_2(r0, r0->nl, r0->nh);
	
	rho = 1.;
	alpha = 1.;
	omeg = 1.;
	
	while ( i<=x->nh && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);
		
		beta = (rho1/rho)*(alpha/omeg);
		
		rho = rho1;
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			p->co[j] = r->co[j] + beta*(p->co[j] - omeg*v->co[j]);
		}
		
		solve_SSOR_preconditioning(w, y, p, Li, Lp, Lx, Ui, Up, Ux);
		product_using_only_upper_diagonal_part(v, y, Ui, Up, Ux);
		
		alpha = rho/product(r0, v);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			s->co[j] = r->co[j] - alpha*v->co[j];
		}
		
		solve_SSOR_preconditioning(w, z, s, Li, Lp, Lx, Ui, Up, Ux);
		product_using_only_upper_diagonal_part(t, z, Ui, Up, Ux);
		
		solve_lower_diagonal_system(tt, t, Ui, Up, Ux);
		solve_lower_diagonal_system(ss, s, Ui, Up, Ux);
		omeg = product(tt, ss)/product(tt, tt);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			x->co[j] += (alpha*y->co[j] + omeg*z->co[j]);
			r->co[j] = s->co[j] - omeg*t->co[j];
		}
		
		i++;
				
	}
	
	
	free_doublevector(r0);
	free_doublevector(r);
	free_doublevector(p);
	free_doublevector(v);
	free_doublevector(s);
	free_doublevector(t);

	free_doublevector(y);
	free_doublevector(z);
	free_doublevector(ss);
	free_doublevector(tt);
	
	return i;
	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/
//FURTHER PERSONALIZATION OF THE FUNCTIONS ABOVE

void product_using_only_strict_lower_diagonal_part(DOUBLEVECTOR *product, DOUBLEVECTOR *x, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){
	
	long c, r, i;
	
	for(i=1;i<=x->nh;i++){
		product->co[i] = 0.0;
	}
	
	c = 1;
	for(i=1;i<=Li->nh;i++){
		
		r = Li->co[i];
		
		if(r > c){
			product->co[c] += Lx->co[i] * (x->co[r] - x->co[c]);
			//printf("c:%ld -> %e i:%ld Lx:%e r:%ld %e c:%ld %e -> %e \n",c,product->co[c],i,Lx->co[i],r,x->co[r],c,x->co[c],Lx->co[i] * (x->co[r] - x->co[c]));

			product->co[r] += Lx->co[i] * (x->co[c] - x->co[r]);
			//printf("r:%ld -> %e i:%ld Lx:%e c:%ld %e r:%ld %e -> %e \n",r,product->co[r],i,Lx->co[i],c,x->co[c],r,x->co[r],Lx->co[i] * (x->co[c] - x->co[r]));

		}else if(r < c){
			printf("r:%ld c:%ld i:%ld Ap[c]:%ld tot:%ld %ld %ld\n",r,c,i,Lp->co[c],Li->nh,Lp->co[x->nh],x->nh);
			t_error("matrix is not L");
		}
		
		if(i<Li->nh){
			while(i >= Lp->co[c]) c++;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(DOUBLEVECTOR *product, DOUBLEVECTOR *x, DOUBLEVECTOR *y,
																	LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){
	
	long i, r, c;
	
	//calculate (A+Iy)*x, A described by its strict lower diagonal part
	
	
	for(i=1;i<=x->nh;i++){
		product->co[i] = y->co[i] * x->co[i];
	}
	
	c = 1;
	for(i=1;i<=Li->nh;i++){
		r = Li->co[i];
		
		if(r > c){
			product->co[c] += Lx->co[i] * (x->co[r] - x->co[c]);
			product->co[r] += Lx->co[i] * (x->co[c] - x->co[r]);
		}else if(r < c){
			printf("r:%ld c:%ld i:%ld Lp[c]:%ld tot:%ld %ld %ld\n",r,c,i,Lp->co[c],Li->nh,Lp->co[x->nh],x->nh);
			t_error("matrix is not L");
		}
		
		if(i<Li->nh){
			while(i >= Lp->co[c]) c++;
		}
	}
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void get_diag_strict_lower_matrix_plus_identity_by_vector(DOUBLEVECTOR *diag, DOUBLEVECTOR *udiag, DOUBLEVECTOR *y,
												   LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){
	
	long i, r, c;
	//find diagonal and upper diagonal of matrix A+Iy, where A is described by its strict lower diagonal part
	
	
	for(i=1;i<=diag->nh;i++){
		diag->co[i] = y->co[i];
		if(i<diag->nh) udiag->co[i] = 0.;
	}
	
	c = 1;
	for(i=1;i<=Li->nh;i++){
		r = Li->co[i];
		
		diag->co[c] -= Lx->co[i];
		diag->co[r] -= Lx->co[i];
			
		if(r == c+1) udiag->co[c] = Lx->co[i];
		
		if(i<Li->nh){
			while(i >= Lp->co[c]) c++;
		}
	}

}	

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x, 
									DOUBLEVECTOR *b, DOUBLEVECTOR *y, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){
	
	//solve sistem (A+Iy)*x = B, find x
	//A M-matrix described by its lower diagonal part
	
	DOUBLEVECTOR *r0, *r, *p, *v, *s, *t, *diag, *udiag, *yy, *z;
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j, maxiter;
	short sux;
	
	r0 = new_doublevector(x->nh);
	r = new_doublevector(x->nh);
	p = new_doublevector(x->nh);
	v = new_doublevector(x->nh);
	s = new_doublevector(x->nh);
	t = new_doublevector(x->nh);
	diag = new_doublevector(x->nh);
	udiag = new_doublevector(x->nh-1);
	yy = new_doublevector(x->nh);
	z = new_doublevector(x->nh);
	
	get_diag_strict_lower_matrix_plus_identity_by_vector(diag, udiag, y, Li, Lp, Lx);

	//product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(r, x, y, Li, Lp, Lx);
	
    for (j=x->nl;j<=x->nh;j++ ) {	
    	//r->co[j] = b->co[j] - r->co[j];
		r->co[j] = b->co[j];
		r0->co[j] = r->co[j];
		p->co[j] = 0.;
		v->co[j] = 0.;
    }
	
	norm_r0 = norm_2(r0, r0->nl, r0->nh);
		
	rho = 1.;
	alpha = 1.;
	omeg = 1.;
	
	maxiter = (long)(x->nh/100.);
	if (maxiter < 100) maxiter = 100;
	
	while ( i<=maxiter && norm_2(r, r->nl, r->nh) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);
		
		beta = (rho1/rho)*(alpha/omeg);
		
		rho = rho1;
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			p->co[j] = r->co[j] + beta*(p->co[j] - omeg*v->co[j]);
		}
		
		sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, p, yy);
		if(sux==0) return(-1);

		product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(v, yy, y, Li, Lp, Lx);
		
		alpha = rho/product(r0, v);
		
		for (j=x->nl;j<=x->nh;j++ ) {		
			s->co[j] = r->co[j] - alpha*v->co[j];
		}
		
		if(norm_2(s,s->nl,s->nh)>1.E-10){
		
			sux=tridiag(0, 0, 0, x->nh, udiag, diag, udiag, s, z);	
			if(sux==0) return(-1);

			product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(t, z, y, Li, Lp, Lx);

			omeg = product(t, s)/product(t, t);
		
			for (j=x->nl;j<=x->nh;j++ ) {		
				x->co[j] += (alpha*yy->co[j] + omeg*z->co[j]);
				r->co[j] = s->co[j] - omeg*t->co[j];
			}
			
		}else{
			
			for (j=x->nl;j<=x->nh;j++ ) {		
				x->co[j] += alpha*yy->co[j];
				r->co[j] = s->co[j];
			}	
			
		}
		
		i++;
				
	}
	
	free_doublevector(r0);
	free_doublevector(r);
	free_doublevector(p);
	free_doublevector(v);
	free_doublevector(s);
	free_doublevector(t);
	free_doublevector(diag);
	free_doublevector(udiag);
	free_doublevector(yy);
	free_doublevector(z);
	
	return i;
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void product_matrix_using_lower_part_by_vector_plus_vector(double k, DOUBLEVECTOR *out, DOUBLEVECTOR *y, DOUBLEVECTOR *x, 
														   LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){
	
	//calculates k*(y + Ax), where k is coefficient, y and x vectors, and A a SPD matrix defined with its lower diagonal part
	
	long i;
	
	product_using_only_strict_lower_diagonal_part(out, x, Li, Lp, Lx);

	for(i=1;i<=x->nh;i++){
		//printf("-> i:%ld k:%e out:%e y:%e out:%e\n",i,k,out->co[i],y->co[i],k * (out->co[i] + y->co[i]));
		out->co[i] = k * (out->co[i] + y->co[i]);
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


