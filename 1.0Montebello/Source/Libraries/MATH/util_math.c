
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */

#include "turtle.h"
#include "util_math.h"
#include "t_utilities.h"
#include "constant.h"

/*----------------------------------------------------------------------------------------------------------*/

void tridiag(short a, long r, long c, long nx, DOUBLEVECTOR *diag_inf, DOUBLEVECTOR *diag, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e)

{

long j;
double bet;
DOUBLEVECTOR *gam;

gam=new_doublevector(nx);

if(diag->co[1]==0.0){
	printf("type=%d r=%ld c=%ld\n",a,r,c);
	t_error("Error 1 in tridiag");
}

bet=diag->co[1];
e->co[1]=b->co[1]/bet;

//Decomposition and forward substitution
for(j=2;j<=nx;j++){
	gam->co[j]=diag_sup->co[j-1]/bet;
	bet=diag->co[j]-diag_inf->co[j-1]*gam->co[j];
	if(bet==0.0){
		printf("type=%d r=%ld c=%ld\n",a,r,c);
		printf("l=%ld diag(l)=%f diag_inf(l-1)=%f diag_sup(l-1)=%f\n",j,diag->co[j],diag_inf->co[j-1],diag_sup->co[j-1]);
		t_error("Error 2 in tridiag");
	}
	e->co[j]=(b->co[j]-diag_inf->co[j-1]*e->co[j-1])/bet;
}

//Backsubstitution
for(j=(nx-1);j>=1;j--){
	e->co[j]-=gam->co[j+1]*e->co[j+1];
}

free_doublevector(gam);

}

/*----------------------------------------------------------------------------------------------------------*/


void tridiag2(short a, long r, long c, long nx, double wdi, DOUBLEVECTOR *diag_inf, double wd, DOUBLEVECTOR *diag, double wds, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e)

//solve A(wdi*diag_inf,wd*diag,wds*diag_sup) * e + b = 0

{
	
	long j;
	double bet;
	DOUBLEVECTOR *gam;
	
	gam=new_doublevector(nx);
	
	if(wd*diag->co[1]==0.0){
		printf("type=%d r=%ld c=%ld\n",a,r,c);
		t_error("Error 1 in tridiag");
	}
	
	bet=wd*diag->co[1];
	e->co[1]=-b->co[1]/bet;
	
	//Decomposition and forward substitution
	for(j=2;j<=nx;j++){
		gam->co[j]=wds*diag_sup->co[j-1]/bet;
		bet=wd*diag->co[j]-wdi*diag_inf->co[j-1]*gam->co[j];
		if(bet==0.0){
			printf("type=%d r=%ld c=%ld\n",a,r,c);
			printf("l=%ld diag(l)=%f diag_inf(l-1)=%f diag_sup(l-1)=%f\n",j,wd*diag->co[j],wdi*diag_inf->co[j-1],wds*diag_sup->co[j-1]);
			t_error("Error 2 in tridiag");
		}
		e->co[j]=(-b->co[j]-wdi*diag_inf->co[j-1]*e->co[j-1])/bet;
	}
	
	//Backsubstitution
	for(j=(nx-1);j>=1;j--){
		e->co[j]-=gam->co[j+1]*e->co[j+1];
	}
	
	free_doublevector(gam);
	
}

/*----------------------------------------------------------------------------------------------------------*/

double norm_2(DOUBLEVECTOR *V, long n){

	long l;
	double N=0.0;
	
	for(l=1;l<=n;l++){
		N+=pow(V->co[l],2.0);
	}
	
	N=sqrt(N);
	return(N);		

}

/*----------------------------------------------------------------------------------------------------------*/

double norm_inf(DOUBLEVECTOR *V, long n){
	
	long l;
	double N=0.0;
	
	for(l=1;l<=n;l++){
		if( fabs(V->co[l]) > N ) N = fabs(V->co[l]);
	}
	
	return(N);		
	
}
/*----------------------------------------------------------------------------------------------------------*/

double sum(DOUBLEVECTOR *V, long n){
	
	long l;
	double N=0.0;
	
	for(l=1;l<=n;l++){
		N += V->co[l];
	}
	
	return(N);		
	
}

/*----------------------------------------------------------------------------------------------------------*/

void Cramer_rule(double A, double B, double C, double D, double E, double F, double *x, double *y){

	/*
	Ax + By = C
	Dx + Ey = F

	x = (CE - FB) / (AE - DB)
	y = (AF - CD) / (AE - DB)
	*/
	
	*x = (C*E - F*B) / (A*E - D*B);
	*y = (A*F - C*D) / (A*E - D*B);
	
}

/*----------------------------------------------------------------------------------------------------------*/


double minimize_merit_function(double res0, double lambda1, double res1, double lambda2, double res2) {
	
	double lambda;
	double a, b, c;			//interpolating polynomial: ax2 + bx + c
	
	//calculate three-point quadratic polynomial interpolating the merit function
	c = res0;
	Cramer_rule(pow(lambda1, 2.0), lambda1, res1-res0, pow(lambda2, 2.0), lambda2, res2-res0, &a, &b);
	
	//minimize ax^2+bx+c
	if(a>0){
		lambda = -b/(2*a);
		if(lambda < lambda1*thmin){
			lambda = lambda1*thmin;
		}else if(lambda > lambda1*thmax){
			lambda = lambda1*thmax;
		}
	}else{
		if(a * lambda1*thmin*lambda1*thmin + b * lambda1*thmin + c < a * lambda1*thmax*lambda1*thmax + b * lambda1*thmax + c){
			lambda = lambda1*thmin;
		}else{
			lambda = lambda1*thmax;
		}
	}
	
	return(lambda);
}

/*----------------------------------------------------------------------------------------------------------*/

double product(DOUBLEVECTOR *a, DOUBLEVECTOR *b){

	double p=0.;
	long i,n=a->nh;
	
	for(i=1;i<=n;i++){
		p += a->co[i] * b->co[i];
	}
	
	return(p);
	
}

/*----------------------------------------------------------------------------------------------------------*/
