
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 - Version 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145
 
 GEOtop 1.145 is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "../fluidturtle/turtle.h"
#include "util_math.h"
#include "../fluidturtle/t_utilities.h"
#include "../../geotop/constants.h"

/*----------------------------------------------------------------------------------------------------------*/

short tridiag(short a, long r, long c, long nx, DOUBLEVECTOR *diag_inf, DOUBLEVECTOR *diag, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e)

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
		printf("Error 2 in tridiag");
		return 0;
	}
	e->co[j]=(b->co[j]-diag_inf->co[j-1]*e->co[j-1])/bet;
}

//Backsubstitution
for(j=(nx-1);j>=1;j--){
	e->co[j]-=gam->co[j+1]*e->co[j+1];
}

free_doublevector(gam);
	
return 1;

}

/*----------------------------------------------------------------------------------------------------------*/


void tridiag2(short a, long r, long c, long nbeg, long nend, DOUBLEVECTOR *ld, DOUBLEVECTOR *d, DOUBLEVECTOR *ud, DOUBLEVECTOR *b, DOUBLEVECTOR *e)

//solve A(ld,d,ud) * e + b = 0

{
	
	long j;
	double bet;
	DOUBLEVECTOR *gam;
	
	gam = new_doublevector(nend);
		
	bet = d->co[nbeg];
	if(bet == 0.0){
		printf("type=%d r=%ld c=%ld\n",a,r,c);
		t_error("Error 1 in tridiag");
	}	
	e->co[nbeg] = -b->co[nbeg]/bet;
	
	//Decomposition and forward substitution
	for(j=nbeg+1; j<=nend; j++){
		gam->co[j] = ud->co[j-1]/bet;
		bet = d->co[j]-ld->co[j-1]*gam->co[j];
		if(bet == 0.0){
			printf("type=%d r=%ld c=%ld\n",a,r,c);
			printf("l=%ld d(l)=%f ld(l-1)=%f ud(l-1)=%f\n",j,d->co[j],ld->co[j-1],ud->co[j-1]);
			t_error("Error 2 in tridiag");
		}
		e->co[j] = (-b->co[j]-ld->co[j-1]*e->co[j-1])/bet;
	}
	
	//Backsubstitution
	for(j=(nend-1); j>=nbeg; j--){
		e->co[j] -= gam->co[j+1]*e->co[j+1];
	}
	
	free_doublevector(gam);
	
}

/*----------------------------------------------------------------------------------------------------------*/

double norm_inf(DOUBLEVECTOR *V, long nbeg, long nend){
	
	long l;
	double N=0.0;
	
	for(l=nbeg;l<=nend;l++){
		if (fabs(V->co[l])> N) N = fabs(V->co[l]);
	}
	
	return(N);		
	
}

/*----------------------------------------------------------------------------------------------------------*/

double norm_2(DOUBLEVECTOR *V, long nbeg, long nend){

	long l;
	double N=0.0;
	
	for(l=nbeg;l<=nend;l++){
		N+=pow(V->co[l],2.0);
	}
	N=sqrt(N);

	return(N);		

}

/*----------------------------------------------------------------------------------------------------------*/

double norm_1(DOUBLEVECTOR *V, long nbeg, long nend){
	
	long l;
	double N=0.0;
	
	for(l=nbeg;l<=nend;l++){
		N += fabs(V->co[l]);
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

//
// Recursive auxiliary function for adaptiveSimpsons() function below
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//                                                                                                 
double adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon,                 
						   double S, double fa, double fb, double fc, int bottom) {                 
	double c = (a + b)/2, h = b - a;                                                                  
	double d = (a + c)/2, e = (c + b)/2;                                                              
	double fd = f(d), fe = f(e);                                                                      
	double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
	double Sright = (h/12)*(fc + 4*fe + fb);                                                          
	double S2 = Sleft + Sright;                                                                       
	if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)                                                    
		return S2 + (S2 - S)/15;                                                                        
	return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
	adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
}         

//
// Adaptive Simpson's Rule
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsons(double (*f)(double),   // ptr to function
						double a, double b,  // interval [a,b]
						double epsilon,  // error tolerance
						int maxRecursionDepth) {   // recursion cap        
	double c = (a + b)/2, h = b - a;                                                                  
	double fa = f(a), fb = f(b), fc = f(c);                                                           
	double S = (h/6)*(fa + 4*fc + fb);                                                                
	return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
}

/*----------------------------------------------------------------------------------------------------------*/

//
// Recursive auxiliary function for adaptiveSimpsons() function below
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//                                                                                                 
double adaptiveSimpsonsAux2(double (*f)(double x, void *p), void *arg, double a, double b, double epsilon, 
							double S, double fa, double fb, double fc, int bottom) {                 
	double c = (a + b)/2, h = b - a;                                                                  
	double d = (a + c)/2, e = (c + b)/2;                                                              
	double fd = f(d, arg), fe = f(e, arg);                                                                      
	double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
	double Sright = (h/12)*(fc + 4*fe + fb);                                                          
	double S2 = Sleft + Sright;                                                                       
	if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)                                                    
		return S2 + (S2 - S)/15;                                                                        
	return adaptiveSimpsonsAux2(f, arg, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
		   adaptiveSimpsonsAux2(f, arg, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
}         

//
// Adaptive Simpson's Rule
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsons2(double (*f)(double x, void *p), void *arg,   // ptr to function
						 double a, double b,  // interval [a,b]
						 double epsilon,  // error tolerance
						 int maxRecursionDepth) {   // recursion cap        
	double c = (a + b)/2, h = b - a;                                                                  
	double fa = f(a, arg), fb = f(b, arg), fc = f(c, arg);                                                           
	double S = (h/6)*(fa + 4*fc + fb);                                                                
	return adaptiveSimpsonsAux2(f, arg, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
}

/*----------------------------------------------------------------------------------------------------------*/

