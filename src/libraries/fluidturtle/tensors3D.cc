#include "tensors3D.h"

/*===============functions copied from utilities.c ===================*/

void stop_execution(void){
	
	char ch;

	printf("\nPRESS RETURN TO CONTINUE\n");
	scanf("%c",&ch);//this generates a warning
	
}


double Fmin(double a, double b){

	double min=a;

	if(b<a) min=b;

	return(min);

}

long Fminlong(long a, long b){

	long min=a;

	if(b<a) min=b;

	return(min);

}


double Fmax(double a, double b){

	double max=a;

	if(b>a) max=b;

	return(max);

}

long Fmaxlong(long a, long b){

	long max=a;

	if(b>a) max=b;

	return(max);

}


/*===============================functions copied from util_math.c========================================*/

short tridiag2(long nbeg, long nend, const GeoVector<double>& ld, const GeoVector<double>& d, const GeoVector<double>& ud, const GeoVector<double>& b, GeoVector<double>& e)
//solve A(ld,d,ud) * e + b = 0
{
	long j;
	double bet;
#ifdef VERYVERBOSE
	size_t lIndex = 0 ;
#endif
    GeoVector<double> gam;

    gam.resize(nend+1);

	bet = d[nbeg];
	if(bet == 0.0){
		return 1;
	}

	e[nbeg] = -b[nbeg]/bet;

	//Decomposition and forward substitution
	for(j=nbeg+1; j<=nend; j++){
		gam[j] = ud[j-1]/bet;
		bet = d[j]-ld[j-1]*gam[j];
		if(bet == 0.0){
			return 1;
		}
		e[j] = (-b[j]-ld[j-1]*e[j-1])/bet;
	}

	//Backsubstitution
	for(j=(nend-1); j>=nbeg; j--){
		e[j] -= gam[j+1]*e[j+1];
	}

#ifdef VERYVERBOSE
	printf("DEBUG_PRINT: nbeg(%ld),nend(%ld)\n", nbeg, nend);
	printf("DEBUG_PRINT: ld(");
	for(lIndex = 0 ; lIndex < ld.size() ; lIndex ++)
	{
		printf("(%lu:%.12g)",lIndex, ld[lIndex]);
	}
	printf(")\n");
	printf("DEBUG_PRINT: d(");
	for(lIndex = 0 ; lIndex < d.size() ; lIndex ++)
	{
		printf("(%lu:%.12g)",lIndex, d[lIndex]);
	}
	printf(")\n");
	printf("DEBUG_PRINT: b(");
	for(lIndex = 0 ; lIndex < b.size() ; lIndex ++)
	{
		printf("(%lu:%.12g)",lIndex, b[lIndex]);
	}
	printf(")\n");
	printf("DEBUG_PRINT: e(");
	for(lIndex = 0 ; lIndex < e.size() ; lIndex ++)
	{
		printf("(%lu:%.12g)",lIndex, e[lIndex]);
	}
	printf(")\n");
#endif
	return 0;
}

/*----------------------------------------------------------------------------------------------------------*/

double norm_inf(const GeoVector<double>& V, long nbeg, long nend){

	long l;
	double N=0.0;

	for(l=nbeg;l<=nend;l++){
		if (fabs(V[l])> N) N = fabs(V[l]);
	}

	return(N);

}

/*----------------------------------------------------------------------------------------------------------*/


double norm_2(const GeoVector<double>& V, long nbeg, long nend){

	long l;
	double N=0.0;

	for(l=nbeg;l<nend;l++){
		  N+=pow(V[l],2.0);
	}
	N=sqrt(N);

	return(N);

}



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
		if(lambda < lambda1*GTConst::thmin){
			lambda = lambda1*GTConst::thmin;
		}else if(lambda > lambda1*GTConst::thmax){
			lambda = lambda1*GTConst::thmax;
		}
	}else{
		if(a * lambda1*GTConst::thmin*lambda1*GTConst::thmin + b * lambda1*GTConst::thmin + c < a * lambda1*GTConst::thmax*lambda1*GTConst::thmax + b * lambda1*GTConst::thmax + c){
			lambda = lambda1*GTConst::thmin;
		}else{
			lambda = lambda1*GTConst::thmax;
		}
	}

	return(lambda);
}


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

