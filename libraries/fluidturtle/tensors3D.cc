#include "tensors3D.h"

/* Note that depth is the first indices and that the indices were perutated
with respect to NR */

/*-----------------------------------------------------------------------*/


double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) t_error("allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) t_error("allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) t_error("allocation failure 3 in d3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	/* return pointer to array of pointers to rows */
	return t;
}


/*-----------------------------------------------------------------------*/


DOUBLETENSOR *new_doubletensor(long ndh,long nrh,long nch)
{
	DOUBLETENSOR *m;
	m=(DOUBLETENSOR *)malloc(sizeof(DOUBLETENSOR));
	if (!m) t_error("allocation failure in new_doubletensor()");
	m->isdynamic=isDynamic;
	m->nrl=NL;
	m->nrh=nrh;
	m->ncl=NL;
	m->nch=nch;
	m->ndl=NL;
	m->ndh=ndh;
	m->co=d3tensor(m->ndl,m->ndh,m->nrl,m->nrh,m->ncl,m->nch);
	return m;
}

DOUBLETENSOR *new_doubletensor0(long ndh,long nrh,long nch)
{
	DOUBLETENSOR *m;	
	m=(DOUBLETENSOR *)malloc(sizeof(DOUBLETENSOR));
	if (!m) t_error("allocation failure in new_doubletensor()");
	m->isdynamic=isDynamic;
	m->nrl=NL;
	m->nrh=nrh;
	m->ncl=NL;
	m->nch=nch;
	m->ndl=0;
	m->ndh=ndh;
	m->co=d3tensor(m->ndl,m->ndh,m->nrl,m->nrh,m->ncl,m->nch);
	return m;
}


/*-----------------------------------------------------------------------*/


void free_d3tensor(double ***t, long nrl, long ncl, long ndl)
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}


/*-----------------------------------------------------------------------*/


void free_doubletensor( DOUBLETENSOR *m)

{
		  if(m==NULL || m->co==NULL){
			  	t_error("This matrix was never allocated");
		}else if(m->isdynamic==1){

			free_d3tensor(m->co,m->ndl,m->nrl,m->ncl);
			m->isdynamic=m->nrl=m->ncl=m->nrh=m->nch=m->ndl=m->ndh=-1;
			free(m);

			return;

		  }else{
			printf("\nWarning::An attemp was made to free a non dynamic tensor\n");
		  }
}



/*---------------------------------------------------------------------------*/
void initialize_doubletensor(DOUBLETENSOR *L, double sign)

{

long i,j,k;

if(L!=NULL){
	if(L->isdynamic==1){
		for(k=L->ndl;k<=L->ndh;k++){
			for(i=L->nrl;i<=L->nrh;i++){
				for(j=L->ncl;j<=L->nch;j++){			
					L->co[k][i][j]=sign;
				}
			}
		}
	}else{
		t_error("This tensor was no properly allocated");
	}
}else{
	t_error("A null tensor was addressed");
}
}


void copy_doubletensor(DOUBLETENSOR *origin,DOUBLETENSOR *destination)
// added by Emanuele e Matteo for netCDF
{
  long i,j, l;
  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){
    t_error("A tensor was not allocated");
  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->ndh <1 || destination->ndh <1 || origin->nrh <1 || destination->nrh <1 ||  origin->nch <1 || destination->nch <1 ){
    t_error("A tensor was not allocated properly");
  }else if( origin->nrh != destination->nrh ||  origin->nch != destination->nch ||  origin->ndh != destination->ndh ){
    t_error("The tensors do not have the same dimensions");
  }

  for(i=1;i<=origin->nrh;i++){
    for(j=1;j<=origin->nch;j++){
    	for(l=origin->ndl;l<=origin->ndh;l++){
    		destination->co[l][i][j]=origin->co[l][i][j];
    	}
    }
  }
}




DOUBLETENSOR *new_doubletensor_flexlayer(long ndl,long ndh,long nrh,long nch)
// added by Emanuele e Matteo for netCDF
{
	DOUBLETENSOR *m;
	m=(DOUBLETENSOR *)malloc(sizeof(DOUBLETENSOR));
	if (!m) t_error("allocation failure in new_doubletensor()");
	m->isdynamic=isDynamic;
	m->nrl=NL;
	m->nrh=nrh;
	m->ncl=NL;
	m->nch=nch;
	m->ndl=ndl;
	m->ndh=ndh;
	m->co=d3tensor(m->ndl,m->ndh,m->nrl,m->nrh,m->ncl,m->nch);
	return m;
}

/*===============functions copied from utilities.c ===================*/

void stop_execution(void)

{

char ch;

printf("\nPRESS RETURN TO CONTINUE\n");
scanf("%c",&ch);

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

//short tridiag2(short a, long r, long c, long nbeg, long nend, DOUBLEVECTOR *ld, DOUBLEVECTOR *d, DOUBLEVECTOR *ud, DOUBLEVECTOR *b, DOUBLEVECTOR *e)
short tridiag2(short a, long r, long c, long nbeg, long nend, const GeoVector<double>& ld, const GeoVector<double>& d, const GeoVector<double>& ud, const GeoVector<double>& b, GeoVector<double>& e)

//solve A(ld,d,ud) * e + b = 0

{
	long j;
	double bet;
	size_t lIndex = 0 ;
	//DOUBLEVECTOR *gam;
    GeoVector<double> gam;

	//gam = new_doublevector(nend);
    gam.resize(nend+1);

	//bet = d->co[nbeg];
	bet = d[nbeg];
	if(bet == 0.0){
		return 1;
	}
	//e->co[nbeg] = -b->co[nbeg]/bet;
	e[nbeg] = -b[nbeg]/bet;

	//Decomposition and forward substitution
	for(j=nbeg+1; j<=nend; j++){
		//gam->co[j] = ud->co[j-1]/bet;
		gam[j] = ud[j-1]/bet;
		//bet = d->co[j]-ld->co[j-1]*gam->co[j];
		bet = d[j]-ld[j-1]*gam[j];
		if(bet == 0.0){
			return 1;
		}
		//e->co[j] = (-b->co[j]-ld->co[j-1]*e->co[j-1])/bet;
		e[j] = (-b[j]-ld[j-1]*e[j-1])/bet;
	}

	//Backsubstitution
	for(j=(nend-1); j>=nbeg; j--){
		//e->co[j] -= gam->co[j+1]*e->co[j+1];
		e[j] -= gam[j+1]*e[j+1];
	}

	//free_doublevector(gam);

#ifdef VERBOSE
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

//double norm_inf(DOUBLEVECTOR *V, long nbeg, long nend){
double norm_inf(const GeoVector<double>& V, long nbeg, long nend){

	long l;
	double N=0.0;

	for(l=nbeg;l<=nend;l++){
	//	if (fabs(V->co[l])> N) N = fabs(V->co[l]);
		if (fabs(V[l])> N) N = fabs(V[l]);
	}

	return(N);

}

/*----------------------------------------------------------------------------------------------------------*/

//double norm_2(DOUBLEVECTOR *V, long nbeg, long nend){
double norm_2(const GeoVector<double>& V, long nbeg, long nend){

	long l;
	double N=0.0;

	for(l=nbeg;l<nend;l++){
		//N+=pow(V->co[l],2.0);
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





