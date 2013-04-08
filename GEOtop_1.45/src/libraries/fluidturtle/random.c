#include "t_random.h"

/*--------------------------------------------------------------------------*/ 
/* Quick random number generator modified from Numerical Recipes */

long urand(long *idum,long range)
{
static long y, maxran, v[98];
long dum;
static long iff=0;
long j=0;
long i;

if(*idum <0 || iff==0){
iff=1; i=2;

maxran=RAND_MAX+1.0;
srand(*idum);
*idum=1;
for(i=1;j<=97;j++) dum=rand();
for(i=1;i<=97;i++) v[i]=range*(long)rand()/maxran;

y=1+97*(long)rand()/maxran;

}
j=y;
y=v[j];
v[j]=range*(long)rand()/maxran;

return y;
}
/*--------------------------------------------------------------------------*/
/* Another random number generator from Numerical Recipes  with some
modifications                                                               */

double ran1(long *idum)
{
static long ix1,ix2,ix3;
static double r[98];
double temp;
static long iff;
long j;

if(*idum<0 || iff==0){
	iff=1;
	ix1=(IC1-(*idum)) %  M1;
	ix1=(IA1*ix1+IC1) % M1;
    ix2=ix1 % M2;
    ix1=(IA1*ix1+IC1) % M1;
    ix3=ix1 % M3;
    for(j=0;j<=96;j++){
    	ix1=(IA1*ix1+IC1) % M1;
    	ix2=(IA2*ix2+IC2) % M2;
    	r[j]=(ix1+ix2*RM2)*RM1;
	}
	*idum=1;
}

ix1=(IA1*ix1+IC1) % M1;
ix2=(IA2*ix2+IC2) % M2;
ix3=(IA3*ix3+IC3) % M3;
j=1+((96*ix3)/M3);
if(j > 96 || j<0) t_error("This cannot happen.");
temp=r[j];
r[j]=(ix1+ix2*RM2)*RM1;
return temp;

}


/* As above a random generator from NR */

double ran2(long *idum)

{


long j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
double temp;


if(*idum <=0){
	if(-(*idum) <1) *idum=1;
	else *idum=-(*idum);
	idum2=(*idum);
	for(j=NTAB+7;j>=0;j--){
		k=(*idum)/IQ1;
		*idum=IAA1*(*idum-k*IQ1)-k*IR1;
		if (*idum <0) *idum +=IM1;
		if ( j < NTAB) iv[j] =* idum;
		
	}
    iy=iv[0];
}

k=(*idum)/IQ1;
*idum=IAA1*(*idum-k*IQ1)-k*IR1;
if(*idum <0) *idum += IM1;
k=idum2/IQ2;
idum2=IAA2*(idum2-k*IQ2)-k*IR2;
if(idum2 < 0) idum += IM2;
j=iy/NDIV;
iy=iv[j]-idum2;
iv[j] = *idum;
if( iy <1 ) iy +=IMM1;

/* printf("rnd-> %f %f\n",RNMX,temp);*/
if ((temp=AM1*iy) > RNMX) return RNMX;
else return temp;

}

double ran3( long *idum)
/* "Minimal" random number generator of Park and Miller...
see Press et al."Numerical Recipes in C", 1992. */
{ int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;

if(*idum<=0 || !iy){
if(-(*idum)<1) *idum=1;
else *idum= -(*idum);
for (j=NTAB+7;j>=0; j--){
k=(*idum)/IQ;
*idum=IA*(*idum-k*IQ)-IR*k;
if(*idum<0) *idum +=IM;
if(j<NTAB) iv[j]=*idum;
}iy=iv[0];
}
k=(*idum)/IQ;
*idum=IA*(*idum-k*IQ)-IR*k;
if (*idum<0) *idum+=IM;
j=iy/NDIV;
iy=iv[j];
iv[j]=*idum;
if((temp=AM*iy)>RNMX)return RNMX;
else return temp;
}

//

/* commentata per  warning: integer overflow in expression
void prand(long *pto, long row, long col,long toggle )

{

long i;
static long vx[98],vy[99];
long dnr;
if(toggle==0){
		srand(time(NULL));
		for(i=0;i<=97;i++) rand();
		for(i=0;i<=97;i++) vx[i]=(row-2)*(long)rand()/(RAND_MAX+1);
		for(i=0;i<=97;i++) vy[i]=(col-2)*(long)rand()/(RAND_MAX+1);
		toggle=1;
		}

dnr=98*(long)rand()/(RAND_MAX+1);
if(dnr>97) {printf("dnr was %d",dnr); t_error("Something is wrong here");}
pto[0]=vx[dnr]+2;
vx[dnr]=(row-2)*(long)rand()/(RAND_MAX+1);
dnr=98*(long)rand()/(RAND_MAX+1);
if(dnr>97) {printf("dnr was %d",dnr); t_error("Something is wrong here");}
pto[1]=vy[dnr]+2;
vy[dnr]=(col-2)*(long)rand()/(RAND_MAX+1);

}


//

// commentata per warning: integer overflow in expression
long rrand(long row )

{
long i;
static double vx[98];
static long toggle;
long pto,dnr;
//char ch;
//Initializing the randon number generator
//commentata per  warning: integer overflow in expression
if(toggle==0){
		srand(time(NULL));
		for(i=0;i<=97;i++) rand();
		for(i=0;i<=97;i++) vx[i]=(double)rand()/(RAND_MAX+1);
		toggle=1;
		}

dnr=98*(long)rand()/(RAND_MAX+1);
if(dnr>97) {printf("dnr was %d",dnr); t_error("merda");}
pto=abs(row*vx[dnr]);
vx[dnr]=(double)rand()/(RAND_MAX+1);

// commentata per  warning: integer overflow in expression
return pto;

}


// generator of random variables with lognormal distribution, average vmed, standard deviation vvar
  */


double uuvel(float vmed,float vvar)
{
static long idum=3;
float umed, uvar,u,v;
double y;
uvar=log(1+pow((vvar/vmed),2));
umed=0.5*log(pow(vmed,2)/(1+pow((vvar/vmed),2)));
y=gasdev(&idum);
u=umed+y*sqrt(uvar);
v=exp(u);
return v;
}



/*--------------------------------------------------------------------------*/
double gasdev(long *idum)

/* returns a normally distributed random deviate with zero mean mean
and unit variance, using rann1(idum) as source of uniform deviates */
{

/* double ran1(long *idum); */
static int iset=0;
static double gset;
double fac,rsq,v1,v2;

if(iset==0){

/* we don't have an extra deviate handy, so pick two uniform numbers in
the square extending from -1 to +1 in each direction, see if they are in
the unit circle and if they are not, try again */

do{
v1=2.0*ran1(idum)-1.0;
v2=2.0*ran1(idum)-1.0;
/*printf("v1=%f,v2=%f",(v1+1.0)/2.0,(v2+1.0)/2.0);*/
rsq=v1*v1+v2*v2;
} while (rsq>=1.0 || rsq==0.0);
fac=sqrt(-2.0*log(rsq)/rsq);

/* now make the Box-Muller transformation to get two normal deviates.
return one and save the other for the next time. */

gset=v1*fac;

/* set flag */

iset=1;
return v2*fac;
}else{

/* we have an extra deviate handy, so unset the flag and return it */

iset=0;
return gset;
}
}



