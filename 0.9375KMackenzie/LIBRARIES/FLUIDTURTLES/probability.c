#include "turtle.h"
#include "t_probability.h"
#define PI 3.141592654
/*From the Numerical Recipes Second Ediction pag.287*/
double expdev(long *idum)
{
float dum;
do{
	dum=ran3(idum);
}while(dum==0.0);
return -log(dum);
}

/*From the Numerical Recipes Second Ediction pag.294*/
double poisdev(float xm, long *idum)
{
static float sq,alxm,G,oldm=-(1.0);
float em,t,y;

if(xm<12.0){
     if(xm!=oldm){
     	oldm=xm;
     	G=exp(-xm);
     }
     em=-1;
     t=1.0;
     do{
     	++em;
      	t*=ran3(idum);
     }while(t>G);
}else{
     if(xm!=oldm){
     	oldm=xm;
     	sq=sqrt(2.0*xm);
     	alxm=log(xm);
     	G=xm*alxm-gammln(xm+1.0);
     }
     do{
     	do{
	     y=tan(PI*ran3(idum));
	     em=sq*y+xm;
     	}while(em<0.0);
     	em=floor(em);
     	t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-G);
     }while(ran3(idum)>t);
}
return em;
}
/*from the Numerical Recepis*/
/* RETURNS AS A FLOATING POINT NUMBER AN INTEGER VALUE THAT IS
A RANDOM DEVIATE DRAWN FROM A BINOMIAL DISTRIBUTION OF n
TRIALS EACH OF PROBABILITY PP - see Numerical Rec. */
/*--------------------------------------------------------------*/
float bnldev(float pp, int n, long *idum)
{

long j;
long nold=(-1);
float am, em,G,angle,p,bnl,sq,t,y;
static float pold=(-1.0),pc,plog,pclog,en,oldg;

p=(pp<=0.5 ? pp : 1.0-pp);
am=n*p;
if(n<25){
    bnl=0.0;
    for(j=1;j<=n;j++)
    if(ran3(idum)<p) ++bnl;
}else if (am<1.0){
    G=exp(-am);
    t=1.0;
    for(j=0;j<=n;j++){
	t *=ran3(idum);
	if(t<G) break;
    }
    bnl=(j<=n ? j : n);
}else{
    if(n!=nold){
	en=n;
	oldg=gammln(en+1.0);
	nold=n;
    }if (p!=pold){
    	pc=1.0-p;
    	plog=log(p);
    	pclog=log(pc);
    	pold=p;
    }
    sq=sqrt(2.0*am*pc);
    do{
    	do{
    	    angle=PI*ran3(idum);
	    y=tan(angle);
    	    em=sq*y+am;
    	} while(em < 0.0 || em>= (en+1.0));
    	em=floor(em);
    	t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
    }while(ran3(idum)>t);
    bnl=em;
}
if(p!=pp) bnl=n-bnl;
return bnl;
}
/*From the Numerical recepis pag.214:
   Returns the value ln[gamma(xx)] for x>0*/
float gammln(float xx)
{
double x,y,tmp,ser;
static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,
  0.1208650973866179e-2,-0.5395239384953e-5};
long j;
y=x=xx;
tmp=x+5.5;
tmp-=(x+0.5)*log(tmp);
ser=1.000000000190015;
for(j=0;j<=5;j++)ser+=cof[j]/++y;
return -tmp+log(2.5066282746310005*ser/x);
}

/* ... Double precision variables matrices' mean estimating function 
         V rappresenta i novalues*/
double meandoublem(DOUBLEMATRIX *m,FLOATVECTOR *V)
{
short segna;
long i,j,k,NN;
double mm=0;
NN=m->nrh*m->nch;
for(i=1;i<=m->nrh;i++){
for(j=1;j<m->nch;j++){
	segna=0;
	for(k=1;k<=V->nh;k++){
	    if(m->co[i][j]==V->co[k])segna=1;
	}
	if(segna==0)mm+=m->co[i][j];
}
}
mm/=(double)NN;
return mm;
}

/* ... Double precision variables matrices' variance estimating function
        V rappresenta i novalues*/
double vardoublem(DOUBLEMATRIX *m,FLOATVECTOR *V)
{
short segna;
long i,j,k;
double mm=0,mn=0,mx;
mx=(m->nch*m->nrh);
for(i=1;i<=m->nrh;i++){
for(j=1;j<=m->nch;j++) {
	segna=0;
	for(k=1;k<=V->nh;k++){
	    if(m->co[i][j]==V->co[k])segna=1;
	}
	if(segna==0){
	    mn+=m->co[i][j]*m->co[i][j];
	    mm+=m->co[i][j];
	}
}
}
mm/=mx;
return mn/(mx-1)-mx*mm*mm/(mx-1);
}
