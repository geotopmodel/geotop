
/* MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file util_math.c

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


/*----------------------------------------------------------------------------------------------------------*/
/*Used to calculate shadow matrix*/
/*----------------------------------------------------------------------------------------------------------*/

#include "turtle.h"
#include "util_math.h"


/*----------------------------------------------------------------------------------------------------------*/

void tridiag(short a, long r, long c, long nx, DOUBLEVECTOR *diag_inf, DOUBLEVECTOR *diag, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e)

{

long j;
double bet;
DOUBLEVECTOR *gam;

gam=new_doublevector(nx);

if(diag->co[1]==0.0){
	printf("type=%ld r=%ld c=%ld\n",a,r,c);
	t_error("Error 1 in tridiag");
}

bet=diag->co[1];
e->co[1]=b->co[1]/bet;

//Decomposition and forward substitution
for(j=2;j<=nx;j++){
	gam->co[j]=diag_sup->co[j-1]/bet;
	bet=diag->co[j]-diag_inf->co[j-1]*gam->co[j];
	if(bet==0.0){
		printf("type=%ld r=%ld c=%ld\n",a,r,c);
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

double norm(DOUBLEVECTOR *V){

	long l;
	double N=0.0;

	for(l=1;l<=V->nh;l++){
		N+=pow(V->co[l],2.0);
	}

	N=pow(N,0.5);
	return(N);

}

/*----------------------------------------------------------------------------------------------------------*/

