
/* MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file geo_statistic.09375.c

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
#include "t_datamanipulation.h"
#include "geo_statistic.09375.h"
#include "linearalgebra.h"

/*void ordi_kriging(DOUBLEMATRIX *pesi,DOUBLEMATRIX *coord,DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,
		DOUBLEVECTOR *V, double *scala_integr1, double *varianza1);

void variogramma(DOUBLEMATRIX *coord,DOUBLEMATRIX *Z0,DOUBLEMATRIX *variogr,DOUBLEVECTOR *U,
		double scala_integr, double varianza);

double gamma1(double r, double scala_integr, double varianza);

int ludcmp(SHORTVECTOR *indx, DOUBLEMATRIX *var);
int lubksb(DOUBLEMATRIX *var, SHORTVECTOR *indx,DOUBLEVECTOR *gam);
*/

/**

Name: ordi_kriging

Synopsis: void ordi_kriging(DOUBLEMATRIX *pesi,DOUBLEMATRIX *coord,DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,
		DOUBLEVECTOR *V, double *scala_integr1, double *varianza1)

General information: Calcola la distribuzione delle pioggia tramite algoritmo di ordinary kriging.
	I pesi vengono calcolati tramite la risoluzione delle equazioni con una decomposizione LU
	(cfr. 2.3 in Numerical Recepies in C)

Authors & Date: Riccardo Rigon, Marco Pegoretti, Giacomo Bertoldi 1998

	Inputs: 1) coord: coord matrice con coordinate stazioni
			2) Z0 matrice con elevazioni
			3) U vettore con dimensioni pixels
			4) V vettore con header elevazioni
			5) scala_integr1: integral spatial scale
			6) varianza1: spatial variance
	Outputs:1) pesi: krigingweights	DOUBLEMATRIX [row*col,n_stazioni]

Needs:  t_io.c, t_alloc.c, t_error.c, LinearAlgebra.c

Example: Pegoretti, Marco, Geomodel, implementazione di un modello scalabile di deflusso e bilancio idrologico di bacino,
               Tesi di Laurea, Relatore R.Rigon, Universita' degli studi di Trento, A.A. 1997-98.

See Also: turtle.dat file

References: Pegoretti, Marco, Geomodel, implementazione di un modello scalabile di deflusso e bilancio idrologico di bacino,
                  Tesi di Laurea, Relatore R.Rigon, Universita' degli studi di Trento, A.A. 1997-98.
Bugs:

*/


void ordi_kriging2(DOUBLEMATRIX *weights, DOUBLEVECTOR *E, DOUBLEVECTOR *N, DOUBLEMATRIX *Z0, T_INIT *UV, double int_scale, double variance)

{

long meno,i,j,k,row,col,n,ii,ncol;
short mode;
double r,x,y,dx,dy,coordx,coordy;
DOUBLEMATRIX *var,*variogr;
DOUBLEVECTOR *gam;
SHORTVECTOR *indx;

/* coordinate vertice basso sx matrice */
/* controllare la convenzione FluidTurtle!
in r.out.turtle mi sembra essere: {dx,dy,S,W}
*/

coordy=UV->U->co[3];
coordx=UV->U->co[4];

mode=0;
meno=0;
row=Z0->nrh;
col=Z0->nch;
ncol=row*col;
dx=UV->U->co[1];
dy=UV->U->co[2];
n=E->nh; /* numero stazioni */


if(n>1){ /* nn>1: piu stazioni: faccio il kriging */
	var=new_doublematrix(n+1,n+1);
	gam=new_doublevector(n+1);
	indx=new_shortvector(n+1);
	for(i=1;i<=n+1;i++){
   		indx->co[i]=0;
	}
	/* Applico la funzione del kriging*/
	/* valuto la coordinata di ogni stazione rispetto alla coordinata dell'angolo
   		sinistro in basso della mia matrice*/
	/* Calcolo la funzione var[i][j] definita dalla distanza tra la stazione[i]
   		rispetto alla stazione j*/
	for(i=1;i<=n;i++){
   		for(j=1;j<=n;j++){
    		x=E->co[i]-E->co[j];
    		y=N->co[i]-N->co[j];
    		r=sqrt(pow(x,(double)2)+pow(y,(double)2)); /* r: distanza */
    		/************* calcolo il variogramma con gamma1(r,int_scale,variance) in geo_statistic.c *********/
    		var->co[i][j]=gamma2(r, int_scale, variance);
    		/**************************************************************/
   		}
	}
	for(i=1;i<=n;i++){
    	var->co[i][n+1]=1;
	}
	for(j=1;j<=n;j++){
    	var->co[n+1][j]=1;
 	}
	var->co[n+1][n+1]=0;
	/*Calcolo le coordinate di ogni pixel rispetto allo stesso punto di riferimento*/
	/* variogr: doublematrix (nch*nrh, n stazioni) */
	variogr=new_doublematrix(ncol,n);
	/************* chiamo variogramma in geo_statistic.c *********/
	/* calcola la variogramma tra i pixel e le stazioni  */
	variogramma2(E, N , Z0, variogr, UV->U, int_scale, variance);
	/**************************************************************/
	if(n>1){
		for(i=1;i<=row;i++){
			for(j=1;j<=col;j++){
				/*Calcolo la distanza di ogni pixel rispetto a tutte le stazioni*/
				ii=(i-1)*col+j;
				for(k=1;k<=n;k++){
					/* gam: vettore con la distanza del pixel dalle stazioni */
					gam->co[k]=variogr->co[ii][k];
				}
				gam->co[n+1]=1;
				/*Le seguenti subroutine mi risolvono il sistema che avra' come matrice gli elementi
					var relativi alle stazioni, come vettore dei termini noti il vettore gam  relativo
					ad ogni pixel. La prima routine mi riporta la matrice var modificata per eseguire una
					sotituzione indietro, questo e' sufficiente eseguirla un unica volta */
				if(i==1 && j==1){
					/************* chiamo ludcmp in linearalgebra.c *********/
					ludcmp(indx,var);
					/**************************************************************/
				}
				/************* chiamo lubksb in linearalgebra.c *********/
				lubksb(var,indx,gam);
				/**************************************************************/
				/* metto i risultati ottenuti in gam nella matrice pesi */
				if(Z0->co[i][j]!=UV->V->co[2]){
					for(k=1;k<=n;k++){
						weights->co[ii][k]=gam->co[k];
					}
				}
			}
		}
	}else{
		for(i=1;i<=row;i++){
			for(j=1;j<=col;j++){
				/*Calcolo la distanza di ogni pixel rispetto a tutte le stazioni*/
				ii=(i-1)*col+j;
				weights->co[ii][1]=variogr->co[ii][2]/(variogr->co[ii][1]+variogr->co[ii][2]);
				weights->co[ii][2]=variogr->co[ii][1]/(variogr->co[ii][1]+variogr->co[ii][2]);
			}
		}
	}
	free_doublevector(gam);
	free_doublematrix(var);
	free_doublematrix(variogr);
	free_shortvector(indx);
}else{ /* n=1: 1 stazione: non faccio il kriging e pongo tutti i pesi=1 */
	//printf("THERE IS ONE STATION\n");
	for(i=1;i<=Z0->nrh;i++){
		for(j=1;j<=Z0->nch;j++){
			ii=(i-1)*col+j;
			if(Z0->co[i][j]!=UV->V->co[2]){
				weights->co[ii][1]=1;
			}
		}
	}
}
/* tolgo eventuali pesi < 0 che non hanno senso */
for(i=1;i<=ncol;i++){
	for(j=1;j<=n;j++){
    	if(weights->co[i][j]<=0) weights->co[i][j]=0;
	}
}
}


/**
Name: variogramma;

void variogramma(DOUBLEMATRIX *coord,DOUBLEMATRIX *Z0,DOUBLEMATRIX *variogr,DOUBLEVECTOR *U,
		double scala_integr, double varianza);


General information: variogramma e' una routine che riporta il variogramma di una matrice di elevazioni
	rispetto ad alcuni punti di cui sono note le coordinate. Questi punti potrebbero essere per esempio
	delle stazioni pluviometriche

Authors & Date: Riccardo Rigon, Paolo Verardo, 1999

Inputs: 1)Z0: matrice delle elevazioni(Z0)
        2)coord: matrice che indica le coordinate dei punti
        3)U: vettore con le dimensioni dei pixel e le coordinate X,Y dell'angolo in basso a sinistra
        4)variogr: matrice con la varianza di ogni pixel rispetto alle stazioni
		4) V vettore con header elevazioni
		5) scala_integr1: integral spatial scale
		6) varianza1: spatial variance

Needs:  t_io.c, t_alloc.c, t_error.c.

Example: Pegoretti, M., Geomodel:implementazione di un modello scalabile di deflusso e bilancio idrologico di bacino,
              Tesi di Laurea, Relatore R.Rigon, Universita' degli studi di Trento, A.A. 1997-98.

See Also: turtle.dat file

References:  Pegoretti, M., Geomodel:implementazione di un modello scalabile di deflusso e bilancio idrologico di bacino,
              Tesi di Laurea, Relatore R.Rigon, Universita' degli studi di Trento, A.A. 1997-98.


Bugs:

*/

void variogramma2(DOUBLEVECTOR *E, DOUBLEVECTOR *N, DOUBLEMATRIX *Z0, DOUBLEMATRIX *variogr, DOUBLEVECTOR *U, double int_scale, double variance)

{
long ii,i,j,k;
double X,Y,rx,ry,rz,h2;

ii=0;
for(i=1;i<=Z0->nrh;i++){
	Y=U->co[3]+(Z0->nrh-i)*U->co[2]+(U->co[2]/(double)2);
	for(j=1;j<=Z0->nch;j++){
	    ii++;
	    for(k=1;k<=E->nh;k++){
			X=U->co[4]+(j-1)*U->co[1]+(U->co[1]/(double)2);
			rx=X-E->co[k];
			ry=Y-N->co[k];
			rz=0;
			h2=sqrt(pow(rx,(double)2)+pow(ry,(double)2)+pow(rz,(double)2));
			if(h2!=0)variogr->co[ii][k]=gamma2(h2,int_scale,variance);
		}
	}
}
}


/**
Name: gamma1

double gamma1(double r, double scala_integr, double varianza);


General information: e` una funzione di correalzione esponenziale

Authors & Date: Riccardo Rigon, Paolo Verardo, Giacomo Bertoldi 1999

Inputs: 1) r distanza [m]
        2) scala_integr scala integrale [m]
        3) varianza

Needs:  t_io.c, t_alloc.c, t_error.c.

Example: Pegoretti, M., Geomodel:implementazione di un modello scalabile di deflusso e bilancio idrologico di bacino,
              Tesi di Laurea, Relatore R.Rigon, Universita' degli studi di Trento, A.A. 1997-98.

See Also: turtle.dat file

References:  Pegoretti, M., Geomodel:implementazione di un modello scalabile di deflusso e bilancio idrologico di bacino,
              Tesi di Laurea, Relatore R.Rigon, Universita' degli studi di Trento, A.A. 1997-98.
Bugs:

*/

double gamma2(double r, double scala_integr, double varianza)

{
double gamma2;

gamma2=varianza*(1.0-exp(-r/scala_integr));

return gamma2;
}

