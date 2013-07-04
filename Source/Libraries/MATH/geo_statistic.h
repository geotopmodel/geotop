
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
void ordi_kriging2(DOUBLEMATRIX *weights, DOUBLEVECTOR *E, DOUBLEVECTOR *N, DOUBLEMATRIX *Z0, T_INIT *UV, double int_scale, double variance);


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
void variogramma2(DOUBLEVECTOR *E, DOUBLEVECTOR *N, DOUBLEMATRIX *Z0, DOUBLEMATRIX *variogr, DOUBLEVECTOR *U, double int_scale, double variance);




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
double gamma2(double r, double scala_integr, double varianza);


//int ludcmp(SHORTVECTOR *indx, DOUBLEMATRIX *var);



//int lubksb(DOUBLEMATRIX *var, SHORTVECTOR *indx,DOUBLEVECTOR *gam);

