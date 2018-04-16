#define NRANSI
#define THRESH 0
#define ITOL 3
#define TOL 0.00001
#define ITMAX 1000
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define EPS 1.0e-14
#define TINY 1.0E-20

/**

Name: convlv

Synopsis: void convlv(FLOATVECTOR *data,unsigned long n,FLOATVECTOR *respns,unsigned long m,int isign,
            double delta_time,FLOATVECTOR *ans);


Description: Convolves or deconvolves a real data set data [1..n] with a responce function respns[1..n]
The responce function must be stored in wrap-around order in the first element of respns, where m is an odd
integer <=n. Wrap-around order means that the first half of the array respns contains the impulse responce
function at positive times, while the second half of the array contains the impulceresponce function at negative
times, counting down from the highest element respns[m]. On input isign is +1 for convolution, -1 for deconvolution.
The answer is returned in the first n componenets of ans. However, ans must be supplied in the calling program with
dimension [1..2*n], for consistency with twofft. n MUST be an integer power of two.

Authors & Date: Angelo Zacchia, Marco Pegoretti, 1998

Inputs: data is a complex array of lenght n

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

Notes: It is a modified version of the routine in Numerical Recipes, second edition

*/

void convlv(FLOATVECTOR *data,unsigned long n,FLOATVECTOR *respns,
            unsigned long m,int isign,
            double delta_time,FLOATVECTOR *ans);


/**

Name: four1

Synopsis: void four1(FLOATVECTOR *data,unsigned long nn,int isign)

Description: Replaces data[1..2*nn] by its discrete Fourier transform,if isign is input as 1;
or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input
as -1, data is a complex array of lenght nn or, equivalently, a real array of lenght 2*nn:
nn MUST be an integer power of 2(this is not checked for).

Authors & Date: Angelo Zacchia, Marco Pegoretti, 1998

Inputs: data is a complex array of lenght n

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

Notes: It is a modified version of the routine in Numerical Recipes, second edition


*/
void four1(FLOATVECTOR *data,unsigned long nn,int isign);


/**

Name: ludcmp

Synopsis: int ludcmp(SHORTVECTOR *indx, DOUBLEMATRIX *var)

Description: Convolves or deconvolves a real data set data [1..n] with a responce function respns[1..n]
The responce function must be stored in wrap-around order in the first element of respns, where m is an odd
integer <=n. Wrap-around order means that the first half of the array respns contains the impulse responce
function at positive times, while the second half of the array contains the impulceresponce function at negative
times, counting down from the highest element respns[m]. On input isign is +1 for convolution, -1 for deconvolution.
The answer is returned in the first n componenets of ans. However, ans must be supplied in the calling program with
dimension [1..2*n], for consistency with twofft. n MUST be an integer power of two.

Authors & Date: Angelo Zacchia, Marco Pegoretti, 1998

Inputs:

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

Notes: It is a modified version of the routine in Numerical Recipes, second edition

Examples:

References: Pegoretti, Marco, Geomodel, implementazione di un modello scalabile di deflusso
e bilancio idrologico di bacino, 1997

*/

void ludcmp(SHORTVECTOR *indx, DOUBLEMATRIX *var);

/**

Name: lubksb

Synopsis: int lubksb(DOUBLEMATRIX *var, SHORTVECTOR *indx,DOUBLEVECTOR *gam)

Description: Convolves or deconvolves a real data set data [1..n] with a responce function respns[1..n]
The responce function must be stored in wrap-around order in the first element of respns, where m is an odd
integer <=n. Wrap-around order means that the first half of the array respns contains the impulse responce
function at positive times, while the second half of the array contains the impulceresponce function at negative
times, counting down from the highest element respns[m]. On input isign is +1 for convolution, -1 for deconvolution.
The answer is returned in the first n componenets of ans. However, ans must be supplied in the calling program with
dimension [1..2*n], for consistency with twofft. n MUST be an integer power of two.

Authors & Date: Angelo Zacchia, Marco Pegoretti, 1998

Inputs: data is a complex array of lenght n

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

Notes: It is a modified version of the routine in Numerical Recipes, second edition

Examples:

References: Pegoretti, Marco, Geomodel, implementazione di un modello scalabile di deflusso
e bilancio idrologico di bacino, 1997


*/

void lubksb(DOUBLEMATRIX *var, SHORTVECTOR *indx,DOUBLEVECTOR *gam);

/**

Name: realft

Synopsis: void realft(FLOATVECTOR *data,unsigned long n,int isign);

Description: Calculates the Fourier transform of a set of n real-value data points. Replaces
this data (which is stored in array data[1..n]) by the positive frequency half of its complex Fourier
transform. The real-valued first and last components of the complex transform are returned as elements
data[1] data[2], respectively, n must be a power of 2. This routine also calculates the inverse transform
of a complex data array if it is the transform of real data (result in this case must be multiplied by 2/n).

Authors & Date: modified from NR by Marco Pegoretti and Angelo Zacchia, 1998

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h



*/

void realft(FLOATVECTOR *data,unsigned long n,int isign);

/**

Name: twofft

Synopsis: voidtwofft(FLOATVECTOR *data1,FLOATVECTOR *data2,FLOATVECTOR *fft1,FLOATVECTOR *fft2,unsigned long n);

Description: Given two real input array data1[1..2n] and data2[1..2n] this routine calls four1
and return two complex output arrays, fft1[1..2n] and fft2[1..2n], each of complex lenght
n, which contain the discrete Fourier transforms of the respective data arrays. n MUST
be an integer power of 2

Authors & Date: For the Numerical Recipes,1998.

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

*/
void twofft(FLOATVECTOR *data1,FLOATVECTOR *data2,FLOATVECTOR *fft1,
            FLOATVECTOR *fft2,unsigned long n);

/**

Name: ris_sistema


Synopsis:

void ris_sistema (double d[], double ds[], double di[], double b[], double x[], int n);

Version:0.9

Description: This solves the specified linear system.
   Questa funzione innesca la soluzione del sistema lineare.
   Ha solo lo scopo di ricevere gli argomenti necessari dal programma
   principale, allocare e deallocare le strutture necessarie e di
   restituire il risultato.

   Definizione dei parametri:
   - THRESH: gli elementi della matrice che in valore assoluto sono minori di
             THRESH vengono trascurati;
   - ITOL: varia da 1 a 4 e determina quale criterio di convergenza adottare
           per la soluzione del sistema (vedi Numerical Recipes in C pag. 86);
   - TOL: tolleranza ammessa nella ricerca della soluzione;
   - ITMAX: numero massimo di iterazioni.

Inputs:
   In input riceve gli elementi delle diagonali principale d[], superiore ds[],
   inferiore di[], il vettore dei termini noti b[], una soluzione di primo
   tentativo x[] e la dimensione del sistema, n.

Return:
  Vengono restituiti gli elementi del vettore soluzione del sistema x[].

   Tutti vettori sopra menzionati devono essere variabili globali, sono quindi
   dichiarati all'interno del programma principale.


 Authors & Date: Angelo Zacchia, June 1998

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h


*/

//commentata per l'errore (dovuto ad una errata definizione di asolve):
//error: conflicting types for `asolve'
//void ris_sistema (double d[], double ds[], double di[], double b[], double x[], int n);


/**


Name: sprsin

Synopsis:

void sprsin(double **,int,float,long,double *,long *);

Description:
  Funzione che ordina una matrice sparsa (in questo caso tridiagonale)
  alla maniera di Num. Rec.
  Questa funzione converte una matrice memorizzata nel modo convenzionale
  in un vettore sa[] che contiene solo i valori non nulli della matrice
  e in un vettore ija[] che permette di individuare la posizione originale
  degli elementi di sa[].

Inputs:
  - **a, un puntatore agli elementi della matrice originale;
  - n, dimensione della matrice;
  - thresh, gli elementi della matrice minori di thresh non vengono
            letti;
  - nmax, la lunghezza dei vettori sa[] e ija[].

Return:


 Authors & Date: Angelo Zacchia, June 1998

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

 References: DA NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press).



 */
void sprsin(double **,int,float,long,double *,long *);

/**

Name: linbcg


Description:
   Questa funzione consente di risolvere di risolvere un sistema lineare del tipo
   A x = b con il metodo iterativo del gradiente coniugato.

Inputs:
   - n: dimensione del sistema;
   - sa[] e ija[]: vettori generati dalla funzione sprssin() che memorizzano la
     matrice;
   - b[]: elementi del vettore dei termini noti;
   - x[]: elementi del vettore soluzione (in ingresso questo vettore deve contenere
     una soluzione di primo tentativo);
   - itol, tol, itmax: parametri gli definiti sopra.

   Oltre alla soluzione la funzione calcola anche il numero di iterazione effetuate
   ( iter ) e l'errore commesso ( err ).


 Authors & Date: Angelo Zacchia, June 1998

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

    References:  NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press). pag 86-88

     Bugs & limitations:

 */

// commentata per l'errore (dovuto ad una errata definizione di asolve):
//error: conflicting types for `asolve'
//void linbcg(long, double *,double *, int, double, int, int *,double *,
//             double *, long *);


/**

Name:       vett_mat

Version:

Synopsis:

Description:

  Questa funzione converte tre vettori in una matrice quadrata tridiagonale.


 Inputs:
    Bisogna passare alla funzione i puntatori agli elementi dei vettori
   - d elementi delle diagonale principale;
   - ds elementi della diagonale superiore;
   - di elementi della diagonale inferiore;
   inoltre bisogna passare n, che e' la dimensione della matrice che
   si vuole generare.


 Return:
    La funzione restituisce il puntatore alla matrice.


  Authors & Date: Angelo Zacchia, June 1998

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

 */

DOUBLEMATRIX *vett_mat (double *d,double *ds,double *di,int n);


/**

Name:     snrm

Version:

Synopsis:   double snrm(long n, double sx[], int itol);


Description:

Questa funzione calcola la norma di un vettore con la modalita' specificata
   dal parametro itol

 Authors & Date: Angelo Zacchia, June 1998

FILE: LIBRARIES/LINEARALGEBRA/linearalgebra.c, LIBRARIES/LINEARALGEBRA/linearalgebra.h

References: DA NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press). pag 88

Bugs & limitations:



*/

double snrm(long n, double sx[], int itol);


/**
Name:  atimes

Version:

Synopsis:
void atimes(long n, double x[], double r[], int itrnsp,double sa[], long ija[]);

Description:

References: DA NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press). pag 88

Notes:

*/

void atimes(long n, double x[], double r[], int itrnsp,double sa[],
            long ija[]);

/**

Name:  asolve

Version:

Synopsis:   void asolve(long n, double b[], double x[], int itrnsp,double sa[], long ija[])


Description:


References: DA NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press). pag 89


*/


void asolve(long n, double b[], double x[],double sa[]);

/**

Name:  dsprsax

Version:

Synopsis:   void dsprsax(double sa[], long ija[], double x[], double b[], long n)


Description:
Questa funzione moltiplica una matrice, memoririzzata alla maniera di N.R.,
   per un vettore x[]. Il risultato e' un vettore b[].


Inputs:

Return:


References: DA NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press). pag 89


*/

void dsprsax(double sa[], long ija[], double x[], double b[], long n);


/**
Name:  dsprstx


Synopsis:   void dsprsax(double sa[], long ija[], double x[], double b[], long n)


Description:
Questa funzione moltiplica la trasposta di una matrice, memoririzzata alla
   maniera di N.R., per un vettore x[]. Il risultato e' un vettore b[].


Inputs:

Return:

Needs:

Related Routines:

See Also:

References: DA NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press). pag 80



*/


void dsprstx(double sa[], long ija[], double x[], double b[], long n);

/**
Name:  integration


Synopsis:   float integration(float (*func)(float), float a, float b, int n);

Description:
Questa funzione calcola l'integrale con il metodo dei trapezi


Inputs:

Return:

Needs:

Related Routines:

See Also:

References: DA NUMERICAL RECIPES IN C. (Second Edition - Cambridge Univ. Press).
pag 137
*/
float integration ( float ( *fun )( float ), float, float, int);
