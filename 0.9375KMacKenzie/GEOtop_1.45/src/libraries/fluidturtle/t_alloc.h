float *vector(long nl,long nh);
int *ivector(long nl, long nh);
long  *lvector(long nl, long nh) ;
short  *svector(long nl, long nh) ;
double *dvector(long nl, long nh) ;
char *cvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl,long nch);
int **imatrix(long nrl, long nrh, long ncl,long nch);
short **smatrix(long nrl, long nrh, long ncl,long nch);
long **lmatrix(long nrl, long nrh, long ncl,long nch);
double **dmatrix(long nrl, long nrh, long ncl,long nch);
char **charmatrix(long nrl, long nrh, long ncl,long nch);

/*-------------------------------------------------------------------------------
You will find for each method (function, routine) acting on a type a replica acting
on the other types. To indicate collectively the replicas we will use "*". Thus
in that context,   * substitutes:
shortvector, intvector, longvector, floatvector, doublevector, charvector,
shortmatrix, intmatrix, longmatrix, floatmatrix, doublematrix, charmatrix,
shortbin, intbin, longbin,doublebin ,stringbin.
---------------------------------------------------------------------------------*/


/**

Name:  new_

Synopsis:

SHORTVECTOR *new_shortvector(long);

INTVECTOR *new_intvector(long);

FLOATVECTOR *new_floatvector(long);

LONGVECTOR *new_longvector(long);

DOUBLEVECTOR *new_doublevector(long);

CHARVECTOR *new_charvector(long);


SHORTMATRIX *new_shortmatrix( long,long);

INTMATRIX *new_intmatrix( long,long);

FLOATMATRIX *new_floatmatrix( long,long);

LONGMATRIX *new_longmatrix( long,long);

DOUBLEMATRIX *new_doublematrix( long,long);


INTBIN *new_intbin(LONGVECTOR *);

STRINGBIN *new_stringbin(LONGVECTOR *);

SHORTBIN *new_shortbin(LONGVECTOR *);

LONGBIN *new_longbin(LONGVECTOR *);

DOUBLEBIN *new_doublebin(LONGVECTOR *);

Version: 1.0

Description:
Allocate vectors, matrixes and bind. Strings length is not known apriori.
Thus, the allocation program is one with the reading routine one finds in
t_io.h.


Authors & Date: Paolo D'Odorico, Riccardo Rigon, 1996-1997

FILE: LIBRARIES/BASICS/t_alloc.h, LIBRARIES/BASICS/alloc.c

Inputs: the number of elements in the  vector, the number of rows and columns matrixes,
 a LONGVECTOR for bins.

Return: a pointer to the allocated structure

See Also: _vector, free_

Keywords: Dynamic memory allocation, pointers.

Examples: 1.example.c, 2.example.c

*/

SHORTVECTOR *new_shortvector(long);

INTVECTOR *new_intvector(long);

FLOATVECTOR *new_floatvector(long);

LONGVECTOR *new_longvector(long);

DOUBLEVECTOR *new_doublevector(long);
DOUBLEVECTOR *new_doublevector0(long);

CHARVECTOR *new_charvector(long);


SHORTMATRIX *new_shortmatrix( long,long);

INTMATRIX *new_intmatrix( long,long);

FLOATMATRIX *new_floatmatrix( long,long);

LONGMATRIX *new_longmatrix( long,long);

DOUBLEMATRIX *new_doublematrix( long,long);

DOUBLEMATRIX *new_doublematrix0(long ,long );

INTBIN *new_intbin(LONGVECTOR *);

STRINGBIN *new_stringbin(LONGVECTOR *);

SHORTBIN *new_shortbin(LONGVECTOR *);

LONGBIN *new_longbin(LONGVECTOR *);

DOUBLEBIN *new_doublebin(LONGVECTOR *);





/**

Name: free_

Description:
This section includes the deallocation routines
for the vector, matrixes and bins


Version: 1.0

Synopsis:
void free_shortvector( SHORTVECTOR *);
void free_intvector( INTVECTOR *);
void free_longvector( LONGVECTOR *);
void free_floatvector( FLOATVECTOR *);
void free_doublevector( DOUBLEVECTOR *);
void free_charvector( CHARVECTOR *);

void free_shortmatrix( SHORTMATRIX *);
void free_intmatrix( INTMATRIX *);
void free_longmatrix( LONGMATRIX *);
void free_floatmatrix( FLOATMATRIX *);
void free_doublematrix( DOUBLEMATRIX *);

void free_intbin( INTBIN *);
void free_stringbin( STRINGBIN *);
void free_shortbin( SHORTBIN *);
void free_longbin(  LONGBIN *);
void free_doublebin( DOUBLEBIN *);

Authors  & date: Paolo D'Odorico, Riccardo Rigon 1996-1997

FILE: LIBRARIES/BASICS/t_alloc.h, LIBRARIES/BASICS/alloc.c

Inputs: the name of the structure to be deallocated

Return: void

Related Routines: free, malloc

Keywords: Dynamic memory allocation, Pointers.

Examples: 1.example.c, 2.example.c


*/

void free_shortvector( SHORTVECTOR *);
void free_intvector( INTVECTOR *);
void free_longvector( LONGVECTOR *);
void free_floatvector( FLOATVECTOR *);
void free_doublevector( DOUBLEVECTOR *);
void free_charvector( CHARVECTOR *);

void free_shortmatrix( SHORTMATRIX *);
void free_intmatrix( INTMATRIX *);
void free_longmatrix( LONGMATRIX *);
void free_floatmatrix( FLOATMATRIX *);
void free_doublematrix( DOUBLEMATRIX *);

void free_intbin( INTBIN *);
void free_stringbin( STRINGBIN *);
void free_shortbin( SHORTBIN *);
void free_longbin(  LONGBIN *);
void free_doublebin( DOUBLEBIN *);

void free_svector(short* v, long nl);
void free_ivector(int* v, long nl);
void free_vector(float * v, long nl);
void free_lvector(long * v, long nl);
void free_dvector(double * v, long nl);
void free_cvector(char * v, long nl);
void free_smatrix(short **m,long nrl,long ncl);
void free_imatrix(int **m,long nrl,long ncl);
void free_matrix(float **m,long nrl,long ncl);
void free_lmatrix(long **m,long nrl,long ncl);
void free_dmatrix(double **m,long nrl,long ncl);
void free_charmatrix(char **m,long nrl,long ncl);



/**

Name:  free_header

Version: 1.0

Synopsis:  void free_header(HEADER )

General information:
This section includes the deallocation routines
for blocks and header. Types HEADER are allocated by read_header.

Authors & Date: Riccardo Rigon 1996-1997

FILE: LIBRARIES/BASICS/t_alloc.h, LIBRARIES/BASICS/alloc.c

Inputs: the name of the structure to be deallocated

Return: void

Examples: 1.example.c, 2.example.c

Keywords: Dynamic memory allocation, Pointers.

*/


void free_header(HEADER );





double ***d3tensor( long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
DOUBLETENSOR *new_doubletensor(long nrh,long nch,long ndh);
DOUBLETENSOR *new_doubletensor0(long ndh,long nrh,long nch);

void free_d3tensor(double ***t, long nrl, long ncl, long ndl);
void free_doubletensor( DOUBLETENSOR *m);

double **dmatrix(long nrl, long nrh, long ncl,long nch);
double *dvector(long nl, long nh);