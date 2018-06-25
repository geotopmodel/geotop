float *vector(long nl,long nh);

float **matrix(long nrl, long nrh, long ncl,long nch);

short **smatrix(long nrl, long nrh, long ncl,long nch);
long **lmatrix(long nrl, long nrh, long ncl,long nch);
double **dmatrix(long nrl, long nrh, long ncl,long nch);

/*-------------------------------------------------------------------------------
You will find for each method (function, routine) acting on a type a replica acting
on the other types. To indicate collectively the replicas we will use "*". Thus
in that context,   * substitutes:
shortmatrix, intmatrix, longmatrix, floatmatrix, doublematrix, charmatrix,
shortbin, intbin, longbin,doublebin ,stringbin.
---------------------------------------------------------------------------------*/


SHORTMATRIX *new_shortmatrix( long,long);

LONGMATRIX *new_longmatrix( long,long);

DOUBLEMATRIX *new_doublematrix( long,long);

DOUBLEMATRIX *new_doublematrix0_(long,long );

DOUBLEMATRIX *new_doublematrix_0(long,long );


void free_shortmatrix( SHORTMATRIX *);

void free_longmatrix( LONGMATRIX *);

void free_doublematrix( DOUBLEMATRIX *);

void free_smatrix(short **m,long nrl,long ncl);

void free_lmatrix(long **m,long nrl,long ncl);
void free_dmatrix(double **m,long nrl,long ncl);


double ***d3tensor( long nrl, long nrh, long ncl, long nch, long ndl,
                    long ndh);
DOUBLETENSOR *new_doubletensor(long nrh,long nch,long ndh);
DOUBLETENSOR *new_doubletensor0(long ndh,long nrh,long nch);

void free_d3tensor(double ***t, long nrl, long ncl, long ndl);
void free_doubletensor( DOUBLETENSOR *m);

double **dmatrix(long nrl, long nrh, long ncl,long nch);
