#ifndef TURTLE_H
#define TURTLE_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#define element co //to inhertit old nomunclature for FluidTurle

#define isDynamic 1         /*This number will be used to mark a dynamic allocates
                              quantity */
#define MAXBUFFERSIZE 50000 /* Many functions uses a buffer to allocate data. A
                             safe habit is to limit its size. */

#define MAXLABELSIZE 1024    /* The same as above */

#define BUFFERINCREMENT 256 /*Buffers memory is increases when necessary of this
							 storage capability till MAXBUFFERSIZE*/

#define LABELINCREMENT 64   /*The same as above */

#define NR_END    1         /* Required by the Numerical Recipes' allocation routines
							*/

#define FREE_ARG  char*     /* The same as above */

#define OK 1

#define NOK -1

#define WOK 0

#define MAX_KEYWORD_LENGTH 10

#define PRINT 1
#define NOPRINT 0         /*With the above PRINT is used to set a printint option
                           in some functions */

#define NL 1              /* Numerical Recipes allocation routines allow to have
							 arbitrary  subscripts for vector and matrixes. The
							 fluid turtle library restrict this freedom by setting
						     their  lower value to NL  */
						     
						     
#define TEST 1

#define NOTEST 0

#define SQRT2  1.414213562373095

/**-------------------------------------------------------------------
VECTORS: essentially the same types as in Numerical Recipes, except
that char and long are not unsigned as there. We have vector for short, 
int,long , char, float, double.  
---------------------------------------------------------------------*/

typedef struct {
	short isdynamic;          /* set to 1 when dynamically allocated */
	const char * name;        /* the name of the data structure      */
	long nl,nh;               /* the lower bound, nl, and the upper, nh */
	float *co;           /* the real stuff                      */	
} FLOATVECTOR;

typedef struct {
	short isdynamic;
	const char * name; 
	long nl,nh;
	int *co;
} INTVECTOR;


typedef struct {
	short isdynamic;
	const char * name; 
	long nl,nh;
	long  *co;	
} LONGVECTOR;


typedef struct {
	short isdynamic;
	const char * name; 
	long nl,nh;
	short *co;
} SHORTVECTOR;


typedef struct {
	short isdynamic;
	const char * name; 
	long nl,nh;
	double *co;
} DOUBLEVECTOR;

typedef struct {
	short isdynamic;
	const char * name; 
	long nl,nh;
	char *co;
} CHARVECTOR;                   


/**-------------------------------------------------------------------
MATRIXES: The same types as for vectors except for char:
SHORTMATRIX, INTMATRIX, LONGMATRIX, FLOATMATRIX, DOUBLEMATRIX
---------------------------------------------------------------------*/

typedef struct {
	short isdynamic;         /* see FLOATVECTOR */
	const char * name; 
	long nrl,nrh,ncl,nch;    /* lower and upper bound for rows: nrl, nrh;
	                            lower and upper bounds for columns: ncl, nch  */
	short **co;
} SHORTMATRIX;

typedef struct {
	short isdynamic;
	const char * name; 
	long nrl,nrh,ncl,nch;
	int **co;
} INTMATRIX;


typedef struct {
	short isdynamic;
	const char * name; 
	long nrl,nrh,ncl,nch;
	float **co;
} FLOATMATRIX;


typedef struct {
	short isdynamic;
	const char * name; 
	long nrl,nrh,ncl,nch;
	double **co;
} DOUBLEMATRIX;


typedef struct {
	short isdynamic;
	const char * name; 
	long nrl,nrh,ncl,nch;
	long **co;
} LONGMATRIX;


/**-------------------------------------------------------------------
BINS: A BIN is a collection of vectors of different lengths. 
We have BINS for short, long and double, a BIN of char vectors is 
here called STRINGBIN. Thus, we have:
SHORTBIN, INTBIN, LONGBIN, DOUBLEBIN, STRING
---------------------------------------------------------------------*/

typedef struct {
	short isdynamic;
	const char * name; 
	LONGVECTOR  * index;      /* The index of the list: see  the
	                             example file for more informations  */
	short **co;
} SHORTBIN;


typedef struct {
	short isdynamic;
	const char * name; 
	LONGVECTOR  *index;
	int **co;
} INTBIN;

typedef struct {
	short isdynamic;
	const char * name; 	
	LONGVECTOR  *index;
	long **co;
} LONGBIN;

typedef struct {
	short isdynamic;
	const char * name; 
	LONGVECTOR  * index;
	float **co;
} FLOATBIN;

typedef struct db{
	short isdynamic;
	const char * name; 	
	LONGVECTOR  * index;
	double **co;
	struct db *next;
} DOUBLEBIN;


typedef struct st{
	short isdynamic;
	const char * name; 
	LONGVECTOR  *index;
	char **co;
	struct st *next;
} STRINGBIN;
/*  Tensor3D */
typedef struct {
	short isdynamic;
	const char * name; 
	long nrl,nrh,ncl,nch,ndl,ndh;
	double ***co;
	
} DOUBLETENSOR;


/**-------------------------------------------------------------------------
Linear linked lists are defined in t_lists.h.  Turtle arrays are vectors written
in a file between enclosing curling braces and separated by commas: i.e.
{1,2,3} is an array of 3 integers. They are stored in memory as vectors of the
appropriate type. Thus, they do no requires a special type to be defined.
Other two special types defined in the fluid turtle library are the HEADER type
and the t_keyword type below. The HEADER contains the description of a 
vector, matrix, array or BIN. See the file turtle.dat for more details.
The structure t_keyword contains the keywords used by the turtle libraries.
--------------------------------------------------------------------------*/


typedef struct{

	long number;                 /* Its position in file */
    short   gender;               /* It can be ascii=1, binary=2 .... */
	short   type;                /* type of the vector to be stored:
	                                 char=1,short=2, int=3, long=4, float=5
	                                 double=6, string=7             */
	short  category;            /* array=1, vector=1, matrix=2,list=3 ,tensor3d=4*/
	long dimensions[4];         /* The length for a vector, the number of
	                                           rows and columns for a matrix, rows,colums,depth for a tensor3D     */
	char *name;                  /* Whatever name you choose or NULL */           
    
} HEADER;

typedef struct {
const char *gender[3];
const char *type[8];
const char *category[6];
const char *symbol[3];
const char *separator[3];
const char *delimiter[3];

}t_keywords;

/*
t_keywords T_KEYWORDS={{"2","ascii","binary"},
					   {"7","char","short","int","long","float","double","string"},
					   {"4","array","vector","matrix","list"},
					   {"2","->","<-"},
					   {"2"," ",","},
					   {"2","{","}"}};							  
*/							   
typedef struct { /*header of maps in fluid turtle format*/
DOUBLEVECTOR *U;  /*dx,dy*/
DOUBLEVECTOR *V;  /*sign of novalue,novalue*/
} T_INIT;

/* #include "nr_util.h" */
#include "t_alloc.h"
#include "t_io.h"
#include "t_list.h"

/* #include "t_random.h" */
/*
#include "t_datamanipulation.h"
#include "t_probability.h"
*/
void t_error(char *error_text);



#endif
