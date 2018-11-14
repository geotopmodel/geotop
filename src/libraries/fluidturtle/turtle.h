#ifndef __turtle_h_
#define __turtle_h_

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cctype>

#include "vector.h"
#include "matrix.h"
#include "tensor.h"

#define isDynamic 1         /*This number will be used to mark a dynamic allocates
                              quantity */

#define MAXBUFFERSIZE 50000 /* Many functions uses a buffer to allocate data. A
                             safe habit is to limit its size. */


#define MAXLABELSIZE 1024    /* The same as above */


#define BUFFERINCREMENT 256 /*Buffers memory is increases when necessary of this
               storage capability till MAXBUFFERSIZE*/


#define LABELINCREMENT 64   /*The same as above */


#define NR_END    1         /* Required by the Numerical Recipes' allocation routines */

#define FREE_ARG  char*     /* The same as above */

#define OK 1

#define NOK -1

#define WOK 0

#define MAX_KEYWORD_LENGTH 10

#define PRINT 1

#define NOPRINT 0         /*With the above PRINT is used to set a printint option
                           in some functions */


#define NL 1   /* Numerical Recipes allocation routines allow to have
               arbitrary  subscripts for vector and matrixes. The
               fluid turtle library restrict this freedom by setting
                 their  lower value to NL  */

#define TEST 1

#define NOTEST 0

#define SQRT2  1.414213562373095


/*  Tensor3D */

typedef struct
{
    short isdynamic;

    const char *name;

    long nrl,nrh,ncl,nch,ndl,ndh;

    double ***co;

} DOUBLETENSOR;



/*

t_keywords T_KEYWORDS={{"2","ascii","binary"},

             {"7","char","short","int","long","float","double","string"},

             {"4","array","vector","matrix","list"},

             {"2","->","<-"},

             {"2"," ",","},

             {"2","{","}"}};

*/
/*header of maps in fluid turtle format*/
struct T_INIT
{
    std::unique_ptr<Vector<double>> U;  /*dx,dy*/
    std::unique_ptr<Vector<double>> V;  /*sign of novalue,novalue*/
};


/* #include "nr_util.h" */

#include "t_alloc.h"

#include "t_io.h"

/*#include "t_list.h"

#include "t_random.h" /

#include "t_datamanipulation.h"

#include "t_probability.h"

*/

void t_error(const char *error_text);

#endif
