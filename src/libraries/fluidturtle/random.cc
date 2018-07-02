#include "turtle.h"
#include "t_random.h"

/* For the Numerical Recipes Random Generator  ran1 */
#define PI 3.141592654
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243200
#define IA3 4561
#define IC3 51349
#define NR_END 1
/* For the Numerical Recipes Random Generator  ran2 */

#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/IM1)
#define IMM1 (IM1-1)
#define IAA1 40014
#define IAA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2  3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX  (1.0 - EPS)


/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/



