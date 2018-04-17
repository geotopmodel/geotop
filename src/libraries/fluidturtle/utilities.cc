#include "turtle.h"
#include "t_utilities.h"
#include "math.h"


/*--------------------------------------------------------------------------------------*/


void meter( long  index, long rows,  short frequence,const char *message,
            const char *separator)

{
  short m;
  if (frequence>rows)
    {
      printf("%s %ld/%ld%s", message, index,rows,separator);
    }
  else
    {
      m=ceil (rows/frequence);
      if ( (index  % m )==0)
        {
          printf("%s %ld/%ld%s", message, index,rows,separator);
        }
    }
}

/*--------------------------------------------------------------------------*/

double Fmin(double a, double b)
{

  double min=a;

  if (b<a) min=b;

  return (min);

}

long Fminlong(long a, long b)
{

  long min=a;

  if (b<a) min=b;

  return (min);

}

/*--------------------------------------------------------------------------*/

double Fmax(double a, double b)
{

  double max=a;

  if (b>a) max=b;

  return (max);

}

long Fmaxlong(long a, long b)
{

  long max=a;

  if (b>a) max=b;

  return (max);

}


/*--------------------------------------------------------------------------*/
