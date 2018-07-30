#include "turtle.h"
#include "t_utilities.h"
#include "math.h"


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
