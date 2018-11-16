
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 1.225-15

 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#include "turtle.h"
#include "util_math.h"
#include "constants.h"

/*----------------------------------------------------------------------------------------------------------*/

short tridiag(short a, long r, long c, long nx, Vector<double> *diag_inf,
              Vector<double> *diag, Vector<double> *diag_sup, Vector<double> *b, Vector<double> *e)

{

  long j;
  double bet;
  std::unique_ptr<Vector<double>> gam {new Vector<double>{nx}};

  /*for (j=1; j<=nx; j++) {
    printf("d[%ld]=%e\n",j,diag->co[j]);
  }
  for (j=1; j<nx; j++) {
    printf("d[%ld]=%e %e\n",j,diag_inf->co[j],diag_sup->co[j]);
  }*/
  if (diag->co[1]==0.0)
    {
      printf("type=%d r=%ld c=%ld\n",a,r,c);
      t_error("Error 1 in tridiag");
    }

  bet=diag->co[1];
  e->co[1]=b->co[1]/bet;

  //Decomposition and forward substitution
  for (j=2; j<=nx; j++)
    {
      gam->co[j]=diag_sup->co[j-1]/bet;
      bet=diag->co[j]-diag_inf->co[j-1]*gam->co[j];
      if (bet==0.0)
        {
          /*printf("type=%d r=%ld c=%ld\n",a,r,c);
          printf("l=%ld diag(l)=%f diag_inf(l-1)=%f diag_sup(l-1)=%f\n",j,diag->co[j],diag_inf->co[j-1],diag_sup->co[j-1]);
          printf("Error 2 in tridiag\n");*/
          return 0;
        }
      e->co[j]=(b->co[j]-diag_inf->co[j-1]*e->co[j-1])/bet;
    }

  //Backsubstitution
  for (j=(nx-1); j>=1; j--)
    {
      e->co[j]-=gam->co[j+1]*e->co[j+1];
    }

  return 1;

}

/*----------------------------------------------------------------------------------------------------------*/


short tridiag2(short  /*a*/, long  /*r*/, long  /*c*/, long nbeg, long nend,
               Vector<double> *ld, Vector<double> *d, Vector<double> *ud, Vector<double> *b,
               Vector<double> *e)

//solve A(ld,d,ud) * e + b = 0

{
  long j;
  double bet;
  std::unique_ptr<Vector<double>> gam{new Vector<double>{nend}};

  bet = d->co[nbeg];
  if (bet == 0.0)
    {
      return 1;
      //printf("type=%d r=%ld c=%ld\n",a,r,c);
      //t_error("Error 1 in tridiag");
    }
  e->co[nbeg] = -b->co[nbeg]/bet;

  //Decomposition and forward substitution
  for (j=nbeg+1; j<=nend; j++)
    {
      gam->co[j] = ud->co[j-1]/bet;
      bet = d->co[j]-ld->co[j-1]*gam->co[j];
      if (bet == 0.0)
        {
          return 1;
          //printf("type=%d r=%ld c=%ld\n",a,r,c);
          //printf("l=%ld d(l)=%f ld(l-1)=%f ud(l-1)=%f\n",j,d->co[j],ld->co[j-1],ud->co[j-1]);
          //t_error("Error 2 in tridiag");
        }
      e->co[j] = (-b->co[j]-ld->co[j-1]*e->co[j-1])/bet;
    }

  //Backsubstitution
  for (j=(nend-1); j>=nbeg; j--)
    {
      e->co[j] -= gam->co[j+1]*e->co[j+1];
    }

  return 0;

}

/*----------------------------------------------------------------------------------------------------------*/

double norm_inf(Vector<double> *V, long nbeg, long nend)
{

  long l;
  double N=0.0;

  for (l=nbeg; l<=nend; l++)
    {
      if (fabs(V->co[l])> N) N = fabs(V->co[l]);
    }

  return (N);

}

/*----------------------------------------------------------------------------------------------------------*/

double norm_2(Vector<double> *V, long nbeg, long nend)
{

  long l;
  double N=0.0;

  for (l=nbeg; l<=nend; l++)
    {
      N+= (V->co[l])*(V->co[l]);
    }
  N=sqrt(N);

  return (N);

}

/*----------------------------------------------------------------------------------------------------------*/

double norm_1(Vector<double> *V, long nbeg, long nend)
{

  long l;
  double N=0.0;

  for (l=nbeg; l<=nend; l++)
    {
      N += fabs(V->co[l]);
    }

  return (N);

}
/*----------------------------------------------------------------------------------------------------------*/


void Cramer_rule(double A, double B, double C, double D, double E, double F,
                 double *x, double *y)
{

  /*
  Ax + By = C
  Dx + Ey = F

  x = (CE - FB) / (AE - DB)
  y = (AF - CD) / (AE - DB)
  */

  *x = (C*E - F*B) / (A*E - D*B);
  *y = (A*F - C*D) / (A*E - D*B);

}

/*----------------------------------------------------------------------------------------------------------*/


double minimize_merit_function(double res0, double lambda1, double res1,
                               double lambda2, double res2)
{

  double lambda;
  double a, b, c;     //interpolating polynomial: ax2 + bx + c

  //calculate three-point quadratic polynomial interpolating the merit function
  c = res0;
  Cramer_rule(pow(lambda1, 2.0), lambda1, res1-res0, pow(lambda2, 2.0), lambda2,
              res2-res0, &a, &b);

  //minimize ax^2+bx+c
  if (a>0)
    {
      lambda = -b/(2*a);
      if (lambda < lambda1*thmin)
        {
          lambda = lambda1*thmin;
        }
      else if (lambda > lambda1*thmax)
        {
          lambda = lambda1*thmax;
        }
    }
  else
    {
      if (a * lambda1*thmin*lambda1*thmin + b * lambda1*thmin + c < a *
          lambda1*thmax*lambda1*thmax + b * lambda1*thmax + c)
        {
          lambda = lambda1*thmin;
        }
      else
        {
          lambda = lambda1*thmax;
        }
    }

  return (lambda);
}

/*----------------------------------------------------------------------------------------------------------*/

double product(Vector<double> *a, Vector<double> *b)
{

  double p=0.;
  long i,n=a->nh;

  for (i=1; i<=n; i++)
    {
      p += a->co[i] * b->co[i];
    }

  return (p);

}

/*----------------------------------------------------------------------------------------------------------*/

//
// Recursive auxiliary function for adaptiveSimpsons() function below
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsonsAux(double (*f)(double), double a, double b,
                           double epsilon,
                           double S, double fa, double fb, double fc, int bottom)
{
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d), fe = f(e);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
    return S2 + (S2 - S)/15;
  return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
         adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

//
// Adaptive Simpson's Rule
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsons(double (*f)(double),   // ptr to function
                        double a, double b,  // interval [a,b]
                        double epsilon,  // error tolerance
                        int maxRecursionDepth)     // recursion cap
{
  double c = (a + b)/2, h = b - a;
  double fa = f(a), fb = f(b), fc = f(c);
  double S = (h/6)*(fa + 4*fc + fb);
  return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc,
                             maxRecursionDepth);
}

/*----------------------------------------------------------------------------------------------------------*/

//
// Recursive auxiliary function for adaptiveSimpsons() function below
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsonsAux2(double (*f)(double x, void *p), void *arg,
                            double a, double b, double epsilon,
                            double S, double fa, double fb, double fc, int bottom)
{
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
  double fd = f(d, arg), fe = f(e, arg);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
    return S2 + (S2 - S)/15;
  return adaptiveSimpsonsAux2(f, arg, a, c, epsilon/2, Sleft,  fa, fc, fd,
                              bottom-1) +
         adaptiveSimpsonsAux2(f, arg, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

//
// Adaptive Simpson's Rule
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsons2(double (*f)(double x, void *p),
                         void *arg,   // ptr to function
                         double a, double b,  // interval [a,b]
                         double epsilon,  // error tolerance
                         int maxRecursionDepth)     // recursion cap
{
  double c = (a + b)/2, h = b - a;
  double fa = f(a, arg), fb = f(b, arg), fc = f(c, arg);
  double S = (h/6)*(fa + 4*fc + fb);
  return adaptiveSimpsonsAux2(f, arg, a, b, epsilon, S, fa, fb, fc,
                              maxRecursionDepth);
}

/*----------------------------------------------------------------------------------------------------------*/

