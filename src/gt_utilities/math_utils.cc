/* STATEMENT:

 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)

 Copyright (c), 2016 - GEOtop Foundation

 This file is part of GEOtop 2.1

 GEOtop 2.1  is a free software and is distributed under GNU General Public
 License v. 3.0 <http://www.gnu.org/licenses/> WITHOUT ANY WARRANTY; without
 even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE

 GEOtop 2.1  is distributed as a free software in the hope to create and support
 a community of developers and users that constructively interact. If you just
 use the code, please give feedback to the authors and the community. Any way
 you use the model, may be the most trivial one, is significantly helpful for
 the future development of the GEOtop model. Any feedback will be highly
 appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#include "math_utils.h"
#include "read_command_line.h"

/*===============functions copied from utilities.c ===================*/

void stop_execution(void)
{
  char ch;

  printf("\nPRESS RETURN TO CONTINUE\n");
  int ret = scanf("%c", &ch);
  // silence warning
  (void)ret;
}

double Fmin(double a, double b)
{
  double min = a;

  if (b < a) min = b;

  return (min);
}

long Fminlong(long a, long b)
{
  long min = a;

  if (b < a) min = b;

  return (min);
}

double Fmax(double a, double b)
{
  double max = a;

  if (b > a) max = b;

  return (max);
}

long Fmaxlong(long a, long b)
{
  long max = a;

  if (b > a) max = b;

  return (max);
}

/*===============================functions copied from
 * util_math.c========================================*/

short tridiag2(long nbeg,
               long nend,
               const GeoVector<double> &ld,
               const GeoVector<double> &d,
               const GeoVector<double> &ud,
               const GeoVector<double> &b,
               GeoVector<double> &e)
// solve A(ld,d,ud) * e + b = 0
{
  long j;
  double bet;
#ifdef VERYVERBOSE
  size_t lIndex = 0;
#endif
  GeoVector<double> gam;

  gam.resize(nend + 1);

  bet = d[nbeg];
  if (bet == 0.0) { return 1; }

  e[nbeg] = -b[nbeg] / bet;

  // Decomposition and forward substitution
  for (j = nbeg + 1; j <= nend; j++)
    {
      gam[j] = ud[j - 1] / bet;
      bet = d[j] - ld[j - 1] * gam[j];
      if (bet == 0.0) { return 1; }
      e[j] = (-b[j] - ld[j - 1] * e[j - 1]) / bet;
    }

  // Backsubstitution
  for (j = (nend - 1); j >= nbeg; j--)
    {
      e[j] -= gam[j + 1] * e[j + 1];
    }

#ifdef VERY_VERBOSE
  printf("DEBUG_PRINT: nbeg(%ld),nend(%ld)\n", nbeg, nend);
  printf("DEBUG_PRINT: ld(");
  for (lIndex = 0; lIndex < ld.size(); lIndex++)
    {
      printf("(%lu:%.12g)", lIndex, ld[lIndex]);
    }
  printf(")\n");
  printf("DEBUG_PRINT: d(");
  for (lIndex = 0; lIndex < d.size(); lIndex++)
    {
      printf("(%lu:%.12g)", lIndex, d[lIndex]);
    }
  printf(")\n");
  printf("DEBUG_PRINT: b(");
  for (lIndex = 0; lIndex < b.size(); lIndex++)
    {
      printf("(%lu:%.12g)", lIndex, b[lIndex]);
    }
  printf(")\n");
  printf("DEBUG_PRINT: e(");
  for (lIndex = 0; lIndex < e.size(); lIndex++)
    {
      printf("(%lu:%.12g)", lIndex, e[lIndex]);
    }
  printf(")\n");
#endif
  return 0;
}

/*----------------------------------------------------------------------------------------------------------*/

double norm_inf(const GeoVector<double> &V, long nbeg, long nend)
{
  long l;
  double N = 0.0;

  for (l = nbeg; l < nend; l++)
    {
      if (fabs(V[l]) > N) N = fabs(V[l]);
    }

  return (N);
}

//-----------------------------------------

double norm_1(const GeoVector<double> &V, long nbeg, long nend)
{
  long l;
  double N = 0.0;
  for (l = nbeg; l < nend; l++)
    {
      N += fabs(V[l]);
    }

  return (N);
}

/*----------------------------------------------------------------------------------------------------------*/

double norm_2(const GeoVector<double> &V, long nbeg, long nend)
{
  long l;
  double N = 0.0;

  for (l = nbeg; l < nend; l++)
    {
      N += V[l]*V[l];
    }
  N = sqrt(N);

  return (N);
}

void Cramer_rule(double A,
                 double B,
                 double C,
                 double D,
                 double E,
                 double F,
                 double *x,
                 double *y)
{
  /*
  Ax + By = C
  Dx + Ey = F
  x = (CE - FB) / (AE - DB)
  y = (AF - CD) / (AE - DB)
  */
  *x = (C * E - F * B) / (A * E - D * B);
  *y = (A * F - C * D) / (A * E - D * B);
}

/*----------------------------------------------------------------------------------------------------------*/

double minimize_merit_function(double res0,
                               double lambda1,
                               double res1,
                               double lambda2,
                               double res2)
{
  double lambda;
  double a, b, c;  // interpolating polynomial: ax2 + bx + c

  // calculate three-point quadratic polynomial interpolating the merit function
  c = res0;
  Cramer_rule(pow(lambda1, 2.0), lambda1, res1 - res0, pow(lambda2, 2.0),
              lambda2, res2 - res0, &a, &b);

  // minimize ax^2+bx+c
  if (a > 0)
    {
      lambda = -b / (2 * a);
      if (lambda < lambda1 * GTConst::thmin)
        {
          lambda = lambda1 * GTConst::thmin;
        }
      else if (lambda > lambda1 * GTConst::thmax)
        {
          lambda = lambda1 * GTConst::thmax;
        }
    }
  else
    {
      if (a * lambda1 * GTConst::thmin * lambda1 * GTConst::thmin +
          b * lambda1 * GTConst::thmin + c <
          a * lambda1 * GTConst::thmax * lambda1 * GTConst::thmax +
          b * lambda1 * GTConst::thmax + c)
        {
          lambda = lambda1 * GTConst::thmin;
        }
      else
        {
          lambda = lambda1 * GTConst::thmax;
        }
    }

  return (lambda);
}

//
// Recursive auxiliary function for adaptiveSimpsons() function below
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsonsAux(double (*f)(double),
                           double a,
                           double b,
                           double epsilon,
                           double S,
                           double fa,
                           double fb,
                           double fc,
                           int bottom)
{
  double c = (a + b) / 2, h = b - a;
  double d = (a + c) / 2, e = (c + b) / 2;
  double fd = f(d), fe = f(e);
  double Sleft = (h / 12) * (fa + 4 * fd + fc);
  double Sright = (h / 12) * (fc + 4 * fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15 * epsilon) return S2 + (S2 - S) / 15;
  return adaptiveSimpsonsAux(f, a, c, epsilon / 2, Sleft, fa, fc, fd,
                             bottom - 1) +
         adaptiveSimpsonsAux(f, c, b, epsilon / 2, Sright, fc, fb, fe,
                             bottom - 1);
}

//
// Adaptive Simpson's Rule
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsons(double (*f)(double),  // ptr to function
                        double a,
                        double b,                 // interval [a,b]
                        double epsilon,           // error tolerance
                        int maxRecursionDepth)    // recursion cap
{
  double c = (a + b) / 2, h = b - a;
  double fa = f(a), fb = f(b), fc = f(c);
  double S = (h / 6) * (fa + 4 * fc + fb);
  return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc,
                             maxRecursionDepth);
}

/*----------------------------------------------------------------------------------------------------------*/

//
// Recursive auxiliary function for adaptiveSimpsons() function below
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsonsAux2(double (*f)(double x, void *p),
                            void *arg,
                            double a,
                            double b,
                            double epsilon,
                            double S,
                            double fa,
                            double fb,
                            double fc,
                            int bottom)
{
  double c = (a + b) / 2, h = b - a;
  double d = (a + c) / 2, e = (c + b) / 2;
  double fd = f(d, arg), fe = f(e, arg);
  double Sleft = (h / 12) * (fa + 4 * fd + fc);
  double Sright = (h / 12) * (fc + 4 * fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 0 || fabs(S2 - S) <= 15 * epsilon) return S2 + (S2 - S) / 15;
  return adaptiveSimpsonsAux2(f, arg, a, c, epsilon / 2, Sleft, fa, fc, fd,
                              bottom - 1) +
         adaptiveSimpsonsAux2(f, arg, c, b, epsilon / 2, Sright, fc, fb, fe,
                              bottom - 1);
}

//
// Adaptive Simpson's Rule
// from http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
//
double adaptiveSimpsons2(double (*f)(double x, void *p),
                         void *arg,  // ptr to function
                         double a,
                         double b,                 // interval [a,b]
                         double epsilon,           // error tolerance
                         int maxRecursionDepth)    // recursion cap
{
  double c = (a + b) / 2, h = b - a;
  double fa = f(a, arg), fb = f(b, arg), fc = f(c, arg);
  double S = (h / 6) * (fa + 4 * fc + fb);
  return adaptiveSimpsonsAux2(f, arg, a, b, epsilon, S, fa, fb, fc,
                              maxRecursionDepth);
}

/////////////////////////////////////////////////////////////////////////////
// These routines coming from sparse_matrix_library.. S.C. 07.12.2013
/////////////////////////////////////////////////////////////////////////////

void product_matrix_using_lower_part_by_vector_plus_vector(
  double k,
  GeoVector<double> &out,
  const GeoVector<double> &y,
  const GeoVector<double> &x,
  const GeoVector<long> &Li,
  const GeoVector<long> &Lp,
  GeoVector<double> &Lx)
{
  /// @brief calculates k*(y + Ax), where k is coefficient, y and x vectors, and
  /// A a SPD matrix defined with its lower diagonal part

  size_t i;

  product_using_only_strict_lower_diagonal_part(out, x, Li, Lp, Lx);

  for (i = 1; i < x.size(); i++)
    {
#ifdef VERY_VERBOSE
      printf(
        "product_matrix_using_lower_part_by_vector_plus_vector-> i:%ld k:%e "
        "out:%e y:%e out:%e\n",
        i, k, out[i], y[i], k * (out[i] + y[i]));
#endif
      out[i] = k * (out[i] + y[i]);
    }
}

//======

void product_using_only_strict_lower_diagonal_part(GeoVector<double> &product,
                                                   const GeoVector<double> &x,
                                                   const GeoVector<long> &Li,
                                                   const GeoVector<long> &Lp,
                                                   GeoVector<double> &Lx)
{
  long c, r;

  for (size_t i = 1; i < x.size(); i++)
    {
      product[i] = 0.0;
    }

  c = 1;
  for (size_t i = 1; i < Li.size(); i++)
    {
      r = Li[i];

      if (r > c)
        {
          product[c] += Lx[i] * (x[r] - x[c]);
          product[r] += Lx[i] * (x[c] - x[r]);

#ifdef VERY_VERBOSE
          printf("c:%ld -> %e i:%ld Lx:%e r:%ld %e c:%ld %e -> %e \n", c,
                 product[c], i, Lx[i], r, x[r], c, x[c], Lx[i] * (x[r] - x[c]));
          printf("r:%ld -> %e i:%ld Lx:%e c:%ld %e r:%ld %e -> %e \n", r,
                 product[r], i, Lx[i], c, x[c], r, x[r], Lx[i] * (x[c] - x[r]));
#endif

        }
      else if (r < c)
        {
#ifdef VERBOSE
          printf(
            "product_using_only_strict_lower_diagonal_part r:%ld c:%ld i:%ld "
            "Ap[c]:%ld tot:%ld %ld %ld\n",
            r, c, i, Lp[c], Li.size(), Lp[x.size() - 1], x.size());
#endif
          t_error(
            "product_using_only_strict_lower_diagonal_part:matrix is not L, see "
            "function: " +
            std::string(__FUNCTION__));
        }

      if (i < Li.size() && (long)Li.size() > 0)
        {
          while (c < (long)Lp.size() && (long)i >= Lp[c])
            c++;
        }
    }
}

//====================================

long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(
  double tol_rel,
  double tol_min,
  double tol_max,
  GeoVector<double> &x,
  const GeoVector<double> &b,
  GeoVector<double> &y,
  const GeoVector<long> &Li,
  const GeoVector<long> &Lp,
  const GeoVector<double> &Lx)
{
  /// @brief solve sistem (A+Iy)*x = B, find x
  // A M-matrix described by its lower diagonal part

  GeoVector<double> r0, r, p, v, s, t, diag, udiag, yy, z;
  double rho, rho1, alpha, omeg, beta, norm_r0;
  size_t i = 0, maxiter;
  short sux;
  size_t j;

  r0.resize(x.size());
  r.resize(x.size());
  p.resize(x.size());
  v.resize(x.size());
  s.resize(x.size());
  t.resize(x.size());
  diag.resize(x.size());
  udiag.resize(x.size() - 1);
  yy.resize(x.size());
  z.resize(x.size());

  get_diag_strict_lower_matrix_plus_identity_by_vector(diag, udiag, y, Li, Lp,
                                                       Lx);

  // product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(r, x,
  // y, Li, Lp, Lx);

  /// printf("==== >> x->nl=%ld x->nh=%ld ",x->nl,x->nh);
  // stop_execution();

  for (j = GTConst::nl; j < x.size(); j++)
    {
      r[j] = b[j];
      r0[j] = r[j];
      p[j] = 0.;
      v[j] = 0.;
    }

  norm_r0 = norm_2(r0, GTConst::nl, r0.size());
  if (norm_r0 != norm_r0)
    {
#ifdef VERBOSE
      printf(" BiCGSTAB_strict norm_r0: %f\n", norm_r0);
#endif
      t_error("Fatal Error! NAN in a variable. See failing report.");
    }
#ifdef VERBOSE
  printf(" BiCGSTAB_strict norm_r0: %f\n", norm_r0);
#endif

  rho = 1.;
  alpha = 1.;
  omeg = 1.;

  maxiter = x.size() / 100;
  if (maxiter < 100) maxiter = 100;

  while (i < maxiter && norm_2(r, GTConst::nl, r.size()) >
         Fmax(tol_min, Fmin(tol_max, tol_rel * norm_r0)))
    {
      rho1 = product(r0, r);

      beta = (rho1 / rho) * (alpha / omeg);

      rho = rho1;

      for (j = GTConst::nl; j < x.size(); j++)
        {
          p[j] = r[j] + beta * (p[j] - omeg * v[j]);
        }

      sux = tridiag(0, 0, 0, x.size(), udiag, diag, udiag, p, yy);

      if (sux == 0) return (-1);

      product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(
        v, yy, y, Li, Lp, Lx);

      alpha = rho / product(r0, v);

      for (j = GTConst::nl; j < x.size(); j++)
        {
          s[j] = r[j] - alpha * v[j];
        }

      if (norm_2(s, GTConst::nl, s.size()) > 1.E-10)
        {
          sux = tridiag(0, 0, 0, x.size(), udiag, diag, udiag, s, z);
          if (sux == 0) return (-1);

          product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(
            t, z, y, Li, Lp, Lx);

          omeg = product(t, s) / product(t, t);

          for (j = GTConst::nl; j < x.size(); j++)
            {
              x[j] += (alpha * yy[j] + omeg * z[j]);
              r[j] = s[j] - omeg * t[j];
            }

        }
      else
        {
          for (j = GTConst::nl; j < x.size(); j++)
            {
              x[j] += alpha * yy[j];
              r[j] = s[j];
            }
        }

      i++;
    }

  return i;
}

//===============================

short tridiag(short a,
              long r,
              long c,
              long nx,
              const GeoVector<double> &diag_inf,
              const GeoVector<double> &diag,
              const GeoVector<double> &diag_sup,
              const GeoVector<double> &b,
              GeoVector<double> &e)
{
  long j;
  double bet;

  GeoVector<double> gam;

  gam.resize(nx);

  if (diag[1] == 0.0)
    {
      printf("type=%d r=%ld c=%ld\n", a, r, c);
      t_error("Error 1 in tridiag");
    }

  bet = diag[1];
  e[1] = b[1] / bet;

  for (j = 2; j < nx; j++)
    {
      gam[j] = diag_sup[j - 1] / bet;
      bet = diag[j] - diag_inf[j - 1] * gam[j];
      if (bet == 0.0)
        {
          printf("type=%d r=%ld c=%ld\n", a, r, c);
          printf("l=%ld diag(l)=%f diag_inf(l-1)=%f diag_sup(l-1)=%f\n", j, diag[j],
                 diag_inf[j - 1], diag_sup[j - 1]);
          printf("Error 2 in tridiag\n");
          return 0;
        }

      e[j] = (b[j] - diag_inf[j - 1] * e[j - 1]) / bet;
    }

  // Backsubstitution
  for (j = (nx - 2); j >= 1; j--)
    {
      e[j] -= gam[j + 1] * e[j + 1];
    }

  // free_doublevector(gam);

  return 1;
}

void get_diag_strict_lower_matrix_plus_identity_by_vector(
  GeoVector<double> &diag,
  GeoVector<double> &udiag,
  const GeoVector<double> &y,
  const GeoVector<long> &Li,
  const GeoVector<long> &Lp,
  const GeoVector<double> &Lx)
{
  long r, c;
  size_t i;

  // find diagonal and upper diagonal of matrix A+Iy, where A is described by
  // its strict lower diagonal part

  for (i = 1; i < diag.size(); i++)
    {
      diag[i] = y[i];
      if (i < diag.size() - 1) udiag[i] = 0.;
    }

  c = 1;
  for (i = 1; i < Li.size(); i++)
    {
      r = Li[i];

      diag[c] -= Lx[i];
      diag[r] -= Lx[i];

      if (r == c + 1) udiag[c] = Lx[i];
      if (i < Li.size() - 1 && (long)i > 0)
        {
          while (c < (long)Lp.size() && (long)i >= Lp[c])
            c++;
        }
    }
}

void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(
  GeoVector<double> &product,
  const GeoVector<double> &x,
  const GeoVector<double> &y,
  const GeoVector<long> &Li,
  const GeoVector<long> &Lp,
  const GeoVector<double> &Lx)
{
  long r, c;
  size_t i;

  // calculate (A+Iy)*x, A described by its strict lower diagonal part

  for (i = 1; i < x.size(); i++)
    {
      product[i] = y[i] * x[i];
    }

  c = 1;

  for (i = 1; i < Li.size(); i++)
    {
      r = Li[i];

      if (r > c)
        {
          product[c] += Lx[i] * (x[r] - x[c]);
          product[r] += Lx[i] * (x[c] - x[r]);
        }
      else if (r < c)
        {
          t_error(
            "product_using_only_strict_lower_diagonal_part_plus_identity_by_vector:"
            " matrix is not L, see function: " +
            std::string(__FUNCTION__));
        }

      if (i < Li.size() && (long)i > 0)
        {
          while (c < (long)Lp.size() && (long)i >= Lp[c])
            c++;
        }
    }
}

double product(const GeoVector<double> &a, const GeoVector<double> &b)
{
  double p = 0.;

  size_t i;

  for (i = 1; i < a.size(); i++)
    {
      p += a[i] * b[i];
    }

  return (p);
}
