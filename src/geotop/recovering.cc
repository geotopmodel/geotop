
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

#include "struct.geotop.h"
#include "rw_maps.h"
#include "tabs.h"

extern T_INIT *UV;
extern long number_novalue;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map(short old, long n, char *name, Matrix<double> *assign,
                          PAR *par, Matrix<double> *Zdistr)
{

  long r, c;
  std::unique_ptr<Matrix<double>> M;
  char *temp, *temp2;

  temp = namefile_i_we2(name, n);

  if (old == 1)
    {
      temp2 = join_strings(temp, ".old");
      free(temp);
      temp = assign_string(temp2);
      free(temp2);
    }

  M.reset(read_map(1, temp, Zdistr, UV, (double)number_novalue));
  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          (*assign)(r,c) = (*M)(r,c);
        }
    }

  free(temp);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map_vector(short old, long n, char *name,
                                 Vector<double> *assign, Matrix<long> *rc, PAR *par, Matrix<double> *Zdistr)
{

  long i, r, c;
  std::unique_ptr<Matrix<double>> M;
  char *temp, *temp2;

  temp = namefile_i_we2(name, n);

  if (old == 1)
    {
      temp2 = join_strings(temp, ".old");
      free(temp);
      temp = assign_string(temp2);
      free(temp2);
    }

  M.reset(read_map(1, temp, Zdistr, UV, (double)number_novalue));
  for (i=1; i<=rc->nrh; i++)
    {
      r = (*rc)(i,1);
      c = (*rc)(i,2);
      assign->co[i] = (*M)(r,c);
    }

  free(temp);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_map_long(short old, long n, char *name,
                               Matrix<long> *assign, PAR *par, Matrix<double> *Zdistr)
{

  long r, c;
  std::unique_ptr<Matrix<double>> M;
  char *temp, *temp2;

  temp = namefile_i_we2(name, n);

  if (old == 1)
    {
      temp2 = join_strings(temp, ".old");
      free(temp);
      temp = assign_string(temp2);
      free(temp2);
    }

  M.reset(read_map(1, temp, Zdistr, UV, (double)number_novalue));
  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          (*assign)(r,c) = (long)(*M)(r,c);
        }
    }
    free(temp);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_tensor(short old, long n, char *name,
                             DOUBLETENSOR *assign, PAR *par, Matrix<double> *Zdistr)
{

  long r, c, l;
  //long i;
  std::unique_ptr<Matrix<double>> M;
  char *temp1, *temp2, *temp3;

  for (l=assign->ndl; l<=assign->ndh; l++)
    {

      temp1 = namefile_i_we2(name, n);
      temp2 = namefile_i_we(temp1, l);

      if (old == 1)
        {
          temp3 = join_strings(temp2, ".old");
          free(temp2);
          temp2 = assign_string(temp3);
          free(temp3);
        }

      M.reset(read_map(1, temp2, Zdistr, UV, (double)number_novalue));

      for (r=1; r<=M->nrh; r++)
        {
          for (c=1; c<=M->nch; c++)
            {
              assign->co[l][r][c] = (*M)(r,c);
            }
        }
        free(temp2);
      free(temp1);

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_tensor_vector(short old, long n, char *name, Matrix<double> *assign,
                                    Matrix<long> *rc, PAR *par, Matrix<double> *Zdistr)
{

  long r, c, i, l;
  std::unique_ptr<Matrix<double>> M;
  char *temp1, *temp2, *temp3;

  for (l=assign->nrl; l<=assign->nrh; l++)
    {

      temp1 = namefile_i_we2(name, n);
      temp2 = namefile_i_we(temp1, l);

      if (old == 1)
        {
          temp3 = join_strings(temp2, ".old");
          free(temp2);
          temp2 = assign_string(temp3);
          free(temp3);
        }

      M.reset(read_map(1, temp2, Zdistr, UV, (double)number_novalue));
      for (i=1; i<=rc->nrh; i++)
        {
          r = (*rc)(i,1);
          c = (*rc)(i,2);
          (*assign)(l,i) = (*M)(r,c);
        }
        free(temp2);
      free(temp1);

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assign_recovered_tensor_channel(short old, long n, char *name,
                                     Matrix<double> *assign, Vector<long> *r, Vector<long> *c, Matrix<double> *Zdistr)
{

  long ch, l;
  std::unique_ptr<Matrix<double>> M;
  char *temp1, *temp2, *temp3;

  for (l=assign->nrl; l<=assign->nrh; l++)
    {

      temp1 = namefile_i_we2(name, n);
      temp2 = namefile_i_we(temp1, l);

      if (old == 1)
        {
          temp3 = join_strings(temp2, ".old");
          free(temp2);
          temp2 = assign_string(temp3);
          free(temp3);
        }

      M.reset(read_map(1, temp2, Zdistr, UV, (double)number_novalue));

      for (ch=1; ch<=r->nh; ch++)
        {
          if (r->co[ch] > 0)
            (*assign)(l,ch) = (*M)(r->co[ch],c->co[ch]);
        }
        free(temp2);
      free(temp1);

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void recover_run_averages(short old, Matrix<double> *A, char *name,
                          Matrix<double> *LC, Matrix<long> *rc, PAR *par, long n)
{

  std::unique_ptr<Matrix<double>> M;
  long j, l;

  M.reset(new Matrix<double>{n, par->total_pixel});
  assign_recovered_tensor_vector(old, par->recover, name, M.get(), rc, par, LC);
  for (j=1; j<=par->total_pixel; j++)
    {
      if ((*par->jplot)(j) > 0)
        {
          for (l=1; l<=n; l++)
            {
              (*A)((*par->jplot)(j),l) = (*M)(l,j);
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void print_run_averages_for_recover(Matrix<double> *A, char *name,
                                    long **j_cont, PAR *par, long n, long nr, long nc)
{

  std::unique_ptr<Matrix<double>> M;
  long j, l;

  M.reset(new Matrix<double>{n, par->total_pixel});
  *M = (double)number_novalue;
  for (j=1; j<=par->total_pixel; j++)
    {
      if ((*par->jplot)(j) > 0)
        {
          for (l=1; l<=n; l++)
            {
              (*M)(l,j) = (*A)((*par->jplot)(j),l);
            }
        }
    }
  for (l=1; l<=n; l++)
    {
      rename_tensorseries(1, l, 0, name);
      write_tensorseries_vector(1, l, 0, name, 0, par->format_out, M.get(), UV,
                                number_novalue, j_cont, nr, nc);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
