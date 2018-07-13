
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

extern long number_novalue;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void i_lrc_cont(Matrix<double> *LC, long ***i, LONGMATRIX *lrc, long nl,
                long nr, long nc)
{

  long cont=0;
  long l, r, c;

  for (r=1; r<=nr; r++)
    {
      for (c=1; c<=nc; c++)
        {
          if ((long)(*LC)(r,c)!=number_novalue)
            {
              for (l=0; l<=nl; l++)
                {
                  cont++;
                  i[l][r][c]=cont;
                  lrc->co[cont][1]=l;
                  lrc->co[cont][2]=r;
                  lrc->co[cont][3]=c;
                }
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void j_rc_cont(Matrix<double> *LC, long **j, Matrix<long> *rc, long nr, long nc)
{

  long cont=0;
  long r, c;
  for (r=1; r<=nr; r++)
    {
      for (c=1; c<=nc; c++)
        {
          if ((long)(*LC)(r,c)!=number_novalue)
            {
              cont++;
              j[r][c]=cont;
              (*rc)(cont,1)=r;
              (*rc)(cont,2)=c;
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void lch3_cont(long **ch3, LONGMATRIX *lch, long nl, long nch)
{

  long cont=0;
  long l, ch;

  for (ch=1; ch<=nch; ch++)
    {
      for (l=0; l<=nl; l++)
        {
          cont++;
          ch3[l][ch]=cont;
          lch->co[cont][1]=l;
          lch->co[cont][2]=ch;
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void cont_nonzero_values_matrix2(long *tot, long *totdiag, CHANNEL *cnet,
                                 Matrix<double> *LC, LONGMATRIX *lrc, long ***i, long n, long nch, long nl)
{

  long j, jj, l, r=0, c=0;
  long cnt=0;
  long N, NCH;

  N = n*(nl+1);
  NCH = nch*(nl+1);

  for (j=1; j<=N+NCH; j++)
    {

      if (j<=N)
        {
          l=lrc->co[j][1];
          r=lrc->co[j][2];
          c=lrc->co[j][3];
        }
      else
        {
          jj=j-N;
          l=cnet->lch->co[jj][1];
        }

      //the cell itself
      //cnt ++;

      //the cell below
      if (l<nl) cnt ++;

      if (j<=N)
        {
          if ((long)(*LC)(r-1,c)!=number_novalue)
            {
              if (i[l][r-1][c]>j) cnt ++;
            }

          if ((long)(*LC)(r+1,c)!=number_novalue)
            {
              if (i[l][r+1][c]>j) cnt ++;
            }

          if ((long)(*LC)(r,c-1)!=number_novalue)
            {
              if (i[l][r][c-1]>j) cnt ++;
            }

          if ((long)(*LC)(r,c+1)!=number_novalue)
            {
              if (i[l][r][c+1]>j) cnt ++;
            }

          if (l>0 && cnet->ch->co[r][c]>0) cnt++;
        }
    }

  *tot = cnt;
  *totdiag = N+NCH;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void cont_nonzero_values_matrix3(Vector<long> *Lp, Vector<long> *Li,
                                 CHANNEL *cnet, Matrix<double> *LC, LONGMATRIX *lrc, long ***i, long n, long nch,
                                 long nl)
{

  //Ai = line index
  //Ap = number of values for each row
  long j,jj,l,r=0,c=0;
  long cnt = 0;
  long N, NCH;

  N = n*(nl+1);
  NCH = nch*(nl+1);

  for (j=1; j<=N+NCH; j++)
    {

      if (j<=N)
        {
          l=lrc->co[j][1];
          r=lrc->co[j][2];
          c=lrc->co[j][3];
        }
      else
        {
          jj=j-N;
          l=cnet->lch->co[jj][1];
        }

      //the cell itself
      //cnt++;
      //Li->co[cnt] = j;

      //the cell below
      if (l<nl)
        {
          cnt++;
          Li->co[cnt] = j+1;
        }

      if (j<=N)
        {
          if ((long)(*LC)(r-1,c)!=number_novalue)
            {
              if (i[l][r-1][c]>j)
                {
                  cnt++;
                  Li->co[cnt] = i[l][r-1][c];
                }
            }

          if ((long)(*LC)(r+1,c)!=number_novalue)
            {
              if (i[l][r+1][c]>j)
                {
                  cnt++;
                  Li->co[cnt] = i[l][r+1][c];
                }
            }

          if ((long)(*LC)(r,c-1)!=number_novalue)
            {
              if (i[l][r][c-1]>j)
                {
                  cnt++;
                  Li->co[cnt] = i[l][r][c-1];
                }
            }

          if ((long)(*LC)(r,c+1)!=number_novalue)
            {
              if (i[l][r][c+1]>j)
                {
                  cnt++;
                  Li->co[cnt] = i[l][r][c+1];
                }
            }

          if (l>0 && cnet->ch->co[r][c]>0)
            {
              cnt++;
              Li->co[cnt] = N + cnet->ch3[l][cnet->ch->co[r][c]];
            }

        }

      (*Lp)(j) = cnt;
    }


}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

