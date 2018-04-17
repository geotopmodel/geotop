
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
#include "rw_maps.h"
#include "import_ascii.h"
#include "write_ascii.h"
#include "extensions.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "tensor3D.h"


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

SHORTMATRIX *copyshort_doublematrix(DOUBLEMATRIX *M)
{

  SHORTMATRIX *S;
  long r, c;

  S=new_shortmatrix(M->nrh,M->nch);
  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          S->co[r][c]=(short)M->co[r][c];
        }
    }

  return (S);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

LONGMATRIX *copylong_doublematrix(DOUBLEMATRIX *M)
{

  LONGMATRIX *L;
  long r, c;

  L=new_longmatrix(M->nrh,M->nch);
  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          L->co[r][c]=(long)M->co[r][c];
        }
    }

  return (L);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *copydouble_shortmatrix(SHORTMATRIX *S)
{

  DOUBLEMATRIX *M;
  long r, c;

  M=new_doublematrix(S->nrh,S->nch);
  for (r=1; r<=S->nrh; r++)
    {
      for (c=1; c<=S->nch; c++)
        {
          M->co[r][c]=(double)S->co[r][c];
        }
    }

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *copydouble_longmatrix(LONGMATRIX *L)
{

  DOUBLEMATRIX *M;
  long r, c;

  M=new_doublematrix(L->nrh,L->nch);
  for (r=1; r<=L->nrh; r++)
    {
      for (c=1; c<=L->nch; c++)
        {
          M->co[r][c]=(double)L->co[r][c];
        }
    }

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *copydoublematrix_const(double c0, DOUBLEMATRIX *Mref,
                                     double NOVALUE)
{

  DOUBLEMATRIX *M;
  long r, c;

  M=new_doublematrix(Mref->nrh,Mref->nch);
  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          if (Mref->co[r][c]==NOVALUE)
            {
              M->co[r][c]=NOVALUE;
            }
          else
            {
              M->co[r][c]=c0;
            }
        }
    }

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *multiplydoublematrix(double f, DOUBLEMATRIX *Mref,
                                   double NOVALUE)
{

  DOUBLEMATRIX *M;
  long r, c;

  M=new_doublematrix(Mref->nrh,Mref->nch);
  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          if (Mref->co[r][c]==NOVALUE)
            {
              M->co[r][c]=NOVALUE;
            }
          else
            {
              M->co[r][c]=f*Mref->co[r][c];
            }
        }
    }

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void build_doubletensor(DOUBLETENSOR *T, DOUBLEMATRIX *M, long l)
{

  long r,c;

  if (l<=0 || l>T->ndh) t_error("Invalid doubletensor construction");
  if (T->nrh!=M->nrh) t_error("Invalid doubletensor construction");
  if (T->nch!=M->nch) t_error("Invalid doubletensor construction");

  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          T->co[l][r][c]=M->co[r][c];
        }
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *extract_doublematrix(DOUBLETENSOR *T, long l)
{

  long r,c;
  DOUBLEMATRIX *M;

  if (l<=0 || l>T->ndh) t_error("Invalid doubletensor extraction");

  M=new_doublematrix(T->nrh,T->nch);

  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          M->co[r][c]=T->co[l][r][c];
        }
    }

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *extract_fromtensor(DOUBLETENSOR *T, long l)
{

  long r, c;
  DOUBLEMATRIX *M;

  if (l<1 || l>T->ndh) t_error("cannot extract a matrix from a tensor");

  M=new_doublematrix(T->nrh,T->nch);
  for (r=1; r<=T->nrh; r++)
    {
      for (c=1; c<=T->nch; c++)
        {
          M->co[r][c]=T->co[l][r][c];
        }
    }

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLETENSOR *build_frommatrix(DOUBLEMATRIX *M, long l, long lmax)
{

  long r,c;
  DOUBLETENSOR *T;

  if (lmax<1) t_error("cannot build a tensor from a matrix");
  if (l>lmax) t_error("cannot build a tensor from a matrix");

  T=new_doubletensor(lmax, M->nrh, M->nch);
  initialize_doubletensor(T, 0.0);
  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          T->co[l][r][c]=M->co[r][c];
        }
    }

  return (T);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_frommatrix(long l, DOUBLEMATRIX *M, DOUBLETENSOR *T)
{

  long r,c;

  if (l<0 || l>T->ndh) t_error("cannot write a tensor from a matrix");
  if (M->nrh!=T->nrh) t_error("cannot write a tensor from a matrix");
  if (M->nch!=T->nch) t_error("cannot write a tensor from a matrix");

  for (r=1; r<=M->nrh; r++)
    {
      for (c=1; c<=M->nch; c++)
        {
          T->co[l][r][c]=M->co[r][c];
        }
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void fmultiplydoublematrix(DOUBLEMATRIX *destination, DOUBLEMATRIX *origin,
                           double f, double novalue)
{

  long r,c;

  if (origin->nrh!=destination->nrh) t_error("cannot copy matrix");
  if (origin->nch!=destination->nch) t_error("cannot copy matrix");
  for (r=1; r<=origin->nrh; r++)
    {
      for (c=1; c<=origin->nch; c++)
        {
          if (origin->co[r][c]!=novalue)
            {
              destination->co[r][c]=f*origin->co[r][c];
            }
          else
            {
              destination->co[r][c]=novalue;
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void assignnovalue(DOUBLEMATRIX *destination, DOUBLEMATRIX *origin,
                   double novalue)
{

  long r,c;

  if (origin->nrh!=destination->nrh) t_error("cannot assign novalue");
  if (origin->nch!=destination->nch) t_error("cannot assign novalue");
  for (r=1; r<=origin->nrh; r++)
    {
      for (c=1; c<=origin->nch; c++)
        {
          if (origin->co[r][c]==novalue) destination->co[r][c]=novalue;
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_suffix(char *suffix, long i, short start)
{

  short m=0, c=0, d=0, u=0;

  if (i>=0 && i<=9)
    {
      m=0;
      c=0;
      d=0;
      u=i;
    }
  else if (i<=99)
    {
      m=0;
      c=0;
      d=(short)(i/10.0);
      u=i-10.0*d;
    }
  else if (i<=999)
    {
      m=0;
      c=(short)(i/100.0);
      d=(short)((i-100.0*c)/10.0);
      u=i-100.0*c-10*d;
    }
  else if (i<=9999)
    {
      m=(short)(i/1000.0);
      c=(short)((i-1000.0*m)/100.0);
      d=(short)((i-1000.0*m-100.0*c)/10.0);
      u=i-1000*m-100.0*c-10*d;
    }
  else
    {
      t_error("Number too high");
    }

  m+=48;
  c+=48;
  d+=48;
  u+=48;

  suffix[start]=m;
  suffix[start+1]=c;
  suffix[start+2]=d;
  suffix[start+3]=u;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

char *namefile_i(char *name, long i)
{

  char SSSS[ ]= {"SSSS"};
  char *name_out;
  char *temp;

  write_suffix(SSSS, i, 0);

  temp=join_strings(name,SSSS);
  name_out=join_strings(temp,textfile);
  free(temp);

  return (name_out);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

char *namefile_i_we(char *name, long i)
{

  char SSSS[ ]= {"LSSSS"};
  char *name_out;

  write_suffix(SSSS, i, 1);

  name_out=join_strings(name,SSSS);

  return (name_out);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

char *namefile_i_we2(char *name, long i)
{

  char SSSS[ ]= {"SSSS"};
  char *name_out;

  write_suffix(SSSS, i, 0);

  name_out=join_strings(name,SSSS);

  return (name_out);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short existing_file(char *name)
{

  //if the file exists gives 1 (fluidturtle), 2(grassascii), 3(esriascii), 0 if the file doesn't exist

  short a=0;
  FILE *f;
  //char *ft;
  char *esri,*grass;

  //ft=join_strings(name,ascii_ft);
  esri=join_strings(name,ascii_esri);
  grass=join_strings(name,ascii_grass);

  //if( (f=fopen(ft,"r"))!=NULL ){
  //  a=1;
  //  fclose(f);
  //}else
  if ( (f=fopen(grass,"r"))!=NULL )
    {
      a=2;
      fclose(f);
    }
  else if ( (f=fopen(esri,"r"))!=NULL )
    {
      a=3;
      fclose(f);
    }

  //free(ft);
  free(esri);
  free(grass);

  return (a);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short existing_file_wext(char *name, char *extension)
{

  //if the file exists gives 1, 0 if the file doesn't exist

  short a=0;
  FILE *f;
  char *temp;

  temp=join_strings(name,extension);

  if ( (f=fopen(temp,"r"))!=NULL )
    {
      a=1;
      fclose(f);
    }

  free(temp);

  return (a);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short existing_file_woext(char *name)
{

  //if the file exists gives 1, 0 if the file doesn't exist

  short a=0;
  FILE *f;

  if ( (f=fopen(name,"r"))!=NULL )
    {
      a=1;
      fclose(f);
    }

  return (a);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *read_map(short a, char *filename, DOUBLEMATRIX *Mref,
                       T_INIT *UVref, double no_value)
{

  //  a=0 non usa Mref, UVref output
  //  a=1 non esegue controllo non values, Mref e UVref input
  //  a=2 esegue controllo novalues, Mref e UVref input

  DOUBLEMATRIX *M=NULL;
  long r=0, c=0, nr=0, nc=0;
  double *header=NULL, *m=NULL;
  double Dxmap=0, Dymap=0, X0map=0, Y0map=0;

  if (a != 0 && a != 1
      && a != 2)
    t_error("Value of flag not supported in the subroutine in read_map");

  if (existing_file(filename)==1)
    {

      t_error("Fluidturtle map format not supported any more");

    }
  else
    {

      header = (double *) malloc(6*sizeof(double));

      if (existing_file(filename)==2) //grass ascii
        {
          m=read_grassascii(header, no_value, filename);
          nr=(long)header[4];
          nc=(long)header[5];
          Dxmap=(header[2]-header[3])/((long)header[5]);
          Dymap=(header[0]-header[1])/((long)header[4]);
          X0map=header[3];
          Y0map=header[1];

        }
      else if (existing_file(filename)==3)  //esri ascii
        {
          m=read_esriascii(header, no_value, filename);
          printf("%s\n",filename);
          nr=(long)header[1];
          nc=(long)header[0];
          Dxmap=header[4];
          Dymap=header[4];
          X0map=header[2];
          Y0map=header[3];

        }
      else
        {

          printf("The file %s doesn't exist\n",filename);
          t_error("Fatal error");

        }

      free(header);
      if (a==1 || a==2)
        {
          //Check header
          if (Dxmap!=UVref->U->co[2])
            {
              printf("Dx in %s file is different from Dx in DTM file! \n",filename);
              t_error("Inconsistent map");
            }
          if (Dymap!=UVref->U->co[1])
            {
              printf("Dy:%f in %s file is different from Dy:%f in DTM file! \n",Dymap,
                     filename,UVref->U->co[1]);
              t_error("Inconsistent map");
            }
          if (X0map!=UVref->U->co[4])
            {
              printf("X0 in %s file is different from X0 in DTM file! \n",filename);
              t_error("Inconsistent map");
            }
          if (Y0map!=UVref->U->co[3])
            {
              printf("Y0 in %s file is different from Y0 in DTM file! \n",filename);
              t_error("Inconsistent map");
            }
          if (nr!=Mref->nrh)
            {
              printf("Number of rows in %s file (%ld) is not consistent with DTM file (%ld)! \n",
                     filename,nr,Mref->nrh);
              t_error("Inconsistent map");
            }
          if (nc!=Mref->nch)
            {
              printf("Number of columns in %s file is not consistent with DTM file! \n",
                     filename);
              t_error("Inconsistent map");
            }

        }
      else if (a==0)
        {

          UVref->U.reset(new Vector<double>{4});
          UVref->V.reset(new Vector<double>{2});
          UVref->U->co[2]=Dxmap;
          UVref->U->co[1]=Dymap;
          UVref->U->co[4]=X0map;
          UVref->U->co[3]=Y0map;

          UVref->V->co[2]=no_value;
          if (no_value<0)
            {
              UVref->V->co[1]=-1;
            }
          else
            {
              UVref->V->co[1]=1;
            }
        }

      //assign values and check novalues
      M=new_doublematrix(nr,nc);
      for (r=1; r<=nr; r++)
        {
          for (c=1; c<=nc; c++)
            {
              M->co[r][c]=m[(r-1)*nc+c-1];
              if (a==2)
                {
                  if (M->co[r][c]==no_value && Mref->co[r][c]!=no_value)
                    {
                      printf("ERROR:: Problem reading map %s, it has novalue where the reference maps has value: %ld %ld ref:%f %f\n",
                             filename,r,c,Mref->co[r][c],M->co[r][c]);
                      t_error("Fatal Error (9)");
                    }
                  if (M->co[r][c]!=no_value && Mref->co[r][c]==no_value) M->co[r][c]=no_value;
                }
            }
        }

      free(m);
    }

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

std::unique_ptr<Vector<double>> read_map_vector(short type, char *namefile, DOUBLEMATRIX *mask,
                                                T_INIT *grid, double no_value, LONGMATRIX *rc)
{

  DOUBLEMATRIX *M;
  long i, n=rc->nrh;
  std::unique_ptr<Vector<double>> V{new Vector<double>{n}};

  M = read_map(type, namefile, mask, grid, no_value);

  for (i=1; i<=n; i++)
    {
      V->co[i] = M->co[rc->co[i][1]][rc->co[i][2]];
    }

  free_doublematrix(M);

  return V;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLEMATRIX *read_mapseries(long i, char *filename, DOUBLEMATRIX *Mref,
                             T_INIT *UVref, double no_value)
{

  char SSSS[ ]= {"LSSSS"};
  char *name;
  DOUBLEMATRIX *M;

  write_suffix(SSSS, i, 1);
  name=join_strings(filename,SSSS);
  M=read_map(2, name, Mref, UVref, no_value);
  free(name);

  return (M);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLETENSOR *read_tensor(long nl, char *filename, DOUBLEMATRIX *Mref,
                          T_INIT *UVref, double no_value)
{
  long l;
  DOUBLEMATRIX *M;
  DOUBLETENSOR *T=NULL;

  if (nl<1) t_error("The dimension of the tensor must be greater than or equal to 1");

  for (l=1; l<=nl; l++)
    {
      M=read_mapseries(l, filename, Mref, UVref, no_value);
      if (l==1)
        {
          T=build_frommatrix(M, l, nl);
        }
      else
        {
          write_frommatrix(l, M, T);
        }
      free_doublematrix(M);
    }

  return (T);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

DOUBLETENSOR *read_maptensor(long i, long lmax, char *filename,
                             DOUBLEMATRIX *Mref, T_INIT *UVref, double no_value)
{

  char SSSSLLLLL[ ]= {"SSSSLLLLL"};
  char *name;
  DOUBLETENSOR *T=NULL;
  DOUBLEMATRIX *M;
  long l;

  if (i<0 || lmax<1)
    t_error("cannot read a tensor with null or negative columns");

  write_suffix(SSSSLLLLL, i, 0);

  for (l=1; l<=lmax; l++)
    {
      write_suffix(SSSSLLLLL, l, 5);
      name=join_strings(filename,SSSSLLLLL);
      M=read_map(2, name, Mref, UVref, no_value);
      free(name);
      if (l==1)
        {
          T=build_frommatrix(M, l, lmax);
        }
      else
        {
          write_frommatrix(l, M, T);
        }
    }

  return (T);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_map(char *filename, short type, short format, DOUBLEMATRIX *M,
               T_INIT *UV, long novalue)
{

  //  type=0  floating point
  //  type=1  integer

  //  format=1 fluidturtle
  //  format=2 grassascii
  //  format=3 esriascii

  if (format==1)
    {
      t_error("The fluidturtle format is not support any more");
    }
  else if (format==2)
    {
      write_grassascii(filename, type, M, UV, novalue);
    }
  else if (format==3)
    {
      write_esriascii(filename, type, M, UV, novalue);
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_map_vector(char *filename, short type, short format,
                      Vector<double> *V, T_INIT *UV, long novalue, long **j, long nr, long nc)
{

  //  type=0  floating point
  //  type=1  integer

  //  format=1 fluidturtle
  //  format=2 grassascii
  //  format=3 esriascii

  if (format==1)
    {
      t_error("The fluidturtle format is not support any more");
    }
  else if (format==2)
    {
      write_grassascii_vector(filename, type, V, j, nr, nc, UV, novalue);
    }
  else if (format==3)
    {
      write_esriascii_vector(filename, type, V, j, nr, nc, UV, novalue);
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_mapseries(long i, char *filename, short type, short format,
                     DOUBLEMATRIX *M, T_INIT *UV, long novalue)
{

  char SSSS[ ]= {"SSSS"};
  char *name;

  write_suffix(SSSS, i, 0);
  name=join_strings(filename, SSSS);
  write_map(name, type, format, M, UV, novalue);
  free(name);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_tensorseries(short a, long l, long i, char *filename, short type,
                        short format, DOUBLETENSOR *T, T_INIT *UV, long novalue)
{

  //  a=0 non include "l" nel suffisso
  //  a=1 include "l" nel suffisso
  //  l:layer
  //  i:temporal step

  char SSSSLLLLL[ ]= {"SSSSLLLLL"};
  char SSSS[ ]= {"SSSS"};
  char *name=NULL;
  long r, c;
  DOUBLEMATRIX *M;

  if (a==0)
    {
      write_suffix(SSSS, i, 0);
      name=join_strings(filename,SSSS);
    }
  else if (a==1)
    {
      write_suffix(SSSSLLLLL, i, 0);
      write_suffix(SSSSLLLLL, l, 5);
      name=join_strings(filename,SSSSLLLLL);
    }
  else
    {
      t_error("Value not admitted");
    }


  M=new_doublematrix(T->nrh,T->nch);
  for (r=1; r<=T->nrh; r++)
    {
      for (c=1; c<=T->nch; c++)
        {
          M->co[r][c]=T->co[l][r][c];
        }
    }

  write_map(name, type, format, M, UV, novalue);

  free_doublematrix(M);
  free(name);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_tensorseries_vector(short a, long l, long i, char *filename,
                               short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J,
                               long nr, long nc)
{

  //  a=0 non include "l" nel suffisso
  //  a=1 include "l" nel suffisso
  //  l:layer
  //  i:temporal step

  char SSSSLLLLL[ ]= {"SSSSLLLLL"};
  char SSSS[ ]= {"SSSS"};
  char *name=NULL;
  long j, npoints=T->nch;
  std::unique_ptr<Vector<double>> V;

  if (a==0)
    {
      write_suffix(SSSS, i, 0);
      name=join_strings(filename,SSSS);
    }
  else if (a==1)
    {
      write_suffix(SSSSLLLLL, i, 0);
      write_suffix(SSSSLLLLL, l, 5);
      name=join_strings(filename,SSSSLLLLL);
    }
  else
    {
      t_error("Value not admitted");
    }

  V.reset(new Vector<double>{npoints});
  for (j=1; j<=npoints; j++)
    {
      V->co[j]=T->co[l][j];
    }

  write_map_vector(name, type, format, V.get(), UV, novalue, J, nr, nc);

  free(name);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void rename_tensorseries(short a, long l, long i, char *filename)
{

  //  a=0 non include "l" nel suffisso
  //  a=1 include "l" nel suffisso
  //  l:layer
  //  i:temporal step

  char SSSSLLLLL[ ]= {"SSSSLLLLL"};
  char SSSS[ ]= {"SSSS"};
  char *name=NULL;

  if (a==0)
    {
      write_suffix(SSSS, i, 0);
      name=join_strings(filename,SSSS);
    }
  else if (a==1)
    {
      write_suffix(SSSSLLLLL, i, 0);
      write_suffix(SSSSLLLLL, l, 5);
      name=join_strings(filename,SSSSLLLLL);
    }
  else
    {
      t_error("Value not admitted");
    }

  rename_map(name);

  free(name);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void rename_map(char *filename)
{

  char *name, *namenw, *temp;

  if (existing_file(filename) == 2)
    {
      name = join_strings(filename, ascii_grass);
      temp = join_strings(filename, ".old");
      namenw = join_strings(temp, ascii_grass);
      rename(name, namenw);
      free(namenw);
      free(name);
      free(temp);
    }
  else if (existing_file(filename) == 3)
    {
      name = join_strings(filename, ascii_esri);
      temp = join_strings(filename, ".old");
      namenw = join_strings(temp, ascii_esri);
      rename(name, namenw);
      free(namenw);
      free(name);
      free(temp);
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_tensorseries2(char *suf, long l, char *filename, short type,
                         short format, DOUBLETENSOR *T, T_INIT *UV, long novalue)
{

  char LLLLL[ ]= {"LLLLL"};
  char *temp1, *temp2;
  long r, c;
  DOUBLEMATRIX *M;

  temp1 = join_strings(LLLLL, suf);
  write_suffix(temp1, l, 1);

  M = new_doublematrix(T->nrh,T->nch);

  for (r=1; r<=T->nrh; r++)
    {
      for (c=1; c<=T->nch; c++)
        {
          M->co[r][c] = T->co[l][r][c];
        }
    }

  temp2 = join_strings(filename, temp1);
  write_map(temp2, type, format, M, UV, novalue);

  free_doublematrix(M);
  free(temp1);
  free(temp2);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_tensorseries2_vector(char *suf, long l, char *filename, short type,
                                short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr,
                                long nc)
{

  char LLLLL[ ]= {"LLLLL"};
  char *temp1, *temp2;
  long i, npoints=T->nch;
  std::unique_ptr<Vector<double>> V;

  temp1 = join_strings(LLLLL, suf);
  write_suffix(temp1, l, 1);

  V.reset(new Vector<double>{npoints});

  for (i=1; i<=npoints; i++)
    {
      V->co[i] = T->co[l][i];
    }

  temp2 = join_strings(filename, temp1);
  write_map_vector(temp2, type, format, V.get(), UV, novalue, J, nr, nc);

  free(temp1);
  free(temp2);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_tensorseries3(char *suffix, char *filename, short type,
                         short format, DOUBLETENSOR *T, T_INIT *UV, long novalue)
{

  long l;
  for (l=T->ndl; l<=T->ndh; l++)
    {
      write_tensorseries2(suffix, l, filename, type, format, T, UV, novalue);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_tensorseries3_vector(char *suffix, char *filename, short type,
                                short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr,
                                long nc)
{

  long l;
  for (l=T->nrl; l<=T->nrh; l++)
    {
      write_tensorseries2_vector(suffix, l, filename, type, format, T, UV, novalue,
                                 J, nr, nc);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
