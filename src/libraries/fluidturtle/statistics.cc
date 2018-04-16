#include "turtle.h"
#include  "t_statistics.h"
#define MAXCOUNTER 3

/*----------------------------------------------------------------------------------*/
float longvector_n_moment(LONGVECTOR *v, float mean,float NN, long novalue)
{
  unsigned i,n=0;
  float  moment=0;


  if (v==NULL || v->co==NULL )
    {
      t_error("this vector was never allocated");
    }
  else if (v->nh <1 ||  v->isdynamic !=1)
    {
      t_error("this vector was not properly allocated");
    }


  if (NN==1)
    {

      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=v->co[i];
              n++;
            }

        }
      moment/=n;

    }
  else if (NN==2)
    {

      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=(v->co[i])*(v->co[i]);
              n++;
            }

        }

      moment=(moment/n-mean*mean);

    }
  else
    {
      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=pow((v->co[i]-mean),NN);
              n++;
            }

        }

      moment/=n;

    }
  return moment;
}

/*----------------------------------------------------------------------------------*/
float floatvector_n_moment(FLOATVECTOR *v, float mean,float NN, float novalue)
{
  unsigned i,n=0;
  float  moment=0;


  if (v==NULL || v->co==NULL )
    {
      t_error("this vector was never allocated");
    }
  else if (v->nh <1 ||  v->isdynamic !=1)
    {
      t_error("this vector was not properly allocated");
    }


  if (NN==1)
    {

      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=v->co[i];
              n++;
            }

        }
      moment/=n;

    }
  else if (NN==2)
    {

      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=(v->co[i])*(v->co[i]);
              n++;
            }

        }

      moment=(moment/n-mean*mean);

    }
  else
    {
      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=pow((v->co[i]-mean),NN);
              n++;
            }

        }

      moment/=n;

    }
  return moment;
}


/*----------------------------------------------------------------------------------*/
double doublevector_n_moment(DOUBLEVECTOR *v, float mean,float NN,
                             double novalue)
{
  unsigned i,n=0;
  double  moment=0;


  if (v==NULL || v->co==NULL )
    {
      t_error("this vector was never allocated");
    }
  else if (v->nh <1 ||  v->isdynamic !=1)
    {
      t_error("this vector was not properly allocated");
    }


  if (NN==1)
    {

      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=v->co[i];
              n++;
            }

        }
      moment/=n;

    }
  else if (NN==2)
    {

      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=(v->co[i])*(v->co[i]);
              n++;
            }

        }

      moment=(moment/n-mean*mean);

    }
  else
    {
      for (i=1; i<v->nh; i++)
        {

          if (v->co[i]!=novalue)
            {
              moment+=pow((v->co[i]-mean),NN);
              n++;
            }

        }

      moment/=n;

    }
  return moment;
}

/*----------------------------------------------------------------------------------*/

/*... Double precision variables' vector n-moment estimation  */
float floatmatrix_n_moment(FLOATMATRIX *m, float mean,float NN, float novalue)
{
  long i,j;
  float  moment=0, n;


  if (m==NULL || m->co==NULL )
    {
      t_error("this matrix was never allocated");
    }
  else if (m->nrh <1 || m->nch <1 || m->isdynamic !=1)
    {
      t_error("this matrix was not properly allocated");
    }

  n=0;

  if (NN==1)
    {
      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=m->co[i][j];
                  n++;
                }

            }
        }

      moment/=n;

    }
  else if (NN==2)
    {

      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=(m->co[i][j])*(m->co[i][j]);
                  n++;
                }

            }
        }

      moment=(moment/n-mean*mean);

    }
  else
    {
      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=pow((m->co[i][j]-mean),NN);
                  n++;
                }
            }
        }

      moment/=n;

    }
  return moment;
}

/*----------------------------------------------------------------------------------*/

double doublematrix_n_moment(DOUBLEMATRIX *m, double mean,double  NN,
                             double novalue)
{
  long i,j;
  double  moment=0, n;


  if (m==NULL || m->co==NULL )
    {
      t_error("this matrix was never allocated");
    }
  else if (m->nrh <1 || m->nch <1 || m->isdynamic !=1)
    {
      t_error("this matrix was not properly allocated");
    }

  n=0;

  if (NN==1)
    {

      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=m->co[i][j];
                  n++;
                }
            }
        }
      moment/=n;

    }
  else if (NN==2)
    {

      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=(m->co[i][j])*(m->co[i][j]);
                  n++;
                }
            }
        }

      moment=(moment/n-mean*mean);

    }
  else
    {
      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=pow((m->co[i][j]-mean),NN);
                  n++;
                }
            }
        }

      moment/=n;

    }
  return moment;
}



/*----------------------------------------------------------------------------------*/

float longmatrix_n_moment(LONGMATRIX *m, double mean,double  NN, long novalue)
{
  long i,j,n=0;
  float  moment=0;


  if (m==NULL || m->co==NULL )
    {
      t_error("this matrix was never allocated");
    }
  else if (m->nrh <1 || m->nch <1 || m->isdynamic !=1)
    {
      t_error("this matrix was not properly allocated");
    }


  if (NN==1)
    {

      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=m->co[i][j];
                  n++;
                }
            }
        }
      moment/=n;

    }
  else if (NN==2)
    {

      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=(m->co[i][j])*(m->co[i][j]);
                  n++;
                }
            }
        }

      moment=(moment/n-mean*mean);

    }
  else
    {
      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue)
                {
                  moment+=pow((m->co[i][j]-mean),NN);
                  n++;
                }
            }
        }

      moment/=n;

    }
  return moment;
}

/*----------------------------------------------------------------------------------*/

float floatmatrix_restricted_n_moment(FLOATMATRIX *m,FLOATMATRIX *s,
                                      float mean,float NN, float novalue,float novalue2)
{
  long i,j;
  float  moment=0, n;


  if (m==NULL || m->co==NULL || s==NULL || s->co==NULL)
    {
      t_error("these matrices was never allocated");
    }
  else if (m->nrh <1 || m->nch <1 || m->isdynamic !=1 || s->nrh <1 || s->nch <1
           || s->isdynamic !=1)
    {
      t_error("thess matrices was not properly allocated");

    }
  else if (m->nrh !=s->nrh || m->nch!=s->nch)
    {
      t_error("the matrixes do not have the same dimensions");
    }

  n=0;

  if (NN==1)
    {
      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue && s->co[i][j]!=novalue2)
                {
                  moment+=m->co[i][j];
                  n++;
                }

            }
        }

      moment/=n;

    }
  else if (NN==2)
    {

      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue && s->co[i][j]!=novalue2)
                {
                  moment+=(m->co[i][j])*(m->co[i][j]);
                  n++;
                }

            }
        }

      moment=(moment/n-mean*mean);

    }
  else
    {
      for (i=1; i<m->nrh; i++)
        {
          for (j=1; j<m->nch; j++)
            {

              if (m->co[i][j]!=novalue && s->co[i][j]!=novalue2)
                {
                  moment+=pow((m->co[i][j]-mean),NN);
                  n++;
                }
            }
        }

      moment/=n;

    }
  return moment;
}

/*----------------------------------------------------------------------------------*/
float longvector_correlation(LONGVECTOR *v,LONGVECTOR *u,float m1,float m2,
                             long r,long novalue)
{
  long i,n=0;
  float correlation=0;


  if (v==NULL || v->co==NULL || u==NULL || u->co==NULL )
    {
      t_error("this vector was never allocated");
    }
  else if (v->nh <1 ||  v->isdynamic !=1 || u->nh <1 ||  u->isdynamic !=1)
    {
      t_error("this vector was not properly allocated");
    }




  for (i=1; i<v->nh-r; i++)
    {

      if (v->co[i]!=novalue && u->co[i]!=novalue)
        {
          correlation+=v->co[i]*u->co[i+r];
          n++;
        }
    }




  return correlation/n -m1*m2;



}


/*----------------------------------------------------------------------------------*/
float floatvector_correlation(FLOATVECTOR *v,FLOATVECTOR *u,float m1,float m2,
                              long r,float novalue)
{
  long i,n=0;
  float correlation=0;


  if (v==NULL || v->co==NULL || u==NULL || u->co==NULL )
    {
      t_error("this vector was never allocated");
    }
  else if (v->nh <1 ||  v->isdynamic !=1 || u->nh <1 ||  u->isdynamic !=1)
    {
      t_error("this vector was not properly allocated");
    }




  for (i=1; i<v->nh-r; i++)
    {

      if (v->co[i]!=novalue && u->co[i]!=novalue)
        {
          correlation+=v->co[i]*u->co[i+r];
          n++;
        }
    }




  return correlation/n -m1*m2;



}


/*----------------------------------------------------------------------------------*/
double doublevector_correlation(DOUBLEVECTOR *v,DOUBLEVECTOR *u,double m1,
                                double m2,long r,double novalue)
{
  long i,n=0;
  float correlation=0;


  if (v==NULL || v->co==NULL || u==NULL || u->co==NULL )
    {
      t_error("this vector was never allocated");
    }
  else if (v->nh <1 ||  v->isdynamic !=1 || u->nh <1 ||  u->isdynamic !=1)
    {
      t_error("this vector was not properly allocated");
    }




  for (i=1; i<v->nh-r; i++)
    {

      if (v->co[i]!=novalue && v->co[i]!=novalue)
        {
          correlation+=v->co[i]*u->co[i+r];
          n++;
        }
    }




  return correlation/n -m1*m2;



}


/*----------------------------------------------------------------------------------*/

double double_n_moment(double *m, long nh,double mean,double  NN,
                       double novalue)
{
  long i;
  double  moment=0, n;


  if (m==NULL  )
    {
      t_error("this matrix was never allocated");
    }
  else if (nh <1)
    {
      t_error("this matrix was not properly allocated");
    }

  n=0;

  if (NN==1)
    {

      for (i=1; i<=nh; i++)
        {

          if (m[i]!=novalue)
            {
              moment+=m[i];
              n++;
              /*printf("(%f,%f,%f)",m[i],moment,n);
              scanf("%hd/n",&ch);*/
            }

        }

      if (n >=1 )
        {
          moment/=n;
          /* printf("%f %f",moment,n);
          scanf("%hd/n",&ch);*/
        }
      else
        {
          printf("\nWarning::No valid data were processed\n");
          printf("\nWarning::setting moment value to zero\n");
          moment=0;
        }

    }
  else if (NN==2)
    {
      for (i=1; i<=nh; i++)
        {
          if (m[i]!=novalue)
            {
              moment+=(m[i])*(m[i]);
              n++;
              /*printf("(%f,%f,%f)",m[i],moment,n);
              scanf("%hd/n",&ch);*/
            }

        }
      if (n >=1 )
        {
          moment=(moment/n-mean*mean);
          /*printf("%f %f %f",moment,n,mean);
          scanf("%hd/n",&ch);*/
        }
      else
        {
          printf("\nWarning::No valid data were processed\n");
          printf("\nWarning::setting moment value to zero\n");
          moment=0;
        }


    }
  else
    {
      for (i=1; i<=nh; i++)
        {
          if (m[i]!=novalue)
            {
              moment+=pow((m[i]-mean),NN);
              n++;
              /*printf("%f %f %f %f %f",m[i],moment,n,mean,NN);
              scanf("%hd/n",&ch);*/
            }

        }
      if (n >=1 )
        {
          moment/=n;
          /*printf("%f %f %f",moment,n);
          scanf("%hd/n",&ch);*/

        }
      else
        {
          printf("\nWarning::No valid data were processed\n");
          printf("\nWarning::setting moment value to zero\n");
          moment=0;
        }

    }

  /* printf("%f\n",moment); */
  return moment;
}


