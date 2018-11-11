#include "turtle.h"
#include "geomorphology.0875.h"
#include "t_datamanipulation.h"
#include "t_random.h"
#include "networks.h"
#include <logger.h>
#include <timer.h>

#define Pi 3.14159265358979     /* P greco */


//***************************************************************************


/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of parts in which one wants to divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti
                                                 */
void sky_view_factor(Matrix<double> *sky, long N, T_INIT *UV, Matrix<double> *input, Matrix<short> *convess,
                     long novalue)
{
    GEOTIMER_PREFIX(__func__);
    long i,j,t,m,n,p,q,h,k,r,s; //counter
    double deltateta; //amplitude of the angles in which the horizon is divided
    std::unique_ptr<Matrix<double>> alfa; //matrices with the angles of the direction
    std::unique_ptr<Vector<double>>
            vv; //vector with the view factor of the current pixel for one of the N parts
    std::unique_ptr<Vector<double>>
            v; //vector with the minimum view factor of the current pixel for one of the N parts
    double vvv; //mean of the sky view for a pixel of the N parts

    if (sky->nrh!=input->nrh)
        t_error("Sky view factor fatal error, number of rows not consistent");
    if (sky->nch!=input->nch)
        t_error("Sky view factor fatal error, number of cols not consistent");

    // Computation of the matrix with the angles of the direction
    alfa.reset(new Matrix<double>{2*input->nrh-1,2*input->nch-1});
    *alfa = (double)novalue; //initialisation with novalue
    for (i=1; i<=2*input->nrh-1; i++)
    {
        for (j=1; j<=2*input->nch-1; j++)
        {
            if (i<=input->nrh && j<input->nch)
            {
                (*alfa)(i,j)=3.0/2.0*Pi+atan(((input->nrh-i)*  (*UV->U)(1))/ ((input->nch-j)*(*UV->U)(1)));
            }
            if (i>input->nrh && j<=input->nch)
            {
                (*alfa)(i,j)=Pi+atan(((input->nch-j)*(*UV->U)(1))/
                                     ((i-input->nrh)*(*UV->U)(1)));
            }
            if (i>=input->nrh && j>input->nch)
            {
                (*alfa)(i,j)=Pi/2.0+atan(((i-input->nrh)*(*UV->U)(1))/
                                         ((j-input->nch)*(*UV->U)(1)));
            }
            if (i<input->nrh && j>=input->nch)
            {
                (*alfa)(i,j)=atan(((j-input->nch)*(*UV->U)(1))/
                                  ((input->nrh-i)*(*UV->U)(1)));
            }
        }
    }

    // Computation of matrix with sky view factor:
    for (i=1; i<=sky->nrh; i++)
    {
        for (j=1; j<=sky->nch; j++)
        {
            (*sky)(i,j)=(double)novalue;
        }
    }

    v.reset(new Vector<double>{N});
    vv.reset(new Vector<double>{N});
    deltateta=2.0*Pi/N;

    for (i=1; i<=input->nrh; i++)
    {
        for (j=1; j<=input->nch; j++)
        {
            if ((long)(*input)(i,j)!=novalue)  //computation only of novalue pixels
            {
                for (t=1; t<=N; t++)
                {
                    v->co[t]=1.0;
                }
                m=input->nrh-i+1;
                n=input->nch-j+1;
                p=m+input->nrh-1;
                q=n+input->nch-1;
                for (h=m; h<=p; h++)
                {
                    for (k=n; k<=q; k++)
                    {
                        for (t=1; t<=N; t++)
                        {
                            if ((*alfa)(h,k)>=(t-1)*deltateta && (*alfa)(h,k)<t*deltateta)
                            {
                                r=h-m+1;
                                s=k-n+1;
                                if ((*convess)(r,s)==1 && sqrt(pow((r-i),2)+pow((s-j),2))!=0)
                                {
                                    vv->co[t]=1-sin(atan(((*input)(r,s)-(*input)(i,j))
                                                         /(sqrt(pow((r-i),2)+pow((s-j),2))*(*UV->U)(1))));
                                    if (vv->co[t]<v->co[t])
                                    {
                                        v->co[t]=vv->co[t];
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                vvv=0.0;
                for (t=1; t<=N; t++)
                {
                    vvv=vvv+v->co[t];
                }
                (*sky)(i,j)=(1.0/N*vvv);
            }
        }
        geolog << "Percentage of the calculation of the sky view factor matrix: %5.2f%%\n"
               << 100.0*(double)i/(double)sky->nrh << std::endl;
    }
}

//***************************************************************************

//***************************************************************************

//Presa da geomorphology099 e modificato der_min
void nablaquadro_mask(Matrix<double> *Z0, Matrix<short> *curv, Vector<double> *U, Vector<double> *V)
{

    short y;
    long i,j,h,rows,cols;
    double grid[9],z[9],derivate2;
    double der_min=0.00001; /*limite per la limite per la planarita'*/

    short v[13][2] = {           { 0, 0},
                                 { 0, 1},
                                 {-1, 1},
                                 {-1, 0},
                                 {-1,-1},
                                 { 0,-1},
                                 { 1,-1},
                                 { 1, 0},
                                 { 1, 1},
                                 { 0, 0},
                                 { 0, 0},
                                 { 0, 0},
                                 { 0, 0}
    };

    grid[0]=0;
    grid[1]=grid[5]=U->co[1];
    grid[3]=grid[7]=U->co[2];
    grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);

    rows=Z0->nrh;
    cols=Z0->nch;

    for (i=2; i<=rows-1; i++)
    {
        for (j=2; j<=cols-1; j++)
        {
            z[0]=(*Z0)(i,j);
            if (z[0]!=V->co[2])
            {
                y=1;
                for (h=1; h<=8; h++)
                {
                    z[h]=(*Z0)(i+v[h][0],j+v[h][1]);
                    if (z[h]==V->co[2])
                    {
                        y=0;
                        break;
                    }
                }
                if (y==0)
                {
                    (*curv)(i,j)=1;
                }
                else
                {
                    derivate2=0.5*((z[1]+z[5]-2*z[0])/(grid[1]*grid[1])+ (z[3]+z[7]-2*z[0])/
                                                                         (grid[3]*grid[3]));
                    derivate2=derivate2+0.5*((z[2]+z[4]+z[6]+z[8]-4*z[0])/(grid[6]*grid[6]));

                    if (fabs(derivate2)<=der_min || derivate2>der_min)  //plane or concave
                    {
                        (*curv)(i,j)=0;
                    }
                    else
                    {
                        (*curv)(i,j)=1; //convex
                    }
                }
            }
        }
    }
}
//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

//***************************************************************************

