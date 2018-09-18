
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


//Author: Stefano Endrizzi
//Date: 13 November 2005
//Contents: Snow subroutines

#include "constants.h"
#include "struct.geotop.h"
#include "snow.h"
#include "extensions.h"
#include "rw_maps.h"
#include "times.h"
#include "output.h"

#include "PBSM.h"

extern long number_novalue, number_absent;

extern T_INIT *UV;
extern const char *WORKING_DIRECTORY;
extern char **files;
extern long Nl, Nr, Nc;
extern long i_sim, i_run;
extern char *FailedRunFile;

#define simpl_snow 1E-1
#define no_snow 1E-6

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//Jordan et al., 1999
double rho_newlyfallensnow(double u, double Tatm, double Tfreez)
{

  double rho,T;

  if (Tatm>260.15-tk)
    {
      T=Tatm;
      if (T>=275.65-tk) T=275.65-tk;
      rho=500.0*(1.0-0.951*exp(-1.4*pow(278.15-tk-T,-1.15)-0.008*pow(u,1.7)));
    }
  else
    {
      rho=500.0*(1.0-0.904*exp(-0.008*pow(u,1.7)));
    }

  /*if(Tatm>Tfreez+2.0){
   rho=50.0+1.7*pow(17,1.5);
   }else if(Tatm<=Tfreez+2.0 && Tatm>Tfreez-15.0){
   rho=50.0+1.7*pow(Tatm-Tfreez+15.0,1.5);
   }else{
   rho=50.0;
   }*/

  return (rho);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snow_compactation(double Dt, long r, long c, long l, STATEVAR_3D *snow,
                       double slope, PAR *par)
{

  long m;
  double theta_i,theta_w,c1,c2,c3,c4,c5,load,eta0,eta,CR1,CR2;

  theta_i=snow->w_ice->co[l][r][c]/(0.001*(*snow->Dzl)(l,r,c)*rho_i);
  theta_w=snow->w_liq->co[l][r][c]/(0.001*(*snow->Dzl)(l,r,c)*rho_w);

  if ( theta_i < par->snow_maxpor)
    {

      //DESTRUCTIVE METAMORPHISM
      //gamma=theta*rho
      if (theta_i*rho_i<=par->snow_density_cutoff)
        {
          c1=par->drysnowdef_rate;
        }
      else
        {
          c1=par->drysnowdef_rate*exp(-0.046*(rho_i*theta_i -
                                              par->snow_density_cutoff));
        }
      if (theta_w>0.001) c1*=par->wetsnowdef_rate;
      c2=2.777E-6; //[s^-1]
      c3=0.04;   //[K^-1]
      CR1=-c1*c2*exp(-c3*(Tfreezing-(*snow->T)(l,r,c)));

      //OVERBURDEN
      eta0=par->snow_viscosity;  //[kg s m^-2]
      c4=0.08;   //[K^-1]
      c5=0.023;    //[m^3 kg^-1]
      load=0.0;    //[kg m^-2]
      for (m=l; m<=(*snow->lnum)(r,c); m++)
        {
          load+=(snow->w_ice->co[m][r][c]+snow->w_liq->co[m][r][c]);
        }
      load*=fabs(cos(slope*Pi/180.));
      eta=eta0*exp(c4*(Tfreezing-(*snow->T)(l,r,c))+c5*(rho_i*theta_i));
      CR2=-load/eta;

      (*snow->Dzl)(l,r,c) *= exp( (CR1 + CR2) * Dt ); // Cuell

    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snow_layer_combination(double a, long r, long c, STATEVAR_3D *snow, double Ta, Vector<long> *inf, double SWEmax_layer,
                            double SWEmax_tot)

{

  long l, linf, n, k, j, max=snow->Dzl->ndh;
  double D, Dnew, Dmin, SWE, SWEnew, ice;
  short occurring;
  FILE *f;

  //check on SWEmax[kg/m2]
  SWE = 0.0;
  l=1;
  while (l <= (*snow->lnum)(r,c))
    {
      SWE += snow->w_ice->co[l][r][c];
      ice = snow->w_ice->co[l][r][c];
      snow->w_ice->co[l][r][c] -= Fmax(0.0, SWE-SWEmax_tot);
      if (snow->w_ice->co[l][r][c] < 0) snow->w_ice->co[l][r][c] = 0.0;
      if (ice > 0)
        {
          (*snow->Dzl)(l,r,c) *= snow->w_ice->co[l][r][c]/ice;
          snow->w_liq->co[l][r][c] *= snow->w_ice->co[l][r][c]/ice;
        }
      l++;
    }

  //D=snow depth(mm)
  D = 0.0;
  SWE = 0.0;
  for (l=1; l<=max; l++)
    {
      D += (*snow->Dzl)(l,r,c);
      SWE += (snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c]);
    }

  //PREPROCESSING
  //1. If the snow depth is too small, then it is reset to 0
  if (SWE < no_snow*SWEmax_layer)
    {
      (*snow->lnum)(r,c)=0;
      (*snow->type)(r,c)=0;
      for (l=1; l<=max; l++)
        {
          (*snow->Dzl)(l,r,c)=0.0;
          snow->w_ice->co[l][r][c]=0.0;
          snow->w_liq->co[l][r][c]=0.0;
          (*snow->T)(l,r,c)=0.0; //Temperatura di inizializzazione
        }


      //2. If D<Dsnow_simpl, we are in the simplified case
    }
  else if ((*snow->lnum)(r,c) > 0 && SWE < simpl_snow*SWEmax_layer)
    {

      (*snow->type)(r,c)=1;
      if ((*snow->lnum)(r,c)>1)
        {
          for (l=(*snow->lnum)(r,c); l>1; l--)
            {
              snowlayer_merging(a, r, c, snow, l, l-1, l-1);
            }
          for (l=2; l<=max; l++)
            {
              (*snow->T)(l,r,c)=0.0;
              (*snow->Dzl)(l,r,c)=0.0;
              snow->w_liq->co[l][r][c]=0.0;
              snow->w_ice->co[l][r][c]=0.0;
            }
          (*snow->lnum)(r,c)=1;
        }

      //3. if z>=Dsnow_simpl, ordinary case
    }
  else if ((*snow->lnum)(r,c) > 0 && SWE >= simpl_snow*SWEmax_layer)
    {

      (*snow->type)(r,c)=2;

      //4. if there is not yet a snow layer and D<Dsnow_simpl, simplified case
    }
  else if ((*snow->lnum)(r,c) == 0 && SWE < simpl_snow*SWEmax_layer)
    {

      (*snow->lnum)(r,c)=1;
      (*snow->type)(r,c)=1;
      (*snow->T)(1,r,c)=Fmin(Ta,-0.1);


      //5. if there is not yet a snow layer and D>=Dmin(max), simplified case
    }
  else if ((*snow->lnum)(r,c) == 0 && SWE >= simpl_snow*SWEmax_layer)
    {

      (*snow->lnum)(r,c)=1;
      (*snow->type)(r,c)=2;
      (*snow->T)(1,r,c)=Fmin(Ta,-0.1);

    }

  // SIMMETRICAL PARAMETERIZATION SCHEME (new)
  if ((*snow->type)(r,c)==2)
    {

      //remove layers < 0.01 mm
      n = 0;
      do
        {
          occurring = 0;
          if ((*snow->lnum)(r,c) > 1)
            {
              for (l=1; l<=(*snow->lnum)(r,c); l++)
                {
                  if (snow->w_ice->co[l][r][c] < simpl_snow*SWEmax_layer)
                    {
                      merge_layers(a, r, c, snow, l);
                      occurring = 1;
                    }
                }
            }
          n ++;
        }
      while (n<=max && occurring==1);

      //add new layer
      if (snow->w_ice->co[(*snow->lnum)(r,c)][r][c] > SWEmax_layer *
          (1.+simpl_snow) )
        {
          occurring = 1;
        }
      else
        {
          occurring = 0;
        }

      if (occurring == 1)
        {
          if ((*snow->lnum)(r,c) == max)
            {
              linf = 0;
              Dmin = 9.E99;
              for (n=1; n<=inf->nh; n++)
                {
                  if (inf->co[n]>0 && inf->co[n]<max)
                    {
                      if (Dmin > snow->w_ice->co[inf->co[n]][r][c] + snow->w_ice->co[inf->co[n]
                          +1][r][c])
                        {
                          Dmin = snow->w_ice->co[inf->co[n]][r][c] + snow->w_ice->co[inf->co[n]
                                                                                     +1][r][c];
                          linf = inf->co[n];
                        }
                    }
                  else if (inf->co[n]>1 && inf->co[n]<=max)
                    {
                      if (Dmin > snow->w_ice->co[inf->co[n]][r][c] + snow->w_ice->co[inf->co[n]
                          -1][r][c])
                        {
                          Dmin = snow->w_ice->co[inf->co[n]][r][c] + snow->w_ice->co[inf->co[n]
                                                                                     -1][r][c];
                          linf = -inf->co[n];
                        }
                    }
                }

              if (linf > 0)
                {
                  snowlayer_merging(a, r, c, snow, linf, linf+1, linf);
                }
              else if (linf < 0)
                {
                  snowlayer_merging(a, r, c, snow, -linf, -linf-1, -linf-1);
                  linf = -linf-1;
                }
              else
                {
                  f = fopen(FailedRunFile, "w");
                  fprintf(f,"r:%ld c:%ld \n",r,c);
                  fprintf(f,
                          "Error in snow combination - Rules to combine layers not applicable\n");
                  fclose(f);
                  t_error("Fatal Error! Geotop is closed. See failing report.");
                }

              for (l=linf+1; l<(*snow->lnum)(r,c); l++)
                {
                  snow->w_ice->co[l][r][c] = snow->w_ice->co[l+1][r][c];
                  snow->w_liq->co[l][r][c] = snow->w_liq->co[l+1][r][c];
                  (*snow->T)(l,r,c) = (*snow->T)(l+1,r,c);
                  (*snow->Dzl)(l,r,c) = (*snow->Dzl)(l+1,r,c);
                }

              initialize_snow(r, c, (*snow->lnum)(r,c), snow);

              (*snow->lnum)(r,c) --;
            }

          snow->w_liq->co[(*snow->lnum)(r,c)+1][r][c] =
            snow->w_liq->co[(*snow->lnum)(r,c)][r][c] *
            (snow->w_ice->co[(*snow->lnum)(r,c)][r][c] - SWEmax_layer) /
            snow->w_ice->co[(*snow->lnum)(r,c)][r][c];
          (*snow->Dzl)((*snow->lnum)(r,c)+1,r,c) =
            (*snow->Dzl)((*snow->lnum)(r,c),r,c) *
            (snow->w_ice->co[(*snow->lnum)(r,c)][r][c] - SWEmax_layer) /
            snow->w_ice->co[(*snow->lnum)(r,c)][r][c];
          snow->w_ice->co[(*snow->lnum)(r,c)+1][r][c] =
            snow->w_ice->co[(*snow->lnum)(r,c)][r][c] - SWEmax_layer;
          (*snow->T)((*snow->lnum)(r,c)+1,r,c) =
            (*snow->T)((*snow->lnum)(r,c),r,c);

          snow->w_ice->co[(*snow->lnum)(r,c)][r][c] -=
            snow->w_ice->co[(*snow->lnum)(r,c)+1][r][c];
          snow->w_liq->co[(*snow->lnum)(r,c)][r][c] -=
            snow->w_liq->co[(*snow->lnum)(r,c)+1][r][c];
          (*snow->Dzl)((*snow->lnum)(r,c),r,c) -=
            (*snow->Dzl)((*snow->lnum)(r,c)+1,r,c);

          (*snow->lnum)(r,c) ++;
        }

      //split layers
      if ((*snow->lnum)(r,c) < max)
        {

          do
            {

              occurring = 0;

              //FROM UP DOWN
              k = inf->nh;

              do
                {

                  l = inf->co[k];

                  if (snow->w_ice->co[l][r][c] > SWEmax_layer*2.)
                    {

                      occurring = 1;

                      for (j=(*snow->lnum)(r,c); j>l; j--)
                        {
                          snow->w_ice->co[j+1][r][c] = snow->w_ice->co[j][r][c];
                          (*snow->Dzl)(j+1,r,c) = (*snow->Dzl)(j,r,c);
                          snow->w_liq->co[j+1][r][c] = snow->w_liq->co[j][r][c];
                          (*snow->T)(j+1,r,c) = (*snow->T)(j,r,c);
                        }

                      (*snow->Dzl)(l+1,r,c) = (*snow->Dzl)(l,r,c) * SWEmax_layer /
                                                 snow->w_ice->co[l][r][c];
                      snow->w_liq->co[l+1][r][c] = snow->w_liq->co[l][r][c] * SWEmax_layer /
                                                   snow->w_ice->co[l][r][c];
                      snow->w_ice->co[l+1][r][c] = SWEmax_layer;

                      (*snow->T)(l+1,r,c) = (*snow->T)(l,r,c);

                      (*snow->Dzl)(l,r,c) -= (*snow->Dzl)(l+1,r,c);
                      snow->w_liq->co[l][r][c] -= snow->w_liq->co[l+1][r][c];
                      snow->w_ice->co[l][r][c] -= snow->w_ice->co[l+1][r][c];

                      (*snow->lnum)(r,c) ++;
                    }

                  k --;

                }
              while (occurring == 0 && k > 0 );

            }
          while ((*snow->lnum)(r,c) != max && occurring != 0);

        }

      //check
      Dnew = 0.0;
      SWEnew = 0.0;
      for (l=1; l<=max; l++)
        {
          Dnew += (*snow->Dzl)(l,r,c);
          SWEnew += (snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c]);
        }

      if (fabs(D-Dnew)>0.001 || fabs(SWE-SWEnew)>0.001)
        {
          f = fopen(FailedRunFile, "w");
          fprintf(f,"r:%ld c:%ld Dold:%f Dnew:%f SWEold:%f SWEnew:%f\n",r,c,D,Dnew,SWE,
                  SWEnew);
          fprintf(f,"Error in snow combination\n");
          fclose(f);
          t_error("Fatal Error! Geotop is closed. See failing report.");
        }

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double DEPTH(long r, long c, Matrix<long> *n, Tensor<double> *Dz)
{

  double d=0.0;
  long l;

  for (l=1; l<=(*n)(r,c); l++)
    {
      d+= (*Dz)(l,r,c);
    }

  return (d);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double get_SWE(long r, long c, Matrix<long> *n, DOUBLETENSOR *w1, DOUBLETENSOR *w2)
{

  double d=0.0;
  long l;

  for (l=1; l<=(*n)(r,c); l++)
    {
      d+=(w1->co[l][r][c]+w2->co[l][r][c]);
    }

  return (d);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snowlayer_merging(double a, long r, long c, STATEVAR_3D *snow, long l1,
                       long l2, long lres)
{

  double h;

  h=internal_energy(snow->w_ice->co[l1][r][c],snow->w_liq->co[l1][r][c],
                    (*snow->T)(l1,r,c))+internal_energy(snow->w_ice->co[l2][r][c],
                                                           snow->w_liq->co[l2][r][c],(*snow->T)(l2,r,c));
  (*snow->Dzl)(lres,r,c)=(*snow->Dzl)(l1,r,c)+(*snow->Dzl)(l2,r,c);
  snow->w_ice->co[lres][r][c]=snow->w_ice->co[l1][r][c]
                              +snow->w_ice->co[l2][r][c];
  snow->w_liq->co[lres][r][c]=snow->w_liq->co[l1][r][c]
                              +snow->w_liq->co[l2][r][c];
  if ((*snow->Dzl)(lres,r,c)<0 || snow->w_ice->co[lres][r][c]<0
      || snow->w_liq->co[lres][r][c]<0)
    {
      printf("ERROR 1 in snow layer merging r:%ld c:%ld l1:%ld l2:%ld lres:%ld\n",r,
             c,l1,l2,lres);
      write_snow_all(r, c, snow);
      t_error("Stop Execution");
    }
  from_internal_energy(a, r+6000,c+6000,h,&(snow->w_ice->co[lres][r][c]),
                       &(snow->w_liq->co[lres][r][c]),&((*snow->T)(lres,r,c)));
  if ((*snow->T)(lres,r,c)>0)
    {
      printf("ERROR 2 in snow layer merging r:%ld c:%ld l1:%ld l2:%ld lres:%ld\n",r,
             c,l1,l2,lres);
      write_snow_all(r, c, snow);
      t_error("Stop Execution");
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double internal_energy(double w_ice, double w_liq, double T)
{

  return (c_ice*w_ice+c_liq*w_liq)*(T-Tfreezing) + Lf*w_liq;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void from_internal_energy(double a, long r, long c, double h, double *w_ice,
                          double *w_liq, double *T)
{

  double SWE=(*w_ice)+(*w_liq);
  double T0;
  double A=-h/(c_ice*SWE), B=c_liq/(c_ice*a*a), C=(Lf*SWE-h)/(c_ice*SWE*a*a);
  long cont=0;

  if (SWE>0)
    {
      *T=0.0;
      do
        {
          T0=*T;
          *T=T0-(pow(T0,3.0)+A*pow(T0,2.0)+B*T0+C)/(3.0*pow(T0,2.0)+2*A*T0+B);
          cont++;
          if (cont>100) printf("%ld %e %e %e %e\n",cont,SWE,h,T0,*T);
        }
      while (fabs(*T-T0)>1.E-10 && cont<100);

      *w_liq=theta_snow(a, 1., *T)*SWE;
      *w_ice=SWE-(*w_liq);
    }
  else
    {
      *T=0.0;
      *w_liq=0.0;
      *w_ice=0.0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_snow(long r, long c, long l, STATEVAR_3D *snow)
{

  printf("r:%ld c:%ld wice(%ld/%ld):%f wliq(%ld):%f T(%ld):%f Dz(%ld):%f\n",r,c,
         l,(*snow->lnum)(r,c),snow->w_ice->co[l][r][c],l,snow->w_liq->co[l][r][c],
         l,(*snow->T)(l,r,c),l,(*snow->Dzl)(l,r,c));

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void write_snow_all(long r, long c, STATEVAR_3D *snow)
{
  long l;
  for (l=1; l<=(*snow->lnum)(r,c); l++)
    {
      write_snow(r,c,l,snow);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short set_snow_min(double a, long r, long c, STATEVAR_3D *snow, long l1,
                   long l2, double Dmin)
{

  double h, f, dwl, dwi, dd;

  if ((*snow->Dzl)(l1,r,c) < Dmin
      && (*snow->Dzl)(l2,r,c) > 1.E-5)  //l1 too shallow and takes mass from l2
    {
      f = Fmin(Dmin - (*snow->Dzl)(l1,r,c),
               (*snow->Dzl)(l2,r,c))/(*snow->Dzl)(l2,r,c);
      h = internal_energy(snow->w_ice->co[l1][r][c], snow->w_liq->co[l1][r][c],
                          (*snow->T)(l1,r,c));
      dd = f*(*snow->Dzl)(l2,r,c);
      (*snow->Dzl)(l1,r,c) += dd;
      (*snow->Dzl)(l2,r,c) -= dd;
      dwl = f*snow->w_liq->co[l2][r][c];
      snow->w_liq->co[l1][r][c] += dwl;
      snow->w_liq->co[l2][r][c] -= dwl;
      dwi = f*snow->w_ice->co[l2][r][c];
      snow->w_ice->co[l1][r][c] += dwi;
      snow->w_ice->co[l2][r][c] -= dwi;
      h += internal_energy(dwi, dwl, (*snow->T)(l2,r,c));
      from_internal_energy(a, r+1000, c+1000, h, &(snow->w_ice->co[l1][r][c]),
                           &(snow->w_liq->co[l1][r][c]), &((*snow->T)(l1,r,c)));
      return 1;
    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short set_snow_max(double a, long r, long c, STATEVAR_3D *snow, long l1,
                   long l2, double Dmax)
{

  double h, f, dwl, dwi, dd;

  if ((*snow->Dzl)(l1,r,c) > Dmax) //l1 too thick and gives mass to l2
    {
      f = ((*snow->Dzl)(l1,r,c) - Dmax)/(*snow->Dzl)(l1,r,c);
      h = internal_energy(snow->w_ice->co[l2][r][c], snow->w_liq->co[l2][r][c], (*snow->T)(l2,r,c));
      dd = f*(*snow->Dzl)(l1,r,c);
      (*snow->Dzl)(l1,r,c) -= dd;
      (*snow->Dzl)(l2,r,c) += dd;
      dwl = f*snow->w_liq->co[l1][r][c];
      snow->w_liq->co[l1][r][c] -= dwl;
      snow->w_liq->co[l2][r][c] += dwl;
      dwi = f*snow->w_ice->co[l1][r][c];
      snow->w_ice->co[l1][r][c] -= dwi;
      snow->w_ice->co[l2][r][c] += dwi;
      h += internal_energy(dwi, dwl, (*snow->T)(l1,r,c));
      from_internal_energy(a, r+2000, c+2000, h, &(snow->w_ice->co[l2][r][c]),
                           &(snow->w_liq->co[l2][r][c]), &((*snow->T)(l2,r,c)));
      return 1;
    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short set_snowice_min(double a, long r, long c, STATEVAR_1D *snow, long l1,
                      long l2, double wicemin)
{

  //double h;
  double f, dwl, dwi, dd;

  if ((*snow->w_ice)(l1) < wicemin
      && (*snow->w_ice)(l2)> 1.E-6) //l1 too shallow and takes mass from l2
    {
      f = Fmin(wicemin - (*snow->w_ice)(l1),
               (*snow->w_ice)(l2))/(*snow->w_ice)(l2);
      //h = internal_energy((*snow->w_ice)(l1), snow->w_liq->co[l1], snow->T->co[l1]);
      dd = f * (*snow->Dzl)(l2);
        (*snow->Dzl)(l1) += dd;
        (*snow->Dzl)(l2) -= dd;
      dwl = f * (*snow->w_liq)(l2);
      (*snow->w_liq)(l1) += dwl;
      (*snow->w_liq)(l2) -= dwl;
      dwi = f*(*snow->w_ice)(l2);
      (*snow->w_ice)(l1) += dwi;
      (*snow->w_ice)(l2)-= dwi;
      //h += internal_energy(dwi, dwl, snow->T->co[l2]);
      //from_internal_energy(a, r+1000, c+1000, h, &(snow->w_ice->co[l1]), &(snow->w_liq->co[l1]), &(snow->T->co[l1]));
      return 1;
    }
  else
    {
      return 0;
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void split_layers(long r, long c, STATEVAR_3D *snow, long l1)
{

  long l;
  FILE *f;

  if (l1>(*snow->lnum)(r,c))
    {
      f = fopen(FailedRunFile, "w");
      fprintf(f,"Error 1 in split_layers\n");
      fclose(f);
      t_error("Fatal Error! Geotop is closed. See failing report.");
    }

  snow->w_ice->co[l1][r][c]*=0.5;
  snow->w_liq->co[l1][r][c]*=0.5;
  (*snow->Dzl)(l1,r,c)*=0.5;

  for (l=(*snow->lnum)(r,c); l>=l1; l--)
    {
      snow->w_ice->co[l+1][r][c]=snow->w_ice->co[l][r][c];
      snow->w_liq->co[l+1][r][c]=snow->w_liq->co[l][r][c];
      (*snow->T)(l+1,r,c)=(*snow->T)(l,r,c);
      (*snow->Dzl)(l+1,r,c)=(*snow->Dzl)(l,r,c);
    }

  (*snow->lnum)(r,c)+=1;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void merge_layers(double a, long r, long c, STATEVAR_3D *snow, long l1)
{

  long l;
  FILE *f;

  if (l1>(*snow->lnum)(r,c))
    {
      f = fopen(FailedRunFile, "w");
      fprintf(f,"Error 1 in merge_layers\n");
      fclose(f);
      t_error("Fatal Error! Geotop is closed. See failing report.");
    }

  if (l1==(*snow->lnum)(r,c))
    {
      snowlayer_merging(a, r, c, snow, l1, l1-1, l1-1);
      initialize_snow(r, c, l1, snow);
    }
  else
    {
      snowlayer_merging(a, r, c, snow, l1, l1+1, l1);
      for (l=l1+1; l<(*snow->lnum)(r,c); l++)
        {
          snow->w_ice->co[l][r][c]=snow->w_ice->co[l+1][r][c];
          snow->w_liq->co[l][r][c]=snow->w_liq->co[l+1][r][c];
          (*snow->T)(l,r,c)=(*snow->T)(l+1,r,c);
          (*snow->Dzl)(l,r,c)=(*snow->Dzl)(l+1,r,c);
        }
      initialize_snow(r, c, (*snow->lnum)(r,c), snow);
    }
  (*snow->lnum)(r,c)-=1;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void min_max_layer(long n, Vector<double>* Dmin, Vector<double>* Dmax,
                   Vector<double>* Dmin2, Vector<double>* Dmax2, long linf)
{

  long l,mup,mdw,N=Dmin->nh;

  if (n==N)
    {

      for (l=1; l<=n; l++)
        {
          Dmin2->co[l]=Dmin->co[l];
          Dmax2->co[l]=Dmax->co[l];
        }

    }
  else
    {

      mup=ceil(n/2.0);
      mdw=floor(n/2.0);

      if (linf<N-mup && linf<=mdw)
        {
          mdw=linf-1;
          mup=n-mdw;
        }
      else if (linf>mdw && linf>=N-mup)
        {
          mup=N-linf;
          mdw=n-mup;
        }

      for (l=1; l<=mdw; l++)
        {
          Dmin2->co[l]=Dmin->co[l];
          Dmax2->co[l]=Dmax->co[l];
        }

      for (l=n; l>n-mup; l--)
        {
          Dmin2->co[l]=Dmin->co[N+l-n];
          Dmax2->co[l]=Dmax->co[N+l-n];
        }

    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void initialize_snow(long r, long c, long l, STATEVAR_3D *snow)
{

  snow->w_ice->co[l][r][c]=0.0;
  snow->w_liq->co[l][r][c]=0.0;
  (*snow->Dzl)(l,r,c)=0.0;
  (*snow->T)(l,r,c)=0.0;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void show_Dminmax(long r, long c, double *Dmin, double *Dmax, long n)
{

  long l;

  printf("n:%ld\n",n);
  for (l=1; l<=n; l++)
    {
      printf("l:%ld ltot:%ld Dmin:%f Dmax:%f\n",l,n,Dmin[l],Dmax[l]);
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_snow_age(double Psnow, double Ts, double Dt, double Prestore,
                     double *tsnow_nondim)
{

  double r1, r2, r3;

  //effect of grain growth due to vapour diffusion
  r1=exp(5000.0*(1.0/tk-1.0/(Ts+tk)));

  //effect melt and refreezing*/
  r2=pow(r1,10);
  if (r2>1.0) r2=1.0;

  //effect of dirt
  r3=0.3;

  //non-dimensional snow age: 10 mm of snow precipitation restore snow age Dt(s)
  *tsnow_nondim=Fmax( 0.0,
                      (*tsnow_nondim+(r1+r2+r3)*Dt*1.0E-6)*(1.0-Psnow/Prestore) );
  if ((*tsnow_nondim)!=(*tsnow_nondim))
    printf("tsnow no value - tausn:%f P:%f Ts:%f r1:%f r2:%f r3:%f\n",
           *tsnow_nondim,Psnow,Ts,r1,r2,r3);

}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double snow_albedo(double ground_alb, double snowD, double AEP,
                   double freshsnow_alb, double C, double tsnow, double cosinc,
                   double ( *F)(double x))
{
  double A, Fage=1.0-1.0/(1.0+tsnow), w;
  A=freshsnow_alb*(1.0-C*Fage);
  A+=0.4*(1.0-A)*(*F)(cosinc);
  if (snowD<AEP)  //if snow is shallow (<AEP), interpolate between snow and ground albedo
    {
      w=(1.0-snowD/AEP)*exp(-snowD*0.5/AEP);
      A=w*ground_alb+(1.0-w)*A;
    }
  return (A);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Fzen(double cosinc)
{
  double f, b=2.0;
  if (cosinc<0.5)
    {
      f=(1.0/b)*((b+1.0)/(1.0+2.0*b*cosinc)-1.0);
    }
  else
    {
      f=0.0;
    }
  return (f);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double k_thermal_snow_Sturm(double density)   //W m^-1 K^-1 (Sturm, 1997)
{

  double kt;
  density*=0.001;
  if (density<0.156)
    {
      kt=0.023+0.234*density;
    }
  else
    {
      kt=0.138-1.01*density+3.233*pow(density,2.0);
    }
  return (kt);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double k_thermal_snow_Yen(double density)   //W m^-1 K^-1 (Yen, 1981)
{

  double kt;
  kt=k_ice*pow((density/rho_w),1.88);
  return (kt);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void non_dimensionalize_snowage(double *snowage, double Ta)
{

  double r1, r2, r3;

  r1=exp(5000.0*(1.0/273.16-1.0/(Ta+273.16)));
  r2=pow(r1,10);
  if (r2>1.0) r2=1.0;
  r3=0.3;

  *snowage*=((r1+r2+r3)*1.0E-6);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void glac2snow(double a, long r, long c, STATEVAR_3D *snow, STATEVAR_3D *glac)
{

  double h;
  long ns=(*snow->lnum)(r,c);

  h = internal_energy(glac->w_ice->co[1][r][c], glac->w_liq->co[1][r][c], (*glac->T)(1,r,c));
  if (ns>0) 
    h += internal_energy(snow->w_ice->co[1][r][c], snow->w_liq->co[1][r][c], (*snow->T)(1,r,c));

  (*snow->Dzl)(1,r,c) +=(*glac->Dzl)(1,r,c);
  snow->w_ice->co[1][r][c] += glac->w_ice->co[1][r][c];
  snow->w_liq->co[1][r][c] += glac->w_liq->co[1][r][c];
  if (ns>0) from_internal_energy(a, r, c, h, &(snow->w_ice->co[1][r][c]),
                                   &(snow->w_liq->co[1][r][c]), &((*snow->T)(1,r,c)));

 (*glac->Dzl)(1,r,c)=0.0;
  glac->w_liq->co[1][r][c]=0.0;
  glac->w_ice->co[1][r][c]=0.0;
  (*glac->T)(1,r,c)=0.0;
  (*glac->lnum)(r,c)=0;
  (*glac->type)(r,c)=0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void snow2glac(double a, long r, long c, STATEVAR_3D *snow, STATEVAR_3D *glac)
{

  double h;
  long ng=(*glac->lnum)(r,c);

  h = internal_energy(snow->w_ice->co[1][r][c], snow->w_liq->co[1][r][c],
                      (*snow->T)(1,r,c));
  if (ng>0) 
    h += internal_energy(glac->w_ice->co[ng][r][c], glac->w_liq->co[ng][r][c], (*glac->T)(ng,r,c));

  glac->Dzl->co[Fmaxlong(ng,1)][r][c] += (*snow->Dzl)(1,r,c);
  glac->w_ice->co[Fmaxlong(ng,1)][r][c] += snow->w_ice->co[1][r][c];
  glac->w_liq->co[Fmaxlong(ng,1)][r][c] += snow->w_liq->co[1][r][c];
  if (ng>0) 
    from_internal_energy(a, r, c, h, &(glac->w_ice->co[ng][r][c]), &(glac->w_liq->co[ng][r][c]), &((*glac->T)(ng,r,c)));

  (*snow->Dzl)(1,r,c)=0.0;
  snow->w_liq->co[1][r][c]=0.0;
  snow->w_ice->co[1][r][c]=0.0;
  (*snow->T)(1,r,c)=0.0;
  (*snow->lnum)(r,c)=0;
  (*snow->type)(r,c)=0;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void WBsnow(double Dt, long ns, long r, long c, STATEVAR_3D *snow,
            double *Melt, double *RainOnSnow, PAR *par, double slope, double Rain,
            ENERGY *E, double Evap)
{

  long l, m;
  double Se, th, thi, Wdt, Edt;

  Wdt = Rain;
  Edt = Evap;

  for (l=(*snow->lnum)(r,c); l>=1; l--)
    {

      if (l>ns)
        {

          snow->w_ice->co[l][r][c] = 0.;
          snow->w_liq->co[l][r][c] = 0.;
          (*snow->Dzl)(l,r,c) = 0.;

        }
      else
        {

          m = ns - l + 1;

          //assign
          if (Edt > E->ice->co[m] - E->deltaw->co[m])
            {

              Edt -= (E->ice->co[m] - E->deltaw->co[m]);
              Wdt += (E->liq->co[m] + E->deltaw->co[m]);
              snow->w_ice->co[l][r][c] = 0.;
              snow->w_liq->co[l][r][c] = 0.;
              (*snow->Dzl)(l,r,c) = 0;

            }
          else
            {

              (*snow->T)(l,r,c) = E->Temp->co[m];
              snow->w_ice->co[l][r][c] = Fmax(0., E->ice->co[m] - E->deltaw->co[m] - Edt);
              snow->w_liq->co[l][r][c] = Fmax(0., E->liq->co[m] + E->deltaw->co[m] + Wdt);
              (*snow->Dzl)(l,r,c) = 1.E3 * E->Dlayer->co[m];

              Edt = 0.;
              Wdt = 0.;

              if (snow->w_ice->co[l][r][c] > simpl_snow*par->max_weq_snow)
                {

                  //ACCOUNT FOR SNOW COMPACTION
                  //a)destructive metamorphism and overburden
                  snow_compactation(Dt, r, c, l, snow, slope, par);

                  //b)melting: snow depth decreases maintaining the same density
                  if (snow->w_ice->co[l][r][c]/E->ice->co[m] < 1) (*snow->Dzl)(l,r,c) *=
                      (snow->w_ice->co[l][r][c]/E->ice->co[m]);

                  //limit on max porosity
                  if (snow->w_ice->co[l][r][c] / (1.E-3*(*snow->Dzl)(l,r,c)*rho_w) >
                      par->snow_maxpor)
                    {
                      (*snow->Dzl)(l,r,c) = 1.E3 * snow->w_ice->co[l][r][c] /
                                               ( rho_w * par->snow_maxpor );
                    }

                  //CALCULATE LIQUID WATER GOING BELOW
                  th = snow->w_liq->co[l][r][c]/(1.0E-3*(*snow->Dzl)(l,r,c)*rho_w);    //[-]
                  thi = snow->w_ice->co[l][r][c]/(1.0E-3*(*snow->Dzl)(l,r,c)*rho_i);   //[-]
                  Se = (th - par->Sr*(1.0-thi))/( (1.0-thi) - par->Sr*(1.0-thi));
                  if (Se<0) Se=0.0;
                  if (Se>1) Se=1.0;
                  if (th>par->Sr*(1.0-thi)) Wdt += Fmin(5.0*pow(Se,3.0)*Dt,
                                                          (th - par->Sr*(1.0-thi)))*(*snow->Dzl)(l,r,c)*1.E-3*rho_w;
                  snow->w_liq->co[l][r][c] -= Wdt;

                }
              else
                {

                  Wdt += snow->w_liq->co[l][r][c];
                  snow->w_liq->co[l][r][c] = 0.;

                }

            }
        }
    }

  (*snow->lnum)(r,c) = ns;

  *Melt = Wdt - Rain;

  if (Wdt < Rain)
    {
      *RainOnSnow = Rain;
    }
  else
    {
      *RainOnSnow = 0.;
    }

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void new_snow(double a, long r, long c, STATEVAR_3D *snow, double P,
              double Dz, double T)
{

  long ns;
  double h;


  if ((*snow->type)(r,c)==0)
    {

      (*snow->Dzl)(1,r,c)+=Dz;
      snow->w_ice->co[1][r][c]+=P;

    }
  else
    {

      ns=(*snow->lnum)(r,c);

      h=internal_energy(snow->w_ice->co[ns][r][c], snow->w_liq->co[ns][r][c], (*snow->T)(ns,r,c));
      h+=(c_ice*P)*(Fmin(T, -0.1) - Tfreezing);

      (*snow->Dzl)(ns,r,c)+=Dz;
      snow->w_ice->co[ns][r][c]+=P;

      from_internal_energy(a, r, c, h, &(snow->w_ice->co[ns][r][c]),
                           &(snow->w_liq->co[ns][r][c]),&((*snow->T)(ns,r,c)));

    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void WBglacier(long ns, long ng, long r, long c, STATEVAR_3D *glac,
               double *Melt, PAR *par, ENERGY *E, double Evap)
{

  long l, m;
  double th, thi, Se, Edt;

  *Melt = 0;

  Edt = Evap;

  for (l=(*glac->lnum)(r,c); l>=1; l--)
    {

      if (l>ng)
        {

          glac->w_ice->co[l][r][c] = 0.;
          glac->w_liq->co[l][r][c] = 0.;
         (*glac->Dzl)(l,r,c) = 0.;

        }
      else
        {

          m = ns + ng - l + 1;

          //assign
          if (Edt > E->ice->co[m] - E->deltaw->co[m])
            {

              Edt -= (E->ice->co[m] - E->deltaw->co[m]);
              *Melt = *Melt +  (E->liq->co[m] + E->deltaw->co[m]);

              glac->w_ice->co[l][r][c] = 0.;
              glac->w_liq->co[l][r][c] = 0.;
             (*glac->Dzl)(l,r,c) = 0;

            }
          else
            {

              (*glac->T)(l,r,c) = E->Temp->co[m];
              glac->w_ice->co[l][r][c] = Fmax(0., E->ice->co[m] - E->deltaw->co[m] - Edt);
              glac->w_liq->co[l][r][c] = Fmax(0., E->liq->co[m] + E->deltaw->co[m]);
             (*glac->Dzl)(l,r,c) = 1.E3 * E->Dlayer->co[m];

              Edt = 0.;

              if (glac->w_ice->co[l][r][c] > simpl_snow*par->max_weq_glac)
                {

                  //COMPACTION
                  //melting: snow depth decreases maintaining the same density
                  if (glac->w_ice->co[l][r][c]/E->ice->co[m] < 1)(*glac->Dzl)(l,r,c) *=
                      (glac->w_ice->co[l][r][c]/E->ice->co[m]);

                  //limit on max porosity
                  if (glac->w_ice->co[l][r][c] / (1.E-3*(*glac->Dzl)(l,r,c)*rho_w) > 0.95)
                    {
                     (*glac->Dzl)(l,r,c) = 1.E3 * glac->w_ice->co[l][r][c] / ( rho_w * 0.95 );
                    }

                  //CALCULATE LIQUID WATER GOING BELOW
                  th = glac->w_liq->co[l][r][c]/(1.0E-3*(*glac->Dzl)(l,r,c)*rho_w);    //[-]
                  thi = glac->w_ice->co[l][r][c]/(1.0E-3*(*glac->Dzl)(l,r,c)*rho_i);   //[-]
                  Se = (th - par->Sr*(1.0-thi))/( (1.0-thi) - par->Sr*(1.0-thi));
                  if (Se<0) Se=0.0;
                  if (Se>1) Se=1.0;
                  if (th>par->Sr*(1.0-thi))
                    {
                      *Melt = *Melt + (th - par->Sr*(1.0-thi))*(*glac->Dzl)(l,r,c)*1.E-3*rho_w;
                      glac->w_liq->co[l][r][c] = par->Sr*(1.0-thi)
                                                 *(*glac->Dzl)(l,r,c)*1.E-3*rho_w;
                    }

                }
              else
                {

                  *Melt = *Melt + glac->w_liq->co[l][r][c];
                  glac->w_liq->co[l][r][c] = 0.;

                }
            }
        }
    }

  (*glac->lnum)(r,c) = ns;

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_SCA(STATEVAR_3D *snow, PAR *par, Matrix<double> *Z, double t)
{

  long l, r, c, cont=0, conttot=0;
  double T, D, SWE, Tmean=0.0, Tsmean=0.0, Dmean=0.0, SWEmean=0.0, SCA;

  double JD, JDfrom0;
  long day, month, year, hour, minute;

  FILE *f;

  char rec[ ]= {"_recNNNN"},crec[ ]= {"_crecNNNN"};
  char *name, *temp;

  JDfrom0 = convert_tfromstart_JDfrom0(t, (*par->init_date)(i_sim));
  convert_JDfrom0_JDandYear(JDfrom0, &JD, &year);
  convert_JDandYear_daymonthhourmin(JD, year, &day, &month, &hour, &minute);

  for (r=1; r<=Nr; r++)
    {
      for (c=1; c<=Nc; c++)
        {
          if ((long)(*Z)(r,c)!=number_novalue)
            {
              D=0.0;
              T=0.0;
              SWE=0.0;
              for (l=1; l<=(*snow->lnum)(r,c); l++)
                {
                  D+=(*snow->Dzl)(l,r,c);
                  SWE+=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c]);
                  T+= (*snow->T)(l,r,c) * (*snow->Dzl)(l,r,c);
                }
              if (D>0)T/=D;

              conttot++;
              if (D>10)
                {
                  cont++;
                  Dmean+=D;
                  SWEmean+=SWE;
                  Tmean+=T;
                  Tsmean+= (*snow->T)((*snow->lnum)(r,c),r,c);
                }

            }
        }
    }

  if (cont>0)
    {
      Dmean/=(double)cont;
      SWEmean/=(double)cont;
      Tmean/=(double)cont;
      Tsmean/=(double)cont;
      SCA=cont/(double)conttot;
    }
  else
    {
      SCA=0.0;
    }

  if (par->recover > 0) write_suffix(rec, par->recover, 4);
  if (par->n_ContRecovery > 0) write_suffix(crec, par->n_ContRecovery, 5);

  if (par->recover>0)
    {
      temp = join_strings(files[fSCA], rec);
      name = join_strings(temp, textfile);
      free(temp);
    }
  else if (par->n_ContRecovery>0)
    {
      temp = join_strings(files[fSCA], crec);
      name = join_strings(temp, textfile);
      free(temp);
    }
  else
    {
      name = join_strings(files[fSCA], textfile);
    }

  f=fopen(name,"a");
  fprintf(f,"%ld/%ld/%ld %ld:%02.0f",day,month,year,hour,(float)minute);
  fprintf(f,",%f,%f,%f",(JDfrom0-(*par->init_date)(i_sim))+(i_run-1)*
          ((*par->end_date)(i_sim)-(*par->init_date)(i_sim)),JDfrom0,JD);
  fprintf(f,",%f,%f,%f,%f,%f,%f\n",Dmean,SWEmean,Tmean,Tsmean,(1.0-SCA)*100.0,
          SCA*100.0);
  fclose(f);

  free(name);
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double max_dtheta_snow(double a, double b) {return (-pow(b/(3*a*a), 0.5));}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double theta_snow(double a, double b, double T)
{

  double th;

  th=1./(b+pow(a*T, 2.));
  if (T>0) th=1./(b+pow(a*0., 2.));
  return (th);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double dtheta_snow(double a, double b, double T)
{

  double dth;

  dth=-2.*a*a*T/pow(b+pow(a*T, 2.), 2.);
  if (T>0) dth=-2.*a*a*max_dtheta_snow(a,b)/pow(b+pow(a*max_dtheta_snow(a,b),
                                                        2.), 2.);
  return (dth);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void allocate_and_initialize_statevar_1D(STATEVAR_1D *V, double nan, long nl)
{

  V->Dzl.reset(new Vector<double>{nl});
  V->T.reset(new Vector<double>{nl});
  *(V->T) = nan;
  V->w_ice.reset(new Vector<double>{nl});
  V->w_liq.reset(new Vector<double>{nl});
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void deallocate_statevar_1D(STATEVAR_1D *V)
{

  free(V);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short copy_statevar_from3D_to1D(long r, long c, STATEVAR_3D *origin,
                                STATEVAR_1D *destination)
{

  long nl, l;

  nl = origin->Dzl->ndh;

  if (r<1 || r>origin->type->nrh) return 0;
  if (c<1 || c>origin->type->nch) return 0;

  if (nl != destination->Dzl->nh) return 0;

  destination->type = (*origin->type)(r,c);
  destination->lnum = (*origin->lnum)(r,c);
  for (l=1; l<=nl; l++)
    {
      (*destination->Dzl)(l) = (*origin->Dzl)(l,r,c);
      (*destination->T)(l) = (*origin->T)(l,r,c);
      (*destination->w_ice)(l) = origin->w_ice->co[l][r][c];
      (*destination->w_liq)(l) = origin->w_liq->co[l][r][c];
    }

  return 1;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double interpolate_snow(long r, long c, double h, long max, Tensor<double> *Dz,
                        Tensor<double> *Q, short k)
{

  double q, z, z0=0.;
  long l;
  short u;

  q = (double)number_novalue;

  if (h>0)  //downwards
    {
      u = -1;
      l = max+1;
    }
  else   //upwards
    {
      u = 1;
      l = 1;
      h = -h;
    }

  do
    {

      if (l == 1)
        {
          z = z0 + (*Dz)(l,r,c)/2.;
        }
      else if (l <= max)
        {
          z = z0 + (*Dz)(l,r,c)/2. + (*Dz)(l-1,r,c)/2.;
        }
      else
        {
          z = z0 + (*Dz)(max,r,c)/2.;
        }

      if (fabs(h) <= z && fabs(h) > z0)
        {
          if (l == 1)
            {
              if (k==0)
                {
                  q = (*Q)(l,r,c);
                }
              else
                {
                  q = (*Q)(l,r,c)/(*Dz)(l,r,c);
                }
            }
          else if (l <= max)
            {
              if (u>0)
                {
                  if (k==0)
                    {
                      q = ( (*Q)(l-1,r,c) * (h-z0) + (*Q)(l,r,c) * (z-h) ) / (z - z0);
                    }
                  else
                    {
                      q = ( (*Q)(l-1,r,c)/(*Dz)(l-1,r,c) * (h-z0) +
                            (*Q)(l,r,c)/(*Dz)(l,r,c) * (z-h) ) / (z - z0);
                    }
                }
              else
                {
                  if (k==0)
                    {
                      q = ( (*Q)(l-1,r,c) * (z-h) + (*Q)(l,r,c) * (h-z0) ) / (z - z0);
                    }
                  else
                    {
                      q = ( (*Q)(l-1,r,c)/(*Dz)(l-1,r,c) * (z-h) +
                            (*Q)(l,r,c)/(*Dz)(l,r,c) * (h-z0) ) / (z - z0);

                    }
                }
            }
          else
            {
              if (k==0)
                {
                  q = (*Q)(max,r,c);
                }
              else
                {
                  q = (*Q)(max,r,c)/(*Dz)(max,r,c);
                }

            }
        }

      z0 = z;

      l += u;

    }
  while ( (long)q == number_novalue && l <= max+1 && l >= 1);

  return q;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void copy_snowvar3D(STATEVAR_3D *from, STATEVAR_3D *to)
{

  long l, r, c;
  long nl=from->Dzl->ndh, nr=from->Dzl->nrh, nc=from->Dzl->nch;

  for (r=1; r<=nr; r++)
    {
      for (c=1; c<=nc; c++)
        {
          (*to->type)(r,c) = (*from->type)(r,c);
          (*to->lnum)(r,c) = (*from->lnum)(r,c);
          for (l=1; l<=nl; l++)
            {
              (*to->Dzl)(l,r,c) = (*from->Dzl)(l,r,c);
              to->w_liq->co[l][r][c] = from->w_liq->co[l][r][c];
              to->w_ice->co[l][r][c] = from->w_ice->co[l][r][c];
              (*to->T)(l,r,c) = (*from->T)(l,r,c);
            }
        }
    }
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


