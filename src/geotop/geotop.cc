
/** STATEMENT:

    Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
    Geotop 3.0.0 - 31 Oct 2013
 
    Copyright (c), 2013 - Stefano Endrizzi
 
    This file is part of Geotop 2.0.0
 
    Geotop 2.0.0  is a free software and is distributed under GNU General Public License 
    v. 3.0 <http://www.gnu.org/licenses/>
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
    PARTICULAR PURPOSE
 
    Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of 
    developers and users that constructively interact.
    If you just use the code, please give feedback to the authors and the community.
    Any way you use the model, may be the most trivial one, is significantly helpful for the future 
    development of the Geotop model. Any feedback will be highly appreciated.
 
    If you have satisfactorily used the code, please acknowledge the authors.

*/

/*--------  1.  Include File, Prototype of the subroutine "time_loop", global variables  -------*/

#include <sys/stat.h>
#include "struct.geotop.h"
#include "input.h"
#include "output.h"
#include "times.h"
#include "constants.h"
#include "energy.balance.h"
#include "meteo.h"
#include "water.balance.h"
#include "snow.h"
#include "blowingsnow.h"
#include "tabs.h"
#include "deallocate.h"
#include "pedo.funct.h"
#include <string>
#include <iostream>
#include <version.h>
#include <logger.h>
#include <fstream>

void time_loop(ALLDATA *A);


/*----------   1. Global variables  ------------*/

long number_novalue;
long number_absent;
char *string_novalue;

std::unique_ptr<T_INIT> UV;

char *logfile = "/dev/null";
char **files;

long Nl,Nr,Nc;
double t_meteo, t_energy, t_water, t_sub, t_sup, t_blowingsnow, t_out;

double **odpnt, * *odp;
long *opnt, nopnt;
short *ipnt, *ibsn;
char **hpnt;

double *odbsn, *odb;
long *obsn, nobsn;
char **hbsn;

long *osnw, nosnw;
char **hsnw;

long *oglc, noglc;
char **hglc;

long *osl, nosl;
char **hsl;

FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav,
        *ffice, *fficeav, *ffsnowT, *ffsnowl, *ffsnowi, *ffsnowd, *ffglac;

long i_sim=0, i_run, i_sim0, i_run0;

char *SuccessfulRunFile, *FailedRunFile;

const char * WORKING_DIRECTORY;

time_t start_time;
double elapsed_time, elapsed_time_start, cum_time, max_time;

/*----------   2.  Begin of main and declaration of its variables (several structs)   ----------*/

int main(int argc,char *argv[])
{
  std::string wd;
  if (!argv[1])
  {
    std::cerr << "Wrong number of arguments. Abort.\n"
              << "Example of usage:\n"
              << "$ ./geotop /path/to/a/folder" << std::endl;
    exit(9);
  }
  else
  {
    wd = argv[1];
    if (wd.back() != '/')
      wd.append("/");
  }

  WORKING_DIRECTORY = wd.c_str();

  std::ofstream olog{wd+"geotop.log"};
  geolog.attach_file_stream(olog);
  geolog << "STATEMENT:\n\n"
         << "Geotop 3.0.0 - 2018\n\n"
         << "Geotop 3.0.0  is a free software and is distributed under GNU General Public License \
v. 3.0 <http://www.gnu.org/licenses/>\n"
         << "WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS \
FOR A PARTICULAR PURPOSE.\n" << std::endl;

  geolog << "WORKING DIRECTORY: " << wd << std::endl;
  geolog << "LOGFILE: " << wd << "geotop.log\n" << std::endl;

  std::unique_ptr<ALLDATA> adt;
  FILE* f;

  // assign novalues
  number_novalue = -9999;
  number_absent = -9998;
  string_novalue = assign_string("none");
  i_sim0 = 1;
  i_run0 = 1;
  cum_time = 0.;
  elapsed_time_start = 0.;

  /* dinamic allocations: */
  UV.reset(new T_INIT{});

  adt.reset(new ALLDATA);

  t_meteo=0.;
  t_energy=0.;
  t_water=0.;
  t_sub=0.;
  t_sup=0.;
  t_out=0.;
  t_blowingsnow=0.;


  /*------------------    3.  Acquisition of input data and initialisation    --------------------*/
  get_all_input(argc, argv, adt->T.get(), adt->S.get(), adt->L.get(), adt->M.get(), adt->W.get(),
                adt->C.get(), adt->P.get(), adt->E.get(), adt->N.get(), adt->G.get(), adt->I.get());

  /*-----------------   4. Time-loop for the balances of water-mass and egy   -----------------*/
  time_loop(adt.get());

  /*--------------------   5.Completion of the output files and deallocaions  --------------------*/
  dealloc_all(adt->T.get(), adt->S.get(), adt->L.get(), adt->W.get(), adt->C.get(), adt->P.get(),
              adt->E.get(), adt->N.get(), adt->G.get(), adt->M.get(), adt->I.get());


  geolog << "End of simulation!" << std::endl;

  f = fopen(SuccessfulRunFile, "w");
  fclose(f);
  free(SuccessfulRunFile);
  free(FailedRunFile);
  return 0;
}


/*----------------   6. The most important subroutine of the main: "time_loop"   ---------------*/
void time_loop(ALLDATA *A)
{
  GEOLOG_PREFIX(__func__);
  clock_t tstart, tend;
  short en=0, wt=0, out;
  long i, sy, r, c, j, l;
  double t, Dt, JD0, JDb, JDe, W, th, th0;
  double Vout, Voutsub, Voutsup, Vbottom, C0, C1;
  FILE* f;

  // double mean;

  std::unique_ptr<STATEVAR_3D> S, G;
  std::unique_ptr<SOIL_STATE> L, C;
  std::unique_ptr<STATE_VEG> V;
  std::unique_ptr<Vector<double>> a, Vsup_ch, Vsub_ch;


  S.reset(new STATEVAR_3D{(double)number_novalue,
                          A->P->max_snow_layers, Nr, Nc});

  if (A->P->max_glac_layers>0)
  {
    G.reset(new STATEVAR_3D{(double)number_novalue,
                            A->P->max_glac_layers, Nr, Nc});
  }
  L.reset(new SOIL_STATE{A->P->total_pixel, Nl});
  C.reset(new SOIL_STATE {A->C->r->nh, Nl});
  V.reset(new STATE_VEG);
  initialize_veg_state(V.get(), A->P->total_pixel);
  a.reset(new Vector<double>{A->P->total_pixel});
  Vsub_ch.reset(new Vector<double>{A->C->r->nh});
  Vsup_ch.reset(new Vector<double>{A->C->r->nh});

  time( &start_time );

  // periods
  i_sim = i_sim0;

  do
  {
    // runs
    A->I->time = A->P->delay_day_recover*86400.; // initialize time
    A->P->delay_day_recover = 0.;

    do
    {
      if ( A->I->time > (A->P->end_date->co[i_sim] - A->P->init_date->co[i_sim])
                        *86400. - 1.E-5)
      {
        // printf("Number of times the simulation #%ld has been run: %ld\n",i_sim,i_run);
        // f=fopen(logfile, "a");
        // fprintf(f,"Number of times the simulation #%ld has been run: %ld\n",i_sim,i_run);
        // fclose(f);

        print_run_average(A->S.get(), A->T.get(), A->P.get());

        i_run++;
        A->I->time = 0.0; // Initialize time

        A->M->line_interp_WEB_LR = 0;
        A->M->line_interp_Bsnow_LR = 0;

        for (i=1; i<=A->M->st->Z->nh; i++)
        {
          A->M->line_interp_WEB[i-1] = 0;
          A->M->line_interp_Bsnow[i-1] = 0;
        }

        if (i_run <= (*A->P->run_times)(i_sim))
        {
          reset_to_zero(A->P.get(), A->S.get(), A->L.get(), A->N.get(), A->G.get(),
                        A->E.get(), A->M.get(), A->W.get());
          init_run(A->S.get(), A->P.get());
        }
      }
      else
      {
        // find time step from file or inpts
        set_time_step(A->P.get(), A->I.get());

        // time at the beginning of the time step
        JD0 = A->P->init_date->co[i_sim]+A->I->time/secinday;

        // time step variables
        t = 0.;
        Dt = A->P->Dt;

        // time step subdivisions
        do
        {
          JDb = A->P->init_date->co[i_sim]+(A->I->time+t)/secinday;

          if (t + Dt > A->P->Dt)
            Dt = A->P->Dt - t;

          // iterations
          do
          {
            JDe = A->P->init_date->co[i_sim]+(A->I->time+t+Dt)/secinday;

            // copy state variables on
            copy_snowvar3D(A->N->S, S.get());
            *a = *(A->N->age);

            if (A->P->max_glac_layers>0) copy_snowvar3D(A->G->G, G.get());
            copy_soil_state(A->S->SS, L.get());
            copy_soil_state(A->C->SS, C.get());
            copy_veg_state(A->S->VS.get(), V.get());

            /* for (j=1; j<=A->W->H1->nh; j++) {
   l=A->T->lrc_cont->co[j][1];
   r=A->T->lrc_cont->co[j][2];
   c=A->T->lrc_cont->co[j][3];
   printf("START %ld %ld %ld %e\n",l,r,c,(*A->S->SS->P)(l,A->T->j_cont[r][c]));
   } */

            // init
            Vout = 0.;
            Voutsub = 0.;
            Voutsup = 0.;
            Vbottom = 0.;
            *Vsub_ch = 0.;
            *Vsup_ch = 0.;

            // meteo
            tstart=clock();
            meteo_distr(A->M->line_interp_WEB, A->M->line_interp_WEB_LR, A->M.get(),
                        A->W.get(), A->T.get(), A->P.get(), JD0, JDb, JDe);
            tend=clock();
            t_meteo+=(tend-tstart)/(double)CLOCKS_PER_SEC;

            if (A->P->en_balance == 1)
            {
              tstart=clock();
              en = EnergyBalance(Dt, JD0, JDb, JDe, L.get(), C.get(), S.get(), G.get(),
                                 V.get(), a.get(), A, &W);
              tend=clock();
              t_energy+=(tend-tstart)/(double)CLOCKS_PER_SEC;
            }

            if (A->P->wat_balance == 1 && en == 0)
            {
              tstart=clock();
              wt = water_balance(Dt, JD0, JDb, JDe, L.get(), C.get(), A, Vsub_ch.get(),
                                 Vsup_ch.get(), &Vout, &Voutsub, &Voutsup, &Vbottom);
              tend=clock();
              t_water+=(tend-tstart)/(double)CLOCKS_PER_SEC;
            }

            if (en != 0 || wt != 0)
            {
              if (Dt > A->P->min_Dt) Dt *= 0.5;
              out = 0;

              if (en != 0)
              {
                geolog << "Energy balance not converging" << std::endl;
              }
              else
              {
                geolog << "Water balance not converging" << std::endl;
              }
              geolog << "Reducing time step to "<< Dt << "s, t: "<< t <<" s" << std::endl;
            }
            else
            {
              out = 1;
            }
            // printf("Dt:%f min:%f\n",Dt,A->P->min_Dt);
          }
          while ( out == 0 && Dt > A->P->min_Dt );

          /* if (en != 0 || wt != 0) {
     f = fopen(FailedRunFile, "w");
     fprintf(f, "Simulation Period:%ld\n",i_sim);
     fprintf(f, "Run Time:%ld\n",i_run);
     fprintf(f, "Number of days after start:%f\n",A->I->time/86400.);

     if (en != 0 && wt == 0) {
     fprintf(f, "ERROR: Energy balance does not converge, Dt:%f\n",Dt);
     }else if (en == 0 && wt != 0) {
     fprintf(f, "ERROR: Water balance does not converge, Dt:%f\n",Dt);
     }else {
     fprintf(f, "ERROR: Water and energy balance do not converge, Dt:%f\n",Dt);
     }

     fclose(f);
     t_error("Fatal Error! Geotop is closed. See failing report.");
     } */

          if (en != 0 || wt != 0)
          {
            // f = fopen(FailedRunFile, "w");

            // fprintf(f, "Simulation Period:%ld\n",i_sim);
            // fprintf(f, "Run Time:%ld\n",i_run);
            // fprintf(f, "Number of days after start:%f\n",A->I->time/86400.);

            if (en != 0 && wt == 0)
            {
              geolog <<"WARNING: Energy balance does not converge, Dt: " << Dt<< std::endl;
            }
            else if (en == 0 && wt != 0)
            {
              geolog <<"WARNING: Water balance does not converge, Dt: " << Dt<< std::endl;
            }
            else
            {
              geolog <<"WARNING: Water and energy balance do not converge, Dt: " << Dt
                     << std::endl;
            }

            // t_error("Fatal Error! Geotop is closed. See failing report.");
          }

          t += Dt;

          if (A->P->state_pixel == 1 && A->P->dUzrun == 1)
          {
            for (j=1; j<=A->P->rc->nrh; j++)
            {
              for (l=1; l<=Nl; l++)
              {
                r = (*A->P->rc)(j,1);  // r (*A->P->rc)(j,1)
                c = (*A->P->rc)(j,2);
                sy = (*A->S->type)(r,c);

                th = theta_from_psi((*A->S->SS->P)(l,A->T->j_cont[r][c]),
                                    (*A->S->SS->thi)(l,A->T->j_cont[r][c]), l,
                                    A->S->pa->co[sy], PsiMin);

                if (th > A->S->pa->co[sy][jsat][l] - (*A->S->SS->thi)(l,A->T->j_cont[r][c]))
                  th = A->S->pa->co[sy][jsat][l] - (*A->S->SS->thi)(l,A->T->j_cont[r][c]);

                C0 = A->S->pa->co[sy][jct][l] * (1. - A->S->pa->co[sy][jsat][l])
                     * A->S->pa->co[sy][jdz][l] + c_ice*(*A->S->SS->thi)(l,A->T->j_cont[r][c]) + c_liq*th;

                th0 = th;

                th = theta_from_psi((*L->P)(l,A->T->j_cont[r][c]),
                                    (*L->thi)(l,A->T->j_cont[r][c]), l,
                                    A->S->pa->co[sy], PsiMin);

                if (th > A->S->pa->co[sy][jsat][l]-(*L->thi)(l,A->T->j_cont[r][c]))
                  th = A->S->pa->co[sy][jsat][l]-(*L->thi)(l,A->T->j_cont[r][c]);

                C1 = A->S->pa->co[sy][jct][l] * (1. - A->S->pa->co[sy][jsat][l])
                     * A->S->pa->co[sy][jdz][l] + c_ice*(*L->thi)(l,A->T->j_cont[r][c]) + c_liq*th;

                (*A->S->dUzrun)(j,l) += 1.E-6
                                          * (0.5*(C0+C1) * ((*L->T)(l,A->T->j_cont[r][c]) - (*A->S->SS->T)(l,A->T->j_cont[r][c]))
                                             + Lf *(th-th0) *A->S->pa->co[sy][jdz][l] );
              }
            }
          }

          // write state variables
          copy_snowvar3D(S.get(), A->N->S);
          *(A->N->age) = *a;

          if (A->P->max_glac_layers>0) copy_snowvar3D(G.get(), A->G->G);
          copy_soil_state(L.get(), A->S->SS);
          copy_soil_state(C.get(), A->C->SS);
          copy_veg_state(V.get(), A->S->VS.get());
          *(A->C->Vsub) += *Vsub_ch;
          *(A->C->Vsup) += *Vsup_ch;
          A->C->Vout += Vout;
          A->W->Voutbottom += Vbottom;
          A->W->Voutlandsub += Voutsub;
          A->W->Voutlandsup += Voutsup;

          // printf("%f\n",A->I->time);

          // record time step
          odb[ootimestep] = Dt * (Dt/A->P->Dtplot_basin->co[i_sim]);

          // write output variables
          fill_output_vectors(Dt, W, A->E.get(), A->N.get(), A->G.get(), A->W.get(),
                              A->M.get(), A->P.get(), A->I.get(), A->T.get(), A->S.get());

          // reset Dt
          if (Dt < A->P->Dt) Dt *= 2.;

        }
        while (t < A->P->Dt);

        if (A->P->blowing_snow==1)
        {
          tstart=clock();
          windtrans_snow(A->N.get(), A->M.get(), A->W.get(), A->L.get(), A->T.get(),
                         A->P.get(), A->I->time);
          tend=clock();
          t_blowingsnow+=(tend-tstart)/(double)CLOCKS_PER_SEC;
        }

        tstart=clock();
        write_output(A->I.get(), A->W.get(), A->C.get(), A->P.get(), A->T.get(), A->L.get(),
                     A->S.get(), A->E.get(), A->N.get(), A->G.get(), A->M.get());
        tend=clock();
        t_out+=(tend-tstart)/(double)CLOCKS_PER_SEC;

        A->I->time += A->P->Dt; // Increase TIME
      }
    }
    while (i_run <= (*A->P->run_times)(i_sim)); // end of time-cycle

    if (A->P->newperiodinit != 0)
      end_period_1D(A->S.get(), A->T.get(), A->P.get());
    if (i_sim < A->P->init_date->nh)
      change_grid(i_sim, i_sim+1, A->P.get(), A->T.get(), A->L.get(), A->W.get(), A->C.get());

    reset_to_zero(A->P.get(), A->S.get(), A->L.get(), A->N.get(), A->G.get(), A->E.get(),
                  A->M.get(), A->W.get());
    init_run(A->S.get(), A->P.get());

    i_sim++;
    i_run0 = 1;
    i_run = i_run0;

  }
  while (i_sim <= A->P->init_date->nh);
}
/*--------------------------------------------------------------------------------------------*/

