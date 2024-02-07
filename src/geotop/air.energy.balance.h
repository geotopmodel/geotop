
#ifndef _GEOTOP_AIR_ENERGY_BALANCE_H
#define _GEOTOP_AIR_ENERGY_BALANCE_H

#include <cmath>
#include "math.optim.h"

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

short air_energy_balance(double Dt, double JD0, double JDb, double JDe,
                    SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, STATEVAR_3D *G, STATE_VEG *V,
                    Vector<double> *snowage,
                    ALLDATA *A, double *W);
                    
short Energy3D(double Dt, SOIL_STATE *L, SOIL_STATE *C, STATEVAR_3D *S, ALLDATA *adt, double *loss, short updateK);

int find_matrix_Kthermal_3D(double Dt, SOIL_STATE *SL, SOIL_STATE *SC,STATEVAR_3D *S,
                     Vector<double> *Lx, ALLDATA *adt, Vector<double> *T);
                     
int find_fthermal_3D(double Dt, Vector<double> *f, ALLDATA *adt, SOIL_STATE *L,
              SOIL_STATE *C, STATEVAR_3D *S,Vector<double> *T);
                     
double k_thermalcm2b(short snow, short a, double th_liq, double th_ice,
                 double th_sat, double k_solid);
                 
double TransferTestConstant(short HeatTransferModel, double Vel, double dp, double kthcm,
                 double Ht_a, double Ht_b, double Ht_n, double ThetaAir1,double printaux);

int find_dfdHthermal_3D(double Dt, Vector<double> *df, ALLDATA *adt, SOIL_STATE *L,
                 SOIL_STATE *C, Vector<double> *H, Matrix<double> *Klat);
                     
                 
#endif
