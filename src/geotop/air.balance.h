#ifndef _GEOTOP_AIR_BALANCE_H
#define _GEOTOP_AIR_BALANCE_H

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

short air_balance(double Dt, double JD0, double JD1, double JD2,
                    SOIL_STATE *L, SOIL_STATE *C,STATEVAR_3D *S, ALLDATA *adt, Vector<double> *Vsub,
                    Vector<double> *Vsup,
                    double *Voutnet, double *Voutlandsub, double *Voutlandsup,
                    double *Voutlandbottom);

short AirRichards3D(double Dt, SOIL_STATE *L, SOIL_STATE *C,STATEVAR_3D *S, ALLDATA *adt, double *loss, Vector<double> *Vsub,
                 double *Vbottom, double *Vlatsub, double *Total_Pnet, short updateK);


int find_matrix_Kair_3D(double Dt, SOIL_STATE *SL, SOIL_STATE *SC,STATEVAR_3D *S,
                     Vector<double> *Lx,Vector<double> *LxB,Vector<short> *FluxDir, ALLDATA *adt, Vector<double> *P1);

int find_matrix_LxHeat_air_3D(double Dt, SOIL_STATE *SL, SOIL_STATE *SC,
                     Vector<double> *LxJair,Vector<short> *FluxDir, Matrix<double> *Klat, Matrix<double> *Kbottom_l,
                     Vector<double> *Kbottom_ch, ALLDATA *adt, Vector<double> *H);
                     
int find_dfdPair_3D(double Dt, Vector<double> *df, ALLDATA *adt, SOIL_STATE *L,
                 SOIL_STATE *C,STATEVAR_3D *S, Vector<double> *P);
                 
int find_fair_3D(double Dt, Vector<double> *f, ALLDATA *adt, SOIL_STATE *L,
              SOIL_STATE *C,STATEVAR_3D *S, Vector<double> *P);


void find_dt_max(short DD, double Courant, RowView<double> &&h, LAND *land, TOPO *top,
                 CHANNEL *cnet, PAR *par, METEO *met, double t, double *dt);


#endif
