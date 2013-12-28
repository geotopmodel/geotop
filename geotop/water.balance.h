
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#ifndef WATER_BALANCE_H
#define WATER_BALANCE_H

#include "constants.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
//#include "../libraries/math/sparse_matrix.h"
//#include "../libraries/math/util_math.h"
#include "meteodata.h"

#include <time.h>

extern long number_novalue, number_absent;

//extern T_INIT *UV;
extern TInit *UV;
extern std::vector<std::string> files;
extern std::string logfile;

extern long Nl, Nr, Nc;
extern double *odb;
extern double t_sub, t_sup;

//extern char *FailedRunFile;
extern std::string FailedRunFile;

extern long i_sim, i_run;

//subsurface flow constants
#define tol_max_GC 1.E+5
#define tol_min_GC 1.E-13
#define max_res_adm 1.E-2
#define MM 1
#define ni_wat 1.E-7
#define maxITER_rec_K 10


//short water_balance(double Dt, double JD0, double JD1, double JD2, SOIL_STATE *L, SOIL_STATE *C, ALLDATA *adt, DOUBLEVECTOR *Vsub, DOUBLEVECTOR *Vsup,
//					double *Voutnet, double *Voutlandsub, double *Voutlandsup, double *Voutlandbottom);

short water_balance(double Dt, double JD0, double JD1, double JD2, SoilState *L, SoilState *C, AllData *adt, GeoVector<double>& Vsub, GeoVector<double>& Vsup,
					double *Voutnet, double *Voutlandsub, double *Voutlandsup, double *Voutlandbottom);


//short Richards3D(double Dt, SOIL_STATE *L, SOIL_STATE *C, ALLDATA *adt, FILE *flog, double *loss, DOUBLEVECTOR *Vsub, double *Vbottom, double *Vlatsub, double *Total_Pnet, short updateK);
short Richards3D(double Dt, SoilState *L, SoilState *C, AllData *adt, FILE *flog, double *loss, GeoVector<double>& Vsub, double *Vbottom, double *Vlatsub, double *Total_Pnet, short updateK);

//short Richards1D(long c, double Dt, SOIL_STATE *L, ALLDATA *adt, FILE *flog, double *loss, double *Vbottom, double *Vlat, double *Total_Pnet, short updateK);
short Richards1D(long c, double Dt, SoilState *L, AllData *adt, FILE *flog, double *loss, double *Vbottom, double *Vlat, double *Total_Pnet, short updateK);

double cm_h(double cm0, double h, double h_thres1, double h_thres2);

//int find_matrix_K_3D(double Dt, SOIL_STATE *SL, SOIL_STATE *SC, DOUBLEVECTOR *Lx, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom_l, DOUBLEVECTOR *Kbottom_ch, ALLDATA *adt, DOUBLEVECTOR *H);
int find_matrix_K_3D(double Dt, SoilState *SL, SoilState *SC, GeoVector<double>& Lx, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom_l, GeoVector<double>& Kbottom_ch, AllData *adt, const GeoVector<double>& H);

//int find_matrix_K_1D(long c, double Dt, SOIL_STATE *L, DOUBLEVECTOR *Lx, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom, ALLDATA *adt, DOUBLEVECTOR *H);

int find_matrix_K_1D(long c, double Dt, SoilState *L, GeoVector<double>& Lx, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom, AllData *adt, const GeoVector<double>& H);

//int find_dfdH_3D(double Dt, DOUBLEVECTOR *df, ALLDATA *adt, SOIL_STATE *L, SOIL_STATE *C, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat);
  int find_dfdH_3D(double Dt, GeoVector<double>& df, AllData *adt, SoilState *L, SoilState *C, const GeoVector<double>& H, GeoMatrix<double>& Klat);

//int find_dfdH_1D(long c, double Dt, SOIL_STATE *L, DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat);
  int find_dfdH_1D(long c, double Dt, SoilState *L, GeoVector<double>& df, AllData *adt, const GeoVector<double>& H, GeoMatrix<double>& Klat);

//int find_f_3D(double Dt, DOUBLEVECTOR *f, ALLDATA *adt, SOIL_STATE *L, SOIL_STATE *C, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom_l, DOUBLEVECTOR *Kbottom_ch);
int find_f_3D(double Dt, GeoVector<double>& f, AllData *adt, SoilState *L, SoilState *C, const GeoVector<double>& H, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom_l, const GeoVector<double>& Kbottom_ch);

//int find_f_1D(long c, double Dt, SoilState *L, DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *Klat, DOUBLEMATRIX *Kbottom);
  int find_f_1D(long c, double Dt, SoilState *L, GeoVector<double>& f, AllData *adt, const GeoVector<double>& H, GeoMatrix<double>& Klat, GeoMatrix<double>& Kbottom);

double find_3Ddistance(double horizontal_distance, double vertical_distance);

//void find_dt_max(double Courant, double *h, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, METEO *met, double t, double *dt);
  void find_dt_max(double Courant, GeoMatrix<double>& h, Land *land, Topo *top, Channel *cnet, Par *par, Meteo *met, double t, double *dt);

//void find_dt_max_chla(double Courant, double *h, double *hch, TOPO *top, CHANNEL *cnet, PAR *par, double t, double *dt);
  void find_dt_max_chla(double Courant, GeoMatrix<double>& h, GeoMatrix<double>& hch, Topo *top, Channel *cnet, Par *par, double t, double *dt);

//void supflow(double Dt, double t, double *h, double *dV, double *hch, double *dhch, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet,
//			 PAR *par, METEO *met, DOUBLEVECTOR *Vsup, double *Voutnet, double *Voutland, FILE *flog);
  void supflow(double Dt, double t, GeoMatrix<double>& h, double *dV, GeoMatrix<double>& hch, double *dhch, Topo *top, Land *land, Water *wat, Channel *cnet,Par *par, Meteo *met, GeoVector<double>& Vsup, double *Voutnet, double *Voutland, FILE *flog, double *mm1, double *mm2, double *mmo);

//void supflow_chla(double Dt, double t, double *h, double *hch, TOPO *top, WATER *wat, CHANNEL *cnet, PAR *par, DOUBLEVECTOR *Vsup, FILE *flog, long *cnt);
  void supflow_chla(double Dt, double t, GeoMatrix<double>& h, GeoMatrix<double>& hch, Topo *top, Water *wat, Channel *cnet, Par *par, GeoVector<double>& Vsup, FILE *flog, long *cnt);

//void channel_flow(double Dt, double t, short DDcomplex, double *h, double *dV, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double *Vout, FILE *f, long *cnt);
  void channel_flow(double Dt, double t, short DDcomplex, GeoMatrix<double>& h, double *dV, Topo *top, Channel *cnet, Par *par, Land *land, double *Vout, FILE *f, long *cnt);

//void find_dt_max_channel(short DDcomplex, double Courant, double *h, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double t, double *dt);
  void find_dt_max_channel(short DDcomplex, double Courant, GeoMatrix<double>& h, Topo *top, Channel *cnet, Par *par, Land *land, double t, double *dt);

//void draining_land(double alpha, long i, TOPO *T, LAND *L, PAR *P, double *h, long *I, double *Q);
  void draining_land(double alpha, long i, Topo *T, Land *L, Par *P, GeoMatrix<double>& h, GeoMatrix<long>& I, GeoMatrix<double>& Q, long row);

//void draining_channel(double alpha, long ch, DOUBLEMATRIX *Z, double *h, CHANNEL *cnet, long *CH);
  void draining_channel(double alpha, long ch, GeoMatrix<double>& Z, GeoMatrix<double>& h, Channel *cnet, long *CH);


#endif

