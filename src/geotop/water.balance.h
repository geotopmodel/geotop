
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

short water_balance(double Dt, double JD0, double JD1, double JD2,
                    SOIL_STATE *L, SOIL_STATE *C, ALLDATA *adt, Vector<double> *Vsub,
                    Vector<double> *Vsup,
                    double *Voutnet, double *Voutlandsub, double *Voutlandsup,
                    double *Voutlandbottom);

short Richards3D(double Dt, SOIL_STATE *L, SOIL_STATE *C, ALLDATA *adt, double *loss, Vector<double> *Vsub,
                 double *Vbottom, double *Vlatsub, double *Total_Pnet, short updateK);

short Richards1D(long c, double Dt, SOIL_STATE *L, ALLDATA *adt, double *loss, double *Vbottom, double *Vlat,
                 double *Total_Pnet, short updateK);

double cm_h(double cm0, double h, double h_thres1, double h_thres2);

int find_matrix_K_3D(double Dt, SOIL_STATE *SL, SOIL_STATE *SC,
                     Vector<double> *Lx, Matrix<double> *Klat, Matrix<double> *Kbottom_l,
                     Vector<double> *Kbottom_ch, ALLDATA *adt, Vector<double> *H);

int find_matrix_K_1D(long c, double Dt, SOIL_STATE *L, Vector<double> *Lx,
                     Matrix<double> *Klat, Matrix<double> *Kbottom, ALLDATA *adt, Vector<double> *H);

int find_dfdH_3D(double Dt, Vector<double> *df, ALLDATA *adt, SOIL_STATE *L,
                 SOIL_STATE *C, Vector<double> *H, Matrix<double> *Klat);

int find_dfdH_1D(long c, double Dt, SOIL_STATE *L, Vector<double> *df,
                 ALLDATA *adt, Vector<double> *H, Matrix<double> *Klat);

int find_f_3D(double Dt, Vector<double> *f, ALLDATA *adt, SOIL_STATE *L,
              SOIL_STATE *C, Vector<double> *H, Matrix<double> *Klat, Matrix<double> *Kbottom_l,
              Vector<double> *Kbottom_ch);

int find_f_1D(long c, double Dt, SOIL_STATE *L, Vector<double> *f, ALLDATA *adt,
              Vector<double> *H, Matrix<double> *Klat, Matrix<double> *Kbottom);

double find_3Ddistance(double horizontal_distance, double vertical_distance);

void find_dt_max(short DD, double Courant, MatrixRow<double> &&h, LAND *land, TOPO *top,
                 CHANNEL *cnet, PAR *par, METEO *met, double t, double *dt);

void find_dt_max_chla(double Courant, MatrixRow<double> &&h, MatrixRow<double> &&hch, TOPO *top,
                      CHANNEL *cnet, PAR *par, double t, double *dt);

void supflow(short DDland, short DDch, double Dt, double t, MatrixRow<double> &&h, double *dV, MatrixRow<double> &&hch,
             double *dhch, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, METEO *met,
             Vector<double> *Vsup, double *Voutnet, double *Voutland, double *mm1, double *mm2, double *mmo);

void supflow_chla(double Dt, double t, MatrixRow<double> &&h, MatrixRow<double> &&hch, TOPO *top, WATER *wat,
                  CHANNEL *cnet, PAR *par,
                  Vector<double> *Vsup, long *cnt);

void channel_flow(double Dt, double t, short DDcomplex, MatrixRow<double> &&h, double *dV, TOPO *top, CHANNEL *cnet,
                  PAR *par,
                  LAND *land, double *Vout, long *cnt);

void find_dt_max_channel(short DDcomplex, double Courant, MatrixRow<double> &&h,
                         TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double t, double *dt);

void draining_land(double alpha, long i, TOPO *T, LAND *L, PAR *P,
                   CHANNEL *cnet, MatrixRow<double> &&h, MatrixRow<long> &&I, MatrixRow<double> &&Q);

void draining_channel(double alpha, long ch, Matrix<double> *Z, MatrixRow<double> &&h,
                      CHANNEL *cnet, long *CH);





