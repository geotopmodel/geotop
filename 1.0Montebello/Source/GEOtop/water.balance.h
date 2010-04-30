
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */


    
void water_balance(ALLDATA *adt);

short Richards(double Dt, DOUBLETENSOR *P, double *loss, long *iter, ALLDATA *adt);

double J_d_water(long i, DOUBLEVECTOR *d, DOUBLEVECTOR *H, double Dt, ALLDATA *adt);

double J_d_water_(long i, DOUBLEVECTOR *d, DOUBLEVECTOR *H, double Dt, void *adt);

double F_water(long i, DOUBLEVECTOR *H, double Dt, ALLDATA *adt);

double cm_h(double cm0, double h, double h_thres);

void supflow(double Dt, double **h, double **dh, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double *Vout);

void find_dt_max(double Courant, double **h, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, double *dt);

//void routing(CHANNEL *cnet);

//int J_water_triplet(UMFPACK_REAL_TRIPLET *triplet, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

int J_water(DOUBLEVECTOR *Lx, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

int J_water_with_transposed(DOUBLEVECTOR *Lx, DOUBLEVECTOR *Ux, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

int find_matrix_K(DOUBLEVECTOR *Lx, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

int find_matrix_K_with_transposed(DOUBLEVECTOR *Lx, DOUBLEVECTOR *Ux, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

int find_dfdH(DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

int find_f(DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

void find_dt_max_channel(double Courant, DOUBLEVECTOR *h, TOPO *top, CHANNEL *cnet, PAR *par, double *dt);

void channel_flow(double Dt, DOUBLEVECTOR *h, DOUBLEVECTOR *dV, TOPO *top, CHANNEL *cnet, PAR *par, double *Vout);
