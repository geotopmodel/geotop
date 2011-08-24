
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

    
void water_balance(ALLDATA *adt);

short Richards(double Dt, double *loss, long *iter, ALLDATA *adt, FILE *f);

double cm_h(double cm0, double h, double h_thres1, double h_thres2);

int find_matrix_K(DOUBLEVECTOR *Lx, DOUBLEMATRIX *Klat, ALLDATA *adt, DOUBLEVECTOR *H, double Dt);

int find_dfdH(DOUBLEVECTOR *df, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *K, double Dt);

int find_f(DOUBLEVECTOR *f, ALLDATA *adt, DOUBLEVECTOR *H, DOUBLEMATRIX *K, double Dt);

void find_dt_Richards(short DDcomplex, double Courant, double **h, LAND *land, TOPO *top, SOIL *sl, PAR *par, double t, double *dt);

double find_hsup(double Psurface, double slope);

double find_dhsup(double slope);

double find_sup_pressure(double hsup, double slope);

double find_3Ddistance(double horizontal_distance, double vertical_distance);

void supflow(double Dt, double t, short DDland, short DDchannel, double **h, double **dh, double *hch, double *dhch, 
			 TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double *Voutnet, double *Voutland, FILE *flog);

void find_dt_max(short DDcomplex, double Courant, double **h, double *hch, LAND *land, TOPO *top, CHANNEL *cnet, PAR *par, double t, double *dt);

void channel_flow(double Dt, double t, short DDcomplex, double *h, double *dV, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double *Vout, FILE *f);

void find_dt_max_channel(short DDcomplex, double Courant, double *h, TOPO *top, CHANNEL *cnet, PAR *par, LAND *land, double t, double *dt);

void draining_land(double alpha, long r, long c, DOUBLEMATRIX *Z, double **h, DOUBLEMATRIX *LC, long *R, long *C);

void draining_channel(double alpha, long ch, DOUBLEMATRIX *Z, double *h, CHANNEL *cnet, long *CH);





