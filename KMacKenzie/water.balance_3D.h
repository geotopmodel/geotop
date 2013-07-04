
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie 

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie. 
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
    
 	
    
void water_balance_3D(ALLDATA *adt);

void Richards_3D(double Dt, DOUBLETENSOR *P, DOUBLEMATRIX *h, DOUBLEVECTOR *h_ch, double *loss, ALLDATA *adt);

double Solve_Richards_3D(long i, DOUBLEVECTOR *H1, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, ALLDATA *adt);

double Solve_Richards_3D_p(long i, DOUBLEVECTOR *H1, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, void *adt);

double Find_b(long i, DOUBLEVECTOR *H0, DOUBLEVECTOR *H00, double Dt, ALLDATA *adt);

void routing3(CHANNEL *cnet);

void output_waterbalance3(WATER *wat, SOIL *sl, PAR *par, DOUBLEMATRIX *Z);

void find_slope_H(long ***I, DOUBLEVECTOR *H, DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, DOUBLEMATRIX *slope);

double find_k_sup(long j1, long j2, LONGMATRIX *rc, DOUBLEMATRIX *slope, DOUBLEMATRIX *Z, DOUBLEVECTOR *H, double cm1, double cm2, double gamma);

void supflow3(double Dt, double Dtmax, DOUBLEMATRIX *h, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par);
