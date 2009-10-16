
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
    
    
    
void water_balance_1D(TOPO *top, SOIL *sl, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double time);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/												
void vertical_water_balance(double Dt, DOUBLEMATRIX *Z, SOIL *sl, WATER *wat, double time, PAR *par);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/													
void supflow(double Dt, double Dtmax, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/													
void subflow(double Dt, LAND *land, TOPO *top, SOIL *sl, PAR *par, WATER *wat, double *DeltaPsiMax, long *rref, long *cref, long *lref, 
	DOUBLETENSOR *P1, DOUBLETENSOR *P2);
		
/*----------------------------------------------------------------------------------------------------------------------------------------------------*/												
void subflow_channel(double Dt, double Dtmax, CHANNEL *cnet, SOIL *sl, PAR *par);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/												
void routing(CHANNEL *cnet);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/												
void output_waterbalance(double Dt, WATER *wat, SOIL *sl, PAR *par, DOUBLEMATRIX *Z);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/													
void check_q_1(DOUBLETENSOR *Pbeg, DOUBLETENSOR *Pend, double *Q, SOIL *sl, double dt, long l, long r, long c, double psimin, double Esoil, double max_dpsi);

void check_q_2(DOUBLETENSOR *Pbeg, DOUBLETENSOR *Pend, double *Q, SOIL *sl, double dt, long l, long r1, long c1, long r2, long c2, 
	double psimin, double Esoil);

void set_psi(DOUBLETENSOR *psi, DOUBLETENSOR *q, SOIL *sl, double dt, long l, long r, long c, double psimin, double Esoil);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/													
void Richards_1D(double Dt, long r, long c, SOIL *sl, DOUBLEVECTOR *psi, double *h, double Pnet, double t, PAR *par, double *masslosscum, 
	DOUBLETENSOR *q_sub);

void find_coeff_Richards_1D(long r, long c, double dt, short *bc, double Inf, double h, DOUBLEVECTOR *psit, DOUBLEVECTOR *e0, DOUBLEVECTOR *adi, 
	DOUBLEVECTOR *ad, DOUBLEVECTOR *ads, DOUBLEVECTOR *b, SOIL *sl, PAR *par, DOUBLETENSOR *q_sub);

void solve_Richards_1D(long r, long c, double dt, double *Inf, double h, DOUBLEVECTOR *psit, DOUBLEVECTOR *e1, double *massloss, SOIL *sl, 
	PAR *par, DOUBLETENSOR *q_sub);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/													