
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion KMackenzie

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 KMackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOtop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/



void water_balance(TOPO *top, SOIL *sl, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par, double time);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void vertical_water_balance (double Dt, DOUBLEMATRIX *Z, SOIL *sl, WATER *wat, double time, PAR *par);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void Richards(double Dt, long r, long c, SOIL *sl, DOUBLEVECTOR *psi, double Pnet, double t, PAR *par, double *masserrorcum);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void supflow(double Dt, double Dtmax, TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void subflow(double Dt, TOPO *top, SOIL *sl, PAR *par, WATER *wat, double *DeltaPsiMax, DOUBLETENSOR *P1);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void subflow_channel(double Dt, double Dtmax, CHANNEL *cnet, SOIL *sl, PAR *par);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void routing(CHANNEL *cnet);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void output_waterbalance(double Dt, WATER *wat, SOIL *sl, PAR *par, DOUBLEMATRIX *Z);

/*----------------------------------------------------------------------------------------------------------------------------------------------------*/
void set_psi(DOUBLETENSOR *P1, SOIL *sl, double *Q, double dt, long l, long r1, long c1, long r2, long c2, double psimin, double Esoil);
void set_psi_single(DOUBLETENSOR *psi, SOIL *sl, double *q, double dt, long l, long r, long c, double psimin, double Esoil);
