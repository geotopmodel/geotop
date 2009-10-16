
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
    
 	
    
void water_balance_3D(ALLDATA *adt, double time);

void Richards_3D(ALLDATA *adt);

void Assign_Psi(DOUBLEVECTOR *psi, SOIL *sl, LAND *land, PAR *par, TOPO *top);

double Solve_Richards_3D(long i, DOUBLEVECTOR *psi, ALLDATA *adt);

double Solve_Richards_3D_p(long i, DOUBLEVECTOR *psi, void *adt);

double Find_b(long i, DOUBLEVECTOR *theta0, ALLDATA *adt);

void supflow2(TOPO *top, LAND *land, WATER *wat, CHANNEL *cnet, PAR *par);

void routing2(CHANNEL *cnet);

void output_waterbalance2(WATER *wat, SOIL *sl, PAR *par, DOUBLEMATRIX *Z);


