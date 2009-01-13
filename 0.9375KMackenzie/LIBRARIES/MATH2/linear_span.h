
/* MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file linear_span.h

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico and Riccardo Rigon

This file is part of MATH2.
 MATH2 is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MATH2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

int no_value_function(double x,DOUBLEVECTOR *V);

DOUBLEMATRIX *extract_a_new_map(DOUBLETENSOR *xtensor, long l,T_INIT *UVref);

DOUBLEVECTOR *prod_doublematvet(DOUBLEMATRIX *m, DOUBLEVECTOR *v);

double prodscal (DOUBLEVECTOR *a, DOUBLEVECTOR *b);

DOUBLEVECTOR *scalxvet (double a, DOUBLEVECTOR *b);

DOUBLETENSOR *linear_span_doubletensor(double c1, double c2, DOUBLETENSOR *T1,DOUBLETENSOR *T2, DOUBLEVECTOR *V);

DOUBLEMATRIX *linear_span_doublematrix(double c1, double c2, DOUBLEMATRIX *M1a,DOUBLEMATRIX *M2a, DOUBLEVECTOR *V);

DOUBLEMATRIX *transpose_doublematrix(DOUBLEMATRIX *M);

DOUBLEVECTOR *extract_a_column_from_doublematrix(long d,DOUBLEMATRIX *M);

DOUBLEVECTOR *extract_a_vertical_column_from_doubletensor(long r,long c, DOUBLETENSOR *T);
