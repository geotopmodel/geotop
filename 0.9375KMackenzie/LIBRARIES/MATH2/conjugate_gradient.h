
/* MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file conjugate_gradient.h

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


long conjugate_gradient_search(long icnt, double epsilon,  DOUBLEVECTOR *x, DOUBLEVECTOR *b, int (* funz)(DOUBLEVECTOR *y,DOUBLEVECTOR *x));
int linear_comb_doublevector(DOUBLEVECTOR *result,DOUBLEVECTOR *a, DOUBLEVECTOR *b, double ca, double cb);
long conjugate_gradient_search_LONG(long epsilon,LONGVECTOR *x,LONGVECTOR *b, int (* funz)(LONGVECTOR *y,LONGVECTOR *x));
double max_doublevector(DOUBLEVECTOR *v);
