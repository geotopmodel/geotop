/* STATEMENT:

C-FORTRAN libraries

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon

 LICENSE:

 This file is part of C-FORTRAN libraries. 
 C-FORTRAN libraries is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
    
    
    
//float *c2f(DOUBLEMATRIX *M);
//void f2c(float *V, DOUBLEMATRIX *M);
int im(int ix, int iy, int nx, int ny);
char *allocate_c(int n);
float *allocatev_f(int n);
float *allocatem_f(int n, int m);
float *allocatet_f(int n, int m, int k);
int *allocatev_int(int n);
int *allocatem_int(int n, int k);
int *allocatet_int(int n, int k, int q);
