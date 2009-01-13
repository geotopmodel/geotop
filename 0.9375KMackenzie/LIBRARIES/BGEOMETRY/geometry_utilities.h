
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry_utilities.h

Copyright, 2009 Emanuele Cordano and Riccardo Rigon

This file is part of BGEOMETRY.
 BGEOMETRY is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BGEOMETRY is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

int signum (long a);

double area2d(POINT *P1,POINT *P2, POINT *P3);

int check_3_numbers(long a, long b, long c);

int check_vertex_intersection_2_lines(long b1,long e1,long b2,long e2);

long query_freq_longvector(LONGVECTOR *lvector,long value);

int segment_intersection_occurence(LINE *L1,LINE *L2);

int signum_double (double a);

double distance2D(POINT *P1,POINT *P2);

long shared_edges(POLYGON *PO1, POLYGON *PO2, long no_intersection, long *l_po1, long *l_po2, double *dist_centroids);
