
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file additional_read_functions.c

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



//long index_pixel(long r, long c,long nrh, long nch);

//LONGMATRIX *indices(long nrh, long nch,long IBASE,long (*t_index)(long r, long c,long nrh, long nch));

LONGBIN *addresses(LONGMATRIX *lmask,long row_shift,long col_shift);

long index_pixel_from_a_bin(long r, long c, LONGVECTOR *s_index);

LONGMATRIX *m_indices_with_novalues(LONGBIN *laddresses, long nch, long novalue, long IBASE, long (*index_pixel_from_a_bin)(long r, long c, LONGVECTOR *s_index));

LONGMATRIX *m_indices_from_mask(DOUBLEMATRIX *mask,long row_shift,long col_shift,long novalue,long IBASE,long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index),DOUBLEVECTOR *V, int (*check_novalues)(double x, DOUBLEVECTOR *V));

POINT *new_point_from_raster(long r,long c, long nrh, long nch, double lx, double ly, double nsres, double ewres, double blc_x, double blc_y, long index);

LINE *new_horizontal_line_from_raster(long r,long c, long nrh, long nch, double nsres, double ewres, double blc_x, double blc_y, long iline, long ivertex1, long ivertex2);

LINE *new_vertical_line_from_raster(long r,long c, long nrh, long nch, double nsres, double ewres, double blc_x, double blc_y, long iline, long ivertex1, long ivertex2);



LINEVECTOR *get_linevector_from_raster_grid(LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,LONGMATRIX *i_vertex, double nsres, double ewres, double blc_x, double blc_y, long novalue);




POLYGON *new_pixel_from_raster(long index,long r, long c ,LINEVECTOR *lines, LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,double nsres, double ewres, double blc_x, double blc_y,short print);


POLYGONVECTOR *get_polygonvector_from_raster(LINEVECTOR *lines,LONGMATRIX *i_pixels,LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,double nsres, double ewres, double blc_x, double blc_y,long novalue, short print);
DOUBLEVECTOR *read_doublevector_from_raster(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref,LONGMATRIX *indices);

DOUBLEVECTOR *get_doublevector_from_doublematrix(LONGMATRIX *indices,DOUBLEMATRIX *M, double novalue);

DOUBLEVECTOR *read_doublevector_from_raster(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref,LONGMATRIX *indices);

int write_raster_from_doublevector(char *filename, DOUBLEVECTOR *v, T_INIT *UVref, LONGMATRIX *indices, DOUBLEMATRIX *Mref);

DOUBLEMATRIX *get_doublematrix_from_mapseries(LONGMATRIX *indices,DOUBLETENSOR *mapseries, T_INIT *UVref);
