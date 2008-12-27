/*
 * g_raster2plvector.h
 *
 *  Created on: Nov 25, 2008
 *      Author: ecor
 */
long index_pixel(long r, long c,long nrh, long nch);

LONGMATRIX *indices(long nrh, long nch,long IBASE,long (*t_index)(long r, long c,long nrh, long nch));

POINT *new_point_from_raster(long r,long c, long nrh, long nch, double lx, double ly, double nsres, double ewres, double blc_x, double blc_y, long index);

LINE *new_horizontal_line_from_raster(long r,long c, long nrh, long nch, double nsres, double ewres, double blc_x, double blc_y, long iline, long ivertex1, long ivertex2);

LINE *new_vertical_line_from_raster(long r,long c, long nrh, long nch, double nsres, double ewres, double blc_x, double blc_y, long iline, long ivertex1, long ivertex2);




LINEVECTOR *get_linevector_from_raster_grid(LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,LONGMATRIX *i_vertex, double nsres, double ewres, double blc_x, double blc_y);


POLYGON *new_pixel_from_raster(long index,long r, long c ,LINEVECTOR *lines, LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,double nsres, double ewres, double blc_x, double blc_y,short print);

POLYGONVECTOR *get_polygonvector_from_raster(LINEVECTOR *lines,LONGMATRIX *i_pixels,LONGMATRIX *i_horizontal,LONGMATRIX *i_vertical,double nsres, double ewres, double blc_x, double blc_y,short print);

DOUBLEVECTOR *read_doublevector_from_raster(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref,LONGMATRIX *indices);

DOUBLEVECTOR *get_doublevector_from_doublematrix(LONGMATRIX *indices,DOUBLEMATRIX *M, double novalue);

DOUBLEVECTOR *read_doublevector_from_raster(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref,LONGMATRIX *indices);

int write_raster_from_doublevector(char *filename, DOUBLEVECTOR *v, T_INIT *UVref, LONGMATRIX *indices, DOUBLEMATRIX *Mref);

DOUBLEMATRIX *get_doublematrix_from_mapseries(LONGMATRIX *indices,DOUBLETENSOR *mapseries, T_INIT *UVref);
