/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
BGEOMETRY Version 0.9375 KMackenzie

file bigcells2.h

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

typedef struct {
	LINEVECTOR *lines;

	POLYGONVECTOR *polygons;

	polygon_connection_attribute_array *links;

	char *file_resume_lines;
	char *file_resume_polygons;
	char *file_resume_connections;

	long boundary_indicator;

} GRID;

typedef struct {

	GRID * grid;

	LONGMATRIX *indices_pixel, *indices_horizontal_lines, *indices_vertical_lines, *indices_vertex;

	long novalue;

	long nhorizontal_lines;


} SQUARE_GRID;


typedef struct {

	//short is_Dynamic;

	long nl,nh;

	LONGMATRIX **element;

} LONGMATRIX_VECTOR;


typedef struct {

	SQUARE_GRID *big;
	SQUARE_GRID *fine;

	LONGMATRIX_VECTOR *small_content_polygon;
	LONGBIN   *small_content_line;

	char *file_resume_c_polygon;
	char *file_resume_c_line;


} DOUBLESQUARE_GRID ;





typedef struct {

	long reference_index_map;
	long nl,nh; /* nl=NL, nh: numbr of layers */

	DOUBLEMATRIX **layer;
	T_INIT *UV;
	int (*check_novalues)(double x, DOUBLEVECTOR *V); /* function which indentigies null values */

} RASTER_MAP;

typedef struct {

	RASTER_MAP *coarse;
	RASTER_MAP *fine;

	char *mapname_coarse;
	char *mapname_fine;

} DOUBLERASTER_MAP;

SQUARE_GRID *get_square_grid(DOUBLEMATRIX *DTM, T_INIT *UV,char *file_resume_lines, char *file_resume_polygons, char *file_resume_connections,long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index), int (*check_novalues)(double x, DOUBLEVECTOR *V),long displacement,short print);

void write_grid(GRID *grid);

void free_longmatrix_vector (LONGMATRIX_VECTOR *lmv);

RASTER_MAP *new_raster_map(long nh);
int add_map_to_raster(long index,char *filename, RASTER_MAP *raster, short a, short print);

DOUBLERASTER_MAP *get_doubleraster_map(long n_coarsemaps, long n_finemaps,char *coarse_mapname, char *fine_mapname,short print);

DOUBLESQUARE_GRID *get_doublesquare_grid(DOUBLERASTER_MAP *draster, char *resume_filenames,long (*index_pixel_from_a_bin_coarse)(long r, long c,LONGVECTOR *s_index),long (*index_pixel_from_a_bin_fine)(long r, long c,LONGVECTOR *s_index),long d_coarse,long d_fine,short print);

LONGMATRIX_VECTOR *new_longmatrix_vector(long nh);
LONGMATRIX_VECTOR *get_fine_indices(LONGMATRIX *fine,LONGMATRIX *coarse,long i_firstcell,long i_lastcell,long novalue,short print);




LONGBIN *get_fine_line_indices (LONGMATRIX *h_fine,LONGMATRIX *h_coarse,LONGMATRIX *v_fine,LONGMATRIX *v_coarse,long n_horizontal,long n_lines, long novalue,short print);

void free_longmatrix_vector (LONGMATRIX_VECTOR *lmv);

int write_longmatrix_vector(LONGMATRIX_VECTOR *lmv,char *filename);
int write_doublesquare_grid(DOUBLESQUARE_GRID *dsq);

void free_grid(GRID *grid);
void free_square_grid(SQUARE_GRID *sq);
void free_raster_map(RASTER_MAP *raster);
void free_doubleraster_map(DOUBLERASTER_MAP *draster);
void free_doublesquare_grid(DOUBLESQUARE_GRID *dsq);
char *strcpy2(char *dest, const char *src);
char *copy_string(char *origin);
