
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry_io.h

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



/* This functions create polygons and lines if the line or vertex indices are written in a string */

LINE *get_line(DOUBLEVECTOR *vertex_x_coord, DOUBLEVECTOR *vertex_y_coord, long line_index, char *number_strings, short print);

LINEVECTOR *get_linevector(DOUBLEVECTOR *vertex_x_coord, DOUBLEVECTOR *vertex_y_coord,STRINGBIN *line_strings , short print);



POLYGON *get_polygon(long index,double x, double y, double z, char *edge_index_string,LINEVECTOR *all_lines, short print);

POLYGONVECTOR *get_polygonvector(DOUBLEVECTOR *centroid_x_coord, DOUBLEVECTOR *centroid_y_coord,STRINGBIN *polygon_strings ,LINEVECTOR *all_lines,short print);

/* end  This functions create polugons and lines if the line or vertex indices are written in a string  */

LINEVECTOR *extract_linvector_from_linevector(LONGVECTOR *nlines, LINEVECTOR *lines);

int write_linevector(char *filename, LINEVECTOR *lines);
int fprint_linevector(char *filename, LINEVECTOR *lines);
int write_polygonvector(char *filename, POLYGONVECTOR *polygons);
int fprint_polygonvector(char *filename, POLYGONVECTOR *polygons);
