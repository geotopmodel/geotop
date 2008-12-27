/*
 * geometry_io.h
 *
 *  Created on: Nov 1, 2008
 *      Author: ecor
 */
/* This functions create polugons and lines if the line or vertex indices are written in a string */

LINE *get_line(DOUBLEVECTOR *vertex_x_coord, DOUBLEVECTOR *vertex_y_coord, long line_index, char *number_strings, short print);

LINEVECTOR *get_linevector(DOUBLEVECTOR *vertex_x_coord, DOUBLEVECTOR *vertex_y_coord,STRINGBIN *line_strings , short print);



POLYGON *get_polygon(long index,double x, double y, double z, char *edge_index_string,LINEVECTOR *all_lines, short print);

POLYGONVECTOR *get_polygonvector(DOUBLEVECTOR *centroid_x_coord, DOUBLEVECTOR *centroid_y_coord,STRINGBIN *polygon_strings ,LINEVECTOR *all_lines,short print);

/* end  This functions create polugons and lines if the line or vertex indices are written in a string  */

LINEVECTOR *extract_linvector_from_linevector(LONGVECTOR *nlines, LINEVECTOR *lines);


int fprint_linevector(char *filename, LINEVECTOR *lines);
int fprint_polygonvector(char *filename, POLYGONVECTOR *polygons);
