/*
 * geometry_attribute.h
 *
 *  Created on: Nov 7, 2008
 *      Author: ecor
 */


typedef void attribute_point;
typedef void attribute_line;
typedef void attribute_polygon;


typedef struct {

	LONGVECTOR *connections;
	DOUBLEVECTOR *d_connections;

} polygon_connection_attributes;


/* array of attributes for each point,lines, polygons */

typedef struct {
	short isdynamic;
	long nh,nl;
	attribute_point **element;
}  attribute_point_array;


typedef struct {
	short isdynamic;
	long nh,nl;
	attribute_line **element;
}  attribute_line_array;


typedef struct {
	short isdynamic;
	long nh,nl;
	polygon_connection_attributes **element;
}  polygon_connection_attribute_array;



/* function header */

polygon_connection_attributes *get_connection(POLYGON *polygon,POLYGONVECTOR *polygons, long boundary, short print);

polygon_connection_attribute_array *new_connection_attributes(long nh);


polygon_connection_attribute_array *get_connection_array(POLYGONVECTOR *polygons, long boundary, short print);

int fprint_polygonconnectionattributearray(char *filename,polygon_connection_attribute_array *pca);
