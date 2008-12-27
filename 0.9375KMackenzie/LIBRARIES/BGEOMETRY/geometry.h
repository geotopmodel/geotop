/*
 * geometry.h
 *
 *  Created on: Oct 11, 2008
 *      Author: ecor
 */
#define MISSING_ELEMENT -99
#define TRUE 1
#define FALSE 0

/* the following lines define the attribute types
 * for a POINT, LINE and POLYGON
 *
 */

typedef struct {
	long index;
	double x;
	double y;
	double z;
	//attribute_point *attributes;

} POINT;

typedef struct {
	long index;
	POINT *begin;
	POINT *end;
	//long begin_point_index;
	//long end_point_index;
	double length2d;
	//attribute_line *attributes;
} LINE;


typedef struct {
	long index;
	double area2D;
	POINT *centroid;
	LONGVECTOR *edge_indices;
	//attribute_polygon *attributes;
} POLYGON;


typedef struct {
	short isdynamic;
	long nh,nl;
	LINE **element;
}  LINEVECTOR;



typedef struct {
	short isdynamic;
	long nh,nl;
	POLYGON **element;
}  POLYGONVECTOR;

/* array of attributes for each point,lines, polygons

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
	attribute_polygon **element;
}  attribute_polygon_array;*/

/* function header */

POINT *new_point(long index,double x,double y,double z);

LINE *new_line_from_points(long index, POINT *P1,POINT *P2);

POLYGON *new_polygon_from_a_linevector(LINEVECTOR *lines,POINT *centroid);

LINEVECTOR *new_linevector(long nh);

POLYGONVECTOR *new_polygonvector(long nh);
