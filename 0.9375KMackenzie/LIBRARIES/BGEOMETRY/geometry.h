

/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry.h

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
