
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry.c

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

#include "turtle.h"
#include "geometry.h"
#include "geometry_utilities.h"
#include "geometry_freememory.h"

POINT *new_point(long index,double x,double y,double z){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date October 2008
	 *
	 *\param index - (long) point index
	 *\param x - (double) x coordinate
	 *\param y - (doulbe) y coordinate
	 *\param z - (double) z coordinate
	 *
	 *\return cretes a new point with fixed coordinates and index (Point Attributes are not allocated!!)
	 *
	 */
	POINT *P;

	P=(POINT *)malloc((sizeof(POINT)));
	if (!P) t_error("point was not allocated");
    P->index=index;
	P->x=x;
	P->y=y;
	P->z=z;
	//P->attributes=(attribute_point *)malloc((size_t)(sizeof(attribute_point)));
	//if (!P) t_error("attribute_point was not allocated");
	//P->attributes=attributes;

  return P;


}



LINE *new_line_from_points(long index, POINT *P1,POINT *P2){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date October 2008
	 *
	 * \parem index - (long) line index
	 * \param P1 - (POINT *) begin point of the line
	 * \param P2 - (POINT *) end   point of the line
	 *
	 * \return creates a line between two point
	 * \warning Once allocated the line, the two points P1 and P2 are dellacated!!!
	 */
	LINE *L;
	L=(LINE *)malloc((sizeof(LINE)));
	if (!L) t_error("line was not allocated");
	L->index=index;

	L->begin=new_point(P1->index,P1->x,P1->y,P1->z);
	L->end=new_point(P2->index,P2->x,P2->y,P2->z);
//	L->begin_point_index=P1->index;
//	L->end_point_index=P2->index;
//	L->attributes=(attribute_line *)malloc(sizeof(attribute_line));
//	if (!L->attributes) t_error("attribute_line was not allocated");
//	L->attributes=attributes;
	L->length2d=distance2D(P1,P2);


	return L;

}

POLYGON *new_polygon_from_a_linevector(LINEVECTOR *lines,POINT *centroid){
	/*
	 * \author Emanuele Cordano, Davide Giacomelli
	 * \date October 2008
	 *
	 *\param lines - (LINEVECTOR *) a vector of lines
	 *\param centrod - (POINT *) the centroid of the polygon
	 *
	 *
	 *\return creates a polygon with edges taken by the line vector (Polygon Attributes are not allocated!!)
	 */
	POLYGON *PO;
	LONGVECTOR *lvertices;
	long i,j,l,s,se,a;
	double area;

	PO=(POLYGON *)malloc(sizeof(POLYGON));
	if (!PO) t_error("polygon was not allocated");

	PO->centroid=new_point(centroid->index,centroid->x,centroid->y,centroid->z);

	PO->index=centroid->index;





	for (i=lines->nl;i<=lines->nh;i++){
		for (j=lines->nl;j<i;j++){
			 if (segment_intersection_occurence(lines->element[i],lines->element[j])!=0) {
				 printf("Warning:the lines %ld and %ld (polyogon-internal numeration) share an internal point, the polygon %ld cannot be created!! \n",i,j,centroid->index);
				 return PO;
			 }
		}
	}

	lvertices=new_longvector(lines->nh*2);


    for (i=lines->nl;i<=lines->nh;i++){
    	lvertices->element[2*i-1]=lines->element[i]->begin->index;
    	lvertices->element[2*i]=lines->element[i]->end->index;
    	if (lines->element[i]->end->index==lines->element[i]->begin->index) {
    		printf("Warning:the line %ld has the same extremes in index %ld, the polygon %ld cannot be created!!",i,lines->element[i]->begin->index,centroid->index);
    	    return PO;
    	}
    }

  //  print_longvector_elements(lvertices,1);

    for (i=lines->nl;i<=lines->nh;i++){
    	if(query_freq_longvector(lvertices,lines->element[i]->begin->index)!=2) {
    		printf("Warning: problem in begin vertex of %ld, (%ld points in the same vertex) the polygon %ld cannot be created!! \n",lines->element[i]->begin->index,query_freq_longvector(lvertices,lines->element[i]->begin->index),centroid->index);

    		return PO;
    	}

    	if(query_freq_longvector(lvertices,lines->element[i]->end->index)!=2) {
    		printf("Warning: problem in end vertex of line %ld, (%ld points in the same vertex) the polygon %ld cannot be created!! \n",lines->element[i]->end->index,query_freq_longvector(lvertices,lines->element[i]->begin->index),centroid->index);

    	    return PO;
    	 }
    }

/* The polygon can be defined */
/* It calculates the area of a polygon */
/* It defines the edge of the polygon */
    area=0.0;

    PO->edge_indices=new_longvector(lines->nh);
    for (i=lines->nl;i<=lines->nh;i++){
    	area=area+area2d(centroid,lines->element[i]->begin,lines->element[i]->end);

    	PO->edge_indices->element[i]=lines->element[i]->index;


    }
    PO->area2D=area;



    free_longvector(lvertices);




	return PO;



}


LINEVECTOR *new_linevector(long nh){
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \return allocates a line vector with nh elements
	 */
	LINEVECTOR *lv;

	lv=(LINEVECTOR *)malloc(sizeof(LINEVECTOR));
	if (!lv) t_error("Linevector struct was not allocated");

	lv->isdynamic=isDynamic;
	lv->nl=NL;
	lv->nh=nh;


	lv->element=(LINE **)malloc((size_t)((nh-NL+1+NR_END)*sizeof(LINE *)));
	if (!lv->element) t_error("Linevector element struct was not allocated");

	return lv;

}


POLYGONVECTOR *new_polygonvector(long nh){
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \return allocates a polygon vector with nh elements
	 */
	POLYGONVECTOR *pv;

	pv=(POLYGONVECTOR *)malloc(sizeof(POLYGONVECTOR));
	if (!pv) t_error("Polygonvector struct was not allocated");

	pv->isdynamic=isDynamic;
	pv->nl=NL;
	pv->nh=nh;


	pv->element=(POLYGON **)malloc((size_t)((nh-NL+1+NR_END)*sizeof(POLYGON *)));
	if (!pv) t_error("Polygonvector element struct was not allocated");

	return pv;

}



