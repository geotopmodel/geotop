
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry_freememory.c

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
#include "t_datamanipulation.h"
#include "geometry_utilities.h"
#include "geometry_attribute.h"
#include "geometry_freememory.h"

void free_point(POINT *P){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 */


	free(P);
///	if (P) printf ("Warning: point %ld was not deallacated!!",P->index);

}

void free_line(LINE *L){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 */
	free_point(L->begin);
	free_point(L->end);
//	if (L->attributes) free_line_attributes(L->attributes);

	free(L);

}

void free_polygon(POLYGON *PO){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 */
	free_point(PO->centroid);
	free_longvector(PO->edge_indices);
	//if (PO->attributes) free_polygon_attributes(PO->attributes);

	free(PO);
}

void free_linevector(LINEVECTOR *lv){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 */
	long i;

	for (i=lv->nl;i<=lv->nh;i++){
		free_line(lv->element[i]);

	}

	free(lv);


}

void free_polygonvector(POLYGONVECTOR *pv){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 */
	long i;

	for (i=pv->nl;i<=pv->nh;i++){
		free_polygon(pv->element[i]);
	}

	free(pv);


}


/* The following routines lines must be rewritten and rebuilt if attribute struct type are changed!! */
void free_point_attributes(attribute_point *point){

	free(point);


}

void free_line_attributes(attribute_line *line){

	free(line);


}

void free_polygon_attributes(attribute_polygon *polygon){

	free(polygon);


}

void free_polygon_connection_attributes (polygon_connection_attributes *pca){

	free_longvector(pca->connections);
	free_doublevector(pca->d_connections);
	free(pca);
}

void  free_polygon_connection_attribute_array(polygon_connection_attribute_array *pcaa){

	long l;
	for (l=pcaa->nl;l<=pcaa->nh;l++){
		free_polygon_connection_attributes (pcaa->element[l]);
	}
	free(pcaa);
}
