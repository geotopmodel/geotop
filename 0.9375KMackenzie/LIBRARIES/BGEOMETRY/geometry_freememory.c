/*
 * attrites_boussinesq1.c
 *
 *  Created on: Nov 3, 2008
 *      Author: ecor
 */

#include "turtle.h"
#include "geometry.h"
#include "geometry_freememory.h"
#include "geometry_attribute.h"
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
