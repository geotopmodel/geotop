
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry_io.c

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
//#include "additional_read_functions.h"
#include "geometry_utilities.h"
#include "geometry_io.h"
#include "geometry_freememory.h"


#define DELIMITERS ";"
#define MAX_POINTS 20
#define NULL_VALUE -99
#define NO_ELEVATION 0.0

#define NO_COLS_LINEVECTOR_DOUBLEMATRIX 8

LINE *get_line(DOUBLEVECTOR *vertex_x_coord, DOUBLEVECTOR *vertex_y_coord, long line_index, char *number_strings, short print){
	/*
	 * \param vertex_x_coord - (DOUBLEVECTOR *) vector of x coordinates
	 * \param vertex_y_coord - (DOUBLEVECTOR *) vector of y coordinates
	 *
	 * \param line_index     - (long int) index of the line
	 * \param number_string  - (char *) string containing the number of two vertex points
	 * \param print          - (short) if activated it prints error or warning messages
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 */
	POINT *P1, *P2;
	long index1,index2,j;
	LONGVECTOR *lvertices;
	LINE *li;

	lvertices=read_longarray_from_string(number_strings,DELIMITERS,MAX_POINTS,0);

    if (lvertices->nh<2) printf("Error:: line %s (%ld) cannot be created, there no vertices!!\n",number_strings,line_index);
    for (j=lvertices->nl;j<=2;j++){
    	if ((lvertices->element[j]>vertex_x_coord->nh) || (lvertices->element[j]>vertex_y_coord->nh)) {
    		printf ("Error:: line %s (%ld) cannot be created, point %ld has no coordinates or no atributes!!\n",number_strings,line_index,j);
    	}
    }

    j=1;
    P1=new_point(lvertices->element[j],vertex_x_coord->element[lvertices->element[j]],vertex_y_coord->element[lvertices->element[j]],NO_ELEVATION);

    j=2;
    P2=new_point(lvertices->element[j],vertex_x_coord->element[lvertices->element[j]],vertex_y_coord->element[lvertices->element[j]],NO_ELEVATION);

    if (print==1) printf ("The two vertices %ld (x=%lf,y=%lf) and %ld (x=%lf,y=%lf) of line %s (%ld) were succesfully created!!\n",P1->index,P1->x,P1->y,P2->index,P2->x,P2->y,number_strings,line_index);

    li=new_line_from_points(line_index, P1,P2);


    free_point(P1);
    free_point(P2);
    free_longvector(lvertices);


    return li;



}



LINEVECTOR *get_linevector(DOUBLEVECTOR *vertex_x_coord, DOUBLEVECTOR *vertex_y_coord,STRINGBIN *line_strings , short print)
{
	/*
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \param vertex_x_coord - (DOUBLEVECTOR *) vector of x coordinates
	 * \param vertex_y_coord - (DOUBLEVECTOR *) vector of y coordinates
	 *
	 * \param line_strings   - (STRINGBIN *) string array containing line information
	 * \param print          - (short int)
	 *
	 *
	 */
	LINEVECTOR *lv;
	long l;


	//if (attribute_line_array->nh!=line_strings->index->nh) printf("Error: the line string array and attribute array has different length!!!\n");

	lv=new_linevector(line_strings->index->nh);

	for(l=lv->nl;l<=lv->nh;l++){
		lv->element[l]=get_line(vertex_x_coord,vertex_y_coord,l,line_strings->element[l]+1,print);

	}




	return lv;

}

LINEVECTOR *extract_linvector_from_linevector(LONGVECTOR *nlines, LINEVECTOR *lines){
/*
 * \author Emanuele Cordano
 * \date November 2008
 *
 * \param nlines - (LONGVECTOR *) numbers of lines to be extracted
 * \param lines  - (LINEVECTOR *) vector of all lines
 *
 * \return a reduced linevoctor with the lines extracted
 */
	LINEVECTOR *lve;

	long j;
//	print_longvector_elements(nlines,1);


    lve=new_linevector(nlines->nh);

    for (j=lve->nl;j<=lve->nh;j++){
    	lve->element[j]=new_line_from_points(nlines->element[j],lines->element[nlines->element[j]]->begin,lines->element[nlines->element[j]]->end);

    }


    return lve;

}

POLYGON *get_polygon(long index,double x, double y, double z, char *edge_index_string,LINEVECTOR *all_lines, short print) {
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \param index - (long) index of the centroid and the polygons
	 * \param x  - (double) x of the centroid
	 * \param y  - (double) y of the centroid
	 * \param z  - (double) z of the centroid
	 *
	 * \param edge_line_index - (char *) string containing the edge indices
	 *
	 * \param all_lines       -  (LINEVECTOR *) vector of all lines
	 * \param print   - (short)
	 *
	 * \return a polygon (Attributes are not allocated!!)
	 *
	 */

	LINEVECTOR *edges;
	POLYGON *PO;
	LONGVECTOR *ledges;
	POINT *centroid;

	centroid=new_point(index,x,y,z);
	ledges=read_longarray_from_string(edge_index_string,DELIMITERS,MAX_POINTS,0);
	edges=extract_linvector_from_linevector(ledges,all_lines);
	PO=new_polygon_from_a_linevector(edges,centroid);

	free_point(centroid);
	free_longvector(ledges);
	free_linevector(edges);

	return PO;


}

POLYGONVECTOR *get_polygonvector(DOUBLEVECTOR *centroid_x_coord, DOUBLEVECTOR *centroid_y_coord,STRINGBIN *polygon_strings ,LINEVECTOR *all_lines,short print){
/*
 * \author Emanuele Cordano
 * \date November 2008
 *
 *\param centroid_x_coord - (DOUBLEVECTOR*) vector of centroid x coordinates
 *\param centroid_y_coord - (DOUBLEVECTOR*) vector of centroid y coordinates
 *\param polygon_strings -  (STRINGBIN *) string containing polygon information
 *\param lines            - (LINEVECTOR *) vector containing all the lines
 *\param print            - (short)
 *
 *
 *
 */
POLYGONVECTOR *pv;
long l;

if ((centroid_x_coord->nh!=centroid_y_coord->nh) && (polygon_strings->index->nh!=centroid_y_coord->nh) && (polygon_strings->index->nh!=centroid_x_coord->nh)) printf ("Error:: Polygons stringarray (STRINGBIN *) and centroid coordinete vectors have a different number of elements!!\n");
pv=new_polygonvector(polygon_strings->index->nh);

for(l=pv->nl;l<=pv->nh;l++){
	pv->element[l]=get_polygon(l,centroid_x_coord->element[l],centroid_y_coord->element[l],NO_ELEVATION,polygon_strings->element[l]+1,all_lines,print);
}

return pv;

}

int write_linevector(char *filename, LINEVECTOR *lines){
	/*
		 *
		 * \author Emanuele Cordano
		 * \data November 2008
		 *
		 * \param name - (char *)name of file where to write the linevector properties
		 * \param linevector - (LINEVECTOR *) the linevector to be printed
		 *
		 *\brief write a linevecto in a fluidturtle formalism.
		 *\brief
		 *\return 0 if the linevector is written correctly, -1 otherwise
		 *
		 */
	FILE *fd;
	long i;
	double xm,ym;

	fd=fopen(filename,"w");
	fprintf(fd,"index{1}\n");
	fprintf(fd,"/**  FILE CONTAINIG NECESSARY INFORMATION FOR LINES \n");
	fprintf(fd,"x    ");
	fprintf(fd,"y    ");

	fprintf(fd,"line_index    ");
	fprintf(fd,"lenght2d    ");
	fprintf(fd,"x_P1    ");
	fprintf(fd,"y_P1    ");
	fprintf(fd,"x_P2   ");
	fprintf(fd,"y_P2   \n ");
	fprintf(fd,"*/ \n");
	fprintf(fd,"1: double matrix lines information  {%ld,%ld}",lines->nh,NO_COLS_LINEVECTOR_DOUBLEMATRIX);
		for (i=lines->nl;i<=lines->nh;i++){
	//		printf("entrato!!! %ld\n",i);
			fprintf(fd,"\n");
			xm=(lines->element[i]->begin->x+lines->element[i]->end->x)/2.0;
			ym=(lines->element[i]->begin->y+lines->element[i]->end->y)/2.0;
			fprintf(fd,"%lf    ",xm);
			fprintf(fd,"%lf    ",ym);
			fprintf(fd,"%ld    ",lines->element[i]->index);
			fprintf(fd,"%lf    ",lines->element[i]->length2d);
				fprintf(fd,"%lf    ",lines->element[i]->begin->x);
				fprintf(fd,"%lf    ",lines->element[i]->begin->y);
				fprintf(fd,"%lf   ",lines->element[i]->end->x);
				fprintf(fd,"%lf    ",lines->element[i]->end->y);
		}


		fclose(fd);
		return 0;
}

int write_polygonvector(char *filename, POLYGONVECTOR *polygons) {
	/*
		 *
		 * \author Emanuele Cordano
		 * \data May 2009
		 *
		 * \param name - (char *)name of file where to write the linevector properties
		 * \param polygons - (POLYGONVECTOR *) the polygonvector to be printed
		 *
		 *\brief This functions writes a polygonvector in an ascii files with fluidturle formalism
		 *
		 *\return 0 if the polygonvector is written correctly, -1 otherwise
		 *
		 */

	FILE *fd;
	long i,nh,nl,l,lval;
	nl=NL;
	nh=NL;
	for (i=polygons->nl;i<=polygons->nh;i++){
		if (nl>polygons->element[i]->edge_indices->nl) nl=polygons->element[i]->edge_indices->nl;
		if (nh<polygons->element[i]->edge_indices->nh) nh=polygons->element[i]->edge_indices->nh;
	}

	fd=fopen(filename,"w");
	/* PRINT HEADER */
	fprintf(fd,"index{%ld}\n",polygons->nh);
	fprintf(fd,"/**  FILE CONTAINIG NECESSARY INFORMATION FOR POLOYGONS  \n");
	fprintf(fd," each polygon  is expressed as a double array coteining the following information \n");
	fprintf(fd,"x     ");
	fprintf(fd,"y    ");

	fprintf(fd,"polygon_index    ");
	fprintf(fd,"area2d    ");
	for (l=nl;l<=nh;l++) {
		fprintf(fd,"L%ld    ",l);
	}
	fprintf(fd,"  */ \n");

	/* end print header */
	for (i=polygons->nl;i<=polygons->nh;i++){
		fprintf(fd,"\n/** %ld block x     ",i);
		fprintf(fd,"y    ");

			fprintf(fd,"polygon_index    ");
			fprintf(fd,"area2d    ");
			for (l=nl;l<=nh;l++) {
				fprintf(fd,"L%ld    ",l);
			}
		fprintf(fd," */ \n");
		fprintf(fd,"%ld: double array {",i);
		fprintf(fd,"%lf,    ",polygons->element[i]->centroid->x);
		fprintf(fd,"%lf,    ",polygons->element[i]->centroid->y);
		fprintf(fd,"%ld,    ",polygons->element[i]->index);
		fprintf(fd,"%lf    ",polygons->element[i]->area2D);
		for (l=nl;l<=nh;l++){
			lval=NULL_VALUE;
			if ((l<=polygons->element[i]->edge_indices->nh) && (l>=polygons->element[i]->edge_indices->nl)) lval=polygons->element[i]->edge_indices->element[l];
			fprintf(fd,",%ld    ",lval);
		}
		fprintf(fd,"}");
	}


	fprintf(fd,"\n");
	fclose(fd);
	return 0;

}

int fprint_linevector(char *filename, LINEVECTOR *lines){
	/*
	 *
	 * \author Emanuele Cordano
	 * \data November 2008
	 *
	 * \param name - (char *)name of file where to write the linevector properties
	 * \param linevector - (LINEVECTOR *) the linevector to be printed
	 *
	 *
	 *\return 0 if the linevector is written correctly, -1 otherwise
	 *
	 */

	FILE *fd;
	long i;
	double xm,ym;

	fd=fopen(filename,"w");
	fprintf(fd,"x    ");
	fprintf(fd,"y    ");

	fprintf(fd,"line_index    ");
	fprintf(fd,"lenght2d    ");
	fprintf(fd,"x_P1    ");
	fprintf(fd,"y_P1    ");
	fprintf(fd,"x_P2   ");
	fprintf(fd,"y_P2   ");


	for (i=lines->nl;i<=lines->nh;i++){
//		printf("entrato!!! %ld\n",i);
		fprintf(fd,"\n");
		xm=(lines->element[i]->begin->x+lines->element[i]->end->x)/2.0;
		ym=(lines->element[i]->begin->y+lines->element[i]->end->y)/2.0;
		fprintf(fd,"%lf    ",xm);
		fprintf(fd,"%lf    ",ym);
		fprintf(fd,"%ld    ",lines->element[i]->index);
		fprintf(fd,"%lf    ",lines->element[i]->length2d);
			fprintf(fd,"%lf    ",lines->element[i]->begin->x);
			fprintf(fd,"%lf    ",lines->element[i]->begin->y);
			fprintf(fd,"%lf   ",lines->element[i]->end->x);
			fprintf(fd,"%lf    ",lines->element[i]->end->y);
	}


	fclose(fd);
	return 0;

}

int fprint_polygonvector(char *filename, POLYGONVECTOR *polygons) {
	/*
		 *
		 * \author Emanuele Cordano
		 * \data November 2008
		 *
		 * \param name - (char *)name of file where to write the linevector properties
		 * \param polygons - (POLYGONVECTOR *) the polygonvector to be printed
		 *
		 *\return 0 if the polygonvector is written correctly, -1 otherwise
		 *
		 */

	FILE *fd;
	long i,nh,nl,l,lval;
	nl=NL;
	nh=NL;
	for (i=polygons->nl;i<=polygons->nh;i++){
		if (nl>polygons->element[i]->edge_indices->nl) nl=polygons->element[i]->edge_indices->nl;
		if (nh<polygons->element[i]->edge_indices->nh) nh=polygons->element[i]->edge_indices->nh;
	}

	fd=fopen(filename,"w");
	/* PRINT HEADER */


	fprintf(fd,"x     ");
	fprintf(fd,"y    ");

	fprintf(fd,"polygon_index    ");
	fprintf(fd,"area2d    ");
	for (l=nl;l<=nh;l++) {
		fprintf(fd,"L%ld    ",l);
	}
	/* end print header */
	for (i=polygons->nl;i<=polygons->nh;i++){
		fprintf(fd,"\n");
		fprintf(fd,"%lf    ",polygons->element[i]->centroid->x);
		fprintf(fd,"%lf    ",polygons->element[i]->centroid->y);
		fprintf(fd,"%ld    ",polygons->element[i]->index);
		fprintf(fd,"%lf    ",polygons->element[i]->area2D);
		for (l=nl;l<=nh;l++){
			lval=NULL_VALUE;
			if ((l<=polygons->element[i]->edge_indices->nh) && (l>=polygons->element[i]->edge_indices->nl)) lval=polygons->element[i]->edge_indices->element[l];
			fprintf(fd,"%ld    ",lval);
		}

	}


	fprintf(fd,"\n");
	fclose(fd);
	return 0;

}
