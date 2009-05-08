
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry_attribute.c

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


//#include "geometry_io.h"
#include "geometry_attribute.h"

#define NO_INTERSECTION -99
#define NULL_VALUE -99

polygon_connection_attributes *get_connection(POLYGON *polygon,POLYGONVECTOR *polygons, long boundary, long displacement ,short print){
	/*
	 *
	 * \param polygon - (POLYGON *) polygon to which link are referred
	 * \param polygons - (POLYGONS *) vector of polygons
	 * \param boundary - (long) long value which indetifies the boundary
	 * \param displacement - (long) displacement in the POLYNGONVECTOR polygons around the index value of polygon where to find the connections
	 * \param print    - (short)
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 *\return a poygon_connection_atributes struct for a given polygon (polygon) within a polygon array (polygons)
	 *
	 */

	polygon_connection_attributes *pca;
	long l,s,l_po1,l_po2,icnt;
	double dist;
	long A_min,A_max; /* extremes of the search interval */

	icnt=1;

	pca=(polygon_connection_attributes *)malloc((sizeof(polygon_connection_attributes)));
	if (!pca) printf("Error: polygon_connection_attributes was not allocated at %ld polygon",polygon->index);
	pca->connections=new_longvector(polygon->edge_indices->nh);
	pca->d_connections=new_doublevector(polygon->edge_indices->nh);
	initialize_longvector(pca->connections,boundary);
	initialize_doublevector(pca->d_connections,NULL_VALUE);
	A_min=fmax(polygon->index-displacement,polygons->nl);
	A_max=fmin(polygon->index+displacement,polygons->nh);

	for (l=A_min;l<=A_max;l++){
		if (l!=polygon->index) {
			dist=1.0;
			l_po1=0;
			l_po2=0;

			s=shared_edges(polygon,polygons->element[l],NO_INTERSECTION,&l_po1,&l_po2,&dist);

			if (s!=NO_INTERSECTION ){
				if (l_po1>polygon->edge_indices->nh) printf ("Error: Line %ld : (%ld for polygon %ld) Not coherent data!!! \n ",s,l_po1,polygon->index);
				pca->connections->element[l_po1]=polygons->element[l]->index;
				pca->d_connections->element[l_po1]=dist;
//				printf(" pca conn. %ld %lf",pca->connections->element[l_po1],pca->d_connections->element[l_po1]);
			}
		}


	}


	return pca;

}






polygon_connection_attribute_array *new_connection_attributes(long nh){
	/*
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 * \return allocates a new polygon connection attributes vector with nh elements
	 */
	polygon_connection_attribute_array *apv;

	apv=(polygon_connection_attribute_array *)malloc(sizeof(polygon_connection_attribute_array));
	if (!apv) t_error("polygon_connection_attribute_array was not allocated");

	apv->isdynamic=isDynamic;
	apv->nl=NL;
	apv->nh=nh;


	apv->element=(polygon_connection_attributes **)malloc((size_t)((nh-NL+1+NR_END)*sizeof(polygon_connection_attributes *)));
	if (!apv) t_error("polygon_connection_attribute_array element struct was not allocated");

	return apv;

}

polygon_connection_attribute_array *get_connection_array(POLYGONVECTOR *polygons, long boundary, long displacement, short print){

	/*

	 * \param polygons - (POLYGONS *) vector of polygons
	 * \param boundary - (long) long value which indetifies the boundary
	 * \param displacement - (long) displacement in the POLYNGONVECTOR polygons around the index value of polygon where to find the connections
	 * \param print    - (short)
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 *\return build a a poygon_connection_atribute_array struct for a given array of polygons
	 *
	 *
	 *
	 */
	polygon_connection_attribute_array *pca;
	long l,j;

	pca=new_connection_attributes(polygons->nh);
	printf ("pca Allocated!\n");

	for(l=pca->nl;l<=pca->nh;l++){
		pca->element[l]=get_connection(polygons->element[l],polygons,boundary,displacement,print);
		if (print==1) printf ("Polygons %ld (%ld connections) of %ld (displacement=%ld)!! \n",l,pca->element[l]->connections->nh,pca->nh,displacement);
	}

	if (print==1) printf (" Finished polygon_connection_array struct filling!! \n ");



	return pca;
}


int write_polygonconnectionattributearray(char *filename,polygon_connection_attribute_array *pca) {
	/*
		 *
		 * \author Emanuele Cordano
		 * \data November 2008
		 *
		 * \param name - (char *)name of file where to write the linevector properties
		 * \param pca - (polygon_connection_attribute_array *) the polygon_connection_attribute_array to be printed
		 *
		 *\brief This functions writes a polygonconnectionattributearray in an ascii files with fluidturle formalism
		 *\return 0 if the polygon_connection_attribute_array is written correctly, -1 otherwise
		 *
		 */

	FILE *fd;
	long i,nh,nl,l,link_val;
	double d_link_val;
	nl=NL;
	nh=NL;
	for (i=pca->nl;i<=pca->nh;i++){
		if (nl>pca->element[i]->connections->nl) nl=pca->element[i]->connections->nl;
		if (nh<pca->element[i]->connections->nh) nh=pca->element[i]->connections->nh;
	}

	fd=fopen(filename,"w");
	/* PRINT HEADER */
	fprintf(fd,"index{%ld}\n",pca->nh);
	fprintf(fd,"/**  FILE CONTAINIG NECESSARY INFORMATION FOR POLOYGON CONNECTION  \n");
	fprintf(fd," each polygon connection array is expressed as a double array containing the following information \n");
	fprintf(fd,"polygon_index    ");

	for (l=nl;l<=nh;l++) {
		fprintf(fd,"P%ld    ",l);
		fprintf(fd,"d_P%ld    ",l);
	}
	fprintf(fd," */  \n  ");

	/* end print header */
	for (i=pca->nl;i<=pca->nh;i++){
		fprintf(fd,"\n");
		fprintf(fd,"%ld: double array connections %ld {",i,i);
		fprintf(fd,"%ld    ",i);

		for (l=nl;l<=nh;l++){
			link_val=NULL_VALUE;
			d_link_val=NULL_VALUE;
			if ((l<=pca->element[i]->connections->nh) && (l>=pca->element[i]->connections->nl)) link_val=pca->element[i]->connections->element[l];
			if ((l<=pca->element[i]->d_connections->nh) && (l>=pca->element[i]->d_connections->nl)) d_link_val=pca->element[i]->d_connections->element[l];
			fprintf(fd,",%ld,   ",link_val);
			fprintf(fd,"%lf    ",d_link_val);
		}
		fprintf(fd,"}");
	}


	fprintf(fd,"\n");
	fclose(fd);
	return 0;

}


int fprint_polygonconnectionattributearray(char *filename,polygon_connection_attribute_array *pca) {
	/*
		 *
		 * \author Emanuele Cordano
		 * \data November 2008
		 *
		 * \param name - (char *)name of file where to write the linevector properties
		 * \param pca - (polygon_connection_attribute_array *) the polygon_connection_attribute_array to be printed
		 *
		 *\return 0 if the polygon_connection_attribute_array is written correctly, -1 otherwise
		 *
		 */

	FILE *fd;
	long i,nh,nl,l,link_val;
	double d_link_val;
	nl=NL;
	nh=NL;
	for (i=pca->nl;i<=pca->nh;i++){
		if (nl>pca->element[i]->connections->nl) nl=pca->element[i]->connections->nl;
		if (nh<pca->element[i]->connections->nh) nh=pca->element[i]->connections->nh;
	}

	fd=fopen(filename,"w");
	/* PRINT HEADER */


	fprintf(fd,"polygon_index    ");

	for (l=nl;l<=nh;l++) {
		fprintf(fd,"P%ld    ",l);
		fprintf(fd,"d_P%ld    ",l);
	}
	/* end print header */
	for (i=pca->nl;i<=pca->nh;i++){
		fprintf(fd,"\n");

		fprintf(fd,"%ld    ",i);

		for (l=nl;l<=nh;l++){
			link_val=NULL_VALUE;
			d_link_val=NULL_VALUE;
			if ((l<=pca->element[i]->connections->nh) && (l>=pca->element[i]->connections->nl)) link_val=pca->element[i]->connections->element[l];
			if ((l<=pca->element[i]->d_connections->nh) && (l>=pca->element[i]->d_connections->nl)) d_link_val=pca->element[i]->d_connections->element[l];
			fprintf(fd,"%ld    ",link_val);
			fprintf(fd,"%lf    ",d_link_val);
		}

	}


	fprintf(fd,"\n");
	fclose(fd);
	return 0;

}


int connections_symmetry(polygon_connection_attribute_array* pca, long boundary) {
	/*
	 *\param  pca - (polygon_connection_attrivute_array *) the attributte array of which the symmetry is verified;
	 *\param boundary - (long) - boundary indicator
	 *
	 *\author Emanuele Cordano
	 *\date  March 2008
	 *
	 *
	 */
	long i,j,kl,c;
	int s=1,sk; /* symmetry is assumed */

	for (i=pca->nl;i<=pca->nh;i++) {
		for (j=pca->element[i]->connections->nl;j<=pca->element[i]->connections->nh;j++) {
			if (pca->element[i]->connections->element[j]!=boundary){

				kl=pca->element[i]->connections->element[j];
				sk=0;
				for (c=pca->element[kl]->connections->nl;c<=pca->element[kl]->connections->nh;c++) {
					if (pca->element[kl]->connections->element[c]==i) sk=1;
				}
				if (sk==0) {
					s=0;
					printf("Error in pca: struct is not symmetric in the calls %ld , %ld      /n",i,kl);
				}
			}
		}
	}

	return s;

}
