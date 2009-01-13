
/* BGEOMETRY BUILDS THE MESH FROM A RASTER FOR THE BOUSSINESQ MODEL
KeyPalette Version 0.9375 KMackenzie

file geometry_utilities.c

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

int signum (long a){
/*
 * \author Emanuele Cordano
 *
 * \date October 2008
 *
 */
if (a==0){
	return 0;
}else if(a>0) {
	return 1;
} else if(a<0){
	return -1;
}
return 1;
}

int signum_double (double a){
/*
 * \author Emanuele Cordano
 *
 * \date October 2008
 *
 */
if (a==0.0){
	return 0;
}else if(a>0) {
	return 1;
} else if(a<0){
	return -1;
}
return 1;
}

int segment_intersection_occurence(LINE *L1,LINE *L2){
/*
 * \author Emanuele Cordano
 * \date October 2008
 * \return 0 if the lines are not connected, 1 if they share only 1 point, 2 if they share more than 1 point

 *
 */
//int l;
double det,det1,det2;
double lambda1,lambda2;

	/* Apply Cramer's Rule (http://en.wikipedia.org/wiki/Cramer%27s_rule)*/
det=-(L1->end->x-L1->begin->x)*(L2->end->y-L2->begin->y)+(L2->end->x-L2->begin->x)*(L1->end->y-L1->begin->y);
det1=-(L2->begin->x-L1->begin->x)*(L2->end->y-L2->begin->y)+(L2->end->x-L2->begin->x)*(L2->begin->y-L1->begin->y);
det2=(L1->end->x-L1->begin->x)*(L2->begin->y-L1->begin->y)-(L1->end->y-L1->begin->y)*(L2->begin->x-L1->begin->x);

if (det==0.0) {
	if ((det1==0) || (det2==0)) {
		return 2;
	}else {
		return 0;
	}
} else {
	lambda1=det1/det;
	lambda2=det2/det;
	if ((lambda1>0.0) && (lambda1<1.0) && (lambda2>0.0) && (lambda2<1.0)) {
		return 1;
	} else{
	return 0;
	}

}
return 0;

}

long query_freq_longvector(LONGVECTOR *lvector,long value){
/*
 *
 * \author Emanuele Cordano
 * \\date October 2008
 *
 */
long j,cnt;
cnt=0;
for (j=lvector->nl;j<=lvector->nh;j++){
	if (lvector->element[j]==value) cnt++;
}
return cnt;
}

int check_vertex_intersection_2_lines(long b1,long e1,long b2,long e2){
	/*
	 * \param (long) b1 - begin point index of line 1
	 * \param (long) e1 - end   point index of line 1
	 * \param (long) b2 - begin point index of line 2
	 * \param (long) e2 - end   point index of line 2
	 *
	 * \return 0 if the lines are not connected, 1 if they share only 1 edge, 2 if they share both the edges
	 * \author Emanuele Cordano
	 * \data October 2008
	 *
	 */


	if (signum(b1-e2)*signum(b1-b2)*signum(e1-e2)*signum(e1-b2)==0) {
		if ((signum(b1-e2)*signum(b1-e2)+signum(e1-b2)*signum(e1-b2))*(signum(b1-b2)*signum(b1-b2)+signum(e1-e2)*signum(e1-e2))==0)   {
			return 2;/* the two lines coincide */
		}else {
			return 1; /* the lines share a vertex */
		}

	}
	/* the lines has no intersection at the edges  */
	return 0;
}

int check_3_numbers(long a, long b, long c) {
   /*
    * \param (long) a,b,c - three long integer values
    * \author Emanuele Cordano
    * \date October 2008
    *
    * \return 1 if a=b=c , 0 otherwise
    */

//if (((a-b)*(a-b)+(a-c)*(a-c))==0) return 1;
	if ((a==b) && (a=c)) return TRUE;
	return FALSE;
}


double area2d(POINT *P1,POINT *P2, POINT *P3) {
	/*
	 * \return calculetes the area of the triangle among three points;
	 * \author Emanuele Cordano
	 * \date October 2008
	 *
	 * \brief area's formulas on http://it.wikipedia.org/wiki/Triangolo
	 */
	double area;


	area=(P1->x*(P3->y-P2->y)+P2->x*(P1->y-P3->y)+P3->x*(P2->y-P1->y))/2.0;
//	printf("x= %lf %lf %lf \n",P1->x,P2->x,P3->x);
//	printf("y= %lf %lf %lf \n",P1->y,P2->y,P3->y);
//	printf ("area= %lf",area);

	if (area<0.0) area=-area;

	return area;

}


double distance2D(POINT *P1,POINT *P2){
/*
 * \return euclidean 2d distance between points P1 and P2
 * \author Emanuele Cordano
 * \date October 2008
 */
	return pow(pow(P1->x-P2->x,2.0)+pow(P1->y-P2->y,2.0),0.5);
}

long shared_edges(POLYGON *PO1, POLYGON *PO2, long no_intersection, long *l_po1, long *l_po2, double *dist_centroids){
/*
 *
 *
 * \author Emanuele Cordano
 * \date November 2008
 *
 * \param PO1 - (POLYGON *) the first polygon
 * \param PO2 - (POLYGON *) the second polygon
 * \param no_intersection - (long) integer vale which is returend in case PO1 and PO2 do not share any edge.
 * \param l_po1 - (double) number of the shared edge respect to the polygon PO1
 * \param l_po2 - (double) number of the shared edge respect to the polygon PO2
 *
 * \param dist_centroids -  (double) distance between the two centroids of the polygons
 *
 * \return the index of the edge shared between the two polygons;
 */

long l,k,l_shared;
double distance;
int s;


l_shared=no_intersection;
s=0;

for(l=PO1->edge_indices->nl;l<=PO1->edge_indices->nh;l++){
	for(k=PO2->edge_indices->nl;k<=PO2->edge_indices->nh;k++){
		if (PO1->edge_indices->element[l]==PO2->edge_indices->element[k]) {
			if (s==0){
				l_shared=PO1->edge_indices->element[l];
				*l_po1=l;
				*l_po2=k;
				s=1;
			}else {
				printf ("Warning:: polygons %ld and %ld share more than one edge (%ld,%ld)!! \n",PO1->index,PO2->index,PO1->edge_indices->element[l],l_shared);
				l_shared=no_intersection;
			}

		}
	}
}

distance=distance2D(PO1->centroid,PO2->centroid);
*dist_centroids=distance;
//printf("... PO1=%ld PO2=%ld dist=%lf",PO1->index,PO2->index,distance);
//stop_execution();
return l_shared;

}


