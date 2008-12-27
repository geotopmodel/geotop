/*
 * geometry_utilities.h
 *
 *  Created on: Oct 22, 2008
 *      Author: ecor
 */

int signum (long a);

double area2d(POINT *P1,POINT *P2, POINT *P3);

int check_3_numbers(long a, long b, long c);

int check_vertex_intersection_2_lines(long b1,long e1,long b2,long e2);

long query_freq_longvector(LONGVECTOR *lvector,long value);

int segment_intersection_occurence(LINE *L1,LINE *L2);

int signum_double (double a);

double distance2D(POINT *P1,POINT *P2);

long shared_edges(POLYGON *PO1, POLYGON *PO2, long no_intersection, long *l_po1, long *l_po2, double *dist_centroids);
