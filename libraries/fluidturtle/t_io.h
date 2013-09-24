#ifndef T_IO_H
#define T_IO_H
#include "turtle.h"
#include "tensors3D.h"
#include <vector>
#include "../../geotop/datastructs.h"

int mkdirp(const char *pathname, mode_t mode);

FILE *t_fopen(const char * ,const char *);FILE *t_fclose(FILE * stream);
char * join_strings(char *,char *);
char *join_strings(const char *, char *);

/**
Names: get_workingdirectory
Synopsis: char *get_workingdirectory(void );
Version: 1.0
Description: It asks for the working directory, i.e. the directory where inputs data are stored.
This directory can be alternatively specified in a file $WorkingPaths to be found in the executable directory
(for Windows or Mac systems where it is useful. In Unix systems it can be also the directory from where
the program containing this routine is called) FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c

Inputs: void
Return:  a pointer to the resulting string (the working directory )
*/

char *get_workingdirectory(void );

/*================functions copied from datamanipulation.c===========*/

/**
Name: initialize_
Version: 1.0
Synopsis:
Description:  It initialize the matrix or vector with  the specified value
Inputs:  1) The pointer to the structure to initalize; 2) the value used for initialization
*/


void initialize_longvector(LONGVECTOR *, long);
void initialize_shortvector(SHORTVECTOR *,short );
void initialize_intvector(INTVECTOR *,int );
void initialize_floatvector(FLOATVECTOR *,float );
void initialize_doublevector(DOUBLEVECTOR *,double );

void initialize_longmatrix(LONGMATRIX *, long);
void initialize_shortmatrix(SHORTMATRIX *,short );
void initialize_intmatrix(INTMATRIX *,int );
void initialize_floatmatrix(FLOATMATRIX *,float );
void initialize_doublematrix(DOUBLEMATRIX *,double );


/**
Name: copy__matrix
Description: copy a matrix into a  new one
Inputs: 1- the destination matrix; 2- the origin matrix
*/

void copy_longmatrix(LONGMATRIX *,LONGMATRIX *);
void copy_shortmatrix(SHORTMATRIX *,SHORTMATRIX *);
void copy_intmatrix(INTMATRIX *,INTMATRIX *);
void copy_floatmatrix(FLOATMATRIX *,FLOATMATRIX *);
void copy_doublematrix(DOUBLEMATRIX *,DOUBLEMATRIX *);
void copy_shortvector(SHORTVECTOR *,SHORTVECTOR *);
void copy_intvector(INTVECTOR *,INTVECTOR *);
void copy_longvector(LONGVECTOR *,LONGVECTOR *);
void copy_floatvector(FLOATVECTOR *,FLOATVECTOR *);
void copy_doublevector(DOUBLEVECTOR *,DOUBLEVECTOR *);
void add_doublevector(DOUBLEVECTOR *, DOUBLEVECTOR *);

void print_doublematrix_elements(DOUBLEMATRIX *,long );

void multipass_topofilter(long ntimes, DOUBLEMATRIX *Zin, DOUBLEMATRIX *Zout, long novalue, long n);
void multipass_topofilter(long ntimes, GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n);
void find_min_max(DOUBLEMATRIX *M, long novalue, double *max, double *min);
void find_min_max(GeoMatrix<double>& M, long novalue, double *max, double *min);

//double norm_1(DOUBLEVECTOR *V, long nbeg, long nend);
double norm_1(const GeoVector<double>& V, long nbeg, long nend);

/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of part in which you want divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti                                                    */
//void sky_view_factor(DOUBLEMATRIX *sky, long N, T_INIT *UV, DOUBLEMATRIX *input, SHORTMATRIX *convess, long novalue);
void sky_view_factor(DOUBLEMATRIX *sky, long N, TInit *UV, DOUBLEMATRIX *input, SHORTMATRIX *convess, long novalue);
//void sky_view_factor(GeoMatrix<double>& sky, long N, T_INIT *UV, GeoMatrix<double>& input, GeoMatrix<short>& convess, long novalue);
void sky_view_factor(GeoMatrix<double>& sky, long N, TInit *UV, GeoMatrix<double>& input, GeoMatrix<short>& convess, long novalue);

DOUBLEMATRIX *find_aspect(DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx, DOUBLEMATRIX *dzdy, long undef);
DOUBLEMATRIX *find_aspect(DOUBLEMATRIX *topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef);
void find_aspect(DOUBLEMATRIX *topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M);
void find_aspect(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M);

void topofilter(DOUBLEMATRIX *Zin, DOUBLEMATRIX *Zout, long novalue, long n);
void topofilter(GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n);

void find_slope(double deltax, double deltay, DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx, DOUBLEMATRIX *dzdy, long undef);
void find_slope(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef);

DOUBLEMATRIX *find_max_slope(DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx, DOUBLEMATRIX *dzdy, long undef);
void find_max_slope(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M);

//void product_matrix_using_lower_part_by_vector_plus_vector(double k, DOUBLEVECTOR *out, DOUBLEVECTOR *y, DOUBLEVECTOR *x, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);
  void product_matrix_using_lower_part_by_vector_plus_vector(double k, GeoVector<double>& out, const GeoVector<double>& y, const GeoVector<double>& x, const GeoVector<long>& Li, const GeoVector<long>& Lp, GeoVector<double>& Lx);

//void product_using_only_strict_lower_diagonal_part(DOUBLEVECTOR *product, DOUBLEVECTOR *x, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);
  void product_using_only_strict_lower_diagonal_part(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<long>& Li, const GeoVector<long>& Lp, GeoVector<double>& Lx);

//long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x,
//									DOUBLEVECTOR *b, DOUBLEVECTOR *y, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);
  long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(double tol_rel, double tol_min, double tol_max, GeoVector<double>& x,
		  	  	  	  	  	  	  	const GeoVector<double>& b, GeoVector<double>& y, const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);

//short tridiag(short a, long r, long c, long nx, DOUBLEVECTOR *diag_inf, DOUBLEVECTOR *diag, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e);
  short tridiag(short a, long r, long c, long nx, const GeoVector<double>& diag_inf, const GeoVector<double>& diag, const GeoVector<double>& diag_sup, const GeoVector<double>& b, GeoVector<double>& e);

//void get_diag_strict_lower_matrix_plus_identity_by_vector(DOUBLEVECTOR *diag, DOUBLEVECTOR *udiag, DOUBLEVECTOR *y,
//												   LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);
  void get_diag_strict_lower_matrix_plus_identity_by_vector(GeoVector<double>& diag, GeoVector<double>& udiag, const GeoVector<double>& y,
													const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);

//void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(DOUBLEVECTOR *product, DOUBLEVECTOR *x, DOUBLEVECTOR *y,
//																LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);

  void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<double>& y,
																const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);

//double product(DOUBLEVECTOR *a, DOUBLEVECTOR *b);
  double product(const GeoVector<double>& a, const GeoVector<double>& b);

//char *copy_stringnames(const char *origin);
#endif
