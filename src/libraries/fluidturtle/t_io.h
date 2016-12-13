#ifndef T_IO_H
#define T_IO_H
#include "turtle.h"
#include "tensors3D.h"
#include <vector>
//#include "../../geotop/datastructs.h"

int mkdirp(const char *pathname, mode_t mode);

FILE *t_fopen(const char * ,const char *);FILE *t_fclose(FILE * stream);
//char * join_strings(char *,char *);
//char *join_strings(char const * const first, char const * const second) ;
//char *join_strings(const char *, char *);

/**
  * @brief It asks for the working directory, i.e. the directory where inputs data are stored.
  * This directory can be alternatively specified in a file $WorkingPaths to be found in the executable directory
  * (for Windows or Mac systems where it is useful. In Unix systems it can be also the directory from where
  * the program containing this routine is called) FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c
  *
  * @return the resulting string (the working directory )
  */
std::string get_workingdirectory();

/*================functions copied from datamanipulation.c===========*/

/**
Name: initialize_
Version: 1.0
Synopsis:
Description:  It initialize the matrix or vector with  the specified value
Inputs:  1) The pointer to the structure to initalize; 2) the value used for initialization
*/



/**
Name: copy__matrix
Description: copy a matrix into a  new one
Inputs: 1- the destination matrix; 2- the origin matrix
*/


void multipass_topofilter(long ntimes, GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n);
void find_min_max(GeoMatrix<double>& M, long novalue, double *max, double *min);

double norm_1(const GeoVector<double>& V, long nbeg, long nend);

/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of part in which you want divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti                                                    */
void sky_view_factor(GeoMatrix<double>& sky, long N, TInit *UV, GeoMatrix<double>& input, GeoMatrix<short>& convess, long novalue);
void find_aspect(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M);

void topofilter(GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n);

void find_slope(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef);

void find_max_slope(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M);

void product_matrix_using_lower_part_by_vector_plus_vector(double k, GeoVector<double>& out, const GeoVector<double>& y, const GeoVector<double>& x, const GeoVector<long>& Li, const GeoVector<long>& Lp, GeoVector<double>& Lx);

void product_using_only_strict_lower_diagonal_part(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<long>& Li, const GeoVector<long>& Lp, GeoVector<double>& Lx);

long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(double tol_rel, double tol_min, double tol_max, GeoVector<double>& x,
		  	  	  	  	  	  	  	const GeoVector<double>& b, GeoVector<double>& y, const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);

 short tridiag(short a, long r, long c, long nx, const GeoVector<double>& diag_inf, const GeoVector<double>& diag, const GeoVector<double>& diag_sup, const GeoVector<double>& b, GeoVector<double>& e);

void get_diag_strict_lower_matrix_plus_identity_by_vector(GeoVector<double>& diag, GeoVector<double>& udiag, const GeoVector<double>& y,const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);
void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<double>& y,const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx);

double product(const GeoVector<double>& a, const GeoVector<double>& b);

#endif
