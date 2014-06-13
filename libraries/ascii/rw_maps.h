
#ifndef RW_MAPS_H
#define RW_MAPS_H

#include <string>
#include "../../geotop/datastructs.h"
#include "../fluidturtle/turtle.h"

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//BASE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void copyshort_doublematrix(GeoMatrix<short>& S, GeoMatrix<double>& M);

void copylong_doublematrix(GeoMatrix<long>& L, GeoMatrix<double>& M);

void copydoublematrix_const(double c0, GeoMatrix<double>& Mref,GeoMatrix<double>& M, double NOVALUE);


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//UTILITITY subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_suffix(std::string &suffix, long i, short start);

std::string namefile_i(std::string name, long i);
std::string namefile_i_we(std::string name, long i);
std::string namefile_i_we2(std::string name, long i);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//READ subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

GeoVector<double> read_map_vector(short type, std::string namefile, GeoMatrix<double>& mask, TInit *grid, double no_value, GeoMatrix<long>& rc);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//WRITE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_map(std::string filename, short type, short format, GeoMatrix<double>& M, TInit *UV, long novalue);
void write_map(std::string filename, short type, short format, GeoMatrix<long>& M, TInit *UV, long novalue);

void write_map_vector(std::string filename, short type, short format, const GeoVector<double>& V, TInit *UV, long novalue, long **j, long nr, long nc);

void write_tensorseries(short a, long l, long i, std::string filename, short type, short format, GeoTensor<double>& T, TInit *UV, long novalue);

void write_tensorseries_vector(short a, long l, long i, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries2_vector(std::string suf, long l, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);

void write_tensorseries3_vector(std::string suffix, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);


/*===================================================================================*/
/*===================functions copied from the file write_ascii.c====================*/
/*===================================================================================*/
void write_esriascii(std::string name, short type, GeoMatrix<double>& DTM, TInit *UV, long novalue);
void write_esriascii(std::string name, short type, GeoMatrix<long>& DTM, TInit *UV, long novalue);
void write_esriascii_vector(std::string name, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue);
/*==================================coppied functions from import_ascii.c============*/

#define max_figures 30

void error_message(short format, long n, long n1, long n2, long n3, char *name);

/*===================functions copied from geomorphology.0875.c======================*/

void nablaquadro_mask(GeoMatrix<double>& Z0,GeoMatrix<short>& curv,GeoVector<double>& U, GeoVector<double>& V);

void curvature(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& c1, GeoMatrix<double>& c2, GeoMatrix<double>& c3, GeoMatrix<double>& c4, long undef);

short is_boundary(long r, long c, GeoMatrix<double>& dem, long novalue);
long row(double N, long nrows, TInit *UV, long novalue);
long col(double E, long ncols, TInit *UV, long novalue);

/*===========copied function from init.h========*/
void initmatrix(double val, GeoMatrix<double>& destination, GeoMatrix<double>& origin, double novalue);

/* ===========copied data from extensions.h==========*/

char const * const ascii_grass = ".grass" ;
char const * const ascii_esri = ".asc" ;
char const * const textfile = ".txt" ;

#endif
