
#ifndef RW_MAPS_H
#define RW_MAPS_H

#include "../fluidturtle/turtle.h"
#include "../fluidturtle/t_io.h"
#include "../fluidturtle/tensors3D.h"

#include <string>
#include "../../meteoio_plugin/meteoioplugin.h"

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//BASE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

//SHORTMATRIX *copyshort_doublematrix(DOUBLEMATRIX *M);
GeoMatrix<short> copyshort_doublematrix(DOUBLEMATRIX *M);
void copyshort_doublematrix(GeoMatrix<short>& S, GeoMatrix<double>& M);

//LONGMATRIX *copylong_doublematrix(DOUBLEMATRIX *M);
GeoMatrix<long> copylong_doublematrix(DOUBLEMATRIX *M);
void copylong_doublematrix(GeoMatrix<long>& L, GeoMatrix<double>& M);

DOUBLEMATRIX *copydouble_shortmatrix(SHORTMATRIX *S);

DOUBLEMATRIX *copydouble_longmatrix(LONGMATRIX *L);

//DOUBLEMATRIX *copydoublematrix_const(double c0, DOUBLEMATRIX *Mref, double NOVALUE);
//DOUBLEMATRIX *copydoublematrix_const(double c0, GeoMatrix<double>& Mref, double NOVALUE);
void copydoublematrix_const(double c0, GeoMatrix<double>& Mref,GeoMatrix<double>& M, double NOVALUE);


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//UTILITITY subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_suffix(char *suffix, long i, short start);
void write_suffix(std::string suffix, long i, short start);

//short existing_file(char *name);
//short existing_file_wext(char *name, char *extension);

char *namefile_i(char *name, long i);
char *namefile_i_we(char *name, long i);
char *namefile_i_we2(char *name, long i);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//READ subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

//DOUBLEMATRIX *read_map(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref, double no_value);
//DOUBLEVECTOR *read_map_vector(short type, char *namefile, DOUBLEMATRIX *mask, T_INIT *grid, double no_value, LONGMATRIX *rc);
GeoVector<double> read_map_vector(short type, char *namefile, GeoMatrix<double>& mask, TInit *grid, double no_value, GeoMatrix<long>& rc);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//WRITE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

//void write_map(char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV, long novalue);
void write_map(char *filename, short type, short format, DOUBLEMATRIX *M, TInit *UV, long novalue);
//void write_map(char *filename, short type, short format, GeoMatrix<double>& M, T_INIT *UV, long novalue);
void write_map(char *filename, short type, short format, GeoMatrix<double>& M, TInit *UV, long novalue);
void write_map(std::string filename, short type, short format, GeoMatrix<double>& M, TInit *UV, long novalue);
//void write_map(char *filename, short type, short format, GeoMatrix<long>& M, T_INIT *UV, long novalue);
void write_map(char *filename, short type, short format, GeoMatrix<long>& M, TInit *UV, long novalue);
void write_map(std::string filename, short type, short format, GeoMatrix<long>& M, TInit *UV, long novalue);

//void write_map_vector(char *filename, short type, short format, DOUBLEVECTOR *V, T_INIT *UV, long novalue, long **j, long nr, long nc);
void write_map_vector(char *filename, short type, short format, const GeoVector<double>& V, TInit *UV, long novalue, long **j, long nr, long nc);
void write_map_vector(std::string filename, short type, short format, const GeoVector<double>& V, TInit *UV, long novalue, long **j, long nr, long nc);

//void write_tensorseries(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);
void write_tensorseries(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, TInit *UV, long novalue);
//void write_tensorseries(short a, long l, long i, char *filename, short type, short format, GeoTensor<double>& T, T_INIT *UV, long novalue);
void write_tensorseries(short a, long l, long i, char *filename, short type, short format, GeoTensor<double>& T, TInit *UV, long novalue);

//void write_tensorseries_vector(short a, long l, long i, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc);
//void write_tensorseries_vector(short a, long l, long i, char *filename, short type, short format, DOUBLEMATRIX *T, TInit *UV, long novalue, long **J, long nr, long nc);
//void write_tensorseries_vector(short a, long l, long i, char *filename, short type, short format, GeoMatrix<double>& T, T_INIT *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries_vector(short a, long l, long i, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);
//void write_tensorseries2(char *suf, long l, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);
void write_tensorseries2(char *suf, long l, char *filename, short type, short format, DOUBLETENSOR *T, TInit *UV, long novalue);
//void write_tensorseries2_vector(char *suf, long l, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries2_vector(char *suf, long l, char *filename, short type, short format, DOUBLEMATRIX *T, TInit *UV, long novalue, long **J, long nr, long nc);
//void write_tensorseries2_vector(char *suf, long l, char *filename, short type, short format, GeoMatrix<double>& T, T_INIT *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries2_vector(char *suf, long l, char *filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries2_vector(std::string suf, long l, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);

//void write_tensorseries3(char *suffix, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);
void write_tensorseries3(char *suffix, char *filename, short type, short format, DOUBLETENSOR *T, TInit *UV, long novalue);
//void write_tensorseries3_vector(char *suffix, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries3_vector(char *suffix, char *filename, short type, short format, DOUBLEMATRIX *T, TInit *UV, long novalue, long **J, long nr, long nc);
//void write_tensorseries3_vector(char *suffix, char *filename, short type, short format, GeoMatrix<double>& T, T_INIT *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries3_vector(char *suffix, char *filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);
void write_tensorseries3_vector(std::string suffix, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc);


/*===================================================================================*/
/*===================functions copied from the file write_ascii.c====================*/
/*===================================================================================*/
//void write_esriascii(char *name, short type, DOUBLEMATRIX *DTM, T_INIT *UV, long novalue);
void write_esriascii(char *name, short type, DOUBLEMATRIX *DTM, TInit *UV, long novalue);
//void write_esriascii(char *name, short type, GeoMatrix<double>& DTM, T_INIT *UV, long novalue);
void write_esriascii(char *name, short type, GeoMatrix<double>& DTM, TInit *UV, long novalue);
void write_esriascii(std::string name, short type, GeoMatrix<double>& DTM, TInit *UV, long novalue);
//void write_esriascii(char *name, short type, GeoMatrix<long>& DTM, T_INIT *UV, long novalue);
void write_esriascii(char *name, short type, GeoMatrix<long>& DTM, TInit *UV, long novalue);
void write_esriascii(std::string name, short type, GeoMatrix<long>& DTM, TInit *UV, long novalue);
//void write_esriascii_vector(char *name, short type, DOUBLEVECTOR *DTM, long **j, long nr, long nc, T_INIT *UV, long novalue);
  void write_esriascii_vector(char *name, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue);
  void write_esriascii_vector(std::string name, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue);
/*==================================coppied functions from import_ascii.c============*/

#define max_figures 30

//double *read_esriascii(double *header, double novalue, char *name);
void error_message(short format, long n, long n1, long n2, long n3, char *name);

/*===================functions copied from geomorphology.0875.c======================*/
//#define	Pi 3.14159265358979			/* P greco */

void nablaquadro_mask(DOUBLEMATRIX *Z0,SHORTMATRIX *curv,DOUBLEVECTOR *U,DOUBLEVECTOR *V);
void nablaquadro_mask(GeoMatrix<double>& Z0,GeoMatrix<short>& curv,GeoVector<double>& U, GeoVector<double>& V);

void curvature(double deltax, double deltay, DOUBLEMATRIX *topo, DOUBLEMATRIX *c1, DOUBLEMATRIX *c2, DOUBLEMATRIX *c3, DOUBLEMATRIX *c4,long undef);
void curvature(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& c1, GeoMatrix<double>& c2, GeoMatrix<double>& c3, GeoMatrix<double>& c4, long undef);

short is_boundary(long r, long c, GeoMatrix<double>& dem, long novalue);
//long row(double N, long nrows, T_INIT *UV, long novalue);
long row(double N, long nrows, TInit *UV, long novalue);
//long col(double E, long ncols, T_INIT *UV, long novalue);
long col(double E, long ncols, TInit *UV, long novalue);

/*===========copied function from init.h========*/
void initmatrix(double val, DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue);
void initmatrix(double val, GeoMatrix<double>& destination, GeoMatrix<double>& origin, double novalue);

/* ===========copied data from extensions.h==========*/
#define ascii_grass ".grass"
#define ascii_esri ".asc"
#define textfile ".txt"

#endif
