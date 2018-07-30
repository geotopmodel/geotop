
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 1.225-15

 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//BASE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

Matrix<short> * copyshort_doublematrix(Matrix<double> *M);

Matrix<long> * copylong_doublematrix(Matrix<double> *M);

Matrix<double> * copydouble_longmatrix(Matrix<long> *L);

Matrix<double> *copydoublematrix_const(double c0, Matrix<double> *Mref, double NOVALUE);

DOUBLETENSOR *build_frommatrix(DOUBLEMATRIX *M, long l, long lmax);

void write_frommatrix(long l, DOUBLEMATRIX *M, DOUBLETENSOR *T);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//UTILITITY subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_suffix(char *suffix, long i, short start);

short existing_file(char *name);
short existing_file_wext(char *name, char *extension);
short existing_file_woext(char *name);

char *namefile_i(char *name, long i);
char *namefile_i_we(char *name, long i);
char *namefile_i_we2(char *name, long i);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//READ subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

Matrix<double> * read_map(short a, char *filename, Matrix<double> *Mref,
                          T_INIT *UVref, double no_value);

std::unique_ptr<Vector<double>> read_map_vector(short type, char *namefile, Matrix<double> *mask,
                                                T_INIT *grid, double no_value, Matrix<long> *rc);

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//WRITE subroutines
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

void write_map(char *filename, short type, short format, Matrix<double> *M,
               T_INIT *UV, long novalue);

void write_map_vector(char *filename, short type, short format,
                      Vector<double> *V, T_INIT *UV, long novalue, long **j, long nr, long nc);

void write_tensorseries(short a, long l, long i, char *filename, short type,
                        short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);

void write_tensorseries_vector(short a, long l, long i, char *filename,
                               short type, short format, Matrix<double> *T, T_INIT *UV, long novalue, long **J,
                               long nr, long nc);

void rename_tensorseries(short a, long l, long i, char *filename);

void rename_map(char *filename);

void write_tensorseries2(char *suf, long l, char *filename, short type,
                         short format, DOUBLETENSOR *T, T_INIT *UV, long novalue);

void write_tensorseries2_vector(char *suf, long l, char *filename, short type,
                                short format, Matrix<double> *T, T_INIT *UV, long novalue, long **J, long nr,
                                long nc);

void write_tensorseries3_vector(char *suffix, char *filename, short type,
                                short format, Matrix<double> *T, T_INIT *UV, long novalue, long **J, long nr,
                                long nc);
