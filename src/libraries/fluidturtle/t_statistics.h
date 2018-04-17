
float floatvector_correlation(FLOATVECTOR *,FLOATVECTOR *u,float,float,long,
                              float );
float longvector_correlation(LONGVECTOR *,LONGVECTOR *,float,float,long,long);

/**

Name: _n_moment

Synopsis:
float longmatrix_n_moment(LONGMATRIX *, double ,double , long );
double doublematrix_n_moment(DOUBLEMATRIX *, double ,double  , double );
float floatmatrix_n_moment(FLOATMATRIX *, float ,float , float );
double doublevector_n_moment(Vector<double>* , float ,float , double );
float floatvector_n_moment(FLOATVECTOR *, float ,float , float );
float longvector_n_moment(LONGVECTOR *, float ,float , long );


Description: It calculates the biased n moment of a vector or a matrix of the specified
type

Inputs:  1) the first vector; 2) the second vector; 3) the first vector mean; 4) the second
vector mean; 5) the value that identify "NO DATA" or "MISSING DATA" (a sequence of mesasured data
usually contains points where the data were not detected or monitored)

Return: biased n moment (the division factor is the length of the vector minus the
number of missing data)


Examples: LIBRARIES/APPLICATIONS/DATA_MANIPULATION/hystogram.c

See Also: doublevector_correlation

Authors & date: Paolo D'Odorico, Riccardo Rigon, February 1998.

FILE: LIBRARIES/BASICMATHSTAT/statistics.c, LIBRARIES/BASICMATHSTAT/t_statistics.h



*/

float longmatrix_n_moment(LONGMATRIX *, double,double, long );
double doublematrix_n_moment(DOUBLEMATRIX *, double,double, double );
float floatmatrix_n_moment(FLOATMATRIX *, float,float, float );

float floatvector_n_moment(FLOATVECTOR *, float,float, float );
float longvector_n_moment(LONGVECTOR *, float,float, long );
float floatmatrix_moment(FLOATMATRIX *, float,float, float );

/**

Name:floatmatrix_restricted _n_moment

Synopsis:
float floatmatrix_n_moment(FLOATMATRIX *first,FLOATMATRIX *second, float ,float , float,float );


Description: It calculates the biased n moment of a vector or a matrix of the specified
type as floatmatrix_n_moment but it uses the subset of points that are not novalue2 in a second
nmatrix

Inputs:  1) the first matrix; 2) The second matrix;  3) the first matrix mean; 4) the second
matrix mean; 5) the value that identify "NO DATA" or "MISSING DATA" in the first matrix (a sequence of mesasured data
usually contains points where the data were not detected or monitored); 6) the vector that identify "NO DATA"
in the second matrix

Return: biased n moment (the division factor is the length of the vector minus the
number of missing data)



Authors & date: Paolo D'Odorico, Riccardo Rigon, February 1998.

FILE: LIBRARIES/BASICMATHSTAT/statistics.c, LIBRARIES/BASICMATHSTAT/t_statistics.h

Examples: LIBRARIES/APPLICATIONS/GEOMORPHOLOGY/RIVER_NETWORKS/sum_downstream.c

*/

float floatmatrix_restricted_n_moment(FLOATMATRIX *,FLOATMATRIX *, float,
                                      float, float, float );

/**


*/
double double_n_moment(double *m, long nh,double mean,double  NN,
                       double novalue);
