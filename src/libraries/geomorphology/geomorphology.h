#ifndef _LIBRARIES_GEOMORPHOLOGY_GEOMORPHOLOGY_H
#define _LIBRARIES_GEOMORPHOLOGY_GEOMORPHOLOGY_H


void find_slope(double deltax, double deltay, Matrix<double> *topo, Matrix<double> *dzdx, Matrix<double> *dzdy,
                long undef);

Matrix<double> * find_max_slope(Matrix<double> *topo, Matrix<double> *dzdx, Matrix<double> *dzdy, long undef);

Matrix<double> * find_aspect(Matrix<double> *topo, Matrix<double> *dzdx, Matrix<double> *dzdy, long undef);

void curvature(double deltax, double deltay, Matrix<double> *topo,
               Matrix<double> *c1, Matrix<double> *c2, Matrix<double> *c3, Matrix<double> *c4,
               long undef);

void topofilter(Matrix<double> *Zin, Matrix<double> *Zout, long novalue, long n);

void order_values(Vector<double>* list, long n);

void multipass_topofilter(long ntimes, Matrix<double> *Zin, Matrix<double> *Zout, long novalue, long n);

short is_boundary(long r, long c, Matrix<double> *dem, long novalue);

void find_min_max(Matrix<double> *M, long novalue, double *max, double *min);

long row(double N, long nrows, T_INIT *UV, long novalue);

long col(double E, long ncols, T_INIT *UV, long novalue);

double topo_from_origin(double **topo, double E, double N, long ncols, long nrows, T_INIT *UV, long novalue);

double interp_value(double E, double N, Matrix<double> *M, Matrix<double> *Z, T_INIT *UV, long novalue);






#endif
