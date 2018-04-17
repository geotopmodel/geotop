void find_slope(double deltax, double deltay, DOUBLEMATRIX *topo,
                DOUBLEMATRIX *dzdx, DOUBLEMATRIX *dzdy, long undef);
DOUBLEMATRIX *find_max_slope(DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx,
                             DOUBLEMATRIX *dzdy, long undef);
DOUBLEMATRIX *find_aspect(DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx,
                          DOUBLEMATRIX *dzdy, long undef);
void curvature(double deltax, double deltay, DOUBLEMATRIX *topo,
               DOUBLEMATRIX *c1, DOUBLEMATRIX *c2, DOUBLEMATRIX *c3, DOUBLEMATRIX *c4,
               long undef);
void topofilter(DOUBLEMATRIX *Zin, DOUBLEMATRIX *Zout, long novalue, long n);
void order_values(Vector<double>* list, long n);
void multipass_topofilter(long ntimes, DOUBLEMATRIX *Zin, DOUBLEMATRIX *Zout,
                          long novalue, long n);
short is_boundary(long r, long c, DOUBLEMATRIX *dem, long novalue);
void find_min_max(DOUBLEMATRIX *M, long novalue, double *max, double *min);
long row(double N, long nrows, T_INIT *UV, long novalue);
long col(double E, long ncols, T_INIT *UV, long novalue);
double topo_from_origin(double **topo, double E, double N, long ncols,
                        long nrows, T_INIT *UV, long novalue);
double interp_value(double E, double N, DOUBLEMATRIX *M, DOUBLEMATRIX *Z,
                    T_INIT *UV, long novalue);





