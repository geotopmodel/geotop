
/* TO DO
 * \author Emanuele Cordano
 */

#ifdef USE_NETCDF_ONGOING

char *copy_stringnames(const char *origin);

int rotate180_y_doublematrix(DOUBLEMATRIX *M);
int rotate180_y_doubletensor(DOUBLETENSOR *M);
int rotate180_y_floatmatrix(FLOATMATRIX *M);



int rotate180_y_longmatrix(LONGMATRIX *M);
int rotate180_y_intmatrix(INTMATRIX *M);
int rotate180_y_shortmatrix(SHORTMATRIX *M);

int invert_order_doublevector(DOUBLEVECTOR *v);
int invert_order_longvector(LONGVECTOR *v);

#endif
