


#include<cstdio>


FILE *t_fopen(const char *, const char *);


FILE *t_fclose(FILE *stream);


long write_shortmatrix_elements(FILE *, SHORTMATRIX *, long);


long write_intmatrix_elements(FILE *, INTMATRIX *, long);


long write_longmatrix_elements(FILE *, LONGMATRIX *, long);


long write_floatmatrix_elements(FILE *, FLOATMATRIX *, long);


long write_doublematrix_elements(FILE *, DOUBLEMATRIX *, long);


char *join_strings(char *, char *);


char *get_workingdirectory(void);


void write_floatarray_elements(FILE *, FLOATVECTOR *, long);


void write_doublearray_elements(FILE *, Vector<double> *, long);


