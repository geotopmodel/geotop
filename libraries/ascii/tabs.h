#ifndef TABS_H
#define TABS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "tabs.h"
#include "../fluidturtle/turtle.h"

#include "../../geotop/constants.h"

#define max_components 200
#define max_string_length 200

extern long number_novalue, number_absent;
extern char *string_novalue;

short readline_par(FILE *f, long comment_char, long sepfield_char, long sepvect_char, long maxcharstring, long maxnumvect, long *key, 
				   long *keylength, long *string, long *stringlength, double *number, long *numberlength, short *endoffile);
double find_number(long *vector, long lengthvector);
char *find_string(long *vector, long lengthvector);
double *find_number_vector(double *vector, long lengthvector);
long *find_string_int(long *vector, long lengthvector);

/*----------------------------------------------------------------------------------------------------------*/

short readline(FILE *f, long comment_char, long sep_char, long **string, long *string_length, long *components, long maxcomponents, long maxstringlength, short *endoffile);
char **readline_of_strings(FILE *f, long comment_char, long sep_char, long *components, short *endoffile, short *success);
double *readline_of_numbers(FILE *f, long comment_char, long sep_char, long *components, short *endoffile, short *success);
//char **ReadHeader(FILE *f, char *filename, long *num_cols);
char **ReadHeader(FILE *f, std::string filename, long *num_cols);
//long *ColumnCoder(char *filename, char **ColDescr, long max_num_cols, char **header, long num_cols_header, FILE *flog);
long *ColumnCoder(std::string filename, char **ColDescr, long max_num_cols, char **header, long num_cols_header, FILE *flog);
//long count_lines(char *meteo_file_name, long comment_char, long sep_char);
long count_lines(std::string meteo_file_name, long comment_char, long sep_char);
//double **read_datamatrix(FILE *f, char *filename, long comment_char, long sep_char, long number_lines, long components_header);
double **read_datamatrix(FILE *f, std::string filename, long comment_char, long sep_char, long number_lines, long components_header);
//double **read_txt_matrix(char *filename, long comment_char, long sep_char, char **Col_Descr, long ncolsCol_Descr, long *nlines, FILE *flog);
double **read_txt_matrix(std::string filename, long comment_char, long sep_char, char **Col_Descr, long ncolsCol_Descr, long *nlines, FILE *flog);
//double **read_txt_matrix_2(char *filename, long comment_char, long sep_char, long ncolsCol_Descr, long *nlines);
double **read_txt_matrix_2(std::string filename, long comment_char, long sep_char, long ncolsCol_Descr, long *nlines);
/*----------------------------------------------------------------------------------------------------------*/

char *assign_string(char *a);

/*----------------------------------------------------------------------------------------------------------*/

void convert_string_in_lower_case(char *s);

#endif
