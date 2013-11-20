
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
    

/*----------------------------------------------------------------------------------------------------------*/

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

char **ReadHeader(FILE *f, char *filename, long *num_cols);

long *ColumnCoder(char *filename, char **ColDescr, long max_num_cols, char **header, long num_cols_header, FILE *flog);

long count_lines(char *meteo_file_name, long comment_char, long sep_char);

double **read_datamatrix(FILE *f, char *filename, long comment_char, long sep_char, long number_lines, long components_header);

double **read_txt_matrix(char *filename, long comment_char, long sep_char, char **Col_Descr, long ncolsCol_Descr, long *nlines, FILE *flog);

double **read_txt_matrix_2(char *filename, long comment_char, long sep_char, long ncolsCol_Descr, long *nlines);

/*----------------------------------------------------------------------------------------------------------*/

char *assign_string(char *a);

/*----------------------------------------------------------------------------------------------------------*/

void convert_string_in_lower_case(char *s);
