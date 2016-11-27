#ifndef TABS_H
#define TABS_H

#include <string>
#include <vector>

/*----------------------------------------------------------------------------------------------------------*/

long count_lines(std::string meteo_file_name, long comment_char, long sep_char);

double **read_txt_matrix(std::string filename, long comment_char, long sep_char,
                         std::vector<std::string> Col_Descr, long ncolsCol_Descr,
                         long *nlines);

double **read_txt_matrix_2(std::string filename, long comment_char,
                           long sep_char, long ncolsCol_Descr, long *nlines);

/*----------------------------------------------------------------------------------------------------------*/
#endif
