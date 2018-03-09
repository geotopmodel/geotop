/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */


#include "tabs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

//#include "../fluidturtle/turtle.h"
#include "../../gt_utilities/read_command_line.h"
#include <boost/algorithm/string.hpp>
#include "../../geotop/inputKeywords.h"

#include <iostream>
#include "../../geotop/global_logger.h"
#include <sstream>

/*============================================================================*/
/*                               Constants                                    */
/*============================================================================*/
const long max_components = 200;
const long max_string_length = 200;

/*============================================================================*/
/*                      Private functions prototypes                          */
/*============================================================================*/
static double find_number(long *vector, long lengthvector);

static std::string find_string(long *vector, long lengthvector);

static std::vector<std::string> ReadHeader(FILE *f, std::string filename,
                                           long *num_cols);

static short readline(FILE *f, long comment_char, long sep_char, long **string,
                      long *string_length, long *components, long maxcomponents,
                      long maxstringlength, short *endoffile);

static std::vector<std::string> readline_of_strings(FILE *f,
                                                    long comment_char,
                                                    long sep_char,
                                                    long *components,
                                                    short *endoffile,
                                                    short *success);

static double* readline_of_numbers(FILE *f, long comment_char,
                                   long sep_char, long *components,
                                   short *endoffile, short *success);

// this function is never defined...
//static long* ColumnCoder(std::string filename, std::vector<std::string> ColDescr,
//                         long max_num_cols, std::vector<std::string> header,
//                         long num_cols_header, FILE *flog);

static double **read_datamatrix(FILE *f, long comment_char, long sep_char,
                                long number_lines, long components_header);

/*============================================================================*/

double find_number(long *vector, long lengthvector)
{

    double N = 0.0, Nexp = 0.0;
    long i, ie, ids, cnt;

    //looking for the E or e
    ie = -1;
    i = 0;

    do
    {
        if(vector[i] == 69 || vector[i] == 101) ie = i;
        i++;
    }
    while (ie == -1 && i < lengthvector);

    if(ie == -1) ie = lengthvector;

    //looking for decimal separator
    ids = -1;
    i = 0;

    do
    {
        if(vector[i] == 46) ids = i;
        i++;
    }
    while (ids == -1 && i < lengthvector);

    if(ids == -1 && ie > -1) ids = ie;
    if(ids == -1) ids = lengthvector;

    //integer part
    cnt = 0;

    for (i = ids - 1; i >= 0; i--)
    {
        if (vector[i] >= 48 && vector[i] <= 57)
        {
            N += (double)(vector[i] - 48) * pow(10.0, (double)cnt);
            cnt++;
        }
    }

    //fractional part
    cnt = -1;

    for (i = ids + 1; i <= ie - 1; i++)
    {
        if (vector[i] >= 48 && vector[i] <= 57)
        {
            N += (double)(vector[i] - 48) * pow(10.0, (double)cnt);
            cnt--;
        }
    }

    //exponential part
    cnt = 0;

    for (i = lengthvector - 1; i >= ie + 1; i--)
    {
        if (vector[i] >= 48 && vector[i] <= 57)
        {
            Nexp += (double)(vector[i] - 48) * pow(10.0, (double)cnt);
            cnt++;
        }
    }

    //sign of the exponential part
    if (ie < lengthvector - 1)
    {
        if (vector[ie + 1] == 45) Nexp *= (-1.);
    }

    //exponential format
    N = N * pow(10.0, Nexp);

    //sign
    if (vector[0] == 45) N *= (-1.);

    return(N);
}


std::string find_string(long *vector, long lengthvector)
{

  std::stringstream stream;
   long i;

   for (i = 0; i < lengthvector; i++)
   {
     stream << char(vector[i]);
   }

   return stream.str();
}


double *find_number_vector(double *vector, long lengthvector)
{

    double *number_vector;
    long i;

    number_vector = (double*)malloc(lengthvector * sizeof(double));

    for (i = 0; i < lengthvector; i++)
    {
        number_vector[i] = vector[i];
    }

    return(number_vector);
}


long *find_string_int(long *vector, long lengthvector)
{

    long *string;
    long i;

    string = (long*)malloc(lengthvector * sizeof(long));

    for (i = 0; i < lengthvector; i++)
    {
        string[i] = vector[i + 1];
    }

    return(string);
}


static short readline(FILE *f, long comment_char, long sep_char, long **string, long *string_length, long *components, long maxcomponents,
                      long maxstringlength, short *endoffile)
{

    long i, j;
    char *c;

    *endoffile = 0;

    //initialization of vector that are already allocated
    for (i = 0; i < maxcomponents; i++)
    {
        string_length[i] = 0;

        for (j = 0; j < maxstringlength; j++)
        {
            string[i][j] = 0;
        }
    }

    //allocate character
    c = (char*)malloc(sizeof(char));

    //read first character
    do
    {
        c[0] = fgetc(f);
        if(c[0] == -1) *endoffile = 1;
    }
    while ( *endoffile == 0 && ( c[0] <= 42 || ( c[0] >= 58 && c[0] <= 64 ) || c[0] == 96 || c[0] >= 123 ) && c[0] != comment_char );

    //end of file reached
    if( *endoffile == 1)
    {

        free(c);
        return -1;

        //first character is the comment tag
    }
    else if( c[0] == comment_char )
    {

        //read until end of line
        do
        {
            c[0] = fgetc(f);
        }
        while (c[0] != 10 && c[0] != -1);

        if(c[0] == -1) *endoffile = 1;

        free(c);
        return -1;

    }
    else
    {

        //read argument
        j = 0;

        do
        {

            i = 0;

            if(j == 0)
            {
                string[j][i] = c[0];
                i++;
            }

            do
            {
                do
                {
                    c[0] = fgetc(f);
                    if(c[0] == -1) *endoffile = 1;
                }
                while ( *endoffile == 0 && ( c[0] <= 42 || ( c[0] >= 58 && c[0] <= 64 ) || c[0] == 96 || c[0] >= 123 ) && c[0] != 10 && c[0] != sep_char);

                if (*endoffile == 0 && c[0] != 10 && c[0] != sep_char)
                {
                    if(i < maxstringlength)
                    {
                        string[j][i] = c[0];
                        i++;
                    }
                }

            }
            while (*endoffile == 0 && c[0] != 10 && c[0] != sep_char);

            string_length[j] = i;

            j++;

        }
        while( *endoffile == 0 && c[0] != 10 && j < maxcomponents );

        *components = j;

        free(c);
        return 1;
    }
}


static std::vector<std::string> readline_of_strings(FILE *f, long comment_char, long sep_char, long *components, short *endoffile, short *success)
{

    long i, n;
    long **string, *string_length;
    std::vector<std::string> line_of_strings;

    n = max_components;
    string_length = (long*)malloc(n * sizeof(long));
    string = (long**)malloc(n * sizeof(long*));

    n = max_string_length;
    for (i = 0; i < max_components; i++)
    {
        string[i] = (long*)malloc(n * sizeof(long));
    }

    *success = readline(f, comment_char, sep_char, string, string_length, components, max_components, max_string_length, endoffile);

    if(*success == 1)
    {
        for (i = 0; i < (*components); i++)
        {
            line_of_strings.push_back(find_string(string[i], string_length[i]));
        }

    }

    for (i = 0; i < max_components; i++)
    {
        free(string[i]);
    }

    free(string_length);
    free(string);

    return line_of_strings;
}


static double *readline_of_numbers(FILE *f, long comment_char, long sep_char, long *components, short *endoffile, short *success)
{

    long i, n;
    long **string, *string_length;
    double *line_of_numbers = NULL;

    n = max_components;
    string_length = (long*)malloc(n * sizeof(long));
    string = (long**)malloc(n * sizeof(long*));

    n = max_string_length;
    for (i = 0; i < max_components; i++)
    {
        string[i] = (long*)malloc(n * sizeof(long));
    }

    *success = readline(f, comment_char, sep_char, string, string_length, components, max_components, max_string_length, endoffile);

    if(*success == 1)
    {

        line_of_numbers = (double*)malloc(*components * sizeof(double));

        for (i = 0; i < (*components); i++)
        {
            line_of_numbers[i] = find_number(string[i], string_length[i]);
        }

    }

    for (i = 0; i < max_components; i++)
    {
        free(string[i]);
    }

    free(string_length);
    free(string);

    return line_of_numbers;
}


static std::vector<std::string> ReadHeader(FILE *f, std::string filename, long *num_cols)
{

    short endoffile, success;
    std::vector<std::string> Header;

    do
    {
        Header = readline_of_strings(f, 33, 44, num_cols, &endoffile, &success);
    }
    while (endoffile == 0 && success != 1);

    if (endoffile == 1)
    {
        printf("Error!! File %s contains only the header\n", filename.c_str());
        t_error("Impossible To Continue");
    }

    return Header;

}


static long* ColumnCoder(std::string filename, std::vector<std::string> ColDescr,
                         long max_num_cols, std::vector<std::string> header,
                         long num_cols_header)
{

    long *coder, i, j;
    std::string lowercaseColDescr;
    geotop::logger::GlobalLogger* lg = geotop::logger::GlobalLogger::getInstance();

    //allocation
    coder = (long*)malloc(max_num_cols * sizeof(long));

    //initialization
    for( i = 0; i < max_num_cols; i++)
    {
        coder[i] = -1;
    }

    for( i = 0; i < max_num_cols; i++)
    {
        lowercaseColDescr = ColDescr[i];
        boost::algorithm::to_lower(lowercaseColDescr);

        for (j = 0; j < num_cols_header; j++)
        {
            boost::algorithm::to_lower(header[j]);

            if (lowercaseColDescr == header[j] && coder[i] == -1 && geotop::input::gStringNoValue != header[j])
            {
              coder[i] = j;
              lg->logf("Column %ld in file %s assigned to %s", j + 1, filename.c_str(), ColDescr[i].c_str());
            }
            else if (lowercaseColDescr == header[j] && coder[i] != -1)
            {
                t_error("Column name DUPLICATED in a comma-separated-value tables!");
            }
        }
    }

    return coder;
}


long count_lines(std::string meteo_file_name, long comment_char, long sep_char)
{
    FILE *f;
    std::vector<std::string> header;
    double *line;
    short success, endoffile;
    long  components, cont;

    f = fopen(meteo_file_name.c_str(), "r");
    if (f == NULL)
    {
        printf("File -0- %s not existing\n", meteo_file_name.c_str());
        t_error("Fatal Error (10)");
    }

    //read header line
    do
    {
        header = readline_of_strings(f, comment_char, sep_char, &components, &endoffile, &success);
    }
    while(success != 1 && endoffile == 0);

    if (endoffile == 1)
    {
        fclose(f);
        return 0;
    }

    //count lines
    cont = 0;
    do
    {
        line = readline_of_numbers(f, comment_char, sep_char, &components, &endoffile, &success);
        if (success == 1)
        {
            free(line);
            cont++;
        }
    }
    while( endoffile == 0);

    fclose(f);
    return cont;
}


static double **read_datamatrix(FILE *f, long comment_char, long sep_char, long number_lines, long components_header)
{

    double **data, *line;
    short endoffile, success;
    long i, cont, components;

    data = (double**)malloc(number_lines * sizeof(double*));

    cont = 0;
    do
    {
        line = readline_of_numbers(f, comment_char, sep_char, &components, &endoffile, &success);

        if (success == 1)
        {
            data[cont] = (double*)malloc(components_header * sizeof(double));
            if (components <= components_header)
            {
                for (i = 0; i < components; i++)
                {
                    data[cont][i] = line[i];
                }
                for (i = components; i < components_header; i++)
                {
                    data[cont][i] = geotop::input::gDoubleNoValue;
                }
            }
            else
            {
                for (i = 0; i < components_header; i++)
                {
                    data[cont][i] = line[i];
                }
            }

            free(line);
            cont++;

        }
    }
    while( cont < number_lines );

    return data;
}


double **read_txt_matrix(std::string filename, long comment_char, long sep_char, std::vector<std::string> Col_Descr, long ncolsCol_Descr, long *nlines)
{

    /*Read header, and create a **double with the same columns as the header. Then fill with geotop::input::gDoubleAbsent the columns
     missing with respect to Col_Descr*/

    FILE *f;
    std::vector<std::string> Header;
    double **Data, **Dataout;
    long *Coder, ncols, i, j;

    *nlines = count_lines(filename, comment_char, sep_char);

    f = fopen(filename.c_str(), "r");
    if (f == NULL)
    {
        printf("File -1- %s not existing\n", filename.c_str());
        t_error("Fatale Error (11)");
    }
    Header = ReadHeader(f, filename, &ncols);

    Coder = ColumnCoder(filename, Col_Descr, ncolsCol_Descr, Header, ncols);
    
    Data = read_datamatrix(f, comment_char, sep_char, *nlines, ncols);
    fclose(f);

    Dataout = (double**)malloc((*nlines) * sizeof(double*));
    for (i = 0; i < (*nlines); i++)
    {
        Dataout[i] = (double*)malloc(ncolsCol_Descr * sizeof(double));
        for (j = 0; j < ncolsCol_Descr; j++)
        {
            if (Coder[j] != -1)
            {
                Dataout[i][j] = Data[i][Coder[j]];
            }
            else
            {
                Dataout[i][j] = geotop::input::gDoubleAbsent;
            }
        }
        free(Data[i]);
    }
    free(Data);
    free(Coder);

    return Dataout;

}


double **read_txt_matrix_2(std::string filename, long comment_char, long sep_char, long ncolsCol_Descr, long *nlines)
{

    /*Read header, and create a **double with the same columns as Col_Descr, filling with geotop::input::gDoubleNoValue the missing columns*/

    FILE *f;
    std::vector<std::string> Header;
    double **Dataout;
    long ncols;

    *nlines = count_lines(filename, comment_char, sep_char);

    f = fopen(filename.c_str(), "r");
    if (f == NULL)
    {
        printf("File -2- %s not existing\n", filename.c_str());
        t_error("Fatale Error (13)");
    }

    Header = ReadHeader(f, filename, &ncols);
    Dataout = read_datamatrix(f, comment_char, sep_char, *nlines, ncolsCol_Descr);
    fclose(f);

    return Dataout;

}

