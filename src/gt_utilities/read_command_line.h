
/* GEOTRIVIALUtilities is an under-construction set of  C functions which
supports i/o interface and other utilities for models like GEOtop
GEOTRIVIALUtilities Version 1.0

file read_command_line.h

Copyright (c), 2011 Emanuele Cordano

This file is part of GEOTRIVIALUtilities.
 GEOTRIVIALUtilities is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOTRIVIALUtilities is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License v. 3.0 for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef READ_COMMANDLINE_H
#define READ_COMMANDLINE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <fstream>
#include <math.h>
#include <time.h>
#include <ctype.h>

const int SUCCESS = 1;
const int NO_SUCCESS = 0;

std::string get_workingdirectory();
void t_error(std::string error_text);

char *read_option_string(int argc,
                         char *argv[],
                         char *option_f,
                         char *no_option_argument,
                         short print);

std::string read_option_string(int argc,
                               char *argv[],
                               std::string option_f,
                               std::string no_option_argument,
                               short print);

double read_option_double(int argc,
                          char *argv[],
                          char *option_f,
                          char *no_option_argument,
                          double default_value,
                          short print);

int read_flag(int argc, char *argv[], char *flag, short print);

/* debugging options */
#define PRINT_FLAG "-print" /*!< flag which prints possible warning message */
#define PRINT_ACTIVATED read_flag(argc, argv, PRINT_FLAG, 0)

/* */
#define MISSING_ARGUMENT "missing_argument"
#endif
