
/* GEOTRIVIALUtilities is an under-construction set of  C functions which supports i/o interface
   and other utilities for models like GEOtop
GEOTRIVIALUtilities Version 1.0

file read_command_line.c

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

#include "read_command_line.h"

char *read_option_string(int argc, char *argv[], char *option_f, char *no_option_argument, short print)
{
    /*
     * \author Emanuele Cordano
     * \date July 2008
     *
     *\param (int) - argc
     *\param (char *) - argv
     *\param (char *) - option_f - a string related to the flag which can appear in the command
     *\param (char *) - no_option_f - a string returned by the routine in case of MISSING option_field
     *\param (short) -  print
     *
     *\return Returns the string followed by the string option_f in the command line
     *
     */
    int i;

    for (i = 1; i < argc; i++)                    /* Skip argv[0] (program name). */
    {
        /*
         * Use the 'strcmp' function to compare the argv values
         * to a string of your choice (here, it's the optional
         * argument "-q").  When strcmp returns 0, it means that the
         * two strings are identical.
         */

        if (!strcmp(argv[i], option_f))
        {
            if (print == 1) printf(" \n READ ARGUMENT of %s : %s \n", argv[i], argv[i + 1]);
            return argv[i + 1];
        }
    }
    /* warning message */
    if (print == 1)  printf("\nWARNING: Option %s is set to default value %s", option_f, no_option_argument);

    return no_option_argument;
}


//overloaded function
std::string read_option_string(int argc, char *argv[], std::string option_f, std::string no_option_argument, short print)
{
    /*
     * \author Noori
     * \date July 2012
     *
     *\param (int) - argc
     *\param (char *) - argv
     *\param (char *) - option_f - a string related to the flag which can appear in the command
     *\param (char *) - no_option_f - a string returned by the routine in case of MISSING option_field
     *\param (short) -  print
     *
     *\return Returns the string followed by the string option_f in the command line
     *
     */

    for (int i = 1; i < argc; i++)                /* Skip argv[0] (program name). */
    {
        /*
         * Use the 'strcmp' function to compare the argv values
         * to a string of your choice (here, it's the optional
         * argument "-q").  When strcmp returns 0, it means that the
         * two strings are identical.
         */

        if (!strcmp(argv[i], option_f.c_str()))
        {
            if (print == 1) printf(" \n READ ARGUMENT of %s : %s \n", argv[i], argv[i + 1]);
            return argv[i + 1];
        }
    }
    /* warning message */
    if (print == 1)  printf("\nWARNING: Option %s is set to default value %s", option_f.c_str(), no_option_argument.c_str());

    return no_option_argument;
}


double read_option_double(int argc, char *argv[], char *option_f, char *no_option_argument, double default_value, short print)
{
    /*
     * \author Emanuele Cordano
     * \date July 2008
     *
     *\param (int) - argc
     *\param (char *) - argv
     *\param (char *) - option_f - a string related to the flag which can appear in the command
     *\param (char *) - no_option_f - a string returned by the routine in case of MISSING option_field
     *\param (double)- default_value - the DOUBKLE value applied if option_f is missing
     *\param (short) -  print
     *
     *\return Return the double value followed by the the string option_f in the command line
     *
     */
    int s = NO_SUCCESS;
    double d;
    std::string argument;

    argument = read_option_string(argc, argv, option_f, no_option_argument, print);

    if (strcmp(argument.c_str(), no_option_argument) != 0)
    {
        s = sscanf(argument.c_str(), "%lf", &d);
    }

    if (s == NO_SUCCESS && print == 1)
    {
        printf("\nWARNING: No real value found for %s option, the default value %lf is set  \n", option_f, default_value);
        return default_value;
    }

    return d;

}


int read_flag(int argc, char *argv[], char *flag, short print)
/*
 * \author Emanuele Cordano
 * \date July 2008
 *
 */
{
    int i, s = NO_SUCCESS;

    for (i = 1; i < argc; i++)                    /* Skip argv[0] (program name). */
    {
        /*
         * Use the 'strcmp' function to compare the argv values
         * to a string of your choice (here, it's the optional
         * argument "-q").  When strcmp returns 0, it means that the
         * two strings are identical.
         */

        if (strcmp(argv[i], flag) == 0) s = SUCCESS;

    }

    if (print == 1)
    {
        if (s == NO_SUCCESS) printf("\nWARNING: flag %s is missing and not activated \n ", flag);
        if (s == SUCCESS) printf("\nWARNING: flag %s is activated \n ", flag);

    }
    return s;
}
