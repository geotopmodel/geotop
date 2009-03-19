
/* KeyPalette MANAGES  THE I/O FILES  OF A MODEL
KeyPalette Version 0.9375 KMackenzie

file read_command_line.c

Copyright, 2009 Emanuele Cordano and Riccardo Rigon

This file is part of KeyPalette.
 KeyPalette is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    KeyPalette is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



//#include <cstdio>
#include <stdio.h>
#include <string.h>
#include "read_command_line.h"

#define SUCCESS 1
#define NO_SUCCESS 0




char *read_option_string(int argc,char *argv[], char *option_f,char *no_option_argument,short print)
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
	int i,s;



	s=NO_SUCCESS;

	for (i = 1; i < argc; i++)  /* Skip argv[0] (program name). */
		{
	        /*
	         * Use the 'strcmp' function to compare the argv values
	         * to a string of your choice (here, it's the optional
	         * argument "-q").  When strcmp returns 0, it means that the
	         * two strings are identical.
	         */

	        if (!strcmp(argv[i],option_f)) {
	        	if (print==1) printf(" \n READ ARGUMENT of %s : %s \n",argv[i],argv[i+1]);
	        	return argv[i+1];
	        }
		}
    /* warning message */
	if (print==1)  printf("\nWARNING: Option %s is set to default value %s",option_f,no_option_argument);




	return no_option_argument;
}


double read_option_double(int argc,char *argv[], char *option_f, char *no_option_argument, double default_value,short print)
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
	int s;
	double d;
	char *argument;


	argument=read_option_string(argc,argv,option_f,no_option_argument,print);

	s=sscanf(argument,"%lf",&d);
	if (s==0 && print==1) {
			printf("\nWARNING: No real value found for %s option, the default value %lf is set  ",option_f,default_value);
			return default_value;
	}
//	d=default_value;
	return d;
}

int read_flag(int argc,char *argv[],char *flag,short print)
/*
	 * \author Emanuele Cordano
	 * \date July 2008
	 *
	 */
{
	int s,i;
	s=NO_SUCCESS;
	for (i = 1; i < argc; i++)  /* Skip argv[0] (program name). */
		{
		        /*
		         * Use the 'strcmp' function to compare the argv values
		         * to a string of your choice (here, it's the optional
		         * argument "-q").  When strcmp returns 0, it means that the
		         * two strings are identical.
		         */

		     if (strcmp(argv[i],flag) == 0) s=SUCCESS;

		}
	if (print==1) {
		if (s==NO_SUCCESS) printf("\nWARNING: flag %s is missing and not activated \n ",flag);
		if (s==SUCCESS) printf("\nWARNING: flag %s is activated \n ",flag);

	}
	return s;
}










