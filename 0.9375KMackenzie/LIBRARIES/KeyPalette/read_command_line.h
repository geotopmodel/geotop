
/* KeyPalette MANAGES  THE I/O FILES  OF A MODEL
KeyPalette Version 0.9375 KMackenzie

file read_command_line.h

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

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


char *read_option_string(int argc,char *argv[], char *option_f,char *no_option_argument,short print);

double read_option_double(int argc,char *argv[], char *option_f,char *no_option_argument,double default_value,short print);

int read_flag(int argc,char *argv[],char *flag,short print);

/* debugging options */
#define PRINT_FLAG   "-print"  /*!< flag which prints possible warning message	*/
#define PRINT_ACTIVATED  read_flag(argc,argv,PRINT_FLAG,0)
