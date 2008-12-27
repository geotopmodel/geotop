/*
 *! GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop Version 0.9375 Lavagna

Copyright, 2008 Riccardo Rigon, Emanuele Cordano, Stefano Endrizzi,  Giacomo Bertoldi, Matteo Dall'Amico, Davide Tamanini, Fabrizio Zanotti, Silvia Simoni, Marco Pegoretti and Paolo Verardo following the footsteps of the FREE SOFTWARE FOUNDATION
LINCENSE

This file is part of GEOtop.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *
 *
 * read_command_line.h
 *
 *  Created on: Jul 19, 2008
 *      Author: ecor
 */

#ifndef READ_COMMAND_LINE_H_
#define READ_COMMAND_LINE_H_


#endif /* READ_COMMAND_LINE_H_ */

char *read_option_string(int argc,char *argv[], char *option_f,char *no_option_argument,short print);

double read_option_double(int argc,char *argv[], char *option_f,char *no_option_argument,double default_value,short print);

int read_flag(int argc,char *argv[],char *flag,short print);

/* debugging options */
#define PRINT_FLAG   "-print"  /*!< flag which prints possible warning message	*/
#define PRINT_ACTIVATED  read_flag(argc,argv,PRINT_FLAG,0)
