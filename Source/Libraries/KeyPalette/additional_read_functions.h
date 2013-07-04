
/* KeyPalette MANAGES  THE I/O FILES  OF A MODEL
KeyPalette Version 0.9375 KMackenzie

file additional_read_functions.h

Copyright (c), 2009 Emanuele Cordano and Riccardo Rigon

This file is part of KeyPalette.
 KeyPalette is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    KeyPalette is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License v. 3.0 for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

DOUBLEVECTOR *read_doublearray_from_string(char *argument,char *delimiters,int max_numbers,short print);

LONGVECTOR *read_longarray_from_string(char *argument,char *delimiters,int max_numbers,short print);
