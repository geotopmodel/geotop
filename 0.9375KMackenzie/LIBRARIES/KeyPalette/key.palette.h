
/* KeyPalette MANAGES  THE I/O FILES  OF A MODEL
KeyPalette Version 0.9375 KMackenzie

file key.palette.h

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


/*!< \file  key.palette.h */


#define EMPTY_VALUE -99

#define NOFILE_NAME "MISSING FILE"

#define NOMODEL_NAME "OFF"



typedef struct {
	STRINGBIN *names;
	STRINGBIN *comments;
} KEYWORDS;

typedef struct {
	short isdynamic;

	const char * name;
	long nh,nl;
	KEYWORDS **element;

} KEYWORDS_LIST;

KEYWORDS *read_keywords (FILE *init, short print);

void free_keywords(KEYWORDS *keywords);

STRINGBIN *read_names(FILE *fd, STRINGBIN *keywords, char *empty_name, short print);

void write_read_filenames(STRINGBIN *filenames, STRINGBIN *comments);

void write_keywords(KEYWORDS *keywords);

KEYWORDS_LIST *read_keywords_list (char *filename, short print);

KEYWORDS_LIST *reorder_keywords_list(KEYWORDS_LIST *written_palette, STRINGBIN *keywords_model, long jread, char *empty_name);

void free_keywords_list(KEYWORDS_LIST* keywords_list);


KEYWORDS_LIST *read_and_reorder_keywords_list(char *filename, STRINGBIN *keywords_model, long jread, char *empty_name,short print);

STRINGBIN *join_path_to_stringbin(char *path, STRINGBIN *stringvector,char *no_joinstring);
