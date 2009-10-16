
/* KeyPalette MANAGES  THE I/O FILES  OF A MODEL
KeyPalette Version 0.9375 KMackenzie

file get_filenames.c

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

#include "turtle.h"
#include "t_io.h"
#include "key.palette.h"
#include "get_filenames.h"
#define INIT ".init"
#define INPTS ".inpts"


#define NOFILE_NAME join_strings(WORKING_DIRECTORY,"MISSING_FILE")

STRINGBIN *get_filenames_from_keys(char *WORKING_DIRECTORY, char *program, short print){
	/*
	 *
	 * \author Emanuele Cordano
	 *
	 * \date Dec 2008
	 *
	 *\param WORKING_DIRECTORY (char *) - working directory where *.init and *.inpts are contained
	 *\param program (char *) - name of the program
	 *\param print (short) - flag for verdose mode
	 *
	 *\return a stringbin with filenames
	 *
	 */

	char *keyfile, *inptsfile;
	KEYWORDS *keywords_file;
	STRINGBIN *files;
	FILE *fkey,*finpts;
	STRINGBIN *names;
	short ik;
	char *pinit,*pinpts, *nofilename;
	nofilename=NOFILE_NAME; /* allocated with joint-strings */
	if (print==1) printf("\nWORKING DIRECTORY: %s",WORKING_DIRECTORY);
	/* READ FILENAMES with KeyPalette */
	pinit=join_strings(program,INIT);
	keyfile=join_strings(WORKING_DIRECTORY,pinit);
	free(pinit);
	pinpts=join_strings(program,INPTS);
	inptsfile=join_strings(WORKING_DIRECTORY,pinpts);
	free(pinpts);
	/* GET KEYWORD FOR EACH FILES */
	fkey=t_fopen(keyfile,"r");
	ik=read_index(fkey,print); /* read file keywords string */
	keywords_file=read_keywords(fkey,print); /* read  keywords for each file */
	t_fclose(fkey);
	if (print==1) write_keywords(keywords_file);



	/* GET FILENAMES */
	finpts=t_fopen(inptsfile,"r");
	ik=read_index(finpts,print);
	names=read_names(finpts,keywords_file->names,nofilename,print);
	files=join_path_to_stringbin(WORKING_DIRECTORY,names,"/");
	free_stringbin(names);

	if (print==1) write_read_filenames(files,keywords_file->comments);

	t_fclose(finpts);
	free(nofilename);
	free(keyfile);
	free(inptsfile);
	free_keywords(keywords_file);

	return files;
}

