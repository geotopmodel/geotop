/*
 * get_filenames.c
 *
 *  Created on: Dec 26, 2008
 *      Author: ecor
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
	short ik;

	if (print==1) printf("\nWORKING DIRECTORY: %s",WORKING_DIRECTORY);
	/* READ FILENAMES with KeyPalette */
	keyfile=join_strings(WORKING_DIRECTORY,join_strings(program,INIT));
	inptsfile=join_strings(WORKING_DIRECTORY,join_strings(program,INPTS));
	/* GET KEYWORD FOR EACH FILES */
	fkey=t_fopen(keyfile,"r");
	ik=read_index(fkey,print); /* read file keywords string */
	keywords_file=read_keywords(fkey,print); /* read  keywords for each file */
	t_fclose(fkey);
	write_keywords(keywords_file);



	/* GET FILENAMES */
	finpts=t_fopen(inptsfile,"r");
	ik=read_index(finpts,print);
	files=join_path_to_stringbin(WORKING_DIRECTORY,read_names(finpts,keywords_file->names,NOFILE_NAME,print),"/");


	write_read_filenames(files,keywords_file->comments);

	t_fclose(finpts);

	free_keywords(keywords_file);

	return files;
}
