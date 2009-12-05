
/* KeyPalette MANAGES  THE I/O FILES  OF A MODEL
KeyPalette Version 0.9375 KMackenzie

file key.palette.c

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



#include "tensor3D.h"
//#include "networks.h"
#include "t_utilities.h"
#include "key.palette.h"

#define REFK "_to_"

KEYWORDS *read_keywords (FILE *init, short print) {

/*!

 *  KEYWORDS *read_keywords (char* namefile_init, short print)
 *
 * \param init   -(* FILE) File pointer
 * \param print  -(short) print as
 *
 *
 *
 * \return a KEYWORDS struct which contains two stringbins, one for the keyword names, the other for the referenced comments.
 *
 *	\brief it reads a KEYWORDS struct from a string array
 *
 * \author Emanuele Cordano
 *
 * \date 30 January 2008
 *
 *
 * \relate FLUIDTURLE_STATICLIBRARY.read_stringarray
 */


KEYWORDS *keywords;



STRINGBIN *stringarray;

LONGVECTOR *vname, *vcomment;

long nrows,kj;


keywords=(KEYWORDS *)malloc(sizeof(KEYWORDS));
if (!keywords) t_error("Keywords struct was not allocated");


stringarray=read_stringarray(init,print);


if (stringarray->index->nh%2==0) {
	nrows=stringarray->index->nh/2;

}else{
	t_error("Missing string (keywords or comment!!) in the *.init or *.inpts file   ");
	return NULL;

}

vname=new_longvector(nrows);
vcomment=new_longvector(nrows);

for (kj=vname->nl;kj<=nrows;kj++){
	vname->element[kj]=stringarray->index->element[2*kj-1];
	vcomment->element[kj]=stringarray->index->element[2*kj];
//	printf("\n vcom= %d \n",vcomment->element[kj]);
}

keywords->names=new_stringbin(vname);
keywords->comments=new_stringbin(vcomment);


for (kj=1;kj<=nrows;kj++){

	strcpy(keywords->comments->element[kj]+1,stringarray->element[2*kj]+1);
	strcpy(keywords->names->element[kj]+1,stringarray->element[2*kj-1]+1);


}

free_stringbin(stringarray);
free_longvector(vname);
free_longvector(vcomment);


if (!keywords) t_error(" Function read_keywords did not work !!");
return keywords;






}


KEYWORDS_LIST *read_keywords_list (char *filename, short print){

/*!
	 * KEYWORDS_LIST *read_keywords_array (char *filename, short print)
	 *
	 * \param filename   -(* char) a filename
     * \param print      -(short)  a short int used by FluidTurle functions
     *
     *
	 *
     * \return a KEYWORDS_LIST struct which contains two or more KEYWORD structs.
     *
     * \brief it reads a KEYWORDS_LIST struct from a
	 *
	 *
     * \author Emanuele Cordano
	 *
	 * \date 13 Fabruary 2008

	* \details The function reads a ascii file written in FluidTurtle Format
	*
	* \relate FLUIDTURLE_STATICLIBRARY.read_stringarray






	 *
	 */

	FILE *fd;
	KEYWORDS_LIST* array;

	long index,j,nl,nh;


	fd=t_fopen(filename,"r");
    index=read_index(fd,print); /*!<   read file keywords string */
    array=(KEYWORDS_LIST *)malloc(sizeof(KEYWORDS_LIST));
    if (!array) t_error("Keywords_list struct was not allocated");

    nl=NL;
	nh=index;
    array->element=(KEYWORDS **)malloc((size_t)((nh-nl+1+NR_END)*sizeof(KEYWORDS *)));
    array->isdynamic=isDynamic;
    array->nh=nh;
    array->nl=nl;

    for (j=nl;j<=nh;j++){
    	array->element[j]=read_keywords(fd,print);
    }

	return array;


	t_fclose(fd);

}

void write_keywords(KEYWORDS *keywords){

/*!
 * void write_keywords(KEYWORDS *keywords)
 *
 * \date 13 February 2006
 *
 * \author Emanuele Cordano
 *
 * \brief This function writes all the elements of a KEYWORD struct. This function was made to test all the functions contains in this .c file.

 */


	long kcount;

	printf ("\n HERE ARE WRITTEN THE FOLLOWING KEYWORDS  :  \n ");

    for (kcount=keywords->names->index->nl;kcount<=keywords->names->index->nh;kcount++){

		  printf ("\nkeyword n %ld : %s  as  %s \n ",kcount,keywords->names->element[kcount]+1,keywords->comments->element[kcount]+1);

    }

	 printf("\n");


}




void write_read_filenames(STRINGBIN *filenames, STRINGBIN *comments) {


	/*!

	 *  void write_read_filenames(STRINGBIN *filenames, STRINGBIN *comments)

	 * \author Emanuele Cordano
	 *
	 * \date 8 Fabruary 2008
	 *
	 *
	 *
	*/





	 long kcount;


	 if (filenames->index->nh!=comments->index->nh) {

		 t_error("Error in write_read_filenames: there is no correspondence between comments and keyword files");
	 }else{

		 printf ("THE PROGRAM HAS JUST READ THE FOLLOWING FILENAMES:  \n ");

		 for (kcount=filenames->index->nl;kcount<=filenames->index->nh;kcount++){
	//	 if (strcmp(filenames->element[kcount]+1,NOFILE_NAME))
			 printf ("\nREAD FILENAME: %s as file %ld  i.e. %s",filenames->element[kcount]+1,kcount,comments->element[kcount]+1);

		 }

	 printf("\n");


	 }

 }


KEYWORDS_LIST *reorder_keywords_list(KEYWORDS_LIST *written_palette, STRINGBIN *keywords_model, long jread, char *empty_name) {


/*!   KEYWORD_LIST *reorder_keyword_list(KEYWOWORD_LIST *written_palette, STRINGBIN *keywords_model,long jread)
tempo troppo alto
		 \param written_palette - (* KEYWORDS_LIST) a KEYWORDS_LIST struct written in a scattered way
		 \param keywords_model -  (* STRINGBIN )  a stringbin which contains the model  keywords in the correct order
		 \param jread          - (long) - the number releted to the element for each KEYWORDS struct contained in keywords_model which contains the model keyword
		 \param empty_name     - (* char) the value which correspons to a missing or inactive model component

		 \return a KEYWORDS_LIST for the model keyword and value in the correct order as reported in written_palette


		 \brief It reorder the keywords_model according to the keyword oreder in written_palette.

		 \details The variable keyword_model contains two or more KEYWORDS struct. In each of these the jread-th element is the keyword of the model component (normally jread is assumed equal to 1), the other elements all possible values of the model component. Here the order of the model component is here modifed according to  the order of the keywords in written_palette.

	 * \relate FluidTurtle librariy

	 * \date 13 February 2008
	 *
	 * \author Emanuele Cordano
	 *
	 */




	KEYWORDS_LIST *palette;

	LONGVECTOR *lennames;
	LONGVECTOR *keyreferences;
	LONGVECTOR *vi,*wi;
	long keycnt,namescnt,j;
	char *key_ref, *key_written;

	palette=(KEYWORDS_LIST *)malloc(sizeof(KEYWORDS_LIST));
    if (!palette) t_error("Keywords_list struct palette was not allocated");

    palette->nl=keywords_model->index->nl;
    palette->nh=keywords_model->index->nh;
    palette->element=(KEYWORDS **)malloc((size_t)((keywords_model->index->nh-keywords_model->index->nl+1+NR_END)*sizeof(KEYWORDS *)));
    palette->isdynamic=isDynamic;

    lennames=new_longvector(keywords_model->index->nh);
    keyreferences=new_longvector(keywords_model->index->nh);
	for (keycnt=keywords_model->index->nl;keycnt<=keywords_model->index->nh;keycnt++){

		key_ref=keywords_model->element[keycnt]+1;

		lennames->element[keycnt]=(long)strlen(empty_name);
		keyreferences->element[keycnt]=EMPTY_VALUE;

		for (namescnt=written_palette->nl;namescnt<=written_palette->nh;namescnt++){

			key_written=written_palette->element[namescnt]->names->element[jread]+1;

			if (!strcmp(key_ref,key_written)) keyreferences->element[keycnt]=namescnt;


		}

	}

	for (keycnt=keywords_model->index->nl;keycnt<=keywords_model->index->nh;keycnt++){

		palette->element[keycnt]=(KEYWORDS *)malloc(sizeof(KEYWORDS *));
		if (!palette->element[keycnt]) t_error("Keywords palette->element[j]  was not allocated");

		if (keyreferences->element[keycnt]==EMPTY_VALUE) {
			vi=new_longvector(2);
			wi=new_longvector(2);

			vi->element[1]=(long)strlen(keywords_model->element[keycnt]+1)+1;
			vi->element[2]=(long)strlen(empty_name)+1;
			wi->element[1]=(long)strlen(empty_name)+1;
			wi->element[2]=(long)strlen(empty_name)+1;

			palette->element[keycnt]->names=new_stringbin(vi);
			palette->element[keycnt]->comments=new_stringbin(wi);

			strcpy(palette->element[keycnt]->names->element[1]+1,keywords_model->element[keycnt]+1);
			strcpy(palette->element[keycnt]->comments->element[1]+1,empty_name);

			strcpy(palette->element[keycnt]->names->element[2]+1,empty_name);
			strcpy(palette->element[keycnt]->comments->element[2]+1,empty_name);

			free_longvector(vi);
			free_longvector(wi);

		}else {

			palette->element[keycnt]->names=new_stringbin(written_palette->element[keyreferences->element[keycnt]]->names->index);
			palette->element[keycnt]->comments=new_stringbin(written_palette->element[keyreferences->element[keycnt]]->comments->index);

			for (j=palette->element[keycnt]->names->index->nl;j<=palette->element[keycnt]->names->index->nh;j++){
				strcpy(palette->element[keycnt]->names->element[j]+1,written_palette->element[keyreferences->element[keycnt]]->names->element[j]+1);
			    strcpy(palette->element[keycnt]->comments->element[j]+1,written_palette->element[keyreferences->element[keycnt]]->comments->element[j]+1);

			}

		}

	}

	free_longvector(keyreferences);
	free_longvector(lennames);



	return palette;


	}



KEYWORDS_LIST *read_and_reorder_keywords_list(char *filename, STRINGBIN *keywords_model, long jread, char *empty_name,short print){


/*!   KEYWORDS_LIST *read_&_reorder_keywords_list(char *filename, STRINGBIN *keywords_model, long jread, char *empty_name,short print)
	 *

	 * \param filename - (* char) a filemane used by read_keywords_list fuction
	 * \param keywords_model - (* STRINGBIN) a string array containg the keyword of the model component
	 * \param jread - (* long) parameteres jread utized in reorder_keyword_list function
	   \param empty_name - (* char) name correspondis to a mising or inactive model component as utilized in reorder_keyword_list function
	   \param print      -(short)  a short int used by FluidTurle library

	   \brief It reads a KEYWORDS_LIST containing the model component values using read_keywords_list function and reorders it according to the STRINGBIN keywords_model using the reorder_keywords_list.
	 *
		\relate read_keywords_list, reorder_keywords_list, FluidTurtle Library

	   \author Emanuele Cordano
	 *
	 * \date 15 February 2008
	 *
	 */



	KEYWORDS_LIST *keylist1, *keylist2;

	keylist1=read_keywords_list(filename,print);

	keylist2=reorder_keywords_list(keylist1,keywords_model,jread,empty_name);

	free_keywords_list(keylist1);

	return keylist2;


}


void free_keywords_list(KEYWORDS_LIST* keywords_list){


	/*!
		 * void free_keywords_array(KEYWORDS_LIST* keywords_list
		 *

		 \param keywords_list - (KEYWORD_LIST *) a generic variable
		\brief It frees the mamory allovcated by keywords_list

		 * \author Emanuele Cordano
		 *
		 * \date 13 Fabruary 2008
		 *
		 */




	long jj;
	for (jj=keywords_list->nl;jj<=keywords_list->nh;jj++) {
		free_keywords(keywords_list->element[jj]);
	}

	free(keywords_list);
}







void free_keywords(KEYWORDS *keywords) {

/*!
	 * void free_keywords(KEYWORDS *keywords)
	 *

	  \param keywords - (KEYWORDS *) a generic variable

	   \brief It frees the mamory allovcated by keywords
	 * \author Emanuele Cordano
	 *
	 * \date 1 Fabruary 2008
	 *




	 */
	free_stringbin(keywords->names);
	free_stringbin(keywords->comments);

	free(keywords);

}




 STRINGBIN *read_names(FILE *fd, STRINGBIN *keywords, char *empty_name,short print)	{



 /*!<
	 	 *  STRINGBIN *read_names(FILE *fd, STRINGBIN *keywords, char *empty_name,short print)	 	 *
	 	 *
		 * \param fd - (* FILE) file pointer
		 * \param keywords - (* STRINGBIN) the string array containig the file keywords
		 * \param empty_name (* char) - a string indicating a missing file
		 *
		 * \return a STRINGNIN containing the filenames order according to STRINGBIN *keywords.
		 *
		 * \brief It reads a string array containg a series of filenames and reorders it according to STRINGBIN *keywords.

		 \details The reed string array is like a KEYWORD struct: filename1, KEYWWORD1, filename22, KEYWORD22. The filenames can be written in a scatterdc order and must be followed by the keyword of the filetype ther are referred to .
		 *
		 \relate FluidTurtle Library,

		 \author Emanuele Cordano
	 	 *
		 *
	 	 * \date 6 Fabruary 2008
	 	 *




	 	 */




		 KEYWORDS *written_names;
	 /*   WARNING: the comment of a keyword a generic comment,
	  * the comment of a filename is a keyword !!! */
	 LONGVECTOR *lennames;
	 LONGVECTOR *keyreferences;
	 char *empty_name_extended;






	 STRINGBIN *names;

	 long keycnt,namescnt;

	 char *key_written, *key_ref;
	 char *empty_name_temp;
	 written_names=read_keywords(fd,print);

	 lennames=new_longvector(keywords->index->nh);
	 keyreferences=new_longvector(keywords->index->nh);
	 empty_name_temp=join_strings(empty_name,REFK);
	 for (keycnt=keywords->index->nl; keycnt<=keywords->index->nh;keycnt++) {

		 key_ref=keywords->element[keycnt]+1;
		 empty_name_extended=join_strings(empty_name_temp,keywords->element[keycnt]+1);
		 lennames->element[keycnt]=(long)strlen(empty_name_extended)+1;
		 free(empty_name_extended);
		 keyreferences->element[keycnt]=EMPTY_VALUE;



		 for (namescnt=written_names->comments->index->nl;namescnt<=written_names->comments->index->nh;namescnt++){

			 key_written=written_names->comments->element[namescnt]+1;
//			 printf ("\n keycnt= %d namescnt= %d  \n ",keycnt,namescnt);
//			 printf(key_ref);
//			 printf("\n and \n");
//			 printf(key_written);
//			 printf("\n");

			 if (!strcmp(key_ref,key_written))   {
//				printf ("\n keycnt= %d namescnt= %d  \n ",keycnt,namescnt);

				lennames->element[keycnt]=written_names->names->index->element[namescnt];
				keyreferences->element[keycnt]=namescnt;

			 }
	     }

	 }
	 names=new_stringbin(lennames);

	 for (keycnt=names->index->nl;keycnt<=names->index->nh;keycnt++){
		 if (keyreferences->element[keycnt]==EMPTY_VALUE) {

			 empty_name_extended=join_strings(empty_name_temp,keywords->element[keycnt]+1);
			 strcpy(names->element[keycnt]+1,empty_name_extended);
			 free(empty_name_extended);
		 }else  {
			 strcpy(names->element[keycnt]+1,written_names->names->element[keyreferences->element[keycnt]]+1);
		 }

	 }

	free_longvector(lennames);
	free_longvector(keyreferences);
	free_keywords(written_names);
	//free(key_ref);
	//free(key_written);

	free(empty_name_temp);
	return names;


 }


 STRINGBIN *join_path_to_stringbin(char *path, STRINGBIN *stringvector,char *no_joinstring) {


 /*
    STRINGBIN *join_path_to_stringbin(char *path, STRINGBIN *stringvector)

 \brief It adds a string to each element of a stringbin
 \author Emanuele Cordano

 \date 5 March 2008

 */

 LONGVECTOR *lenstring;

 STRINGBIN *joinedstring;
 long i;
 char *str;

 lenstring=new_longvector(stringvector->index->nh);
 for (i=lenstring->nl;i<=lenstring->nh;i++) {
#ifdef USE_NETCDF_MAP
	 //calculate exact lenght
	 if (strncmp(stringvector->element[i]+1,"#",1) == 0){
 		//if current line start with # -> netcdf variable definition
		 lenstring->element[i]=(long)strlen(stringvector->element[i]+1)+1;
	}else if (!strcmp(stringvector->element[i]+1,no_joinstring)) {
#else
 	if (!strcmp(stringvector->element[i]+1,no_joinstring)) {
#endif
 		lenstring->element[i]=(long)strlen(stringvector->element[i]+1)+1;
 	} else {
 		str=join_strings(path,stringvector->element[i]+1);
 		lenstring->element[i]=(long)strlen(str)+1;
 		free(str);
 	}
	//printf("%ld-%s\n",i,stringvector->element[i]+1);
 }
 joinedstring=new_stringbin(lenstring);
 for (i=lenstring->nl;i<=lenstring->nh;i++) {
#ifdef USE_NETCDF_MAP
	 if (strncmp(stringvector->element[i]+1,"#",1) == 0){
 		//don't add path to lines starting with # -> netcdf variable definition
	 	strcpy(joinedstring->element[i]+1,stringvector->element[i]+1);
	}else if (!strcmp(stringvector->element[i]+1,no_joinstring)) {
#else
	if (!strcmp(stringvector->element[i]+1,no_joinstring)) {
#endif
 		strcpy(joinedstring->element[i]+1,stringvector->element[i]+1);
 	} else {
 		str=join_strings(path,stringvector->element[i]+1);
 		strcpy(joinedstring->element[i]+1,str);
 		free(str);
 	}
	//printf("%ld-%s-%s\n",i,stringvector->element[i]+1,joinedstring->element[i]+1);
 }
 free_longvector(lenstring);

 return joinedstring;


 }
