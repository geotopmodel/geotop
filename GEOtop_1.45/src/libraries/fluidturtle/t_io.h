/*-------------------------------------------------------------------------------



 You will find for each method (function, routine) acting on a type a replica acting



 on the other types. To indicate collectively the replicas we will use "*". Thus



 in that context,   * substitutes:



 shortvector, intvector, longvector, floatvector, doublevector, charvector,



 shortmatrix, intmatrix, longmatrix, floatmatrix, doublematrix, charmatrix,



 shortbin, intbin, longbin,doublebin ,stringbin.



 "**" instead substitutes:



 shortarray,intarray,longarray,chararray,floatarray,doublearray,stringarray





 Any of the routines is commented in a standard way as follows:



 ---------------------------------------------------------------------------------*/

/**







 Name: t_fopen, t_fclose



 Synopsis: FILE *t_fopen(const char * ,const char *); FILE *t_fclose(FILE * stream);



 Version: 1.0



 Description:



 Safe open  and close file. t_fopen() used in writing mode checks if an file with the existing name already exists.



 In this case it moves it to filename.old. Future inplementations will implement the check of the mode of opening.



 



 Authors & Date:  Riccardo Rigon, 1996



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: name of the file and mode of opening





 Return: a pointer to the opened file or a NULL pointer in the case of t_fclose









 See Also: fopen, close, t_error



 Keywords: streams



 Examples: 1.example.c





 */

FILE *t_fopen(const char *, const char *);

FILE *t_fclose(FILE * stream);

/**







 Name: read_, write__elements, binarywrite__elements.



 Version: 1.0



 Synopsis:



 long read_shortvector_elements(FILE *  ,SHORTVECTOR *,char *);



 long write_shortvector_elements(FILE *,SHORTVECTOR *,long );



 long binarywrite_shortvector_elements(FILE *,SHORTVECTOR *);



 void print_shortvector_elements(SHORTVECTOR *,long );







 Description: For each data type four functions are given: 1) read_* to

 read (both in ascii or binary format) the elements of a vector, matrix or bin,

 from the specified FILE stream; 2) write_* to write it on the specified file;

 3) binarywrite_* to write in binary format and 4) print_* to print it on the screen.

 Using "long" as a return value obviosly limits the size of the vector and matrixes used.

 A better choice would have been size_t but this prevents to have returned

 negative numbers that are used to handle the errors. Thus, use of size_t would require

 a rewriting of the error handling part of the routines. Also, one could observe that, in the

 case of matrixes, the total number of elements (rows times columns) must be limited to

 the maximum of long, say MAX_LONG, in order to have correct return values even if one could theoretically

 allocate matrixes of size MAX_LONG*MAX_LONG.







 Authors & Date: Riccardo Rigon, February 1997,  April 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: read_ : The pointer to the opened file stream FILE *,  the pointer to the structure



 were the data are going to be allocated and the mode in which the data are read ("a"



 for ascii and "b" for binary).



 write_: The pointer to the file were the data are going to be write, the pointer to



 the structure were they are actually allocated, the number of columns of the output.



 binarywrite_*:The pointer to the file were the data are going to be write, the pointer to



 the structure were they are actually allocated.



 print_ : the pointer to the file and the maximum number of data to be printed



 on the same row.







 Return: read_, write_, binarywrite_: a long indicating the number of elements read/written. A negative long is



 returned if some error is encountered in reading/writing. print_ returns a void.



 In the case of reading a STRINGBIN or array, a pointer to STRINGBIN



 is returned.



 Examples: 1.examples.c, 2.examples.c





 References:





 */

long read_shortvector_elements(FILE *, SHORTVECTOR *, char *);

long write_shortvector_elements(FILE *, SHORTVECTOR *, long);

long binarywrite_shortvector_elements(FILE *, SHORTVECTOR *);

void print_shortvector_elements(SHORTVECTOR *, long);

long read_intvector_elements(FILE *, INTVECTOR *, char *);

long write_intvector_elements(FILE *, INTVECTOR *, long);

long binarywrite_intvector_elements(FILE *, INTVECTOR *);

void print_intvector_elements(INTVECTOR *, long);

long read_longvector_elements(FILE *, LONGVECTOR *, char *);

long write_longvector_elements(FILE *, LONGVECTOR *, long);

long binarywrite_longvector_elements(FILE *, LONGVECTOR *);

void print_longvector_elements(LONGVECTOR *, long);

long read_floatvector_elements(FILE *, FLOATVECTOR *, char *);

long write_floatvector_elements(FILE *, FLOATVECTOR *, long);

long binarywrite_floatvector_elements(FILE *, FLOATVECTOR *);

void print_floatvector_elements(FLOATVECTOR *, long);

long read_doublevector_elements(FILE *, DOUBLEVECTOR *, char *);

long write_doublevector_elements(FILE *, DOUBLEVECTOR *, long);

long binarywrite_doublevector_elements(FILE *, DOUBLEVECTOR *);

void print_doublevector_elements(DOUBLEVECTOR *, long);

long read_charvector_elements(FILE *, CHARVECTOR *, char *);

long write_charvector_elements(FILE *, CHARVECTOR *, long);

long binarywrite_charvector_elements(FILE *, CHARVECTOR *);

void print_charvector_elements(CHARVECTOR *, long);

long read_shortmatrix_elements(FILE *, SHORTMATRIX *, char *);

long write_shortmatrix_elements(FILE *, SHORTMATRIX *, long);

long binarywrite_shortmatrix_elements(FILE *, SHORTMATRIX *);

void print_shortmatrix_elements(SHORTMATRIX *, long);

long read_intmatrix_elements(FILE *, INTMATRIX *, char *);

long write_intmatrix_elements(FILE *, INTMATRIX *, long);

long binarywrite_intmatrix_elements(FILE *, INTMATRIX *);

void print_intmatrix_elements(INTMATRIX *, long);

long read_longmatrix_elements(FILE *, LONGMATRIX *, char *);

long write_longmatrix_elements(FILE *, LONGMATRIX *, long);

long binarywrite_longmatrix_elements(FILE *, LONGMATRIX *);

void print_longmatrix_elements(LONGMATRIX *, long);

long read_floatmatrix_elements(FILE *, FLOATMATRIX *, char *);

long write_floatmatrix_elements(FILE *, FLOATMATRIX *, long);

long binarywrite_floatmatrix_elements(FILE *, FLOATMATRIX *);

void print_floatmatrix_elements(FLOATMATRIX *, long);

long read_doublematrix_elements(FILE *, DOUBLEMATRIX *, char *);

long write_doublematrix_elements(FILE *, DOUBLEMATRIX *, long);

long binarywrite_doublematrix_elements(FILE *, DOUBLEMATRIX *);

void print_doublematrix_elements(DOUBLEMATRIX *, long);

long read_intbin_elements(FILE *, INTBIN *, char *);

long write_intbin_elements(FILE *, INTBIN *, long);

long binarywrite_intbin_elements(FILE * output, INTBIN *l);

void print_intbin_elements(INTBIN *, long);

STRINGBIN *read_plane_strings(FILE *, long, long);

long write_stringbin_elements(FILE *, STRINGBIN *, long);

void print_stringbin_elements(STRINGBIN *, long);

long read_shortbin_elements(FILE *, SHORTBIN *, char *);

long write_shortbin_elements(FILE *, SHORTBIN *, long);

long binarywrite_shortbin_elements(FILE *, SHORTBIN *);

void print_shortbin_elements(SHORTBIN *, long);

long read_longbin_elements(FILE *, LONGBIN *, char *);

long write_longbin_elements(FILE *, LONGBIN *, long);

long binarywrite_longbin_elements(FILE *, LONGBIN *);

void print_longbin_elements(LONGBIN *, long);

long read_doublebin_elements(FILE *, DOUBLEBIN *, char *);

long write_doublebin_elements(FILE *, DOUBLEBIN *, long);

long binarywrite_doublebin_elements(FILE *, DOUBLEBIN *);

void print_doublebin_elements(DOUBLEBIN *, long);

/**











 Name: copy_buffer_into_stringbin



 Synopsis: int copy_buffer_into_stringbin(char *,STRINGBIN *);



 Version: 1.0



 Description: In order to read a fixed number of strings whose



 length is unknown,  a buffer is inizialized and filled with the strings



 reallocating each time is necessary a certain amount of memory



 specified in the global variable BUFFER_INCREMENT. A buffer is termined



 by a null character. Each string is separated by another by a whitespace



 (in the sense of iswhitespace() - e.g. Harbison and Steele, 1995). While



 reading them is kept track of the length of each string and the suitable



 vector index for the string bin is inizialized (note: the number of string



 must be known in advance). Finally one has to copy the content of the buffer



 into its stringbin.





 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c





 Authors & Date: Riccardo Rigon, August 1997





 Inputs: read_: A pointer to the buffer and a pointer to the stringbin.

 Return: 1 if successful (otherwise it aborts).



 See Also: iswhitespace, read_stringbin, print_stringbin_elements

 Keywords: buffer, strings



 Examples: 3.examples.c







 */

int copy_buffer_into_stringbin(char *, STRINGBIN *);

/**







 Name: read_comment, readandstore_comment



 Synopsis:

 int read_comment(FILE *,int, long, short );

 char * readandstore_comment(FILE *,int, long);++++



 Version: 1.0





 Description: reads a comment from a turtle file. read_comment prints it on the screen, readandstore_comment

 save it in a dynamically allocated buffer pointed by the return value.







 Authors & Date: Riccardo Rigon, February, September, November 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: a pointer to the file to be read, 1 or 0, the maximum length



 of the comment to be printed. If 1 the text read is checked to be a comment



 (i.e. the program looks for the string ). If 0 the check is suppressed.



 







 Return: read_comment: 1 if successful 0 otherwise; readandstore_comment a pointer to a buffer that

 contains the data read















 See Also: iscomment, skip_whitespaces







 Examples: read_





 References:



 */

int read_comment(FILE *, int, long, short);

char * readandstore_comment(FILE *, int, long);

/**



 

 Names: read_buffer_from stdio



 Synopsis: char *read_buffer_from stdio(long maxbufferzsize)



 Version: 1.0



 Description: reads from the standard input a comment. It is opened by  and closed the same way

 this comment is.

 

 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: maxbuffersize is the length of the  maximum number of characters



 Return:  a pointer to the buffer stored





 Examples: matrixconverter.c



 Authors & Date: Riccardo Rigon, October 1997







 */

char *read_buffer_from_stdio(long);

/**







 Name: iscomment



 

 Synopsis: int iscomment(char *,FILE *);



 Version: 1.0



 Description: inspect a file a look for the string that starts



 a comment (/**)







 Authors & Date: Riccardo Rigon, February 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: a buffer where to store the comment and the pointer to the opened



 file







 Return: 1 if successful 0 otherwise





 See Also: read_comment, skip_whitespaces







 */

int iscomment(char *, FILE *);

/**



 

 Names: join_strings



 Synopsis: char * join_strings(char *,char *);



 Version: 1.0



 Description: joins two strings and store the results in a third string.

 

 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: the two strings to be joined.



 Return:  a pointer to the resulting string





 Examples: 1.examples.c



 Authors & Date: Riccardo Rigon, October 1997







 */

char * join_strings(char *, char *);

/**







 Name: get_filename



 Synopsis: char *get_filename(char *WORKING_DIRECTORY,char *program);



 Version: 0.8



 Description: get_filename returns the full name of a file to work with. It is intended to be a platform

 independent way to store input/output file names in a file. Thus avoiding to specify them on the standard input every

 time the program is executed (this is particularly useful when the same program with the same inputs is executed several times).

 The working directory (the variable: WORKING_DIRECTORY )is specified in the file

 $PathFile that has to be placed in the directory where the programs are (Windows95/NT or Mac OS) or in the

 current directory (Unixes). The file names are instead stored in a file whose name is $program where "program" is the

 name of the program being executed.



 Authors & Date : Riccardo Rigon, April 1998



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs:

 1) The pointer to working directory name

 2) The pointer to the name of the program being executed







 Return: the requested name





 See Also: simplefind, get_parameter, get_string





 */

char *get_filename(char *, char *);

STRINGBIN *read_filenames(char *working_directory, char *program,
		char *extension, char *position);

DOUBLEVECTOR *read_parameters(char *working_directory, char *program,
		char *extension, char *position);

/**







 Name: get_parameter



 Synopsis: double *get_parameter(char *WORKING_DIRECTORY,char *program);



 Version: 0.8



 Description: get_parameter returns a requested parameter. It is intended to be a platform

 independent way to store program parameter values in a file without having to specify them on the standard input every

 time the program is executed. The working directory (the variable: WORKING_DIRECTORY )is specified in the file

 $PathFile that has to be placed in the directory where the programs are (Windows95/NT or Mac OS) or in the

 current directory (Unixes). The file names are instead stored in a file whose name is $program where "program" is the

 name of the program being executed.  Every parameter is first read as a double but assignement to a different type

 of variable is possible



 Authors & Date: Riccardo Rigon, April 1998



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs:

 1) The pointer to working directory name

 2) The pointer to the name of the program being executed







 Return: the requested variable







 See Also: simplefind, get_filename, get_string











 */

double get_parameter(char *, char *);

/**







 Name: char *get_string



 Synopsis: char *get_string(char *WORKING_DIRECTORY,char *program);



 Version: 0.8



 Description: get_string returns a string needed by a program. It is intended to be a platform

 independent way to store strings in a file without having to specify them on the standard input every

 time the program is executed. The strings are  stored in a file whose name is $program where "program" is the

 name of the program being executed.



 Authors & Date: Riccardo Rigon, April 1998



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs:

 1) The pointer to working directory name

 2) The pointer to the name of the program being executed





 Return: the requested string



 See Also: simplefind, get_filename, get_parameter







 */

char * get_strings(char *, char *);

/**



 

 Names: get_workingdirectory



 Synopsis: char *get_workingdirectory(void );



 Version: 1.0



 Description: It asks for the working directory, i.e. the directory where inputs data are stored.

 This directory can be alternatively specified in a file $WorkingPaths to be found in the executable directory

 (for Windows or Mac systems where it is useful. In Unix systems it can be also the directory from where

 the program containing this routine is called)





 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: void



 Return:  a pointer to the resulting string (the working directory )



 Examples: 1.examples.c



 Authors & Date: Riccardo Rigon, October 1997



 */

char *get_workingdirectory(void);

/**







 Name: skip_whitespaces



 Version: 1.0



 Synopsis: void skip_whitespaces(FILE *);



 Description: skips the whitespaces in a file and exits when a



 normal character is found. White spaces are defined by the standard C routine isspace



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Authors & Date: Riccardo Rigon, February 1997





 Inputs: the pointer to the opened file



 Return: void



 See Also: read_comment, iscomment, goto_EOF







 */

void skip_whitespaces(FILE *);

/**







 Name: goto_EOF



 Synopsis: void goto_EOF(FILE *);



 Version: 1.0



 Description: skips the whitespaces in a file and exits when a



 normal character or EOF is found. It is very similar to skip_whitespaces.



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Authors & Date: Riccardo Rigon, February 1997







 Inputs: the pointer to the opened file



 Return: void



 See Also: read_comment, iscomment, goto_EOF







 */

void goto_EOF(FILE *);

/**







 Name: read_header, print_header



 

 Synopsis:

 HEADER read_header(FILE *,const char *);

 void print_header(HEADER *);



 Version: 1.0



 Description: read in a turtle file the header



 of the stored data. See the file "turtle.dat" for



 further details. print_header prints the header



 # below stands for _vectorheader, matrixheader or



 header







 Authors & Date: Riccardo Rigon, February, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: read_#: the pointer to the opened file. print_#



 the HEADER to be printed







 Return: a HEADER





 Keywords: turtle file



 Examples: 1.example.c





 */

HEADER read_header(FILE *, const char *);

void print_header(HEADER *);

/**







 Name: query_for_token, query_for_label, get_phrase



 







 Description: ask the file to receive either a known word (query_for_token)



 or an unknown sequence of words terminated by a '{' (query for label) or



 a sequence of words terminated by a given separator (get_phrase).







 Authors & Date: Riccardo Rigon, February, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: a pointer to a file, the word to be searched; a pointer to a file.



 get_phrase second input is the separator, for instance: ';'.







 Return: 1 if succesfull 0 otherwise ;  the read word or NOLABEL in case



 of failure.





 See Also: read_comment, read_header







 Keywords: turtle file













 */

int query_for_token(FILE *, const char *);

char *query_for_label(FILE *);

char *get_phrase(FILE *, const char);

/**







 Name: read_



 

 Synopsis:

 INTMATRIX *  read_intmatrix(FILE *, char *,short );



 Version: 1.0



 Description: read_* reads a data set of the specified type after having allocated an appropriate structure







 Authors & Date: Riccardo Rigon, February, September 1997





 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: a pointer to a file (FILE) where the data are,  the mode of reading and 1 if

 the header is going to be printed , 0 otherwise







 Return: the pointer to the appropriate structure where data are stored;





 See Also: read_comment, read_header.



 Keywords: turtle file







 Examples: 1.example.c





 */

INTMATRIX * read_intmatrix(FILE *, char *, short);

SHORTMATRIX * read_shortmatrix(FILE *, char *, short);

LONGMATRIX * read_longmatrix(FILE *, char *, short);

FLOATMATRIX * read_floatmatrix(FILE *, char *, short);

DOUBLEMATRIX * read_doublematrix(FILE *, char *, short);

CHARVECTOR * read_charvector(FILE *, char *, short);

INTVECTOR * read_intvector(FILE *, char *, short);

SHORTVECTOR * read_shortvector(FILE *, char *, short);

LONGVECTOR * read_longvector(FILE *, char *, short);

FLOATVECTOR * read_floatvector(FILE *, char *, short);

DOUBLEVECTOR * read_doublevector(FILE *, char *, short);

SHORTBIN * read_shortbin(FILE *, char *, short);

INTBIN * read_intbin(FILE *, char *, short);

LONGBIN * read_longbin(FILE *, char *, short);

DOUBLEBIN * read_doublebin(FILE *, char *, short);

/**







 Name: read_index



 Synopsis: long read_index(FILE *);



 Version: 1.0



 Description: reads the number of data blocks in a file. the format is:



 index{a,Nm} where a is the number of blocks and Nm is the type of file. Nm is optional.







 Authors & Date: Riccardo Rigon, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: a pointer to a file







 Return: the number of blocks or 0 if unsuccessful in reading



 See Also: read_comment, read_header



 Keywords: turtle file



 Examples: 5.example.c 6.example.c





 */

long read_index(FILE *, short);

/**







 Name: read__





 Synopsis: INTVECTOR * read_intarray(FILE *,short );



 Version: 1.0



 Description: read_** reads a data set of the specified type;





 Authors & Date: Riccardo Rigon, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: a pointer to a file





 Return: the structure read





 See Also: read_comment, read_header.





 Keywords: turtle file





 Examples: 1.example.c







 */

INTVECTOR * read_intarray(FILE *, short);

SHORTVECTOR * read_shortarray(FILE *, short);

FLOATVECTOR * read_floatarray(FILE *, short);

LONGVECTOR * read_longarray(FILE *, short);

CHARVECTOR * read_chararray(FILE *, short);

DOUBLEVECTOR * read_doublearray(FILE *, short);

STRINGBIN * read_stringarray(FILE *, short);

/**







 Name: justread_floatarray, justread_chararray, justread_longarray





 Synopsis: void justread_floatarray(FILE *,FLOATVECTOR *,short );



 Version: 1.0



 Description: works as read_** but assumes that the vector containing

 the data has been already allocated





 Authors & Date: Riccardo Rigon, November 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c





 Inputs: a pointer to a file, a pointer to the allocated strucure, PRINT or NOPRINT

 to indicate if the header specification are printed on the standard output







 Return: void





 See Also: read_comment, read_header, read__.







 Keywords: turtle file







 Examples: 5.example.c 6.example.c







 Notes:  only long, char and float array reading is implemented



 */

void justread_floatarray(FILE *, FLOATVECTOR *, short);

void justread_floatmatrix(FILE *, FLOATMATRIX *, char *mode, short print);

void justread_longarray(FILE *, LONGVECTOR *, short);

void justread_chararray(FILE *, CHARVECTOR *, short);

/**







 Name: write__elements



 Version: 1.0



 Synopsis:

 void write_floatarray_elements(FILE *,FLOATVECTOR *, long );







 Descriptiontion: write an array of the specified type on



 a specified file







 Authors & Date: Riccardo Rigon, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: the pointer to the file and the pointer to the data stored,



 the length of the line



 Return: void



 See Also: read__



 Keywords: turtle file



 Examples: 1.example.c





 */

void write_floatarray_elements(FILE *, FLOATVECTOR *, long);

void write_doublearray_elements(FILE *, DOUBLEVECTOR *, long);

void write_shortarray_elements(FILE *, SHORTVECTOR *, long);

void write_intarray_elements(FILE *, INTVECTOR *, long);

void write_longarray_elements(FILE *, LONGVECTOR *, long);

void write_chararray_elements(FILE *, CHARVECTOR *, long);

void write_stringsarray_elements(FILE *, STRINGBIN *);

/**







 Name: headercmp



 Synopsis: long headercmp(HEADER *,HEADER *);



 Version: 1.0



 Descriptiontion: compares two HEADERS. They



 are said to be equal if they have the same type,



 gender, category and name (they possibly differ



 for the number that indicate their position in



 the file.







 Authors & Date: Riccardo Rigon, September 1997



 Inputs: the pointers to the two header to compare





 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c





 Return: 0 if the HEADER are equal, a suitable number



 if some difference is found. This number is the



 result of four comparison operation. Each one returns



 zero if successful. Otherwise it is returned: 1



 if the two headers differ for type, 10 if they differ



 for gender, 1000 if they differ for category,



 10000 if they differ for name. If more that one comparison



 is unsuccessful the returned value is the sum of



 the unseccesful comparison value, i.e. 11 means that



 the two headers differs for type and gender.



 See Also: read_comment, read_header.



 Keywords: turtle file



 Examples: 1.example.c







 */

long headercmp(HEADER *, HEADER *);

/**







 Name: longvectorcmp



 Synopsis: long longvectorcmp(LONGVECTOR *, LONGVECTOR *);



 Version: 1.0



 Description: compares two lonvector.  This routine is used

 in read_*bin





 Authors & Date: Riccardo Rigon, September 1997



 Inputs: the pointers to the two header to compare



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Return: 0 if the longvector are equal, a suitable number



 if some difference is found.



 See Also: read_comment, read_header, header_cmp.



 Keywords: turtle file



 Examples: 1.example.c







 */

long longvectorcmp(LONGVECTOR *, LONGVECTOR *);

/**







 Name: read_header, read_matrix_header,read_vector_header



 Version: 1.0



 Synopsis:



 Description: read a vector, matrix or generic header



 Authors & Date: Riccardo Rigon, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: the pointer to the file where the header is stored



 and the pointer to the HEADER variable  where to write read information









 Keywords: turtle file







 Examples: 1.example.c



 See also: ../HANDMADE/TurtleFile





 */

void read_matrixheader(FILE *, HEADER *);

void read_vectorheader(FILE *, HEADER *);

/**







 Name: header_scan



 Version: 1.0



 Description: reads the first part of a header







 Authors & Date: Riccardo Rigon, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: the pointer to the file where the header is stored



 and the pointer to the HEADER variable where to write read information







 Return: 1 if successfull





 See Also: read_header.



 Keywords: turtle file



 Examples: 1.example.c





 */

int header_scan(FILE *, HEADER *);

short header_name_scan(FILE *, HEADER);

/**







 Name: write_turtle



 Synopsis: void write_turtle(FILE * ,char *,char *);



 Version: 1.0



 Description: writes: This is a turtle file and the date on a file. It also write:

 1) the name of the programs that created the file and 2) the files whose elaboration generated the file







 Authors & Date: Riccardo Rigon, September 1997, November 1987



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: the pointer to the file where to write and 1) and 2) above







 Return: void





 See Also: read_comment







 Keywords: turtle file







 Examples: 1.example.c.









 */

void write_turtle(FILE *, char *, char *);

/**







 Name: write_comment



 Synopsis: void write_comment(FILE *,const char *,long);



 Version: 1.0



 Description: It writes: "This is a turtle file" and the date on a file





 Authors & Date: Riccardo Rigon, September 1997



 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c



 Inputs: the pointer to the file where to write, the string to write and the line\



 length







 Return: void





 See Also: read_comment







 Keywords: turtle file







 Examples:1.example.c





 */

void write_comment(FILE *, const char *, long);

/**







 Name: write_header_header



 Synopsis: void write_header_header(FILE *,HEADER *H);



 Version: 1.0



 Description: writes a HEADER CONTENT to a file





 Authors & Date: Riccardo Rigon, September 1997





 FILE: LIBRARIES/BASICS/t_io.h, LIBRARIES/BASICS/t_io.c





 Inputs: the pointer to the file where to write and a pointer to the header



 where the information are stored







 Return: void





 See Also: read_header





 Keywords: turtle file



 Examples: 1.example.c



 */

void write_header_header(FILE *, HEADER *H);

long write_doubletensor_elements(FILE * output, DOUBLETENSOR *m, long maxcols);

void print_doubletensor_elements(DOUBLETENSOR *m, long maxcols);

DOUBLETENSOR *read_doubletensor(FILE *inputfile, char *mode, short print);

long read_doubletensor_elements(FILE *input, DOUBLETENSOR *m, char *mode);

void read_tensorheader(FILE *inputfile, HEADER *h);

