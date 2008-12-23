#define   FL     512

/**



Name: search_array

Synopsis: long search_array(FILE * input,char *type, char *category)

Version: 0.8

Description: It searches array looks for the next array in file
having the type  and the category specified. Here as "array" is meant also
the header of a matrix or a vector

Authors & date: Riccardo Rigon, February 1998

FILE: LIBRAIRES/BASICS/t_utilities.h, LIBRARIES/BASICS/utilities.c

Inputs: 
1) The pointer to the opened file; 
2) The type of the array
3) The category of the array



Return: the position in file of the searched array

See Also: search_named_array, simple_find
 

*/

long search_array(FILE *,char *, char *);

/**



Name: search_named_array

Synopsis: long search_named_array(FILE * input,char *name)

Version: 0.8

Description: search named array looks for the next array in file
having the name specified. Here as "array" is meant also
the header of a matrix or a vector


Authors & date: Riccardo Rigon, February 1998

FILE: LIBRAIRES/BASICS/t_utilities.h, LIBRARIES/BASICS/utilities.c

Inputs: 
1) The pointer to the opened file; 
2) The name of the array



Return: the position in file of the searched array


See Also: search_array, simple_find


*/

long search_named_array(FILE *,char* );

/**



Name: long simplefind

Synopsis: long simplefind(FILE *input,const char *string);

Version: 0.8

Description: simplefind look for the next occurence in input file  of the "string"

FILE: LIBRAIRES/BASICS/t_utilities.h, LIBRARIES/BASICS/utilities.c

Authors & date: Riccardo Rigon, February 1998

FILE: LIBRAIRES/BASICS/t_utilities.h, LIBRARIES/BASICS/utilities.c


Inputs: 
1) The pointer to the opened file
2) The pointer to the searched string



Return: the position in input file of the found array


See ALso: search_array, search_named_array



*/

long simplefind(FILE *,const char *);

/**

Name: file_copy

Version: 0.9

Synopsis: short file_copy(FILE *destination,FILE *origin)


Description: It copies the file "origin" in file "destination"

Authors & Date: Riccardo Rigon, June 1998

FILE: LIBRAIRES/BASICS/t_utilities.h, LIBRARIES/BASICS/utilities.c

Inputs: A pointer to the  destination FILE, a pointer to the origin FILE 

*/

short file_copy(FILE *,FILE *);

/**

Name: ssimplefindkeyword

Version: 0.95

Synopsis: long ssimplefindkeyword(char *buffer,long bufferlength, const char *string)

Description:  It finds a keyword contained in "string" in a "buffer" of length "bufferlength". It also
mark the position of the buffer preceding the returned position with a null character so that finally "buffer" is
subdivided in many substrings.

Authors & Date: Riccardo Rigon, June 1998, December 1999.

Inputs: A pointer to the buffer, the buffer length, the string to be searched

Return: the position in buffer of the searched keyword 

File: LIBRARIES/UTILITIES/utilities.c,   LIBRARIES/UTILITIES/t_utilities.h

Examples: make_doc

See Also:  ../HANDMADE/documentation.html



*/

long ssimplefindkeyword(char *,long,const char *);


/*--------------------------------------------------------------------------------------*/

/**

Name: meter

Synopsis: void meter( long  index, long rows,  short frequence,const char* message, const char* separator);


Version: 1.0

Description:  It can be used to print to the video a message saying how much of a cycle is already done


Authors & Date: Riccardo Rigon, 1999

FILE: LIBRARIES/BASICS/t_utilities.h, LIBRARIES/BASICS/utilities.c

Inputs:   1- the variable  containing the cycle position; 2- the number of rows of a matrix or the total number of elements to be
parsed; 3- How many messages are requested to be output; 4- the message to be print; 5- the separator between successive messages.
i.e \n, \t and so on

Examples: tca

*/

void meter( long  index, long rows,  short frequence,const char* message, const char* separator);

/**

Name: join_strings_into, join_3strings_into

Synopsis: 
short join_strings_into(char *str, char *first,char *second);
short join_3strings_into(char *str, char *first,char *second,char* third);

Version: 1.0

Description:  joins two, or three strings -first and second ... and third -  into str

Authors & Date: Riccardo Rigon, 1999

FILE: LIBRARIES/BASICS/t_utilities.h, LIBRARIES/BASICS/utilities.c

Inputs:   1- a pointer vector opf characters containing the final string; 2- the first string; 3-the second string;
4 - the third string

Notes: Strings str is assumed to have the rigth dimensions to contain the other strings.

Examples: make_doc

*/

short join_strings_into(char *string, char *first,char *second);


short join_3strings_into(char *string, char *first,char *second,char* third);

/**


Name:   stop_execution

Synopsis:  
void stop_execution(void);

Version: 0.96

Description: top the execution of a program an wait for a key to bey pressed. It can be used in
old style debugging to avoid to declare unuseful variables in routines


Authors &  Date: Riccardo Rigon, June 2000


Inputs: void

Return: void

FIle: LIBRARIES/BASICS/utilities.c

See Also: variance_doublematrix_column


*/

void stop_execution(void);

/*--------------------------------------------------------------------------*/

void time2date(float time, long *giulian, long *year, long *month, long *day, long *hour, long *min, float *sec);

/* given a inputs 
	1:the time in second 
	2:date (giulian day, year, month, day, hour, min, sec)
   return as outputs the date updated for the time given
   	4:date (giulian day, year, month, day, hour, min, sec)
   bug: time have to be less than 1 year */

/*--------------------------------------------------------------------------*/

void giulian2day(long giulian, long year, long *month, long *day);

/* given a inputs 
	giulian day, year
   return as outputs the date
   	month, day */


/*--------------------------------------------------------------------------*/

void day2giulian(long year, long month, long day, long *giulian);

/* given a inputs 
	day, year, month
   return as outputs the julian day */

