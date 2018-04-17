#define   FL     512


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

void meter( long  index, long rows,  short frequence,const char *message,
            const char *separator);


/* given a inputs
  1:the time in second
  2:date (giulian day, year, month, day, hour, min, sec)
   return as outputs the date updated for the time given
    4:date (giulian day, year, month, day, hour, min, sec)
   bug: time have to be less than 1 year */

/* given a inputs
  giulian day, year
   return as outputs the date
    month, day */


/* given a inputs
  day, year, month
   return as outputs the julian day */

double Fmin(double a, double b);
long Fminlong(long a, long b);
double Fmax(double a, double b);
long Fmaxlong(long a, long b);
