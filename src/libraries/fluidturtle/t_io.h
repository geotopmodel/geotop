#ifndef _LIBRARIES_FLUIDTURTLE_T_IO_H
#define _LIBRARIES_FLUIDTURTLE_T_IO_H




#include<cstdio>


FILE *t_fopen(const char *, const char *);


FILE *t_fclose(FILE *stream);


char *join_strings(const char *, const char *);



#endif
