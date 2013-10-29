#ifndef __GEOTOP_COMMON__
#define __GEOTOP_COMMON__

/*----------   Global variables  ------------*/

#include "keywords.h"	//contains the definition of char** keywords_num and char** keywords_char

#include "times.h"
#include <string>

extern const double TZ ;

extern long number_novalue;
extern long number_absent;
extern char *string_novalue;

// T_INIT *UV;
extern TInit *UV;

//char *logfile;
extern std::string logfile;
extern char **files;

extern long Nl ;
extern long Nr ;
extern long Nc ;

extern double t_meteo ;
extern double t_energy ;
extern double t_water ;
extern double t_sub ;
extern double t_sup ;
extern double t_blowingsnow ;
extern double t_out ;

extern double **odpnt ;
extern double **odp;
extern long *opnt ;
extern long nopnt;
extern short *ipnt ;
extern short *ibsn;
extern char **hpnt;

extern double *odbsn ;
extern double *odb;
extern long *obsn ;
extern long nobsn;
extern char **hbsn;

extern long *osnw ;
extern long nosnw;
extern char **hsnw;

extern long *oglc ;
extern long noglc;
extern char **hglc;

extern long *osl ;
extern long nosl;
extern char **hsl;

extern FILE *ffbas ;
extern FILE *ffpoint ;
extern FILE *ffT ;
extern FILE *ffTav ;
extern FILE *ffpsi ;
extern FILE *ffpsitot ;
extern FILE *ffliq ;
extern FILE *ffliqav ;
extern FILE *ffice ;
extern FILE *fficeav ;
extern FILE *ffsnow ;
extern FILE *ffsnowT ;
extern FILE *ffsnowl ;
extern FILE *ffsnowi ;
extern FILE *ffsnowd ;
extern FILE *ffglac ;

extern long i_sim ;
extern long i_run ;
extern long i_sim0 ;
extern long i_run0 ;

//char *SuccessfulRunFile, *FailedRunFile;

extern std::string SuccessfulRunFile ;
extern std::string FailedRunFile ;

extern time_t start_time ; 
extern double elapsed_time ;
extern double elapsed_time_start ;
extern double cum_time ;
extern double max_time ;

#endif //__GEOTOP_COMMON__
