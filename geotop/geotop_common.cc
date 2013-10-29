#include <geotop_common.h>

const double TZ = 1;

long number_novalue;
long number_absent;
char *string_novalue;

// T_INIT *UV;
TInit *UV;

//char *logfile;
std::string logfile;
char **files;

long Nl = 0 ;
long Nr = 0 ;
long Nc = 0 ;

double t_meteo = 0 ;
double t_energy = 0 ;
double t_water = 0 ;
double t_sub = 0 ;
double t_sup = 0 ;
double t_blowingsnow = 0 ;
double t_out = 0 ;

double **odpnt = NULL ;
double **odp = NULL ;

long *opnt = NULL ;
long nopnt = 0;
short *ipnt = NULL ;
short *ibsn =  NULL ;
char **hpnt = NULL ;

double *odbsn = NULL ;
double *odb = NULL ;
long *obsn = NULL ;
long nobsn = 0 ;
char **hbsn = NULL;

long *osnw = NULL ;
long nosnw = 0 ;
char **hsnw = NULL ;

long *oglc = NULL ;
long noglc = 0 ;
char **hglc = NULL ;

long *osl = NULL ;
long nosl = 0 ;
char **hsl = NULL ;

FILE *ffbas = NULL ;
FILE *ffpoint = NULL ;
FILE *ffT = NULL ;
FILE *ffTav = NULL ;
FILE *ffpsi = NULL ;
FILE *ffpsitot = NULL ;
FILE *ffliq = NULL ;
FILE *ffliqav = NULL ;
FILE *ffice = NULL ;
FILE *fficeav = NULL ;
FILE *ffsnowT = NULL ;
FILE *ffsnowl = NULL ;
FILE *ffsnow = NULL ;
FILE *ffsnowi = NULL ;
FILE *ffsnowd = NULL ;
FILE *ffglac = NULL ;

long i_sim = 0 ;
long i_run = 0 ;
long i_sim0 = 0;
long i_run0 = 0 ;

//char *SuccessfulRunFile, *FailedRunFile;

std::string SuccessfulRunFile ;
std::string FailedRunFile ;

time_t start_time ; 
double elapsed_time = 0 ;
double elapsed_time_start = 0 ;
double cum_time = 0 ;
double max_time = 0 ;

