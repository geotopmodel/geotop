#include <geotop_common.h>

std::string geotop::common::Variables::WORKING_DIRECTORY ;

std::vector<std::string> geotop::common::Variables::hpnt ;
std::vector<std::string> geotop::common::Variables::hbsn ;
std::vector<std::string> geotop::common::Variables::hsnw ;
std::vector<std::string> geotop::common::Variables::hglc ;
std::vector<std::string> geotop::common::Variables::hsl ;

const double geotop::common::Variables::TZ = 1;

TInit *geotop::common::Variables::UV;

std::string geotop::common::Variables::logfile;
std::vector<std::string> geotop::common::Variables::files;

long geotop::common::Variables::Nl = 0 ;
long geotop::common::Variables::Nr = 0 ;
long geotop::common::Variables::Nc = 0 ;

double geotop::common::Variables::t_meteo = 0 ;
double geotop::common::Variables::t_energy = 0 ;
double geotop::common::Variables::t_water = 0 ;
double geotop::common::Variables::t_sub = 0 ;
double geotop::common::Variables::t_sup = 0 ;
double geotop::common::Variables::t_blowingsnow = 0 ;
double geotop::common::Variables::t_out = 0 ;

double **geotop::common::Variables::odpnt = NULL ;
double **geotop::common::Variables::odp = NULL ;

long *geotop::common::Variables::opnt = NULL ;
long geotop::common::Variables::nopnt = 0;
short *geotop::common::Variables::ipnt = NULL ;
short *geotop::common::Variables::ibsn =  NULL ;

double *geotop::common::Variables::odbsn = NULL ;
double *geotop::common::Variables::odb = NULL ;
long *geotop::common::Variables::obsn = NULL ;
long geotop::common::Variables::nobsn = 0 ;

long *geotop::common::Variables::osnw = NULL ;
long geotop::common::Variables::nosnw = 0 ;

long *geotop::common::Variables::oglc = NULL ;
long geotop::common::Variables::noglc = 0 ;

long *geotop::common::Variables::osl = NULL ;
long geotop::common::Variables::nosl = 0 ;

FILE *geotop::common::Variables::ffbas = NULL ;
FILE *geotop::common::Variables::ffpoint = NULL ;
FILE *geotop::common::Variables::ffT = NULL ;
FILE *geotop::common::Variables::ffTav = NULL ;
FILE *geotop::common::Variables::ffpsi = NULL ;
FILE *geotop::common::Variables::ffpsitot = NULL ;
FILE *geotop::common::Variables::ffliq = NULL ;
FILE *geotop::common::Variables::ffliqav = NULL ;
FILE *geotop::common::Variables::ffice = NULL ;
FILE *geotop::common::Variables::fficeav = NULL ;
FILE *geotop::common::Variables::ffsnowT = NULL ;
FILE *geotop::common::Variables::ffsnowl = NULL ;
FILE *geotop::common::Variables::ffsnow = NULL ;
FILE *geotop::common::Variables::ffsnowi = NULL ;
FILE *geotop::common::Variables::ffsnowd = NULL ;
FILE *geotop::common::Variables::ffglac = NULL ;

long geotop::common::Variables::i_sim = 0 ;
long geotop::common::Variables::i_run = 0 ;
long geotop::common::Variables::i_sim0 = 0;
long geotop::common::Variables::i_run0 = 0 ;

//char *SuccessfulRunFile, *FailedRunFile;

std::string geotop::common::Variables::SuccessfulRunFile ;
std::string geotop::common::Variables::FailedRunFile ;

time_t geotop::common::Variables::start_time ;
double geotop::common::Variables::elapsed_time = 0 ;
double geotop::common::Variables::elapsed_time_start = 0 ;
double geotop::common::Variables::cum_time = 0 ;
double geotop::common::Variables::max_time = 0 ;
